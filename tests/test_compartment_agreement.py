"""Acceptance permits per-view raw signs/rounding, while preserving biological direction."""
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from microc2compare.acceptance import compare_compartments
from microc2compare.contracts import ContractError, atomic_json, read_table, write_table


def fixture(path, flip=False):
    rows = []
    for i in range(6):
        region = "left" if i < 3 else "right"
        factors = (-1, 1, -1) if region == "left" else (1, -1, -1)
        components = [float(i+1), float(3-i), float(i-2)]
        rows.append(dict(chrom="chr1", start=i*100, end=(i+1)*100,
                         **{f"E{j+1}": x*factors[j]+1e-15 if flip else x for j,x in enumerate(components)},
                         oriented_score=float(3-i)+1e-15 if flip else float(3-i), selected_component="E2"))
    write_table(path, rows, list(rows[0]))
    views = [dict(chrom="chr1",start=0,end=300,name="left",eigval1=8.,eigval2=4.,eigval3=-2.),
             dict(chrom="chr1",start=300,end=600,name="right",eigval1=7.,eigval2=3.,eigval3=-1.)]
    write_table(path.with_suffix(".eigenvalues.tsv"), views, list(views[0]))
    atomic_json(path.with_suffix(".phasing.json"), [
        dict(region="left",component="E2",correlation=.8,status="oriented"),
        dict(region="right",component="E2",correlation=-.7 if flip else .7,status="oriented")])
    path.with_suffix(".bedgraph").write_text("".join(
        f"chr1\t{row['start']}\t{row['end']}\t{row['oriented_score']}\n" for row in rows))


@pytest.fixture
def pair(tmp_path):
    full, selected = tmp_path/"full.tsv", tmp_path/"selected.tsv"
    fixture(full); fixture(selected, True)
    return full, selected


def change(path, column, value, index=0):
    rows = read_table(path)
    rows[index][column] = value
    write_table(path, rows, list(rows[0]))


def test_raw_signs_are_allowed_independently_per_view_but_not_per_bin(pair):
    full, selected = pair
    before = {str(p):p.read_bytes() for p in full.parent.iterdir() if p.is_file()}
    result = compare_compartments(full, selected)
    assert result["status"] == "numerically_equivalent" and result["rows"] == result["browser_rows"] == 6
    assert result["views"][0]["raw_component_signs"] == dict(E1=-1, E2=1, E3=-1)
    assert result["views"][1]["raw_component_signs"] == dict(E1=1, E2=-1, E3=-1)
    assert max(row["max_oriented_score_error"] for row in result["views"]) < 2e-15
    assert before == {str(p):p.read_bytes() for p in full.parent.iterdir() if p.is_file()}
    change(selected, "E1", 2., index=1)
    with pytest.raises(ContractError, match="numerical disagreement"):
        compare_compartments(full, selected)


@pytest.mark.parametrize("column,value,message", [
    ("start", 1, "start differs"), ("chrom", "chr2", "chrom differs"),
    ("selected_component", "E1", "selected_component differs"),
    ("E1", -1.001, "numerical disagreement"),
    ("oriented_score", -3., "numerical disagreement"),
    ("oriented_score", 3.001, "numerical disagreement"),
    ("E2", "", "missing-value pattern"), ("oriented_score", "", "missing-value pattern"),
    ("E3", "inf", "infinite value"), ("E3", "garbage", "invalid number")])
def test_semantic_and_numeric_changes_still_fail(pair, column, value, message):
    full, selected = pair
    change(selected, column, value)
    with pytest.raises(ContractError, match=message): compare_compartments(full, selected)


@pytest.mark.parametrize("kind", ["rows", "headers", "order", "eigenvalue", "view", "phasing", "correlation", "browser"])
def test_bundle_structure_eigenvalues_phasing_and_tracks_are_checked(pair, kind):
    full, selected = pair
    if kind in ("rows", "order"):
        rows = read_table(selected)
        rows = rows[:-1] if kind == "rows" else list(reversed(rows))
        write_table(selected, rows, list(rows[0]))
    elif kind == "headers": selected.write_text("chrom\tstart\tend\nchr1\t0\t100\n")
    elif kind == "eigenvalue": change(selected.with_suffix(".eigenvalues.tsv"), "eigval1", -8.)
    elif kind == "view": change(selected.with_suffix(".eigenvalues.tsv"), "end", 400)
    elif kind in ("phasing", "correlation"):
        path = selected.with_suffix(".phasing.json"); qc = json.loads(path.read_text())
        qc[1]["status" if kind == "phasing" else "correlation"] = "ambiguous" if kind == "phasing" else .7
        atomic_json(path, qc)
    else:
        path = selected.with_suffix(".bedgraph")
        path.write_text(path.read_text().replace("3.000000000000001", "-3.000000000000001"))
    with pytest.raises(ContractError): compare_compartments(full, selected)


def test_matching_nan_masks_and_empty_tracks_are_valid(pair):
    full, selected = pair
    for path in pair:
        change(path, "E3", "")
        rows = read_table(path)
        for row in rows: row.update(oriented_score="", selected_component="")
        write_table(path, rows, list(rows[0]))
        qc = json.loads(path.with_suffix(".phasing.json").read_text())
        for row in qc: row.update(component=None, correlation=None, status="ambiguous")
        atomic_json(path.with_suffix(".phasing.json"), qc)
        path.with_suffix(".bedgraph").write_text("# No finite measurements\n")
    result = compare_compartments(full, selected)
    assert result["browser_rows"] == 0 and result["status"] == "numerically_equivalent"


def test_real_cli_under_bootstrap_without_scientific_imports(pair):
    full, selected = pair
    script = Path(__file__).resolve().parents[1]/"scripts/compare_compartments.py"
    wrapper = '''
import builtins,runpy,sys
original=builtins.__import__
def stdlib_import(name,*args,**kwargs):
    if name.split('.')[0] in {'numpy','pandas','cooler','cooltools','h5py','bioframe','matplotlib'}:
        raise AssertionError('Unexpected scientific import: '+name)
    return original(name,*args,**kwargs)
builtins.__import__=stdlib_import
sys.argv=sys.argv[1:]
runpy.run_path(sys.argv[0],run_name='__main__')
'''
    args = [sys.executable,"-I","-S","-B","-c",wrapper,str(script),"--full",str(full),"--selected",str(selected)]
    result = subprocess.run(args,capture_output=True,text=True)
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)["status"] == "numerically_equivalent"
    change(selected, "oriented_score", -3.)
    result = subprocess.run(args,capture_output=True,text=True)
    assert result.returncode == 2 and "numerical disagreement" in result.stderr


def test_acceptance_bash_calls_numerical_checker_and_preserves_inputs(tmp_path):
    # Coordination test: pipeline/validator stubs do not establish science.
    root = Path(__file__).resolve().parents[1]
    def shell_path(path):
        return '/'+path.drive[0].lower()+'/'+path.relative_to(path.anchor).as_posix() if os.name=='nt' else str(path)
    scripts = tmp_path/"source/scripts"; scripts.mkdir(parents=True)
    (scripts/"validate_module_demo.sh").write_bytes((root/"scripts/validate_module_demo.sh").read_bytes().replace(b"\r\n",b"\n"))
    (scripts/"microc2compare.sh").write_text("#!/usr/bin/env bash\nexit 0\n",newline="\n")
    (scripts/"validate_outputs.py").write_text("# Coordination stub; no scientific validation.\n")
    (scripts/"compare_compartments.py").write_text(
        "import runpy\nrunpy.run_path("+repr(str(root/"scripts/compare_compartments.py"))+",run_name='__main__')\n")
    site = tmp_path/"site.conf"
    site.write_text('PYTHON_BOOTSTRAP="'+shell_path(Path(sys.executable))+'"\n',newline="\n")
    destination = tmp_path/"private"; destination.mkdir()
    (destination/"project.conf").write_text('source "'+shell_path(site)+'"\n',newline="\n")
    identity = destination/"original.mcool"
    identity.write_bytes(b"Identity-only fixture, not a scientific contact map")
    import hashlib
    (destination/"input-sha256.txt").write_text(hashlib.sha256(identity.read_bytes()).hexdigest()+'  '+shell_path(identity)+'\n',newline="\n")
    base = destination/"results/demo/synthetic_demo"
    for run, flipped in (("run001",False),("compartments-only",True)):
        for sid in ("C1","C2","C3","T1","T2","T3"):
            path = base/run/f"03_features/compartments/{sid}.tsv"
            path.parent.mkdir(parents=True,exist_ok=True); fixture(path,flipped)
    for run in ("run001","insulation-only","loops-only"):
        for sid in ("C1","C2","C3","T1","T2","T3"):
            for kind, name in (("insulation",sid+".tsv"),("loops",sid+".10000.bedpe.tsv")):
                path=base/run/"03_features"/kind/name
                path.parent.mkdir(parents=True,exist_ok=True)
                path.write_text("Coordination-only identical text fixture\n")
    before=identity.read_bytes()
    bash='C:/Program Files/Git/bin/bash.exe' if os.name=='nt' else 'bash'
    result=subprocess.run([bash,str(scripts/"validate_module_demo.sh"),"--site",str(site),
                           "--out",str(destination),"--resume"],capture_output=True,text=True)
    assert result.returncode==0,result.stdout+result.stderr
    assert result.stdout.count('"status": "numerically_equivalent"')==6
    assert 'Full/selected/variant acceptance passed' in result.stdout
    assert identity.read_bytes()==before
