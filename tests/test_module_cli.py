import importlib.util
import os
import subprocess
import sys
from pathlib import Path

from microc2compare.contracts import file_hash, read_json
from microc2compare.state import complete

ROOT = Path(__file__).resolve().parents[1]


def test_fresh_named_qc_without_comparison_or_phasing_assets_and_resume(tmp_path):
    spec=importlib.util.spec_from_file_location("generator", ROOT/"scripts/generate_demo.py")
    generator=importlib.util.module_from_spec(spec); spec.loader.exec_module(generator)
    data=tmp_path/"data"; generator.generate(data, development=True, raw_only=True)
    (data/"comparisons.tsv").unlink(); (data/"phasing.bedgraph").unlink()
    before={p.name:file_hash(p) for p in data.glob("*.mcool")}
    config=data/"project.conf"
    with config.open("a") as handle:
        handle.write('\nTOTAL_CPUS="2"\nMAX_MEMORY_MB="2000"\n')
        for kind in ("QC", "PREFLIGHT", "REPORT"):
            handle.write(f'THREADS_{kind}="1"\nMEMORY_MB_{kind}="1000"\n')
    bash=r"C:/Program Files/Git/bin/bash.exe" if os.name=="nt" else "bash"
    env=dict(os.environ, PYTHON_BOOTSTRAP=sys.executable, PYTHONDONTWRITEBYTECODE="1")
    command=[bash, str(ROOT/"scripts/microc2compare.sh"), "-c", str(config), "--analysis", "qc"]
    first=subprocess.run(command, env=env, capture_output=True, text=True)
    assert first.returncode==0, first.stdout+first.stderr
    run=data/"results/demo/synthetic_demo/run001"; c=read_json(run/"00_metadata/config.json")
    assert complete(c,"report","run") and not c["RUN_HICREP"]
    assert not (run/"03_features/loops").exists()
    second=subprocess.run(command+["--resume"], env=env, capture_output=True, text=True)
    assert second.returncode==0, second.stdout+second.stderr
    assert "REUSED qc/C1" in (run/"logs/qc.C1.log").read_text()
    assert before=={p.name:file_hash(p) for p in data.glob("*.mcool")}
    html=(run/"08_reports/index.html").read_text(encoding="utf-8")
    assert "Registered input files" in html and "Analysis settings" in html
