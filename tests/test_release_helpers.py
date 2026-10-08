import hashlib
import importlib.util
import json
import zipfile
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]


def helper(name):
    spec=importlib.util.spec_from_file_location(name,ROOT/"scripts"/(name+".py"))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def archive(path, name="scripts/microc2compare.sh", corrupt=False):
    contents=b"#!/bin/bash\necho test\n"
    manifest={"version":"0.2.0","files":{name:hashlib.sha256(contents).hexdigest()}}
    with zipfile.ZipFile(path,"w") as z:
        z.writestr("microc2compare-0.2.0/"+name,contents+b"corrupt" if corrupt else contents)
        z.writestr("microc2compare-0.2.0/RELEASE_MANIFEST.json",json.dumps(manifest))
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_verified_new_release_and_no_overwrite(tmp_path):
    module=helper("install_release"); p=tmp_path/"release.zip"; sha=archive(p)
    base=tmp_path/"software"; source=module.install(p,sha,base,"candidate")
    assert (source/"scripts/microc2compare.sh").read_bytes().startswith(b"#!/bin/bash")
    with pytest.raises(ValueError,match="refusing overwrite"):module.install(p,sha,base,"candidate")


@pytest.mark.parametrize("name,corrupt,message",[("../escape.py",False,"Unsafe"),("scripts/a.py",True,"checksum")])
def test_rejected_archive_does_not_write(tmp_path,name,corrupt,message):
    module=helper("install_release");p=tmp_path/"release.zip";sha=archive(p,name,corrupt)
    base=tmp_path/"software"
    with pytest.raises(ValueError,match=message):module.install(p,sha,base,"candidate")
    assert not base.exists()


def test_wrong_sha_and_bad_release_id_rejected(tmp_path):
    module=helper("install_release");p=tmp_path/"release.zip";archive(p)
    with pytest.raises(ValueError,match="SHA256"):module.install(p,"wrong",tmp_path/"base","candidate")
    with pytest.raises(ValueError,match="release ID"):module.install(p,"wrong",tmp_path/"base","../outside")


def test_variant_literals_and_relative_paths(tmp_path):
    import subprocess
    import sys
    source=tmp_path/"project.conf";source.write_text('SAMPLESHEET="samples.csv"\n')
    target=tmp_path/"variant.conf"
    result=subprocess.run([sys.executable,str(ROOT/"scripts/create_variant.py"),"--config",str(source),"--out",str(target),"--run-id","new","--set","INSULATION_BOUNDARY_THRESHOLD=0.5"],capture_output=True,text=True)
    assert result.returncode==0,result.stderr
    assert "RUN_ID=new" in target.read_text() and "INSULATION_BOUNDARY_THRESHOLD=0.5" in target.read_text()
    assert subprocess.run([sys.executable,str(ROOT/"scripts/create_variant.py"),"--config",str(source),"--out",str(target),"--run-id","other"],capture_output=True).returncode!=0
