#!/usr/bin/env python3
"""Verify an archive completely before installing a new shared source release.

No environment, project, previous source release or command link is modified.
"""
import argparse
import hashlib
import json
import re
import stat
import sys
import zipfile
from pathlib import Path, PurePosixPath


def install(archive_path, expected_sha256, base, release_id):
    archive_path=Path(archive_path).resolve(); base=Path(base).resolve()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*",release_id): raise ValueError("Unsafe release ID")
    if hashlib.sha256(archive_path.read_bytes()).hexdigest()!=expected_sha256: raise ValueError("Archive SHA256 mismatch")
    with zipfile.ZipFile(archive_path) as archive:
        names=archive.namelist()
        if len(names)!=len(set(names)): raise ValueError("Duplicate archive member")
        manifests=[n for n in names if n.endswith("/RELEASE_MANIFEST.json")]
        if len(manifests)!=1: raise ValueError("Expected one release manifest")
        prefix=manifests[0].split("/")[0]
        manifest=json.loads(archive.read(manifests[0]))
        if prefix!="microc2compare-"+manifest["version"]: raise ValueError("Release version/root mismatch")
        expected={prefix+"/"+name for name in manifest["files"]} | {manifests[0]}
        if set(names)!=expected: raise ValueError("Archive content differs from allowlisted manifest")
        contents={}
        for info in archive.infolist():
            path=PurePosixPath(info.filename)
            if path.is_absolute() or ".." in path.parts or "\\" in info.filename or ":" in info.filename: raise ValueError("Unsafe archive path")
            if stat.S_ISLNK(info.external_attr>>16): raise ValueError("Archive symlinks are unsupported")
            data=archive.read(info.filename); relative=str(path.relative_to(prefix))
            if relative!="RELEASE_MANIFEST.json" and hashlib.sha256(data).hexdigest()!=manifest["files"][relative]: raise ValueError("Source checksum mismatch: "+relative)
            contents[relative]=data
    target=base/"releases"/release_id/prefix
    if not target.resolve().is_relative_to(base): raise ValueError("Resolved release path escapes the shared base")
    if (base/"releases"/release_id).exists(): raise ValueError("Release ID exists; refusing overwrite")
    # All paths/content were checked before the first write.
    target.mkdir(parents=True,mode=0o755)
    for relative,data in contents.items():
        path=target/relative; path.parent.mkdir(parents=True,exist_ok=True,mode=0o755)
        path.write_bytes(data); path.chmod(0o644)
    for directory in [target.parent,target,*[p for p in target.rglob("*") if p.is_dir()]]: directory.chmod(0o755)
    print(f"Verified and installed {len(manifest['files'])} files\nSource: {target}")
    return target


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive",required=True);p.add_argument("--sha256",required=True)
    p.add_argument("--base",default="/opt/conda_envs/microc2compare");p.add_argument("--release-id",required=True)
    a=p.parse_args()
    try:install(a.archive,a.sha256,a.base,a.release_id)
    except (ValueError,OSError,KeyError,zipfile.BadZipFile) as error:print(f"ERROR: {error}",file=sys.stderr);sys.exit(2)
