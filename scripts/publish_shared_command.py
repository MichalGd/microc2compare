#!/usr/bin/env python3
"""Publish a versioned launcher; activate only with an explicit --activate flag."""
import argparse
import hashlib
import json
import os
import re
import shlex
import stat
import sys
from pathlib import Path


def publish(source, base, activate=False, site_profile=None):
    source=Path(source).resolve();base=Path(base).resolve()
    if not source.is_relative_to(base/"releases"):raise ValueError("Source must be a release below the shared base")
    manifest=json.loads((source/"RELEASE_MANIFEST.json").read_text())
    for name,expected in manifest["files"].items():
        path=(source/name).resolve()
        if not path.is_relative_to(source) or hashlib.sha256(path.read_bytes()).hexdigest()!=expected:raise ValueError("Source manifest mismatch: "+name)
    name="microc2compare-"+manifest["version"]+"-"+source.parent.name
    if not re.fullmatch(r"[A-Za-z0-9_.-]+",name):raise ValueError("Unsafe launcher name")
    binding=''
    if site_profile is not None:
        profile=Path(site_profile).resolve()
        if not profile.is_relative_to((base/"site").resolve()) or not profile.is_file():
            raise ValueError("Shared software profile must be a file below the shared site directory")
        if os.name!="nt" and any(p.stat().st_mode & (stat.S_IWGRP|stat.S_IWOTH) for p in (profile,profile.parent)):
            raise ValueError("Shared software profile/directory is writable by non-maintainers")
        checksum=hashlib.sha256(profile.read_bytes()).hexdigest()
        binding='export M2C_DEFAULT_SOFTWARE_PROFILE='+shlex.quote(profile.as_posix())+'\n'
        binding+='export M2C_DEFAULT_SOFTWARE_PROFILE_SHA256='+shlex.quote(checksum)+'\n'
    directory=base/"bin";directory.mkdir(exist_ok=True,mode=0o755)
    if os.name!="nt" and directory.stat().st_mode & (stat.S_IWGRP|stat.S_IWOTH):raise ValueError("Shared bin directory is writable by non-maintainers")
    launcher=directory/name
    text='#!/usr/bin/env bash\nset -Eeuo pipefail\nexport PYTHONDONTWRITEBYTECODE=1\n'+binding+'exec bash '+shlex.quote((source/"scripts/microc2compare.sh").as_posix())+' "$@"\n'
    if launcher.exists():
        if launcher.is_symlink() or launcher.read_text()!=text:raise ValueError("Existing launcher differs; refusing overwrite")
    else:
        with launcher.open("x",newline="\n") as handle:handle.write(text)
        launcher.chmod(0o755)
    if activate:
        alias=directory/"microc2compare"
        if alias.exists() and not alias.is_symlink():raise ValueError("Alias is not a symlink; refusing replacement")
        pending=directory/".microc2compare.next"
        if pending.exists() or pending.is_symlink():raise ValueError("Pending alias exists; inspect it before retrying")
        pending.symlink_to(name);os.replace(pending,alias)
    print(f"Versioned command: {launcher}\nShared alias activated: {activate}")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--source",required=True)
    p.add_argument("--base",default="/opt/conda_envs/microc2compare");p.add_argument("--activate",action="store_true")
    p.add_argument("--site-profile",help="Pin this maintainer-owned software profile into the versioned launcher")
    a=p.parse_args()
    try:publish(a.source,a.base,a.activate,a.site_profile)
    except (ValueError,OSError,KeyError) as error:print(f"ERROR: {error}",file=sys.stderr);sys.exit(2)
