#!/usr/bin/env python3
"""Create an explicitly allowlisted public source archive; does not publish."""
import argparse
import hashlib
import json
import zipfile
import tomllib
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
VERSION=tomllib.loads((ROOT/'pyproject.toml').read_text())['project']['version']
ARCHIVE_ROOT='microc2compare-'+VERSION
DIRECTORIES=('src','scripts','schemas','docs','config','envs','tests','.github')
FILES=('README.md','CONTRIBUTING.md','CHANGELOG.md','LICENSE.md','CITATION.cff',
       'pyproject.toml','environment.yml','.gitignore','.gitattributes',
       'examples/tiny/README.md',
       'examples/upstream-handoff/README.md','examples/upstream-handoff/samplesheet.csv',
       'validation/calibration-protocol.md','validation/calibration_registry.json',
       'validation/datasets.tsv','validation/local-verification.md','validation/environment-solves.json',
       'validation/performance-033-local.json','validation/general-purpose-040-local.json')
EXTENSIONS={'.py','.R','.sh','.json','.yml','.yaml','.md','.conf','.csv','.tsv','.txt','.cff','.toml'}
PRIVATE_WORKSPACE_TESTS = {'test_browser_loop_tracks11oct2026.py'}

def selected():
    result=[ROOT/p for p in FILES if (ROOT/p).is_file()]
    for directory in DIRECTORIES:
        for path in (ROOT/directory).rglob('*'):
            if not path.is_file() or path.suffix not in EXTENSIONS:continue
            if any(part in ('__pycache__','.pytest_cache') for part in path.parts):continue
            # Historical cohort-export tests import private Analysis scripts.
            # General browser tests shipped below cover the public exporter.
            if directory == 'tests' and path.name in PRIVATE_WORKSPACE_TESTS:continue
            if directory=='config' and (path.name=='project.conf' or '.local.' in path.name):continue
            # Config publication is deliberately limited to templates and examples.
            if directory=='config' and not (path.name=='config_template.conf' or path.name.startswith('resources_') or '.example.' in path.name):continue
            result.append(path)
    return sorted(set(result))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',required=True);p.add_argument('--list',action='store_true');a=p.parse_args()
    paths=selected()
    if a.list:
        print('\n'.join(str(p.relative_to(ROOT)) for p in paths));return
    output=Path(a.out).resolve()
    if output.exists():raise SystemExit('Refusing to overwrite an existing release archive')
    output.parent.mkdir(parents=True,exist_ok=True)
    def public_bytes(path):
        # Normalize repository text for reproducible Linux-facing extraction.
        return path.read_bytes().replace(b'\r\n',b'\n')
    hashes={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(public_bytes(p)).hexdigest() for p in paths}
    with zipfile.ZipFile(output,'w',compression=zipfile.ZIP_DEFLATED) as archive:
        for path in paths:archive.writestr(ARCHIVE_ROOT+'/'+str(path.relative_to(ROOT)).replace('\\','/'),public_bytes(path))
        archive.writestr(ARCHIVE_ROOT+'/RELEASE_MANIFEST.json',json.dumps({'status':'validation_candidate','version':VERSION,'files':hashes},indent=2)+'\n')
    print(f'{output}: {len(paths)} source/documentation files; review before publication')

if __name__=='__main__':main()
