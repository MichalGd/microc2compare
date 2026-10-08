# Windows development verification

Windows/Git Bash is a development environment; the production target remains Linux. This workspace used a project-local embedded Python 3.11.9 because the WindowsApps Python alias was unusable. The runtime and downloaded tools stay under ignored `.devtools/` and are excluded from the public archive.

Use PowerShell with an actual Python executable and Git Bash path. These example commands match the local verification setup:

```powershell
Set-Location 'D:\path\to\microc2compare'
& '.\.devtools\python\python.exe' -m pytest -q
$env:PYTHON_BOOTSTRAP = (Resolve-Path '.\.devtools\python\python.exe').Path
& 'C:\Program Files\Git\bin\bash.exe' scripts/ci_checks.sh
& 'C:\Program Files\Git\bin\bash.exe' scripts/microc2compare.sh -c examples/tiny/generated/project.conf --resume
```

For a fresh empty fixture destination:

```powershell
& '.\.devtools\python\python.exe' scripts/generate_demo.py --out examples/tiny/new-demo --development --raw-only
& 'C:\Program Files\Git\bin\bash.exe' scripts/microc2compare.sh -c examples/tiny/new-demo/project.conf
```

The generator embeds the actual Python path in the local development configuration. It creates synthetic maps with known test weights; this is not production balancing. The raw-only demonstration exercises QC, HiCRep and reports without claiming omitted cooltools/R modules passed.

Cooler/HiCRep numerical dependencies installed and ran locally. cooltools 0.7.1 source installation reached the extension build and failed because Microsoft Visual C++ 14+ was unavailable. At initial verification, Rscript and a working Linux runtime were unavailable; WSL was disabled and the Docker Linux daemon was stopped. No system service/compiler installation was performed. A Windows micromamba cross-solve verified Linux dependency availability, not Linux binaries executing.

For the 0.2.4 loop correction, an isolated ignored `.devtools/r-adapter-test` runtime was created inside this workspace: R 4.4.3, edgeR 4.4.2, limma 3.62.2, statmod 1.5.2 and jsonlite 2.0.0. Actual loop/insulation adapter tests and base-R uncertainty regressions passed. Set `M2C_TEST_RSCRIPT` to the accepted Rscript executable and ensure that runtime's DLL directories are on the process PATH to run those tests. This does not change any server or upstream environment. It does not establish Linux acceptance, calibration, or compatibility of untested package versions; see [the correction guide](24-r-loop-compatibility.md).

Do not translate Windows numerical success into Linux process-group/R/resource certification. Use the separate Linux CI/full fixture, R numerical tests and server pilot described in the other guides.
