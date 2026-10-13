#!/usr/bin/env python3
import argparse
import json
import sys
import tempfile
from pathlib import Path
sys.stdout.reconfigure(newline="\n")
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from microc2compare.contracts import validate_metadata
from microc2compare.state import inventory

p=argparse.ArgumentParser();p.add_argument("kind",choices=("required","core","hicrep","loops","stripes","browser"));p.add_argument("--config",required=True);a=p.parse_args()
c=json.loads(Path(a.config).read_text())
if a.kind=="required":
    kinds=["core"]
    if c["RUN_CHROMOSIGHT"] or c.get("RUN_CHROMOSIGHT_DISCOVERY") or c["RUN_LOOPS"] and c["LOOP_CALLER"]=="mustache":kinds.append("loops")
    if c.get("RUN_BROWSER_EXPORT") and c.get("EXPORT_HIGLASS") and c.get("BROWSER_ENV_PREFIX"):kinds.append("browser")
    if c.get("RUN_STRIPES"):kinds.append("stripes")
    if c["RUN_HICREP"]:kinds.append("hicrep")
    if c["RUN_DIFFERENTIAL"] and any(x["inference"]!="descriptive" for x in validate_metadata(c)["comparisons"]):kinds.append("diff_r")
    print("\n".join(kinds))
elif a.kind=="browser":
    from microc2compare.browser_runtime import probe
    print(json.dumps(probe()))
elif a.kind=="stripes":
    from microc2compare.stripes import probe
    print(json.dumps(probe()))
elif a.kind=="loops":
    from microc2compare.tooling import probe_command
    import subprocess
    # Already entered the loop prefix by bootstrap; avoid a redundant nested Conda run.
    command=probe_command(dict(c,ENVIRONMENT_MODE="development"))
    subprocess.run(command,check=True)
else:
    import cooler,h5py,numpy as np,pandas as pd
    with tempfile.TemporaryDirectory(prefix="microc2compare-probe-") as temporary:
        path=Path(temporary)/"probe.cool"
        bins=pd.DataFrame({"chrom":["chrTest"]*3,"start":[0,10000,20000],"end":[10000,20000,30000]})
        pixels=pd.DataFrame({"bin1_id":[0,0,1],"bin2_id":[1,2,2],"count":[2,3,4]})
        cooler.create_cooler(str(path),bins,pixels,ordered=True)
        assert cooler.Cooler(str(path)).matrix(balance=False)[0:1,1:2][0,0]==2
        with h5py.File(path,"r+") as handle:handle["bins"].create_dataset("weight",data=[2.,3.,4.])
        assert cooler.Cooler(str(path)).matrix(balance=True)[0:1,1:2][0,0]==12
        if a.kind=="hicrep":
            from hicrep import hicrepSCC
            score=hicrepSCC(cooler.Cooler(str(path)),cooler.Cooler(str(path)),0,30000,False)
            assert np.isfinite(score).all() and np.allclose(score,1)
        elif c["RUN_COMPARTMENTS"] or c["RUN_INSULATION"] or c["RUN_LOOPS"]:
            import cooltools
            assert all(hasattr(cooltools,name) for name in ("eigs_cis","insulation","dots"))
    print(json.dumps(inventory(),indent=2))
