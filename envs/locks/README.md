# Exact environment records

No tested Linux explicit locks are claimed yet. The checked-in YAMLs are pinned candidate requirements, not complete binary locks. Windows numerical dependencies have been exercised; cooltools needs a compiler here, and Linux/R remain separate validation gates.

After a successful Linux solve and fixture run, export `conda list -p PREFIX --explicit` to `core-linux-64.explicit.txt` (similarly `hicrep`/`diff-r`), record `conda list --json`, Python `pip freeze --all`, R `sessionInfo()` and package-library versions. For pip dependencies, download the exact wheel with `pip download --no-deps`, record SHA-256 and install using `--require-hashes`; a Conda explicit export alone does not reconstruct pip HiCRep. Store the companion inventory JSON with `package_digest`, and register the actual platform and verification results. `ENVIRONMENT_MODE=locked` fails without these records; it never falls back to a fresh solve.

The HiCRep explicit lock also needs `LOCK.pip.txt`, using entries such as `hicrep==0.2.6 --hash=sha256:ACTUAL_VERIFIED_WHEEL_HASH`. Record every pip-only transitive dependency with its own version/hash; `--no-deps` prevents hidden resolution. Make artifact URLs/cache files accessible on the target host. Do not copy the example placeholder hash into a real lock.

Generate the inventory after verifying a recreated prefix with the same registered Conda executable:

```bash
conda list -p /owned/verified/prefix --json > /owned/locks/packages.json
python3 scripts/lock_inventory.py --packages /owned/locks/packages.json \
  --out /owned/locks/core-linux-64.explicit.txt.inventory.json
```

This utility records a digest only; it cannot certify numerical tests. Keep the passing Linux fixture/CI evidence alongside the inventory. A Conda command-format change can change the inventory digest and requires an explicit reviewed update.
