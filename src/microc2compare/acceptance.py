"""Read-only numerical agreement checks for full/selected compartment products.

Raw eigenvectors may change sign independently within each declared view.
Oriented scores and browser tracks retain their direction. No products are edited.
"""
import csv
import math
from pathlib import Path

from .contracts import ContractError, read_json, require

ATOL = 1e-10
RTOL = 1e-8
COMPONENTS = ("E1", "E2", "E3")
AXES = ("chrom", "start", "end")


def table(path, required):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        require(fields and len(fields) == len(set(fields)) and set(required) <= set(fields),
                f"{path}: missing/duplicate required columns")
        rows = list(reader)
    require(all(None not in row and all(value is not None for value in row.values()) for row in rows),
            f"{path}: malformed row")
    return fields, rows


def number(value, label):
    try: result = float(value) if value not in (None, "", "NA") else math.nan
    except (TypeError, ValueError) as error: raise ContractError(f"{label}: invalid number {value!r}") from error
    require(not math.isinf(result), f"{label}: infinite value")
    return result


def close(a, b):
    return (math.isnan(a) and math.isnan(b)) or math.isclose(a, b, rel_tol=RTOL, abs_tol=ATOL)


def numeric_pair(left, right, label, allow_sign=False):
    require(len(left) == len(right), f"{label}: value count differs")
    a = [number(x, label) for x in left]
    b = [number(x, label) for x in right]
    require(all(math.isnan(x) == math.isnan(y) for x, y in zip(a, b)), f"{label}: missing-value pattern differs")
    signs = (1, -1) if allow_sign else (1,)
    for sign in signs:
        if all(close(sign*x, y) for x, y in zip(a, b)):
            errors = [abs(sign*x-y) for x,y in zip(a,b) if not math.isnan(x)]
            return sign, max(errors, default=0.)
    error = max((abs(x-y) for x,y in zip(a,b) if not math.isnan(x)), default=0.)
    raise ContractError(f"{label}: numerical disagreement (max direct difference {error:.6g}; atol={ATOL:g}, rtol={RTOL:g})")


def paired_tables(full, selected, required):
    fields, a = table(full, required)
    other_fields, b = table(selected, required)
    require(fields == other_fields, f"{full.name}: column names/order differ")
    require(len(a) == len(b), f"{full.name}: row count differs")
    return fields, a, b


def compare_compartments(full, selected):
    """Compare vectors, view eigenvalues, phasing QC and oriented browser tracks."""
    full, selected = Path(full), Path(selected)
    fields, a, b = paired_tables(full, selected, (*AXES, *COMPONENTS, "oriented_score", "selected_component"))
    for index, (left, right) in enumerate(zip(a, b), 2):
        for field in set(fields)-set(COMPONENTS)-{"oriented_score"}:
            require(left[field] == right[field], f"{full.name}:{index}: {field} differs")
    ef, ea, eb = paired_tables(full.with_suffix(".eigenvalues.tsv"), selected.with_suffix(".eigenvalues.tsv"),
                                (*AXES, "name", "eigval1", "eigval2", "eigval3"))
    require(ea and len({row["name"] for row in ea}) == len(ea), "Compartment views must be nonempty with unique names")
    for left, right in zip(ea, eb):
        for field in set(ef)-{"eigval1", "eigval2", "eigval3"}:
            require(left[field] == right[field], f"Compartment eigenvalue view {field} differs")
    for field in ("eigval1", "eigval2", "eigval3"):
        numeric_pair([row[field] for row in ea], [row[field] for row in eb], field)
    memberships = [[] for _ in ea]
    for index, row in enumerate(a):
        matches = [i for i, view in enumerate(ea) if row["chrom"] == view["chrom"]
                   and int(view["start"]) <= int(row["start"]) < int(row["end"]) <= int(view["end"])]
        require(len(matches) == 1, f"{full.name}: bin has missing/overlapping view membership")
        memberships[matches[0]].append(index)
    summary = {"vectors": str(full), "selected_vectors": str(selected), "rows": len(a),
               "atol": ATOL, "rtol": RTOL, "views": []}
    signs = {}
    for view, indexes in zip(ea, memberships):
        require(indexes, f"Compartment view {view['name']} has no bins")
        result = {"name": view["name"], "raw_component_signs": {}, "max_component_error_after_sign": {}}
        for component in COMPONENTS:
            sign, error = numeric_pair([a[i][component] for i in indexes], [b[i][component] for i in indexes],
                                       f"{view['name']}/{component}", allow_sign=True)
            result["raw_component_signs"][component] = sign
            result["max_component_error_after_sign"][component] = error
            signs[view["name"], component] = sign
        _, error = numeric_pair([a[i]["oriented_score"] for i in indexes], [b[i]["oriented_score"] for i in indexes],
                                f"{view['name']}/oriented_score")
        result["max_oriented_score_error"] = error
        summary["views"].append(result)
    qa, qb = read_json(full.with_suffix(".phasing.json")), read_json(selected.with_suffix(".phasing.json"))
    require(isinstance(qa, list) and isinstance(qb, list) and len(qa) == len(qb) == len(ea), "Phasing QC view count differs")
    require([q["region"] for q in qa] == [view["name"] for view in ea], "Phasing QC view identities/order differ")
    for left, right in zip(qa, qb):
        require(set(left) == set(right), "Phasing QC fields differ")
        for field in set(left)-{"correlation"}:
            require(left[field] == right[field], f"Phasing QC {field} differs")
        component = left["component"]
        if component is None:
            require(left["correlation"] is None and right["correlation"] is None, "Unselected component has correlation")
        else:
            require(component in COMPONENTS, "Unknown phasing component")
            sign = signs[left["region"], component]
            x, y = number(left["correlation"], "phasing correlation"), number(right["correlation"], "phasing correlation")
            require(close(sign*x, y), "Phasing correlation disagrees with raw eigenvector sign")
    def browser(path):
        with path.open(encoding="utf-8") as handle:
            rows = [line.rstrip("\r\n").split("\t") for line in handle if line.strip() and not line.startswith("#")]
        require(all(len(row) == 4 for row in rows), f"{path}: malformed bedGraph")
        return rows
    ba, bb = browser(full.with_suffix(".bedgraph")), browser(selected.with_suffix(".bedgraph"))
    require(len(ba) == len(bb) and all(x[:3] == y[:3] for x,y in zip(ba,bb)), "Compartment browser coordinates differ")
    require(all(math.isfinite(number(row[3], "browser score")) for row in ba+bb), "Browser track has missing scores")
    numeric_pair([row[3] for row in ba], [row[3] for row in bb], "oriented browser scores")
    summary["browser_rows"] = len(ba)
    summary["status"] = "numerically_equivalent"
    return summary
