"""Explicit chromosome classes and a single analysis-scope policy."""
from .contracts import digest, require


def chromosome_selection(c, ref, available, viewed):
    available, viewed = list(available), list(dict.fromkeys(viewed))
    classes = ref.get("chromosome_classes")
    declared = classes is not None
    if classes is None:
        assembly = ref.get("assembly", ref["reference_id"])
        count = 19 if assembly in ("mm10", "mm39") else 22 if assembly in ("hg19", "hg38", "GRCh37", "GRCh38") else None
        require(count is not None or c.get("CHROMOSOME_SCOPE") == "explicit",
                "Custom references require chromosome_classes or CHROMOSOME_SCOPE=explicit")
        classes = {"autosomes": [f"chr{i}" for i in range(1, (count or 0)+1)],
                   "sex": [x for x in ("chrX", "chrY") if x in available]}
    flat = [x for names in classes.values() for x in names]
    require(all(isinstance(names, list) and all(isinstance(x, str) for x in names) for names in classes.values()),
            "Chromosome classes must contain lists of exact names")
    require(len(flat) == len(set(flat)), "Duplicate/overlapping chromosome classes")
    if declared:
        require(set(flat) <= set(available), "Chromosome class names absent from chromosome sizes")
    mode = c.get("CHROMOSOME_SCOPE", "autosomes")
    sex = [x for x in c.get("SEX_CHROMOSOMES", "").split(",") if x]
    explicit = [x for x in c.get("ANALYSIS_CHROMOSOMES", "").split(",") if x]
    if mode == "explicit":
        require(explicit and not sex, "Explicit scope requires ANALYSIS_CHROMOSOMES and an empty SEX_CHROMOSOMES")
        selected = explicit
    else:
        require(not explicit, "ANALYSIS_CHROMOSOMES is only allowed with explicit scope")
        selected = list(classes.get("autosomes", []))
        if mode == "autosomes_plus_sex":
            require(sex and set(sex) <= set(classes.get("sex", [])), "Select exact registered SEX_CHROMOSOMES (X, Y or both)")
            selected += sex
        else:
            require(mode == "autosomes" and not sex, "Autosome scope requires empty SEX_CHROMOSOMES")
    require(selected and len(selected) == len(set(selected)), "Empty or duplicate chromosome selection")
    require(set(selected) <= set(available) and set(selected) <= set(viewed),
            "A requested chromosome is absent from registered sizes/views; no silent scope reduction")
    selected = [x for x in available if x in selected]
    identity = {"mode": mode, "chromosomes": selected, "reference_id": ref["reference_id"]}
    return dict(identity, digest=digest(identity), excluded=[x for x in available if x not in selected])
