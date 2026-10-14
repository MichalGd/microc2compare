"""Explicit culture observations preserve donor identity and inference limits."""
import pytest

from microc2compare.contracts import ContractError, SAMPLE_FIELDS, validate_metadata, write_table


def metadata(tmp_path, observations, inference="descriptive", design="~condition"):
    samples = []
    for index, observation in enumerate(observations):
        sid = f"S{index}"
        sample = dict(sample_id=sid, project="fresh", condition="control", assay="microc",
                      biological_replicate=str(index), subject_id=sid, reference_id="customAssembly",
                      matrix_path=f"{sid}.cool", provenance_path=f"{sid}.json")
        samples.append(dict(sample, **observation))
    fields = list(SAMPLE_FIELDS) + ["experimental_unit_id", "culture_id"]
    sample_path = tmp_path / "samples.tsv"
    write_table(sample_path, samples, fields)
    membership_path = tmp_path / "membership.tsv"
    write_table(membership_path, [dict(sample_id=s["sample_id"], library_id=f"L{i}", source_id=f"source{i}",
                                      source_kind="independent_library", dedup_group_id=f"group{i}",
                                      dedup_scope="joint_library", source_checksum=f"checksum{i}")
                                   for i, s in enumerate(samples)],
                ["sample_id", "library_id", "source_id", "source_kind", "dedup_group_id", "dedup_scope", "source_checksum"])
    comparisons_path = tmp_path / "comparisons.tsv"
    write_table(comparisons_path, [dict(comparison_id="effect", case_condition="treated", reference_condition="control",
                                       design=design, inference=inference, modules="loops")],
                ["comparison_id", "case_condition", "reference_condition", "design", "inference", "modules"])
    return dict(SAMPLESHEET=str(sample_path), TECHNICAL_MEMBERSHIP=str(membership_path),
                COMPARISONS=str(comparisons_path), PROJECT="fresh", REFERENCE_ID="customAssembly",
                ANALYSIS_MODE="condition_comparison", RUN_DIFFERENTIAL=True, RUN_TADS=False,
                FORMAL_MIN_UNITS_PER_CONDITION=3, PRELIMINARY_MIN_UNITS_PER_CONDITION=2)


def test_descriptive_independent_cultures_keep_shared_donor(tmp_path):
    c = metadata(tmp_path, [dict(subject_id="donorA", experimental_unit_id="cultureA", condition="control"),
                            dict(subject_id="donorA", experimental_unit_id="cultureB", condition="control"),
                            dict(subject_id="donorB", culture_id="cultureC", condition="treated")])
    result = validate_metadata(c)
    assert [s["subject_id"] for s in result["samples"]] == ["donorA", "donorA", "donorB"]
    assert result["comparisons"][0]["inference"] == "descriptive"


@pytest.mark.parametrize("units", [{}, {"experimental_unit_id": "sameCulture"}])
def test_repeated_donor_condition_needs_unique_explicit_units(tmp_path, units):
    c = metadata(tmp_path, [dict(subject_id="donorA", **units), dict(subject_id="donorA", **units),
                            dict(subject_id="donorB", condition="treated")])
    with pytest.raises(ContractError, match="unique explicit"):
        validate_metadata(c)


def test_one_missing_culture_identity_cannot_infer_another_unit(tmp_path):
    c = metadata(tmp_path, [dict(subject_id="donorA", culture_id="cultureA"), dict(subject_id="donorA"),
                            dict(subject_id="donorB", condition="treated")])
    with pytest.raises(ContractError, match="unique explicit"):
        validate_metadata(c)


def test_unit_has_one_donor_owner_and_aliases_agree(tmp_path):
    c = metadata(tmp_path, [dict(subject_id="donorA", culture_id="shared"),
                            dict(subject_id="donorB", culture_id="shared", condition="treated")])
    with pytest.raises(ContractError, match="multiple subjects"):
        validate_metadata(c)
    c = metadata(tmp_path, [dict(experimental_unit_id="A", culture_id="B"), dict(condition="treated")])
    with pytest.raises(ContractError, match="must agree"):
        validate_metadata(c)


@pytest.mark.parametrize("inference", ["preliminary", "formal"])
def test_nested_cultures_reject_inference_even_with_enough_donors(tmp_path, inference):
    observations = [dict(subject_id=f"donor{donor}", experimental_unit_id=f"culture{donor}_{culture}", condition=condition)
                    for donor in range(3) for culture in range(2) for condition in ("control", "treated")]
    c = metadata(tmp_path, observations, inference=inference, design="~subject_id+condition")
    with pytest.raises(ContractError, match="nested donor/culture model"):
        validate_metadata(c)


def test_existing_paired_donors_and_explicit_unit_split_remain_supported(tmp_path):
    observations = [dict(subject_id=f"donor{donor}", culture_id=f"culture{donor}", condition=condition)
                    for donor in range(3) for condition in ("control", "treated")]
    result = validate_metadata(metadata(tmp_path, observations, inference="formal", design="~subject_id+condition"))
    assert len(result["samples"]) == 6
    assert len({s["subject_id"] for s in result["samples"]}) == 3


def test_independent_subject_minimum_still_applies(tmp_path):
    observations = [dict(subject_id=f"{condition}{donor}", condition=condition)
                    for donor in range(2) for condition in ("control", "treated")]
    c = metadata(tmp_path, observations, inference="formal")
    with pytest.raises(ContractError, match="3 distinct independent subjects"):
        validate_metadata(c)
