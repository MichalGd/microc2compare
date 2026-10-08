"""Compartment orientation and export, independent of the optional eigensolver."""
import json
import sys
from multiprocessing.pool import ThreadPool
from types import SimpleNamespace

import numpy as np
import pandas as pd

from microc2compare import science


def test_compartment_export_preserves_orientation_missing_and_ambiguous_bins(tmp_path, monkeypatch):
    bins = pd.DataFrame({"chrom": ["chr1"] * 8,
                         "start": np.arange(8) * 100,
                         "end": (np.arange(8) + 1) * 100})
    views = pd.DataFrame({"chrom": ["chr1", "chr1"], "start": [0, 400],
                          "end": [400, 800], "name": ["left", "right"]})
    phase = bins.assign(value=[0, 1, 2, 3] * 2)
    phase_path = tmp_path / "phase.bedgraph"
    phase.to_csv(phase_path, sep="\t", header=False, index=False)
    vectors = bins.assign(E1=[1, -1, -1, 1] * 2,
                          E2=[3, 2, np.nan, 0, 1, 1, 1, 1], E3=1.0)
    eigenvalues = pd.DataFrame({"region": ["left", "right"], "eigval1": [2.0, 2.0]})
    clr = SimpleNamespace(extent=lambda region: (region[1] // 100, region[2] // 100))
    monkeypatch.setattr(science, "context", lambda *args: (clr, bins, None, None, views, None))
    monkeypatch.setattr(science, "reference", lambda c: ({"phasing_track": str(phase_path)}, None, None))

    def eigs_cis(actual_clr, **kwargs):
        assert actual_clr is clr
        pd.testing.assert_frame_equal(kwargs["view_df"], views)
        assert kwargs["clr_weight_name"] == "weight"
        assert kwargs["n_eigs"] == 3
        return eigenvalues.copy(), vectors.copy()

    monkeypatch.setitem(sys.modules, "cooltools", SimpleNamespace(eigs_cis=eigs_cis))
    # The stub eigensolver does not dispatch tasks; no child processes are needed.
    import multiprocess
    monkeypatch.setattr(multiprocess, "Pool", ThreadPool)
    c = {"run_dir": str(tmp_path / "results"), "COMPARTMENT_RESOLUTION": 100,
         "MAX_DENSE_BINS": 10, "THREADS_COMPARTMENTS": 1,
         "MEMORY_MB_COMPARTMENTS": 1, "INSULATION_IGNORE_DIAGS": 2,
         "MIN_PHASING_CORRELATION": 0.1}
    outputs = science.compartments(c, "C1")
    assert len(outputs) == 4 and all(path.is_file() for path in outputs)
    table = pd.read_csv(outputs[0], sep="\t")
    np.testing.assert_allclose(table.oriented_score[:4], [-3, -2, np.nan, 0], equal_nan=True)
    assert table.selected_component[:4].eq("E2").all()
    assert table.oriented_score[4:].isna().all()
    assert table.selected_component[4:].isna().all()
    assert "value" not in table
    qc = json.loads(outputs[2].read_text())
    assert qc[0]["status"] == "oriented" and qc[0]["component"] == "E2"
    assert qc[0]["correlation"] < -0.99 and qc[1]["status"] == "ambiguous"
    browser = pd.read_csv(outputs[3], sep="\t", header=None)
    assert list(browser.iloc[:, 1]) == [0, 100, 300]
    assert list(browser.iloc[:, 3]) == [-3, -2, 0]
    pd.testing.assert_frame_equal(pd.read_csv(outputs[1], sep="\t"), eigenvalues)
