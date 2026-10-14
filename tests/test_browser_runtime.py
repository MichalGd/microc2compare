"""Optional prefix admission requires real converter APIs and smoke evidence."""
from pathlib import Path
import subprocess
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from microc2compare.browser_runtime import probe


def test_missing_converters_are_not_an_accepted_browser_prefix():
    with patch("microc2compare.browser_runtime.shutil.which", return_value=None):
        with pytest.raises(ValueError, match="requires bedGraphToBigWig and clodius"):
            probe()


@pytest.mark.parametrize("old_cli", [False, True], ids=["supported-cli", "unsupported-cli"])
def test_prefix_probe_checks_commands_and_decoded_smoke(tmp_path, old_cli):
    tools = {}
    for name in ("bedGraphToBigWig", "clodius"):
        path = tmp_path / name
        path.write_bytes(b"controlled converter executable fixture")
        tools[name] = str(path)
    records = {"probe_chr1": [(0, 10, 0.0), (20, 30, -2.5)], "probe_chr2": [(10, 20, 1.25)]}
    reader = SimpleNamespace(isBigWig=lambda: True, chroms=lambda: {"probe_chr2": 80, "probe_chr1": 100},
                             intervals=lambda chrom: records[chrom], close=lambda: None)
    commands = []
    def execute(command, **kwargs):
        commands.append(command)
        if "--help" in command:
            help_text = "old CLI" if old_cli else "--chromsizes-filename --output-file --delimiter --chr1-col --from1-col --to1-col --chr2-col --from2-col --to2-col"
            return subprocess.CompletedProcess(command, 0, stdout=help_text, stderr="")
        Path(command[-1]).write_bytes(b"\x26\xfc\x8f\x88" + b"controlled binary fixture")
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")
    with patch("microc2compare.browser_runtime.shutil.which", side_effect=tools.get), patch("microc2compare.browser_runtime.importlib.metadata.version", return_value="fixture-version"), patch("microc2compare.browser_runtime.subprocess.run", side_effect=execute), patch.dict("sys.modules", {"pyBigWig": SimpleNamespace(open=lambda path: reader)}):
        if old_cli:
            with pytest.raises(ValueError, match="lacks required exporter options"):
                probe()
            assert len(commands) == 1
        else:
            evidence = probe()
            assert evidence["status"] == "passed"
            assert evidence["bigwig_smoke"]["zero_and_negative_values_checked"]
            assert evidence["bigwig_smoke"]["chromosome_order"] == ["probe_chr2", "probe_chr1"]
            assert not evidence["higlass_ingested"] and not evidence["higlass_render_tested"]
            assert set(evidence["clodius_cli"]) == {"bedfile", "bedpe"}
            assert len(commands) == 3
