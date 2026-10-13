"""Generic project initialization and common measurement-universe admission."""
import importlib.util
from pathlib import Path

import pytest

from microc2compare.contracts import ContractError, normalize, read_json


def environment(config):
    return {key: str(value).lower() if isinstance(value, bool) else str(value)
            for key, value in config.items()}


@pytest.mark.parametrize("bound", ["minimum", "maximum"])
def test_independent_discovery_cannot_bypass_the_measured_distance_universe(config, tmp_path, bound):
    values = environment(config)
    values['RUN_CHROMOSIGHT_DISCOVERY'] = 'true'
    key = 'CHROMOSIGHT_DISCOVERY_MIN_DISTANCE_BP' if bound == 'minimum' else 'CHROMOSIGHT_DISCOVERY_MAX_DISTANCE_BP'
    values[key] = str(config['LOOP_MIN_DISTANCE_BP'] - 1 if bound == 'minimum' else config['LOOP_MAX_DISTANCE_BP'] + 1)
    with pytest.raises(ContractError, match='common .* measurement range'):
        normalize(values, tmp_path)


def test_raw_only_subkilobase_quantification_needs_explicit_opt_in(config, tmp_path):
    values = environment(config)
    for key in ('RUN_LOOPS', 'RUN_CHROMOSIGHT_DISCOVERY', 'RUN_CHROMOSIGHT', 'RUN_ANNOTATION',
                'RUN_AGGREGATES', 'RUN_COMPARTMENTS', 'RUN_INSULATION', 'RUN_DIFFERENTIAL'):
        values[key] = 'false'
    values.update(M2C_ANALYSES='quantify', INDEPENDENT_FEATURES_REGISTRY='independent-features.json',
                  LOOP_RESOLUTION='400', LOOP_QUANTIFICATION_RESOLUTIONS='400',
                  ALLOW_SUBKILOBASE_ANALYSIS='false')
    with pytest.raises(ContractError, match='Sub-kilobase analysis requires explicit'):
        normalize(values, tmp_path)


@pytest.mark.parametrize('preset', ['existing', 'hsieh_8grid', 'short_range_4kb'])
@pytest.mark.parametrize('assembly', ['hg38', 'newAssemblyV2'])
def test_new_project_preset_has_no_historical_sample_or_reference_dependency(tmp_path, preset, assembly):
    script = Path(__file__).resolve().parents[1]/'scripts/project_tools.py'
    spec = importlib.util.spec_from_file_location('fresh_project_tools', script)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    destination = tmp_path/'new-project/project.conf'
    module.initialize(destination, preset, assembly)
    text = destination.read_text()
    assert f'REFERENCE_ID="{assembly}"' in text and f'use_analysis_profile {preset}' in text
    assert 'RUN_BROWSER_EXPORT="true"' in text
    registry = read_json(destination.parent/'references.json')
    assert registry['reference_id'] == assembly and registry['chromsizes'] == ''
    assert 'basal_PMKs' not in text and 'mESCs' not in text
    with pytest.raises(ContractError, match='exists'):
        module.initialize(destination, preset, assembly)
