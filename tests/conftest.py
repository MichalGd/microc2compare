import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))

import pytest
from microc2compare.contracts import read_json


@pytest.fixture
def config():
    schema=read_json(Path(__file__).resolve().parents[1]/"schemas/config.schema.json")
    return {k:v["default"] for k,v in schema["properties"].items()}
