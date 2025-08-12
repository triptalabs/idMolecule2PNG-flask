import sys
from pathlib import Path

import pytest

# Ensure project root is on path so 'backend' package can be imported
sys.path.append(str(Path(__file__).resolve().parents[2]))
from backend.app import app


@pytest.fixture
def client():
    """Flask test client for the application."""
    app.config.update(TESTING=True)
    with app.test_client() as client:
        yield client


def test_render_returns_png(client):
    """Valid identifier should return PNG image."""
    res = client.get("/render", query_string={"identifier": "C"})
    assert res.status_code == 200
    assert res.mimetype == "image/png"
    assert len(res.data) > 0


def test_render_missing_identifier(client):
    """Missing identifier should return 400 with error message."""
    res = client.get("/render")
    assert res.status_code == 400
    assert "Falta parámetro" in res.get_data(as_text=True)


def test_render_invalid_identifier(client):
    """Invalid identifier should return 400 with error message."""
    res = client.get("/render", query_string={"identifier": "??"})
    assert res.status_code == 400
    assert "Identificador no soportado" in res.get_data(as_text=True)
