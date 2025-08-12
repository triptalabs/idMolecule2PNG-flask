import pathlib
import sys

import pytest
import requests

sys.path.append(str(pathlib.Path(__file__).resolve().parents[2]))
from backend.resolver import fetch_smiles, ResolveError


def test_pubchem_network_error(monkeypatch):
    """Se lanza ResolveError si la solicitud a PubChem falla."""

    def fake_get(*args, **kwargs):
        raise requests.RequestException("network down")

    monkeypatch.setattr("backend.resolver.requests.get", fake_get)

    with pytest.raises(ResolveError, match="PubChem"):
        fetch_smiles("50-00-0")  # CAS de formaldehído


def test_pubchem_invalid_response(monkeypatch):
    """Se lanza ResolveError si la respuesta de PubChem es incompleta."""

    class FakeResponse:
        def raise_for_status(self):
            pass

        def json(self):
            return {}

    monkeypatch.setattr("backend.resolver.requests.get", lambda *a, **k: FakeResponse())

    with pytest.raises(ResolveError, match="Respuesta de PubChem inválida"):
        fetch_smiles("50-00-0")


def test_chembl_network_error(monkeypatch):
    """Se lanza ResolveError si la solicitud a ChEMBL falla."""

    def fake_get(*args, **kwargs):
        raise requests.RequestException("boom")

    monkeypatch.setattr("backend.resolver.requests.get", fake_get)

    with pytest.raises(ResolveError, match="ChEMBL"):
        fetch_smiles("ChEMBL1")

