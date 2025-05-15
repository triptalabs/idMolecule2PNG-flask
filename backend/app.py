"""app.py - Flask API + CLI for mol2png project"""
from __future__ import annotations

import sys
from io import BytesIO
from pathlib import Path

from flask import Flask, abort, request, send_file
from rdkit import Chem
from rdkit.Chem import Draw

from backend.config import IMG_SIZE
from backend.resolver import ResolveError, fetch_smiles

app = Flask(__name__)


@app.get("/render")
def render():
    """HTTP endpoint that converts an identifier to a PNG image."""
    identifier: str = request.args.get("identifier", "").strip()
    if not identifier:
        abort(400, "Falta parámetro 'identifier'.")

    try:
        smiles = fetch_smiles(identifier)
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            abort(400, "SMILES inválido.")

        return _mol_to_response(mol)
    except ResolveError as e:
        abort(400, str(e))
    except Exception as e:  # pylint: disable=broad-except
        abort(500, f"Error interno: {e}")


# ---------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------

def _mol_to_image(mol: Chem.Mol):
    """Return a PIL.Image instance for a molecule."""
    return Draw.MolToImage(mol, size=IMG_SIZE)


def _mol_to_response(mol: Chem.Mol):
    """Return a Flask response with PNG image of a molecule."""
    img = _mol_to_image(mol)
    buf = BytesIO()
    img.save(buf, format="PNG")
    buf.seek(0)
    return send_file(buf, mimetype="image/png")


def save_png(identifier: str, outfile: str = "mol.png") -> Path:
    """Generate a PNG from an identifier and save it to *outfile*.

    This is used when the module is executed as a CLI:
        python -m backend.app "CHEMBL263881" lsd.png
    """
    smiles = fetch_smiles(identifier)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("SMILES inválido")

    img = _mol_to_image(mol)
    path = Path(outfile)
    img.save(path)
    print(
        f"✔  PNG guardado en {path.resolve()}  "
        f"[{identifier} → {smiles}]"
    )
    return path


# ---------------------------------------------------------------------
# Entrypoint
# ---------------------------------------------------------------------

def _run_cli() -> None:
    """CLI mode: convert identifier to PNG without starting the server."""
    if len(sys.argv) < 2:
        sys.exit(
            "Uso:\n"
            "  python -m backend.app IDENTIFICADOR [archivo.png]\n"
            "Ejemplos:\n"
            "  python -m backend.app 50-37-3 lsd.png\n"
            "  python -m backend.app CHEMBL263881\n"
        )

    identifier = sys.argv[1]
    outfile = sys.argv[2] if len(sys.argv) > 2 else "mol.png"
    try:
        save_png(identifier, outfile)
    except ResolveError as exc:
        sys.exit(f"✖  {exc}")


if __name__ == "__main__":
    # If arguments are supplied, run in CLI mode; otherwise start Flask.
    if len(sys.argv) == 1:
        app.run(host="0.0.0.0", port=5000, debug=True)
    else:
        _run_cli()
