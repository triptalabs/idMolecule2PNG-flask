from io import BytesIO
from flask import Flask, request, send_file, abort
from rdkit import Chem
from rdkit.Chem import Draw
from backend.resolver import fetch_smiles, ResolveError
from backend.config import IMG_SIZE

app = Flask(__name__)

@app.get("/render")
def render():
    identifier = request.args.get("identifier", "").strip()
    if not identifier:
        abort(400, "Falta parámetro 'identifier'.")

    try:
        smiles = fetch_smiles(identifier)
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            abort(400, "SMILES inválido.")

        img = Draw.MolToImage(mol, size=IMG_SIZE)
        buf = BytesIO()
        img.save(buf, format="PNG")
        buf.seek(0)
        return send_file(buf, mimetype="image/png")
    except ResolveError as e:
        abort(400, str(e))
    except Exception as e:
        abort(500, f"Error interno: {e}")

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000, debug=True)
