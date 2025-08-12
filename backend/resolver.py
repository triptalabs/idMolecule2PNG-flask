import re
import requests
from backend.config import HTTP_TIMEOUT

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CHEMBL  = "https://www.ebi.ac.uk/chembl/api/data/molecule/{}?format=json"

CAS_RX     = re.compile(r"^\d{2,7}-\d{2}-\d$")
SMILES_RX  = re.compile(r"^[A-Za-z0-9@+\-\[\]\(\)\\\/%=#$]+$")

class ResolveError(Exception):
    """Error de resolución de identificador químico."""

def fetch_smiles(identifier: str) -> str:
    """Convierte CAS o ChEMBL a SMILES; si nada coincide, asume SMILES."""
    # 1️⃣ CAS
    if CAS_RX.fullmatch(identifier):
        url = f"{PUBCHEM}/compound/name/{identifier}/property/IsomericSMILES/JSON"
        try:
            r = requests.get(url, timeout=HTTP_TIMEOUT)
            r.raise_for_status()
        except requests.RequestException as exc:
            raise ResolveError(f"Error al consultar PubChem: {exc}") from exc
        data = r.json()
        smiles = data.get("PropertyTable", {}).get("Properties", [{}])[0].get("IsomericSMILES")
        if not smiles:
            raise ResolveError("Respuesta de PubChem inválida")
        return smiles

    # 2️⃣ ChEMBL
    if identifier.lower().startswith("chembl"):
        try:
            r = requests.get(CHEMBL.format(identifier.upper()), timeout=HTTP_TIMEOUT)
            r.raise_for_status()
        except requests.RequestException as exc:
            raise ResolveError(f"Error al consultar ChEMBL: {exc}") from exc
        data = r.json()
        smiles = data.get("molecule_structures", {}).get("canonical_smiles")
        if not smiles:
            raise ResolveError("Respuesta de ChEMBL inválida")
        return smiles

    # 3️⃣ SMILES (por descarte)
    if SMILES_RX.fullmatch(identifier):
        return identifier

    raise ResolveError("Identificador no soportado")
