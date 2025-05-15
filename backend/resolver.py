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
    """Convierte CAS o ChEMBL a SMILES.  
    Si ya es SMILES lo devuelve tal cual."""
    if SMILES_RX.fullmatch(identifier):
        return identifier

    if CAS_RX.fullmatch(identifier):
        url = f"{PUBCHEM}/compound/name/{identifier}/property/IsomericSMILES/JSON"
        r = requests.get(url, timeout=HTTP_TIMEOUT)
        r.raise_for_status()
        return r.json()["PropertyTable"]["Properties"][0]["IsomericSMILES"]

    if identifier.lower().startswith("chembl"):
        r = requests.get(CHEMBL.format(identifier.upper()), timeout=HTTP_TIMEOUT)
        r.raise_for_status()
        return r.json()["molecule_structures"]["canonical_smiles"]

    raise ResolveError("Identificador no soportado")
