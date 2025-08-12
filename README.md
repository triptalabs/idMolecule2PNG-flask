# idMolecule2PNG-flask

Aplicación full‑stack que transforma identificadores químicos (CAS, ChEMBL o SMILES) en imágenes PNG de las moléculas correspondientes. El backend usa **Flask** y **RDKit** para resolver el identificador y dibujar la estructura; el frontend ofrece una interfaz React para consultar el servicio y descargar el resultado.

## Estructura del repositorio

| Ruta | Descripción |
|------|-------------|
| `backend/` | API de Flask y CLI para generar imágenes. |
| `frontend/` | Aplicación React + Vite que consume la API. |
| `render` | archivo temporal generado al probar la API (puede ignorarse). |

## Backend

El servicio HTTP expone un único endpoint `GET /render` que recibe el parámetro `identifier` y devuelve una imagen PNG. Internamente utiliza RDKit para generar la imagen y la función `fetch_smiles` para resolver los identificadores admitidos.

```python
@app.get("/render")
def render():
    identifier: str = request.args.get("identifier", "").strip()
    ...
    smiles = fetch_smiles(identifier)
    mol = Chem.MolFromSmiles(smiles)
    return _mol_to_response(mol)
```

La resolución admite tres formatos:

1. **CAS**: consulta a la API de PubChem.
2. **ChEMBL**: consulta a la API de ChEMBL.
3. **SMILES**: si no coincide con las anteriores se interpreta directamente como SMILES.

Estos pasos están implementados en `backend/resolver.py`.

El módulo también puede ejecutarse como CLI para generar un PNG sin iniciar el servidor:

```bash
python -m backend.app IDENTIFICADOR [archivo.png]
```

### Dependencias

Instalación recomendada (crear entorno virtual si se desea):

```bash
pip install -r backend/requirements.txt
```

## Frontend

Interfaz escrita en React, Vite y Tailwind CSS. Ofrece un formulario donde el usuario introduce el identificador, envía la consulta al endpoint `/render` y muestra la imagen resultante con opción de descarga.

### Comandos útiles

```bash
cd frontend
npm install       # instalar dependencias una sola vez
npm run dev       # ejecutar servidor de desarrollo
npm run build     # generar versión de producción en dist/
```

## Ejecución conjunta

1. Arrancar el backend (`python -m backend.app`).
2. Arrancar el frontend (`npm run dev`) y abrir el URL que indica Vite.
3. Introducir un identificador químico (p.ej. `50-37-3`) y obtener la representación PNG.

## Pruebas

El directorio `backend/tests/` contiene la estructura para pruebas con `pytest`. Para ejecutarlas:

```bash
pytest
```

El frontend incluye reglas de linting ejecutables con:

```bash
cd frontend
npm run lint
```

## Licencia

Este proyecto se ofrece sin licencia explícita; verifique antes de usar en producción.
