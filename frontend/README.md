# Frontend de idMolecule2PNG

Aplicación React que permite convertir identificadores químicos en imágenes PNG consumiendo la API del backend.

## Desarrollo

```bash
npm install       # instalar dependencias
npm run dev       # servidor de desarrollo en modo HMR
```

El formulario envía una petición `GET /render` con el parámetro `identifier` y muestra la imagen devuelta. También ofrece la opción de descargar el archivo PNG con un nombre definido por el usuario.

## Construcción

```bash
npm run build     # genera artefactos en dist/
```

## Lint

```bash
npm run lint
```
