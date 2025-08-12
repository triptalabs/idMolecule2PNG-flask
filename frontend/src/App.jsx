import { useState, useEffect } from 'react'
import axios from 'axios'

export default function App() {
  const [identifier, setIdentifier] = useState('')
  const [filename, setFilename] = useState('mol.png')
  const [imgUrl, setImgUrl] = useState(null)
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState('')

  useEffect(() => {
    return () => {
      if (imgUrl) URL.revokeObjectURL(imgUrl)
    }
  }, [imgUrl])

  const renderMol = async (e) => {
    e.preventDefault()
    if (!identifier.trim()) return
    setLoading(true)
    setError('')
    setImgUrl(null)
    try {
      const { data } = await axios.get('/render', {
        params: { identifier: identifier.trim() },
        responseType: 'blob'
      })
      const url = URL.createObjectURL(data)
      setImgUrl(url)
    } catch (err) {
      if (!err.response) setError('Error de conexión')
      else setError(err.response.data || 'Error al renderizar')
    } finally {
      setLoading(false)
    }
  }

  const download = () => {
    if (!imgUrl) return
    const link = document.createElement('a')
    const safeName = filename.trim().replace(/[^a-z0-9_\-.]/gi, '_') || 'mol.png'
    link.href = imgUrl
    link.download = safeName
    document.body.appendChild(link)
    link.click()
    document.body.removeChild(link)
  }

  return (
    <div className="min-h-screen bg-gray-50 flex flex-col items-center p-6">
      <h1 className="text-3xl font-bold mb-6">mol2png v1</h1>

      <form onSubmit={renderMol} className="w-full max-w-md space-y-4">
        <div>
          <label htmlFor="identifier" className="block text-sm font-medium mb-1">
            Identificador (SMILES, CAS o ChEMBL)
          </label>
          <input
            id="identifier"
            type="text"
            value={identifier}
            onChange={e => setIdentifier(e.target.value)}
            className="w-full border rounded px-3 py-2"
            placeholder="p.ej. 50-37-3"
          />
        </div>

        <div>
          <label htmlFor="filename" className="block text-sm font-medium mb-1">
            Nombre de archivo PNG
          </label>
          <input
            id="filename"
            type="text"
            value={filename}
            onChange={e => setFilename(e.target.value)}
            className="w-full border rounded px-3 py-2"
            placeholder="mol.png"
          />
        </div>

        <button
          type="submit"
          disabled={loading}
          className="w-full bg-blue-600 text-white py-2 rounded hover:bg-blue-700 disabled:opacity-50"
        >
          {loading ? 'Renderizando...' : 'Renderizar y Mostrar'}
        </button>
      </form>

      {error && <p className="mt-4 text-red-600">{error}</p>}

      {imgUrl && (
        <div className="mt-6 flex flex-col items-center space-y-2">
          <div className="border p-2 bg-white">
            <img src={imgUrl} alt="molécula" className="max-w-xs" />
          </div>
          <button
            onClick={download}
            disabled={loading}
            className="bg-green-600 text-white px-4 py-2 rounded hover:bg-green-700 disabled:opacity-50"
          >
            Descargar PNG
          </button>
        </div>
      )}
    </div>
  )
}
