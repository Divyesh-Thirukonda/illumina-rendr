import React, { useRef, useState } from 'react';

type Sphere = {
  x: number;
  y: number;
  z: number;
  radius: number;
  r: number;
  g: number;
  b: number;
};

function App() {
  const canvasRef = useRef<HTMLCanvasElement>(null);

  const [sphere, setSphere] = useState<Sphere>({
    x: 0.0,
    y: 0.0,
    z: 0.0,
    radius: 0.3,
    r: 255,
    g: 0,
    b: 0,
  });

  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const { name, value } = e.target;
    setSphere((prev) => ({ ...prev, [name]: parseFloat(value) }));
  };

  const addSphere = async () => {
    setLoading(true);
    setError(null);

    try {
      await fetch('http://localhost:3001/add_sphere', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(sphere),
      });

      const res = await fetch('http://localhost:3001/render');
      const buffer = new Uint8Array(await res.arrayBuffer());

      // --- Robust PPM header parsing ---
      const textChunk = new TextDecoder('ascii').decode(buffer.slice(0, 1024));
      const parts = textChunk.split(/\s+/).filter(Boolean);

      if (parts[0] !== 'P6') throw new Error('Not a P6 PPM file');
      const width = parseInt(parts[1], 10);
      const height = parseInt(parts[2], 10);
      const maxval = parseInt(parts[3], 10);
      if (maxval !== 255) throw new Error('Unsupported maxval');

      const headerEnd = textChunk.indexOf('255') + 3;
      const pixelStart = textChunk.indexOf('\n', headerEnd) + 1;

      const pixelData = buffer.slice(pixelStart);


      // --- Create RGBA buffer ---
      const rgba = new Uint8ClampedArray(width * height * 4);
      for (let i = 0, j = 0; i < pixelData.length; i += 3, j += 4) {
        rgba[j] = pixelData[i];       // R
        rgba[j + 1] = pixelData[i + 1]; // G
        rgba[j + 2] = pixelData[i + 2]; // B
        rgba[j + 3] = 255;             // A
      }

      // --- Draw to canvas ---
      const canvas = canvasRef.current;
      if (canvas) {
        canvas.width = width;
        canvas.height = height;
        const ctx = canvas.getContext('2d');
        if (ctx) {
          const imageData = new ImageData(rgba, width, height);
          ctx.putImageData(imageData, 0, 0);
        }
      }

      setLoading(false);
    } catch (err) {
      console.error('Error:', err);
      setError('Failed to render');
      setLoading(false);
    }
  };

  return (
    <div style={{ padding: '2rem', fontFamily: 'monospace' }}>
      <h1>Raycaster UI</h1>

      <div style={{ display: 'grid', gridTemplateColumns: 'auto auto', gap: '1rem', maxWidth: '400px' }}>
        {Object.entries(sphere).map(([key, value]) => (
          <label key={key}>
            {key}: <input type="number" step="0.1" name={key} value={value} onChange={handleChange} />
          </label>
        ))}
      </div>

      <button onClick={addSphere} style={{ marginTop: '1rem' }}>
        Add Sphere + Render
      </button>

      {loading && <p>Rendering...</p>}
      {error && <p style={{ color: 'red' }}>{error}</p>}

      <canvas ref={canvasRef} style={{ marginTop: '1rem', border: '1px solid #ccc' }} />
    </div>
  );
}

export default App;
