import React, { useState, useRef } from 'react';

function App() {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const [sceneJSON, setSceneJSON] = useState<string>('{\n  "width": 400,\n  "height": 300,\n  "spheres": []\n}');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const handleRender = async () => {
    setLoading(true);
    setError(null);

    try {
      // Validate JSON before sending
      const parsed = JSON.parse(sceneJSON);

      // Send JSON to server
      await fetch('http://localhost:3001/scene', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(parsed),
      });

      // Request rendered image
      const res = await fetch('http://localhost:3001/render');
const text = await res.text();

// --- Parse P3 header ---
const lines = text.trim().split(/\s+/);
if (lines[0] !== 'P3') throw new Error('Not a P3 PPM file');

const width = parseInt(lines[1]);
const height = parseInt(lines[2]);
const maxVal = parseInt(lines[3]);

const pixelValues = lines.slice(4).map(Number);
if (pixelValues.length !== width * height * 3) throw new Error('Invalid pixel data length');

const rgba = new Uint8ClampedArray(width * height * 4);
for (let i = 0, j = 0; i < pixelValues.length; i += 3, j += 4) {
  rgba[j] = (pixelValues[i] / maxVal) * 255;
  rgba[j + 1] = (pixelValues[i + 1] / maxVal) * 255;
  rgba[j + 2] = (pixelValues[i + 2] / maxVal) * 255;
  rgba[j + 3] = 255;
}

// Draw to canvas
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
      <h1>Raytracer Scene Renderer</h1>
      <textarea
        value={sceneJSON}
        onChange={(e) => setSceneJSON(e.target.value)}
        rows={20}
        cols={80}
        style={{ fontFamily: 'monospace', fontSize: '14px', width: '100%' }}
      />
      <br />
      <button onClick={handleRender} style={{ marginTop: '1rem' }}>
        Render Scene
      </button>
      {loading && <p>Rendering...</p>}
      {error && <p style={{ color: 'red' }}>{error}</p>}
      <canvas ref={canvasRef} style={{ marginTop: '2rem', border: '1px solid #ccc' }} />
    </div>
  );
}

export default App;
