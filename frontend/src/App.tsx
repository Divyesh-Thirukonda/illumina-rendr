import React, { useState, useRef, useEffect } from 'react';

import './App.css';

type Vec3 = [number, number, number];

const defaultScene = {
  eye: [0, 0, 5] as Vec3,
  viewDir: [0, 0, -1] as Vec3,
  upDir: [0, 1, 0] as Vec3,
  vfov: 45,
  width: 800,
  height: 600,
  bkgcolor: {
    rgb: [0.5, 0.5, 0.5] as Vec3,
    eta: 1.0
  },
  parallel: false,
  frustumHeight: 2.0,
  depthCueingEnabled: false,
  depthCueing: {
    color: [0.5, 0.5, 0.5] as Vec3,
    alphaMin: 0.1,
    alphaMax: 1.0,
    distMin: 0,
    distMax: 10
  },
  spheres: [
  {
    center: [0, 0, 0] as Vec3,
    radius: 1,
    diffuse: [0.8, 0.1, 0.1] as Vec3,
    specular: [1.0, 1.0, 1.0] as Vec3,
    intensity: [0.1, 0.8, 0.5] as Vec3,
    shininess: 32,
    alpha: 1.0,
    eta: 1.0,
    hasAlphaAt: false
  }
  ],
  lights: [
  {
    position: [5, 10, 5] as Vec3,
    w: 1,
    intensity: 1.0,
    useAttenuation: false,
    attenuation: [1, 0, 0] as Vec3
  }
  ],
  triangles: [],
  cubes: [
  {
    position: [0, -2.5, 0],
    scale: 2,
    rotation: [0, 45, 0], // in degrees
    material: {
      diffuse: [0.5, 0.8, 0.2],
      specular: [1, 1, 1],
      intensity: [0.1, 0.3, 0.6],
      shininess: 32,
      alpha: 1.0,
      eta: 1.0,
      hasAlphaAt: false
    }
  }
  ]

};

function App() {
  const canvasRef = useRef<HTMLCanvasElement>(null);

  const [scene, setScene] = useState(defaultScene);
  const [showSceneSettings, setShowSceneSettings] = useState(true);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const updateScene = (key: string, value: any) => {
    setScene((prev) => ({ ...prev, [key]: value }));
  };

  const updateVec3 = (key: keyof typeof defaultScene, index: number, value: number) => {
    setScene((prev) => ({
      ...prev,
      [key]: [...(prev[key] as Vec3).slice(0, index), value, ...(prev[key] as Vec3).slice(index + 1)]
    }));
  };

  const updateNestedVec3 = (
    outer: 'bkgcolor' | 'depthCueing',
    key: string,
    index: number,
    value: number
  ) => {
    setScene((prev) => ({
      ...prev,
      [outer]: {
        ...prev[outer],
        [key]: [
          ...(prev[outer] as any)[key].slice(0, index),
          value,
          ...(prev[outer] as any)[key].slice(index + 1)
        ]
      }
    }));
  };


  const cubeToTriangles = (cube: {
    position: number[];
    scale: number;
    rotation: number[];
    material: {
      diffuse: number[];
      specular: number[];
      intensity: number[];
      shininess: number;
      alpha: number;
      eta: number;
      hasAlphaAt: boolean;
      alphaR?: number;
      alphaG?: number;
      alphaB?: number;
    };
  }) => {
    const [cx, cy, cz] = cube.position;
    const [rx, ry, rz] = cube.rotation.map(r => r * Math.PI / 180); // radians
    const s = cube.scale / 2;

    let verts: [number, number, number][] = [
      [-s, -s, -s], [s, -s, -s], [s, s, -s], [-s, s, -s],
      [-s, -s, s], [s, -s, s], [s, s, s], [-s, s, s],
    ];

    // Rotation matrices
    const rotate = ([x, y, z]: [number, number, number]): [number, number, number] => {
      let [sx, cx_] = [Math.sin(rx), Math.cos(rx)];
      let [sy, cy_] = [Math.sin(ry), Math.cos(ry)];
      let [sz, cz_] = [Math.sin(rz), Math.cos(rz)];

      // rotate x
      let y1 = y * cx_ - z * sx;
      let z1 = y * sx + z * cx_;
      y = y1; z = z1;

      // rotate y
      let x1 = x * cy_ + z * sy;
      z1 = -x * sy + z * cy_;
      x = x1; z = z1;

      // rotate z
      x1 = x * cz_ - y * sz;
      y1 = x * sz + y * cz_;
      return [x1 + cx, y1 + cy, z + cz];
    };

    verts = verts.map(rotate);

    const f = [
      [0, 1, 2, 3], [5, 4, 7, 6],
      [4, 0, 3, 7], [1, 5, 6, 2],
      [3, 2, 6, 7], [4, 5, 1, 0],
    ];

    const m = cube.material;
    const tri = (a: number, b: number, c: number) => ({
      v0: verts[a], v1: verts[b], v2: verts[c],
      vt0: [0, 0], vt1: [1, 0], vt2: [0, 1],
      diffuse: m.diffuse,
      specular: m.specular,
      intensity: m.intensity,
      shininess: m.shininess,
      alpha: m.alpha,
      eta: m.eta,
      hasAlphaAt: m.hasAlphaAt,
      ...(m.hasAlphaAt ? {
        alphaR: m.alphaR ?? m.alpha,
        alphaG: m.alphaG ?? m.alpha,
        alphaB: m.alphaB ?? m.alpha,
      } : {})
    });

    return f.flatMap(([a, b, c, d]) => [tri(a, b, c), tri(a, c, d)]);
  };



  const handleRender = async () => {
    setLoading(true);
    setError(null);

    try {
      // Assemble JSON for backend
      const sceneToSend: any = {
        width: scene.width,
        height: scene.height,
        eye: scene.eye,
        viewDir: scene.viewDir,
        upDir: scene.upDir,
        vfov: scene.vfov,
        bkgcolor: {
          rgb: scene.bkgcolor.rgb,
          eta: scene.bkgcolor.eta
        },
        parallel: scene.parallel,
        ...(scene.parallel && { frustumHeight: scene.frustumHeight }),
        ...(scene.depthCueingEnabled && {
          depthCueing: {
            color: scene.depthCueing.color,
            alphaMin: scene.depthCueing.alphaMin,
            alphaMax: scene.depthCueing.alphaMax,
            distMin: scene.depthCueing.distMin,
            distMax: scene.depthCueing.distMax
          }
        }),
        spheres: scene.spheres.map(sphere => ({
          center: sphere.center,
          radius: sphere.radius,
          diffuse: sphere.diffuse,
          specular: sphere.specular,
          intensity: sphere.intensity,
          shininess: sphere.shininess,
          alpha: sphere.alpha,
          eta: sphere.eta,
          hasAlphaAt: sphere.hasAlphaAt
        })),
        lights: scene.lights.map(light => {
          const base = {
            position: light.position,
            w: light.w,
            intensity: light.intensity
          };
          if (light.useAttenuation) {
            return { ...base, attenuation: light.attenuation };
          }
          return base;
        }),
        triangles: [
          ...scene.triangles,
          ...scene.cubes.flatMap(cubeToTriangles)
        ]


      };

      const BACKEND_URL = process.env.REACT_APP_BACKEND_URL || 'http://localhost:3001';

      // POST scene to backend
      const res = await fetch(BACKEND_URL + '/render', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(sceneToSend),
      });

      const text = await res.text();

      // Parse P3 PPM
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


  useEffect(() => {
    handleRender();
  }, [scene]);


  // all spheres visibility state should be true by default
  const [visibleSpheres, setVisibleSpheres] = useState(() => {
    const initialVisibility: { [key: number]: boolean } = {};
    for (let i = 0; i < scene.spheres.length; i++) {
      initialVisibility[i] = false;
    }
    return initialVisibility;
  });
  const toggleSphereVisibility = (idx: number) => {
    setVisibleSpheres(prev => ({
      ...prev,
      [idx]: !prev[idx],
    }));
  };

  const [visibleLights, setVisibleLights] = useState(() => {
    const initialVisibility: { [key: number]: boolean } = {};
    for (let i = 0; i < scene.lights.length; i++) {
      initialVisibility[i] = false;
    }
    return initialVisibility;
  });
  const toggleLightVisibility = (idx: number) => {
    setVisibleLights(prev => ({
      ...prev,
      [idx]: !prev[idx],
    }));
  };

  const [visibleCubes, setVisibleCubes] = useState(() => {
    const initialVisibility: { [key: number]: boolean } = {};
    for (let i = 0; i < scene.cubes.length; i++) {
      initialVisibility[i] = false;
    }
    return initialVisibility;
  });
  const toggleCubeVisibility = (idx: number) => {
    setVisibleCubes(prev => ({
      ...prev,
      [idx]: !prev[idx],
    }));
  };












  return (
    
  <div className="app-container">
    <div className="canvas-container">
      <canvas ref={canvasRef} style={{ border: '1px solid #ccc' }} />
    </div>
    <div className="settings-container">

      <button onClick={() => setShowSceneSettings(!showSceneSettings)} style={{ marginBottom: '1rem' }}>
        {showSceneSettings ? 'Hide' : 'Show'} Scene Settings
      </button>

      {showSceneSettings && (
        <div style={{ border: '1px solid #ccc', padding: '1rem', marginBottom: '1rem' }}>
          <h3>Camera</h3>
          <label>eye: {scene.eye.map((v, i) => (
            <input key={i} type="number" value={v} onChange={(e) => updateVec3("eye", i, parseFloat(e.target.value))} />
          ))}</label><br />
          <label>viewDir: {scene.viewDir.map((v, i) => (
            <input key={i} type="number" value={v} onChange={(e) => updateVec3("viewDir", i, parseFloat(e.target.value))} />
          ))}</label><br />
          <label>upDir: {scene.upDir.map((v, i) => (
            <input key={i} type="number" value={v} onChange={(e) => updateVec3("upDir", i, parseFloat(e.target.value))} />
          ))}</label><br />
          <label>vfov: <input type="number" value={scene.vfov} onChange={(e) => updateScene("vfov", parseFloat(e.target.value))} /></label><br />
          <label>width: <input type="number" value={scene.width} onChange={(e) => updateScene("width", parseFloat(e.target.value))} /></label><br />
          <label>height: <input type="number" value={scene.height} onChange={(e) => updateScene("height", parseFloat(e.target.value))} /></label><br />

          <h3>Background</h3>
          <label>rgb: {scene.bkgcolor.rgb.map((v, i) => (
            <input key={i} type="number" step="0.01" value={v} onChange={(e) =>
              updateNestedVec3("bkgcolor", "rgb", i, parseFloat(e.target.value))
            } />
          ))}</label><br />
          <label>eta: <input type="number" value={scene.bkgcolor.eta}
            onChange={(e) => updateScene("bkgcolor", { ...scene.bkgcolor, eta: parseFloat(e.target.value) })} /></label><br />

          <h3>Projection</h3>
          <label>
            <input type="checkbox" checked={scene.parallel} onChange={() => updateScene("parallel", !scene.parallel)} />
            Parallel Projection
          </label><br />
          {scene.parallel && (
            <label>frustumHeight: <input type="number" value={scene.frustumHeight}
              onChange={(e) => updateScene("frustumHeight", parseFloat(e.target.value))} /></label>
          )}

          <h3>Depth Cueing</h3>
          <label>
            <input type="checkbox" checked={scene.depthCueingEnabled}
              onChange={() => updateScene("depthCueingEnabled", !scene.depthCueingEnabled)} />
            Enable Depth Cueing
          </label><br />
          {scene.depthCueingEnabled && (
            <>
              <label>color: {scene.depthCueing.color.map((v, i) => (
                <input key={i} type="number" value={v} onChange={(e) =>
                  updateScene("depthCueing", {
                    ...scene.depthCueing,
                    color: [...scene.depthCueing.color.slice(0, i), parseFloat(e.target.value), ...scene.depthCueing.color.slice(i + 1)]
                  })} />
              ))}</label><br />
              <label>alphaMin: <input type="number" value={scene.depthCueing.alphaMin}
                onChange={(e) => updateScene("depthCueing", { ...scene.depthCueing, alphaMin: parseFloat(e.target.value) })} /></label><br />
              <label>alphaMax: <input type="number" value={scene.depthCueing.alphaMax}
                onChange={(e) => updateScene("depthCueing", { ...scene.depthCueing, alphaMax: parseFloat(e.target.value) })} /></label><br />
              <label>distMin: <input type="number" value={scene.depthCueing.distMin}
                onChange={(e) => updateScene("depthCueing", { ...scene.depthCueing, distMin: parseFloat(e.target.value) })} /></label><br />
              <label>distMax: <input type="number" value={scene.depthCueing.distMax}
                onChange={(e) => updateScene("depthCueing", { ...scene.depthCueing, distMax: parseFloat(e.target.value) })} /></label><br />
            </>
          )}
        </div>
      )}



      <h2>Spheres</h2>
      {scene.spheres.map((s, idx) => (
        <div key={idx} style={{ border: '1px solid #ccc', padding: '1rem', marginBottom: '1rem' }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
            <h4>Sphere {idx + 1}</h4>
            <button
              onClick={() => toggleSphereVisibility(idx)}
              style={{ color: 'gray' }}
            >
              {visibleSpheres[idx] === false ? '+' : '-'}
            </button>
          </div>

          {visibleSpheres[idx] !== false && (
            <>
              <label>Center: {s.center.map((v, i) => (
                <input key={i} type="number" value={v}
                  onChange={(e) => {
                    const updated = [...scene.spheres];
                    updated[idx].center[i] = parseFloat(e.target.value);
                    updateScene("spheres", updated);
                  }} />
              ))}</label><br />
              <label>Radius: <input type="number" value={s.radius}
                onChange={(e) => {
                  const updated = [...scene.spheres];
                  updated[idx].radius = parseFloat(e.target.value);
                  updateScene("spheres", updated);
                }} /></label><br />
              <label>Diffuse: {s.diffuse.map((v, i) => (
                <input key={i} type="number" value={v}
                  onChange={(e) => {
                    const updated = [...scene.spheres];
                    updated[idx].diffuse[i] = parseFloat(e.target.value);
                    updateScene("spheres", updated);
                  }} />
              ))}</label><br />
              <label>Specular: {s.specular.map((v, i) => (
                <input key={i} type="number" value={v}
                  onChange={(e) => {
                    const updated = [...scene.spheres];
                    updated[idx].specular[i] = parseFloat(e.target.value);
                    updateScene("spheres", updated);
                  }} />
              ))}</label><br />
              <label>Intensity: {s.intensity.map((v, i) => (
                <input key={i} type="number" value={v}
                  onChange={(e) => {
                    const updated = [...scene.spheres];
                    updated[idx].intensity[i] = parseFloat(e.target.value);
                    updateScene("spheres", updated);
                  }} />
              ))}</label><br />
              <label>Shininess: <input type="number" value={s.shininess}
                onChange={(e) => {
                  const updated = [...scene.spheres];
                  updated[idx].shininess = parseFloat(e.target.value);
                  updateScene("spheres", updated);
                }} /></label><br />
              <label>Alpha: <input type="number" value={s.alpha}
                onChange={(e) => {
                  const updated = [...scene.spheres];
                  updated[idx].alpha = parseFloat(e.target.value);
                  updateScene("spheres", updated);
                }} /></label><br />
              <label>Eta: <input type="number" value={s.eta}
                onChange={(e) => {
                  const updated = [...scene.spheres];
                  updated[idx].eta = parseFloat(e.target.value);
                  updateScene("spheres", updated);
                }} /></label><br />
            </>
          )}

          <button onClick={() => {
            const updated = scene.spheres.filter((_, i) => i !== idx);
            updateScene("spheres", updated);
            // make sure to update visibility state
            setVisibleSpheres(prev => {
              const newVisibility = { ...prev };
              delete newVisibility[idx];
              return newVisibility;
            });
          }} style={{ marginTop: '0.5rem', color: 'red' }}>Delete Sphere</button>
        </div>
      ))}


      <button onClick={() => {
        const newSphere = {
          center: [0, 0, 0] as Vec3,
          radius: 1,
          diffuse: [0.8, 0.1, 0.1] as Vec3,
          specular: [1, 1, 1] as Vec3,
          intensity: [0.1, 0.8, 0.5] as Vec3,
          shininess: 32,
          alpha: 1.0,
          eta: 1.0,
          hasAlphaAt: false
        };
        updateScene("spheres", [...scene.spheres, newSphere]);
        // make sure to update visibility state
        setVisibleSpheres(prev => ({
          ...prev,
          [scene.spheres.length]: true
        }));
      }} style={{ marginBottom: '2rem' }}>+ Add Sphere</button>

      {/* ---------- Lights ---------- */}
      <h2>Lights</h2>
      {scene.lights.map((light, idx) => (
        <div key={idx} style={{ border: '1px solid #ccc', padding: '1rem', marginBottom: '1rem' }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
            <h4>Light {idx + 1}</h4>
            <button
              onClick={() => toggleLightVisibility(idx)}
              style={{ color: 'gray' }}
            >
              {visibleLights[idx] === false ? '+' : '-'}
            </button>
          </div>

          {visibleLights[idx] !== false && (
            <>
              <label>Position: {light.position.map((v, i) => (
                <input key={i} type="number" value={v}
                  onChange={(e) => {
                    const updated = [...scene.lights];
                    updated[idx].position[i] = parseFloat(e.target.value);
                    updateScene("lights", updated);
                  }} />
              ))}</label><br />

              <label>w: <input type="number" value={light.w}
                onChange={(e) => {
                  const updated = [...scene.lights];
                  updated[idx].w = parseFloat(e.target.value);
                  updateScene("lights", updated);
                }} /></label><br />

              <label>intensity: <input type="number" value={light.intensity}
                onChange={(e) => {
                  const updated = [...scene.lights];
                  updated[idx].intensity = parseFloat(e.target.value);
                  updateScene("lights", updated);
                }} /></label><br />

              <label>
                <input type="checkbox" checked={light.useAttenuation}
                  onChange={() => {
                    const updated = [...scene.lights];
                    updated[idx].useAttenuation = !updated[idx].useAttenuation;
                    updateScene("lights", updated);
                  }} />
                Use Attenuation
              </label><br />

              {light.useAttenuation && (
                <label>Attenuation: {light.attenuation.map((v, i) => (
                  <input key={i} type="number" value={v}
                    onChange={(e) => {
                      const updated = [...scene.lights];
                      updated[idx].attenuation[i] = parseFloat(e.target.value);
                      updateScene("lights", updated);
                    }} />
                ))}</label>
              )}<br />
            </>
          )}

          <button onClick={() => {
            const updated = scene.lights.filter((_, i) => i !== idx);
                      // make sure to update visibility state
                      setVisibleSpheres(prev => {
            const newVisibility = { ...prev };
            delete newVisibility[idx];
            return newVisibility;
          });
            updateScene("lights", updated);
          }} style={{ marginTop: '0.5rem', color: 'red' }}>Delete Light</button>
        </div>
      ))}

      <button onClick={() => {
        const newLight = {
          position: [0, 10, 0] as Vec3,
          w: 1,
          intensity: 1.0,
          useAttenuation: false,
          attenuation: [1, 0, 0] as Vec3
        };
        // make sure to update visibility state
        setVisibleLights(prev => ({
          ...prev,
          [scene.lights.length]: true
        }));
        updateScene("lights", [...scene.lights, newLight]);
      }} style={{ marginBottom: '2rem' }}>+ Add Light</button>




      {/* ---------- Cubes ---------- */}
      <h2>Cubes</h2>
      {scene.cubes.map((cube, idx) => (
        <div key={idx} style={{ border: '1px solid #ccc', padding: '1rem', marginBottom: '1rem' }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
            <h4>Cube {idx + 1}</h4>
            <button
              onClick={() => toggleCubeVisibility(idx)}
              style={{ color: 'gray' }}
            >
              {visibleCubes[idx] === false ? '+' : '-'}
            </button>
          </div>

          {visibleCubes[idx] !== false && (
            <>
              <label>Position: {cube.position.map((v, i) => (
                <input key={i} type="number" value={v}
                  onChange={(e) => {
                    const updated = [...scene.cubes];
                    updated[idx].position[i] = parseFloat(e.target.value);
                    updateScene("cubes", updated);
                  }} />
              ))}</label><br />

              <label>Scale: <input type="number" value={cube.scale}
                onChange={(e) => {
                  const updated = [...scene.cubes];
                  updated[idx].scale = parseFloat(e.target.value);
                  updateScene("cubes", updated);
                }} /></label><br />

              <label>Rotation (degrees): {cube.rotation.map((v, i) => (
                <input key={i} type="number" value={v}
                  onChange={(e) => {
                    const updated = [...scene.cubes];
                    updated[idx].rotation[i] = parseFloat(e.target.value);
                    updateScene("cubes", updated);
                  }} />
              ))}</label><br />

              <h5>Material</h5>
              <label>Diffuse: {cube.material.diffuse.map((v, i) => (
                <input key={i} type="number" value={v}
                  onChange={(e) => {
                    const updated = [...scene.cubes];
                    updated[idx].material.diffuse[i] = parseFloat(e.target.value);
                    updateScene("cubes", updated);
                  }} />
              ))}</label><br />

              <label>Specular: {cube.material.specular.map((v, i) => (
                <input key={i} type="number" value={v}
                  onChange={(e) => {
                    const updated = [...scene.cubes];
                    updated[idx].material.specular[i] = parseFloat(e.target.value);
                    updateScene("cubes", updated);
                  }} />
              ))}</label><br />

              <label>Intensity: {cube.material.intensity.map((v, i) => (
                <input key={i} type="number" value={v}
                  onChange={(e) => {
                    const updated = [...scene.cubes];
                    updated[idx].material.intensity[i] = parseFloat(e.target.value);
                    updateScene("cubes", updated);
                  }} />
              ))}</label><br />

              <label>Shininess: <input type="number" value={cube.material.shininess}
                onChange={(e) => {
                  const updated = [...scene.cubes];
                  updated[idx].material.shininess = parseFloat(e.target.value);
                  updateScene("cubes", updated);
                }} /></label><br />
            </>
          )}

          <button onClick={() => {
            const updated = scene.cubes.filter((_, i) => i !== idx);
            // make sure to update visibility state
            setVisibleCubes(prev => {
              const newVisibility = { ...prev };
              delete newVisibility[idx];
              return newVisibility;
            })
            updateScene("cubes", updated);
          }} style={{ color: 'red', marginTop: '0.5rem' }}>Delete Cube</button>
        </div>
      ))}


      <button onClick={() => {
        const newCube = {
          position: [0, -2.5, 0],
          scale: 2,
          rotation: [0, 45, 0],
          material: {
            diffuse: [0.5, 0.8, 0.2],
            specular: [1, 1, 1],
            intensity: [0.1, 0.3, 0.6],
            shininess: 32,
            alpha: 1.0,
            eta: 1.0,
            hasAlphaAt: false
          }
        };
        // make sure to update visibility state
        setVisibleCubes(prev => ({
          ...prev,
          [scene.cubes.length]: true
        }));
        updateScene("cubes", [...scene.cubes, newCube]);
      }} style={{ marginBottom: '2rem' }}>
        + Add Cube
      </button>




      <br></br>
      <button onClick={handleRender} style={{ marginTop: '1rem' }}>
        (Re)render Scene
      </button>
    </div>
    </div>
  );
}

export default App;
