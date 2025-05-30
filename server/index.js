const express = require('express');
const fs = require('fs');
const { exec } = require('child_process');
const cors = require('cors');
const path = require('path');

const app = express();
const PORT = 3001;

let scene = {
  width: 400,
  height: 300,
  spheres: []
};

app.use(cors()); // Allow requests from React dev server
app.use(express.json()); // Parse JSON bodies

// POST /add_sphere
app.post('/add_sphere', (req, res) => {
  scene.spheres.push(req.body);
  const scenePath = path.resolve(__dirname, '../scene.json');
    fs.writeFileSync(scenePath, JSON.stringify(scene));

    const exePath = path.resolve(__dirname, '../backend/build/raycaster.exe');
    const outputPath = path.resolve(__dirname, '../backend/output.ppm');
    exec(`"${exePath}" "${scenePath}" "${outputPath}"`, (err) => {


        if (err) {
        console.error(err);
        return res.status(500).send('Render failed');
        }
        res.send('Sphere added and rendered');
    });
});

// GET /render
app.get('/render', (req, res) => {
  res.setHeader('Content-Type', 'image/x-portable-pixmap');
  const outputPath = path.resolve(__dirname, '../backend/output.ppm');
fs.createReadStream(outputPath).pipe(res);

});

// Start server
app.listen(PORT, () => {
  console.log(`Server running at http://localhost:${PORT}`);
});
