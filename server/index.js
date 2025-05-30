const express = require('express');
const cors = require('cors');
const fs = require('fs');
const path = require('path');
const { exec } = require('child_process');

const app = express();
const PORT = 3001;

app.use(cors());
app.use(express.json());

const scenePath = path.resolve(__dirname, '../scene.json');
const outputPath = path.resolve(__dirname, '../backend/output.ppm');
const exePath = path.resolve(__dirname, '../backend/build/raycaster.exe');

// POST /scene
app.post('/scene', (req, res) => {
  try {
    fs.writeFileSync(scenePath, JSON.stringify(req.body, null, 2));
    res.send('Scene received');
  } catch (err) {
    console.error(err);
    res.status(500).send('Failed to write scene');
  }
});

// GET /render
app.get('/render', (req, res) => {
  exec(`"${exePath}" "${scenePath}" "${outputPath}"`, (err, stdout, stderr) => {
    if (err) {
      console.error(stderr);
      return res.status(500).send('Render failed');
    }

    // fs.readFile(outputPath, (err, data) => {
    //     if (err) return res.status(500).send('Cannot read output file');

    //     console.log('File starts with:', data.slice(0, 10).toString());
    //     res.setHeader('Content-Type', 'image/x-portable-pixmap');
    //     res.send(data);
    // });


    fs.access(outputPath, fs.constants.F_OK, (err) => {
      if (err) return res.status(404).send('Output not found');
      res.setHeader('Content-Type', 'image/x-portable-pixmap');
      fs.createReadStream(outputPath).pipe(res);
    });
  });
});

app.listen(PORT, () => {
  console.log(`Server running at http://localhost:${PORT}`);
});
