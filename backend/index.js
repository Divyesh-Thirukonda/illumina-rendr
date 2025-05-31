const express = require('express');
const { spawn } = require('child_process');
const cors = require('cors');
const path = require('path');

const app = express();
const PORT = 3001;

app.use(cors());
app.use(express.json({ limit: '5mb' }));

app.post('/render', async (req, res) => {
  const sceneData = JSON.stringify(req.body);

  // Since index.js is in backend/, and so is build/
  const raycasterPath = path.join(__dirname, 'build', 'raycaster.exe');

  const raycaster = spawn(raycasterPath, [], {
    cwd: __dirname,
    stdio: ['pipe', 'pipe', 'inherit'],
  });

  let outputData = '';

  raycaster.stdout.on('data', (chunk) => {
    outputData += chunk.toString();
  });

  raycaster.on('error', (err) => {
    console.error('Spawn failed:', err);
    res.status(500).send('Failed to run raycaster');
  });

  raycaster.on('close', (code) => {
    if (code !== 0) {
      return res.status(500).send(`Raycaster exited with code ${code}`);
    }
    res.setHeader('Content-Type', 'text/plain');
    res.send(outputData);
  });

  raycaster.stdin.write(sceneData);
  raycaster.stdin.end();
});

app.listen(PORT, () => {
  console.log(`Server running at http://localhost:${PORT}`);
});
