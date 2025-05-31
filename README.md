# Illumina Rendr 🌟

A browser-based raytracer built with a React frontend and a C++ backend compiled to run via Node.js.

## 🚀 Getting Started (Local Setup)

This guide assumes you're running on a Unix-like environment (macOS/Linux/WSL). Windows users should use WSL or a Unix-compatible terminal for best compatibility.

---

## 🧱 Project Structure

```
illumina-rendr/
├── backend/
│   ├── build/                # Contains the compiled C++ raytracer (raycaster executable)
│   ├── index.js              # Node.js backend server
│   ├── raycaster_pro.cpp     # Scene rendering logic
│   ├── cli.cpp               # C++ CLI entry point
│   └── json.hpp              # nlohmann JSON library
├── frontend/
│   ├── src/
│   │   └── App.tsx           # Main React component
│   └── ...
└── README.md                 # This file
```

---

## ✅ Prerequisites

- Node.js (v16 or later recommended)
- npm
- A C++ compiler (g++, clang, or MSVC)
- Git (optional)

---

## ⚙️ Backend Setup

1. **Open Terminal #1** and navigate to the `backend/` directory:

   ```bash
   cd backend
   ```

2. **Install dependencies**:

   ```bash
   npm install
   ```

3. **Compile your C++ raytracer** (make sure `raycaster.cpp` and dependencies are up-to-date):

   ```bash
   g++ cli.cpp -o build/raycaster
   ```

   > This will generate the executable at `backend/build/raycaster`

4. **Start the Node.js server**:

   ```bash
   node index.js
   ```

   You should see:
   ```
   Server running at http://localhost:3001
   ```

---

## 💻 Frontend Setup

1. **Open Terminal #2** and navigate to the `frontend/` directory:

   ```bash
   cd frontend
   ```

2. **Install dependencies**:

   ```bash
   npm install
   ```

3. **Start the React development server**:

   ```bash
   npm run dev
   ```

   It should automatically open `http://localhost:3000` (or similar) in your browser.

---

## 🧪 Usage

- Modify scene objects (spheres, lights, cubes) through the UI.
- Click **"(Re)render Scene"** to send your scene to the backend.
- The raytracer renders your scene and returns a PPM image, which is drawn directly onto the canvas.

---

## 🤝 Contributing

Pull requests welcome! Let’s make simple raytracing more fun and accessible for everyone.

---

## 📄 License

MIT © 2025
