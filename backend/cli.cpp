#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdint>
#include "json.hpp"

using json = nlohmann::json;

struct Sphere {
    float x, y, z, radius;
    uint8_t r, g, b;
};

// Declare functions from raycaster.cpp
void addSphere(const Sphere& s);
void clearScene();
void render(const std::string& outputFile, int width, int height);

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: raycaster <scene.json path> <output.ppm path>\n";
        return 1;
    }

    std::string sceneFile = argv[1];
    std::string outputFile = argv[2];

    std::ifstream file(sceneFile);
    if (!file.is_open()) {
        std::cerr << "Failed to open " << sceneFile << "\n";
        return 1;
    }

    json scene;
    try {
        file >> scene;
    } catch (...) {
        std::cerr << "Invalid JSON in " << sceneFile << "\n";
        return 1;
    }

    clearScene();
    for (const auto& obj : scene["spheres"]) {
        Sphere s = {
            obj["x"], obj["y"], obj["z"], obj["radius"],
            obj["r"], obj["g"], obj["b"]
        };
        addSphere(s);
    }

    render(outputFile, scene["width"], scene["height"]);
    return 0;
}
