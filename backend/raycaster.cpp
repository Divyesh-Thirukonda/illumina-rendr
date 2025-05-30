#include <vector>
#include <fstream>
#include <cmath>
#include <cstdint>

struct Sphere {
    float x, y, z, radius;
    uint8_t r, g, b;
};

static std::vector<Sphere> spheres;

void addSphere(const Sphere& s) {
    spheres.push_back(s);
}

void clearScene() {
    spheres.clear();
}


void render(const std::string& outputFile, int width, int height) {
    std::ofstream ofs(outputFile, std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";

    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            float cx = (i - width / 2.0f) / (float)width;
            float cy = (j - height / 2.0f) / (float)height;
            uint8_t r = 0, g = 0, b = 0;
            for (const auto& s : spheres) {
                float dx = cx - s.x, dy = cy - s.y;
                if (std::sqrt(dx * dx + dy * dy) < s.radius) {
                    r = s.r; g = s.g; b = s.b;
                }
            }
            ofs << r << g << b;
        }
    }
}
