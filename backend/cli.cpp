#include <iostream>
#include <fstream>
#include <string>
#include "json.hpp" // using nlohmann::json for JSON parsing
#include "raycaster_pro.cpp"

using json = nlohmann::json;
using namespace std;

Vec3 toVec3(const json& j) {
    return Vec3(j[0], j[1], j[2]);
}

Vec2 toVec2(const json& j) {
    return Vec2(j[0], j[1]);
}

int main(int argc, char** argv) {
    // Read all stdin into a string
    stringstream buffer;
    buffer << cin.rdbuf();

    json jscene;
    try {
        jscene = json::parse(buffer.str());
    } catch (const json::parse_error& e) {
        cerr << "Invalid JSON: " << e.what() << endl;
        return 1;
    }


    Scene scene;
    scene.width = jscene["width"];
    scene.height = jscene["height"];
    if (jscene.contains("eye")) scene.eye = toVec3(jscene["eye"]);
    if (jscene.contains("viewDir")) scene.viewDir = toVec3(jscene["viewDir"]);
    if (jscene.contains("upDir")) scene.upDir = toVec3(jscene["upDir"]);
    if (jscene.contains("vfov")) scene.vfov = jscene["vfov"];
    if (jscene.contains("bkgcolor")) {
        scene.bgColor = toVec3(jscene["bkgcolor"]["rgb"]);
        scene.bgEta = jscene["bkgcolor"]["eta"];
    }

    // Lights
    for (auto& l : jscene["lights"]) {
        Light light;
        light.position = toVec3(l["position"]);
        light.w = l["w"];
        light.intensity = l["intensity"];
        if (l.contains("attenuation")) {
            scene.hasAt = true;
            light.fattc1 = l["attenuation"][0];
            light.fattc2 = l["attenuation"][1];
            light.fattc3 = l["attenuation"][2];
        }
        scene.lights.push_back(light);
    }

    // Spheres
    for (auto& s : jscene["spheres"]) {
        Sphere sphere;
        sphere.center = toVec3(s["center"]);
        sphere.radius = s["radius"];
        sphere.diffuseColor = toVec3(s["diffuse"]);
        sphere.specularColor = toVec3(s["specular"]);
        sphere.colorIntensity = toVec3(s["intensity"]);
        sphere.specIntensity = s["shininess"];
        sphere.alpha = s["alpha"];
        sphere.eta = s["eta"];
        sphere.hasAlphaAt = s.value("hasAlphaAt", false);
        if (sphere.hasAlphaAt) {
            sphere.alphaR = s["alphaR"];
            sphere.alphaG = s["alphaG"];
            sphere.alphaB = s["alphaB"];
        }
        scene.spheres.push_back(sphere);
    }

    // Triangles
    for (auto& t : jscene["triangles"]) {
        Triangle tri;
        tri.v0 = toVec3(t["v0"]);
        tri.v1 = toVec3(t["v1"]);
        tri.v2 = toVec3(t["v2"]);
        if (t.contains("normals")) {
            tri.hasNormals = true;
            tri.n0 = toVec3(t["normals"][0]);
            tri.n1 = toVec3(t["normals"][1]);
            tri.n2 = toVec3(t["normals"][2]);
        }
        tri.vt1 = toVec2(t["vt0"]);
        tri.vt2 = toVec2(t["vt1"]);
        tri.vt3 = toVec2(t["vt2"]);
        tri.diffuseColor = toVec3(t["diffuse"]);
        tri.specularColor = toVec3(t["specular"]);
        tri.colorIntensity = toVec3(t["intensity"]);
        tri.specIntensity = t["shininess"];
        tri.alpha = t["alpha"];
        tri.eta = t["eta"];
        tri.hasAlphaAt = t.value("hasAlphaAt", false);
        if (tri.hasAlphaAt) {
            tri.alphaR = t["alphaR"];
            tri.alphaG = t["alphaG"];
            tri.alphaB = t["alphaB"];
        }
        scene.triangles.push_back(tri);
    }

    // Depth cueing
    if (jscene.contains("depthCueing")) {
        scene.hasDc = true;
        const auto& dc = jscene["depthCueing"];
        scene.dcCol = toVec3(dc["color"]);
        scene.dcAlphaMin = dc["alphaMin"];
        scene.dcAlphaMax = dc["alphaMax"];
        scene.dcDistMin = dc["distMin"];
        scene.dcDistMax = dc["distMax"];
    }

    // BVH bounds
    if (jscene.contains("bvh")) {
        for (int i = 0; i < 3; i++) {
            scene.bvh.push_back(Vec2(jscene["bvh"][i][0], jscene["bvh"][i][1]));
        }
    } else {
        scene.bvh = { Vec2(-10000, 10000), Vec2(-10000, 10000), Vec2(-10000, 10000) };
    }

    if (jscene.contains("parallel")) {
        scene.parallel = jscene["parallel"];
        if (scene.parallel && jscene.contains("frustumHeight")) {
            scene.frustumHeight = jscene["frustumHeight"];
        }
    }

    // Optional: vertices, normals, etc.
    for (auto& v : jscene.value("vertices", json::array())) {
        scene.vertices.push_back(toVec3(v));
    }
    for (auto& n : jscene.value("normals", json::array())) {
        scene.normals.push_back(toVec3(n));
    }

    render(scene, std::cout);
    return 0;
}
