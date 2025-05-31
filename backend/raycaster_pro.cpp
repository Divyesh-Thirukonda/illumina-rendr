#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <limits>
#include <random>
#include <stack>
#include <cstdint>

using namespace std;

const double pi = 3.14159265358979323846;

struct Vec3 {
    // Vector of 3 elems
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double x, double y, double z) : x(x), y(y), z(z) {}

    // standard vector operations
    Vec3 operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator*(double s) const { return Vec3(x * s, y * s, z * s); }
    Vec3 operator/(double s) const { return Vec3(x / s, y / s, z / s); }
    
    Vec3 clone() const {
        return Vec3(x, y, z);
    }    

    // dot product
    double dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }

    // cross product
    Vec3 cross(const Vec3& v) const { return Vec3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x); }

    Vec3 normalize() const {
        // normalize function
        double len = sqrt(x * x + y * y + z * z);
        return Vec3(x / len, y / len, z / len);
    }

    Vec3 normalizeTo1() const {
        // normalize function
        double sum = x + y + z;
        return Vec3(x / sum, y / sum, z / sum);
    }

    double norm() const {
        return sqrt(x * x + y * y + z * z);
    }

    // helper function for clamp
    double clampValue(double val, double a, double b) const {
        if (val < a) return a;
        if (val > b) return b;
        return val;
    }

    Vec3 clamp(double a, double b) const {
        return Vec3(clampValue(x, a, b), clampValue(y, a, b), clampValue(z, a, b));
    }
};


struct Vec2 {
    // Vector of 2 elems
    double x, y;
    Vec2() : x(0), y(0) {}
    Vec2(double x, double y) : x(x), y(y) {}
    Vec2 operator*(double s) const { return Vec2(x * s, y * s); }
    Vec2 operator+(const Vec2& v) const { return Vec2(x + v.x, y + v.y); }
    Vec2 operator-(const Vec2& v) const { return Vec2(x - v.x, y - v.y); }
};

struct Triangle {
    Vec3 v0, v1, v2;
    Vec3 n0, n1, n2;
    Vec2 vt1, vt2, vt3;

    Vec3 diffuseColor;
    Vec3 specularColor;
    Vec3 colorIntensity;
    int specIntensity;
    double alpha;
    double eta;

    bool hasNormals = false;

    vector<vector<Vec3>> texture;

    bool hasAlphaAt;
    float alphaR;
    float alphaG;
    float alphaB;
    

    // Triangle(Vec3 _v0, Vec3 _v1, Vec3 _v2, Vec2 _vt1, Vec2 _vt2, Vec2 _vt3, Vec3 _diffuseColor, Vec3 _specularColor, Vec3 _colorIntensity, int _specIntensity, double _alpha, double _eta, vector<vector<Vec3>> _texture) :
    //     v0(_v0), v1(_v1), v2(_v2),
    //     vt1(_vt1), vt2(_vt2), vt3(_vt3),
    //     diffuseColor(_diffuseColor), specularColor(_specularColor), colorIntensity(_colorIntensity), specIntensity(_specIntensity),
    //     hasNormals(false), alpha(_alpha), eta(_eta), texture(_texture) {}

    // Triangle(Vec3 _v0, Vec3 _v1, Vec3 _v2, Vec3 _n0, Vec3 _n1, Vec3 _n2, Vec2 _vt1, Vec2 _vt2, Vec2 _vt3, Vec3 _diffuseColor, Vec3 _specularColor, Vec3 _colorIntensity, int _specIntensity, double _alpha, double _eta, vector<vector<Vec3>> _texture) :
    //     v0(_v0), v1(_v1), v2(_v2), 
    //     n0(_n0), n1(_n1), n2(_n2),
    //     vt1(_vt1), vt2(_vt2), vt3(_vt3),
    //     diffuseColor(_diffuseColor), specularColor(_specularColor), colorIntensity(_colorIntensity), specIntensity(_specIntensity),
    //     hasNormals(true), alpha(_alpha), eta(_eta),texture(_texture) {}

    Triangle() : v0(Vec3()), v1(Vec3()), v2(Vec3()), 
        n0(Vec3()), n1(Vec3()), n2(Vec3()), 
        vt1(Vec2()), vt2(Vec2()), vt3(Vec2()),
        diffuseColor(Vec3()), specularColor(Vec3()), 
        colorIntensity(Vec3()), specIntensity(0), 
        hasNormals(false), alpha(1.0), eta(1.0), 
        hasAlphaAt(false), alphaR(0), alphaG(0), alphaB(0) {}
    
    Triangle(Vec3 _v0, Vec3 _v1, Vec3 _v2, Vec2 _vt1, Vec2 _vt2, Vec2 _vt3, 
        Vec3 _diffuseColor, Vec3 _specularColor, Vec3 _colorIntensity, 
        int _specIntensity, double _alpha, double _eta, 
        vector<vector<Vec3>> _texture, bool _hasAlphaAt, 
        double _alphaR, double _alphaG, double _alphaB) 
   : v0(_v0), v1(_v1), v2(_v2), 
     vt1(_vt1), vt2(_vt2), vt3(_vt3), 
     diffuseColor(_diffuseColor), specularColor(_specularColor), 
     colorIntensity(_colorIntensity), specIntensity(_specIntensity), 
     hasNormals(false), alpha(_alpha), eta(_eta), 
     texture(_texture), hasAlphaAt(_hasAlphaAt), 
     alphaR(_alphaR), alphaG(_alphaG), alphaB(_alphaB) {}

    Triangle(Vec3 _v0, Vec3 _v1, Vec3 _v2, Vec3 _n0, Vec3 _n1, Vec3 _n2, 
            Vec2 _vt1, Vec2 _vt2, Vec2 _vt3, Vec3 _diffuseColor, 
            Vec3 _specularColor, Vec3 _colorIntensity, int _specIntensity, 
            double _alpha, double _eta, vector<vector<Vec3>> _texture, 
            bool _hasAlphaAt, double _alphaR, double _alphaG, double _alphaB) 
    : v0(_v0), v1(_v1), v2(_v2), 
        n0(_n0), n1(_n1), n2(_n2), 
        vt1(_vt1), vt2(_vt2), vt3(_vt3), 
        diffuseColor(_diffuseColor), specularColor(_specularColor), 
        colorIntensity(_colorIntensity), specIntensity(_specIntensity), 
        hasNormals(true), alpha(_alpha), eta(_eta), 
        texture(_texture), hasAlphaAt(_hasAlphaAt), 
        alphaR(_alphaR), alphaG(_alphaG), alphaB(_alphaB) {}
};

struct Sphere {
    Vec3 center;
    double radius;

    Vec3 diffuseColor; // Odr Odg Odb
    Vec3 specularColor; // Osr Osg Osb
    Vec3 colorIntensity; // ka kd ks
    int specIntensity; // n
    double alpha; // opacity
    double eta; // index of refraction

    vector<vector<Vec3>> texture;

    bool hasAlphaAt;
    double alphaR;
    double alphaG;
    double alphaB;
};

struct Ray {
    Vec3 origin;
    Vec3 direction;
};

struct Light {
    Vec3 position; 
    double w; // 1 for point, 0 for directional
    double intensity;
    double fattc1, fattc2, fattc3;
};

struct Scene {
    Vec3 eye;
    Vec3 viewDir;
    Vec3 upDir;

    double vfov;
    int width;
    int height;

    Vec3 bgColor;
    double bgEta;

    vector<Sphere> spheres;
    vector<Light> lights;

    vector<Vec3> vertices;
    vector<Vec3> normals;
    vector<Triangle> triangles;

    Vec3 dcCol;
    double dcAlphaMin, dcAlphaMax, dcDistMin, dcDistMax;

    bool hasDc = false;
    bool hasAt = false;

    bool parallel = false;
    double frustumHeight = 1.0;
    
    vector<Vec2> bvh;
    
};

bool rayTriangleIntersect(const Ray& ray, const Triangle& tri, double& t, Vec3& normal, Vec2& texCoord) {
    Vec3 edge1 = tri.v1 - tri.v0;
    Vec3 edge2 = tri.v2 - tri.v0;
    Vec3 h = ray.direction.cross(edge2);
    double a = edge1.dot(h);

    if (fabs(a) < 1e-6) return false;  // Ray is parallel

    double f = 1.0 / a;
    Vec3 s = ray.origin - tri.v0;
    double u = f * s.dot(h);

    if (u < 0.0 || u > 1.0) return false;

    Vec3 q = s.cross(edge1);
    double v = f * ray.direction.dot(q);

    if (v < 0.0 || u + v > 1.0) return false;

    t = f * edge2.dot(q);

    if (t <= 1e-6) return false;

    // Compute barycentric coordinates
    double w = 1.0 - u - v;

    if (tri.hasNormals) {
        double alpha = 1 - u - v;
        double beta = u;
        double gamma = v;

        // Interpolate the normals using the barycentric coordinates
        normal = tri.n0*alpha + tri.n1*beta + tri.n2*gamma;
        normal = normal.normalize();
    } else {
        // flat shading = use face normal
        normal = (edge1.cross(edge2)).normalize(); // (tri.v1 - tri.v0).cross(tri.v2 - tri.v0).normalize()
    }

    texCoord = tri.vt1*(1.0 - u - v) + tri.vt2*u + tri.vt3*v;

    return true;
}

bool raySphereIntersect(const Ray& ray, const Sphere& sphere, double& t) {
    // this entire function is just the quadratic formula method we covered in class

    Vec3 oc = ray.origin - sphere.center;
    double a = ray.direction.dot(ray.direction);
    double b = 2.0 * oc.dot(ray.direction);
    double c = oc.dot(oc) - sphere.radius * sphere.radius;
    double discriminant = b * b - 4 * a * c;
    if (discriminant < 0) return false;
    t = (-b - sqrt(discriminant)) / (2.0 * a);
    return t > 0; // >= 0 ?
}



double randomFloat(double min, double max) {
    return min + (rand() / (RAND_MAX + 1.0)) * (max - min);
}

// ec func to calculate the shadow factor on jittered rays: NO LONGER USED, JUST RETURNS 1 NOW.
double computeShadowFactor(const Scene& scene, const Vec3& intersectionPoint, const Vec3& normal, const Vec3& lightPosition, bool isPointLight) {
    int numRays = 1; // num jittered rays
    int numHits = 0;  // num rays that hit something

    for (int i = 0; i < numRays; ++i) {
        Vec3 jitteredLightPosition = lightPosition;
        if (isPointLight) {
            jitteredLightPosition.x += randomFloat(-0.1, 0.1);
            jitteredLightPosition.y += randomFloat(-0.1, 0.1);
            jitteredLightPosition.z += randomFloat(-0.1, 0.1);
        }

        Vec3 lightDir = ((jitteredLightPosition - intersectionPoint)).normalize();
        Ray shadowRay = {intersectionPoint, lightDir};

        bool inShadow = false;

        // check intersections (except itself)
        for (const auto& sphere : scene.spheres) {
            // if (sphere.center.x == intersectionPoint.x && sphere.center.y == intersectionPoint.y && sphere.center.z == intersectionPoint.z)
            //     continue;

            double tShadow;
            if (raySphereIntersect(shadowRay, sphere, tShadow) && tShadow > 0) {
                inShadow = true;
                break;
            }
        }

        for (const auto& triangle : scene.triangles) {
            double tShadow;
            Vec3 normal;
            Vec2 texCoord;
            if (rayTriangleIntersect(shadowRay, triangle, tShadow, normal, texCoord) && tShadow > 0) {
                inShadow = true;
                break;
            }
        }

        if (inShadow) {
            numHits++;
        }
    }

    double shadowFactor = (numRays - numHits) / double(numRays);
    

    // fraction of rays that are not in shadow (shadow factor)
    return 1; // * lastObjHit.alpha;
}



Vec3 getTextureColorAt(const Vec2& texCoord, const vector<vector<Vec3>>& texture) {
    int x = texCoord.x * texture[0].size();
    if (x < 0) x = 0;

    int y = texCoord.y * texture.size();
    if (y < 0) y = 0;
    
    x = x % texture[0].size();
    y = y % texture.size();

    return texture[y][x];
}


Vec3 shadeRaySphere(const Scene& scene, const Ray& ray, const Sphere& sphere, double t, int countBounces); // predefine for rec call

Vec3 shadeRayTriangle(const Scene& scene, const Ray& ray, const Triangle& tri, double t, const Vec3& normal, const Vec2& texCoord, int countBounces) {
    Vec3 intersectionPoint = ray.origin + ray.direction * t;
    Vec3 viewDir = ray.direction.normalize() * -1;

    Vec3 newDiffuseColor = tri.diffuseColor;
    if (tri.texture.size() > 0) {
        Vec3 textureColor = getTextureColorAt(texCoord, tri.texture);
        newDiffuseColor = textureColor;
    }
    // Ambient light
    Vec3 ambientColor = newDiffuseColor * tri.colorIntensity.x;
    Vec3 finalColor = ambientColor;

    for (const auto& light : scene.lights) {
        // Direction to the light
        Vec3 lightDir;
        bool isPointLight = light.w == 1;
        if (isPointLight) {
            lightDir = (light.position - intersectionPoint).normalize();
        } else {
            lightDir = light.position.normalize() * -1;
        }

        // Distance to light and attenuation
        double distanceToLight = (light.position - intersectionPoint).dot(light.position - intersectionPoint);
        distanceToLight = (light.position - intersectionPoint).norm();
        double attenuationFactor = 1;
        if (scene.hasAt) {
            attenuationFactor = 1.0f / (light.fattc1 + light.fattc2 * distanceToLight + light.fattc3 * distanceToLight * distanceToLight);
            if (attenuationFactor >= std::numeric_limits<double>::max()) {
                attenuationFactor = 1;
            }
        }

        // Shadow factor computation (checks if the point is in shadow)
        // double shadowFactor = computeShadowFactor(scene, intersectionPoint, normal, light.position, isPointLight);
        
        // For each light, we will cast a shadow ray to check if the point is in shadow
        double shadowFactor = 1;
        for (const auto& sphere : scene.spheres) {
            double tShadow;
            if (raySphereIntersect({intersectionPoint, lightDir}, sphere, tShadow) && tShadow > 0) {
                shadowFactor *= (1-sphere.alpha);
                break;
            }
        }
        for (const auto& triangle : scene.triangles) {
            double tShadow;
            Vec3 _normal;
            Vec2 texCoord;
            if (rayTriangleIntersect({intersectionPoint + normal * (1e-4), lightDir}, triangle, tShadow, _normal, texCoord) && tShadow > 0) {
                shadowFactor *= (1-triangle.alpha);
                break;
            }
        }

        // Only compute lighting if not in shadow
        if (shadowFactor > 0.0) {
            // Diffuse color
            double diffuseFactor = max(normal.dot(lightDir), 0.0);
            Vec3 diffuseColor = newDiffuseColor * tri.colorIntensity.y * diffuseFactor * light.intensity * shadowFactor;

            // Specular color
            Vec3 halfwayDir = (lightDir + viewDir).normalize();
            double specFactor = pow(max(normal.dot(halfwayDir), 0.0), tri.specIntensity);
            Vec3 specularColor = tri.specularColor * tri.colorIntensity.z * specFactor * light.intensity * shadowFactor;

            // Adding diffuse and specular components
            Vec3 lightingColor = diffuseColor + specularColor;

            // Apply attenuation
            lightingColor = lightingColor * attenuationFactor;

            // Final color update
            finalColor = finalColor + lightingColor;
        }
    }

    // Apply depth cueing if enabled
    Vec3 Icurr = finalColor.clamp(0, 1);
    
    if (scene.hasDc) {
        double dcAlpha = 1;

        if (t <= scene.dcDistMin) {
            dcAlpha = scene.dcAlphaMax;
        } else if (t >= scene.dcDistMax) {
            dcAlpha = scene.dcAlphaMin;
        } else {
            dcAlpha = scene.dcAlphaMin + (scene.dcAlphaMax - scene.dcAlphaMin) * ((scene.dcDistMax - t) / (scene.dcDistMax - scene.dcDistMin));
        }

        // Apply depth cueing color adjustment
        Icurr = (Icurr * dcAlpha) + (scene.dcCol * (1 - dcAlpha));
    }




    
    if (countBounces > 10) {
        return Icurr;
    }

    if (tri.colorIntensity.z > 0) {
        
        Vec3 I = viewDir;
        Vec3 N = normal;
        bool inside = I.dot(N) < 0;

        double eta_i = tri.eta;
        double eta_t = scene.bgEta;

        if (inside) {
            eta_i = scene.bgEta;
            eta_t = tri.eta;
            N = N * -1;
        }
        double F0 = pow((eta_i - eta_t) / (eta_i + eta_t), 2);
        double F_r = F0 + (1 - F0) * pow(1 - N.dot(I), 5);

        Vec3 R = (N * 2 * N.dot(I) - I).normalize();
        Ray reflectionRay = {intersectionPoint + N * 1e-4, R};
        Vec3 reflectionColor = scene.bgColor;
        
        double t_clos = numeric_limits<double>::max();        
        for (const auto& _sphere : scene.spheres) {
            double t_s;
            if (raySphereIntersect(reflectionRay, _sphere, t_s) && t_s > 0 && t_s < t_clos) {
                t_clos = t_s;
                // BVH (EC)
                // if intersectionPoint in rect bounds
                if (intersectionPoint.x < scene.bvh[0].y && intersectionPoint.x > scene.bvh[0].x && intersectionPoint.y < scene.bvh[1].y && intersectionPoint.y > scene.bvh[1].x && intersectionPoint.z < scene.bvh[2].y && intersectionPoint.z > scene.bvh[2].x) {
                    // and t moves past rect bounds
                    Vec3 hitPoint = intersectionPoint + R*t_s;
                    if (hitPoint.x < scene.bvh[0].y && hitPoint.x > scene.bvh[0].x && hitPoint.y < scene.bvh[1].y && hitPoint.y > scene.bvh[1].x && hitPoint.z < scene.bvh[2].y && hitPoint.z > scene.bvh[2].x) {
                        reflectionColor = shadeRaySphere(scene, reflectionRay, _sphere, t_s, countBounces + 1);
                    }
                }
            }
        }
        for (const auto& _triangle : scene.triangles) {
            double t_tr;
            Vec3 normal;
            Vec2 texCoord;
            if (rayTriangleIntersect(reflectionRay, _triangle, t_tr, normal, texCoord) && t_tr > 0 && t_tr < t_clos) {
                t_clos = t_tr;
                // BVH (EC)
                // if intersectionPoint in rect bounds
                if (intersectionPoint.x < scene.bvh[0].y && intersectionPoint.x > scene.bvh[0].x && intersectionPoint.y < scene.bvh[1].y && intersectionPoint.y > scene.bvh[1].x && intersectionPoint.z < scene.bvh[2].y && intersectionPoint.z > scene.bvh[2].x) {
                    // and t moves past rect bounds
                    Vec3 hitPoint = intersectionPoint + R*t_tr;
                    if (hitPoint.x < scene.bvh[0].y && hitPoint.x > scene.bvh[0].x && hitPoint.y < scene.bvh[1].y && hitPoint.y > scene.bvh[1].x && hitPoint.z < scene.bvh[2].y && hitPoint.z > scene.bvh[2].x) {
                        reflectionColor = shadeRayTriangle(scene, reflectionRay, _triangle, t_tr, normal, texCoord, countBounces + 1);
                    }
                }
            }
        }
        
        Icurr = (Icurr + reflectionColor * F_r);
    }

    if (tri.alpha < 1) {
        std::stack<double> eta_stack;
        eta_stack.push(scene.bgEta);
        
        Vec3 I = viewDir;
        Vec3 N = normal;
        bool inside = I.dot(N) < 0;

        double eta_i = abs(1-tri.eta);
        double eta_t = eta_stack.top();

        if (inside) {
            eta_i = eta_stack.top();
            eta_t = abs(1-tri.eta);
            N = N * -1;
        }
        double F0 = pow((eta_i - eta_t) / (eta_i + eta_t), 2);
        double F_r = F0 + (1 - F0) * pow(1 - N.dot(I), 5);

        // Compute refraction
        Vec3 refractionColor = scene.bgColor;
        double etaRatio = eta_i / eta_t;
        double cosThetaI = N.dot(I);
        double sin2ThetaT = etaRatio * etaRatio * (1 - cosThetaI * cosThetaI);

        if (sin2ThetaT <= 1.0) { // No total internal reflection
            double cosThetaT = sqrt(1 - sin2ThetaT);
            // Vec3 T = (N * -1) * cosThetaT + ((N * cosThetaI) - I)*etaRatio;
            Vec3 T = ((N * -1) * sqrt(1 - (etaRatio * etaRatio) * (1 - (N.dot(I)) * (N.dot(I))))) + (((N*(N.dot(I))) - I)*etaRatio);
            // Vec3 T = (I * etaRatio + N * (etaRatio * cosThetaI - cosThetaT));
            
            T = T.normalize();

            Ray refractionRay = {intersectionPoint - N * 1e-4, T};

            if (inside) {
                eta_stack.push(eta_t);
            } else {
                eta_stack.pop();
            }

            double t_clos = numeric_limits<double>::max();
            for (const auto& _sphere : scene.spheres) {
                double t_s;
                if (raySphereIntersect(refractionRay, _sphere, t_s) && t_s > 0 && t_s < t_clos) {
                    t_clos = t_s;
                    refractionColor = shadeRaySphere(scene, refractionRay, _sphere, t_s, countBounces + 1);
                }
            }
            for (const auto& _triangle : scene.triangles) {
                double t_tr;
                Vec3 normal;
                Vec2 texCoord;
                if (rayTriangleIntersect(refractionRay, _triangle, t_tr, normal, texCoord) && t_tr > 0 && t_tr < t_clos) {
                    t_clos = t_tr;
                    refractionColor = shadeRayTriangle(scene, refractionRay, _triangle, t_tr, normal, texCoord, countBounces + 1);
                }
            }
        }

        

        // final color contribution
        
        // I did the extra credit part for the attenuation lighting
        if (tri.hasAlphaAt) {
            Vec3 alphaAtValRGB = Vec3( pow(2.7183, -tri.alphaR*t), pow(2.7183, -tri.alphaG*t), pow(2.7183, -tri.alphaB*t) );
            float Vr = (1 - F_r) * alphaAtValRGB.x * refractionColor.x;
            float Vg = (1 - F_r) * alphaAtValRGB.y * refractionColor.y;
            float Vb = (1 - F_r) * alphaAtValRGB.z * refractionColor.z;
            Icurr = (Icurr + Vec3(Vr, Vg, Vb));
        } else {
            Icurr = (Icurr + refractionColor * (1 - F_r) * (pow(2.7183, -tri.alpha*t)));

            // Or:
            // Icurr = (Icurr + refractionColor * (1 - F_r) * (1 - tri.alpha));
        }
    }


    return Icurr;
}


Vec3 getTextureColor(const Vec3& intersectionPoint, const Sphere& sphere) {
    Vec3 intersectionPointMinusCenter = intersectionPoint - sphere.center;
    
    // Convert 3D point to 2D texture coordinates (u, v)
    double u = 0.5 + (atan2(intersectionPointMinusCenter.y, intersectionPointMinusCenter.x) / (2 * pi));
    double v = 0.5 - (asin(intersectionPointMinusCenter.z / sphere.radius) / pi);

    // Map u, v to pixel coordinates (nearest neighbor)
    int i = round(u * (sphere.texture[0].size() - 1));
    int j = round(v * (sphere.texture.size() - 1));

    // Clamp values to ensure they stay within valid range
    i = std::max(0, std::min(i, (int)sphere.texture[0].size() - 1));
    j = std::max(0, std::min(j, (int)sphere.texture.size() - 1));

    Vec3 myCol = sphere.texture[j][i];

    return myCol;  // Return the nearest pixel color
}



Vec3 shadeRaySphere(const Scene& scene, const Ray& ray, const Sphere& sphere, double t, int countBounces) {
    Vec3 intersectionPoint = ray.origin + ray.direction * t;
    Vec3 normal = (intersectionPoint - sphere.center).normalize();

    Vec3 viewDir = ray.direction.normalize() * -1; // direction towards the eye

    Vec3 newDiffuseColor = sphere.diffuseColor;

    if (sphere.texture.size() > 0) {
        Vec3 textureColor = getTextureColor(intersectionPoint, sphere);
        newDiffuseColor = textureColor;
    }

    Vec3 ambientColor = newDiffuseColor * sphere.colorIntensity.x; // ambient
    Vec3 finalColor = ambientColor;

    for (Light light : scene.lights) {
        Vec3 lightDir;
        bool isPointLight = light.w == 1;
        if (isPointLight) {
            lightDir = (light.position - intersectionPoint).normalize();
        } else {
            lightDir = light.position.normalize() * -1;
        }

        double distanceToLight = (light.position - intersectionPoint).dot(light.position - intersectionPoint);
        distanceToLight = (light.position - intersectionPoint).norm();
        double attenuationFactor = 1;
        if (scene.hasAt) {
            attenuationFactor = 1.0f / (light.fattc1 + light.fattc2 * distanceToLight + light.fattc3 * distanceToLight * distanceToLight);
            if (attenuationFactor >= std::numeric_limits<double>::max()) {
                attenuationFactor = 1;
            }
        }
        
        // double shadowFactor = computeShadowFactor(scene, intersectionPoint, normal, light.position, isPointLight);
        double shadowFactor = 1;
        for (const auto& sphere : scene.spheres) {
            double tShadow;
            if (raySphereIntersect({intersectionPoint, lightDir}, sphere, tShadow) && tShadow > 0) {
                shadowFactor *= (1-sphere.alpha);
                // break;
            }
        }
        for (const auto& triangle : scene.triangles) {
            double tShadow;
            Vec3 normal;
            Vec2 texCoord;
            if (rayTriangleIntersect({intersectionPoint, lightDir}, triangle, tShadow, normal, texCoord) && tShadow > 0) {
                shadowFactor *= (1-triangle.alpha);
                // break;
            }
        }

        // only compute lighting if not in shadow
        if (shadowFactor > 0.0) {
            Vec3 diffColor = newDiffuseColor * sphere.colorIntensity.y * max(normal.dot(lightDir), 0.0) * light.intensity * shadowFactor;

            Vec3 halfwayDir = (lightDir + viewDir).normalize(); // viewDir = I
            double specFactor = pow(max(normal.dot(halfwayDir), 0.0), sphere.specIntensity);
            Vec3 specularColor = sphere.specularColor * sphere.colorIntensity.z * specFactor * light.intensity * shadowFactor;

            Vec3 tsColor = (diffColor + specularColor) * attenuationFactor;
            finalColor = finalColor + tsColor;
        }
    }

    Vec3 Icurr = finalColor.clamp(0, 1);

    if (scene.hasDc) {
        double dcAlpha = 1;

        if (t <= scene.dcDistMin) {
            dcAlpha = scene.dcAlphaMax;
        } else if (t >= scene.dcDistMax) {
            dcAlpha = scene.dcAlphaMin;
        } else {
            dcAlpha = scene.dcAlphaMin + (scene.dcAlphaMax - scene.dcAlphaMin) * ((scene.dcDistMax - t) / (scene.dcDistMax - scene.dcDistMin));
        }
        Icurr = (Icurr * dcAlpha) + (scene.dcCol * (1-dcAlpha));
    }
    



    if (countBounces > 10) {
        return Icurr;
    }

    if (sphere.colorIntensity.z > 0) {
        // Reflections
        Vec3 I = viewDir;
        Vec3 N = normal;
        bool inside = I.dot(N) < 0;

        double eta_i = sphere.eta;
        double eta_t = scene.bgEta;

        if (inside) { //  if we are starting inside, invert N and swap eta
            eta_i = scene.bgEta;
            eta_t = sphere.eta;
            N = N * -1;
        }

        // F_r = F0 + (1-F0) (1- I*N)^5
        // F0 = ((n-1)/(n+1)) ^2
        double F0 = pow((eta_i - eta_t) / (eta_i + eta_t), 2);
        double F_r = F0 + (1 - F0) * pow(1 - N.dot(I), 5);

        // R = A+S = 2(N*I)N - I
        // R_lambda = color at point in direction R
        Vec3 R = (N * 2 * N.dot(I) - I).normalize();
        Ray reflectionRay = {intersectionPoint + N * 1e-4, R};
        Vec3 reflectionColor = scene.bgColor;
        
        double t_clos = numeric_limits<double>::max();        
        for (const auto& _sphere : scene.spheres) {
            double t_s;
            if (raySphereIntersect(reflectionRay, _sphere, t_s) && t_s > 0 && t_s < t_clos) {
                t_clos = t_s;
                // BVH (EC)
                // if intersectionPoint in rect bounds
                if (intersectionPoint.x < scene.bvh[0].y && intersectionPoint.x > scene.bvh[0].x && intersectionPoint.y < scene.bvh[1].y && intersectionPoint.y > scene.bvh[1].x && intersectionPoint.z < scene.bvh[2].y && intersectionPoint.z > scene.bvh[2].x) {
                    // and t moves past rect bounds
                    Vec3 hitPoint = intersectionPoint + R*t_s;
                    if (hitPoint.x < scene.bvh[0].y && hitPoint.x > scene.bvh[0].x && hitPoint.y < scene.bvh[1].y && hitPoint.y > scene.bvh[1].x && hitPoint.z < scene.bvh[2].y && hitPoint.z > scene.bvh[2].x) {
                        reflectionColor = shadeRaySphere(scene, reflectionRay, _sphere, t_s, countBounces + 1);
                    }
                }
            }
        }
        for (const auto& _triangle : scene.triangles) {
            double t_tr;
            Vec3 normal;
            Vec2 texCoord;
            if (rayTriangleIntersect(reflectionRay, _triangle, t_tr, normal, texCoord) && t_tr > 0 && t_tr < t_clos) {
                t_clos = t_tr;
                // BVH (EC)
                // if intersectionPoint in rect bounds
                if (intersectionPoint.x < scene.bvh[0].y && intersectionPoint.x > scene.bvh[0].x && intersectionPoint.y < scene.bvh[1].y && intersectionPoint.y > scene.bvh[1].x && intersectionPoint.z < scene.bvh[2].y && intersectionPoint.z > scene.bvh[2].x) {
                    // and t moves past rect bounds
                    Vec3 hitPoint = intersectionPoint + R*t_tr;
                    if (hitPoint.x < scene.bvh[0].y && hitPoint.x > scene.bvh[0].x && hitPoint.y < scene.bvh[1].y && hitPoint.y > scene.bvh[1].x && hitPoint.z < scene.bvh[2].y && hitPoint.z > scene.bvh[2].x) {
                        reflectionColor = shadeRayTriangle(scene, reflectionRay, _triangle, t_tr, normal, texCoord, countBounces + 1);
                    }
                }
            }
        }
        
        Icurr = (Icurr + reflectionColor * F_r);
    }

    
    if (sphere.alpha < 1) {
        // Refractions
        std::stack<double> eta_stack; // stack to keep track of etas
        eta_stack.push(scene.bgEta);
        
        Vec3 I = viewDir;
        Vec3 N = normal;
        bool inside = I.dot(N) < 0;

        double eta_i = abs(1-sphere.eta);
        double eta_t = eta_stack.top();

        if (inside) {
            eta_i = eta_stack.top();
            eta_t = abs(1-sphere.eta);
            N = N * -1;
        }
        double F0 = pow((eta_i - eta_t) / (eta_i + eta_t), 2);
        double F_r = F0 + (1 - F0) * pow(1 - N.dot(I), 5);

        // refraction variable computation
        Vec3 refractionColor = scene.bgColor;
        double etaRatio = eta_i / eta_t;
        double cosThetaI = N.dot(I) * -1;
        double sin2ThetaT = etaRatio * etaRatio * (1 - cosThetaI * cosThetaI);

        if (sin2ThetaT <= 1.0) { // No total internal reflection
            double cosThetaT = sqrt(1 - sin2ThetaT);
            // Vec3 T = (N * -1) * cosThetaT + ((N * cosThetaI) - I)*etaRatio;
            Vec3 T = ((N * -1) * sqrt(1 - (etaRatio * etaRatio) * (1 - (N.dot(I)) * (N.dot(I))))) + (((N*(N.dot(I))) - I)*etaRatio);
            // Vec3 T = (I * etaRatio + N * (etaRatio * cosThetaI - cosThetaT));
            T = T.normalize();

            Ray refractionRay = {intersectionPoint - N * 1e-4, T};

            if (inside) { // add eta_t if we are entering a new medium
                eta_stack.push(eta_t);
            } else { // exited medium means pop last eta
                eta_stack.pop();
            }

            // Set refraction color by shooting ray and checking if it hits any objects
            double t_clos = numeric_limits<double>::max();
            for (const auto& _sphere : scene.spheres) {
                double t_s;
                if (raySphereIntersect(refractionRay, _sphere, t_s) && t_s > 0 && t_s < t_clos) {
                    t_clos = t_s;
                    refractionColor = shadeRaySphere(scene, refractionRay, _sphere, t_s, countBounces + 1);
                }
            }
            for (const auto& _triangle : scene.triangles) {
                double t_tr;
                Vec3 normal;
                Vec2 texCoord;
                if (rayTriangleIntersect(refractionRay, _triangle, t_tr, normal, texCoord) && t_tr > 0 && t_tr < t_clos) {
                    t_clos = t_tr;
                    refractionColor = shadeRayTriangle(scene, refractionRay, _triangle, t_tr, normal, texCoord, countBounces + 1);
                }
            }
        }
        
        // final color contribution
        // I did the extra credit part for the attenuation lighting

        if (sphere.hasAlphaAt) {
            Vec3 alphaAtValRGB = Vec3( pow(2.7183, -sphere.alphaR*t), pow(2.7183, -sphere.alphaG*t), pow(2.7183, -sphere.alphaB*t) );
            float Vr = (1 - F_r) * alphaAtValRGB.x * refractionColor.x;
            float Vg = (1 - F_r) * alphaAtValRGB.y * refractionColor.y;
            float Vb = (1 - F_r) * alphaAtValRGB.z * refractionColor.z;
            Icurr = (Icurr + Vec3(Vr, Vg, Vb));
        } else {
            Icurr = (Icurr + refractionColor * (1 - F_r) * (pow(2.7183, -sphere.alpha*t))); // or (1 -sphere.alpha)
        }
        
    }

    return Icurr;
}

void addSphere(float x, float y, float z, float radius, const Vec3& diffuseColor, const Vec3& specularColor, const Vec3& colorIntensity, int specIntensity, double eta, double alpha, const vector<vector<Vec3>>& texture, Scene& scene) {
    Sphere sphere;
    sphere.center = Vec3(x, y, z);
    sphere.radius = radius;
    sphere.diffuseColor = diffuseColor;
    sphere.specularColor = specularColor;
    sphere.colorIntensity = colorIntensity;
    sphere.specIntensity = specIntensity;
    sphere.eta = eta;
    sphere.alpha = alpha;
    sphere.texture = texture;

    scene.spheres.push_back(sphere);
}

void render(const Scene& scene, std::ostream& outFile) {
    outFile << "P3\n" << scene.width << " " << scene.height << "\n255\n";

    // this is all math either from lecture slides or self explanatory
    Vec3 right = scene.viewDir.cross(scene.upDir).normalize();
    Vec3 up = right.cross(scene.viewDir).normalize();

    double aspectRatio = (double)scene.width / scene.height;
    // ec functionality: if we are doing a parallel scene, then set height to frustum
    double viewHeight = scene.parallel ? scene.frustumHeight : 2 * tan((scene.vfov * pi / 180) / 2);
    double viewWidth = viewHeight * aspectRatio;

    Vec3 center = scene.eye + scene.viewDir.normalize(); // moving in view dir from eye

    Vec3 pixel00Loc = center - (right * (viewWidth / 2)) + (up * (viewHeight / 2)); // top-left pixel

    // change in pixel position
    Vec3 pixelDeltaU = right * (viewWidth / scene.width);
    Vec3 pixelDeltaV = up * (viewHeight / scene.height);

    for (int j = 0; j < scene.height; ++j) {
        for (int i = 0; i < scene.width; ++i) {
            // center of the current pixel in world space
            Vec3 pixelCenter = pixel00Loc + pixelDeltaU * (i + 0.5) - pixelDeltaV * (j + 0.5);

            // direction of the ray from the eye to the pixel center
            // if it is a parallel scene, just use viewDir
            Vec3 rayDir =  scene.parallel ? scene.viewDir.normalize() : (pixelCenter - scene.eye).normalize();

            // we either want to get rays going from the eyes or parallely outwards (based on pixelCenter)
            Vec3 rayOrigin = scene.parallel ? pixelCenter : scene.eye;
            Ray ray = {rayOrigin, rayDir};

            // initially, pixel is the bg color
            Vec3 pixelColor = scene.bgColor;

            // find the closest intersecting sphere
            double closestT = numeric_limits<double>::max();
            for (const auto& sphere : scene.spheres) {
                double t;
                if (raySphereIntersect(ray, sphere, t) && t < closestT) {
                    closestT = t;
                    pixelColor = shadeRaySphere(scene, ray, sphere, t, 0);
                }
            }
            for (const auto& triangle : scene.triangles) {
                double t;
                Vec3 normal;
                Vec2 texCoord;
                if (rayTriangleIntersect(ray, triangle, t, normal, texCoord) && t < closestT) {
                    closestT = t;
                    pixelColor = shadeRayTriangle(scene, ray, triangle, t, normal, texCoord, 0);
                }
            }
            

            // write computed pixel color converting from [0,1] to [0,255]
            outFile << (int)(pixelColor.x * 255) << " " << (int)(pixelColor.y * 255) << " " << (int)(pixelColor.z * 255) << " ";
        }
        outFile << "\n";
    }
}

vector<vector<Vec3>> loadTexture(const string& textureFilename) {
    ifstream file(textureFilename);
    string ppmType;
    int ppmWidth, ppmHeight, ppmMaxColor;

    // Read PPM header
    file >> ppmType >> ppmWidth >> ppmHeight >> ppmMaxColor;

    // Initialize the texture 2D array
    vector<vector<Vec3>> texture(ppmHeight, vector<Vec3>(ppmWidth));

    // Read pixel data
    for (int i = 0; i < ppmHeight; ++i) {
        for (int j = 0; j < ppmWidth; ++j) {
            int r, g, b;
            file >> r >> g >> b;
            // Normalize color values based on max color value
            texture[i][j] = Vec3((double)(r / (double)ppmMaxColor), (double)(g / (double)ppmMaxColor), (double)(b / (double)ppmMaxColor));
        }
    }

    return texture;
}
