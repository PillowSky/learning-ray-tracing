#include <iostream>
#include <random>
#include <glm/glm.hpp>
#include <Magick++.h>

using namespace std;
using namespace glm;

#define EPS 1e-4
default_random_engine rng;
uniform_real_distribution<double> uniform(0, 1);

enum struct Material {
	Diffuse,
	Specular,
	Refraction
};

struct Ray {
	dvec3 origin;
	dvec3 direction;
	
	Ray(dvec3 _origin, dvec3 _direction) : origin(_origin), direction(_direction) {};
};

struct Sphere {
	double radius;
	dvec3 position;
	dvec3 emission;
	dvec3 color;
	Material material;

	Sphere(double _radius, dvec3 _position, dvec3 _emission, dvec3 _color, Material _material) : radius(_radius), position(_position), emission(_emission), color(_color), material(_material) {};

	/*
	* Solve t^2*D.D + 2*t*(O-P).D + (O-P).(O-P)-R^2 = 0
	*
	* return distance to the incomming ray
	* return negative if nohit
	*/
	double intersect(const Ray& ray) const {
		dvec3 op = position - ray.origin;
		double leg = dot(op, ray.direction);
		double delta = leg * leg - dot(op, op) + radius * radius;
		if (delta < 0) {
			return delta;
		} else {
			delta = sqrt(delta);
			double distance;
			return (distance = leg - delta) > EPS ? distance : ((distance = leg + delta) > EPS ? distance : -1);
		}
	}
};

/*
* Scene: radius, position, emission, color, material
*/
Sphere spheres[] = {
	Sphere(1e5, dvec3( 1e5+1,40.8,81.6), dvec3(),dvec3(.75,.25,.25), Material::Diffuse),//Left
	Sphere(1e5, dvec3(-1e5+99,40.8,81.6),dvec3(),dvec3(.25,.25,.75), Material::Diffuse),//Rght
	Sphere(1e5, dvec3(50,40.8, 1e5),     dvec3(),dvec3(.75,.75,.75), Material::Diffuse),//Back
	Sphere(1e5, dvec3(50,40.8,-1e5+170), dvec3(),dvec3(),            Material::Diffuse),//Frnt
	Sphere(1e5, dvec3(50, 1e5, 81.6),    dvec3(),dvec3(.75,.75,.75), Material::Diffuse),//Botm
	Sphere(1e5, dvec3(50,-1e5+81.6,81.6),dvec3(),dvec3(.75,.75,.75), Material::Diffuse),//Top
	Sphere(16.5,dvec3(27,16.5,47),       dvec3(),dvec3(1,1,1)*.999,  Material::Specular),//Mirr
	Sphere(16.5,dvec3(73,16.5,78),       dvec3(),dvec3(1,1,1)*.999,  Material::Refraction),//Glas
	Sphere(600, dvec3(50,681.6-.27,81.6),dvec3(12,12,12),  dvec3(),  Material::Diffuse) //Lite
};
int spheres_n = sizeof(spheres) / sizeof(Sphere);

inline bool trace(const Ray& ray, double& distance, int& id) {
	double d = 0;
	distance = INFINITY;
	for (int i = spheres_n; i--;) {
		if ((d = spheres[i].intersect(ray)) >= 0 && d < distance) {
			distance = d;
			id = i;
		}
	}
	return distance < INFINITY;
}

dvec3 radiance(const Ray& ray, int depth) {
	double distance = 0;
	int id = 0;
	if (!trace(ray, distance, id)) return dvec3();
	
	const Sphere& obj = spheres[id];
	dvec3 hit = ray.origin + ray.direction * distance;
	dvec3 normal = normalize(hit - obj.position);
	bool into = dot(normal, ray.direction) < 0;
	dvec3 normal_forward = into ? normal : normal * -1.0;
	dvec3 color = obj.color;
	double max_intensity = glm::max(glm::max(color.r, color.g), color.b);
	if (++depth > 5) {
		if (uniform(rng) < max_intensity) {
			color *= 1 / max_intensity;
		} else {
			return obj.emission;
		}
	}

	if (obj.material == Material::Diffuse) {
		double angle = 2 * M_PI * uniform(rng);
		double distance_squared = uniform(rng);
		double distance = sqrt(distance_squared);
		dvec3 w = normal_forward;
		dvec3 u = normalize(cross(abs(w.x) > 0.1 ? dvec3(0, 1, 0) : dvec3(1, 0, 0), w));
		dvec3 v = cross(w, u);
		dvec3 direction = normalize(u * cos(angle) * distance + v * sin(angle) * distance + w * sqrt(1 - distance_squared));
		return obj.emission + color * radiance(Ray(hit, direction), depth);
	} else if (obj.material == Material::Specular) {
		return obj.emission + color * radiance(Ray(hit, reflect(ray.direction, normal)), depth);
	}

	Ray reflect_ray(hit, reflect(ray.direction, normal));
	double incident_normal = dot(ray.direction, normal_forward);
	double n_vaccum = 1;
	double n_glass = 1.5;
	double eta = into ? n_vaccum/n_glass : n_glass/n_vaccum;
	double cos2t;
	if ((cos2t = 1 - eta*eta * (1 - incident_normal*incident_normal)) < 0) {
		return obj.emission + color * radiance(reflect_ray, depth);
	} else {
		dvec3 refract_direction = normalize(ray.direction * eta - normal * ((into ? 1 : -1) * (incident_normal*eta + sqrt(cos2t))));
		double a = n_glass - n_vaccum;
		double b = n_glass + n_vaccum;
		double R0 = a * a / (b * b);
		double c = 1 - (into ? -incident_normal : dot(refract_direction, normal));
		double Re = R0 + (1 - R0) * c*c*c*c*c;
		double Tr = 1 - Re;
		double P = 0.25 + 0.5 * Re;
		double RP = Re / P;
		double TP = Tr / (1 - P);
		return obj.emission + color * (depth > 2 ? (uniform(rng) < P ?
			radiance(reflect_ray, depth) * RP : radiance(Ray(hit, refract_direction) , depth) * TP):
			radiance(reflect_ray, depth) * Re + radiance(Ray(hit, refract_direction) , depth) * Tr);
	}
}

int main(int argc, char* argv[]) {
	int width = 800;
	int height = 600;
	int sample = 1;

	Ray camera(dvec3(50,50,295.6), normalize(dvec3(0,-0.042612,-1)));
	dvec3 cx(width * .5135 / height, 0, 0);
	dvec3 cy = normalize(cross(cx, camera.direction)) * 0.5135;
	dvec3* buffer = new dvec3[width * height];
	dvec3 r;

	#pragma omp parallel for schedule(dynamic, 1) private(r)
	for (int y = 0; y < height; y++) {
		fprintf(stderr, "\rRendering (%d spp) %5.2f%%", sample * 4, 100.0 * y / (height - 1));
		for (int x = 0; x < width; x++) {
			for (int sy = 0, i = (height-y-1) * width + x; sy < 2; sy++) {
				for (int sx = 0; sx < 2; sx++, r = dvec3()) {
					for (int s = 0; s < sample; s++) {
						double r1 = 2 * uniform(rng);
						double r2 = 2 * uniform(rng);
						double dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
						double dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
						dvec3 d = cx*( ( (sx + 0.5 + dx)/2 + x)/width - 0.5) +
								cy*( ( (sy + 0.5 + dy)/2 + y)/height - 0.5) + camera.direction;
						r += radiance(Ray(camera.origin + d * 140.0, normalize(d)), 0) * (1.0 / sample);
					}
					buffer[i] += clamp(r, dvec3(0, 0, 0), dvec3(1, 1, 1)) * 0.25;
				}
			}
		}
	}
	fprintf(stderr, "\n");

	Magick::Image image(width, height, "RGB", Magick::DoublePixel, buffer);
	image.gamma(1.5);
	image.write("image.png");

	return EXIT_SUCCESS;
}
