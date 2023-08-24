/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Kevin Hopkins
 * *************************
*/

#ifdef WIN32
#include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
#include <GL/gl.h>
#include <GL/glut.h>
#elif defined(__APPLE__)
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <string.h>
#define _USE_MATH_DEFINES
#include <cmath>
#ifdef WIN32
#define strcasecmp _stricmp
#endif

#include <imageIO.h>

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char* filename = NULL;

// Different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

// You may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

// The field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
    double position[3];
    double color_diffuse[3];
    double color_specular[3];
    double normal[3];
    double shininess;
};

struct Triangle
{
    Vertex v[3];
};

struct Sphere
{
    double position[3];
    double color_diffuse[3];
    double color_specular[3];
    double shininess;
    double radius;
};

struct Light
{
    double position[3];
    double color[3];
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

int antialiased = 0;

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b);

// Helper function to normalize a vector
void normalize(double vec[3])
{
    double length = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
    vec[0] /= length;
    vec[1] /= length;
    vec[2] /= length;
}

// Helper function to compute the dot product of two vectors
double dot_product(double vec1[3], double vec2[3])
{
    return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

// Helper functions
// Calculate the length of a vector
double length(double vec[3])
{
    return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

// Subtract two vectors
void subtract(double vec1[3], double vec2[3], double result[3])
{
    result[0] = vec1[0] - vec2[0];
    result[1] = vec1[1] - vec2[1];
    result[2] = vec1[2] - vec2[2];
}

// Multiply a vector by a scalar
void multiply(double vec[3], double scalar, double result[3])
{
    result[0] = vec[0] * scalar;
    result[1] = vec[1] * scalar;
    result[2] = vec[2] * scalar;
}

// Add two vectors
void add(double vec1[3], double vec2[3], double result[3])
{
    result[0] = vec1[0] + vec2[0];
    result[1] = vec1[1] + vec2[1];
    result[2] = vec1[2] + vec2[2];
}

// Calculate the cross product of two vectors
void cross_product(double vec1[3], double vec2[3], double result[3])
{
    result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}

// Determine the intersection between a ray and a sphere
bool ray_sphere_intersection(double ray_origin[3], double ray_dir[3], Sphere sphere, double& t)
{
    double oc[3];
    subtract(ray_origin, sphere.position, oc);

    double a = dot_product(ray_dir, ray_dir);
    double b = 2.0 * dot_product(oc, ray_dir);
    double c = dot_product(oc, oc) - sphere.radius * sphere.radius;
    double discriminant = b * b - 4 * a * c;

    if (discriminant < 0)
    {
        return false;
    }
    else
    {
        t = (-b - sqrt(discriminant)) / (2.0 * a);
        return true;
    }
}

// Determine the intersection between a ray and a triangle
bool ray_triangle_intersection(double ray_origin[3], double ray_dir[3], Triangle triangle, double& t)
{
    double EPSILON = 1e-8;
    double e1[3], e2[3], h[3], s[3], q[3];
    subtract(triangle.v[1].position, triangle.v[0].position, e1);
    subtract(triangle.v[2].position, triangle.v[0].position, e2);
    cross_product(ray_dir, e2, h);
    double a = dot_product(e1, h);

    if (a > -EPSILON && a < EPSILON)
        return false;

    double f = 1 / a;
    subtract(ray_origin, triangle.v[0].position, s);
    double u = f * dot_product(s, h);
    if (u < 0.0 || u > 1.0)
        return false;

    cross_product(s, e1, q);
    double v = f * dot_product(ray_dir, q);
    if (v < 0.0 || u + v > 1.0)
        return false;

    t = f * dot_product(e2, q);
    if (t > EPSILON)
        return true;

    return false;
}

// Clamp function
double clamp(double value, double min_value, double max_value) {
    if (value < min_value) {
        return min_value;
    }
    else if (value > max_value) {
        return max_value;
    }
    else {
        return value;
    }
}

bool is_in_shadow(Vertex& vertex, Light& light)
{
    // Compute the direction from vertex position to the light source
    double light_dir[3] = { light.position[0] - vertex.position[0],
                            light.position[1] - vertex.position[1],
                            light.position[2] - vertex.position[2] };
    // Create variable light_dist that is a copy of light_dir
    double light_dist[3] = { light.position[0] - vertex.position[0],
                             light.position[1] - vertex.position[1],
                             light.position[2] - vertex.position[2] };
    double light_distance = length(light_dist);
    normalize(light_dist);

    // Launch a shadow ray from the vertex position towards the light source
    double min_t = INFINITY;
    int intersected_triangle = -1;
    int intersected_sphere = -1;

    // Check for triangle intersections
    for (int i = 0; i < num_triangles; i++)
    {
        double t;
        if (ray_triangle_intersection(vertex.position, light_dist, triangles[i], t))
        {
            if (t < min_t && t > 1e-8)
            {
                min_t = t;
                intersected_triangle = i;
            }
        }
    }

    // Check for sphere intersections
    for (int i = 0; i < num_spheres; i++)
    {
        double t;
        if (ray_sphere_intersection(vertex.position, light_dist, spheres[i], t))
        {
            if (t < min_t && t > 1e-8)
            {
                min_t = t;
                intersected_sphere = i;
            }
        }
    }

    // If the shadow ray hits an object before the light source, the point is in shadow
    if ((intersected_triangle != -1 || intersected_sphere != -1) && min_t < light_distance)
    {
        return true;
    }

    return false;
}


// Compute Phong shading for a point
void phong_shading(const Vertex& vertex, const Light& light,
    const double camera_pos[3], double light_color[3])
{
    // Calculate L (light direction), N (normal), R (reflected light direction), and V (view direction)
    double L[3] = { light.position[0] - vertex.position[0],
                    light.position[1] - vertex.position[1],
                    light.position[2] - vertex.position[2] };
    normalize(L);

    double N[3] = { vertex.normal[0],
                    vertex.normal[1],
                    vertex.normal[2] };
    normalize(N);

    double R[3] = { 2.0 * dot_product(N, L) * N[0] - L[0],
                    2.0 * dot_product(N, L) * N[1] - L[1],
                    2.0 * dot_product(N, L) * N[2] - L[2] };
    normalize(R);

    double V[3] = { camera_pos[0] - vertex.position[0],
                    camera_pos[1] - vertex.position[1],
                    camera_pos[2] - vertex.position[2] };
    normalize(V);

    // Calculate Phong shading components: diffuse and specular
    double diffuse = clamp(dot_product(L, N), 0.0, 1.0);
    double specular = clamp(pow(clamp(dot_product(R, V), 0.0, 1.0), vertex.shininess), 0.0, 1.0);

    // Calculate Phong shading for each color channel
    for (int i = 0; i < 3; ++i)
    {
        light_color[i] = light.color[i] * (vertex.color_diffuse[i] * diffuse +
            vertex.color_specular[i] * specular);
    }
}

void barycentric_coordinates(Triangle& triangle, double point[3], double& u, double& v)
{
    double v0[3], v1[3], v2[3];
    subtract(triangle.v[1].position, triangle.v[0].position, v0);
    subtract(triangle.v[2].position, triangle.v[0].position, v1);
    double temp[3];
    subtract(point, triangle.v[0].position, temp);
    double d00 = dot_product(v0, v0);
    double d01 = dot_product(v0, v1);
    double d11 = dot_product(v1, v1);
    double d20 = dot_product(temp, v0);
    double d21 = dot_product(temp, v1);
    double denom = d00 * d11 - d01 * d01;
    u = (d11 * d20 - d01 * d21) / denom;
    v = (d00 * d21 - d01 * d20) / denom;
}

void interpolate_color(const Triangle& triangle, double u, double v, double out_color[3], bool is_specular)
{
    const double* color_a = is_specular ? triangle.v[0].color_specular : triangle.v[0].color_diffuse;
    const double* color_b = is_specular ? triangle.v[1].color_specular : triangle.v[1].color_diffuse;
    const double* color_c = is_specular ? triangle.v[2].color_specular : triangle.v[2].color_diffuse;

    for (int i = 0; i < 3; i++)
    {
        out_color[i] = color_a[i] * (1 - u - v) + color_b[i] * u + color_c[i] * v;
    }
}

double interpolate_shininess(const Triangle& triangle, double u, double v)
{
    return triangle.v[0].shininess * (1 - u - v) + triangle.v[1].shininess * u + triangle.v[2].shininess * v;
}

void interpolate_normal(const Triangle& triangle, double u, double v, double out_normal[3])
{
    for (int i = 0; i < 3; i++)
    {
        out_normal[i] = triangle.v[0].normal[i] * (1 - u - v) + triangle.v[1].normal[i] * u + triangle.v[2].normal[i] * v;
    }
    normalize(out_normal);
}

void draw_scene()
{
    double aspect_ratio = (double)WIDTH / (double)HEIGHT; // Aspect ratio: this is never actually used!
    double camera_pos[3] = { 0, 0, 0 }; // Camera position
    double plane_dist = (0.5 * HEIGHT) / tan(fov * 0.5 * (M_PI / 180.0)); // Distance from camera to image plane

    for (unsigned int x = 0; x < WIDTH; x++)
    {
        for (unsigned int y = 0; y < HEIGHT; y++)
        {
            // Calculate the ray direction
            double ray_dir[3] = { ((double)x - 0.5 * (double)WIDTH),
                     ((double)y - 0.5 * (double)HEIGHT),
                     -plane_dist };
            normalize(ray_dir);

            // Implement raytracing steps:
            // 1. Find intersections with triangles and spheres
            double min_t = INFINITY;
            int intersected_triangle = -1;
            int intersected_sphere = -1;

            // Check for triangle intersections
            for (int i = 0; i < num_triangles; i++)
            {
                double t;
                if (ray_triangle_intersection(camera_pos, ray_dir, triangles[i], t))
                {
                    if (t < min_t)
                    {
                        min_t = t;
                        intersected_triangle = i;
                    }
                }
            }

            // Check for sphere intersections
            for (int i = 0; i < num_spheres; i++)
            {
                double t;
                if (ray_sphere_intersection(camera_pos, ray_dir, spheres[i], t))
                {
                    if (t < min_t)
                    {
                        min_t = t;
                        intersected_sphere = i;
                    }
                }
            }

            // Check if there is an intersection
            if (intersected_triangle != -1 || intersected_sphere != -1)
            {
                double final_color[3] = { 0, 0, 0 };

                Vertex intersection;
                intersection.position[0] = camera_pos[0] + min_t * ray_dir[0];
                intersection.position[1] = camera_pos[1] + min_t * ray_dir[1];
                intersection.position[2] = camera_pos[2] + min_t * ray_dir[2];
                intersection.color_diffuse[0] = 0;
                intersection.color_diffuse[1] = 0;
                intersection.color_diffuse[2] = 0;
                intersection.color_specular[0] = 0;
                intersection.color_specular[1] = 0;
                intersection.color_specular[2] = 0;
                intersection.normal[0] = 0;
                intersection.normal[1] = 0;
                intersection.normal[2] = 0;
                intersection.shininess = 0;

                if (intersected_sphere != -1) {
                    Sphere sphere = spheres[intersected_sphere];
                    // Set color properties
                    memcpy(intersection.color_diffuse, sphere.color_diffuse, 3 * sizeof(double));
                    memcpy(intersection.color_specular, sphere.color_specular, 3 * sizeof(double));
                    intersection.shininess = sphere.shininess;

                    // Calculate intersection normal for the sphere
                    subtract(intersection.position, sphere.position, intersection.normal);
                    normalize(intersection.normal);
                }
                else if (intersected_triangle != -1) {

                    Triangle triangle = triangles[intersected_triangle];
                    // Find barycentric coordinates of intersection point
                    double u, v;
                    barycentric_coordinates(triangle, intersection.position, u, v);

                    // Set color properties
                    interpolate_color(triangle, u, v, intersection.color_diffuse, false);
                    interpolate_color(triangle, u, v, intersection.color_specular, true);
                    intersection.shininess = interpolate_shininess(triangle, u, v);

                    // Calculate intersection normal for the triangle
                    interpolate_normal(triangle, u, v, intersection.normal);
                }


                // Add ambient light
                final_color[0] += ambient_light[0] * intersection.color_diffuse[0];
                final_color[1] += ambient_light[1] * intersection.color_diffuse[1];
                final_color[2] += ambient_light[2] * intersection.color_diffuse[2];

                // 2. Determine if intersection point is in shadow and calculate Phong shading
                for (int light_idx = 0; light_idx < num_lights; ++light_idx)
                {
                    Light light = lights[light_idx];
                    if (!is_in_shadow(intersection, light))
                    {
                        double light_color[3] = { 0, 0, 0 };

                        // Calculate Phong shading for the point if not in shadow
                        phong_shading(intersection, light, camera_pos, light_color);

                        // Add light contribution to the final color
                        final_color[0] += light_color[0];
                        final_color[1] += light_color[1];
                        final_color[2] += light_color[2];
                    }
                }

                // 3. Clamp color to [0, 1]
                unsigned char r = clamp(final_color[0], 0.0, 1.0) * 255;
                unsigned char g = clamp(final_color[1], 0.0, 1.0) * 255;
                unsigned char b = clamp(final_color[2], 0.0, 1.0) * 255;

                plot_pixel(x, y, r, g, b);
            }
            else
            {
                // Set background color if no intersection
                unsigned char r = 255;
                unsigned char g = 255;
                unsigned char b = 255;

                plot_pixel(x, y, r, g, b);
            }
        }
    }
    glEnd();
    glFlush();
}

#include <random>

// Generate a random double between 0 and 1
double random_double()
{
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    static std::mt19937 generator;
    return distribution(generator);
}

void draw_antialiased_scene()
{
    double aspect_ratio = (double)WIDTH / (double)HEIGHT; // Aspect ratio: this is never actually used!
    double camera_pos[3] = { 0, 0, 0 }; // Camera position
    double plane_dist = (0.5 * HEIGHT) / tan(fov * 0.5 * (M_PI / 180.0)); // Distance from camera to image plane
    const int SUPERSAMPLING_FACTOR = 2; // Change this value to increase/decrease supersampling grid size

    for (unsigned int x = 0; x < WIDTH; x++)
    {
        for (unsigned int y = 0; y < HEIGHT; y++)
        {
            double final_color_accum[3] = { 0, 0, 0 };

            for (int sx = 0; sx < SUPERSAMPLING_FACTOR; sx++) {

                for (int sy = 0; sy < SUPERSAMPLING_FACTOR; sy++) {

                    double x_supersampled = x + (double)sx / SUPERSAMPLING_FACTOR;
                    double y_supersampled = y + (double)sy / SUPERSAMPLING_FACTOR;

                    // Calculate the ray direction
                    double ray_dir[3] = { ((double)x_supersampled - 0.5 * (double)WIDTH),
                                          ((double)y_supersampled - 0.5 * (double)HEIGHT),
                                          -plane_dist };
                    normalize(ray_dir);

                    // Implement raytracing steps:
                    // 1. Find intersections with triangles and spheres
                    double min_t = INFINITY;
                    int intersected_triangle = -1;
                    int intersected_sphere = -1;

                    // Check for triangle intersections
                    for (int i = 0; i < num_triangles; i++)
                    {
                        double t;
                        if (ray_triangle_intersection(camera_pos, ray_dir, triangles[i], t))
                        {
                            if (t < min_t)
                            {
                                min_t = t;
                                intersected_triangle = i;
                            }
                        }
                    }

                    // Check for sphere intersections
                    for (int i = 0; i < num_spheres; i++)
                    {
                        double t;
                        if (ray_sphere_intersection(camera_pos, ray_dir, spheres[i], t))
                        {
                            if (t < min_t)
                            {
                                min_t = t;
                                intersected_sphere = i;
                            }
                        }
                    }

                    // Check if there is an intersection
                    if (intersected_triangle != -1 || intersected_sphere != -1)
                    {
                        double final_color[3] = { 0, 0, 0 };

                        Vertex intersection;
                        intersection.position[0] = camera_pos[0] + min_t * ray_dir[0];
                        intersection.position[1] = camera_pos[1] + min_t * ray_dir[1];
                        intersection.position[2] = camera_pos[2] + min_t * ray_dir[2];
                        intersection.color_diffuse[0] = 0;
                        intersection.color_diffuse[1] = 0;
                        intersection.color_diffuse[2] = 0;
                        intersection.color_specular[0] = 0;
                        intersection.color_specular[1] = 0;
                        intersection.color_specular[2] = 0;
                        intersection.normal[0] = 0;
                        intersection.normal[1] = 0;
                        intersection.normal[2] = 0;
                        intersection.shininess = 0;

                        if (intersected_sphere != -1) {
                            Sphere sphere = spheres[intersected_sphere];
                            // Set color properties
                            memcpy(intersection.color_diffuse, sphere.color_diffuse, 3 * sizeof(double));
                            memcpy(intersection.color_specular, sphere.color_specular, 3 * sizeof(double));
                            intersection.shininess = sphere.shininess;

                            // Calculate intersection normal for the sphere
                            subtract(intersection.position, sphere.position, intersection.normal);
                            normalize(intersection.normal);
                        }
                        else if (intersected_triangle != -1) {

                            Triangle triangle = triangles[intersected_triangle];
                            // Find barycentric coordinates of intersection point
                            double u, v;
                            barycentric_coordinates(triangle, intersection.position, u, v);

                            // Set color properties
                            interpolate_color(triangle, u, v, intersection.color_diffuse, false);
                            interpolate_color(triangle, u, v, intersection.color_specular, true);
                            intersection.shininess = interpolate_shininess(triangle, u, v);

                            // Calculate intersection normal for the triangle
                            interpolate_normal(triangle, u, v, intersection.normal);
                        }


                        // Add ambient light
                        final_color[0] += ambient_light[0] * intersection.color_diffuse[0];
                        final_color[1] += ambient_light[1] * intersection.color_diffuse[1];
                        final_color[2] += ambient_light[2] * intersection.color_diffuse[2];

                        // 2. Determine if intersection point is in shadow and calculate Phong shading
                        for (int light_idx = 0; light_idx < num_lights; ++light_idx)
                        {
                            Light light = lights[light_idx];
                            if (!is_in_shadow(intersection, light))
                            {
                                double light_color[3] = { 0, 0, 0 };

                                // Calculate Phong shading for the point if not in shadow
                                phong_shading(intersection, light, camera_pos, light_color);

                                // Add light contribution to the final color
                                final_color[0] += light_color[0];
                                final_color[1] += light_color[1];
                                final_color[2] += light_color[2];
                            }
                        }

                        // 3. Clamp color to [0, 1]
                        unsigned char r = clamp(final_color[0], 0.0, 1.0) * 255;
                        unsigned char g = clamp(final_color[1], 0.0, 1.0) * 255;
                        unsigned char b = clamp(final_color[2], 0.0, 1.0) * 255;

                        //plot_pixel(x, y, r, g, b);

                        // Add the calculated color to the accumulated color
                        final_color_accum[0] += final_color[0];
                        final_color_accum[1] += final_color[1];
                        final_color_accum[2] += final_color[2];
                    }
                    else
                    {
                        // Set background color if no intersection
                        final_color_accum[0] += 1;
                        final_color_accum[1] += 1;
                        final_color_accum[2] += 1;
                    }
                }
            }
            // Calculate the average color
            final_color_accum[0] /= (SUPERSAMPLING_FACTOR * SUPERSAMPLING_FACTOR);
            final_color_accum[1] /= (SUPERSAMPLING_FACTOR * SUPERSAMPLING_FACTOR);
            final_color_accum[2] /= (SUPERSAMPLING_FACTOR * SUPERSAMPLING_FACTOR);

            // Plot the pixel with the averaged color
            plot_pixel(x, y, final_color_accum[0] * 255, final_color_accum[1] * 255, final_color_accum[2] * 255);
        }
    }
    glEnd();
    glFlush();
}

void draw_monte_carlo_scene()
{
    double aspect_ratio = (double)WIDTH / (double)HEIGHT;
    double camera_pos[3] = { 0, 0, 0 };
    double plane_dist = (0.5 * HEIGHT) / tan(fov * 0.5 * (M_PI / 180.0));
    const int NUM_SAMPLES = 4; // Change this value to increase/decrease the number of samples

    for (unsigned int x = 0; x < WIDTH; x++)
    {
        for (unsigned int y = 0; y < HEIGHT; y++)
        {
            double final_color_accum[3] = { 0, 0, 0 };

            for (int sample = 0; sample < NUM_SAMPLES; sample++)
            {
                double x_random = x + random_double();
                double y_random = y + random_double();

                // Calculate the ray direction
                double ray_dir[3] = { ((double)x_random - 0.5 * (double)WIDTH),
                                      ((double)y_random - 0.5 * (double)HEIGHT),
                                      -plane_dist };
                normalize(ray_dir);

                // Implement raytracing steps:
                    // 1. Find intersections with triangles and spheres
                double min_t = INFINITY;
                int intersected_triangle = -1;
                int intersected_sphere = -1;

                // Check for triangle intersections
                for (int i = 0; i < num_triangles; i++)
                {
                    double t;
                    if (ray_triangle_intersection(camera_pos, ray_dir, triangles[i], t))
                    {
                        if (t < min_t)
                        {
                            min_t = t;
                            intersected_triangle = i;
                        }
                    }
                }

                // Check for sphere intersections
                for (int i = 0; i < num_spheres; i++)
                {
                    double t;
                    if (ray_sphere_intersection(camera_pos, ray_dir, spheres[i], t))
                    {
                        if (t < min_t)
                        {
                            min_t = t;
                            intersected_sphere = i;
                        }
                    }
                }

                // Check if there is an intersection
                if (intersected_triangle != -1 || intersected_sphere != -1)
                {
                    double final_color[3] = { 0, 0, 0 };

                    Vertex intersection;
                    intersection.position[0] = camera_pos[0] + min_t * ray_dir[0];
                    intersection.position[1] = camera_pos[1] + min_t * ray_dir[1];
                    intersection.position[2] = camera_pos[2] + min_t * ray_dir[2];
                    intersection.color_diffuse[0] = 0;
                    intersection.color_diffuse[1] = 0;
                    intersection.color_diffuse[2] = 0;
                    intersection.color_specular[0] = 0;
                    intersection.color_specular[1] = 0;
                    intersection.color_specular[2] = 0;
                    intersection.normal[0] = 0;
                    intersection.normal[1] = 0;
                    intersection.normal[2] = 0;
                    intersection.shininess = 0;

                    if (intersected_sphere != -1) {
                        Sphere sphere = spheres[intersected_sphere];
                        // Set color properties
                        memcpy(intersection.color_diffuse, sphere.color_diffuse, 3 * sizeof(double));
                        memcpy(intersection.color_specular, sphere.color_specular, 3 * sizeof(double));
                        intersection.shininess = sphere.shininess;

                        // Calculate intersection normal for the sphere
                        subtract(intersection.position, sphere.position, intersection.normal);
                        normalize(intersection.normal);
                    }
                    else if (intersected_triangle != -1) {

                        Triangle triangle = triangles[intersected_triangle];
                        // Find barycentric coordinates of intersection point
                        double u, v;
                        barycentric_coordinates(triangle, intersection.position, u, v);

                        // Set color properties
                        interpolate_color(triangle, u, v, intersection.color_diffuse, false);
                        interpolate_color(triangle, u, v, intersection.color_specular, true);
                        intersection.shininess = interpolate_shininess(triangle, u, v);

                        // Calculate intersection normal for the triangle
                        interpolate_normal(triangle, u, v, intersection.normal);
                    }


                    // Add ambient light
                    final_color[0] += ambient_light[0] * intersection.color_diffuse[0];
                    final_color[1] += ambient_light[1] * intersection.color_diffuse[1];
                    final_color[2] += ambient_light[2] * intersection.color_diffuse[2];

                    // 2. Determine if intersection point is in shadow and calculate Phong shading
                    for (int light_idx = 0; light_idx < num_lights; ++light_idx)
                    {
                        Light light = lights[light_idx];
                        if (!is_in_shadow(intersection, light))
                        {
                            double light_color[3] = { 0, 0, 0 };

                            // Calculate Phong shading for the point if not in shadow
                            phong_shading(intersection, light, camera_pos, light_color);

                            // Add light contribution to the final color
                            final_color[0] += light_color[0];
                            final_color[1] += light_color[1];
                            final_color[2] += light_color[2];
                        }
                    }

                    // 3. Clamp color to [0, 1]
                    unsigned char r = clamp(final_color[0], 0.0, 1.0) * 255;
                    unsigned char g = clamp(final_color[1], 0.0, 1.0) * 255;
                    unsigned char b = clamp(final_color[2], 0.0, 1.0) * 255;

                    //plot_pixel(x, y, r, g, b);

                    // Add the calculated color to the accumulated color
                    final_color_accum[0] += final_color[0];
                    final_color_accum[1] += final_color[1];
                    final_color_accum[2] += final_color[2];
                }
                else
                {
                    // Set background color if no intersection
                    final_color_accum[0] += 1;
                    final_color_accum[1] += 1;
                    final_color_accum[2] += 1;
                }
            }

            // Calculate the average color
            final_color_accum[0] /= NUM_SAMPLES;
            final_color_accum[1] /= NUM_SAMPLES;
            final_color_accum[2] /= NUM_SAMPLES;

            // Plot the pixel with the averaged color
            plot_pixel(x, y, final_color_accum[0] * 255, final_color_accum[1] * 255, final_color_accum[2] * 255);
        }
    }
    glEnd();
    glFlush();
}


void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
    glColor3f(((double)r) / 256.f, ((double)g) / 256.f, ((double)b) / 256.f);
    glBegin(GL_POINTS);
    glVertex2i(x, y);
    glEnd();
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
    buffer[y][x][0] = r;
    buffer[y][x][1] = g;
    buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
    plot_pixel_display(x, y, r, g, b);
    if (mode == MODE_JPEG)
        plot_pixel_jpeg(x, y, r, g, b);
}

void save_jpg()
{
    printf("Saving JPEG file: %s\n", filename);

    ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
    if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
        printf("Error in Saving\n");
    else
        printf("File saved Successfully\n");
}

void parse_check(const char* expected, char* found)
{
    if (strcasecmp(expected, found))
    {
        printf("Expected '%s ' found '%s '\n", expected, found);
        printf("Parse error, abnormal abortion\n");
        exit(0);
    }
}

void parse_doubles(FILE* file, const char* check, double p[3])
{
    char str[100];
    fscanf(file, "%s", str);
    parse_check(check, str);
    fscanf(file, "%lf %lf %lf", &p[0], &p[1], &p[2]);
    printf("%s %lf %lf %lf\n", check, p[0], p[1], p[2]);
}

void parse_rad(FILE* file, double* r)
{
    char str[100];
    fscanf(file, "%s", str);
    parse_check("rad:", str);
    fscanf(file, "%lf", r);
    printf("rad: %f\n", *r);
}

void parse_shi(FILE* file, double* shi)
{
    char s[100];
    fscanf(file, "%s", s);
    parse_check("shi:", s);
    fscanf(file, "%lf", shi);
    printf("shi: %f\n", *shi);
}

int loadScene(char* argv)
{
    FILE* file = fopen(argv, "r");
    int number_of_objects;
    char type[50];
    Triangle t;
    Sphere s;
    Light l;
    fscanf(file, "%i", &number_of_objects);

    printf("number of objects: %i\n", number_of_objects);

    parse_doubles(file, "amb:", ambient_light);

    for (int i = 0; i < number_of_objects; i++)
    {
        fscanf(file, "%s\n", type);
        printf("%s\n", type);
        if (strcasecmp(type, "triangle") == 0)
        {
            printf("found triangle\n");
            for (int j = 0; j < 3; j++)
            {
                parse_doubles(file, "pos:", t.v[j].position);
                parse_doubles(file, "nor:", t.v[j].normal);
                parse_doubles(file, "dif:", t.v[j].color_diffuse);
                parse_doubles(file, "spe:", t.v[j].color_specular);
                parse_shi(file, &t.v[j].shininess);
            }

            if (num_triangles == MAX_TRIANGLES)
            {
                printf("too many triangles, you should increase MAX_TRIANGLES!\n");
                exit(0);
            }
            triangles[num_triangles++] = t;
        }
        else if (strcasecmp(type, "sphere") == 0)
        {
            printf("found sphere\n");

            parse_doubles(file, "pos:", s.position);
            parse_rad(file, &s.radius);
            parse_doubles(file, "dif:", s.color_diffuse);
            parse_doubles(file, "spe:", s.color_specular);
            parse_shi(file, &s.shininess);

            if (num_spheres == MAX_SPHERES)
            {
                printf("too many spheres, you should increase MAX_SPHERES!\n");
                exit(0);
            }
            spheres[num_spheres++] = s;
        }
        else if (strcasecmp(type, "light") == 0)
        {
            printf("found light\n");
            parse_doubles(file, "pos:", l.position);
            parse_doubles(file, "col:", l.color);

            if (num_lights == MAX_LIGHTS)
            {
                printf("too many lights, you should increase MAX_LIGHTS!\n");
                exit(0);
            }
            lights[num_lights++] = l;
        }
        else
        {
            printf("unknown type in scene description:\n%s\n", type);
            exit(0);
        }
    }
    return 0;
}

void display()
{
}

void init()
{
    glMatrixMode(GL_PROJECTION);
    glOrtho(0, WIDTH, 0, HEIGHT, 1, -1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
    //hack to make it only draw once
    static int once = 0;
    if (!once)
    {
        //use printf to print the value of antialiased
        printf("Anti-aliasing: %d\n", antialiased);
        if (antialiased == 0) {
            printf("Mode: Standard\n");
            draw_scene();
        }
        else if (antialiased == 1) {
            printf("Mode: Anti-aliasing with supersampling!\n");
            draw_antialiased_scene();
        }
        else if (antialiased == 2){
            printf("Mode: Monte-Carlo sampling!\n");
            draw_monte_carlo_scene();
        }
        if (mode == MODE_JPEG)
            save_jpg();
    }
    once = 1;
}

int main(int argc, char** argv)
{
    if ((argc < 2) || (argc > 4))
    {
        printf("Usage: %s <input scenefile> [output jpegname] [optional: 1 if anti-aliased]\n", argv[0]);
        exit(0);
    }
    if (argc >= 3)
    {
        mode = MODE_JPEG;
        filename = argv[2];
        if (argc == 4) {
            // if the third argument is 1 or 2, set antialiased to 1 or 2 respectively
            antialiased = atoi(argv[3]);
        }
    }
    else if (argc == 2)
        mode = MODE_DISPLAY;

    glutInit(&argc, argv);
    loadScene(argv[1]);

    glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(WIDTH, HEIGHT);
    int window = glutCreateWindow("Ray Tracer");
#ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
#endif
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    init();
    glutMainLoop();
}


