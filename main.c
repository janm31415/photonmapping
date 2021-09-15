#include <stdio.h>
#include <math.h>

typedef struct vec3
{
    float x, y, z;
    
} vec3;

float dot(vec3 v1, vec3 v2)
{
    float dp = v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
    return dp;
}

vec3 cross(vec3 v1, vec3 v2)
{
    vec3 cp;
    cp.x = v1.y*v2.z - v1.z*v2.y;
    cp.y = v1.x*v2.z - v1.z*v2.x;
    cp.z = v1.x*v2.y - v1.y*v2.x;
    return cp;
}

vec3 add(vec3 v1, vec3 v2)
{
    vec3 res;
    res.x = v1.x + v2.x;
    res.y = v1.y + v2.y;
    res.z = v1.z + v2.z;
    return res;
}

vec3 sub(vec3 v1, vec3 v2)
{
    vec3 res;
    res.x = v1.x - v2.x;
    res.y = v1.y - v2.y;
    res.z = v1.z - v2.z;
    return res;
}

vec3 mulv(vec3 v1, vec3 v2)
{
    vec3 res;
    res.x = v1.x * v2.x;
    res.y = v1.y * v2.y;
    res.z = v1.z * v2.z;
    return res;
}

vec3 muls(vec3 v1, float a)
{
    vec3 res;
    res.x = v1.x * a;
    res.y = v1.y * a;
    res.z = v1.z * a;
    return res;
}

vec3 smul(float a, vec3 v1)
{
    vec3 res;
    res.x = v1.x * a;
    res.y = v1.y * a;
    res.z = v1.z * a;
    return res;
}

vec3 normalize(vec3 v)
{
    float L = sqrt(dot(v,v));
    return muls(v, 1.f/L);
}


int image_size = 512;
int nr_types = 2; // 2 object types (sphere = 0, plane = 1)
int nr_objects[2] = {2, 5}; // 2 spheres, 5 planes
float g_ambient = 0.1f;
vec3 g_origin = {0,0,0};
vec3 light = {0.f, 1.2f, 3.75f}; // point light source position
float spheres[2][4] = {{1.f, 0.f, 4.f, 0.5f},{-0.6f,-1.f,4.5f,0.5f}}; // sphere centers & radius
float planes[5][2] = {{0,1.5f},{1,-1.5f},{0,-1.5f},{1,1.5f},{2,5.f}}; // plane axis & distance-to-origin
int g_intersect = 0; // global set by latest ray tracing call if anything was intersected
int g_type; // type of the intersected object (sphere or plane)
int g_index; // index of the intersected object
float g_sq_dist = -1.f; // squared distane from ray origin to intersection
float g_dist = -1.f; // distance from ray origin to intersection
vec3 g_point = {0,0,0}; // point where the ray intersected the object
int light_photons = 0;


void check_distance(float dist, int p, int i)
{
    if (dist < g_dist && dist > 0)
    {
        g_type = p;
        g_index = i;
        g_dist = dist;
        g_intersect = 1;
    }
}

void ray_sphere( int idx, vec3 ray, vec3 origin)
{
    vec3 so = {spheres[idx][0], spheres[idx][1], spheres[idx][2]};
    vec3 s = sub(so, origin);
    float radius = spheres[idx][3];
    
    float A = dot(ray, ray);
    float B = -2.f * dot(s, ray);
    float C = dot(s, s) - radius*radius;
    float D = B*B - 4.f*A*C;
    
    if (D > 0)
    {
        float sign = C < -0.00001f ? 1 : -1;
        float dist = (-B + sign*sqrt(D))/(2.f*A);
        check_distance(dist,0,idx);
    }
}

void ray_plane( int idx, float* ray, float* origin)
{
    int axis = (int)planes[idx][0];
    if (ray[axis]!=0.f)
    {
        float dist = (planes[idx][1] - origin[axis]) / ray[axis];
        check_distance(dist, 1, idx);
    }
}

void ray_object(int type, int idx, vec3 ray, vec3 origin)
{
    type == 0 ? ray_sphere(idx, ray, origin) : ray_plane(idx, &ray.x,
                                                         &origin.x);
}

void raytrace(vec3 ray, vec3 origin)
{
    g_intersect = 0;
    g_dist = 99999999.f;
    for (int t = 0; t < nr_types; ++t)
    {
        for (int i = 0; i < nr_objects[t]; ++i)
        {
            ray_object(t,i,ray,origin);
        }
    }
}

vec3 sphere_normal(int idx, vec3 P)
{
    vec3 so = {spheres[idx][0], spheres[idx][1], spheres[idx][2]};
    return normalize(sub(P, so));
}

vec3 plane_normal(int idx, vec3 P, float* O)
{
    int axis = (int)planes[idx][0];
    float N[3] = {0,0,0};
    N[axis] = O[axis] - planes[idx][1];
    vec3 n = {N[0], N[1], N[2]};
    return normalize(n);
}

vec3 surface_normal(int type, int index, vec3 P, vec3 Inside)
{
    if (type == 0)
        return sphere_normal(index, P);
    else
        return plane_normal(index, P, &Inside.x);
}

vec3 reflect(vec3 ray, vec3 from_point)
{
    vec3 N = surface_normal(g_type, g_index, g_point, from_point);
    return normalize(sub(ray, muls(N, 2.f*dot(ray,N))));
}

float light_diffuse(vec3 N, vec3 P)
{
    vec3 L = normalize(sub(light, P));
    return dot(N, L);
}

float minimum(float a, float b)
{
    return a < b ? a : b;
}

float maximum(float a, float b)
{
    return a < b?b:a;
}

float light_object(int type, int idx, vec3 P, float ambient)
{
    float i = light_diffuse(surface_normal(type, idx, P, light), P);
    return minimum(1.f, maximum(i, ambient));
}

vec3 filter_color(vec3 rgb, float r, float g, float b)
{
    vec3 out = {r, g, b};
    out.x = minimum(out.x, rgb.x);
    out.y = minimum(out.y, rgb.y);
    out.z = minimum(out.z, rgb.z);
        
    return out;
}

vec3 get_color(vec3 rgb, int type, int index)
{
    if (type == 1 && index == 0)
        return filter_color(rgb, 0.f, 1.f, 0.1);
    else if (type == 1 && index == 2)
        return filter_color(rgb, 1.f, 0.f, 0.f);
    else
        return filter_color(rgb, 1.f, 1.f, 1.f);
}

vec3 compute_pixel_color(int x, int y)
{
    vec3 clr = {0,0,0};
    vec3 ray = {(float)x/(float)image_size - 0.5f, -((float)y/(float)image_size - 0.5f), 1.f};
    
    raytrace(ray, g_origin);
    
    if (g_intersect)
    {
        g_point = muls(ray,g_dist);
        if (g_type == 0  && g_index == 1) // mirror surface
        {
            ray = reflect(ray, g_origin);
            raytrace(ray, g_point);
            if (g_intersect)
            {
                g_point = add(muls(ray, g_dist), g_point);
            }
        }
        if (light_photons)
        {
            
        } else
        {
            int type = g_type;
            int index = g_index;
            float i = g_ambient;
            raytrace(sub(g_point, light), light);
            if (type == g_type && index == g_index)
                i = light_object(g_type, g_index, g_point, g_ambient);
            clr.x = i;
            clr.y = i;
            clr.z = i;
            clr = get_color(clr, type, index);
        }
    }
    
    return clr;
}

int main()
{
    printf("P6 %d %d 255\n", image_size, image_size);
    for (int y = 0; y < image_size; ++y)
    {
        for (int x = 0; x < image_size; ++x)
        {
            vec3 clr = muls(compute_pixel_color(x,y), 255.f);
            printf("%c%c%c", (int)clr.x, (int)clr.y, (int)clr.z);
        }
    }
    return 0;
}
