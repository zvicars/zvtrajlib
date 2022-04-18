//implementation of the Moller Trombore algorithm modified from wikipedia
//https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
#include "../meshtools.hpp"
#include "../stlmath.hpp"
#define EPSILON 0.0000001
bool MollerTrombore(std::array<double,3> rayOrigin, 
                           std::array<double,3> rayVector, 
                           std::array<std::array<double,3>, 3> inTriangle,
                           std::array<double,3>& outIntersectionPoint)
{
    rayVector = rayVector * (1.0/norm2(rayVector));
    std::array<double,3> vertex0 = inTriangle[0];
    std::array<double,3> vertex1 = inTriangle[1];  
    std::array<double,3> vertex2 = inTriangle[2];
    std::array<double,3> edge1, edge2, h, s, q;
    double a,f,u,v;
    edge1 = vertex1 - vertex0;
    edge2 = vertex2 - vertex0;
    h = cross(rayVector, edge2);
    a = dot(edge1, h);
    if (a > -EPSILON && a < EPSILON)
        return 0;    // This ray is parallel to this triangle.
    f = 1.0/a;
    s = rayOrigin - vertex0;
    u = f * dot(s,h);
    if (u < 0.0 || u > 1.0)
        return 0;
    q = cross(s,edge1);
    v = f * dot(rayVector, q);
    if (v < 0.0 || u + v > 1.0)
        return 0;
    // At this stage we can compute t to find out where the intersection point is on the line.
    double t = f * dot(edge2, q);
    if (t > EPSILON) // ray intersection
    {
        outIntersectionPoint = rayOrigin + rayVector*t;
        return 1;
    }
    else // This means that there is a line intersection but not a ray intersection.
        return 0;
}