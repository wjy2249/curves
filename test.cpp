#include<cstdio>
#include<vector>
#include"NURBS3D/NURBS3D.hpp"
#define eps 1e-7

int main() 
{
    std::vector<Point3D> controlPoints = {
        { 0.0, 0.0, 0.0 },
        { 1.0, 2.0, 3.0 },
        { 3.0, 2.0, 1.0 },
        { 4.0, 0.0, 0.0 }
    };

    std::vector<double> knots = { 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 };

    std::vector<double> weights = { 1.0, 2.0, 2.0, 1.0 };

    NURBSCurve3D curve(2, controlPoints, knots, weights);
    for (double i = 0; i < 1 - eps; i += 0.01)
    {
        Point3D p = curve.evaluate(i);
        printf("%f %f %f\n", p.x, p.y, p.z);
    }
    return 0;
}