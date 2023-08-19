#include<cmath>
#include<vector>
struct Point3D 
{
    double x;
    double y;
    double z;
};

class BezierCurve3D
{
public:
    explicit BezierCurve3D(const std::vector<Point3D>& controlPoints) : controlPoints_(controlPoints)
    {
        n_ = controlPoints.size() - 1;
    }

    inline Point3D evaluate(double t) const
    {
        Point3D result = { 0.0, 0.0, 0.0 };

        for (int i = 0; i <= n_; i++)
        {
            double factor = Bernstein(n_, i, t);
            result.x += controlPoints_[i].x * factor;
            result.y += controlPoints_[i].y * factor;
            result.z += controlPoints_[i].z * factor;
        }
        return result;
    }

private:
    inline int binomialCoefficient(int n, int k) const
    {
        if (k > n || n <= 0)
            return 0;
        if (k == 0 || k == n)
            return 1;
        else
            return binomialCoefficient(n - 1, k - 1) + binomialCoefficient(n - 1, k);
    }

    inline double Bernstein(int n, int i, double t) const
    {
        return (double)binomialCoefficient(n, i) * pow(t, i) * pow(1 - t, n - i);
    }

private:
    std::vector<Point3D> controlPoints_;
    int n_; // ´ÎÊý,control points - 1  
};