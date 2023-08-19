#include<cmath>
#include<vector>
struct Point2D
{
    double x;
    double y;
};

class BezierCurve2D
{
public:
    explicit BezierCurve2D(const std::vector<Point2D>& controlPoints) : controlPoints_(controlPoints)
    {
        n_ = controlPoints.size() - 1;
    }

    inline Point2D evaluate(double t) const
    {
        Point2D result = { 0.0, 0.0 };

        for (int i = 0; i <= n_; i++)
        {
            double factor = Bernstein(n_, i, t);
            result.x += controlPoints_[i].x * factor;
            result.y += controlPoints_[i].y * factor;
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
    std::vector<Point2D> controlPoints_;
    int n_; // ´ÎÊý,control points - 1  
};