#include<vector>
#define eps 1e-7

struct Point2D
{
    double x;
    double y;
    double z;
};

class BSplineCurve2D
{
public:
    explicit BSplineCurve2D(const int& degree, const std::vector<Point2D>& controlPoints, const std::vector<double>& knots) : degree_(degree), controlPoints_(controlPoints), knots_(knots)
    {
        n_ = controlPoints.size() - 1;
    }

    inline Point2D evaluate(double t) const
    {
        int span = findSpan(t);
        std::vector<double> basis = computeBasis(span, t);
        Point2D result = { 0.0, 0.0 };

        for (int i = 0; i <= degree_; ++i)
        {
            double basisMultiplier = basis[i];
            result.x += controlPoints_[span - degree_ + i].x * basisMultiplier;
            result.y += controlPoints_[span - degree_ + i].y * basisMultiplier;
        }

        return result;
    }

private:
    inline int findSpan(double t) const
    {
        int low = degree_;
        int high = n_ + 1;
        int mid;

        while (high - low > 1)
        {
            mid = (high + low) >> 1;

            if (t < knots_[mid] - eps)
                high = mid;
            else
                low = mid;
        }

        return low;
    }

    inline std::vector<double> computeBasis(int span, double t) const
    {
        std::vector<double> basis(degree_ + 1, 0.0);
        std::vector<double> left(degree_ + 1, 0.0);
        std::vector<double> right(degree_ + 1, 0.0);

        basis[0] = 1.0;

        for (int j = 1; j <= degree_; ++j)
        {
            left[j] = t - knots_[span + 1 - j];
            right[j] = knots_[span + j] - t;

            double saved = 0.0;

            for (int r = 0; r < j; ++r)
            {
                double basisMultiplier = basis[r] / (right[r + 1] + left[j - r]);
                basis[r] = saved + right[r + 1] * basisMultiplier;
                saved = left[j - r] * basisMultiplier;
            }

            basis[j] = saved;
        }

        return basis;
    }

private:
    std::vector<Point2D> controlPoints_;
    int n_; // control points - 1  
    int degree_;
    std::vector<double> knots_;
};
