#include<vector>
#define eps 1e-7

struct Point3D
{
    double x;
    double y;
    double z;
};

class BSplineCurve3D
{
public:
    explicit BSplineCurve3D(const int& degree, const std::vector<Point3D>& controlPoints, const std::vector<double>& knots) : degree_(degree), controlPoints_(controlPoints), knots_(knots)
    {
        n_ = controlPoints.size() - 1;
    }

    inline Point3D evaluate(double t) const
    {
        int span = findSpan(t);
        std::vector<double> basis = computeBasis(span, t);
        Point3D result = { 0.0, 0.0, 0.0 };

        for (int i = 0; i <= degree_; ++i)
        {
            double basisMultiplier = basis[i];
            result.x += controlPoints_[span - degree_ + i].x * basisMultiplier;
            result.y += controlPoints_[span - degree_ + i].y * basisMultiplier;
            result.z += controlPoints_[span - degree_ + i].z * basisMultiplier;
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
    std::vector<Point3D> controlPoints_;
    int n_; // control points - 1  
    int degree_;
    std::vector<double> knots_;
};

class BSplineSurface3D
{
public:
    explicit BSplineSurface3D(const int& degreeU, const int& degreeV, const std::vector<std::vector<Point3D>>& controlPoints, const std::vector<double>& knotsU, const std::vector<double>& knotsV)
        : degreeU_(degreeU), degreeV_(degreeV), controlPoints(controlPoints), knotsU_(knotsU), knotsV_(knotsV)
    {
        nU_ = controlPoints.size() - 1;
        nV_ = controlPoints[0].size() - 1;
    }

    inline Point3D evaluate(double u, double v) const
    {
        int spanU = findSpan(u, nU_, degreeU_, knotsU_);
        int spanV = findSpan(v, nV_, degreeV_, knotsV_);
        std::vector<double> basisU = computeBasis(spanU, u, nU_, degreeU_, knotsU_);
        std::vector<double> basisV = computeBasis(spanV, v, nV_, degreeV_, knotsV_);
        Point3D result = { 0.0, 0.0, 0.0 };

        for (int i = 0; i <= degreeU_; ++i)
        {
            for (int j = 0; j <= degreeV_; ++j)
            {
                double basisMultiplier = basisU[i] * basisV[j];
                result.x += controlPoints[spanU - degreeU_ + i][spanV - degreeV_ + j].x * basisMultiplier;
                result.y += controlPoints[spanU - degreeU_ + i][spanV - degreeV_ + j].y * basisMultiplier;
                result.z += controlPoints[spanU - degreeU_ + i][spanV - degreeV_ + j].z * basisMultiplier;
            }
        }

        return result;
    }

private:
    int findSpan(double t, int n, int degree, const std::vector<double>& knots) const
    {
        int low = degree;
        int high = n + 1;
        int mid;

        while (high - low > 1)
        {
            mid = (high + low) >> 1;

            if (t < knots[mid] - eps)
                high = mid;
            else
                low = mid;
        }

        return low;
    }

    std::vector<double> computeBasis(int span, double t, int n, int degree, const std::vector<double>& knots) const
    {
        std::vector<double> basis(degree + 1, 0.0);
        std::vector<double> left(degree + 1, 0.0);
        std::vector<double> right(degree + 1, 0.0);

        basis[0] = 1.0;

        for (int j = 1; j <= degree; ++j)
        {
            left[j] = t - knots[span + 1 - j];
            right[j] = knots[span + j] - t;

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
    std::vector<std::vector<Point3D>> controlPoints;
    int nU_, nV_; // control points - 1    
    int degreeU_, degreeV_;
    std::vector<double> knotsU_;
    std::vector<double> knotsV_;
};
