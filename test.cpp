#include<cstdio>
#include<vector>
#include<iostream>
#include"Bspline/Bspline3D.hpp"
#define eps 1e-7

int main() 
{
 //   std::vector<Point3D> controlPoints = {
 //       {0.25212955, 0.20192656, 0.73066454},
 //       {0.07040837, 0.04212473, 0.14618082},
 //       {0.80285129, 0.90044722, 0.02047026},
 //       {0.54550146, 0.31983326, 0.44662081},
 //       {0.91162146, 0.9554301,  0.72781836},
 //       {0.0466265,  0.8317827,  0.25163542},
 //       {0.58953454, 0.80531284, 0.4129208},
 //       {0.50378757, 0.55654583, 0.46313412},
 //       {0.43309917, 0.5984493,  0.49776652},
 //       {0.67352757, 0.68633497, 0.47716614}
 //   };

 //   std::vector<double> knots = { 0.         ,0.       ,  0.        , 0.     ,    0.1634004  ,0.29624581,
 //0.30784826, 0.61157454 ,0.62200473 ,0.83257783 ,1.,        1.,
 //1. ,        1. };

 //   BSplineCurve3D curve(3, controlPoints, knots);
 //   for (double i = 0; i <= 1-eps ; i += 0.01)
 //   {
 //       Point3D p = curve.evaluate(i);
 //       printf("%f %f %f\n", p.x, p.y, p.z);
 //   }

    std::vector<std::vector<Point3D>> controlPoints = {
    {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}},
    {{0, 1, 0}, {1, 1, 1}, {2, 1, 0}},
    {{0, 2, 0}, {1, 2, 0}, {2, 2, 0}}
    };

    // 定义节点向量  
    std::vector<double> knotsU = { 0, 0, 0, 1, 1, 1 };
    std::vector<double> knotsV = { 0, 0, 0, 1, 1, 1 };

    // 创建 B 样条曲面对象  
    BSplineSurface3D surface(2, 2, controlPoints, knotsU, knotsV);

    // 测试曲面上的一些点  
    Point3D point1 = surface.evaluate(0.5, 0.5);
    Point3D point2 = surface.evaluate(0.3, 0.7);
    Point3D point3 = surface.evaluate(1.0, 0.0);

    // 打印结果  
    std::cout << "Point 1: " << point1.x << ", " << point1.y << ", " << point1.z << std::endl;
    std::cout << "Point 2: " << point2.x << ", " << point2.y << ", " << point2.z << std::endl;
    std::cout << "Point 3: " << point3.x << ", " << point3.y << ", " << point3.z << std::endl;

    for (double i = 0; i <= 1; i += 0.1)
    {
        for (double j = 0; j <= 1; j += 0.1)
        {
            Point3D p = surface.evaluate(i, j);
            std::cout << "i: " << i << " j: " << j <<" 坐标: " << p.x << ", " << p.y << ", " << p.z << std::endl;
        }
    }

    return 0;
}