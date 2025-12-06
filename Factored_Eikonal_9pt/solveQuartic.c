#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <time.h>

// 求解二次方程
int solveQuadratic(double a, double b, double c, double** roots) {
    if (fabs(a) < DBL_EPSILON) {
        if (fabs(b) < DBL_EPSILON) {
            *roots = NULL;
            return 0;
        }
        *roots = (double*)malloc(sizeof(double));
        if (*roots == NULL) {
            return -1;
        }
        (*roots)[0] = -c / b;
        return 1;
    }
    
    double discriminant = b * b - 4 * a * c;
    
    if (discriminant < -DBL_EPSILON) {
        *roots = NULL;
        return 0;
    } else if (fabs(discriminant) < DBL_EPSILON) {
        *roots = (double*)malloc(sizeof(double));
        if (*roots == NULL) {
            return -1;
        }
        (*roots)[0] = -b / (2 * a);
        return 1;
    } else {
        *roots = (double*)malloc(2 * sizeof(double));
        if (*roots == NULL) {
            return -1;
        }
        double sqrtDiscriminant = sqrt(discriminant);
        (*roots)[0] = (-b - sqrtDiscriminant) / (2 * a);
        (*roots)[1] = (-b + sqrtDiscriminant) / (2 * a);
        return 2;
    }
}

// 改进的三次方程求解器 - 专门处理重根情况
int solveCubic(double a, double b, double c, double d, double** roots) {
    if (fabs(a) < DBL_EPSILON) {
        return solveQuadratic(b, c, d, roots);
    }
    
    // 归一化系数
    double b_norm = b / a;
    double c_norm = c / a;
    double d_norm = d / a;
    
    double q = (3.0 * c_norm - b_norm * b_norm) / 9.0;
    double r = (9.0 * b_norm * c_norm - 27.0 * d_norm - 2.0 * b_norm * b_norm * b_norm) / 54.0;
    double discriminant = q * q * q + r * r;
    
    // 使用更大的容差来检测重根
    const double tolerance = 1e-6;
    
    if (discriminant > tolerance) {
        *roots = (double*)malloc(sizeof(double));
        if (*roots == NULL) return -1;
        double s = r + sqrt(discriminant);
        double t = r - sqrt(discriminant);
        s = (s >= 0) ? pow(s, 1.0/3.0) : -pow(-s, 1.0/3.0);
        t = (t >= 0) ? pow(t, 1.0/3.0) : -pow(-t, 1.0/3.0);
        (*roots)[0] = s + t - b_norm / 3.0;
        return 1;
    } else if (fabs(discriminant) <= tolerance) {
        // 处理重根情况
        double sqrtDiscriminant = cbrt(fabs(r));
        if (r < 0) sqrtDiscriminant = -sqrtDiscriminant;
        
        *roots = (double*)malloc(3 * sizeof(double));
        if (*roots == NULL) return -1;
        
        (*roots)[0] = 2.0 * sqrtDiscriminant - b_norm / 3.0;
        (*roots)[1] = -sqrtDiscriminant - b_norm / 3.0;
        (*roots)[2] = -sqrtDiscriminant - b_norm / 3.0; // 重根
        
        return 3;
    } else {
        *roots = (double*)malloc(3 * sizeof(double));
        if (*roots == NULL) return -1;
        double theta = acos(r / sqrt(-q * q * q));
        double sqrtQ = sqrt(-q);
        (*roots)[0] = 2.0 * sqrtQ * cos(theta / 3.0) - b_norm / 3.0;
        (*roots)[1] = 2.0 * sqrtQ * cos((theta + 2.0 * M_PI) / 3.0) - b_norm / 3.0;
        (*roots)[2] = 2.0 * sqrtQ * cos((theta + 4.0 * M_PI) / 3.0) - b_norm / 3.0;
        return 3;
    }
}

// 求解四次方程的主函数
void solveQuartic(double a, double b, double c, double d, double e, double* roots, int* num) {
    for (int i = 0; i < 4; i++) roots[i] = 0.0;
    *num = 0;
    
    if (fabs(a) < DBL_EPSILON) {
        double* cubicRoots = NULL;
        int count = solveCubic(b, c, d, e, &cubicRoots);
        if (count > 0 && cubicRoots != NULL) {
            for (int i = 0; i < count && i < 4; i++) {
                roots[*num] = cubicRoots[i];
                (*num)++;
            }
            free(cubicRoots);
        }
        return;
    }
    
    // 归一化系数
    double normB = b / a;
    double normC = c / a;
    double normD = d / a;
    double normE = e / a;
    
    // 求解辅助三次方程
    double yCoeffs[4] = {
        1.0,
        -normC,
        normB * normD - 4 * normE,
        -normB * normB * normE + 4 * normC * normE - normD * normD
    };
    
    double* yRoots = NULL;
    int yCount = solveCubic(yCoeffs[0], yCoeffs[1], yCoeffs[2], yCoeffs[3], &yRoots);
    
   
    
    if (yCount <= 0 || yRoots == NULL) {
        if (yRoots != NULL) free(yRoots);
        return;
    }
    
    // 选择最大的实数y根
    double y = yRoots[0];
    for (int i = 1; i < yCount; i++) {
        if (yRoots[i] > y) {
            y = yRoots[i];
        }
    }

    
    // 计算r和term2
    double r_squared = normB * normB / 4.0 - normC + y;
    double r = (r_squared >= 0) ? sqrt(r_squared) : 0;
    
    double term2 = normB * y / 2.0 - normD;
    
    // 处理r接近0的情况
    if (fabs(r) < 1e-10) {
        // 当r接近0时，使用替代方法
        r = 0;
        term2 = 0;
    } else {
        term2 /= r;
    }
    
  
    
    // 构造并求解两个二次方程
    double* roots1 = NULL;
    double* roots2 = NULL;
    
    int count1 = solveQuadratic(1.0, normB / 2.0 + r, y / 2.0 + term2 / 2.0, &roots1);
    int count2 = solveQuadratic(1.0, normB / 2.0 - r, y / 2.0 - term2 / 2.0, &roots2);
    
    // 合并结果
    if (count1 > 0 && roots1 != NULL) {
        for (int i = 0; i < count1 && *num < 4; i++) {
            roots[*num] = roots1[i];
            (*num)++;
        }
        free(roots1);
    }
    
    if (count2 > 0 && roots2 != NULL) {
        for (int i = 0; i < count2 && *num < 4; i++) {
            roots[*num] = roots2[i];
            (*num)++;
        }
        free(roots2);
    }
    
    if (yRoots != NULL) free(yRoots);
}
