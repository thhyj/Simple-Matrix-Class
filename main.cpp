#include <cmath>
#include <iostream>
#include <bits/stdc++.h>




/**
 * 矩阵类
 */
struct Matrix {
    constexpr static const float eps = 1e-6;
    constexpr const static float pi = 3.1415926;
    std::vector<std::vector<float>>mat;
    int n, m;

    Matrix() = default;
    /**
     * 生成n*m的空矩阵
     * @param n 行数
     * @param m 列数
     */
    Matrix(int n, int m):n(n), m(m){
        mat.resize(n);
        for(int i = 0; i < n; ++i)
            mat[i].resize(m);
    }
    /**
     * 以行优先的方式输入参数构造矩阵
     * @tparam T 可变参数模板
     * @param n 行数
     * @param m 列数
     * @param a 矩阵元素，行优先
     */
    template<typename... T>
    Matrix(int n, int m, T... a): n(n), m(m){
        float parameter[16], *pointer;
        pointer = parameter;
        for(auto temp :{a...}) {
            *pointer = temp;
            ++pointer;
        }
        *this = Matrix(n, m, parameter);
    }
    /**
     * 以行优先的一维float数组构造矩阵
     * @param n 行数
     * @param m 列数
     * @param a 行优先的一维数组
     */
    Matrix(int n, int m, float *a): n(n), m(m){
        mat.resize(n);

        for(int i = 0; i < n; ++i) {
            mat[i].resize(m);
            for(int j = 0; j < m; ++j) {
                mat[i][j] = *a;
                ++a;
            }
        }
    }
    /**
     * 矩阵取负
     * @return
     */
    const Matrix operator -() const {
        Matrix matrix = Matrix(n, m);
        for(int i =0 ;i  < n; ++i) {
            for(int j = 0; j < m; ++j) {
                matrix.mat[i][j] = -mat[i][j];
            }
        }
        return matrix;
    }
    /**
     * 矩阵减法
     * @param rhs
     * @return
     */
    const Matrix operator -(const Matrix& rhs) const {
        if(n != rhs.n || m != rhs.m) {
            puts("invalid");
            return Matrix(1,1);
        }
        Matrix matrix = Matrix(n, m);
        for(int i =0 ;i  < n; ++i) {
            for(int j = 0; j < m; ++j) {
                matrix.mat[i][j] = mat[i][j] - rhs.mat[i][j];
            }
        }
        return matrix;
    }

    /**
     * 矩阵加法
     * @param rhs
     * @return
     */
    const Matrix operator +(const Matrix&  rhs) const {
        if(n != rhs.n || m != rhs.m) {
            puts("invalid");
            return Matrix(1,1);
        }
        Matrix matrix = Matrix(n, m);
        for(int i =0 ;i  < n; ++i) {
            for(int j = 0; j < m; ++j) {
                matrix.mat[i][j] = mat[i][j] - rhs.mat[i][j];
            }
        }
        return matrix;
    }
    /**
     * 矩阵乘法
     * 采用了基本的ikj取指优化和稀疏矩阵优化
     * @param rhs
     * @return 相乘得到的矩阵
     */
    const Matrix operator * (const Matrix& rhs) const {
        if(m != rhs.n) {
            puts("两个矩阵无法相乘");
            return Matrix(1,1);
        }
        const std::vector<std::vector<float>> &a = mat;
        const std::vector<std::vector<float>> &b = (rhs.mat);
        Matrix temp = Matrix(n, rhs.m);
        float tempValue;
        std::vector<float>::const_iterator st;


        for(int i = 0; i < n; ++i) {
            for(int k = 0; k < m; ++k) {
                if(fabs(a[i][k]) > eps) {
                    tempValue = a[i][k];
                    st = b[k].begin();
                    for(int j = 0; j < rhs.m; ++j) {
                        temp.mat[i][j] += tempValue * (*st);
                        ++st;
                    }
                }
            }
        }
        return temp;
    }

    /**
     * 打印矩阵
     */
    void print() const {
        puts("--------------------------print Matrix--------------------------");
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                printf("%f ",  mat[i][j]);
            }
            puts("");
        }
        puts("---------------------------------------------------------------");
    }
    /**
     * 将矩阵转换成opengl格式的一维列优先数组
     * @param target 用于接收矩阵的float数组
     */
    void getOpenglModeMatrix(float * target) const {
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                *target = mat[j][i];
                ++target;
            }
        }
    }
    /**
     * 矩阵转置
     * @return 返回转置的矩阵
     */
    const Matrix transpose() const {
        Matrix matrix = Matrix(m, n);
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                matrix.mat[j][i] = mat[i][j];
            }
        }
        return matrix;
    }

    /**
     * 生成单位矩阵
     * @param n 单位矩阵的维数
     * @return 生成的单位矩阵
     */

    static const Matrix getIdentityMatrix(int n) {
        Matrix matrix = Matrix(n, n);
        for(int i = 0; i < n; ++i) {
            matrix.mat[i][i] = 1;
        }
        return matrix;
    }

    /**
     * 生成三维坐标的位移矩阵，由于采用齐次坐标，实际生成四维矩阵
     * @param x x方向位移
     * @param y y方向位移
     * @param z z方向位移
     * @return 位移矩阵
     */

    static const Matrix getTranslationMatrix3(float x, float y, float z) {
        Matrix matrix = getIdentityMatrix(4);
        matrix.mat[0][3] = x;
        matrix.mat[1][3] = y;
        matrix.mat[2][3] = z;
        return matrix;
    }

    /**
     * 生成三维坐标的位移矩阵，由于采用齐次坐标，实际生成四维矩阵
     * 每次能且仅能绕一个轴旋转
     * @param angle 旋转角度,采用角度制
     * @param x 为1则绕x轴旋转
     * @param y 为1则绕y轴旋转
     * @param z 为1则绕z轴旋转
     * @return 旋转矩阵
     */
    static const Matrix getRotationMatrix3(float angle, int x, int y, int z) {
        if((x == 1) + (y == 1) + (z == 1) != 1) {
            puts("旋转轴不合法");
            return Matrix(1,1);
        }
        Matrix matrix = getIdentityMatrix(4);
        angle = angle / 180 * pi;
        if(z == 1) {
            matrix.mat[0][0] = std::cos(angle);
            matrix.mat[1][0] = std::sin(angle);
            matrix.mat[0][1] = -std::sin(angle);
            matrix.mat[1][1] = std::cos(angle);
        } else {
            if(y == 1) {
                matrix.mat[0][0] = std::cos(angle);
                matrix.mat[2][0] = -std::sin(angle);
                matrix.mat[0][2] = std::sin(angle);
                matrix.mat[2][2] = std::cos(angle);
            } else {
                if(x == 1) {
                    matrix.mat[1][1] = std::cos(angle);
                    matrix.mat[2][1] = std::sin(angle);
                    matrix.mat[1][2] = -std::sin(angle);
                    matrix.mat[2][2] = std::cos(angle);
                }
            }
        }
        return matrix;
    }

    /**
     * 生成投影矩阵
     * 由于opengl中只显示NDC立方体范围内的内容，且opengl坐标系为右手系，而NDC坐标系为左手系
     * 所以这个矩阵做了投影，坐标轴转换，归一化三件事，
     * 因此这个矩阵和课件上的完全不一样
     * 为了搞清楚opengl的几种坐标系的转换，最后推出这个矩阵，还是耗费了不少时间的
     * 毕竟课件上只讲了将一个物体投影到一个二维平面上的方式，而opengl的NDC坐标系是一个以(0,0,0)为中心的长为1的立方体
     * 虽然可以投影到二维平面后用opengl自带的正交投影工具，但是我还是想了解更底层的原理
     * 这个矩阵假定了视点在原点，视景体在z轴负方向上
     * @param n 视景体近平面离原点的距离
     * @param f 视景体远平面离原点的距离
     * @param l 视景体近平面x坐标最小值
     * @param r 视景体近平面x坐标最大值
     * @param b 视景体近平面y坐标最小值
     * @param t 视景体近平面y坐标最大值
     * @return 返回投影矩阵
     */

    static const Matrix getProjectionMatrix3(float l, float r, float b, float t, float n, float f) {
        Matrix matrix = getIdentityMatrix(4);
        matrix.mat[0][0] = (2 * n) / (r - l);
        matrix.mat[1][1] = (2 * n) / (t - b);
        matrix.mat[0][2] = (r + l) / (r - l);
        matrix.mat[1][2] = (t + b) / (t - b);
        matrix.mat[2][2] = -(f + n) / (f - n);
        matrix.mat[2][3] = -(f * n * 2) / (f - n);
        matrix.mat[3][2] = -1;
        return matrix;
    }
};


int main() {
    Matrix::getIdentityMatrix(5).print();
    Matrix::getTranslationMatrix3(3,4,5).print();
    Matrix cube = Matrix(8, 4,
                         0.0f, 0.0f, 0.0f,1.0f,
                         0.0f, 30.0f, 0.0f,1.0f,
                         0.0f,  0.0f, -25.0f,1.0f,
                         0.0f, 30.0f, -25.0f,1.0f,
                         20.0f, 0.0f, 0.0f,1.0f,
                         20.0f, 30.0f, 0.0f,1.0f,
                         20.0f,  0.0f, -25.0f,1.0f,
                         20.0f, 30.0f, -25.0f,1.0f);
    cube = cube.transpose();
    cube.print();
    cube = Matrix::getTranslationMatrix3(3,4,5) * cube;
    cube.print();
    Matrix::getRotationMatrix3(45, 1,0,0).print();
    Matrix::getRotationMatrix3(45, 0,1,0).print();
    (Matrix::getRotationMatrix3(45, 0,0,1) + Matrix::getRotationMatrix3(45, 0,1,0)).print();
    return 0;
}
