#include <cmath>
#include <iostream>
#include <bits/stdc++.h>




struct Matrix {
    constexpr static const float eps = 1e-6;
    constexpr const static float pi = 3.1415926;
   // float mat[MatrixSize][MatrixSize] = {0};
    std::vector<std::vector<float>>mat;
    int n, m;

    Matrix() = default;

    Matrix(int n, int m):n(n), m(m){
        mat.resize(n);
        for(int i = 0; i < n; ++i)
            mat[i].resize(m);
    }

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

    Matrix operator * (Matrix &another) {
        if(m != another.n) {
            puts("两个矩阵无法相乘");
            return Matrix(1,1);
        }
        std::vector<std::vector<float>> &a = mat;
        std::vector<std::vector<float>> &b = (another.mat);
        Matrix temp = Matrix(n, another.m);
        float tempValue;
        std::vector<float>::iterator st;


        for(int i = 0; i < n; ++i) {
            for(int k = 0; k < m; ++k) {
                if(fabs(a[i][k]) > eps) {
                    tempValue = a[i][k];
                    st = b[k].begin();
                     for(int j = 0; j < another.m; ++j) {
                        temp.mat[i][j] += tempValue * (*st);
                        ++st;
                    }
                }
            }
        }
        return temp;
    }


    void print() {
        puts("--------------------------print Matrix--------------------------");
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                printf("%f ",  mat[i][j]);
            }
            puts("");
        }
        puts("---------------------------------------------------------------");
    }

    void getOpenglModeMatrix(float * target) {
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                *target = mat[j][i];
                ++target;
            }
        }
    }

    Matrix transpose() {
        Matrix matrix = Matrix(m, n);
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                matrix.mat[j][i] = mat[i][j];
            }
        }
        return matrix;
    }

    static Matrix getIdentityMatrix(int n) {
        Matrix matrix = Matrix(n, n);
        for(int i = 0; i < n; ++i) {
            matrix.mat[i][i] = 1;
        }
        return matrix;
    }

    static Matrix getTranslationMatrix3(float x, float y, float z) {
        Matrix matrix = getIdentityMatrix(4);
        matrix.mat[0][3] = x;
        matrix.mat[1][3] = y;
        matrix.mat[2][3] = z;
        return matrix;
    }

    static Matrix getRotationMatrix3(float angle, int x, int y, int z) {
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
                matrix.mat[2][0] = std::sin(angle);
                matrix.mat[0][2] = -std::sin(angle);
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
    Matrix::getRotationMatrix3(45, 0,0,1).print();
    return 0;
}
