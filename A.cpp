// Alexey Tkachenko
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <algorithm>

using namespace std;

class Matrix {
private:
    vector<vector<double>> data;

public:
    const int n, m;

    Matrix(int n, int m) : n(n), m(m) {
        data.resize(n, vector<double>(m, 0));
    }

    double &operator()(int i, int j) {
        return data[i][j];
    }

    double operator()(int i, int j) const {
        return data[i][j];
    }

    friend istream& operator>>(istream& in, Matrix &mat) {
        for (int i = 0; i < mat.n; i++)
            for (int j = 0; j < mat.m; j++)
                in >> mat(i, j);

        return in;
    }

    friend ostream& operator<<(ostream& out, Matrix &mat) {
        for (int i = 0; i < mat.n; i++) {
            for (int j = 0; j < mat.m; j++) {
                out << mat(i, j);
                if (j != mat.m - 1) out << " ";
            }
            out << endl;
        }
        return out;
    }

    friend Matrix operator*(const Matrix &A, const Matrix &B) {
        int n = A.n;
        int m = A.m;
        int p = B.m;
        Matrix C(n, p);

        for (int i = 0; i < n; ++i)
            for (int j = 0; j < p; ++j)
                for (int k = 0; k < m; ++k)
                    C(i,j) += A(i,k) * B(k,j);

        return C;
    }

    Matrix transpose() {
        Matrix B(m, n);

        for (int i = 0; i < m; ++i)
            for (int j = 0; j < n; ++j)
                B(i,j) = (*this)(j,i);

        return B;
    }

    Matrix inverse() {
        Matrix A = (*this);
        Matrix I(n,n);

        for (int i = 0; i < n; ++i) {
            I(i,i) = 1;
        }

        for (int i = 0; i < n; ++i) {
            double pivot = A(i,i);
            for (int j = 0; j < n; ++j) {
                A(i,j) /= pivot;
                I(i,j) /= pivot;
            }

            for (int j = 0; j < n; ++j) {
                if (i == j) continue;
                double factor = A(j,i);
                for (int k = 0; k < n; ++k) {
                    A(j,k) -= factor * A(i,k);
                    I(j,k) -= factor * I(i,k);
                }
            }
        }

        return I;
    }
};


int main() {
    int m, n;
    cin >> m;
    vector<double> t(m), b(m);
    for (int i = 0; i < m; ++i) {
        cin >> t[i] >> b[i];
    }
    cin >> n;

    Matrix A(m, n + 1);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j <= n; ++j) {
            A(i,j) = pow(t[i], j);
        }
    }

    Matrix A_T = A.transpose();
    Matrix A_T_A = A_T * A;
    Matrix A_T_A_inv = A_T_A.inverse();

    Matrix A_T_b(m, 1);
    for (int i = 0; i < m; ++i)
        A_T_b(i,0) = b[i];

    Matrix A_T_b_product = A_T * A_T_b;
    Matrix x = A_T_A_inv * A_T_b_product;

    cout << setprecision(4) << fixed;
    cout << "A:" << endl << A;
    cout << "A_T*A:" << endl << A_T_A;
    cout << "(A_T*A)^-1:" << endl << A_T_A_inv;
    cout << "A_T*b:" << endl << A_T_b_product;
    cout << "x~:" << endl << x;

    return 0;
}
