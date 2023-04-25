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

    int min_x = *min_element(t.begin(), t.end()) - 5;
    int max_x = *max_element(t.begin(), t.end()) + 5;
    int min_y = *min_element(b.begin(), b.end()) - 5;
    int max_y = *max_element(b.begin(), b.end()) + 5;

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

    // GNUPlot
    FILE *pipe = popen("gnuplot -persist", "w");

    if (pipe == NULL) {
        cout << "Pipe failed" << endl;
        return 0;
    }

    // Setup GNUPlot
    fprintf(pipe, "set xrange [%d:%d]\n", min_x, max_x);
    fprintf(pipe, "set yrange [%d:%d]\n", min_y, max_y);
    fprintf(pipe, "set xtics 1\n");
    fprintf(pipe, "set ytics 1\n");
    fprintf(pipe, "set grid xtics ytics\n");

    // Function print
    string fx = "";
    for (int i = 0; i < x.n; i++) {
        fx += to_string(x(i, 0)) + " * x**" + to_string(i);
        if (i != x.n - 1)
            fx += " + ";
    }

    cout << fx << endl;

    fprintf(pipe, "plot %s with lines,", fx.c_str());

    // Print points
    fprintf(pipe, " '-' with points pointtype 7 pointsize 1\n");
    for (int i = 0; i < t.size(); i++)
        fprintf(pipe, "%f %f\n", t[i], b[i]);

    fflush(pipe);
    pclose(pipe);


    return 0;
}
