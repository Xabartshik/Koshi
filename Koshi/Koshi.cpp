#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double f(double x, double y)
{
    return sin(x) * cos(y);
    //    return (cos(6*x) * pow(sin(4*y), 2));
}


void result_print(vector<vector<double>> result)
{
    if (result[0].size() == 3)
    {
    for (int i = 0; i < result.size(); i++) {
        cout.width(2);
        cout << result[i][0] << " ";
        cout.width(10);
        cout << result[i][1] << " ";
        cout.width(10);
        cout << result[i][2] << endl;
    }
    }
    else if (result[0].size() == 5)
    {
        for (int i = 0; i < result.size(); i++) {
            cout.width(2);
            cout << result[i][0] << " ";
            cout.width(10);
            cout << result[i][1] << " ";
            cout.width(10);
            cout << result[i][2] << " ";
            cout.width(10);
            cout << result[i][3] << " ";
            cout.width(10);
            cout << result[i][4] << endl;
        }
    }
}



vector<vector<double>> Euler(double xo, double yo, double h, int k) {
    //Массив под результаты
    vector<vector<double>> result(k, vector<double>(3));
    vector<double> xx(k);
    vector<double> yy(k);
    xx[0] = xo;
    yy[0] = yo;
    for (int n = 0; n < k; n++) {
        result[n][0] = n;
        result[n][1] = xx[n];
        result[n][2] = yy[n];
        if (n < k - 1) {
            xx[n + 1] = xx[n] + h;
            yy[n + 1] = yy[n] + h * f(xx[n], yy[n]);
        }
    }
    return result;
}

vector<vector<double>> EulerMod(double xo, double yo, double h, int k) {
    //Массив под результаты
    vector<vector<double>> result(k, vector<double>(3));
    vector<double> xx(k);
    vector<double> yy(k);
    xx[0] = xo;
    yy[0] = yo;
    for (int n = 0; n < k; n++) {
        result[n][0] = n;
        result[n][1] = xx[n];
        result[n][2] = yy[n];
        if (n < k - 1) {
            xx[n + 1] = xx[n] + h;
            yy[n + 1] = yy[n] + h * f(xx[n] + h / 2, yy[n] + h / 2 * f(xx[n], yy[n]));
        }
    }
    return result;
}

vector<vector<double>> RungeKutt(double xo, double yo, double h, int k) {
    //Массив под результаты
    vector<vector<double>> result(k, vector<double>(3));
    vector<double> xx(k);
    vector<double> yy(k);
    double a, b, c, d;
    xx[0] = xo;
    yy[0] = yo;
    for (int n = 0; n < k; n++) {
        result[n][0] = n;
        result[n][1] = xx[n];
        result[n][2] = yy[n];
        if (n < k - 1) {
            xx[n + 1] = xx[n] + h;
             a = h * f(xx[n], yy[n]);
             b = h * f(xx[n] + h / 2, yy[n] + a / 2);
             c = h * f(xx[n] + h / 2, yy[n] + b / 2);
             d = h * f(xx[n] + h, yy[n] + c);
            yy[n + 1] = yy[n] + (a + 2 * b + 2 * c + d) / 6;
        }
    }
    return result;
}


vector<vector<double>> Adams(double xo, double yo, double h, int k) {
    //Массив под результаты
    vector<vector<double>> result(k, vector<double>(5));
    vector<double> xx(k + 1);
    vector<double> yy(k + 1);
    xx[0] = xo;
    yy[0] = yo;
    double a, b, c, d, aa, bb, cc, dd;
    // m = 1
    for (int n = 0; n < k; n++) {
        result[n][0] = n;
        result[n][1] = xx[n];
        result[n][2] = yy[n];
        if (n < 1) {
            xx[n + 1] = xx[n] + h;
             a = h * f(xx[n], yy[n]);
             b = h * f(xx[n] + h / 2, yy[n] + a / 2);
             c = h * f(xx[n] + h / 2, yy[n] + b / 2);
             d = h * f(xx[n] + h, yy[n] + c);
            yy[n + 1] = yy[n] + (a + 2 * b + 2 * c + d) / 6;
        }
        else {
            if (n < k) {
                xx[n + 1] = xx[n] + h;
                 aa = f(xx[n], yy[n]);
                 bb = f(xx[n - 1], yy[n - 1]);
                yy[n + 1] = yy[n] + (h / 2) * (3 * aa - bb);
            }
        }
    }

    // m = 2
    for (int n = 0; n < k; n++) {
        result[n][0] = n;
        result[n][1] = xx[n];
        result[n][3] = yy[n];
        if (n < 2) {
            xx[n + 1] = xx[n] + h;
             a = h * f(xx[n], yy[n]);
             b = h * f(xx[n] + h / 2, yy[n] + a / 2);
             c = h * f(xx[n] + h / 2, yy[n] + b / 2);
             d = h * f(xx[n] + h, yy[n] + c);
            yy[n + 1] = yy[n] + (a + 2 * b + 2 * c + d) / 6;
        }
        else {
            if (n < k) {
                xx[n + 1] = xx[n] + h;
                 aa = f(xx[n], yy[n]);
                 bb = f(xx[n - 1], yy[n - 1]);
                 cc = f(xx[n - 2], yy[n - 2]);
                yy[n + 1] = yy[n] + (h / 12) * (23 * aa - 16 * bb + 5 * cc);
            }
        }
    }

    // m = 3
    for (int n = 0; n < k; n++) {
        result[n][0] = n;
        result[n][1] = xx[n];
        result[n][4] = yy[n];
        if (n < 3) {
            xx[n + 1] = xx[n] + h;
             a = h * f(xx[n], yy[n]);
             b = h * f(xx[n] + h / 2, yy[n] + a / 2);
             c = h * f(xx[n] + h / 2, yy[n] + b / 2);
             d = h * f(xx[n] + h, yy[n] + c);
            yy[n + 1] = yy[n] + (a + 2 * b + 2 * c + d) / 6;
        }
        else {
            if (n < k) {
                xx[n + 1] = xx[n] + h;
                 aa = f(xx[n], yy[n]);
                 bb = f(xx[n - 1], yy[n - 1]);
                 cc = f(xx[n - 2], yy[n - 2]);
                 dd = f(xx[n - 3], yy[n - 3]);
                yy[n + 1] = yy[n] +
                    (h / 24) * (55 * aa - 59 * bb + 37 * cc - 9 * dd);
            }
        }
    }
    return result;
}

vector<vector<double>> Miln(double xo, double yo, double h, double eps, int k) {
    //Массив под результаты
    vector<vector<double>> result(k, vector<double>(3));
    vector<double> x(k);
    vector<double> y(k);
    x[0] = xo;
    y[0] = yo;
    double a, b, c, d, aa, bb, cc, dd, y1, y2;
    result[0][0] = 0;
    result[0][1] = x[0];
    result[0][2] = y[0];

    bool isDone = false;
    while (!isDone) {
        for (int n = 1; n < k; n++) {
            if (n <= 3) {
                x[n] = x[n - 1] + h;
                 a = h * f(x[n - 1], y[n - 1]);
                 b = h * f(x[n - 1] + h / 2, y[n - 1] + a / 2);
                 c = h * f(x[n - 1] + h / 2, y[n - 1] + b / 2);
                 d = h * f(x[n - 1] + h, y[n - 1] + c);
                y[n] = y[n - 1] + (a + 2 * b + 2 * c + d) / 6;
            }
            else {
                if (n <= k) {
                    x[n] = x[n - 1] + h;
                     aa = f(x[n - 1], y[n - 1]);
                     bb = f(x[n - 2], y[n - 2]);
                     cc = f(x[n - 3], y[n - 3]);
                     dd = f(x[n - 4], y[n - 4]);
                     y1 = y[n - 4] + (4 * h / 3) * (2 * cc - bb + 2 * aa);
                     y2 = y[n - 2] + (h / 3) * (bb + 4 * aa + y1);
                    if (abs(y2 - y1) / 29 < eps) {
                        y[n] = y2;
                    }
                    else {
                        h = h / 2;
                        cout << "Пришлось уменьшить шаг..." << endl;
                        break;
                    }
                }
            }
            if (n == (k - 1)) {
                isDone = true;
            }
            result[n][0] = n;
            result[n][1] = x[n];
            result[n][2] = y[n];
        }
    }
    return result;
}


int main() {
    double xo = 0, yo = 1, h = 0.75, eps = 0.001;
    //    double xo = 0, yo = 1, h = 0.1, eps = 0.001;  //Для Милна
    //    double xo = 0, yo = 0.250, h = 0.1, eps = 0.001;
    int k = 11;
    //cout << "Enter initial conditions (xo, yo): ";
    //cin >> xo >> yo;
    //cout << "Enter step size (h): ";
    //cin >> h;
    //cout << "Enter number of steps (k): ";
    //cin >> k;

    vector<vector<double>> result_euler(k, vector<double>(3));
    vector<vector<double>> result_em(k, vector<double>(3));
    vector<vector<double>> result_rk(k, vector<double>(3));
    vector<vector<double>> result_adams(k, vector<double>(5));
    vector<vector<double>> result_miln(k, vector<double>(3));
    result_euler = Euler(xo, yo, h, k);
    result_em = EulerMod(xo, yo, h, k);
    result_rk = RungeKutt(xo, yo, h, k);
    result_adams = Adams(xo, yo, h, k);
    result_miln = Miln(xo, yo, h, eps, k);
    cout << "Results Euler: " << endl;
    result_print(result_euler);
    cout << "Results Euler Mod: " << endl;
    result_print(result_em);
    cout << "Results Runge-Kutt: " << endl;
    result_print(result_rk);
    cout << "Results Adams: " << endl;
    result_print(result_adams);
    cout << "Results Miln: " << endl;
    result_print(result_miln);
    return 0;
}