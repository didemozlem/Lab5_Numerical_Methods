#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace std;


double PodintegralnayaFunc(const double& x)
{

    return pow((x + x * x * x), 0.5);
}

double CubeFunction(const double& x, const double& y)
{
    return x * x + 2 * y;
}

double SimpsonMethod(const double& a, const double& b, const double& eps, double function(const double&))
{
    int N = 2;
    double oddAmount = 0;
    double evenAmount = 0;
    double h = (b - a) / N;
    double I = function(a) + function(b);
    double I2;

    do
    {
        I2 = I;
        N *= 2;
        h = (b - a) / N;
        double xi = a + h;
        oddAmount = 0;
        evenAmount = 0;
        for (int i = 1; i < N; i++) {
            if (i % 2 == 0)
            {
                evenAmount += function(xi);
            }
            else
            {
                oddAmount += function(xi);
            }
            xi += h;
        }

        evenAmount *= 2;
        oddAmount *= 4;
        I += (evenAmount + oddAmount);
        I *= (h / 3);
    } while (fabs(I - I2) > 15 * eps);

    return I;
}


double SimpsonCubeMethod(const double& a, const double& b, const double& c, const double& d,
                         double function(const double&, const double&))
{

    int m = 4;
    int n = 4;
    double hx = (b - a) / (2 * n);
    double hy = (d - c) / (2 * m);
    double sum = 0.0;
    double I = 0.0;

    vector<double> Xi;
    for (double xi = a; xi <= b; xi += hx)
    {
        Xi.push_back(xi);
    }

    vector<double> Yi;
    for (double yi = c; yi <= d; yi += hy)
    {
        Yi.push_back(yi);
    }


    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            sum += function(Xi[2 * i], Yi[2 * j]);
            sum += 4 * function(Xi[2 * i + 1], Yi[2 * j]);
            sum += function(Xi[2 * i + 2], Yi[2 * j]);
            sum += 4 * function(Xi[2 * i], Yi[2 * j + 1]);
            sum += 16 * function(Xi[2 * i + 1], Yi[2 * j + 1]);
            sum += 4 * function(Xi[2 * i + 2], Yi[2 * j + 1]);
            sum += function(Xi[2 * i], Yi[2 * j + 2]);
            sum += 4 * function(Xi[2 * i + 1], Yi[2 * j + 2]);
            sum += function(Xi[2 * i + 2], Yi[2 * j + 2]);
        }
    }
    I += sum;
    I *= ((hx * hy) / 9);
    return I;
}


double TrapezoidalMethod(const double& a, const double& b, const double& eps, double function(const double&))
{
    int n = 1;
    double h = (b - a);
    double xi = a;
    double sum = 0;
    double I = function(a) + function(b);
    I *= h;
    double I2;
    do
    {
        n *= 2;
        I2 = I;
        h = (b - a) / n;
        xi = a;
        sum = 0;
        for (int i = 0; i < n - 1; i++)
        {
            xi += h;
            sum += function(xi);
        }
        I = (function(a) + function(b) + 2 * sum) * (h / 2);
    } while (fabs(I - I2) >= 3 * eps);

    return I;
}



int main()
{
    const double eps1 = 1e-4;
    const double eps2 = 1e-5;
    double a = 0.6;
    double b = 1.724;
    cout << "Trapezoidal Method (eps = 10^-4): " << endl;
    cout << setprecision(20) << "I: " << TrapezoidalMethod(a, b, eps1, PodintegralnayaFunc) << endl;
    cout << "Trapezoidal Method (eps = 10^-5): " << endl;
    cout << setprecision(20) << "I: " << TrapezoidalMethod(a, b, eps2, PodintegralnayaFunc) << endl;
    cout << "Simpson Method (eps = 10^-4): " << endl;
    cout << setprecision(20) << "I: " << SimpsonMethod(a, b, eps1, PodintegralnayaFunc) << endl;
    cout << "Simpson Method (eps = 10^-5): " << endl;
    cout << setprecision(20) << "I: " << SimpsonMethod(a, b, eps2, PodintegralnayaFunc) << endl;

    //Variant 31 for simpson cube method
    double a_ = 0;
    double b_ = 2.0;
    double c_ = 0;
    double d_ = 1.0;
    cout << "Simson Cube Method = " << SimpsonCubeMethod(a_, b_, c_, d_, CubeFunction);
    cout << endl;


    return 0;
}
