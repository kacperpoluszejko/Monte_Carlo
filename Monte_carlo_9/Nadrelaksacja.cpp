#include <iostream>
#include <fstream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <random>


using namespace std;

int main()
{
    int nx = 30, ny = 30, VL = 1, VT = -1, VB = -1, eps = 1, pmax = 1;
    double Delta = 0.1;
    double xmax = Delta*nx, ymax = Delta*ny, sigma_p = xmax/10;
    double omega = 1.8, tol = 0.000001, F_old = 0, F_new = 0;

    double V[nx+1][ny+1] = {0};
    double p[nx+1][ny+1];

    int itmax = 10000;

    for (int i = 0; i<=nx; i++)
    {
        for (int j = 0; j<=ny; j++)
        {
            p[i][j] = pmax*exp(-( (i*Delta - xmax/2)*(i*Delta - xmax/2) + (j*Delta - ymax/2)*(j*Delta - ymax/2) )/(2*sigma_p*sigma_p));
        }
    }

    for (int j = 0; j<=ny; j++)
    {
        V[0][j] = VL * sin(M_PI*Delta*j/ymax);
        for (int i = 0; i<=nx; i++)
        {
            V[i][0] = VB*sin(M_PI*Delta*i/xmax);
            V[i][ny] = VT*sin(M_PI*Delta*i/xmax);
        }
    }

    for (int it = 0; it<itmax; it++)
    {
        for (int i = 1; i<nx; i++)
        {
            for (int j = 1; j<ny; j++)
            V[i][j] = (1-omega)*V[i][j] + (omega/4)*(V[i+1][j] + V[i-1][j] + V[i][j+1] + V[i][j-1] + Delta*Delta*p[i][j]/eps);
        }
        for(int j=1; j<ny;j++) V[nx][j] = V[nx-1][j];

        F_old = F_new;
        F_new = 0;
        for (int i = 1; i<nx; i++)
        {
            for (int j = 1; j<ny; j++)
            {
                double Ex = (V[i+1][j] - V[i-1][j])/(2*Delta);
                double Ey = (V[i][j+1] - V[i][j-1])/(2*Delta);
                F_new = F_new + (Ex*Ex + Ey*Ey)/2 - p[i][j]*V[i][j];
            }
        }
        if(abs((F_new - F_old)/(F_new))< tol) break;
        cout<<abs((F_new - F_old)/(F_new))<<endl;
    }
    

    return 0;
}