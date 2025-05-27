#include <iostream>
#include <fstream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <random>


using namespace std;


class RandomNumberGenerator 
{
    private:
        std::random_device rd;
        std::mt19937 gen;
        std::uniform_real_distribution<double> dist;
    
    public:
        RandomNumberGenerator()
            : gen(rd()), dist(0.0, 1.0) {}

        double get_random() { 
            return dist(gen);
        }
};


int main()
{
    RandomNumberGenerator rn;
    int nx = 30, ny = 30, VL = 1, VT = -1, VB = -1, eps = 1, pmax = 1;
    double Delta = 0.1;
    double xmax = Delta*nx, ymax = Delta*ny, sigma_p = xmax/10;
    double omega = 1.8, tol = 0.000001, F_old = 0, F_new = 0;

    double V[nx+1][ny+1] = {0};
    double sigma_V[nx+1][ny+1] = {0};
    double p[nx+1][ny+1];
    double B[nx+1][ny+1] = {0};
    double S[nx+1][ny+1] = {0};

    //Parametry do zmiany dla 3 przypadków z polecenia
    int N_chains = 100, n_length = 100;

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


    for (int i0 = 1; i0<nx; i0++)
        {
            for (int j0 = 1; j0<ny; j0++)
            {
                double sum_V1 = 0, sum_V2 = 0;
                int k_chains = 0;
                for (int N = 1; N<=N_chains; N++)
                {
                    int i = i0, j = j0;
                    double g = 0;
                    for (int n = 1; n<=n_length; n++)
                    {
                        double U = rn.get_random();
                        int m = floor(rn.get_random()*U);
                        if (m == 0) i--;
                        else if (m == 1) i++;
                        else if (m == 2) j--;
                        else if (m == 3) j++;
                        if (i==(nx+1)) i = nx - 1;
                        if (B[i][j] == 1)
                        {
                            double dV = V[i][j] + g;
                            sum_V1 += dV;
                            sum_V2 += dV*dV;
                            k_chains += 1;
                            break;
                        }
                        g = g + p[i][j]*Delta*Delta/(4*eps);
                    } //n
                }//N
                double V1 = sum_V1/k_chains;
                double V2 = sum_V2/k_chains;
                V[i0][j0] = V1;
                sigma_V[i0][j0] = sqrt((V2 - V1*V1)/(k_chains));
                B[i0][j0] = 0; // PRZYPADEK 1
                S[i0][j0] = k_chains/N_chains;
            }
        }

    return 0;
}