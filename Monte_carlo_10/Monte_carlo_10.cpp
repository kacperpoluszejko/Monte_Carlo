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

const double L = 0.25 * pow(10, -6);
const double C = 100 * pow(10, -12);
const double R = 12.5;
const double G = 0.5 * pow(10, -3);
const double l = 2;
const double R_l = 12.5;
const double R_g = 75;
const double v = 1 * pow(10, 9);
const double t0 = 7.5 * pow(10, -9);
const double sigma = 0.75 * pow(10, -9);
const double c = 1/(sqrt(L*C));
const double mu = G/C;
const double R_0 = sqrt(L/C);
const double lambda = 0.5*(R/L - G/C);
const double tau = (R_0)/(R_0 + R_g);
const double tau_g = (R_g - R_0)/(R_0 + R_g);
const double tau_l = (R_l - R_0)/(R_0 + R_l);

double V_g(double t)
{
    return sin(M_PI*2*v*t)*exp(-(t - t0)*(t - t0)/(2*sigma*sigma));
}

int main()
{
    RandomNumberGenerator rn;


    int npaths = 100000;
    double t_start = 35 * pow(10, -9);
    int Nx = 1000;
    double x_table [Nx] = {0};
    double dx = l/Nx;
    double F_table [Nx] = {0};
    double B_table [Nx] = {0};
    double fxt = 0, bxt = 0;

    for (int nx = 0; nx<Nx; nx++)
    {
        double x_start = nx*dx;
        x_table[nx] = x_start;
        for (int i = 1; i<=2; i++)
        {
            double suma = 0;
            for (int n = 0; n<npaths; n++)
            {
                double x = x_start;
                double t = t_start;
                double eta = 1;
                int sign = pow(-1, i);

                while (t>0)
                {
                    double s = -1*(log(rn.get_random()))/(lambda + mu);
                    if (sign == (-1))
                    {
                        if ((x-c*s)>0)
                        {
                            eta = eta*(lambda)/(lambda+mu);
                        }
                        else 
                        {
                            s = x/c;
                            suma += eta*tau*V_g(t-s);
                            eta = eta * tau_g;
                        }
                        x = x - c*s;
                        t = t - s;
                    }
                    else if (sign == 1)
                    {
                        if ((x+c*s)<l)
                        {
                            eta = eta*(lambda)/(lambda + mu);
                        }
                        else 
                        {
                            s = (l-x)/c;
                            eta = eta*tau_l;
                        }
                        x = x + c*s;
                        t = t-s;
                    }
                    sign = -sign;
                }
            }

            if (i == 1)
            {
                fxt = (suma/npaths);
            } 
            else if (i == 2)
            {
                 bxt = suma/npaths;
            }
        }
        F_table[nx] = fxt;
        B_table[nx] = bxt;

    }

    double u_table[Nx] = {0};
    double i_table[Nx] = {0};

    ofstream plik1;
    plik1.open("C:\\Users\\kacpe\\OneDrive\\Pulpit\\C_plus\\Monte_Carlo\\Monte_carlo_10\\monte_carlo_10_u4_3.txt");
    // ofstream plik2;
    // plik2.open("C:\\Users\\kacpe\\OneDrive\\Pulpit\\C_plus\\Monte_Carlo\\Monte_carlo_10\\monte_carlo_10_i5_3.txt");
    for (int nx = 0; nx<Nx; nx++)
    {
        u_table[nx] = F_table[nx]+B_table[nx];
        // i_table[nx] = (F_table[nx] - B_table[nx])/R_0;
        plik1<<u_table[nx]<<" "<<x_table[nx]<<endl;
        // plik2<<i_table[nx]<<" "<<x_table[nx]<<endl;
    }

    


    return 0;
}