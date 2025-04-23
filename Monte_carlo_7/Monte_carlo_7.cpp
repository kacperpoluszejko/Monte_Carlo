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

    double k1 = 1, k2 = 1, k3 = 0.001, k4 = 0.01;
    double x1_0 = 120, x2_0 = 80, x3_0 = 1;
    int t_max = 200, N = 50, P_max = 5, m = 0;
    double dt = t_max/N;
    double h0[N] = {0}, h1[N] = {0}, h2[N] = {0};
    int ncount[N] = {0};

    for (int p=1; p<=P_max; p++)
    {
        double t = 0;
        double x1 = x1_0, x2 = x2_0, x3 = x3_0;
        h0[N] = {0};
        ncount[N] = {0};

        while (t < t_max)
        {
            double g1 = k1, g2 = k2, g3 = k3*x1*x2, g4 = k4*x3;
            double g_max = g1 + g2 + g3 + g4;
            double U1 = rn.get_random(), U2 = rn.get_random();
            double delta_t = (-1)*log(U1)/g_max;

            if(U2 <= g1/g_max)
            {
                m=1;
                x1 +=1;
            }
            else if (U2 <= (g1 + g2)/g_max) 
            {
                m=2;
                x2 +=1;
            }
            else if (U2 <= (g1 + g2 + g3)/g_max)
            {
                m=3;
                x1 = x1 - 1;
                x2 = x2 - 1;
                x3 = x3 + 1; 
            }
            else
            {
                m=4;
                x3 = x3 - 1;
            }
            t = t + delta_t;

            // plik1<<t<<" "<<x1<<" "<<x2<<" "<<x3<<endl;

            //Wkład pojedynczej ścieżki
            int l = floor(t/dt);
            h0[l] += x3;
            ncount[l]++;

        }

        //HISTOGRAM Uśrednianie po wielu ścieżkach
        for(int l=0; l<N; l++)
        {
            double X3 = h0[l]/ncount[l];
            h1[l] += X3;
            h2[l] += X3*X3;
        }
    }
    ofstream plik1;
    plik1.open("C:\\Users\\kacpe\\OneDrive\\Pulpit\\C_plus\\Monte_Carlo\\Monte_carlo_7\\monte_carlo_7_zad2_2.txt");

    //HISTOGRAM - przetwarzanie z wszystkich ścieżek
    for(int l=0; l<N; l++)
    {
        double x3_av= h1[l]/P_max;
        double x3_av2 = h2[l]/P_max;
        double sigma = sqrt((x3_av2 - x3_av*x3_av)/P_max);
        double t = (l + 0.5)*dt;
        plik1<<t<<" "<<x3_av<<" "<<sigma<<endl;
    }



    return 0;
}