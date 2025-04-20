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

double oblicz_srednia (double * tab, int rozmiar)
{
    double suma = 0;
    for (int i=0; i<rozmiar; i++)
    {
        suma += tab[i];
    }
    return suma/rozmiar;
}

double oblicz_srednia_kwadrat (double * tab1, double * tab2, int rozmiar)
{
    double suma = 0;
    for (int i=0; i<rozmiar; i++)
    {
        suma += tab1[i]*tab2[i];
    }
    return suma/rozmiar;
}

int main(){
    RandomNumberGenerator rn;

    int D = 1, N_max = 100000;

    double X[N_max] = {0}; // Tablica położeń x
    double Y[N_max] = {0}; // Tablica położeń y
    double dt = 0.1, t_max = 100;
    int Nt = t_max/dt; // liczba kroków czasowych
    double sigma = sqrt(2*D*dt); // odchylenie rozkładu nromalnego

    double Dxx[Nt] = {0};
    double Dyy[Nt] = {0};
    double Dxy[Nt] = {0};

    for (int t=0; t<Nt; t++)
    {
        for (int i=0; i<N_max; i++)
        {
            double U1 = rn.get_random();
            double U2 = rn.get_random();
            double dx = sqrt(-2*log(1-U1))*cos(2*M_PI*U2)*sigma;
            double dy = sqrt(-2*log(1-U1))*sin(2*M_PI*U2)*sigma;
            X[i] = X[i] + dx;
            Y[i] = Y[i] + dy;
        }
        double av_xx = oblicz_srednia(X, N_max);
        double av_yy = oblicz_srednia(Y, N_max);
        double av_xx_2 = oblicz_srednia_kwadrat(X,X, N_max);
        double av_yy_2 = oblicz_srednia_kwadrat(Y,Y, N_max);
        double av_xy = oblicz_srednia_kwadrat(X,Y, N_max);
        Dxx[t] = (av_xx_2 - av_xx*av_xx)/(2*t*dt);
        Dyy[t] = (av_yy_2 - av_yy*av_yy)/(2*t*dt);
        Dxy[t] = (av_xy - av_xx*av_yy)/(2*t*dt);

    }
    ofstream plik1;
    plik1.open("C:\\Users\\kacpe\\OneDrive\\Pulpit\\C_plus\\Monte_Carlo\\Monte_carlo_6\\monte_carlo_6_10^5.txt");

    for (int i = 0; i <Nt; i++)
    {
        plik1<<i*dt<<" "<<Dxx[i]<<" "<<Dyy[i]<<" "<<Dxy[i]<<endl;
    }




    return 0;
}