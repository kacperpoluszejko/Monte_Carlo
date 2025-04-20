#include <iostream>
#include <fstream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <conio.h>
#include <math.h>
#include <vector>
#include <ctime>
#include <algorithm>

using namespace std;

double uniform ()
{
    return static_cast<double>(std::rand()) / RAND_MAX;
}

double fun (double x)
{

    return 0.8*(1+x-x*x*x);
}

double p_acc(double f1, double f2)
{
    double a = 1;
    return min(f2/f1, a);
}

int main ()
{

    double delta = 0.05;
    int N = 1000000;

    vector<double> tab(N);

    // double f0 = fun(uniform());
    double x0 = 0.05;
    tab[0] = x0;

    for (int i=1; i<N; i++)
    {
        double u1 = uniform();
        double u2 = uniform();
        double x_new = tab[i-1] + delta*(2*u1-1);
        double f_x1 = fun(tab[i-1]);
        double f_x2 = fun(x_new);
        double p = p_acc(f_x1, f_x2);
        if(x_new >=0 && x_new <= 1 && u2<=p)
        {
            tab[i] = x_new;
        } 
        else tab[i] = tab[i-1];
    }


    ofstream plik1;
    plik1.open("C:\\Users\\kacpe\\OneDrive\\Pulpit\\C_plus\\FUZ\\monte_carlo_2_1.txt");
    
    for (int i = 0; i<N; i++)
    {
        plik1<<tab[i]<<endl;

    }

    cout<<"koniec";
    return 0;






}