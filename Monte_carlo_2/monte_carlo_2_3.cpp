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

double g_uniform()
{

    return 1.15*(static_cast<double>(std::rand()) / RAND_MAX);
}

double fun (double x)
{

    return 0.8*(1+x-x*x*x);
}


int main ()
{
    int N = 1000000;
    vector<double> tab(N);

    for (int  i = 0; i<N;)
    {
        double g2 = g_uniform();
        double u1 = uniform();
        if(g2 <= fun(u1))
        {
            tab[i] = u1;
            i++;
        } 

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