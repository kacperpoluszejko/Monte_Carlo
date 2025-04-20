#include <iostream>
#include <fstream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <conio.h>
#include <math.h>
#include <vector>
#include <ctime>

using namespace std;

double uniform ()
{
    return static_cast<double>(std::rand()) / RAND_MAX;
}


int main ()
{
    int N = 1000000;
    double g1 = 0.8, g2 = 0.2;
    std::vector<double> tab(N, 0);
    std::srand(static_cast<unsigned>(std::time(0)));
    // double a = uniform();

    for (int i = 0; i<N; i++)
    {
        double a = uniform();
        double b = uniform();
        // cout<<a<<"  "<<b<<endl;
        if (a<=g1) tab[i] = b;
        if (a>g1) tab[i] = sqrt(1-sqrt(1-b));
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
