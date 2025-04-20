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


int main()
{
    std::srand(static_cast<unsigned>(std::time(0)));
    double R_a = 2;
    double R_b = sqrt(2)*R_a;
    double x_a = R_b + 0.5*R_a;
    int N = 1000000;

    ofstream plik1;
    plik1.open("C:\\Users\\kacpe\\OneDrive\\Pulpit\\C_plus\\Monte_Carlo\\Monte_carlo_4\\monte_carlo_4_5.txt");

    vector<double> X(N);
    vector<double> Y(N);
    double srednia_b = 0;
    int k =2;

   for (int i = 0; i<=N; i++)
   {
        double U1 = uniform();
        double U2 = uniform();
        double q = sqrt(uniform());

        X[i] = sqrt(-2*log(U1))*sin(2*M_PI*U2);
        Y[i] = sqrt(-2*log(U1))*cos(2*M_PI*U2);

        double R = sqrt(X[i]*X[i]+Y[i]*Y[i]);
        X[i] = X[i]/R;
        Y[i] = Y[i]/R;

        X[i] = q*X[i]*R_a + x_a;
        Y[i] = q*Y[i]*R_a;

        if (((X[i]-x_a)*(X[i]-x_a) + Y[i]*Y[i] <= R_a*R_a) && (X[i]*X[i] + Y[i]*Y[i] <= R_b*R_b))
        {
            // cout<<i<<endl;
            srednia_b += M_PI * R_a*R_a;
        }

        if( i == pow(10, k) )
        {
            cout<<k;
            double srednia = srednia_b/i;
            double wariancja = sqrt((M_PI*R_a*R_a*srednia - srednia*srednia)/i);
            plik1<<srednia<<" "<<wariancja<<" "<<i<<endl;
            k++;

        }



   }

   ofstream plik2;
   plik2.open("C:\\Users\\kacpe\\OneDrive\\Pulpit\\C_plus\\Monte_Carlo\\Monte_carlo_4\\monte_carlo_4_6.txt");

   X.clear();
   Y.clear();
   double srednia_d = 0, wariancja_d=0;
   k=2;

  for (int i = 0; i<=N; i++)
  {
       double U1 = uniform();
       double U2 = uniform();
       double q = sqrt(uniform());

       X[i] = sqrt(-2*log(U1))*sin(2*M_PI*U2);
       Y[i] = sqrt(-2*log(U1))*cos(2*M_PI*U2);

       double R = sqrt(X[i]*X[i]+Y[i]*Y[i]);
       X[i] = X[i]/R;
       Y[i] = Y[i]/R;

       X[i] = q*X[i]*R_b;
       Y[i] = q*Y[i]*R_b;

    //    plik2<<X[i]<<" "<<Y[i]<<endl;

    if (((X[i]-x_a)*(X[i]-x_a) + Y[i]*Y[i] <= R_a*R_a) && (X[i]*X[i] + Y[i]*Y[i] <= R_b*R_b))
    {
        // cout<<i<<endl;
         srednia_d += M_PI * R_b*R_b;
    }


       if( i == pow(10, k) )
       {
           cout<<k;
           double srednia = srednia_d/i;
           double wariancja = sqrt((M_PI*R_b*R_b*srednia - srednia*srednia)/i);
           plik2<<srednia<<" "<<wariancja<<" "<<i<<endl;
           k++;

       }

  }



}