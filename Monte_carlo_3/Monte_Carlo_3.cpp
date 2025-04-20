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
    int n = 10000;

    std::srand(static_cast<unsigned>(std::time(0)));

    vector<double> X(n);
    vector<double> Y(n);



    ofstream plik1;
    plik1.open("C:\\Users\\kacpe\\OneDrive\\Pulpit\\C_plus\\Monte_Carlo\\Monte_carlo_3\\monte_carlo_gauss.txt");

   for (int i = 0; i<n; i++)
   {
        double U1 = uniform();
        double U2 = uniform();

        X[i] = sqrt(-2*log(1-U1))*cos(2*M_PI*U2);
        Y[i] = sqrt(-2*log(1-U1))*sin(2*M_PI*U2);

        plik1<<X[i]<<" "<<Y[i]<<endl;

   }

   vector<double> X2(n);
   vector<double> Y2(n);

   ofstream plik2;
   plik2.open("C:\\Users\\kacpe\\OneDrive\\Pulpit\\C_plus\\Monte_Carlo\\Monte_carlo_3\\monte_carlo_gauss_2.txt");

   for (int i = 0; i<n; i++)
   {
        X2[i] = X[i]/(sqrt(X[i]*X[i]+Y[i]*Y[i]));
        Y2[i] = Y[i]/(sqrt(X[i]*X[i]+Y[i]*Y[i]));

        plik2<<X2[i]<<" "<<Y2[i]<<endl;

   }

   vector<double> X3(n);
   vector<double> Y3(n);

   ofstream plik3;
   plik3.open("C:\\Users\\kacpe\\OneDrive\\Pulpit\\C_plus\\Monte_Carlo\\Monte_carlo_3\\monte_carlo_gauss_3.txt");

   for (int i = 0; i<n; i++)
   {
        double U1 = uniform();
        double R = sqrt(U1);

        X3[i] = X2[i]*R;
        Y3[i] = Y2[i]*R;

        plik3<<X3[i]<<" "<<Y3[i]<<endl;

   }

   std::vector<std::vector<int>> vec(n, std::vector<int>(2, 0));

   double b1 = 1, b2 = 0.2, alpha = M_PI/4;
   double r1[2];
   double r2[2];

   vector<double> X4(n);
   vector<double> Y4(n);

   r1[0] = b1*cos(alpha);
   r1[1] = b1*sin(alpha);
   r2[0] = (-1)*b2*sin(alpha);
   r2[1] = b2*cos(alpha);

   ofstream plik4;
   plik4.open("C:\\Users\\kacpe\\OneDrive\\Pulpit\\C_plus\\Monte_Carlo\\Monte_carlo_3\\monte_carlo_gauss_4.txt");



   for (int i = 0; i<n; i++)
   {
        double U1 = uniform();
        double R = sqrt(U1);

        X4[i] = r1[0]*X[i] + r2[0]*Y[i];
        Y4[i] = r1[1]*X[i] + r2[1]*Y[i];

        plik4<<X4[i]<<" "<<Y4[i]<<endl;

   }

   double x_bar=0, y_bar=0, x_bar2=0, y_bar2=0, xy_bar=0, sigma_x=0, sigma_y=0, sigma_xy=0;


   for (int i = 0; i<n; i++)
   {
       x_bar += X4[i];
       y_bar += Y4[i];
       x_bar2 += X4[i]*X4[i];
       y_bar2 += Y4[i]*Y4[i];
       xy_bar += Y4[i]*X4[i]; 
   }

   x_bar = x_bar/n;
   y_bar = y_bar/n;
   x_bar2 = x_bar2/n;
   y_bar2 = y_bar2/n;
   xy_bar = xy_bar/n;

//    cout<<x_bar<<"  "<<y_bar<<"  "<<x_bar2<<"  "<<y_bar2<<"  "<<xy_bar<<endl;

   sigma_x = x_bar2 - x_bar*x_bar;
   sigma_y = y_bar2 - y_bar*y_bar;
   sigma_xy = xy_bar -x_bar*y_bar;

   cout<<"Elementy macierzy kowariancji: "<<sigma_x<<"  "<< sigma_y<<"  "<<sigma_xy<<endl;

   cout<<"Elementy macierzy A: "<<r1[0]<<" "<<r1[1]<<" "<<r2[0]<<" "<<r2[1]<<endl;

   double r_xy = sigma_xy/(sqrt(sigma_x*sigma_y));

   cout<<"Współczynnik korelacji = "<<r_xy;


    return 0;

}