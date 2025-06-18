#include <iostream>
#include <fstream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <random>
#include<cmath>
#include<vector>
#include<stdlib.h>
#include<stdio.h>

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

struct Result {
    double mean;
    double variance;
};

double calc_e_loc(double a, double c, double r)
{
    double e_loc = (-a*a*c*r*r + (-a*a + 4*a*c -2*c)*r + 2*a -2*c -2)/(2*c*r*r + 2*r);
    return e_loc;
};

double calc_p(double r, double a, double c)
{
    double exp_val = exp(-2 * a * r); 
    double psi_sq = (1 + c * r) * (1 + c * r); 
    return r * r * psi_sq * exp_val;
}



Result integral(double a, double c, double dr, int N)
{
    RandomNumberGenerator rn;

    double sum = 0;
    double sum_sq = 0;
    double ri = 1.0; 

    for (int i = 0; i < N; i++)
    {   
        double U1 = rn.get_random();
        double U2 = rn.get_random();
        double r_new = ri + dr * (2 * U1 - 1);

        if (r_new > 0)
        {
            double p_new = calc_p(r_new, a, c);
            double p_old = calc_p(ri, a, c);
            double p_acc = std::min(p_new / p_old, 1.0);

            if (U2 <= p_acc)
                ri = r_new;
        }

        double e = calc_e_loc(a, c, ri);
        sum += e;
        sum_sq += e * e;
    }

    double mean = sum / N;
    double variance = (sum_sq / N) - (mean * mean);

    return {mean, variance};
}




int main()
{
    int N = 10000;
    double dr = 0.1, da = 0.02, dc = 0.02;

    std::ofstream outfile("C:\\Users\\kacpe\\OneDrive\\Pulpit\\C_plus\\Monte_Carlo\\Monte_carlo_14\\energia.txt");
    std::ofstream outfile2("C:\\Users\\kacpe\\OneDrive\\Pulpit\\C_plus\\Monte_Carlo\\Monte_carlo_14\\wariancja.txt");

    for (double c = -0.7; c <= 0.3 + 1e-9; c += dc){
         for (double a = 0.3; a <= 1.2 + 1e-9; a += da) 
        {
            Result res = integral(a, c, dr, N);
            outfile << res.mean << " ";
            outfile2 << res.variance << " ";
        }
        outfile << endl;
        outfile2 << endl;
    }

    outfile.close();
    return 0;


    return 0;
}