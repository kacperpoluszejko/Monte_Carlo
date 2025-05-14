#include <iostream>
#include <fstream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <random>


using namespace std;

// Do losowania
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

struct atom
{
    double x;
    double y;
    double z;
    double r;
    double theta;
    double fi;
};

const double c0 = 19, d0 = 2.5, a0 = 0.011304, delta = 0.80460, S = 1.29;
const double R0 = 1.315, R1 = 1.7, R2 = 2.0;
const int  N=60;


void convertToSpherical(atom &a)
{
    a.r = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);

    if (a.r != 0) 
        a.theta = acos(a.z / a.r);
    else 
        a.theta = 0;

    a.fi = atan2(a.y, a.x);
}

double calc_g_theta_ijk(atom i, atom j, atom k)
{
    atom r_ij = {j.x - i.x, j.y - i.y, j.z - i.z, 0, 0, 0};
    atom r_ik = {k.x - i.x, k.y - i.y, k.z - i.z, 0, 0, 0};
    convertToSpherical(r_ij);
    convertToSpherical(r_ik);

    double cos_theta_ijk = (r_ij.x*r_ik.x + r_ij.y*r_ik.y + r_ij.z*r_ik.z)/(r_ij.r*r_ik.r);

    double g_theta_ijk = a0*(1 + (c0*c0)/(d0*d0) - (c0*c0)/(d0*d0 + pow(1+cos_theta_ijk,2)) );

    return g_theta_ijk;
}

double calc_f_cut(double r)
{
    if(r<=R1) return 1;
    else if(r>R1 && r<=R2) return 0.5*(1+cos(M_PI*((r-R1)/(R2-R1))));
    else return 0;

}

double calc_B_ij(atom atoms[N], atom i, atom j)
{
    double zeta_ij = 0;
    for (int it = 0; it<N; it++)
    {   atom k = atoms[it];
        atom r_ik = {k.x - i.x, k.y - i.y, k.z - i.z, 0, 0, 0};
        convertToSpherical(r_ik);
        zeta_ij += calc_f_cut(r_ik.r)*calc_g_theta_ijk(i, j, k);
    }

    double B_ij = pow(1+zeta_ij, -delta);
    return B_ij;
}

int main()
{
    RandomNumberGenerator rn;
    atom atoms[N];


    ifstream plik1;
    plik1.open("C:\\Users\\kacpe\\OneDrive\\Pulpit\\C_plus\\Monte_Carlo\\Monte_carlo_8\\c60.txt");
    if (!plik1.is_open())
    {
        cerr << "Nie udało się otworzyć pliku c60.txt!" << endl;
        return 1;
    }

    for (int i = 0; i<60; i++)
    {
        plik1>>atoms[i].x>>atoms[i].y>>atoms[i].z;
        convertToSpherical(atoms[i]);
        
    }

    plik1.close();

    for (int i = 0; i < 60; i++)
    {
        cout << "Atom " << i << ": x=" << atoms[i].x << ", y=" << atoms[i].y << ", z=" << atoms[i].z << endl;
    }

    return 0;
}