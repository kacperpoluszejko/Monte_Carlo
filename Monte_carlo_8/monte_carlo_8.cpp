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

const double angstrom_to_meter = 1e-10;
const double eV_to_joule = 1.602e-19;

const double c0 = 19, d0 = 2.5, a0 = 0.011304, delta = 0.80460, S = 1.29;
const double R0 = 1.315 * angstrom_to_meter;
const double R1 = 1.70 * angstrom_to_meter;
const double R2 = 2.0 * angstrom_to_meter;
const double De = 6.325 * eV_to_joule;
const double lambda = 1.5 / angstrom_to_meter;
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

double calc_bar_B_ij(atom atoms[N], atom i, atom j)
{
    double zeta_ij = 0;
    for (int it = 0; it<N; it++)
    {   atom k = atoms[it];
        if (fabs(i.r - k.r) > 1e-10 && fabs(j.r - k.r) > 1e-10)
        {
            atom r_ik = {k.x - i.x, k.y - i.y, k.z - i.z, 0, 0, 0};
            convertToSpherical(r_ik);
            zeta_ij += calc_f_cut(r_ik.r)*calc_g_theta_ijk(i, j, k);
        }
    }

    double zeta_ji = 0;
    for (int it = 0; it<N; it++)
    {   atom k = atoms[it];
        if (fabs(i.r - k.r) > 1e-10 && fabs(j.r - k.r) > 1e-10)
        {
            atom r_jk = {k.x - j.x, k.y - j.y, k.z - j.z, 0, 0, 0};
            convertToSpherical(r_jk);
            zeta_ji += calc_f_cut(r_jk.r)*calc_g_theta_ijk(j, i, k);
        }
    }

    double B_ij = pow(1+zeta_ij, -delta);
    double B_ji = pow(1+zeta_ji, -delta);

    return (B_ij + B_ji)/2;
}

double calc_V_R(double r)
{
    double VR = De/(S-1)*exp(-sqrt(2*S)*(r-R0));
    return VR;
}

double calc_V_A(double r)
{
    double VA = De*S/(S-1)*exp(-sqrt(2/S)*(r-R0));
    return VA;
}

double calc_V_i(atom atoms[N], atom i, int atom_index)
{
    double Vi = 0;
    for (int it = 0; it<N; it++)
    {
        if (atom_index != it)
        {
            atom j = atoms[it];
            double rij = sqrt(pow(j.x - i.x, 2) + pow(j.y - i.y, 2) + pow(j.z - i.z, 2));
            cout<<rij<<" "<<Vi<<" "<<calc_f_cut(rij)<<" "<<calc_V_R(rij)<<" "<<calc_bar_B_ij(atoms, i, j)<<" "<<calc_V_A(rij)<<endl;
            Vi += calc_f_cut(rij)*(calc_V_R(rij) - calc_bar_B_ij(atoms, i, j)*calc_V_A(rij));
        }
    }

    return Vi;
}

double calc_total_V(atom atoms[N])
{
    double V_total = 0;
     for (int it = 0; it<N; it++)
    {
        atom i = atoms[it];
        V_total += calc_V_i(atoms, i, it);
    }
    return V_total/2;
}

void f_corelation(atom atoms[N], double *pcf, int M)
{
    double r_sr = 0;
    for (int i = 0; i < N; i++)
    {
        r_sr += atoms[i].r;
    }
    r_sr /= N;
    double r_max = 2.5*r_sr;
    double dr = r_max/M;

    for (int it = 0; it<N; it++)
    {
        for (int it2 = it + 1; it2<N; it2++)
        {
            atom i = atoms[it];
            atom j = atoms[it2];
            double rij = sqrt(pow(j.x - i.x, 2) + pow(j.y - i.y, 2) + pow(j.z - i.z, 2));
            int m = floor(rij/dr);
            if(m<M) pcf[m] += (8*M_PI*r_sr*r_sr)/(N*N*2*M_PI*rij*dr);
        }
    }
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

    for (int i = 0; i<N; i++)
    {
        plik1>>atoms[i].x>>atoms[i].y>>atoms[i].z;
        atoms[i].x *= angstrom_to_meter;
        atoms[i].y *= angstrom_to_meter;
        atoms[i].z *= angstrom_to_meter;
        convertToSpherical(atoms[i]);  
    }

    plik1.close();

    for (int i = 0; i < N; i++)
    {
        cout << "Atom " << i << ": x=" << atoms[i].x << ", y=" << atoms[i].y << ", z=" << atoms[i].z << endl;
    }

    double test = calc_total_V(atoms);
    cout<<test/eV_to_joule;

    //Funkcja korealcji par
    ofstream plik2;
    plik2.open("C:\\Users\\kacpe\\OneDrive\\Pulpit\\C_plus\\Monte_Carlo\\Monte_carlo_8\\pcf_1.txt");
    int M = 100;
    double pcf[M];
    f_corelation(atoms, pcf, M);
    for (int i = 0; i<M; i++)
    {
        plik2<<i<<" "<<pcf[i]<<endl;
    }


    return 0;
}