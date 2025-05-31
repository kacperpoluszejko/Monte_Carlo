#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include<vector>
#include<complex>
#include <cmath>
#include <fstream>
using namespace std;


/*****************************************************************************
 *  liczymy rozwiazanie dokladne dla impulsu generowanego przez zrodlo
 *  oznaczenia:
 * 			x -polozenie, t - aktualny czas, t0 - maksimum sygnalu zrodla, Omega=2*Pi*f - (f to czestotliwosc zrodla)
 *                sigma - szerokosc impulsu,  R - opor linii, G - konduktancja linii, C-pojemnosc linii, Rg - opor zrodla
 *                Rl-opor cewki w x=l,length-dlugosc linii, number_nodes - liczba wezlow w calkowaniu po czestosci (omedze),
 *			n_sum_terms - maksymalna liczba wyrazow uwzgledniana w sumowaniu wkladow ui (i=0,1,...,n_sum_terms)
 * 
 * 
 *****************************************************************************/


double u_xt_exact(double x, double t, double t0,  double freq, double sigma, double R, double G, double L, double C, 
			    double Rg, double Rl, double length, int number_nodes, int n_sum_terms){
	
	double omega, domega, omega_min,omega_max,c_speed,aj,Omega;
	double sigma_w=1./sigma;
	double uxt=0.;
	complex<double> vg_m, I, z0, zl, ksi_omega, Gamma_l, Gamma_g,k_omega, u0_x_omega,u_sum_omega, multiplier,single_term, uxt_cmplx;
	int zakres;
	
	Omega=2.*M_PI*freq;
	zakres=4.;
	
	c_speed=1/sqrt(L*C); // predkosc swiatla
	I=0.0+1.0i;
	
	uxt_cmplx=0.+0.i;
	
	//wklady scentrowane w:  +/- Omega
	for(int m=-1;m<=1;m+=2){
		
		omega_min=m*Omega-zakres*sigma_w;
		omega_max=m*Omega+zakres*sigma_w;
		domega=(omega_max-omega_min)/number_nodes;
		
		
		//calkowanie po czestosci -(mala omega)
		for(int j=0;j<=number_nodes;j++){
			omega=omega_min+domega*j;
			vg_m=(-1.0*m)*I*sigma*sqrt(M_PI/2.)*exp(-pow((omega-m*Omega)/sigma_w,2)/2-I*(omega-m*Omega)*t0);
			z0=sqrt((R+I*omega*L)/(G+I*omega*C));
			ksi_omega=z0/(z0+Rg);
			
			zl=Rl+I*omega*0.; //tylko opor na wyjsciu
			
			Gamma_l=(zl-z0)/(zl+z0);
			Gamma_g=(Rg-z0)/(Rg+z0);
			k_omega=-I/c_speed*sqrt((I*omega+R/L)*(I*omega+G/C));
			u0_x_omega=ksi_omega*vg_m*( exp(-I*k_omega*x)+Gamma_l*exp(-I*k_omega*(2*length-x)) );
			
			//sumowanie wkladow: ui_x_omega
			u_sum_omega=0.+0.i;
			multiplier=Gamma_l*Gamma_g*exp(-2.*I*k_omega*length);
			single_term=1.;
			for(int ii=0;ii<=n_sum_terms;ii++){
				if(ii>0)single_term*=multiplier;
				u_sum_omega=u_sum_omega+u0_x_omega*single_term;
			}
			
			if(j==0 || j==number_nodes){
				aj=0.5;
			}else{
				aj=1.0;
			}	
			uxt_cmplx+=u_sum_omega/2./M_PI*exp(I*omega*t)*domega*aj;
			
		}//j
	}//m
	
	
	uxt=real(uxt_cmplx);
	
	
	return uxt;
}

const double L = 0.25 * pow(10, -6);
const double C = 100 * pow(10, -12);
const double R_0 = 12.5;
const double G = 0.5 * pow(10, -3);
const double l = 2;
const double R_l = 12.5;
const double R_g = 75;
const double v = 1 * pow(10, 9);
const double t0 = 7.5 * pow(10, -9);
const double sigma = 0.75 * pow(10, -9);
const double c = 1/(sqrt(L*C));
const double mu = G/C;
const double lambda = 0.5*(R_0/L + G/C);
const double tau = (R_0)/(R_0 + R_g);
const double tau_g = (R_g - R_0)/(R_0 + R_g);
const double tau_l = (R_l - R_0)/(R_0 + R_l);

int main()
{


    int npaths = 1000;
    double t_start = 50 * pow(10, -9);
    int Nx = 1000;
    double x_table [Nx] = {0};
    double dx = l/Nx;
    double fxt = 0, bxt = 0;
    double u_table[Nx] = {0};

    for (int nx = 0; nx<Nx; nx++)
    {
        double x_start = nx*dx;
        x_table[nx] = x_start;
        u_table[nx] =  u_xt_exact(x_start, t_start, t0, v, sigma, R_0, G, L, C, R_g, R_l, l, 1000, 100);

    }



    ofstream plik1;
    plik1.open("C:\\Users\\kacpe\\OneDrive\\Pulpit\\C_plus\\Monte_Carlo\\Monte_carlo_10\\monte_carlo_10_u5_teo.txt");
    for (int nx = 0; nx<Nx; nx++)
    {
        plik1<<u_table[nx]<<" "<<x_table[nx]<<endl;
    }

    


    return 0;
}

