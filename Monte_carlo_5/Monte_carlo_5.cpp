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


struct Wynik
 {
    double srednia;
    double wariancja;
    double error;
};


double g_1(double a, double b)
{
    RandomNumberGenerator rn;
    double U1 = a + (b-a)*rn.get_random();
    return 1+tanh(U1);
}


double g_2(double a, double b)
{
    RandomNumberGenerator rn;
    double U1 = a + (b-a)*rn.get_random();
    return 1/(1+U1*U1);
}

double g_3(double a, double b)
{
    RandomNumberGenerator rn;
    double U1 = a + (b-a)*rn.get_random();
    return pow(cos(M_PI*U1), 10);
}


Wynik calc_C1(double x_min, double x_max, int N)
{
    double a=-3, b=3;
    double suma=0, suma_2=0;
    for (int i=0; i<N; i++)
    {
        double g = g_1(x_min, x_max);
        suma += g;
        suma_2 += g*g;
    }
    double srednia = suma*(b-a)/N;
    double srednia_2 = suma_2*(b-a)*(b-a)/N;
    double wariancja = (srednia_2 - srednia*srednia)/N;
    double error = sqrt(wariancja)/srednia;

    return {srednia, wariancja, error};
}

Wynik calc_C2(double x_min, double x_max, int N)
{
    double a=0, b=10;
    double suma=0, suma_2=0;
    for (int i=0; i<N; i++)
    {
        double g = g_2(x_min, x_max);
        suma += g;
        suma_2 += g*g;
    }
    double srednia = suma*(b-a)/N;
    double srednia_2 = suma_2*(b-a)*(b-a)/N;
    double wariancja = (srednia_2 - srednia*srednia)/N;
    double error = sqrt(wariancja)/srednia;

    return {srednia, wariancja, error};
}

Wynik calc_C3(double x_min, double x_max, int N)
{
    double a=0, b=1;
    double suma=0, suma_2=0;
    for (int i=0; i<N; i++)
    {
        double g = g_3(x_min, x_max);
        suma += g;
        suma_2 += g*g;
    }
    double srednia = suma*(b-a)/N;
    double srednia_2 = suma_2*(b-a)*(b-a)/N;
    double wariancja = (srednia_2 - srednia*srednia)/N;
    double error = sqrt(wariancja)/srednia;

    return {srednia, wariancja, error};
}

int main()
{
    
    RandomNumberGenerator rn;
    int N = 1000;
    double a1 = -3, b1 = 3, a2 = 0, b2 = 10, a3 = 0, b3 = 1;

    // // Metoda podstawowa
    // Wynik C1 = calc_C1(a1, b1, N);
    // Wynik C2 = calc_C2(a2, b2, N);
    // Wynik C3 = calc_C3(a3, b3, N);
    // cout<<C1.srednia<<" "<< sqrt(C1.wariancja)<<" "<<C1.error<<endl;
    // cout<<C2.srednia<<" "<< sqrt(C2.wariancja)<<" "<<C2.error<<endl;
    // cout<<C3.srednia<<" "<< sqrt(C3.wariancja)<<" "<<C3.error<<endl;

    // double random[N];
    
    // ofstream plik1;
    // plik1.open("C:\\Users\\kacpe\\OneDrive\\Pulpit\\C_plus\\Monte_Carlo\\Monte_carlo_5\\monte_carlo_5_1.txt");
    // for (int  i = 0 ; i<N; i++)
    // {
    //     random[i] =a1 + (b1-a1)*rn.get_random();
    //     plik1<<random[i]<<" ";
    // }


    //Metoda systematyczna
    int M = 10;
    double pm = 1.0/M;
    double Nm = pm*N;

    double dx1 = (b1 - a1)/M;
    double dx2 = (b2 - a2)/M;
    double dx3 = (b3 - a3)/M;

    // Wynik C1 = {0,0,0};
    // Wynik C2 = {0,0,0};
    // Wynik C3 = {0,0,0};

    // for (int i=0; i<M; i++)
    // {
    //     Wynik C_1m = calc_C1(a1+dx1*i, a1+(i+1)*dx1, Nm);
    //     Wynik C_2m = calc_C2(a2+dx2*i, a2+(i+1)*dx2, Nm);
    //     Wynik C_3m = calc_C3(a3+dx3*i, a3+(i+1)*dx3, Nm);

    //     C1.srednia += C_1m.srednia*pm;
    //     C1.wariancja +=C_1m.wariancja*pm*pm;
    //     C2.srednia += C_2m.srednia*pm;
    //     C2.wariancja +=C_2m.wariancja*pm*pm;
    //     C3.srednia += C_3m.srednia*pm;
    //     C3.wariancja +=C_3m.wariancja*pm*pm;
    // }

    // C1.error = sqrt(C1.wariancja)/C1.srednia;
    // C2.error = sqrt(C2.wariancja)/C2.srednia;
    // C3.error = sqrt(C3.wariancja)/C3.srednia;

    // cout<<C1.srednia<<" "<< sqrt(C1.wariancja)<<" "<<C1.error<<endl;
    // cout<<C2.srednia<<" "<< sqrt(C2.wariancja)<<" "<<C2.error<<endl;
    // cout<<C3.srednia<<" "<< sqrt(C3.wariancja)<<" "<<C3.error<<endl;

    // vector<double> random(N);

    // ofstream plik2;
    // plik2.open("C:\\Users\\kacpe\\OneDrive\\Pulpit\\C_plus\\Monte_Carlo\\Monte_carlo_5\\monte_carlo_5_2.txt");

    // for (int i = 0; i<M; i++)
    // {
    //     for (int j = 0; j<Nm; j++)
    //     {
    //         double value = a1+dx1*i + dx1*rn.get_random();
    //         plik2<<value<<" ";
    //     }
    // }



    //Metoda warstwowa

    double sigma_m_table[M];
    double N_m_table[M];

    double suma = 0;
    for (int i=0; i<M; i++)
    {
        sigma_m_table[i] = sqrt(calc_C1(a1+dx1*i, a1+(i+1)*dx1, Nm).wariancja*Nm);
        suma += sigma_m_table[i];
    }

    double suma2 = 0;
    for (int i=0; i<M; i++)
    {
        N_m_table[i] = sigma_m_table[i]*N/suma;
        cout<<N_m_table[i]<<endl;
        suma2 += N_m_table[i];
    }
    cout<<suma2<<endl;

    Wynik C1 = {0,0,0};
    Wynik C2 = {0,0,0};
    Wynik C3 = {0,0,0};

    for (int i=0; i<M; i++)
    {
        if (N_m_table[i]<1) N_m_table[i]=1;
        Wynik C_1m = calc_C1(a1+dx1*i, a1+(i+1)*dx1, N_m_table[i]);
        Wynik C_2m = calc_C2(a2+dx2*i, a2+(i+1)*dx2, N_m_table[i]);
        Wynik C_3m = calc_C3(a3+dx3*i, a3+(i+1)*dx3, N_m_table[i]);

        C1.srednia += C_1m.srednia*pm;
        C1.wariancja +=C_1m.wariancja*pm*pm;
        C2.srednia += C_2m.srednia*pm;
        C2.wariancja +=C_2m.wariancja*pm*pm;
        C3.srednia += C_3m.srednia*pm;
        C3.wariancja +=C_3m.wariancja*pm*pm;
    }

    C1.error = sqrt(C1.wariancja)/C1.srednia;
    C2.error = sqrt(C2.wariancja)/C2.srednia;
    C3.error = sqrt(C3.wariancja)/C3.srednia;

    cout<<C1.srednia<<" "<< sqrt(C1.wariancja)<<" "<<C1.error<<endl;
    cout<<C2.srednia<<" "<< sqrt(C2.wariancja)<<" "<<C2.error<<endl;
    cout<<C3.srednia<<" "<< sqrt(C3.wariancja)<<" "<<C3.error<<endl;


    ofstream plik3;
    plik3.open("C:\\Users\\kacpe\\OneDrive\\Pulpit\\C_plus\\Monte_Carlo\\Monte_carlo_5\\monte_carlo_5_3.txt");

    for (int i = 0; i<M; i++)
    {
        for (int j = 0; j<N_m_table[i]; j++)
        {
            double value = a1+dx1*i + dx1*rn.get_random();
            plik3<<value<<" ";
        }
    }


    return 0;
} 
       