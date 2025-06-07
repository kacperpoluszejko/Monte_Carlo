#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <random>

using namespace std;

class PHOTON_DIFFUSION_2D {
private:
    class BEAM {
    public:
        double w;
        double x, y;
        double x_new, y_new;
        double rx, ry;
        int layer;
        bool alive = true;
        int length;
        vector<double> path;
    };

public:
    int nlayers;
    int nx, ny;
    double xmax, ymax, dx, dy;
    double x_source, dx_source;
    double x_detect, dx_detect;
    double p_min;
    double w_min;
    double rx0, ry0;

    int write_all_paths;
    int write_source_detection_paths;

    BEAM beam;
    double abs_specular;

    vector<vector<double>> absorption;
    vector<double> reflectance;
    vector<double> transmittance;
    vector<vector<double>> layers_data;

    PHOTON_DIFFUSION_2D();
    void init();
    double uniform();
    void single_path();
    void roulette();
    void calculate_new_position();
    void scatter_in_layer();
    void scatter_up_down_boundary();
    double sign(double);
    void segment_intersection(double, double, double, double, double, double, double, double,
                               double&, double&, int&);
    void write_paths_to_file();
    void save_absorption_to_file(const string& filename) const;
};

PHOTON_DIFFUSION_2D::PHOTON_DIFFUSION_2D() {
    write_all_paths = 0;
    write_source_detection_paths = 0;
    nlayers = 0;
    nx = 0;
    ny = 0;
    xmax = 0.;
    ymax = 0.;
    dx = 0;
    dy = 0;
    p_min = 0.1;
    w_min = 1.0E-4;
    abs_specular = 0.;
    rx0 = 0.;
    ry0 = 1.;

    int M = 20;
    layers_data.resize(M, vector<double>(M, 0.));
    for (int i = 0; i < M; i++) {
        layers_data[i][3] = 0.0;
        layers_data[i][4] = 1.0;
    }
}

void PHOTON_DIFFUSION_2D::init() {
    absorption.resize(nx + 1, vector<double>(ny + 1, 0));
    reflectance.resize(nx + 1, 0.);
    transmittance.resize(nx + 1, 0.);

    dx = xmax / nx;
    ymax = 0.;
    for (int i = 1; i <= nlayers; i++) ymax += layers_data[i][2];
    dy = ymax / ny;

    for (int i = 1; i <= nlayers + 1; i++) {
        layers_data[i][5] = layers_data[i - 1][6];
        layers_data[i][6] += layers_data[i][5] + layers_data[i][2];
    }
}

double PHOTON_DIFFUSION_2D::uniform() {
    return (double)rand() / RAND_MAX;
}

void PHOTON_DIFFUSION_2D::single_path() {
    if (dx_source >= 0. && x_source >= 0 && x_source <= xmax) {
        beam.x = x_source + dx_source / 2 * (2 * uniform() - 1);
        beam.y = 1.0E-10;
        beam.rx = rx0;
        beam.ry = ry0;
        beam.w = 1.0;
        beam.alive = true;
        beam.layer = 1;
        beam.path.clear();
        beam.path.push_back(beam.x);
        beam.path.push_back(beam.y);
        beam.path.push_back(beam.w);
    }

    while (beam.alive == true) {
        calculate_new_position();
        scatter_in_layer();
        scatter_up_down_boundary();
        roulette();
        beam.path.push_back(beam.x);
        beam.path.push_back(beam.y);
        beam.path.push_back(beam.w);
    }

    write_paths_to_file();
}

void PHOTON_DIFFUSION_2D::write_paths_to_file() {
    if (write_all_paths == 1) {
        FILE* fp = fopen("all_paths.dat", "a");
        fprintf(fp, "\n");
        for (int i = 0; i < beam.path.size(); i += 3) {
            double x = beam.path[i];
            double y = beam.path[i + 1];
            double w = beam.path[i + 2];
            fprintf(fp, "%15.5E  %15.5E   %15.5E \n", x, y, w);
        }
        fclose(fp);
    }
}

void PHOTON_DIFFUSION_2D::save_absorption_to_file(const string& filename) const {
    ofstream file(filename);
    for (int j = ny; j >= 0; --j) {
        for (int i = 0; i <= nx; ++i) {
            file << absorption[i][j] << " ";
        }
        file << "\n";
    }
    file.close();
}

int main() {
    PHOTON_DIFFUSION_2D ob;
    ob.xmax = 0.2;
    ob.x_source = 0.1;
    ob.dx_source = 0.0;
    ob.x_detect = 0.15;
    ob.dx_detect = 0.01;
    ob.nx = 100;
    ob.ny = 100;
    ob.rx0 = 0.0;
    ob.ry0 = 1.0;
    ob.nlayers = 3;

    ob.layers_data[1][0] = 1.;
    ob.layers_data[1][1] = 10.;
    ob.layers_data[1][2] = 0.02;
    ob.layers_data[1][3] = 0.75;
    ob.layers_data[1][4] = 1.3;

    ob.layers_data[2][0] = 1.;
    ob.layers_data[2][1] = 190.;
    ob.layers_data[2][2] = 0.02;
    ob.layers_data[2][3] = 0.75;
    ob.layers_data[2][4] = 1.0;

    ob.layers_data[3][0] = 10;
    ob.layers_data[3][1] = 90.;
    ob.layers_data[3][2] = 0.02;
    ob.layers_data[3][3] = 0.95;
    ob.layers_data[3][4] = 1.0;

    ob.init();
    int N = 10000;  // zwiększ liczbę fotonów
    ob.write_all_paths = 0;
    ob.write_source_detection_paths = 0;
    for (int k = 0; k < N; k++) {
        ob.single_path();
    }

    ob.save_absorption_to_file("absorption_map.dat");

    return 0;
}
