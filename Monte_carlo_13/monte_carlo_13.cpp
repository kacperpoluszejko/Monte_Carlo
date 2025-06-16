#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <random>
#include "dsmc_2d.cpp"  


int main() {
    DSMC_2D ob;
    ob.read("i.dat");
    ob.init();
    ob.write_position_velocity("rv0.dat");
    ob.nthreads = 8;
    ob.icol = 1;
    ob.evolution(0.0, 2000);
    ob.write_position_velocity("rv.dat");
    return 0;
}
