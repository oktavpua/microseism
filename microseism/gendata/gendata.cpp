#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

double rand_float(double a, double b) {
    return a + (b - a) * rand() / RAND_MAX;
}

//---------------------------------------------------------------------------
int main(int argc, char *argv[]) {
    int nx = (argc < 2 ? 20 : atoi(argv[1]));
    int ny = nx;
    int nz = nx;

    double dx = 50;
    double dy = 50;
    double dz = 50;

    // Геометрия.
    {
	ofstream f("geometry.txt");

	f << nx << " " << ny << " " << nz << endl;

	for(int i = 0; i <= nx; i++)
	    f << i * dx << " ";
	f << endl;

	for(int i = 0; i <= ny; i++)
	    f << i * dy << " ";
	f << endl;

	for(int i = 0; i <= nz; i++)
	    f << i * dz << " ";
	f << endl;
    }

    // Свойства.
    {
	ofstream f("property.txt");

	int beg = 3;
	int end = nx - 8;
	for(int k = 0; k < end - beg; k++) {
	    for(int j = beg; j < end; j++) {
		for(int i = beg; i < end; i++) {
		    f << i << " " << j << " " << k << " "
		      << i + 1 << " " << j + 1 << " " << k + 1 << " "
		      << rand_float(5e9, 1.5e10) << " "
		      << rand_float(2e-1, 4e-1) << " "
		      << rand_float(2e3, 4e3) << " "
		      //<< rand_float(7e2, 9e2) << " "
		      << rand_float(1e-3, 5e-3) << " "
		      << rand_float(1e-3, 1e-2) << " "
		      << 1 << endl;
		      //<< rand_float(4.5e-3, 5.5e-3) << " "
		      //<< rand_float(3e-1, 4e-1) << endl;
		}
	    }
	}
    }

    // Удары.
    {
	ofstream f("power.txt");

	f << "0.001 0.001 " << (nx + 1) * (ny + 1) << endl;

	for(int j = 0; j <= ny; j++)
	    for(int i = 0; i <= nx; i++)
		f << i << " " << j << " 0 ";
	f << endl;

	for(int j = 0; j <= ny; j++)
	    for(int i = 0; i <= nx; i++)
		f << "0 0 12000 ";
	f << endl;
    }
}
