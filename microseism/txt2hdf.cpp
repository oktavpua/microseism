#include <fstream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <H5Cpp.h>
#include "util.h"
using namespace std;
using namespace H5;

int main(int argc, char *argv[]) {
    if (argc < 5) {
	std::cerr << "Usage: " << argv[0]
	    << " <geometry> <properties> <power> <out.h5>"
	    << std::endl;
	return 1;
    }

    try {
	H5File hdf(argv[4], H5F_ACC_TRUNC);

	hsize_t nx, ny, nz;

	/* Geometry */
	{
	    ifstream f(argv[1]);
	    precondition(f, "Failed to open geometry file");
	    
	    precondition(
		    f >> nx >> ny >> nz,
		    "Wrong geometry format (nx, ny, nz)"
		    ); 

	    vector<double> x(++nx);
	    vector<double> y(++ny);
	    vector<double> z(++nz);

	    for(int i = 0; i < nx; i++)
		precondition(f >> x[i], "Wrong geometry format (x)");

	    for(int i = 0; i < ny; i++)
		precondition(f >> y[i], "Wrong geometry format (y)");

	    for(int i = 0; i < nz; i++)
		precondition(f >> z[i], "Wrong geometry format (z)");

	    hdf.createDataSet(
		    "x", PredType::NATIVE_DOUBLE, DataSpace(1, &nx)
		    ).write(x.data(), PredType::NATIVE_DOUBLE);

	    hdf.createDataSet(
		    "y", PredType::NATIVE_DOUBLE, DataSpace(1, &ny)
		    ).write(y.data(), PredType::NATIVE_DOUBLE);

	    hdf.createDataSet(
		    "z", PredType::NATIVE_DOUBLE, DataSpace(1, &nz)
		    ).write(z.data(), PredType::NATIVE_DOUBLE);
	}

	/* Properties */
	{
	    ifstream f(argv[2]);
	    precondition(f, "Failed to open properties file");

	    vector<int>    box;
	    vector<double> E;
	    vector<double> mu;
	    vector<double> rho;
	    vector<double> alpha;
	    vector<double> beta;
	    vector<double> visc;

	    while(f) {
		int lo[3], hi[3];
		struct {
		    double E;
		    double mu;
		    double rho;
		    double alpha;
		    double beta;
		    double visc;
		} tmp;

		f >> lo[0] >> lo[1] >> lo[2]
		  >> hi[0] >> hi[1] >> hi[2]
		  >> tmp.E >> tmp.mu >> tmp.rho
		  >> tmp.alpha >> tmp.beta >> tmp.visc;

		if (f) {
		    box.push_back(lo[0]);
		    box.push_back(lo[1]);
		    box.push_back(lo[2]);
		    box.push_back(hi[0]);
		    box.push_back(hi[1]);
		    box.push_back(hi[2]);

		    E.push_back(tmp.E);
		    mu.push_back(tmp.mu);
		    rho.push_back(tmp.rho);
		    alpha.push_back(tmp.alpha);
		    beta.push_back(tmp.beta);
		    visc.push_back(tmp.visc);
		}
	    }

	    precondition(box.size() && (box.size() % 6 == 0),
		    "Wrong properties format (boxes)"
		    );

	    hdf.createGroup("prop");

	    hsize_t count[] = {E.size(), 6};

	    hdf.createDataSet(
		    "prop/box", PredType::NATIVE_INT32, DataSpace(2, count)
		    ).write(box.data(), PredType::NATIVE_INT32);

	    hdf.createDataSet(
		    "prop/E", PredType::NATIVE_DOUBLE, DataSpace(1, count)
		    ).write(E.data(), PredType::NATIVE_DOUBLE);

	    hdf.createDataSet(
		    "prop/mu", PredType::NATIVE_DOUBLE, DataSpace(1, count)
		    ).write(mu.data(), PredType::NATIVE_DOUBLE);

	    hdf.createDataSet(
		    "prop/rho", PredType::NATIVE_DOUBLE, DataSpace(1, count)
		    ).write(rho.data(), PredType::NATIVE_DOUBLE);

	    hdf.createDataSet(
		    "prop/alpha", PredType::NATIVE_DOUBLE, DataSpace(1, count)
		    ).write(alpha.data(), PredType::NATIVE_DOUBLE);

	    hdf.createDataSet(
		    "prop/beta", PredType::NATIVE_DOUBLE, DataSpace(1, count)
		    ).write(beta.data(), PredType::NATIVE_DOUBLE);

	    hdf.createDataSet(
		    "prop/visc", PredType::NATIVE_DOUBLE, DataSpace(1, count)
		    ).write(visc.data(), PredType::NATIVE_DOUBLE);
	}

	/* Shocks */
	{
	    ifstream f(argv[3]);
	    precondition(f, "Failed to open properties file");

	    double tmax, period;
	    size_t n;

	    precondition(
		    f >> tmax >> period >> n,
		    "Wrong format for power file (tmax, period, n)"
		    );

	    vector<int> coord(3 * n);

	    for(auto c = coord.begin(); c != coord.end(); c += 3)
		precondition(
			f >> c[0] >> c[1] >> c[2],
			"Wrong format for power file (coord)"
			);

	    istream_iterator<double> ifb(f);
	    istream_iterator<double> ife;
	    vector<double> power(ifb, ife);

	    precondition(power.size() && (power.size() % (3 * n) == 0),
		    "Wrong format for power file (power)"
		    );

	    vector<double> time(power.size() / (3 * n));
	    for(int i = 0; i < time.size(); i++)
		time[i] = i * period;

	    hdf.createGroup("shock");

	    {
		hsize_t count[] = {n, 3};
		hdf.createDataSet(
			"shock/coord", PredType::NATIVE_INT32, DataSpace(2, count)
			).write(coord.data(), PredType::NATIVE_INT32);
	    }

	    {
		hsize_t count[] = {time.size()};
		hdf.createDataSet(
			"shock/time", PredType::NATIVE_DOUBLE, DataSpace(1, count)
			).write(time.data(), PredType::NATIVE_DOUBLE);
	    }

	    {
		hsize_t count[] = {time.size(), n, 3};
		hdf.createDataSet(
			"shock/power", PredType::NATIVE_DOUBLE, DataSpace(3, count)
			).write(power.data(), PredType::NATIVE_DOUBLE);
	    }
	}
    } catch(const std::exception &e) {
	std::cerr << "Error: " << e.what() << std::endl;
    }
}
