#include <fstream>
#include <H5Cpp.h>
#include <string.h>
#include "shock.h"
#include "util.h"
#include "config.h"

using namespace H5;

//---------------------------------------------------------------------------
Shock::Shock(const Grid &grid, const Partition &part, const std::string &power) {
    if (hdf_input) {
	H5File hdf(power, H5F_ACC_RDONLY);

	{
	    DataSet ds = hdf.openDataSet("shock/time");
	    DataSpace space = ds.getSpace();

	    precondition(space.getSimpleExtentNdims() == 1,
		    std::string("Wrong format in file ") + std::string(power));

	    hsize_t n;
	    space.getSimpleExtentDims(&n);

	    time.resize(n);
	    ds.read(time.data(), PredType::NATIVE_DOUBLE);
	}

	std::vector<char> skip;

	{
	    DataSet ds = hdf.openDataSet("shock/coord");
	    DataSpace space = ds.getSpace();

	    precondition(space.getSimpleExtentNdims() == 2,
		    std::string("Wrong format in file ") + std::string(power));

	    hsize_t dim[2];
	    space.getSimpleExtentDims(dim);

	    precondition(dim[1] == 3,
		    std::string("Wrong format in file ") + std::string(power));

	    std::vector<std::array<int,3>> c(dim[0]);
	    ds.read(c.data(), PredType::NATIVE_INT32);

	    coord.reserve(dim[0]);
	    skip.resize(dim[0], false);

	    for(uint i = 0; i < c.size(); i++) {
		uint id = grid.idx(c[i][0], c[i][1], c[i][2]);
		if (part[id] == mpi_rank && !grid.boundary(id))
		    coord.push_back(id);
		else
		    skip[i] = true;
	    }
	}

	{
	    DataSet ds = hdf.openDataSet("shock/power");
	    DataSpace space = ds.getSpace();

	    precondition(space.getSimpleExtentNdims() == 3,
		    std::string("Wrong format in file ") + std::string(power));

	    hsize_t dim[3];
	    space.getSimpleExtentDims(dim);

	    precondition(time.size() == dim[0],
		    std::string("Wrong format in file ") + std::string(power));

	    precondition(dim[2] == 3,
		    std::string("Wrong format in file ") + std::string(power));

	    std::vector<double> f(dim[0] * dim[1] * dim[2]);
	    ds.read(f.data(), PredType::NATIVE_DOUBLE);

	    F.resize(time.size());
	    for(uint i = 0, idx = 0; i < time.size(); i++) {
		F[i].resize(coord.size());
		for(uint j = 0, k = 0; j < dim[1]; j++) {
		    if (skip[j]) {
			idx += 3;
		    } else {
			F[i][k][0] = f[idx++];
			F[i][k][1] = f[idx++];
			F[i][k][2] = f[idx++];

			k++;
		    }
		}
	    }
	}
    } else {
	uint n;

	std::ifstream fpower(power);

	precondition(fpower, "Failed to open power file.");

	double tmax, period;

	precondition(fpower >> tmax >> period >> n,
		"Wrong format for power file");

	coord.reserve(n);

	std::vector<char> skip(n, false);

	for(uint i = 0; i < n; i++) {
	    uint cx, cy, cz;
	    precondition(fpower >> cx >> cy >> cz, "Wrong format for power file");

	    uint id = grid.idx(cx, cy, cz);
	    if (part[id] == mpi_rank && !grid.boundary(id))
		coord.push_back(id);
	    else
		skip[i] = true;
	}

	for(double t = 0; t < tmax; t += period) {
	    time.emplace_back(t);
	    F.emplace_back(coord.size());

	    real fx, fy, fz;

	    for(uint i = 0, j = 0; i < n; i++) {
		precondition(fpower >> fx >> fy >> fz,
			"Wrong format for power file (not enough lines)");

		if (skip[i]) continue;

		F.back()[j][0] = fx;
		F.back()[j][1] = fy;
		F.back()[j][2] = fz;

		j++;
	    }
	}
    }
}
