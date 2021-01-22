#include <fstream>
#include <algorithm>
#include <H5Cpp.h>

#include "grid.h"
#include "util.h"
#include "config.h"

using namespace H5;

/*
 * [1] Голованов, Бережной. "Метод конечных элементов в механике деформируемых
 *     твердых тел"
 * [2] G. R. Liu, S. S. Quek. "The Finite Element Method. A Practical Course".
 */

static const float abscoef[100][2] = {
    {0.00059f, 0.006f},
    {0.00068f, 0.011f},
    {0.00077f, 0.016f},
    {0.00086f, 0.021f},
    {0.00095f, 0.026f},
    {0.00104f, 0.031f},
    {0.00113f, 0.036f},
    {0.00122f, 0.041f},
    {0.00131f, 0.046f},
    {0.0014f, 0.051f},
    {0.00149f, 0.056f},
    {0.00158f, 0.061f},
    {0.00167f, 0.066f},
    {0.00176f, 0.071f},
    {0.00185f, 0.076f},
    {0.00194f, 0.081f},
    {0.00203f, 0.086f},
    {0.00212f, 0.091f},
    {0.00221f, 0.096f},
    {0.0023f, 0.101f},
    {0.00239f, 0.106f},
    {0.00248f, 0.111f},
    {0.00257f, 0.116f},
    {0.00266f, 0.121f},
    {0.00275f, 0.126f},
    {0.00284f, 0.131f},
    {0.00293f, 0.231f},
    {0.00302f, 0.331f},
    {0.00311f, 0.431f},
    {0.0032f, 0.531f},
    {0.00329f, 0.631f},
    {0.00338f, 0.731f},
    {0.00347f, 0.831f},
    {0.00356f, 0.931f},
    {0.00365f, 1.031f},
    {0.00374f, 1.131f},
    {0.00383f, 1.231f},
    {0.00392f, 1.331f},
    {0.00401f, 1.431f},
    {0.0041f, 1.531f},
    {0.00419f, 1.631f},
    {0.00428f, 1.731f},
    {0.00437f, 1.831f},
    {0.00446f, 1.931f},
    {0.00455f, 2.031f},
    {0.00464f, 2.131f},
    {0.00473f, 2.231f},
    {0.00482f, 2.331f},
    {0.00491f, 2.431f},
    {0.005f, 2.531f},
    {0.00509f, 2.631f},
    {0.00518f, 3.231f},
    {0.00527f, 3.831f},
    {0.00536f, 4.431f},
    {0.00545f, 5.031f},
    {0.00554f, 5.631f},
    {0.00563f, 6.231f},
    {0.00572f, 6.831f},
    {0.00581f, 7.431f},
    {0.0059f, 8.031f},
    {0.00599f, 8.631f},
    {0.00608f, 9.231f},
    {0.00617f, 9.831f},
    {0.00626f, 10.431f},
    {0.00635f, 11.031f},
    {0.00644f, 11.631f},
    {0.00653f, 12.231f},
    {0.00662f, 12.831f},
    {0.00671f, 13.431f},
    {0.0068f, 14.031f},
    {0.00689f, 14.631f},
    {0.00698f, 15.231f},
    {0.00707f, 15.831f},
    {0.00716f, 16.431f},
    {0.00725f, 17.031f},
    {0.00734f, 17.631f},
    {0.00743f, 19.631f},
    {0.00752f, 21.631f},
    {0.00761f, 23.631f},
    {0.0077f, 25.631f},
    {0.00779f, 27.631f},
    {0.00788f, 29.631f},
    {0.00797f, 31.631f},
    {0.00806f, 33.631f},
    {0.00815f, 35.631f},
    {0.00824f, 37.631f},
    {0.00833f, 39.631f},
    {0.00842f, 41.631f},
    {0.00851f, 43.631f},
    {0.0086f, 45.631f},
    {0.00869f, 47.631f},
    {0.00878f, 53.631f},
    {0.00887f, 58.613f},
    {0.00896f, 63.595f},
    {0.00905f, 68.577f},
    {0.00914f, 73.559f},
    {0.00923f, 78.541f},
    {0.00932f, 83.523f},
    {0.00941f, 88.505f},
    {0.0095f, 93.487f}
};

//---------------------------------------------------------------------------
Grid::Grid(const std::string &geometry, const std::string &properties)
    : nx(0), ny(0), nz(0)
{
    // Геометрия.
    {
	std::ifstream fgeom(geometry);

	precondition(fgeom, "Failed to open geometry file");

	precondition(fgeom >> nx >> ny >> nz,
		"Wrong format for geometry file");

	x.resize(nx + 1);
	y.resize(ny + 1);
	z.resize(nz + 1);

	for(uint i = 0; i <= nx; i++)
	    precondition(fgeom >> x[i], "Wrong format for geometry file");

	for(uint i = 0; i <= ny; i++)
	    precondition(fgeom >> y[i], "Wrong format for geometry file");

	for(uint i = 0; i <= nz; i++)
	    precondition(fgeom >> z[i], "Wrong format for geometry file");
    }

    // Свойства.
    {
	std::ifstream fprop(properties);

	precondition(fprop, "Failed to open properties file");

	while(fprop) {
	    Domain d;

	    if (fprop >> d.bbox.lo[0] >> d.bbox.lo[1] >> d.bbox.lo[2]
		      >> d.bbox.hi[0] >> d.bbox.hi[1] >> d.bbox.hi[2]
		      >> d.prop.young
		      >> d.prop.poisson
		      >> d.prop.density
		      >> d.prop.stiffness
		      >> d.prop.mass
		      >> d.prop.viscosity
	       )
            {
                domain.push_back(d);
            }
	}

	precondition(domain.size(), "Wrong format for properties file");
    }

    nx_1 = nx + 1;
    nx_ny_1 = (nx + 1) * (ny + 1);

    nnodes = (nx + 1) * (ny + 1) * (nz + 1);
    nelems = nx * ny * nz;

    fill_properties();
}

//---------------------------------------------------------------------------
Grid::Grid(const std::string &h5name) : nx(0), ny(0), nz(0) {
    H5File hdf(h5name, H5F_ACC_RDONLY);

    {
	DataSet ds = hdf.openDataSet("x");
	DataSpace space = ds.getSpace();

	precondition(space.getSimpleExtentNdims() == 1,
		std::string("Wrong format in file ") + std::string(h5name));

	hsize_t dim;
	space.getSimpleExtentDims(&dim);

	nx = dim - 1;

	x.resize(dim);
	ds.read(x.data(), PredType::NATIVE_DOUBLE);
    }

    {
	DataSet ds = hdf.openDataSet("y");
	DataSpace space = ds.getSpace();

	precondition(space.getSimpleExtentNdims() == 1,
		std::string("Wrong format in file ") + std::string(h5name));

	hsize_t dim;
	space.getSimpleExtentDims(&dim);

	ny = dim - 1;

	y.resize(dim);
	ds.read(y.data(), PredType::NATIVE_DOUBLE);
    }

    {
	DataSet ds = hdf.openDataSet("z");
	DataSpace space = ds.getSpace();

	precondition(space.getSimpleExtentNdims() == 1,
		std::string("Wrong format in file ") + std::string(h5name));

	hsize_t dim;
	space.getSimpleExtentDims(&dim);

	nz = dim - 1;

	z.resize(dim);
	ds.read(z.data(), PredType::NATIVE_DOUBLE);
    }

    {
	DataSet ds = hdf.openDataSet("prop/box");
	DataSpace space = ds.getSpace();

	precondition(space.getSimpleExtentNdims() == 2,
		std::string("Wrong format in file ") + std::string(h5name));

	hsize_t dim[2];
	space.getSimpleExtentDims(dim);

	precondition(dim[1] == 6,
		std::string("Wrong format in file ") + std::string(h5name));

	std::vector<int> box(dim[0] * dim[1]);
	ds.read(box.data(), PredType::NATIVE_INT32);

	domain.resize(dim[0]);

	for(uint i= 0; i < dim[0]; i++) {
	    domain[i].bbox.lo[0] = box[i * 6];
	    domain[i].bbox.lo[1] = box[i * 6 + 1];
	    domain[i].bbox.lo[2] = box[i * 6 + 2];
	    domain[i].bbox.hi[0] = box[i * 6 + 3];
	    domain[i].bbox.hi[1] = box[i * 6 + 4];
	    domain[i].bbox.hi[2] = box[i * 6 + 5];
	}
    }

#define READ_PROP(name, field) \
    { \
	DataSet ds = hdf.openDataSet(name); \
	DataSpace space = ds.getSpace(); \
	precondition(space.getSimpleExtentNdims() == 1, \
		std::string("Wrong format in file ") + std::string(h5name)); \
	hsize_t dim; \
	space.getSimpleExtentDims(&dim); \
	precondition(dim == domain.size(), \
		std::string("Wrong format in file ") + std::string(h5name)); \
	std::vector<double> prop(dim); \
	ds.read(prop.data(), PredType::NATIVE_DOUBLE); \
	for(uint i = 0; i < dim; i++) domain[i].prop.field = prop[i]; \
    }

    READ_PROP("prop/E",     young);
    READ_PROP("prop/mu",    poisson);
    READ_PROP("prop/rho",   density);
    READ_PROP("prop/alpha", stiffness);
    READ_PROP("prop/beta",  mass);
    READ_PROP("prop/visc",  viscosity);

    nx_1 = nx + 1;
    nx_ny_1 = (nx + 1) * (ny + 1);

    nnodes = (nx + 1) * (ny + 1) * (nz + 1);
    nelems = nx * ny * nz;

    fill_properties();
}

//---------------------------------------------------------------------------
Grid::Properties Grid::properties(uint i, uint j, uint k) const {
    std::array<int,3> p = {{static_cast<int>(i), static_cast<int>(j), static_cast<int>(k)}};
    BoundingBox b(p);

    if (absorb.overlaps(b)) {
	return P[i][j][k];
    } else {
	std::array<int,3> ref = p;
	int radius = std::numeric_limits<int>::max();

	for(int d = 0; d < 3; d++) {
	    if (p[d] < absorb.lo[d]) {
		ref[d] = absorb.lo[d];
	    } else if (p[d] >= absorb.hi[d]) {
		ref[d] = absorb.hi[d] - 1;
	    }
	}

	radius = std::abs(ref[0] - int(i));
	radius = std::max<int>(radius, std::abs(ref[1] - int(j)));
	radius = std::max<int>(radius, std::abs(ref[2] - int(k)));
	radius = std::min<int>(radius, 99);

	Properties prop = P[ref[0]][ref[1]][ref[2]];

	prop.stiffness += radius * (0.0095 - prop.stiffness) / 99;
	prop.mass      += abscoef[radius][1] - 5e-3;

	return prop;
    }
}

//---------------------------------------------------------------------------
std::array<uint,8> Grid::node(uint i, uint j, uint k) const {
    // Локальная нумерация узлов в соответствии с [2], стр. 209.
    std::array<uint,8> id = {{
	idx(i,     j,     k    ),
	idx(i + 1, j,     k    ),
	idx(i + 1, j + 1, k    ),
	idx(i,     j + 1, k    ),
	idx(i,     j,     k + 1),
	idx(i + 1, j,     k + 1),
	idx(i + 1, j + 1, k + 1),
	idx(i,     j + 1, k + 1)
    }};

    return id;
}

//---------------------------------------------------------------------------
void Grid::fill_properties() {
    P.resize(boost::extents[nx][ny][nz]);

    for(const auto &d : domain)
        for(uint i = d.bbox.lo[0]; i < d.bbox.hi[0]; ++i)
            for(uint j = d.bbox.lo[1]; j < d.bbox.hi[1]; ++j)
                for(uint k = d.bbox.lo[2]; k < d.bbox.hi[2]; ++k)
                    if (i < nx && j < ny && k < nz)
                        P[i][j][k] = d.prop;
}
