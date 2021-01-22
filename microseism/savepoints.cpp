#include <fstream>
#include <H5Cpp.h>
#include <string.h>
#include "savepoints.h"
#include "util.h"

using namespace H5;

//---------------------------------------------------------------------------
void read_save_points(
        const std::string &fname,
        std::vector< int > &id,
        std::vector< std::array<int, 3> > &c
	)
{
    std::string wrong_format = "Wrong format in file " + fname;

    if (fname.find(".h5") + 3 == fname.size()) {
	H5File hdf(fname, H5F_ACC_RDONLY);

	{
	    DataSet ds = hdf.openDataSet("savepoints/pids");
	    DataSpace space = ds.getSpace();

	    precondition(space.getSimpleExtentNdims() == 1, wrong_format);

	    hsize_t dim;
	    space.getSimpleExtentDims(&dim);

	    id.resize(dim);
	    ds.read(id.data(), PredType::NATIVE_INT32);
	}

	{
	    DataSet ds = hdf.openDataSet("savepoints/coord");
	    DataSpace space = ds.getSpace();

	    precondition(space.getSimpleExtentNdims() == 2, wrong_format);

	    hsize_t dim[2];
	    space.getSimpleExtentDims(dim);

	    precondition(dim[1] == 3, wrong_format);

	    c.resize(dim[0]);
	    ds.read(c.data(), PredType::NATIVE_INT32);
	}

        precondition(id.size() == c.size(), wrong_format);
    } else {
	uint n;

	std::ifstream fp(fname);

	precondition(fp, "Failed to open file " + fname);

	precondition(fp >> n, wrong_format);

        id.resize(n);
	c.resize(n);

	for(uint i = 0; i < n; i++)
            precondition(fp >> id[i], wrong_format);

	for(uint i = 0; i < n; i++)
	    precondition(fp >> c[i][0] >> c[i][1] >> c[i][2], wrong_format);
    }
}
