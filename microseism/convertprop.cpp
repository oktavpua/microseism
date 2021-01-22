#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <libconfig.h++>
#include <H5Cpp.h>

#include "grid.h"
#include "partition.h"
#include "util.h"

using namespace std;
using namespace H5;

//---------------------------------------------------------------------------
void save_prop(const Grid &grid, const BoundingBox &sbox)
{
    H5File hdf("prop.h5", H5F_ACC_TRUNC);

    FloatType datatype(PredType::NATIVE_DOUBLE);
    datatype.setOrder(H5T_ORDER_LE);

    hsize_t nx = sbox.hi[0] - sbox.lo[0];
    hsize_t ny = sbox.hi[1] - sbox.lo[1];
    hsize_t nz = sbox.hi[2] - sbox.lo[2];

    {
	DataSpace space(1, &nx);
	hdf.createDataSet("x", datatype, space).write(&grid.x[sbox.lo[0]], datatype);
    }
    {
	DataSpace space(1, &ny);
	hdf.createDataSet("y", datatype, space).write(&grid.y[sbox.lo[1]], datatype);
    }
    {
	DataSpace space(1, &nz);
	hdf.createDataSet("z", datatype, space).write(&grid.z[sbox.lo[2]], datatype);
    }

    {
	hsize_t dim[3] = {nz - 1, ny - 1, nx - 1};

	DataSpace space(3, dim);

	vector<double> buf(dim[0] * dim[1] * dim[2]);

	ofstream xmf("prop.xmf");

	xmf << "<?xml version=\"1.0\" ?>\n"
	    << "<!DOCTYPE Xdmf SYSTEM \"Sdmf.dtd\" []>\n"
	    << "<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.2\">\n"
	    << "<Domain>\n"
	    << "<Grid Name=\"Prop\" GridType=\"Uniform\">\n"
	    << "<Topology TopologyType=\"3DRectMesh\" NumberOfElements=\""
	    << nz << " " << ny << " " << nx << "\"/>\n"
	    << "<Geometry GeometryType=\"VXVYVZ\">\n"
	    << "<DataItem Dimensions=\"" << nx << "\" NumberType=\"Float\" Precision=\""
	    << sizeof(double) << "\" Format=\"HDF\">prop.h5:/x</DataItem>\n"
	    << "<DataItem Dimensions=\"" << ny << "\" NumberType=\"Float\" Precision=\""
	    << sizeof(double) << "\" Format=\"HDF\">prop.h5:/y</DataItem>\n"
	    << "<DataItem Dimensions=\"" << nz << "\" NumberType=\"Float\" Precision=\""
	    << sizeof(double) << "\" Format=\"HDF\">prop.h5:/z</DataItem>\n"
	    << "</Geometry>\n";

#define WRITE_PROP(name, field) \
	for(uint k = sbox.lo[2], idx = 0; k < sbox.lo[2] + dim[0]; k++) \
	    for(uint j = sbox.lo[1]; j < sbox.lo[1] + dim[1]; j++) \
		for(uint i = sbox.lo[0]; i < sbox.lo[0] + dim[2]; i++, idx++) \
		    buf[idx] = grid.properties(i, j, k).field; \
	hdf.createDataSet(name, datatype, space).write(buf.data(), datatype); \
	xmf << "<Attribute Name=\"" << name \
	    << "\" AttributeType=\"Scalar\" Center=\"Cell\">\n" \
	    << "<DataItem Dimensions=\"" << nz << " " << ny << " " << nx \
	    << "\" NumberType=\"Float\" Precision=\"" << sizeof(double) \
	    << "\" Format=\"HDF\">prop.h5:/" << name << "</DataItem>\n" \
	    << "</Attribute>\n"

	WRITE_PROP("young",   young);
	WRITE_PROP("poisson", poisson);
	WRITE_PROP("density", density);
	WRITE_PROP("alpha",   stiffness);
	WRITE_PROP("beta",    mass);

	xmf << "</Grid>\n</Domain>\n</Xdmf>\n";
    }
}

//---------------------------------------------------------------------------
int main(int argc, char *argv[]) {
    try {
	if (argc < 5) {
	    cout << "Usage: " << argv[0]
		 << " <geometry> <properies> <power> <config>"
		 << endl;
	    return 1;
	}

	// Чтение сетки.
	cout << "Reading in \"" << argv[1] << "\" and \""
	     << argv[2] << "\"" << endl;
	Grid grid(argv[1], argv[2]);

	// Чтение конфигурации.
	BoundingBox sbox;
	sbox.lo[0] = 0;
	sbox.lo[1] = 0;
	sbox.lo[2] = 0;
	sbox.hi[0] = grid.nx + 1;
	sbox.hi[1] = grid.ny + 1;
	sbox.hi[2] = grid.nz + 1;

	grid.absorb = sbox;

	cout << "Parsing \"" << argv[4] << "\"" << endl;
	try {
	    libconfig::Config cfg;
	    cfg.setAutoConvert(true);
	    cfg.readFile(argv[4]);

	    if(cfg.exists("savebox")) {
		libconfig::Setting &bb = cfg.lookup("savebox");

		sbox.lo[0] = bb[0];
		sbox.lo[1] = bb[1];
		sbox.lo[2] = bb[2];
		sbox.hi[0] = bb[3];
		sbox.hi[1] = bb[4];
		sbox.hi[2] = bb[5];
	    }

	    if(cfg.exists("absorb")) {
		libconfig::Setting &bb = cfg.lookup("absorb");

		grid.absorb.lo[0] = bb[0];
		grid.absorb.lo[1] = bb[1];
		grid.absorb.lo[2] = bb[2];
		grid.absorb.hi[0] = bb[3];
		grid.absorb.hi[1] = bb[4];
		grid.absorb.hi[2] = bb[5];
	    }
	} catch(const libconfig::FileIOException&) {
	    precondition(false, "Failed to read configuration file");
	}
	cout << "  savebox    = " << sbox        << endl
	     << "  absorb     = " << grid.absorb << endl;

	// Сохранение свойств.
	save_prop(grid, sbox);
    } catch (const std::exception &e) {
	cout << "Error: " << e.what() << endl;
	return 1;
    }

    return 0;
}
