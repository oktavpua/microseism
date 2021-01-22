#include <iostream>
#include <iomanip>
#include "partition.h"
#include "util.h"
#include "config.h"

//---------------------------------------------------------------------------
Partition::Partition(const Grid &grid)
    : nx(grid.nx + 1), ny(grid.ny + 1), nz(grid.nz + 1)
{
    bbox.reserve(mpi_size);

    BoundingBox bb;

    bb.lo[0] = 0;
    bb.lo[1] = 0;
    bb.lo[2] = 0;

    bb.hi[0] = grid.nx + 1;
    bb.hi[1] = grid.ny + 1;
    bb.hi[2] = grid.nz + 1;

    split(bb, mpi_size);

    for(uint i = 0; i < bbox.size(); i++)
	rtree.insert(bbox[i], i);
}

//---------------------------------------------------------------------------
void Partition::split(BoundingBox bb, int nparts) {
    if (nparts == 1) {
	bbox.push_back(bb);
	return;
    }

    int wdim = 0;
    for(int i = 1; i < 3; i++)
	if (dim_weight[i] * (bb.hi[i] - bb.lo[i]) >
		dim_weight[wdim] * (bb.hi[wdim] - bb.lo[wdim])
	   ) wdim = i;

    int middle = bb.lo[wdim] +
	(bb.hi[wdim] - bb.lo[wdim]) * (nparts / 2) / nparts;

    BoundingBox bb1 = bb;
    BoundingBox bb2 = bb;

    bb1.hi[wdim] = middle;
    bb2.lo[wdim] = middle;

    split(bb1, nparts / 2);
    split(bb2, nparts - nparts / 2);
}

//---------------------------------------------------------------------------
uint Partition::size(int p) const {
    return bbox[p].area();
}

//---------------------------------------------------------------------------
int Partition::operator[](uint id) const {
    std::array<int,3> p = {{
	static_cast<int>(id % nx),
	static_cast<int>((id / nx) % ny),
	static_cast<int>(id / (nx * ny))
    }};
    return rtree.pick(p);
}

//---------------------------------------------------------------------------
BoundingBox Partition::box(int p) const {
    return bbox[p];
}

//---------------------------------------------------------------------------
uint Partition::g2l(uint id) const {
    BoundingBox bb = bbox[mpi_rank];

    uint i = id % nx;
    uint j = (id / nx) % ny;
    uint k = id / (nx * ny);

    i -= bb.lo[0];
    j -= bb.lo[1];
    k -= bb.lo[2];

    uint dx = bb.hi[0] - bb.lo[0];
    uint dy = bb.hi[1] - bb.lo[1];

    return i + j * dx + k * dx * dy;
}

//---------------------------------------------------------------------------
uint Partition::l2g(uint id) const {
    BoundingBox bb = bbox[mpi_rank];

    uint dx = bb.hi[0] - bb.lo[0];
    uint dy = bb.hi[1] - bb.lo[1];

    uint i = id % dx;
    uint j = (id / dx) % dy;
    uint k = id / (dx * dy);

    i += bb.lo[0];
    j += bb.lo[1];
    k += bb.lo[2];

    return i + j * nx + k * nx * ny;
}

//---------------------------------------------------------------------------
std::ostream& operator<<(std::ostream &os, const Partition &part) {
    for(auto bb = part.bbox.begin(); bb != part.bbox.end(); bb++)
	os << *bb << std::endl;
    return os;
}

