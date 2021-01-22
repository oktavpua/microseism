#include <iostream>
#include <unordered_map>
#include <cassert>
#include "buildvoigt.h"

std::hash<uint> VoigtModel::HashRow::hi;
std::hash<real>   VoigtModel::HashRow::hd;

const VoigtModel::HashRow::mask_type VoigtModel::HashRow::mask = hash_mask;
const real   VoigtModel::HashRow::eps  = hash_eps;
const real   VoigtModel::CompRow::eps  = hash_eps;

//---------------------------------------------------------------------------
VoigtModel::VoigtModel(const Grid &grid, const Partition &part,
	int mpi_rank, real tau, int first_plane, bool skip_first
	)
    : grid(grid), part(part), mpi_rank(mpi_rank), tau(tau),
      top_ready(!skip_first), unique_idx(1024, HashRow(), CompRow())
{
    auto bbox = part.box(mpi_rank);

    uint nx = bbox.hi[0] - bbox.lo[0];
    uint ny = bbox.hi[1] - bbox.lo[1];
    uint nz = bbox.hi[2] - bbox.lo[2];

    plane_size  = 3 * nx * ny;
    plane_start = first_plane * plane_size;

    top.resize(plane_size);
    bot.resize(plane_size);
    ptr.resize(plane_size * nz);
}

//---------------------------------------------------------------------------
void VoigtModel::add(uint i, real m) {
    assert(i >= plane_start && i < plane_start + 2 * plane_size);

    i -= plane_start;

    if (i < plane_size) {
	top[i].add(m);
    } else {
	bot[i - plane_size].add(m);
    }
}

//---------------------------------------------------------------------------
void VoigtModel::add(uint i, int j, real a1, real a2) {
    assert(i >= plane_start && i < plane_start + 2 * plane_size);

    i -= plane_start;

    if (i < plane_size) {
	top[i].add(j, a1, a2);
    } else {
	bot[i - plane_size].add(j, a1, a2);
    }
}

//---------------------------------------------------------------------------
void VoigtModel::shift_plane() {
    // Добавляем строки из верхней плоскости в набор уникальных строк.
    if (top_ready) {
	for(uint i = 0; i < plane_size; i++) {
	    top[i].preprocess(tau, grid.boundary(part.l2g((plane_start + i) / 3)));

	    auto found = unique_idx.find(&top[i]);

	    if (found == unique_idx.end()) {
		unique_rows.push_back(top[i]);
		unique_idx.insert(&unique_rows.back());

		ptr[plane_start + i] = &unique_rows.back();
	    } else {
		ptr[plane_start + i] = *found;
	    }
	}
    }

    std::swap(top, bot);

    for(auto r = bot.begin(); r != bot.end(); r++) {
	r->M   = 0;
	r->len = 0;
    }

    plane_start += plane_size;
    top_ready = true;
}

//---------------------------------------------------------------------------
VoigtModel::Matrices VoigtModel::matrices() const {
    Matrices mtx;

    mtx.idx.resize(ptr.size());

    mtx.row.resize(unique_rows.size() + 1);
    mtx.M.resize(unique_rows.size());

    uint nnz = 0;
    for(auto r = unique_rows.begin(); r != unique_rows.end(); r++)
	nnz += r->len;

    mtx.col.reserve(nnz);
    mtx.A1.reserve(nnz);
    mtx.A2.reserve(nnz);

    std::unordered_map<const FullRow*, uint> ptr2idx(2 * unique_rows.size());

    for(uint i = 0; i < unique_rows.size(); i++)
	ptr2idx[&unique_rows[i]] = i;

    for(uint i = 0; i < ptr.size(); i++)
	mtx.idx[i] = ptr2idx[ptr[i]];

    mtx.row[0] = 0;
    {
	uint i = 0;
	auto   r = unique_rows.begin();

	for(; r != unique_rows.end(); r++, i++) {
	    mtx.M[i] = r->M;

	    for(int j = 0; j < r->len; j++) {
		mtx.col.push_back(r->col[j]);
		mtx.A1.push_back(r->A1[j]);
		mtx.A2.push_back(r->A2[j]);
	    }

	    mtx.row[i + 1] = mtx.row[i] + r->len;
	}
    }

    return mtx;
}

//---------------------------------------------------------------------------
VoigtModel::Matrices VoigtModel::remote_matrices() const {
    Matrices mtx;

    mtx.row.resize(ptr.size() + 1);

    uint nnz = 0;

    for(uint i = 0; i < ptr.size(); i++)
	nnz += ptr[i]->len;

    mtx.col.reserve(nnz);
    mtx.A1.reserve(nnz);
    mtx.A2.reserve(nnz);

    mtx.row[0] = 0;
    for(uint i = 0; i < ptr.size(); i++) {
	for(int j = 0; j < ptr[i]->len; j++) {
	    mtx.col.push_back(ptr[i]->col[j] + 3 * part.l2g(i/3) + (i % 3));

	    mtx.A1.push_back(ptr[i]->A1[j]);
	    mtx.A2.push_back(ptr[i]->A2[j]);
	}

	mtx.row[i + 1] = mtx.row[i] + ptr[i]->len;
    }

    return mtx;
}
