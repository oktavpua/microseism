#include <fstream>
#include <sstream>
#include <algorithm>
#include <boost/filesystem.hpp>

#include "saver.h"

using namespace H5;

//---------------------------------------------------------------------------
HDFSaver::HDFSaver(const Grid &grid, const Partition &part,
        const std::vector<int> &sp_ids,
        const std::vector<std::array<int,3>> &sp_coords,
	const std::string &fmask
	)
    : part(part), hdfname(mpi_size), do_savebox(mpi_size), do_savepnt(mpi_size),
#ifdef DOUBLE_PRECISION
      datatype(PredType::NATIVE_DOUBLE),
#else
      datatype(PredType::NATIVE_FLOAT),
#endif
      nx(mpi_size), ny(mpi_size), nz(mpi_size), np(mpi_size),
      step(0), last_checkpoint(std::chrono::high_resolution_clock::now())
{
    std::vector<int>  pids;
    std::vector<real> pcoo;
    std::vector<std::array<int,3>> pi;

    for(int i = 0; i < mpi_size; i++) {
        // File name
        std::ostringstream fname;
        fname << fmask << "-" << i << ".h5";
        hdfname[i] = fname.str();


        // Data in a save_box
        BoundingBox sb = part.box(i);

        do_savebox[i] = sb.overlaps(save_box) && save_box.area();

        np[i] = std::count_if(sp_coords.begin(), sp_coords.end(),
                [&sb](const std::array<int, 3> &p){ return sb.owns(p); }
                );

        do_savepnt[i] = (np[i] > 0);

        sb.lo[0] = std::max(sb.lo[0], save_box.lo[0]);
        sb.lo[1] = std::max(sb.lo[1], save_box.lo[1]);
        sb.lo[2] = std::max(sb.lo[2], save_box.lo[2]);

        sb.hi[0] = std::min(sb.hi[0], save_box.hi[0]);
        sb.hi[1] = std::min(sb.hi[1], save_box.hi[1]);
        sb.hi[2] = std::min(sb.hi[2], save_box.hi[2]);

        nx[i] = sb.hi[0] - sb.lo[0];
        ny[i] = sb.hi[1] - sb.lo[1];
        nz[i] = sb.hi[2] - sb.lo[2];

        if (i == mpi_rank) sbox = sb;
    }

    if (do_savepnt[mpi_rank]) {
        BoundingBox bb = part.box(mpi_rank);

        for(size_t i = 0; i < sp_coords.size(); ++i) {
            auto p = sp_coords[i];
            if (bb.owns(p)) {
                pcoo.push_back( grid.x[p[0]] );
                pcoo.push_back( grid.y[p[1]] );
                pcoo.push_back( grid.z[p[2]] );

                p[0] -= bb.lo[0];
                p[1] -= bb.lo[1];
                p[2] -= bb.lo[2];

                pi.push_back( p );

                pids.push_back(sp_ids[i]);
            }
        }
    }

    xmf_name = fmask + ".xmf";
    std::string full_hdfname = scratch   + "/" + hdfname[mpi_rank];
    std::string chkpt_name   = chkpt_dir + "/" + hdfname[mpi_rank];

    if (do_savebox[mpi_rank]) {
        {
            hsize_t count[] = {nz[mpi_rank], ny[mpi_rank], nx[mpi_rank], 3};
            fbspace = DataSpace(4, count);
        }

        if (compress) {
            hsize_t chunk_size[] = {
                std::min<uint>(32, nz[mpi_rank]),
                std::min<uint>(32, ny[mpi_rank]),
                std::min<uint>(32, nx[mpi_rank]),
                3};
            prop.setChunk(4, chunk_size);
            prop.setDeflate(9);
        }

        BoundingBox bb = part.box(mpi_rank);

        hsize_t count[] = {
            static_cast<hsize_t>(bb.hi[2] - bb.lo[2]),
            static_cast<hsize_t>(bb.hi[1] - bb.lo[1]),
            static_cast<hsize_t>(bb.hi[0] - bb.lo[0]),
            3};

        mbspace = DataSpace(4, count);

        hsize_t start[] = {
            static_cast<hsize_t>(sbox.lo[2] - part.box(mpi_rank).lo[2]),
            static_cast<hsize_t>(sbox.lo[1] - part.box(mpi_rank).lo[1]),
            static_cast<hsize_t>(sbox.lo[0] - part.box(mpi_rank).lo[0]),
            0};

        count[0] = nz[mpi_rank];
        count[1] = ny[mpi_rank];
        count[2] = nx[mpi_rank];

        mbspace.selectHyperslab(H5S_SELECT_SET, count, start);
    }

    if (do_savepnt[mpi_rank]) {
        {
            BoundingBox bb = part.box(mpi_rank);

            hsize_t count[] = {
                static_cast<hsize_t>(bb.hi[2] - bb.lo[2]),
                static_cast<hsize_t>(bb.hi[1] - bb.lo[1]),
                static_cast<hsize_t>(bb.hi[0] - bb.lo[0]),
                3};

            mpspace = DataSpace(4, count);

            std::vector<hsize_t> coords;
            coords.reserve(pi.size() * 12);
            for(const auto &p : pi) {
                for(int i = 0; i < 3; ++i) {
                    coords.push_back(p[2]);
                    coords.push_back(p[1]);
                    coords.push_back(p[0]);
                    coords.push_back(i);
                }
            }

            mpspace.selectElements(H5S_SELECT_SET, pi.size() * 3, coords.data());
        }

        {
            hsize_t count[] = {pi.size(), 3};
            fpspace = DataSpace(2, count);
        }
    }

    if (!boost::filesystem::exists(chkpt_dir))
        boost::filesystem::create_directory(chkpt_dir);

    if (
            boost::filesystem::exists(xmf_name)     &&
            boost::filesystem::exists(full_hdfname) &&
            boost::filesystem::exists(chkpt_name)
       )
    {
        hdf       = H5File(full_hdfname, H5F_ACC_RDWR);
        hdf_chkpt = H5File(chkpt_name,   H5F_ACC_RDWR);

        chkpt_step = hdf_chkpt.openDataSet("checkpoint/step");
        chkpt_time = hdf_chkpt.openDataSet("checkpoint/time");
        chkpt_q0   = hdf_chkpt.openDataSet("checkpoint/q0");
        chkpt_q1   = hdf_chkpt.openDataSet("checkpoint/q1");
        time_ds    = hdf.openDataSet("time");

        xmf.load_file(xmf_name.c_str());
        time_collection = xmf.child("Xdmf").child("Domain").find_child_by_attribute("Grid", "Name", "Time");

        int last_step;
        chkpt_step.read(&last_step, PredType::NATIVE_INT32);

        std::vector<pugi::xml_node> nodes_to_remove;

        for(auto grid = time_collection.child("Grid"); grid; grid = grid.next_sibling("Grid"))
            if (++step > (uint)last_step)
                nodes_to_remove.push_back(grid);

        step = last_step;

        for(auto n = nodes_to_remove.begin(); n != nodes_to_remove.end(); ++n)
            time_collection.remove_child(*n);

        have_chkpt = true;
    } else {
        if (do_savebox[mpi_rank] || do_savepnt[mpi_rank])
            hdf = H5File(full_hdfname, H5F_ACC_TRUNC);

        if (do_savebox[mpi_rank]) {
            {
                DataSpace space(1, &nx[mpi_rank]);
                DataSet   ds = hdf.createDataSet("x", datatype, space);
                ds.write(&grid.x[sbox.lo[0]], datatype);
            }
            {
                DataSpace space(1, &ny[mpi_rank]);
                DataSet   ds = hdf.createDataSet("y", datatype, space);
                ds.write(&grid.y[sbox.lo[1]], datatype);
            }
            {
                DataSpace space(1, &nz[mpi_rank]);
                DataSet   ds = hdf.createDataSet("z", datatype, space);
                ds.write(&grid.z[sbox.lo[2]], datatype);
            }

            {
                if (nx[mpi_rank] > 1 && ny[mpi_rank] > 1 && nz[mpi_rank] > 1) {
                    hsize_t start[3] = {
                        static_cast<hsize_t>(sbox.lo[2]),
                        static_cast<hsize_t>(sbox.lo[1]),
                        static_cast<hsize_t>(sbox.lo[0])
                    };

                    hsize_t dim[3] = {
                        nz[mpi_rank] - 1,
                        ny[mpi_rank] - 1,
                        nx[mpi_rank] - 1
                    };

                    DataSpace space(3, dim);

                    std::vector<real> buf(dim[0] * dim[1] * dim[2]);
#define WRITE_PROP(name, field) \
                    for(uint k = start[0], idx = 0; k < start[0] + dim[0]; k++) \
                    for(uint j = start[1]; j < start[1] + dim[1]; j++) \
                    for(uint i = start[2]; i < start[2] + dim[2]; i++, idx++) \
                    buf[idx] = grid.properties(i, j, k).field; \
                    hdf.createDataSet(name, datatype, space).write(buf.data(), datatype)

                    WRITE_PROP("young",   young);
                    WRITE_PROP("poisson", poisson);
                    WRITE_PROP("density", density);
                    WRITE_PROP("alpha",   stiffness);
                    WRITE_PROP("beta",    mass);
#undef WRITE_PROP
                }
            }

        }

        if (do_savepnt[mpi_rank]) {
            hsize_t count[] = {pi.size(), 3};

            {
                DataSpace space(2, count);
                DataSet   ds = hdf.createDataSet("pcoords", datatype, space);
                ds.write(pcoo.data(), datatype);
            }

            {
                DataSpace space(1, count);
                DataSet   ds = hdf.createDataSet("pids", PredType::NATIVE_INT32, space);
                ds.write(pids.data(), PredType::NATIVE_INT32);
            }
        }

        if (mpi_rank == 0) {
            xmf.append_child(pugi::node_doctype).set_value("Xdmf SYSTEM \"Sdmf.dtd\" []");

            auto xdmf = xmf.append_child("Xdmf");
            xdmf.append_attribute("xmlns:xi").set_value("http://www.w3.org/2003/XInclude");
            xdmf.append_attribute("Version" ).set_value("2.2");
            auto domain = xdmf.append_child("Domain");

            auto pgrid = domain.append_child("Grid");
            pgrid.append_attribute("Name"          ).set_value("Properties");
            pgrid.append_attribute("GridType"      ).set_value("Collection");
            pgrid.append_attribute("CollectionType").set_value("Spatial");

            for(int p = 0; p < mpi_size; p++) {
                if (!do_savebox[p] || nx[p] <= 1 || ny[p] <= 1 || nz[p] <= 1) continue;

                auto  prop = pgrid.append_child("Grid");
                prop.append_attribute("Name"    ).set_value((std::string("Prop_" + std::to_string(p))).c_str());
                prop.append_attribute("GridType").set_value("Uniform");

                auto topo = prop.append_child("Topology");
                topo.append_attribute("TopologyType"    ).set_value("3DRectMesh");
                topo.append_attribute("NumberOfElements").set_value( (
                            std::to_string(nz[p]) + " " +
                            std::to_string(ny[p]) + " " +
                            std::to_string(nx[p])
                            ).c_str());

                auto geom = prop.append_child("Geometry");
                geom.append_attribute("GeometryType").set_value("VXVYVZ");

                auto add_dimension = [&](size_t size, const std::string &name) {
                    auto data = geom.append_child("DataItem");
                    data.append_attribute("Dimensions").set_value(std::to_string(size).c_str());
                    data.append_attribute("NumberType").set_value("Float");
                    data.append_attribute("Precision" ).set_value(std::to_string(sizeof(real)).c_str());
                    data.append_attribute("Format"    ).set_value("HDF");

                    data.append_child(pugi::node_pcdata).set_value(
                            (hdfname[p] + ":/" + name).c_str()
                            );
                };

                add_dimension(nx[p], "x");
                add_dimension(ny[p], "y");
                add_dimension(nz[p], "z");

                auto add_attr = [&](const std::string &name) {
                    auto attr = prop.append_child("Attribute");
                    attr.append_attribute("Name"         ).set_value(name.c_str());
                    attr.append_attribute("AttributeType").set_value("Scalar");
                    attr.append_attribute("Center"       ).set_value("Cell");

                    auto data = attr.append_child("DataItem");
                    data.append_attribute("Dimensions").set_value( (
                                std::to_string(nz[p]) + " " +
                                std::to_string(ny[p]) + " " +
                                std::to_string(nx[p])
                                ).c_str());
                    data.append_attribute("NumberType").set_value("Float");
                    data.append_attribute("Precision" ).set_value(std::to_string(sizeof(real)).c_str());
                    data.append_attribute("Format"    ).set_value("HDF");

                    data.append_child(pugi::node_pcdata).set_value(
                            (hdfname[p] + ":/" + name).c_str()
                            );
                };

                add_attr("young");
                add_attr("poisson");
                add_attr("density");
                add_attr("alpha");
                add_attr("beta");
            }

            time_collection = domain.append_child("Grid");
            time_collection.append_attribute("Name"          ).set_value("Time");
            time_collection.append_attribute("GridType"      ).set_value("Collection");
            time_collection.append_attribute("CollectionType").set_value("Temporal");
        }

        if (do_savebox[mpi_rank] || do_savepnt[mpi_rank]) {
            hsize_t zero  = 0;
            hsize_t unlim = H5S_UNLIMITED;
            DataSpace space(1, &zero, &unlim);

            hsize_t chunk_size = 1024;
            H5::DSetCreatPropList p;
            p.setChunk(1, &chunk_size);

            time_ds = hdf.createDataSet("time", datatype, space, p);
        }

        {
            hsize_t n = 3 * part.size(mpi_rank);
            std::vector<real> zero(3 * n, 0);

            if (chkpt != 0) {
                hdf_chkpt = H5File(chkpt_name, H5F_ACC_TRUNC);
                hsize_t one = 1;

                hdf_chkpt.createGroup("checkpoint");
                chkpt_step = hdf_chkpt.createDataSet("checkpoint/step", PredType::NATIVE_INT32, DataSpace(1, &one));
                chkpt_time = hdf_chkpt.createDataSet("checkpoint/time", datatype, DataSpace(1, &one));
                chkpt_q0   = hdf_chkpt.createDataSet("checkpoint/q0",   datatype, DataSpace(1, &n));
                chkpt_q1   = hdf_chkpt.createDataSet("checkpoint/q1",   datatype, DataSpace(1, &n));

                save_checkpoint(0, zero, zero);

                have_chkpt = true;
            } else {
                have_chkpt = false;
            }

            add_step(0, zero, zero, zero);
        }

        if (do_savebox[mpi_rank] || do_savepnt[mpi_rank])
            hdf.flush(H5F_SCOPE_GLOBAL);
    }
}

//---------------------------------------------------------------------------
HDFSaver::~HDFSaver() {
    if (mpi_rank == 0)
	xmf.save_file(xmf_name.c_str(), "  ");
}

//---------------------------------------------------------------------------
void HDFSaver::add_step(real time,
        const std::vector<real> &q0,
        const std::vector<real> &q1,
        const std::vector<real> &d
        )
{
    step++;

    if (do_savebox[mpi_rank] || do_savepnt[mpi_rank]) {
        hsize_t cur_size;
        time_ds.getSpace().getSimpleExtentDims(&cur_size);

        if (step >= cur_size)
            time_ds.extend(&step);

        hsize_t start = step - 1;
        hsize_t count = 1;

        H5::DataSpace mspace = DataSpace(1, &count);
        H5::DataSpace fspace = time_ds.getSpace();

        fspace.selectHyperslab(H5S_SELECT_SET, &count, &start);

        time_ds.write(&time, datatype, mspace, fspace);
    }

    if (do_savebox[mpi_rank]) {
	std::ostringstream dsname;
	dsname << "Displ-" << step;

	DataSet dataset;

        if (H5Lexists(hdf.openGroup("/").getId(), dsname.str().c_str(), H5P_DEFAULT)) {
            dataset = hdf.openDataSet(dsname.str());
        } else {
            dataset = hdf.createDataSet(
                    dsname.str(), datatype, fbspace, prop);
        }

	dataset.write(d.data(), datatype, mbspace, fbspace);
    }

    if (do_savepnt[mpi_rank]) {
	std::ostringstream dsname;
	dsname << "points-" << step;

	DataSet dataset;

        if (H5Lexists(hdf.openGroup("/").getId(), dsname.str().c_str(), H5P_DEFAULT)) {
            dataset = hdf.openDataSet(dsname.str());
        } else {
            dataset = hdf.createDataSet(dsname.str(), datatype, fpspace);
        }

	dataset.write(d.data(), datatype, mpspace, fpspace);
    }

    if (mpi_rank == 0) {
        auto collection = time_collection.append_child("Grid");
        collection.append_attribute("Name"          ).set_value("Space");
        collection.append_attribute("GridType"      ).set_value("Collection");
        collection.append_attribute("CollectionType").set_value("Spatial");

        collection.append_child("Time").append_attribute("Value").set_value(std::to_string(time).c_str());

	for(int p = 0; p < mpi_size; p++) {
	    if (do_savebox[p]) {
                auto grid = collection.append_child("Grid");
                grid.append_attribute("Name"    ).set_value((std::string("Box_") + std::to_string(p)).c_str());
                grid.append_attribute("GridType").set_value("Uniform");

                auto topo = grid.append_child("Topology");
                topo.append_attribute("TopologyType"    ).set_value("3DRectMesh");
                topo.append_attribute("NumberOfElements").set_value( (
                            std::to_string(nz[p]) + " " +
                            std::to_string(ny[p]) + " " +
                            std::to_string(nx[p])
                            ).c_str());

                auto geom = grid.append_child("Geometry");
                geom.append_attribute("GeometryType").set_value("VXVYVZ");

                auto add_dimension = [&](size_t size, const std::string &name) {
                    auto data = geom.append_child("DataItem");
                    data.append_attribute("Dimensions").set_value(std::to_string(size).c_str());
                    data.append_attribute("NumberType").set_value("Float");
                    data.append_attribute("Precision" ).set_value(std::to_string(sizeof(real)).c_str());
                    data.append_attribute("Format"    ).set_value("HDF");

                    data.append_child(pugi::node_pcdata).set_value(
                            (hdfname[p] + ":/" + name).c_str()
                            );
                };

                add_dimension(nx[p], "x");
                add_dimension(ny[p], "y");
                add_dimension(nz[p], "z");

                auto attr = grid.append_child("Attribute");

                attr.append_attribute("Name"         ).set_value("Displ");
                attr.append_attribute("AttributeType").set_value("Vector");
                attr.append_attribute("Center"       ).set_value("Node");

                auto data = attr.append_child("DataItem");
                data.append_attribute("Dimensions").set_value( (
                            std::to_string(nz[p]) + " " +
                            std::to_string(ny[p]) + " " +
                            std::to_string(nx[p]) + " 3"
                            ).c_str());
                data.append_attribute("NumberType").set_value("Float");
                data.append_attribute("Precision" ).set_value(std::to_string(sizeof(real)).c_str());
                data.append_attribute("Format"    ).set_value("HDF");

                data.append_child(pugi::node_pcdata).set_value( (
                            hdfname[p] + ":/Displ-" + std::to_string(step)
                            ).c_str());
            }

	    if (do_savepnt[p]) {
                auto grid = collection.append_child("Grid");
                grid.append_attribute("Name"    ).set_value((std::string("Points_") + std::to_string(p)).c_str());
                grid.append_attribute("GridType").set_value("Uniform");

                auto topo = grid.append_child("Topology");
                topo.append_attribute("TopologyType"    ).set_value("Polyvertex");
                topo.append_attribute("NumberOfElements").set_value( std::to_string(np[p]).c_str() );

                auto geom = grid.append_child("Geometry");
                geom.append_attribute("GeometryType").set_value("XYZ");

                {
                    auto data = geom.append_child("DataItem");
                    data.append_attribute("Dimensions").set_value((std::to_string(np[p]) + " 3").c_str());
                    data.append_attribute("NumberType").set_value("Float");
                    data.append_attribute("Precision" ).set_value(std::to_string(sizeof(real)).c_str());
                    data.append_attribute("Format"    ).set_value("HDF");

                    data.append_child(pugi::node_pcdata).set_value(
                            (hdfname[p] + ":/pcoords").c_str()
                            );
                }

                {
                    auto attr = grid.append_child("Attribute");

                    attr.append_attribute("Name"         ).set_value("Displ");
                    attr.append_attribute("AttributeType").set_value("Vector");
                    attr.append_attribute("Center"       ).set_value("Node");

                    auto data = attr.append_child("DataItem");
                    data.append_attribute("Dimensions").set_value( (
                                std::to_string(np[p]) + " 3"
                                ).c_str());
                    data.append_attribute("NumberType").set_value("Float");
                    data.append_attribute("Precision" ).set_value(std::to_string(sizeof(real)).c_str());
                    data.append_attribute("Format"    ).set_value("HDF");

                    data.append_child(pugi::node_pcdata).set_value( (
                                hdfname[p] + ":/points-" + std::to_string(step)
                                ).c_str());
                }
            }
	}

        xmf.save_file(xmf_name.c_str(), "  ");
    }

    if (time_to_checkpoint()) {
        save_checkpoint(time, q0, q1);
    }

    if (do_savebox[mpi_rank] || do_savepnt[mpi_rank])
        hdf.flush(H5F_SCOPE_GLOBAL);
}

//---------------------------------------------------------------------------
bool HDFSaver::time_to_checkpoint() const {
    if (chkpt == 0) return false;

    return
        std::chrono::minutes(chkpt) <
        std::chrono::high_resolution_clock::now() - last_checkpoint;
}

//---------------------------------------------------------------------------
void HDFSaver::save_checkpoint(
        real t,
        const std::vector<real> &q0,
        const std::vector<real> &q1
        )
{
    int s = step;
    chkpt_step.write(&s, PredType::NATIVE_INT32);
    chkpt_time.write(&t, datatype);

    chkpt_q0.write(q0.data(), datatype);
    chkpt_q1.write(q1.data(), datatype);

    last_checkpoint = std::chrono::high_resolution_clock::now();

    hdf_chkpt.flush(H5F_SCOPE_GLOBAL);
}

//---------------------------------------------------------------------------
real HDFSaver::read_checkpoint(std::vector<real> &q0, std::vector<real> &q1) {
    if (have_chkpt) {
        chkpt_q0.read(q0.data(), datatype);
        chkpt_q1.read(q1.data(), datatype);

        int  s;
        real t;
        chkpt_step.read(&s, PredType::NATIVE_INT32);
        chkpt_time.read(&t, datatype);

        return t;
    } else {
        std::fill(q0.begin(), q0.end(), static_cast<real>(0));
        std::fill(q1.begin(), q1.end(), static_cast<real>(0));
        return 0;
    }
}

//---------------------------------------------------------------------------
void save_matrix(
        const VoigtModel::Matrices &local,
        const VoigtModel::Matrices &remote,
        const std::vector<int>     &recv_nbr,
        const std::vector<uint>    &recv_idx,
        const std::vector<int>     &send_nbr,
        const std::vector<uint>    &send_idx,
        const std::vector<uint>    &send_col
        )
{
    namespace fs = boost::filesystem;
    fs::create_directory(matrix_dir);

    H5File hdf(
            matrix_dir + "/process-" + std::to_string(mpi_rank) + ".h5",
            H5F_ACC_TRUNC
            );

#ifdef DOUBLE_PRECISION
      FloatType ftype(PredType::NATIVE_DOUBLE);
#else
      FloatType ftype(PredType::NATIVE_FLOAT);
#endif

      auto write_real = [&](const std::string &name, const std::vector<real> &v) {
          hsize_t n = v.size();
          DataSpace space(1, &n);
          hdf.createDataSet(name, ftype, space).write(v.data(), ftype);
      };

      auto write_uint = [&](const std::string &name, const std::vector<uint> &v) {
          hsize_t n = v.size();
          DataSpace space(1, &n);
          hdf.createDataSet(name, PredType::NATIVE_UINT32, space).write(
                  v.data(), PredType::NATIVE_UINT32);
      };

      auto write_int = [&](const std::string &name, const std::vector<int> &v) {
          hsize_t n = v.size();
          DataSpace space(1, &n);
          hdf.createDataSet(name, PredType::NATIVE_INT32, space).write(
                  v.data(), PredType::NATIVE_INT32);
      };

      hdf.createGroup("local");
      write_real("local/M",   local.M);
      write_uint("local/idx", local.idx);
      write_uint("local/row", local.row);
      write_int ("local/col", local.col);
      write_real("local/A1",  local.A1);
      write_real("local/A2",  local.A2);

      if (mpi_size > 1) {
          hdf.createGroup("remote");
          write_uint("remote/row", remote.row);
          write_int ("remote/col", remote.col);
          write_real("remote/A1",  remote.A1);
          write_real("remote/A2",  remote.A2);

          hdf.createGroup("recv");
          write_int ("recv/nbr", recv_nbr);
          write_uint("recv/idx", recv_idx);

          hdf.createGroup("send");
          write_int ("send/nbr", send_nbr);
          write_uint("send/idx", send_idx);
          write_uint("send/col", send_col);
      }
}

//---------------------------------------------------------------------------
void read_matrix(
        VoigtModel::Matrices &local,
        VoigtModel::Matrices &remote,
        std::vector<int>     &recv_nbr,
        std::vector<uint>    &recv_idx,
        std::vector<int>     &send_nbr,
        std::vector<uint>    &send_idx,
        std::vector<uint>    &send_col
        )
{
    H5File hdf(
            matrix_dir + "/process-" + std::to_string(mpi_rank) + ".h5",
            H5F_ACC_RDONLY
            );

#ifdef DOUBLE_PRECISION
      FloatType ftype(PredType::NATIVE_DOUBLE);
#else
      FloatType ftype(PredType::NATIVE_FLOAT);
#endif

    auto read_real = [&](const std::string &name, std::vector<real> &v) {
        DataSet ds = hdf.openDataSet(name);

        hsize_t n;
        ds.getSpace().getSimpleExtentDims(&n);

        v.resize(n);
	if (n) ds.read(v.data(), ftype);
    };

    auto read_uint = [&](const std::string &name, std::vector<uint> &v) {
        DataSet ds = hdf.openDataSet(name);

        hsize_t n;
        ds.getSpace().getSimpleExtentDims(&n);

        v.resize(n);
	if (n) ds.read(v.data(), PredType::NATIVE_UINT32);
    };

    auto read_int = [&](const std::string &name, std::vector<int> &v) {
        DataSet ds = hdf.openDataSet(name);

        hsize_t n;
        ds.getSpace().getSimpleExtentDims(&n);

        v.resize(n);
	if (n) ds.read(v.data(), PredType::NATIVE_INT32);
    };

      read_real("local/M",   local.M);
      read_uint("local/idx", local.idx);
      read_uint("local/row", local.row);
      read_int ("local/col", local.col);
      read_real("local/A1",  local.A1);
      read_real("local/A2",  local.A2);

      if (mpi_size > 1) {
          read_uint("remote/row", remote.row);
          read_int ("remote/col", remote.col);
          read_real("remote/A1",  remote.A1);
          read_real("remote/A2",  remote.A2);

          read_int ("recv/nbr", recv_nbr);
          read_uint("recv/idx", recv_idx);

          read_int ("send/nbr", send_nbr);
          read_uint("send/idx", send_idx);
          read_uint("send/col", send_col);
      } else {
          send_idx.push_back(0);
          recv_idx.push_back(0);
      }
}

