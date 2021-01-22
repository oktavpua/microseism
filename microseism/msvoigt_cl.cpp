#include <iostream>
#include <utility>
#include <sstream>
#include <vexcl/util.hpp>
#include "msvoigt_cl.h"
#include "util.h"

//---------------------------------------------------------------------------
MSVoigtCL::MSVoigtCL(
        const std::vector<real> &p_M,
        const std::vector<uint> &p_idx,
        const std::vector<uint> &p_loc_row,
        const std::vector<int>  &p_loc_col,
        const std::vector<real> &p_loc_A1,
        const std::vector<real> &p_loc_A2,

        const std::vector<uint> &p_rem_row,
        const std::vector<int>  &p_rem_col,
        const std::vector<real> &p_rem_A1,
        const std::vector<real> &p_rem_A2,

        const std::vector<uint> &p_send_col,
        uint                    p_recv_size,
        const std::vector<std::vector<real>> &p_q0,
        const std::vector<std::vector<real>> &p_q1
        )
    : local_size(p_idx.size()), send_size(p_send_col.size()),
      recv_size(p_recv_size), ns(p_q0.size()),
      ctx(vex::Filter::Exclusive(
                vex::Filter::DoublePrecision  &&
                vex::Filter::Env              &&
                vex::Filter::Count(1)
                )),
      M(ctx, p_M), idx(ctx, p_idx),
      loc_row(ctx, p_loc_row), loc_col(ctx, p_loc_col),
      loc_A1(ctx, p_loc_A1), loc_A2(ctx, p_loc_A2),
      q0(ns), q1(ns), q2(ns)
{
    precondition(!ctx.empty(), "No OpenCL devices found.");

    if (recv_size) {
        q1_remote.resize(ns);
        q2_remote.resize(ns);

        rem_row.resize(ctx, p_rem_row);
        rem_col.resize(ctx, p_rem_col);
        rem_A1.resize(ctx, p_rem_A1);
        rem_A2.resize(ctx, p_rem_A2);
    }

    if (send_size) {
        send_val.resize(ns);

        send_col.resize(ctx, p_send_col);
    }

    for(unsigned s = 0; s < ns; ++s) {
        q0[s].resize(ctx, p_q0[s]);
        q1[s].resize(ctx, p_q1[s]);
        q2[s].resize(ctx, local_size);

        q2[s] = 0;

        if (recv_size) {
            q1_remote[s].resize(ctx, p_recv_size);
            q2_remote[s].resize(ctx, p_recv_size);

            q2_remote[s] = 0;
        }

        if (send_size) {
            send_val[s].resize(ctx, send_size);
        }
    }


    std::string preamble =
#ifdef DOUBLE_PRECISION
        "typedef double real;\n"
#else
        "typedef float real;\n"
#endif
        ;

    using vex::global_ptr;

    {
        vex::backend::source_generator src(ctx.queue(0));
        src << preamble;

        src.begin_kernel("derivative");
        src.begin_kernel_parameters();
        src.parameter<cl_uint>("n");
        src.parameter<real   >("tau");

        for(unsigned s = 0; s < ns; ++s) {
            std::string suffix = "_" + std::to_string(s);

            src.parameter< global_ptr<const real> >("q1" + suffix);
            src.parameter< global_ptr<const real> >("q2" + suffix);
            src.parameter< global_ptr<      real> >("d"  + suffix);
        }
        src.end_kernel_parameters();

        src.grid_stride_loop("i").open("{");

        for(unsigned s = 0; s < ns; ++s)
            src.new_line() << "d_" << s << "[i] = "
                "(q1_" << s << "[i] - q2_" << s << "[i]) / tau;";

        src.close("}");
        src.end_kernel();

        deriv = vex::backend::kernel(ctx.queue(0), src.str(), "derivative");
    }

    {
        vex::backend::source_generator src(ctx.queue(0));
        src << preamble;

        src.begin_kernel("advance_local");
        src.begin_kernel_parameters();
        src.parameter< cl_uint                   >("n");
        src.parameter< global_ptr<const cl_uint> >("idx");
        src.parameter< global_ptr<const cl_uint> >("row");
        src.parameter< global_ptr<const cl_int > >("col");
        src.parameter< global_ptr<const real   > >("A1");
        src.parameter< global_ptr<const real   > >("A2");

        for(unsigned s = 0; s < ns; ++s) {
            std::string suffix = "_" + std::to_string(s);

            src.parameter< global_ptr<const real> >("q1" + suffix);
            src.parameter< global_ptr<const real> >("q2" + suffix);
            src.parameter< global_ptr<      real> >("q0" + suffix);
        }

        src.end_kernel_parameters();
        src.grid_stride_loop("i").open("{");
        src.new_line() << "uint ii        = idx[i];";
        src.new_line() << "uint row_start = row[ii];";
        src.new_line() << "uint row_end   = row[ii + 1];";

        for(unsigned s = 0; s < ns; ++s)
            src.new_line() << "real new_q_" << s << " = 0;";

        src.new_line() << "for(uint j = row_start; j < row_end; ++j)";
        src.open("{");

        src.new_line() << "int  c  = col[j] + i;";
        src.new_line() << "real a1 = A1[j];";
        src.new_line() << "real a2 = A2[j];";

        for(unsigned s = 0; s < ns; ++s) {
            src.new_line() << "new_q_" << s << " += "
                "a1 * q1_" << s << "[c] + a2 * q2_" << s << "[c];";
        }

        src.close("}");

        for(unsigned s = 0; s < ns; ++s)
            src.new_line() << "q0_" << s << "[i] = new_q_" << s << ";";

        src.close("}");
        src.end_kernel();

        adv_loc = vex::backend::kernel(ctx.queue(0), src.str(), "advance_local");
    }

    {
        vex::backend::source_generator src(ctx.queue(0));
        src << preamble;

        src.begin_kernel("advance_remote");
        src.begin_kernel_parameters();
        src.parameter< cl_uint                   >("n");
        src.parameter< global_ptr<const cl_uint> >("row");
        src.parameter< global_ptr<const cl_int > >("col");
        src.parameter< global_ptr<const real   > >("A1");
        src.parameter< global_ptr<const real   > >("A2");

        for(unsigned s = 0; s < ns; ++s) {
            std::string suffix = "_" + std::to_string(s);

            src.parameter< global_ptr<const real> >("q1" + suffix);
            src.parameter< global_ptr<const real> >("q2" + suffix);
            src.parameter< global_ptr<      real> >("q0" + suffix);
        }

        src.end_kernel_parameters();
        src.grid_stride_loop("i").open("{");

        src.new_line() << "uint row_start = row[i];";
        src.new_line() << "uint row_end   = row[i + 1];";

        for(unsigned s = 0; s < ns; ++s)
            src.new_line() << "real new_q_" << s << " = 0;";

        src.new_line() << "for(uint j = row_start; j < row_end; ++j)";
        src.open("{");

        src.new_line() << "int  c  = col[j];";
        src.new_line() << "real a1 = A1[j];";
        src.new_line() << "real a2 = A2[j];";

        for(unsigned s = 0; s < ns; ++s) {
            src.new_line() << "new_q_" << s << " += "
                "a1 * q1_" << s << "[c] + a2 * q2_" << s << "[c];";
        }

        src.close("}");

        for(unsigned s = 0; s < ns; ++s)
            src.new_line() << "q0_" << s << "[i] += new_q_" << s << ";";

        src.close("}");
        src.end_kernel();

        adv_rem = vex::backend::kernel(ctx.queue(0), src.str(), "advance_remote");
    }

    {
        vex::backend::source_generator src(ctx.queue(0));
        src << preamble;

        src.begin_kernel("gather");
        src.begin_kernel_parameters();
        src.parameter< cl_uint                   >("n");
        src.parameter< global_ptr<const cl_uint> >("send_col");

        for(unsigned s = 0; s < ns; ++s) {
            std::string suffix = "_" + std::to_string(s);

            src.parameter< global_ptr<const real> >("q" + suffix);
            src.parameter< global_ptr<      real> >("v" + suffix);
        }

        src.end_kernel_parameters();
        src.grid_stride_loop("i").open("{");

        src.new_line() << "uint c = send_col[i];";

        for(unsigned s = 0; s < ns; ++s)
            src.new_line() << "v_" << s << "[i] = q_" << s << "[c];";

        src.close("}");
        src.end_kernel();

        gather = vex::backend::kernel(ctx.queue(0), src.str(), "gather");
    }
}

//---------------------------------------------------------------------------
std::string MSVoigtCL::device() const {
    std::ostringstream buf;
    buf << ctx.queue(0);
    return buf.str();
}

//---------------------------------------------------------------------------
void MSVoigtCL::get_q0(unsigned s, std::vector<real> &q) const {
    vex::copy(q0[s], q);
}

//---------------------------------------------------------------------------
void MSVoigtCL::get_q1(unsigned s, std::vector<real> &q) const {
    vex::copy(q1[s], q);
}

//---------------------------------------------------------------------------
void MSVoigtCL::set_q0(unsigned s, const std::vector<real> &q) {
    vex::copy(q, q0[s]);
}

//---------------------------------------------------------------------------
void MSVoigtCL::get_local_data(std::vector< std::vector<real> > &q) {
    if (send_size) {
        gather.push_arg(send_size);
        gather.push_arg(send_col());

        for(unsigned s = 0; s < ns; ++s) {
            gather.push_arg(q0[s]());
            gather.push_arg(send_val[s]());
        }

        gather(ctx.queue(0));

        for(unsigned s = 0; s < ns; ++s)
            vex::copy(send_val[s], q[s]);
    }
}

//---------------------------------------------------------------------------
void MSVoigtCL::set_remote_data(const std::vector< std::vector<real> > &q) {
    for(unsigned s = 0; s < ns; ++s)
        if (q[s].size()) vex::copy(q[s], q1_remote[s]);
}

//---------------------------------------------------------------------------
void MSVoigtCL::get_deriv(std::vector< std::vector<real> > &d, real tau) {
    deriv.push_arg(local_size);
    deriv.push_arg(tau);

    for(unsigned s = 0; s < ns; ++s) {
        deriv.push_arg(q0[s]());
        deriv.push_arg(q1[s]());
        deriv.push_arg(q2[s]());
    }

    deriv(ctx.queue(0));

    for(unsigned s = 0; s < ns; ++s)
        vex::copy(q2[s], d[s]);
}

//---------------------------------------------------------------------------
void MSVoigtCL::step_local() {
    std::swap(q2, q1);
    std::swap(q1, q0);

    adv_loc.push_arg(local_size);
    adv_loc.push_arg(idx());
    adv_loc.push_arg(loc_row());
    adv_loc.push_arg(loc_col());
    adv_loc.push_arg(loc_A1());
    adv_loc.push_arg(loc_A2());

    for(unsigned s = 0; s < ns; ++s) {
        adv_loc.push_arg(q1[s]());
        adv_loc.push_arg(q2[s]());
        adv_loc.push_arg(q0[s]());
    }

    adv_loc(ctx.queue(0));
}

//---------------------------------------------------------------------------
void MSVoigtCL::step_remote() {
    if (recv_size) {
        adv_rem.push_arg(local_size);
        adv_rem.push_arg(rem_row());
        adv_rem.push_arg(rem_col());
        adv_rem.push_arg(rem_A1());
        adv_rem.push_arg(rem_A2());

        for(unsigned s = 0; s < ns; ++s) {
            adv_rem.push_arg(q1_remote[s]());
            adv_rem.push_arg(q2_remote[s]());
            adv_rem.push_arg(q0[s]());
        }

        adv_rem(ctx.queue(0));

        std::swap(q1_remote, q2_remote);
    }
}
