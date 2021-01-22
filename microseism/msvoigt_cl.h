#ifndef MSVOIGT_CL_H
#define MSVOIGT_CL_H

#include <vector>
#include <vexcl/vexcl.hpp>
#include "config.h"

class MSVoigtCL {
    public:
	MSVoigtCL(
		const std::vector<real>      &M,
		const std::vector<uint>      &idx,
		const std::vector<uint>      &loc_row,
		const std::vector<int>       &loc_col,
		const std::vector<real>	     &loc_A1,
		const std::vector<real>	     &loc_A2,

		const std::vector<uint>      &rem_row,
		const std::vector<int>       &rem_col,
		const std::vector<real>	     &rem_A1,
		const std::vector<real>	     &rem_A2,

		const std::vector<uint>      &send_col,
		uint                         remote_size,
                const std::vector<std::vector<real>> &q0,
                const std::vector<std::vector<real>> &q1
		);

	void get_q0(unsigned s, std::vector<real> &q) const;
	void get_q1(unsigned s, std::vector<real> &q) const;
	void set_q0(unsigned s, const std::vector<real> &q);

	void get_local_data (std::vector< std::vector<real> > &q);
	void set_remote_data(const std::vector< std::vector<real> > &q);

	void get_deriv(std::vector< std::vector<real> > &d, real tau);

	void step_local();
	void step_remote();

        std::string device() const;
    private:
	uint local_size, send_size, recv_size, ns;

        vex::Context            ctx;
        vex::vector<real>       M;
	vex::vector<uint>       idx;
	vex::vector<uint>       loc_row;
	vex::vector<int>        loc_col;
	vex::vector<real>       loc_A1;
	vex::vector<real>       loc_A2;
	vex::vector<uint>       rem_row;
	vex::vector<int>        rem_col;
	vex::vector<real>       rem_A1;
	vex::vector<real>       rem_A2;
	vex::vector<uint>       send_col;

        std::vector< vex::vector<real> > q0;
	std::vector< vex::vector<real> > q1;
	std::vector< vex::vector<real> > q2;
	std::vector< vex::vector<real> > q1_remote;
	std::vector< vex::vector<real> > q2_remote;
	std::vector< vex::vector<real> > send_val;

        vex::backend::kernel              adv_loc;
        vex::backend::kernel              adv_rem;
        vex::backend::kernel              gather;
        vex::backend::kernel              deriv;
};

#endif
