/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef BOND_CLASS

BondStyle(Fraclangevin, BondFRACLANGEVIN)

#else

#ifndef LMP_BOND_FRACLANGEVIN_H
#define LMP_BOND_FRACLANGEVIN_H

#include <stdio.h>
#include "bond.h"

namespace LAMMPS_NS {

	class BondFRACLANGEVIN : public Bond {
	public:
		BondFRACLANGEVIN(class LAMMPS*);
		virtual ~BondFRACLANGEVIN();
		virtual void compute(int, int);
		virtual void coeff(int, char**);
		void init_style();
		double equilibrium_distance(int);
		void write_restart(FILE*);
		void read_restart(FILE*);
		void write_data(FILE*);
		double single(int, double, int, int, double&);

	protected:
		double* b, * N, *critical_r_Nb;

		virtual void allocate();
	};

}

#endif
#endif

/* ERROR/WARNING messages:

W: Langevin bond too long: %ld %d %d %g

A Langevin bond has stretched dangerously far.  It's interaction strength
will be truncated to attempt to prevent the bond from blowing up.

E: Bad Langevin bond

Two atoms in a Langevin bond have become so far apart that the bond cannot
be computed.

E: Incorrect args for bond coefficients

Self-explanatory.  Check the input script or data file.

*/
