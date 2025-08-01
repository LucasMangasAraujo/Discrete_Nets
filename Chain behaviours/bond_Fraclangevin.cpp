/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include "bond_Fraclangevin.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "update.h"
#include "force.h"
#include "memory.h"
#include "error.h"

#include <iostream>
using namespace std;

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondFRACLANGEVIN::BondFRACLANGEVIN(LAMMPS* lmp) : Bond(lmp) {}

/* ---------------------------------------------------------------------- */

BondFRACLANGEVIN::~BondFRACLANGEVIN()
{
    if (allocated && !copymode) {
        memory->destroy(setflag);
        memory->destroy(b);
        memory->destroy(N);
        memory->destroy(critical_r_Nb);
    }
}

/* ---------------------------------------------------------------------- */

void BondFRACLANGEVIN::compute(int eflag, int vflag)
{
    int i1, i2, n, type, m;
    double delx, dely, delz, ebond, fbond;
    double r, rsq, Lc, beta, lambda;

    ebond = 0.0;
    ev_init(eflag, vflag);

    double** x = atom->x;
    double** f = atom->f;
    int** bondlist = neighbor->bondlist;
    int nbondlist = neighbor->nbondlist;
    int nlocal = atom->nlocal;
    int newton_bond = force->newton_bond;


    for (n = 0; n < nbondlist; n++) {

        // Check if bond is broken
        if (bondlist[n][2] <= 0) continue;

        i1 = bondlist[n][0];
        i2 = bondlist[n][1];
        type = bondlist[n][2];

        delx = x[i1][0] - x[i2][0];
        dely = x[i1][1] - x[i2][1];
        delz = x[i1][2] - x[i2][2];

        Lc = N[type] * b[type]; //contour length

        rsq = delx * delx + dely * dely + delz * delz;
        r = sqrt(rsq);
        lambda = r / Lc;


        //Pade approximation of the inverse Langevin function
        beta = lambda * (3. - pow(lambda, 2)) / (1. - pow(lambda, 2));
        
        // Check if broken
        
        if (lambda >= 1) {
            cout << "--------------Bond " << n << "---------------------" << "\n";
            cout << "The ration r/Lc is: " << lambda << "\n";
            cout << "-----------------------------------------------" << "\n";
            char str[128];
            sprintf(str, "Langevin bond too long: " BIGINT_FORMAT " "
                TAGINT_FORMAT " " TAGINT_FORMAT " %g",
                update->ntimestep, atom->tag[i1], atom->tag[i2], sqrt(rsq));
            error->warning(FLERR, str, 0);
            lambda = 0.99999999;
        }
        

        if (r > 0.0) fbond = -beta / (r * b[type]);
        else fbond = 0.0;

        // energy
        if (eflag) {
            if (r > 0) ebond = N[type] * (lambda * beta + log(beta / sinh(beta)));
            else ebond = 0;
        }

        // apply force to each of 2 atoms
        if (newton_bond || i1 < nlocal) {
            f[i1][0] += delx * fbond;
            f[i1][1] += dely * fbond;
            f[i1][2] += delz * fbond;
        }

        if (newton_bond || i2 < nlocal) {
            f[i2][0] -= delx * fbond;
            f[i2][1] -= dely * fbond;
            f[i2][2] -= delz * fbond;
        }

        if (evflag) ev_tally(i1, i2, nlocal, newton_bond, ebond, fbond, delx, dely, delz);
    } //end loop on bonds
}

/* ---------------------------------------------------------------------- */

void BondFRACLANGEVIN::allocate()
{
    allocated = 1;
    int n = atom->nbondtypes;
    memory->create(b, n + 1, "bond:b");
    memory->create(N, n + 1, "bond:N");
    memory->create(critical_r_Nb, n + 1, "bond:critical_r_Nb");
    memory->create(setflag, n + 1, "bond:setflag");
    for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void BondFRACLANGEVIN::coeff(int narg, char** arg)
{
    if (narg != 4) error->all(FLERR, "Incorrect args for bond coefficients");
    if (!allocated) allocate();

    int ilo, ihi;
    utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

    double b_one = utils::numeric(FLERR, arg[1], false, lmp);
    double N_one = utils::numeric(FLERR, arg[2], false, lmp);
    double critical_r_Nb_one = utils::numeric(FLERR, arg[3], false, lmp);

    int count = 0;
    for (int i = ilo; i <= ihi; i++) {
        b[i] = b_one;
        N[i] = N_one;
        critical_r_Nb[i] = critical_r_Nb_one;
        setflag[i] = 1;
        count++;
    }

    if (count == 0) error->all(FLERR, "Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   check if special_bond settings are valid
------------------------------------------------------------------------- */

void BondFRACLANGEVIN::init_style()
{   
    return;
}

/* ---------------------------------------------------------------------- */

double BondFRACLANGEVIN::equilibrium_distance(int i)
{
    return 0.0;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void BondFRACLANGEVIN::write_restart(FILE* fp)
{
    fwrite(&b[1], sizeof(double), atom->nbondtypes, fp);
    fwrite(&N[1], sizeof(double), atom->nbondtypes, fp);
    fwrite(&critical_r_Nb[1], sizeof(double), atom->nbondtypes, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void BondFRACLANGEVIN::read_restart(FILE* fp)
{
    allocate();

    if (comm->me == 0) {
        fread(&b[1], sizeof(double), atom->nbondtypes, fp);
        fread(&N[1], sizeof(double), atom->nbondtypes, fp);
        fread(&critical_r_Nb[1], sizeof(double), atom->nbondtypes, fp);
    }
    MPI_Bcast(&b[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
    MPI_Bcast(&N[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
    MPI_Bcast(&critical_r_Nb[1], atom->nbondtypes, MPI_DOUBLE, 0, world);

    for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondFRACLANGEVIN::write_data(FILE* fp)
{
    for (int i = 1; i <= atom->nbondtypes; i++)
        fprintf(fp, "%d %g %g %g\n", i, b[i], N[i], critical_r_Nb[i]);
}

/* ---------------------------------------------------------------------- */

double BondFRACLANGEVIN::single(int type, double rsq, int i, int j,
    double& fforce)
{
    double r, lambda, Lc, beta, eng;
    r = sqrt(rsq);
    Lc = b[type] * N[type];
    lambda = r / Lc;

    if (type <= 0) return 0.0;
    

    //if (lambda >= 1){
    //  error->all(FLERR,"Bad Langevin bond");
    //}

    //if (lambda > 0.99999999) {
    //  char str[128];
    //  sprintf(str,"Langevin bond too long: " BIGINT_FORMAT " %g",
    //  update->ntimestep,r);
    //  error->warning(FLERR,str,0);
    //  lambda = 0.99999999;
    //}


    beta = lambda * (3. - pow(lambda, 2)) / (1. - pow(lambda, 2));

    if (r > 0) {
        fforce = -beta / (r * b[type]);
        eng = N[type] * (lambda * beta + log(beta / sinh(beta)));
    }
    else {
        fforce = 0.0;
        eng = 0.0;
    }

    return eng;
}
