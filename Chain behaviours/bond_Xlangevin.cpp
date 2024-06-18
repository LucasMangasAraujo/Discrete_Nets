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

#include "bond_Xlangevin.h"

#include <math.h>
#include <stdlib.h>
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "update.h"
#include "force.h"
#include "memory.h"
#include "error.h"

// For debuging
#include <iostream>
using namespace std;

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondXLANGEVIN::BondXLANGEVIN(LAMMPS* lmp) : Bond(lmp) {}

/* ---------------------------------------------------------------------- */

BondXLANGEVIN::~BondXLANGEVIN()
{
    if (allocated && !copymode) {
        memory->destroy(setflag);
        memory->destroy(b);
        memory->destroy(N);
        memory->destroy(Eb); // Bonds Stiffness
        memory->destroy(RuptureEnergy); // Cleaveage Energy

    }
}

/* ---------------------------------------------------------------------- */

void BondXLANGEVIN::compute(int eflag, int vflag)
{
    int i1, i2, n, type, m;
    double delx, dely, delz, ebond, fbond;
    double r, rsq, Lc, beta, lambda;
    double bond_stretch, A, B, C, D, R, S, T, dummy, internal_energy, temp, entropic;

    double pi = M_PI;

    ebond = 0.0;
    ev_init(eflag, vflag);

    double** x = atom->x;
    double** f = atom->f;
    int** bondlist = neighbor->bondlist;
    int nbondlist = neighbor->nbondlist;
    int nlocal = atom->nlocal;
    int newton_bond = force->newton_bond;


    for (n = 0; n < nbondlist; n++) {

        // If the bond is broken skip#
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

        

        // Calculate the Bond Stretch
        A = 1.; B = -(1. + lambda); C = lambda; D = (-1. / Eb[type]) * lambda; // Polynomium coefficients

        R = ((3. * A * C) - pow(B, 2)) / (3. * pow(A, 2));
        S = (2. * pow(B, 3) - 9. * A * B * C + 27. * pow(A, 2) * D) / (27. * pow(A, 3)); // Auxiliar variables

        temp = (3. * S / (2. * R)) * sqrt(-3. / R);

        // Check for double precision errors
        if (temp >= 1) {
            temp = 0.;
            T = (1. / 3.) * temp;
        }
        else if (temp <= -1) {
            temp = pi;
            T = (1. / 3.) * temp;
        }
        else {
            //T = (1. / 3.) * acos((3. * S / (2. * R)) * sqrt(-3. / R));
            T = (1. / 3.) * acos(temp);
        }
    

        bond_stretch = (2. * sqrt(-R / 3.) * cos(T)) - (B / (3. * A)); // Final solution for the bond stretch
        internal_energy = (1. / 2.) * Eb[type] * pow(bond_stretch - 1., 2); // Internal Energy Contribution

        //Pade approximation of the inverse Langevin function accounting for the bond stretch
        dummy = lambda / bond_stretch;
        beta = dummy * (3. - pow(dummy, 2)) / (1. - pow(dummy, 2));
        entropic = ((lambda / bond_stretch) * beta + log(beta / sinh(beta)));

        // if bond breaks, set type to 0
        // both in temporary bondlist and permanent bond_type
        // if this proc owns both atoms,
        // negate bond_type twice if other atom stores it
        // if other proc owns 2nd atom, other proc will also break bond

        // Decide if bond has broken
        if (internal_energy >= RuptureEnergy[type] || isinf(entropic)) {
        //if (lambda >= 1 || isinf(entropic)) {
        //if (isinf(entropic) || lambda >= 0.95) {
            bondlist[n][2] = 0;
            cout << "--------------Bond " << n << "---------------------" << "\n";
            cout << "lambda_b is: " << bond_stretch << " entropic part is: " << entropic << "\n";
            cout << "The end-to-end distance is: " << r << " while the contour length is: "<< Lc << "\n";
            cout << "The ration r/Lc is: " << lambda << "\n";
            cout << "-----------------------------------------------" << "\n";
            for (m = 0; m < atom->num_bond[i1]; m++)
                if (atom->bond_atom[i1][m] == atom->tag[i2])
                    atom->bond_type[i1][m] = 0;
            if (i2 < atom->nlocal)
                for (m = 0; m < atom->num_bond[i2]; m++)
                    if (atom->bond_atom[i2][m] == atom->tag[i1])
                        atom->bond_type[i2][m] = 0;
            continue;
        }
        

        //if (r > 0.0) fbond = -beta / (r * b[type]);
        if (r > 0.0) fbond = ( -beta / bond_stretch) / (r * b[type]);
        else fbond = 0.0;
        
        // energy
        if (eflag) {

            //if (r > 0) ebond = N[type] * internal_energy + N[type] * ((lambda / bond_stretch) * beta + log(beta / sinh(beta)));
            if (r > 0) ebond = N[type] * internal_energy + N[type] * entropic;
            else ebond = 0., fbond = 0.;
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

        // subtract out pairwise contribution from 2 atoms via pair->single()
        // required since special_bond = 1,1,1
        // tally energy/virial in pair, using newton_bond as newton flag

        //itype = atom->type[i1];
        //jtype = atom->type[i2];




    } //end loop on bonds
}

/* ---------------------------------------------------------------------- */

void BondXLANGEVIN::allocate()
{
    allocated = 1;
    int n = atom->nbondtypes;

    memory->create(b, n + 1, "bond:b");
    memory->create(N, n + 1, "bond:N");
    memory->create(Eb, n + 1, "bond:Eb");
    memory->create(RuptureEnergy, n + 1, "bond:RuptureEnergy");
    memory->create(setflag, n + 1, "bond:setflag");
    for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void BondXLANGEVIN::coeff(int narg, char** arg)
{
    if (narg != 5) error->all(FLERR, "Incorrect args for bond coefficients");
    if (!allocated) allocate();

    int ilo, ihi;
    utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

    double b_one = utils::numeric(FLERR, arg[1], false, lmp);
    double N_one = utils::numeric(FLERR, arg[2], false, lmp);
    double Eb_one = utils::numeric(FLERR, arg[3], false, lmp);
    double RuptureEnergy_one = utils::numeric(FLERR, arg[4], false, lmp);
    
    int count = 0;
    for (int i = ilo; i <= ihi; i++) {
        b[i] = b_one;
        N[i] = N_one;
        Eb[i] = Eb_one;
        RuptureEnergy[i] = RuptureEnergy_one;
        setflag[i] = 1;
        count++;
    }

    if (count == 0) error->all(FLERR, "Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   check if special_bond settings are valid
------------------------------------------------------------------------- */

void BondXLANGEVIN::init_style()
{
    return;
}

/* ---------------------------------------------------------------------- */

double BondXLANGEVIN::equilibrium_distance(int i)
{
    return 0.0;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void BondXLANGEVIN::write_restart(FILE* fp)
{
    fwrite(&b[1], sizeof(double), atom->nbondtypes, fp);
    fwrite(&N[1], sizeof(double), atom->nbondtypes, fp);
    fwrite(&Eb[1], sizeof(double), atom->nbondtypes, fp);
    fwrite(&RuptureEnergy[1], sizeof(double), atom->nbondtypes, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void BondXLANGEVIN::read_restart(FILE* fp)
{
    allocate();

    if (comm->me == 0) {
        fread(&b[1], sizeof(double), atom->nbondtypes, fp);
        fread(&N[1], sizeof(double), atom->nbondtypes, fp);
        fread(&Eb[1], sizeof(double), atom->nbondtypes, fp);
        fread(&RuptureEnergy[1], sizeof(double), atom->nbondtypes, fp);
    }
    MPI_Bcast(&b[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
    MPI_Bcast(&N[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
    MPI_Bcast(&Eb[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
    MPI_Bcast(&RuptureEnergy[1], atom->nbondtypes, MPI_DOUBLE, 0, world);

    for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondXLANGEVIN::write_data(FILE* fp)
{
    for (int i = 1; i <= atom->nbondtypes; i++)
        fprintf(fp, "%d %g %g %g %g\n", i, b[i], N[i],Eb[i], RuptureEnergy[i]);
}

/* ---------------------------------------------------------------------- */

double BondXLANGEVIN::single(int type, double rsq, int i, int j,
    double& fforce)
{
    double r, lambda, Lc, beta, eng;
    double bond_stretch, A, B, C, D, R, S, T, dummy, internal_energy, temp;
    double pi = M_PI;
    r = sqrt(rsq);
    Lc = b[type] * N[type];
    lambda = r / Lc;

    if (type <= 0) return 0.0;

    // Calculate the Bond Stretch
    A = 1.; B = -(1. + lambda); C = lambda; D = (-1. / Eb[type]) * lambda; // Polynomium coefficients

    R = ((3. * A * C) - pow(B, 2)) / (3. * pow(A, 2));
    S = (2. * pow(B, 3) - 9. * A * B * C + 27. * pow(A, 2) * D) / (27. * pow(A, 3)); // Auxiliar variables

    temp = (3. * S / (2. * R)) * sqrt(-3. / R);

    // Check for double precision errors
    if (temp >= 1) {
        temp = 0.;
        T = (1. / 3.) * temp;
    }
    else if (temp <= -1) {
        temp = pi;
        T = (1. / 3.) * temp;
    }
    else {
        //T = (1. / 3.) * acos((3. * S / (2. * R)) * sqrt(-3. / R));
        T = (1. / 3.) * acos(temp);
    }


    bond_stretch = (2. * sqrt(-R / 3.) * cos(T)) - (B / (3. * A)); // Final solution for the bond stretch



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


    //beta = lambda * (3. - pow(lambda, 2)) / (1. - pow(lambda, 2));
    dummy = lambda / bond_stretch;
    beta = dummy * (3. - pow(dummy, 2)) / (1. - pow(dummy, 2));


    // Calculate the energetic contribution for the free-enregy
    internal_energy = (1. / 2.) * Eb[type] * pow(bond_stretch - 1., 2);

    if (r > 0 && internal_energy < RuptureEnergy[type]) {

        eng = N[type] * internal_energy + N[type] * ((lambda / bond_stretch) * beta + log(beta / sinh(beta)));

        //fforce = -beta / (r * b[type]);
        fforce = (-beta / bond_stretch) / (r * b[type]);

        
    }
    else {
        fforce = 0.0;
        eng = 0.0;
    }

    return eng;
}