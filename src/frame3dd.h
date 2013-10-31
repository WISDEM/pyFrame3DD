/*
 This file is part of FRAME3DD:
 Static and dynamic structural analysis of 2D and 3D frames and trusses with
 elastic and geometric stiffness.
 ---------------------------------------------------------------------------
 http://frame3dd.sourceforge.net/
 ---------------------------------------------------------------------------
 Copyright (C) 1992-2010  Henri P. Gavin
 
    FRAME3DD is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FRAME3DD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FRAME3DD.  If not, see <http://www.gnu.org/licenses/>.
*//**
	@file
	Main functions of the FRAME3DD solver API
*/

#ifndef FRAME_FRAME_H
#define FRAME_FRAME_H

/* for Micro-Stran compatability, structure for cartesian vectors */
#include "microstran/vec3.h"

/* maximum number of load cases */
#define _NL_ 32


/** form the global stiffness matrix */
void assemble_K(
	double **K,		/**< stiffness matrix			*/
	int DoF,		/**< number of degrees of freedom	*/
	int nE,			/**< number of frame elements		*/
	vec3 *xyz,		/**< XYZ locations of every node	*/
	float *r,		/**< rigid radius of every node	*/
	double *L, double *Le,	/**< length of each frame element, effective */
	int *N1, int *N2,	/**< node connectivity			*/
	float *Ax, float *Asy, float *Asz,	/**< section areas	*/
	float *Jx, float *Iy, float *Iz,	/**< section inertias	*/
	float *E, float *G,	/**< elastic and shear moduli		*/
	float *p,		/**< roll angle, radians		*/
	int shear,		/**< 1: include shear deformation, 0: don't */
	int geom,		/**< 1: include goemetric stiffness, 0: don't */
	double **Q,		/**< frame element end forces		*/
	int debug		/**< 1: write element stiffness matrices*/
);


/* compute_reaction_forces --- comput [K(r,q)] * {D(q)} + [K(r,r)] * {D(r)} */
 
void compute_reaction_forces( 
	double *F,	/**< vector of external loads and reaction forces  */
	double **K,	/**< stiffness matrix				*/
	double *D,	/**< displacement vector to be solved		*/
	int DoF,	/**< number of structural coordinates		*/
	int *r		/**< 0: not a reaction; 1: a reaction coordinate */
);


/** solve {F} =   [K]{D} via L D L' decomposition */
void solve_system(
	double **K,	/**< stiffness matrix for the restrained frame	*/
	double *D,	/**< displacement vector to be solved		*/
	double *F,	/**< load vector				*/
	int DoF,	/**< number of degrees of freedom		*/
	int *q,		/**< 1: not a reaction; 0: a reaction coordinate */
	int *r,		/**< 0: not a reaction; 1: a reaction coordinate */
	int *ok,	/**< indicates positive definite stiffness matrix */
	int verbose,	/**< 1: copious screen output; 0: none		*/
        double *rms_resid /**< the RMS error of the solution residual */
);


/** compute {dF} = {F} - [K]{D} and return ||dF|| / ||F||*/
double equilibrium_error(
        double *dF,	/**< equilibrium error  {dF} = {F} - [K]{D}	*/
        double *F,	/**< load vector                                */
        double **K,	/**< stiffness matrix for the restrained frame  */
        double *D,	/**< displacement vector to be solved           */
        int DoF,	/**< number of degrees of freedom               */
        int *q		/**< 1: not a reaction; 0: a reaction coordinate */
);


/** evaluate the member end forces for every member */
void element_end_forces(
	double **Q,	/**< frame element end forces			*/
	int nE,		/**< number of frame elements			*/
	vec3 *xyz,	/** XYZ locations of each node			*/
	double *L, double *Le,	/**< length of each frame element, effective */
	int *N1, int *N2,	/**< node connectivity			*/
	float *Ax, float *Asy, float *Asz,	/**< section areas	*/
	float *Jx, float *Iy, float *Iz,	/**< section area inertias */
	float *E, float *G,	/**< elastic and shear moduli		*/
	float *p,		/**< roll angle, radians		*/
	double *D,	/**< displacement vector			*/
	int shear,	/**< 1: include shear deformation, 0: don't */
	int geom	/**< 1: include goemetric stiffness, 0: don't */
);


/** add fixed end forces to internal element forces */
void add_feF(	
	vec3 *xyz,	/**< XYZ locations of each node		*/
	double *L,	/**< length of each frame element, effective	*/
	int *N1, int *N2, /**< node connectivity			*/
	float *p,	/**< roll angle, radians			*/
	double **Q,	/**< frame element end forces			*/
	double **feF_temp, /**< temp. fixed end forces for every frame element*/
	double **feF_mech, /**< mech. fixed end forces for every frame element*/
	int nE,		/**< number of frame elements			*/
	int DoF,	/**< number of degrees of freedom		*/
	double *F,	/**< vector of external loads and reaction forces  */
	int *r,		/**< 0: not a reaction; 1: a reaction coordinate */
	int verbose	/**< 1: copious screen output; 0: none		*/
);


/** assemble global mass matrix from element mass & inertia */
void assemble_M(
	double **M,	/**< mass matrix				*/
	int DoF,	/**< number of degrees of freedom		*/
	int nN, int nE,	/**< number of nodes, number of frame elements	*/
	vec3 *xyz,	/** XYZ locations of each node			*/
	float *r,	/**< rigid radius of every node		*/
	double *L,	/**< length of each frame element, effective	*/
	int *N1, int *N2, /**< node connectivity			*/
	float *Ax,	/**< node connectivity				*/
	float *Jx, float *Iy, float *Iz,	/**< section area inertias*/
	float *p,	/**< roll angle, radians			*/
	float *d,	/**< frame element density			*/
	float *EMs,	/**< extra frame element mass			*/
	float *NMs,	/**< node mass					*/
	float *NMx, float *NMy, float *NMz,	/**< node inertias	*/
	int lump,	/**< 1: lumped mass matrix, 0: consistent mass	*/
	int debug	/**< 1: write element mass matrices	 	*/
);


/** static condensation of stiffness matrix from NxN to nxn */
void condense(
	double **A,	/**< a square matrix				*/
	int N,		/**< the dimension of the matrix		*/
	int *q,		/**< list of matrix indices to retain		*/
	int n,		/**< the dimension of the condensed matrix	*/
	double **Ac,	/**< the condensed matrix			*/
	int verbose	/**< 1: copious screen output; 0: none		*/
);


/**
	generalized Guyan reduction of mass and stiffness matrices
	matches the response at a particular frequency, sqrt(L)/2/pi
	Guyan, Robert J., "Reduction of Stiffness and Mass Matrices",
	AIAA Journal, Vol. 3, No. 2 (1965) p 380.
*/
void guyan(
	double **M, double **K,	/**< mass and stiffness matrices	*/
	int N,			/**< dimension of the matrices, DoF	*/
	int *q,			/**< list of degrees of freedom to retain */
	int n,			/**< dimension of the condensed matrices */
	double **Mc, double **Kc,	/**< the condensed matrices	*/
	double w2,		/**< matched value of frequency squared	*/
	int verbose	/**< 1: copious screen output; 0: none		*/
);


/**
	dynamic condensation of mass and stiffness matrices
	matches the response at a set of frequencies

	@NOTE Kc and Mc may be ill-conditioned, and xyzsibly non-positive def.
*/
void dyn_conden(
	double **M, double **K,	/**< mass and stiffness matrices	*/
	int N,			/**< dimension of the matrices, DoF	*/
	int *R,		/**< R[i]=1: DoF i is fixed, R[i]=0: DoF i is free */
	int *p,		/**< list of primary degrees of freedom		*/
	int n,		/**< the dimension of the condensed matrix	*/
	double **Mc, double **Kc,	/**< the condensed matrices	*/
	double **V, double *f,	/**< mode shapes and natural frequencies*/
	int *m,		/**< list of modes to match in the condensed model */
	int verbose	/**< 1: copious screen output; 0: none		*/
);


/**
	release allocated memory
*/
void deallocate( 
	int nN, int nE, int nL, int *nF, int *nU, int *nW, int *nP, int *nT, int DoF,
	int modes,
	vec3 *xyz, float *rj, double *L, double *Le,
	int *N1, int *N2, int *q, int *r,
	float *Ax, float *Asy, float *Asz,
	float *Jx, float *Iy, float *Iz,
	float *E, float *G,
	float *p,
	float ***U, float ***W, float ***P, float ***T,
	float **Dp,
	double **F_mech, double **F_temp,
	double ***feF_mech, double ***feF_temp, double **F, double *dF, 
	double **K, double **Q,
	double *D, double *dD,
	float *d, float *EMs,
	float *NMs, float *NMx, float *NMy, float *NMz,
	double **M, double *f, double **V, 
	int *c, int *m
);


#endif /* FRAME_FRAME_H */

