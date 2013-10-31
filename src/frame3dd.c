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
*/

/* NOTE main 'driver' routine is now moved to main.c */

#include <math.h>
#include <assert.h>

#include "frame3dd.h"
#include "common.h"
#include "coordtrans.h"
#include "eig.h"
#include "HPGmatrix.h"
#include "NRutil.h"


/* #define MATRIX_DEBUG */

/* forward declarations */

static void elastic_K(
	double **k, vec3 *xyz, float *r,
	double L, double Le,
	int n1, int n2,
	float Ax, float Asy, float Asz,
	float Jx, float Iy, float Iz, float E, float G, float p,
	int shear
);

static void geometric_K(
	double **k, vec3 *xyz, float *r,
	double L, double Le, 
	int n1, int n2, float
	Ax, float Asy, float Asz, float Jx, float Iy, float Iz,
	float E, float G, float p, double T, 
	int shear
);

static void member_force(
	double *s, int M, vec3 *xyz, double L, double Le,
	int n1, int n2,
	float Ax, float Asy, float Asz, float Jx, float Iy, float Iz,
	float E, float G, float p, double *D,
	int shear, int geom, double *axial_strain
);

static void lumped_M(
	double **m, vec3 *xyz,
	double L, int n1, int n2,
	float Ax, float Jx, float Iy, float Iz, float p, 
	float d, float EMs
);

static void consistent_M(
	double **m, vec3 *xyz, float *r, double L,
	int n1, int n2,
	float Ax, float Jx, float Iy, float Iz, float p, 
	float d, float EMs
);


/*------------------------------------------------------------------------------
ASSEMBLE_K  -  assemble global stiffness matrix from individual elements 23feb94
------------------------------------------------------------------------------*/
void assemble_K(
	double **K,
	int DoF, int nE,
	vec3 *xyz, float *r, double *L, double *Le,
	int *N1, int *N2,
	float *Ax, float *Asy, float *Asz,
	float *Jx, float *Iy, float *Iz,
	float *E, float *G, float *p,
	int shear, int geom, double **Q, int debug
){
	double	**k;		/* element stiffness matrix in global coord */
	int	**ind,		/* member-structure DoF index table	*/
		res=0,
		i, j, ii, jj, l, ll;
	char	stiffness_fn[FILENMAX];

	for (i=1; i<=DoF; i++)	for (j=1; j<=DoF; j++)	K[i][j] = 0.0;

	k   =  dmatrix(1,12,1,12);
	ind = imatrix(1,12,1,nE);


	for ( i=1; i<= nE; i++ ) {
		ind[1][i] = 6*N1[i] - 5;	ind[7][i]  = 6*N2[i] - 5;
		ind[2][i] = ind[1][i] + 1;	ind[8][i]  = ind[7][i] + 1;
		ind[3][i] = ind[1][i] + 2;	ind[9][i]  = ind[7][i] + 2;
		ind[4][i] = ind[1][i] + 3;	ind[10][i] = ind[7][i] + 3;
		ind[5][i] = ind[1][i] + 4;	ind[11][i] = ind[7][i] + 4;
		ind[6][i] = ind[1][i] + 5;	ind[12][i] = ind[7][i] + 5;
	}

	for ( i = 1; i <= nE; i++ ) {

		elastic_K ( k, xyz, r, L[i], Le[i], N1[i], N2[i],
		Ax[i],Asy[i],Asz[i], Jx[i],Iy[i],Iz[i], E[i],G[i], p[i], shear);

		if (geom)
		 geometric_K( k, xyz, r, L[i], Le[i], N1[i], N2[i],
		           Ax[i], Asy[i],Asz[i], 
                           Jx[i], Iy[i], Iz[i], 
                           E[i],G[i], p[i], -Q[i][1], shear);

		if (debug) {
			res = sprintf(stiffness_fn,"k_%03d",i);
			save_dmatrix(stiffness_fn,k,1,12,1,12,0, "w");
		}

		for ( l=1; l <= 12; l++ ) {
			ii = ind[l][i];
			for ( ll=1; ll <= 12; ll++ ) {
				jj = ind[ll][i];
				K[ii][jj] += k[l][ll];
			}
		}
	}
	free_dmatrix ( k,1,12,1,12);
	free_imatrix(ind,1,12,1,nE);
	return;
}


/*------------------------------------------------------------------------------
ELASTIC_K - space frame elastic stiffness matrix in global coordnates	22oct02
------------------------------------------------------------------------------*/
void elastic_K(
	double **k, vec3 *xyz, float *r,
	double L, double Le,
	int n1, int n2,
	float Ax, float Asy, float Asz,
	float J, float Iy, float Iz, float E, float G, float p,
	int shear
){
	double   t1, t2, t3, t4, t5, t6, t7, t8, t9,     /* coord Xformn */
		Ksy, Ksz;		/* shear deformatn coefficients	*/
	int     i, j;

	coord_trans ( xyz, L, n1, n2,
				&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p );

	for (i=1;i<=12;i++)	for (j=1;j<=12;j++)	k[i][j] = 0.0;

	if ( shear ) {
		Ksy = 12.*E*Iz / (G*Asy*Le*Le);
		Ksz = 12.*E*Iy / (G*Asz*Le*Le);
	} else	Ksy = Ksz = 0.0;

	k[1][1]  = k[7][7]   = E*Ax / Le;
	k[2][2]  = k[8][8]   = 12.*E*Iz / ( Le*Le*Le*(1.+Ksy) );
	k[3][3]  = k[9][9]   = 12.*E*Iy / ( Le*Le*Le*(1.+Ksz) );
	k[4][4]  = k[10][10] = G*J / Le;
	k[5][5]  = k[11][11] = (4.+Ksz)*E*Iy / ( Le*(1.+Ksz) );
	k[6][6]  = k[12][12] = (4.+Ksy)*E*Iz / ( Le*(1.+Ksy) );

	k[5][3]  = k[3][5]   = -6.*E*Iy / ( Le*Le*(1.+Ksz) );
	k[6][2]  = k[2][6]   =  6.*E*Iz / ( Le*Le*(1.+Ksy) );
	k[7][1]  = k[1][7]   = -k[1][1];

	k[12][8] = k[8][12]  =  k[8][6] = k[6][8] = -k[6][2];
	k[11][9] = k[9][11]  =  k[9][5] = k[5][9] = -k[5][3];
	k[10][4] = k[4][10]  = -k[4][4];
	k[11][3] = k[3][11]  =  k[5][3];
	k[12][2] = k[2][12]  =  k[6][2];

	k[8][2]  = k[2][8]   = -k[2][2];
	k[9][3]  = k[3][9]   = -k[3][3];
	k[11][5] = k[5][11]  = (2.-Ksz)*E*Iy / ( Le*(1.+Ksz) );
	k[12][6] = k[6][12]  = (2.-Ksy)*E*Iz / ( Le*(1.+Ksy) );

#ifdef MATRIX_DEBUG
	save_dmatrix ( "ke", k, 1,12, 1,12, 0, "w" ); /* element elastic stiffness matrix */
#endif

	atma(t1,t2,t3,t4,t5,t6,t7,t8,t9, k, r[n1],r[n2]);	/* globalize */

	/* check and enforce symmetry of elastic element stiffness matrix */

	for(i=1; i<=12; i++){
	    for (j=i+1; j<=12; j++){
			if( k[i][j] != k[j][i] ) {
				if(fabs(k[i][j]/k[j][i]-1.0) > 1.0e-6 
					&&(fabs(k[i][j]/k[i][i]) > 1e-6 
						|| fabs(k[j][i]/k[i][i]) > 1e-6
					)
				){
					fprintf(stderr,"elastic_K: element stiffness matrix not symetric ...\n" ); 
					fprintf(stderr," ... k[%d][%d] = %15.6e \n",i,j,k[i][j] ); 
					fprintf(stderr," ... k[%d][%d] = %15.6e   ",j,i,k[j][i] ); 
					fprintf(stderr," ... relative error = %e \n",  fabs(k[i][j]/k[j][i]-1.0) ); 
					fprintf(stderr," ... element matrix saved in file 'kt'\n");
					save_dmatrix ( "kt", k, 1,12, 1,12, 0, "w" ); 
				}

				k[i][j] = k[j][i] = 0.5 * ( k[i][j] + k[j][i] );
			}
		}
	}
#ifdef MATRIX_DEBUG
	save_dmatrix ( "ket", k, 1,12, 1,12, 0, "w" ); /* transformed element matx */
#endif
}


/*------------------------------------------------------------------------------
GEOMETRIC_K - space frame geometric stiffness matrix, global coordnates 20dec07
------------------------------------------------------------------------------*/
void geometric_K(
	double **k, vec3 *xyz, float *r,
	double L, double Le,
	int n1, int n2, float Ax, float Asy, float Asz, float J,
	float Iy, float Iz, float E, float G, float p, double T,
	int shear
){
	double t1, t2, t3, t4, t5, t6, t7, t8, t9; /* coord Xformn */
	double **kg;
	double Ksy,Ksz,Dsy,Dsz; /* shear deformation coefficients	*/
	int i, j;

	coord_trans(
		xyz, L, n1, n2, &t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p 
	);

	kg = dmatrix(1,12,1,12);
	for (i=1;i<=12;i++)	for (j=1;j<=12;j++)	kg[i][j] = 0.0;

	if ( shear ) {
		Ksy = 12.*E*Iz / (G*Asy*Le*Le);
		Ksz = 12.*E*Iy / (G*Asz*Le*Le);
		Dsy = (1+Ksy)*(1+Ksy);
		Dsz = (1+Ksz)*(1+Ksz);
	} else{
		Ksy = Ksz = 0.0;
		Dsy = Dsz = 1.0;
	}


        kg[1][1]  = kg[7][7]   = 0.0; // T/L;

	kg[2][2]  = kg[8][8]   = T/L*(1.2+2.0*Ksy+Ksy*Ksy)/Dsy;
	kg[3][3]  = kg[9][9]   = T/L*(1.2+2.0*Ksz+Ksz*Ksz)/Dsz;
	kg[4][4]  = kg[10][10] = T/L*J/Ax;
	kg[5][5]  = kg[11][11] = T*L*(2.0/15.0+Ksz/6.0+Ksz*Ksz/12.0)/Dsz;
	kg[6][6]  = kg[12][12] = T*L*(2.0/15.0+Ksy/6.0+Ksy*Ksy/12.0)/Dsy;

        kg[1][7]  = kg[7][1]   = 0.0; // -T/L;
 
	kg[5][3]  = kg[3][5]   =  kg[11][3] = kg[3][11] = -T/10.0/Dsz;
	kg[9][5]  = kg[5][9]   =  kg[11][9] = kg[9][11] =  T/10.0/Dsz;
	kg[6][2]  = kg[2][6]   =  kg[12][2] = kg[2][12] =  T/10.0/Dsy;
	kg[8][6]  = kg[6][8]   =  kg[12][8] = kg[8][12] = -T/10.0/Dsy;

        kg[4][10] = kg[10][4]  = -kg[4][4];

	kg[8][2]  = kg[2][8]   = -T/L*(1.2+2.0*Ksy+Ksy*Ksy)/Dsy;
	kg[9][3]  = kg[3][9]   = -T/L*(1.2+2.0*Ksz+Ksz*Ksz)/Dsz;

	kg[11][5] = kg[5][11]  = -T*L*(1.0/30.0+Ksz/6.0+Ksz*Ksz/12.0)/Dsz;
	kg[12][6] = kg[6][12]  = -T*L*(1.0/30.0+Ksy/6.0+Ksy*Ksy/12.0)/Dsy;

#ifdef MATRIX_DEBUG
	save_dmatrix ( "kg", kg, 1,12, 1,12, 0, "w" ); /* element geom. stiffness matrix */
#endif

	atma(t1,t2,t3,t4,t5,t6,t7,t8,t9, kg, r[n1],r[n2]);	/* globalize */

	/* check and enforce symmetry of geometric element stiffness matrix */ 

	for(i=1; i<=12; i++){
	    for (j=i+1; j<=12; j++){
			if( kg[i][j] != kg[j][i] ) {
				if(fabs(kg[i][j]/kg[j][i]-1.0) > 1.0e-6 
					&&(fabs(kg[i][j]/kg[i][i]) > 1e-6 
						|| fabs(kg[j][i]/kg[i][i]) > 1e-6
					)
				){
					fprintf(stderr,"geometric_K element stiffness matrix not symetric ...\n" ); 
					fprintf(stderr," ... kg[%d][%d] = %15.6e \n",i,j,kg[i][j] ); 
					fprintf(stderr," ... kg[%d][%d] = %15.6e   ",j,i,kg[j][i] ); 
					fprintf(stderr," ... relative error = %e \n",  fabs(kg[i][j]/kg[j][i]-1.0) ); 
					fprintf(stderr," ... element matrix saved in file 'kg'\n");
					save_dmatrix ( "kg", kg, 1,12, 1,12, 0, "w" ); 
				}

				kg[i][j] = kg[j][i] = 0.5 * ( kg[i][j] + kg[j][i] );
			}
	    }
	}

#ifdef MATRIX_DEBUG
	save_dmatrix ( "kgt", kg, 1,12, 1,12, 0, "w" );   /* transformed element matx */
#endif

	/* add geometric stiffness matrix to elastic stiffness matrix ... */

	for (i=1; i<=12; i++)   for (j=1; j<=12; j++)	k[i][j] += kg[i][j];

	free_dmatrix(kg,1,12,1,12);
}


#if 0 /* DISUSED CODE */
/*------------------------------------------------------------------------------
END_RELEASE - apply matrix condensation for one member end force release 20nov04
------------------------------------------------------------------------------*/
void end_release ( X, r )
double	**X;
int	r;
{
	int	i,j;

	for (i=1; i<=12; i++)  if ( i != r) 
		for (j=1; j<=12; j++) if ( j != r) 
			X[i][j]  =  X[i][j]  -  X[i][r] * X[r][j] / X[r][r];
	
	for (i=1; i<=12; i++) 
		X[r][i] = X[i][r] = 0.0;

	return;
}
#endif


/*------------------------------------------------------------------------------
COMPUTE_REACTION_FORCES -  12oct12
Compute  [K(r,q)] * {D(q)} + [K(r,r)] * {D(r)}
The load vector modified for prescribed displacements Dp is returned as F
Prescribed displacements are "mechanican loads" not "temperature loads"  
------------------------------------------------------------------------------*/
void compute_reaction_forces( double *F, double **K, double *D, int DoF, int *r)
{
	int	i,j;

	for (i=1; i<=DoF; i++) {
		if (r[i]) {
			F[i] = 0.0;
			for (j=1; j<=DoF; j++)	F[i] += K[i][j]*D[j];
		}
	}

}


/*----------------------------------------------------------------------------
SOLVE_SYSTEM  -  solve {F} =   [K]{D} via L D L' decomposition        27dec01
Prescribed displacements are "mechanical loads" not "temperature loads"  
----------------------------------------------------------------------------*/
void solve_system(
	double **K, double *D, double *F, int DoF, int *q, int *r,
	int *ok, int verbose, double *rms_resid
){
	double	*diag;		/* diagonal vector of the L D L' decomp. */
	int	i;

	verbose = 0;		/* suppress verbose output		*/

	diag = dvector ( 1, DoF );

	/*  L D L' decomp */
	ldl_dcmp_pm ( K, DoF, diag, F, D, q,r, 1, 0, ok );
	if ( *ok < 0 ) {
	 	fprintf(stderr," Make sure that all six");
		fprintf(stderr," rigid body translations are restrained!\n");
		/* exit(31); */
	} else {				/* back substitute for D */
		/* LDL'  back-sub */
		ldl_dcmp_pm ( K, DoF, diag, F, D, q,r, 0, 1, ok );
		if ( verbose ) fprintf(stdout,"    LDL' RMS residual:");
		*rms_resid = *ok = 1;
		do {					/* improve solution */
			ldl_mprove_pm ( K, DoF, diag, F,D, q,r, rms_resid, ok );
			if ( verbose ) fprintf(stdout,"%9.2e", *rms_resid );
		} while ( *ok );
	        if ( verbose ) fprintf(stdout,"\n");
	}
	
	free_dvector( diag, 1, DoF );
}


/*----------------------------------------------------------------------------
EQUILIBRIUM_ERROR -  compute {dF} =   {F} - [K]{D}  and return ||dF||/||F||
----------------------------------------------------------------------------*/
double equilibrium_error( double *dF, double *F, double **K, double *D, int DoF, int *q )
{
	double	ss_dF = 0.0,	//  sum of squares of dF
		ss_F  = 0.0;	//  sum of squares of F	
	int	i,j;

	for (i=1; i<=DoF; i++) { // compute equilibrium error
		dF[i] = F[i];
		for (j=1; j<=DoF; j++) {
			if ( q[i] && K[i][j] != 0.0 && D[j] != 0.0 )
				dF[i] -= K[i][j]*D[j];
		}
	}       

	for (i=1; i<=DoF; i++) if (q[i]) ss_dF += ( dF[i] * dF[i] );
	for (i=1; i<=DoF; i++) if (q[i]) ss_F  += ( F[i]  * F[i] );

	return ( sqrt(ss_dF) / sqrt(ss_F) );	// convergence criterion
}


/*------------------------------------------------------------------------------
END_FORCES  -  evaluate the member end forces for every member		23feb94
------------------------------------------------------------------------------*/
void element_end_forces(
	double **Q, int nE, vec3 *xyz,
	double *L, double *Le,
	int *N1, int *N2,
	float *Ax, float *Asy, float *Asz,
	float *Jx, float *Iy, float *Iz, float *E, float *G, float *p,
	double *D, int shear, int geom
){
	double	*s, axial_strain;
	int	i,j;

	s = dvector(1,12);

	for(i=1; i <= nE; i++) {

     		member_force ( s, i, xyz, L[i], Le[i], N1[i], N2[i],
			Ax[i], Asy[i], Asz[i], Jx[i], Iy[i], Iz[i],
			E[i], G[i], p[i], D, shear, geom, &axial_strain );

		for(j=1; j<=12; j++)	Q[i][j] = s[j];

		if ( fabs(axial_strain > 0.001) )
		 fprintf(stderr," Warning! Frame element %2d has an average axial strain of %8.6f\n", i, axial_strain ); 

	}

	free_dvector(s,1,12);
}


/*------------------------------------------------------------------------------
MEMBER_FORCE  -  evaluate the end forces for a member			12nov02
------------------------------------------------------------------------------*/
void member_force(
	double *s, int M, vec3 *xyz, double L, double Le,
	int n1, int n2, float Ax, float Asy, float Asz, float J,
	float Iy, float Iz, float E, float G, float p, double *D,
	int shear, int geom, double *axial_strain
){
	double	t1, t2, t3, t4, t5, t6, t7, t8, t9, /* coord Xformn	*/
		d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12,
		x1, y1, z1, x2, y2, z2,	/* node coordinates	*/
		Ls,			/* stretched length of element */
		delta=0.0,		/* stretch in the frame element */
		Ksy, Ksz, Dsy, Dsz,	/* shear deformation coeff's	*/
		T = 0.0;		/* axial force for geometric stiffness */

	coord_trans ( xyz, L, n1, n2,
			&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p );

	x1 = xyz[n1].x;	y1 = xyz[n1].y;	z1 = xyz[n1].z;
	x2 = xyz[n2].x;	y2 = xyz[n2].y;	z2 = xyz[n2].z;

	n1 = 6*(n1-1);	n2 = 6*(n2-1);

	d1  = D[n1+1];	d2  = D[n1+2];	d3  = D[n1+3];
	d4  = D[n1+4];	d5  = D[n1+5];	d6  = D[n1+6];
	d7  = D[n2+1];	d8  = D[n2+2];	d9  = D[n2+3];
	d10 = D[n2+4];	d11 = D[n2+5];	d12 = D[n2+6];

	if ( shear ) {
		Ksy = 12.*E*Iz / (G*Asy*Le*Le);
		Ksz = 12.*E*Iy / (G*Asz*Le*Le);
		Dsy = (1+Ksy)*(1+Ksy);
		Dsz = (1+Ksz)*(1+Ksz);
	} else {
		Ksy = Ksz = 0.0;
		Dsy = Dsz = 1.0;
	}


	/* infinitessimal strain ... */
	delta = (d7-d1)*t1 + (d8-d2)*t2 + (d9-d3)*t3; 

	/* finite strain ... (not consistent with 2nd order formulation) */
/*  
 	delta += ( pow(((d7-d1)*t4 + (d8-d2)*t5 + (d9-d3)*t6),2.0) + 
 		   pow(((d7-d1)*t7 + (d8-d2)*t8 + (d9-d3)*t9),2.0) )/(2.0*L);
*/

	/* true strain ... (not appropriate for structural materials) */
/* 
  	Ls =	pow((x2+d7-x1-d1),2.0) + 
  		pow((y2+d8-y1-d2),2.0) + 
  		pow((z2+d9-z1-d3),2.0);
  	Ls = sqrt(Ls) + Le - L;

	delta = Le*log(Ls/Le);
*/

	*axial_strain = delta / Le;	// log(Ls/Le);

	if ( geom )	 T = Ax*E/Le * delta;
			 // T  = Ax*E*log(Ls/Le); 	/* true strain */

	if ( geom )
		s[1] = -T;
	else
		s[1]  =  -(Ax*E/Le)*( (d7-d1)*t1 + (d8-d2)*t2 + (d9-d3)*t3 );

	s[2]  = -( 12.*E*Iz/(Le*Le*Le*(1.+Ksy)) + 
		   T/L*(1.2+2.0*Ksy+Ksy*Ksy)/Dsy ) *
				( (d7-d1)*t4 + (d8-d2)*t5 + (d9-d3)*t6 )
		+ (6.*E*Iz/(Le*Le*(1.+Ksy)) + T/10.0/Dsy) *
				( (d4+d10)*t7 + (d5+d11)*t8 + (d6+d12)*t9 );
	s[3]  = -(12.*E*Iy/(Le*Le*Le*(1.+Ksz)) + 
		  T/L*(1.2+2.0*Ksz+Ksz*Ksz)/Dsz ) * 
				( (d7-d1)*t7  + (d8-d2)*t8  + (d9-d3)*t9 )
		- ( 6.*E*Iy/(Le*Le*(1.+Ksz)) + T/10.0/Dsz ) *
				( (d4+d10)*t4 + (d5+d11)*t5 + (d6+d12)*t6 );
	s[4]  =   -(G*J/Le) * ( (d10-d4)*t1 + (d11-d5)*t2 + (d12-d6)*t3 );
	s[5]  =   (6.*E*Iy/(Le*Le*(1.+Ksz)) + T/10.0/Dsz ) * 
				( (d7-d1)*t7 + (d8-d2)*t8 + (d9-d3)*t9 )
		+ ( (4.+Ksz)*E*Iy/(Le*(1.+Ksz)) + 
		    T*L*(2.0/15.0+Ksz/6.0+Ksz*Ksz/12.0)/Dsz ) *
				(d4 *t4 + d5 *t5 + d6 *t6)
		+ ((2.-Ksz)*E*Iy/(Le*(1.+Ksz)) -  
		    T*L*(1.0/30.0+Ksz/6.0+Ksz*Ksz/12.0)/Dsz ) *
				(d10*t4 + d11*t5 + d12*t6);
	s[6]  =  -( 6.*E*Iz/(Le*Le*(1.+Ksy)) + T/10.0/Dsy ) *
				( (d7-d1)*t4 + (d8-d2)*t5 + (d9-d3)*t6 )
		+ ((4.+Ksy)*E*Iz/(Le*(1.+Ksy)) + 
		    T*L*(2.0/15.0+Ksy/6.0+Ksy*Ksy/12.0)/Dsy ) *
				( d4 *t7 + d5 *t8 + d6 *t9 )
		+ ((2.-Ksy)*E*Iz/(Le*(1.+Ksy)) - 
		    T*L*(1.0/30.0+Ksy/6.0+Ksy*Ksy/12.0)/Dsy ) * 
				( d10*t7 + d11*t8 + d12*t9 );
	s[7]  = -s[1];
	s[8]  = -s[2]; 
	s[9]  = -s[3]; 
	s[10] = -s[4]; 

	s[11] =   ( 6.*E*Iy/(Le*Le*(1.+Ksz)) + T/10.0/Dsz ) * 
				( (d7-d1)*t7 + (d8-d2)*t8 + (d9-d3)*t9 )
		+ ((4.+Ksz)*E*Iy/(Le*(1.+Ksz)) + 
		    T*L*(2.0/15.0+Ksz/6.0+Ksz*Ksz/12.0)/Dsz ) * 
				( d10*t4 + d11*t5 + d12*t6 )
		+ ((2.-Ksz)*E*Iy/(Le*(1.+Ksz)) - 
		    T*L*(1.0/30.0+Ksz/6.0+Ksz*Ksz/12.0)/Dsz ) *
				( d4 *t4 + d5 *t5 + d6 *t6 );
	s[12] =  -(6.*E*Iz/(Le*Le*(1.+Ksy)) + T/10.0/Dsy ) *
				( (d7-d1)*t4 + (d8-d2)*t5 + (d9-d3)*t6 )
		+ ((4.+Ksy)*E*Iz/(Le*(1.+Ksy)) + 
		    T*L*(2.0/15.0+Ksy/6.0+Ksy*Ksy/12.0)/Dsy ) * 
				( d10*t7 + d11*t8 + d12*t9 )
		+ ((2.-Ksy)*E*Iz/(Le*(1.+Ksy)) - 
		    T*L*(1.0/30.0+Ksy/6.0+Ksy*Ksy/12.0)/Dsy ) * 
				( d4 *t7 + d5 *t8 + d6 *t9 );
}


/*----------------------------------------------------------------------------- 
ADD_FEF -  add fixed end forces to reactions and internal element forces 18oct12
------------------------------------------------------------------------------*/
void add_feF(	
	vec3 *xyz,
	double *L, int *N1, int *N2, float *p,
	double **Q, double **f_t, double **f_m, int nE, int DoF, 
	double *F, int *r,
	int verbose
){
	double  t1, t2, t3, t4, t5, t6, t7, t8, t9,	/* 3D coord Xformn */
		f1=0, f2=0, f3=0, f4=0,  f5=0,  f6=0, 
		f7=0, f8=0, f9=0, f10=0, f11=0, f12=0;
	int	m, n1, n2, i, i1, i2; //, J, x;

	for (m=1; m <= nE; m++) {	/* loop over all frame elements */

		n1 = N1[m];	n2 = N2[m];

		/* add fixed-end forces to reaction forces */
		for (i=1; i<=6; i++) {
			i1 = 6*(n1-1) + i;
			if (r[i1])
				F[i1] -= ( f_t[m][i] + f_m[m][i] );
		}
		for (i=1; i<=6; i++) {
			i2 = 6*(n2-1) + i;
			if (r[i2])
				F[i2] -= ( f_t[m][i+6] + f_m[m][i+6] );
		}

		coord_trans ( xyz, L[m], n1, n2,
			&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p[m] );

		// n1 = 6*(n1-1);	n2 = 6*(n2-1);	// ??

		// break out temperature fixed-end-forces to variables f1-f12
		f1  = f_t[m][1];   f2  = f_t[m][2];   f3  = f_t[m][3];
		f4  = f_t[m][4];   f5  = f_t[m][5];   f6  = f_t[m][6];
		f7  = f_t[m][7];   f8  = f_t[m][8];   f9  = f_t[m][9];
		f10 = f_t[m][10];  f11 = f_t[m][11];  f12 = f_t[m][12];

		// add mechanical load fixed-end-forces to variables f1-f12
		f1  += f_m[m][1];  f2  += f_m[m][2];  f3  += f_m[m][3];
		f4  += f_m[m][4];  f5  += f_m[m][5];  f6  += f_m[m][6];
		f7  += f_m[m][7];  f8  += f_m[m][8];  f9  += f_m[m][9];
		f10 += f_m[m][10]; f11 += f_m[m][11]; f12 += f_m[m][12];

		// add fixed end forces (-equivalent loads) to internal loads 
		// {Q} = [T]{F}
		Q[m][1]  -= ( f1 *t1 + f2 *t2 + f3 *t3 );    
		Q[m][2]  -= ( f1 *t4 + f2 *t5 + f3 *t6 );
		Q[m][3]  -= ( f1 *t7 + f2 *t8 + f3 *t9 );
		Q[m][4]  -= ( f4 *t1 + f5 *t2 + f6 *t3 );
		Q[m][5]  -= ( f4 *t4 + f5 *t5 + f6 *t6 );
		Q[m][6]  -= ( f4 *t7 + f5 *t8 + f6 *t9 );

		Q[m][7]  -= ( f7 *t1 + f8 *t2 + f9 *t3 );
		Q[m][8]  -= ( f7 *t4 + f8 *t5 + f9 *t6 );
		Q[m][9]  -= ( f7 *t7 + f8 *t8 + f9 *t9 );
		Q[m][10] -= ( f10*t1 + f11*t2 + f12*t3 );
		Q[m][11] -= ( f10*t4 + f11*t5 + f12*t6 );
		Q[m][12] -= ( f10*t7 + f11*t8 + f12*t9 );

	}

}


/*------------------------------------------------------------------------------
ASSEMBLE_M  -  assemble global mass matrix from element mass & inertia  24nov98
------------------------------------------------------------------------------*/
void assemble_M(
	double **M, int DoF, int nN, int nE,
	vec3 *xyz, float *r, double *L,
	int *N1, int *N2,
	float *Ax, float *Jx, float *Iy, float *Iz, float *p,
	float *d, float *EMs,
	float *NMs, float *NMx, float *NMy, float *NMz,
	int lump, int debug
){
	double  **m,	    /* element mass matrix in global coord */
		**dmatrix();
	int     **ind,	  /* member-structure DoF index table     */
		**imatrix(),
		res=0, 
		i, j, ii, jj, l, ll;
	char	mass_fn[FILENMAX];

	for (i=1; i<=DoF; i++)  for (j=1; j<=DoF; j++)  M[i][j] = 0.0;

	m      = dmatrix(1,12,1,12);
	ind    = imatrix(1,12,1,nE);


	for ( i=1; i<= nE; i++ ) {
		ind[1][i] = 6*N1[i] - 5;	ind[7][i]  = 6*N2[i] - 5;
		ind[2][i] = ind[1][i] + 1;      ind[8][i]  = ind[7][i] + 1;
		ind[3][i] = ind[1][i] + 2;      ind[9][i]  = ind[7][i] + 2;
		ind[4][i] = ind[1][i] + 3;      ind[10][i] = ind[7][i] + 3;
		ind[5][i] = ind[1][i] + 4;      ind[11][i] = ind[7][i] + 4;
		ind[6][i] = ind[1][i] + 5;      ind[12][i] = ind[7][i] + 5;
	}

	for ( i = 1; i <= nE; i++ ) {

		if ( lump )	lumped_M ( m, xyz, L[i], N1[i], N2[i],
				Ax[i], Jx[i], Iy[i], Iz[i], p[i], d[i], EMs[i]);
		else		consistent_M ( m, xyz,r,L[i], N1[i], N2[i],
				Ax[i], Jx[i], Iy[i], Iz[i], p[i], d[i], EMs[i]);

		if (debug) {
			res = sprintf(mass_fn,"m_%03d",i);
			save_dmatrix(mass_fn, m, 1,12, 1,12, 0, "w");
		}

		for ( l=1; l <= 12; l++ ) {
			ii = ind[l][i];
			for ( ll=1; ll <= 12; ll++ ) {
				jj = ind[ll][i];
				M[ii][jj] += m[l][ll];
			}
		}
	}

	for ( j = 1; j <= nN; j++ ) {		// add extra node mass
		i = 6*(j-1);
		M[i+1][i+1] += NMs[j];
		M[i+2][i+2] += NMs[j];
		M[i+3][i+3] += NMs[j];
		M[i+4][i+4] += NMx[j];
		M[i+5][i+5] += NMy[j];
		M[i+6][i+6] += NMz[j];
	}

	for (i=1; i<= DoF; i++) {
		if ( M[i][i] <= 0.0 ) {
			fprintf(stderr,"  error: Non pos-def mass matrix\n");
			fprintf(stderr,"  M[%d][%d] = %lf\n", i,i, M[i][i] );
		}
	}
	free_dmatrix ( m, 1,12,1,12);
	free_imatrix( ind,1,12,1,nE);
}


/*------------------------------------------------------------------------------
LUMPED_M  -  space frame element lumped mass matrix in global coordnates 7apr94
------------------------------------------------------------------------------*/
static void lumped_M(
	double **m, vec3 *xyz, double L,
	int n1, int n2, float Ax, float J, float Iy, float Iz, float p,
	float d, float EMs
){
	double   t1, t2, t3, t4, t5, t6, t7, t8, t9,     /* coord Xformn */
		t, ry,rz, po;	/* translational, rotational & polar inertia */
	int     i, j;

	coord_trans ( xyz, L, n1, n2,
				&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p );

			/* rotatory inertia of extra mass is neglected */

	t = ( d*Ax*L + EMs ) / 2.0;
	ry = d*Iy*L / 2.0;
	rz = d*Iz*L / 2.0;
	po = d*L*J / 2.0;		/* assumes simple cross-section	*/

	for (i=1;i<=12;i++)	for (j=1;j<=12;j++)	m[i][j] = 0.0;

	m[1][1] = m[2][2] = m[3][3] = m[7][7] = m[8][8] = m[9][9] = t;

	m[4][4] = m[10][10] = po*t1*t1 + ry*t4*t4 + rz*t7*t7;
	m[5][5] = m[11][11] = po*t2*t2 + ry*t5*t5 + rz*t8*t8;
	m[6][6] = m[12][12] = po*t3*t3 + ry*t6*t6 + rz*t9*t9;

	m[4][5] = m[5][4] = m[10][11] = m[11][10] =po*t1*t2 +ry*t4*t5 +rz*t7*t8;
	m[4][6] = m[6][4] = m[10][12] = m[12][10] =po*t1*t3 +ry*t4*t6 +rz*t7*t9;
	m[5][6] = m[6][5] = m[11][12] = m[12][11] =po*t2*t3 +ry*t5*t6 +rz*t8*t9;
}


/*------------------------------------------------------------------------------
CONSISTENT_M  -  space frame consistent mass matrix in global coordnates 2oct97
		 does not include shear deformations
------------------------------------------------------------------------------*/
void consistent_M(
	double **m, vec3 *xyz, float *r, double L,
	int n1, int n2, float Ax, float J, float Iy, float Iz, float p, 
	float d, float EMs
){
	double	t1, t2, t3, t4, t5, t6, t7, t8, t9,     /* coord Xformn */
		t, ry, rz, po;	/* translational, rotational & polar inertia */
	int     i, j;

	coord_trans ( xyz, L, n1, n2,
				&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p );

	t  =  d*Ax*L;	
	ry =  d*Iy;
	rz =  d*Iz;
	po =  d*J*L;

	for (i=1;i<=12;i++)	for (j=1;j<=12;j++)	m[i][j] = 0.0;

	m[1][1]  = m[7][7]   = t/3.;
	m[2][2]  = m[8][8]   = 13.*t/35. + 6.*rz/(5.*L);
	m[3][3]  = m[9][9]   = 13.*t/35. + 6.*ry/(5.*L);
	m[4][4]  = m[10][10] = po/3.;
	m[5][5]  = m[11][11] = t*L*L/105. + 2.*L*ry/15.;
	m[6][6]  = m[12][12] = t*L*L/105. + 2.*L*rz/15.;

	m[5][3]  = m[3][5]   = -11.*t*L/210. - ry/10.;
	m[6][2]  = m[2][6]   =  11.*t*L/210. + rz/10.;
	m[7][1]  = m[1][7]   =  t/6.;

	m[8][6]  = m[6][8]   =  13.*t*L/420. - rz/10.;
	m[9][5]  = m[5][9]   = -13.*t*L/420. + ry/10.;
	m[10][4] = m[4][10]  =  po/6.; 
	m[11][3] = m[3][11]  =  13.*t*L/420. - ry/10.;
	m[12][2] = m[2][12]  = -13.*t*L/420. + rz/10.;

	m[11][9] = m[9][11]  =  11.*t*L/210. + ry/10.;
	m[12][8] = m[8][12]  = -11.*t*L/210. - rz/10.;

	m[8][2]  = m[2][8]   =  9.*t/70. - 6.*rz/(5.*L);
	m[9][3]  = m[3][9]   =  9.*t/70. - 6.*ry/(5.*L);
	m[11][5] = m[5][11]  = -L*L*t/140. - ry*L/30.;
	m[12][6] = m[6][12]  = -L*L*t/140. - rz*L/30.;

			/* rotatory inertia of extra beam mass is neglected */

	for (i=1; i<=3; i++)	m[i][i] += 0.5*EMs;
	for (i=7; i<=9; i++)	m[i][i] += 0.5*EMs;

#ifdef MATRIX_DEBUG
	save_dmatrix ( "mo", m, 1,12, 1,12, 0, "w" ); /* element mass matrix */
#endif

	atma(t1,t2,t3,t4,t5,t6,t7,t8,t9, m, r[n1],r[n2]);	/* globalize */


	/* check and enforce symmetry of consistent element mass matrix */ 

	for(i=1; i<=12; i++){
	    for (j=i+1; j<=12; j++){
			if( m[i][j] != m[j][i] ) {
				if(fabs(m[i][j]/m[j][i]-1.0) > 1.0e-6 
					&&(fabs(m[i][j]/m[i][i]) > 1e-6 
						|| fabs(m[j][i]/m[i][i]) > 1e-6
					)
				){
					fprintf(stderr,"consistent_M: element mass matrix not symetric ...\n" ); 
					fprintf(stderr," ... m[%d][%d] = %15.6e \n",i,j,m[i][j] ); 
					fprintf(stderr," ... m[%d][%d] = %15.6e   ",j,i,m[j][i] ); 
					fprintf(stderr," ... relative error = %e \n",  fabs(m[i][j]/m[j][i]-1.0) ); 
					fprintf(stderr," ... element matrix saved in file 'mc'\n");
					save_dmatrix ( "mc", m, 1,12, 1,12, 0, "w" ); 
				}

				m[i][j] = m[j][i] = 0.5 * ( m[i][j] + m[j][i] );
			}
		}
	}

#ifdef MATRIX_DEBUG
	save_dmatrix ( "mt", m, 1,12, 1,12, 0, "w" );/* transformed matrix */
#endif
}



/*------------------------------------------------------------------------------
CONDENSE - static condensation of stiffness matrix from NxN to nxn    30aug01
------------------------------------------------------------------------------*/
void condense(
	double **A, int N, int *c, int n, double **Ac, int verbose
){
	double	**Arr, **Arc;
	int	i,j,k, ri,rj,ci,cj, ok, 
		*r;

	r    = ivector(1,N-n);
	Arr  = dmatrix(1,N-n,1,N-n);
	Arc  = dmatrix(1,N-n,1,n);

	k = 1;
	for (i=1; i<=N; i++) {
		ok = 1;
		for (j=1; j<=n; j++) {
			if ( c[j] == i ) {
				ok = 0;
				break;
			}
		}
		if ( ok )	r[k++] = i;
	}

	for (i=1; i<=N-n; i++) {
		for (j=i; j<=N-n; j++) { /* use only upper triangle of A */
			ri = r[i];
			rj = r[j];
			if ( ri <= rj )	Arr[j][i] = Arr[i][j] = A[ri][rj];
		}
	}

	for (i=1; i<=N-n; i++) {
		for (j=1; j<=n; j++) {	/* use only upper triangle of A */
			ri = r[i];
			cj = c[j];
			if ( ri < cj )	Arc[i][j] = A[ri][cj];
			else		Arc[i][j] = A[cj][ri];
		}
	}

	xtinvAy( Arc, Arr, Arc, N-n, n, Ac, verbose );

	for (i=1; i<=n; i++) {
		for (j=i; j<=n; j++) { /* use only upper triangle of A */
			ci = c[i];
			cj = c[j];
			if ( ci <= cj ) Ac[j][i]=Ac[i][j] = A[ci][cj]-Ac[i][j];
		}
	}

	free_ivector ( r,   1,N-n );
	free_dmatrix  ( Arr, 1,N-n,1,N-n );
	free_dmatrix  ( Arc, 1,N-n,1,n );
}


/*---------------------------------------------------------------------------- 
GUYAN  -   generalized Guyan reduction of mass and stiffness matrices    6jun07
           matches the response at a particular frequency, sqrt(L)/2/pi
           Guyan, Robert J., ``Reduction of Stiffness and Mass Matrices,''
           AIAA Journal, Vol. 3, No. 2 (1965) p 380.
-----------------------------------------------------------------------------*/
void guyan(
	double **M, double **K, int N,
	int *c, int n,
	double **Mc, double **Kc, double w2, 
	int verbose
){
	double	**Drr, **Drc, **invDrrDrc, **T;
	int	i,j,k, ri,rj,cj, ok, 
		*r;
	
	assert(M!=NULL);

	r   = ivector(1,N-n);
	Drr =  dmatrix(1,N-n,1,N-n);
	Drc =  dmatrix(1,N-n,1,n);
	invDrrDrc = dmatrix(1,N-n,1,n);	/* inv(Drr) * Drc	*/
	T   = dmatrix(1,N,1,n);	/* coordinate transformation matrix	*/

	w2 = 4.0 * PI * PI * w2 * w2;	/* eigen-value ... omega^2 	*/

	/* find "remaining" (r) degrees of freedom, not "condensed" (c)	*/
	k = 1;
	for (i=1; i<=N; i++) {
		ok = 1;
		for (j=1; j<=n; j++) {
			if ( c[j] == i ) {
				ok = 0;
				break;
			}
		}
		if ( ok )	r[k++] = i;
	}

	for (i=1; i<=N-n; i++) {
		for (j=1; j<=N-n; j++) { /* use only upper triangle of K,M */
			ri = r[i];
			rj = r[j];
			if ( ri <= rj )	
				Drr[j][i] = Drr[i][j] = K[ri][rj]-w2*M[ri][rj];
			else	Drr[j][i] = Drr[i][j] = K[rj][ri]-w2*M[rj][ri];
		}
	}

	for (i=1; i<=N-n; i++) {
		for (j=1; j<=n; j++) {	/* use only upper triangle of K,M */
			ri = r[i];
			cj = c[j];
			if ( ri < cj )	Drc[i][j] = K[ri][cj] - w2*M[ri][cj];
			else		Drc[i][j] = K[cj][ri] - w2*M[cj][ri];
		}
	}

	invAB(Drr, Drc, N-n, n, invDrrDrc, &ok, verbose); /* inv(Drr) * Drc */

	/* coordinate transformation matrix	*/	
	for (i=1; i<=n; i++) {
		for (j=1; j<=n; j++)	T[c[i]][j] =  0.0;
		T[c[i]][i] = 1.0;
	}	
	for (i=1; i<=N-n; i++) 
		for (j=1; j<=n; j++)	T[r[i]][j] = -invDrrDrc[i][j];

	xtAx ( K, T, Kc, N, n );		/* Kc = T' * K * T	*/

	xtAx ( M, T, Mc, N, n );		/* Mc = T' * M * T	*/

	free_ivector ( r,   1, N-n );
	free_dmatrix ( Drr, 1,N-n,1,N-n );
	free_dmatrix ( Drc, 1,N-n,1,n );
	free_dmatrix ( invDrrDrc, 1,N-n,1,N-n );
	free_dmatrix ( T, 1,N-n,1,n );
}


/*---------------------------------------------------------------------------- 
DYN_CONDEN - dynamic condensation of mass and stiffness matrices    8oct01
	     matches the response at a set of frequencies
WARNING: Kc and Mc may be ill-conditioned, and xyzsibly non-positive def.
-----------------------------------------------------------------------------*/
void dyn_conden(
	double **M, double **K, int N, int *R, int *p, int n,
	double **Mc, double **Kc, double **V, double *f, int *m,
	int verbose
){
	double	**P, **invP,
		traceM = 0, traceMc = 0, 
		Aij;		/* temporary storage for matrix mult. */
	int	i,j,k; 

	P    =  dmatrix(1,n,1,n);
	invP =  dmatrix(1,n,1,n);

	for (i=1; i<=n; i++)	/* first n modal vectors at primary DoF's */
		for (j=1; j<=n; j++)
			P[i][j] = V[p[i]][m[j]];

	pseudo_inv ( P, invP, n, n, 1e-9, verbose );

	for (i=1; i<=N; i++) if ( !R[i] ) traceM += M[i][i];

	for (i=1; i<=n; i++) { 		/* compute inv(P)' * I * inv(P)	*/
	    for (j=1; j<=n; j++) {
		Aij = 0.0;
	        for (k=1; k<=n; k++)
			Aij += invP[k][i] * invP[k][j];
		Mc[i][j] = Aij;
	    }
	}

	for (i=1; i<=n; i++) traceMc += Mc[i][i];

	for (i=1; i<=n; i++) { 		/* compute inv(P)' * W^2 * inv(P) */
	    for (j=1; j<=n; j++) {
		Aij = 0.0;
	        for (k=1; k<=n; k++) 
		    Aij += invP[k][i] * 4.0*PI*PI*f[m[k]]*f[m[k]] * invP[k][j];
		Kc[i][j] = Aij;
	    }
	}

	for (i=1; i<=n; i++)
	       for (j=1; j<=n; j++)
		       Mc[i][j] *= (traceM / traceMc);

	for (i=1; i<=n; i++)
	       for (j=1; j<=n; j++)
		       Kc[i][j] *= (traceM / traceMc);

	free_dmatrix  ( P,    1,n,1,n);
	free_dmatrix  ( invP, 1,n,1,n);
}


/*------------------------------------------------------------------------------
DEALLOCATE  -  release allocated memory					9sep08
------------------------------------------------------------------------------*/
void deallocate( 
	int nN, int nE, int nL, int *nF, int *nU, int *nW, int *nP, int *nT,
	int DoF, int nM,
	vec3 *xyz, float *rj, double *L, double *Le,
	int *N1, int *N2, int *q, int *r,
	float *Ax, float *Asy, float *Asz, float *J, float *Iy, float *Iz,
	float *E, float *G, float *p,
	float ***U, float ***W, float ***P, float ***T,
	float **Dp,
	double **F_mech, double **F_temp, 
	double ***feF_mech, double ***feF_temp, double **F, double *dF,
	double **K, double **Q,
	double *D, double *dD,
	float *d, float *EMs, float *NMs, float *NMx, float *NMy, float *NMz,
	double **M, double *f, double **V,
	int *c, int *m
){

	void	free();

	free(xyz);

	free_vector(rj,1,nN);
	free_dvector(L,1,nE);
	free_dvector(Le,1,nE);

// printf("..B\n"); /* debug */
	free_ivector(N1,1,nE);
	free_ivector(N2,1,nE);
	free_ivector(q,1,DoF);
	free_ivector(r,1,DoF);

// printf("..C\n"); /* debug */
	free_vector(Ax,1,nE);
	free_vector(Asy,1,nE);
	free_vector(Asz,1,nE);
	free_vector(J,1,nE);
	free_vector(Iy,1,nE);
	free_vector(Iz,1,nE);
	free_vector(E,1,nE);
	free_vector(G,1,nE);
	free_vector(p,1,nE);

// printf("..D\n"); /* debug */
	free_D3matrix(U,1,nL,1,nE,1,4);
	free_D3matrix(W,1,nL,1,10*nE,1,13);
	free_D3matrix(P,1,nL,1,10*nE,1,5);
	free_D3matrix(T,1,nL,1,nE,1,8);
	free_matrix(Dp,1,nL,1,DoF);

// printf("..E\n"); /* debug */
	free_dmatrix(F_mech,1,nL,1,DoF);
	free_dmatrix(F_temp,1,nL,1,DoF);

// printf("..F\n"); /* debug */
	free_D3dmatrix(feF_mech,1,nL,1,nE,1,12);
	free_D3dmatrix(feF_temp,1,nL,1,nE,1,12);

// printf("..G\n"); /* debug */
	free_dmatrix(F,1,nL,1,DoF);
	free_dvector(dF,1,DoF);

// printf("..H\n"); /* debug */
	free_dmatrix(K,1,DoF,1,DoF);
	free_dmatrix(Q,1,nE,1,12);

// printf("..I\n"); /* debug */
	free_dvector(D,1,DoF);
	free_dvector(dD,1,DoF);

// printf("..J\n"); /* debug */
	free_vector(d,1,nE);
	free_vector(EMs,1,nE);
	free_vector(NMs,1,nN);
	free_vector(NMx,1,nN);
	free_vector(NMy,1,nN);
	free_vector(NMz,1,nN);

// printf("..K\n"); /* debug */
	free_ivector(c,1,DoF);
	free_ivector(m,1,DoF);

// printf("..L\n"); /* debug */
	if ( nM > 0 ) {
		free_dmatrix(M,1,DoF,1,DoF);
		free_dvector(f,1,nM);
		free_dmatrix(V,1,DoF,1,DoF);
	}
}

/* itoa moved to frame3dd_io.c */

/* removed strcat -- it's in <string.h> in the standard C library */

/* removed strcpy -- it's in <string.h> in the standard C library */

/* dots moved to frame3dd_io.c */

