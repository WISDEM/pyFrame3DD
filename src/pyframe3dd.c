
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "frame3dd.h"
#include "frame3dd_io.h"
#include "eig.h"
#include "HPGmatrix.h"
#include "HPGutil.h"
#include "NRutil.h"



// ----------------
// Structs
// ----------------

typedef struct {

    int nN;
    int *N;
    double *x, *y, *z, *r;

} Nodes;


typedef struct {

    int nR;
    int* N;
    int *Rx, *Ry, *Rz, *Rxx, *Ryy, *Rzz;

} Reactions;


typedef struct {

    int nE;
    int *EL, *N1, *N2;
    double *Ax, *Asy, *Asz, *Jx, *Iy, *Iz, *E, *G, *roll, *density;

} Elements;


typedef struct {
    int nF;
    int *N;
    double *Fx, *Fy, *Fz, *Mxx, *Myy, *Mzz;

} PointLoads;


typedef struct {
    int nU;
    int *EL;
    double *Ux, *Uy, *Uz;

} UniformLoads;


typedef struct {
    int nD;
    int *N;
    double *Dx, *Dy, *Dz, *Dxx, *Dyy, *Dzz;

} PrescribedDisplacements;


typedef struct{

    double gx, gy, gz;
    PointLoads pointLoads;
    UniformLoads uniformLoads;
    PrescribedDisplacements prescribedDisplacements;

} LoadCase;


// return objects

typedef struct{

    int *node;
    double *x, *y, *z, *xrot, *yrot, *zrot;

} Displacements;


typedef struct{

    int *element;
    int *node;
    double *Nx, *Vy, *Vz, *Txx, *Myy, *Mzz;

} Forces;


typedef struct{

    int *node;
    double *Fx, *Fy, *Fz, *Mxx, *Myy, *Mzz;

} ReactionForces;


void readinputs(Nodes* nodes, Reactions* reactions, Elements* elements,
    int shear, int geom, double exagg_static, float dx,
    int nL, LoadCase* loadcases,
    Displacements* displacements, Forces* forces, ReactionForces* reactionForces){

    int i, j, l, n;
    char errMsg[MAXL];


    // ---- parse node data ---------

    int nN = nodes->nN;
    vec3    *xyz;       // X,Y,Z node coordinates (global)
    float *rj  =  vector(1, nN);        /* rigid radius around each node */
    xyz = (vec3 *)malloc(sizeof(vec3)*(1+nN));  /* node coordinates */

    for (i = 0; i < nN; i++) {
        j = nodes->N[i];
        if ( j <= 0 || j > nN ) {
            sprintf(errMsg,"\nERROR: in node coordinate data, node number out of range\n(node number %d is <= 0 or > %d)\n", j, nN);
            errorMsg(errMsg);
            exit(41);
        }
        xyz[j].x = nodes->x[i];
        xyz[j].y = nodes->y[i];
        xyz[j].z = nodes->z[i];
        rj[j] = fabs(nodes->r[i]);
    }

    int DoF = 6*nN;     /* total number of degrees of freedom   */

    // -----------------------------

    // ----- parse reaction data ---------

    int* q   = ivector(1, DoF);   /* allocate memory for reaction data ... */
    int* r   = ivector(1, DoF);   /* allocate memory for reaction data ... */

    for (i=1; i<=DoF; i++)  r[i] = 0;

    int nR = reactions->nR;

    if ( nR < 0 || nR > DoF/6 ) {
        fprintf(stderr," number of nodes with reactions ");
        dots(stderr,21);
        fprintf(stderr," nR = %3d ", nR );
        sprintf(errMsg,"\n  error: valid ranges for nR is 0 ... %d \n", DoF/6 );
        errorMsg(errMsg);
        exit(80);
    }


    int sumR;
    for (i = 0; i < nR; i++) {
        j = reactions->N[i];

        r[6*j-5] = reactions->Rx[i];
        r[6*j-4] = reactions->Ry[i];
        r[6*j-3] = reactions->Rz[i];
        r[6*j-2] = reactions->Rxx[i];
        r[6*j-1] = reactions->Ryy[i];
        r[6*j] = reactions->Rzz[i];

        if ( j > nN ) {
            sprintf(errMsg,"\n  error in reaction data: node number %d is greater than the number of nodes, %d \n", j, nN );
            errorMsg(errMsg);
            exit(81);
        }

        for (l=5; l >=0; l--) {
            if ( r[6*j-l] != 0 && r[6*j-l] != 1 ) {
                sprintf(errMsg,"\n  error in reaction data: Reaction data must be 0 or 1\n   Data for node %d, DoF %d is %d\n", j, 6-l, r[6*j-l] );
                errorMsg(errMsg);
                exit(82);
            }
        }

        sumR = 0;
        for (l=5; l >=0; l--)   sumR += r[6*j-l];
        if ( sumR == 0 ) {
            sprintf(errMsg,"\n  error: node %3d has no reactions\n   Remove node %3d from the list of reactions\n   and set nR to %3d \n", j, j, nR-1 );
            errorMsg(errMsg);
            exit(83);
        }
    }

    sumR=0;    for (i=1;i<=DoF;i++)    sumR += r[i];
    if ( sumR < 4 ) {
        sprintf(errMsg,"\n  Warning:  un-restrained structure   %d imposed reactions.\n  At least 4 reactions are required to support static loads.\n", sumR );
        errorMsg(errMsg);
        /*  exit(84); */
    }
    if ( sumR >= DoF ) {
        sprintf(errMsg,"\n  error in reaction data:  Fully restrained structure\n   %d imposed reactions >= %d degrees of freedom\n", sumR, DoF );
        errorMsg(errMsg);
        exit(85);
    }

    for (i=1; i<=DoF; i++)  if (r[i]) q[i] = 0; else q[i] = 1;

    // ---------------



    // ----- parse element data --------
    int nE = elements->nE;

    if (nN > nE + 1) { /* not enough elements */
        fprintf(stderr,"\n  warning: %d nodes and %d members...", nN, nE );
        fprintf(stderr," not enough elements to connect all nodes.\n");
    }

    /* allocate memory for frame elements ... */
    double *L   = dvector(1,nE);    /* length of each element       */
    double *Le  = dvector(1,nE);    /* effective length of each element */

    int *N1  = ivector(1,nE);    /* node #1 of each element      */
    int *N2  = ivector(1,nE);    /* node #2 of each element      */

    float *Ax  =  vector(1,nE);    /* cross section area of each element   */
    float *Asy =  vector(1,nE);    /* shear area in local y direction  */
    float *Asz =  vector(1,nE);    /* shear area in local z direction  */
    float *Jx  =  vector(1,nE);    /* torsional moment of inertia      */
    float *Iy  =  vector(1,nE);    /* bending moment of inertia about y-axis */
    float *Iz  =  vector(1,nE);    /* bending moment of inertia about z-axis */

    float *E   =  vector(1,nE);    /* frame element Young's modulus    */
    float *G   =  vector(1,nE);    /* frame element shear modulus      */
    float *p   =  vector(1,nE);    /* member rotation angle about local x axis */
    float *d   =  vector(1,nE);    /* member rotation angle about local x axis */



    int n1, n2;
    int epn0 = 0;  /* vector of elements per node */
    int *epn = ivector(1,nN);
    int b;

    for (n=1;n<=nN;n++) epn[n] = 0;

    for (i=0;i<nE;i++) {       /* read frame element properties */
        b = elements->EL[i];

        if ( b <= 0 || b > nE ) {
            sprintf(errMsg,"\n  error in frame element property data: Element number out of range  \n Frame element number: %d  \n", b);
            errorMsg(errMsg);
            exit(51);
        }

        N1[b] = elements->N1[i];
        N2[b] = elements->N2[i];

        epn[N1[b]] += 1;
        epn[N2[b]] += 1;

        if ( N1[b] <= 0 || N1[b] > nN || N2[b] <= 0 || N2[b] > nN ) {
            sprintf(errMsg,"\n  error in frame element property data: node number out of range  \n Frame element number: %d \n", b);
            errorMsg(errMsg);
            exit(52);
        }
        Ax[b] = elements->Ax[i];
        Asy[b] = elements->Asy[i];
        Asz[b] = elements->Asz[i];
        Jx[b] = elements->Jx[i];
        Iy[b] = elements->Iy[i];
        Iz[b] = elements->Iz[i];
        E[b] = elements->E[i];
        G[b] = elements->G[i];
        p[b] = elements->roll[i];
        d[b] = elements->density[i];

        p[b] = p[b]*PI/180.0;   /* convert from degrees to radians */

        if ( Ax[b] < 0 || Asy[b] < 0 || Asz[b] < 0 ||
             Jx[b] < 0 ||  Iy[b] < 0 ||  Iz[b] < 0  ) {
            sprintf(errMsg,"\n  error in frame element property data: section property < 0 \n  Frame element number: %d  \n", b);
            errorMsg(errMsg);
            exit(53);
        }
        if ( Ax[b] == 0 ) {
            sprintf(errMsg,"\n  error in frame element property data: cross section area is zero   \n  Frame element number: %d  \n", b);
            errorMsg(errMsg);
            exit(54);
        }
        if ( (Asy[b] == 0 || Asz[b] == 0) && G[b] == 0 ) {
            sprintf(errMsg,"\n  error in frame element property data: a shear area and shear modulus are zero   \n  Frame element number: %d  \n", b);
            errorMsg(errMsg);
            exit(55);
        }
        if ( Jx[b] == 0 ) {
            sprintf(errMsg,"\n  error in frame element property data: torsional moment of inertia is zero   \n  Frame element number: %d  \n", b);
            errorMsg(errMsg);
            exit(56);
        }
        if ( Iy[b] == 0 || Iz[b] == 0 ) {
            sprintf(errMsg,"\n  error: cross section bending moment of inertia is zero   \n  Frame element number : %d  \n", b);
            errorMsg(errMsg);
            exit(57);
        }
        if ( E[b] <= 0 || G[b] <= 0 ) {
            sprintf(errMsg,"\n  error : material elastic modulus E or G is not positive   \n  Frame element number: %d  \n", b);
            errorMsg(errMsg);
            exit(58);
        }
        if ( d[b] <= 0 ) {
            sprintf(errMsg,"\n  error : mass density d is not positive   \n  Frame element number: %d  \n", b);
            errorMsg(errMsg);
            exit(59);
        }
    }

    for (b=1;b<=nE;b++) {       /* calculate frame element lengths */
        n1 = N1[b];
        n2 = N2[b];

#define SQ(X) ((X)*(X))
        L[b] =  SQ( xyz[n2].x - xyz[n1].x ) +
            SQ( xyz[n2].y - xyz[n1].y ) +
            SQ( xyz[n2].z - xyz[n1].z );
#undef SQ

        L[b] = sqrt( L[b] );
        Le[b] = L[b] - r[n1] - r[n2];
        if ( n1 == n2 || L[b] == 0.0 ) {
            sprintf(errMsg,
              " Frame elements must start and stop at different nodes\n  frame element %d  N1= %d N2= %d L= %e\n   Perhaps frame element number %d has not been specified.\n  or perhaps the Input Data file is missing expected data.\n",
              b, n1,n2, L[b], i );
            errorMsg(errMsg);
            exit(60);
        }
        if ( Le[b] <= 0.0 ) {
            sprintf(errMsg, " Node  radii are too large.\n  frame element %d  N1= %d N2= %d L= %e \n  r1= %e r2= %e Le= %e \n",
              b, n1,n2, L[b], r[n1], r[n2], Le[b] );
            errorMsg(errMsg);
            exit(61);
        }
    }

    for ( n=1; n<=nN; n++ ) {
        if ( epn[n] == 0 ) {
            sprintf(errMsg,"node or frame element property data:\n     node number %3d is unconnected. \n", n);
            sferr(errMsg);
            epn0 += 1;
        }
    }

    free_ivector(epn,1,nN);

    if ( epn0 > 0 ) exit(42);

    // -------------------------------


    // -------- parse run data ---------------

    if (shear != 0 && shear != 1) {
        errorMsg(" Rember to specify shear deformations with a 0 or a 1 \n after the frame element property info.\n");
        exit(71);
    }

    if (geom != 0 && geom != 1) {
        errorMsg(" Rember to specify geometric stiffness with a 0 or a 1 \n after the frame element property info.\n");
        exit(72);
    }

    if ( exagg_static < 0.0 ) {
        errorMsg(" Remember to specify an exageration factor greater than zero.\n");
        exit(73);
    }

    if ( dx <= 0.0 && dx != -1 ) {
        errorMsg(" Remember to specify a frame element increment greater than zero.\n");
        exit(74);
    }

    //  -------------------------------


    // -------- parse load cases -------------

    if ( nL < 1 ) { /* not enough load cases */
        errorMsg("\n ERROR: the number of load cases must be at least 1\n");
        exit(101);
    }
    if ( nL >= _NL_ ) { /* too many load cases */
        sprintf(errMsg,"\n ERROR: maximum of %d load cases allowed\n", _NL_-1);
        errorMsg(errMsg);
        exit(102);
    }

    /* allocate memory for loads ... */
    float ***U   =  D3matrix(1,nL,1,nE,1,4);    /* uniform load on each member */
    float ***W   =  D3matrix(1,nL,1,10*nE,1,13);/* trapezoidal load on each member */
    float ***P   =  D3matrix(1,nL,1,10*nE,1,5); /* internal point load each member */
    float ***T   =  D3matrix(1,nL,1,nE,1,8);    /* internal temp change each member*/
    float **Dp  =  matrix(1,nL,1,DoF); /* prescribed displacement of each node */

    double **F_mech  = dmatrix(1,nL,1,DoF);  /* mechanical load vector   */
    double **F_temp  = dmatrix(1,nL,1,DoF);  /* temperature load vector  */
    double **F       = dmatrix(1,nL,1,DoF);  /* external load vectors    */
    double *dF  = dvector(1,DoF);   /* equilibrium error {F} - [K]{D} */

    double ***feF_mech =  D3dmatrix(1,nL,1,nE,1,12); /* feF due to mech loads */
    double ***feF_temp =  D3dmatrix(1,nL,1,nE,1,12); /* feF due to temp loads */

    double **K   = dmatrix(1,DoF,1,DoF); /* global stiffness matrix  */
    double **Q   = dmatrix(1,nE,1,12);   /* end forces for each member   */

    double *D   = dvector(1,DoF);   /* displacments of each node        */
    double *dD  = dvector(1,DoF);   /* incremental displ. of each node  */

    float *EMs =  vector(1,nE);    /* lumped mass for each frame element   */
    float *NMs =  vector(1,nN);    /* node mass for each node      */
    float *NMx =  vector(1,nN);    /* node inertia about global X axis */
    float *NMy =  vector(1,nN);    /* node inertia about global Y axis */
    float *NMz =  vector(1,nN);    /* node inertia about global Z axis */

    int *c = ivector(1,DoF);     /* vector of condensed degrees of freedom */
    int *m = ivector(1,DoF);     /* vector of condensed mode numbers */

    float gX[_NL_], gY[_NL_], gZ[_NL_];   // gravitational acceleration in global


    double **Fo = F;
    float   hy, hz;         /* section dimensions in local coords */
    float   x1,x2, w1,w2;
    double  Ln, R1o, R2o, f01, f02;
    double  Nx1, Vy1, Vz1, Mx1=0.0, My1=0.0, Mz1=0.0, /* fixed end forces */
        Nx2, Vy2, Vz2, Mx2=0.0, My2=0.0, Mz2=0.0,
        Ksy, Ksz,       /* shear deformatn coefficients */
        // a, b,           /* point load locations */
        t1, t2, t3, t4, t5, t6, t7, t8, t9; /* 3D coord Xfrm coeffs */
    int lc;

    for (j=1; j<=DoF; j++)
        for (lc=1; lc <= nL; lc++)
            Fo[lc][j] = F_temp[lc][j] = F_mech[lc][j] = 0.0;
    for (i=1; i<=12; i++)
        for (n=1; n<=nE; n++)
            for (lc=1; lc <= nL; lc++)
                feF_mech[lc][n][i] = feF_temp[lc][n][i] = 0.0;

    for (i=1; i<=DoF; i++)  for (lc=1; lc<=nL; lc++) Dp[lc][i] = 0.0;

    for (i=1;i<=nE;i++) for(j=1;j<=12;j++)  Q[i][j] = 0.0;

    LoadCase lcase;
    PointLoads pL;
    UniformLoads uL;
    PrescribedDisplacements pD;

    int *J1 = N1;
    int *J2 = N2;
    int nD[_NL_],   // number of prescribed nodal displ'nts
        nF[_NL_],   // number of loaded nodes
        nU[_NL_],   // number of members w/ unifm dist loads
        nW[_NL_],   // number of members w/ trapz dist loads
        nP[_NL_],   // number of members w/ conc point loads
        nT[_NL_],   // number of members w/ temp. changes
        nI=0,       // number of nodes w/ extra inertia
        nX=0,       // number of elemts w/ extra mass
        nC=0;       // number of condensed nodes

    for (lc = 1; lc <= nL; lc++) {      /* begin load-case loop */

        lcase = loadcases[lc-1];
        pL = lcase.pointLoads;
        uL = lcase.uniformLoads;
        pD = lcase.prescribedDisplacements;

        /* gravity loads applied uniformly to all frame elements  */
        gX[lc] = lcase.gx;
        gY[lc] = lcase.gy;
        gZ[lc] = lcase.gz;


        for (n=1; n<=nE; n++) {

            n1 = J1[n]; n2 = J2[n];

            coord_trans ( xyz, L[n], n1, n2,
                &t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p[n] );

            feF_mech[lc][n][1]  = d[n]*Ax[n]*L[n]*gX[lc] / 2.0;
            feF_mech[lc][n][2]  = d[n]*Ax[n]*L[n]*gY[lc] / 2.0;
            feF_mech[lc][n][3]  = d[n]*Ax[n]*L[n]*gZ[lc] / 2.0;

            feF_mech[lc][n][4]  = d[n]*Ax[n]*L[n]*L[n] / 12.0 *
                ( (-t4*t8+t5*t7)*gY[lc] + (-t4*t9+t6*t7)*gZ[lc] );
            feF_mech[lc][n][5]  = d[n]*Ax[n]*L[n]*L[n] / 12.0 *
                ( (-t5*t7+t4*t8)*gX[lc] + (-t5*t9+t6*t8)*gZ[lc] );
            feF_mech[lc][n][6]  = d[n]*Ax[n]*L[n]*L[n] / 12.0 *
                ( (-t6*t7+t4*t9)*gX[lc] + (-t6*t8+t5*t9)*gY[lc] );

            feF_mech[lc][n][7]  = d[n]*Ax[n]*L[n]*gX[lc] / 2.0;
            feF_mech[lc][n][8]  = d[n]*Ax[n]*L[n]*gY[lc] / 2.0;
            feF_mech[lc][n][9]  = d[n]*Ax[n]*L[n]*gZ[lc] / 2.0;

            feF_mech[lc][n][10] = d[n]*Ax[n]*L[n]*L[n] / 12.0 *
                ( ( t4*t8-t5*t7)*gY[lc] + ( t4*t9-t6*t7)*gZ[lc] );
            feF_mech[lc][n][11] = d[n]*Ax[n]*L[n]*L[n] / 12.0 *
                ( ( t5*t7-t4*t8)*gX[lc] + ( t5*t9-t6*t8)*gZ[lc] );
            feF_mech[lc][n][12] = d[n]*Ax[n]*L[n]*L[n] / 12.0 *
                ( ( t6*t7-t4*t9)*gX[lc] + ( t6*t8-t5*t9)*gY[lc] );

        }                 /* end gravity loads */


        nF[lc] = pL.nF;

        for (i=1; i <= nF[lc]; i++) { /* ! global structural coordinates ! */

            j = pL.N[i-1];

            if ( j < 1 || j > nN ) {
                sprintf(errMsg,"\n  error in node load data: node number out of range ... Node : %d\n   Perhaps you did not specify %d node loads \n  or perhaps the Input Data file is missing expected data.\n", j, nF[lc] );
                errorMsg(errMsg);
                exit(121);
            }

            F_mech[lc][6*j-5] = pL.Fx[i-1];
            F_mech[lc][6*j-4] = pL.Fy[i-1];
            F_mech[lc][6*j-3] = pL.Fz[i-1];
            F_mech[lc][6*j-2] = pL.Mxx[i-1];
            F_mech[lc][6*j-1] = pL.Myy[i-1];
            F_mech[lc][6*j] = pL.Mzz[i-1];

            if ( F_mech[lc][6*j-5]==0 && F_mech[lc][6*j-4]==0 && F_mech[lc][6*j-3]==0 && F_mech[lc][6*j-2]==0 && F_mech[lc][6*j-1]==0 && F_mech[lc][6*j]==0 )
                fprintf(stderr,"\n   Warning: All node loads applied at node %d  are zero\n", j );
        }         /* end node point loads  */


        nU[lc] = uL.nU;

        if ( nU[lc] < 0 || nU[lc] > nE ) {
            fprintf(stderr,"  number of uniformly distributed loads ");
            dots(stderr,13);
            fprintf(stderr," nU = %3d\n", nU[lc]);
            sprintf(errMsg,"\n  error: valid ranges for nU is 0 ... %d \n", nE );
            errorMsg(errMsg);
            exit(131);
        }
        for (i=1; i <= nU[lc]; i++) { /* ! local element coordinates ! */

            n = uL.EL[i-1];

            if ( n < 1 || n > nE ) {
                sprintf(errMsg,"\n  error in uniform distributed loads: element number %d is out of range\n",n);
                errorMsg(errMsg);
                exit(132);
            }
            U[lc][i][1] = (double) n;
            U[lc][i][2] = uL.Ux[i-1];
            U[lc][i][3] = uL.Uy[i-1];
            U[lc][i][4] = uL.Uz[i-1];

            if ( U[lc][i][2]==0 && U[lc][i][3]==0 && U[lc][i][4]==0 )
                fprintf(stderr,"\n   Warning: All distributed loads applied to frame element %d  are zero\n", n );

            Nx1 = Nx2 = U[lc][i][2]*Le[n] / 2.0;
            Vy1 = Vy2 = U[lc][i][3]*Le[n] / 2.0;
            Vz1 = Vz2 = U[lc][i][4]*Le[n] / 2.0;
            Mx1 = Mx2 = 0.0;
            My1 = -U[lc][i][4]*Le[n]*Le[n] / 12.0;  My2 = -My1;
            Mz1 =  U[lc][i][3]*Le[n]*Le[n] / 12.0;  Mz2 = -Mz1;

            n1 = J1[n]; n2 = J2[n];

            coord_trans ( xyz, L[n], n1, n2,
                &t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p[n] );

            /* {F} = [T]'{Q} */
            feF_mech[lc][n][1]  += ( Nx1*t1 + Vy1*t4 + Vz1*t7 );
            feF_mech[lc][n][2]  += ( Nx1*t2 + Vy1*t5 + Vz1*t8 );
            feF_mech[lc][n][3]  += ( Nx1*t3 + Vy1*t6 + Vz1*t9 );
            feF_mech[lc][n][4]  += ( Mx1*t1 + My1*t4 + Mz1*t7 );
            feF_mech[lc][n][5]  += ( Mx1*t2 + My1*t5 + Mz1*t8 );
            feF_mech[lc][n][6]  += ( Mx1*t3 + My1*t6 + Mz1*t9 );

            feF_mech[lc][n][7]  += ( Nx2*t1 + Vy2*t4 + Vz2*t7 );
            feF_mech[lc][n][8]  += ( Nx2*t2 + Vy2*t5 + Vz2*t8 );
            feF_mech[lc][n][9]  += ( Nx2*t3 + Vy2*t6 + Vz2*t9 );
            feF_mech[lc][n][10] += ( Mx2*t1 + My2*t4 + Mz2*t7 );
            feF_mech[lc][n][11] += ( Mx2*t2 + My2*t5 + Mz2*t8 );
            feF_mech[lc][n][12] += ( Mx2*t3 + My2*t6 + Mz2*t9 );


        }             /* end uniformly distributed loads */


      // sfrv=fscanf(fp,"%d", &nW[lc] ); /* trapezoidally distributed loads */
      // if (sfrv != 1) sferr("nW value in load data");
      // if ( verbose ) {
      //   fprintf(stdout,"  number of trapezoidally distributed loads ");
      //   dots(stdout,9); fprintf(stdout," nW = %3d\n", nW[lc]);
      // }
      // if ( nW[lc] < 0 || nW[lc] > 10*nE ) {
      //   sprintf(errMsg,"\n  error: valid ranges for nW is 0 ... %d \n", 10*nE );
      //   errorMsg(errMsg);
      //   exit(140);
      // }
      // for (i=1; i <= nW[lc]; i++) { /* ! local element coordinates ! */
      //   sfrv=fscanf(fp,"%d", &n );
      //   if (sfrv != 1) sferr("frame element number in trapezoidal load data");
      //   if ( n < 1 || n > nE ) {
      //       sprintf(errMsg,"\n  error in trapezoidally-distributed loads: element number %d is out of range\n",n);
      //       errorMsg(errMsg);
      //       exit(141);
      //   }
      //   W[lc][i][1] = (double) n;
      //   for (l=2; l<=13; l++) {
      //       sfrv=fscanf(fp,"%f", &W[lc][i][l] );
      //       if (sfrv != 1) sferr("value in trapezoidal load data");
      //   }

      //   Ln = L[n];

      //   /* error checking */

      //   if ( W[lc][i][ 4]==0 && W[lc][i][ 5]==0 &&
      //        W[lc][i][ 8]==0 && W[lc][i][ 9]==0 &&
      //        W[lc][i][12]==0 && W[lc][i][13]==0 ) {
      //     fprintf(stderr,"\n   Warning: All trapezoidal loads applied to frame element %d  are zero\n", n );
      //     fprintf(stderr,"     load case: %d , element %d , load %d\n ", lc, n, i );
      //   }

      //   if ( W[lc][i][ 2] < 0 ) {
      //     sprintf(errMsg,"\n   error in x-axis trapezoidal loads, load case: %d , element %d , load %d\n  starting location = %f < 0\n",
      //     lc, n, i , W[lc][i][2]);
      //     errorMsg(errMsg);
      //     exit(142);
      //   }
      //   if ( W[lc][i][ 2] > W[lc][i][3] ) {
      //     sprintf(errMsg,"\n   error in x-axis trapezoidal loads, load case: %d , element %d , load %d\n  starting location = %f > ending location = %f \n",
      //     lc, n, i , W[lc][i][2], W[lc][i][3] );
      //     errorMsg(errMsg);
      //     exit(143);
      //   }
      //   if ( W[lc][i][ 3] > Ln ) {
      //     sprintf(errMsg,"\n   error in x-axis trapezoidal loads, load case: %d , element %d , load %d\n ending location = %f > L (%f) \n",
      //     lc, n, i, W[lc][i][3], Ln );
      //     errorMsg(errMsg);
      //     exit(144);
      //   }
      //   if ( W[lc][i][ 6] < 0 ) {
      //     sprintf(errMsg,"\n   error in y-axis trapezoidal loads, load case: %d , element %d , load %d\n starting location = %f < 0\n",
      //     lc, n, i, W[lc][i][6]);
      //     errorMsg(errMsg);
      //     exit(142);
      //   }
      //   if ( W[lc][i][ 6] > W[lc][i][7] ) {
      //     sprintf(errMsg,"\n   error in y-axis trapezoidal loads, load case: %d , element %d , load %d\n starting location = %f > ending location = %f \n",
      //     lc, n, i, W[lc][i][6], W[lc][i][7] );
      //     errorMsg(errMsg);
      //     exit(143);
      //   }
      //   if ( W[lc][i][ 7] > Ln ) {
      //     sprintf(errMsg,"\n   error in y-axis trapezoidal loads, load case: %d , element %d , load %d\n ending location = %f > L (%f) \n",
      //     lc, n, i, W[lc][i][7],Ln );
      //     errorMsg(errMsg);
      //     exit(144);
      //   }
      //   if ( W[lc][i][10] < 0 ) {
      //     sprintf(errMsg,"\n   error in z-axis trapezoidal loads, load case: %d , element %d , load %d\n starting location = %f < 0\n",
      //     lc, n, i, W[lc][i][10]);
      //     errorMsg(errMsg);
      //     exit(142);
      //   }
      //   if ( W[lc][i][10] > W[lc][i][11] ) {
      //     sprintf(errMsg,"\n   error in z-axis trapezoidal loads, load case: %d , element %d , load %d\n starting location = %f > ending location = %f \n",
      //     lc, n, i, W[lc][i][10], W[lc][i][11] );
      //     errorMsg(errMsg);
      //     exit(143);
      //   }
      //   if ( W[lc][i][11] > Ln ) {
      //     sprintf(errMsg,"\n   error in z-axis trapezoidal loads, load case: %d , element %d , load %d\n ending location = %f > L (%f) \n",lc, n, i, W[lc][i][11], Ln );
      //     errorMsg(errMsg);
      //     exit(144);
      //   }

      //   if ( shear ) {
      //       Ksy = (12.0*E[n]*Iz[n]) / (G[n]*Asy[n]*Le[n]*Le[n]);
      //       Ksz = (12.0*E[n]*Iy[n]) / (G[n]*Asz[n]*Le[n]*Le[n]);
      //   } else  Ksy = Ksz = 0.0;

      //   /* x-axis trapezoidal loads (along the frame element length) */
      //   x1 =  W[lc][i][2]; x2 =  W[lc][i][3];
      //   w1 =  W[lc][i][4]; w2 =  W[lc][i][5];

      //   Nx1 = ( 3.0*(w1+w2)*Ln*(x2-x1) - (2.0*w2+w1)*x2*x2 + (w2-w1)*x2*x1 + (2.0*w1+w2)*x1*x1 ) / (6.0*Ln);
      //   Nx2 = ( -(2.0*w1+w2)*x1*x1 + (2.0*w2+w1)*x2*x2  - (w2-w1)*x1*x2 ) / ( 6.0*Ln );

      //   /* y-axis trapezoidal loads (across the frame element length) */
      //   x1 =  W[lc][i][6];  x2 = W[lc][i][7];
      //   w1 =  W[lc][i][8]; w2 =  W[lc][i][9];

      //   R1o = ( (2.0*w1+w2)*x1*x1 - (w1+2.0*w2)*x2*x2 +
      //        3.0*(w1+w2)*Ln*(x2-x1) - (w1-w2)*x1*x2 ) / (6.0*Ln);
      //   R2o = ( (w1+2.0*w2)*x2*x2 + (w1-w2)*x1*x2 -
      //       (2.0*w1+w2)*x1*x1 ) / (6.0*Ln);

      //   f01 = (  3.0*(w2+4.0*w1)*x1*x1*x1*x1 -  3.0*(w1+4.0*w2)*x2*x2*x2*x2
      //         - 15.0*(w2+3.0*w1)*Ln*x1*x1*x1 + 15.0*(w1+3.0*w2)*Ln*x2*x2*x2
      //         -  3.0*(w1-w2)*x1*x2*(x1*x1 + x2*x2)
      //         + 20.0*(w2+2.0*w1)*Ln*Ln*x1*x1 - 20.0*(w1+2.0*w2)*Ln*Ln*x2*x2
      //         + 15.0*(w1-w2)*Ln*x1*x2*(x1+x2)
      //         -  3.0*(w1-w2)*x1*x1*x2*x2 - 20.0*(w1-w2)*Ln*Ln*x1*x2 ) / 360.0;

      //   f02 = (  3.0*(w2+4.0*w1)*x1*x1*x1*x1 - 3.0*(w1+4.0*w2)*x2*x2*x2*x2
      //         -  3.0*(w1-w2)*x1*x2*(x1*x1+x2*x2)
      //         - 10.0*(w2+2.0*w1)*Ln*Ln*x1*x1 + 10.0*(w1+2.0*w2)*Ln*Ln*x2*x2
      //         -  3.0*(w1-w2)*x1*x1*x2*x2 + 10.0*(w1-w2)*Ln*Ln*x1*x2 ) / 360.0;

      //   Mz1 = -( 4.0*f01 + 2.0*f02 + Ksy*(f01 - f02) ) / ( Ln*Ln*(1.0+Ksy) );
      //   Mz2 = -( 2.0*f01 + 4.0*f02 - Ksy*(f01 - f02) ) / ( Ln*Ln*(1.0+Ksy) );

      //   Vy1 =  R1o + Mz1/Ln + Mz2/Ln;
      //   Vy2 =  R2o - Mz1/Ln - Mz2/Ln;

      //   /* z-axis trapezoidal loads (across the frame element length) */
      //   x1 =  W[lc][i][10]; x2 =  W[lc][i][11];
      //   w1 =  W[lc][i][12]; w2 =  W[lc][i][13];

      //   R1o = ( (2.0*w1+w2)*x1*x1 - (w1+2.0*w2)*x2*x2 +
      //        3.0*(w1+w2)*Ln*(x2-x1) - (w1-w2)*x1*x2 ) / (6.0*Ln);
      //   R2o = ( (w1+2.0*w2)*x2*x2 + (w1-w2)*x1*x2 -
      //       (2.0*w1+w2)*x1*x1 ) / (6.0*Ln);

      //   f01 = (  3.0*(w2+4.0*w1)*x1*x1*x1*x1 -  3.0*(w1+4.0*w2)*x2*x2*x2*x2
      //         - 15.0*(w2+3.0*w1)*Ln*x1*x1*x1 + 15.0*(w1+3.0*w2)*Ln*x2*x2*x2
      //         -  3.0*(w1-w2)*x1*x2*(x1*x1 + x2*x2)
      //         + 20.0*(w2+2.0*w1)*Ln*Ln*x1*x1 - 20.0*(w1+2.0*w2)*Ln*Ln*x2*x2
      //         + 15.0*(w1-w2)*Ln*x1*x2*(x1+x2)
      //         -  3.0*(w1-w2)*x1*x1*x2*x2 - 20.0*(w1-w2)*Ln*Ln*x1*x2 ) / 360.0;

      //   f02 = (  3.0*(w2+4.0*w1)*x1*x1*x1*x1 - 3.0*(w1+4.0*w2)*x2*x2*x2*x2
      //         -  3.0*(w1-w2)*x1*x2*(x1*x1+x2*x2)
      //         - 10.0*(w2+2.0*w1)*Ln*Ln*x1*x1 + 10.0*(w1+2.0*w2)*Ln*Ln*x2*x2
      //         -  3.0*(w1-w2)*x1*x1*x2*x2 + 10.0*(w1-w2)*Ln*Ln*x1*x2 ) / 360.0;

      //   My1 = ( 4.0*f01 + 2.0*f02 + Ksz*(f01 - f02) ) / ( Ln*Ln*(1.0+Ksz) );
      //   My2 = ( 2.0*f01 + 4.0*f02 - Ksz*(f01 - f02) ) / ( Ln*Ln*(1.0+Ksz) );

      //   Vz1 =  R1o - My1/Ln - My2/Ln;
      //   Vz2 =  R2o + My1/Ln + My2/Ln;

      //   /* debugging ... check internal force values
      //   printf("n=%d\n Nx1=%9.3f\n Nx2=%9.3f\n Vy1=%9.3f\n Vy2=%9.3f\n Vz1=%9.3f\n Vz2=%9.3f\n My1=%9.3f\n My2=%9.3f\n Mz1=%9.3f\n Mz2=%9.3f\n",
      //           n, Nx1,Nx2,Vy1,Vy2,Vz1,Vz2, My1,My2,Mz1,Mz2 );
      //   */

      //   n1 = J1[n]; n2 = J2[n];

      //   coord_trans ( xyz, Ln, n1, n2,
      //       &t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p[n] );

      //    debugging ... check coordinate transformation coefficients
      //   printf("t1=%5.2f t2=%5.2f t3=%5.2f \n", t1, t2, t3 );
      //   printf("t4=%5.2f t5=%5.2f t6=%5.2f \n", t4, t5, t6 );
      //   printf("t7=%5.2f t8=%5.2f t9=%5.2f \n", t7, t8, t9 );


      //   /* {F} = [T]'{Q} */
      //   feF_mech[lc][n][1]  += ( Nx1*t1 + Vy1*t4 + Vz1*t7 );
      //   feF_mech[lc][n][2]  += ( Nx1*t2 + Vy1*t5 + Vz1*t8 );
      //   feF_mech[lc][n][3]  += ( Nx1*t3 + Vy1*t6 + Vz1*t9 );
      //   feF_mech[lc][n][4]  += ( Mx1*t1 + My1*t4 + Mz1*t7 );
      //   feF_mech[lc][n][5]  += ( Mx1*t2 + My1*t5 + Mz1*t8 );
      //   feF_mech[lc][n][6]  += ( Mx1*t3 + My1*t6 + Mz1*t9 );

      //   feF_mech[lc][n][7]  += ( Nx2*t1 + Vy2*t4 + Vz2*t7 );
      //   feF_mech[lc][n][8]  += ( Nx2*t2 + Vy2*t5 + Vz2*t8 );
      //   feF_mech[lc][n][9]  += ( Nx2*t3 + Vy2*t6 + Vz2*t9 );
      //   feF_mech[lc][n][10] += ( Mx2*t1 + My2*t4 + Mz2*t7 );
      //   feF_mech[lc][n][11] += ( Mx2*t2 + My2*t5 + Mz2*t8 );
      //   feF_mech[lc][n][12] += ( Mx2*t3 + My2*t6 + Mz2*t9 );

      //   /* debugging ... check feF data
      //   for (l=1;l<=13;l++) printf(" %9.2e ", W[lc][i][l] );
      //   printf("\n");
      //   printf("n=%d ", n);
      //   for (l=1;l<=12;l++) {
      //       if (feF_mech[lc][n][l] != 0)
      //          printf(" feF %d = %9.3f ", l, feF_mech[lc][n][l] );
      //   }
      //   printf("\n");
      //   */
      // }         /* end trapezoidally distributed loads */

      // sfrv=fscanf(fp,"%d", &nP[lc] );   /* element point loads  */
      // if (sfrv != 1) sferr("nP value load data");
      // if ( verbose ) {
      //   fprintf(stdout,"  number of concentrated frame element point loads ");
      //   dots(stdout,2); fprintf(stdout," nP = %3d\n", nP[lc]);
      // }
      // if ( nP[lc] < 0 || nP[lc] > 10*nE ) {
      //   fprintf(stderr,"  number of concentrated frame element point loads ");
      //   dots(stderr,3);
      //   fprintf(stderr," nP = %3d\n", nP[lc]);
      //   sprintf(errMsg,"\n  error: valid ranges for nP is 0 ... %d \n", 10*nE );
      //   errorMsg(errMsg);
      //   exit(150);
      // }
      // for (i=1; i <= nP[lc]; i++) { /* ! local element coordinates ! */
      //   sfrv=fscanf(fp,"%d", &n );
      //   if (sfrv != 1) sferr("frame element number value point load data");
      //   if ( n < 1 || n > nE ) {
      //       sprintf(errMsg,"\n   error in internal point loads: frame element number %d is out of range\n",n);
      //       errorMsg(errMsg);
      //       exit(151);
      //   }
      //   P[lc][i][1] = (double) n;
      //   for (l=2; l<=5; l++) {
      //       sfrv=fscanf(fp,"%f", &P[lc][i][l] );
      //       if (sfrv != 1) sferr("value in point load data");
      //   }
      //   a = P[lc][i][5];    b = L[n] - a;

      //   if ( a < 0 || L[n] < a || b < 0 || L[n] < b ) {
      //       sprintf(errMsg,"\n  error in point load data: Point load coord. out of range\n   Frame element number: %d  L: %lf  load coord.: %lf\n",
      //       n, L[n], P[lc][i][5] );
      //       errorMsg(errMsg);
      //       exit(152);
      //   }

      //   if ( shear ) {
      //       Ksy = (12.0*E[n]*Iz[n]) / (G[n]*Asy[n]*Le[n]*Le[n]);
      //       Ksz = (12.0*E[n]*Iy[n]) / (G[n]*Asz[n]*Le[n]*Le[n]);
      //   } else  Ksy = Ksz = 0.0;

      //   Ln = L[n];

      //   Nx1 = P[lc][i][2]*a/Ln;
      //   Nx2 = P[lc][i][2]*b/Ln;

      //   Vy1 = (1./(1.+Ksz))    * P[lc][i][3]*b*b*(3.*a + b) / ( Ln*Ln*Ln ) +
      //       (Ksz/(1.+Ksz)) * P[lc][i][3]*b/Ln;
      //   Vy2 = (1./(1.+Ksz))    * P[lc][i][3]*a*a*(3.*b + a) / ( Ln*Ln*Ln ) +
      //       (Ksz/(1.+Ksz)) * P[lc][i][3]*a/Ln;

      //   Vz1 = (1./(1.+Ksy))    * P[lc][i][4]*b*b*(3.*a + b) / ( Ln*Ln*Ln ) +
      //       (Ksy/(1.+Ksy)) * P[lc][i][4]*b/Ln;
      //   Vz2 = (1./(1.+Ksy))    * P[lc][i][4]*a*a*(3.*b + a) / ( Ln*Ln*Ln ) +
      //       (Ksy/(1.+Ksy)) * P[lc][i][4]*a/Ln;

      //   Mx1 = Mx2 = 0.0;

      //   My1 = -(1./(1.+Ksy))  * P[lc][i][4]*a*b*b / ( Ln*Ln ) -
      //       (Ksy/(1.+Ksy))* P[lc][i][4]*a*b   / (2.*Ln);
      //   My2 =  (1./(1.+Ksy))  * P[lc][i][4]*a*a*b / ( Ln*Ln ) +
      //       (Ksy/(1.+Ksy))* P[lc][i][4]*a*b   / (2.*Ln);

      //   Mz1 =  (1./(1.+Ksz))  * P[lc][i][3]*a*b*b / ( Ln*Ln ) +
      //       (Ksz/(1.+Ksz))* P[lc][i][3]*a*b   / (2.*Ln);
      //   Mz2 = -(1./(1.+Ksz))  * P[lc][i][3]*a*a*b / ( Ln*Ln ) -
      //       (Ksz/(1.+Ksz))* P[lc][i][3]*a*b   / (2.*Ln);

      //   n1 = J1[n]; n2 = J2[n];

      //   coord_trans ( xyz, Ln, n1, n2,
      //       &t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p[n] );

      //   /* {F} = [T]'{Q} */
      //   feF_mech[lc][n][1]  += ( Nx1*t1 + Vy1*t4 + Vz1*t7 );
      //   feF_mech[lc][n][2]  += ( Nx1*t2 + Vy1*t5 + Vz1*t8 );
      //   feF_mech[lc][n][3]  += ( Nx1*t3 + Vy1*t6 + Vz1*t9 );
      //   feF_mech[lc][n][4]  += ( Mx1*t1 + My1*t4 + Mz1*t7 );
      //   feF_mech[lc][n][5]  += ( Mx1*t2 + My1*t5 + Mz1*t8 );
      //   feF_mech[lc][n][6]  += ( Mx1*t3 + My1*t6 + Mz1*t9 );

      //   feF_mech[lc][n][7]  += ( Nx2*t1 + Vy2*t4 + Vz2*t7 );
      //   feF_mech[lc][n][8]  += ( Nx2*t2 + Vy2*t5 + Vz2*t8 );
      //   feF_mech[lc][n][9]  += ( Nx2*t3 + Vy2*t6 + Vz2*t9 );
      //   feF_mech[lc][n][10] += ( Mx2*t1 + My2*t4 + Mz2*t7 );
      //   feF_mech[lc][n][11] += ( Mx2*t2 + My2*t5 + Mz2*t8 );
      //   feF_mech[lc][n][12] += ( Mx2*t3 + My2*t6 + Mz2*t9 );
      // }                 /* end element point loads */

      // sfrv=fscanf(fp,"%d", &nT[lc] );   /* thermal loads    */
      // if (sfrv != 1) sferr("nT value in load data");
      // if ( verbose ) {
      //   fprintf(stdout,"  number of temperature changes ");
      //   dots(stdout,21); fprintf(stdout," nT = %3d\n", nT[lc] );
      // }
      // if ( nT[lc] < 0 || nT[lc] > nE ) {
      //   fprintf(stderr,"  number of temperature changes ");
      //   dots(stderr,21);
      //   fprintf(stderr," nT = %3d\n", nT[lc] );
      //   sprintf(errMsg,"\n  error: valid ranges for nT is 0 ... %d \n", nE );
      //   errorMsg(errMsg);
      //   exit(160);
      // }
      // for (i=1; i <= nT[lc]; i++) { /* ! local element coordinates ! */
      //   sfrv=fscanf(fp,"%d", &n );
      //   if (sfrv != 1) sferr("frame element number in temperature load data");
      //   if ( n < 1 || n > nE ) {
      //       sprintf(errMsg,"\n  error in temperature loads: frame element number %d is out of range\n",n);
      //       errorMsg(errMsg);
      //       exit(161);
      //   }
      //   T[lc][i][1] = (double) n;
      //   for (l=2; l<=8; l++) {
      //       sfrv=fscanf(fp,"%f", &T[lc][i][l] );
      //       if (sfrv != 1) sferr("value in temperature load data");
      //   }
      //   a  = T[lc][i][2];
      //   hy = T[lc][i][3];
      //   hz = T[lc][i][4];

      //   if ( hy < 0 || hz < 0 ) {
      //       sprintf(errMsg,"\n  error in thermal load data: section dimension < 0\n   Frame element number: %d  hy: %f  hz: %f\n", n,hy,hz);
      //       errorMsg(errMsg);
      //       exit(162);
      //   }

      //   Nx2 = (a/4.0)*( T[lc][i][5]+T[lc][i][6]+T[lc][i][7]+T[lc][i][8])*E[n]*Ax[n];
      //   Nx1 = -Nx2;
      //   Vy1 = Vy2 = Vz1 = Vz2 = 0.0;
      //   Mx1 = Mx2 = 0.0;
      //   My1 =  (a/hz)*(T[lc][i][8]-T[lc][i][7])*E[n]*Iy[n];
      //   My2 = -My1;
      //   Mz1 =  (a/hy)*(T[lc][i][5]-T[lc][i][6])*E[n]*Iz[n];
      //   Mz2 = -Mz1;

      //   n1 = J1[n]; n2 = J2[n];

      //   coord_trans ( xyz, L[n], n1, n2,
      //       &t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p[n] );

      //   /* {F} = [T]'{Q} */
      //   feF_temp[lc][n][1]  += ( Nx1*t1 + Vy1*t4 + Vz1*t7 );
      //   feF_temp[lc][n][2]  += ( Nx1*t2 + Vy1*t5 + Vz1*t8 );
      //   feF_temp[lc][n][3]  += ( Nx1*t3 + Vy1*t6 + Vz1*t9 );
      //   feF_temp[lc][n][4]  += ( Mx1*t1 + My1*t4 + Mz1*t7 );
      //   feF_temp[lc][n][5]  += ( Mx1*t2 + My1*t5 + Mz1*t8 );
      //   feF_temp[lc][n][6]  += ( Mx1*t3 + My1*t6 + Mz1*t9 );

      //   feF_temp[lc][n][7]  += ( Nx2*t1 + Vy2*t4 + Vz2*t7 );
      //   feF_temp[lc][n][8]  += ( Nx2*t2 + Vy2*t5 + Vz2*t8 );
      //   feF_temp[lc][n][9]  += ( Nx2*t3 + Vy2*t6 + Vz2*t9 );
      //   feF_temp[lc][n][10] += ( Mx2*t1 + My2*t4 + Mz2*t7 );
      //   feF_temp[lc][n][11] += ( Mx2*t2 + My2*t5 + Mz2*t8 );
      //   feF_temp[lc][n][12] += ( Mx2*t3 + My2*t6 + Mz2*t9 );
      // }             /* end thermal loads    */

      /* debugging ...  check feF's prior to asembly
      for (n=1; n<=nE; n++) {
        printf("n=%d ", n);
        for (l=1;l<=12;l++) {
            if (feF_mech[lc][n][l] != 0)
               printf(" feF %d = %9.2e ", l, feF_mech[lc][n][l] );
        }
        printf("\n");
      }
      */

        for (n=1; n<=nE; n++) {
            n1 = J1[n];    n2 = J2[n];
            for (i=1; i<= 6; i++) F_mech[lc][6*n1- 6+i] += feF_mech[lc][n][i];
            for (i=7; i<=12; i++) F_mech[lc][6*n2-12+i] += feF_mech[lc][n][i];
            for (i=1; i<= 6; i++) F_temp[lc][6*n1- 6+i] += feF_temp[lc][n][i];
            for (i=7; i<=12; i++) F_temp[lc][6*n2-12+i] += feF_temp[lc][n][i];
        }

        nD[lc] = pD.nD;

        for (i=1; i <= nD[lc]; i++) {
            j = pD.N[i-1];

            Dp[lc][6*j-5] = pD.Dx[i-1];
            Dp[lc][6*j-4] = pD.Dy[i-1];
            Dp[lc][6*j-3] = pD.Dz[i-1];
            Dp[lc][6*j-2] = pD.Dxx[i-1];
            Dp[lc][6*j-1] = pD.Dyy[i-1];
            Dp[lc][6*j] = pD.Dzz[i-1];

            for (l=5; l >=0; l--) {
                if ( r[6*j-l] == 0 && Dp[lc][6*j-l] != 0.0 ) {
                    sprintf(errMsg," Initial displacements can be prescribed only at restrained coordinates\n  node: %d  dof: %d  r: %d\n",
                    j, 6-l, r[6*j-l] );
                    errorMsg(errMsg);
                    exit(171);
                }
            }
        }

    }                   /* end load-case loop */




    /* solve the problem    */
    int debug=0;    // 1: debugging screen output, 0: none
    int verbose = 0;  // 1: copious screen output, 0: none
    int ok=1;       // number of (-ve) diag. terms of L D L'
    int iter=0;     // number of iterations
    double error = 1.0;    // rms equilibrium error and reactions
    double tol = 1.0e-9;   // tolerance for modal convergence
    double rms_resid=1.0;  // root mean square of residual displ. error

    srand(time(NULL));
    for (lc=1; lc<=nL; lc++) { /* begin load case analysis loop */

        /*  initialize displacements and displ. increment to {0}  */
        for (i=1; i<=DoF; i++)  D[i] = dD[i] = 0.0;

        /*  initialize internal forces to {0}   */
        for (i=1; i<=nE; i++)   for (j=1;j<=12;j++) Q[i][j] = 0.0;

        /*  assemble stiffness matrix [K({D}^(i))], {D}^(0)={0} (i=0) */
        assemble_K ( K, DoF, nE, xyz, rj, L, Le, N1, N2,
                    Ax, Asy, Asz, Jx,Iy,Iz, E, G, p,
                    shear, geom, Q, debug );

#ifdef MATRIX_DEBUG
        save_dmatrix ( "Ku", K, 1,DoF, 1,DoF, 0, "w" ); // unloaded stiffness matrix
#endif

        /* first apply temperature loads only, if there are any ... */
        if (nT[lc] > 0) {
            // if ( verbose )
            //     fprintf(stdout," Linear Elastic Analysis ... Temperature Loads\n");

            /*  solve {F_t} = [K({D=0})] * {D_t} */
            solve_system(K,dD,F_temp[lc],DoF,q,r,&ok,verbose,&rms_resid);

            /* increment {D_t} = {0} + {D_t} temp.-induced displ */
            for (i=1; i<=DoF; i++)  if (q[i]) D[i] += dD[i];

            if (geom) {
             /* compute   {Q}={Q_t} ... temp.-induced forces     */
             element_end_forces ( Q, nE, xyz, L, Le, N1,N2,
                Ax, Asy,Asz, Jx,Iy,Iz, E,G, p, D, shear, geom );

             /* assemble temp.-stressed stiffness [K({D_t})]     */
             assemble_K ( K, DoF, nE, xyz, rj, L, Le, N1, N2,
                        Ax,Asy,Asz, Jx,Iy,Iz, E, G, p,
                        shear,geom, Q, debug );
            }
        }

        /* ... then apply mechanical loads only, if there are any ... */
        if ( nF[lc]>0 || nU[lc]>0 || nW[lc]>0 || nP[lc]>0 || nD[lc]>0 ||
             gX[lc] != 0 || gY[lc] != 0 || gZ[lc] != 0 ) {

            for (i=1; i<=DoF; i++)  if (r[i]) dD[i] = Dp[lc][i];

            /*  solve {F_m} = [K({D_t})] * {D_m}    */
            solve_system(K,dD,F_mech[lc],DoF,q,r,&ok,verbose,&rms_resid);

            /* combine {D} = {D_t} + {D_m}  */
            for (i=1; i<=DoF; i++) {
                if (q[i])   D[i] += dD[i];
                else {      D[i] = Dp[lc][i]; dD[i] = 0.0; }
            }
        }

        /*  combine {F} = {F_t} + {F_m} */
        for (i=1; i<=DoF; i++)  F[lc][i] = F_temp[lc][i] + F_mech[lc][i];

        /*  element forces {Q} for displacements {D}    */
        element_end_forces ( Q, nE, xyz, L, Le, N1,N2,
                Ax, Asy,Asz, Jx,Iy,Iz, E,G, p, D, shear, geom );


        /* initialize Broyden secant stiffness matrix, Ks */
/*
        if ( geom ) {
            Ks  = dmatrix( 1, DoF, 1, DoF );
            for (i=1;i<=DoF;i++) {
                for(j=i;j<=DoF;j++) {
                    Ks[i][j]=Ks[j][i]=K[i][j];
                }
            }
        }
*/

        /* quasi Newton-Raphson iteration for geometric nonlinearity */
        ok = 0; iter = 0; error = 1.0;  /* re-initialize */
        while ( geom && error > tol && iter < 500 && ok >= 0) {

            ++iter;

            /*  assemble stiffness matrix [K({D}^(i))]  */
            assemble_K ( K, DoF, nE, xyz, rj, L, Le, N1, N2,
                Ax,Asy,Asz, Jx,Iy,Iz, E, G, p,
                shear,geom, Q, debug );

            /*  compute equilibrium error, {dF}, at iteration i */
            /*  {dF}^(i) = {F} - [K({D}^(i))]*{D}^(i) */
            /*  convergence criteria = || {dF}^(i) ||  /  || F || */
            error = equilibrium_error ( dF, F[lc], K, D, DoF, q );


            /*  solve {dF}^(i) = [K({D}^(i))] * {dD}^(i) */
            solve_system(K,dD,dF,DoF,q,r,&ok,verbose,&rms_resid);

            if ( ok < 0 ) { /*  K is not pos.def.  */
                fprintf(stderr,"   The stiffness matrix is not pos-def. \n");
                fprintf(stderr,"   Reduce loads and re-run the analysis.\n");
                break;
            }

            /*  increment {D}^(i+1} = {D}^{i} + {dD}^(i) */
            for (i=1; i<=DoF; i++)  if (q[i])   D[i] += dD[i];

            /*  compute element forces {Q} for displacements {D}^(i) */
            element_end_forces ( Q, nE, xyz, L, Le, N1,N2,
                Ax, Asy,Asz, Jx,Iy,Iz, E,G, p, D, shear, geom );

        }           /* end quasi Newton-Raphson iteration */

        compute_reaction_forces( F[lc], K, D, DoF, r );

        add_feF ( xyz, L, N1,N2, p, Q, feF_temp[lc], feF_mech[lc],
                nE, DoF, F[lc], r, verbose );




        if ( ok < 0 ) {
            printf("  * The Stiffness Matrix is not positive-definite *\n");
            printf("    Check that all six rigid-body translations are restrained\n");
            printf("    If geometric stiffness is included, reduce the loads.\n");
            /*   return; */
        }

        printf("\nL O A D   C A S E   %d   O F   %d  ... \n\n", lc, nL);


        for (j=1; j<= nN; j++) {
            displacements[lc-1].node[j-1] = j;
            displacements[lc-1].x[j-1] = D[6*j-5];
            displacements[lc-1].y[j-1] = D[6*j-4];
            displacements[lc-1].z[j-1] = D[6*j-3];
            displacements[lc-1].xrot[j-1] = D[6*j-2];
            displacements[lc-1].yrot[j-1] = D[6*j-1];
            displacements[lc-1].zrot[j-1] = D[6*j];
        }

        printf("here1\n");


        for (n=1; n<= nE; n++) {
            forces[lc-1].element[2*n-2] = n;
            forces[lc-1].node[2*n-2] = J1[n];
            forces[lc-1].Nx[2*n-2] = Q[n][1];
            forces[lc-1].Vy[2*n-2] = Q[n][2];
            forces[lc-1].Vz[2*n-2] = Q[n][3];
            forces[lc-1].Txx[2*n-2] = Q[n][4];
            forces[lc-1].Myy[2*n-2] = Q[n][5];
            forces[lc-1].Mzz[2*n-2] = Q[n][6];

            forces[lc-1].element[2*n-1] = n;
            forces[lc-1].node[2*n-1] = J2[n];
            forces[lc-1].Nx[2*n-1] = Q[n][7];
            forces[lc-1].Vy[2*n-1] = Q[n][8];
            forces[lc-1].Vz[2*n-1] = Q[n][9];
            forces[lc-1].Txx[2*n-1] = Q[n][10];
            forces[lc-1].Myy[2*n-1] = Q[n][11];
            forces[lc-1].Mzz[2*n-1] = Q[n][12];

        }


        printf("here1\n");

        for (j=1; j<=nN; j++) {

            reactionForces[lc-1].node[j-1] = j;
            reactionForces[lc-1].Fx[j-1] = F[lc][6*j-5];
            reactionForces[lc-1].Fy[j-1] = F[lc][6*j-4];
            reactionForces[lc-1].Fz[j-1] = F[lc][6*j-3];
            reactionForces[lc-1].Mxx[j-1] = F[lc][6*j-2];
            reactionForces[lc-1].Myy[j-1] = F[lc][6*j-1];
            reactionForces[lc-1].Mzz[j-1] = F[lc][6*j];

        }
        printf("R M S    R E L A T I V E    E Q U I L I B R I U M    E R R O R: %9.3e\n", rms_resid );


        // write_static_results ( fp, nN,nE,nL, lc, DoF, N1,N2,
        //         F[lc], D,r,Q, rms_resid, ok, axial_sign );

        // if ( filetype == 1 ) {      // .CSV format output
        //     write_static_csv(OUT_file, title,
        //         nN,nE,nL,lc, DoF, N1,N2, F[lc], D,r,Q, error, ok );
        // }

        // if ( filetype == 2 ) {      // .m matlab format output
        //     write_static_mfile (OUT_file, title, nN,nE,nL,lc, DoF,
        //             N1,N2, F[lc], D,r,Q, error, ok );
        // }

/*
 *      if ( verbose )
 *       printf("\n   If the program pauses here for very long,"
 *       " hit CTRL-C to stop execution, \n"
 *       "    reduce exagg_static in the Input Data,"
 *       " and re-run the analysis. \n");
 */

        // TODO:
        // write_internal_forces ( infcpath, lc, nL, title, dx, xyz,
        //             Q, nN, nE, L, N1, N2,
        //             Ax, Asy, Asz, Jx, Iy, Iz, E, G, p,
        //             d, gX[lc], gY[lc], gZ[lc],
        //             nU[lc],U[lc],nW[lc],W[lc],nP[lc],P[lc],
        //             D, shear, error );

        // static_mesh ( IN_file, infcpath, meshpath, plotpath, title,
        //             nN, nE, nL, lc, DoF,
        //             xyz, L, N1,N2, p, D,
        //             exagg_static, D3_flag, anlyz, dx );

    } /* end load case loop */


//     if (nM > 0) { /* carry out modal analysis */

//         if(verbose & anlyz) fprintf(stdout,"\n\n Modal Analysis ...\n");

//         nM_calc = (nM+8)<(2*nM) ? nM+8 : 2*nM;      /* Bathe */

//         M   = dmatrix(1,DoF,1,DoF);
//         f   = dvector(1,nM_calc);
//         V   = dmatrix(1,DoF,1,nM_calc);

//         assemble_M ( M, DoF, nN, nE, xyz, rj, L, N1,N2,
//                 Ax, Jx,Iy,Iz, p, d, EMs, NMs, NMx, NMy, NMz,
//                 lump, debug );

// #ifdef MATRIX_DEBUG
//         save_dmatrix ( "Mf", M, 1,DoF, 1,DoF, 0, "w" ); /* free mass matrix */
// #endif

//         for (j=1; j<=DoF; j++) { /*  compute traceK and traceM */
//             if ( !r[j] ) {
//                 traceK += K[j][j];
//                 traceM += M[j][j];
//             }
//         }
//         for (i=1; i<=DoF; i++) { /*  modify K and M for reactions    */
//             if ( r[i] ) {   /* apply reactions to upper triangle */
//                 K[i][i] = traceK * 1e4;
//                 M[i][i] = traceM;
//                 for (j=i+1; j<=DoF; j++)
//                     K[j][i]=K[i][j]=M[j][i]=M[i][j] = 0.0;
//             }
//         }

//         if ( write_matrix ) {   /* write Kd and Md matrices */
//             save_ut_dmatrix ( "Kd", K, DoF, "w" );/* dynamic stff matx */
//             save_ut_dmatrix ( "Md", M, DoF, "w" );/* dynamic mass matx */
//         }

//         if ( anlyz ) {  /* subspace or stodola methods */
//             if( Mmethod == 1 )
//                 subspace( K, M, DoF, nM_calc, f, V, tol,shift,&iter,&ok, verbose );
//             if( Mmethod == 2 )
//                 stodola ( K, M, DoF, nM_calc, f, V, tol,shift,&iter,&ok, verbose );

//             for (j=1; j<=nM_calc; j++) f[j] = sqrt(f[j])/(2.0*PI);

//             write_modal_results ( fp, nN,nE,nI, DoF, M,f,V,
//                     total_mass, struct_mass,
//                     iter, sumR, nM, shift, lump, tol, ok );
//         }
//     }

    // fprintf(fp,"\n");
    // fclose (fp);

    // if ( nM > 0 && anlyz ) {    /* write modal analysis results */

    //     modal_mesh ( IN_file, meshpath, modepath, plotpath, title,
    //             nN,nE, DoF, nM, xyz, L, N1,N2, p,
    //             M, f, V, exagg_modal, D3_flag, anlyz );

    //     animate ( IN_file, meshpath, modepath, plotpath, title,anim,
    //             nN,nE, DoF, nM, xyz, L, p, N1,N2, f,
    //             V, exagg_modal, D3_flag, pan );
    // }

    // if ( nC > 0 ) {     /* matrix condensation of stiffness and mass */

    //     if ( verbose ) fprintf(stdout,"\n Matrix Condensation ...\n");

    //     if(Cdof > nM && Cmethod == 3){
    //         fprintf(stderr,"  Cdof > nM ... Cdof = %d  nM = %d \n",
    //              Cdof, nM );
    //         fprintf(stderr,"  The number of condensed degrees of freedom");
    //         fprintf(stderr," may not exceed the number of computed modes");
    //         fprintf(stderr," when using dynamic condensation.\n");
    //         exit(94);
    //     }

    //     Kc = dmatrix(1,Cdof,1,Cdof);
    //     Mc = dmatrix(1,Cdof,1,Cdof);

    //     if ( m[1] > 0 && nM > 0 )   Cfreq = f[m[1]];

    //     if ( Cmethod == 1 ) {   /* static condensation only */
    //         condense(K, DoF, c, Cdof, Kc, verbose );
    //         if ( verbose )
    //             fprintf(stdout,"   static condensation of K complete\n");
    //     }
    //     if ( Cmethod == 2 ) {
    //         guyan(M, K, DoF, c, Cdof, Mc,Kc, Cfreq, verbose );
    //         if ( verbose ) {
    //             fprintf(stdout,"   Guyan condensation of K and M complete");
    //             fprintf(stdout," ... dynamics matched at %f Hz.\n", Cfreq );
    //         }
    //     }
    //     if ( Cmethod == 3 && nM > 0 ) {
    //         dyn_conden(M,K, DoF, r, c, Cdof, Mc,Kc, V,f, m, verbose );
    //         if ( verbose )
    //             fprintf(stdout,"   dynamic condensation of K and M complete\n");
    //     }
    //     save_dmatrix("Kc", Kc, 1,Cdof, 1,Cdof, 0, "w" );
    //     save_dmatrix("Mc", Mc, 1,Cdof, 1,Cdof, 0, "w" );

    //     free_dmatrix(Kc, 1,Cdof,1,Cdof );
    //     free_dmatrix(Mc, 1,Cdof,1,Cdof );
    // }


    /* deallocate memory used for each frame analysis variable */
    // TODO
    // deallocate ( nN, nE, nL, nF, nU, nW, nP, nT, DoF,
    //         xyz, rj, L, Le, N1, N2, q,r,
    //         Ax, Asy, Asz, Jx, Iy, Iz, E, G, p,
    //         U,W,P,T, Dp, F_mech, F_temp,
    //         feF_mech, feF_temp, F, dF,
    //         K, Q, D, dD,
    //         d,EMs,NMs,NMx,NMy,NMz, c, m
    // );

    // if ( verbose ) fprintf(stdout,"\n");

    // if ( argc == 1 ) { /* wait for keyboard entry to close the terminal */
    //    fprintf(stderr," The Output Data was appended to %s \n", OUT_file );
    //    fprintf(stderr," A Gnuplot script was written to %s \n", plotpath );
    //    fprintf(stderr," Press the 'Enter' key to close.\n");
    //    (void) getchar();    // clear the buffer ??
    //    while( !getchar() ) ;    // wait for the Enter key to be hit
    // }
    // color(0);

    return;


}



