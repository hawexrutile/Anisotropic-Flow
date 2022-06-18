/* art1f.f -- translated by f2c (version 20160102).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    real r__[450003]	/* was [3][150001] */;
} aa_;

#define aa_1 aa_

struct {
    real p[450003]	/* was [3][150001] */;
} bb_;

#define bb_1 bb_

struct {
    real e[150001];
} cc_;

#define cc_1 cc_

struct {
    real rho[82369]	/* was [41][41][49] */, rhop[82369]	/* was [41][
	    41][49] */, rhon[82369]	/* was [41][41][49] */;
} dd_;

#define dd_1 dd_

struct {
    integer id[150001], lb[150001];
} ee_;

#define ee_1 ee_

struct {
    real proper[150001];
} hh_;

#define hh_1 hh_

struct {
    real f[2342277]	/* was [9][9][17][9][9][21] */;
} ff_;

#define ff_1 ff_

struct {
    real dx, dy, dz, dpx, dpy, dpz;
} gg_;

#define gg_1 gg_

struct {
    integer nstar, ndirct;
    real dir;
} input_;

#define input_1 input_

struct {
    real prho[2009]	/* was [41][49] */;
} pp_;

#define pp_1 pp_

struct {
    real phrho[2401]	/* was [49][49] */;
} qq_;

#define qq_1 qq_

struct {
    integer massr[2];
} rr_;

#define rr_1 rr_

struct {
    integer inout[20];
} ss_;

#define ss_1 ss_

struct {
    integer zta, zpr;
} zz_;

#define zz_1 zz_

struct {
    integer num;
} run_;

#define run_1 run_

struct {
    real tkaon[7], ekaon[14007]	/* was [7][2001] */;
} kkk_;

#define kkk_1 kkk_

struct {
    real ak[5400]	/* was [3][50][36] */, speck[12600]	/* was [50][
	    36][7] */;
    integer mf;
} kaon_;

#define kaon_1 kaon_

struct {
    real xarray[1001], earray[1001];
} table_;

#define table_1 table_

struct {
    integer masspr, massta, iseed, iavoid;
    real dt;
} input1_;

#define input1_1 input1_

struct {
    real pirho[82369]	/* was [41][41][49] */;
} ddpi_;

#define ddpi_1 ddpi_

struct {
    real pel[82369]	/* was [41][41][49] */, rxy[82369]	/* was [41][
	    41][49] */;
} tt_;

#define tt_1 tt_

struct {
    integer ilab, manyb, ntmax, icoll, insys, ipot, mode, imomen, nfreq, 
	    icflow, icrho, icou, kpoten, kmul;
} input2_;

#define input2_1 input2_

struct {
    real plab, elab, zeropt, b0, bi, bm, dencut, cycbox;
} input3_;

#define input3_1 input3_

struct {
    real arpar1[100];
    integer iapar2[50];
    real arint1[100];
    integer iaint2[50];
} arprnt_;

#define arprnt_1 arprnt_

struct {
    real pro1[150001]	/* was [150001][1] */;
} arercp_;

#define arercp_1 arercp_

struct {
    integer multi1[1];
} arerc1_;

#define arerc1_1 arerc1_

struct {
    integer ityp1[150001]	/* was [150001][1] */;
    real gx1[150001]	/* was [150001][1] */, gy1[150001]	/* was [
	    150001][1] */, gz1[150001]	/* was [150001][1] */, ft1[150001]	
	    /* was [150001][1] */, px1[150001]	/* was [150001][1] */, py1[
	    150001]	/* was [150001][1] */, pz1[150001]	/* was [
	    150001][1] */, ee1[150001]	/* was [150001][1] */, xm1[150001]	
	    /* was [150001][1] */;
} arprc1_;

#define arprc1_1 arprc1_

struct {
    integer itimeh;
    real bimp;
} lastt_;

#define lastt_1 lastt_

struct {
    real efrm;
    integer npart1, npart2;
    real epsipz, epsipt, pzproj, pztarg;
} snn_;

#define snn_1 snn_

struct {
    integer lblast[150001];
    real xlast[600004]	/* was [4][150001] */, plast[600004]	/* was [4][
	    150001] */;
    integer nlast;
} hbt_;

#define hbt_1 hbt_

struct {
    integer nsav, iksdcy;
} resdcy_;

#define resdcy_1 resdcy_

struct {
    integer nseed;
} rndf77_;

#define rndf77_1 rndf77_

struct {
    real ftsv[150001], ftsvt[150001]	/* was [150001][1] */;
} ftmax_;

#define ftmax_1 ftmax_

struct {
    real dpertt[150001]	/* was [150001][1] */, dpertp[150001], dplast[150001],
	     dpdcy[150001], dpdpi[150001]	/* was [150001][1] */, dpt[
	    150001]	/* was [150001][1] */, dpp1[150001]	/* was [
	    150001][1] */, dppion[150001]	/* was [150001][1] */;
} dpert_;

#define dpert_1 dpert_

struct {
    real hipr1[100];
    integer ihpr2[50];
    real hint1[100];
    integer ihnt2[50];
} hparnt_;

#define hparnt_1 hparnt_

struct {
    integer nnn;
} nn_;

#define nn_1 nn_

struct {
    real betax, betay, betaz, gamma;
} bg_;

#define bg_1 bg_

struct {
    real rpion[450003]	/* was [3][150001][1] */;
} pa_;

#define pa_1 pa_

struct {
    real ppion[450003]	/* was [3][150001][1] */;
} pb_;

#define pb_1 pb_

struct {
    real epion[150001]	/* was [150001][1] */;
} pc_;

#define pc_1 pc_

struct {
    integer lpion[150001]	/* was [150001][1] */;
} pd_;

#define pd_1 pd_

struct {
    real propi[150001]	/* was [150001][1] */;
} pe_;

#define pe_1 pe_

struct {
    integer lb1;
    real px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, 
	    dp1n;
} leadng_;

#define leadng_1 leadng_

struct {
    real tfdcy[150001], tfdpi[150001]	/* was [150001][1] */, tft[150001];
} tdecay_;

#define tdecay_1 tdecay_

struct ppbmas_1_ {
    integer niso[15], nstate;
    real ppbm[30]	/* was [15][2] */, thresh[15], weight[15];
};

#define ppbmas_1 (*(struct ppbmas_1_ *) &ppbmas_)

struct ppb1_1_ {
    real ene, factr2[6], fsum, ppinnb, s, wtot;
};

#define ppb1_1 (*(struct ppb1_1_ *) &ppb1_)

struct {
    real pprr, ppee, pppe, rpre, xopoe, rree;
} ppmm_;

#define ppmm_1 ppmm_

struct {
    real em2;
    integer lb2;
} dpi_;

#define dpi_1 dpi_

struct {
    integer iphidcy;
    real pttrig;
    integer ntrig, maxmiss, ipi0dcy;
} phidcy_;

#define phidcy_1 phidcy_

struct {
    integer idpert, npertd, idxsec;
} para8_;

#define para8_1 para8_

struct {
    real bxx[82369]	/* was [41][41][49] */, byy[82369]	/* was [41][
	    41][49] */, bzz[82369]	/* was [41][41][49] */;
} bbb_;

#define bbb_1 bbb_

struct {
    integer iaevt, iarun, miss;
} arevt_;

#define arevt_1 arevt_

struct {
    integer lbnn1, lbnn2, lbnd1, lbnd2, lbns1, lbns2, lbnp1, lbnp2, lbdd1, 
	    lbdd2, lbds1, lbds2, lbdp1, lbdp2, lbss1, lbss2, lbsp1, lbsp2, 
	    lbpp1, lbpp2;
} dpifsl_;

#define dpifsl_1 dpifsl_

struct {
    real xmnn1, xmnn2, xmnd1, xmnd2, xmns1, xmns2, xmnp1, xmnp2, xmdd1, xmdd2,
	     xmds1, xmds2, xmdp1, xmdp2, xmss1, xmss2, xmsp1, xmsp2, xmpp1, 
	    xmpp2;
} dpifsm_;

#define dpifsm_1 dpifsm_

struct {
    real sdmel, sdmnn, sdmnd, sdmns, sdmnp, sdmdd, sdmds, sdmdp, sdmss, sdmsp,
	     sdmpp;
} dpisig_;

#define dpisig_1 dpisig_

/* Initialized data */

struct {
    integer e_1[15];
    integer fill_2[1];
    real e_3[45];
    integer fill_4[15];
    } ppbmas_ = { 1, 2, 1, 16, 16, 4, 4, 64, 4, 4, 32, 32, 4, 8, 4, {0}, 
	    .93828f, .93828f, .939457f, .93828f, .939457f, .93828f, .939457f, 
	    1.232f, .93828f, .939457f, 1.232f, 1.232f, 1.44f, 1.44f, 1.535f, 
	    .93828f, .939457f, .939457f, 1.232f, 1.232f, 1.44f, 1.44f, 1.232f,
	     1.535f, 1.535f, 1.44f, 1.535f, 1.44f, 1.535f, 1.535f, 1.87656f, 
	    1.877737f, 1.878914f, 2.17028f, 2.171457f, 2.37828f, 2.379457f, 
	    2.464f, 2.47328f, 2.474457f, 2.672f, 2.767f, 2.88f, 2.975f, 3.07f 
	    };

struct {
    integer fill_1[1];
    real e_2[6];
    integer fill_3[4];
    } ppb1_ = { {0}, 0.f, 1.f, .117f, .00327f, 3.58e-5f, 1.93e-7f };


/* Table of constant values */

static doublereal c_b5 = .33333333333333331;
static integer c__1 = 1;
static integer c__0 = 0;
static integer c__2 = 2;
static real c_b106 = 5.f;
static integer c__9 = 9;
static integer c__3 = 3;
static real c_b122 = .6f;
static real c_b123 = .3f;
static integer c__5 = 5;
static real c_b182 = 0.f;
static real c_b183 = 1.f;
static real c_b195 = 1.232f;
static real c_b197 = 1.08f;
static real c_b201 = 1.44f;
static integer c__4 = 4;
static real c_b254 = 2.f;
static real c_b255 = -1.f;
static doublereal c_b397 = .16666666666666666;
static doublereal c_b398 = -.333;
static doublereal c_b407 = .67;
static doublereal c_b413 = .333;
static real c_b725 = .77f;
static real c_b727 = .28f;
static doublereal c_b783 = .141;
static doublereal c_b784 = 1.9;
static real c_b842 = 2.01f;
static doublereal c_b927 = .7;
static doublereal c_b932 = .048383999999999983;
static doublereal c_b933 = 1.5;
static doublereal c_b941 = -.519;
static doublereal c_b942 = -.534;
static doublereal c_b943 = -.304;
static doublereal c_b944 = -.786;
static doublereal c_b945 = -.277;
static doublereal c_b1071 = 3.3;

/* ....................art1f.f */
/* ************************************* */

/*                           PROGRAM ART1.0 */

/*        A relativistic transport (ART) model for heavy-ion collisions */

/*   sp/01/04/2002 */
/*   calculates K+K- from phi decay, dimuons from phi decay */
/*   has finite baryon density & possibilites of varying Kaon */
/*   in-medium mass in phiproduction-annhilation channel only. */


/* RELEASING DATE: JAN., 1997 */
/* ************************************** */

/* Bao-An Li & Che Ming Ko */
/* Cyclotron Institute, Texas A&M University. */
/* Phone: (409) 845-1411 */
/* e-mail: Bali@comp.tamu.edu & Ko@comp.tamu.edu */
/* http://wwwcyc.tamu.edu/~bali */
/* ************************************** */
/* Speical notice on the limitation of the code: */

/* (1) ART is a hadronic transport model */

/* (2) E_beam/A <= 15 GeV */

/* (3) The mass of the colliding system is limited by the dimensions of arrays */
/*    which can be extended purposely. Presently the dimensions are large enough */
/*     for running Au+Au at 15 GeV/A. */

/* (4) The production and absorption of antiparticles (e.g., ki-, anti-nucleons, */
/*     etc) are not fully included in this version of the model. They, however, */
/*     have essentially no effect on the reaction dynamics and observables */
/*     related to nucleons, pions and kaons (K+) at and below AGS energies. */

/* (5) Bose enhancement for mesons and Pauli blocking for fermions are */
/*     turned off. */

/* ******************************** */

/* USEFUL REFERENCES ON PHYSICS AND NUMERICS OF NUCLEAR TRANSPORT MODELS: */
/*     G.F. BERTSCH AND DAS GUPTA, PHYS. REP. 160 (1988) 189. */
/*     B.A. LI AND W. BAUER, PHYS. REV. C44 (1991) 450. */
/*     B.A. LI, W. BAUER AND G.F. BERTSCH, PHYS. REV. C44 (1991) 2095. */
/*     P. DANIELEWICZ AND G.F. BERTSCH, NUCL. PHYS. A533 (1991) 712. */

/* MAIN REFERENCES ON THIS VERSION OF ART MODEL: */
/*     B.A. LI AND C.M. KO, PHYS. REV. C52 (1995) 2037; */
/*                          NUCL. PHYS. A601 (1996) 457. */

/* ********************************* */
/* ********************************* */
/*  VARIABLES IN INPUT-SECTION:                                               * */
/*                                                                      * */
/*  1) TARGET-RELATED QUANTITIES                                        * */
/*       MASSTA, ZTA -  TARGET MASS IN AMU, TARGET CHARGE  (INTEGER)    * */
/*                                                                      * */
/*  2) PROJECTILE-RELATED QUANTITIES                                    * */
/*       MASSPR, ZPR -  PROJECTILE MASS IN AMU, PROJ. CHARGE(INTEGER)   * */
/*       ELAB     -  BEAM ENERGY IN [MEV/NUCLEON]               (REAL)  * */
/*       ZEROPT   -  DISPLACEMENT OF THE SYSTEM IN Z-DIREC. [FM](REAL)  * */
/*       B        -  IMPACT PARAMETER [FM]                      (REAL)  * */
/*                                                                      * */
/*  3) PROGRAM-CONTROL PARAMETERS                                       * */
/*       ISEED    -  SEED FOR RANDOM NUMBER GENERATOR        (INTEGER)  * */
/*       DT       -  TIME-STEP-SIZE [FM/C]                      (REAL)  * */
/*       NTMAX    -  TOTAL NUMBER OF TIMESTEPS               (INTEGER)  * */
/*       ICOLL    -  (= 1 -> MEAN FIELD ONLY,                           * */
/*                -   =-1 -> CACADE ONLY, ELSE FULL ART)     (INTEGER)  * */
/*       NUM      -  NUMBER OF TESTPARTICLES PER NUCLEON     (INTEGER)  * */
/*       INSYS    -  (=0 -> LAB-SYSTEM, ELSE C.M. SYSTEM)    (INTEGER)  * */
/*       IPOT     -  1 -> SIGMA=2; 2 -> SIGMA=4/3; 3 -> SIGMA=7/6       * */
/*                   IN MEAN FIELD POTENTIAL                 (INTEGER)  * */
/*       MODE     -  (=1 -> interpolation for pauli-blocking,           * */
/*                    =2 -> local lookup, other -> unblocked)(integer)  * */
/*       DX,DY,DZ -  widths of cell for paulat in coor. sp. [fm](real)  * */
/*       DPX,DPY,DPZ-widths of cell for paulat in mom. sp.[GeV/c](real) * */
/*       IAVOID   -  (=1 -> AVOID FIRST COLL. WITHIN SAME NUCL.         * */
/*                    =0 -> ALLOW THEM)                      (INTEGER)  * */
/*       IMOMEN   -  FLAG FOR CHOICE OF INITIAL MOMENTUM DISTRIBUTION   * */
/*                   (=1 -> WOODS-SAXON DENSITY AND LOCAL THOMAS-FERMI  * */
/*                    =2 -> NUCLEAR MATTER DEN. AND LOCAL THOMAS-FERMI  * */
/*                    =3 -> COHERENT BOOST IN Z-DIRECTION)   (INTEGER)  * */
/*  4) CONTROL-PRINTOUT OPTIONS                                         * */
/*       NFREQ    -  NUMBER OF TIMSTEPS AFTER WHICH PRINTOUT            * */
/*                   IS REQUIRED OR ON-LINE ANALYSIS IS PERFORMED       * */
/*       ICFLOW      =1 PERFORM ON-LINE FLOW ANALYSIS EVERY NFREQ STEPS * */
/*       ICRHO       =1 PRINT OUT THE BARYON,PION AND ENERGY MATRIX IN  * */
/*                      THE REACTION PLANE EVERY NFREQ TIME-STEPS       * */
/*  5) */
/*       CYCBOX   -  ne.0 => cyclic boundary conditions;boxsize CYCBOX  * */

/* ********************************* */
/*               Lables of particles used in this code                     * */
/* ********************************* */

/*         LB(I) IS USED TO LABEL PARTICLE'S CHARGE STATE */

/*         LB(I)   = */
/* lin-11/07/00: */
/*                -30 K*- */
/* lin-8/29/00 */
/*                -13 anti-N*(+1)(1535),s_11 */
/*                -12 anti-N*0(1535),s_11 */
/*                 -11 anti-N*(+1)(1440),p_11 */
/*                 -10 anti-N*0(1440), p_11 */
/*                  -9 anti-DELTA+2 */
/*                  -8 anti-DELTA+1 */
/*                  -7 anti-DELTA0 */
/*                  -6 anti-DELTA-1 */
/* lin-8/29/00-end */
/* bali2/7/99 */
/*                  -2 antineutron */
/*                             -1       antiproton */
/* bali2/7/99 end */
/*                   0 eta */
/*                        1 PROTON */
/*                   2 NUETRON */
/*                   3 PION- */
/*                   4 PION0 */
/*                   5 PION+ */
/*                   6 DELTA-1 */
/*                   7 DELTA0 */
/*                   8 DELTA+1 */
/*                   9 DELTA+2 */
/*                   10 N*0(1440), p_11 */
/*                   11 N*(+1)(1440),p_11 */
/*                  12 N*0(1535),s_11 */
/*                  13 N*(+1)(1535),s_11 */
/*                  14 LAMBDA */
/*                   15 sigma-, since we used isospin averaged xsection for */
/*                   16 sigma0  sigma associated K+ production, sigma0 and */
/*                   17 sigma+  sigma+ are counted as sigma- */
/*                   21 kaon- */
/*                   23 KAON+ */
/*                   24 kaon0 */
/*                   25 rho- */
/*                         26 rho0 */
/*                   27 rho+ */
/*                   28 omega meson */
/*                   29 phi */
/* lin-11/07/00: */
/*                  30 K*+ */
/* sp01/03/01 */
/*                 -14 LAMBDA(bar) */
/*                  -15 sigma-(bar) */
/*                  -16 sigma0(bar) */
/*                  -17 sigma+(bar) */
/*                   31 eta-prime */
/*                   40 cascade- */
/*                  -40 cascade-(bar) */
/*                   41 cascade0 */
/*                  -41 cascade0(bar) */
/*                   45 Omega baryon */
/*                  -45 Omega baryon(bar) */
/* sp01/03/01 end */
/* lin-5/2008: */
/*                   42 Deuteron (same in ampt.dat) */
/*                  -42 anti-Deuteron (same in ampt.dat) */

/*                   ++  ------- SEE BAO-AN LI'S NOTE BOOK */
/* ********************************* */
/* bz11/16/98 */
/*      PROGRAM ART */
/* Subroutine */ int artmn_(void)
{
    /* Initialized data */

    static real zet[91] = { 1.f,0.f,0.f,0.f,0.f,1.f,0.f,0.f,0.f,0.f,0.f,0.f,
	    0.f,0.f,0.f,-1.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,
	    -1.f,0.f,1.f,0.f,-1.f,0.f,-1.f,0.f,-2.f,-1.f,0.f,1.f,0.f,0.f,0.f,
	    0.f,-1.f,0.f,1.f,0.f,-1.f,0.f,1.f,-1.f,0.f,1.f,2.f,0.f,1.f,0.f,
	    1.f,0.f,-1.f,0.f,1.f,0.f,0.f,0.f,-1.f,0.f,1.f,0.f,-1.f,0.f,1.f,
	    0.f,0.f,1.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,-1.f,0.f,0.f,0.f,
	    0.f,-1.f };

    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6;
    real r__1, r__2, r__3, r__4;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);
    integer s_wsfe(cilist *), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    double sqrt(doublereal), log(doublereal), exp(doublereal);
    integer i_nint(real *);

    /* Local variables */
    static real b;
    static integer i__, j;
    static real s, b2;
    static integer i0, l0, ia, ld, le, ic, ie, il, jl;
    static real ct, et[150001]	/* was [150001][1] */, bx, by;
    static integer ir, lt[150001]	/* was [150001][1] */, ix, iy, iz;
    static real pt[450003]	/* was [3][150001][1] */, rt[450003]	/* 
	    was [3][150001][1] */;
    static integer is, np, nt;
    static real sp, sk;
    static integer lp, ln;
    static real tz;
    static integer ib, ii, lb1, ld1, ld2, ld3, ld4, ln1, lp1, lp2, lp3;
    static real pr0, pr1;
    static integer ln2, ln5, np1;
    static real xx1, xx2, add, ald, ale, akg;
    static integer ldd;
    static real apd, bkg, den, eta, rdd, arh, aom, akn;
    static integer lpd, kkk;
    static real app;
    static integer lkn;
    static real rpd;
    static integer npi[1];
    static real epr, rkn, asy, rpp, rpn;
    static integer lpp, lpn;
    static real apn;
    static integer iso;
    static real alp, aln, ddx, ddy, ddz, grp, spt, spz, drr, udt;
    static integer iss, nsh;
    static real ska0, ska1, aln5, addk, facl, cden;
    static integer lddk;
    static real andk;
    static integer mean, ncen;
    static real rddk, bmax;
    extern /* Subroutine */ int dens_(integer *, integer *, integer *, 
	    integer *);
    static integer lndk;
    static real ares, adou, rndk, appk, annk;
    static integer lnnk, mass;
    extern /* Subroutine */ int init_(integer *, integer *, integer *, real *,
	     real *, real *, real *, real *, integer *, integer *, integer *);
    static real temp[450003]	/* was [3][150001] */;
    static integer irun;
    static real pzta, prot[150001]	/* was [150001][1] */, rppk, rnnk, 
	    rres;
    static integer lppk;
    static real pzpr;
    static integer lrho, lres, ldou, mrun;
    static real rnsg, ecor;
    static integer idir;
    static real ekin, rads, zras, vols, ecms0, engs;
    extern /* Subroutine */ int flow_(integer *);
    static integer lrun;
    static real betac, acnnd, gammc, acndn, akaon, radta;
    static integer lbloc, lcnne, ikaon;
    static real esbin, radpr;
    extern /* Subroutine */ int gradu_(integer *, integer *, integer *, 
	    integer *, real *, real *, real *);
    static real skaon[7];
    static integer imany;
    static real gradx, grady, gradz, rcnne, rcnnd, zdist, rcndn, rcoll, rbloc,
	     rdirt;
    static integer lkaon, numnt, lcoll, lcnnd, lcndn, ldirt, lnnom;
    static real adirt, acoll, annom, rdiff;
    static integer nlost;
    static real denst, radut;
    static integer ipart;
    static real adecay, betata;
    static integer ldecay;
    static real addrho;
    extern /* Subroutine */ int tablem_(void);
    static integer lomega;
    static real gammta, rdecay, alkaon, gammas;
    static integer lddrho;
    static real betapr;
    extern /* Subroutine */ int graduk_(integer *, integer *, integer *, real 
	    *, real *, real *);
    static real sekaon[14007]	/* was [7][2001] */;
    extern integer iarflv_(integer *);
    static real gammpr, psqare;
    static integer ntotal;
    extern doublereal ranart_(integer *);
    extern integer invflv_(integer *);
    static integer outpar;
    extern /* Subroutine */ int coulin_(integer *, integer *, integer *);
    static integer lkaons;
    extern /* Subroutine */ int relcol_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, real *, real *, real *
	    );
    static integer lnnrho;
    static real annrho, sumene, etotal, atotal;
    static integer nquark, nbaryn;
    static real edenst, gradxk, gradyk, gradzk, gradxn, gradyn, gradzn, 
	    gradxp, gradyp, gradzp;
    extern /* Subroutine */ int gradup_(integer *, integer *, integer *, real 
	    *, real *, real *), gradun_(integer *, integer *, integer *, real 
	    *, real *, real *);
    static real ctlong;
    extern /* Subroutine */ int hbtout_(integer *, integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___14 = { 0, 12, 0, "(//10X,'**** FATAL ERROR: TOO MANY "
	    "TEST PART. **** ')", 0 };


/* bz11/16/98end */
/* ********************************* */
/* PARAMETERS:                                                           * */
/*  MAXPAR     - MAXIMUM NUMBER OF PARTICLES      PROGRAM CAN HANDLE     * */
/*  MAXP       - MAXIMUM NUMBER OF CREATED MESONS PROGRAM CAN HANDLE     * */
/*  MAXR       - MAXIMUM NUMBER OF EVENTS AT EACH IMPACT PARAMETER       * */
/*  MAXX       - NUMBER OF MESHPOINTS IN X AND Y DIRECTION = 2 MAXX + 1  * */
/*  MAXZ       - NUMBER OF MESHPOINTS IN Z DIRECTION       = 2 MAXZ + 1  * */
/*  AMU        - 1 ATOMIC MASS UNIT "GEV/C**2"                           * */
/*  MX,MY,MZ   - MESH SIZES IN COORDINATE SPACE [FM] FOR PAULI LATTICE   * */
/*  MPX,MPY,MPZ- MESH SIZES IN MOMENTUM SPACE [GEV/C] FOR PAULI LATTICE  * */
/* ---------------------------------------------------------------------- * */
/* lin      PARAMETER     (maxpar=200000,MAXR=50,AMU= 0.9383, */
/* lin      PARAMETER (MAXP = 14000) */
/* ----------------------------------------------------------------------* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /DD/ */
/* c      SAVE /EE/ */
/* c      SAVE /HH/ */
/* c      SAVE /ff/ */
/* c      SAVE /gg/ */
/* c      SAVE /INPUT/ */
/* c      SAVE /RR/ */
/* c      SAVE /ss/ */
/* c      SAVE /zz/ */
/* c      SAVE /RUN/ */
/* lin-4/2008: */
/*      COMMON  /KKK/     TKAON(7),EKAON(7,0:200) */
/* c      SAVE /KKK/ */
/* c      SAVE /KAON/ */
/* c      SAVE /TABLE/ */
/* c      SAVE /input1/ */
/* c      SAVE /DDpi/ */
/* c      SAVE /tt/ */
/* lin-4/2008: */
/*      DIMENSION TEMP(3,MAXSTR),SKAON(7),SEKAON(7,0:200) */
/* bz12/2/98 */
/* c      SAVE /INPUT2/ */
/* c      SAVE /INPUT3/ */
/* bz12/2/98end */
/* bz11/16/98 */
/* c      SAVE /ARPRNT/ */
/* .....note in the below, since a common block in ART is called EE, */
/* .....the variable EE in /ARPRC/is changed to PEAR. */
/* lin-9/29/03 changed name in order to distinguish from /prec2/ */
/*        COMMON /ARPRC/ ITYPAR(MAXSTR), */
/*     &       GXAR(MAXSTR), GYAR(MAXSTR), GZAR(MAXSTR), FTAR(MAXSTR), */
/*     &       PXAR(MAXSTR), PYAR(MAXSTR), PZAR(MAXSTR), PEAR(MAXSTR), */
/*     &       XMAR(MAXSTR) */
/* c      SAVE /ARPRC/ */
/* lin-9/29/03-end */
/* c      SAVE /ARERCP/ */
/* c      SAVE /ARERC1/ */
/* c      SAVE /ARPRC1/ */

/* bz11/16/98end */
/* c      SAVE /lastt/ */
/* c      SAVE /snn/ */
/* c      SAVE /hbt/ */
/* c      SAVE /resdcy/ */
/* c      SAVE /RNDF77/ */
/* lin-4/2008 zet() expanded to avoid out-of-bound errors: */
    hbt_1.nlast = 0;
    for (i__ = 1; i__ <= 150001; ++i__) {
	ftmax_1.ftsv[i__ - 1] = 0.f;
	for (irun = 1; irun <= 1; ++irun) {
	    ftmax_1.ftsvt[i__ + irun * 150001 - 150002] = 0.f;
/* L1101: */
	}
	hbt_1.lblast[i__ - 1] = 999;
	for (j = 1; j <= 4; ++j) {
/* lin-4/2008 bugs pointed out by Vander Molen & Westfall: */
/*            xlast(i,j)=0. */
/*            plast(i,j)=0. */
	    hbt_1.xlast[j + (i__ << 2) - 5] = 0.f;
	    hbt_1.plast[j + (i__ << 2) - 5] = 0.f;
/* L1001: */
	}
/* L1002: */
    }
/* -------------------------------------------------------------------* */
/* Input information about the reaction system and contral parameters* */
/* -------------------------------------------------------------------* */
/*              input section starts here                           * */
/* -------------------------------------------------------------------* */
/* bz12/2/98 */
/* .....input section is moved to subroutine ARTSET */
/* bz12/2/98end */
/* -----------------------------------------------------------------------* */
/*                   input section ends here                            * */
/* -----------------------------------------------------------------------* */
/* read in the table for gengrating the transverse momentum */
/* IN THE NN-->DDP PROCESS */
    tablem_();
/* several control parameters, keep them fixed in this code. */
    ikaon = 1;
    input_1.nstar = 1;
    input_1.ndirct = 0;
    input_1.dir = .02f;
    asy = .032f;
    esbin = .04f;
    kaon_1.mf = 36;
/* ----------------------------------------------------------------------* */
/*      CALL FRONT(12,MASSTA,MASSPR,ELAB) */
/* ----------------------------------------------------------------------* */
    d__1 = (doublereal) ((real) input1_1.massta);
    radta = pow_dd(&d__1, &c_b5) * 1.124f;
    d__1 = (doublereal) ((real) input1_1.masspr);
    radpr = pow_dd(&d__1, &c_b5) * 1.124f;
    zdist = radta + radpr;
/*      if ( cycbox.ne.0 ) zdist=0 */
    bmax = radta + radpr;
    mass = input1_1.massta + input1_1.masspr;
    ntotal = run_1.num * mass;

    if (ntotal > 150001) {
	s_wsfe(&io___14);
	e_wsfe();
	s_stop("", (ftnlen)0);
    }

/* ----------------------------------------------------------------------- */
/*       RELATIVISTIC KINEMATICS */

/*       1) LABSYSTEM */

    eta = (real) input1_1.massta * .9383f;
    pzta = 0.f;
    betata = 0.f;
    gammta = 1.f;

    epr = (real) input1_1.masspr * (input3_1.elab * .001f + .9383f);
/* Computing 2nd power */
    r__1 = epr;
/* Computing 2nd power */
    r__2 = (real) input1_1.masspr * .9383f;
    pzpr = sqrt(r__1 * r__1 - r__2 * r__2);
    betapr = pzpr / epr;
/* Computing 2nd power */
    r__1 = betapr;
    gammpr = 1.f / sqrt(1.f - r__1 * r__1);

/* BETAC AND GAMMAC OF THE C.M. OBSERVED IN THE LAB. FRAME */
    betac = (pzpr + pzta) / (epr + eta);
/* Computing 2nd power */
    r__1 = betac;
    gammc = 1.f / sqrt(1.f - r__1 * r__1);

/*      WRITE(12,'(/10x,''****    KINEMATICAL PARAMETERS    ****''/)') */
/*      WRITE(12,'(10x,''1) LAB-FRAME:        TARGET PROJECTILE'')') */
/*      WRITE(12,'(10x,''   ETOTAL "GEV" '',2F11.4)') ETA, EPR */
/*      WRITE(12,'(10x,''   P "GEV/C"    '',2F11.4)') PZTA, PZPR */
/*      WRITE(12,'(10x,''   BETA         '',2F11.4)') BETATA, BETAPR */
/*      WRITE(12,'(10x,''   GAMMA        '',2F11.4)') GAMMTA, GAMMPR */
    if (input2_1.insys != 0) {

/*       2) C.M. SYSTEM */

/* Computing 2nd power */
	r__1 = epr + eta;
/* Computing 2nd power */
	r__2 = pzpr;
	s = r__1 * r__1 - r__2 * r__2;
	xx1 = log((real) input1_1.massta) * 4.f;
	xx2 = log((real) input1_1.masspr) * 4.f;
	xx1 = exp(xx1);
	xx2 = exp(xx2);
/* Computing 2nd power */
	r__1 = s;
/* Computing 2nd power */
	i__1 = input1_1.massta;
/* Computing 2nd power */
	i__2 = input1_1.masspr;
/* Computing 2nd power */
	i__3 = input1_1.massta;
/* Computing 2nd power */
	i__4 = input1_1.masspr;
	psqare = (r__1 * r__1 + (xx1 + xx2) * .77511629195947218f - s * 2.f * 
		.88040689000000005f * (real) (i__1 * i__1 + i__2 * i__2) - (
		real) (i__3 * i__3 * (i__4 * i__4)) * 2.f * 
		.77511629195947218f) / (s * 4.f);

/* Computing 2nd power */
	r__1 = (real) input1_1.massta * .9383f;
	eta = sqrt(psqare + r__1 * r__1);
	pzta = -sqrt(psqare);
	betata = pzta / eta;
/* Computing 2nd power */
	r__1 = betata;
	gammta = 1.f / sqrt(1.f - r__1 * r__1);

/* Computing 2nd power */
	r__1 = (real) input1_1.masspr * .9383f;
	epr = sqrt(psqare + r__1 * r__1);
	pzpr = sqrt(psqare);
	betapr = pzpr / epr;
/* Computing 2nd power */
	r__1 = betapr;
	gammpr = 1.f / sqrt(1.f - r__1 * r__1);

/*        WRITE(12,'(10x,''2) C.M.-FRAME:  '')') */
/*        WRITE(12,'(10x,''   ETOTAL "GEV" '',2F11.4)') ETA, EPR */
/*        WRITE(12,'(10x,''   P "GEV/C"    '',2F11.4)') PZTA, PZPR */
/*        WRITE(12,'(10x,''   BETA         '',2F11.4)') BETATA, BETAPR */
/*        WRITE(12,'(10x,''   GAMMA        '',2F11.4)') GAMMTA, GAMMPR */
/*        WRITE(12,'(10x,''S "GEV**2"      '',F11.4)')  S */
/*        WRITE(12,'(10x,''PSQARE "GEV/C"2 '',E14.3)')  PSQARE */
/*        WRITE(12,'(/10x,''*** CALCULATION DONE IN CM-FRAME ***''/)') */
    } else {
/*        WRITE(12,'(/10x,''*** CALCULATION DONE IN LAB-FRAME ***''/)') */
    }
/* MOMENTUM PER PARTICLE */
    pzta /= (real) input1_1.massta;
    pzpr /= (real) input1_1.masspr;
/* total initial energy in the N-N cms frame */
    ecms0 = eta + epr;
/* ----------------------------------------------------------------------- */

/* Start loop over many runs of different impact parameters */
/* IF MANYB=1, RUN AT A FIXED IMPACT PARAMETER B0, OTHERWISE GENERATE */
/* MINIMUM BIAS EVENTS WITHIN THE IMPACT PARAMETER RANGE OF B_MIN AND B_MAX */
    i__1 = input2_1.manyb;
    for (imany = 1; imany <= i__1; ++imany) {
/* ------------------------------------------------------------------------ */
/* Initialize the impact parameter B */
	if (input2_1.manyb > 1) {
L111:
	    bx = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	    by = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	    b2 = bx * bx + by * by;
	    if (b2 > 1.f) {
		goto L111;
	    }
	    b = sqrt(b2) * (input3_1.bm - input3_1.bi) + input3_1.bi;
	} else {
	    b = input3_1.b0;
	}
/*      WRITE(12,'(///10X,''RUN NUMBER:'',I6)') IMANY */
/*      WRITE(12,'(//10X,''IMPACT PARAMETER B FOR THIS RUN:'', */
/*     &             F9.3,'' FM''/10X,49(''*'')/)') B */

/* ----------------------------------------------------------------------- */
/*       INITIALIZATION */
/* 1 INITIALIZATION IN ISOSPIN SPACE FOR BOTH THE PROJECTILE AND TARGET */
	coulin_(&input1_1.masspr, &input1_1.massta, &run_1.num);
/* 2 INITIALIZATION IN PHASE SPACE FOR THE TARGET */
	r__1 = b / 2.f;
	r__2 = input3_1.zeropt + zdist / 2.f;
	init_(&c__1, &input1_1.massta, &run_1.num, &radta, &r__1, &r__2, &
		pzta, &gammta, &input1_1.iseed, &mass, &input2_1.imomen);
/* 3.1 INITIALIZATION IN PHASE SPACE FOR THE PROJECTILE */
	i__2 = input1_1.massta + 1;
	r__1 = -b / 2.f;
	r__2 = input3_1.zeropt - zdist / 2.f;
	init_(&i__2, &mass, &run_1.num, &radpr, &r__1, &r__2, &pzpr, &gammpr, 
		&input1_1.iseed, &mass, &input2_1.imomen);
/* 3.2 OUTPAR IS THE NO. OF ESCAPED PARTICLES */
	outpar = 0;
/* 3.3 INITIALIZATION FOR THE NO. OF PARTICLES IN EACH SAMPLE */
/*    THIS IS NEEDED DUE TO THE FACT THAT PIONS CAN BE PRODUCED OR ABSORBED */
	rr_1.massr[0] = 0;
	i__2 = run_1.num;
	for (ir = 1; ir <= i__2; ++ir) {
	    rr_1.massr[ir] = mass;
/* L1003: */
	}
/* 3.4 INITIALIZation FOR THE KAON SPECTRUM */
/*      CALL KSPEC0(BETAC,GAMMC) */
/* calculate the local baryon density matrix */
	dens_(&input2_1.ipot, &mass, &run_1.num, &outpar);

/* ----------------------------------------------------------------------- */
/*       CONTROL PRINTOUT OF INITIAL CONFIGURATION */

/*      WRITE(12,'(''**********  INITIAL CONFIGURATION  **********''/)') */

/* print out the INITIAL density matrix in the reaction plane */
/*       do ix=-10,10 */
/*       do iz=-10,10 */
/*       write(1053,992)ix,iz,rho(ix,0,iz)/0.168 */
/*       end do */
/*       end do */
/* ----------------------------------------------------------------------- */
/*       CALCULATE MOMENTA FOR T = 0.5 * DT */
/*       (TO OBTAIN 2ND DEGREE ACCURACY!) */
/*       "Reference: J. AICHELIN ET AL., PHYS. REV. C31, 1730 (1985)" */

	if (input2_1.icoll != -1) {
	    i__2 = ntotal;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		ix = i_nint(&aa_1.r__[i__ * 3 - 3]);
		iy = i_nint(&aa_1.r__[i__ * 3 - 2]);
		iz = i_nint(&aa_1.r__[i__ * 3 - 1]);
/* lin-4/2008 check bounds: */
		if (ix >= 20 || iy >= 20 || iz >= 24 || ix <= -20 || iy <= 
			-20 || iz <= -24) {
		    goto L700;
		}
		gradu_(&input2_1.ipot, &ix, &iy, &iz, &gradx, &grady, &gradz);
		bb_1.p[i__ * 3 - 3] -= input1_1.dt * .5f * gradx;
		bb_1.p[i__ * 3 - 2] -= input1_1.dt * .5f * grady;
		bb_1.p[i__ * 3 - 1] -= input1_1.dt * .5f * gradz;
L700:
		;
	    }
	}
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* 4 INITIALIZATION OF TIME-LOOP VARIABLES */
/* 4.1 COLLISION NUMBER COUNTERS */
/* lin 51      RCNNE  = 0 */
	rcnne = 0.f;
	rdd = 0.f;
	rpp = 0.f;
	rppk = 0.f;
	rpn = 0.f;
	rpd = 0.f;
	rkn = 0.f;
	rnnk = 0.f;
	rddk = 0.f;
	rndk = 0.f;
	rcnnd = 0.f;
	rcndn = 0.f;
	rcoll = 0.f;
	rbloc = 0.f;
	rdirt = 0.f;
	rdecay = 0.f;
	rres = 0.f;
/* 4.11 KAON PRODUCTION PROBABILITY COUNTER FOR PERTURBATIVE CALCULATIONS ONLY */
	for (kkk = 1; kkk <= 5; ++kkk) {
	    skaon[kkk - 1] = 0.f;
	    for (is = 1; is <= 2000; ++is) {
		sekaon[kkk + is * 7 - 1] = 0.f;
/* L1004: */
	    }
/* L1005: */
	}
/* 4.12 anti-proton and anti-kaon counters */
	pr0 = 0.f;
	pr1 = 0.f;
	ska0 = 0.f;
	ska1 = 0.f;
/*       ============== LOOP OVER ALL TIME STEPS ================       * */
/*                             STARTS HERE                              * */
/*       ========================================================       * */
/* bz11/16/98 */
	if (arprnt_1.iapar2[0] != 1) {
	    for (i__ = 1; i__ <= 150001; ++i__) {
		for (j = 1; j <= 3; ++j) {
		    aa_1.r__[j + i__ * 3 - 4] = 0.f;
		    bb_1.p[j + i__ * 3 - 4] = 0.f;
/* L1015: */
		}
		cc_1.e[i__ - 1] = 0.f;
		ee_1.lb[i__ - 1] = 0;
/* bz3/25/00 */
		ee_1.id[i__ - 1] = 0;
/*     sp 12/19/00 */
		hh_1.proper[i__ - 1] = 1.f;
/* L1016: */
	    }
	    mass = 0;
/* bz12/22/98 */
/*         MASSR(1) = 0 */
/*         NP = 0 */
/*         NPI = 1 */
	    np = 0;
	    i__2 = run_1.num;
	    for (j = 1; j <= i__2; ++j) {
		rr_1.massr[j] = 0;
		npi[j - 1] = 1;
/* L1017: */
	    }
	    for (i__ = 1; i__ <= 1; ++i__) {
		for (j = 1; j <= 150001; ++j) {
		    rt[(j + i__ * 150001) * 3 - 450006] = 0.f;
		    rt[(j + i__ * 150001) * 3 - 450005] = 0.f;
		    rt[(j + i__ * 150001) * 3 - 450004] = 0.f;
		    pt[(j + i__ * 150001) * 3 - 450006] = 0.f;
		    pt[(j + i__ * 150001) * 3 - 450005] = 0.f;
		    pt[(j + i__ * 150001) * 3 - 450004] = 0.f;
		    et[j + i__ * 150001 - 150002] = 0.f;
		    lt[j + i__ * 150001 - 150002] = 0;
/*     sp 12/19/00 */
		    prot[j + i__ * 150001 - 150002] = 1.f;
/* L1018: */
		}
/* L1019: */
	    }
/* bz12/22/98end */
	}
/* bz11/16/98end */
	i__2 = input2_1.ntmax;
	for (nt = 1; nt <= i__2; ++nt) {
/* TEMPORARY PARTICLE COUNTERS */
/* 4.2 PION COUNTERS : LP1,LP2 AND LP3 ARE THE NO. OF P+,P0 AND P- */
	    lp1 = 0;
	    lp2 = 0;
	    lp3 = 0;
/* 4.3 DELTA COUNTERS : LD1,LD2,LD3 AND LD4 ARE THE NO. OF D++,D+,D0 AND D- */
	    ld1 = 0;
	    ld2 = 0;
	    ld3 = 0;
	    ld4 = 0;
/* 4.4 N*(1440) COUNTERS : LN1 AND LN2 ARE THE NO. OF N*+ AND N*0 */
	    ln1 = 0;
	    ln2 = 0;
/* 4.5 N*(1535) counters */
	    ln5 = 0;
/* 4.6 ETA COUNTERS */
	    le = 0;
/* 4.7 KAON COUNTERS */
	    lkaon = 0;
/* lin-11/09/00: */
/* KAON* COUNTERS */
	    lkaons = 0;
/* ----------------------------------------------------------------------- */
	    if (input2_1.icoll != 1) {
/* STUDYING BINARY COLLISIONS AMONG PARTICLES DURING THIS TIME INTERVAL * */
/* lin-10/25/02 get rid of argument usage mismatch in relcol(.nt.): */
		numnt = nt;
		relcol_(&lcoll, &lbloc, &lcnne, &ldd, &lpp, &lppk, &lpn, &lpd,
			 &lrho, &lomega, &lkn, &lnnk, &lddk, &lndk, &lcnnd, &
			lcndn, &ldirt, &ldecay, &lres, &ldou, &lddrho, &
			lnnrho, &lnnom, &numnt, &input2_1.ntmax, &sp, &akaon, 
			&sk);
/*     &    LNNOM,NT,ntmax,sp,akaon,sk) */
/* lin-10/25/02-end */
/* ----------------------------------------------------------------------- */
/* dilepton production from Dalitz decay */
/* of pi0 at final time */
/*      if(nt .eq. ntmax) call dalitz_pi(nt,ntmax) */
/*                                                                      * */
/* ********************************* */
/*                Lables of collision channels                             * */
/* ********************************* */
/*         LCOLL   - NUMBER OF COLLISIONS              (INTEGER,OUTPUT) * */
/*         LBLOC   - NUMBER OF PULI-BLOCKED COLLISIONS (INTEGER,OUTPUT) * */
/*         LCNNE   - NUMBER OF ELASTIC COLLISION       (INTEGER,OUTPUT) * */
/*         LCNND   - NUMBER OF N+N->N+DELTA REACTION   (INTEGER,OUTPUT) * */
/*         LCNDN   - NUMBER OF N+DELTA->N+N REACTION   (INTEGER,OUTPUT) * */
/*         LDD     - NUMBER OF RESONANCE+RESONANCE COLLISIONS */
/*         LPP     - NUMBER OF PION+PION elastic COLIISIONS */
/*         lppk    - number of pion(RHO,OMEGA)+pion(RHO,OMEGA) */
/*                 -->K+K- collisions */
/*         LPN     - NUMBER OF PION+N-->KAON+X */
/*         lpd     - number of pion+n-->delta+pion */
/*         lrho    - number of pion+n-->Delta+rho */
/*         lomega  - number of pion+n-->Delta+omega */
/*         LKN     - NUMBER OF KAON RESCATTERINGS */
/*         LNNK    - NUMBER OF bb-->kAON PROCESS */
/*         LDDK    - NUMBER OF DD-->KAON PROCESS */
/*         LNDK    - NUMBER OF ND-->KAON PROCESS */
/* ********************************** */
/* TIME-INTEGRATED COLLISIONS NUMBERS OF VARIOUS PROCESSES */
		rcoll += (real) lcoll / run_1.num;
		rbloc += (real) lbloc / run_1.num;
		rcnne += (real) lcnne / run_1.num;
		rdd += (real) ldd / run_1.num;
		rpp += (real) lpp / run_1.num;
		rppk += (real) lppk / run_1.num;
		rpn += (real) lpn / run_1.num;
		rpd += (real) lpd / run_1.num;
		rkn += (real) lkn / run_1.num;
		rnnk += (real) lnnk / run_1.num;
		rddk += (real) lddk / run_1.num;
		rndk += (real) lndk / run_1.num;
		rcnnd += (real) lcnnd / run_1.num;
		rcndn += (real) lcndn / run_1.num;
		rdirt += (real) ldirt / run_1.num;
		rdecay += (real) ldecay / run_1.num;
		rres += (real) lres / run_1.num;
/* AVERAGE RATES OF VARIOUS COLLISIONS IN THE CURRENT TIME STEP */
		adirt = ldirt / input1_1.dt / run_1.num;
		acoll = (lcoll - lbloc) / input1_1.dt / run_1.num;
		acnnd = lcnnd / input1_1.dt / run_1.num;
		acndn = lcndn / input1_1.dt / run_1.num;
		adecay = ldecay / input1_1.dt / run_1.num;
		ares = lres / input1_1.dt / run_1.num;
		adou = ldou / input1_1.dt / run_1.num;
		addrho = lddrho / input1_1.dt / run_1.num;
		annrho = lnnrho / input1_1.dt / run_1.num;
		annom = lnnom / input1_1.dt / run_1.num;
		add = ldd / input1_1.dt / run_1.num;
		app = lpp / input1_1.dt / run_1.num;
		appk = lppk / input1_1.dt / run_1.num;
		apn = lpn / input1_1.dt / run_1.num;
		apd = lpd / input1_1.dt / run_1.num;
		arh = lrho / input1_1.dt / run_1.num;
		aom = lomega / input1_1.dt / run_1.num;
		akn = lkn / input1_1.dt / run_1.num;
		annk = lnnk / input1_1.dt / run_1.num;
		addk = lddk / input1_1.dt / run_1.num;
		andk = lndk / input1_1.dt / run_1.num;
/* PRINT OUT THE VARIOUS COLLISION RATES */
/* (1)N-N COLLISIONS */
/*       WRITE(1010,9991)NT*DT,ACNND,ADOU,ADIRT,ADDRHO,ANNRHO+ANNOM */
/* 9991       FORMAT(6(E10.3,2X)) */
/* (2)PION-N COLLISIONS */
/*       WRITE(1011,'(5(E10.3,2X))')NT*DT,apd,ARH,AOM,APN */
/* (3)KAON PRODUCTION CHANNELS */
/*        WRITE(1012,9993)NT*DT,ANNK,ADDK,ANDK,APN,Appk */
/* (4)D(N*)+D(N*) COLLISION */
/*       WRITE(1013,'(4(E10.3,2X))')NT*DT,ADDK,ADD,ADD+ADDK */
/* (5)MESON+MESON */
/*       WRITE(1014,'(4(E10.3,2X))')NT*DT,APPK,APP,APP+APPK */
/* (6)DECAY AND RESONANCE */
/*       WRITE(1016,'(3(E10.3,2X))')NT*DT,ARES,ADECAY */
/* (7)N+D(N*) */
/*       WRITE(1017,'(4(E10.3,2X))')NT*DT,ACNDN,ANDK,ACNDN+ANDK */
/* 9992    FORMAT(5(E10.3,2X)) */
/* 9993    FORMAT(6(E10.3,2X)) */
/* PRINT OUT TIME-INTEGRATED COLLISION INFORMATION */
/* bz12/28/98 */
/*        write(1018,'(5(e10.3,2x),/, 4(e10.3,2x))') */
/*     &           RCNNE,RCNND,RCNDN,RDIRT,rpd, */
/*     &           RDECAY,RRES,RDD,RPP */
/*        write(1018,'(6(e10.3,2x),/, 5(e10.3,2x))') */
/*     &           NT*DT,RCNNE,RCNND,RCNDN,RDIRT,rpd, */
/*     &           NT*DT,RDECAY,RRES,RDD,RPP */
/* bz12/18/98end */
/* PRINT OUT TIME-INTEGRATED KAON MULTIPLICITIES FROM DIFFERENT CHANNELS */
/*       WRITE(1019,'(7(E10.3,2X))')NT*DT,RNNK,RDDK,RNDK,RPN,Rppk, */
/*     &                           RNNK+RDDK+RNDK+RPN+Rppk */
/*                                                                      * */
	    }

/*       UPDATE BARYON DENSITY */

	    dens_(&input2_1.ipot, &mass, &run_1.num, &outpar);

/*       UPDATE POSITIONS FOR ALL THE PARTICLES PRESENT AT THIS TIME */

	    sumene = 0.f;
	    iso = 0;
	    i__3 = run_1.num;
	    for (mrun = 1; mrun <= i__3; ++mrun) {
		iso += rr_1.massr[mrun - 1];
		i__4 = rr_1.massr[mrun];
		for (i0 = 1; i0 <= i__4; ++i0) {
		    i__ = i0 + iso;
/* Computing 2nd power */
		    r__1 = cc_1.e[i__ - 1];
/* Computing 2nd power */
		    r__2 = bb_1.p[i__ * 3 - 3];
/* Computing 2nd power */
		    r__3 = bb_1.p[i__ * 3 - 2];
/* Computing 2nd power */
		    r__4 = bb_1.p[i__ * 3 - 1];
		    etotal = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + 
			    r__4 * r__4);
		    sumene += etotal;
/* for kaons, if there is a potential */
/* CALCULATE THE ENERGY OF THE KAON ACCORDING TO THE IMPULSE APPROXIMATION */
/* REFERENCE: B.A. LI AND C.M. KO, PHYS. REV. C 54 (1996) 3283. */
		    if (input2_1.kpoten != 0 && ee_1.lb[i__ - 1] == 23) {
			den = 0.f;
			ix = i_nint(&aa_1.r__[i__ * 3 - 3]);
			iy = i_nint(&aa_1.r__[i__ * 3 - 2]);
			iz = i_nint(&aa_1.r__[i__ * 3 - 1]);
/* lin-4/2008: */
/*       IF (ABS(IX) .LT. MAXX .AND. ABS(IY) .LT. MAXX .AND. */
/*     & ABS(IZ) .LT. MAXZ) den=rho(ix,iy,iz) */
			if (ix < 20 && iy < 20 && iz < 24 && ix > -20 && iy > 
				-20 && iz > -24) {
			    den = dd_1.rho[ix + (iy + iz * 41) * 41 + 41184];
			}
/*         ecor=0.1973**2*0.255*kmul*4*3.14159*(1.+0.4396/0.938) */
/*         etotal=sqrt(etotal**2+ecor*den) */
/* ** G.Q Li potential form with n_s = n_b and pot(n_0)=29 MeV, m^*=m */
/*     GeV^2 fm^3 */
			akg = .1727f;
/*     GeV fm^3 */
			bkg = .333f;
			rnsg = den;
/* Computing 2nd power */
			r__1 = bkg * den;
			ecor = -akg * rnsg + r__1 * r__1;
/* Computing 2nd power */
			r__1 = etotal;
			etotal = sqrt(r__1 * r__1 + ecor);
		    }

		    if (input2_1.kpoten != 0 && ee_1.lb[i__ - 1] == 21) {
			den = 0.f;
			ix = i_nint(&aa_1.r__[i__ * 3 - 3]);
			iy = i_nint(&aa_1.r__[i__ * 3 - 2]);
			iz = i_nint(&aa_1.r__[i__ * 3 - 1]);
/* lin-4/2008: */
/*       IF (ABS(IX) .LT. MAXX .AND. ABS(IY) .LT. MAXX .AND. */
/*     & ABS(IZ) .LT. MAXZ) den=rho(ix,iy,iz) */
			if (ix < 20 && iy < 20 && iz < 24 && ix > -20 && iy > 
				-20 && iz > -24) {
			    den = dd_1.rho[ix + (iy + iz * 41) * 41 + 41184];
			}
/* * for song potential no effect on position */
/* ** G.Q Li potential form with n_s = n_b and pot(n_0)=29 MeV, m^*=m */
/*     GeV^2 fm^3 */
			akg = .1727f;
/*     GeV fm^3 */
			bkg = .333f;
			rnsg = den;
/* Computing 2nd power */
			r__1 = bkg * den;
			ecor = -akg * rnsg + r__1 * r__1;
/* Computing 2nd power */
			r__1 = etotal;
			etotal = sqrt(r__1 * r__1 + ecor);
		    }

/* UPDATE POSITIONS */
		    aa_1.r__[i__ * 3 - 3] += input1_1.dt * bb_1.p[i__ * 3 - 3]
			     / etotal;
		    aa_1.r__[i__ * 3 - 2] += input1_1.dt * bb_1.p[i__ * 3 - 2]
			     / etotal;
		    aa_1.r__[i__ * 3 - 1] += input1_1.dt * bb_1.p[i__ * 3 - 1]
			     / etotal;
/* use cyclic boundary conitions */
		    if (input3_1.cycbox != 0.f) {
			if (aa_1.r__[i__ * 3 - 3] > input3_1.cycbox / 2) {
			    aa_1.r__[i__ * 3 - 3] -= input3_1.cycbox;
			}
			if (aa_1.r__[i__ * 3 - 3] <= -input3_1.cycbox / 2) {
			    aa_1.r__[i__ * 3 - 3] += input3_1.cycbox;
			}
			if (aa_1.r__[i__ * 3 - 2] > input3_1.cycbox / 2) {
			    aa_1.r__[i__ * 3 - 2] -= input3_1.cycbox;
			}
			if (aa_1.r__[i__ * 3 - 2] <= -input3_1.cycbox / 2) {
			    aa_1.r__[i__ * 3 - 2] += input3_1.cycbox;
			}
			if (aa_1.r__[i__ * 3 - 1] > input3_1.cycbox / 2) {
			    aa_1.r__[i__ * 3 - 1] -= input3_1.cycbox;
			}
			if (aa_1.r__[i__ * 3 - 1] <= -input3_1.cycbox / 2) {
			    aa_1.r__[i__ * 3 - 1] += input3_1.cycbox;
			}
		    }
/* UPDATE THE DELTA, N* AND PION COUNTERS */
		    lb1 = ee_1.lb[i__ - 1];
/* 1. FOR DELTA++ */
		    if (lb1 == 9) {
			++ld1;
		    }
/* 2. FOR DELTA+ */
		    if (lb1 == 8) {
			++ld2;
		    }
/* 3. FOR DELTA0 */
		    if (lb1 == 7) {
			++ld3;
		    }
/* 4. FOR DELTA- */
		    if (lb1 == 6) {
			++ld4;
		    }
/* 5. FOR N*+(1440) */
		    if (lb1 == 11) {
			++ln1;
		    }
/* 6. FOR N*0(1440) */
		    if (lb1 == 10) {
			++ln2;
		    }
/* 6.1 FOR N*(1535) */
		    if (lb1 == 13 || lb1 == 12) {
			++ln5;
		    }
/* 6.2 FOR ETA */
		    if (lb1 == 0) {
			++le;
		    }
/* 6.3 FOR KAONS */
		    if (lb1 == 23) {
			++lkaon;
		    }
/* lin-11/09/00: FOR KAON* */
		    if (lb1 == 30) {
			++lkaons;
		    }
/* UPDATE PION COUNTER */
/* 7. FOR PION+ */
		    if (lb1 == 5) {
			++lp1;
		    }
/* 8. FOR PION0 */
		    if (lb1 == 4) {
			++lp2;
		    }
/* 9. FOR PION- */
		    if (lb1 == 3) {
			++lp3;
		    }
/* L201: */
		}
	    }
	    lp = lp1 + lp2 + lp3;
	    ld = ld1 + ld2 + ld3 + ld4;
	    ln = ln1 + ln2;
	    alp = (real) lp / (real) run_1.num;
	    ald = (real) ld / (real) run_1.num;
	    aln = (real) ln / (real) run_1.num;
	    aln5 = (real) ln5 / (real) run_1.num;
	    atotal = alp + ald + aln + aln5 * .5f;
	    ale = (real) le / (real) run_1.num;
	    alkaon = (real) lkaon / (real) run_1.num;
/* UPDATE MOMENTUM DUE TO COULOMB INTERACTION */
	    if (input2_1.icou == 1) {
/*       with Coulomb interaction */
		iso = 0;
		i__4 = run_1.num;
		for (irun = 1; irun <= i__4; ++irun) {
		    iso += rr_1.massr[irun - 1];
		    i__3 = rr_1.massr[irun];
		    for (il = 1; il <= i__3; ++il) {
			temp[il * 3 - 3] = 0.f;
			temp[il * 3 - 2] = 0.f;
			temp[il * 3 - 1] = 0.f;
/* L1021: */
		    }
		    i__3 = rr_1.massr[irun];
		    for (il = 1; il <= i__3; ++il) {
			i__ = iso + il;
			if (zet[ee_1.lb[i__ - 1] + 45] != 0.f) {
			    i__5 = il - 1;
			    for (jl = 1; jl <= i__5; ++jl) {
				j = iso + jl;
				if (zet[ee_1.lb[j - 1] + 45] != 0.f) {
				    ddx = aa_1.r__[i__ * 3 - 3] - aa_1.r__[j *
					     3 - 3];
				    ddy = aa_1.r__[i__ * 3 - 2] - aa_1.r__[j *
					     3 - 2];
				    ddz = aa_1.r__[i__ * 3 - 1] - aa_1.r__[j *
					     3 - 1];
/* Computing 2nd power */
				    r__1 = ddx;
/* Computing 2nd power */
				    r__2 = ddy;
/* Computing 2nd power */
				    r__3 = ddz;
				    rdiff = sqrt(r__1 * r__1 + r__2 * r__2 + 
					    r__3 * r__3);
				    if (rdiff <= 1.f) {
					rdiff = 1.f;
				    }
/* Computing 3rd power */
				    r__1 = rdiff;
				    grp = zet[ee_1.lb[i__ - 1] + 45] * zet[
					    ee_1.lb[j - 1] + 45] / (r__1 * (
					    r__1 * r__1));
				    ddx *= grp;
				    ddy *= grp;
				    ddz *= grp;
				    temp[il * 3 - 3] += ddx;
				    temp[il * 3 - 2] += ddy;
				    temp[il * 3 - 1] += ddz;
				    temp[jl * 3 - 3] -= ddx;
				    temp[jl * 3 - 2] -= ddy;
				    temp[jl * 3 - 1] -= ddz;
				}
/* L1022: */
			    }
			}
/* L1023: */
		    }
		    i__3 = rr_1.massr[irun];
		    for (il = 1; il <= i__3; ++il) {
			i__ = iso + il;
			if (zet[ee_1.lb[i__ - 1] + 45] != 0.f) {
			    for (idir = 1; idir <= 3; ++idir) {
				bb_1.p[idir + i__ * 3 - 4] += temp[idir + il *
					 3 - 4] * input1_1.dt * .00144f;
/* L1024: */
			    }
			}
/* L1025: */
		    }
/* L1026: */
		}
	    }
/*       In the following, we shall: */
/*       (1) UPDATE MOMENTA DUE TO THE MEAN FIELD FOR BARYONS AND KAONS, */
/*       (2) calculate the thermalization, temperature in a sphere of */
/*           radius 2.0 fm AROUND THE CM */
/*       (3) AND CALCULATE THE NUMBER OF PARTICLES IN THE HIGH DENSITY REGION */
	    spt = 0.f;
	    spz = 0.f;
	    ncen = 0;
	    ekin = 0.f;
	    nlost = 0;
	    mean = 0;
	    nquark = 0;
	    nbaryn = 0;
/* sp06/18/01 */
	    rads = 2.f;
	    zras = .1f;
	    denst = 0.f;
	    edenst = 0.f;
/* sp06/18/01 end */
	    i__4 = run_1.num;
	    for (irun = 1; irun <= i__4; ++irun) {
		mean += rr_1.massr[irun - 1];
		i__3 = rr_1.massr[irun];
		for (j = 1; j <= i__3; ++j) {
		    i__ = j + mean;

/* sp06/18/01 */
/* Computing 2nd power */
		    r__1 = aa_1.r__[i__ * 3 - 3];
/* Computing 2nd power */
		    r__2 = aa_1.r__[i__ * 3 - 2];
		    radut = sqrt(r__1 * r__1 + r__2 * r__2);
		    if (radut <= rads) {
			if ((r__1 = aa_1.r__[i__ * 3 - 1], dabs(r__1)) <= 
				zras * nt * input1_1.dt) {
/*         vols = 3.14159*radut**2*abs(r(3,i))      ! cylinder pi*r^2*l */
/*     cylinder pi*r^2*l */
/* Computing 2nd power */
			    r__1 = rads;
			    vols = r__1 * r__1 * 3.14159f * zras;
/* Computing 2nd power */
			    r__1 = bb_1.p[i__ * 3 - 3];
/* Computing 2nd power */
			    r__2 = bb_1.p[i__ * 3 - 2];
/* Computing 2nd power */
			    r__3 = bb_1.p[i__ * 3 - 1];
/* Computing 2nd power */
			    r__4 = cc_1.e[i__ - 1];
			    engs = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * 
				    r__3 + r__4 * r__4);
			    gammas = 1.f;
			    if (cc_1.e[i__ - 1] != 0.f) {
				gammas = engs / cc_1.e[i__ - 1];
			    }
/*     rho */
			    denst += 1.f / gammas / vols;
/*     energy density */
			    edenst += engs / gammas / gammas / vols;
			}
		    }
/* sp06/18/01 end */

/* Computing 2nd power */
		    r__1 = aa_1.r__[i__ * 3 - 3];
/* Computing 2nd power */
		    r__2 = aa_1.r__[i__ * 3 - 2];
/* Computing 2nd power */
		    r__3 = aa_1.r__[i__ * 3 - 1];
		    drr = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
		    if (drr <= 2.f) {
/* Computing 2nd power */
			r__1 = bb_1.p[i__ * 3 - 3];
/* Computing 2nd power */
			r__2 = bb_1.p[i__ * 3 - 2];
			spt = spt + r__1 * r__1 + r__2 * r__2;
/* Computing 2nd power */
			r__1 = bb_1.p[i__ * 3 - 1];
			spz += r__1 * r__1;
			++ncen;
/* Computing 2nd power */
			r__1 = bb_1.p[i__ * 3 - 3];
/* Computing 2nd power */
			r__2 = bb_1.p[i__ * 3 - 2];
/* Computing 2nd power */
			r__3 = bb_1.p[i__ * 3 - 1];
/* Computing 2nd power */
			r__4 = cc_1.e[i__ - 1];
			ekin = ekin + sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * 
				r__3 + r__4 * r__4) - cc_1.e[i__ - 1];
		    }
		    ix = i_nint(&aa_1.r__[i__ * 3 - 3]);
		    iy = i_nint(&aa_1.r__[i__ * 3 - 2]);
		    iz = i_nint(&aa_1.r__[i__ * 3 - 1]);
/* calculate the No. of particles in the high density region */
/* lin-4/2008: */
/*              IF (ABS(IX) .LT. MAXX .AND. ABS(IY) .LT. MAXX .AND. */
/*     & ABS(IZ) .LT. MAXZ) THEN */
		    if (ix < 20 && iy < 20 && iz < 24 && ix > -20 && iy > -20 
			    && iz > -24) {
			if (dd_1.rho[ix + (iy + iz * 41) * 41 + 41184] / 
				.168f > input3_1.dencut) {
			    goto L5800;
			}
			if (dd_1.rho[ix + (iy + iz * 41) * 41 + 41184] / 
				.168f > 5.f && cc_1.e[i__ - 1] > .9f) {
			    ++nbaryn;
			}
			if (tt_1.pel[ix + (iy + iz * 41) * 41 + 41184] > 2.f) 
				{
			    ++nquark;
			}
		    }
/* * */
/* If there is a kaon potential, propogating kaons */
		    if (input2_1.kpoten != 0 && ee_1.lb[i__ - 1] == 23) {
			den = 0.f;
/* lin-4/2008: */
/*       IF (ABS(IX) .LT. MAXX .AND. ABS(IY) .LT. MAXX .AND. */
/*     & ABS(IZ) .LT. MAXZ)then */
			if (ix < 20 && iy < 20 && iz < 24 && ix > -20 && iy > 
				-20 && iz > -24) {
			    den = dd_1.rho[ix + (iy + iz * 41) * 41 + 41184];
/*        ecor=0.1973**2*0.255*kmul*4*3.14159*(1.+0.4396/0.938) */
/*       etotal=sqrt(P(1,i)**2+p(2,I)**2+p(3,i)**2+e(i)**2+ecor*den) */
/* ** for G.Q Li potential form with n_s = n_b and pot(n_0)=29 MeV */
/*     !! GeV^2 fm^3 */
			    akg = .1727f;
/*     !! GeV fm^3 */
			    bkg = .333f;
			    rnsg = den;
/* Computing 2nd power */
			    r__1 = bkg * den;
			    ecor = -akg * rnsg + r__1 * r__1;
/* Computing 2nd power */
			    r__1 = bb_1.p[i__ * 3 - 3];
/* Computing 2nd power */
			    r__2 = bb_1.p[i__ * 3 - 2];
/* Computing 2nd power */
			    r__3 = bb_1.p[i__ * 3 - 1];
/* Computing 2nd power */
			    r__4 = cc_1.e[i__ - 1];
			    etotal = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * 
				    r__3 + r__4 * r__4 + ecor);
/* Computing 2nd power */
			    r__1 = bkg;
			    ecor = -akg + r__1 * r__1 * 2.f * den + bkg * 2.f 
				    * etotal;
/* ** G.Q. Li potential (END) */
			    graduk_(&ix, &iy, &iz, &gradxk, &gradyk, &gradzk);
			    bb_1.p[i__ * 3 - 3] -= input1_1.dt * gradxk * 
				    ecor / (etotal * 2.f);
			    bb_1.p[i__ * 3 - 2] -= input1_1.dt * gradyk * 
				    ecor / (etotal * 2.f);
			    bb_1.p[i__ * 3 - 1] -= input1_1.dt * gradzk * 
				    ecor / (etotal * 2.f);
			}
		    }

		    if (input2_1.kpoten != 0 && ee_1.lb[i__ - 1] == 21) {
			den = 0.f;
/* lin-4/2008: */
/*           IF (ABS(IX) .LT. MAXX .AND. ABS(IY) .LT. MAXX .AND. */
/*     &        ABS(IZ) .LT. MAXZ)then */
			if (ix < 20 && iy < 20 && iz < 24 && ix > -20 && iy > 
				-20 && iz > -24) {
			    den = dd_1.rho[ix + (iy + iz * 41) * 41 + 41184];
			    graduk_(&ix, &iy, &iz, &gradxk, &gradyk, &gradzk);
/*        P(1,I) = P(1,I) - DT * GRADXk*(-0.12/0.168)    !! song potential */
/*        P(2,I) = P(2,I) - DT * GRADYk*(-0.12/0.168) */
/*        P(3,I) = P(3,I) - DT * GRADZk*(-0.12/0.168) */
/* ** for G.Q Li potential form with n_s = n_b and pot(n_0)=29 MeV */
/*    !! GeV^2 fm^3 */
			    akg = .1727f;
/*     !! GeV fm^3 */
			    bkg = .333f;
			    rnsg = den;
/* Computing 2nd power */
			    r__1 = bkg * den;
			    ecor = -akg * rnsg + r__1 * r__1;
/* Computing 2nd power */
			    r__1 = bb_1.p[i__ * 3 - 3];
/* Computing 2nd power */
			    r__2 = bb_1.p[i__ * 3 - 2];
/* Computing 2nd power */
			    r__3 = bb_1.p[i__ * 3 - 1];
/* Computing 2nd power */
			    r__4 = cc_1.e[i__ - 1];
			    etotal = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * 
				    r__3 + r__4 * r__4 + ecor);
/* Computing 2nd power */
			    r__1 = bkg;
			    ecor = -akg + r__1 * r__1 * 2.f * den - bkg * 2.f 
				    * etotal;
			    bb_1.p[i__ * 3 - 3] -= input1_1.dt * gradxk * 
				    ecor / (etotal * 2.f);
			    bb_1.p[i__ * 3 - 2] -= input1_1.dt * gradyk * 
				    ecor / (etotal * 2.f);
			    bb_1.p[i__ * 3 - 1] -= input1_1.dt * gradzk * 
				    ecor / (etotal * 2.f);
/* ** G.Q. Li potential (END) */
			}
		    }

/* for other mesons, there is no potential */
		    if (j > mass) {
			goto L5800;
		    }
/*  with mean field interaction for baryons   (open endif below) !!sp05 */
/* *      if( (iabs(lb(i)).eq.1.or.iabs(lb(i)).eq.2) .or. */
/* *    &     (iabs(lb(i)).ge.6.and.iabs(lb(i)).le.17) .or. */
/* *    &      iabs(lb(i)).eq.40.or.iabs(lb(i)).eq.41 )then */
		    if (input2_1.icoll != -1) {
/* check if the baryon has run off the lattice */
/*             IX0=NINT(R(1,I)/DX) */
/*             IY0=NINT(R(2,I)/DY) */
/*             IZ0=NINT(R(3,I)/DZ) */
/*             IPX0=NINT(P(1,I)/DPX) */
/*             IPY0=NINT(P(2,I)/DPY) */
/*             IPZ0=NINT(P(3,I)/DPZ) */
/*      if ( (abs(ix0).gt.mx) .or. (abs(iy0).gt.my) .or. (abs(iz0).gt.mz) */
/*     &  .or. (abs(ipx0).gt.mpx) .or. (abs(ipy0) */
/*     &  .or. (ipz0.lt.-mpz) .or. (ipz0.gt.mpzp)) NLOST=NLOST+1 */
/* lin-4/2008: */
/*              IF (ABS(IX) .LT. MAXX .AND. ABS(IY) .LT. MAXX .AND. */
/*     &                                    ABS(IZ) .LT. MAXZ     ) THEN */
			if (ix < 20 && iy < 20 && iz < 24 && ix > -20 && iy > 
				-20 && iz > -24) {
			    gradu_(&input2_1.ipot, &ix, &iy, &iz, &gradx, &
				    grady, &gradz);
			    tz = 0.f;
			    gradxn = 0.f;
			    gradyn = 0.f;
			    gradzn = 0.f;
			    gradxp = 0.f;
			    gradyp = 0.f;
			    gradzp = 0.f;
			    if (input2_1.icou == 1) {
				gradup_(&ix, &iy, &iz, &gradxp, &gradyp, &
					gradzp);
				gradun_(&ix, &iy, &iz, &gradxn, &gradyn, &
					gradzn);
				if (zet[ee_1.lb[i__ - 1] + 45] != 0.f) {
				    tz = -1.f;
				}
				if (zet[ee_1.lb[i__ - 1] + 45] == 0.f) {
				    tz = 1.f;
				}
			    }
			    if ((i__5 = ee_1.lb[i__ - 1], abs(i__5)) >= 14 && 
				    (i__6 = ee_1.lb[i__ - 1], abs(i__6)) <= 
				    17) {
				facl = .66666666666666663f;
			    } else if ((i__5 = ee_1.lb[i__ - 1], abs(i__5)) ==
				     40 || (i__6 = ee_1.lb[i__ - 1], abs(i__6)
				    ) == 41) {
				facl = .33333333333333331f;
			    } else {
				facl = 1.f;
			    }
			    bb_1.p[i__ * 3 - 3] -= facl * input1_1.dt * (
				    gradx + asy * (gradxn - gradxp) * tz);
			    bb_1.p[i__ * 3 - 2] -= facl * input1_1.dt * (
				    grady + asy * (gradyn - gradyp) * tz);
			    bb_1.p[i__ * 3 - 1] -= facl * input1_1.dt * (
				    gradz + asy * (gradzn - gradzp) * tz);
			}
		    }
/* *          endif          !!sp05 */
L5800:
		    ;
		}
/* L6000: */
	    }
/* print out the average no. of particles in regions where the local */
/* baryon density is higher than 5*rho0 */
/*       write(1072,'(e10.3,2x,e10.3)')nt*dt,float(nbaryn)/float(num) */
/* print out the average no. of particles in regions where the local */
/* energy density is higher than 2 GeV/fm^3. */
/*       write(1073,'(e10.3,2x,e10.3)')nt*dt,float(nquark)/float(num) */
/* print out the no. of particles that have run off the lattice */
/*          IF (NLOST .NE. 0 .AND. (NT/NFREQ)*NFREQ .EQ. NT) THEN */
/*            WRITE(12,'(5X,''***'',I7,'' TESTPARTICLES LOST AFTER '', */
/*     &                   ''TIME STEP NUMBER'',I4)') NLOST, NT */
/*         END IF */

/*       update phase space density */
/*        call platin(mode,mass,num,dx,dy,dz,dpx,dpy,dpz,fnorm) */

/*       CONTROL-PRINTOUT OF CONFIGURATION (IF REQUIRED) */

/*        if (inout(5) .eq. 2) CALL ENERGY(NT,IPOT,NUM,MASS,EMIN,EMAX) */


/* print out central baryon density as a function of time */
	    cden = dd_1.rho[41184] / .168f;
/* c        WRITE(1002,990)FLOAT(NT)*DT,CDEN */
/*        WRITE(1002,1990)FLOAT(NT)*DT,CDEN,denst/real(num) */
/* print out the central energy density as a function of time */
/* c        WRITE(1003,990)FLOAT(NT)*DT,PEL(0,0,0) */
/*        WRITE(1003,1990)FLOAT(NT)*DT,PEL(0,0,0),edenst/real(num) */
/* print out the no. of pion-like particles as a function of time */
/*        WRITE(1004,9999)FLOAT(NT)*DT,ALD,ALN,ALP,ALN5, */
/*     &               ALD+ALN+ALP+0.5*ALN5 */
/* print out the no. of eta-like particles as a function of time */
/*        WRITE(1005,991)FLOAT(NT)*DT,ALN5,ALE,ALE+0.5*ALN5 */
/* 990       FORMAT(E10.3,2X,E10.3) */
/* 1990       FORMAT(E10.3,2X,E10.3,2X,E10.3) */
/* 991       FORMAT(E10.3,2X,E10.3,2X,E10.3,2X,E10.3) */
/* 9999    FORMAT(e10.3,2X,e10.3,2X,E10.3,2X,E10.3,2X, */
/*     1  E10.3,2X,E10.3) */
/* THE FOLLOWING OUTPUTS CAN BE TURNED ON/OFF by setting icflow and icrho=0 */
/* print out the baryon and meson density matrix in the reaction plane */
	    if (nt / input2_1.nfreq * input2_1.nfreq == nt) {
		if (input2_1.icflow == 1) {
		    flow_(&nt);
		}
/* bz11/18/98 */
/*       if(icrho.ne.1)go to 10000 */
/*       if (icrho .eq. 1) then */
/* bz11/18/98end */
/*       do ix=-10,10 */
/*       do iz=-10,10 */
/*       write(1053,992)ix,iz,rho(ix,0,iz)/0.168 */
/*       write(1054,992)ix,iz,pirho(ix,0,iz)/0.168 */
/*       write(1055,992)ix,iz,pel(ix,0,iz) */
/*       end do */
/*       end do */
/* bz11/18/98 */
/*        end if */
/* bz11/18/98end */
/* 992       format(i3,i3,e11.4) */
	    }
/* print out the ENERGY density matrix in the reaction plane */
/* CHECK LOCAL MOMENTUM EQUILIBRIUM IN EACH CELL, */
/* AND PERFORM ON-LINE FLOW ANALYSIS AT A FREQUENCY OF NFREQ */
/*        IF ((NT/NFREQ)*NFREQ .EQ. NT ) THEN */
/*       call flow(nt) */
/*       call equ(ipot,mass,num,outpar) */
/*       do ix=-10,10 */
/*       do iz=-10,10 */
/*       write(1055,992)ix,iz,pel(ix,0,iz) */
/*       write(1056,992)ix,iz,rxy(ix,0,iz) */
/*       end do */
/*       end do */
/*       endif */
/* calculate the volume of high BARYON AND ENERGY density */
/* matter as a function of time */
/*       vbrho=0. */
/*       verho=0. */
/*       do ix=-20,20 */
/*       do iy=-20,20 */
/*       do iz=-20,20 */
/*       if(rho(ix,iy,iz)/0.168.gt.5.)vbrho=vbrho+1. */
/*       if(pel(ix,iy,iz).gt.2.)verho=verho+1. */
/*       end do */
/*       end do */
/*       end do */
/*       write(1081,993)dt*nt,vbrho */
/*       write(1082,993)dt*nt,verho */
/* 993       format(e11.4,2x,e11.4) */
/* ----------------------------------------------------------------------- */
/* bz11/16/98 */
/* .....for read-in initial conditions produce particles from read-in */
/* .....common block. */
/* .....note that this part is only for cascade with number of test particles */
/* .....NUM = 1. */
	    if (arprnt_1.iapar2[0] != 1) {
		ct = nt * input1_1.dt;
/* bz12/22/98 */
/*         NP = MASSR(1) */
/*         DO WHILE (FTAR(NPI) .GT. CT - DT .AND. FTAR(NPI) .LE. CT) */
/*            NP = NP + 1 */
/*            R(1, NP) = GXAR(NPI) + PXAR(NPI) / PEAR(NPI) * (CT - FTAR(NPI)) */
/*            R(2, NP) = GYAR(NPI) + PYAR(NPI) / PEAR(NPI) * (CT - FTAR(NPI)) */
/*            R(3, NP) = GZAR(NPI) + PZAR(NPI) / PEAR(NPI) * (CT - FTAR(NPI)) */
/*            P(1, NP) = PXAR(NPI) */
/*            P(2, NP) = PYAR(NPI) */
/*            P(3, NP) = PZAR(NPI) */
/*            E(NP) = XMAR(NPI) */
/*            LB(NP) = IARFLV(ITYPAR(NPI)) */
/*            NPI = NPI + 1 */
/*         END DO */
/*         MASSR(1) = NP */
		ia = 0;
		i__4 = run_1.num;
		for (irun = 1; irun <= i__4; ++irun) {
		    i__3 = rr_1.massr[irun];
		    for (ic = 1; ic <= i__3; ++ic) {
			ie = ia + ic;
			rt[(ic + irun * 150001) * 3 - 450006] = aa_1.r__[ie * 
				3 - 3];
			rt[(ic + irun * 150001) * 3 - 450005] = aa_1.r__[ie * 
				3 - 2];
			rt[(ic + irun * 150001) * 3 - 450004] = aa_1.r__[ie * 
				3 - 1];
			pt[(ic + irun * 150001) * 3 - 450006] = bb_1.p[ie * 3 
				- 3];
			pt[(ic + irun * 150001) * 3 - 450005] = bb_1.p[ie * 3 
				- 2];
			pt[(ic + irun * 150001) * 3 - 450004] = bb_1.p[ie * 3 
				- 1];
			et[ic + irun * 150001 - 150002] = cc_1.e[ie - 1];
			lt[ic + irun * 150001 - 150002] = ee_1.lb[ie - 1];
/*         !! sp 12/19/00 */
			prot[ic + irun * 150001 - 150002] = hh_1.proper[ie - 
				1];
/* lin-5/2008: */
			dpert_1.dpertt[ic + irun * 150001 - 150002] = 
				dpert_1.dpertp[ie - 1];
/* L1027: */
		    }
		    np = rr_1.massr[irun];
		    np1 = npi[irun - 1];
/* bz10/05/99 */
/*            DO WHILE (FT1(NP1, IRUN) .GT. CT - DT .AND. */
/*     &           FT1(NP1, IRUN) .LE. CT) */
/* bz10/06/99 */
/*            DO WHILE (NPI(IRUN).LE.MULTI1(IRUN).AND. */
/* bz10/06/99 end */
/* lin-11/13/00 finally read in all unformed particles and do the decays in ART: */
/*           DO WHILE (NP1.LE.MULTI1(IRUN).AND. */
/*    &           FT1(NP1, IRUN) .GT. CT - DT .AND. */
/*    &           FT1(NP1, IRUN) .LE. CT) */

		    ctlong = ct;
		    if (nt == input2_1.ntmax - 1) {
			ctlong = 1e30f;
		    } else if (nt == input2_1.ntmax) {
			goto L1111;
		    }

		    while(np1 <= arerc1_1.multi1[irun - 1] && arprc1_1.ft1[
			    np1 + irun * 150001 - 150002] > (nt - 1) * 
			    input1_1.dt && arprc1_1.ft1[np1 + irun * 150001 - 
			    150002] <= ctlong) {
/* lin-ma-5/2016 changed the following to 2nd line above to avoid bug */
/*     that leads to loss of hadrons inside ART due to finite accuracy */
/*     [which results in (ct-dt) + dt != ct exactly]: */
/*     &           FT1(NP1, IRUN) .GT. (CT - DT) .AND. */
			++np;
			udt = (ct - arprc1_1.ft1[np1 + irun * 150001 - 150002]
				) / arprc1_1.ee1[np1 + irun * 150001 - 150002]
				;
/* lin-10/28/03 since all unformed hadrons at time ct are read in at nt=ntmax-1, */
/*     their positions should not be propagated to time ct: */
			if (nt == input2_1.ntmax - 1) {
			    ftmax_1.ftsvt[np + irun * 150001 - 150002] = 
				    arprc1_1.ft1[np1 + irun * 150001 - 150002]
				    ;
			    if (arprc1_1.ft1[np1 + irun * 150001 - 150002] > 
				    ct) {
				udt = 0.f;
			    }
			}
			rt[(np + irun * 150001) * 3 - 450006] = arprc1_1.gx1[
				np1 + irun * 150001 - 150002] + arprc1_1.px1[
				np1 + irun * 150001 - 150002] * udt;
			rt[(np + irun * 150001) * 3 - 450005] = arprc1_1.gy1[
				np1 + irun * 150001 - 150002] + arprc1_1.py1[
				np1 + irun * 150001 - 150002] * udt;
			rt[(np + irun * 150001) * 3 - 450004] = arprc1_1.gz1[
				np1 + irun * 150001 - 150002] + arprc1_1.pz1[
				np1 + irun * 150001 - 150002] * udt;
			pt[(np + irun * 150001) * 3 - 450006] = arprc1_1.px1[
				np1 + irun * 150001 - 150002];
			pt[(np + irun * 150001) * 3 - 450005] = arprc1_1.py1[
				np1 + irun * 150001 - 150002];
			pt[(np + irun * 150001) * 3 - 450004] = arprc1_1.pz1[
				np1 + irun * 150001 - 150002];
			et[np + irun * 150001 - 150002] = arprc1_1.xm1[np1 + 
				irun * 150001 - 150002];
			lt[np + irun * 150001 - 150002] = iarflv_(&
				arprc1_1.ityp1[np1 + irun * 150001 - 150002]);
/* lin-5/2008: */
			dpert_1.dpertt[np + irun * 150001 - 150002] = 
				dpert_1.dpp1[np1 + irun * 150001 - 150002];
/* lin-4/30/03 ctest off */
/*     record initial phi,K*,Lambda(1520) resonances formed during the timestep: */
/*               if(LT(NP, IRUN).eq.29.or.iabs(LT(NP, IRUN)).eq.30) */
/*     1              write(17,112) 'formed',LT(NP, IRUN),PX1(NP1, IRUN), */
/*     2 PY1(NP1, IRUN),PZ1(NP1, IRUN),XM1(NP1, IRUN),nt */
/* 112           format(a10,1x,I4,4(1x,f9.3),1x,I4) */

			++np1;
/*     !! sp 12/19/00 */
			prot[np + irun * 150001 - 150002] = 1.f;
		    }

L1111:
		    npi[irun - 1] = np1;
		    ia += rr_1.massr[irun];
		    rr_1.massr[irun] = np;
/* L1028: */
		}
		ia = 0;
		i__4 = run_1.num;
		for (irun = 1; irun <= i__4; ++irun) {
		    ia += rr_1.massr[irun - 1];
		    i__3 = rr_1.massr[irun];
		    for (ic = 1; ic <= i__3; ++ic) {
			ie = ia + ic;
			aa_1.r__[ie * 3 - 3] = rt[(ic + irun * 150001) * 3 - 
				450006];
			aa_1.r__[ie * 3 - 2] = rt[(ic + irun * 150001) * 3 - 
				450005];
			aa_1.r__[ie * 3 - 1] = rt[(ic + irun * 150001) * 3 - 
				450004];
			bb_1.p[ie * 3 - 3] = pt[(ic + irun * 150001) * 3 - 
				450006];
			bb_1.p[ie * 3 - 2] = pt[(ic + irun * 150001) * 3 - 
				450005];
			bb_1.p[ie * 3 - 1] = pt[(ic + irun * 150001) * 3 - 
				450004];
			cc_1.e[ie - 1] = et[ic + irun * 150001 - 150002];
			ee_1.lb[ie - 1] = lt[ic + irun * 150001 - 150002];
/*     !! sp 12/19/00 */
			hh_1.proper[ie - 1] = prot[ic + irun * 150001 - 
				150002];
			if (nt == input2_1.ntmax - 1) {
			    ftmax_1.ftsv[ie - 1] = ftmax_1.ftsvt[ic + irun * 
				    150001 - 150002];
			}
/* lin-5/2008: */
			dpert_1.dpertp[ie - 1] = dpert_1.dpertt[ic + irun * 
				150001 - 150002];
/* L1029: */
		    }
/* lin-3/2009 Moved here to better take care of freezeout spacetime: */
		    hbtout_(&rr_1.massr[irun], &nt, &input2_1.ntmax);
/* L1030: */
		}
/* bz12/22/98end */
	    }
/* bz11/16/98end */
/* lin-5/2009 ctest off: */
/*      call flowh(ct) */
/* L10000: */
	}
/*                                                                      * */
/*       ==============  END OF TIME STEP LOOP   ================       * */
/* *********************************** */
/*     WRITE OUT particle's MOMENTA ,and/OR COORDINATES , */
/*     label and/or their local baryon density in the final state */
	iss = 0;
	i__2 = run_1.num;
	for (lrun = 1; lrun <= i__2; ++lrun) {
	    iss += rr_1.massr[lrun - 1];
	    i__4 = rr_1.massr[lrun];
	    for (l0 = 1; l0 <= i__4; ++l0) {
		ipart = iss + l0;
/* L1031: */
	    }
/* L1032: */
	}
/* bz11/16/98 */
	if (arprnt_1.iapar2[0] != 1) {
/* bz12/22/98 */
/*        NSH = MASSR(1) - NPI + 1 */
/*        IAINT2(1) = IAINT2(1) + NSH */
/* .....to shift the unformed particles to the end of the common block */
/*        IF (NSH .GT. 0) THEN */
/*           IB = IAINT2(1) */
/*           IE = MASSR(1) + 1 */
/*           II = -1 */
/*        ELSE IF (NSH .LT. 0) THEN */
/*           IB = MASSR(1) + 1 */
/*           IE = IAINT2(1) */
/*           II = 1 */
/*        END IF */
/*        IF (NSH .NE. 0) THEN */
/*           DO I = IB, IE, II */
/*              J = I - NSH */
/*              ITYPAR(I) = ITYPAR(J) */
/*              GXAR(I) = GXAR(J) */
/*              GYAR(I) = GYAR(J) */
/*              GZAR(I) = GZAR(J) */
/*              FTAR(I) = FTAR(J) */
/*              PXAR(I) = PXAR(J) */
/*              PYAR(I) = PYAR(J) */
/*              PZAR(I) = PZAR(J) */
/*              PEAR(I) = PEAR(J) */
/*              XMAR(I) = XMAR(J) */
/*           END DO */
/*        END IF */
/* .....to copy ART particle info to COMMON /ARPRC/ */
/*        DO I = 1, MASSR(1) */
/*           ITYPAR(I) = INVFLV(LB(I)) */
/*           GXAR(I) = R(1, I) */
/*           GYAR(I) = R(2, I) */
/*           GZAR(I) = R(3, I) */
/*           FTAR(I) = CT */
/*           PXAR(I) = P(1, I) */
/*           PYAR(I) = P(2, I) */
/*           PZAR(I) = P(3, I) */
/*           XMAR(I) = E(I) */
/*           PEAR(I) = SQRT(PXAR(I) ** 2 + PYAR(I) ** 2 + PZAR(I) ** 2 */
/*     &        + XMAR(I) ** 2) */
/*        END DO */
	    ia = 0;
	    i__2 = run_1.num;
	    for (irun = 1; irun <= i__2; ++irun) {
		ia += rr_1.massr[irun - 1];
		np1 = npi[irun - 1];
		nsh = rr_1.massr[irun] - np1 + 1;
		arerc1_1.multi1[irun - 1] += nsh;
/* .....to shift the unformed particles to the end of the common block */
		if (nsh > 0) {
		    ib = arerc1_1.multi1[irun - 1];
		    ie = rr_1.massr[irun] + 1;
		    ii = -1;
		} else if (nsh < 0) {
		    ib = rr_1.massr[irun] + 1;
		    ie = arerc1_1.multi1[irun - 1];
		    ii = 1;
		}
		if (nsh != 0) {
		    i__4 = ie;
		    i__3 = ii;
		    for (i__ = ib; i__3 < 0 ? i__ >= i__4 : i__ <= i__4; i__ 
			    += i__3) {
			j = i__ - nsh;
			arprc1_1.ityp1[i__ + irun * 150001 - 150002] = 
				arprc1_1.ityp1[j + irun * 150001 - 150002];
			arprc1_1.gx1[i__ + irun * 150001 - 150002] = 
				arprc1_1.gx1[j + irun * 150001 - 150002];
			arprc1_1.gy1[i__ + irun * 150001 - 150002] = 
				arprc1_1.gy1[j + irun * 150001 - 150002];
			arprc1_1.gz1[i__ + irun * 150001 - 150002] = 
				arprc1_1.gz1[j + irun * 150001 - 150002];
			arprc1_1.ft1[i__ + irun * 150001 - 150002] = 
				arprc1_1.ft1[j + irun * 150001 - 150002];
			arprc1_1.px1[i__ + irun * 150001 - 150002] = 
				arprc1_1.px1[j + irun * 150001 - 150002];
			arprc1_1.py1[i__ + irun * 150001 - 150002] = 
				arprc1_1.py1[j + irun * 150001 - 150002];
			arprc1_1.pz1[i__ + irun * 150001 - 150002] = 
				arprc1_1.pz1[j + irun * 150001 - 150002];
			arprc1_1.ee1[i__ + irun * 150001 - 150002] = 
				arprc1_1.ee1[j + irun * 150001 - 150002];
			arprc1_1.xm1[i__ + irun * 150001 - 150002] = 
				arprc1_1.xm1[j + irun * 150001 - 150002];
/*     !! sp 12/19/00 */
			arercp_1.pro1[i__ + irun * 150001 - 150002] = 
				arercp_1.pro1[j + irun * 150001 - 150002];
/* lin-5/2008: */
			dpert_1.dpp1[i__ + irun * 150001 - 150002] = 
				dpert_1.dpp1[j + irun * 150001 - 150002];
/* L1033: */
		    }
		}
/* .....to copy ART particle info to COMMON /ARPRC1/ */
		i__3 = rr_1.massr[irun];
		for (i__ = 1; i__ <= i__3; ++i__) {
		    ib = ia + i__;
		    arprc1_1.ityp1[i__ + irun * 150001 - 150002] = invflv_(&
			    ee_1.lb[ib - 1]);
		    arprc1_1.gx1[i__ + irun * 150001 - 150002] = aa_1.r__[ib *
			     3 - 3];
		    arprc1_1.gy1[i__ + irun * 150001 - 150002] = aa_1.r__[ib *
			     3 - 2];
		    arprc1_1.gz1[i__ + irun * 150001 - 150002] = aa_1.r__[ib *
			     3 - 1];
/* lin-10/28/03: */
/* since all unformed hadrons at time ct are read in at nt=ntmax-1, */
/* their formation time ft1 should be kept to determine their freezeout(x,t): */
/*              FT1(I, IRUN) = CT */
		    if (arprc1_1.ft1[i__ + irun * 150001 - 150002] < ct) {
			arprc1_1.ft1[i__ + irun * 150001 - 150002] = ct;
		    }
		    arprc1_1.px1[i__ + irun * 150001 - 150002] = bb_1.p[ib * 
			    3 - 3];
		    arprc1_1.py1[i__ + irun * 150001 - 150002] = bb_1.p[ib * 
			    3 - 2];
		    arprc1_1.pz1[i__ + irun * 150001 - 150002] = bb_1.p[ib * 
			    3 - 1];
		    arprc1_1.xm1[i__ + irun * 150001 - 150002] = cc_1.e[ib - 
			    1];
/* Computing 2nd power */
		    r__1 = arprc1_1.px1[i__ + irun * 150001 - 150002];
/* Computing 2nd power */
		    r__2 = arprc1_1.py1[i__ + irun * 150001 - 150002];
/* Computing 2nd power */
		    r__3 = arprc1_1.pz1[i__ + irun * 150001 - 150002];
/* Computing 2nd power */
		    r__4 = arprc1_1.xm1[i__ + irun * 150001 - 150002];
		    arprc1_1.ee1[i__ + irun * 150001 - 150002] = sqrt(r__1 * 
			    r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
/*     !! sp 12/19/00 */
		    arercp_1.pro1[i__ + irun * 150001 - 150002] = hh_1.proper[
			    ib - 1];
/* L1034: */
		}
/* L1035: */
	    }
/* bz12/22/98end */
	}
/* bz11/16/98end */

/* ********************************* */
/*                                                                      * */
/*       ======= END OF MANY LOOPS OVER IMPACT PARAMETERS ==========    * */
/*                                                               * */
/* ********************************* */
/* L50000: */
    }

/* ----------------------------------------------------------------------- */
/*                       ==== ART COMPLETED ==== */
/* ----------------------------------------------------------------------- */
/* bz11/16/98 */
/*      STOP */
    return 0;
/* bz11/16/98end */
} /* artmn_ */

/* ********************************* */
/* Subroutine */ int coulin_(integer *masspr, integer *massta, integer *num)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, mass, irun;

/*                                                                      * */
/*     purpose:   initialization of array zet() and lb() for all runs  * */
/*                lb(i) = 1   =>  proton                               * */
/*                lb(i) = 2   =>  neutron                              * */
/* ********************************* */
/* c      SAVE /EE/ */
/* c      SAVE /zz/ */
    mass = *massta + *masspr;
    i__1 = *num;
    for (irun = 1; irun <= i__1; ++irun) {
	i__2 = zz_1.zta + (irun - 1) * mass;
	for (i__ = (irun - 1) * mass + 1; i__ <= i__2; ++i__) {
	    ee_1.lb[i__ - 1] = 1;
/* L100: */
	}
	i__2 = *massta + (irun - 1) * mass;
	for (i__ = zz_1.zta + 1 + (irun - 1) * mass; i__ <= i__2; ++i__) {
	    ee_1.lb[i__ - 1] = 2;
/* L200: */
	}
	i__2 = *massta + zz_1.zpr + (irun - 1) * mass;
	for (i__ = *massta + 1 + (irun - 1) * mass; i__ <= i__2; ++i__) {
	    ee_1.lb[i__ - 1] = 1;
/* L300: */
	}
	i__2 = *massta + *masspr + (irun - 1) * mass;
	for (i__ = *massta + zz_1.zpr + 1 + (irun - 1) * mass; i__ <= i__2; 
		++i__) {
	    ee_1.lb[i__ - 1] = 2;
/* L400: */
	}
/* L500: */
    }
    return 0;
} /* coulin_ */

/* ********************************* */
/*                                                                      * */
/* Subroutine */ int relcol_(integer *lcoll, integer *lbloc, integer *lcnne, 
	integer *ldd, integer *lpp, integer *lppk, integer *lpn, integer *lpd,
	 integer *lrho, integer *lomega, integer *lkn, integer *lnnk, integer 
	*lddk, integer *lndk, integer *lcnnd, integer *lcndn, integer *ldirt, 
	integer *ldecay, integer *lres, integer *ldou, integer *lddrho, 
	integer *lnnrho, integer *lnnom, integer *nt, integer *ntmax, real *
	sp, real *akaon, real *sk)
{
    /* Initialized data */

    static real zet[91] = { 1.f,0.f,0.f,0.f,0.f,1.f,0.f,0.f,0.f,0.f,0.f,0.f,
	    0.f,0.f,0.f,-1.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,
	    -1.f,0.f,1.f,0.f,-1.f,0.f,-1.f,0.f,-2.f,-1.f,0.f,1.f,0.f,0.f,0.f,
	    0.f,-1.f,0.f,1.f,0.f,-1.f,0.f,1.f,-1.f,0.f,1.f,2.f,0.f,1.f,0.f,
	    1.f,0.f,-1.f,0.f,1.f,0.f,0.f,0.f,-1.f,0.f,1.f,0.f,-1.f,0.f,1.f,
	    0.f,0.f,1.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,-1.f,0.f,0.f,0.f,
	    0.f,-1.f };

    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    real r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);
    integer i_nint(real *), s_wsle(cilist *), do_lio(integer *, integer *, 
	    char *, ftnlen), e_wsle(void);
    double log(doublereal);

    /* Local variables */
    static integer i__, j, n;
    static real e2;
    static integer i0, i1, j1, j2, i2, n0;
    static real t0, x1, y1, z1, x2, y2, z2, ec;
    static integer ia, j10, ic, ib, ie, ig, il, im;
    static real ds, et[150001];
    static integer kp, lt[150001], is;
    static real pt[450003]	/* was [3][150001] */, rt[450003]	/* 
	    was [3][150001] */;
    static integer in;
    static real xx, yy, ec0;
    extern /* Subroutine */ int sdbelastic_(real *, real *);
    static integer id1;
    static real am1;
    static integer id2;
    static real am2, er1, er2, pk0;
    extern doublereal pp1_(real *);
    static integer ix1, iy1, iz1, ix2, iy2, iz2;
    extern doublereal pp2_(real *);
    static real px2, py2, pz2, xx0, ece, sdb;
    static integer lbm;
    static real dse;
    extern doublereal w1440_(real *);
    extern /* Subroutine */ int cms_(integer *, integer *, real *, real *, 
	    real *, real *);
    static real wid;
    extern doublereal w1535_(real *);
    static real sig, pcx, pcy, pcz, srt;
    static integer iss;
    extern /* Subroutine */ int xnd_(real *, real *, real *, real *, integer *
	    , integer *, real *, real *, real *, real *, real *, real *, real 
	    *);
    static real xmm;
    extern doublereal ppt_(real *), xpp_(real *), xnp_(real *);
    static real dsr, sdm;
    static integer ipp;
    static real ert;
    static integer ikk;
    static real e1cm;
    static integer ilb1, ilb2, lb1i, lb2i;
    static real em1i, em2i, e2cm;
    static integer lbp1;
    static real emm1, emm2;
    static integer lbp2, ipx1, ipy1, ipz1;
    static real px1i, py1i, pz1i, px2i, py2i, pz2i;
    static integer ipx2, ipy2, ipz2;
    static real xsk1, xsk2, xsk3, xsk4, xsk5, xsk6, xsk7, xky1, xky2;
    extern doublereal reab_(integer *, integer *, real *, integer *);
    static real xky3, xky4;
    extern /* Subroutine */ int crdd_(integer *, real *, real *, real *, real 
	    *, integer *, integer *, integer *, integer *, real *, real *, 
	    integer *, integer *);
    static real xky5;
    extern /* Subroutine */ int crhb_(real *, real *, real *, real *, integer 
	    *, integer *, integer *);
    static real xky6, xky7, eini, brel;
    extern /* Subroutine */ int crnd_(integer *, real *, real *, real *, real 
	    *, integer *, integer *, integer *, real *, real *, real *, real *
	    , real *, real *, real *, real *, integer *, integer *), cren_(
	    real *, real *, real *, real *, integer *, integer *, integer *);
    static integer ntag;
    extern /* Subroutine */ int crrd_(real *, real *, real *, real *, integer 
	    *, integer *, integer *, real *, real *, real *, real *), crpd_(
	    real *, real *, real *, real *, integer *, integer *, integer *, 
	    real *, real *, real *, real *);
    static real sigk;
    static integer lpdr;
    static real prot[150001];
    static integer msum, irun, mass;
    static real xeta;
    extern doublereal xn1535_(integer *, integer *, integer *);
    static real xphi, xmax;
    extern /* Subroutine */ int wida1_(real *, real *, real *, integer *);
    extern doublereal xnpi_(integer *, integer *, integer *, real *);
    static real xres;
    extern /* Subroutine */ int crpn_(real *, real *, real *, real *, integer 
	    *, integer *, integer *, real *, real *, real *, real *), crnn_(
	    integer *, real *, real *, real *, real *, integer *, integer *, 
	    integer *, integer *, real *, real *, integer *, integer *);
    static real ppel;
    extern /* Subroutine */ int ppxs_(integer *, integer *, real *, real *, 
	    real *, integer *);
    static real ppin, dspp;
    extern /* Subroutine */ int crpp_(real *, real *, real *, real *, integer 
	    *, integer *, integer *, real *, real *, real *, integer *);
    static real dskk, dskk0, sigp, pzrt;
    static integer ikkg;
    static real sigr0;
    static integer ikkl;
    static real dshn, xky8, px1cm, py1cm, pz1cm;
    extern doublereal pipp1_(real *);
    static real xky9, xky10, xky11, xky12, xky13, xky14, xky15, xky16, xky17;
    static integer ikmp;
    static real dskn;
    extern /* Subroutine */ int crkn_(real *, real *, real *, real *, integer 
	    *, integer *, integer *);
    static real pt1i1, pt2i1, pt3i1, pt1i2, pt2i2, pt3i2;
    static integer icase;
    extern /* Subroutine */ int decay_(integer *, integer *, integer *, 
	    integer *, real *, integer *), sbbdm_(real *, real *, integer *, 
	    integer *, real *, real *), sdmbb_(real *, real *, integer *);
    extern doublereal akpel_(real *), width_(real *);
    static real drmax;
    static integer lrhor;
    static real rhomp, pxini, pyini, pzini, spipi, bmass;
    extern /* Subroutine */ int decay2_(integer *, integer *, integer *, 
	    integer *, real *, integer *);
    static real pkaon;
    extern doublereal aknel_(real *);
    static real brsgm, brsig;
    static integer nchrg;
    extern /* Subroutine */ int newka_(integer *, integer *, integer *, real *
	    , integer *, integer *, integer *, integer *, real *, real *, 
	    real *, real *, integer *);
    static integer ictrl;
    extern doublereal pnlka_(real *);
    static real sigma0;
    extern doublereal pnska_(real *);
    static real xkaon, xphin, xmaxn;
    extern /* Subroutine */ int distc0_(real *, real *, real *, integer *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *);
    extern doublereal dirct1_(real *), dirct2_(real *), dirct3_(real *);
    static real deltr0, dr0max, xnpin, xnpid, xreab;
    extern doublereal dpion_(real *, real *, integer *, integer *, real *);
    static real xkaon0;
    extern /* Subroutine */ int crdir_(real *, real *, real *, real *, 
	    integer *, integer *, integer *);
    static real spika;
    extern doublereal erhon_(real *, real *, integer *, integer *, real *);
    static real signn, xinel;
    static integer ianti, ipert1;
    static real signn0;
    extern /* Subroutine */ int xddin_(real *, real *, real *, real *, 
	    integer *, integer *, real *, real *, real *, real *, real *, 
	    real *, real *);
    static real ppsig, ppink, xmaxn1, xnpin1;
    extern doublereal pipik_(real *);
    extern /* Subroutine */ int spprr_(integer *, integer *, real *), sppee_(
	    integer *, integer *, real *), spppe_(integer *, integer *, real *
	    ), srpre_(integer *, integer *, real *), sopoe_(integer *, 
	    integer *, real *), srree_(integer *, integer *, real *);
    static real dsppr;
    static integer icheck;
    static real dsppb;
    extern /* Subroutine */ int crlaba_(real *, real *, real *, real *, real *
	    , real *, integer *, integer *, integer *, integer *, integer *, 
	    integer *), crdmbb_(real *, real *, real *, real *, integer *, 
	    integer *, integer *, integer *, real *, integer *, integer *), 
	    xphib_(integer *, integer *, real *, real *, real *, real *, real 
	    *, real *, real *, real *, real *), crdbel_(real *, real *, real *
	    , real *, integer *, integer *, integer *, integer *, real *, 
	    integer *, integer *);
    static real dshnr, sigkp;
    static integer idecay;
    extern /* Subroutine */ int lambar_(integer *, integer *, real *, real *);
    static integer ichann;
    static real siglab;
    static integer iblock;
    static real pdecay, gfactr;
    extern /* Subroutine */ int resdec_(integer *, integer *, integer *, real 
	    *, integer *, integer *);
    static real sigela;
    extern doublereal akplam_(real *), aknlam_(real *);
    static real xdecay;
    static integer inewka;
    extern /* Subroutine */ int inidcy_(void);
    static real dptemp[150001], fttemp[150001];
    static integer massrn[2];
    static real ftpisv[150001]	/* was [150001][1] */, resona;
    static integer nodelt, lomgar;
    extern doublereal ranart_(integer *);
    static real pnstar;
    static integer nnnini;
    static real rppmax, rsqare;
    static integer ifirst;
    static real sigsgm;
    extern doublereal akpsgm_(real *), aknsgm_(real *);
    static real deltre;
    extern /* Subroutine */ int distce_(integer *, integer *, real *, real *, 
	    real *, real *, real *, integer *, real *, real *, real *), 
	    pibphi_(real *, integer *, integer *, real *, real *, real *, 
	    real *);
    extern doublereal pionpp_(real *);
    static real sumsrt, xdirct, xnelas, deltar;
    extern /* Subroutine */ int dreson_(integer *, integer *), crkpla_(real *,
	     real *, real *, real *, real *, real *, real *, real *, integer *
	    , integer *, integer *, integer *, integer *, real *), ksreso_(
	    integer *, integer *);
    static real xelstc, cutoff, sdprod, pfinal, dspert, spprho;
    extern /* Subroutine */ int getnst_(real *);
    extern doublereal ppbbar_(real *), prbbar_(real *), rrbbar_(real *), 
	    pobbar_(real *), robbar_(real *), oobbar_(real *), xppbar_(real *)
	    ;
    static real dsppbr;
    extern /* Subroutine */ int crppba_(real *, real *, real *, real *, 
	    integer *, integer *, integer *), pertur_(real *, real *, real *, 
	    real *, integer *, integer *, integer *, integer *, integer *, 
	    integer *);
    static integer icontp;
    extern /* Subroutine */ int crphib_(real *, real *, real *, real *, 
	    integer *, integer *, real *, real *, real *, real *, real *, 
	    real *, integer *), phimes_(integer *, integer *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *);
    static real sigphi;
    extern /* Subroutine */ int crphim_(real *, real *, real *, real *, 
	    integer *, integer *, real *, real *, real *, real *, real *, 
	    real *, real *, integer *, integer *, integer *), xkhype_(integer 
	    *, integer *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *), crkhyp_(real *, real *, 
	    real *, real *, integer *, integer *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, integer *,
	     integer *), crlan_(real *, real *, real *, real *, integer *, 
	    integer *, integer *), crkphi_(real *, real *, real *, real *, 
	    real *, integer *, real *, real *, integer *, integer *, integer *
	    , integer *, integer *, integer *, real *, real *), crksph_(real *
	    , real *, real *, real *, real *, real *, real *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, real *);
    static real dsknr, p1beta, transf;
    static integer ipion, ipdflag;
    static real dsrpert;

    /* Fortran I/O blocks */
    static cilist io___379 = { 0, 6, 0, 0, 0 };


/*                                                                      * */
/*       PURPOSE:    CHECK CONDITIONS AND CALCULATE THE KINEMATICS      * */
/*                   FOR BINARY COLLISIONS AMONG PARTICLES              * */
/*                                 - RELATIVISTIC FORMULA USED          * */
/*                                                                      * */
/*       REFERENCES: HAGEDORN, RELATIVISTIC KINEMATICS (1963)           * */
/*                                                                      * */
/*       VARIABLES:                                                     * */
/*         MASSPR  - NUMBER OF NUCLEONS IN PROJECTILE   (INTEGER,INPUT) * */
/*         MASSTA  - NUMBER OF NUCLEONS IN TARGET       (INTEGER,INPUT) * */
/*         NUM     - NUMBER OF TESTPARTICLES PER NUCLEON(INTEGER,INPUT) * */
/*         ISEED   - SEED FOR RANDOM NUMBER GENERATOR   (INTEGER,INPUT) * */
/*         IAVOID  - (= 1 => AVOID FIRST CLLISIONS WITHIN THE SAME      * */
/*                   NUCLEUS, ELSE ALL COLLISIONS)      (INTEGER,INPUT) * */
/*         DELTAR  - MAXIMUM SPATIAL DISTANCE FOR WHICH A COLLISION     * */
/*                   STILL CAN OCCUR                       (REAL,INPUT) * */
/*         DT      - TIME STEP SIZE                        (REAL,INPUT) * */
/*         LCOLL   - NUMBER OF COLLISIONS              (INTEGER,OUTPUT) * */
/*         LBLOC   - NUMBER OF PULI-BLOCKED COLLISIONS (INTEGER,OUTPUT) * */
/*         LCNNE   - NUMBER OF ELASTIC COLLISION       (INTEGER,OUTPUT) * */
/*         LCNND   - NUMBER OF N+N->N+DELTA REACTION   (INTEGER,OUTPUT) * */
/*         LCNDN   - NUMBER OF N+DELTA->N+N REACTION   (INTEGER,OUTPUT) * */
/*         LDD     - NUMBER OF RESONANCE+RESONANCE COLLISIONS */
/*         LPP     - NUMBER OF PION+PION elastic COLIISIONS */
/*         lppk    - number of pion(RHO,OMEGA)+pion(RHO,OMEGA) */
/*                   -->K+K- collisions */
/*         LPN     - NUMBER OF PION+N-->KAON+X */
/*         lpd     - number of pion+n-->delta+pion */
/*         lrho    - number of pion+n-->Delta+rho */
/*         lomega  - number of pion+n-->Delta+omega */
/*         LKN     - NUMBER OF KAON RESCATTERINGS */
/*         LNNK    - NUMBER OF bb-->kAON PROCESS */
/*         LDDK    - NUMBER OF DD-->KAON PROCESS */
/*         LNDK    - NUMBER OF ND-->KAON PROCESS */
/*         LB(I) IS USED TO LABEL PARTICLE'S CHARGE STATE */
/*         LB(I)   = */
/* bali2/7/99 */
/*                 -45 Omega baryon(bar) */
/*                 -41 cascade0(bar) */
/*                 -40 cascade-(bar) */
/* lin-11/07/00: */
/*                 -30 K*- */
/*                 -17 sigma+(bar) */
/*                 -16 sigma0(bar) */
/*                 -15 sigma-(bar) */
/*                 -14 LAMBDA(bar) */
/* lin-8/29/00 */
/*                 -13 anti-N*(+1)(1535),s_11 */
/*                 -12 anti-N*0(1535),s_11 */
/*                 -11 anti-N*(+1)(1440),p_11 */
/*                 -10 anti-N*0(1440), p_11 */
/*                  -9 anti-DELTA+2 */
/*                  -8 anti-DELTA+1 */
/*                  -7 anti-DELTA0 */
/*                  -6 anti-DELTA-1 */

/*                  -2 antineutron */
/*                  -1 antiproton */
/* bali2/7/99end */
/*                   0 eta */
/*                   1 PROTON */
/*                   2 NUETRON */
/*                   3 PION- */
/*                   4 PION0 */
/*                   5 PION+ */
/*                   6 DELTA-1 */
/*                   7 DELTA0 */
/*                   8 DELTA+1 */
/*                   9 DELTA+2 */
/*                   10 N*0(1440), p_11 */
/*                   11 N*(+1)(1440),p_11 */
/*                  12 N*0(1535),s_11 */
/*                  13 N*(+1)(1535),s_11 */
/*                  14 LAMBDA */
/*                   15 sigma- */
/*                   16 sigma0 */
/*                   17 sigma+ */
/*                   21 kaon- */
/* lin-2/23/03        22 Kaon0Long (converted at the last timestep) */
/*                   23 KAON+ */
/*                   24 Kaon0short (converted at the last timestep then decay) */
/*                   25 rho- */
/*                   26 rho0 */
/*                   27 rho+ */
/*                   28 omega meson */
/*                   29 phi */
/*                   30 K*+ */
/* sp01/03/01 */
/*                   31 eta-prime */
/*                   40 cascade- */
/*                   41 cascade0 */
/*                   45 Omega baryon */
/* sp01/03/01 end */

/*                   ++  ------- SEE NOTE BOOK */
/*         NSTAR=1 INCLUDING N* RESORANCE */
/*         ELSE DELTA RESORANCE ONLY */
/*         NDIRCT=1 INCLUDING DIRECT PROCESS,ELSE NOT */
/*         DIR - PERCENTAGE OF DIRECT PION PRODUCTION PROCESS */
/* ********************************* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /DD/ */
/* c      SAVE /EE/ */
/* c      SAVE /HH/ */
/* c      SAVE /ff/ */
/* c      SAVE /gg/ */
/* c      SAVE /INPUT/ */
/* c      SAVE /NN/ */
/* c      SAVE /RR/ */
/* c      SAVE /ss/ */
/* c      SAVE /BG/ */
/* c      SAVE /RUN/ */
/* c      SAVE /PA/ */
/* c      SAVE /PB/ */
/* c      SAVE /PC/ */
/* c      SAVE /PD/ */
/* c      SAVE /PE/ */
/* c      SAVE /KKK/ */
/* c      SAVE /KAON/ */
/* c      SAVE /TABLE/ */
/* c      SAVE /input1/ */
/* c      SAVE /leadng/ */
/* c      SAVE /tdecay/ */
/* c      SAVE /lastt/ */

/* c      SAVE /ppbmas/ */
/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
/* c      SAVE /hbt/ */
/* c      SAVE /resdcy/ */
/* c      SAVE /RNDF77/ */
/* lin-5/2008: */

/* lin-2/19/03 initialize n and nsav for resonance decay at each timestep */
/*     in order to prevent integer overflow: */
    inidcy_();
/* OFF skip ART collisions to reproduce HJ: */
/* c       if(nt.ne.ntmax) return */
/* lin-11/07/00 rrkk is assumed to be 0.6mb(default) for mm->KKbar */
/*     with m=rho or omega, estimated from Ko's paper: */
/*      rrkk=0.6 */
/* prkk: cross section of pi (rho or omega) -> K* Kbar (AND) K*bar K: */
/*      prkk=0.3 */
/*     cross section in mb for (rho or omega) K* -> pi K: */
/*      srhoks=5. */
/* lin-11/07/00-end */
/*      ESBIN=0.04 */
    resona = 5.f;
/* ----------------------------------------------------------------------- */
/*     INITIALIZATION OF COUNTING VARIABLES */
    nodelt = 0;
    sumsrt = 0.f;
    *lcoll = 0;
    *lbloc = 0;
    *lcnne = 0;
    *ldd = 0;
    *lpp = 0;
    *lpd = 0;
    lpdr = 0;
    *lrho = 0;
    lrhor = 0;
    *lomega = 0;
    lomgar = 0;
    *lpn = 0;
    *lkn = 0;
    *lnnk = 0;
    *lddk = 0;
    *lndk = 0;
    *lppk = 0;
    *lcnnd = 0;
    *lcndn = 0;
    *ldirt = 0;
    *ldecay = 0;
    *lres = 0;
    *ldou = 0;
    *lddrho = 0;
    *lnnrho = 0;
    *lnnom = 0;
    msum = 0;
    massrn[0] = 0;
/* COM: MSUM IS USED TO COUNT THE TOTAL NO. OF PARTICLES */
/*      IN PREVIOUS IRUN-1 RUNS */
/* KAON COUNTERS */
    for (il = 1; il <= 5; ++il) {
	kkk_1.tkaon[il - 1] = 0.f;
	for (is = 1; is <= 2000; ++is) {
	    kkk_1.ekaon[il + is * 7 - 1] = 0.f;
/* L1001: */
	}
/* L1002: */
    }
/* sp 12/19/00 */
    i__1 = run_1.num;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 150001; ++j) {
	    pe_1.propi[j + i__ * 150001 - 150002] = 1.f;
/* L1003: */
	}
/* L1004: */
    }
    for (i__ = 1; i__ <= 150001; ++i__) {
	fttemp[i__ - 1] = 0.f;
	for (irun = 1; irun <= 1; ++irun) {
	    ftpisv[i__ + irun * 150001 - 150002] = 0.f;
/* L1101: */
	}
/* L1102: */
    }
/* sp 12/19/00 end */
    *sp = 0.f;
/* antikaon counters */
    *akaon = 0.f;
    *sk = 0.f;
/* ----------------------------------------------------------------------- */
/*     LOOP OVER ALL PARALLEL RUNS */
/* bz11/17/98 */
/*      MASS=MASSPR+MASSTA */
    mass = 0;
/* bz11/17/98end */
    i__1 = run_1.num;
    for (irun = 1; irun <= i__1; ++irun) {
	nn_1.nnn = 0;
	msum += rr_1.massr[irun - 1];
/*     LOOP OVER ALL PSEUDOPARTICLES 1 IN THE SAME RUN */
	j10 = 2;
	if (*nt == *ntmax) {
	    j10 = 1;
	}

/* test off skips the check of energy conservation after each timestep: */
/*         enetot=0. */
/*         do ip=1,MASSR(IRUN) */
/*            if(e(ip).ne.0.or.lb(ip).eq.10022) enetot=enetot */
/*     1           +sqrt(p(1,ip)**2+p(2,ip)**2+p(3,ip)**2+e(ip)**2) */
/*         enddo */
/*         write(91,*) 'A:',nt,enetot,massr(irun),bimp */
	i__2 = rr_1.massr[irun];
	for (j1 = j10; j1 <= i__2; ++j1) {
	    i1 = j1 + msum;
/* E(I)=0 are for pions having been absorbed or photons which do not enter here: */
/* lin-4/2012 option of pi0 decays: */
/*            IF(E(I1).EQ.0.)GO TO 800 */
	    if (cc_1.e[i1 - 1] == 0.f) {
		goto L798;
	    }
/*     To include anti-(Delta,N*1440 and N*1535): */
/*          IF ((LB(I1) .LT. -13 .OR. LB(I1) .GT. 28) */
/*     1         .and.iabs(LB(I1)) .ne. 30 ) GOTO 800 */
/* lin-4/2012 option of pi0 decays: */
/*            IF (LB(I1) .LT. -45 .OR. LB(I1) .GT. 45) GOTO 800 */
	    if (ee_1.lb[i1 - 1] < -45 || ee_1.lb[i1 - 1] > 45) {
		goto L798;
	    }
	    x1 = aa_1.r__[i1 * 3 - 3];
	    y1 = aa_1.r__[i1 * 3 - 2];
	    z1 = aa_1.r__[i1 * 3 - 1];
	    leadng_1.px1 = bb_1.p[i1 * 3 - 3];
	    leadng_1.py1 = bb_1.p[i1 * 3 - 2];
	    leadng_1.pz1 = bb_1.p[i1 * 3 - 1];
	    leadng_1.em1 = cc_1.e[i1 - 1];
	    am1 = leadng_1.em1;
/* Computing 2nd power */
	    r__1 = leadng_1.em1;
/* Computing 2nd power */
	    r__2 = leadng_1.px1;
/* Computing 2nd power */
	    r__3 = leadng_1.py1;
/* Computing 2nd power */
	    r__4 = leadng_1.pz1;
	    leadng_1.e1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 
		    * r__4);
	    id1 = ee_1.id[i1 - 1];
	    leadng_1.lb1 = ee_1.lb[i1 - 1];
/*     generate k0short and k0long from K+ and K- at the last timestep: */
	    if (*nt == *ntmax && (leadng_1.lb1 == 21 || leadng_1.lb1 == 23)) {
		pk0 = ranart_(&rndf77_1.nseed);
		if (pk0 < .25f) {
		    ee_1.lb[i1 - 1] = 22;
		} else if (pk0 < .5f) {
		    ee_1.lb[i1 - 1] = 24;
		}
		leadng_1.lb1 = ee_1.lb[i1 - 1];
	    }
/* lin-8/07/02 these particles don't decay strongly, so skip decay routines: */
/*            IF( (lb1.ge.-2.and.lb1.le.5) .OR. lb1.eq.31 .OR. */
/*     &           (iabs(lb1).ge.14.and.iabs(lb1).le.24) .OR. */
/*     &           (iabs(lb1).ge.40.and.iabs(lb1).le.45) .or. */
/*     &           lb1.eq.31)GO TO 1 */
/*     only decay K0short when iksdcy=1: */
	    if (leadng_1.lb1 == 0 || leadng_1.lb1 == 25 || leadng_1.lb1 == 26 
		    || leadng_1.lb1 == 27 || leadng_1.lb1 == 28 || 
		    leadng_1.lb1 == 29 || abs(leadng_1.lb1) == 30 || abs(
		    leadng_1.lb1) >= 6 && abs(leadng_1.lb1) <= 13 || 
		    resdcy_1.iksdcy == 1 && leadng_1.lb1 == 24 || abs(
		    leadng_1.lb1) == 16 || phidcy_1.ipi0dcy == 1 && *nt == *
		    ntmax && leadng_1.lb1 == 4) {
/* lin-4/2012-above for option of pi0 decay: */
/*     &           .or.iabs(lb1).eq.16) then */
	    } else {
		goto L1;
	    }
/* IF I1 IS A RESONANCE, CHECK WHETHER IT DECAYS DURING THIS TIME STEP */
	    if (leadng_1.lb1 >= 25 && leadng_1.lb1 <= 27) {
		wid = .151f;
	    } else if (leadng_1.lb1 == 28) {
		wid = .00841f;
	    } else if (leadng_1.lb1 == 29) {
		wid = .00443f;
	    } else if (abs(leadng_1.lb1) == 30) {
		wid = .051f;
	    } else if (leadng_1.lb1 == 0) {
		wid = 1.18e-6f;
/*     to give K0short ct0=2.676cm: */
	    } else if (resdcy_1.iksdcy == 1 && leadng_1.lb1 == 24) {
		wid = 7.36e-15f;
/* lin-4/29/03 add Sigma0 decay to Lambda, ct0=2.22E-11m: */
	    } else if (abs(leadng_1.lb1) == 16) {
		wid = 8.87e-6f;
/* sp-07/25/01 test a1 resonance: */
/* c          ELSEIF(LB1.EQ.32) then */
/* c             WID=0.40 */
	    } else if (leadng_1.lb1 == 32) {
		wida1_(&leadng_1.em1, &rhomp, &wid, &input1_1.iseed);
	    } else if (abs(leadng_1.lb1) >= 6 && abs(leadng_1.lb1) <= 9) {
		wid = width_(&leadng_1.em1);
	    } else if (abs(leadng_1.lb1) == 10 || abs(leadng_1.lb1) == 11) {
		wid = w1440_(&leadng_1.em1);
	    } else if (abs(leadng_1.lb1) == 12 || abs(leadng_1.lb1) == 13) {
		wid = w1535_(&leadng_1.em1);
/* lin-4/2012 for option of pi0 decay: */
	    } else if (phidcy_1.ipi0dcy == 1 && *nt == *ntmax && leadng_1.lb1 
		    == 4) {
		wid = 7.85e-9f;
	    }
/* if it is the last time step, FORCE all resonance to strong-decay */
/* and go out of the loop */
	    if (*nt == *ntmax) {
		pdecay = 1.1f;
/* lin-5b/2008 forbid phi decay at the end of hadronic cascade: */
		if (phidcy_1.iphidcy == 0 && abs(leadng_1.lb1) == 29) {
		    pdecay = 0.f;
		}
/* test off clin-9/2012 forbid long-time decays (eta,omega,K*,Sigma0) */
/*     at the end of hadronic cascade to analyze freezeout time: */
/*             if(LB1.eq.0.or.LB1.eq.28.or.iabs(LB1).eq.30 */
/*     1            .or.iabs(LB1).eq.16) pdecay=0. */
	    } else {
		t0 = .19733f / wid;
		gfactr = leadng_1.e1 / leadng_1.em1;
		t0 *= gfactr;
		if (t0 > 0.f) {
		    pdecay = 1.f - exp(-input1_1.dt / t0);
		} else {
		    pdecay = 0.f;
		}
	    }
	    xdecay = ranart_(&rndf77_1.nseed);
/* c dilepton production from rho0, omega, phi decay */
/* c        if(lb1.eq.26 .or. lb1.eq.28 .or. lb1.eq.29) */
/* c     &   call dec_ceres(nt,ntmax,irun,i1) */
/* c */
	    if (xdecay < pdecay) {
/* lin-10/25/02 get rid of argument usage mismatch in rhocay(): */
		idecay = irun;
		leadng_1.tfnl = *nt * input1_1.dt;
/* lin-10/28/03 keep formation time of hadrons unformed at nt=ntmax-1: */
		if (*nt == *ntmax && ftmax_1.ftsv[i1 - 1] > (*ntmax - 1) * 
			input1_1.dt) {
		    leadng_1.tfnl = ftmax_1.ftsv[i1 - 1];
		}
		leadng_1.xfnl = x1;
		leadng_1.yfnl = y1;
		leadng_1.zfnl = z1;
/* use PYTHIA to perform decays of eta,rho,omega,phi,K*,(K0s) and Delta: */
		if (leadng_1.lb1 == 0 || leadng_1.lb1 == 25 || leadng_1.lb1 ==
			 26 || leadng_1.lb1 == 27 || leadng_1.lb1 == 28 || 
			leadng_1.lb1 == 29 || abs(leadng_1.lb1) == 30 || abs(
			leadng_1.lb1) >= 6 && abs(leadng_1.lb1) <= 9 || 
			resdcy_1.iksdcy == 1 && leadng_1.lb1 == 24 || abs(
			leadng_1.lb1) == 16 || phidcy_1.ipi0dcy == 1 && *nt ==
			 *ntmax && leadng_1.lb1 == 4) {
/* lin-4/2012 Above for option of pi0 decay: */
/*     &           .or.iabs(lb1).eq.16) then */
/*     previous rho decay performed in rhodecay(): */
/*                nnn=nnn+1 */
/*                call rhodecay(idecay,i1,nnn,iseed) */

/* test off record decays of phi,K*,Lambda(1520) resonances: */
/*                if(lb1.eq.29.or.iabs(lb1).eq.30) */
/*     1               write(18,112) 'decay',lb1,px1,py1,pz1,am1,nt */

/* lin-4/2012 option of pi0 decays: */
/*                call resdec(i1,nt,nnn,wid,idecay) */
		    resdec_(&i1, nt, &nn_1.nnn, &wid, &idecay, &c__0);
		    bb_1.p[i1 * 3 - 3] = leadng_1.px1n;
		    bb_1.p[i1 * 3 - 2] = leadng_1.py1n;
		    bb_1.p[i1 * 3 - 1] = leadng_1.pz1n;
/* lin-5/2008: */
		    dpert_1.dpertp[i1 - 1] = leadng_1.dp1n;
/*     add decay time to freezeout positions & time at the last timestep: */
		    if (*nt == *ntmax) {
			aa_1.r__[i1 * 3 - 3] = leadng_1.xfnl;
			aa_1.r__[i1 * 3 - 2] = leadng_1.yfnl;
			aa_1.r__[i1 * 3 - 1] = leadng_1.zfnl;
			tdecay_1.tfdcy[i1 - 1] = leadng_1.tfnl;
		    }

/* decay number for baryon resonance or L/S decay */
		    if (abs(leadng_1.lb1) >= 6 && abs(leadng_1.lb1) <= 9) {
			++(*ldecay);
		    }
/* for a1 decay */
/*             elseif(lb1.eq.32)then */
/*                NNN=NNN+1 */
/*                call a1decay(idecay,i1,nnn,iseed,rhomp) */
/* FOR N*(1440) */
		} else if (abs(leadng_1.lb1) == 10 || abs(leadng_1.lb1) == 11)
			 {
		    ++nn_1.nnn;
		    ++(*ldecay);
		    pnstar = 1.f;
		    if (cc_1.e[i1 - 1] > 1.22f) {
			pnstar = .6f;
		    }
		    if (ranart_(&rndf77_1.nseed) <= pnstar) {
/* (1) DECAY TO SINGLE PION+NUCLEON */
			decay_(&idecay, &i1, &nn_1.nnn, &input1_1.iseed, &wid,
				 nt);
		    } else {
/* (2) DECAY TO TWO PIONS + NUCLEON */
			decay2_(&idecay, &i1, &nn_1.nnn, &input1_1.iseed, &
				wid, nt);
			++nn_1.nnn;
		    }
/* for N*(1535) decay */
		} else if (abs(leadng_1.lb1) == 12 || abs(leadng_1.lb1) == 13)
			 {
		    ++nn_1.nnn;
		    decay_(&idecay, &i1, &nn_1.nnn, &input1_1.iseed, &wid, nt)
			    ;
		    ++(*ldecay);
		}

/* COM: AT HIGH ENERGIES WE USE VERY SHORT TIME STEPS, */
/*     IN ORDER TO TAKE INTO ACCOUNT THE FINITE FORMATIOM TIME, WE */
/*     DO NOT ALLOW PARTICLES FROM THE DECAY OF RESONANCE TO INTERACT */
/*     WITH OTHERS IN THE SAME TIME STEP. CHANGE 9000 TO REVERSE THIS */
/*     ASSUMPTION. EFFECTS OF THIS ASSUMPTION CAN BE STUDIED BY CHANGING */
/*     THE STATEMENT OF 9000. See notebook for discussions on effects of */
/*     changing statement 9000. */

/*     kaons from K* decay are converted to k0short (and k0long), */
/*     phi decay may produce rho, K0S or eta, N*(1535) decay may produce eta, */
/*     and these decay daughters need to decay again if at the last timestep: */
/*     (note: these daughters have been assigned to lb(i1) only, not to lpion) */
/*             if(nt.eq.ntmax.and.(lb1.eq.29.or.iabs(lb1).eq.30 */
/*     1            .iabs(lb1).eq.12.or.iabs(lb1).eq.13)) then */
		if (*nt == *ntmax) {
		    if (ee_1.lb[i1 - 1] == 25 || ee_1.lb[i1 - 1] == 26 || 
			    ee_1.lb[i1 - 1] == 27) {
			wid = .151f;
		    } else if (ee_1.lb[i1 - 1] == 0) {
			wid = 1.18e-6f;
		    } else if (ee_1.lb[i1 - 1] == 24 && resdcy_1.iksdcy == 1) 
			    {
/* lin-4/2012 corrected K0s decay width: */
/*                   wid=7.36e-17 */
			wid = 7.36e-15f;
/* lin-4/2012 option of pi0 decays: */
		    } else if (phidcy_1.ipi0dcy == 1 && ee_1.lb[i1 - 1] == 4) 
			    {
			wid = 7.85e-9f;
		    } else {
			goto L9000;
		    }
		    leadng_1.lb1 = ee_1.lb[i1 - 1];
		    leadng_1.px1 = bb_1.p[i1 * 3 - 3];
		    leadng_1.py1 = bb_1.p[i1 * 3 - 2];
		    leadng_1.pz1 = bb_1.p[i1 * 3 - 1];
		    leadng_1.em1 = cc_1.e[i1 - 1];
/* Computing 2nd power */
		    r__1 = leadng_1.em1;
/* Computing 2nd power */
		    r__2 = leadng_1.px1;
/* Computing 2nd power */
		    r__3 = leadng_1.py1;
/* Computing 2nd power */
		    r__4 = leadng_1.pz1;
		    leadng_1.e1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * 
			    r__3 + r__4 * r__4);
/* lin-4/2012 option of pi0 decays: */
/*                call resdec(i1,nt,nnn,wid,idecay) */
		    resdec_(&i1, nt, &nn_1.nnn, &wid, &idecay, &c__0);
		    bb_1.p[i1 * 3 - 3] = leadng_1.px1n;
		    bb_1.p[i1 * 3 - 2] = leadng_1.py1n;
		    bb_1.p[i1 * 3 - 1] = leadng_1.pz1n;
		    aa_1.r__[i1 * 3 - 3] = leadng_1.xfnl;
		    aa_1.r__[i1 * 3 - 2] = leadng_1.yfnl;
		    aa_1.r__[i1 * 3 - 1] = leadng_1.zfnl;
		    tdecay_1.tfdcy[i1 - 1] = leadng_1.tfnl;
/* lin-5/2008: */
		    dpert_1.dpertp[i1 - 1] = leadng_1.dp1n;
		}
/*     Decay daughter of the above decay in lb(i1) may be a pi0: */
		if (*nt == *ntmax && phidcy_1.ipi0dcy == 1 && ee_1.lb[i1 - 1] 
			== 4) {
		    wid = 7.85e-9f;
		    leadng_1.lb1 = ee_1.lb[i1 - 1];
		    leadng_1.px1 = bb_1.p[i1 * 3 - 3];
		    leadng_1.py1 = bb_1.p[i1 * 3 - 2];
		    leadng_1.pz1 = bb_1.p[i1 * 3 - 1];
		    leadng_1.em1 = cc_1.e[i1 - 1];
/* Computing 2nd power */
		    r__1 = leadng_1.em1;
/* Computing 2nd power */
		    r__2 = leadng_1.px1;
/* Computing 2nd power */
		    r__3 = leadng_1.py1;
/* Computing 2nd power */
		    r__4 = leadng_1.pz1;
		    leadng_1.e1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * 
			    r__3 + r__4 * r__4);
		    resdec_(&i1, nt, &nn_1.nnn, &wid, &idecay, &c__0);
		    bb_1.p[i1 * 3 - 3] = leadng_1.px1n;
		    bb_1.p[i1 * 3 - 2] = leadng_1.py1n;
		    bb_1.p[i1 * 3 - 1] = leadng_1.pz1n;
		    aa_1.r__[i1 * 3 - 3] = leadng_1.xfnl;
		    aa_1.r__[i1 * 3 - 2] = leadng_1.yfnl;
		    aa_1.r__[i1 * 3 - 1] = leadng_1.zfnl;
		    tdecay_1.tfdcy[i1 - 1] = leadng_1.tfnl;
		    dpert_1.dpertp[i1 - 1] = leadng_1.dp1n;
		}
/* negelecting the Pauli blocking at high energies */
/* lin-4/2012 option of pi0 decays: */
/* 9000        go to 800 */
L9000:
		goto L798;
	    }
/* LOOP OVER ALL PSEUDOPARTICLES 2 IN THE SAME RUN */
/* SAVE ALL THE COORDINATES FOR POSSIBLE CHANGE IN THE FOLLOWING COLLISION */
/* lin-4/2012 option of pi0 decays: */
/* 1        if(nt.eq.ntmax)go to 800 */
L1:
	    if (*nt == *ntmax) {
		goto L798;
	    }
	    x1 = aa_1.r__[i1 * 3 - 3];
	    y1 = aa_1.r__[i1 * 3 - 2];
	    z1 = aa_1.r__[i1 * 3 - 1];

	    i__3 = j1 - 1;
	    for (j2 = 1; j2 <= i__3; ++j2) {
		i2 = j2 + msum;
/* IF I2 IS A MESON BEING ABSORBED, THEN GO OUT OF THE LOOP */
		if (cc_1.e[i2 - 1] == 0.f) {
		    goto L600;
		}
/* lin-5/2008 in case the first particle is already destroyed: */
		if (cc_1.e[i1 - 1] == 0.f) {
		    goto L800;
		}
/* lin-4/2012 option of pi0 decays: */
		if (ee_1.lb[i2 - 1] < -45 || ee_1.lb[i2 - 1] > 45) {
		    goto L600;
		}
/* lin-7/26/03 improve speed */
		x2 = aa_1.r__[i2 * 3 - 3];
		y2 = aa_1.r__[i2 * 3 - 2];
		z2 = aa_1.r__[i2 * 3 - 1];
		dr0max = 5.f;
/* lin-9/2008 deuteron+nucleon elastic cross sections could reach ~2810mb: */
		ilb1 = (i__4 = ee_1.lb[i1 - 1], abs(i__4));
		ilb2 = (i__4 = ee_1.lb[i2 - 1], abs(i__4));
		if (ilb1 == 42 || ilb2 == 42) {
		    if (ilb1 >= 1 && ilb1 <= 2 || ilb1 >= 6 && ilb1 <= 13 || 
			    ilb2 >= 1 && ilb2 <= 2 || ilb2 >= 6 && ilb2 <= 13)
			     {
			if (ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1] > 0) {
			    dr0max = 10.f;
			}
		    }
		}

/* Computing 2nd power */
		r__1 = x1 - x2;
/* Computing 2nd power */
		r__2 = y1 - y2;
/* Computing 2nd power */
		r__3 = z1 - z2;
/* Computing 2nd power */
		r__4 = dr0max;
		if (r__1 * r__1 + r__2 * r__2 + r__3 * r__3 > r__4 * r__4) {
		    goto L600;
		}
		if (ee_1.id[i1 - 1] * ee_1.id[i2 - 1] == input1_1.iavoid) {
		    goto L400;
		}
		id1 = ee_1.id[i1 - 1];
		id2 = ee_1.id[i2 - 1];

		r__1 = x1 / gg_1.dx;
		ix1 = i_nint(&r__1);
		r__1 = y1 / gg_1.dy;
		iy1 = i_nint(&r__1);
		r__1 = z1 / gg_1.dz;
		iz1 = i_nint(&r__1);
		leadng_1.px1 = bb_1.p[i1 * 3 - 3];
		leadng_1.py1 = bb_1.p[i1 * 3 - 2];
		leadng_1.pz1 = bb_1.p[i1 * 3 - 1];
		leadng_1.em1 = cc_1.e[i1 - 1];
		am1 = leadng_1.em1;
		leadng_1.lb1 = ee_1.lb[i1 - 1];
/* Computing 2nd power */
		r__1 = leadng_1.em1;
/* Computing 2nd power */
		r__2 = leadng_1.px1;
/* Computing 2nd power */
		r__3 = leadng_1.py1;
/* Computing 2nd power */
		r__4 = leadng_1.pz1;
		leadng_1.e1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + 
			r__4 * r__4);
		r__1 = leadng_1.px1 / gg_1.dpx;
		ipx1 = i_nint(&r__1);
		r__1 = leadng_1.py1 / gg_1.dpy;
		ipy1 = i_nint(&r__1);
		r__1 = leadng_1.pz1 / gg_1.dpz;
		ipz1 = i_nint(&r__1);
		dpi_1.lb2 = ee_1.lb[i2 - 1];
		px2 = bb_1.p[i2 * 3 - 3];
		py2 = bb_1.p[i2 * 3 - 2];
		pz2 = bb_1.p[i2 * 3 - 1];
		dpi_1.em2 = cc_1.e[i2 - 1];
		am2 = dpi_1.em2;
		lb1i = ee_1.lb[i1 - 1];
		lb2i = ee_1.lb[i2 - 1];
		px1i = bb_1.p[i1 * 3 - 3];
		py1i = bb_1.p[i1 * 3 - 2];
		pz1i = bb_1.p[i1 * 3 - 1];
		em1i = cc_1.e[i1 - 1];
		px2i = bb_1.p[i2 * 3 - 3];
		py2i = bb_1.p[i2 * 3 - 2];
		pz2i = bb_1.p[i2 * 3 - 1];
		em2i = cc_1.e[i2 - 1];
/* lin-2/26/03 ctest off check energy conservation after each binary search: */
/* Computing 2nd power */
		r__1 = cc_1.e[i1 - 1];
/* Computing 2nd power */
		r__2 = bb_1.p[i1 * 3 - 3];
/* Computing 2nd power */
		r__3 = bb_1.p[i1 * 3 - 2];
/* Computing 2nd power */
		r__4 = bb_1.p[i1 * 3 - 1];
/* Computing 2nd power */
		r__5 = cc_1.e[i2 - 1];
/* Computing 2nd power */
		r__6 = bb_1.p[i2 * 3 - 3];
/* Computing 2nd power */
		r__7 = bb_1.p[i2 * 3 - 2];
/* Computing 2nd power */
		r__8 = bb_1.p[i2 * 3 - 1];
		eini = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * 
			r__4) + sqrt(r__5 * r__5 + r__6 * r__6 + r__7 * r__7 
			+ r__8 * r__8);
		pxini = bb_1.p[i1 * 3 - 3] + bb_1.p[i2 * 3 - 3];
		pyini = bb_1.p[i1 * 3 - 2] + bb_1.p[i2 * 3 - 2];
		pzini = bb_1.p[i1 * 3 - 1] + bb_1.p[i2 * 3 - 1];
		nnnini = nn_1.nnn;

/* lin-4/30/03 initialize value: */
		iblock = 0;

/* TO SAVE COMPUTING TIME we do the following */
/* (1) make a ROUGH estimate to see whether particle i2 will collide with */
/* particle I1, and (2) skip the particle pairs for which collisions are */
/* not modeled in the code. */
/* FOR MESON-BARYON AND MESON-MESON COLLISIONS, we use a maximum */
/* interaction distance DELTR0=2.6 */
/* for ppbar production from meson (pi rho omega) interactions: */

		deltr0 = 3.f;
		if (abs(leadng_1.lb1) >= 14 && abs(leadng_1.lb1) <= 17 || abs(
			leadng_1.lb1) >= 30 && abs(leadng_1.lb1) <= 45) {
		    deltr0 = 5.f;
		}
		if (abs(dpi_1.lb2) >= 14 && abs(dpi_1.lb2) <= 17 || abs(
			dpi_1.lb2) >= 30 && abs(dpi_1.lb2) <= 45) {
		    deltr0 = 5.f;
		}
		if (leadng_1.lb1 == 28 && dpi_1.lb2 == 28) {
		    deltr0 = 4.84f;
		}
/* lin-10/08/00 to include pi pi -> rho rho: */
		if (leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5 && (dpi_1.lb2 >= 3 
			&& dpi_1.lb2 <= 5)) {
/* Computing 2nd power */
		    r__1 = dpi_1.em2;
/* Computing 2nd power */
		    r__2 = px2;
/* Computing 2nd power */
		    r__3 = py2;
/* Computing 2nd power */
		    r__4 = pz2;
		    e2 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 *
			     r__4);
/* Computing 2nd power */
		    r__1 = leadng_1.e1 + e2;
/* Computing 2nd power */
		    r__2 = leadng_1.px1 + px2;
/* Computing 2nd power */
		    r__3 = leadng_1.py1 + py2;
/* Computing 2nd power */
		    r__4 = leadng_1.pz1 + pz2;
		    spipi = r__1 * r__1 - r__2 * r__2 - r__3 * r__3 - r__4 * 
			    r__4;
		    if (spipi >= 2.3715999999999999f) {
			deltr0 = 3.5f;
		    }
		}
/* khyperon */
		if (leadng_1.lb1 == 23 && (dpi_1.lb2 >= 14 && dpi_1.lb2 <= 17)
			) {
		    goto L3699;
		}
		if (dpi_1.lb2 == 23 && (leadng_1.lb1 >= 14 && leadng_1.lb1 <= 
			17)) {
		    goto L3699;
		}
/* K(K*) + Kbar(K*bar) scattering including */
/*     K(K*) + Kbar(K*bar) --> phi + pi(rho,omega) and pi pi(rho,omega) */
		if (leadng_1.lb1 == 21 && dpi_1.lb2 == 23) {
		    goto L3699;
		}
		if (dpi_1.lb2 == 21 && leadng_1.lb1 == 23) {
		    goto L3699;
		}
		if (leadng_1.lb1 == 30 && dpi_1.lb2 == 21) {
		    goto L3699;
		}
		if (dpi_1.lb2 == 30 && leadng_1.lb1 == 21) {
		    goto L3699;
		}
		if (leadng_1.lb1 == -30 && dpi_1.lb2 == 23) {
		    goto L3699;
		}
		if (dpi_1.lb2 == -30 && leadng_1.lb1 == 23) {
		    goto L3699;
		}
		if (leadng_1.lb1 == -30 && dpi_1.lb2 == 30) {
		    goto L3699;
		}
		if (dpi_1.lb2 == -30 && leadng_1.lb1 == 30) {
		    goto L3699;
		}

/* lin-12/15/00 */
/*     kaon+rho(omega,eta) collisions: */
		if (leadng_1.lb1 == 21 || leadng_1.lb1 == 23) {
		    if (dpi_1.lb2 == 0 || dpi_1.lb2 >= 25 && dpi_1.lb2 <= 28) 
			    {
			goto L3699;
		    }
		} else if (dpi_1.lb2 == 21 || dpi_1.lb2 == 23) {
		    if (leadng_1.lb1 == 0 || leadng_1.lb1 >= 25 && 
			    leadng_1.lb1 <= 28) {
			goto L3699;
		    }
		}
/* lin-8/14/02 K* (pi, rho, omega, eta) collisions: */
		if (abs(leadng_1.lb1) == 30 && (dpi_1.lb2 == 0 || dpi_1.lb2 >=
			 25 && dpi_1.lb2 <= 28 || dpi_1.lb2 >= 3 && dpi_1.lb2 
			<= 5)) {
		    goto L3699;
		} else if (abs(dpi_1.lb2) == 30 && (leadng_1.lb1 == 0 || 
			leadng_1.lb1 >= 25 && leadng_1.lb1 <= 28 || 
			leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5)) {
		    goto L3699;
/* lin-8/14/02-end */
/* K+/K*-bar + baryon/antibaryon collisions: */
		} else if (abs(leadng_1.lb1) == 30 && (abs(dpi_1.lb2) == 1 || 
			abs(dpi_1.lb2) == 2 || abs(dpi_1.lb2) >= 6 && abs(
			dpi_1.lb2) <= 13)) {
		    goto L3699;
		}
		if (abs(dpi_1.lb2) == 30 && (abs(leadng_1.lb1) == 1 || abs(
			leadng_1.lb1) == 2 || abs(leadng_1.lb1) >= 6 && abs(
			leadng_1.lb1) <= 13)) {
		    goto L3699;
		}
/* K^+ baryons and antibaryons: */
/* ** K+ + B-bar  --> La(Si)-bar + pi */
/* K^- and antibaryons, note K^- and baryons are included in newka(): */
/* note that we fail to satisfy charge conjugation for these cross sections: */
		if ((leadng_1.lb1 == 23 || leadng_1.lb1 == 21) && (abs(
			dpi_1.lb2) == 1 || abs(dpi_1.lb2) == 2 || abs(
			dpi_1.lb2) >= 6 && abs(dpi_1.lb2) <= 13)) {
		    goto L3699;
		} else if ((dpi_1.lb2 == 23 || dpi_1.lb2 == 21) && (abs(
			leadng_1.lb1) == 1 || abs(leadng_1.lb1) == 2 || abs(
			leadng_1.lb1) >= 6 && abs(leadng_1.lb1) <= 13)) {
		    goto L3699;
		}

/* For anti-nucleons annihilations: */
/* Assumptions: */
/* (1) for collisions involving a p_bar or n_bar, */
/* we allow only collisions between a p_bar and a baryon or a baryon */
/* resonance (as well as a n_bar and a baryon or a baryon resonance), */
/* we skip all other reactions involving a p_bar or n_bar, */
/* such as collisions between p_bar (n_bar) and mesons, */
/* and collisions between two p_bar's (n_bar's). */
/* (2) we introduce a new parameter rppmax: the maximum interaction */
/* distance to make the quick collision check,rppmax=3.57 fm */
/* corresponding to a cutoff of annihilation xsection= 400mb which is */
/* also used consistently in the actual annihilation xsection to be */
/* used in the following as given in the subroutine xppbar(srt) */
		rppmax = 3.57f;
/* anti-baryon on baryons */
		if ((leadng_1.lb1 == -1 || leadng_1.lb1 == -2 || leadng_1.lb1 
			>= -13 && leadng_1.lb1 <= -6) && (dpi_1.lb2 == 1 || 
			dpi_1.lb2 == 2 || dpi_1.lb2 >= 6 && dpi_1.lb2 <= 13)) 
			{
		    deltr0 = rppmax;
		    goto L2699;
		} else if ((dpi_1.lb2 == -1 || dpi_1.lb2 == -2 || dpi_1.lb2 >=
			 -13 && dpi_1.lb2 <= -6) && (leadng_1.lb1 == 1 || 
			leadng_1.lb1 == 2 || leadng_1.lb1 >= 6 && 
			leadng_1.lb1 <= 13)) {
		    deltr0 = rppmax;
		    goto L2699;
		}
/* *  ((anti) lambda, cascade, omega  should not be rejected) */
		if (abs(leadng_1.lb1) >= 14 && abs(leadng_1.lb1) <= 17 || abs(
			dpi_1.lb2) >= 14 && abs(dpi_1.lb2) <= 17) {
		    goto L3699;
		}

/* lin-9/2008 maximum sigma~2810mb for deuteron+nucleon elastic collisions: */
		if (abs(leadng_1.lb1) == 42 || abs(dpi_1.lb2) == 42) {
		    ilb1 = abs(leadng_1.lb1);
		    ilb2 = abs(dpi_1.lb2);
		    if (ilb1 >= 1 && ilb1 <= 2 || ilb1 >= 6 && ilb1 <= 13 || 
			    ilb2 >= 1 && ilb2 <= 2 || ilb2 >= 6 && ilb2 <= 13)
			     {
			if (leadng_1.lb1 * dpi_1.lb2 > 0) {
			    deltr0 = 9.5f;
			}
		    }
		}

		if (abs(leadng_1.lb1) >= 40 && abs(leadng_1.lb1) <= 45 || abs(
			dpi_1.lb2) >= 40 && abs(dpi_1.lb2) <= 45) {
		    goto L3699;
		}

/* * phi channel --> elastic + inelastic scatt. */
		if (leadng_1.lb1 == 29 && (dpi_1.lb2 >= 1 && dpi_1.lb2 <= 13 
			|| dpi_1.lb2 >= 21 && dpi_1.lb2 <= 28 || abs(
			dpi_1.lb2) == 30) || dpi_1.lb2 == 29 && (leadng_1.lb1 
			>= 1 && leadng_1.lb1 <= 13 || leadng_1.lb1 >= 21 && 
			leadng_1.lb1 <= 28 || abs(leadng_1.lb1) == 30)) {
		    deltr0 = 3.f;
		    goto L3699;
		}

/*  La/Si, Cas, Om (bar)-meson elastic colln */
/* pion vs. La & Ca (bar) coll. are treated in resp. subroutines */
/* SKIP all other K* RESCATTERINGS */
		if (abs(leadng_1.lb1) == 30 || abs(dpi_1.lb2) == 30) {
		    goto L400;
		}
/* SKIP KAON(+) RESCATTERINGS WITH particles other than pions and baryons */
		if (leadng_1.lb1 == 23 && (dpi_1.lb2 < 1 || dpi_1.lb2 > 17)) {
		    goto L400;
		}
		if (dpi_1.lb2 == 23 && (leadng_1.lb1 < 1 || leadng_1.lb1 > 17)
			) {
		    goto L400;
		}

/* anti-baryon proccess: B-bar+M, N-bar+R-bar, N-bar+N-bar, R-bar+R-bar */
/*  R = (D,N*) */
		if (leadng_1.lb1 <= -1 && leadng_1.lb1 >= -13 && (dpi_1.lb2 ==
			 0 || dpi_1.lb2 >= 3 && dpi_1.lb2 <= 5 || dpi_1.lb2 >=
			 25 && dpi_1.lb2 <= 28) || dpi_1.lb2 <= -1 && 
			dpi_1.lb2 >= -13 && (leadng_1.lb1 == 0 || 
			leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5 || 
			leadng_1.lb1 >= 25 && leadng_1.lb1 <= 28)) {
		} else if ((leadng_1.lb1 == -1 || leadng_1.lb1 == -2) && (
			dpi_1.lb2 < -5 && dpi_1.lb2 >= -13) || (dpi_1.lb2 == 
			-1 || dpi_1.lb2 == -2) && (leadng_1.lb1 < -5 && 
			leadng_1.lb1 >= -13)) {
		} else if ((leadng_1.lb1 == -1 || leadng_1.lb1 == -2) && (
			dpi_1.lb2 == -1 || dpi_1.lb2 == -2)) {
		} else if (leadng_1.lb1 < -5 && leadng_1.lb1 >= -13 && (
			dpi_1.lb2 < -5 && dpi_1.lb2 >= -13)) {
/*        elseif((lb1.lt.0).or.(lb2.lt.0)) then */
/*         go to 400 */
		}
L2699:
/* for baryon-baryon collisions */
		if (leadng_1.lb1 == 1 || leadng_1.lb1 == 2 || leadng_1.lb1 >= 
			6 && leadng_1.lb1 <= 17) {
		    if (dpi_1.lb2 == 1 || dpi_1.lb2 == 2 || dpi_1.lb2 >= 6 && 
			    dpi_1.lb2 <= 17) {
			deltr0 = 2.f;
		    }
		}

L3699:
/* Computing 2nd power */
		r__1 = x1 - x2;
/* Computing 2nd power */
		r__2 = y1 - y2;
/* Computing 2nd power */
		r__3 = z1 - z2;
		rsqare = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
/* Computing 2nd power */
		r__1 = deltr0;
		if (rsqare > r__1 * r__1) {
		    goto L400;
		}
/* NOW PARTICLES ARE CLOSE ENOUGH TO EACH OTHER ! */
/* KEEP ALL COORDINATES FOR POSSIBLE PHASE SPACE CHANGE */
		r__1 = x2 / gg_1.dx;
		ix2 = i_nint(&r__1);
		r__1 = y2 / gg_1.dy;
		iy2 = i_nint(&r__1);
		r__1 = z2 / gg_1.dz;
		iz2 = i_nint(&r__1);
		r__1 = px2 / gg_1.dpx;
		ipx2 = i_nint(&r__1);
		r__1 = py2 / gg_1.dpy;
		ipy2 = i_nint(&r__1);
		r__1 = pz2 / gg_1.dpz;
		ipz2 = i_nint(&r__1);
/* FIND MOMENTA OF PARTICLES IN THE CMS OF THE TWO COLLIDING PARTICLES */
/* AND THE CMS ENERGY SRT */
		cms_(&i1, &i2, &pcx, &pcy, &pcz, &srt);
/* lin-7/26/03 improve speed */
		drmax = dr0max;
		distc0_(&drmax, &deltr0, &input1_1.dt, &ifirst, &pcx, &pcy, &
			pcz, &x1, &y1, &z1, &leadng_1.px1, &leadng_1.py1, &
			leadng_1.pz1, &leadng_1.em1, &x2, &y2, &z2, &px2, &
			py2, &pz2, &dpi_1.em2);
		if (ifirst == -1) {
		    goto L400;
		}
		r__1 = srt / .04f;
		iss = i_nint(&r__1);
/* lin-4/2008 use last bin if ISS is out of EKAON's upper bound of 2000: */
		if (iss > 2000) {
		    iss = 2000;
		}
/* Sort collisions */

/* lin-8/2008 Deuteron+Meson->B+B; */
/*     meson=(pi,rho,omega,eta), B=(n,p,Delta,N*1440,N*1535): */
		if (abs(leadng_1.lb1) == 42 || abs(dpi_1.lb2) == 42) {
		    ilb1 = abs(leadng_1.lb1);
		    ilb2 = abs(dpi_1.lb2);
		    if (leadng_1.lb1 == 0 || leadng_1.lb1 >= 3 && 
			    leadng_1.lb1 <= 5 || leadng_1.lb1 >= 25 && 
			    leadng_1.lb1 <= 28 || dpi_1.lb2 == 0 || dpi_1.lb2 
			    >= 3 && dpi_1.lb2 <= 5 || dpi_1.lb2 >= 25 && 
			    dpi_1.lb2 <= 28) {
			goto L505;
/* lin-9/2008 Deuteron+Baryon or antiDeuteron+antiBaryon elastic collisions: */
		    } else if ((ilb1 >= 1 && ilb1 <= 2 || ilb1 >= 6 && ilb1 <=
			     13 || ilb2 >= 1 && ilb2 <= 2 || ilb2 >= 6 && 
			    ilb2 <= 13) && leadng_1.lb1 * dpi_1.lb2 > 0) {
			goto L506;
		    } else {
			goto L400;
		    }
		}

/* K+ + (N,N*,D)-bar --> L/S-bar + pi */
		if ((leadng_1.lb1 == 23 || leadng_1.lb1 == 30) && (dpi_1.lb2 
			== -1 || dpi_1.lb2 == -2 || dpi_1.lb2 >= -13 && 
			dpi_1.lb2 <= -6) || (dpi_1.lb2 == 23 || dpi_1.lb2 == 
			30) && (leadng_1.lb1 == -1 || leadng_1.lb1 == -2 || 
			leadng_1.lb1 >= -13 && leadng_1.lb1 <= -6)) {
		    bmass = .938f;
		    if (srt <= bmass + .498f) {
			pkaon = 0.f;
		    } else {
/* Computing 2nd power */
			r__2 = srt;
/* Computing 2nd power */
			r__3 = bmass;
/* Computing 2nd power */
			r__1 = (r__2 * r__2 - (r__3 * r__3 + .248004f)) / 2.f 
				/ bmass;
			pkaon = sqrt(r__1 * r__1 - .248004f);
		    }
/* lin-10/31/02 cross sections are isospin-averaged, same as those in newka */
/*     for K- + (N,N*,D) --> L/S + pi: */
		    sigela = (akpel_(&pkaon) + aknel_(&pkaon)) * .5f;
		    sigsgm = akpsgm_(&pkaon) * 1.5f + aknsgm_(&pkaon);
		    sig = sigela + sigsgm + akplam_(&pkaon);
		    if (sig > 1e-7f) {
/*     ! K+ + N-bar reactions */
			icase = 3;
			brel = sigela / sig;
			brsgm = sigsgm / sig;
			brsig = sig;
			nchrg = 1;
			goto L3555;
		    }
		    goto L400;
		}


/*  meson + hyperon-bar -> K+ + N-bar */
		if (leadng_1.lb1 >= -17 && leadng_1.lb1 <= -14 && (dpi_1.lb2 
			>= 3 && dpi_1.lb2 <= 5) || dpi_1.lb2 >= -17 && 
			dpi_1.lb2 <= -14 && (leadng_1.lb1 >= 3 && 
			leadng_1.lb1 <= 5)) {
		    nchrg = -100;
/* *       first classify the reactions due to total charge. */
		    if (leadng_1.lb1 == -15 && (dpi_1.lb2 == 5 || dpi_1.lb2 ==
			     27) || dpi_1.lb2 == -15 && (leadng_1.lb1 == 5 || 
			    leadng_1.lb1 == 27)) {
			nchrg = -2;
/*     ! D-(bar) */
			bmass = 1.232f;
			goto L110;
		    }
		    if (leadng_1.lb1 == -15 && (dpi_1.lb2 == 0 || dpi_1.lb2 ==
			     4 || dpi_1.lb2 == 26 || dpi_1.lb2 == 28) || 
			    dpi_1.lb2 == -15 && (leadng_1.lb1 == 0 || 
			    leadng_1.lb1 == 4 || leadng_1.lb1 == 26 || 
			    leadng_1.lb1 == 28) || (leadng_1.lb1 == -14 || 
			    leadng_1.lb1 == -16) && (dpi_1.lb2 == 5 || 
			    dpi_1.lb2 == 27) || (dpi_1.lb2 == -14 || 
			    dpi_1.lb2 == -16) && (leadng_1.lb1 == 5 || 
			    leadng_1.lb1 == 27)) {
			nchrg = -1;
/*     ! n-bar */
			bmass = .938f;
			goto L110;
		    }
		    if (leadng_1.lb1 == -15 && (dpi_1.lb2 == 3 || dpi_1.lb2 ==
			     25) || dpi_1.lb2 == -15 && (leadng_1.lb1 == 3 || 
			    leadng_1.lb1 == 25) || leadng_1.lb1 == -17 && (
			    dpi_1.lb2 == 5 || dpi_1.lb2 == 27) || dpi_1.lb2 ==
			     -17 && (leadng_1.lb1 == 5 || leadng_1.lb1 == 27) 
			    || (leadng_1.lb1 == -14 || leadng_1.lb1 == -16) &&
			     (dpi_1.lb2 == 0 || dpi_1.lb2 == 4 || dpi_1.lb2 ==
			     26 || dpi_1.lb2 == 28) || (dpi_1.lb2 == -14 || 
			    dpi_1.lb2 == -16) && (leadng_1.lb1 == 0 || 
			    leadng_1.lb1 == 4 || leadng_1.lb1 == 26 || 
			    leadng_1.lb1 == 28)) {
			nchrg = 0;
/*     ! p-bar */
			bmass = .938f;
			goto L110;
		    }
		    if (leadng_1.lb1 == -17 && (dpi_1.lb2 == 0 || dpi_1.lb2 ==
			     4 || dpi_1.lb2 == 26 || dpi_1.lb2 == 28) || 
			    dpi_1.lb2 == -17 && (leadng_1.lb1 == 0 || 
			    leadng_1.lb1 == 4 || leadng_1.lb1 == 26 || 
			    leadng_1.lb1 == 28) || (leadng_1.lb1 == -14 || 
			    leadng_1.lb1 == -16) && (dpi_1.lb2 == 3 || 
			    dpi_1.lb2 == 25) || (dpi_1.lb2 == -14 || 
			    dpi_1.lb2 == -16) && (leadng_1.lb1 == 3 || 
			    leadng_1.lb1 == 25)) {
			nchrg = 1;
/*     ! D++(bar) */
			bmass = 1.232f;
		    }

/* 110     if(nchrg.ne.-100.and.srt.ge.(aka+bmass))then !! for elastic */
L110:
		    sig = 0.f;
/* !! for elastic */
		    if (nchrg != -100 && srt >= bmass + .498f) {
/* c110        if(nchrg.eq.-100.or.srt.lt.(aka+bmass)) go to 400 */
/*             ! PI + La(Si)-bar => K+ + N-bar reactions */
			icase = 4;
/* c       pkaon=sqrt(((srt**2-(aka**2+bmass**2))/2./bmass)**2-aka**2) */
/* Computing 2nd power */
			r__2 = srt;
/* Computing 2nd power */
			r__1 = (r__2 * r__2 - 1.1278479999999997f) / 2.f / 
				.938f;
			pkaon = sqrt(r__1 * r__1 - .248004f);
/* ! lambda-bar + Pi */
			if (leadng_1.lb1 == -14 || dpi_1.lb2 == -14) {
			    if (nchrg >= 0) {
				sigma0 = akplam_(&pkaon);
			    }
			    if (nchrg < 0) {
				sigma0 = aknlam_(&pkaon);
			    }
/*                ! sigma-bar + pi */
			} else {
/* !K-p or K-D++ */
			    if (nchrg >= 0) {
				sigma0 = akpsgm_(&pkaon);
			    }
/* !K-n or K-D- */
			    if (nchrg < 0) {
				sigma0 = aknsgm_(&pkaon);
			    }
			    sigma0 = akpsgm_(&pkaon) * 1.5f + aknsgm_(&pkaon);
			}
/* Computing 2nd power */
			r__1 = srt;
/* Computing 2nd power */
			r__2 = bmass + .498f;
/* Computing 2nd power */
			r__3 = srt;
/* Computing 2nd power */
			r__4 = .498f - bmass;
/* Computing 2nd power */
			r__5 = srt;
/* Computing 2nd power */
			r__6 = leadng_1.em1 + dpi_1.em2;
/* Computing 2nd power */
			r__7 = srt;
/* Computing 2nd power */
			r__8 = leadng_1.em1 - dpi_1.em2;
			sig = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - 
				r__4 * r__4) / (r__5 * r__5 - r__6 * r__6) / (
				r__7 * r__7 - r__8 * r__8) * sigma0;
/* ! K0barD++, K-D- */
			if (nchrg == -2 || nchrg == 2) {
			    sig *= 2.f;
			}
/* *     the factor 2 comes from spin of delta, which is 3/2 */
/* *     detailed balance. copy from Page 423 of N.P. A614 1997 */
			if (leadng_1.lb1 == -14 || dpi_1.lb2 == -14) {
			    sig *= 1.3333333333333333f;
			} else if (nchrg == -2 || nchrg == 2) {
			    sig *= .88888888888888884f;
			} else {
			    sig *= .44444444444444442f;
			}
/* c        brel=0. */
/* c        brsgm=0. */
/* c        brsig = sig */
/* c          if(sig.lt.1.e-7) go to 400 */
/* - */
		    }
/*                ! PI + La(Si)-bar => elastic included */
		    icase = 4;
		    sigela = 10.f;
		    sig += sigela;
		    brel = sigela / sig;
		    brsgm = 0.f;
		    brsig = sig;
/* - */
		    goto L3555;
		}
/* * MULTISTRANGE PARTICLE (Cas,Omega -bar) PRODUCTION - (NON)PERTURBATIVE */
/* K-/K*0bar + La/Si --> cascade + pi/eta */
		if ((leadng_1.lb1 == 21 || leadng_1.lb1 == -30) && (dpi_1.lb2 
			>= 14 && dpi_1.lb2 <= 17) || (dpi_1.lb2 == 21 || 
			dpi_1.lb2 == -30) && (leadng_1.lb1 >= 14 && 
			leadng_1.lb1 <= 17)) {
		    kp = 0;
		    goto L3455;
		}
/* K+/K*0 + La/Si(bar) --> cascade-bar + pi/eta */
		if ((leadng_1.lb1 == 23 || leadng_1.lb1 == 30) && (dpi_1.lb2 
			<= -14 && dpi_1.lb2 >= -17) || (dpi_1.lb2 == 23 || 
			dpi_1.lb2 == 30) && (leadng_1.lb1 <= -14 && 
			leadng_1.lb1 >= -17)) {
		    kp = 1;
		    goto L3455;
		}
/* K-/K*0bar + cascade --> omega + pi */
		if ((leadng_1.lb1 == 21 || leadng_1.lb1 == -30) && (dpi_1.lb2 
			== 40 || dpi_1.lb2 == 41) || (dpi_1.lb2 == 21 || 
			dpi_1.lb2 == -30) && (leadng_1.lb1 == 40 || 
			leadng_1.lb1 == 41)) {
		    kp = 0;
		    goto L3455;
		}
/* K+/K*0 + cascade-bar --> omega-bar + pi */
		if ((leadng_1.lb1 == 23 || leadng_1.lb1 == 30) && (dpi_1.lb2 
			== -40 || dpi_1.lb2 == -41) || (dpi_1.lb2 == 23 || 
			dpi_1.lb2 == 30) && (leadng_1.lb1 == -40 || 
			leadng_1.lb1 == -41)) {
		    kp = 1;
		    goto L3455;
		}
/* Omega + Omega --> Di-Omega + photon(eta) */
/* c        if( lb1.eq.45.and.lb2.eq.45 ) go to 3455 */
/* annhilation of cascade(bar), omega(bar) */
		kp = 3;
/* K- + L/S <-- cascade(bar) + pi/eta */
		if ((leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5 || leadng_1.lb1 ==
			 0) && (abs(dpi_1.lb2) == 40 || abs(dpi_1.lb2) == 41) 
			|| (dpi_1.lb2 >= 3 && dpi_1.lb2 <= 5 || dpi_1.lb2 == 
			0) && (abs(leadng_1.lb1) == 40 || abs(leadng_1.lb1) ==
			 41)) {
		    goto L3455;
		}
/* K- + cascade(bar) <-- omega(bar) + pi */
/*         if(  (lb1.eq.0.and.iabs(lb2).eq.45) */
/*    &       .OR. (lb2.eq.0.and.iabs(lb1).eq.45) )go to 3455 */
		if (leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5 && abs(dpi_1.lb2) 
			== 45 || dpi_1.lb2 >= 3 && dpi_1.lb2 <= 5 && abs(
			leadng_1.lb1) == 45) {
		    goto L3455;
		}

/* **  MULTISTRANGE PARTICLE PRODUCTION  (END) */
/* * K+ + La(Si) --> Meson + B */
		if (leadng_1.lb1 == 23 && (dpi_1.lb2 >= 14 && dpi_1.lb2 <= 17)
			) {
		    goto L5699;
		}
		if (dpi_1.lb2 == 23 && (leadng_1.lb1 >= 14 && leadng_1.lb1 <= 
			17)) {
		    goto L5699;
		}
/* * K- + La(Si)-bar --> Meson + B-bar */
		if (leadng_1.lb1 == 21 && (dpi_1.lb2 >= -17 && dpi_1.lb2 <= 
			-14)) {
		    goto L5699;
		}
		if (dpi_1.lb2 == 21 && (leadng_1.lb1 >= -17 && leadng_1.lb1 <=
			 -14)) {
		    goto L5699;
		}
/* La/Si-bar + B --> pi + K+ */
		if ((leadng_1.lb1 == 1 || leadng_1.lb1 == 2 || leadng_1.lb1 >=
			 6 && leadng_1.lb1 <= 13) && (dpi_1.lb2 >= -17 && 
			dpi_1.lb2 <= -14) || (dpi_1.lb2 == 1 || dpi_1.lb2 == 
			2 || dpi_1.lb2 >= 6 && dpi_1.lb2 <= 13) && (
			leadng_1.lb1 >= -17 && leadng_1.lb1 <= -14)) {
		    goto L5999;
		}
/* La/Si + B-bar --> pi + K- */
		if ((leadng_1.lb1 == -1 || leadng_1.lb1 == -2 || leadng_1.lb1 
			<= -6 && leadng_1.lb1 >= -13) && (dpi_1.lb2 >= 14 && 
			dpi_1.lb2 <= 17) || (dpi_1.lb2 == -1 || dpi_1.lb2 == 
			-2 || dpi_1.lb2 <= -6 && dpi_1.lb2 >= -13) && (
			leadng_1.lb1 >= 14 && leadng_1.lb1 <= 17)) {
		    goto L5999;
		}


/* K(K*) + Kbar(K*bar) --> phi + pi(rho,omega), M + M (M=pi,rho,omega,eta) */
		if (leadng_1.lb1 == 21 && dpi_1.lb2 == 23) {
		    goto L8699;
		}
		if (dpi_1.lb2 == 21 && leadng_1.lb1 == 23) {
		    goto L8699;
		}
		if (leadng_1.lb1 == 30 && dpi_1.lb2 == 21) {
		    goto L8699;
		}
		if (dpi_1.lb2 == 30 && leadng_1.lb1 == 21) {
		    goto L8699;
		}
		if (leadng_1.lb1 == -30 && dpi_1.lb2 == 23) {
		    goto L8699;
		}
		if (dpi_1.lb2 == -30 && leadng_1.lb1 == 23) {
		    goto L8699;
		}
		if (leadng_1.lb1 == -30 && dpi_1.lb2 == 30) {
		    goto L8699;
		}
		if (dpi_1.lb2 == -30 && leadng_1.lb1 == 30) {
		    goto L8699;
		}
/* * (K,K*)-bar + rho(omega) --> phi +(K,K*)-bar, piK and elastic */
		if ((leadng_1.lb1 == 23 || leadng_1.lb1 == 21 || abs(
			leadng_1.lb1) == 30) && (dpi_1.lb2 >= 25 && dpi_1.lb2 
			<= 28) || (dpi_1.lb2 == 23 || dpi_1.lb2 == 21 || abs(
			dpi_1.lb2) == 30) && (leadng_1.lb1 >= 25 && 
			leadng_1.lb1 <= 28)) {
		    goto L8799;
		}

/* * K*(-bar) + pi --> phi + (K,K*)-bar */
		if (abs(leadng_1.lb1) == 30 && (dpi_1.lb2 >= 3 && dpi_1.lb2 <=
			 5) || abs(dpi_1.lb2) == 30 && (leadng_1.lb1 >= 3 && 
			leadng_1.lb1 <= 5)) {
		    goto L8799;
		}


/* * phi + N --> pi+N(D),  rho+N(D),  K+ +La */
/* * phi + D --> pi+N(D),  rho+N(D) */
		if (leadng_1.lb1 == 29 && (dpi_1.lb2 == 1 || dpi_1.lb2 == 2 ||
			 dpi_1.lb2 >= 6 && dpi_1.lb2 <= 9) || dpi_1.lb2 == 29 
			&& (leadng_1.lb1 == 1 || leadng_1.lb1 == 2 || 
			leadng_1.lb1 >= 6 && leadng_1.lb1 <= 9)) {
		    goto L7222;
		}

/* * phi + (pi,rho,ome,K,K*-bar) --> K+K, K+K*, K*+K*, (pi,rho,omega)+(K,K*-bar) */
		if (leadng_1.lb1 == 29 && (dpi_1.lb2 >= 3 && dpi_1.lb2 <= 5 ||
			 dpi_1.lb2 >= 21 && dpi_1.lb2 <= 28 || abs(dpi_1.lb2) 
			== 30) || dpi_1.lb2 == 29 && (leadng_1.lb1 >= 3 && 
			leadng_1.lb1 <= 5 || leadng_1.lb1 >= 21 && 
			leadng_1.lb1 <= 28 || abs(leadng_1.lb1) == 30)) {
		    goto L7444;
		}


/* La/Si, Cas, Om (bar)-(rho,omega,phi) elastic colln */
/* pion vs. La, Ca, Omega-(bar) elastic coll. treated in resp. subroutines */
		if ((abs(leadng_1.lb1) >= 14 && abs(leadng_1.lb1) <= 17 || 
			abs(leadng_1.lb1) >= 40) && (dpi_1.lb2 >= 25 && 
			dpi_1.lb2 <= 29 || dpi_1.lb2 == 0)) {
		    goto L888;
		}
		if ((abs(dpi_1.lb2) >= 14 && abs(dpi_1.lb2) <= 17 || abs(
			dpi_1.lb2) >= 40) && (leadng_1.lb1 >= 25 && 
			leadng_1.lb1 <= 29 || leadng_1.lb1 == 0)) {
		    goto L888;
		}

/* K+/K* (N,R)  OR   K-/K*- (N,R)-bar  elastic scatt */
		if ((leadng_1.lb1 == 23 || leadng_1.lb1 == 30) && (dpi_1.lb2 
			== 1 || dpi_1.lb2 == 2 || dpi_1.lb2 >= 6 && dpi_1.lb2 
			<= 13) || (dpi_1.lb2 == 23 || dpi_1.lb2 == 30) && (
			leadng_1.lb1 == 1 || leadng_1.lb1 == 2 || 
			leadng_1.lb1 >= 6 && leadng_1.lb1 <= 13)) {
		    goto L888;
		}
		if ((leadng_1.lb1 == 21 || leadng_1.lb1 == -30) && (dpi_1.lb2 
			== -1 || dpi_1.lb2 == -2 || dpi_1.lb2 >= -13 && 
			dpi_1.lb2 <= -6) || (dpi_1.lb2 == 21 || dpi_1.lb2 == 
			-30) && (leadng_1.lb1 == -1 || leadng_1.lb1 == -2 || 
			leadng_1.lb1 >= -13 && leadng_1.lb1 <= -6)) {
		    goto L888;
		}

/* L/S-baryon elastic collision */
		if (leadng_1.lb1 >= 14 && leadng_1.lb1 <= 17 && (dpi_1.lb2 >= 
			6 && dpi_1.lb2 <= 13) || dpi_1.lb2 >= 14 && dpi_1.lb2 
			<= 17 && (leadng_1.lb1 >= 6 && leadng_1.lb1 <= 13)) {
		    goto L7799;
		}
		if (leadng_1.lb1 <= -14 && leadng_1.lb1 >= -17 && (dpi_1.lb2 
			<= -6 && dpi_1.lb2 >= -13) || dpi_1.lb2 <= -14 && 
			dpi_1.lb2 >= -17 && (leadng_1.lb1 <= -6 && 
			leadng_1.lb1 >= -13)) {
		    goto L7799;
		}

/* skip other collns with perturbative particles or hyperon-bar */
		if (abs(leadng_1.lb1) >= 40 || abs(dpi_1.lb2) >= 40 || 
			leadng_1.lb1 <= -14 && leadng_1.lb1 >= -17 || 
			dpi_1.lb2 <= -14 && dpi_1.lb2 >= -17) {
		    goto L400;
		}


/* anti-baryon on baryon resonaces */
		if ((leadng_1.lb1 == -1 || leadng_1.lb1 == -2 || leadng_1.lb1 
			>= -13 && leadng_1.lb1 <= -6) && (dpi_1.lb2 == 1 || 
			dpi_1.lb2 == 2 || dpi_1.lb2 >= 6 && dpi_1.lb2 <= 13)) 
			{
		    goto L2799;
		} else if ((dpi_1.lb2 == -1 || dpi_1.lb2 == -2 || dpi_1.lb2 >=
			 -13 && dpi_1.lb2 <= -6) && (leadng_1.lb1 == 1 || 
			leadng_1.lb1 == 2 || leadng_1.lb1 >= 6 && 
			leadng_1.lb1 <= 13)) {
		    goto L2799;
		}

/* lin-10/25/02 get rid of argument usage mismatch in newka(): */
		inewka = irun;
/*        call newka(icase,irun,iseed,dt,nt, */
/* lin-5/01/03 set iblock value in art1f.f, necessary for resonance studies: */
/*        call newka(icase,inewka,iseed,dt,nt, */
/*     &                  ictrl,i1,i2,srt,pcx,pcy,pcz) */
		newka_(&icase, &inewka, &input1_1.iseed, &input1_1.dt, nt, &
			ictrl, &i1, &i2, &srt, &pcx, &pcy, &pcz, &iblock);
/* lin-10/25/02-end */
		if (ictrl == 1) {
		    goto L400;
		}

/* SEPARATE NUCLEON+NUCLEON( BARYON RESONANCE+ BARYON RESONANCE ELASTIC */
/* COLLISION), BARYON RESONANCE+NUCLEON AND BARYON-PION */
/* COLLISIONS INTO THREE PARTS TO CHECK IF THEY ARE GOING TO SCATTER, */
/* WE only allow L/S to COLLIDE elastically with a nucleon and meson */
		if (abs(leadng_1.lb1) >= 14 && abs(leadng_1.lb1) <= 17 || abs(
			dpi_1.lb2) >= 14 && abs(dpi_1.lb2) <= 17) {
		    goto L400;
		}
/* IF PION+PION COLLISIONS GO TO 777 */
/* if pion+eta, eta+eta to create kaons go to 777 */
		if (leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5 && (dpi_1.lb2 >= 3 
			&& dpi_1.lb2 <= 5)) {
		    goto L777;
		}
		if (leadng_1.lb1 == 0 && (dpi_1.lb2 >= 3 && dpi_1.lb2 <= 5)) {
		    goto L777;
		}
		if (dpi_1.lb2 == 0 && (leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5)
			) {
		    goto L777;
		}
		if (leadng_1.lb1 == 0 && dpi_1.lb2 == 0) {
		    goto L777;
		}
/* we assume that rho and omega behave the same way as pions in */
/* kaon production */
/* (1) rho(omega)+rho(omega) */
		if (leadng_1.lb1 >= 25 && leadng_1.lb1 <= 28 && (dpi_1.lb2 >= 
			25 && dpi_1.lb2 <= 28)) {
		    goto L777;
		}
/* (2) rho(omega)+pion */
		if (leadng_1.lb1 >= 25 && leadng_1.lb1 <= 28 && (dpi_1.lb2 >= 
			3 && dpi_1.lb2 <= 5)) {
		    goto L777;
		}
		if (dpi_1.lb2 >= 25 && dpi_1.lb2 <= 28 && (leadng_1.lb1 >= 3 
			&& leadng_1.lb1 <= 5)) {
		    goto L777;
		}
/* (3) rho(omega)+eta */
		if (leadng_1.lb1 >= 25 && leadng_1.lb1 <= 28 && dpi_1.lb2 == 
			0) {
		    goto L777;
		}
		if (dpi_1.lb2 >= 25 && dpi_1.lb2 <= 28 && leadng_1.lb1 == 0) {
		    goto L777;
		}

/* if kaon+pion collisions go to 889 */
		if ((leadng_1.lb1 == 23 || leadng_1.lb1 == 21) && (dpi_1.lb2 
			>= 3 && dpi_1.lb2 <= 5)) {
		    goto L889;
		}
		if ((dpi_1.lb2 == 23 || dpi_1.lb2 == 21) && (leadng_1.lb1 >= 
			3 && leadng_1.lb1 <= 5)) {
		    goto L889;
		}

/* lin-2/06/03 skip all other (K K* Kbar K*bar) channels: */
/* SKIP all other K and K* RESCATTERINGS */
		if (abs(leadng_1.lb1) == 30 || abs(dpi_1.lb2) == 30) {
		    goto L400;
		}
		if (leadng_1.lb1 == 21 || dpi_1.lb2 == 21) {
		    goto L400;
		}
		if (leadng_1.lb1 == 23 || dpi_1.lb2 == 23) {
		    goto L400;
		}

/* IF PION+baryon COLLISION GO TO 3 */
		if (leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5 && (abs(dpi_1.lb2) 
			== 1 || abs(dpi_1.lb2) == 2 || abs(dpi_1.lb2) >= 6 && 
			abs(dpi_1.lb2) <= 13)) {
		    goto L3;
		}
		if (dpi_1.lb2 >= 3 && dpi_1.lb2 <= 5 && (abs(leadng_1.lb1) == 
			1 || abs(leadng_1.lb1) == 2 || abs(leadng_1.lb1) >= 6 
			&& abs(leadng_1.lb1) <= 13)) {
		    goto L3;
		}

/* IF rho(omega)+NUCLEON (baryon resonance) COLLISION GO TO 33 */
		if (leadng_1.lb1 >= 25 && leadng_1.lb1 <= 28 && (abs(
			dpi_1.lb2) == 1 || abs(dpi_1.lb2) == 2 || abs(
			dpi_1.lb2) >= 6 && abs(dpi_1.lb2) <= 13)) {
		    goto L33;
		}
		if (dpi_1.lb2 >= 25 && dpi_1.lb2 <= 28 && (abs(leadng_1.lb1) 
			== 1 || abs(leadng_1.lb1) == 2 || abs(leadng_1.lb1) >=
			 6 && abs(leadng_1.lb1) <= 13)) {
		    goto L33;
		}

/* IF ETA+NUCLEON (baryon resonance) COLLISIONS GO TO 547 */
		if (leadng_1.lb1 == 0 && (abs(dpi_1.lb2) == 1 || abs(
			dpi_1.lb2) == 2 || abs(dpi_1.lb2) >= 6 && abs(
			dpi_1.lb2) <= 13)) {
		    goto L547;
		}
		if (dpi_1.lb2 == 0 && (abs(leadng_1.lb1) == 1 || abs(
			leadng_1.lb1) == 2 || abs(leadng_1.lb1) >= 6 && abs(
			leadng_1.lb1) <= 13)) {
		    goto L547;
		}

/* IF NUCLEON+BARYON RESONANCE COLLISION GO TO 44 */
		if ((leadng_1.lb1 == 1 || leadng_1.lb1 == 2) && (dpi_1.lb2 > 
			5 && dpi_1.lb2 <= 13)) {
		    goto L44;
		}
		if ((dpi_1.lb2 == 1 || dpi_1.lb2 == 2) && (leadng_1.lb1 > 5 &&
			 leadng_1.lb1 <= 13)) {
		    goto L44;
		}
		if ((leadng_1.lb1 == -1 || leadng_1.lb1 == -2) && (dpi_1.lb2 <
			 -5 && dpi_1.lb2 >= -13)) {
		    goto L44;
		}
		if ((dpi_1.lb2 == -1 || dpi_1.lb2 == -2) && (leadng_1.lb1 < 
			-5 && leadng_1.lb1 >= -13)) {
		    goto L44;
		}

/* IF NUCLEON+NUCLEON COLLISION GO TO 4 */
		if ((leadng_1.lb1 == 1 || leadng_1.lb1 == 2) && (dpi_1.lb2 == 
			1 || dpi_1.lb2 == 2)) {
		    goto L4;
		}
		if ((leadng_1.lb1 == -1 || leadng_1.lb1 == -2) && (dpi_1.lb2 
			== -1 || dpi_1.lb2 == -2)) {
		    goto L4;
		}

/* IF BARYON RESONANCE+BARYON RESONANCE COLLISION GO TO 444 */
		if (leadng_1.lb1 > 5 && leadng_1.lb1 <= 13 && (dpi_1.lb2 > 5 
			&& dpi_1.lb2 <= 13)) {
		    goto L444;
		}
		if (leadng_1.lb1 < -5 && leadng_1.lb1 >= -13 && (dpi_1.lb2 < 
			-5 && dpi_1.lb2 >= -13)) {
		    goto L444;
		}

/* if L/S+L/S or L/s+nucleon go to 400 */
/* otherwise, develop a model for their collisions */
		if (leadng_1.lb1 < 3 && (dpi_1.lb2 >= 14 && dpi_1.lb2 <= 17)) 
			{
		    goto L400;
		}
		if (dpi_1.lb2 < 3 && (leadng_1.lb1 >= 14 && leadng_1.lb1 <= 
			17)) {
		    goto L400;
		}
		if (leadng_1.lb1 >= 14 && leadng_1.lb1 <= 17 && (dpi_1.lb2 >= 
			14 && dpi_1.lb2 <= 17)) {
		    goto L400;
		}

/* otherwise, go out of the loop */
		goto L400;


L547:
		if (leadng_1.lb1 * dpi_1.lb2 == 0) {
/* (1) FOR ETA+NUCLEON SYSTEM, we allow both elastic collision, */
/*     i.e. N*(1535) formation and kaon production */
/*     the total kaon production cross section is */
/*     ASSUMED to be THE SAME AS PION+NUCLEON COLLISIONS */
/* (2) for eta+baryon resonance we only allow kaon production */
/* Computing 2nd power */
		    r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
		    ece = r__1 * r__1;
		    xkaon0 = 0.f;
		    if (srt >= 1.63f && srt <= 1.7f) {
			xkaon0 = pnlka_(&srt);
		    }
		    if (srt > 1.7f) {
			xkaon0 = pnlka_(&srt) + pnska_(&srt);
		    }
/* bz3/7/99 neutralk */
		    xkaon0 *= 2.f;
/* bz3/7/99 neutralk end */
/* Here we negelect eta+n inelastic collisions other than the */
/* kaon production, therefore the total inelastic cross section */
/* xkaon equals to the xkaon0 (kaon production cross section) */
		    xkaon = xkaon0;
/* note here the xkaon is in unit of fm**2 */
		    xeta = xn1535_(&i1, &i2, &c__0);
		    if ((i__4 = ee_1.lb[i1 - 1], abs(i__4)) >= 6 && (i__5 = 
			    ee_1.lb[i1 - 1], abs(i__5)) <= 13 || (i__6 = 
			    ee_1.lb[i2 - 1], abs(i__6)) >= 6 && (i__7 = 
			    ee_1.lb[i2 - 1], abs(i__7)) <= 13) {
			xeta = 0.f;
		    }
		    if (xeta + xkaon <= 1e-6f) {
			goto L400;
		    }
		    dse = sqrt((xeta + xkaon) / 3.1415926f);
		    deltre = dse + .1f;
		    px1cm = pcx;
		    py1cm = pcy;
		    pz1cm = pcz;
/* CHECK IF N*(1535) resonance CAN BE FORMED */
		    distce_(&i1, &i2, &deltre, &dse, &input1_1.dt, &ece, &srt,
			     &ic, &pcx, &pcy, &pcz);
		    if (ic == -1) {
			goto L400;
		    }
		    kkk_1.ekaon[iss * 7 + 3] += 1;
		    if (xkaon0 / (xkaon + xeta) > ranart_(&rndf77_1.nseed)) {
/* kaon production, USE CREN TO CALCULATE THE MOMENTUM OF L/S K+ */
			cren_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock)
				;
/* kaon production */
			if (iblock == 7) {
			    ++(*lpn);
			} else if (iblock == -7) {
			}

			leadng_1.em1 = cc_1.e[i1 - 1];
			dpi_1.em2 = cc_1.e[i2 - 1];
			goto L440;
		    }
/* N*(1535) FORMATION */
		    resona = 1.f;
		    goto L98;
		}
/* IF PION+NUCLEON (baryon resonance) COLLISION THEN */
L3:
		px1cm = pcx;
		py1cm = pcy;
		pz1cm = pcz;
/* the total kaon production cross section for pion+baryon (resonance) is */
/* assumed to be the same as in pion+nucleon */
		xkaon0 = 0.f;
		if (srt >= 1.63f && srt <= 1.7f) {
		    xkaon0 = pnlka_(&srt);
		}
		if (srt > 1.7f) {
		    xkaon0 = pnlka_(&srt) + pnska_(&srt);
		}
		xkaon0 *= 2.f;

/* sp11/21/01  phi production: pi +N(D) -> phi + N(D) */
		xphi = 0.f;
		if ((leadng_1.lb1 >= 1 && leadng_1.lb1 <= 2 || leadng_1.lb1 >=
			 6 && leadng_1.lb1 <= 9 || (dpi_1.lb2 >= 1 && 
			dpi_1.lb2 <= 2 || dpi_1.lb2 >= 6 && dpi_1.lb2 <= 9)) 
			&& srt > 1.958f) {
		    pibphi_(&srt, &leadng_1.lb1, &dpi_1.lb2, &leadng_1.em1, &
			    dpi_1.em2, &xphi, &xphin);
		}
/* !! in fm^2 above */
/* if a pion collide with a baryon resonance, */
/* we only allow kaon production AND the reabsorption */
/* processes: Delta+pion-->N+pion, N*+pion-->N+pion */
/* Later put in pion+baryon resonance elastic */
/* cross through forming higher resonances implicitly. */
/*          If(em1.gt.1.or.em2.gt.1.)go to 31 */
		if ((i__4 = ee_1.lb[i1 - 1], abs(i__4)) >= 6 && (i__5 = 
			ee_1.lb[i1 - 1], abs(i__5)) <= 13 || (i__6 = ee_1.lb[
			i2 - 1], abs(i__6)) >= 6 && (i__7 = ee_1.lb[i2 - 1], 
			abs(i__7)) <= 13) {
		    goto L31;
		}
/* For pion+nucleon collisions: */
/* using the experimental pion+nucleon inelastic cross section, we assume it */
/* is exhausted by the Delta+pion, Delta+rho and Delta+omega production */
/* and kaon production. In the following we first check whether */
/* inelastic pion+n collision can happen or not, then determine in */
/* crpn whether it is through pion production or through kaon production */
/* note that the xkaon0 is the kaon production cross section */
/* Note in particular that: */
/* xkaon in the following is the total pion+nucleon inelastic cross section */
/* note here the xkaon is in unit of fm**2, xnpi is also in unit of fm**2 */
/* FOR PION+NUCLEON SYSTEM, THE MINIMUM S IS 1.2056 the minimum srt for */
/* elastic scattering, and it is 1.60 for pion production, 1.63 for LAMBDA+kaon */
/* production and 1.7 FOR SIGMA+KAON */
/* (EC = PION MASS+NUCLEON MASS+20MEV)**2 */
/* Computing 2nd power */
		r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
		ec = r__1 * r__1;
		xkaon = 0.f;
		if (srt > 1.23f) {
		    xkaon = (pionpp_(&srt) + pipp1_(&srt)) / 2.f;
		}
/* pion+nucleon elastic cross section is divided into two parts: */
/* (1) forming D(1232)+N*(1440) +N*(1535) */
/* (2) cross sections forming higher resonances are calculated as */
/*     the difference between the total elastic and (1), this part is */
/*     treated as direct process since we do not explicitLY include */
/*     higher resonances. */
/* the following is the resonance formation cross sections. */
/* 1. PION(+)+PROTON-->DELTA++,PION(-)+NEUTRON-->DELTA(-) */
		if (leadng_1.lb1 * dpi_1.lb2 == 5 || leadng_1.lb1 * dpi_1.lb2 
			== 6 && (leadng_1.lb1 == 3 || dpi_1.lb2 == 3) || (
			leadng_1.lb1 * dpi_1.lb2 == -3 || leadng_1.lb1 * 
			dpi_1.lb2 == -10 && (leadng_1.lb1 == 5 || dpi_1.lb2 ==
			 5))) {
		    xmax = 190.f;
		    xmaxn = 0.f;
		    xmaxn1 = 0.f;
		    xdirct = dirct1_(&srt);
		    goto L678;
		}
/* 2. PION(-)+PROTON-->DELTA0,PION(+)+NEUTRON-->DELTA+ */
/*   or N*(+)(1440) or N*(+)(1535) */
/* note the factor 2/3 is from the isospin consideration and */
/* the factor 0.6 or 0.5 is the branching ratio for the resonance to decay */
/* into pion+nucleon */
		if (leadng_1.lb1 * dpi_1.lb2 == 3 || leadng_1.lb1 * dpi_1.lb2 
			== 10 && (leadng_1.lb1 == 5 || dpi_1.lb2 == 5) || (
			leadng_1.lb1 * dpi_1.lb2 == -5 || leadng_1.lb1 * 
			dpi_1.lb2 == -6 && (leadng_1.lb1 == 3 || dpi_1.lb2 == 
			3))) {
		    xmax = 27.f;
		    xmaxn = 9.9999999999999982f;
		    xmaxn1 = 13.333333333333332f;
		    xdirct = dirct2_(&srt);
		    goto L678;
		}
/* 3. PION0+PROTON-->DELTA+,PION0+NEUTRON-->DELTA0, or N*(0)(1440) or N*(0)(1535) */
		if ((leadng_1.lb1 == 4 || dpi_1.lb2 == 4) && ((i__4 = 
			leadng_1.lb1 * dpi_1.lb2, abs(i__4)) == 4 || (i__5 = 
			leadng_1.lb1 * dpi_1.lb2, abs(i__5)) == 8)) {
		    xmax = 50.f;
		    xmaxn = 4.9999999999999991f;
		    xmaxn1 = 6.6666666666666661f;
		    xdirct = dirct3_(&srt);
		    goto L678;
		}
L678:
		xnpin1 = 0.f;
		xnpin = 0.f;
		xnpid = xnpi_(&i1, &i2, &c__1, &xmax);
		if (xmaxn1 != 0.f) {
		    xnpin1 = xnpi_(&i1, &i2, &c__2, &xmaxn1);
		}
		if (xmaxn != 0.f) {
		    xnpin = xnpi_(&i1, &i2, &c__0, &xmaxn);
		}
/* the following */
		xres = xnpid + xnpin + xnpin1;
		xnelas = xres + xdirct;
		icheck = 1;
		goto L34;
/* For pion + baryon resonance the reabsorption */
/* cross section is calculated from the detailed balance */
/* using reab(i1,i2,srt,ictrl), ictrl=1, 2 and 3 */
/* for pion, rho and omega + baryon resonance */
L31:
/* Computing 2nd power */
		r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
		ec = r__1 * r__1;
		xreab = reab_(&i1, &i2, &srt, &c__1);
/* lin-12/02/00 to satisfy detailed balance, forbid N* absorptions: */
		if (abs(leadng_1.lb1) >= 10 && abs(leadng_1.lb1) <= 13 || abs(
			dpi_1.lb2) >= 10 && abs(dpi_1.lb2) <= 13) {
		    xreab = 0.f;
		}
		xkaon = xkaon0 + xreab;
/* a constant of 10 mb IS USED FOR PION + N* RESONANCE, */
		if (abs(leadng_1.lb1) > 9 && abs(leadng_1.lb1) <= 13 || abs(
			dpi_1.lb2) > 9 && abs(dpi_1.lb2) <= 13) {
		    xnelas = 1.f;
		} else {
		    xnelas = dpion_(&leadng_1.em1, &dpi_1.em2, &leadng_1.lb1, 
			    &dpi_1.lb2, &srt);
		}
		icheck = 2;
L34:
		if (xnelas + xkaon + xphi <= 1e-6f) {
		    goto L400;
		}
		ds = sqrt((xnelas + xkaon + xphi) / 3.1415926f);
/* sp09/20/01 */
/*           totcr = xnelas+xkaon */
/*           if(srt .gt. 3.5)totcr = max1(totcr,3.) */
/*           DS=SQRT(totcr/PI) */
/* sp09/20/01 end */
		deltar = ds + .1f;
		distce_(&i1, &i2, &deltar, &ds, &input1_1.dt, &ec, &srt, &ic, 
			&pcx, &pcy, &pcz);
		if (ic == -1) {
		    goto L400;
		}
		kkk_1.ekaon[iss * 7 + 3] += 1;
/* *** */
/* check what kind of collision has happened */
/* (1) pion+baryon resonance */
/* if direct elastic process */
		if (icheck == 2) {
/*  !!sp11/21/01 */
		    if (xnelas / (xnelas + xkaon + xphi) >= ranart_(&
			    rndf77_1.nseed)) {
/*               call Crdir(PX1CM,PY1CM,PZ1CM,SRT,I1,I2) */
			crdir_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &
				iblock);
			goto L440;
		    } else {
/* for inelastic process, go to 96 to check */
/* kaon production and pion reabsorption : pion+D(N*)-->pion+N */
			goto L96;
		    }
		}
/* (2) pion+n */
/* CHECK IF inELASTIC COLLISION IS POSSIBLE FOR PION+N COLLISIONS */
/* lin-8/17/00 typo corrected, many other occurences: */
/*        IF(XKAON/(XKAON+Xnelas).GT.RANART(NSEED))GO TO 95 */
		if ((xkaon + xphi) / (xkaon + xphi + xnelas) > ranart_(&
			rndf77_1.nseed)) {
		    goto L95;
		}
/* direct process */
		if (xdirct / xnelas >= ranart_(&rndf77_1.nseed)) {
/*               call Crdir(PX1CM,PY1CM,PZ1CM,SRT,I1,I2) */
		    crdir_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock);
		    goto L440;
		}
/* now resonance formation or direct process (higher resonances) */
		if (leadng_1.lb1 * dpi_1.lb2 == 5 || leadng_1.lb1 * dpi_1.lb2 
			== 6 && (leadng_1.lb1 == 3 || dpi_1.lb2 == 3) || (
			leadng_1.lb1 * dpi_1.lb2 == -3 || leadng_1.lb1 * 
			dpi_1.lb2 == -10 && (leadng_1.lb1 == 5 || dpi_1.lb2 ==
			 5))) {

/* ONLY DELTA RESONANCE IS POSSIBLE, go to 99 */
		    goto L99;
		} else {
/* NOW BOTH DELTA AND N* RESORANCE ARE POSSIBLE */
/* DETERMINE THE RESORANT STATE BY USING THE MONTRE CARLO METHOD */
		    xx = (xnpin + xnpin1) / xres;
		    if (ranart_(&rndf77_1.nseed) < xx) {
/* N* RESONANCE IS SELECTED */
/* decide N*(1440) or N*(1535) formation */
			xx0 = xnpin / (xnpin + xnpin1);
			if (ranart_(&rndf77_1.nseed) < xx0) {
			    resona = 0.f;
/* N*(1440) formation */
			    goto L97;
			} else {
/* N*(1535) formation */
			    resona = 1.f;
			    goto L98;
			}
		    } else {
/* DELTA RESONANCE IS SELECTED */
			goto L99;
		    }
		}
L97:
		if (resona == 0.f) {
/* N*(1440) IS PRODUCED,WE DETERMINE THE CHARGE STATE OF THE PRODUCED N* */
		    i__ = i1;
		    if (leadng_1.em1 < .6f) {
			i__ = i2;
		    }
/* (0.1) n+pion(+)-->N*(+) */
		    if (leadng_1.lb1 * dpi_1.lb2 == 10 && (leadng_1.lb1 == 5 
			    || dpi_1.lb2 == 5) || leadng_1.lb1 * dpi_1.lb2 == 
			    -6 && (leadng_1.lb1 == 3 || dpi_1.lb2 == 3)) {
			ee_1.lb[i__ - 1] = 11;
			goto L303;
		    }
/* (0.2) p+pion(0)-->N*(+) */
/*            IF(LB(I1)*LB(I2).EQ.4.AND.(LB(I1).EQ.1.OR.LB(I2).EQ.1))THEN */
		    if ((i__4 = ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1], abs(i__4)) 
			    == 4 && (ee_1.lb[i1 - 1] == 4 || ee_1.lb[i2 - 1] 
			    == 4)) {
			ee_1.lb[i__ - 1] = 11;
			goto L303;
		    }
/* (0.3) n+pion(0)-->N*(0) */
/*            IF(LB(I1)*LB(I2).EQ.8.AND.(LB(I1).EQ.2.OR.LB(I2).EQ.2))THEN */
		    if ((i__4 = ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1], abs(i__4)) 
			    == 8 && (ee_1.lb[i1 - 1] == 4 || ee_1.lb[i2 - 1] 
			    == 4)) {
			ee_1.lb[i__ - 1] = 10;
			goto L303;
		    }
/* (0.4) p+pion(-)-->N*(0) */
/*            IF(LB(I1)*LB(I2).EQ.3)THEN */
		    if (ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1] == 3 || ee_1.lb[i1 
			    - 1] * ee_1.lb[i2 - 1] == -5) {
			ee_1.lb[i__ - 1] = 10;
		    }
L303:
		    dreson_(&i1, &i2);
		    if (leadng_1.lb1 < 0 || dpi_1.lb2 < 0) {
			ee_1.lb[i__ - 1] = -ee_1.lb[i__ - 1];
		    }
		    ++(*lres);
		    goto L101;
/* COM: GO TO 101 TO CHANGE THE PHASE SPACE DENSITY OF THE NUCLEON */
		}
L98:
		if (resona == 1.f) {
/* N*(1535) IS PRODUCED, WE DETERMINE THE CHARGE STATE OF THE PRODUCED N* */
		    i__ = i1;
		    if (leadng_1.em1 < .6f) {
			i__ = i2;
		    }
/* note: this condition applies to both eta and pion */
/* (0.1) n+pion(+)-->N*(+) */
/*            IF(LB1*LB2.EQ.10.AND.(LB1.EQ.2.OR.LB2.EQ.2))THEN */
		    if (leadng_1.lb1 * dpi_1.lb2 == 10 && (leadng_1.lb1 == 5 
			    || dpi_1.lb2 == 5) || leadng_1.lb1 * dpi_1.lb2 == 
			    -6 && (leadng_1.lb1 == 3 || dpi_1.lb2 == 3)) {
			ee_1.lb[i__ - 1] = 13;
			goto L304;
		    }
/* (0.2) p+pion(0)-->N*(+) */
/*            IF(LB(I1)*LB(I2).EQ.4.AND.(LB(I1).EQ.1.OR.LB(I2).EQ.1))THEN */
		    if ((i__4 = ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1], abs(i__4)) 
			    == 4 && (ee_1.lb[i1 - 1] == 4 || ee_1.lb[i2 - 1] 
			    == 4)) {
			ee_1.lb[i__ - 1] = 13;
			goto L304;
		    }
/* (0.3) n+pion(0)-->N*(0) */
/*            IF(LB(I1)*LB(I2).EQ.8.AND.(LB(I1).EQ.2.OR.LB(I2).EQ.2))THEN */
		    if ((i__4 = ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1], abs(i__4)) 
			    == 8 && (ee_1.lb[i1 - 1] == 4 || ee_1.lb[i2 - 1] 
			    == 4)) {
			ee_1.lb[i__ - 1] = 12;
			goto L304;
		    }
/* (0.4) p+pion(-)-->N*(0) */
/*            IF(LB(I1)*LB(I2).EQ.3)THEN */
		    if (ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1] == 3 || ee_1.lb[i1 
			    - 1] * ee_1.lb[i2 - 1] == -5) {
			ee_1.lb[i__ - 1] = 12;
			goto L304;
		    }
/* (0.5) p+eta-->N*(+)(1535),n+eta-->N*(0)(1535) */
		    if (ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1] == 0) {
/*            if((lb(i1).eq.1).or.(lb(i2).eq.1))then */
			if ((i__4 = ee_1.lb[i1 - 1], abs(i__4)) == 1 || (i__5 
				= ee_1.lb[i2 - 1], abs(i__5)) == 1) {
			    ee_1.lb[i__ - 1] = 13;
			    goto L304;
			} else {
			    ee_1.lb[i__ - 1] = 12;
			}
		    }
L304:
		    dreson_(&i1, &i2);
		    if (leadng_1.lb1 < 0 || dpi_1.lb2 < 0) {
			ee_1.lb[i__ - 1] = -ee_1.lb[i__ - 1];
		    }
		    ++(*lres);
		    goto L101;
/* COM: GO TO 101 TO CHANGE THE PHASE SPACE DENSITY OF THE NUCLEON */
		}
/* DELTA IS PRODUCED,IN THE FOLLOWING WE DETERMINE THE */
/* CHARGE STATE OF THE PRODUCED DELTA */
L99:
		++(*lres);
		i__ = i1;
		if (leadng_1.em1 <= .6f) {
		    i__ = i2;
		}
/* (1) p+pion(+)-->DELTA(++) */
/*        IF(LB(I1)*LB(I2).EQ.5)THEN */
		if (ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1] == 5 || ee_1.lb[i1 - 1] 
			* ee_1.lb[i2 - 1] == -3) {
		    ee_1.lb[i__ - 1] = 9;
		    goto L305;
		}
/* (2) p+pion(0)-->delta(+) */
/*        IF(LB(I1)*LB(I2).EQ.4.AND.(LB(I1).EQ.1.OR.LB(I2).EQ.1))then */
		if ((i__4 = ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1], abs(i__4)) == 
			4 && (ee_1.lb[i1 - 1] == 4 || ee_1.lb[i2 - 1] == 4)) {
		    ee_1.lb[i__ - 1] = 8;
		    goto L305;
		}
/* (3) n+pion(+)-->delta(+) */
/*        IF(LB(I1)*LB(I2).EQ.10.AND.(LB(I1).EQ.2.OR.LB(I2).EQ.2))THEN */
		if (ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1] == 10 && (ee_1.lb[i1 - 
			1] == 5 || ee_1.lb[i2 - 1] == 5) || ee_1.lb[i1 - 1] * 
			ee_1.lb[i2 - 1] == -6 && (ee_1.lb[i1 - 1] == 3 || 
			ee_1.lb[i2 - 1] == 3)) {
		    ee_1.lb[i__ - 1] = 8;
		    goto L305;
		}
/* (4) n+pion(0)-->delta(0) */
/*        IF(LB(I1)*LB(I2).EQ.8.AND.(LB(I1).EQ.2.OR.LB(I2).EQ.2))THEN */
		if ((i__4 = ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1], abs(i__4)) == 
			8 && (ee_1.lb[i1 - 1] == 4 || ee_1.lb[i2 - 1] == 4)) {
		    ee_1.lb[i__ - 1] = 7;
		    goto L305;
		}
/* (5) p+pion(-)-->delta(0) */
/*        IF(LB(I1)*LB(I2).EQ.3)THEN */
		if (ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1] == 3 || ee_1.lb[i1 - 1] 
			* ee_1.lb[i2 - 1] == -5) {
		    ee_1.lb[i__ - 1] = 7;
		    goto L305;
		}
/* (6) n+pion(-)-->delta(-) */
/*        IF(LB(I1)*LB(I2).EQ.6.AND.(LB(I1).EQ.2.OR.LB(I2).EQ.2))THEN */
		if (ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1] == 6 && (ee_1.lb[i1 - 1]
			 == 3 || ee_1.lb[i2 - 1] == 3) || ee_1.lb[i1 - 1] * 
			ee_1.lb[i2 - 1] == -10 && (ee_1.lb[i1 - 1] == 5 || 
			ee_1.lb[i2 - 1] == 5)) {
		    ee_1.lb[i__ - 1] = 6;
		}
L305:
		dreson_(&i1, &i2);
		if (leadng_1.lb1 < 0 || dpi_1.lb2 < 0) {
		    ee_1.lb[i__ - 1] = -ee_1.lb[i__ - 1];
		}
		goto L101;
/* sp-11/08/01 K* */
/* FOR kaON+pion COLLISIONS, form K* (bar) or */
/* La/Si-bar + N <-- pi + K+ */
/* La/Si + N-bar <-- pi + K- */
/* phi + K <-- pi + K */
/* lin (rho,omega) + K* <-- pi + K */
L889:
		px1cm = pcx;
		py1cm = pcy;
		pz1cm = pcz;
/* Computing 2nd power */
		r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
		ec = r__1 * r__1;
/* the cross section is from C.M. Ko, PRC 23, 2760 (1981). */
/* Computing 2nd power */
		r__1 = srt - .895f;
		spika = 60.f / (r__1 * r__1 * 4.f / .0025000000000000005f + 
			1.f);

/* c       if(lb(i1).eq.23.or.lb(i2).eq.23)then   !! block  K- + pi->La + B-bar */
		crkpla_(&px1cm, &py1cm, &pz1cm, &ec, &srt, &spika, &emm1, &
			emm2, &lbp1, &lbp2, &i1, &i2, &icase, &c_b106);
/* c */
/* * only K* or K*bar formation */
/*       else */
/*      DSkn=SQRT(spika/PI/10.) */
/*      dsknr=dskn+0.1 */
/*      CALL DISTCE(I1,I2,dsknr,DSkn,DT,EC,SRT,IC, */
/*    1     PX1CM,PY1CM,PZ1CM) */
/*        IF(IC.EQ.-1) GO TO 400 */
/*       icase = 1 */
/*      endif */

		if (icase == 0) {
		    iblock = 0;
		    goto L400;
		}
		if (icase == 1) {
		    ksreso_(&i1, &i2);
/* lin-4/30/03 give non-zero iblock for resonance selections: */
		    iblock = 171;
/* test off for resonance (phi, K*) studies: */
/*             if(iabs(lb(i1)).eq.30) then */
/*             write(17,112) 'ks',lb(i1),p(1,i1),p(2,i1),p(3,i1),e(i1),nt */
/*             elseif(iabs(lb(i2)).eq.30) then */
/*             write(17,112) 'ks',lb(i2),p(1,i2),p(2,i2),p(3,i2),e(i2),nt */
/*             endif */

		    ++(*lres);
		    goto L101;
		} else if (icase == 2) {
		    iblock = 71;

/* La/Si (bar) formation */
		} else if (abs(icase) == 5) {
		    iblock = 88;
		} else {

/* phi formation */
		    iblock = 222;
		}
		ee_1.lb[i1 - 1] = lbp1;
		ee_1.lb[i2 - 1] = lbp2;
		cc_1.e[i1 - 1] = emm1;
		cc_1.e[i2 - 1] = emm2;
		leadng_1.em1 = cc_1.e[i1 - 1];
		dpi_1.em2 = cc_1.e[i2 - 1];
		ntag = 0;
		goto L440;

L33:
		leadng_1.em1 = cc_1.e[i1 - 1];
		dpi_1.em2 = cc_1.e[i2 - 1];
/* (1) if rho or omega collide with a nucleon we allow both elastic */
/*     scattering and kaon production to happen if collision conditions */
/*     are satisfied. */
/* (2) if rho or omega collide with a baryon resonance we allow */
/*     kaon production, pion reabsorption: rho(omega)+D(N*)-->pion+N */
/*     and NO elastic scattering to happen */
		xelstc = 0.f;
		if (leadng_1.lb1 >= 25 && leadng_1.lb1 <= 28 && (abs(
			dpi_1.lb2) == 1 || abs(dpi_1.lb2) == 2)) {
		    xelstc = erhon_(&leadng_1.em1, &dpi_1.em2, &leadng_1.lb1, 
			    &dpi_1.lb2, &srt);
		}
		if (dpi_1.lb2 >= 25 && dpi_1.lb2 <= 28 && (abs(leadng_1.lb1) 
			== 1 || abs(leadng_1.lb1) == 2)) {
		    xelstc = erhon_(&leadng_1.em1, &dpi_1.em2, &leadng_1.lb1, 
			    &dpi_1.lb2, &srt);
		}
/* Computing 2nd power */
		r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
		ec = r__1 * r__1;
/* the kaon production cross section is */
		xkaon0 = 0.f;
		if (srt >= 1.63f && srt <= 1.7f) {
		    xkaon0 = pnlka_(&srt);
		}
		if (srt > 1.7f) {
		    xkaon0 = pnlka_(&srt) + pnska_(&srt);
		}
		if (xkaon0 < 0.f) {
		    xkaon0 = 0.f;
		}
/* bz3/7/99 neutralk */
		xkaon0 *= 2.f;
/* bz3/7/99 neutralk end */
/* the total inelastic cross section for rho(omega)+N is */
		xkaon = xkaon0;
		ichann = 0;
/* the total inelastic cross section for rho (omega)+D(N*) is */
/* xkaon=xkaon0+reab(**) */
/* sp11/21/01  phi production: rho + N(D) -> phi + N(D) */
		xphi = 0.f;
		if (((leadng_1.lb1 >= 1 && leadng_1.lb1 <= 2 || leadng_1.lb1 
			>= 6 && leadng_1.lb1 <= 9) && (dpi_1.lb2 >= 25 && 
			dpi_1.lb2 <= 27) || (dpi_1.lb2 >= 1 && dpi_1.lb2 <= 2 
			|| dpi_1.lb2 >= 6 && dpi_1.lb2 <= 9) && (leadng_1.lb1 
			>= 25 && leadng_1.lb1 <= 27)) && srt > 1.958f) {
		    pibphi_(&srt, &leadng_1.lb1, &dpi_1.lb2, &leadng_1.em1, &
			    dpi_1.em2, &xphi, &xphin);
		}
/* !! in fm^2 above */

		if (abs(leadng_1.lb1) >= 6 && dpi_1.lb2 >= 25 || leadng_1.lb1 
			>= 25 && abs(dpi_1.lb2) >= 6) {
		    ichann = 1;
		    ictrl = 2;
		    if (leadng_1.lb1 == 28 || dpi_1.lb2 == 28) {
			ictrl = 3;
		    }
		    xreab = reab_(&i1, &i2, &srt, &ictrl);
/* lin-12/02/00 to satisfy detailed balance, forbid N* absorptions: */
		    if (abs(leadng_1.lb1) >= 10 && abs(leadng_1.lb1) <= 13 || 
			    abs(dpi_1.lb2) >= 10 && abs(dpi_1.lb2) <= 13) {
			xreab = 0.f;
		    }
		    if (xreab < 0.f) {
			xreab = 1e-6f;
		    }
		    xkaon = xkaon0 + xreab;
		    xelstc = 1.f;
		}
		ds = sqrt((xkaon + xphi + xelstc) / 3.1415926f);

/* sp09/20/01 */
/*           totcr = xelstc+xkaon */
/*           if(srt .gt. 3.5)totcr = max1(totcr,3.) */
/*           DS=SQRT(totcr/PI) */
/* sp09/20/01 end */

		deltar = ds + .1f;
		px1cm = pcx;
		py1cm = pcy;
		pz1cm = pcz;
/* CHECK IF the collision can happen */
		distce_(&i1, &i2, &deltar, &ds, &input1_1.dt, &ec, &srt, &ic, 
			&pcx, &pcy, &pcz);
		if (ic == -1) {
		    goto L400;
		}
		kkk_1.ekaon[iss * 7 + 3] += 1;
/* * */
/* NOW rho(omega)+N or D(N*) COLLISION IS POSSIBLE */
/* (1) check elastic collision */
		if (xelstc / (xelstc + xkaon + xphi) > ranart_(&
			rndf77_1.nseed)) {
/*       call crdir(px1CM,py1CM,pz1CM,srt,I1,i2) */
		    crdir_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock);
		    goto L440;
		}
/* (2) check pion absorption or kaon production */
		crrd_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock, &
			xkaon0, &xkaon, &xphi, &xphin);
/* kaon production */
/* sp05/16/01 */
		if (iblock == 7) {
		    ++(*lpn);
		} else if (iblock == -7) {
		}
/* sp05/16/01 end */
/* rho obsorption */
		if (iblock == 81) {
		    ++lrhor;
		}
/* omega obsorption */
		if (iblock == 82) {
		    ++lomgar;
		}
		leadng_1.em1 = cc_1.e[i1 - 1];
		dpi_1.em2 = cc_1.e[i2 - 1];
		goto L440;
/* for pion+n now using the subroutine crpn to change */
/* the particle label and set the new momentum of L/S+K final state */
L95:
/* NOW PION+N INELASTIC COLLISION IS POSSIBLE */
/* check pion production or kaon production */
		crpn_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock, &
			xkaon0, &xkaon, &xphi, &xphin);
/* kaon production */
/* sp05/16/01 */
		if (iblock == 7) {
		    ++(*lpn);
		} else if (iblock == -7) {
		}
/* sp05/16/01 end */
/* pion production */
		if (iblock == 77) {
		    ++(*lpd);
		}
/* rho production */
		if (iblock == 78) {
		    ++(*lrho);
		}
/* omega production */
		if (iblock == 79) {
		    ++(*lomega);
		}
		leadng_1.em1 = cc_1.e[i1 - 1];
		dpi_1.em2 = cc_1.e[i2 - 1];
		goto L440;
/* for pion+D(N*) now using the subroutine crpd to */
/* (1) check kaon production or pion reabsorption */
/* (2) change the particle label and set the new */
/*     momentum of L/S+K final state */
L96:
		crpd_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock, &
			xkaon0, &xkaon, &xphi, &xphin);
/* kaon production */
/* sp05/16/01 */
		if (iblock == 7) {
		    ++(*lpn);
		} else if (iblock == -7) {
		}
/* sp05/16/01 end */
/* pion obserption */
		if (iblock == 80) {
		    ++lpdr;
		}
		leadng_1.em1 = cc_1.e[i1 - 1];
		dpi_1.em2 = cc_1.e[i2 - 1];
		goto L440;
/* CALCULATE KAON PRODUCTION PROBABILITY FROM PION + N COLLISIONS */
/*        IF(SRT.GT.1.615)THEN */
/*        CALL PKAON(SRT,XXp,PK) */
/*        TKAON(7)=TKAON(7)+PK */
/*        EKAON(7,ISS)=EKAON(7,ISS)+1 */
/*        CALL KSPEC1(SRT,PK) */
/*        call LK(3,srt,iseed,pk) */
/*        ENDIF */
/* negelecting the pauli blocking at high energies */
L101:
		if (cc_1.e[i2 - 1] == 0.f) {
		    goto L600;
		}
		if (cc_1.e[i1 - 1] == 0.f) {
		    goto L800;
		}
/* IF NUCLEON+BARYON RESONANCE COLLISIONS */
L44:
/* CALCULATE THE TOTAL CROSS SECTION OF NUCLEON+ BARYON RESONANCE COLLISION */
/* WE ASSUME THAT THE ELASTIC CROSS SECTION IS THE SAME AS NUCLEON+NUCLEON */
/* COM: WE USE THE PARAMETERISATION BY CUGNON FOR LOW ENERGIES */
/*      AND THE PARAMETERIZATIONS FROM CERN DATA BOOK FOR HIGHER */
/*      ENERGIES. THE CUTOFF FOR THE TOTAL CROSS SECTION IS 55 MB */
		cutoff = leadng_1.em1 + dpi_1.em2 + .02f;
		if (srt <= cutoff) {
		    goto L400;
		}
		if (srt > 2.245f) {
		    signn = pp2_(&srt);
		} else {
		    signn = 35.f / ((srt - cutoff) * 100.f + 1.f) + 20.f;
		}
		xnd_(&pcx, &pcy, &pcz, &srt, &i1, &i2, &xinel, &sigk, &xsk1, &
			xsk2, &xsk3, &xsk4, &xsk5);
		sig = signn + xinel;
/* For nucleon+baryon resonance collision, the minimum cms**2 energy is */
/* Computing 2nd power */
		r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
		ec = r__1 * r__1;
/* CHECK THE DISTENCE BETWEEN THE TWO PARTICLES */
		px1cm = pcx;
		py1cm = pcy;
		pz1cm = pcz;
/* lin-6/2008 Deuteron production: */
		ianti = 0;
		if (ee_1.lb[i1 - 1] < 0 && ee_1.lb[i2 - 1] < 0) {
		    ianti = 1;
		}
		sbbdm_(&srt, &sdprod, &ianti, &lbm, &xmm, &pfinal);
		sig += sdprod;
/* lin-6/2008 perturbative treatment of deuterons: */
		ipdflag = 0;
		if (para8_1.idpert == 1) {
		    ipert1 = 1;
		    sigr0 = sig;
		    dspert = sqrt(sigr0 / 3.1415926f / 10.f);
		    dsrpert = dspert + .1f;
		    distce_(&i1, &i2, &dsrpert, &dspert, &input1_1.dt, &ec, &
			    srt, &ic, &px1cm, &py1cm, &pz1cm);
		    if (ic == -1) {
			goto L363;
		    }
		    signn0 = 0.f;
		    crnd_(&irun, &px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &
			    iblock, &signn0, &sigr0, &sigk, &xsk1, &xsk2, &
			    xsk3, &xsk4, &xsk5, nt, &ipert1);
/*     &  IBLOCK,SIGNN,SIG,sigk,xsk1,xsk2,xsk3,xsk4,xsk5) */
		    ipdflag = 1;
L363:
		    ipert1 = 0;
		}
		if (para8_1.idpert == 2) {
		    ipert1 = 1;
		}

		ds = sqrt(sig / 31.415925999999999f);
		deltar = ds + .1f;
		distce_(&i1, &i2, &deltar, &ds, &input1_1.dt, &ec, &srt, &ic, 
			&px1cm, &py1cm, &pz1cm);
/*        IF(IC.EQ.-1)GO TO 400 */
		if (ic == -1) {
		    if (ipdflag == 1) {
			iblock = 501;
		    }
		    goto L400;
		}
		kkk_1.ekaon[iss * 7 + 2] += 1;
/* CALCULATE KAON PRODUCTION PROBABILITY FROM NUCLEON + BARYON RESONANCE */
/* COLLISIONS */
		goto L361;
/* CHECK WHAT KIND OF COLLISION HAS HAPPENED */
L361:
		crnd_(&irun, &px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock, 
			&signn, &sig, &sigk, &xsk1, &xsk2, &xsk3, &xsk4, &
			xsk5, nt, &ipert1);
/*     &  IBLOCK,SIGNN,SIG,sigk,xsk1,xsk2,xsk3,xsk4,xsk5) */
		if (iblock == 0 && ipdflag == 1) {
		    iblock = 501;
		}
		if (iblock == 11) {
		    ++(*lndk);
		    goto L400;
/*        elseIF(IBLOCK.EQ.-11) then */
		} else if (iblock == -11 || iblock == 501) {
		    goto L400;
		}
		if (iblock == 222) {
/*    !! sp12/17/01 */
		    goto L400;
		}
		leadng_1.em1 = cc_1.e[i1 - 1];
		dpi_1.em2 = cc_1.e[i2 - 1];
		goto L440;
/* IF NUCLEON+NUCLEON OR BARYON RESONANCE+BARYON RESONANCE COLLISIONS */
L4:
/* PREPARE THE EALSTIC CROSS SECTION FOR BARYON+BARYON COLLISIONS */
/* COM: WE USE THE PARAMETERISATION BY CUGNON FOR SRT LEQ 2.0 GEV */
/*      AND THE PARAMETERIZATIONS FROM CERN DATA BOOK FOR HIGHER */
/*      ENERGIES. THE CUTOFF FOR THE TOTAL CROSS SECTION IS 55 MB */
/*      WITH LOW-ENERGY-CUTOFF */
		cutoff = leadng_1.em1 + dpi_1.em2 + .14f;
/* AT HIGH ENERGIES THE ISOSPIN DEPENDENCE IS NEGLIGIBLE */
/* THE TOTAL CROSS SECTION IS TAKEN AS THAT OF THE PP */
/* ABOVE E_KIN=800 MEV, WE USE THE ISOSPIN INDEPENDNET XSECTION */
		if (srt > 2.245f) {
		    sig = ppt_(&srt);
		    signn = sig - pp1_(&srt);
		} else {
/* AT LOW ENERGIES THE ISOSPIN DEPENDENCE FOR NN COLLISION IS STRONG */
		    sig = xpp_(&srt);
		    if (zet[ee_1.lb[i1 - 1] + 45] * zet[ee_1.lb[i2 - 1] + 45] 
			    <= 0.f) {
			sig = xnp_(&srt);
		    }
		    if (zet[ee_1.lb[i1 - 1] + 45] * zet[ee_1.lb[i2 - 1] + 45] 
			    > 0.f) {
			sig = xpp_(&srt);
		    }
		    if (zet[ee_1.lb[i1 - 1] + 45] == 0.f && zet[ee_1.lb[i2 - 
			    1] + 45] == 0.f) {
			sig = xpp_(&srt);
		    }
		    if (ee_1.lb[i1 - 1] == -1 && ee_1.lb[i2 - 1] == -2 || 
			    ee_1.lb[i2 - 1] == -1 && ee_1.lb[i1 - 1] == -2) {
			sig = xnp_(&srt);
		    }
/*     WITH LOW-ENERGY-CUTOFF */
		    if (srt < 1.897f) {
			signn = sig;
		    } else {
			signn = 35.f / ((srt - 1.897f) * 100.f + 1.f) + 20.f;
		    }
		}
		px1cm = pcx;
		py1cm = pcy;
		pz1cm = pcz;
/* lin-5/2008 Deuteron production cross sections were not included */
/*     in the previous parameterized inelastic cross section of NN collisions */
/*     (SIGinel=SIG-SIGNN), so they are added here: */
		ianti = 0;
		if (ee_1.lb[i1 - 1] < 0 && ee_1.lb[i2 - 1] < 0) {
		    ianti = 1;
		}
		sbbdm_(&srt, &sdprod, &ianti, &lbm, &xmm, &pfinal);
		sig += sdprod;

/* lin-5/2008 perturbative treatment of deuterons: */
		ipdflag = 0;
		if (para8_1.idpert == 1) {
/*     For idpert=1: ipert1=1 means we will first treat deuteron perturbatively, */
/*     then we set ipert1=0 to treat regular NN or NbarNbar collisions including */
/*     the regular deuteron productions. */
/*     ipdflag=1 means perturbative deuterons are produced here: */
		    ipert1 = 1;
		    ec = 4.0481439999999997f;
/*     Use the same cross section for NN/NNBAR collisions */
/*     to trigger perturbative production */
		    sigr0 = sig;
/*     One can also trigger with X*sbbdm() so the weight will not be too small; */
/*     but make sure to limit the maximum trigger Xsec: */
/*           sigr0=sdprod*25. */
/*           if(sigr0.ge.100.) sigr0=100. */
		    dspert = sqrt(sigr0 / 3.1415926f / 10.f);
		    dsrpert = dspert + .1f;
		    distce_(&i1, &i2, &dsrpert, &dspert, &input1_1.dt, &ec, &
			    srt, &ic, &px1cm, &py1cm, &pz1cm);
		    if (ic == -1) {
			goto L365;
		    }
		    signn0 = 0.f;
		    crnn_(&irun, &px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &
			    iblock, &ntag, &signn0, &sigr0, nt, &ipert1);
		    ipdflag = 1;
L365:
		    ipert1 = 0;
		}
		if (para8_1.idpert == 2) {
		    ipert1 = 1;
		}

/* lin-5/2008 in case perturbative deuterons are produced for idpert=1: */
/*        IF(SIGNN.LE.0)GO TO 400 */
		if (signn <= 0.f) {
		    if (ipdflag == 1) {
			iblock = 501;
		    }
		    goto L400;
		}

		ec = 3.59709f;
		ds = sqrt(sig / 3.1415926f / 10.f);
		dsr = ds + .1f;
		if (cc_1.e[i1 - 1] >= 1.f && cc_1.e[i2 - 1] >= 1.f) {
		    ec = 4.75f;
		}
		distce_(&i1, &i2, &dsr, &ds, &input1_1.dt, &ec, &srt, &ic, &
			px1cm, &py1cm, &pz1cm);
/* lin-5/2008 in case perturbative deuterons are produced above: */
/*        IF(IC.EQ.-1) GO TO 400 */
		if (ic == -1) {
		    if (ipdflag == 1) {
			iblock = 501;
		    }
		    goto L400;
		}

/* CALCULATE KAON PRODUCTION PROBABILITY FROM NUCLEON+NUCLEON OR */
/* RESONANCE+RESONANCE COLLISIONS */
		goto L362;
/* CHECK WHAT KIND OF COLLISION HAS HAPPENED */
L362:
		kkk_1.ekaon[iss * 7] += 1;
		crnn_(&irun, &px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock, 
			&ntag, &signn, &sig, nt, &ipert1);
/* lin-5/2008 give iblock # in case pert deuterons are produced for idpert=1: */
		if (iblock == 0 && ipdflag == 1) {
		    iblock = 501;
		}
/* lin-5/2008 add iblock # for deuteron formation: */
/*        IF(IBLOCK.EQ.4.OR.IBLOCK.Eq.9.or.iblock.ge.44.OR.IBLOCK.EQ.-9 */
/*     &       .or.iblock.eq.222)THEN */
		if (iblock == 4 || iblock == 9 || iblock >= 44 || iblock == 
			-9 || iblock == 222 || iblock == 501) {

/*     !! sp12/17/01 above */
/* momentum of the three particles in the final state have been calculated */
/* in the crnn, go out of the loop */
		    ++(*lcoll);
		    if (iblock == 4) {
			++(*ldirt);
		    } else if (iblock == 44) {
			++(*lddrho);
		    } else if (iblock == 45) {
			++(*lnnrho);
		    } else if (iblock == 46) {
			++(*lnnom);
		    } else if (iblock == 222) {
		    } else if (iblock == 9) {
			++(*lnnk);
		    } else if (iblock == -9) {
		    }
		    goto L400;
		}
		leadng_1.em1 = cc_1.e[i1 - 1];
		dpi_1.em2 = cc_1.e[i2 - 1];
		goto L440;
/* lin-8/2008 B+B->Deuteron+Meson over */

/* lin-8/2008 Deuteron+Meson->B+B collisions: */
L505:
		ianti = 0;
		if (ee_1.lb[i1 - 1] < 0 || ee_1.lb[i2 - 1] < 0) {
		    ianti = 1;
		}
		sdmbb_(&srt, &sdm, &ianti);
		px1cm = pcx;
		py1cm = pcy;
		pz1cm = pcz;
/*     minimum srt**2, note a 2.012GeV lower cutoff is used in N+N->Deuteron+pi: */
		ec = 4.0481439999999997f;
		ds = sqrt(sdm / 31.4f);
		dsr = ds + .1f;
		distce_(&i1, &i2, &dsr, &ds, &input1_1.dt, &ec, &srt, &ic, &
			px1cm, &py1cm, &pz1cm);
		if (ic == -1) {
		    goto L400;
		}
		crdmbb_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock, &
			ntag, &sdm, nt, &ianti);
		++(*lcoll);
		goto L400;
/* lin-8/2008 Deuteron+Meson->B+B collisions over */

/* lin-9/2008 Deuteron+Baryon elastic collisions: */
L506:
		ianti = 0;
		if (ee_1.lb[i1 - 1] < 0 || ee_1.lb[i2 - 1] < 0) {
		    ianti = 1;
		}
		sdbelastic_(&srt, &sdb);
		px1cm = pcx;
		py1cm = pcy;
		pz1cm = pcz;
/*     minimum srt**2, note a 2.012GeV lower cutoff is used in N+N->Deuteron+pi: */
		ec = 4.0481439999999997f;
		ds = sqrt(sdb / 31.4f);
		dsr = ds + .1f;
		distce_(&i1, &i2, &dsr, &ds, &input1_1.dt, &ec, &srt, &ic, &
			px1cm, &py1cm, &pz1cm);
		if (ic == -1) {
		    goto L400;
		}
		crdbel_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock, &
			ntag, &sdb, nt, &ianti);
		++(*lcoll);
		goto L400;
/* lin-9/2008 Deuteron+Baryon elastic collisions over */

/* IF BARYON RESONANCE+BARYON RESONANCE COLLISIONS */
L444:
/* PREPARE THE EALSTIC CROSS SECTION FOR BARYON+BARYON COLLISIONS */
		cutoff = leadng_1.em1 + dpi_1.em2 + .02f;
/* AT HIGH ENERGIES THE ISOSPIN DEPENDENCE IS NEGLIGIBLE */
/* THE TOTAL CROSS SECTION IS TAKEN AS THAT OF THE PP */
		if (srt <= cutoff) {
		    goto L400;
		}
		if (srt > 2.245f) {
		    signn = pp2_(&srt);
		} else {
		    signn = 35.f / ((srt - cutoff) * 100.f + 1.f) + 20.f;
		}
		if (signn <= 0.f) {
		    goto L400;
		}
		xddin_(&pcx, &pcy, &pcz, &srt, &i1, &i2, &xinel, &sigk, &xsk1,
			 &xsk2, &xsk3, &xsk4, &xsk5);
		sig = signn + xinel;
/* Computing 2nd power */
		r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
		ec = r__1 * r__1;
		px1cm = pcx;
		py1cm = pcy;
		pz1cm = pcz;
/* lin-6/2008 Deuteron production: */
		ianti = 0;
		if (ee_1.lb[i1 - 1] < 0 && ee_1.lb[i2 - 1] < 0) {
		    ianti = 1;
		}
		sbbdm_(&srt, &sdprod, &ianti, &lbm, &xmm, &pfinal);
		sig += sdprod;
/* lin-6/2008 perturbative treatment of deuterons: */
		ipdflag = 0;
		if (para8_1.idpert == 1) {
		    ipert1 = 1;
		    sigr0 = sig;
		    dspert = sqrt(sigr0 / 3.1415926f / 10.f);
		    dsrpert = dspert + .1f;
		    distce_(&i1, &i2, &dsrpert, &dspert, &input1_1.dt, &ec, &
			    srt, &ic, &px1cm, &py1cm, &pz1cm);
		    if (ic == -1) {
			goto L367;
		    }
		    signn0 = 0.f;
		    crdd_(&irun, &px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &
			    iblock, &ntag, &signn0, &sigr0, nt, &ipert1);
/*     1          IBLOCK,NTAG,SIGNN,SIG) */
		    ipdflag = 1;
L367:
		    ipert1 = 0;
		}
		if (para8_1.idpert == 2) {
		    ipert1 = 1;
		}

		ds = sqrt(sig / 31.4f);
		dsr = ds + .1f;
		distce_(&i1, &i2, &dsr, &ds, &input1_1.dt, &ec, &srt, &ic, &
			px1cm, &py1cm, &pz1cm);
/*        IF(IC.EQ.-1) GO TO 400 */
		if (ic == -1) {
		    if (ipdflag == 1) {
			iblock = 501;
		    }
		    goto L400;
		}
/* CALCULATE KAON PRODUCTION PROBABILITY FROM NUCLEON+NUCLEON OR */
/* RESONANCE+RESONANCE COLLISIONS */
		goto L364;
/* CHECK WHAT KIND OF COLLISION HAS HAPPENED */
L364:
		kkk_1.ekaon[iss * 7 + 1] += 1;
/* for resonance+resonance */
/* lin-6/2008: */
		crdd_(&irun, &px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock, 
			&ntag, &signn, &sig, nt, &ipert1);
/*     1  IBLOCK,NTAG,SIGNN,SIG) */
		if (iblock == 0 && ipdflag == 1) {
		    iblock = 501;
		}

		if (abs(iblock) == 10) {
/* momentum of the three particles in the final state have been calculated */
/* in the crnn, go out of the loop */
		    ++(*lcoll);
		    if (iblock == 10) {
			++(*lddk);
		    } else if (iblock == -10) {
		    }
		    goto L400;
		}
/* lin-6/2008 */
/*        if(iblock .eq. 222)then */
		if (iblock == 222 || iblock == 501) {
/*    !! sp12/17/01 */
		    goto L400;
		}
		leadng_1.em1 = cc_1.e[i1 - 1];
		dpi_1.em2 = cc_1.e[i2 - 1];
		goto L440;
/* FOR PION+PION,pion+eta, eta+eta and rho(omega)+pion(rho,omega) or eta */
L777:
		px1cm = pcx;
		py1cm = pcy;
		pz1cm = pcz;
/* energy thresh for collisions */
		ec0 = leadng_1.em1 + dpi_1.em2 + .02f;
		if (srt <= ec0) {
		    goto L400;
		}
/* Computing 2nd power */
		r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
		ec = r__1 * r__1;
/* we negelect the elastic collision between mesons except that betwen */
/* two pions because of the lack of information about these collisions */
/* However, we do let them to collide inelastically to produce kaons */
/* lin-8/15/02       ppel=1.e-09 */
		ppel = 20.f;
		ipp = 1;
		if (leadng_1.lb1 < 3 || leadng_1.lb1 > 5 || dpi_1.lb2 < 3 || 
			dpi_1.lb2 > 5) {
		    goto L778;
		}
		ppxs_(&leadng_1.lb1, &dpi_1.lb2, &srt, &ppsig, &spprho, &ipp);
		ppel = ppsig;
L778:
		ppink = pipik_(&srt);
/* pi+eta and eta+eta are assumed to be the same as pipik( for pi+pi -> K+K-) */
/* estimated from Ko's paper: */
		ppink *= 2.f;
		if (leadng_1.lb1 >= 25 && dpi_1.lb2 >= 25) {
		    ppink = .6f;
		}
/* lin-2/13/03 include omega the same as rho, eta the same as pi: */
/*        if(((lb1.ge.3.and.lb1.le.5).and.(lb2.ge.25.and.lb2.le.27)) */
/*     1  .or.((lb2.ge.3.and.lb2.le.5).and.(lb1.ge.25.and.lb1.le.27))) */
		if ((leadng_1.lb1 == 0 || leadng_1.lb1 >= 3 && leadng_1.lb1 <=
			 5) && (dpi_1.lb2 >= 25 && dpi_1.lb2 <= 28) || (
			dpi_1.lb2 == 0 || dpi_1.lb2 >= 3 && dpi_1.lb2 <= 5) &&
			 (leadng_1.lb1 >= 25 && leadng_1.lb1 <= 28)) {
		    ppink = 0.f;
		    if (srt >= 1.393f) {
			ppink = .3f;
		    }
		}
/* pi pi <-> rho rho: */
		spprr_(&leadng_1.lb1, &dpi_1.lb2, &srt);
/* lin-4/03/02 pi pi <-> eta eta: */
		sppee_(&leadng_1.lb1, &dpi_1.lb2, &srt);
/* lin-4/03/02 pi pi <-> pi eta: */
		spppe_(&leadng_1.lb1, &dpi_1.lb2, &srt);
/* lin-4/03/02 rho pi <-> rho eta: */
		srpre_(&leadng_1.lb1, &dpi_1.lb2, &srt);
/* lin-4/03/02 omega pi <-> omega eta: */
		sopoe_(&leadng_1.lb1, &dpi_1.lb2, &srt);
/* lin-4/03/02 rho rho <-> eta eta: */
		srree_(&leadng_1.lb1, &dpi_1.lb2, &srt);
		ppb1_1.ppinnb = 0.f;
		if (srt > ppbmas_1.thresh[0]) {
		    getnst_(&srt);
		    if (leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5 && dpi_1.lb2 >=
			     3 && dpi_1.lb2 <= 5) {
			ppb1_1.ppinnb = ppbbar_(&srt);
		    } else if (leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5 && 
			    dpi_1.lb2 >= 25 && dpi_1.lb2 <= 27 || dpi_1.lb2 >=
			     3 && dpi_1.lb2 <= 5 && leadng_1.lb1 >= 25 && 
			    leadng_1.lb1 <= 27) {
			ppb1_1.ppinnb = prbbar_(&srt);
		    } else if (leadng_1.lb1 >= 25 && leadng_1.lb1 <= 27 && 
			    dpi_1.lb2 >= 25 && dpi_1.lb2 <= 27) {
			ppb1_1.ppinnb = rrbbar_(&srt);
		    } else if (leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5 && 
			    dpi_1.lb2 == 28 || dpi_1.lb2 >= 3 && dpi_1.lb2 <= 
			    5 && leadng_1.lb1 == 28) {
			ppb1_1.ppinnb = pobbar_(&srt);
		    } else if (leadng_1.lb1 >= 25 && leadng_1.lb1 <= 27 && 
			    dpi_1.lb2 == 28 || dpi_1.lb2 >= 25 && dpi_1.lb2 <=
			     27 && leadng_1.lb1 == 28) {
			ppb1_1.ppinnb = robbar_(&srt);
		    } else if (leadng_1.lb1 == 28 && dpi_1.lb2 == 28) {
			ppb1_1.ppinnb = oobbar_(&srt);
		    } else {
			if (leadng_1.lb1 != 0 && dpi_1.lb2 != 0) {
			    s_wsle(&io___379);
			    do_lio(&c__9, &c__1, "missed MM lb1,lb2=", (
				    ftnlen)18);
			    do_lio(&c__3, &c__1, (char *)&leadng_1.lb1, (
				    ftnlen)sizeof(integer));
			    do_lio(&c__3, &c__1, (char *)&dpi_1.lb2, (ftnlen)
				    sizeof(integer));
			    e_wsle();
			}
		    }
		}
		ppin = ppink + ppb1_1.ppinnb + ppmm_1.pprr + ppmm_1.ppee + 
			ppmm_1.pppe + ppmm_1.rpre + ppmm_1.xopoe + 
			ppmm_1.rree;
/* check if a collision can happen */
		if (ppel + ppin <= .01f) {
		    goto L400;
		}
		dspp = sqrt((ppel + ppin) / 31.4f);
		dsppr = dspp + .1f;
		distce_(&i1, &i2, &dsppr, &dspp, &input1_1.dt, &ec, &srt, &ic,
			 &px1cm, &py1cm, &pz1cm);
		if (ic == -1) {
		    goto L400;
		}
		if (ppel == 0.f) {
		    goto L400;
		}
/* the collision can happen */
/* check what kind collision has happened */
		kkk_1.ekaon[iss * 7 + 4] += 1;
		crpp_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock, &ppel, 
			&ppin, &spprho, &ipp);
/* rho formation, go to 400 */
/*       if(iblock.eq.666)go to 600 */
		if (iblock == 666) {
		    goto L555;
		}
		if (iblock == 6) {
		    ++(*lpp);
		}
		if (iblock == 66) {
		    ++(*lppk);
		} else if (iblock == 366) {
		    ++(*lppk);
		} else if (iblock == 367) {
		    ++(*lppk);
		}
		leadng_1.em1 = cc_1.e[i1 - 1];
		dpi_1.em2 = cc_1.e[i2 - 1];
		goto L440;
/* In this block we treat annihilations of */
/* lin-9/28/00* an anti-nucleon and a baryon or baryon resonance */
/* an anti-baryon and a baryon (including resonances) */
L2799:
		px1cm = pcx;
		py1cm = pcy;
		pz1cm = pcz;
/* Computing 2nd power */
		r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
		ec = r__1 * r__1;
/* lin assume the same cross section (as a function of sqrt s) as for PPbar: */
/* lin-ctest annih maximum */
/*        DSppb=SQRT(amin1(xppbar(srt),30.)/PI/10.) */
		dsppb = sqrt(xppbar_(&srt) / 3.1415926f / 10.f);
		dsppbr = dsppb + .1f;
		distce_(&i1, &i2, &dsppbr, &dsppb, &input1_1.dt, &ec, &srt, &
			ic, &px1cm, &py1cm, &pz1cm);
		if (ic == -1) {
		    goto L400;
		}
		crppba_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock);
		leadng_1.em1 = cc_1.e[i1 - 1];
		dpi_1.em2 = cc_1.e[i2 - 1];
		goto L440;

L3555:
		px1cm = pcx;
		py1cm = pcy;
		pz1cm = pcz;
/* Computing 2nd power */
		r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
		ec = r__1 * r__1;
		dskk = sqrt(sig / 3.1415926f / 10.f);
		dskk0 = dskk + .1f;
		distce_(&i1, &i2, &dskk0, &dskk, &input1_1.dt, &ec, &srt, &ic,
			 &px1cm, &py1cm, &pz1cm);
		if (ic == -1) {
		    goto L400;
		}
		crlaba_(&px1cm, &py1cm, &pz1cm, &srt, &brel, &brsgm, &i1, &i2,
			 nt, &iblock, &nchrg, &icase);
		leadng_1.em1 = cc_1.e[i1 - 1];
		dpi_1.em2 = cc_1.e[i2 - 1];
		goto L440;

/* perturbative production of cascade and omega */
L3455:
		px1cm = pcx;
		py1cm = pcy;
		pz1cm = pcz;
		pertur_(&px1cm, &py1cm, &pz1cm, &srt, &irun, &i1, &i2, nt, &
			kp, &icontp);
		if (icontp == 0) {
/*     inelastic collisions: */
		    leadng_1.em1 = cc_1.e[i1 - 1];
		    dpi_1.em2 = cc_1.e[i2 - 1];
		    iblock = 727;
		    goto L440;
		}
/*     elastic collisions: */
		if (cc_1.e[i1 - 1] == 0.f) {
		    goto L800;
		}
		if (cc_1.e[i2 - 1] == 0.f) {
		    goto L600;
		}
		goto L400;

/* * phi + N --> pi+N(D),  N(D,N*)+N(D,N*),  K+ +La */
/* * phi + D --> pi+N(D) */
L7222:
		px1cm = pcx;
		py1cm = pcy;
		pz1cm = pcz;
/* Computing 2nd power */
		r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
		ec = r__1 * r__1;
		xphib_(&leadng_1.lb1, &dpi_1.lb2, &leadng_1.em1, &dpi_1.em2, &
			srt, &xsk1, &xsk2, &xsk3, &xsk4, &xsk5, &sigp);
		dskk = sqrt(sigp / 3.1415926f / 10.f);
		dskk0 = dskk + .1f;
		distce_(&i1, &i2, &dskk0, &dskk, &input1_1.dt, &ec, &srt, &ic,
			 &px1cm, &py1cm, &pz1cm);
		if (ic == -1) {
		    goto L400;
		}
		crphib_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &xsk1, &xsk2, 
			&xsk3, &xsk4, &xsk5, &sigp, &iblock);
		leadng_1.em1 = cc_1.e[i1 - 1];
		dpi_1.em2 = cc_1.e[i2 - 1];
		goto L440;

/* * phi + M --> K+ + K* ..... */
L7444:
		px1cm = pcx;
		py1cm = pcy;
		pz1cm = pcz;
/* Computing 2nd power */
		r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
		ec = r__1 * r__1;
		phimes_(&i1, &i2, &srt, &xsk1, &xsk2, &xsk3, &xsk4, &xsk5, &
			xsk6, &xsk7, &sigphi);
		dskk = sqrt(sigphi / 3.1415926f / 10.f);
		dskk0 = dskk + .1f;
		distce_(&i1, &i2, &dskk0, &dskk, &input1_1.dt, &ec, &srt, &ic,
			 &px1cm, &py1cm, &pz1cm);
		if (ic == -1) {
		    goto L400;
		}
/* *--- */
		pzrt = bb_1.p[i1 * 3 - 1] + bb_1.p[i2 * 3 - 1];
/* Computing 2nd power */
		r__1 = bb_1.p[i1 * 3 - 3];
/* Computing 2nd power */
		r__2 = bb_1.p[i1 * 3 - 2];
/* Computing 2nd power */
		r__3 = bb_1.p[i1 * 3 - 1];
/* Computing 2nd power */
		r__4 = cc_1.e[i1 - 1];
		er1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * 
			r__4);
/* Computing 2nd power */
		r__1 = bb_1.p[i2 * 3 - 3];
/* Computing 2nd power */
		r__2 = bb_1.p[i2 * 3 - 2];
/* Computing 2nd power */
		r__3 = bb_1.p[i2 * 3 - 1];
/* Computing 2nd power */
		r__4 = cc_1.e[i2 - 1];
		er2 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * 
			r__4);
		ert = er1 + er2;
		yy = log((ert + pzrt) / (ert - pzrt)) * .5f;
/* *------ */
		crphim_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &xsk1, &xsk2, 
			&xsk3, &xsk4, &xsk5, &xsk6, &sigphi, &ikkg, &ikkl, &
			iblock);
		leadng_1.em1 = cc_1.e[i1 - 1];
		dpi_1.em2 = cc_1.e[i2 - 1];
		goto L440;

/* lambda-N elastic xsection, Li & Ko, PRC 54(1996)1897. */
L7799:
		px1cm = pcx;
		py1cm = pcy;
		pz1cm = pcz;
/* Computing 2nd power */
		r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
		ec = r__1 * r__1;
		lambar_(&i1, &i2, &srt, &siglab);
		dshn = sqrt(siglab / 3.1415926f / 10.f);
		dshnr = dshn + .1f;
		distce_(&i1, &i2, &dshnr, &dshn, &input1_1.dt, &ec, &srt, &ic,
			 &px1cm, &py1cm, &pz1cm);
		if (ic == -1) {
		    goto L400;
		}
		crhb_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock);
		leadng_1.em1 = cc_1.e[i1 - 1];
		dpi_1.em2 = cc_1.e[i2 - 1];
		goto L440;

/* * K+ + La(Si) --> Meson + B */
/* * K- + La(Si)-bar --> Meson + B-bar */
L5699:
		px1cm = pcx;
		py1cm = pcy;
		pz1cm = pcz;
/* Computing 2nd power */
		r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
		ec = r__1 * r__1;
		xkhype_(&i1, &i2, &srt, &xky1, &xky2, &xky3, &xky4, &xky5, &
			xky6, &xky7, &xky8, &xky9, &xky10, &xky11, &xky12, &
			xky13, &xky14, &xky15, &xky16, &xky17, &sigk);
		dskk = sqrt(sigk / 3.1415926f);
		dskk0 = dskk + .1f;
		distce_(&i1, &i2, &dskk0, &dskk, &input1_1.dt, &ec, &srt, &ic,
			 &px1cm, &py1cm, &pz1cm);
		if (ic == -1) {
		    goto L400;
		}

		if (ee_1.lb[i1 - 1] == 23 || ee_1.lb[i2 - 1] == 23) {
		    ikmp = 1;
		} else {
		    ikmp = -1;
		}
		crkhyp_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &xky1, &xky2, 
			&xky3, &xky4, &xky5, &xky6, &xky7, &xky8, &xky9, &
			xky10, &xky11, &xky12, &xky13, &xky14, &xky15, &xky16,
			 &xky17, &sigk, &ikmp, &iblock);
		leadng_1.em1 = cc_1.e[i1 - 1];
		dpi_1.em2 = cc_1.e[i2 - 1];
		goto L440;
/* khyperon end */

/* sp11/03/01 La/Si-bar + N --> pi + K+ */
/*  La/Si + N-bar --> pi + K- */
L5999:
		px1cm = pcx;
		py1cm = pcy;
		pz1cm = pcz;
/* Computing 2nd power */
		r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
		ec = r__1 * r__1;
		sigkp = 15.f;
/*      if((lb1.ge.14.and.lb1.le.17) */
/*     &    .or.(lb2.ge.14.and.lb2.le.17))sigkp=10. */
		dskk = sqrt(sigkp / 3.1415926f / 10.f);
		dskk0 = dskk + .1f;
		distce_(&i1, &i2, &dskk0, &dskk, &input1_1.dt, &ec, &srt, &ic,
			 &px1cm, &py1cm, &pz1cm);
		if (ic == -1) {
		    goto L400;
		}

		crlan_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock);
		leadng_1.em1 = cc_1.e[i1 - 1];
		dpi_1.em2 = cc_1.e[i2 - 1];
		goto L440;

/* * */
/* K(K*) + K(K*) --> phi + pi(rho,omega) */
L8699:
		px1cm = pcx;
		py1cm = pcy;
		pz1cm = pcz;
/* Computing 2nd power */
		r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
		ec = r__1 * r__1;
/*  CALL CROSSKKPHI(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,IBLOCK)  used for KK*->phi+rho */
		crkphi_(&px1cm, &py1cm, &pz1cm, &ec, &srt, &iblock, &emm1, &
			emm2, &lbp1, &lbp2, &i1, &i2, &ikk, &icase, &c_b122, &
			c_b123);
		if (icase == 0) {
		    iblock = 0;
		    goto L400;
		}
/* *--- */
		if (lbp1 == 29 || lbp2 == 29) {
		    pzrt = bb_1.p[i1 * 3 - 1] + bb_1.p[i2 * 3 - 1];
/* Computing 2nd power */
		    r__1 = bb_1.p[i1 * 3 - 3];
/* Computing 2nd power */
		    r__2 = bb_1.p[i1 * 3 - 2];
/* Computing 2nd power */
		    r__3 = bb_1.p[i1 * 3 - 1];
/* Computing 2nd power */
		    r__4 = cc_1.e[i1 - 1];
		    er1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 
			    * r__4);
/* Computing 2nd power */
		    r__1 = bb_1.p[i2 * 3 - 3];
/* Computing 2nd power */
		    r__2 = bb_1.p[i2 * 3 - 2];
/* Computing 2nd power */
		    r__3 = bb_1.p[i2 * 3 - 1];
/* Computing 2nd power */
		    r__4 = cc_1.e[i2 - 1];
		    er2 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 
			    * r__4);
		    ert = er1 + er2;
		    yy = log((ert + pzrt) / (ert - pzrt)) * .5f;
/* *------ */
		    iblock = 222;
		    ntag = 0;
		}
		ee_1.lb[i1 - 1] = lbp1;
		ee_1.lb[i2 - 1] = lbp2;
		cc_1.e[i1 - 1] = emm1;
		cc_1.e[i2 - 1] = emm2;
		leadng_1.em1 = cc_1.e[i1 - 1];
		dpi_1.em2 = cc_1.e[i2 - 1];
		goto L440;
/* * */
/* rho(omega) + K(K*)  --> phi + K(K*) */
L8799:
		px1cm = pcx;
		py1cm = pcy;
		pz1cm = pcz;
/* Computing 2nd power */
		r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
		ec = r__1 * r__1;
/*  CALL CROSSKKPHI(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,IBLOCK)  used for KK*->phi+rho */
		crksph_(&px1cm, &py1cm, &pz1cm, &ec, &srt, &emm1, &emm2, &
			lbp1, &lbp2, &i1, &i2, &ikkg, &ikkl, &iblock, &icase, 
			&c_b106);
		if (icase == 0) {
		    iblock = 0;
		    goto L400;
		}

		if (lbp1 == 29 || lbp2 == 20) {
/* *--- */
		    pzrt = bb_1.p[i1 * 3 - 1] + bb_1.p[i2 * 3 - 1];
/* Computing 2nd power */
		    r__1 = bb_1.p[i1 * 3 - 3];
/* Computing 2nd power */
		    r__2 = bb_1.p[i1 * 3 - 2];
/* Computing 2nd power */
		    r__3 = bb_1.p[i1 * 3 - 1];
/* Computing 2nd power */
		    r__4 = cc_1.e[i1 - 1];
		    er1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 
			    * r__4);
/* Computing 2nd power */
		    r__1 = bb_1.p[i2 * 3 - 3];
/* Computing 2nd power */
		    r__2 = bb_1.p[i2 * 3 - 2];
/* Computing 2nd power */
		    r__3 = bb_1.p[i2 * 3 - 1];
/* Computing 2nd power */
		    r__4 = cc_1.e[i2 - 1];
		    er2 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 
			    * r__4);
		    ert = er1 + er2;
		    yy = log((ert + pzrt) / (ert - pzrt)) * .5f;
		}
		ee_1.lb[i1 - 1] = lbp1;
		ee_1.lb[i2 - 1] = lbp2;
		cc_1.e[i1 - 1] = emm1;
		cc_1.e[i2 - 1] = emm2;
		leadng_1.em1 = cc_1.e[i1 - 1];
		dpi_1.em2 = cc_1.e[i2 - 1];
		goto L440;
/* for kaon+baryon scattering, using a constant xsection of 10 mb. */
L888:
		px1cm = pcx;
		py1cm = pcy;
		pz1cm = pcz;
/* Computing 2nd power */
		r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
		ec = r__1 * r__1;
		sig = 10.f;
		if (abs(leadng_1.lb1) == 14 || abs(dpi_1.lb2) == 14 || abs(
			leadng_1.lb1) == 30 || abs(dpi_1.lb2) == 30) {
		    sig = 20.f;
		}
		if (leadng_1.lb1 == 29 || dpi_1.lb2 == 29) {
		    sig = 5.f;
		}
		dskn = sqrt(sig / 3.1415926f / 10.f);
		dsknr = dskn + .1f;
		distce_(&i1, &i2, &dsknr, &dskn, &input1_1.dt, &ec, &srt, &ic,
			 &px1cm, &py1cm, &pz1cm);
		if (ic == -1) {
		    goto L400;
		}
		crkn_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock);
		leadng_1.em1 = cc_1.e[i1 - 1];
		dpi_1.em2 = cc_1.e[i2 - 1];
		goto L440;
/* ** */
L440:
/*                IBLOCK = 0 ; NOTHING HAS HAPPENED */
/*                IBLOCK = 1 ; ELASTIC N-N COLLISION */
/*                IBLOCK = 2 ; N + N -> N + DELTA */
/*                IBLOCK = 3 ; N + DELTA -> N + N */
/*                IBLOCK = 4 ; N + N -> d + d + PION,DIRECT PROCESS */
/*               IBLOCK = 5 ; D(N*)+D(N*) COLLISIONS */
/*                IBLOCK = 6 ; PION+PION COLLISIONS */
/*                iblock = 7 ; pion+nucleon-->l/s+kaon */
/*               iblock =77;  pion+nucleon-->delta+pion */
/*               iblock = 8 ; kaon+baryon rescattering */
/*                IBLOCK = 9 ; NN-->KAON+X */
/*                IBLOCK = 10; DD-->KAON+X */
/*               IBLOCK = 11; ND-->KAON+X */
/* bali2/1/99 */

/*           iblock   - 1902 annihilation-->pion(+)+pion(-)   (2 pion) */
/*           iblock   - 1903 annihilation-->pion(+)+rho(-)    (3 pion) */
/*           iblock   - 1904 annihilation-->rho(+)+rho(-)     (4 pion) */
/*           iblock   - 1905 annihilation-->rho(0)+omega      (5 pion) */
/*           iblock   - 1906 annihilation-->omega+omega       (6 pion) */
/* bali3/5/99 */
/*           iblock   - 1907 K+K- to pi+pi- */
/* bali3/5/99 end */
/* bz3/9/99 khyperon */
/*           iblock   - 1908 K+Y -> piN */
/* bz3/9/99 khyperon end */
/* bali2/1/99end */
/* lin-9/28/00 Processes: m(pi rho omega)+m(pi rho omega) */
/*     to anti-(p n D N*1 N*2)+(p n D N*1 N*2): */
/*           iblock   - 1801  mm -->pbar p */
/*           iblock   - 18021 mm -->pbar n */
/*           iblock   - 18022 mm -->nbar p */
/*           iblock   - 1803  mm -->nbar n */
/*           iblock   - 18041 mm -->pbar Delta */
/*           iblock   - 18042 mm -->anti-Delta p */
/*           iblock   - 18051 mm -->nbar Delta */
/*           iblock   - 18052 mm -->anti-Delta n */
/*           iblock   - 18061 mm -->pbar N*(1400) */
/*           iblock   - 18062 mm -->anti-N*(1400) p */
/*           iblock   - 18071 mm -->nbar N*(1400) */
/*           iblock   - 18072 mm -->anti-N*(1400) n */
/*           iblock   - 1808  mm -->anti-Delta Delta */
/*           iblock   - 18091 mm -->pbar N*(1535) */
/*           iblock   - 18092 mm -->anti-N*(1535) p */
/*           iblock   - 18101 mm -->nbar N*(1535) */
/*           iblock   - 18102 mm -->anti-N*(1535) n */
/*           iblock   - 18111 mm -->anti-Delta N*(1440) */
/*           iblock   - 18112 mm -->anti-N*(1440) Delta */
/*           iblock   - 18121 mm -->anti-Delta N*(1535) */
/*           iblock   - 18122 mm -->anti-N*(1535) Delta */
/*           iblock   - 1813  mm -->anti-N*(1440) N*(1440) */
/*           iblock   - 18141 mm -->anti-N*(1440) N*(1535) */
/*           iblock   - 18142 mm -->anti-N*(1535) N*(1440) */
/*           iblock   - 1815  mm -->anti-N*(1535) N*(1535) */
/* lin-9/28/00-end */
/* lin-10/08/00 Processes: pi pi <-> rho rho */
/*           iblock   - 1850  pi pi -> rho rho */
/*           iblock   - 1851  rho rho -> pi pi */
/* lin-10/08/00-end */
/* lin-08/14/02 Processes: pi pi <-> eta eta */
/*           iblock   - 1860  pi pi -> eta eta */
/*           iblock   - 1861  eta eta -> pi pi */
/* Processes: pi pi <-> pi eta */
/*           iblock   - 1870  pi pi -> pi eta */
/*           iblock   - 1871  pi eta -> pi pi */
/* Processes: rho pi <-> rho eta */
/*           iblock   - 1880  pi pi -> pi eta */
/*           iblock   - 1881  pi eta -> pi pi */
/* Processes: omega pi <-> omega eta */
/*           iblock   - 1890  pi pi -> pi eta */
/*           iblock   - 1891  pi eta -> pi pi */
/* Processes: rho rho <-> eta eta */
/*           iblock   - 1895  rho rho -> eta eta */
/*           iblock   - 1896  eta eta -> rho rho */
/* lin-08/14/02-end */
/* lin-11/07/00 Processes: */
/*           iblock   - 366  pi rho -> K* Kbar or K*bar K */
/*           iblock   - 466  pi rho <- K* Kbar or K*bar K */
/* lin-9/2008 Deuteron: */
/*           iblock   - 501  B+B -> Deuteron+Meson */
/*           iblock   - 502  Deuteron+Meson -> B+B */
/*           iblock   - 503  Deuteron+Baryon elastic */
/*           iblock   - 504  Deuteron+Meson elastic */

		if (iblock == 0) {
		    goto L400;
		}
/* COM: FOR DIRECT PROCESS WE HAVE TREATED THE PAULI BLOCKING AND FIND */
/*     THE MOMENTUM OF PARTICLES IN THE ''LAB'' FRAME. SO GO TO 400 */
/* A COLLISION HAS TAKEN PLACE !! */
		++(*lcoll);
/* WAS COLLISION PAULI-FORBIDEN? IF YES, NTAG = -1 */
		ntag = 0;

/*             LORENTZ-TRANSFORMATION INTO CMS FRAME */
/* Computing 2nd power */
		r__1 = leadng_1.em1;
/* Computing 2nd power */
		r__2 = px1cm;
/* Computing 2nd power */
		r__3 = py1cm;
/* Computing 2nd power */
		r__4 = pz1cm;
		e1cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * 
			r__4);
		p1beta = px1cm * bg_1.betax + py1cm * bg_1.betay + pz1cm * 
			bg_1.betaz;
		transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) 
			+ e1cm);
		pt1i1 = bg_1.betax * transf + px1cm;
		pt2i1 = bg_1.betay * transf + py1cm;
		pt3i1 = bg_1.betaz * transf + pz1cm;
/* negelect the pauli blocking at high energies */
		goto L90002;
/* lin-10/25/02-comment out following, since there is no path to it: */
/* *CHECK IF PARTICLE #1 IS PAULI BLOCKED */
/*              CALL PAULat(I1,occup) */
/*              if (RANART(NSEED) .lt. occup) then */
/*                ntag = -1 */
/*              else */
/*                ntag = 0 */
/*              end if */
/* lin-10/25/02-end */
L90002:
/* IF PARTICLE #1 IS NOT PAULI BLOCKED */
/*              IF (NTAG .NE. -1) THEN */
/* Computing 2nd power */
		r__1 = dpi_1.em2;
/* Computing 2nd power */
		r__2 = px1cm;
/* Computing 2nd power */
		r__3 = py1cm;
/* Computing 2nd power */
		r__4 = pz1cm;
		e2cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * 
			r__4);
		transf = bg_1.gamma * (-bg_1.gamma * p1beta / (bg_1.gamma + 
			1.f) + e2cm);
		pt1i2 = bg_1.betax * transf - px1cm;
		pt2i2 = bg_1.betay * transf - py1cm;
		pt3i2 = bg_1.betaz * transf - pz1cm;
		goto L90003;
/* lin-10/25/02-comment out following, since there is no path to it: */
/* *CHECK IF PARTICLE #2 IS PAULI BLOCKED */
/*                CALL PAULat(I2,occup) */
/*                if (RANART(NSEED) .lt. occup) then */
/*                  ntag = -1 */
/*                else */
/*                  ntag = 0 */
/*                end if */
/* c              END IF */
/* * IF COLLISION IS BLOCKED,RESTORE THE MOMENTUM,MASSES */
/* * AND LABELS OF I1 AND I2 */
/* c             IF (NTAG .EQ. -1) THEN */
/*                LBLOC  = LBLOC + 1 */
/*                P(1,I1) = PX1 */
/*                P(2,I1) = PY1 */
/*                P(3,I1) = PZ1 */
/*                P(1,I2) = PX2 */
/*                P(2,I2) = PY2 */
/*                P(3,I2) = PZ2 */
/*                E(I1)   = EM1 */
/*                E(I2)   = EM2 */
/*                LB(I1)  = LB1 */
/*                LB(I2)  = LB2 */
/* c              ELSE */
/* lin-10/25/02-end */
L90003:
		if (iblock == 1) {
		    ++(*lcnne);
		}
		if (iblock == 5) {
		    ++(*ldd);
		}
		if (iblock == 2) {
		    ++(*lcnnd);
		}
		if (iblock == 8) {
		    ++(*lkn);
		}
		if (iblock == 43) {
		    ++(*ldou);
		}
/*                IF(IBLOCK.EQ.2) THEN */
/* CALCULATE THE AVERAGE SRT FOR N + N---> N + DELTA PROCESS */
/*                NODELT=NODELT+1 */
/*                SUMSRT=SUMSRT+SRT */
/*                ENDIF */
		if (iblock == 3) {
		    ++(*lcndn);
		}
/* assign final momenta to particles while keep the leadng particle */
/* behaviour */
/*              if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0)then */
		bb_1.p[i1 * 3 - 3] = pt1i1;
		bb_1.p[i1 * 3 - 2] = pt2i1;
		bb_1.p[i1 * 3 - 1] = pt3i1;
		bb_1.p[i2 * 3 - 3] = pt1i2;
		bb_1.p[i2 * 3 - 2] = pt2i2;
		bb_1.p[i2 * 3 - 1] = pt3i2;
/*              else */
/*              p(1,i1)=pt1i2 */
/*              p(2,i1)=pt2i2 */
/*              p(3,i1)=pt3i2 */
/*              p(1,i2)=pt1i1 */
/*              p(2,i2)=pt2i1 */
/*              p(3,i2)=pt3i1 */
/*              endif */
		leadng_1.px1 = bb_1.p[i1 * 3 - 3];
		leadng_1.py1 = bb_1.p[i1 * 3 - 2];
		leadng_1.pz1 = bb_1.p[i1 * 3 - 1];
		leadng_1.em1 = cc_1.e[i1 - 1];
		dpi_1.em2 = cc_1.e[i2 - 1];
		leadng_1.lb1 = ee_1.lb[i1 - 1];
		dpi_1.lb2 = ee_1.lb[i2 - 1];
		ee_1.id[i1 - 1] = 2;
		ee_1.id[i2 - 1] = 2;
/* Computing 2nd power */
		r__1 = leadng_1.em1;
/* Computing 2nd power */
		r__2 = leadng_1.px1;
/* Computing 2nd power */
		r__3 = leadng_1.py1;
/* Computing 2nd power */
		r__4 = leadng_1.pz1;
		leadng_1.e1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + 
			r__4 * r__4);
		id1 = ee_1.id[i1 - 1];
		goto L90004;
/* lin-10/25/02-comment out following, since there is no path to it: */
/* * change phase space density FOR NUCLEONS INVOLVED : */
/* * NOTE THAT f is the phase space distribution function for nucleons only */
/*                if ((abs(ix1).le.mx) .and. (abs(iy1).le.my) .and. */
/*     &              (abs(iz1).le.mz)) then */
/*                  ipx1p = nint(p(1,i1)/dpx) */
/*                  ipy1p = nint(p(2,i1)/dpy) */
/*                  ipz1p = nint(p(3,i1)/dpz) */
/*                  if ((ipx1p.ne.ipx1) .or. (ipy1p.ne.ipy1) .or. */
/*     &                (ipz1p.ne.ipz1)) then */
/*                    if ((abs(ipx1).le.mpx) .and. (abs(ipy1).le.my) */
/*     &                .and. (ipz1.ge.-mpz) .and. (ipz1.le.mpzp) */
/*     &                .AND. (AM1.LT.1.)) */
/*     &                f(ix1,iy1,iz1,ipx1,ipy1,ipz1) = */
/*     &                f(ix1,iy1,iz1,ipx1,ipy1,ipz1) - 1. */
/*                    if ((abs(ipx1p).le.mpx) .and. (abs(ipy1p).le.my) */
/*     &                .and. (ipz1p.ge.-mpz).and. (ipz1p.le.mpzp) */
/*     &                .AND. (EM1.LT.1.)) */
/*     &                f(ix1,iy1,iz1,ipx1p,ipy1p,ipz1p) = */
/*     &                f(ix1,iy1,iz1,ipx1p,ipy1p,ipz1p) + 1. */
/*                  end if */
/*                end if */
/*                if ((abs(ix2).le.mx) .and. (abs(iy2).le.my) .and. */
/*     &              (abs(iz2).le.mz)) then */
/*                  ipx2p = nint(p(1,i2)/dpx) */
/*                  ipy2p = nint(p(2,i2)/dpy) */
/*                  ipz2p = nint(p(3,i2)/dpz) */
/*                  if ((ipx2p.ne.ipx2) .or. (ipy2p.ne.ipy2) .or. */
/*     &                (ipz2p.ne.ipz2)) then */
/*                    if ((abs(ipx2).le.mpx) .and. (abs(ipy2).le.my) */
/*     &                .and. (ipz2.ge.-mpz) .and. (ipz2.le.mpzp) */
/*     &                .AND. (AM2.LT.1.)) */
/*     &                f(ix2,iy2,iz2,ipx2,ipy2,ipz2) = */
/*     &                f(ix2,iy2,iz2,ipx2,ipy2,ipz2) - 1. */
/*                    if ((abs(ipx2p).le.mpx) .and. (abs(ipy2p).le.my) */
/*     &                .and. (ipz2p.ge.-mpz) .and. (ipz2p.le.mpzp) */
/*     &                .AND. (EM2.LT.1.)) */
/*     &                f(ix2,iy2,iz2,ipx2p,ipy2p,ipz2p) = */
/*     &                f(ix2,iy2,iz2,ipx2p,ipy2p,ipz2p) + 1. */
/*                  end if */
/*                end if */
/* lin-10/25/02-end */
L90004:
		am1 = leadng_1.em1;
		am2 = dpi_1.em2;
/*            END IF */
L400:

/* lin-6/10/03 skips the info output on resonance creations: */
/*            goto 550 */
/* clin-4/30/03 study phi,K*,Lambda(1520) resonances at creation: */
/* c     note that no decays give these particles, so don't need to consider nnn: */
/*            if(iblock.ne.0.and.(lb(i1).eq.29.or.iabs(lb(i1)).eq.30 */
/*     1           .or.lb(i2).eq.29.or.iabs(lb(i2)).eq.30 */
/*     2           .or.lb1i.eq.29.or.iabs(lb1i).eq.30 */
/*     3           .or.lb2i.eq.29.or.iabs(lb2i).eq.30)) then */
/*               lb1now=lb(i1) */
/*               lb2now=lb(i2) */
/* c */
/*               nphi0=0 */
/*               nksp0=0 */
/*               nksm0=0 */
/* c               nlar0=0 */
/* c               nlarbar0=0 */
/*               if(lb1i.eq.29) then */
/*                  nphi0=nphi0+1 */
/*               elseif(lb1i.eq.30) then */
/*                  nksp0=nksp0+1 */
/*               elseif(lb1i.eq.-30) then */
/*                  nksm0=nksm0+1 */
/*               endif */
/*               if(lb2i.eq.29) then */
/*                  nphi0=nphi0+1 */
/*               elseif(lb2i.eq.30) then */
/*                  nksp0=nksp0+1 */
/*               elseif(lb2i.eq.-30) then */
/*                  nksm0=nksm0+1 */
/*               endif */
/* c */
/*               nphi=0 */
/*               nksp=0 */
/*               nksm=0 */
/*               nlar=0 */
/*               nlarbar=0 */
/*               if(lb1now.eq.29) then */
/*                  nphi=nphi+1 */
/*               elseif(lb1now.eq.30) then */
/*                  nksp=nksp+1 */
/*               elseif(lb1now.eq.-30) then */
/*                  nksm=nksm+1 */
/*               endif */
/*               if(lb2now.eq.29) then */
/*                  nphi=nphi+1 */
/*               elseif(lb2now.eq.30) then */
/*                  nksp=nksp+1 */
/*               elseif(lb2now.eq.-30) then */
/*                  nksm=nksm+1 */
/*               endif */
/* c */
/*               if(nphi.eq.2.or.nksp.eq.2.or.nksm.eq.2) then */
/*                  write(91,*) '2 same resonances in one reaction!' */
/*                  write(91,*) nphi,nksp,nksm,iblock */
/*               endif */

/* c     All reactions create or destroy no more than 1 these resonance, */
/* c     otherwise file "fort.91" warns us: */
/*               do 222 ires=1,3 */
/*                  if(ires.eq.1.and.nphi.ne.nphi0) then */
/*                     idr=29 */
/*                  elseif(ires.eq.2.and.nksp.ne.nksp0) then */
/*                     idr=30 */
/*                  elseif(ires.eq.3.and.nksm.ne.nksm0) then */
/*                     idr=-30 */
/*                  else */
/*                     goto 222 */
/*                  endif */
/* ctest off for resonance (phi, K*) studies: */
/* c               if(lb1now.eq.idr) then */
/* c       write(17,112) 'collision',lb1now,P(1,I1),P(2,I1),P(3,I1),e(I1),nt */
/* c               elseif(lb2now.eq.idr) then */
/* c       write(17,112) 'collision',lb2now,P(1,I2),P(2,I2),P(3,I2),e(I2),nt */
/* c               elseif(lb1i.eq.idr) then */
/* c       write(18,112) 'collision',lb1i,px1i,py1i,pz1i,em1i,nt */
/* c               elseif(lb2i.eq.idr) then */
/* c       write(18,112) 'collision',lb2i,px2i,py2i,pz2i,em2i,nt */
/* c               endif */
/* 222           continue */

/*            else */
/*            endif */
/* c 112        format(a10,I4,4(1x,f9.3),1x,I4) */

/* lin-2/26/03 skips the check of energy conservation after each binary search: */
/* 550        goto 555 */
/*            pxfin=0 */
/*            pyfin=0 */
/*            pzfin=0 */
/*            efin=0 */
/*            if(e(i1).ne.0.or.lb(i1).eq.10022) then */
/*               efin=efin+SQRT(E(I1)**2+P(1,I1)**2+P(2,I1)**2+P(3,I1)**2) */
/*               pxfin=pxfin+P(1,I1) */
/*               pyfin=pyfin+P(2,I1) */
/*               pzfin=pzfin+P(3,I1) */
/*            endif */
/*            if(e(i2).ne.0.or.lb(i2).eq.10022) then */
/*               efin=efin+SQRT(E(I2)**2+P(1,I2)**2+P(2,I2)**2+P(3,I2)**2) */
/*               pxfin=pxfin+P(1,I2) */
/*               pyfin=pyfin+P(2,I2) */
/*               pzfin=pzfin+P(3,I2) */
/*            endif */
/*            if((nnn-nnnini).ge.1) then */
/*               do imore=nnnini+1,nnn */
/*                  if(EPION(imore,IRUN).ne.0) then */
/*                     efin=efin+SQRT(EPION(imore,IRUN)**2 */
/*     1                    +PPION(1,imore,IRUN)**2+PPION(2,imore,IRUN)**2 */
/*     2                    +PPION(3,imore,IRUN)**2) */
/*                     pxfin=pxfin+PPION(1,imore,IRUN) */
/*                     pyfin=pyfin+PPION(2,imore,IRUN) */
/*                     pzfin=pzfin+PPION(3,imore,IRUN) */
/*                  endif */
/*               enddo */
/*            endif */
/*            devio=sqrt((pxfin-pxini)**2+(pyfin-pyini)**2 */
/*     1           +(pzfin-pzini)**2+(efin-eini)**2) */
/* c */
/*            if(devio.ge.0.1) then */
/*               write(92,'a20,5(1x,i6),2(1x,f8.3)') 'iblock,lb,npi=', */
/*     1              iblock,lb1i,lb2i,lb(i1),lb(i2),e(i1),e(i2) */
/*               do imore=nnnini+1,nnn */
/*                  if(EPION(imore,IRUN).ne.0) then */
/*                     write(92,'a10,2(1x,i6)') 'ipi,lbm=', */
/*     1                    imore,LPION(imore,IRUN) */
/*                  endif */
/*               enddo */
/*               write(92,'a3,4(1x,f8.3)') 'I:',eini,pxini,pyini,pzini */
/*               write(92,'a3,5(1x,f8.3)') */
/*     1              'F:',efin,pxfin,pyfin,pzfin,devio */
/*            endif */

L555:
/* test off only one collision for the same 2 particles in the same timestep: */
/*            if(iblock.ne.0) then */
/*               goto 800 */
/*            endif */
/* test off collisions history: */
/*            if(iblock.ne.0) then */
/*               write(10,*) nt,i1,i2,iblock,x1,z1,x2,z2 */
/*            endif */
L600:
		;
	    }
/* lin-4/2012 option of pi0 decays: */
/*     particles in lpion() may be a pi0, and when ipi0dcy=1 */
/*     we need to decay them at nt=ntmax after all lb(i1) decays are done: */
L798:
	    if (*nt == *ntmax && phidcy_1.ipi0dcy == 1 && i1 == rr_1.massr[
		    irun] + msum) {
		i__3 = nn_1.nnn;
		for (ipion = 1; ipion <= i__3; ++ipion) {
		    if (pd_1.lpion[ipion + irun * 150001 - 150002] == 4) {
			wid = 7.85e-9f;
			resdec_(&i1, nt, &nn_1.nnn, &wid, &idecay, &ipion);
		    }
		}
	    }
/* test off */
/*          if(nt.eq.ntmax.and.i1.eq.(MASSR(IRUN)+MSUM)) then */
/*             do ip=1,i1 */
/*                write(98,*) lb(ip),e(ip),ip */
/*             enddo */
/*          endif */
/* lin-4/2012 option of pi0 decays-end */
L800:
	    ;
	}
/* RELABLE MESONS LEFT IN THIS RUN EXCLUDING THOSE BEING CREATED DURING */
/* THIS TIME STEP AND COUNT THE TOTAL NO. OF PARTICLES IN THIS RUN */
/* note that the first mass=mta+mpr particles are baryons */
/*        write(*,*)'I: NNN,massr ', nnn,massr(irun) */
	n0 = mass + msum;
	i__2 = rr_1.massr[irun] + msum;
	for (n = n0 + 1; n <= i__2; ++n) {
/* bz11/25/98 */
/* lin-2/19/03 lb>5000: keep particles with no LB codes in ART(photon,lepton,..): */
/*        IF(E(N).GT.0.)THEN */
	    if (cc_1.e[n - 1] > 0.f || ee_1.lb[n - 1] > 5000) {
/* bz11/25/98end */
		++nn_1.nnn;
		pa_1.rpion[(nn_1.nnn + irun * 150001) * 3 - 450006] = 
			aa_1.r__[n * 3 - 3];
		pa_1.rpion[(nn_1.nnn + irun * 150001) * 3 - 450005] = 
			aa_1.r__[n * 3 - 2];
		pa_1.rpion[(nn_1.nnn + irun * 150001) * 3 - 450004] = 
			aa_1.r__[n * 3 - 1];
/* lin-10/28/03: */
		if (*nt == *ntmax) {
		    ftpisv[nn_1.nnn + irun * 150001 - 150002] = ftmax_1.ftsv[
			    n - 1];
		    tdecay_1.tfdpi[nn_1.nnn + irun * 150001 - 150002] = 
			    tdecay_1.tfdcy[n - 1];
		}

		pb_1.ppion[(nn_1.nnn + irun * 150001) * 3 - 450006] = bb_1.p[
			n * 3 - 3];
		pb_1.ppion[(nn_1.nnn + irun * 150001) * 3 - 450005] = bb_1.p[
			n * 3 - 2];
		pb_1.ppion[(nn_1.nnn + irun * 150001) * 3 - 450004] = bb_1.p[
			n * 3 - 1];
		pc_1.epion[nn_1.nnn + irun * 150001 - 150002] = cc_1.e[n - 1];
		pd_1.lpion[nn_1.nnn + irun * 150001 - 150002] = ee_1.lb[n - 1]
			;
/*       !! sp 12/19/00 */
		pe_1.propi[nn_1.nnn + irun * 150001 - 150002] = hh_1.proper[n 
			- 1];
/* lin-5/2008: */
		dpert_1.dppion[nn_1.nnn + irun * 150001 - 150002] = 
			dpert_1.dpertp[n - 1];
/*        if(lb(n) .eq. 45) */
/*    &   write(*,*)'IN-1  NT,NNN,LB,P ',nt,NNN,lb(n),proper(n) */
	    }
/* L1005: */
	}
	massrn[irun] = nn_1.nnn + mass;
/*        write(*,*)'F: NNN,massrn ', nnn,massrn(irun) */
/* L1000: */
    }
/* CALCULATE THE AVERAGE SRT FOR N + N--->N +DELTA PROCESSES */
/*        IF(NODELT.NE.0)THEN */
/*        AVSRT=SUMSRT/FLOAT(NODELT) */
/*        ELSE */
/*        AVSRT=0. */
/*        ENDIF */
/*        WRITE(1097,'(F8.2,2X,E10.3)')FLOAT(NT)*DT,AVSRT */
/* RELABLE ALL THE PARTICLES EXISTING AFTER THIS TIME STEP */
    ia = 0;
    ib = 0;
    i__1 = run_1.num;
    for (irun = 1; irun <= i__1; ++irun) {
	ia += rr_1.massr[irun - 1];
	ib += massrn[irun - 1];
	i__2 = massrn[irun];
	for (ic = 1; ic <= i__2; ++ic) {
	    ie = ia + ic;
	    ig = ib + ic;
	    if (ic <= mass) {
		rt[ig * 3 - 3] = aa_1.r__[ie * 3 - 3];
		rt[ig * 3 - 2] = aa_1.r__[ie * 3 - 2];
		rt[ig * 3 - 1] = aa_1.r__[ie * 3 - 1];
/* lin-10/28/03: */
		if (*nt == *ntmax) {
		    fttemp[ig - 1] = ftmax_1.ftsv[ie - 1];
		    tdecay_1.tft[ig - 1] = tdecay_1.tfdcy[ie - 1];
		}

		pt[ig * 3 - 3] = bb_1.p[ie * 3 - 3];
		pt[ig * 3 - 2] = bb_1.p[ie * 3 - 2];
		pt[ig * 3 - 1] = bb_1.p[ie * 3 - 1];
		et[ig - 1] = cc_1.e[ie - 1];
		lt[ig - 1] = ee_1.lb[ie - 1];
		prot[ig - 1] = hh_1.proper[ie - 1];
/* lin-5/2008: */
		dptemp[ig - 1] = dpert_1.dpertp[ie - 1];
	    } else {
		i0 = ic - mass;
		rt[ig * 3 - 3] = pa_1.rpion[(i0 + irun * 150001) * 3 - 450006]
			;
		rt[ig * 3 - 2] = pa_1.rpion[(i0 + irun * 150001) * 3 - 450005]
			;
		rt[ig * 3 - 1] = pa_1.rpion[(i0 + irun * 150001) * 3 - 450004]
			;
/* lin-10/28/03: */
		if (*nt == *ntmax) {
		    fttemp[ig - 1] = ftpisv[i0 + irun * 150001 - 150002];
		    tdecay_1.tft[ig - 1] = tdecay_1.tfdpi[i0 + irun * 150001 
			    - 150002];
		}

		pt[ig * 3 - 3] = pb_1.ppion[(i0 + irun * 150001) * 3 - 450006]
			;
		pt[ig * 3 - 2] = pb_1.ppion[(i0 + irun * 150001) * 3 - 450005]
			;
		pt[ig * 3 - 1] = pb_1.ppion[(i0 + irun * 150001) * 3 - 450004]
			;
		et[ig - 1] = pc_1.epion[i0 + irun * 150001 - 150002];
		lt[ig - 1] = pd_1.lpion[i0 + irun * 150001 - 150002];
		prot[ig - 1] = pe_1.propi[i0 + irun * 150001 - 150002];
/* lin-5/2008: */
		dptemp[ig - 1] = dpert_1.dppion[i0 + irun * 150001 - 150002];
	    }
/* L10001: */
	}
    }

    il = 0;
/* lin-10/26/01-hbt: */
/*        DO 10002 IRUN=1,NUM */
    i__2 = run_1.num;
    for (irun = 1; irun <= i__2; ++irun) {
	rr_1.massr[irun] = massrn[irun];
	il += rr_1.massr[irun - 1];
	i__1 = rr_1.massr[irun];
	for (im = 1; im <= i__1; ++im) {
	    in = il + im;
	    aa_1.r__[in * 3 - 3] = rt[in * 3 - 3];
	    aa_1.r__[in * 3 - 2] = rt[in * 3 - 2];
	    aa_1.r__[in * 3 - 1] = rt[in * 3 - 1];
/* lin-10/28/03: */
	    if (*nt == *ntmax) {
		ftmax_1.ftsv[in - 1] = fttemp[in - 1];
		tdecay_1.tfdcy[in - 1] = tdecay_1.tft[in - 1];
	    }
	    bb_1.p[in * 3 - 3] = pt[in * 3 - 3];
	    bb_1.p[in * 3 - 2] = pt[in * 3 - 2];
	    bb_1.p[in * 3 - 1] = pt[in * 3 - 1];
	    cc_1.e[in - 1] = et[in - 1];
	    ee_1.lb[in - 1] = lt[in - 1];
	    hh_1.proper[in - 1] = prot[in - 1];
/* lin-5/2008: */
	    dpert_1.dpertp[in - 1] = dptemp[in - 1];
	    if (ee_1.lb[in - 1] < 1 || ee_1.lb[in - 1] > 2) {
		ee_1.id[in - 1] = 0;
	    }
/* L10002: */
	}
/* lin-ctest off check energy conservation after each timestep */
/*         enetot=0. */
/*         do ip=1,MASSR(IRUN) */
/*            if(e(ip).ne.0.or.lb(ip).eq.10022) enetot=enetot */
/*     1           +sqrt(p(1,ip)**2+p(2,ip)**2+p(3,ip)**2+e(ip)**2) */
/*         enddo */
/*         write(91,*) 'B:',nt,enetot,massr(irun),bimp */
/* lin-3/2009 move to the end of a timestep to take care of freezeout spacetime: */
/*        call hbtout(MASSR(IRUN),nt,ntmax) */
/* L10003: */
    }

    return 0;
} /* relcol_ */

/* lin-9/2012: use double precision for S in CMS(): to avoid crash */
/*     (segmentation fault due to s<0, which happened at high energies */
/*     such as LHC with large NTMAX for two almost-comoving hadrons */
/*     that have small Pt but large |Pz|): */
/* *************************************** */
/*            SUBROUTINE CMS(I1,I2,PX1CM,PY1CM,PZ1CM,SRT) */
/* PURPOSE : FIND THE MOMENTA OF PARTICLES IN THE CMS OF THE */
/*          TWO COLLIDING PARTICLES */
/* VARIABLES : */
/* **************************************** */
/*            PARAMETER (MAXSTR=150001) */
/*            COMMON   /AA/  R(3,MAXSTR) */
/* cc      SAVE /AA/ */
/*            COMMON   /BB/  P(3,MAXSTR) */
/* cc      SAVE /BB/ */
/*            COMMON   /CC/  E(MAXSTR) */
/* cc      SAVE /CC/ */
/*            COMMON   /BG/  BETAX,BETAY,BETAZ,GAMMA */
/* cc      SAVE /BG/ */
/*            SAVE */
/*            PX1=P(1,I1) */
/*            PY1=P(2,I1) */
/*            PZ1=P(3,I1) */
/*            PX2=P(1,I2) */
/*            PY2=P(2,I2) */
/*            PZ2=P(3,I2) */
/*            EM1=E(I1) */
/*            EM2=E(I2) */
/*            E1=SQRT(EM1**2+PX1**2+PY1**2+PZ1**2) */
/*            E2=SQRT(EM2**2 + PX2**2 + PY2**2 + PZ2**2 ) */
/*            S=(E1+E2)**2-(PX1+PX2)**2-(PY1+PY2)**2-(PZ1+PZ2)**2 */
/*            SRT=SQRT(S) */
/* *LORENTZ-TRANSFORMATION IN I1-I2-C.M. SYSTEM */
/*              ETOTAL = E1 + E2 */
/*              BETAX  = (PX1+PX2) / ETOTAL */
/*              BETAY  = (PY1+PY2) / ETOTAL */
/*              BETAZ  = (PZ1+PZ2) / ETOTAL */
/*              GAMMA  = 1.0 / SQRT(1.0-BETAX**2-BETAY**2-BETAZ**2) */
/* *TRANSFORMATION OF MOMENTA (PX1CM = - PX2CM) */
/*              P1BETA = PX1*BETAX + PY1*BETAY + PZ1 * BETAZ */
/*              TRANSF = GAMMA * ( GAMMA * P1BETA / (GAMMA + 1) - E1 ) */
/*              PX1CM  = BETAX * TRANSF + PX1 */
/*              PY1CM  = BETAY * TRANSF + PY1 */
/*              PZ1CM  = BETAZ * TRANSF + PZ1 */
/*              RETURN */
/*              END */
/* *************************************** */
/* Subroutine */ int cms_(integer *i1, integer *i2, real *px1cm, real *py1cm, 
	real *pz1cm, real *srt)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static doublereal s, e1, e2, em1, em2, px1, py1, pz1, px2, py2, pz2, 
	    p1beta, dgamma, scheck, dbetax, dbetay, dbetaz, etotal, transf;

    /* Fortran I/O blocks */
    static cilist io___466 = { 0, 99, 0, 0, 0 };


/* PURPOSE : FIND THE MOMENTA OF PARTICLES IN THE CMS OF THE */
/*          TWO COLLIDING PARTICLES */
/* VARIABLES : */
/* **************************************** */
    px1 = (doublereal) bb_1.p[*i1 * 3 - 3];
    py1 = (doublereal) bb_1.p[*i1 * 3 - 2];
    pz1 = (doublereal) bb_1.p[*i1 * 3 - 1];
    px2 = (doublereal) bb_1.p[*i2 * 3 - 3];
    py2 = (doublereal) bb_1.p[*i2 * 3 - 2];
    pz2 = (doublereal) bb_1.p[*i2 * 3 - 1];
    em1 = (doublereal) cc_1.e[*i1 - 1];
    em2 = (doublereal) cc_1.e[*i2 - 1];
/* Computing 2nd power */
    d__1 = em1;
/* Computing 2nd power */
    d__2 = px1;
/* Computing 2nd power */
    d__3 = py1;
/* Computing 2nd power */
    d__4 = pz1;
    e1 = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4);
/* Computing 2nd power */
    d__1 = em2;
/* Computing 2nd power */
    d__2 = px2;
/* Computing 2nd power */
    d__3 = py2;
/* Computing 2nd power */
    d__4 = pz2;
    e2 = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4);
/* Computing 2nd power */
    d__1 = e1 + e2;
/* Computing 2nd power */
    d__2 = px1 + px2;
/* Computing 2nd power */
    d__3 = py1 + py2;
/* Computing 2nd power */
    d__4 = pz1 + pz2;
    s = d__1 * d__1 - d__2 * d__2 - d__3 * d__3 - d__4 * d__4;
    if (s <= 0.) {
	s = 0.;
    }
    *srt = (real) sqrt(s);
/* LORENTZ-TRANSFORMATION IN I1-I2-C.M. SYSTEM */
    etotal = e1 + e2;
    dbetax = (px1 + px2) / etotal;
    dbetay = (py1 + py2) / etotal;
    dbetaz = (pz1 + pz2) / etotal;
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    d__1 = dbetax;
/* Computing 2nd power */
    d__2 = dbetay;
/* Computing 2nd power */
    d__3 = dbetaz;
    scheck = 1. - d__1 * d__1 - d__2 * d__2 - d__3 * d__3;
    if (scheck <= 0.) {
	s_wsle(&io___466);
	do_lio(&c__9, &c__1, "scheck1: ", (ftnlen)9);
	do_lio(&c__5, &c__1, (char *)&scheck, (ftnlen)sizeof(doublereal));
	e_wsle();
	s_stop("", (ftnlen)0);
    }
    dgamma = 1. / sqrt(scheck);
/* TRANSFORMATION OF MOMENTA (PX1CM = - PX2CM) */
    p1beta = px1 * dbetax + py1 * dbetay + pz1 * dbetaz;
    transf = dgamma * (dgamma * p1beta / (dgamma + 1.) - e1);
    *px1cm = (real) (dbetax * transf + px1);
    *py1cm = (real) (dbetay * transf + py1);
    *pz1cm = (real) (dbetaz * transf + pz1);
    bg_1.betax = (real) dbetax;
    bg_1.betay = (real) dbetay;
    bg_1.betaz = (real) dbetaz;
    bg_1.gamma = (real) dgamma;
    return 0;
} /* cms_ */

/* lin-9/2012-end */
/* ************************************** */
/* Subroutine */ int distce_(integer *i1, integer *i2, real *deltar, real *ds,
	 real *dt, real *ec, real *srt, integer *ic, real *px1cm, real *py1cm,
	 real *pz1cm)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real s, e2, x1, y1, z1, x2, y2, z2, px2, py2, pz2, bbb, ddd, dzz, 
	    drcm, dxcm, dycm, dzcm, prcm, p1beta, drbeta, relvel, rsqare, 
	    transf;

/* PURPOSE : CHECK IF THE COLLISION BETWEEN TWO PARTICLES CAN HAPPEN */
/*           BY CHECKING */
/*                      (1) IF THE DISTANCE BETWEEN THEM IS SMALLER */
/*           THAN THE MAXIMUM DISTANCE DETERMINED FROM THE CROSS SECTION. */
/*                      (2) IF PARTICLE WILL PASS EACH OTHER WITHIN */
/*           TWO HARD CORE RADIUS. */
/*                      (3) IF PARTICLES WILL GET CLOSER. */
/* VARIABLES : */
/*           IC=1 COLLISION HAPPENED */
/*           IC=-1 COLLISION CAN NOT HAPPEN */
/* **************************************** */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /BG/ */
    *ic = 0;
    x1 = aa_1.r__[*i1 * 3 - 3];
    y1 = aa_1.r__[*i1 * 3 - 2];
    z1 = aa_1.r__[*i1 * 3 - 1];
    leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
    leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
    leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
    x2 = aa_1.r__[*i2 * 3 - 3];
    y2 = aa_1.r__[*i2 * 3 - 2];
    z2 = aa_1.r__[*i2 * 3 - 1];
    px2 = bb_1.p[*i2 * 3 - 3];
    py2 = bb_1.p[*i2 * 3 - 2];
    pz2 = bb_1.p[*i2 * 3 - 1];
    leadng_1.em1 = cc_1.e[*i1 - 1];
    dpi_1.em2 = cc_1.e[*i2 - 1];
/* Computing 2nd power */
    r__1 = leadng_1.em1;
/* Computing 2nd power */
    r__2 = leadng_1.px1;
/* Computing 2nd power */
    r__3 = leadng_1.py1;
/* Computing 2nd power */
    r__4 = leadng_1.pz1;
    leadng_1.e1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
/*            IF (ABS(X1-X2) .GT. DELTAR) GO TO 400 */
/*            IF (ABS(Y1-Y2) .GT. DELTAR) GO TO 400 */
/*            IF (ABS(Z1-Z2) .GT. DELTAR) GO TO 400 */
/* Computing 2nd power */
    r__1 = x1 - x2;
/* Computing 2nd power */
    r__2 = y1 - y2;
/* Computing 2nd power */
    r__3 = z1 - z2;
    rsqare = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
/* Computing 2nd power */
    r__1 = *deltar;
    if (rsqare > r__1 * r__1) {
	goto L400;
    }
/* NOW PARTICLES ARE CLOSE ENOUGH TO EACH OTHER ! */
/* Computing 2nd power */
    r__1 = dpi_1.em2;
/* Computing 2nd power */
    r__2 = px2;
/* Computing 2nd power */
    r__3 = py2;
/* Computing 2nd power */
    r__4 = pz2;
    e2 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    s = *srt * *srt;
    if (s < *ec) {
	goto L400;
    }
/* NOW THERE IS ENOUGH ENERGY AVAILABLE ! */
/* LORENTZ-TRANSFORMATION IN I1-I2-C.M. SYSTEM */
/* BETAX, BETAY, BETAZ AND GAMMA HAVE BEEN GIVEN IN THE SUBROUTINE CMS */
/* TRANSFORMATION OF MOMENTA (PX1CM = - PX2CM) */
    p1beta = leadng_1.px1 * bg_1.betax + leadng_1.py1 * bg_1.betay + 
	    leadng_1.pz1 * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) - 
	    leadng_1.e1);
/* Computing 2nd power */
    r__1 = *px1cm;
/* Computing 2nd power */
    r__2 = *py1cm;
/* Computing 2nd power */
    r__3 = *pz1cm;
    prcm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
    if (prcm <= 1e-5f) {
	goto L400;
    }
/* TRANSFORMATION OF SPATIAL DISTANCE */
    drbeta = bg_1.betax * (x1 - x2) + bg_1.betay * (y1 - y2) + bg_1.betaz * (
	    z1 - z2);
    transf = bg_1.gamma * bg_1.gamma * drbeta / (bg_1.gamma + 1);
    dxcm = bg_1.betax * transf + x1 - x2;
    dycm = bg_1.betay * transf + y1 - y2;
    dzcm = bg_1.betaz * transf + z1 - z2;
/* DETERMINING IF THIS IS THE POINT OF CLOSEST APPROACH */
/* Computing 2nd power */
    r__1 = dxcm;
/* Computing 2nd power */
    r__2 = dycm;
/* Computing 2nd power */
    r__3 = dzcm;
    drcm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
    dzz = (*px1cm * dxcm + *py1cm * dycm + *pz1cm * dzcm) / prcm;
/* Computing 2nd power */
    r__1 = drcm;
/* Computing 2nd power */
    r__2 = dzz;
    if (r__1 * r__1 - r__2 * r__2 <= 0.f) {
	bbb = 0.f;
    } else {
/* Computing 2nd power */
	r__1 = drcm;
/* Computing 2nd power */
	r__2 = dzz;
	bbb = sqrt(r__1 * r__1 - r__2 * r__2);
    }
/* WILL PARTICLE PASS EACH OTHER WITHIN 2 * HARD CORE RADIUS ? */
    if (bbb > *ds) {
	goto L400;
    }
    relvel = prcm * (1.f / leadng_1.e1 + 1.f / e2);
    ddd = relvel * *dt * .5f;
/* WILL PARTICLES GET CLOSER ? */
    if (dabs(ddd) < dabs(dzz)) {
	goto L400;
    }
    *ic = 1;
    goto L500;
L400:
    *ic = -1;
L500:
    return 0;
} /* distce_ */

/* *************************************** */
/*                                                                      * */
/*                                                                      * */
/* Subroutine */ int crnn_(integer *irun, real *px, real *py, real *pz, real *
	srt, integer *i1, integer *i2, integer *iblock, integer *ntag, real *
	signn, real *sig, integer *nt, integer *ipert1)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    real r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal), log(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double atan2(doublereal, doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    extern /* Subroutine */ int bbdangle_(real *, real *, real *, integer *, 
	    integer *, integer *, integer *, real *, real *, integer *);
    static real a, x, c1, c2, s1, t1, t2, s2, x1;
    static integer ic;
    static real al;
    static integer m12, n12;
    static real dm, fm, as, ta, es, pr, ss, cc1, ak0;
    extern doublereal fd5_(real *, real *, real *);
    static real dm1, dm2, dm3, dm4;
    static integer id1;
    static real ct1, ct2, pr2, st1, px3, py3, pz3, px4, py4, pz4, pz2, st2, 
	    ada;
    extern doublereal fde_(real *, real *, real *);
    static real ana;
    static integer lbd;
    extern doublereal ang_(real *, integer *);
    static integer lbm;
    extern /* Subroutine */ int n1535_(integer *, integer *, real *, real *);
    static real akp, ppd[30000]	/* was [3][10000] */, x1535;
    extern doublereal fns_(real *, real *, real *);
    static real pxd, pyd, xmm, pzd;
    extern doublereal ptr_(real *, integer *);
    static real ppx, ppy, ppz, e1cm, e2cm;
    static integer lbi1, lbi2;
    extern /* Subroutine */ int ddp2_(real *, integer *, real *, real *, real 
	    *, real *, real *, real *, real *, real *, real *, real *, real *,
	     integer *);
    static real dm1n, dm2n, sig1, sig3, sig4, sig2, eti1, eti2;
    extern doublereal ppk0_(real *), ppk1_(real *);
    static real s4pi, pxi1;
    extern doublereal x2pi_(real *), x3pi_(real *), x4pi_(real *);
    static real pyi1, xsk1, xsk2, xsk3, xsk4, xsk5, pzi1, pxi2, pyi2, pzi2;
    static integer lbpd[10000];
    static real epcm, dmin__, dmax__, arho, sigk, pt1i1, pt2i1, pt3i1, pt1i2;
    extern doublereal x33pi_(real *);
    static real srho, xdir, pt2i2, pt3i2;
    extern doublereal xrho_(real *);
    static real e1dcm, xptr, t1dlk, t2dlk;
    static integer icou1;
    static real t1dsk, t2dsk, t1nlk, t1nsk, t2nsk;
    static integer ntry1, ntry2;
    extern /* Subroutine */ int sbbdm_(real *, real *, integer *, integer *, 
	    real *, real *);
    extern doublereal omega_(real *), sigma_(real *, integer *, integer *, 
	    integer *);
    extern /* Subroutine */ int ddrho_(real *, integer *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, integer *);
    static integer ianti;
    static real pmdlk, signd, amrho, pmdsk;
    extern doublereal pplpk_(real *);
    static real pmnsk;
    extern /* Subroutine */ int pprho_(real *, integer *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, integer *);
    static real xmass, p1beta, p2beta, e2picm, dprob1, pmdlk2, pmdsk2, pmnsk2,
	     aomega;
    extern /* Subroutine */ int bbkaon_(integer *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, integer *);
    static real scheck, pfinal;
    extern /* Subroutine */ int rmasdd_(real *, real *, real *, real *, real *
	    , integer *, integer *, real *, real *);
    static real somega, ppbeta, xfinal;
    extern /* Subroutine */ int ppomga_(real *, integer *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, integer *);
    extern doublereal ranart_(integer *);
    static real sdprod, xdmass;
    extern /* Subroutine */ int rotate_(real *, real *, real *, real *, real *
	    , real *);
    static real transf;
    static integer ndloop, idloop, ipertd;
    static real p1dbeta, p2pibeta;

    /* Fortran I/O blocks */
    static cilist io___560 = { 0, 99, 0, 0, 0 };
    static cilist io___620 = { 0, 99, 0, 0, 0 };


/*     PURPOSE:                                                         * */
/*             DEALING WITH NUCLEON-NUCLEON COLLISIONS                    * */
/*     NOTE   :                                                         * */
/*     QUANTITIES:                                                 * */
/*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME* */
/*           SRT      - SQRT OF S                                       * */
/*           NSTAR =1 INCLUDING N* RESORANCE,ELSE NOT                   * */
/*           NDIRCT=1 INCLUDING DIRECT PION PRODUCTION PROCESS         * */
/*           IBLOCK   - THE INFORMATION BACK                            * */
/*                      0-> COLLISION CANNOT HAPPEN                     * */
/*                      1-> N-N ELASTIC COLLISION                       * */
/*                      2-> N+N->N+DELTA,OR N+N->N+N* REACTION          * */
/*                      3-> N+DELTA->N+N OR N+N*->N+N REACTION          * */
/*                      4-> N+N->D+D+pion reaction */
/*                     43->N+N->D(N*)+D(N*) reaction */
/*                     44->N+N->D+D+rho reaction */
/*                     45->N+N->N+N+rho */
/*                     46->N+N->N+N+omega */
/*           N12       - IS USED TO SPECIFY BARYON-BARYON REACTION      * */
/*                      CHANNELS. M12 IS THE REVERSAL CHANNEL OF N12    * */
/*                      N12,                                            * */
/*                      M12=1 FOR p+n-->delta(+)+ n                     * */
/*                          2     p+n-->delta(0)+ p                     * */
/*                          3     p+p-->delta(++)+n                     * */
/*                          4     p+p-->delta(+)+p                      * */
/*                          5     n+n-->delta(0)+n                      * */
/*                          6     n+n-->delta(-)+p                      * */
/*                          7     n+p-->N*(0)(1440)+p                   * */
/*                          8     n+p-->N*(+)(1440)+n                   * */
/*                        9     p+p-->N*(+)(1535)+p                     * */
/*                        10    n+n-->N*(0)(1535)+n                     * */
/*                         11    n+p-->N*(+)(1535)+n                     * */
/*                        12    n+p-->N*(0)(1535)+p */
/*                        13    D(++)+D(-)-->N*(+)(1440)+n */
/*                         14    D(++)+D(-)-->N*(0)(1440)+p */
/*                        15    D(+)+D(0)--->N*(+)(1440)+n */
/*                        16    D(+)+D(0)--->N*(0)(1440)+p */
/*                        17    D(++)+D(0)-->N*(+)(1535)+p */
/*                        18    D(++)+D(-)-->N*(0)(1535)+p */
/*                        19    D(++)+D(-)-->N*(+)(1535)+n */
/*                        20    D(+)+D(+)-->N*(+)(1535)+p */
/*                        21    D(+)+D(0)-->N*(+)(1535)+n */
/*                        22    D(+)+D(0)-->N*(0)(1535)+p */
/*                        23    D(+)+D(-)-->N*(0)(1535)+n */
/*                        24    D(0)+D(0)-->N*(0)(1535)+n */
/*                          25    N*(+)(14)+N*(+)(14)-->N*(+)(15)+p */
/*                          26    N*(0)(14)+N*(0)(14)-->N*(0)(15)+n */
/*                          27    N*(+)(14)+N*(0)(14)-->N*(+)(15)+n */
/*                        28    N*(+)(14)+N*(0)(14)-->N*(0)(15)+p */
/*                        29    N*(+)(14)+D+-->N*(+)(15)+p */
/*                        30    N*(+)(14)+D0-->N*(+)(15)+n */
/*                        31    N*(+)(14)+D(-)-->N*(0)(1535)+n */
/*                        32    N*(0)(14)+D++--->N*(+)(15)+p */
/*                        33    N*(0)(14)+D+--->N*(+)(15)+n */
/*                        34    N*(0)(14)+D+--->N*(0)(15)+p */
/*                        35    N*(0)(14)+D0-->N*(0)(15)+n */
/*                        36    N*(+)(14)+D0--->N*(0)(15)+p */
/*                        ++    see the note book for more listing */


/*     NOTE ABOUT N*(1440) RESORANCE IN Nucleon+NUCLEON COLLISION:      * */
/*     As it has been discussed in VerWest's paper,I= 1(initial isospin)* */
/*     channel can all be attributed to delta resorance while I= 0      * */
/*     channel can all be  attribured to N* resorance.Only in n+p       * */
/*     one can have I=0 channel so is the N*(1440) resonance            * */
/*                                                                      * */
/*                             REFERENCES:                            * */
/*                    J. CUGNON ET AL., NUCL. PHYS. A352, 505 (1981)    * */
/*                    Y. KITAZOE ET AL., PHYS. LETT. 166B, 35 (1986)    * */
/*                    B. VerWest el al., PHYS. PRV. C25 (1982)1979      * */
/*                    Gy. Wolf  et al, Nucl Phys A517 (1990) 615;       * */
/*                                     Nucl phys A552 (1993) 349.       * */
/* ********************************* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /ff/ */
/* c      SAVE /gg/ */
/* c      SAVE /INPUT/ */
/* c      SAVE /NN/ */
/* c      SAVE /BG/ */
/* c      SAVE /RUN/ */
/* c      SAVE /PA/ */
/* c      SAVE /PB/ */
/* c      SAVE /PC/ */
/* c      SAVE /PD/ */
/* c      SAVE /TABLE/ */
/* c      SAVE /input1/ */
/* c      SAVE /leadng/ */
/* c      SAVE /RNDF77/ */
/* ----------------------------------------------------------------------- */
    n12 = 0;
    m12 = 0;
    *iblock = 0;
    *ntag = 0;
    leadng_1.em1 = cc_1.e[*i1 - 1];
    dpi_1.em2 = cc_1.e[*i2 - 1];
/* Computing 2nd power */
    r__1 = *px;
/* Computing 2nd power */
    r__2 = *py;
/* Computing 2nd power */
    r__3 = *pz;
    pr = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
    c2 = *pz / pr;
    x1 = ranart_(&rndf77_1.nseed);
    ianti = 0;
    if (ee_1.lb[*i1 - 1] < 0 && ee_1.lb[*i2 - 1] < 0) {
	ianti = 1;
    }
    sbbdm_(srt, &sdprod, &ianti, &lbm, &xmm, &pfinal);
/* lin-5/2008 Production of perturbative deuterons for idpert=1: */
    if (para8_1.idpert == 1 && *ipert1 == 1) {
	if (*srt < 2.012f) {
	    return 0;
	}
	if (((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1 || (i__2 = ee_1.lb[*i1 
		- 1], abs(i__2)) == 2) && ((i__3 = ee_1.lb[*i2 - 1], abs(i__3)
		) == 1 || (i__4 = ee_1.lb[*i2 - 1], abs(i__4)) == 2)) {
	    goto L108;
	} else {
	    return 0;
	}
    }

/* ----------------------------------------------------------------------- */
/* COM: TEST FOR ELASTIC SCATTERING (EITHER N-N OR DELTA-DELTA 0R */
/*      N-DELTA OR N*-N* or N*-Delta) */
/*      IF (X1 .LE. SIGNN/SIG) THEN */
    if (x1 <= *signn / *sig) {
/* COM:  PARAMETRISATION IS TAKEN FROM THE CUGNON-PAPER */
/* Computing 6th power */
	r__1 = (*srt - 1.8766f) * 3.65f, r__1 *= r__1;
	as = r__1 * (r__1 * r__1);
	a = as * 6.f / (as + 1.f);
/* Computing 2nd power */
	r__1 = pr;
	ta = r__1 * r__1 * -2.f;
	x = ranart_(&rndf77_1.nseed);
/* lin-10/24/02        T1  = DLOG( (1-X) * DEXP(dble(A)*dble(TA)) + X )  /  A */
	t1 = (real) log((doublereal) (1.f - x) * exp((doublereal) a * (
		doublereal) ta) + (doublereal) x) / a;
	c1 = 1.f - t1 / ta;
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	*iblock = 1;
	goto L107;
    } else {
/* COM: TEST FOR INELASTIC SCATTERING */
/*     IF THE AVAILABLE ENERGY IS LESS THAN THE PION-MASS, NOTHING */
/*     CAN HAPPEN ANY MORE ==> RETURN (2.012 = 2*AVMASS + PI-MASS) */
/* lin-5/2008: Mdeuteron+Mpi=2.0106 to 2.0152 GeV/c2, so we can still use this: */
	if (*srt < 2.012f) {
	    return 0;
	}
/*     calculate the N*(1535) production cross section in N+N collisions */
/*     note that the cross sections in this subroutine are in units of mb */
/*     as only ratios of the cross sections are used to determine the */
/*     reaction channels */
	i__3 = (i__1 = ee_1.lb[*i1 - 1], abs(i__1));
	i__4 = (i__2 = ee_1.lb[*i2 - 1], abs(i__2));
	n1535_(&i__3, &i__4, srt, &x1535);
/* COM: HERE WE HAVE A PROCESS N+N ==> N+DELTA,OR N+N==>N+N*(144) or N*(1535) */
/*     OR */
/* 3 pi channel : N+N==>d1+d2+PION */
	sig3 = (x3pi_(srt) + x33pi_(srt)) * 3.f;
/* 2 pi channel : N+N==>d1+d2+d1*n*+n*n* */
	sig4 = x2pi_(srt) * 4.f;
/* 4 pi channel : N+N==>d1+d2+rho */
	s4pi = x4pi_(srt);
/* N+N-->NN+rho channel */
	srho = xrho_(srt);
/* N+N-->NN+omega */
	somega = omega_(srt);
/* CROSS SECTION FOR KAON PRODUCTION from the four channels */
/* for NLK channel */
	akp = .498f;
	ak0 = .498f;
	ana = .94f;
	ada = 1.232f;
	al = 1.1157f;
	as = 1.1197f;
	xsk1 = 0.f;
	xsk2 = 0.f;
	xsk3 = 0.f;
	xsk4 = 0.f;
	xsk5 = 0.f;
	t1nlk = ana + al + akp;
	if (*srt <= t1nlk) {
	    goto L222;
	}
	xsk1 = pplpk_(srt) * 1.5f;
/* for DLK channel */
	t1dlk = ada + al + akp;
	t2dlk = ada + al - akp;
	if (*srt <= t1dlk) {
	    goto L222;
	}
	es = *srt;
/* Computing 2nd power */
	r__1 = es;
/* Computing 2nd power */
	r__2 = t1dlk;
/* Computing 2nd power */
	r__3 = es;
/* Computing 2nd power */
	r__4 = t2dlk;
/* Computing 2nd power */
	r__5 = es;
	pmdlk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (
		r__5 * r__5 * 4.f);
	pmdlk = sqrt(pmdlk2);
	xsk3 = pplpk_(srt) * 1.5f;
/* for NSK channel */
	t1nsk = ana + as + akp;
	t2nsk = ana + as - akp;
	if (*srt <= t1nsk) {
	    goto L222;
	}
/* Computing 2nd power */
	r__1 = es;
/* Computing 2nd power */
	r__2 = t1nsk;
/* Computing 2nd power */
	r__3 = es;
/* Computing 2nd power */
	r__4 = t2nsk;
/* Computing 2nd power */
	r__5 = es;
	pmnsk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (
		r__5 * r__5 * 4.f);
	pmnsk = sqrt(pmnsk2);
	xsk2 = (ppk1_(srt) + ppk0_(srt)) * 1.5f;
/* for DSK channel */
	t1dsk = ada + as + akp;
	t2dsk = ada + as - akp;
	if (*srt <= t1dsk) {
	    goto L222;
	}
/* Computing 2nd power */
	r__1 = es;
/* Computing 2nd power */
	r__2 = t1dsk;
/* Computing 2nd power */
	r__3 = es;
/* Computing 2nd power */
	r__4 = t2dsk;
/* Computing 2nd power */
	r__5 = es;
	pmdsk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (
		r__5 * r__5 * 4.f);
	pmdsk = sqrt(pmdsk2);
	xsk4 = (ppk1_(srt) + ppk0_(srt)) * 1.5f;
/* sp11/21/01 */
/* phi production */
	if (*srt <= 2.898914f) {
	    goto L222;
	}
/*  !! mb put the correct form */
	xsk5 = 1e-4f;
/* sp11/21/01 end */

/* THE TOTAL KAON+ PRODUCTION CROSS SECTION IS THEN */
L222:
	sigk = xsk1 + xsk2 + xsk3 + xsk4;
/* bz3/7/99 neutralk */
	xsk1 *= 2.f;
	xsk2 *= 2.f;
	xsk3 *= 2.f;
	xsk4 *= 2.f;
	sigk = sigk * 2.f + xsk5;
/* bz3/7/99 neutralk end */

/* * FOR P+P or L/S+L/S COLLISION: */
/*       lb1=lb(i1) */
/*       lb2=lb(i2) */
	leadng_1.lb1 = (i__1 = ee_1.lb[*i1 - 1], abs(i__1));
	dpi_1.lb2 = (i__1 = ee_1.lb[*i2 - 1], abs(i__1));
	if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 1 || leadng_1.lb1 <= 17 && 
		leadng_1.lb1 >= 14 && (dpi_1.lb2 <= 17 && dpi_1.lb2 >= 14) || 
		leadng_1.lb1 <= 2 && (dpi_1.lb2 <= 17 && dpi_1.lb2 >= 14) || 
		dpi_1.lb2 <= 2 && (leadng_1.lb1 <= 17 && leadng_1.lb1 >= 14)) 
		{
/* lin-8/2008 PP->d+meson here: */
	    if (x1 <= (*signn + sdprod) / *sig) {
		goto L108;
	    }
	    sig1 = sigma_(srt, &c__1, &c__1, &c__0) + sigma_(srt, &c__1, &
		    c__1, &c__1) * .5f;
	    sig2 = sigma_(srt, &c__1, &c__1, &c__1) * 1.5f;
	    signd = sig1 + sig2 + sig3 + sig4 + x1535 + sigk + s4pi + srho + 
		    somega;
/* lin-5/2008: */
/*           IF (X1.GT.(SIGNN+SIGND)/SIG)RETURN */
	    if (x1 > (*signn + signd + sdprod) / *sig) {
		return 0;
	    }
	    input_1.dir = sig3 / signd;
	    if (ranart_(&rndf77_1.nseed) <= input_1.dir) {
		goto L106;
	    }
	    if (ranart_(&rndf77_1.nseed) <= sigk / (sigk + x1535 + sig4 + 
		    sig2 + sig1 + s4pi + srho + somega)) {
		goto L306;
	    }
	    if (ranart_(&rndf77_1.nseed) <= s4pi / (x1535 + sig4 + sig2 + 
		    sig1 + s4pi + srho + somega)) {
		goto L307;
	    }
	    if (ranart_(&rndf77_1.nseed) <= srho / (x1535 + sig4 + sig2 + 
		    sig1 + srho + somega)) {
		goto L308;
	    }
	    if (ranart_(&rndf77_1.nseed) <= somega / (x1535 + sig4 + sig2 + 
		    sig1 + somega)) {
		goto L309;
	    }
	    if (ranart_(&rndf77_1.nseed) <= x1535 / (sig1 + sig2 + sig4 + 
		    x1535)) {
/* N*(1535) production */
		n12 = 9;
	    } else {
		if (ranart_(&rndf77_1.nseed) <= sig4 / (sig1 + sig2 + sig4)) {
/* DOUBLE DELTA PRODUCTION */
		    n12 = 66;
		    goto L1012;
		} else {
/* DELTA PRODUCTION */
		    n12 = 3;
		    if (ranart_(&rndf77_1.nseed) > sig1 / (sig1 + sig2)) {
			n12 = 4;
		    }
		}
	    }
	    goto L1011;
	}
/* * FOR N+N COLLISION: */
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 && (i__2 = ee_1.lb[*i2 
		- 1], abs(i__2)) == 2) {
/* lin-8/2008 NN->d+meson here: */
	    if (x1 <= (*signn + sdprod) / *sig) {
		goto L108;
	    }
	    sig1 = sigma_(srt, &c__1, &c__1, &c__0) + sigma_(srt, &c__1, &
		    c__1, &c__1) * .5f;
	    sig2 = sigma_(srt, &c__1, &c__1, &c__1) * 1.5f;
	    signd = sig1 + sig2 + x1535 + sig3 + sig4 + sigk + s4pi + srho + 
		    somega;
/* lin-5/2008: */
/*           IF (X1.GT.(SIGNN+SIGND)/SIG)RETURN */
	    if (x1 > (*signn + signd + sdprod) / *sig) {
		return 0;
	    }
	    input_1.dir = sig3 / signd;
	    if (ranart_(&rndf77_1.nseed) <= input_1.dir) {
		goto L106;
	    }
	    if (ranart_(&rndf77_1.nseed) <= sigk / (sigk + x1535 + sig4 + 
		    sig2 + sig1 + s4pi + srho + somega)) {
		goto L306;
	    }
	    if (ranart_(&rndf77_1.nseed) <= s4pi / (x1535 + sig4 + sig2 + 
		    sig1 + s4pi + srho + somega)) {
		goto L307;
	    }
	    if (ranart_(&rndf77_1.nseed) <= srho / (x1535 + sig4 + sig2 + 
		    sig1 + srho + somega)) {
		goto L308;
	    }
	    if (ranart_(&rndf77_1.nseed) <= somega / (x1535 + sig4 + sig2 + 
		    sig1 + somega)) {
		goto L309;
	    }
	    if (ranart_(&rndf77_1.nseed) <= x1535 / (x1535 + sig1 + sig2 + 
		    sig4)) {
/* N*(1535) PRODUCTION */
		n12 = 10;
	    } else {
		if (ranart_(&rndf77_1.nseed) <= sig4 / (sig1 + sig2 + sig4)) {
/* double delta production */
		    n12 = 67;
		    goto L1013;
		} else {
/* DELTA PRODUCTION */
		    n12 = 6;
		    if (ranart_(&rndf77_1.nseed) > sig1 / (sig1 + sig2)) {
			n12 = 5;
		    }
		}
	    }
	    goto L1011;
	}
/* * FOR N+P COLLISION */
	if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 2) {
/* lin-5/2008 NP->d+meson here: */
	    if (x1 <= (*signn + sdprod) / *sig) {
		goto L108;
	    }
	    sig1 = sigma_(srt, &c__1, &c__1, &c__1) * .5f + sigma_(srt, &c__1,
		     &c__1, &c__0) * .25f;
	    if (input_1.nstar == 1) {
		sig2 = sigma_(srt, &c__2, &c__0, &c__1) * .75f;
	    } else {
		sig2 = 0.f;
	    }
	    signd = (sig1 + sig2 + x1535) * 2.f + sig3 + sig4 + sigk + s4pi + 
		    srho + somega;
/* lin-5/2008: */
/*           IF (X1.GT.(SIGNN+SIGND)/SIG)RETURN */
	    if (x1 > (*signn + signd + sdprod) / *sig) {
		return 0;
	    }
	    input_1.dir = sig3 / signd;
	    if (ranart_(&rndf77_1.nseed) <= input_1.dir) {
		goto L106;
	    }
	    if (ranart_(&rndf77_1.nseed) <= sigk / (signd - sig3)) {
		goto L306;
	    }
	    if (ranart_(&rndf77_1.nseed) <= s4pi / (signd - sig3 - sigk)) {
		goto L307;
	    }
	    if (ranart_(&rndf77_1.nseed) <= srho / (signd - sig3 - sigk - 
		    s4pi)) {
		goto L308;
	    }
	    if (ranart_(&rndf77_1.nseed) <= somega / (signd - sig3 - sigk - 
		    s4pi - srho)) {
		goto L309;
	    }
	    if (ranart_(&rndf77_1.nseed) < x1535 / (sig1 + sig2 + x1535 + 
		    sig4 * .5f)) {
/* N*(1535) PRODUCTION */
		n12 = 11;
		if (ranart_(&rndf77_1.nseed) <= .5f) {
		    n12 = 12;
		}
	    } else {
		if (ranart_(&rndf77_1.nseed) <= sig4 / (sig4 + (sig1 + sig2) *
			 2.f)) {
/* double resonance production */
		    n12 = 68;
		    goto L1014;
		} else {
		    if (ranart_(&rndf77_1.nseed) <= sig1 / (sig1 + sig2)) {
/* DELTA PRODUCTION */
			n12 = 2;
			if (ranart_(&rndf77_1.nseed) >= .5f) {
			    n12 = 1;
			}
		    } else {
/* N*(1440) PRODUCTION */
			n12 = 8;
			if (ranart_(&rndf77_1.nseed) >= .5f) {
			    n12 = 7;
			}
		    }
		}
	    }
	}
L1011:
	*iblock = 2;
/* PARAMETRIZATION OF THE SHAPE OF THE DELTA RESONANCE ACCORDING */
/*     TO kitazoe's or J.D.JACKSON'S MASS FORMULA AND BREIT WIGNER */
/*     FORMULA FOR N* RESORANCE */
/*     DETERMINE DELTA MASS VIA REJECTION METHOD. */
	dmax__ = *srt - .9383f - .005f;
	dmax__ = *srt - .9383f - .005f;
	dmin__ = 1.078f;
	if (n12 < 7) {
/* Delta(1232) production */
	    if (dmax__ < 1.232f) {
		fm = fde_(&dmax__, srt, &c_b182);
	    } else {
/* lin-10/25/02 get rid of argument usage mismatch in FDE(): */
		xdmass = 1.232f;
/*          FM=FDE(1.232,SRT,1.) */
		fm = fde_(&xdmass, srt, &c_b183);
/* lin-10/25/02-end */
	    }
	    if (fm == 0.f) {
		fm = 1e-9f;
	    }
	    ntry1 = 0;
L10:
	    dm = ranart_(&rndf77_1.nseed) * (dmax__ - dmin__) + dmin__;
	    ++ntry1;
	    if (ranart_(&rndf77_1.nseed) > fde_(&dm, srt, &c_b183) / fm && 
		    ntry1 <= 30) {
		goto L10;
	    }
/* lin-2/26/03 limit the Delta mass below a certain value */
/*     (here taken as its central value + 2* B-W fullwidth): */
	    if (dm > 1.47f) {
		goto L10;
	    }
	    goto L13;
	}
	if (n12 == 7 || n12 == 8) {
/* N*(1440) production */
	    if (dmax__ < 1.44f) {
		fm = fns_(&dmax__, srt, &c_b182);
	    } else {
/* lin-10/25/02 get rid of argument usage mismatch in FNS(): */
		xdmass = 1.44f;
/*          FM=FNS(1.44,SRT,1.) */
		fm = fns_(&xdmass, srt, &c_b183);
/* lin-10/25/02-end */
	    }
	    if (fm == 0.f) {
		fm = 1e-9f;
	    }
	    ntry2 = 0;
L11:
	    dm = ranart_(&rndf77_1.nseed) * (dmax__ - dmin__) + dmin__;
	    ++ntry2;
	    if (ranart_(&rndf77_1.nseed) > fns_(&dm, srt, &c_b183) / fm && 
		    ntry2 <= 10) {
		goto L11;
	    }
/* lin-2/26/03 limit the N* mass below a certain value */
/*     (here taken as its central value + 2* B-W fullwidth): */
	    if (dm > 2.14f) {
		goto L11;
	    }
	    goto L13;
	}
	if (n12 >= 17) {
/* N*(1535) production */
	    if (dmax__ < 1.535f) {
		fm = fd5_(&dmax__, srt, &c_b182);
	    } else {
/* lin-10/25/02 get rid of argument usage mismatch in FNS(): */
		xdmass = 1.535f;
/*          FM=FD5(1.535,SRT,1.) */
		fm = fd5_(&xdmass, srt, &c_b183);
/* lin-10/25/02-end */
	    }
	    if (fm == 0.f) {
		fm = 1e-9f;
	    }
	    ntry1 = 0;
L12:
	    dm = ranart_(&rndf77_1.nseed) * (dmax__ - dmin__) + dmin__;
	    ++ntry1;
	    if (ranart_(&rndf77_1.nseed) > fd5_(&dm, srt, &c_b183) / fm && 
		    ntry1 <= 10) {
		goto L12;
	    }
/* lin-2/26/03 limit the N* mass below a certain value */
/*     (here taken as its central value + 2* B-W fullwidth): */
	    if (dm > 1.84f) {
		goto L12;
	    }
	    goto L13;
	}
/* CALCULATE THE MASSES OF BARYON RESONANCES IN THE DOUBLE RESONANCE */
/* PRODUCTION PROCESS AND RELABLE THE PARTICLES */
L1012:
	*iblock = 43;
	rmasdd_(srt, &c_b195, &c_b195, &c_b197, &c_b197, &input1_1.iseed, &
		c__1, &dm1, &dm2);
	rmasdd_(srt, &c_b195, &c_b201, &c_b197, &c_b197, &input1_1.iseed, &
		c__3, &dm1n, &dm2n);
	if (n12 == 66) {
/* (1) PP-->DOUBLE RESONANCES */
/* DETERMINE THE FINAL STATE */
	    xfinal = ranart_(&rndf77_1.nseed);
	    if (xfinal <= .25f) {
/* (1.1) D+++D0 */
		ee_1.lb[*i1 - 1] = 9;
		ee_1.lb[*i2 - 1] = 7;
		cc_1.e[*i1 - 1] = dm1;
		cc_1.e[*i2 - 1] = dm2;
		goto L200;
/* go to 200 to set the new momentum */
	    }
	    if (xfinal > .25f && xfinal <= .5f) {
/* (1.2) D++D+ */
		ee_1.lb[*i1 - 1] = 8;
		ee_1.lb[*i2 - 1] = 8;
		cc_1.e[*i1 - 1] = dm1;
		cc_1.e[*i2 - 1] = dm2;
		goto L200;
/* go to 200 to set the new momentum */
	    }
	    if (xfinal > .5f && xfinal <= .75f) {
/* (1.3) D+++N*0 */
		ee_1.lb[*i1 - 1] = 9;
		ee_1.lb[*i2 - 1] = 10;
		cc_1.e[*i1 - 1] = dm1n;
		cc_1.e[*i2 - 1] = dm2n;
		goto L200;
/* go to 200 to set the new momentum */
	    }
	    if (xfinal > .75f) {
/* (1.4) D++N*+ */
		ee_1.lb[*i1 - 1] = 8;
		ee_1.lb[*i2 - 1] = 11;
		cc_1.e[*i1 - 1] = dm1n;
		cc_1.e[*i2 - 1] = dm2n;
		goto L200;
/* go to 200 to set the new momentum */
	    }
	}
L1013:
	*iblock = 43;
	rmasdd_(srt, &c_b195, &c_b195, &c_b197, &c_b197, &input1_1.iseed, &
		c__1, &dm1, &dm2);
	rmasdd_(srt, &c_b195, &c_b201, &c_b197, &c_b197, &input1_1.iseed, &
		c__3, &dm1n, &dm2n);
	if (n12 == 67) {
/* (2) NN-->DOUBLE RESONANCES */
/* DETERMINE THE FINAL STATE */
	    xfinal = ranart_(&rndf77_1.nseed);
	    if (xfinal <= .25f) {
/* (2.1) D0+D0 */
		ee_1.lb[*i1 - 1] = 7;
		ee_1.lb[*i2 - 1] = 7;
		cc_1.e[*i1 - 1] = dm1;
		cc_1.e[*i2 - 1] = dm2;
		goto L200;
/* go to 200 to set the new momentum */
	    }
	    if (xfinal > .25f && xfinal <= .5f) {
/* (2.2) D++D+ */
		ee_1.lb[*i1 - 1] = 6;
		ee_1.lb[*i2 - 1] = 8;
		cc_1.e[*i1 - 1] = dm1;
		cc_1.e[*i2 - 1] = dm2;
		goto L200;
/* go to 200 to set the new momentum */
	    }
	    if (xfinal > .5f && xfinal <= .75f) {
/* (2.3) D0+N*0 */
		ee_1.lb[*i1 - 1] = 7;
		ee_1.lb[*i2 - 1] = 10;
		cc_1.e[*i1 - 1] = dm1n;
		cc_1.e[*i2 - 1] = dm2n;
		goto L200;
/* go to 200 to set the new momentum */
	    }
	    if (xfinal > .75f) {
/* (2.4) D++N*+ */
		ee_1.lb[*i1 - 1] = 8;
		ee_1.lb[*i2 - 1] = 11;
		cc_1.e[*i1 - 1] = dm1n;
		cc_1.e[*i2 - 1] = dm2n;
		goto L200;
/* go to 200 to set the new momentum */
	    }
	}
L1014:
	*iblock = 43;
	rmasdd_(srt, &c_b195, &c_b195, &c_b197, &c_b197, &input1_1.iseed, &
		c__1, &dm1, &dm2);
	rmasdd_(srt, &c_b195, &c_b201, &c_b197, &c_b197, &input1_1.iseed, &
		c__3, &dm1n, &dm2n);
	if (n12 == 68) {
/* (3) NP-->DOUBLE RESONANCES */
/* DETERMINE THE FINAL STATE */
	    xfinal = ranart_(&rndf77_1.nseed);
	    if (xfinal <= .25f) {
/* (3.1) D0+D+ */
		ee_1.lb[*i1 - 1] = 7;
		ee_1.lb[*i2 - 1] = 8;
		cc_1.e[*i1 - 1] = dm1;
		cc_1.e[*i2 - 1] = dm2;
		goto L200;
/* go to 200 to set the new momentum */
	    }
	    if (xfinal > .25f && xfinal <= .5f) {
/* (3.2) D+++D- */
		ee_1.lb[*i1 - 1] = 9;
		ee_1.lb[*i2 - 1] = 6;
		cc_1.e[*i1 - 1] = dm1;
		cc_1.e[*i2 - 1] = dm2;
		goto L200;
/* go to 200 to set the new momentum */
	    }
	    if (xfinal > .5f && xfinal <= .75f) {
/* (3.3) D0+N*+ */
		ee_1.lb[*i1 - 1] = 7;
		ee_1.lb[*i2 - 1] = 11;
		cc_1.e[*i1 - 1] = dm1n;
		cc_1.e[*i2 - 1] = dm2n;
		goto L200;
/* go to 200 to set the new momentum */
	    }
	    if (xfinal > .75f) {
/* (3.4) D++N*0 */
		ee_1.lb[*i1 - 1] = 8;
		ee_1.lb[*i2 - 1] = 10;
		cc_1.e[*i1 - 1] = dm1n;
		cc_1.e[*i2 - 1] = dm2n;
		goto L200;
/* go to 200 to set the new momentum */
	    }
	}
L13:
/* ------------------------------------------------------- */
/* RELABLE BARYON I1 AND I2 */
/* 1. p+n-->delta(+)+n */
	if (n12 == 1) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1) {
		ee_1.lb[*i2 - 1] = 2;
		ee_1.lb[*i1 - 1] = 8;
		cc_1.e[*i1 - 1] = dm;
	    } else {
		ee_1.lb[*i1 - 1] = 2;
		ee_1.lb[*i2 - 1] = 8;
		cc_1.e[*i2 - 1] = dm;
	    }
	    goto L200;
	}
/* 2 p+n-->delta(0)+p */
	if (n12 == 2) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2) {
		ee_1.lb[*i2 - 1] = 1;
		ee_1.lb[*i1 - 1] = 7;
		cc_1.e[*i1 - 1] = dm;
	    } else {
		ee_1.lb[*i1 - 1] = 1;
		ee_1.lb[*i2 - 1] = 7;
		cc_1.e[*i2 - 1] = dm;
	    }
	    goto L200;
	}
/* 3 p+p-->delta(++)+n */
	if (n12 == 3) {
	    ee_1.lb[*i1 - 1] = 9;
	    cc_1.e[*i1 - 1] = dm;
	    ee_1.lb[*i2 - 1] = 2;
	    cc_1.e[*i2 - 1] = .939457f;
	    goto L200;
	}
/* 4 p+p-->delta(+)+p */
	if (n12 == 4) {
	    ee_1.lb[*i2 - 1] = 1;
	    ee_1.lb[*i1 - 1] = 8;
	    cc_1.e[*i1 - 1] = dm;
	    goto L200;
	}
/* 5 n+n--> delta(0)+n */
	if (n12 == 5) {
	    ee_1.lb[*i2 - 1] = 2;
	    ee_1.lb[*i1 - 1] = 7;
	    cc_1.e[*i1 - 1] = dm;
	    goto L200;
	}
/* 6 n+n--> delta(-)+p */
	if (n12 == 6) {
	    ee_1.lb[*i1 - 1] = 6;
	    cc_1.e[*i1 - 1] = dm;
	    ee_1.lb[*i2 - 1] = 1;
	    cc_1.e[*i2 - 1] = .93828f;
	    goto L200;
	}
/* 7 n+p--> N*(0)+p */
	if (n12 == 7) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1) {
		ee_1.lb[*i1 - 1] = 1;
		ee_1.lb[*i2 - 1] = 10;
		cc_1.e[*i2 - 1] = dm;
	    } else {
		ee_1.lb[*i2 - 1] = 1;
		ee_1.lb[*i1 - 1] = 10;
		cc_1.e[*i1 - 1] = dm;
	    }
	    goto L200;
	}
/* 8 n+p--> N*(+)+n */
	if (n12 == 8) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1) {
		ee_1.lb[*i2 - 1] = 2;
		ee_1.lb[*i1 - 1] = 11;
		cc_1.e[*i1 - 1] = dm;
	    } else {
		ee_1.lb[*i1 - 1] = 2;
		ee_1.lb[*i2 - 1] = 11;
		cc_1.e[*i2 - 1] = dm;
	    }
	    goto L200;
	}
/* 9 p+p--> N*(+)(1535)+p */
	if (n12 == 9) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 1;
		ee_1.lb[*i1 - 1] = 13;
		cc_1.e[*i1 - 1] = dm;
	    } else {
		ee_1.lb[*i1 - 1] = 1;
		ee_1.lb[*i2 - 1] = 13;
		cc_1.e[*i2 - 1] = dm;
	    }
	    goto L200;
	}
/* 10 n+n--> N*(0)(1535)+n */
	if (n12 == 10) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 2;
		ee_1.lb[*i1 - 1] = 12;
		cc_1.e[*i1 - 1] = dm;
	    } else {
		ee_1.lb[*i1 - 1] = 2;
		ee_1.lb[*i2 - 1] = 12;
		cc_1.e[*i2 - 1] = dm;
	    }
	    goto L200;
	}
/* 11 n+p--> N*(+)(1535)+n */
	if (n12 == 11) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2) {
		ee_1.lb[*i1 - 1] = 2;
		ee_1.lb[*i2 - 1] = 13;
		cc_1.e[*i2 - 1] = dm;
	    } else {
		ee_1.lb[*i2 - 1] = 2;
		ee_1.lb[*i1 - 1] = 13;
		cc_1.e[*i1 - 1] = dm;
	    }
	    goto L200;
	}
/* 12 n+p--> N*(0)(1535)+p */
	if (n12 == 12) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1) {
		ee_1.lb[*i1 - 1] = 1;
		ee_1.lb[*i2 - 1] = 12;
		cc_1.e[*i2 - 1] = dm;
	    } else {
		ee_1.lb[*i2 - 1] = 1;
		ee_1.lb[*i1 - 1] = 12;
		cc_1.e[*i1 - 1] = dm;
	    }
	}
    }
/* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH */
/* ENERGY CONSERVATION */
L200:
    leadng_1.em1 = cc_1.e[*i1 - 1];
    dpi_1.em2 = cc_1.e[*i2 - 1];
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = leadng_1.em1;
/* Computing 2nd power */
    r__4 = dpi_1.em2;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = leadng_1.em1 * dpi_1.em2;
    pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
    if (pr2 <= 0.f) {
	pr2 = 1e-9f;
    }
    pr = sqrt(pr2) / (*srt * 2.f);
    if (*srt <= 2.14f) {
	c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
    }
    if (*srt > 2.14f && *srt <= 2.4f) {
	c1 = ang_(srt, &input1_1.iseed);
    }
    if (*srt > 2.4f) {
/* lin-10/25/02 get rid of argument usage mismatch in PTR(): */
	xptr = pr * .33f;
/*         cc1=ptr(0.33*pr,iseed) */
	cc1 = ptr_(&xptr, &input1_1.iseed);
/* lin-10/25/02-end */
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
	r__1 = pr;
/* Computing 2nd power */
	r__2 = cc1;
	scheck = r__1 * r__1 - r__2 * r__2;
	if (scheck < 0.f) {
	    s_wsle(&io___560);
	    do_lio(&c__9, &c__1, "scheck2: ", (ftnlen)9);
	    do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	    e_wsle();
	    scheck = 0.f;
	}
	c1 = sqrt(scheck) / pr;
/*             c1=sqrt(pr**2-cc1**2)/pr */
    }
    t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
    if (ianti == 1 && ee_1.lb[*i1 - 1] >= 1 && ee_1.lb[*i2 - 1] >= 1) {
	ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
	ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
    }
    goto L107;
/* FOR THE NN-->D1+D2+PI PROCESS, FIND MOMENTUM OF THE FINAL TWO */
/* DELTAS AND PION IN THE NUCLEUS-NUCLEUS CMS. */
L106:
    ntry1 = 0;
L123:
    ddp2_(srt, &input1_1.iseed, &px3, &py3, &pz3, &dm3, &px4, &py4, &pz4, &
	    dm4, &ppx, &ppy, &ppz, &icou1);
    ++ntry1;
    if (icou1 < 0 && ntry1 <= 40) {
	goto L123;
    }
/*       if(icou1.lt.0)return */
/* ROTATE THE MOMENTA OF PARTICLES IN THE CMS OF P1+P2 */
    rotate_(px, py, pz, &px3, &py3, &pz3);
    rotate_(px, py, pz, &px4, &py4, &pz4);
    rotate_(px, py, pz, &ppx, &ppy, &ppz);
    ++nn_1.nnn;
/* DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE */
/* (1) FOR P+P */
    xdir = ranart_(&rndf77_1.nseed);
    if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 1) {
	if (xdir <= .2f) {
/* (1.1)P+P-->D+++D0+PION(0) */
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 4;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13496f;
	    ee_1.lb[*i1 - 1] = 9;
	    ee_1.lb[*i2 - 1] = 7;
	    goto L205;
	}
/* (1.2)P+P -->D++D+PION(0) */
	if (xdir <= .4f && xdir > .2f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 4;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13496f;
	    ee_1.lb[*i1 - 1] = 8;
	    ee_1.lb[*i2 - 1] = 8;
	    goto L205;
	}
/* (1.3)P+P-->D+++D+PION(-) */
	if (xdir <= .6f && xdir > .4f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 3;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13957f;
	    ee_1.lb[*i1 - 1] = 9;
	    ee_1.lb[*i2 - 1] = 8;
	    goto L205;
	}
	if (xdir <= .8f && xdir > .6f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 5;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13957f;
	    ee_1.lb[*i1 - 1] = 9;
	    ee_1.lb[*i2 - 1] = 6;
	    goto L205;
	}
	if (xdir > .8f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 5;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13957f;
	    ee_1.lb[*i1 - 1] = 7;
	    ee_1.lb[*i2 - 1] = 8;
	    goto L205;
	}
    }
/* (2)FOR N+N */
    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 && (i__2 = ee_1.lb[*i2 - 1],
	     abs(i__2)) == 2) {
	if (xdir <= .2f) {
/* (2.1)N+N-->D++D-+PION(0) */
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 4;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13496f;
	    ee_1.lb[*i1 - 1] = 6;
	    ee_1.lb[*i2 - 1] = 7;
	    goto L205;
	}
/* (2.2)N+N -->D+++D-+PION(-) */
	if (xdir <= .4f && xdir > .2f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 3;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13957f;
	    ee_1.lb[*i1 - 1] = 6;
	    ee_1.lb[*i2 - 1] = 9;
	    goto L205;
	}
/* (2.3)P+P-->D0+D-+PION(+) */
	if (xdir > .4f && xdir <= .6f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 5;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13957f;
	    ee_1.lb[*i1 - 1] = 9;
	    ee_1.lb[*i2 - 1] = 8;
	    goto L205;
	}
/* (2.4)P+P-->D0+D0+PION(0) */
	if (xdir > .6f && xdir <= .8f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 4;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13496f;
	    ee_1.lb[*i1 - 1] = 7;
	    ee_1.lb[*i2 - 1] = 7;
	    goto L205;
	}
/* (2.5)P+P-->D0+D++PION(-) */
	if (xdir > .8f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 3;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13957f;
	    ee_1.lb[*i1 - 1] = 7;
	    ee_1.lb[*i2 - 1] = 8;
	    goto L205;
	}
    }
/* (3)FOR N+P */
    if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 2) {
	if (xdir <= .17f) {
/* (3.1)N+P-->D+++D-+PION(0) */
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 4;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13496f;
	    ee_1.lb[*i1 - 1] = 6;
	    ee_1.lb[*i2 - 1] = 9;
	    goto L205;
	}
/* (3.2)N+P -->D+++D0+PION(-) */
	if (xdir <= .34f && xdir > .17f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 3;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13957f;
	    ee_1.lb[*i1 - 1] = 7;
	    ee_1.lb[*i2 - 1] = 9;
	    goto L205;
	}
/* (3.3)N+P-->D++D-+PION(+) */
	if (xdir > .34f && xdir <= .51f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 5;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13957f;
	    ee_1.lb[*i1 - 1] = 7;
	    ee_1.lb[*i2 - 1] = 8;
	    goto L205;
	}
/* (3.4)N+P-->D++D++PION(-) */
	if (xdir > .51f && xdir <= .68f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 3;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13957f;
	    ee_1.lb[*i1 - 1] = 8;
	    ee_1.lb[*i2 - 1] = 8;
	    goto L205;
	}
/* (3.5)N+P-->D0+D++PION(0) */
	if (xdir > .68f && xdir <= .85f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 4;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13957f;
	    ee_1.lb[*i1 - 1] = 7;
	    ee_1.lb[*i2 - 1] = 8;
	    goto L205;
	}
/* (3.6)N+P-->D0+D0+PION(+) */
	if (xdir > .85f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 5;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13957f;
	    ee_1.lb[*i1 - 1] = 7;
	    ee_1.lb[*i2 - 1] = 7;
	}
    }
/* FIND THE MOMENTUM OF PARTICLES IN THE FINAL STATE IN THE NUCLEUS- */
/* NUCLEUS CMS. FRAME */
/*             LORENTZ-TRANSFORMATION INTO LAB FRAME FOR DELTA1 */
L205:
/* Computing 2nd power */
    r__1 = dm3;
/* Computing 2nd power */
    r__2 = px3;
/* Computing 2nd power */
    r__3 = py3;
/* Computing 2nd power */
    r__4 = pz3;
    e1cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    p1beta = px3 * bg_1.betax + py3 * bg_1.betay + pz3 * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) + e1cm);
    pt1i1 = bg_1.betax * transf + px3;
    pt2i1 = bg_1.betay * transf + py3;
    pt3i1 = bg_1.betaz * transf + pz3;
    eti1 = dm3;

    if (ianti == 1 && ee_1.lb[*i1 - 1] >= 1 && ee_1.lb[*i2 - 1] >= 1) {
	ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
	ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
	if (pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] == 3) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 5;
	} else if (pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] == 5) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 3;
	}
    }

    leadng_1.lb1 = ee_1.lb[*i1 - 1];
/* FOR DELTA2 */
/* Computing 2nd power */
    r__1 = dm4;
/* Computing 2nd power */
    r__2 = px4;
/* Computing 2nd power */
    r__3 = py4;
/* Computing 2nd power */
    r__4 = pz4;
    e2cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    p2beta = px4 * bg_1.betax + py4 * bg_1.betay + pz4 * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * p2beta / (bg_1.gamma + 1.f) + e2cm);
    pt1i2 = bg_1.betax * transf + px4;
    pt2i2 = bg_1.betay * transf + py4;
    pt3i2 = bg_1.betaz * transf + pz4;
    eti2 = dm4;
    dpi_1.lb2 = ee_1.lb[*i2 - 1];
/* assign delta1 and delta2 to i1 or i2 to keep the leadng particle */
/* behaviour */
/*              if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0)then */
    bb_1.p[*i1 * 3 - 3] = pt1i1;
    bb_1.p[*i1 * 3 - 2] = pt2i1;
    bb_1.p[*i1 * 3 - 1] = pt3i1;
    cc_1.e[*i1 - 1] = eti1;
    ee_1.lb[*i1 - 1] = leadng_1.lb1;
    bb_1.p[*i2 * 3 - 3] = pt1i2;
    bb_1.p[*i2 * 3 - 2] = pt2i2;
    bb_1.p[*i2 * 3 - 1] = pt3i2;
    cc_1.e[*i2 - 1] = eti2;
    ee_1.lb[*i2 - 1] = dpi_1.lb2;
    leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
    leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
    leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
    leadng_1.em1 = cc_1.e[*i1 - 1];
    ee_1.id[*i1 - 1] = 2;
    ee_1.id[*i2 - 1] = 2;
    id1 = ee_1.id[*i1 - 1];
    *iblock = 4;
/* GET PION'S MOMENTUM AND COORDINATES IN NUCLEUS-NUCLEUS CMS. FRAME */
/* Computing 2nd power */
    r__1 = pc_1.epion[nn_1.nnn + *irun * 150001 - 150002];
/* Computing 2nd power */
    r__2 = ppx;
/* Computing 2nd power */
    r__3 = ppy;
/* Computing 2nd power */
    r__4 = ppz;
    epcm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    ppbeta = ppx * bg_1.betax + ppy * bg_1.betay + ppz * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * ppbeta / (bg_1.gamma + 1.f) + epcm);
    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = bg_1.betax * 
	    transf + ppx;
    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = bg_1.betay * 
	    transf + ppy;
    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = bg_1.betaz * 
	    transf + ppz;
/* lin-5/2008: */
    dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = dpert_1.dpertp[*i1 - 
	    1] * dpert_1.dpertp[*i2 - 1];
/* lin-5/2008 do not allow smearing in position of produced particles */
/*     to avoid immediate reinteraction with the particle I1, I2 or themselves: */
/* 2002        X01 = 1.0 - 2.0 * RANART(NSEED) */
/*            Y01 = 1.0 - 2.0 * RANART(NSEED) */
/*            Z01 = 1.0 - 2.0 * RANART(NSEED) */
/*        IF ((X01*X01+Y01*Y01+Z01*Z01) .GT. 1.0) GOTO 2002 */
/*                RPION(1,NNN,IRUN)=R(1,I1)+0.5*x01 */
/*                RPION(2,NNN,IRUN)=R(2,I1)+0.5*y01 */
/*                RPION(3,NNN,IRUN)=R(3,I1)+0.5*z01 */
    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = aa_1.r__[*i1 * 3 - 
	    3];
    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = aa_1.r__[*i1 * 3 - 
	    2];
    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = aa_1.r__[*i1 * 3 - 
	    1];

    goto L90005;
/* lin-5/2008 N+N->Deuteron+pi: */
/*     FIND MOMENTUM OF THE FINAL PARTICLES IN THE NUCLEUS-NUCLEUS CMS. */
L108:
    if (para8_1.idpert == 1 && *ipert1 == 1 && para8_1.npertd >= 1) {
/*     For idpert=1: we produce npertd pert deuterons: */
	ndloop = para8_1.npertd;
    } else if (para8_1.idpert == 2 && para8_1.npertd >= 1) {
/*     For idpert=2: we first save information for npertd pert deuterons; */
/*     at the last ndloop we create the regular deuteron+pi */
/*     and those pert deuterons: */
	ndloop = para8_1.npertd + 1;
    } else {
/*     Just create the regular deuteron+pi: */
	ndloop = 1;
    }

    dprob1 = sdprod / *sig / (real) para8_1.npertd;
    i__1 = ndloop;
    for (idloop = 1; idloop <= i__1; ++idloop) {
	bbdangle_(&pxd, &pyd, &pzd, nt, ipert1, &ianti, &idloop, &pfinal, &
		dprob1, &lbm);
	rotate_(px, py, pz, &pxd, &pyd, &pzd);
/*     LORENTZ-TRANSFORMATION OF THE MOMENTUM OF PARTICLES IN THE FINAL STATE */
/*     FROM THE NN CMS FRAME INTO THE GLOBAL CMS FRAME: */
/*     For the Deuteron: */
	xmass = 1.8756f;
/* Computing 2nd power */
	r__1 = xmass;
/* Computing 2nd power */
	r__2 = pxd;
/* Computing 2nd power */
	r__3 = pyd;
/* Computing 2nd power */
	r__4 = pzd;
	e1dcm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p1dbeta = pxd * bg_1.betax + pyd * bg_1.betay + pzd * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p1dbeta / (bg_1.gamma + 1.f) + 
		e1dcm);
	pxi1 = bg_1.betax * transf + pxd;
	pyi1 = bg_1.betay * transf + pyd;
	pzi1 = bg_1.betaz * transf + pzd;
	if (ianti == 0) {
	    lbd = 42;
	} else {
	    lbd = -42;
	}
	if (para8_1.idpert == 1 && *ipert1 == 1 && para8_1.npertd >= 1) {
/* ccc  Perturbative production for idpert=1: */
	    ++nn_1.nnn;
	    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = pxi1;
	    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = pyi1;
	    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = pzi1;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = 1.8756f;
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = lbd;
	    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = aa_1.r__[*
		    i1 * 3 - 3];
	    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = aa_1.r__[*
		    i1 * 3 - 2];
	    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = aa_1.r__[*
		    i1 * 3 - 1];
/* lin-5/2008 assign the perturbative probability: */
	    dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = sdprod / *
		    sig / (real) para8_1.npertd;
	} else if (para8_1.idpert == 2 && idloop <= para8_1.npertd) {
/* lin-5/2008 For idpert=2, we produce NPERTD perturbative (anti)deuterons */
/*     only when a regular (anti)deuteron+pi is produced in NN collisions. */
/*     First save the info for the perturbative deuterons: */
	    ppd[idloop * 3 - 3] = pxi1;
	    ppd[idloop * 3 - 2] = pyi1;
	    ppd[idloop * 3 - 1] = pzi1;
	    lbpd[idloop - 1] = lbd;
	} else {
/* ccc  Regular production: */
/*     For the regular pion: do LORENTZ-TRANSFORMATION: */
	    cc_1.e[*i1 - 1] = xmm;
/* Computing 2nd power */
	    r__1 = xmm;
/* Computing 2nd power */
	    r__2 = pxd;
/* Computing 2nd power */
	    r__3 = pyd;
/* Computing 2nd power */
	    r__4 = pzd;
	    e2picm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * 
		    r__4);
	    p2pibeta = -pxd * bg_1.betax - pyd * bg_1.betay - pzd * 
		    bg_1.betaz;
	    transf = bg_1.gamma * (bg_1.gamma * p2pibeta / (bg_1.gamma + 1.f) 
		    + e2picm);
	    pxi2 = bg_1.betax * transf - pxd;
	    pyi2 = bg_1.betay * transf - pyd;
	    pzi2 = bg_1.betaz * transf - pzd;
	    bb_1.p[*i1 * 3 - 3] = pxi2;
	    bb_1.p[*i1 * 3 - 2] = pyi2;
	    bb_1.p[*i1 * 3 - 1] = pzi2;
/*     Remove regular pion to check the equivalence */
/*     between the perturbative and regular deuteron results: */
/*                 E(i1)=0. */

	    ee_1.lb[*i1 - 1] = lbm;
	    leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
	    leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
	    leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
	    leadng_1.em1 = cc_1.e[*i1 - 1];
	    ee_1.id[*i1 - 1] = 2;
	    id1 = ee_1.id[*i1 - 1];
/* Computing 2nd power */
	    r__1 = leadng_1.em1;
/* Computing 2nd power */
	    r__2 = leadng_1.px1;
/* Computing 2nd power */
	    r__3 = leadng_1.py1;
/* Computing 2nd power */
	    r__4 = leadng_1.pz1;
	    leadng_1.e1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 
		    * r__4);
	    leadng_1.lb1 = ee_1.lb[*i1 - 1];
/*     For the regular deuteron: */
	    bb_1.p[*i2 * 3 - 3] = pxi1;
	    bb_1.p[*i2 * 3 - 2] = pyi1;
	    bb_1.p[*i2 * 3 - 1] = pzi1;
	    ee_1.lb[*i2 - 1] = lbd;
	    dpi_1.lb2 = ee_1.lb[*i2 - 1];
	    cc_1.e[*i2 - 1] = 1.8756f;
	    eti2 = cc_1.e[*i2 - 1];
	    ee_1.id[*i2 - 1] = 2;
/*     For idpert=2: create the perturbative deuterons: */
	    if (para8_1.idpert == 2 && idloop == ndloop) {
		i__2 = para8_1.npertd;
		for (ipertd = 1; ipertd <= i__2; ++ipertd) {
		    ++nn_1.nnn;
		    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = 
			    ppd[ipertd * 3 - 3];
		    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = 
			    ppd[ipertd * 3 - 2];
		    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = 
			    ppd[ipertd * 3 - 1];
		    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = 1.8756f;
		    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = lbpd[
			    ipertd - 1];
		    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = 
			    aa_1.r__[*i1 * 3 - 3];
		    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = 
			    aa_1.r__[*i1 * 3 - 2];
		    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = 
			    aa_1.r__[*i1 * 3 - 1];
/* lin-5/2008 assign the perturbative probability: */
		    dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = 1.f /
			     (real) para8_1.npertd;
		}
	    }
	}
    }
    *iblock = 501;
    goto L90005;
/* lin-5/2008 N+N->Deuteron+pi over */
/* FOR THE NN-->KAON+X PROCESS, FIND MOMENTUM OF THE FINAL PARTICLES IN */
/* THE NUCLEUS-NUCLEUS CMS. */
L306:
/* sp11/21/01 phi production */
    if (xsk5 / sigk > ranart_(&rndf77_1.nseed)) {
	leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
	pz2 = bb_1.p[*i2 * 3 - 1];
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	++nn_1.nnn;
	pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 29;
	pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = 1.02f;
	*iblock = 222;
	goto L208;
    }

    *iblock = 9;
    if (ianti == 1) {
	*iblock = -9;
    }

    leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
    pz2 = bb_1.p[*i2 * 3 - 1];
/* DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE */
    ++nn_1.nnn;
    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 23;
    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .498f;
    if (*srt <= 2.63f) {
/* only lambda production is possible */
/* (1.1)P+P-->p+L+kaon+ */
	ic = 1;
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	ee_1.lb[*i2 - 1] = 14;
	goto L208;
    }
    if (*srt <= 2.74f && *srt > 2.63f) {
/* both Lambda and sigma production are possible */
	if (xsk1 / (xsk1 + xsk2) > ranart_(&rndf77_1.nseed)) {
/* lambda production */
	    ic = 1;
	    ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	    ee_1.lb[*i2 - 1] = 14;
	} else {
/* sigma production */
	    ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	    ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 15;
	    ic = 2;
	}
	goto L208;
    }
    if (*srt <= 2.77f && *srt > 2.74f) {
/* then pp-->Delta lamda kaon can happen */
	if (xsk1 / (xsk1 + xsk2 + xsk3) > ranart_(&rndf77_1.nseed)) {
/* * (1.1)P+P-->p+L+kaon+ */
	    ic = 1;
	    ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	    ee_1.lb[*i2 - 1] = 14;
	    goto L208;
	} else {
	    if (xsk2 / (xsk2 + xsk3) > ranart_(&rndf77_1.nseed)) {
/* pp-->psk */
		ic = 2;
		ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 
			1;
		ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 
			15;
	    } else {
/* pp-->D+l+k */
		ic = 3;
		ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 4) + 
			6;
		ee_1.lb[*i2 - 1] = 14;
	    }
	    goto L208;
	}
    }
    if (*srt > 2.77f) {
/* all four channels are possible */
	if (xsk1 / (xsk1 + xsk2 + xsk3 + xsk4) > ranart_(&rndf77_1.nseed)) {
/* p lambda k production */
	    ic = 1;
	    ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	    ee_1.lb[*i2 - 1] = 14;
	    goto L208;
	} else {
	    if (xsk3 / (xsk2 + xsk3 + xsk4) > ranart_(&rndf77_1.nseed)) {
/* delta l K production */
		ic = 3;
		ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 4) + 
			6;
		ee_1.lb[*i2 - 1] = 14;
		goto L208;
	    } else {
		if (xsk2 / (xsk2 + xsk4) > ranart_(&rndf77_1.nseed)) {
/* n sigma k production */
		    ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 
			    2) + 1;
		    ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 
			    3) + 15;
		    ic = 2;
		} else {
		    ic = 4;
		    ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 
			    4) + 6;
		    ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 
			    3) + 15;
		}
		goto L208;
	    }
	}
    }
L208:
    if (ianti == 1 && ee_1.lb[*i1 - 1] >= 1 && ee_1.lb[*i2 - 1] >= 1) {
	ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
	ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
	if (pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] == 23) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 21;
	}
    }
/* KEEP ALL COORDINATES OF PARTICLE 2 FOR POSSIBLE PHASE SPACE CHANGE */
    ntry1 = 0;
L127:
    bbkaon_(&ic, srt, &px3, &py3, &pz3, &dm3, &px4, &py4, &pz4, &dm4, &ppx, &
	    ppy, &ppz, &icou1);
    ++ntry1;
    if (icou1 < 0 && ntry1 <= 20) {
	goto L127;
    }
/*       if(icou1.lt.0)return */
/* ROTATE THE MOMENTA OF PARTICLES IN THE CMS OF P1+P2 */
    rotate_(px, py, pz, &px3, &py3, &pz3);
    rotate_(px, py, pz, &px4, &py4, &pz4);
    rotate_(px, py, pz, &ppx, &ppy, &ppz);
/* FIND THE MOMENTUM OF PARTICLES IN THE FINAL STATE IN THE NUCLEUS- */
/* NUCLEUS CMS. FRAME */
/* (1) for the necleon/delta */
/*             LORENTZ-TRANSFORMATION INTO LAB FRAME FOR DELTA1 */
/* Computing 2nd power */
    r__1 = dm3;
/* Computing 2nd power */
    r__2 = px3;
/* Computing 2nd power */
    r__3 = py3;
/* Computing 2nd power */
    r__4 = pz3;
    e1cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    p1beta = px3 * bg_1.betax + py3 * bg_1.betay + pz3 * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) + e1cm);
    pt1i1 = bg_1.betax * transf + px3;
    pt2i1 = bg_1.betay * transf + py3;
    pt3i1 = bg_1.betaz * transf + pz3;
    eti1 = dm3;
    lbi1 = ee_1.lb[*i1 - 1];
/* (2) for the lambda/sigma */
/* Computing 2nd power */
    r__1 = dm4;
/* Computing 2nd power */
    r__2 = px4;
/* Computing 2nd power */
    r__3 = py4;
/* Computing 2nd power */
    r__4 = pz4;
    e2cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    p2beta = px4 * bg_1.betax + py4 * bg_1.betay + pz4 * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * p2beta / (bg_1.gamma + 1.f) + e2cm);
    pt1i2 = bg_1.betax * transf + px4;
    pt2i2 = bg_1.betay * transf + py4;
    pt3i2 = bg_1.betaz * transf + pz4;
    eti2 = dm4;
    lbi2 = ee_1.lb[*i2 - 1];
/* GET the kaon'S MOMENTUM AND COORDINATES IN NUCLEUS-NUCLEUS CMS. FRAME */
/* Computing 2nd power */
    r__1 = ppx;
/* Computing 2nd power */
    r__2 = ppy;
/* Computing 2nd power */
    r__3 = ppz;
    epcm = sqrt(r__1 * r__1 + .248004f + r__2 * r__2 + r__3 * r__3);
    ppbeta = ppx * bg_1.betax + ppy * bg_1.betay + ppz * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * ppbeta / (bg_1.gamma + 1.f) + epcm);
    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = bg_1.betax * 
	    transf + ppx;
    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = bg_1.betay * 
	    transf + ppy;
    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = bg_1.betaz * 
	    transf + ppz;
/* lin-5/2008 */
    dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = dpert_1.dpertp[*i1 - 
	    1] * dpert_1.dpertp[*i2 - 1];
/* lin-5/2008 */
/* 2003        X01 = 1.0 - 2.0 * RANART(NSEED) */
/*            Y01 = 1.0 - 2.0 * RANART(NSEED) */
/*            Z01 = 1.0 - 2.0 * RANART(NSEED) */
/*        IF ((X01*X01+Y01*Y01+Z01*Z01) .GT. 1.0) GOTO 2003 */
/*                RPION(1,NNN,IRUN)=R(1,I1)+0.5*x01 */
/*                RPION(2,NNN,IRUN)=R(2,I1)+0.5*y01 */
/*                RPION(3,NNN,IRUN)=R(3,I1)+0.5*z01 */
    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = aa_1.r__[*i1 * 3 - 
	    3];
    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = aa_1.r__[*i1 * 3 - 
	    2];
    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = aa_1.r__[*i1 * 3 - 
	    1];

/* assign the nucleon/delta and lambda/sigma to i1 or i2 to keep the */
/* leadng particle behaviour */
/*              if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0)then */
    bb_1.p[*i1 * 3 - 3] = pt1i1;
    bb_1.p[*i1 * 3 - 2] = pt2i1;
    bb_1.p[*i1 * 3 - 1] = pt3i1;
    cc_1.e[*i1 - 1] = eti1;
    ee_1.lb[*i1 - 1] = lbi1;
    bb_1.p[*i2 * 3 - 3] = pt1i2;
    bb_1.p[*i2 * 3 - 2] = pt2i2;
    bb_1.p[*i2 * 3 - 1] = pt3i2;
    cc_1.e[*i2 - 1] = eti2;
    ee_1.lb[*i2 - 1] = lbi2;
    leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
    leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
    leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
    leadng_1.em1 = cc_1.e[*i1 - 1];
    ee_1.id[*i1 - 1] = 2;
    ee_1.id[*i2 - 1] = 2;
    id1 = ee_1.id[*i1 - 1];
    goto L90005;
/* FOR THE NN-->Delta+Delta+rho PROCESS, FIND MOMENTUM OF THE FINAL */
/* PARTICLES IN THE NUCLEUS-NUCLEUS CMS. */
L307:
    ntry1 = 0;
L125:
    ddrho_(srt, &input1_1.iseed, &px3, &py3, &pz3, &dm3, &px4, &py4, &pz4, &
	    dm4, &ppx, &ppy, &ppz, &amrho, &icou1);
    ++ntry1;
    if (icou1 < 0 && ntry1 <= 20) {
	goto L125;
    }
/*       if(icou1.lt.0)return */
/* ROTATE THE MOMENTA OF PARTICLES IN THE CMS OF P1+P2 */
    rotate_(px, py, pz, &px3, &py3, &pz3);
    rotate_(px, py, pz, &px4, &py4, &pz4);
    rotate_(px, py, pz, &ppx, &ppy, &ppz);
    ++nn_1.nnn;
    arho = amrho;
/* DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE */
/* (1) FOR P+P */
    xdir = ranart_(&rndf77_1.nseed);
    if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 1) {
	if (xdir <= .2f) {
/* (1.1)P+P-->D+++D0+rho(0) */
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 26;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
	    ee_1.lb[*i1 - 1] = 9;
	    ee_1.lb[*i2 - 1] = 7;
	    goto L2051;
	}
/* (1.2)P+P -->D++D+rho(0) */
	if (xdir <= .4f && xdir > .2f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 26;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
	    ee_1.lb[*i1 - 1] = 8;
	    ee_1.lb[*i2 - 1] = 8;
	    goto L2051;
	}
/* (1.3)P+P-->D+++D+arho(-) */
	if (xdir <= .6f && xdir > .4f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 25;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
	    ee_1.lb[*i1 - 1] = 9;
	    ee_1.lb[*i2 - 1] = 8;
	    goto L2051;
	}
	if (xdir <= .8f && xdir > .6f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 27;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
	    ee_1.lb[*i1 - 1] = 9;
	    ee_1.lb[*i2 - 1] = 6;
	    goto L2051;
	}
	if (xdir > .8f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 27;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
	    ee_1.lb[*i1 - 1] = 7;
	    ee_1.lb[*i2 - 1] = 8;
	    goto L2051;
	}
    }
/* (2)FOR N+N */
    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 && (i__2 = ee_1.lb[*i2 - 1],
	     abs(i__2)) == 2) {
	if (xdir <= .2f) {
/* (2.1)N+N-->D++D-+rho(0) */
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 26;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
	    ee_1.lb[*i1 - 1] = 6;
	    ee_1.lb[*i2 - 1] = 7;
	    goto L2051;
	}
/* (2.2)N+N -->D+++D-+rho(-) */
	if (xdir <= .4f && xdir > .2f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 25;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
	    ee_1.lb[*i1 - 1] = 6;
	    ee_1.lb[*i2 - 1] = 9;
	    goto L2051;
	}
/* (2.3)P+P-->D0+D-+rho(+) */
	if (xdir > .4f && xdir <= .6f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 27;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
	    ee_1.lb[*i1 - 1] = 9;
	    ee_1.lb[*i2 - 1] = 8;
	    goto L2051;
	}
/* (2.4)P+P-->D0+D0+rho(0) */
	if (xdir > .6f && xdir <= .8f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 26;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
	    ee_1.lb[*i1 - 1] = 7;
	    ee_1.lb[*i2 - 1] = 7;
	    goto L2051;
	}
/* (2.5)P+P-->D0+D++rho(-) */
	if (xdir > .8f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 25;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
	    ee_1.lb[*i1 - 1] = 7;
	    ee_1.lb[*i2 - 1] = 8;
	    goto L2051;
	}
    }
/* (3)FOR N+P */
    if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 2) {
	if (xdir <= .17f) {
/* (3.1)N+P-->D+++D-+rho(0) */
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 25;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
	    ee_1.lb[*i1 - 1] = 6;
	    ee_1.lb[*i2 - 1] = 9;
	    goto L2051;
	}
/* (3.2)N+P -->D+++D0+rho(-) */
	if (xdir <= .34f && xdir > .17f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 25;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
	    ee_1.lb[*i1 - 1] = 7;
	    ee_1.lb[*i2 - 1] = 9;
	    goto L2051;
	}
/* (3.3)N+P-->D++D-+rho(+) */
	if (xdir > .34f && xdir <= .51f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 27;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
	    ee_1.lb[*i1 - 1] = 7;
	    ee_1.lb[*i2 - 1] = 8;
	    goto L2051;
	}
/* (3.4)N+P-->D++D++rho(-) */
	if (xdir > .51f && xdir <= .68f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 25;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
	    ee_1.lb[*i1 - 1] = 8;
	    ee_1.lb[*i2 - 1] = 8;
	    goto L2051;
	}
/* (3.5)N+P-->D0+D++rho(0) */
	if (xdir > .68f && xdir <= .85f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 26;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
	    ee_1.lb[*i1 - 1] = 7;
	    ee_1.lb[*i2 - 1] = 8;
	    goto L2051;
	}
/* (3.6)N+P-->D0+D0+rho(+) */
	if (xdir > .85f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 27;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
	    ee_1.lb[*i1 - 1] = 7;
	    ee_1.lb[*i2 - 1] = 7;
	}
    }
/* FIND THE MOMENTUM OF PARTICLES IN THE FINAL STATE IN THE NUCLEUS- */
/* NUCLEUS CMS. FRAME */
/*             LORENTZ-TRANSFORMATION INTO LAB FRAME FOR DELTA1 */
L2051:
/* Computing 2nd power */
    r__1 = dm3;
/* Computing 2nd power */
    r__2 = px3;
/* Computing 2nd power */
    r__3 = py3;
/* Computing 2nd power */
    r__4 = pz3;
    e1cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    p1beta = px3 * bg_1.betax + py3 * bg_1.betay + pz3 * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) + e1cm);
    pt1i1 = bg_1.betax * transf + px3;
    pt2i1 = bg_1.betay * transf + py3;
    pt3i1 = bg_1.betaz * transf + pz3;
    eti1 = dm3;

    if (ianti == 1 && ee_1.lb[*i1 - 1] >= 1 && ee_1.lb[*i2 - 1] >= 1) {
	ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
	ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
	if (pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] == 25) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 27;
	} else if (pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] == 27) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 25;
	}
    }

    leadng_1.lb1 = ee_1.lb[*i1 - 1];
/* FOR DELTA2 */
/* Computing 2nd power */
    r__1 = dm4;
/* Computing 2nd power */
    r__2 = px4;
/* Computing 2nd power */
    r__3 = py4;
/* Computing 2nd power */
    r__4 = pz4;
    e2cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    p2beta = px4 * bg_1.betax + py4 * bg_1.betay + pz4 * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * p2beta / (bg_1.gamma + 1.f) + e2cm);
    pt1i2 = bg_1.betax * transf + px4;
    pt2i2 = bg_1.betay * transf + py4;
    pt3i2 = bg_1.betaz * transf + pz4;
    eti2 = dm4;
    dpi_1.lb2 = ee_1.lb[*i2 - 1];
/* assign delta1 and delta2 to i1 or i2 to keep the leadng particle */
/* behaviour */
/*              if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0)then */
    bb_1.p[*i1 * 3 - 3] = pt1i1;
    bb_1.p[*i1 * 3 - 2] = pt2i1;
    bb_1.p[*i1 * 3 - 1] = pt3i1;
    cc_1.e[*i1 - 1] = eti1;
    ee_1.lb[*i1 - 1] = leadng_1.lb1;
    bb_1.p[*i2 * 3 - 3] = pt1i2;
    bb_1.p[*i2 * 3 - 2] = pt2i2;
    bb_1.p[*i2 * 3 - 1] = pt3i2;
    cc_1.e[*i2 - 1] = eti2;
    ee_1.lb[*i2 - 1] = dpi_1.lb2;
    leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
    leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
    leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
    leadng_1.em1 = cc_1.e[*i1 - 1];
    ee_1.id[*i1 - 1] = 2;
    ee_1.id[*i2 - 1] = 2;
    id1 = ee_1.id[*i1 - 1];
    *iblock = 44;
/* GET rho'S MOMENTUM AND COORDINATES IN NUCLEUS-NUCLEUS CMS. FRAME */
/* Computing 2nd power */
    r__1 = pc_1.epion[nn_1.nnn + *irun * 150001 - 150002];
/* Computing 2nd power */
    r__2 = ppx;
/* Computing 2nd power */
    r__3 = ppy;
/* Computing 2nd power */
    r__4 = ppz;
    epcm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    ppbeta = ppx * bg_1.betax + ppy * bg_1.betay + ppz * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * ppbeta / (bg_1.gamma + 1.f) + epcm);
    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = bg_1.betax * 
	    transf + ppx;
    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = bg_1.betay * 
	    transf + ppy;
    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = bg_1.betaz * 
	    transf + ppz;
/* lin-5/2008: */
    dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = dpert_1.dpertp[*i1 - 
	    1] * dpert_1.dpertp[*i2 - 1];
/* lin-5/2008: */
/* 2004        X01 = 1.0 - 2.0 * RANART(NSEED) */
/*            Y01 = 1.0 - 2.0 * RANART(NSEED) */
/*            Z01 = 1.0 - 2.0 * RANART(NSEED) */
/*        IF ((X01*X01+Y01*Y01+Z01*Z01) .GT. 1.0) GOTO 2004 */
/*                RPION(1,NNN,IRUN)=R(1,I1)+0.5*x01 */
/*                RPION(2,NNN,IRUN)=R(2,I1)+0.5*y01 */
/*                RPION(3,NNN,IRUN)=R(3,I1)+0.5*z01 */
    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = aa_1.r__[*i1 * 3 - 
	    3];
    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = aa_1.r__[*i1 * 3 - 
	    2];
    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = aa_1.r__[*i1 * 3 - 
	    1];

    goto L90005;
/* FOR THE NN-->N+N+rho PROCESS, FIND MOMENTUM OF THE FINAL */
/* PARTICLES IN THE NUCLEUS-NUCLEUS CMS. */
L308:
    ntry1 = 0;
L126:
    pprho_(srt, &input1_1.iseed, &px3, &py3, &pz3, &dm3, &px4, &py4, &pz4, &
	    dm4, &ppx, &ppy, &ppz, &amrho, &icou1);
    ++ntry1;
    if (icou1 < 0 && ntry1 <= 20) {
	goto L126;
    }
/*       if(icou1.lt.0)return */
/* ROTATE THE MOMENTA OF PARTICLES IN THE CMS OF P1+P2 */
    rotate_(px, py, pz, &px3, &py3, &pz3);
    rotate_(px, py, pz, &px4, &py4, &pz4);
    rotate_(px, py, pz, &ppx, &ppy, &ppz);
    ++nn_1.nnn;
    arho = amrho;
/* DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE */
/* (1) FOR P+P */
    xdir = ranart_(&rndf77_1.nseed);
    if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 1) {
	if (xdir <= .5f) {
/* (1.1)P+P-->P+P+rho(0) */
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 26;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
	    ee_1.lb[*i1 - 1] = 1;
	    ee_1.lb[*i2 - 1] = 1;
	    goto L2052;
	} else {
/* (1.2)P+P -->p+n+rho(+) */
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 27;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
	    ee_1.lb[*i1 - 1] = 1;
	    ee_1.lb[*i2 - 1] = 2;
	    goto L2052;
	}
    }
/* (2)FOR N+N */
    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 && (i__2 = ee_1.lb[*i2 - 1],
	     abs(i__2)) == 2) {
	if (xdir <= .5f) {
/* (2.1)N+N-->N+N+rho(0) */
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 26;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
	    ee_1.lb[*i1 - 1] = 2;
	    ee_1.lb[*i2 - 1] = 2;
	    goto L2052;
	} else {
/* (2.2)N+N -->N+P+rho(-) */
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 25;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
	    ee_1.lb[*i1 - 1] = 1;
	    ee_1.lb[*i2 - 1] = 2;
	    goto L2052;
	}
    }
/* (3)FOR N+P */
    if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 2) {
	if (xdir <= .33f) {
/* (3.1)N+P-->N+P+rho(0) */
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 26;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
	    ee_1.lb[*i1 - 1] = 1;
	    ee_1.lb[*i2 - 1] = 2;
	    goto L2052;
/* (3.2)N+P -->P+P+rho(-) */
	} else if (xdir <= .67f && xdir > .34f) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 25;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
	    ee_1.lb[*i1 - 1] = 1;
	    ee_1.lb[*i2 - 1] = 1;
	    goto L2052;
	} else {
/* (3.3)N+P-->N+N+rho(+) */
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 27;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
	    ee_1.lb[*i1 - 1] = 2;
	    ee_1.lb[*i2 - 1] = 2;
	    goto L2052;
	}
    }
/* FIND THE MOMENTUM OF PARTICLES IN THE FINAL STATE IN THE NUCLEUS- */
/* NUCLEUS CMS. FRAME */
/*             LORENTZ-TRANSFORMATION INTO LAB FRAME FOR DELTA1 */
L2052:
/* Computing 2nd power */
    r__1 = dm3;
/* Computing 2nd power */
    r__2 = px3;
/* Computing 2nd power */
    r__3 = py3;
/* Computing 2nd power */
    r__4 = pz3;
    e1cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    p1beta = px3 * bg_1.betax + py3 * bg_1.betay + pz3 * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) + e1cm);
    pt1i1 = bg_1.betax * transf + px3;
    pt2i1 = bg_1.betay * transf + py3;
    pt3i1 = bg_1.betaz * transf + pz3;
    eti1 = dm3;

    if (ianti == 1 && ee_1.lb[*i1 - 1] >= 1 && ee_1.lb[*i2 - 1] >= 1) {
	ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
	ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
	if (pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] == 25) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 27;
	} else if (pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] == 27) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 25;
	}
    }

    leadng_1.lb1 = ee_1.lb[*i1 - 1];
/* FOR p2 */
/* Computing 2nd power */
    r__1 = dm4;
/* Computing 2nd power */
    r__2 = px4;
/* Computing 2nd power */
    r__3 = py4;
/* Computing 2nd power */
    r__4 = pz4;
    e2cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    p2beta = px4 * bg_1.betax + py4 * bg_1.betay + pz4 * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * p2beta / (bg_1.gamma + 1.f) + e2cm);
    pt1i2 = bg_1.betax * transf + px4;
    pt2i2 = bg_1.betay * transf + py4;
    pt3i2 = bg_1.betaz * transf + pz4;
    eti2 = dm4;
    dpi_1.lb2 = ee_1.lb[*i2 - 1];
/* assign p1 and p2 to i1 or i2 to keep the leadng particle */
/* behaviour */
/*              if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0)then */
    bb_1.p[*i1 * 3 - 3] = pt1i1;
    bb_1.p[*i1 * 3 - 2] = pt2i1;
    bb_1.p[*i1 * 3 - 1] = pt3i1;
    cc_1.e[*i1 - 1] = eti1;
    ee_1.lb[*i1 - 1] = leadng_1.lb1;
    bb_1.p[*i2 * 3 - 3] = pt1i2;
    bb_1.p[*i2 * 3 - 2] = pt2i2;
    bb_1.p[*i2 * 3 - 1] = pt3i2;
    cc_1.e[*i2 - 1] = eti2;
    ee_1.lb[*i2 - 1] = dpi_1.lb2;
    leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
    leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
    leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
    leadng_1.em1 = cc_1.e[*i1 - 1];
    ee_1.id[*i1 - 1] = 2;
    ee_1.id[*i2 - 1] = 2;
    id1 = ee_1.id[*i1 - 1];
    *iblock = 45;
/* GET rho'S MOMENTUM AND COORDINATES IN NUCLEUS-NUCLEUS CMS. FRAME */
/* Computing 2nd power */
    r__1 = pc_1.epion[nn_1.nnn + *irun * 150001 - 150002];
/* Computing 2nd power */
    r__2 = ppx;
/* Computing 2nd power */
    r__3 = ppy;
/* Computing 2nd power */
    r__4 = ppz;
    epcm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    ppbeta = ppx * bg_1.betax + ppy * bg_1.betay + ppz * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * ppbeta / (bg_1.gamma + 1.f) + epcm);
    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = bg_1.betax * 
	    transf + ppx;
    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = bg_1.betay * 
	    transf + ppy;
    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = bg_1.betaz * 
	    transf + ppz;
/* lin-5/2008: */
    dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = dpert_1.dpertp[*i1 - 
	    1] * dpert_1.dpertp[*i2 - 1];
/* lin-5/2008: */
/* 2005        X01 = 1.0 - 2.0 * RANART(NSEED) */
/*            Y01 = 1.0 - 2.0 * RANART(NSEED) */
/*            Z01 = 1.0 - 2.0 * RANART(NSEED) */
/*        IF ((X01*X01+Y01*Y01+Z01*Z01) .GT. 1.0) GOTO 2005 */
/*                RPION(1,NNN,IRUN)=R(1,I1)+0.5*x01 */
/*                RPION(2,NNN,IRUN)=R(2,I1)+0.5*y01 */
/*                RPION(3,NNN,IRUN)=R(3,I1)+0.5*z01 */
    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = aa_1.r__[*i1 * 3 - 
	    3];
    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = aa_1.r__[*i1 * 3 - 
	    2];
    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = aa_1.r__[*i1 * 3 - 
	    1];

    goto L90005;
/* FOR THE NN-->p+p+omega PROCESS, FIND MOMENTUM OF THE FINAL */
/* PARTICLES IN THE NUCLEUS-NUCLEUS CMS. */
L309:
    ntry1 = 0;
L138:
    ppomga_(srt, &input1_1.iseed, &px3, &py3, &pz3, &dm3, &px4, &py4, &pz4, &
	    dm4, &ppx, &ppy, &ppz, &icou1);
    ++ntry1;
    if (icou1 < 0 && ntry1 <= 20) {
	goto L138;
    }
/*       if(icou1.lt.0)return */
/* ROTATE THE MOMENTA OF PARTICLES IN THE CMS OF P1+P2 */
    rotate_(px, py, pz, &px3, &py3, &pz3);
    rotate_(px, py, pz, &px4, &py4, &pz4);
    rotate_(px, py, pz, &ppx, &ppy, &ppz);
    ++nn_1.nnn;
    aomega = .782f;
/* DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE */
/* (1) FOR P+P */
    if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 1) {
/* (1.1)P+P-->P+P+omega(0) */
	pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 28;
	pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = aomega;
	ee_1.lb[*i1 - 1] = 1;
	ee_1.lb[*i2 - 1] = 1;
	goto L2053;
    }
/* (2)FOR N+N */
    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 && (i__2 = ee_1.lb[*i2 - 1],
	     abs(i__2)) == 2) {
/* (2.1)N+N-->N+N+omega(0) */
	pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 28;
	pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = aomega;
	ee_1.lb[*i1 - 1] = 2;
	ee_1.lb[*i2 - 1] = 2;
	goto L2053;
    }
/* (3)FOR N+P */
    if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 2) {
/* (3.1)N+P-->N+P+omega(0) */
	pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 28;
	pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = aomega;
	ee_1.lb[*i1 - 1] = 1;
	ee_1.lb[*i2 - 1] = 2;
	goto L2053;
    }
/* FIND THE MOMENTUM OF PARTICLES IN THE FINAL STATE IN THE NUCLEUS- */
/* NUCLEUS CMS. FRAME */
/*             LORENTZ-TRANSFORMATION INTO LAB FRAME FOR DELTA1 */
L2053:
/* Computing 2nd power */
    r__1 = dm3;
/* Computing 2nd power */
    r__2 = px3;
/* Computing 2nd power */
    r__3 = py3;
/* Computing 2nd power */
    r__4 = pz3;
    e1cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    p1beta = px3 * bg_1.betax + py3 * bg_1.betay + pz3 * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) + e1cm);
    pt1i1 = bg_1.betax * transf + px3;
    pt2i1 = bg_1.betay * transf + py3;
    pt3i1 = bg_1.betaz * transf + pz3;
    eti1 = dm3;
    if (ianti == 1 && ee_1.lb[*i1 - 1] >= 1 && ee_1.lb[*i2 - 1] >= 1) {
	ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
	ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
    }
    leadng_1.lb1 = ee_1.lb[*i1 - 1];
/* FOR DELTA2 */
/* Computing 2nd power */
    r__1 = dm4;
/* Computing 2nd power */
    r__2 = px4;
/* Computing 2nd power */
    r__3 = py4;
/* Computing 2nd power */
    r__4 = pz4;
    e2cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    p2beta = px4 * bg_1.betax + py4 * bg_1.betay + pz4 * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * p2beta / (bg_1.gamma + 1.f) + e2cm);
    pt1i2 = bg_1.betax * transf + px4;
    pt2i2 = bg_1.betay * transf + py4;
    pt3i2 = bg_1.betaz * transf + pz4;
    eti2 = dm4;
    dpi_1.lb2 = ee_1.lb[*i2 - 1];
/* assign delta1 and delta2 to i1 or i2 to keep the leadng particle */
/* behaviour */
/*              if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0)then */
    bb_1.p[*i1 * 3 - 3] = pt1i1;
    bb_1.p[*i1 * 3 - 2] = pt2i1;
    bb_1.p[*i1 * 3 - 1] = pt3i1;
    cc_1.e[*i1 - 1] = eti1;
    ee_1.lb[*i1 - 1] = leadng_1.lb1;
    bb_1.p[*i2 * 3 - 3] = pt1i2;
    bb_1.p[*i2 * 3 - 2] = pt2i2;
    bb_1.p[*i2 * 3 - 1] = pt3i2;
    cc_1.e[*i2 - 1] = eti2;
    ee_1.lb[*i2 - 1] = dpi_1.lb2;
    leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
    leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
    leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
    leadng_1.em1 = cc_1.e[*i1 - 1];
    ee_1.id[*i1 - 1] = 2;
    ee_1.id[*i2 - 1] = 2;
    id1 = ee_1.id[*i1 - 1];
    *iblock = 46;
/* GET omega'S MOMENTUM AND COORDINATES IN NUCLEUS-NUCLEUS CMS. FRAME */
/* Computing 2nd power */
    r__1 = pc_1.epion[nn_1.nnn + *irun * 150001 - 150002];
/* Computing 2nd power */
    r__2 = ppx;
/* Computing 2nd power */
    r__3 = ppy;
/* Computing 2nd power */
    r__4 = ppz;
    epcm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    ppbeta = ppx * bg_1.betax + ppy * bg_1.betay + ppz * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * ppbeta / (bg_1.gamma + 1.f) + epcm);
    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = bg_1.betax * 
	    transf + ppx;
    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = bg_1.betay * 
	    transf + ppy;
    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = bg_1.betaz * 
	    transf + ppz;
/* lin-5/2008: */
    dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = dpert_1.dpertp[*i1 - 
	    1] * dpert_1.dpertp[*i2 - 1];
/* lin-5/2008: */
/* 2006        X01 = 1.0 - 2.0 * RANART(NSEED) */
/*            Y01 = 1.0 - 2.0 * RANART(NSEED) */
/*            Z01 = 1.0 - 2.0 * RANART(NSEED) */
/*        IF ((X01*X01+Y01*Y01+Z01*Z01) .GT. 1.0) GOTO 2006 */
/*                RPION(1,NNN,IRUN)=R(1,I1)+0.5*x01 */
/*                RPION(2,NNN,IRUN)=R(2,I1)+0.5*y01 */
/*                RPION(3,NNN,IRUN)=R(3,I1)+0.5*z01 */
    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = aa_1.r__[*i1 * 3 - 
	    3];
    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = aa_1.r__[*i1 * 3 - 
	    2];
    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = aa_1.r__[*i1 * 3 - 
	    1];

    goto L90005;
/* change phase space density FOR NUCLEONS AFTER THE PROCESS */
/* lin-10/25/02-comment out following, since there is no path to it: */
/* lin-8/16/02 used before set */
/*     IX1,IY1,IZ1,IPX1,IPY1,IPZ1, IX2,IY2,IZ2,IPX2,IPY2,IPZ2: */
/*                if ((abs(ix1).le.mx) .and. (abs(iy1).le.my) .and. */
/*     &              (abs(iz1).le.mz)) then */
/*                  ipx1p = nint(p(1,i1)/dpx) */
/*                  ipy1p = nint(p(2,i1)/dpy) */
/*                  ipz1p = nint(p(3,i1)/dpz) */
/*                  if ((ipx1p.ne.ipx1) .or. (ipy1p.ne.ipy1) .or. */
/*     &                (ipz1p.ne.ipz1)) then */
/*                    if ((abs(ipx1).le.mpx) .and. (abs(ipy1).le.my) */
/*     &                .and. (ipz1.ge.-mpz) .and. (ipz1.le.mpzp)) */
/*     &                f(ix1,iy1,iz1,ipx1,ipy1,ipz1) = */
/*     &                f(ix1,iy1,iz1,ipx1,ipy1,ipz1) - 1. */
/*                    if ((abs(ipx1p).le.mpx) .and. (abs(ipy1p).le.my) */
/*     &                .and. (ipz1p.ge.-mpz).and. (ipz1p.le.mpzp)) */
/*     &                f(ix1,iy1,iz1,ipx1p,ipy1p,ipz1p) = */
/*     &                f(ix1,iy1,iz1,ipx1p,ipy1p,ipz1p) + 1. */
/*                  end if */
/*                end if */
/*                if ((abs(ix2).le.mx) .and. (abs(iy2).le.my) .and. */
/*     &              (abs(iz2).le.mz)) then */
/*                  ipx2p = nint(p(1,i2)/dpx) */
/*                  ipy2p = nint(p(2,i2)/dpy) */
/*                  ipz2p = nint(p(3,i2)/dpz) */
/*                  if ((ipx2p.ne.ipx2) .or. (ipy2p.ne.ipy2) .or. */
/*     &                (ipz2p.ne.ipz2)) then */
/*                    if ((abs(ipx2).le.mpx) .and. (abs(ipy2).le.my) */
/*     &                .and. (ipz2.ge.-mpz) .and. (ipz2.le.mpzp)) */
/*     &                f(ix2,iy2,iz2,ipx2,ipy2,ipz2) = */
/*     &                f(ix2,iy2,iz2,ipx2,ipy2,ipz2) - 1. */
/*                    if ((abs(ipx2p).le.mpx) .and. (abs(ipy2p).le.my) */
/*     &                .and. (ipz2p.ge.-mpz) .and. (ipz2p.le.mpzp)) */
/*     &                f(ix2,iy2,iz2,ipx2p,ipy2p,ipz2p) = */
/*     &                f(ix2,iy2,iz2,ipx2p,ipy2p,ipz2p) + 1. */
/*                  end if */
/*                end if */
/* lin-10/25/02-end */
L90005:
    return 0;
/* ----------------------------------------------------------------------- */
/* COM: SET THE NEW MOMENTUM COORDINATES */
L107:
    if (*px == 0.f && *py == 0.f) {
	t2 = 0.f;
    } else {
	t2 = atan2(*py, *px);
    }
/* Computing 2nd power */
    r__1 = c1;
    s1 = 1.f - r__1 * r__1;
    if (s1 <= 0.f) {
	s1 = 0.f;
    }
    s1 = sqrt(s1);
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__1 = c2;
    scheck = 1.f - r__1 * r__1;
    if (scheck < 0.f) {
	s_wsle(&io___620);
	do_lio(&c__9, &c__1, "scheck3: ", (ftnlen)9);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	scheck = 0.f;
    }
    s2 = sqrt(scheck);
/*       S2  =  SQRT( 1.0 - C2**2 ) */
    ct1 = cos(t1);
    st1 = sin(t1);
    ct2 = cos(t2);
    st2 = sin(t2);
    *pz = pr * (c1 * c2 - s1 * s2 * ct1);
    ss = c2 * s1 * ct1 + s2 * c1;
    *px = pr * (ss * ct2 - s1 * st1 * st2);
    *py = pr * (ss * st2 + s1 * st1 * ct2);
    return 0;
} /* crnn_ */

/* lin-5/2008 CRNN over */
/* ********************************* */
/* ********************************* */
/*                                                                      * */
/*                                                                      * */

/* Subroutine */ int crpp_(real *px, real *py, real *pz, real *srt, integer *
	i1, integer *i2, integer *iblock, real *ppel, real *ppin, real *
	spprho, integer *ipp)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    double sqrt(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static real c1, s1, t1, pr, ei1, ei2;
    static integer lb1, lb2;
    static real em1, em2, ct1, pr2, px0, py0, pz0, st1;
    static integer lbb1, lbb2, lb1i, lb2i, ntag;
    extern /* Subroutine */ int pi2et2_(integer *, integer *, integer *, 
	    integer *, real *, real *, integer *, integer *);
    static real ranpi;
    extern /* Subroutine */ int pi2ro2_(integer *, integer *, integer *, 
	    integer *, real *, real *, integer *, integer *), ro2et2_(integer 
	    *, integer *, integer *, integer *, real *, real *, integer *, 
	    integer *), pi3eta_(integer *, integer *, integer *, integer *, 
	    real *, real *, integer *, integer *), bbarfs_(integer *, integer 
	    *, real *, real *, integer *, integer *);
    extern doublereal ranart_(integer *);
    extern /* Subroutine */ int rotate_(real *, real *, real *, real *, real *
	    , real *), opioet_(integer *, integer *, integer *, integer *, 
	    real *, real *, integer *, integer *), rhores_(integer *, integer 
	    *), rpiret_(integer *, integer *, integer *, integer *, real *, 
	    real *, integer *, integer *);

/*     PURPOSE:                                                         * */
/*             DEALING WITH PION-PION COLLISIONS                        * */
/*     NOTE   :                                                         * */
/*           VALID ONLY FOR PION-PION-DISTANCES LESS THAN 2.5 FM        * */
/*     QUANTITIES:                                                 * */
/*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME* */
/*           SRT      - SQRT OF S                                       * */
/*           IBLOCK   - THE INFORMATION BACK                            * */
/*                     6-> Meson+Meson elastic */
/*                     66-> Meson+meson-->K+K- */
/* ********************************* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /input1/ */
/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
/* c      SAVE /RNDF77/ */
    lb1i = ee_1.lb[*i1 - 1];
    lb2i = ee_1.lb[*i2 - 1];
    px0 = *px;
    py0 = *py;
    pz0 = *pz;
    *iblock = 1;
/* ----------------------------------------------------------------------- */
/* check Meson+Meson inelastic collisions */
/* lin-9/28/00 */
/*        if((srt.gt.1.).and.(ppin/(ppin+ppel).gt.RANART(NSEED)))then */
/*        iblock=66 */
/*        e(i1)=0.498 */
/*        e(i2)=0.498 */
/*        lb(i1)=21 */
/*        lb(i2)=23 */
/*        go to 10 */
/* lin-11/07/00 */
/*        if(srt.gt.1.and.(ppin/(ppin+ppel)).gt.RANART(NSEED)) then */
/* lin-4/03/02 */
    if (*srt > .996f && *ppin / (*ppin + *ppel) > ranart_(&rndf77_1.nseed)) {
/*        if(ppin/(ppin+ppel).gt.RANART(NSEED)) then */
/* lin-10/08/00 */
	ranpi = ranart_(&rndf77_1.nseed);
	if (ppmm_1.pprr / *ppin >= ranpi) {
/*     1) pi pi <-> rho rho: */
	    pi2ro2_(i1, i2, &lbb1, &lbb2, &ei1, &ei2, iblock, &input1_1.iseed)
		    ;
/* lin-4/03/02 eta equilibration: */
	} else if ((ppmm_1.pprr + ppmm_1.ppee) / *ppin >= ranpi) {
/*     4) pi pi <-> eta eta: */
	    pi2et2_(i1, i2, &lbb1, &lbb2, &ei1, &ei2, iblock, &input1_1.iseed)
		    ;
	} else if ((ppmm_1.pprr + ppmm_1.ppee + ppmm_1.pppe) / *ppin >= ranpi)
		 {
/*     5) pi pi <-> pi eta: */
	    pi3eta_(i1, i2, &lbb1, &lbb2, &ei1, &ei2, iblock, &input1_1.iseed)
		    ;
	} else if ((ppmm_1.pprr + ppmm_1.ppee + ppmm_1.pppe + ppmm_1.rpre) / *
		ppin >= ranpi) {
/*     6) rho pi <-> pi eta: */
	    rpiret_(i1, i2, &lbb1, &lbb2, &ei1, &ei2, iblock, &input1_1.iseed)
		    ;
	} else if ((ppmm_1.pprr + ppmm_1.ppee + ppmm_1.pppe + ppmm_1.rpre + 
		ppmm_1.xopoe) / *ppin >= ranpi) {
/*     7) omega pi <-> omega eta: */
	    opioet_(i1, i2, &lbb1, &lbb2, &ei1, &ei2, iblock, &input1_1.iseed)
		    ;
	} else if ((ppmm_1.pprr + ppmm_1.ppee + ppmm_1.pppe + ppmm_1.rpre + 
		ppmm_1.xopoe + ppmm_1.rree) / *ppin >= ranpi) {
/*     8) rho rho <-> eta eta: */
	    ro2et2_(i1, i2, &lbb1, &lbb2, &ei1, &ei2, iblock, &input1_1.iseed)
		    ;
/* lin-4/03/02-end */
/*     2) BBbar production: */
	} else if ((ppmm_1.pprr + ppmm_1.ppee + ppmm_1.pppe + ppmm_1.rpre + 
		ppmm_1.xopoe + ppmm_1.rree + ppb1_1.ppinnb) / *ppin >= ranpi) 
		{
	    bbarfs_(&lbb1, &lbb2, &ei1, &ei2, iblock, &input1_1.iseed);
/*     3) KKbar production: */
	} else {
	    *iblock = 66;
	    ei1 = .498f;
	    ei2 = .498f;
	    lbb1 = 21;
	    lbb2 = 23;
/* lin-11/07/00 pi rho -> K* Kbar and K*bar K productions: */
	    lb1 = ee_1.lb[*i1 - 1];
	    lb2 = ee_1.lb[*i2 - 1];
/* lin-2/13/03 include omega the same as rho, eta the same as pi: */
/*        if(((lb1.ge.3.and.lb1.le.5).and.(lb2.ge.25.and.lb2.le.27)) */
/*     1  .or.((lb2.ge.3.and.lb2.le.5).and.(lb1.ge.25.and.lb1.le.27))) */
	    if ((lb1 == 0 || lb1 >= 3 && lb1 <= 5) && (lb2 >= 25 && lb2 <= 28)
		     || (lb2 == 0 || lb2 >= 3 && lb2 <= 5) && (lb1 >= 25 && 
		    lb1 <= 28)) {
		ei1 = .895f;
		ei2 = .498f;
		if (ranart_(&rndf77_1.nseed) >= .5f) {
		    *iblock = 366;
		    lbb1 = 30;
		    lbb2 = 21;
		} else {
		    *iblock = 367;
		    lbb1 = -30;
		    lbb2 = 23;
		}
	    }
/* lin-11/07/00-end */
	}
/* lin-ppbar-8/25/00 */
	cc_1.e[*i1 - 1] = ei1;
	cc_1.e[*i2 - 1] = ei2;
	ee_1.lb[*i1 - 1] = lbb1;
	ee_1.lb[*i2 - 1] = lbb2;
/* lin-10/08/00-end */
    } else {
/* bzdbg10/15/99 */
/* .....for meson+meson elastic srt.le.2Mk, if not pi+pi collision return */
	if ((ee_1.lb[*i1 - 1] < 3 || ee_1.lb[*i1 - 1] > 5) && (ee_1.lb[*i2 - 
		1] < 3 || ee_1.lb[*i2 - 1] > 5)) {
	    return 0;
	}
/* bzdbg10/15/99 end */
/* check Meson+Meson elastic collisions */
	*iblock = 6;
/* direct process */
	if (*ipp == 1 || *ipp == 4 || *ipp == 6) {
	    goto L10;
	}
	if (*spprho / *ppel > ranart_(&rndf77_1.nseed)) {
	    goto L20;
	}
    }
L10:
    ntag = 0;
    em1 = cc_1.e[*i1 - 1];
    em2 = cc_1.e[*i2 - 1];
/* ----------------------------------------------------------------------- */
/* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH */
/* ENERGY CONSERVATION */
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = em1;
/* Computing 2nd power */
    r__4 = em2;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = em1 * em2;
    pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
    if (pr2 <= 0.f) {
	pr2 = 1e-9f;
    }
    pr = sqrt(pr2) / (*srt * 2.f);
    c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
    t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
/* Computing 2nd power */
    r__1 = c1;
    s1 = sqrt(1.f - r__1 * r__1);
    ct1 = cos(t1);
    st1 = sin(t1);
    *pz = pr * c1;
    *px = pr * s1 * ct1;
    *py = pr * s1 * st1;
/* for isotropic distribution no need to ROTATE THE MOMENTUM */
/* ROTATE IT */
    rotate_(&px0, &py0, &pz0, px, py, pz);
    return 0;
L20:
    *iblock = 666;
/* treat rho formation in pion+pion collisions */
/* calculate the mass and momentum of rho in the nucleus-nucleus frame */
    rhores_(i1, i2);
    if (*ipp == 2) {
	ee_1.lb[*i1 - 1] = 27;
    }
    if (*ipp == 3) {
	ee_1.lb[*i1 - 1] = 26;
    }
    if (*ipp == 5) {
	ee_1.lb[*i1 - 1] = 25;
    }
    return 0;
} /* crpp_ */

/* ********************************* */
/* ********************************* */
/*                                                                      * */
/*                                                                      * */
/* Subroutine */ int crnd_(integer *irun, real *px, real *py, real *pz, real *
	srt, integer *i1, integer *i2, integer *iblock, real *signn, real *
	sig, real *sigk, real *xsk1, real *xsk2, real *xsk3, real *xsk4, real 
	*xsk5, integer *nt, integer *ipert1)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8;
    real r__1, r__2, r__3, r__4;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal), log(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double atan2(doublereal, doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    extern /* Subroutine */ int bbdangle_(real *, real *, real *, integer *, 
	    integer *, integer *, integer *, real *, real *, integer *);
    static real a, x, c1, c2, s1, t1, t2, s2, x1;
    static integer ic, m12, n12;
    static real dm, fm, as, ta, pr, ss, cc1;
    extern doublereal fd5_(real *, real *, real *);
    static real dm3, dm4, pf2, ct1, ct2;
    static integer id1;
    static real am1, am2, st1, st2, pz2, px3, py3, pz3, px4, py4, pz4;
    static integer lbd;
    extern doublereal ang_(real *, integer *);
    static integer lbm;
    extern /* Subroutine */ int m1535_(integer *, integer *, real *, real *);
    static real x1440, ppd[30000]	/* was [3][10000] */, x1535;
    extern doublereal fns_(real *, real *, real *);
    static real prf, xmm, pxd, pyd;
    extern doublereal ptr_(real *, integer *);
    static real ppx, ppy, ppz, pzd, e1cm, e2cm;
    static integer lbi1, lbi2;
    static real eti1, eti2, pxi1, pyi1, pzi1, pxi2, pyi2, pzi2;
    static integer lbpd[10000];
    static real epcm, dmin__;
    static integer ntag;
    static real dmax__, pt1i1, pt2i1, pt3i1, pt1i2, pt2i2, pt3i2, e1dcm, xptr;
    static integer icou1, ntry1, ntry2;
    extern /* Subroutine */ int sbbdm_(real *, real *, integer *, integer *, 
	    real *, real *);
    extern doublereal sigma_(real *, integer *, integer *, integer *), denom_(
	    real *, real *);
    static integer ianti;
    static real signd, sigdn, renom, xmass, p1beta, p2beta, e2picm, dprob1, 
	    renom1;
    extern /* Subroutine */ int bbkaon_(integer *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, integer *);
    static real scheck, deltam, pfinal, ppbeta;
    extern doublereal ranart_(integer *);
    static real sdprod, renomn, xdmass;
    extern /* Subroutine */ int rotate_(real *, real *, real *, real *, real *
	    , real *);
    static real transf;
    static integer ndloop, idloop, ipertd;
    static real p1dbeta, p2pibeta;

    /* Fortran I/O blocks */
    static cilist io___686 = { 0, 99, 0, 0, 0 };
    static cilist io___688 = { 0, 99, 0, 0, 0 };
    static cilist io___690 = { 0, 99, 0, 0, 0 };


/*     PURPOSE:                                                         * */
/*             DEALING WITH NUCLEON-BARYON RESONANCE COLLISIONS         * */
/*     NOTE   :                                                         * */
/*           VALID ONLY FOR BARYON-BARYON-DISTANCES LESS THAN 1.32 FM   * */
/*           (1.32 = 2 * HARD-CORE-RADIUS [HRC] )                       * */
/*     QUANTITIES:                                                 * */
/*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME* */
/*           SRT      - SQRT OF S                                       * */
/*           NSTAR =1 INCLUDING N* RESORANCE,ELSE NOT                   * */
/*           NDIRCT=1 INCLUDING DIRECT PION PRODUCTION PROCESS         * */
/*           IBLOCK   - THE INFORMATION BACK                            * */
/*                      0-> COLLISION CANNOT HAPPEN                     * */
/*                      1-> N-N ELASTIC COLLISION                       * */
/*                      2-> N+N->N+DELTA,OR N+N->N+N* REACTION          * */
/*                      3-> N+DELTA->N+N OR N+N*->N+N REACTION          * */
/*                      4-> N+N->N+N+PION,DIRTCT PROCESS                * */
/*           N12       - IS USED TO SPECIFY BARYON-BARYON REACTION      * */
/*                      CHANNELS. M12 IS THE REVERSAL CHANNEL OF N12    * */
/*                      N12,                                            * */
/*                      M12=1 FOR p+n-->delta(+)+ n                     * */
/*                          2     p+n-->delta(0)+ p                     * */
/*                          3     p+p-->delta(++)+n                     * */
/*                          4     p+p-->delta(+)+p                      * */
/*                          5     n+n-->delta(0)+n                      * */
/*                          6     n+n-->delta(-)+p                      * */
/*                          7     n+p-->N*(0)(1440)+p                   * */
/*                          8     n+p-->N*(+)(1440)+n                   * */
/*                        9     p+p-->N*(+)(1535)+p                     * */
/*                        10    n+n-->N*(0)(1535)+n                     * */
/*                         11    n+p-->N*(+)(1535)+n                     * */
/*                        12    n+p-->N*(0)(1535)+p */
/*                        13    D(++)+D(-)-->N*(+)(1440)+n */
/*                         14    D(++)+D(-)-->N*(0)(1440)+p */
/*                        15    D(+)+D(0)--->N*(+)(1440)+n */
/*                        16    D(+)+D(0)--->N*(0)(1440)+p */
/*                        17    D(++)+D(0)-->N*(+)(1535)+p */
/*                        18    D(++)+D(-)-->N*(0)(1535)+p */
/*                        19    D(++)+D(-)-->N*(+)(1535)+n */
/*                        20    D(+)+D(+)-->N*(+)(1535)+p */
/*                        21    D(+)+D(0)-->N*(+)(1535)+n */
/*                        22    D(+)+D(0)-->N*(0)(1535)+p */
/*                        23    D(+)+D(-)-->N*(0)(1535)+n */
/*                        24    D(0)+D(0)-->N*(0)(1535)+n */
/*                          25    N*(+)(14)+N*(+)(14)-->N*(+)(15)+p */
/*                          26    N*(0)(14)+N*(0)(14)-->N*(0)(15)+n */
/*                          27    N*(+)(14)+N*(0)(14)-->N*(+)(15)+n */
/*                        28    N*(+)(14)+N*(0)(14)-->N*(0)(15)+p */
/*                        29    N*(+)(14)+D+-->N*(+)(15)+p */
/*                        30    N*(+)(14)+D0-->N*(+)(15)+n */
/*                        31    N*(+)(14)+D(-)-->N*(0)(1535)+n */
/*                        32    N*(0)(14)+D++--->N*(+)(15)+p */
/*                        33    N*(0)(14)+D+--->N*(+)(15)+n */
/*                        34    N*(0)(14)+D+--->N*(0)(15)+p */
/*                        35    N*(0)(14)+D0-->N*(0)(15)+n */
/*                        36    N*(+)(14)+D0--->N*(0)(15)+p */
/*                        ++    see the note book for more listing */
/* ********************************* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /ff/ */
/* c      SAVE /gg/ */
/* c      SAVE /INPUT/ */
/* c      SAVE /NN/ */
/* c      SAVE /BG/ */
/* c      SAVE /RUN/ */
/* c      SAVE /PA/ */
/* c      SAVE /PB/ */
/* c      SAVE /PC/ */
/* c      SAVE /PD/ */
/* c      SAVE /input1/ */
/* c      SAVE /leadng/ */
/* c      SAVE /RNDF77/ */
/* ----------------------------------------------------------------------- */
    n12 = 0;
    m12 = 0;
    *iblock = 0;
    ntag = 0;
    leadng_1.em1 = cc_1.e[*i1 - 1];
    dpi_1.em2 = cc_1.e[*i2 - 1];
/* Computing 2nd power */
    r__1 = *px;
/* Computing 2nd power */
    r__2 = *py;
/* Computing 2nd power */
    r__3 = *pz;
    pr = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
    c2 = *pz / pr;
    x1 = ranart_(&rndf77_1.nseed);
    ianti = 0;
    if (ee_1.lb[*i1 - 1] < 0 && ee_1.lb[*i2 - 1] < 0) {
	ianti = 1;
    }
/* lin-6/2008 Production of perturbative deuterons for idpert=1: */
    sbbdm_(srt, &sdprod, &ianti, &lbm, &xmm, &pfinal);
    if (para8_1.idpert == 1 && *ipert1 == 1) {
	if (*srt < 2.012f) {
	    return 0;
	}
	if (((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1 || (i__2 = ee_1.lb[*i1 
		- 1], abs(i__2)) == 2) && ((i__3 = ee_1.lb[*i2 - 1], abs(i__3)
		) >= 6 && (i__4 = ee_1.lb[*i2 - 1], abs(i__4)) <= 13)) {
	    goto L108;
	} else if (((i__1 = ee_1.lb[*i2 - 1], abs(i__1)) == 1 || (i__2 = 
		ee_1.lb[*i2 - 1], abs(i__2)) == 2) && ((i__3 = ee_1.lb[*i1 - 
		1], abs(i__3)) >= 6 && (i__4 = ee_1.lb[*i1 - 1], abs(i__4)) <=
		 13)) {
	    goto L108;
	} else {
	    return 0;
	}
    }
/* ----------------------------------------------------------------------- */
/* COM: TEST FOR ELASTIC SCATTERING (EITHER N-N OR DELTA-DELTA 0R */
/*      N-DELTA OR N*-N* or N*-Delta) */
    if (x1 <= *signn / *sig) {
/* COM:  PARAMETRISATION IS TAKEN FROM THE CUGNON-PAPER */
/* Computing 6th power */
	r__1 = (*srt - 1.8766f) * 3.65f, r__1 *= r__1;
	as = r__1 * (r__1 * r__1);
	a = as * 6.f / (as + 1.f);
/* Computing 2nd power */
	r__1 = pr;
	ta = r__1 * r__1 * -2.f;
	x = ranart_(&rndf77_1.nseed);
/* lin-10/24/02        T1  = ALOG( (1-X) * DEXP(dble(A)*dble(TA)) + X )  /  A */
	t1 = (real) log((doublereal) (1.f - x) * exp((doublereal) a * (
		doublereal) ta) + (doublereal) x) / a;
	c1 = 1.f - t1 / ta;
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	*iblock = 1;
	goto L107;
    } else {
/* COM: TEST FOR INELASTIC SCATTERING */
/*     IF THE AVAILABLE ENERGY IS LESS THAN THE PION-MASS, NOTHING */
/*     CAN HAPPEN ANY MORE ==> RETURN (2.04 = 2*AVMASS + PI-MASS+0.02) */
	if (*srt < 2.04f) {
	    return 0;
	}
/* lin-6/2008 add d+meson production for n*N*(0)(1440) and p*N*(+)(1440) channels */
/*     (they did not have any inelastic reactions before): */
	if (((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 || (i__2 = ee_1.lb[*i2 
		- 1], abs(i__2)) == 2) && ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] 
		== 20 || ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 13) {
	    if (x1 <= (*signn + sdprod) / *sig) {
		goto L108;
	    }
	}

/* Resonance absorption or Delta + N-->N*(1440), N*(1535) */
/* COM: TEST FOR DELTA OR N* ABSORPTION */
/*      IN THE PROCESS DELTA+N-->NN, N*+N-->NN */
/* Computing 2nd power */
	r__1 = *srt;
	prf = sqrt(r__1 * r__1 * .25f - .88040689000000005f);
	if (leadng_1.em1 > 1.f) {
	    deltam = leadng_1.em1;
	} else {
	    deltam = dpi_1.em2;
	}
/* Computing 2nd power */
	r__1 = prf;
	renom = deltam * (r__1 * r__1) / denom_(srt, &c_b183) / pr;
/* Computing 2nd power */
	r__1 = prf;
	renomn = deltam * (r__1 * r__1) / denom_(srt, &c_b254) / pr;
/* Computing 2nd power */
	r__1 = prf;
	renom1 = deltam * (r__1 * r__1) / denom_(srt, &c_b255) / pr;
/* avoid the inelastic collisions between n+delta- -->N+N */
/*       and p+delta++ -->N+N due to charge conservation, */
/*       but they can scatter to produce kaons */
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 && (i__2 = ee_1.lb[*i2 
		- 1], abs(i__2)) == 6) {
	    renom = 0.f;
	}
	if ((i__1 = ee_1.lb[*i2 - 1], abs(i__1)) == 2 && (i__2 = ee_1.lb[*i1 
		- 1], abs(i__2)) == 6) {
	    renom = 0.f;
	}
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1 && (i__2 = ee_1.lb[*i2 
		- 1], abs(i__2)) == 9) {
	    renom = 0.f;
	}
	if ((i__1 = ee_1.lb[*i2 - 1], abs(i__1)) == 1 && (i__2 = ee_1.lb[*i1 
		- 1], abs(i__2)) == 9) {
	    renom = 0.f;
	}
	i__3 = (i__1 = ee_1.lb[*i1 - 1], abs(i__1));
	i__4 = (i__2 = ee_1.lb[*i2 - 1], abs(i__2));
	m1535_(&i__3, &i__4, srt, &x1535);
	x1440 = sigma_(srt, &c__2, &c__0, &c__1) * .75f;
/* CROSS SECTION FOR KAON PRODUCTION from the four channels */
/* for NLK channel */
/* avoid the inelastic collisions between n+delta- -->N+N */
/*       and p+delta++ -->N+N due to charge conservation, */
/*       but they can scatter to produce kaons */
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 && (i__2 = ee_1.lb[*i2 
		- 1], abs(i__2)) == 6 || (i__3 = ee_1.lb[*i2 - 1], abs(i__3)) 
		== 2 && (i__4 = ee_1.lb[*i1 - 1], abs(i__4)) == 6 || (i__5 = 
		ee_1.lb[*i1 - 1], abs(i__5)) == 1 && (i__6 = ee_1.lb[*i2 - 1],
		 abs(i__6)) == 9 || (i__7 = ee_1.lb[*i2 - 1], abs(i__7)) == 1 
		&& (i__8 = ee_1.lb[*i1 - 1], abs(i__8)) == 9) {
/* lin-6/2008 */
	    if (x1 <= (*signn + sdprod) / *sig) {
		goto L108;
	    }
/*          IF((SIGK+SIGNN)/SIG.GE.X1)GO TO 306 */
	    if ((*sigk + *signn + sdprod) / *sig >= x1) {
		goto L306;
	    }

	}
/* WE DETERMINE THE REACTION CHANNELS IN THE FOLLOWING */
/* FOR n+delta(++)-->p+p or n+delta(++)-->n+N*(+)(1440),n+N*(+)(1535) */
/* REABSORPTION OR N*(1535) PRODUCTION LIKE IN P+P OR N*(1440) LIKE PN, */
	if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 18 && ((i__1 = ee_1.lb[*i1 
		- 1], abs(i__1)) == 2 || (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) 
		== 2)) {
	    signd = sigma_(srt, &c__1, &c__1, &c__0) + sigma_(srt, &c__1, &
		    c__1, &c__1) * .5f;
	    sigdn = signd * .25f * renom;
/* lin-6/2008 */
	    if (x1 <= (*signn + sdprod) / *sig) {
		goto L108;
	    }
/*        IF(X1.GT.(SIGNN+SIGDN+X1440+X1535+SIGK)/SIG)RETURN */
	    if (x1 > (*signn + sigdn + x1440 + x1535 + *sigk + sdprod) / *sig)
		     {
		return 0;
	    }

	    if (*sigk / (*sigk + sigdn + x1440 + x1535) > ranart_(&
		    rndf77_1.nseed)) {
		goto L306;
	    }
/* REABSORPTION: */
	    if (ranart_(&rndf77_1.nseed) < sigdn / (sigdn + x1440 + x1535)) {
		m12 = 3;
		goto L206;
	    } else {
/* N* PRODUCTION */
		if (ranart_(&rndf77_1.nseed) < x1440 / (x1440 + x1535)) {
/* N*(1440) */
		    m12 = 37;
		} else {
/* N*(1535)       M12=38 */
/* lin-2/26/03 why is the above commented out? leads to M12=0 but */
/*     particle mass is changed after 204 (causes energy violation). */
/*     replace by elastic process (return): */
		    return 0;
		}
		goto L204;
	    }
	}
/* FOR p+delta(-)-->n+n or p+delta(-)-->n+N*(0)(1440),n+N*(0)(1535) */
/* REABSORPTION OR N*(1535) PRODUCTION LIKE IN P+P OR N*(1440) LIKE PN, */
	if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 6 && ((i__1 = ee_1.lb[*i1 
		- 1], abs(i__1)) == 1 || (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) 
		== 1)) {
	    signd = sigma_(srt, &c__1, &c__1, &c__0) + sigma_(srt, &c__1, &
		    c__1, &c__1) * .5f;
	    sigdn = signd * .25f * renom;
/* lin-6/2008 */
	    if (x1 <= (*signn + sdprod) / *sig) {
		goto L108;
	    }
/*        IF (X1.GT.(SIGNN+SIGDN+X1440+X1535+SIGK)/SIG)RETURN */
	    if (x1 > (*signn + sigdn + x1440 + x1535 + *sigk + sdprod) / *sig)
		     {
		return 0;
	    }

	    if (*sigk / (*sigk + sigdn + x1440 + x1535) > ranart_(&
		    rndf77_1.nseed)) {
		goto L306;
	    }
/* REABSORPTION: */
	    if (ranart_(&rndf77_1.nseed) < sigdn / (sigdn + x1440 + x1535)) {
		m12 = 6;
		goto L206;
	    } else {
/* N* PRODUCTION */
		if (ranart_(&rndf77_1.nseed) < x1440 / (x1440 + x1535)) {
/* N*(1440) */
		    m12 = 47;
		} else {
/* N*(1535)       M12=48 */
/* lin-2/26/03 causes energy violation, replace by elastic process (return): */
		    return 0;
		}
		goto L204;
	    }
	}
/* FOR p+delta(+)-->p+p, N*(+)(144)+p, N*(+)(1535)+p */
	if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 8 && ((i__1 = ee_1.lb[*i1 
		- 1], abs(i__1)) == 1 || (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) 
		== 1)) {
	    signd = sigma_(srt, &c__1, &c__1, &c__1) * 1.5f;
	    sigdn = signd * .25f * renom;
/* lin-6/2008 */
	    if (x1 <= (*signn + sdprod) / *sig) {
		goto L108;
	    }
/*        IF(X1.GT.(SIGNN+SIGDN+x1440+x1535+SIGK)/SIG)RETURN */
	    if (x1 > (*signn + sigdn + x1440 + x1535 + *sigk + sdprod) / *sig)
		     {
		return 0;
	    }

	    if (*sigk / (*sigk + sigdn + x1440 + x1535) > ranart_(&
		    rndf77_1.nseed)) {
		goto L306;
	    }
	    if (ranart_(&rndf77_1.nseed) < sigdn / (sigdn + x1440 + x1535)) {
		m12 = 4;
		goto L206;
	    } else {
		if (ranart_(&rndf77_1.nseed) < x1440 / (x1440 + x1535)) {
/* N*(144) */
		    m12 = 39;
		} else {
		    m12 = 40;
		}
		goto L204;
	    }
	}
/* FOR n+delta(0)-->n+n, N*(0)(144)+n, N*(0)(1535)+n */
	if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 14 && ((i__1 = ee_1.lb[*i1 
		- 1], abs(i__1)) == 2 || (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) 
		== 2)) {
	    signd = sigma_(srt, &c__1, &c__1, &c__1) * 1.5f;
	    sigdn = signd * .25f * renom;
/* lin-6/2008 */
	    if (x1 <= (*signn + sdprod) / *sig) {
		goto L108;
	    }
/*        IF(X1.GT.(SIGNN+SIGDN+x1440+x1535+SIGK)/SIG)RETURN */
	    if (x1 > (*signn + sigdn + x1440 + x1535 + *sigk + sdprod) / *sig)
		     {
		return 0;
	    }

	    if (*sigk / (*sigk + sigdn + x1440 + x1535) > ranart_(&
		    rndf77_1.nseed)) {
		goto L306;
	    }
	    if (ranart_(&rndf77_1.nseed) < sigdn / (sigdn + x1440 + x1535)) {
		m12 = 5;
		goto L206;
	    } else {
		if (ranart_(&rndf77_1.nseed) < x1440 / (x1440 + x1535)) {
/* N*(144) */
		    m12 = 48;
		} else {
		    m12 = 49;
		}
		goto L204;
	    }
	}
/* FOR n+delta(+)-->n+p, N*(+)(1440)+n,N*(0)(1440)+p, */
/*                       N*(+)(1535)+n,N*(0)(1535)+p */
	if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 16 && ((i__1 = ee_1.lb[*i1 
		- 1], abs(i__1)) == 2 || (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) 
		== 2)) {
	    signd = sigma_(srt, &c__1, &c__1, &c__1) * .5f + sigma_(srt, &
		    c__1, &c__1, &c__0) * .25f;
	    sigdn = signd * .5f * renom;
/* lin-6/2008 */
	    if (x1 <= (*signn + sdprod) / *sig) {
		goto L108;
	    }
/*        IF(X1.GT.(SIGNN+SIGDN+2.*x1440+2.*x1535+SIGK)/SIG)RETURN */
	    if (x1 > (*signn + sigdn + x1440 * 2.f + x1535 * 2.f + *sigk + 
		    sdprod) / *sig) {
		return 0;
	    }

	    if (*sigk / (*sigk + sigdn + x1440 * 2 + x1535 * 2) > ranart_(&
		    rndf77_1.nseed)) {
		goto L306;
	    }
	    if (ranart_(&rndf77_1.nseed) < sigdn / (sigdn + x1440 * 2.f + 
		    x1535 * 2.f)) {
		m12 = 1;
		goto L206;
	    } else {
		if (ranart_(&rndf77_1.nseed) < x1440 / (x1440 + x1535)) {
		    m12 = 41;
		    if (ranart_(&rndf77_1.nseed) <= .5f) {
			m12 = 43;
		    }
		} else {
		    m12 = 42;
		    if (ranart_(&rndf77_1.nseed) <= .5f) {
			m12 = 44;
		    }
		}
		goto L204;
	    }
	}
/* FOR p+delta(0)-->n+p, N*(+)(1440)+n,N*(0)(1440)+p, */
/*                       N*(+)(1535)+n,N*(0)(1535)+p */
	if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 7) {
	    signd = sigma_(srt, &c__1, &c__1, &c__1) * .5f + sigma_(srt, &
		    c__1, &c__1, &c__0) * .25f;
	    sigdn = signd * .5f * renom;
/* lin-6/2008 */
	    if (x1 <= (*signn + sdprod) / *sig) {
		goto L108;
	    }
/*        IF(X1.GT.(SIGNN+SIGDN+2.*x1440+2.*x1535+SIGK)/SIG)RETURN */
	    if (x1 > (*signn + sigdn + x1440 * 2.f + x1535 * 2.f + *sigk + 
		    sdprod) / *sig) {
		return 0;
	    }

	    if (*sigk / (*sigk + sigdn + x1440 * 2 + x1535 * 2) > ranart_(&
		    rndf77_1.nseed)) {
		goto L306;
	    }
	    if (ranart_(&rndf77_1.nseed) < sigdn / (sigdn + x1440 * 2.f + 
		    x1535 * 2.f)) {
		m12 = 2;
		goto L206;
	    } else {
		if (ranart_(&rndf77_1.nseed) < x1440 / (x1440 + x1535)) {
		    m12 = 50;
		    if (ranart_(&rndf77_1.nseed) <= .5f) {
			m12 = 51;
		    }
		} else {
		    m12 = 52;
		    if (ranart_(&rndf77_1.nseed) <= .5f) {
			m12 = 53;
		    }
		}
		goto L204;
	    }
	}
/* FOR p+N*(0)(14)-->p+n, N*(+)(1535)+n,N*(0)(1535)+p */
/* OR  P+N*(0)(14)-->D(+)+N, D(0)+P, */
	if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 10 && ((i__1 = ee_1.lb[*i1 
		- 1], abs(i__1)) == 1 || (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) 
		== 1)) {
	    signd = sigma_(srt, &c__2, &c__0, &c__1) * .75f;
	    sigdn = signd * renomn;
/* lin-6/2008 */
	    if (x1 <= (*signn + sdprod) / *sig) {
		goto L108;
	    }
/*        IF(X1.GT.(SIGNN+SIGDN+X1535+SIGK)/SIG)RETURN */
	    if (x1 > (*signn + sigdn + x1535 + *sigk + sdprod) / *sig) {
		return 0;
	    }

	    if (*sigk / (*sigk + sigdn + x1535) > ranart_(&rndf77_1.nseed)) {
		goto L306;
	    }
	    if (ranart_(&rndf77_1.nseed) < sigdn / (sigdn + x1535)) {
		m12 = 7;
		goto L206;
	    } else {
		m12 = 54;
		if (ranart_(&rndf77_1.nseed) <= .5f) {
		    m12 = 55;
		}
	    }
	    goto L204;
	}
/* FOR n+N*(+)-->p+n, N*(+)(1535)+n,N*(0)(1535)+p */
	if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 22 && ((i__1 = ee_1.lb[*i1 
		- 1], abs(i__1)) == 2 || (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) 
		== 2)) {
	    signd = sigma_(srt, &c__2, &c__0, &c__1) * .75f;
	    sigdn = signd * renomn;
/* lin-6/2008 */
	    if (x1 <= (*signn + sdprod) / *sig) {
		goto L108;
	    }
/*        IF(X1.GT.(SIGNN+SIGDN+X1535+SIGK)/SIG)RETURN */
	    if (x1 > (*signn + sigdn + x1535 + *sigk + sdprod) / *sig) {
		return 0;
	    }

	    if (*sigk / (*sigk + sigdn + x1535) > ranart_(&rndf77_1.nseed)) {
		goto L306;
	    }
	    if (ranart_(&rndf77_1.nseed) < sigdn / (sigdn + x1535)) {
		m12 = 8;
		goto L206;
	    } else {
		m12 = 56;
		if (ranart_(&rndf77_1.nseed) <= .5f) {
		    m12 = 57;
		}
	    }
	    goto L204;
	}
/* FOR N*(1535)+N-->N+N COLLISIONS */
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 12 || (i__2 = ee_1.lb[*i1 
		- 1], abs(i__2)) == 13 || (i__3 = ee_1.lb[*i2 - 1], abs(i__3))
		 == 12 || (i__4 = ee_1.lb[*i2 - 1], abs(i__4)) == 13) {
	    signd = x1535;
	    sigdn = signd * renom1;
/* lin-6/2008 */
	    if (x1 <= (*signn + sdprod) / *sig) {
		goto L108;
	    }
/*        IF(X1.GT.(SIGNN+SIGDN+SIGK)/SIG)RETURN */
	    if (x1 > (*signn + sigdn + *sigk + sdprod) / *sig) {
		return 0;
	    }

	    if (*sigk / (*sigk + sigdn) > ranart_(&rndf77_1.nseed)) {
		goto L306;
	    }
	    if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 24) {
		m12 = 10;
	    }
	    if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 12) {
		m12 = 12;
	    }
	    if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 26) {
		m12 = 11;
	    }
	    if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 13) {
		m12 = 9;
	    }
	    goto L206;
	}
L204:
/* (1) GENERATE THE MASS FOR THE N*(1440) AND N*(1535) */
/* (2) CALCULATE THE FINAL MOMENTUM OF THE n+N* SYSTEM */
/* (3) RELABLE THE FINAL STATE PARTICLES */
/* PARAMETRIZATION OF THE SHAPE OF THE N* RESONANCE ACCORDING */
/*     TO kitazoe's or J.D.JACKSON'S MASS FORMULA AND BREIT WIGNER */
/*     FORMULA FOR N* RESORANCE */
/*     DETERMINE DELTA MASS VIA REJECTION METHOD. */
	dmax__ = *srt - .9383f - .005f;
	dmin__ = 1.078f;
	if (m12 == 37 || m12 == 39 || m12 == 41 || m12 == 43 || m12 == 46 || 
		m12 == 48 || m12 == 50 || m12 == 51) {
/* N*(1440) production */
	    if (dmax__ < 1.44f) {
		fm = fns_(&dmax__, srt, &c_b182);
	    } else {
/* lin-10/25/02 get rid of argument usage mismatch in FNS(): */
		xdmass = 1.44f;
/*          FM=FNS(1.44,SRT,1.) */
		fm = fns_(&xdmass, srt, &c_b183);
/* lin-10/25/02-end */
	    }
	    if (fm == 0.f) {
		fm = 1e-9f;
	    }
	    ntry2 = 0;
L11:
	    dm = ranart_(&rndf77_1.nseed) * (dmax__ - dmin__) + dmin__;
	    ++ntry2;
	    if (ranart_(&rndf77_1.nseed) > fns_(&dm, srt, &c_b183) / fm && 
		    ntry2 <= 10) {
		goto L11;
	    }
/* lin-2/26/03 limit the N* mass below a certain value */
/*     (here taken as its central value + 2* B-W fullwidth): */
	    if (dm > 2.14f) {
		goto L11;
	    }
	    goto L13;
	} else {
/* N*(1535) production */
	    if (dmax__ < 1.535f) {
		fm = fd5_(&dmax__, srt, &c_b182);
	    } else {
/* lin-10/25/02 get rid of argument usage mismatch in FNS(): */
		xdmass = 1.535f;
/*          FM=FD5(1.535,SRT,1.) */
		fm = fd5_(&xdmass, srt, &c_b183);
/* lin-10/25/02-end */
	    }
	    if (fm == 0.f) {
		fm = 1e-9f;
	    }
	    ntry1 = 0;
L12:
	    dm = ranart_(&rndf77_1.nseed) * (dmax__ - dmin__) + dmin__;
	    ++ntry1;
	    if (ranart_(&rndf77_1.nseed) > fd5_(&dm, srt, &c_b183) / fm && 
		    ntry1 <= 10) {
		goto L12;
	    }
/* lin-2/26/03 limit the N* mass below a certain value */
/*     (here taken as its central value + 2* B-W fullwidth): */
	    if (dm > 1.84f) {
		goto L12;
	    }
	}
L13:
/* (2) DETERMINE THE FINAL MOMENTUM */
	prf = 0.f;
/* Computing 2nd power */
	r__2 = *srt;
/* Computing 2nd power */
	r__3 = dm;
/* Computing 2nd power */
	r__1 = (r__2 * r__2 - r__3 * r__3 + .88040689000000005f) / (*srt * 
		2.f);
	pf2 = r__1 * r__1 - .88040689000000005f;
	if (pf2 > 0.f) {
	    prf = sqrt(pf2);
	}
/* (3) RELABLE FINAL STATE PARTICLES */
/* 37 D(++)+n-->N*(+)(14)+p */
	if (m12 == 37) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 9) {
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
		ee_1.lb[*i2 - 1] = 11;
		cc_1.e[*i2 - 1] = dm;
	    } else {
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .93828f;
		ee_1.lb[*i1 - 1] = 11;
		cc_1.e[*i1 - 1] = dm;
	    }
	    goto L207;
	}
/* 38 D(++)+n-->N*(+)(15)+p */
	if (m12 == 38) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 9) {
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
		ee_1.lb[*i2 - 1] = 13;
		cc_1.e[*i2 - 1] = dm;
	    } else {
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .93828f;
		ee_1.lb[*i1 - 1] = 13;
		cc_1.e[*i1 - 1] = dm;
	    }
	    goto L207;
	}
/* 39 D(+)+P-->N*(+)(14)+p */
	if (m12 == 39) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8) {
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
		ee_1.lb[*i2 - 1] = 11;
		cc_1.e[*i2 - 1] = dm;
	    } else {
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .93828f;
		ee_1.lb[*i1 - 1] = 11;
		cc_1.e[*i1 - 1] = dm;
	    }
	    goto L207;
	}
/* 40 D(+)+P-->N*(+)(15)+p */
	if (m12 == 40) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8) {
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
		ee_1.lb[*i2 - 1] = 13;
		cc_1.e[*i2 - 1] = dm;
	    } else {
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .93828f;
		ee_1.lb[*i1 - 1] = 13;
		cc_1.e[*i1 - 1] = dm;
	    }
	    goto L207;
	}
/* 41 D(+)+N-->N*(+)(14)+N */
	if (m12 == 41) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8) {
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
		ee_1.lb[*i2 - 1] = 11;
		cc_1.e[*i2 - 1] = dm;
	    } else {
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
		ee_1.lb[*i1 - 1] = 11;
		cc_1.e[*i1 - 1] = dm;
	    }
	    goto L207;
	}
/* 42 D(+)+N-->N*(+)(15)+N */
	if (m12 == 42) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8) {
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
		ee_1.lb[*i2 - 1] = 13;
		cc_1.e[*i2 - 1] = dm;
	    } else {
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
		ee_1.lb[*i1 - 1] = 13;
		cc_1.e[*i1 - 1] = dm;
	    }
	    goto L207;
	}
/* 43 D(+)+N-->N*(0)(14)+P */
	if (m12 == 43) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8) {
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
		ee_1.lb[*i2 - 1] = 10;
		cc_1.e[*i2 - 1] = dm;
	    } else {
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .93828f;
		ee_1.lb[*i1 - 1] = 10;
		cc_1.e[*i1 - 1] = dm;
	    }
	    goto L207;
	}
/* 44 D(+)+N-->N*(0)(15)+P */
	if (m12 == 44) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8) {
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
		ee_1.lb[*i2 - 1] = 12;
		cc_1.e[*i2 - 1] = dm;
	    } else {
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .93828f;
		ee_1.lb[*i1 - 1] = 12;
		cc_1.e[*i1 - 1] = dm;
	    }
	    goto L207;
	}
/* 46 D(-)+P-->N*(0)(14)+N */
	if (m12 == 46) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 6) {
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
		ee_1.lb[*i2 - 1] = 10;
		cc_1.e[*i2 - 1] = dm;
	    } else {
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
		ee_1.lb[*i1 - 1] = 10;
		cc_1.e[*i1 - 1] = dm;
	    }
	    goto L207;
	}
/* 47 D(-)+P-->N*(0)(15)+N */
	if (m12 == 47) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 6) {
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
		ee_1.lb[*i2 - 1] = 12;
		cc_1.e[*i2 - 1] = dm;
	    } else {
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
		ee_1.lb[*i1 - 1] = 12;
		cc_1.e[*i1 - 1] = dm;
	    }
	    goto L207;
	}
/* 48 D(0)+N-->N*(0)(14)+N */
	if (m12 == 48) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7) {
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
		ee_1.lb[*i2 - 1] = 11;
		cc_1.e[*i2 - 1] = dm;
	    } else {
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
		ee_1.lb[*i1 - 1] = 11;
		cc_1.e[*i1 - 1] = dm;
	    }
	    goto L207;
	}
/* 49 D(0)+N-->N*(0)(15)+N */
	if (m12 == 49) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7) {
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
		ee_1.lb[*i2 - 1] = 12;
		cc_1.e[*i2 - 1] = dm;
	    } else {
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
		ee_1.lb[*i1 - 1] = 12;
		cc_1.e[*i1 - 1] = dm;
	    }
	    goto L207;
	}
/* 50 D(0)+P-->N*(0)(14)+P */
	if (m12 == 50) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7) {
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
		ee_1.lb[*i2 - 1] = 10;
		cc_1.e[*i2 - 1] = dm;
	    } else {
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .93828f;
		ee_1.lb[*i1 - 1] = 10;
		cc_1.e[*i1 - 1] = dm;
	    }
	    goto L207;
	}
/* 51 D(0)+P-->N*(+)(14)+N */
	if (m12 == 51) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7) {
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
		ee_1.lb[*i2 - 1] = 11;
		cc_1.e[*i2 - 1] = dm;
	    } else {
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
		ee_1.lb[*i1 - 1] = 11;
		cc_1.e[*i1 - 1] = dm;
	    }
	    goto L207;
	}
/* 52 D(0)+P-->N*(0)(15)+P */
	if (m12 == 52) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7) {
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
		ee_1.lb[*i2 - 1] = 12;
		cc_1.e[*i2 - 1] = dm;
	    } else {
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .93828f;
		ee_1.lb[*i1 - 1] = 12;
		cc_1.e[*i1 - 1] = dm;
	    }
	    goto L207;
	}
/* 53 D(0)+P-->N*(+)(15)+N */
	if (m12 == 53) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7) {
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
		ee_1.lb[*i2 - 1] = 13;
		cc_1.e[*i2 - 1] = dm;
	    } else {
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
		ee_1.lb[*i1 - 1] = 13;
		cc_1.e[*i1 - 1] = dm;
	    }
	    goto L207;
	}
/* 54 N*(0)(14)+P-->N*(+)(15)+N */
	if (m12 == 54) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 10) {
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
		ee_1.lb[*i2 - 1] = 13;
		cc_1.e[*i2 - 1] = dm;
	    } else {
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
		ee_1.lb[*i1 - 1] = 13;
		cc_1.e[*i1 - 1] = dm;
	    }
	    goto L207;
	}
/* 55 N*(0)(14)+P-->N*(0)(15)+P */
	if (m12 == 55) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 10) {
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
		ee_1.lb[*i2 - 1] = 12;
		cc_1.e[*i2 - 1] = dm;
	    } else {
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .93828f;
		ee_1.lb[*i1 - 1] = 12;
		cc_1.e[*i1 - 1] = dm;
	    }
	    goto L207;
	}
/* 56 N*(+)(14)+N-->N*(+)(15)+N */
	if (m12 == 56) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 11) {
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
		ee_1.lb[*i2 - 1] = 13;
		cc_1.e[*i2 - 1] = dm;
	    } else {
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
		ee_1.lb[*i1 - 1] = 13;
		cc_1.e[*i1 - 1] = dm;
	    }
	    goto L207;
	}
/* 57 N*(+)(14)+N-->N*(0)(15)+P */
	if (m12 == 57) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 11) {
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
		ee_1.lb[*i2 - 1] = 12;
		cc_1.e[*i2 - 1] = dm;
	    } else {
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .93828f;
		ee_1.lb[*i1 - 1] = 12;
		cc_1.e[*i1 - 1] = dm;
	    }
	}
	goto L207;
/* ------------------------------------------------ */
/* RELABLE NUCLEONS AFTER DELTA OR N* BEING ABSORBED */
/* (1) n+delta(+)-->n+p */
L206:
	if (m12 == 1) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8) {
		ee_1.lb[*i2 - 1] = 2;
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
	    } else {
		ee_1.lb[*i1 - 1] = 2;
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .93828f;
	    }
	    goto L207;
	}
/* (2) p+delta(0)-->p+n */
	if (m12 == 2) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7) {
		ee_1.lb[*i2 - 1] = 1;
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
	    } else {
		ee_1.lb[*i1 - 1] = 1;
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
	    }
	    goto L207;
	}
/* (3) n+delta(++)-->p+p */
	if (m12 == 3) {
	    ee_1.lb[*i1 - 1] = 1;
	    ee_1.lb[*i2 - 1] = 1;
	    cc_1.e[*i1 - 1] = .93828f;
	    cc_1.e[*i2 - 1] = .93828f;
	    goto L207;
	}
/* (4) p+delta(+)-->p+p */
	if (m12 == 4) {
	    ee_1.lb[*i1 - 1] = 1;
	    ee_1.lb[*i2 - 1] = 1;
	    cc_1.e[*i1 - 1] = .93828f;
	    cc_1.e[*i2 - 1] = .93828f;
	    goto L207;
	}
/* (5) n+delta(0)-->n+n */
	if (m12 == 5) {
	    ee_1.lb[*i1 - 1] = 2;
	    ee_1.lb[*i2 - 1] = 2;
	    cc_1.e[*i1 - 1] = .939457f;
	    cc_1.e[*i2 - 1] = .939457f;
	    goto L207;
	}
/* (6) p+delta(-)-->n+n */
	if (m12 == 6) {
	    ee_1.lb[*i1 - 1] = 2;
	    ee_1.lb[*i2 - 1] = 2;
	    cc_1.e[*i1 - 1] = .939457f;
	    cc_1.e[*i2 - 1] = .939457f;
	    goto L207;
	}
/* (7) p+N*(0)-->n+p */
	if (m12 == 7) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1) {
		ee_1.lb[*i1 - 1] = 1;
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i1 - 1] = .93828f;
		cc_1.e[*i2 - 1] = .939457f;
	    } else {
		ee_1.lb[*i1 - 1] = 2;
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i1 - 1] = .939457f;
		cc_1.e[*i2 - 1] = .93828f;
	    }
	    goto L207;
	}
/* (8) n+N*(+)-->n+p */
	if (m12 == 8) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2) {
		ee_1.lb[*i1 - 1] = 2;
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i1 - 1] = .939457f;
		cc_1.e[*i2 - 1] = .93828f;
	    } else {
		ee_1.lb[*i1 - 1] = 1;
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i1 - 1] = .93828f;
		cc_1.e[*i2 - 1] = .939457f;
	    }
	    goto L207;
	}
/* lin-6/2008 */
/* *(9) N*(+)p-->pp */
/* (9) N*(+)(1535) p-->pp */
	if (m12 == 9) {
	    ee_1.lb[*i1 - 1] = 1;
	    ee_1.lb[*i2 - 1] = 1;
	    cc_1.e[*i1 - 1] = .93828f;
	    cc_1.e[*i2 - 1] = .93828f;
	    goto L207;
	}
/* (12) N*(0)P-->nP */
	if (m12 == 12) {
	    ee_1.lb[*i1 - 1] = 2;
	    ee_1.lb[*i2 - 1] = 1;
	    cc_1.e[*i1 - 1] = .939457f;
	    cc_1.e[*i2 - 1] = .93828f;
	    goto L207;
	}
/* (11) N*(+)n-->nP */
	if (m12 == 11) {
	    ee_1.lb[*i1 - 1] = 2;
	    ee_1.lb[*i2 - 1] = 1;
	    cc_1.e[*i1 - 1] = .939457f;
	    cc_1.e[*i2 - 1] = .93828f;
	    goto L207;
	}
/* lin-6/2008 */
/* *(12) N*(0)p-->Np */
/* (12) N*(0)(1535) p-->Np */
	if (m12 == 12) {
	    ee_1.lb[*i1 - 1] = 1;
	    ee_1.lb[*i2 - 1] = 2;
	    cc_1.e[*i1 - 1] = .93828f;
	    cc_1.e[*i2 - 1] = .939457f;
	}
/* ---------------------------------------------- */
L207:
	pr = prf;
	c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	if (*srt <= 2.14f) {
	    c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	}
	if (*srt > 2.14f && *srt <= 2.4f) {
	    c1 = ang_(srt, &input1_1.iseed);
	}
	if (*srt > 2.4f) {
/* lin-10/25/02 get rid of argument usage mismatch in PTR(): */
	    xptr = pr * .33f;
/*         cc1=ptr(0.33*pr,iseed) */
	    cc1 = ptr_(&xptr, &input1_1.iseed);
/* lin-10/25/02-end */
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
	    r__1 = pr;
/* Computing 2nd power */
	    r__2 = cc1;
	    scheck = r__1 * r__1 - r__2 * r__2;
	    if (scheck < 0.f) {
		s_wsle(&io___686);
		do_lio(&c__9, &c__1, "scheck4: ", (ftnlen)9);
		do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
		e_wsle();
		scheck = 0.f;
	    }
	    c1 = sqrt(scheck) / pr;
/*         c1=sqrt(pr**2-cc1**2)/pr */
	}
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	*iblock = 3;
    }
    if (ianti == 1 && ee_1.lb[*i1 - 1] >= 1 && ee_1.lb[*i2 - 1] >= 1) {
	ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
	ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
    }
/* ----------------------------------------------------------------------- */
/* COM: SET THE NEW MOMENTUM COORDINATES */
L107:
    if (*px == 0.f && *py == 0.f) {
	t2 = 0.f;
    } else {
	t2 = atan2(*py, *px);
    }
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__1 = c1;
    scheck = 1.f - r__1 * r__1;
    if (scheck < 0.f) {
	s_wsle(&io___688);
	do_lio(&c__9, &c__1, "scheck5: ", (ftnlen)9);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	scheck = 0.f;
    }
    s1 = sqrt(scheck);
/*      S1   = SQRT( 1.0 - C1**2 ) */
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__1 = c2;
    scheck = 1.f - r__1 * r__1;
    if (scheck < 0.f) {
	s_wsle(&io___690);
	do_lio(&c__9, &c__1, "scheck6: ", (ftnlen)9);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	scheck = 0.f;
    }
    s2 = sqrt(scheck);
/*      S2  =  SQRT( 1.0 - C2**2 ) */
    ct1 = cos(t1);
    st1 = sin(t1);
    ct2 = cos(t2);
    st2 = sin(t2);
    *pz = pr * (c1 * c2 - s1 * s2 * ct1);
    ss = c2 * s1 * ct1 + s2 * c1;
    *px = pr * (ss * ct2 - s1 * st1 * st2);
    *py = pr * (ss * st2 + s1 * st1 * ct2);
    return 0;
/* FOR THE NN-->KAON+X PROCESS, FIND MOMENTUM OF THE FINAL PARTICLES IN */
/* THE NUCLEUS-NUCLEUS CMS. */
L306:
/* sp11/21/01 phi production */
    if (*xsk5 / *sigk > ranart_(&rndf77_1.nseed)) {
	leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
	pz2 = bb_1.p[*i2 * 3 - 1];
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	++nn_1.nnn;
	pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 29;
	pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = 1.02f;
	*iblock = 222;
	goto L208;
    }
/* sp11/21/01 end */
    *iblock = 11;
    if (ianti == 1) {
	*iblock = -11;
    }

    leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
    pz2 = bb_1.p[*i2 * 3 - 1];
/* DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE */
    ++nn_1.nnn;
    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 23;
    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .498f;
    if (*srt <= 2.63f) {
/* only lambda production is possible */
/* (1.1)P+P-->p+L+kaon+ */
	ic = 1;
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	ee_1.lb[*i2 - 1] = 14;
	goto L208;
    }
    if (*srt <= 2.74f && *srt > 2.63f) {
/* both Lambda and sigma production are possible */
	if (*xsk1 / (*xsk1 + *xsk2) > ranart_(&rndf77_1.nseed)) {
/* lambda production */
	    ic = 1;
	    ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	    ee_1.lb[*i2 - 1] = 14;
	} else {
/* sigma production */
	    ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	    ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 15;
	    ic = 2;
	}
	goto L208;
    }
    if (*srt <= 2.77f && *srt > 2.74f) {
/* then pp-->Delta lamda kaon can happen */
	if (*xsk1 / (*xsk1 + *xsk2 + *xsk3) > ranart_(&rndf77_1.nseed)) {
/* * (1.1)P+P-->p+L+kaon+ */
	    ic = 1;
	    ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	    ee_1.lb[*i2 - 1] = 14;
	    goto L208;
	} else {
	    if (*xsk2 / (*xsk2 + *xsk3) > ranart_(&rndf77_1.nseed)) {
/* pp-->psk */
		ic = 2;
		ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 
			1;
		ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 
			15;
	    } else {
/* pp-->D+l+k */
		ic = 3;
		ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 4) + 
			6;
		ee_1.lb[*i2 - 1] = 14;
	    }
	    goto L208;
	}
    }
    if (*srt > 2.77f) {
/* all four channels are possible */
	if (*xsk1 / (*xsk1 + *xsk2 + *xsk3 + *xsk4) > ranart_(&rndf77_1.nseed)
		) {
/* p lambda k production */
	    ic = 1;
	    ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	    ee_1.lb[*i2 - 1] = 14;
	    goto L208;
	} else {
	    if (*xsk3 / (*xsk2 + *xsk3 + *xsk4) > ranart_(&rndf77_1.nseed)) {
/* delta l K production */
		ic = 3;
		ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 4) + 
			6;
		ee_1.lb[*i2 - 1] = 14;
		goto L208;
	    } else {
		if (*xsk2 / (*xsk2 + *xsk4) > ranart_(&rndf77_1.nseed)) {
/* n sigma k production */
		    ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 
			    2) + 1;
		    ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 
			    3) + 15;
		    ic = 2;
		} else {
		    ic = 4;
		    ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 
			    4) + 6;
		    ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 
			    3) + 15;
		}
		goto L208;
	    }
	}
    }
L208:
    if (ianti == 1 && ee_1.lb[*i1 - 1] >= 1 && ee_1.lb[*i2 - 1] >= 1) {
	ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
	ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
	if (pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] == 23) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 21;
	}
    }
    lbi1 = ee_1.lb[*i1 - 1];
    lbi2 = ee_1.lb[*i2 - 1];
/* KEEP ALL COORDINATES OF PARTICLE 2 FOR POSSIBLE PHASE SPACE CHANGE */
    ntry1 = 0;
L128:
    bbkaon_(&ic, srt, &px3, &py3, &pz3, &dm3, &px4, &py4, &pz4, &dm4, &ppx, &
	    ppy, &ppz, &icou1);
    ++ntry1;
    if (icou1 < 0 && ntry1 <= 20) {
	goto L128;
    }
/*       if(icou1.lt.0)return */
/* ROTATE THE MOMENTA OF PARTICLES IN THE CMS OF P1+P2 */
    rotate_(px, py, pz, &px3, &py3, &pz3);
    rotate_(px, py, pz, &px4, &py4, &pz4);
    rotate_(px, py, pz, &ppx, &ppy, &ppz);
/* FIND THE MOMENTUM OF PARTICLES IN THE FINAL STATE IN THE NUCLEUS- */
/* NUCLEUS CMS. FRAME */
/* (1) for the necleon/delta */
/*             LORENTZ-TRANSFORMATION INTO LAB FRAME FOR DELTA1 */
/* Computing 2nd power */
    r__1 = dm3;
/* Computing 2nd power */
    r__2 = px3;
/* Computing 2nd power */
    r__3 = py3;
/* Computing 2nd power */
    r__4 = pz3;
    e1cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    p1beta = px3 * bg_1.betax + py3 * bg_1.betay + pz3 * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) + e1cm);
    pt1i1 = bg_1.betax * transf + px3;
    pt2i1 = bg_1.betay * transf + py3;
    pt3i1 = bg_1.betaz * transf + pz3;
    eti1 = dm3;
/* (2) for the lambda/sigma */
/* Computing 2nd power */
    r__1 = dm4;
/* Computing 2nd power */
    r__2 = px4;
/* Computing 2nd power */
    r__3 = py4;
/* Computing 2nd power */
    r__4 = pz4;
    e2cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    p2beta = px4 * bg_1.betax + py4 * bg_1.betay + pz4 * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * p2beta / (bg_1.gamma + 1.f) + e2cm);
    pt1i2 = bg_1.betax * transf + px4;
    pt2i2 = bg_1.betay * transf + py4;
    pt3i2 = bg_1.betaz * transf + pz4;
    eti2 = dm4;
/* GET the kaon'S MOMENTUM AND COORDINATES IN NUCLEUS-NUCLEUS CMS. FRAME */
/* Computing 2nd power */
    r__1 = ppx;
/* Computing 2nd power */
    r__2 = ppy;
/* Computing 2nd power */
    r__3 = ppz;
    epcm = sqrt(r__1 * r__1 + .248004f + r__2 * r__2 + r__3 * r__3);
    ppbeta = ppx * bg_1.betax + ppy * bg_1.betay + ppz * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * ppbeta / (bg_1.gamma + 1.f) + epcm);
    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = bg_1.betax * 
	    transf + ppx;
    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = bg_1.betay * 
	    transf + ppy;
    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = bg_1.betaz * 
	    transf + ppz;
/* lin-5/2008: */
    dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = dpert_1.dpertp[*i1 - 
	    1] * dpert_1.dpertp[*i2 - 1];
/* lin-5/2008: */
/* 2008        X01 = 1.0 - 2.0 * RANART(NSEED) */
/*            Y01 = 1.0 - 2.0 * RANART(NSEED) */
/*            Z01 = 1.0 - 2.0 * RANART(NSEED) */
/*        IF ((X01*X01+Y01*Y01+Z01*Z01) .GT. 1.0) GOTO 2008 */
/*                RPION(1,NNN,IRUN)=R(1,I1)+0.5*x01 */
/*                RPION(2,NNN,IRUN)=R(2,I1)+0.5*y01 */
/*                RPION(3,NNN,IRUN)=R(3,I1)+0.5*z01 */
    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = aa_1.r__[*i1 * 3 - 
	    3];
    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = aa_1.r__[*i1 * 3 - 
	    2];
    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = aa_1.r__[*i1 * 3 - 
	    1];

/* assign the nucleon/delta and lambda/sigma to i1 or i2 to keep the */
/* leadng particle behaviour */
/*              if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0)then */
    bb_1.p[*i1 * 3 - 3] = pt1i1;
    bb_1.p[*i1 * 3 - 2] = pt2i1;
    bb_1.p[*i1 * 3 - 1] = pt3i1;
    cc_1.e[*i1 - 1] = eti1;
    ee_1.lb[*i1 - 1] = lbi1;
    bb_1.p[*i2 * 3 - 3] = pt1i2;
    bb_1.p[*i2 * 3 - 2] = pt2i2;
    bb_1.p[*i2 * 3 - 1] = pt3i2;
    cc_1.e[*i2 - 1] = eti2;
    ee_1.lb[*i2 - 1] = lbi2;
    leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
    leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
    leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
    leadng_1.em1 = cc_1.e[*i1 - 1];
    ee_1.id[*i1 - 1] = 2;
    ee_1.id[*i2 - 1] = 2;
    id1 = ee_1.id[*i1 - 1];
    if (pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] != 29) {
	*iblock = 11;
    }
    leadng_1.lb1 = ee_1.lb[*i1 - 1];
    dpi_1.lb2 = ee_1.lb[*i2 - 1];
    am1 = leadng_1.em1;
    am2 = dpi_1.em2;
/* Computing 2nd power */
    r__1 = leadng_1.em1;
/* Computing 2nd power */
    r__2 = leadng_1.px1;
/* Computing 2nd power */
    r__3 = leadng_1.py1;
/* Computing 2nd power */
    r__4 = leadng_1.pz1;
    leadng_1.e1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    return 0;
/* lin-6/2008 N+D->Deuteron+pi: */
/*     FIND MOMENTUM OF THE FINAL PARTICLES IN THE NUCLEUS-NUCLEUS CMS. */
L108:
    if (para8_1.idpert == 1 && *ipert1 == 1 && para8_1.npertd >= 1) {
/*     For idpert=1: we produce npertd pert deuterons: */
	ndloop = para8_1.npertd;
    } else if (para8_1.idpert == 2 && para8_1.npertd >= 1) {
/*     For idpert=2: we first save information for npertd pert deuterons; */
/*     at the last ndloop we create the regular deuteron+pi */
/*     and those pert deuterons: */
	ndloop = para8_1.npertd + 1;
    } else {
/*     Just create the regular deuteron+pi: */
	ndloop = 1;
    }

    dprob1 = sdprod / *sig / (real) para8_1.npertd;
    i__1 = ndloop;
    for (idloop = 1; idloop <= i__1; ++idloop) {
	bbdangle_(&pxd, &pyd, &pzd, nt, ipert1, &ianti, &idloop, &pfinal, &
		dprob1, &lbm);
	rotate_(px, py, pz, &pxd, &pyd, &pzd);
/*     LORENTZ-TRANSFORMATION OF THE MOMENTUM OF PARTICLES IN THE FINAL STATE */
/*     FROM THE NN CMS FRAME INTO THE GLOBAL CMS FRAME: */
/*     For the Deuteron: */
	xmass = 1.8756f;
/* Computing 2nd power */
	r__1 = xmass;
/* Computing 2nd power */
	r__2 = pxd;
/* Computing 2nd power */
	r__3 = pyd;
/* Computing 2nd power */
	r__4 = pzd;
	e1dcm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p1dbeta = pxd * bg_1.betax + pyd * bg_1.betay + pzd * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p1dbeta / (bg_1.gamma + 1.f) + 
		e1dcm);
	pxi1 = bg_1.betax * transf + pxd;
	pyi1 = bg_1.betay * transf + pyd;
	pzi1 = bg_1.betaz * transf + pzd;
	if (ianti == 0) {
	    lbd = 42;
	} else {
	    lbd = -42;
	}
	if (para8_1.idpert == 1 && *ipert1 == 1 && para8_1.npertd >= 1) {
/* ccc  Perturbative production for idpert=1: */
	    ++nn_1.nnn;
	    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = pxi1;
	    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = pyi1;
	    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = pzi1;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = 1.8756f;
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = lbd;
	    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = aa_1.r__[*
		    i1 * 3 - 3];
	    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = aa_1.r__[*
		    i1 * 3 - 2];
	    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = aa_1.r__[*
		    i1 * 3 - 1];
/* lin-6/2008 assign the perturbative probability: */
	    dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = sdprod / *
		    sig / (real) para8_1.npertd;
	} else if (para8_1.idpert == 2 && idloop <= para8_1.npertd) {
/* lin-6/2008 For idpert=2, we produce NPERTD perturbative (anti)deuterons */
/*     only when a regular (anti)deuteron+pi is produced in NN collisions. */
/*     First save the info for the perturbative deuterons: */
	    ppd[idloop * 3 - 3] = pxi1;
	    ppd[idloop * 3 - 2] = pyi1;
	    ppd[idloop * 3 - 1] = pzi1;
	    lbpd[idloop - 1] = lbd;
	} else {
/* ccc  Regular production: */
/*     For the regular pion: do LORENTZ-TRANSFORMATION: */
	    cc_1.e[*i1 - 1] = xmm;
/* Computing 2nd power */
	    r__1 = xmm;
/* Computing 2nd power */
	    r__2 = pxd;
/* Computing 2nd power */
	    r__3 = pyd;
/* Computing 2nd power */
	    r__4 = pzd;
	    e2picm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * 
		    r__4);
	    p2pibeta = -pxd * bg_1.betax - pyd * bg_1.betay - pzd * 
		    bg_1.betaz;
	    transf = bg_1.gamma * (bg_1.gamma * p2pibeta / (bg_1.gamma + 1.f) 
		    + e2picm);
	    pxi2 = bg_1.betax * transf - pxd;
	    pyi2 = bg_1.betay * transf - pyd;
	    pzi2 = bg_1.betaz * transf - pzd;
	    bb_1.p[*i1 * 3 - 3] = pxi2;
	    bb_1.p[*i1 * 3 - 2] = pyi2;
	    bb_1.p[*i1 * 3 - 1] = pzi2;
/*     Remove regular pion to check the equivalence */
/*     between the perturbative and regular deuteron results: */
/*                 E(i1)=0. */

	    ee_1.lb[*i1 - 1] = lbm;
	    leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
	    leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
	    leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
	    leadng_1.em1 = cc_1.e[*i1 - 1];
	    ee_1.id[*i1 - 1] = 2;
	    id1 = ee_1.id[*i1 - 1];
/* Computing 2nd power */
	    r__1 = leadng_1.em1;
/* Computing 2nd power */
	    r__2 = leadng_1.px1;
/* Computing 2nd power */
	    r__3 = leadng_1.py1;
/* Computing 2nd power */
	    r__4 = leadng_1.pz1;
	    leadng_1.e1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 
		    * r__4);
	    leadng_1.lb1 = ee_1.lb[*i1 - 1];
/*     For the regular deuteron: */
	    bb_1.p[*i2 * 3 - 3] = pxi1;
	    bb_1.p[*i2 * 3 - 2] = pyi1;
	    bb_1.p[*i2 * 3 - 1] = pzi1;
	    ee_1.lb[*i2 - 1] = lbd;
	    dpi_1.lb2 = ee_1.lb[*i2 - 1];
	    cc_1.e[*i2 - 1] = 1.8756f;
	    eti2 = cc_1.e[*i2 - 1];
	    ee_1.id[*i2 - 1] = 2;
/*     For idpert=2: create the perturbative deuterons: */
	    if (para8_1.idpert == 2 && idloop == ndloop) {
		i__2 = para8_1.npertd;
		for (ipertd = 1; ipertd <= i__2; ++ipertd) {
		    ++nn_1.nnn;
		    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = 
			    ppd[ipertd * 3 - 3];
		    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = 
			    ppd[ipertd * 3 - 2];
		    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = 
			    ppd[ipertd * 3 - 1];
		    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = 1.8756f;
		    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = lbpd[
			    ipertd - 1];
		    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = 
			    aa_1.r__[*i1 * 3 - 3];
		    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = 
			    aa_1.r__[*i1 * 3 - 2];
		    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = 
			    aa_1.r__[*i1 * 3 - 1];
/* lin-6/2008 assign the perturbative probability: */
		    dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = 1.f /
			     (real) para8_1.npertd;
		}
	    }
	}
    }
    *iblock = 501;
    return 0;
/* lin-6/2008 N+D->Deuteron+pi over */
} /* crnd_ */

/* ********************************* */
/*                                                                      * */
/*                                                                      * */
/* Subroutine */ int crdd_(integer *irun, real *px, real *py, real *pz, real *
	srt, integer *i1, integer *i2, integer *iblock, integer *ntag, real *
	signn, real *sig, integer *nt, integer *ipert1)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6;
    real r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    double sqrt(doublereal), atan2(doublereal, doublereal), exp(doublereal), 
	    log(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    extern /* Subroutine */ int bbdangle_(real *, real *, real *, integer *, 
	    integer *, integer *, integer *, real *, real *, integer *);
    static real a, x, c1, c2, s1, t1, t2, s2, x1;
    static integer ic;
    static real al;
    static integer m12, n12;
    static real dm, fm, as, ta, es, pr, ss, cc1, ak0;
    extern doublereal fd5_(real *, real *, real *);
    static real dm3, dm4, ct1, s2d, ct2;
    static integer id1;
    static real am1, am2, pr2, st1, st2, pz2, px3, py3, pz3, px4, py4, pz4, 
	    ada, ana;
    static integer idd, lbd, ich;
    extern doublereal ang_(real *, integer *);
    static integer lbm;
    extern /* Subroutine */ int n1535_(integer *, integer *, real *, real *);
    static real akp, ppd[30000]	/* was [3][10000] */, x1535;
    extern doublereal fns_(real *, real *, real *);
    static real pxd, xmm, pyd, pzd;
    extern doublereal ptr_(real *, integer *);
    static real ppx, ppy, ppz, e1cm, e2cm;
    static integer lbi1, lbi2;
    static real eti1, sig2, eti2;
    extern doublereal ppk0_(real *), ppk1_(real *);
    static real pxi1, pyi1, pzi1, pxi2, pyi2, xsk1, xsk2, xsk3, xsk4, xsk5, 
	    pzi2;
    static integer lbpd[10000];
    static real epcm, dmin__, dmax__, sigk, pt1i1, pt2i1, pt3i1, pt1i2, pt2i2,
	     pt3i2, e1dcm, xptr, t1dlk, t2dlk;
    static integer icou1;
    static real t1dsk, t2dsk, t1nlk, t1nsk, t2nsk;
    static integer ntry1, ntry2;
    extern /* Subroutine */ int sbbdm_(real *, real *, integer *, integer *, 
	    real *, real *);
    extern doublereal sigma_(real *, integer *, integer *, integer *);
    static integer ianti;
    static real pmdlk, signd, pmdsk;
    extern doublereal pplpk_(real *);
    static real pmnsk, xmass;
    extern doublereal reab2d_(integer *, integer *, real *);
    static real p1beta, p2beta, e2picm, dprob1, pmdlk2, pmdsk2, pmnsk2;
    extern /* Subroutine */ int bbkaon_(integer *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, integer *);
    static real scheck, pfinal, ppbeta;
    extern doublereal ranart_(integer *);
    static real sdprod, xdmass;
    extern /* Subroutine */ int rotate_(real *, real *, real *, real *, real *
	    , real *);
    static real transf;
    static integer ndloop, idloop, ipertd;
    static real p1dbeta, p2pibeta;

    /* Fortran I/O blocks */
    static cilist io___811 = { 0, 99, 0, 0, 0 };
    static cilist io___812 = { 0, 99, 0, 0, 0 };
    static cilist io___814 = { 0, 99, 0, 0, 0 };


/*     1NTAG,SIGNN,SIG) */
/*     PURPOSE:                                                         * */
/*             DEALING WITH BARYON RESONANCE-BARYON RESONANCE COLLISIONS* */
/*     NOTE   :                                                         * */
/*     QUANTITIES:                                                 * */
/*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME* */
/*           SRT      - SQRT OF S                                       * */
/*           NSTAR =1 INCLUDING N* RESORANCE,ELSE NOT                   * */
/*           NDIRCT=1 INCLUDING DIRECT PION PRODUCTION PROCESS         * */
/*           IBLOCK   - THE INFORMATION BACK                            * */
/*                      0-> COLLISION CANNOT HAPPEN                     * */
/*                      1-> N-N ELASTIC COLLISION                       * */
/*                      2-> N+N->N+DELTA,OR N+N->N+N* REACTION          * */
/*                      3-> N+DELTA->N+N OR N+N*->N+N REACTION          * */
/*                      4-> N+N->N+N+PION,DIRTCT PROCESS                * */
/*                     5-> DELTA(N*)+DELTA(N*)   TOTAL   COLLISIONS    * */
/*           N12       - IS USED TO SPECIFY BARYON-BARYON REACTION      * */
/*                      CHANNELS. M12 IS THE REVERSAL CHANNEL OF N12    * */
/*                      N12,                                            * */
/*                      M12=1 FOR p+n-->delta(+)+ n                     * */
/*                          2     p+n-->delta(0)+ p                     * */
/*                          3     p+p-->delta(++)+n                     * */
/*                          4     p+p-->delta(+)+p                      * */
/*                          5     n+n-->delta(0)+n                      * */
/*                          6     n+n-->delta(-)+p                      * */
/*                          7     n+p-->N*(0)(1440)+p                   * */
/*                          8     n+p-->N*(+)(1440)+n                   * */
/*                        9     p+p-->N*(+)(1535)+p                     * */
/*                        10    n+n-->N*(0)(1535)+n                     * */
/*                         11    n+p-->N*(+)(1535)+n                     * */
/*                        12    n+p-->N*(0)(1535)+p */
/*                        13    D(++)+D(-)-->N*(+)(1440)+n */
/*                         14    D(++)+D(-)-->N*(0)(1440)+p */
/*                        15    D(+)+D(0)--->N*(+)(1440)+n */
/*                        16    D(+)+D(0)--->N*(0)(1440)+p */
/*                        17    D(++)+D(0)-->N*(+)(1535)+p */
/*                        18    D(++)+D(-)-->N*(0)(1535)+p */
/*                        19    D(++)+D(-)-->N*(+)(1535)+n */
/*                        20    D(+)+D(+)-->N*(+)(1535)+p */
/*                        21    D(+)+D(0)-->N*(+)(1535)+n */
/*                        22    D(+)+D(0)-->N*(0)(1535)+p */
/*                        23    D(+)+D(-)-->N*(0)(1535)+n */
/*                        24    D(0)+D(0)-->N*(0)(1535)+n */
/*                          25    N*(+)(14)+N*(+)(14)-->N*(+)(15)+p */
/*                          26    N*(0)(14)+N*(0)(14)-->N*(0)(15)+n */
/*                          27    N*(+)(14)+N*(0)(14)-->N*(+)(15)+n */
/*                        28    N*(+)(14)+N*(0)(14)-->N*(0)(15)+p */
/*                        29    N*(+)(14)+D+-->N*(+)(15)+p */
/*                        30    N*(+)(14)+D0-->N*(+)(15)+n */
/*                        31    N*(+)(14)+D(-)-->N*(0)(1535)+n */
/*                        32    N*(0)(14)+D++--->N*(+)(15)+p */
/*                        33    N*(0)(14)+D+--->N*(+)(15)+n */
/*                        34    N*(0)(14)+D+--->N*(0)(15)+p */
/*                        35    N*(0)(14)+D0-->N*(0)(15)+n */
/*                        36    N*(+)(14)+D0--->N*(0)(15)+p */
/*                        +++ */
/*               AND MORE CHANNELS AS LISTED IN THE NOTE BOOK */

/* NOTE ABOUT N*(1440) RESORANCE:                                       * */
/*     As it has been discussed in VerWest's paper,I= 1 (initial isospin) */
/*     channel can all be attributed to delta resorance while I= 0      * */
/*     channel can all be  attribured to N* resorance.Only in n+p       * */
/*     one can have I=0 channel so is the N*(1440) resorance            * */
/* REFERENCES:    J. CUGNON ET AL., NUCL. PHYS. A352, 505 (1981)        * */
/*                    Y. KITAZOE ET AL., PHYS. LETT. 166B, 35 (1986)    * */
/*                    B. VerWest el al., PHYS. PRV. C25 (1982)1979      * */
/*                    Gy. Wolf  et al, Nucl Phys A517 (1990) 615        * */
/*                    CUTOFF = 2 * AVMASS + 20 MEV                      * */
/*                                                                      * */
/*       for N*(1535) we use the parameterization by Gy. Wolf et al     * */
/*       Nucl phys A552 (1993) 349, added May 18, 1994                  * */
/* ********************************* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /ff/ */
/* c      SAVE /gg/ */
/* c      SAVE /INPUT/ */
/* c      SAVE /NN/ */
/* c      SAVE /BG/ */
/* c      SAVE /RUN/ */
/* c      SAVE /PA/ */
/* c      SAVE /PB/ */
/* c      SAVE /PC/ */
/* c      SAVE /PD/ */
/* c      SAVE /input1/ */
/* c      SAVE /leadng/ */
/* c      SAVE /RNDF77/ */
/* ----------------------------------------------------------------------- */
    n12 = 0;
    m12 = 0;
    *iblock = 0;
    *ntag = 0;
    leadng_1.em1 = cc_1.e[*i1 - 1];
    dpi_1.em2 = cc_1.e[*i2 - 1];
/* Computing 2nd power */
    r__1 = *px;
/* Computing 2nd power */
    r__2 = *py;
/* Computing 2nd power */
    r__3 = *pz;
    pr = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
    c2 = *pz / pr;
    if (*px == 0.f && *py == 0.f) {
	t2 = 0.f;
    } else {
	t2 = atan2(*py, *px);
    }
    x1 = ranart_(&rndf77_1.nseed);
    ianti = 0;
    if (ee_1.lb[*i1 - 1] < 0 && ee_1.lb[*i2 - 1] < 0) {
	ianti = 1;
    }
/* lin-6/2008 Production of perturbative deuterons for idpert=1: */
    sbbdm_(srt, &sdprod, &ianti, &lbm, &xmm, &pfinal);
    if (para8_1.idpert == 1 && *ipert1 == 1) {
	if (*srt < 2.012f) {
	    return 0;
	}
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) >= 6 && (i__2 = ee_1.lb[*i1 
		- 1], abs(i__2)) <= 13 && ((i__3 = ee_1.lb[*i2 - 1], abs(i__3)
		) >= 6 && (i__4 = ee_1.lb[*i2 - 1], abs(i__4)) <= 13)) {
	    goto L108;
	} else {
	    return 0;
	}
    }
/* ----------------------------------------------------------------------- */
/* COM: TEST FOR ELASTIC SCATTERING (EITHER N-N OR DELTA-DELTA 0R */
/*      N-DELTA OR N*-N* or N*-Delta) */
    if (x1 <= *signn / *sig) {
/* COM:  PARAMETRISATION IS TAKEN FROM THE CUGNON-PAPER */
/* Computing 6th power */
	r__1 = (*srt - 1.8766f) * 3.65f, r__1 *= r__1;
	as = r__1 * (r__1 * r__1);
	a = as * 6.f / (as + 1.f);
/* Computing 2nd power */
	r__1 = pr;
	ta = r__1 * r__1 * -2.f;
	x = ranart_(&rndf77_1.nseed);
/* lin-10/24/02        T1  = DLOG( (1-X) * DEXP(dble(A)*dble(TA)) + X )  /  A */
	t1 = (real) log((doublereal) (1.f - x) * exp((doublereal) a * (
		doublereal) ta) + (doublereal) x) / a;
	c1 = 1.f - t1 / ta;
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	*iblock = 20;
	goto L107;
    } else {
/* COM: TEST FOR INELASTIC SCATTERING */
/*     IF THE AVAILABLE ENERGY IS LESS THAN THE PION-MASS, NOTHING */
/*     CAN HAPPEN ANY MORE ==> RETURN (2.15 = 2*AVMASS +2*PI-MASS) */
	if (*srt < 2.15f) {
	    return 0;
	}
/*     IF THERE WERE 2 N*(1535) AND THEY DIDN'T SCATT. ELAST., */
/*     ALLOW THEM TO PRODUCE KAONS. NO OTHER INELASTIC CHANNELS */
/*     ARE KNOWN */
/*       if((lb(i1).ge.12).and.(lb(i2).ge.12))return */
/*     ALL the inelastic collisions between N*(1535) and Delta as well */
/*     as N*(1440) TO PRODUCE KAONS, NO OTHER CHANNELS ARE KNOWN */
/*       if((lb(i1).ge.12).and.(lb(i2).ge.3))return */
/*       if((lb(i2).ge.12).and.(lb(i1).ge.3))return */
/*     calculate the N*(1535) production cross section in I1+I2 collisions */
	i__3 = (i__1 = ee_1.lb[*i1 - 1], abs(i__1));
	i__4 = (i__2 = ee_1.lb[*i2 - 1], abs(i__2));
	n1535_(&i__3, &i__4, srt, &x1535);
/* for Delta+Delta-->N*(1440 OR 1535)+N AND N*(1440)+N*(1440)-->N*(1535)+X */
/*     AND DELTA+N*(1440)-->N*(1535)+X */
/* WE ASSUME THEY HAVE THE SAME CROSS SECTIONS as CORRESPONDING N+N COLLISION): */
/* FOR D++D0, D+D+,D+D-,D0D0,N*+N*+,N*0N*0,N*(+)D+,N*(+)D(-),N*(0)D(0) */
/* N*(1535) production, kaon production and reabsorption through */
/* D(N*)+D(N*)-->NN are ALLOWED. */
/* CROSS SECTION FOR KAON PRODUCTION from the four channels are */
/* for NLK channel */
	akp = .498f;
	ak0 = .498f;
	ana = .938f;
	ada = 1.232f;
	al = 1.1157f;
	as = 1.1197f;
	xsk1 = 0.f;
	xsk2 = 0.f;
	xsk3 = 0.f;
	xsk4 = 0.f;
	xsk5 = 0.f;
	t1nlk = ana + al + akp;
	if (*srt <= t1nlk) {
	    goto L222;
	}
	xsk1 = pplpk_(srt) * 1.5f;
/* for DLK channel */
	t1dlk = ada + al + akp;
	t2dlk = ada + al - akp;
	if (*srt <= t1dlk) {
	    goto L222;
	}
	es = *srt;
/* Computing 2nd power */
	r__1 = es;
/* Computing 2nd power */
	r__2 = t1dlk;
/* Computing 2nd power */
	r__3 = es;
/* Computing 2nd power */
	r__4 = t2dlk;
/* Computing 2nd power */
	r__5 = es;
	pmdlk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (
		r__5 * r__5 * 4.f);
	pmdlk = sqrt(pmdlk2);
	xsk3 = pplpk_(srt) * 1.5f;
/* for NSK channel */
	t1nsk = ana + as + akp;
	t2nsk = ana + as - akp;
	if (*srt <= t1nsk) {
	    goto L222;
	}
/* Computing 2nd power */
	r__1 = es;
/* Computing 2nd power */
	r__2 = t1nsk;
/* Computing 2nd power */
	r__3 = es;
/* Computing 2nd power */
	r__4 = t2nsk;
/* Computing 2nd power */
	r__5 = es;
	pmnsk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (
		r__5 * r__5 * 4.f);
	pmnsk = sqrt(pmnsk2);
	xsk2 = (ppk1_(srt) + ppk0_(srt)) * 1.5f;
/* for DSK channel */
	t1dsk = ada + as + akp;
	t2dsk = ada + as - akp;
	if (*srt <= t1dsk) {
	    goto L222;
	}
/* Computing 2nd power */
	r__1 = es;
/* Computing 2nd power */
	r__2 = t1dsk;
/* Computing 2nd power */
	r__3 = es;
/* Computing 2nd power */
	r__4 = t2dsk;
/* Computing 2nd power */
	r__5 = es;
	pmdsk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (
		r__5 * r__5 * 4.f);
	pmdsk = sqrt(pmdsk2);
	xsk4 = (ppk1_(srt) + ppk0_(srt)) * 1.5f;
/* sp11/21/01 */
/* phi production */
	if (*srt <= 2.898914f) {
	    goto L222;
	}
/*  !! mb put the correct form */
	xsk5 = 1e-4f;
/* sp11/21/01 end */
/* THE TOTAL KAON+ PRODUCTION CROSS SECTION IS THEN */
L222:
	sigk = xsk1 + xsk2 + xsk3 + xsk4;
/* bz3/7/99 neutralk */
	xsk1 *= 2.f;
	xsk2 *= 2.f;
	xsk3 *= 2.f;
	xsk4 *= 2.f;
	sigk = sigk * 2.f + xsk5;
/* bz3/7/99 neutralk end */
/* The reabsorption cross section for the process */
/* D(N*)D(N*)-->NN is */
	s2d = reab2d_(i1, i2, srt);
/* bz3/16/99 pion */
	s2d = 0.f;
/* bz3/16/99 pion end */
/* (1) N*(1535)+D(N*(1440)) reactions */
/*    we allow kaon production and reabsorption only */
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) >= 12 && (i__2 = ee_1.lb[*i2 
		- 1], abs(i__2)) >= 12 || (i__3 = ee_1.lb[*i1 - 1], abs(i__3))
		 >= 12 && (i__4 = ee_1.lb[*i2 - 1], abs(i__4)) >= 6 || (i__5 =
		 ee_1.lb[*i2 - 1], abs(i__5)) >= 12 && (i__6 = ee_1.lb[*i1 - 
		1], abs(i__6)) >= 6) {
	    signd = sigk + s2d;
/* lin-6/2008 */
	    if (x1 <= (*signn + sdprod) / *sig) {
		goto L108;
	    }
/*       if(x1.gt.(signd+signn)/sig)return */
	    if (x1 > (signd + *signn + sdprod) / *sig) {
		return 0;
	    }

/* if kaon production */
/* lin-6/2008 */
/*       IF(SIGK/SIG.GE.RANART(NSEED))GO TO 306 */
	    if ((sigk + sdprod) / *sig >= ranart_(&rndf77_1.nseed)) {
		goto L306;
	    }

/* if reabsorption */
	    goto L1012;
	}
	idd = (i__1 = ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1], abs(i__1));
/* channels have the same charge as pp */
	if (idd == 63 || idd == 64 || idd == 48 || idd == 49 || idd == 121 || 
		idd == 100 || idd == 88 || idd == 66 || idd == 90 || idd == 
		70) {
	    signd = x1535 + sigk + s2d;
/* lin-6/2008 */
	    if (x1 <= (*signn + sdprod) / *sig) {
		goto L108;
	    }
/*        IF (X1.GT.(SIGNN+SIGND)/SIG)RETURN */
	    if (x1 > (*signn + signd + sdprod) / *sig) {
		return 0;
	    }

/* if kaon production */
	    if (sigk / signd > ranart_(&rndf77_1.nseed)) {
		goto L306;
	    }
/* if reabsorption */
	    if (s2d / (x1535 + s2d) > ranart_(&rndf77_1.nseed)) {
		goto L1012;
	    }
/* if N*(1535) production */
	    if (idd == 63) {
		n12 = 17;
	    }
	    if (idd == 64) {
		n12 = 20;
	    }
	    if (idd == 48) {
		n12 = 23;
	    }
	    if (idd == 49) {
		n12 = 24;
	    }
	    if (idd == 121) {
		n12 = 25;
	    }
	    if (idd == 100) {
		n12 = 26;
	    }
	    if (idd == 88) {
		n12 = 29;
	    }
	    if (idd == 66) {
		n12 = 31;
	    }
	    if (idd == 90) {
		n12 = 32;
	    }
	    if (idd == 70) {
		n12 = 35;
	    }
	    goto L1011;
	}
/* IN DELTA+N*(1440) and N*(1440)+N*(1440) COLLISIONS, */
/* N*(1535), kaon production and reabsorption are ALLOWED */
/* IN N*(1440)+N*(1440) COLLISIONS, ONLY N*(1535) IS ALLOWED */
	if (idd == 110 || idd == 77 || idd == 80) {
/* lin-6/2008 */
	    if (x1 <= (*signn + sdprod) / *sig) {
		goto L108;
	    }
/*       IF(X1.GT.(SIGNN+X1535+SIGK+s2d)/SIG)RETURN */
	    if (x1 > (*signn + x1535 + sigk + s2d + sdprod) / *sig) {
		return 0;
	    }

	    if (sigk / (x1535 + sigk + s2d) > ranart_(&rndf77_1.nseed)) {
		goto L306;
	    }
	    if (s2d / (x1535 + s2d) > ranart_(&rndf77_1.nseed)) {
		goto L1012;
	    }
	    if (idd == 77) {
		n12 = 30;
	    }
	    if (idd == 77 && ranart_(&rndf77_1.nseed) <= .5f) {
		n12 = 36;
	    }
	    if (idd == 80) {
		n12 = 34;
	    }
	    if (idd == 80 && ranart_(&rndf77_1.nseed) <= .5f) {
		n12 = 35;
	    }
	    if (idd == 110) {
		n12 = 27;
	    }
	    if (idd == 110 && ranart_(&rndf77_1.nseed) <= .5f) {
		n12 = 28;
	    }
	    goto L1011;
	}
	if (idd == 54 || idd == 56) {
/* LIKE FOR N+P COLLISION, */
/* IN DELTA+DELTA COLLISIONS BOTH N*(1440) AND N*(1535) CAN BE PRODUCED */
	    sig2 = sigma_(srt, &c__2, &c__0, &c__1) * .75f;
	    signd = (sig2 + x1535) * 2.f + sigk + s2d;
/* lin-6/2008 */
	    if (x1 <= (*signn + sdprod) / *sig) {
		goto L108;
	    }
/*        IF(X1.GT.(SIGNN+SIGND)/SIG)RETURN */
	    if (x1 > (*signn + signd + sdprod) / *sig) {
		return 0;
	    }

	    if (sigk / signd > ranart_(&rndf77_1.nseed)) {
		goto L306;
	    }
	    if (s2d / ((sig2 + x1535) * 2.f + s2d) > ranart_(&rndf77_1.nseed))
		     {
		goto L1012;
	    }
	    if (ranart_(&rndf77_1.nseed) < x1535 / (sig2 + x1535)) {
/* N*(1535) PRODUCTION */
		if (idd == 54) {
		    n12 = 18;
		}
		if (idd == 54 && ranart_(&rndf77_1.nseed) <= .5f) {
		    n12 = 19;
		}
		if (idd == 56) {
		    n12 = 21;
		}
		if (idd == 56 && ranart_(&rndf77_1.nseed) <= .5f) {
		    n12 = 22;
		}
	    } else {
/* N*(144) PRODUCTION */
		if (idd == 54) {
		    n12 = 13;
		}
		if (idd == 54 && ranart_(&rndf77_1.nseed) <= .5f) {
		    n12 = 14;
		}
		if (idd == 56) {
		    n12 = 15;
		}
		if (idd == 56 && ranart_(&rndf77_1.nseed) <= .5f) {
		    n12 = 16;
		}
	    }
	}
L1011:
	*iblock = 5;
/* PARAMETRIZATION OF THE SHAPE OF THE N*(1440) AND N*(1535) */
/* RESONANCE ACCORDING */
/*     TO kitazoe's or J.D.JACKSON'S MASS FORMULA AND BREIT WIGNER */
/*     FORMULA FOR N* RESORANCE */
/*     DETERMINE DELTA MASS VIA REJECTION METHOD. */
	dmax__ = *srt - .9383f - .005f;
	dmin__ = 1.078f;
	if (n12 >= 13 && n12 <= 16) {
/* N*(1440) production */
	    if (dmax__ < 1.44f) {
		fm = fns_(&dmax__, srt, &c_b182);
	    } else {
/* lin-10/25/02 get rid of argument usage mismatch in FNS(): */
		xdmass = 1.44f;
/*          FM=FNS(1.44,SRT,1.) */
		fm = fns_(&xdmass, srt, &c_b183);
/* lin-10/25/02-end */
	    }
	    if (fm == 0.f) {
		fm = 1e-9f;
	    }
	    ntry2 = 0;
L11:
	    dm = ranart_(&rndf77_1.nseed) * (dmax__ - dmin__) + dmin__;
	    ++ntry2;
	    if (ranart_(&rndf77_1.nseed) > fns_(&dm, srt, &c_b183) / fm && 
		    ntry2 <= 10) {
		goto L11;
	    }
/* lin-2/26/03 limit the N* mass below a certain value */
/*     (here taken as its central value + 2* B-W fullwidth): */
	    if (dm > 2.14f) {
		goto L11;
	    }
	    goto L13;
	}
	if (n12 >= 17 && n12 <= 36) {
/* N*(1535) production */
	    if (dmax__ < 1.535f) {
		fm = fd5_(&dmax__, srt, &c_b182);
	    } else {
/* lin-10/25/02 get rid of argument usage mismatch in FNS(): */
		xdmass = 1.535f;
/*          FM=FD5(1.535,SRT,1.) */
		fm = fd5_(&xdmass, srt, &c_b183);
/* lin-10/25/02-end */
	    }
	    if (fm == 0.f) {
		fm = 1e-9f;
	    }
	    ntry1 = 0;
L12:
	    dm = ranart_(&rndf77_1.nseed) * (dmax__ - dmin__) + dmin__;
	    ++ntry1;
	    if (ranart_(&rndf77_1.nseed) > fd5_(&dm, srt, &c_b183) / fm && 
		    ntry1 <= 10) {
		goto L12;
	    }
/* lin-2/26/03 limit the N* mass below a certain value */
/*     (here taken as its central value + 2* B-W fullwidth): */
	    if (dm > 1.84f) {
		goto L12;
	    }
	}
L13:
/* ------------------------------------------------------- */
/* RELABLE BARYON I1 AND I2 */
/* 13 D(++)+D(-)--> N*(+)(14)+n */
	if (n12 == 13) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 11;
		cc_1.e[*i2 - 1] = dm;
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
	    } else {
		ee_1.lb[*i1 - 1] = 11;
		cc_1.e[*i1 - 1] = dm;
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
	    }
	    goto L200;
	}
/* 14 D(++)+D(-)--> N*(0)(14)+P */
	if (n12 == 14) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 10;
		cc_1.e[*i2 - 1] = dm;
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
	    } else {
		ee_1.lb[*i1 - 1] = 10;
		cc_1.e[*i1 - 1] = dm;
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .93828f;
	    }
	    goto L200;
	}
/* 15 D(+)+D(0)--> N*(+)(14)+n */
	if (n12 == 15) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 11;
		cc_1.e[*i2 - 1] = dm;
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
	    } else {
		ee_1.lb[*i1 - 1] = 11;
		cc_1.e[*i1 - 1] = dm;
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
	    }
	    goto L200;
	}
/* 16 D(+)+D(0)--> N*(0)(14)+P */
	if (n12 == 16) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 10;
		cc_1.e[*i2 - 1] = dm;
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
	    } else {
		ee_1.lb[*i1 - 1] = 10;
		cc_1.e[*i1 - 1] = dm;
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .93828f;
	    }
	    goto L200;
	}
/* 17 D(++)+D(0)--> N*(+)(14)+P */
	if (n12 == 17) {
	    ee_1.lb[*i2 - 1] = 13;
	    cc_1.e[*i2 - 1] = dm;
	    ee_1.lb[*i1 - 1] = 1;
	    cc_1.e[*i1 - 1] = .93828f;
	    goto L200;
	}
/* 18 D(++)+D(-)--> N*(0)(15)+P */
	if (n12 == 18) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 12;
		cc_1.e[*i2 - 1] = dm;
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
	    } else {
		ee_1.lb[*i1 - 1] = 12;
		cc_1.e[*i1 - 1] = dm;
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .93828f;
	    }
	    goto L200;
	}
/* 19 D(++)+D(-)--> N*(+)(15)+N */
	if (n12 == 19) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 13;
		cc_1.e[*i2 - 1] = dm;
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
	    } else {
		ee_1.lb[*i1 - 1] = 13;
		cc_1.e[*i1 - 1] = dm;
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
	    }
	    goto L200;
	}
/* 20 D(+)+D(+)--> N*(+)(15)+P */
	if (n12 == 20) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 13;
		cc_1.e[*i2 - 1] = dm;
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
	    } else {
		ee_1.lb[*i1 - 1] = 13;
		cc_1.e[*i1 - 1] = dm;
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .93828f;
	    }
	    goto L200;
	}
/* 21 D(+)+D(0)--> N*(+)(15)+N */
	if (n12 == 21) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 13;
		cc_1.e[*i2 - 1] = dm;
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
	    } else {
		ee_1.lb[*i1 - 1] = 13;
		cc_1.e[*i1 - 1] = dm;
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
	    }
	    goto L200;
	}
/* 22 D(+)+D(0)--> N*(0)(15)+P */
	if (n12 == 22) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 12;
		cc_1.e[*i2 - 1] = dm;
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
	    } else {
		ee_1.lb[*i1 - 1] = 12;
		cc_1.e[*i1 - 1] = dm;
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .93828f;
	    }
	    goto L200;
	}
/* 23 D(+)+D(-)--> N*(0)(15)+N */
	if (n12 == 23) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 12;
		cc_1.e[*i2 - 1] = dm;
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
	    } else {
		ee_1.lb[*i1 - 1] = 12;
		cc_1.e[*i1 - 1] = dm;
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
	    }
	    goto L200;
	}
/* 24 D(0)+D(0)--> N*(0)(15)+N */
	if (n12 == 24) {
	    ee_1.lb[*i2 - 1] = 12;
	    cc_1.e[*i2 - 1] = dm;
	    ee_1.lb[*i1 - 1] = 2;
	    cc_1.e[*i1 - 1] = .939457f;
	    goto L200;
	}
/* 25 N*(+)+N*(+)--> N*(0)(15)+P */
	if (n12 == 25) {
	    ee_1.lb[*i2 - 1] = 12;
	    cc_1.e[*i2 - 1] = dm;
	    ee_1.lb[*i1 - 1] = 1;
	    cc_1.e[*i1 - 1] = .93828f;
	    goto L200;
	}
/* 26 N*(0)+N*(0)--> N*(0)(15)+N */
	if (n12 == 26) {
	    ee_1.lb[*i2 - 1] = 12;
	    cc_1.e[*i2 - 1] = dm;
	    ee_1.lb[*i1 - 1] = 2;
	    cc_1.e[*i1 - 1] = .939457f;
	    goto L200;
	}
/* 27 N*(+)+N*(0)--> N*(+)(15)+N */
	if (n12 == 27) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 13;
		cc_1.e[*i2 - 1] = dm;
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
	    } else {
		ee_1.lb[*i1 - 1] = 13;
		cc_1.e[*i1 - 1] = dm;
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
	    }
	    goto L200;
	}
/* 28 N*(+)+N*(0)--> N*(0)(15)+P */
	if (n12 == 28) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 12;
		cc_1.e[*i2 - 1] = dm;
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
	    } else {
		ee_1.lb[*i1 - 1] = 12;
		cc_1.e[*i1 - 1] = dm;
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .93828f;
	    }
	    goto L200;
	}
/* 27 N*(+)+N*(0)--> N*(+)(15)+N */
	if (n12 == 27) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 13;
		cc_1.e[*i2 - 1] = dm;
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
	    } else {
		ee_1.lb[*i1 - 1] = 13;
		cc_1.e[*i1 - 1] = dm;
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
	    }
	    goto L200;
	}
/* 29 N*(+)+D(+)--> N*(+)(15)+P */
	if (n12 == 29) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 13;
		cc_1.e[*i2 - 1] = dm;
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
	    } else {
		ee_1.lb[*i1 - 1] = 13;
		cc_1.e[*i1 - 1] = dm;
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .93828f;
	    }
	    goto L200;
	}
/* 30 N*(+)+D(0)--> N*(+)(15)+N */
	if (n12 == 30) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 13;
		cc_1.e[*i2 - 1] = dm;
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
	    } else {
		ee_1.lb[*i1 - 1] = 13;
		cc_1.e[*i1 - 1] = dm;
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
	    }
	    goto L200;
	}
/* 31 N*(+)+D(-)--> N*(0)(15)+N */
	if (n12 == 31) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 12;
		cc_1.e[*i2 - 1] = dm;
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
	    } else {
		ee_1.lb[*i1 - 1] = 12;
		cc_1.e[*i1 - 1] = dm;
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
	    }
	    goto L200;
	}
/* 32 N*(0)+D(++)--> N*(+)(15)+P */
	if (n12 == 32) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 13;
		cc_1.e[*i2 - 1] = dm;
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
	    } else {
		ee_1.lb[*i1 - 1] = 13;
		cc_1.e[*i1 - 1] = dm;
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .93828f;
	    }
	    goto L200;
	}
/* 33 N*(0)+D(+)--> N*(+)(15)+N */
	if (n12 == 33) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 13;
		cc_1.e[*i2 - 1] = dm;
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
	    } else {
		ee_1.lb[*i1 - 1] = 13;
		cc_1.e[*i1 - 1] = dm;
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
	    }
	    goto L200;
	}
/* 34 N*(0)+D(+)--> N*(0)(15)+P */
	if (n12 == 34) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 12;
		cc_1.e[*i2 - 1] = dm;
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
	    } else {
		ee_1.lb[*i1 - 1] = 12;
		cc_1.e[*i1 - 1] = dm;
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .93828f;
	    }
	    goto L200;
	}
/* 35 N*(0)+D(0)--> N*(0)(15)+N */
	if (n12 == 35) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 12;
		cc_1.e[*i2 - 1] = dm;
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
	    } else {
		ee_1.lb[*i1 - 1] = 12;
		cc_1.e[*i1 - 1] = dm;
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
	    }
	    goto L200;
	}
/* 36 N*(+)+D(0)--> N*(0)(15)+P */
	if (n12 == 36) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 12;
		cc_1.e[*i2 - 1] = dm;
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
	    } else {
		ee_1.lb[*i1 - 1] = 12;
		cc_1.e[*i1 - 1] = dm;
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .93828f;
	    }
	    goto L200;
	}
L1012:
	*iblock = 55;
	leadng_1.lb1 = ee_1.lb[*i1 - 1];
	dpi_1.lb2 = ee_1.lb[*i2 - 1];
	ich = (i__1 = leadng_1.lb1 * dpi_1.lb2, abs(i__1));
/* ------------------------------------------------------- */
/* RELABLE BARYON I1 AND I2 in the reabsorption processes */
/* 37 D(++)+D(-)--> n+p */
	if (ich == 54) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .93828f;
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
	    } else {
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
	    }
	    goto L200;
	}
/* 38 D(+)+D(0)--> n+p */
	if (ich == 56) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .93828f;
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
	    } else {
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
	    }
	    goto L200;
	}
/* 39 D(++)+D(0)--> p+p */
	if (ich == 63) {
	    ee_1.lb[*i2 - 1] = 1;
	    cc_1.e[*i2 - 1] = .93828f;
	    ee_1.lb[*i1 - 1] = 1;
	    cc_1.e[*i1 - 1] = .93828f;
	    goto L200;
	}
/* 40 D(+)+D(+)--> p+p */
	if (ich == 64) {
	    ee_1.lb[*i2 - 1] = 1;
	    cc_1.e[*i2 - 1] = .93828f;
	    ee_1.lb[*i1 - 1] = 1;
	    cc_1.e[*i1 - 1] = .93828f;
	    goto L200;
	}
/* 41 D(+)+D(-)--> n+n */
	if (ich == 48) {
	    ee_1.lb[*i2 - 1] = 2;
	    cc_1.e[*i2 - 1] = .939457f;
	    ee_1.lb[*i1 - 1] = 2;
	    cc_1.e[*i1 - 1] = .939457f;
	    goto L200;
	}
/* 42 D(0)+D(0)--> n+n */
	if (ich == 36) {
	    ee_1.lb[*i2 - 1] = 2;
	    cc_1.e[*i2 - 1] = .939457f;
	    ee_1.lb[*i1 - 1] = 2;
	    cc_1.e[*i1 - 1] = .939457f;
	    goto L200;
	}
/* 43 N*(+)+N*(+)--> p+p */
	if (ich == 121 || ich == 169 || ich == 143) {
	    ee_1.lb[*i2 - 1] = 1;
	    cc_1.e[*i2 - 1] = .93828f;
	    ee_1.lb[*i1 - 1] = 1;
	    cc_1.e[*i1 - 1] = .93828f;
	    goto L200;
	}
/* 44 N*(0)(1440)+N*(0)--> n+n */
	if (ich == 100 || ich == 144 || ich == 120) {
	    ee_1.lb[*i2 - 1] = 2;
	    cc_1.e[*i2 - 1] = .939457f;
	    ee_1.lb[*i1 - 1] = 2;
	    cc_1.e[*i1 - 1] = .939457f;
	    goto L200;
	}
/* 45 N*(+)+N*(0)--> n+p */
	if (ich == 110 || ich == 156 || ich == 130 || ich == 132) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .93828f;
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
	    } else {
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
	    }
	    goto L200;
	}
/* 46 N*(+)+D(+)--> p+p */
	if (ich == 88 || ich == 104) {
	    ee_1.lb[*i2 - 1] = 1;
	    cc_1.e[*i2 - 1] = .93828f;
	    ee_1.lb[*i1 - 1] = 1;
	    cc_1.e[*i1 - 1] = .93828f;
	    goto L200;
	}
/* 47 N*(+)+D(0)--> n+p */
	if (ich == 77 || ich == 91) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .93828f;
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
	    } else {
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
	    }
	    goto L200;
	}
/* 48 N*(+)+D(-)--> n+n */
	if (ich == 66 || ich == 78) {
	    ee_1.lb[*i2 - 1] = 2;
	    cc_1.e[*i2 - 1] = .939457f;
	    ee_1.lb[*i1 - 1] = 2;
	    cc_1.e[*i1 - 1] = .939457f;
	    goto L200;
	}
/* 49 N*(0)+D(++)--> p+p */
	if (ich == 90 || ich == 108) {
	    ee_1.lb[*i2 - 1] = 1;
	    cc_1.e[*i2 - 1] = .93828f;
	    ee_1.lb[*i1 - 1] = 1;
	    cc_1.e[*i1 - 1] = .93828f;
	    goto L200;
	}
/* 50 N*(0)+D(0)--> n+n */
	if (ich == 70 || ich == 84) {
	    ee_1.lb[*i2 - 1] = 2;
	    cc_1.e[*i2 - 1] = .939457f;
	    ee_1.lb[*i1 - 1] = 2;
	    cc_1.e[*i1 - 1] = .939457f;
	    goto L200;
	}
/* 51 N*(0)+D(+)--> n+p */
	if (ich == 80 || ich == 96) {
	    if (ranart_(&rndf77_1.nseed) <= .5f) {
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
	    } else {
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .93828f;
	    }
	    goto L200;
	}
	ee_1.lb[*i1 - 1] = 1;
	cc_1.e[*i1 - 1] = .93828f;
	ee_1.lb[*i2 - 1] = 2;
	cc_1.e[*i2 - 1] = .939457f;
/* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH */
/* ENERGY CONSERVATION */
/* resonance production or absorption in resonance+resonance collisions is */
/* assumed to have the same pt distribution as pp */
L200:
	leadng_1.em1 = cc_1.e[*i1 - 1];
	dpi_1.em2 = cc_1.e[*i2 - 1];
/* Computing 2nd power */
	r__2 = *srt;
/* Computing 2nd power */
	r__3 = leadng_1.em1;
/* Computing 2nd power */
	r__4 = dpi_1.em2;
/* Computing 2nd power */
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
	r__5 = leadng_1.em1 * dpi_1.em2;
	pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
	if (pr2 <= 0.f) {
	    pr2 = 1e-9f;
	}
	pr = sqrt(pr2) / (*srt * 2.f);
	if (*srt <= 2.14f) {
	    c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	}
	if (*srt > 2.14f && *srt <= 2.4f) {
	    c1 = ang_(srt, &input1_1.iseed);
	}
	if (*srt > 2.4f) {
/* lin-10/25/02 get rid of argument usage mismatch in PTR(): */
	    xptr = pr * .33f;
/*         cc1=ptr(0.33*pr,iseed) */
	    cc1 = ptr_(&xptr, &input1_1.iseed);
/* lin-10/25/02-end */
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
	    r__1 = pr;
/* Computing 2nd power */
	    r__2 = cc1;
	    scheck = r__1 * r__1 - r__2 * r__2;
	    if (scheck < 0.f) {
		s_wsle(&io___811);
		do_lio(&c__9, &c__1, "scheck7: ", (ftnlen)9);
		do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
		e_wsle();
		scheck = 0.f;
	    }
	    c1 = sqrt(scheck) / pr;
/*         c1=sqrt(pr**2-cc1**2)/pr */
	}
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	if (ianti == 1 && ee_1.lb[*i1 - 1] >= 1 && ee_1.lb[*i2 - 1] >= 1) {
	    ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
	    ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
	}
    }
/* COM: SET THE NEW MOMENTUM COORDINATES */
/* lin-9/2012: check argument in sqrt(): */
L107:
/* Computing 2nd power */
    r__1 = c1;
    scheck = 1.f - r__1 * r__1;
    if (scheck < 0.f) {
	s_wsle(&io___812);
	do_lio(&c__9, &c__1, "scheck8: ", (ftnlen)9);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	scheck = 0.f;
    }
    s1 = sqrt(scheck);
/* 107   S1   = SQRT( 1.0 - C1**2 ) */
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__1 = c2;
    scheck = 1.f - r__1 * r__1;
    if (scheck < 0.f) {
	s_wsle(&io___814);
	do_lio(&c__9, &c__1, "scheck9: ", (ftnlen)9);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	scheck = 0.f;
    }
    s2 = sqrt(scheck);
/*      S2  =  SQRT( 1.0 - C2**2 ) */
    ct1 = cos(t1);
    st1 = sin(t1);
    ct2 = cos(t2);
    st2 = sin(t2);
    *pz = pr * (c1 * c2 - s1 * s2 * ct1);
    ss = c2 * s1 * ct1 + s2 * c1;
    *px = pr * (ss * ct2 - s1 * st1 * st2);
    *py = pr * (ss * st2 + s1 * st1 * ct2);
    return 0;
/* FOR THE DD-->KAON+X PROCESS, FIND MOMENTUM OF THE FINAL PARTICLES IN */
/* THE NUCLEUS-NUCLEUS CMS. */
L306:
/* sp11/21/01 phi production */
    if (xsk5 / sigk > ranart_(&rndf77_1.nseed)) {
	leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
	pz2 = bb_1.p[*i2 * 3 - 1];
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	++nn_1.nnn;
	pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 29;
	pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = 1.02f;
	*iblock = 222;
	goto L208;
    }
    *iblock = 10;
    if (ianti == 1) {
	*iblock = -10;
    }
    leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
    pz2 = bb_1.p[*i2 * 3 - 1];
/* DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE */
    ++nn_1.nnn;
    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 23;
    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .498f;
    if (*srt <= 2.63f) {
/* only lambda production is possible */
/* (1.1)P+P-->p+L+kaon+ */
	ic = 1;
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	ee_1.lb[*i2 - 1] = 14;
	goto L208;
    }
    if (*srt <= 2.74f && *srt > 2.63f) {
/* both Lambda and sigma production are possible */
	if (xsk1 / (xsk1 + xsk2) > ranart_(&rndf77_1.nseed)) {
/* lambda production */
	    ic = 1;
	    ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	    ee_1.lb[*i2 - 1] = 14;
	} else {
/* sigma production */
	    ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	    ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 15;
	    ic = 2;
	}
	goto L208;
    }
    if (*srt <= 2.77f && *srt > 2.74f) {
/* then pp-->Delta lamda kaon can happen */
	if (xsk1 / (xsk1 + xsk2 + xsk3) > ranart_(&rndf77_1.nseed)) {
/* * (1.1)P+P-->p+L+kaon+ */
	    ic = 1;
	    ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	    ee_1.lb[*i2 - 1] = 14;
	    goto L208;
	} else {
	    if (xsk2 / (xsk2 + xsk3) > ranart_(&rndf77_1.nseed)) {
/* pp-->psk */
		ic = 2;
		ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 
			1;
		ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 
			15;
	    } else {
/* pp-->D+l+k */
		ic = 3;
		ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 4) + 
			6;
		ee_1.lb[*i2 - 1] = 14;
	    }
	    goto L208;
	}
    }
    if (*srt > 2.77f) {
/* all four channels are possible */
	if (xsk1 / (xsk1 + xsk2 + xsk3 + xsk4) > ranart_(&rndf77_1.nseed)) {
/* p lambda k production */
	    ic = 1;
	    ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	    ee_1.lb[*i2 - 1] = 14;
	    goto L208;
	} else {
	    if (xsk3 / (xsk2 + xsk3 + xsk4) > ranart_(&rndf77_1.nseed)) {
/* delta l K production */
		ic = 3;
		ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 4) + 
			6;
		ee_1.lb[*i2 - 1] = 14;
		goto L208;
	    } else {
		if (xsk2 / (xsk2 + xsk4) > ranart_(&rndf77_1.nseed)) {
/* n sigma k production */
		    ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 
			    2) + 1;
		    ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 
			    3) + 15;
		    ic = 2;
		} else {
/* D sigma K */
		    ic = 4;
		    ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 
			    4) + 6;
		    ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 
			    3) + 15;
		}
		goto L208;
	    }
	}
    }
L208:
    if (ianti == 1 && ee_1.lb[*i1 - 1] >= 1 && ee_1.lb[*i2 - 1] >= 1) {
	ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
	ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
	if (pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] == 23) {
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 21;
	}
    }
    lbi1 = ee_1.lb[*i1 - 1];
    lbi2 = ee_1.lb[*i2 - 1];
/* KEEP ALL COORDINATES OF PARTICLE 2 FOR POSSIBLE PHASE SPACE CHANGE */
    ntry1 = 0;
L129:
    bbkaon_(&ic, srt, &px3, &py3, &pz3, &dm3, &px4, &py4, &pz4, &dm4, &ppx, &
	    ppy, &ppz, &icou1);
    ++ntry1;
    if (icou1 < 0 && ntry1 <= 20) {
	goto L129;
    }
/*       if(icou1.lt.0)return */
/* ROTATE THE MOMENTA OF PARTICLES IN THE CMS OF P1+P2 */
    rotate_(px, py, pz, &px3, &py3, &pz3);
    rotate_(px, py, pz, &px4, &py4, &pz4);
    rotate_(px, py, pz, &ppx, &ppy, &ppz);
/* FIND THE MOMENTUM OF PARTICLES IN THE FINAL STATE IN THE NUCLEUS- */
/* NUCLEUS CMS. FRAME */
/* (1) for the necleon/delta */
/*             LORENTZ-TRANSFORMATION INTO LAB FRAME FOR DELTA1 */
/* Computing 2nd power */
    r__1 = dm3;
/* Computing 2nd power */
    r__2 = px3;
/* Computing 2nd power */
    r__3 = py3;
/* Computing 2nd power */
    r__4 = pz3;
    e1cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    p1beta = px3 * bg_1.betax + py3 * bg_1.betay + pz3 * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) + e1cm);
    pt1i1 = bg_1.betax * transf + px3;
    pt2i1 = bg_1.betay * transf + py3;
    pt3i1 = bg_1.betaz * transf + pz3;
    eti1 = dm3;
/* (2) for the lambda/sigma */
/* Computing 2nd power */
    r__1 = dm4;
/* Computing 2nd power */
    r__2 = px4;
/* Computing 2nd power */
    r__3 = py4;
/* Computing 2nd power */
    r__4 = pz4;
    e2cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    p2beta = px4 * bg_1.betax + py4 * bg_1.betay + pz4 * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * p2beta / (bg_1.gamma + 1.f) + e2cm);
    pt1i2 = bg_1.betax * transf + px4;
    pt2i2 = bg_1.betay * transf + py4;
    pt3i2 = bg_1.betaz * transf + pz4;
    eti2 = dm4;
/* GET the kaon'S MOMENTUM AND COORDINATES IN NUCLEUS-NUCLEUS CMS. FRAME */
/* Computing 2nd power */
    r__1 = ppx;
/* Computing 2nd power */
    r__2 = ppy;
/* Computing 2nd power */
    r__3 = ppz;
    epcm = sqrt(r__1 * r__1 + .248004f + r__2 * r__2 + r__3 * r__3);
    ppbeta = ppx * bg_1.betax + ppy * bg_1.betay + ppz * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * ppbeta / (bg_1.gamma + 1.f) + epcm);
    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = bg_1.betax * 
	    transf + ppx;
    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = bg_1.betay * 
	    transf + ppy;
    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = bg_1.betaz * 
	    transf + ppz;
/* lin-5/2008: */
    dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = dpert_1.dpertp[*i1 - 
	    1] * dpert_1.dpertp[*i2 - 1];
/* lin-5/2008: */
/* 2007        X01 = 1.0 - 2.0 * RANART(NSEED) */
/*            Y01 = 1.0 - 2.0 * RANART(NSEED) */
/*            Z01 = 1.0 - 2.0 * RANART(NSEED) */
/*        IF ((X01*X01+Y01*Y01+Z01*Z01) .GT. 1.0) GOTO 2007 */
/*                RPION(1,NNN,IRUN)=R(1,I1)+0.5*x01 */
/*                RPION(2,NNN,IRUN)=R(2,I1)+0.5*y01 */
/*                RPION(3,NNN,IRUN)=R(3,I1)+0.5*z01 */
    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = aa_1.r__[*i1 * 3 - 
	    3];
    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = aa_1.r__[*i1 * 3 - 
	    2];
    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = aa_1.r__[*i1 * 3 - 
	    1];

/* assign the nucleon/delta and lambda/sigma to i1 or i2 to keep the */
/* leadng particle behaviour */
/*              if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0)then */
    bb_1.p[*i1 * 3 - 3] = pt1i1;
    bb_1.p[*i1 * 3 - 2] = pt2i1;
    bb_1.p[*i1 * 3 - 1] = pt3i1;
    cc_1.e[*i1 - 1] = eti1;
    ee_1.lb[*i1 - 1] = lbi1;
    bb_1.p[*i2 * 3 - 3] = pt1i2;
    bb_1.p[*i2 * 3 - 2] = pt2i2;
    bb_1.p[*i2 * 3 - 1] = pt3i2;
    cc_1.e[*i2 - 1] = eti2;
    ee_1.lb[*i2 - 1] = lbi2;
    leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
    leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
    leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
    leadng_1.em1 = cc_1.e[*i1 - 1];
    ee_1.id[*i1 - 1] = 2;
    ee_1.id[*i2 - 1] = 2;
    id1 = ee_1.id[*i1 - 1];
    leadng_1.lb1 = ee_1.lb[*i1 - 1];
    dpi_1.lb2 = ee_1.lb[*i2 - 1];
    am1 = leadng_1.em1;
    am2 = dpi_1.em2;
/* Computing 2nd power */
    r__1 = leadng_1.em1;
/* Computing 2nd power */
    r__2 = leadng_1.px1;
/* Computing 2nd power */
    r__3 = leadng_1.py1;
/* Computing 2nd power */
    r__4 = leadng_1.pz1;
    leadng_1.e1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    return 0;
/* lin-6/2008 D+D->Deuteron+pi: */
/*     FIND MOMENTUM OF THE FINAL PARTICLES IN THE NUCLEUS-NUCLEUS CMS. */
L108:
    if (para8_1.idpert == 1 && *ipert1 == 1 && para8_1.npertd >= 1) {
/*     For idpert=1: we produce npertd pert deuterons: */
	ndloop = para8_1.npertd;
    } else if (para8_1.idpert == 2 && para8_1.npertd >= 1) {
/*     For idpert=2: we first save information for npertd pert deuterons; */
/*     at the last ndloop we create the regular deuteron+pi */
/*     and those pert deuterons: */
	ndloop = para8_1.npertd + 1;
    } else {
/*     Just create the regular deuteron+pi: */
	ndloop = 1;
    }

    dprob1 = sdprod / *sig / (real) para8_1.npertd;
    i__1 = ndloop;
    for (idloop = 1; idloop <= i__1; ++idloop) {
	bbdangle_(&pxd, &pyd, &pzd, nt, ipert1, &ianti, &idloop, &pfinal, &
		dprob1, &lbm);
	rotate_(px, py, pz, &pxd, &pyd, &pzd);
/*     LORENTZ-TRANSFORMATION OF THE MOMENTUM OF PARTICLES IN THE FINAL STATE */
/*     FROM THE NN CMS FRAME INTO THE GLOBAL CMS FRAME: */
/*     For the Deuteron: */
	xmass = 1.8756f;
/* Computing 2nd power */
	r__1 = xmass;
/* Computing 2nd power */
	r__2 = pxd;
/* Computing 2nd power */
	r__3 = pyd;
/* Computing 2nd power */
	r__4 = pzd;
	e1dcm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p1dbeta = pxd * bg_1.betax + pyd * bg_1.betay + pzd * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p1dbeta / (bg_1.gamma + 1.f) + 
		e1dcm);
	pxi1 = bg_1.betax * transf + pxd;
	pyi1 = bg_1.betay * transf + pyd;
	pzi1 = bg_1.betaz * transf + pzd;
	if (ianti == 0) {
	    lbd = 42;
	} else {
	    lbd = -42;
	}
	if (para8_1.idpert == 1 && *ipert1 == 1 && para8_1.npertd >= 1) {
/* ccc  Perturbative production for idpert=1: */
	    ++nn_1.nnn;
	    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = pxi1;
	    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = pyi1;
	    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = pzi1;
	    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = 1.8756f;
	    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = lbd;
	    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = aa_1.r__[*
		    i1 * 3 - 3];
	    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = aa_1.r__[*
		    i1 * 3 - 2];
	    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = aa_1.r__[*
		    i1 * 3 - 1];
/* lin-6/2008 assign the perturbative probability: */
	    dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = sdprod / *
		    sig / (real) para8_1.npertd;
	} else if (para8_1.idpert == 2 && idloop <= para8_1.npertd) {
/* lin-6/2008 For idpert=2, we produce NPERTD perturbative (anti)deuterons */
/*     only when a regular (anti)deuteron+pi is produced in NN collisions. */
/*     First save the info for the perturbative deuterons: */
	    ppd[idloop * 3 - 3] = pxi1;
	    ppd[idloop * 3 - 2] = pyi1;
	    ppd[idloop * 3 - 1] = pzi1;
	    lbpd[idloop - 1] = lbd;
	} else {
/* ccc  Regular production: */
/*     For the regular pion: do LORENTZ-TRANSFORMATION: */
	    cc_1.e[*i1 - 1] = xmm;
/* Computing 2nd power */
	    r__1 = xmm;
/* Computing 2nd power */
	    r__2 = pxd;
/* Computing 2nd power */
	    r__3 = pyd;
/* Computing 2nd power */
	    r__4 = pzd;
	    e2picm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * 
		    r__4);
	    p2pibeta = -pxd * bg_1.betax - pyd * bg_1.betay - pzd * 
		    bg_1.betaz;
	    transf = bg_1.gamma * (bg_1.gamma * p2pibeta / (bg_1.gamma + 1.f) 
		    + e2picm);
	    pxi2 = bg_1.betax * transf - pxd;
	    pyi2 = bg_1.betay * transf - pyd;
	    pzi2 = bg_1.betaz * transf - pzd;
	    bb_1.p[*i1 * 3 - 3] = pxi2;
	    bb_1.p[*i1 * 3 - 2] = pyi2;
	    bb_1.p[*i1 * 3 - 1] = pzi2;
/*     Remove regular pion to check the equivalence */
/*     between the perturbative and regular deuteron results: */
/*                 E(i1)=0. */

	    ee_1.lb[*i1 - 1] = lbm;
	    leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
	    leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
	    leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
	    leadng_1.em1 = cc_1.e[*i1 - 1];
	    ee_1.id[*i1 - 1] = 2;
	    id1 = ee_1.id[*i1 - 1];
/* Computing 2nd power */
	    r__1 = leadng_1.em1;
/* Computing 2nd power */
	    r__2 = leadng_1.px1;
/* Computing 2nd power */
	    r__3 = leadng_1.py1;
/* Computing 2nd power */
	    r__4 = leadng_1.pz1;
	    leadng_1.e1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 
		    * r__4);
	    leadng_1.lb1 = ee_1.lb[*i1 - 1];
/*     For the regular deuteron: */
	    bb_1.p[*i2 * 3 - 3] = pxi1;
	    bb_1.p[*i2 * 3 - 2] = pyi1;
	    bb_1.p[*i2 * 3 - 1] = pzi1;
	    ee_1.lb[*i2 - 1] = lbd;
	    dpi_1.lb2 = ee_1.lb[*i2 - 1];
	    cc_1.e[*i2 - 1] = 1.8756f;
	    eti2 = cc_1.e[*i2 - 1];
	    ee_1.id[*i2 - 1] = 2;
/*     For idpert=2: create the perturbative deuterons: */
	    if (para8_1.idpert == 2 && idloop == ndloop) {
		i__2 = para8_1.npertd;
		for (ipertd = 1; ipertd <= i__2; ++ipertd) {
		    ++nn_1.nnn;
		    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = 
			    ppd[ipertd * 3 - 3];
		    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = 
			    ppd[ipertd * 3 - 2];
		    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = 
			    ppd[ipertd * 3 - 1];
		    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = 1.8756f;
		    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = lbpd[
			    ipertd - 1];
		    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = 
			    aa_1.r__[*i1 * 3 - 3];
		    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = 
			    aa_1.r__[*i1 * 3 - 2];
		    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = 
			    aa_1.r__[*i1 * 3 - 1];
/* lin-6/2008 assign the perturbative probability: */
		    dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = 1.f /
			     (real) para8_1.npertd;
		}
	    }
	}
    }
    *iblock = 501;
    return 0;
/* lin-6/2008 D+D->Deuteron+pi over */
} /* crdd_ */

/* ********************************* */
/* ********************************* */
/*                                                                      * */
/* Subroutine */ int init_(integer *minnum, integer *maxnum, integer *num, 
	real *radius, real *x0, real *z0, real *p0, real *gamma, integer *
	iseed, integer *mass, integer *iopt)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2, r__3;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double sqrt(doublereal), exp(doublereal), pow_dd(doublereal *, doublereal 
	    *);

    /* Local variables */
    static integer i__;
    static real x, y, z__, px, py, pz, beta;
    static integer idir;
    static real sign;
    static integer irun;
    static real ptot[3], rhow0, epart;
    static integer idnum, npart;
    static real rdist, rhows, scheck, pfermi;
    extern doublereal ranart_(integer *);

    /* Fortran I/O blocks */
    static cilist io___878 = { 0, 99, 0, 0, 0 };


/*                                                                      * */
/*       PURPOSE:     PROVIDING INITIAL CONDITIONS FOR PHASE-SPACE      * */
/*                    DISTRIBUTION OF TESTPARTICLES                     * */
/*       VARIABLES:   (ALL INPUT)                                       * */
/*         MINNUM  - FIRST TESTPARTICLE TREATED IN ONE RUN    (INTEGER) * */
/*         MAXNUM  - LAST TESTPARTICLE TREATED IN ONE RUN     (INTEGER) * */
/*         NUM     - NUMBER OF TESTPARTICLES PER NUCLEON      (INTEGER) * */
/*         RADIUS  - RADIUS OF NUCLEUS "FM"                      (REAL) * */
/*         X0,Z0   - DISPLACEMENT OF CENTER OF NUCLEUS IN X,Z-          * */
/*                   DIRECTION "FM"                              (REAL) * */
/*         P0      - MOMENTUM-BOOST IN C.M. FRAME "GEV/C"        (REAL) * */
/*         GAMMA   - RELATIVISTIC GAMMA-FACTOR                   (REAL) * */
/*         ISEED   - SEED FOR RANDOM-NUMBER GENERATOR         (INTEGER) * */
/*         MASS    - TOTAL MASS OF THE SYSTEM                 (INTEGER) * */
/*         IOPT    - OPTION FOR DIFFERENT OCCUPATION OF MOMENTUM        * */
/*                   SPACE                                    (INTEGER) * */
/*                                                                      * */
/* ********************************* */

/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /DD/ */
/* c      SAVE /EE/ */
/* c      SAVE /ss/ */
/* c      SAVE /RNDF77/ */
/* ---------------------------------------------------------------------- */
/*     PREPARATION FOR LORENTZ-TRANSFORMATIONS */

    if (*p0 != 0.f) {
	sign = *p0 / dabs(*p0);
    } else {
	sign = 0.f;
    }
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__1 = *gamma;
    scheck = r__1 * r__1 - 1.f;
    if (scheck < 0.f) {
	s_wsle(&io___878);
	do_lio(&c__9, &c__1, "scheck10: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	scheck = 0.f;
    }
    beta = sign * sqrt(scheck) / *gamma;
/*      BETA = SIGN * SQRT(GAMMA**2-1.)/GAMMA */
/* ----------------------------------------------------------------------- */
/*     TARGET-ID = 1 AND PROJECTILE-ID = -1 */

    if (*minnum == 1) {
	idnum = 1;
    } else {
	idnum = -1;
    }
/* ----------------------------------------------------------------------- */
/*     IDENTIFICATION OF TESTPARTICLES AND ASSIGMENT OF RESTMASS */

/*     LOOP OVER ALL PARALLEL RUNS: */
    i__1 = *num;
    for (irun = 1; irun <= i__1; ++irun) {
	i__2 = *maxnum + (irun - 1) * *mass;
	for (i__ = *minnum + (irun - 1) * *mass; i__ <= i__2; ++i__) {
	    ee_1.id[i__ - 1] = idnum;
	    cc_1.e[i__ - 1] = .9383f;
/* L100: */
	}
/* ----------------------------------------------------------------------- */
/*       OCCUPATION OF COORDINATE-SPACE */

	i__2 = *maxnum + (irun - 1) * *mass;
	for (i__ = *minnum + (irun - 1) * *mass; i__ <= i__2; ++i__) {
L200:
	    x = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	    y = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	    z__ = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	    if (x * x + y * y + z__ * z__ > 1.f) {
		goto L200;
	    }
	    aa_1.r__[i__ * 3 - 3] = x * *radius;
	    aa_1.r__[i__ * 3 - 2] = y * *radius;
	    aa_1.r__[i__ * 3 - 1] = z__ * *radius;
/* L300: */
	}
/* L400: */
    }
/* ======================================================================= */
    if (*iopt != 3) {
/* ----- */
/*     OPTION 1: USE WOODS-SAXON PARAMETRIZATION FOR DENSITY AND */
/* -----          CALCULATE LOCAL FERMI-MOMENTUM */

	rhow0 = .168f;
	i__1 = *num;
	for (irun = 1; irun <= i__1; ++irun) {
	    i__2 = *maxnum + (irun - 1) * *mass;
	    for (i__ = *minnum + (irun - 1) * *mass; i__ <= i__2; ++i__) {
L500:
		px = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
		py = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
		pz = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
		if (px * px + py * py + pz * pz > 1.f) {
		    goto L500;
		}
/* Computing 2nd power */
		r__1 = aa_1.r__[i__ * 3 - 3];
/* Computing 2nd power */
		r__2 = aa_1.r__[i__ * 3 - 2];
/* Computing 2nd power */
		r__3 = aa_1.r__[i__ * 3 - 1];
		rdist = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
		rhows = rhow0 / (exp((rdist - *radius) / .55f) + 1.f);
		d__1 = (doublereal) (rhows * 14.804406096562142f);
		pfermi = pow_dd(&d__1, &c_b5) * .197f;
/* ----- */
/*     OPTION 2: NUCLEAR MATTER CASE */
		if (*iopt == 2) {
		    pfermi = .27f;
		}
		if (*iopt == 4) {
		    pfermi = 0.f;
		}
/* ----- */
		bb_1.p[i__ * 3 - 3] = pfermi * px;
		bb_1.p[i__ * 3 - 2] = pfermi * py;
		bb_1.p[i__ * 3 - 1] = pfermi * pz;
/* L600: */
	    }

/*         SET TOTAL MOMENTUM TO 0 IN REST FRAME AND BOOST */

	    for (idir = 1; idir <= 3; ++idir) {
		ptot[idir - 1] = 0.f;
/* L700: */
	    }
	    npart = 0;
	    i__2 = *maxnum + (irun - 1) * *mass;
	    for (i__ = *minnum + (irun - 1) * *mass; i__ <= i__2; ++i__) {
		++npart;
		for (idir = 1; idir <= 3; ++idir) {
		    ptot[idir - 1] += bb_1.p[idir + i__ * 3 - 4];
/* L800: */
		}
/* L900: */
	    }
	    i__2 = *maxnum + (irun - 1) * *mass;
	    for (i__ = *minnum + (irun - 1) * *mass; i__ <= i__2; ++i__) {
		for (idir = 1; idir <= 3; ++idir) {
		    bb_1.p[idir + i__ * 3 - 4] -= ptot[idir - 1] / (real) 
			    npart;
/* L925: */
		}
/*           BOOST */
		if (*iopt == 1 || *iopt == 2) {
/* Computing 2nd power */
		    r__1 = bb_1.p[i__ * 3 - 3];
/* Computing 2nd power */
		    r__2 = bb_1.p[i__ * 3 - 2];
/* Computing 2nd power */
		    r__3 = bb_1.p[i__ * 3 - 1];
		    epart = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + 
			    .88040689000000005f);
		    bb_1.p[i__ * 3 - 1] = *gamma * (bb_1.p[i__ * 3 - 1] + 
			    beta * epart);
		} else {
		    bb_1.p[i__ * 3 - 1] += *p0;
		}
/* L950: */
	    }
/* L1000: */
	}
/* ----- */
    } else {
/* ----- */
/*     OPTION 3: GIVE ALL NUCLEONS JUST A Z-MOMENTUM ACCORDING TO */
/*               THE BOOST OF THE NUCLEI */

	i__1 = *num;
	for (irun = 1; irun <= i__1; ++irun) {
	    i__2 = *maxnum + (irun - 1) * *mass;
	    for (i__ = *minnum + (irun - 1) * *mass; i__ <= i__2; ++i__) {
		bb_1.p[i__ * 3 - 3] = 0.f;
		bb_1.p[i__ * 3 - 2] = 0.f;
		bb_1.p[i__ * 3 - 1] = *p0;
/* L1100: */
	    }
/* L1200: */
	}
/* ----- */
    }
/* ======================================================================= */
/*     PUT PARTICLES IN THEIR POSITION IN COORDINATE-SPACE */
/*     (SHIFT AND RELATIVISTIC CONTRACTION) */

    i__1 = *num;
    for (irun = 1; irun <= i__1; ++irun) {
	i__2 = *maxnum + (irun - 1) * *mass;
	for (i__ = *minnum + (irun - 1) * *mass; i__ <= i__2; ++i__) {
	    aa_1.r__[i__ * 3 - 3] += *x0;
/* two nuclei in touch after contraction */
	    aa_1.r__[i__ * 3 - 1] = (aa_1.r__[i__ * 3 - 1] + *z0) / *gamma;
/* two nuclei in touch before contraction */
/*          R(3,I) = R(3,I) / GAMMA + Z0 */
/* L1300: */
	}
/* L1400: */
    }

    return 0;
} /* init_ */

/* ********************************* */
/*                                                                      * */
/* Subroutine */ int dens_(integer *ipot, integer *mass, integer *num, 
	integer *nesc)
{
    /* Initialized data */

    static real zet[91] = { 1.f,0.f,0.f,0.f,0.f,1.f,0.f,0.f,0.f,0.f,0.f,0.f,
	    0.f,0.f,0.f,-1.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,
	    -1.f,0.f,1.f,0.f,-1.f,0.f,-1.f,0.f,-2.f,-1.f,0.f,1.f,0.f,0.f,0.f,
	    0.f,-1.f,0.f,1.f,0.f,-1.f,0.f,1.f,-1.f,0.f,1.f,2.f,0.f,1.f,0.f,
	    1.f,0.f,-1.f,0.f,1.f,0.f,0.f,0.f,-1.f,0.f,1.f,0.f,-1.f,0.f,1.f,
	    0.f,0.f,1.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,-1.f,0.f,0.f,0.f,
	    0.f,-1.f };

    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2, r__3, r__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer i_nint(real *);
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static real a, b;
    static integer i__, j;
    static real s, u;
    static integer ix, iy, iz;
    static real big, pxl[82369]	/* was [41][41][49] */, pyl[82369]	/* 
	    was [41][41][49] */, pzl[82369]	/* was [41][41][49] */, rho0, 
	    denr;
    static integer irun, msum;
    static real gamma, small, smass, smass2;

/*                                                                      * */
/*       PURPOSE:     CALCULATION OF LOCAL BARYON, MESON AND ENERGY     * */
/*                    DENSITY FROM SPATIAL DISTRIBUTION OF TESTPARTICLES* */
/*                                                                      * */
/*       VARIABLES (ALL INPUT, ALL INTEGER)                             * */
/*         MASS    -  MASS NUMBER OF THE SYSTEM                         * */
/*         NUM     -  NUMBER OF TESTPARTICLES PER NUCLEON               * */
/*                                                                      * */
/*         NESC    -  NUMBER OF ESCAPED PARTICLES      (INTEGER,OUTPUT) * */
/*                                                                      * */
/* ********************************* */

/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /DD/ */
/* c      SAVE /DDpi/ */
/* c      SAVE /EE/ */
/* c      SAVE /ss/ */
/* c      SAVE /RR/ */
/* c      SAVE /tt/ */

    for (iz = -24; iz <= 24; ++iz) {
	for (iy = -20; iy <= 20; ++iy) {
	    for (ix = -20; ix <= 20; ++ix) {
		dd_1.rho[ix + (iy + iz * 41) * 41 + 41184] = 0.f;
		dd_1.rhon[ix + (iy + iz * 41) * 41 + 41184] = 0.f;
		dd_1.rhop[ix + (iy + iz * 41) * 41 + 41184] = 0.f;
		ddpi_1.pirho[ix + (iy + iz * 41) * 41 + 41184] = 0.f;
		pxl[ix + (iy + iz * 41) * 41 + 41184] = 0.f;
		pyl[ix + (iy + iz * 41) * 41 + 41184] = 0.f;
		pzl[ix + (iy + iz * 41) * 41 + 41184] = 0.f;
		tt_1.pel[ix + (iy + iz * 41) * 41 + 41184] = 0.f;
		bbb_1.bxx[ix + (iy + iz * 41) * 41 + 41184] = 0.f;
		bbb_1.byy[ix + (iy + iz * 41) * 41 + 41184] = 0.f;
		bbb_1.bzz[ix + (iy + iz * 41) * 41 + 41184] = 0.f;
/* L100: */
	    }
/* L200: */
	}
/* L300: */
    }

    *nesc = 0;
    big = 1.f / ((real) (*num) * 3.f);
    small = 1.f / ((real) (*num) * 9.f);

    msum = 0;
    i__1 = *num;
    for (irun = 1; irun <= i__1; ++irun) {
	msum += rr_1.massr[irun - 1];
	i__2 = rr_1.massr[irun];
	for (j = 1; j <= i__2; ++j) {
	    i__ = j + msum;
	    ix = i_nint(&aa_1.r__[i__ * 3 - 3]);
	    iy = i_nint(&aa_1.r__[i__ * 3 - 2]);
	    iz = i_nint(&aa_1.r__[i__ * 3 - 1]);
	    if (ix <= -20 || ix >= 20 || iy <= -20 || iy >= 20 || iz <= -24 ||
		     iz >= 24) {
		++(*nesc);
	    } else {

/* sp01/04/02 include baryon density */
		if (j > *mass) {
		    goto L30;
		}
/*         if( (lb(i).eq.1.or.lb(i).eq.2) .or. */
/*    &    (lb(i).ge.6.and.lb(i).le.17) )then */
/* (1) baryon density */
		dd_1.rho[ix + (iy + iz * 41) * 41 + 41184] += big;
		dd_1.rho[ix + 1 + (iy + iz * 41) * 41 + 41184] += small;
		dd_1.rho[ix - 1 + (iy + iz * 41) * 41 + 41184] += small;
		dd_1.rho[ix + (iy + 1 + iz * 41) * 41 + 41184] += small;
		dd_1.rho[ix + (iy - 1 + iz * 41) * 41 + 41184] += small;
		dd_1.rho[ix + (iy + (iz + 1) * 41) * 41 + 41184] += small;
		dd_1.rho[ix + (iy + (iz - 1) * 41) * 41 + 41184] += small;
/* (2) CALCULATE THE PROTON DENSITY */
		if (zet[ee_1.lb[i__ - 1] + 45] != 0.f) {
		    dd_1.rhop[ix + (iy + iz * 41) * 41 + 41184] += big;
		    dd_1.rhop[ix + 1 + (iy + iz * 41) * 41 + 41184] += small;
		    dd_1.rhop[ix - 1 + (iy + iz * 41) * 41 + 41184] += small;
		    dd_1.rhop[ix + (iy + 1 + iz * 41) * 41 + 41184] += small;
		    dd_1.rhop[ix + (iy - 1 + iz * 41) * 41 + 41184] += small;
		    dd_1.rhop[ix + (iy + (iz + 1) * 41) * 41 + 41184] += 
			    small;
		    dd_1.rhop[ix + (iy + (iz - 1) * 41) * 41 + 41184] += 
			    small;
		    goto L40;
		}
/* (3) CALCULATE THE NEUTRON DENSITY */
		if (zet[ee_1.lb[i__ - 1] + 45] == 0.f) {
		    dd_1.rhon[ix + (iy + iz * 41) * 41 + 41184] += big;
		    dd_1.rhon[ix + 1 + (iy + iz * 41) * 41 + 41184] += small;
		    dd_1.rhon[ix - 1 + (iy + iz * 41) * 41 + 41184] += small;
		    dd_1.rhon[ix + (iy + 1 + iz * 41) * 41 + 41184] += small;
		    dd_1.rhon[ix + (iy - 1 + iz * 41) * 41 + 41184] += small;
		    dd_1.rhon[ix + (iy + (iz + 1) * 41) * 41 + 41184] += 
			    small;
		    dd_1.rhon[ix + (iy + (iz - 1) * 41) * 41 + 41184] += 
			    small;
		    goto L40;
		}
/*           else    !! sp01/04/02 */
/* (4) meson density */
L30:
		ddpi_1.pirho[ix + (iy + iz * 41) * 41 + 41184] += big;
		ddpi_1.pirho[ix + 1 + (iy + iz * 41) * 41 + 41184] += small;
		ddpi_1.pirho[ix - 1 + (iy + iz * 41) * 41 + 41184] += small;
		ddpi_1.pirho[ix + (iy + 1 + iz * 41) * 41 + 41184] += small;
		ddpi_1.pirho[ix + (iy - 1 + iz * 41) * 41 + 41184] += small;
		ddpi_1.pirho[ix + (iy + (iz + 1) * 41) * 41 + 41184] += small;
		ddpi_1.pirho[ix + (iy + (iz - 1) * 41) * 41 + 41184] += small;
/*           endif    !! sp01/04/02 */
/* to calculate the Gamma factor in each cell */
/* (1) PX */
L40:
		pxl[ix + (iy + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 - 3] *
			 big;
		pxl[ix + 1 + (iy + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 - 
			3] * small;
		pxl[ix - 1 + (iy + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 - 
			3] * small;
		pxl[ix + (iy + 1 + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 - 
			3] * small;
		pxl[ix + (iy - 1 + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 - 
			3] * small;
		pxl[ix + (iy + (iz + 1) * 41) * 41 + 41184] += bb_1.p[i__ * 3 
			- 3] * small;
		pxl[ix + (iy + (iz - 1) * 41) * 41 + 41184] += bb_1.p[i__ * 3 
			- 3] * small;
/* (2) PY */
		pyl[ix + (iy + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 - 2] *
			 big;
		pyl[ix + 1 + (iy + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 - 
			2] * small;
		pyl[ix - 1 + (iy + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 - 
			2] * small;
		pyl[ix + (iy + 1 + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 - 
			2] * small;
		pyl[ix + (iy - 1 + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 - 
			2] * small;
		pyl[ix + (iy + (iz + 1) * 41) * 41 + 41184] += bb_1.p[i__ * 3 
			- 2] * small;
		pyl[ix + (iy + (iz - 1) * 41) * 41 + 41184] += bb_1.p[i__ * 3 
			- 2] * small;
/* (3) PZ */
		pzl[ix + (iy + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 - 1] *
			 big;
		pzl[ix + 1 + (iy + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 - 
			1] * small;
		pzl[ix - 1 + (iy + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 - 
			1] * small;
		pzl[ix + (iy + 1 + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 - 
			1] * small;
		pzl[ix + (iy - 1 + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 - 
			1] * small;
		pzl[ix + (iy + (iz + 1) * 41) * 41 + 41184] += bb_1.p[i__ * 3 
			- 1] * small;
		pzl[ix + (iy + (iz - 1) * 41) * 41 + 41184] += bb_1.p[i__ * 3 
			- 1] * small;
/* (4) ENERGY */
/* Computing 2nd power */
		r__1 = cc_1.e[i__ - 1];
/* Computing 2nd power */
		r__2 = bb_1.p[i__ * 3 - 3];
/* Computing 2nd power */
		r__3 = bb_1.p[i__ * 3 - 2];
/* Computing 2nd power */
		r__4 = bb_1.p[i__ * 3 - 1];
		tt_1.pel[ix + (iy + iz * 41) * 41 + 41184] += sqrt(r__1 * 
			r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4) * big;
/* Computing 2nd power */
		r__1 = cc_1.e[i__ - 1];
/* Computing 2nd power */
		r__2 = bb_1.p[i__ * 3 - 3];
/* Computing 2nd power */
		r__3 = bb_1.p[i__ * 3 - 2];
/* Computing 2nd power */
		r__4 = bb_1.p[i__ * 3 - 1];
		tt_1.pel[ix + 1 + (iy + iz * 41) * 41 + 41184] += sqrt(r__1 * 
			r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4) * 
			small;
/* Computing 2nd power */
		r__1 = cc_1.e[i__ - 1];
/* Computing 2nd power */
		r__2 = bb_1.p[i__ * 3 - 3];
/* Computing 2nd power */
		r__3 = bb_1.p[i__ * 3 - 2];
/* Computing 2nd power */
		r__4 = bb_1.p[i__ * 3 - 1];
		tt_1.pel[ix - 1 + (iy + iz * 41) * 41 + 41184] += sqrt(r__1 * 
			r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4) * 
			small;
/* Computing 2nd power */
		r__1 = cc_1.e[i__ - 1];
/* Computing 2nd power */
		r__2 = bb_1.p[i__ * 3 - 3];
/* Computing 2nd power */
		r__3 = bb_1.p[i__ * 3 - 2];
/* Computing 2nd power */
		r__4 = bb_1.p[i__ * 3 - 1];
		tt_1.pel[ix + (iy + 1 + iz * 41) * 41 + 41184] += sqrt(r__1 * 
			r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4) * 
			small;
/* Computing 2nd power */
		r__1 = cc_1.e[i__ - 1];
/* Computing 2nd power */
		r__2 = bb_1.p[i__ * 3 - 3];
/* Computing 2nd power */
		r__3 = bb_1.p[i__ * 3 - 2];
/* Computing 2nd power */
		r__4 = bb_1.p[i__ * 3 - 1];
		tt_1.pel[ix + (iy - 1 + iz * 41) * 41 + 41184] += sqrt(r__1 * 
			r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4) * 
			small;
/* Computing 2nd power */
		r__1 = cc_1.e[i__ - 1];
/* Computing 2nd power */
		r__2 = bb_1.p[i__ * 3 - 3];
/* Computing 2nd power */
		r__3 = bb_1.p[i__ * 3 - 2];
/* Computing 2nd power */
		r__4 = bb_1.p[i__ * 3 - 1];
		tt_1.pel[ix + (iy + (iz + 1) * 41) * 41 + 41184] += sqrt(r__1 
			* r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4) * 
			small;
/* Computing 2nd power */
		r__1 = cc_1.e[i__ - 1];
/* Computing 2nd power */
		r__2 = bb_1.p[i__ * 3 - 3];
/* Computing 2nd power */
		r__3 = bb_1.p[i__ * 3 - 2];
/* Computing 2nd power */
		r__4 = bb_1.p[i__ * 3 - 1];
		tt_1.pel[ix + (iy + (iz - 1) * 41) * 41 + 41184] += sqrt(r__1 
			* r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4) * 
			small;
	    }
/* L400: */
	}
    }

    for (iz = -24; iz <= 24; ++iz) {
	for (iy = -20; iy <= 20; ++iy) {
	    for (ix = -20; ix <= 20; ++ix) {
		if (dd_1.rho[ix + (iy + iz * 41) * 41 + 41184] == 0.f || 
			tt_1.pel[ix + (iy + iz * 41) * 41 + 41184] == 0.f) {
		    goto L101;
		}
/* Computing 2nd power */
		r__1 = tt_1.pel[ix + (iy + iz * 41) * 41 + 41184];
/* Computing 2nd power */
		r__2 = pxl[ix + (iy + iz * 41) * 41 + 41184];
/* Computing 2nd power */
		r__3 = pyl[ix + (iy + iz * 41) * 41 + 41184];
/* Computing 2nd power */
		r__4 = pzl[ix + (iy + iz * 41) * 41 + 41184];
		smass2 = r__1 * r__1 - r__2 * r__2 - r__3 * r__3 - r__4 * 
			r__4;
		if (smass2 <= 0.f) {
		    smass2 = 1e-6f;
		}
		smass = sqrt(smass2);
		if (smass == 0.f) {
		    smass = 1e-6f;
		}
		gamma = tt_1.pel[ix + (iy + iz * 41) * 41 + 41184] / smass;
		if (gamma == 0.f) {
		    goto L101;
		}
		bbb_1.bxx[ix + (iy + iz * 41) * 41 + 41184] = pxl[ix + (iy + 
			iz * 41) * 41 + 41184] / tt_1.pel[ix + (iy + iz * 41) 
			* 41 + 41184];
		bbb_1.byy[ix + (iy + iz * 41) * 41 + 41184] = pyl[ix + (iy + 
			iz * 41) * 41 + 41184] / tt_1.pel[ix + (iy + iz * 41) 
			* 41 + 41184];
		bbb_1.bzz[ix + (iy + iz * 41) * 41 + 41184] = pzl[ix + (iy + 
			iz * 41) * 41 + 41184] / tt_1.pel[ix + (iy + iz * 41) 
			* 41 + 41184];
		dd_1.rho[ix + (iy + iz * 41) * 41 + 41184] /= gamma;
		dd_1.rhon[ix + (iy + iz * 41) * 41 + 41184] /= gamma;
		dd_1.rhop[ix + (iy + iz * 41) * 41 + 41184] /= gamma;
		ddpi_1.pirho[ix + (iy + iz * 41) * 41 + 41184] /= gamma;
/* Computing 2nd power */
		r__1 = gamma;
		tt_1.pel[ix + (iy + iz * 41) * 41 + 41184] /= r__1 * r__1;
		rho0 = .163f;
		if (*ipot == 0) {
		    u = 0.f;
		    goto L70;
		}
		if (*ipot == 1 || *ipot == 6) {
		    a = -.1236f;
		    b = .0704f;
		    s = 2.f;
		    goto L60;
		}
		if (*ipot == 2 || *ipot == 7) {
		    a = -.218f;
		    b = .164f;
		    s = 1.3333333333333333f;
		}
		if (*ipot == 3) {
		    a = -.3581f;
		    b = .3048f;
		    s = 1.167f;
		    goto L60;
		}
		if (*ipot == 4) {
		    denr = dd_1.rho[ix + (iy + iz * 41) * 41 + 41184] / rho0;
		    b = .3048f;
		    s = 1.167f;
		    if (denr <= 4.f || denr > 7.f) {
			a = -.3581f;
		    } else {
			d__1 = (doublereal) denr;
			d__2 = (doublereal) denr;
			a = -b * pow_dd(&d__1, &c_b397) - pow_dd(&d__2, &
				c_b398) * .023999999999999997f;
		    }
		    goto L60;
		}
L60:
/* Computing 2nd power */
		r__1 = dd_1.rho[ix + (iy + iz * 41) * 41 + 41184];
		d__1 = (doublereal) (dd_1.rho[ix + (iy + iz * 41) * 41 + 
			41184] / rho0);
		d__2 = (doublereal) s;
		u = a * .5f * (r__1 * r__1) / rho0 + b / (s + 1) * pow_dd(&
			d__1, &d__2) * dd_1.rho[ix + (iy + iz * 41) * 41 + 
			41184];
L70:
		tt_1.pel[ix + (iy + iz * 41) * 41 + 41184] += u;
L101:
		;
	    }
/* L201: */
	}
/* L301: */
    }
    return 0;
} /* dens_ */

/* ********************************* */
/*                                                                      * */
/* Subroutine */ int gradu_(integer *iopt, integer *ix, integer *iy, integer *
	iz, real *gradx, real *grady, real *gradz)
{
    /* System generated locals */
    real r__1, r__2;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static real ef, eh, cf0, den0, ene0, denr, eqgp, acoef, expnt, acoef1, 
	    acoef2, expnt2, rxmins, rymins, rzmins, rxplus, ryplus, rzplus;

/*                                                                      * */
/*       PURPOSE:     DETERMINE GRAD(U(RHO(X,Y,Z)))                     * */
/*       VARIABLES:                                                     * */
/*         IOPT                - METHOD FOR EVALUATING THE GRADIENT     * */
/*                                                      (INTEGER,INPUT) * */
/*         IX, IY, IZ          - COORDINATES OF POINT   (INTEGER,INPUT) * */
/*         GRADX, GRADY, GRADZ - GRADIENT OF U            (REAL,OUTPUT) * */
/*                                                                      * */
/* ********************************* */

/* c      SAVE /DD/ */
/* c      SAVE /ss/ */
/* c      SAVE /tt/ */

    rxplus = dd_1.rho[*ix + 1 + (*iy + *iz * 41) * 41 + 41184] / .167f;
    rxmins = dd_1.rho[*ix - 1 + (*iy + *iz * 41) * 41 + 41184] / .167f;
    ryplus = dd_1.rho[*ix + (*iy + 1 + *iz * 41) * 41 + 41184] / .167f;
    rymins = dd_1.rho[*ix + (*iy - 1 + *iz * 41) * 41 + 41184] / .167f;
    rzplus = dd_1.rho[*ix + (*iy + (*iz + 1) * 41) * 41 + 41184] / .167f;
    rzmins = dd_1.rho[*ix + (*iy + (*iz - 1) * 41) * 41 + 41184] / .167f;
    den0 = dd_1.rho[*ix + (*iy + *iz * 41) * 41 + 41184] / .167f;
    ene0 = tt_1.pel[*ix + (*iy + *iz * 41) * 41 + 41184];
/* ----------------------------------------------------------------------- */
    switch (*iopt) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L3;
	case 4:  goto L4;
	case 5:  goto L5;
    }
    if (*iopt == 6) {
	goto L6;
    }
    if (*iopt == 7) {
	goto L7;
    }

L1:
/*       POTENTIAL USED IN 1) (STIFF): */
/*       U = -.124 * RHO/RHO0 + .0705 (RHO/RHO0)**2 GEV */

/* Computing 2nd power */
    r__1 = rxplus;
/* Computing 2nd power */
    r__2 = rxmins;
    *gradx = (rxplus - rxmins) * -.062f + (r__1 * r__1 - r__2 * r__2) * 
	    .03525f;
/* Computing 2nd power */
    r__1 = ryplus;
/* Computing 2nd power */
    r__2 = rymins;
    *grady = (ryplus - rymins) * -.062f + (r__1 * r__1 - r__2 * r__2) * 
	    .03525f;
/* Computing 2nd power */
    r__1 = rzplus;
/* Computing 2nd power */
    r__2 = rzmins;
    *gradz = (rzplus - rzmins) * -.062f + (r__1 * r__1 - r__2 * r__2) * 
	    .03525f;
    return 0;

L2:
/*       POTENTIAL USED IN 2): */
/*       U = -.218 * RHO/RHO0 + .164 (RHO/RHO0)**(4/3) GEV */

    expnt = 1.3333333f;
    d__1 = (doublereal) rxplus;
    d__2 = (doublereal) expnt;
    d__3 = (doublereal) rxmins;
    d__4 = (doublereal) expnt;
    *gradx = (rxplus - rxmins) * -.109f + (pow_dd(&d__1, &d__2) - pow_dd(&
	    d__3, &d__4)) * .082f;
    d__1 = (doublereal) ryplus;
    d__2 = (doublereal) expnt;
    d__3 = (doublereal) rymins;
    d__4 = (doublereal) expnt;
    *grady = (ryplus - rymins) * -.109f + (pow_dd(&d__1, &d__2) - pow_dd(&
	    d__3, &d__4)) * .082f;
    d__1 = (doublereal) rzplus;
    d__2 = (doublereal) expnt;
    d__3 = (doublereal) rzmins;
    d__4 = (doublereal) expnt;
    *gradz = (rzplus - rzmins) * -.109f + (pow_dd(&d__1, &d__2) - pow_dd(&
	    d__3, &d__4)) * .082f;
    return 0;

L3:
/*       POTENTIAL USED IN 3) (SOFT): */
/*       U = -.356 * RHO/RHO0 + .303 * (RHO/RHO0)**(7/6)  GEV */

    expnt = 1.1666667f;
    acoef = .178f;
    d__1 = (doublereal) rxplus;
    d__2 = (doublereal) expnt;
    d__3 = (doublereal) rxmins;
    d__4 = (doublereal) expnt;
    *gradx = -acoef * (rxplus - rxmins) + (pow_dd(&d__1, &d__2) - pow_dd(&
	    d__3, &d__4)) * .1515f;
    d__1 = (doublereal) ryplus;
    d__2 = (doublereal) expnt;
    d__3 = (doublereal) rymins;
    d__4 = (doublereal) expnt;
    *grady = -acoef * (ryplus - rymins) + (pow_dd(&d__1, &d__2) - pow_dd(&
	    d__3, &d__4)) * .1515f;
    d__1 = (doublereal) rzplus;
    d__2 = (doublereal) expnt;
    d__3 = (doublereal) rzmins;
    d__4 = (doublereal) expnt;
    *gradz = -acoef * (rzplus - rzmins) + (pow_dd(&d__1, &d__2) - pow_dd(&
	    d__3, &d__4)) * .1515f;
    return 0;


L4:
/*       POTENTIAL USED IN 4) (super-soft in the mixed phase of 4 < rho/rho <7): */
/*       U1 = -.356 * RHO/RHO0 + .303 * (RHO/RHO0)**(7/6)  GEV */
/*       normal phase, soft eos of iopt=3 */
/*       U2 = -.02 * (RHO/RHO0)**(2/3) -0.0253 * (RHO/RHO0)**(7/6)  GEV */

    eh = 4.f;
    eqgp = 7.f;
    acoef = .178f;
    expnt = 1.1666667f;
    denr = dd_1.rho[*ix + (*iy + *iz * 41) * 41 + 41184] / .167f;
    if (denr <= eh || denr >= eqgp) {
	d__1 = (doublereal) rxplus;
	d__2 = (doublereal) expnt;
	d__3 = (doublereal) rxmins;
	d__4 = (doublereal) expnt;
	*gradx = -acoef * (rxplus - rxmins) + (pow_dd(&d__1, &d__2) - pow_dd(&
		d__3, &d__4)) * .1515f;
	d__1 = (doublereal) ryplus;
	d__2 = (doublereal) expnt;
	d__3 = (doublereal) rymins;
	d__4 = (doublereal) expnt;
	*grady = -acoef * (ryplus - rymins) + (pow_dd(&d__1, &d__2) - pow_dd(&
		d__3, &d__4)) * .1515f;
	d__1 = (doublereal) rzplus;
	d__2 = (doublereal) expnt;
	d__3 = (doublereal) rzmins;
	d__4 = (doublereal) expnt;
	*gradz = -acoef * (rzplus - rzmins) + (pow_dd(&d__1, &d__2) - pow_dd(&
		d__3, &d__4)) * .1515f;
    } else {
	acoef1 = .178f;
	acoef2 = 0.f;
	expnt2 = .66666666666666663f;
	d__1 = (doublereal) rxplus;
	d__2 = (doublereal) expnt;
	d__3 = (doublereal) rxmins;
	d__4 = (doublereal) expnt;
	d__5 = (doublereal) rxplus;
	d__6 = (doublereal) expnt2;
	d__7 = (doublereal) rxmins;
	d__8 = (doublereal) expnt2;
	*gradx = -acoef1 * (pow_dd(&d__1, &d__2) - pow_dd(&d__3, &d__4)) - 
		acoef2 * (pow_dd(&d__5, &d__6) - pow_dd(&d__7, &d__8));
	d__1 = (doublereal) ryplus;
	d__2 = (doublereal) expnt;
	d__3 = (doublereal) rymins;
	d__4 = (doublereal) expnt;
	d__5 = (doublereal) ryplus;
	d__6 = (doublereal) expnt2;
	d__7 = (doublereal) rymins;
	d__8 = (doublereal) expnt2;
	*grady = -acoef1 * (pow_dd(&d__1, &d__2) - pow_dd(&d__3, &d__4)) - 
		acoef2 * (pow_dd(&d__5, &d__6) - pow_dd(&d__7, &d__8));
	d__1 = (doublereal) rzplus;
	d__2 = (doublereal) expnt;
	d__3 = (doublereal) rzmins;
	d__4 = (doublereal) expnt;
	d__5 = (doublereal) rzplus;
	d__6 = (doublereal) expnt2;
	d__7 = (doublereal) rzmins;
	d__8 = (doublereal) expnt2;
	*gradz = -acoef1 * (pow_dd(&d__1, &d__2) - pow_dd(&d__3, &d__4)) - 
		acoef2 * (pow_dd(&d__5, &d__6) - pow_dd(&d__7, &d__8));
    }
    return 0;

L5:
/*       POTENTIAL USED IN 5) (SUPER STIFF): */
/*       U = -.10322 * RHO/RHO0 + .04956 * (RHO/RHO0)**(2.77)  GEV */

    expnt = 2.77f;
    d__1 = (doublereal) rxplus;
    d__2 = (doublereal) expnt;
    d__3 = (doublereal) rxmins;
    d__4 = (doublereal) expnt;
    *gradx = (rxplus - rxmins) * -.0516f + (pow_dd(&d__1, &d__2) - pow_dd(&
	    d__3, &d__4)) * .02498f;
    d__1 = (doublereal) ryplus;
    d__2 = (doublereal) expnt;
    d__3 = (doublereal) rymins;
    d__4 = (doublereal) expnt;
    *grady = (ryplus - rymins) * -.0516f + (pow_dd(&d__1, &d__2) - pow_dd(&
	    d__3, &d__4)) * .02498f;
    d__1 = (doublereal) rzplus;
    d__2 = (doublereal) expnt;
    d__3 = (doublereal) rzmins;
    d__4 = (doublereal) expnt;
    *gradz = (rzplus - rzmins) * -.0516f + (pow_dd(&d__1, &d__2) - pow_dd(&
	    d__3, &d__4)) * .02498f;
    return 0;

L6:
/*       POTENTIAL USED IN 6) (STIFF-qgp): */
/*       U = -.124 * RHO/RHO0 + .0705 (RHO/RHO0)**2 GEV */

    if (ene0 <= .5f) {
/* Computing 2nd power */
	r__1 = rxplus;
/* Computing 2nd power */
	r__2 = rxmins;
	*gradx = (rxplus - rxmins) * -.062f + (r__1 * r__1 - r__2 * r__2) * 
		.03525f;
/* Computing 2nd power */
	r__1 = ryplus;
/* Computing 2nd power */
	r__2 = rymins;
	*grady = (ryplus - rymins) * -.062f + (r__1 * r__1 - r__2 * r__2) * 
		.03525f;
/* Computing 2nd power */
	r__1 = rzplus;
/* Computing 2nd power */
	r__2 = rzmins;
	*gradz = (rzplus - rzmins) * -.062f + (r__1 * r__1 - r__2 * r__2) * 
		.03525f;
	return 0;
    }
    if (ene0 > .5f && ene0 <= 1.5f) {
/*       U=c1-ef*rho/rho0**2/3 */
	ef = .035999999999999997f;
	d__1 = (doublereal) rxplus;
	d__2 = (doublereal) rxmins;
	*gradx = ef * -.5f * (pow_dd(&d__1, &c_b407) - pow_dd(&d__2, &c_b407))
		;
	d__1 = (doublereal) ryplus;
	d__2 = (doublereal) rymins;
	*grady = ef * -.5f * (pow_dd(&d__1, &c_b407) - pow_dd(&d__2, &c_b407))
		;
	d__1 = (doublereal) rzplus;
	d__2 = (doublereal) rzmins;
	*gradz = ef * -.5f * (pow_dd(&d__1, &c_b407) - pow_dd(&d__2, &c_b407))
		;
	return 0;
    }
    if (ene0 > 1.5f) {
/* U=800*(rho/rho0)**1/3.-Ef*(rho/rho0)**2/3.-c2 */
	ef = .035999999999999997f;
	cf0 = .8f;
	d__1 = (doublereal) rxplus;
	d__2 = (doublereal) rxmins;
	d__3 = (doublereal) rxplus;
	d__4 = (doublereal) rxmins;
	*gradx = cf0 * .5f * (pow_dd(&d__1, &c_b413) - pow_dd(&d__2, &c_b413))
		 - ef * .5f * (pow_dd(&d__3, &c_b407) - pow_dd(&d__4, &c_b407)
		);
	d__1 = (doublereal) ryplus;
	d__2 = (doublereal) rymins;
	d__3 = (doublereal) ryplus;
	d__4 = (doublereal) rymins;
	*grady = cf0 * .5f * (pow_dd(&d__1, &c_b413) - pow_dd(&d__2, &c_b413))
		 - ef * .5f * (pow_dd(&d__3, &c_b407) - pow_dd(&d__4, &c_b407)
		);
	d__1 = (doublereal) rzplus;
	d__2 = (doublereal) rzmins;
	d__3 = (doublereal) rzplus;
	d__4 = (doublereal) rzmins;
	*gradz = cf0 * .5f * (pow_dd(&d__1, &c_b413) - pow_dd(&d__2, &c_b413))
		 - ef * .5f * (pow_dd(&d__3, &c_b407) - pow_dd(&d__4, &c_b407)
		);
	return 0;
    }

L7:
/*       POTENTIAL USED IN 7) (Soft-qgp): */
    if (den0 <= 4.5f) {
/*       POTENTIAL USED is the same as IN 3) (SOFT): */
/*       U = -.356 * RHO/RHO0 + .303 * (RHO/RHO0)**(7/6)  GEV */

	expnt = 1.1666667f;
	acoef = .178f;
	d__1 = (doublereal) rxplus;
	d__2 = (doublereal) expnt;
	d__3 = (doublereal) rxmins;
	d__4 = (doublereal) expnt;
	*gradx = -acoef * (rxplus - rxmins) + (pow_dd(&d__1, &d__2) - pow_dd(&
		d__3, &d__4)) * .1515f;
	d__1 = (doublereal) ryplus;
	d__2 = (doublereal) expnt;
	d__3 = (doublereal) rymins;
	d__4 = (doublereal) expnt;
	*grady = -acoef * (ryplus - rymins) + (pow_dd(&d__1, &d__2) - pow_dd(&
		d__3, &d__4)) * .1515f;
	d__1 = (doublereal) rzplus;
	d__2 = (doublereal) expnt;
	d__3 = (doublereal) rzmins;
	d__4 = (doublereal) expnt;
	*gradz = -acoef * (rzplus - rzmins) + (pow_dd(&d__1, &d__2) - pow_dd(&
		d__3, &d__4)) * .1515f;
	return 0;
    }
    if (den0 > 4.5f && den0 <= 5.1f) {
/*       U=c1-ef*rho/rho0**2/3 */
	ef = .035999999999999997f;
	d__1 = (doublereal) rxplus;
	d__2 = (doublereal) rxmins;
	*gradx = ef * -.5f * (pow_dd(&d__1, &c_b407) - pow_dd(&d__2, &c_b407))
		;
	d__1 = (doublereal) ryplus;
	d__2 = (doublereal) rymins;
	*grady = ef * -.5f * (pow_dd(&d__1, &c_b407) - pow_dd(&d__2, &c_b407))
		;
	d__1 = (doublereal) rzplus;
	d__2 = (doublereal) rzmins;
	*gradz = ef * -.5f * (pow_dd(&d__1, &c_b407) - pow_dd(&d__2, &c_b407))
		;
	return 0;
    }
    if (den0 > 5.1f) {
/* U=800*(rho/rho0)**1/3.-Ef*(rho/rho0)**2/3.-c2 */
	ef = .035999999999999997f;
	cf0 = .8f;
	d__1 = (doublereal) rxplus;
	d__2 = (doublereal) rxmins;
	d__3 = (doublereal) rxplus;
	d__4 = (doublereal) rxmins;
	*gradx = cf0 * .5f * (pow_dd(&d__1, &c_b413) - pow_dd(&d__2, &c_b413))
		 - ef * .5f * (pow_dd(&d__3, &c_b407) - pow_dd(&d__4, &c_b407)
		);
	d__1 = (doublereal) ryplus;
	d__2 = (doublereal) rymins;
	d__3 = (doublereal) ryplus;
	d__4 = (doublereal) rymins;
	*grady = cf0 * .5f * (pow_dd(&d__1, &c_b413) - pow_dd(&d__2, &c_b413))
		 - ef * .5f * (pow_dd(&d__3, &c_b407) - pow_dd(&d__4, &c_b407)
		);
	d__1 = (doublereal) rzplus;
	d__2 = (doublereal) rzmins;
	d__3 = (doublereal) rzplus;
	d__4 = (doublereal) rzmins;
	*gradz = cf0 * .5f * (pow_dd(&d__1, &c_b413) - pow_dd(&d__2, &c_b413))
		 - ef * .5f * (pow_dd(&d__3, &c_b407) - pow_dd(&d__4, &c_b407)
		);
	return 0;
    }
    return 0;
} /* gradu_ */

/* ********************************* */
/*                                                                      * */
/* Subroutine */ int graduk_(integer *ix, integer *iy, integer *iz, real *
	gradxk, real *gradyk, real *gradzk)
{
    static real rxmins, rymins, rzmins, rxplus, ryplus, rzplus;

/*                                                                      * */
/*       PURPOSE:     DETERMINE the baryon density gradient for         * */
/*                    proporgating kaons in a mean field caused by      * */
/*                    surrounding baryons                               * */
/*       VARIABLES:                                                     * */
/*         IX, IY, IZ          - COORDINATES OF POINT   (INTEGER,INPUT) * */
/*         GRADXk, GRADYk, GRADZk                       (REAL,OUTPUT)   * */
/*                                                                      * */
/* ********************************* */

/* c      SAVE /DD/ */
/* c      SAVE /ss/ */

    rxplus = dd_1.rho[*ix + 1 + (*iy + *iz * 41) * 41 + 41184];
    rxmins = dd_1.rho[*ix - 1 + (*iy + *iz * 41) * 41 + 41184];
    ryplus = dd_1.rho[*ix + (*iy + 1 + *iz * 41) * 41 + 41184];
    rymins = dd_1.rho[*ix + (*iy - 1 + *iz * 41) * 41 + 41184];
    rzplus = dd_1.rho[*ix + (*iy + (*iz + 1) * 41) * 41 + 41184];
    rzmins = dd_1.rho[*ix + (*iy + (*iz - 1) * 41) * 41 + 41184];
    *gradxk = (rxplus - rxmins) / 2.f;
    *gradyk = (ryplus - rymins) / 2.f;
    *gradzk = (rzplus - rzmins) / 2.f;
    return 0;
} /* graduk_ */

/* ----------------------------------------------------------------------- */
/* Subroutine */ int gradup_(integer *ix, integer *iy, integer *iz, real *
	gradxp, real *gradyp, real *gradzp)
{
    static real rxmins, rymins, rzmins, rxplus, ryplus, rzplus;

/*                                                                      * */
/*       PURPOSE:     DETERMINE THE GRADIENT OF THE PROTON DENSITY      * */
/*       VARIABLES:                                                     * */
/*                                                                           * */
/*         IX, IY, IZ          - COORDINATES OF POINT   (INTEGER,INPUT) * */
/*         GRADXP, GRADYP, GRADZP - GRADIENT OF THE PROTON              * */
/*                                  DENSITY(REAL,OUTPUT)                * */
/*                                                                      * */
/* ********************************* */

/* c      SAVE /DD/ */
/* c      SAVE /ss/ */

    rxplus = dd_1.rhop[*ix + 1 + (*iy + *iz * 41) * 41 + 41184] / .168f;
    rxmins = dd_1.rhop[*ix - 1 + (*iy + *iz * 41) * 41 + 41184] / .168f;
    ryplus = dd_1.rhop[*ix + (*iy + 1 + *iz * 41) * 41 + 41184] / .168f;
    rymins = dd_1.rhop[*ix + (*iy - 1 + *iz * 41) * 41 + 41184] / .168f;
    rzplus = dd_1.rhop[*ix + (*iy + (*iz + 1) * 41) * 41 + 41184] / .168f;
    rzmins = dd_1.rhop[*ix + (*iy + (*iz - 1) * 41) * 41 + 41184] / .168f;
/* ----------------------------------------------------------------------- */

    *gradxp = (rxplus - rxmins) / 2.f;
    *gradyp = (ryplus - rymins) / 2.f;
    *gradzp = (rzplus - rzmins) / 2.f;
    return 0;
} /* gradup_ */

/* ----------------------------------------------------------------------- */
/* Subroutine */ int gradun_(integer *ix, integer *iy, integer *iz, real *
	gradxn, real *gradyn, real *gradzn)
{
    static real rxmins, rymins, rzmins, rxplus, ryplus, rzplus;

/*                                                                      * */
/*       PURPOSE:     DETERMINE THE GRADIENT OF THE NEUTRON DENSITY     * */
/*       VARIABLES:                                                     * */
/*                                                                           * */
/*         IX, IY, IZ          - COORDINATES OF POINT   (INTEGER,INPUT) * */
/*         GRADXN, GRADYN, GRADZN - GRADIENT OF THE NEUTRON             * */
/*                                  DENSITY(REAL,OUTPUT)                * */
/*                                                                      * */
/* ********************************* */

/* c      SAVE /DD/ */
/* c      SAVE /ss/ */

    rxplus = dd_1.rhon[*ix + 1 + (*iy + *iz * 41) * 41 + 41184] / .168f;
    rxmins = dd_1.rhon[*ix - 1 + (*iy + *iz * 41) * 41 + 41184] / .168f;
    ryplus = dd_1.rhon[*ix + (*iy + 1 + *iz * 41) * 41 + 41184] / .168f;
    rymins = dd_1.rhon[*ix + (*iy - 1 + *iz * 41) * 41 + 41184] / .168f;
    rzplus = dd_1.rhon[*ix + (*iy + (*iz + 1) * 41) * 41 + 41184] / .168f;
    rzmins = dd_1.rhon[*ix + (*iy + (*iz - 1) * 41) * 41 + 41184] / .168f;
/* ----------------------------------------------------------------------- */

    *gradxn = (rxplus - rxmins) / 2.f;
    *gradyn = (ryplus - rymins) / 2.f;
    *gradzn = (rzplus - rzmins) / 2.f;
    return 0;
} /* gradun_ */

/* ----------------------------------------------------------------------------- */
/* FUNCTION FDE(DMASS) GIVES DELTA MASS DISTRIBUTION BY USING OF */
/* KITAZOE'S FORMULA */
doublereal fde_(real *dmass, real *srt, real *con)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4, r__5, r__6;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real p1, fd, p11, am0, amn, avpi;
    extern doublereal width_(real *);

    amn = .938869f;
    avpi = .13803333f;
    am0 = 1.232f;
/* Computing 2nd power */
    r__1 = am0;
/* Computing 2nd power */
    r__3 = *dmass;
/* Computing 2nd power */
    r__2 = r__3 * r__3 - 1.5178240000000001f;
/* Computing 2nd power */
    r__4 = am0;
/* Computing 2nd power */
    r__5 = width_(dmass);
    fd = r__1 * r__1 * 4.f * width_(dmass) / (r__2 * r__2 + r__4 * r__4 * (
	    r__5 * r__5));
    if (*con == 1.f) {
/* Computing 2nd power */
	r__2 = *srt;
/* Computing 2nd power */
	r__3 = *dmass;
/* Computing 2nd power */
	r__4 = amn;
/* Computing 2nd power */
	r__1 = r__2 * r__2 + r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
	r__5 = *srt;
/* Computing 2nd power */
	r__6 = *dmass;
	p11 = r__1 * r__1 / (r__5 * r__5 * 4.f) - r__6 * r__6;
	if (p11 <= 0.f) {
	    p11 = 1e-6f;
	}
	p1 = sqrt(p11);
    } else {
	*dmass = amn + avpi;
/* Computing 2nd power */
	r__2 = *srt;
/* Computing 2nd power */
	r__3 = *dmass;
/* Computing 2nd power */
	r__4 = amn;
/* Computing 2nd power */
	r__1 = r__2 * r__2 + r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
	r__5 = *srt;
/* Computing 2nd power */
	r__6 = *dmass;
	p11 = r__1 * r__1 / (r__5 * r__5 * 4.f) - r__6 * r__6;
	if (p11 <= 0.f) {
	    p11 = 1e-6f;
	}
	p1 = sqrt(p11);
    }
    ret_val = fd * p1 * *dmass;
    return ret_val;
} /* fde_ */

/* ------------------------------------------------------------- */
/* FUNCTION FDE(DMASS) GIVES N*(1535) MASS DISTRIBUTION BY USING OF */
/* KITAZOE'S FORMULA */
doublereal fd5_(real *dmass, real *srt, real *con)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4, r__5, r__6;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double sqrt(doublereal);

    /* Local variables */
    static real p1, fd, am0, amn;
    extern doublereal w1535_(real *);
    static real avpi, scheck;

    /* Fortran I/O blocks */
    static cilist io___966 = { 0, 99, 0, 0, 0 };
    static cilist io___968 = { 0, 99, 0, 0, 0 };


    amn = .938869f;
    avpi = .13803333f;
    am0 = 1.535f;
/* Computing 2nd power */
    r__1 = am0;
/* Computing 2nd power */
    r__3 = *dmass;
/* Computing 2nd power */
    r__2 = r__3 * r__3 - 2.3562249999999998f;
/* Computing 2nd power */
    r__4 = am0;
/* Computing 2nd power */
    r__5 = w1535_(dmass);
    fd = r__1 * r__1 * 4.f * w1535_(dmass) / (r__2 * r__2 + r__4 * r__4 * (
	    r__5 * r__5));
    if (*con == 1.f) {
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
	r__2 = *srt;
/* Computing 2nd power */
	r__3 = *dmass;
/* Computing 2nd power */
	r__4 = amn;
/* Computing 2nd power */
	r__1 = r__2 * r__2 + r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
	r__5 = *srt;
/* Computing 2nd power */
	r__6 = *dmass;
	scheck = r__1 * r__1 / (r__5 * r__5 * 4.f) - r__6 * r__6;
	if (scheck < 0.f) {
	    s_wsle(&io___966);
	    do_lio(&c__9, &c__1, "scheck11: ", (ftnlen)10);
	    do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	    e_wsle();
	    scheck = 0.f;
	}
	p1 = sqrt(scheck);
/*           P1=SQRT((SRT**2+DMASS**2-AMN**2)**2 */
/*     1          /(4.*SRT**2)-DMASS**2) */
    } else {
	*dmass = amn + avpi;
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
	r__2 = *srt;
/* Computing 2nd power */
	r__3 = *dmass;
/* Computing 2nd power */
	r__4 = amn;
/* Computing 2nd power */
	r__1 = r__2 * r__2 + r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
	r__5 = *srt;
/* Computing 2nd power */
	r__6 = *dmass;
	scheck = r__1 * r__1 / (r__5 * r__5 * 4.f) - r__6 * r__6;
	if (scheck < 0.f) {
	    s_wsle(&io___968);
	    do_lio(&c__9, &c__1, "scheck12: ", (ftnlen)10);
	    do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	    e_wsle();
	    scheck = 0.f;
	}
	p1 = sqrt(scheck);
/*        P1=SQRT((SRT**2+DMASS**2-AMN**2)**2 */
/*     1  /(4.*SRT**2)-DMASS**2) */
    }
    ret_val = fd * p1 * *dmass;
    return ret_val;
} /* fd5_ */

/* -------------------------------------------------------------------------- */
/* FUNCTION FNS(DMASS) GIVES N* MASS DISTRIBUTION */
/*     BY USING OF BREIT-WIGNER FORMULA */
doublereal fns_(real *dmass, real *srt, real *con)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4, r__5, r__6;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double sqrt(doublereal);

    /* Local variables */
    static real p1, fn, an0, amn, avpi, width, scheck;

    /* Fortran I/O blocks */
    static cilist io___975 = { 0, 99, 0, 0, 0 };
    static cilist io___977 = { 0, 99, 0, 0, 0 };


    width = .2f;
    amn = .938869f;
    avpi = .13803333f;
    an0 = 1.43f;
/* Computing 2nd power */
    r__1 = an0;
/* Computing 2nd power */
    r__3 = *dmass;
/* Computing 2nd power */
    r__2 = r__3 * r__3 - 2.0735999999999999f;
/* Computing 2nd power */
    r__4 = an0;
/* Computing 2nd power */
    r__5 = width;
    fn = r__1 * r__1 * 4.f * width / (r__2 * r__2 + r__4 * r__4 * (r__5 * 
	    r__5));
    if (*con == 1.f) {
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
	r__2 = *srt;
/* Computing 2nd power */
	r__3 = *dmass;
/* Computing 2nd power */
	r__4 = amn;
/* Computing 2nd power */
	r__1 = r__2 * r__2 + r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
	r__5 = *srt;
/* Computing 2nd power */
	r__6 = *dmass;
	scheck = r__1 * r__1 / (r__5 * r__5 * 4.f) - r__6 * r__6;
	if (scheck < 0.f) {
	    s_wsle(&io___975);
	    do_lio(&c__9, &c__1, "scheck13: ", (ftnlen)10);
	    do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	    e_wsle();
	    scheck = 0.f;
	}
	p1 = sqrt(scheck);
/*        P1=SQRT((SRT**2+DMASS**2-AMN**2)**2 */
/*     1  /(4.*SRT**2)-DMASS**2) */
    } else {
	*dmass = amn + avpi;
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
	r__2 = *srt;
/* Computing 2nd power */
	r__3 = *dmass;
/* Computing 2nd power */
	r__4 = amn;
/* Computing 2nd power */
	r__1 = r__2 * r__2 + r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
	r__5 = *srt;
/* Computing 2nd power */
	r__6 = *dmass;
	scheck = r__1 * r__1 / (r__5 * r__5 * 4.f) - r__6 * r__6;
	if (scheck < 0.f) {
	    s_wsle(&io___977);
	    do_lio(&c__9, &c__1, "scheck14: ", (ftnlen)10);
	    do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	    e_wsle();
	    scheck = 0.f;
	}
	p1 = sqrt(scheck);
/*        P1=SQRT((SRT**2+DMASS**2-AMN**2)**2 */
/*     1  /(4.*SRT**2)-DMASS**2) */
    }
    ret_val = fn * p1 * *dmass;
    return ret_val;
} /* fns_ */

/* ----------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------------- */
/* PURPOSE:1. SORT N*(1440) and N*(1535) 2-body DECAY PRODUCTS */
/*         2. DETERMINE THE MOMENTUM AND COORDINATES OF NUCLEON AND PION */
/*            AFTER THE DELTA OR N* DECAYING */
/* DATE   : JAN. 24,1990, MODIFIED ON MAY 17, 1994 TO INCLUDE ETA */
/* Subroutine */ int decay_(integer *irun, integer *i__, integer *nnn, 
	integer *iseed, real *wid, integer *nt)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static real x3, x4, x5, x6, x8, dm;
    static integer lbi, lbm, nlab, nalb;
    static real ctrl;
    extern /* Subroutine */ int dkine_(integer *, integer *, integer *, 
	    integer *, integer *, real *, integer *);
    static integer lbanti, lbsave;
    static real dpsave;
    extern doublereal ranart_(integer *);
    static real xmsave, pxsave, pysave, pzsave;

/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /RUN/ */
/* c      SAVE /PA/ */
/* c      SAVE /PB/ */
/* c      SAVE /PC/ */
/* c      SAVE /PD/ */
/* c      SAVE /INPUT2/ */
/* c      SAVE /RNDF77/ */
    lbanti = ee_1.lb[*i__ - 1];

    dm = cc_1.e[*i__ - 1];
/* 1. FOR N*+(1440) DECAY */
    if ((i__1 = ee_1.lb[*i__ - 1], abs(i__1)) == 11) {
	x3 = ranart_(&rndf77_1.nseed);
	if (x3 > .33333333333333331f) {
	    ee_1.lb[*i__ - 1] = 2;
	    nlab = 2;
	    pd_1.lpion[*nnn + *irun * 150001 - 150002] = 5;
	    pc_1.epion[*nnn + *irun * 150001 - 150002] = .13957f;
	} else {
	    ee_1.lb[*i__ - 1] = 1;
	    nlab = 1;
	    pd_1.lpion[*nnn + *irun * 150001 - 150002] = 4;
	    pc_1.epion[*nnn + *irun * 150001 - 150002] = .13496f;
	}
/* 2. FOR N*0(1440) DECAY */
    } else if ((i__1 = ee_1.lb[*i__ - 1], abs(i__1)) == 10) {
	x4 = ranart_(&rndf77_1.nseed);
	if (x4 > .33333333333333331f) {
	    ee_1.lb[*i__ - 1] = 1;
	    nlab = 1;
	    pd_1.lpion[*nnn + *irun * 150001 - 150002] = 3;
	    pc_1.epion[*nnn + *irun * 150001 - 150002] = .13957f;
	} else {
	    ee_1.lb[*i__ - 1] = 2;
	    nalb = 2;
	    pd_1.lpion[*nnn + *irun * 150001 - 150002] = 4;
	    pc_1.epion[*nnn + *irun * 150001 - 150002] = .13496f;
	}
/* N*(1535) CAN DECAY TO A PION OR AN ETA IF DM > 1.49 GeV */
/* 3 N*(0)(1535) DECAY */
    } else if ((i__1 = ee_1.lb[*i__ - 1], abs(i__1)) == 12) {
	ctrl = .65f;
	if (dm <= 1.49f) {
	    ctrl = -1.f;
	}
	x5 = ranart_(&rndf77_1.nseed);
	if (x5 >= ctrl) {
/* DECAY TO PION+NUCLEON */
	    x6 = ranart_(&rndf77_1.nseed);
	    if (x6 > .33333333333333331f) {
		ee_1.lb[*i__ - 1] = 1;
		nlab = 1;
		pd_1.lpion[*nnn + *irun * 150001 - 150002] = 3;
		pc_1.epion[*nnn + *irun * 150001 - 150002] = .13957f;
	    } else {
		ee_1.lb[*i__ - 1] = 2;
		nalb = 2;
		pd_1.lpion[*nnn + *irun * 150001 - 150002] = 4;
		pc_1.epion[*nnn + *irun * 150001 - 150002] = .13496f;
	    }
	} else {
/* DECAY TO ETA+NEUTRON */
	    ee_1.lb[*i__ - 1] = 2;
	    nlab = 2;
	    pd_1.lpion[*nnn + *irun * 150001 - 150002] = 0;
	    pc_1.epion[*nnn + *irun * 150001 - 150002] = .5475f;
	}
/* 4. FOR N*+(1535) DECAY */
    } else if ((i__1 = ee_1.lb[*i__ - 1], abs(i__1)) == 13) {
	ctrl = .65f;
	if (dm <= 1.49f) {
	    ctrl = -1.f;
	}
	x5 = ranart_(&rndf77_1.nseed);
	if (x5 >= ctrl) {
/* DECAY TO PION+NUCLEON */
	    x8 = ranart_(&rndf77_1.nseed);
	    if (x8 > .33333333333333331f) {
		ee_1.lb[*i__ - 1] = 2;
		nlab = 2;
		pd_1.lpion[*nnn + *irun * 150001 - 150002] = 5;
		pc_1.epion[*nnn + *irun * 150001 - 150002] = .13957f;
	    } else {
		ee_1.lb[*i__ - 1] = 1;
		nlab = 1;
		pd_1.lpion[*nnn + *irun * 150001 - 150002] = 4;
		pc_1.epion[*nnn + *irun * 150001 - 150002] = .13496f;
	    }
	} else {
/* DECAY TO ETA+NUCLEON */
	    ee_1.lb[*i__ - 1] = 1;
	    nlab = 1;
	    pd_1.lpion[*nnn + *irun * 150001 - 150002] = 0;
	    pc_1.epion[*nnn + *irun * 150001 - 150002] = .5475f;
	}
    }

    dkine_(irun, i__, nnn, &nlab, iseed, wid, nt);

/*     anti-particle ID for anti-N* decays: */
    if (lbanti < 0) {
	lbi = ee_1.lb[*i__ - 1];
	if (lbi == 1 || lbi == 2) {
	    lbi = -lbi;
	} else if (lbi == 3) {
	    lbi = 5;
	} else if (lbi == 5) {
	    lbi = 3;
	}
	ee_1.lb[*i__ - 1] = lbi;

	lbi = pd_1.lpion[*nnn + *irun * 150001 - 150002];
	if (lbi == 3) {
	    lbi = 5;
	} else if (lbi == 5) {
	    lbi = 3;
	} else if (lbi == 1 || lbi == 2) {
	    lbi = -lbi;
	}
	pd_1.lpion[*nnn + *irun * 150001 - 150002] = lbi;
    }

    if (*nt == input2_1.ntmax) {
/*     at the last timestep, assign rho or eta (decay daughter) */
/*     to lb(i1) only (not to lpion) in order to decay them again: */
	lbm = pd_1.lpion[*nnn + *irun * 150001 - 150002];
	if (lbm == 0 || lbm == 25 || lbm == 26 || lbm == 27) {
/*     switch rho or eta with baryon, positions are the same (no change needed): */
	    lbsave = lbm;
	    xmsave = pc_1.epion[*nnn + *irun * 150001 - 150002];
	    pxsave = pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450006];
	    pysave = pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450005];
	    pzsave = pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450004];
/* lin-5/2008: */
	    dpsave = dpert_1.dppion[*nnn + *irun * 150001 - 150002];
	    pd_1.lpion[*nnn + *irun * 150001 - 150002] = ee_1.lb[*i__ - 1];
	    pc_1.epion[*nnn + *irun * 150001 - 150002] = cc_1.e[*i__ - 1];
	    pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450006] = bb_1.p[*i__ * 
		    3 - 3];
	    pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450005] = bb_1.p[*i__ * 
		    3 - 2];
	    pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450004] = bb_1.p[*i__ * 
		    3 - 1];
/* lin-5/2008: */
	    dpert_1.dppion[*nnn + *irun * 150001 - 150002] = dpert_1.dpertp[*
		    i__ - 1];
	    ee_1.lb[*i__ - 1] = lbsave;
	    cc_1.e[*i__ - 1] = xmsave;
	    bb_1.p[*i__ * 3 - 3] = pxsave;
	    bb_1.p[*i__ * 3 - 2] = pysave;
	    bb_1.p[*i__ * 3 - 1] = pzsave;
/* lin-5/2008: */
	    dpert_1.dpertp[*i__ - 1] = dpsave;
	}
    }
    return 0;
} /* decay_ */

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* PURPOSE: */
/*         CALCULATE THE MOMENTUM OF NUCLEON AND PION (OR ETA) */
/*         IN THE LAB. FRAME AFTER DELTA OR N* DECAY */
/* DATE   : JAN. 24,1990, MODIFIED ON MAY 17, 1994 TO INCLUDE ETA PRODUCTION */
/* Subroutine */ int dkine_(integer *irun, integer *i__, integer *nnn, 
	integer *nlab, integer *iseed, real *wid, integer *nt)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal);

    /* Local variables */
    static real q, q2, gd, am, dm, en, ep, pm, qs, px, py, pz, rx, ry, rz, qx,
	     qy, qz, fgd, bdx, bdy, bdz, bpp, bpn, pxn, pyn, pxp, pyp, pzp, 
	    pzn, tau0, devio, edelta;
    extern doublereal ranart_(integer *);
    static real taudcy;

/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /RUN/ */
/* c      SAVE /PA/ */
/* c      SAVE /PB/ */
/* c      SAVE /PC/ */
/* c      SAVE /PD/ */
/* c      SAVE /leadng/ */
/* c      SAVE /tdecay/ */
/* c      SAVE /INPUT2/ */
/* c      SAVE /RNDF77/ */
/* READ IN THE COORDINATES OF DELTA OR N* UNDERGOING DECAY */
    px = bb_1.p[*i__ * 3 - 3];
    py = bb_1.p[*i__ * 3 - 2];
    pz = bb_1.p[*i__ * 3 - 1];
    rx = aa_1.r__[*i__ * 3 - 3];
    ry = aa_1.r__[*i__ * 3 - 2];
    rz = aa_1.r__[*i__ * 3 - 1];
    dm = cc_1.e[*i__ - 1];
/* Computing 2nd power */
    r__1 = dm;
/* Computing 2nd power */
    r__2 = px;
/* Computing 2nd power */
    r__3 = py;
/* Computing 2nd power */
    r__4 = pz;
    edelta = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    pm = pc_1.epion[*nnn + *irun * 150001 - 150002];
    am = .93828f;
    if (*nlab == 2) {
	am = .939457f;
    }
/* FIND OUT THE MOMENTUM AND ENERGY OF PION AND NUCLEON IN DELTA REST FRAME */
/* THE MAGNITUDE OF MOMENTUM IS DETERMINED BY ENERGY CONSERVATION ,THE FORMULA */
/* CAN BE FOUND ON PAGE 716,W BAUER P.R.C40,1989 */
/* THE DIRECTION OF THE MOMENTUM IS ASSUMED ISOTROPIC. NOTE THAT P(PION)=-P(N) */
/* Computing 2nd power */
    r__2 = dm;
/* Computing 2nd power */
    r__3 = am;
/* Computing 2nd power */
    r__4 = pm;
/* Computing 2nd power */
    r__1 = (r__2 * r__2 - r__3 * r__3 + r__4 * r__4) / (dm * 2.f);
/* Computing 2nd power */
    r__5 = pm;
    q2 = r__1 * r__1 - r__5 * r__5;
    if (q2 <= 0.f) {
	q2 = 1e-9f;
    }
    q = sqrt(q2);
L11:
    qx = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
    qy = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
    qz = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
/* Computing 2nd power */
    r__1 = qx;
/* Computing 2nd power */
    r__2 = qy;
/* Computing 2nd power */
    r__3 = qz;
    qs = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
    if (qs > 1.f) {
	goto L11;
    }
    pxp = q * qx / sqrt(qs);
    pyp = q * qy / sqrt(qs);
    pzp = q * qz / sqrt(qs);
/* Computing 2nd power */
    r__1 = q;
/* Computing 2nd power */
    r__2 = pm;
    ep = sqrt(r__1 * r__1 + r__2 * r__2);
    pxn = -pxp;
    pyn = -pyp;
    pzn = -pzp;
/* Computing 2nd power */
    r__1 = q;
/* Computing 2nd power */
    r__2 = am;
    en = sqrt(r__1 * r__1 + r__2 * r__2);
/* TRANSFORM INTO THE LAB. FRAME. THE GENERAL LORENTZ TRANSFORMATION CAN */
/* BE FOUND ON PAGE 34 OF R. HAGEDORN " RELATIVISTIC KINEMATICS" */
    gd = edelta / dm;
    fgd = gd / (gd + 1.f);
    bdx = px / edelta;
    bdy = py / edelta;
    bdz = pz / edelta;
    bpp = bdx * pxp + bdy * pyp + bdz * pzp;
    bpn = bdx * pxn + bdy * pyn + bdz * pzn;
    bb_1.p[*i__ * 3 - 3] = pxn + bdx * gd * (fgd * bpn + en);
    bb_1.p[*i__ * 3 - 2] = pyn + bdy * gd * (fgd * bpn + en);
    bb_1.p[*i__ * 3 - 1] = pzn + bdz * gd * (fgd * bpn + en);
    cc_1.e[*i__ - 1] = am;
/* WE ASSUME THAT THE SPACIAL COORDINATE OF THE NUCLEON */
/* IS THAT OF THE DELTA */
    pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450006] = pxp + bdx * gd * (fgd *
	     bpp + ep);
    pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450005] = pyp + bdy * gd * (fgd *
	     bpp + ep);
    pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450004] = pzp + bdz * gd * (fgd *
	     bpp + ep);
/* lin-5/2008: */
    dpert_1.dppion[*nnn + *irun * 150001 - 150002] = dpert_1.dpertp[*i__ - 1];
/* WE ASSUME THE PION OR ETA COMING FROM DELTA DECAY IS LOCATED ON THE SPHERE */
/* OF RADIUS 0.5FM AROUND DELTA, THIS POINT NEED TO BE CHECKED */
/* AND OTHER CRIERTION MAY BE TRIED */
/* lin-2/20/03 no additional smearing for position of decay daughters: */
/* 200         X0 = 1.0 - 2.0 * RANART(NSEED) */
/*            Y0 = 1.0 - 2.0 * RANART(NSEED) */
/*            Z0 = 1.0 - 2.0 * RANART(NSEED) */
/*        IF ((X0*X0+Y0*Y0+Z0*Z0) .GT. 1.0) GOTO 200 */
/*        RPION(1,NNN,IRUN)=R(1,I)+0.5*x0 */
/*        RPION(2,NNN,IRUN)=R(2,I)+0.5*y0 */
/*        RPION(3,NNN,IRUN)=R(3,I)+0.5*z0 */
    pa_1.rpion[(*nnn + *irun * 150001) * 3 - 450006] = aa_1.r__[*i__ * 3 - 3];
    pa_1.rpion[(*nnn + *irun * 150001) * 3 - 450005] = aa_1.r__[*i__ * 3 - 2];
    pa_1.rpion[(*nnn + *irun * 150001) * 3 - 450004] = aa_1.r__[*i__ * 3 - 1];

/* Computing 2nd power */
    r__1 = pc_1.epion[*nnn + *irun * 150001 - 150002];
/* Computing 2nd power */
    r__2 = pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450006];
/* Computing 2nd power */
    r__3 = pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450005];
/* Computing 2nd power */
    r__4 = pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450004];
/* Computing 2nd power */
    r__5 = cc_1.e[*i__ - 1];
/* Computing 2nd power */
    r__6 = bb_1.p[*i__ * 3 - 3];
/* Computing 2nd power */
    r__7 = bb_1.p[*i__ * 3 - 2];
/* Computing 2nd power */
    r__8 = bb_1.p[*i__ * 3 - 1];
    devio = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4) + 
	    sqrt(r__5 * r__5 + r__6 * r__6 + r__7 * r__7 + r__8 * r__8) - 
	    leadng_1.e1;
/*        if(abs(devio).gt.0.02) write(93,*) 'decay(): nt=',nt,devio,lb1 */
/*     add decay time to daughter's formation time at the last timestep: */
    if (*nt == input2_1.ntmax) {
	tau0 = .19733f / *wid;
	taudcy = tau0 * -1.f * log(1.f - ranart_(&rndf77_1.nseed));
/*     lorentz boost: */
	taudcy = taudcy * leadng_1.e1 / leadng_1.em1;
	leadng_1.tfnl += taudcy;
	leadng_1.xfnl += leadng_1.px1 / leadng_1.e1 * taudcy;
	leadng_1.yfnl += leadng_1.py1 / leadng_1.e1 * taudcy;
	leadng_1.zfnl += leadng_1.pz1 / leadng_1.e1 * taudcy;
	aa_1.r__[*i__ * 3 - 3] = leadng_1.xfnl;
	aa_1.r__[*i__ * 3 - 2] = leadng_1.yfnl;
	aa_1.r__[*i__ * 3 - 1] = leadng_1.zfnl;
	tdecay_1.tfdcy[*i__ - 1] = leadng_1.tfnl;
	pa_1.rpion[(*nnn + *irun * 150001) * 3 - 450006] = leadng_1.xfnl;
	pa_1.rpion[(*nnn + *irun * 150001) * 3 - 450005] = leadng_1.yfnl;
	pa_1.rpion[(*nnn + *irun * 150001) * 3 - 450004] = leadng_1.zfnl;
	tdecay_1.tfdpi[*nnn + *irun * 150001 - 150002] = leadng_1.tfnl;
    }
/* L200: */
/* L210: */
/* L220: */
    return 0;
} /* dkine_ */

/* ----------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------------- */
/* PURPOSE:1. N*-->N+PION+PION  DECAY PRODUCTS */
/*         2. DETERMINE THE MOMENTUM AND COORDINATES OF NUCLEON AND PION */
/*            AFTER THE DELTA OR N* DECAYING */
/* DATE   : NOV.7,1994 */
/* ---------------------------------------------------------------------------- */
/* Subroutine */ int decay2_(integer *irun, integer *i__, integer *nnn, 
	integer *iseed, real *wid, integer *nt)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static real x3, dm;
    static integer lbi, nlab;
    extern /* Subroutine */ int dkine2_(integer *, integer *, integer *, 
	    integer *, integer *, real *, integer *);
    static integer lbanti;
    extern doublereal ranart_(integer *);

/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /RUN/ */
/* c      SAVE /PA/ */
/* c      SAVE /PB/ */
/* c      SAVE /PC/ */
/* c      SAVE /PD/ */
/* c      SAVE /RNDF77/ */
    lbanti = ee_1.lb[*i__ - 1];

    dm = cc_1.e[*i__ - 1];
/* DETERMINE THE DECAY PRODUCTS */
/* FOR N*+(1440) DECAY */
    if ((i__1 = ee_1.lb[*i__ - 1], abs(i__1)) == 11) {
	x3 = ranart_(&rndf77_1.nseed);
	if (x3 < .33333333333333331f) {
	    ee_1.lb[*i__ - 1] = 2;
	    nlab = 2;
	    pd_1.lpion[*nnn + *irun * 150001 - 150002] = 5;
	    pc_1.epion[*nnn + *irun * 150001 - 150002] = .13957f;
	    pd_1.lpion[*nnn + 1 + *irun * 150001 - 150002] = 4;
	    pc_1.epion[*nnn + 1 + *irun * 150001 - 150002] = .13496f;
	} else if (x3 < .66666666666666663f && x3 > .33333333333333331f) {
	    ee_1.lb[*i__ - 1] = 1;
	    nlab = 1;
	    pd_1.lpion[*nnn + *irun * 150001 - 150002] = 5;
	    pc_1.epion[*nnn + *irun * 150001 - 150002] = .13957f;
	    pd_1.lpion[*nnn + 1 + *irun * 150001 - 150002] = 3;
	    pc_1.epion[*nnn + 1 + *irun * 150001 - 150002] = .13957f;
	} else {
	    ee_1.lb[*i__ - 1] = 1;
	    nlab = 1;
	    pd_1.lpion[*nnn + *irun * 150001 - 150002] = 4;
	    pc_1.epion[*nnn + *irun * 150001 - 150002] = .13496f;
	    pd_1.lpion[*nnn + 1 + *irun * 150001 - 150002] = 4;
	    pc_1.epion[*nnn + 1 + *irun * 150001 - 150002] = .13496f;
	}
/* FOR N*0(1440) DECAY */
    } else if ((i__1 = ee_1.lb[*i__ - 1], abs(i__1)) == 10) {
	x3 = ranart_(&rndf77_1.nseed);
	if (x3 < .33333333333333331f) {
	    ee_1.lb[*i__ - 1] = 2;
	    nlab = 2;
	    pd_1.lpion[*nnn + *irun * 150001 - 150002] = 4;
	    pc_1.epion[*nnn + *irun * 150001 - 150002] = .13496f;
	    pd_1.lpion[*nnn + 1 + *irun * 150001 - 150002] = 4;
	    pc_1.epion[*nnn + 1 + *irun * 150001 - 150002] = .13496f;
	} else if (x3 < .66666666666666663f && x3 > .33333333333333331f) {
	    ee_1.lb[*i__ - 1] = 1;
	    nlab = 1;
	    pd_1.lpion[*nnn + *irun * 150001 - 150002] = 3;
	    pc_1.epion[*nnn + *irun * 150001 - 150002] = .13957f;
	    pd_1.lpion[*nnn + 1 + *irun * 150001 - 150002] = 4;
	    pc_1.epion[*nnn + 1 + *irun * 150001 - 150002] = .13496f;
	} else {
	    ee_1.lb[*i__ - 1] = 2;
	    nlab = 2;
	    pd_1.lpion[*nnn + *irun * 150001 - 150002] = 5;
	    pc_1.epion[*nnn + *irun * 150001 - 150002] = .13957f;
	    pd_1.lpion[*nnn + 1 + *irun * 150001 - 150002] = 3;
	    pc_1.epion[*nnn + 1 + *irun * 150001 - 150002] = .13957f;
	}
    }
    dkine2_(irun, i__, nnn, &nlab, iseed, wid, nt);

/*     anti-particle ID for anti-N* decays: */
    if (lbanti < 0) {
	lbi = ee_1.lb[*i__ - 1];
	if (lbi == 1 || lbi == 2) {
	    lbi = -lbi;
	} else if (lbi == 3) {
	    lbi = 5;
	} else if (lbi == 5) {
	    lbi = 3;
	}
	ee_1.lb[*i__ - 1] = lbi;

	lbi = pd_1.lpion[*nnn + *irun * 150001 - 150002];
	if (lbi == 3) {
	    lbi = 5;
	} else if (lbi == 5) {
	    lbi = 3;
	} else if (lbi == 1 || lbi == 2) {
	    lbi = -lbi;
	}
	pd_1.lpion[*nnn + *irun * 150001 - 150002] = lbi;

	lbi = pd_1.lpion[*nnn + 1 + *irun * 150001 - 150002];
	if (lbi == 3) {
	    lbi = 5;
	} else if (lbi == 5) {
	    lbi = 3;
	} else if (lbi == 1 || lbi == 2) {
	    lbi = -lbi;
	}
	pd_1.lpion[*nnn + 1 + *irun * 150001 - 150002] = lbi;
    }

    return 0;
} /* decay2_ */

/* ------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/*         CALCULATE THE MOMENTUM OF NUCLEON AND PION (OR ETA) */
/*         IN THE LAB. FRAME AFTER DELTA OR N* DECAY */
/* DATE   : JAN. 24,1990, MODIFIED ON MAY 17, 1994 TO INCLUDE ETA PRODUCTION */
/* -------------------------------------------------------------------------- */
/* Subroutine */ int dkine2_(integer *irun, integer *i__, integer *nnn, 
	integer *nlab, integer *iseed, real *wid, integer *nt)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8, r__9, r__10, r__11, 
	    r__12;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double cos(doublereal), sin(doublereal);
    /* Subroutine */ int s_stop(char *, ftnlen);
    double log(doublereal);

    /* Local variables */
    static real q, q2, gd, am, dm, en, ep, qs, px, py, pz, rx, ry, rz, qx, qy,
	     qz, gd1, ep0, bp0, pm1, pm2, p1m, p2m, p3m, p1p, p2p, p3p, px0, 
	    py0, pz0, fai, fgd, bdx, bdy, bdz, bpp, epn, epp, bpn, css, pxn, 
	    pyn, pxp, sss, pyp, pzp, pzn, fgd1, bpn1, bpp1, tau0, pmax, pmax2,
	     betax, betay, betaz, enucl, devio, epion1, epion2, edelta, 
	    scheck;
    extern doublereal ranart_(integer *);
    static real taudcy;

    /* Fortran I/O blocks */
    static cilist io___1048 = { 0, 99, 0, 0, 0 };
    static cilist io___1060 = { 0, 99, 0, 0, 0 };


/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /RUN/ */
/* c      SAVE /PA/ */
/* c      SAVE /PB/ */
/* c      SAVE /PC/ */
/* c      SAVE /PD/ */
/* c      SAVE /leadng/ */
/* c      SAVE /tdecay/ */
/* c      SAVE /INPUT2/ */
/* c      SAVE /RNDF77/ */
/* READ IN THE COORDINATES OF THE N*(1440) UNDERGOING DECAY */
    px = bb_1.p[*i__ * 3 - 3];
    py = bb_1.p[*i__ * 3 - 2];
    pz = bb_1.p[*i__ * 3 - 1];
    rx = aa_1.r__[*i__ * 3 - 3];
    ry = aa_1.r__[*i__ * 3 - 2];
    rz = aa_1.r__[*i__ * 3 - 1];
    dm = cc_1.e[*i__ - 1];
/* Computing 2nd power */
    r__1 = dm;
/* Computing 2nd power */
    r__2 = px;
/* Computing 2nd power */
    r__3 = py;
/* Computing 2nd power */
    r__4 = pz;
    edelta = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    pm1 = pc_1.epion[*nnn + *irun * 150001 - 150002];
    pm2 = pc_1.epion[*nnn + 1 + *irun * 150001 - 150002];
    am = .939457f;
    if (*nlab == 1) {
	am = .93828f;
    }
/* THE MAXIMUM MOMENTUM OF THE NUCLEON FROM THE DECAY OF A N* */
/* Computing 2nd power */
    r__1 = dm;
/* Computing 2nd power */
    r__2 = am + pm1 + pm2;
/* Computing 2nd power */
    r__3 = dm;
/* Computing 2nd power */
    r__4 = am - pm1 - pm2;
/* Computing 2nd power */
    r__5 = dm;
    pmax2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / 4 / (
	    r__5 * r__5);
/* lin-9/2012: check argument in sqrt(): */
    scheck = pmax2;
    if (scheck < 0.f) {
	s_wsle(&io___1048);
	do_lio(&c__9, &c__1, "scheck15: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	scheck = 0.f;
    }
    pmax = sqrt(scheck);
/*       PMAX=SQRT(PMAX2) */
/* GENERATE THE MOMENTUM OF THE NUCLEON IN THE N* REST FRAME */
    css = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
/* Computing 2nd power */
    r__1 = css;
    sss = sqrt(1 - r__1 * r__1);
    fai = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
    px0 = pmax * sss * cos(fai);
    py0 = pmax * sss * sin(fai);
    pz0 = pmax * css;
/* Computing 2nd power */
    r__1 = px0;
/* Computing 2nd power */
    r__2 = py0;
/* Computing 2nd power */
    r__3 = pz0;
/* Computing 2nd power */
    r__4 = am;
    ep0 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
/* lin-5/23/01 bug: P0 for pion0 is equal to PMAX, leaving pion+ and pion- */
/*     without no relative momentum, thus producing them with equal momenta, */
/* BETA AND GAMMA OF THE CMS OF PION+-PION- */
    betax = -px0 / (dm - ep0);
    betay = -py0 / (dm - ep0);
    betaz = -pz0 / (dm - ep0);
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__1 = betax;
/* Computing 2nd power */
    r__2 = betay;
/* Computing 2nd power */
    r__3 = betaz;
    scheck = 1 - r__1 * r__1 - r__2 * r__2 - r__3 * r__3;
    if (scheck <= 0.f) {
	s_wsle(&io___1060);
	do_lio(&c__9, &c__1, "scheck16: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	s_stop("", (ftnlen)0);
    }
    gd1 = 1.f / sqrt(scheck);
/*       GD1=1./SQRT(1-BETAX**2-BETAY**2-BETAZ**2) */
    fgd1 = gd1 / (gd1 + 1);
/* GENERATE THE MOMENTA OF PIONS IN THE CMS OF PION+PION- */
/* Computing 2nd power */
    r__1 = (dm - ep0) / (gd1 * 2.f);
/* Computing 2nd power */
    r__2 = pm1;
    q2 = r__1 * r__1 - r__2 * r__2;
    if (q2 <= 0.f) {
	q2 = 1e-9f;
    }
    q = sqrt(q2);
L11:
    qx = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
    qy = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
    qz = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
/* Computing 2nd power */
    r__1 = qx;
/* Computing 2nd power */
    r__2 = qy;
/* Computing 2nd power */
    r__3 = qz;
    qs = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
    if (qs > 1.f) {
	goto L11;
    }
    pxp = q * qx / sqrt(qs);
    pyp = q * qy / sqrt(qs);
    pzp = q * qz / sqrt(qs);
/* Computing 2nd power */
    r__1 = q;
/* Computing 2nd power */
    r__2 = pm1;
    ep = sqrt(r__1 * r__1 + r__2 * r__2);
    pxn = -pxp;
    pyn = -pyp;
    pzn = -pzp;
/* Computing 2nd power */
    r__1 = q;
/* Computing 2nd power */
    r__2 = pm2;
    en = sqrt(r__1 * r__1 + r__2 * r__2);
/* TRANSFORM THE MOMENTA OF PION+PION- INTO THE N* REST FRAME */
    bpp1 = betax * pxp + betay * pyp + betaz * pzp;
    bpn1 = betax * pxn + betay * pyn + betaz * pzn;
/* FOR PION- */
    p1m = pxn + betax * gd1 * (fgd1 * bpn1 + en);
    p2m = pyn + betay * gd1 * (fgd1 * bpn1 + en);
    p3m = pzn + betaz * gd1 * (fgd1 * bpn1 + en);
/* Computing 2nd power */
    r__1 = p1m;
/* Computing 2nd power */
    r__2 = p2m;
/* Computing 2nd power */
    r__3 = p3m;
/* Computing 2nd power */
    r__4 = pm2;
    epn = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
/* FOR PION+ */
    p1p = pxp + betax * gd1 * (fgd1 * bpp1 + ep);
    p2p = pyp + betay * gd1 * (fgd1 * bpp1 + ep);
    p3p = pzp + betaz * gd1 * (fgd1 * bpp1 + ep);
/* Computing 2nd power */
    r__1 = p1p;
/* Computing 2nd power */
    r__2 = p2p;
/* Computing 2nd power */
    r__3 = p3p;
/* Computing 2nd power */
    r__4 = pm1;
    epp = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
/* TRANSFORM MOMENTA OF THE THREE PIONS INTO THE */
/* THE NUCLEUS-NUCLEUS CENTER OF MASS  FRAME. */
/* THE GENERAL LORENTZ TRANSFORMATION CAN */
/* BE FOUND ON PAGE 34 OF R. HAGEDORN " RELATIVISTIC KINEMATICS" */
    gd = edelta / dm;
    fgd = gd / (gd + 1.f);
    bdx = px / edelta;
    bdy = py / edelta;
    bdz = pz / edelta;
    bp0 = bdx * px0 + bdy * py0 + bdz * pz0;
    bpp = bdx * p1p + bdy * p2p + bdz * p3p;
    bpn = bdx * p1m + bdy * p2m + bdz * p3m;
/* FOR THE NUCLEON */
    bb_1.p[*i__ * 3 - 3] = px0 + bdx * gd * (fgd * bp0 + ep0);
    bb_1.p[*i__ * 3 - 2] = py0 + bdy * gd * (fgd * bp0 + ep0);
    bb_1.p[*i__ * 3 - 1] = pz0 + bdz * gd * (fgd * bp0 + ep0);
    cc_1.e[*i__ - 1] = am;
    ee_1.id[*i__ - 1] = 0;
/* Computing 2nd power */
    r__1 = bb_1.p[*i__ * 3 - 3];
/* Computing 2nd power */
    r__2 = bb_1.p[*i__ * 3 - 2];
/* Computing 2nd power */
    r__3 = bb_1.p[*i__ * 3 - 1];
/* Computing 2nd power */
    r__4 = cc_1.e[*i__ - 1];
    enucl = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
/* WE ASSUME THAT THE SPACIAL COORDINATE OF THE PION0 */
/* IS in a sphere of radius 0.5 fm around N* */
/* FOR PION+ */
    pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450006] = p1p + bdx * gd * (fgd *
	     bpp + epp);
    pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450005] = p2p + bdy * gd * (fgd *
	     bpp + epp);
    pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450004] = p3p + bdz * gd * (fgd *
	     bpp + epp);
/* Computing 2nd power */
    r__1 = pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450006];
/* Computing 2nd power */
    r__2 = pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450005];
/* Computing 2nd power */
    r__3 = pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450004];
/* Computing 2nd power */
    r__4 = pc_1.epion[*nnn + *irun * 150001 - 150002];
    epion1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
/* lin-2/20/03 no additional smearing for position of decay daughters: */
/* 200         X0 = 1.0 - 2.0 * RANART(NSEED) */
/*            Y0 = 1.0 - 2.0 * RANART(NSEED) */
/*            Z0 = 1.0 - 2.0 * RANART(NSEED) */
/*        IF ((X0*X0+Y0*Y0+Z0*Z0) .GT. 1.0) GOTO 200 */
/*        RPION(1,NNN,IRUN)=R(1,I)+0.5*x0 */
/*        RPION(2,NNN,IRUN)=R(2,I)+0.5*y0 */
/*        RPION(3,NNN,IRUN)=R(3,I)+0.5*z0 */
    pa_1.rpion[(*nnn + *irun * 150001) * 3 - 450006] = aa_1.r__[*i__ * 3 - 3];
    pa_1.rpion[(*nnn + *irun * 150001) * 3 - 450005] = aa_1.r__[*i__ * 3 - 2];
    pa_1.rpion[(*nnn + *irun * 150001) * 3 - 450004] = aa_1.r__[*i__ * 3 - 1];
/* FOR PION- */
    pb_1.ppion[(*nnn + 1 + *irun * 150001) * 3 - 450006] = p1m + bdx * gd * (
	    fgd * bpn + epn);
    pb_1.ppion[(*nnn + 1 + *irun * 150001) * 3 - 450005] = p2m + bdy * gd * (
	    fgd * bpn + epn);
    pb_1.ppion[(*nnn + 1 + *irun * 150001) * 3 - 450004] = p3m + bdz * gd * (
	    fgd * bpn + epn);
/* lin-5/2008: */
    dpert_1.dppion[*nnn + *irun * 150001 - 150002] = dpert_1.dpertp[*i__ - 1];
    dpert_1.dppion[*nnn + 1 + *irun * 150001 - 150002] = dpert_1.dpertp[*i__ 
	    - 1];

/* Computing 2nd power */
    r__1 = pb_1.ppion[(*nnn + 1 + *irun * 150001) * 3 - 450006];
/* Computing 2nd power */
    r__2 = pb_1.ppion[(*nnn + 1 + *irun * 150001) * 3 - 450005];
/* Computing 2nd power */
    r__3 = pb_1.ppion[(*nnn + 1 + *irun * 150001) * 3 - 450004];
/* Computing 2nd power */
    r__4 = pc_1.epion[*nnn + 1 + *irun * 150001 - 150002];
    epion2 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
/* lin-2/20/03 no additional smearing for position of decay daughters: */
/* 300         X0 = 1.0 - 2.0 * RANART(NSEED) */
/*            Y0 = 1.0 - 2.0 * RANART(NSEED) */
/*            Z0 = 1.0 - 2.0 * RANART(NSEED) */
/*        IF ((X0*X0+Y0*Y0+Z0*Z0) .GT. 1.0) GOTO 300 */
/*        RPION(1,NNN+1,IRUN)=R(1,I)+0.5*x0 */
/*        RPION(2,NNN+1,IRUN)=R(2,I)+0.5*y0 */
/*        RPION(3,NNN+1,IRUN)=R(3,I)+0.5*z0 */
    pa_1.rpion[(*nnn + 1 + *irun * 150001) * 3 - 450006] = aa_1.r__[*i__ * 3 
	    - 3];
    pa_1.rpion[(*nnn + 1 + *irun * 150001) * 3 - 450005] = aa_1.r__[*i__ * 3 
	    - 2];
    pa_1.rpion[(*nnn + 1 + *irun * 150001) * 3 - 450004] = aa_1.r__[*i__ * 3 
	    - 1];

/* check energy conservation in the decay */
/*       efinal=enucl+epion1+epion2 */
/*       DEEE=(EDELTA-EFINAL)/EDELTA */
/*       IF(ABS(DEEE).GE.1.E-03)write(6,*)1,edelta,efinal */
/* Computing 2nd power */
    r__1 = pc_1.epion[*nnn + *irun * 150001 - 150002];
/* Computing 2nd power */
    r__2 = pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450006];
/* Computing 2nd power */
    r__3 = pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450005];
/* Computing 2nd power */
    r__4 = pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450004];
/* Computing 2nd power */
    r__5 = cc_1.e[*i__ - 1];
/* Computing 2nd power */
    r__6 = bb_1.p[*i__ * 3 - 3];
/* Computing 2nd power */
    r__7 = bb_1.p[*i__ * 3 - 2];
/* Computing 2nd power */
    r__8 = bb_1.p[*i__ * 3 - 1];
/* Computing 2nd power */
    r__9 = pc_1.epion[*nnn + 1 + *irun * 150001 - 150002];
/* Computing 2nd power */
    r__10 = pb_1.ppion[(*nnn + 1 + *irun * 150001) * 3 - 450006];
/* Computing 2nd power */
    r__11 = pb_1.ppion[(*nnn + 1 + *irun * 150001) * 3 - 450005];
/* Computing 2nd power */
    r__12 = pb_1.ppion[(*nnn + 1 + *irun * 150001) * 3 - 450004];
    devio = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4) + 
	    sqrt(r__5 * r__5 + r__6 * r__6 + r__7 * r__7 + r__8 * r__8) + 
	    sqrt(r__9 * r__9 + r__10 * r__10 + r__11 * r__11 + r__12 * r__12) 
	    - leadng_1.e1;
/*        if(abs(devio).gt.0.02) write(93,*) 'decay2(): nt=',nt,devio,lb1 */
/*     add decay time to daughter's formation time at the last timestep: */
    if (*nt == input2_1.ntmax) {
	tau0 = .19733f / *wid;
	taudcy = tau0 * -1.f * log(1.f - ranart_(&rndf77_1.nseed));
/*     lorentz boost: */
	taudcy = taudcy * leadng_1.e1 / leadng_1.em1;
	leadng_1.tfnl += taudcy;
	leadng_1.xfnl += leadng_1.px1 / leadng_1.e1 * taudcy;
	leadng_1.yfnl += leadng_1.py1 / leadng_1.e1 * taudcy;
	leadng_1.zfnl += leadng_1.pz1 / leadng_1.e1 * taudcy;
	aa_1.r__[*i__ * 3 - 3] = leadng_1.xfnl;
	aa_1.r__[*i__ * 3 - 2] = leadng_1.yfnl;
	aa_1.r__[*i__ * 3 - 1] = leadng_1.zfnl;
	tdecay_1.tfdcy[*i__ - 1] = leadng_1.tfnl;
	pa_1.rpion[(*nnn + *irun * 150001) * 3 - 450006] = leadng_1.xfnl;
	pa_1.rpion[(*nnn + *irun * 150001) * 3 - 450005] = leadng_1.yfnl;
	pa_1.rpion[(*nnn + *irun * 150001) * 3 - 450004] = leadng_1.zfnl;
	tdecay_1.tfdpi[*nnn + *irun * 150001 - 150002] = leadng_1.tfnl;
	pa_1.rpion[(*nnn + 1 + *irun * 150001) * 3 - 450006] = leadng_1.xfnl;
	pa_1.rpion[(*nnn + 1 + *irun * 150001) * 3 - 450005] = leadng_1.yfnl;
	pa_1.rpion[(*nnn + 1 + *irun * 150001) * 3 - 450004] = leadng_1.zfnl;
	tdecay_1.tfdpi[*nnn + 1 + *irun * 150001 - 150002] = leadng_1.tfnl;
    }
/* L200: */
/* L210: */
/* L220: */
    return 0;
} /* dkine2_ */

/* --------------------------------------------------------------------------- */
/* --------------------------------------------------------------------------- */
/* PURPOSE : CALCULATE THE MASS AND MOMENTUM OF BARYON RESONANCE */
/*           AFTER PION OR ETA BEING ABSORBED BY A NUCLEON */
/* NOTE    : */

/* DATE    : JAN.29,1990 */
/* Subroutine */ int dreson_(integer *i1, integer *i2)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static integer i__;
    static doublereal p1, p2, p3, e10, e20;
    static real dm;
    static doublereal scheck;

    /* Fortran I/O blocks */
    static cilist io___1108 = { 0, 99, 0, 0, 0 };
    static cilist io___1109 = { 0, 99, 0, 0, 0 };
    static cilist io___1110 = { 0, 99, 0, 0, 0 };
    static cilist io___1111 = { 0, 99, 0, 0, 0 };


/* lin-9/2012: improve precision for argument in sqrt(): */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /RUN/ */
/* c      SAVE /PA/ */
/* c      SAVE /PB/ */
/* c      SAVE /PC/ */
/* c      SAVE /PD/ */
/* 1. DETERMINE THE MOMENTUM COMPONENT OF DELTA/N* IN THE LAB. FRAME */
/* lin-9/2012: improve precision for argument in sqrt(): */
/*        E10=SQRT(E(I1)**2+P(1,I1)**2+P(2,I1)**2+P(3,I1)**2) */
/*        E20=SQRT(E(I2)**2+P(1,I2)**2+P(2,I2)**2+P(3,I2)**2) */
/* Computing 2nd power */
    d__1 = (doublereal) cc_1.e[*i1 - 1];
/* Computing 2nd power */
    d__2 = (doublereal) bb_1.p[*i1 * 3 - 3];
/* Computing 2nd power */
    d__3 = (doublereal) bb_1.p[*i1 * 3 - 2];
/* Computing 2nd power */
    d__4 = (doublereal) bb_1.p[*i1 * 3 - 1];
    e10 = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4);
/* Computing 2nd power */
    d__1 = (doublereal) cc_1.e[*i2 - 1];
/* Computing 2nd power */
    d__2 = (doublereal) bb_1.p[*i2 * 3 - 3];
/* Computing 2nd power */
    d__3 = (doublereal) bb_1.p[*i2 * 3 - 2];
/* Computing 2nd power */
    d__4 = (doublereal) bb_1.p[*i2 * 3 - 1];
    e20 = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4);
    p1 = (doublereal) bb_1.p[*i1 * 3 - 3] + (doublereal) bb_1.p[*i2 * 3 - 3];
    p2 = (doublereal) bb_1.p[*i1 * 3 - 2] + (doublereal) bb_1.p[*i2 * 3 - 2];
    p3 = (doublereal) bb_1.p[*i1 * 3 - 1] + (doublereal) bb_1.p[*i2 * 3 - 1];
    if ((i__1 = ee_1.lb[*i2 - 1], abs(i__1)) == 1 || (i__2 = ee_1.lb[*i2 - 1],
	     abs(i__2)) == 2 || (i__3 = ee_1.lb[*i2 - 1], abs(i__3)) >= 6 && (
	    i__4 = ee_1.lb[*i2 - 1], abs(i__4)) <= 17) {
	cc_1.e[*i1 - 1] = 0.f;
	i__ = *i2;
    } else {
	cc_1.e[*i2 - 1] = 0.f;
	i__ = *i1;
    }
    bb_1.p[i__ * 3 - 3] = bb_1.p[*i1 * 3 - 3] + bb_1.p[*i2 * 3 - 3];
    bb_1.p[i__ * 3 - 2] = bb_1.p[*i1 * 3 - 2] + bb_1.p[*i2 * 3 - 2];
    bb_1.p[i__ * 3 - 1] = bb_1.p[*i1 * 3 - 1] + bb_1.p[*i2 * 3 - 1];
/* 2. DETERMINE THE MASS OF DELTA/N* BY USING THE REACTION KINEMATICS */
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    d__1 = e10 + e20;
/* Computing 2nd power */
    d__2 = p1;
/* Computing 2nd power */
    d__3 = p2;
/* Computing 2nd power */
    d__4 = p3;
    scheck = d__1 * d__1 - d__2 * d__2 - d__3 * d__3 - d__4 * d__4;
    if (scheck < 0.) {
	s_wsle(&io___1108);
	do_lio(&c__9, &c__1, "scheck17: ", (ftnlen)10);
	do_lio(&c__5, &c__1, (char *)&scheck, (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___1109);
	do_lio(&c__9, &c__1, "scheck17", (ftnlen)8);
	do_lio(&c__5, &c__1, (char *)&scheck, (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&e10, (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&e20, (ftnlen)sizeof(doublereal));
	do_lio(&c__4, &c__1, (char *)&bb_1.p[i__ * 3 - 3], (ftnlen)sizeof(
		real));
	do_lio(&c__4, &c__1, (char *)&bb_1.p[i__ * 3 - 2], (ftnlen)sizeof(
		real));
	do_lio(&c__4, &c__1, (char *)&bb_1.p[i__ * 3 - 1], (ftnlen)sizeof(
		real));
	e_wsle();
	s_wsle(&io___1110);
	do_lio(&c__9, &c__1, "scheck17-1", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&cc_1.e[*i1 - 1], (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&bb_1.p[*i1 * 3 - 3], (ftnlen)sizeof(
		real));
	do_lio(&c__4, &c__1, (char *)&bb_1.p[*i1 * 3 - 2], (ftnlen)sizeof(
		real));
	do_lio(&c__4, &c__1, (char *)&bb_1.p[*i1 * 3 - 1], (ftnlen)sizeof(
		real));
	e_wsle();
	s_wsle(&io___1111);
	do_lio(&c__9, &c__1, "scheck17-2", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&cc_1.e[*i2 - 1], (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&bb_1.p[*i2 * 3 - 3], (ftnlen)sizeof(
		real));
	do_lio(&c__4, &c__1, (char *)&bb_1.p[*i2 * 3 - 2], (ftnlen)sizeof(
		real));
	do_lio(&c__4, &c__1, (char *)&bb_1.p[*i2 * 3 - 1], (ftnlen)sizeof(
		real));
	e_wsle();
	scheck = 0.;
    }
    dm = sqrt((real) scheck);
/*        DM=SQRT((E10+E20)**2-P(1,I)**2-P(2,I)**2-P(3,I)**2) */
    cc_1.e[i__ - 1] = dm;
    return 0;
} /* dreson_ */

/* --------------------------------------------------------------------------- */
/* PURPOSE : CALCULATE THE MASS AND MOMENTUM OF RHO RESONANCE */
/*           AFTER PION + PION COLLISION */
/* DATE    : NOV. 30,1994 */
/* Subroutine */ int rhores_(integer *i1, integer *i2)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static doublereal p1, p2, p3, e10, e20;
    static real dm;
    static doublereal scheck;

    /* Fortran I/O blocks */
    static cilist io___1119 = { 0, 99, 0, 0, 0 };


/* lin-9/2012: improve precision for argument in sqrt(): */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /RUN/ */
/* c      SAVE /PA/ */
/* c      SAVE /PB/ */
/* c      SAVE /PC/ */
/* c      SAVE /PD/ */
/* 1. DETERMINE THE MOMENTUM COMPONENT OF THE RHO IN THE CMS OF NN FRAME */
/*    WE LET I1 TO BE THE RHO AND ABSORB I2 */
/* lin-9/2012: improve precision for argument in sqrt(): */
/*        E10=SQRT(E(I1)**2+P(1,I1)**2+P(2,I1)**2+P(3,I1)**2) */
/*        E20=SQRT(E(I2)**2+P(1,I2)**2+P(2,I2)**2+P(3,I2)**2) */
/* Computing 2nd power */
    d__1 = (doublereal) cc_1.e[*i1 - 1];
/* Computing 2nd power */
    d__2 = (doublereal) bb_1.p[*i1 * 3 - 3];
/* Computing 2nd power */
    d__3 = (doublereal) bb_1.p[*i1 * 3 - 2];
/* Computing 2nd power */
    d__4 = (doublereal) bb_1.p[*i1 * 3 - 1];
    e10 = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4);
/* Computing 2nd power */
    d__1 = (doublereal) cc_1.e[*i2 - 1];
/* Computing 2nd power */
    d__2 = (doublereal) bb_1.p[*i2 * 3 - 3];
/* Computing 2nd power */
    d__3 = (doublereal) bb_1.p[*i2 * 3 - 2];
/* Computing 2nd power */
    d__4 = (doublereal) bb_1.p[*i2 * 3 - 1];
    e20 = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4);
    p1 = (doublereal) bb_1.p[*i1 * 3 - 3] + (doublereal) bb_1.p[*i2 * 3 - 3];
    p2 = (doublereal) bb_1.p[*i1 * 3 - 2] + (doublereal) bb_1.p[*i2 * 3 - 2];
    p3 = (doublereal) bb_1.p[*i1 * 3 - 1] + (doublereal) bb_1.p[*i2 * 3 - 1];
    bb_1.p[*i1 * 3 - 3] += bb_1.p[*i2 * 3 - 3];
    bb_1.p[*i1 * 3 - 2] += bb_1.p[*i2 * 3 - 2];
    bb_1.p[*i1 * 3 - 1] += bb_1.p[*i2 * 3 - 1];
/* 2. DETERMINE THE MASS OF THE RHO BY USING THE REACTION KINEMATICS */
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    d__1 = e10 + e20;
/* Computing 2nd power */
    d__2 = p1;
/* Computing 2nd power */
    d__3 = p2;
/* Computing 2nd power */
    d__4 = p3;
    scheck = d__1 * d__1 - d__2 * d__2 - d__3 * d__3 - d__4 * d__4;
    if (scheck < 0.) {
	s_wsle(&io___1119);
	do_lio(&c__9, &c__1, "scheck18: ", (ftnlen)10);
	do_lio(&c__5, &c__1, (char *)&scheck, (ftnlen)sizeof(doublereal));
	e_wsle();
	scheck = 0.;
    }
    dm = sqrt((real) scheck);
/*        DM=SQRT((E10+E20)**2-P(1,I1)**2-P(2,I1)**2-P(3,I1)**2) */
    cc_1.e[*i1 - 1] = dm;
    cc_1.e[*i2 - 1] = 0.f;
    return 0;
} /* rhores_ */

/* --------------------------------------------------------------------------- */
/* PURPOSE : CALCULATE THE PION+NUCLEON CROSS SECTION ACCORDING TO THE */
/*           BREIT-WIGNER FORMULA/(p*)**2 */
/* VARIABLE : LA = 1 FOR DELTA RESONANCE */
/*            LA = 0 FOR N*(1440) RESONANCE */
/*            LA = 2 FRO N*(1535) RESONANCE */
/* DATE    : JAN.29,1990 */
doublereal xnpi_(integer *i1, integer *i2, integer *la, real *xmax)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4, r__5;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static real f1;
    static doublereal p1, p2, p3, e10, e20;
    static real dm, gam;
    extern doublereal w1440_(real *), w1535_(real *);
    static real avpi;
    extern doublereal width_(real *);
    static real pdelt2, pstar2;
    static doublereal scheck;
    static real avmass;

    /* Fortran I/O blocks */
    static cilist io___1129 = { 0, 99, 0, 0, 0 };


/* lin-9/2012: improve precision for argument in sqrt(): */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /RUN/ */
/* c      SAVE /PA/ */
/* c      SAVE /PB/ */
/* c      SAVE /PC/ */
/* c      SAVE /PD/ */
    avmass = .93886849999999999f;
    avpi = .13803333333333334f;
/* 1. DETERMINE THE MOMENTUM COMPONENT OF DELTA IN THE LAB. FRAME */
/* lin-9/2012: improve precision for argument in sqrt(): */
/*        E10=SQRT(E(I1)**2+P(1,I1)**2+P(2,I1)**2+P(3,I1)**2) */
/*        E20=SQRT(E(I2)**2+P(1,I2)**2+P(2,I2)**2+P(3,I2)**2) */
/* Computing 2nd power */
    d__1 = (doublereal) cc_1.e[*i1 - 1];
/* Computing 2nd power */
    d__2 = (doublereal) bb_1.p[*i1 * 3 - 3];
/* Computing 2nd power */
    d__3 = (doublereal) bb_1.p[*i1 * 3 - 2];
/* Computing 2nd power */
    d__4 = (doublereal) bb_1.p[*i1 * 3 - 1];
    e10 = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4);
/* Computing 2nd power */
    d__1 = (doublereal) cc_1.e[*i2 - 1];
/* Computing 2nd power */
    d__2 = (doublereal) bb_1.p[*i2 * 3 - 3];
/* Computing 2nd power */
    d__3 = (doublereal) bb_1.p[*i2 * 3 - 2];
/* Computing 2nd power */
    d__4 = (doublereal) bb_1.p[*i2 * 3 - 1];
    e20 = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4);
/*        P1=P(1,I1)+P(1,I2) */
/*        P2=P(2,I1)+P(2,I2) */
/*        P3=P(3,I1)+P(3,I2) */
    p1 = (doublereal) bb_1.p[*i1 * 3 - 3] + (doublereal) bb_1.p[*i2 * 3 - 3];
    p2 = (doublereal) bb_1.p[*i1 * 3 - 2] + (doublereal) bb_1.p[*i2 * 3 - 2];
    p3 = (doublereal) bb_1.p[*i1 * 3 - 1] + (doublereal) bb_1.p[*i2 * 3 - 1];
/* 2. DETERMINE THE MASS OF DELTA BY USING OF THE REACTION KINEMATICS */
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    d__1 = e10 + e20;
/* Computing 2nd power */
    d__2 = p1;
/* Computing 2nd power */
    d__3 = p2;
/* Computing 2nd power */
    d__4 = p3;
    scheck = d__1 * d__1 - d__2 * d__2 - d__3 * d__3 - d__4 * d__4;
    if (scheck < 0.) {
	s_wsle(&io___1129);
	do_lio(&c__9, &c__1, "scheck19: ", (ftnlen)10);
	do_lio(&c__5, &c__1, (char *)&scheck, (ftnlen)sizeof(doublereal));
	e_wsle();
	scheck = 0.;
    }
    dm = sqrt((real) scheck);
/*        DM=SQRT((E10+E20)**2-P1**2-P2**2-P3**2) */
    if (dm <= 1.1f) {
	ret_val = 1e-9f;
	return ret_val;
    }
/* 3. DETERMINE THE PION+NUCLEON CROSS SECTION ACCORDING TO THE */
/*    BREIT-WIGNER FORMULA IN UNIT OF FM**2 */
    if (*la == 1) {
	gam = width_(&dm);
/* Computing 2nd power */
	r__1 = gam;
/* Computing 2nd power */
	r__2 = gam;
/* Computing 2nd power */
	r__3 = dm - 1.232f;
	f1 = r__1 * r__1 * .25f / (r__2 * r__2 * .25f + r__3 * r__3);
	pdelt2 = .051622f;
	goto L10;
    }
    if (*la == 0) {
	gam = w1440_(&dm);
/* Computing 2nd power */
	r__1 = gam;
/* Computing 2nd power */
	r__2 = gam;
/* Computing 2nd power */
	r__3 = dm - 1.44f;
	f1 = r__1 * r__1 * .25f / (r__2 * r__2 * .25f + r__3 * r__3);
	pdelt2 = .157897f;
	goto L10;
    }
    if (*la == 2) {
	gam = w1535_(&dm);
/* Computing 2nd power */
	r__1 = gam;
/* Computing 2nd power */
	r__2 = gam;
/* Computing 2nd power */
	r__3 = dm - 1.535f;
	f1 = r__1 * r__1 * .25f / (r__2 * r__2 * .25f + r__3 * r__3);
	pdelt2 = .2181f;
    }
L10:
/* Computing 2nd power */
    r__2 = dm;
/* Computing 2nd power */
    r__3 = avmass;
/* Computing 2nd power */
    r__4 = avpi;
/* Computing 2nd power */
    r__1 = (r__2 * r__2 - r__3 * r__3 + r__4 * r__4) / (dm * 2.f);
/* Computing 2nd power */
    r__5 = avpi;
    pstar2 = r__1 * r__1 - r__5 * r__5;
    if (pstar2 <= 0.f) {
	ret_val = 1e-9f;
    } else {
/* give the cross section in unit of fm**2 */
	ret_val = f1 * (pdelt2 / pstar2) * *xmax / 10.f;
    }
    return ret_val;
} /* xnpi_ */

/* ------------------------------------------------------------------------------ */
/* **************************************** */
doublereal sigma_(real *srt, integer *id, integer *ioi, integer *iof)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4, r__5, r__6, r__7;
    doublereal d__1, d__2;

    /* Builtin functions */
    double atan(doublereal), log(doublereal), sqrt(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static real q, s, t, p0, q0, p2, s0, t0, q2, p02, q02, pr, ss, am0, pr2, 
	    ss0, alfa, beta, deln, amass, zplus, amass0, scheck, zminus;

    /* Fortran I/O blocks */
    static cilist io___1158 = { 0, 99, 0, 0, 0 };


/* PURPOSE : THIS IS THE PROGRAM TO CALCULATE THE ISOSPIN DECOMPOSED CROSS */
/*       SECTION BY USING OF B.J.VerWEST AND R.A.ARNDT'S PARAMETERIZATION */
/* REFERENCE: PHYS. REV. C25(1982)1979 */
/* QUANTITIES: IOI -- INITIAL ISOSPIN OF THE TWO NUCLEON SYSTEM */
/*            IOF -- FINAL   ISOSPIN ------------------------- */
/*            ID -- =1 FOR DELTA RESORANCE */
/*                  =2 FOR N*    RESORANCE */
/* DATE : MAY 15,1990 */
/* **************************************** */
    if (*id == 1) {
	amass0 = 1.22f;
	t0 = .12f;
    } else {
	amass0 = 1.43f;
	t0 = .2f;
    }
    if (*ioi == 1 && *iof == 1) {
	alfa = 3.772f;
	beta = 1.262f;
	am0 = 1.188f;
	t = .09902f;
    }
    if (*ioi == 1 && *iof == 0) {
	alfa = 15.28f;
	beta = 0.f;
	am0 = 1.245f;
	t = .1374f;
    }
    if (*ioi == 0 && *iof == 1) {
	alfa = 146.3f;
	beta = 0.f;
	am0 = 1.472f;
	t = .02649f;
    }
    zplus = (*srt - .9383f - amass0) * 2.f / t0;
    zminus = (1.0767f - amass0) * 2.f / t0;
    deln = atan(zplus) - atan(zminus);
    if (deln == 0.f) {
	deln = 1e-6f;
    }
/* Computing 2nd power */
    r__1 = zplus;
/* Computing 2nd power */
    r__2 = zminus;
    amass = amass0 + t0 / 4.f * log((r__1 * r__1 + 1.f) / (r__2 * r__2 + 1.f))
	     / deln;
/* Computing 2nd power */
    r__1 = *srt;
    s = r__1 * r__1;
    p2 = s / 4.f - .88040689000000005f;
/* Computing 2nd power */
    r__1 = am0 + .9383f;
    s0 = r__1 * r__1;
    p02 = s0 / 4.f - .88040689000000005f;
    p0 = sqrt(p02);
/* Computing 2nd power */
    r__1 = .9383f - amass;
/* Computing 2nd power */
    r__2 = amass + .9383f;
    pr2 = (s - r__1 * r__1) * (s - r__2 * r__2) / (s * 4.f);
    if (pr2 > 1e-6f) {
	pr = sqrt(pr2);
    } else {
	pr = 0.f;
	ret_val = 1e-6f;
	return ret_val;
    }
/* Computing 2nd power */
    r__1 = amass;
    ss = r__1 * r__1;
    q2 = (ss - .63984001000000013f) * (ss - 1.1592828900000001f) / (ss * 4.f);
    if (q2 > 1e-6f) {
	q = sqrt(q2);
    } else {
	q = 0.f;
	ret_val = 1e-6f;
	return ret_val;
    }
/* Computing 2nd power */
    r__1 = am0;
    ss0 = r__1 * r__1;
    q02 = (ss0 - .63984001000000013f) * (ss0 - 1.1592828900000001f) / (ss0 * 
	    4.f);
/* lin-9/2012: check argument in sqrt(): */
    scheck = q02;
    if (scheck < 0.f) {
	s_wsle(&io___1158);
	do_lio(&c__9, &c__1, "scheck20: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	scheck = 0.f;
    }
    q0 = sqrt(scheck);
/*        Q0=SQRT(Q02) */
    d__1 = (doublereal) (pr / p0);
    d__2 = (doublereal) beta;
/* Computing 2nd power */
    r__1 = am0;
/* Computing 2nd power */
    r__2 = t;
/* Computing 3rd power */
    r__3 = q / q0;
/* Computing 2nd power */
    r__5 = am0;
/* Computing 2nd power */
    r__4 = ss - r__5 * r__5;
/* Computing 2nd power */
    r__6 = am0;
/* Computing 2nd power */
    r__7 = t;
    ret_val = .12233087920268615f / (p2 * 2.f) * alfa * pow_dd(&d__1, &d__2) *
	     (r__1 * r__1) * (r__2 * r__2) * (r__3 * (r__3 * r__3)) / (r__4 * 
	    r__4 + r__6 * r__6 * (r__7 * r__7));
    ret_val *= 10.f;
    if (ret_val == 0.f) {
	ret_val = 1e-6f;
    }
    return ret_val;
} /* sigma_ */

/* **************************** */
doublereal denom_(real *srt, real *con)
{
    /* System generated locals */
    integer i__1;
    real ret_val, r__1, r__2, r__3, r__4, r__5, r__6;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real f;
    static integer i__;
    static real q, s, a1, p0, p1, q2, dm, tq, am0, amn, amp, sum, amin, amax, 
	    avpi;
    static integer nmax;
    static real dmass;

/* NOTE: CON=1 FOR DELTA RESONANCE, CON=2 FOR N*(1440) RESONANCE */
/*       con=-1 for N*(1535) */
/* PURPOSE : CALCULATE THE INTEGRAL IN THE DETAILED BALANCE */

/* DATE : NOV. 15, 1991 */
/* ****************************** */
    avpi = .13803333333333334f;
    am0 = 1.232f;
    amn = .9383f;
    amp = avpi;
    amax = *srt - .9383f;
    amin = avpi + .9383f;
    nmax = 200;
    dmass = (amax - amin) / (real) nmax;
    sum = 0.f;
    i__1 = nmax + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dm = amin + (real) (i__ - 1) * dmass;
	if (*con == 1.f) {
/* Computing 2nd power */
	    r__2 = dm;
/* Computing 2nd power */
	    r__3 = amn;
/* Computing 2nd power */
	    r__4 = amp;
/* Computing 2nd power */
	    r__1 = (r__2 * r__2 - r__3 * r__3 + r__4 * r__4) / (dm * 2.f);
/* Computing 2nd power */
	    r__5 = amp;
	    q2 = r__1 * r__1 - r__5 * r__5;
	    if (q2 > 0.f) {
		q = sqrt(q2);
	    } else {
		q = 1e-6f;
	    }
/* Computing 3rd power */
	    r__1 = q;
/* Computing 2nd power */
	    r__2 = amp;
/* Computing 2nd power */
	    r__3 = q / amp;
	    tq = r__1 * (r__1 * r__1) * .47f / (r__2 * r__2 * (r__3 * r__3 * 
		    .6f + 1.f));
	} else if (*con == 2.f) {
	    tq = .2f;
	    am0 = 1.44f;
	} else if (*con == -1.f) {
	    tq = .1f;
	    am0 = 1.535f;
	}
/* Computing 2nd power */
	r__1 = am0;
/* Computing 2nd power */
	r__2 = am0;
/* Computing 2nd power */
	r__3 = tq;
/* Computing 2nd power */
	r__5 = dm;
/* Computing 2nd power */
	r__6 = am0;
/* Computing 2nd power */
	r__4 = r__5 * r__5 - r__6 * r__6;
	a1 = tq * 4.f * (r__1 * r__1) / (r__2 * r__2 * (r__3 * r__3) + r__4 * 
		r__4);
/* Computing 2nd power */
	r__1 = *srt;
	s = r__1 * r__1;
/* Computing 2nd power */
	r__2 = dm;
/* Computing 2nd power */
	r__3 = amn;
/* Computing 2nd power */
	r__1 = s + r__2 * r__2 - r__3 * r__3;
/* Computing 2nd power */
	r__4 = dm;
	p0 = r__1 * r__1 / (s * 4.f) - r__4 * r__4;
	if (p0 <= 0.f) {
	    p1 = 1e-6f;
	} else {
	    p1 = sqrt(p0);
	}
	f = dm * a1 * p1;
	if (i__ == 1 || i__ == nmax + 1) {
	    sum += f * .5f;
	} else {
	    sum += f;
	}
/* L10: */
    }
    ret_val = sum * dmass / 6.2831852000000001f;
    return ret_val;
} /* denom_ */

/* ********************************* */
/* subroutine : ang.FOR */
/* PURPOSE : Calculate the angular distribution of Delta production process */
/* DATE    : Nov. 19, 1992 */
/* REFERENCE: G. WOLF ET. AL., NUCL. PHYS. A517 (1990) 615 */
/* Note: this function applies when srt is larger than 2.14 GeV, */
/* for less energetic reactions, we assume the angular distribution */
/* is isotropic. */
/* ********************************** */
doublereal ang_(real *srt, integer *iseed)
{
    /* System generated locals */
    real ret_val, r__1, r__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static real p, q, x, b1s, b2s, ang1, ang2;
    extern doublereal ranart_(integer *);

/* c      SAVE /RNDF77/ */
/*        if(srt.le.2.14)then */
/*       b1s=0.5 */
/*       b2s=0. */
/*      endif */
    if (*srt > 2.14f && *srt <= 2.4f) {
/* Computing 2nd power */
	r__1 = *srt;
	b1s = 29.03f - *srt * 23.75f + r__1 * r__1 * 4.865f;
/* Computing 2nd power */
	r__1 = *srt;
	b2s = *srt * 25.53f - 30.33f - r__1 * r__1 * 5.301f;
    }
    if (*srt > 2.4f) {
	b1s = .06f;
	b2s = .4f;
    }
    x = ranart_(&rndf77_1.nseed);
    p = b1s / b2s;
    q = (x * 2.f - 1.f) * (b1s + b2s) / b2s;
/* Computing 2nd power */
    r__1 = q / 2.f;
/* Computing 3rd power */
    r__2 = p / 3.f;
    if (-q / 2.f + sqrt(r__1 * r__1 + r__2 * (r__2 * r__2)) >= 0.f) {
/* Computing 2nd power */
	r__1 = q / 2.f;
/* Computing 3rd power */
	r__2 = p / 3.f;
	d__1 = (doublereal) (-q / 2.f + sqrt(r__1 * r__1 + r__2 * (r__2 * 
		r__2)));
	ang1 = pow_dd(&d__1, &c_b5);
    } else {
/* Computing 2nd power */
	r__1 = q / 2.f;
/* Computing 3rd power */
	r__2 = p / 3.f;
	d__1 = (doublereal) (q / 2.f - sqrt(r__1 * r__1 + r__2 * (r__2 * r__2)
		));
	ang1 = -pow_dd(&d__1, &c_b5);
    }
/* Computing 2nd power */
    r__1 = q / 2.f;
/* Computing 3rd power */
    r__2 = p / 3.f;
    if (-q / 2.f - sqrt(r__1 * r__1 + r__2 * (r__2 * r__2)) >= 0.f) {
/* Computing 2nd power */
	r__1 = q / 2.f;
/* Computing 3rd power */
	r__2 = p / 3.f;
	d__1 = (doublereal) (-q / 2.f - sqrt(r__1 * r__1 + r__2 * (r__2 * 
		r__2)));
	ang2 = pow_dd(&d__1, &c_b5);
    } else {
/* Computing 2nd power */
	r__1 = q / 2.f;
/* Computing 3rd power */
	r__2 = p / 3.f;
	d__1 = (doublereal) (q / 2.f + sqrt(r__1 * r__1 + r__2 * (r__2 * r__2)
		));
	ang2 = -pow_dd(&d__1, &c_b5);
    }
    ret_val = ang1 + ang2;
    return ret_val;
} /* ang_ */

/* -------------------------------------------------------------------------- */
/* ****subprogram * kaon production from pi+B collisions ******************* */
doublereal pnlka_(real *srt)
{
    /* System generated locals */
    real ret_val;

    /* Local variables */
    static real t1, aka, ala, ana, sbbk;

/* units: fm**2 */
/* **********************************C */
    ala = 1.116f;
    aka = .498f;
    ana = .939f;
    t1 = ala + aka;
    if (*srt <= t1) {
	ret_val = 0.f;
    } else {
	if (*srt < 1.7f) {
	    sbbk = (*srt - t1) * 9.8901098901098905f;
	}
	if (*srt >= 1.7f) {
	    sbbk = .09f / (*srt - 1.6f);
	}
	ret_val = sbbk * .25f;
/* give the cross section in units of fm**2 */
	ret_val /= 10.f;
    }
    return ret_val;
} /* pnlka_ */

/* ------------------------------------------------------------------------- */
/* ****subprogram * kaon production from pi+B collisions ******************* */
doublereal pnska_(real *srt)
{
    /* System generated locals */
    real ret_val;

    /* Local variables */
    static real t1, aka, ala, ana, asa, sbb1, sbb2;

/* ********************************** */
    if (*srt > 3.f) {
	ret_val = 0.f;
	return ret_val;
    }
    ala = 1.116f;
    aka = .498f;
    ana = .939f;
    asa = 1.197f;
    t1 = asa + aka;
    if (*srt <= t1) {
	ret_val = 0.f;
	return ret_val;
    }
    if (*srt < 1.9f) {
	sbb1 = (*srt - t1) * 3.2110091743119265f;
    }
    if (*srt >= 1.9f) {
	sbb1 = .14f / (*srt - 1.7f);
    }
    sbb2 = 0.f;
    if (*srt > 1.682f) {
	sbb2 = (1.f - (*srt - 1.682f) * .75f) * .5f;
    }
    ret_val = (sbb1 + sbb2) * .25f;
/* give the cross section in fm**2 */
    ret_val /= 10.f;
    return ret_val;
} /* pnska_ */

/* ******************************* */

/*       Kaon momentum distribution in baryon-baryon-->N lamda K process */

/*       NOTE: dsima/dp is prototional to (1-p/p_max)(p/p_max)^2 */
/*              we use rejection method to generate kaon momentum */

/*       Variables: Fkaon = F(p)/F_max */
/*                 srt   = cms energy of the colliding pair, */
/*                          used to calculate the P_max */
/*       Date: Feb. 8, 1994 */

/*       Reference: C. M. Ko et al. */
/* ******************************* */
doublereal fkaon_(real *p, real *pmax)
{
    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real fmax;

    fmax = .148f;
    if (*pmax == 0.f) {
	*pmax = 1e-6f;
    }
/* Computing 2nd power */
    r__1 = *p / *pmax;
    ret_val = (1.f - *p / *pmax) * (r__1 * r__1);
    if (ret_val > fmax) {
	ret_val = fmax;
    }
    ret_val /= fmax;
    return ret_val;
} /* fkaon_ */

/* ************************ */
/* cross section for N*(1535) production in ND OR NN* collisions */
/* VARIABLES: */
/* LB1,LB2 ARE THE LABLES OF THE TWO COLLIDING PARTICLES */
/* SRT IS THE CMS ENERGY */
/* X1535 IS THE N*(1535) PRODUCTION CROSS SECTION */
/* NOTE THAT THE N*(1535) PRODUCTION CROSS SECTION IS 2 TIMES THE ETA */
/* PRODUCTION CROSS SECTION */
/* DATE: MAY 18, 1994 */
/* *********************** */
/* Subroutine */ int m1535_(integer *lb1, integer *lb2, real *srt, real *
	x1535)
{
    /* System generated locals */
    real r__1;

    /* Local variables */
    static real s0, sigma;

    s0 = 2.424f;
    *x1535 = 0.f;
    if (*srt <= s0) {
	return 0;
    }
/* Computing 2nd power */
    r__1 = *srt - s0;
    sigma = (*srt - s0) * .20399999999999999f / (r__1 * r__1 + .058f);
/* I N*(1535) PRODUCTION IN NUCLEON-DELTA COLLISIONS */
/* (1) nD(++)->pN*(+)(1535), pD(-)->nN*(0)(1535),pD(+)-->N*(+)p */
/* bz11/25/98 */
/*       IF((LB1*LB2.EQ.18).OR.(LB1*LB2.EQ.6). */
/*     1  or.(lb1*lb2).eq.8)then */
    if (*lb1 * *lb2 == 18 && (*lb1 == 2 || *lb2 == 2) || *lb1 * *lb2 == 6 && (
	    *lb1 == 1 || *lb2 == 1) || *lb1 * *lb2 == 8 && (*lb1 == 1 || *lb2 
	    == 1)) {
/* bz11/25/98end */
	*x1535 = sigma;
	return 0;
    }
/* (2) pD(0)->pN*(0)(1535),pD(0)->nN*(+)(1535) */
    if (*lb1 * *lb2 == 7) {
	*x1535 = sigma * 3.f;
	return 0;
    }
/* II N*(1535) PRODUCTION IN N*(1440)+NUCLEON REACTIONS */
/* (3) N*(+)(1440)p->N*(0+)(1535)p, N*(0)(1440)n->N*(0)(1535) */
/* bz11/25/98 */
/*       IF((LB1*LB2.EQ.11).OR.(LB1*LB2.EQ.20))THEN */
    if (*lb1 * *lb2 == 11 || *lb1 * *lb2 == 20 && (*lb1 == 2 || *lb2 == 2)) {
/* bz11/25/98end */
	*x1535 = sigma;
	return 0;
    }
/* (4) N*(0)(1440)p->N*(0+) or N*(+)(1440)n->N*(0+)(1535) */
/* bz11/25/98 */
/*       IF((LB1*LB2.EQ.10).OR.(LB1*LB2.EQ.22))X1535=3.*SIGMA */
    if (*lb1 * *lb2 == 10 && (*lb1 == 1 || *lb2 == 1) || *lb1 * *lb2 == 22 && 
	    (*lb1 == 2 || *lb2 == 2)) {
	*x1535 = sigma * 3.f;
    }
/* bz11/25/98end */
    return 0;
} /* m1535_ */

/* ************************ */
/* cross section for N*(1535) production in NN collisions */
/* VARIABLES: */
/* LB1,LB2 ARE THE LABLES OF THE TWO COLLIDING PARTICLES */
/* SRT IS THE CMS ENERGY */
/* X1535 IS THE N*(1535) PRODUCTION CROSS SECTION */
/* NOTE THAT THE N*(1535) PRODUCTION CROSS SECTION IS 2 TIMES THE ETA */
/* PRODUCTION CROSS SECTION */
/* DATE: MAY 18, 1994 */
/* *********************** */
/* Subroutine */ int n1535_(integer *lb1, integer *lb2, real *srt, real *
	x1535)
{
    /* System generated locals */
    real r__1;

    /* Local variables */
    static real s0, sigma;

    s0 = 2.424f;
    *x1535 = 0.f;
    if (*srt <= s0) {
	return 0;
    }
/* Computing 2nd power */
    r__1 = *srt - s0;
    sigma = (*srt - s0) * .20399999999999999f / (r__1 * r__1 + .058f);
/* I N*(1535) PRODUCTION IN NUCLEON-NUCLEON COLLISIONS */
/* (1) pp->pN*(+)(1535), nn->nN*(0)(1535) */
/* bdbg11/25/98 */
/*       IF((LB1*LB2.EQ.1).OR.(LB1*LB2.EQ.4))then */
    if (*lb1 * *lb2 == 1 || *lb1 == 2 && *lb2 == 2) {
/* bz11/25/98end */
	*x1535 = sigma;
	return 0;
    }
/* (2) pn->pN*(0)(1535),pn->nN*(+)(1535) */
    if (*lb1 * *lb2 == 2) {
	*x1535 = sigma * 3.f;
	return 0;
    }
/* III N*(1535) PRODUCTION IN DELTA+DELTA REACTIONS */
/* (5) D(++)+D(0), D(+)+D(+),D(+)+D(-),D(0)+D(0) */
/* bz11/25/98 */
/*       IF((LB1*LB2.EQ.63).OR.(LB1*LB2.EQ.64).OR.(LB1*LB2.EQ.48). */
/*     1  OR.(LB1*LB2.EQ.49))then */
    if (*lb1 * *lb2 == 63 && (*lb1 == 7 || *lb2 == 7) || *lb1 * *lb2 == 64 && 
	    (*lb1 == 8 || *lb2 == 8) || *lb1 * *lb2 == 48 && (*lb1 == 6 || *
	    lb2 == 6) || *lb1 * *lb2 == 49 && (*lb1 == 7 || *lb2 == 7)) {
/* bz11/25/98end */
	*x1535 = sigma;
	return 0;
    }
/* (6) D(++)+D(-),D(+)+D(0) */
/* bz11/25/98 */
/*       IF((LB1*LB2.EQ.54).OR.(LB1*LB2.EQ.56))then */
    if (*lb1 * *lb2 == 54 && (*lb1 == 6 || *lb2 == 6) || *lb1 * *lb2 == 56 && 
	    (*lb1 == 7 || *lb2 == 7)) {
/* bz11/25/98end */
	*x1535 = sigma * 3.f;
	return 0;
    }
/* IV N*(1535) PRODUCTION IN N*(1440)+N*(1440) REACTIONS */
/* bz11/25/98 */
/*       IF((LB1*LB2.EQ.100).OR.(LB1*LB2.EQ.11*11))X1535=SIGMA */
    if (*lb1 == 10 && *lb2 == 10 || *lb1 == 11 && *lb2 == 11) {
	*x1535 = sigma;
    }
/*       IF(LB1*LB2.EQ.110)X1535=3.*SIGMA */
    if (*lb1 * *lb2 == 110 && (*lb1 == 10 || *lb2 == 10)) {
	*x1535 = sigma * 3.f;
    }
/* bdbg11/25/98end */
    return 0;
} /* n1535_ */

/* *********************************** */
/* FUNCTION WA1(DMASS) GIVES THE A1 DECAY WIDTH */
/* Subroutine */ int wida1_(real *dmass, real *rhomp, real *wa1, integer *
	iseed)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real epi, qqp, qqp2, erho, coupa, epirho;
    extern doublereal rhomas_(real *, integer *);
    static real pimass, rhomax;
    static integer icount;


    pimass = .137265f;
    coupa = 14.8f;

    rhomax = *dmass - pimass - .02f;
    if (rhomax <= 0.f) {
	*rhomp = 0.f;
/*   !! no decay */
	*wa1 = -10.f;
    }
    icount = 0;
L711:
    *rhomp = rhomas_(&rhomax, iseed);
    ++icount;
    if (*dmass <= pimass + *rhomp) {
	if (icount <= 100) {
	    goto L711;
	} else {
	    *rhomp = 0.f;
/*   !! no decay */
	    *wa1 = -10.f;
	    return 0;
	}
    }
/* Computing 2nd power */
    r__1 = *dmass;
/* Computing 2nd power */
    r__2 = *rhomp + pimass;
/* Computing 2nd power */
    r__3 = *dmass;
/* Computing 2nd power */
    r__4 = *rhomp - pimass;
    qqp2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4);
    qqp = sqrt(qqp2) / (*dmass * 2.f);
/* Computing 2nd power */
    r__1 = pimass;
/* Computing 2nd power */
    r__2 = qqp;
    epi = sqrt(r__1 * r__1 + r__2 * r__2);
/* Computing 2nd power */
    r__1 = *rhomp;
/* Computing 2nd power */
    r__2 = qqp;
    erho = sqrt(r__1 * r__1 + r__2 * r__2);
/* Computing 2nd power */
    r__2 = qqp;
/* Computing 2nd power */
    r__1 = epi * erho + r__2 * r__2;
/* Computing 2nd power */
    r__3 = *rhomp;
/* Computing 2nd power */
    r__4 = epi;
    epirho = r__1 * r__1 * 2.f + r__3 * r__3 * (r__4 * r__4);
/* Computing 2nd power */
    r__1 = coupa;
/* Computing 2nd power */
    r__2 = *dmass;
    *wa1 = r__1 * r__1 * qqp * epirho / (r__2 * r__2 * 75.398399999999995f);
    return 0;
} /* wida1_ */

/* *********************************** */
/* FUNCTION W1535(DMASS) GIVES THE N*(1535) DECAY WIDTH */
/*     FOR A GIVEN N*(1535) MASS */
/* HERE THE FORMULA GIVEN BY KITAZOE IS USED */
doublereal w1535_(real *dmass)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real aux, qavail, avmass, pimass;

    avmass = .938868f;
    pimass = .137265f;
/* Computing 2nd power */
    r__2 = *dmass;
/* Computing 2nd power */
    r__3 = avmass;
/* Computing 2nd power */
    r__4 = pimass;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = avmass * pimass;
    aux = r__1 * r__1 * .25f - r__5 * r__5;
    if (aux > 0.f) {
/* Computing 2nd power */
	r__1 = *dmass;
	qavail = sqrt(aux / (r__1 * r__1));
    } else {
	qavail = 1e-6f;
    }
    ret_val = qavail * .15f / .467f;
/*       W1535=0.15 */
    return ret_val;
} /* w1535_ */

/* *********************************** */
/* FUNCTION W1440(DMASS) GIVES THE N*(1440) DECAY WIDTH */
/*     FOR A GIVEN N*(1535) MASS */
/* HERE THE FORMULA GIVEN BY KITAZOE IS USED */
doublereal w1440_(real *dmass)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real aux, qavail, avmass, pimass;

    avmass = .938868f;
    pimass = .137265f;
/* Computing 2nd power */
    r__2 = *dmass;
/* Computing 2nd power */
    r__3 = avmass;
/* Computing 2nd power */
    r__4 = pimass;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = avmass * pimass;
    aux = r__1 * r__1 * .25f - r__5 * r__5;
    if (aux > 0.f) {
	qavail = sqrt(aux) / *dmass;
    } else {
	qavail = 1e-6f;
    }
/*              w1440=0.2 */
/* Computing 3rd power */
    r__1 = qavail / .397f;
    ret_val = r__1 * (r__1 * r__1) * .2f;
    return ret_val;
} /* w1440_ */

/* *************** */
/* PURPOSE : CALCULATE THE PION(ETA)+NUCLEON CROSS SECTION */
/*           ACCORDING TO THE BREIT-WIGNER FORMULA, */
/*           NOTE THAT N*(1535) IS S_11 */
/* VARIABLE : LA = 1 FOR PI+N */
/*            LA = 0 FOR ETA+N */
/* DATE    : MAY 16, 1994 */
/* *************** */
doublereal xn1535_(integer *i1, integer *i2, integer *la)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static real f1;
    static doublereal p1, p2, p3, e10, e20;
    static real dm, gam;
    extern doublereal w1535_(real *);
    static real gam0, avpi, xmax;
    static doublereal scheck;
    static real avmass;

    /* Fortran I/O blocks */
    static cilist io___1228 = { 0, 99, 0, 0, 0 };


/* lin-9/2012: improve precision for argument in sqrt(): */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /RUN/ */
/* c      SAVE /PA/ */
/* c      SAVE /PB/ */
/* c      SAVE /PC/ */
/* c      SAVE /PD/ */
    avmass = .93886849999999999f;
    avpi = .13803333333333334f;
/* 1. DETERMINE THE MOMENTUM COMPONENT OF N*(1535) IN THE LAB. FRAME */
/* lin-9/2012: improve precision for argument in sqrt(): */
/*        E10=SQRT(E(I1)**2+P(1,I1)**2+P(2,I1)**2+P(3,I1)**2) */
/*        E20=SQRT(E(I2)**2+P(1,I2)**2+P(2,I2)**2+P(3,I2)**2) */
/* Computing 2nd power */
    d__1 = (doublereal) cc_1.e[*i1 - 1];
/* Computing 2nd power */
    d__2 = (doublereal) bb_1.p[*i1 * 3 - 3];
/* Computing 2nd power */
    d__3 = (doublereal) bb_1.p[*i1 * 3 - 2];
/* Computing 2nd power */
    d__4 = (doublereal) bb_1.p[*i1 * 3 - 1];
    e10 = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4);
/* Computing 2nd power */
    d__1 = (doublereal) cc_1.e[*i2 - 1];
/* Computing 2nd power */
    d__2 = (doublereal) bb_1.p[*i2 * 3 - 3];
/* Computing 2nd power */
    d__3 = (doublereal) bb_1.p[*i2 * 3 - 2];
/* Computing 2nd power */
    d__4 = (doublereal) bb_1.p[*i2 * 3 - 1];
    e20 = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4);
/*        P1=P(1,I1)+P(1,I2) */
/*        P2=P(2,I1)+P(2,I2) */
/*        P3=P(3,I1)+P(3,I2) */
    p1 = (doublereal) bb_1.p[*i1 * 3 - 3] + (doublereal) bb_1.p[*i2 * 3 - 3];
    p2 = (doublereal) bb_1.p[*i1 * 3 - 2] + (doublereal) bb_1.p[*i2 * 3 - 2];
    p3 = (doublereal) bb_1.p[*i1 * 3 - 1] + (doublereal) bb_1.p[*i2 * 3 - 1];
/* 2. DETERMINE THE MASS OF DELTA BY USING OF THE REACTION KINEMATICS */
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    d__1 = e10 + e20;
/* Computing 2nd power */
    d__2 = p1;
/* Computing 2nd power */
    d__3 = p2;
/* Computing 2nd power */
    d__4 = p3;
    scheck = d__1 * d__1 - d__2 * d__2 - d__3 * d__3 - d__4 * d__4;
    if (scheck < 0.) {
	s_wsle(&io___1228);
	do_lio(&c__9, &c__1, "scheck21: ", (ftnlen)10);
	do_lio(&c__5, &c__1, (char *)&scheck, (ftnlen)sizeof(doublereal));
	e_wsle();
	scheck = 0.;
    }
    dm = sqrt((real) scheck);
/*        DM=SQRT((E10+E20)**2-P1**2-P2**2-P3**2) */
    if (dm <= 1.1f) {
	ret_val = 1e-6f;
	return ret_val;
    }
/* 3. DETERMINE THE PION(ETA)+NUCLEON->N*(1535) CROSS SECTION ACCORDING TO THE */
/*    BREIT-WIGNER FORMULA IN UNIT OF FM**2 */
    gam = w1535_(&dm);
    gam0 = .15f;
/* Computing 2nd power */
    r__1 = gam0;
/* Computing 2nd power */
    r__2 = gam;
/* Computing 2nd power */
    r__3 = dm - 1.535f;
    f1 = r__1 * r__1 * .25f / (r__2 * r__2 * .25f + r__3 * r__3);
    if (*la == 1) {
	xmax = 11.3f;
    } else {
	xmax = 74.f;
    }
    ret_val = f1 * xmax / 10.f;
    return ret_val;
} /* xn1535_ */

/* **************************8 */
/* FUNCTION FDE(DMASS) GIVES DELTA MASS DISTRIBUTION BY USING OF */
/* KITAZOE'S FORMULA */
doublereal fdelta_(real *dmass)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3;

    /* Local variables */
    static real fd, am0, amn, avpi;
    extern doublereal width_(real *);

    amn = .938869f;
    avpi = .13803333f;
    am0 = 1.232f;
/* Computing 2nd power */
    r__1 = width_(dmass);
/* Computing 2nd power */
    r__2 = *dmass - 1.232f;
/* Computing 2nd power */
    r__3 = width_(dmass);
    fd = r__1 * r__1 * .25f / (r__2 * r__2 + r__3 * r__3 * .25f);
    ret_val = fd;
    return ret_val;
} /* fdelta_ */

/* FUNCTION WIDTH(DMASS) GIVES THE DELTA DECAY WIDTH FOR A GIVEN DELTA MASS */
/* HERE THE FORMULA GIVEN BY KITAZOE IS USED */
doublereal width_(real *dmass)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real aux, qavail, avmass, pimass;

    avmass = .938868f;
    pimass = .137265f;
/* Computing 2nd power */
    r__2 = *dmass;
/* Computing 2nd power */
    r__3 = avmass;
/* Computing 2nd power */
    r__4 = pimass;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = avmass * pimass;
    aux = r__1 * r__1 * .25f - r__5 * r__5;
    if (aux > 0.f) {
/* Computing 2nd power */
	r__1 = *dmass;
	qavail = sqrt(aux / (r__1 * r__1));
    } else {
	qavail = 1e-6f;
    }
/* Computing 3rd power */
    r__1 = qavail;
/* Computing 2nd power */
    r__2 = pimass;
/* Computing 2nd power */
    r__3 = qavail / pimass;
    ret_val = r__1 * (r__1 * r__1) * .47f / (r__2 * r__2 * (r__3 * r__3 * .6f 
	    + 1.f));
/*       width=0.115 */
    return ret_val;
} /* width_ */

/* *********************************** */
/* Subroutine */ int ddp2_(real *srt, integer *iseed, real *px, real *py, 
	real *pz, real *dm1, real *pnx, real *pny, real *pnz, real *dm2, real 
	*ppx, real *ppy, real *ppz, integer *icou1)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5, r__6;

    /* Builtin functions */
    double sqrt(doublereal), cos(doublereal), sin(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static real v, w, x, f00, ga, ek, en, ep, x00, pi, bx, by, bz, pn, pt, 
	    pn2, fai, amn, amp, eln, sig, pnt;
    extern doublereal ptr_(real *, integer *);
    static real srt1, fain, elnc, xmax, xptr, fmax00, pbeta, ptmax, pzmax;
    static integer ntrym, ntryx;
    static real trans0, ptmax2, pzmax2, scheck;
    extern /* Subroutine */ int rmasdd_(real *, real *, real *, real *, real *
	    , integer *, integer *, real *, real *);
    extern doublereal ranart_(integer *);
    static real xratio;

    /* Fortran I/O blocks */
    static cilist io___1268 = { 0, 99, 0, 0, 0 };
    static cilist io___1277 = { 0, 99, 0, 0, 0 };


/* PURPOSE : CALCULATE MOMENTUM OF PARTICLES IN THE FINAL SATAT FROM */
/* THE PROCESS N+N--->D1+D2+PION */
/*       DATE : July 25, 1994 */
/* Generate the masses and momentum for particles in the NN-->DDpi process */
/* for a given center of mass energy srt, the momenta are given in the center */
/* of mass of the NN */
/* **************************************** */
/* c      SAVE /TABLE/ */
/* c      SAVE /RNDF77/ */
    *icou1 = 0;
    pi = 3.1415926f;
    amn = .93892500000000001f;
    amp = .137265f;
/* (1) GENGRATE THE MASS OF DELTA1 AND DELTA2 USING */
    srt1 = *srt - amp - .02f;
    ntrym = 0;
L8:
    rmasdd_(&srt1, &c_b195, &c_b195, &c_b197, &c_b197, iseed, &c__1, dm1, dm2)
	    ;
    ++ntrym;
/* CONSTANTS FOR GENERATING THE LONGITUDINAL MOMENTUM */
/* FOR ONE OF THE RESONANCES */
    v = .43f;
    w = -.84f;
/* (2) Generate the transverse momentum */
/*     OF DELTA1 */
/* (2.1) estimate the maximum transverse momentum */
/* Computing 2nd power */
    r__1 = *srt;
/* Computing 2nd power */
    r__2 = *dm1 + *dm2 + amp;
/* Computing 2nd power */
    r__3 = *srt;
/* Computing 2nd power */
    r__4 = *dm1 - amp - *dm2;
/* Computing 2nd power */
    r__5 = *srt;
    ptmax2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / 4.f /
	     (r__5 * r__5);
    if (ptmax2 <= 0.f) {
	goto L8;
    }
    ptmax = sqrt(ptmax2) * 1.f / 3.f;
L7:
    pt = ptr_(&ptmax, iseed);
/* (3.1) THE MAXIMUM LONGITUDINAL MOMENTUM IS */
/* Computing 2nd power */
    r__1 = *srt;
/* Computing 2nd power */
    r__2 = *dm1 + *dm2 + amp;
/* Computing 2nd power */
    r__3 = *srt;
/* Computing 2nd power */
    r__4 = *dm1 - amp - *dm2;
/* Computing 2nd power */
    r__5 = *srt;
/* Computing 2nd power */
    r__6 = pt;
    pzmax2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / 4.f /
	     (r__5 * r__5) - r__6 * r__6;
    if (pzmax2 < 0.f && ntrym <= 100) {
	goto L7;
    } else {
	pzmax2 = 1e-9f;
    }
    pzmax = sqrt(pzmax2);
    xmax = pzmax * 2.f / *srt;
/* (3.2) THE GENERATED X IS */
/* THE DSTRIBUTION HAS A MAXIMUM AT X0=-V/(2*w), f(X0)=1.056 */
    ntryx = 0;
    fmax00 = 1.056f;
    x00 = .26f;
    if (dabs(xmax) > .26f) {
	f00 = fmax00;
    } else {
/* Computing 2nd power */
	r__1 = xmax;
	f00 = v * dabs(xmax) + 1.f + w * (r__1 * r__1);
    }
L9:
    x = xmax * (1.f - ranart_(&rndf77_1.nseed) * 2.f);
    ++ntryx;
/* Computing 2nd power */
    r__1 = x;
    xratio = (v * dabs(x) + 1.f + w * (r__1 * r__1)) / f00;
/* lin-8/17/00       IF(xratio.LT.RANART(NSEED).and.ntryx.le.50)GO TO 9 */
    if (xratio < ranart_(&rndf77_1.nseed) && ntryx <= 50) {
	goto L9;
    }
/* (3.5) THE PZ IS */
    *pz = *srt * .5f * x;
/* The x and y components of the deltA1 */
    fai = pi * 2.f * ranart_(&rndf77_1.nseed);
    *px = pt * cos(fai);
    *py = pt * sin(fai);
/* find the momentum of delta2 and pion */
/* the energy of the delta1 */
/* Computing 2nd power */
    r__1 = *dm1;
/* Computing 2nd power */
    r__2 = pt;
/* Computing 2nd power */
    r__3 = *pz;
    ek = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
/* (1) Generate the momentum of the delta2 in the cms of delta2 and pion */
/*     the energy of the cms of DP */
    eln = *srt - ek;
    if (eln <= 0.f) {
	*icou1 = -1;
	return 0;
    }
/* beta and gamma of the cms of delta2+pion */
    bx = -(*px) / eln;
    by = -(*py) / eln;
    bz = -(*pz) / eln;
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__1 = bx;
/* Computing 2nd power */
    r__2 = by;
/* Computing 2nd power */
    r__3 = bz;
    scheck = 1.f - r__1 * r__1 - r__2 * r__2 - r__3 * r__3;
    if (scheck <= 0.f) {
	s_wsle(&io___1268);
	do_lio(&c__9, &c__1, "scheck22: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	s_stop("", (ftnlen)0);
    }
    ga = 1.f / sqrt(scheck);
/*       ga=1./sqrt(1.-bx**2-by**2-bz**2) */
/* the momentum of delta2 and pion in their cms frame */
    elnc = eln / ga;
/* Computing 2nd power */
    r__2 = elnc;
/* Computing 2nd power */
    r__3 = *dm2;
/* Computing 2nd power */
    r__4 = amp;
/* Computing 2nd power */
    r__1 = (r__2 * r__2 + r__3 * r__3 - r__4 * r__4) / (elnc * 2.f);
/* Computing 2nd power */
    r__5 = *dm2;
    pn2 = r__1 * r__1 - r__5 * r__5;
    if (pn2 <= 0.f) {
	*icou1 = -1;
	return 0;
    }
    pn = sqrt(pn2);
/* lin-10/25/02 get rid of argument usage mismatch in PTR(): */
    xptr = pn * .33f;
/*       PNT=PTR(0.33*PN,ISEED) */
    pnt = ptr_(&xptr, iseed);
/* lin-10/25/02-end */
    fain = pi * 2.f * ranart_(&rndf77_1.nseed);
    *pnx = pnt * cos(fain);
    *pny = pnt * sin(fain);
    sig = 1.f;
    if (x > 0.f) {
	sig = -1.f;
    }
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__1 = pn;
/* Computing 2nd power */
    r__2 = pnt;
    scheck = r__1 * r__1 - r__2 * r__2;
    if (scheck < 0.f) {
	s_wsle(&io___1277);
	do_lio(&c__9, &c__1, "scheck23: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	scheck = 0.f;
    }
    *pnz = sig * sqrt(scheck);
/*       pnz=SIG*SQRT(pn**2-PNT**2) */
/* Computing 2nd power */
    r__1 = *dm2;
/* Computing 2nd power */
    r__2 = *pnx;
/* Computing 2nd power */
    r__3 = *pny;
/* Computing 2nd power */
    r__4 = *pnz;
    en = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
/* (2) the momentum for the pion */
    *ppx = -(*pnx);
    *ppy = -(*pny);
    *ppz = -(*pnz);
/* Computing 2nd power */
    r__1 = amp;
/* Computing 2nd power */
    r__2 = *ppx;
/* Computing 2nd power */
    r__3 = *ppy;
/* Computing 2nd power */
    r__4 = *ppz;
    ep = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
/* (3) for the delta2, LORENTZ-TRANSFORMATION INTO nn cms FRAME */
    pbeta = *pnx * bx + *pny * by + *pnz * bz;
    trans0 = ga * (ga * pbeta / (ga + 1.f) + en);
    *pnx = bx * trans0 + *pnx;
    *pny = by * trans0 + *pny;
    *pnz = bz * trans0 + *pnz;
/* (4) for the pion, LORENTZ-TRANSFORMATION INTO nn cms FRAME */
    if (ep == 0.f) {
	ep = 1e-9f;
    }
    pbeta = *ppx * bx + *ppy * by + *ppz * bz;
    trans0 = ga * (ga * pbeta / (ga + 1.f) + ep);
    *ppx = bx * trans0 + *ppx;
    *ppy = by * trans0 + *ppy;
    *ppz = bz * trans0 + *ppz;
    return 0;
} /* ddp2_ */

/* *************************************** */
/* Subroutine */ int ddrho_(real *srt, integer *iseed, real *px, real *py, 
	real *pz, real *dm1, real *pnx, real *pny, real *pnz, real *dm2, real 
	*ppx, real *ppy, real *ppz, real *amp, integer *icou1)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5, r__6;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double sqrt(doublereal), cos(doublereal), sin(doublereal);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static real v, w, x, f00, ga, ek, en, ep, x00, pi, bx, by, bz, pn, pt, 
	    pn2, fai, amn, eln, sig, pnt;
    extern doublereal ptr_(real *, integer *);
    static real srt1, fain, elnc, xmax, xptr, fmax00, pbeta, ptmax, pzmax;
    static integer ntrym, ntryx;
    static real trans0, ptmax2, pzmax2, scheck;
    extern /* Subroutine */ int rmasdd_(real *, real *, real *, real *, real *
	    , integer *, integer *, real *, real *);
    extern doublereal ranart_(integer *), rhomas_(real *, integer *);
    static real rhomax, xratio;

    /* Fortran I/O blocks */
    static cilist io___1291 = { 0, 99, 0, 0, 0 };
    static cilist io___1309 = { 0, 99, 0, 0, 0 };
    static cilist io___1318 = { 0, 99, 0, 0, 0 };


/* PURPOSE : CALCULATE MOMENTUM OF PARTICLES IN THE FINAL SATAT FROM */
/* THE PROCESS N+N--->D1+D2+rho */
/*       DATE : Nov.5, 1994 */
/* Generate the masses and momentum for particles in the NN-->DDrho process */
/* for a given center of mass energy srt, the momenta are given in the center */
/* of mass of the NN */
/* **************************************** */
/* c      SAVE /TABLE/ */
/* c      SAVE /RNDF77/ */
    *icou1 = 0;
    pi = 3.1415926f;
    amn = .93892500000000001f;
    *amp = .77000000000000002f;
/* (1) GENGRATE THE MASS OF DELTA1 AND DELTA2 USING */
    srt1 = *srt - *amp - .02f;
    ntrym = 0;
L8:
    rmasdd_(&srt1, &c_b195, &c_b195, &c_b197, &c_b197, iseed, &c__1, dm1, dm2)
	    ;
    ++ntrym;
/* GENERATE THE MASS FOR THE RHO */
    rhomax = *srt - *dm1 - *dm2 - .02f;
    if (rhomax <= 0.f && ntrym <= 20) {
	goto L8;
    }
    *amp = rhomas_(&rhomax, iseed);
/* CONSTANTS FOR GENERATING THE LONGITUDINAL MOMENTUM */
/* FOR ONE OF THE RESONANCES */
    v = .43f;
    w = -.84f;
/* (2) Generate the transverse momentum */
/*     OF DELTA1 */
/* (2.1) estimate the maximum transverse momentum */
/* Computing 2nd power */
    r__1 = *srt;
/* Computing 2nd power */
    r__2 = *dm1 + *dm2 + *amp;
/* Computing 2nd power */
    r__3 = *srt;
/* Computing 2nd power */
    r__4 = *dm1 - *amp - *dm2;
/* Computing 2nd power */
    r__5 = *srt;
    ptmax2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / 4.f /
	     (r__5 * r__5);
/* lin-9/2012: check argument in sqrt(): */
    scheck = ptmax2;
    if (scheck < 0.f) {
	s_wsle(&io___1291);
	do_lio(&c__9, &c__1, "scheck24: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	scheck = 0.f;
    }
    ptmax = sqrt(scheck) * 1.f / 3.f;
/*       PTMAX=SQRT(PTMAX2)*1./3. */
L7:
    pt = ptr_(&ptmax, iseed);
/* (3) GENGRATE THE LONGITUDINAL MOMENTUM FOR DM1 */
/*     USING THE GIVEN DISTRIBUTION */
/* (3.1) THE MAXIMUM LONGITUDINAL MOMENTUM IS */
/* Computing 2nd power */
    r__1 = *srt;
/* Computing 2nd power */
    r__2 = *dm1 + *dm2 + *amp;
/* Computing 2nd power */
    r__3 = *srt;
/* Computing 2nd power */
    r__4 = *dm1 - *amp - *dm2;
/* Computing 2nd power */
    r__5 = *srt;
/* Computing 2nd power */
    r__6 = pt;
    pzmax2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / 4.f /
	     (r__5 * r__5) - r__6 * r__6;
    if (pzmax2 < 0.f && ntrym <= 100) {
	goto L7;
    } else {
	pzmax2 = 1e-6f;
    }
    pzmax = sqrt(pzmax2);
    xmax = pzmax * 2.f / *srt;
/* (3.2) THE GENERATED X IS */
/* THE DSTRIBUTION HAS A MAXIMUM AT X0=-V/(2*w), f(X0)=1.056 */
    ntryx = 0;
    fmax00 = 1.056f;
    x00 = .26f;
    if (dabs(xmax) > .26f) {
	f00 = fmax00;
    } else {
/* Computing 2nd power */
	r__1 = xmax;
	f00 = v * dabs(xmax) + 1.f + w * (r__1 * r__1);
    }
L9:
    x = xmax * (1.f - ranart_(&rndf77_1.nseed) * 2.f);
    ++ntryx;
/* Computing 2nd power */
    r__1 = x;
    xratio = (v * dabs(x) + 1.f + w * (r__1 * r__1)) / f00;
/* lin-8/17/00       IF(xratio.LT.RANART(NSEED).and.ntryx.le.50)GO TO 9 */
    if (xratio < ranart_(&rndf77_1.nseed) && ntryx <= 50) {
	goto L9;
    }
/* (3.5) THE PZ IS */
    *pz = *srt * .5f * x;
/* The x and y components of the delta1 */
    fai = pi * 2.f * ranart_(&rndf77_1.nseed);
    *px = pt * cos(fai);
    *py = pt * sin(fai);
/* find the momentum of delta2 and rho */
/* the energy of the delta1 */
/* Computing 2nd power */
    r__1 = *dm1;
/* Computing 2nd power */
    r__2 = pt;
/* Computing 2nd power */
    r__3 = *pz;
    ek = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
/* (1) Generate the momentum of the delta2 in the cms of delta2 and rho */
/*     the energy of the cms of Drho */
    eln = *srt - ek;
    if (eln <= 0.f) {
	*icou1 = -1;
	return 0;
    }
/* beta and gamma of the cms of delta2 and rho */
    bx = -(*px) / eln;
    by = -(*py) / eln;
    bz = -(*pz) / eln;
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__1 = bx;
/* Computing 2nd power */
    r__2 = by;
/* Computing 2nd power */
    r__3 = bz;
    scheck = 1.f - r__1 * r__1 - r__2 * r__2 - r__3 * r__3;
    if (scheck <= 0.f) {
	s_wsle(&io___1309);
	do_lio(&c__9, &c__1, "scheck25: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	s_stop("", (ftnlen)0);
    }
    ga = 1.f / sqrt(scheck);
/*       ga=1./sqrt(1.-bx**2-by**2-bz**2) */
    elnc = eln / ga;
/* Computing 2nd power */
    r__2 = elnc;
/* Computing 2nd power */
    r__3 = *dm2;
/* Computing 2nd power */
    r__4 = *amp;
/* Computing 2nd power */
    r__1 = (r__2 * r__2 + r__3 * r__3 - r__4 * r__4) / (elnc * 2.f);
/* Computing 2nd power */
    r__5 = *dm2;
    pn2 = r__1 * r__1 - r__5 * r__5;
    if (pn2 <= 0.f) {
	*icou1 = -1;
	return 0;
    }
    pn = sqrt(pn2);
/* lin-10/25/02 get rid of argument usage mismatch in PTR(): */
    xptr = pn * .33f;
/*       PNT=PTR(0.33*PN,ISEED) */
    pnt = ptr_(&xptr, iseed);
/* lin-10/25/02-end */
    fain = pi * 2.f * ranart_(&rndf77_1.nseed);
    *pnx = pnt * cos(fain);
    *pny = pnt * sin(fain);
    sig = 1.f;
    if (x > 0.f) {
	sig = -1.f;
    }
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__1 = pn;
/* Computing 2nd power */
    r__2 = pnt;
    scheck = r__1 * r__1 - r__2 * r__2;
    if (scheck < 0.f) {
	s_wsle(&io___1318);
	do_lio(&c__9, &c__1, "scheck26: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	scheck = 0.f;
    }
    *pnz = sig * sqrt(scheck);
/*       pnz=SIG*SQRT(pn**2-PNT**2) */
/* Computing 2nd power */
    r__1 = *dm2;
/* Computing 2nd power */
    r__2 = *pnx;
/* Computing 2nd power */
    r__3 = *pny;
/* Computing 2nd power */
    r__4 = *pnz;
    en = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
/* (2) the momentum for the rho */
    *ppx = -(*pnx);
    *ppy = -(*pny);
    *ppz = -(*pnz);
/* Computing 2nd power */
    r__1 = *amp;
/* Computing 2nd power */
    r__2 = *ppx;
/* Computing 2nd power */
    r__3 = *ppy;
/* Computing 2nd power */
    r__4 = *ppz;
    ep = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
/* (3) for the delta2, LORENTZ-TRANSFORMATION INTO nn cms FRAME */
    pbeta = *pnx * bx + *pny * by + *pnz * bz;
    trans0 = ga * (ga * pbeta / (ga + 1.f) + en);
    *pnx = bx * trans0 + *pnx;
    *pny = by * trans0 + *pny;
    *pnz = bz * trans0 + *pnz;
/* (4) for the rho, LORENTZ-TRANSFORMATION INTO nn cms FRAME */
    if (ep == 0.f) {
	ep = 1e-9f;
    }
    pbeta = *ppx * bx + *ppy * by + *ppz * bz;
    trans0 = ga * (ga * pbeta / (ga + 1.f) + ep);
    *ppx = bx * trans0 + *ppx;
    *ppy = by * trans0 + *ppy;
    *ppz = bz * trans0 + *ppz;
    return 0;
} /* ddrho_ */

/* *************************************** */
/* Subroutine */ int pprho_(real *srt, integer *iseed, real *px, real *py, 
	real *pz, real *dm1, real *pnx, real *pny, real *pnz, real *dm2, real 
	*ppx, real *ppy, real *ppz, real *amp, integer *icou1)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5, r__6;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double sqrt(doublereal), cos(doublereal), sin(doublereal);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static real v, w, x, f00, ga, ek, en, ep, x00, pi, bx, by, bz, pn, pt, 
	    pn2, fai, amn, eln, sig, pnt;
    extern doublereal ptr_(real *, integer *);
    static real fain, elnc;
    static integer icou;
    static real xmax, xptr, fmax00, pbeta, ptmax, pzmax;
    static integer ntrym, ntryx;
    static real trans0, ptmax2, pzmax2, scheck;
    extern doublereal ranart_(integer *), rhomas_(real *, integer *);
    static real rhomax, xratio;

    /* Fortran I/O blocks */
    static cilist io___1332 = { 0, 99, 0, 0, 0 };
    static cilist io___1350 = { 0, 99, 0, 0, 0 };
    static cilist io___1359 = { 0, 99, 0, 0, 0 };


/* PURPOSE : CALCULATE MOMENTUM OF PARTICLES IN THE FINAL SATAT FROM */
/* THE PROCESS N+N--->N1+N2+rho */
/*       DATE : Nov.5, 1994 */
/* Generate the masses and momentum for particles in the NN--> process */
/* for a given center of mass energy srt, the momenta are given in the center */
/* of mass of the NN */
/* **************************************** */
/* c      SAVE /TABLE/ */
/* c      SAVE /RNDF77/ */
    ntrym = 0;
    *icou1 = 0;
    pi = 3.1415926f;
    amn = .93892500000000001f;
/*        AMP=770./1000. */
    *dm1 = amn;
    *dm2 = amn;
/* GENERATE THE MASS FOR THE RHO */
    rhomax = *srt - *dm1 - *dm2 - .02f;
    if (rhomax <= 0.f) {
	icou = -1;
	return 0;
    }
    *amp = rhomas_(&rhomax, iseed);
/* CONSTANTS FOR GENERATING THE LONGITUDINAL MOMENTUM */
/* FOR ONE OF THE nucleons */
    v = .43f;
    w = -.84f;
/* (2) Generate the transverse momentum */
/*     OF p1 */
/* (2.1) estimate the maximum transverse momentum */
/* Computing 2nd power */
    r__1 = *srt;
/* Computing 2nd power */
    r__2 = *dm1 + *dm2 + *amp;
/* Computing 2nd power */
    r__3 = *srt;
/* Computing 2nd power */
    r__4 = *dm1 - *amp - *dm2;
/* Computing 2nd power */
    r__5 = *srt;
    ptmax2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / 4.f /
	     (r__5 * r__5);
/* lin-9/2012: check argument in sqrt(): */
    scheck = ptmax2;
    if (scheck < 0.f) {
	s_wsle(&io___1332);
	do_lio(&c__9, &c__1, "scheck27: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	scheck = 0.f;
    }
    ptmax = sqrt(scheck) * 1.f / 3.f;
/*       PTMAX=SQRT(PTMAX2)*1./3. */
L7:
    pt = ptr_(&ptmax, iseed);
/* (3) GENGRATE THE LONGITUDINAL MOMENTUM FOR DM1 */
/*     USING THE GIVEN DISTRIBUTION */
/* (3.1) THE MAXIMUM LONGITUDINAL MOMENTUM IS */
/* Computing 2nd power */
    r__1 = *srt;
/* Computing 2nd power */
    r__2 = *dm1 + *dm2 + *amp;
/* Computing 2nd power */
    r__3 = *srt;
/* Computing 2nd power */
    r__4 = *dm1 - *amp - *dm2;
/* Computing 2nd power */
    r__5 = *srt;
/* Computing 2nd power */
    r__6 = pt;
    pzmax2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / 4.f /
	     (r__5 * r__5) - r__6 * r__6;
    ++ntrym;
    if (pzmax2 < 0.f && ntrym <= 100) {
	goto L7;
    } else {
	pzmax2 = 1e-6f;
    }
    pzmax = sqrt(pzmax2);
    xmax = pzmax * 2.f / *srt;
/* (3.2) THE GENERATED X IS */
/* THE DSTRIBUTION HAS A MAXIMUM AT X0=-V/(2*w), f(X0)=1.056 */
    ntryx = 0;
    fmax00 = 1.056f;
    x00 = .26f;
    if (dabs(xmax) > .26f) {
	f00 = fmax00;
    } else {
/* Computing 2nd power */
	r__1 = xmax;
	f00 = v * dabs(xmax) + 1.f + w * (r__1 * r__1);
    }
L9:
    x = xmax * (1.f - ranart_(&rndf77_1.nseed) * 2.f);
    ++ntryx;
/* Computing 2nd power */
    r__1 = x;
    xratio = (v * dabs(x) + 1.f + w * (r__1 * r__1)) / f00;
/* lin-8/17/00       IF(xratio.LT.RANART(NSEED).and.ntryx.le.50)GO TO 9 */
    if (xratio < ranart_(&rndf77_1.nseed) && ntryx <= 50) {
	goto L9;
    }
/* (3.5) THE PZ IS */
    *pz = *srt * .5f * x;
/* The x and y components of the delta1 */
    fai = pi * 2.f * ranart_(&rndf77_1.nseed);
    *px = pt * cos(fai);
    *py = pt * sin(fai);
/* find the momentum of delta2 and rho */
/* the energy of the delta1 */
/* Computing 2nd power */
    r__1 = *dm1;
/* Computing 2nd power */
    r__2 = pt;
/* Computing 2nd power */
    r__3 = *pz;
    ek = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
/* (1) Generate the momentum of the delta2 in the cms of delta2 and rho */
/*     the energy of the cms of Drho */
    eln = *srt - ek;
    if (eln <= 0.f) {
	*icou1 = -1;
	return 0;
    }
/* beta and gamma of the cms of the two partciles */
    bx = -(*px) / eln;
    by = -(*py) / eln;
    bz = -(*pz) / eln;
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__1 = bx;
/* Computing 2nd power */
    r__2 = by;
/* Computing 2nd power */
    r__3 = bz;
    scheck = 1.f - r__1 * r__1 - r__2 * r__2 - r__3 * r__3;
    if (scheck <= 0.f) {
	s_wsle(&io___1350);
	do_lio(&c__9, &c__1, "scheck28: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	s_stop("", (ftnlen)0);
    }
    ga = 1.f / sqrt(scheck);
/*       ga=1./sqrt(1.-bx**2-by**2-bz**2) */
    elnc = eln / ga;
/* Computing 2nd power */
    r__2 = elnc;
/* Computing 2nd power */
    r__3 = *dm2;
/* Computing 2nd power */
    r__4 = *amp;
/* Computing 2nd power */
    r__1 = (r__2 * r__2 + r__3 * r__3 - r__4 * r__4) / (elnc * 2.f);
/* Computing 2nd power */
    r__5 = *dm2;
    pn2 = r__1 * r__1 - r__5 * r__5;
    if (pn2 <= 0.f) {
	*icou1 = -1;
	return 0;
    }
    pn = sqrt(pn2);
/* lin-10/25/02 get rid of argument usage mismatch in PTR(): */
    xptr = pn * .33f;
/*       PNT=PTR(0.33*PN,ISEED) */
    pnt = ptr_(&xptr, iseed);
/* lin-10/25/02-end */
    fain = pi * 2.f * ranart_(&rndf77_1.nseed);
    *pnx = pnt * cos(fain);
    *pny = pnt * sin(fain);
    sig = 1.f;
    if (x > 0.f) {
	sig = -1.f;
    }
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__1 = pn;
/* Computing 2nd power */
    r__2 = pnt;
    scheck = r__1 * r__1 - r__2 * r__2;
    if (scheck < 0.f) {
	s_wsle(&io___1359);
	do_lio(&c__9, &c__1, "scheck29: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	scheck = 0.f;
    }
    *pnz = sig * sqrt(scheck);
/*       pnz=SIG*SQRT(pn**2-PNT**2) */
/* Computing 2nd power */
    r__1 = *dm2;
/* Computing 2nd power */
    r__2 = *pnx;
/* Computing 2nd power */
    r__3 = *pny;
/* Computing 2nd power */
    r__4 = *pnz;
    en = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
/* (2) the momentum for the rho */
    *ppx = -(*pnx);
    *ppy = -(*pny);
    *ppz = -(*pnz);
/* Computing 2nd power */
    r__1 = *amp;
/* Computing 2nd power */
    r__2 = *ppx;
/* Computing 2nd power */
    r__3 = *ppy;
/* Computing 2nd power */
    r__4 = *ppz;
    ep = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
/* (3) for the delta2, LORENTZ-TRANSFORMATION INTO nn cms FRAME */
    pbeta = *pnx * bx + *pny * by + *pnz * bz;
    trans0 = ga * (ga * pbeta / (ga + 1.f) + en);
    *pnx = bx * trans0 + *pnx;
    *pny = by * trans0 + *pny;
    *pnz = bz * trans0 + *pnz;
/* (4) for the rho, LORENTZ-TRANSFORMATION INTO nn cms FRAME */
    if (ep == 0.f) {
	ep = 1e-9f;
    }
    pbeta = *ppx * bx + *ppy * by + *ppz * bz;
    trans0 = ga * (ga * pbeta / (ga + 1.f) + ep);
    *ppx = bx * trans0 + *ppx;
    *ppy = by * trans0 + *ppy;
    *ppz = bz * trans0 + *ppz;
    return 0;
} /* pprho_ */

/* **************************8 */
/* *************************************** */
/* Subroutine */ int ppomga_(real *srt, integer *iseed, real *px, real *py, 
	real *pz, real *dm1, real *pnx, real *pny, real *pnz, real *dm2, real 
	*ppx, real *ppy, real *ppz, integer *icou1)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5, r__6;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double sqrt(doublereal), cos(doublereal), sin(doublereal);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static real v, w, x, f00, ga, ek, en, ep, x00, pi, bx, by, bz, pn, pt, 
	    pn2, fai, amn, amp, eln, sig, pnt;
    extern doublereal ptr_(real *, integer *);
    static real fain, elnc, xmax, xptr, fmax00, pbeta, ptmax, pzmax;
    static integer ntrym, ntryx;
    static real trans0, ptmax2, pzmax2, scheck;
    extern doublereal ranart_(integer *);
    static real xratio;

    /* Fortran I/O blocks */
    static cilist io___1372 = { 0, 99, 0, 0, 0 };
    static cilist io___1390 = { 0, 99, 0, 0, 0 };
    static cilist io___1399 = { 0, 99, 0, 0, 0 };


/* PURPOSE : CALCULATE MOMENTUM OF PARTICLES IN THE FINAL SATAT FROM */
/* THE PROCESS N+N--->N1+N2+OMEGA */
/*       DATE : Nov.5, 1994 */
/* Generate the masses and momentum for particles in the NN--> process */
/* for a given center of mass energy srt, the momenta are given in the center */
/* of mass of the NN */
/* **************************************** */
/* c      SAVE /TABLE/ */
/* c      SAVE /RNDF77/ */
    ntrym = 0;
    *icou1 = 0;
    pi = 3.1415926f;
    amn = .93892500000000001f;
    amp = .78200000000000003f;
    *dm1 = amn;
    *dm2 = amn;
/* CONSTANTS FOR GENERATING THE LONGITUDINAL MOMENTUM */
/* FOR ONE OF THE nucleons */
    v = .43f;
    w = -.84f;
/* (2) Generate the transverse momentum */
/*     OF p1 */
/* (2.1) estimate the maximum transverse momentum */
/* Computing 2nd power */
    r__1 = *srt;
/* Computing 2nd power */
    r__2 = *dm1 + *dm2 + amp;
/* Computing 2nd power */
    r__3 = *srt;
/* Computing 2nd power */
    r__4 = *dm1 - amp - *dm2;
/* Computing 2nd power */
    r__5 = *srt;
    ptmax2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / 4.f /
	     (r__5 * r__5);
/* lin-9/2012: check argument in sqrt(): */
    scheck = ptmax2;
    if (scheck < 0.f) {
	s_wsle(&io___1372);
	do_lio(&c__9, &c__1, "scheck30: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	scheck = 0.f;
    }
    ptmax = sqrt(scheck) * 1.f / 3.f;
/*       PTMAX=SQRT(PTMAX2)*1./3. */
L7:
    pt = ptr_(&ptmax, iseed);
/* (3) GENGRATE THE LONGITUDINAL MOMENTUM FOR DM1 */
/*     USING THE GIVEN DISTRIBUTION */
/* (3.1) THE MAXIMUM LONGITUDINAL MOMENTUM IS */
/* Computing 2nd power */
    r__1 = *srt;
/* Computing 2nd power */
    r__2 = *dm1 + *dm2 + amp;
/* Computing 2nd power */
    r__3 = *srt;
/* Computing 2nd power */
    r__4 = *dm1 - amp - *dm2;
/* Computing 2nd power */
    r__5 = *srt;
/* Computing 2nd power */
    r__6 = pt;
    pzmax2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / 4.f /
	     (r__5 * r__5) - r__6 * r__6;
    ++ntrym;
    if (pzmax2 < 0.f && ntrym <= 100) {
	goto L7;
    } else {
	pzmax2 = 1e-9f;
    }
    pzmax = sqrt(pzmax2);
    xmax = pzmax * 2.f / *srt;
/* (3.2) THE GENERATED X IS */
/* THE DSTRIBUTION HAS A MAXIMUM AT X0=-V/(2*w), f(X0)=1.056 */
    ntryx = 0;
    fmax00 = 1.056f;
    x00 = .26f;
    if (dabs(xmax) > .26f) {
	f00 = fmax00;
    } else {
/* Computing 2nd power */
	r__1 = xmax;
	f00 = v * dabs(xmax) + 1.f + w * (r__1 * r__1);
    }
L9:
    x = xmax * (1.f - ranart_(&rndf77_1.nseed) * 2.f);
    ++ntryx;
/* Computing 2nd power */
    r__1 = x;
    xratio = (v * dabs(x) + 1.f + w * (r__1 * r__1)) / f00;
/* lin-8/17/00       IF(xratio.LT.RANART(NSEED).and.ntryx.le.50)GO TO 9 */
    if (xratio < ranart_(&rndf77_1.nseed) && ntryx <= 50) {
	goto L9;
    }
/* (3.5) THE PZ IS */
    *pz = *srt * .5f * x;
/* The x and y components of the delta1 */
    fai = pi * 2.f * ranart_(&rndf77_1.nseed);
    *px = pt * cos(fai);
    *py = pt * sin(fai);
/* find the momentum of delta2 and rho */
/* the energy of the delta1 */
/* Computing 2nd power */
    r__1 = *dm1;
/* Computing 2nd power */
    r__2 = pt;
/* Computing 2nd power */
    r__3 = *pz;
    ek = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
/* (1) Generate the momentum of the delta2 in the cms of delta2 and rho */
/*     the energy of the cms of Drho */
    eln = *srt - ek;
    if (eln <= 0.f) {
	*icou1 = -1;
	return 0;
    }
    bx = -(*px) / eln;
    by = -(*py) / eln;
    bz = -(*pz) / eln;
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__1 = bx;
/* Computing 2nd power */
    r__2 = by;
/* Computing 2nd power */
    r__3 = bz;
    scheck = 1.f - r__1 * r__1 - r__2 * r__2 - r__3 * r__3;
    if (scheck <= 0.f) {
	s_wsle(&io___1390);
	do_lio(&c__9, &c__1, "scheck31: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	s_stop("", (ftnlen)0);
    }
    ga = 1.f / sqrt(scheck);
/*       ga=1./sqrt(1.-bx**2-by**2-bz**2) */
    elnc = eln / ga;
/* Computing 2nd power */
    r__2 = elnc;
/* Computing 2nd power */
    r__3 = *dm2;
/* Computing 2nd power */
    r__4 = amp;
/* Computing 2nd power */
    r__1 = (r__2 * r__2 + r__3 * r__3 - r__4 * r__4) / (elnc * 2.f);
/* Computing 2nd power */
    r__5 = *dm2;
    pn2 = r__1 * r__1 - r__5 * r__5;
    if (pn2 <= 0.f) {
	*icou1 = -1;
	return 0;
    }
    pn = sqrt(pn2);
/* lin-10/25/02 get rid of argument usage mismatch in PTR(): */
    xptr = pn * .33f;
/*       PNT=PTR(0.33*PN,ISEED) */
    pnt = ptr_(&xptr, iseed);
/* lin-10/25/02-end */
    fain = pi * 2.f * ranart_(&rndf77_1.nseed);
    *pnx = pnt * cos(fain);
    *pny = pnt * sin(fain);
    sig = 1.f;
    if (x > 0.f) {
	sig = -1.f;
    }
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__1 = pn;
/* Computing 2nd power */
    r__2 = pnt;
    scheck = r__1 * r__1 - r__2 * r__2;
    if (scheck < 0.f) {
	s_wsle(&io___1399);
	do_lio(&c__9, &c__1, "scheck32: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	scheck = 0.f;
    }
    *pnz = sig * sqrt(scheck);
/*       pnz=SIG*SQRT(pn**2-PNT**2) */
/* Computing 2nd power */
    r__1 = *dm2;
/* Computing 2nd power */
    r__2 = *pnx;
/* Computing 2nd power */
    r__3 = *pny;
/* Computing 2nd power */
    r__4 = *pnz;
    en = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
/* (2) the momentum for the rho */
    *ppx = -(*pnx);
    *ppy = -(*pny);
    *ppz = -(*pnz);
/* Computing 2nd power */
    r__1 = amp;
/* Computing 2nd power */
    r__2 = *ppx;
/* Computing 2nd power */
    r__3 = *ppy;
/* Computing 2nd power */
    r__4 = *ppz;
    ep = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
/* (3) for the delta2, LORENTZ-TRANSFORMATION INTO nn cms FRAME */
    pbeta = *pnx * bx + *pny * by + *pnz * bz;
    trans0 = ga * (ga * pbeta / (ga + 1.f) + en);
    *pnx = bx * trans0 + *pnx;
    *pny = by * trans0 + *pny;
    *pnz = bz * trans0 + *pnz;
/* (4) for the rho, LORENTZ-TRANSFORMATION INTO nn cms FRAME */
    if (ep == 0.f) {
	ep = 1e-9f;
    }
    pbeta = *ppx * bx + *ppy * by + *ppz * bz;
    trans0 = ga * (ga * pbeta / (ga + 1.f) + ep);
    *ppx = bx * trans0 + *ppx;
    *ppy = by * trans0 + *ppy;
    *ppz = bz * trans0 + *ppz;
    return 0;
} /* ppomga_ */

/* **************************8 */
/* **************************8 */
/*   DELTA MASS GENERATOR */
doublereal rmass_(real *dmax__, integer *iseed)
{
    /* System generated locals */
    real ret_val;

    /* Local variables */
    static real dm, fm, dmin__;
    static integer ntry1;
    extern doublereal fdelta_(real *), ranart_(integer *);

/* c      SAVE /RNDF77/ */
/* THE MINIMUM MASS FOR DELTA */
    dmin__ = 1.078f;
/* Delta(1232) production */
    if (*dmax__ < 1.232f) {
	fm = fdelta_(dmax__);
    } else {
	fm = 1.f;
    }
    if (fm == 0.f) {
	fm = 1e-6f;
    }
    ntry1 = 0;
L10:
    dm = ranart_(&rndf77_1.nseed) * (*dmax__ - dmin__) + dmin__;
    ++ntry1;
    if (ranart_(&rndf77_1.nseed) > fdelta_(&dm) / fm && ntry1 <= 10) {
	goto L10;
    }
/* lin-2/26/03 sometimes Delta mass can reach very high values (e.g. 15.GeV), */
/*     thus violating the thresh of the collision which produces it */
/*     and leads to large violation of energy conservation. */
/*     To limit the above, limit the Delta mass below a certain value */
/*     (here taken as its central value + 2* B-W fullwidth): */
    if (dm > 1.47f) {
	goto L10;
    }
    ret_val = dm;
    return ret_val;
} /* rmass_ */

/* ------------------------------------------------------------------ */
/* THE Breit Wigner FORMULA */
doublereal frho_(real *dmass)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3;

    /* Local variables */
    static real fd, am0, wid;

    am0 = .77f;
    wid = .153f;
/* Computing 2nd power */
    r__1 = wid;
/* Computing 2nd power */
    r__2 = *dmass - am0;
/* Computing 2nd power */
    r__3 = wid;
    fd = r__1 * r__1 * .25f / (r__2 * r__2 + r__3 * r__3 * .25f);
    ret_val = fd;
    return ret_val;
} /* frho_ */

/* **************************8 */
/*   RHO MASS GENERATOR */
doublereal rhomas_(real *dmax__, integer *iseed)
{
    /* System generated locals */
    real ret_val;

    /* Local variables */
    static real dm, fm, dmin__;
    extern doublereal frho_(real *);
    static integer ntry1;
    extern doublereal ranart_(integer *);

/* c      SAVE /RNDF77/ */
/* THE MINIMUM MASS FOR DELTA */
    dmin__ = .28f;
/* RHO(770) production */
    if (*dmax__ < .77f) {
	fm = frho_(dmax__);
    } else {
	fm = 1.f;
    }
    if (fm == 0.f) {
	fm = 1e-6f;
    }
    ntry1 = 0;
L10:
    dm = ranart_(&rndf77_1.nseed) * (*dmax__ - dmin__) + dmin__;
    ++ntry1;
    if (ranart_(&rndf77_1.nseed) > frho_(&dm) / fm && ntry1 <= 10) {
	goto L10;
    }
/* lin-2/26/03 limit the rho mass below a certain value */
/*     (here taken as its central value + 2* B-W fullwidth): */
    if (dm > 1.07f) {
	goto L10;
    }
    ret_val = dm;
    return ret_val;
} /* rhomas_ */

/* ***************************************** */
/* for pp-->pp+2pi */
/*      real*4 function X2pi(srt) */
doublereal x2pi_(real *srt)
{
    /* Initialized data */

    static real earray[15] = { 2.23f,2.81f,3.67f,4.f,4.95f,5.52f,5.97f,6.04f,
	    6.6f,6.9f,7.87f,8.11f,10.01f,16.f,19.f };
    static real xarray[15] = { 1.22f,2.51f,2.67f,2.95f,2.96f,2.84f,2.8f,3.2f,
	    2.7f,3.f,2.54f,2.46f,2.4f,1.66f,1.5f };

    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal), exp(doublereal);

    /* Local variables */
    static integer ie;
    static real plab, xmin, ymin, xmax, ymax, pmass;

/*  This function contains the experimental */
/*     total pp-pp+pi(+)pi(-) Xsections    * */
/*  srt    = DSQRT(s) in GeV                                                  * */
/*  xsec   = production cross section in mb                                   * */
/*  earray = EXPerimental table with proton momentum in GeV/c                 * */
/*  xarray = EXPerimental table with cross sections in mb (curve to guide eye)* */
/*                                                                            * */
/* ***************************************** */
/*      real*4   xarray(15), earray(15) */
    pmass = .9383f;
/* 1.Calculate p(lab)  from srt [GeV] */
/*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1) */
/*      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.) */
    ret_val = 1e-6f;
    if (*srt <= 2.2f) {
	return ret_val;
    }
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = pmass;
/* Computing 2nd power */
    r__1 = (r__2 * r__2 - r__3 * r__3 * 2.f) / (pmass * 2.f);
/* Computing 2nd power */
    r__4 = pmass;
    plab = sqrt(r__1 * r__1 - r__4 * r__4);
    if (plab < earray[0]) {
	ret_val = xarray[0];
	return ret_val;
    }

/* 2.Interpolate double logarithmically to find sigma(srt) */

    for (ie = 1; ie <= 15; ++ie) {
	if (earray[ie - 1] == plab) {
	    ret_val = xarray[ie - 1];
	    return ret_val;
	} else if (earray[ie - 1] > plab) {
	    ymin = log(xarray[ie - 2]);
	    ymax = log(xarray[ie - 1]);
	    xmin = log(earray[ie - 2]);
	    xmax = log(earray[ie - 1]);
	    ret_val = exp(ymin + (log(plab) - xmin) * (ymax - ymin) / (xmax - 
		    xmin));
	    return ret_val;
	}
/* L1001: */
    }
    return ret_val;
} /* x2pi_ */

/* ***************************************** */
/* for pp-->pn+pi(+)pi(+)pi(-) */
/*      real*4 function X3pi(srt) */
doublereal x3pi_(real *srt)
{
    /* Initialized data */

    static real xarray[12] = { .02f,.4f,1.15f,1.6f,2.19f,2.85f,2.3f,3.1f,
	    2.47f,2.6f,2.4f,1.7f };
    static real earray[12] = { 2.23f,2.81f,3.67f,4.f,4.95f,5.52f,5.97f,6.04f,
	    6.6f,6.9f,10.01f,19.f };

    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal), exp(doublereal);

    /* Local variables */
    static integer ie;
    static real plab, xmin, ymin, xmax, ymax, pmass;

/*  This function contains the experimental pp->pp+3pi cross sections          * */
/*  srt    = DSQRT(s) in GeV                                                   * */
/*  xsec   = production cross section in mb                                    * */
/*  earray = EXPerimental table with proton energies in MeV                    * */
/*  xarray = EXPerimental table with cross sections in mb (curve to guide eye) * */
/*                                                                             * */
/* ***************************************** */
/*      real*4   xarray(12), earray(12) */
    pmass = .9383f;
/* 1.Calculate p(lab)  from srt [GeV] */
/*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1) */
/*      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.) */
    ret_val = 1e-6f;
    if (*srt <= 2.3f) {
	return ret_val;
    }
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = pmass;
/* Computing 2nd power */
    r__1 = (r__2 * r__2 - r__3 * r__3 * 2.f) / (pmass * 2.f);
/* Computing 2nd power */
    r__4 = pmass;
    plab = sqrt(r__1 * r__1 - r__4 * r__4);
    if (plab < earray[0]) {
	ret_val = xarray[0];
	return ret_val;
    }

/* 2.Interpolate double logarithmically to find sigma(srt) */

    for (ie = 1; ie <= 12; ++ie) {
	if (earray[ie - 1] == plab) {
	    ret_val = xarray[ie - 1];
	    return ret_val;
	} else if (earray[ie - 1] > plab) {
	    ymin = log(xarray[ie - 2]);
	    ymax = log(xarray[ie - 1]);
	    xmin = log(earray[ie - 2]);
	    xmax = log(earray[ie - 1]);
	    ret_val = exp(ymin + (log(plab) - xmin) * (ymax - ymin) / (xmax - 
		    xmin));
	    return ret_val;
	}
/* L1001: */
    }
    return ret_val;
} /* x3pi_ */

/* ***************************************** */
/* ***************************************** */
/* for pp-->pp+pi(+)pi(-)pi(0) */
/*      real*4 function X33pi(srt) */
doublereal x33pi_(real *srt)
{
    /* Initialized data */

    static real xarray[12] = { .02f,.22f,.74f,1.1f,1.76f,1.84f,2.2f,2.4f,
	    2.15f,2.6f,2.3f,1.7f };
    static real earray[12] = { 2.23f,2.81f,3.67f,4.f,4.95f,5.52f,5.97f,6.04f,
	    6.6f,6.9f,10.01f,19.f };

    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal), exp(doublereal);

    /* Local variables */
    static integer ie;
    static real plab, xmin, ymin, xmax, ymax, pmass;

/*  This function contains the experimental pp->pp+3pi cross sections          * */
/*  srt    = DSQRT(s) in GeV                                                   * */
/*  xsec   = production cross section in mb                                    * */
/*  earray = EXPerimental table with proton energies in MeV                    * */
/*  xarray = EXPerimental table with cross sections in mb (curve to guide eye) * */
/*                                                                             * */
/* ***************************************** */
/*      real*4   xarray(12), earray(12) */
    pmass = .9383f;
    ret_val = 1e-6f;
    if (*srt <= 2.3f) {
	return ret_val;
    }
/* 1.Calculate p(lab)  from srt [GeV] */
/*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1) */
/*      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.) */
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = pmass;
/* Computing 2nd power */
    r__1 = (r__2 * r__2 - r__3 * r__3 * 2.f) / (pmass * 2.f);
/* Computing 2nd power */
    r__4 = pmass;
    plab = sqrt(r__1 * r__1 - r__4 * r__4);
    if (plab < earray[0]) {
	ret_val = xarray[0];
	return ret_val;
    }

/* 2.Interpolate double logarithmically to find sigma(srt) */

    for (ie = 1; ie <= 12; ++ie) {
	if (earray[ie - 1] == plab) {
	    ret_val = xarray[ie - 1];
	    return ret_val;
	} else if (earray[ie - 1] > plab) {
	    ymin = log(xarray[ie - 2]);
	    ymax = log(xarray[ie - 1]);
	    xmin = log(earray[ie - 2]);
	    xmax = log(earray[ie - 1]);
	    ret_val = exp(ymin + (log(plab) - xmin) * (ymax - ymin) / (xmax - 
		    xmin));
	    return ret_val;
	}
/* L1001: */
    }
    return ret_val;
} /* x33pi_ */

/* ***************************************** */
/*       REAL*4 FUNCTION X4pi(SRT) */
doublereal x4pi_(real *srt)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real al, as, es, ak0;
    extern doublereal pp1_(real *);
    static real xk1, xk2, xk3, xk4, ada, ana, akp;
    extern doublereal s1535_(real *), ppk0_(real *), ppk1_(real *), x2pi_(
	    real *), x3pi_(real *);
    static real pps1, pps2;
    extern doublereal x33pi_(real *);
    static real t1dlk, t2dlk, t1dsk, t2dsk, t1nlk, t2nlk, t1nsk, t2nsk;
    extern doublereal sigma_(real *, integer *, integer *, integer *);
    static real pmdlk, pmdsk, xkaon, pmnlk, pmass;
    extern doublereal pplpk_(real *);
    static real pmnsk, pmdlk2, pmdsk2, pmnlk2, pmnsk2, xpp2pi, xpp3pi, ppsngl;

/*       CROSS SECTION FOR NN-->DD+rho PROCESS */
/* ***************************** */
    akp = .498f;
    ak0 = .498f;
    ana = .94f;
    ada = 1.232f;
    al = 1.1157f;
    as = 1.1197f;
    pmass = .9383f;
    es = *srt;
    if (es <= 4.f) {
	ret_val = 0.f;
    } else {
/* cross section for two resonance pp-->DD+DN*+N*N* */
	xpp2pi = x2pi_(&es) * 4.f;
/* cross section for pp-->pp+spi */
	xpp3pi = (x3pi_(&es) + x33pi_(&es)) * 3.f;
/* cross section for pp-->pD+ and nD++ */
	pps1 = sigma_(&es, &c__1, &c__1, &c__0) + sigma_(&es, &c__1, &c__1, &
		c__1) * .5f;
	pps2 = sigma_(&es, &c__1, &c__1, &c__1) * 1.5f;
	ppsngl = pps1 + pps2 + s1535_(&es);
/* CROSS SECTION FOR KAON PRODUCTION from the four channels */
/* for NLK channel */
	xk1 = 0.f;
	xk2 = 0.f;
	xk3 = 0.f;
	xk4 = 0.f;
	t1nlk = ana + al + akp;
	t2nlk = ana + al - akp;
	if (es <= t1nlk) {
	    goto L333;
	}
/* Computing 2nd power */
	r__1 = es;
/* Computing 2nd power */
	r__2 = t1nlk;
/* Computing 2nd power */
	r__3 = es;
/* Computing 2nd power */
	r__4 = t2nlk;
/* Computing 2nd power */
	r__5 = es;
	pmnlk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (
		r__5 * r__5 * 4.f);
	pmnlk = sqrt(pmnlk2);
	xk1 = pplpk_(&es);
/* for DLK channel */
	t1dlk = ada + al + akp;
	t2dlk = ada + al - akp;
	if (es <= t1dlk) {
	    goto L333;
	}
/* Computing 2nd power */
	r__1 = es;
/* Computing 2nd power */
	r__2 = t1dlk;
/* Computing 2nd power */
	r__3 = es;
/* Computing 2nd power */
	r__4 = t2dlk;
/* Computing 2nd power */
	r__5 = es;
	pmdlk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (
		r__5 * r__5 * 4.f);
	pmdlk = sqrt(pmdlk2);
	xk3 = pplpk_(&es);
/* for NSK channel */
	t1nsk = ana + as + akp;
	t2nsk = ana + as - akp;
	if (es <= t1nsk) {
	    goto L333;
	}
/* Computing 2nd power */
	r__1 = es;
/* Computing 2nd power */
	r__2 = t1nsk;
/* Computing 2nd power */
	r__3 = es;
/* Computing 2nd power */
	r__4 = t2nsk;
/* Computing 2nd power */
	r__5 = es;
	pmnsk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (
		r__5 * r__5 * 4.f);
	pmnsk = sqrt(pmnsk2);
	xk2 = ppk1_(&es) + ppk0_(&es);
/* for DSK channel */
	t1dsk = ada + as + akp;
	t2dsk = ada + as - akp;
	if (es <= t1dsk) {
	    goto L333;
	}
/* Computing 2nd power */
	r__1 = es;
/* Computing 2nd power */
	r__2 = t1dsk;
/* Computing 2nd power */
	r__3 = es;
/* Computing 2nd power */
	r__4 = t2dsk;
/* Computing 2nd power */
	r__5 = es;
	pmdsk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (
		r__5 * r__5 * 4.f);
	pmdsk = sqrt(pmdsk2);
	xk4 = ppk1_(&es) + ppk0_(&es);
/* THE TOTAL KAON+ AND KAON0 PRODUCTION CROSS SECTION IS THEN */
L333:
	xkaon = (xk1 + xk2 + xk3 + xk4) * 3.f;
/* cross section for pp-->DD+rho */
	ret_val = pp1_(&es) - ppsngl - xpp2pi - xpp3pi - xkaon;
	if (ret_val <= 0.f) {
	    ret_val = 1e-6f;
	}
    }
    return ret_val;
} /* x4pi_ */

/* ***************************************** */
/* for pp-->inelastic */
/*      real*4 function pp1(srt) */
doublereal pp1_(real *srt)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal), pow_dd(doublereal *, doublereal 
	    *);

    /* Local variables */
    static real a, b, c__, d__, an, plab, pmin, pmax, plab2, pmass;

/*  srt    = DSQRT(s) in GeV                                                   * */
/*  xsec   = production cross section in mb                                    * */
/*  earray = EXPerimental table with proton energies in MeV                    * */
/*  xarray = EXPerimental table with cross sections in mb (curve to guide eye) * */
/*                                                                             * */
/* ***************************************** */
    pmass = .9383f;
    ret_val = 0.f;
/* 1.Calculate p(lab)  from srt [GeV] */
/*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1) */
/*      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.) */
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = pmass;
/* Computing 2nd power */
    r__1 = (r__2 * r__2 - r__3 * r__3 * 2.f) / (pmass * 2.f);
/* Computing 2nd power */
    r__4 = pmass;
    plab2 = r__1 * r__1 - r__4 * r__4;
    if (plab2 <= 0.f) {
	return ret_val;
    }
    plab = sqrt(plab2);
    pmin = .968f;
    pmax = 2080.f;
    if (plab < pmin || plab > pmax) {
	ret_val = 0.f;
	return ret_val;
    }
/* * fit parameters */
    a = 30.9f;
    b = -28.9f;
    c__ = .192f;
    d__ = -.835f;
    an = -2.46f;
    d__1 = (doublereal) plab;
    d__2 = (doublereal) an;
/* Computing 2nd power */
    r__1 = log(plab);
    ret_val = a + b * pow_dd(&d__1, &d__2) + c__ * (r__1 * r__1);
    if (ret_val <= 0.f) {
	ret_val = 0.f;
    }
    return ret_val;
} /* pp1_ */

/* ***************************************** */
/* for pp-->elastic */
/*      real*4 function pp2(srt) */
doublereal pp2_(real *srt)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double sqrt(doublereal), log(doublereal), pow_dd(doublereal *, doublereal 
	    *);

    /* Local variables */
    static real a, b, c__, d__, an, plab, pmin, pmax, pmass, scheck;

    /* Fortran I/O blocks */
    static cilist io___1488 = { 0, 99, 0, 0, 0 };


/*  srt    = DSQRT(s) in GeV                                                   * */
/*  xsec   = production cross section in mb                                    * */
/*  earray = EXPerimental table with proton energies in MeV                    * */
/*  xarray = EXPerimental table with cross sections in mb (curve to guide eye) * */
/*                                                                             * */
/* ***************************************** */
    pmass = .9383f;
/* 1.Calculate p(lab)  from srt [GeV] */
/*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1) */
/*      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.) */
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = pmass;
/* Computing 2nd power */
    r__1 = (r__2 * r__2 - r__3 * r__3 * 2.f) / (pmass * 2.f);
/* Computing 2nd power */
    r__4 = pmass;
    scheck = r__1 * r__1 - r__4 * r__4;
    if (scheck < 0.f) {
	s_wsle(&io___1488);
	do_lio(&c__9, &c__1, "scheck33: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	scheck = 0.f;
    }
    plab = sqrt(scheck);
/*      plab=sqrt(((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2) */
    pmin = 2.f;
    pmax = 2050.f;
    if (plab > pmax) {
	ret_val = 8.f;
	return ret_val;
    }
    if (plab < pmin) {
	ret_val = 25.f;
	return ret_val;
    }
/* * fit parameters */
    a = 11.2f;
    b = 25.5f;
    c__ = .151f;
    d__ = -1.62f;
    an = -1.12f;
    d__1 = (doublereal) plab;
    d__2 = (doublereal) an;
/* Computing 2nd power */
    r__1 = log(plab);
    ret_val = a + b * pow_dd(&d__1, &d__2) + c__ * (r__1 * r__1) + d__ * log(
	    plab);
    if (ret_val <= 0.f) {
	ret_val = 0.f;
    }
    return ret_val;
} /* pp2_ */

/* ***************************************** */
/* for pp-->total */
/*      real*4 function ppt(srt) */
doublereal ppt_(real *srt)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double sqrt(doublereal), log(doublereal), pow_dd(doublereal *, doublereal 
	    *);

    /* Local variables */
    static real a, b, c__, d__, an, plab, pmin, pmax, pmass, scheck;

    /* Fortran I/O blocks */
    static cilist io___1499 = { 0, 99, 0, 0, 0 };


/*  srt    = DSQRT(s) in GeV                                                   * */
/*  xsec   = production cross section in mb                                    * */
/*  earray = EXPerimental table with proton energies in MeV                    * */
/*  xarray = EXPerimental table with cross sections in mb (curve to guide eye) * */
/*                                                                             * */
/* ***************************************** */
    pmass = .9383f;
/* 1.Calculate p(lab)  from srt [GeV] */
/*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1) */
/*      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.) */
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = pmass;
/* Computing 2nd power */
    r__1 = (r__2 * r__2 - r__3 * r__3 * 2.f) / (pmass * 2.f);
/* Computing 2nd power */
    r__4 = pmass;
    scheck = r__1 * r__1 - r__4 * r__4;
    if (scheck < 0.f) {
	s_wsle(&io___1499);
	do_lio(&c__9, &c__1, "scheck34: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	scheck = 0.f;
    }
    plab = sqrt(scheck);
/*      plab=sqrt(((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2) */
    pmin = 3.f;
    pmax = 2100.f;
    if (plab < pmin || plab > pmax) {
	ret_val = 55.f;
	return ret_val;
    }
/* * fit parameters */
    a = 45.6f;
    b = 219.f;
    c__ = .41f;
    d__ = -3.41f;
    an = -4.23f;
    d__1 = (doublereal) plab;
    d__2 = (doublereal) an;
/* Computing 2nd power */
    r__1 = log(plab);
    ret_val = a + b * pow_dd(&d__1, &d__2) + c__ * (r__1 * r__1) + d__ * log(
	    plab);
    if (ret_val <= 0.f) {
	ret_val = 0.f;
    }
    return ret_val;
} /* ppt_ */

/* ************************ */
/* cross section for N*(1535) production in PP collisions */
/* VARIABLES: */
/* LB1,LB2 ARE THE LABLES OF THE TWO COLLIDING PARTICLES */
/* SRT IS THE CMS ENERGY */
/* X1535 IS THE N*(1535) PRODUCTION CROSS SECTION */
/* NOTE THAT THE N*(1535) PRODUCTION CROSS SECTION IS 2 TIMES THE ETA */
/* PRODUCTION CROSS SECTION */
/* DATE: Aug. 1 , 1994 */
/* ******************************** */
doublereal s1535_(real *srt)
{
    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real s0;

    s0 = 2.424f;
    ret_val = 0.f;
    if (*srt <= s0) {
	return ret_val;
    }
/* Computing 2nd power */
    r__1 = *srt - s0;
    ret_val = (*srt - s0) * .20399999999999999f / (r__1 * r__1 + .058f);
    return ret_val;
} /* s1535_ */

/* *************************************** */
/* generate a table for pt distribution for */
/* Subroutine */ int tablem_(void)
{
    static integer l;
    static real x, rr, anorm;
    extern doublereal ptdis_(real *);
    static real ptmax;

/* THE PROCESS N+N--->N+N+PION */
/*       DATE : July 11, 1994 */
/* **************************************** */
/* c      SAVE /TABLE/ */
    ptmax = 2.01f;
    anorm = ptdis_(&ptmax);
    for (l = 0; l <= 200; ++l) {
	x = (real) (l + 1) * .01f;
	rr = ptdis_(&x) / anorm;
	table_1.earray[l] = rr;
	table_1.xarray[l] = x;
/* L10: */
    }
    return 0;
} /* tablem_ */

/* ******************************** */
doublereal ptdis_(real *x)
{
    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static real b, c__, d__;

/* NUCLEON TRANSVERSE MOMENTUM DISTRIBUTION AT HIGH ENERGIES */
/* DATE: Aug. 11, 1994 */
/* ******************************** */
    b = 3.78f;
    c__ = .47f;
    d__ = 3.6f;
/*       b=b*3 */
/*       d=d*3 */
/* Computing 2nd power */
    r__1 = *x;
/* Computing 2nd power */
    r__2 = d__;
    ret_val = 1.f / (b * 2.f) * (1.f - exp(-b * (r__1 * r__1))) - c__ / d__ * 
	    *x * exp(-d__ * *x) - c__ / (r__2 * r__2) * (exp(-d__ * *x) - 1.f)
	    ;
    return ret_val;
} /* ptdis_ */

/* **************************** */
/* Subroutine */ int ppxs_(integer *lb1, integer *lb2, real *srt, real *ppsig,
	 real *spprho, integer *ipp)
{
    /* System generated locals */
    real r__1, r__2;

    /* Builtin functions */
    double sqrt(doublereal), atan(doublereal), sin(doublereal);

    /* Local variables */
    static real q, s0, s1, s2, d00, d11, d20, erh, esi, erho, trho, esigma, 
	    tsigma;

/* purpose: this subroutine gives the cross section for pion+pion */
/*          elastic collision */
/* variables: */
/*       input: lb1,lb2 and srt are the labels and srt for I1 and I2 */
/*       output: ppsig: pp xsection */
/*               ipp: label for the pion+pion channel */
/*               Ipp=0 NOTHING HAPPEND */
/*                  1 for Pi(+)+PI(+) DIRECT */
/*                   2     PI(+)+PI(0) FORMING RHO(+) */
/*                  3     PI(+)+PI(-) FORMING RHO(0) */
/*                   4     PI(0)+PI(O) DIRECT */
/*                  5     PI(0)+PI(-) FORMING RHO(-) */
/*                  6     PI(-)+PI(-) DIRECT */
/* reference: G.F. Bertsch, Phys. Rev. D37 (1988) 1202. */
/* date : Aug 29, 1994 */
/* **************************** */
    *ppsig = 0.f;
/* bzdbg10/15/99 */
    *spprho = 0.f;
/* bzdbg10/15/99 end */
    *ipp = 0;
    if (*srt <= .3f) {
	return 0;
    }
/* Computing 2nd power */
    r__1 = *srt / 2;
    q = sqrt(r__1 * r__1 - .019600000000000003f);
    esigma = .81200000000000006f;
    tsigma = q * 2.06f;
    erho = .77f;
/* Computing 2nd power */
    r__2 = q / erho;
/* Computing 2nd power */
    r__1 = q / .14f / (r__2 * r__2 + 1.f);
    trho = q * .095f * (r__1 * r__1);
    esi = esigma - *srt;
    if (esi == 0.f) {
	d00 = 1.5707963f;
	goto L10;
    }
    d00 = atan(tsigma / 2.f / esi);
L10:
    erh = erho - *srt;
    if (erh == 0.f) {
	d11 = 1.5707963f;
	goto L20;
    }
    d11 = atan(trho / 2.f / erh);
L20:
    d20 = q * -.12f / .14f;
/* Computing 2nd power */
    r__1 = sin(d00);
/* Computing 2nd power */
    r__2 = q;
    s0 = r__1 * r__1 * 25.132740800000001f / (r__2 * r__2);
/* Computing 2nd power */
    r__1 = sin(d11);
/* Computing 2nd power */
    r__2 = q;
    s1 = r__1 * r__1 * 75.398222400000009f / (r__2 * r__2);
/* Computing 2nd power */
    r__1 = sin(d20);
/* Computing 2nd power */
    r__2 = q;
    s2 = r__1 * r__1 * 125.663704f / (r__2 * r__2);
/*    !! GeV^-2 to mb */
    s0 = s0 * .038809000000000003f * 10.f;
    s1 = s1 * .038809000000000003f * 10.f;
    s2 = s2 * .038809000000000003f * 10.f;
/*       ppXS=s0/9.+s1/3.+s2*0.56 */
/*       if(ppxs.le.0)ppxs=0.00001 */
    *spprho = s1 / 2.f;
/* (1) PI(+)+PI(+) */
    if (*lb1 == 5 && *lb2 == 5) {
	*ipp = 1;
	*ppsig = s2;
	return 0;
    }
/* (2) PI(+)+PI(0) */
    if (*lb1 == 5 && *lb2 == 4 || *lb1 == 4 && *lb2 == 5) {
	*ipp = 2;
	*ppsig = s2 / 2.f + s1 / 2.f;
	return 0;
    }
/* (3) PI(+)+PI(-) */
    if (*lb1 == 5 && *lb2 == 3 || *lb1 == 3 && *lb2 == 5) {
	*ipp = 3;
	*ppsig = s2 / 6.f + s1 / 2.f + s0 / 3.f;
	return 0;
    }
/* (4) PI(0)+PI(0) */
    if (*lb1 == 4 && *lb2 == 4) {
	*ipp = 4;
	*ppsig = s2 * 2 / 3.f + s0 / 3.f;
	return 0;
    }
/* (5) PI(0)+PI(-) */
    if (*lb1 == 4 && *lb2 == 3 || *lb1 == 3 && *lb2 == 4) {
	*ipp = 5;
	*ppsig = s2 / 2.f + s1 / 2.f;
	return 0;
    }
/* (6) PI(-)+PI(-) */
    if (*lb1 == 3 && *lb2 == 3) {
	*ipp = 6;
	*ppsig = s2;
    }
    return 0;
} /* ppxs_ */

/* ********************************* */
/* elementary kaon production cross sections */
/*  from the CERN data book */
/*  date: Sept.2, 1994 */
/*  for pp-->pLK+ */
/*      real*4 function pplpk(srt) */
doublereal pplpk_(real *srt)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double sqrt(doublereal), log(doublereal), pow_dd(doublereal *, doublereal 
	    *);

    /* Local variables */
    static real a, b, c__, an, plab, pmin, pmax, pmass, scheck;

    /* Fortran I/O blocks */
    static cilist io___1532 = { 0, 99, 0, 0, 0 };


/*  srt    = DSQRT(s) in GeV                                                   * */
/*  xsec   = production cross section in mb                                    * */
/*  earray = EXPerimental table with proton energies in MeV                    * */
/*  xarray = EXPerimental table with cross sections in mb (curve to guide eye) * */
/*                                                                             * */
/* ***************************************** */
    pmass = .9383f;
/* 1.Calculate p(lab)  from srt [GeV] */
/*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1) */
/*   find the center of mass energy corresponding to the given pm as */
/*   if Lambda+N+K are produced */
    ret_val = 0.f;
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = pmass;
/* Computing 2nd power */
    r__1 = (r__2 * r__2 - r__3 * r__3 * 2.f) / (pmass * 2.f);
/* Computing 2nd power */
    r__4 = pmass;
    scheck = r__1 * r__1 - r__4 * r__4;
    if (scheck < 0.f) {
	s_wsle(&io___1532);
	do_lio(&c__9, &c__1, "scheck35: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	scheck = 0.f;
    }
    plab = sqrt(scheck);
/*        plab=sqrt(((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2) */
    pmin = 2.82f;
    pmax = 25.f;
    if (plab > pmax) {
	ret_val = .036f;
	return ret_val;
    }
    if (plab < pmin) {
	ret_val = 0.f;
	return ret_val;
    }
/* * fit parameters */
    a = .0654f;
    b = -3.16f;
    c__ = -.0029f;
    an = -4.14f;
    d__1 = (doublereal) plab;
    d__2 = (doublereal) an;
/* Computing 2nd power */
    r__1 = log(plab);
    ret_val = a + b * pow_dd(&d__1, &d__2) + c__ * (r__1 * r__1);
    if (ret_val <= 0.f) {
	ret_val = 0.f;
    }
    return ret_val;
} /* pplpk_ */

/* ***************************************** */
/* for pp-->pSigma+K0 */
/*      real*4 function ppk0(srt) */
doublereal ppk0_(real *srt)
{
    /* Initialized data */

    static real xarray[7] = { .03f,.025f,.025f,.026f,.02f,.014f,.06f };
    static real earray[7] = { 3.67f,4.95f,5.52f,6.05f,6.92f,7.87f,10.f };

    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal), exp(doublereal);

    /* Local variables */
    static integer ie;
    static real plab, xmin, ymin, xmax, ymax, pmass;

/*  srt    = DSQRT(s) in GeV                                                   * */
/*  xsec   = production cross section in mb                                    * */
/*                                                                             * */
/* ***************************************** */
/*      real*4   xarray(7), earray(7) */
    pmass = .9383f;
/* 1.Calculate p(lab)  from srt [GeV] */
/*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1) */
/*      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.) */
    ret_val = 0.f;
    if (*srt <= 2.63f) {
	return ret_val;
    }
    if (*srt > 4.54f) {
	ret_val = .037f;
	return ret_val;
    }
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = pmass;
/* Computing 2nd power */
    r__1 = (r__2 * r__2 - r__3 * r__3 * 2.f) / (pmass * 2.f);
/* Computing 2nd power */
    r__4 = pmass;
    plab = sqrt(r__1 * r__1 - r__4 * r__4);
    if (plab < earray[0]) {
	ret_val = xarray[0];
	return ret_val;
    }

/* 2.Interpolate double logarithmically to find sigma(srt) */

    for (ie = 1; ie <= 7; ++ie) {
	if (earray[ie - 1] == plab) {
	    ret_val = xarray[ie - 1];
	    goto L10;
	} else if (earray[ie - 1] > plab) {
	    ymin = log(xarray[ie - 2]);
	    ymax = log(xarray[ie - 1]);
	    xmin = log(earray[ie - 2]);
	    xmax = log(earray[ie - 1]);
	    ret_val = exp(ymin + (log(plab) - xmin) * (ymax - ymin) / (xmax - 
		    xmin));
	    goto L10;
	}
/* L1001: */
    }
L10:
    return ret_val;
} /* ppk0_ */

/* ***************************************** */
/* for pp-->pSigma0K+ */
/*      real*4 function ppk1(srt) */
doublereal ppk1_(real *srt)
{
    /* Initialized data */

    static real xarray[7] = { .013f,.025f,.016f,.012f,.017f,.029f,.025f };
    static real earray[7] = { 3.67f,4.95f,5.52f,5.97f,6.05f,6.92f,7.87f };

    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal), exp(doublereal);

    /* Local variables */
    static integer ie;
    static real plab, xmin, ymin, xmax, ymax, pmass;

/*  srt    = DSQRT(s) in GeV                                                   * */
/*  xsec   = production cross section in mb                                    * */
/*                                                                             * */
/* ***************************************** */
/*      real*4   xarray(7), earray(7) */
    pmass = .9383f;
/* 1.Calculate p(lab)  from srt [GeV] */
/*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1) */
/*      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.) */
    ret_val = 0.f;
    if (*srt <= 2.63f) {
	return ret_val;
    }
    if (*srt > 4.08f) {
	ret_val = .025f;
	return ret_val;
    }
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = pmass;
/* Computing 2nd power */
    r__1 = (r__2 * r__2 - r__3 * r__3 * 2.f) / (pmass * 2.f);
/* Computing 2nd power */
    r__4 = pmass;
    plab = sqrt(r__1 * r__1 - r__4 * r__4);
    if (plab < earray[0]) {
	ret_val = xarray[0];
	return ret_val;
    }

/* 2.Interpolate double logarithmically to find sigma(srt) */

    for (ie = 1; ie <= 7; ++ie) {
	if (earray[ie - 1] == plab) {
	    ret_val = xarray[ie - 1];
	    goto L10;
	} else if (earray[ie - 1] > plab) {
	    ymin = log(xarray[ie - 2]);
	    ymax = log(xarray[ie - 1]);
	    xmin = log(earray[ie - 2]);
	    xmax = log(earray[ie - 1]);
	    ret_val = exp(ymin + (log(plab) - xmin) * (ymax - ymin) / (xmax - 
		    xmin));
	    goto L10;
	}
/* L1001: */
    }
L10:
    return ret_val;
} /* ppk1_ */

/* ********************************* */
/*                                                                      * */
/*                                                                      * */
/* Subroutine */ int crpn_(real *px, real *py, real *pz, real *srt, integer *
	i1, integer *i2, integer *iblock, real *xkaon0, real *xkaon, real *
	xphi, real *xphin)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static real c1, s1, t1, x1, x2, dm;
    static integer ii, jj;
    static real pr, cc1, em1, em2, ct1, pr2, px0, py0, pz0, st1;
    static integer ipi;
    extern doublereal ptr_(real *, integer *);
    static integer ntag;
    static real dmax__, arho, xptr;
    static integer kaonc, ianti;
    extern doublereal pnlka_(real *), pnska_(real *), rmass_(real *, integer *
	    ), twopi_(real *);
    static real aomega, scheck;
    extern /* Subroutine */ int rmasdd_(real *, real *, real *, real *, real *
	    , integer *, integer *, real *, real *);
    static real ameson;
    extern doublereal ranart_(integer *), threpi_(real *);
    extern /* Subroutine */ int rotate_(real *, real *, real *, real *, real *
	    , real *);
    extern doublereal fourpi_(real *);

    /* Fortran I/O blocks */
    static cilist io___1581 = { 0, 99, 0, 0, 0 };


/*     PURPOSE:                                                         * */
/*           DEALING WITH PION+N-->L/S+KAON PROCESS AND PION PRODUCTION * */
/*     NOTE   :                                                         * */

/*     QUANTITIES:                                                 * */
/*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME* */
/*           SRT      - SQRT OF S                                       * */
/*           IBLOCK   - THE INFORMATION BACK                            * */
/*                     7  PION+N-->L/S+KAON */
/*           iblock   - 77 pion+N-->Delta+pion */
/*           iblock   - 78 pion+N-->Delta+RHO */
/*           iblock   - 79 pion+N-->Delta+OMEGA */
/*           iblock   - 222 pion+N-->Phi */
/* ********************************* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /input1/ */
/* c      SAVE /RNDF77/ */
    px0 = *px;
    py0 = *py;
    pz0 = *pz;
    *iblock = 1;
    x1 = ranart_(&rndf77_1.nseed);
    ianti = 0;
    if (ee_1.lb[*i1 - 1] < 0 || ee_1.lb[*i2 - 1] < 0) {
	ianti = 1;
    }
    if (*xkaon0 / (*xkaon + *xphi) >= x1) {
/* kaon production */
/* ----------------------------------------------------------------------- */
	*iblock = 7;
	if (ianti == 1) {
	    *iblock = -7;
	}
	ntag = 0;
/* RELABLE PARTICLES FOR THE PROCESS PION+n-->LAMBDA K OR SIGMA k */
/* DECIDE LAMBDA OR SIGMA PRODUCTION, AND TO CALCULATE THE NEW */
/* MOMENTA FOR PARTICLES IN THE FINAL STATE. */
	kaonc = 0;
	if (pnlka_(srt) / (pnlka_(srt) + pnska_(srt)) > ranart_(&
		rndf77_1.nseed)) {
	    kaonc = 1;
	}
	if (cc_1.e[*i1 - 1] <= .2f) {
	    ee_1.lb[*i1 - 1] = 23;
	    cc_1.e[*i1 - 1] = .498f;
	    if (kaonc == 1) {
		ee_1.lb[*i2 - 1] = 14;
		cc_1.e[*i2 - 1] = 1.1157f;
	    } else {
		ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 
			15;
		cc_1.e[*i2 - 1] = 1.1974f;
	    }
	    if (ianti == 1) {
		ee_1.lb[*i1 - 1] = 21;
		ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
	    }
	} else {
	    ee_1.lb[*i2 - 1] = 23;
	    cc_1.e[*i2 - 1] = .498f;
	    if (kaonc == 1) {
		ee_1.lb[*i1 - 1] = 14;
		cc_1.e[*i1 - 1] = 1.1157f;
	    } else {
		ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 
			15;
		cc_1.e[*i1 - 1] = 1.1974f;
	    }
	    if (ianti == 1) {
		ee_1.lb[*i2 - 1] = 21;
		ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
	    }
	}
	em1 = cc_1.e[*i1 - 1];
	em2 = cc_1.e[*i2 - 1];
	goto L50;
/* to gererate the momentum for the kaon and L/S */
    } else if (*xphi / (*xkaon + *xphi) >= x1) {
	*iblock = 222;
	if (*xphin / *xphi >= ranart_(&rndf77_1.nseed)) {
	    ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	    cc_1.e[*i1 - 1] = .939457f;
	} else {
	    ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 4) + 6;
	    cc_1.e[*i1 - 1] = 1.232f;
	}
/*  !! at present only baryon */
	if (ianti == 1) {
	    ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
	}
	ee_1.lb[*i2 - 1] = 29;
	cc_1.e[*i2 - 1] = 1.02f;
	em1 = cc_1.e[*i1 - 1];
	em2 = cc_1.e[*i2 - 1];
	goto L50;
    } else {
/* CHECK WHAT KIND OF PION PRODUCTION PROCESS HAS HAPPENED */
	if (ranart_(&rndf77_1.nseed) <= twopi_(srt) / (twopi_(srt) + threpi_(
		srt) + fourpi_(srt))) {
	    *iblock = 77;
	} else {
	    if (threpi_(srt) / (threpi_(srt) + fourpi_(srt)) > ranart_(&
		    rndf77_1.nseed)) {
		*iblock = 78;
	    } else {
		*iblock = 79;
	    }
	}
	ntag = 0;
/* pion production (Delta+pion/rho/omega in the final state) */
/* generate the mass of the delta resonance */
	x2 = ranart_(&rndf77_1.nseed);
/* relable the particles */
	if (*iblock == 77) {
/* GENERATE THE DELTA MASS */
	    dmax__ = *srt - .13496f - .02f;
	    dm = rmass_(&dmax__, &input1_1.iseed);
/* pion+baryon-->pion+delta */
/* Relable particles, I1 is assigned to the Delta and I2 is assigned to the */
/* meson */
/* (1) for pi(+)+p-->D(+)+P(+) OR D(++)+p(0) */
	    if (ee_1.lb[*i1 - 1] == 1 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 
		    - 1] == 5 && ee_1.lb[*i2 - 1] == 1 || (ee_1.lb[*i1 - 1] ==
		     -1 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && 
		    ee_1.lb[*i2 - 1] == -1)) {
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1) {
		    ii = *i1;
		    if (x2 <= .5f) {
			ee_1.lb[*i1 - 1] = 8;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 5;
			cc_1.e[*i2 - 1] = .13496f;
			goto L40;
		    } else {
			ee_1.lb[*i1 - 1] = 9;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 4;
			ipi = 4;
			cc_1.e[*i2 - 1] = .13496f;
			goto L40;
		    }
		} else {
		    ii = *i2;
		    if (x2 <= .5f) {
			ee_1.lb[*i2 - 1] = 8;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 5;
			cc_1.e[*i1 - 1] = .13496f;
			goto L40;
		    } else {
			ee_1.lb[*i2 - 1] = 9;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 4;
			cc_1.e[*i1 - 1] = .13496f;
			goto L40;
		    }
		}
	    }
/* (2) for pi(-)+p-->D(0)+P(0) OR D(+)+p(-),or D(-)+p(+) */
	    if (ee_1.lb[*i1 - 1] == 1 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 
		    - 1] == 3 && ee_1.lb[*i2 - 1] == 1 || (ee_1.lb[*i1 - 1] ==
		     -1 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && 
		    ee_1.lb[*i2 - 1] == -1)) {
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1) {
		    ii = *i1;
		    if (x2 <= .33f) {
			ee_1.lb[*i1 - 1] = 6;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 5;
			cc_1.e[*i2 - 1] = .13496f;
			goto L40;
		    }
		    if (x2 > .33f && x2 <= .67f) {
			ee_1.lb[*i1 - 1] = 7;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 4;
			cc_1.e[*i2 - 1] = .13496f;
			goto L40;
		    }
		    if (x2 > .67f) {
			ee_1.lb[*i1 - 1] = 8;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 3;
			cc_1.e[*i2 - 1] = .13496f;
			goto L40;
		    }
		} else {
		    ii = *i2;
		    if (x2 <= .33f) {
			ee_1.lb[*i2 - 1] = 6;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 5;
			cc_1.e[*i1 - 1] = .13496f;
			goto L40;
		    }
		    if (x2 > .33f && x2 <= .67f) {
			ee_1.lb[*i2 - 1] = 7;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 4;
			cc_1.e[*i1 - 1] = .13496f;
			goto L40;
		    }
		    if (x2 > .67f) {
			ee_1.lb[*i2 - 1] = 8;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 3;
			cc_1.e[*i1 - 1] = .13496f;
			goto L40;
		    }
		}
	    }
/* (3) for pi(+)+n-->D(+)+Pi(0) OR D(++)+p(-) or D(0)+pi(+) */
	    if (ee_1.lb[*i1 - 1] == 2 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 
		    - 1] == 5 && ee_1.lb[*i2 - 1] == 2 || (ee_1.lb[*i1 - 1] ==
		     -2 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && 
		    ee_1.lb[*i2 - 1] == -2)) {
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2) {
		    ii = *i1;
		    if (x2 <= .33f) {
			ee_1.lb[*i1 - 1] = 8;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 4;
			cc_1.e[*i2 - 1] = .13496f;
			goto L40;
		    }
		    if (x2 > .33f && x2 <= .67f) {
			ee_1.lb[*i1 - 1] = 7;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 5;
			cc_1.e[*i2 - 1] = .13496f;
			goto L40;
		    }
		    if (x2 > .67f) {
			ee_1.lb[*i1 - 1] = 9;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 3;
			cc_1.e[*i2 - 1] = .13496f;
			goto L40;
		    }
		} else {
		    ii = *i2;
		    if (x2 <= .33f) {
			ee_1.lb[*i2 - 1] = 8;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 4;
			cc_1.e[*i1 - 1] = .13496f;
			goto L40;
		    }
		    if (x2 > .33f && x2 <= .67f) {
			ee_1.lb[*i2 - 1] = 7;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 5;
			cc_1.e[*i1 - 1] = .13496f;
			goto L40;
		    }
		    if (x2 > .67f) {
			ee_1.lb[*i2 - 1] = 9;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 3;
			cc_1.e[*i1 - 1] = .13496f;
			goto L40;
		    }
		}
	    }
/* (4) for pi(0)+p-->D(+)+Pi(0) OR D(++)+p(-) or D(0)+pi(+) */
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1 && ee_1.lb[*i2 - 1] 
		    == 4 || ee_1.lb[*i1 - 1] == 4 && (i__2 = ee_1.lb[*i2 - 1],
		     abs(i__2)) == 1) {
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1) {
		    ii = *i1;
		    if (x2 <= .33f) {
			ee_1.lb[*i1 - 1] = 8;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 4;
			cc_1.e[*i2 - 1] = .13496f;
			goto L40;
		    }
		    if (x2 > .33f && x2 <= .67f) {
			ee_1.lb[*i1 - 1] = 7;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 5;
			cc_1.e[*i2 - 1] = .13496f;
			goto L40;
		    }
		    if (x2 > .67f) {
			ee_1.lb[*i1 - 1] = 9;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 3;
			cc_1.e[*i2 - 1] = .13496f;
			goto L40;
		    }
		} else {
		    ii = *i2;
		    if (x2 <= .33f) {
			ee_1.lb[*i2 - 1] = 8;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 4;
			cc_1.e[*i1 - 1] = .13496f;
			goto L40;
		    }
		    if (x2 > .33f && x2 <= .67f) {
			ee_1.lb[*i2 - 1] = 7;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 5;
			cc_1.e[*i1 - 1] = .13496f;
			goto L40;
		    }
		    if (x2 > .67f) {
			ee_1.lb[*i2 - 1] = 9;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 3;
			cc_1.e[*i1 - 1] = .13496f;
			goto L40;
		    }
		}
	    }
/* (5) for pi(-)+n-->D(-)+P(0) OR D(0)+p(-) */
	    if (ee_1.lb[*i1 - 1] == 2 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 
		    - 1] == 3 && ee_1.lb[*i2 - 1] == 2 || (ee_1.lb[*i1 - 1] ==
		     -2 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && 
		    ee_1.lb[*i2 - 1] == -2)) {
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2) {
		    ii = *i1;
		    if (x2 <= .5f) {
			ee_1.lb[*i1 - 1] = 6;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 4;
			cc_1.e[*i2 - 1] = .13496f;
			goto L40;
		    } else {
			ee_1.lb[*i1 - 1] = 7;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 3;
			cc_1.e[*i2 - 1] = .13496f;
			goto L40;
		    }
		} else {
		    ii = *i2;
		    if (x2 <= .5f) {
			ee_1.lb[*i2 - 1] = 6;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 4;
			cc_1.e[*i1 - 1] = .13496f;
			goto L40;
		    } else {
			ee_1.lb[*i2 - 1] = 7;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 3;
			cc_1.e[*i1 - 1] = .13496f;
			goto L40;
		    }
		}
	    }
/* (6) for pi(0)+n-->D(0)+P(0), D(-)+p(+) or D(+)+p(-) */
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 && ee_1.lb[*i2 - 1] 
		    == 4 || ee_1.lb[*i1 - 1] == 4 && (i__2 = ee_1.lb[*i2 - 1],
		     abs(i__2)) == 2) {
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2) {
		    ii = *i1;
		    if (x2 <= .33f) {
			ee_1.lb[*i1 - 1] = 7;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 4;
			cc_1.e[*i2 - 1] = .13496f;
			goto L40;
		    }
		    if (x2 <= .67f && x2 > .33f) {
			ee_1.lb[*i1 - 1] = 6;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 5;
			cc_1.e[*i2 - 1] = .13496f;
			goto L40;
		    }
		    if (x2 > .67f) {
			ee_1.lb[*i1 - 1] = 8;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 3;
			cc_1.e[*i2 - 1] = .13496f;
			goto L40;
		    }
		} else {
		    ii = *i2;
		    if (x2 <= .33f) {
			ee_1.lb[*i2 - 1] = 7;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 4;
			cc_1.e[*i1 - 1] = .13496f;
			goto L40;
		    }
		    if (x2 <= .67f && x2 > .33f) {
			ee_1.lb[*i2 - 1] = 6;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 5;
			cc_1.e[*i1 - 1] = .13496f;
			goto L40;
		    }
		    if (x2 > .67f) {
			ee_1.lb[*i2 - 1] = 8;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 3;
			cc_1.e[*i1 - 1] = .13496f;
			goto L40;
		    }
		}
	    }
	}
	if (*iblock == 78) {
	    rmasdd_(srt, &c_b195, &c_b725, &c_b197, &c_b727, &input1_1.iseed, 
		    &c__4, &dm, &ameson);
	    arho = ameson;
/* pion+baryon-->Rho+delta */
/* (1) for pi(+)+p-->D(+)+rho(+) OR D(++)+rho(0) */
	    if (ee_1.lb[*i1 - 1] == 1 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 
		    - 1] == 5 && ee_1.lb[*i2 - 1] == 1 || (ee_1.lb[*i1 - 1] ==
		     -1 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && 
		    ee_1.lb[*i2 - 1] == -1)) {
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1) {
		    ii = *i1;
		    if (x2 <= .5f) {
			ee_1.lb[*i1 - 1] = 8;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 27;
			cc_1.e[*i2 - 1] = arho;
			goto L40;
		    } else {
			ee_1.lb[*i1 - 1] = 9;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 26;
			cc_1.e[*i2 - 1] = arho;
			goto L40;
		    }
		} else {
		    ii = *i2;
		    if (x2 <= .5f) {
			ee_1.lb[*i2 - 1] = 8;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 27;
			cc_1.e[*i1 - 1] = arho;
			goto L40;
		    } else {
			ee_1.lb[*i2 - 1] = 9;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 26;
			cc_1.e[*i1 - 1] = arho;
			goto L40;
		    }
		}
	    }
/* (2) for pi(-)+p-->D(+)+rho(-) OR D(0)+rho(0) or D(-)+rho(+) */
	    if (ee_1.lb[*i1 - 1] == 1 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 
		    - 1] == 3 && ee_1.lb[*i2 - 1] == 1 || (ee_1.lb[*i1 - 1] ==
		     -1 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && 
		    ee_1.lb[*i2 - 1] == -1)) {
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1) {
		    ii = *i1;
		    if (x2 <= .33f) {
			ee_1.lb[*i1 - 1] = 6;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 27;
			cc_1.e[*i2 - 1] = arho;
			goto L40;
		    }
		    if (x2 > .33f && x2 <= .67f) {
			ee_1.lb[*i1 - 1] = 7;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 26;
			cc_1.e[*i2 - 1] = arho;
			goto L40;
		    }
		    if (x2 > .67f) {
			ee_1.lb[*i1 - 1] = 8;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 25;
			cc_1.e[*i2 - 1] = arho;
			goto L40;
		    }
		} else {
		    ii = *i2;
		    if (x2 <= .33f) {
			ee_1.lb[*i2 - 1] = 6;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 27;
			cc_1.e[*i1 - 1] = arho;
			goto L40;
		    }
		    if (x2 > .33f && x2 <= .67f) {
			ee_1.lb[*i2 - 1] = 7;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 26;
			cc_1.e[*i1 - 1] = arho;
			goto L40;
		    }
		    if (x2 > .67f) {
			ee_1.lb[*i2 - 1] = 8;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 25;
			cc_1.e[*i1 - 1] = arho;
			goto L40;
		    }
		}
	    }
/* (3) for pi(+)+n-->D(+)+rho(0) OR D(++)+rho(-) or D(0)+rho(+) */
	    if (ee_1.lb[*i1 - 1] == 2 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 
		    - 1] == 5 && ee_1.lb[*i2 - 1] == 2 || (ee_1.lb[*i1 - 1] ==
		     -2 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && 
		    ee_1.lb[*i2 - 1] == -2)) {
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2) {
		    ii = *i1;
		    if (x2 <= .33f) {
			ee_1.lb[*i1 - 1] = 8;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 26;
			cc_1.e[*i2 - 1] = arho;
			goto L40;
		    }
		    if (x2 > .33f && x2 <= .67f) {
			ee_1.lb[*i1 - 1] = 7;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 27;
			cc_1.e[*i2 - 1] = arho;
			goto L40;
		    }
		    if (x2 > .67f) {
			ee_1.lb[*i1 - 1] = 9;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 25;
			cc_1.e[*i2 - 1] = arho;
			goto L40;
		    }
		} else {
		    ii = *i2;
		    if (x2 <= .33f) {
			ee_1.lb[*i2 - 1] = 8;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 26;
			cc_1.e[*i1 - 1] = arho;
			goto L40;
		    }
		    if (x2 > .33f && x2 <= .67f) {
			ee_1.lb[*i2 - 1] = 7;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 27;
			cc_1.e[*i1 - 1] = arho;
			goto L40;
		    }
		    if (x2 > .67f) {
			ee_1.lb[*i2 - 1] = 9;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 25;
			cc_1.e[*i1 - 1] = arho;
			goto L40;
		    }
		}
	    }
/* (4) for pi(0)+p-->D(+)+rho(0) OR D(++)+rho(-) or D(0)+rho(+) */
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1 && ee_1.lb[*i2 - 1] 
		    == 4 || ee_1.lb[*i1 - 1] == 4 && (i__2 = ee_1.lb[*i2 - 1],
		     abs(i__2)) == 1) {
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1) {
		    ii = *i1;
		    if (x2 <= .33f) {
			ee_1.lb[*i1 - 1] = 7;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 27;
			cc_1.e[*i2 - 1] = arho;
			goto L40;
		    }
		    if (x2 > .33f && x2 <= .67f) {
			ee_1.lb[*i1 - 1] = 8;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 26;
			cc_1.e[*i2 - 1] = arho;
			goto L40;
		    }
		    if (x2 > .67f) {
			ee_1.lb[*i1 - 1] = 9;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 25;
			cc_1.e[*i2 - 1] = arho;
			goto L40;
		    }
		} else {
		    ii = *i2;
		    if (x2 <= .33f) {
			ee_1.lb[*i2 - 1] = 7;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 27;
			cc_1.e[*i1 - 1] = arho;
			goto L40;
		    }
		    if (x2 > .33f && x2 <= .67f) {
			ee_1.lb[*i2 - 1] = 8;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 26;
			cc_1.e[*i1 - 1] = arho;
			goto L40;
		    }
		    if (x2 > .67f) {
			ee_1.lb[*i2 - 1] = 9;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 25;
			cc_1.e[*i1 - 1] = arho;
			goto L40;
		    }
		}
	    }
/* (5) for pi(-)+n-->D(-)+rho(0) OR D(0)+rho(-) */
	    if (ee_1.lb[*i1 - 1] == 2 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 
		    - 1] == 3 && ee_1.lb[*i2 - 1] == 2 || (ee_1.lb[*i1 - 1] ==
		     -2 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && 
		    ee_1.lb[*i2 - 1] == -2)) {
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2) {
		    ii = *i1;
		    if (x2 <= .5f) {
			ee_1.lb[*i1 - 1] = 6;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 26;
			cc_1.e[*i2 - 1] = arho;
			goto L40;
		    } else {
			ee_1.lb[*i1 - 1] = 7;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 25;
			cc_1.e[*i2 - 1] = arho;
			goto L40;
		    }
		} else {
		    ii = *i2;
		    if (x2 <= .5f) {
			ee_1.lb[*i2 - 1] = 6;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 26;
			cc_1.e[*i1 - 1] = arho;
			goto L40;
		    } else {
			ee_1.lb[*i2 - 1] = 7;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 25;
			cc_1.e[*i1 - 1] = arho;
			goto L40;
		    }
		}
	    }
/* (6) for pi(0)+n-->D(0)+rho(0), D(-)+rho(+) and D(+)+rho(-) */
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 && ee_1.lb[*i2 - 1] 
		    == 4 || ee_1.lb[*i1 - 1] == 4 && (i__2 = ee_1.lb[*i2 - 1],
		     abs(i__2)) == 2) {
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2) {
		    ii = *i1;
		    if (x2 <= .33f) {
			ee_1.lb[*i1 - 1] = 7;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 26;
			cc_1.e[*i2 - 1] = arho;
			goto L40;
		    }
		    if (x2 > .33f && x2 <= .67f) {
			ee_1.lb[*i1 - 1] = 6;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 27;
			cc_1.e[*i2 - 1] = arho;
			goto L40;
		    }
		    if (x2 > .67f) {
			ee_1.lb[*i1 - 1] = 8;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 25;
			cc_1.e[*i2 - 1] = arho;
		    }
		} else {
		    ii = *i2;
		    if (x2 <= .33f) {
			ee_1.lb[*i2 - 1] = 7;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 26;
			cc_1.e[*i1 - 1] = arho;
			goto L40;
		    }
		    if (x2 <= .67f && x2 > .33f) {
			ee_1.lb[*i2 - 1] = 6;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 27;
			cc_1.e[*i1 - 1] = arho;
			goto L40;
		    }
		    if (x2 > .67f) {
			ee_1.lb[*i2 - 1] = 8;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 25;
			cc_1.e[*i1 - 1] = arho;
		    }
		}
	    }
	}
	if (*iblock == 79) {
	    aomega = .782f;
/* GENERATE THE DELTA MASS */
	    dmax__ = *srt - .782f - .02f;
	    dm = rmass_(&dmax__, &input1_1.iseed);
/* pion+baryon-->omega+delta */
/* (1) for pi(+)+p-->D(++)+omega(0) */
	    if (ee_1.lb[*i1 - 1] == 1 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 
		    - 1] == 5 && ee_1.lb[*i2 - 1] == 1 || (ee_1.lb[*i1 - 1] ==
		     -1 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && 
		    ee_1.lb[*i2 - 1] == -1)) {
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1) {
		    ii = *i1;
		    ee_1.lb[*i1 - 1] = 9;
		    cc_1.e[*i1 - 1] = dm;
		    ee_1.lb[*i2 - 1] = 28;
		    cc_1.e[*i2 - 1] = aomega;
		    goto L40;
		} else {
		    ii = *i2;
		    ee_1.lb[*i2 - 1] = 9;
		    cc_1.e[*i2 - 1] = dm;
		    ee_1.lb[*i1 - 1] = 28;
		    cc_1.e[*i1 - 1] = aomega;
		    goto L40;
		}
	    }
/* (2) for pi(-)+p-->D(0)+omega(0) */
	    if (ee_1.lb[*i1 - 1] == 1 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 
		    - 1] == 3 && ee_1.lb[*i2 - 1] == 1 || (ee_1.lb[*i1 - 1] ==
		     -1 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && 
		    ee_1.lb[*i2 - 1] == -1)) {
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1) {
		    ii = *i1;
		    ee_1.lb[*i1 - 1] = 7;
		    cc_1.e[*i1 - 1] = dm;
		    ee_1.lb[*i2 - 1] = 28;
		    cc_1.e[*i2 - 1] = aomega;
		    goto L40;
		} else {
		    ii = *i2;
		    ee_1.lb[*i2 - 1] = 7;
		    cc_1.e[*i2 - 1] = dm;
		    ee_1.lb[*i1 - 1] = 28;
		    cc_1.e[*i1 - 1] = aomega;
		    goto L40;
		}
	    }
/* (3) for pi(+)+n-->D(+)+omega(0) */
	    if (ee_1.lb[*i1 - 1] == 2 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 
		    - 1] == 5 && ee_1.lb[*i2 - 1] == 2 || (ee_1.lb[*i1 - 1] ==
		     -2 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && 
		    ee_1.lb[*i2 - 1] == -2)) {
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2) {
		    ii = *i1;
		    ee_1.lb[*i1 - 1] = 8;
		    cc_1.e[*i1 - 1] = dm;
		    ee_1.lb[*i2 - 1] = 28;
		    cc_1.e[*i2 - 1] = aomega;
		    goto L40;
		} else {
		    ii = *i2;
		    ee_1.lb[*i2 - 1] = 8;
		    cc_1.e[*i2 - 1] = dm;
		    ee_1.lb[*i1 - 1] = 28;
		    cc_1.e[*i1 - 1] = aomega;
		    goto L40;
		}
	    }
/* (4) for pi(0)+p-->D(+)+omega(0) */
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1 && ee_1.lb[*i2 - 1] 
		    == 4 || ee_1.lb[*i1 - 1] == 4 && (i__2 = ee_1.lb[*i2 - 1],
		     abs(i__2)) == 1) {
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1) {
		    ii = *i1;
		    ee_1.lb[*i1 - 1] = 8;
		    cc_1.e[*i1 - 1] = dm;
		    ee_1.lb[*i2 - 1] = 28;
		    cc_1.e[*i2 - 1] = aomega;
		    goto L40;
		} else {
		    ii = *i2;
		    ee_1.lb[*i2 - 1] = 8;
		    cc_1.e[*i2 - 1] = dm;
		    ee_1.lb[*i1 - 1] = 28;
		    cc_1.e[*i1 - 1] = aomega;
		    goto L40;
		}
	    }
/* (5) for pi(-)+n-->D(-)+omega(0) */
	    if (ee_1.lb[*i1 - 1] == 2 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 
		    - 1] == 3 && ee_1.lb[*i2 - 1] == 2 || (ee_1.lb[*i1 - 1] ==
		     -2 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && 
		    ee_1.lb[*i2 - 1] == -2)) {
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2) {
		    ii = *i1;
		    ee_1.lb[*i1 - 1] = 6;
		    cc_1.e[*i1 - 1] = dm;
		    ee_1.lb[*i2 - 1] = 28;
		    cc_1.e[*i2 - 1] = aomega;
		    goto L40;
		} else {
		    ii = *i2;
		    ee_1.lb[*i2 - 1] = 6;
		    cc_1.e[*i2 - 1] = dm;
		    ee_1.lb[*i1 - 1] = 28;
		    cc_1.e[*i1 - 1] = aomega;
		}
	    }
/* (6) for pi(0)+n-->D(0)+omega(0) */
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 && ee_1.lb[*i2 - 1] 
		    == 4 || ee_1.lb[*i1 - 1] == 4 && (i__2 = ee_1.lb[*i2 - 1],
		     abs(i__2)) == 2) {
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2) {
		    ii = *i1;
		    ee_1.lb[*i1 - 1] = 7;
		    cc_1.e[*i1 - 1] = dm;
		    ee_1.lb[*i2 - 1] = 28;
		    cc_1.e[*i2 - 1] = aomega;
		    goto L40;
		} else {
		    ii = *i2;
		    ee_1.lb[*i2 - 1] = 7;
		    cc_1.e[*i2 - 1] = dm;
		    ee_1.lb[*i1 - 1] = 26;
		    cc_1.e[*i1 - 1] = arho;
		    goto L40;
		}
	    }
	}
L40:
	em1 = cc_1.e[*i1 - 1];
	em2 = cc_1.e[*i2 - 1];
	if (ianti == 1 && ee_1.lb[*i1 - 1] >= 1 && ee_1.lb[*i2 - 1] >= 1) {
	    ee_1.lb[ii - 1] = -ee_1.lb[ii - 1];
	    jj = *i2;
	    if (ii == *i2) {
		jj = *i1;
	    }
	    if (*iblock == 77) {
		if (ee_1.lb[jj - 1] == 3) {
		    ee_1.lb[jj - 1] = 5;
		} else if (ee_1.lb[jj - 1] == 5) {
		    ee_1.lb[jj - 1] = 3;
		}
	    } else if (*iblock == 78) {
		if (ee_1.lb[jj - 1] == 25) {
		    ee_1.lb[jj - 1] = 27;
		} else if (ee_1.lb[jj - 1] == 27) {
		    ee_1.lb[jj - 1] = 25;
		}
	    }
	}
    }
/* ----------------------------------------------------------------------- */
/* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH */
/* ENERGY CONSERVATION */
L50:
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = em1;
/* Computing 2nd power */
    r__4 = em2;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = em1 * em2;
    pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
    if (pr2 <= 0.f) {
	pr2 = 1e-8f;
    }
    pr = sqrt(pr2) / (*srt * 2.f);
/* here we use the same transverse momentum distribution as for */
/* pp collisions, it might be necessary to use a different distribution */
/* lin-10/25/02 get rid of argument usage mismatch in PTR(): */
    xptr = pr * .33f;
/*         cc1=ptr(0.33*pr,iseed) */
    cc1 = ptr_(&xptr, &input1_1.iseed);
/* lin-10/25/02-end */
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__1 = pr;
/* Computing 2nd power */
    r__2 = cc1;
    scheck = r__1 * r__1 - r__2 * r__2;
    if (scheck < 0.f) {
	s_wsle(&io___1581);
	do_lio(&c__9, &c__1, "scheck36: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	scheck = 0.f;
    }
    c1 = sqrt(scheck) / pr;
/*         c1=sqrt(pr**2-cc1**2)/pr */
/*          C1   = 1.0 - 2.0 * RANART(NSEED) */
    t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
/* Computing 2nd power */
    r__1 = c1;
    s1 = sqrt(1.f - r__1 * r__1);
    ct1 = cos(t1);
    st1 = sin(t1);
/* THE MOMENTUM IN THE CMS IN THE FINAL STATE */
    *pz = pr * c1;
    *px = pr * s1 * ct1;
    *py = pr * s1 * st1;
/* ROTATE IT */
    rotate_(&px0, &py0, &pz0, px, py, pz);
    return 0;
} /* crpn_ */

/* ********************************* */
/*                                                                      * */
/*                                                                      * */
/* Subroutine */ int cren_(real *px, real *py, real *pz, real *srt, integer *
	i1, integer *i2, integer *iblock)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    double sqrt(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static real c1, s1, t1, pr, em1, em2, ct1, pr2, px0, py0, pz0, st1;
    static integer ntag, kaonc, ianti;
    extern doublereal pnlka_(real *), pnska_(real *), ranart_(integer *);

/*     PURPOSE:                                                         * */
/*             DEALING WITH ETA+N-->L/S+KAON PROCESS                   * */
/*     NOTE   :                                                         * */

/*     QUANTITIES:                                                 * */
/*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME* */
/*           SRT      - SQRT OF S                                       * */
/*           IBLOCK   - THE INFORMATION BACK                            * */
/*                     7  ETA+N-->L/S+KAON */
/* ********************************* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /input1/ */
/* c      SAVE /RNDF77/ */
    px0 = *px;
    py0 = *py;
    pz0 = *pz;
    ntag = 0;
    *iblock = 7;
    ianti = 0;
    if (ee_1.lb[*i1 - 1] < 0 || ee_1.lb[*i2 - 1] < 0) {
	ianti = 1;
	*iblock = -7;
    }
/* RELABLE PARTICLES FOR THE PROCESS eta+n-->LAMBDA K OR SIGMA k */
/* DECIDE LAMBDA OR SIGMA PRODUCTION, AND TO CALCULATE THE NEW */
/* MOMENTA FOR PARTICLES IN THE FINAL STATE. */
    kaonc = 0;
    if (pnlka_(srt) / (pnlka_(srt) + pnska_(srt)) > ranart_(&rndf77_1.nseed)) 
	    {
	kaonc = 1;
    }
    if (cc_1.e[*i1 - 1] <= .6f) {
	ee_1.lb[*i1 - 1] = 23;
	cc_1.e[*i1 - 1] = .498f;
	if (kaonc == 1) {
	    ee_1.lb[*i2 - 1] = 14;
	    cc_1.e[*i2 - 1] = 1.1157f;
	} else {
	    ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 15;
	    cc_1.e[*i2 - 1] = 1.1974f;
	}
	if (ianti == 1) {
	    ee_1.lb[*i1 - 1] = 21;
	    ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
	}
    } else {
	ee_1.lb[*i2 - 1] = 23;
	cc_1.e[*i2 - 1] = .498f;
	if (kaonc == 1) {
	    ee_1.lb[*i1 - 1] = 14;
	    cc_1.e[*i1 - 1] = 1.1157f;
	} else {
	    ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 15;
	    cc_1.e[*i1 - 1] = 1.1974f;
	}
	if (ianti == 1) {
	    ee_1.lb[*i2 - 1] = 21;
	    ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
	}
    }
    em1 = cc_1.e[*i1 - 1];
    em2 = cc_1.e[*i2 - 1];
/* ----------------------------------------------------------------------- */
/* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH */
/* ENERGY CONSERVATION */
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = em1;
/* Computing 2nd power */
    r__4 = em2;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = em1 * em2;
    pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
    if (pr2 <= 0.f) {
	pr2 = 1e-9f;
    }
    pr = sqrt(pr2) / (*srt * 2.f);
    c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
    t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
/* Computing 2nd power */
    r__1 = c1;
    s1 = sqrt(1.f - r__1 * r__1);
    ct1 = cos(t1);
    st1 = sin(t1);
/* THE MOMENTUM IN THE CMS IN THE FINAL STATE */
    *pz = pr * c1;
    *px = pr * s1 * ct1;
    *py = pr * s1 * st1;
/* FOR THE ISOTROPIC DISTRIBUTION THERE IS NO NEED TO ROTATE */
    return 0;
} /* cren_ */

/* ********************************* */
/*                                                                      * */
/*                                                                      * */
/*      SUBROUTINE Crdir(PX,PY,PZ,SRT,I1,I2) */
/* Subroutine */ int crdir_(real *px, real *py, real *pz, real *srt, integer *
	i1, integer *i2, integer *iblock)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static real c1, s1, t1, pr, cc1, em1, em2, ct1, pr2, px0, py0, pz0, st1;
    extern doublereal ptr_(real *, integer *);
    static integer ntag;
    static real xptr, scheck;
    extern doublereal ranart_(integer *);
    extern /* Subroutine */ int rotate_(real *, real *, real *, real *, real *
	    , real *);

    /* Fortran I/O blocks */
    static cilist io___1613 = { 0, 99, 0, 0, 0 };


/*     PURPOSE:                                                         * */
/*             DEALING WITH pion+N-->pion+N PROCESS                   * */
/*     NOTE   :                                                         * */

/*     QUANTITIES:                                                 * */
/*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME* */
/*           SRT      - SQRT OF S                                       * */
/*           IBLOCK   - THE INFORMATION BACK                            * */

/* ********************************* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /input1/ */
/* c      SAVE /RNDF77/ */
    px0 = *px;
    py0 = *py;
    pz0 = *pz;
    *iblock = 999;
    ntag = 0;
    em1 = cc_1.e[*i1 - 1];
    em2 = cc_1.e[*i2 - 1];
/* ----------------------------------------------------------------------- */
/* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH */
/* ENERGY CONSERVATION */
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = em1;
/* Computing 2nd power */
    r__4 = em2;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = em1 * em2;
    pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
    if (pr2 <= 0.f) {
	pr2 = 1e-9f;
    }
    pr = sqrt(pr2) / (*srt * 2.f);
/* lin-10/25/02 get rid of argument usage mismatch in PTR(): */
    xptr = pr * .33f;
/*         cc1=ptr(0.33*pr,iseed) */
    cc1 = ptr_(&xptr, &input1_1.iseed);
/* lin-10/25/02-end */
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__1 = pr;
/* Computing 2nd power */
    r__2 = cc1;
    scheck = r__1 * r__1 - r__2 * r__2;
    if (scheck < 0.f) {
	s_wsle(&io___1613);
	do_lio(&c__9, &c__1, "scheck37: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	scheck = 0.f;
    }
    c1 = sqrt(scheck) / pr;
/*         c1=sqrt(pr**2-cc1**2)/pr */
    t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
/* Computing 2nd power */
    r__1 = c1;
    s1 = sqrt(1.f - r__1 * r__1);
    ct1 = cos(t1);
    st1 = sin(t1);
/* THE MOMENTUM IN THE CMS IN THE FINAL STATE */
    *pz = pr * c1;
    *px = pr * s1 * ct1;
    *py = pr * s1 * st1;
/* ROTATE the momentum */
    rotate_(&px0, &py0, &pz0, px, py, pz);
    return 0;
} /* crdir_ */

/* ********************************* */
/*                                                                      * */
/*                                                                      * */
/* Subroutine */ int crpd_(real *px, real *py, real *pz, real *srt, integer *
	i1, integer *i2, integer *iblock, real *xkaon0, real *xkaon, real *
	xphi, real *xphin)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    real r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static real c1, s1, t1, x1, x2;
    static integer ii, jj;
    static real pr, cc1, em1, em2, ct1, pr2, px0, py0, pz0, st1;
    extern doublereal ptr_(real *, integer *);
    static integer ntag;
    static real xptr;
    static integer kaonc, ianti;
    extern doublereal pnlka_(real *), pnska_(real *);
    static real scheck;
    extern doublereal ranart_(integer *);
    extern /* Subroutine */ int rotate_(real *, real *, real *, real *, real *
	    , real *);

    /* Fortran I/O blocks */
    static cilist io___1636 = { 0, 99, 0, 0, 0 };


/*     PURPOSE:                                                         * */
/*     DEALING WITH PION+D(N*)-->PION +N OR */
/*                                             L/S+KAON PROCESS         * */
/*     NOTE   :                                                         * */

/*     QUANTITIES:                                                 * */
/*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME* */
/*           SRT      - SQRT OF S                                       * */
/*           IBLOCK   - THE INFORMATION BACK                            * */
/*                     7  PION+D(N*)-->L/S+KAON */
/*           iblock   - 80 pion+D(N*)-->pion+N */
/*           iblock   - 81 RHO+D(N*)-->PION+N */
/*           iblock   - 82 OMEGA+D(N*)-->PION+N */
/*                     222  PION+D --> PHI */
/* ********************************* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /input1/ */
/* c      SAVE /RNDF77/ */
    px0 = *px;
    py0 = *py;
    pz0 = *pz;
    *iblock = 1;
    x1 = ranart_(&rndf77_1.nseed);
    ianti = 0;
    if (ee_1.lb[*i1 - 1] < 0 || ee_1.lb[*i2 - 1] < 0) {
	ianti = 1;
    }
    if (*xkaon0 / (*xkaon + *xphi) >= x1) {
/* kaon production */
/* ----------------------------------------------------------------------- */
	*iblock = 7;
	if (ianti == 1) {
	    *iblock = -7;
	}
	ntag = 0;
/* RELABLE PARTICLES FOR THE PROCESS PION+n-->LAMBDA K OR SIGMA k */
/* DECIDE LAMBDA OR SIGMA PRODUCTION, AND TO CALCULATE THE NEW */
/* MOMENTA FOR PARTICLES IN THE FINAL STATE. */
	kaonc = 0;
	if (pnlka_(srt) / (pnlka_(srt) + pnska_(srt)) > ranart_(&
		rndf77_1.nseed)) {
	    kaonc = 1;
	}
/* lin-8/17/00     & +PNSKA(SRT)).GT.RANART(NSEED))KAONC=1 */
	if (cc_1.e[*i1 - 1] <= .2f) {
	    ee_1.lb[*i1 - 1] = 23;
	    cc_1.e[*i1 - 1] = .498f;
	    if (kaonc == 1) {
		ee_1.lb[*i2 - 1] = 14;
		cc_1.e[*i2 - 1] = 1.1157f;
	    } else {
		ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 
			15;
		cc_1.e[*i2 - 1] = 1.1974f;
	    }
	    if (ianti == 1) {
		ee_1.lb[*i1 - 1] = 21;
		ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
	    }
	} else {
	    ee_1.lb[*i2 - 1] = 23;
	    cc_1.e[*i2 - 1] = .498f;
	    if (kaonc == 1) {
		ee_1.lb[*i1 - 1] = 14;
		cc_1.e[*i1 - 1] = 1.1157f;
	    } else {
		ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 
			15;
		cc_1.e[*i1 - 1] = 1.1974f;
	    }
	    if (ianti == 1) {
		ee_1.lb[*i2 - 1] = 21;
		ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
	    }
	}
	em1 = cc_1.e[*i1 - 1];
	em2 = cc_1.e[*i2 - 1];
	goto L50;
/* to gererate the momentum for the kaon and L/S */

/* * Phi production */
    } else if (*xphi / (*xkaon + *xphi) >= x1) {
	*iblock = 222;
	if (*xphin / *xphi >= ranart_(&rndf77_1.nseed)) {
	    ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	    cc_1.e[*i1 - 1] = .939457f;
	} else {
	    ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 4) + 6;
	    cc_1.e[*i1 - 1] = 1.232f;
	}
/*   !! at present only baryon */
	if (ianti == 1) {
	    ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
	}
	ee_1.lb[*i2 - 1] = 29;
	cc_1.e[*i2 - 1] = 1.02f;
	em1 = cc_1.e[*i1 - 1];
	em2 = cc_1.e[*i2 - 1];
	goto L50;
    } else {
/* PION REABSORPTION HAS HAPPENED */
	x2 = ranart_(&rndf77_1.nseed);
	*iblock = 80;
	ntag = 0;
/* Relable particles, I1 is assigned to the nucleon */
/* and I2 is assigned to the pion */
/* for the reverse of the following process */
/* (1) for D(+)+P(+)-->p+pion(+) */
	if (ee_1.lb[*i1 - 1] == 8 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1]
		 == 5 && ee_1.lb[*i2 - 1] == 8 || (ee_1.lb[*i1 - 1] == -8 && 
		ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 
		- 1] == -8)) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8) {
		ii = *i1;
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .939457f;
		ee_1.lb[*i2 - 1] = 5;
		cc_1.e[*i2 - 1] = .13496f;
		goto L40;
	    } else {
		ii = *i2;
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .939457f;
		ee_1.lb[*i1 - 1] = 5;
		cc_1.e[*i1 - 1] = .13496f;
		goto L40;
	    }
	}

/* (2) for D(0)+P(0)-->n+pi(0) or p+pi(-) */
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7 && ee_1.lb[*i2 - 1] == 
		4 || ee_1.lb[*i1 - 1] == 4 && (i__2 = ee_1.lb[*i2 - 1], abs(
		i__2)) == 7) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7) {
		ii = *i1;
		if (x2 <= .5f) {
		    ee_1.lb[*i1 - 1] = 2;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 4;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i1 - 1] = 1;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 3;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		}
	    } else {
		ii = *i2;
		if (x2 <= .5f) {
		    ee_1.lb[*i2 - 1] = 2;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 4;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i2 - 1] = 1;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 3;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		}
	    }
	}
/* (3) for D(+)+Pi(0)-->pi(+)+n or pi(0)+p */
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8 && ee_1.lb[*i2 - 1] == 
		4 || ee_1.lb[*i1 - 1] == 4 && (i__2 = ee_1.lb[*i2 - 1], abs(
		i__2)) == 8) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8) {
		ii = *i1;
		if (x2 <= .5f) {
		    ee_1.lb[*i1 - 1] = 2;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 5;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i1 - 1] = 1;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 4;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		}
	    } else {
		ii = *i2;
		if (x2 <= .5f) {
		    ee_1.lb[*i2 - 1] = 2;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 5;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i2 - 1] = 1;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 4;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		}
	    }
	}
/* (4) for D(-)+Pi(0)-->n+pi(-) */
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 6 && ee_1.lb[*i2 - 1] == 
		4 || ee_1.lb[*i1 - 1] == 4 && (i__2 = ee_1.lb[*i2 - 1], abs(
		i__2)) == 6) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 6) {
		ii = *i1;
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
		ee_1.lb[*i2 - 1] = 3;
		cc_1.e[*i2 - 1] = .13496f;
		goto L40;
	    } else {
		ii = *i2;
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
		ee_1.lb[*i1 - 1] = 3;
		cc_1.e[*i1 - 1] = .13496f;
		goto L40;
	    }
	}
/* (5) for D(+)+Pi(-)-->pi(0)+n or pi(-)+p */
	if (ee_1.lb[*i1 - 1] == 8 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1]
		 == 3 && ee_1.lb[*i2 - 1] == 8 || (ee_1.lb[*i1 - 1] == -8 && 
		ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 
		- 1] == -8)) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8) {
		ii = *i1;
		if (x2 <= .5f) {
		    ee_1.lb[*i1 - 1] = 2;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 4;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i1 - 1] = 1;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 3;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		}
	    } else {
		ii = *i2;
		if (x2 <= .5f) {
		    ee_1.lb[*i2 - 1] = 2;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 4;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i2 - 1] = 1;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 3;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		}
	    }
	}
/* (6) D(0)+P(+)-->n+pi(+) or p+pi(0) */
	if (ee_1.lb[*i1 - 1] == 7 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1]
		 == 5 && ee_1.lb[*i2 - 1] == 7 || (ee_1.lb[*i1 - 1] == -7 && 
		ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 
		- 1] == -7)) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7) {
		ii = *i1;
		if (x2 <= .5f) {
		    ee_1.lb[*i1 - 1] = 2;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 5;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i1 - 1] = 1;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 4;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		}
	    } else {
		ii = *i2;
		if (x2 <= .5f) {
		    ee_1.lb[*i2 - 1] = 2;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 5;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i2 - 1] = 1;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 4;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		}
	    }
	}
/* (7) for D(0)+Pi(-)-->n+pi(-) */
	if (ee_1.lb[*i1 - 1] == 7 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1]
		 == 3 && ee_1.lb[*i2 - 1] == 7 || (ee_1.lb[*i1 - 1] == -7 && 
		ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 
		- 1] == -7)) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7) {
		ii = *i1;
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
		ee_1.lb[*i2 - 1] = 3;
		cc_1.e[*i2 - 1] = .13496f;
		goto L40;
	    } else {
		ii = *i2;
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
		ee_1.lb[*i1 - 1] = 3;
		cc_1.e[*i1 - 1] = .13496f;
		goto L40;
	    }
	}
/* (8) D(-)+P(+)-->n+pi(0) or p+pi(-) */
	if (ee_1.lb[*i1 - 1] == 6 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1]
		 == 5 && ee_1.lb[*i2 - 1] == 6 || (ee_1.lb[*i1 - 1] == -6 && 
		ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 
		- 1] == -6)) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 6) {
		ii = *i1;
		if (x2 <= .5f) {
		    ee_1.lb[*i1 - 1] = 2;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 4;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i1 - 1] = 1;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 3;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		}
	    } else {
		ii = *i2;
		if (x2 <= .5f) {
		    ee_1.lb[*i2 - 1] = 2;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 4;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i2 - 1] = 1;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 3;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		}
	    }
	}

/* (9) D(++)+P(-)-->n+pi(+) or p+pi(0) */
	if (ee_1.lb[*i1 - 1] == 9 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1]
		 == 3 && ee_1.lb[*i2 - 1] == 9 || (ee_1.lb[*i1 - 1] == -9 && 
		ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 
		- 1] == -9)) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 9) {
		ii = *i1;
		if (x2 <= .5f) {
		    ee_1.lb[*i1 - 1] = 2;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 5;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i1 - 1] = 1;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 4;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		}
	    } else {
		ii = *i2;
		if (x2 <= .5f) {
		    ee_1.lb[*i2 - 1] = 2;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 5;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i2 - 1] = 1;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 4;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		}
	    }
	}
/* (10) for D(++)+Pi(0)-->p+pi(+) */
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 9 && ee_1.lb[*i2 - 1] == 
		4 || ee_1.lb[*i1 - 1] == 4 && (i__2 = ee_1.lb[*i2 - 1], abs(
		i__2)) == 9) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 9) {
		ii = *i1;
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .939457f;
		ee_1.lb[*i2 - 1] = 5;
		cc_1.e[*i2 - 1] = .13496f;
		goto L40;
	    } else {
		ii = *i2;
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .939457f;
		ee_1.lb[*i1 - 1] = 5;
		cc_1.e[*i1 - 1] = .13496f;
		goto L40;
	    }
	}
/* (11) for N*(1440)(+)or N*(1535)(+)+P(+)-->p+pion(+) */
	if (ee_1.lb[*i1 - 1] == 11 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 
		1] == 5 && ee_1.lb[*i2 - 1] == 11 || ee_1.lb[*i1 - 1] == 13 &&
		 ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*
		i2 - 1] == 13 || (ee_1.lb[*i1 - 1] == -11 && ee_1.lb[*i2 - 1] 
		== 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == -11 || 
		ee_1.lb[*i1 - 1] == -13 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*
		i1 - 1] == 3 && ee_1.lb[*i2 - 1] == -13)) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 11 || (i__2 = ee_1.lb[
		    *i1 - 1], abs(i__2)) == 13) {
		ii = *i1;
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .939457f;
		ee_1.lb[*i2 - 1] = 5;
		cc_1.e[*i2 - 1] = .13496f;
		goto L40;
	    } else {
		ii = *i2;
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .939457f;
		ee_1.lb[*i1 - 1] = 5;
		cc_1.e[*i1 - 1] = .13496f;
		goto L40;
	    }
	}
/* (12) for N*(1440) or N*(1535)(0)+P(0)-->n+pi(0) or p+pi(-) */
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 10 && ee_1.lb[*i2 - 1] == 
		4 || ee_1.lb[*i1 - 1] == 4 && (i__2 = ee_1.lb[*i2 - 1], abs(
		i__2)) == 10 || ee_1.lb[*i1 - 1] == 4 && (i__3 = ee_1.lb[*i2 
		- 1], abs(i__3)) == 12 || ee_1.lb[*i2 - 1] == 4 && (i__4 = 
		ee_1.lb[*i1 - 1], abs(i__4)) == 12) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 10 || (i__2 = ee_1.lb[
		    *i1 - 1], abs(i__2)) == 12) {
		ii = *i1;
		if (x2 <= .5f) {
		    ee_1.lb[*i1 - 1] = 2;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 4;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i1 - 1] = 1;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 3;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		}
	    } else {
		ii = *i2;
		if (x2 <= .5f) {
		    ee_1.lb[*i2 - 1] = 2;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 4;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i2 - 1] = 1;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 3;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		}
	    }
	}
/* (13) for N*(1440) or N*(1535)(+)+Pi(0)-->pi(+)+n or pi(0)+p */
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 11 && ee_1.lb[*i2 - 1] == 
		4 || ee_1.lb[*i1 - 1] == 4 && (i__2 = ee_1.lb[*i2 - 1], abs(
		i__2)) == 11 || ee_1.lb[*i1 - 1] == 4 && (i__3 = ee_1.lb[*i2 
		- 1], abs(i__3)) == 13 || ee_1.lb[*i2 - 1] == 4 && (i__4 = 
		ee_1.lb[*i1 - 1], abs(i__4)) == 13) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 11 || (i__2 = ee_1.lb[
		    *i1 - 1], abs(i__2)) == 13) {
		ii = *i1;
		if (x2 <= .5f) {
		    ee_1.lb[*i1 - 1] = 2;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 5;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i1 - 1] = 1;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 4;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		}
	    } else {
		ii = *i2;
		if (x2 <= .5f) {
		    ee_1.lb[*i2 - 1] = 2;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 5;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i2 - 1] = 1;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 4;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		}
	    }
	}
/* (14) for N*(1440) or N*(1535)(+)+Pi(-)-->pi(0)+n or pi(-)+p */
	if (ee_1.lb[*i1 - 1] == 11 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 
		1] == 3 && ee_1.lb[*i2 - 1] == 11 || ee_1.lb[*i1 - 1] == 3 && 
		ee_1.lb[*i2 - 1] == 13 || ee_1.lb[*i2 - 1] == 3 && ee_1.lb[*
		i1 - 1] == 13 || (ee_1.lb[*i1 - 1] == -11 && ee_1.lb[*i2 - 1] 
		== 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == -11 || 
		ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == -13 || ee_1.lb[*
		i2 - 1] == 5 && ee_1.lb[*i1 - 1] == -13)) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 11 || (i__2 = ee_1.lb[
		    *i1 - 1], abs(i__2)) == 13) {
		ii = *i1;
		if (x2 <= .5f) {
		    ee_1.lb[*i1 - 1] = 2;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 4;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i1 - 1] = 1;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 3;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		}
	    } else {
		ii = *i2;
		if (x2 <= .5f) {
		    ee_1.lb[*i2 - 1] = 2;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 4;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i2 - 1] = 1;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 3;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		}
	    }
	}
/* (15) N*(1440) or N*(1535)(0)+P(+)-->n+pi(+) or p+pi(0) */
	if (ee_1.lb[*i1 - 1] == 10 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 
		1] == 5 && ee_1.lb[*i2 - 1] == 10 || ee_1.lb[*i1 - 1] == 12 &&
		 ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*
		i2 - 1] == 12 || (ee_1.lb[*i1 - 1] == -10 && ee_1.lb[*i2 - 1] 
		== 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == -10 || 
		ee_1.lb[*i1 - 1] == -12 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*
		i1 - 1] == 3 && ee_1.lb[*i2 - 1] == -12)) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 10 || (i__2 = ee_1.lb[
		    *i1 - 1], abs(i__2)) == 12) {
		ii = *i1;
		if (x2 <= .5f) {
		    ee_1.lb[*i1 - 1] = 2;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 5;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i1 - 1] = 1;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 4;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		}
	    } else {
		ii = *i2;
		if (x2 <= .5f) {
		    ee_1.lb[*i2 - 1] = 2;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 5;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i2 - 1] = 1;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 4;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		}
	    }
	}
/* (16) for N*(1440) or N*(1535) (0)+Pi(-)-->n+pi(-) */
	if (ee_1.lb[*i1 - 1] == 10 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 
		1] == 3 && ee_1.lb[*i2 - 1] == 10 || ee_1.lb[*i1 - 1] == 3 && 
		ee_1.lb[*i2 - 1] == 12 || ee_1.lb[*i1 - 1] == 12 && ee_1.lb[*
		i2 - 1] == 3 || (ee_1.lb[*i1 - 1] == -10 && ee_1.lb[*i2 - 1] 
		== 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == -10 || 
		ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == -12 || ee_1.lb[*
		i1 - 1] == -12 && ee_1.lb[*i2 - 1] == 5)) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 10 || (i__2 = ee_1.lb[
		    *i1 - 1], abs(i__2)) == 12) {
		ii = *i1;
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
		ee_1.lb[*i2 - 1] = 3;
		cc_1.e[*i2 - 1] = .13496f;
		goto L40;
	    } else {
		ii = *i2;
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
		ee_1.lb[*i1 - 1] = 3;
		cc_1.e[*i1 - 1] = .13496f;
		goto L40;
	    }
	}
L40:
	em1 = cc_1.e[*i1 - 1];
	em2 = cc_1.e[*i2 - 1];
	if (ianti == 1 && ee_1.lb[*i1 - 1] >= 1 && ee_1.lb[*i2 - 1] >= 1) {
	    ee_1.lb[ii - 1] = -ee_1.lb[ii - 1];
	    jj = *i2;
	    if (ii == *i2) {
		jj = *i1;
	    }
	    if (ee_1.lb[jj - 1] == 3) {
		ee_1.lb[jj - 1] = 5;
	    } else if (ee_1.lb[jj - 1] == 5) {
		ee_1.lb[jj - 1] = 3;
	    }
	}
    }
/* ----------------------------------------------------------------------- */
/* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH */
/* ENERGY CONSERVATION */
L50:
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = em1;
/* Computing 2nd power */
    r__4 = em2;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = em1 * em2;
    pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
    if (pr2 <= 0.f) {
	pr2 = 1e-9f;
    }
    pr = sqrt(pr2) / (*srt * 2.f);
/* lin-10/25/02 get rid of argument usage mismatch in PTR(): */
    xptr = pr * .33f;
/*         cc1=ptr(0.33*pr,iseed) */
    cc1 = ptr_(&xptr, &input1_1.iseed);
/* lin-10/25/02-end */
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__1 = pr;
/* Computing 2nd power */
    r__2 = cc1;
    scheck = r__1 * r__1 - r__2 * r__2;
    if (scheck < 0.f) {
	s_wsle(&io___1636);
	do_lio(&c__9, &c__1, "scheck38: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	scheck = 0.f;
    }
    c1 = sqrt(scheck) / pr;
/*         c1=sqrt(pr**2-cc1**2)/pr */
/*         C1   = 1.0 - 2.0 * RANART(NSEED) */
    t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
/* Computing 2nd power */
    r__1 = c1;
    s1 = sqrt(1.f - r__1 * r__1);
    ct1 = cos(t1);
    st1 = sin(t1);
    *pz = pr * c1;
    *px = pr * s1 * ct1;
    *py = pr * s1 * st1;
/* rotate the momentum */
    rotate_(&px0, &py0, &pz0, px, py, pz);
    return 0;
} /* crpd_ */

/* ********************************* */
/*                                                                      * */
/*                                                                      * */
/* Subroutine */ int crrd_(real *px, real *py, real *pz, real *srt, integer *
	i1, integer *i2, integer *iblock, real *xkaon0, real *xkaon, real *
	xphi, real *xphin)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    real r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static real c1, s1, t1, x1, x2;
    static integer ii, jj;
    static real pr, cc1, em1, em2, ct1, pr2, px0, py0, pz0, st1;
    extern doublereal ptr_(real *, integer *);
    static integer ntag;
    static real xptr;
    static integer kaonc, ianti;
    extern doublereal pnlka_(real *), pnska_(real *);
    static real scheck;
    extern doublereal ranart_(integer *);
    extern /* Subroutine */ int rotate_(real *, real *, real *, real *, real *
	    , real *);

    /* Fortran I/O blocks */
    static cilist io___1659 = { 0, 99, 0, 0, 0 };


/*     PURPOSE:                                                         * */
/*     DEALING WITH rho(omega)+N or D(N*)-->PION +N OR */
/*                                             L/S+KAON PROCESS         * */
/*     NOTE   :                                                         * */

/*     QUANTITIES:                                                 * */
/*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME* */
/*           SRT      - SQRT OF S                                       * */
/*           IBLOCK   - THE INFORMATION BACK                            * */
/*                     7  rho(omega)+N or D(N*)-->L/S+KAON */
/*           iblock   - 80 pion+D(N*)-->pion+N */
/*           iblock   - 81 RHO+D(N*)-->PION+N */
/*           iblock   - 82 OMEGA+D(N*)-->PION+N */
/*           iblock   - 222 pion+N-->Phi */
/* ********************************* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /input1/ */
/* c      SAVE /RNDF77/ */
    px0 = *px;
    py0 = *py;
    pz0 = *pz;
    *iblock = 1;
    ianti = 0;
    if (ee_1.lb[*i1 - 1] < 0 || ee_1.lb[*i2 - 1] < 0) {
	ianti = 1;
    }
    x1 = ranart_(&rndf77_1.nseed);
    if (*xkaon0 / (*xkaon + *xphi) >= x1) {
/* kaon production */
/* ----------------------------------------------------------------------- */
	*iblock = 7;
	if (ianti == 1) {
	    *iblock = -7;
	}
	ntag = 0;
/* RELABLE PARTICLES FOR THE PROCESS PION+n-->LAMBDA K OR SIGMA k */
/* DECIDE LAMBDA OR SIGMA PRODUCTION, AND TO CALCULATE THE NEW */
/* MOMENTA FOR PARTICLES IN THE FINAL STATE. */
	kaonc = 0;
	if (pnlka_(srt) / (pnlka_(srt) + pnska_(srt)) > ranart_(&
		rndf77_1.nseed)) {
	    kaonc = 1;
	}
/* lin-8/17/00     & +PNSKA(SRT)).GT.RANART(NSEED))KAONC=1 */
	if (cc_1.e[*i1 - 1] <= .92f) {
	    ee_1.lb[*i1 - 1] = 23;
	    cc_1.e[*i1 - 1] = .498f;
	    if (kaonc == 1) {
		ee_1.lb[*i2 - 1] = 14;
		cc_1.e[*i2 - 1] = 1.1157f;
	    } else {
		ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 
			15;
		cc_1.e[*i2 - 1] = 1.1974f;
	    }
	    if (ianti == 1) {
		ee_1.lb[*i1 - 1] = 21;
		ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
	    }
	} else {
	    ee_1.lb[*i2 - 1] = 23;
	    cc_1.e[*i2 - 1] = .498f;
	    if (kaonc == 1) {
		ee_1.lb[*i1 - 1] = 14;
		cc_1.e[*i1 - 1] = 1.1157f;
	    } else {
		ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 
			15;
		cc_1.e[*i1 - 1] = 1.1974f;
	    }
	    if (ianti == 1) {
		ee_1.lb[*i2 - 1] = 21;
		ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
	    }
	}
	em1 = cc_1.e[*i1 - 1];
	em2 = cc_1.e[*i2 - 1];
	goto L50;
/* to gererate the momentum for the kaon and L/S */

/* * Phi production */
    } else if (*xphi / (*xkaon + *xphi) >= x1) {
	*iblock = 222;
	if (*xphin / *xphi >= ranart_(&rndf77_1.nseed)) {
	    ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	    cc_1.e[*i1 - 1] = .939457f;
	} else {
	    ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 4) + 6;
	    cc_1.e[*i1 - 1] = 1.232f;
	}
/*   !! at present only baryon */
	if (ianti == 1) {
	    ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
	}
	ee_1.lb[*i2 - 1] = 29;
	cc_1.e[*i2 - 1] = 1.02f;
	em1 = cc_1.e[*i1 - 1];
	em2 = cc_1.e[*i2 - 1];
	goto L50;
    } else {
/* rho(omega) REABSORPTION HAS HAPPENED */
	x2 = ranart_(&rndf77_1.nseed);
	*iblock = 81;
	ntag = 0;
	if (ee_1.lb[*i1 - 1] == 28 || ee_1.lb[*i2 - 1] == 28) {
	    goto L60;
	}
/* we treat Rho reabsorption in the following */
/* Relable particles, I1 is assigned to the Delta */
/* and I2 is assigned to the meson */
/* for the reverse of the following process */
/* (1) for D(+)+rho(+)-->p+pion(+) */
	if (ee_1.lb[*i1 - 1] == 8 && ee_1.lb[*i2 - 1] == 27 || ee_1.lb[*i1 - 
		1] == 27 && ee_1.lb[*i2 - 1] == 8 || (ee_1.lb[*i1 - 1] == -8 
		&& ee_1.lb[*i2 - 1] == 25 || ee_1.lb[*i1 - 1] == 25 && 
		ee_1.lb[*i2 - 1] == -8)) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8) {
		ii = *i1;
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .939457f;
		ee_1.lb[*i2 - 1] = 5;
		cc_1.e[*i2 - 1] = .13496f;
		goto L40;
	    } else {
		ii = *i2;
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .939457f;
		ee_1.lb[*i1 - 1] = 5;
		cc_1.e[*i1 - 1] = .13496f;
		goto L40;
	    }
	}
/* (2) for D(0)+rho(0)-->n+pi(0) or p+pi(-) */
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7 && ee_1.lb[*i2 - 1] == 
		26 || ee_1.lb[*i1 - 1] == 26 && (i__2 = ee_1.lb[*i2 - 1], abs(
		i__2)) == 7) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7) {
		ii = *i1;
		if (x2 <= .5f) {
		    ee_1.lb[*i1 - 1] = 2;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 4;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i1 - 1] = 1;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 3;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		}
	    } else {
		ii = *i2;
		if (x2 <= .5f) {
		    ee_1.lb[*i2 - 1] = 2;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 4;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i2 - 1] = 1;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 3;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		}
	    }
	}
/* (3) for D(+)+rho(0)-->pi(+)+n or pi(0)+p */
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8 && ee_1.lb[*i2 - 1] == 
		26 || ee_1.lb[*i1 - 1] == 26 && (i__2 = ee_1.lb[*i2 - 1], abs(
		i__2)) == 8) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8) {
		ii = *i1;
		if (x2 <= .5f) {
		    ee_1.lb[*i1 - 1] = 2;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 5;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i1 - 1] = 1;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 4;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		}
	    } else {
		ii = *i2;
		if (x2 <= .5f) {
		    ee_1.lb[*i2 - 1] = 2;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 5;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i2 - 1] = 1;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 4;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		}
	    }
	}
/* (4) for D(-)+rho(0)-->n+pi(-) */
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 6 && ee_1.lb[*i2 - 1] == 
		26 || ee_1.lb[*i1 - 1] == 26 && (i__2 = ee_1.lb[*i2 - 1], abs(
		i__2)) == 6) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 6) {
		ii = *i1;
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
		ee_1.lb[*i2 - 1] = 3;
		cc_1.e[*i2 - 1] = .13496f;
		goto L40;
	    } else {
		ii = *i2;
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
		ee_1.lb[*i1 - 1] = 3;
		cc_1.e[*i1 - 1] = .13496f;
		goto L40;
	    }
	}
/* (5) for D(+)+rho(-)-->pi(0)+n or pi(-)+p */
	if (ee_1.lb[*i1 - 1] == 8 && ee_1.lb[*i2 - 1] == 25 || ee_1.lb[*i1 - 
		1] == 25 && ee_1.lb[*i2 - 1] == 8 || (ee_1.lb[*i1 - 1] == -8 
		&& ee_1.lb[*i2 - 1] == 27 || ee_1.lb[*i1 - 1] == 27 && 
		ee_1.lb[*i2 - 1] == -8)) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8) {
		ii = *i1;
		if (x2 <= .5f) {
		    ee_1.lb[*i1 - 1] = 2;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 4;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i1 - 1] = 1;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 3;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		}
	    } else {
		ii = *i2;
		if (x2 <= .5f) {
		    ee_1.lb[*i2 - 1] = 2;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 4;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i2 - 1] = 1;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 3;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		}
	    }
	}
/* (6) D(0)+rho(+)-->n+pi(+) or p+pi(0) */
	if (ee_1.lb[*i1 - 1] == 7 && ee_1.lb[*i2 - 1] == 27 || ee_1.lb[*i1 - 
		1] == 27 && ee_1.lb[*i2 - 1] == 7 || (ee_1.lb[*i1 - 1] == -7 
		&& ee_1.lb[*i2 - 1] == 25 || ee_1.lb[*i1 - 1] == 25 && 
		ee_1.lb[*i2 - 1] == -7)) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7) {
		ii = *i1;
		if (x2 <= .5f) {
		    ee_1.lb[*i1 - 1] = 2;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 5;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i1 - 1] = 1;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 4;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		}
	    } else {
		ii = *i2;
		if (x2 <= .5f) {
		    ee_1.lb[*i2 - 1] = 2;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 5;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i2 - 1] = 1;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 4;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		}
	    }
	}
/* (7) for D(0)+rho(-)-->n+pi(-) */
	if (ee_1.lb[*i1 - 1] == 7 && ee_1.lb[*i2 - 1] == 25 || ee_1.lb[*i1 - 
		1] == 25 && ee_1.lb[*i2 - 1] == 7 || (ee_1.lb[*i1 - 1] == -7 
		&& ee_1.lb[*i2 - 1] == 27 || ee_1.lb[*i1 - 1] == 27 && 
		ee_1.lb[*i2 - 1] == -7)) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7) {
		ii = *i1;
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
		ee_1.lb[*i2 - 1] = 3;
		cc_1.e[*i2 - 1] = .13496f;
		goto L40;
	    } else {
		ii = *i2;
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
		ee_1.lb[*i1 - 1] = 3;
		cc_1.e[*i1 - 1] = .13496f;
		goto L40;
	    }
	}
/* (8) D(-)+rho(+)-->n+pi(0) or p+pi(-) */
	if (ee_1.lb[*i1 - 1] == 6 && ee_1.lb[*i2 - 1] == 27 || ee_1.lb[*i1 - 
		1] == 27 && ee_1.lb[*i2 - 1] == 6 || (ee_1.lb[*i1 - 1] == -6 
		&& ee_1.lb[*i2 - 1] == 25 || ee_1.lb[*i1 - 1] == 25 && 
		ee_1.lb[*i2 - 1] == -6)) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 6) {
		ii = *i1;
		if (x2 <= .5f) {
		    ee_1.lb[*i1 - 1] = 2;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 4;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i1 - 1] = 1;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 3;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		}
	    } else {
		ii = *i2;
		if (x2 <= .5f) {
		    ee_1.lb[*i2 - 1] = 2;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 4;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i2 - 1] = 1;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 3;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		}
	    }
	}
/* (9) D(++)+rho(-)-->n+pi(+) or p+pi(0) */
	if (ee_1.lb[*i1 - 1] == 9 && ee_1.lb[*i2 - 1] == 25 || ee_1.lb[*i1 - 
		1] == 25 && ee_1.lb[*i2 - 1] == 9 || (ee_1.lb[*i1 - 1] == -9 
		&& ee_1.lb[*i2 - 1] == 27 || ee_1.lb[*i1 - 1] == 27 && 
		ee_1.lb[*i2 - 1] == -9)) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 9) {
		ii = *i1;
		if (x2 <= .5f) {
		    ee_1.lb[*i1 - 1] = 2;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 5;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i1 - 1] = 1;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 4;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		}
	    } else {
		ii = *i2;
		if (x2 <= .5f) {
		    ee_1.lb[*i2 - 1] = 2;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 5;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i2 - 1] = 1;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 4;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		}
	    }
	}
/* (10) for D(++)+rho(0)-->p+pi(+) */
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 9 && ee_1.lb[*i2 - 1] == 
		26 || ee_1.lb[*i1 - 1] == 26 && (i__2 = ee_1.lb[*i2 - 1], abs(
		i__2)) == 9) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 9) {
		ii = *i1;
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .939457f;
		ee_1.lb[*i2 - 1] = 5;
		cc_1.e[*i2 - 1] = .13496f;
		goto L40;
	    } else {
		ii = *i2;
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .939457f;
		ee_1.lb[*i1 - 1] = 5;
		cc_1.e[*i1 - 1] = .13496f;
		goto L40;
	    }
	}
/* (11) for N*(1440)(+)or N*(1535)(+)+rho(+)-->p+pion(+) */
	if (ee_1.lb[*i1 - 1] == 11 && ee_1.lb[*i2 - 1] == 27 || ee_1.lb[*i1 - 
		1] == 27 && ee_1.lb[*i2 - 1] == 11 || ee_1.lb[*i1 - 1] == 13 
		&& ee_1.lb[*i2 - 1] == 27 || ee_1.lb[*i1 - 1] == 27 && 
		ee_1.lb[*i2 - 1] == 13 || (ee_1.lb[*i1 - 1] == -11 && ee_1.lb[
		*i2 - 1] == 25 || ee_1.lb[*i1 - 1] == 25 && ee_1.lb[*i2 - 1] 
		== -11 || ee_1.lb[*i1 - 1] == -13 && ee_1.lb[*i2 - 1] == 25 ||
		 ee_1.lb[*i1 - 1] == 25 && ee_1.lb[*i2 - 1] == -13)) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 11 || (i__2 = ee_1.lb[
		    *i1 - 1], abs(i__2)) == 13) {
		ii = *i1;
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .939457f;
		ee_1.lb[*i2 - 1] = 5;
		cc_1.e[*i2 - 1] = .13496f;
		goto L40;
	    } else {
		ii = *i2;
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .939457f;
		ee_1.lb[*i1 - 1] = 5;
		cc_1.e[*i1 - 1] = .13496f;
		goto L40;
	    }
	}
/* (12) for N*(1440) or N*(1535)(0)+rho(0)-->n+pi(0) or p+pi(-) */
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 10 && ee_1.lb[*i2 - 1] == 
		26 || ee_1.lb[*i1 - 1] == 26 && (i__2 = ee_1.lb[*i2 - 1], abs(
		i__2)) == 10 || ee_1.lb[*i1 - 1] == 26 && (i__3 = ee_1.lb[*i2 
		- 1], abs(i__3)) == 12 || ee_1.lb[*i2 - 1] == 26 && (i__4 = 
		ee_1.lb[*i1 - 1], abs(i__4)) == 12) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 10 || (i__2 = ee_1.lb[
		    *i1 - 1], abs(i__2)) == 12) {
		ii = *i1;
		if (x2 <= .5f) {
		    ee_1.lb[*i1 - 1] = 2;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 4;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i1 - 1] = 1;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 3;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		}
	    } else {
		ii = *i2;
		if (x2 <= .5f) {
		    ee_1.lb[*i2 - 1] = 2;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 4;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i2 - 1] = 1;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 3;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		}
	    }
	}
/* (13) for N*(1440) or N*(1535)(+)+rho(0)-->pi(+)+n or pi(0)+p */
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 11 && ee_1.lb[*i2 - 1] == 
		26 || ee_1.lb[*i1 - 1] == 26 && (i__2 = ee_1.lb[*i2 - 1], abs(
		i__2)) == 11 || ee_1.lb[*i1 - 1] == 26 && (i__3 = ee_1.lb[*i2 
		- 1], abs(i__3)) == 13 || ee_1.lb[*i2 - 1] == 26 && (i__4 = 
		ee_1.lb[*i1 - 1], abs(i__4)) == 13) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 11 || (i__2 = ee_1.lb[
		    *i1 - 1], abs(i__2)) == 13) {
		ii = *i1;
		if (x2 <= .5f) {
		    ee_1.lb[*i1 - 1] = 2;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 5;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i1 - 1] = 1;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 4;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		}
	    } else {
		ii = *i2;
		if (x2 <= .5f) {
		    ee_1.lb[*i2 - 1] = 2;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 5;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i2 - 1] = 1;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 4;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		}
	    }
	}
/* (14) for N*(1440) or N*(1535)(+)+rho(-)-->pi(0)+n or pi(-)+p */
	if (ee_1.lb[*i1 - 1] == 11 && ee_1.lb[*i2 - 1] == 25 || ee_1.lb[*i1 - 
		1] == 25 && ee_1.lb[*i2 - 1] == 11 || ee_1.lb[*i1 - 1] == 25 
		&& ee_1.lb[*i2 - 1] == 13 || ee_1.lb[*i2 - 1] == 25 && 
		ee_1.lb[*i1 - 1] == 13 || (ee_1.lb[*i1 - 1] == -11 && ee_1.lb[
		*i2 - 1] == 27 || ee_1.lb[*i1 - 1] == 27 && ee_1.lb[*i2 - 1] 
		== -11 || ee_1.lb[*i1 - 1] == 27 && ee_1.lb[*i2 - 1] == -13 ||
		 ee_1.lb[*i2 - 1] == 27 && ee_1.lb[*i1 - 1] == -13)) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 11 || (i__2 = ee_1.lb[
		    *i1 - 1], abs(i__2)) == 13) {
		ii = *i1;
		if (x2 <= .5f) {
		    ee_1.lb[*i1 - 1] = 2;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 4;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i1 - 1] = 1;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 3;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		}
	    } else {
		ii = *i2;
		if (x2 <= .5f) {
		    ee_1.lb[*i2 - 1] = 2;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 4;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i2 - 1] = 1;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 3;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		}
	    }
	}
/* (15) N*(1440) or N*(1535)(0)+rho(+)-->n+pi(+) or p+pi(0) */
	if (ee_1.lb[*i1 - 1] == 10 && ee_1.lb[*i2 - 1] == 27 || ee_1.lb[*i1 - 
		1] == 27 && ee_1.lb[*i2 - 1] == 10 || ee_1.lb[*i1 - 1] == 12 
		&& ee_1.lb[*i2 - 1] == 27 || ee_1.lb[*i1 - 1] == 27 && 
		ee_1.lb[*i2 - 1] == 12 || (ee_1.lb[*i1 - 1] == -10 && ee_1.lb[
		*i2 - 1] == 25 || ee_1.lb[*i1 - 1] == 25 && ee_1.lb[*i2 - 1] 
		== -10 || ee_1.lb[*i1 - 1] == -12 && ee_1.lb[*i2 - 1] == 25 ||
		 ee_1.lb[*i1 - 1] == 25 && ee_1.lb[*i2 - 1] == -12)) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 10 || (i__2 = ee_1.lb[
		    *i1 - 1], abs(i__2)) == 12) {
		ii = *i1;
		if (x2 <= .5f) {
		    ee_1.lb[*i1 - 1] = 2;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 5;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i1 - 1] = 1;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 4;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		}
	    } else {
		ii = *i2;
		if (x2 <= .5f) {
		    ee_1.lb[*i2 - 1] = 2;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 5;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i2 - 1] = 1;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 4;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		}
	    }
	}
/* (16) for N*(1440) or N*(1535) (0)+rho(-)-->n+pi(-) */
	if (ee_1.lb[*i1 - 1] == 10 && ee_1.lb[*i2 - 1] == 25 || ee_1.lb[*i1 - 
		1] == 25 && ee_1.lb[*i2 - 1] == 10 || ee_1.lb[*i1 - 1] == 25 
		&& ee_1.lb[*i2 - 1] == 12 || ee_1.lb[*i1 - 1] == 12 && 
		ee_1.lb[*i2 - 1] == 25 || (ee_1.lb[*i1 - 1] == -10 && ee_1.lb[
		*i2 - 1] == 27 || ee_1.lb[*i1 - 1] == 27 && ee_1.lb[*i2 - 1] 
		== -10 || ee_1.lb[*i1 - 1] == 27 && ee_1.lb[*i2 - 1] == -12 ||
		 ee_1.lb[*i1 - 1] == -12 && ee_1.lb[*i2 - 1] == 27)) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 10 || (i__2 = ee_1.lb[
		    *i1 - 1], abs(i__2)) == 12) {
		ii = *i1;
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
		ee_1.lb[*i2 - 1] = 3;
		cc_1.e[*i2 - 1] = .13496f;
		goto L40;
	    } else {
		ii = *i2;
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
		ee_1.lb[*i1 - 1] = 3;
		cc_1.e[*i1 - 1] = .13496f;
		goto L40;
	    }
	}
L60:
	*iblock = 82;
/* FOR OMEGA REABSORPTION */
/* Relable particles, I1 is assigned to the Delta */
/* and I2 is assigned to the meson */
/* for the reverse of the following process */
/* (1) for D(0)+OMEGA(0)-->n+pi(0) or p+pi(-) */
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7 && ee_1.lb[*i2 - 1] == 
		28 || ee_1.lb[*i1 - 1] == 28 && (i__2 = ee_1.lb[*i2 - 1], abs(
		i__2)) == 7) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7) {
		ii = *i1;
		if (x2 <= .5f) {
		    ee_1.lb[*i1 - 1] = 2;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 4;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i1 - 1] = 1;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 3;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		}
	    } else {
		ii = *i2;
		if (x2 <= .5f) {
		    ee_1.lb[*i2 - 1] = 2;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 4;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i2 - 1] = 1;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 3;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		}
	    }
	}
/* (2) for D(+)+OMEGA(0)-->pi(+)+n or pi(0)+p */
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8 && ee_1.lb[*i2 - 1] == 
		28 || ee_1.lb[*i1 - 1] == 28 && (i__2 = ee_1.lb[*i2 - 1], abs(
		i__2)) == 8) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8) {
		ii = *i1;
		if (x2 <= .5f) {
		    ee_1.lb[*i1 - 1] = 2;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 5;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i1 - 1] = 1;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 4;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		}
	    } else {
		ii = *i2;
		if (x2 <= .5f) {
		    ee_1.lb[*i2 - 1] = 2;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 5;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i2 - 1] = 1;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 4;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		}
	    }
	}
/* (3) for D(-)+OMEGA(0)-->n+pi(-) */
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 6 && ee_1.lb[*i2 - 1] == 
		28 || ee_1.lb[*i1 - 1] == 28 && (i__2 = ee_1.lb[*i2 - 1], abs(
		i__2)) == 6) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 6) {
		ii = *i1;
		ee_1.lb[*i1 - 1] = 2;
		cc_1.e[*i1 - 1] = .939457f;
		ee_1.lb[*i2 - 1] = 3;
		cc_1.e[*i2 - 1] = .13496f;
		goto L40;
	    } else {
		ii = *i2;
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
		ee_1.lb[*i1 - 1] = 3;
		cc_1.e[*i1 - 1] = .13496f;
		goto L40;
	    }
	}
/* (4) for D(++)+OMEGA(0)-->p+pi(+) */
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 9 && ee_1.lb[*i2 - 1] == 
		28 || ee_1.lb[*i1 - 1] == 28 && (i__2 = ee_1.lb[*i2 - 1], abs(
		i__2)) == 9) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 9) {
		ii = *i1;
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .939457f;
		ee_1.lb[*i2 - 1] = 5;
		cc_1.e[*i2 - 1] = .13496f;
		goto L40;
	    } else {
		ii = *i2;
		ee_1.lb[*i2 - 1] = 1;
		cc_1.e[*i2 - 1] = .939457f;
		ee_1.lb[*i1 - 1] = 5;
		cc_1.e[*i1 - 1] = .13496f;
		goto L40;
	    }
	}
/* (5) for N*(1440) or N*(1535)(0)+omega(0)-->n+pi(0) or p+pi(-) */
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 10 && ee_1.lb[*i2 - 1] == 
		28 || ee_1.lb[*i1 - 1] == 28 && (i__2 = ee_1.lb[*i2 - 1], abs(
		i__2)) == 10 || ee_1.lb[*i1 - 1] == 28 && (i__3 = ee_1.lb[*i2 
		- 1], abs(i__3)) == 12 || ee_1.lb[*i2 - 1] == 28 && (i__4 = 
		ee_1.lb[*i1 - 1], abs(i__4)) == 12) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 10 || (i__2 = ee_1.lb[
		    *i1 - 1], abs(i__2)) == 12) {
		ii = *i1;
		if (x2 <= .5f) {
		    ee_1.lb[*i1 - 1] = 2;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 4;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i1 - 1] = 1;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 3;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		}
	    } else {
		ii = *i2;
		if (x2 <= .5f) {
		    ee_1.lb[*i2 - 1] = 2;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 4;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i2 - 1] = 1;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 3;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		}
	    }
	}
/* (6) for N*(1440) or N*(1535)(+)+omega(0)-->pi(+)+n or pi(0)+p */
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 11 && ee_1.lb[*i2 - 1] == 
		28 || ee_1.lb[*i1 - 1] == 28 && (i__2 = ee_1.lb[*i2 - 1], abs(
		i__2)) == 11 || ee_1.lb[*i1 - 1] == 28 && (i__3 = ee_1.lb[*i2 
		- 1], abs(i__3)) == 13 || ee_1.lb[*i2 - 1] == 28 && (i__4 = 
		ee_1.lb[*i1 - 1], abs(i__4)) == 13) {
	    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 11 || (i__2 = ee_1.lb[
		    *i1 - 1], abs(i__2)) == 13) {
		ii = *i1;
		if (x2 <= .5f) {
		    ee_1.lb[*i1 - 1] = 2;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 5;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i1 - 1] = 1;
		    cc_1.e[*i1 - 1] = .939457f;
		    ee_1.lb[*i2 - 1] = 4;
		    cc_1.e[*i2 - 1] = .13496f;
		    goto L40;
		}
	    } else {
		ii = *i2;
		if (x2 <= .5f) {
		    ee_1.lb[*i2 - 1] = 2;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 5;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		} else {
		    ee_1.lb[*i2 - 1] = 1;
		    cc_1.e[*i2 - 1] = .939457f;
		    ee_1.lb[*i1 - 1] = 4;
		    cc_1.e[*i1 - 1] = .13496f;
		    goto L40;
		}
	    }
	}
L40:
	em1 = cc_1.e[*i1 - 1];
	em2 = cc_1.e[*i2 - 1];
	if (ianti == 1 && ee_1.lb[*i1 - 1] >= 1 && ee_1.lb[*i2 - 1] >= 1) {
	    ee_1.lb[ii - 1] = -ee_1.lb[ii - 1];
	    jj = *i2;
	    if (ii == *i2) {
		jj = *i1;
	    }
	    if (ee_1.lb[jj - 1] == 3) {
		ee_1.lb[jj - 1] = 5;
	    } else if (ee_1.lb[jj - 1] == 5) {
		ee_1.lb[jj - 1] = 3;
	    }
	}
    }
/* ----------------------------------------------------------------------- */
/* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH */
/* ENERGY CONSERVATION */
L50:
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = em1;
/* Computing 2nd power */
    r__4 = em2;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = em1 * em2;
    pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
    if (pr2 <= 0.f) {
	pr2 = 1e-9f;
    }
    pr = sqrt(pr2) / (*srt * 2.f);
/*          C1   = 1.0 - 2.0 * RANART(NSEED) */
/* lin-10/25/02 get rid of argument usage mismatch in PTR(): */
    xptr = pr * .33f;
/*         cc1=ptr(0.33*pr,iseed) */
    cc1 = ptr_(&xptr, &input1_1.iseed);
/* lin-10/25/02-end */
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__1 = pr;
/* Computing 2nd power */
    r__2 = cc1;
    scheck = r__1 * r__1 - r__2 * r__2;
    if (scheck < 0.f) {
	s_wsle(&io___1659);
	do_lio(&c__9, &c__1, "scheck39: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	scheck = 0.f;
    }
    c1 = sqrt(scheck) / pr;
/*         c1=sqrt(pr**2-cc1**2)/pr */
    t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
/* Computing 2nd power */
    r__1 = c1;
    s1 = sqrt(1.f - r__1 * r__1);
    ct1 = cos(t1);
    st1 = sin(t1);
    *pz = pr * c1;
    *px = pr * s1 * ct1;
    *py = pr * s1 * st1;
/* ROTATE THE MOMENTUM */
    rotate_(&px0, &py0, &pz0, px, py, pz);
    return 0;
} /* crrd_ */

/* ********************************* */
/* sp 03/19/01                                                          * */
/*                                                                      * */
/* Subroutine */ int crlaba_(real *px, real *py, real *pz, real *srt, real *
	brel, real *brsgm, integer *i1, integer *i2, integer *nt, integer *
	iblock, integer *nchrg, integer *icase)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    double sqrt(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static real c1, s1, t1, pr, em1, em2, ct1, pr2, px0, py0, pz0, st1, rrr;
    extern doublereal ranart_(integer *);
    extern /* Subroutine */ int rotate_(real *, real *, real *, real *, real *
	    , real *);

/*     PURPOSE:                                                         * */
/*            DEALING WITH   K+ + N(D,N*)-bar <-->  La(Si)-bar + pi     * */
/*     NOTE   :                                                         * */
/*                                                                      * */
/*     QUANTITIES:                                                 * */
/*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME* */
/*           SRT      - SQRT OF S                                       * */
/*           IBLOCK   - THE INFORMATION BACK                            * */
/*                     8-> elastic scatt                               * */
/*                     100-> K+ + N-bar -> Sigma-bar + PI */
/*                     102-> PI + Sigma(Lambda)-bar -> K+ + N-bar */
/* ********************************* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /input1/ */
/* c      SAVE /RNDF77/ */

    px0 = *px;
    py0 = *py;
    pz0 = *pz;

    if (*icase == 3) {
	rrr = ranart_(&rndf77_1.nseed);
	if (rrr < *brel) {
/*            !! elastic scat.  (avoid in reverse process) */
	    *iblock = 8;
	} else {
	    *iblock = 100;
	    if (rrr < *brel + *brsgm) {
/* *    K+ + N-bar -> Sigma-bar + PI */
		ee_1.lb[*i1 - 1] = -15 - (integer) (ranart_(&rndf77_1.nseed) *
			 3);
		cc_1.e[*i1 - 1] = 1.1974f;
	    } else {
/* *    K+ + N-bar -> Lambda-bar + PI */
		ee_1.lb[*i1 - 1] = -14;
		cc_1.e[*i1 - 1] = 1.1157f;
	    }
	    ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	    cc_1.e[*i2 - 1] = .138f;
	}
    }


    if (*icase == 4) {
	rrr = ranart_(&rndf77_1.nseed);
	if (rrr < *brel) {
/*            !! elastic scat. */
	    *iblock = 8;
	} else {
	    *iblock = 102;
/*    PI + Sigma(Lambda)-bar -> K+ + N-bar */
/*         ! K+ */
	    ee_1.lb[*i1 - 1] = 23;
	    ee_1.lb[*i2 - 1] = -1 - (integer) (ranart_(&rndf77_1.nseed) * 2);
	    if (*nchrg == -2) {
		ee_1.lb[*i2 - 1] = -6;
	    }
	    if (*nchrg == 1) {
		ee_1.lb[*i2 - 1] = -9;
	    }
	    cc_1.e[*i1 - 1] = .498f;
	    cc_1.e[*i2 - 1] = .938f;
	    if (*nchrg == -2 || *nchrg == 1) {
		cc_1.e[*i2 - 1] = 1.232f;
	    }
	}
    }

    em1 = cc_1.e[*i1 - 1];
    em2 = cc_1.e[*i2 - 1];
/* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH */
/* ENERGY CONSERVATION */
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = em1;
/* Computing 2nd power */
    r__4 = em2;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = em1 * em2;
    pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
    if (pr2 <= 0.f) {
	pr2 = 1e-9f;
    }
    pr = sqrt(pr2) / (*srt * 2.f);
    c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
    t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
/* Computing 2nd power */
    r__1 = c1;
    s1 = sqrt(1.f - r__1 * r__1);
    ct1 = cos(t1);
    st1 = sin(t1);
    *pz = pr * c1;
    *px = pr * s1 * ct1;
    *py = pr * s1 * st1;
/* ROTATE IT */
    rotate_(&px0, &py0, &pz0, px, py, pz);
    return 0;
} /* crlaba_ */

/* ********************************* */
/*                                                                      * */
/*                                                                      * */
/* Subroutine */ int crkn_(real *px, real *py, real *pz, real *srt, integer *
	i1, integer *i2, integer *iblock)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    double sqrt(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static real c1, s1, t1, pr, em1, em2, ct1, pr2, px0, py0, pz0, st1;
    static integer ntag;
    extern doublereal ranart_(integer *);

/*     PURPOSE:                                                         * */
/*             DEALING WITH kaON+N/pi-->KAON +N/pi elastic PROCESS      * */
/*     NOTE   :                                                         * */

/*     QUANTITIES:                                                 * */
/*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME* */
/*           SRT      - SQRT OF S                                       * */
/*           IBLOCK   - THE INFORMATION BACK                            * */
/*                     8-> PION+N-->L/S+KAON */
/* ********************************* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /input1/ */
/* c      SAVE /RNDF77/ */
    px0 = *px;
    py0 = *py;
    pz0 = *pz;
/* ----------------------------------------------------------------------- */
    *iblock = 8;
    ntag = 0;
    em1 = cc_1.e[*i1 - 1];
    em2 = cc_1.e[*i2 - 1];
/* ----------------------------------------------------------------------- */
/* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH */
/* ENERGY CONSERVATION */
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = em1;
/* Computing 2nd power */
    r__4 = em2;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = em1 * em2;
    pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
    if (pr2 <= 0.f) {
	pr2 = 1e-9f;
    }
    pr = sqrt(pr2) / (*srt * 2.f);
    c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
    t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
/* Computing 2nd power */
    r__1 = c1;
    s1 = sqrt(1.f - r__1 * r__1);
    ct1 = cos(t1);
    st1 = sin(t1);
    *pz = pr * c1;
    *px = pr * s1 * ct1;
    *py = pr * s1 * st1;
    return 0;
} /* crkn_ */

/* ********************************* */
/*                                                                      * */
/*                                                                      * */
/* Subroutine */ int crppba_(real *px, real *py, real *pz, real *srt, integer 
	*i1, integer *i2, integer *iblock)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    double sqrt(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static real c1, s1, t1, pr, em1, em2, ct1, pr2, px0, py0, pz0, st1;
    static integer npion, nchrg1, nchrg2;
    static real pmass1, pmass2;
    extern /* Subroutine */ int pbarfs_(real *, integer *, integer *);
    extern doublereal ranart_(integer *);
    extern /* Subroutine */ int rotate_(real *, real *, real *, real *, real *
	    , real *);

/*     PURPOSE:                                                         * */
/* lin-8/29/00*             DEALING WITH anti-nucleon annihilation with */
/*             DEALING WITH anti-baryon annihilation with */
/*             nucleons or baryon resonances */
/*             Determine:                                               * */
/*             (1) no. of pions in the final state */
/*             (2) relable particles in the final state */
/*             (3) new momenta of final state particles                 * */

/*     QUANTITIES:                                                      * */
/*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME* */
/*           SRT      - SQRT OF S                                       * */
/*           IBLOCK   - INFORMATION about the reaction channel          * */

/*           iblock   - 1902 annihilation-->pion(+)+pion(-)   (2 pion) */
/*           iblock   - 1903 annihilation-->pion(+)+rho(-)    (3 pion) */
/*           iblock   - 1904 annihilation-->rho(+)+rho(-)     (4 pion) */
/*           iblock   - 1905 annihilation-->rho(0)+omega      (5 pion) */
/*           iblock   - 1906 annihilation-->omega+omega       (6 pion) */
/*       charge conservation is enforced in relabling particles */
/*       in the final state (note: at the momentum we don't check the */
/*       initial charges while dealing with annihilation, since some */
/*       annihilation channels between antinucleons and nucleons (baryon */
/*       resonances) might be forbiden by charge conservation, this effect */
/*       should be small, but keep it in mind. */
/* ********************************* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /input1/ */
/* c      SAVE /RNDF77/ */
    px0 = *px;
    py0 = *py;
    pz0 = *pz;
/* determine the no. of pions in the final state using a */
/* statistical model */
    pbarfs_(srt, &npion, &input1_1.iseed);
/* find the masses of the final state particles before calculate */
/* their momenta, and relable them. The masses of rho and omega */
/* will be generated according to the Breit Wigner formula       (NOTE!!! */
/* NOT DONE YET, AT THE MOMENT LET US USE FIXED RHO AND OMEGA MAEES) */
/* bali2/22/99 */
/* Here we generate two stes of integer random numbers (3,4,5) */
/* one or both of them are used directly as the lables of pions */
/* similarly, 22+nchrg1 and 22+nchrg2 are used directly */
/* to label rhos */
    nchrg1 = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
    nchrg2 = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
/* the corresponding masses of pions */
    pmass1 = .13496f;
    pmass2 = .13496f;
    if (nchrg1 == 3 || nchrg1 == 5) {
	pmass1 = .13957f;
    }
    if (nchrg2 == 3 || nchrg2 == 5) {
	pmass2 = .13957f;
    }
/* (1) for 2 pion production */
    if (npion == 2) {
	*iblock = 1902;
/* randomly generate the charges of final state particles, */
	ee_1.lb[*i1 - 1] = nchrg1;
	cc_1.e[*i1 - 1] = pmass1;
	ee_1.lb[*i2 - 1] = nchrg2;
	cc_1.e[*i2 - 1] = pmass2;
/* TO CALCULATE THE FINAL MOMENTA */
	goto L50;
    }
/* (2) FOR 3 PION PRODUCTION */
    if (npion == 3) {
	*iblock = 1903;
	ee_1.lb[*i1 - 1] = nchrg1;
	cc_1.e[*i1 - 1] = pmass1;
	ee_1.lb[*i2 - 1] = nchrg2 + 22;
	cc_1.e[*i2 - 1] = .769f;
	goto L50;
    }
/* (3) FOR 4 PION PRODUCTION */
/* we allow both rho+rho and pi+omega with 50-50% probability */
    if (npion == 4) {
	*iblock = 1904;
/* determine rho+rho or pi+omega */
	if (ranart_(&rndf77_1.nseed) >= .5f) {
/* rho+rho */
	    ee_1.lb[*i1 - 1] = nchrg1 + 22;
	    cc_1.e[*i1 - 1] = .769f;
	    ee_1.lb[*i2 - 1] = nchrg2 + 22;
	    cc_1.e[*i2 - 1] = .769f;
	} else {
/* pion+omega */
	    ee_1.lb[*i1 - 1] = nchrg1;
	    cc_1.e[*i1 - 1] = pmass1;
	    ee_1.lb[*i2 - 1] = 28;
	    cc_1.e[*i2 - 1] = .782f;
	}
	goto L50;
    }
/* (4) FOR 5 PION PRODUCTION */
    if (npion == 5) {
	*iblock = 1905;
/* RHO AND OMEGA */
	ee_1.lb[*i1 - 1] = nchrg1 + 22;
	cc_1.e[*i1 - 1] = .769f;
	ee_1.lb[*i2 - 1] = 28;
	cc_1.e[*i2 - 1] = .782f;
	goto L50;
    }
/* (5) FOR 6 PION PRODUCTION */
    if (npion == 6) {
	*iblock = 1906;
/* OMEGA AND OMEGA */
	ee_1.lb[*i1 - 1] = 28;
	cc_1.e[*i1 - 1] = .782f;
	ee_1.lb[*i2 - 1] = 28;
	cc_1.e[*i2 - 1] = .782f;
    }
/* bali2/22/99 */
L50:
    em1 = cc_1.e[*i1 - 1];
    em2 = cc_1.e[*i2 - 1];
/* ----------------------------------------------------------------------- */
/* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH */
/* ENERGY CONSERVATION */
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = em1;
/* Computing 2nd power */
    r__4 = em2;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = em1 * em2;
    pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
    if (pr2 <= 0.f) {
	pr2 = 1e-8f;
    }
    pr = sqrt(pr2) / (*srt * 2.f);
/* WE ASSUME AN ISOTROPIC ANGULAR DISTRIBUTION IN THE CMS */
    c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
    t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
/* Computing 2nd power */
    r__1 = c1;
    s1 = sqrt(1.f - r__1 * r__1);
    ct1 = cos(t1);
    st1 = sin(t1);
/* THE MOMENTUM IN THE CMS IN THE FINAL STATE */
    *pz = pr * c1;
    *px = pr * s1 * ct1;
    *py = pr * s1 * st1;
/* ROTATE IT */
    rotate_(&px0, &py0, &pz0, px, py, pz);
    return 0;
} /* crppba_ */

/* bali2/7/99end */
/* bali3/5/99 */
/* ********************************* */
/*     PURPOSE:                                                         * */
/*     assign final states for K+K- --> light mesons */

/* Subroutine */ int crkkpi_(integer *i1, integer *i2, real *xsk1, real *xsk2,
	 real *xsk3, real *xsk4, real *xsk5, real *xsk6, real *xsk7, real *
	xsk8, real *xsk9, real *xsk10, real *xsk11, real *sigk, integer *
	iblock, integer *lbp1, integer *lbp2, real *emm1, real *emm2)
{
    static real x1;
    extern doublereal ranart_(integer *);
    extern /* Subroutine */ int rhores_(integer *, integer *);


/*     QUANTITIES:                                                     * */
/*           IBLOCK   - INFORMATION about the reaction channel          * */

/*             iblock   - 1907 */
/* ********************************* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /input1/ */
/* c      SAVE /RNDF77/ */
    *iblock = 1907;
    x1 = ranart_(&rndf77_1.nseed) * *sigk;
    *xsk2 = *xsk1 + *xsk2;
    *xsk3 = *xsk2 + *xsk3;
    *xsk4 = *xsk3 + *xsk4;
    *xsk5 = *xsk4 + *xsk5;
    *xsk6 = *xsk5 + *xsk6;
    *xsk7 = *xsk6 + *xsk7;
    *xsk8 = *xsk7 + *xsk8;
    *xsk9 = *xsk8 + *xsk9;
    *xsk10 = *xsk9 + *xsk10;
    if (x1 <= *xsk1) {
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	cc_1.e[*i1 - 1] = .13957f;
	cc_1.e[*i2 - 1] = .13957f;
	goto L100;
    } else if (x1 <= *xsk2) {
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 25;
	cc_1.e[*i1 - 1] = .13957f;
	cc_1.e[*i2 - 1] = .769f;
	goto L100;
    } else if (x1 <= *xsk3) {
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	ee_1.lb[*i2 - 1] = 28;
	cc_1.e[*i1 - 1] = .13957f;
	cc_1.e[*i2 - 1] = .782f;
	goto L100;
    } else if (x1 <= *xsk4) {
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	ee_1.lb[*i2 - 1] = 0;
	cc_1.e[*i1 - 1] = .13957f;
	cc_1.e[*i2 - 1] = .5473f;
	goto L100;
    } else if (x1 <= *xsk5) {
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 25;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 25;
	cc_1.e[*i1 - 1] = .769f;
	cc_1.e[*i2 - 1] = .769f;
	goto L100;
    } else if (x1 <= *xsk6) {
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 25;
	ee_1.lb[*i2 - 1] = 28;
	cc_1.e[*i1 - 1] = .769f;
	cc_1.e[*i2 - 1] = .782f;
	goto L100;
    } else if (x1 <= *xsk7) {
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 25;
	ee_1.lb[*i2 - 1] = 0;
	cc_1.e[*i1 - 1] = .769f;
	cc_1.e[*i2 - 1] = .5473f;
	goto L100;
    } else if (x1 <= *xsk8) {
	ee_1.lb[*i1 - 1] = 28;
	ee_1.lb[*i2 - 1] = 28;
	cc_1.e[*i1 - 1] = .782f;
	cc_1.e[*i2 - 1] = .782f;
	goto L100;
    } else if (x1 <= *xsk9) {
	ee_1.lb[*i1 - 1] = 28;
	ee_1.lb[*i2 - 1] = 0;
	cc_1.e[*i1 - 1] = .782f;
	cc_1.e[*i2 - 1] = .5473f;
	goto L100;
    } else if (x1 <= *xsk10) {
	ee_1.lb[*i1 - 1] = 0;
	ee_1.lb[*i2 - 1] = 0;
	cc_1.e[*i1 - 1] = .5473f;
	cc_1.e[*i2 - 1] = .5473f;
    } else {
	*iblock = 222;
	rhores_(i1, i2);
/*     !! phi */
	ee_1.lb[*i1 - 1] = 29;
/*          return */
	cc_1.e[*i2 - 1] = 0.f;
    }
L100:
    *lbp1 = ee_1.lb[*i1 - 1];
    *lbp2 = ee_1.lb[*i2 - 1];
    *emm1 = cc_1.e[*i1 - 1];
    *emm2 = cc_1.e[*i2 - 1];
    return 0;
} /* crkkpi_ */

/* ********************************* */
/*     PURPOSE:                                                         * */
/*             DEALING WITH K+Y -> piN scattering */

/* Subroutine */ int crkhyp_(real *px, real *py, real *pz, real *srt, integer 
	*i1, integer *i2, real *xky1, real *xky2, real *xky3, real *xky4, 
	real *xky5, real *xky6, real *xky7, real *xky8, real *xky9, real *
	xky10, real *xky11, real *xky12, real *xky13, real *xky14, real *
	xky15, real *xky16, real *xky17, real *sigk, integer *ikmp, integer *
	iblock)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    double sqrt(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static real c1, s1, t1, x1, pr, em1, em2, ct1, pr2, px0, py0, pz0, st1;
    extern doublereal ranart_(integer *);
    extern /* Subroutine */ int rotate_(real *, real *, real *, real *, real *
	    , real *);


/*             Determine:                                               * */
/*             (1) relable particles in the final state                 * */
/*             (2) new momenta of final state particles                 * */
/*                                                                        * */
/*     QUANTITIES:                                                    * */
/*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME* */
/*           SRT      - SQRT OF S                                       * */
/*           IBLOCK   - INFORMATION about the reaction channel          * */
/*                                                                     * */
/*             iblock   - 1908                                          * */
/*             iblock   - 222   !! phi                                  * */
/* ********************************* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /input1/ */
/* c      SAVE /RNDF77/ */
    px0 = *px;
    py0 = *py;
    pz0 = *pz;
    *iblock = 1908;

    x1 = ranart_(&rndf77_1.nseed) * *sigk;
    *xky2 = *xky1 + *xky2;
    *xky3 = *xky2 + *xky3;
    *xky4 = *xky3 + *xky4;
    *xky5 = *xky4 + *xky5;
    *xky6 = *xky5 + *xky6;
    *xky7 = *xky6 + *xky7;
    *xky8 = *xky7 + *xky8;
    *xky9 = *xky8 + *xky9;
    *xky10 = *xky9 + *xky10;
    *xky11 = *xky10 + *xky11;
    *xky12 = *xky11 + *xky12;
    *xky13 = *xky12 + *xky13;
    *xky14 = *xky13 + *xky14;
    *xky15 = *xky14 + *xky15;
    *xky16 = *xky15 + *xky16;
    if (x1 <= *xky1) {
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	cc_1.e[*i1 - 1] = .14f;
	cc_1.e[*i2 - 1] = .93828f;
	goto L100;
    } else if (x1 <= *xky2) {
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 4) + 6;
	cc_1.e[*i1 - 1] = .14f;
	cc_1.e[*i2 - 1] = 1.232f;
	goto L100;
    } else if (x1 <= *xky3) {
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 10;
	cc_1.e[*i1 - 1] = .14f;
	cc_1.e[*i2 - 1] = 1.44f;
	goto L100;
    } else if (x1 <= *xky4) {
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 12;
	cc_1.e[*i1 - 1] = .14f;
	cc_1.e[*i2 - 1] = 1.535f;
	goto L100;
    } else if (x1 <= *xky5) {
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 25;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	cc_1.e[*i1 - 1] = .769f;
	cc_1.e[*i2 - 1] = .93828f;
	goto L100;
    } else if (x1 <= *xky6) {
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 25;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 4) + 6;
	cc_1.e[*i1 - 1] = .769f;
	cc_1.e[*i2 - 1] = 1.232f;
	goto L100;
    } else if (x1 <= *xky7) {
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 25;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 10;
	cc_1.e[*i1 - 1] = .769f;
	cc_1.e[*i2 - 1] = 1.44f;
	goto L100;
    } else if (x1 <= *xky8) {
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 25;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 12;
	cc_1.e[*i1 - 1] = .769f;
	cc_1.e[*i2 - 1] = 1.535f;
	goto L100;
    } else if (x1 <= *xky9) {
	ee_1.lb[*i1 - 1] = 28;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	cc_1.e[*i1 - 1] = .782f;
	cc_1.e[*i2 - 1] = .93828f;
	goto L100;
    } else if (x1 <= *xky10) {
	ee_1.lb[*i1 - 1] = 28;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 4) + 6;
	cc_1.e[*i1 - 1] = .782f;
	cc_1.e[*i2 - 1] = 1.232f;
	goto L100;
    } else if (x1 <= *xky11) {
	ee_1.lb[*i1 - 1] = 28;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 10;
	cc_1.e[*i1 - 1] = .782f;
	cc_1.e[*i2 - 1] = 1.44f;
	goto L100;
    } else if (x1 <= *xky12) {
	ee_1.lb[*i1 - 1] = 28;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 12;
	cc_1.e[*i1 - 1] = .782f;
	cc_1.e[*i2 - 1] = 1.535f;
	goto L100;
    } else if (x1 <= *xky13) {
	ee_1.lb[*i1 - 1] = 0;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	cc_1.e[*i1 - 1] = .5473f;
	cc_1.e[*i2 - 1] = .93828f;
	goto L100;
    } else if (x1 <= *xky14) {
	ee_1.lb[*i1 - 1] = 0;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 4) + 6;
	cc_1.e[*i1 - 1] = .5473f;
	cc_1.e[*i2 - 1] = 1.232f;
	goto L100;
    } else if (x1 <= *xky15) {
	ee_1.lb[*i1 - 1] = 0;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 10;
	cc_1.e[*i1 - 1] = .5473f;
	cc_1.e[*i2 - 1] = 1.44f;
	goto L100;
    } else if (x1 <= *xky16) {
	ee_1.lb[*i1 - 1] = 0;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 12;
	cc_1.e[*i1 - 1] = .5473f;
	cc_1.e[*i2 - 1] = 1.535f;
	goto L100;
    } else {
	ee_1.lb[*i1 - 1] = 29;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	cc_1.e[*i1 - 1] = 1.02f;
	cc_1.e[*i2 - 1] = .939457f;
	*iblock = 222;
	goto L100;
    }
L100:
    if (*ikmp == -1) {
	ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
    }
    em1 = cc_1.e[*i1 - 1];
    em2 = cc_1.e[*i2 - 1];
/* ----------------------------------------------------------------------- */
/* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH */
/* ENERGY CONSERVATION */
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = em1;
/* Computing 2nd power */
    r__4 = em2;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = em1 * em2;
    pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
    if (pr2 <= 0.f) {
	pr2 = 1e-8f;
    }
    pr = sqrt(pr2) / (*srt * 2.f);
/* WE ASSUME AN ISOTROPIC ANGULAR DISTRIBUTION IN THE CMS */
    c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
    t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
/* Computing 2nd power */
    r__1 = c1;
    s1 = sqrt(1.f - r__1 * r__1);
    ct1 = cos(t1);
    st1 = sin(t1);
/* THE MOMENTUM IN THE CMS IN THE FINAL STATE */
    *pz = pr * c1;
    *px = pr * s1 * ct1;
    *py = pr * s1 * st1;
/* ROTATE IT */
    rotate_(&px0, &py0, &pz0, px, py, pz);
    return 0;
} /* crkhyp_ */

/* ********************************* */
/*                                                                      * */
/*                                                                      * */
/* Subroutine */ int crlan_(real *px, real *py, real *pz, real *srt, integer *
	i1, integer *i2, integer *iblock)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    double sqrt(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static real c1, s1, t1, pr, em1, em2, ct1, pr2, px0, py0, pz0, st1;
    static integer ntag;
    extern doublereal ranart_(integer *);

/*     PURPOSE:                                                         * */
/*      DEALING WITH La/Si-bar + N --> K+ + pi PROCESS                  * */
/*                   La/Si + N-bar --> K- + pi                          * */
/*     NOTE   :                                                         * */

/*     QUANTITIES:                                                      * */
/*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME* */
/*           SRT      - SQRT OF S                                       * */
/*           IBLOCK   - THE INFORMATION BACK                            * */
/*                      71 */
/* ********************************* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /input1/ */
/* c      SAVE /RNDF77/ */
    px0 = *px;
    py0 = *py;
    pz0 = *pz;
    *iblock = 71;
    ntag = 0;
    if (ee_1.lb[*i1 - 1] >= 14 && ee_1.lb[*i1 - 1] <= 17 || ee_1.lb[*i2 - 1] 
	    >= 14 && ee_1.lb[*i2 - 1] <= 17) {
	ee_1.lb[*i1 - 1] = 21;
    } else {
	ee_1.lb[*i1 - 1] = 23;
    }
    ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
    cc_1.e[*i1 - 1] = .498f;
    cc_1.e[*i2 - 1] = .138f;
    em1 = cc_1.e[*i1 - 1];
    em2 = cc_1.e[*i2 - 1];
/* ----------------------------------------------------------------------- */
/* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH */
/* ENERGY CONSERVATION */
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = em1;
/* Computing 2nd power */
    r__4 = em2;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = em1 * em2;
    pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
    if (pr2 <= 0.f) {
	pr2 = 1e-9f;
    }
    pr = sqrt(pr2) / (*srt * 2.f);
    c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
    t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
/* Computing 2nd power */
    r__1 = c1;
    s1 = sqrt(1.f - r__1 * r__1);
    ct1 = cos(t1);
    st1 = sin(t1);
/* THE MOMENTUM IN THE CMS IN THE FINAL STATE */
    *pz = pr * c1;
    *px = pr * s1 * ct1;
    *py = pr * s1 * st1;
/* FOR THE ISOTROPIC DISTRIBUTION THERE IS NO NEED TO ROTATE */
    return 0;
} /* crlan_ */

/* sp11/03/01 end */
/* ********************************* */
/* ********************************* */
/*                                                                      * */
/*                                                                      * */
/* Subroutine */ int crkpla_(real *px, real *py, real *pz, real *ec, real *
	srt, real *spika, real *emm1, real *emm2, integer *lbp1, integer *
	lbp2, integer *i1, integer *i2, integer *icase, real *srhoks)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5, r__6;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static real c1, s1, t1;
    static integer ic;
    static real pr, ct1, pr2, px0, py0, pz0, st1, pdd, pff, sig1, sig2, xkp0, 
	    xkp1, xkp2, xkp3, xkp4, xkp5, xkp6, xkp7, xkp8, xkp9, sigm, dskn, 
	    xkp10, randu, sigkp, dsknr, scheck;
    extern /* Subroutine */ int distce_(integer *, integer *, real *, real *, 
	    real *, real *, real *, integer *, real *, real *, real *);
    static real sigpik;
    extern doublereal ranart_(integer *);

    /* Fortran I/O blocks */
    static cilist io___1752 = { 0, 99, 0, 0, 0 };
    static cilist io___1753 = { 0, 99, 0, 0, 0 };


/*     PURPOSE:                                                         * */
/*     DEALING WITH  K+ + Pi ---> La/Si-bar + B, phi+K, phi+K* OR  K* * */
/*                   K- + Pi ---> La/Si + B-bar  OR   K*-bar          * */
/*     NOTE   :                                                         * */

/*     QUANTITIES:                                                      * */
/*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME* */
/*           SRT      - SQRT OF S                                       * */
/*           IBLOCK   - THE INFORMATION BACK                            * */
/*                      71 */
/* ********************************* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /input1/ */
/* c      SAVE /RNDF77/ */
    *emm1 = 0.f;
    *emm2 = 0.f;
    *lbp1 = 0;
    *lbp2 = 0;
    xkp0 = *spika;
    xkp1 = 0.f;
    xkp2 = 0.f;
    xkp3 = 0.f;
    xkp4 = 0.f;
    xkp5 = 0.f;
    xkp6 = 0.f;
    xkp7 = 0.f;
    xkp8 = 0.f;
    xkp9 = 0.f;
    xkp10 = 0.f;
    sigm = 15.f;
/*         if(lb(i1).eq.21.or.lb(i2).eq.21)sigm=10. */
/* Computing 2nd power */
    r__1 = *srt;
/* Computing 2nd power */
    r__2 = *srt;
    pdd = (r__1 * r__1 - .40063836159999994f) * (r__2 * r__2 - 
	    .13179804160000003f);

    if (*srt < 2.0551569999999999f) {
	goto L70;
    }
/* Computing 2nd power */
    r__1 = *srt;
/* Computing 2nd power */
    r__2 = *srt;
    xkp1 = sigm * 1.3333333333333333f * (r__1 * r__1 - 4.2236702946489997f) * 
	    (r__2 * r__2 - .031061595048999975f) / pdd;
    if (*srt > 2.3476999999999997f) {
/* Computing 2nd power */
	r__1 = *srt;
/* Computing 2nd power */
	r__2 = *srt;
	xkp2 = sigm * 5.333333333333333f * (r__1 * r__1 - 5.5116952899999987f)
		 * (r__2 * r__2 - .013525690000000016f) / pdd;
    }
    if (*srt > 2.5556999999999999f) {
/* Computing 2nd power */
	r__1 = *srt;
/* Computing 2nd power */
	r__2 = *srt;
	xkp3 = sigm * 1.3333333333333333f * (r__1 * r__1 - 
		6.5316024899999992f) * (r__2 * r__2 - .10517049000000002f) / 
		pdd;
    }
    if (*srt > 2.6506999999999996f) {
/* Computing 2nd power */
	r__1 = *srt;
/* Computing 2nd power */
	r__2 = *srt;
	xkp4 = sigm * 1.3333333333333333f * (r__1 * r__1 - 
		7.0262104899999978f) * (r__2 * r__2 - .17581249000000002f) / 
		pdd;
    }

    if (*srt > 2.136857f) {
/* Computing 2nd power */
	r__1 = *srt;
/* Computing 2nd power */
	r__2 = *srt;
	xkp5 = sigm * 4.f * (r__1 * r__1 - 4.5661578384489996f) * (r__2 * 
		r__2 - .066534591249000019f) / pdd;
    }
    if (*srt > 2.4294000000000002f) {
/* Computing 2nd power */
	r__1 = *srt;
/* Computing 2nd power */
	r__2 = *srt;
	xkp6 = sigm * 16.f * (r__1 * r__1 - 5.901984360000001f) * (r__2 * 
		r__2 - .0011971599999999975f) / pdd;
    }
    if (*srt > 2.6374f) {
/* Computing 2nd power */
	r__1 = *srt;
/* Computing 2nd power */
	r__2 = *srt;
	xkp7 = sigm * 4.f * (r__1 * r__1 - 6.95587876f) * (r__2 * r__2 - 
		.058854759999999964f) / pdd;
    }
    if (*srt > 2.7324000000000002f) {
/* Computing 2nd power */
	r__1 = *srt;
/* Computing 2nd power */
	r__2 = *srt;
	xkp8 = sigm * 4.f * (r__1 * r__1 - 7.4660097600000013f) * (r__2 * 
		r__2 - .11397375999999994f) / pdd;
    }
L70:
    sig1 = 195.639f;
    sig2 = 372.378f;
    if (*srt > 1.518f) {
/* Computing 2nd power */
	r__1 = *srt;
/* Computing 2nd power */
	r__2 = *srt;
	pff = sqrt((r__1 * r__1 - 2.3043240000000003f) * (r__2 * r__2 - 
		.272484f));
/* lin-9/2012: check argument in sqrt(): */
	scheck = pdd;
	if (scheck <= 0.f) {
	    s_wsle(&io___1752);
	    do_lio(&c__9, &c__1, "scheck40: ", (ftnlen)10);
	    do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	    e_wsle();
	    s_stop("", (ftnlen)0);
	}
/* Computing 2nd power */
	r__1 = *srt;
	xkp9 = sig1 * pff / sqrt(pdd) * 1.f / 32.f / 3.1415926f / (r__1 * 
		r__1);
	if (*srt > 1.915f) {
/* Computing 2nd power */
	    r__1 = *srt;
/* Computing 2nd power */
	    r__2 = *srt;
	    pff = sqrt((r__1 * r__1 - 3.6672250000000002f) * (r__2 * r__2 - 
		    .015625f));
/* lin-9/2012: check argument in sqrt(): */
	    scheck = pdd;
	    if (scheck <= 0.f) {
		s_wsle(&io___1753);
		do_lio(&c__9, &c__1, "scheck41: ", (ftnlen)10);
		do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
		e_wsle();
		s_stop("", (ftnlen)0);
	    }
/* Computing 2nd power */
	    r__1 = *srt;
	    xkp10 = sig2 * pff / sqrt(pdd) * 3.f / 32.f / 3.1415926f / (r__1 *
		     r__1);
	}
    }
/* lin-8/15/02 K pi -> K* (rho omega), from detailed balance, */
/* neglect rho and omega mass difference for now: */
    sigpik = 0.f;
    if (*srt > 1.6640000000000001f) {
/* Computing 2nd power */
	r__1 = *srt;
/* Computing 2nd power */
	r__2 = *srt;
/* Computing 2nd power */
	r__3 = *srt;
/* Computing 2nd power */
	r__4 = *px;
/* Computing 2nd power */
	r__5 = *py;
/* Computing 2nd power */
	r__6 = *pz;
	sigpik = *srhoks * 9.f * (r__1 * r__1 - .015625f) * (r__2 * r__2 - 
		2.7722250000000002f) / 4 / (r__3 * r__3) / (r__4 * r__4 + 
		r__5 * r__5 + r__6 * r__6);
	if (*srt > 1.677f) {
	    sigpik = sigpik * 12.f / 9.f;
	}
    }

    sigkp = xkp0 + xkp1 + xkp2 + xkp3 + xkp4 + xkp5 + xkp6 + xkp7 + xkp8 + 
	    xkp9 + xkp10 + sigpik;
    *icase = 0;
    dskn = sqrt(sigkp / 3.1415926f / 10.f);
    dsknr = dskn + .1f;
    distce_(i1, i2, &dsknr, &dskn, &input1_1.dt, ec, srt, &ic, px, py, pz);
    if (ic == -1) {
	return 0;
    }

    randu = ranart_(&rndf77_1.nseed) * sigkp;
    xkp1 = xkp0 + xkp1;
    xkp2 = xkp1 + xkp2;
    xkp3 = xkp2 + xkp3;
    xkp4 = xkp3 + xkp4;
    xkp5 = xkp4 + xkp5;
    xkp6 = xkp5 + xkp6;
    xkp7 = xkp6 + xkp7;
    xkp8 = xkp7 + xkp8;
    xkp9 = xkp8 + xkp9;
    xkp10 = xkp9 + xkp10;

/*   !! K* formation */
    if (randu <= xkp0) {
	*icase = 1;
	return 0;
    } else {
/* La/Si-bar + B formation */
	*icase = 2;
	if (randu <= xkp1) {
	    *lbp1 = -14;
	    *lbp2 = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	    *emm1 = 1.1157f;
	    *emm2 = .939457f;
	    goto L60;
	} else if (randu <= xkp2) {
	    *lbp1 = -14;
	    *lbp2 = (integer) (ranart_(&rndf77_1.nseed) * 4) + 6;
	    *emm1 = 1.1157f;
	    *emm2 = 1.232f;
	    goto L60;
	} else if (randu <= xkp3) {
	    *lbp1 = -14;
	    *lbp2 = (integer) (ranart_(&rndf77_1.nseed) * 2) + 10;
	    *emm1 = 1.1157f;
	    *emm2 = 1.44f;
	    goto L60;
	} else if (randu <= xkp4) {
	    *lbp1 = -14;
	    *lbp2 = (integer) (ranart_(&rndf77_1.nseed) * 2) + 12;
	    *emm1 = 1.1157f;
	    *emm2 = 1.535f;
	    goto L60;
	} else if (randu <= xkp5) {
	    *lbp1 = -15 - (integer) (ranart_(&rndf77_1.nseed) * 3);
	    *lbp2 = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	    *emm1 = 1.1974f;
	    *emm2 = .939457f;
	    goto L60;
	} else if (randu <= xkp6) {
	    *lbp1 = -15 - (integer) (ranart_(&rndf77_1.nseed) * 3);
	    *lbp2 = (integer) (ranart_(&rndf77_1.nseed) * 4) + 6;
	    *emm1 = 1.1974f;
	    *emm2 = 1.232f;
	    goto L60;
	} else if (randu < xkp7) {
	    *lbp1 = -15 - (integer) (ranart_(&rndf77_1.nseed) * 3);
	    *lbp2 = (integer) (ranart_(&rndf77_1.nseed) * 2) + 10;
	    *emm1 = 1.1974f;
	    *emm2 = 1.44f;
	    goto L60;
	} else if (randu < xkp8) {
	    *lbp1 = -15 - (integer) (ranart_(&rndf77_1.nseed) * 3);
	    *lbp2 = (integer) (ranart_(&rndf77_1.nseed) * 2) + 12;
	    *emm1 = 1.1974f;
	    *emm2 = 1.535f;
	    goto L60;
	} else if (randu < xkp9) {
/*       !! phi +K  formation (iblock=224) */
	    *icase = 3;
	    *lbp1 = 29;
	    *lbp2 = 23;
	    *emm1 = 1.02f;
	    *emm2 = .498f;
	    if (ee_1.lb[*i1 - 1] == 21 || ee_1.lb[*i2 - 1] == 21) {
/*         !! phi +K-bar  formation (iblock=124) */
		*lbp2 = 21;
		*icase = -3;
	    }
	    goto L60;
	} else if (randu < xkp10) {
/*       !! phi +K* formation (iblock=226) */
	    *icase = 4;
	    *lbp1 = 29;
	    *lbp2 = 30;
	    *emm1 = 1.02f;
	    *emm2 = .895f;
	    if (ee_1.lb[*i1 - 1] == 21 || ee_1.lb[*i2 - 1] == 21) {
		*lbp2 = -30;
		*icase = -4;
	    }
	    goto L60;
	} else {
/*       !! (rho,omega) +K* formation (iblock=88) */
	    *icase = 5;
	    *lbp1 = (integer) (ranart_(&rndf77_1.nseed) * 3) + 25;
	    *lbp2 = 30;
	    *emm1 = .769f;
	    *emm2 = .895f;
	    if (*srt > 1.677f && ranart_(&rndf77_1.nseed) < .25f) {
		*lbp1 = 28;
		*emm1 = .782f;
	    }
	    if (ee_1.lb[*i1 - 1] == 21 || ee_1.lb[*i2 - 1] == 21) {
		*lbp2 = -30;
		*icase = -5;
	    }
	}
    }

L60:
    if (*icase == 2 && (ee_1.lb[*i1 - 1] == 21 || ee_1.lb[*i2 - 1] == 21)) {
	*lbp1 = -(*lbp1);
	*lbp2 = -(*lbp2);
    }
    px0 = *px;
    py0 = *py;
    pz0 = *pz;
/* ----------------------------------------------------------------------- */
/* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH */
/* ENERGY CONSERVATION */
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = *emm1;
/* Computing 2nd power */
    r__4 = *emm2;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = *emm1 * *emm2;
    pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
    if (pr2 <= 0.f) {
	pr2 = 1e-9f;
    }
    pr = sqrt(pr2) / (*srt * 2.f);
    c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
    t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
/* Computing 2nd power */
    r__1 = c1;
    s1 = sqrt(1.f - r__1 * r__1);
    ct1 = cos(t1);
    st1 = sin(t1);
/* THE MOMENTUM IN THE CMS IN THE FINAL STATE */
    *pz = pr * c1;
    *px = pr * s1 * ct1;
    *py = pr * s1 * st1;
/* FOR THE ISOTROPIC DISTRIBUTION THERE IS NO NEED TO ROTATE */
    return 0;
} /* crkpla_ */

/* ********************************* */
/*                                                                      * */
/*                                                                      * */
/* Subroutine */ int crkphi_(real *px, real *py, real *pz, real *ec, real *
	srt, integer *iblock, real *emm1, real *emm2, integer *lbp1, integer *
	lbp2, integer *i1, integer *i2, integer *ikk, integer *icase, real *
	rrkk, real *prkk)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal), cos(
	    doublereal), sin(doublereal);

    /* Local variables */
    static real c1, s1, t1;
    static integer ic;
    static real pr;
    static integer lb1, lb2;
    static real ct1, pr2, px0, py0, pz0, st1, pii, sig, dnr, sig1, sig2, sig3,
	     xsk1, srr1, srr2, srr3, xsk2, xsk3, xsk4, xsk5, xsk6, xsk7, xsk8,
	     xsk9, sigm, dskn, xsk10, xsk11, ranx, srri, srrt, sigm0, prkk0, 
	    rrkk0, sigks, dsknr, sigks1, sigks2, sigks3, sigks4;
    extern /* Subroutine */ int distce_(integer *, integer *, real *, real *, 
	    real *, real *, real *, integer *, real *, real *, real *), 
	    crkkpi_(integer *, integer *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    integer *, integer *, integer *, real *, real *);
    extern doublereal ranart_(integer *);
    extern /* Subroutine */ int xkkann_(real *, real *, real *, real *, real *
	    , real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *), crkspi_(integer *, integer *, real *, real *, real *, 
	    real *, real *, integer *, integer *, integer *, real *, real *), 
	    xkksan_(integer *, integer *, real *, real *, real *, real *, 
	    real *, real *, real *);

/*     PURPOSE:                                                         * */
/*     DEALING WITH   KKbar, KK*bar, KbarK*, K*K*bar --> Phi + pi(rho,omega) */
/*     and KKbar --> (pi eta) (pi eta), (rho omega) (rho omega) */
/*     and KK*bar or Kbar K* --> (pi eta) (rho omega) */

/*     NOTE   :                                                         * */

/*     QUANTITIES:                                                      * */
/*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME* */
/*           SRT      - SQRT OF S                                       * */
/*           IBLOCK   - THE INFORMATION BACK                            * */
/*                      222 */
/* ********************************* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /input1/ */
/* c      SAVE /RNDF77/ */
    lb1 = ee_1.lb[*i1 - 1];
    lb2 = ee_1.lb[*i2 - 1];
    *icase = 0;
/*        if(srt .lt. aphi+ap1)return */
/* c        if(srt .lt. aphi+ap1) then */
    if (*srt < 1.15496f) {
	sig1 = 0.f;
	sig2 = 0.f;
	sig3 = 0.f;
    } else {

	if (lb1 == 23 && lb2 == 21 || lb2 == 23 && lb1 == 21) {
	    dnr = 4.f;
	    *ikk = 2;
	} else if (lb1 == 21 && lb2 == 30 || lb2 == 21 && lb1 == 30 || lb1 == 
		23 && lb2 == -30 || lb2 == 23 && lb1 == -30) {
	    dnr = 12.f;
	    *ikk = 1;
	} else {
	    dnr = 36.f;
	    *ikk = 0;
	}
	sig1 = 0.f;
	sig2 = 0.f;
	sig3 = 0.f;
	srri = cc_1.e[*i1 - 1] + cc_1.e[*i2 - 1];
	srr1 = 1.15496f;
	srr2 = 1.8019000000000001f;
	srr3 = 1.79f;

/* Computing 2nd power */
	r__1 = *srt;
/* Computing 2nd power */
	r__2 = cc_1.e[*i1 - 1] + cc_1.e[*i2 - 1];
/* Computing 2nd power */
	r__3 = *srt;
/* Computing 2nd power */
	r__4 = cc_1.e[*i1 - 1] - cc_1.e[*i2 - 1];
	pii = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4);
	srrt = *srt - dmax(srri,srr1);
/* c   to avoid divergent/negative values at small srrt: */
/*          if(srrt .lt. 0.3)then */
	if (srrt < .3f && srrt > .01f) {
	    d__1 = (doublereal) srrt;
	    sig = 1.69f / (pow_dd(&d__1, &c_b783) - .407f);
	} else {
	    d__1 = (doublereal) srrt;
	    sig = pow_dd(&d__1, &c_b784) * .008f + 3.74f;
	}
/* Computing 2nd power */
	r__1 = *srt;
/* Computing 2nd power */
	r__2 = *srt;
	sig1 = sig * (9.f / dnr) * (r__1 * r__1 - 1.3339326015999999f) * (
		r__2 * r__2 - .78329580160000012f) / pii;
	if (*srt > 1.8019000000000001f) {
	    srrt = *srt - dmax(srri,srr2);
/* c         if(srrt .lt. 0.3)then */
	    if (srrt < .3f && srrt > .01f) {
		d__1 = (doublereal) srrt;
		sig = 1.69f / (pow_dd(&d__1, &c_b783) - .407f);
	    } else {
		d__1 = (doublereal) srrt;
		sig = pow_dd(&d__1, &c_b784) * .008f + 3.74f;
	    }
/* Computing 2nd power */
	    r__1 = *srt;
/* Computing 2nd power */
	    r__2 = *srt;
	    sig2 = sig * (9.f / dnr) * (r__1 * r__1 - 3.24684361f) * (r__2 * 
		    r__2 - .05669160999999999f) / pii;
	}
	if (*srt > 1.79f) {
	    srrt = *srt - dmax(srri,srr3);
/* c         if(srrt .lt. 0.3)then */
	    if (srrt < .3f && srrt > .01f) {
		d__1 = (doublereal) srrt;
		sig = 1.69f / (pow_dd(&d__1, &c_b783) - .407f);
	    } else {
		d__1 = (doublereal) srrt;
		sig = pow_dd(&d__1, &c_b784) * .008f + 3.74f;
	    }
/* Computing 2nd power */
	    r__1 = *srt;
/* Computing 2nd power */
	    r__2 = *srt;
	    sig3 = sig * (27.f / dnr) * (r__1 * r__1 - 3.2040999999999999f) * 
		    (r__2 * r__2 - .0625f) / pii;
	}
/*         sig1 = amin1(20.,sig1) */
/*         sig2 = amin1(20.,sig2) */
/*         sig3 = amin1(20.,sig3) */
    }
    rrkk0 = *rrkk;
    prkk0 = *prkk;
    sigm = 0.f;
    if (lb1 == 23 && lb2 == 21 || lb2 == 23 && lb1 == 21) {
	xkkann_(srt, &xsk1, &xsk2, &xsk3, &xsk4, &xsk5, &xsk6, &xsk7, &xsk8, &
		xsk9, &xsk10, &xsk11, &sigm, &rrkk0);
    } else if (lb1 == 21 && lb2 == 30 || lb2 == 21 && lb1 == 30 || lb1 == 23 
	    && lb2 == -30 || lb2 == 23 && lb1 == -30) {
	xkksan_(i1, i2, srt, &sigks1, &sigks2, &sigks3, &sigks4, &sigm, &
		prkk0);
    } else {
    }

/*         sigks = sig1 + sig2 + sig3 */
    sigm0 = sigm;
    sigks = sig1 + sig2 + sig3 + sigm;
    dskn = sqrt(sigks / 3.1415926f / 10.f);
    dsknr = dskn + .1f;
    distce_(i1, i2, &dsknr, &dskn, &input1_1.dt, ec, srt, &ic, px, py, pz);
    if (ic == -1) {
	return 0;
    }
    *icase = 1;
    ranx = ranart_(&rndf77_1.nseed);
    *lbp1 = 29;
    *emm1 = 1.02f;
    if (ranx <= sig1 / sigks) {
	*lbp2 = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	*emm2 = .13496f;
    } else if (ranx <= (sig1 + sig2) / sigks) {
	*lbp2 = 28;
	*emm2 = .7819f;
    } else if (ranx <= (sig1 + sig2 + sig3) / sigks) {
	*lbp2 = (integer) (ranart_(&rndf77_1.nseed) * 3) + 25;
	*emm2 = .77f;
    } else {
	if (lb1 == 23 && lb2 == 21 || lb2 == 23 && lb1 == 21) {
	    crkkpi_(i1, i2, &xsk1, &xsk2, &xsk3, &xsk4, &xsk5, &xsk6, &xsk7, &
		    xsk8, &xsk9, &xsk10, &xsk11, &sigm0, iblock, lbp1, lbp2, 
		    emm1, emm2);
	} else if (lb1 == 21 && lb2 == 30 || lb2 == 21 && lb1 == 30 || lb1 == 
		23 && lb2 == -30 || lb2 == 23 && lb1 == -30) {
	    crkspi_(i1, i2, &sigks1, &sigks2, &sigks3, &sigks4, &sigm0, 
		    iblock, lbp1, lbp2, emm1, emm2);
	} else {
	}
    }

    px0 = *px;
    py0 = *py;
    pz0 = *pz;
/* ----------------------------------------------------------------------- */
/* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH */
/* ENERGY CONSERVATION */
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = *emm1;
/* Computing 2nd power */
    r__4 = *emm2;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = *emm1 * *emm2;
    pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
    if (pr2 <= 0.f) {
	pr2 = 1e-9f;
    }
    pr = sqrt(pr2) / (*srt * 2.f);
    c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
    t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
/* Computing 2nd power */
    r__1 = c1;
    s1 = sqrt(1.f - r__1 * r__1);
    ct1 = cos(t1);
    st1 = sin(t1);
/* THE MOMENTUM IN THE CMS IN THE FINAL STATE */
    *pz = pr * c1;
    *px = pr * s1 * ct1;
    *py = pr * s1 * st1;
/* FOR THE ISOTROPIC DISTRIBUTION THERE IS NO NEED TO ROTATE */
    return 0;
} /* crkphi_ */

/* sp11/21/01 end */
/* ********************************* */
/*                                                                      * */
/*                                                                      * */
/* Subroutine */ int crksph_(real *px, real *py, real *pz, real *ec, real *
	srt, real *emm1, real *emm2, integer *lbp1, integer *lbp2, integer *
	i1, integer *i2, integer *ikkg, integer *ikkl, integer *iblock, 
	integer *icase, real *srhoks)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    double sqrt(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static real c1, s1, t1;
    static integer ic;
    static real pr;
    static integer lb1, lb2;
    static real ct1, pr2, px0, py0, pz0, st1, pff, pii, dnr, sig1, sig2, 
	    sig11, sig22, dskn, ranx, sigkm, sigks, dsknr, scheck, sigela;
    extern /* Subroutine */ int distce_(integer *, integer *, real *, real *, 
	    real *, real *, real *, integer *, real *, real *, real *);
    extern doublereal ranart_(integer *);

    /* Fortran I/O blocks */
    static cilist io___1827 = { 0, 99, 0, 0, 0 };
    static cilist io___1829 = { 0, 99, 0, 0, 0 };


/*     PURPOSE:                                                         * */
/*     DEALING WITH   K + rho(omega) or K* + pi(rho,omega) */
/*                    --> Phi + K(K*), pi + K* or pi + K, and elastic */
/*     NOTE   :                                                         * */

/*     QUANTITIES:                                                      * */
/*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME* */
/*           SRT      - SQRT OF S                                       * */
/*           IBLOCK   - THE INFORMATION BACK                            * */
/*                      222 */
/*                      223 --> phi + pi(rho,omega) */
/*                      224 --> phi + K <-> K + pi(rho,omega) */
/*                      225 --> phi + K <-> K* + pi(rho,omega) */
/*                      226 --> phi + K* <-> K + pi(rho,omega) */
/*                      227 --> phi + K* <-> K* + pi(rho,omega) */
/* ********************************* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /input1/ */
/* c      SAVE /RNDF77/ */
    lb1 = ee_1.lb[*i1 - 1];
    lb2 = ee_1.lb[*i2 - 1];
    *icase = 0;
    sigela = 10.f;
    sigkm = 0.f;
/*     K(K*) + rho(omega) -> pi K*(K) */
    if (lb1 >= 25 && lb1 <= 28 || lb2 >= 25 && lb2 <= 28) {
	if (abs(lb1) == 30 || abs(lb2) == 30) {
	    sigkm = *srhoks;
/* lin-2/26/03 check whether (rho K) is above the (pi K*) thresh: */
	} else if ((lb1 == 23 || lb1 == 21 || lb2 == 23 || lb2 == 21) && *srt 
		> 1.03457f) {
	    sigkm = *srhoks;
	}
    }
/*        if(srt .lt. aphi+aka)return */
    if (*srt < 1.518f) {
	sig11 = 0.f;
	sig22 = 0.f;
    } else {
/* K*-bar +pi --> phi + (K,K*)-bar */
	if (abs(lb1) == 30 && (lb2 >= 3 && lb2 <= 5) || abs(lb2) == 30 && (
		lb1 >= 3 && lb1 <= 5)) {
	    dnr = 18.f;
	    *ikkl = 0;
	    *iblock = 225;
/*               sig1 = 15.0 */
/*               sig2 = 30.0 */
/* lin-2/06/03 these large values reduces to ~10 mb for sig11 or sig22 */
/*     due to the factors of ~1/(32*pi*s)~1/200: */
	    sig1 = 2047.042f;
	    sig2 = 1496.692f;
/* K(-bar)+rho --> phi + (K,K*)-bar */
	} else if (lb1 == 23 || lb1 == 21 && (lb2 >= 25 && lb2 <= 27) || (lb2 
		== 23 || lb2 == 21 && (lb1 >= 25 && lb1 <= 27))) {
	    dnr = 18.f;
	    *ikkl = 1;
	    *iblock = 224;
/*               sig1 = 3.5 */
/*               sig2 = 9.0 */
	    sig1 = 526.702f;
	    sig2 = 1313.96f;
/* K*(-bar) +rho */
	} else if (abs(lb1) == 30 && (lb2 >= 25 && lb2 <= 27) || abs(lb2) == 
		30 && (lb1 >= 25 && lb1 <= 27)) {
	    dnr = 54.f;
	    *ikkl = 0;
	    *iblock = 225;
/*               sig1 = 3.5 */
/*               sig2 = 9.0 */
	    sig1 = 1371.257f;
	    sig2 = 6999.84f;
/* K(-bar) + omega */
	} else if ((lb1 == 23 || lb1 == 21) && lb2 == 28 || (lb2 == 23 || lb2 
		== 21) && lb1 == 28) {
	    dnr = 6.f;
	    *ikkl = 1;
	    *iblock = 224;
/*               sig1 = 3.5 */
/*               sig2 = 6.5 */
	    sig1 = 355.429f;
	    sig2 = 440.558f;
/* K*(-bar) +omega */
	} else {
	    dnr = 18.f;
	    *ikkl = 0;
	    *iblock = 225;
/*               sig1 = 3.5 */
/*               sig2 = 15.0 */
	    sig1 = 482.292f;
	    sig2 = 1698.903f;
	}
	sig11 = 0.f;
	sig22 = 0.f;
/*         sig11=sig1*(6./dnr)*(srt**2-(aphi+aka)**2)* */
/*    &           (srt**2-(aphi-aka)**2)/(srt**2-(e(i1)+e(i2))**2)/ */
/*    &           (srt**2-(e(i1)-e(i2))**2) */
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
	r__1 = *srt;
/* Computing 2nd power */
	r__2 = cc_1.e[*i1 - 1] + cc_1.e[*i2 - 1];
/* Computing 2nd power */
	r__3 = *srt;
/* Computing 2nd power */
	r__4 = cc_1.e[*i1 - 1] - cc_1.e[*i2 - 1];
	scheck = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4);
	if (scheck <= 0.f) {
	    s_wsle(&io___1827);
	    do_lio(&c__9, &c__1, "scheck42: ", (ftnlen)10);
	    do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	    e_wsle();
	    s_stop("", (ftnlen)0);
	}
	pii = sqrt(scheck);
/*        pii = sqrt((srt**2-(e(i1)+e(i2))**2)*(srt**2-(e(i1)-e(i2))**2)) */
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
	r__1 = *srt;
/* Computing 2nd power */
	r__2 = *srt;
	scheck = (r__1 * r__1 - 2.3043240000000003f) * (r__2 * r__2 - 
		.272484f);
	if (scheck < 0.f) {
	    s_wsle(&io___1829);
	    do_lio(&c__9, &c__1, "scheck43: ", (ftnlen)10);
	    do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	    e_wsle();
	    scheck = 0.f;
	}
	pff = sqrt(scheck);
/*        pff = sqrt((srt**2-(aphi+aka)**2)*(srt**2-(aphi-aka)**2)) */
/* Computing 2nd power */
	r__1 = *srt;
	sig11 = sig1 * pff / pii * 6.f / dnr / 32.f / 3.1415926f / (r__1 * 
		r__1);

	if (*srt > 1.915f) {
/*         sig22=sig2*(18./dnr)*(srt**2-(aphi+aks)**2)* */
/*    &           (srt**2-(aphi-aks)**2)/(srt**2-(e(i1)+e(i2))**2)/ */
/*    &           (srt**2-(e(i1)-e(i2))**2) */
/* Computing 2nd power */
	    r__1 = *srt;
/* Computing 2nd power */
	    r__2 = *srt;
	    pff = sqrt((r__1 * r__1 - 3.6672250000000002f) * (r__2 * r__2 - 
		    .015625f));
/* Computing 2nd power */
	    r__1 = *srt;
	    sig22 = sig2 * pff / pii * 18.f / dnr / 32.f / 3.1415926f / (r__1 
		    * r__1);
	}
/*         sig11 = amin1(20.,sig11) */
/*         sig22 = amin1(20.,sig22) */

    }
/*         sigks = sig11 + sig22 */
    sigks = sig11 + sig22 + sigela + sigkm;

    dskn = sqrt(sigks / 3.1415926f / 10.f);
    dsknr = dskn + .1f;
    distce_(i1, i2, &dsknr, &dskn, &input1_1.dt, ec, srt, &ic, px, py, pz);
    if (ic == -1) {
	return 0;
    }
    *icase = 1;
    ranx = ranart_(&rndf77_1.nseed);
    if (ranx <= sigela / sigks) {
	*lbp1 = lb1;
	*emm1 = cc_1.e[*i1 - 1];
	*lbp2 = lb2;
	*emm2 = cc_1.e[*i2 - 1];
	*iblock = 111;
    } else if (ranx <= (sigela + sigkm) / sigks) {
	*lbp1 = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	*emm1 = .14f;
	if (lb1 == 23 || lb2 == 23) {
	    *lbp2 = 30;
	    *emm2 = .895f;
	} else if (lb1 == 21 || lb2 == 21) {
	    *lbp2 = -30;
	    *emm2 = .895f;
	} else if (lb1 == 30 || lb2 == 30) {
	    *lbp2 = 23;
	    *emm2 = .498f;
	} else {
	    *lbp2 = 21;
	    *emm2 = .498f;
	}
	*iblock = 112;
    } else if (ranx <= (sigela + sigkm + sig11) / sigks) {
	*lbp2 = 23;
	*emm2 = .498f;
	*ikkg = 1;
	if (lb1 == 21 || lb2 == 21 || lb1 == -30 || lb2 == -30) {
	    *lbp2 = 21;
	    *iblock += -100;
	}
	*lbp1 = 29;
	*emm1 = 1.02f;
    } else {
	*lbp2 = 30;
	*emm2 = .895f;
	*ikkg = 0;
	*iblock += 2;
	if (lb1 == 21 || lb2 == 21 || lb1 == -30 || lb2 == -30) {
	    *lbp2 = -30;
	    *iblock += -100;
	}
	*lbp1 = 29;
	*emm1 = 1.02f;
    }

    px0 = *px;
    py0 = *py;
    pz0 = *pz;
/* ----------------------------------------------------------------------- */
/* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH */
/* ENERGY CONSERVATION */
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = *emm1;
/* Computing 2nd power */
    r__4 = *emm2;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = *emm1 * *emm2;
    pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
    if (pr2 <= 0.f) {
	pr2 = 1e-9f;
    }
    pr = sqrt(pr2) / (*srt * 2.f);
    c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
    t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
/* Computing 2nd power */
    r__1 = c1;
    s1 = sqrt(1.f - r__1 * r__1);
    ct1 = cos(t1);
    st1 = sin(t1);
/* THE MOMENTUM IN THE CMS IN THE FINAL STATE */
    *pz = pr * c1;
    *px = pr * s1 * ct1;
    *py = pr * s1 * st1;
/* FOR THE ISOTROPIC DISTRIBUTION THERE IS NO NEED TO ROTATE */
    return 0;
} /* crksph_ */

/* sp11/21/01 end */
/* ********************************* */
/* ********************************* */
/* Subroutine */ int bbkaon_(integer *ic, real *srt, real *px, real *py, real 
	*pz, real *ana, real *plx, real *ply, real *plz, real *ala, real *pkx,
	 real *pky, real *pkz, integer *icou1)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    double sqrt(doublereal), cos(doublereal), sin(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static real t1, t2, ga, ek, el, en, cs, pi, bx, pk, by, bz, pn, ss, dm1, 
	    pn2, aka, fai, eln, csn, ssn, fain, elnc, dmax__, prob, pmax;
    static integer ntry;
    static real pbeta;
    extern doublereal fkaon_(real *, real *), rmass_(real *, integer *);
    static real trans0, scheck;
    extern doublereal ranart_(integer *);

    /* Fortran I/O blocks */
    static cilist io___1865 = { 0, 99, 0, 0, 0 };


/* purpose: generate the momenta for kaon,lambda/sigma and nucleon/delta */
/*          in the BB-->nlk process */
/* date: Sept. 9, 1994 */

/* c      SAVE /input1/ */
/* c      SAVE /RNDF77/ */
    pi = 3.1415962f;
    *icou1 = 0;
    aka = .498f;
    *ala = 1.116f;
    if (*ic == 2 || *ic == 4) {
	*ala = 1.197f;
    }
    *ana = .939f;
/* generate the mass of the delta */
    if (*ic > 2) {
	dmax__ = *srt - aka - *ala - .02f;
	dm1 = rmass_(&dmax__, &input1_1.iseed);
	*ana = dm1;
    }
    t1 = aka + *ana + *ala;
    t2 = *ana + *ala - aka;
    if (*srt <= t1) {
	*icou1 = -1;
	return 0;
    }
/* Computing 2nd power */
    r__1 = *srt;
/* Computing 2nd power */
    r__2 = t1;
/* Computing 2nd power */
    r__3 = *srt;
/* Computing 2nd power */
    r__4 = t2;
    pmax = sqrt((r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4)) / (
	    *srt * 2.f);
    if (pmax == 0.f) {
	pmax = 1e-9f;
    }
/* (1) Generate the momentum of the kaon according to the distribution Fkaon */
/*     and assume that the angular distribution is isotropic */
/*     in the cms of the colliding pair */
    ntry = 0;
L1:
    pk = pmax * ranart_(&rndf77_1.nseed);
    ++ntry;
    prob = fkaon_(&pk, &pmax);
    if (prob < ranart_(&rndf77_1.nseed) && ntry <= 40) {
	goto L1;
    }
    cs = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
/* Computing 2nd power */
    r__1 = cs;
    ss = sqrt(1.f - r__1 * r__1);
    fai = ranart_(&rndf77_1.nseed) * 6.2800000000000002f;
    *pkx = pk * ss * cos(fai);
    *pky = pk * ss * sin(fai);
    *pkz = pk * cs;
/* the energy of the kaon */
/* Computing 2nd power */
    r__1 = aka;
/* Computing 2nd power */
    r__2 = pk;
    ek = sqrt(r__1 * r__1 + r__2 * r__2);
/* (2) Generate the momentum of the nucleon/delta in the cms of N/delta */
/*     and lamda/sigma */
/*  the energy of the cms of NL */
    eln = *srt - ek;
    if (eln <= 0.f) {
	*icou1 = -1;
	return 0;
    }
/* beta and gamma of the cms of L/S+N */
    bx = -(*pkx) / eln;
    by = -(*pky) / eln;
    bz = -(*pkz) / eln;
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__1 = bx;
/* Computing 2nd power */
    r__2 = by;
/* Computing 2nd power */
    r__3 = bz;
    scheck = 1.f - r__1 * r__1 - r__2 * r__2 - r__3 * r__3;
    if (scheck <= 0.f) {
	s_wsle(&io___1865);
	do_lio(&c__9, &c__1, "scheck44: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	s_stop("", (ftnlen)0);
    }
    ga = 1.f / sqrt(scheck);
/*       ga=1./sqrt(1.-bx**2-by**2-bz**2) */
    elnc = eln / ga;
/* Computing 2nd power */
    r__2 = elnc;
/* Computing 2nd power */
    r__3 = *ana;
/* Computing 2nd power */
    r__4 = *ala;
/* Computing 2nd power */
    r__1 = (r__2 * r__2 + r__3 * r__3 - r__4 * r__4) / (elnc * 2.f);
/* Computing 2nd power */
    r__5 = *ana;
    pn2 = r__1 * r__1 - r__5 * r__5;
    if (pn2 <= 0.f) {
	pn2 = 1e-9f;
    }
    pn = sqrt(pn2);
    csn = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
/* Computing 2nd power */
    r__1 = csn;
    ssn = sqrt(1.f - r__1 * r__1);
    fain = ranart_(&rndf77_1.nseed) * 6.2800000000000002f;
    *px = pn * ssn * cos(fain);
    *py = pn * ssn * sin(fain);
    *pz = pn * csn;
/* Computing 2nd power */
    r__1 = *ana;
    en = sqrt(r__1 * r__1 + pn2);
/* the momentum of the lambda/sigma in the n-l cms frame is */
    *plx = -(*px);
    *ply = -(*py);
    *plz = -(*pz);
/* (3) LORENTZ-TRANSFORMATION INTO nn cms FRAME for the neutron/delta */
    pbeta = *px * bx + *py * by + *pz * bz;
    trans0 = ga * (ga * pbeta / (ga + 1.f) + en);
    *px = bx * trans0 + *px;
    *py = by * trans0 + *py;
    *pz = bz * trans0 + *pz;
/* (4) Lorentz-transformation for the lambda/sigma */
/* Computing 2nd power */
    r__1 = *ala;
/* Computing 2nd power */
    r__2 = *plx;
/* Computing 2nd power */
    r__3 = *ply;
/* Computing 2nd power */
    r__4 = *plz;
    el = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    pbeta = *plx * bx + *ply * by + *plz * bz;
    trans0 = ga * (ga * pbeta / (ga + 1.f) + el);
    *plx = bx * trans0 + *plx;
    *ply = by * trans0 + *ply;
    *plz = bz * trans0 + *plz;
    return 0;
} /* bbkaon_ */

/* ***************************************** */
/* for pion+pion-->K+K- */
/*      real*4 function pipik(srt) */
doublereal pipik_(real *srt)
{
    /* Initialized data */

    static real xarray[5] = { .001f,.7f,1.5f,1.7f,2.f };
    static real earray[5] = { 1.f,1.2f,1.6f,2.f,2.4f };

    /* System generated locals */
    real ret_val;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static integer ie;
    static real xmin, ymin, xmax, ymax, pmass;

/*  srt    = DSQRT(s) in GeV                                                   * */
/*  xsec   = production cross section in mb                                    * */
/*  NOTE: DEVIDE THE CROSS SECTION TO OBTAIN K+ PRODUCTION                     * */
/* ***************************************** */
/*      real*4   xarray(5), earray(5) */
    pmass = .9383f;
/* 1.Calculate p(lab)  from srt [GeV] */
/*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1) */
/*      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.) */
    ret_val = 0.f;
    if (*srt <= 1.f) {
	return ret_val;
    }
    if (*srt > 2.4f) {
	ret_val = 1.f;
	return ret_val;
    }
    if (*srt < earray[0]) {
	ret_val = xarray[0] / 2.f;
	return ret_val;
    }

/* 2.Interpolate double logarithmically to find sigma(srt) */

    for (ie = 1; ie <= 5; ++ie) {
	if (earray[ie - 1] == *srt) {
	    ret_val = xarray[ie - 1];
	    goto L10;
	} else if (earray[ie - 1] > *srt) {
	    ymin = log(xarray[ie - 2]);
	    ymax = log(xarray[ie - 1]);
	    xmin = log(earray[ie - 2]);
	    xmax = log(earray[ie - 1]);
	    ret_val = exp(ymin + (log(*srt) - xmin) * (ymax - ymin) / (xmax - 
		    xmin));
	    goto L10;
	}
/* L1001: */
    }
L10:
    ret_val /= 2.f;
    return ret_val;
} /* pipik_ */

/* ********************************* */
/* TOTAL PION-P INELASTIC CROSS SECTION */
/*  from the CERN data book */
/*  date: Sept.2, 1994 */
/*  for pion++p-->Delta+pion */
/*      real*4 function pionpp(srt) */
doublereal pionpp_(real *srt)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4, r__5;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal), pow_dd(doublereal *, doublereal 
	    *);

    /* Local variables */
    static real a, b, c__, d__, an, plab, pmin, pmax, pmass, pmass1;

/*  srt    = DSQRT(s) in GeV                                                   * */
/*  xsec   = production cross section in fm**2                                 * */
/*  earray = EXPerimental table with proton energies in MeV                    * */
/*  xarray = EXPerimental table with cross sections in mb (curve to guide eye) * */
/*                                                                             * */
/* ***************************************** */
    pmass = .14f;
    pmass1 = .938f;
    ret_val = 1e-5f;
    if (*srt <= 1.22f) {
	return ret_val;
    }
/* 1.Calculate p(lab)  from srt [GeV] */
/*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1) */
/*      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.) */
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = pmass;
/* Computing 2nd power */
    r__4 = pmass1;
/* Computing 2nd power */
    r__1 = (r__2 * r__2 - r__3 * r__3 - r__4 * r__4) / (pmass1 * 2.f);
/* Computing 2nd power */
    r__5 = pmass;
    plab = sqrt(r__1 * r__1 - r__5 * r__5);
    pmin = .3f;
    pmax = 25.f;
    if (plab > pmax) {
	ret_val = 2.f;
	return ret_val;
    }
    if (plab < pmin) {
	ret_val = 0.f;
	return ret_val;
    }
/* * fit parameters */
    a = 24.3f;
    b = -12.3f;
    c__ = .324f;
    an = -1.91f;
    d__ = -2.44f;
    d__1 = (doublereal) plab;
    d__2 = (doublereal) an;
/* Computing 2nd power */
    r__1 = log(plab);
    ret_val = a + b * pow_dd(&d__1, &d__2) + c__ * (r__1 * r__1) + d__ * log(
	    plab);
    if (ret_val <= 0.f) {
	ret_val = 0.f;
    }
    ret_val /= 10.f;
    return ret_val;
} /* pionpp_ */

/* ********************************* */
/* elementary cross sections */
/*  from the CERN data book */
/*  date: Sept.2, 1994 */
/*  for pion-+p-->INELASTIC */
/*      real*4 function pipp1(srt) */
doublereal pipp1_(real *srt)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4, r__5;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal), pow_dd(doublereal *, doublereal 
	    *);

    /* Local variables */
    static real a, b, c__, d__, an, plab, pmin, pmax, pmass, pmass1;

/*  srt    = DSQRT(s) in GeV                                                   * */
/*  xsec   = production cross section in fm**2                                 * */
/*  earray = EXPerimental table with proton energies in MeV                    * */
/*  xarray = EXPerimental table with cross sections in mb (curve to guide eye) * */
/*  UNITS: FM**2 */
/* ***************************************** */
    pmass = .14f;
    pmass1 = .938f;
    ret_val = 1e-4f;
    if (*srt <= 1.22f) {
	return ret_val;
    }
/* 1.Calculate p(lab)  from srt [GeV] */
/*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1) */
/*      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.) */
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = pmass;
/* Computing 2nd power */
    r__4 = pmass1;
/* Computing 2nd power */
    r__1 = (r__2 * r__2 - r__3 * r__3 - r__4 * r__4) / (pmass1 * 2.f);
/* Computing 2nd power */
    r__5 = pmass;
    plab = sqrt(r__1 * r__1 - r__5 * r__5);
    pmin = .3f;
    pmax = 25.f;
    if (plab > pmax) {
	ret_val = 2.f;
	return ret_val;
    }
    if (plab < pmin) {
	ret_val = 0.f;
	return ret_val;
    }
/* * fit parameters */
    a = 26.6f;
    b = -7.18f;
    c__ = .327f;
    an = -1.86f;
    d__ = -2.81f;
    d__1 = (doublereal) plab;
    d__2 = (doublereal) an;
/* Computing 2nd power */
    r__1 = log(plab);
    ret_val = a + b * pow_dd(&d__1, &d__2) + c__ * (r__1 * r__1) + d__ * log(
	    plab);
    if (ret_val <= 0.f) {
	ret_val = 0.f;
    }
    ret_val /= 10.f;
    return ret_val;
} /* pipp1_ */

/* ***************************** */
/*       real*4 function xrho(srt) */
doublereal xrho_(real *srt)
{
    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real es, trho, xrho0, esmin, pmass, rmass;

/*       xsection for pp-->pp+rho */
/* ***************************** */
    pmass = .9383f;
    rmass = .77f;
    trho = .151f;
    ret_val = 1e-9f;
    if (*srt <= 2.67f) {
	return ret_val;
    }
    esmin = rmass + 1.8766f - trho / 2.f;
    es = *srt;
/* the cross section for tho0 production is */
/* Computing 2nd power */
    r__1 = es - esmin;
    xrho0 = (es - esmin) * .24f / (r__1 * r__1 + 1.4f);
    ret_val = xrho0 * 3.f;
    return ret_val;
} /* xrho_ */

/* ***************************** */
/*       real*4 function omega(srt) */
doublereal omega_(real *srt)
{
    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real es, esmin, omass, pmass, tomega;

/*       xsection for pp-->pp+omega */
/* ***************************** */
    pmass = .9383f;
    omass = .782f;
    tomega = .0084f;
    ret_val = 1e-8f;
    if (*srt <= 2.68f) {
	return ret_val;
    }
    esmin = omass + 1.8766f - tomega / 2.f;
    es = *srt;
/* Computing 2nd power */
    r__1 = es - esmin;
    ret_val = (es - esmin) * .36f / (r__1 * r__1 + 1.25f);
    return ret_val;
} /* omega_ */

/* ***************************************** */
/* for ppi(+)-->DELTA+pi */
/*      real*4 function TWOPI(srt) */
doublereal twopi_(real *srt)
{
    /* Initialized data */

    static real xarray[19] = { 3e-6f,1.87f,11.f,14.9f,9.35f,7.65f,4.62f,3.45f,
	    2.41f,1.85f,1.65f,1.5f,1.32f,1.17f,1.16f,1.f,.856f,.745f,3e-6f };
    static real earray[19] = { 1.22f,1.47f,1.72f,1.97f,2.22f,2.47f,2.72f,
	    2.97f,3.22f,3.47f,3.72f,3.97f,4.22f,4.47f,4.72f,4.97f,5.22f,5.47f,
	    5.72f };

    /* System generated locals */
    real ret_val;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static integer ie;
    static real plab, xmin, ymin, xmax, ymax, pmass, pmass1;

/*  This function contains the experimental pi+p-->DELTA+PION cross sections   * */
/*  srt    = DSQRT(s) in GeV                                                   * */
/*  xsec   = production cross section in mb                                    * */
/*  earray = EXPerimental table with proton energies in MeV                    * */
/*  xarray = EXPerimental table with cross sections in mb (curve to guide eye) * */
/*                                                                             * */
/* ***************************************** */
/*      real*4   xarray(19), earray(19) */
    pmass = .14f;
    pmass1 = .938f;
    ret_val = 1e-6f;
    if (*srt <= 1.22f) {
	return ret_val;
    }
/* 1.Calculate p(lab)  from srt [GeV] */
/*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1) */
    plab = *srt;
    if (plab < earray[0]) {
	ret_val = 1e-5f;
	return ret_val;
    }

/* 2.Interpolate double logarithmically to find sigma(srt) */

    for (ie = 1; ie <= 19; ++ie) {
	if (earray[ie - 1] == plab) {
	    ret_val = xarray[ie - 1];
	    return ret_val;
	} else if (earray[ie - 1] > plab) {
	    ymin = log(xarray[ie - 2]);
	    ymax = log(xarray[ie - 1]);
	    xmin = log(earray[ie - 2]);
	    xmax = log(earray[ie - 1]);
	    ret_val = exp(ymin + (log(plab) - xmin) * (ymax - ymin) / (xmax - 
		    xmin));
	    return ret_val;
	}
/* L1001: */
    }
    return ret_val;
} /* twopi_ */

/* ***************************************** */
/* ***************************************** */
/* for ppi(+)-->DELTA+RHO */
/*      real*4 function THREPI(srt) */
doublereal threpi_(real *srt)
{
    /* Initialized data */

    static real xarray[15] = { 8e-6f,6.1999999e-5f,1.88194f,5.02569f,
	    11.80154f,13.92114f,15.07308f,11.79571f,11.53772f,10.01197f,
	    9.792673f,9.465264f,8.97049f,7.944254f,6.88632f };
    static real earray[15] = { 1.22f,1.47f,1.72f,1.97f,2.22f,2.47f,2.72f,
	    2.97f,3.22f,3.47f,3.72f,3.97f,4.22f,4.47f,4.72f };

    /* System generated locals */
    real ret_val;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static integer ie;
    static real plab, xmin, ymin, xmax, ymax, pmass, pmass1;

/*  This function contains the experimental pi+p-->DELTA + rho cross sections  * */
/*  srt    = DSQRT(s) in GeV                                                   * */
/*  xsec   = production cross section in mb                                    * */
/*  earray = EXPerimental table with proton energies in MeV                    * */
/*  xarray = EXPerimental table with cross sections in mb (curve to guide eye) * */
/*                                                                             * */
/* ***************************************** */
/*      real*4   xarray(15), earray(15) */
    pmass = .14f;
    pmass1 = .938f;
    ret_val = 1e-6f;
    if (*srt <= 1.36f) {
	return ret_val;
    }
/* 1.Calculate p(lab)  from srt [GeV] */
/*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1) */
    plab = *srt;
    if (plab < earray[0]) {
	ret_val = 1e-5f;
	return ret_val;
    }

/* 2.Interpolate double logarithmically to find sigma(srt) */

    for (ie = 1; ie <= 15; ++ie) {
	if (earray[ie - 1] == plab) {
	    ret_val = xarray[ie - 1];
	    return ret_val;
	} else if (earray[ie - 1] > plab) {
	    ymin = log(xarray[ie - 2]);
	    ymax = log(xarray[ie - 1]);
	    xmin = log(earray[ie - 2]);
	    xmax = log(earray[ie - 1]);
	    ret_val = exp(ymin + (log(plab) - xmin) * (ymax - ymin) / (xmax - 
		    xmin));
	    return ret_val;
	}
/* L1001: */
    }
    return ret_val;
} /* threpi_ */

/* ***************************************** */
/* ***************************************** */
/* for ppi(+)-->DELTA+omega */
/*      real*4 function FOURPI(srt) */
doublereal fourpi_(real *srt)
{
    /* Initialized data */

    static real xarray[10] = { 1e-4f,1.986597f,6.411932f,7.636956f,9.598362f,
	    9.88974f,10.24317f,10.80138f,11.86988f,12.83925f };
    static real earray[10] = { 2.468f,2.718f,2.968f,3.22f,3.47f,3.72f,3.97f,
	    4.22f,4.47f,4.72f };

    /* System generated locals */
    real ret_val;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static integer ie;
    static real plab, xmin, ymin, xmax, ymax, pmass, pmass1;

/*  This function contains the experimental pi+p-->DELTA+PION cross sections   * */
/*  srt    = DSQRT(s) in GeV                                                   * */
/*  xsec   = production cross section in mb                                    * */
/*  earray = EXPerimental table with proton energies in MeV                    * */
/*  xarray = EXPerimental table with cross sections in mb (curve to guide eye) * */
/*                                                                             * */
/* ***************************************** */
/*      real*4   xarray(10), earray(10) */
    pmass = .14f;
    pmass1 = .938f;
    ret_val = 1e-6f;
    if (*srt <= 1.52f) {
	return ret_val;
    }
/* 1.Calculate p(lab)  from srt [GeV] */
/*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1) */
    plab = *srt;
    if (plab < earray[0]) {
	ret_val = 1e-5f;
	return ret_val;
    }

/* 2.Interpolate double logarithmically to find sigma(srt) */

    for (ie = 1; ie <= 10; ++ie) {
	if (earray[ie - 1] == plab) {
	    ret_val = xarray[ie - 1];
	    return ret_val;
	} else if (earray[ie - 1] > plab) {
	    ymin = log(xarray[ie - 2]);
	    ymax = log(xarray[ie - 1]);
	    xmin = log(earray[ie - 2]);
	    xmax = log(earray[ie - 1]);
	    ret_val = exp(ymin + (log(plab) - xmin) * (ymax - ymin) / (xmax - 
		    xmin));
	    return ret_val;
	}
/* L1001: */
    }
    return ret_val;
} /* fourpi_ */

/* ***************************************** */
/* ***************************************** */
/* for pion (rho or omega)+baryon resonance collisions */
/*      real*4 function reab(i1,i2,srt,ictrl) */
doublereal reab_(integer *i1, integer *i2, real *srt, integer *ictrl)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4, r__5;

    /* Local variables */
    static real ed;
    static integer lb1, lb2;
    static real pin2, xpro, arho1, pout2;
    extern doublereal twopi_(real *);
    static real factor;
    extern doublereal threpi_(real *), fourpi_(real *);

/*  This function calculates the cross section for */
/*  pi+Delta(N*)-->N+PION process                                              * */
/*  srt    = DSQRT(s) in GeV                                                   * */
/*  reab   = cross section in fm**2                                            * */
/*  ictrl=1,2,3 for pion, rho and omega+D(N*) */
/* *************************************** */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /DD/ */
/* c      SAVE /EE/ */
    lb1 = ee_1.lb[*i1 - 1];
    lb2 = ee_1.lb[*i2 - 1];
    ret_val = 0.f;
    if (*ictrl == 1 && *srt <= 1.238f) {
	return ret_val;
    }
    if (*ictrl == 3 && *srt <= 1.8799999999999999f) {
	return ret_val;
    }
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__1 = (r__2 * r__2 + .019600000000000003f - .87984399999999985f) / (*srt 
	    * 2.f);
    pin2 = r__1 * r__1 - .019600000000000003f;
    if (pin2 <= 0.f) {
	return ret_val;
    }
/* for pion+D(N*)-->pion+N */
    if (*ictrl == 1) {
	if (cc_1.e[*i1 - 1] > 1.f) {
	    ed = cc_1.e[*i1 - 1];
	} else {
	    ed = cc_1.e[*i2 - 1];
	}
/* Computing 2nd power */
	r__2 = *srt;
/* Computing 2nd power */
	r__3 = ed;
/* Computing 2nd power */
	r__1 = (r__2 * r__2 + .019600000000000003f - r__3 * r__3) / (*srt * 
		2.f);
	pout2 = r__1 * r__1 - .019600000000000003f;
	if (pout2 <= 0.f) {
	    return ret_val;
	}
	xpro = twopi_(srt) / 10.f;
	factor = .33333333333333331f;
	if (lb1 == 8 && lb2 == 5 || lb1 == 5 && lb2 == 8 || (lb1 == -8 && lb2 
		== 3 || lb1 == 3 && lb2 == -8)) {
	    factor = .25f;
	}
	if (abs(lb1) >= 10 && abs(lb1) <= 13 || abs(lb2) >= 10 && abs(lb2) <= 
		13) {
	    factor = 1.f;
	}
	ret_val = factor * pin2 / pout2 * xpro;
	return ret_val;
    }
/* for rho reabsorption */
    if (*ictrl == 2) {
	if (ee_1.lb[*i2 - 1] >= 25) {
	    ed = cc_1.e[*i1 - 1];
	    arho1 = cc_1.e[*i2 - 1];
	} else {
	    ed = cc_1.e[*i2 - 1];
	    arho1 = cc_1.e[*i1 - 1];
	}
	if (*srt <= arho1 + 1.0779999999999998f + .02f) {
	    return ret_val;
	}
/* Computing 2nd power */
	r__2 = *srt;
/* Computing 2nd power */
	r__3 = arho1;
/* Computing 2nd power */
	r__4 = ed;
/* Computing 2nd power */
	r__1 = (r__2 * r__2 + r__3 * r__3 - r__4 * r__4) / (*srt * 2.f);
/* Computing 2nd power */
	r__5 = arho1;
	pout2 = r__1 * r__1 - r__5 * r__5;
	if (pout2 <= 0.f) {
	    return ret_val;
	}
	xpro = threpi_(srt) / 10.f;
	factor = .33333333333333331f;
	if (lb1 == 8 && lb2 == 27 || lb1 == 27 && lb2 == 8 || (lb1 == -8 && 
		lb2 == 25 || lb1 == 25 && lb2 == -8)) {
	    factor = .25f;
	}
	if (abs(lb1) >= 10 && abs(lb1) <= 13 || abs(lb2) >= 10 && abs(lb2) <= 
		13) {
	    factor = 1.f;
	}
	ret_val = factor * pin2 / pout2 * xpro;
	return ret_val;
    }
/* for omega reabsorption */
    if (*ictrl == 3) {
	if (cc_1.e[*i1 - 1] > 1.f) {
	    ed = cc_1.e[*i1 - 1];
	}
	if (cc_1.e[*i2 - 1] > 1.f) {
	    ed = cc_1.e[*i2 - 1];
	}
/* Computing 2nd power */
	r__2 = *srt;
/* Computing 2nd power */
	r__3 = ed;
/* Computing 2nd power */
	r__1 = (r__2 * r__2 + .61152400000000007f - r__3 * r__3) / (*srt * 
		2.f);
	pout2 = r__1 * r__1 - .61152400000000007f;
	if (pout2 <= 0.f) {
	    return ret_val;
	}
	xpro = fourpi_(srt) / 10.f;
	factor = .16666666666666666f;
	if (abs(lb1) >= 10 && abs(lb1) <= 13 || abs(lb2) >= 10 && abs(lb2) <= 
		13) {
	    factor = .33333333333333331f;
	}
	ret_val = factor * pin2 / pout2 * xpro;
    }
    return ret_val;
} /* reab_ */

/* ***************************************** */
/* for the reabsorption of two resonances */
/* This function calculates the cross section for */
/* DD-->NN, N*N*-->NN and DN*-->NN */
/*      real*4 function reab2d(i1,i2,srt) */
doublereal reab2d_(integer *i1, integer *i2, real *srt)
{
    /* System generated locals */
    integer i__1;
    real ret_val, r__1, r__2, r__3, r__4, r__5;

    /* Local variables */
    static real ed1, ed2;
    static integer lb1, lb2;
    static real pin2;
    extern doublereal x2pi_(real *);
    static real xpro, pout2, factor;

/*  srt    = DSQRT(s) in GeV                                                   * */
/*  reab   = cross section in mb */
/* *************************************** */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /DD/ */
/* c      SAVE /EE/ */
    ret_val = 0.f;
    lb1 = (i__1 = ee_1.lb[*i1 - 1], abs(i__1));
    lb2 = (i__1 = ee_1.lb[*i2 - 1], abs(i__1));
    ed1 = cc_1.e[*i1 - 1];
    ed2 = cc_1.e[*i2 - 1];
/* Computing 2nd power */
    r__1 = *srt / 2.f;
    pin2 = r__1 * r__1 - .87984399999999985f;
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = ed1;
/* Computing 2nd power */
    r__4 = ed2;
/* Computing 2nd power */
    r__1 = (r__2 * r__2 + r__3 * r__3 - r__4 * r__4) / (*srt * 2.f);
/* Computing 2nd power */
    r__5 = ed1;
    pout2 = r__1 * r__1 - r__5 * r__5;
    if (pout2 <= 0.f) {
	return ret_val;
    }
    xpro = x2pi_(srt);
    factor = .25f;
    if (lb1 >= 10 && lb1 <= 13 && (lb2 >= 10 && lb2 <= 13)) {
	factor = 1.f;
    }
    if (lb1 >= 6 && lb1 <= 9 && (lb2 > 10 && lb2 <= 13)) {
	factor = .5f;
    }
    if (lb2 >= 6 && lb2 <= 9 && (lb1 > 10 && lb1 <= 13)) {
	factor = .5f;
    }
    ret_val = factor * pin2 / pout2 * xpro;
    return ret_val;
} /* reab2d_ */

/* ************************************** */
/* Subroutine */ int rotate_(real *px0, real *py0, real *pz0, real *px, real *
	py, real *pz)
{
    /* System generated locals */
    real r__1, r__2, r__3;

    /* Builtin functions */
    double sqrt(doublereal), atan2(doublereal, doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static real c1, c2, s1, s2, t2, t1, pr, ss, ct1, ct2, pr0, st2, st1, 
	    scheck;

    /* Fortran I/O blocks */
    static cilist io___1966 = { 0, 99, 0, 0, 0 };
    static cilist io___1973 = { 0, 99, 0, 0, 0 };


/* purpose: rotate the momentum of a particle in the CMS of p1+p2 such that */
/* the x' y' and z' in the cms of p1+p2 is the same as the fixed x y and z */
/* quantities: */
/*            px0,py0 and pz0 are the cms momentum of the incoming colliding */
/*            particles */
/*            px, py and pz are the cms momentum of any one of the particles */
/*            after the collision to be rotated */
/* ************************************** */
/* the momentum, polar and azimuthal angles of the incoming momentm */
/* Computing 2nd power */
    r__1 = *px0;
/* Computing 2nd power */
    r__2 = *py0;
/* Computing 2nd power */
    r__3 = *pz0;
    pr0 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
    if (pr0 == 0.f) {
	pr0 = 1e-8f;
    }
    c2 = *pz0 / pr0;
    if (*px0 == 0.f && *py0 == 0.f) {
	t2 = 0.f;
    } else {
	t2 = atan2(*py0, *px0);
    }
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__1 = c2;
    scheck = 1.f - r__1 * r__1;
    if (scheck < 0.f) {
	s_wsle(&io___1966);
	do_lio(&c__9, &c__1, "scheck45: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	scheck = 0.f;
    }
    s2 = sqrt(scheck);
/*      S2  =  SQRT( 1.0 - C2**2 ) */
    ct2 = cos(t2);
    st2 = sin(t2);
/* the momentum, polar and azimuthal angles of the momentum to be rotated */
/* Computing 2nd power */
    r__1 = *px;
/* Computing 2nd power */
    r__2 = *py;
/* Computing 2nd power */
    r__3 = *pz;
    pr = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
    if (pr == 0.f) {
	pr = 1e-7f;
    }
    c1 = *pz / pr;
    if (*px == 0.f && *py == 0.f) {
	t1 = 0.f;
    } else {
	t1 = atan2(*py, *px);
    }
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__1 = c1;
    scheck = 1.f - r__1 * r__1;
    if (scheck < 0.f) {
	s_wsle(&io___1973);
	do_lio(&c__9, &c__1, "scheck46: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	scheck = 0.f;
    }
    s1 = sqrt(scheck);
/*      S1   = SQRT( 1.0 - C1**2 ) */
    ct1 = cos(t1);
    st1 = sin(t1);
    ss = c2 * s1 * ct1 + s2 * c1;
/* THE MOMENTUM AFTER ROTATION */
    *px = pr * (ss * ct2 - s1 * st1 * st2);
    *py = pr * (ss * st2 + s1 * st1 * ct2);
    *pz = pr * (c1 * c2 - s1 * s2 * ct1);
    return 0;
} /* rotate_ */

/* ***************************************** */
/*      real*4 function Xpp(srt) */
doublereal xpp_(real *srt)
{
    /* Initialized data */

    static real earray[14] = { 20.f,30.f,40.f,60.f,80.f,100.f,170.f,250.f,
	    310.f,350.f,460.f,560.f,660.f,800.f };
    static real xarray[14] = { 150.f,90.f,80.6f,48.f,36.6f,31.6f,25.9f,24.f,
	    23.1f,24.f,28.3f,33.6f,41.5f,47.f };

    /* System generated locals */
    real ret_val, r__1;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static integer ie;
    static real ekin, xmin, ymin, xmax, ymax, pmass;

/*  This function contains the experimental total n-p cross sections           * */
/*  srt    = DSQRT(s) in GeV                                                   * */
/*  xsec   = production cross section in mb                                    * */
/*  earray = EXPerimental table with proton energies in MeV                    * */
/*  xarray = EXPerimental table with cross sections in mb (curve to guide eye) * */
/*  WITH A CUTOFF AT 55MB                                                      * */
/* ***************************************** */
/*      real*4   xarray(14), earray(14) */
    pmass = .9383f;
/* 1.Calculate E_kin(lab) [MeV] from srt [GeV] */
/*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1) */
/* Computing 2nd power */
    r__1 = *srt / (pmass * 2.f);
    ekin = pmass * 2e3f * (r__1 * r__1 - 1.f);
    if (ekin < earray[0]) {
	ret_val = xarray[0];
	if (ret_val > 55.f) {
	    ret_val = 55.f;
	}
	return ret_val;
    }
    if (ekin > earray[13]) {
	ret_val = xarray[13];
	return ret_val;
    }


/* 2.Interpolate double logarithmically to find sigma(srt) */

    for (ie = 1; ie <= 14; ++ie) {
	if (earray[ie - 1] == ekin) {
	    ret_val = xarray[ie - 1];
	    if (ret_val > 55.f) {
		ret_val = 55.f;
	    }
	    return ret_val;
	}
	if (earray[ie - 1] > ekin) {
	    ymin = log(xarray[ie - 2]);
	    ymax = log(xarray[ie - 1]);
	    xmin = log(earray[ie - 2]);
	    xmax = log(earray[ie - 1]);
	    ret_val = exp(ymin + (log(ekin) - xmin) * (ymax - ymin) / (xmax - 
		    xmin));
	    if (ret_val > 55.f) {
		ret_val = 55.f;
	    }
	    goto L50;
	}
/* L1001: */
    }
L50:
    return ret_val;
} /* xpp_ */

/* ***************************************** */
doublereal xnp_(real *srt)
{
    /* Initialized data */

    static real earray[11] = { 20.f,30.f,40.f,60.f,90.f,135.f,200.f,300.f,
	    400.f,600.f,800.f };
    static real xarray[11] = { 410.f,270.f,214.5f,130.f,78.f,53.5f,41.6f,
	    35.9f,34.2f,34.3f,34.9f };

    /* System generated locals */
    real ret_val, r__1;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static integer ie;
    static real ekin, xmin, ymin, xmax, ymax, pmass;

/*  This function contains the experimental total n-p cross sections           * */
/*  srt    = DSQRT(s) in GeV                                                   * */
/*  xsec   = production cross section in mb                                    * */
/*  earray = EXPerimental table with proton energies in MeV                    * */
/*  xarray = EXPerimental table with cross sections in mb (curve to guide eye) * */
/*  WITH  A CUTOFF AT 55MB                                                * */
/* ***************************************** */
/*      real*4   xarray(11), earray(11) */
    pmass = .9383f;
/* 1.Calculate E_kin(lab) [MeV] from srt [GeV] */
/*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1) */
/* Computing 2nd power */
    r__1 = *srt / (pmass * 2.f);
    ekin = pmass * 2e3f * (r__1 * r__1 - 1.f);
    if (ekin < earray[0]) {
	ret_val = xarray[0];
	if (ret_val > 55.f) {
	    ret_val = 55.f;
	}
	return ret_val;
    }
    if (ekin > earray[10]) {
	ret_val = xarray[10];
	return ret_val;
    }

/* Interpolate double logarithmically to find sigma(srt) */

    for (ie = 1; ie <= 11; ++ie) {
	if (earray[ie - 1] == ekin) {
	    ret_val = xarray[ie - 1];
	    if (ret_val > 55.f) {
		ret_val = 55.f;
	    }
	    return ret_val;
	}
	if (earray[ie - 1] > ekin) {
	    ymin = log(xarray[ie - 2]);
	    ymax = log(xarray[ie - 1]);
	    xmin = log(earray[ie - 2]);
	    xmax = log(earray[ie - 1]);
	    ret_val = exp(ymin + (log(ekin) - xmin) * (ymax - ymin) / (xmax - 
		    xmin));
	    if (ret_val > 55.f) {
		ret_val = 55.f;
	    }
	    goto L50;
	}
/* L1001: */
    }
L50:
    return ret_val;
} /* xnp_ */

/* ****************************** */
doublereal ptr_(real *ptmax, integer *iseed)
{
    /* System generated locals */
    real ret_val;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static integer ie;
    static real xt, xmin, ymin, xmax, ymax;
    extern doublereal ptdis_(real *), ranart_(integer *);
    static real tryial;

/* (2) Generate the transverse momentum */
/*     OF nucleons */
/* ****************************** */
/* c      SAVE /TABLE/ */
/* c      SAVE /RNDF77/ */
    ret_val = 0.f;
    if (*ptmax <= .01f) {
	ret_val = *ptmax;
	return ret_val;
    }
    if (*ptmax > 2.01f) {
	*ptmax = 2.01f;
    }
    tryial = ptdis_(ptmax) / ptdis_(&c_b842);
    xt = ranart_(&rndf77_1.nseed) * tryial;
/* look up the table and */
/* Interpolate double logarithmically to find pt */
    for (ie = 1; ie <= 200; ++ie) {
	if (table_1.earray[ie] == xt) {
	    ret_val = table_1.xarray[ie];
	    return ret_val;
	}
	if (table_1.xarray[ie - 1] <= 1e-5f) {
	    goto L50;
	}
	if (table_1.xarray[ie] <= 1e-5f) {
	    goto L50;
	}
	if (table_1.earray[ie - 1] <= 1e-5f) {
	    goto L50;
	}
	if (table_1.earray[ie] <= 1e-5f) {
	    goto L50;
	}
	if (table_1.earray[ie] > xt) {
	    ymin = log(table_1.xarray[ie - 1]);
	    ymax = log(table_1.xarray[ie]);
	    xmin = log(table_1.earray[ie - 1]);
	    xmax = log(table_1.earray[ie]);
	    ret_val = exp(ymin + (log(xt) - xmin) * (ymax - ymin) / (xmax - 
		    xmin));
	    if (ret_val > *ptmax) {
		ret_val = *ptmax;
	    }
	    return ret_val;
	}
L50:
	;
    }
    return ret_val;
} /* ptr_ */

/* ********************************* */
/* ********************************* */
/*                                                                      * */
/*                                                                      * */
/* Subroutine */ int xnd_(real *px, real *py, real *pz, real *srt, integer *
	i1, integer *i2, real *xinel, real *sigk, real *xsk1, real *xsk2, 
	real *xsk3, real *xsk4, real *xsk5)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8;
    real r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real al, as, es, pr, ak0, em1, em2, ada, ana;
    extern /* Subroutine */ int m1535_(integer *, integer *, real *, real *);
    static real akp, x1440, x1535, prf;
    extern doublereal ppk0_(real *), ppk1_(real *);
    static real t1dlk, t2dlk, t1dsk, t2dsk, t1nlk, t1nsk, t2nsk;
    extern doublereal sigma_(real *, integer *, integer *, integer *), denom_(
	    real *, real *);
    static real signd, pmdlk, sigdn, pmdsk, renom;
    extern doublereal pplpk_(real *);
    static real pmnsk, pmdlk2, pmdsk2, renom1, pmnsk2, deltam, renomn;

/*     PURPOSE:                                                         * */
/*             calculate NUCLEON-BARYON RESONANCE inelatic Xsection     * */
/*     NOTE   :                                                         * */
/*     QUANTITIES:                                                 * */
/*                      CHANNELS. M12 IS THE REVERSAL CHANNEL OF N12    * */
/*                      N12,                                            * */
/*                      M12=1 FOR p+n-->delta(+)+ n                     * */
/*                          2     p+n-->delta(0)+ p                     * */
/*                          3     p+p-->delta(++)+n                     * */
/*                          4     p+p-->delta(+)+p                      * */
/*                          5     n+n-->delta(0)+n                      * */
/*                          6     n+n-->delta(-)+p                      * */
/*                          7     n+p-->N*(0)(1440)+p                   * */
/*                          8     n+p-->N*(+)(1440)+n                   * */
/*                        9     p+p-->N*(+)(1535)+p                     * */
/*                        10    n+n-->N*(0)(1535)+n                     * */
/*                         11    n+p-->N*(+)(1535)+n                     * */
/*                        12    n+p-->N*(0)(1535)+p */
/*                        13    D(++)+D(-)-->N*(+)(1440)+n */
/*                         14    D(++)+D(-)-->N*(0)(1440)+p */
/*                        15    D(+)+D(0)--->N*(+)(1440)+n */
/*                        16    D(+)+D(0)--->N*(0)(1440)+p */
/*                        17    D(++)+D(0)-->N*(+)(1535)+p */
/*                        18    D(++)+D(-)-->N*(0)(1535)+p */
/*                        19    D(++)+D(-)-->N*(+)(1535)+n */
/*                        20    D(+)+D(+)-->N*(+)(1535)+p */
/*                        21    D(+)+D(0)-->N*(+)(1535)+n */
/*                        22    D(+)+D(0)-->N*(0)(1535)+p */
/*                        23    D(+)+D(-)-->N*(0)(1535)+n */
/*                        24    D(0)+D(0)-->N*(0)(1535)+n */
/*                          25    N*(+)(14)+N*(+)(14)-->N*(+)(15)+p */
/*                          26    N*(0)(14)+N*(0)(14)-->N*(0)(15)+n */
/*                          27    N*(+)(14)+N*(0)(14)-->N*(+)(15)+n */
/*                        28    N*(+)(14)+N*(0)(14)-->N*(0)(15)+p */
/*                        29    N*(+)(14)+D+-->N*(+)(15)+p */
/*                        30    N*(+)(14)+D0-->N*(+)(15)+n */
/*                        31    N*(+)(14)+D(-)-->N*(0)(1535)+n */
/*                        32    N*(0)(14)+D++--->N*(+)(15)+p */
/*                        33    N*(0)(14)+D+--->N*(+)(15)+n */
/*                        34    N*(0)(14)+D+--->N*(0)(15)+p */
/*                        35    N*(0)(14)+D0-->N*(0)(15)+n */
/*                        36    N*(+)(14)+D0--->N*(0)(15)+p */
/*                            and more */
/* ********************************** */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /ff/ */
/* c      SAVE /gg/ */
/* c      SAVE /INPUT/ */
/* c      SAVE /NN/ */
/* c      SAVE /BG/ */
/* c      SAVE /RUN/ */
/* c      SAVE /PA/ */
/* c      SAVE /PB/ */
/* c      SAVE /PC/ */
/* c      SAVE /PD/ */
/* c      SAVE /input1/ */
/* ----------------------------------------------------------------------- */
    *xinel = 0.f;
    *sigk = 0.f;
    *xsk1 = 0.f;
    *xsk2 = 0.f;
    *xsk3 = 0.f;
    *xsk4 = 0.f;
    *xsk5 = 0.f;
    em1 = cc_1.e[*i1 - 1];
    em2 = cc_1.e[*i2 - 1];
/* Computing 2nd power */
    r__1 = *px;
/* Computing 2nd power */
    r__2 = *py;
/* Computing 2nd power */
    r__3 = *pz;
    pr = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
/*     CAN HAPPEN ANY MORE ==> RETURN (2.04 = 2*AVMASS + PI-MASS+0.02) */
    if (*srt < 2.04f) {
	return 0;
    }
/* Resonance absorption or Delta + N-->N*(1440), N*(1535) */
/* COM: TEST FOR DELTA OR N* ABSORPTION */
/*      IN THE PROCESS DELTA+N-->NN, N*+N-->NN */
/* Computing 2nd power */
    r__1 = *srt;
    prf = sqrt(r__1 * r__1 * .25f - .88040689000000005f);
    if (em1 > 1.f) {
	deltam = em1;
    } else {
	deltam = em2;
    }
/* Computing 2nd power */
    r__1 = prf;
    renom = deltam * (r__1 * r__1) / denom_(srt, &c_b183) / pr;
/* Computing 2nd power */
    r__1 = prf;
    renomn = deltam * (r__1 * r__1) / denom_(srt, &c_b254) / pr;
/* Computing 2nd power */
    r__1 = prf;
    renom1 = deltam * (r__1 * r__1) / denom_(srt, &c_b255) / pr;
/* avoid the inelastic collisions between n+delta- -->N+N */
/*       and p+delta++ -->N+N due to charge conservation, */
/*       but they can scatter to produce kaons */
    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 && (i__2 = ee_1.lb[*i2 - 1],
	     abs(i__2)) == 6) {
	renom = 0.f;
    }
    if ((i__1 = ee_1.lb[*i2 - 1], abs(i__1)) == 2 && (i__2 = ee_1.lb[*i1 - 1],
	     abs(i__2)) == 6) {
	renom = 0.f;
    }
    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1 && (i__2 = ee_1.lb[*i2 - 1],
	     abs(i__2)) == 9) {
	renom = 0.f;
    }
    if ((i__1 = ee_1.lb[*i2 - 1], abs(i__1)) == 1 && (i__2 = ee_1.lb[*i1 - 1],
	     abs(i__2)) == 9) {
	renom = 0.f;
    }
    i__3 = (i__1 = ee_1.lb[*i1 - 1], abs(i__1));
    i__4 = (i__2 = ee_1.lb[*i2 - 1], abs(i__2));
    m1535_(&i__3, &i__4, srt, &x1535);
    x1440 = sigma_(srt, &c__2, &c__0, &c__1) * .75f;
/* CROSS SECTION FOR KAON PRODUCTION from the four channels */
/* for NLK channel */
    akp = .498f;
    ak0 = .498f;
    ana = .94f;
    ada = 1.232f;
    al = 1.1157f;
    as = 1.1197f;
    *xsk1 = 0.f;
    *xsk2 = 0.f;
    *xsk3 = 0.f;
    *xsk4 = 0.f;
/*      !! phi production */
    *xsk5 = 0.f;
    t1nlk = ana + al + akp;
    if (*srt <= t1nlk) {
	goto L222;
    }
    *xsk1 = pplpk_(srt) * 1.5f;
/* for DLK channel */
    t1dlk = ada + al + akp;
    t2dlk = ada + al - akp;
    if (*srt <= t1dlk) {
	goto L222;
    }
    es = *srt;
/* Computing 2nd power */
    r__1 = es;
/* Computing 2nd power */
    r__2 = t1dlk;
/* Computing 2nd power */
    r__3 = es;
/* Computing 2nd power */
    r__4 = t2dlk;
/* Computing 2nd power */
    r__5 = es;
    pmdlk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (
	    r__5 * r__5 * 4.f);
    pmdlk = sqrt(pmdlk2);
    *xsk3 = pplpk_(srt) * 1.5f;
/* for NSK channel */
    t1nsk = ana + as + akp;
    t2nsk = ana + as - akp;
    if (*srt <= t1nsk) {
	goto L222;
    }
/* Computing 2nd power */
    r__1 = es;
/* Computing 2nd power */
    r__2 = t1nsk;
/* Computing 2nd power */
    r__3 = es;
/* Computing 2nd power */
    r__4 = t2nsk;
/* Computing 2nd power */
    r__5 = es;
    pmnsk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (
	    r__5 * r__5 * 4.f);
    pmnsk = sqrt(pmnsk2);
    *xsk2 = (ppk1_(srt) + ppk0_(srt)) * 1.5f;
/* for DSK channel */
    t1dsk = ada + as + akp;
    t2dsk = ada + as - akp;
    if (*srt <= t1dsk) {
	goto L222;
    }
/* Computing 2nd power */
    r__1 = es;
/* Computing 2nd power */
    r__2 = t1dsk;
/* Computing 2nd power */
    r__3 = es;
/* Computing 2nd power */
    r__4 = t2dsk;
/* Computing 2nd power */
    r__5 = es;
    pmdsk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (
	    r__5 * r__5 * 4.f);
    pmdsk = sqrt(pmdsk2);
    *xsk4 = (ppk1_(srt) + ppk0_(srt)) * 1.5f;
/* sp11/21/01 */
/* phi production */
    if (*srt <= 2.898914f) {
	goto L222;
    }
/*  !! mb put the correct form */
    *xsk5 = 1e-4f;
/* sp11/21/01 end */
/* THE TOTAL KAON+ PRODUCTION CROSS SECTION IS THEN */
L222:
    *sigk = *xsk1 + *xsk2 + *xsk3 + *xsk4;
/* bz3/7/99 neutralk */
    *xsk1 *= 2.f;
    *xsk2 *= 2.f;
    *xsk3 *= 2.f;
    *xsk4 *= 2.f;
    *sigk = *sigk * 2.f + *xsk5;
/* bz3/7/99 neutralk end */
/* avoid the inelastic collisions between n+delta- -->N+N */
/*       and p+delta++ -->N+N due to charge conservation, */
/*       but they can scatter to produce kaons */
    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 && (i__2 = ee_1.lb[*i2 - 1],
	     abs(i__2)) == 6 || (i__3 = ee_1.lb[*i2 - 1], abs(i__3)) == 2 && (
	    i__4 = ee_1.lb[*i1 - 1], abs(i__4)) == 6 || (i__5 = ee_1.lb[*i1 - 
	    1], abs(i__5)) == 1 && (i__6 = ee_1.lb[*i2 - 1], abs(i__6)) == 9 
	    || (i__7 = ee_1.lb[*i2 - 1], abs(i__7)) == 1 && (i__8 = ee_1.lb[*
	    i1 - 1], abs(i__8)) == 9) {
	*xinel = *sigk;
	return 0;
    }
/* WE DETERMINE THE REACTION CHANNELS IN THE FOLLOWING */
/* FOR n+delta(++)-->p+p or n+delta(++)-->n+N*(+)(1440),n+N*(+)(1535) */
/* REABSORPTION OR N*(1535) PRODUCTION LIKE IN P+P OR N*(1440) LIKE PN, */
    if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 18 && ((i__1 = ee_1.lb[*i1 - 1]
	    , abs(i__1)) == 2 || (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 2)) {
	signd = sigma_(srt, &c__1, &c__1, &c__0) + sigma_(srt, &c__1, &c__1, &
		c__1) * .5f;
	sigdn = signd * .25f * renom;
	*xinel = sigdn + x1440 + x1535 + *sigk;
	return 0;
    }
/* FOR p+delta(-)-->n+n or p+delta(-)-->n+N*(0)(1440),n+N*(0)(1535) */
/* REABSORPTION OR N*(1535) PRODUCTION LIKE IN P+P OR N*(1440) LIKE PN, */
    if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 6 && ((i__1 = ee_1.lb[*i1 - 1],
	     abs(i__1)) == 1 || (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 1)) {
	signd = sigma_(srt, &c__1, &c__1, &c__0) + sigma_(srt, &c__1, &c__1, &
		c__1) * .5f;
	sigdn = signd * .25f * renom;
	*xinel = sigdn + x1440 + x1535 + *sigk;
	return 0;
    }
/* FOR p+delta(+)-->p+p, N*(+)(144)+p, N*(+)(1535)+p */
/* bz11/25/98 */
    if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 8 && ((i__1 = ee_1.lb[*i1 - 1],
	     abs(i__1)) == 1 || (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 1)) {
	signd = sigma_(srt, &c__1, &c__1, &c__1) * 1.5f;
	sigdn = signd * .25f * renom;
	*xinel = sigdn + x1440 + x1535 + *sigk;
	return 0;
    }
/* FOR n+delta(0)-->n+n, N*(0)(144)+n, N*(0)(1535)+n */
    if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 14 && ((i__1 = ee_1.lb[*i1 - 1]
	    , abs(i__1)) == 2 && (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 2)) {
	signd = sigma_(srt, &c__1, &c__1, &c__1) * 1.5f;
	sigdn = signd * .25f * renom;
	*xinel = sigdn + x1440 + x1535 + *sigk;
	return 0;
    }
/* FOR n+delta(+)-->n+p, N*(+)(1440)+n,N*(0)(1440)+p, */
/*                       N*(+)(1535)+n,N*(0)(1535)+p */
    if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 16 && ((i__1 = ee_1.lb[*i1 - 1]
	    , abs(i__1)) == 2 || (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 2)) {
	signd = sigma_(srt, &c__1, &c__1, &c__1) * .5f + sigma_(srt, &c__1, &
		c__1, &c__0) * .25f;
	sigdn = signd * .5f * renom;
	*xinel = sigdn + x1440 * 2.f + x1535 * 2.f + *sigk;
	return 0;
    }
/* FOR p+delta(0)-->n+p, N*(+)(1440)+n,N*(0)(1440)+p, */
/*                       N*(+)(1535)+n,N*(0)(1535)+p */
    if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 7) {
	signd = sigma_(srt, &c__1, &c__1, &c__1) * .5f + sigma_(srt, &c__1, &
		c__1, &c__0) * .25f;
	sigdn = signd * .5f * renom;
	*xinel = sigdn + x1440 * 2.f + x1535 * 2.f + *sigk;
	return 0;
    }
/* FOR p+N*(0)(14)-->p+n, N*(+)(1535)+n,N*(0)(1535)+p */
/* OR  P+N*(0)(14)-->D(+)+N, D(0)+P, */
    if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 10 && ((i__1 = ee_1.lb[*i1 - 1]
	    , abs(i__1)) == 1 || (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 1)) {
	signd = sigma_(srt, &c__2, &c__0, &c__1) * .75f;
	sigdn = signd * renomn;
	*xinel = sigdn + x1535 + *sigk;
	return 0;
    }
/* FOR n+N*(+)-->p+n, N*(+)(1535)+n,N*(0)(1535)+p */
    if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 22 && ((i__1 = ee_1.lb[*i1 - 1]
	    , abs(i__1)) == 2 || (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 2)) {
	signd = sigma_(srt, &c__2, &c__0, &c__1) * .75f;
	sigdn = signd * renomn;
	*xinel = sigdn + x1535 + *sigk;
	return 0;
    }
/* FOR N*(1535)+N-->N+N COLLISIONS */
    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 12 || (i__2 = ee_1.lb[*i1 - 1]
	    , abs(i__2)) == 13 || (i__3 = ee_1.lb[*i2 - 1], abs(i__3)) == 12 
	    || (i__4 = ee_1.lb[*i2 - 1], abs(i__4)) == 13) {
	signd = x1535;
	sigdn = signd * renom1;
	*xinel = sigdn + *sigk;
	return 0;
    }
    return 0;
} /* xnd_ */

/* ********************************* */
/*                                                                      * */
/*                                                                      * */
/* Subroutine */ int xddin_(real *px, real *py, real *pz, real *srt, integer *
	i1, integer *i2, real *xinel, real *sigk, real *xsk1, real *xsk2, 
	real *xsk3, real *xsk4, real *xsk5)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6;
    real r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real al, as, es, pr, ak0, em1, em2, s2d, ada, ana;
    static integer idd;
    extern /* Subroutine */ int n1535_(integer *, integer *, real *, real *);
    static real akp, x1535, sig2;
    extern doublereal ppk0_(real *), ppk1_(real *);
    static real t1dlk, t2dlk, t1dsk, t2dsk, t1nlk, t1nsk, t2nsk;
    extern doublereal sigma_(real *, integer *, integer *, integer *);
    static real pmdlk, pmdsk;
    extern doublereal pplpk_(real *);
    static real pmnsk;
    extern doublereal reab2d_(integer *, integer *, real *);
    static real pmdlk2, pmdsk2, pmnsk2;

/*     PURPOSE:                                                         * */
/*             DEALING WITH BARYON RESONANCE-BARYON RESONANCE COLLISIONS* */
/*     NOTE   :                                                         * */
/*           VALID ONLY FOR BARYON-BARYON-DISTANCES LESS THAN 1.32 FM   * */
/*           (1.32 = 2 * HARD-CORE-RADIUS [HRC] )                       * */
/*     QUANTITIES:                                                 * */
/*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME* */
/*           SRT      - SQRT OF S                                       * */
/*           NSTAR =1 INCLUDING N* RESORANCE,ELSE NOT                   * */
/*           NDIRCT=1 INCLUDING DIRECT PION PRODUCTION PROCESS         * */
/*           IBLOCK   - THE INFORMATION BACK                            * */
/*                      0-> COLLISION CANNOT HAPPEN                     * */
/*                      1-> N-N ELASTIC COLLISION                       * */
/*                      2-> N+N->N+DELTA,OR N+N->N+N* REACTION          * */
/*                      3-> N+DELTA->N+N OR N+N*->N+N REACTION          * */
/*                      4-> N+N->N+N+PION,DIRTCT PROCESS                * */
/*                     5-> DELTA(N*)+DELTA(N*)   TOTAL   COLLISIONS    * */
/*           N12       - IS USED TO SPECIFY BARYON-BARYON REACTION      * */
/*                      CHANNELS. M12 IS THE REVERSAL CHANNEL OF N12    * */
/*                      N12,                                            * */
/*                      M12=1 FOR p+n-->delta(+)+ n                     * */
/*                          2     p+n-->delta(0)+ p                     * */
/*                          3     p+p-->delta(++)+n                     * */
/*                          4     p+p-->delta(+)+p                      * */
/*                          5     n+n-->delta(0)+n                      * */
/*                          6     n+n-->delta(-)+p                      * */
/*                          7     n+p-->N*(0)(1440)+p                   * */
/*                          8     n+p-->N*(+)(1440)+n                   * */
/*                        9     p+p-->N*(+)(1535)+p                     * */
/*                        10    n+n-->N*(0)(1535)+n                     * */
/*                         11    n+p-->N*(+)(1535)+n                     * */
/*                        12    n+p-->N*(0)(1535)+p */
/*                        13    D(++)+D(-)-->N*(+)(1440)+n */
/*                         14    D(++)+D(-)-->N*(0)(1440)+p */
/*                        15    D(+)+D(0)--->N*(+)(1440)+n */
/*                        16    D(+)+D(0)--->N*(0)(1440)+p */
/*                        17    D(++)+D(0)-->N*(+)(1535)+p */
/*                        18    D(++)+D(-)-->N*(0)(1535)+p */
/*                        19    D(++)+D(-)-->N*(+)(1535)+n */
/*                        20    D(+)+D(+)-->N*(+)(1535)+p */
/*                        21    D(+)+D(0)-->N*(+)(1535)+n */
/*                        22    D(+)+D(0)-->N*(0)(1535)+p */
/*                        23    D(+)+D(-)-->N*(0)(1535)+n */
/*                        24    D(0)+D(0)-->N*(0)(1535)+n */
/*                          25    N*(+)(14)+N*(+)(14)-->N*(+)(15)+p */
/*                          26    N*(0)(14)+N*(0)(14)-->N*(0)(15)+n */
/*                          27    N*(+)(14)+N*(0)(14)-->N*(+)(15)+n */
/*                        28    N*(+)(14)+N*(0)(14)-->N*(0)(15)+p */
/*                        29    N*(+)(14)+D+-->N*(+)(15)+p */
/*                        30    N*(+)(14)+D0-->N*(+)(15)+n */
/*                        31    N*(+)(14)+D(-)-->N*(0)(1535)+n */
/*                        32    N*(0)(14)+D++--->N*(+)(15)+p */
/*                        33    N*(0)(14)+D+--->N*(+)(15)+n */
/*                        34    N*(0)(14)+D+--->N*(0)(15)+p */
/*                        35    N*(0)(14)+D0-->N*(0)(15)+n */
/*                        36    N*(+)(14)+D0--->N*(0)(15)+p */
/*                        +++ */
/*               AND MORE CHANNELS AS LISTED IN THE NOTE BOOK */

/* NOTE ABOUT N*(1440) RESORANCE:                                       * */
/*     As it has been discussed in VerWest's paper,I= 1 (initial isospin) */
/*     channel can all be attributed to delta resorance while I= 0      * */
/*     channel can all be  attribured to N* resorance.Only in n+p       * */
/*     one can have I=0 channel so is the N*(1440) resorance            * */
/* REFERENCES:    J. CUGNON ET AL., NUCL. PHYS. A352, 505 (1981)        * */
/*                    Y. KITAZOE ET AL., PHYS. LETT. 166B, 35 (1986)    * */
/*                    B. VerWest el al., PHYS. PRV. C25 (1982)1979      * */
/*                    Gy. Wolf  et al, Nucl Phys A517 (1990) 615        * */
/*                    CUTOFF = 2 * AVMASS + 20 MEV                      * */
/*                                                                      * */
/*       for N*(1535) we use the parameterization by Gy. Wolf et al     * */
/*       Nucl phys A552 (1993) 349, added May 18, 1994                  * */
/* ********************************* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /ff/ */
/* c      SAVE /gg/ */
/* c      SAVE /INPUT/ */
/* c      SAVE /NN/ */
/* c      SAVE /BG/ */
/* c      SAVE /RUN/ */
/* c      SAVE /PA/ */
/* c      SAVE /PB/ */
/* c      SAVE /PC/ */
/* c      SAVE /PD/ */
/* c      SAVE /input1/ */
/* ----------------------------------------------------------------------- */
    *xinel = 0.f;
    *sigk = 0.f;
    *xsk1 = 0.f;
    *xsk2 = 0.f;
    *xsk3 = 0.f;
    *xsk4 = 0.f;
    *xsk5 = 0.f;
    em1 = cc_1.e[*i1 - 1];
    em2 = cc_1.e[*i2 - 1];
/* Computing 2nd power */
    r__1 = *px;
/* Computing 2nd power */
    r__2 = *py;
/* Computing 2nd power */
    r__3 = *pz;
    pr = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
/*     IF THERE WERE 2 N*(1535) AND THEY DIDN'T SCATT. ELAST., */
/*     ALLOW THEM TO PRODUCE KAONS. NO OTHER INELASTIC CHANNELS */
/*     ARE KNOWN */
/*       if((lb(i1).ge.12).and.(lb(i2).ge.12))return */
/*     ALL the inelastic collisions between N*(1535) and Delta as well */
/*     as N*(1440) TO PRODUCE KAONS, NO OTHER CHANNELS ARE KNOWN */
/*       if((lb(i1).ge.12).and.(lb(i2).ge.3))return */
/*       if((lb(i2).ge.12).and.(lb(i1).ge.3))return */
/*     calculate the N*(1535) production cross section in I1+I2 collisions */
    i__3 = (i__1 = ee_1.lb[*i1 - 1], abs(i__1));
    i__4 = (i__2 = ee_1.lb[*i2 - 1], abs(i__2));
    n1535_(&i__3, &i__4, srt, &x1535);

/* for Delta+Delta-->N*(1440 OR 1535)+N AND N*(1440)+N*(1440)-->N*(1535)+X */
/*     AND DELTA+N*(1440)-->N*(1535)+X */
/* WE ASSUME THEY HAVE THE SAME CROSS SECTIONS as CORRESPONDING N+N COLLISION): */
/* FOR D++D0, D+D+,D+D-,D0D0,N*+N*+,N*0N*0,N*(+)D+,N*(+)D(-),N*(0)D(0) */
/* N*(1535) production, kaon production and reabsorption through */
/* D(N*)+D(N*)-->NN are ALLOWED. */
/* CROSS SECTION FOR KAON PRODUCTION from the four channels are */
/* for NLK channel */
    akp = .498f;
    ak0 = .498f;
    ana = .94f;
    ada = 1.232f;
    al = 1.1157f;
    as = 1.1197f;
    *xsk1 = 0.f;
    *xsk2 = 0.f;
    *xsk3 = 0.f;
    *xsk4 = 0.f;
    t1nlk = ana + al + akp;
    if (*srt <= t1nlk) {
	goto L222;
    }
    *xsk1 = pplpk_(srt) * 1.5f;
/* for DLK channel */
    t1dlk = ada + al + akp;
    t2dlk = ada + al - akp;
    if (*srt <= t1dlk) {
	goto L222;
    }
    es = *srt;
/* Computing 2nd power */
    r__1 = es;
/* Computing 2nd power */
    r__2 = t1dlk;
/* Computing 2nd power */
    r__3 = es;
/* Computing 2nd power */
    r__4 = t2dlk;
/* Computing 2nd power */
    r__5 = es;
    pmdlk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (
	    r__5 * r__5 * 4.f);
    pmdlk = sqrt(pmdlk2);
    *xsk3 = pplpk_(srt) * 1.5f;
/* for NSK channel */
    t1nsk = ana + as + akp;
    t2nsk = ana + as - akp;
    if (*srt <= t1nsk) {
	goto L222;
    }
/* Computing 2nd power */
    r__1 = es;
/* Computing 2nd power */
    r__2 = t1nsk;
/* Computing 2nd power */
    r__3 = es;
/* Computing 2nd power */
    r__4 = t2nsk;
/* Computing 2nd power */
    r__5 = es;
    pmnsk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (
	    r__5 * r__5 * 4.f);
    pmnsk = sqrt(pmnsk2);
    *xsk2 = (ppk1_(srt) + ppk0_(srt)) * 1.5f;
/* for DSK channel */
    t1dsk = ada + as + akp;
    t2dsk = ada + as - akp;
    if (*srt <= t1dsk) {
	goto L222;
    }
/* Computing 2nd power */
    r__1 = es;
/* Computing 2nd power */
    r__2 = t1dsk;
/* Computing 2nd power */
    r__3 = es;
/* Computing 2nd power */
    r__4 = t2dsk;
/* Computing 2nd power */
    r__5 = es;
    pmdsk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (
	    r__5 * r__5 * 4.f);
    pmdsk = sqrt(pmdsk2);
    *xsk4 = (ppk1_(srt) + ppk0_(srt)) * 1.5f;
/* sp11/21/01 */
/* phi production */
    if (*srt <= 2.898914f) {
	goto L222;
    }
/*  !! mb put the correct form */
    *xsk5 = 1e-4f;
/* sp11/21/01 end */
/* THE TOTAL KAON+ PRODUCTION CROSS SECTION IS THEN */
L222:
    *sigk = *xsk1 + *xsk2 + *xsk3 + *xsk4;
/* bz3/7/99 neutralk */
    *xsk1 *= 2.f;
    *xsk2 *= 2.f;
    *xsk3 *= 2.f;
    *xsk4 *= 2.f;
    *sigk = *sigk * 2.f + *xsk5;
/* bz3/7/99 neutralk end */
    idd = (i__1 = ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1], abs(i__1));
/* The reabsorption cross section for the process */
/* D(N*)D(N*)-->NN is */
    s2d = reab2d_(i1, i2, srt);
/* bz3/16/99 pion */
    s2d = 0.f;
/* bz3/16/99 pion end */
/* (1) N*(1535)+D(N*(1440)) reactions */
/*    we allow kaon production and reabsorption only */
    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) >= 12 && (i__2 = ee_1.lb[*i2 - 1]
	    , abs(i__2)) >= 12 || (i__3 = ee_1.lb[*i1 - 1], abs(i__3)) >= 12 
	    && (i__4 = ee_1.lb[*i2 - 1], abs(i__4)) >= 6 || (i__5 = ee_1.lb[*
	    i2 - 1], abs(i__5)) >= 12 && (i__6 = ee_1.lb[*i1 - 1], abs(i__6)) 
	    >= 6) {
	*xinel = *sigk + s2d;
	return 0;
    }
/* channels have the same charge as pp */
    if (idd == 63 || idd == 64 || idd == 48 || idd == 49 || idd == 121 || idd 
	    == 100 || idd == 88 || idd == 66 || idd == 90 || idd == 70) {
	*xinel = x1535 + *sigk + s2d;
	return 0;
    }
/* IN DELTA+N*(1440) and N*(1440)+N*(1440) COLLISIONS, */
/* N*(1535), kaon production and reabsorption are ALLOWED */
/* IN N*(1440)+N*(1440) COLLISIONS, ONLY N*(1535) IS ALLOWED */
    if (idd == 110 || idd == 77 || idd == 80) {
	*xinel = x1535 + *sigk + s2d;
	return 0;
    }
    if (idd == 54 || idd == 56) {
/* LIKE FOR N+P COLLISION, */
/* IN DELTA+DELTA COLLISIONS BOTH N*(1440) AND N*(1535) CAN BE PRODUCED */
	sig2 = sigma_(srt, &c__2, &c__0, &c__1) * .75f;
	*xinel = (sig2 + x1535) * 2.f + *sigk + s2d;
	return 0;
    }
    return 0;
} /* xddin_ */

/* ***************************************** */
doublereal dirct1_(real *srt)
{
    /* Initialized data */

    static real earray[122] = { 1.5683f,1.5783f,1.5883f,1.5983f,1.6083f,
	    1.6183f,1.6283f,1.6383f,1.6483f,1.6583f,1.6683f,1.6783f,1.6883f,
	    1.6983f,1.7083f,1.7183f,1.7283f,1.7383f,1.7483f,1.7583f,1.7683f,
	    1.7783f,1.7883f,1.7983f,1.8083f,1.8183f,1.8283f,1.8383f,1.8483f,
	    1.8583f,1.8683f,1.8783f,1.8883f,1.8983f,1.9083f,1.9183f,1.9283f,
	    1.9383f,1.9483f,1.9583f,1.9683f,1.9783f,1.9883f,1.9983f,2.0083f,
	    2.0183f,2.0283f,2.0383f,2.0483f,2.0583f,2.0683f,2.0783f,2.0883f,
	    2.0983f,2.1083f,2.1183f,2.1283f,2.1383f,2.1483f,2.1583f,2.1683f,
	    2.1783f,2.1883f,2.1983f,2.2083f,2.2183f,2.2283f,2.2383f,2.2483f,
	    2.2583f,2.2683f,2.2783f,2.2883f,2.2983f,2.3083f,2.3183f,2.3283f,
	    2.3383f,2.3483f,2.3583f,2.3683f,2.3783f,2.3883f,2.3983f,2.4083f,
	    2.4183f,2.4283f,2.4383f,2.4483f,2.4583f,2.4683f,2.4783f,2.4883f,
	    2.4983f,2.5083f,2.5183f,2.5283f,2.5383f,2.5483f,2.5583f,2.5683f,
	    2.5783f,2.5883f,2.5983f,2.6083f,2.6183f,2.6283f,2.6383f,2.6483f,
	    2.6583f,2.6683f,2.6783f,2.6883f,2.6983f,2.7083f,2.7183f,2.7283f,
	    2.7383f,2.7483f,2.7583f,2.7683f,2.7783f };
    static real xarray[122] = { .017764091f,.5643668f,.8150568f,1.045565f,
	    2.133695f,3.327922f,4.206488f,3.471242f,4.486876f,5.542213f,
	    6.800052f,7.192446f,6.829848f,6.580306f,6.86841f,8.527946f,
	    10.1572f,9.716511f,9.298335f,8.90131f,10.31213f,10.52185f,
	    11.1763f,11.61639f,12.05577f,12.71596f,13.46036f,14.2206f,
	    14.65449f,14.94775f,14.9331f,15.32907f,16.56481f,16.29422f,
	    15.18548f,14.12658f,13.72544f,13.24488f,13.31003f,14.4268f,
	    12.84423f,12.49025f,12.14858f,11.8187f,11.18993f,11.35816f,
	    11.09447f,10.83873f,10.61592f,10.53754f,9.425521f,8.195912f,
	    9.661075f,9.696192f,9.200142f,8.953734f,8.715461f,8.484999f,
	    8.320765f,8.255512f,8.190969f,8.127125f,8.079508f,8.073004f,
	    8.010611f,7.948909f,7.887895f,7.761005f,7.62629f,7.494696f,
	    7.366132f,7.530178f,8.392097f,9.046881f,8.962544f,8.879403f,
	    8.797427f,8.716601f,8.636904f,8.558312f,8.404368f,8.328978f,
	    8.254617f,8.181265f,8.108907f,8.037527f,7.9671f,7.897617f,
	    7.829057f,7.761405f,7.694647f,7.628764f,7.563742f,7.49957f,
	    7.387562f,7.273281f,7.161334f,6.973375f,6.529592f,6.280323f,
	    6.293136f,6.305725f,6.318097f,6.330258f,6.342214f,6.353968f,
	    6.365528f,6.376895f,6.388079f,6.399081f,6.409906f,6.42056f,
	    6.431045f,6.441367f,6.451529f,6.461533f,6.471386f,6.481091f,
	    6.49065f,6.476413f,6.297259f,6.097826f };

    /* System generated locals */
    real ret_val;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static integer ie;
    static real xmin, ymin, xmax, ymax;

/*  This function contains the experimental, direct pion(+) + p cross sections * */
/*  srt    = DSQRT(s) in GeV                                                   * */
/*  dirct1  = cross section in fm**2                                     * */
/*  earray = EXPerimental table with the srt */
/*  xarray = EXPerimental table with cross sections in mb (curve to guide eye) * */
/* ***************************************** */
/*      real*4   xarray(122), earray(122) */
    if (*srt < earray[0]) {
	ret_val = 1e-5f;
	return ret_val;
    }
    if (*srt > earray[121]) {
	ret_val = xarray[121];
	ret_val /= 10.f;
	return ret_val;
    }

/* Interpolate double logarithmically to find xdirct2(srt) */

    for (ie = 1; ie <= 122; ++ie) {
	if (earray[ie - 1] == *srt) {
	    ret_val = xarray[ie - 1];
	    ret_val /= 10.f;
	    return ret_val;
	}
	if (earray[ie - 1] > *srt) {
	    ymin = log(xarray[ie - 2]);
	    ymax = log(xarray[ie - 1]);
	    xmin = log(earray[ie - 2]);
	    xmax = log(earray[ie - 1]);
	    ret_val = exp(ymin + (log(*srt) - xmin) * (ymax - ymin) / (xmax - 
		    xmin));
	    ret_val /= 10.f;
	    goto L50;
	}
/* L1001: */
    }
L50:
    return ret_val;
} /* dirct1_ */

/* ****************************** */
/* ***************************************** */
doublereal dirct2_(real *srt)
{
    /* Initialized data */

    static real earray[122] = { 1.5683f,1.5783f,1.5883f,1.5983f,1.6083f,
	    1.6183f,1.6283f,1.6383f,1.6483f,1.6583f,1.6683f,1.6783f,1.6883f,
	    1.6983f,1.7083f,1.7183f,1.7283f,1.7383f,1.7483f,1.7583f,1.7683f,
	    1.7783f,1.7883f,1.7983f,1.8083f,1.8183f,1.8283f,1.8383f,1.8483f,
	    1.8583f,1.8683f,1.8783f,1.8883f,1.8983f,1.9083f,1.9183f,1.9283f,
	    1.9383f,1.9483f,1.9583f,1.9683f,1.9783f,1.9883f,1.9983f,2.0083f,
	    2.0183f,2.0283f,2.0383f,2.0483f,2.0583f,2.0683f,2.0783f,2.0883f,
	    2.0983f,2.1083f,2.1183f,2.1283f,2.1383f,2.1483f,2.1583f,2.1683f,
	    2.1783f,2.1883f,2.1983f,2.2083f,2.2183f,2.2283f,2.2383f,2.2483f,
	    2.2583f,2.2683f,2.2783f,2.2883f,2.2983f,2.3083f,2.3183f,2.3283f,
	    2.3383f,2.3483f,2.3583f,2.3683f,2.3783f,2.3883f,2.3983f,2.4083f,
	    2.4183f,2.4283f,2.4383f,2.4483f,2.4583f,2.4683f,2.4783f,2.4883f,
	    2.4983f,2.5083f,2.5183f,2.5283f,2.5383f,2.5483f,2.5583f,2.5683f,
	    2.5783f,2.5883f,2.5983f,2.6083f,2.6183f,2.6283f,2.6383f,2.6483f,
	    2.6583f,2.6683f,2.6783f,2.6883f,2.6983f,2.7083f,2.7183f,2.7283f,
	    2.7383f,2.7483f,2.7583f,2.7683f,2.7783f };
    static real xarray[122] = { .5773182f,1.404156f,2.578629f,3.832013f,
	    4.906011f,9.076963f,13.10492f,10.65975f,15.31156f,19.77611f,
	    19.92874f,18.68979f,19.80114f,18.39536f,14.34269f,13.35353f,
	    13.58822f,14.57031f,10.24686f,11.23386f,9.764803f,10.35652f,
	    10.53539f,10.07524f,9.582198f,9.596469f,9.818489f,9.012848f,
	    9.378012f,9.529244f,9.529698f,8.835624f,6.671396f,8.797758f,
	    8.133437f,7.866227f,7.823946f,7.808504f,7.791755f,7.502062f,
	    7.417275f,7.592349f,7.752028f,7.910585f,8.068122f,8.224736f,
	    8.075289f,7.895902f,7.721359f,7.551512f,7.386224f,7.225343f,
	    7.068739f,6.916284f,6.767842f,6.623294f,6.48252f,6.345404f,
	    6.211833f,7.33951f,7.531462f,7.724824f,7.91962f,7.848021f,
	    7.639856f,7.571083f,7.508881f,7.447474f,7.386855f,7.327011f,
	    7.164454f,7.001266f,6.842526f,6.688094f,6.537823f,6.391583f,
	    6.249249f,6.110689f,5.97579f,5.8942f,5.959503f,6.024602f,
	    6.089505f,6.154224f,6.21876f,6.283128f,6.347331f,6.297411f,
	    6.120248f,5.948606f,6.494864f,6.357106f,6.222824f,6.09191f,
	    5.964267f,5.839795f,5.718402f,5.599994f,5.499146f,5.451325f,
	    5.404156f,5.357625f,5.311721f,5.266435f,5.301964f,5.343963f,
	    5.385833f,5.427577f,5.4692f,5.510702f,5.552088f,5.593359f,
	    5.63452f,5.67557f,5.716515f,5.757356f,5.798093f,5.838732f,
	    5.879272f,5.919717f,5.960068f,5.980941f };

    /* System generated locals */
    real ret_val;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static integer ie;
    static real xmin, ymin, xmax, ymax;

/*  This function contains the experimental, direct pion(-) + p cross sections * */
/*  srt    = DSQRT(s) in GeV                                                   * */
/*  dirct2 = cross section in fm**2 */
/*  earray = EXPerimental table with the srt */
/*  xarray = EXPerimental table with cross sections in mb (curve to guide eye) * */
/* ***************************************** */
/*      real*4   xarray(122), earray(122) */
    if (*srt < earray[0]) {
	ret_val = 1e-5f;
	return ret_val;
    }
    if (*srt > earray[121]) {
	ret_val = xarray[121];
	ret_val /= 10.f;
	return ret_val;
    }

/* Interpolate double logarithmically to find xdirct2(srt) */

    for (ie = 1; ie <= 122; ++ie) {
	if (earray[ie - 1] == *srt) {
	    ret_val = xarray[ie - 1];
	    ret_val /= 10.f;
	    return ret_val;
	}
	if (earray[ie - 1] > *srt) {
	    ymin = log(xarray[ie - 2]);
	    ymax = log(xarray[ie - 1]);
	    xmin = log(earray[ie - 2]);
	    xmax = log(earray[ie - 1]);
	    ret_val = exp(ymin + (log(*srt) - xmin) * (ymax - ymin) / (xmax - 
		    xmin));
	    ret_val /= 10.f;
	    goto L50;
	}
/* L1001: */
    }
L50:
    return ret_val;
} /* dirct2_ */

/* ****************************** */
/* ***************************** */
/* this program calculates the elastic cross section for rho+nucleon */
/* through higher resonances */
/*       real*4 function ErhoN(em1,em2,lb1,lb2,srt) */
doublereal erhon_(real *em1, real *em2, integer *lb1, integer *lb2, real *srt)
{
    /* Initialized data */

    static real arrayj[19] = { .5f,1.5f,.5f,.5f,2.5f,2.5f,1.5f,.5f,1.5f,3.5f,
	    1.5f,.5f,1.5f,.5f,2.5f,.5f,1.5f,2.5f,3.5f };
    static real arrayl[19] = { 1.f,2.f,0.f,0.f,2.f,3.f,2.f,1.f,1.f,3.f,1.f,
	    0.f,2.f,0.f,3.f,1.f,1.f,2.f,3.f };
    static real arraym[19] = { 1.44f,1.52f,1.535f,1.65f,1.675f,1.68f,1.7f,
	    1.71f,1.72f,1.99f,1.6f,1.62f,1.7f,1.9f,1.905f,1.91f,1.86f,1.93f,
	    1.95f };
    static real arrayw[19] = { .2f,.125f,.15f,.15f,.155f,.125f,.1f,.11f,.2f,
	    .29f,.25f,.16f,.28f,.15f,.3f,.22f,.25f,.25f,.24f };
    static real arrayb[19] = { .15f,.2f,.05f,.175f,.025f,.125f,.1f,.2f,.53f,
	    .34f,.05f,.07f,.15f,.45f,.45f,.058f,.08f,.12f,.08f };

    /* System generated locals */
    integer i__1, i__2;
    real ret_val;

    /* Local variables */
    static real pi;
    static integer ir;
    static real xs, xs0;
    extern doublereal fdr_(real *, real *, real *, real *, real *, real *, 
	    real *, real *);
    static real branch;

/* date : Dec. 19, 1994 */
/* **************************** */
/*       implicit real*4 (a-h,o-z) */
/* the minimum energy for pion+delta collision */
    pi = 3.1415926f;
    xs = 0.f;
/* include contribution from each resonance */
    for (ir = 1; ir <= 19; ++ir) {
/* bz11/25/98 */
	if (ir <= 8) {
/*       if(lb1*lb2.eq.27.OR.LB1*LB2.EQ.25*2)branch=0. */
/*       if(lb1*lb2.eq.26.OR.LB1*LB2.EQ.26*2)branch=1./3. */
/*       if(lb1*lb2.eq.27*2.OR.LB1*LB2.EQ.25)branch=2./3. */
/*       ELSE */
/*       if(lb1*lb2.eq.27.OR.LB1*LB2.EQ.25*2)branch=1. */
/*       if(lb1*lb2.eq.26.OR.LB1*LB2.EQ.26*2)branch=2./3. */
/*       if(lb1*lb2.eq.27*2.OR.LB1*LB2.EQ.25)branch=1./3. */
/*       ENDIF */
	    if (*lb1 * *lb2 == 27 && (*lb1 == 1 || *lb2 == 1) || *lb1 * *lb2 
		    == 50 && (*lb1 == 2 || *lb2 == 2) || (*lb1 * *lb2 == -25 
		    && (*lb1 == -1 || *lb2 == -1) || *lb1 * *lb2 == -54 && (*
		    lb1 == -2 || *lb2 == -2))) {
		branch = 0.f;
	    }
	    if ((i__1 = *lb1 * *lb2, abs(i__1)) == 26 && (abs(*lb1) == 1 || 
		    abs(*lb2) == 1) || (i__2 = *lb1 * *lb2, abs(i__2)) == 52 
		    && (abs(*lb1) == 2 || abs(*lb2) == 2)) {
		branch = .33333333333333331f;
	    }
	    if (*lb1 * *lb2 == 54 && (*lb1 == 2 || *lb2 == 2) || *lb1 * *lb2 
		    == 25 && (*lb1 == 1 || *lb2 == 1) || (*lb1 * *lb2 == -50 
		    && (*lb1 == -2 || *lb2 == -2) || *lb1 * *lb2 == -27 && (*
		    lb1 == -1 || *lb2 == -1))) {
		branch = .66666666666666663f;
	    }
	} else {
	    if (*lb1 * *lb2 == 27 && (*lb1 == 1 || *lb2 == 1) || *lb1 * *lb2 
		    == 50 && (*lb1 == 2 || *lb2 == 2) || (*lb1 * *lb2 == -25 
		    && (*lb1 == -1 || *lb2 == -1) || *lb1 * *lb2 == -54 && (*
		    lb1 == -2 || *lb2 == -2))) {
		branch = 1.f;
	    }
	    if ((i__1 = *lb1 * *lb2, abs(i__1)) == 26 && (abs(*lb1) == 1 || 
		    abs(*lb2) == 1) || (i__2 = *lb1 * *lb2, abs(i__2)) == 52 
		    && (abs(*lb1) == 2 || abs(*lb2) == 2)) {
		branch = .66666666666666663f;
	    }
	    if (*lb1 * *lb2 == 54 && (*lb1 == 2 || *lb2 == 2) || *lb1 * *lb2 
		    == 25 && (*lb1 == 1 || *lb2 == 1) || (*lb1 * *lb2 == -50 
		    && (*lb1 == -2 || *lb2 == -2) || *lb1 * *lb2 == -27 && (*
		    lb1 == -1 || *lb2 == -1))) {
		branch = .33333333333333331f;
	    }
	}
/* bz11/25/98end */
	xs0 = fdr_(&arraym[ir - 1], &arrayj[ir - 1], &arrayl[ir - 1], &arrayw[
		ir - 1], &arrayb[ir - 1], srt, em1, em2);
	xs += pi * 1.3f * branch * xs0 * .038927290000000003f;
/* L1001: */
    }
    ret_val = xs;
    return ret_val;
} /* erhon_ */

/* **************************8 */
/* FUNCTION FDE(DMASS) GIVES DELTA MASS DISTRIBUTION BY USING OF */
/* KITAZOE'S FORMULA */
/*        REAL*4 FUNCTION FDR(DMASS,aj,al,width,widb0,srt,em1,em2) */
doublereal fdr_(real *dmass, real *aj, real *al, real *width, real *widb0, 
	real *srt, real *em1, real *em2)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4, r__5;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static real b, q, q0, ak2, ak02, amd, amp;

    amd = *em1;
    amp = *em2;
/* Computing 2nd power */
    r__2 = *dmass;
/* Computing 2nd power */
    r__3 = amd;
/* Computing 2nd power */
    r__4 = amp;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = amp * amd;
    ak02 = r__1 * r__1 * .25f - r__5 * r__5;
    if (ak02 > 0.f) {
	q0 = sqrt(ak02 / *dmass);
    } else {
	q0 = 0.f;
	ret_val = 0.f;
	return ret_val;
    }
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = amd;
/* Computing 2nd power */
    r__4 = amp;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = amp * amd;
    ak2 = r__1 * r__1 * .25f - r__5 * r__5;
    if (ak2 > 0.f) {
	q = sqrt(ak2 / *dmass);
    } else {
	q = 0.f;
	ret_val = 0.f;
	return ret_val;
    }
    d__1 = (doublereal) (q / q0);
    d__2 = (doublereal) (*al * 2.f + 1);
    d__3 = (doublereal) (q / q0);
    d__4 = (doublereal) (*al * 2);
    b = *widb0 * 1.2f * *dmass / *srt * pow_dd(&d__1, &d__2) / (pow_dd(&d__3, 
	    &d__4) * .2f + 1.f);
/* Computing 2nd power */
    r__1 = *width;
/* Computing 2nd power */
    r__2 = *srt - *dmass;
/* Computing 2nd power */
    r__3 = *width;
/* Computing 2nd power */
    r__4 = q;
    ret_val = (*aj * 2.f + 1) * (r__1 * r__1) * b / (r__2 * r__2 + r__3 * 
	    r__3 * .25f) / (r__4 * r__4 * 6.f);
    return ret_val;
} /* fdr_ */

/* ***************************** */
/* this program calculates the elastic cross section for pion+delta */
/* through higher resonances */
/*       REAL*4 FUNCTION DIRCT3(SRT) */
doublereal dirct3_(real *srt)
{
    /* Initialized data */

    static real arrayj[17] = { 1.5f,.5f,2.5f,2.5f,1.5f,.5f,1.5f,3.5f,1.5f,.5f,
	    1.5f,.5f,2.5f,.5f,1.5f,2.5f,3.5f };
    static real arrayl[17] = { 2.f,0.f,2.f,3.f,2.f,1.f,1.f,3.f,1.f,0.f,2.f,
	    0.f,3.f,1.f,1.f,2.f,3.f };
    static real arraym[17] = { 1.52f,1.65f,1.675f,1.68f,1.7f,1.71f,1.72f,
	    1.99f,1.6f,1.62f,1.7f,1.9f,1.905f,1.91f,1.86f,1.93f,1.95f };
    static real arrayw[17] = { .125f,.15f,.155f,.125f,.1f,.11f,.2f,.29f,.25f,
	    .16f,.28f,.15f,.3f,.22f,.25f,.25f,.24f };
    static real arrayb[17] = { .55f,.6f,.375f,.6f,.1f,.15f,.15f,.05f,.35f,.3f,
	    .15f,.1f,.1f,.22f,.2f,.09f,.4f };

    /* System generated locals */
    real ret_val;

    /* Local variables */
    static real pi;
    static integer ir;
    static real xs;
    extern doublereal fd1_(real *, real *, real *, real *, real *, real *);
    static real xs0, amn, amp, branch;

/* date : Dec. 19, 1994 */
/* **************************** */
/*     implicit real*4 (a-h,o-z) */
/* the minimum energy for pion+delta collision */
    pi = 3.1415926f;
    amn = .938f;
    amp = .138f;
    xs = 0.f;
/* include contribution from each resonance */
    branch = .33333333333333331f;
    for (ir = 1; ir <= 17; ++ir) {
	if (ir > 8) {
	    branch = .66666666666666663f;
	}
	xs0 = fd1_(&arraym[ir - 1], &arrayj[ir - 1], &arrayl[ir - 1], &arrayw[
		ir - 1], &arrayb[ir - 1], srt);
	xs += pi * 1.3f * branch * xs0 * .038927290000000003f;
/* L1001: */
    }
    ret_val = xs;
    return ret_val;
} /* dirct3_ */

/* **************************8 */
/* FUNCTION FDE(DMASS) GIVES DELTA MASS DISTRIBUTION BY USING OF */
/* KITAZOE'S FORMULA */
/*        REAL*4 FUNCTION FD1(DMASS,aj,al,width,widb0,srt) */
doublereal fd1_(real *dmass, real *aj, real *al, real *width, real *widb0, 
	real *srt)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4, r__5;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static real b, q, q0, ak2, ak02, amd, amn, amp;

    amn = .938f;
    amp = .138f;
    amd = amn;
/* Computing 2nd power */
    r__2 = *dmass;
/* Computing 2nd power */
    r__3 = amd;
/* Computing 2nd power */
    r__4 = amp;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = amp * amd;
    ak02 = r__1 * r__1 * .25f - r__5 * r__5;
    if (ak02 > 0.f) {
	q0 = sqrt(ak02 / *dmass);
    } else {
	q0 = 0.f;
	ret_val = 0.f;
	return ret_val;
    }
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = amd;
/* Computing 2nd power */
    r__4 = amp;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = amp * amd;
    ak2 = r__1 * r__1 * .25f - r__5 * r__5;
    if (ak2 > 0.f) {
	q = sqrt(ak2 / *dmass);
    } else {
	q = 0.f;
	ret_val = 0.f;
	return ret_val;
    }
    d__1 = (doublereal) (q / q0);
    d__2 = (doublereal) (*al * 2.f + 1);
    d__3 = (doublereal) (q / q0);
    d__4 = (doublereal) (*al * 2);
    b = *widb0 * 1.2f * *dmass / *srt * pow_dd(&d__1, &d__2) / (pow_dd(&d__3, 
	    &d__4) * .2f + 1.f);
/* Computing 2nd power */
    r__1 = *width;
/* Computing 2nd power */
    r__2 = *srt - *dmass;
/* Computing 2nd power */
    r__3 = *width;
/* Computing 2nd power */
    r__4 = q;
    ret_val = (*aj * 2.f + 1) * (r__1 * r__1) * b / (r__2 * r__2 + r__3 * 
	    r__3 * .25f) / (r__4 * r__4 * 2.f);
    return ret_val;
} /* fd1_ */

/* ***************************** */
/* this program calculates the elastic cross section for pion+delta */
/* through higher resonances */
/*       REAL*4 FUNCTION DPION(EM1,EM2,LB1,LB2,SRT) */
doublereal dpion_(real *em1, real *em2, integer *lb1, integer *lb2, real *srt)
{
    /* Initialized data */

    static real arrayj[19] = { .5f,1.5f,.5f,.5f,2.5f,2.5f,1.5f,.5f,1.5f,3.5f,
	    1.5f,.5f,1.5f,.5f,2.5f,.5f,1.5f,2.5f,3.5f };
    static real arrayl[19] = { 1.f,2.f,0.f,0.f,2.f,3.f,2.f,1.f,1.f,3.f,1.f,
	    0.f,2.f,0.f,3.f,1.f,1.f,2.f,3.f };
    static real arraym[19] = { 1.44f,1.52f,1.535f,1.65f,1.675f,1.68f,1.7f,
	    1.71f,1.72f,1.99f,1.6f,1.62f,1.7f,1.9f,1.905f,1.91f,1.86f,1.93f,
	    1.95f };
    static real arrayw[19] = { .2f,.125f,.15f,.15f,.155f,.125f,.1f,.11f,.2f,
	    .29f,.25f,.16f,.28f,.15f,.3f,.22f,.25f,.25f,.24f };
    static real arrayb[19] = { .15f,.25f,0.f,.05f,.575f,.125f,.379f,.1f,.1f,
	    .062f,.45f,.6f,.6984f,.05f,.25f,.089f,.19f,.2f,.13f };

    /* System generated locals */
    integer i__1, i__2;
    real ret_val;

    /* Local variables */
    static real pi;
    static integer ir;
    static real xs;
    extern doublereal fd2_(real *, real *, real *, real *, real *, real *, 
	    real *, real *);
    static real xs0, amn, amp, branch;

/* date : Dec. 19, 1994 */
/* **************************** */
/*     implicit real*4 (a-h,o-z) */
/* the minimum energy for pion+delta collision */
    pi = 3.1415926f;
    amn = .94f;
    amp = .14f;
    xs = 0.f;
/* include contribution from each resonance */
    for (ir = 1; ir <= 19; ++ir) {
	branch = 0.f;
/* bz11/25/98 */
	if (ir <= 8) {
/*       IF(LB1*LB2.EQ.5*7.OR.LB1*LB2.EQ.3*8)branch=1./6. */
/*       IF(LB1*LB2.EQ.4*7.OR.LB1*LB2.EQ.4*8)branch=1./3. */
/*       IF(LB1*LB2.EQ.5*6.OR.LB1*LB2.EQ.3*9)branch=1./2. */
/*       ELSE */
/*       IF(LB1*LB2.EQ.5*8.OR.LB1*LB2.EQ.5*6)branch=2./5. */
/*       IF(LB1*LB2.EQ.3*9.OR.LB1*LB2.EQ.3*7)branch=2./5. */
/*       IF(LB1*LB2.EQ.5*7.OR.LB1*LB2.EQ.3*8)branch=8./15. */
/*       IF(LB1*LB2.EQ.4*7.OR.LB1*LB2.EQ.4*8)branch=1./15. */
/*       IF(LB1*LB2.EQ.4*9.OR.LB1*LB2.EQ.4*6)branch=3./5. */
/*       ENDIF */
	    if (*lb1 * *lb2 == 35 && (*lb1 == 5 || *lb2 == 5) || *lb1 * *lb2 
		    == 24 && (*lb1 == 3 || *lb2 == 3) || (*lb1 * *lb2 == -21 
		    && (*lb1 == 3 || *lb2 == 3) || *lb1 * *lb2 == -40 && (*
		    lb1 == 5 || *lb2 == 5))) {
		branch = .16666666666666666f;
	    }
	    if ((i__1 = *lb1 * *lb2, abs(i__1)) == 28 && (*lb1 == 4 || *lb2 ==
		     4) || (i__2 = *lb1 * *lb2, abs(i__2)) == 32 && (*lb1 == 
		    4 || *lb2 == 4)) {
		branch = .33333333333333331f;
	    }
	    if (*lb1 * *lb2 == 30 && (*lb1 == 5 || *lb2 == 5) || *lb1 * *lb2 
		    == 27 && (*lb1 == 3 || *lb2 == 3) || (*lb1 * *lb2 == -18 
		    && (*lb1 == 3 || *lb2 == 3) || *lb1 * *lb2 == -45 && (*
		    lb1 == 5 || *lb2 == 5))) {
		branch = .5f;
	    }
	} else {
	    if (*lb1 * *lb2 == 40 && (*lb1 == 5 || *lb2 == 5) || *lb1 * *lb2 
		    == 30 && (*lb1 == 5 || *lb2 == 5) || (*lb1 * *lb2 == -24 
		    && (*lb1 == 3 || *lb2 == 3) || *lb1 * *lb2 == -18 && (*
		    lb1 == 3 || *lb2 == 3))) {
		branch = .40000000000000002f;
	    }
	    if (*lb1 * *lb2 == 27 && (*lb1 == 3 || *lb2 == 3) || *lb1 * *lb2 
		    == 21 && (*lb1 == 3 || *lb2 == 3) || (*lb1 * *lb2 == -45 
		    && (*lb1 == 5 || *lb2 == 5) || *lb1 * *lb2 == -35 && (*
		    lb1 == 5 || *lb2 == 5))) {
		branch = .40000000000000002f;
	    }
	    if (*lb1 * *lb2 == 35 && (*lb1 == 5 || *lb2 == 5) || *lb1 * *lb2 
		    == 24 && (*lb1 == 3 || *lb2 == 3) || (*lb1 * *lb2 == -21 
		    && (*lb1 == 3 || *lb2 == 3) || *lb1 * *lb2 == -40 && (*
		    lb1 == 5 || *lb2 == 5))) {
		branch = .53333333333333333f;
	    }
	    if ((i__1 = *lb1 * *lb2, abs(i__1)) == 28 && (*lb1 == 4 || *lb2 ==
		     4) || (i__2 = *lb1 * *lb2, abs(i__2)) == 32 && (*lb1 == 
		    4 || *lb2 == 4)) {
		branch = .066666666666666666f;
	    }
	    if ((i__1 = *lb1 * *lb2, abs(i__1)) == 36 && (*lb1 == 4 || *lb2 ==
		     4) || (i__2 = *lb1 * *lb2, abs(i__2)) == 24 && (*lb1 == 
		    4 || *lb2 == 4)) {
		branch = .59999999999999998f;
	    }
	}
/* bz11/25/98end */
	xs0 = fd2_(&arraym[ir - 1], &arrayj[ir - 1], &arrayl[ir - 1], &arrayw[
		ir - 1], &arrayb[ir - 1], em1, em2, srt);
	xs += pi * 1.3f * branch * xs0 * .038927290000000003f;
/* L1001: */
    }
    ret_val = xs;
    return ret_val;
} /* dpion_ */

/* **************************8 */
/* FUNCTION FDE(DMASS) GIVES DELTA MASS DISTRIBUTION BY USING OF */
/* KITAZOE'S FORMULA */
/*        REAL*4 FUNCTION FD2(DMASS,aj,al,width,widb0,EM1,EM2,srt) */
doublereal fd2_(real *dmass, real *aj, real *al, real *width, real *widb0, 
	real *em1, real *em2, real *srt)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4, r__5;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static real b, q, q0, ak2, ak02, amd, amp;

    amp = *em1;
    amd = *em2;
/* Computing 2nd power */
    r__2 = *dmass;
/* Computing 2nd power */
    r__3 = amd;
/* Computing 2nd power */
    r__4 = amp;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = amp * amd;
    ak02 = r__1 * r__1 * .25f - r__5 * r__5;
    if (ak02 > 0.f) {
	q0 = sqrt(ak02 / *dmass);
    } else {
	q0 = 0.f;
	ret_val = 0.f;
	return ret_val;
    }
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = amd;
/* Computing 2nd power */
    r__4 = amp;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = amp * amd;
    ak2 = r__1 * r__1 * .25f - r__5 * r__5;
    if (ak2 > 0.f) {
	q = sqrt(ak2 / *dmass);
    } else {
	q = 0.f;
	ret_val = 0.f;
	return ret_val;
    }
    d__1 = (doublereal) (q / q0);
    d__2 = (doublereal) (*al * 2.f + 1);
    d__3 = (doublereal) (q / q0);
    d__4 = (doublereal) (*al * 2);
    b = *widb0 * 1.2f * *dmass / *srt * pow_dd(&d__1, &d__2) / (pow_dd(&d__3, 
	    &d__4) * .2f + 1.f);
/* Computing 2nd power */
    r__1 = *width;
/* Computing 2nd power */
    r__2 = *srt - *dmass;
/* Computing 2nd power */
    r__3 = *width;
/* Computing 2nd power */
    r__4 = q;
    ret_val = (*aj * 2.f + 1) * (r__1 * r__1) * b / (r__2 * r__2 + r__3 * 
	    r__3 * .25f) / (r__4 * r__4 * 4.f);
    return ret_val;
} /* fd2_ */

/* **************************8 */
/*   MASS GENERATOR for two resonances simultaneously */
/* Subroutine */ int rmasdd_(real *srt, real *am10, real *am20, real *dmin1, 
	real *dmin2, integer *iseed, integer *ic, real *dm1, real *dm2)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5, r__6;

    /* Local variables */
    static real q2, fm1, fm2, fff, amn, amp, prob;
    static integer ntry;
    static real dmax1, dmax2, prob0;
    static integer ntry1, ntry2, ictrl;
    extern doublereal fmassd_(real *), ranart_(integer *), fmassn_(real *), 
	    fmassr_(real *);

/* c      SAVE /RNDF77/ */
    amn = .94f;
    amp = .14f;
/* the maximum mass for resonance 1 */
    dmax1 = *srt - *dmin2;
/* generate the mass for the first resonance */
L5:
    ntry1 = 0;
    ntry2 = 0;
    ntry = 0;
    ictrl = 0;
L10:
    *dm1 = ranart_(&rndf77_1.nseed) * (dmax1 - *dmin1) + *dmin1;
    ++ntry1;
/* the maximum mass for resonance 2 */
    if (ictrl == 0) {
	dmax2 = *srt - *dm1;
    }
/* generate the mass for the second resonance */
L20:
    *dm2 = ranart_(&rndf77_1.nseed) * (dmax2 - *dmin2) + *dmin2;
    ++ntry2;
/* check the energy-momentum conservation with two masses */
/* q2 in the following is q**2*4*srt**2 */
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = *dm1;
/* Computing 2nd power */
    r__4 = *dm2;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = *dm1;
/* Computing 2nd power */
    r__6 = *dm2;
    q2 = r__1 * r__1 - r__5 * r__5 * 4.f * (r__6 * r__6);
    if (q2 <= 0.f) {
	dmax2 = *dm2 - .01f;
/*         dmax1=dm1-0.01 */
	ictrl = 1;
	goto L20;
    }
/* determine the weight of the mass pair */
    if (dmax1 < *am10) {
	if (*ic == 1) {
	    fm1 = fmassd_(&dmax1);
	}
	if (*ic == 2) {
	    fm1 = fmassn_(&dmax1);
	}
	if (*ic == 3) {
	    fm1 = fmassd_(&dmax1);
	}
	if (*ic == 4) {
	    fm1 = fmassd_(&dmax1);
	}
    } else {
	if (*ic == 1) {
	    fm1 = fmassd_(am10);
	}
	if (*ic == 2) {
	    fm1 = fmassn_(am10);
	}
	if (*ic == 3) {
	    fm1 = fmassd_(am10);
	}
	if (*ic == 4) {
	    fm1 = fmassd_(am10);
	}
    }
    if (dmax2 < *am20) {
	if (*ic == 1) {
	    fm2 = fmassd_(&dmax2);
	}
	if (*ic == 2) {
	    fm2 = fmassn_(&dmax2);
	}
	if (*ic == 3) {
	    fm2 = fmassn_(&dmax2);
	}
	if (*ic == 4) {
	    fm2 = fmassr_(&dmax2);
	}
    } else {
	if (*ic == 1) {
	    fm2 = fmassd_(am20);
	}
	if (*ic == 2) {
	    fm2 = fmassn_(am20);
	}
	if (*ic == 3) {
	    fm2 = fmassn_(am20);
	}
	if (*ic == 4) {
	    fm2 = fmassr_(am20);
	}
    }
    if (fm1 == 0.f) {
	fm1 = 1e-4f;
    }
    if (fm2 == 0.f) {
	fm2 = 1e-4f;
    }
    prob0 = fm1 * fm2;
    if (*ic == 1) {
	prob = fmassd_(dm1) * fmassd_(dm2);
    }
    if (*ic == 2) {
	prob = fmassn_(dm1) * fmassn_(dm2);
    }
    if (*ic == 3) {
	prob = fmassd_(dm1) * fmassn_(dm2);
    }
    if (*ic == 4) {
	prob = fmassd_(dm1) * fmassr_(dm2);
    }
    if (prob <= 1e-6f) {
	prob = 1e-6f;
    }
    fff = prob / prob0;
    ++ntry;
    if (ranart_(&rndf77_1.nseed) > fff && ntry <= 20) {
	goto L10;
    }
/* lin-2/26/03 limit the mass of (rho,Delta,N*1440) below a certain value */
/*     (here taken as its central value + 2* B-W fullwidth): */
    if ((r__1 = *am10 - .77f, dabs(r__1)) <= .01f && *dm1 > 1.07f || (r__2 = *
	    am10 - 1.232f, dabs(r__2)) <= .01f && *dm1 > 1.47f || (r__3 = *
	    am10 - 1.44f, dabs(r__3)) <= .01f && *dm1 > 2.14f) {
	goto L5;
    }
    if ((r__1 = *am20 - .77f, dabs(r__1)) <= .01f && *dm2 > 1.07f || (r__2 = *
	    am20 - 1.232f, dabs(r__2)) <= .01f && *dm2 > 1.47f || (r__3 = *
	    am20 - 1.44f, dabs(r__3)) <= .01f && *dm2 > 2.14f) {
	goto L5;
    }
    return 0;
} /* rmasdd_ */

/* FUNCTION Fmassd(DMASS) GIVES the delta MASS DISTRIBUTION */
doublereal fmassd_(real *dmass)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4, r__5;

    /* Local variables */
    static real am0;
    extern doublereal width_(real *);

    am0 = 1.232f;
/* Computing 2nd power */
    r__2 = *dmass;
/* Computing 2nd power */
    r__3 = am0;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3;
/* Computing 2nd power */
    r__4 = am0;
/* Computing 2nd power */
    r__5 = width_(dmass);
    ret_val = am0 * width_(dmass) / (r__1 * r__1 + r__4 * r__4 * (r__5 * r__5)
	    );
    return ret_val;
} /* fmassd_ */

/* FUNCTION Fmassn(DMASS) GIVES the N* MASS DISTRIBUTION */
doublereal fmassn_(real *dmass)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4, r__5;

    /* Local variables */
    static real am0;
    extern doublereal w1440_(real *);

    am0 = 1.44f;
/* Computing 2nd power */
    r__2 = *dmass;
/* Computing 2nd power */
    r__3 = am0;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3;
/* Computing 2nd power */
    r__4 = am0;
/* Computing 2nd power */
    r__5 = w1440_(dmass);
    ret_val = am0 * w1440_(dmass) / (r__1 * r__1 + r__4 * r__4 * (r__5 * r__5)
	    );
    return ret_val;
} /* fmassn_ */

/* FUNCTION Fmassr(DMASS) GIVES the rho MASS DISTRIBUTION */
doublereal fmassr_(real *dmass)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4, r__5;

    /* Local variables */
    static real am0, wid;

    am0 = .77f;
    wid = .153f;
/* Computing 2nd power */
    r__2 = *dmass;
/* Computing 2nd power */
    r__3 = am0;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3;
/* Computing 2nd power */
    r__4 = am0;
/* Computing 2nd power */
    r__5 = wid;
    ret_val = am0 * wid / (r__1 * r__1 + r__4 * r__4 * (r__5 * r__5));
    return ret_val;
} /* fmassr_ */

/* ********************************* */
/* PURPOSE : flow analysis */
/* DATE : Feb. 1, 1995 */
/* ********************************** */
/* Subroutine */ int flow_(integer *nt)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2, r__3, r__4;

    /* Builtin functions */
    integer i_nint(real *);
    double sqrt(doublereal), log(doublereal);

    /* Local variables */
    static integer i__, j, m;
    static real e00;
    static integer kk;
    static real y00;
    static integer is;
    static real dy;
    static integer iy, ly, npr, npt;
    static real ypr[161], dnuc, dypr;
    static integer nrun;
    static real ycut1, ycut2, dnuck;
    static integer nkaon;
    static real dnucp, ykaon[161];
    static integer npion;
    static real ypion[161], pxpro[161], dykaon, pxkaon[161], dypion, pxpion[
	    161];

/*       IMPLICIT REAL*4 (A-H,O-Z) */
/* ----------------------------------------------------------------------* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /RR/ */
/* c      SAVE /RUN/ */
/* c      SAVE /input1/ */
/* ----------------------------------------------------------------------* */
    ycut1 = -2.6f;
    ycut2 = 2.6f;
    dy = .2f;
    r__1 = (ycut2 - ycut1) / dy;
    ly = i_nint(&r__1);
/* ********************************** */
/* initialize the transverse momentum counters */
    for (kk = -80; kk <= 80; ++kk) {
	pxpion[kk + 80] = 0.f;
	pxpro[kk + 80] = 0.f;
	pxkaon[kk + 80] = 0.f;
/* L11: */
    }
    i__1 = ly;
    for (j = -ly; j <= i__1; ++j) {
	ypion[j + 80] = 0.f;
	ykaon[j + 80] = 0.f;
	ypr[j + 80] = 0.f;
/* L701: */
    }
    nkaon = 0;
    npr = 0;
    npion = 0;
    is = 0;
    i__1 = run_1.num;
    for (nrun = 1; nrun <= i__1; ++nrun) {
	is += rr_1.massr[nrun - 1];
	i__2 = rr_1.massr[nrun];
	for (j = 1; j <= i__2; ++j) {
	    i__ = j + is;
/* for protons go to 200 to calculate its rapidity and transvese momentum */
/* distributions */
/* Computing 2nd power */
	    r__1 = bb_1.p[i__ * 3 - 3];
/* Computing 2nd power */
	    r__2 = bb_1.p[i__ * 3 - 2];
/* Computing 2nd power */
	    r__3 = bb_1.p[i__ * 3 - 1];
/* Computing 2nd power */
	    r__4 = cc_1.e[i__ - 1];
	    e00 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	    y00 = log((e00 + bb_1.p[i__ * 3 - 1]) / (e00 - bb_1.p[i__ * 3 - 1]
		    )) * .5f;
	    if (dabs(y00) >= ycut2) {
		goto L20;
	    }
	    r__1 = y00 / dy;
	    iy = i_nint(&r__1);
	    if (abs(iy) >= 80) {
		goto L20;
	    }
	    if (cc_1.e[i__ - 1] == 0.f) {
		goto L20;
	    }
	    if (ee_1.lb[i__ - 1] >= 25) {
		goto L20;
	    }
	    if (ee_1.lb[i__ - 1] <= 5 && ee_1.lb[i__ - 1] >= 3) {
		goto L50;
	    }
	    if (ee_1.lb[i__ - 1] == 1 || ee_1.lb[i__ - 1] == 2) {
		goto L200;
	    }
/* bz3/10/99 */
/*       if(lb(i).ge.6.and.lb(i).le.15)go to 200 */
	    if (ee_1.lb[i__ - 1] >= 6 && ee_1.lb[i__ - 1] <= 17) {
		goto L200;
	    }
/* bz3/10/99 end */
	    if (ee_1.lb[i__ - 1] == 23) {
		goto L400;
	    }
	    goto L20;
/* calculate rapidity and transverse momentum distribution for pions */
L50:
	    ++npion;
/* (2) rapidity distribution in the cms frame */
	    ypion[iy + 80] += 1;
	    pxpion[iy + 80] += bb_1.p[i__ * 3 - 3] / cc_1.e[i__ - 1];
	    goto L20;
/* calculate rapidity and transverse energy distribution for baryons */
L200:
	    ++npr;
	    pxpro[iy + 80] += bb_1.p[i__ * 3 - 3] / cc_1.e[i__ - 1];
	    ypr[iy + 80] += 1.f;
	    goto L20;
L400:
	    ++nkaon;
	    ykaon[iy + 80] += 1.f;
	    pxkaon[iy + 80] += bb_1.p[i__ * 3 - 3] / cc_1.e[i__ - 1];
L20:
	    ;
	}
    }
/* PRINT OUT NUCLEON'S TRANSVERSE MOMENTUM distribution */
/*       write(1041,*)Nt */
/*       write(1042,*)Nt */
/*       write(1043,*)Nt */
/*       write(1090,*)Nt */
/*       write(1091,*)Nt */
/*       write(1092,*)Nt */
    for (npt = -10; npt <= 10; ++npt) {
	if (ypr[npt + 80] == 0.f) {
	    goto L101;
	}
	pxpro[npt + 80] = -pxpro[npt + 80] / ypr[npt + 80];
	dnuc = pxpro[npt + 80] / sqrt(ypr[npt + 80]);
/*       WRITE(1041,*)NPT*DY,Pxpro(NPT),DNUC */
/* print pion's transverse momentum distribution */
L101:
	if (ypion[npt + 80] == 0.f) {
	    goto L102;
	}
	pxpion[npt + 80] = -pxpion[npt + 80] / ypion[npt + 80];
	dnucp = pxpion[npt + 80] / sqrt(ypion[npt + 80]);
/*       WRITE(1042,*)NPT*DY,Pxpion(NPT),DNUCp */
/* kaons */
L102:
	if (ykaon[npt + 80] == 0.f) {
	    goto L3;
	}
	pxkaon[npt + 80] = -pxkaon[npt + 80] / ykaon[npt + 80];
	dnuck = pxkaon[npt + 80] / sqrt(ykaon[npt + 80]);
/*       WRITE(1043,*)NPT*DY,Pxkaon(NPT),DNUCk */
L3:
	;
    }
/* ******************************* */
/* OUTPUT PION AND PROTON RAPIDITY DISTRIBUTIONS */
    i__2 = ly;
    for (m = -ly; m <= i__2; ++m) {
/* PROTONS */
	dypr = 0.f;
	if (ypr[m + 80] != 0.f) {
	    dypr = sqrt(ypr[m + 80]) / (real) nrun / dy;
	}
	ypr[m + 80] = ypr[m + 80] / (real) nrun / dy;
/*       WRITE(1090,'(E11.3,2X,E11.3,2X,E11.3)')m*DY,YPR(M),DYPR */
/* PIONS */
	dypion = 0.f;
	if (ypion[m + 80] != 0.f) {
	    dypion = sqrt(ypion[m + 80]) / (real) nrun / dy;
	}
	ypion[m + 80] = ypion[m + 80] / (real) nrun / dy;
/*       WRITE(1091,'(E11.3,2X,E11.3,2X,E11.3)')m*DY,YPION(M),DYPION */
/* KAONS */
	dykaon = 0.f;
	if (ykaon[m + 80] != 0.f) {
	    dykaon = sqrt(ykaon[m + 80]) / (real) nrun / dy;
	}
	ykaon[m + 80] = ykaon[m + 80] / (real) nrun / dy;
/*       WRITE(1092,'(E11.3,2X,E11.3,2X,E11.3)')m*DY,YKAON(M),DYKAON */
/* L1001: */
    }
    return 0;
} /* flow_ */

/* bali1/16/99 */
/* ******************************************* */
/* Purpose: pp_bar annihilation cross section as a functon of their cms energy */
/*      real*4 function xppbar(srt) */
doublereal xppbar_(real *srt)
{
    /* System generated locals */
    real ret_val, r__1, r__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static real plab, plab2;

/*  srt    = DSQRT(s) in GeV                                                   * */
/*  xppbar = pp_bar annihilation cross section in mb                           * */

/*  Reference: G.J. Wang, R. Bellwied, C. Pruneau and G. Welke */
/*             Proc. of the 14th Winter Workshop on Nuclear Dynamics, */
/*             Snowbird, Utah 31, Eds. W. Bauer and H.G. Ritter */
/*             (Plenum Publishing, 1998)                             * */

/* ***************************************** */
/* Note: */
/* (1) we introduce a new parameter xmax=400 mb: */
/*     the maximum annihilation xsection */
/* there are shadowing effects in pp_bar annihilation, with this parameter */
/* we can probably look at these effects */
/* (2) Calculate p(lab) from srt [GeV], since the formular in the */
/* reference applies only to the case of a p_bar on a proton at rest */
/* Formula used: srt**2=2.*pmass*(pmass+sqrt(pmass**2+plab**2)) */
    ret_val = 1e-6f;
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__1 = r__2 * r__2 / 1.8766f - .9383f;
    plab2 = r__1 * r__1 - .88040689000000005f;
    if (plab2 > 0.f) {
	plab = sqrt(plab2);
	d__1 = (doublereal) plab;
	ret_val = 67.f / pow_dd(&d__1, &c_b927);
	if (ret_val > 400.f) {
	    ret_val = 400.f;
	}
    }
    return ret_val;
} /* xppbar_ */

/* bali1/16/99 end */
/* ********************************* */
/* bali2/6/99 */
/* ******************************************* */
/* Purpose: To generate randomly the no. of pions in the final */
/*          state of pp_bar annihilation according to a statistical */
/*          model by using of the rejection method. */
/* bz2/25/99 */
/*      real*4 function pbarfs(srt,npion,iseed) */
/* Subroutine */ int pbarfs_(real *srt, integer *npion, integer *iseed)
{
    /* System generated locals */
    real r__1;

    /* Builtin functions */
    double pow_ri(real *, integer *);

    /* Local variables */
    static integer n;
    static real ene, pmax, pnpi[6];
    static integer ntry;
    static real thisp, factor[6];
    extern doublereal ranart_(integer *);

/* bz2/25/99end */
/* Quantities: */
/*  srt: DSQRT(s) in GeV                                                    * */
/*  npion: No. of pions produced in the annihilation of ppbar at srt        * */
/*  nmax=6, cutoff of the maximum no. of n the code can handle */

/*  Reference: C.M. Ko and R. Yuan, Phys. Lett. B192 (1987) 31      * */

/* ***************************************** */
/* c      SAVE /RNDF77/ */
/* the factorial coefficients in the pion no. distribution */
/* from n=2 to 6 calculated use the formula in the reference */
    factor[1] = 1.f;
    factor[2] = .117f;
    factor[3] = .00327f;
    factor[4] = 3.58e-5f;
    factor[5] = 1.93e-7f;
/* Computing 3rd power */
    r__1 = *srt / .14f;
    ene = r__1 * (r__1 * r__1) / 59.217624386248559f;
/* the relative probability from n=2 to 6 */
    for (n = 2; n <= 6; ++n) {
	pnpi[n - 1] = pow_ri(&ene, &n) * factor[n - 1];
/* L1001: */
    }
/* find the maximum of the probabilities, I checked a */
/* Fortan manual: max() returns the maximum value of */
/* the same type as in the argument list */
/* Computing MAX */
    r__1 = max(pnpi[1],pnpi[2]), r__1 = max(r__1,pnpi[3]), r__1 = max(r__1,
	    pnpi[4]);
    pmax = dmax(r__1,pnpi[5]);
/* randomly generate n between 2 and 6 */
    ntry = 0;
L10:
    *npion = (integer) (ranart_(&rndf77_1.nseed) * 5) + 2;
/* lin-4/2008 check bounds: */
    if (*npion > 6) {
	goto L10;
    }
    thisp = pnpi[*npion - 1] / pmax;
    ++ntry;
/* decide whether to take this npion according to the distribution */
/* using rejection method. */
    if (thisp < ranart_(&rndf77_1.nseed) && ntry <= 20) {
	goto L10;
    }
/* now take the last generated npion and return */
    return 0;
} /* pbarfs_ */

/* ********************************* */
/* bali2/6/99 end */
/* bz3/9/99 kkbar */
/* bali3/5/99 */
/* ***************************************** */
/* purpose: Xsection for K+ K- to pi+ pi- */
/*      real*4 function xkkpi(srt) */
/*  srt    = DSQRT(s) in GeV                                  * */
/*  xkkpi   = xsection in mb obtained from */
/*           the detailed balance                             * */
/* ****************************************** */
/*          parameter (pimass=0.140,aka=0.498) */
/*       xkkpi=1.e-08 */
/*       ppi2=(srt/2)**2-pimass**2 */
/*       pk2=(srt/2)**2-aka**2 */
/*       if(ppi2.le.0.or.pk2.le.0)return */
/* bz3/9/99 kkbar */
/*       xkkpi=ppi2/pk2*pipik(srt) */
/*       xkkpi=9.0 / 4.0 * ppi2/pk2*pipik(srt) */
/*        xkkpi = 2.0 * xkkpi */
/* bz3/9/99 kkbar end */
/* bz3/9/99 kkbar */
/*       end */
/*       return */
/*        END */
/* bz3/9/99 kkbar end */
/* bali3/5/99 end */
/* bz3/9/99 kkbar end */
/* bz3/9/99 kkbar */
/* **************************** */
/* purpose: Xsection for K+ K- to pi+ pi- */
/* Subroutine */ int xkkann_(real *srt, real *xsk1, real *xsk2, real *xsk3, 
	real *xsk4, real *xsk5, real *xsk6, real *xsk7, real *xsk8, real *
	xsk9, real *xsk10, real *xsk11, real *sigk, real *rrkk)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    static real s, pf2, pi2, xm1, xm2, fwdp, pkaon;
    extern doublereal pipik_(real *);
    static real xpion0, scheck;

    /* Fortran I/O blocks */
    static cilist io___2196 = { 0, 99, 0, 0, 0 };


/*  srt    = DSQRT(s) in GeV                                       * */
/*  xsk1   = annihilation into pi pi                               * */
/*  xsk2   = annihilation into pi rho (shifted to XKKSAN)         * */
/*  xsk3   = annihilation into pi omega (shifted to XKKSAN)       * */
/*  xsk4   = annihilation into pi eta                              * */
/*  xsk5   = annihilation into rho rho                             * */
/*  xsk6   = annihilation into rho omega                           * */
/*  xsk7   = annihilation into rho eta (shifted to XKKSAN)        * */
/*  xsk8   = annihilation into omega omega                         * */
/*  xsk9   = annihilation into omega eta (shifted to XKKSAN)      * */
/*  xsk10  = annihilation into eta eta                             * */
/*  sigk   = xsection in mb obtained from                          * */
/*           the detailed balance                                  * */
/* *************************** */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /EE/ */
/* c      SAVE /DD/ */
/* Computing 2nd power */
    r__1 = *srt;
    s = r__1 * r__1;
    *sigk = 1e-8f;
    *xsk1 = 0.f;
    *xsk2 = 0.f;
    *xsk3 = 0.f;
    *xsk4 = 0.f;
    *xsk5 = 0.f;
    *xsk6 = 0.f;
    *xsk7 = 0.f;
    *xsk8 = 0.f;
    *xsk9 = 0.f;
    *xsk10 = 0.f;
    *xsk11 = 0.f;
    xpion0 = pipik_(srt);
/* .....take into account both K+ and K0 */
    xpion0 *= 2.f;
    pi2 = s * (s - .99201600000000001f);
    if (pi2 <= 0.f) {
	return 0;
    }
    xm1 = .14f;
    xm2 = .14f;
/* Computing 2nd power */
    r__1 = xm1 + xm2;
/* Computing 2nd power */
    r__2 = xm1 - xm2;
    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (pf2 > 0.f) {
	*xsk1 = pf2 * 2.25f / pi2 * xpion0;
    }
/* lin-8/28/00 (pi eta) eta -> K+K- is assumed the same as pi pi -> K+K-: */
    xm1 = .14f;
    xm2 = .5473f;
/* Computing 2nd power */
    r__1 = xm1 + xm2;
/* Computing 2nd power */
    r__2 = xm1 - xm2;
    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (pf2 > 0.f) {
	*xsk4 = pf2 * .75f / pi2 * xpion0;
    }
    xm1 = .5473f;
    xm2 = .5473f;
/* Computing 2nd power */
    r__1 = xm1 + xm2;
/* Computing 2nd power */
    r__2 = xm1 - xm2;
    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (pf2 > 0.f) {
	*xsk10 = pf2 * .25f / pi2 * xpion0;
    }
    xpion0 = *rrkk;
/* lin-11/07/00: (pi eta) (rho omega) -> K* Kbar (or K*bar K) instead to K Kbar: */
/*        XM1 = PIMASS */
/*        XM2 = RHOM */
/*        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2) */
/*        IF (PF2 .GT. 0.0) THEN */
/*           XSK2 = 27.0 / 4.0 * PF2 / PI2 * XPION0 */
/*        END IF */
/*        XM1 = PIMASS */
/*        XM2 = OMEGAM */
/*        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2) */
/*        IF (PF2 .GT. 0.0) THEN */
/*           XSK3 = 9.0 / 4.0 * PF2 / PI2 * XPION0 */
/*        END IF */
    xm1 = .77f;
    xm2 = .77f;
/* Computing 2nd power */
    r__1 = xm1 + xm2;
/* Computing 2nd power */
    r__2 = xm1 - xm2;
    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (pf2 > 0.f) {
	*xsk5 = pf2 * 20.25f / pi2 * xpion0;
    }
    xm1 = .77f;
    xm2 = .7819f;
/* Computing 2nd power */
    r__1 = xm1 + xm2;
/* Computing 2nd power */
    r__2 = xm1 - xm2;
    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (pf2 > 0.f) {
	*xsk6 = pf2 * 6.75f / pi2 * xpion0;
    }
/*        XM1 = RHOM */
/*        XM2 = ETAM */
/*        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2) */
/*        IF (PF2 .GT. 0.0) THEN */
/*           XSK7 = 9.0 / 4.0 * PF2 / PI2 * XPION0 */
/*        END IF */
    xm1 = .7819f;
    xm2 = .7819f;
/* Computing 2nd power */
    r__1 = xm1 + xm2;
/* Computing 2nd power */
    r__2 = xm1 - xm2;
    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (pf2 > 0.f) {
	*xsk8 = pf2 * 2.25f / pi2 * xpion0;
    }
/*        XM1 = OMEGAM */
/*        XM2 = ETAM */
/*        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2) */
/*        IF (PF2 .GT. 0.0) THEN */
/*           XSK9 = 3.0 / 4.0 * PF2 / PI2 * XPION0 */
/*        END IF */
/* * K+ + K- --> phi */
    fwdp = pow_dd(&c_b932, &c_b933) * 1.68f / 6.f / 1.02f / 1.02f;
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__1 = *srt;
    scheck = r__1 * r__1 - .99201600000000001f;
    if (scheck <= 0.f) {
	s_wsle(&io___2196);
	do_lio(&c__9, &c__1, "scheck47: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	s_stop("", (ftnlen)0);
    }
    pkaon = sqrt(scheck) * .5f;
/*          pkaon=0.5*sqrt(srt**2-4.0*aka**2) */
/* Computing 2nd power */
    r__1 = fwdp * 1.02f;
/* Computing 2nd power */
    r__3 = *srt;
/* Computing 2nd power */
    r__2 = r__3 * r__3 - 1.0404f;
/* Computing 2nd power */
    r__4 = fwdp * 1.02f;
/* Computing 2nd power */
    r__5 = pkaon;
    *xsk11 = r__1 * r__1 * 3.6688075497330002f / (r__2 * r__2 + r__4 * r__4) /
	     (r__5 * r__5);

    *sigk = *xsk1 + *xsk2 + *xsk3 + *xsk4 + *xsk5 + *xsk6 + *xsk7 + *xsk8 + *
	    xsk9 + *xsk10 + *xsk11;
    return 0;
} /* xkkann_ */

/* bz3/9/99 kkbar end */
/* **************************** */
/* purpose: Xsection for Phi + B */
/* Subroutine */ int xphib_(integer *lb1, integer *lb2, real *em1, real *em2, 
	real *srt, real *xsk1, real *xsk2, real *xsk3, real *xsk4, real *xsk5,
	 real *sigp)
{
    /* System generated locals */
    real r__1;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static real xsk6, srrt;


/* *************************** */
    *sigp = 1e-8f;
    *xsk1 = 0.f;
    *xsk2 = 0.f;
    *xsk3 = 0.f;
    *xsk4 = 0.f;
    *xsk5 = 0.f;
    xsk6 = 0.f;
    srrt = *srt - (*em1 + *em2);
/* * phi + N(D) -> elastic scattering */
/*            XSK1 = 0.56  !! mb */
/*  !! mb  (photo-production xsecn used) */
    *xsk1 = 8.f;

/* * phi + N(D) -> pi + N */
    if (*srt > 1.074417f) {
	d__1 = (doublereal) srrt;
	*xsk2 = pow_dd(&d__1, &c_b941) * .0235f;
    }

/* * phi + N(D) -> pi + D */
    if (*srt > 1.36696f) {
	if (srrt < .7f) {
	    d__1 = (doublereal) srrt;
	    *xsk3 = pow_dd(&d__1, &c_b942) * .0119f;
	} else {
	    d__1 = (doublereal) srrt;
	    *xsk3 = pow_dd(&d__1, &c_b943) * .013f;
	}
    }

/* * phi + N(D) -> rho + N */
    if (*srt > 1.709457f) {
	if (srrt < .7f) {
	    d__1 = (doublereal) srrt;
	    *xsk4 = pow_dd(&d__1, &c_b944) * .0166f;
	} else {
	    d__1 = (doublereal) srrt;
	    *xsk4 = pow_dd(&d__1, &c_b945) * .0189f;
	}
    }

/* * phi + N(D) -> rho + D   (same as pi + D) */
    if (*srt > 2.0019999999999998f) {
	if (srrt < .7f) {
	    d__1 = (doublereal) srrt;
	    *xsk5 = pow_dd(&d__1, &c_b942) * .0119f;
	} else {
	    d__1 = (doublereal) srrt;
	    *xsk5 = pow_dd(&d__1, &c_b943) * .013f;
	}
    }

/* * phi + N -> K+ + La */
    if (*lb1 >= 1 && *lb1 <= 2 || *lb2 >= 1 && *lb2 <= 2) {
	if (*srt > 1.6136999999999999f) {
/* Computing 2nd power */
	    r__1 = srrt + 3.508f;
	    xsk6 = 1.715f / (r__1 * r__1 - 12.138f);
	}
    }
    *sigp = *xsk1 + *xsk2 + *xsk3 + *xsk4 + *xsk5 + xsk6;
    return 0;
} /* xphib_ */


/* ********************************* */

/* Subroutine */ int crphib_(real *px, real *py, real *pz, real *srt, integer 
	*i1, integer *i2, real *xsk1, real *xsk2, real *xsk3, real *xsk4, 
	real *xsk5, real *sigp, integer *iblock)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    double sqrt(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static real c1, s1, t1, x1, pr, em1, em2, ct1, pr2, px0, py0, pz0, st1;
    extern doublereal ranart_(integer *);
    extern /* Subroutine */ int rotate_(real *, real *, real *, real *, real *
	    , real *);


/*     PURPOSE:                                                         * */
/*             DEALING WITH PHI + N(D) --> pi+N(D), rho+N(D),  K+ + La */
/*     QUANTITIES:                                                      * */
/*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME* */
/*           SRT      - SQRT OF S                                       * */
/*           IBLOCK   - INFORMATION about the reaction channel          * */

/*             iblock   - 20  elastic */
/*             iblock   - 221  K+ formation */
/*             iblock   - 223  others */
/* ********************************* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /input1/ */
/* c      SAVE /RNDF77/ */

    px0 = *px;
    py0 = *py;
    pz0 = *pz;
    *iblock = 223;

    x1 = ranart_(&rndf77_1.nseed) * *sigp;
    *xsk2 = *xsk1 + *xsk2;
    *xsk3 = *xsk2 + *xsk3;
    *xsk4 = *xsk3 + *xsk4;
    *xsk5 = *xsk4 + *xsk5;

/*  !! elastic scatt. */
    if (x1 <= *xsk1) {
	*iblock = 20;
	goto L100;
    } else if (x1 <= *xsk2) {
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	cc_1.e[*i1 - 1] = .13496f;
	cc_1.e[*i2 - 1] = .939457f;
	goto L100;
    } else if (x1 <= *xsk3) {
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 4) + 6;
	cc_1.e[*i1 - 1] = .13496f;
	cc_1.e[*i2 - 1] = 1.232f;
	goto L100;
    } else if (x1 <= *xsk4) {
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 25;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
	cc_1.e[*i1 - 1] = .77f;
	cc_1.e[*i2 - 1] = .939457f;
	goto L100;
    } else if (x1 <= *xsk5) {
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 25;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 4) + 6;
	cc_1.e[*i1 - 1] = .77f;
	cc_1.e[*i2 - 1] = 1.232f;
	goto L100;
    } else {
	ee_1.lb[*i1 - 1] = 23;
	ee_1.lb[*i2 - 1] = 14;
	cc_1.e[*i1 - 1] = .498f;
	cc_1.e[*i2 - 1] = 1.1157f;
	*iblock = 221;
    }
L100:
    em1 = cc_1.e[*i1 - 1];
    em2 = cc_1.e[*i2 - 1];
/* ----------------------------------------------------------------------- */
/* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH */
/* ENERGY CONSERVATION */
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = em1;
/* Computing 2nd power */
    r__4 = em2;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = em1 * em2;
    pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
    if (pr2 <= 0.f) {
	pr2 = 1e-8f;
    }
    pr = sqrt(pr2) / (*srt * 2.f);
/* WE ASSUME AN ISOTROPIC ANGULAR DISTRIBUTION IN THE CMS */
    c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
    t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
/* Computing 2nd power */
    r__1 = c1;
    s1 = sqrt(1.f - r__1 * r__1);
    ct1 = cos(t1);
    st1 = sin(t1);
/* THE MOMENTUM IN THE CMS IN THE FINAL STATE */
    *pz = pr * c1;
    *px = pr * s1 * ct1;
    *py = pr * s1 * st1;
/* ROTATE IT */
    rotate_(&px0, &py0, &pz0, px, py, pz);
    return 0;
} /* crphib_ */


/* **************************** */
/* purpose: Xsection for Phi + B */
/* !! in fm^2 */
/* Subroutine */ int pibphi_(real *srt, integer *lb1, integer *lb2, real *em1,
	 real *em2, real *xphi, real *xphin)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5, r__6;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static real sig, srrt, xphid;


/*      phi + N(D) <- pi + N */
/*      phi + N(D) <- pi + D */
/*      phi + N(D) <- rho + N */
/*      phi + N(D) <- rho + D   (same as pi + D) */

/* *************************** */
    *xphi = 0.f;
    *xphin = 0.f;
    xphid = 0.f;

    if (*lb1 >= 3 && *lb1 <= 5 || *lb2 >= 3 && *lb2 <= 5) {

	if (abs(*lb1) >= 1 && abs(*lb1) <= 2 || abs(*lb2) >= 1 && abs(*lb2) <=
		 2) {
/* * phi + N <- pi + N */
	    if (*srt > 1.959457f) {
		srrt = *srt - 1.959457f;
		d__1 = (doublereal) srrt;
		sig = pow_dd(&d__1, &c_b941) * .0235f;
/* Computing 2nd power */
		r__1 = *srt;
/* Computing 2nd power */
		r__2 = *srt;
/* Computing 2nd power */
		r__3 = *srt;
/* Computing 2nd power */
		r__4 = *em1 + *em2;
/* Computing 2nd power */
		r__5 = *srt;
/* Computing 2nd power */
		r__6 = *em1 - *em2;
		*xphin = sig * 1.f * (r__1 * r__1 - 3.839471734849f) * (r__2 *
			 r__2 - .0064871748490000049f) / (r__3 * r__3 - r__4 *
			 r__4) / (r__5 * r__5 - r__6 * r__6);
	    }
/* * phi + D <- pi + N */
	    if (*srt > 2.2519999999999998f) {
		srrt = *srt - 2.2519999999999998f;
		d__1 = (doublereal) srrt;
		sig = pow_dd(&d__1, &c_b941) * .0235f;
/* Computing 2nd power */
		r__1 = *srt;
/* Computing 2nd power */
		r__2 = *srt;
/* Computing 2nd power */
		r__3 = *srt;
/* Computing 2nd power */
		r__4 = *em1 + *em2;
/* Computing 2nd power */
		r__5 = *srt;
/* Computing 2nd power */
		r__6 = *em1 - *em2;
		xphid = sig * 4.f * (r__1 * r__1 - 5.0715039999999991f) * (
			r__2 * r__2 - .044943999999999984f) / (r__3 * r__3 - 
			r__4 * r__4) / (r__5 * r__5 - r__6 * r__6);
	    }
	} else {
/* * phi + N <- pi + D */
	    if (*srt > 1.959457f) {
		srrt = *srt - 1.959457f;
		if (srrt < .7f) {
		    d__1 = (doublereal) srrt;
		    sig = pow_dd(&d__1, &c_b942) * .0119f;
		} else {
		    d__1 = (doublereal) srrt;
		    sig = pow_dd(&d__1, &c_b943) * .013f;
		}
/* Computing 2nd power */
		r__1 = *srt;
/* Computing 2nd power */
		r__2 = *srt;
/* Computing 2nd power */
		r__3 = *srt;
/* Computing 2nd power */
		r__4 = *em1 + *em2;
/* Computing 2nd power */
		r__5 = *srt;
/* Computing 2nd power */
		r__6 = *em1 - *em2;
		*xphin = sig * .25f * (r__1 * r__1 - 3.839471734849f) * (r__2 
			* r__2 - .0064871748490000049f) / (r__3 * r__3 - r__4 
			* r__4) / (r__5 * r__5 - r__6 * r__6);
	    }
/* * phi + D <- pi + D */
	    if (*srt > 2.2519999999999998f) {
		srrt = *srt - 2.2519999999999998f;
		if (srrt < .7f) {
		    d__1 = (doublereal) srrt;
		    sig = pow_dd(&d__1, &c_b942) * .0119f;
		} else {
		    d__1 = (doublereal) srrt;
		    sig = pow_dd(&d__1, &c_b943) * .013f;
		}
/* Computing 2nd power */
		r__1 = *srt;
/* Computing 2nd power */
		r__2 = *srt;
/* Computing 2nd power */
		r__3 = *srt;
/* Computing 2nd power */
		r__4 = *em1 + *em2;
/* Computing 2nd power */
		r__5 = *srt;
/* Computing 2nd power */
		r__6 = *em1 - *em2;
		xphid = sig * 1.f * (r__1 * r__1 - 5.0715039999999991f) * (
			r__2 * r__2 - .044943999999999984f) / (r__3 * r__3 - 
			r__4 * r__4) / (r__5 * r__5 - r__6 * r__6);
	    }
	}


/* ** for rho + N(D) colln */

    } else {

	if (abs(*lb1) >= 1 && abs(*lb1) <= 2 || abs(*lb2) >= 1 && abs(*lb2) <=
		 2) {

/* * phi + N <- rho + N */
	    if (*srt > 1.959457f) {
		srrt = *srt - 1.959457f;
		if (srrt < .7f) {
		    d__1 = (doublereal) srrt;
		    sig = pow_dd(&d__1, &c_b944) * .0166f;
		} else {
		    d__1 = (doublereal) srrt;
		    sig = pow_dd(&d__1, &c_b945) * .0189f;
		}
/* Computing 2nd power */
		r__1 = *srt;
/* Computing 2nd power */
		r__2 = *srt;
/* Computing 2nd power */
		r__3 = *srt;
/* Computing 2nd power */
		r__4 = *em1 + *em2;
/* Computing 2nd power */
		r__5 = *srt;
/* Computing 2nd power */
		r__6 = *em1 - *em2;
		*xphin = sig * .33333333333333331f * (r__1 * r__1 - 
			3.839471734849f) * (r__2 * r__2 - 
			.0064871748490000049f) / (r__3 * r__3 - r__4 * r__4) /
			 (r__5 * r__5 - r__6 * r__6);
	    }
/* * phi + D <- rho + N */
	    if (*srt > 2.2519999999999998f) {
		srrt = *srt - 2.2519999999999998f;
		if (srrt < .7f) {
		    d__1 = (doublereal) srrt;
		    sig = pow_dd(&d__1, &c_b944) * .0166f;
		} else {
		    d__1 = (doublereal) srrt;
		    sig = pow_dd(&d__1, &c_b945) * .0189f;
		}
/* Computing 2nd power */
		r__1 = *srt;
/* Computing 2nd power */
		r__2 = *srt;
/* Computing 2nd power */
		r__3 = *srt;
/* Computing 2nd power */
		r__4 = *em1 + *em2;
/* Computing 2nd power */
		r__5 = *srt;
/* Computing 2nd power */
		r__6 = *em1 - *em2;
		xphid = sig * 1.3333333333333333f * (r__1 * r__1 - 
			5.0715039999999991f) * (r__2 * r__2 - 
			.044943999999999984f) / (r__3 * r__3 - r__4 * r__4) / 
			(r__5 * r__5 - r__6 * r__6);
	    }
	} else {
/* * phi + N <- rho + D  (same as pi+D->phi+N) */
	    if (*srt > 1.959457f) {
		srrt = *srt - 1.959457f;
		if (srrt < .7f) {
		    d__1 = (doublereal) srrt;
		    sig = pow_dd(&d__1, &c_b942) * .0119f;
		} else {
		    d__1 = (doublereal) srrt;
		    sig = pow_dd(&d__1, &c_b943) * .013f;
		}
/* Computing 2nd power */
		r__1 = *srt;
/* Computing 2nd power */
		r__2 = *srt;
/* Computing 2nd power */
		r__3 = *srt;
/* Computing 2nd power */
		r__4 = *em1 + *em2;
/* Computing 2nd power */
		r__5 = *srt;
/* Computing 2nd power */
		r__6 = *em1 - *em2;
		*xphin = sig * .083333333333333329f * (r__1 * r__1 - 
			3.839471734849f) * (r__2 * r__2 - 
			.0064871748490000049f) / (r__3 * r__3 - r__4 * r__4) /
			 (r__5 * r__5 - r__6 * r__6);
	    }
/* * phi + D <- rho + D  (same as pi+D->phi+D) */
	    if (*srt > 2.2519999999999998f) {
		srrt = *srt - 2.2519999999999998f;
		if (srrt < .7f) {
		    d__1 = (doublereal) srrt;
		    sig = pow_dd(&d__1, &c_b942) * .0119f;
		} else {
		    d__1 = (doublereal) srrt;
		    sig = pow_dd(&d__1, &c_b943) * .013f;
		}
/* Computing 2nd power */
		r__1 = *srt;
/* Computing 2nd power */
		r__2 = *srt;
/* Computing 2nd power */
		r__3 = *srt;
/* Computing 2nd power */
		r__4 = *em1 + *em2;
/* Computing 2nd power */
		r__5 = *srt;
/* Computing 2nd power */
		r__6 = *em1 - *em2;
		xphid = sig * .33333333333333331f * (r__1 * r__1 - 
			5.0715039999999991f) * (r__2 * r__2 - 
			.044943999999999984f) / (r__3 * r__3 - r__4 * r__4) / 
			(r__5 * r__5 - r__6 * r__6);
	    }
	}
    }
/*   !! in fm^2 */
    *xphin /= 10.f;
/*   !! in fm^2 */
    xphid /= 10.f;
    *xphi = *xphin + xphid;
    return 0;
} /* pibphi_ */


/* **************************** */
/* purpose: Xsection for phi +M to K+K etc */
/* Subroutine */ int phimes_(integer *i1, integer *i2, real *srt, real *xsk1, 
	real *xsk2, real *xsk3, real *xsk4, real *xsk5, real *xsk6, real *
	xsk7, real *sigphi)
{
    /* System generated locals */
    real r__1, r__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static real s;
    static integer lb1, lb2;
    static real em1, em2, pff, pii, srr, srr1, srr2, akap, srrt, scheck;

    /* Fortran I/O blocks */
    static cilist io___2223 = { 0, 99, 0, 0, 0 };


/*     QUANTITIES:                                                      * */
/*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME* */
/*           SRT      - SQRT OF S                                       * */
/*           IBLOCK   - THE INFORMATION BACK                            * */
/*                      223 --> phi destruction */
/*                      20 -->  elastic */
/* ********************************* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /DD/ */
/* c      SAVE /EE/ */
/* Computing 2nd power */
    r__1 = *srt;
    s = r__1 * r__1;
    *sigphi = 1e-8f;
    *xsk1 = 0.f;
    *xsk2 = 0.f;
    *xsk3 = 0.f;
    *xsk4 = 0.f;
    *xsk5 = 0.f;
    *xsk6 = 0.f;
    *xsk7 = 0.f;
    em1 = cc_1.e[*i1 - 1];
    em2 = cc_1.e[*i2 - 1];
    lb1 = ee_1.lb[*i1 - 1];
    lb2 = ee_1.lb[*i2 - 1];
    akap = .498f;
/* ****** */

/*   !! mb, elastic */
    *xsk1 = 5.f;
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__1 = em1 + em2;
/* Computing 2nd power */
    r__2 = em1 - em2;
    scheck = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (scheck <= 0.f) {
	s_wsle(&io___2223);
	do_lio(&c__9, &c__1, "scheck48: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	s_stop("", (ftnlen)0);
    }
    pii = sqrt(scheck);
/*           pii = sqrt((S-(em1+em2)**2)*(S-(em1-em2)**2)) */
/* phi + K(-bar) channel */
    if (lb1 == 23 || lb2 == 23 || lb1 == 21 || lb2 == 21) {
	if (*srt > akap + .13496f) {
/*             XSK2 = 2.5 */
/* Computing 2nd power */
	    r__1 = akap + .13496f;
/* Computing 2nd power */
	    r__2 = .13496f - akap;
	    pff = sqrt((s - r__1 * r__1) * (s - r__2 * r__2));
	    *xsk2 = pff * 195.639f / pii / 32.f / 3.1415926f / s;
	}
	if (*srt > akap + .77f) {
/*              XSK3 = 3.5 */
/* Computing 2nd power */
	    r__1 = akap + .77f;
/* Computing 2nd power */
	    r__2 = .77f - akap;
	    pff = sqrt((s - r__1 * r__1) * (s - r__2 * r__2));
	    *xsk3 = pff * 526.702f / pii / 32.f / 3.1415926f / s;
	}
	if (*srt > akap + .7819f) {
/*               XSK4 = 3.5 */
/* Computing 2nd power */
	    r__1 = akap + .7819f;
/* Computing 2nd power */
	    r__2 = .7819f - akap;
	    pff = sqrt((s - r__1 * r__1) * (s - r__2 * r__2));
	    *xsk4 = pff * 355.429f / pii / 32.f / 3.1415926f / s;
	}
	if (*srt > 1.02996f) {
/*           XSK5 = 15.0 */
	    pff = sqrt((s - 1.0608176015999999f) * (s - .57766080160000011f));
	    *xsk5 = pff * 2047.042f / pii / 32.f / 3.1415926f / s;
	}
	if (*srt > 1.665f) {
/*            XSK6 = 3.5 */
	    pff = sqrt((s - 2.7722250000000002f) * (s - .015625f));
	    *xsk6 = pff * 1371.257f / pii / 32.f / 3.1415926f / s;
	}
	if (*srt > 1.6769000000000001f) {
/*            XSK7 = 3.5 */
	    pff = sqrt((s - 2.81199361f) * (s - .012791609999999995f));
	    *xsk7 = pff * 482.292f / pii / 32.f / 3.1415926f / s;
	}

    } else if (abs(lb1) == 30 || abs(lb2) == 30) {
/* phi + K*(-bar) channel */

	if (*srt > akap + .13496f) {
/*             XSK2 = 3.5 */
/* Computing 2nd power */
	    r__1 = akap + .13496f;
/* Computing 2nd power */
	    r__2 = .13496f - akap;
	    pff = sqrt((s - r__1 * r__1) * (s - r__2 * r__2));
	    *xsk2 = pff * 372.378f / pii / 32.f / 3.1415926f / s;
	}
	if (*srt > akap + .77f) {
/*              XSK3 = 9.0 */
/* Computing 2nd power */
	    r__1 = akap + .77f;
/* Computing 2nd power */
	    r__2 = .77f - akap;
	    pff = sqrt((s - r__1 * r__1) * (s - r__2 * r__2));
	    *xsk3 = pff * 1313.96f / pii / 32.f / 3.1415926f / s;
	}
	if (*srt > akap + .7819f) {
/*               XSK4 = 6.5 */
/* Computing 2nd power */
	    r__1 = akap + .7819f;
/* Computing 2nd power */
	    r__2 = .7819f - akap;
	    pff = sqrt((s - r__1 * r__1) * (s - r__2 * r__2));
	    *xsk4 = pff * 440.558f / pii / 32.f / 3.1415926f / s;
	}
	if (*srt > 1.02996f) {
/*           XSK5 = 30.0 !wrong */
	    pff = sqrt((s - 1.0608176015999999f) * (s - .57766080160000011f));
	    *xsk5 = pff * 1496.692f / pii / 32.f / 3.1415926f / s;
	}
	if (*srt > 1.665f) {
/*            XSK6 = 9.0 */
	    pff = sqrt((s - 2.7722250000000002f) * (s - .015625f));
	    *xsk6 = pff * 6999.84f / pii / 32.f / 3.1415926f / s;
	}
	if (*srt > 1.6769000000000001f) {
/*            XSK7 = 15.0 */
	    pff = sqrt((s - 2.81199361f) * (s - .012791609999999995f));
	    *xsk7 = pff * 1698.903f / pii / 32.f / 3.1415926f / s;
	}
    } else {

/* phi + rho(pi,omega) channel */

	srr1 = em1 + em2;
	if (*srt > akap + akap) {
	    srrt = *srt - srr1;
/* c          if(srrt .lt. 0.3)then */
	    if (srrt < .3f && srrt > .01f) {
		d__1 = (doublereal) srrt;
		*xsk2 = 1.69f / (pow_dd(&d__1, &c_b783) - .407f);
	    } else {
		d__1 = (doublereal) srrt;
		*xsk2 = pow_dd(&d__1, &c_b784) * .008f + 3.74f;
	    }
	}
	if (*srt > akap + .895f) {
	    srr2 = akap + .895f;
	    srr = dmax(srr1,srr2);
	    srrt = *srt - srr;
/* c          if(srrt .lt. 0.3)then */
	    if (srrt < .3f && srrt > .01f) {
		d__1 = (doublereal) srrt;
		*xsk3 = 1.69f / (pow_dd(&d__1, &c_b783) - .407f);
	    } else {
		d__1 = (doublereal) srrt;
		*xsk3 = pow_dd(&d__1, &c_b784) * .008f + 3.74f;
	    }
	}
	if (*srt > 1.79f) {
	    srr2 = 1.79f;
	    srr = dmax(srr1,srr2);
	    srrt = *srt - srr;
/* c          if(srrt .lt. 0.3)then */
	    if (srrt < .3f && srrt > .01f) {
		d__1 = (doublereal) srrt;
		*xsk4 = 1.69f / (pow_dd(&d__1, &c_b783) - .407f);
	    } else {
		d__1 = (doublereal) srrt;
		*xsk4 = pow_dd(&d__1, &c_b784) * .008f + 3.74f;
	    }
	}
/*          xsk2 = amin1(20.,xsk2) */
/*          xsk3 = amin1(20.,xsk3) */
/*          xsk4 = amin1(20.,xsk4) */
    }
    *sigphi = *xsk1 + *xsk2 + *xsk3 + *xsk4 + *xsk5 + *xsk6 + *xsk7;
    return 0;
} /* phimes_ */

/* ********************************* */
/*     PURPOSE:                                                         * */
/*             DEALING WITH phi+M  scatt. */

/* Subroutine */ int crphim_(real *px, real *py, real *pz, real *srt, integer 
	*i1, integer *i2, real *xsk1, real *xsk2, real *xsk3, real *xsk4, 
	real *xsk5, real *xsk6, real *sigphi, integer *ikkg, integer *ikkl, 
	integer *iblock)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    double sqrt(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static real c1, s1, t1, x1, pr;
    static integer lb1, lb2;
    static real em1, em2, ct1, pr2, px0, py0, pz0, st1;
    static integer iad1, iad2;
    extern doublereal ranart_(integer *);
    extern /* Subroutine */ int rotate_(real *, real *, real *, real *, real *
	    , real *);


/*     QUANTITIES:                                                      * */
/*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME* */
/*           SRT      - SQRT OF S                                       * */
/*           IBLOCK   - THE INFORMATION BACK                            * */
/*                      20 -->  elastic */
/*                      223 --> phi + pi(rho,omega) */
/*                      224 --> phi + K -> K + pi(rho,omega) */
/*                      225 --> phi + K -> K* + pi(rho,omega) */
/*                      226 --> phi + K* -> K + pi(rho,omega) */
/*                      227 --> phi + K* -> K* + pi(rho,omega) */
/* ********************************* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /input1/ */
/* c      SAVE /RNDF77/ */

    px0 = *px;
    py0 = *py;
    pz0 = *pz;
    lb1 = ee_1.lb[*i1 - 1];
    lb2 = ee_1.lb[*i2 - 1];
    x1 = ranart_(&rndf77_1.nseed) * *sigphi;
    *xsk2 = *xsk1 + *xsk2;
    *xsk3 = *xsk2 + *xsk3;
    *xsk4 = *xsk3 + *xsk4;
    *xsk5 = *xsk4 + *xsk5;
    *xsk6 = *xsk5 + *xsk6;
    if (x1 <= *xsk1) {
/*        !! elastic scatt */
	*iblock = 20;
	goto L100;
    } else {

/* phi + (K,K*)-bar */
	if (lb1 == 23 || lb1 == 21 || abs(lb1) == 30 || lb2 == 23 || lb2 == 
		21 || abs(lb2) == 30) {

	    if (lb1 == 23 || lb2 == 23) {
		*ikkl = 1;
		*iblock = 224;
		iad1 = 23;
		iad2 = 30;
	    } else if (lb1 == 30 || lb2 == 30) {
		*ikkl = 0;
		*iblock = 226;
		iad1 = 23;
		iad2 = 30;
	    } else if (lb1 == 21 || lb2 == 21) {
		*ikkl = 1;
		*iblock = 124;
		iad1 = 21;
		iad2 = -30;
/*         !! -30 */
	    } else {
		*ikkl = 0;
		*iblock = 126;
		iad1 = 21;
		iad2 = -30;
	    }
	    if (x1 <= *xsk2) {
		ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 
			3;
		ee_1.lb[*i2 - 1] = iad1;
		cc_1.e[*i1 - 1] = .13496f;
		cc_1.e[*i2 - 1] = .498f;
		*ikkg = 1;
		goto L100;
	    } else if (x1 <= *xsk3) {
		ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 
			25;
		ee_1.lb[*i2 - 1] = iad1;
		cc_1.e[*i1 - 1] = .77f;
		cc_1.e[*i2 - 1] = .498f;
		*ikkg = 1;
		goto L100;
	    } else if (x1 <= *xsk4) {
		ee_1.lb[*i1 - 1] = 28;
		ee_1.lb[*i2 - 1] = iad1;
		cc_1.e[*i1 - 1] = .7819f;
		cc_1.e[*i2 - 1] = .498f;
		*ikkg = 1;
		goto L100;
	    } else if (x1 <= *xsk5) {
		ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 
			3;
		ee_1.lb[*i2 - 1] = iad2;
		cc_1.e[*i1 - 1] = .13496f;
		cc_1.e[*i2 - 1] = .895f;
		*ikkg = 0;
		++(*iblock);
		goto L100;
	    } else if (x1 <= *xsk6) {
		ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 
			25;
		ee_1.lb[*i2 - 1] = iad2;
		cc_1.e[*i1 - 1] = .77f;
		cc_1.e[*i2 - 1] = .895f;
		*ikkg = 0;
		++(*iblock);
		goto L100;
	    } else {
		ee_1.lb[*i1 - 1] = 28;
		ee_1.lb[*i2 - 1] = iad2;
		cc_1.e[*i1 - 1] = .7819f;
		cc_1.e[*i2 - 1] = .895f;
		*ikkg = 0;
		++(*iblock);
		goto L100;
	    }
	} else {
/*      !! phi destruction via (pi,rho,omega) */
	    *iblock = 223;
/* phi + pi(rho,omega) */
	    if (x1 <= *xsk2) {
		ee_1.lb[*i1 - 1] = 23;
		ee_1.lb[*i2 - 1] = 21;
		cc_1.e[*i1 - 1] = .498f;
		cc_1.e[*i2 - 1] = .498f;
		*ikkg = 2;
		*ikkl = 0;
		goto L100;
	    } else if (x1 <= *xsk3) {
		ee_1.lb[*i1 - 1] = 23;
/*           LB(I2) = 30 */
		ee_1.lb[*i2 - 1] = -30;
/* lin-2/10/03 currently take XSK3 to be the sum of KK*bar & KbarK*: */
		if (ranart_(&rndf77_1.nseed) <= .5f) {
		    ee_1.lb[*i1 - 1] = 21;
		    ee_1.lb[*i2 - 1] = 30;
		}
		cc_1.e[*i1 - 1] = .498f;
		cc_1.e[*i2 - 1] = .895f;
		*ikkg = 1;
		*ikkl = 0;
		goto L100;
	    } else if (x1 <= *xsk4) {
		ee_1.lb[*i1 - 1] = 30;
/*           LB(I2) = 30 */
		ee_1.lb[*i2 - 1] = -30;
		cc_1.e[*i1 - 1] = .895f;
		cc_1.e[*i2 - 1] = .895f;
		*ikkg = 0;
		*ikkl = 0;
		goto L100;
	    }
	}
    }

L100:
    em1 = cc_1.e[*i1 - 1];
    em2 = cc_1.e[*i2 - 1];
/* ----------------------------------------------------------------------- */
/* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH */
/* ENERGY CONSERVATION */
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = em1;
/* Computing 2nd power */
    r__4 = em2;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = em1 * em2;
    pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
    if (pr2 <= 0.f) {
	pr2 = 1e-8f;
    }
    pr = sqrt(pr2) / (*srt * 2.f);
/* WE ASSUME AN ISOTROPIC ANGULAR DISTRIBUTION IN THE CMS */
    c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
    t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
/* Computing 2nd power */
    r__1 = c1;
    s1 = sqrt(1.f - r__1 * r__1);
    ct1 = cos(t1);
    st1 = sin(t1);
/* THE MOMENTUM IN THE CMS IN THE FINAL STATE */
    *pz = pr * c1;
    *px = pr * s1 * ct1;
    *py = pr * s1 * st1;
/* ROTATE IT */
    rotate_(&px0, &py0, &pz0, px, py, pz);
    return 0;
} /* crphim_ */

/* ********************************* */
/* ********************************* */
/* bz3/9/99 khyperon */
/* ************************************ */
/* purpose: Xsection for K+Y ->  piN                                       * */
/*          Xsection for K+Y-bar ->  piN-bar   !! sp03/29/01               * */

/* Subroutine */ int xkhype_(integer *i1, integer *i2, real *srt, real *xky1, 
	real *xky2, real *xky3, real *xky4, real *xky5, real *xky6, real *
	xky7, real *xky8, real *xky9, real *xky10, real *xky11, real *xky12, 
	real *xky13, real *xky14, real *xky15, real *xky16, real *xky17, real 
	*sigk)
{
    /* System generated locals */
    real r__1, r__2;

    /* Local variables */
    static real s;
    static integer lb1, lb2;
    static real pf2, pi2, xm1, xm2, ddf, sig, srrt;
    extern doublereal pnlka_(real *), pnska_(real *);
    static real xkaon0;

/*      subroutine xkhype(i1, i2, srt, sigk) */
/*  srt    = DSQRT(s) in GeV                                               * */
/*  xkkpi   = xsection in mb obtained from                                 * */
/*           the detailed balance                                          * */
/* *********************************** */
/* c      SAVE /EE/ */
/* Computing 2nd power */
    r__1 = *srt;
    s = r__1 * r__1;
    *sigk = 1e-8f;
    *xky1 = 0.f;
    *xky2 = 0.f;
    *xky3 = 0.f;
    *xky4 = 0.f;
    *xky5 = 0.f;
    *xky6 = 0.f;
    *xky7 = 0.f;
    *xky8 = 0.f;
    *xky9 = 0.f;
    *xky10 = 0.f;
    *xky11 = 0.f;
    *xky12 = 0.f;
    *xky13 = 0.f;
    *xky14 = 0.f;
    *xky15 = 0.f;
    *xky16 = 0.f;
    *xky17 = 0.f;
    lb1 = ee_1.lb[*i1 - 1];
    lb2 = ee_1.lb[*i2 - 1];
    if (abs(lb1) == 14 || abs(lb2) == 14) {
	xkaon0 = pnlka_(srt);
	xkaon0 *= 2.f;
	pi2 = (s - 2.6049960000000003f) * (s - .38192400000000015f);
    } else {
	xkaon0 = pnska_(srt);
	xkaon0 *= 2.f;
	pi2 = (s - 2.8594810000000002f) * (s - .48302500000000009f);
    }
    if (pi2 <= 0.f) {
	return 0;
    }
    xm1 = .14f;
    xm2 = .93828f;
/* Computing 2nd power */
    r__1 = xm1 + xm2;
/* Computing 2nd power */
    r__2 = xm1 - xm2;
    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (pf2 > 0.f) {
	*xky1 = pf2 * 3.f / pi2 * xkaon0;
    }
    xm1 = .14f;
    xm2 = 1.232f;
/* Computing 2nd power */
    r__1 = xm1 + xm2;
/* Computing 2nd power */
    r__2 = xm1 - xm2;
    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (pf2 > 0.f) {
	*xky2 = pf2 * 12.f / pi2 * xkaon0;
    }
    xm1 = .14f;
    xm2 = 1.44f;
/* Computing 2nd power */
    r__1 = xm1 + xm2;
/* Computing 2nd power */
    r__2 = xm1 - xm2;
    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (pf2 > 0.f) {
	*xky3 = pf2 * 3.f / pi2 * xkaon0;
    }
    xm1 = .14f;
    xm2 = 1.535f;
/* Computing 2nd power */
    r__1 = xm1 + xm2;
/* Computing 2nd power */
    r__2 = xm1 - xm2;
    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (pf2 > 0.f) {
	*xky4 = pf2 * 3.f / pi2 * xkaon0;
    }
    xm1 = .769f;
    xm2 = .93828f;
/* Computing 2nd power */
    r__1 = xm1 + xm2;
/* Computing 2nd power */
    r__2 = xm1 - xm2;
    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (pf2 > 0.f) {
	*xky5 = pf2 * 9.f / pi2 * xkaon0;
    }
    xm1 = .769f;
    xm2 = 1.232f;
/* Computing 2nd power */
    r__1 = xm1 + xm2;
/* Computing 2nd power */
    r__2 = xm1 - xm2;
    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (pf2 > 0.f) {
	*xky6 = pf2 * 36.f / pi2 * xkaon0;
    }
    xm1 = .769f;
    xm2 = 1.44f;
/* Computing 2nd power */
    r__1 = xm1 + xm2;
/* Computing 2nd power */
    r__2 = xm1 - xm2;
    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (pf2 > 0.f) {
	*xky7 = pf2 * 9.f / pi2 * xkaon0;
    }
    xm1 = .769f;
    xm2 = 1.535f;
/* Computing 2nd power */
    r__1 = xm1 + xm2;
/* Computing 2nd power */
    r__2 = xm1 - xm2;
    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (pf2 > 0.f) {
	*xky8 = pf2 * 9.f / pi2 * xkaon0;
    }
    xm1 = .782f;
    xm2 = .93828f;
/* Computing 2nd power */
    r__1 = xm1 + xm2;
/* Computing 2nd power */
    r__2 = xm1 - xm2;
    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (pf2 > 0.f) {
	*xky9 = pf2 * 3.f / pi2 * xkaon0;
    }
    xm1 = .782f;
    xm2 = 1.232f;
/* Computing 2nd power */
    r__1 = xm1 + xm2;
/* Computing 2nd power */
    r__2 = xm1 - xm2;
    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (pf2 > 0.f) {
	*xky10 = pf2 * 12.f / pi2 * xkaon0;
    }
    xm1 = .782f;
    xm2 = 1.44f;
/* Computing 2nd power */
    r__1 = xm1 + xm2;
/* Computing 2nd power */
    r__2 = xm1 - xm2;
    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (pf2 > 0.f) {
	*xky11 = pf2 * 3.f / pi2 * xkaon0;
    }
    xm1 = .782f;
    xm2 = 1.535f;
/* Computing 2nd power */
    r__1 = xm1 + xm2;
/* Computing 2nd power */
    r__2 = xm1 - xm2;
    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (pf2 > 0.f) {
	*xky12 = pf2 * 3.f / pi2 * xkaon0;
    }
    xm1 = .5473f;
    xm2 = .93828f;
/* Computing 2nd power */
    r__1 = xm1 + xm2;
/* Computing 2nd power */
    r__2 = xm1 - xm2;
    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (pf2 > 0.f) {
	*xky13 = pf2 * 1.f / pi2 * xkaon0;
    }
    xm1 = .5473f;
    xm2 = 1.232f;
/* Computing 2nd power */
    r__1 = xm1 + xm2;
/* Computing 2nd power */
    r__2 = xm1 - xm2;
    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (pf2 > 0.f) {
	*xky14 = pf2 * 4.f / pi2 * xkaon0;
    }
    xm1 = .5473f;
    xm2 = 1.44f;
/* Computing 2nd power */
    r__1 = xm1 + xm2;
/* Computing 2nd power */
    r__2 = xm1 - xm2;
    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (pf2 > 0.f) {
	*xky15 = pf2 * 1.f / pi2 * xkaon0;
    }
    xm1 = .5473f;
    xm2 = 1.535f;
/* Computing 2nd power */
    r__1 = xm1 + xm2;
/* Computing 2nd power */
    r__2 = xm1 - xm2;
    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (pf2 > 0.f) {
	*xky16 = pf2 * 1.f / pi2 * xkaon0;
    }
/* sp11/21/01  K+ + La --> phi + N */
    if (lb1 == 14 || lb2 == 14) {
	if (*srt > 1.959457f) {
	    srrt = *srt - 1.959457f;
/* Computing 2nd power */
	    r__1 = srrt + 3.508f;
	    sig = 1.715f / (r__1 * r__1 - 12.138f);
	    xm1 = .939457f;
	    xm2 = 1.02f;
/* Computing 2nd power */
	    r__1 = xm1 + xm2;
/* Computing 2nd power */
	    r__2 = xm1 - xm2;
	    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
/*     ! fm^-1 */
	    *xky17 = pf2 * 3.f / pi2 * sig / 10.f;
	}
    }
/* sp11/21/01  end */

    if (abs(lb1) >= 15 && abs(lb1) <= 17 || abs(lb2) >= 15 && abs(lb2) <= 17) 
	    {
	ddf = 3.f;
	*xky1 /= ddf;
	*xky2 /= ddf;
	*xky3 /= ddf;
	*xky4 /= ddf;
	*xky5 /= ddf;
	*xky6 /= ddf;
	*xky7 /= ddf;
	*xky8 /= ddf;
	*xky9 /= ddf;
	*xky10 /= ddf;
	*xky11 /= ddf;
	*xky12 /= ddf;
	*xky13 /= ddf;
	*xky14 /= ddf;
	*xky15 /= ddf;
	*xky16 /= ddf;
    }
    *sigk = *xky1 + *xky2 + *xky3 + *xky4 + *xky5 + *xky6 + *xky7 + *xky8 + *
	    xky9 + *xky10 + *xky11 + *xky12 + *xky13 + *xky14 + *xky15 + *
	    xky16 + *xky17;
    return 0;
} /* xkhype_ */

/* ******************************* */
/* Subroutine */ int ppbdat_(void)
{
    return 0;
} /* ppbdat_ */

/*     to give default values to parameters for BbarB production from mesons */
/* c      SAVE /ppbmas/ */
/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
/*     thresh(i) gives the mass thresh for final channel i: */
/*     ppbm(i,j=1,2) gives masses for the two final baryons of channel i, */
/*     with j=1 for the lighter baryon: */
/*     factr2(i) gives weights for producing i pions from ppbar annihilation: */
/*     niso(i) gives the degeneracy factor for final channel i: */

/* **************************************** */
/* get the number of BbarB states available for mm collisions of energy srt */
/* Subroutine */ int getnst_(real *srt)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2, r__3;

    /* Local variables */
    static integer i__;
    static real pf2;

/*  srt    = DSQRT(s) in GeV                                                   * */
/* **************************************** */
/* c      SAVE /ppbmas/ */
/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
/* Computing 2nd power */
    r__1 = *srt;
    ppb1_1.s = r__1 * r__1;
    ppbmas_1.nstate = 0;
    ppb1_1.wtot = 0.f;
    if (*srt <= ppbmas_1.thresh[0]) {
	return 0;
    }
    for (i__ = 1; i__ <= 15; ++i__) {
	ppbmas_1.weight[i__ - 1] = 0.f;
	if (*srt > ppbmas_1.thresh[i__ - 1]) {
	    ppbmas_1.nstate = i__;
	}
/* L1001: */
    }
    i__1 = ppbmas_1.nstate;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	r__1 = ppbmas_1.ppbm[i__ - 1] + ppbmas_1.ppbm[i__ + 14];
/* Computing 2nd power */
	r__2 = ppbmas_1.ppbm[i__ - 1] - ppbmas_1.ppbm[i__ + 14];
	pf2 = (ppb1_1.s - r__1 * r__1) * (ppb1_1.s - r__2 * r__2) / 4 / 
		ppb1_1.s;
	ppbmas_1.weight[i__ - 1] = pf2 * ppbmas_1.niso[i__ - 1];
	ppb1_1.wtot += ppbmas_1.weight[i__ - 1];
/* L1002: */
    }
/* Computing 3rd power */
    r__1 = *srt / .14f;
    ppb1_1.ene = r__1 * (r__1 * r__1) / 59.217624386248559f;
/* Computing 2nd power */
    r__1 = ppb1_1.ene;
/* Computing 3rd power */
    r__2 = ppb1_1.ene;
/* Computing 4th power */
    r__3 = ppb1_1.ene, r__3 *= r__3;
    ppb1_1.fsum = ppb1_1.factr2[1] + ppb1_1.factr2[2] * ppb1_1.ene + 
	    ppb1_1.factr2[3] * (r__1 * r__1) + ppb1_1.factr2[4] * (r__2 * (
	    r__2 * r__2)) + ppb1_1.factr2[5] * (r__3 * r__3);
    return 0;
} /* getnst_ */

/* **************************************** */
/* for pion+pion-->Bbar B                                                      * */
/*      real*4 function ppbbar(srt) */
doublereal ppbbar_(real *srt)
{
    /* System generated locals */
    real ret_val;

    /* Local variables */
    static real pi2, sppb2p;
    extern doublereal xppbar_(real *);

/* **************************************** */
/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
    sppb2p = xppbar_(srt) * ppb1_1.factr2[1] / ppb1_1.fsum;
    pi2 = (ppb1_1.s - .078400000000000011f) / 4;
    ret_val = sppb2p * .44444444444444442f / pi2 * ppb1_1.wtot;
    return ret_val;
} /* ppbbar_ */

/* **************************************** */
/* for pion+rho-->Bbar B                                                      * */
/*      real*4 function prbbar(srt) */
doublereal prbbar_(real *srt)
{
    /* System generated locals */
    real ret_val;

    /* Local variables */
    static real pi2, sppb3p;
    extern doublereal xppbar_(real *);

/* **************************************** */
/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
    sppb3p = xppbar_(srt) * ppb1_1.factr2[2] * ppb1_1.ene / ppb1_1.fsum;
    pi2 = (ppb1_1.s - .82810000000000006f) * (ppb1_1.s - .39690000000000003f) 
	    / 4 / ppb1_1.s;
    ret_val = sppb3p * .14814814814814814f / pi2 * ppb1_1.wtot;
    return ret_val;
} /* prbbar_ */

/* **************************************** */
/* for rho+rho-->Bbar B                                                      * */
/*      real*4 function rrbbar(srt) */
doublereal rrbbar_(real *srt)
{
    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real pi2, sppb4p;
    extern doublereal xppbar_(real *);

/* **************************************** */
/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
/* Computing 2nd power */
    r__1 = ppb1_1.ene;
    sppb4p = xppbar_(srt) * ppb1_1.factr2[3] * (r__1 * r__1) / ppb1_1.fsum;
    pi2 = (ppb1_1.s - 2.3715999999999999f) / 4;
    ret_val = sppb4p / 2 * .049382716049382713f / pi2 * ppb1_1.wtot;
    return ret_val;
} /* rrbbar_ */

/* **************************************** */
/* for pi+omega-->Bbar B                                                      * */
/*      real*4 function pobbar(srt) */
doublereal pobbar_(real *srt)
{
    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real pi2, sppb4p;
    extern doublereal xppbar_(real *);

/* **************************************** */
/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
/* Computing 2nd power */
    r__1 = ppb1_1.ene;
    sppb4p = xppbar_(srt) * ppb1_1.factr2[3] * (r__1 * r__1) / ppb1_1.fsum;
    pi2 = (ppb1_1.s - .85008400000000006f) * (ppb1_1.s - .41216400000000003f) 
	    / 4 / ppb1_1.s;
    ret_val = sppb4p / 2 * .44444444444444442f / pi2 * ppb1_1.wtot;
    return ret_val;
} /* pobbar_ */

/* **************************************** */
/* for rho+omega-->Bbar B                                                      * */
/*      real*4 function robbar(srt) */
doublereal robbar_(real *srt)
{
    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real pi2, sppb5p;
    extern doublereal xppbar_(real *);

/* **************************************** */
/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
/* Computing 3rd power */
    r__1 = ppb1_1.ene;
    sppb5p = xppbar_(srt) * ppb1_1.factr2[4] * (r__1 * (r__1 * r__1)) / 
	    ppb1_1.fsum;
    pi2 = (ppb1_1.s - 2.4087040000000002f) * (ppb1_1.s - 
	    1.4400000000000025e-4f) / 4 / ppb1_1.s;
    ret_val = sppb5p * .14814814814814814f / pi2 * ppb1_1.wtot;
    return ret_val;
} /* robbar_ */

/* **************************************** */
/* for omega+omega-->Bbar B                                                    * */
/*      real*4 function oobbar(srt) */
doublereal oobbar_(real *srt)
{
    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real pi2, sppb6p;
    extern doublereal xppbar_(real *);

/* **************************************** */
/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
/* Computing 4th power */
    r__1 = ppb1_1.ene, r__1 *= r__1;
    sppb6p = xppbar_(srt) * ppb1_1.factr2[5] * (r__1 * r__1) / ppb1_1.fsum;
    pi2 = (ppb1_1.s - 2.4460960000000003f) / 4;
    ret_val = sppb6p * .44444444444444442f / pi2 * ppb1_1.wtot;
    return ret_val;
} /* oobbar_ */

/* **************************************** */
/* Generate final states for mm-->Bbar B                                       * */
/* Subroutine */ int bbarfs_(integer *lbb1, integer *lbb2, real *ei1, real *
	ei2, integer *iblock, integer *iseed)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static real rd, rd1, rd2;
    static integer ifs;
    static real wsum;
    extern doublereal ranart_(integer *);

/* **************************************** */
/* c      SAVE /ppbmas/ */
/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
/* c      SAVE /RNDF77/ */
/*     determine which final BbarB channel occurs: */
    rd = ranart_(&rndf77_1.nseed);
    wsum = 0.f;
    i__1 = ppbmas_1.nstate;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wsum += ppbmas_1.weight[i__ - 1];
	if (rd <= wsum / ppb1_1.wtot) {
	    ifs = i__;
	    *ei1 = ppbmas_1.ppbm[i__ - 1];
	    *ei2 = ppbmas_1.ppbm[i__ + 14];
	    goto L10;
	}
/* L1001: */
    }
L10:
/* 1    pbar p */
    if (ifs == 1) {
	*iblock = 1801;
	*lbb1 = -1;
	*lbb2 = 1;
    } else if (ifs == 2) {
/* 2    pbar n */
	if (ranart_(&rndf77_1.nseed) <= .5f) {
	    *iblock = 18021;
	    *lbb1 = -1;
	    *lbb2 = 2;
/* 2    nbar p */
	} else {
	    *iblock = 18022;
	    *lbb1 = 1;
	    *lbb2 = -2;
	}
/* 3    nbar n */
    } else if (ifs == 3) {
	*iblock = 1803;
	*lbb1 = -2;
	*lbb2 = 2;
/* 4&5  (pbar nbar) Delta, (p n) anti-Delta */
    } else if (ifs == 4 || ifs == 5) {
	rd = ranart_(&rndf77_1.nseed);
	if (rd <= .5f) {
/*     (pbar nbar) Delta */
	    if (ifs == 4) {
		*iblock = 18041;
		*lbb1 = -1;
	    } else {
		*iblock = 18051;
		*lbb1 = -2;
	    }
	    rd2 = ranart_(&rndf77_1.nseed);
	    if (rd2 <= .25f) {
		*lbb2 = 6;
	    } else if (rd2 <= .5f) {
		*lbb2 = 7;
	    } else if (rd2 <= .75f) {
		*lbb2 = 8;
	    } else {
		*lbb2 = 9;
	    }
	} else {
/*     (p n) anti-Delta */
	    if (ifs == 4) {
		*iblock = 18042;
		*lbb1 = 1;
	    } else {
		*iblock = 18052;
		*lbb1 = 2;
	    }
	    rd2 = ranart_(&rndf77_1.nseed);
	    if (rd2 <= .25f) {
		*lbb2 = -6;
	    } else if (rd2 <= .5f) {
		*lbb2 = -7;
	    } else if (rd2 <= .75f) {
		*lbb2 = -8;
	    } else {
		*lbb2 = -9;
	    }
	}
/* 6&7  (pbar nbar) N*(1440), (p n) anti-N*(1440) */
    } else if (ifs == 6 || ifs == 7) {
	rd = ranart_(&rndf77_1.nseed);
	if (rd <= .5f) {
/*     (pbar nbar) N*(1440) */
	    if (ifs == 6) {
		*iblock = 18061;
		*lbb1 = -1;
	    } else {
		*iblock = 18071;
		*lbb1 = -2;
	    }
	    rd2 = ranart_(&rndf77_1.nseed);
	    if (rd2 <= .5f) {
		*lbb2 = 10;
	    } else {
		*lbb2 = 11;
	    }
	} else {
/*     (p n) anti-N*(1440) */
	    if (ifs == 6) {
		*iblock = 18062;
		*lbb1 = 1;
	    } else {
		*iblock = 18072;
		*lbb1 = 2;
	    }
	    rd2 = ranart_(&rndf77_1.nseed);
	    if (rd2 <= .5f) {
		*lbb2 = -10;
	    } else {
		*lbb2 = -11;
	    }
	}
/* 8    Delta anti-Delta */
    } else if (ifs == 8) {
	*iblock = 1808;
	rd1 = ranart_(&rndf77_1.nseed);
	if (rd1 <= .25f) {
	    *lbb1 = 6;
	} else if (rd1 <= .5f) {
	    *lbb1 = 7;
	} else if (rd1 <= .75f) {
	    *lbb1 = 8;
	} else {
	    *lbb1 = 9;
	}
	rd2 = ranart_(&rndf77_1.nseed);
	if (rd2 <= .25f) {
	    *lbb2 = -6;
	} else if (rd2 <= .5f) {
	    *lbb2 = -7;
	} else if (rd2 <= .75f) {
	    *lbb2 = -8;
	} else {
	    *lbb2 = -9;
	}
/* 9&10 (pbar nbar) N*(1535), (p n) anti-N*(1535) */
    } else if (ifs == 9 || ifs == 10) {
	rd = ranart_(&rndf77_1.nseed);
	if (rd <= .5f) {
/*     (pbar nbar) N*(1440) */
	    if (ifs == 9) {
		*iblock = 18091;
		*lbb1 = -1;
	    } else {
		*iblock = 18101;
		*lbb1 = -2;
	    }
	    rd2 = ranart_(&rndf77_1.nseed);
	    if (rd2 <= .5f) {
		*lbb2 = 12;
	    } else {
		*lbb2 = 13;
	    }
	} else {
/*     (p n) anti-N*(1535) */
	    if (ifs == 9) {
		*iblock = 18092;
		*lbb1 = 1;
	    } else {
		*iblock = 18102;
		*lbb1 = 2;
	    }
	    rd2 = ranart_(&rndf77_1.nseed);
	    if (rd2 <= .5f) {
		*lbb2 = -12;
	    } else {
		*lbb2 = -13;
	    }
	}
/* 11&12 anti-Delta N*, Delta anti-N* */
    } else if (ifs == 11 || ifs == 12) {
	rd = ranart_(&rndf77_1.nseed);
	if (rd <= .5f) {
/*     anti-Delta N* */
	    rd1 = ranart_(&rndf77_1.nseed);
	    if (rd1 <= .25f) {
		*lbb1 = -6;
	    } else if (rd1 <= .5f) {
		*lbb1 = -7;
	    } else if (rd1 <= .75f) {
		*lbb1 = -8;
	    } else {
		*lbb1 = -9;
	    }
	    if (ifs == 11) {
		*iblock = 18111;
		rd2 = ranart_(&rndf77_1.nseed);
		if (rd2 <= .5f) {
		    *lbb2 = 10;
		} else {
		    *lbb2 = 11;
		}
	    } else {
		*iblock = 18121;
		rd2 = ranart_(&rndf77_1.nseed);
		if (rd2 <= .5f) {
		    *lbb2 = 12;
		} else {
		    *lbb2 = 13;
		}
	    }
	} else {
/*     Delta anti-N* */
	    rd1 = ranart_(&rndf77_1.nseed);
	    if (rd1 <= .25f) {
		*lbb1 = 6;
	    } else if (rd1 <= .5f) {
		*lbb1 = 7;
	    } else if (rd1 <= .75f) {
		*lbb1 = 8;
	    } else {
		*lbb1 = 9;
	    }
	    if (ifs == 11) {
		*iblock = 18112;
		rd2 = ranart_(&rndf77_1.nseed);
		if (rd2 <= .5f) {
		    *lbb2 = -10;
		} else {
		    *lbb2 = -11;
		}
	    } else {
		*iblock = 18122;
		rd2 = ranart_(&rndf77_1.nseed);
		if (rd2 <= .5f) {
		    *lbb2 = -12;
		} else {
		    *lbb2 = -13;
		}
	    }
	}
/* 13   N*(1440) anti-N*(1440) */
    } else if (ifs == 13) {
	*iblock = 1813;
	rd1 = ranart_(&rndf77_1.nseed);
	if (rd1 <= .5f) {
	    *lbb1 = 10;
	} else {
	    *lbb1 = 11;
	}
	rd2 = ranart_(&rndf77_1.nseed);
	if (rd2 <= .5f) {
	    *lbb2 = -10;
	} else {
	    *lbb2 = -11;
	}
/* 14   anti-N*(1440) N*(1535), N*(1440) anti-N*(1535) */
    } else if (ifs == 14) {
	rd = ranart_(&rndf77_1.nseed);
	if (rd <= .5f) {
/*     anti-N*(1440) N*(1535) */
	    *iblock = 18141;
	    rd1 = ranart_(&rndf77_1.nseed);
	    if (rd1 <= .5f) {
		*lbb1 = -10;
	    } else {
		*lbb1 = -11;
	    }
	    rd2 = ranart_(&rndf77_1.nseed);
	    if (rd2 <= .5f) {
		*lbb2 = 12;
	    } else {
		*lbb2 = 13;
	    }
	} else {
/*     N*(1440) anti-N*(1535) */
	    *iblock = 18142;
	    rd1 = ranart_(&rndf77_1.nseed);
	    if (rd1 <= .5f) {
		*lbb1 = 10;
	    } else {
		*lbb1 = 11;
	    }
	    rd2 = ranart_(&rndf77_1.nseed);
	    if (rd2 <= .5f) {
		*lbb2 = -12;
	    } else {
		*lbb2 = -13;
	    }
	}
/* 15   N*(1535) anti-N*(1535) */
    } else if (ifs == 15) {
	*iblock = 1815;
	rd1 = ranart_(&rndf77_1.nseed);
	if (rd1 <= .5f) {
	    *lbb1 = 12;
	} else {
	    *lbb1 = 13;
	}
	rd2 = ranart_(&rndf77_1.nseed);
	if (rd2 <= .5f) {
	    *lbb2 = -12;
	} else {
	    *lbb2 = -13;
	}
    } else {
    }
    return 0;
} /* bbarfs_ */

/* **************************************** */
/* for pi pi <-> rho rho cross sections */
/* Subroutine */ int spprr_(integer *lb1, integer *lb2, real *srt)
{
    extern doublereal ptor_(real *), rtop_(real *);

/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
    ppmm_1.pprr = 0.f;
    if (*lb1 >= 3 && *lb1 <= 5 && (*lb2 >= 3 && *lb2 <= 5)) {
/*     for now, rho mass taken to be the central value in these two processes */
	if (*srt > 1.54f) {
	    ppmm_1.pprr = ptor_(srt);
	}
    } else if (*lb1 >= 25 && *lb1 <= 27 && (*lb2 >= 25 && *lb2 <= 27)) {
	ppmm_1.pprr = rtop_(srt);
    }

    return 0;
} /* spprr_ */

/* **************************************** */
/* for pi pi -> rho rho, determined from detailed balance */
doublereal ptor_(real *srt)
{
    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real s2;
    extern doublereal rtop_(real *);

/* **************************************** */
/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
/* Computing 2nd power */
    r__1 = *srt;
    s2 = r__1 * r__1;
    ret_val = (s2 - 2.3715999999999999f) * 9 / (s2 - .078400000000000011f) * 
	    rtop_(srt);
    return ret_val;
} /* ptor_ */

/* **************************************** */
/* for rho rho -> pi pi, assumed a constant cross section (in mb) */
doublereal rtop_(real *srt)
{
    /* System generated locals */
    real ret_val;

/* **************************************** */
    ret_val = 5.f;
    return ret_val;
} /* rtop_ */

/* **************************************** */
/* for pi pi <-> rho rho final states */
/* Subroutine */ int pi2ro2_(integer *i1, integer *i2, integer *lbb1, integer 
	*lbb2, real *ei1, real *ei2, integer *iblock, integer *iseed)
{
    extern doublereal ranart_(integer *);

/* c      SAVE /EE/ */
/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
/* c      SAVE /RNDF77/ */
    if (ee_1.lb[*i1 - 1] >= 3 && ee_1.lb[*i1 - 1] <= 5 && (ee_1.lb[*i2 - 1] >=
	     3 && ee_1.lb[*i2 - 1] <= 5)) {
	*iblock = 1850;
	*ei1 = .77f;
	*ei2 = .77f;
/*     for now, we don't check isospin states(allowing pi+pi+ & pi0pi0 -> 2rho) */
/*     thus the cross sections used are considered as the isospin-averaged ones. */
	*lbb1 = (integer) (ranart_(&rndf77_1.nseed) * 3) + 25;
	*lbb2 = (integer) (ranart_(&rndf77_1.nseed) * 3) + 25;
    } else if (ee_1.lb[*i1 - 1] >= 25 && ee_1.lb[*i1 - 1] <= 27 && (ee_1.lb[*
	    i2 - 1] >= 25 && ee_1.lb[*i2 - 1] <= 27)) {
	*iblock = 1851;
	*lbb1 = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	*lbb2 = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	*ei1 = .13957f;
	*ei2 = .13957f;
	if (*lbb1 == 4) {
	    *ei1 = .13496f;
	}
	if (*lbb2 == 4) {
	    *ei2 = .13496f;
	}
    }
    return 0;
} /* pi2ro2_ */

/* **************************************** */
/* for pi pi <-> eta eta cross sections */
/* Subroutine */ int sppee_(integer *lb1, integer *lb2, real *srt)
{
    extern doublereal ptoe_(real *), etop_(real *);

/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
    ppmm_1.ppee = 0.f;
    if (*lb1 >= 3 && *lb1 <= 5 && (*lb2 >= 3 && *lb2 <= 5)) {
	if (*srt > 1.095f) {
	    ppmm_1.ppee = ptoe_(srt);
	}
    } else if (*lb1 == 0 && *lb2 == 0) {
	ppmm_1.ppee = etop_(srt);
    }
    return 0;
} /* sppee_ */

/* **************************************** */
/* for pi pi -> eta eta, determined from detailed balance, spin-isospin averaged */
doublereal ptoe_(real *srt)
{
    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real s2;
    extern doublereal etop_(real *);

/* **************************************** */
/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
/* Computing 2nd power */
    r__1 = *srt;
    s2 = r__1 * r__1;
    ret_val = (s2 - 1.199025f) * .1111111111111111f / (s2 - 
	    .078400000000000011f) * etop_(srt);
    return ret_val;
} /* ptoe_ */

/* **************************************** */
/* for eta eta -> pi pi, assumed a constant cross section (in mb) */
doublereal etop_(real *srt)
{
    /* System generated locals */
    real ret_val;

/* **************************************** */
/*     eta equilibration: */
/*     most important channel is found to be pi pi <-> pi eta, then */
/*     rho pi <-> rho eta. */
    ret_val = 5.f;
    return ret_val;
} /* etop_ */

/* **************************************** */
/* for pi pi <-> eta eta final states */
/* Subroutine */ int pi2et2_(integer *i1, integer *i2, integer *lbb1, integer 
	*lbb2, real *ei1, real *ei2, integer *iblock, integer *iseed)
{
    extern doublereal ranart_(integer *);

/* c      SAVE /EE/ */
/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
/* c      SAVE /RNDF77/ */
    if (ee_1.lb[*i1 - 1] >= 3 && ee_1.lb[*i1 - 1] <= 5 && (ee_1.lb[*i2 - 1] >=
	     3 && ee_1.lb[*i2 - 1] <= 5)) {
	*iblock = 1860;
	*ei1 = .5475f;
	*ei2 = .5475f;
/*     for now, we don't check isospin states(allowing pi+pi+ & pi0pi0 -> 2rho) */
/*     thus the cross sections used are considered as the isospin-averaged ones. */
	*lbb1 = 0;
	*lbb2 = 0;
    } else if (ee_1.lb[*i1 - 1] == 0 && ee_1.lb[*i2 - 1] == 0) {
	*iblock = 1861;
	*lbb1 = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	*lbb2 = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	*ei1 = .13957f;
	*ei2 = .13957f;
	if (*lbb1 == 4) {
	    *ei1 = .13496f;
	}
	if (*lbb2 == 4) {
	    *ei2 = .13496f;
	}
    }
    return 0;
} /* pi2et2_ */

/* **************************************** */
/* for pi pi <-> pi eta cross sections */
/* Subroutine */ int spppe_(integer *lb1, integer *lb2, real *srt)
{
    extern doublereal pptope_(real *), petopp_(real *);

/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
    ppmm_1.pppe = 0.f;
    if (*lb1 >= 3 && *lb1 <= 5 && (*lb2 >= 3 && *lb2 <= 5)) {
	if (*srt > .6875f) {
	    ppmm_1.pppe = pptope_(srt);
	}
    } else if (*lb1 >= 3 && *lb1 <= 5 && *lb2 == 0) {
	ppmm_1.pppe = petopp_(srt);
    } else if (*lb2 >= 3 && *lb2 <= 5 && *lb1 == 0) {
	ppmm_1.pppe = petopp_(srt);
    }
    return 0;
} /* spppe_ */

/* **************************************** */
/* for pi pi -> pi eta, determined from detailed balance, spin-isospin averaged */
doublereal pptope_(real *srt)
{
    /* System generated locals */
    real ret_val, r__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real s2, pf2, pi2;
    extern doublereal petopp_(real *);

/* **************************************** */
/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
/* Computing 2nd power */
    r__1 = *srt;
    s2 = r__1 * r__1;
    pf2 = (s2 - .47265625f) * (s2 - .16605624999999999f) / 2 / sqrt(s2);
    pi2 = (s2 - .078400000000000011f) * s2 / 2 / sqrt(s2);
    ret_val = pf2 * .33333333333333331f / pi2 * petopp_(srt);
    return ret_val;
} /* pptope_ */

/* **************************************** */
/* for pi eta -> pi pi, assumed a constant cross section (in mb) */
doublereal petopp_(real *srt)
{
    /* System generated locals */
    real ret_val;

/* **************************************** */
/*     eta equilibration: */
    ret_val = 5.f;
    return ret_val;
} /* petopp_ */

/* **************************************** */
/* for pi pi <-> pi eta final states */
/* Subroutine */ int pi3eta_(integer *i1, integer *i2, integer *lbb1, integer 
	*lbb2, real *ei1, real *ei2, integer *iblock, integer *iseed)
{
    extern doublereal ranart_(integer *);

/* c      SAVE /EE/ */
/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
/* c      SAVE /RNDF77/ */
    if (ee_1.lb[*i1 - 1] >= 3 && ee_1.lb[*i1 - 1] <= 5 && (ee_1.lb[*i2 - 1] >=
	     3 && ee_1.lb[*i2 - 1] <= 5)) {
	*iblock = 1870;
	*ei1 = .13957f;
	*ei2 = .5475f;
/*     for now, we don't check isospin states(allowing pi+pi+ & pi0pi0 -> 2rho) */
/*     thus the cross sections used are considered as the isospin-averaged ones. */
	*lbb1 = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	if (*lbb1 == 4) {
	    *ei1 = .13496f;
	}
	*lbb2 = 0;
    } else if (ee_1.lb[*i1 - 1] >= 3 && ee_1.lb[*i1 - 1] <= 5 && ee_1.lb[*i2 
	    - 1] == 0 || ee_1.lb[*i2 - 1] >= 3 && ee_1.lb[*i2 - 1] <= 5 && 
	    ee_1.lb[*i1 - 1] == 0) {
	*iblock = 1871;
	*lbb1 = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	*lbb2 = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	*ei1 = .13957f;
	*ei2 = .13957f;
	if (*lbb1 == 4) {
	    *ei1 = .13496f;
	}
	if (*lbb2 == 4) {
	    *ei2 = .13496f;
	}
    }
    return 0;
} /* pi3eta_ */

/* **************************************** */
/* for rho pi <-> rho eta cross sections */
/* Subroutine */ int srpre_(integer *lb1, integer *lb2, real *srt)
{
    extern doublereal rptore_(real *), retorp_(real *);

/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
    ppmm_1.rpre = 0.f;
    if (*lb1 >= 25 && *lb1 <= 27 && *lb2 >= 3 && *lb2 <= 5) {
	if (*srt > 1.3174999999999999f) {
	    ppmm_1.rpre = rptore_(srt);
	}
    } else if (*lb2 >= 25 && *lb2 <= 27 && *lb1 >= 3 && *lb1 <= 5) {
	if (*srt > 1.3174999999999999f) {
	    ppmm_1.rpre = rptore_(srt);
	}
    } else if (*lb1 >= 25 && *lb1 <= 27 && *lb2 == 0) {
	if (*srt > .91000000000000003f) {
	    ppmm_1.rpre = retorp_(srt);
	}
    } else if (*lb2 >= 25 && *lb2 <= 27 && *lb1 == 0) {
	if (*srt > .91000000000000003f) {
	    ppmm_1.rpre = retorp_(srt);
	}
    }
    return 0;
} /* srpre_ */

/* **************************************** */
/* for rho pi->rho eta, determined from detailed balance, spin-isospin averaged */
doublereal rptore_(real *srt)
{
    /* System generated locals */
    real ret_val, r__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real s2, pf2, pi2;
    extern doublereal retorp_(real *);

/* **************************************** */
/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
/* Computing 2nd power */
    r__1 = *srt;
    s2 = r__1 * r__1;
    pf2 = (s2 - 1.7358062499999998f) * (s2 - .049506250000000016f) / 2 / sqrt(
	    s2);
    pi2 = (s2 - .82810000000000006f) * (s2 - .39690000000000003f) / 2 / sqrt(
	    s2);
    ret_val = pf2 * .33333333333333331f / pi2 * retorp_(srt);
    return ret_val;
} /* rptore_ */

/* **************************************** */
/* for rho eta -> rho pi, assumed a constant cross section (in mb) */
doublereal retorp_(real *srt)
{
    /* System generated locals */
    real ret_val;

/* **************************************** */
/*     eta equilibration: */
    ret_val = 5.f;
    return ret_val;
} /* retorp_ */

/* **************************************** */
/* for rho pi <-> rho eta final states */
/* Subroutine */ int rpiret_(integer *i1, integer *i2, integer *lbb1, integer 
	*lbb2, real *ei1, real *ei2, integer *iblock, integer *iseed)
{
    extern doublereal ranart_(integer *);

/* c      SAVE /EE/ */
/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
/* c      SAVE /RNDF77/ */
    if (ee_1.lb[*i1 - 1] >= 25 && ee_1.lb[*i1 - 1] <= 27 && ee_1.lb[*i2 - 1] 
	    >= 3 && ee_1.lb[*i2 - 1] <= 5 || ee_1.lb[*i1 - 1] >= 3 && ee_1.lb[
	    *i1 - 1] <= 5 && ee_1.lb[*i2 - 1] >= 25 && ee_1.lb[*i2 - 1] <= 27)
	     {
	*iblock = 1880;
	*ei1 = .77f;
	*ei2 = .5475f;
/*     for now, we don't check isospin states(allowing pi+pi+ & pi0pi0 -> 2rho) */
/*     thus the cross sections used are considered as the isospin-averaged ones. */
	*lbb1 = (integer) (ranart_(&rndf77_1.nseed) * 3) + 25;
	*lbb2 = 0;
    } else if (ee_1.lb[*i1 - 1] >= 25 && ee_1.lb[*i1 - 1] <= 27 && ee_1.lb[*
	    i2 - 1] == 0 || ee_1.lb[*i2 - 1] >= 25 && ee_1.lb[*i2 - 1] <= 27 
	    && ee_1.lb[*i1 - 1] == 0) {
	*iblock = 1881;
	*lbb1 = (integer) (ranart_(&rndf77_1.nseed) * 3) + 25;
	*lbb2 = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	*ei1 = .77f;
	*ei2 = .13957f;
	if (*lbb2 == 4) {
	    *ei2 = .13496f;
	}
    }
    return 0;
} /* rpiret_ */

/* **************************************** */
/* for omega pi <-> omega eta cross sections */
/* Subroutine */ int sopoe_(integer *lb1, integer *lb2, real *srt)
{
    extern doublereal xop2oe_(real *), xoe2op_(real *);

/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
    ppmm_1.xopoe = 0.f;
    if (*lb1 == 28 && *lb2 >= 3 && *lb2 <= 5 || *lb2 == 28 && *lb1 >= 3 && *
	    lb1 <= 5) {
	if (*srt > 1.3294999999999999f) {
	    ppmm_1.xopoe = xop2oe_(srt);
	}
    } else if (*lb1 == 28 && *lb2 == 0 || *lb1 == 0 && *lb2 == 28) {
	if (*srt > 1.3294999999999999f) {
	    ppmm_1.xopoe = xoe2op_(srt);
	}
    }
    return 0;
} /* sopoe_ */

/* **************************************** */
/* for omega pi -> omega eta, */
/*     determined from detailed balance, spin-isospin averaged */
doublereal xop2oe_(real *srt)
{
    /* System generated locals */
    real ret_val, r__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real s2, pf2, pi2;
    extern doublereal xoe2op_(real *);

/* **************************************** */
/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
/* Computing 2nd power */
    r__1 = *srt;
    s2 = r__1 * r__1;
    pf2 = (s2 - 1.7675702499999997f) * (s2 - .054990250000000018f) / 2 / sqrt(
	    s2);
    pi2 = (s2 - .85008400000000006f) * (s2 - .41216400000000003f) / 2 / sqrt(
	    s2);
    ret_val = pf2 * .33333333333333331f / pi2 * xoe2op_(srt);
    return ret_val;
} /* xop2oe_ */

/* **************************************** */
/* for omega eta -> omega pi, assumed a constant cross section (in mb) */
doublereal xoe2op_(real *srt)
{
    /* System generated locals */
    real ret_val;

/* **************************************** */
/*     eta equilibration: */
    ret_val = 5.f;
    return ret_val;
} /* xoe2op_ */

/* **************************************** */
/* for omega pi <-> omega eta final states */
/* Subroutine */ int opioet_(integer *i1, integer *i2, integer *lbb1, integer 
	*lbb2, real *ei1, real *ei2, integer *iblock, integer *iseed)
{
    extern doublereal ranart_(integer *);

/* c      SAVE /EE/ */
/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
/* c      SAVE /RNDF77/ */
    if (ee_1.lb[*i1 - 1] >= 3 && ee_1.lb[*i1 - 1] <= 5 && ee_1.lb[*i2 - 1] == 
	    28 || ee_1.lb[*i2 - 1] >= 3 && ee_1.lb[*i2 - 1] <= 5 && ee_1.lb[*
	    i1 - 1] == 28) {
	*iblock = 1890;
	*ei1 = .782f;
	*ei2 = .5475f;
/*     for now, we don't check isospin states(allowing pi+pi+ & pi0pi0 -> 2rho) */
/*     thus the cross sections used are considered as the isospin-averaged ones. */
	*lbb1 = 28;
	*lbb2 = 0;
    } else if (ee_1.lb[*i1 - 1] == 28 && ee_1.lb[*i2 - 1] == 0 || ee_1.lb[*i1 
	    - 1] == 0 && ee_1.lb[*i2 - 1] == 28) {
	*iblock = 1891;
	*lbb1 = 28;
	*lbb2 = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	*ei1 = .782f;
	*ei2 = .13957f;
	if (*lbb2 == 4) {
	    *ei2 = .13496f;
	}
    }
    return 0;
} /* opioet_ */

/* **************************************** */
/* for rho rho <-> eta eta cross sections */
/* Subroutine */ int srree_(integer *lb1, integer *lb2, real *srt)
{
    extern doublereal rrtoee_(real *), eetorr_(real *);

/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
    ppmm_1.rree = 0.f;
    if (*lb1 >= 25 && *lb1 <= 27 && *lb2 >= 25 && *lb2 <= 27) {
	if (*srt > 1.095f) {
	    ppmm_1.rree = rrtoee_(srt);
	}
    } else if (*lb1 == 0 && *lb2 == 0) {
	if (*srt > 1.54f) {
	    ppmm_1.rree = eetorr_(srt);
	}
    }
    return 0;
} /* srree_ */

/* **************************************** */
/* for eta eta -> rho rho */
/*     determined from detailed balance, spin-isospin averaged */
doublereal eetorr_(real *srt)
{
    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real s2;
    extern doublereal rrtoee_(real *);

/* **************************************** */
/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
/* Computing 2nd power */
    r__1 = *srt;
    s2 = r__1 * r__1;
    ret_val = (s2 - 2.3715999999999999f) * 81.f / (s2 - 1.199025f) * rrtoee_(
	    srt);
    return ret_val;
} /* eetorr_ */

/* **************************************** */
/* for rho rho -> eta eta, assumed a constant cross section (in mb) */
doublereal rrtoee_(real *srt)
{
    /* System generated locals */
    real ret_val;

/* **************************************** */
/*     eta equilibration: */
    ret_val = 5.f;
    return ret_val;
} /* rrtoee_ */

/* **************************************** */
/* for rho rho <-> eta eta final states */
/* Subroutine */ int ro2et2_(integer *i1, integer *i2, integer *lbb1, integer 
	*lbb2, real *ei1, real *ei2, integer *iblock, integer *iseed)
{
    extern doublereal ranart_(integer *);

/* c      SAVE /EE/ */
/* c      SAVE /ppb1/ */
/* c      SAVE /ppmm/ */
/* c      SAVE /RNDF77/ */
    if (ee_1.lb[*i1 - 1] >= 25 && ee_1.lb[*i1 - 1] <= 27 && ee_1.lb[*i2 - 1] 
	    >= 25 && ee_1.lb[*i2 - 1] <= 27) {
	*iblock = 1895;
	*ei1 = .5475f;
	*ei2 = .5475f;
/*     for now, we don't check isospin states(allowing pi+pi+ & pi0pi0 -> 2rho) */
/*     thus the cross sections used are considered as the isospin-averaged ones. */
	*lbb1 = 0;
	*lbb2 = 0;
    } else if (ee_1.lb[*i1 - 1] == 0 && ee_1.lb[*i2 - 1] == 0) {
	*iblock = 1896;
	*lbb1 = (integer) (ranart_(&rndf77_1.nseed) * 3) + 25;
	*lbb2 = (integer) (ranart_(&rndf77_1.nseed) * 3) + 25;
	*ei1 = .77f;
	*ei2 = .77f;
    }
    return 0;
} /* ro2et2_ */

/* **************************** */
/* purpose: Xsection for K* Kbar or K*bar K to pi(eta) rho(omega) */
/* Subroutine */ int xkksan_(integer *i1, integer *i2, real *srt, real *
	sigks1, real *sigks2, real *sigks3, real *sigks4, real *sigk, real *
	prkk)
{
    /* System generated locals */
    real r__1, r__2;

    /* Local variables */
    static real s, pf2, pi2, xm1, xm2, xpion0;

/*  srt    = DSQRT(s) in GeV                                       * */
/*  sigk   = xsection in mb obtained from                          * */
/*           the detailed balance                                  * */
/* *************************** */
/* c      SAVE /CC/ */
/* Computing 2nd power */
    r__1 = *srt;
    s = r__1 * r__1;
    *sigks1 = 1e-8f;
    *sigks2 = 1e-8f;
    *sigks3 = 1e-8f;
    *sigks4 = 1e-8f;
    xpion0 = *prkk;
/* lin note that prkk is for pi (rho omega) -> K* Kbar (AND!) K*bar K: */
    xpion0 /= 2;
/* c */
/*        PI2 = (S - (aks + AKA) ** 2) * (S - (aks - AKA) ** 2) */
/* Computing 2nd power */
    r__1 = cc_1.e[*i1 - 1] + cc_1.e[*i2 - 1];
/* Computing 2nd power */
    r__2 = cc_1.e[*i1 - 1] - cc_1.e[*i2 - 1];
    pi2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    *sigk = 1e-8f;
    if (pi2 <= 0.f) {
	return 0;
    }
    xm1 = .14f;
    xm2 = .77f;
/* Computing 2nd power */
    r__1 = xm1 + xm2;
/* Computing 2nd power */
    r__2 = xm1 - xm2;
    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (pi2 > 0.f && pf2 > 0.f) {
	*sigks1 = pf2 * 6.75f / pi2 * xpion0;
    }
    xm1 = .14f;
    xm2 = .7819f;
/* Computing 2nd power */
    r__1 = xm1 + xm2;
/* Computing 2nd power */
    r__2 = xm1 - xm2;
    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (pi2 > 0.f && pf2 > 0.f) {
	*sigks2 = pf2 * 2.25f / pi2 * xpion0;
    }
    xm1 = .77f;
    xm2 = .5473f;
/* Computing 2nd power */
    r__1 = xm1 + xm2;
/* Computing 2nd power */
    r__2 = xm1 - xm2;
    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (pf2 > 0.f) {
	*sigks3 = pf2 * 2.25f / pi2 * xpion0;
    }
    xm1 = .7819f;
    xm2 = .5473f;
/* Computing 2nd power */
    r__1 = xm1 + xm2;
/* Computing 2nd power */
    r__2 = xm1 - xm2;
    pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (pf2 > 0.f) {
	*sigks4 = pf2 * .75f / pi2 * xpion0;
    }
    *sigk = *sigks1 + *sigks2 + *sigks3 + *sigks4;
    return 0;
} /* xkksan_ */

/* ********************************* */
/*     PURPOSE:                                                         * */
/*     assign final states for KK*bar or K*Kbar --> light mesons */

/*      SUBROUTINE Crkspi(PX,PY,PZ,SRT,I1,I2,IBLOCK) */
/* Subroutine */ int crkspi_(integer *i1, integer *i2, real *xsk1, real *xsk2,
	 real *xsk3, real *xsk4, real *sigk, integer *iblock, integer *lbp1, 
	integer *lbp2, real *emm1, real *emm2)
{
    static real x1;
    extern doublereal ranart_(integer *);

/*             iblock   - 466 */
/* ********************************* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /input1/ */
/* c      SAVE /RNDF77/ */
    *iblock = 466;
/* charges of final state mesons: */
    x1 = ranart_(&rndf77_1.nseed) * *sigk;
    *xsk2 = *xsk1 + *xsk2;
    *xsk3 = *xsk2 + *xsk3;
    *xsk4 = *xsk3 + *xsk4;
    if (x1 <= *xsk1) {
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 25;
	cc_1.e[*i1 - 1] = .13957f;
	cc_1.e[*i2 - 1] = .77f;
    } else if (x1 <= *xsk2) {
	ee_1.lb[*i1 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	ee_1.lb[*i2 - 1] = 28;
	cc_1.e[*i1 - 1] = .13957f;
	cc_1.e[*i2 - 1] = .782f;
    } else if (x1 <= *xsk3) {
	ee_1.lb[*i1 - 1] = 0;
	ee_1.lb[*i2 - 1] = (integer) (ranart_(&rndf77_1.nseed) * 3) + 25;
	cc_1.e[*i1 - 1] = .548f;
	cc_1.e[*i2 - 1] = .77f;
    } else {
	ee_1.lb[*i1 - 1] = 0;
	ee_1.lb[*i2 - 1] = 28;
	cc_1.e[*i1 - 1] = .548f;
	cc_1.e[*i2 - 1] = .782f;
    }
    if (ee_1.lb[*i1 - 1] == 4) {
	cc_1.e[*i1 - 1] = .13496f;
    }
    *lbp1 = ee_1.lb[*i1 - 1];
    *lbp2 = ee_1.lb[*i2 - 1];
    *emm1 = cc_1.e[*i1 - 1];
    *emm2 = cc_1.e[*i2 - 1];
    return 0;
} /* crkspi_ */

/* --------------------------------------------------------------------------- */
/* PURPOSE : CALCULATE THE MASS AND MOMENTUM OF K* RESONANCE */
/*           AFTER PION + KAON COLLISION */
/* clin only here the K* mass may be different from aks=0.895 */
/* Subroutine */ int ksreso_(integer *i1, integer *i2)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static integer i__;
    static doublereal p1, p2, p3, e10, e20;
    static real dm;
    static doublereal scheck;

    /* Fortran I/O blocks */
    static cilist io___2304 = { 0, 99, 0, 0, 0 };
    static cilist io___2305 = { 0, 99, 0, 0, 0 };
    static cilist io___2306 = { 0, 99, 0, 0, 0 };
    static cilist io___2307 = { 0, 99, 0, 0, 0 };


/* lin-9/2012: improve precision for argument in sqrt(): */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /RUN/ */
/* c      SAVE /PA/ */
/* c      SAVE /PB/ */
/* c      SAVE /PC/ */
/* c      SAVE /PD/ */
/* 1. DETERMINE THE MOMENTUM COMPONENT OF THE K* IN THE CMS OF PI-K FRAME */
/*    WE LET I1 TO BE THE K* AND ABSORB I2 */
/* lin-9/2012: improve precision for argument in sqrt(): */
/*        E10=SQRT(E(I1)**2+P(1,I1)**2+P(2,I1)**2+P(3,I1)**2) */
/*        E20=SQRT(E(I2)**2+P(1,I2)**2+P(2,I2)**2+P(3,I2)**2) */
/* Computing 2nd power */
    d__1 = (doublereal) cc_1.e[*i1 - 1];
/* Computing 2nd power */
    d__2 = (doublereal) bb_1.p[*i1 * 3 - 3];
/* Computing 2nd power */
    d__3 = (doublereal) bb_1.p[*i1 * 3 - 2];
/* Computing 2nd power */
    d__4 = (doublereal) bb_1.p[*i1 * 3 - 1];
    e10 = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4);
/* Computing 2nd power */
    d__1 = (doublereal) cc_1.e[*i2 - 1];
/* Computing 2nd power */
    d__2 = (doublereal) bb_1.p[*i2 * 3 - 3];
/* Computing 2nd power */
    d__3 = (doublereal) bb_1.p[*i2 * 3 - 2];
/* Computing 2nd power */
    d__4 = (doublereal) bb_1.p[*i2 * 3 - 1];
    e20 = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4);
    p1 = (doublereal) bb_1.p[*i1 * 3 - 3] + (doublereal) bb_1.p[*i2 * 3 - 3];
    p2 = (doublereal) bb_1.p[*i1 * 3 - 2] + (doublereal) bb_1.p[*i2 * 3 - 2];
    p3 = (doublereal) bb_1.p[*i1 * 3 - 1] + (doublereal) bb_1.p[*i2 * 3 - 1];
    if (ee_1.lb[*i2 - 1] == 21 || ee_1.lb[*i2 - 1] == 23) {
	cc_1.e[*i1 - 1] = 0.f;
	i__ = *i2;
    } else {
	cc_1.e[*i2 - 1] = 0.f;
	i__ = *i1;
    }
    if (ee_1.lb[i__ - 1] == 23) {
	ee_1.lb[i__ - 1] = 30;
    } else if (ee_1.lb[i__ - 1] == 21) {
	ee_1.lb[i__ - 1] = -30;
    }
    bb_1.p[i__ * 3 - 3] = bb_1.p[*i1 * 3 - 3] + bb_1.p[*i2 * 3 - 3];
    bb_1.p[i__ * 3 - 2] = bb_1.p[*i1 * 3 - 2] + bb_1.p[*i2 * 3 - 2];
    bb_1.p[i__ * 3 - 1] = bb_1.p[*i1 * 3 - 1] + bb_1.p[*i2 * 3 - 1];
/* 2. DETERMINE THE MASS OF K* BY USING THE REACTION KINEMATICS */
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    d__1 = e10 + e20;
/* Computing 2nd power */
    d__2 = p1;
/* Computing 2nd power */
    d__3 = p2;
/* Computing 2nd power */
    d__4 = p3;
    scheck = d__1 * d__1 - d__2 * d__2 - d__3 * d__3 - d__4 * d__4;
    if (scheck < 0.) {
	s_wsle(&io___2304);
	do_lio(&c__9, &c__1, "scheck49: ", (ftnlen)10);
	do_lio(&c__5, &c__1, (char *)&scheck, (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___2305);
	do_lio(&c__9, &c__1, "scheck49", (ftnlen)8);
	do_lio(&c__5, &c__1, (char *)&scheck, (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&e10, (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&e20, (ftnlen)sizeof(doublereal));
	do_lio(&c__4, &c__1, (char *)&bb_1.p[i__ * 3 - 3], (ftnlen)sizeof(
		real));
	do_lio(&c__4, &c__1, (char *)&bb_1.p[i__ * 3 - 2], (ftnlen)sizeof(
		real));
	do_lio(&c__4, &c__1, (char *)&bb_1.p[i__ * 3 - 1], (ftnlen)sizeof(
		real));
	e_wsle();
	s_wsle(&io___2306);
	do_lio(&c__9, &c__1, "scheck49-1", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&cc_1.e[*i1 - 1], (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&bb_1.p[*i1 * 3 - 3], (ftnlen)sizeof(
		real));
	do_lio(&c__4, &c__1, (char *)&bb_1.p[*i1 * 3 - 2], (ftnlen)sizeof(
		real));
	do_lio(&c__4, &c__1, (char *)&bb_1.p[*i1 * 3 - 1], (ftnlen)sizeof(
		real));
	e_wsle();
	s_wsle(&io___2307);
	do_lio(&c__9, &c__1, "scheck49-2", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&cc_1.e[*i2 - 1], (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&bb_1.p[*i2 * 3 - 3], (ftnlen)sizeof(
		real));
	do_lio(&c__4, &c__1, (char *)&bb_1.p[*i2 * 3 - 2], (ftnlen)sizeof(
		real));
	do_lio(&c__4, &c__1, (char *)&bb_1.p[*i2 * 3 - 1], (ftnlen)sizeof(
		real));
	e_wsle();
    }
    dm = sqrt((real) scheck);
/*        DM=SQRT((E10+E20)**2-P(1,I)**2-P(2,I)**2-P(3,I)**2) */
    cc_1.e[i__ - 1] = dm;
    return 0;
} /* ksreso_ */

/* -------------------------------------------------------- */
/* ************************************ */
/*                                                                         * */
/* Subroutine */ int pertur_(real *px, real *py, real *pz, real *srt, integer 
	*irun, integer *i1, integer *i2, integer *nt, integer *kp, integer *
	icont)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8;

    /* Builtin functions */
    double sqrt(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static real c1, s1, t1, x1, y1, z1, x2, y2, z2, ec;
    static integer ic;
    static real ds, pr;
    static integer lb1, lb2;
    static real em1, em2, ct1, pr2, px0, py0, pz0, st1;
    static integer idn, idp;
    static real app, pii, pff, sig, dfr, dsr, xpt, ypt, zpt, e1cm, e2cm, acap,
	     akal, akap, alas, asap, cmat, ames, aomp, brpp, ppt11, ppt12, 
	    ppt13, ppt21, ppt22, ppt23, srrt;
    static integer lbpp1, lbpp2;
    static real empp1, prob1, prob2, empp2, sigca, pkaon, sigpe, sigpi, xrand,
	     sigom, p1beta;
    static integer icsbel;
    static real sigcal, sigcas, sigeta;
    extern /* Subroutine */ int distce_(integer *, integer *, real *, real *, 
	    real *, real *, real *, integer *, real *, real *, real *);
    extern doublereal aknpsg_(real *), ranart_(integer *);
    static real sigomm, transf;
    extern /* Subroutine */ int rotate_(real *, real *, real *, real *, real *
	    , real *);

/*                                                                         * */
/*       PURPOSE:   TO PRODUCE CASCADE AND OMEGA PERTURBATIVELY            * */
/* sp 01/03/01 */
/*                   40 cascade- */
/*                  -40 cascade-(bar) */
/*                   41 cascade0 */
/*                  -41 cascade0(bar) */
/*                   45 Omega baryon */
/*                  -45 Omega baryon(bar) */
/*                   44 Di-Omega */
/* ********************************* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /HH/ */
/* c      SAVE /ff/ */
/* c      SAVE /gg/ */
/* c      SAVE /INPUT/ */
/* c      SAVE /NN/ */
/* c      SAVE /PA/ */
/* c      SAVE /PB/ */
/* c      SAVE /PC/ */
/* c      SAVE /PD/ */
/* c      SAVE /PE/ */
/* c      SAVE /RR/ */
/* c      SAVE /BG/ */
/* c      SAVE /input1/ */
/*     perturbative method is disabled: */
/*      common /imulst/ iperts */

/* c      SAVE /RNDF77/ */
    px0 = *px;
    py0 = *py;
    pz0 = *pz;
    lb1 = ee_1.lb[*i1 - 1];
    em1 = cc_1.e[*i1 - 1];
    x1 = aa_1.r__[*i1 * 3 - 3];
    y1 = aa_1.r__[*i1 * 3 - 2];
    z1 = aa_1.r__[*i1 * 3 - 1];
    prob1 = hh_1.proper[*i1 - 1];

    lb2 = ee_1.lb[*i2 - 1];
    em2 = cc_1.e[*i2 - 1];
    x2 = aa_1.r__[*i2 * 3 - 3];
    y2 = aa_1.r__[*i2 * 3 - 2];
    z2 = aa_1.r__[*i2 * 3 - 1];
    prob2 = hh_1.proper[*i2 - 1];

/*                 !! flag for real 2-body process (1/0=no/yes) */
    *icont = 1;
/*                !! flag for elastic scatt only (-1=no) */
    icsbel = -1;
/* K-/K*0bar + La/Si --> cascade + pi */
/* K+/K*0 + La/Si (bar) --> cascade-bar + pi */
    if ((lb1 == 21 || lb1 == 23 || abs(lb1) == 30) && (abs(lb2) >= 14 && abs(
	    lb2) <= 17)) {
	goto L60;
    }
    if ((lb2 == 21 || lb2 == 23 || abs(lb2) == 30) && (abs(lb1) >= 14 && abs(
	    lb1) <= 17)) {
	goto L60;
    }
/* K-/K*0bar + cascade --> omega + pi */
/* K+/K*0 + cascade-bar --> omega-bar + pi */
    if ((lb1 == 21 || lb1 == 23 || abs(lb1) == 30) && (abs(lb2) == 40 || abs(
	    lb2) == 41)) {
	goto L70;
    }
    if ((lb2 == 21 || lb2 == 23 || abs(lb2) == 30) && (abs(lb1) == 40 || abs(
	    lb1) == 41)) {
	goto L70;
    }

/* annhilation of cascade,cascade-bar, omega,omega-bar */

/* K- + La/Si <-- cascade + pi(eta,rho,omega) */
/* K+ + La/Si(bar) <-- cascade-bar + pi(eta,rho,omega) */
    if ((lb1 >= 3 && lb1 <= 5 || lb1 == 0) && (abs(lb2) == 40 || abs(lb2) == 
	    41) || (lb2 >= 3 && lb2 <= 5 || lb2 == 0) && (abs(lb1) == 40 || 
	    abs(lb1) == 41)) {
	goto L90;
    }
/* K- + cascade <-- omega + pi */
/* K+ + cascade-bar <-- omega-bar + pi */
/*         if( (lb1.eq.0.and.iabs(lb2).eq.45) */
/*    &    .OR. (lb2.eq.0.and.iabs(lb1).eq.45) ) go to 110 */
    if (lb1 >= 3 && lb1 <= 5 && abs(lb2) == 45 || lb2 >= 3 && lb2 <= 5 && abs(
	    lb1) == 45) {
	goto L110;
    }

/* ---------------------------------------------------- */
/*  for process:  K-bar + L(S) --> Ca + pi */

L60:
    if (abs(lb1) >= 14 && abs(lb1) <= 17) {
	asap = cc_1.e[*i1 - 1];
	akap = cc_1.e[*i2 - 1];
	idp = *i1;
    } else {
	asap = cc_1.e[*i2 - 1];
	akap = cc_1.e[*i1 - 1];
	idp = *i2;
    }
    app = .138f;
    if (*srt < app + 1.3213f) {
	return 0;
    }
    srrt = *srt - (app + 1.3213f) + (akap + .939457f);
/* Computing 2nd power */
    r__2 = srrt;
/* Computing 2nd power */
    r__3 = akap;
/* Computing 2nd power */
    r__1 = (r__2 * r__2 - (r__3 * r__3 + .88257945484900002f)) / 2.f / 
	    .939457f;
/* Computing 2nd power */
    r__4 = akap;
    pkaon = sqrt(r__1 * r__1 - r__4 * r__4);
    sigca = (aknpsg_(&pkaon) + aknpsg_(&pkaon)) * 1.5f;
/* lin pii & pff should be each divided by (4*srt**2), */
/*     but these two factors cancel out in the ratio pii/pff: */
/* Computing 2nd power */
    r__1 = *srt;
/* Computing 2nd power */
    r__2 = akap + .939457f;
/* Computing 2nd power */
    r__3 = *srt;
/* Computing 2nd power */
    r__4 = .939457f - akap;
    pii = sqrt((r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4));
/* Computing 2nd power */
    r__1 = *srt;
/* Computing 2nd power */
    r__2 = asap + app;
/* Computing 2nd power */
    r__3 = *srt;
/* Computing 2nd power */
    r__4 = asap - app;
    pff = sqrt((r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4));
    cmat = sigca * pii / pff;
/* Computing 2nd power */
    r__1 = *srt;
/* Computing 2nd power */
    r__2 = app + 1.3213f;
/* Computing 2nd power */
    r__3 = *srt;
/* Computing 2nd power */
    r__4 = 1.3213f - app;
/* Computing 2nd power */
    r__5 = *srt;
/* Computing 2nd power */
    r__6 = asap + akap;
/* Computing 2nd power */
    r__7 = *srt;
/* Computing 2nd power */
    r__8 = asap - akap;
    sigpi = cmat * sqrt((r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * 
	    r__4)) / sqrt((r__5 * r__5 - r__6 * r__6) * (r__7 * r__7 - r__8 * 
	    r__8));

    sigeta = 0.f;
    if (*srt > 1.8693f) {
	srrt = *srt - 1.8693f + (akap + .939457f);
/* Computing 2nd power */
	r__2 = srrt;
/* Computing 2nd power */
	r__3 = akap;
/* Computing 2nd power */
	r__1 = (r__2 * r__2 - (r__3 * r__3 + .88257945484900002f)) / 2.f / 
		.939457f;
/* Computing 2nd power */
	r__4 = akap;
	pkaon = sqrt(r__1 * r__1 - r__4 * r__4);
	sigca = (aknpsg_(&pkaon) + aknpsg_(&pkaon)) * 1.5f;
	cmat = sigca * pii / pff;
/* Computing 2nd power */
	r__1 = *srt;
/* Computing 2nd power */
	r__2 = *srt;
/* Computing 2nd power */
	r__3 = *srt;
/* Computing 2nd power */
	r__4 = asap + akap;
/* Computing 2nd power */
	r__5 = *srt;
/* Computing 2nd power */
	r__6 = asap - akap;
	sigeta = cmat * sqrt((r__1 * r__1 - 3.4942824899999998f) * (r__2 * 
		r__2 - .59799288999999978f)) / sqrt((r__3 * r__3 - r__4 * 
		r__4) * (r__5 * r__5 - r__6 * r__6));
    }

    sigca = sigpi + sigeta;
    sigpe = 0.f;
/* lin-2/25/03 disable the perturb option: */
/*        if(iperts .eq. 1) sigpe = 40.   !! perturbative xsecn */
    sig = dmax(sigpe,sigca);
    ds = sqrt(sig / 31.4f);
    dsr = ds + .1f;
/* Computing 2nd power */
    r__1 = em1 + em2 + .02f;
    ec = r__1 * r__1;
    distce_(i1, i2, &dsr, &ds, &input1_1.dt, &ec, srt, &ic, px, py, pz);
    if (ic == -1) {
	return 0;
    }
    brpp = sigca / sig;

/* else particle production */
    if (lb1 >= 14 && lb1 <= 17 || lb2 >= 14 && lb2 <= 17) {
/*   !! cascade- or cascde0 */
	lbpp1 = (integer) (ranart_(&rndf77_1.nseed) * 2) + 40;
    } else {
/* elseif(lb1 .eq. -14 .or. lb2 .eq. -14) */
/*     !! cascade-bar- or cascde0 -bar */
	lbpp1 = -40 - (integer) (ranart_(&rndf77_1.nseed) * 2);
    }
    empp1 = 1.3213f;
    if (ranart_(&rndf77_1.nseed) < sigpi / sigca) {
/*    !! pion */
	lbpp2 = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	empp2 = .138f;
    } else {
/*    !! eta */
	lbpp2 = 0;
	empp2 = .548f;
    }
/* * check real process of cascade(bar) and pion formation */
    if (ranart_(&rndf77_1.nseed) < brpp) {
/*       !! real process flag */
	*icont = 0;
	ee_1.lb[*i1 - 1] = lbpp1;
	cc_1.e[*i1 - 1] = empp1;
/*  !! cascade formed with prob Gam */
	hh_1.proper[*i1 - 1] = brpp;
	ee_1.lb[*i2 - 1] = lbpp2;
	cc_1.e[*i2 - 1] = empp2;
/*         !! pion/eta formed with prob 1. */
	hh_1.proper[*i2 - 1] = 1.f;
    }
/* else only cascade(bar) formed perturbatively */
    goto L700;
/* ---------------------------------------------------- */
/*  for process:  Cas(bar) + K_bar(K) --> Om(bar) + pi  !! eta */

L70:
    if (abs(lb1) == 40 || abs(lb1) == 41) {
	acap = cc_1.e[*i1 - 1];
	akap = cc_1.e[*i2 - 1];
	idp = *i1;
    } else {
	acap = cc_1.e[*i2 - 1];
	akap = cc_1.e[*i1 - 1];
	idp = *i2;
    }
    app = .138f;
/*         ames = aeta */
/*  !! only pion */
    ames = .138f;
    if (*srt < ames + 1.6724f) {
	return 0;
    }
    srrt = *srt - (ames + 1.6724f) + (akap + .939457f);
/* Computing 2nd power */
    r__2 = srrt;
/* Computing 2nd power */
    r__3 = akap;
/* Computing 2nd power */
    r__1 = (r__2 * r__2 - (r__3 * r__3 + .88257945484900002f)) / 2.f / 
	    .939457f;
/* Computing 2nd power */
    r__4 = akap;
    pkaon = sqrt(r__1 * r__1 - r__4 * r__4);
/* use K(bar) + Ca --> Om + eta  xsecn same as  K(bar) + N --> Si + Pi */
/*  as Omega have no resonances */
/* ** using same matrix elements as K-bar + N -> Si + pi */
    sigomm = (aknpsg_(&pkaon) + aknpsg_(&pkaon)) * 1.5f;
/* Computing 2nd power */
    r__1 = *srt;
/* Computing 2nd power */
    r__2 = akap + .939457f;
/* Computing 2nd power */
    r__3 = *srt;
/* Computing 2nd power */
    r__4 = .939457f - akap;
/* Computing 2nd power */
    r__5 = *srt;
/* Computing 2nd power */
    r__6 = app + 1.1974f;
/* Computing 2nd power */
    r__7 = *srt;
/* Computing 2nd power */
    r__8 = 1.1974f - app;
    cmat = sigomm * sqrt((r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * 
	    r__4)) / sqrt((r__5 * r__5 - r__6 * r__6) * (r__7 * r__7 - r__8 * 
	    r__8));
/* Computing 2nd power */
    r__1 = *srt;
/* Computing 2nd power */
    r__2 = ames + 1.6724f;
/* Computing 2nd power */
    r__3 = *srt;
/* Computing 2nd power */
    r__4 = 1.6724f - ames;
/* Computing 2nd power */
    r__5 = *srt;
/* Computing 2nd power */
    r__6 = acap + akap;
/* Computing 2nd power */
    r__7 = *srt;
/* Computing 2nd power */
    r__8 = acap - akap;
    sigom = cmat * sqrt((r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * 
	    r__4)) / sqrt((r__5 * r__5 - r__6 * r__6) * (r__7 * r__7 - r__8 * 
	    r__8));
    sigpe = 0.f;
/* lin-2/25/03 disable the perturb option: */
/*         if(iperts .eq. 1) sigpe = 40.   !! perturbative xsecn */
    sig = dmax(sigpe,sigom);
    ds = sqrt(sig / 31.4f);
    dsr = ds + .1f;
/* Computing 2nd power */
    r__1 = em1 + em2 + .02f;
    ec = r__1 * r__1;
    distce_(i1, i2, &dsr, &ds, &input1_1.dt, &ec, srt, &ic, px, py, pz);
    if (ic == -1) {
	return 0;
    }
    brpp = sigom / sig;

/* else particle production */
    if (lb1 >= 40 && lb1 <= 41 || lb2 >= 40 && lb2 <= 41) {
/*    !! omega */
	lbpp1 = 45;
    } else {
/* elseif(lb1 .eq. -40 .or. lb2 .eq. -40) */
/*    !! omega-bar */
	lbpp1 = -45;
    }
    empp1 = 1.6724f;
/*           lbpp2 = 0    !! eta */
/*    !! pion */
    lbpp2 = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
    empp2 = ames;

/* * check real process of omega(bar) and pion formation */
    xrand = ranart_(&rndf77_1.nseed);
    if (xrand < hh_1.proper[idp - 1] * brpp) {
/*       !! real process flag */
	*icont = 0;
	ee_1.lb[*i1 - 1] = lbpp1;
	cc_1.e[*i1 - 1] = empp1;
/*  !! P_Om = P_Cas*Gam */
	hh_1.proper[*i1 - 1] = hh_1.proper[idp - 1] * brpp;
	ee_1.lb[*i2 - 1] = lbpp2;
	cc_1.e[*i2 - 1] = empp2;
/*   !! pion formed with prob 1. */
	hh_1.proper[*i2 - 1] = 1.f;
    } else if (xrand < brpp) {
/* else omega(bar) formed perturbatively and cascade destroyed */
	cc_1.e[idp - 1] = 0.f;
    }
    goto L700;
/* ----------------------------------------------------------- */
/*  for process:  Ca + pi/eta --> K-bar + L(S) */

L90:
    if (abs(lb1) == 40 || abs(lb1) == 41) {
	acap = cc_1.e[*i1 - 1];
	app = cc_1.e[*i2 - 1];
	idp = *i1;
	idn = *i2;
    } else {
	acap = cc_1.e[*i2 - 1];
	app = cc_1.e[*i1 - 1];
	idp = *i2;
	idn = *i1;
    }
/*            akal = (aka+aks)/2.  !! average of K and K* taken */
/*  !! using K only */
    akal = .498f;

    alas = 1.1157f;
    if (*srt <= alas + .498f) {
	return 0;
    }
    srrt = *srt - (acap + app) + 1.437457f;
/* Computing 2nd power */
    r__2 = srrt;
/* Computing 2nd power */
    r__1 = (r__2 * r__2 - 1.1305834548489999f) / 2.f / .939457f;
    pkaon = sqrt(r__1 * r__1 - .248004f);
/* ** using same matrix elements as K-bar + N -> La/Si + pi */
    sigca = (aknpsg_(&pkaon) + aknpsg_(&pkaon)) * 1.5f;
/* Computing 2nd power */
    r__1 = *srt;
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = *srt;
/* Computing 2nd power */
    r__4 = alas + .138f;
/* Computing 2nd power */
    r__5 = *srt;
/* Computing 2nd power */
    r__6 = alas - .138f;
    cmat = sigca * sqrt((r__1 * r__1 - 2.066282626849f) * (r__2 * r__2 - 
	    .194884282849f)) / sqrt((r__3 * r__3 - r__4 * r__4) * (r__5 * 
	    r__5 - r__6 * r__6));
/* Computing 2nd power */
    r__1 = *srt;
/* Computing 2nd power */
    r__2 = acap + app;
/* Computing 2nd power */
    r__3 = *srt;
/* Computing 2nd power */
    r__4 = acap - app;
/* Computing 2nd power */
    r__5 = *srt;
/* Computing 2nd power */
    r__6 = alas + .498f;
/* Computing 2nd power */
    r__7 = *srt;
/* Computing 2nd power */
    r__8 = alas - .498f;
    sigca = cmat * sqrt((r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * 
	    r__4)) / sqrt((r__5 * r__5 - r__6 * r__6) * (r__7 * r__7 - r__8 * 
	    r__8));
/*    !! pi */
    dfr = .33333333333333331f;
/*       !! eta */
    if (ee_1.lb[idn - 1] == 0) {
	dfr = 1.f;
    }
/* Computing 2nd power */
    r__1 = *srt;
/* Computing 2nd power */
    r__2 = alas + .498f;
/* Computing 2nd power */
    r__3 = *srt;
/* Computing 2nd power */
    r__4 = alas - .498f;
/* Computing 2nd power */
    r__5 = *srt;
/* Computing 2nd power */
    r__6 = acap + app;
/* Computing 2nd power */
    r__7 = *srt;
/* Computing 2nd power */
    r__8 = acap - app;
    sigcal = sigca * dfr * (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 *
	     r__4) / (r__5 * r__5 - r__6 * r__6) / (r__7 * r__7 - r__8 * r__8)
	    ;

    alas = 1.1974f;
    if (*srt <= alas + .498f) {
	sigcas = 0.f;
    } else {
	srrt = *srt - (acap + app) + 1.437457f;
/* Computing 2nd power */
	r__2 = srrt;
/* Computing 2nd power */
	r__1 = (r__2 * r__2 - 1.1305834548489999f) / 2.f / .939457f;
	pkaon = sqrt(r__1 * r__1 - .248004f);
/* use K(bar) + La/Si --> Ca + Pi  xsecn same as  K(bar) + N --> Si + Pi */
/* ** using same matrix elements as K-bar + N -> La/Si + pi */
	sigca = (aknpsg_(&pkaon) + aknpsg_(&pkaon)) * 1.5f;
/* Computing 2nd power */
	r__1 = *srt;
/* Computing 2nd power */
	r__2 = *srt;
/* Computing 2nd power */
	r__3 = *srt;
/* Computing 2nd power */
	r__4 = alas + .138f;
/* Computing 2nd power */
	r__5 = *srt;
/* Computing 2nd power */
	r__6 = alas - .138f;
	cmat = sigca * sqrt((r__1 * r__1 - 2.066282626849f) * (r__2 * r__2 - 
		.194884282849f)) / sqrt((r__3 * r__3 - r__4 * r__4) * (r__5 * 
		r__5 - r__6 * r__6));
/* Computing 2nd power */
	r__1 = *srt;
/* Computing 2nd power */
	r__2 = acap + app;
/* Computing 2nd power */
	r__3 = *srt;
/* Computing 2nd power */
	r__4 = acap - app;
/* Computing 2nd power */
	r__5 = *srt;
/* Computing 2nd power */
	r__6 = alas + .498f;
/* Computing 2nd power */
	r__7 = *srt;
/* Computing 2nd power */
	r__8 = alas - .498f;
	sigca = cmat * sqrt((r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 
		* r__4)) / sqrt((r__5 * r__5 - r__6 * r__6) * (r__7 * r__7 - 
		r__8 * r__8));
/*    !! pi */
	dfr = 1.f;
/*    !! eta */
	if (ee_1.lb[idn - 1] == 0) {
	    dfr = 3.f;
	}
/* Computing 2nd power */
	r__1 = *srt;
/* Computing 2nd power */
	r__2 = alas + .498f;
/* Computing 2nd power */
	r__3 = *srt;
/* Computing 2nd power */
	r__4 = alas - .498f;
/* Computing 2nd power */
	r__5 = *srt;
/* Computing 2nd power */
	r__6 = acap + app;
/* Computing 2nd power */
	r__7 = *srt;
/* Computing 2nd power */
	r__8 = acap - app;
	sigcas = sigca * dfr * (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - 
		r__4 * r__4) / (r__5 * r__5 - r__6 * r__6) / (r__7 * r__7 - 
		r__8 * r__8);
    }

    sig = sigcal + sigcas;
    brpp = 1.f;
    ds = sqrt(sig / 31.4f);
    dsr = ds + .1f;
/* Computing 2nd power */
    r__1 = em1 + em2 + .02f;
    ec = r__1 * r__1;
    distce_(i1, i2, &dsr, &ds, &input1_1.dt, &ec, srt, &ic, px, py, pz);

/* lin-2/25/03: checking elastic scatt after failure of inelastic scatt gives */
/*     conditional probability (in general incorrect), tell Pal to correct: */
    if (ic == -1) {
/* check for elastic scatt, no particle annhilation */
/*  !! elastic cross section of 20 mb */
	ds = sqrt(.63694267515923575f);
	dsr = ds + .1f;
	distce_(i1, i2, &dsr, &ds, &input1_1.dt, &ec, srt, &icsbel, px, py, 
		pz);
	if (icsbel == -1) {
	    return 0;
	}
	empp1 = em1;
	empp2 = em2;
	goto L700;
    }

/* else pert. produced cascade(bar) is annhilated OR real process */

/* DECIDE LAMBDA OR SIGMA PRODUCTION */

    if (sigcal / sig > ranart_(&rndf77_1.nseed)) {
	if (lb1 == 40 || lb1 == 41 || lb2 == 40 || lb2 == 41) {
	    lbpp1 = 21;
	    lbpp2 = 14;
	} else {
	    lbpp1 = 23;
	    lbpp2 = -14;
	}
	alas = 1.1157f;
    } else {
	if (lb1 == 40 || lb1 == 41 || lb2 == 40 || lb2 == 41) {
	    lbpp1 = 21;
	    lbpp2 = (integer) (ranart_(&rndf77_1.nseed) * 3) + 15;
	} else {
	    lbpp1 = 23;
	    lbpp2 = -15 - (integer) (ranart_(&rndf77_1.nseed) * 3);
	}
	alas = 1.1974f;
    }
    empp1 = .498f;
    empp2 = alas;

/* check for real process for L/S(bar) and K(bar) formation */
    if (ranart_(&rndf77_1.nseed) < hh_1.proper[idp - 1]) {
/* real process */
/*       !! real process flag */
	*icont = 0;
	ee_1.lb[*i1 - 1] = lbpp1;
	cc_1.e[*i1 - 1] = empp1;
/*   !! K(bar) formed with prob 1. */
	hh_1.proper[*i1 - 1] = 1.f;
	ee_1.lb[*i2 - 1] = lbpp2;
	cc_1.e[*i2 - 1] = empp2;
/*   !! L/S(bar) formed with prob 1. */
	hh_1.proper[*i2 - 1] = 1.f;
	goto L700;
    } else {
/* else only cascade(bar) annhilation & go out */
	cc_1.e[idp - 1] = 0.f;
    }
    return 0;

/* ---------------------------------------------------- */
/*  for process:  Om(bar) + pi --> Cas(bar) + K_bar(K) */

L110:
    if (lb1 == 45 || lb1 == -45) {
	aomp = cc_1.e[*i1 - 1];
	app = cc_1.e[*i2 - 1];
	idp = *i1;
	idn = *i2;
    } else {
	aomp = cc_1.e[*i2 - 1];
	app = cc_1.e[*i1 - 1];
	idp = *i2;
	idn = *i1;
    }
/*            akal = (aka+aks)/2.  !! average of K and K* taken */
/*  !! using K only */
    akal = .498f;
    if (*srt <= 1.8192999999999999f) {
	return 0;
    }
    srrt = *srt - (app + 1.6724f) + 1.437457f;
/* Computing 2nd power */
    r__2 = srrt;
/* Computing 2nd power */
    r__1 = (r__2 * r__2 - 1.1305834548489999f) / 2.f / .939457f;
    pkaon = sqrt(r__1 * r__1 - .248004f);
/* use K(bar) + Ca --> Om + eta  xsecn same as  K(bar) + N --> Si + Pi */
/* ** using same matrix elements as K-bar + N -> La/Si + pi */
    sigca = (aknpsg_(&pkaon) + aknpsg_(&pkaon)) * 1.5f;
/* Computing 2nd power */
    r__1 = *srt;
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = *srt;
/* Computing 2nd power */
    r__4 = *srt;
    cmat = sigca * sqrt((r__1 * r__1 - 2.066282626849f) * (r__2 * r__2 - 
	    .194884282849f)) / sqrt((r__3 * r__3 - 1.7832931599999997f) * (
	    r__4 * r__4 - 1.1223283600000002f));
/* Computing 2nd power */
    r__1 = *srt;
/* Computing 2nd power */
    r__2 = aomp + app;
/* Computing 2nd power */
    r__3 = *srt;
/* Computing 2nd power */
    r__4 = aomp - app;
/* Computing 2nd power */
    r__5 = *srt;
/* Computing 2nd power */
    r__6 = *srt;
    sigom = cmat * sqrt((r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * 
	    r__4)) / sqrt((r__5 * r__5 - 3.3098524899999995f) * (r__6 * r__6 
	    - .67782288999999984f));
/*            dfr = 2.    !! eta */
/*    !! pion */
    dfr = .66666666666666663f;
/* Computing 2nd power */
    r__1 = *srt;
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = *srt;
/* Computing 2nd power */
    r__4 = aomp + app;
/* Computing 2nd power */
    r__5 = *srt;
/* Computing 2nd power */
    r__6 = aomp - app;
    sigom = sigom * dfr * (r__1 * r__1 - 3.3098524899999995f) * (r__2 * r__2 
	    - .67782288999999984f) / (r__3 * r__3 - r__4 * r__4) / (r__5 * 
	    r__5 - r__6 * r__6);

    brpp = 1.f;
    ds = sqrt(sigom / 31.4f);
    dsr = ds + .1f;
/* Computing 2nd power */
    r__1 = em1 + em2 + .02f;
    ec = r__1 * r__1;
    distce_(i1, i2, &dsr, &ds, &input1_1.dt, &ec, srt, &ic, px, py, pz);

/* lin-2/25/03: checking elastic scatt after failure of inelastic scatt gives */
/*     conditional probability (in general incorrect), tell Pal to correct: */
    if (ic == -1) {
/* check for elastic scatt, no particle annhilation */
/*  !! elastic cross section of 20 mb */
	ds = sqrt(.63694267515923575f);
	dsr = ds + .1f;
	distce_(i1, i2, &dsr, &ds, &input1_1.dt, &ec, srt, &icsbel, px, py, 
		pz);
	if (icsbel == -1) {
	    return 0;
	}
	empp1 = em1;
	empp2 = em2;
	goto L700;
    }

/* else pert. produced omega(bar) annhilated  OR real process */
/* annhilate only pert. omega, rest from hijing go out WITHOUT annhil. */
    if (lb1 == 45 || lb2 == 45) {
/*  !! Ca */
	lbpp1 = (integer) (ranart_(&rndf77_1.nseed) * 2) + 40;
/*   !! K- */
	lbpp2 = 21;
    } else {
/* elseif(lb1 .eq. -45 .or. lb2 .eq. -45) */
/*    !! Ca-bar */
	lbpp1 = -40 - (integer) (ranart_(&rndf77_1.nseed) * 2);
/*      !! K+ */
	lbpp2 = 23;
    }
    empp1 = 1.3213f;
    empp2 = .498f;

/* check for real process for Cas(bar) and K(bar) formation */
    if (ranart_(&rndf77_1.nseed) < hh_1.proper[idp - 1]) {
/*       !! real process flag */
	*icont = 0;
	ee_1.lb[*i1 - 1] = lbpp1;
	cc_1.e[*i1 - 1] = empp1;
/*   !! P_Cas(bar) = P_Om(bar) */
	hh_1.proper[*i1 - 1] = hh_1.proper[idp - 1];
	ee_1.lb[*i2 - 1] = lbpp2;
	cc_1.e[*i2 - 1] = empp2;
/*   !! K(bar) formed with prob 1. */
	hh_1.proper[*i2 - 1] = 1.f;

    } else {
/* else Cascade(bar)  produced and Omega(bar) annhilated */
	cc_1.e[idp - 1] = 0.f;
    }
/*   !! for produced particles */
    goto L700;

/* ----------------------------------------------------------- */
L700:
/* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH */
/* ENERGY CONSERVATION */
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = empp1;
/* Computing 2nd power */
    r__4 = empp2;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = empp1 * empp2;
    pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
    if (pr2 <= 0.f) {
	pr2 = 1e-8f;
    }
    pr = sqrt(pr2) / (*srt * 2.f);
/* using isotropic */
    c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
    t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
/* Computing 2nd power */
    r__1 = c1;
    s1 = sqrt(1.f - r__1 * r__1);
    ct1 = cos(t1);
    st1 = sin(t1);
/* THE MOMENTUM IN THE CMS IN THE FINAL STATE */
    *pz = pr * c1;
    *px = pr * s1 * ct1;
    *py = pr * s1 * st1;
/* ROTATE IT */
    rotate_(&px0, &py0, &pz0, px, py, pz);
    if (*icont == 0) {
	return 0;
    }

/* LORENTZ-TRANSFORMATION INTO CMS FRAME */
/* Computing 2nd power */
    r__1 = empp1;
/* Computing 2nd power */
    r__2 = *px;
/* Computing 2nd power */
    r__3 = *py;
/* Computing 2nd power */
    r__4 = *pz;
    e1cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    p1beta = *px * bg_1.betax + *py * bg_1.betay + *pz * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) + e1cm);
    ppt11 = bg_1.betax * transf + *px;
    ppt12 = bg_1.betay * transf + *py;
    ppt13 = bg_1.betaz * transf + *pz;

/* c** for elastic scattering update the momentum of pertb particles */
    if (icsbel != -1) {
/*            if(EMpp1 .gt. 0.9)then */
	bb_1.p[*i1 * 3 - 3] = ppt11;
	bb_1.p[*i1 * 3 - 2] = ppt12;
	bb_1.p[*i1 * 3 - 1] = ppt13;
/*            else */
/* Computing 2nd power */
	r__1 = empp2;
/* Computing 2nd power */
	r__2 = *px;
/* Computing 2nd power */
	r__3 = *py;
/* Computing 2nd power */
	r__4 = *pz;
	e2cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	transf = bg_1.gamma * (-bg_1.gamma * p1beta / (bg_1.gamma + 1) + e2cm)
		;
	ppt21 = bg_1.betax * transf - *px;
	ppt22 = bg_1.betay * transf - *py;
	ppt23 = bg_1.betaz * transf - *pz;
	bb_1.p[*i2 * 3 - 3] = ppt21;
	bb_1.p[*i2 * 3 - 2] = ppt22;
	bb_1.p[*i2 * 3 - 1] = ppt23;
/*            endif */
	return 0;
    }
/* lin-5/2008: */
/* 2008        X01 = 1.0 - 2.0 * RANART(NSEED) */
/*            Y01 = 1.0 - 2.0 * RANART(NSEED) */
/*            Z01 = 1.0 - 2.0 * RANART(NSEED) */
/*        IF ((X01*X01+Y01*Y01+Z01*Z01) .GT. 1.0) GOTO 2008 */
/*                Xpt=X1+0.5*x01 */
/*                Ypt=Y1+0.5*y01 */
/*                Zpt=Z1+0.5*z01 */
    xpt = x1;
    ypt = y1;
    zpt = z1;


/*          if(lbpp1 .eq. 45)then */
/*           write(*,*)'II lb1,lb2,lbpp1,empp1,proper(idp),brpp' */
/*           write(*,*)lb1,lb2,lbpp1,empp1,proper(idp),brpp */
/*          endif */

    ++nn_1.nnn;
    pe_1.propi[nn_1.nnn + *irun * 150001 - 150002] = hh_1.proper[idp - 1] * 
	    brpp;
    pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = lbpp1;
    pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = empp1;
    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = xpt;
    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = ypt;
    pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = zpt;
    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = ppt11;
    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = ppt12;
    pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = ppt13;
/* lin-5/2008: */
    dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = dpert_1.dpertp[*i1 - 
	    1] * dpert_1.dpertp[*i2 - 1];
    return 0;
} /* pertur_ */

/* ********************************* */
/*  sp 12/08/00                                                         * */
/* Subroutine */ int crhb_(real *px, real *py, real *pz, real *srt, integer *
	i1, integer *i2, integer *iblock)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    double sqrt(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static real c1, s1, t1, pr, em1, em2, ct1, pr2, px0, py0, pz0, st1;
    static integer ntag;
    extern doublereal ranart_(integer *);

/*     PURPOSE:                                                         * */
/*        DEALING WITH hyperon+N(D,N*)->hyp+N(D,N*) elastic PROCESS     * */
/*     NOTE   :                                                         * */

/*     QUANTITIES:                                                 * */
/*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME* */
/*           SRT      - SQRT OF S                                       * */
/*           IBLOCK   - THE INFORMATION BACK                            * */
/*                     144-> hyp+N(D,N*)->hyp+N(D,N*) */
/* ********************************* */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
/* c      SAVE /input1/ */
/* c      SAVE /RNDF77/ */
    px0 = *px;
    py0 = *py;
    pz0 = *pz;
/* ----------------------------------------------------------------------- */
    *iblock = 144;
    ntag = 0;
    em1 = cc_1.e[*i1 - 1];
    em2 = cc_1.e[*i2 - 1];
/* ----------------------------------------------------------------------- */
/* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH */
/* ENERGY CONSERVATION */
/* Computing 2nd power */
    r__2 = *srt;
/* Computing 2nd power */
    r__3 = em1;
/* Computing 2nd power */
    r__4 = em2;
/* Computing 2nd power */
    r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
/* Computing 2nd power */
    r__5 = em1 * em2;
    pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
    if (pr2 <= 0.f) {
	pr2 = 1e-9f;
    }
    pr = sqrt(pr2) / (*srt * 2.f);
    c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
    t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
/* Computing 2nd power */
    r__1 = c1;
    s1 = sqrt(1.f - r__1 * r__1);
    ct1 = cos(t1);
    st1 = sin(t1);
    *pz = pr * c1;
    *px = pr * s1 * ct1;
    *py = pr * s1 * st1;
    return 0;
} /* crhb_ */

/* *************************************** */
/* sp 04/05/01 */
/* Purpose: lambda-baryon elastic xsection as a functon of their cms energy */
/* Subroutine */ int lambar_(integer *i1, integer *i2, real *srt, real *
	siglab)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2, r__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static real emb, eml, plab, pthr, plab2;

/*  srt    = DSQRT(s) in GeV                                               * */
/*  siglab = lambda-nuclar elastic cross section in mb */
/*         = 12 + 0.43/p_lab**3.3 (mb) */

/* (2) Calculate p(lab) from srt [GeV], since the formular in the */
/* reference applies only to the case of a p_bar on a proton at rest */
/* Formula used: srt**2=2.*pmass*(pmass+sqrt(pmass**2+plab**2)) */
/* **************************** */
/* c      SAVE /AA/ */
/* c      SAVE /BB/ */
/* c      SAVE /CC/ */
/* c      SAVE /EE/ */
    *siglab = 1e-6f;
    if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) >= 14 && (i__2 = ee_1.lb[*i1 - 1]
	    , abs(i__2)) <= 17) {
	eml = cc_1.e[*i1 - 1];
	emb = cc_1.e[*i2 - 1];
    } else {
	eml = cc_1.e[*i2 - 1];
	emb = cc_1.e[*i1 - 1];
    }
/* Computing 2nd power */
    r__1 = *srt;
/* Computing 2nd power */
    r__2 = eml;
/* Computing 2nd power */
    r__3 = emb;
    pthr = r__1 * r__1 - r__2 * r__2 - r__3 * r__3;
    if (pthr > 0.f) {
/* Computing 2nd power */
	r__1 = pthr / 2.f / emb;
/* Computing 2nd power */
	r__2 = eml;
	plab2 = r__1 * r__1 - r__2 * r__2;
	if (plab2 > 0.f) {
	    plab = sqrt(plab2);
	    d__1 = (doublereal) plab;
	    *siglab = .43f / pow_dd(&d__1, &c_b1071) + 12.f;
	    if (*siglab > 200.f) {
		*siglab = 200.f;
	    }
	}
    }
    return 0;
} /* lambar_ */

/* ------------------------------------------------------------------ */
/* lin-7/26/03 improve speed */
/* ************************************** */
/* Subroutine */ int distc0_(real *drmax, real *deltr0, real *dt, integer *
	ifirst, real *px1cm, real *py1cm, real *pz1cm, real *x1, real *y1, 
	real *z1, real *px1, real *py1, real *pz1, real *em1, real *x2, real *
	y2, real *z2, real *px2, real *py2, real *pz2, real *em2)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real e1, e2, bbb, ddd, dzz, drcm, dxcm, dycm, dzcm, prcm, p1beta, 
	    drbeta, relvel, transf;

/* PURPOSE : CHECK IF THE COLLISION BETWEEN TWO PARTICLES CAN HAPPEN */
/*           BY CHECKING */
/*                      (2) IF PARTICLE WILL PASS EACH OTHER WITHIN */
/*           TWO HARD CORE RADIUS. */
/*                      (3) IF PARTICLES WILL GET CLOSER. */
/* VARIABLES : */
/*           Ifirst=1 COLLISION may HAPPENED */
/*           Ifirst=-1 COLLISION CAN NOT HAPPEN */
/* **************************************** */
/* c      SAVE /BG/ */
    *ifirst = -1;
/* Computing 2nd power */
    r__1 = *em1;
/* Computing 2nd power */
    r__2 = *px1;
/* Computing 2nd power */
    r__3 = *py1;
/* Computing 2nd power */
    r__4 = *pz1;
    e1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
/* NOW PARTICLES ARE CLOSE ENOUGH TO EACH OTHER ! */
/* Computing 2nd power */
    r__1 = *em2;
/* Computing 2nd power */
    r__2 = *px2;
/* Computing 2nd power */
    r__3 = *py2;
/* Computing 2nd power */
    r__4 = *pz2;
    e2 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
/* NOW THERE IS ENOUGH ENERGY AVAILABLE ! */
/* LORENTZ-TRANSFORMATION IN I1-I2-C.M. SYSTEM */
/* BETAX, BETAY, BETAZ AND GAMMA HAVE BEEN GIVEN IN THE SUBROUTINE CMS */
/* TRANSFORMATION OF MOMENTA (PX1CM = - PX2CM) */
    p1beta = *px1 * bg_1.betax + *py1 * bg_1.betay + *pz1 * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) - e1);
/* Computing 2nd power */
    r__1 = *px1cm;
/* Computing 2nd power */
    r__2 = *py1cm;
/* Computing 2nd power */
    r__3 = *pz1cm;
    prcm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
    if (prcm <= 1e-5f) {
	return 0;
    }
/* TRANSFORMATION OF SPATIAL DISTANCE */
    drbeta = bg_1.betax * (*x1 - *x2) + bg_1.betay * (*y1 - *y2) + bg_1.betaz 
	    * (*z1 - *z2);
    transf = bg_1.gamma * bg_1.gamma * drbeta / (bg_1.gamma + 1);
    dxcm = bg_1.betax * transf + *x1 - *x2;
    dycm = bg_1.betay * transf + *y1 - *y2;
    dzcm = bg_1.betaz * transf + *z1 - *z2;
/* DETERMINING IF THIS IS THE POINT OF CLOSEST APPROACH */
/* Computing 2nd power */
    r__1 = dxcm;
/* Computing 2nd power */
    r__2 = dycm;
/* Computing 2nd power */
    r__3 = dzcm;
    drcm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
    dzz = (*px1cm * dxcm + *py1cm * dycm + *pz1cm * dzcm) / prcm;
/* Computing 2nd power */
    r__1 = drcm;
/* Computing 2nd power */
    r__2 = dzz;
    if (r__1 * r__1 - r__2 * r__2 <= 0.f) {
	bbb = 0.f;
    } else {
/* Computing 2nd power */
	r__1 = drcm;
/* Computing 2nd power */
	r__2 = dzz;
	bbb = sqrt(r__1 * r__1 - r__2 * r__2);
    }
/* WILL PARTICLE PASS EACH OTHER WITHIN 2 * HARD CORE RADIUS ? */
    if (bbb > *drmax) {
	return 0;
    }
    relvel = prcm * (1.f / e1 + 1.f / e2);
    ddd = relvel * *dt * .5f;
/* WILL PARTICLES GET CLOSER ? */
    if (dabs(ddd) < dabs(dzz)) {
	return 0;
    }
    *ifirst = 1;
    return 0;
} /* distc0_ */

/* --------------------------------------------------------------------------- */

/* lin-8/2008 B+B->Deuteron+Meson cross section in mb: */
/* Subroutine */ int sbbdm_(real *srt, real *sdprod, integer *ianti, integer *
	lbm, real *xmm, real *pfinal)
{
    /* System generated locals */
    real r__1, r__2;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    static real pifactor, pinitial, s, sbbdomega, x1, threshold, fs;
    static integer ilb1, ilb2;
    static real snew, scheck, sbbdpi;
    extern doublereal fnndpi_(real *), ranart_(integer *);
    static real sbbdeta, sbbdrho;

    /* Fortran I/O blocks */
    static cilist io___2420 = { 0, 99, 0, 0, 0 };



    *sdprod = 0.f;
    sbbdpi = 0.f;
    sbbdrho = 0.f;
    sbbdomega = 0.f;
    sbbdeta = 0.f;
    if (*srt <= leadng_1.em1 + dpi_1.em2) {
	return 0;
    }

    ilb1 = abs(leadng_1.lb1);
    ilb2 = abs(dpi_1.lb2);
/* test off check Xsec using fixed mass for resonances: */
/*      if(ilb1.ge.6.and.ilb1.le.9) then */
/*         em1=1.232 */
/*      elseif(ilb1.ge.10.and.ilb1.le.11) then */
/*         em1=1.44 */
/*      elseif(ilb1.ge.12.and.ilb1.le.13) then */
/*         em1=1.535 */
/*      endif */
/*      if(ilb2.ge.6.and.ilb2.le.9) then */
/*         em2=1.232 */
/*      elseif(ilb2.ge.10.and.ilb2.le.11) then */
/*         em2=1.44 */
/*      elseif(ilb2.ge.12.and.ilb2.le.13) then */
/*         em2=1.535 */
/*      endif */

/* Computing 2nd power */
    r__1 = *srt;
    s = r__1 * r__1;
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__1 = leadng_1.em1 + dpi_1.em2;
/* Computing 2nd power */
    r__2 = leadng_1.em1 - dpi_1.em2;
    scheck = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (scheck <= 0.f) {
	s_wsle(&io___2420);
	do_lio(&c__9, &c__1, "scheck50: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	s_stop("", (ftnlen)0);
    }
    pinitial = sqrt(scheck) / 2.f / *srt;
/*      pinitial=sqrt((s-(em1+em2)**2)*(s-(em1-em2)**2))/2./srt */
    fs = fnndpi_(&s);
/*     Determine isospin and spin factors for the ratio between */
/*     BB->Deuteron+Meson and Deuteron+Meson->BB cross sections: */
    if (para8_1.idxsec == 1 || para8_1.idxsec == 2) {
/*     Assume B+B -> d+Meson has the same cross sections as N+N -> d+pi: */
    } else {
/*     Assume d+Meson -> B+B has the same cross sections as d+pi -> N+N, */
/*     then determine B+B -> d+Meson cross sections: */
	if (ilb1 >= 1 && ilb1 <= 2 && ilb2 >= 1 && ilb2 <= 2) {
	    pifactor = 1.125f;
	} else if (ilb1 >= 1 && ilb1 <= 2 && ilb2 >= 6 && ilb2 <= 9 || ilb2 >=
		 1 && ilb2 <= 2 && ilb1 >= 6 && ilb1 <= 9) {
	    pifactor = .140625f;
	} else if (ilb1 >= 1 && ilb1 <= 2 && ilb2 >= 10 && ilb2 <= 13 || ilb2 
		>= 1 && ilb2 <= 2 && ilb1 >= 10 && ilb1 <= 13) {
	    pifactor = .5625f;
	} else if (ilb1 >= 6 && ilb1 <= 9 && ilb2 >= 6 && ilb2 <= 9) {
	    pifactor = .0703125f;
	} else if (ilb1 >= 6 && ilb1 <= 9 && ilb2 >= 10 && ilb2 <= 13 || ilb2 
		>= 6 && ilb2 <= 9 && ilb1 >= 10 && ilb1 <= 13) {
	    pifactor = .140625f;
	} else if (ilb1 >= 10 && ilb1 <= 11 && ilb2 >= 10 && ilb2 <= 11 || 
		ilb2 >= 12 && ilb2 <= 13 && ilb1 >= 12 && ilb1 <= 13) {
	    pifactor = 1.125f;
	} else if (ilb1 >= 10 && ilb1 <= 11 && ilb2 >= 12 && ilb2 <= 13 || 
		ilb2 >= 10 && ilb2 <= 11 && ilb1 >= 12 && ilb1 <= 13) {
	    pifactor = .5625f;
	}
    }
/*     d pi: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE */
/*     (1) FOR P+P->Deuteron+pi+: */
    if (ilb1 * ilb2 == 1) {
	*lbm = 5;
	if (*ianti == 1) {
	    *lbm = 3;
	}
	*xmm = .13957f;
/*     (2)FOR N+N->Deuteron+pi-: */
    } else if (ilb1 == 2 && ilb2 == 2) {
	*lbm = 3;
	if (*ianti == 1) {
	    *lbm = 5;
	}
	*xmm = .13957f;
/*     (3)FOR N+P->Deuteron+pi0: */
    } else if (ilb1 * ilb2 == 2) {
	*lbm = 4;
	*xmm = .13496f;
    } else {
/*     For baryon resonances, use isospin-averaged cross sections: */
	*lbm = (integer) (ranart_(&rndf77_1.nseed) * 3) + 3;
	if (*lbm == 4) {
	    *xmm = .13496f;
	} else {
	    *xmm = .13957f;
	}
    }

    if (*srt >= *xmm + 1.8756f) {
/* Computing 2nd power */
	r__1 = *xmm + 1.8756f;
/* Computing 2nd power */
	r__2 = 1.8756f - *xmm;
	*pfinal = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
	if (ilb1 == 1 && ilb2 == 1 || ilb1 == 2 && ilb2 == 2) {
/*     for pp or nn initial states: */
	    sbbdpi = fs * *pfinal / pinitial / 4.f;
	} else if (ilb1 == 1 && ilb2 == 2 || ilb1 == 2 && ilb2 == 1) {
/*     factor of 1/2 for pn or np initial states: */
	    sbbdpi = fs * *pfinal / pinitial / 4.f / 2.f;
	} else {
/*     for other BB initial states (spin- and isospin averaged): */
	    if (para8_1.idxsec == 1) {
/*     1: assume the same |matrix element|**2/s (after averaging over initial */
/*     spins and isospins) for B+B -> deuteron+meson at the same sqrt(s); */
		sbbdpi = fs * *pfinal / pinitial * 3.f / 16.f;
	    } else if (para8_1.idxsec == 2 || para8_1.idxsec == 4) {
/* Computing MAX */
		r__1 = *xmm + 1.8756f, r__2 = leadng_1.em1 + dpi_1.em2;
		threshold = dmax(r__1,r__2);
/* Computing 2nd power */
		r__1 = *srt - threshold + 2.012f;
		snew = r__1 * r__1;
		if (para8_1.idxsec == 2) {
/*     2: assume the same |matrix element|**2/s for B+B -> deuteron+meson */
/*     at the same sqrt(s)-threshold: */
		    sbbdpi = fnndpi_(&snew) * *pfinal / pinitial * 3.f / 16.f;
		} else if (para8_1.idxsec == 4) {
/*     4: assume the same |matrix element|**2/s for B+B <- deuteron+meson */
/*     at the same sqrt(s)-threshold: */
		    sbbdpi = fnndpi_(&snew) * *pfinal / pinitial / 6.f * 
			    pifactor;
		}
	    } else if (para8_1.idxsec == 3) {
/*     3: assume the same |matrix element|**2/s for B+B <- deuteron+meson */
/*     at the same sqrt(s): */
		sbbdpi = fs * *pfinal / pinitial / 6.f * pifactor;
	    }

	}
    }

/*     d rho: DETERMINE THE CROSS SECTION TO THIS FINAL STATE: */
    if (*srt > 2.6456f) {
	*pfinal = sqrt((s - 6.9991993599999995f) * (s - 1.2223513599999998f)) 
		/ 2.f / *srt;
	if (para8_1.idxsec == 1) {
	    sbbdrho = fs * *pfinal / pinitial * 3.f / 16.f;
	} else if (para8_1.idxsec == 2 || para8_1.idxsec == 4) {
/* Computing MAX */
	    r__1 = 2.6456f, r__2 = leadng_1.em1 + dpi_1.em2;
	    threshold = dmax(r__1,r__2);
/* Computing 2nd power */
	    r__1 = *srt - threshold + 2.012f;
	    snew = r__1 * r__1;
	    if (para8_1.idxsec == 2) {
		sbbdrho = fnndpi_(&snew) * *pfinal / pinitial * 3.f / 16.f;
	    } else if (para8_1.idxsec == 4) {
/*     The spin- and isospin-averaged factor is 3-times larger for rho: */
		sbbdrho = fnndpi_(&snew) * *pfinal / pinitial / 6.f * (
			pifactor * 3.f);
	    }
	} else if (para8_1.idxsec == 3) {
	    sbbdrho = fs * *pfinal / pinitial / 6.f * (pifactor * 3.f);
	}
    }

/*     d omega: DETERMINE THE CROSS SECTION TO THIS FINAL STATE: */
    if (*srt > 2.6576f) {
	*pfinal = sqrt((s - 7.0628377599999999f) * (s - 1.1959609599999999f)) 
		/ 2.f / *srt;
	if (para8_1.idxsec == 1) {
	    sbbdomega = fs * *pfinal / pinitial * 3.f / 16.f;
	} else if (para8_1.idxsec == 2 || para8_1.idxsec == 4) {
/* Computing MAX */
	    r__1 = 2.6576f, r__2 = leadng_1.em1 + dpi_1.em2;
	    threshold = dmax(r__1,r__2);
/* Computing 2nd power */
	    r__1 = *srt - threshold + 2.012f;
	    snew = r__1 * r__1;
	    if (para8_1.idxsec == 2) {
		sbbdomega = fnndpi_(&snew) * *pfinal / pinitial * 3.f / 16.f;
	    } else if (para8_1.idxsec == 4) {
		sbbdomega = fnndpi_(&snew) * *pfinal / pinitial / 6.f * 
			pifactor;
	    }
	} else if (para8_1.idxsec == 3) {
	    sbbdomega = fs * *pfinal / pinitial / 6.f * pifactor;
	}
    }

/*     d eta: DETERMINE THE CROSS SECTION TO THIS FINAL STATE: */
    if (*srt > 2.4236f) {
	*pfinal = sqrt((s - 5.8738369600000002f) * (s - 1.7625217599999996f)) 
		/ 2.f / *srt;
	if (para8_1.idxsec == 1) {
	    sbbdeta = fs * *pfinal / pinitial * 3.f / 16.f;
	} else if (para8_1.idxsec == 2 || para8_1.idxsec == 4) {
/* Computing MAX */
	    r__1 = 2.4236f, r__2 = leadng_1.em1 + dpi_1.em2;
	    threshold = dmax(r__1,r__2);
/* Computing 2nd power */
	    r__1 = *srt - threshold + 2.012f;
	    snew = r__1 * r__1;
	    if (para8_1.idxsec == 2) {
		sbbdeta = fnndpi_(&snew) * *pfinal / pinitial * 3.f / 16.f;
	    } else if (para8_1.idxsec == 4) {
		sbbdeta = fnndpi_(&snew) * *pfinal / pinitial / 6.f * (
			pifactor / 3.f);
	    }
	} else if (para8_1.idxsec == 3) {
	    sbbdeta = fs * *pfinal / pinitial / 6.f * (pifactor / 3.f);
	}
    }

    *sdprod = sbbdpi + sbbdrho + sbbdomega + sbbdeta;
/* test off */
/*      write(99,111) srt,sbbdpi,sbbdrho,sbbdomega,sbbdeta,sdprod */
/* 111  format(6(f8.2,1x)) */

    if (*sdprod <= 0.f) {
	return 0;
    }

/*     choose final state and assign masses here: */
    x1 = ranart_(&rndf77_1.nseed);
    if (x1 <= sbbdpi / *sdprod) {
/*     use the above-determined lbm and xmm. */
    } else if (x1 <= (sbbdpi + sbbdrho) / *sdprod) {
	*lbm = (integer) (ranart_(&rndf77_1.nseed) * 3) + 25;
	*xmm = .77f;
    } else if (x1 <= (sbbdpi + sbbdrho + sbbdomega) / *sdprod) {
	*lbm = 28;
	*xmm = .782f;
    } else {
	*lbm = 0;
	*xmm = .548f;
    }

    return 0;
} /* sbbdm_ */


/*     Generate angular distribution of Deuteron in the CMS frame: */
/* Subroutine */ int bbdangle_(real *pxd, real *pyd, real *pzd, integer *nt, 
	integer *ipert1, integer *ianti, integer *idloop, real *pfinal, real *
	dprob1, integer *lbm)
{
    /* System generated locals */
    real r__1;

    /* Builtin functions */
    double sqrt(doublereal), cos(doublereal), sin(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static real c1, s1, t1, ct1, st1, dprob;
    extern doublereal ranart_(integer *);

    /* Fortran I/O blocks */
    static cilist io___2433 = { 0, 91, 0, 0, 0 };
    static cilist io___2434 = { 0, 91, 0, 0, 0 };
    static cilist io___2435 = { 0, 91, 0, 0, 0 };
    static cilist io___2436 = { 0, 91, 0, 0, 0 };


/*     take isotropic distribution for now: */
    c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
    t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
/* Computing 2nd power */
    r__1 = c1;
    s1 = sqrt(1.f - r__1 * r__1);
    ct1 = cos(t1);
    st1 = sin(t1);
/* THE MOMENTUM IN THE CMS IN THE FINAL STATE */
    *pzd = *pfinal * c1;
    *pxd = *pfinal * s1 * ct1;
    *pyd = *pfinal * s1 * st1;
/* lin-5/2008 track the number of produced deuterons: */
    if (para8_1.idpert == 1 && para8_1.npertd >= 1) {
	dprob = *dprob1;
    } else if (para8_1.idpert == 2 && para8_1.npertd >= 1) {
	dprob = 1.f / (real) para8_1.npertd;
    }
    if (*ianti == 0) {
	if (para8_1.idpert == 0 || para8_1.idpert == 1 && *ipert1 == 0 || 
		para8_1.idpert == 2 && *idloop == para8_1.npertd + 1) {
	    s_wsle(&io___2433);
	    do_lio(&c__3, &c__1, (char *)&leadng_1.lb1, (ftnlen)sizeof(
		    integer));
	    do_lio(&c__9, &c__1, " *", (ftnlen)2);
	    do_lio(&c__3, &c__1, (char *)&dpi_1.lb2, (ftnlen)sizeof(integer));
	    do_lio(&c__9, &c__1, " ->d+", (ftnlen)5);
	    do_lio(&c__3, &c__1, (char *)&(*lbm), (ftnlen)sizeof(integer));
	    do_lio(&c__9, &c__1, " (regular d prodn)    @evt#", (ftnlen)27);
	    do_lio(&c__3, &c__1, (char *)&arevt_1.iaevt, (ftnlen)sizeof(
		    integer));
	    do_lio(&c__9, &c__1, " @nt=", (ftnlen)5);
	    do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(integer));
	    e_wsle();
	} else if ((para8_1.idpert == 1 || para8_1.idpert == 2) && *idloop == 
		para8_1.npertd) {
	    s_wsle(&io___2434);
	    do_lio(&c__3, &c__1, (char *)&leadng_1.lb1, (ftnlen)sizeof(
		    integer));
	    do_lio(&c__9, &c__1, " *", (ftnlen)2);
	    do_lio(&c__3, &c__1, (char *)&dpi_1.lb2, (ftnlen)sizeof(integer));
	    do_lio(&c__9, &c__1, " ->d+", (ftnlen)5);
	    do_lio(&c__3, &c__1, (char *)&(*lbm), (ftnlen)sizeof(integer));
	    do_lio(&c__9, &c__1, " (pert d prodn)       @evt#", (ftnlen)27);
	    do_lio(&c__3, &c__1, (char *)&arevt_1.iaevt, (ftnlen)sizeof(
		    integer));
	    do_lio(&c__9, &c__1, " @nt=", (ftnlen)5);
	    do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(integer));
	    do_lio(&c__9, &c__1, " @prob=", (ftnlen)7);
	    do_lio(&c__4, &c__1, (char *)&dprob, (ftnlen)sizeof(real));
	    e_wsle();
	}
    } else {
	if (para8_1.idpert == 0 || para8_1.idpert == 1 && *ipert1 == 0 || 
		para8_1.idpert == 2 && *idloop == para8_1.npertd + 1) {
	    s_wsle(&io___2435);
	    do_lio(&c__3, &c__1, (char *)&leadng_1.lb1, (ftnlen)sizeof(
		    integer));
	    do_lio(&c__9, &c__1, " *", (ftnlen)2);
	    do_lio(&c__3, &c__1, (char *)&dpi_1.lb2, (ftnlen)sizeof(integer));
	    do_lio(&c__9, &c__1, " ->d+", (ftnlen)5);
	    do_lio(&c__3, &c__1, (char *)&(*lbm), (ftnlen)sizeof(integer));
	    do_lio(&c__9, &c__1, " (regular dbar prodn) @evt#", (ftnlen)27);
	    do_lio(&c__3, &c__1, (char *)&arevt_1.iaevt, (ftnlen)sizeof(
		    integer));
	    do_lio(&c__9, &c__1, " @nt=", (ftnlen)5);
	    do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(integer));
	    e_wsle();
	} else if ((para8_1.idpert == 1 || para8_1.idpert == 2) && *idloop == 
		para8_1.npertd) {
	    s_wsle(&io___2436);
	    do_lio(&c__3, &c__1, (char *)&leadng_1.lb1, (ftnlen)sizeof(
		    integer));
	    do_lio(&c__9, &c__1, " *", (ftnlen)2);
	    do_lio(&c__3, &c__1, (char *)&dpi_1.lb2, (ftnlen)sizeof(integer));
	    do_lio(&c__9, &c__1, " ->d+", (ftnlen)5);
	    do_lio(&c__3, &c__1, (char *)&(*lbm), (ftnlen)sizeof(integer));
	    do_lio(&c__9, &c__1, " (pert dbar prodn)    @evt#", (ftnlen)27);
	    do_lio(&c__3, &c__1, (char *)&arevt_1.iaevt, (ftnlen)sizeof(
		    integer));
	    do_lio(&c__9, &c__1, " @nt=", (ftnlen)5);
	    do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(integer));
	    do_lio(&c__9, &c__1, " @prob=", (ftnlen)7);
	    do_lio(&c__4, &c__1, (char *)&dprob, (ftnlen)sizeof(real));
	    e_wsle();
	}
    }

    return 0;
} /* bbdangle_ */


/*     Deuteron+Meson->B+B cross section (in mb) */
/* Subroutine */ int sdmbb_(real *srt, real *sdm, integer *ianti)
{
    /* System generated locals */
    real r__1, r__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real pinitial, s, threshold, xnnfactor, fs, snew;
    extern doublereal fdpiel_(real *);
    static real pfinal;
    extern doublereal fnndpi_(real *), ranart_(integer *);


    *sdm = 0.f;
    dpisig_1.sdmel = 0.f;
    dpisig_1.sdmnn = 0.f;
    dpisig_1.sdmnd = 0.f;
    dpisig_1.sdmns = 0.f;
    dpisig_1.sdmnp = 0.f;
    dpisig_1.sdmdd = 0.f;
    dpisig_1.sdmds = 0.f;
    dpisig_1.sdmdp = 0.f;
    dpisig_1.sdmss = 0.f;
    dpisig_1.sdmsp = 0.f;
    dpisig_1.sdmpp = 0.f;
/* test off check Xsec using fixed mass for resonances: */
/*      if(lb1.ge.25.and.lb1.le.27) then */
/*         em1=0.776 */
/*      elseif(lb1.eq.28) then */
/*         em1=0.783 */
/*      elseif(lb1.eq.0) then */
/*         em1=0.548 */
/*      endif */
/*      if(lb2.ge.25.and.lb2.le.27) then */
/*         em2=0.776 */
/*      elseif(lb2.eq.28) then */
/*         em2=0.783 */
/*      elseif(lb2.eq.0) then */
/*         em2=0.548 */
/*      endif */

    if (*srt <= leadng_1.em1 + dpi_1.em2) {
	return 0;
    }
/* Computing 2nd power */
    r__1 = *srt;
    s = r__1 * r__1;
/* Computing 2nd power */
    r__1 = leadng_1.em1 + dpi_1.em2;
/* Computing 2nd power */
    r__2 = leadng_1.em1 - dpi_1.em2;
    pinitial = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
    fs = fnndpi_(&s);
/*     Determine isospin and spin factors for the ratio between */
/*     Deuteron+Meson->BB and BB->Deuteron+Meson cross sections: */
    if (para8_1.idxsec == 1 || para8_1.idxsec == 2) {
/*     Assume B+B -> d+Meson has the same cross sections as N+N -> d+pi, */
/*     then determine d+Meson -> B+B cross sections: */
	if (leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5 || dpi_1.lb2 >= 3 && 
		dpi_1.lb2 <= 5) {
	    xnnfactor = .88888888888888884f;
	} else if (leadng_1.lb1 >= 25 && leadng_1.lb1 <= 27 || dpi_1.lb2 >= 
		25 && dpi_1.lb2 <= 27) {
	    xnnfactor = .29629629629629628f;
	} else if (leadng_1.lb1 == 28 || dpi_1.lb2 == 28) {
	    xnnfactor = .88888888888888884f;
	} else if (leadng_1.lb1 == 0 || dpi_1.lb2 == 0) {
	    xnnfactor = 2.6666666666666665f;
	}
    } else {
/*     Assume d+Meson -> B+B has the same cross sections as d+pi -> N+N: */
    }
/* lin-9/2008 For elastic collisions: */
    if (para8_1.idxsec == 1 || para8_1.idxsec == 3) {
/*     1/3: assume the same |matrix element|**2/s (after averaging over initial */
/*     spins and isospins) for d+Meson elastic at the same sqrt(s); */
	dpisig_1.sdmel = fdpiel_(&s);
    } else if (para8_1.idxsec == 2 || para8_1.idxsec == 4) {
/*     2/4: assume the same |matrix element|**2/s (after averaging over initial */
/*     spins and isospins) for d+Meson elastic at the same sqrt(s)-threshold: */
	threshold = leadng_1.em1 + dpi_1.em2;
/* Computing 2nd power */
	r__1 = *srt - threshold + 2.012f;
	snew = r__1 * r__1;
	dpisig_1.sdmel = fdpiel_(&snew);
    }

/*     NN: DETERMINE THE CHARGE STATES OF PARTICLESIN THE FINAL STATE */
    if ((leadng_1.lb1 == 5 || dpi_1.lb2 == 5 || leadng_1.lb1 == 27 || 
	    dpi_1.lb2 == 27) && *ianti == 0 || (leadng_1.lb1 == 3 || 
	    dpi_1.lb2 == 3 || leadng_1.lb1 == 25 || dpi_1.lb2 == 25) && *
	    ianti == 1) {
/*     (1) FOR Deuteron+(pi+,rho+) -> P+P or DeuteronBar+(pi-,rho-)-> PBar+PBar: */
	dpifsl_1.lbnn1 = 1;
	dpifsl_1.lbnn2 = 1;
	dpifsm_1.xmnn1 = .93828f;
	dpifsm_1.xmnn2 = .93828f;
    } else if (leadng_1.lb1 == 3 || dpi_1.lb2 == 3 || leadng_1.lb1 == 26 || 
	    dpi_1.lb2 == 26 || leadng_1.lb1 == 28 || dpi_1.lb2 == 28 || 
	    leadng_1.lb1 == 0 || dpi_1.lb2 == 0) {
/*     (2) FOR Deuteron+(pi0,rho0,omega,eta) -> N+P */
/*     or DeuteronBar+(pi0,rho0,omega,eta) ->NBar+PBar: */
	dpifsl_1.lbnn1 = 2;
	dpifsl_1.lbnn2 = 1;
	dpifsm_1.xmnn1 = .939457f;
	dpifsm_1.xmnn2 = .93828f;
    } else {
/*     (3) FOR Deuteron+(pi-,rho-) -> N+N or DeuteronBar+(pi+,rho+)-> NBar+NBar: */
	dpifsl_1.lbnn1 = 2;
	dpifsl_1.lbnn2 = 2;
	dpifsm_1.xmnn1 = .939457f;
	dpifsm_1.xmnn2 = .939457f;
    }
    if (*srt > dpifsm_1.xmnn1 + dpifsm_1.xmnn2) {
/* Computing 2nd power */
	r__1 = dpifsm_1.xmnn1 + dpifsm_1.xmnn2;
/* Computing 2nd power */
	r__2 = dpifsm_1.xmnn1 - dpifsm_1.xmnn2;
	pfinal = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
	if (para8_1.idxsec == 1) {
/*     1: assume the same |matrix element|**2/s (after averaging over initial */
/*     spins and isospins) for B+B -> deuteron+meson at the same sqrt(s); */
	    dpisig_1.sdmnn = fs * pfinal / pinitial * 3.f / 16.f * xnnfactor;
	} else if (para8_1.idxsec == 2 || para8_1.idxsec == 4) {
/* Computing MAX */
	    r__1 = dpifsm_1.xmnn1 + dpifsm_1.xmnn2, r__2 = leadng_1.em1 + 
		    dpi_1.em2;
	    threshold = dmax(r__1,r__2);
/* Computing 2nd power */
	    r__1 = *srt - threshold + 2.012f;
	    snew = r__1 * r__1;
	    if (para8_1.idxsec == 2) {
/*     2: assume the same |matrix element|**2/s for B+B -> deuteron+meson */
/*     at the same sqrt(s)-threshold: */
		dpisig_1.sdmnn = fnndpi_(&snew) * pfinal / pinitial * 3.f / 
			16.f * xnnfactor;
	    } else if (para8_1.idxsec == 4) {
/*     4: assume the same |matrix element|**2/s for B+B <- deuteron+meson */
/*     at the same sqrt(s)-threshold: */
		dpisig_1.sdmnn = fnndpi_(&snew) * pfinal / pinitial / 6.f;
	    }
	} else if (para8_1.idxsec == 3) {
/*     3: assume the same |matrix element|**2/s for B+B <- deuteron+meson */
/*     at the same sqrt(s): */
	    dpisig_1.sdmnn = fs * pfinal / pinitial / 6.f;
	}
    }

/*     ND: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE */
    dpifsl_1.lbnd1 = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
    dpifsl_1.lbnd2 = (integer) (ranart_(&rndf77_1.nseed) * 4) + 6;
    if (dpifsl_1.lbnd1 == 1) {
	dpifsm_1.xmnd1 = .93828f;
    } else if (dpifsl_1.lbnd1 == 2) {
	dpifsm_1.xmnd1 = .939457f;
    }
    dpifsm_1.xmnd2 = 1.232f;
    if (*srt > dpifsm_1.xmnd1 + dpifsm_1.xmnd2) {
/* Computing 2nd power */
	r__1 = dpifsm_1.xmnd1 + dpifsm_1.xmnd2;
/* Computing 2nd power */
	r__2 = dpifsm_1.xmnd1 - dpifsm_1.xmnd2;
	pfinal = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
	if (para8_1.idxsec == 1) {
/*     The spin- and isospin-averaged factor is 8-times larger for ND: */
	    dpisig_1.sdmnd = fs * pfinal / pinitial * 3.f / 16.f * (xnnfactor 
		    * 8.f);
	} else if (para8_1.idxsec == 2 || para8_1.idxsec == 4) {
/* Computing MAX */
	    r__1 = dpifsm_1.xmnd1 + dpifsm_1.xmnd2, r__2 = leadng_1.em1 + 
		    dpi_1.em2;
	    threshold = dmax(r__1,r__2);
/* Computing 2nd power */
	    r__1 = *srt - threshold + 2.012f;
	    snew = r__1 * r__1;
	    if (para8_1.idxsec == 2) {
		dpisig_1.sdmnd = fnndpi_(&snew) * pfinal / pinitial * 3.f / 
			16.f * (xnnfactor * 8.f);
	    } else if (para8_1.idxsec == 4) {
		dpisig_1.sdmnd = fnndpi_(&snew) * pfinal / pinitial / 6.f;
	    }
	} else if (para8_1.idxsec == 3) {
	    dpisig_1.sdmnd = fs * pfinal / pinitial / 6.f;
	}
    }

/*     NS: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE */
    dpifsl_1.lbns1 = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
    dpifsl_1.lbns2 = (integer) (ranart_(&rndf77_1.nseed) * 2) + 10;
    if (dpifsl_1.lbns1 == 1) {
	dpifsm_1.xmns1 = .93828f;
    } else if (dpifsl_1.lbns1 == 2) {
	dpifsm_1.xmns1 = .939457f;
    }
    dpifsm_1.xmns2 = 1.44f;
    if (*srt > dpifsm_1.xmns1 + dpifsm_1.xmns2) {
/* Computing 2nd power */
	r__1 = dpifsm_1.xmns1 + dpifsm_1.xmns2;
/* Computing 2nd power */
	r__2 = dpifsm_1.xmns1 - dpifsm_1.xmns2;
	pfinal = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
	if (para8_1.idxsec == 1) {
	    dpisig_1.sdmns = fs * pfinal / pinitial * 3.f / 16.f * (xnnfactor 
		    * 2.f);
	} else if (para8_1.idxsec == 2 || para8_1.idxsec == 4) {
/* Computing MAX */
	    r__1 = dpifsm_1.xmns1 + dpifsm_1.xmns2, r__2 = leadng_1.em1 + 
		    dpi_1.em2;
	    threshold = dmax(r__1,r__2);
/* Computing 2nd power */
	    r__1 = *srt - threshold + 2.012f;
	    snew = r__1 * r__1;
	    if (para8_1.idxsec == 2) {
		dpisig_1.sdmns = fnndpi_(&snew) * pfinal / pinitial * 3.f / 
			16.f * (xnnfactor * 2.f);
	    } else if (para8_1.idxsec == 4) {
		dpisig_1.sdmns = fnndpi_(&snew) * pfinal / pinitial / 6.f;
	    }
	} else if (para8_1.idxsec == 3) {
	    dpisig_1.sdmns = fs * pfinal / pinitial / 6.f;
	}
    }

/*     NP: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE */
    dpifsl_1.lbnp1 = (integer) (ranart_(&rndf77_1.nseed) * 2) + 1;
    dpifsl_1.lbnp2 = (integer) (ranart_(&rndf77_1.nseed) * 2) + 12;
    if (dpifsl_1.lbnp1 == 1) {
	dpifsm_1.xmnp1 = .93828f;
    } else if (dpifsl_1.lbnp1 == 2) {
	dpifsm_1.xmnp1 = .939457f;
    }
    dpifsm_1.xmnp2 = 1.535f;
    if (*srt > dpifsm_1.xmnp1 + dpifsm_1.xmnp2) {
/* Computing 2nd power */
	r__1 = dpifsm_1.xmnp1 + dpifsm_1.xmnp2;
/* Computing 2nd power */
	r__2 = dpifsm_1.xmnp1 - dpifsm_1.xmnp2;
	pfinal = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
	if (para8_1.idxsec == 1) {
	    dpisig_1.sdmnp = fs * pfinal / pinitial * 3.f / 16.f * (xnnfactor 
		    * 2.f);
	} else if (para8_1.idxsec == 2 || para8_1.idxsec == 4) {
/* Computing MAX */
	    r__1 = dpifsm_1.xmnp1 + dpifsm_1.xmnp2, r__2 = leadng_1.em1 + 
		    dpi_1.em2;
	    threshold = dmax(r__1,r__2);
/* Computing 2nd power */
	    r__1 = *srt - threshold + 2.012f;
	    snew = r__1 * r__1;
	    if (para8_1.idxsec == 2) {
		dpisig_1.sdmnp = fnndpi_(&snew) * pfinal / pinitial * 3.f / 
			16.f * (xnnfactor * 2.f);
	    } else if (para8_1.idxsec == 4) {
		dpisig_1.sdmnp = fnndpi_(&snew) * pfinal / pinitial / 6.f;
	    }
	} else if (para8_1.idxsec == 3) {
	    dpisig_1.sdmnp = fs * pfinal / pinitial / 6.f;
	}
    }

/*     DD: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE */
    dpifsl_1.lbdd1 = (integer) (ranart_(&rndf77_1.nseed) * 4) + 6;
    dpifsl_1.lbdd2 = (integer) (ranart_(&rndf77_1.nseed) * 4) + 6;
    dpifsm_1.xmdd1 = 1.232f;
    dpifsm_1.xmdd2 = 1.232f;
    if (*srt > dpifsm_1.xmdd1 + dpifsm_1.xmdd2) {
/* Computing 2nd power */
	r__1 = dpifsm_1.xmdd1 + dpifsm_1.xmdd2;
/* Computing 2nd power */
	r__2 = dpifsm_1.xmdd1 - dpifsm_1.xmdd2;
	pfinal = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
	if (para8_1.idxsec == 1) {
	    dpisig_1.sdmdd = fs * pfinal / pinitial * 3.f / 16.f * (xnnfactor 
		    * 16.f);
	} else if (para8_1.idxsec == 2 || para8_1.idxsec == 4) {
/* Computing MAX */
	    r__1 = dpifsm_1.xmdd1 + dpifsm_1.xmdd2, r__2 = leadng_1.em1 + 
		    dpi_1.em2;
	    threshold = dmax(r__1,r__2);
/* Computing 2nd power */
	    r__1 = *srt - threshold + 2.012f;
	    snew = r__1 * r__1;
	    if (para8_1.idxsec == 2) {
		dpisig_1.sdmdd = fnndpi_(&snew) * pfinal / pinitial * 3.f / 
			16.f * (xnnfactor * 16.f);
	    } else if (para8_1.idxsec == 4) {
		dpisig_1.sdmdd = fnndpi_(&snew) * pfinal / pinitial / 6.f;
	    }
	} else if (para8_1.idxsec == 3) {
	    dpisig_1.sdmdd = fs * pfinal / pinitial / 6.f;
	}
    }

/*     DS: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE */
    dpifsl_1.lbds1 = (integer) (ranart_(&rndf77_1.nseed) * 4) + 6;
    dpifsl_1.lbds2 = (integer) (ranart_(&rndf77_1.nseed) * 2) + 10;
    dpifsm_1.xmds1 = 1.232f;
    dpifsm_1.xmds2 = 1.44f;
    if (*srt > dpifsm_1.xmds1 + dpifsm_1.xmds2) {
/* Computing 2nd power */
	r__1 = dpifsm_1.xmds1 + dpifsm_1.xmds2;
/* Computing 2nd power */
	r__2 = dpifsm_1.xmds1 - dpifsm_1.xmds2;
	pfinal = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
	if (para8_1.idxsec == 1) {
	    dpisig_1.sdmds = fs * pfinal / pinitial * 3.f / 16.f * (xnnfactor 
		    * 8.f);
	} else if (para8_1.idxsec == 2 || para8_1.idxsec == 4) {
/* Computing MAX */
	    r__1 = dpifsm_1.xmds1 + dpifsm_1.xmds2, r__2 = leadng_1.em1 + 
		    dpi_1.em2;
	    threshold = dmax(r__1,r__2);
/* Computing 2nd power */
	    r__1 = *srt - threshold + 2.012f;
	    snew = r__1 * r__1;
	    if (para8_1.idxsec == 2) {
		dpisig_1.sdmds = fnndpi_(&snew) * pfinal / pinitial * 3.f / 
			16.f * (xnnfactor * 8.f);
	    } else if (para8_1.idxsec == 4) {
		dpisig_1.sdmds = fnndpi_(&snew) * pfinal / pinitial / 6.f;
	    }
	} else if (para8_1.idxsec == 3) {
	    dpisig_1.sdmds = fs * pfinal / pinitial / 6.f;
	}
    }

/*     DP: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE */
    dpifsl_1.lbdp1 = (integer) (ranart_(&rndf77_1.nseed) * 4) + 6;
    dpifsl_1.lbdp2 = (integer) (ranart_(&rndf77_1.nseed) * 2) + 12;
    dpifsm_1.xmdp1 = 1.232f;
    dpifsm_1.xmdp2 = 1.535f;
    if (*srt > dpifsm_1.xmdp1 + dpifsm_1.xmdp2) {
/* Computing 2nd power */
	r__1 = dpifsm_1.xmdp1 + dpifsm_1.xmdp2;
/* Computing 2nd power */
	r__2 = dpifsm_1.xmdp1 - dpifsm_1.xmdp2;
	pfinal = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
	if (para8_1.idxsec == 1) {
	    dpisig_1.sdmdp = fs * pfinal / pinitial * 3.f / 16.f * (xnnfactor 
		    * 8.f);
	} else if (para8_1.idxsec == 2 || para8_1.idxsec == 4) {
/* Computing MAX */
	    r__1 = dpifsm_1.xmdp1 + dpifsm_1.xmdp2, r__2 = leadng_1.em1 + 
		    dpi_1.em2;
	    threshold = dmax(r__1,r__2);
/* Computing 2nd power */
	    r__1 = *srt - threshold + 2.012f;
	    snew = r__1 * r__1;
	    if (para8_1.idxsec == 2) {
		dpisig_1.sdmdp = fnndpi_(&snew) * pfinal / pinitial * 3.f / 
			16.f * (xnnfactor * 8.f);
	    } else if (para8_1.idxsec == 4) {
		dpisig_1.sdmdp = fnndpi_(&snew) * pfinal / pinitial / 6.f;
	    }
	} else if (para8_1.idxsec == 3) {
	    dpisig_1.sdmdp = fs * pfinal / pinitial / 6.f;
	}
    }

/*     SS: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE */
    dpifsl_1.lbss1 = (integer) (ranart_(&rndf77_1.nseed) * 2) + 10;
    dpifsl_1.lbss2 = (integer) (ranart_(&rndf77_1.nseed) * 2) + 10;
    dpifsm_1.xmss1 = 1.44f;
    dpifsm_1.xmss2 = 1.44f;
    if (*srt > dpifsm_1.xmss1 + dpifsm_1.xmss2) {
/* Computing 2nd power */
	r__1 = dpifsm_1.xmss1 + dpifsm_1.xmss2;
/* Computing 2nd power */
	r__2 = dpifsm_1.xmss1 - dpifsm_1.xmss2;
	pfinal = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
	if (para8_1.idxsec == 1) {
	    dpisig_1.sdmss = fs * pfinal / pinitial * 3.f / 16.f * xnnfactor;
	} else if (para8_1.idxsec == 2 || para8_1.idxsec == 4) {
/* Computing MAX */
	    r__1 = dpifsm_1.xmss1 + dpifsm_1.xmss2, r__2 = leadng_1.em1 + 
		    dpi_1.em2;
	    threshold = dmax(r__1,r__2);
/* Computing 2nd power */
	    r__1 = *srt - threshold + 2.012f;
	    snew = r__1 * r__1;
	    if (para8_1.idxsec == 2) {
		dpisig_1.sdmss = fnndpi_(&snew) * pfinal / pinitial * 3.f / 
			16.f * xnnfactor;
	    } else if (para8_1.idxsec == 4) {
		dpisig_1.sdmss = fnndpi_(&snew) * pfinal / pinitial / 6.f;
	    }
	} else if (para8_1.idxsec == 3) {
	    dpisig_1.sdmns = fs * pfinal / pinitial / 6.f;
	}
    }

/*     SP: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE */
    dpifsl_1.lbsp1 = (integer) (ranart_(&rndf77_1.nseed) * 2) + 10;
    dpifsl_1.lbsp2 = (integer) (ranart_(&rndf77_1.nseed) * 2) + 12;
    dpifsm_1.xmsp1 = 1.44f;
    dpifsm_1.xmsp2 = 1.535f;
    if (*srt > dpifsm_1.xmsp1 + dpifsm_1.xmsp2) {
/* Computing 2nd power */
	r__1 = dpifsm_1.xmsp1 + dpifsm_1.xmsp2;
/* Computing 2nd power */
	r__2 = dpifsm_1.xmsp1 - dpifsm_1.xmsp2;
	pfinal = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
	if (para8_1.idxsec == 1) {
	    dpisig_1.sdmsp = fs * pfinal / pinitial * 3.f / 16.f * (xnnfactor 
		    * 2.f);
	} else if (para8_1.idxsec == 2 || para8_1.idxsec == 4) {
/* Computing MAX */
	    r__1 = dpifsm_1.xmsp1 + dpifsm_1.xmsp2, r__2 = leadng_1.em1 + 
		    dpi_1.em2;
	    threshold = dmax(r__1,r__2);
/* Computing 2nd power */
	    r__1 = *srt - threshold + 2.012f;
	    snew = r__1 * r__1;
	    if (para8_1.idxsec == 2) {
		dpisig_1.sdmsp = fnndpi_(&snew) * pfinal / pinitial * 3.f / 
			16.f * (xnnfactor * 2.f);
	    } else if (para8_1.idxsec == 4) {
		dpisig_1.sdmsp = fnndpi_(&snew) * pfinal / pinitial / 6.f;
	    }
	} else if (para8_1.idxsec == 3) {
	    dpisig_1.sdmsp = fs * pfinal / pinitial / 6.f;
	}
    }

/*     PP: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE */
    dpifsl_1.lbpp1 = (integer) (ranart_(&rndf77_1.nseed) * 2) + 12;
    dpifsl_1.lbpp2 = (integer) (ranart_(&rndf77_1.nseed) * 2) + 12;
    dpifsm_1.xmpp1 = 1.535f;
    dpifsm_1.xmpp2 = 1.535f;
    if (*srt > dpifsm_1.xmpp1 + dpifsm_1.xmpp2) {
/* Computing 2nd power */
	r__1 = dpifsm_1.xmpp1 + dpifsm_1.xmpp2;
/* Computing 2nd power */
	r__2 = dpifsm_1.xmpp1 - dpifsm_1.xmpp2;
	pfinal = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
	if (para8_1.idxsec == 1) {
	    dpisig_1.sdmpp = fs * pfinal / pinitial * 3.f / 16.f * xnnfactor;
	} else if (para8_1.idxsec == 2 || para8_1.idxsec == 4) {
/* Computing MAX */
	    r__1 = dpifsm_1.xmpp1 + dpifsm_1.xmpp2, r__2 = leadng_1.em1 + 
		    dpi_1.em2;
	    threshold = dmax(r__1,r__2);
/* Computing 2nd power */
	    r__1 = *srt - threshold + 2.012f;
	    snew = r__1 * r__1;
	    if (para8_1.idxsec == 2) {
		dpisig_1.sdmpp = fnndpi_(&snew) * pfinal / pinitial * 3.f / 
			16.f * xnnfactor;
	    } else if (para8_1.idxsec == 4) {
		dpisig_1.sdmpp = fnndpi_(&snew) * pfinal / pinitial / 6.f;
	    }
	} else if (para8_1.idxsec == 3) {
	    dpisig_1.sdmpp = fs * pfinal / pinitial / 6.f;
	}
    }

    *sdm = dpisig_1.sdmel + dpisig_1.sdmnn + dpisig_1.sdmnd + dpisig_1.sdmns 
	    + dpisig_1.sdmnp + dpisig_1.sdmdd + dpisig_1.sdmds + 
	    dpisig_1.sdmdp + dpisig_1.sdmss + dpisig_1.sdmsp + dpisig_1.sdmpp;
    if (*ianti == 1) {
	dpifsl_1.lbnn1 = -dpifsl_1.lbnn1;
	dpifsl_1.lbnn2 = -dpifsl_1.lbnn2;
	dpifsl_1.lbnd1 = -dpifsl_1.lbnd1;
	dpifsl_1.lbnd2 = -dpifsl_1.lbnd2;
	dpifsl_1.lbns1 = -dpifsl_1.lbns1;
	dpifsl_1.lbns2 = -dpifsl_1.lbns2;
	dpifsl_1.lbnp1 = -dpifsl_1.lbnp1;
	dpifsl_1.lbnp2 = -dpifsl_1.lbnp2;
	dpifsl_1.lbdd1 = -dpifsl_1.lbdd1;
	dpifsl_1.lbdd2 = -dpifsl_1.lbdd2;
	dpifsl_1.lbds1 = -dpifsl_1.lbds1;
	dpifsl_1.lbds2 = -dpifsl_1.lbds2;
	dpifsl_1.lbdp1 = -dpifsl_1.lbdp1;
	dpifsl_1.lbdp2 = -dpifsl_1.lbdp2;
	dpifsl_1.lbss1 = -dpifsl_1.lbss1;
	dpifsl_1.lbss2 = -dpifsl_1.lbss2;
	dpifsl_1.lbsp1 = -dpifsl_1.lbsp1;
	dpifsl_1.lbsp2 = -dpifsl_1.lbsp2;
	dpifsl_1.lbpp1 = -dpifsl_1.lbpp1;
	dpifsl_1.lbpp2 = -dpifsl_1.lbpp2;
    }
/* test off */
/*      write(98,100) srt,sdmnn,sdmnd,sdmns,sdmnp,sdmdd,sdmds,sdmdp, */
/*     1     sdmss,sdmsp,sdmpp,sdm */
/* 100  format(f5.2,11(1x,f5.1)) */

    return 0;
} /* sdmbb_ */


/* lin-9/2008 Deuteron+Meson ->B+B and elastic collisions */
/* Subroutine */ int crdmbb_(real *px, real *py, real *pz, real *srt, integer 
	*i1, integer *i2, integer *iblock, integer *ntag, real *sig, integer *
	nt, integer *ianti)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double sqrt(doublereal);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static real s;
    extern /* Subroutine */ int dmelangle_(real *, real *, real *, real *);
    static real x1;
    static integer idm, lbm;
    static real pxn, pyn, pzn;
    static integer lbb1, lbb2;
    static real e1cm, e2cm, xmb1, pt1d, pt2d, pt3d, xmb2, edcm, pt1i1, pt2i1, 
	    pt3i1, pt1i2, pt2i2, pt3i2;
    static integer ideut;
    static real p1beta, p2beta, pdbeta, scheck, pfinal;
    extern doublereal ranart_(integer *);
    static real transf;
    extern /* Subroutine */ int rotate_(real *, real *, real *, real *, real *
	    , real *), dmangle_(real *, real *, real *, integer *, integer *, 
	    real *, integer *);

    /* Fortran I/O blocks */
    static cilist io___2449 = { 0, 91, 0, 0, 0 };
    static cilist io___2450 = { 0, 91, 0, 0, 0 };
    static cilist io___2452 = { 0, 99, 0, 0, 0 };
    static cilist io___2463 = { 0, 91, 0, 0, 0 };
    static cilist io___2464 = { 0, 91, 0, 0, 0 };
    static cilist io___2469 = { 0, 99, 0, 0, 0 };
    static cilist io___2470 = { 0, 91, 0, 0, 0 };
    static cilist io___2471 = { 0, 91, 0, 0, 0 };
    static cilist io___2472 = { 0, 6, 0, 0, 0 };


/* ----------------------------------------------------------------------- */
    *iblock = 0;
    *ntag = 0;
    leadng_1.em1 = cc_1.e[*i1 - 1];
    dpi_1.em2 = cc_1.e[*i2 - 1];
/* Computing 2nd power */
    r__1 = *srt;
    s = r__1 * r__1;
    if (*sig <= 0.f) {
	return 0;
    }

    if (abs(leadng_1.lb1) == 42) {
	ideut = *i1;
	lbm = dpi_1.lb2;
	idm = *i2;
    } else {
	ideut = *i2;
	lbm = leadng_1.lb1;
	idm = *i1;
    }
/* ccc  Elastic collision or destruction of perturbatively-produced deuterons: */
    if ((para8_1.idpert == 1 || para8_1.idpert == 2) && dpert_1.dpertp[ideut 
	    - 1] != 1.f) {
/*     choose reaction channels: */
	x1 = ranart_(&rndf77_1.nseed);
	if (x1 <= dpisig_1.sdmel / *sig) {
/*     Elastic collisions: */
	    if (*ianti == 0) {
		s_wsle(&io___2449);
		do_lio(&c__9, &c__1, "  d+", (ftnlen)4);
		do_lio(&c__3, &c__1, (char *)&lbm, (ftnlen)sizeof(integer));
		do_lio(&c__9, &c__1, " (pert d M elastic) @nt=", (ftnlen)24);
		do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(integer));
		do_lio(&c__9, &c__1, " @prob=", (ftnlen)7);
		do_lio(&c__4, &c__1, (char *)&dpert_1.dpertp[ideut - 1], (
			ftnlen)sizeof(real));
		e_wsle();
	    } else {
		s_wsle(&io___2450);
		do_lio(&c__9, &c__1, "  d+", (ftnlen)4);
		do_lio(&c__3, &c__1, (char *)&lbm, (ftnlen)sizeof(integer));
		do_lio(&c__9, &c__1, " (pert dbar M elastic) @nt=", (ftnlen)
			27);
		do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(integer));
		do_lio(&c__9, &c__1, " @prob=", (ftnlen)7);
		do_lio(&c__4, &c__1, (char *)&dpert_1.dpertp[ideut - 1], (
			ftnlen)sizeof(real));
		e_wsle();
	    }
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
	    r__1 = leadng_1.em1 + dpi_1.em2;
/* Computing 2nd power */
	    r__2 = leadng_1.em1 - dpi_1.em2;
	    scheck = (s - r__1 * r__1) * (s - r__2 * r__2);
	    if (scheck < 0.f) {
		s_wsle(&io___2452);
		do_lio(&c__9, &c__1, "scheck51: ", (ftnlen)10);
		do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
		e_wsle();
		scheck = 0.f;
	    }
	    pfinal = sqrt(scheck) / 2.f / *srt;
/*            pfinal=sqrt((s-(em1+em2)**2)*(s-(em1-em2)**2))/2./srt */
	    dmelangle_(&pxn, &pyn, &pzn, &pfinal);
	    rotate_(px, py, pz, &pxn, &pyn, &pzn);
/* Computing 2nd power */
	    r__1 = cc_1.e[ideut - 1];
/* Computing 2nd power */
	    r__2 = pxn;
/* Computing 2nd power */
	    r__3 = pyn;
/* Computing 2nd power */
	    r__4 = pzn;
	    edcm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4)
		    ;
	    pdbeta = pxn * bg_1.betax + pyn * bg_1.betay + pzn * bg_1.betaz;
	    transf = bg_1.gamma * (bg_1.gamma * pdbeta / (bg_1.gamma + 1.f) + 
		    edcm);
	    pt1d = bg_1.betax * transf + pxn;
	    pt2d = bg_1.betay * transf + pyn;
	    pt3d = bg_1.betaz * transf + pzn;
	    bb_1.p[ideut * 3 - 3] = pt1d;
	    bb_1.p[ideut * 3 - 2] = pt2d;
	    bb_1.p[ideut * 3 - 1] = pt3d;
	    *iblock = 504;
	    leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
	    leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
	    leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
	    ee_1.id[*i1 - 1] = 2;
	    ee_1.id[*i2 - 1] = 2;
/*     Change the position of the perturbative deuteron to that of */
/*     the meson to avoid consecutive collisions between them: */
	    aa_1.r__[ideut * 3 - 3] = aa_1.r__[idm * 3 - 3];
	    aa_1.r__[ideut * 3 - 2] = aa_1.r__[idm * 3 - 2];
	    aa_1.r__[ideut * 3 - 1] = aa_1.r__[idm * 3 - 1];
	} else {
/*     Destruction of deuterons: */
	    if (*ianti == 0) {
		s_wsle(&io___2463);
		do_lio(&c__9, &c__1, "  d+", (ftnlen)4);
		do_lio(&c__3, &c__1, (char *)&lbm, (ftnlen)sizeof(integer));
		do_lio(&c__9, &c__1, " ->BB (pert d destrn) @nt=", (ftnlen)26)
			;
		do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(integer));
		do_lio(&c__9, &c__1, " @prob=", (ftnlen)7);
		do_lio(&c__4, &c__1, (char *)&dpert_1.dpertp[ideut - 1], (
			ftnlen)sizeof(real));
		e_wsle();
	    } else {
		s_wsle(&io___2464);
		do_lio(&c__9, &c__1, "  d+", (ftnlen)4);
		do_lio(&c__3, &c__1, (char *)&lbm, (ftnlen)sizeof(integer));
		do_lio(&c__9, &c__1, " ->BB (pert dbar destrn) @nt=", (ftnlen)
			29);
		do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(integer));
		do_lio(&c__9, &c__1, " @prob=", (ftnlen)7);
		do_lio(&c__4, &c__1, (char *)&dpert_1.dpertp[ideut - 1], (
			ftnlen)sizeof(real));
		e_wsle();
	    }
	    cc_1.e[ideut - 1] = 0.f;
	    *iblock = 502;
	}
	return 0;
    }

/* ccc  Destruction of regularly-produced deuterons: */
    *iblock = 502;
/*     choose final state and assign masses here: */
    x1 = ranart_(&rndf77_1.nseed);
    if (x1 <= dpisig_1.sdmnn / *sig) {
	lbb1 = dpifsl_1.lbnn1;
	lbb2 = dpifsl_1.lbnn2;
	xmb1 = dpifsm_1.xmnn1;
	xmb2 = dpifsm_1.xmnn2;
    } else if (x1 <= (dpisig_1.sdmnn + dpisig_1.sdmnd) / *sig) {
	lbb1 = dpifsl_1.lbnd1;
	lbb2 = dpifsl_1.lbnd2;
	xmb1 = dpifsm_1.xmnd1;
	xmb2 = dpifsm_1.xmnd2;
    } else if (x1 <= (dpisig_1.sdmnn + dpisig_1.sdmnd + dpisig_1.sdmns) / *
	    sig) {
	lbb1 = dpifsl_1.lbns1;
	lbb2 = dpifsl_1.lbns2;
	xmb1 = dpifsm_1.xmns1;
	xmb2 = dpifsm_1.xmns2;
    } else if (x1 <= (dpisig_1.sdmnn + dpisig_1.sdmnd + dpisig_1.sdmns + 
	    dpisig_1.sdmnp) / *sig) {
	lbb1 = dpifsl_1.lbnp1;
	lbb2 = dpifsl_1.lbnp2;
	xmb1 = dpifsm_1.xmnp1;
	xmb2 = dpifsm_1.xmnp2;
    } else if (x1 <= (dpisig_1.sdmnn + dpisig_1.sdmnd + dpisig_1.sdmns + 
	    dpisig_1.sdmnp + dpisig_1.sdmdd) / *sig) {
	lbb1 = dpifsl_1.lbdd1;
	lbb2 = dpifsl_1.lbdd2;
	xmb1 = dpifsm_1.xmdd1;
	xmb2 = dpifsm_1.xmdd2;
    } else if (x1 <= (dpisig_1.sdmnn + dpisig_1.sdmnd + dpisig_1.sdmns + 
	    dpisig_1.sdmnp + dpisig_1.sdmdd + dpisig_1.sdmds) / *sig) {
	lbb1 = dpifsl_1.lbds1;
	lbb2 = dpifsl_1.lbds2;
	xmb1 = dpifsm_1.xmds1;
	xmb2 = dpifsm_1.xmds2;
    } else if (x1 <= (dpisig_1.sdmnn + dpisig_1.sdmnd + dpisig_1.sdmns + 
	    dpisig_1.sdmnp + dpisig_1.sdmdd + dpisig_1.sdmds + dpisig_1.sdmdp)
	     / *sig) {
	lbb1 = dpifsl_1.lbdp1;
	lbb2 = dpifsl_1.lbdp2;
	xmb1 = dpifsm_1.xmdp1;
	xmb2 = dpifsm_1.xmdp2;
    } else if (x1 <= (dpisig_1.sdmnn + dpisig_1.sdmnd + dpisig_1.sdmns + 
	    dpisig_1.sdmnp + dpisig_1.sdmdd + dpisig_1.sdmds + dpisig_1.sdmdp 
	    + dpisig_1.sdmss) / *sig) {
	lbb1 = dpifsl_1.lbss1;
	lbb2 = dpifsl_1.lbss2;
	xmb1 = dpifsm_1.xmss1;
	xmb2 = dpifsm_1.xmss2;
    } else if (x1 <= (dpisig_1.sdmnn + dpisig_1.sdmnd + dpisig_1.sdmns + 
	    dpisig_1.sdmnp + dpisig_1.sdmdd + dpisig_1.sdmds + dpisig_1.sdmdp 
	    + dpisig_1.sdmss + dpisig_1.sdmsp) / *sig) {
	lbb1 = dpifsl_1.lbsp1;
	lbb2 = dpifsl_1.lbsp2;
	xmb1 = dpifsm_1.xmsp1;
	xmb2 = dpifsm_1.xmsp2;
    } else if (x1 <= (dpisig_1.sdmnn + dpisig_1.sdmnd + dpisig_1.sdmns + 
	    dpisig_1.sdmnp + dpisig_1.sdmdd + dpisig_1.sdmds + dpisig_1.sdmdp 
	    + dpisig_1.sdmss + dpisig_1.sdmsp + dpisig_1.sdmpp) / *sig) {
	lbb1 = dpifsl_1.lbpp1;
	lbb2 = dpifsl_1.lbpp2;
	xmb1 = dpifsm_1.xmpp1;
	xmb2 = dpifsm_1.xmpp2;
    } else {
/*     Elastic collision: */
	lbb1 = leadng_1.lb1;
	lbb2 = dpi_1.lb2;
	xmb1 = leadng_1.em1;
	xmb2 = dpi_1.em2;
	*iblock = 504;
    }
    ee_1.lb[*i1 - 1] = lbb1;
    cc_1.e[*i1 - 1] = xmb1;
    ee_1.lb[*i2 - 1] = lbb2;
    cc_1.e[*i2 - 1] = xmb2;
    leadng_1.lb1 = ee_1.lb[*i1 - 1];
    dpi_1.lb2 = ee_1.lb[*i2 - 1];
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__1 = xmb1 + xmb2;
/* Computing 2nd power */
    r__2 = xmb1 - xmb2;
    scheck = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (scheck < 0.f) {
	s_wsle(&io___2469);
	do_lio(&c__9, &c__1, "scheck52: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	scheck = 0.f;
    }
    pfinal = sqrt(scheck) / 2.f / *srt;
/*      pfinal=sqrt((s-(xmb1+xmb2)**2)*(s-(xmb1-xmb2)**2))/2./srt */
    if (*iblock == 502) {
	dmangle_(&pxn, &pyn, &pzn, nt, ianti, &pfinal, &lbm);
    } else if (*iblock == 504) {
	if (*ianti == 0) {
	    s_wsle(&io___2470);
	    do_lio(&c__9, &c__1, " d+", (ftnlen)3);
	    do_lio(&c__3, &c__1, (char *)&lbm, (ftnlen)sizeof(integer));
	    do_lio(&c__9, &c__1, " (regular d M elastic) @evt#", (ftnlen)28);
	    do_lio(&c__3, &c__1, (char *)&arevt_1.iaevt, (ftnlen)sizeof(
		    integer));
	    do_lio(&c__9, &c__1, " @nt=", (ftnlen)5);
	    do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(integer));
	    do_lio(&c__9, &c__1, " lb1,2=", (ftnlen)7);
	    do_lio(&c__3, &c__1, (char *)&leadng_1.lb1, (ftnlen)sizeof(
		    integer));
	    do_lio(&c__3, &c__1, (char *)&dpi_1.lb2, (ftnlen)sizeof(integer));
	    e_wsle();
	} else {
	    s_wsle(&io___2471);
	    do_lio(&c__9, &c__1, " d+", (ftnlen)3);
	    do_lio(&c__3, &c__1, (char *)&lbm, (ftnlen)sizeof(integer));
	    do_lio(&c__9, &c__1, " (regular dbar M elastic) @evt#", (ftnlen)
		    31);
	    do_lio(&c__3, &c__1, (char *)&arevt_1.iaevt, (ftnlen)sizeof(
		    integer));
	    do_lio(&c__9, &c__1, " @nt=", (ftnlen)5);
	    do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(integer));
	    do_lio(&c__9, &c__1, " lb1,2=", (ftnlen)7);
	    do_lio(&c__3, &c__1, (char *)&leadng_1.lb1, (ftnlen)sizeof(
		    integer));
	    do_lio(&c__3, &c__1, (char *)&dpi_1.lb2, (ftnlen)sizeof(integer));
	    e_wsle();
	}
	dmelangle_(&pxn, &pyn, &pzn, &pfinal);
    } else {
	s_wsle(&io___2472);
	do_lio(&c__9, &c__1, "Wrong iblock number in crdmbb()", (ftnlen)31);
	e_wsle();
	s_stop("", (ftnlen)0);
    }
/*     ROTATE THE MOMENTA OF PARTICLES IN THE CMS OF P1+P2 */
/*     (This is not needed for isotropic distributions) */
    rotate_(px, py, pz, &pxn, &pyn, &pzn);
/*     LORENTZ-TRANSFORMATION OF THE MOMENTUM OF PARTICLES IN THE FINAL STATE */
/*     FROM THE NUCLEUS-NUCLEUS CMS. FRAME INTO LAB FRAME: */
/*     For the 1st baryon: */
/* Computing 2nd power */
    r__1 = cc_1.e[*i1 - 1];
/* Computing 2nd power */
    r__2 = pxn;
/* Computing 2nd power */
    r__3 = pyn;
/* Computing 2nd power */
    r__4 = pzn;
    e1cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    p1beta = pxn * bg_1.betax + pyn * bg_1.betay + pzn * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1.f) + e1cm);
    pt1i1 = bg_1.betax * transf + pxn;
    pt2i1 = bg_1.betay * transf + pyn;
    pt3i1 = bg_1.betaz * transf + pzn;

    bb_1.p[*i1 * 3 - 3] = pt1i1;
    bb_1.p[*i1 * 3 - 2] = pt2i1;
    bb_1.p[*i1 * 3 - 1] = pt3i1;
/*     For the 2nd baryon: */
/* Computing 2nd power */
    r__1 = cc_1.e[*i2 - 1];
/* Computing 2nd power */
    r__2 = pxn;
/* Computing 2nd power */
    r__3 = pyn;
/* Computing 2nd power */
    r__4 = pzn;
    e2cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    p2beta = -pxn * bg_1.betax - pyn * bg_1.betay - pzn * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * p2beta / (bg_1.gamma + 1.f) + e2cm);
    pt1i2 = bg_1.betax * transf - pxn;
    pt2i2 = bg_1.betay * transf - pyn;
    pt3i2 = bg_1.betaz * transf - pzn;

    bb_1.p[*i2 * 3 - 3] = pt1i2;
    bb_1.p[*i2 * 3 - 2] = pt2i2;
    bb_1.p[*i2 * 3 - 1] = pt3i2;

    leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
    leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
    leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
    leadng_1.em1 = cc_1.e[*i1 - 1];
    dpi_1.em2 = cc_1.e[*i2 - 1];
    ee_1.id[*i1 - 1] = 2;
    ee_1.id[*i2 - 1] = 2;
    return 0;
} /* crdmbb_ */


/*     Generate angular distribution of BB from d+meson in the CMS frame: */
/* Subroutine */ int dmangle_(real *pxn, real *pyn, real *pzn, integer *nt, 
	integer *ianti, real *pfinal, integer *lbm)
{
    /* System generated locals */
    real r__1;

    /* Builtin functions */
    double sqrt(doublereal), cos(doublereal), sin(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static real c1, s1, t1, ct1, st1;
    extern doublereal ranart_(integer *);

    /* Fortran I/O blocks */
    static cilist io___2488 = { 0, 91, 0, 0, 0 };
    static cilist io___2489 = { 0, 91, 0, 0, 0 };


/*     take isotropic distribution for now: */
    c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
    t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
/* Computing 2nd power */
    r__1 = c1;
    s1 = sqrt(1.f - r__1 * r__1);
    ct1 = cos(t1);
    st1 = sin(t1);
/* THE MOMENTUM IN THE CMS IN THE FINAL STATE */
    *pzn = *pfinal * c1;
    *pxn = *pfinal * s1 * ct1;
    *pyn = *pfinal * s1 * st1;
/* lin-5/2008 track the number of regularly-destructed deuterons: */
    if (*ianti == 0) {
	s_wsle(&io___2488);
	do_lio(&c__9, &c__1, " d+", (ftnlen)3);
	do_lio(&c__3, &c__1, (char *)&(*lbm), (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, " ->BB (regular d destrn) @evt#", (ftnlen)30);
	do_lio(&c__3, &c__1, (char *)&arevt_1.iaevt, (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, " @nt=", (ftnlen)5);
	do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, " lb1,2=", (ftnlen)7);
	do_lio(&c__3, &c__1, (char *)&leadng_1.lb1, (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&dpi_1.lb2, (ftnlen)sizeof(integer));
	e_wsle();
    } else {
	s_wsle(&io___2489);
	do_lio(&c__9, &c__1, " d+", (ftnlen)3);
	do_lio(&c__3, &c__1, (char *)&(*lbm), (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, " ->BB (regular dbar destrn) @evt#", (ftnlen)33);
	do_lio(&c__3, &c__1, (char *)&arevt_1.iaevt, (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, " @nt=", (ftnlen)5);
	do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, " lb1,2=", (ftnlen)7);
	do_lio(&c__3, &c__1, (char *)&leadng_1.lb1, (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&dpi_1.lb2, (ftnlen)sizeof(integer));
	e_wsle();
    }

    return 0;
} /* dmangle_ */


/*     Angular distribution of d+meson elastic collisions in the CMS frame: */
/* Subroutine */ int dmelangle_(real *pxn, real *pyn, real *pzn, real *pfinal)
{
    /* System generated locals */
    real r__1;

    /* Builtin functions */
    double sqrt(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static real c1, s1, t1, ct1, st1;
    extern doublereal ranart_(integer *);

/*     take isotropic distribution for now: */
    c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
    t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
/* Computing 2nd power */
    r__1 = c1;
    s1 = sqrt(1.f - r__1 * r__1);
    ct1 = cos(t1);
    st1 = sin(t1);
/* THE MOMENTUM IN THE CMS IN THE FINAL STATE */
    *pzn = *pfinal * c1;
    *pxn = *pfinal * s1 * ct1;
    *pyn = *pfinal * s1 * st1;
    return 0;
} /* dmelangle_ */


/* lin-9/2008 Deuteron+Baryon elastic cross section (in mb) */
/* Subroutine */ int sdbelastic_(real *srt, real *sdb)
{
    /* System generated locals */
    real r__1;

    /* Local variables */
    static real s, threshold, snew;
    extern doublereal fdbel_(real *);
    static real sdbel;


    *sdb = 0.f;
    sdbel = 0.f;
    if (*srt <= leadng_1.em1 + dpi_1.em2) {
	return 0;
    }
/* Computing 2nd power */
    r__1 = *srt;
    s = r__1 * r__1;
/*     For elastic collisions: */
    if (para8_1.idxsec == 1 || para8_1.idxsec == 3) {
/*     1/3: assume the same |matrix element|**2/s (after averaging over initial */
/*     spins and isospins) for d+Baryon elastic at the same sqrt(s); */
	sdbel = fdbel_(&s);
    } else if (para8_1.idxsec == 2 || para8_1.idxsec == 4) {
/*     2/4: assume the same |matrix element|**2/s (after averaging over initial */
/*     spins and isospins) for d+Baryon elastic at the same sqrt(s)-threshold: */
	threshold = leadng_1.em1 + dpi_1.em2;
/* Computing 2nd power */
	r__1 = *srt - threshold + 2.012f;
	snew = r__1 * r__1;
	sdbel = fdbel_(&snew);
    }
    *sdb = sdbel;
    return 0;
} /* sdbelastic_ */

/* lin-9/2008 Deuteron+Baryon elastic collisions */
/* Subroutine */ int crdbel_(real *px, real *py, real *pz, real *srt, integer 
	*i1, integer *i2, integer *iblock, integer *ntag, real *sig, integer *
	nt, integer *ianti)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double sqrt(doublereal);

    /* Local variables */
    static real s;
    extern /* Subroutine */ int dbelangle_(real *, real *, real *, real *);
    static integer idb, lbb;
    static real pxn, pyn, pzn, e1cm, e2cm, pt1d, pt2d, pt3d, edcm, pt1i1, 
	    pt2i1, pt3i1, pt1i2, pt2i2, pt3i2;
    static integer ideut;
    static real p1beta, p2beta, pdbeta, scheck, pfinal, transf;
    extern /* Subroutine */ int rotate_(real *, real *, real *, real *, real *
	    , real *);

    /* Fortran I/O blocks */
    static cilist io___2503 = { 0, 91, 0, 0, 0 };
    static cilist io___2504 = { 0, 91, 0, 0, 0 };
    static cilist io___2506 = { 0, 99, 0, 0, 0 };
    static cilist io___2517 = { 0, 91, 0, 0, 0 };
    static cilist io___2518 = { 0, 91, 0, 0, 0 };
    static cilist io___2519 = { 0, 99, 0, 0, 0 };


/* ----------------------------------------------------------------------- */
    *iblock = 0;
    *ntag = 0;
    leadng_1.em1 = cc_1.e[*i1 - 1];
    dpi_1.em2 = cc_1.e[*i2 - 1];
/* Computing 2nd power */
    r__1 = *srt;
    s = r__1 * r__1;
    if (*sig <= 0.f) {
	return 0;
    }
    *iblock = 503;

    if (abs(leadng_1.lb1) == 42) {
	ideut = *i1;
	lbb = dpi_1.lb2;
	idb = *i2;
    } else {
	ideut = *i2;
	lbb = leadng_1.lb1;
	idb = *i1;
    }
/* ccc  Elastic collision of perturbatively-produced deuterons: */
    if ((para8_1.idpert == 1 || para8_1.idpert == 2) && dpert_1.dpertp[ideut 
	    - 1] != 1.f) {
	if (*ianti == 0) {
	    s_wsle(&io___2503);
	    do_lio(&c__9, &c__1, "  d+", (ftnlen)4);
	    do_lio(&c__3, &c__1, (char *)&lbb, (ftnlen)sizeof(integer));
	    do_lio(&c__9, &c__1, " (pert d B elastic) @nt=", (ftnlen)24);
	    do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(integer));
	    do_lio(&c__9, &c__1, " @prob=", (ftnlen)7);
	    do_lio(&c__4, &c__1, (char *)&dpert_1.dpertp[ideut - 1], (ftnlen)
		    sizeof(real));
	    do_lio(&c__4, &c__1, (char *)&bb_1.p[idb * 3 - 3], (ftnlen)sizeof(
		    real));
	    do_lio(&c__4, &c__1, (char *)&bb_1.p[idb * 3 - 2], (ftnlen)sizeof(
		    real));
	    do_lio(&c__4, &c__1, (char *)&bb_1.p[ideut * 3 - 3], (ftnlen)
		    sizeof(real));
	    do_lio(&c__4, &c__1, (char *)&bb_1.p[ideut * 3 - 2], (ftnlen)
		    sizeof(real));
	    e_wsle();
	} else {
	    s_wsle(&io___2504);
	    do_lio(&c__9, &c__1, "  d+", (ftnlen)4);
	    do_lio(&c__3, &c__1, (char *)&lbb, (ftnlen)sizeof(integer));
	    do_lio(&c__9, &c__1, " (pert dbar Bbar elastic) @nt=", (ftnlen)30)
		    ;
	    do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(integer));
	    do_lio(&c__9, &c__1, " @prob=", (ftnlen)7);
	    do_lio(&c__4, &c__1, (char *)&dpert_1.dpertp[ideut - 1], (ftnlen)
		    sizeof(real));
	    do_lio(&c__4, &c__1, (char *)&bb_1.p[idb * 3 - 3], (ftnlen)sizeof(
		    real));
	    do_lio(&c__4, &c__1, (char *)&bb_1.p[idb * 3 - 2], (ftnlen)sizeof(
		    real));
	    do_lio(&c__4, &c__1, (char *)&bb_1.p[ideut * 3 - 3], (ftnlen)
		    sizeof(real));
	    do_lio(&c__4, &c__1, (char *)&bb_1.p[ideut * 3 - 2], (ftnlen)
		    sizeof(real));
	    e_wsle();
	}
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
	r__1 = leadng_1.em1 + dpi_1.em2;
/* Computing 2nd power */
	r__2 = leadng_1.em1 - dpi_1.em2;
	scheck = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (scheck < 0.f) {
	    s_wsle(&io___2506);
	    do_lio(&c__9, &c__1, "scheck53: ", (ftnlen)10);
	    do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	    e_wsle();
	    scheck = 0.f;
	}
	pfinal = sqrt(scheck) / 2.f / *srt;
/*         pfinal=sqrt((s-(em1+em2)**2)*(s-(em1-em2)**2))/2./srt */
	dbelangle_(&pxn, &pyn, &pzn, &pfinal);
	rotate_(px, py, pz, &pxn, &pyn, &pzn);
/* Computing 2nd power */
	r__1 = cc_1.e[ideut - 1];
/* Computing 2nd power */
	r__2 = pxn;
/* Computing 2nd power */
	r__3 = pyn;
/* Computing 2nd power */
	r__4 = pzn;
	edcm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	pdbeta = pxn * bg_1.betax + pyn * bg_1.betay + pzn * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * pdbeta / (bg_1.gamma + 1.f) + 
		edcm);
	pt1d = bg_1.betax * transf + pxn;
	pt2d = bg_1.betay * transf + pyn;
	pt3d = bg_1.betaz * transf + pzn;
	bb_1.p[ideut * 3 - 3] = pt1d;
	bb_1.p[ideut * 3 - 2] = pt2d;
	bb_1.p[ideut * 3 - 1] = pt3d;
	leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
	leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
	leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
	ee_1.id[*i1 - 1] = 2;
	ee_1.id[*i2 - 1] = 2;
/*     Change the position of the perturbative deuteron to that of */
/*     the baryon to avoid consecutive collisions between them: */
	aa_1.r__[ideut * 3 - 3] = aa_1.r__[idb * 3 - 3];
	aa_1.r__[ideut * 3 - 2] = aa_1.r__[idb * 3 - 2];
	aa_1.r__[ideut * 3 - 1] = aa_1.r__[idb * 3 - 1];
	return 0;
    }

/*     Elastic collision of regularly-produced deuterons: */
    if (*ianti == 0) {
	s_wsle(&io___2517);
	do_lio(&c__9, &c__1, " d+", (ftnlen)3);
	do_lio(&c__3, &c__1, (char *)&lbb, (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, " (regular d B elastic) @evt#", (ftnlen)28);
	do_lio(&c__3, &c__1, (char *)&arevt_1.iaevt, (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, " @nt=", (ftnlen)5);
	do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, " lb1,2=", (ftnlen)7);
	do_lio(&c__3, &c__1, (char *)&leadng_1.lb1, (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&dpi_1.lb2, (ftnlen)sizeof(integer));
	e_wsle();
    } else {
	s_wsle(&io___2518);
	do_lio(&c__9, &c__1, " d+", (ftnlen)3);
	do_lio(&c__3, &c__1, (char *)&lbb, (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, " (regular dbar Bbar elastic) @evt#", (ftnlen)34)
		;
	do_lio(&c__3, &c__1, (char *)&arevt_1.iaevt, (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, " @nt=", (ftnlen)5);
	do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, " lb1,2=", (ftnlen)7);
	do_lio(&c__3, &c__1, (char *)&leadng_1.lb1, (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&dpi_1.lb2, (ftnlen)sizeof(integer));
	e_wsle();
    }
/* lin-9/2012: check argument in sqrt(): */
/* Computing 2nd power */
    r__1 = leadng_1.em1 + dpi_1.em2;
/* Computing 2nd power */
    r__2 = leadng_1.em1 - dpi_1.em2;
    scheck = (s - r__1 * r__1) * (s - r__2 * r__2);
    if (scheck < 0.f) {
	s_wsle(&io___2519);
	do_lio(&c__9, &c__1, "scheck54: ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&scheck, (ftnlen)sizeof(real));
	e_wsle();
	scheck = 0.f;
    }
    pfinal = sqrt(scheck) / 2.f / *srt;
/*      pfinal=sqrt((s-(em1+em2)**2)*(s-(em1-em2)**2))/2./srt */
    dbelangle_(&pxn, &pyn, &pzn, &pfinal);
/*     ROTATE THE MOMENTA OF PARTICLES IN THE CMS OF P1+P2 */
/*     (This is not needed for isotropic distributions) */
    rotate_(px, py, pz, &pxn, &pyn, &pzn);
/*     LORENTZ-TRANSFORMATION OF THE MOMENTUM OF PARTICLES IN THE FINAL STATE */
/*     FROM THE NUCLEUS-NUCLEUS CMS. FRAME INTO LAB FRAME: */
/*     For the 1st baryon: */
/* Computing 2nd power */
    r__1 = cc_1.e[*i1 - 1];
/* Computing 2nd power */
    r__2 = pxn;
/* Computing 2nd power */
    r__3 = pyn;
/* Computing 2nd power */
    r__4 = pzn;
    e1cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    p1beta = pxn * bg_1.betax + pyn * bg_1.betay + pzn * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1.f) + e1cm);
    pt1i1 = bg_1.betax * transf + pxn;
    pt2i1 = bg_1.betay * transf + pyn;
    pt3i1 = bg_1.betaz * transf + pzn;

    bb_1.p[*i1 * 3 - 3] = pt1i1;
    bb_1.p[*i1 * 3 - 2] = pt2i1;
    bb_1.p[*i1 * 3 - 1] = pt3i1;
/*     For the 2nd baryon: */
/* Computing 2nd power */
    r__1 = cc_1.e[*i2 - 1];
/* Computing 2nd power */
    r__2 = pxn;
/* Computing 2nd power */
    r__3 = pyn;
/* Computing 2nd power */
    r__4 = pzn;
    e2cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
    p2beta = -pxn * bg_1.betax - pyn * bg_1.betay - pzn * bg_1.betaz;
    transf = bg_1.gamma * (bg_1.gamma * p2beta / (bg_1.gamma + 1.f) + e2cm);
    pt1i2 = bg_1.betax * transf - pxn;
    pt2i2 = bg_1.betay * transf - pyn;
    pt3i2 = bg_1.betaz * transf - pzn;

    bb_1.p[*i2 * 3 - 3] = pt1i2;
    bb_1.p[*i2 * 3 - 2] = pt2i2;
    bb_1.p[*i2 * 3 - 1] = pt3i2;

    leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
    leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
    leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
    leadng_1.em1 = cc_1.e[*i1 - 1];
    dpi_1.em2 = cc_1.e[*i2 - 1];
    ee_1.id[*i1 - 1] = 2;
    ee_1.id[*i2 - 1] = 2;
    return 0;
} /* crdbel_ */


/*     Part of the cross section function of NN->Deuteron+Pi (in mb): */
doublereal fnndpi_(real *s)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3;

    /* Builtin functions */
    double exp(doublereal);

    if (*s <= 4.0481439999999997f) {
	ret_val = 0.f;
    } else {
/* Computing 2nd power */
	r__1 = *s - 4.65f;
/* Computing 2nd power */
	r__2 = *s - 4.65f;
/* Computing 2nd power */
	r__3 = *s - 6.f;
	ret_val = exp(-(r__1 * r__1) / .1f) * 26.f + exp(-(r__2 * r__2) / 2.f)
		 * 4.f + exp(-(r__3 * r__3) / 10.f) * .28f;
    }
    return ret_val;
} /* fnndpi_ */


/*     Angular distribution of d+baryon elastic collisions in the CMS frame: */
/* Subroutine */ int dbelangle_(real *pxn, real *pyn, real *pzn, real *pfinal)
{
    /* System generated locals */
    real r__1;

    /* Builtin functions */
    double sqrt(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static real c1, s1, t1, ct1, st1;
    extern doublereal ranart_(integer *);

/*     take isotropic distribution for now: */
    c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
    t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
/* Computing 2nd power */
    r__1 = c1;
    s1 = sqrt(1.f - r__1 * r__1);
    ct1 = cos(t1);
    st1 = sin(t1);
/* THE MOMENTUM IN THE CMS IN THE FINAL STATE */
    *pzn = *pfinal * c1;
    *pxn = *pfinal * s1 * ct1;
    *pyn = *pfinal * s1 * st1;
    return 0;
} /* dbelangle_ */


/*     Cross section of Deuteron+Pi elastic (in mb): */
doublereal fdpiel_(real *s)
{
    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double exp(doublereal);

    if (*s <= 4.0481439999999997f) {
	ret_val = 0.f;
    } else {
/* Computing 2nd power */
	r__1 = *s - 4.67f;
/* Computing 2nd power */
	r__2 = *s - 6.25f;
	ret_val = exp(-(r__1 * r__1) / .15f) * 63.f + exp(-(r__2 * r__2) / 
		.3f) * 15.f;
    }
    return ret_val;
} /* fdpiel_ */


/*     Cross section of Deuteron+N elastic (in mb): */
doublereal fdbel_(real *s)
{
    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double exp(doublereal);

    if (*s <= 4.0481439999999997f) {
	ret_val = 0.f;
    } else {
/* Computing 2nd power */
	r__1 = *s - 7.93f;
/* Computing 2nd power */
	r__2 = *s - 7.93f;
	ret_val = exp(-(r__1 * r__1) / .003f) * 2500.f + exp(-(r__2 * r__2) / 
		.1f) * 300.f + 10.f;
    }
    return ret_val;
} /* fdbel_ */

