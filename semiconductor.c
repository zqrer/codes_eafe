#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "phg.h"

/* NOTE:
   1) The folloing codes implement the FVSG method,
      Zlamal FEM, Inverse Average Mixed FEM, and Exoponential 
      Basis FEM.
   2) The variables that we use are carrier concentrations.
   3) Only Gummel iteration is implemented.
   4) SRH Recombination model is implemented.
   5) The devices we simulate include PN diode, NPN bipolar junction
      transistor, and nMOSFET. The NPN transistor case is not carefully
      test, and may have problems.
*/

static const char *poisson_names[] = { "fvm", "fem", NULL };
static const char *continuity_names[] =
    { "box", "mixed", "zlamal", "exp", "eafe", NULL };
static int poisson_index = 0;
static int continuity_index = 0;

static enum { PN = 0, NPN = 1, PNP = 2, MOS = 3 } device = PN;

typedef enum { POISSON = 0, ELECTRON = 1, HOLE = 2 } EQUATION;

typedef enum { REGION_N = 1, REGION_P = 2,	/* PN DIODE */
    REGION_E = 1, REGION_B = 2, REGION_C = 3,	/* PNP/NPN BJT */
    REGION_S = 1, REGION_D = 2, CHANNEL = 3, REGION_ST = 4, OXIDE = 5,	/* MOSFET */
} REGION;

typedef enum { CATHODE = BDRY_USER1, ANODE = BDRY_USER2,
    JUNCTION_PN = BDRY_USER3,	/* PN DIODE */
    EMITTER = BDRY_USER1, COLLECTOR = BDRY_USER2, BASE = BDRY_USER3,
    JUNCTION_EB = BDRY_USER4, JUNCTION_CB = BDRY_USER5,	/* BJT */
    SOURCE = BDRY_USER1, DRAIN = BDRY_USER2, SUBSTRATE = BDRY_USER4,
    JUNCTION_SS = BDRY_USER5, JUNCTION_DS = BDRY_USER6, JUNCTION_OS = BDRY_USER7,	/* MOSFET */
    GATE = BDRY_USER8
} CONTACT;

typedef enum { MAT_SEMICONDUCTOR = 1, MAT_OXIDE = 2 } MATERIAL;

/*---- constants from Sentaurus-Device manual ----*/
static const FLOAT KB = 8.617343e-5 /* eV/K */ ;	/* Boltzmann constant */
static const FLOAT Q = 1.602176487e-19 /* C */ ;	/* charge of a proton */
static const FLOAT EPSILONE0 = 8.8541878176e-12 /* F/m */ ;	/* vacumm permittivity */
static const FLOAT CM2MICRON = 1.0e4;
static const FLOAT T = 300.0 /* K */ ;	/* absolute temperture */

/* silicon relative permittivity */
static const FLOAT EPSILONE = 11.9;
/* silicon band structure */
static const FLOAT NCF = 1.5;
static const FLOAT NVF = 1.5;
static const FLOAT NC300 = 2.89e+19 /* 1/cm^3 */ ;
static const FLOAT NV300 = 3.14e+19 /* 1/cm^3 */ ;
static const FLOAT EG0 = 1.16964 /* eV */ ;
static const FLOAT EGALPHA = 4.73e-4 /* eV/K */ ;
static const FLOAT EGBETA = 636.0 /* K */ ;
/* oxide */
static const FLOAT EPSILONE_O = 3.9;	/* oxide relative permittivity */
/* metal */
//static const FLOAT WORKFUNC = 4.17 /* eV */ ;	/* metal work function */
static const FLOAT WORKFUNC = 0.0 /* eV */ ;	/* metal work function */
/* Shockley-Read-Hall(SRH) Recombination*/
static const FLOAT TAUN0 = 1.0e-7 /* s */ ;	/* electron life time */
static const FLOAT TAUP0 = 1.0e-7 /* s */ ;	/* hole life time */
/* constant mobility model */
static const FLOAT MUN = 1417.0 /* cm^2/(V*s) */ ;
static const FLOAT MUP = 470.5 /* cm^2/(V*s) */ ;
static const FLOAT TMUN = 2.5;
static const FLOAT TMUP = 2.2;
static FLOAT nie /* 1/cm^3 */ ;	/* intrinsic carrier concentration */
static FLOAT Dn /* cm^2/s */ ;
static FLOAT Dp /* cm^2/s */ ;

/*---- device specifications ----*/
static BTYPE ohmic_contact;
/** PN diode **/
static FLOAT doping_pregion = 1.0e19 /* 1/cm^3 */ ;	/* P region doping */
static FLOAT doping_nregion = 1.0e16 /* 1/cm^3 */ ;	/* N region doping */
static FLOAT bias_anode = 0.05, bias_cathode = 0. /* V */ ;
static FLOAT anode = 0., cathode = 0. /* V */ ;	/* runtime */
static FLOAT h = 0.05 /* V */ ;	/* bias increment */
static FLOAT resistor_anode = 0. /* Omega */ ;	/* contact resistor */
static FLOAT resistor_cathode = 0. /* Omega */ ;	/* contact resistor */
/** NPN/PNP **/
static FLOAT doping_e = 1.0e18 /* 1/cm^3 */ ;	/* BJT emitter region doping */
static FLOAT doping_c = 1.0e14 /* 1/cm^3 */ ;	/* BJT collector region doping */
static FLOAT doping_b = 1.0e16 /* 1/cm^3 */ ;	/* BJT base region doping */
static FLOAT bias_base = 0., bias_emitter = 0., bias_collector = 0. /* V */ ;
static FLOAT base = 0., emitter = 0., collector = 0. /* V */ ;	/* runtime */
static FLOAT h_b = 0.05 /* V */ ;	/* base bias increment */
static FLOAT h_c = 0.05 /* V */ ;	/* collector bias increment */
static FLOAT h_e = 0.05 /* V */ ;	/* emitter bias increment */
static FLOAT resistor_emitter = 0. /* Omega */ ;	/* contact resistor */
static FLOAT resistor_collector = 0. /* Omega */ ;	/* contact resistor */
static FLOAT resistor_base = 0. /* Omega */ ;	/* contact resistor */
/** MOS **/
static FLOAT t_ox = 2.5e-6 /* cm */ ;	/* oxide thickness */
static FLOAT bias_source = 0., bias_drain = 1, bias_gate =
    8., bias_substrate = 0. /* V */ ;
static FLOAT source, drain, gate, substrate /* V */ ;	/* runtime */
static FLOAT doping_s = 1.0e19;	/* SOURCE */
static FLOAT doping_d = 1.0e19;	/* DRAIN */
static FLOAT doping_st = 1.0e16;	/* SUBSTRATE */
static FLOAT h_g = 0.05 /* V */ ;	/* gate bias increment */
static FLOAT h_d = 0.05 /* V */ ;	/* drain bias increment */
//static FLOAT Not = 1.0e12;
static SHORT *verts_region = NULL;
static const FLOAT SCALE = 1.0e14;

static int recombination = 1;	/* 1 ==> with SRH recombination */
static int V_d_count = 0, V_g_count = 0;
static int version_eafe = 1;
static void
InitConstant(void)
{
    /* the intrinsic carrier concentration,
     * bandgap narrowing effect neglected */
    FLOAT NC = Pow(T / 300.0, NCF) * NC300;
    FLOAT NV = Pow(T / 300.0, NVF) * NV300;
    FLOAT EG = EG0 - EGALPHA * T * T / (T + EGBETA);
    nie = Sqrt(NC * NV) * Exp(-EG / (2.0 * KB * T)) /* 1/cm^3 */ ;

    /* constant mobility model from sentaurus.pdf p296 */
    FLOAT mu_n = MUN * Pow(T / 300.0, -TMUN);
    FLOAT mu_p = MUP * Pow(T / 300.0, -TMUP);
    Dn = KB * T * mu_n /* cm^2/s */ ;
    Dp = KB * T * mu_p /* cm^2/s */ ;
    //phgPrintf("Dn = %le Dp = %le\n",Dn, Dp);
}

static void
SRH(FLOAT v, FLOAT p, FLOAT n, FLOAT nie, EQUATION eqn, FLOAT *lv, FLOAT *rv)
/* recombination model */
{
    FLOAT c;

    c = TAUP0 * (n + nie/SCALE) + TAUN0 * (p + nie/SCALE);
    if (rv != NULL)
	*rv = (nie*nie/SCALE/SCALE - p * n) / c;
    if (lv != NULL)
	*lv = eqn == ELECTRON ? p / c : n / c;
}

static void
Doping(int region_mark, FLOAT *ND, FLOAT *NA)
{
    if (device == PN) {
	switch (region_mark) {
	    case REGION_P:
		*ND = 0.;
		*NA = doping_pregion;
		return;
	    case REGION_N:
		*ND = doping_nregion;
		*NA = 0.;
		return;
	}
    }
    else if (device == NPN) {
	switch (region_mark) {
	    case REGION_E:
		*ND = doping_e;
		*NA = 0.;
		return;
	    case REGION_C:
		*ND = doping_c;
		*NA = 0.;
		return;
	    case REGION_B:
		*ND = 0.;
		*NA = doping_b;
		return;
	}
    }
    else if (device == PNP) {
	switch (region_mark) {
	    case REGION_E:
		*ND = 0.;
		*NA = doping_e;
		return;
	    case REGION_C:
		*ND = 0.;
		*NA = doping_c;
		return;
	    case REGION_B:
		*ND = doping_b;
		*NA = 0.;
		return;
	}
    }
    else if (device == MOS) {
	switch (region_mark) {
	    case REGION_S:
		*ND = doping_s;
		*NA = 0.;
		return;
	    case REGION_D:
		*ND = doping_d;
		*NA = 0.;
		return;
	    case REGION_ST:
		*ND = 0.;
		*NA = doping_st;
		return;
	}
    }
}

static void
PrintDopingProfile(void)
{
    if (device == PN) {
	phgPrintf("\n  Doping profile: P region %0.4le, N region %0.4le\n",
		  (double)doping_pregion, (double)doping_nregion);
    }
    else if (device == NPN) {
	phgPrintf
	    ("\n  Doping profile: emitter %0.4le, base %0.4le, collector %0.4le\n",
	     (double)doping_e, (double)doping_b, (double)doping_c);
    }
    else if (device == PNP) {
	phgPrintf
	    ("\n  Doping profile: emitter %0.4le, base %0.4le, collector %0.4le\n",
	     (double)doping_e, (double)doping_b, (double)doping_c);
    }
    else if (device == MOS) {
	phgPrintf
	    ("\n  Doping profile: source %0.4le, drain %0.4le, substrate %0.4le\n",
	     (double)doping_s, (double)doping_d, (double)doping_st);
    }
}

static BOOLEAN
BiasSweepPN(void)
{
    /* cathode grounded, anode ramps up/down */
    static BOOLEAN initilized = FALSE;

    if (!initilized) {
	initilized = TRUE;
    }

    if (Fabs(anode) + 1.0e-6 > Fabs(bias_anode))
	return TRUE;

    if (bias_anode >= bias_cathode)
	anode += h;
    else
	anode -= h;

    phgPrintf("\n  ================ BIAS: %0.4lf ================\n",
	      (double)(anode - bias_cathode));

    return FALSE;
}

static BOOLEAN
BiasSweepNPN(void)
{
    static BOOLEAN stage1 = TRUE;
    static BOOLEAN stage2 = FALSE;
    static BOOLEAN initialized = FALSE;
    FLOAT eps = 1.0e-6;
#if 0
    /* emitter grounded, first base goes up to bias_base, then
       collector goes up to bias_collector */
    if (stage1) {
	if (!initialized) {
	    initialized = TRUE;
	}

	if (base + eps > bias_base) {
	    stage1 = FALSE;
	    stage2 = TRUE;
	    initialized = FALSE;
	}
	else
	    base += h_b;
    }

    if (stage2) {
	if (!initialized) {
	    initialized = TRUE;
	}

	if (collector + eps > bias_collector)
	    return TRUE;

	collector += h_c;
    }
#else
    /* emitter grounded, first collector goes up to bias_collector, then
       base goes up to bias_base */
    if (stage1) {
	if (!initialized) {
	    initialized = TRUE;
	}

	if (collector + eps > bias_collector) {
	    stage1 = FALSE;
	    stage2 = TRUE;
	    initialized = FALSE;
	}
	else
	    collector += h_c;
    }

    if (stage2) {
	if (!initialized) {
	    initialized = TRUE;
	}

	if (base + eps > bias_base)
	    return TRUE;

	base += h_b;
    }
#endif
#if 0
    /* base grounded, first collector goes up to bias_collector, then
       emitter goes up to bias_emitter */
    if (stage1) {
	if (!initialized) {
	    initialized = TRUE;
	}

	if (collector + eps > bias_collector) {
	    stage1 = FALSE;
	    stage2 = TRUE;
	    initialized = FALSE;
	}
	else
	    collector += h_c;
    }

    if (stage2) {
	if (!initialized) {
	    initialized = TRUE;
	}

	if (Fabs(emitter) + eps > Fabs(bias_emitter))
	    return TRUE;

	emitter -= h_e;
    }
#endif
    phgPrintf
	("\n  ============ V_b: %0.4lf, V_c: %0.4lf, V_e: %0.4lf ============\n",
	 (double)base, (double)collector, (double)emitter);

    return FALSE;
}

static BOOLEAN
BiasSweepPNP(void)
{
    static BOOLEAN stage1 = TRUE;
    static BOOLEAN stage2 = FALSE;
    static BOOLEAN initialized = FALSE;
    FLOAT eps = 1.0e-6;
#if 0
    /* emitter grounded, first base goes up to bias_base, then
       collector goes up to bias_collector */
    if (stage1) {
	if (!initialized) {
	    initialized = TRUE;
	}

	if (Fabs(base) + eps > Fabs(bias_base)) {
	    stage1 = FALSE;
	    stage2 = TRUE;
	    initialized = FALSE;
	}
	else
	    base -= h_b;
    }

    if (stage2) {
	if (!initialized) {
	    initialized = TRUE;
	}

	if (Fabs(collector) + eps > Fabs(bias_collector))
	    return TRUE;

	collector -= h_c;
    }
#else 
    /* emitter grounded, first collector goes up to bias_collector, then
       base goes up to bias_base */
    if (stage1) {
	if (!initialized) {
	    initialized = TRUE;
	}

	if (Fabs(collector) + eps > Fabs(bias_collector)) {
	    stage1 = FALSE;
	    stage2 = TRUE;
	    initialized = FALSE;
	}
	else
	    collector -= h_c;
    }

    if (stage2) {
	if (!initialized) {
	    initialized = TRUE;
	}

	if (Fabs(base) + eps > Fabs(bias_base))
	    return TRUE;

	base -= h_b;
    }
#endif
#if 0
    /* base grounded, first collector goes up to bias_collector, then
       emitter goes up to bias_emitter */
    if (stage1) {
	if (!initialized) {
	    initialized = TRUE;
	}

	if (Fabs(collector) + eps > Fabs(bias_collector)) {
	    stage1 = FALSE;
	    stage2 = TRUE;
	    initialized = FALSE;
	}
	else
	    collector -= h_c;
    }

    if (stage2) {
	if (!initialized) {
	    initialized = TRUE;
	}

	if (emitter + eps > bias_emitter)
	    return TRUE;

	emitter += h_e;
    }
#endif
    phgPrintf
	("\n  ============ V_b: %0.4lf, V_c: %0.4lf, V_e: %0.4lf ============\n",
	 (double)base, (double)collector, (double)emitter);

    return FALSE;
}

static BOOLEAN
BiasSweepMOS(void)
{
    static BOOLEAN stage1 = TRUE;
    static BOOLEAN stage2 = FALSE;
    static BOOLEAN initialized = FALSE;
    FLOAT eps = 1.0e-6;

#if 1
    /* substrate and source grounded, first drain goes up to bias_drain, then
       gate goes up to bias_gate */
    if (stage1) {
	if (!initialized) {
	    initialized = TRUE;
	}

	if (Fabs(drain) + eps > Fabs(bias_drain)) {
	    stage1 = FALSE;
	    stage2 = TRUE;
	    initialized = FALSE;
	}
	else
	    drain += h_d;
    }

    if (stage2) {
	if (!initialized) {
	    initialized = TRUE;
	}

 //  if(Fabs(gate) + eps > Fabs(bias_gate))
	if (gate + eps > bias_gate)
	    return TRUE;

	gate += h_g;
    }
#else
    /* substrate and source grounded, first gate goes up to bias_gate, then
       drain goes up to bias_drain */
    if (stage1) {
	if (!initialized) {
	    initialized = TRUE;
	}

//	if (Fabs(gate) + eps > Fabs(bias_gate)) {
	if (gate + eps > bias_gate){
	    stage1 = FALSE;
	    stage2 = TRUE;
	    initialized = FALSE;
	}
	else
	    gate += h_g;
    }

    if (stage2) {
	if (!initialized) {
	    initialized = TRUE;
	}

	if (drain + eps > bias_drain)
	    return TRUE;

	drain += h_d;
    }
#endif

    phgPrintf
	("\n  =========== V_g: %0.4lf, V_s: %0.4lf, V_d: %0.4lf, V_st: %0.4lf ===========\n",
	 (double)gate, (double)source, (double)drain, (double)substrate);

    return FALSE;
}

static BOOLEAN
BiasSweep(void)
{
    if (device == PN)
	return BiasSweepPN();
    else if (device == NPN)
	return BiasSweepNPN();
    else if (device == PNP)
	return BiasSweepPNP();
    else if (device == MOS)
	return BiasSweepMOS();

    return TRUE;
}


static void
SetupVertsRegion(GRID *g, BOOLEAN flag)
{
    SIMPLEX *e;
    int i;

    if (flag == FALSE) {
	phgFree(verts_region);
	verts_region = NULL;
	return;
    }

    assert(verts_region == NULL);
    verts_region = phgCalloc(g->nvert, sizeof(*verts_region));

    ForAllElements(g, e) {
	for (i = 0; i < NVert; i++) {
	    if (e->region_mark == OXIDE)
		*(verts_region + e->verts[i]) |= MAT_OXIDE;
	    else		/* now there are only two types of material */
		*(verts_region + e->verts[i]) |= MAT_SEMICONDUCTOR;

	    if (g->types_vert[e->verts[i]] & JUNCTION_OS) {
		/* to remedy the possibility that JUNCTION_OS is also the 
		   interface between processes */
		*(verts_region + e->verts[i]) |= MAT_OXIDE;
		*(verts_region + e->verts[i]) |= MAT_SEMICONDUCTOR;
	    }
	}
    }
}

#define B(x) (Fabs(x)>1.0e-4? (x)/(Exp(x)-1): \
	(((-1./720 *(x)*(x) + 1./12)*(x) -1./2)*(x) +1.))

static void
exponential(GRID *g, SIMPLEX *e, DOF *u, EQUATION eqn, QUAD *quad, int n,
	    const FLOAT *a, FLOAT bas[], FLOAT flux[][Dim])
/* cf: Three-dimentional exponentially fitted conforming tetrahedral 
   finite elements for the semiconductor continuity equations */
{
    int i, j;
    FLOAT x[Dim];
    FLOAT lm[NVert][Dim];
    FLOAT sm[NVert], bs[NVert], nbs[NVert];
    FLOAT bas_val[NVert];
    FLOAT A[Dim][Dim], b[Dim][NVert];
    BOOLEAN flag;
    const FLOAT *lambda = quad->points + n * (Dim + 1);
    FLOAT sum;

    phgGeomLambda2XYZ(g, e, quad->points + n * (Dim + 1), &(x[0]), &(x[1]),
		      &(x[2]));

    for (i = 0; i < NVert; i++) {
	sm[i] = 0;
	for (j = 0; j < Dim; j++) {
	    lm[i][j] = x[j] - g->verts[e->verts[i]][j];
	    if (eqn == ELECTRON)	/* electrons */
		sm[i] += a[j] * lm[i][j];
	    else if (eqn == HOLE)	/* holes */
		sm[i] += -a[j] * lm[i][j];
	}
	bs[i] = B(sm[i]);
	nbs[i] = B(-sm[i]);
    }

    sum = 0.;
    for (i = 0; i < NVert; i++)
	sum += lambda[i] * bs[i];

    for (i = 0; i < NVert; i++)
	bas_val[i] = lambda[i] * nbs[i] / sum;

    if (bas != NULL) {
	for (i = 0; i < NVert; i++)
	    bas[i] = bas_val[i];
    }

    if (flux == NULL)
	return;

    bzero(&(b[0][0]), Dim * NVert * sizeof(b[0][0]));
    for (i = 0; i < Dim; i++) {
	for (j = 0; j < Dim; j++)
	    A[i][j] = lm[i][j];
	for (j = 0; j < NVert; j++)
	    b[i][j] = bs[i] * bas_val[j];
	b[i][i] += -nbs[i];
    }

    flag = phgSolverDenseSolver(Dim, NVert, &(A[0][0]), &(b[0][0]));
    if (!flag)
	phgError(1, "%s:%d,unexpected error.\n", __FILE__, __LINE__);

    for (i = 0; i < NVert; i++)
	for (j = 0; j < Dim; j++)
	    flux[i][j] = b[j][i];
}

static int
cmp_float(const void *p0, const void *p1)
{
    FLOAT f = *((FLOAT *)p0) - *((FLOAT *)p1);
    return (f < 0.) ? -1 : ((f > 0.) ? 1 : 0);
}



static FLOAT *
JBJ(DOF *V, SIMPLEX *e, int n, EQUATION eqn)
{
    GRID *g = V->g;
    static FLOAT JT[Dim][Dim];
    FLOAT IJT[Dim][Dim];
    FLOAT psi0, psi1, b;
    int order[NVert];
    int i, j;

    for (i = 0; i < NVert; i++)
	order[i] = (n + i) % NVert;

    psi0 = *(DofVertexData(V, e->verts[order[0]]));
    for (i = 1; i < NVert; i++) {
	psi1 = *(DofVertexData(V, e->verts[order[i]]));
	b = eqn == ELECTRON ? B(psi0 - psi1) : B(psi1 - psi0);
	for (j = 0; j < Dim; j++) {
	    JT[i - 1][j] =
		g->verts[e->verts[order[i]]][j] -
		g->verts[e->verts[order[0]]][j];
	    IJT[i - 1][j] = JT[i - 1][j];
	    JT[i - 1][j] *= b;	/* B x JT */
	}
    }

    /* JT^-1 x B x JT */
    phgSolverDenseSolver(Dim, Dim, &(IJT[0][0]), &(JT[0][0]));

    return JT[0];
}

static void
exp_quad_dof_eval(SIMPLEX *e, DOF *u, EQUATION eqn, DOF *gradv,
		  QUAD *quad, FLOAT *val)
{
    int i, n;
    FLOAT bas_val[NVert];
    const FLOAT *a;

    a = DofElementData(gradv, e->index);
    for (n = 0; n < quad->npoints; n++) {
	exponential(u->g, e, u, eqn, quad, n, a, bas_val, NULL);
	val[n] = 0.;
	for (i = 0; i < NVert; i++)
	    val[n] += *(DofVertexData(u, e->verts[i])) * bas_val[i];
    }
}

static void
InitialGuess(DOF *V, DOF *P, DOF *N)
/* use values of potential and carrier concentrations at thermal
   equilibrium states as initial values */
{
    GRID *g = V->g;
    SIMPLEX *e;
    FLOAT ND, NA, phi, v, p, n, D;
    BTYPE interface, btype;
    int i;
    INT ind;
    VEC *vecV, *vecP, *vecN;
    VEC *diagV, *diagP, *diagN;
    MAP *map;

    map = phgMapCreate(V, NULL);
    vecV = phgMapCreateVec(map, 1);
    vecP = phgMapCreateVec(map, 1);
    vecN = phgMapCreateVec(map, 1);
    diagV = phgMapCreateVec(map, 1);
    diagP = phgMapCreateVec(map, 1);
    diagN = phgMapCreateVec(map, 1);
    phgVecDisassemble(vecV);
    phgVecDisassemble(vecP);
    phgVecDisassemble(vecN);
    phgVecDisassemble(diagV);
    phgVecDisassemble(diagP);
    phgVecDisassemble(diagN);

    if (device == PN) {
	assert(cathode == anode);
	phi = anode / (KB * T);
	interface = JUNCTION_PN;
    }
    else if (device == NPN) {
   #if 1
    /*emitter grounded*/
	assert(emitter == base && emitter == collector);
	phi = emitter / (KB * T);
   #else
   /*base grounded*/
	assert(base == emitter && base == collector);
	phi = base / (KB * T);
   #endif
	interface = JUNCTION_EB | JUNCTION_CB;
    }
    else if (device == PNP) {
   #if 1
    /*emitter grounded*/
	assert(emitter == base && emitter == collector);
	phi = emitter / (KB * T);
   #else
   /*base grounded*/
	assert(base == emitter && base == collector);
	phi = base / (KB * T);
   #endif
	interface = JUNCTION_EB | JUNCTION_CB;
    }
    else if (device == MOS) {
//	assert(source == substrate);
	assert(source == substrate && drain == substrate);
//	assert(gate == substrate);
	phi = substrate / (KB * T);
	interface = JUNCTION_SS | JUNCTION_DS;
    }

    ForAllElements(g, e) {
	for (i = 0; i < NVert; i++) {
	    ind = phgMapE2L(map, 0, e, i);
	    if (e->region_mark == OXIDE) {
		if (g->types_vert[e->verts[i]] & JUNCTION_OS)
		    continue;
		v = (gate - WORKFUNC) / (KB * T);
		phgVecAddEntry(vecV, 0, ind, v);
		phgVecAddEntry(diagV, 0, ind, 1.);
	    }
	    else {
		btype = phgDofGetElementBoundaryType(V, e, i);
		if (btype & interface) {
		    v = phi;
		    p = n = nie / SCALE;
		}
		else {
		    Doping(e->region_mark, &ND, &NA);
		    D = (ND - NA) / nie;
		    if (D > 0.)
			v = phi + Log(0.5 * (D + Sqrt(D * D + 4.0)));
		    else
			v = phi - Log(0.5 * (-D + Sqrt(D * D + 4.0)));

		    p = nie * Exp(phi - v) / SCALE;
		    n = nie * Exp(v - phi) / SCALE;
		}
		phgVecAddEntry(vecV, 0, ind, v);
		phgVecAddEntry(diagV, 0, ind, 1.);
		phgVecAddEntry(vecP, 0, ind, p);
		phgVecAddEntry(diagP, 0, ind, 1.);
		phgVecAddEntry(vecN, 0, ind, n);
		phgVecAddEntry(diagN, 0, ind, 1.);
	    }
	}
    }

    phgVecAssemble(vecV);
    phgVecAssemble(vecP);
    phgVecAssemble(vecN);
    phgVecAssemble(diagV);
    phgVecAssemble(diagP);
    phgVecAssemble(diagN);

    for (i = 0; i < map->nlocal; i++) {
	*(vecV->data + i) /= *(diagV->data + i);
	if (*(diagP->data + i) > 0.5)
	    *(vecP->data + i) /= *(diagP->data + i);
	if (*(diagN->data + i) > 0.5)
	    *(vecN->data + i) /= *(diagN->data + i);
    }

    phgVecDestroy(&diagV);
    phgVecDestroy(&diagP);
    phgVecDestroy(&diagN);

    phgMapVecToDofArrays(map, vecV, FALSE, &V, NULL);
    phgMapVecToDofArrays(map, vecP, FALSE, &P, NULL);
    phgMapVecToDofArrays(map, vecN, FALSE, &N, NULL);

    phgMapDestroy(&map);
    phgVecDestroy(&vecV);
    phgVecDestroy(&vecP);
    phgVecDestroy(&vecN);
}

static void
TakeBoundary(DOF *V, DOF *P, DOF *N, int who)
    /* who == 0 ==> take bdry for V,
     * who == 1 ==> take bdry for P and N,
     * who == 2 ==> take bdry for V, P, and N */
{
    GRID *g = V->g;
    SIMPLEX *e;
    BTYPE btype;
    FLOAT ND, NA;
    FLOAT psi, phi, bias;
    FLOAT D;
    int i;

    ForAllElements(g, e) {
	if (e->region_mark == OXIDE) {
	    for (i = 0; i < NVert; i++) {
		btype = phgDofGetElementBoundaryType(V, e, i);
		if (btype & GATE) {
		    psi = (gate - WORKFUNC) / (KB * T);
		    if (who == 0 || who == 2) {
			*(V->data + phgDofMapE2D(V, e, i)) = psi;
		    }
		}
	    }
	    continue;
	}

	for (i = 0; i < NVert; i++) {
	    if (phgDofGetElementBoundaryType(V, e, i) & ohmic_contact)
		break;
	}
	if (i >= NVert)		/* no bdry d.o.f */
	    continue;

	Doping(e->region_mark, &ND, &NA);

	for (i = 0; i < NVert; i++) {
	    btype = phgDofGetElementBoundaryType(V, e, i);
	    if (!(btype & ohmic_contact))
		continue;

	    if (device == PN) {
		bias = (btype & CATHODE) ? cathode : anode;
	    }
	    else if (device == NPN) {
		bias = (btype & EMITTER) ? emitter :
		    (btype & COLLECTOR ? collector : base);
	    }
	    else if (device == PNP) {
		bias = (btype & EMITTER) ? emitter :
		    (btype & COLLECTOR ? collector : base);
	    }
	    else if (device == MOS) {
		bias = (btype & SOURCE) ? source :
		    ((btype & DRAIN) ? drain : substrate);
	    }

	    D = (ND - NA) / nie;
	    phi = bias / (KB * T);

	    if (D > 0)
		psi = phi + Log(0.5 * (D + Sqrt(D * D + 4.0)));
	    else
		psi = phi - Log(0.5 * (-D + Sqrt(D * D + 4.0)));

	    if (who == 1 || who == 2) {
		*(P->data + phgDofMapE2D(P, e, i)) =
		    nie * Exp(phi - psi) / SCALE;
		*(N->data + phgDofMapE2D(N, e, i)) =
		    nie * Exp(psi - phi) / SCALE;
	    }
	    if (who == 0 || who == 2)
		*(V->data + phgDofMapE2D(V, e, i)) = psi;
	}
    }
}

FLOAT Average2(FLOAT x)
{
   if(Fabs(x)>1.0e4)
   {
      if(x < 0.0)
         return 0.0;
      else 
         return 1.0;
   }
      return Exp(x) / (1 + Exp(x));
}

FLOAT Average3_0(FLOAT x)
{
   FLOAT alpha = 0.5;
#if 1
   if(Fabs(x)<1.0e-4)
      alpha = 0.5;
   else
   {
      if(x > 0)
         alpha = 1. / x;
      else
         alpha = 1 + 1. /x;
   }
#endif
      return Exp(x) / (1 - alpha + alpha * Exp(x));
}
FLOAT Average3_1(FLOAT x)
{
   FLOAT alpha = 0.5;
#if 1
   if(Fabs(x)<1.0e-4)
      alpha = 0.5;
   else
   {
      if(x > 0)
         alpha = 1. / x;
      else
         alpha = 1 + 1. / x;
   }
#endif
      return 1.0 / (1 - alpha + alpha * Exp(x));
}

static void
average(FLOAT local_v[4], EQUATION eqn, FLOAT average_coeff[2], int psi0, int psi1)
{
    FLOAT s0, s1;
    FLOAT p[4];
    FLOAT p1, p2, p3, p4;
    FLOAT d41, d34, d31, d23, d24, d21;
    FLOAT v0, v1;
    FLOAT eps = 1.0e-3;
    FLOAT r;
    int i;

    for (i = 0; i < 4; i++) {
	if (eqn == ELECTRON)
	    p[i] = -local_v[i];
	else
	    p[i] = +local_v[i];
    }

    s0 = p[psi0];
    s1 = p[psi1];

    qsort(p, 4, sizeof(p[0]), cmp_float);
    p1 = p[0];
    p2 = p[1];
    p3 = p[2];
    p4 = p[3];

    d41 = p4 - p1;
    d34 = p3 - p4;
    d31 = p3 - p1;
    d23 = p2 - p3;
    d24 = p2 - p4;
    d21 = p2 - p1;

#define IB(x) (Fabs(x)>eps ? (Exp(x)-1.)/(x): \
	(((1./24 *(x) + 1./6)*(x) + 1./2)*(x) +1.))
#define W(x,y) (1./2 + 1./6*(x+y)+1./24*(x*x+x*y+y*y)+ \
	1./120*(x+y)*(x*x+y*y))
#define T(x,y,z) (1./6 + 1./24*(x+y+z) + 1./120*(x*x+y*y+z*z+x*y+x*z+y*z) + \
	1./720*(x*x*x+y*y*y+z*z*z+(y+z)*x*x+(y*y+z*y+z*z)*x+z*y*y+z*z*y))

    if (Fabs(d34) > eps)
	v0 = 1. / d34 * (Exp(d34) * IB(d23) - IB(d24));
    else
	v0 = IB(d23) * IB(d34) - W(d23, d24);

    if (Fabs(d31) > eps)
	v1 = 1. / d31 * (Exp(d31) * IB(d23) - IB(d21));
    else
	v1 = IB(d23) * IB(d31) - W(d23, d21);

    if (Fabs(d41) > eps) {
	average_coeff[0] = 6. / d41 * (Exp(p4 - s0) * v0 - Exp(p1 - s0) * v1);
	average_coeff[1] = 6. / d41 * (Exp(p4 - s1) * v0 - Exp(p1 - s1) * v1);
    }
    else {
	r = IB(d41) * v0;
	r -= IB(d23) * W(d34, d31);
	r += T(d24, d21, d23);
	average_coeff[0] = 6. * r * Exp(p1 - s0);
	average_coeff[1] = 6. * r * Exp(p1 - s1);
    }

#undef IB
#undef W
#undef T
}
static void
current_eafe(DOF *V, DOF *P, DOF *N, CONTACT cont[], int ncont,
	       FLOAT currents[][3])
    /* currents[][0]: total currents,
       currents[][1]: electron currents, 
       currents[][2]: hole currents */
{
    GRID *g = V->g;
    SIMPLEX *e;
    BYTE flag[ncont];
    int i, j, k, c, order = 2;
    FLOAT psi0, psi1, vp[NVert], vn[NVert], local_v[NVert];
    FLOAT coeff_p, coeff_n;
    FLOAT rp, rn;
    FLOAT Bp_ij, Bp_ii, Bn_ij, Bn_ii;
    FLOAT Ap[NVert][NVert], An[NVert][NVert];
    FLOAT Ap_bas, An_bas;
    FLOAT averagep_coeff[2], averagen_coeff[2];


    bzero(currents[0], ncont * 3 * sizeof(currents[0][0]));
    
    ForAllElements(g, e) {

	if (e->region_mark == OXIDE)
	    continue;

	bzero(flag, ncont * sizeof(flag[0]));
	for (i = 0; i < NVert; i++) {
	    for (j = 0; j < ncont; j++)
		if (g->types_vert[e->verts[i]] & cont[j])
		    flag[j] = 1;
	}

	for (i = 0; i < NVert; i++) {
	    local_v[i] = *(DofVertexData(V, e->verts[i]));
	    vp[i] = *(DofVertexData(P, e->verts[i]));
	    vn[i] = *(DofVertexData(N, e->verts[i]));
	}

	coeff_n = Dn * SCALE;
	coeff_p = Dp * SCALE;

	for (c = 0; c < ncont; c++) {
	    if (flag[c] == 0)
		continue;

	    for (i = 0; i < NVert; i++) {
		if (!(g->types_vert[e->verts[i]] & cont[c]))
		    continue;

      psi0 = local_v[i];
      Ap[i][i] = An[i][i] = 0.0;
   if(version_eafe == 1){
      for(j = 0; j < NVert; j++)
      {
         if(j != i)
         {
            psi1 = local_v[j];
            Bp_ij = B(psi0 - psi1);
            Bp_ii = B(psi1 - psi0);
            Bn_ij = B(psi1 - psi0);
            Bn_ii = B(psi0 - psi1);
            Ap_bas = phgQuadGradBasDotGradBas(e, P, j, P, i, 3);
            An_bas = phgQuadGradBasDotGradBas(e, N, j, N, i, 3);
            Ap[i][j] = coeff_p * Bp_ij * Ap_bas;
            Ap[i][i] -= coeff_p * Bp_ii * Ap_bas;
            An[i][j] = coeff_n * Bn_ij * An_bas;
            An[i][i] -= coeff_n * Bn_ii * An_bas;
         }
      }
   }
   else if(version_eafe == 2){
      for(j = 0; j < NVert; j++)
      {
         if(j != i)
         {
            psi1 = local_v[j];
            Bp_ii = Average2(psi0 - psi1);
            Bp_ij = 1.0 - Bp_ii;
            Bn_ii = Average2(psi1 - psi0);
            Bn_ij = 1.0 - Bn_ii;
            Ap_bas = phgQuadGradBasDotGradBas(e, P, j, P, i, 3);
            An_bas = phgQuadGradBasDotGradBas(e, N, j, N, i, 3);
            Ap[i][j] = coeff_p * 2 * Bp_ij * Ap_bas;
            Ap[i][i] -= coeff_p * 2 * Bp_ii * Ap_bas;
            An[i][j] = coeff_n * 2 * Bn_ij * An_bas;
            An[i][i] -= coeff_n * 2 * Bn_ii * An_bas;
         }
      }
   }
   else if(version_eafe == 3){
      for(j = 0; j < NVert; j++)
      {
         if(j != i)
         {
            psi1 = local_v[j];
            Bp_ij = Exp(psi1 - psi0);
            Bp_ii = Exp(psi0 - psi1);
            Bn_ij = Exp(psi0 - psi1);
            Bn_ii = Exp(psi1 - psi0);
            Ap_bas = phgQuadGradBasDotGradBas(e, P, j, P, i, 3);
            An_bas = phgQuadGradBasDotGradBas(e, N, j, N, i, 3);
            Ap[i][j] = coeff_p / 2 * (1 + Bp_ij) * Ap_bas;
            Ap[i][i] -= coeff_p / 2 * (1 + Bp_ii) * Ap_bas;
            An[i][j] = coeff_n / 2 * (1 + Bn_ij) * An_bas;
            An[i][i] -= coeff_n / 2 * (1 + Bn_ii) * An_bas;
         }
      }
   }
   else if(version_eafe == 4){
      for(j = 0; j < NVert; j++)
      {
         if(j != i)
         {
	         bzero(averagep_coeff, 2 * sizeof(averagep_coeff[0]));
	         bzero(averagen_coeff, 2 * sizeof(averagen_coeff[0]));
            average(local_v, HOLE, averagep_coeff, i, j);
            average(local_v, ELECTRON, averagen_coeff, i, j);
            Ap_bas = phgQuadGradBasDotGradBas(e, P, j, P, i, 3);
            An_bas = phgQuadGradBasDotGradBas(e, N, j, N, i, 3);
            Ap[i][j] = coeff_p / averagep_coeff[1] * Ap_bas;
            Ap[i][i] -= coeff_p / averagep_coeff[0] * Ap_bas;
            An[i][j] = coeff_n / averagen_coeff[1] * An_bas;
            An[i][i] -= coeff_n / averagen_coeff[0] * An_bas;
         }
      }
   }
   else if(version_eafe == 5)
   {
      for(j = 0; j < NVert; j++)
      {
         if(j != i)
         {
            psi1 = local_v[j];
            Bp_ii = Average3_0(psi0 - psi1);
            Bp_ij = Average3_1(psi0 - psi1);
            Bn_ii = Average3_0(psi1 - psi0);
            Bn_ij = Average3_1(psi1 - psi0);
            Ap_bas = phgQuadGradBasDotGradBas(e, P, j, P, i, 3);
            An_bas = phgQuadGradBasDotGradBas(e, N, j, N, i, 3);
            Ap[i][j] = coeff_p * Bp_ij * Ap_bas;
            Ap[i][i] -= coeff_p * Bp_ii * Ap_bas;
            An[i][j] = coeff_n * Bn_ij * An_bas;
            An[i][i] -= coeff_n * Bn_ii * An_bas;
         }
   }
   }
      rp = rn = 0.0;
      for(j = 0; j < NVert; j++)
      {
         rp += Ap[i][j] * vp[j];
         rn += An[i][j] * vn[j];
      }
		    currents[c][1] += rn;
		    currents[c][2] -= rp;
	    }/*loop on NVert*/
	}/*loop on ncont*/
    }/*loop on element*/

    for (c = 0; c < ncont; c++) {
	/* handle units */
	for (i = 1; i < 3; i++)
	    currents[c][i] /= CM2MICRON;
	currents[c][0] = currents[c][1] + currents[c][2];
    }

#if USE_MPI
    if (g->nprocs > 1) {
	FLOAT currents0[ncont][3];
	memcpy(currents0[0], currents[0],
	       ncont * 3 * sizeof(currents0[0][0]));
	MPI_Allreduce(currents0[0], currents[0], ncont * 3, PHG_MPI_FLOAT,
		      MPI_SUM, g->comm);
    }
#endif
}

static void
current_zlamal(DOF *V, DOF *P, DOF *N, CONTACT cont[], int ncont,
	       FLOAT currents[][3])
    /* currents[][0]: total currents,
       currents[][1]: electron currents, 
       currents[][2]: hole currents */
{
    GRID *g = V->g;
    SIMPLEX *e;
    BYTE flag[ncont];
    int i, j, k, n, c, order = 2;
    QUAD *quad = phgQuadGetQuad3D(order);
    int v0, v1;
    const FLOAT *w, *dphi[NVert];
    FLOAT psi[NVert], vp[NVert], vn[NVert];
    FLOAT expdp[NVert], expdn[NVert];
    FLOAT coeff_p, coeff_n;
    FLOAT dp[Dim], dn[Dim], fp[Dim], fn[Dim];
    FLOAT *jbj_p, *jbj_n;
    FLOAT rp, rn;

    bzero(currents[0], ncont * 3 * sizeof(currents[0][0]));
    jbj_p = (FLOAT *)malloc(Dim * Dim * sizeof(*jbj_p));
    jbj_n = (FLOAT *)malloc(Dim * Dim * sizeof(*jbj_n));

    ForAllElements(g, e) {
	FLOAT vol = phgGeomGetVolume(g, e);

	if (e->region_mark == OXIDE)
	    continue;

	bzero(flag, ncont * sizeof(flag[0]));
	for (i = 0; i < NVert; i++) {
	    for (j = 0; j < ncont; j++)
		if (g->types_vert[e->verts[i]] & cont[j])
		    flag[j] = 1;
	}

	for (i = 0; i < NVert; i++) {
	    dphi[i] = phgQuadGetBasisGradient(e, N, i, quad);
	    psi[i] = *(DofVertexData(V, e->verts[i]));
	    vp[i] = *(DofVertexData(P, e->verts[i]));
	    vn[i] = *(DofVertexData(N, e->verts[i]));
	}

	coeff_n = Dn * SCALE;
	coeff_p = Dp * SCALE;

	for (c = 0; c < ncont; c++) {
	    if (flag[c] == 0)
		continue;

	    for (i = 0; i < NVert; i++) {
		if (!(g->types_vert[e->verts[i]] & cont[c]))
		    continue;

		memcpy(jbj_p, JBJ(V, e, i, HOLE), Dim * Dim * sizeof(*jbj_p));
		memcpy(jbj_n, JBJ(V, e, i, ELECTRON),Dim * Dim * sizeof(*jbj_n));
		FLOAT (*Jp)[Dim] = (void *)jbj_p;
		FLOAT (*Jn)[Dim] = (void *)jbj_n;

		for (j = 0; j < NVert; j++) {
		    expdn[j] = Exp(psi[i] - psi[j]);
		    expdp[j] = Exp(psi[j] - psi[i]);
		}

		w = quad->weights;
		for (n = 0; n < quad->npoints; n++) {
		    dp[0] = dp[1] = dp[2] = 0.;
		    dn[0] = dn[1] = dn[2] = 0.;
		    for (j = 0; j < NVert; j++) {
			for (k = 0; k < Dim; k++) {
			    dp[k] += vp[j] * expdp[j] * dphi[j][n * Dim + k];
			    dn[k] += vn[j] * expdn[j] * dphi[j][n * Dim + k];
			}
		    }

		    rp = rn = 0.;
		    for (k = 0; k < Dim; k++) {
			fp[k] =
			    Jp[k][0] * dp[0] + Jp[k][1] * dp[1] +
			    Jp[k][2] * dp[2];
			rp += fp[k] * dphi[i][n * Dim + k];
			fn[k] =
			    Jn[k][0] * dn[0] + Jn[k][1] * dn[1] +
			    Jn[k][2] * dn[2];
			rn += fn[k] * dphi[i][n * Dim + k];
		    }
		    rp *= *w * vol;
		    rn *= *w * vol;

		    currents[c][1] += coeff_n * rn;
		    currents[c][2] -= coeff_p * rp;

		    w++;
		}
	    }
	}
    }

    for (c = 0; c < ncont; c++) {
	/* handle units */
	for (i = 1; i < 3; i++)
	    currents[c][i] /= CM2MICRON;
	currents[c][0] = currents[c][1] + currents[c][2];
    }

#if USE_MPI
    if (g->nprocs > 1) {
	FLOAT currents0[ncont][3];
	memcpy(currents0[0], currents[0],
	       ncont * 3 * sizeof(currents0[0][0]));
	MPI_Allreduce(currents0[0], currents[0], ncont * 3, PHG_MPI_FLOAT,
		      MPI_SUM, g->comm);
    }
#endif

    free(jbj_p);
    free(jbj_n);
}

static void
build_linear_system_poisson_fem(SOLVER *solver, DOF *V, DOF *P, DOF *N)
{
    int i, j, n, f, order;
    GRID *g = V->g;
    SIMPLEX *e;
    DOF *gradv = phgDofGradient(V, NULL, NULL, NULL);
    FLOAT A[NVert][NVert], rhs[NVert];
    FLOAT oxide_A[NVert][NVert], oxide_rhs[NVert];
    FLOAT vol, ck, cr, qc, coeff, expp, expn;
    FLOAT NA, ND;
    FLOAT *valp, *valn, *valu;
    BTYPE btype;
    QUAD *quad;
    const FLOAT *w, *phi[NVert];
    INT I[NVert];
    DOF *NOT;
    order = 2;			/* order 1 too less, order 3 too more, 2011.5.7 */
    quad = phgQuadGetQuad3D(order);

    valp = phgAlloc(quad->npoints * 3 * sizeof(*valp));
    valn = valp + quad->npoints;
    valu = valp + 2 * quad->npoints;
    
    NOT = phgDofNew(g,DOF_DEFAULT,1,"NOT",DofInterpolation);
    phgDofSetDataByValue(NOT,0.0);

    ForAllElements(g, e) {
	bzero(A[0], NVert * NVert * sizeof(A[0][0]));
	bzero(rhs, NVert * sizeof(rhs[0]));

	if (e->region_mark != OXIDE) {
	    Doping(e->region_mark, &ND, &NA);

	    for (i = 0; i < NVert; i++)
		phi[i] = phgQuadGetBasisValues(e, V, i, quad);

	    exp_quad_dof_eval(e, P, HOLE, gradv, quad, valp);
	    exp_quad_dof_eval(e, N, ELECTRON, gradv, quad, valn);

	    memcpy(valu, phgQuadGetDofValues(e, V, quad),
		   quad->npoints * sizeof(*valu));
	    vol = phgGeomGetVolume(g, e);

	    qc = 100.0 * Q / (EPSILONE0 * KB * T) /* cm */ ;
	    qc /= CM2MICRON * CM2MICRON * CM2MICRON;	/* handle units */

	    w = quad->weights;
	    for (n = 0; n < quad->npoints; n++) {
		ck = qc * (valp[n] + valn[n]) * SCALE;
		cr = qc * ((valp[n] - valn[n]) * SCALE + (ND - NA));

		for (i = 0; i < NVert; i++) {
		    rhs[i] += cr * phi[i][n] * (*w);
		    for (j = 0; j < NVert; j++)
			A[i][j] += ck * phi[i][n] * phi[j][n] * (*w);
		}
		w++;
	    }

	    for (i = 0; i < NVert; i++) {
		rhs[i] *= vol;
		for (j = 0; j < NVert; j++)
		    A[i][j] *= vol;
	    }
	}


	coeff = e->region_mark == OXIDE ? EPSILONE_O : EPSILONE;
	coeff /= CM2MICRON;	/* handle units */
	for (i = 0; i < NVert; i++) {
            if(e->region_mark != OXIDE){
	    rhs[i] -=
		coeff * phgQuadDofDotGradBas(e, gradv, V, i, QUAD_DEFAULT);
		}
	    else {
		rhs[i] -=
      coeff * phgQuadDofDotGradBas(e, gradv, V, i, QUAD_DEFAULT);// + qc*phgQuadDofDotBas(e,NOT,V,i,QUAD_DEFAULT);

		}
	    for (j = 0; j < NVert; j++)
		A[i][j] += coeff * phgQuadGradBasDotGradBas
		    (e, V, i, V, j, QUAD_DEFAULT);
	}

	/* handle bdry */
	for (i = 0; i < NVert; i++) {
	    btype = phgDofGetElementBoundaryType(V, e, i);
	    if (btype & (ohmic_contact | GATE)) {
		for (j = 0; j < NVert; j++)
		    A[i][j] = A[j][i] = 0.0;
		rhs[i] = 0.0;
		A[i][i] = 1.0;
	    }
       /*interface*/
    if (device == MOS) {	/* source and substrate grounded */
         for(j = 0; j < NFace; j++)
         {
            if(e->bound_type[j] & JUNCTION_OS)
            {
               rhs[i] += qc * CM2MICRON * phgQuadFaceDofDotBas(e, j, NOT, DOF_PROJ_NONE, V, i, QUAD_DEFAULT);
            }
         }
    }
	}

	/* add entries */
	for (i = 0; i < NVert; i++)
	    I[i] = phgSolverMapE2L(solver, 0, e, i);
	phgSolverAddMatrixEntries(solver, NVert, I, NVert, I, A[0]);
	phgSolverAddRHSEntries(solver, NVert, I, rhs);
    }

    phgFree(valp);
    phgDofFree(&gradv);
    phgDofFree(&NOT);
}

static void
build_linear_system_continuity_eafe(SOLVER *solver, DOF *V, DOF *P, DOF *N, EQUATION eqn)
{
    GRID *g = V->g;
    SIMPLEX *e;
    DOF *dof = eqn == ELECTRON ? N : P;
    FLOAT A[NVert][NVert], rhs[NVert], A_bas;
    FLOAT u[NVert];
    INT I[NVert];
    int i, j, k, n;
    int order = 2;
    QUAD *quad = phgQuadGetQuad3D(order);
    const FLOAT* bas_val[NVert];
    FLOAT ND, NA, vol, r, coeff, psi0, psi1;
    const FLOAT *w, *dphi[NVert];
    FLOAT B_ij, B_ii;
    FLOAT local_v[NVert];
    FLOAT average_coeff[2];

   ForAllElements(g, e) {
	bzero(A[0], NVert * NVert * sizeof(A[0][0]));
	bzero(local_v, NVert * sizeof(local_v[0]));
	bzero(rhs, NVert * sizeof(rhs[0]));

	if (device == MOS && e->region_mark == OXIDE) {
	    for (i = 0; i < NVert; i++) {
		rhs[i] = 0.;
		if (!(verts_region[e->verts[i]] & MAT_SEMICONDUCTOR))
		    A[i][i] = 1.;
		else
		    A[i][i] = 0.;
	    }
	    goto add_entries;
	}

	Doping(e->region_mark, &ND, &NA);
	vol = phgGeomGetVolume(g, e);
      for (i = 0; i < NVert; i++)
      {
         u[i] = *(DofVertexData(dof, e->verts[i]));
	      bas_val[i] = phgQuadGetBasisValues(e, dof, i, quad);
         local_v[i] = *(DofVertexData(V, e->verts[i]));     
      }

	coeff = eqn == ELECTRON ? Dn : Dp;
	coeff /= CM2MICRON;	/* handle units */

if(version_eafe == 1){
      for(i = 0; i < NVert; i++)
      {
         //psi0 = *(DofVertexData(V, e->verts[i]));
         psi0 = local_v[i];
         A[i][i] = 0.0;
         for(j = 0; j < NVert; j++)
         {
            if(j != i)
            {
            //   psi1 = *(DofVertexData(V, e->verts[j]));
               psi1 = local_v[j];
               B_ij = eqn == HOLE ? B(psi0 - psi1) : B(psi1 - psi0);
               B_ii = eqn == HOLE ? B(psi1 - psi0) : B(psi0 - psi1);
               A_bas = phgQuadGradBasDotGradBas(e, dof, j, dof, i, 3);
               A[i][j] =  coeff * B_ij * A_bas;
               A[i][i] -= coeff * B_ii * A_bas;
            }
         }
      }
   }
   else if(version_eafe == 2){
      for(i = 0; i < NVert; i++)
      {
         //psi0 = *(DofVertexData(V, e->verts[i]));
         psi0 = local_v[i];
         A[i][i] = 0.0;
         for(j = 0; j < NVert; j++)
         {
            if(j != i)
            {
               //psi1 = *(DofVertexData(V, e->verts[j]));
               psi1 = local_v[j];
               B_ii = eqn == HOLE ? Average2(psi0 - psi1) : Average2(psi1 - psi0);
               B_ij = 1.0 - B_ii;
               A_bas = phgQuadGradBasDotGradBas(e, dof, j, dof, i, 3);
               A[i][j] =  coeff * 2 * B_ij * A_bas;
               A[i][i] -= coeff * 2 * B_ii * A_bas;
            }
         }
      }
   }
else if(version_eafe == 3){
      for(i = 0; i < NVert; i++)
      {
         //psi0 = *(DofVertexData(V, e->verts[i]));
         psi0 = local_v[i];
         A[i][i] = 0.0;
         for(j = 0; j < NVert; j++)
         {
            if(j != i)
            {
              // psi1 = *(DofVertexData(V, e->verts[j]));
               psi1 = local_v[j];
               B_ij = eqn == HOLE ? Exp(psi1 - psi0) : Exp(psi0 - psi1);
               B_ii = eqn == HOLE ? Exp(psi0 - psi1) : Exp(psi1 - psi0);
               A_bas = phgQuadGradBasDotGradBas(e, dof, j, dof, i, 3);
               A[i][j] =  coeff / 2 * (1 + B_ij) * A_bas;
               A[i][i] -= coeff / 2 * (1 + B_ii) * A_bas;
            }
         }
      }
   }
else if(version_eafe == 4){
      for(i = 0; i < NVert; i++)
      {
         A[i][i] = 0.0;
            for(j = 0; j < NVert; j++)
            {
               if(j != i)
               {
	               bzero(average_coeff, 2 * sizeof(average_coeff[0]));
                  average(local_v, eqn, average_coeff, i, j);
                  A_bas = phgQuadGradBasDotGradBas(e, dof, j, dof, i, 3);
                  A[i][j] = coeff / average_coeff[1] * A_bas;
                  A[i][i] -= coeff / average_coeff[0] * A_bas;
               }
            }
      }
   }
else if(version_eafe == 5)
{
      for(i = 0; i < NVert; i++)
      {
         psi0 = local_v[i];
         A[i][i] = 0.0;
         for(j = 0; j < NVert; j++)
         {
            if(j != i)
            {
               psi1 = local_v[j];
               B_ii = eqn == HOLE ? Average3_0(psi0 - psi1) : Average3_0(psi1 - psi0);
               B_ij = eqn == HOLE ? Average3_1(psi0 - psi1) : Average3_1(psi1 - psi0);
               A_bas = phgQuadGradBasDotGradBas(e, dof, j, dof, i, 3);
               A[i][j] =  coeff *  B_ij * A_bas;
               A[i][i] -= coeff *  B_ii * A_bas;
            }
         }
      }
}
      for(i = 0; i < NVert; i++)
      {
         r = 0.0;
         for(j = 0; j < NVert; j++)
            r += A[i][j] * u[j]; 
         rhs[i] = -1.0 * r;
      }

	/* handle bdry */
	for (i = 0; i < NVert; i++) {
	    BTYPE btype = phgDofGetElementBoundaryType(dof, e, i);
	    if (btype & ohmic_contact) {
		   for (j = 0; j < NVert; j++)
		       A[i][j] =  A[j][i] = 0.0;
		   rhs[i] = 0.0;
		   A[i][i] = 1.0;
	    }
	}

      add_entries:
	/* add entries */
	for (i = 0; i < NVert; i++)
	    I[i] = phgSolverMapE2L(solver, 0, e, i);

	phgSolverAddMatrixEntries(solver, NVert, I, NVert, I, A[0]);
	phgSolverAddRHSEntries(solver, NVert, I, rhs);
   }
}

static void
build_linear_system_continuity_zlamal(SOLVER *solver, DOF *V, DOF *P, DOF *N,
				      EQUATION eqn)
{
    GRID *g = V->g;
    SIMPLEX *e;
    DOF *dof = eqn == ELECTRON ? N : P;
    DOF *gradv = phgDofGradient(V, NULL, NULL, NULL);
    FLOAT A[NVert][NVert], rhs[NVert], expd[NVert];
    FLOAT u[NVert], psi[NVert], du[Dim], f[Dim];
    INT I[NVert];
    int i, j, k, n;
    int order = 1;
    QUAD *quad = phgQuadGetQuad3D(order);
    FLOAT fr[quad->npoints], fl[quad->npoints];
    const FLOAT* bas_val[NVert];
    FLOAT ND, NA, vol, r, coeff, psi0, psi1;
    const FLOAT *w, *dphi[NVert];

    ForAllElements(g, e) {
	bzero(A[0], NVert * NVert * sizeof(A[0][0]));
	bzero(rhs, NVert * sizeof(rhs[0]));

	if (device == MOS && e->region_mark == OXIDE) {
	    for (i = 0; i < NVert; i++) {
		rhs[i] = 0.;
		if (!(verts_region[e->verts[i]] & MAT_SEMICONDUCTOR))
		    A[i][i] = 1.;
		else
		    A[i][i] = 0.;
	    }
	    goto add_entries;
	}

	Doping(e->region_mark, &ND, &NA);
	vol = phgGeomGetVolume(g, e);
	for (i = 0; i < NVert; i++) {
	    dphi[i] = phgQuadGetBasisGradient(e, dof, i, quad);
	    u[i] = *(DofVertexData(dof, e->verts[i]));
	    psi[i] = *(DofVertexData(V, e->verts[i]));
	    bas_val[i] = phgQuadGetBasisValues(e, V, i, quad);
	}

	if (recombination) {
	    FLOAT valp[quad->npoints], valn[quad->npoints];
	    exp_quad_dof_eval(e, P, HOLE, gradv, quad, valp);
	    exp_quad_dof_eval(e, N, ELECTRON, gradv, quad, valn);
	    //memcpy(valp, phgQuadGetDofValues(e, P, quad), quad->npoints * sizeof(valp[0]));
	    //memcpy(valn, phgQuadGetDofValues(e, N, quad), quad->npoints * sizeof(valn[0]));
	    const FLOAT *v = phgQuadGetDofValues(e, V, quad);
	    for (i = 0; i < quad->npoints; i++)
		SRH(v[i], valp[i], valn[i], nie, eqn, &(fl[i]), &(fr[i]));
	}

	coeff = eqn == ELECTRON ? Dn : Dp;
	coeff /= CM2MICRON;	/* handle units */

	for (i = 0; i < NVert; i++) {
	    FLOAT (*J)[Dim] = (void *)JBJ(V, e, i, eqn);

	    for (j = 0; j < NVert; j++) {
		expd[j] = eqn == ELECTRON ?
		    Exp(psi[i] - psi[j]) : Exp(psi[j] - psi[i]);
	    }

	    w = quad->weights;
	    for (n = 0; n < quad->npoints; n++) {
		du[0] = du[1] = du[2] = 0.;
		for (j = 0; j < NVert; j++) {
		    for (k = 0; k < Dim; k++) {
			du[k] += u[j] * expd[j] * dphi[j][n * Dim + k];
		    }
		}

		r = 0.;
		for (k = 0; k < Dim; k++) {
		    f[k] =
			J[k][0] * du[0] + J[k][1] * du[1] + J[k][2] * du[2];
		    r += f[k] * dphi[i][n * Dim + k];
		}
		r *= *w * vol;
		rhs[i] -= coeff * r;	/* residual */

		for (j = 0; j < NVert; j++) {
		    for (k = 0; k < Dim; k++) {
			f[k] =
			    J[k][0] * dphi[j][n * Dim + 0] +
			    J[k][1] * dphi[j][n * Dim + 1] +
			    J[k][2] * dphi[j][n * Dim + 2];
		    }

		    r = f[0] * dphi[i][n * Dim + 0] +
			f[1] * dphi[i][n * Dim + 1] +
			f[2] * dphi[i][n * Dim + 2];
		    r *= *w * vol;
		    A[i][j] += coeff * expd[j] * r;
		}

		if(recombination){
		    FLOAT l = 1./(CM2MICRON * CM2MICRON * CM2MICRON);
		    for(j=0;j<NVert;j++){
			A[i][j] +=
			    l * fl[n] * bas_val[i][n] * bas_val[j][n] * (*w) * vol;
		    }
		    rhs[i] += l * fr[n] * bas_val[i][n] * (*w) * vol;
		}

		w++;
	    }
	}

	/* handle bdry */
	for (i = 0; i < NVert; i++) {
	    BTYPE btype = phgDofGetElementBoundaryType(V, e, i);
	    if (btype & ohmic_contact) {
		for (j = 0; j < NVert; j++)
		    A[i][j] = A[j][i] = 0.0;
		rhs[i] = 0.0;
		A[i][i] = 1.0;
	    }
	}

      add_entries:
	/* add entries */
	for (i = 0; i < NVert; i++)
	    I[i] = phgSolverMapE2L(solver, 0, e, i);

	phgSolverAddMatrixEntries(solver, NVert, I, NVert, I, A[0]);
	phgSolverAddRHSEntries(solver, NVert, I, rhs);
    }

    phgDofFree(&gradv);
}

RelativeError(DOF *V, DOF *P, DOF *N, BOOLEAN flag)
{
    GRID *g = V->g;
    static DOF *V_old, *P_old, *N_old;
    static FLOAT err[3];
    INT i;

    if (!flag) {
	V_old = phgDofCopy(V, NULL, NULL, NULL);
	P_old = phgDofCopy(P, NULL, NULL, NULL);
	N_old = phgDofCopy(N, NULL, NULL, NULL);
	return NULL;
    }
    else {
	FLOAT v_old, v_new, p_old, p_new, n_old, n_new, r, rv, rp, rn;
	FLOAT C;

#if 0
	C = 1.e14 / SCALE;
#else
	C = 1.;
#endif

#define max(a,b) ((a)>(b)? (a):(b))

	rv = rp = rn = 0.;
	for (i = 0; i < g->nvert; i++) {
	    if (g->types_vert[i] == UNREFERENCED)
		continue;

	    v_old = *(V_old->data + i);
	    v_new = *(V->data + i);
	    r = Fabs(v_new - v_old) / max(Fabs(v_old), 1.);
	    if (rv < r)
		rv = r;

	    if (!(verts_region[i] & MAT_SEMICONDUCTOR))
		continue;

	    p_old = *(P_old->data + i);
	    p_new = *(P->data + i);
	    n_old = *(N_old->data + i);
	    n_new = *(N->data + i);

	    if (p_old < 1. / SCALE)
		p_old = 1. / SCALE;
	    if (p_new < 1. / SCALE)
		p_new = 1. / SCALE;
	    if (n_old < 1. / SCALE)
		n_old = 1. / SCALE;
	    if (n_new < 1. / SCALE)
		n_new = 1. / SCALE;

	    /* concentration to fermi level */
	    p_old = v_old + Log(p_old * SCALE / nie);
	    n_old = v_old - Log(n_old * SCALE / nie);
	    p_new = v_new + Log(p_new * SCALE / nie);
	    n_new = v_new - Log(n_new * SCALE / nie);

	    r = Fabs(p_new - p_old) / max(Fabs(p_old), C);
	    if (rp < r)
		rp = r;
	    r = Fabs(n_new - n_old) / max(Fabs(n_old), C);
	    if (rn < r)
		rn = r;
	}
#undef max

	err[0] = rv;
	err[1] = rp;
	err[2] = rn;

#if USE_MPI
	if (g->nprocs > 1) {
	    FLOAT err0[3];
	    memcpy(err0, err, 3 * sizeof(err0[0]));
	    MPI_Allreduce(err0, err, 3, PHG_MPI_FLOAT, MPI_MAX, g->comm);
	}
#endif

	phgDofFree(&V_old);
	phgDofFree(&P_old);
	phgDofFree(&N_old);

	return err;
    }
}

static double
elapsed_time(GRID *g, BOOLEAN flag, SOLVER *solver)
/* returns and prints elapsed time since last call to this function */
{
    static double t0 = 0.0;
    double et, tt[3];
    size_t mem;

    phgGetTime(tt);
    et = tt[2] - t0;
    t0 = tt[2];

    mem = phgMemoryUsage(g, NULL);

    if (flag) {
	if (solver != NULL) {
	    phgPrintf("  Linear systems solved: [nits %d, res %0.2le ",
		      solver->nits, (double)solver->residual);
	}
	phgPrintf("%0.4lgMB %0.4lfs]\n", mem / (1024.0 * 1024.0), et);
    }

    return et;
}

static void
Update(DOF *V, DOF *P, DOF *N, DOF *delta_v, FLOAT damp)
/* to keep the Quasi-Fermi potential constant when solving Poisson's equation */
{
    INT i;
    GRID *g = V->g;
    FLOAT dv;

    phgDofAXPBY(damp, delta_v, 1.0, &V);

    for (i = 0; i < g->nvert; i++) {
	if (g->types_vert[i] == UNREFERENCED ||
	    !(verts_region[i] & MAT_SEMICONDUCTOR))
	    continue;

	dv = *(DofVertexData(delta_v, i)) * damp;
	*(P->data + i) *= Exp(-dv);
	*(N->data + i) *= Exp(dv);
    }
}

static void
Current(DOF *V, DOF *P, DOF *N, CONTACT cont[], int ncont,
	FLOAT currents[][3])
{
    if (!strcmp(continuity_names[continuity_index], "zlamal"))
	current_zlamal(V, P, N, cont, ncont, currents);
    else if (!strcmp(continuity_names[continuity_index], "eafe"))
	   current_eafe(V, P, N, cont, ncont, currents);
    else
	phgError(1, "%s:%d,unexpected error.\n", __FILE__, __LINE__);
}

static void
SolvePoisson(DOF *V, DOF *P, DOF *N, FLOAT damp)
{
    GRID *g = V->g;
    DOF *delta;
    SOLVER *solver;
    int iter_newton;
    int maxit_newton;
    FLOAT r, rtol_v = 1.0e-5;

    if (damp < 0.5)
	maxit_newton = 10000;	/* it can take many iterations to get to the 
				   equilibrium state for MOS */
    else
	maxit_newton = 100;

    delta = phgDofNew(g, V->type, V->dim, "delta", DofNoAction);

    phgOptionsPush();
#if 0
    //poisson
    phgOptionsSetOptions("-solver  mumps -mumps_symmetry spd");
#elif 1				/* iterative method */
    phgOptionsSetOptions
	("-solver hypre -hypre_solver pcg -hypre_pc boomeramg");
    phgOptionsSetOptions
	("-solver_rtol 1.0e-15 -solver_atol 1.0e-15 -solver_btol 1.0e-15");
    phgOptionsSetOptions("-solver_maxit 200");
#else /* __float128 */
   phgOptionsSetOptions
	("-solver pcg -pcg_pc_type solver -pcg_pc_opts \"-solver mumps -mumps_symmetry spd\"");
    phgOptionsSetOptions
	("-solver_rtol 1.0e-30 -solver_atol 1.0e-30 -solver_btol 1.0e-30");
    phgOptionsSetOptions("-solver_maxit 5");
#endif

    for (iter_newton = 1; iter_newton <= maxit_newton; iter_newton++) {
	solver = phgSolverCreate(SOLVER_DEFAULT, delta, NULL);
	solver->mat->handle_bdry_eqns = FALSE;

	if (!strcmp(poisson_names[poisson_index], "fem"))
	    build_linear_system_poisson_fem(solver, V, P, N);
	else
	    phgError(1, "%s:%d,unexpected error.\n", __FILE__, __LINE__);

	phgDofSetDataByValue(delta, 0.0);
	elapsed_time(g, FALSE, NULL);
	phgSolverSolve(solver, TRUE, delta, NULL);
	elapsed_time(g, TRUE, solver);
	phgSolverDestroy(&solver);
	Update(V, P, N, delta, damp);

	r = phgDofNormInftyVec(V);
	r = r > 1.0 ? r : 1.0;
	r = phgDofNormInftyVec(delta) / r;
	phgPrintf("  rtol %0.6le\n", (double)r);
	if (r < rtol_v) {
	    phgPrintf("  Converged.\n");
	    break;
	}
    }

    if (iter_newton > maxit_newton)
	phgError(1, "%s:%d,unexpected error.\n", __FILE__, __LINE__);

    phgOptionsPop();
    phgDofFree(&delta);
}

static void
SolveContinuity(DOF *V, DOF *P, DOF *N, EQUATION eqn)
{
    GRID *g = V->g;
    SOLVER *solver;
    DOF *u;
    int iter_newton;
    int maxit_newton = 1;
    FLOAT r, rtol_v = 1.0e-5;

    assert(eqn == ELECTRON || eqn == HOLE);
    u = eqn == ELECTRON ? N : P;

    phgOptionsPush();
//#if 0

//    phgOptionsSetOptions("-solver mumps -mumps_symmetry unsym");
//#elif 1				/* iterative method */
  //  phgOptionsSetOptions("-solver petsc");
  //  phgOptionsSetOptions("-oem_options \"-ksp_type gmres\"");
 phgOptionsSetOptions
	("-petsc_pc_opts \"-solver hypre -hypre_solver boomeramg \
	    -hypre_amg_coarsen_type falgout -hypre_pc none\"");
    phgOptionsSetOptions("-solver_maxit 5000");
    phgOptionsSetOptions
	("-solver_rtol 1.0e-15 -solver_atol 1.0e-15 -solver_btol 1.0e-15");
//#elif 0
///    phgOptionsSetOptions("-solver petsc");
//   phgOptionsSetOptions("-oem_options \"-ksp_type gmres\"");
//   phgOptionsSetOptions("-oem_options \"-pc_asm_overlap 1\"");
//   phgOptionsSetOptions
//	("-oem_options \"-sub_pc_factor_mat_solver_package mumps\"");
//    phgOptionsSetOptions("-oem_options \"-mat_mumps_sym 0\"");
//    phgOptionsSetOptions("-oem_options \"-mat_mumps_icntl_4 0\"");
//   phgOptionsSetOptions("-solver_maxit 1000");
//   phgOptionsSetOptions
//	("-solver_rtol 1.0e-12 -solver_atol 1.0e-12 -solver_btol 1.0e-12");
//#else /* __float128 */
//   phgOptionsSetOptions
//	("-solver gmres -gmres_pc_type solver -gmres_pc_opts \"-solver mumps -mumps_symmetry unsym\"");
//    phgOptionsSetOptions
//	("-solver_rtol 1.0e-30 -solver_atol 1.0e-30 -solver_btol 1.0e-30");
//    phgOptionsSetOptions("-solver_maxit 5");
//#endif
    DOF *delta = phgDofNew(g, u->type, u->dim, "increment", DofNoAction);

    for (iter_newton = 1; iter_newton <= maxit_newton; iter_newton++) {
	solver = phgSolverCreate(SOLVER_DEFAULT, delta, NULL);
	solver->mat->handle_bdry_eqns = FALSE;
	if (!strcmp(continuity_names[continuity_index], "zlamal"))
	    build_linear_system_continuity_zlamal(solver, V, P, N, eqn);
   else if (!strcmp(continuity_names[continuity_index], "eafe"))
	    build_linear_system_continuity_eafe(solver, V, P, N, eqn);
	else
	    phgError(1, "%s:%d,unexpected error.\n", __FILE__, __LINE__);


	phgDofSetDataByValue(delta, 0.0);
	elapsed_time(g, FALSE, NULL);
	phgSolverSolve(solver, TRUE, delta, NULL);
        phgPrintf("$$$$$$$");
	elapsed_time(g, TRUE, solver);
	phgSolverDestroy(&solver);
#if 0
   FLOAT udamp;
   if(anode - bias_cathode < 1.2)
      udamp = eqn == HOLE ? 1.0 : 0.5;
   else
      udamp = eqn == HOLE ? 0.5 : 0.1;

	phgDofAXPBY(udamp, delta, 1.0, &u);
#elif 1
	phgDofAXPBY(1.0, delta, 1.0, &u);
#endif
    }

    phgDofFree(&delta);
    phgOptionsPop();
}

#undef r2
#undef B

static int
bc_map(int bctype)
{
    if (device == PN) {
	switch (bctype) {
	    case 1:
		return CATHODE;
	    case 2:
		return ANODE;
	    case 4:
		return JUNCTION_PN;
	    default:
		return UNDEFINED;
	}
    }
    else if (device == NPN) {
	switch (bctype) {
	    case 2:
		return EMITTER;
	    case 3:
		return COLLECTOR;
	    case 4:
		return BASE;
	    case 8:
		return JUNCTION_CB;
	    case 9:
		return JUNCTION_EB;
	    default:
		return UNDEFINED;
	}
    }
    else if (device == PNP) {
	switch (bctype) {
	    case 2:
		return EMITTER;
	    case 3:
		return COLLECTOR;
	    case 4:
		return BASE;
	    case 8:
		return JUNCTION_CB;
	    case 9:
		return JUNCTION_EB;
	    default:
		return UNDEFINED;
	}
    }
    else if (device == MOS) {
	switch (bctype) {
	    case 2:
		return SOURCE;
	    case 3:
		return DRAIN;
	    case 4:
		return SUBSTRATE;
	    case 5:
		return JUNCTION_OS;
	    case 6:
		return GATE;
	    case 8:
		return JUNCTION_SS;
	    case 9:
		return JUNCTION_DS;
	    default:
		return UNDEFINED;
	}
    }

    return UNDEFINED;
}

int
main(int argc, char *argv[])
{
    GRID *g;
    DOF *V, *P, *N;
    size_t mem, mem_peak;
    const char *fn = NULL;
    int pre_refines = 0;

    int iter_fp, maxit_fp = 10000;
    const FLOAT *err;
    FLOAT rtol = 1.0e-4, rtol_v = 1.0e-4;

    phgOptionsRegisterKeyword("poisson",
			      "discretization methods for Poisson's equations",
			      poisson_names, &poisson_index);
    phgOptionsRegisterKeyword("continuity",
			      "discretization methods for continuity equations",
			      continuity_names, &continuity_index);
    phgOptionsRegisterFilename("mesh_file", "Mesh file", (char **)&fn);
    phgOptionsRegisterInt("pre_refines", "Pre-refines", &pre_refines);
    phgOptionsRegisterInt("recombination", "recombination", &recombination);
    phgOptionsRegisterInt("version_eafe", "eafe version", &version_eafe);
    phgOptionsRegisterInt("device", "device", &device);

    /* PN */
    phgOptionsRegisterFloat("bias_cathode", "bias_cathode", &bias_cathode);
    phgOptionsRegisterFloat("bias_anode", "bias_anode", &bias_anode);
    phgOptionsRegisterFloat("doping_pregion", "doping_pregion",
			    &doping_pregion);
    phgOptionsRegisterFloat("doping_nregion", "doping_nregion",
			    &doping_nregion);
    phgOptionsRegisterFloat("h", "bias increment", &h);

    /* NPN/PNP */
    phgOptionsRegisterFloat("bias_emitter", "bias_emitter", &bias_emitter);
    phgOptionsRegisterFloat("bias_collector", "bias_collector",
			    &bias_collector);
    phgOptionsRegisterFloat("bias_base", "bias_base", &bias_base);
    phgOptionsRegisterFloat("h_b", "bias increment", &h_b);
    phgOptionsRegisterFloat("h_c", "bias increment", &h_c);
    phgOptionsRegisterFloat("h_e", "bias increment", &h_e);
    phgOptionsRegisterFloat("doping_e", "doping_e", &doping_e);
    phgOptionsRegisterFloat("doping_c", "doping_c", &doping_c);
    phgOptionsRegisterFloat("doping_b", "doping_b", &doping_b);

    /* MOS */
    phgOptionsRegisterFloat("bias_source", "bias_source", &bias_source);
    phgOptionsRegisterFloat("bias_drain", "bias_drain", &bias_drain);
    phgOptionsRegisterFloat("bias_gate", "bias_gate", &bias_gate);
    phgOptionsRegisterFloat("bias_substrate", "bias_substrate",
			    &bias_substrate);
    phgOptionsRegisterFloat("h_d", "bias increment", &h_d);
    phgOptionsRegisterFloat("h_g", "bias increment", &h_g);
    phgOptionsRegisterFloat("doping_s", "doping_s", &doping_s);
    phgOptionsRegisterFloat("doping_d", "doping_d", &doping_d);
    phgOptionsRegisterFloat("doping_st", "doping_st", &doping_st);

    /** pre-set options **/
#if 0
    phgOptionsPreset("-device 0");	/* 0 --> PN */
    phgOptionsPreset("-doping_pregion 1.0e17 -doping_nregion 1.0e17");
    phgOptionsPreset("-mesh_file ./mesh/pn10w.mesh");
    phgOptionsPreset("-bias_anode -5.0 -h 0.05");
    phgOptionsPreset("-pre_refines 0");
    phgOptionsPreset("-recombination 0");
    phgOptionsPreset("-poisson fem");
  //  phgOptionsPreset("-continuity zlamal");
    phgOptionsPreset("-continuity eafe");
    phgOptionsPreset("-version_eafe 1");
#endif
#if 1
    phgOptionsPreset("-device 3");	/* 3 --> MOS */
    phgOptionsPreset("-doping_s 1.0e18 -doping_d 1.0e18 -doping_st 1.0e16");
    phgOptionsPreset("-bias_gate 2.0 -bias_drain 0.5");
    phgOptionsPreset("-h_g 0.05 -h_d 0.05");
    phgOptionsPreset("-mesh_file ./mesh/mos8w.mesh");
    phgOptionsPreset("-pre_refines 0");
    phgOptionsPreset("-recombination 0");
    phgOptionsPreset("-poisson fem");
   // phgOptionsPreset("-continuity zlamal");
    phgOptionsPreset("-continuity eafe");
    phgOptionsPreset("-version_eafe 1");
#endif
#if 0
    phgOptionsPreset("-device 1");	/* 1 --> NPN */
    phgOptionsPreset("-doping_e 1.0e17 -doping_c 2.0e16 -doping_b 2.0e14");
    phgOptionsPreset("-bias_collector 2.5 -bias_base 6.0");/*emitter grounded*/
    phgOptionsPreset("-h_b 0.05 -h_c 0.05");/*emitter grounded*/
    //phgOptionsPreset("-bias_collector 3.0 -bias_emitter -10.0");/*base grounded*/
    //phgOptionsPreset("-h_e 0.05 -h_c 0.05");/*base grounded*/
    phgOptionsPreset("-mesh_file ./mesh/GLPNP_nOXIDE.mesh");
    phgOptionsPreset("-pre_refines 0");
    phgOptionsPreset("-recombination 0");
    phgOptionsPreset("-poisson fem");
  //  phgOptionsPreset("-continuity zlamal");
    phgOptionsPreset("-continuity eafe");
    phgOptionsPreset("-version_eafe 1");
#endif
#if 0
    phgOptionsPreset("-device 2");	/* 2 --> PNP */
    //phgOptionsPreset("-doping_e 1.0e18 -doping_c 1.0e17 -doping_b 1.0e14");
    phgOptionsPreset("-doping_e 1.0e17 -doping_c 2.0e16 -doping_b 2.0e14");
    //phgOptionsPreset("-doping_e 1.0e16 -doping_c 2.0e15 -doping_b 2.0e13");
    phgOptionsPreset("-bias_collector -1.0 -bias_base -1.0");/*emitter grounded*/
    phgOptionsPreset("-h_b 0.05 -h_c 0.05");/*emitter grounded*/
    //phgOptionsPreset("-bias_collector -0.5 -bias_emitter 5.0");/*base grounded*/
    //phgOptionsPreset("-h_e 0.05 -h_c 0.05");/*base grounded*/
    phgOptionsPreset("-mesh_file ./mesh/GLPNP_nOXIDE.mesh");
    phgOptionsPreset("-pre_refines 0");
    phgOptionsPreset("-recombination 0");
    phgOptionsPreset("-poisson fem");
    //phgOptionsPreset("-continuity zlamal");
    phgOptionsPreset("-continuity eafe");
    phgOptionsPreset("-version_eafe 1");
#endif
    
phgOptionsPreset("-solver petsc -solver_maxit 100 -solver_rtol 1e-2 "
		     "-oem_options \"-ksp_type gmres -ksp_gmres_restart 50\" "
		     /*"-oem_options \"-ksp_monitor\""*/);
//  phgOptionsPreset("-solver petsc  -oem_options {-ksp_type gmres -ksp_gmres_restart 50 -ksp_pc_side right -pc_type asm -pc_asm_type restrict -pc_asm_overlap 2 -sub_ksp_type preonly -sub_pc_type lu }");
    phgInit(&argc, &argv);
    g = phgNewGrid(-1);		/* WARNING: no idle processes are allowed */
    phgImportSetBdryMapFunc(bc_map);
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);
    phgRefineAllElements(g, pre_refines > 6 ? 6 : pre_refines);
    phgBalanceGrid(g, 1.2, -1, NULL, 0.);
    if (pre_refines > 6) {
	phgRefineAllElements(g, pre_refines - 6);
	phgBalanceGrid(g, 1.2, -1, NULL, 0.);
    }

    InitConstant();
    PrintDopingProfile();

    phgPrintf
	("  %d vertices, %d elements, %d submeshes, load imbalance: %lg\n",
	 g->nvert_global, g->nleaf_global, g->nprocs, (double)g->lif);

    V = phgDofNew(g, DOF_P1, 1, "electric potential", DofInterpolation);
    P = phgDofNew(g, DOF_P1, 1, "hole concentration", DofInterpolation);
    N = phgDofNew(g, DOF_P1, 1, "electron concentration", DofInterpolation);

 //   Center(g, &tet_center, &tri_center, &volume);
    SetupVertsRegion(g, TRUE);
    /* initialize the device */
    if (device == PN) {
	anode = cathode = bias_cathode;	/* bias_cathode grounded */
	ohmic_contact = CATHODE | ANODE;
	phgPrintf("\n *********** \n");
    }
    else if (device == NPN) {
	collector = base = emitter = bias_emitter;	/* emitter grounded */
	//collector = base = emitter = bias_base;	/* base grounded */
	ohmic_contact = EMITTER | COLLECTOR | BASE;
    }
    else if (device == PNP) {
	collector = base = emitter = bias_emitter;	/* emitter grounded */
	//collector = base = emitter = bias_base;	/* base grounded */
	ohmic_contact = EMITTER | COLLECTOR | BASE;
    }
    else if (device == MOS) {	/* source and substrate grounded */
	gate = substrate;
   gate = -2.0;
	source = drain = substrate = bias_substrate;
//	source = substrate = bias_substrate;
//	drain = bias_drain;
	ohmic_contact = SOURCE | DRAIN | SUBSTRATE;
    }

    InitialGuess(V, P, N);

    phgPrintf("\n  ================ EQUILIBRIUM STATE ================\n\n");
    SolvePoisson(V, P, N, 0.2);	/* to get thermal equilibrium state */

#if 0
    /* rescale the data for output */
    INT i;
    for (i = 0; i < g->nvert; i++) {
	if (g->types_vert[i] == UNREFERENCED)
	    continue;
	*(V->data + i) *= KB * T;

	if (!(verts_region[i] & MAT_SEMICONDUCTOR))
	    continue;

	*(P->data + i) *= SCALE;
	*(N->data + i) *= SCALE;
    }

   char vtk_name[40];
   sprintf(vtk_name, "pn%0.2f.vtk", anode);
//    phgPrintf("\n  \"%s\" created.\n", phgExportVTK(g, "pn.vtk", V, P, N, NULL));
    phgPrintf("\n  \"%s\" created.\n", phgExportVTK(g, vtk_name, V, P, N, NULL));

    /* rescale the data for input */
    for (i = 0; i < g->nvert; i++) {
	if (g->types_vert[i] == UNREFERENCED)
	    continue;
	*(V->data + i) /= KB * T;

	if (!(verts_region[i] & MAT_SEMICONDUCTOR))
	    continue;

	*(P->data + i) /= SCALE;
	*(N->data + i) /= SCALE;
    }
#endif
    while (TRUE) {		/* bias ramp */
	if (BiasSweep())
	    break;

	TakeBoundary(V, P, N, 0);

	/* Gummel iteration for each bias */
	for (iter_fp = 1; iter_fp <= maxit_fp; iter_fp++) {
	    RelativeError(V, P, N, FALSE);

	    phgPrintf("\n  ---- The %d-th Gummel iteration ----\n", iter_fp);

	    phgPrintf("\n  The Poisson equation:\n");
	     SolvePoisson(V, P, N, 1.0);

//	    phgPrintf("\n  The electron continuity equation:\n");
//	    SolveContinuity(V, P, N, ELECTRON);

	    phgPrintf("\n  The hole continuity equation:\n");
	    SolveContinuity(V, P, N, HOLE);

            phgPrintf("\n  The electron continuity equation:\n");
            SolveContinuity(V, P, N, ELECTRON);

	    /* convergence criteria for gummel iterations */
	    err = RelativeError(V, P, N, TRUE);
	    phgPrintf("\n  rtol: V %0.6le, P %0.6le, N %0.6le\n",
		      (double)err[0], (double)err[1], (double)err[2]);

	    if (err[0] < rtol_v && err[1] < rtol && err[2] < rtol) {
		phgPrintf("  Converged!\n");
		break;
	    }
	}

	if (iter_fp > maxit_fp)
	    phgError(1, "  Fail to converge.\n");

	/* compute terminal currents */
	if (device == PN) {
	    CONTACT cont[2] = { ANODE, CATHODE };
	    FLOAT currents[2][3];
	    Current(V, P, N, cont, 2, currents);
	    phgPrintf("\n  Currents (V_a: %0.4lfV, V_c: %0.4lfV): ",
		      (double)anode, (double)cathode);
	    phgPrintf("anode: %0.6le, cathode %0.6le\n",
		      (double)(Q * currents[0][0]),
		      (double)(Q * currents[1][0]));
	}
	else if (device == NPN) {
	    CONTACT cont[3] = { EMITTER, COLLECTOR, BASE };
	    FLOAT currents[3][3];
	    Current(V, P, N, cont, 3, currents);
	    phgPrintf
		("\n  Currents (V_c: %0.4lf V_b: %0.4lf, V_e: %0.4lf): ",
		 (double)collector, (double)base, (double)emitter);
	    phgPrintf("emitter: %0.4le, collector: %0.4le, base %0.4le\n",
		      (double)(Q * currents[0][0]),
		      (double)(Q * currents[1][0]),
		      (double)(Q * currents[2][0]));
	}
	else if (device == PNP) {
	    CONTACT cont[3] = { EMITTER, COLLECTOR, BASE };
	    FLOAT currents[3][3];
	    Current(V, P, N, cont, 3, currents);
	    phgPrintf
		("\n  Currents (V_c: %0.4lf V_b: %0.4lf, V_e: %0.4lf): ",
		 (double)collector, (double)base, (double)emitter);
	    phgPrintf("emitter: %0.4le, collector: %0.4le, base %0.4le\n",
		      (double)(Q * currents[0][0]),
		      (double)(Q * currents[1][0]),
		      (double)(Q * currents[2][0]));
	}
	else if (device == MOS) {
	    CONTACT cont[3] = { SOURCE, DRAIN, SUBSTRATE };
	    FLOAT currents[3][3];
	    Current(V, P, N, cont, 3, currents);
	    phgPrintf
		("\n  Currents (V_g: %0.4lf V_d: %0.4lf): ", (double)gate,
		 (double)drain);
	    phgPrintf("source: %0.4le, drain: %0.4le, substrate: %0.4le\n",
		      (double)(Q * currents[0][0]),
		      (double)(Q * currents[1][0]),
		      (double)(Q * currents[2][0]));


	}

}				/* end bias ramp */

    INT i;
#if 0 
//if (gate + 1.0e-6 > bias_gate) {
//if (Fabs(drain) + 1.0e-6 > Fabs(bias_drain)) {
    /* rescale the data for output */
    for (i = 0; i < g->nvert; i++) {
	if (g->types_vert[i] == UNREFERENCED)
	    continue;
	*(V->data + i) *= KB * T;

	if (!(verts_region[i] & MAT_SEMICONDUCTOR))
	    continue;

	*(P->data + i) *= SCALE;
	*(N->data + i) *= SCALE;
    }

//#if 0
   char vtk_name[40];
   if(device == PN)
      sprintf(vtk_name, "pn%0.2f.vtk", anode);
   else if(device == MOS)
      sprintf(vtk_name, "mos_g%0.2f_d%0.2f_NOT5e11.vtk", gate, drain);
   else if(device == PNP)
      sprintf(vtk_name, "pnp_c%0.2f_b%0.2f_e%0.2f.vtk",collector, base, emitter);
    phgPrintf("\n  \"%s\" created.\n", phgExportVTK(g, vtk_name, V, P, N, NULL));
//#endif
    /* rescale the data for input */
   // INT i;
    for (i = 0; i < g->nvert; i++) {
	if (g->types_vert[i] == UNREFERENCED)
	    continue;
	*(V->data + i) /= KB * T;

	if (!(verts_region[i] & MAT_SEMICONDUCTOR))
	    continue;

	*(P->data + i) /= SCALE;
	*(N->data + i) /= SCALE;
    }
//}
#endif

//    phgDofFree(&tet_center);
//    phgDofFree(&tri_center);
//    phgDofFree(&volume);
//    SetupVertsRegion(g, FALSE);

    mem = phgMemoryUsage(g, &mem_peak);
    phgPrintf("\n  Memory usage: current %0.4lgMB, peak %0.4lgMB\n",
	      (double)mem / (1024.0 * 1024.0),
	      (double)mem_peak / (1024.0 * 1024.0));

    phgDofFree(&V);
    phgDofFree(&P);
    phgDofFree(&N);

    phgFreeGrid(&g);
    phgFinalize();

    return 0;
}
