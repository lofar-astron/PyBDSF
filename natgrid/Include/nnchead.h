#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef __APPLE__
#include <malloc.h>
#endif

#define SQ(x)   (x) * (x)
#define BIGNUM  1E37
#define EPSILON 0.00001
#define RANGE   10
#define EQ      ==
#define NE      !=
#define AND     &&
#define OR      ||

extern  double  **points, **joints, wbit,
                horilap, vertlap, bI, bJ, nuldat,
                xstart, ystart, xend, yend,
                maxhoriz, aaa, bbb, ccc, det,
                work3[3][3], xx, sumx, sumy, sumz,
                sumx2, sumy2, sumxy, sumxz, sumyz,
                asum, nn_pi, piby2, piby32, rad2deg,
                bigtri[3][3], horilap_save, vertlap_save;

extern  double  magx, magy, magz, magx_orig, magy_orig, magz_orig,
                maxxy[2][3], magx_auto, magy_auto, magz_auto;

extern  int     igrad, non_neg, densi, sdip, rads, southhemi,
                extrap, adf, nndup;

extern  int     datcnt, datcnt3, numtri, imag, numnei, iscale,
                ext, *jndx, neicnt, optim, goodflag, updir,
                scor[3][2], auto_scale,
                single_point, first_single, asflag,
                error_status;

extern  char    tri_file[256], error_file[256], emsg[256];

extern  FILE    *filee;

extern void     Terminate(void);
extern void     ErrorHnd(int, char *, FILE *, char *);

void            FindNeigh(int);
void            TriNeigh(void);
void            Gradient(void);
void            FindProp(double, double);
double          Surface(void);
double          Meld(double, double, double);
void            TooSteep(void);
void            TooShallow(void);
void            TooNarrow(void);
struct datum    *IMakeDatum(void);
struct simp     *IMakeSimp(void);
struct temp     *IMakeTemp(void);
struct neig     *IMakeNeig(void);
int             *IntVect(int ncols);
void            FreeVecti(int *vectptr);
double          *DoubleVect(int ncols);
void            FreeVectd(double *vectptr);
int             **IntMatrix(int nrows, int ncols);
void            FreeMatrixi(int **matptr);
float           **FloatMatrix(int nrows, int ncols);
void            FreeMatrixf(float **matptr);
double          **DoubleMatrix(int nrows, int ncols);
void            FreeMatrixd(double **matptr);
