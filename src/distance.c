#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

// extends the code from stats/src/distance.c
// to handle auto- and cross-distances.
//
// note: the runtime in the symmetric case is
//       not always optimal.
//
// ceeboo 2007, 2008, 2009, 2014, 2016

#define both_non_NA(a,b) (!ISNAN(a) && !ISNAN(b))
#define both_FINITE(a,b) (R_FINITE(a) && R_FINITE(b))

typedef double (* DFUN)(double *, double *, int, int, int);

static double dfp = 1;

static double minkowski(double *x, double *y, int nx, int ny, int nc)
{
    double dev, dist;
    int count, j;

    count = 0;
    dist  = 0;
    for (j = 0; j < nc; j++) {
        if (both_non_NA(*x, *y)) {
            dev = (*x - *y);
            if (!ISNAN(dev)) {
                dist += R_pow(fabs(dev), dfp);
                count++;
            }
        }
        x += nx;
        y += ny;
    }
    if (count == 0) return NA_REAL;
    if (count != nc) dist /= ((double)count/nc);
    return R_pow(dist, 1.0/dfp);
}

static double euclidean(double *x, double *y, int nx, int ny, int nc)
{
    double dev, dist;
    int count, j;

    count = 0;
    dist  = 0;
    for (j = 0; j < nc; j++) {
	if (both_non_NA(*x, *y)) {
	    dev = (*x - *y);
            if (!ISNAN(dev)) {
		dist += dev * dev;
		count++;
            }
        }
        x += nx;
        y += ny;
    }
    if (count == 0) return NA_REAL;
    if (count != nc) dist /= ((double)count/nc);
    return sqrt(dist);
}

static double maximum(double *x, double *y, int nx, int ny, int nc)
{
    double dev, dist;
    int count, j;

    count = 0;
    dist = -DBL_MAX;
    for (j = 0 ; j < nc ; j++) {
        if (both_non_NA(*x, *y)) {
            dev = fabs(*x - *y);
            if (!ISNAN(dev)) {
                if (dev > dist)
                    dist = dev;
                count++;
            }
        }
        x += nx;
        y += ny;
    }
    if (count == 0) return NA_REAL;
    //if (count != nc) dist /= ((double)count/nc);
    return dist;
}

static double manhattan(double *x, double *y, int nx, int ny, int nc)
{
    double dev, dist;
    int count, j;

    count = 0;
    dist  = 0;
    for (j = 0; j < nc; j++) {
        if (both_non_NA(*x, *y)) {
            dev = fabs(*x - *y);
            if (!ISNAN(dev)) {
                dist += dev;
                count++;
            }
        }
        x += nx;
        y += ny;
    }
    if (count == 0) return NA_REAL;
    if (count != nc) dist /= ((double)count/nc);
    return dist;
}

// fixme: NA for two all-zero vectors

static double canberra(double *x, double *y, int nx, int ny, int nc)
{
    double dev, dist, sum, diff;
    int count, j;

    count = 0;
    dist  = 0;
    for (j = 0 ;j < nc; j++) {
        if (both_non_NA(*x, *y)) {
            sum = fabs(*x + *y);
            diff = fabs(*x - *y);
            if (sum > DBL_MIN || diff > DBL_MIN) {
                dev = diff/sum;
                if (!ISNAN(dev) ||
                   (!R_FINITE(diff) && diff == sum &&
                    /* use Inf = lim x -> oo */ (dev = 1.))) {
                    dist += dev;
                    count++;
                }
            }
        }
        x += nx;
        y += ny;
    }
    if (count == 0) return NA_REAL;
    if (count != nc) dist /= ((double)count/nc);
    return dist;
}

// FIXME why treat not both finite as NA?
static double binary(double *x, double *y, int nx, int ny, int nc)
{
    int total, count, dist;
    int j;

    total = count = dist = 0;
    for (j = 0; j < nc; j++) {
	if (both_non_NA(*x, *y)) {
	    if (*x || *y) {
		count++;
		if (!(*x && *y))
		    dist++;
            }
            total++;
        }
        x += nx;
        y += ny;
    }
    if (total == 0) return NA_REAL;
    if (count == 0) return 0;
    return (double) dist / count;
}

static double matching(double *x, double *y, int nx, int ny, int nc)
{
    int total, count;
    int j;

    total = count = 0;
    for (j = 0; j < nc; j++) {
	if (both_non_NA(*x, *y)) {
	    if (*x != *y) 
		count++;
            total++;
        }
        x += nx;
        y += ny;
    }
    if (total == 0) return NA_REAL;
    if (count == 0) return 0;
    return (double) count / total;
}

static double fuzzy(double *x, double *y, int nx, int ny, int nc)
{
    double dist, smax, smin;
    int count, j;

    count = 0;
    smax = smin = 0;
    for (j = 0; j < nc; j++) {
	if (both_FINITE(*x, *y)) {
	    if (*x > *y) {
		smax += *x;
		smin += *y;
	    }
	    else {
		smax += *y;
		smin += *x;
	    }
	    count++;
        }
        x += nx;
        y += ny;
    }
    if (count == 0) return NA_REAL;
    if (!R_FINITE(smin)) return NA_REAL;
    dist = smin / smax;
    if (ISNAN(dist)) return 0;
    return 1-dist;
}

static double mutual(double *x, double *y, int nx, int ny, int nc)
{
    double dist;
    int total, count, cx, cy, j;

    total = count = cx = cy = 0;
    for (j = 0; j < nc; j++) {
	if (both_non_NA(*x, *y)) {
	    if (*x && *y) 
		count++;
	    cx += (*x && 1);
	    cy += (*y && 1);
            total++;
        }
        x += nx;
        y += ny;
    }
    if (total == 0) return NA_REAL;
    if (cx == 0 || cy == 0 || cx == total || cy == total) return 0;
    dist = 0;
    if (count > 0)
	dist += (double) count / total * log((double) count / cx / cy * total);
    cy = total - cy;
    count = cx - count;
    if (count > 0) 
	dist += (double) count / total * log((double) count / cx / cy * total);
    cx = total - cx;
    count = cy - count;
    if (count > 0)
	dist += (double) count / total * log((double) count / cx / cy * total);
    cy = total - cy;
    count = cx - count;
    if (count > 0)
	dist += (double) count / total * log((double) count / cx / cy * total); 
    if (total != nc) dist /= ((double)total/nc);
    return dist ;
}


// wrapper

static SEXP dists(SEXP R_x, SEXP R_y, SEXP R_d, DFUN f, SEXP R_p) {
    if (!isMatrix(R_x))
	error("'x' not of class matrix");
    if (!isNull(R_y) && !isMatrix(R_y))
	error("'y' not of class matrix");
    if (TYPEOF(R_d) != LGLSXP)
	error("'d' not of type logical");
    int i, j, n, nx, ny, nc, m = 0;
    SEXP x = R_x, y = R_y, r;

    if (!isNull(R_p))			    // fixme: check?
	dfp = *REAL(R_p);
    
    if (isNull(y)) 
	y = x;
    else
    if (LOGICAL(R_d)[0] == TRUE)
	m = 2;
    else
	m = 1;				    // return matrix
    
    nc = INTEGER(GET_DIM(x))[1];
    
    if (INTEGER(GET_DIM(y))[1] != nc)
	error("invalid number of columns");

    nx = INTEGER(GET_DIM(x))[0];
    ny = INTEGER(GET_DIM(y))[0];
   
    if (m == 2 && nx != ny)
	error("invalid number of rows for pairwise mode");

    if (TYPEOF(x) != REALSXP) {
	PROTECT(x = coerceVector(R_x, REALSXP));

	if (isNull(R_y) || R_x == R_y)
	    y = x;
    }

    if (TYPEOF(y) != REALSXP) 
	PROTECT(y = coerceVector(R_y, REALSXP));

    if (m == 0) {
	SEXP d;
	
	PROTECT(r = allocVector(REALSXP, nx*(nx-1)/2));
	setAttrib(r, install("Size"), ScalarInteger(nx));
	
	if (!isNull(d = getAttrib(x, R_DimNamesSymbol)))
	    setAttrib(r, install("Labels"), VECTOR_ELT(d, 0));
	// fixme: package?
	setAttrib(r, R_ClassSymbol, mkString("dist"));
    } else
    if (m == 1) {
	SEXP d1, d2;
	
	PROTECT(r = allocMatrix(REALSXP, nx, ny));

	d1 = getAttrib(x, R_DimNamesSymbol);
	d2 = getAttrib(y, R_DimNamesSymbol);
	if (!isNull(d1) || !isNull(d2)) {
	    SEXP d;
	    
	    setAttrib(r, R_DimNamesSymbol, (d = allocVector(VECSXP, 2)));
	    SET_VECTOR_ELT(d, 0, isNull(d1) ? d1 : VECTOR_ELT(d1, 0));
	    SET_VECTOR_ELT(d, 1, isNull(d2) ? d2 : VECTOR_ELT(d2, 0));
	}
    } else
	PROTECT(r = allocVector(REALSXP, nx));
   
    n = 0;
    for (j = 0; j < ny; j++) {
	if (m == 2)
	    REAL(r)[n++] = f(REAL(x)+j, REAL(y)+j, nx, ny, nc);
	else
	    for (i = (m == 0) ? j+1 : 0; i < nx; i++)
		REAL(r)[n++] = f(REAL(x)+i, REAL(y)+j, nx, ny, nc);
	R_CheckUserInterrupt();
    }

    UNPROTECT(1);
    if (x != R_x)
	UNPROTECT(1);
    if (!isNull(R_y) && R_y != R_x && y != R_y)
	UNPROTECT(1);

    return r;
}

// R wrappers

SEXP R_minkowski_dist(SEXP x, SEXP y, SEXP d, SEXP p) {
    if (isNull(p))
	error("'p' invalid");
    return dists(x, y, d, minkowski, p);
}

SEXP R_euclidean_dist(SEXP x, SEXP y, SEXP d) {
    return dists(x, y, d, euclidean, R_NilValue);
}

SEXP R_maximum_dist(SEXP x, SEXP y, SEXP d) {
    return dists(x, y, d, maximum, R_NilValue);
}

SEXP R_manhattan_dist(SEXP x, SEXP y, SEXP d) {
    return dists(x, y, d, manhattan, R_NilValue);
}

SEXP R_canberra_dist(SEXP x, SEXP y, SEXP d) {
    return dists(x, y, d, canberra, R_NilValue);
}

SEXP R_binary_dist(SEXP x, SEXP y, SEXP d) {
    return dists(x, y, d, binary, R_NilValue);
}

SEXP R_matching_dist(SEXP x, SEXP y, SEXP d) {
    return dists(x, y, d, matching, R_NilValue);
}

SEXP R_fuzzy_dist(SEXP x, SEXP y, SEXP d) {
    return dists(x, y, d, fuzzy, R_NilValue);
}

SEXP R_mutual_dist(SEXP x, SEXP y, SEXP d) {
    return dists(x, y, d, mutual, R_NilValue);
}

// R optimized

// compute jaccard similarities given logical data.
// for implicit binarization use one minus binary
// distances. 

SEXP R_bjaccard(SEXP R_x, SEXP R_y, SEXP R_d) {
    if (!isMatrix(R_x) || TYPEOF(R_x) != LGLSXP)
	error("'x' invalid object or mode");
    if (!isNull(R_y) && (!isMatrix(R_y) || TYPEOF(R_y) != LGLSXP))
	error("'y' invalid object or mode");
    if (TYPEOF(R_d) != LGLSXP)
	error("'d' invalid mode");
    int nc, nx, ny, nz; 
    int i, j, k, l, t, n, m = 0;		// matrix flag
    int *x, *y, *s;
    double z;
    SEXP r;
    
    if (isNull(R_y))
	R_y = R_x;
    else
    if (LOGICAL(R_d)[0] == TRUE)
	m = 2;
    else
	m = 1;
    
    nc = INTEGER(GET_DIM(R_x))[1];
    
    if (INTEGER(GET_DIM(R_y))[1] != nc)
	error("tha number of columns of 'x' and 'y' do not conform");

    nz =
    nx = INTEGER(GET_DIM(R_x))[0];
    ny = INTEGER(GET_DIM(R_y))[0];

    if (m == 0) {
	SEXP d;
	
	PROTECT(r = allocVector(REALSXP, nx*(nx-1)/2));
	setAttrib(r, install("Size"), ScalarInteger(nx));
	
	if (!isNull(d = getAttrib(R_x, R_DimNamesSymbol)))
	    setAttrib(r, install("Labels"), VECTOR_ELT(d, 0));
	// fixme: package?
	setAttrib(r, R_ClassSymbol, mkString("dist"));
    } else
    if (m == 1) {
	SEXP d1, d2;
	
	PROTECT(r = allocMatrix(REALSXP, nx, ny));

	d1 = getAttrib(R_x, R_DimNamesSymbol);
	d2 = getAttrib(R_y, R_DimNamesSymbol);
	if (!isNull(d1) || !isNull(d2)) {
	    SEXP d;
	    
	    setAttrib(r, R_DimNamesSymbol, (d = allocVector(VECSXP, 2)));
	    SET_VECTOR_ELT(d, 0, isNull(d1) ? d1 : VECTOR_ELT(d1, 0));
	    SET_VECTOR_ELT(d, 1, isNull(d2) ? d2 : VECTOR_ELT(d2, 0));
	}
    } else {
	if (nx != ny)
	    error("the number of rows of 'x' and 'y' do not conform");

	PROTECT(r = allocVector(REALSXP, nx));
    }

    x = INTEGER(R_x);
    y = INTEGER(R_y);
    
    s = INTEGER(PROTECT(allocVector(INTSXP, nx)));
    memset(s, 0, sizeof(int)*nx);
    
    for (i = 0; i < nx; i++) {
	t = 0;
	for (k = 0; k < nc; k++) {
	    if (x[i+k*nx] == NA_LOGICAL)
		continue;
	    t += x[i+k*nx] == TRUE;
	}
	s[i] = t;	
    }

    n = 0;
    for (j = 0; j < ny; j++) {
	if (m == 0) {
	    t = s[j];
	    i = j + 1;
	}
	else {
	    t = 0;
	    for (k = 0; k < nc; k++) {
		if (y[j+k*ny] == NA_LOGICAL)
		    continue;
		t += y[j+k*ny] == TRUE;
	    }
	    if (m == 1) 
		i  = 0;
	    else {
		i  = j;
		nz = j + 1;
	    }
	}
	for (; i < nz; i++) {
	    l = 0;
	    for (k = 0; k < nc; k++) {
		if (x[i+k*nx] == NA_LOGICAL || 
		    y[j+k*ny] == NA_LOGICAL)
		    continue;
		l += (x[i+k*nx] == TRUE) & 
		     (y[j+k*ny] == TRUE);
	    }
	    z = (double) l / (t + s[i] - l);
	    if (ISNAN(z))			    /* division by zero */
		REAL(r)[n++] = 1;		    /* but be compatible */
	    else
		REAL(r)[n++] = z;
	}
	R_CheckUserInterrupt();
    }
   
    UNPROTECT(2);

    return r;
}

/* calculate Jaccard similarities extended to real-
 * valued data, i.e the scalar product divided by the
 * squared Euclidean distance plus the scalar product.
 */

SEXP R_ejaccard(SEXP R_x, SEXP R_y, SEXP R_d) {
    if (!isMatrix(R_x))
	error("'x' not of class matrix");
    if (!isNull(R_y) && !isMatrix(R_x))
	error("'y' not of class matrix");
    if (TYPEOF(R_d) != LGLSXP)
	error("'d' not of type logical");
    int nc, nx, ny, nz; 
    int i, j, k, l, n, m = 0;
    double t, z;
    double *x, *y, *s;
    SEXP r_x = R_x, r_y = R_y, r;
    
    if (isNull(R_y))
	R_y = R_x;
    else
    if (LOGICAL(R_d)[0] == TRUE)
	m = 2;
    else
	m = 1;
		    
    nc = INTEGER(GET_DIM(R_x))[1];

    if (INTEGER(GET_DIM(R_y))[1] != nc)
	error("the number of columns of 'x' and 'y' do not conform");

    nz =
    nx = INTEGER(GET_DIM(R_x))[0];
    ny = INTEGER(GET_DIM(R_y))[0];

    if (m == 2 && nx != ny)
	error("the number f rows of 'x' and 'y' do not conform");

    if (TYPEOF(R_x) != REALSXP) {
	PROTECT(R_x = coerceVector(r_x, REALSXP));

	if (isNull(r_y) || r_x == r_y)
	    R_y = R_x;
    }

    if (TYPEOF(R_y) != REALSXP) 
	PROTECT(R_y = coerceVector(r_y, REALSXP));

    if (m == 0) {
	SEXP d;
	
	PROTECT(r = allocVector(REALSXP, nx*(nx-1)/2));
	setAttrib(r, install("Size"), ScalarInteger(nx));
	
	if (!isNull(d = getAttrib(R_x, R_DimNamesSymbol)))
	    setAttrib(r, install("Labels"), VECTOR_ELT(d, 0));
	// fixme: package?
	setAttrib(r, R_ClassSymbol, mkString("dist"));
    } else
    if (m == 1) {
	SEXP d1, d2;
	
	PROTECT(r = allocMatrix(REALSXP, nx, ny));

	d1 = getAttrib(R_x, R_DimNamesSymbol);
	d2 = getAttrib(R_y, R_DimNamesSymbol);
	if (!isNull(d1) || !isNull(d2)) {
	    SEXP d;
	    
	    setAttrib(r, R_DimNamesSymbol, (d = allocVector(VECSXP, 2)));
	    SET_VECTOR_ELT(d, 0, isNull(d1) ? d1 : VECTOR_ELT(d1, 0));
	    SET_VECTOR_ELT(d, 1, isNull(d2) ? d2 : VECTOR_ELT(d2, 0));
	}
    } else 
	PROTECT(r = allocVector(REALSXP, nx));

    x = REAL(R_x);
    y = REAL(R_y);
	
    s = REAL(PROTECT(allocVector(REALSXP, nx)));
    memset(s, 0, sizeof(double)*nx);
    
    for (i = 0; i < nx; i++) {
	z = 0;
	l = 0;
	for (k = 0; k < nc; k++) {
	    if (!R_FINITE(x[i+k*nx]))
		continue;
	    l++;
	    z+= pow(x[i+k*nx], 2);
	}
	s[i] = (l > 0) ? z : NA_REAL;
    }

    n = 0; 
    for (j = 0; j < ny; j++) {
	if (m == 0) {
	    t = s[j];
	    i = j + 1;
	}
	else {
	    z = 0;
	    l = 0;
	    for (k = 0; k < nc; k++) {
		if (!R_FINITE(y[j+k*ny]))
		    continue;
		l++;
		z+= pow(y[j+k*ny], 2);
	    }
	    t = (l > 0) ? z : NA_REAL;
	    if (m == 1)
		i = 0;
	    else {
		i  = j;
		nz = j + 1;
	    }
	}
	for (; i < nz; i++) {
	    if (!R_FINITE(t) || !R_FINITE(s[i])) {
		REAL(r)[n++] = NA_REAL;
		continue;
	    }
	    l = 0;
	    z = 0;
	    for (k = 0; k < nc; k++) {
	        if (!R_FINITE(x[i+k*nx]) || !R_FINITE(y[j+k*ny]))
		    continue;
		l++;
		z+= x[i+k*nx] * y[j+k*ny];
	    }
	    if (l > 0) {
	        z = z / (t + s[i] - z);
	        if (ISNAN(z))
		    REAL(r)[n++] = 1;	    /* be compatible */
		else
		    REAL(r)[n++] = z;
	    }
	    else
		REAL(r)[n++] = NA_REAL;
	}
	R_CheckUserInterrupt();
    }

    UNPROTECT(2);
    if (R_x != r_x)
	UNPROTECT(1);
    if (!isNull(r_y) && r_y != r_x && R_y != r_y)
	UNPROTECT(1);

    return r;
}

/* calculate Dice similarities extended to real-
 * valued data, i.e. twice the scalar product divided by the
 * squared Euclidean distance plus twice the scalar product.
 */

SEXP R_edice(SEXP R_x, SEXP R_y, SEXP R_d) {
    if (!isMatrix(R_x))
	error("'x' not of class matrix");
    if (!isNull(R_y) && !isMatrix(R_x))
	error("'y' not of class matrix");
    if (TYPEOF(R_d) != LGLSXP)
	error("'d' not of type logical");
    int nc, nx, ny, nz; 
    int i, j, k, l, n, m = 0;
    double t, z;
    double *x, *y, *s;
    SEXP r_x = R_x, r_y = R_y, r;
    
    if (isNull(R_y))
	R_y = R_x;
    else
    if (LOGICAL(R_d)[0] == TRUE)
	m = 2;
    else
	m = 1;
		    
    nc = INTEGER(GET_DIM(R_x))[1];

    if (INTEGER(GET_DIM(R_y))[1] != nc)
	error("the number of columns of 'x' and 'y' do not conform");

    nz =
    nx = INTEGER(GET_DIM(R_x))[0];
    ny = INTEGER(GET_DIM(R_y))[0];

    if (m == 2 && nx != ny)
	error("the number f rows of 'x' and 'y' do not conform");

    if (TYPEOF(R_x) != REALSXP) {
	PROTECT(R_x = coerceVector(r_x, REALSXP));

	if (isNull(r_y) || r_x == r_y)
	    R_y = R_x;
    }

    if (TYPEOF(R_y) != REALSXP) 
	PROTECT(R_y = coerceVector(r_y, REALSXP));

    if (m == 0) {
	SEXP d;
	
	PROTECT(r = allocVector(REALSXP, nx*(nx-1)/2));
	setAttrib(r, install("Size"), ScalarInteger(nx));
	
	if (!isNull(d = getAttrib(R_x, R_DimNamesSymbol)))
	    setAttrib(r, install("Labels"), VECTOR_ELT(d, 0));
	// fixme: package?
	setAttrib(r, R_ClassSymbol, mkString("dist"));
    } else
    if (m == 1) {
	SEXP d1, d2;
	
	PROTECT(r = allocMatrix(REALSXP, nx, ny));

	d1 = getAttrib(R_x, R_DimNamesSymbol);
	d2 = getAttrib(R_y, R_DimNamesSymbol);
	if (!isNull(d1) || !isNull(d2)) {
	    SEXP d;
	    
	    setAttrib(r, R_DimNamesSymbol, (d = allocVector(VECSXP, 2)));
	    SET_VECTOR_ELT(d, 0, isNull(d1) ? d1 : VECTOR_ELT(d1, 0));
	    SET_VECTOR_ELT(d, 1, isNull(d2) ? d2 : VECTOR_ELT(d2, 0));
	}
    } else 
	PROTECT(r = allocVector(REALSXP, nx));

    x = REAL(R_x);
    y = REAL(R_y);
	
    s = REAL(PROTECT(allocVector(REALSXP, nx)));
    memset(s, 0, sizeof(double)*nx);
    
    for (i = 0; i < nx; i++) {
	z = 0;
	l = 0;
	for (k = 0; k < nc; k++) {
	    if (!R_FINITE(x[i+k*nx]))
		continue;
	    l++;
	    z+= pow(x[i+k*nx], 2);
	}
	s[i] = (l > 0) ? z : NA_REAL;
    }

    n = 0; 
    for (j = 0; j < ny; j++) {
	if (m == 0) {
	    t = s[j];
	    i = j + 1;
	}
	else {
	    z = 0;
	    l = 0;
	    for (k = 0; k < nc; k++) {
		if (!R_FINITE(y[j+k*ny]))
		    continue;
		l++;
		z+= pow(y[j+k*ny], 2);
	    }
	    t = (l > 0) ? z : NA_REAL;
	    if (m == 1)
		i = 0;
	    else {
		i  = j;
		nz = j + 1;
	    }
	}
	for (; i < nz; i++) {
	    if (!R_FINITE(t) || !R_FINITE(s[i])) {
		REAL(r)[n++] = NA_REAL;
		continue;
	    }
	    l = 0;
	    z = 0;
	    for (k = 0; k < nc; k++) {
	        if (!R_FINITE(x[i+k*nx]) || !R_FINITE(y[j+k*ny]))
		    continue;
		l++;
		z+= x[i+k*nx] * y[j+k*ny];
	    }
	    if (l > 0) {
	        z = 2 * z / (t + s[i]);
	        if (ISNAN(z))
		    REAL(r)[n++] = 1;	    /* be compatible */
		else
		    REAL(r)[n++] = z;
	    }
	    else
		REAL(r)[n++] = NA_REAL;
	}
	R_CheckUserInterrupt();
    }

    UNPROTECT(2);
    if (R_x != r_x)
	UNPROTECT(1);
    if (!isNull(r_y) && r_y != r_x && R_y != r_y)
	UNPROTECT(1);

    return r;
}

/* calculate cosine similarities.
 */

SEXP R_cosine(SEXP R_x, SEXP R_y, SEXP R_d) {
    if (!isMatrix(R_x))
	error("'x' not of class matrix");
    if (!isNull(R_y) && !isMatrix(R_y))
	error("'y' not of class matrix");
    if (TYPEOF(R_d) != LGLSXP)
	error("'d' not of type logical");
    int nc, nx, ny, nz; 
    int i, j, k, l, n, m = 0;
    double t, z;
    double *x, *y, *s;
    SEXP r_x = R_x, r_y = R_y, r;
	        
    if (isNull(R_y))
	R_y = R_x;
    else
    if (LOGICAL(R_d)[0] == TRUE)
	m = 2;
    else
	m = 1;
		    
    nc = INTEGER(GET_DIM(R_x))[1];

    if (INTEGER(GET_DIM(R_y))[1] != nc)
       error("the number of columns of 'x' and 'y' do not conform");

    nz =
    nx = INTEGER(GET_DIM(R_x))[0];
    ny = INTEGER(GET_DIM(R_y))[0];

    if (m == 2 && nx != ny)
	error("the number of rows of 'x' and 'y' do not conform");

    if (TYPEOF(R_x) != REALSXP) {
	PROTECT(R_x = coerceVector(r_x, REALSXP));

	if (isNull(r_y) || r_x == r_y)
	    R_y = R_x;
    }

    if (TYPEOF(R_y) != REALSXP) 
	PROTECT(R_y = coerceVector(r_y, REALSXP));

    if (m == 0) {
	SEXP d;
	
	PROTECT(r = allocVector(REALSXP, nx*(nx-1)/2));
	setAttrib(r, install("Size"), ScalarInteger(nx));
	
	if (!isNull(d = getAttrib(R_x, R_DimNamesSymbol)))
	    setAttrib(r, install("Labels"), VECTOR_ELT(d, 0));
	// fixme: package?
	setAttrib(r, R_ClassSymbol, mkString("dist"));
    } else
    if (m == 1) {
	SEXP d1, d2;
	
	PROTECT(r = allocMatrix(REALSXP, nx, ny));

	d1 = getAttrib(R_x, R_DimNamesSymbol);
	d2 = getAttrib(R_y, R_DimNamesSymbol);
	if (!isNull(d1) || !isNull(d2)) {
	    SEXP d;
	    
	    setAttrib(r, R_DimNamesSymbol, (d = allocVector(VECSXP, 2)));
	    SET_VECTOR_ELT(d, 0, isNull(d1) ? d1 : VECTOR_ELT(d1, 0));
	    SET_VECTOR_ELT(d, 1, isNull(d2) ? d2 : VECTOR_ELT(d2, 0));
	}
	setAttrib(r, R_ClassSymbol, mkString("crossdist"));
    } else 
	PROTECT(r = allocVector(REALSXP, nx));

    x = REAL(R_x);
    y = REAL(R_y);
	
    s = REAL(PROTECT(allocVector(REALSXP, nx)));
    memset(s, 0, sizeof(double)*nx);
    
    for (i = 0; i < nx; i++) {
	z = 0;
	l = 0;
	for (k = 0; k < nc; k++) {
	    if (!R_FINITE(x[i+k*nx]))
		continue;
	    l++;
	    z+= pow(x[i+k*nx], 2);
	}
	s[i] = (l > 0) ? sqrt(z) : NA_REAL;
    }

    n = 0; 
    for (j = 0; j < ny; j++) {
	if (m == 0) {
	    t = s[j];
	    i = j + 1;
	}
	else {
	    z = 0;
	    l = 0;
	    for (k = 0; k < nc; k++) {
		if (!R_FINITE(y[j+k*ny]))
		    continue;
		l++;
		z+= pow(y[j+k*ny], 2);
	    }
	    t = (l > 0) ? sqrt(z) : NA_REAL;
	    if (m == 1)
		i = 0;
	    else {
		i  = j;
		nz = j + 1;
	    }
	}
	for (; i < nz; i++) {
	    if (!R_FINITE(t) || !R_FINITE(s[i])) {
		REAL(r)[n++] = NA_REAL;
		continue;
	    }
	    l = 0;
	    z = 0;
	    for (k = 0; k < nc; k++) {
	        if (!R_FINITE(x[i+k*nx]) || !R_FINITE(y[j+k*ny]))
		    continue;
		l++;
		z+= x[i+k*nx] * y[j+k*ny];
	    }
	    if (l > 0) {
	        z =  z / t / s[i];
	        if (ISNAN(z)) {
		    if (t < DBL_MIN && s[i] < DBL_MIN)
			REAL(r)[n++] = 1;	    /* be compatible */
		    else
			REAL(r)[n++] = 0;
		}
		else
		    REAL(r)[n++] = z;
	    }
	    else
		REAL(r)[n++] = NA_REAL;
	}
	R_CheckUserInterrupt();
    }

    UNPROTECT(2);
    if (R_x != r_x)
	UNPROTECT(1);
    if (!isNull(r_y) && r_x != r_y && R_y != r_y)
	UNPROTECT(1);

    return r;
}

//
