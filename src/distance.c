#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

// extends the code from stats/src/distance.c
// to handle auto- and cross-distances.
//
// note: the runtime in the symmetric case is
//       not always optimal.
//
// ceeboo 2007

#define both_non_NA(a,b) (!ISNAN(a) && !ISNAN(b))
#define both_FINITE(a,b) (R_FINITE(a) && R_FINITE(b))

typedef double (* DFUN)(double *, double *, int, int, int);

static double dfp = 1;

static double minkowski(double *x, double *y, int nx, int ny, int nc)
{
    if (x == y) return 0;
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
    if (x == y) return 0;
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
    if (x == y) return 0;
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
    if (x == y) return 0;
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
    if (x == y) return 0;
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

static double binary(double *x, double *y, int nx, int ny, int nc)
{
    if (x == y) return 0;
    int total, count, dist;
    int j;

    total = count = dist = 0;
    for (j = 0; j < nc; j++) {
	if (both_FINITE(*x, *y)) {
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
    if (x == y) return 0;
    int total, count;
    int j;

    total = count = 0;
    for (j = 0; j < nc; j++) {
	if (both_FINITE(*x, *y)) {
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
    if (x == y) return 0;
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
    dist = smin / smax;
    if (ISNAN(dist)) return 0;
    return 1-dist;
}

static double mutual(double *x, double *y, int nx, int ny, int nc)
{
    if (x == y) return 0;
    double dist;
    int total, count, cx, cy, j;

    total = count = cx = cy = 0;
    for (j = 0; j < nc; j++) {
	if (both_FINITE(*x, *y)) {
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

static SEXP dists(SEXP R_x, SEXP R_y, DFUN f, SEXP R_p) {
    if (!isMatrix(R_x))
	error("'x' not of class matrix");
    if (!isNull(R_y) && !isMatrix(R_y))
	error("'y' not of class matrix");
    int i, j, n, nx, ny, nc, m = 0;
    SEXP x = R_x, y = R_y, r;

    if (!isNull(R_p))			    // fixme: check?
	dfp = *REAL(R_p);
    
    nc = INTEGER(GET_DIM(x))[1];
    
    if (isNull(y)) 
	y = x;
    else
	m = 1;				    // return matrix
    
    if (INTEGER(GET_DIM(y))[1] != nc)
	error("invalid number of columns");

    if (TYPEOF(x) != REALSXP) {
	PROTECT(x = coerceVector(R_x, REALSXP));

	if (isNull(R_y) || R_x == R_y)
	    y = x;
    }

    if (TYPEOF(y) != REALSXP) 
	PROTECT(y = coerceVector(R_y, REALSXP));

    nx = INTEGER(GET_DIM(x))[0];
    ny = INTEGER(GET_DIM(y))[0];
   
    if (!m && x == y) {
	SEXP d;
	
	PROTECT(r = allocVector(REALSXP, nx*(nx-1)/2));
	setAttrib(r, install("Size"), ScalarInteger(nx));
	
	if (!isNull(d = getAttrib(x, R_DimNamesSymbol)))
	    setAttrib(r, install("Labels"), VECTOR_ELT(d, 0));
	// fixme: package?
	setAttrib(r, R_ClassSymbol, mkString("dist"));
	
    } else {
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
    }
   
    n = 0;
    for (j = 0; j < ny; j++) {
	for (i = (!m && x == y) ? j+1 : 0; i < nx; i++)
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

// R wrapper for testing

SEXP R_dists(SEXP x, SEXP y, SEXP f, SEXP p) {
    static DFUN dfun[] = {
	minkowski,      //  1 
	euclidean,      //  2
	maximum,	//  3
	manhattan,	//  4
	canberra,	//  5
	binary,		//  6
	matching,	//  7
	fuzzy,		//  8 similarity
	mutual		//  9 similarity
    };
    int c;

    c = INTEGER(f)[0];
    if (c < 0 || c > sizeof(dfun) / sizeof(*dfun))
	error("'f' invalid index");
     
    return dists(x, y, dfun[c-1], p);
}

// R wrappers

SEXP R_minkowski_dist(SEXP x, SEXP y, SEXP p) {
    return dists(x, y, minkowski, p);
}

SEXP R_euclidean_dist(SEXP x, SEXP y, SEXP p) {
    return dists(x, y, euclidean, R_NilValue);
}

SEXP R_maximum_dist(SEXP x, SEXP y, SEXP p) {
    return dists(x, y, maximum, R_NilValue);
}

SEXP R_manhattan_dist(SEXP x, SEXP y, SEXP p) {
    return dists(x, y, manhattan, R_NilValue);
}

SEXP R_canberra_dist(SEXP x, SEXP y, SEXP p) {
    return dists(x, y, canberra, R_NilValue);
}

SEXP R_binary_dist(SEXP x, SEXP y, SEXP p) {
    return dists(x, y, binary, R_NilValue);
}

SEXP R_matching_dist(SEXP x, SEXP y, SEXP p) {
    return dists(x, y, matching, R_NilValue);
}

SEXP R_fuzzy_dist(SEXP x, SEXP y, SEXP p) {
    return dists(x, y, fuzzy, R_NilValue);
}

SEXP R_mutual_dist(SEXP x, SEXP y, SEXP p) {
    return dists(x, y, mutual, R_NilValue);
}

// R optimized

// compute jaccard similarities given logical data.
// for implicit binarization use one minus binary
// distances. 

SEXP R_bjaccard(SEXP R_x, SEXP R_y) {
    if (!isMatrix(R_x) || TYPEOF(R_x) != LGLSXP)
	error("'x' invalid object or mode");
    if (!isNull(R_y) && (!isMatrix(R_y) || TYPEOF(R_y) != LGLSXP))
	error("'y' invalid object or mode");
    int nc, nx, ny; 
    int i, j, k, l, t, n, m = 0;		// matrix flag
    int *x, *y, *s;
    double z;
    SEXP r;
    
    nc = INTEGER(GET_DIM(R_x))[1];
    
    if (isNull(R_y))
	R_y = R_x;
    else
	m = 1;
    
    if (INTEGER(GET_DIM(R_y))[1] != nc)
	error("tha number of columns of 'x' and 'y' do not conform");

    nx = INTEGER(GET_DIM(R_x))[0];
    ny = INTEGER(GET_DIM(R_y))[0];

    if (!m && R_x == R_y) {
	SEXP d;
	
	PROTECT(r = allocVector(REALSXP, nx*(nx-1)/2));
	setAttrib(r, install("Size"), ScalarInteger(nx));
	
	if (!isNull(d = getAttrib(R_x, R_DimNamesSymbol)))
	    setAttrib(r, install("Labels"), VECTOR_ELT(d, 0));
	// fixme: package?
	setAttrib(r, R_ClassSymbol, mkString("dist"));
    } else {
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
	if (!m && R_x == R_y) {
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
	    i = 0;
	}
	for (i = i; i < nx; i++) {
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

/* calculate extended Jaccard similarities, i.e the
 * scalar product divided by the squared Euclidean
 * distance minus the scalar product. 
 */

SEXP R_ejaccard(SEXP R_x, SEXP R_y) {
    if (!isMatrix(R_x))
	error("'x' not of class matrix");
    if (!isNull(R_y) && !isMatrix(R_x))
	error("'y' not of class matrix");
    int nc, nx, ny; 
    int i, j, k, l, n, m = 0;
    double t, z;
    double *x, *y, *s;
    SEXP r_x = R_x, r_y = R_y, r;
    
    if (isNull(R_y))
	R_y = R_x;
    else
	m = 1;
		    
    nc = INTEGER(GET_DIM(R_x))[1];

    if (INTEGER(GET_DIM(R_y))[1] != nc)
	error("the number of columns of 'x' and 'y' do not conform");
	
    if (TYPEOF(R_x) != REALSXP) {
	PROTECT(R_x = coerceVector(r_x, REALSXP));

	if (isNull(r_y) || r_x == r_y)
	    R_y = R_x;
    }

    if (TYPEOF(R_y) != REALSXP) 
	PROTECT(R_y = coerceVector(r_y, REALSXP));

    nx = INTEGER(GET_DIM(R_x))[0];
    ny = INTEGER(GET_DIM(R_y))[0];

    if (!m && R_x == R_y) {
	SEXP d;
	
	PROTECT(r = allocVector(REALSXP, nx*(nx-1)/2));
	setAttrib(r, install("Size"), ScalarInteger(nx));
	
	if (!isNull(d = getAttrib(R_x, R_DimNamesSymbol)))
	    setAttrib(r, install("Labels"), VECTOR_ELT(d, 0));
	// fixme: package?
	setAttrib(r, R_ClassSymbol, mkString("dist"));
    } else {
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
    }

    x = REAL(R_x);
    y = REAL(R_y);
	
    s = REAL(PROTECT(allocVector(REALSXP, nx)));
    memset(s, 0, sizeof(double)*nx);
    
    for (i = 0; i < nx; i++) {
	z = 0;
	l = 0;
	for (k = 0; k < nc; k++) {
	    if (ISNAN(x[i+k*nx]))
		continue;
	    l++;
	    z+= pow(x[i+k*nx], 2);
	}
	s[i] = (l > 0) ? z : NA_REAL;
    }

    n = 0; 
    for (j = 0; j < ny; j++) {
	if (!m && R_x == R_y) {
	    t = s[j];
	    i = j + 1;
	}
	else {
	    z = 0;
	    l = 0;
	    for (k = 0; k < nc; k++) {
		if (ISNAN(y[j+k*ny]))
		    continue;
		l++;
		z+= pow(y[j+k*ny], 2);
	    }
	    t = (l > 0) ? z : NA_REAL;
	    i = 0;
	}
	for (i = i; i < nx; i++) {
	    if (ISNAN(t) || ISNAN(s[i])) {
		REAL(r)[n++] = NA_REAL;
		continue;
	    }
	    l = 0;
	    z = 0;
	    for (k = 0; k < nc; k++) {
	        if (ISNAN(x[i+k*nx]) || ISNAN(y[j+k*ny]))
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

/* calculate cosine similarities.
 */

SEXP R_cosine(SEXP R_x, SEXP R_y) {
    if (!isMatrix(R_x))
	error("'x' not of class matrix");
    if (!isNull(R_y) && !isMatrix(R_y))
	error("'y' not of class matrix");
    int nc, nx, ny; 
    int i, j, k, l, n, m = 0;
    double t, z;
    double *x, *y, *s;
    SEXP r_x = R_x, r_y = R_y, r;
	        
    if (isNull(R_y))
	R_y = R_x;
    else
	m = 1;
		    
    nc = INTEGER(GET_DIM(R_x))[1];

    if (INTEGER(GET_DIM(R_y))[1] != nc)
       error("the number of columns of 'x' and 'y' do not conform");
	
    if (TYPEOF(R_x) != REALSXP) {
	PROTECT(R_x = coerceVector(r_x, REALSXP));

	if (isNull(r_y) || r_x == r_y)
	    R_y = R_x;
    }

    if (TYPEOF(R_y) != REALSXP) 
	PROTECT(R_y = coerceVector(r_y, REALSXP));

    nx = INTEGER(GET_DIM(R_x))[0];
    ny = INTEGER(GET_DIM(R_y))[0];

    if (!m && R_x == R_y) {
	SEXP d;
	
	PROTECT(r = allocVector(REALSXP, nx*(nx-1)/2));
	setAttrib(r, install("Size"), ScalarInteger(nx));
	
	if (!isNull(d = getAttrib(R_x, R_DimNamesSymbol)))
	    setAttrib(r, install("Labels"), VECTOR_ELT(d, 0));
	// fixme: package?
	setAttrib(r, R_ClassSymbol, mkString("dist"));
    } else {
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
    }

    x = REAL(R_x);
    y = REAL(R_y);
	
    s = REAL(PROTECT(allocVector(REALSXP, nx)));
    memset(s, 0, sizeof(double)*nx);
    
    for (i = 0; i < nx; i++) {
	z = 0;
	l = 0;
	for (k = 0; k < nc; k++) {
	    if (ISNAN(x[i+k*nx]))
		continue;
	    l++;
	    z+= pow(x[i+k*nx], 2);
	}
	s[i] = (l > 0) ? sqrt(z) : NA_REAL;
    }

    n = 0; 
    for (j = 0; j < ny; j++) {
	if (!m && R_x == R_y) {
	    t = s[j];
	    i = j + 1;
	}
	else {
	    z = 0;
	    l = 0;
	    for (k = 0; k < nc; k++) {
		if (ISNAN(y[j+k*ny]))
		    continue;
		l++;
		z+= pow(y[j+k*ny], 2);
	    }
	    t = (l > 0) ? sqrt(z) : NA_REAL;
	    i = 0;
	}
	for (i = i; i < nx; i++) {
	    if (ISNAN(t) || ISNAN(s[i])) {
		REAL(r)[n++] = NA_REAL;
		continue;
	    }
	    l = 0;
	    z = 0;
	    for (k = 0; k < nc; k++) {
	        if (ISNAN(x[i+k*nx]) || ISNAN(y[j+k*ny]))
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
