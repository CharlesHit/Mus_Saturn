/*            PURPOSE : Matrix memory allocation and mathematical operations

		PREREQUISITES : None

*/
#pragma once

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
	double **m;
	int l, c;
} dmatrix_t;

static dmatrix_t C;

void error(const char error_text[])
{
	void exit(int code);

	fprintf(stderr, "Run-time error...\n");
	fprintf(stderr, "%s\n", error_text);
	exit(1);
}

void write_dmatrix(dmatrix_t *M)

{
	int i, j;

	//printf("row:%d; column:%d\n",M->l,M->c);

	for (i = 1; i <= (*M).l; i++) {
		for (j = 1; j <= (*M).c; j++) {
			printf("%7.4f ", (*M).m[i][j]);
		}
		printf("\n");
	}
}


double **dmatrix(int nrl, int nrh, int ncl, int nch)

{
	int i;
	double **m;

	m = (double **)malloc((unsigned)(nrh - nrl + 1) * sizeof(double));
	if (!m) {
		error("MATRIX.H: allocation failure");
	}
	m -= nrl;

	for (i = nrl; i <= nrh; i++) {
		m[i] = (double *)malloc((unsigned)(nch - ncl + 1) * sizeof(double));
		if (!m[i]) {
			error("MATRIX.H: allocation failure");
		}
		m[i] -= ncl;
	}
	return m;
}

void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)

{ int i ;

    for (i = nrh ; i >= nrl ; i--) {
        free((char *) (m[i] + ncl)) ;
    }
    free((char *) (m + nrl)) ;
}


void delete_dmatrix(dmatrix_t *A)

{ free_dmatrix(A->m,1,A->l,1,A->c) ;
    //free(A) ;
}

void dmat_alloc(dmatrix_t *A, int l, int c)

{
	(*A).m = dmatrix(1, l, 1, c);
	(*A).l = l;
	(*A).c = c;
}


dmatrix_t *dmat_duplicate(dmatrix_t *A)

{
	int i, j;
	dmatrix_t *B;
	B = (dmatrix_t *)malloc(sizeof(dmatrix_t));
	dmat_alloc(B, (*A).l, (*A).c);

	for (i = 1; i <= (*A).l; i++) {
		for (j = 1; j <= (*A).c; j++) {
			(*B).m[i][j] = (*A).m[i][j];
		}
	}
	return B;
}


dmatrix_t *dmat_init(dmatrix_t *A, double a)

{
	int i, j;

	for (i = 1; i <= (*A).l; i++) {
		for (j = 1; j <= (*A).c; j++) {
			(*A).m[i][j] = a;
		}
	}
	return A;
}

/* to make a matrix like 1 2 3 4; 1 2 3 4; 1 2 3 4 */
dmatrix_t *dmat_init_2(dmatrix_t *A)

{
	int i, j;

	for (i = 1; i <= (*A).l; i++) {
		for (j = 1; j <= (*A).c; j++) {
			(*A).m[i][j] = j;
		}
	}
	return A;
}

dmatrix_t *dmat_identity(dmatrix_t *A)

{
	int i, j;

	for (i = 1; i <= (*A).l; i++) {
		for (j = 1; j <= (*A).c; j++) {
			if (i == j) {
				(*A).m[i][j] = 1.0;
			}
			else {
				(*A).m[i][j] = 0.0;
			}
		}
	}
	return A;
}


dmatrix_t *dmat_scalar_mult(dmatrix_t *A, double a)

{
	dmatrix_t *B;
	int i, j;

	B = (dmatrix_t *)malloc(sizeof(dmatrix_t));
	dmat_alloc(B, (*A).l, (*A).c);

	for (i = 1; i <= (*A).l; i++) {
		for (j = 1; j <= (*A).c; j++) {
			(*B).m[i][j] = (*A).m[i][j] * a;
		}
	}
	return B;
}


dmatrix_t *dmat_mult(dmatrix_t *A, dmatrix_t *B)

{
	dmatrix_t *C;
	double s;
	int i, j, k;

	if ((*A).c != (*B).l) {
		error("MATRIX.H: incompatible matrix sizes");
	}
	C = (dmatrix_t *)malloc(sizeof(dmatrix_t));
	dmat_alloc(C, (*A).l, (*B).c);

	for (i = 1; i <= (*C).l; i++) {
		for (j = 1; j <= (*C).c; j++) {
			for (s = 0.0, k = 1; k <= (*A).c; k++) {
				s += (*A).m[i][k] * (*B).m[k][j];
			}
			(*C).m[i][j] = s;
		}
	}
	return C;
}


dmatrix_t *dmat_add(dmatrix_t *A, dmatrix_t *B)

{
	dmatrix_t *C;
	int i, j;

	if ((*A).l != (*B).l || (*A).c != (*B).c) {
		error("MATRIX.H: incompatible matrix sizes");
	}
	C = (dmatrix_t *)malloc(sizeof(dmatrix_t));
	dmat_alloc(C, (*A).l, (*A).c);

	for (i = 1; i <= (*C).l; i++) {
		for (j = 1; j <= (*C).c; j++) {
			(*C).m[i][j] = (*A).m[i][j] + (*B).m[i][j];
		}
	}
	return C;
}


dmatrix_t *dmat_sub(dmatrix_t *A, dmatrix_t *B)

{
	dmatrix_t *C;
	//double s ;
	int i, j;

	if ((*A).l != (*B).l || (*A).c != (*B).c) {
		error("MATRIX.H: incompatible matrix sizes");
	}
	C = (dmatrix_t *)malloc(sizeof(dmatrix_t));
	dmat_alloc(C, (*A).l, (*A).c);

	for (i = 1; i <= (*C).l; i++) {
		for (j = 1; j <= (*C).c; j++) {
			(*C).m[i][j] = (*A).m[i][j] - (*B).m[i][j];
		}
	}
	return C;
}

bool dmat_compare(dmatrix_t *A, dmatrix_t *B)

{
    //double s ;
    int i, j;

    if ((*A).l != (*B).l || (*A).c != (*B).c) {
        error("MATRIX.H_compare: incompatible matrix sizes");
    }

    for (i = 1; i <= (*A).l; i++) {
        for (j = 1; j <= (*A).c; j++) {
            if((*A).m[i][j] - (*B).m[i][j] < 0.00001) return false;
        }
    }
    return true;
}


double dmat_norm(dmatrix_t *A)

{
	double s;
	int i;

	if ((*A).l != 1 && (*A).c != 1) {
		s = 0;
		error("MATRIX.H: incompatible matrix sizes");
	}
	else {
		if ((*A).l == 1) {
			for (s = 0.0, i = 1; i <= (*A).c; i++) {
				s += pow((*A).m[1][i], 2.0);
			}
		}
		else {
			for (s = 0.0, i = 1; i <= (*A).l; i++) {
				s += pow((*A).m[i][1], 2.0);
			}
		}
	}
	return sqrt(s);
}


dmatrix_t *dmat_normalize(dmatrix_t *A)

{
	return dmat_scalar_mult(A, 1.0 / dmat_norm(A));
}


dmatrix_t *dmat_transpose(dmatrix_t *A)

{
	dmatrix_t *B;
	int i, j;

	B = (dmatrix_t *)malloc(sizeof(dmatrix_t));
	dmat_alloc(B, (*A).c, (*A).l);

	for (i = 1; i <= (*A).l; i++) {
		for (j = 1; j <= (*A).c; j++) {
			(*B).m[j][i] = (*A).m[i][j];
		}
	}
	return B;
}


dmatrix_t *dmat_abs(dmatrix_t *A)

{
    int i, j;

    for (i = 1; i <= (*A).l; i++) {
        for (j = 1; j <= (*A).c; j++) {
            if( (*A).m[i][j] < 0)(*A).m[i][j] = - (*A).m[i][j];
        }
    }
    return A;
}


double ddot_product(dmatrix_t *A, dmatrix_t *B)

{
	dmatrix_t *C;

	C = (dmatrix_t *)malloc(sizeof(dmatrix_t));

	if ((*A).c == (*B).c && (*A).l == 1 && (*B).l == 1) {
		C = dmat_mult(A, dmat_transpose(B));
	}
	else if ((*A).c == (*B).l && (*A).l == 1 && (*B).c == 1) {
		C = dmat_mult(A, B);
	}
	else if ((*A).l == (*B).c && (*A).c == 1 && (*B).l == 1) {
		C = dmat_mult(B, A);
	}
	else if ((*A).l == (*B).l && (*A).c == 1 && (*B).c == 1) {
		C = dmat_mult(dmat_transpose(A), B);
	}
	else error("MATRIX.H.dot: Incompatible matrix sizes");
	return (*C).m[1][1];
}


dmatrix_t *cross_product_matrix(dmatrix_t *A, dmatrix_t *B)

{
	int i, j;
	dmatrix_t *C;

	C = (dmatrix_t *)malloc(sizeof(dmatrix_t));
	dmat_alloc(C, (*A).l, (*A).l);

	for (i = 1; i <= (*A).l; i++) {
		for (j = 1; j <= (*A).l; j++) {
			if (i == 1) {
				(*C).m[i][j] = 1.0;
			}
			else if (i == 2) {
				(*C).m[i][j] = (*A).m[j][1];
			}
			else {
				(*C).m[i][j] = (*B).m[j][1];
			}
		}
	}
	return C;
}

dmatrix_t *sub_matrix(dmatrix_t *A, int r, int c)

{
	int i, j, k, l;
	dmatrix_t *B;

	if (r < 1 || r >(*A).l || c < 1 || c >(*A).c || (*A).c < 2 || (*A).l < 2) {
		error("MATRIX.H: erroneous indices");
	}
	B = (dmatrix_t *)malloc(sizeof(dmatrix_t));
	dmat_alloc(B, (*A).l - 1, (*A).c - 1);

	for (i = 1, k = 1; i <= (*A).l; i++) {
		for (j = 1, l = 1; j <= (*A).c; j++) {
			if (j != c && i != r) {
				(*B).m[k][l] = (*A).m[i][j];
			}
			if (j != c) l++;
		}
		if (i != r) k++;
	}
	return B;
}


double determinant(dmatrix_t *A)

{
	int i;
	double det = 0;

	if ((*A).l < 1 || (*A).c < 1) {
		error("MATRIX.H: erroneous matrix size");
	}
	else if ((*A).l != (*A).c) {
		error("MATRIX.H: not a square matrix");
	}
	else if ((*A).l == 1) {
		det = (*A).m[1][1];
	}
	else {
		det = 0.0;
		for (i = 1; i <= (*A).c; i++) {
			det += pow(-1.0, i + 1.0)*(*A).m[1][i] * determinant(sub_matrix(A, 1, i));
		}
	}
	return det;
}


dmatrix_t *cofactor(dmatrix_t *A)

{
	int i, j;
	dmatrix_t *B;

	B = (dmatrix_t *)malloc(sizeof(dmatrix_t));
	dmat_alloc(B, (*A).l, (*A).c);

	for (i = 1; i <= (*A).l; i++) {
		for (j = 1; j <= (*A).c; j++) {
			(*B).m[i][j] = pow(-1.0, i + j)*determinant(sub_matrix(A, i, j));
		}
	}
	return B;
}


dmatrix_t *dmat_inverse(dmatrix_t *A)

{
	return dmat_scalar_mult(dmat_transpose(cofactor(A)), 1.0 / determinant(A));
}


dmatrix_t *dcross_product(dmatrix_t *A, dmatrix_t *B)

{
	int i;
	dmatrix_t *C;

	if ((*A).l != 3 || (*A).c != 1 || (*B).l != 3 || (*B).c != 1) {
		error("MATRIX.H: Incompatible matrix sizes");
	}

	C = (dmatrix_t *)malloc(sizeof(dmatrix_t));
	dmat_alloc(C, (*A).l, (*A).c);

	for (i = 1; i <= (*A).l; i++) {
		(*C).m[i][1] = pow(-1.0, i + 1)*determinant(sub_matrix(cross_product_matrix(A, B), 1, i));
	}
	return C;
}


dmatrix_t *to_homogeneous(dmatrix_t *A, double l)

{
	int i, j;
	dmatrix_t *B;

	if ((*A).l <= 0 || (*A).c <= 0) {
		error("MATRIX.H: erroneous matrix size");
	}

	B = (dmatrix_t *)malloc(sizeof(dmatrix_t));

	if ((*A).c == 1) {
		dmat_alloc(B, (*A).l + 1, 1);
		for (i = 1; i < (*B).l; i++) {
			(*B).m[i][1] = (*A).m[i][1];
		}
		(*B).m[(*B).l][1] = l;
	}
	else if ((*A).l == 1) {
		dmat_alloc(B, 1, (*A).c + 1);
		for (i = 1; i < (*B).c; i++) {
			(*B).m[1][i] = (*A).m[1][i];
		}
		(*B).m[1][(*B).c] = l;
	}
	else {
		dmat_alloc(B, (*A).l + 1, (*A).c + 1);
		B = dmat_init(B, 0.0);
		for (i = 1; i < (*B).l; i++) {
			for (j = 1; j < (*B).c; j++) {
				(*B).m[i][j] = (*A).m[i][j];
			}
		}
		(*B).m[(*B).l][(*B).c] = l;
	}
	return B;
}


dmatrix_t *from_homogeneous(dmatrix_t *A)

{
	int i, j;
	dmatrix_t *B;

	if ((*A).l < 1 || (*A).c < 1) {
		error("MATRIX.H: erroneous matrix size");
	}

	B = (dmatrix_t *)malloc(sizeof(dmatrix_t));

	if ((*A).c == 1) {
		dmat_alloc(B, (*A).l - 1, 1);
		for (i = 1; i < (*A).l; i++) {
			(*B).m[i][1] = (*A).m[i][1];
		}
	}
	else if ((*A).l == 1) {
		dmat_alloc(B, 1, (*A).c - 1);
		for (i = 1; i < (*A).c; i++) {
			(*B).m[1][i] = (*A).m[1][i];
		}
	}
	else {
		dmat_alloc(B, (*A).l - 1, (*A).c - 1);
		for (i = 1; i < (*A).l; i++) {
			for (j = 1; j < (*A).c; j++) {
				(*B).m[i][j] = (*A).m[i][j];
			}
		}
	}
	return B;
}

#endif
