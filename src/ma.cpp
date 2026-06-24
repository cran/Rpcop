// Contains the original-space dimension, the matrix, and the xo needed to
// transform the points to the original coordinate system

extern "C" {
#include <stdlib.h>
}

#include "ma.h"

M_a::M_a(int d, int p, float **M, float *x) {
  Dim = d;
  profundidad = p;
  Ma = M;
  xa = x;
}

M_a::~M_a() {
  int i;
  for (i = 0; i < Dim; i++)
    delete[] Ma[i];
  delete[] Ma;
  delete[] xa;
}

float *M_a::aplicar_Ma_punt(float *punt) {
  float *p2;
  float *p3;

  p2 = Mxv(Ma, punt - profundidad); // points are represented in Dim+depth form
                                    // with leading depth coordinates set to 0
  p3 = sum_v(p2, xa);

  delete[] p2;
  return p3;
}

float *M_a::aplicar_Ma_vect(float *vect) {
  float *v3;

  v3 = Mxv(Ma, vect - profundidad); // vectors are stored with depth offset

  return v3; // vector in original coordinates
}

M_a *M_a::donar_M_a(float **Mbopt, float *xo) {
  int i, j, n_prof;
  float **n_Ma, **n_Ma2;
  float *n_xa;

  /* new depth */
  n_prof = profundidad + 1;

  /* new Ma */
  n_Ma = new float *[Dim];
  for (i = 0; i < Dim; i++)
    n_Ma[i] = new float[Dim]();

  for (i = 0; i < profundidad; i++)
    n_Ma[i][i] = 1;

  for (i = 0; i < Dim - profundidad; i++)
    for (j = 0; j < Dim - profundidad; j++)
      n_Ma[i + profundidad][j + profundidad] = Mbopt[i][j];

  n_Ma2 = MxM(Ma, n_Ma);

  for (i = 0; i < Dim; i++)
    delete[] n_Ma[i];
  delete[] n_Ma;

  /* new xa */
  n_xa = aplicar_Ma_punt(xo);

  return new M_a(Dim, n_prof, n_Ma2,
                 n_xa); // pass updated subspace depth
}

// private
// vector ops

float *M_a::Mxv(float **M1, float *v) {
  // vxM, we work with row vectors.
  int i, j;
  float sum;
  float *v3 = new float[Dim]();

  if (!M1 || !v) {
    return v3;
  }

  for (i = 0; i < Dim; i++) {
    sum = 0;
    for (j = 0; j < Dim; j++) {
      if (M1[j]) {
        sum += v[j] * M1[j][i];
      }
    }
    v3[i] = sum;
  }
  return v3;
}

float **M_a::MxM(float **M1, float **M2) {
  // vxM, works with row vectors.
  int i, ii, j;
  float sum;
  float **M3;

  M3 = new float *[Dim];

  for (i = 0; i < Dim; i++)
    M3[i] = new float[Dim]();

  if (!M1 || !M2) {
    return M3;
  }

  for (i = 0; i < Dim; i++) {
    for (ii = 0; ii < Dim; ii++) {
      sum = 0;
      for (j = 0; j < Dim; j++) {
        float lhs = M1[i] ? M1[i][j] : 0.0f;
        float rhs = M2[j] ? M2[j][ii] : 0.0f;
        sum += lhs * rhs;
      }
      M3[i][ii] = sum;
    }
  }
  return M3;
}

float *M_a::sum_v(float *v1, float *v2) {
  int i;
  float *v3;

  v3 = new float[Dim];
  if (!v1 || !v2) {
    for (i = 0; i < Dim; i++) {
      float lhs = v1 ? v1[i] : 0.0f;
      float rhs = v2 ? v2[i] : 0.0f;
      v3[i] = lhs + rhs;
    }
    return v3;
  }
  for (i = 0; i < Dim; i++)
    v3[i] = v1[i] + v2[i];
  return v3;
}
