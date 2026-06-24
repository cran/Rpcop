#include "mb.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
extern "C" {
#include <memory.h>
}

M_b::M_b(int d, float *b) {
  int i;
  int j;
  int pivot;
  float *v_dif, *v_acum1, *v_acum2, *mesc;

  Dim = d;
  xo = NULL;
  /* we create and initialize Mb and MId */
  Mb = new float *[Dim];
  MId = new float *[Dim];
  MInv = NULL;

  Mb[0] = b; // this constructor takes ownership of the input vector.
  pivot = 0;
  while (pivot < Dim && !b[pivot]) {
    pivot++;
  }
  if (pivot == Dim) {
    pivot = 0;
  }

  for (i = 1, j = 0; i < Dim; i++, j++) {
    if (j == pivot) {
      j++;
    }
    Mb[i] = new float[Dim]();
    if (j < Dim) {
      Mb[i][j] = 1;
    }
  }

  for (i = 0; i < Dim; i++) {
    MId[i] = new float[Dim]();
  }
  for (i = 0; i < Dim; i++) {
    MId[i][i] = 1;
  }

  /* Gram-Schmidt orthonormalization process */

  v_acum1 = new float[Dim]();
  v_dif = norma_v(Mb[0]);
  delete[] Mb[0];
  Mb[0] = v_dif;
  for (i = 1; i < Dim; i++) {
    delete[] v_acum1;
    v_acum1 = new float[Dim]();
    for (j = 0; j < i; j++) {
      mesc = mult_esc(mult_v(Mb[i], Mb[j]), Mb[j]);
      v_acum2 = sum_v(v_acum1, mesc);
      delete[] mesc;
      delete[] v_acum1;
      v_acum1 = v_acum2;
    }
    v_dif = dif_v(Mb[i], v_acum1);
    v_acum2 = norma_v(v_dif);
    delete[] v_dif;
    delete[] Mb[i];
    Mb[i] = v_acum2;
  }
  delete[] v_acum1;
}

M_b::M_b(int d, float **n_M, float *n_xo) {
  int i;
  Dim = d;
  Mb = n_M;
  xo = n_xo;

  /* we create and initialize MId */
  MId = new float *[Dim];
  MInv = NULL;

  for (i = 0; i < Dim; i++) {
    MId[i] = new float[Dim]();
  }
  for (i = 0; i < Dim; i++) {
    MId[i][i] = 1;
  }
}

M_b::~M_b() {
  int i;

  for (i = 0; i < Dim; i++) {
    delete[] Mb[i];
    delete[] MId[i];
  }
  if (MInv) {
    for (i = 0; i < Dim; i++) {
      delete[] MInv[i];
    }
    delete[] MInv;
  }

  delete[] Mb;
  delete[] MId;

  // xo ownership is managed by espai (curve-direction routines).
}

M_b *M_b::girar(int eix, float angle) {
  // 0 < axis < Dim: rotate b with respect to this axis in positive direction.

  float **Mbaux;

  // insert the rotation into the identity matrix
  MId[0][0] = cos(angle);
  MId[0][eix] = sin(angle);
  MId[eix][eix] = cos(angle);
  MId[eix][0] = -1 * sin(angle);

  Mbaux = MxM(Mb, MId);

  // restore the identity matrix
  MId[0][0] = 1;
  MId[0][eix] = 0;
  MId[eix][eix] = 1;
  MId[eix][0] = 0;

  return new M_b(Dim, Mbaux, xo);
}

M_b *M_b::replicar() {
  float **c_Mb;
  int i;

  c_Mb = new float *[Dim];
  for (i = 0; i < Dim; i++)
    c_Mb[i] = new float[Dim];

  for (i = 0; i < Dim; i++)
    memmove(c_Mb[i], Mb[i], Dim * sizeof(float));
  // for (j=0;j<Dim;j++) c_Mb[i][j] = Mb[i][j];

  return new M_b(Dim, c_Mb, xo);
}

void M_b::calcular_la_inversa() {
  int i;
  float **c_Mb;
  // compute inverse after Mb rotations and before applying Mb to points.
  // If Mb is optimal, this may be recomputed more than once.

  if (MInv) {
    for (i = 0; i < Dim; i++)
      delete[] MInv[i];
    delete[] MInv;
  }

  /* copy Mb */ // inv() modifies its input

  c_Mb = new float *[Dim];
  for (i = 0; i < Dim; i++)
    c_Mb[i] = new float[Dim];

  for (i = 0; i < Dim; i++)
    memmove(c_Mb[i], Mb[i], Dim * sizeof(float));
  // for (j=0;j<Dim;j++) c_Mb[i][j] = Mb[i][j];

  /* calculate inverse */

  MInv = inv(c_Mb);

  /* delete c_Mb */
  for (i = 0; i < Dim; i++)
    delete[] c_Mb[i];
  delete[] c_Mb;
}

float *M_b::aplicar(
    float *punt) { /* apply Mb and return coordinates in the lower-dimensional subspace */
  float *p2;
  float *p3;
  p2 = dif_v(punt, xo);
  p3 = Mxv(MInv, p2); /* apply inverse */

  // First coordinate is signed offset from the plane; remaining coordinates are
  // projected point coordinates in the lower plane.

  delete[] p2;
  return p3;
}

float *M_b::desaplicar(float *punt) { /* inverse operation of apply */
  float *p2;
  float *p3;

  p2 = Mxv(Mb, punt);
  p3 = sum_v(p2, xo);

  // First coordinate is signed offset from the plane; remaining coordinates are
  // projected point coordinates in the lower plane.

  delete[] p2;
  return p3;
}

void M_b::rebre_xo(float *punt) {
  // It will be executed once for advancement_cluster() and every time the pop
  // candidate crosses the previous pop hyperplane.
  // free xo; // no need to free here. Either b is deleted when creating the new
  // xo, or xo is shared with the previous pop matrix.
  xo = punt;
}

/* Ma is the matrix lower-dimensional subspaces need to map points/vectors */
/* back to original coordinates */

M_a *M_b::donar_M_a(M_a *Ma) {
  // Executed at most once if this is the optimal M_b.

  return Ma->donar_M_a(
      Mb, xo); // xo must have size Dim+prof at this point.
}

float *M_b::donar_bopt() { return Mb[0]; }

// PRIVATE //////////////////////////////////////////////////////////////////

float **M_b::inv(float **M) {
  int i, j, ii;
  float **Inv, Mji;

  Inv = new float *[Dim];
  for (i = 0; i < Dim; i++)
    Inv[i] = new float[Dim]();
  for (i = 0; i < Dim; i++)
    Inv[i][i] = 1;

  for (i = 0; i < Dim; i++) {
    for (j = (i + 1) % Dim; j != i; j = (j + 1) % Dim) {
      Mji = M[j][i];
      for (ii = 0; ii < Dim; ii++) {
        Inv[j][ii] = Inv[j][ii] * M[i][i] - Inv[i][ii] * Mji;
        M[j][ii] = M[j][ii] * M[i][i] - M[i][ii] * Mji;
      }
    }
  }

  /* Legacy alternative implementation was kept here historically; removed for clarity. */
  for (j = 0; j < Dim; j++) {
    for (ii = 0; ii < Dim; ii++)
      Inv[j][ii] = Inv[j][ii] / M[j][j];
  }
  return Inv;
}

float *M_b::Mxv(float **M1, float *v) {
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

float **M_b::MxM(float **M1, float **M2) {
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

float *M_b::mult_esc(float e, float *v) {
  int i;
  float *v3;

  v3 = new float[Dim];
  for (i = 0; i < Dim; i++)
    v3[i] = v[i] * e;
  return v3;
}

float M_b::mult_v(float *v1, float *v2) {
  int i;
  float sum = 0.0;
  if (!v1 || !v2) {
    return sum;
  }
  for (i = 0; i < Dim; i++) {
    sum += v1[i] * v2[i];
  }
  return sum;
}

float *M_b::sum_v(float *v1, float *v2) {
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

float *M_b::dif_v(float *v2, float *v1) {
  int i;
  float *v3;

  v3 = new float[Dim];
  if (!v2 || !v1) {
    for (i = 0; i < Dim; i++) {
      float lhs = v2 ? v2[i] : 0.0f;
      float rhs = v1 ? v1[i] : 0.0f;
      v3[i] = lhs - rhs;
    }
    return v3;
  }
  for (i = 0; i < Dim; i++)
    v3[i] = v2[i] - v1[i];
  return v3;
}

float *M_b::norma_v(float *v) {
  // returns the normalized vector
  int i;
  float nrm = 0.0;
  float *v3;

  v3 = new float[Dim]();
  if (!v) {
    return v3;
  }
  for (i = 0; i < Dim; i++)
    nrm += pow(v[i], 2);
  if (nrm == 0.0f) {
    return v3;
  }
  nrm = sqrt(nrm);
  for (i = 0; i < Dim; i++)
    v3[i] = v[i] / nrm;
  return v3;
}
