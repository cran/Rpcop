#include "ma.h"

class M_b {
private:
  int Dim;
  float *xo;
  float **Mb;
  float **MId;
  float **MInv;

  // vector ops

  float **inv(float **M);
  float *Mxv(float **M, float *v);
  float **MxM(float **M1, float **M2);
  float *mult_esc(float e, float *v);
  float mult_v(float *v1, float *v2);
  float *sum_v(float *v1, float *v2);
  float *dif_v(float *v2, float *v1);
  float *norma_v(float *v);

public:
  // constructors

  M_b(int Dim, float *b); // we reuse the vector b passed by parameter.
  M_b(int Dim, float **n_M, float *n_xo);
  ~M_b();

  M_b *girar(int eix,
             float angle); // apply a rotation to Mb for the given axis/dimension.
                           // rotation angle is what differentiates candidate
                           // matrices from one another.
  M_b *replicar();         // create a copy of M_b.
  void calcular_la_inversa(); // compute the inverse after Mb rotations and before
                              // applying Mb to points. If Mb is optimal this may
                              // be computed more than once, otherwise only once.
  M_a *donar_M_a(M_a *Ma);    // Ma is the matrix lower-dimensional subspaces need
                              // to map points back to original coordinates
  void rebre_xo(float *punt);
  float *
  aplicar(float *punt); // apply Mb to a point and return lower-space coordinates
  float *desaplicar(float *punt); // map back to original coordinates (inverse of aplicar)
  float *donar_bopt();
};
