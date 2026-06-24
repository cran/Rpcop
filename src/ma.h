// Contains the original-space dimension, the matrix, and the xo needed to
// transform the points to the original coordinate system

class M_a {
private:
  int Dim;
  int profundidad;
  float **Ma;
  float *xa;

  // vector ops

  float *Mxv(float **M1, float *v);
  float **MxM(float **M1, float **M2);
  float *sum_v(float *v1, float *v2);

public:
  // constructors
  M_a(int d, int p, float **M, float *x);
  ~M_a();

  float *aplicar_Ma_punt(float *punt);
  float *aplicar_Ma_vect(float *vect);
  M_a *donar_M_a(float **Mbopt,
                 float *xo); // builds M_a for the new subspace, mapping points
                             // to original-space coordinates
};
