
#define FALSE 0
#define TRUE 1
#define ESQUERRA 0
#define DRETA 1
#define X 0
#define INF 9999

extern "C" {
#include <math.h>
#include <memory.h>
#include <stdlib.h>
#include <string.h>
}
#include "ll_q.h"
#include "pila.h"

class ll_p {

private:
  int Dim;
  float dmax;
  int orcluster; // candidate scan starts on the right side of xo, then on the
                 // left side along the X coordinate
  int numcl;     // keeps track of cluster traversal rounds
  float sum_w;   // sum of point weights
  float suma_d;

  typedef struct node {
    float *coord;
    int marca;
    node *seg[2];
    void *noin[2]; // 1:link nodes outside the spanning tree. 2:link
                   // satellites of the node
  } node;

  typedef struct node_satelit {
    node *ptnode;
    node_satelit *seg;
  } node_satelit;

  node *xorig;   // middle point of the curve
  node *xoant;   // cluster point close to a given xo
  node *semilla; // node from which we obtain its satellites to use them
                 // as candidates
  node_satelit *candidat; // possible cluster point sent to validate

  pila p_n;

  int vn_punts;
  node *topright;
  node *topleft;
  float *min;
  float *max;
  float *x_mean; // to compute xmid
  bool owns_x_mean;

  // modifiers

  void
  mstinsertar(node *pt); // mark as inside the MST and remove from `noin`
  void add_satelit(int, node *, node *); // link satellites

  // accessors
  int mstinsertat(node *pt); // check whether it is already in the MST

  // initialization
  // float *calculate_xomig_curve(); // search the curve origin point
  // closest to xmean
  void calcular_max_min_cluster();

  void obtener_quartiles(
      ll_q *ll_qt); // compute quartiles on distances obtained
                    // from the minimum spanning tree (MST) of points.
  float *obtener_satelites();

  // vector ops

  float *mult_esc(float e, float *v);
  float distancia(float *pnt1, float *pnt2);
  float *sum_v(float *v1, float *v2);

public:
  // constructors
  ll_p(int d);
  ~ll_p();

  void add_ordX_principal(float *vect);

  // initialization
  float inicialitzacio_principal();
  void inicialitzacio_final();
  void tornar_a_xomig(); // return to the origin point to process in opposite direction
  // for final spaces, this returns weighted xmean (h_tail weighting)

  // modifiers
  void
  trobar_primer_candidat_clt(float *xo); // search for the first candidate in the cluster.
  float *canviar_orientacio_clt(); // starts the search in the opposite direction, returns
                                   // the 1st candidate

  // accessors

  float *primer_candidat_clt(); // returns the first candidate in the cluster.
  float *seguent_candidat_clt(int validacio); // search the next candidate and
                                              // pass the previous validation
  int n_punts();                              // number of points in the cluster
  void donar_max_min_xomig(float **mx, float **mn, float **xm, float *s_d,
                           bool *owns_xm);

  // pointer-based queries

  node *resetpt() { return topleft->seg[DRETA]; }
  void *noend(node *pt);
  float *llpt(node *pt);
  void advpt(node **pt);
};
