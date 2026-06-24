#include "ll_flt.h"
#include "ll_pnt.h"
// #include "Ma.h"
#include "ll_p.h"
#include "mb.h"

extern "C" {
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
}

// #define PROF_REQ 1 // 0: compute curve on STV; 1: on subspace curves;
// 2: on subspace surfaces...

#define C_EPS 0.05 // Without recursion, the value 0.01 worked fine
#define NPTMIN 50  // minimum number of points to compute the curve
#define LD 0.5

#define PINZA_MAX PI / 4
#define PINZA_MIN PI / ((4 * NPARTS) + 1) // 90/4

class espai {
private:
  int Dim;
  int profundidad;
  // ### M_a *Ma;
  ll_p *ll_pt;

  float suma_d; // sum of distances in the MST.
  float h_tail; // distance threshold (half-thickness bandwidth) around the plane
  float delta;  // step size along the curve
  float *eps_x; // minimum distance required between cluster xo and xmean
                // for validity; shared by all clusters in the current space.
  float *xomig; // required for obtain_STV() and to obtain initial xo and b_opt.
  int xomig_offset; // leading coordinates owned before xomig; 0 means raw local coords.
  bool xomig_shared; // true when xomig aliases memory still owned elsewhere.
  // When computing the curve in the opposite direction, we do not necessarily
  // start from this xomig; we can start from its corresponding pop.

  int bficorba; // flag indicating whether a new point appears in the next pop.

  typedef struct m_d_s {
    float *xmean;
    float span;
    float density;
  } m_d_s;

  struct opt {
    float VTG;
    M_b *Mb;
    M_b *Mb_ant;
    espai *espai_;
    m_d_s mds;
  } optims;

  struct x {
    float *act;
    float *ant;
  } xo;

  struct pop {
    int depth;
    float *alpha;
    float I;
    float *b_ast;
    float var_k;
    float span;
    float density;
    espai *espai_;

    ~pop();
  };

  float Var_PC;
  float Var_res;
  float GTV;

  int dist_al_pla(float *n_punt);
  void calcular_htail_delta_xomig_epsx();
  espai *obtenir_cluster(M_b *Mb, m_d_s *mds);
  float *treure_coord(float *n_pnt); // project onto the plane before passing
                                     // to the subspace of a lower dimension
  int fi_corba(float *n_pnt);
  int no_creua_corba(float *pop);
  void calcular_Mb(int ejegir, M_b *Mb, float porcion_pinza);
  float calcular_corba_en_un_sentit();
  float calcular_corba_en_sentit_contrari();
  float finalitzacio(); /* returns the curve GTV */
  float *allargar(float *bopt);
  float obtenir_STV();
  float *obtenir_bo_inicial(float *alfa);
  float Bmst();
  float gammln(float xx);
  float kernel(float d);
  int major(float *v1, float *v2);
  void eigens(float *A, float *RR, float *E,
              int N); // Copyright 1973, 1991 by Stephen L. Moshier

  // vector ops
  float distancia(float *pnt1, float *pnt2);
  float *mult_esc(float e, float *v);
  float mult_v(float *v1, float *v2);
  float *sum_v(float *v1, float *v2);
  float *dif_v(float *v2, float *v1);
  float *norma_v(float *v);

  // Modified 2002-04-16: declare as static so values are shared across spaces.
  static int PROF_REQ;
  static int NPARTS;
  static float C_H;
  static float C_D; // always less than 0.5

public:
  espai(ll_p *ll_punts, int d, int p);
  ~espai();
  float obtenir_VTG(float **xm);
  void rebre_M_a(M_a *n_Ma); // pass the new Ma to the subspace.

  M_a *Ma;             // ###
  ll_pnt<pop> *ll_pop; // ###

  // Modified 2002-04-16: initialize tuning parameters once at the first space.
  void inicializar_nparts_ch_cd(int profreq, int nparts, float c_h, float c_d);
  void obtenir_data(float *, int *, int *, int);
};
