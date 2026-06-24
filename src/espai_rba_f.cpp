// ### r: cluster rectangular.
// ### b: initial bopt / previous bopt.
// ### a: a pop cannot go back.

#include "espai.h"
#include <Rcpp.h>
#include <cmath>
#include <cstdlib>
#include <vector>

#define PI M_PI

espai::espai(ll_p *ll_punts, int d, int p) {
  ll_pt = ll_punts;
  Dim = d;
  profundidad = p;
  suma_d = 0.0f;
  h_tail = 0.0f;
  delta = 0.0f;
  eps_x = NULL;
  xomig = NULL;
  xomig_offset = 0;
  xomig_shared = false;
  bficorba = FALSE;
  optims.VTG = 0.0f;
  optims.Mb = NULL;
  optims.Mb_ant = NULL;
  optims.espai_ = NULL;
  optims.mds.xmean = NULL;
  optims.mds.span = 0.0f;
  optims.mds.density = 0.0f;
  xo.act = NULL;
  xo.ant = NULL;
  Var_PC = 0.0f;
  Var_res = 0.0f;
  GTV = 0.0f;
  Ma = NULL;
  ll_pop = NULL;
}

espai::pop::~pop() {
  delete[] (alpha ? alpha - depth : NULL);
  delete[] (b_ast ? b_ast - depth : NULL);
  delete espai_;
}

espai::~espai() {

  delete ll_pt;
  ll_pt = NULL;

  delete Ma;
  Ma = NULL;

  delete[] eps_x;
  eps_x = NULL;

  if (!xomig_shared && xomig) {
    delete[] (xomig - xomig_offset);
  }
  xomig = NULL;
  xomig_offset = 0;
  xomig_shared = false;

  delete optims.Mb;
  optims.Mb = NULL;
  delete optims.Mb_ant;
  optims.Mb_ant = NULL;
  delete optims.espai_;
  optims.espai_ = NULL;
  delete[] optims.mds.xmean;
  optims.mds.xmean = NULL;
  delete[] xo.act;
  xo.act = NULL;
  delete[] xo.ant;
  xo.ant = NULL;

  // Ownership/lifetime notes:
  // ll_pt, xomig, eps_x, mds.xmean, and optims.mds.xmean are typically
  // released along direct/recursive finalization paths or opposite-direction
  // routines.
  // optims.Mb and optims.espai_ may remain referenced by output lists.
  // xo.act and xo.ant are released when curve-direction routines exit.

  /* clear subspace lists */

  /* delete output lists */

  delete ll_pop;
  ll_pop = NULL;
}

float espai::obtenir_VTG(float **xm) {
  if ((profundidad == PROF_REQ) || (Dim == 1) ||
      (ll_pt->n_punts() < NPTMIN * Dim)) { // verification missing: check whether
                                           // point count is sufficient for Dim.

    // Initialization to obtain the STV
    ll_pt->inicialitzacio_final();
    calcular_htail_delta_xomig_epsx();

    GTV = obtenir_STV();
    // printf(" %f \n",GTV); //###
    delete ll_pop;
    ll_pop = NULL;
  } else {

    ll_pop = new ll_pnt<pop>();

    // Initialization to obtain the GTV
    ll_pt->inicialitzacio_principal();
    calcular_htail_delta_xomig_epsx();

    static_cast<void>(
        calcular_corba_en_un_sentit()); // each subspace-curve computation consumes
                                        // and frees the passed ll_pt as needed
    /* compute in the opposite direction */
    static_cast<void>(
        calcular_corba_en_sentit_contrari()); // may not produce any valid pop.

    GTV = finalitzacio(); /* obtain curve GTV */
  }
  if (!xomig) {
    *xm = NULL;
  } else {
    float *returned_xomig = new float[Dim + 1]();
    if (xomig_offset == 0) {
      memmove(returned_xomig + 1, xomig, Dim * sizeof(float));
    } else {
      memmove(returned_xomig, xomig - 1, (Dim + 1) * sizeof(float));
    }
    *xm = returned_xomig;
  }
  // Keep ll_pt alive until obtenir_data(); xomig may alias a point it owns.
  delete[] eps_x;
  eps_x = NULL;
  return GTV;
}

espai *espai::obtenir_cluster(M_b *Mb, m_d_s *mds) {
  int validacio = 1;
  float *v_pnt;
  float *n_pnt = NULL;
  ll_p *n_ll_pt;
  float sum_w, w;

  v_pnt = ll_pt->primer_candidat_clt();
  if (!v_pnt) {
    return NULL;
  }
  n_pnt = Mb->aplicar(v_pnt); // distance check to the plane is deferred here;
                              // this path assumes dmax and xmean displacement
                              // are within the h_tail threshold.
  if (!n_pnt) {
    return NULL;
  }

  /* calculate cluster densities */

  w = kernel(n_pnt[X] / h_tail) *
      (*(v_pnt - 1)); // always abs(n_pnt[X]/h_tail) < 1
  if (w <= 0) {       // ###
    w = 0;
  } // ###
  sum_w = w;    // needed to compute smooth_mean and density
  n_pnt[X] = w; // In the height coordinate we store the weight.
  /* create and insert first point to the cluster */

  n_pnt = treure_coord(n_pnt);
  n_ll_pt = new ll_p(Dim - 1);
  n_ll_pt->add_ordX_principal(n_pnt);

  v_pnt = ll_pt->seguent_candidat_clt(validacio);
  while (v_pnt && fi_corba(v_pnt)) { // keep iterating until fi_corba rejects v_pnt;
                                     // v_pnt itself may be outside this cluster.
    n_pnt = Mb->aplicar(v_pnt);
    if ((validacio = dist_al_pla(n_pnt))) {

      /* calculate cluster density */

      w = kernel(n_pnt[X] / h_tail) * (*(v_pnt - 1));
      sum_w += w;   // needed to compute smooth_mean and density
      n_pnt[X] = w; // In the height coordinate we store the weight.

      /* add point to cluster */

      n_pnt = treure_coord(n_pnt);
      n_ll_pt->add_ordX_principal(n_pnt);
    } else {
      delete[] n_pnt;
    }
    v_pnt = ll_pt->seguent_candidat_clt(validacio);
  }
  if (v_pnt) {
    while (v_pnt) {
      n_pnt = Mb->aplicar(v_pnt);
      if ((validacio = dist_al_pla(n_pnt))) {

        /* calculate cluster density */
        w = kernel(n_pnt[X] / h_tail) * (*(v_pnt - 1));
        sum_w += w;   // needed to compute smooth_mean and density
        n_pnt[X] = w; // In the height coordinate we store the weight.

        /* insert new point to cluster */
        n_pnt = treure_coord(n_pnt);
        n_ll_pt->add_ordX_principal(n_pnt);
      } else {
        delete[] n_pnt;
      }
      v_pnt = ll_pt->seguent_candidat_clt(
          validacio); // termination is controlled by candidate traversal state.
    }

    /* final calculation of densities */
    mds->span = (float)n_ll_pt->n_punts() / (float)ll_pt->n_punts();
    mds->density = sum_w / (ll_pt->n_punts() * h_tail);

    return (new espai(n_ll_pt, Dim - 1, profundidad + 1));
  } else {
    delete n_ll_pt;
    return NULL; // all the candidates have previously been inserted in other
                 // clusters, there are no new points to deal with.
  }
}

int espai::dist_al_pla(
    float *n_punt) { // check if the point is within distance h_tail of the plane
  return (fabs(n_punt[X]) < h_tail);
}

float *
espai::treure_coord(float *n_pnt) { // project onto the plane before passing
                                    // to the lower-dimensional subspace
  return n_pnt + 1; // the extra leading float is released by ll_p
                    // when points are deleted; also used by obtain_STV().
}

void espai::rebre_M_a(M_a *n_Ma) { // pass M_a to the lower subspace
                                   // after obtaining the best space
  Ma = n_Ma;
}

int espai::fi_corba(float *v_pnt) {
  float *Mba_pnt;
  if (bficorba) {
    Mba_pnt = optims.Mb_ant->aplicar(v_pnt);
    if (Mba_pnt[X] > 2 * delta)
      bficorba = FALSE; // uses 2*delta threshold instead of h_tail.
    //if (Mba_pnt[X]>h_tail) bficorba = FALSE;
    else {
      delete[] Mba_pnt;
      return TRUE;
    }
    delete[] Mba_pnt;
  }
  return FALSE;
}

int espai::no_creua_corba(float *ncand) {
  // Crossings depend on delta. Smaller delta gives better precision, both
  // when defining the curve and when finalizing it.
  float *punt;

  auto pt = ll_pop->resetpt();
  while (pt->seg->seg) {
    if ((distancia(((pop *)ll_pop->llpt(pt))->alpha, ncand)) <
        delta) { // if the curve is a circle, the last point is at distance delta
                 // from xomig.
      punt = optims.Mb_ant->aplicar(((pop *)ll_pop->llpt(pt))->alpha);
      if (punt[X] > 0) {
        delete[] punt;
        return FALSE;
      }
      delete[] punt;
    }
    ll_pop->advpt(&pt);
  }
  return TRUE;
}

void espai::calcular_Mb(int ejegir, M_b *Mb, float porcion_pinza) {
  int i; /* ejegir selects the fixed axis; porcion_pinza is constant through recursion */
  float VTG;
  m_d_s mds;
  espai *espai;
  if (ejegir) {
    /* recursive case */
    for (i = -NPARTS / 2; i < 0; i++) {
      calcular_Mb(ejegir - 1, Mb->girar(ejegir, i * porcion_pinza),
                  porcion_pinza);
    }
    for (i = 1; i <= NPARTS / 2; i++) {
      calcular_Mb(ejegir - 1, Mb->girar(ejegir, i * porcion_pinza),
                  porcion_pinza);
    }
    calcular_Mb(ejegir - 1, Mb, porcion_pinza);
  }
  else {
    /* direct case */
    mds.xmean = NULL; // in case we are at the end of the curve.
    Mb->calcular_la_inversa();
    espai = obtenir_cluster(Mb, &mds); // check whether the curve has ended
                                       // in this direction.
    if (espai) {
      espai->rebre_M_a(
          Mb->donar_M_a(Ma)); // pass Ma positioned in original coordinates
                              // for the Mb plane.
      VTG = espai->obtenir_VTG(&mds.xmean);
      if (VTG < optims.VTG) {
        delete optims.Mb; // delete previous non-optimal Mb candidate.
        delete optims.espai_;
        optims.VTG = VTG;
        optims.Mb = Mb;
        optims.espai_ = espai;
        delete[] optims.mds.xmean;
        optims.mds.xmean = Mb->desaplicar(mds.xmean); // recover candidate pop
        optims.mds.span = mds.span;
        optims.mds.density = mds.density;

      } else {
        delete Mb;
        delete espai;
      }
    } // Even if this cluster has no more points, other clusters may still.
    else {
      delete Mb;
    }
    delete[] mds.xmean;
  }
}

float espai::calcular_corba_en_un_sentit() {
  // Compared with calcular_corba_en_sentit_contrari(), this differs only in
  // initialization and insertion order (`add` vs `addrev`).
  // `b_ast` keeps eigenvector orientation; arc-length distances have opposite
  // signs across the two traversal directions.
  pop *n_pop;
  float pinza;
  float *n_bo, *b_opt;
  float *v_xact2xm = NULL;
  float *v_xant2xm = NULL;
  int naux_delta;
  float sum = 0;
  float lambda;
  float alfa;
  // char c; //###
  // int i; //file output

  /* find the first principal component of the points */

  optims.Mb = new M_b(Dim, obtenir_bo_inicial(&alfa)); // ###
  // optims.Mb = new M_b(Dim,mult_esc(-1,get_bo_initial()));

  /* compute h_tail threshold (plane bandwidth) and delta step */

  lambda = suma_d / (pow(ll_pt->n_punts(), ((float)(Dim - 1)) / Dim) * Bmst()) *
           pow(((float)(Dim - 1)) / (1 - ((1 - LD) * alfa + LD)),
               ((float)(Dim - 1)) / (2 * Dim)) *
           (1 / sqrt(2 * PI)) * pow(((float)(Dim - 1)) / Dim, Dim / 2.);

  h_tail = C_H * 2.214 * pow(4. / 3, 1. / 5) * lambda *
           pow(ll_pt->n_punts(),
               -1. / 5); // 2.214 comes from the Epanechnikov kernel normalization.
  delta = C_D * h_tail;  // Always less than h_tail. This keeps the new xo
                         // within distance dmax of an ll_p point.

  /* find the midpoint pop of the curve */

  optims.Mb_ant = optims.Mb->replicar(); // Mb_ant is required for fi_corba()
                                         // checks until direction changes.
  optims.Mb_ant->calcular_la_inversa();
  n_bo = mult_esc(
      -2 * h_tail,
      optims.Mb_ant->donar_bopt()); // Mb_ant must receive xo.ant for fi_corba().
  xo.ant = sum_v(xomig, n_bo);
  optims.Mb_ant->rebre_xo(xo.ant);
  delete[] n_bo;

  optims.espai_ = NULL;
  xo.act = NULL;
  optims.mds.xmean = new float[Dim];
  memmove(optims.mds.xmean, xomig,
          Dim * sizeof(float)); // copy current xomig (needed by ll_p)
  do {
    delete optims.espai_;
    delete[] v_xact2xm;
    delete[] xo.act;
    xo.act = optims.mds.xmean;
    optims.Mb->rebre_xo(
        xo.act); // this Mb is immediately deleted (VTG=INF), but it
                 // propagates xo_act across candidate matrices per rotation.

    /* compute optimal GTV/VTG for the given xo */

    optims.VTG = INF;        // optims.Mb starts as the previous cluster optimum
                             // and is replaced in calcular_Mb when a better
                             // candidate is found.
    optims.espai_ = NULL;    // prior optimal space/Mb may still be referenced
                             // by result lists, so they are not deleted here.
    optims.mds.xmean = NULL; // xo passed to Mb cannot be freed until replaced.
    pinza = PI / 2;
    while (pinza > PI / 16) {
      calcular_Mb(Dim - 1, optims.Mb->replicar(),
                  pinza / NPARTS); // use one subdivision of current angle window
      pinza = pinza / NPARTS;      // next window equals one previous subdivision
    }
    v_xact2xm = dif_v(optims.mds.xmean, xo.act);
  } while (major(v_xact2xm, eps_x));

  if (!optims.mds.xmean) {
    delete optims.espai_;
    optims.espai_ = NULL;
    delete optims.Mb_ant;
    optims.Mb_ant = NULL;
    delete optims.Mb;
    optims.Mb = NULL;
    delete[] v_xact2xm;
    v_xact2xm = NULL;
    delete[] xo.act;
    xo.act = NULL;
    delete[] xo.ant;
    xo.ant = NULL;
    return sum;
  }

  delete[] xo.act; // delete previous candidate; copied xomig is now consumed.
  xo.act = optims.mds.xmean;

  optims.mds.xmean =
      allargar(optims.mds.xmean); // extend optimal xmean for storage in outputs.
  xomig = optims.mds.xmean;       // new midpoint pop of the curve; do not lose it.
  xomig_offset = profundidad;
  xomig_shared = true;

  optims.Mb->rebre_xo(optims.mds.xmean); // pass extended xo to optimal Mb;
                                         // Ma is computed from it.
  // optims.espai_->rebre_M_a(optims.Mb->donar_M_a(Ma));
  // Pass Ma positioned in original coordinates for this Mb plane.

  n_pop = new pop;
  n_pop->depth = profundidad;
  sum = 0;
  n_pop->I = 0;
  delete[] xo.ant; // discard xo of previous Mb_ant.
  xo.ant = xo.act; // xo.ant will be released after the next cluster step.

  b_opt = optims.Mb->donar_bopt(); // used to advance the cluster
  delete optims.Mb_ant;
  optims.Mb_ant = optims.Mb;
  optims.Mb = optims.Mb->replicar(); // keep a copy as initial matrix
                                     // for the next cluster step.

  n_pop->alpha = optims.mds.xmean; // output points/vectors are Dim+depth.
  n_pop->b_ast = allargar(b_opt);  // store extended b_opt (non-normalized is expected).
  n_pop->var_k = optims.VTG;
  n_pop->span = optims.mds.span;
  n_pop->density = optims.mds.density;
  if (PROF_REQ == 2) {
    n_pop->espai_ = optims.espai_;
  } else {
    delete optims.espai_;
    n_pop->espai_ = NULL;
  }
  optims.espai_ = NULL;

  ll_pop->add(n_pop);
  /* update cluster */

  n_bo = mult_esc(delta, b_opt);
  optims.mds.xmean =
      sum_v(xo.act, n_bo); // preload xmean; it becomes xo.act at next iteration.

  xo.act = NULL; // xo.ant now owns this pointer and will free it later.
  delete[] n_bo;

  /* find subsequent pops in the forward direction of cluster b */

  xo.act = NULL;
  naux_delta = 1;
  // Cap at n_punts: a principal curve cannot meaningfully have more oriented
  // points than the dataset has points. Without this bound, near-zero delta
  // values (which can arise from extreme scale anisotropy causing alfa -> 1 in
  // the bandwidth formula) would cause the loop to run for an unbounded number
  // of steps before the natural termination conditions fire.
  int max_steps = ll_pt->n_punts();
  int n_steps = 0;
  while (no_creua_corba(xo.ant) && n_steps < max_steps) {
    n_steps++;
    /* search for pop */
    bficorba = TRUE; // remains true until no further point can be processed.
    delete[] v_xact2xm;
    delete[] v_xant2xm;
    delete[] xo.act;           // reset xo.act before assigning current xmean.
    xo.act = optims.mds.xmean; // keep current xmean as active xo candidate.
    optims.Mb->rebre_xo(
        xo.act); // this Mb is immediately deleted (VTG=INF), but it
                 // propagates xo_act across candidate matrices per rotation.
    ll_pt->trobar_primer_candidat_clt(
        optims.mds.xmean); // search for the first candidate in the cluster.

    /* compute optimal GTV/VTG for the given xo */

    optims.VTG = INF;        // same initialization logic as above.
    optims.espai_ = NULL;    // preserve previously referenced space objects.
    optims.mds.xmean = NULL; // xo ownership moves only after replacement.
    pinza = PINZA_MAX;
    while (pinza > PINZA_MIN) {
      calcular_Mb(Dim - 1, optims.Mb->replicar(),
                  pinza / NPARTS); // evaluate one subdivision of the angle window
      if (optims.VTG == INF) { // all hyperplanes were already processed in
                               // previous clusters; no points remain.
        delete optims.Mb;
        optims.Mb = NULL;
        delete optims.Mb_ant;
        optims.Mb_ant = NULL;
        delete[] xo.act;
        xo.act = NULL;
        delete[] optims.mds.xmean;
        optims.mds.xmean = NULL;
        // xo.act is freed in the opposite-direction routine.
        delete[] xo.ant;
        xo.ant = NULL;
        return (sum); // exit when there are no more points to process
      }
      pinza = pinza / NPARTS; // shrink to one subdivision of previous window
    }
    // if (depth == 0){
    // c = getchar(); //###
    // }
    v_xact2xm = dif_v(optims.mds.xmean, xo.act);
    v_xant2xm = dif_v(optims.mds.xmean, xo.ant);
    if ((mult_v(optims.Mb->donar_bopt(), optims.Mb_ant->donar_bopt()) < 0) ||
        (mult_v(v_xant2xm, optims.Mb_ant->donar_bopt()) <
         0)) { // If angle < 0, we are changing curve direction
               // (cosine theorem). We must avoid this.
      /* Advance original xo slightly more than at initial advance. */
      naux_delta++;
      delete[] optims.mds.xmean;
      n_bo = mult_esc(naux_delta * delta, optims.Mb_ant->donar_bopt());
      optims.mds.xmean =
          sum_v(xo.ant, n_bo); // increase step when candidate would switch to opposite direction.
      delete[] n_bo;
      delete optims.Mb;
      optims.Mb = optims.Mb_ant->replicar(); // keep original Mb_ant (stored in b_ast);
                                             // assign updated xo to this copy.
    } else if (major(v_xact2xm,
                     eps_x)) { // check that xo and xmean are not too far apart.
                               // If xmean was moved, this condition should still
                               // hold due to naux_delta; remove v_xant2xm.
      /* Values are invalid: search for another cluster starting from xmean. */
      delete optims.espai_;
      optims.espai_ = NULL;
    } else {
      /* Values are valid. */
      /* update lists */
      delete[] xo.act;
      xo.act = optims.mds.xmean;
      optims.mds.xmean =
          allargar(optims.mds.xmean); // store extended optimal xmean in output list.

      optims.Mb->rebre_xo(
          optims.mds.xmean); // pass extended xo to optimal Mb for Ma computation.
      // optims.espai_->rebre_M_a(optims.Mb->donar_M_a(Ma));
      // Pass Ma positioned in original coordinates for this Mb plane.

      n_pop = new pop;
      n_pop->depth = profundidad;
      sum += distancia(xo.act, xo.ant);
      n_pop->I = sum;
      delete[] xo.ant;
      xo.ant = xo.act; // no separate copy needed; previous value is no longer used.

      b_opt = optims.Mb->donar_bopt(); // direction vector used to advance cluster.
      delete optims.Mb_ant; // also releases previous pop b_opt.
      optims.Mb_ant = optims.Mb;
      optims.Mb = optims.Mb->replicar(); // keep copy as initial matrix
                                         // for next cluster step.

      n_pop->alpha = optims.mds.xmean; // output points/vectors are Dim+depth.
      n_pop->b_ast = allargar(b_opt);  // keep extended (possibly non-normalized) b_opt.
      n_pop->var_k = optims.VTG;
      n_pop->span = optims.mds.span;
      n_pop->density = optims.mds.density;
      if (PROF_REQ == 2) {
        n_pop->espai_ = optims.espai_;
      } else {
        delete optims.espai_;
        n_pop->espai_ = NULL;
      }
      optims.espai_ = NULL;

      ll_pop->add(n_pop);

      /* update cluster */

      n_bo = mult_esc(delta, b_opt);
      // delete optims.mds.xmean; // point to the pop
      // saved in alpha. It cannot be deleted.
      optims.mds.xmean = sum_v(xo.act,
                               n_bo); // preload xmean for next loop iteration.

      xo.act = NULL; // xo.ant now owns this pointer and will free it later.
      // delete b_opt; it is deleted together with its Mb.
      // Safe because b_ast stores an extended copy.
      delete[] n_bo;
      naux_delta = 1;
    }
  }
  delete optims.Mb;
  optims.Mb = NULL;
  delete optims.Mb_ant;
  optims.Mb_ant = NULL;
  delete[] optims.mds.xmean;
  optims.mds.xmean = NULL;
  // xo.act is deleted in the opposite-direction routine.
  delete[] xo.ant;
  xo.ant = NULL;
  return (sum); // exit when there is a crossing with the curve itself
}

float espai::calcular_corba_en_sentit_contrari() {
  // Differs from the previous function only in initialization and
  // insertion into output lists via addrev(). b_ast vectors remain
  // eigenvector-oriented and I distances are negative.
  pop *n_pop;
  float pinza;
  float *n_bo, *b_opt;
  float *v_xact2xm = NULL;
  float *v_xant2xm = NULL;
  int naux_delta;
  // Same cap as in the forward direction; see comment there.
  int max_steps = ll_pt->n_punts();
  int n_steps = 0;
  float sum = 0;

  ll_pt->tornar_a_xomig();

  delete optims.Mb;
  optims.Mb = NULL;
  delete[] optims.mds.xmean;
  optims.mds.xmean = NULL;
  delete[] xo.act;
  xo.act = NULL;

  auto pt = ll_pop ? ll_pop->resetpt() : nullptr;
  if (!pt || !ll_pop->llpt(pt) || !((pop *)ll_pop->llpt(pt))->b_ast) {
    return sum;
  }

  optims.Mb = new M_b(
      Dim,
      mult_esc(-1,
               ((pop *)ll_pop->llpt(pt))->b_ast)); // matrix for previous pop, with
                                                   // opposite orientation
  optims.Mb_ant = optims.Mb->replicar();
  optims.Mb_ant->calcular_la_inversa();
  xo.ant = new float[Dim];
  memmove(xo.ant, xomig,
          Dim * sizeof(float)); // keep explicit copy: alpha-owned memory cannot
                                // be reused directly for distance updates.
  optims.Mb_ant->rebre_xo(xo.ant);

  /* advance cluster */
  n_bo = mult_esc(delta,
                  optims.Mb->donar_bopt()); // advance in opposite direction; start
                                            // from the neighbor of the xomig pop.
  optims.mds.xmean =
      sum_v(xo.ant,
            n_bo); // preload xmean; becomes xo.act at start of next iteration.

  // xo.act = NULL; the first xo.act to be freed is
  // the last one from the opposite-direction curve.
  delete[] n_bo;
  naux_delta = 1;
  n_steps = 0;

  while (no_creua_corba(xo.ant) && n_steps < max_steps) {
    n_steps++;
    /* search for pop */
    bficorba = TRUE; // remains true until no further point can be processed.
    delete[] v_xact2xm;
    delete[] v_xant2xm;
    delete[] xo.act;           // reset xo.act before assigning current xmean.
    xo.act = optims.mds.xmean; // keep current xmean as active xo candidate.
    optims.Mb->rebre_xo(
        xo.act); // this Mb is immediately deleted (VTG=INF), but it
                 // propagates xo_act across candidate matrices per rotation.
    ll_pt->trobar_primer_candidat_clt(
        optims.mds.xmean); // search for the first candidate in the cluster.

    /* compute optimal GTV/VTG for the given xo */

    optims.VTG = INF;        // same initialization logic as forward direction.
    optims.espai_ = NULL;    // preserve previously referenced space objects.
    optims.mds.xmean = NULL; // xo ownership moves only after replacement.
    pinza = PINZA_MAX;
    while (pinza > PINZA_MIN) {
      calcular_Mb(Dim - 1, optims.Mb->replicar(),
                  pinza / NPARTS); // evaluate one subdivision of the angle window
      if (optims.VTG == INF) { // for all possible hyperplanes, cluster points
                               // were already processed in previous clusters;
                               // we must continue in opposite direction, else
                               // we would exit at the first window
        delete optims.Mb;
        optims.Mb = NULL;
        delete optims.Mb_ant;
        optims.Mb_ant = NULL;
        delete[] xo.act;
        xo.act = NULL;
        delete[] optims.mds.xmean;
        optims.mds.xmean = NULL;
        delete[] xo.ant;
        xo.ant = NULL;
        return (sum); // exit when there are no more points to process
      }
      pinza = pinza / NPARTS; // shrink to one subdivision of previous window
    }

    v_xact2xm = dif_v(optims.mds.xmean, xo.act);
    v_xant2xm = dif_v(optims.mds.xmean, xo.ant);
    if (mult_v(optims.Mb->donar_bopt(), optims.Mb_ant->donar_bopt()) < 0 ||
        mult_v(v_xant2xm, optims.Mb_ant->donar_bopt()) <
            0) { // If angle < 0, we are changing curve direction. Avoid this.
      /* Advance original xo a bit more than in the previous iteration. */
      naux_delta++;
      delete[] optims.mds.xmean;
      n_bo = mult_esc(naux_delta * delta, optims.Mb_ant->donar_bopt());
      optims.mds.xmean = sum_v(xo.ant, n_bo);
      delete[] n_bo;
      delete optims.Mb;
      optims.Mb = optims.Mb_ant->replicar(); // keep original Mb_ant (stored in b_ast);
                                             // assign updated xo to this copy.
    } else if (major(v_xact2xm,
                     eps_x)) { // check that xo and xmean are not too far apart.
                               // If xmean was moved, this condition should still
                               // hold due to naux_delta; remove v_xant2xm.
      /* Values are invalid: search for another cluster starting from xmean. */
      delete optims.espai_;
      optims.espai_ = NULL;
    } else {
      /* Values are valid. */
      /* update lists */
      delete[] xo.act;
      xo.act = optims.mds.xmean;
      optims.mds.xmean =
          allargar(optims.mds.xmean); // store extended optimal xmean in output list.

      optims.Mb->rebre_xo(
          optims.mds.xmean); // pass extended xo to optimal Mb for Ma computation.
      // optims.espai_->rebre_M_a(optims.Mb->donar_M_a(Ma));
      // Pass Ma positioned in original coordinates for this Mb plane.

      b_opt = optims.Mb->donar_bopt(); // direction vector used to advance cluster.
      n_bo = mult_esc(-1, b_opt);      // all vectors in b_ast
                                       // they will have the same orientation.
      delete optims.Mb_ant;            // also releases previous pop b_opt.
      optims.Mb_ant = optims.Mb;
      optims.Mb = optims.Mb->replicar(); // keep copy as initial matrix
                                         // for next cluster step.

      n_pop = new pop;
      n_pop->depth = profundidad;
      sum -= distancia(xo.act, xo.ant); // negative distances in opposite direction.
      n_pop->I = sum;
      delete[] xo.ant;
      xo.ant = xo.act; // no separate copy needed; previous value is no longer used.

      n_pop->alpha = optims.mds.xmean; // output points/vectors are Dim+depth.
      n_pop->b_ast =
          allargar(n_bo); // output points/vectors are Dim+depth; non-normalized
                          // b_opt is preserved with opposite orientation.
      delete[] n_bo;
      n_pop->var_k = optims.VTG;
      n_pop->span = optims.mds.span;
      n_pop->density = optims.mds.density;
      if (PROF_REQ == 2) {
        n_pop->espai_ = optims.espai_;
      } else {
        delete optims.espai_;
        n_pop->espai_ = NULL;
      }
      optims.espai_ = NULL;

      ll_pop->addrev(n_pop);

      /* update cluster */

      n_bo = mult_esc(delta, b_opt);
      // delete optims.mds.xmean; // point to the pop
      // saved in alpha. It cannot be deleted.
      optims.mds.xmean = sum_v(xo.act,
                               n_bo); // preload xmean for next loop iteration.

      xo.act = NULL; // xo.ant now owns this pointer and will free it later.
      // b_opt is deleted together with its Mb; b_ast stores an extended copy.
      delete[] n_bo;
      naux_delta = 1;
    }
  }
  delete optims.Mb;
  optims.Mb = NULL;
  delete optims.Mb_ant;
  optims.Mb_ant = NULL;
  delete[] optims.mds.xmean;
  optims.mds.xmean = NULL;
  delete[] xo.act;
  xo.act = NULL;
  delete[] xo.ant;
  xo.ant = NULL;
  return (sum); // exit when there is a crossing with the curve itself
}

float *espai::allargar(float *b_opt) {
  int i;
  float *n_b;
  n_b = new float[(Dim + profundidad)]();

  for (i = 0; i < profundidad; i++)
    n_b[i] = 0.0;

  if (b_opt) {
    memmove(n_b + profundidad, b_opt,
            Dim * sizeof(float)); // copy original Dim coordinates
  }

  // for(i=depth;i<Dim+depth;i++)
  // n_b[i]=b_opt[i-profundidad];

  return n_b + profundidad;
}

float *espai::obtenir_bo_inicial(float *alfa) {

  float *A = new float[Dim * (Dim + 1) / 2];
  float *EV = new float[Dim * Dim]();
  float *E = new float[Dim];
  float *b_op = new float[Dim];
  float *punt;
  int np = ll_pt->n_punts();
  int i, j;
  float VarTot;
  std::vector<double> m(Dim, 0.0);
  std::vector<double> S(Dim * Dim, 0.0);

  auto pt = ll_pt->resetpt();
  while (ll_pt->noend(pt)) {
    punt = ((float *)ll_pt->llpt(pt));
    for (i = 0; i < Dim; i++) {
      m[i] += punt[i];
    }
    ll_pt->advpt(&pt);
  }
  for (i = 0; i < Dim; i++) {
    m[i] /= np;
  }

  pt = ll_pt->resetpt();
  while (ll_pt->noend(pt)) {
    punt = ((float *)ll_pt->llpt(pt));
    for (i = 0; i < Dim; i++) {
      double centered_i = punt[i] - m[i];
      for (j = 0; j <= i; j++) {
        S[(i * Dim) + j] += centered_i * (punt[j] - m[j]);
      }
    }
    ll_pt->advpt(&pt);
  }

  for (i = 0; i < Dim; i++) {
    for (j = 0; j < Dim; j++) {
      double value = (i >= j) ? S[(i * Dim) + j] : S[(j * Dim) + i];
      A[((i * i + i) / 2) + j] =
          static_cast<float>(value / np); /* lower-triangular packed storage order */
    }
  }

  eigens(A, EV, E, Dim);

  i = 0;
  VarTot = 0;
  for (j = 0; j < Dim; j++) {
    VarTot += E[j];
    if (E[i] < E[j])
      i = j; // we look for the eigenvector with the highest eigenvalue
  }
  for (j = 0; j < Dim; j++)
    b_op[j] = EV[Dim * i + j]; /* EV[i][j] */

  *alfa = (VarTot > 0.0f) ? E[i] / VarTot : 0.0f;

  delete[] A;
  delete[] EV;
  delete[] E;

  return b_op;
}

float espai::Bmst() {
  float cd, S, fact;
  // float ubeta, lbeta;
  double par_f, res;
  int i;

  par_f = (Dim / 2.) + 1;
  cd = pow(PI, Dim / 2.) / exp(gammln(par_f));
  par_f = 1. / Dim;
  // lbeta = exp(gammln(par_f))/(Dim*pow(cd,1./Dim));
  // ubeta = pow(2,1./Dim)*lbeta;

  S = 0;
  fact = 1.;
  for (i = 1; i <= 30; i++) {
    fact = fact * i;
    par_f = i + (1. / Dim) - 1;
    S += exp(gammln(par_f)) / (fact * pow(i, (1. / Dim) + 1));
  }

  res = S / (Dim * pow(cd, 1. / Dim));
  return res;
}

float espai::gammln(float xx) {
  double x, y, tmp, ser;
  static double cof[6] = {76.18009172947146,     -86.50532032941677,
                          24.01409824083091,     -1.231739572450155,
                          0.1208650973866179e-2, -0.5395239384953e-5};
  int j;
  y = x = xx;
  tmp = x + 5.5;
  tmp -= (x + 0.5) * log(tmp);
  ser = 1.00000000190015;
  for (j = 0; j <= 5; j++)
    ser += cof[j] / ++y;
  return -tmp + log(2.5066282746310005 * ser / x);
}

float espai::obtenir_STV() {
  // When curve generation fails, STV is computed around center point xomig.

  int j;
  float *punt;
  float *c_xomig;
  float sum_w = 0.0;
  float stv = 0.0;

  c_xomig = new float[(Dim + 1)];
  c_xomig[0] = 0;
  c_xomig++; // keep +1 offset convention when passing to upper space
  memmove(c_xomig, xomig, Dim * sizeof(float));
  delete[] xomig;
  xomig = c_xomig;
  xomig_offset = 1;
  xomig_shared = false;

  auto pt = ll_pt->resetpt();
  while (ll_pt->noend(pt)) {
    punt = ll_pt->llpt(pt);
    // *(punt-1) stores precomputed weight.
    sum_w += *(punt - 1);
    for (j = 0; j < Dim; j++)
      stv += pow(punt[j] - xomig[j], 2) *
             (*(punt - 1)); // weight already encodes h_tail/kernel scaling
    ll_pt->advpt(&pt);
  }
  return stv / sum_w;
}

void espai::calcular_htail_delta_xomig_epsx() {
  float *range, *min, *max;
  bool owns_xomig;

  ll_pt->donar_max_min_xomig(&max, &min, &xomig, &suma_d, &owns_xomig);
  xomig_offset = 0;
  xomig_shared = !owns_xomig;

  range = dif_v(max, min);
  eps_x = mult_esc(C_EPS, range);

  delete[] range;
  delete[] max;
  delete[] min;
}

float espai::kernel(float d) { return (0.75 * (1 - pow(d, 2))); }

int espai::major(float *v1, float *v2) {
  int i = 0;

  while (i < Dim && fabs(v1[i]) <= fabs(v2[i]))
    i++; // greater or equal
  return (i != Dim);
}

float espai::finalitzacio() {

  ll_flt w_s;
  float *xomig_0 = NULL, *xomig_1 = NULL;
  float Iant, I2ant, dant, wsact;
  float sum = 0;
  pop *p0;

  if (!ll_pop) {
    return INF;
  }

  auto pt = ll_pop->resetpt();
  if (!pt || !(p0 = (pop *)ll_pop->llpt(pt))) {
    return INF;
  }

  // Degenerate curve with a single POP: keep a valid xomig and GTV.
  if (!ll_pop->noend(pt)) {
    xomig = p0->alpha;
    xomig_offset = profundidad;
    xomig_shared = true;
    p0->I = 0.0;
    p0->density = 1.0;
    Var_PC = 0.0;
    Var_res = p0->var_k;
    return (Var_PC + Var_res);
  }

  /* w_s :: ws(n) = (I(n+1)-I(n-1))*density(n) */
  Iant = ((pop *)ll_pop->llpt(pt))->I;
  dant = ((pop *)ll_pop->llpt(pt))->density;
  ll_pop->advpt(&pt);
  wsact = (((pop *)ll_pop->llpt(pt))->I - Iant) * 2 * dant;
  sum += wsact;
  w_s.add(wsact);
  while (pt->seg->seg) {
    I2ant = Iant;
    Iant = ((pop *)ll_pop->llpt(pt))->I;
    dant = ((pop *)ll_pop->llpt(pt))->density;
    ll_pop->advpt(&pt);
    wsact = (((pop *)ll_pop->llpt(pt))->I - I2ant) * dant;
    sum += wsact;
    w_s.add(wsact);
  }
  wsact = (((pop *)ll_pop->llpt(pt))->I - Iant) * 2 *
          ((pop *)ll_pop->llpt(pt))->density;
  sum += wsact;
  w_s.add(wsact);

  /* w_s :: w_s(n)/sum(w_s) */

  auto ptws = w_s.resetpt();
  while (w_s.noend(ptws)) {
    w_s.modpt(ptws, w_s.llpt(ptws) / sum);
    w_s.advpt(&ptws);
  }

  /* density :: density(n) = 2*w_s(n)/(I(n+1)-I(n-1)) */

  sum = 0;
  pt = ll_pop->resetpt();
  ptws = w_s.resetpt();
  auto ptprev = pt;
  auto ptnext = pt;
  ll_pop->advpt(&ptnext);
  auto cur = (pop *)ll_pop->llpt(pt);
  auto next = (pop *)ll_pop->llpt(ptnext);
  sum += cur->I * w_s.llpt(ptws);
  cur->density = (w_s.llpt(ptws)) / (next->I - cur->I);

  ptprev = pt;
  pt = ptnext;
  w_s.advpt(&ptws);
  while (pt->seg->seg) {
    ptnext = pt;
    ll_pop->advpt(&ptnext);
    auto prev = (pop *)ll_pop->llpt(ptprev);
    cur = (pop *)ll_pop->llpt(pt);
    next = (pop *)ll_pop->llpt(ptnext);
    sum += cur->I * w_s.llpt(ptws);
    cur->density = (2 * w_s.llpt(ptws)) / (next->I - prev->I);
    ptprev = pt;
    pt = ptnext;
    w_s.advpt(&ptws);
  }
  auto prev = (pop *)ll_pop->llpt(ptprev);
  cur = (pop *)ll_pop->llpt(pt);
  sum += cur->I * w_s.llpt(ptws);
  cur->density = (w_s.llpt(ptws)) / (cur->I - prev->I);

  /* I :: I(n) = I(n) - sum(I*ws) */

  /* Var_PC :: principal curve variance */
  /* Var_res :: residual variance */

  Var_PC = 0.0;
  Var_res = 0.0;
  ptws = w_s.resetpt();
  pt = ll_pop->resetpt();
  ((pop *)ll_pop->llpt(pt))->I = ((pop *)ll_pop->llpt(pt))->I - sum;

  auto ptant = ll_pop->resetpt();
  while (((pop *)ll_pop->llpt(pt))->I <
         0) { // there will always be a change from positive to negative
    Var_PC += pow(((pop *)ll_pop->llpt(pt))->I, 2) * w_s.llpt(ptws);
    Var_res += ((pop *)ll_pop->llpt(pt))->var_k * w_s.llpt(ptws);
    ptant = pt;
    ll_pop->advpt(&pt);
    w_s.advpt(&ptws);
    ((pop *)ll_pop->llpt(pt))->I = ((pop *)ll_pop->llpt(pt))->I - sum;
  }
  /* assign curve-center pop to pass to the upper space */

  if (!((pop *)ll_pop->llpt(pt))->I) {
    float *center = mult_esc(((pop *)ll_pop->llpt(pt))->I,
                             ((pop *)ll_pop->llpt(ptant))->alpha);
    xomig_0 = mult_esc(((pop *)ll_pop->llpt(ptant))->I,
                       ((pop *)ll_pop->llpt(pt))->alpha);
    xomig_1 = sum_v(center, xomig_0);
    delete[] xomig_0;
    delete[] center;
    center =
        mult_esc(((pop *)ll_pop->llpt(pt))->I * ((pop *)ll_pop->llpt(ptant))->I,
                 xomig_1);
    xomig = allargar(center);
    xomig_offset = profundidad;
    xomig_shared = false;
    delete[] xomig_1;
    delete[] center;
  } else {
    xomig = ((pop *)ll_pop->llpt(pt))->alpha;
    xomig_offset = profundidad;
    xomig_shared = true;
  }

  ll_pop->advpt(&pt);
  w_s.advpt(&ptws);
  while (w_s.noend(ptws)) {
    ((pop *)ll_pop->llpt(pt))->I = ((pop *)ll_pop->llpt(pt))->I - sum;
    Var_PC += pow(((pop *)ll_pop->llpt(pt))->I, 2) * w_s.llpt(ptws);
    Var_res += ((pop *)ll_pop->llpt(pt))->var_k * w_s.llpt(ptws);
    ll_pop->advpt(&pt);
    w_s.advpt(&ptws);
  }

  /* GTV:: global total variance */

  // printf("var_PC: %f ; var_res: %f \n", Var_PC, Var_res);
  return (Var_PC + Var_res); // GTV
}

void espai::obtenir_data(float *result, int *ncol, int *nrow,
                         int max_rows) {

  int i;
  *ncol = Dim * 2 + 5;
  (*nrow) = 0;
  espai *sespai;
  float *auxa;
  float *auxb;
  auto ensure_row_capacity = [&]() {
    if (*nrow >= max_rows) {
      Rcpp::stop("Internal Rpcop output buffer is too small.\n");
    }
  };
  auto singleton_variance = [&]() {
    if (std::isfinite(static_cast<double>(GTV)) && GTV >= 0.0f &&
        GTV != INF) {
      return GTV;
    }
    if (!ll_pt || !xomig) {
      return 0.0f;
    }

    double sum_w = 0.0;
    double stv = 0.0;
    auto pt_var = ll_pt->resetpt();
    while (ll_pt->noend(pt_var)) {
      auto punt = ll_pt->llpt(pt_var);
      double w = static_cast<double>(*(punt - 1));
      sum_w += w;
      for (int j = 0; j < Dim; j++) {
        double diff = static_cast<double>(punt[j]) -
                      static_cast<double>(xomig[j]);
        stv += diff * diff * w;
      }
      ll_pt->advpt(&pt_var);
    }
    if (sum_w <= 0.0) {
      return 0.0f;
    }
    double variance = stv / sum_w;
    if (!std::isfinite(variance) || variance < 0.0) {
      return 0.0f;
    }

    GTV = static_cast<float>(variance);
    return GTV;
  };
  auto write_singleton_row = [&]() {
    if (result == NULL) {
      Rcpp::stop("result is null in espai::obtenir data.\n");
    }
    ensure_row_capacity();
    *(result++) = 0;
    *(result++) = 0;
    *(result++) = 1;
    *(result++) = 1;
    *(result++) = singleton_variance();
    for (i = 0; i < Dim; i++)
      *(result++) = xomig ? xomig[i] : 0.0f;
    for (i = 0; i < Dim; i++)
      *(result++) = (Dim == 1 && i == 0) ? 1.0f : 0.0f;
    (*nrow) = 1;
  };
  if (ll_pop == NULL) {
    write_singleton_row();
    return;
  }

  auto pt = ll_pop->resetpt();
  if (pt == NULL) {
    Rcpp::stop("pt is null in espai::obtenir data.\n");
  }
  if (ll_pop->llpt(pt) == NULL) {
    write_singleton_row();
    return;
  }
  while (ll_pop->noend(pt)) {
    if (result == NULL) {
      Rcpp::stop("result is null in espai::obtenir data.\n");
    }
    if (pt == NULL) {
      Rcpp::stop("pt is null in espai::obtenir data.\n");
    }
    auto cur = (pop *)ll_pop->llpt(pt);
    if (!cur) {
      return;
    }
    ensure_row_capacity();
    *(result++) = 0;
    if (result == NULL) {
      Rcpp::stop("result is null in espai::obtenir data.\n");
    }
    *(result++) = cur->I;
    *(result++) = cur->density;
    *(result++) = cur->span;
    *(result++) = cur->var_k;

    for (i = 0; i < Dim; i++)
      *(result++) = cur->alpha[i];
    for (i = 0; i < Dim; i++)
      *(result++) = cur->b_ast[i];
    (*nrow)++;

    sespai = cur->espai_;
    auto sll_pop = sespai ? sespai->ll_pop : NULL;

    // subspace //
    if (PROF_REQ == 2 && sll_pop) {
      auto pt2 = sll_pop->resetpt();
      while (sll_pop->noend(pt2)) {
        ensure_row_capacity();
        *(result++) = 1;
        *(result++) = ((pop *)sll_pop->llpt(pt2))->I;
        *(result++) = ((pop *)sll_pop->llpt(pt2))->density;
        *(result++) = ((pop *)sll_pop->llpt(pt2))->span;
        *(result++) = ((pop *)sll_pop->llpt(pt2))->var_k;

        auxb = sespai->Ma->aplicar_Ma_vect(((pop *)sll_pop->llpt(pt2))->b_ast);
        auxa = sespai->Ma->aplicar_Ma_punt(((pop *)sll_pop->llpt(pt2))->alpha);
        for (i = 0; i < Dim; i++)
          *(result++) = auxa[i];
        for (i = 0; i < Dim; i++)
          *(result++) = auxb[i];

        (*nrow)++;
        delete[] auxb;
        delete[] auxa;

        sll_pop->advpt(&pt2);
      }
      ensure_row_capacity();
      *(result++) = 1;
      *(result++) = ((pop *)sll_pop->llpt(pt2))->I;
      *(result++) = ((pop *)sll_pop->llpt(pt2))->density;
      *(result++) = ((pop *)sll_pop->llpt(pt2))->span;
      *(result++) = ((pop *)sll_pop->llpt(pt2))->var_k;

      auxb = sespai->Ma->aplicar_Ma_vect(((pop *)sll_pop->llpt(pt2))->b_ast);
      auxa = sespai->Ma->aplicar_Ma_punt(((pop *)sll_pop->llpt(pt2))->alpha);
      for (i = 0; i < Dim; i++)
        *(result++) = auxa[i];
      for (i = 0; i < Dim; i++)
        *(result++) = auxb[i];
      (*nrow)++;
      delete[] auxb;
      delete[] auxa;
    }
    // end subspace

    ll_pop->advpt(&pt);
  }
  auto last = (pop *)ll_pop->llpt(pt);
  if (!last) {
    return;
  }
  ensure_row_capacity();
  *(result++) = 0;
  *(result++) = last->I;
  *(result++) = last->density;
  *(result++) = last->span;
  *(result++) = last->var_k;

  for (i = 0; i < Dim; i++)
    *(result++) = last->alpha[i];
  for (i = 0; i < Dim; i++)
    *(result++) = last->b_ast[i];

  (*nrow)++;
}

// Modified 2002-04-16: initialize static tuning parameters.
int espai ::PROF_REQ;
int espai ::NPARTS;
float espai ::C_H;
float espai ::C_D;

void espai::inicializar_nparts_ch_cd(int profreq, int nparts, float c_h,
                                     float c_d) {
  if (profreq <= 0)
    PROF_REQ = 1;
  else
    PROF_REQ = profreq;
  if (nparts < 3 || nparts > 6)
    NPARTS = 4;
  else
    NPARTS = nparts;
  if (c_h < 0.5 || c_h > 1.5)
    C_H = 0.75;
  else
    C_H = c_h;
  if (c_d < 0.25 || c_d > 0.5)
    C_D = 0.3; // always less than 0.5
  else
    C_D = c_d;
}

/*
 * Eigenvalues and eigenvectors of a real symmetric matrix
 *
 * SYNOPSIS:
 *
 * int n;
 * double A[n*(n+1)/2], EV[n*n], E[n];
 * void eigens(A, EV, E, n);
 *
 * DESCRIPTION:
 *
 * The algorithm is due to J. vonNeumann.
 *
 * A[] is a symmetric matrix stored in lower triangular form.
 * That is, A[row, column] = A[(row*row+row)/2 + column]
 * or equivalently with row and column interchanged. The
 * indices row and column run from 0 through n-1.
 *
 * EV[] is the output matrix of eigenvectors stored columnwise.
 * That is, the elements of each eigenvector appear in sequential
 * memory order. The jth element of the ith eigenvector is
 * EV[n*i+j] = EV[i][j].
 *
 * E[] is the output matrix of eigenvalues. The ith element
 * of E corresponds to the ith eigenvector (the ith row of EV).
 *
 * On output, the matrix A will have been diagonalized and its
 * original contents are destroyed.
 *
 * ACCURACY:
 *
 * The error is controlled by an internal parameter called RANGE
 * which is set to 1e-10. After diagonalization, the
 * off-diagonal elements of A will have been reduced by
 * this factor.
 *
 * ERROR MESSAGES:
 *
 * None.
 */
/* Copyright 1973, 1991 by Stephen L. Moshier Copyleft version. */

void espai::eigens(float *A, float *RR, float *E, int N)
// double A[], RR[], E[];
// int N;
{
  int IND, L, LL, LM, M, MM, MQ, I, J, IA, LQ;
  int IQ, IM, IL, NLI, NMI;
  double ANORM, ANORMX, AIA, THR, ALM, ALL, AMM, Xv, Y;
  double SINX, SINX2, COSX, COSX2, SINCS, AIL, AIM;
  double RLI, RMI;
  // double sqrt(), fabs();
  static double RANGE = 1.0e-10; /*3.0517578e-5;*/

  /* Initialize identity matrix in RR[] */
  // for( J=0; J<N*N; J++ ) // is already initialized.
  //RR[J] = 0.0;

  MM = 0;
  for (J = 0; J < N; J++) {
    RR[MM + J] = 1.0;
    MM += N;
  }

  ANORM = 0.0;
  for (I = 0; I < N; I++) {
    for (J = 0; J < N; J++) {
      if (I != J) {
        IA = I + (J * J + J) / 2;
        AIA = A[IA];
        ANORM += AIA * AIA;
      }
    }
  }
  if (ANORM <= 0.0)
    goto done;
  ANORM = sqrt(ANORM + ANORM);
  ANORMX = ANORM * RANGE / N;
  THR = ANORM;

  while (THR > ANORMX) {
    THR = THR / N;

    do { /* while IND != 0 */
      IND = 0;

      for (L = 0; L < N - 1; L++) {

        for (M = L + 1; M < N; M++) {
          MQ = (M * M + M) / 2;
          LM = L + MQ;
          ALM = A[LM];
          if (fabs(ALM) < THR)
            continue;

          IND = 1;
          LQ = (L * L + L) / 2;
          LL = L + LQ;
          MM = M + MQ;
          ALL = A[LL];
          AMM = A[MM];
          Xv = (ALL - AMM) / 2.0;
          Y = -ALM / sqrt(ALM * ALM + Xv * Xv);
          if (Xv < 0.0)
            Y = -Y;
          SINX = Y / sqrt(2.0 * (1.0 + sqrt(1.0 - Y * Y)));
          SINX2 = SINX * SINX;
          COSX = sqrt(1.0 - SINX2);
          COSX2 = COSX * COSX;
          SINCS = SINX * COSX;

          /* ROTATE L AND M COLUMNS */
          for (I = 0; I < N; I++) {
            IQ = (I * I + I) / 2;
            if ((I != M) && (I != L)) // we are not interested in the eigenvalues,
                                      // only the eigenvectors
            {
              if (I > M)
                IM = M + IQ;
              else
                IM = I + MQ;
              if (I >= L)
                IL = L + IQ;
              else
                IL = I + LQ;
              AIL = A[IL];
              AIM = A[IM];
              Xv = AIL * COSX - AIM * SINX;
              A[IM] = AIL * SINX + AIM * COSX;
              A[IL] = Xv;
            }
            NLI = N * L + I;
            NMI = N * M + I;
            RLI = RR[NLI];
            RMI = RR[NMI];
            RR[NLI] = RLI * COSX - RMI * SINX;
            RR[NMI] = RLI * SINX + RMI * COSX;
          }

          Xv = 2.0 * ALM * SINCS; // we are not interested in the eigenvalues, only
                                  // the eigenvectors
          A[LL] = ALL * COSX2 + AMM * SINX2 - Xv;
          A[MM] = ALL * SINX2 + AMM * COSX2 + Xv;
          A[LM] = (ALL - AMM) * SINCS + ALM * (COSX2 - SINX2);
        } /* for M=L+1 to N-1 */
      } /* for L=0 to N-2 */

    } while (IND != 0);

  } /* while THR > ANORMX */

done:;

  /* Extract eigenvalues from the reduced matrix */
  L = 0;
  for (J = 1; J <= N; J++) // we are not interested in the eigenvalues.
  {
    L = L + J;
    E[J - 1] = A[L - 1];
  }
}

// vector ops /////////////////

float espai::distancia(float *pnt1, float *pnt2) {
  int i;
  float sum = 0.0;
  for (i = 0; i < Dim; i++)
    sum += pow(pnt1[i] - pnt2[i], 2);
  return sqrt(sum);
}

float *espai::mult_esc(float e, float *v) {
  int i;
  float *v3;

  v3 = new float[Dim];
  for (i = 0; i < Dim; i++)
    v3[i] = v[i] * e;
  return v3;
}

float espai::mult_v(float *v1, float *v2) {
  int i;
  float sum = 0.0;
  for (i = 0; i < Dim; i++) {
    sum += v1[i] * v2[i];
  }
  return sum;
}

float *espai::sum_v(float *v1, float *v2) {
  int i;
  float *v3;

  v3 = new float[Dim];
  for (i = 0; i < Dim; i++)
    v3[i] = v1[i] + v2[i];
  return v3;
}

float *espai::dif_v(float *v2, float *v1) {
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

float *espai::norma_v(float *v) {
  // returns the normalized vector
  int i;
  float nrm = 0.0;
  float *v3;

  v3 = new float[Dim];
  for (i = 0; i < Dim; i++)
    nrm += pow(v[i], 2);
  nrm = sqrt(nrm);
  for (i = 0; i < Dim; i++)
    v3[i] = v[i] / nrm;
  return v3;
}
