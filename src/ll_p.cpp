#include "ll_p.h"

ll_p::ll_p(int d) {
  int i;
  Dim = d;
  orcluster = ESQUERRA;
  numcl = 0;
  sum_w = 0; // sum of point weights
  suma_d = 0;
  xorig = NULL;
  xoant = NULL;
  semilla = NULL;
  candidat = NULL;

  topleft = new node();
  topleft->coord = new float[Dim + 1]();
  topleft->coord++;
  topleft->coord[X] = -1 * INF;
  topright = new node();
  topright->coord = new float[Dim + 1]();
  topright->coord++; // Dim +1 to include the weights
  topright->coord[X] = INF;

  topleft->seg[DRETA] = topright;
  topleft->noin[DRETA] = (node *)topright;
  topright->seg[ESQUERRA] = topleft;
  topright->noin[ESQUERRA] = (node *)topleft;

  topright->seg[DRETA] = NULL;
  topright->noin[DRETA] = NULL;
  topleft->seg[ESQUERRA] = NULL;
  topleft->noin[ESQUERRA] = NULL;

  topright->marca = -1; // marker for nodes already inserted in the spanning tree
  topleft->marca = -1;
  vn_punts = 0;
  max = new float[Dim];
  for (i = 0; i < Dim; i++)
    max[i] = -1 * INF;
  min = new float[Dim];
  for (i = 0; i < Dim; i++)
    min[i] = INF;
  x_mean = new float[Dim + 1](); // for depth-0 spaces this accumulates xmean;
                                 // otherwise it is the local origin passed down.
  owns_x_mean = true;
}

ll_p::~ll_p() {
  node *pt = topleft;
  node *auxpt;
  node_satelit *ptst;
  node_satelit *auxptst;
  if (pt->seg[DRETA] == pt->noin[DRETA])
    while (pt) {
      auxpt = pt;
      pt = pt->seg[DRETA];
      delete[] (
          (auxpt->coord) -
          1); // adjust pointer offset used when descending to Dim-1
      delete auxpt;
    }
  else
    while (pt) {
      ptst = (node_satelit *)pt->noin[DRETA];
      while (ptst) {
        auxptst = ptst;
        ptst = ptst->seg;
        delete auxptst;
      }
      ptst = (node_satelit *)pt->noin[ESQUERRA];
      while (ptst) {
        auxptst = ptst;
        ptst = ptst->seg;
        delete auxptst;
      }
      auxpt = pt;
      pt = pt->seg[DRETA];
      delete[] (
          (auxpt->coord) -
          1); // adjust pointer offset used when descending to Dim-1
      delete auxpt;
    }

  delete[] max;
  delete[] min;
  if (owns_x_mean) {
    delete[] x_mean;
  }
}

void ll_p::add_ordX_principal(float *vect) {
  node *act_node, *new_node;
  float *x_mean02, *x_mean03;

  int i;

  act_node = topright->seg[ESQUERRA];

  if (vect[X] >
      topleft->seg[DRETA]->coord[X] +
          (0.5 *
           (act_node->coord[X] -
            topleft->seg[DRETA]->coord[X]))) { // if greater than midpoint
                                               // start from the right side

    while (vect[X] < act_node->coord[X]) {
      act_node = act_node->seg[ESQUERRA];
    }
    new_node = new node;
    new_node->coord = vect;
    new_node->marca = 0;
    new_node->seg[ESQUERRA] = act_node;
    new_node->noin[ESQUERRA] = (node *)act_node;
    new_node->seg[DRETA] = act_node->seg[DRETA];
    new_node->noin[DRETA] = (node *)act_node->noin[DRETA];
    act_node->seg[DRETA]->seg[ESQUERRA] = new_node;
    act_node->seg[DRETA]->noin[ESQUERRA] = (node *)new_node;
    act_node->seg[DRETA] = new_node;
    act_node->noin[DRETA] = (node *)new_node;
    vn_punts++;
  } else {
    act_node = topleft->seg[DRETA];
    while (vect[X] > act_node->coord[X]) {
      act_node = act_node->seg[DRETA];
    }
    new_node = new node;
    new_node->coord = vect;
    new_node->marca = 0;
    new_node->seg[DRETA] = act_node;
    new_node->noin[DRETA] = (node *)act_node;
    new_node->seg[ESQUERRA] = act_node->seg[ESQUERRA];
    new_node->noin[ESQUERRA] = (node *)act_node->noin[ESQUERRA];
    act_node->seg[ESQUERRA]->seg[DRETA] = new_node;
    act_node->seg[ESQUERRA]->noin[DRETA] = (node *)new_node;
    act_node->seg[ESQUERRA] = new_node;
    act_node->noin[ESQUERRA] = (node *)new_node;
    vn_punts++;
  }

  /* update max, min, xmean */ // since ll_p represents one full cluster, these are
                               // the cluster-wide extrema.
  for (i = 0; i < Dim; i++) {
    if (max[i] < vect[i])
      max[i] = vect[i];
    if (min[i] > vect[i])
      min[i] = vect[i];
  }
  /* calculate xmean */
  sum_w += *(vect - 1);
  x_mean02 = mult_esc(*(vect - 1), vect);
  x_mean03 = sum_v(x_mean, x_mean02);
  delete[] x_mean;
  delete[] x_mean02;
  x_mean = x_mean03;
}

void ll_p::donar_max_min_xomig(float **mx, float **mn, float **xm, float *s_d,
                               bool *owns_xm) {

  *mx = max;
  *mn = min;
  *xm = x_mean; // weighted xmean (or final-space representative origin).
  *owns_xm = owns_x_mean;
  *s_d = suma_d;
  max = NULL;
  min = NULL;
  x_mean = NULL;
  owns_x_mean = false;
}

void ll_p::calcular_max_min_cluster() {
  node_satelit *s_act;
  int i;

  numcl++; // keep different from sentinel mark values

  for (i = 0; i < Dim; i++) {
    if (max[i] < xorig->coord[i])
      max[i] = xorig->coord[i];
    if (min[i] > xorig->coord[i])
      min[i] = xorig->coord[i];
  }

  for (orcluster = ESQUERRA; orcluster <= DRETA; orcluster++) {
    s_act = (node_satelit *)xorig->noin[orcluster];

    for (i = 0; i < Dim; i++) {
      if (max[i] < s_act->ptnode->coord[i])
        max[i] = s_act->ptnode->coord[i];
      if (min[i] > s_act->ptnode->coord[i])
        min[i] = s_act->ptnode->coord[i];
    }
    s_act->ptnode->marca = numcl;
    p_n.apilar(s_act->ptnode);

    while (s_act->seg) {
      s_act = s_act->seg;
      if (s_act->ptnode->marca != numcl) {
        for (i = 0; i < Dim; i++) {
          if (max[i] < s_act->ptnode->coord[i])
            max[i] = s_act->ptnode->coord[i];
          if (min[i] > s_act->ptnode->coord[i])
            min[i] = s_act->ptnode->coord[i];
        }
        s_act->ptnode->marca = numcl;
        p_n.apilar(s_act->ptnode);
      }
    }

    while (!p_n.pila_buida()) {
      while (
          !p_n.pila_buida() &&
          !(s_act =
                ((node_satelit *)((node *)p_n.desapilar())->noin[orcluster])))
        ; // we have to consider points situated at the limit of the cluster on the
          // X coordinate
      if (s_act) {
        if (s_act->ptnode->marca != numcl) {
          for (i = 0; i < Dim; i++) {
            if (max[i] < s_act->ptnode->coord[i])
              max[i] = s_act->ptnode->coord[i];
            if (min[i] > s_act->ptnode->coord[i])
              min[i] = s_act->ptnode->coord[i];
          }
          s_act->ptnode->marca = numcl;
          p_n.apilar(s_act->ptnode);
        }
        while (s_act->seg) {
          s_act = s_act->seg;
          if (s_act->ptnode->marca != numcl) {
            for (i = 0; i < Dim; i++) {
              if (max[i] < s_act->ptnode->coord[i])
                max[i] = s_act->ptnode->coord[i];
              if (min[i] > s_act->ptnode->coord[i])
                min[i] = s_act->ptnode->coord[i];
            }
            s_act->ptnode->marca = numcl;
            p_n.apilar(s_act->ptnode);
          }
        }
      }
    }
  }
}

float ll_p::inicialitzacio_principal() {
  ll_q *quartiles;
  float *auxx_mean;

  auxx_mean = x_mean; // x_mean theoretical.
  x_mean = mult_esc(1.0 / sum_w, auxx_mean);
  delete[] auxx_mean;
  quartiles = new ll_q(vn_punts); // number of points in this space
  obtener_quartiles(quartiles);   // compute MST distances to get quartiles
  dmax = quartiles->dmax();
  x_mean = obtener_satelites(); // use dmax to obtain satellites
                                // and locate xorig.
  delete quartiles;
  return dmax; // we only need dmax to separate the different curves
}

void ll_p::inicialitzacio_final() {
  /* No MST/satellites here: final spaces do not compute additional curves. */
  float *auxx_mean;

  auxx_mean = x_mean; // x_mean theoretical.
  x_mean = mult_esc(
      1.0 / sum_w,
      auxx_mean); // normalize xmean by total point weight.
  delete[] auxx_mean;
}

void ll_p::mstinsertar(node *pt) {
  pt->marca = -1; // will be reset to 0 in calcular_satellites() candidate traversal.

  ((node *)pt->noin[ESQUERRA])->noin[DRETA] = ((node *)pt->noin[DRETA]);
  ((node *)pt->noin[DRETA])->noin[ESQUERRA] = ((node *)pt->noin[ESQUERRA]);

  // clear links used by satellite traversal

  pt->noin[DRETA] = NULL;
  pt->noin[ESQUERRA] = NULL;
}

int ll_p::mstinsertat(node *pt) { return pt->marca; }

void ll_p::obtener_quartiles(
    ll_q *llqt) { // compute the minimum spanning tree and store distances
                  // used to obtain quartiles
  int aux_npunts = vn_punts - 1;
  node *xact;
  node *xpost;
  node *xpt;
  node *nxtnomstin;
  float dpost, dant, d; /* dant: min distance to a node of the tree */
  /* dpost: min distance to a node outside the tree */
  xpt = topleft->seg[DRETA];

  /* preprocessing step */ // locate the second point to insert in the MST

  xact = xpt->seg[DRETA];
  xpost =
      (node *)xpt->noin[DRETA]; /* for the case: Xright = dpost> Xleft */
  dpost = distancia(xpost->coord,
                    xpt->coord); /* when xact=topright, dpost = Xact-Xpt */

  while (xact->coord[X] - xpt->coord[X] < dpost) {
    if ((d = distancia(xact->coord, xpt->coord)) <= dpost)
      dpost = d;
    xact = xact->seg[DRETA]; // all points are initially outside the MST
  }
  mstinsertar(xpt); // insert previous xpost into MST, without unlinking it
                    // from the `noin` list
  xpt = xpost;      // xpost is the next point to insert
  dant = dpost;
  aux_npunts--; // we subtract one too many because the last insertion in ll_q is
                // performed in the post-process

  /* main body */

  while (aux_npunts) {

    xact = xpt->seg[DRETA];
    nxtnomstin = (node *)xpt->noin[DRETA];
    xpost = nxtnomstin; /* for the case: Xright = dpost> Xleft */
    dpost = distancia(xpost->coord,
                      xpt->coord); /* when xact=topright, dpost = Xact-Xpt */

    /* right-side scan */

    while (xact->coord[X] - xpt->coord[X] < dant &&
           xact->coord[X] - xpt->coord[X] < dpost) {
      if (!mstinsertat(xact)) {
        if ((d = distancia(xact->coord, xpt->coord)) < dpost) {
          dpost = d;
          xpost = xact;
        }
        nxtnomstin = (node *)xact->noin[DRETA];
      } else if ((d = distancia(xact->coord, xpt->coord)) < dant)
        dant = d;
      xact = xact->seg[DRETA];
    }

    if (xact->coord[X] - xpt->coord[X] >= dant) { // dant bound reached
      xact = nxtnomstin; // continue from last non-MST candidate
      while (xact->coord[X] - xpt->coord[X] < dpost &&
             xact != topright) { // if dpost is also bounded, this exits early
        if ((d = distancia(xact->coord, xpt->coord)) <= dpost) {
          dpost = d;
          xpost = xact;
        }
        xact = (node *)xact->noin[DRETA];
      }
    } else { // dpost bound reached
      while (xact->coord[X] - xpt->coord[X] < dant) {
        if ((mstinsertat(xact)) &&
            ((d = distancia(xact->coord, xpt->coord)) < dant))
          dant = d;
        xact = xact->seg[DRETA];
      }
    }

    /* left-side scan */
    xact = xpt->seg[ESQUERRA];
    nxtnomstin = (node *)xpt->noin[ESQUERRA];
    while (xpt->coord[X] - xact->coord[X] < dant &&
           xpt->coord[X] - xact->coord[X] < dpost) {
      if (!mstinsertat(xact)) {
        if ((d = distancia(xact->coord, xpt->coord)) < dpost) {
          dpost = d;
          xpost = xact;
        }
        nxtnomstin = (node *)xact->noin[ESQUERRA];
      } else if ((d = distancia(xact->coord, xpt->coord)) < dant)
        dant = d;
      xact = xact->seg[ESQUERRA];
    }

    if (xpt->coord[X] - xact->coord[X] >= dant) { // dant bound reached
      xact = nxtnomstin; // continue from last non-MST candidate
      while (
          xpt->coord[X] - xact->coord[X] <
          dpost) { // if dpost is also bounded, this exits early
        if ((d = distancia(xact->coord, xpt->coord)) < dpost) {
          dpost = d;
          xpost = xact;
        }
        xact = (node *)xact->noin[ESQUERRA];
      }
    } else { // dpost bound reached
      while (xpt->coord[X] - xact->coord[X] < dant) {
        if ((mstinsertat(xact)) &&
            ((d = distancia(xact->coord, xpt->coord)) < dant))
          dant = d;
        xact = xact->seg[ESQUERRA];
      }
    }

    /* final-iteration treatment */

    llqt->add_ord(dant);
    suma_d += dant;
    mstinsertar(xpt); // insert previous xpost into MST, but do not unlink it
                      // from the `noin` list
    xpt = xpost;
    dant = dpost;
    aux_npunts--;
  }

  /* post-processing */ // verify there is no smaller dant for the final xpost
                        // to insert into the MST

  xact = xpt->seg[DRETA];
  while (xact->coord[X] - xpt->coord[X] < dant) {
    if ((d = distancia(xact->coord, xpt->coord)) <
        dant) // only points already inserted remain here
      dant = d;
    xact = xact->seg[DRETA];
  }

  xact = xpt->seg[ESQUERRA];
  while (xpt->coord[X] - xact->coord[X] < dant) {
    if ((d = distancia(xact->coord, xpt->coord)) <
        dant) // only points already inserted remain here
      dant = d;
    xact = xact->seg[ESQUERRA];
  }

  llqt->add_ord(dant);
  suma_d += dant;
  // clear links used by satellite traversal
  xpt->noin[DRETA] = NULL;
  xpt->noin[ESQUERRA] = NULL;

  topleft->noin[DRETA] = NULL;
  topright->noin[ESQUERRA] = NULL;
}

float *ll_p::obtener_satelites() {
  float d, dpost = INF;
  node *xpt = topleft->seg[DRETA];
  node *xact = xpt->seg[DRETA];
  node *fallback = (xpt && xpt != topright) ? xpt : NULL;
  xorig = NULL;

  while (xpt->seg[DRETA]) {
    while (xact->coord[X] - xpt->coord[X] < dmax) {
      if (distancia(xact->coord, xpt->coord) < dmax) {
        add_satelit(DRETA, xpt, xact);
        add_satelit(ESQUERRA, xact, xpt);
        if ((xpt->noin[ESQUERRA]) &&
            ((d = distancia(xpt->coord, x_mean)) < dpost)) {
          dpost = d;
          xorig = xpt;
        }; // search xorig. It must have satellites on both left and right.
      };
      xact = xact->seg[DRETA];
    }

    xpt = xpt->seg[DRETA];
    xact = xpt->seg[DRETA];
  }
  if (!xorig || xorig == topright) {
    xorig = fallback;
  }
  xoant = xorig;

  if (xoant && xoant->coord) {
    delete[] x_mean; // replace theoretical xmean with practical xorig
    owns_x_mean = false;
    return xoant->coord;
  }
  return x_mean;
}

void ll_p::add_satelit(int or_, node *ptor, node *ptdsti) {
  node_satelit *new_nst;

  new_nst = new node_satelit;
  new_nst->ptnode = ptdsti;
  new_nst->seg = (node_satelit *)ptor->noin[or_]; // are prepended
  ptor->noin[or_] = new_nst;
}

void ll_p::tornar_a_xomig() { xoant = xorig; }

void ll_p::trobar_primer_candidat_clt(float *xm) {
  float d, dtop;
  node *xtop;
  node *xodmax;

  if (!xm || !xoant || !xoant->coord) {
    return;
  }

  orcluster = (xm[X] > xoant->coord[X]); // true ->right, false ->left
  while (xoant->seg[orcluster] && fabs(xm[X] - xoant->coord[X]) > dmax)
    xoant = xoant->seg[orcluster];
  if (!xoant || !xoant->coord) {
    return;
  }
  xodmax = xoant;

  while (xoant->seg[orcluster] && xoant->marca < 1)
    xoant = xoant->seg[orcluster]; // we do not require connectivity on both
                                   // sides; one side is enough
  if (!xoant || !xoant->coord) {
    return;
  }
  dtop = distancia(xm, xoant->coord);
  xtop = xoant;
  while (xoant->seg[orcluster] && fabs(xoant->coord[X] - xm[X]) < dtop) {
    xoant = xoant->seg[orcluster];
    if (!xoant || !xoant->coord) {
      break;
    }
    if (((d = distancia(xm, xoant->coord)) < dtop) &&
        (xoant->marca > 0)) { // skip previously visited marks (xoant->marca>0)
      dtop = d;
      xtop = xoant;
    }
  }
  if (dtop > dmax) { // if both candidate sides are farther than dmax from xm
    xoant = xodmax;
    orcluster = (orcluster + 1) % 2;
    while (xoant && xoant->coord && xoant->seg[orcluster] &&
           fabs(xoant->coord[X] - xm[X]) < dtop) {
      xoant = xoant->seg[orcluster];
      if (!xoant || !xoant->coord) {
        break;
      }
      if (((d = distancia(xm, xoant->coord)) < dtop) &&
          (xoant->marca > 0)) { // skip previously visited marks (xoant->marca>0)
        dtop = d;
        xtop = xoant;
      }
    }
  }
  if (xtop) {
    xoant = xtop;
  }
}

float *ll_p::primer_candidat_clt() {
  numcl++;
  semilla = NULL;
  candidat = NULL;

  if (!xoant) {
    return NULL;
  }

  if (!(candidat = (node_satelit *)xoant->noin[orcluster])) {
    orcluster = (orcluster + 1) % 2;
    candidat = (node_satelit *)
                   xoant->noin[orcluster]; // it should have at least one satellite;
                                           // otherwise it would not belong to the cluster.
  }
  if (!candidat || !candidat->ptnode || !candidat->ptnode->coord) {
    return NULL;
  }
  semilla = xoant;
  return candidat->ptnode->coord; // ### known issue: first satellite may be skipped
                                  // while xoant may be inserted twice
}

float *ll_p::seguent_candidat_clt(int validacio) {

  candidat->ptnode->marca = numcl;
  if (validacio) // The space object validated the cluster candidate
    p_n.apilar(candidat->ptnode);
  /* Only search points inside the candidate cluster along X in orcluster direction. */
  do {
    if (candidat->seg)
      candidat = candidat->seg;
    else if (semilla) { // We can come from trying the 2nd dir. or the 1st.
      orcluster = (orcluster + 1) % 2;
      if (!(candidat = (node_satelit *)((node *)semilla)->noin[orcluster])) {
        do {
          if (p_n.pila_buida())
            return NULL;
          else
            semilla = (node *)p_n.desapilar();
          if (!(candidat =
                    (node_satelit *)((node *)semilla)->noin[orcluster])) {
            orcluster = (orcluster + 1) % 2;
            candidat = (node_satelit *)((node *)semilla)->noin[orcluster];
            semilla = NULL;
          }
        } while (!candidat); // if first-direction candidate is valid,
                             // candidate may be evaluated twice.
      } else
        semilla = NULL; // second direction still pending
    } else
      do {
        if (p_n.pila_buida())
          return NULL;
        else
          semilla = (node *)p_n.desapilar();
        if (!(candidat = (node_satelit *)((node *)semilla)->noin[orcluster])) {
          orcluster = (orcluster + 1) % 2;
          candidat = (node_satelit *)((node *)semilla)->noin[orcluster];
          semilla = NULL;
        }
      } while (!candidat); // if first-direction candidate is valid,
                           // candidate may be evaluated twice.
  } while (candidat->ptnode->marca == numcl);

  return candidat->ptnode->coord;
}

float *ll_p::canviar_orientacio_clt() {
  orcluster = (orcluster + 1) % 2; // switch direction without changing xoant
  // candidate = (node_satellite *)xoant->noin[orcluster]; // we take the first
  // satellite in opposite direction to the 1st candidate
  // seed = xoant;
  // if (!candidate) return NULL; // the xoant may be at the end of the
  // cluster at the X coordinate
  // if (candidate->ptnode->marca == numcl) return NULL; // check both directions
  // via semilla; return candidate->ptnode->coord; since this is opposite
  // direction, prior-insertion check is unnecessary
  return NULL;
}

int ll_p::n_punts() { return vn_punts; }

// list access helpers
// needed to compute final STV at maximum depth; kept here with list internals.

void *ll_p::noend(node *pt) {
  return (((node *)pt)->seg[DRETA]); // if NULL, pt is the sentinel
}

float *ll_p::llpt(node *pt) { return pt->coord; }

void ll_p::advpt(node **pt) { *pt = (*pt)->seg[DRETA]; }

// operations needed by main algorithm

// vector ops

float ll_p::distancia(float *pnt1, float *pnt2) {
  int i;
  float sum = 0.0;
  for (i = 0; i < Dim; i++)
    sum += pow(pnt1[i] - pnt2[i], 2);
  return sqrt(sum);
}

float *ll_p::mult_esc(float e, float *v) {
  int i;
  float *v3;

  v3 = new float[Dim];
  for (i = 0; i < Dim; i++)
    v3[i] = v[i] * e;
  return v3;
}

float *ll_p::sum_v(float *v1, float *v2) {
  int i;
  float *v3;

  v3 = new float[Dim];
  for (i = 0; i < Dim; i++)
    v3[i] = v1[i] + v2[i];
  return v3;
}
