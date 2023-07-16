#include "ll_q.h"
extern "C" {
#include <stdlib.h>
}

ll_q::ll_q(int np) {
  npunts = np;
  Topleft = new node(); // does not contain info
}

ll_q::~ll_q() {
  node *pt = Topleft;
  node *auxpt;

  while (pt != NULL) {
    auxpt = pt;
    pt = pt->seg;
    delete auxpt;
  }
}

void ll_q::add_ord(float info) {
  node *pt = Topleft;
  node *newpt;

  while (pt->seg && pt->seg->info > info) {
    pt = pt->seg;
  }
  newpt = new node;
  newpt->info = info;
  newpt->seg = pt->seg;
  pt->seg = newpt;
}

float ll_q::dmax() {
  int limite, i;
  node *nd;

  limite = npunts * 0.00; // a un dmax ms petit -> calculs ms rapids, pero un
                          // cluster ms restrictiu.

  nd = Topleft->seg;
  for (i = 1; i < limite; i++) {
    nd = nd->seg;
  }
  return (nd->info);
}
