extern "C" {
#include "stdlib.h"
}

#include "pila.h"

pila::pila() { top = NULL; }

pila::~pila() {
  while (top) {
    n_top = top->seg;
    delete top;
    top = n_top;
  }
}

void pila::apilar(void *pt) {
  n_top = new node;

  n_top->pt = pt;
  n_top->seg = top;
  top = n_top;
}

void *pila::desapilar() {
  void *tp =
      top->pt; // Caller is expected to check that the stack is non-empty.

  n_top = top->seg;
  delete top;
  top = n_top;

  return tp;
}

int pila::pila_buida() { return !top; }
