
extern "C" {
#include <stdio.h>
#include <stdlib.h>
}

#define CQ 3 // constant to define dmax from the interquartile range

class ll_q {
private:
  typedef struct node {
    float info;
    node *seg;
  } node;

  int npunts;

  node *Topleft;

public:
  // constructors
  ll_q(int np);
  ~ll_q();
  void add_ord(float info);

  // accessors

  float dmax();
};
