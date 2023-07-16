#include <cmath>
#include <cstdlib>

template <typename T> class ll_pnt {

private:
  typedef struct node {
    T *info;
    node *seg;
  } node;

  node *Topleft;
  node *newTopleft;
  node *Topright;

public:
  ll_pnt() {
    Topleft = new node();
    Topright = Topleft;
  }

  ~ll_pnt() {
    auto act_pt = Topleft;
    while (act_pt) {
      auto aux_pt = act_pt->seg;
      delete act_pt->info;
      delete act_pt;
      act_pt = aux_pt;
    }
  }

  void add(T *info) {
    Topright->info = info;
    Topright->seg = new node();
    Topright = Topright->seg;
  }

  void addrev(T *info) {
    newTopleft = new node;
    newTopleft->info = info;
    newTopleft->seg = Topleft;
    Topleft = newTopleft;
  }

  node* resetpt() { return Topleft; }

  node *noend(node *pt) {
    return (pt->seg->seg); // a l'ultim elem. contindr info, s'haur
                                     // de tractar apart.
    // si pt->seg == NULL el camp info de pt estara buit
  }

  T *llpt(node *pt) { return (pt->info); }

  void advpt(node **pt) { *pt = (*pt)->seg; }

  void modpt(node *pt, T *info) { pt->info = info; }
};
