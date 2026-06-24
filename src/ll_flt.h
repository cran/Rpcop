struct node {
  float info;
  node *seg;
};

class ll_flt {
private:
  node *Topleft;
  node *Topright;

public:
  ll_flt() {
    Topleft = new node();
    Topleft->seg = nullptr;
    Topright = Topleft;
  }

  ~ll_flt() {
    node *act_pt, *aux_pt;

    act_pt = Topleft;
    while (act_pt) {
      aux_pt = act_pt->seg;
      delete act_pt;
      act_pt = aux_pt;
    }
  }

  void add(float info) {
    Topright->info = info;
    Topright->seg = new node();
    Topright = Topright->seg;
  }

  node *resetpt() { return Topleft; }

  node *noend(node *pt) {
    return pt->seg; // the final sentinel has seg == NULL
                    // and does not carry a valid info value
  }

  float llpt(node *pt) { return pt->info; }

  void advpt(node **pt) { *pt = (*pt)->seg; }

  void modpt(node *pt, float info) { (pt)->info = info; }
};
