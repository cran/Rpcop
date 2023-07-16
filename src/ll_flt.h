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
        Topleft = new node;
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

    node* resetpt() { return Topleft; }

    node *noend(node *pt) {
        return pt->seg; // a l'ultim sera NULL perque creem amb calloc.
        // si pt->seg == NULL el camp info de pt estara buit
    }

    float llpt(node *pt) { return pt->info; }

    void advpt(node **pt) { *pt = (*pt)->seg; }

    void modpt(node *pt, float info) { (pt)->info = info; }
};
