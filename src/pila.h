class pila {

private:
  typedef struct node {
    void *pt;
    node *seg;
  } node;
  node *top;
  node *n_top;

public:
  // constructor
  pila();
  ~pila();
  void apilar(void *pt);
  void *desapilar();

  // query method
  int pila_buida();
};
