
extern "C" {
#include <stdlib.h>
#include <stdio.h>
}

#define CQ 3 // constante para definir dmax a prtir del rango intercuartilico

class ll_q{
	private :

		typedef struct node{
			float info;
			node  *seg;
		} node;

		int npunts;

		node *Topleft;

	public :

		// constructoras
		ll_q(int np);
		~ll_q();
		void add_ord(float info);

		// consultoras

		float dmax();

};
