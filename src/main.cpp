extern "C" {
#include <stdlib.h>
}
#include "espai.h"
#include <Rcpp.h>

//' @title pcop_backend
//' @description Computes a principal curve as defined in Delicado (2001). DO NOT use this function unless you know what you are doing. Use `pcop()` instead.
//' <doi: 10.1007/s001800300145>
//' @param x    See `pcop()`
//' @param c_h  See `pcop()`
//' @param c_d  See `pcop()`
//' @return A numeric matrix to be parsed by `pcop()`.
// [[Rcpp::export]]
Rcpp::NumericMatrix pcop_backend(const Rcpp::NumericMatrix &x, float c_d, float c_h){
	int profreq = 1;
	int nparts = 4;
	espai *psp;
	float **Ma;
	int i;
	ll_p *ll_pt;
	float vtg;
	float *mx;

	// inicialización valores por defecto, despues en el fichero de setup puede que cambien
	// modificado 16/4/2002
	//PROFREQ =1;
	//NPARTs = 4;
	//C_h = 0.75;
	//C_d = 0.4;		// siempre menor que 0.5

	int Dim = x.ncol();
	ll_pt = new ll_p(Dim);
	for (i = 0; i < x.nrow(); i++) {
		float *d_punt = (float*)malloc((Dim + 1) * sizeof(float));
		d_punt[0] = 1; d_punt++;// coordenada -1 pel pes.
		for (int j = 0; j < x.ncol(); j++) {
			d_punt[j] = x(i, j);
		}
		ll_pt->add_ordX_principal(d_punt); //###
	}
	if (ll_pt->n_punts()< NPTMIN*Dim) {
		Rcpp::stop("Warning: Not enough points in data matrix. At least %d points are needed for dimension %d.\n", NPTMIN * Dim, Dim);
	}
	Ma = (float **) malloc(Dim*sizeof(float *));
	if (Ma == NULL) {
		Rcpp::stop("Could not allocate Ma.\n");
	}

	for (i=0;i<Dim;i++){
		Ma[i] = (float *)calloc(Dim,sizeof(float));
		Ma[i][i]=1;
	}
	mx = (float *) calloc(Dim,sizeof(float));
	if (mx == NULL) {
		Rcpp::stop("Could not allocate mx.\n");
	}
	psp = new espai(ll_pt,Dim,0);
	//modificacion 16/4/2002
	psp->inicializar_nparts_ch_cd(profreq,nparts,c_h,c_d);
	//fin
	psp->rebre_M_a(new M_a(Dim,0,Ma,mx));
	vtg = psp->obtenir_VTG(&mx);
	int nrow, ncol;
	float *out = (float*)malloc((2 * Dim + 5) * x.nrow() * sizeof(float));
	if (out == NULL) {
		Rcpp::stop("Could not allocate out.\n");
	}
	psp->obtenir_data(out, &nrow, &ncol);

	delete psp;
	Rcpp::NumericMatrix proy(nrow, ncol, out);
	proy = transpose(proy);
	return proy;

}
