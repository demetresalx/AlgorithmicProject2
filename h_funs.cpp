#include "h_funs.h"

/*
template <class T>
int our_mod(T a, T b){ //returns remainder as it should
	return (a % b + b) % b;
}
*/


template <class T>
h_funs<T>::h_funs(int k, int dimens, double w_to_be) {
    M = (long int) floor(pow(2,(double)(32/k)));
    //std::cout << "eftiaksa" << M <<"\n";
    dimensions = dimens;
    w = w_to_be;
    //now si's
    //http://www.cplusplus.com/reference/random/uniform_real_distribution/
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_real_distribution<double>  distr(0, w);
    for (int i=0; i<dimensions; ++i) {
        //sis[i].push_back(distr(generator);
        sis.push_back(distr(generator));
        //fprintf(stderr, "%f\n", sis[i]);
    }
}

template <class T>
h_funs<T>::~h_funs() {}


template <class T>
long int h_funs<T>::individual_comp(long int ai, int expon){//expon o ek8eths tou m
  //std::cout << "eftiaksa" << M <<"\n";
  int component_1 = mod_pow(m, expon, M);
  int component_2 = our_mod(ai, M);
  long int result = our_mod(component_1*component_2, M); //de xreiazetai exponentiation, kai ta 2 components mikrotera tou M
  return result;
}



template <class T>
long int h_funs<T>::actual_h_function(my_vector<T> x) {
    std::vector<T> the_v = x.get_v();

    long int result_part = 0; //h ontothta m^d * ai mod M, sto telos 8a a8roistoun auta gia to teliko apotelesma ths h
    long int ai = 0;
    long int result = 0;
    for(int i=dimensions-1; i>=0; i--){
        ai = (long int)floor((the_v[i] - sis[i])/w);
        result_part = individual_comp(ai, dimensions-1-i);
        result += result_part;
    }
    result = our_mod(result, M);

   return result;
}


long int our_mod (long int a, long long int b){ //returns remainder as it should
	return (a % b + b) % b;
}


long int mod_pow(long int b, int e, long long int m) { //kaluptei thn akraia periptwsh k=1 => M = 2^32
	long int c = 1;
	if (m==1) {
		return 0;
	}
	else {
		for (int i = 0; i < e; i++) {
			c = our_mod(c*b, m);
		}

		return c;
	}
}


template class h_funs<float>;
template class h_funs<int>;
template class h_funs<double>;
