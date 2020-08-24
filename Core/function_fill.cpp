#include <iostream>
#include <vector>
#include <cmath>
#include <multigrid.h>

void Multigrid::fill(std::vector<double> &array, double (*func)(std::vector<double> &)){
	int i, ii;
	std::vector<double> x;
	x.resize(ndim);

	for(i = ii = 0; i < usize; ii = ++i){
		for(int dm = ndim-1; dm >= 0; dm--){
			if(ii%(nn+1)){
				x[dm] += 1.0/(1.0*nn);
				break;
			}
			x[dm] = 0.0;
			ii /= (nn+1);
		}
		array[i] = (*func)(x);
	}
	return;
}

double Multigrid::compare_lhs(double (*lhs_func)(std::vector<double> &)){
	int i, ii;
	double diff, max = 0.0;
	std::vector<double> x;
	x.resize(ndim);

	for(i = ii = 0; i < usize; ii = ++i){
		for(int dm = ndim-1; dm >= 0; dm--){
			if(ii%(nn+1)){
				x[dm] += 1.0/(1.0*nn);
				break;
			}
			x[dm] = 0.0;
			ii /= (nn+1);
		}
		diff = fabs(u[i] - (*lhs_func)(x));
		if(max < diff) max = diff;
	}
	return max;
}