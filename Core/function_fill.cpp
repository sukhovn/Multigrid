#include <iostream>
#include <vector>
#include <cmath>
#include <gauss_seidel.h>

void Gauss_Seidel::fill_lhs(double (*lhs_func)(std::vector<double> &)){
	int i, ii;
	std::vector<double> x;
	x.resize(ndim);

	for(i = ii = 0; i < usize; ii = ++i){
		for(int dm = ndim-1; dm >= 0; dm--){
			if(ii%(nn[dm]+1)){
				x[dm] += 1.0/(1.0*nn[dm]);
				break;
			}
			x[dm] = 0.0;
			ii /= (nn[dm]+1);
		}
		u[i] = (*lhs_func)(x);
	}
	return;
}

void Gauss_Seidel::fill_rhs(double (*rhs_func)(std::vector<double> &)){
	int i, ii;
	std::vector<double> x;
	x.resize(ndim);

	for(i = ii = 0; i < usize; ii = ++i){
		for(int dm = ndim-1; dm >= 0; dm--){
			if(ii%(nn[dm]+1)){
				x[dm] += 1.0/(1.0*nn[dm]);
				break;
			}
			x[dm] = 0.0;
			ii /= (nn[dm]+1);
		}
		f[i] = (*rhs_func)(x);
	}
	return;
}

double Gauss_Seidel::compare_lhs(double (*lhs_func)(std::vector<double> &)){
	int i, ii;
	double diff, max = 0.0;
	std::vector<double> x;
	x.resize(ndim);

	for(i = ii = 0; i < usize; ii = ++i){
		for(int dm = ndim-1; dm >= 0; dm--){
			if(ii%(nn[dm]+1)){
				x[dm] += 1.0/(1.0*nn[dm]);
				break;
			}
			x[dm] = 0.0;
			ii /= (nn[dm]+1);
		}
		diff = fabs(u[i] - (*lhs_func)(x));
		if(max < diff) max = diff;
	}
	return max;
}