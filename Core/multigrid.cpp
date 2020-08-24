#include <iostream>
#include <vector>
#include <cmath>
#include <multigrid.h>

#define MAX_STEP 10000
#define ERROR_THRS 1.0e-12

inline int root(int input, int n){
	return std::round(std::pow(input, 1./n));
}

void Multigrid::print_index(int i){
	int urem = usize, ii = i;
	for(int dm = 0; dm < ndim; dm++){
		urem /= nn+1;
		std::cout << ii/urem << " ";
		ii = ii%urem;
	}
	std::cout << "   " << i << "\n";
	return;
}

double Multigrid::gs_step(std::vector<double> &arr, std::vector<double> &rhs){
	int n = root(arr.size(), ndim)-1;
	double dx2 = std::pow(1.0/(1.0*n), 2.0);
	double diff, max = 0.0;
	int ii, urem, ishift = 1;
	int ist = 0, iend;
	for(int dm = ndim - 1; dm >= 0; dm--){
		ist += ishift; 
		ishift *= n+1;
	}
	int ifinal = arr.size() - ist - 1;

	while(ist < ifinal){
		iend = ist + n - 2;
		for(int i = ist; i <= iend; i++){//Calculation main body
			diff = - 2.0 * ndim * arr[i];
			ishift = 1;
			for(int dm = ndim-1; dm >= 0; dm--){
				diff += arr[i+ishift] + arr[i-ishift];
				ishift *= n + 1;
			}
			diff -= dx2 * rhs[i];
			if(max < fabs(diff)) max = fabs(diff);
			arr[i] += diff/(2.0*ndim);
		}

		ist = (ii = iend) + 1;
		ishift = 2;
		for(int dm = ndim-1; dm > 0; dm--){
			if(ii%(n+1) != n-1) break;
			ii /= n+1;
			ist += ishift;
			ishift *= n+1; 
		}
	}

	return max;
}

int Multigrid::gauss_seidel(void){
	int i = 0;
	double error;
	while(i < MAX_STEP){
		error = gs_step(u, f);
		i++;
		if(error < ERROR_THRS) break;
	}

	return i;
}