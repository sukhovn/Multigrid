#include <iostream>
#include <vector>
#include <cmath>
#include <multigrid.h>

#define MAX_STEP 10000
#define ERROR_THRS 1.0e-12

void Multigrid::print_index(int i){
	int urem = usize, ii = i;
	for(int dm = 0; dm < ndim; dm++){
		urem /= nn[dm]+1;
		std::cout << ii/urem << " ";
		ii = ii%urem;
	}
	std::cout << "   " << i << "\n";
	return;
}

double Multigrid::gs_step(double omega){
	double diff, max = 0.0;
	int ii, urem, ishift = 1;
	int ist = 0, iend;
	for(int dm = ndim - 1; dm >= 0; dm--){
		ist += ishift; 
		ishift *= nn[dm]+1;
	}
	int ifinal = usize - ist - 1;

	while(ist < ifinal){
		iend = ist + nn[ndim-1] - 2;
		for(int i = ist; i <= iend; i++){//Calculation main body
			diff = - 2.0 * ndim * u[i];
			ishift = 1;
			for(int dm = ndim-1; dm >= 0; dm--){
				diff += u[i+ishift] + u[i-ishift];
				ishift *= nn[dm] + 1;
			}
			diff -= dx2 * f[i];
			if(max < fabs(diff)) max = fabs(diff);
			u[i] += omega*diff/(2.0*ndim);
		}

		ist = (ii = iend) + 1;
		ishift = 2;
		for(int dm = ndim-1; dm > 0; dm--){
			if(ii%(nn[dm]+1) != nn[dm]-1) break;
			ii /= nn[dm]+1;
			ist += ishift;
			ishift *= nn[dm]+1; 
		}
	}

	return max;
}

int Multigrid::gauss_seidel(void){
	int i = 0;
	double error;
	while(i < MAX_STEP){
		error = gs_step(1.0);
		i++;
		if(error < ERROR_THRS) break;
	}

	return i;
}