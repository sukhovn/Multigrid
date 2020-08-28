#include <iostream>
#include <vector>
#include <cmath>
#include <multigrid.h>

#define MAX_STEP 10000
#define ERROR_THRS 1.0e-12

void Multigrid::print_index(int n, int i){
	int urem = ipow(n+1, ndim);
	int ii = i;
	for(int dm = 0; dm < ndim; dm++){
		urem /= n+1;
		std::cout << ii/urem << " ";
		ii = ii%urem;
	}
	std::cout << "   " << i << "\n";
	return;
}

double Multigrid::gs_step(std::vector<double> &arr, std::vector<double> &rhs){
	int n = iroot(arr.size(), ndim);
	double dx2 = std::pow(1.0/(1.0*(n-1)), 2.0);
	double diff, max = 0.0;
	int ii, ishift = 1;
	int ist = 0, iend;
	for(int dm = ndim - 1; dm >= 0; dm--){
		ist += ishift; 
		ishift *= n;
	}
	int ifinal = arr.size() - ist - 1;

	while(ist <= ifinal){
		iend = ist + n - 3;
		for(int i = ist; i <= iend; i++){//Calculation main body
			//i runs through all the points except boundary points
			diff = - 2.0 * ndim * arr[i];
			ishift = 1;
			for(int dm = ndim-1; dm >= 0; dm--){
				diff += arr[i+ishift] + arr[i-ishift];
				ishift *= n;
			}
			diff -= dx2 * rhs[i];
			if(max < fabs(diff)) max = fabs(diff);
			arr[i] += diff/(2.0*ndim);
		}

		ist = (ii = iend) + 1;
		ishift = 2;
		for(int dm = ndim-1; dm > 0; dm--){
			if(ii%n != n-2) break;
			ii /= n;
			ist += ishift;
			ishift *= n; 
		}
	}

	return max;
}

//Performes prolongation (linear interpolation) of uc and adds it to uf
void Multigrid::addint(std::vector<double> &uf, std::vector<double> &uc){
	int fsize = uf.size();
	int fn = iroot(fsize, ndim);
	int cn = fn/2+1;
	bool *parity = new bool[ndim];
	int coef;

	int ii, itmp, ishift = 1, j;
	int ist = 0, iend;
	for(int dm = ndim - 1; dm >= 0; dm--){
		ist += ishift; 
		ishift *= fn;
	}
	int ifinal = fsize - ist - 1;

	while(ist <= ifinal){
		iend = ist + fn - 3;
		for(int i = ist; i <= iend; i++){//Calculation main body
			//i runs through all the points except boundary points

			ii = i;	j = 0; ishift = 1; coef = 1;
			for(int dm = 0; dm < ndim; dm++){
				itmp = ii%fn;
				if((parity[dm] = itmp%2)) coef <<= 1;
				itmp >>= 1;
				j += ishift*itmp;
				ishift *= cn;
				ii /= fn;
			}
			
			for(int k = 0; k < coef; k++){
				//j runs through even neighboring points
				//main computational body
				uf[i] += uc[j]/(1.0*coef);

				ii = k; ishift = 1;
				for(int dm = 0; dm < ndim; dm++){
					if(parity[dm]){
						if(ii%2) j-= ishift;
						else{j += ishift; break;}
						ii >>= 1;
					}
					ishift *= cn;
				}
			}
		}

		ist = (ii = iend) + 1;
		ishift = 2;
		for(int dm = ndim-1; dm > 0; dm--){
			if(ii%(fn) != fn-2) break;
			ii /= fn;
			ist += ishift;
			ishift *= fn; 
		}
	}

	return;
}

//Half-weighting restriction, uc boundary not affected
void Multigrid::rstrct_half(std::vector<double> &uc, std::vector<double> &uf){
	int cn = iroot(uf.size(), ndim)/2;
	int csize = ipow(cn+1, ndim);
	double coef = 1.0/(4.0*ndim);

	int ii, ishift = 1, jshift = 1;
	int ist = 0, iend, i;
	int jst = 0, j;
	for(int dm = ndim - 1; dm >= 0; dm--){
		ist += ishift;
		jst += 2*jshift;
		ishift *= cn+1;
		jshift *= 2*cn+1;
	}
	int ifinal = csize - ist - 1;

	while(ist <= ifinal){
		iend = ist + cn - 2;
		for(i = ist, j = jst; i <= iend; i++, j+=2){//Calculation main body
			//i runs through all the points except boundary points
			//if i is (a1, a2 ...), j is (2*a1, 2*a2 ...)
			uc[i] = 0.5*uf[j];
			jshift = 1;
			for(int dm = ndim - 1; dm >= 0; dm--){
				uc[i] += coef * (uf[j+jshift] + uf[j-jshift]);
				jshift *= 2*cn+1;
			}
		}

		ist = (ii = iend) + 1;
		jst = j - 1;
		ishift = 2;
		jshift = 1;
		for(int dm = ndim-1; dm > 0; dm--){
			if(ii%(cn+1) != cn-1) break;
			ii /= cn+1;
			ist += ishift;
			jst += 4*jshift;
			ishift *= cn+1;
			jshift *= 2*cn+1; 
		}
		jst += jshift;
	}
	return;
}

//Full-weight restriction, uc boundary not affected
void Multigrid::rstrct_full(std::vector<double> &uc, std::vector<double> &uf){
	int fsize = uf.size();
	int fn = iroot(fsize, ndim);
	int cn = fn/2+1;
	bool *parity = new bool[ndim];
	int jiter;
	double coef, factor = ipow(2, ndim);

	int ii, itmp, ishift = 1, j;
	int ist = 0, iend;
	for(int dm = ndim - 1; dm >= 0; dm--){
		ist += ishift; 
		ishift *= fn;
	}
	int ifinal = fsize - ist - 1;

	while(ist <= ifinal){
		iend = ist + fn - 3;
		for(int i = ist; i <= iend; i++){//Calculation main body
			//i runs through all the points except boundary points on a fine grid
			ii = i;	j = 0; ishift = 1; coef = 1.0; jiter = 1;
			for(int dm = 0; dm < ndim; dm++){
				itmp = ii%fn;
				parity[dm] = itmp%2;
				itmp >>= 1;
				if(parity[dm]){
					coef *= 2.0;
					if(itmp == 0){parity[dm] = 0; j += ishift;}
					else if(itmp == cn-2) parity[dm] = 0;
					else jiter <<= 1;
				}
				j += ishift*itmp;
				ishift *= cn;
				ii /= fn;
			}
			
			for(int k = 0; k < jiter; k++){
				//j runs through even neighboring points on a coarse grid
				//main computational body
				uc[j] += uf[i]/(factor*coef);

				ii = k; ishift = 1;
				for(int dm = 0; dm < ndim; dm++){
					if(parity[dm]){
						if(ii%2) j-= ishift;
						else{j += ishift; break;}
						ii >>= 1;
					}
					ishift *= cn;
				}
			}
		}

		ist = (ii = iend) + 1;
		ishift = 2;
		for(int dm = ndim-1; dm > 0; dm--){
			if(ii%(fn) != fn-2) break;
			ii /= fn;
			ist += ishift;
			ishift *= fn; 
		}
	}

	return;
}

//The routine finds residual and performs a restriction of it
//Full-weight restriction, uc boundary not affected
double Multigrid::res_rstrct_full(std::vector<double> &uc, std::vector<double> &uf, std::vector<double> &rhs){
	int fsize = uf.size();
	int fn = iroot(fsize, ndim);
	int cn = fn/2+1;
	bool *parity = new bool[ndim];
	int jiter;
	double coef, factor = ipow(2, ndim); //Factor is a common dimension dependent factor, coef depend on the number of neighboring points to j

	double dx2 = std::pow(1.0/(1.0*(fn-1)), 2.0);
	double diff, max = 0.0;

	int ii, itmp, ishift = 1, j;
	int ist = 0, iend;
	for(int dm = ndim - 1; dm >= 0; dm--){
		ist += ishift; 
		ishift *= fn;
	}
	int ifinal = fsize - ist - 1;

	while(ist <= ifinal){
		iend = ist + fn - 3;
		for(int i = ist; i <= iend; i++){//Calculation main body
			//i runs through all the points except boundary points on a fine grid
			//Gauss-Seidel part of the routine
			diff = - 2.0 * ndim * uf[i];
			ishift = 1;
			for(int dm = ndim-1; dm >= 0; dm--){
				diff += uf[i+ishift] + uf[i-ishift];
				ishift *= fn;
			}
			diff -= dx2 * rhs[i];
			if(max < fabs(diff)) max = fabs(diff);

			//Restriction part of the routine
			ii = i;	j = 0; ishift = 1; coef = 1.0; jiter = 1;
			for(int dm = 0; dm < ndim; dm++){
				itmp = ii%fn;
				parity[dm] = itmp%2;
				itmp >>= 1;
				if(parity[dm]){
					coef *= 2.0;
					if(itmp == 0){parity[dm] = 0; j += ishift;}
					else if(itmp == cn-2) parity[dm] = 0;
					else jiter <<= 1;
				}
				j += ishift*itmp;
				ishift *= cn;
				ii /= fn;
			}
			
			for(int k = 0; k < jiter; k++){
				//j runs through even neighboring points on a coarse grid
				//main computational body
				uc[j] += (-1.0)*diff/(factor*coef*dx2);
				ii = k; ishift = 1;
				for(int dm = 0; dm < ndim; dm++){
					if(parity[dm]){
						if(ii%2) j-= ishift;
						else{j += ishift; break;}
						ii >>= 1;
					}
					ishift *= cn;
				}
			}
		}

		ist = (ii = iend) + 1;
		ishift = 2;
		for(int dm = ndim-1; dm > 0; dm--){
			if(ii%(fn) != fn-2) break;
			ii /= fn;
			ist += ishift;
			ishift *= fn; 
		}
	}

	return max;
}

void Multigrid::slvsml(std::vector<double> &arr, std::vector<double> &rhs){
	int i = 0, ishift = 1;
	for(int dm = ndim-1; dm >= 0; dm--){
		i += ishift;
		ishift *= 3;
	}
	double diff = - 0.25 * rhs[i];
	
	ishift = 1;
	for(int dm = ndim-1; dm >= 0; dm--){
		diff += arr[i+ishift] + arr[i-ishift];
		ishift *= 3;
	}
			
	arr[i] = diff/(2.0*ndim);
	return;
}

double Multigrid::multigrid_iteration(int j, std::vector<double> &uj, std::vector<double> &rhsj){
	double err;
	if(j == 0){
		slvsml(uj, rhsj);
	}
	else{
		int len = (1 << j) + 1;
		int size = ipow(len, ndim);
		std::vector<double> v, res;
		v.assign(size, 0.0);
		res.assign(size, 0.0);

		for(int i = 0; i < npre; i++)
			gs_step(uj, rhsj);
		res_rstrct_full(res, uj, rhsj);
		multigrid_iteration(j-1, v, res);
		addint(uj, v);
		for(int i = 0; i < npost; i++)
			err = gs_step(uj, rhsj);	
	}
	return err;
}

double Multigrid::full_multigrid(int ncycle){
	double err;
	int n = 3;
	int size = ipow(n, ndim);
	std::vector<double> ui, ui1;
	ui.assign(size, 0.0);
	slvsml(ui, rhs[0]);

	for(int i = 1; i < ng; i++){
		n = 2*n-1;
		size = ipow(n, ndim);
		ui1 = ui;
		ui.assign(size, 0.0);
		addint(ui, ui1);
		for(int j = 0; j < ncycle; j++)
			err = multigrid_iteration(i, ui, rhs[i]);
	}
	u = ui;
	return err;
}