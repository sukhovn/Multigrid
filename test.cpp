#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <multigrid.h>
#include <string.h>
#include <time.h>

//Simple routine to check if n is a power of 2
int check_pow_2(int n){
	if(n < 2) return 1;
	while(n > 1){
		if(n % 2) return 1;
		n /= 2;
	}
	return 0;
}

double test_function(std::vector<double> &x){
	double val = 1.0;
	for(int dm = 0; dm < x.size(); dm++)
		val *= std::sin(x[dm]*M_PI);
	return val;
}

double test_function_laplacian(std::vector<double> &x){
	double val = (-1.0) * x.size() * std::pow(M_PI, 2.0);
	for(int dm = 0; dm < x.size(); dm++)
		val *= std::sin(x[dm]*M_PI);
	return val;
}

void test(int argc, char const *argv[]){
	int ndim = 2, nx = 32, steps = 100;
	
	bool ifsave = false;
	char folder_name[80];
	for(int i = 1; i < argc; i++){
		if(strcmp(argv[i], "-ndim") == 0){
			ndim = atoi(argv[i+1]);
			i++;
		}
		if(strcmp(argv[i], "-nx") == 0){
			nx = atoi(argv[i+1]);
			i++;
		}
		if(strcmp(argv[i], "-out") == 0){
			strcpy(folder_name, argv[i+1]);
			ifsave = true;
			i++;
		}	
	}

	Multigrid gs(ndim, nx);
	
	// gs.fill_lhs(test_function);
	gs.fill_rhs(test_function_laplacian);
	
	clock_t t; 
    t = clock(); 
	steps = gs.gauss_seidel();
	t = clock() - t;
	double duration = ((double) t)/CLOCKS_PER_SEC;
	
	std::cout << "The number of steps is " << steps << "\n";
	std::cout << "The calculation took " << duration << " seconds\n";
	std::cout << "The difference between the calculated and the exact result is: ";
	std::cout << std::scientific << gs.compare_lhs(test_function) << "\n";

	if(ifsave) gs.save_to_file(folder_name);
	
	return;
}

inline int ipow(int input, int n){
	return std::round(std::pow(input, n));
}

void test_rstrct(int argc, char const *argv[]){
	int ndim = 2, nx = 32, steps = 100;
	int nxc = 4, jf = 12;
	
	bool ifsave = false;
	char folder_name[80];
	for(int i = 1; i < argc; i++){
		if(strcmp(argv[i], "-ndim") == 0){
			ndim = atoi(argv[i+1]);
			i++;
		}
		if(strcmp(argv[i], "-nx") == 0){
			nx = atoi(argv[i+1]);
			i++;
		}
		if(strcmp(argv[i], "-nxc") == 0){
			nxc = atoi(argv[i+1]);
			if(check_pow_2(nxc)) throw("nxc is not a power of 2");
			i++;
		}
		if(strcmp(argv[i], "-j") == 0){
			jf = atoi(argv[i+1]);
			i++;
		}
	}

	Multigrid gs(ndim, nx);
	std::vector<double> uc, uf;
	uc.assign(ipow(nxc+1, ndim), 0.0);
	uf.assign(ipow(2*nxc+1, ndim), 0.0);
	uf[jf] = 1.0;

	// gs.rstrct(uc, uf);
	// gs.print_index(2*nxc, jf);
	// for(int i = 0; i < uc.size(); i++){
	// 	if(uc[i] != 0.0){
	// 		gs.print_index(nxc, i);
	// 		std::cout << uc[i] << "\n";
	// 	}
	// }
	gs.addint(uf, uc);

	return;
}

int main(int argc, char const *argv[]){
	// test(argc, argv);
	test_rstrct(argc, argv);
	return 0;
}