#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <multigrid.h>
#include <string.h>
#include <time.h>

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

void test_v_cycle(int argc, char const *argv[]){
	int ndim = 2, np = 5;
	int nrelax = 1, ncycles = 1;
	double err;

	bool ifsave = false;
	char folder_name[80];
	for(int i = 1; i < argc; i++){
		if(strcmp(argv[i], "-ndim") == 0){
			ndim = atoi(argv[i+1]);
			i++;
		}
		if(strcmp(argv[i], "-np") == 0){
			np = atoi(argv[i+1]);
			i++;
		}
		if(strcmp(argv[i], "-nrelax") == 0){
			nrelax = atoi(argv[i+1]);
			i++;
		}
		if(strcmp(argv[i], "-ncycles") == 0){
			ncycles = atoi(argv[i+1]);
			i++;
		}
		if(strcmp(argv[i], "-out") == 0){
			strcpy(folder_name, argv[i+1]);
			ifsave = true;
			i++;
		}	
	}

	Multigrid gs(ndim, ipow(2, np));
	
	// gs.fill_lhs(test_function);
	gs.fill_rhs(test_function_laplacian);
 	gs.set_relax(nrelax);
 	gs.initialize_lhs();

	clock_t t; 
    t = clock(); 
    for(int i = 0; i < ncycles; i++)
   		err = gs.v_cycle();
	t = clock() - t;
	double duration = ((double) t)/CLOCKS_PER_SEC;
	
	std::cout << "The calculation took " << duration << " seconds\n";
	std::cout << "The error estimate is " << err << "\n";
	std::cout << "The difference between the calculated and the exact result is: ";
	std::cout << std::scientific << gs.compare_lhs(test_function) << "\n";

	if(ifsave) gs.save_to_file(folder_name);
	
	return;
}

void test_full_multigrid(int argc, char const *argv[]){
	int ndim = 2, np = 5;
	int nrelax = 1, ncycles = 1;
	double err;

	bool ifsave = false;
	char folder_name[80];
	for(int i = 1; i < argc; i++){
		if(strcmp(argv[i], "-ndim") == 0){
			ndim = atoi(argv[i+1]);
			i++;
		}
		if(strcmp(argv[i], "-np") == 0){
			np = atoi(argv[i+1]);
			i++;
		}
		if(strcmp(argv[i], "-nrelax") == 0){
			nrelax = atoi(argv[i+1]);
			i++;
		}
		if(strcmp(argv[i], "-ncycles") == 0){
			ncycles = atoi(argv[i+1]);
			i++;
		}
		if(strcmp(argv[i], "-out") == 0){
			strcpy(folder_name, argv[i+1]);
			ifsave = true;
			i++;
		}	
	}

	Multigrid gs(ndim, ipow(2, np));
	
	// gs.fill_lhs(test_function);
	gs.fill_rhs(test_function_laplacian);
 	gs.set_relax(nrelax);

	clock_t t; 
    t = clock(); 
    err = gs.full_multigrid(ncycles);
	t = clock() - t;
	double duration = ((double) t)/CLOCKS_PER_SEC;
	
	std::cout << "The calculation took " << duration << " seconds\n";
	std::cout << "The error estimate is " << err << "\n";
	std::cout << "The difference between the calculated and the exact result is: ";
	std::cout << std::scientific << gs.compare_lhs(test_function) << "\n";

	if(ifsave) gs.save_to_file(folder_name);
	
	return;
}

int main(int argc, char const *argv[]){
	for(int i = 1; i < argc; i++){
		if(strcmp(argv[i], "-v_cycle") == 0){
			test_v_cycle(argc, argv);			
			return 0;
		}
		if(strcmp(argv[i], "-full") == 0){	
			test_full_multigrid(argc, argv);
			return 0;
		}
	}
	return 0;
}