#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <gauss_seidel.h>
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
		// if(strcmp(argv[i], "-steps") == 0){
		// 	steps = atoi(argv[i+1]);
		// 	i++;
		// }
		if(strcmp(argv[i], "-out") == 0){
			strcpy(folder_name, argv[i+1]);
			ifsave = true;
			i++;
		}	
	}
	std::vector<int> nn;
	nn.assign(ndim, nx);

	Gauss_Seidel gs(nn);
	
	gs.fill_lhs(test_function);
	gs.fill_rhs(test_function_laplacian);
	
	clock_t t; 
    t = clock(); 
	steps = gs.sor();
	t = clock() - t;
	double duration = ((double) t)/CLOCKS_PER_SEC;
	
	std::cout << "The number of steps is " << steps << "\n";
	std::cout << "The calculation took " << duration << " seconds\n";
	std::cout << "The difference between the calculated and the exact result is: ";
	std::cout << std::scientific << gs.compare_lhs(test_function) << "\n";

	if(ifsave) gs.save_to_file(folder_name);
	
	return;
}

int main(int argc, char const *argv[]){
	test(argc, argv);
	return 0;
}