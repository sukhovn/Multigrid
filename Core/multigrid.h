#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

#ifndef GAUSS_SEIDEL_H
#define GAUSS_SEIDEL_H

class Multigrid{
	private:
		int ndim;
		int nn;
		int ng;

		int usize;
		std::vector<double> u;
		std::vector<double> f;
	public:
		Multigrid(int ndimi, int nni) : nn(nni), ndim(ndimi), usize(1), ng(0){
			while(nni >>= 1) ng++;
			if(nn != (1 << ng)) throw("nn must be a power of 2 in multigrid\n");
			for(int i = 0; i < ndim; i++) usize *= nn+1;
			u.assign(usize, 0.0);
			f.assign(usize, 0.0);
		}

		//Defined in file_operations.cpp
		void save_to_file(const char *folder_name);

		//Defined in function_fill.cpp
		void fill(std::vector<double> &array, double (*func)(std::vector<double> &));
		void fill_lhs(double (*lhs_func)(std::vector<double> &)){fill(u, lhs_func);}
		void fill_rhs(double (*rhs_func)(std::vector<double> &)){fill(f, rhs_func);}

		double compare_lhs(double (*lhs_func)(std::vector<double> &));

		void print_index(int i); //Technical function that decomposes index i to dimensions
		//Gauss-Seidel step, works only if number of points is equal in all dimensions
		double gs_step(std::vector<double> &arr, std::vector<double> &rhs);

		int gauss_seidel(void);
};

#endif
