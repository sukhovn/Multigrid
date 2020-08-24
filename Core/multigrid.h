#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

#ifndef GAUSS_SEIDEL_H
#define GAUSS_SEIDEL_H

class Multigrid{
	private:
		int ndim;
		std::vector<int> nn;
		double dx2;

		int usize;
		std::vector<double> u;
		std::vector<double> f;
	public:
		Multigrid(std::vector<int> &nni) : nn(nni), ndim(nni.size()), usize(1){
			dx2 = std::pow(1.0/(1.0*nn[0]), 2.0);
			for(int i = 0; i < ndim; i++) usize *= nn[i]+1;

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
		double gs_step(double omega);

		int gauss_seidel(void);
};

#endif
