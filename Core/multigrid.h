#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

#ifndef GAUSS_SEIDEL_H
#define GAUSS_SEIDEL_H

inline int ipow(int input, int n){
	return std::round(std::pow(input, n));
}

inline int iroot(int input, int n){
	return std::round(std::pow(input, 1./n));
}

class Multigrid{
	private:
		int ndim;//Number of dimensions
		int nn;//Number of points along the dimension - 1
		int ng;//nn = 2 ** ng

		int usize;//Total size of the array u
		std::vector<double> u;//Array containing the solution
		std::vector<std::vector<double> > rhs;//Array containing right hand side with the restrictions

		int npre; //Number of relaxation preiterations
		int npost; //Number of relaxation postiterations
	public:
		Multigrid(int ndimi, int nni) : nn(nni), ndim(ndimi), usize(1), ng(0), npre(1), npost(1){
			while(nni >>= 1) ng++;
			if(nn != (1 << ng)) throw("nn must be a power of 2 in multigrid\n");
			for(int i = 0; i < ndim; i++) usize *= nn+1;
		}

		void initialize_lhs(void){u.assign(usize, 0.0);}

		//Defined in file_operations.cpp
		void save_to_file(const char *folder_name);
		void save_array(std::vector<double> &arr, const char *folder_name);

		//Defined in function_fill.cpp
		void fill(std::vector<double> &array, double (*func)(std::vector<double> &));
		void fill_lhs(double (*lhs_func)(std::vector<double> &)){fill(u, lhs_func);}
		void fill_rhs(double (*rhs_func)(std::vector<double> &));

		double compare_lhs(double (*lhs_func)(std::vector<double> &));

		void print_index(int n, int i); //Technical function that decomposes index i to dimensions
		//Gauss-Seidel step, works only if number of points is equal in all dimensions
		double gs_step(std::vector<double> &arr, std::vector<double> &rhs);

		int gauss_seidel(void);
		
		//Multigrid routines
		void slvsml(std::vector<double> &arr, std::vector<double> &rhs);
		void addint(std::vector<double> &uf, std::vector<double> &uc);
		void rstrct_half(std::vector<double> &uc, std::vector<double> &uf);
		void rstrct_full(std::vector<double> &uc, std::vector<double> &uf);
		double res_rstrct_full(std::vector<double> &uc, std::vector<double> &uf, std::vector<double> &rhs);

		//V-cycle routines
		void set_relax(int n){npre = npost = n;}
		double multigrid_iteration(int j, std::vector<double> &uj, std::vector<double> &rhsj);
		double v_cycle(void){return multigrid_iteration(ng-1, u, rhs[ng-1]);}
		double full_multigrid(int ncycle);
};

#endif