#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <string.h>
#include <multigrid.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

void Multigrid::save_to_file(const char *folder_name){
	//Creating output folder
	struct stat sttus = {0};
	if (stat(folder_name, &sttus) == -1) {
    	mkdir(folder_name, 0700);
	}

	char file_name[80];
	std::ofstream out;
	
	sprintf(file_name, "%s/x.dat", folder_name);
	out.open(file_name);

	for(int dm = 0; dm < ndim; dm++){
		for(int i = 0; i <= nn; i++)
			out << 1.0*i/(1.0*nn) << " ";
		out << "\n";
	}
	
	out.close();

	sprintf(file_name, "%s/u.dat", folder_name);
	out.open(file_name);

	for(int i = 0; i < usize; i++) out << u[i] << " ";
	
	out.close();
	
	return;
}

void Multigrid::save_array(std::vector<double> &arr, const char *folder_name){
	//Creating output folder
	struct stat sttus = {0};
	if (stat(folder_name, &sttus) == -1) {
    	mkdir(folder_name, 0700);
	}

	char file_name[80];
	std::ofstream out;
	
	sprintf(file_name, "%s/x.dat", folder_name);
	out.open(file_name);

	int n = iroot(arr.size(), ndim)-1;
	for(int dm = 0; dm < ndim; dm++){
		for(int i = 0; i <= n; i++)
			out << 1.0*i/(1.0*n) << " ";
		out << "\n";
	}
	
	out.close();

	std::cout << arr.size() << " " << n << "\n";
	sprintf(file_name, "%s/u.dat", folder_name);
	out.open(file_name);

	for(int i = 0; i < arr.size(); i++) out << arr[i] << " ";
	
	out.close();
	
	return;
}