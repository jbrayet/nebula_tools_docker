/*************************************************************************
Copyright (c) 2013, Haitham ASHOOR

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses

>>> END OF LICENSE >>>
*************************************************************************/
#include <iostream>
#include <vector>
#include<iterator>
#include<algorithm>
#include <string>
#include <iostream>
#include <iomanip>
using namespace std;


#include "Profiler.h"
#include"Reader.h"
#include "HMCan.h"
#include "Parser.h"
const double HMCAN_VERSION = 1.20;

static void print_version()
{
    std::ostringstream ostr;
     ostr << "HMCan v" << std::fixed << std::setprecision(2) << HMCAN_VERSION << " : detection of histone modifications in normal and cancer ChIP-seq datasets\n";
     std::cout << ostr.str();
}


static void check_files(Parser p){
	ifstream test;

	//check GC index file
	test.open(p.GC_index.c_str());
	if (!test){
		cerr<<"Error: Can not find GC index file at "<<p.GC_index<<endl;
		exit(1);
	}
}


int main(int argc, char * argv[]){
	
	print_version();
	if (argc<5){
		cout<<"Not enough arguments"<<endl<<"Usage: ./HMCan <TargetFile> <ControlFile> <configuration file> <Name>"<<endl;
		cout<<"For more details please see README.md"<<endl;
		return 1;
	}
	Parser p = Parser();
	p.parse(argv[3]);

	check_files(p);

	Reader data(p.t,10,true), control(p.t,10,true);
	vector<vector<float> > emission,transition;
	vector<observation_seq> obs;
	cout<<"Reading chip...."<<endl;
	map<string,vector<DNA_fragment> > target_data = data.Read(argv[1]);
	cout<<"Reading control...."<<endl;
	map<string,vector<DNA_fragment> > control_data = control.Read(argv[2]);
	string name = argv[4];


	vector<BedEntry> black_list = data.Read_blacklist(p.blacklist.c_str());

	int step = (int) ceil(float(p.max_length)/p.bin_size);

	Profiler test_gc = 	Profiler(p.GC_index.c_str(),p.min_length,p.median_length,p.max_length ,p.bin_size,p.large_bin_size,p.merge_dist,p.pvalue_threshold);
	test_gc.build_profile(target_data,control_data,p.path,black_list);
	

	target_data.clear();
	control_data.clear();
	emission = test_gc.get_emission();
	transition = test_gc.get_transition();
	obs = test_gc.get_obs();
	int min = test_gc.get_minObs();
	HMCan a = HMCan(obs,transition,emission,p.iter_threshold,p.final_threshold,step,min,p.bin_size,p.posterior_threshold);
	a.run(p.merge_dist,p.max_iter);

	if (p.wig)
		test_gc.print_wig(name);
	a.print(name,p.median_length);
	
	if (p.posterior){
		a.print_posterior(name);
	}

	return 0;

}
