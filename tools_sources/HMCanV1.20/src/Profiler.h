#ifndef PROFILER_H_
#define PROFILER_H_

#define LARGE_BIN_SIZE 100000

#include <map>
#include <string>
#include<vector>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <numeric>
#include <iostream>
#include <cstdlib>
#include "utils.h"

class Profiler {
public:
	Profiler(const char* ,int , int, int ,int ,int ,int,float);
	virtual ~Profiler();
	void build_profile(std::map<std::string,std::vector<DNA_fragment > >& , std::map<std::string,std::vector<DNA_fragment > >&,
			std::string&, std::vector<BedEntry>&);
	std::vector<std::vector<float> > get_transition();
	std::vector<std::vector<float> > get_emission();
	std::vector<observation_seq>  get_obs();
	int get_minObs();
	void print_wig(std::string);

private:

	int med,bin_length,left,right,large_bin_size,merge_dist,min_obs,max_obs;
	float h,slope1,slope2,b1,b2,left_area,right_area,reads_ratio,pvalue_threshold;
	std::vector<std::vector<float> > transition,emission;
	std::map<std::string,std::vector<float> > GC_profile,medians,notNprofile;
	std::map<std::string,int > sizes,sampled_sizes;
	std::map<std::string,float*> sampled_target_density,sampled_control_density,sampled_GC,peaks;
	std::vector<observation_seq> observations;



	void build_single_profile(std::vector<DNA_fragment >& , std::string, bool);
	void read_gc_profile(std::ifstream& );
	void calculate_nonNs_ratio(std::string chr, std::string& chr_seq);
	void call_freec(std::map<std::string,std::vector<DNA_fragment > >&);
	void extend_reads(std::vector<DNA_fragment >&,int);
	void calculate_density_coefs();
	void Normalize_GC();
	void calculate_GC_sampled_bins(std::string,std::string&);
	inline float calculate_GC(std::string&,int,int);
	void print_values();
	void print_vector(float *, int, const char *);
	void derive_transition_probabilities();
	//void derive_emision_probabilities();
	inline void count_emissions();

	void generate_observations();
	void call_peaks();
	void peaks_ideal();
	void clean_reads(std::map<std::string,std::vector<DNA_fragment > >&,
			std::map<std::string,std::vector<DNA_fragment > >&);
	void remove_blacklist(std::vector<BedEntry>&,std::map<std::string,float*>& );



};

#endif /* PROFILER_H_ */
