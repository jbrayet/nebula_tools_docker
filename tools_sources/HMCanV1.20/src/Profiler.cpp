/*************************************************************************
Copyright (c) 2013, Haitham ASHOOR

>>> SOURCE LICENSE >>>
r
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
#include "Profiler.h"
#include "linreg.h"
#include "segmentation.h"
#include "freec.h"
#include "poissondistr.h"
#include "utils.h"

#include <sstream>
#include <set>
using namespace std;


Profiler::Profiler(const char * GC_index, int left, int med, int right,int bin_length,int large_bin_size,int merge_dist,float pvalue_threshold) {

	ifstream gc;
	string chr;
	gc.open(GC_index);

	if(!gc){
		cerr<<"Error:can not open GC information file"<<endl;
		exit(1);
	}
	read_gc_profile(gc);
	gc.close();
	this->right = right;
	this->left = left;
	this->med = med;
	this->bin_length = bin_length;
	this->large_bin_size = large_bin_size;
	this->merge_dist = merge_dist;
	this->pvalue_threshold = pvalue_threshold;
	calculate_density_coefs();
	min_obs = 0;
	max_obs = 0;
}




Profiler::~Profiler() {

	map<string,float*>::iterator chr_it;
	for(chr_it=sampled_target_density.begin();chr_it!=sampled_target_density.begin();++chr_it)
		delete (*chr_it).second;

	for(chr_it=sampled_control_density.begin();chr_it!=sampled_control_density.begin();++chr_it)
			delete (*chr_it).second;

	for(chr_it=sampled_GC.begin();chr_it!=sampled_GC.begin();++chr_it)
			delete (*chr_it).second;

     for(chr_it=peaks.begin();chr_it!=peaks.end();++chr_it)
    	 delete (*chr_it).second;
}

void Profiler::build_single_profile(vector<DNA_fragment >& tags, string chr, bool targetORcontrol){

	// allocate memory for density
	if (targetORcontrol)
		sampled_target_density[chr] =  new float[sampled_sizes[chr]];
	else
		sampled_control_density[chr] =  new float[sampled_sizes[chr]];

	for(int i=0;i<sampled_sizes[chr];i++)
		if (targetORcontrol)
			sampled_target_density[chr][i] = 0;
		else
			sampled_control_density[chr][i] = 0;


	for(vector<DNA_fragment>::iterator fragment_it=tags.begin(); fragment_it !=tags.end();++fragment_it ){


		int fragment_start = (*fragment_it).start;
		int fragment_left  = (*fragment_it).start+left-1>sizes[chr] ? sizes[chr]-1:(*fragment_it).start+left-1;
		int fragment_med = 	(*fragment_it).start+med-1>sizes[chr]?sizes[chr]-1:(*fragment_it).start+med-1;
		int fragment_right = (*fragment_it).start+right-1>sizes[chr]? sizes[chr]-1:(*fragment_it).start+right-1;


		//start at multiple of bins position
		int start_point = fragment_start;
		while ((start_point+1)%bin_length!=0)
			start_point++;

		//calculate density using triangular distribution
		//check findpeaks program for more info
		for(int i=start_point;i<=fragment_right;i+=bin_length){
			int index = ((i+1)/bin_length)-1;
			if(index>=sampled_sizes[chr]){
				cerr<<"Incorrect reference genome, please check your reference genome"<<endl;
				exit(0);
			}
			int x = i-fragment_start;
				if (i<=fragment_left)
					if(targetORcontrol)
						sampled_target_density[chr][index]+= 1;
					else
						sampled_control_density[chr][index]+=1;

				else if (i<=fragment_med)
				{
					int hx = slope1 * x + b1;
					if (targetORcontrol)
						sampled_target_density[chr][index]+= 1 - ((hx * (x-left)) /2);
					else
						sampled_control_density[chr][index]+=1 - ((hx * (x-left)) /2);
				}
				else{
					int hx = slope2 * x + b2;
					float point_density = (1 - left_area)-(right_area - (hx * (float)(right-x) /2));
					if(targetORcontrol)
						sampled_target_density[chr][index]+=point_density;
					else
						sampled_control_density[chr][index]+=point_density;
				}


		}

	}

		//correct for copy number
		for(int i=0;i<sampled_sizes[chr];i++){
			int index = i*bin_length+bin_length;
			int bin_index = index/large_bin_size;
			float median = medians[chr][bin_index];
			if(median>0.3){ // check for low copy number if found assign unknown label for that
				if (targetORcontrol)
					sampled_target_density[chr][i] =sampled_target_density[chr][i]/median;
				else
					sampled_control_density[chr][i]=sampled_control_density[chr][i]/median;
			}

			/*
			else{
				if (targetORcontrol)
					sampled_target_density[chr][i] =-1;
				else
					sampled_control_density[chr][i]=-1;
			}*/
		}



}

void Profiler::build_profile(map<string,vector<DNA_fragment > >& data, map<string,vector<DNA_fragment > >& control, string& path,
		vector<BedEntry>& blacklist){

	float N=0,M=0;
	map<string,vector<DNA_fragment > >::iterator chr_it;



	clean_reads(data,control);
	for(chr_it = data.begin();chr_it!=data.end();++chr_it){
		N+=(*chr_it).second.size();

	}


	for(chr_it = control.begin();chr_it!=control.end();++chr_it)
		M+=(*chr_it).second.size();

	reads_ratio = M/N;

	cout<<"Reads ratio is: "<<reads_ratio<<endl;
	for(chr_it = data.begin();chr_it!=data.end();++chr_it){
		if ((*chr_it).first.compare("chrM")==0)
			continue;
		string chromosome_path;
		if (path[path.size()-1] == '/')
			chromosome_path = path+(*chr_it).first+".fa";
		else
			chromosome_path = path+"/"+(*chr_it).first+".fa";


		string chr_seq = Read_chr(chromosome_path.c_str());
		sizes[(*chr_it).first] = chr_seq.size();
		sampled_sizes[(*chr_it).first]=sizes[(*chr_it).first]/bin_length+1;
		extend_reads((*chr_it).second,sizes[(*chr_it).first]);
		extend_reads(control[(*chr_it).first],sizes[(*chr_it).first]);

		calculate_GC_sampled_bins((*chr_it).first,chr_seq);
	}

	cout<<"calculating copy number alternations........"<<endl;
	call_freec(control);


	for (chr_it = data.begin();chr_it!=data.end();++chr_it){
		cout<<"Building Density profile for chromosome "<<(*chr_it).first<<"......."<<endl;
		build_single_profile((*chr_it).second,(*chr_it).first,true);
		build_single_profile(control[(*chr_it).first],(*chr_it).first,false);
	}


	remove_blacklist(blacklist,sampled_target_density);
	remove_blacklist(blacklist,sampled_target_density);



	for (chr_it = data.begin();chr_it!=data.end();++chr_it){
		for(int i=0;i<sampled_sizes[(*chr_it).first];i++)
			if (sampled_target_density[(*chr_it).first][i]!=-1)
				sampled_target_density[(*chr_it).first][i] = sampled_target_density[(*chr_it).first][i]*reads_ratio;
		peaks[(*chr_it).first] = new float[sampled_sizes[(*chr_it).first]];

	}





	cout<<"Normalize for GC content and noise ratio........."<<endl;
	call_peaks();
	Normalize_GC();
	call_peaks();
	derive_transition_probabilities();
	//derive_emision_probabilities();
	generate_observations();
	count_emissions();

}


void Profiler::read_gc_profile(ifstream& gc_profile){
	string line;
	while(getline(gc_profile,line)){
		vector<string> tokens;
		istringstream stream(line);
		copy(istream_iterator<string>(stream),
		     istream_iterator<string>(),
			 back_inserter<vector<string> >(tokens));
		if (tokens.size()<5){
			cerr<<"Warning: the line: "<<line<<" ,does not follow format. Line will be ignored"<<endl;
                        continue;
		}
		//string key = "chr"+tokens[0];
		string key = tokens[0];
		float content = atof(tokens[2].c_str());
		float mapability = atof(tokens[4].c_str());
		GC_profile[key].push_back(content);
		notNprofile[key].push_back(mapability);

		}
}


void Profiler::call_freec(map<string,vector<DNA_fragment > >& tags){
	map <string, vector <float> > read_count;
	map <string, vector <float> > ratio_profile;
	map <string, vector <int> > breakpoints;
	int read_counts = 0;
	for (map<string,vector<DNA_fragment > >::iterator chr_it = tags.begin();chr_it!=tags.end();++chr_it){
		read_counts+=(*chr_it).second.size();
			int bins_count =GC_profile[(*chr_it).first].size();
			read_count[(*chr_it).first] = vector<float>(bins_count,0);
			for(vector<DNA_fragment>::iterator v = (*chr_it).second.begin();v!=(*chr_it).second.end();++v){
				int bin_index = (*v).start/large_bin_size;
				if (bin_index >= bins_count){
					cerr<<"Error:" <<"There is a problem with "<<(*chr_it).first<<" size please check your GC Index file!"<<endl;
					exit(1);
				}
				read_count[(*chr_it).first][bin_index]+=1;
			}
		}
	cout<<"Total reads count read by FREEC is: "<<read_counts<<endl;
	int ploidy = 4; // to cover all cases; should not affect negatively cases with ploidy==2 or 3
	recalculateRatioUsingCG ( ploidy, read_count, GC_profile, notNprofile,ratio_profile);
	calculateBreakpoints(ratio_profile, breakpoints);
	calculateCopyNumberMedians(ploidy,breakpoints,medians,ratio_profile);

}



void Profiler:: extend_reads(std::vector<DNA_fragment >& tags, int chr_size){
	for(vector<DNA_fragment>::iterator tag_it = tags.begin(); tag_it!=tags.end();++tag_it)
		if ((*tag_it).strand == 0){
			if ((*tag_it).start>=chr_size){
				cerr<<"Warning: a Read is out side chromosome boundaries; read is ignored"<<endl;
				continue;
			}
		}
		else if ((*tag_it).strand == 1) {
			if ((*tag_it).start-right-1<0)
				(*tag_it).start = 0;
			else
				(*tag_it).start = (*tag_it).start-right+1;
			}
		sort(tags.begin(),tags.end(),compare); // resort to take extension into consideration

}


void Profiler::calculate_density_coefs(){
	h = 2/(float)(right-left);
	slope1 = h/(float)(med-left);
	slope2 = -h/(float)(right-med);
	b1 = -slope1*left;
	b2 = -slope2*right;
	left_area = (h*(float)(med-left))/2;
	right_area = (h*(float)(right-med)/2);

}




void Profiler::Normalize_GC(){

	float chip_accum_density[101],chip_windows_count[101],control_accum_density[101],control_windows_count[101],
			chip_percent[101],control_percent[101],GC_stratum=0,chip_lambda=0,control_lambda=0,max_chip=0,
			max_control=0;
	int startum_windows;
	vector<vector<float> > density;
	map<string,float*>::iterator chr_it;
	vector<float> chip,control;


	//normalize target
	for(int i=0;i<101;i++){
		chip_accum_density[i]=0;
		chip_windows_count[i]=0;
		density.push_back(vector<float>());
	}


	for(chr_it=sampled_target_density.begin();chr_it!=sampled_target_density.end();++chr_it){
		for(int i=0;i<sampled_sizes[(*chr_it).first];i++)
			if (peaks[(*chr_it).first][i]==0&& (*chr_it).second[i]!=-1 && sampled_GC[(*chr_it).first][i]!=-1){
				int index = iround(sampled_GC[(*chr_it).first][i]); //(int)(sampled_GC[(*chr_it).first][i]*100);
				chip_windows_count[index]++;
				density[index].push_back((*chr_it).second[i]);
			}
	}



	for(int i=0;i<101;i++)
			if(chip_windows_count[i]>0){
				chip_windows_count[i]=0;
				int ten_percent = density[i].size()/10;
				sort(density[i].begin(),density[i].end());
				float sum=0;
				for(unsigned int j=ten_percent;j<density[i].size()-ten_percent;j++){
					sum+=density[i][j];
					chip_windows_count[i]++;
				}
				chip_accum_density[i] = sum;
			}

		startum_windows=0;
		for(int i=0;i<27;i++){
			GC_stratum+=chip_accum_density[i];
			startum_windows+=chip_windows_count[i];
		}

		for(int i=0;i<27;i++)
			chip_percent[i] = GC_stratum/startum_windows;


		for(int i=28;i<35;i+=2){
			GC_stratum=(chip_accum_density[i]+chip_accum_density[i-1])/(chip_windows_count[i]+chip_windows_count[i-1]);
			chip_percent[i] = chip_percent[i-1] = GC_stratum;
		}


		for(int i=35;i<66;i++){
			chip_percent[i]=chip_accum_density[i]/chip_windows_count[i];
		}



		for(int i = 66;i<75;i+=2){
			GC_stratum=(chip_accum_density[i]+chip_accum_density[i-1])/(chip_windows_count[i]+chip_windows_count[i-1]);
			chip_percent[i] = chip_percent[i-1] = GC_stratum;
		}

		GC_stratum=0;
		startum_windows=0;
		for(int i=75;i<101;i++){
			GC_stratum+=chip_accum_density[i];
			startum_windows+=chip_windows_count[i];
		}

		for(int i=75;i<101;i++)
			chip_percent[i] = GC_stratum/startum_windows;

	density.clear();



								//normalize for control data
	//*****************************************************************************************************************************

	for(int i=0;i<101;i++){
		control_accum_density[i]=0;
		control_windows_count[i]=0;
		density.push_back(vector<float>());
	}

	for(chr_it=sampled_control_density.begin();chr_it!=sampled_control_density.end();++chr_it){
		for(int i=0;i<sampled_sizes[(*chr_it).first];i++)
			if ((*chr_it).second[i]!=-1 && sampled_GC[(*chr_it).first][i]!=-1){
				int index = iround(sampled_GC[(*chr_it).first][i]);// (int)(sampled_GC[(*chr_it).first][i]*100);
				control_windows_count[index]++;
				density[index].push_back((*chr_it).second[i]);
		}
	}

	for(int i=0;i<101;i++)
		if(control_windows_count[i]>0){
			control_windows_count[i]=0;
			int ten_percent = density[i].size()/10;
			sort(density[i].begin(),density[i].end());
			float sum=0;
			for(unsigned int j=ten_percent;j<density[i].size()-ten_percent;j++){
				sum+=density[i][j];
				control_windows_count[i]++;
			}
			control_accum_density[i] = sum;
		}


	startum_windows=0;
	for(int i=0;i<21;i++){
		GC_stratum+=control_accum_density[i];
		startum_windows+=control_windows_count[i];
	}

	for(int i=0;i<21;i++)
		control_percent[i] = GC_stratum/startum_windows;


	for(int i=22;i<35;i+=2){
		GC_stratum=(control_accum_density[i]+control_accum_density[i-1])/(control_windows_count[i]+control_windows_count[i-1]);
		control_percent[i] = control_percent[i-1] = GC_stratum;
	}

	for(int i=35;i<66;i++){
		control_percent[i]=control_accum_density[i]/control_windows_count[i];
	}


	for(int i = 66;i<75;i+=2){
		GC_stratum=(control_accum_density[i]+control_accum_density[i-1])/(control_windows_count[i]+control_windows_count[i-1]);
		control_percent[i] = control_percent[i-1] = GC_stratum;

	}


	GC_stratum=0;
	startum_windows=0;
	for(int i=75;i<101;i++){
		GC_stratum+=control_accum_density[i];
		startum_windows+=control_windows_count[i];
	}

	for(int i=75;i<101;i++)
		control_percent[i] = GC_stratum/startum_windows;


	for (int i=0;i<101;i++){
		if (chip_percent[i]>max_chip)
			max_chip = chip_percent[i];

		if (control_percent[i]>max_control)
			max_control = control_percent[i];
	}

	int i ,chip_lower,chip_upper,control_lower,control_upper;
	i =0;
	while(chip_percent[i]<0.1*max_chip)
		i++;
	chip_lower =i;

	i=100;
	while(chip_percent[i]<0.1*max_chip)
		i--;
	chip_upper=i;

	i=0;
	while(control_percent[i]<0.1*max_control)
		i++;
	control_lower =i;

	i=100;
	while(control_percent[i]<0.1*max_control)
		i--;
	control_upper=i;

	int final_lower = max(control_lower,chip_lower);
	int final_upper = min(control_upper,chip_upper);

	int chip_total_windows=0,control_total_windows=0;
	for (int i=final_lower;i<=final_upper;i++){
		chip_total_windows+=chip_windows_count[i];
		control_total_windows+=control_windows_count[i];
	}

	// calculate the expectations
	for (int i=final_lower;i<=final_upper;i++){
		chip_lambda+=chip_percent[i]*chip_windows_count[i]/chip_total_windows;
		control_lambda+=control_percent[i]*control_windows_count[i]/control_total_windows;
	}


//correcting here
	for(chr_it=sampled_target_density.begin();chr_it!=sampled_target_density.end();++chr_it){
				for(int i=0;i<sampled_sizes[(*chr_it).first];i++)
					if ((*chr_it).second[i] !=-1 && sampled_GC[(*chr_it).first][i]!=-1){
						int index = (int)(sampled_GC[(*chr_it).first][i]*100);
						if (index>=final_lower && index<=final_upper){
						if (chip_percent[index]!=0)
							(*chr_it).second[i] = (*chr_it).second[i]*chip_lambda/chip_percent[index];
					}
					//else
					//	(*chr_it).second[i]=-1;
			}
	}


	for(chr_it=sampled_control_density.begin();chr_it!=sampled_control_density.end();++chr_it){
			for(int i=0;i< sampled_sizes[(*chr_it).first];i++)
				if ((*chr_it).second[i] !=-1 && sampled_GC[(*chr_it).first][i]!=-1){
					int index = (int)(sampled_GC[(*chr_it).first][i]*100);
					if (index>=final_lower && index<=final_upper){
						if (control_percent[index]!=0)
							(*chr_it).second[i] = (*chr_it).second[i]*control_lambda/control_percent[index];
					}
					//else
					//	(*chr_it).second[i]=-1;
			}
	}



// nosie ratio
	cout<<"chip lambda is: "<<chip_lambda<<endl;
	cout<<"control lambda is "<<control_lambda<<endl;
	float noise_ratio = chip_lambda/control_lambda;
	if (noise_ratio>1)
		noise_ratio = 1;

	cout<<"noise ratio is: "<<noise_ratio<<endl;
	for(chr_it=sampled_control_density.begin();chr_it!=sampled_control_density.end();++chr_it){
		for(int i=0;i< sampled_sizes[(*chr_it).first];i++)
			if ((*chr_it).second[i]!=-1)
				(*chr_it).second[i]=(*chr_it).second[i]*noise_ratio;

			}
}

void Profiler::calculate_GC_sampled_bins(string chr,string& chr_seq){
	sampled_GC[chr] = new float[sampled_sizes[chr]];

	for(int i=0;i<sampled_sizes[chr];i++){
		int position = i*bin_length+bin_length;
		int start = position-med < 0 ? 0:position-med;
		int length = position+med < sizes[chr]? 2*med:sizes[chr]-position+1;
		float GC = calculate_GC(chr_seq,start,length);
		sampled_GC[chr][i] = GC;
	}

}



inline float Profiler::calculate_GC(string& chr_seq,int start,int length){
	//string window = chr_seq.substr(start,length);
	float G=0,C=0;
	for(int j =start;j<start+length;j++)
		if (chr_seq[j]!='N'){
			if (chr_seq[j] =='G' || chr_seq[j] == 'g')
				G++;
			if (chr_seq[j] == 'C' || chr_seq[j] == 'c')
				C++;
		}
		else{
			return -1;
		}

		return (G+C)/length;
}


void Profiler::derive_transition_probabilities(){

	map<string,float *>::iterator chr_it;
	float transition_prob[2][2];
	transition_prob[0][0]=0;
	transition_prob[0][1]=0;
	transition_prob[1][0]=0;
	transition_prob[1][1]=0;
	for(chr_it=sampled_target_density.begin();chr_it!=sampled_target_density.end();++chr_it){
			int sampled_size = sampled_sizes[(*chr_it).first];
			for (int i=0;i<sampled_size-1;i++){
				if(peaks[(*chr_it).first][i]==0 && peaks[(*chr_it).first][i+1] == 0)
					transition_prob[0][0]+=1;
				else if(peaks[(*chr_it).first][i]==0 && peaks[(*chr_it).first][i+1] == 1)
					transition_prob[0][1]+=1;
				else if(peaks[(*chr_it).first][i]==1 && peaks[(*chr_it).first][i+1] == 0)
					transition_prob[1][0]+=1;
				else if(peaks[(*chr_it).first][i]==1 && peaks[(*chr_it).first][i+1] == 1)
					transition_prob[1][1]+=1;
			}

	}

	for(int i=0;i<2;i++)
		for(int j=0;j<2;j++){
			vector<float> row;
			row.push_back(i);
			row.push_back(j);
			row.push_back(transition_prob[i][j]);
			transition.push_back(row);
		}
}






void Profiler::generate_observations(){
	map<string,float *>::iterator chr_it;
	observation_seq temp;


	for(chr_it=sampled_target_density.begin();chr_it!=sampled_target_density.end();++chr_it){
		temp.chr = (*chr_it).first;
		bool new_line=true;
		for(int i=0;i<sampled_sizes[(*chr_it).first];i++){
			if ((*chr_it).second[i]!=-1 && sampled_control_density[(*chr_it).first][i]!=-1){
				new_line = false;
				int obs = round((*chr_it).second[i]-sampled_control_density[(*chr_it).first][i]); //difference
				temp.values.push_back(obs);

			}
			else{
				if (!new_line){
					temp.start = i-temp.values.size();
					new_line=true;
					if (temp.values.size()>0){
						observations.push_back(temp);
						temp.values.clear();

					}

				}
			}
	}
		temp.start = sampled_sizes[(*chr_it).first]-temp.values.size();
		new_line=true;
		if (temp.values.size()>0){
			observations.push_back(temp);
			temp.values.clear();

		}

	}
}

vector<observation_seq>  Profiler::get_obs(){
	return observations;
}



void Profiler::call_peaks(){


	map<string,float*>::iterator chr_it;

	int min_dist = merge_dist/bin_length+1;
	for(chr_it=sampled_target_density.begin();chr_it!=sampled_target_density.end();++chr_it){
		int sampled_size = sampled_sizes[(*chr_it).first];
		for(int i =0;i<sampled_size;i++)
				peaks[(*chr_it).first][i]=0;
		for(int i =0;i<sampled_size;i++)
			if((*chr_it).second[i]!=-1){
				float max_dense = sampled_control_density[(*chr_it).first][i]>0 ? sampled_control_density[(*chr_it).first][i]:1;
				float pvalue=poissoncdistribution(round((*chr_it).second[i])-1, max_dense);
				if (pvalue<pvalue_threshold)
					peaks[(*chr_it).first][i]=1;
			}
			else
				peaks[(*chr_it).first][i]=-1;

	}



	//remove single noise points
	for(chr_it=sampled_target_density.begin();chr_it!=sampled_target_density.end();++chr_it){
		int sampled_size = sampled_sizes[(*chr_it).first];
				for(int i =2;i<sampled_size-2;i++)
					if(peaks[(*chr_it).first][i]==1 && peaks[(*chr_it).first][i-1]==0  &&
							peaks[(*chr_it).first][i+1]==0)// && peaks[(*chr_it).first][i+2]==0)&& peaks[(*chr_it).first][i-2]==0
						peaks[(*chr_it).first][i]=0;
	}




	for(chr_it=peaks.begin();chr_it!=peaks.end();++chr_it){
		int i=0;
		while(i<sampled_sizes[(*chr_it).first]-1){
			if((*chr_it).second[i]==1 && (*chr_it).second[i+1]==0){
				int count =0;
				int point =i+1;
				while((*chr_it).second[point]==0 && point<sampled_sizes[(*chr_it).first]){
					count++;
					point++;
				}
				if (count<=min_dist){
					for(int j=i+1;j<point;j++)
						peaks[(*chr_it).first][j]=1;
				}
				i=point;
			}
			else
				i++;

		}


	}



}


void Profiler::print_wig(string name){
    map<string,float*>::iterator chr_it;
    ofstream wig_file;
    ofstream pvalue_wig;
    wig_file.open((name+".wig").c_str());

    if(!wig_file){
    	cerr<<"Error: can not open WIG file for writing"<<endl;
    	exit(1);
    }
    wig_file<<"track name="<<"\""<<name<<"\" type=wig visibility=2"<<endl;
    for(vector<observation_seq>::iterator ii=observations.begin();ii!=observations.end();++ii){
    	wig_file<<"fixedStep chrom="<<(*ii).chr<<" start="<<(*ii).start*bin_length+bin_length<<" step="<<bin_length<<endl;

    	for(unsigned int i=0;i<(*ii).values.size();i++){
    		wig_file<<(*ii).values[i]<<endl;
    	}
    	wig_file<<endl;
    }
}


inline void Profiler::count_emissions(){
	int range,peak_index,obs_index;
	string chr;


	//restrict  observations to be in some specific range

	for (unsigned int i=0;i<observations.size();i++){
			for (unsigned int j=0;j<observations[i].values.size();j++){
				if (observations[i].values[j]<-10)
					observations[i].values[j] = -10;
				//if (observations[i].values[j]>100)
				//	observations[i].values[j] = 100;
			}
	}


	//get max and min obs
	for (unsigned int i=0;i<observations.size();i++){
		for (unsigned int j=0;j<observations[i].values.size();j++){
			if (observations[i].values[j]>max_obs   )
				max_obs = observations[i].values[j];
			if (observations[i].values[j]<min_obs)
				min_obs = observations[i].values[j];
		}
	}





	range = max_obs-min_obs+1;
	emission.push_back(vector<float>(range,0));
	emission.push_back(vector<float>(range,0));


	for (unsigned int i=0;i<observations.size();i++){
		chr = observations[i].chr;
		peak_index = observations[i].start;
		for (unsigned int j=0;j<observations[i].values.size();j++){
			obs_index = observations[i].values[j]-min_obs;
			if (peaks[chr][obs_index] == 1)
				emission[1][obs_index]++;
			else if (peaks[chr][peak_index]==0){
				emission[0][obs_index]++;
			}
			peak_index++;
		}
	}
	/*
	if (peaks[chr][position] == 1)
		emission[1][index]++;
	else
		emission[0][index]++;
	*/



}

int Profiler::get_minObs(){
	return min_obs;
}

vector<vector<float> > Profiler:: get_transition(){
	return transition;
}
vector<vector<float> > Profiler:: get_emission(){
	return emission;
}

void Profiler::clean_reads(std::map<std::string,std::vector<DNA_fragment > >& data,
			std::map<std::string,std::vector<DNA_fragment > >& control){

	map<string,vector<DNA_fragment > >::iterator chr_it;
	set<string> ToDelete;

	for (chr_it=data.begin();chr_it!=data.end();++chr_it){
		if(GC_profile.find((*chr_it).first) == GC_profile.end()){
			ToDelete.insert((*chr_it).first);
			cerr<<"Warning: "<<(*chr_it).first<<
					" can not be found on GC index, it will be ignored in further analysis"<<endl;
		}
	}

	for (chr_it=control.begin();chr_it!=control.end();++chr_it){
		if(GC_profile.find((*chr_it).first) == GC_profile.end()){
			ToDelete.insert((*chr_it).first);
			cerr<<"Warning: "<<(*chr_it).first<<
					" can not be found on GC index, it will be ignored in further analysis"<<endl;
		}
	}

	for (set<string>::iterator it=ToDelete.begin();it!=ToDelete.end();++it){
		data.erase((*it));
		control.erase((*it));
	}

	ToDelete.clear();
	for(chr_it=control.begin();chr_it!=control.end();++chr_it)
		if(data.find((*chr_it).first)==data.end()){
			ToDelete.insert((*chr_it).first);
	}

	for (set<string>::iterator it=ToDelete.begin();it!=ToDelete.end();++it){
			control.erase((*it));
	}

	return;
}



void Profiler::remove_blacklist(vector<BedEntry>& blacklist, std::map<std::string,float*>& profile){
	vector<BedEntry>::iterator chr_it;
	for (chr_it = blacklist.begin();chr_it!=blacklist.end();++chr_it){
		string chr = (*chr_it).chr;
		if(profile.find(chr) == profile.end())
			continue;
		int start_point = (*chr_it).start;
		while((start_point+1)%bin_length!=0)
			start_point++;

		for (int i= start_point;i<(*chr_it).end;i+=bin_length){
			int index = (i+1)/bin_length-1;
			if(index >= sampled_sizes[chr]){
				cerr<<"Blacklisted Value out of chromosome bounds. Exiting..."<<endl;
				exit(1);
			}
			profile[chr][index] = -1;

		}

	}
}

