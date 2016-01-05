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

#include "HMCan.h"
#include <fstream>
#include <iomanip>
#include <map>
#include <cmath>
#include <algorithm>
#include <ctime>
using namespace std;
HMCan::HMCan(vector<observation_seq>& obs_, vector<vector<float> >& transitions,vector<vector<float> >& emissions,
			float threshold1,float threshold2, int step,int min_obs, int bin_size, float posterior_threshold): obs(obs_){
	this->transitions = transitions;
	this->emissions = emissions;
	this->threshold1 = threshold1;
	this->threshold2 = threshold2;
	this->step = step;
	this->min_obs = min_obs;
	this->bin_size = bin_size;
	signal =0;
	background =0;
	this->posterior_threshold = posterior_threshold;
}

HMCan::~HMCan() {
	// TODO Auto-generated destructor stub
}

int HMCan::myrandom (int i) { return rand()%i;}
void HMCan:: estimate_hmm_parameters(vector<vector<int> >& raw_peaks, vector<int>& non_peaks,int total_bins){
	map<int,float> emsision_peak,emision_non_peak;
	map<int,float>::iterator key;
	int total_peaks=0, peak_bins=0;
	float transition_prob[2][2];
	transition_prob[0][0]=0;
	transition_prob[0][1]=0;
	transition_prob[1][0]=0;
	transition_prob[1][1]=0;

	for (unsigned int i=0;i<peaks.size();i++){
		if(peaks[i].score > threshold1){
			peak_bins+=raw_peaks[i].size();
			total_peaks++;
		}
		else{

			for(unsigned int j=0;j<raw_peaks[i].size();j++){
				non_peaks.push_back(raw_peaks[i][j]);
			}
		}
	}


	/* transition re estimation*/
	transition_prob[0][1] = transition_prob[1][0]= total_peaks;
	transition_prob[1][1] = peak_bins-total_peaks;
	transition_prob[0][0] = total_bins - transition_prob[0][1]-transition_prob[0][1]-transition_prob[1][1];
	transitions.clear();
	for(int i=0;i<2;i++)
		for(int j=0;j<2;j++){
			vector<float> row;
			row.push_back(i);
			row.push_back(j);
			row.push_back(transition_prob[i][j]);
			transitions.push_back(row);
	}

	/*re estimation of emission*/

	//clear emissions
	for(int i=0;i<2;i++)
		for(unsigned int j=0;j<emissions[i].size();j++)
			emissions[i][j] = 0;

		for(unsigned int i=0;i<raw_peaks.size();i++){
			if (peaks[i].score>threshold1)
				for(unsigned int j=0;j<raw_peaks[i].size();j++){

					int index = raw_peaks[i][j]-min_obs;
					//cout<<raw_peaks[i][j]<<" "<<emissions[0].size()<<endl;
					emissions[1][index]++;

			}
		}


	for (unsigned int i=0;i<non_peaks.size();i++){
		int index = non_peaks[i]-min_obs;
		emissions[0][index]++;

	}
}


void HMCan:: merge_peaks(int dist){
	map<string,vector<BedEntry> > chrs;
	BedEntry temp;
	map<string,vector<BedEntry> >::iterator chr_it;
	for(vector<BedEntry>::iterator it = peaks.begin();it!=peaks.end();++it)
		chrs[(*it).chr].push_back((*it));
	bool pushed = false;
	for(chr_it=chrs.begin();chr_it!=chrs.end();++chr_it){
		sort((*chr_it).second.begin(),(*chr_it).second.end(),BedStartCompare);
		temp.chr = (*chr_it).first;
		temp.start = (*chr_it).second[0].start;
		temp.end = (*chr_it).second[0].end;
		temp.score = (*chr_it).second[0].score;
		for (unsigned int i=1;i<(*chr_it).second.size();i++){
			if ((*chr_it).second[i].start-temp.end+2<dist){
				temp.end = (*chr_it).second[i].end;
				temp.score = (*chr_it).second[i].score>temp.score ? (*chr_it).second[i].score:temp.score;
				pushed = false;
			}
			else{
				pushed = true;
				regions.push_back(temp);
				temp.start = (*chr_it).second[i].start;
				temp.end = (*chr_it).second[i].end;
				temp.score = (*chr_it).second[i].score;
			}
		}
	}
	if (!pushed)
			regions.push_back(temp);
}



void HMCan::run(int dist, int max_iter){


	vector<vector<int> > raw_peaks;
	vector<int> non_peaks;
	estimate_transitions();
	estimate_emissions();
	//distribution_info();
	
	
	posterior_decoding(raw_peaks,non_peaks);
	int peaks_count=0;
	int total_bins =0;
	for(unsigned int i=0;i<obs.size();i++)
		total_bins+=obs[i].values.size();

	for(unsigned int i=0;i<peaks.size();i++)
		if (peaks[i].score>threshold1)
			peaks_count++;	
	int iterations = 1;

	cout<<"Calling peaks...."<<endl;
	cout<<"Iteration # 1"<<endl;
	int current_count =0;
	do{
		if(iterations>2)
			peaks_count=current_count;
		cout<<"Iteration #"<<iterations+1<<endl;
		sort(peaks.begin(),peaks.end(),BedCompare);
		estimate_hmm_parameters(raw_peaks, non_peaks,total_bins);
		estimate_transitions();
		estimate_emissions();
		//distribution_info();
		peaks.clear();
		non_peaks.clear();
		raw_peaks.clear();
		posterior_decoding(raw_peaks,non_peaks);
		current_count=0;
		for(unsigned int i=0;i<peaks.size();i++)
			if (peaks[i].score>threshold1){
				current_count++;
			}
		iterations++;
		
	}while( iterations<max_iter && abs(peaks_count-current_count)>100);
	cout<<"Done!"<<endl;
//	posterior_decoding(raw_peaks,non_peaks);

	//print_probabilities(prob_log,20);
	merge_peaks(dist);
	sort(peaks.begin(),peaks.end(),BedCompare);
	sort(regions.begin(),regions.end(),BedCompare);
}


void HMCan::print(string name, int fragment_length){
	ofstream out;
	cout.setf(ios::fixed);
	out.open((name+"_peaks.narrowPeak").c_str());
	for(unsigned int i=0;i<peaks.size();i++)
		if (peaks[i].score>threshold2 && (peaks[i].end - peaks[i].start) >= fragment_length)
			out<<peaks[i].chr<<'\t'<<peaks[i].start<<'\t'<<peaks[i].end<<'\t'<<"peak"<<i+1<<"\t"<<peaks[i].score<<"\t.\t"<<
			peaks[i].max_value<<"\t-1\t-1\t"<<peaks[i].max<<endl;
	out.close();

	out.open((name+"_regions.bed").c_str());
	for(unsigned int i=0;i<regions.size();i++)
			if (regions[i].score>threshold2  && (regions[i].end - regions[i].start) >= fragment_length)
				out<<regions[i].chr<<'\t'<<regions[i].start<<'\t'<<regions[i].end<<'\t'<<"peak"<<i+1<<'\t'<<regions[i].score<<"\t."<<endl;
		out.close();


}

void HMCan::print_probabilities(ofstream& prob_log,int iteration){

	for (unsigned int i=0;i<emissions.size();i++)
		for (int j=0;j<emissions[i].size();j++)
			prob_log<<iteration<<'\t'<<i<<'\t'<<j+min_obs<<'\t'<<emissions[i][j]<<endl;

	return;
}

void HMCan::posterior_decoding(vector<vector<int> >& raw_peaks, vector<int>& non_peaks){
	for(unsigned int i=0;i<obs.size();i++){
		obs[i].posterior_prob = vector<float>(obs[i].values.size(),0);
	}

	for(unsigned int i=0;i<obs.size();i++){
		states = vector<int>(obs[i].values.size(),0);
		forward_backward(obs[i].values,obs[i].posterior_prob);
		get_peaks(obs[i],raw_peaks,non_peaks,i);
	}
}



void HMCan::print_posterior(string name){
	ofstream wig_file;
    wig_file.open((name+"_posterior.wig").c_str());

    wig_file<<"track name="<<"\""<<name<<"\" type=wig visibility=2"<<endl;
    for(vector<observation_seq>::iterator ii=obs.begin();ii!=obs.end();++ii){
    	wig_file<<"fixedStep chrom="<<(*ii).chr<<" start="<<(*ii).start*bin_size+bin_size<<" step="<<bin_size<<endl;

    	for(unsigned int i=0;i<(*ii).posterior_prob.size();i++){
    		float p = (*ii).posterior_prob[i];
    		wig_file<<std::setprecision(6)<<p<<endl;

    	}
    	wig_file<<endl;
    }
}


void HMCan::forward_backward(vector<int>& s_obs, vector<float>& posterior){

	vector<vector<double> > alphas, betas,post;
	double sum,sum1;
	int index, index2;

	for(unsigned int i=0;i<s_obs.size();i++){
		alphas.push_back(vector<double> (2,0));
		betas.push_back(vector<double> (2,0.5));
		post.push_back(vector<double> (2,0));
	}


	index = 0-min_obs;
	alphas[0][0] = background*emissions[0][index];
	alphas[0][1] = signal*emissions[1][index];

	//forward algorithm
	for (unsigned int i=1;i<s_obs.size();i++){
		sum1 = 0;
		for(int j=0;j<2;j++){
			sum = 0;
			for(int k=0;k<2;k++){
				index2 = search_transition(k,j);
				sum+=alphas[i-1][k]*transitions[index2][2];
			}		
			index = s_obs[i]-min_obs;
			alphas[i][j] = emissions[j][index]*sum;
			sum1+=alphas[i][j];
		}
		for (int j=0;j<2;j++)
			alphas[i][j]/=sum1;


	}		
	
	//backward algorithm
	for(int i=s_obs.size()-2;i>=0;i--){
		sum1=0;
		for(int j=0;j<2;j++){
			sum=0;
			for (int k=0;k<2;k++){
				index = s_obs[i+1]-min_obs;
				index2 = search_transition(j,k);
				sum+=betas[i+1][k]*transitions[index2][2]*emissions[k][index];
			}
			betas[i][j] = sum;
			sum1+=betas[i][j];
		}
	
		for (int j=0;j<2;j++)
			betas[i][j]/=sum1;
	}


	//calculate posterior probability
	for(unsigned int i=0;i<s_obs.size();i++){
		sum1= 0;
		for(int j=0;j<2;j++){
			post[i][j] = alphas[i][j]*betas[i][j];
			sum1+=post[i][j];
		}

		for(int j=0;j<2;j++)
			post[i][j]/=sum1;
	}

	for(unsigned int i=0;i<s_obs.size();i++){
		if (post[i][1] > 0.99999){
			posterior[i] = 0.99999;
		}
		else if  (post[i][1] < 1e-16)
			posterior[i] = 1e-16;
		else
			posterior[i] = post[i][1];


		if (posterior[i]>=posterior_threshold)
			states[i]=1;
		else
			states[i]=0;
	}

}

int HMCan::search_emission(int state, int obs){

	for (unsigned int i=0;i<emissions.size();i++){
		if (emissions[i][0] == state && emissions[i][1] == obs)
			return i;

	}
	return -1;
}

int HMCan::search_transition(int state1, int state2){
	for(unsigned int i=0;i<transitions.size();i++)
		if (transitions[i][0] == state1 && transitions[i][1] == state2)
			return i;

	return -1;
}



void HMCan::estimate_emissions(){
	int pseudo_count = 1;
	int range = emissions[0].size();
	int smoothing_factor = pseudo_count*range;
	int peaks_sum=0,non_peaks_sum=0;

	for (unsigned int i =0;i<emissions[0].size();i++){
		non_peaks_sum+=emissions[0][i];
		peaks_sum+=emissions[1][i];
	}
		//get emission probabilities
	for (unsigned int i=0;i<emissions[0].size();i++){
		emissions[0][i] = (emissions[0][i]+pseudo_count)/(non_peaks_sum+smoothing_factor);
		emissions[1][i] = (emissions[1][i]+pseudo_count)/(peaks_sum+smoothing_factor);
	}


}


void HMCan::estimate_transitions(){
	float row1=0,row2=0;
	int sum = 0;
	for(unsigned int i=0;i<transitions.size();i++){
			if (transitions[i][0]==0)
				row1+=transitions[i][2];
			if (transitions[i][0]==1)
				row2+=transitions[i][2];
			sum+=transitions[i][2];
	}


	for(unsigned int i=0;i<transitions.size();i++)
		if (transitions[i][0]==0)
			transitions[i][2]=transitions[i][2]/row1;
		else
			transitions[i][2]=transitions[i][2]/row2;
	signal = row2/sum;
	background = row1/sum;

}

void HMCan::get_peaks(observation_seq& s_obs,vector<vector<int> >& raw_peaks,vector<int>& non_peaks, int obs_index){
	BedEntry temp;
	temp.chr = s_obs.chr;
	temp.score = 0;
	int offset = s_obs.start;
	vector<int> peak;
	vector <float> posterior;
	if (states[0]==1){
		temp.start =(0+offset)*bin_size+bin_size-1;
	}

	bool flag = false;
	for(unsigned int i=0;i<states.size()-1;i++){
		if (states[i]==0 && states[i+1]==1){
			temp.start = (i+1+offset)*bin_size+bin_size-1;
			flag = true;
		}
		else if (states[i] == 1 && states[i+1] == 1){
			peak.push_back(s_obs.values[i]);
			posterior.push_back(s_obs.posterior_prob[i]);
		}
		else if (states[i] == 1 && states[i+1]==0){
			peak.push_back(s_obs.values[i]);
			posterior.push_back(s_obs.posterior_prob[i]);
			temp.end = (i+offset)*bin_size+bin_size;
			score_peak(temp,peak,posterior);
			raw_peaks.push_back(peak);
			peak.clear();
			posterior.clear();
			peaks.push_back(temp);
		}

		if (states[i] == 0)
			non_peaks.push_back(s_obs.values[i]);
	}


	if (states[states.size()-1]==1){
		if (states[states.size()-2]==1 && peak.size()>0)
		{
			temp.end = (states.size()-1+offset)*bin_size+bin_size;
			peak.push_back(s_obs.values[states.size()-1]);
			posterior.push_back(s_obs.posterior_prob[states.size()-1]);
			raw_peaks.push_back(peak);
			score_peak(temp,peak,posterior);
			//temp.max_index = temp.start+bin_size*temp.max_index;
			peak.clear();
			posterior.clear();
			peaks.push_back(temp);
		}

		else if (states[states.size()-2]==0 && peak.size()>0){
			temp.end = (states.size()-2+offset)*bin_size+bin_size;
			raw_peaks.push_back(peak);
			score_peak(temp,peak,posterior);
			peak.clear();
			posterior.clear();
			peaks.push_back(temp);
		}
		else
			non_peaks.push_back(s_obs.values[states.size()-1]);
	}





}

void HMCan::score_peak(BedEntry& entry,vector<int>& peak,vector<float>& posterior){
	unsigned int max_index =0;
	int max = peak[0];
	float score=0;
	for(unsigned int i=1;i<peak.size();i++)
	if (peak[i]>max){
		max_index =i;
		max = peak[i];
	}
	int sample_point = max_index;
	while (sample_point>=0){
		score+=log10(posterior[sample_point])-log10(1-posterior[sample_point]);
		sample_point-=step;
	}
	sample_point = max_index+step;

	while(sample_point<peak.size()){
		score+=log10(posterior[sample_point])-log10(1-posterior[sample_point]);
		sample_point+=step;
	}
	entry.score = score;
	entry.max = max_index;
	entry.max_value = max;
}

