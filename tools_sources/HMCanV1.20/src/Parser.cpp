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

#include "Parser.h"
#include <fstream>
#include <iostream>
#include<algorithm>
using namespace std;
#include <cstdlib>

Parser::Parser() {
	GC_index="";
	path ="";
	blacklist = "";
	min_length=0;
	median_length=0;
	max_length=0;
	large_bin_size=0;
	bin_size=0;
	pvalue_threshold =-1;
	merge_dist = -1;
	iter_threshold = -1;
	final_threshold = -1;
	max_iter = 0;
	posterior_threshold =-1;
	wig = false;
	posterior = false;
	t = BAM;

}

Parser::~Parser() {
	// TODO Auto-generated destructor stub
}

void Parser::parse(char * config_file){

	ifstream conf;
	conf.open(config_file);
	string item,tempStr;
	if (!conf){
		cerr<<"Can not open the file "<<config_file<<endl;
		exit(1);
	}
		//geting the configuration
		while(!conf.eof()){

			conf>>item>>tempStr;
			if (item.compare("GCIndex")==0)
				GC_index = tempStr;
			else if (item.compare("genomePath")==0)
				path = tempStr;
			else if (item.compare("minLength")==0)
				min_length = atoi(tempStr.c_str());
			else if(item.compare("medLength")==0)
				median_length = atoi(tempStr.c_str());
			else if (item.compare("maxLength")==0)
				max_length= atoi(tempStr.c_str());
			else if (item.compare("smallBinLength")==0)
				bin_size = atoi(tempStr.c_str());
			else if (item.compare("largeBinLength")==0)
				large_bin_size = atoi(tempStr.c_str());
			else if (item.compare("pvalueThreshold")==0)
				pvalue_threshold = atof(tempStr.c_str());
			else if (item.compare("mergeDistance")==0)
				merge_dist = atoi(tempStr.c_str());
			else if(item.compare("iterationThreshold")==0)
				iter_threshold = atof(tempStr.c_str());
			else if (item.compare("finalThreshold")==0)
				final_threshold = atof(tempStr.c_str());
			else if (item.compare("maxIter")==0)
				max_iter = atoi(tempStr.c_str());
			else if (item.compare("format")==0){
				transform(tempStr.begin(), tempStr.end(),tempStr.begin(), ::tolower);
				if (tempStr.compare("bed")==0)
					t = BED;
				else if (tempStr.compare("sam")==0)
					t = SAM;
				else if (tempStr.compare("bam")==0)
					t = BAM;
				else{
					cerr<<"Please provide a file in SAM or BED format."<<endl;
					exit(EXIT_FAILURE);
				}
			}
			else if (item.compare("PrintWig")==0){
				transform(tempStr.begin(), tempStr.end(),tempStr.begin(), ::tolower);
				if (tempStr.compare("true")==0)
					wig = true;
				else if (tempStr.compare("false")==0)
					wig = false;
				else{
					wig = false;
					cerr<<"Warning: "<<item<<" is not a valid option for PrintWig; it has to be TRUE or FALSE"<<endl;
					cerr<<"Warning: PrintWig is set to FALSE"<<endl;
				}
			}

			else if (item.compare("PrintPosterior")==0){
	                        transform(tempStr.begin(), tempStr.end(),tempStr.begin(), ::tolower);
	                        if (tempStr.compare("true")==0)
	                                posterior = true;
	                        else if (tempStr.compare("false")==0)
	                                posterior = false;
	                        else{
	                                posterior = false;
	                                cerr<<"Warning: "<<item<<" is not a valid option for printPosterior; it has to be TRUE or FALSE"<<endl;
	                                cerr<<"Warning: printPosterior is set to FALSE"<<endl;
	                        }
	                }
			else if (item.compare("PosteriorProb")==0)
				posterior_threshold = atof(tempStr.c_str());
			else if(item.compare("blackListFile")== 0)
				blacklist = tempStr;
			else
				cerr<<"Warning: "<<item<<" is not a parameter for HMCan"<<endl;
			}

		// validating parameters
		if (GC_index=="")
		{
			cerr<<"Please provide a valid file name for GC content information file"<<endl;
			exit(EXIT_FAILURE);
		}

		if (path=="")
		{
			cerr<<"Please provide a valid path for your target genome"<<endl;
			exit(EXIT_FAILURE);
		}

		if(min_length<1){
			min_length = 145;
			cerr<<"Warning: Invalid parameter, minimum fragment length is set to its default value 145"<<endl;
		}

		if(median_length<1){
			median_length = 150;
			cerr<<"Warning: Invalid parameter, median fragment length is set to its default value 150"<<endl;
		}

		if(max_length<1){
			max_length = 155;
			cerr<<"Warning: Invalid or missed parameter, maximum fragment length is set to its default value 155"<<endl;
		}

		if(bin_size<1){
			bin_size = 50;
			cerr<<"Warning: Invalid or missed parameter, bin length length is set to its default value 50"<<endl;
		}

		if(large_bin_size<1){
			large_bin_size = 100000;
			cerr<<"Warning: Invalid or missed parameter, large bin size length is set to its default value 100000"<<endl;
		}

		if(pvalue_threshold<0 || pvalue_threshold>1){
			pvalue_threshold = 0.01;
			cerr<<"Warning: Invalid or missed parameter, P-value threshold is set to its default value 0.01"<<endl;
		}

		if(merge_dist<0){
			 merge_dist= 200;
			cerr<<"Warning: Invalid or missed parameter, merge distance is set to its default value 200"<<endl;
		}

		if(iter_threshold<0){
			iter_threshold = 5;
			cerr<<"Warning: Invalid or missed parameter, iterations score threshold is set to its default value 5"<<endl;
		}

		if(final_threshold<0){
			final_threshold = 0;
			cerr<<"Warning: Invalid or missed parameter, final score threshold is set to its default value 0"<<endl;
		}

		if (max_iter<1){
			max_iter = 10;
			cerr<<"Warning: Invalid or missed parameter, maximum iterations is set to its default value 10"<<endl;
		}
		if(posterior_threshold == -1){
			posterior_threshold = 0.7;
			cerr<<"Warning: Invalid or missed parameter, PosteriorProb is set to its default value 0.7"<<endl;
		}



}


void Parser::print(){
	cout<<GC_index<<endl;
	cout<<path<<endl;
	cout<<min_length<<endl;
	cout<<median_length<<endl;
	cout<<max_length<<endl;
	cout<<large_bin_size<<endl;
	cout<<bin_size<<endl;
	cout<<pvalue_threshold<<endl;
	cout<<merge_dist<<endl;
	cout<<iter_threshold <<endl;
	cout<<final_threshold <<endl;
	cout<<max_iter<<endl;
	cout<<blacklist<<endl;
}
