# Dockerfile
#
# Version:          1.0
# Software:         Filter Peaks
# Software Version: 1.0
# Description:      Filters out peaks present both in ChIP and Control dataset.
# Website:          http://nebula.curie.fr
# Tags:             NGS: Peak Calling
# Provides:         Nebula
# Build Cmd:        docker build --rm -t institutcuriengsintegration/filterpeaks:1.0 .
# Pull Cmd:         docker pull institutcuriengsintegration/filterpeaks:1.0
# Run Cmd:          docker run --rm -it institutcuriengsintegration/filterpeaks:1.0 <options> <files>
#################################################################

# Set the base image to Ubuntu
FROM ubuntu:14.04

#############

RUN apt-get update && \
    apt-get install -y \
      build-essential \
      wget \
      samtools \
      default-jre && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir /usr/bin/filterpeaks
WORKDIR /usr/bin/filterpeaks
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_sources/FilterPeaks.class -P /usr/bin/filterpeaks
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/filterpeaks/tool_shed_galaxy_script/peak2bed.pl -P /usr/bin/filterpeaks
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_data/createAnnotation/create_annotation_files.sh -P /usr/bin
WORKDIR /usr/bin/filterpeaks

##################### INSTALLATION END #####################

# File Author / Maintainer
MAINTAINER Jocelyn Brayet <jocelyn.brayet@curie.fr>
