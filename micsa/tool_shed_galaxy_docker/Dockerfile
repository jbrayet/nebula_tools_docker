# Dockerfile
#
# Version:          1.0
# Software:         MICSA
# Software Version: 1.0
# Description:      Motif Identification for Chip-Seq Analysis.
# Website:          http://nebula.curie.fr
# Tags:             NGS: Peak Calling
# Provides:         Nebula
# Build Cmd:        docker build --rm -t institutcuriengsintegration/micsa:1.0 .
# Pull Cmd:         docker pull institutcuriengsintegration/micsa:1.0
# Run Cmd:          docker run --rm -it institutcuriengsintegration/micsa:1.0 <options> <files>
#################################################################

# Set the base image to Ubuntu and ruby 2.0
FROM ubuntu:14.04
FROM ruby:2.0

#############

RUN apt-get update && \
    apt-get install -y \
      build-essential \
      wget \
      samtools \
      default-jre \
      csh && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir /usr/bin/micsa
WORKDIR /usr/bin/micsa
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_sources/MICSA.tar.gz -P /usr/bin/micsa
RUN tar -zxvf MICSA.tar.gz
RUN rm -f MICSA.tar.gz
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_sources/ChIPMunk_6.0a.tar.gz -P /usr/bin/micsa
RUN tar -zxvf ChIPMunk_6.0a.tar.gz
RUN rm -f ChIPMunk_6.0a.tar.gz
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_sources/VancouverPackage-4.0.9.2.tar.gz -P /usr/bin/micsa
RUN tar -zxvf VancouverPackage-4.0.9.2.tar.gz
RUN rm -f VancouverPackage-4.0.9.2.tar.gz
RUN wget http://meme-suite.org/meme-software/4.9.1/meme_4.9.1_2.tar.gz -P /usr/bin/micsa
RUN tar -zxvf meme_4.9.1_2.tar.gz
RUN rm -f meme_4.9.1_2.tar.gz
RUN mkdir /usr/bin/meme
WORKDIR /usr/bin/micsa/meme_4.9.1
RUN ./configure --prefix=/usr/bin/meme
RUN make
RUN make install
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_data/createAnnotation/create_annotation_files.sh -P /usr/bin
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_data/createAnnotation/splitChr.pl -P /usr/bin
WORKDIR /usr/bin/micsa

##################### INSTALLATION END #####################

# File Author / Maintainer
MAINTAINER Jocelyn Brayet <jocelyn.brayet@curie.fr>
