# Dockerfile
#
# Version:          1
# Software:         MakeTSSdist_histones
# Software Version: 1.0
# Description:      This tool creates a .png file with distribution of peaks around gene TSS
# Website:          http://nebula.curie.fr
# Tags:             NGS: Peak Annotation
# Provides:         Nebula
# Build Cmd:        docker build -t institutcuriengsintegration/maketssdisthistones:1.0 -f Dockerfile .
# Pull Cmd:         docker pull institutcuriengsintegration/maketssdisthistones:1.0
# Run Cmd:          docker run --rm -it institutcuriengsintegration/maketssdisthistones:1.0 <options> <files>
#################################################################

# Set the base image to Ubuntu and R 3.1.2
FROM ubuntu:14.04
FROM r-base:3.1.2

#############

RUN apt-get update && \
    apt-get install -y \
      wget \
      curl \
      build-essential \
      gfortran \
      libxt-dev \
      xorg-dev \
      libpango1.0-dev \
      libcairo-dev && \   
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir /usr/bin/maketssdist_histones
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/maketssdist_histones/tool_shed_galaxy_script.tar.gz -P /usr/bin/maketssdist_histones
WORKDIR /usr/bin/maketssdist_histones
RUN tar -xvf tool_shed_galaxy_script.tar.gz && rm tool_shed_galaxy_script.tar.gz
RUN mv tool_shed_galaxy_script/* . && rm -rf tool_shed_galaxy_script makeTSSdist_histones_wrapper* tool_dependencies.xml

RUN echo "install.packages(\"Hmisc\", repos=\"http://cran.rstudio.com/\", dependencies=TRUE)" > install_Hmisc.R
RUN Rscript --slave install_Hmisc.R
RUN rm install_Hmisc.R

WORKDIR /usr/bin/maketssdist_histones

##################### INSTALLATION END #####################

# File Author / Maintainer
MAINTAINER Jocelyn Brayet <jocelyn.brayet@curie.fr>
