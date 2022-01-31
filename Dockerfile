FROM ubuntu:18.04

MAINTAINER Eric Goubault and Sylvie Putot

RUN apt-get -y update
RUN apt-get install build-essential -y --no-install-recommends
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get install -y git  libblas-dev liblapack-dev libgmp3-dev libmpfr-dev libmpfr-doc libmpfr6 gsl-bin libgsl0-dev bison flex gnuplot-x11 libglpk-dev libyaml-cpp-dev m4 python3-pip emacs autoconf automake libtool libmpfi-dev curl wget 
RUN apt install -y python3-dev
RUN apt-get -y update
RUN apt install -y python3-matplotlib
RUN pip3 install pip --upgrade

# protobuf-compiler
#COPY ./protobuf-3.15.3 /home/protobuf
#RUN cd /home/protobuf && ./autogen.sh && ./configure && make && make install && ldconfig

RUN cd /home && wget http://www.fadbad.com/download/FADBAD++-2.1.tar.gz && tar -xzvf FADBAD++-2.1.tar.gz && rm FADBAD++-2.1.tar.gz
#COPY ./FADBAD++ /home/FADBAD++

RUN cd /home && curl -O http://www2.math.uni-wuppertal.de/wrswt/software/filib++/filibsrc-3.0.2.tar.gz && tar -xzvf filibsrc-3.0.2.tar.gz && rm filibsrc-3.0.2.tar.gz
#COPY ./filibsrc /home/filibsrc

COPY . /home/RINO

RUN cd /home/filibsrc && ./configure && make && make install

#RUN cd /home/RINO/sherlock_2_reduct && make

RUN touch /root/.bashrc.sh
RUN echo "#!/bin/bash" >> /root/.bashrc.sh
RUN echo "export PATH=/usr/local/include:${PATH}" >> /root/.bashrc.sh
RUN echo "export FADBADHOME=/home/FADBAD++" >> /root/.bashrc.sh
RUN echo "export FILIBHOME=/home/filibsrc" >> /root/.bashrc.sh
RUN chmod a+x /root/.bashrc.sh

RUN cd /home/RINO/aaflib-0.1 && make static

RUN . /root/.bashrc.sh && cd /home/RINO && for f in *.cpp; do iconv -c -f UTF-8 -t US-ASCII $f -o $f; done && for f in *.h; do iconv -c -f UTF-8 -t US-ASCII $f -o $f; done && make
