FROM centos:centos7

MAINTAINER Bertrand Neron <bneron@pasteur.fr>

USER root

RUN yum clean all && \
    yum install -y epel-release &&\
    yum install -y make gcc gcc-c++ python34-devel zlib-devel python34-tkinter && \
    yum clean all

WORKDIR /tmp
RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py &&\
    python3 get-pip.py
RUN pip3 install psutil>=4.0.0
RUN pip3 install pysam==0.9.1.4
RUN pip3 install numpy==1.11.2
RUN pip3 install cython && pip3 install pandas==0.17.1
RUN pip3 install matplotlib==1.5.3
RUN pip3 install pillow==3.4.2

CMD ["/bin/bash"]

