FROM centos:centos7

MAINTAINER Bertrand Neron <bneron@pasteur.fr>

USER root

RUN yum clean all && \
    yum install -y epel-release &&\
    yum install -y make gcc python34-devel zlib-devel && \
    yum clean all

WORKDIR /tmp
RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py &&\
    python3 get-pip.py &&\
    pip3 install pysam==0.9.1.4

CMD ["/bin/bash"]