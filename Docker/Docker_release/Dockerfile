FROM registry-gitlab.pasteur.fr/bneron/craw:centos-dep-0.6

MAINTAINER Bertrand Neron <bneron@pasteur.fr>

USER root

RUN pip3 install sphinx
RUN yum install -y texlive-collection-latexrecommended

CMD ["/bin/bash"]

