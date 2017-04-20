FROM registry-gitlab.pasteur.fr/bneron/craw:centos-run-dep

MAINTAINER Bertrand Neron <bneron@pasteur.fr>

USER root

RUN pip3 install sphinx
RUN yum install -y texlive-collection-latexrecommended

CMD ["/bin/bash"]

