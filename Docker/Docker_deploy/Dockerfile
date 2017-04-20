FROM registry-gitlab.pasteur.fr/bneron/craw:centos-dep-0.6

MAINTAINER Bertrand Neron <bneron@pasteur.fr>

USER root

RUN pip3 install sphinx
RUN yum install -y texlive-collection-latexrecommended
RUN pip3 install coverage

CMD ["/bin/bash"]

