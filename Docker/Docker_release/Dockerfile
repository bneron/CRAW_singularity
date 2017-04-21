FROM registry-gitlab.pasteur.fr/bneron/craw:centos-run-dep

MAINTAINER Bertrand Neron <bneron@pasteur.fr>

USER root

RUN pip3 install sphinx
RUN yum install -y texlive-collection-latexrecommended \
                   texlive-titlesec \
                   texlive-framed \
                   texlive-threeparttable \
                   texlive-wrapfig \
                   texlive-multirow \
                   texlive-upquote

COPY tabulary.sty fncychap.sty capt-of.sty \
     eqparbox.sty environ.sty trimspaces.sty \
     /usr/lib64/python3.4/site-packages/sphinx/texinputs/

CMD ["/bin/bash"]

