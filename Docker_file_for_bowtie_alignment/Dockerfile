FROM nfcore/base:1.14
RUN apt-get update
ENV DEBIAN_FRONTEND noninteractive
RUN apt install sudo -y
RUN apt install curl -y
RUN apt-get install zip -y
RUN sudo apt install gzip
RUN mkdir data
RUN mkdir Results
COPY /index_file.1.bt2.gz /data/
COPY /index_file.2.bt2.gz /data/
COPY /index_file.3.bt2.gz /data/
COPY /index_file.4.bt2.gz /data/
COPY /index_file.rev.1.bt2.gz  /data/
COPY /index_file.rev.2.bt2.gz /data/
COPY /environment.yml /environment.yml
RUN conda env create --quiet -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/bowtie/bin:$PATH            #Add conda installation dir to PATH 
RUN conda env export --name bowtie > /Results/bowtie.yml


