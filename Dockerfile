FROM maven:3.3.9-jdk-8

MAINTAINER Myron Peto <peto@ohsu.edu>

USER root
ENV PATH /opt/bin:$PATH

RUN apt-get update && \
    apt-get install --yes \
    build-essential \
    git \
    python \
    && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /opt

COPY ./resource /opt/resource

# RUN git clone https://github.com/PathwayAndDataAnalysis/mutex.git && \
RUN cd resource && \
    mvn clean compile && \
    mvn assembly:single && \
    mv target/maf2matrix.jar /home/

COPY ./maf2matrix.py /opt/bin/
COPY ./0024_filter.maf /home/
RUN chmod +x /opt/bin/maf2matrix.py
RUN chmod +x /home/maf2matrix.jar
ENV PATH=$PATH:/opt/

LABEL version="1.0" description="datamatrix"

WORKDIR /home/
# CMD ["matrix.py"]
