# Base Image
FROM ubuntu:22.04

LABEL MAINTAINER="VDJServer <vdjserver@utsouthwestern.edu>"

# Install OS Dependencies
RUN DEBIAN_FRONTEND='noninteractive' apt-get update && DEBIAN_FRONTEND='noninteractive' apt-get install -y \
    make \
    gcc g++ \
    sendmail-bin \
    supervisor \
    wget \
    curl \
    xz-utils \
    python3 \
    python3-pip \
    libyaml-dev

RUN pip3 install \
    'connexion[swagger-ui,uvicorn]>=3.1.0'

##################
##################

# setup vdj user
RUN echo "vdj:x:816290:803419:VDJServer,,,:/home/vdj:/bin/bash" >> /etc/passwd
RUN echo "G-803419:x:803419:vdj" >> /etc/group
RUN mkdir /home/vdj
RUN chown vdj /home/vdj
RUN chgrp G-803419 /home/vdj

# Setup supervisor
COPY docker/start_supervisor.sh /root/start_supervisor.sh
COPY docker/supervisor.conf /etc/supervisor/conf.d/

# image server with connexion
COPY docker/image_server.py /root/image_server.py

##################
##################

# Copy project source
RUN mkdir /repertoire_zoo
COPY . /repertoire_zoo

# build and install
RUN cd /repertoire_zoo && pip3 install .

WORKDIR /root
CMD ["/root/start_supervisor.sh"]
