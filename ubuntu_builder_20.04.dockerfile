FROM ubuntu:20.04
ENV DEBIAN_FRONTEND noninteractive
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get -qq update \
&&  apt-get -qq --no-install-recommends install software-properties-common \
&&  add-apt-repository -y ppa:deadsnakes \
&&  add-apt-repository -y ppa:mhier/libboost-latest \
&&  apt-get -qq update

RUN apt-get -qq --no-install-recommends install libpython3.7 python3.7 python3-pip libpython3.7-dev pybind11-dev \
&&  apt-get -qq --no-install-recommends install qtbase5-dev libqt5charts5-dev libqt5opengl5-dev qtwebengine5-dev libopengl0 libeigen3-dev libglew-dev zlib1g-dev \
&&  apt-get -qq --no-install-recommends install libboost1.67-dev libboost-system1.67-dev libboost-filesystem1.67-dev libboost-program-options1.67-dev libboost-thread1.67-dev \
&&  apt-get -qq --no-install-recommends install libmkl-dev libmkl-interface-dev libmkl-threading-dev libmkl-computational-dev \
&&  apt-get -qq --no-install-recommends install g++ cmake ccache ninja-build curl unzip git \
&&  python3.7 -m pip install numpy \
&&  apt-get clean


WORKDIR /opt
ENV SOFA_ROOT=/opt/sofa
ENV PYTHONPATH=/opt/sofa/plugins/SofaPython3/lib/python3/site-packages