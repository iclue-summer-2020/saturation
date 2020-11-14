FROM ubuntu:16.04

WORKDIR /iclue-summer-2020/saturation

RUN apt-get -y update    \
 && apt-get -y install   \
      bash               \
      build-essential    \
      git                \
      python3            \
      python3-dev        \
      python3-pip        \
      time               \
 && python3 -m pip install --upgrade pip \
 && python3 -m pip install   \
      cmake==3.17.3

COPY . .
