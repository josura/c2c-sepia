
FROM ubuntu:latest


RUN mkdir /opt/sepia
WORKDIR /opt/sepia

RUN apt-get update 

RUN apt-get -y install build-essential

RUN apt-get -y install make cmake 
RUN apt-get -y install libopenblas-dev liblapack-dev

RUN apt-get -y install wget

RUN apt-get -y install libboost-all-dev

RUN wget https://sourceforge.net/projects/arma/files/armadillo-12.0.1.tar.xz

RUN apt-get -y install 

ADD ./src ./src

ADD CMakeLists.txt CMakeLists.txt

ADD ./resources ./resources

ADD ./data ./data

ADD ./buildarmadillo.sh ./buildarmadillo.sh

RUN sh buildarmadillo.sh

CMD [ "bash" ]

#ENTRYPOINT [ "sbt" , "run" ]

# ENTRYPOINT [ "/bin/sh" , "entrypoint.sh" , $*]