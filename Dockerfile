FROM ubuntu:xenial
## install dependencies
RUN apt-get update \
 && apt-get -y install software-properties-common \
 && add-apt-repository ppa:remik-ziemlinski/nccmp --update \
 && apt-get install -y libnetcdf-dev libnetcdff-dev netcdf-bin gfortran bats nccmp autoconf
## copy repo into container
COPY . .
## run tests
ENTRYPOINT ["build.sh"]
