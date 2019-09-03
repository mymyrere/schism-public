FROM ubuntu:14.04

RUN apt-get -y update &&\
    apt-get install -qy \
        gfortran \
        git \
        wget \
        ipython \
        libffi-dev \
        libfreetype6-dev \
        libgeos-3.4.2 \
        libgeos-dev \
        libnetcdf-dev \
        libpng12-dev \
        libssl-dev \
        pkg-config \
        python-dev \
        python-gdal \
        python-netcdf \
        python-nose \
        python-numpy \
        python-pip \
        python-pyproj \
        python-pyside \
        python-pysqlite2 \
        python-scipy \
        python-sqlalchemy \
        python-qt4 \
        python-yaml \
        python-matplotlib \
        python-mpltoolkits.basemap \
        texinfo \
        subversion

RUN pip install numexpr

#------------- spatialindex
RUN  cd /tmp  &&\
     wget http://download.osgeo.org/libspatialindex/spatialindex-src-1.8.0.tar.gz &&\
     tar -xvf spatialindex-src-1.8.0.tar.gz &&\
     cd spatialindex-src-1.8.0 &&\
     ./configure --prefix=/usr/local &&\
     make &&\
     sudo make install &&\
     cd /tmp &&\
     rm -rf /tmp/spatialindex-src-1.8.0.tar.gz /tmp/spatialindex-src-1.8.0  &&\
     cp /usr/local/lib/libspatialindex_c.so.3 /usr/lib/ &&\
     ldconfig

#------------- rtree
RUN cd /tmp &&\
    wget https://pypi.python.org/packages/16/68/c79bed9d7fa4479b37d7a898b71d1e89099a8c15fac5c5d94328cb0a1ac9/Rtree-0.8.2.tar.gz &&\
    tar -xvf Rtree-0.8.2.tar.gz &&\
    cd Rtree-0.8.2 &&\
    sudo python setup.py install &&\
    cd /tmp &&\
    rm -rf /tmp/Rtree-0.8.2.tar.gz /tmp/Rtree-0.8.2
    
USER $USER
ENV HOME /home/$USER
    