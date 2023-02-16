#!/bin/bash

set -e
set -x

if ! command -v python > /dev/null 2>&1; then
    echo "can not find python"
    exit 1
fi

pyversion='^Python 3.8.[0-9]*$'

if [[ ! "$(python --version)" =~ ${pyversion} ]]; then
    echo "we need python 3.8"
    exit 1 
fi

pip3_install()
{
    pip3 install --user joblib==1.0.1
    pip3 install --user pandas==1.4.1
    pip3 install --user scikit-learn==0.24.2
    pip3 install --user xgboost==1.6.1
    pip3 install --user pyarrow
}

conda_install()
{
    echo "Use const to install required"
    conda install -y joblib=1.0.1
    conda install -y pandas=1.4.1
    conda install -y scikit-learn=0.24.2
    conda install -y pyarrow
    pip3 install --user xgboost==1.6.1
}

htslib_install()
{
    echo "Install htslib"
    git clone git@github.com:samtools/htslib.git
    cd htslib/
    git submodule update --init --recursive
    mkdir build
    autoreconf -i
    ./configure --prefix=$(pwd)/build --disable-lzma
    make
    make install
}

if command -v conda > /dev/null 2>&1; then
    is_active=$(conda info | sed -n "2p" | cut -d":" -f2)
    if [ "${is_active}" = " None" ]; then
        pip3_install
    else
        conda_install
    fi
else
    pip3_install
fi

mkdir build
cd build
echo "Enter build"

htslib_install

cd ..

cmake -DHTS_PREFIX=$(pwd)/htslib/build -DCMAKE_INSTALL_PREFIX=$(pwd) ..
make -j4
make install

echo "Install success"
