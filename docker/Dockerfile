FROM ubuntu:16.04

ENV CONDA_DIR="/opt/conda"
ENV PATH="$CONDA_DIR/bin:$PATH"

# Install conda
RUN CONDA_VERSION="4.5.11" && \
    CONDA_MD5_CHECKSUM="e1045ee415162f944b6aebfe560b8fee" && \
    \
    apt update && \
    apt install --no-install-recommends -y wget tar ca-certificates bzip2 && \
    apt clean &&\
    rm -rf /var/lib/apt/lists/* && \
    \
    mkdir -p "$CONDA_DIR" && \
    wget "http://repo.continuum.io/miniconda/Miniconda3-${CONDA_VERSION}-Linux-x86_64.sh" -O miniconda.sh && \
    echo "$CONDA_MD5_CHECKSUM  miniconda.sh" | md5sum -c && \
    bash miniconda.sh -f -b -p "$CONDA_DIR" && \
    echo "export PATH=$CONDA_DIR/bin:\$PATH" > /etc/profile.d/conda.sh && \
    rm miniconda.sh && \
    \
    conda update --all --yes && \
    conda config --set auto_update_conda False && \
    conda install nomkl && \
    rm -r "$CONDA_DIR/pkgs/" && \
    \
    \
    mkdir -p "$CONDA_DIR/locks" && \
    chmod 777 "$CONDA_DIR/locks"

# Install our software stack
RUN apt update && \ 
    apt install --no-install-recommends -y git && \
    apt clean &&\
    rm -rf /var/lib/apt/lists/* && \ 
    conda install -y --only-deps -c omnia -c cmelab -c mosdef mbuild foyer python=3.6 && \
    conda install -y -c glotzer freud signac signac-flow gsd && \
    conda install -y -c omnia openmm=7.2 && \
    conda clean -tipsy && \
    pip install --no-cache-dir --upgrade git+https://bitbucket.org/cmelab/cme_utils.git git+https://github.com/mikemhenry/foyer.git@mike/plankton git+https://github.com/mikemhenry/mbuild.git@mike/plankton signac signac-flow

# Install HOOMD
# Runtime
ENV CUDA_VERSION=8.0.61
ENV CUDA_PKG_VERSION=8-0=$CUDA_VERSION-1
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates apt-transport-https gnupg-curl && \
    rm -rf /var/lib/apt/lists/* && \
    NVIDIA_GPGKEY_SUM=d1be581509378368edeec8c1eb2958702feedf3bc3d17011adbf24efacce4ab5 && \
    NVIDIA_GPGKEY_FPR=ae09fe4bbd223a84b2ccfce3f60f4b3d7fa2af80 && \
    \
    apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1604/x86_64/7fa2af80.pub && \
    apt-key adv --export --no-emit-version -a $NVIDIA_GPGKEY_FPR | tail -n +5 > cudasign.pub && \
    echo "$NVIDIA_GPGKEY_SUM  cudasign.pub" | sha256sum -c --strict - && rm cudasign.pub  && \
    echo "deb https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1604/x86_64 /" > /etc/apt/sources.list.d/cuda.list  && \
    echo "deb https://developer.download.nvidia.com/compute/machine-learning/repos/ubuntu1604/x86_64 /" > /etc/apt/sources.list.d/nvidia-ml.list && \
    \
    apt-get update && apt-get install -y --no-install-recommends \
          cuda-nvrtc-$CUDA_PKG_VERSION \
          cuda-nvgraph-$CUDA_PKG_VERSION \
          cuda-cusolver-$CUDA_PKG_VERSION \
          cuda-cublas-8-0=8.0.61.2-1 \
          cuda-cufft-$CUDA_PKG_VERSION \
          cuda-curand-$CUDA_PKG_VERSION \
          cuda-cusparse-$CUDA_PKG_VERSION \
          cuda-npp-$CUDA_PKG_VERSION \
          cuda-cudart-$CUDA_PKG_VERSION && \
    \
    ln -s cuda-8.0 /usr/local/cuda && \
    rm -rf /var/lib/apt/lists/* && \
    \ 
    echo "/usr/local/nvidia/lib" >> /etc/ld.so.conf.d/nvidia.conf && \
    echo "/usr/local/nvidia/lib64" >> /etc/ld.so.conf.d/nvidia.conf && \
    apt clean &&\
    rm -rf /var/lib/apt/lists/*

# Devel
RUN apt-get update && apt-get install -y --no-install-recommends \
    cuda-core-$CUDA_PKG_VERSION \
    cuda-misc-headers-$CUDA_PKG_VERSION \
    cuda-command-line-tools-$CUDA_PKG_VERSION \
    cuda-nvrtc-dev-$CUDA_PKG_VERSION \
    cuda-nvml-dev-$CUDA_PKG_VERSION \
    cuda-nvgraph-dev-$CUDA_PKG_VERSION \
    cuda-cusolver-dev-$CUDA_PKG_VERSION \
    cuda-cublas-dev-8-0=8.0.61.2-1 \
    cuda-cufft-dev-$CUDA_PKG_VERSION \
    cuda-curand-dev-$CUDA_PKG_VERSION \
    cuda-cusparse-dev-$CUDA_PKG_VERSION \
    cuda-npp-dev-$CUDA_PKG_VERSION \
    cuda-cudart-dev-$CUDA_PKG_VERSION \
    cuda-driver-dev-$CUDA_PKG_VERSION && \
    \
    apt clean && \
    rm -rf /var/lib/apt/lists/*

#Hoomd
ENV SOFTWARE_ROOT=/opt/hoomd
ENV PYTHONPATH=$PYTHONPATH:${SOFTWARE_ROOT}/lib/python
RUN apt-get update && apt-get install -y --no-install-recommends git cmake && \
    apt clean && \
    rm -rf /var/lib/apt/lists/* && \
    pip --no-cache-dir install numpy && \
    export HOOMD_TAG=v2.3.5 && \
    git clone --recursive https://bitbucket.org/glotzer/hoomd-blue && \
    cd hoomd-blue && \
    git checkout  $HOOMD_TAG && \
    mkdir build && \
    cd build && \
    export CXX="$(command -v g++)" && \
    export CC="$(command -v gcc)" && \
    cmake ../ -DCMAKE_INSTALL_PREFIX=${SOFTWARE_ROOT}/lib/python \
              -DENABLE_CUDA=ON \
              -DDISABLE_SQLITE=ON \
              -DSINGLE_PRECISION=ON && \
    make -j3 && \
    make install

# mount points for filesystems on clusters
RUN mkdir -p /nfs \
    mkdir -p /oasis \
    mkdir -p /scratch \
    mkdir -p /work \
    mkdir -p /projects \
    mkdir -p /home1
