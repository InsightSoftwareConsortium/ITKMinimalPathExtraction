FROM debian:8
MAINTAINER Insight Software Consortium <community@itk.org>

RUN apt-get update && apt-get install -y \
  build-essential \
  curl \
  cmake \
  git \
  libexpat1-dev \
  libhdf5-dev \
  libjpeg-dev \
  libpng12-dev \
  libpython3-dev \
  libtiff5-dev \
  python \
  ninja-build \
  wget \
  vim \
  zlib1g-dev

RUN mkdir -p /usr/src/ITKMinimalPathExtraction-build
WORKDIR /usr/src

# master 2016-07-06
ENV ITK_GIT_TAG 71a779150ee5797e6ca4bef01906d5c12a18188b
RUN git clone git://itk.org/ITK.git && \
  cd ITK && \
  git checkout ${ITK_GIT_TAG} && \
  cd ../ && \
  mkdir ITK-build && \
  cd ITK-build && \
  cmake \
    -G Ninja \
    -DCMAKE_INSTALL_PREFIX:PATH=/usr \
    -DBUILD_EXAMPLES:BOOL=OFF \
    -DBUILD_TESTING:BOOL=OFF \
    -DBUILD_SHARED_LIBS:BOOL=ON \
    -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON \
    -DITK_LEGACY_REMOVE:BOOL=ON \
    -DITK_BUILD_DEFAULT_MODULES:BOOL=OFF \
    -DITK_USE_SYSTEM_LIBRARIES:BOOL=ON \
    -DModule_ITKOptimizers:BOOL=ON \
    -DModule_ITKPath:BOOL=ON \
    -DModule_ITKFastMarching:BOOL=ON \
    -DModule_ITKTestKernel:BOOL=ON \
    -DModule_ITKIOSpatialObjects:BOOL=ON \
    ../ITK && \
  ninja && \
  find . -name '*.o' -delete

