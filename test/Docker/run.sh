#!/bin/sh

script_dir="`cd $(dirname $0); pwd`"

docker run \
  --rm \
  -v $script_dir/../..:/usr/src/ITKMinimalPathExtraction \
    insighttoolkit/minimalpathextraction-test \
      /usr/src/ITKMinimalPathExtraction/test/Docker/test.sh
