#!/bin/bash
set -euo pipefail

cd $1



ls | sort

find . -name "*.metrics" -not -name "neat_5x_EX_hg19_picard.downsample.metrics" -xtype f -exec sh -c "cat {} | md5sum" \;
