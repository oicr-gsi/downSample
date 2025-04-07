#!/bin/bash
set -euo pipefail

cd $1



ls | sort

find . -name *.metrics -xtype f -exec sh -c "cat {} | md5sum" \;
