#!/bin/bash
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
awk 'BEGIN {x=0}{if ($0 == "### Pipeline parameters") x=1; if (x==1) print "# " $0}' config.r > config-custom.r
