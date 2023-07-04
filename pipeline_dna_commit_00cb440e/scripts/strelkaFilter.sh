#!/bin/bash

grep '^#' $1 > $2
awk '($7 == "PASS")'  $1 >> $2

