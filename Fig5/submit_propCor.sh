#!/bin/bash

njobs=2000

qsub -q w-bigmem.q -e iotrash/ -o iotrash/ -l h='biostat-b17|biostat-b18|biostat-b16' -l h_vmem=10G -t 1-$njobs -tc 50 ./call_propCor.sh $1
