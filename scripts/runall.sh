#!/bin/bash

batch=1
for i in `ls simulated_showers/merged/*.root` ; do  
    ./python/showerAnalyze.py $i $batch
done
