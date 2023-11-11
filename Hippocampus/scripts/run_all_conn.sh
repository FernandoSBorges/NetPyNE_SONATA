#!/bin/bash

for i in {0..18}; do 
   echo $i | time -f "Elapsed Time = %E, CPU Usage = %P, CPU Time= %S, Job = $i" -o time.txt --append ipython connfullca1.py $i &
   sleep 1
done

