#!/bin/bash

#extracts columns 6 to 7 from the file 'Tco2.fqfiles' and writes the output to a file named 'srr.list'
cut -f6-7 Tco2.fqfiles >srr.list
#resd each line from 'srr.list'assigning the values in columns 6 and 7 to the variables 'old' and 'new'
#and add the fastqc command into 'Para.sh'file
while read old new
do echo "fastqc -t 4 -o ../fastqc $old $new" >Para.sh
done<srr.list



