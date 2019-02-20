#!/bin/bash
while IFS='' read -r line || [[ -n "$line" ]]; do
    java -jar /home/poxley/bin/trimmomatic.jar PE -phred33 $line ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done < "$1"
