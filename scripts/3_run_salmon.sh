#!/bin/bash
while IFS='' read -r line || [[ -n "$line" ]]; do
    salmon quant -i /home/poxley/Documents/bioinformatics/indices/mmusculus/mmus_index \
		-l A \
		$line \
		-p 20 
done < "$1"
