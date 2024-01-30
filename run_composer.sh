#!/bin/bash

for file in /tmp/mj/sbs_files/*.sbs
do
    filename=$(basename "$file" .sbs)
    echo $filename
    docker run --rm -e DB_HOST="172.17.0.1" -v /tmp/:/tmp/ --entrypoint ./rnacomposer ext-rnacomposer-engine -i/tmp/mj/RNAtools/$filename.dot -usp/tmp/mj/sbs_files/$filename.sbs
done