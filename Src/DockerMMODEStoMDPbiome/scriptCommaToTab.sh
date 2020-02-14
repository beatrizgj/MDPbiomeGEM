#!/bin/bash
for filename in *.tsv; do
    sed -i "s/,/\t/g" $filename 
done
