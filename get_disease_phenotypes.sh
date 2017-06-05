#!/bin/bash

awk '{if ($0 ~ /CLNSIG=5/) print $0}' GRCh38_latest_clinvar.vcf| egrep -o 'CLNDBN[^:]+' | awk -F";" '{print $1}' | sort | uniq -c | sort -nr | grep -v "not_provided" | tr '|' '\n' | tr '\' '\n' | sed 's/CLNDBN[=]//g' | sort -nr  | head -100 > top100_disease_phenotypes
