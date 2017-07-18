#!/usr/bin/bash

## Generate SEQs

###for i in `cat snp_lists`; do echo $i; bash anmol_robot2.sh "$i"; done

####### Come up with a clever line here which checks for when there has been no size change to the SNP_SEQs file for 60 seconds

# makeblastdb

makeblastdb -dbtype nucl -in SNP_SEQs.foo

rm SNP_SEQs.foo

## Do BLAST searches of whatever SRA you are interested in

for i in `cat SRR_list`; do magicBLAST -db SNP_SEQs.foo -sra $i; done  

## Parse the output


