#!/bin/bash

# This script is for use with the sphericalHelix command. Takes the last.*.pdb files and concatenates them, first one chain, then the other, in order of residue number.

rm DA.pdb;
grep -h DA last.?.pdb > DA.pdb
grep -h DA last.??.pdb >> DA.pdb
grep -h DA last.???.pdb >> DA.pdb
grep -h DA last.????.pdb >> DA.pdb
rm DT.pdb
grep -h DT last.?.pdb > DT.pdb
grep -h DT last.??.pdb >> DT.pdb
grep -h DT last.???.pdb >> DT.pdb
grep -h DT last.????.pdb >> DT.pdb

cat DA.pdb > DNA-spiral.pdb
tac DT.pdb >> DNA-spiral.pdb

sed -i 's/H5T/P  /g' DNA-spiral.pdb 
awk '{if ($3 == "P"){print $5"P"$2,$0} else {print $5"Z"$2,$0}}' DNA-spiral.pdb > temp.2.pdb
sort -k1 temp.2.pdb > temp.3.pdb
cut -b 13-90 temp.3.pdb > DNA-spiral.pdb
rm temp.2.pdb temp.3.pdb DA.pdb DT.pdb
