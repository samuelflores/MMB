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
