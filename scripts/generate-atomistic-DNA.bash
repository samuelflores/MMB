#!/bin/bash

# Short procedure to pad residue numbers with zeroes. Otherwise we will find it difficult to sort properly.
pad_residue_numbers () {
counter=1
until [ $counter -gt 9999 ]
do
    if test -f "last.$counter.pdb"; then
        echo "padding last.$counter.pdb residue numbers with zeroes"
        cp last.$counter.pdb temp.1.pdb
        cat temp.1.pdb | awk '{if (substr($0,23,1) == " ") {print substr($0,1,22)"0"substr($0,24)}else {print}}' | awk '{if (substr($0,24,1) == " ") {print substr($0,1,23)"0"substr($0,25)}else {print}}' | awk '{if (substr($0,25,1) == " ") {print substr($0,1,24)"0"substr($0,26)}else {print}}' > temp.2.pdb
	mv -f temp.2.pdb last.$counter.pdb
    fi
    ((counter++))
done
}

rm -f  last.*.pdb; 	docker run  -v $(pwd):/work --rm -it samuelflores/mmb-ubuntu:4.0.0 MMB -C 	commands.spiral.dat ; pad_residue_numbers ; cat last.?.pdb last.??.pdb last.???.pdb last.????.pdb > temp.1.pdb; sed -i 's/H5T/P\ \ /g' ./temp.1.pdb; awk '{if ($3 == "P"){print $5"P"$2,$0} else {print $5"Z"$2,$0}}' temp.1.pdb > temp.2.pdb; sort -k1 temp.2.pdb > temp.3.pdb; awk '{print substr($0,13)}' temp.3.pdb | grep -v REMARK | grep -v H3T > DNA-spiral.pdb; 

echo "all done."
exit 0
