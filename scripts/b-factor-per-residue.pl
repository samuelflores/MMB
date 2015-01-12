#!/usr/bin/perl

# This program takes input pdb file and averages the b-factor for each residue, then outputs this value. 
# Usage:  ./renumber.pl [input file name]  > [output file name]

#open FILE, "input.pdb" or die $!;
open FILE, $ARGV[0]    or die $!;
$averageBFactor = 0.0;
$numAtomsThisResidue = 0;
$oldResNum     = "ZZZZ";
$currentResNum = "XXXX";
$renumberedResNum = $ARGV[1]-1;
while (<FILE>) {
    if (substr($_,0,4) == "ATOM") {
        $currentResNum = substr($_,22,5);
        #print "currentResNum =  >",$currentResNum,"<\n";
        #print "oldResNum =  >",$oldResNum,"<\n";
        if ($currentResNum ne $oldResNum) { 
            if ($numAtomsThisResidue>0) {
                #print "Sum of B-factors = $averageBFactor.  Number of atoms this residue = $numAtomsThisResidue.  Average B-factor this residue = ";
                print ($averageBFactor/$numAtomsThisResidue),"\n";
                print "\n";
                }
            $averageBFactor = 0.0;
            $numAtomsThisResidue = 0;
            $renumberedResNum++;
        }
        #print "current B-factor is : ",substr($_,60,6),"\n";
        $averageBFactor += substr($_,60,6);
        $numAtomsThisResidue ++;
        $paddedRenumberedResNum= sprintf("%04d", $renumberedResNum);
        #print "paddedRenumberedResNum =  >",$paddedRenumberedResNum,"<\n";

        substr($_,22,4) = $paddedRenumberedResNum;
        substr($_,26,1) = " ";                    
        #print $_;
        $oldResNum = $currentResNum;
    }


}
if ($numAtomsThisResidue>0) {
    #print "Sum of B-factors = $averageBFactor.  Number of atoms this residue = $numAtomsThisResidue.  Average B-factor this residue = ";
    print ($averageBFactor/$numAtomsThisResidue),"\n";
    print "\n";
    }

