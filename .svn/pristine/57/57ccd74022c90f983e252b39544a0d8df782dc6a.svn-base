#!/usr/bin/perl

# This program renumbers an input pdb file with sequential integers, starting at a user-specified integer.  Gaps in numbering disappear, and inserted regions are renumbered as integers.
# Usage:  ./renumber.pl [input file name] [new first residue number] > [output file name]
print "REMARK Usage: ./renumber.pl  [input file name] [new first residue number] > [output file name]\n";
#open FILE, "input.pdb" or die $!;
open FILE, $ARGV[0]    or die $!;
$oldResNum     = "ZZZZ";
$currentResNum = "XXXX";
$renumberedResNum = $ARGV[1]-1;
#print "REMARK renumberedResNum =$renumberedResNum \n";
while (<FILE>) {
    if ((substr($_,0,4) eq "ATOM") || (substr($_,0,6) eq "HETATM")) {
        #print substr($_,0,4),"\n";
        #print $_;     
        $currentResNum = substr($_,22,5);
        #print "currentResNum =  >",$currentResNum,"<\n";
        #print "oldResNum =  >",$oldResNum,"<\n";
        if ($currentResNum ne $oldResNum) { 
            $renumberedResNum++;
        }
        #print "renumberedResNum =  >",$renumberedResNum,"<\n";
        $paddedRenumberedResNum= sprintf("%04d", $renumberedResNum);
        #print "paddedRenumberedResNum =  >",$paddedRenumberedResNum,"<\n";

        substr($_,22,4) = $paddedRenumberedResNum;
        substr($_,26,1) = " ";                    
        print $_;
        $oldResNum = $currentResNum;
    }


}

