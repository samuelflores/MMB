#!/usr/bin/perl

print "REMARK Usage: ./renumber-loop.pl  [input file name] \n";
print "REMARK Expecting an input file with a single residue, with chain ID A, and another single residue with chain ID B\n";
#$oldResNum     = "ZZZZ";
#$currentResNum = "XXXX";
#$renumberedResNum = $ARGV[1]-1;
$maxResidueNumber = 9999;
my @nArray = (1..2000);
for $n (@nArray) {
    #print "REMARK working on $n \n";
    open FILE, $ARGV[0]    or die $!;
    open OUTFILE, ">base-pair.$n.pdb"    or die $!;
    while (<FILE>) {
        if ((substr($_,0,4) eq "ATOM") || (substr($_,0,6) eq "HETATM")) {
            $currentChain = substr($_,21,1);
            if ($currentChain eq "A"){
                print ">A< >$currentChain<   > $_ ";
                printf OUTFILE  "%s%04d%s",substr($_,0,22),                   $n, substr($_,26  );   
            } elsif ($currentChain eq "B"){
                print "B $currentChain> $_ ";
                printf OUTFILE  "%s%04d%s",substr($_,0,22), $maxResidueNumber-$n, substr($_,26  );   
            } else {
                print "Was expecting chain ID either A or B\n";
                exit(1);
            }
        }
    }
    close OUTFILE;
    close FILE;
}

