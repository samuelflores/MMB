open RESTRAINTS, $ARGV[0]    or die $!;
@restraints = <RESTRAINTS> ;
#print '@restraints', @restraints,"\n";
#print 'restraints[0]', $restraints[0],"\n";
#print 'restraints[1]', $restraints[1],"\n";
close (RESTRAINTS);
int r;
#for (r = 0; r < 1; r++) 
for ($r = 0; $r < scalar(@restraints); $r++) 
{
    #print "restraints[$r] = $restraints[$r] \n";
    int $i ; $i =  substr($restraints[$r],0,3);
    int $j ; $j =  substr($restraints[$r],10,3); 
    int $dist ; $dist =  substr($restraints[$r],20,5); 
    #print "i = ",$i, " j = ", $j, "\n";
    $atomName1 = "*";    
    $atomName2 = "*";    
    int $residueNumber1; $residueNumber1 = -1;
    int $residueNumber2; $residueNumber2 = -1;
    $chain1 = "A";        
    $chain2 = "A";         
    open PDB, $ARGV[1]    or die $!;
    while (<PDB>) {
        int $PDBatomNumber; $PDBatomNumber = substr($_,6,5);
        #print "PDBatomNumber = $PDBatomNumber\n";
        if ($PDBatomNumber == $i)        {
            $atomName1      = substr($_,12,4);   
            $residueNumber1 = substr($_,22,4);   
            $chain1         = substr($_,21,1);   
            #print '$chain1 ',$chain2,"\n";
        }
        if ($PDBatomNumber == $j)        {
            $atomName2      = substr($_,12,4);   
            $residueNumber2 = substr($_,22,4);   
            $chain2         = substr($_,21,1);   
        }
        
    }
    print "atomTether $chain1 $residueNumber1 $atomName1  $chain2 $residueNumber2 $atomName2 $dist \n";

}

