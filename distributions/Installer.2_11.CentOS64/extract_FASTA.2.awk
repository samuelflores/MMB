BEGIN {
ORS = ""
CHAIN = "CHAIN-NOT-SET"
FIRSTRESIDUE = -1111
BIOPOLYMER = "NOT-SET"
BIOPOLYMER = "protein";
}

(FNR == 1) {print "\n\n";}

( ($1 =="ATOM")  ) {
#((substr($1,1,6) =="HETATM") || ($1 =="ATOM")  ) {
if (substr($0,22,1) != CHAIN  ) {
#print FIRSTRESIDUE ,"!= -1111\n";
if (FIRSTRESIDUE != -1111) {
#print "\nmobilizer Rigid ",CHAIN,"\n ";
#print "constrainToGround ",CHAIN," ",FIRSTRESIDUE,"\n";
}
#print "\n#Chain ID: ",substr($0,22,1) ,"\n";
CHAIN = substr($0,22,1)
FIRSTRESIDUE = substr($0,23,4);
LASTRESIDUE = FIRSTRESIDUE - 1;
#print "# First residue # ",FIRSTRESIDUE,"\n";
print "\n",BIOPOLYMER," ",CHAIN," ",FIRSTRESIDUE," ";
if ( ($3 == "C")) {
#BIOPOLYMER = "protein";
#print "Biopolymer type = ",BIOPOLYMER,"\n";
}
}
}

(((substr($1,1,6) =="HETATM") || ($1 =="ATOM")) && ($3 == "C")) {
BIOPOLYMER = "protein";
LASTRESIDUE =LASTRESIDUE +1;
if (substr($0,18,3) =="ALA") {print "A"}
if (substr($0,18,3) =="ARG") {print "R"}
if (substr($0,18,3) =="ASN") {print "N"}
if (substr($0,18,3) =="ASP") {print "D"}
if (substr($0,18,3) =="CYS") {print "C"}
if (substr($0,18,3) =="GLU") {print "E"}
if (substr($0,18,3) =="GLN") {print "Q"}
if (substr($0,18,3) =="GLY") {print "G"}
if (substr($0,18,3) =="HIS") {print "H"}
if (substr($0,18,3) =="ILE") {print "I"}
if (substr($0,18,3) =="LEU") {print "L"}
if (substr($0,18,3) =="LYS") {print "K"}
if (substr($0,18,3) =="MET") {print "M"}
if (substr($0,18,3) =="PHE") {print "F"}
if (substr($0,18,3) =="PRO") {print "P"}
if (substr($0,18,3) =="SER") {print "S"}
if (substr($0,18,3) =="THR") {print "T"}
if (substr($0,18,3) =="TRP") {print "W"}
if (substr($0,18,3) =="TYR") {print "Y"}
if (substr($0,18,3) =="VAL") {print "V"}

#RNA section 
#if (substr($0,18,3) =="A  ") {print "A"}
#if (substr($0,18,3) =="U  ") {print "U"}
#if (substr($0,18,3) =="G  ") {print "G"}
#if (substr($0,18,3) =="C  ") {print "C"}
#if (substr($0,18,3) =="  A") {print "A"}
#if (substr($0,18,3) =="  U") {print "U"}
#if (substr($0,18,3) =="  G") {print "G"}
#if (substr($0,18,3) =="  C") {print "C"}

}
(((substr($1,1,6) =="HETATM") ||($1 =="ATOM")) && (($3 == "C3'") || ($3 == "C3*"))) {
    BIOPOLYMER = "Nucleic Acid";
    #print $0"\n";
    TEMP  =substr($0,18,3);
    gsub(/ */,"", TEMP);
    print TEMP ;
    
}
END   {

if (FIRSTRESIDUE != -1111) {
#print "\nmobilizer ",CHAIN," ",FIRSTRESIDUE," Rigid ",CHAIN," ",LASTRESIDUE,"\n";
#print "constrainToGround ",CHAIN," ",FIRSTRESIDUE,"\n";
}
ORS = "\n"; print "\n\n"

}
