BEGIN {
ORS = ""
CHAIN = "CHAIN-NOT-SET"
FIRSTRESIDUE = -1111
BIOPOLYMER = "NOT-SET"
BIOPOLYMER = "protein";
}

(FNR == 1) {print "\n\n";}

((substr($1,1,6) =="HETATM") || ($1 =="ATOM")  ) {

sub ( /C1\*/ , "C1\'" , $0);
sub ( /C2\*/ , "C2\'" , $0);
sub ( /C3\*/ , "C3\'" , $0);
sub ( /C4\*/ , "C4\'" , $0);
sub ( /C5\*/ , "C5\'" , $0);

sub ( /O2\*/ , "O2\'" , $0);
sub ( /O3\*/ , "O3\'" , $0);
sub ( /O4\*/ , "O4\'" , $0);
sub ( /O5\*/ , "O5\'" , $0);



CHAIN = substr($0,22,1)
FIRSTRESIDUE = substr($0,23,4);
LASTRESIDUE = FIRSTRESIDUE - 1;
#print BIOPOLYMER," ",CHAIN," ",FIRSTRESIDUE," ";


if (substr($0,22,1) != CHAIN  ) {
#if (FIRSTRESIDUE != -1111) {
#    if ($1 != "TER"){
#        print "TER \n";
#    }
#}

if ( ($3 == "C")) {
}
}
# get atom name:
FIRSTATOMLETTER = substr($0,12,4 );
RESNAME         = substr($0,18,3 );
#print FIRSTATOMLETTER, " >1\n";
# get rid of whitespace and numbers:
gsub ( / / , "" , FIRSTATOMLETTER);
ATOMNAMENOWHITESPACE = FIRSTATOMLETTER;
# molmodel uses CD rather than the PDB-standard CD1
#print ATOMNAMENOWHITESPACE,", ",RESNAME,"\n";
if ((ATOMNAMENOWHITESPACE == "CD") && (RESNAME == "ILE")) {
    sub ( /CD./ , "CD1" , $0 );
}
gsub ( /[0-9]/ , "" , FIRSTATOMLETTER);
#print FIRSTATOMLETTER, " >2 \n";
#sub ( / / , "" , FIRSTATOMLETTER);
#sub ( / / , "" , FIRSTATOMLETTER);
#print FIRSTATOMLETTER, "\n";
# first remaining letter should be element type
FIRSTATOMLETTER = substr(FIRSTATOMLETTER ,0,1  );
#print FIRSTATOMLETTER, "\n";
if (FIRSTATOMLETTER != "H")
    print $0"\n";
}

(((substr($1,1,6) =="HETATM") || ($1 =="ATOM")) && ($3 == "C")) {
BIOPOLYMER = "protein";
LASTRESIDUE =LASTRESIDUE +1;
    #print $0,"\n";

}
(((substr($1,1,6) =="HETATM") ||($1 =="ATOM")) && (($3 == "C3'") || ($3 == "C3*"))) { # else this is RNA, strip out whitespaces and use single letter code:
    BIOPOLYMER = "Nucleic Acid";
    #print $0"\n";
    TEMP  =substr($0,18,3);
    gsub(/ */,"", TEMP);
    #print TEMP ;
    #print $0,"\n";
}

((substr($1,1,6) !="HETATM") && ($1 !="ATOM") && (substr($1,1,6) !="REMARK")  ) {
print $0,"\n";
}

END   {

#if (FIRSTRESIDUE != -1111) {
#    if ($1 != "TER") {
#        print "TER  \n";
#    }
#}
ORS = "\n"; print "\n\n"
print "END\n";
}
