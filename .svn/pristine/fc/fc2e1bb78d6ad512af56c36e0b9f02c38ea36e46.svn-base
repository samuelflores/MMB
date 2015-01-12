 ((substr($1,1,6) =="HETATM") || ($1 =="ATOM"))  {
print $0;
print "REMARK-SIMTK-COORDS ",substr($0,31,8)"00000000000000"," "substr($0,39,8)"00000000000000"," " substr($0,47,8)"00000000000000";
}

