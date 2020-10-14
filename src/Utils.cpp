/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#include "Utils.h"

//#include <boost/lexical_cast.hpp>


#include  <sstream> 
using namespace std;
using namespace SimTK;


int myMkdir(std::string directoryPath){
    MMBLOG_FILE_FUNC_LINE(INFO, " You are asking to create the directory  "<<directoryPath<<std::endl);
    if (!(opendir(directoryPath.c_str()))){
        MMBLOG_FILE_FUNC_LINE(INFO, " opendir failed to open directory "<<directoryPath<<" . Will now create this directory.  " <<std::endl);
        const int dir_err = mkdir(directoryPath.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
        MMBLOG_FILE_FUNC_LINE(INFO, ": mkdir returned : "<<dir_err <<std::endl);
        if (-1 == dir_err)
        {
            printf("Error creating directory!n");
            MMBLOG_FILE_FUNC_LINE(INFO, " Failed to  create directory "<<directoryPath<<"  " <<std::endl);
            exit(1);
        } else if (0 == dir_err) {
            MMBLOG_FILE_FUNC_LINE(INFO, " Successfully created directory "<<directoryPath<<" " <<std::endl);
        } else {
            MMBLOG_FILE_FUNC_LINE(INFO, " An unexpected error occurred when creating the directory "<<directoryPath<<" " <<std::endl);
        }
    }
    if (access((directoryPath ).c_str(), R_OK) == 0) {
        MMBLOG_FILE_FUNC_LINE(INFO, " Found that we have read access to a directory called "<< directoryPath   <<" . So far so good."<<std::endl);
        return 0;
    } else {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, " Heads up! Found that we do NOT have read access to a directory called "<<  directoryPath  <<"  "<<std::endl);
        return 1;
    }
}

int myChdir(std::string directoryPath){
    MMBLOG_FILE_FUNC_LINE(INFO, " About to attempt changing directory to "<<directoryPath<<" . "<<std::endl);
    if (chdir(directoryPath.c_str())){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, " Unable to change directory to "<<directoryPath<<" . Exiting now."<<std::endl);
        return 1;
    } else {
        MMBLOG_FILE_FUNC_LINE(INFO, " Was able to successfully change directory to "<<directoryPath<<" . "<<std::endl);
        return 0;
    }
}

void closingMessage() {
    std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
    std::cout<<__FILE__<<":"<<__LINE__<<" Questions or problems? Please use the public forum at https://simtk.org/forums/viewforum.php?f=359&amp;sid=770b2f0d333ecb740d8c2f9e7e80e51c  .. or email me."<<std::endl;
    std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
    std::cout<<__FILE__<<":"<<__LINE__<<" Please support continued development. "<<std::endl;
    //std::cout<<__FILE__<<":"<<__LINE__<<" Suggested donation: "<<std::endl;
    //std::cout<<__FILE__<<":"<<__LINE__<<" Industry (per user): 1000 EUR "<<std::endl;
    //std::cout<<__FILE__<<":"<<__LINE__<<" By bank transfer to IBAN: SE0750000000053680279418 , SWIFT: ESSESESS "<<std::endl;
    std::cout<<__FILE__<<":"<<__LINE__<<" Industrial and other inquiries: samuel.flores@scilifelab.se "<<std::endl;
    std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
};

CheckFile::CheckFile(const String myFileName){
    fileName = myFileName;
    //struct stat st;
    stat(  fileName.c_str(), &st);
}

void CheckFile::validateNonZeroSize(){
    MMBLOG_FILE_FUNC_LINE(INFO, " About to check that "<<fileName<<" has nonzero size.."<<std::endl);
    if ( st.st_size == 0){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Apparently "<<fileName<<" has size "<<st.st_size <<" . Dying now."<<std::endl);
    } else {
        MMBLOG_FILE_FUNC_LINE(INFO, " Apparently "<<fileName<<" has size "<<st.st_size <<" . This seems OK."<<std::endl);
    }
}

void CheckFile::validateExists(){
    if(stat(fileName.c_str(), &st) != 0){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, " ERROR: stat says no file exists. Dying now."<<std::endl);
    } else {
        MMBLOG_FILE_FUNC_LINE(INFO, " stat says file exists. All is good. "<<std::endl);
    }
}



std::string   trim(const std::string& str, 
                 const std::string& whitespace// = " \t"
                 )   
{ 
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1; 

    return str.substr(strBegin, strRange);
    //return std::string ("hello");
}

bool vectorCompare(String myString, vector<String> & comparisonStringVector) {
    if (comparisonStringVector.size() == 0) { //cout<<", returning TRUE "<<endl ;
        return true;} // If we are comparing to an empty vector, return true.  This is in case no partner chains have been specified, in which case any chain will pass.
    for (int i = 0; i < comparisonStringVector.size(); i++) {
        //cout<<__FILE__<<":"<<__LINE__<<" comparing "<<myString<< " to comparisonStringVector["<<i<<"] : "<<comparisonStringVector[i];
        if (comparisonStringVector[i].compare(myString ) == 0) { //cout<<", returning TRUE "<<endl ;  
            return true; } else {//cout<<", returning FALSE";
            }
        //cout<<endl;
    }
    return false; // If no String in comparisonStringVector is the same as myString
};

BondMobility::Mobility stringToBondMobility(String bondMobilityString) {
       String myBondMobilityString =   bondMobilityString;
       myBondMobilityString.toUpper();
       BondMobility::Mobility myBondMobility;
       // Remember  Free = 1, Torsion = 2, Rigid = 3
       if ((myBondMobilityString).compare("RIGID") == 0)              {
           MMBLOG_FILE_FUNC_LINE(INFO, " Detected Rigid : >"<<bondMobilityString<<"< or >"<<myBondMobilityString<<"< "<<std::endl);
           myBondMobility = SimTK::BondMobility::Rigid;
           //std::cout<<__FILE__<<":"<<__LINE__<<" returning myBondMobility = >"<<myBondMobility<<"< "<<std::endl;
       }
       else if ((myBondMobilityString).compare("TORSION") == 0)       {
           MMBLOG_FILE_FUNC_LINE(INFO, " Detected Torsion :"<<myBondMobilityString<<std::endl);
           myBondMobility = SimTK::BondMobility::Torsion;}
       else if ((myBondMobilityString).compare("DEFAULT") == 0)       {
           MMBLOG_FILE_FUNC_LINE(INFO, " Detected Default :"<<myBondMobilityString<<std::endl);
           myBondMobility = SimTK::BondMobility::Default;}
       else if ((myBondMobilityString).compare("FREE")  == 0)         {
           MMBLOG_FILE_FUNC_LINE(INFO, " Detected Free :"<<myBondMobilityString<<std::endl);
           myBondMobility = SimTK::BondMobility::Free ;}
       else {
           MMBLOG_FILE_FUNC_LINE(CRITICAL, " At this time only Default, Free,Torsion, and Rigid bondMobilities are supported. You are attempting to apply \"" << myBondMobilityString <<"\". "<<std::endl);
       }
       return myBondMobility;
}


void InterfaceContainer::addInterface(vector<String> myChains,vector<String> partnerChains,  double myDepth ,  String myMobilizerString ){
    Interface myInterface; 
    myInterface.Chains.clear(); myInterface.PartnerChains.clear();
    for (int i = 0; i < myChains.size(); i++) {myInterface.Chains.push_back( myChains[i]);}
    for (int i = 0; i < partnerChains.size(); i++) {myInterface.PartnerChains.push_back( partnerChains[i]);}
    myInterface.Depth = myDepth; 
    if (! (myMobilizerString.compare("NONE")==0) ||
	  (myMobilizerString.compare("Rigid")==0) ||
	  (myMobilizerString.compare("Default")==0) || 
	  (myMobilizerString.compare("Torsion")==0) ||
	  (myMobilizerString.compare("Free")==0) 
	) {
       MMBLOG_FILE_FUNC_LINE(CRITICAL, " Expected a mobilizer type (Default, Rigid, Torsion, Free), but got : >"<< myMobilizerString<<"< "<<std::endl);
    }
    myInterface.MobilizerString = myMobilizerString; 
    interfaceVector.push_back(myInterface); 
};


#ifdef USE_OPENMM
vector<TwoAtomClass> InterfaceContainer::retrieveCloseContactPairs(vector<MMBAtomInfo> & concatenatedAtomInfoVector ){
//vector<TwoAtomClass> InterfaceContainer::retrieveCloseContactPairs(  BiopolymerClassContainer & myBiopolymerClassContainer){
        OpenMM::NeighborList neighborList;
        openmmVecType boxSize = openmmVecType(10000,10000,10000);
        //vector<MMBAtomInfo> concatenatedAtomInfoVector = myBiopolymerClassContainer.getConcatenatedAtomInfoVector();
        vector<TwoAtomClass> contactingAtomInfoPairVector;
        contactingAtomInfoPairVector.clear();
        vector<openmmVecType> particleList(concatenatedAtomInfoVector.size());
        vector<set<int> > exclusions( particleList.size() );
        for (int i = 0; i < concatenatedAtomInfoVector.size() ; i++) {
            particleList[i] = concatenatedAtomInfoVector[i].position;
        }

        MMBLOG_FILE_FUNC_LINE(INFO, " neighborList size is : "<<neighborList.size()<<endl);
        for (int h = 0 ; h < numInterfaces(); h++ ){ // loop through interfaceContainer interfaces ..
            vector<String> referenceChains = getInterface(h).getChains();  
            vector<String> partnerChains = getInterface(h).getPartnerChains();  
            double         radius        = getInterface(h).getDepth();  
            MMBLOG_FILE_FUNC_LINE(INFO, "Now turning interface "<< h << " to individual constraints between pairs of atoms."<<endl);
            getInterface(h).print(); 
            cout<<__FILE__<<":"<<__LINE__<<endl;
            computeNeighborListVoxelHash(neighborList, particleList.size() , particleList, exclusions, &boxSize, false, radius  , 0.0);
            for ( int j = 0 ; j < neighborList.size(); j++) {
                if ((((( vectorCompare(concatenatedAtomInfoVector[neighborList[j].first].chain , (referenceChains))) == 1) &&
                     (( vectorCompare(concatenatedAtomInfoVector[neighborList[j].second].chain ,(  partnerChains))) == 1))  || 
//Use an XOR here. This means if the 'partnerChains' evaluation is later set to return 1 when partnerChains is empty, this will still work.
                    ((( vectorCompare(concatenatedAtomInfoVector[neighborList[j].second].chain ,(referenceChains))) == 1) &&
                     (( vectorCompare(concatenatedAtomInfoVector[neighborList[j].first].chain  ,(  partnerChains))) == 1))     //Make sure that exactly one residue is in the 'referenceChains', and the other residue is in the 'partnerChains' .. thus only the desired interface is included
                                                                                                                         ) 
                     && (concatenatedAtomInfoVector[neighborList[j].first].chain.compare(concatenatedAtomInfoVector[neighborList[j].second].chain) != 0 )
                   ) // lastly,make sure the two atoms are not in the same chain.

                {   
                    TwoAtomClass myTwoAtomClass(
			(concatenatedAtomInfoVector[neighborList[j].first].chain),
			(concatenatedAtomInfoVector[neighborList[j].first].residueID),
			(concatenatedAtomInfoVector[neighborList[j].first].atomName),
			(concatenatedAtomInfoVector[neighborList[j].second].chain),
			(concatenatedAtomInfoVector[neighborList[j].second].residueID),
			(concatenatedAtomInfoVector[neighborList[j].second].atomName)//,
                        //(concatenatedAtomInfoVector[neighborList[j].first ].position - concatenatedAtomInfoVector[neighborList[j].second].position) // later, compute distance
                    );
                    contactingAtomInfoPairVector.push_back(myTwoAtomClass );
                    MMBLOG_FILE_FUNC_LINE(INFO, " Detected contact : ");
                    myTwoAtomClass.print();
                }
                else {
                    // Do nothing.
                }
            }
        }
        return contactingAtomInfoPairVector;  
};
#endif

ConstraintClass::ConstraintClass(){
        chain1 = ""; residueID1 = ResidueID(); atomName1 = ""; 
        chain2 = ""; residueID2 = ResidueID(); atomName2 = ""; 
        constraintType = WeldToGround ; 
        };  
ConstraintClass::ConstraintClass(String myChain, ResidueID inputResidueID,String myAtomName) {
        residueID1 = (inputResidueID);
        atomName1 = myAtomName;
        chain1 = myChain;
        residueID2 = ResidueID();
        atomName2 = "" ;    
        chain2 = ""; 
        constraintType = ( WeldToGround );
        //toGround = true ;
    }; 

ConstraintClass::ConstraintClass(String myChain, ResidueID inputResidueID,String myAtomName,String myChain2, ResidueID inputResidueID2,String myAtomName2, ConstraintType myConstraintType) {
        residueID1 = (inputResidueID);
        //residueID1.setInsertionCode ( residueID1.getInsertionCode());
        atomName1 = myAtomName;
        chain1 = myChain;
        residueID2 = (inputResidueID2);
        atomName2 = myAtomName2;
        chain2 = myChain2;
        setConstraintType(myConstraintType);
        //toGround = false;
    };
/*
Array_<MobilizedBodyIndex> ConstraintClass::fetchMobilizedBodyIndexArray_(BiopolymerClassContainer myBiopolymerClassContainer,SimbodyMatterSubsystem & matter ) {
        Array_< MobilizedBodyIndex >    coordMobod(2);
        coordMobod[0] =  myBiopolymerClassContainer.updBiopolymerClass(chain1).getAtomMobilizedBodyIndex(matter,residueID1    ,atomName1    );
        coordMobod[1] =  myBiopolymerClassContainer.updBiopolymerClass(chain2).getAtomMobilizedBodyIndex(matter,residueID2    ,atomName2    );
        return coordMobod;
    };
Array_<MobilizerQIndex> ConstraintClass::fetchMobilizerQIndexArray_(BiopolymerClassContainer myBiopolymerClassContainer,SimbodyMatterSubsystem & matter, State & state) {
        Array_< MobilizerQIndex >       coordQIndex(2);
        coordQIndex[0] =  MobilizerQIndex(0); //mobilizedBody1.getMobilizerQIndex(state);
        coordQIndex[1] =  MobilizerQIndex(0); // // mobilizedBody2.getFirstQIndex(state);
        return coordQIndex;
    }
*/

void ConstraintClass::setConstraintType (ConstraintType myConstraintType) 
{
    constraintType = myConstraintType;
}

ConstraintType ConstraintClass::getConstraintType () const
{
    return constraintType ;
}

String ConstraintClass::constraintTypeString () const  
{  // const promises not to change the object, i.e. ConstraintClass
        if (constraintType == WeldToAtom) { return  "WeldToAtom" ;}
        else if (constraintType == WeldToGround) {  return  "WeldToGround" ;}
        else if (constraintType == CoupledCoordinate) {  return  "CoupledCoordinate" ;}
        else if (constraintType == Undefined) { return "Undefined" ;}
        else return " ERROR! ";
}


void  ConstraintClass::print() const {
        std::cout<<__FILE__<<":"<<__LINE__   // to here is fine
          <<" : Chain ID : "      <<getChain1()
          <<" Residue    ID: "    <<getResidueID1().outString()
          <<" atom name : "       <<getAtomName1()
          <<" : Chain ID2 : "     <<getChain2()
          <<" Residue    ID2: "   <<getResidueID2().outString()
          <<" atom name2 : "      <<getAtomName2()
          //<<" to Ground: "        <<getToGround()
          <<" constraintType : " << constraintTypeString()
          <<endl;
    };  




int IntLen(const char* cstr)
{
  int    k, n = 0;
  if (cstr)
  {
    n = strspn(cstr,spaces);
    cstr += n;
    if (*cstr == '-' || *cstr == '+')
        ++cstr, ++n;
    k = strspn(cstr,digits);
    n = k?n+k:0;
  }
  return n;
}
/// <int>::[spaces][+|-]<digits>[garbage]
bool isNumber(string      line)
{
  //std::string test("1234.56");
  std::istringstream inpStream((line));
  double inpValue = 0.0;
  if (inpStream >> inpValue)
  {
    std::cout<<__FILE__<<":"<<__LINE__<<" Decided that "<<line<<" IS a number."<<std::endl; return 1;
    // ... Success!!  test is a number.
  }
  else
  {
    std::cout<<__FILE__<<":"<<__LINE__<<" Decided that "<<line<<" is not a number."<<std::endl; return 0;
    // ... Failure!!  test is not a number.
  }
}


bool isFixed (const String putativeFixedFloat) { // This checks that the string represents a floating point number in fixed format .. no scientific notation or other stray characters.
	int dotCount = 0;
	for (size_t i = 0 ; i < putativeFixedFloat.length() ; i ++) {
                //std::cout<<__FILE__<<":"<<__LINE__<<" putativeFixedFloat[i] = >"<<putativeFixedFloat[i]<<"< "<<std::endl;
		if (!((String(putativeFixedFloat[i]).compare("0") == 0) || 
		   (String(putativeFixedFloat[i]).compare("1") == 0) || 
		   (String(putativeFixedFloat[i]).compare("2") == 0) || 
		   (String(putativeFixedFloat[i]).compare("3") == 0) || 
		   (String(putativeFixedFloat[i]).compare("4") == 0) || 
		   (String(putativeFixedFloat[i]).compare("5") == 0) || 
		   (String(putativeFixedFloat[i]).compare("6") == 0) || 
		   (String(putativeFixedFloat[i]).compare("7") == 0) || 
		   (String(putativeFixedFloat[i]).compare("8") == 0) || 
		   (String(putativeFixedFloat[i]).compare("9") == 0) || 
		   (String(putativeFixedFloat[i]).compare(".") == 0) || 
		   (String(putativeFixedFloat[i]).compare("+") == 0) || 
		   (String(putativeFixedFloat[i]).compare("-") == 0)))  {
	 	        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Found a character : >"<<putativeFixedFloat[i]<< "< which is not of [0-9].+- .. this is not a fixed/float!"<<endl);
			return false; // actually we shouldn't get to this line
		}
		if (((String(putativeFixedFloat[i]).compare("+") == 0) || (String(putativeFixedFloat[i]).compare("-") == 0)) && (i > 0)) {
	 	        MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have tried to use '+' or '-' somewhere other than at the beginning of the number string.  This is not allowed!"<<endl);
			return false; // actually we shouldn't get to this line
		}
		if  (String(putativeFixedFloat[i]).compare(".") == 0) {
			dotCount++;
			if (dotCount > 1) {
			    MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have tried to use more than one '.' .. This is not allowed!"<<endl);
			    return false; // actually we shouldn't get to this line
			}
		} 
	} // of for
	return true;
}



    // a recursive algorithm for reading an integer from a String.  This String may contain ints, user variables (begin with @), +, and -.  No whitespaces or additional characters should be in the String.
    int   myAtoI(  map<const String,double> myUserVariables,  const char* value){
        MMBLOG_FILE_FUNC_LINE(INFO, ""<<endl);
        size_t plusPosition  = String(value).find_last_of('+');
        size_t minusPosition = String(value).find_last_of('-');
        if ((plusPosition > minusPosition) && (plusPosition  != String::npos) )  minusPosition = String::npos;
        if ((plusPosition < minusPosition) && (minusPosition != String::npos) )  plusPosition  = String::npos;
        String baseIntegerString ;
        int          increment = -1111;
        //int          decrement = -1111;
        size_t lastPlusOrMinus = min(plusPosition,minusPosition);
        if ((lastPlusOrMinus != String::npos) && (lastPlusOrMinus != 0)) {
            baseIntegerString = String(value).substr(0, (lastPlusOrMinus + 0) ); // Put everything to the left of the last +/- into baseIntegerString
            String incrementString = (String(value).substr(lastPlusOrMinus+0,1000)); //  NOT adding 1 to lastPlusOrMinus means that the sign at lastPlusOrMinus goes with incrementString.
            MMBLOG_FILE_FUNC_LINE(INFO, " At this stage, we are adding >"<<baseIntegerString<<"< and >"<<incrementString<<"<"<<endl);
            increment = myAtoI(myUserVariables, incrementString.c_str() );
            MMBLOG_FILE_FUNC_LINE(INFO, incrementString<<" was interpreted as >"<<increment<<"<"<<endl);
        } else if (lastPlusOrMinus == 0) { // There is a leading + or -, and  this is the only +/- in the whole expression.
            baseIntegerString = String(value).substr(1, 1000); // Put everything to from position 1 onwards into baseIntegerString
            if (plusPosition == 0) {
                MMBLOG_FILE_FUNC_LINE(INFO, " Detected that the base string : >"<<String(value) <<"< has a leading \'+\'. "<<endl);
                // Trim the leading '+' and return the rest    
                return myAtoI(myUserVariables, baseIntegerString.c_str() );
                //cout<<__FILE__<<":"<<__LINE__<<" Interpreted >"<<baseIntegerString<<"< as "<<increment<<endl; 
	    } else if (minusPosition == 0) {
                MMBLOG_FILE_FUNC_LINE(INFO, " Detected that the base integer string : >"<<baseIntegerString<<"< has a leading \'-\'.  Inverting sign."<<endl);
                // Trim the leading '-' and return the negative
                return -myAtoI(myUserVariables, baseIntegerString.c_str() );
                //cout<<__FILE__<<":"<<__LINE__<<" Interpreted >"<<baseIntegerString<<"< as "<<increment<<endl; 
            }
        }
        else { // no + or - found.
            MMBLOG_FILE_FUNC_LINE(INFO, ""<<endl);
            if (!((increment == -1111 ))){// && (decrement == -1111 )  )) {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, " Unexplained error!"<<endl);
            }
            MMBLOG_FILE_FUNC_LINE(INFO, ""<<endl);
            baseIntegerString = String(value);
            increment = 0;
            //decrement = 0;
            int baseInteger;
                if ((baseIntegerString.substr(0,1)).compare("@") ==0) {
                    cout<<__FILE__<<":"<<__LINE__<<""<<endl;
                    if (myUserVariables.find(baseIntegerString.c_str()) == myUserVariables.end())
                    {
                        MMBLOG_FILE_FUNC_LINE(CRITICAL, ": Undefined user variable "<<value<<endl);
                    }
                    MMBLOG_FILE_FUNC_LINE(INFO, ""<<endl);
                    double  intCast   = double(int(myUserVariables[baseIntegerString.c_str()]));
                    MMBLOG_FILE_FUNC_LINE(INFO, ""<<endl);
                    double  doubleCast = double(myUserVariables[baseIntegerString.c_str()]);
                    MMBLOG_FILE_FUNC_LINE(INFO, ""<<endl);
                    MMBLOG_FILE_FUNC_LINE(INFO, " Read user variable "<<baseIntegerString.c_str()<<"  which is set to : "<<myUserVariables[baseIntegerString.c_str()]<<endl);
                    SimTK_ERRCHK_ALWAYS(( (intCast) == doubleCast  ) ,"[ParameterReader.cpp]","Expected an int and got a non-integer");
                    MMBLOG_FILE_FUNC_LINE(INFO, ""<<endl);
                    baseInteger = int(myUserVariables[baseIntegerString.c_str()]);
                    MMBLOG_FILE_FUNC_LINE(INFO, ""<<endl);
                }
                else if (isNumber(baseIntegerString.c_str()))
                {
                    double  intCast   = double(int(atof(baseIntegerString.c_str())));
                    double  doubleCast = double(atof(baseIntegerString.c_str()));
                    SimTK_ERRCHK_ALWAYS(( (intCast) == doubleCast  ) ,"[ParameterReader.cpp]","Expected an int and got a non-integer");
                    baseInteger = (atoi(baseIntegerString.c_str()));
                } else {
                    MMBLOG_FILE_FUNC_LINE(CRITICAL, " : What you have entered: >"<<baseIntegerString<<"< is neither a variable (starting with @) nor an explicit number."<<endl);
                }
            //cout<<__FILE__<<":"<<__LINE__<<" : Result of "<<value<<" is : " <<  baseInteger <<endl;
            return baseInteger;
        }

        int baseInteger = myAtoI(myUserVariables,baseIntegerString.c_str() ) ; 

        int finalInteger = baseInteger + increment ;//- decrement;
        MMBLOG_FILE_FUNC_LINE(INFO, " : Result of "<< value  <<" is : " << finalInteger <<endl);
        return finalInteger;
    }   




/*void printBiopolymerSequenceInfo(const Biopolymer & myBiopolymer) {
    for (int i = 0; i < myBiopolymer.getNumResidues(); i++) {
        cout<<__FILE__<<":"<<__LINE__<<" Residue type, number, and insertion code: "<<myBiopolymer.getResidue(ResidueInfo::Index(i)).getOneLetterCode() <<", "<<myBiopolymer.getResidue(ResidueInfo::Index(i)).getPdbResidueNumber()<<", "<<myBiopolymer.getResidue(ResidueInfo::Index(i)).getPdbInsertionCode()<<endl;
    }    
};*/


         String intToString(int i) {
		/*
                // if boost doesn't work, go back to the stringstream method:
                String s ; 
                std::stringstream out;
                out << i;
                s = out.str();
                return s;
		*/
		//return  boost::lexical_cast<String>(i);      
                return SimTK::String(i);
            };  

Vec3 ValidateVec3(Vec3 myVec3){

    if (! myVec3.isFinite()){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, " This Vec3 vector is infinite or not a number : "<<myVec3<<endl);
    }
    return myVec3;
    /*if      (std::isnan(myVec3[0])) {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" x-component "<<myVec3[0] <<" of Vec3 "<<myVec3 <<" is invalid. "<<endl; ErrorManager::instance.treatError();}
    else if (std::isnan(myVec3[1])) {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" y-component "<<myVec3[1] <<" of Vec3 "<<myVec3 <<" is invalid. "<<endl; ErrorManager::instance.treatError();}
    else if (std::isnan(myVec3[2])) {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" z-component "<<myVec3[2] <<" of Vec3 "<<myVec3 <<" is invalid. "<<endl; ErrorManager::instance.treatError();}
    else return myVec3;
*/

}



int ValidateInt (const int myInt) {
	
                /*if      (std::isnan(myInt))      {
			ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" The integer is not valid : "<<myInt<<endl;
			ErrorManager::instance.treatError();
		}
                else if      (std::isinf(myInt)) {
			ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" The integer is not valid : "<<myInt<<endl;
			ErrorManager::instance.treatError();
		}
		else*/ return myInt;
}



int ValidateNonNegativeInt (const int myInt) {
	
        if (!(myInt >= 0)) {
			MMBLOG_FILE_FUNC_LINE(CRITICAL, " Expected a nonnegative integer. "<<endl);
		}
		else return myInt;
}




double ValidateNonNegativeDouble(const double myDouble) {
	
    if (std::isnan(myDouble)) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, " Not a number! "<<endl);
    }
    if (std::isinf(myDouble)) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, " Not a number! "<<endl);
    }
    if (!(myDouble>= 0)) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, " Expected a nonnegative Double   . "<<endl);
    }
    return myDouble;
}

double ValidateDouble(const double myDouble) {
	
    if (std::isnan(myDouble) || std::isinf(myDouble)) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, " Not a number! "<<endl);
    }
    return myDouble;
}


Real DotProduct(const Vec3 vecA, const Vec3 vecB) {
	ValidateVec3(vecA);
	ValidateVec3(vecB);
	return (vecA[0]*vecB[0] + vecA[1]*vecB[1] + vecA[2]*vecB[2] );
}


vector<String> readAndParseLine   (ifstream & inFile) {
        stringstream u;
	String inString;
 	getline (inFile,inString);
	//u.str(inString);
	//cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" read: "<<inString<<endl;
	vector<String> mystring;
	
        istringstream iss(inString);

    	do
    	{
        	String sub;
        	iss >> sub;
        	//cout << "Substring: " << sub << endl;
		mystring.push_back(sub);
    	} while (iss);

	return mystring;
}

String removeAllWhite (String &str)
{
    String temp;
    for (int i = 0; i < str.length(); i++)
        if (std::string(str)[i] != ' ') temp += std::string(str)[i];
    str = temp;
    return str;
}

vector<String> readAndParseOnColWidth   (ifstream & inFile, int columnWidth) {
        stringstream u;
	String inString;
 	getline (inFile,inString);
	//u.str(inString);
	//cout<<__FILE__<<":"<<__LINE__<<" read: "<<inString<<endl;
	vector<String> mystring;
	
        istringstream iss(inString);
    	//for (int j = 0; j < numColumns; j++)
        	String sub;
	int j = 0;
	do
    	{
		sub = inString.substr(j * columnWidth,columnWidth);
		removeAllWhite(sub);
        	//cout << "Substring: >" << sub <<"<"<< endl;
		if (sub.length()>0)  mystring.push_back(sub);			
		j++;
    	} while (sub.length()>0);

	return mystring;
}


    /*ParameterStringClass::ParameterStringClass( const String & paramsLine ){
        char * params = strdup( paramsLine.c_str() );
        char * token = strtok( params, " " );

        clear();
        while( token ){
            add( token );
            token = strtok( NULL, " " );
            std::cout<<__FILE__<<":"<<__LINE__<<" Added token : >"<<token<<"< "<<std::endl;
        }
        free( params );
    }*/

    void ParameterStringClass::validateNumFields(int correctNumFields) const{ // make sure we have the right number of parameters
        std::cout<<__FILE__<<":"<<__LINE__<<" This line contains "<<size()<< " elements. comparing to "<<correctNumFields<<"."<<endl;
        if ( size() < correctNumFields ){
            MMBLOG_FILE_FUNC_LINE(CRITICAL, ": You have not specified enough parameters for this command."<<endl);
        } else if ( size() > correctNumFields ) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, ": You have specified too many parameters for this command."<<endl);
        };  
    };  
    void ParameterStringClass::print() const {
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        for (int i = 0 ; i < size(); i++){
            MMBLOG_FILE_FUNC_LINE(INFO, stringVector[i]<<" ");
            //std::cout<<__FILE__<<":"<<__LINE__<<" "<<i<<" >"<<stringVector[i]<<"< "<<std::endl;
        };
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
    };

    String ParameterStringClass::getString() const {
        std::stringstream ss;
        for (int i = 0 ; i < size(); i++){
            ss <<" "<<stringVector[i];
        };

        return ss.str();
    }

