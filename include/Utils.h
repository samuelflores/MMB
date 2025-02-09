// vim: set sw=4 ts=4 sts=4 expandtab :
/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#ifndef Utils_H_
#define Utils_H_                 
#include <sys/stat.h>

#include <cstring>
#include <fstream>
#include <sstream>
#include <ctype.h>

#include "ExportMacros.h"
#include <cstddef>
#include "SimTKsimbody.h"
#include "SimTKmolmodel.h"
#include "MMBLogger.h"
#include <openmm/reference/RealVec.h>
#include <openmm/reference/ReferenceNeighborList.h>
#include <seqan/align.h>
#ifdef LEPTON_ENABLED
#include "Lepton.h"
#endif

typedef OpenMM::Vec3 openmmVecType ;

class BiopolymerClass;

using namespace SimTK;
using namespace std;  
// some functions to determine whether string is a number (http://www.tek-tips.com/viewthread.cfm?qid=1024751)
int IntLen(const char* cstr);
inline bool isInt(const char* cstr);
inline bool isInt(const string& s);
//bool isNumber(const char* cstr);
bool isNumber(std::string     );

bool isFixed (const String putativeFixedFloat) ; // This checks that the string represents a floating point number in fixed format .. no scientific notation or other stray characters.


struct LessThanComparator {
    template <typename T>
    static bool compare(const T &a, const T &b) {
        return a < b;
    }
    template <typename T>
    static bool invCompare(const T &a, const T &b) {
        return a > b;
    }
};

struct GreaterThanComparator {
    template <typename T>
    static bool compare(const T &a, const T &b) {
        return a > b;
    }
    template <typename T>
    static bool invCompare(const T &a, const T &b) {
        return a < b;
    }
};

// Usage: instantiate with CheckFile (fileName)
// Retrieve file info with accessors
class CheckFile {
private: 
    String fileName;
    struct stat st;
public:
    CheckFile(const String & myFileName);
    bool isDirectory();
    bool ownerCanRead();
    bool ownerCanWrite();
    void validateNonZeroSize();
    void validateExists();
    //void validateReadable();
};

int MMB_EXPORT checkOrCreateDirectory(const std::string & directoryPath);
int MMB_EXPORT myMkdir(const std::string & directoryPath);
int MMB_EXPORT myChdir(const std::string & directoryPath);
int MMB_EXPORT mySystemCall(const std::string & command);

enum class CopyFileResult {
    Success,
    Io_Error,
    No_Space,
};

CopyFileResult MMB_EXPORT mmbCopyFile(const std::string &sourceFileName, const std::string &destinationFileName);

void MMB_EXPORT closingMessage() ;


          String intToString(int i) ;

    // a recursive algorithm for reading an integer from a String.  This String may contain ints, user variables (begin with @), +, and -.  No whitespaces or additional characters should be in the String.
int   myAtoI(  map<const String,double> myUserVariables,  const char* value);


//static double myAtoF(map<const String,double> myUserVariables,const char* value );
// a recursive algorithm for reading a double from a String.  This String may contain ints, user variables (begin with @), +, and -.  No whitespaces or additional characters should be in the String.
inline double   myAtoF(  map<const String,double> myUserVariables,  const char* value){
    MMBLOG_FILE_FUNC_LINE(INFO, "inside myAtoF. converting string : >"<<value<<"<"<<endl);
#ifdef LEPTON_ENABLED

//#ifdef Lepton_USAGE
    MMBLOG_FILE_FUNC_LINE(DEBUG, " Using Lepton : "<<endl);
    map<string,double> leptonFormatUserVariables; // Wish this were not necessary. but userVariables uses the signature const SimTK::String,double . Lepton uses string, double.
    leptonFormatUserVariables.clear();
    for (auto  myUserVariablesIterator = myUserVariables.begin() ; myUserVariablesIterator !=myUserVariables.end(); myUserVariablesIterator++) {
        leptonFormatUserVariables[myUserVariablesIterator->first] = myUserVariablesIterator->second;
    }
    double leptonResult = Lepton::Parser::parse(std::string(value)).evaluate(leptonFormatUserVariables);
    MMBLOG_FILE_FUNC_LINE(INFO, " Lepton evaluation = >"<< leptonResult<<"< "<<std::endl);
    if (isnan(leptonResult)) MMBLOG_FILE_FUNC_LINE(CRITICAL, "The provided string >"<<value <<"< Evaluates to NaN! "<<endl);
    return leptonResult; // if Lepton_USAGE is defined, then we return here and the rest of the procedure is not used. 
#endif
    // If  LEPTON_ENABLED is NOT defined, then we parse the formula the old dumb way, as follows.
    MMBLOG_FILE_FUNC_LINE(DEBUG, " NOT using Lepton : "<<endl);

    size_t plusPosition  = String(value).find_last_of('+'); // returns the position of the last '+' in value
    size_t minusPosition = String(value).find_last_of('-'); // ditto for '-'
    if ((plusPosition > minusPosition) && (plusPosition  != String::npos) )  minusPosition = String::npos; // If the plus sign is closer to the end of the string, pretend we didn't find any '-'
    if ((plusPosition < minusPosition) && (minusPosition != String::npos) )  plusPosition  = String::npos; // Conversely, if the '-' is closer to the end, pretend we didn't find any '+' .. actually we might not have found any '+' anyway.
    String baseDoubleString ;
    double          increment = -1111;
    double          decrement = -1111;
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    if (plusPosition != String::npos) { // We have a '+' to deal with
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        baseDoubleString = String(value).substr(0, (plusPosition + 0) );
        String incrementString = (String(value).substr(plusPosition+1,1000)); // the second parameter is ridiculously large, but will be truncated at the end of the input String.
        //MMBLOG_FILE_FUNC_LINE(" The increment String is : "<<incrementString<<endl;
        stringstream incrementStringStream(incrementString);
        increment = myAtoF(myUserVariables, incrementString.c_str() );
        decrement = 0;
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
    } else if (minusPosition != String::npos ){ // we have a '-' to deal with
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        if (minusPosition== 0) {
            // If this is just a leading '-' sign, then put a zero to the left of that minus.
            MMBLOG_FILE_FUNC_LINE(INFO, "Detected a leading \'-\' sign. Will insert a zero to the left of the \'-\'."<<endl);
            baseDoubleString = "0.0";
        }
        else {
            // Otherwise, parse whatever is to the left of the minus sign:
            baseDoubleString = String(value).substr(0, (minusPosition + 0) ); }
        MMBLOG_FILE_FUNC_LINE(INFO, "baseDoubleString =  >"<<baseDoubleString  <<"< "<<endl);

        String decrementString = (String(value).substr(minusPosition+1,1000)); // the second parameter is ridiculously large, but will be truncated at the end of the input String.
        stringstream decrementStringStream(decrementString);
        MMBLOG_FILE_FUNC_LINE(INFO, "About to extract numerical decrement from the string >"<<decrementString<<"< "<<endl);
        decrement = myAtoF(myUserVariables, decrementString.c_str() );
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        increment = 0;
        //MMBLOG_FILE_FUNC_LINE(endl;
    } else { // no + or - found. This means we can return a result without further recursion.. i.e. we are at a leaf of the recursion tree.
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        if (!((increment == -1111 ) && (decrement == -1111 )  )) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unexplained error!"<<endl);
        }
        baseDoubleString = String(value);
        increment = 0;
        decrement = 0;
        double baseDouble;
        {
            MMBLOG_FILE_FUNC_LINE(INFO, endl);
            if ((baseDoubleString.substr(0,1)).compare("@") ==0) {
                MMBLOG_FILE_FUNC_LINE(INFO, endl);
                if (myUserVariables.find(baseDoubleString.c_str()) == myUserVariables.end())
                {
                    MMBLOG_FILE_FUNC_LINE(CRITICAL, "Undefined user variable "<<value<<endl);
                }
                MMBLOG_FILE_FUNC_LINE(INFO, "Read user variable "<<baseDoubleString.c_str()<<"  which is set to : "<<myUserVariables[baseDoubleString.c_str()]<<endl);
                baseDouble = double(myUserVariables[baseDoubleString.c_str()]);
            }
            else {
                // This has no '+', '-', or leading '@'.  However it's still possible that the user gave scientific notation, so e.g. 1e-5 leads here to baseDoubleString = '1e'.  So we need to make sure there is nothing but [0-9],'+','-' in this string:
        if (isFixed(baseDoubleString)) {
                    baseDouble = (atof(baseDoubleString.c_str()));
                        MMBLOG_FILE_FUNC_LINE(INFO, "We appear to be in a leaf of the recursion tree for myAtoF. Parsed >"<<baseDoubleString<< "< as : "<<baseDouble<<endl);
        } else {
                        MMBLOG_FILE_FUNC_LINE(CRITICAL, "There was an error processing a putative floating point number."<<endl);
        }
            }
        }
        //MMBLOG_FILE_FUNC_LINE(endl;
	if (isnan(baseDouble)) MMBLOG_FILE_FUNC_LINE(CRITICAL, "The provided string >"<<value<<"< Evaluates to NaN! "<<endl);
        return baseDouble;
    }

    //MMBLOG_FILE_FUNC_LINE(endl;
    double baseDouble = myAtoF(myUserVariables,baseDoubleString.c_str() ) ;

    double finalDouble = baseDouble + increment - decrement;
    MMBLOG_FILE_FUNC_LINE(INFO, "Result of >"<< value  <<"< is : " << finalDouble <<endl);
    return finalDouble;
}



//static bool aToBool( const char* value );
inline bool aToBool(  const char* value ) {

    String upperValue(value);
    for(int i=0;i<(int)upperValue.length();i++)  {
        upperValue[i] = toupper(value[i]);
    }

    if (( upperValue ==  "TRUE" ) ||( upperValue ==  "1")) {
        MMBLOG_FILE_FUNC_LINE(DEBUG, "TRUE"<<endl);
        return true;
    }
    else if (( upperValue ==  "FALSE" ) ||( upperValue ==  "0")){
        MMBLOG_FILE_FUNC_LINE(DEBUG, "FALSE"<<endl);
        return false;
    }
    else {
        MMBLOG_FILE_FUNC_LINE(INFO, "Error -- you have specified"<<value<<endl);
        SimTK_ERRCHK_ALWAYS((upperValue == "TRUE" || upperValue == "FALSE" || upperValue == "1"  || upperValue == "0") ,"[ParameterReader.cpp]"," requires either True or False but was set to something else");
        return false;
    }

}
//static bool aToBool( const String& name, const char* value );
inline bool aToBool( const String& name, const char* value ) {
    return aToBool(value);
}

inline bool compareUpper( const String& param, const char* symbol ) {

    String upperParam(param);
    String upperSym(symbol);

    if( upperParam.length() != upperSym.length() ) return false;


    for(int i=0;i<(int)upperParam.length();i++)  {
        upperParam[i] = toupper(param[i]);
        upperSym[i] = toupper(symbol[i]);
    }

    if( upperParam ==  upperSym )
        return true;
    else
        return false;
}



namespace { // unnamed namespace means it is translation unit local.
    const char* spaces = " \t";
    const char* digits = "0123456789";
}



std::string   trim(const std::string& str,
                 const std::string& whitespace = " \t\n\r"
                 );

bool MMB_EXPORT vectorCompare(String myString, vector<String> & comparisonStringVector) ;



BondMobility::Mobility stringToBondMobility(const String bondMobilityString);

class ConstraintFunction : public Function {
     Real calcValue(const Vector& x) const {
         return x[0]-x[1];
     }   
     Real calcDerivative(const Array_<int>& derivComponents,const Vector& x) const {
         if (derivComponents.size() == 1)
             return derivComponents[0] == 0 ? 1 : -1; 
             return 0;
         }   
     int getArgumentSize() const {
         return 2;
     }   
     int getMaxDerivativeOrder() const {
         return std::numeric_limits<int>::max();
     }   
};


class ResidueID {
public:
    int ResidueNumber;
    char InsertionCode;
public:
    
    ResidueID() { ResidueNumber = 0;  InsertionCode = ' ' ;};
    ResidueID(int residueNumber) { ResidueNumber = residueNumber; InsertionCode = ' '          ;};
    ResidueID(int residueNumber, char insertionCode) { ResidueNumber = residueNumber; InsertionCode = insertionCode;};
    ResidueID( map<const String,double> myUserVariables,  const char* value) { // Allow the user to provide a user-configured variable, starting with '@'
	MMBLOG_FILE_LINE(DEBUG, std::endl);
        ResidueNumber = myAtoI( myUserVariables, value );
        InsertionCode = ' ';    
    }
    ResidueID(String inputString/*, bool validate = true*/) {
        char insertionCode;
        int residueNumber;
        int stringLength = inputString.length();
        //char endChar = *(inputString.substr(stringLength-1, 1)).c_str();
        //if (validate)
            for (int i = 1; i < stringLength-1 ; i++) { // check all characters except first and last.
                if (
                    isdigit(inputString.substr(i,1).c_str()[0]) //|| // if the character is a digit
                    //((((inputString.substr(i,1).c_str()))[0]) ==     ('.')  ) // if the character is a decimal point
                    
                    ) { // everything is OK!  
                    } else {
                        MMBLOG_FILE_LINE(CRITICAL, " Error!  You have an unallowed character \'"<<inputString.substr(i,1).c_str()[0]<< "\' in your residue ID : \'"<<inputString<<endl);
                }
            }
        if (
            isdigit((inputString.substr(stringLength-1, 1)).c_str()[0] )   // if last character in the string is 0-9.
            ){
            insertionCode =  ' ';
            residueNumber =   atoi((inputString ).c_str())   ; 
        }
        else {  // last character is not 0-9, but rather another character.
                if (stringLength < 2) {  // if the last char is an insertion code, the string should be at least two characters long -- at least one digit for the residue number, even if it's 0.
                    MMBLOG_FILE_LINE(CRITICAL, " : Something is wrong!  The last character looks like an insertion code, but there are no digits preceding it!  Valid ResidueID's look like e.g. 1 , 10B , 121P.  Not 10 B , Y , etc."<<std::endl);
                }

        if (
            !isdigit((inputString.substr(stringLength-2, 1)).c_str()[0] )   // if last character in the string is NOT 0-9.
            ){
                    MMBLOG_FILE_LINE(CRITICAL, " : Something is wrong!  The last character looks like an insertion code, but the penultimate one is not a digit!  Valid ResidueID's look like e.g. 1 , 10B , 121P.  Not 10 B , Y , etc."<<std::endl);
                }
            insertionCode =  *(inputString.substr(stringLength-1, 1).c_str());
            residueNumber =   atoi( ((inputString.substr(0, stringLength-1) ).c_str()));
        } 
            ResidueNumber = residueNumber;
            InsertionCode = insertionCode;   
        //return ResidueID(residueNumber,insertionCode);
    };
    void setResidueNumber(int myResidueNumber){
        ResidueNumber = myResidueNumber;
    };

    void setInsertionCode(char myInsertionCode){
        InsertionCode = myInsertionCode;
    };

    int getResidueNumber() const{
        return ResidueNumber;
    };

    char getInsertionCode() const {
        return InsertionCode;
    };
    const String outString( ) const {
        stringstream  myStringStream;
        if (InsertionCode != ' ')
            myStringStream << ResidueNumber<<InsertionCode;
        else 
            myStringStream << ResidueNumber;
        return myStringStream.str();
    };
    // This polymorphism of outString takes a chain ID (expected to be 1 character long) and concatenates it with the ResidueID to create a 6-character string suitable to be fed into columns 22-27 (chain ID, residue number, insertion code) of a PDB atom record 
    const String chainIDResidueID(const String &chainID) const {
        stringstream  myStringStream; myStringStream.clear();
        int totalWidth = 6;
        // Columns 23-26 are for the residue number, 27 is for the insertion code.  So 5 characters.
        // setw only affects the NEXT value (ResidueNumber).  hence the 4 here:
        myStringStream << std::setw(totalWidth-2)<< ResidueNumber<<InsertionCode;
        if ((myStringStream.str()).length() != (totalWidth - 1)) {
                MMBLOG_FILE_LINE(CRITICAL, " : Something is wrong!  The residue number + insertion code should be padded to be exactly 5 characters long. You have: >"<<myStringStream.str() <<"< " <<std::endl);
        }
        String outString = chainID + myStringStream.str();
        if ((outString).length() != totalWidth) {
                MMBLOG_FILE_LINE(CRITICAL, " : Something is wrong!  The chain ID + residue number + insertion code string should be padded to be exactly 6 characters long. You have: >"<<outString <<"< " <<std::endl);
        }
        return outString;
        /*else 
            myStringStream << ResidueNumber;*/
        //return myStringStream.str();
    };

    bool operator==(const ResidueID &other) const {
        return (ResidueNumber == other.ResidueNumber) && (InsertionCode == other.InsertionCode);
    }

    bool operator!=(const ResidueID &other) const {
        return !(*this == other);
    }
    // Inequalities don't work , because we no longer require residue ID's to be monotonically increasing, especially wrt insertion codes.
    
    bool operator > (const ResidueID &other) const {
    //...  // Compare the values, and return a bool result.
    if ((ResidueNumber > other.ResidueNumber) )
        {return true;}
    if ( ResidueNumber == other.ResidueNumber ) 
        if ( InsertionCode > other.InsertionCode ){return true;}
    return false;  
    };

    bool operator < (const ResidueID &other) const {
    //...  // Compare the values, and return a bool result.
    if ((ResidueNumber < other.ResidueNumber) )
        {return true;}
    if ( ResidueNumber == other.ResidueNumber ) 
        if ( InsertionCode < other.InsertionCode ){return true;}
    return false;  
    };

    bool operator >= (const ResidueID &other) const {
    if ((*this == other) || (*this > other)) {return true;}
    else {return false;}
    }

    bool operator <= (const ResidueID &other) const {
    if ((*this == other) || (*this < other)) {return true;}
    else {return false;}
    }

};

class BiopolymerModification {
private:
    String    chainToModify;           // Biopolymer chain to be modified.
    ResidueID residueToModify;         // Residue number of residue to be modified.
    String    atomToModifyOnBiopolymer;// Side chain atom to which compound will be attached. Should have one free bond for this. 
    String    freeBondOnBiopolymer;    // Free bond on biopolymer to which the compound will be attached
    String    chainToAdd;           // Compound to add to the biopolymer. Should have one free bond.
    String    atomOnAddedCompound;     // This is the atom on the added compound, which will receive a covalent bond from the biopolymer. expects a distinctive atom name, with no residue number, e.g. SG or 1HB.
    String    freeBondOnAddedCompound; // This is the bond name (e.g. bond1, bond2, ..) on the atom on the added compound which will receive the covalent bond from the biopolymer. 
public:
    void setChainToModify (String myChainToModify) {chainToModify = myChainToModify; }   
    void setResidueToModify (ResidueID myResidueToModify) {residueToModify = myResidueToModify; }   
    void setAtomToModifyOnBiopolymer (String myAtomToModifyOnBiopolymer) {atomToModifyOnBiopolymer = myAtomToModifyOnBiopolymer; }   
    void setFreeBondOnBiopolymer (String myFreeBondOnBiopolymer) {freeBondOnBiopolymer = myFreeBondOnBiopolymer; }   
    void setChainToAdd (String   myChainToAdd) {chainToAdd = myChainToAdd; }   
    void setAtomOnAddedCompound (String myAtomOnAddedCompound) {atomOnAddedCompound = myAtomOnAddedCompound; }   
    void setFreeBondOnAddedCompound (String myFreeBondOnAddedCompound) {freeBondOnAddedCompound = myFreeBondOnAddedCompound; }   
    String    getChainToModify ()      const {return chainToModify; }   
    ResidueID getResidueToModify ()    const {return residueToModify; }   
    String    getChainToAdd ()         const {return chainToAdd; }   
    String    getAtomOnAddedCompound() const {return atomOnAddedCompound; }   
};

struct ChainResidueToGround {
    String chainID;
    ResidueID residueID;
    bool toGround;
};

struct strCmp {
    bool operator()( const char* s1, const char* s2 ) const {
        return strcmp( s1, s2 ) < 0;
    }
};

class ParameterStringClass {
private:
    vector <String> stringVector;
    vector <String> getStringVector(){return stringVector;};
public:
    ParameterStringClass( const String & paramsLine ){
        //std::cout<<__FILE__<<":"<<__LINE__<<" paramsLine prior to trimming : >"<<paramsLine<<"<"<<std::endl; 
        //const String myNewParamsLine = myTrim( paramsLine);
        //trim ( paramsLine);
        //std::cout<<__FILE__<<":"<<__LINE__<<" paramsLine after trimming : >"<<paramsLine<<"<"<<std::endl; 
        char * params = strdup( paramsLine.c_str() );
        char * token = strtok( params, " \t\n\r" ); 

        clear();
        while( token ){
            add( token );
            token = strtok( NULL, " \t\n\r" ); 
            //std::cout<<__FILE__<<":"<<__LINE__<<" Added token : >"<<token<<"< "<<std::endl;
        }
        free( params );
    }

    void clear(){stringVector.clear();};
    int size() const {return stringVector.size();}

   // ParameterStringClass(const ParameterStringClass&){}; // copy 
    //ParameterStringClass& operator=(const ParameterStringClass&){}; //  copy.
    void print() const;
    String getString() const;
    String getString(int stringIndex) const {
        if ( stringIndex >= size() ) {
            return "";
        }
        return stringVector[stringIndex];
    }
    void add(String tempString){
        stringVector.push_back(tempString); 
    }
    void padStringVector(int numFields) {
        clear();
        for (int w = 0;w<numFields;w++) stringVector.push_back("");
    };
    void validateNumFields(int correctNumFields) const; 

};

class  ChainResidueIndex {
public:
    int residueIndex;
    int chainIndex;
    ChainResidueIndex(int myChainIndex,  int myResidueIndex) ;
};
struct SingleBondMobility {
    String chain1;
    ResidueID    residue1;
    String atom1   ;
    String mobility;
    String chain2  ;
    ResidueID    residue2;
    String atom2   ;
};


struct twoIndexCmp {
    bool operator()( const ChainResidueIndex ti1, const ChainResidueIndex ti2 ) const {
        if (ti1.chainIndex < ti2.chainIndex) return 1;
        else if ((ti1.chainIndex == ti2.chainIndex) && (ti1.residueIndex < ti2.residueIndex)) return  1;
        else return 0;
    }
};

struct Displacement {
    String   chain;
    Vec3     displacement;
    Rotation rotation;      
};

struct BasePairPartner   {
    String BPEdge;
    ResidueID    BPPartner;
    String BPPartnerChain;
    int    rigidSegmentNumber;
    ResidueID    rigidSegmentEndNumber;
    ResidueID    rigidSegmentStartNumber;
};



struct BasePair   {
               String FirstBPEdge           ;
               String SecondBPEdge          ;
               String FirstBPChain          ;
               String SecondBPChain         ;
               ResidueID    FirstBPResidue        ;
               ResidueID    SecondBPResidue       ;
               String OrientationBP         ;
               String BasePairIsTwoTransformForce;
               int    BasePairPriority      ;
               int    BasePairTemporary     ;
               Rotation rotationCorrection1 ;
               Rotation rotationCorrection2 ;
               Vec3   translationCorrection1;
               Vec3   translationCorrection2;
               bool   basePairSatisfied     ;
               int    leontisWesthofBondRowIndex;
               
    };

struct BaseInteraction   {
               String FirstBPEdge           ;
               String SecondBPEdge          ;
               String FirstBPChain          ;
               String SecondBPChain         ;
               ResidueID    FirstBPResidue        ;
               ResidueID    SecondBPResidue       ;
               String OrientationBP         ;
               //String BasePairIsTwoTransformForce;
               //int    BasePairPriority      ;
               //int    BasePairTemporary     ;
               Rotation rotationCorrection1 ;
               Rotation rotationCorrection2 ;
               Vec3   translationCorrection1;
               Vec3   translationCorrection2;
               bool   basePairSatisfied     ;
               int    leontisWesthofBondRowIndex;
               
    };

class     NTC_Classes {
    public:
    // NTC
               String NtC_FirstBPChain;
               ResidueID FirstBPResidue;
               ResidueID SecondBPResidue;
               int NtC_step_ID;
               String NtC_Class_String; 
               int    NtC_INDEX;
               double Confalparam;
               Rotation rotationCorrection1;
               Rotation rotationCorrection2;
               Vec3   translationCorrection1;
               Vec3   translationCorrection2;
               int    NTC_PAR_BondRowIndex;
               double weight,weight2;
               int    meta = 0;
               int    count = 0;
               BiopolymerClass *bp;
               std::array<Compound::AtomIndex, 4> atomIndices;

    void print(){
        std::stringstream ss{};
        for (const auto &i : atomIndices) {
            ss << i << ", ";
        }
        ss << std::endl;

        MMBLOG_PLAIN(INFO, " Printing parameters for NtC:"<<std::endl
        << " NtC_FirstBPChain = >"<<NtC_FirstBPChain<<"<"<<std::endl
        << " FirstBPResidue = >"<<FirstBPResidue.outString()<<"<"<<std::endl
        << " SecondBPResidue = >"<<SecondBPResidue.outString()<<"<"<<std::endl
        << " NtC_step_ID = >"<<NtC_step_ID<<"<"<<std::endl
        << " NtC_Class_String = >"<<NtC_Class_String<<"<"<<std::endl
        << " NtC_INDEX = >"<<NtC_INDEX<<"<"<<std::endl
        << " rotationCorrection1 = >"<<rotationCorrection1<<"<"<<std::endl
        << " rotationCorrection2 = >"<<rotationCorrection2<<"<"<<std::endl
        << " translationCorrection1 = >"<<translationCorrection1<<"<"<<std::endl
        << " translationCorrection2 = >"<<translationCorrection2<<"<"<<std::endl
        << " NTC_PAR_BondRowIndex = >"<<NTC_PAR_BondRowIndex<<"<"<<std::endl
        << " meta = >"<<meta<<"<"<<std::endl
        << " count = >"<<count<<"<"<<std::endl
        << " Molmodel atom indices" << ss.str() << std::endl
        << std::endl);
    }
               
};
               
struct IncludeIntraChainInterface {
               String Chain;
               double Depth;
};


class Interface {
    public:
        vector<String> Chains ;
               vector<String> PartnerChains ;
               double    Depth       ;
               String MobilizerString;
               vector<String> getChains() {return Chains;};
               vector<String> getPartnerChains() {return PartnerChains;};
               double         getDepth() {return Depth;};
               void print (){
                   MMBLOG_FILE_LINE(INFO, "For interface, printing reference chains : "<<std::endl);
                   for (unsigned int i = 0; i < Chains.size(); i ++) {
                        MMBLOG_FILE_LINE(INFO, ">"<<Chains[i]<<"<, "<<std::endl);
                   };
                   MMBLOG_FILE_LINE(INFO, " Partner chains : "<<std::endl);
                   for (unsigned int i = 0; i < Chains.size(); i ++) {
                        MMBLOG_FILE_LINE(INFO, ">" << PartnerChains[i]<<"<, "<<std::endl);
                   }
                   MMBLOG_FILE_LINE(INFO, " Depth : "<<Depth<<std::endl);
               };
};

class ResidueStretch   {
    private:
                SimTK::String chain          ;
                ResidueID    startResidue        ;
                ResidueID    endResidue       ;
    public:
                SimTK::String getChain() const {return chain;};
                ResidueID getStartResidue()const {return startResidue;};
                ResidueID getEndResidue()const {return endResidue;};
                void setChain(SimTK::String myChain) {chain = myChain;};
                void setStartResidue(ResidueID myStartResidue) {startResidue = myStartResidue;};
                void setStartResidueNumber(int myStartResidueNumber){startResidue.setResidueNumber(myStartResidueNumber); }
                void setEndResidueNumber  (int myEndResidueNumber)  {  endResidue.setResidueNumber(  myEndResidueNumber); }
                void setEndResidue(ResidueID myEndResidue){endResidue = myEndResidue;};
                bool sameParameters(ResidueStretch myResidueStretch) { // compares all parameters except startResidue and startResidue, returns False if any differ.  This should be overridden by daughter classes, particularly if they have additional parameters which may differ.
                    if (getChain().compare(myResidueStretch.getChain() ) == 0 ) return true;
                    else return false;
                };
                void printStretch() const {
                    MMBLOG_FILE_LINE(INFO, "Stretch chain ="<<getChain()<<", first residue ="<<getStartResidue().outString()<<", last residue ="<<getEndResidue().outString()<<std::endl);
                }
                ResidueStretch() {
                    setStartResidue(ResidueID("-1111" ) );
                    setEndResidue(ResidueID("-1111" ) );
                    setChain( SimTK::String (" ") );
                }

                ResidueStretch(const SimTK::String & myChain, const ResidueID & myStartResidue, const ResidueID & myEndResidue) { 
                    chain = myChain; startResidue = myStartResidue; endResidue= myEndResidue;
		    printStretch();
                }

                ResidueStretch(const SimTK::String & myChain, const ResidueID & myResidue)
                    { 
                    chain = myChain; 
		    startResidue = myResidue; 
		    endResidue= myResidue; }; // This constructor sets endResidue = startResidue
                bool contains(SimTK::String myChain, ResidueID resID) const
                {
                    return chain == myChain && resID >= startResidue && resID <= endResidue;
                }
        bool operator == (ResidueStretch & a){
            if ((this->startResidue == a.startResidue ) &&
                (this->endResidue == a.endResidue ) &&
                (this->chain.compare(a.chain) == 0 ) )
            {return true;}
            else return false;
        }
        // Will need to replace these obsolete operators (ResidueID does not take inequalities any more) with new error traps

        template <typename Comparator>
        bool compareOp (const ResidueStretch & a) const {
            if (this->startResidue > this->endResidue) {
                MMBLOG_FILE_FUNC_LINE(
                    CRITICAL,
                    "The current residue stretch has a start residue : "<<this->startResidue.outString()
                    <<" which is greater than its end residue: "<<this->endResidue.outString()<<std::endl
                );
            }

            if (a.startResidue > a.endResidue) {
                MMBLOG_FILE_FUNC_LINE(
                    CRITICAL,
                    "The current residue stretch has a start residue : "<<a.startResidue.outString()
                    << " which is greater than its end residue: "<<a.endResidue.outString()<<std::endl
                );
            }

            if (Comparator::compare(this->chain, a.chain))
                return true;
            else if (Comparator::invCompare(this->chain, a.chain))
                return false;
            else
                return Comparator::compare(this->startResidue, a.endResidue);
        }

        bool operator > (const ResidueStretch & a) const {
            return compareOp<GreaterThanComparator>(a);
        }

        bool operator < (const ResidueStretch & a) const {
            return compareOp<LessThanComparator>(a);
        }
    };

class SingleResidue : public ResidueStretch  {
        public:
                ResidueID getEndResidue(){
                        MMBLOG_FILE_LINE(CRITICAL, " getEndResidue() is not compatible with SingleResidue! "<<std::endl);
                        return ResidueID(-1111,' ');
                }
                ResidueID  getResidue() const // The second const promises that the method will not change 'this' object
                    { return getStartResidue();} 
                //void setStartResidue(ResidueID myResidue) {startResidue = myResidue; endResidue = myResidue;};
                void setResidue(ResidueID myResidue) {setStartResidue ( myResidue);setEndResidue ( myResidue); }; // This sets both startResidue and endResidue to myResidue, even though endResidue will never be retrieved.
                void setResidueNumber(int       myResidueNumber) {   setStartResidueNumber ( myResidueNumber);   setEndResidueNumber ( myResidueNumber); }; // This sets both startResidue and endResidue to myResidue, even though endResidue will never be retrieved.
                void printStretch() {MMBLOG_FILE_LINE(INFO, "Stretch chain ="<<getChain()<<", first (and only) residue ="<<getStartResidue().outString()<<std::endl);  };
		bool operator == (SingleResidue & a){ // operators appear not to be inherited, perhaps related to inability to access private members of ResidueStretch
		    if ((this->getResidue() == a.getResidue() ) &&
			(this->getChain().compare(a.getChain()) == 0 ) )
		    {return true;}
		    else return false;
		}
		bool operator == (const SingleResidue & a) const { // operators appear not to be inherited, perhaps related to inability to access private members of ResidueStretch
		    if ((this->getResidue() == a.getResidue() ) &&
			(this->getChain().compare(a.getChain()) == 0 ) )
		    {return true;}
		    else return false;
		}
                
		bool operator != (const SingleResidue & a) const { // operators appear not to be inherited, perhaps related to inability to access private members of ResidueStretch
		    if ((this->getResidue() == a.getResidue() ) &&
			(this->getChain().compare(a.getChain()) == 0 ) )
		    {return false;}
		    else return true;
		}
    };



class MMB_EXPORT MobilizerStretch : public ResidueStretch {
        private:
               SimTK::BondMobility::Mobility BondMobility;
               String BondMobilityString;
        public:
               SimTK::BondMobility::Mobility getBondMobility() const {
                   return BondMobility;
               };
	       bool bondMobilityIsRigid() const {
                   return BondMobility == SimTK::BondMobility::Mobility::Rigid;
	       }
               SimTK::BondMobility::Mobility setBondMobility(SimTK::String myBondMobilityString) {
                   MMBLOG_FILE_LINE(INFO, " About to set BondMobility to >"<<myBondMobilityString<<"< "<<std::endl);
                   BondMobilityString = myBondMobilityString;
                   BondMobility = stringToBondMobility(myBondMobilityString);

                   return BondMobility;
               };
               MobilizerStretch(){setStartResidue ( ResidueID()); setEndResidue ( ResidueID()); setChain ( ""); setBondMobility ("Default");};

               MobilizerStretch(ResidueStretch myResidueStretch, SimTK::String myBondMobilityString = "Default"){
                   setStartResidue(myResidueStretch.getStartResidue());
                   setEndResidue(myResidueStretch.getEndResidue());
                   setChain(myResidueStretch.getChain());
                   setBondMobility( myBondMobilityString);
               };

               MobilizerStretch(const SimTK::String & myChain, const ResidueID & myStartResidue, const ResidueID & myEndResidue, const SimTK::String & myBondMobilityString = "") {
                        setChain ( myChain); setStartResidue ( myStartResidue), setEndResidue  (myEndResidue);
                        setBondMobility(myBondMobilityString);

               };
            bool operator == (const MobilizerStretch & a)
            {
                return ((this->getStartResidue() == a.getStartResidue() ) &&
                   (this->getEndResidue() == a.getEndResidue() ) &&
                   (this->getBondMobility() == a.getBondMobility() ) &&
                   (this->getChain().compare(a.getChain()) == 0 ) );
            }
            bool operator == (MobilizerStretch & a)
            {
                return ((this->getStartResidue() == a.getStartResidue() ) &&
                   (this->getEndResidue() == a.getEndResidue() ) &&
                   (this->getBondMobility() == a.getBondMobility() ) &&
                   (this->getChain().compare(a.getChain()) == 0 ) );
            }
        const String & getBondMobilityString() const {
            return BondMobilityString;
        };

        void print() const {
            MMBLOG_FILE_LINE(INFO, " Printing MobilizerStretch. getBondMobilityString() : >"<<getBondMobilityString()<< "< getBondMobility() >"<<getBondMobilityString()<<"<"<<std::endl);
            printStretch();
        } 

    };

typedef SingleResidue IncludeAllNonBondAtomsInResidue;
typedef SingleResidue IncludeResidue;
/*struct IncludeAllNonBondAtomsInResidue {
    String chain;
    ResidueID    residue;

    bool operator==(const IncludeAllNonBondAtomsInResidue &other) const {
        // Compare the values, and return a bool result.
        return chain == other.chain && residue == other.residue;
    };
};*/

class MMB_EXPORT AllResiduesWithin: public  SingleResidue {
    public:
       //String chain;
       //ResidueID    residue;
       double radius;
       AllResiduesWithin () {
	    setResidue(ResidueID("-1111" ) );
	    setChain( SimTK::String (" ") );
            setRadius (0.0);
       };
       AllResiduesWithin (SimTK::String myChain, ResidueID myResidue, double myRadius) {
	    setResidue(myResidue );
	    setChain( myChain );
            setRadius (myRadius);
       };
       void setRadius(double myRadius) {radius = myRadius;}
       double getRadius() const {return radius;}
        void print () const { MMBLOG_FILE_LINE(DEBUG, " I am an AllResiduesWithin object. chain, residue, radius = "<<getChain()<<", "<<getResidue().outString()<<", "<<getRadius()<<std::endl);};

};

class MobilizerWithin: public AllResiduesWithin {
    private:   
               String BondMobilityString           ;
    public:
        void   setBondMobilityString(String myBondMobilityString){BondMobilityString = myBondMobilityString;}
        String getBondMobilityString(){return BondMobilityString;}
};

struct MobilizerDomainsInterface {
    ResidueStretch domain1;
    ResidueStretch domain2;
    double    range;
    String MobilizerString;
    bool rigidBackbone;
};

class TwoAtomClass       {
               //String ConstraintScheme           ;
               //String Chain1          ;
    protected:
               ResidueID    residueID1        ;
               String chain1;
               String atomName1;
               ResidueID    residueID2        ;
               String chain2;
               String atomName2;
               double distance;
    public:
    TwoAtomClass(){
        chain1 = ""; residueID1 = ResidueID(); atomName1 = ""; 
        chain2 = ""; residueID2 = ResidueID(); atomName2 = ""; distance = -1111;};
    TwoAtomClass(String myChain, ResidueID inputResidueID,String myAtomName) {
        residueID1 = (inputResidueID);
        atomName1 = myAtomName;
        chain1 = myChain;
        residueID2 = ResidueID();
        atomName2 = "" ;        
        chain2 = "";
        distance = -1111;
    };
    TwoAtomClass(String myChain, ResidueID inputResidueID,String myAtomName,String myChain2, ResidueID inputResidueID2,String myAtomName2, double myDistance = -1111) {
        residueID1 = (inputResidueID);
        //residueID1.setInsertionCode ( residueID1.getInsertionCode());
        atomName1 = myAtomName;
        chain1 = myChain;
        residueID2 = (inputResidueID2);
        atomName2 = myAtomName2;
        chain2 = myChain2;
        distance = myDistance;
    };
    void  setChain1(String myChain) {chain1     = myChain;}
    void  setResidueID1(ResidueID myResidueID) {residueID1 = myResidueID;}
    void  setAtomName1(String myAtomName    ) {atomName1 = myAtomName;}

    void  setChain2(String myChain) {chain2     = myChain;}
    void  setResidueID2(ResidueID myResidueID) {residueID2 = myResidueID;}
    void  setAtomName2(String myAtomName    ) {atomName2 = myAtomName;}
    const ResidueID getResidueID1() const {return residueID1;};
    const String getAtomName1() const {return atomName1;};
    const String getChain1() const {return chain1;};
    const ResidueID getResidueID2() const {return residueID2;};
    const String getAtomName2() const {return atomName2;};
    const String getChain2() const {return chain2;};
    const double getDistance() const {return distance;}
    void  print() const {
        MMBLOG_FILE_LINE(INFO,   // to here is fine
          " : Chain ID 1 : "      <<getChain1()
          <<" Residue    ID 1 : "    <<getResidueID1().outString()
          <<" atom name 1 : "       <<getAtomName1()
          <<" : Chain ID 2 : "     <<getChain2()
          <<" Residue ID 2: "   <<getResidueID2().outString()
          <<" atom name 2 : "      <<getAtomName2()
          <<" distance : " <<getDistance()
          <<endl);
    };
    };

class CovalentBondClass:public  TwoAtomClass       {
    private:
    SimTK::BondMobility::Mobility mobility; // = SimTK::BondMobility::Rigid;
        SimTK::String bondCenterName1;    
        SimTK::String bondCenterName2;    
    public:
    CovalentBondClass(){
        chain1 = ""; residueID1 = ResidueID(); atomName1 = ""; 
        chain2 = ""; residueID2 = ResidueID(); atomName2 = ""; 
        mobility = SimTK::BondMobility::Rigid;};
    CovalentBondClass(String myChain, ResidueID inputResidueID,String myAtomName) {
        residueID1 = (inputResidueID);
        atomName1 = myAtomName;
        chain1 = myChain;
        residueID2 = ResidueID();
        atomName2 = "" ;        
        chain2 = "";
        mobility = SimTK::BondMobility::Rigid;
    };
    void  print() const {
        MMBLOG_FILE_LINE(INFO,   // to here is fine
          " : Chain ID : "      <<getChain1()
          <<" Residue    ID: "    <<getResidueID1().outString()
          <<" atom name : "       <<getAtomName1()
          <<" : Chain ID2 : "     <<getChain2()
          <<" Residue    ID2: "   <<getResidueID2().outString()
          <<" atom name2 : "      <<getAtomName2()
          <<" bond mobility : "        <<mobility
          <<endl);
    };
    void  setBondCenterName1(String myBondCenterName) {bondCenterName1     = myBondCenterName;}
    void  setBondCenterName2(String myBondCenterName) {bondCenterName2     = myBondCenterName;}
    String  getBondCenterName1() {return bondCenterName1;}
    String  getBondCenterName2() {return bondCenterName2;}
    SimTK::BondMobility::Mobility getBondMobility(){ 
    return mobility;}
    };

enum ConstraintType {WeldToAtom, WeldToGround, CoupledCoordinate, Undefined};
class ConstraintClass : public  TwoAtomClass       {
    private:
    ConstraintType constraintType;
        //bool toGround;
    public:
    ConstraintClass();
    ConstraintClass(String myChain, ResidueID inputResidueID,String myAtomName); 
    ConstraintClass(String myChain, ResidueID inputResidueID,String myAtomName,String myChain2, ResidueID inputResidueID2,String myAtomName2, ConstraintType myConstraintType); 
    void setConstraintType (ConstraintType myConstraintType);
    ConstraintType getConstraintType () const; // {return constraintType ;}
//enum ConstraintType {WeldToAtom, WeldToGround, CoupledCoordinate, Undefined};
    String constraintTypeString () const;/* {
    (constraintType == WeldToAtom)? cout<< "WeldToAtom" :
    (constraintType == WeldToGround)? cout<< "WeldToGround" :
    (constraintType == CoupledCoordinate)? cout<< "CoupledCoordinate" :
    (constraintType == Undefined)? cout<< "Undefined" ;
   
    }*/
    //void  setToGround(bool myToGround) {toGround = myToGround;}
    //const bool   getToGround() const {return toGround;};
    void  print() const;/* {
        std::cout<<__FILE__<<":"<<__LINE__   // to here is fine
          <<" : Chain ID : "      <<getChain1()
          <<" Residue    ID: "    <<getResidueID1().outString()
          <<" atom name : "       <<getAtomName1()
          <<" : Chain ID2 : "     <<getChain2()
          <<" Residue    ID2: "   <<getResidueID2().outString()
          <<" atom name2 : "      <<getAtomName2()
          //<<" to Ground: "        <<getToGround()
          <<" constraintType : " << printConstraintType()
          <<endl;
    };*/


    };

class MMB_EXPORT ContactStretch  : public ResidueStretch  {
    public:
        using ResidueStretch::ResidueStretch;

               SimTK::String ContactScheme           ;
               SimTK::String getContactScheme() {return ContactScheme;};
               //SimTK::String Chain          ;
               int    leontisWesthofBondRowIndex;
               bool sameParameters (ContactStretch contactStretch) {
                    if ((getChain().compare(contactStretch.getChain() ) == 0 )  &&
                        (getContactScheme().compare(contactStretch.getContactScheme() ) == 0 ))  
                        return true;
                    else return false;
               };
};

enum SecondaryStructureType {Alpha, ParallelBeta, AntiParallelBeta};
class SecondaryStructureStretch  : public ResidueStretch  {
    private:
               SecondaryStructureType mySecondaryStructureType;
    public:
           void setSecondaryStructureType(const String &inputSecondaryStructureType) {
            if ((inputSecondaryStructureType.compare("Alpha")) == 0) { 
                mySecondaryStructureType = Alpha;}
            else if ((inputSecondaryStructureType.compare("ParallelBeta")) == 0) { 
                mySecondaryStructureType = ParallelBeta;}
            else if ((inputSecondaryStructureType.compare("AntiParallelBeta")) == 0) { 
                mySecondaryStructureType = AntiParallelBeta;}
            else {
                MMBLOG_FILE_LINE(CRITICAL, " Error!  The only permitted secondary structure types are Alpha, ParallelBeta and AntiParallelBeta."<<endl);
            }
        }
    SecondaryStructureType getSecondaryStructureType() const  {return mySecondaryStructureType ;};
};

class  DensityStretch : public ResidueStretch   {
    public:	
	DensityStretch() : ResidueStretch{}{} // Derived class default constructor calls base class default constructor with (obviously) no arguments. Really it should not be necessary to specify this.
        DensityStretch(const ResidueStretch myResidueStretch){
            setChain(myResidueStretch.getChain());
            setStartResidue(myResidueStretch.getStartResidue());	
            setEndResidue(myResidueStretch.getEndResidue());	
        }	
    };

struct ContactWithin {
    String ContactScheme           ;
    String Chain          ;
    ResidueID    Residue        ;
    double Radius;
    bool operator==(const ContactWithin &other) const {
        // Compare the values, and return a bool result.
        return Chain == other.Chain && Residue == other.Residue;
    }; 
};

//struct AtomSpecification {
//  String chainID;
//  ResidueID residueID;

struct AtomSpring {
    private:	
    public:	
       // Really the parameters should be private.
       String  atom1Name             ;
       String  atom2Name             ;
       ResidueID    atom1Residue          ;
       ResidueID    atom2Residue          ;
       String  atom1Chain;
       String  atom2Chain;
       bool    toGround;
       bool    tether;
       bool    groundLocationIsRelativeToAtom1Location;
       Vec3    groundLocation;
       //AtomSpecification AtomLocationInGround; // This is the chain ID, residue number, and name of an atom. Once the system has been realized to the position stage, we will be able to query that atom and put its location in groundLocation.
       double  forceConstant ;
       double  deadLength    ;
       bool    deadLengthIsFractionOfInitialLength; 
       double  deadLengthFraction;
       /*template <MMBLogger::Severity S>
       void print(){
            MMBLOG_FILE_LINE(S, " Printing AtomSpring  contents:"<<std::endl
                <<" atom1Chain    : "<<atom1Chain                <<std::endl  
                <<" atom2Chain    : "<<atom2Chain                <<std::endl  
                <<" atom1Residue  : "<<atom1Residue.outString()  <<std::endl  
                <<" atom2Residue  : "<<atom2Residue.outString()  <<std::endl  
                <<" atom1Name     : "<<atom1Name                 <<std::endl  
                <<" atom2Name     : "<<atom2Name                 <<std::endl  
                <<" toGround      : "<< toGround                 <<std::endl     
                <<" tether        : "<< tether                   <<std::endl
                <<" groundLocationIsRelativeToAtom1Location : "<< groundLocationIsRelativeToAtom1Location <<std::endl
                <<" groundLocation: "<< groundLocation           <<std::endl
                <<""<<std::endl);
       }*/
       // Note that this works only with loggingSeverity DEBUG:
       void print(){//enum MMBLogger::Severity severity = MMBLogger::Severity::INFO){
	       
            MMBLOG_FILE_LINE(DEBUG, " Printing AtomSpring  contents:"<<std::endl
                <<" atom1Chain    : "<<atom1Chain                  
                <<" atom2Chain    : "<<atom2Chain                  
                <<" atom1Residue  : "<<atom1Residue.outString()    
                <<" atom2Residue  : "<<atom2Residue.outString()    
                <<" atom1Name     : "<<atom1Name                   
                <<" atom2Name     : "<<atom2Name                   
                <<" toGround      : "<< toGround                      
                <<" tether        : "<< tether                   
                <<" deadLength    : "<< deadLength               
                <<" deadLengthIsFractionOfInitialLength : "<< deadLengthIsFractionOfInitialLength  
                <<" deadLengthFraction : "<< deadLengthFraction  
                <<" groundLocationIsRelativeToAtom1Location : "<< groundLocationIsRelativeToAtom1Location 
                <<" groundLocation: "<< groundLocation           
                <<" forceConstant : "<< forceConstant            
                <<""<<std::endl);
       } 
       /*
       void printDebug(){
            MMBLOG_FILE_LINE(DEBUG, " Printing AtomSpring  contents:"<<std::endl
                <<" atom1Chain    : "<<atom1Chain                <<std::endl  
                <<" atom2Chain    : "<<atom2Chain                <<std::endl  
                <<" atom1Residue  : "<<atom1Residue.outString()  <<std::endl  
                <<" atom2Residue  : "<<atom2Residue.outString()  <<std::endl  
                <<" atom1Name     : "<<atom1Name                 <<std::endl  
                <<" atom2Name     : "<<atom2Name                 <<std::endl  
                <<" toGround      : "<< toGround                 <<std::endl     
                <<" tether        : "<< tether                   <<std::endl
                <<" deadLength    : "<< deadLength               <<std::endl
                <<" groundLocationIsRelativeToAtom1Location : "<< groundLocationIsRelativeToAtom1Location <<std::endl
                <<" groundLocation: "<< groundLocation           <<std::endl
                <<""<<std::endl);
       }*/


       void initialize(){
	    MMBLOG_FILE_FUNC_LINE(DEBUG, " initializing atomSpring  "<<endl);
            atom1Name = "XXXX"           ;   
            atom2Name = "XXXX"           ;   
            atom1Residue =  ResidueID(-1111, ' ')        ;   
            atom2Residue =  ResidueID(-1111, ' ')        ;   
            atom1Chain = "XXXX";
            atom2Chain = "XXXX";
            toGround   = false;
            tether     = false;
	    groundLocationIsRelativeToAtom1Location = false;
            groundLocation = Vec3(0);
            forceConstant  = 0.0 ;
            deadLength     = 0.0 ; 
            deadLengthFraction     = 0.0 ; 
	    forceConstant  = 129790.8 ; // base the default value on the spring constant of a carbon-carbon single bond
	    print(); // Only acts if loggingSeverity DEBUG
       }
       AtomSpring(String chain1, ResidueID res1, String name1,
                  String chain2, ResidueID res2, String name2,
                  double constant,
                  double deadLength = 0.0,
                  double deadLengthFraction = 0.0,
		  bool deadLengthIsFractionOfInitialLength = 0,
                  Vec3 location = Vec3(0),
                  bool toGround = false,
                  bool tether = false,
		  bool   groundLocationIsRelativeToAtom1Location = 0
                 )
       {
            atom1Chain = chain1;
            atom2Chain = chain2;
            atom1Name = name1;
            atom2Name = name2;
            atom1Residue = res1;
            atom2Residue = res2;
            forceConstant = constant;
            this->toGround = toGround;
            this->tether = tether;
            this->groundLocationIsRelativeToAtom1Location = groundLocationIsRelativeToAtom1Location;
            this->groundLocation = location;
            this->deadLength = deadLength;
            this->deadLengthFraction = deadLengthFraction;
	    this->deadLengthIsFractionOfInitialLength = deadLengthIsFractionOfInitialLength;
       }

       AtomSpring()
       {
	    initialize();
            /*atom1Name = "XXXX"           ;   
            atom2Name = "XXXX"           ;   
            atom1Residue =  ResidueID(-1111, ' ')        ;   
            atom2Residue =  ResidueID(-1111, ' ')        ;   
            atom1Chain = "XXXX";
            atom2Chain = "XXXX";
            toGround   = false;
            tether     = false;
	    groundLocationIsRelativeToAtom1Location = false;
            groundLocation = Vec3(0);
            forceConstant  = 0.0 ;
            deadLength     = 0.0 ; 
            deadLengthFraction     = 0.0 ; */
       }
};


struct IncludeNonBondAtomInBiopolymerStruct {
    String chain;
    ResidueID residue;
    String atomName;        
};

struct WaterDropletAboutResidueStruct {
    String biopolymerChainID;
    ResidueID residue;
    String waterDropletChainID;
    double radius;
    double tetherStrength;
};

class MMBAtomInfo {
    private:
        String chain;
        ResidueID residueID;
        ResidueInfo::Index  residueIndex;
        String atomName;
    public:
        MobilizedBody mobilizedBody;
        MobilizedBodyIndex mobilizedBodyIndex;
        Compound::AtomIndex compoundAtomIndex;
        double  mass;
        int  atomicNumber;
        openmmVecType position;
        double partialCharge;
        std::vector<MMBAtomInfo*> neighbors;
        //ChargedAtomType chargedAtomType;
        void setAtomName(const String myAtomName ){ atomName = myAtomName;}
        const String & getAtomName() const { return atomName; }
        void setResidueID(const ResidueID myResidueID ){ residueID = myResidueID;}
        void setChain(const String myChain ){ chain = myChain;}
        String getChain( ){ return chain ;}
        ResidueInfo::Index getResidueIndex( ){ return  residueIndex;}
        ResidueID          getResidueID   ( ){ return  residueID   ;}
        void setResidueIndex(ResidueInfo::Index  myResidueIndex ){ residueIndex = myResidueIndex;}
        MMBAtomInfo(){};

        MMBAtomInfo(String myChain, ResidueID myResidueID, String myAtomName) :
            chain(std::move(myChain)),
            residueID(std::move(myResidueID)),
            atomName(std::move(myAtomName))
        {}

        MMBAtomInfo(String myChain, ResidueID myResidueID, ResidueInfo::Index myResidueIndex, String myAtomName) :
            chain(std::move(myChain)),
            residueID(std::move(myResidueID)),
            residueIndex(std::move(myResidueIndex)),
            atomName(std::move(myAtomName))
        {}

        bool operator == (MMBAtomInfo & a){
        if (this->compoundAtomIndex == a.compoundAtomIndex ) {return true;}
        else return false;
        }
        void print () {
            MMBLOG_FILE_LINE(INFO, " Printing MMBAtomInfo contents:"<<std::endl
                <<" mobilizedBodyIndex : "<<mobilizedBodyIndex
                <<" compoundAtomIndex : "<< compoundAtomIndex
                <<" position : "<< position
                <<" residueID: "<< residueID.outString()
                <<" chain: "<<chain
                <<" atomName : "<<atomName
                <<" partial charge: "<<partialCharge<<std::endl);
            //std::cout<<__FILE__<<":"<<__LINE__<< mobilizedBodyIndex<<","<< compoundAtomIndex<<","<<  position<<","<< residueID.outString()<<","<< chain <<", "<< atomName<<", "<<partialCharge<<std::endl;
            
        }

        double distance(const MMBAtomInfo & atom2) const
        {
            double x = (atom2.position[0] - position[0]);
            double y = (atom2.position[1] - position[1]);
            double z = (atom2.position[2] - position[2]);
            return sqrt( (x*x) + (y*y) +(z*z));
        }
        void clearNeighbors()
        {
            neighbors.clear();
        }
        void addNeighbor(MMBAtomInfo* neighbor)
        {
            neighbors.push_back(neighbor);
        }
};




/*class MMB_EXPORT BiopolymerClassContainer{
    public:
	vector<MMBAtomInfo> getConcatenatedAtomInfoVector(const State & state,bool activeChainsOnly=false);
	vector<MMBAtomInfo> getConcatenatedAtomInfoVector(bool activeChainsOnly=false);
};*/

class InterfaceContainer    { 
    private:
        vector <Interface> interfaceVector;
     public:
	void clear() {interfaceVector.clear() ;};
	InterfaceContainer() {clear() ;};
        //Interface getInterface(int interfaceIndex) {return interfaceVector[interfaceIndex];};
        vector<String> getChains(int interfaceIndex) {return getInterface(interfaceIndex).Chains;};
        vector<String> getPartnerChains(int interfaceIndex) {return getInterface(interfaceIndex).PartnerChains;};
        double getDepth(int interfaceIndex) {return interfaceVector[interfaceIndex].Depth;};
        String getMobilizerString(int interfaceIndex) {return interfaceVector[interfaceIndex].MobilizerString;};
        void addInterface(String myChain, double myDepth ,  String myMobilizerString = "NONE")
            {Interface myInterface; myInterface.Chains.push_back( myChain);  myInterface.Depth = myDepth; myInterface.MobilizerString = myMobilizerString; interfaceVector.push_back(myInterface); };
        void addInterface(vector<String> myChains, double myDepth = 0.0 ,  String myMobilizerString = "NONE"){Interface myInterface; 
                        myInterface.Chains.clear(); myInterface.PartnerChains.clear();
            for (unsigned int i = 0; i < myChains.size(); i++) {myInterface.Chains.push_back( myChains[i]);}  myInterface.Depth = myDepth; myInterface.MobilizerString = myMobilizerString; interfaceVector.push_back(myInterface); };
        void addInterface(vector<String> myChains,vector<String> partnerChains,  double myDepth ,  String myMobilizerString = "NONE");
        vector<TwoAtomClass> retrieveCloseContactPairs(vector<MMBAtomInfo> & concatenatedAtomInfoVector );
        Interface getInterface(int interfaceIndex) {return  interfaceVector[interfaceIndex];};
        unsigned int numInterfaces() {return interfaceVector.size();};
        void print(){
            for (unsigned int i = 0 ; i < numInterfaces(); i++){
                 getInterface(i).print();
            } 
        };
};

namespace BiopolymerType {
    enum BiopolymerTypeEnum {
        RNA,
        Protein,
        DNA,
        Unassigned
    };
}

//double ValidateDouble(double myDouble);

Vec3 ValidateVec3(const Vec3 myVec3);
int ValidateNonNegativeInt (const int myInt);
int ValidateInt (const int myInt);
double ValidateNonNegativeDouble(const double myDouble) ;
double DoubleFromStringNoSymbols(const String inString);
double ValidateDouble(const double myDouble) ;
Real DotProduct(const Vec3 vecA, const Vec3 vecB);
vector<String> readAndParseLine   (ifstream & inFile );
String removeAllWhite (String &str);
vector<String> readAndParseOnColWidth   (ifstream & inFile, int columnWidth);


#endif
