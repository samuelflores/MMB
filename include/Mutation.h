// Mutation.h
#ifndef Mutation_H_
#define Mutation_H_
#include "Utils.h"
#include <string>
#include <iostream>
#define FOLDXSEPARATOR ","
using namespace SimTK;
class Mutation : public SingleResidue {
private:
        //std::string chain;
        // class SingleResidue uses the startResidue variable to hold its one and only ResidueID   
        //ResidueID startResidue;//residueID;
        std::string substitutedResidueType;
        std::string wildTypeResidueType;
public:
        Mutation() {setChain ( " "); setResidue (ResidueID( -1111, ' ')); 
            //startResidue.setInsertionCode(' '); 
            setWildTypeResidueType ( "RESIDUE-TYPE-NOT-SET" );
            substitutedResidueType = "RESIDUE-TYPE-NOT-SET"; 
        }
        // Converts the current Mutation to an AllResiduesWithin, adding a value for the radius member
        AllResiduesWithin allResiduesWithin(double radius){AllResiduesWithin myAllResiduesWithin((*this).getChain(), (*this).getResidue(), radius); return myAllResiduesWithin; };
        void print(){
            std::cout<<__FILE__<<":"<<__LINE__<<" Mutation chain = >"<<getChain()<<"<, residue = >"<<  getResidue().outString()<<"<, substitution = >"<<substitutedResidueType<<"<"<<std::endl; }
        void setSubstitutedResidueType(std::string myResidueType) {substitutedResidueType = myResidueType;};
        std::string getSubstitutedResidueType() {
                std::cout<<__FILE__<<":"<<__LINE__<<" About to return substitutedResidueType = "<<std::endl;
                std::cout<<__FILE__<<":"<<__LINE__<<"  >"<<substitutedResidueType<<"<"<<std::endl;
                return substitutedResidueType;
        };
        int getResidueNumber() {return   getResidue().getResidueNumber() ;}
        const std::string getInsertionCode() const { std::string myInsertionCode; myInsertionCode = (  getResidue().getInsertionCode()) ; return myInsertionCode;}
        std::string getResidueIDAsString() {return   getResidue().outString();}
        std::string getMutationAsString() {return (getChain() +  getResidueIDAsString() + getSubstitutedResidueType());}
        std::string getMutationAsFoldxString() {return  getWildTypeResidueType() + getChain() +  getResidueIDAsString() + getSubstitutedResidueType();}
        void validate() {
                if (getResidueNumber() < 0) {
                        std::cout<<__FILE__<<":"<<__LINE__<<" WARNING! The residueNumber is < 0!"<<std::endl;
                        // Downgraded from error to warning, 24 nov 2019. seems to work fine.
                        //exit(1);
                }
                if ((substitutedResidueType.length() == 0) || (substitutedResidueType.length() > 1)) {
                        std::cout<<__FILE__<<":"<<__LINE__<<" The substitutedResidueType.length() is either 0 or >1."<<std::endl;
                        exit(1);
                }
                if ((getInsertionCode().length() == 0) || (getInsertionCode().length() > 1)) {
                        std::cout<<__FILE__<<":"<<__LINE__<<" getInsertionCode().length() returned either 0 or >1."<<std::endl;
                        exit(1);
                }
                std::cout<<__FILE__<<":"<<__LINE__<<" Validated: "<<getMutationAsString()<<std::endl;
        }

        Mutation(String myChainID, ResidueID myResidueID, String mySubstitutedResidueType ) {
		setChain (myChainID);  
		setResidue(myResidueID); 
		setSubstitutedResidueType(mySubstitutedResidueType);
		validate();}

       // This function expects either a breeder formatted (C-NNNI-S) or SKEMPI formatted (WCNNIS) single mutation string. Based on the presence or absence of MUTATIONMINORSEPARATOR (currently '-'), it determines which is being used.
        std::string getWildTypeResidueType( ){ 
            std::cout<<__FILE__<<":"<<__LINE__<<" Returning wildTypeResidueType = >"<<wildTypeResidueType<<"< "<<std::endl;
            return wildTypeResidueType;}
        void setWildTypeResidueType(std::string myWildTypeResidueType ){
            wildTypeResidueType = myWildTypeResidueType;
        }
        void setChainSubstitutionFromSingleMutationString(std::string mySingleMutationString){
            std::cout<<__FILE__<<":"<<__LINE__<<" processing : >"<<mySingleMutationString<<"< "<<std::endl; 
            size_t minorSeparatorPosition1 = mySingleMutationString.find(MUTATIONMINORSEPARATOR);
            std::cout<<__FILE__<<":"<<__LINE__<<std::endl; 
            if (minorSeparatorPosition1 == std::string::npos){ // This is the case the mutation string is in SKEMPI format.
                std::cout<<__FILE__<<":"<<__LINE__<<std::endl; 
                std::string myChain = mySingleMutationString.substr(1,1); // The second character must be the chain ID
                std::cout<<__FILE__<<":"<<__LINE__<<std::endl; 
                setChain (myChain);
                std::cout<<__FILE__<<":"<<__LINE__<<std::endl; 
                std::string myResidueIDString = mySingleMutationString.substr(2,mySingleMutationString.length()-3); // From the third to the penultimate character, must be the residue number and insertion string.
                std::cout<<__FILE__<<":"<<__LINE__<<std::endl; 
                setResidue(ResidueID(myResidueIDString));
                std::cout<<__FILE__<<":"<<__LINE__<<" mySingleMutationString.length() = "<<mySingleMutationString.length()<< std::endl; 
                std::string mySubstitution = mySingleMutationString.substr(mySingleMutationString.length()-1,1); // The last character, must be the substituted residue type.
                std::cout<<__FILE__<<":"<<__LINE__<<std::endl; 
                setSubstitutedResidueType(mySubstitution);
                std::cout<<__FILE__<<":"<<__LINE__<<std::endl; 
            } else { // This is the case the mutation string is in breeder format.
                std::cout<<__FILE__<<":"<<__LINE__<<std::endl; 
                if ((minorSeparatorPosition1 == 0) || (minorSeparatorPosition1 == (mySingleMutationString.length()-1))){
                    std::cout<<__FILE__<<":"<<__LINE__<<" Invalid minor separator position : "<<minorSeparatorPosition1<<std::endl; exit(1);
                }
                size_t minorSeparatorPosition2 = mySingleMutationString.find(MUTATIONMINORSEPARATOR, (minorSeparatorPosition1 + 1));
                std::cout<<__FILE__<<":"<<__LINE__<<std::endl; 
                setChain(mySingleMutationString.substr(0,minorSeparatorPosition1));
                std::cout<<__FILE__<<":"<<__LINE__<<std::endl; 
                std::string residueIDString = mySingleMutationString.substr(minorSeparatorPosition1+1,( minorSeparatorPosition2 - minorSeparatorPosition1 - 1));
                std::cout<<__FILE__<<":"<<__LINE__<<std::endl; 
                setResidue(ResidueID(residueIDString));
                std::cout<<__FILE__<<":"<<__LINE__<<"..."<<std::endl; 
                
                std::string mySubstitution = mySingleMutationString.substr(minorSeparatorPosition2+1,mySingleMutationString.length()-1 ); // The last character, must be the substituted residue type.
                if (mySubstitution.length() != 1){
                    std::cout<<__FILE__<<":"<<__LINE__<<" Invalid substituted residue type : >"<< mySubstitution <<"< " <<std::endl; exit(1);
                }
                std::cout<<__FILE__<<":"<<__LINE__<<std::endl; 
                setSubstitutedResidueType(mySubstitution);
            }
        } // of setChainSubstitutionFromSingleMutationString 
};
#endif
