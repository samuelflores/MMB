// Mutation.h
#ifndef Mutation_H_
#define Mutation_H_
#include "Utils.h"
#include <string>
#include <iostream>

class Mutation {
private:
        std::string chainID;
        ResidueID residueID;
        std::string substitutedResidueType;
public:
        Mutation() {chainID = " "; residueID.setResidueNumber ( -1111); residueID.setInsertionCode(' '); substitutedResidueType = "RESIDUE-TYPE-NOT-SET";}
        void setResidueNumber(int myResidueNumber) {residueID.setResidueNumber  (myResidueNumber);}
        void setResidueID    (ResidueID myResidueID) {residueID = myResidueID;}
        void setInsertionCode(std::string myInsertionCode) {
                if (myInsertionCode.length() > 1) {
                        std::cout<<__FILE__<<":"<<__LINE__<<" The insertion code you provided is too long! Max is 1 char."<<std::endl; exit(1);
                }
                std::cout<<__FILE__<<":"<<__LINE__<<" About to set InsertionCode to first char of = >"<<myInsertionCode<<"< "<<std::endl;
                residueID.setInsertionCode(myInsertionCode[0]);
	}
        void setSubstitutedResidueType(std::string myResidueType) {substitutedResidueType = myResidueType;};
        void setChainID(std::string myChainID ){
                if (myChainID.length() >1) {
                        std::cout<<__FILE__<<":"<<__LINE__<<" Invalid chain ID: "<<myChainID<<std::endl; exit(1);
                }
                chainID = myChainID;};
        std::string getChainID(){return chainID;};
        std::string getSubstitutedResidueType() {
                std::cout<<__FILE__<<":"<<__LINE__<<" About to return substitutedResidueType = "<<std::endl;
                std::cout<<__FILE__<<":"<<__LINE__<<"  >"<<substitutedResidueType<<"<"<<std::endl;
                return substitutedResidueType;
        };
        ResidueID getResidueID() {return residueID ;}
        int getResidueNumber() {return residueID.getResidueNumber() ;}
        std::string getInsertionCode() { std::string myInsertionCode; myInsertionCode = (residueID.getInsertionCode()) ; return myInsertionCode;}
        std::string getResidueIDAsString() {return residueID.outString();}
        std::string getMutationAsString() {return (getChainID() +  getResidueIDAsString() + getSubstitutedResidueType());}
        void validate() {
                if (getResidueNumber() < 0) {
                        std::cout<<__FILE__<<":"<<__LINE__<<" The residueNumber is < 0!"<<std::endl;
                        exit(1);
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
		setChainID (myChainID);  
		setResidueID(myResidueID); 
		setSubstitutedResidueType(mySubstitutedResidueType);
		validate();}
};
#endif
