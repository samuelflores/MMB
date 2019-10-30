/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#include "SimTKmolmodel.h"
#include "BiopolymerClass.h"
#include "Utils.h"          
#include "ResidueStretchContainer.h"
#include "NtC_Class_Container.h"
#include "NTC_PARAMETER_READER.h"

#include <sstream>

void NTC_Class_Container::clear(){
    myNTC_Class_Vector.clear();
};

void NTC_Class_Container::add_NTC_Class(BiopolymerClassContainer & myBiopolymerClassContainer, 
                                    const NTC_PAR_Class & ntc_par_class, NTC_Classes & NTC, bool NTC_add){
    
    
    String dihedraltype;
    
    dihedraltype = "delta";
    
    validate_NTC_Class(myBiopolymerClassContainer, ntc_par_class,NTC, dihedraltype);
    myNTC_Class_Vector.push_back(NTC);
    
    dihedraltype = "epsilon";

    validate_NTC_Class(myBiopolymerClassContainer, ntc_par_class,NTC, dihedraltype);
    myNTC_Class_Vector.push_back(NTC);

    dihedraltype = "zeta";

    validate_NTC_Class(myBiopolymerClassContainer, ntc_par_class,NTC, dihedraltype);
    myNTC_Class_Vector.push_back(NTC);    
    
    dihedraltype = "alpha1";

    validate_NTC_Class(myBiopolymerClassContainer, ntc_par_class,NTC, dihedraltype);
    myNTC_Class_Vector.push_back(NTC);   

    dihedraltype = "beta1";

    validate_NTC_Class(myBiopolymerClassContainer, ntc_par_class,NTC, dihedraltype);
    myNTC_Class_Vector.push_back(NTC);       
    
    dihedraltype = "gamma1";

    validate_NTC_Class(myBiopolymerClassContainer, ntc_par_class,NTC, dihedraltype);
    myNTC_Class_Vector.push_back(NTC);     
    
    dihedraltype = "delta1";

    validate_NTC_Class(myBiopolymerClassContainer, ntc_par_class,NTC, dihedraltype);
    myNTC_Class_Vector.push_back(NTC);     
    
    dihedraltype = "chi";

    validate_NTC_Class(myBiopolymerClassContainer, ntc_par_class,NTC, dihedraltype);
    myNTC_Class_Vector.push_back(NTC);        
    
    dihedraltype = "chi1";

    validate_NTC_Class(myBiopolymerClassContainer, ntc_par_class,NTC, dihedraltype);
    myNTC_Class_Vector.push_back(NTC);     

    dihedraltype = "nccn";

    validate_NTC_Class(myBiopolymerClassContainer, ntc_par_class,NTC, dihedraltype);
    myNTC_Class_Vector.push_back(NTC);     

    dihedraltype = "nn";

    validate_NTC_Class(myBiopolymerClassContainer, ntc_par_class,NTC, dihedraltype);
    myNTC_Class_Vector.push_back(NTC);     

    dihedraltype = "cc";

    validate_NTC_Class(myBiopolymerClassContainer, ntc_par_class,NTC, dihedraltype);
    myNTC_Class_Vector.push_back(NTC);  

    dihedraltype = "tau0a";

    validate_NTC_Class(myBiopolymerClassContainer, ntc_par_class,NTC, dihedraltype);
    myNTC_Class_Vector.push_back(NTC);     

    dihedraltype = "tau1a";

    validate_NTC_Class(myBiopolymerClassContainer, ntc_par_class,NTC, dihedraltype);
    myNTC_Class_Vector.push_back(NTC);     

    dihedraltype = "tau2a";

    validate_NTC_Class(myBiopolymerClassContainer, ntc_par_class,NTC, dihedraltype);
    myNTC_Class_Vector.push_back(NTC);     

    dihedraltype = "tau3a";

    validate_NTC_Class(myBiopolymerClassContainer, ntc_par_class,NTC, dihedraltype);
    myNTC_Class_Vector.push_back(NTC);     

    dihedraltype = "tau4a";

    validate_NTC_Class(myBiopolymerClassContainer, ntc_par_class,NTC, dihedraltype);
    myNTC_Class_Vector.push_back(NTC);     

    dihedraltype = "tau0b";

    validate_NTC_Class(myBiopolymerClassContainer, ntc_par_class,NTC, dihedraltype);
    myNTC_Class_Vector.push_back(NTC);     

    dihedraltype = "tau1b";

    validate_NTC_Class(myBiopolymerClassContainer, ntc_par_class,NTC, dihedraltype);
    myNTC_Class_Vector.push_back(NTC);     

    dihedraltype = "tau2b";

    validate_NTC_Class(myBiopolymerClassContainer, ntc_par_class,NTC, dihedraltype);
    myNTC_Class_Vector.push_back(NTC);     

    dihedraltype = "tau3b";

    validate_NTC_Class(myBiopolymerClassContainer, ntc_par_class,NTC, dihedraltype);
    myNTC_Class_Vector.push_back(NTC);     

    dihedraltype = "tau4b";

    validate_NTC_Class(myBiopolymerClassContainer, ntc_par_class,NTC, dihedraltype);
    myNTC_Class_Vector.push_back(NTC);     

}

void NTC_Class_Container::validate_NTC_Class(BiopolymerClassContainer & myBiopolymerClassContainer, 
                                         const NTC_PAR_Class & ntc_par_class, NTC_Classes & NTC,
                                         String dihedraltype){
    
  //  NTC_PAR_BondRow ntc_par_bondrow;
    
   //cout << NTC.NtC_FirstBPChain << " " << NTC.NtC_step_ID << " " << NTC.NtC_Class_String << endl; 
    
    String ntc2 = to_string((stoi(NTC.NtC_step_ID))+1);
    
    String resName1 = myBiopolymerClassContainer.getPdbResidueName(NTC.NtC_FirstBPChain, NTC.NtC_step_ID);
    String resName2 = myBiopolymerClassContainer.getPdbResidueName(NTC.NtC_FirstBPChain, ntc2);
    
    NTC_PAR_BondRow br = ntc_par_class.getNTC_PAR_BondRow(NTC.NtC_step_ID,ntc2,resName1,NTC.NtC_Class_String,resName2,NTC.NtC_Class_String,dihedraltype,"ntcstep");
 
    NTC.NTC_PAR_BondRowIndex  = ntc_par_class.getNTC_PAR_BondRowIndex(resName1,resName2,NTC.NtC_Class_String,dihedraltype,"ntcstep",NTC ); 
    
    cout << NTC.NTC_PAR_BondRowIndex << " BOND ROW ";
   
 //   NTC.NtC_INDEX = ntc_par_class.getNTC_PAR_BondRow
    
}

int NTC_Class_Container::numNTC_Torsions() {
    return myNTC_Class_Vector.size();
}

const NTC_Classes & NTC_Class_Container::getNTC_Class(int NTC_Class_Index) {
    return myNTC_Class_Vector[NTC_Class_Index];
}


