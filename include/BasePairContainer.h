/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#ifndef BasePairContainer_H_
#define BasePairContainer_H_

#include "MMB_config.h"

#include "SimTKmolmodel.h"
#include "BiopolymerClass.h"
#include "Utils.h"          
#include "ResidueStretchContainer.h"
#include "BaseInteractionParameterReader.h"
#include "NTC_PARAMETER_READER.h"

class MMB_EXPORT BasePairContainer {

public:
	void		clear() ;
	    		BasePairContainer(){clear();} ;
	void		addBasePair(BiopolymerClassContainer & myBiopolymerClassContainer, 
                            const LeontisWesthofClass & lhClass, 
                            BaseInteraction myBasePair, bool helicalStacking=false);
	void		deleteBasePair(int basePairIndex );
    void        updateBasePair(int index, 
                               String ch1, int res1, String edge1, 
                               String ch2, int res2, String edge2, 
                               String orient,
                               BiopolymerClassContainer& myBiopolymerClassContainer,
                               const LeontisWesthofClass& lhClass, bool helicalStacking=false);

	void		validateBasePair(BiopolymerClassContainer & myBiopolymerClassContainer, 
                                 const LeontisWesthofClass & lhClass, 
                                 BaseInteraction & myBasePair, bool helicalStacking=false);
	const BaseInteraction & getBasePair(int basePairIndex) ;
	//BaseInteraction & 	updBasePair(int);
	int 		numBasePairs() ;
#ifdef MMB_NTC_ENABLED
	void		addHelicalStacking(BiopolymerClassContainer & myBiopolymerClassContainer, const LeontisWesthofClass & lhClass,const  NTC_PAR_Class & ntc_par_class, NTC_Class_Container  & my_ntc_class_container);
#else
	void		addHelicalStacking(BiopolymerClassContainer & myBiopolymerClassContainer, const LeontisWesthofClass & lhClass);
#endif // MMB_NTC_ENABLED
    vector<BaseInteraction>	myBasePairVector;	    
    void 		printBasePairs();
    void		setBasePairSatisfied(int,bool);

    vector<int> getSatisfiedBasePairs();

    String getBasePairsStrings();


private:
	
	bool 		hasWatsonCrickCisPair(String chainID ,ResidueID residueNumber) ; 
	const BaseInteraction 	getWatsonCrickCisPair(String chainID ,ResidueID residueNumber) ; 
	const String 	getWatsonCrickCisPairingChain(String chainID ,ResidueID residueNumber) ; 
	const ResidueID 	getWatsonCrickCisPairingResidue(String chainID ,ResidueID residueNumber) ; 
	const ResidueID		getLastWatsonCrickCisPairingResidueOfRun(String chainID ,ResidueID firstResidueNumberInStack, BiopolymerClassContainer & myBiopolymerClassContainer) ; 
	void		generateHelicalStackingInteractions(String chainID,ResidueID  firstResidue, ResidueID lastResidue,BiopolymerClassContainer & myBiopolymerClassContainer, const LeontisWesthofClass & lhClass);
}; // of class

#endif
