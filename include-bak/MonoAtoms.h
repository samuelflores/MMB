/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#ifndef MonoAtoms_H_
#define MonoAtoms_H_

#include <string>
#include "SimTKmolmodel.h"
//#include "DuMMForceFieldSubsystem.h"


using namespace std;
using namespace SimTK;

class MonoAtoms {
	private:
		String 	chainID;
    		ResidueID firstResidueID    ;
		int 	numAtoms;
		String 	atomName;
		vector<Compound> compoundVector; 
	public:
		MonoAtoms();
		MonoAtoms(String, ResidueID, int, String);
                void initialize(CompoundSystem & system,  bool readPreviousFrameFile,  String previousFrameFileName, bool matchExact, bool matchIdealized);
		int 		validate();  //return 0 if fine, return 1 and exit if not fine.
		void		validateResidue(ResidueID) ; //die if residue number out of range, or if atom names mismatch.
		const String	getChainID();
		const ResidueID	getLastResidueID ();
		const ResidueID	getFirstResidueID ();
                const int       getResidueIndex(ResidueID myResidueID);
		const int 	getNumAtoms();
		const String	getAtomName();  
		Compound 	getSingleCompound(ResidueID); 
		MobilizedBodyIndex getMobilizedBodyIndex(ResidueID);
		MobilizedBody 	updMobilizedBody(ResidueID, SimbodyMatterSubsystem &) ;
		const bool      hasAtom(ResidueID residueNumber);
		const Compound::AtomIndex getAtomIndex(ResidueID residueNumber);
		const Vec3	getAtomLocationInMobilizedBodyFrame(ResidueID residueNumber);
		void		adoptCompounds(CompoundSystem & mySystem, bool readPreviousFrameFile);
		void		matchDefaultConfiguration(SimTK::PdbStructure,bool matchExact, bool matchIdealized);
                void            renumberPdbResidues(ResidueID firstResidueID);
                void            setPdbChainId(String chainID);
		String		getAtomPathName(ResidueID);
		void 		includeAllAtoms(DuMMForceFieldSubsystem & dumm);
};
	
class MonoAtomsContainer  {
	private:
		map <String,MonoAtoms> monoAtomsMap;
	public:
		MonoAtomsContainer();
                void initialize(CompoundSystem &  , bool , String ,bool matchExact, bool matchIdealized  );
		//const int       getChainID(String);
		const bool	hasChainID(String) ; 
		MonoAtoms 	getMonoAtoms(String) ; 
		void 		addMonoAtoms(MonoAtoms) ; 
		void		adoptCompounds(CompoundSystem & mySystem, bool readPreviousFrameFile);
		void 		matchDefaultConfiguration(SimTK::PdbStructure ,bool matchExact, bool matchIdealized );
		String		getAtomPathName(String,ResidueID);
		void 		clear();
		void 		includeAllAtoms(DuMMForceFieldSubsystem & dumm);
};
#endif

