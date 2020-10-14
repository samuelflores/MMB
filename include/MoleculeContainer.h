
#ifndef MoleculeContainer_H_
#define MoleculeContainer_H_

#include <vector>
#include "SimTKmolmodel.h"
#include "ExportMacros.h"
#include "ConstraintContainer.h"
#include "Utils.h"
using namespace SimTK;
using namespace std;  

//struct CompoundObject {
//    String   compoundName;
//    Compound compound;
//}

//vector <String> moleculeBuildCommand;

class CompoundObjectMapContainer {
public:
/**
 * \brief This map should ideally contain all the descendants of Compound.  We will use it to translate strings to corresponding Compound objects. Note that this does NOT contain actual molecules for a given simulation, it is a static catalog of Molmodel Compound objects for internal MMB use.
 */
    void     clear() {compoundMap.clear(); /*singleAtomMap.clear();*/}
    void     loadCompoundMap () ;
    void     printCompoundMap() ;

    //void     loadSingleAtomMap () ;
    void     loadBiotypeMap () ;
    CompoundObjectMapContainer () {clear(); loadCompoundMap(); /*loadSingleAtomMap ()*/ ; loadBiotypeMap(); }
    //Compound::SingleAtom fetchSingleAtom(const String);
    Compound::SingleAtom fetchSingleAtom(const String className,  Compound::AtomName& atomName, String elementName , Angle );
    const Element fetchElement (String elementName);
    Compound fetchCompound(const String);
    Biotype fetchBiotype(const String);
private:
    map <const String , Compound> compoundMap;
    //map <const String , Compound::SingleAtom> singleAtomMap;
    map <const String , Biotype> biotypeMap;
};



/**
 * \brief Make a custom molecule
 */
class  CustomMolecule : public Molecule {
public:
    CustomMolecule(){};
    CustomMolecule(const vector <vector <String> > moleculeBuildCommandVector, DuMMForceFieldSubsystem & dumm) ;
};

class MMB_EXPORT MoleculeClass{
public:
    MoleculeClass () {moleculeBuildCommandVector.clear();}
    MoleculeClass (const vector < vector <String> > myMoleculeBuildCommandVector ) {moleculeBuildCommandVector.clear(); moleculeBuildCommandVector=myMoleculeBuildCommandVector; } // Doesn't actually build the molecule.
    void setChainID(String myChainID) {chainID = myChainID;};
    String getChainID() {return chainID;};
    void setResidueName(String myResidueName) {residueName = myResidueName;};
    String getResidueName() {return residueName;};
    void setPdbChainID() {molecule.setPdbChainId(chainID);};
    void setPdbResidueName() ;
    void addOneCommand(vector <String> command) {moleculeBuildCommandVector.push_back(command);  
    }
    
    void includeAllAtoms(DuMMForceFieldSubsystem & dumm);
    void addRingClosingBond(CovalentBondClass myBond);
    vector < vector <String> > moleculeBuildCommandVector;
    CustomMolecule molecule;
private:
    String chainID;
    String residueName;
};


class MMB_EXPORT MoleculeClassContainer {
public :
    map <const String, MoleculeClass> moleculeClassMap;
    void clear() {moleculeClassMap.clear();}
    MoleculeClassContainer() {clear();}
    MoleculeClass &   updMoleculeClass(String myChainID);
    void validateChainID(String myChainID);
    void add(String chainID, String myResidueName, MoleculeClass & myMoleculeClass);
    void initializeCompounds(DuMMForceFieldSubsystem & dumm);
    void adoptCompounds(SimTK::CompoundSystem & mySystem);
    void matchDefaultConfiguration(bool readPreviousFrameFile, String pdbFileName,bool matchExact, bool matchIdealized);
    void includeAllAtoms(DuMMForceFieldSubsystem & dumm);
    bool hasChainID(String chain);
    //void MoleculeClassContainer::matchDefaultConfiguration(bool readPreviousFrameFile, String pdbFileName,bool matchExact, bool matchIdealized)
    void addConstraintToGround(map<const String,double> myUserVariables,  const String chain, const String atomName, ConstraintToGroundContainer & constraintToGroundContainer);

};

#endif
