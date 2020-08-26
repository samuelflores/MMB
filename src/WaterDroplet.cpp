/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#include "WaterDroplet.h"
#include "BiopolymerClass.h"
//#include "Utils.h"
using namespace std;

WaterDroplet::WaterDroplet() {
    //chainID = "D";
    firstResidueNumber = ResidueID(1,' ');
    tetherType = "ToGround";
    tetherStrength = 3.0;
    defaultTransformVector.clear();
}


double WaterDroplet::getRadius(){
    return radius;
}
double WaterDroplet::setRadius(double myRadius ){
    radius = myRadius;
    return radius;
}


void WaterDroplet::printPDB(){
    ErrorManager::instance <<__FILE__<<":"<<__LINE__<<"This doesn't work yet!"<<endl; ErrorManager::instance.treatError();
};

int WaterDroplet::validate(){
    float minRadius = 0;
    if (getRadius() < minRadius){
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" The droplet radius must be greater than "<<minRadius<<".  You have specified: "<<getRadius()<<endl;
        ErrorManager::instance.treatError();}
    if (std::isnan(center[0])) {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" x-component "<<center[0] <<" of Vec3 "<<center <<" is invalid. "<<endl; 
        ErrorManager::instance.treatError();
    }
    else if (std::isnan(center[1])) {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" y-component "<<center[1] <<" of Vec3 "<<center <<" is invalid. "<<endl; 
        ErrorManager::instance.treatError();
    }
    else if (std::isnan(center[2])) {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" z-component "<<center[2] <<" of Vec3 "<<center <<" is invalid. "<<endl; 
        ErrorManager::instance.treatError();
    }

    /*if (chainID.length() > 1){
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" The droplet chain ID must be a single character in length. You have specified: "<<chainID<<endl;
        ErrorManager::instance.treatError();
    }*/
    validateWaterVector();
    return 0;
};

void WaterDropletContainer::clear   (){
    waterDropletVector.clear();
}

void WaterDropletContainer::printPDB(){
	for (int i = 0; i< (int)waterDropletVector.size(); i++){
		waterDropletVector[i].printPDB();
	}
};

void WaterDropletContainer::add( WaterDroplet & waterDroplet) {
    waterDroplet.validate();
    if (hasChainID(waterDroplet.chainID)) {ErrorManager::instance <<__FILE__<<":"<<__LINE__<<"There is already a water droplet with chain ID "<<waterDroplet.chainID<<endl; ErrorManager::instance.treatError();}
    waterDropletVector.push_back(waterDroplet);
};


void WaterDropletContainer::multiplySmallGroupInertia(double smallGroupInertiaMultiplier, CompoundSystem & system, SimbodyMatterSubsystem& matter, State & state) {
	for (int i = 0; i< (int)waterDropletVector.size(); i++){
		waterDropletVector[i].multiplySmallGroupInertia(smallGroupInertiaMultiplier, system,   matter,state);
	}
};

void WaterDropletContainer::addWaterMolecules( CompoundSystem & system, DuMMForceFieldSubsystem &dumm,  SimbodyMatterSubsystem& matter, BiopolymerClassContainer & myBiopolymerClassContainer) {
	for (int i = 0; i< (int)waterDropletVector.size(); i++){
		waterDropletVector[i].addWaterMolecules(  system, dumm,  matter, myBiopolymerClassContainer);
	}
};

void WaterDropletContainer::addTethers(AtomSpringContainer & atomSpringContainer) {
	for (int i = 0; i< (int)waterDropletVector.size(); i++){
		waterDropletVector[i].addTethers(   atomSpringContainer);
	}
};


WaterDroplet WaterDropletContainer::getWaterDroplet(string chainID){
    for (int i = 0; i < (int)waterDropletVector.size(); i++) {
        if (waterDropletVector[i].chainID.compare("chainID") == 0){
            return waterDropletVector[i];
        }
    }
    ErrorManager::instance <<__FILE__<<" : "<<__LINE__<<" Failed to find water droplet with chain ID "<<chainID<<endl;
    ErrorManager::instance.treatError();
}


WaterDroplet & WaterDropletContainer::updWaterDroplet(string chainID){
    for (int i = 0; i < (int)waterDropletVector.size(); i++) {
        if (waterDropletVector[i].chainID.compare(chainID) == 0){
            waterDropletVector[i].validateWaterVector();
            return waterDropletVector[i];
        }
    } 
    ErrorManager::instance <<__FILE__<<" : "<<__LINE__<<" Failed to find water droplet with chain ID "<<chainID<<endl;
    ErrorManager::instance.treatError();
}

bool WaterDropletContainer::hasChainID     (string chainID){
    bool hasChain = false;
    for (int i = 0; i < (int)waterDropletVector.size(); i++) {
       
        if ((waterDropletVector[i].chainID).compare(chainID) == 0){
            hasChain = true;
        }
    }
    return hasChain;
}

MobilizedBody & WaterDroplet::updAtomMobilizedBody(SimbodyMatterSubsystem & matter, ResidueID myResidueNumber, string myAtomName){
    Compound::AtomIndex myAtomIndex = waterVector[myResidueNumber.getResidueNumber()-firstResidueNumber.getResidueNumber()].getAtomIndex(myAtomName);
    MobilizedBodyIndex myAtomMobilizedBodyIndex = waterVector[myResidueNumber.getResidueNumber()-firstResidueNumber.getResidueNumber()].getAtomMobilizedBodyIndex(myAtomIndex);
    return matter.updMobilizedBody(myAtomMobilizedBodyIndex); 
}

void WaterDroplet::multiplySmallGroupInertia(double smallGroupInertiaMultiplier, CompoundSystem & system, SimbodyMatterSubsystem& matter, State & state) {
	for (int i = 0; i< (int)waterVector.size(); i++){
            MobilizedBody myBody = updAtomMobilizedBody(matter,ResidueID(i+firstResidueNumber.getResidueNumber(),' '),string("OW"));
            MassProperties myBodyMassProperties = myBody.getBodyMassProperties(state);
            if (myBodyMassProperties.getMass() < 40){
                //cout<<__FILE__<<__LINE__<<" : "<<myBody.getBodyMassProperties(state).getInertia()<<endl;
                //cout<<__FILE__<<__LINE__<<" : "<< smallGroupInertiaMultiplier <<endl;
                myBody.setDefaultMassProperties (MassProperties(myBodyMassProperties.getMass(), myBodyMassProperties.getMassCenter(),  myBodyMassProperties.getInertia ()* smallGroupInertiaMultiplier));
                state = system.realizeTopology();
                system.realize(state,Stage::Position);
                //cout<<__FILE__<<__LINE__<<" : "<<myBody.getBodyMassProperties(state).getInertia()<<endl;
            } // of if 
	} // of for
}; // of multiplySmallGroupInertia method

void WaterDroplet::addTethers(AtomSpringContainer & atomSpringContainer) {
    if (tetherType.compare("ToGround") == 0) for (int i = 0; i < (int)waterVector.size(); i++) {
        AtomSpring myAtomSpring;
        myAtomSpring.tether = true ;

                // format:  springToGround atom1Chain atom1Residue atom1Name  groundX groundY groundZ       deadLength forceConstant
        myAtomSpring.atom1Chain   = chainID;    
        myAtomSpring.atom1Residue =  ResidueID((i+1),' ');
        myAtomSpring.atom1Name    = "OW";        
        myAtomSpring.atom2Chain   = string("");
        myAtomSpring.atom2Residue = ResidueID();//0.0;        
        myAtomSpring.atom2Name    = string("");
        myAtomSpring.groundLocation[0]    = center[0]; // used to convert Å to nm, no longer done.
        myAtomSpring.groundLocation[1]    = center[1];
        myAtomSpring.groundLocation[2]    = center[2];
        myAtomSpring.toGround     = true ;
        myAtomSpring.deadLength   = 1.1*getRadius(); // used to convert Å to nm, no longer done.
        myAtomSpring.forceConstant= tetherStrength;
        //cout<<__FILE__<<__LINE__<<" : "<<myAtomSpring.groundLocation<<endl;
        atomSpringContainer.add(myAtomSpring);        
    } // of for int i
} // of addTethers


MobilizedBodyIndex WaterDroplet::getOxygenMobilizedBodyIndex(ResidueID residueNumber) {
    if (! waterVector[residueNumber.getResidueNumber()-firstResidueNumber.getResidueNumber()].hasAtom("OW")) {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<"The requested atom was not found!"<<endl; ErrorManager::instance.treatError();
    }
    Compound::AtomIndex myAtomIndex = waterVector[residueNumber.getResidueNumber()-firstResidueNumber.getResidueNumber()].getAtomIndex("OW" );
    //cout<<myAtomIndex<<endl;
    return waterVector[residueNumber.getResidueNumber()-firstResidueNumber.getResidueNumber()].getAtomMobilizedBodyIndex( myAtomIndex);
}



MobilizedBody  WaterDroplet::getOxygenMobilizedBody(SimbodyMatterSubsystem & matter, ResidueID   residueNumber){
    MobilizedBodyIndex myAtomMobilizedBodyIndex = getOxygenMobilizedBodyIndex(residueNumber);
    return matter.updMobilizedBody(myAtomMobilizedBodyIndex);
}


MobilizedBody  & WaterDroplet::updOxygenMobilizedBody(SimbodyMatterSubsystem & matter, ResidueID   residueNumber){
    //cout<<__FILE__<<" : "<<__LINE__<<endl;
    validateWaterVector();
    MobilizedBodyIndex myAtomMobilizedBodyIndex = getOxygenMobilizedBodyIndex(residueNumber);
    return matter.updMobilizedBody(myAtomMobilizedBodyIndex);
}


Vec3 WaterDroplet::getOxygenLocationInMobilizedBodyFrame(ResidueID residueNumber) {
   return Vec3(0,0,0); 
}

void WaterDroplet::includeAllAtoms(DuMMForceFieldSubsystem & dumm){
    for (int i = 0; i < (int)waterVector.size(); i++) {
        Compound::AtomIndex myAtomIndex = waterVector[i].getAtomIndex("OW");
        DuMM::AtomIndex myDuMMAtomIndex = waterVector[i].getDuMMAtomIndex(myAtomIndex);
        dumm.includeNonbondAtom(myDuMMAtomIndex);

        myAtomIndex = waterVector[i].getAtomIndex("HW0");
        myDuMMAtomIndex = waterVector[i].getDuMMAtomIndex(myAtomIndex);
        dumm.includeNonbondAtom(myDuMMAtomIndex);

        myAtomIndex = waterVector[i].getAtomIndex("HW1");
        myDuMMAtomIndex = waterVector[i].getDuMMAtomIndex(myAtomIndex);
        dumm.includeNonbondAtom(myDuMMAtomIndex);
    }
}

void WaterDroplet::validateWaterVector() {
    for (int i =0 ; i < (int)waterVector.size() ; i++) {
        cout<<__FILE__<<" : "<<__LINE__<<" : Validating water #"<<i<<", chain "<<chainID<<" : "<<waterVector[i].hasAtom("OW")<<endl;;
        if (! waterVector[i].hasAtom("OW")) {

            cout<<__FILE__<<" : "<<__LINE__<<" : Failed to find atom OW in waterVector element "<<i<<endl;;
        }
    }
}

void WaterDropletContainer::includeAllAtoms(DuMMForceFieldSubsystem & dumm){
    for (int i = 0; i< (int)waterDropletVector.size(); i++){
        waterDropletVector[i].includeAllAtoms(dumm);
    }
}

void WaterDropletContainer::validateWaterVectors(){
    for (int i = 0; i< (int)waterDropletVector.size(); i++){
        waterDropletVector[i].validateWaterVector();
    }
}

void WaterDroplet::addWaterMolecules( CompoundSystem & system, DuMMForceFieldSubsystem &dumm,  SimbodyMatterSubsystem& matter, BiopolymerClassContainer & myBiopolymerClassContainer) {

	waterVector.clear();
    	Vec3 atomLocation;
	//Compound::AtomPathName myAtomName;
        double dropletRadius = getRadius() ;  // used to convert Ångstroms to nanometers
        double arista = .310433899; 
	int cellsAcross = (int)(getRadius()*2/arista+1); // number of cells across.  add 1 b.c. water molecules are at corners of cells.

        cout<<__FILE__<<":"<<__LINE__<<" Droplet radius = "<<dropletRadius<<" nm, cell width (arista) = "<<arista<<" nm"<<endl;
	std::cout<<"cellsAcross = "<<cellsAcross<<std::endl;

        Vec3 dropletCenter = Vec3 (center[0], center[1], center[2]); // used to convert Ångstroms to nanometers

	int totalAtoms =0;

        vector< vector < vector<int> > > myWaterVecOccupation( cellsAcross, vector < vector<int> >(cellsAcross, vector<int>(cellsAcross) ) );
	int numWater =0;//lsAcross*cellsAcross*cellsAcross;

	double myRadius=0.;
        for (int xi = 0; xi < cellsAcross; xi++)
        	for (int yi = 0; yi < cellsAcross; yi++)
        		for (int zi = 0; zi < cellsAcross; zi++)
                	{
				myRadius = pow((pow((xi*arista-dropletRadius),2.)+ pow((yi*arista-dropletRadius),2.)+ pow((zi*arista-dropletRadius),2.)),.5 );
				if (myRadius <= dropletRadius) //(pow((pow((dropletCenter[0]-(xi*arista+dropletRadius)),2)+ pow((dropletCenter[1]-(yi*arista+dropletRadius)),2)+ pow((dropletCenter[2]-(zi*arista+dropletRadius)),2)),.5 ) <= dropletRadius)  
					{
						
					numWater++;myWaterVecOccupation[xi][yi][zi] = 1;
					}
				else
		  		    {//cout<<"outside of water droplet radius; no water to be placed here."<<endl; 
                                    myWaterVecOccupation[xi][yi][zi] = 0;}       
			}	
	//Count the atoms in the CompoundSystem
	
        for (int i =0 ; i< (int)system.getNumCompounds(); i++)
	{
		cout<<"counting atoms in compound # "<<i<<endl;	
		totalAtoms += system.updCompound(CompoundSystem::CompoundIndex(i)).getNumAtoms();
	}	
	cout<<"counted "<<totalAtoms<<" total atoms."<<endl;
	//Create an array to hold their locations
	vector<Vec3> atomLocationArray(totalAtoms);
 
        vector<MMBAtomInfo> myAtomInfoVector = myBiopolymerClassContainer.getConcatenatedAtomInfoVector();
	for (int i =0 ; i< (int)system.getNumCompounds(); i++)
	{
	    for (int j = 0; j< system.updCompound(CompoundSystem::CompoundIndex(i)).getNumAtoms(); j++) 
                if (! (myBiopolymerClassContainer.hasChainID(system.updCompound(CompoundSystem::CompoundIndex(i)).getPdbChainId())))
	        {
                        MMBAtomInfo myAtomInfo  ;
			Compound & myCompound = system.updCompound(CompoundSystem::CompoundIndex(i));
			myAtomInfo.atomName = myCompound.getAtomName(Compound::AtomIndex(j));
			myAtomInfo.residueID = ResidueID(myCompound.getPdbResidueNumber(),' ');
                        cout<<__FILE__<<":"<<__LINE__<<" loading myAtomInfoVector with : "<<endl;
                        myAtomInfo.print();
                        Vec3 tempVec3 = myCompound.calcDefaultAtomLocationInGroundFrame(myAtomInfo.atomName);
			myAtomInfo.position = openmmVecType(tempVec3[0], tempVec3[1], tempVec3[2]);
                        myAtomInfoVector.push_back(myAtomInfo);
                }
        }
   
        if (totalAtoms != myAtomInfoVector.size()) {
	    cout<<"system had "<<totalAtoms<<" total atoms. myAtomInfoVector has "<<myAtomInfoVector.size()<<" elements. "<<endl;
	    ErrorManager::instance<<__FILE__<<":"<<__LINE__<<": Error! wrong number of elements in myAtomInfoVector! "<<endl;
	    ErrorManager::instance.treatError();
        }
	//int myAtomIndex = 0;
        for (int i = 0; i < myAtomInfoVector.size(); i++) 
	//for (int i =0 ; i< (int)system.getNumCompounds(); i++)
	{
	    //for (int j = 0; j< system.updCompound(CompoundSystem::CompoundIndex(i)).getNumAtoms(); j++)
                //if (system.updCompound(CompoundSystem::CompoundIndex(i)).getPdbChainId())    
		//{
                        MMBAtomInfo myAtomInfo = myAtomInfoVector[i]  ;
                        
			cout << "my Atom Name= "<<myAtomInfo.atomName <<" atomLocation= "<<myAtomInfo.position<<endl;
			int myxi = (int)((myAtomInfo.position[0]-dropletCenter[0]+dropletRadius)/arista+.5);
			int myyi = (int)((myAtomInfo.position[1]-dropletCenter[1]+dropletRadius)/arista+.5);
			int myzi = (int)((myAtomInfo.position[2]-dropletCenter[2]+dropletRadius)/arista+.5);
			if (((myxi < cellsAcross) && (myxi >=0))
				&& ((myyi < cellsAcross) && (myyi >=0))
				&& ((myzi < cellsAcross) && (myzi >=0))
				&& (myWaterVecOccupation[myxi][myyi][myzi] == 1))
			{
					numWater--;
					myWaterVecOccupation[myxi][myyi][myzi] = 0;			
			}
			//myAtomIndex++;
		//}
	}
      
        waterVector.clear();// this is redundant and should be removed; is done in parameterreader.cpp
        for (int k =0;k<numWater;k++)	{
		Water myWaterMolecule(dumm);
		(myWaterMolecule).setPdbResidueName("H2O");
		(myWaterMolecule).setPdbResidueNumber(k+1  );
		(myWaterMolecule).setPdbChainId(chainID); 
                waterVector.push_back(myWaterMolecule);
                cout<<__FILE__<<":"<<__LINE__<<" Just added water with residue number: "<<(k+1)<<endl;
	}
        validateWaterVector();
        cout<<__FILE__<<":"<<__LINE__<<endl;
	int myWaterIndex =-1;
        for (int xi = 0; xi < cellsAcross; xi++) {
        	for (int yi = 0; yi < cellsAcross; yi++)
                    {
        		for (int zi = 0; zi < cellsAcross; zi++)
			{
				cout<<"checking occupation for "<<xi<<","<<yi<<","<<zi<<endl;
				if (myWaterVecOccupation[xi][yi][zi]) 
				{
					myWaterIndex++;
					if (myWaterIndex > numWater) {ErrorManager::instance<< "Attempted to adopt more water than exists in array"<<endl; ErrorManager::instance.treatError();}
				  	//cout<<"adopting water at "<< xi*arista-dropletRadius+dropletCenter[0]<<","<<yi*arista-dropletRadius+dropletCenter[1]<<"=,"<<zi*arista-dropletRadius+dropletCenter[2]<<endl;
					Rotation myRotation(45.*Deg2Rad, YAxis);
					myRotation =  myRotation * Rotation(45.*Deg2Rad, XAxis);
					Transform myTransform(myRotation,Vec3(xi*arista-dropletRadius+dropletCenter[0],yi*arista-dropletRadius+dropletCenter[1],zi*arista-dropletRadius+dropletCenter[2]));
                                        //defaultTransformVector.push_back(myTransform);
                                        waterVector[myWaterIndex].setTopLevelTransform(myTransform);
				        //system.adoptCompound(waterVector[myWaterIndex],myTransform); // This should be moved to its own method

				}	
			} // of for zi

                    } // of for yi
     
		    } // of for xi
	cout<<__FILE__<<":"<<__LINE__<<endl;	   
	

	};

void WaterDroplet::adopt( CompoundSystem & system, bool readPreviousFrameFile) {		      
    /*if (!(readPreviousFrameFile)) {
        for (int i = 0 ; i < waterVector.size(); i++) {
            cout <<__FILE__<<":"<<__LINE__<<waterVector[i].getTopLevelTransform();
            system.adoptCompound(waterVector[i],defaultTransformVector[i]); // This should be moved to its own method
            cout <<__FILE__<<":"<<__LINE__<<waterVector[i].getTopLevelTransform();}
    } else {*/
        for (int i = 0 ; i < waterVector.size(); i++){
            cout <<__FILE__<<":"<<__LINE__<<waterVector[i].getTopLevelTransform();
            system.adoptCompound(waterVector[i]); // This should be moved to its own method
            cout <<__FILE__<<":"<<__LINE__<<waterVector[i].getTopLevelTransform();}
    //}
}

void WaterDropletContainer::matchDefaultConfiguration(bool readPreviousFrameFile, String pdbFileName,bool matchExact, bool matchIdealized)
{
    if (readPreviousFrameFile) {
        PdbStructure pdbStructure;
        
        //============================================ Read in PDB or CIF
        if ( pdbFileName.length() > 4 )
        {
            if ( pdbFileName.substr ( pdbFileName.length() - 4, pdbFileName.length() - 1) == ".pdb" )
            {
                std::ifstream inputFile               ( pdbFileName.c_str(), ifstream::in );
                
                if ( !inputFile.good() )
                {
                    ErrorManager::instance << "!!! Error !!! The file " << pdbFileName << " could not be opened. If this is not the file you wanted to open, please supply the requested file name after the loadSequencesFromPdb command. Note that the supported file extensions currently are \".pdb\", \".cif\" and \".cif.gz\"." << std::endl;
                    ErrorManager::instance.treatError ( );
                }
                else
                {
                    pdbStructure                      = PdbStructure (inputFile);
                }
            }
            else if ( pdbFileName.substr ( pdbFileName.length() - 4, pdbFileName.length() - 1) == ".cif" )
            {
#ifdef GEMMI_USAGE
                std::ifstream testOpen                ( pdbFileName.c_str() );
                if ( testOpen.good() )
                {
                    pdbStructure                      = PdbStructure (pdbFileName);
                }
                else
                {
                    std::string pdbFileHlp            = pdbFileName;
                    pdbFileHlp.append                 ( ".gz" );
                    
                    std::ifstream testOpen2           ( pdbFileHlp.c_str() );
                    if ( testOpen2.good() )
                    {
                        pdbStructure                  = PdbStructure ( pdbFileHlp );
                    }
                    else
                    {
                        ErrorManager::instance << "!!! Error !!! The file " << pdbFileName << " could not be opened. If this is not the file you wanted to open, please supply the requested file name after the loadSequencesFromPdb command. Note that the supported file extensions currently are \".pdb\", \".cif\" and \".cif.gz\"." << std::endl;
                        ErrorManager::instance.treatError ( );
                    }
                    testOpen2.close                   ( );
                }
                testOpen.close                        ( );
#else
                ErrorManager::instance << "!!! Error !!! MMB was not compiled with the Gemmi library required for mmCIF support. Cannot proceed, if you want to use mmCIF files, please re-compile with the Gemmi library option allowed." << std::endl;
                ErrorManager::instance.treatError     ( );
#endif
            }
            else if ( pdbFileName.length() > 7 )
            {
                if ( pdbFileName.substr ( pdbFileName.length() - 7, pdbFileName.length() - 1) == ".cif.gz" )
                {
#ifdef GEMMI_USAGE
                pdbStructure                          = PdbStructure (pdbFileName);
#else
                ErrorManager::instance << "!!! Error !!! MMB was not compiled with the Gemmi library required for mmCIF support. Cannot proceed, if you want to use mmCIF files, please re-compile with the Gemmi library option allowed." << std::endl;
                ErrorManager::instance.treatError     ( );
#endif
                }
                else
                {
                    ErrorManager::instance << "!!! Error !!! The file " << pdbFileName << " could not be opened. If this is not the file you wanted to open, please supply the requested file name after the loadSequencesFromPdb command. Note that the supported file extensions currently are \".pdb\", \".cif\" and \".cif.gz\"." << std::endl;
                    ErrorManager::instance.treatError ( );
                }
            }
            else
            {
                ErrorManager::instance << "!!! Error !!! The file " << pdbFileName << " could not be opened. If this is not the file you wanted to open, please supply the requested file name after the loadSequencesFromPdb command. Note that the supported file extensions currently are \".pdb\", \".cif\" and \".cif.gz\"." << std::endl;
                ErrorManager::instance.treatError     ( );
            }
        }
        else
        {
            ErrorManager::instance << "!!! Error !!! The file " << pdbFileName << " could not be opened. If this is not the file you wanted to open, please supply the requested file name after the loadSequencesFromPdb command. Note that the supported file extensions currently are \".pdb\", \".cif\" and \".cif.gz\"." << std::endl;
            ErrorManager::instance.treatError         ( );
        }
        
        cout <<__FILE__<<":"<<__LINE__<<endl;
	for (int i = 0; i< (int)waterDropletVector.size(); i++){
  
            ////////////////////////////////////////////////
            // First, look in the structure file .. if there are more waters in there than expected, expand the waterDropletVector[i].waterVector accordingly:
            ////////////////////////////////////////////////
            bool keepGoing = true;
            int k = waterDropletVector[i].waterVector.size()-1;
	    Water myExtraWater = waterDropletVector[i].waterVector[k]; // Note that we are making a COPY of the last water in the vector;
            while (keepGoing) {
                cout <<__FILE__<<":"<<__LINE__<<" extra water PDB residue number = "<<myExtraWater.getPdbResidueNumber()<<endl;
		myExtraWater.setPdbResidueNumber(myExtraWater.getPdbResidueNumber()+1);    
                cout <<__FILE__<<":"<<__LINE__<<" extra water PDB residue number = "<<myExtraWater.getPdbResidueNumber()<<endl;
                cout <<__FILE__<<":"<<__LINE__<<" myExtraWater.createAtomTargets(pdbStructure).size() = "<< myExtraWater.createAtomTargets(pdbStructure).size() <<endl;
		if (myExtraWater.createAtomTargets(pdbStructure).size() > 0) {
		    waterDropletVector[i].waterVector.push_back(myExtraWater);} else keepGoing = false;
            }
            ////////////////////////////////////////////////
	    for (int j = 0 ; j <   waterDropletVector[i].waterVector.size(); j++) {
               cout <<__FILE__<<":"<<__LINE__<<endl;
	       Water & myWater = waterDropletVector[i].waterVector[j]; 
	       cout <<__FILE__<<":"<<__LINE__<<" About to create atom targets from file "<<pdbFileName<<endl;
	       cout <<__FILE__<<":"<<__LINE__<<myWater.getPdbChainId()<<endl;
	       Compound::AtomTargetLocations atomTargets = myWater.createAtomTargets(pdbStructure);
               if ((j == 0) && (atomTargets.size() == 0)) {
	           cout <<__FILE__<<":"<<__LINE__<<" No water molecules found for   droplet "<<i<<" in "<<pdbFileName<<" .. using computed coordinates."<<endl;
                   return; // There is no initial structure. Just leave the initial atom positions as they are (probably in a regular lattice), and quit this function.
               } else if ((j > 0) && (atomTargets.size() == 0)) { // This is the case that vector has a water molecule which is not in the structure file.  deleting this and any sbusequent water molecules.
	           cout <<__FILE__<<":"<<__LINE__<<" No coordinates found for droplet "<<i<<" molecule "<<j<<" in "<<pdbFileName<<" .. deleting this and subsequent waters."<<endl;
                   waterDropletVector[i].waterVector.erase(waterDropletVector[i].waterVector.begin()+j,waterDropletVector[i].waterVector.end()) ;
                   return;                  
               }
               cout <<__FILE__<<":"<<__LINE__<<endl;

               for (map<Compound::AtomIndex, Vec3> ::iterator it= atomTargets.begin(); it !=atomTargets.end(); it++) {   
                   cout <<__FILE__<<":"<<__LINE__<<it->second<<endl;
               }
	       if (matchExact)
			{
			//cout<<__FILE__<<":"<<__LINE__<<" Top level transform before fitting: "<<myWater.getTopLevelTransform()<<endl;
			myWater.matchDefaultConfiguration(atomTargets,   Compound::Match_Exact );
			//cout<<__FILE__<<":"<<__LINE__<<" Top level transform after fitting: "<<waterDropletVector[i].getTopLevelTransform()<<endl;
			}
	       if (matchIdealized)
			{myWater.matchDefaultConfiguration(atomTargets,   Compound::Match_Idealized );} //planarity tolerance is in Radians, according to Sherm
               cout <<__FILE__<<":"<<__LINE__<<myWater.getTopLevelTransform();
	    }
            cout <<__FILE__<<":"<<__LINE__<<endl;
    }
    cout <<__FILE__<<":"<<__LINE__<<endl;
}  
}

void WaterDropletContainer::adopt(CompoundSystem & system,bool readPreviousFrameFile) {		      
    for (int i = 0; i< (int)waterDropletVector.size(); i++){
        waterDropletVector[i].adopt(system,readPreviousFrameFile);
    }
}

