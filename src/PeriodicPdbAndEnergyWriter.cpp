/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *         Alex Tek                                                           *  
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */
#include <cstdio>
#include <iostream>
#include <iomanip>

#include "PeriodicPdbAndEnergyWriter.h"
#include "SimTKsimbody.h"

#ifdef CPP4_MAPS_USAGE
  #include <mmdb2/mmdb_manager.h>
  #include <mmdb2/mmdb_cifdefs.h>
#endif

using namespace std;

std::vector<double> PeriodicPdbAndEnergyWriter::myEnergies;

SimTK::PeriodicPdbAndEnergyWriter::PeriodicPdbAndEnergyWriter( 	const CompoundSystem& system,
    															const DuMMForceFieldSubsystem& dumm  , 
    															std::ostream& outputStream,
    															double interval,
    															ParameterReader & myParameterReader,
    															BiopolymerClassContainer & myBiopolymerClassContainer
    															//vector<MagnesiumIon> myMagnesiumIonVec
    														 ) : PeriodicEventHandler(interval),
																 system(system), 
																 dumm(dumm), // the second dumm comes from the parameters fed to the constructor.  the first is the private member.
																 outputStream(outputStream) ,
																 myParameterReader(myParameterReader),
                                                                                                                                 myBiopolymerClassContainer(myBiopolymerClassContainer)
																 //myMagnesiumIonVec(myMagnesiumIonVec)
{
#ifdef CPP4_MAPS_USAGE
    if ( myParameterReader.useCIFFileFormat )
    {
        //============================================ Initialise MMDB2 Manager object for CIF writing, if need be
        mmdb2Manager                                  = new mmdb::Manager ( );
    }
    else
    {
        //============================================ Leave MMDB2 Manager empty, as it will not be used for this run
        mmdb2Manager                                  = nullptr;
    }
#endif
}

SimTK::PeriodicPdbAndEnergyWriter::~PeriodicPdbAndEnergyWriter ( )
{
#ifdef CPP4_MAPS_USAGE
    if ( mmdb2Manager != nullptr )
    {
        //============================================ Delete MMDB2 Manager object for CIF writing, if need be
        delete mmdb2Manager;
    }
#endif
}

void SimTK::PeriodicPdbAndEnergyWriter::handleEvent(State& state, Real accuracy, bool& shouldTerminate) const  {

    static int modelNumber = 1; // increments by one at each reporting step
    int compoundNumber                                = 1;
    int requiredPrecision                             = 12;
    
    system.realize(state, Stage::Dynamics);
#ifdef CPP4_MAPS_USAGE
    mmdb::Model *mmdb2Model;
#endif
    if ( myParameterReader.useCIFFileFormat )
    {
#ifdef CPP4_MAPS_USAGE
        mmdb2Model                                    = new mmdb::Model ( mmdb2Manager, modelNumber );
        mmdb2Manager->MakeBiomolecule                 ( 1, modelNumber );
#else
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Error! Requested mmCIF file output, but did not compile with the MMDB2 library. Cannot proceed, if you want to use mmCIF files, please re-compile with the MMDB2 library option allowed." <<endl;
        ErrorManager::instance.treatError             ( );
#endif
    }
    else
    {
        outputStream << "MODEL     " << std::setw(4) << modelNumber << std::endl;
    }
    PdbAtom::setWriteFullPrecisionLocation(myParameterReader.writeDoublePrecisionTrajectories);// false by default to save a little disk space here

    if (myParameterReader.contactInterfaceContainer.numInterfaces() > 0) {
	vector<MMBAtomInfo> concatenatedAtomInfoVector = myBiopolymerClassContainer.getConcatenatedAtomInfoVector(state);
        #ifdef USE_OPENMM  
	vector<TwoAtomClass> myTwoAtomClassVector = myParameterReader.contactInterfaceContainer.retrieveCloseContactPairs(concatenatedAtomInfoVector); // Just calling this method is enough to get the close contacts written to stdout.
        #endif  
    }

    for (SimTK::CompoundSystem::CompoundIndex c(0); c < system.getNumCompounds(); ++c)
    {
        if ( myParameterReader.useCIFFileFormat )
        {
#ifdef CPP4_MAPS_USAGE
            //======================================== Re-open file for writing
            myParameterReader.cifTrajectoryFile.assign ( myParameterReader.outTrajectoryFileName.c_str ( ) );
            if ( !myParameterReader.cifTrajectoryFile.rewrite ( ) )
            {
                ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Failed to open the file "<< myParameterReader.outTrajectoryFileName << " for writing." << std::endl;
                ErrorManager::instance.treatError     ( );
            }
            
            //======================================== Build the MMDB2 structure
            (system.getCompound(c)).buildCif          ( state, mmdb2Model, Transform(Vec3(0) ) );
            
            //======================================== Add the model to the data
            for ( int noModel = 0; noModel < mmdb2Model->GetNumberOfChains(); noModel++ )
            {
                if ( mmdb2Model->GetChain(noModel) )
                {
                    for ( int noRes = 0; noRes < mmdb2Model->GetChain(noModel)->GetNumberOfResidues(); noRes++ )
                    {
                        if ( mmdb2Model->GetChain(noModel)->GetResidue(noRes) )
                        {
                            for ( int noAt = 0; noAt < mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetNumberOfAtoms(); noAt++ )
                            {
                                if ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt) )
                                {
                                    //==================== Initialise variables
                                    mmdb::mmcif::PLoop Loop;
                                    mmdb::AtomName AtName;
                                    mmdb::Element el;
                                    char N[10];
                                    int i,j,RC;
                                    mmdb::PChain chain  = mmdb2Model->GetChain(noModel)->GetResidue(noRes)->chain;
                                    mmdb::PModel model  = mmdb::PModel ( mmdb2Model->GetChain(noModel)->GetModel() );
                                    
                                    //================ Initialise the Loop object
                                    RC                = myParameterReader.cifData.AddLoop ( mmdb::CIFCAT_ATOM_SITE, Loop );
                                    if ( RC != mmdb::mmcif::CIFRC_Ok )
                                    {
                                        //============ The category was (re)created, provide tags
                                        Loop->AddLoopTag ( mmdb::CIFTAG_GROUP_PDB          ); // ATOM, TER etc.
                                        Loop->AddLoopTag ( mmdb::CIFTAG_ID                 ); // serial number

                                        Loop->AddLoopTag ( mmdb::CIFTAG_TYPE_SYMBOL        ); // element symbol
                                        Loop->AddLoopTag ( mmdb::CIFTAG_LABEL_ATOM_ID      ); // atom name
                                        Loop->AddLoopTag ( mmdb::CIFTAG_LABEL_ALT_ID       ); // alt location
                                        Loop->AddLoopTag ( mmdb::CIFTAG_LABEL_COMP_ID      ); // residue name
                                        Loop->AddLoopTag ( mmdb::CIFTAG_LABEL_ASYM_ID      ); // chain ID
                                        Loop->AddLoopTag ( mmdb::CIFTAG_LABEL_ENTITY_ID    ); // entity ID
                                        Loop->AddLoopTag ( mmdb::CIFTAG_LABEL_SEQ_ID       ); // res seq number
                                        Loop->AddLoopTag ( mmdb::CIFTAG_PDBX_PDB_INS_CODE  ); // insertion code
                                        Loop->AddLoopTag ( mmdb::CIFTAG_SEGMENT_ID         ); // segment ID

                                        Loop->AddLoopTag ( mmdb::CIFTAG_CARTN_X            ); // x-coordinate
                                        Loop->AddLoopTag ( mmdb::CIFTAG_CARTN_Y            ); // y-coordinate
                                        Loop->AddLoopTag ( mmdb::CIFTAG_CARTN_Z            ); // z-coordinate
                                        Loop->AddLoopTag ( mmdb::CIFTAG_OCCUPANCY          ); // occupancy
                                        Loop->AddLoopTag ( mmdb::CIFTAG_B_ISO_OR_EQUIV     ); // temp factor

                                        Loop->AddLoopTag ( mmdb::CIFTAG_CARTN_X_ESD        ); // x-sigma
                                        Loop->AddLoopTag ( mmdb::CIFTAG_CARTN_Y_ESD        ); // y-sigma
                                        Loop->AddLoopTag ( mmdb::CIFTAG_CARTN_Z_ESD        ); // z-sigma
                                        Loop->AddLoopTag ( mmdb::CIFTAG_OCCUPANCY_ESD      ); // occupancy-sigma
                                        Loop->AddLoopTag ( mmdb::CIFTAG_B_ISO_OR_EQUIV_ESD ); // t-factor-sigma

                                        Loop->AddLoopTag ( mmdb::CIFTAG_PDBX_FORMAL_CHARGE ); // charge on atom

                                        Loop->AddLoopTag ( mmdb::CIFTAG_AUTH_SEQ_ID        ); // res seq number
                                        Loop->AddLoopTag ( mmdb::CIFTAG_AUTH_COMP_ID       ); // residue name
                                        Loop->AddLoopTag ( mmdb::CIFTAG_AUTH_ASYM_ID       ); // chain id
                                        Loop->AddLoopTag ( mmdb::CIFTAG_AUTH_ATOM_ID       ); // atom name

                                        Loop->AddLoopTag ( mmdb::CIFTAG_PDBX_PDB_MODEL_NUM ); // model number
                                    }
                                    
                                    //================ Is this a normal atom record?
                                    if ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->WhatIsSet & ( mmdb::ASET_Coordinates | mmdb::ASET_CoordSigma))
                                    {
                                        //============ Yes!
                                        
                                        // group_PDB field
                                        if ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->Het ) { Loop->AddString ( pstr ( "HETATM" ) ); }
                                        else                                                                        { Loop->AddString ( pstr( "ATOM" ) ); }
                                  
                                        // id field
                                        if ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->serNum > 0 ) { Loop->AddInteger ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->serNum ); }
                                        else                                                                               { Loop->AddInteger ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->GetIndex() ); }
                                        
                                        if ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->WhatIsSet & mmdb::ASET_Coordinates)
                                        {
                                            // type_symbol field
                                            mmdb::strcpy_css ( el, mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->element );
                                            Loop->AddString ( el, true );
                                            
                                            // label_atom_id field
                                            mmdb::strcpy_css ( AtName, mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->label_atom_id );
                                            Loop->AddString ( AtName );
                                            
                                            // label_alt_id field
                                            Loop->AddString  ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->altLoc, true );

                                            // label_comp_id field
                                            Loop->AddString ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->label_comp_id );
                                            
                                            // label_asym_id field
                                            Loop->AddString ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->label_asym_id );
                                            
                                            // label_entity_id field
                                            if ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->label_entity_id > 0 ) { Loop->AddInteger ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->label_entity_id ); }
                                            else                                                                         { Loop->AddNoData  ( mmdb::mmcif::CIF_NODATA_DOT ); }
                                            
                                            // label_seq_id field
                                            if ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->label_seq_id > mmdb::MinInt4 ) { Loop->AddInteger ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->label_seq_id ); }
                                            else                                                                                  { Loop->AddNoData ( mmdb::mmcif::CIF_NODATA_DOT ); }
                                            
                                            // pdbx_PDB_ins_code field
                                            Loop->AddString ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->insCode, true );

                                            // segment_id field
                                            Loop->AddString ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->segID, true );
                                            
                                            // Cartn_x, Cartn_y, Cartn_z fields
                                            Loop->AddReal ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->x, requiredPrecision );
                                            Loop->AddReal ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->y, requiredPrecision );
                                            Loop->AddReal ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->z, requiredPrecision );
                                            
                                            // occupancy field
                                            if ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->WhatIsSet & mmdb::ASET_Occupancy ) { Loop->AddReal ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->occupancy, requiredPrecision ); }
                                            else { Loop->AddNoData ( mmdb::mmcif::CIF_NODATA_QUESTION ); }
                                            
                                            // B_iso_or_equiv field
                                            if ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->WhatIsSet & mmdb::ASET_tempFactor ) { Loop->AddReal ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->tempFactor, requiredPrecision ); }
                                            else { Loop->AddNoData ( mmdb::mmcif::CIF_NODATA_QUESTION ); }

                                            // cartn_x_esd, cartn_y_esd, cartn_z_esd fields
                                            if ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->WhatIsSet & mmdb::ASET_CoordSigma )
                                            {
                                              Loop->AddReal ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->sigX, requiredPrecision );
                                              Loop->AddReal ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->sigY, requiredPrecision );
                                              Loop->AddReal ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->sigZ, requiredPrecision );
                                            } else
                                            {
                                              Loop->AddNoData ( mmdb::mmcif::CIF_NODATA_QUESTION );
                                              Loop->AddNoData ( mmdb::mmcif::CIF_NODATA_QUESTION );
                                              Loop->AddNoData ( mmdb::mmcif::CIF_NODATA_QUESTION );
                                            }
                                            
                                            // occupancy_esd field
                                            if ( ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->WhatIsSet & mmdb::ASET_OccSigma) &&
                                                 ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->WhatIsSet & mmdb::ASET_Occupancy) )
                                            {
                                                  Loop->AddReal ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->sigOcc, requiredPrecision );
                                            }
                                            else { Loop->AddNoData ( mmdb::mmcif::CIF_NODATA_QUESTION ); }
                                            
                                          // B_iso_or_equiv_esd field
                                          if ( ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->WhatIsSet & mmdb::ASET_tFacSigma) &&
                                               ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->WhatIsSet & mmdb::ASET_tempFactor) )
                                          {
                                                Loop->AddReal ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->sigTemp, requiredPrecision );
                                          }
                                          else { Loop->AddNoData ( mmdb::mmcif::CIF_NODATA_QUESTION ); }
                                        }
                                        else
                                        {
                                          for ( int iter = 0; iter < 18; iter++ )
                                          {
                                            Loop->AddNoData ( mmdb::mmcif::CIF_NODATA_QUESTION );
                                          }
                                        }
                                        
                                        // pdbx_formal_charge field
                                        if ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->WhatIsSet & mmdb::ASET_Charge)
                                        {
                                          sprintf ( N, "%+2i", mmdb::mround ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->charge ) );
                                          Loop->AddString ( N, true );
                                        }
                                        else { Loop->AddNoData ( mmdb::mmcif::CIF_NODATA_QUESTION ); }

                                        // auth_seq_id field
                                        if ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->seqNum > mmdb::MinInt4 ) { Loop->AddInteger ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->seqNum ); }
                                        else                                                                            { Loop->AddNoData ( mmdb::mmcif::CIF_NODATA_DOT ); }
                                        
                                        // auth_comp_id field
                                        Loop->AddString ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->name );

                                        // auth_asym_id field
                                        if ( mmdb2Model->GetChain(noModel) ) { Loop->AddString ( mmdb2Model->GetChain(noModel)->GetChainID(), true ); }
                                        else { Loop->AddNoData ( mmdb::mmcif::CIF_NODATA_DOT ); }
                                        
                                        // auth_atom_id field
                                        mmdb::strcpy_css ( AtName, mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->name );
                                        Loop->AddString  ( AtName );
                                        
                                        // pdbx_PDB_model_num field
                                        if ( mmdb2Model->GetSerNum() > 0) { Loop->AddInteger ( mmdb2Model->GetSerNum() ); }
                                        else                              { Loop->AddNoData ( mmdb::mmcif::CIF_NODATA_QUESTION ); }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            
            //==================================== Write out the CIF file. This is basically a copy of the WriteMMCIF function, but since this function is hardcoded to print the data_ line as a second line and this is not compatible with the current mmCIF formatting, I had to copy and make the single change...
            myParameterReader.cifTrajectoryFile.Write     ( pstr("data_") );
            myParameterReader.cifTrajectoryFile.WriteLine ( myParameterReader.cifData.GetDataName() );
            
            for ( int i = 0; i < myParameterReader.cifData.GetNumberOfCategories(); i++ )
            {
                mmdb::mmcif::Category* cat        = myParameterReader.cifData.GetCategory(i);
                if ( cat )
                {
                    cat->WriteMMCIF               ( myParameterReader.cifTrajectoryFile );
                }
            }
    
            //======================================== Write the entity poly seq loop
            (system.getCompound(c)).writeEntityPolySeqLoop ( state, &myParameterReader.cifTrajectoryFile, compoundNumber );
            
            //======================================== And shut the file
            myParameterReader.cifTrajectoryFile.shut  ( );
            
            //======================================== Prepare for next compound
            ++compoundNumber;
#endif
        }
        else
        {
            (system.getCompound(c)).writePdb(state, outputStream,Transform(Vec3(0)));//, nextAtomSerialNumber);
        }
    }
    
    filebuf fb;
    fb.open("frame.pdb",ios::out);
    std::ostream  fbstream (&fb);
    PdbAtom::setWriteFullPrecisionLocation(true); // get higher precision from file that might be reused for reading.
    for (SimTK::CompoundSystem::CompoundIndex c(0); c < system.getNumCompounds(); ++c)
        (system.getCompound(c)).writePdb(state, fbstream,Transform(Vec3(0)));

    //scf added time reporting 
    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );

    if ( myParameterReader.useCIFFileFormat )
    {
#ifdef CPP4_MAPS_USAGE
        //======================================== Re-open file for writing
        myParameterReader.cifTrajectoryFile.assign ( myParameterReader.outTrajectoryFileName.c_str ( ) );
        if ( !myParameterReader.cifTrajectoryFile.rewrite ( ) )
        {
            ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Failed to open the file "<< myParameterReader.outTrajectoryFileName << " for writing." << std::endl;
            ErrorManager::instance.treatError     ( );
        }
        
        //============================================ Initialise variables
        std::stringstream hlpRem;
        hlpRem.clear();
        hlpRem.str("");

        //============================================ Create the closing remarks
        hlpRem << "REMARK seconds since January 1st 1970: " << time ( NULL );
        mmdb::Remark closRem1                         = mmdb::Remark ( hlpRem.str().c_str() );
        closRem1.remarkNum                            = myParameterReader.mmbRemarkNum;
        closRem1.MakeCIF                              ( &myParameterReader.cifData, myParameterReader.mmbRemarkCounter );
        myParameterReader.mmbRemarkCounter++;
        
        hlpRem.clear();
        hlpRem.str("");
        hlpRem << "REMARK Current time is: " << asctime (timeinfo);
        mmdb::Remark closRem2                         = mmdb::Remark ( hlpRem.str().c_str() );
        closRem2.remarkNum                            = myParameterReader.mmbRemarkNum;
        closRem2.MakeCIF                              ( &myParameterReader.cifData, myParameterReader.mmbRemarkCounter );
        myParameterReader.mmbRemarkCounter++;
        
        hlpRem.clear();
        hlpRem.str("");
        hlpRem << "REMARK elapsed time: " << (clock()/CLOCKS_PER_SEC);
        mmdb::Remark closRem3                         = mmdb::Remark ( hlpRem.str().c_str() );
        closRem3.remarkNum                            = myParameterReader.mmbRemarkNum;
        closRem3.MakeCIF                              ( &myParameterReader.cifData, myParameterReader.mmbRemarkCounter );
        myParameterReader.mmbRemarkCounter++;
        
        //============================================ Write out energies, if computed
        if (myParameterReader.calcEnergy)
        {
            //======================================== Initialise variables
            double myPotentialEnergy                  = system.calcPotentialEnergy(state);
            Real dummPotentialEnergy                  = dumm.calcPotentialEnergy(state);
            double myKineticEnergy                    = system.calcKineticEnergy(state);
            double myEnergy                           = system.calcEnergy(state);
            
            //======================================== Save results
            myParameterReader.kineticEnergy           = myKineticEnergy;
            myParameterReader.totalEnergy             = myEnergy;
            
            //======================================== Write remarks
            hlpRem.clear();
            hlpRem.str("");
            hlpRem << "REMARK Total Potential Energy = " << myPotentialEnergy << " kJ/mol, " << myPotentialEnergy/4.184 << " kcal/mol";
            mmdb::Remark closRem4                     = mmdb::Remark ( hlpRem.str().c_str() );
            closRem4.remarkNum                        = myParameterReader.mmbRemarkNum;
            closRem4.MakeCIF                          ( &myParameterReader.cifData, myParameterReader.mmbRemarkCounter );
            myParameterReader.mmbRemarkCounter++;
            
            hlpRem.clear();
            hlpRem.str("");
            hlpRem << "REMARK DuMM Molecular Dynamics Potential Energy = " << dummPotentialEnergy << " kJ/mol, " << dummPotentialEnergy/4.184 << " kcal/mol";
            mmdb::Remark closRem5                     = mmdb::Remark ( hlpRem.str().c_str() );
            closRem5.remarkNum                        = myParameterReader.mmbRemarkNum;
            closRem5.MakeCIF                          ( &myParameterReader.cifData, myParameterReader.mmbRemarkCounter );
            myParameterReader.mmbRemarkCounter++;
            
            hlpRem.clear();
            hlpRem.str("");
            hlpRem << "REMARK Kinetic Energy = " << myKineticEnergy << " kJ/mol, " << myKineticEnergy/4.184 << " kcal/mol";
            mmdb::Remark closRem6                     = mmdb::Remark ( hlpRem.str().c_str() );
            closRem6.remarkNum                        = myParameterReader.mmbRemarkNum;
            closRem6.MakeCIF                          ( &myParameterReader.cifData, myParameterReader.mmbRemarkCounter );
            myParameterReader.mmbRemarkCounter++;
            
            hlpRem.clear();
            hlpRem.str("");
            hlpRem << "REMARK Energy = " << myKineticEnergy << " kJ/mol, " << myKineticEnergy/4.184 << " kcal/mol";
            mmdb::Remark closRem7                     = mmdb::Remark ( hlpRem.str().c_str() );
            closRem7.remarkNum                        = myParameterReader.mmbRemarkNum;
            closRem7.MakeCIF                          ( &myParameterReader.cifData, myParameterReader.mmbRemarkCounter );
            myParameterReader.mmbRemarkCounter++;
            
            //======================================== Write energies to output log
            cout <<"Total Potential Energy =                   "<<myPotentialEnergy <<" kJ/mol, "<<myPotentialEnergy/4.184<<" kcal/mol"<<std::endl;
            cout <<"DuMM Molecular Dynamics Potential Energy = "<<dummPotentialEnergy <<" kJ/mol, "<<dummPotentialEnergy/4.184<<" kcal/mol"<<std::endl;
            cout <<"Kinetic Energy = "<< myKineticEnergy <<" kJ/mol, "<<myKineticEnergy/4.184<<" kcal/mol"<<std::endl;
            cout <<"REMARK Energy = "<< myEnergy <<" kJ/mol "<<std::endl;
            
            //======================================== Save results
            myEnergies.push_back                      ( myEnergy );
        }
        
        //============================================ Further closing remarks
        hlpRem.clear();
        hlpRem.str("");
        hlpRem << "REMARK Angular, Linear Momentum = " << system.calcSystemRigidBodyMomentum(state);
        mmdb::Remark closRem8                         = mmdb::Remark ( hlpRem.str().c_str() );
        closRem8.remarkNum                            = myParameterReader.mmbRemarkNum;
        closRem8.MakeCIF                              ( &myParameterReader.cifData, myParameterReader.mmbRemarkCounter );
        myParameterReader.mmbRemarkCounter++;
        
        hlpRem.clear();
        hlpRem.str("");
        hlpRem << "REMARK [" << __FILE__ << "] state.getNU()    = " << state.getNU();
        mmdb::Remark closRem9                         = mmdb::Remark ( hlpRem.str().c_str() );
        closRem9.remarkNum                            = myParameterReader.mmbRemarkNum;
        closRem9.MakeCIF                              ( &myParameterReader.cifData, myParameterReader.mmbRemarkCounter );
        myParameterReader.mmbRemarkCounter++;
        
        //==================================== Write out the CIF file. This is basically a copy of the WriteMMCIF function, but since this function is hardcoded to print the data_ line as a second line and this is not compatible with the current mmCIF formatting, I had to copy and make the single change...
        myParameterReader.cifTrajectoryFile.Write     ( pstr("data_") );
        myParameterReader.cifTrajectoryFile.WriteLine ( myParameterReader.cifData.GetDataName() );
        
        for ( int i = 0; i < myParameterReader.cifData.GetNumberOfCategories(); i++ )
        {
            mmdb::mmcif::Category* cat                = myParameterReader.cifData.GetCategory(i);
            if ( cat )
            {
                cat->WriteMMCIF                       ( myParameterReader.cifTrajectoryFile );
            }
        }

        //============================================ Add the poly_seq loop again (it is not written in the Data, but only in the file, so it needs to be do again here.)
        compoundNumber                                = 1;
        for (SimTK::CompoundSystem::CompoundIndex c(0); c < system.getNumCompounds(); ++c) { (system.getCompound(c)).writeEntityPolySeqLoop ( state, &myParameterReader.cifTrajectoryFile, compoundNumber ); ++compoundNumber; }
        myParameterReader.cifTrajectoryFile.shut ( );
#endif
    }
    else
    {
        outputStream << "ENDMDL" << std::endl;
        
        outputStream <<"REMARK seconds since January 1st 1970: "<<time ( NULL     )<<std::endl; //<<"REMARK elapsed time: "<<(clock()/CLOCKS_PER_SEC)<<std::endl;
        outputStream <<"REMARK Current time is: "<<asctime (timeinfo) <<"REMARK elapsed time: "<<(clock()/CLOCKS_PER_SEC)<<std::endl;
        outputStream.setf(ios::fixed, ios::floatfield); // set output to fixed rather than scientific format
        if (myParameterReader.calcEnergy) {
            double myPotentialEnergy = system.calcPotentialEnergy(state);
             //DuMMForceFieldSubsystem& tempDumm = dumm;
            Real dummPotentialEnergy =     dumm.calcPotentialEnergy(state);
            outputStream <<"REMARK Total Potential Energy = "<<myPotentialEnergy <<" kJ/mol, "<<myPotentialEnergy/4.184<<" kcal/mol"<<std::endl;
            outputStream <<"REMARK DuMM Molecular Dynamics Potential Energy = "<<dummPotentialEnergy <<" kJ/mol, "<<dummPotentialEnergy/4.184<<" kcal/mol"<<std::endl;
            cout <<"Total Potential Energy =                   "<<myPotentialEnergy <<" kJ/mol, "<<myPotentialEnergy/4.184<<" kcal/mol"<<std::endl;
            cout <<"DuMM Molecular Dynamics Potential Energy = "<<dummPotentialEnergy <<" kJ/mol, "<<dummPotentialEnergy/4.184<<" kcal/mol"<<std::endl;

            double myKineticEnergy = system.calcKineticEnergy(state);
            myParameterReader.kineticEnergy = myKineticEnergy;
            outputStream <<"REMARK Kinetic Energy = "<< myKineticEnergy <<" kJ/mol, "<<myKineticEnergy/4.184<<" kcal/mol"<<std::endl;
            cout <<"Kinetic Energy = "<< myKineticEnergy <<" kJ/mol, "<<myKineticEnergy/4.184<<" kcal/mol"<<std::endl;

            double myEnergy = system.calcEnergy(state);
            myParameterReader.totalEnergy = myEnergy;
            outputStream <<"REMARK Energy = "<< myEnergy <<" kJ/mol "<<std::endl;
            cout <<"REMARK Energy = "<< myEnergy <<" kJ/mol "<<std::endl;
            myEnergies.push_back(myEnergy);
        }
        outputStream <<"REMARK Angular, Linear Momentum = "<<system.calcSystemRigidBodyMomentum(state)<<endl;

        //cout<<__FILE__<<" : "<<__LINE__<<" "<<dumm.getVdwMixingRuleName (dumm.getVdwMixingRule())<<endl;

        outputStream <<"REMARK ["<< __FILE__<<"] state.getNU()    = "<<state.getNU()            <<std::endl;
        //outputStream <<"REMARK ["<< __FILE__<<"]Satisfied contacts : "<<myParameterReader.satisfiedBasePairs<<endl;
        //outputStream <<"REMARK ["<< __FILE__<<"]Unsatisfied contacts : "<<myParameterReader.unSatisfiedBasePairs<<endl;
    }

    //cout<<__FILE__<<" : "<<__LINE__<<" "<<dumm.getVdwMixingRuleName (dumm.getVdwMixingRule())<<endl;

    outputStream <<"REMARK ["<< __FILE__<<"] state.getNU()    = "<<state.getNU()            <<std::endl;
    //outputStream <<"REMARK ["<< __FILE__<<"]Satisfied contacts : "<<myParameterReader.satisfiedBasePairs<<endl;
    //outputStream <<"REMARK ["<< __FILE__<<"]Unsatisfied contacts : "<<myParameterReader.unSatisfiedBasePairs<<endl;

    cout<<"Just wrote structure for reporting interval # "<<modelNumber<<std::endl; 
    //cout <<"Satisfied base pairs : "<<myParameterReader.satisfiedBasePairs<<" out of "<<myParameterReader.satisfiedBasePairs+myParameterReader.unSatisfiedBasePairs<<endl;
    //cout <<"Unsatisfied contacts : "<<myParameterReader.unSatisfiedBasePairs<<endl;
    ++modelNumber;

    // Check if converged or not
    if (myParameterReader.detectConvergence) {
        if(myEnergies.size() > myParameterReader.convergenceTimeout) {
            double lastEnergy = myEnergies.back();
            double energyDiffMean = 0.0;
            vector<double>::iterator it;
            for(it = myEnergies.end()-myParameterReader.convergenceTimeout; it != myEnergies.end(); it++)
            {
                double diffEnergy = fabs(*it - *(it-1));
                // cout << diffEnergy << " " << endl;
                energyDiffMean += diffEnergy;
            }
            energyDiffMean /= myParameterReader.convergenceTimeout;
            // cout << "Last " << myParameterReader.convergenceTimeout << " energies diff mean: " << energyDiffMean << endl;
            if(energyDiffMean < myParameterReader.convergenceEpsilon)
            {
                myParameterReader.converged = true;   
                cout << "Converged! Energy difference between two reporting intervals has been less than "<< myParameterReader.convergenceEpsilon << " kJ/mol for the last " << myParameterReader.convergenceTimeout << " frames." << endl;
            }
        }
    }
}
