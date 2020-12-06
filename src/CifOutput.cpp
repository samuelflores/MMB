#include "CifOutput.h"

#ifdef GEMMI_USAGE
    #define GEMMI_WRITE_IMPLEMENTATION
    #include <gemmi/to_mmcif.hpp>
    #include <gemmi/to_cif.hpp>
    #include <gemmi/gz.hpp>

void SimTK::CIFOut::writeOutCif ( gemmi::Structure outStruct, std::string fileName, std::vector < std::pair < std::string, std::string > > remarks )
{
    //================================================ Create document from the structure with models
    gemmi::cif::Document outDocument                  = gemmi::make_mmcif_document ( outStruct );
    
    //================================================ Add remarks
    for ( unsigned int blIter = 0; blIter < static_cast<unsigned int> ( outDocument.blocks.size() ); blIter++ )
    {
        if ( remarks.size() != 0 )
        {
            //======================================== Generate single string with the remark
            std::string remText                       = std::string ( "\n; " );
            for ( unsigned int noRem = 0; noRem < static_cast<unsigned int> ( remarks.size() ); noRem++ )
            {
                remText.append                        ( remarks.at(noRem).second );
                remText.append                        ( "\n" );
            }
            remText.append                            ( "\n;" );
            
            //======================================== Write the remark
            outDocument.blocks.at(blIter).init_loop   ( "_database_PDB_remark.", {"id", "text"} ).add_row ( { "3", remText } );
        }
    }
    
    //================================================ Update all blocks
    for ( unsigned int blIter = 0; blIter < static_cast<unsigned int> ( outDocument.blocks.size() ); blIter++ )
    {
        gemmi::update_mmcif_block                       ( outStruct, outDocument.blocks.at(blIter), true );
    }
    
    //================================================ Open output stream
    std::ofstream outputFile;
    outputFile.open                                   ( fileName );
    if ( !outputFile.is_open ( ) )
    {
        std::cerr << "!!! ERROR !!! Failed to open ostream to file " << fileName << ". Maybe issue with write priviledges to the requested location?" << std::endl << std::flush;
        exit                                          ( EXIT_FAILURE );
    }
    
    //================================================ Write out mmCIF file
    gemmi::cif::write_cif_to_stream                   ( outputFile, outDocument, gemmi::cif::Style::Simple );
    outputFile.close                                  ( );
    
    //================================================ Done
    return ;
    
}

void SimTK::CIFOut::reWriteOutCif ( gemmi::Model gModel, std::string modelName, std::string fileName, ParameterReader& myParameterReader, const CompoundSystem& system, bool firstInStage )
{
    //================================================ Try to read in the file
    if ( firstInStage )
    {
        //============================================ File does not exist, create it.
        gemmi::Structure myTrajectoryOutputFile;
        
        //============================================ Add model name
        myTrajectoryOutputFile.name                   = modelName;
        
        //============================================ Save model to structure
        myTrajectoryOutputFile.models.push_back       ( gModel );

        //============================================ Set Gemmi structure internal values based on model information
        gemmi::setup_entities                         ( myTrajectoryOutputFile );
        gemmi::assign_label_seq_id                    ( myTrajectoryOutputFile, true );
        gemmi::assign_subchains                       ( myTrajectoryOutputFile, true );
        
        //============================================ Add sequence to entities
        int compoundNumber                            = 1;
        for (SimTK::CompoundSystem::CompoundIndex c(0); c < system.getNumCompounds(); ++c)
        {
            //======================================== Print log
            std::cout <<__FILE__<<":"<<__LINE__<<" c = "<<c<< " compoundNumber = "<<compoundNumber <<std::endl;
            
            //======================================== Get sequence for compound
            std::string compSeq                       = (myParameterReader.myBiopolymerClassContainer.getBiopolymerClassMap ())[myParameterReader.myBiopolymerClassContainer.updBiopolymerClass(compoundNumber-1).getChainID()].getSequence();
            
            //======================================== For each structure entity
            for ( unsigned int enIt = 0; enIt < static_cast<unsigned int> ( myTrajectoryOutputFile.entities.size() ); enIt++ )
            {
                //==================================== If entity name = MMB chain ID
                if ( myTrajectoryOutputFile.entities.at(enIt).name == myParameterReader.myBiopolymerClassContainer.updBiopolymerClass(compoundNumber-1).getChainID() )
                {
                    //================================ Copy sequence
                    myTrajectoryOutputFile.entities.at(enIt).full_sequence.clear();
                    for ( unsigned int sqIt = 0; sqIt < static_cast<unsigned int> ( compSeq.length() ); sqIt++ )
                    {
                        myTrajectoryOutputFile.entities.at(enIt).full_sequence.emplace_back ( std::to_string ( compSeq[sqIt] ) );
                    }
                }
            }
            
            //======================================== Update compound number
            compoundNumber++;
        }
        
        //============================================ Write out the file
        writeOutCif ( myTrajectoryOutputFile, fileName, myParameterReader.trajectoryFileRemarks );
        
    }
    else
    {
        //============================================ File does exist, read it, update it and re-write it
        gemmi::cif::Document doc                      = gemmi::cif::read ( gemmi::MaybeGzipped ( fileName ) );
        gemmi::Structure myTrajectoryOutputFile       = gemmi::make_structure ( doc );
        
        //============================================ Add model to structure
        myTrajectoryOutputFile.models.push_back       ( gModel );

        //============================================ Set Gemmi structure internal values based on model information
        gemmi::setup_entities                         ( myTrajectoryOutputFile );
        gemmi::assign_label_seq_id                    ( myTrajectoryOutputFile, true );
        gemmi::assign_subchains                       ( myTrajectoryOutputFile, true );
        
        //============================================ Write out the file
        writeOutCif ( myTrajectoryOutputFile, fileName, myParameterReader.trajectoryFileRemarks );
    }
    
    //============================================ Done
    return ;
    
}

#endif
