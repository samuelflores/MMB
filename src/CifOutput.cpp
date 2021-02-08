#include "CifOutput.h"

#ifdef GEMMI_USAGE
    #define GEMMI_WRITE_IMPLEMENTATION
    #include <gemmi/to_mmcif.hpp>
    #include <gemmi/to_cif.hpp>
    #include <gemmi/gz.hpp>

void SimTK::CIFOut::writeOutCif ( const gemmi::Structure &outStruct, const std::string &fileName, const std::vector < std::pair < std::string, std::string > > &remarks )
{
    //================================================ Create document from the structure with models
    gemmi::cif::Document outDocument                  = gemmi::make_mmcif_document ( outStruct );

    //================================================ Add remarks
    for ( decltype(outDocument.blocks)::size_type blIter = 0; blIter < outDocument.blocks.size(); blIter++ )
    {
        if ( remarks.size() != 0 )
        {
            //======================================== Generate single string with the remark
            std::string remText                       = ( "\n; " );
            for ( size_t noRem = 0; noRem < remarks.size(); noRem++ )
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
    for ( decltype(outDocument.blocks)::size_type blIter = 0; blIter < outDocument.blocks.size(); blIter++ )
    {
        gemmi::update_mmcif_block                       ( outStruct, outDocument.blocks.at(blIter), gemmi::MmcifOutputGroups( true ) );
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

void SimTK::CIFOut::reWriteOutCif ( const gemmi::Model& gModel, const std::string& modelName, const std::string& fileName, ParameterReader& myParameterReader, const CompoundSystem& system, bool firstInStage )
{
    //================================================ Try to read in the file
    if ( firstInStage )
    {
        //============================================ File does not exist, create it.
        gemmi::Structure myTrajectoryOutputFile;
        
        //============================================ Add model name
        myTrajectoryOutputFile.name                   = modelName;
        
        //============================================ Save model to structure
        myTrajectoryOutputFile.models.emplace_back    ( gModel );

        //============================================ Set Gemmi structure internal values based on model information
        gemmi::setup_entities                         ( myTrajectoryOutputFile );
        gemmi::assign_label_seq_id                    ( myTrajectoryOutputFile, true );
        gemmi::assign_subchains                       ( myTrajectoryOutputFile, true );
        
        //============================================ Add sequence to entities
        int compoundNumber                            = 1;

        auto biopolymerMap                            = myParameterReader.myBiopolymerClassContainer.getBiopolymerClassMap();
        for ( SimTK::CompoundSystem::CompoundIndex c(0); c < system.getNumCompounds(); ++c )
        {
            //======================================== Print log
            std::cout <<__FILE__<<":"<<__LINE__<<" c = "<<c<< " compoundNumber = "<<compoundNumber <<std::endl;
            
            //======================================== Get sequence for compound
            const auto chainID                        = myParameterReader.myBiopolymerClassContainer.updBiopolymerClass(compoundNumber-1).getChainID();
            auto& compound                            = biopolymerMap.at(chainID);
            const auto residues                       = compound.getResidueIdVector();

            //======================================== For each structure entity
            for ( decltype(myTrajectoryOutputFile.entities)::size_type enIt = 0; enIt < myTrajectoryOutputFile.entities.size(); enIt++ )
            {
                auto& entity = myTrajectoryOutputFile.entities[enIt];
                //==================================== If entity name = MMB chain ID
                if ( entity.name == chainID )
                {
                    //================================ Copy sequence
                    entity.full_sequence.resize( residues.size() );
                    for ( size_t idx = 0; idx < residues.size(); idx++ )
                    {
                        entity.full_sequence[idx] = compound.getPdbResidueName(residues[idx]);
                    }
                    break;
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
