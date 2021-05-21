#include "CifOutput.h"
#include "molmodel/internal/CompoundSystem.h"
#include <BiopolymerClass.h>
#include <gemmi/metadata.hpp>
#include <gemmi/model.hpp>

#include <cassert>

#define GEMMI_WRITE_IMPLEMENTATION
#include <gemmi/to_mmcif.hpp>
#include <gemmi/to_cif.hpp>
#include <gemmi/gz.hpp>

static
gemmi::EntityType getEntityType( const std::string& chainID, const SimTK::CIFOut::Data &data )
{
    auto bpIt = data.biopolymers.find(chainID);
    if (bpIt != data.biopolymers.cend())
    {
        const auto bioType   = bpIt->second.getBiopolymerType();
        const bool isPolymer = ( bioType == BiopolymerType::RNA ) || ( bioType == BiopolymerType::DNA ) || ( bioType == BiopolymerType::Protein );
        return isPolymer ? gemmi::EntityType::Polymer : gemmi::EntityType::NonPolymer;
    }

    return gemmi::EntityType::NonPolymer;
}

void SimTK::CIFOut::buildModel ( const State& state, gemmi::Model& gModel, const Data& data, const CompoundSystem& system, int precision )
{
    assert(precision > 0);

    for (SimTK::CompoundSystem::CompoundIndex c(0); c < system.getNumCompounds(); ++c)
    {
        const auto& sysCompound                    = system.getCompound(c);
        const auto  entityType                     = getEntityType(sysCompound.getPdbChainId(), data);

        //======================================== Build the Gemmi model from molmodel data
        sysCompound.buildCif                       ( state, &gModel, entityType, precision, Transform( Vec3 ( 0 ) ) );
    }
}

void SimTK::CIFOut::writeOutCif ( const gemmi::Structure& outStruct, const std::string& fileName, const std::vector < std::pair < std::string, std::string > >& remarks )
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
        
        //============================================ Write out the file
        writeOutCif                                   ( myTrajectoryOutputFile, fileName, myParameterReader.trajectoryFileRemarks );
        
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
}
