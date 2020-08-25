#include "CifOutput.h"

#ifdef GEMMI_USAGE
    #define GEMMI_WRITE_IMPLEMENTATION
    #include <gemmi/to_mmcif.hpp>
    #include <gemmi/to_cif.hpp>

void  SimTK::CIFOut::writeOutCif ( gemmi::Structure outStruct, std::string fileName, std::vector < std::pair < std::string, std::string > > remarks )
{
    //================================================ Create document from the structure with models
    gemmi::cif::Document outDocument                  = gemmi::make_mmcif_document ( outStruct );
    
    //================================================ Add remarks
    for ( unsigned int blIter = 0; blIter < static_cast<unsigned int> ( outDocument.blocks.size() ); blIter++ )
    {
        if ( remarks.size() != 0 )
        {
            gemmi::cif::Loop &remarkLoop              = outDocument.blocks.at(blIter).init_loop ( "_database_PDB_remark.", {"id", "text"} );
            remarkLoop.add_row                        ( { "3", "OTHER REFINEMENT REMARKS:" } );
            
            for ( unsigned int noRem = 0; noRem < static_cast<unsigned int> ( remarks.size() ); noRem++ )
            {
                remarkLoop.add_row                    ( { remarks.at(noRem).first, remarks.at(noRem).second } );
            }
        }
    }
    
    //================================================ Update all blocks
    for ( unsigned int blIter = 0; blIter < static_cast<unsigned int> ( outDocument.blocks.size() ); blIter++ )
    {
        gemmi::update_cif_block                       ( outStruct, outDocument.blocks.at(blIter), true );
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
#endif
