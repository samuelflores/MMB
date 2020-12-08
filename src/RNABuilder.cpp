/* vim: set ts=4 sw=4 sts=4 expandtab */

/* -------------------------------------------------------------------------- *
                          MMB (MacroMoleculeBuilder)                       
   -------------------------------------------------------------------------- */
                                                                           

#include <cstdlib>
#include <fstream>
#include <ios>
#include <iostream>
#include <vector>

#include <getopt.h>
#include <signal.h>

#include "MMBLogger.h"
#include "SimTKmolmodel.h"
 //#include "SimTKsimbody_aux.h"
#include "ParameterReader.h"
#include "Repel.h"
#include "Utils.h"
//#include "RNANoHydrogens.h"
//#include "PeriodicPdbAndEnergyWriter.h"
#include "ProgressWriter.h"
#define _DEBUG_FLAGS_ON_

#define PARAM_DOUT 257
#define PARAM_PROG 256

static
void sig_handler(int signum) {
    if (signum != SIGTERM)
        return;

    GlobalProgressWriter::get().update(ProgressWriter::State::FINISHED);
    GlobalProgressWriter::close();

    std::cerr << "SIGTERM-ing" << signum << std::endl;
    std::quick_exit(EXIT_SUCCESS);
}

static struct option long_opts[] = {
    {"commands",  required_argument, 0,        'C'},
    {"directory", required_argument, 0,        'D'},
    {"HELP",      no_argument,       0,        'H'},
    {"output",    required_argument, 0, PARAM_DOUT},
    {"progress",  required_argument, 0, PARAM_PROG},
    {0,           0,                 0,          0}
};

static
void printUsage() {
    std::cout<<std::endl;
    std::cout << "Usage: MMB [options] \n" << endl;
    std::cout << ".. your MMB executable name will vary depending on platform and release." << endl;
    std::cout << "Options: " << std::endl;
    std::cout << " -H  HELP                 Display this information " << std::endl;
    std::cout << " -C  commands             Set name of contacts file " << std::endl;
    std::cout << " -output                  Name of diagnostic output file (stdout of not specified)" << std::endl;
    std::cout << " -progress                Name of progress report file " << std::endl;
    //std::cout << " -D  directory         Set working directory " << std::endl<<std::endl;
    std::cout << "Last compiled on "<<__DATE__<<" at "<<__TIME__<< std::endl<<std::endl;
    std::cout << "MMB units are nm, kJ/mol, ps, and daltons (g/mol). In MMB 2.10 and earlier, we took some lengths and ground locations in Å (atomSpring, springToGround, atomTether, applyContactsWithin, applyMobilizersWithin, etc.).  As of MMB 2.11 all such lengths and locations are in nm.  Please update your older scripts if you plan to reuse them in MMB 2.11 and later." <<std::endl <<std::endl;
    //std::cout<<"Debug flag : "<<__DEBUG__<<std::endl;
    std::cout<<std::endl;

}

int main(int num_args, char *args[]){  //int argc, char *argv[]) {

    printUsage();

    String option ="";
    String arg ="";
    String parameterFile = "commands.dat";
    String outputDir = "./";
    String progressFile = "";
    String diagOutputFile = "";
    std::ofstream diagOutputStm;


    MMBLOG_FILE_FUNC_LINE(INFO, " Current working directory: "<<Pathname::getCurrentWorkingDirectory()<<endl);

    bool useCurrentDir = true;

    int oc, opt_idx = 0;
    while ((oc = getopt_long_only(num_args, args, "C:D:H", long_opts, &opt_idx)) != -1) {
        switch (oc) {
        case 'C':
            parameterFile = optarg;
            break;
        case 'D':
            outputDir = optarg;
            break;
        case 'H':
            printUsage();
            return EXIT_SUCCESS;
        case PARAM_DOUT:
            diagOutputFile = optarg;
            break;
        case PARAM_PROG:
            progressFile = optarg;
            break;
        case ':':
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Parameter -" << char(optopt) << " requires an argument" << std::endl);
            return EXIT_FAILURE;
        case '?':
            if (optopt == 0)
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unknown parameter " << args[optind - 1] << std::endl);
            else
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "Parameter -" << char(optopt) << " requires an argument" << std::endl);
            return EXIT_FAILURE;
        default:
            printUsage();
            return EXIT_FAILURE;
       }
    }

    /* determine the directory that the RNABuilder executable is in */
    if( useCurrentDir ) {
        String dir = args[0];
        int slashPos = dir.find_last_of("/");
        if( slashPos < 0 ) {
            outputDir = "./";
        } else {
            outputDir  = dir.substr(0,slashPos+1 );
        }
    }

    cout << " output directory = " << outputDir << endl;

    GlobalProgressWriter::initialize(progressFile);
    if (signal(SIGTERM, sig_handler) == SIG_ERR)
	    MMBLOG_PLAIN(CRITICAL, "Failed to install SIGTERM handler");

    try 
    {

        //int stdReportingIntervals = 20000;
        //int stdReportingIntervals =  50;
        //int startHere = 0;
        //int firsti = 1;
        //int lasti  = 3;
        stringstream ss2b;
        stringstream ss3;
        map<const char*, int, strCmp> firstResidueNumbers;

        if (diagOutputFile != "") {
            diagOutputStm.open(diagOutputFile, std::ios_base::out | std::ios_base::trunc);
            if (!diagOutputStm.is_open())
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "Cannot open diagnostic output file");
            MMBLogger::instance().setOutput(&diagOutputStm);
        }

        ParameterReader myParameterReader;
        
        MMBLOG_FILE_FUNC_LINE(INFO, "About to read "<<parameterFile<<" to set firstStage and lastStage"<<endl);
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        myParameterReader.setFirstAndLastStageAndUseCifFiles(parameterFile.c_str());

        MMBLOG_FILE_FUNC_LINE(INFO, "lastStage = "<<myParameterReader.lastStage<<endl);
        MMBLOG_FILE_FUNC_LINE(INFO, "firstStage = "<<myParameterReader.firstStage<<endl);
        MMBLOG_FILE_FUNC_LINE(INFO, "useCIFFileFormat = "<<myParameterReader.useCIFFileFormat<<endl);
        if ((myParameterReader.lastStage < myParameterReader.firstStage))   
        {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "stage < 1 error!  Most likely you have failed to specify the command file, (currently "<< parameterFile<<"), or it was not found"<<endl);
        }

        for (int i = myParameterReader.firstStage; i<=  myParameterReader.lastStage; i++) {
            myParameterReader.initializeDefaults();

            myParameterReader.currentStage = i;
            if (myParameterReader.currentStage<1) {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "stage < 1 error!  "<<endl); //Most likely you have failed to specify the command file, (currently "<< parameterFile<<"), or it was not found"<<endl;
            }

            if  (myParameterReader.lastFrameFileName == "NOT-SET"){
                stringstream ss2;
                ss2.clear();
                ss2.str("");
                //ss2<<"./last."<<i<<".cif"; 
                ss2<<"./last."<<i<<".pdb"; 
                myParameterReader.lastFrameFileName = ss2.str();
            } 
            if  (myParameterReader.previousFrameFileName == "NOT-SET"){
                stringstream ss4;
                ss4.clear();
                ss4.str("");
                ss4<<"./last."<<(i-1)<<".pdb"; 
                //ss4<<"./last."<<(i-1)<<".cif"; 
                myParameterReader.previousFrameFileName = ss4.str();
            }
            // SCF moved from here
            //printBiopolymerSequenceInfo(myParameterReader.myBiopolymerClassContainer.updBiopolymerClass("g").myBiopolymer);
            // end new way
            
            std::stringstream ss2;
            ss2.clear(); ss2.str("");
            ss2 << "./last." << i << ".pdb";
            MMBLOG_FILE_FUNC_LINE(INFO, "myParameterReader.lastFrameFileName = " << myParameterReader.lastFrameFileName << " , myParameterReader.useCIFFileFormat = "<<myParameterReader.useCIFFileFormat <<endl);
            if  (myParameterReader.lastFrameFileName == ss2.str())
            {
                if ( myParameterReader.useCIFFileFormat )
                {
                    ss2.clear(); ss2.str("");
                    ss2 << "./last." << i << ".cif";
                    myParameterReader.lastFrameFileName = ss2.str();
                }
                else
                {
                    //================================ Do nothing, PDB already set
                }
            }
            
            std::stringstream ss4; ss4.clear(); ss4.str("");
            ss4 << "./last." << ( i - 1 ) << ".pdb";
            if  (myParameterReader.previousFrameFileName == ss4.str())
            {
                if ( myParameterReader.useCIFFileFormat )
                {
                    ss4.clear(); ss4.str("");
                    ss4 << "./last." << (i-1) << ".cif";
                    myParameterReader.previousFrameFileName = ss4.str();
                }
                else
                {
                    //================================ Do nothing, PDB already set
                }
            }
	    // SCF moved to here
            MMBLOG_FILE_FUNC_LINE(INFO, "about to make sure we aren't trying to read from BOTH .pdb and QVector files, but exactly one."<<endl);
            if (i == 1)  {myParameterReader.readPreviousFrameFile = 0;} //else myParameterReader.readPreviousFrameFile = 1;
            if (i == 1)  {myParameterReader.readInQVector         = 0;} //else myParameterReader.readPreviousFrameFile = 1;

            myParameterReader.initializeFromFileOnly(parameterFile.c_str());
            //
            myParameterReader.postInitialize();
            //scf added .. there was an issue with counting Coulomb forces
            myParameterReader.removeNonPriorityBasePairs(myParameterReader.currentStage);//i);
            MMBLOG_FILE_FUNC_LINE(INFO, "myParameterReader.basePairContainer.numBasePairs()"<<  myParameterReader.basePairContainer.numBasePairs()<<endl);

            MMBLOG_FILE_FUNC_LINE(INFO, "myParameterReader.basePairContainer.numBasePairs()"<<  myParameterReader.basePairContainer.numBasePairs()<<endl);
            
            ConstrainedDynamics  myConstrainedDynamics(&myParameterReader);
            myConstrainedDynamics.initializeDumm();

            MMBLOG_FILE_FUNC_LINE(INFO, "myParameterReader.firstStage="<<myParameterReader.firstStage<<endl);
            ss2b.clear();
            ss2b.str("");
            ss2b<<"./last.qVector."<<i<<".dat"; 
            myParameterReader.outQVectorFileName = (ss2b.str());
            ss2b.clear();
            ss2b.str("");
            ss2b<<"./last.qVector."<<i-1<<".dat"; 
            myParameterReader.inQVectorFileName = (ss2b.str());
            //myConstrainedDynamics.readInQVector =0;

            if ( myParameterReader.useCIFFileFormat )
            {
                //==================================== Using mmCIF format
                ss3.clear();
                ss3.str("");
                ss3<<"./trajectory."<<i<<".cif";
                myParameterReader.outTrajectoryFileName = ss3.str();
            }
            else
            {
                //==================================== Using PDB format
                ss3.clear();
                ss3.str("");
                ss3<<"./trajectory."<<i<<".pdb";
                myParameterReader.outTrajectoryFileName = ss3.str();
            }
            
            stringstream ss6;
            ss6.clear();
            ss6.str("");
            ss6<<"./trajectory.MonteCarlo."<<i<<".pdb";

            ofstream output;
            if ( myParameterReader.useCIFFileFormat )
            {
#ifdef GEMMI_USAGE
                myParameterReader.lastFileRemarks.push_back ( std::pair < std::string, std::string > ( "3", "Contents of user input file " + std::string( parameterFile ) + " :" ) );
                ifstream inFile                       = ifstream ( parameterFile.c_str(),ios_base::in );
                while ( inFile.good() )
                {
                    String tempString;
                    getline                           ( inFile, tempString );
                    myParameterReader.lastFileRemarks.push_back ( std::pair < std::string, std::string > ( "3", tempString ) );
                }
                myParameterReader.lastFileRemarks.push_back ( std::pair < std::string, std::string > ( "3", "End of user input file." ) );
#endif
            }
            else
            {
                output                                = ofstream ( myParameterReader.outTrajectoryFileName.c_str() );
                ifstream inFile                       = ifstream ( parameterFile.c_str(),ios_base::in );
                output << endl << "REMARK contents of user input file " << parameterFile << " : " << endl << endl;
                while ( inFile.good() )
                {
                    String tempString;
                    getline                           ( inFile, tempString );
                    output << String ( "REMARK " ) + tempString << endl;
                }
                output << endl << "REMARK end of user input file " << endl << endl;
            }


            myParameterReader.outMonteCarloFileName = ss6.str();

            MMBLOG_FILE_FUNC_LINE(INFO, endl);
            myParameterReader.removeNonPriorityBasePairs(myParameterReader.currentStage);//i);
            MMBLOG_FILE_FUNC_LINE(INFO, "BASE PAIRS before removing from rigid stretches:"<<endl);
            myParameterReader.basePairContainer.printBasePairs();
            if (myParameterReader.setRemoveBasePairsInRigidStretch) myParameterReader.removeBasePairsInRigidStretch();
            MMBLOG_FILE_FUNC_LINE(INFO, "  1"<<endl);
            MMBLOG_FILE_FUNC_LINE(INFO, "  stage ="<<i<<endl);
            if (myParameterReader.currentStage<1) {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "stage < 1 error!  Most likely you have failed to specify the command file, (currently "<< parameterFile<<"), or it was not found"<<endl);
            }
            myParameterReader.printAllSettings(std::cout,String("") );
            if ( myParameterReader.useCIFFileFormat )
            {
#ifdef GEMMI_USAGE
                myParameterReader.lastFileRemarks.push_back ( std::pair < std::string, std::string > ( "3", "About to call myParameterReader.printAllSettings." ) );
                myParameterReader.printAllSettingsToMMCIF ( myParameterReader.lastFileRemarks );
                myParameterReader.lastFileRemarks.push_back ( std::pair < std::string, std::string > ( "3", "Done with call to myParameterReader.printAllSettings." ) );
#endif
            }
            else
            {
                output<<endl<<"REMARK "<< __FILE__<<":"<<__LINE__<<  " about to call myParameterReader.printAllSettings "<<endl<<endl;
                myParameterReader.printAllSettings(output ,String("REMARK ")   );
                output<<endl<<"REMARK "<< __FILE__<<":"<<__LINE__<<  " done with call to myParameterReader.printAllSettings "<<endl<<endl;
            }

            //printBiopolymerSequenceInfo(myParameterReader.myBiopolymerClassContainer);
            MMBLOG_FILE_FUNC_LINE(INFO, endl);
            myConstrainedDynamics.runDynamics();
            myParameterReader.atomSpringContainer.calcKabschRmsd(myConstrainedDynamics.getCurrentState(),myParameterReader.myBiopolymerClassContainer);
            MMBLOG_FILE_FUNC_LINE(INFO, endl);
            closingMessage();
            MMBLOG_FILE_FUNC_LINE(INFO, endl);
            //exit(0); //hoping to avoid the corrupted double linked list issue
        }

    }
    catch (const std::exception & e) {
        MMBLOG_FILE_FUNC_LINE(INFO, e.what() <<endl); //FIXME: This will log the error twice!
    }

}

