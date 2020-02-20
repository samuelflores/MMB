
/* -------------------------------------------------------------------------- *
                          MMB (MacroMoleculeBuilder)                       
   -------------------------------------------------------------------------- 
                                                                           
Copyright (c) 2011-16 by the author                                        
Author: Samuel Flores                                                      
                                                                           
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to use the Software for non-profit academic research only. Any commercial, for-profit, or industrial use requires the payment of usage fees.                  

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.            
                            
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                         
-------------------------------------------------------------------------- */


/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the author                                        *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to use the Software for non-profit academic research only. Any commercial, *
 * for-profit, or industrial use requires the payment of usage fees.        *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include <fstream>
#include <ios>
#include <iostream>
#include <vector>

 #include "SimTKmolmodel.h"
 //#include "SimTKsimbody_aux.h"
#include "ErrorManager.h"
#include "ParameterReader.h"
#include "Repel.h"
#include "Utils.h"
//#include "RNANoHydrogens.h"
//#include "PeriodicPdbAndEnergyWriter.h"
#define _DEBUG_FLAGS_ON_


void printUsage() {
    std::cout<<std::endl;
    std::cout << "Usage: MMB [options] \n" << endl;
    std::cout << ".. your MMB executable name will vary depending on platform and release." << endl;
    std::cout << "Options: " << std::endl;
    std::cout << " -help                 Display this information " << std::endl;
    std::cout << " -c  contactsFile      Set name of contacts file " << std::endl;
    //std::cout << " -d  directory         Set working directory " << std::endl<<std::endl;
    std::cout << "Commercial use of this software requires payment of usage fees " << std::endl<<std::endl;
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


    std::cout << __FILE__<<":"<<__LINE__<<" Current working directory: "<<Pathname::getCurrentWorkingDirectory()<<std::endl;

    bool useCurrentDir = true;

    for( int i=1;i<num_args;i++) {

        arg = String( args[i]);
        option.resize(arg.length());

        for(int j=0;j<(int)arg.length();j++) {
            option[j] = toupper(arg[j]);
        }

        if( (option == "-HELP" )  || (option == "-H") ) {
            printUsage();
            // user specified name of input file as a command line argument 
        } else if ( option == "-C")   {
            if( num_args < i+2 ) {
                std::cout << "Error missing name of contacts file "  << std::endl;
            } else  {
                parameterFile = args[++i];
            }
            // user specified working directory as command line argument 
        } else if ( option == "-D")  {
            if( num_args < i+2 ) {
                std::cout << "Error missing name of output directory "  << std::endl;
            } else  {
                outputDir = args[++i];
                useCurrentDir = false;
            }

        } else {
            std::cout << "Error Unrecognized option:" << args[i] << std::endl;
            printUsage();
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

        ParameterReader myParameterReader;
        
        cout<<__FILE__<<":"<<__LINE__<<" About to read "<<parameterFile<<" to set firstStage and lastStage"<<endl; 
        cout<<__FILE__<<":"<<__LINE__<<endl; 
        myParameterReader.setFirstAndLastStage(parameterFile.c_str());

        cout<<__FILE__<<":"<<__LINE__<<"  lastStage = "<<myParameterReader.lastStage<<endl; 
        cout<<__FILE__<<":"<<__LINE__<<"  firstStage = "<<myParameterReader.firstStage<<endl; 
        if ((myParameterReader.lastStage < myParameterReader.firstStage))   
        {
            ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" stage < 1 error!  Most likely you have failed to specify the command file, (currently "<< parameterFile<<"), or it was not found"<<endl;
            ErrorManager::instance.treatError();
        }

        for (int i = myParameterReader.firstStage; i<=  myParameterReader.lastStage; i++) {
            myParameterReader.initializeDefaults();

            myParameterReader.currentStage = i;
            if   (myParameterReader.currentStage<1) {
                ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" stage < 1 error!  "<<endl; //Most likely you have failed to specify the command file, (currently "<< parameterFile<<"), or it was not found"<<endl;
                ErrorManager::instance.treatError();
            }

            if  (myParameterReader.lastFrameFileName == "NOT-SET"){
                stringstream ss2;
                ss2.clear();
                ss2.str("");
                ss2<<"./last."<<i<<".pdb"; 
                myParameterReader.lastFrameFileName = ss2.str();
            } 
            if  (myParameterReader.previousFrameFileName == "NOT-SET"){
                stringstream ss4;
                ss4.clear();
                ss4.str("");
                ss4<<"./last."<<(i-1)<<".pdb"; 
                myParameterReader.previousFrameFileName = ss4.str();
            }

            cout<<__FILE__<<":"<<__LINE__<<"  about to make sure we aren't trying to read from BOTH .pdb and QVector files, but exactly one."<<endl;
            if (i == 1)  {myParameterReader.readPreviousFrameFile = 0;} //else myParameterReader.readPreviousFrameFile = 1;
            if (i == 1)  {myParameterReader.readInQVector         = 0;} //else myParameterReader.readPreviousFrameFile = 1;

            myParameterReader.initializeFromFileOnly(parameterFile.c_str());
            //printBiopolymerSequenceInfo(myParameterReader.myBiopolymerClassContainer.updBiopolymerClass("g").myBiopolymer);
            // end new way
            myParameterReader.postInitialize();
            //scf added .. there was an issue with counting Coulomb forces
            myParameterReader.removeNonPriorityBasePairs(myParameterReader.currentStage);//i);
            cout<<__FILE__<<":"<<__LINE__<<"  myParameterReader.basePairContainer.numBasePairs()"<<  myParameterReader.basePairContainer.numBasePairs()<<        endl;

            cout<<__FILE__<<":"<<__LINE__<<"  myParameterReader.basePairContainer.numBasePairs()"<<  myParameterReader.basePairContainer.numBasePairs()<<        endl;
            
            ConstrainedDynamics  myConstrainedDynamics(&myParameterReader);
            myConstrainedDynamics.initializeDumm();

            cout<<__FILE__<<":"<<__LINE__<<"  myParameterReader.firstStage="<<myParameterReader.firstStage<<endl;
            ss2b.clear();
            ss2b.str("");
            ss2b<<"./last.qVector."<<i<<".dat"; 
            myParameterReader.outQVectorFileName = (ss2b.str());
            ss2b.clear();
            ss2b.str("");
            ss2b<<"./last.qVector."<<i-1<<".dat"; 
            myParameterReader.inQVectorFileName = (ss2b.str());
            //myConstrainedDynamics.readInQVector =0;
            ss3.clear();
            ss3.str("");
            ss3<<"./trajectory."<<i<<".pdb"; 
            stringstream ss6;
            ss6.clear();
            ss6.str("");
            ss6<<"./trajectory.MonteCarlo."<<i<<".pdb"; 
            myParameterReader.outTrajectoryFileName = ss3.str();

            ofstream output(myParameterReader.outTrajectoryFileName.c_str());
            ifstream inFile(parameterFile.c_str(),ios_base::in);
            output<<endl<<"REMARK contents of user input file "<<parameterFile<<" : "<<endl<<endl;
            while (inFile.good()) {
                String tempString;
                getline(inFile, tempString);
                output<<String("REMARK ") + tempString << endl;
            }  
            output<<endl<<"REMARK end of user input file "<<endl<<endl;


            myParameterReader.outMonteCarloFileName = ss6.str();

            cout <<__FILE__<<":"<<__LINE__<<" "<<endl;
            myParameterReader.removeNonPriorityBasePairs(myParameterReader.currentStage);//i);
            cout<<__FILE__<<":"<<__LINE__<<"  BASE PAIRS before removing from rigid stretches:"<<endl;
            myParameterReader.basePairContainer.printBasePairs();
            if (myParameterReader.setRemoveBasePairsInRigidStretch) myParameterReader.removeBasePairsInRigidStretch();
            cout <<__FILE__<<":"<<__LINE__<<"  1"<<endl;
            cout<<__FILE__<<":"<<__LINE__<<"  stage ="<<i<<endl;
            if (myParameterReader.currentStage<1) {
                ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" stage < 1 error!  Most likely you have failed to specify the command file, (currently "<< parameterFile<<"), or it was not found"<<endl;
                ErrorManager::instance.treatError();
            }
            myParameterReader.printAllSettings(std::cout,String("") );
            output<<endl<<"REMARK "<< __FILE__<<":"<<__LINE__<<  " about to call myParameterReader.printAllSettings "<<endl<<endl;
            myParameterReader.printAllSettings(output ,String("REMARK ")   );
            output<<endl<<"REMARK "<< __FILE__<<":"<<__LINE__<<  " done with call to myParameterReader.printAllSettings "<<endl<<endl;

            //printBiopolymerSequenceInfo(myParameterReader.myBiopolymerClassContainer);
            cout <<__FILE__<<":"<<__LINE__<<endl;
            myConstrainedDynamics.runDynamics();
            myParameterReader.atomSpringContainer.calcKabschRmsd(myConstrainedDynamics.getCurrentState(),myParameterReader.myBiopolymerClassContainer);
            cout <<__FILE__<<":"<<__LINE__<<endl;
            closingMessage();
            cout <<__FILE__<<":"<<__LINE__<<endl;
            //exit(0); //hoping to avoid the corrupted double linked list issue
        }

    }
    catch (const std::exception & e) {
        std::cout<<__FILE__<<":"<<__LINE__<<e.what() <<endl;
    }

}

