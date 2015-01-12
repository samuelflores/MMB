#ifndef RandomizeRNACoordinates_H_
#define RandomizeRNACoordinates_H_
//using namespace std;
//using namespace SimTK;

static int RandomizeRNACoordinates  (CompoundSystem & system,  State & state,SimbodyMatterSubsystem &  matter , ParameterReader & myParameterReader,  RNA & myChain)

{
    float initialTemperature = myParameterReader.monteCarloTemperature;  //400;
    float temperatureIncrement = myParameterReader.monteCarloTemperatureIncrement;//.02;
    Vec3        newTranslation;
    float myRandx,myRandy,myRandz,myRandA;
    ofstream myoutfile(myParameterReader.outMonteCarloFileName.c_str());
    system.realizeTopology();
    system.realize(state,Stage::Dynamics);
    float oldEnergy = system.calcPotentialEnergy(state);
    float newEnergy = system.calcPotentialEnergy(state);
    State oldState = state;
    srand(time(NULL)*time(NULL));
    //cout<<"[RandomizeRNACoordinates.h] check 1 myRand :"<<myRand<<endl;
    //cout<<"[RandomizeRNACoordinates.h] check 2 myRand :"<<myRand<<endl;
    float temperature = initialTemperature;
    int counter = 0;
    myoutfile << "MODEL "<<temperature<<endl;
    for (SimTK::CompoundSystem::CompoundIndex c(0); c < system.getNumCompounds(); ++c)
        (system.getCompound(c)).writePdb(state, myoutfile   ,Transform(Vec3(0)));

    //myChain.writePdb(state,myoutfile);
    myoutfile << "ENDMDL"<<endl;
    while (temperature>initialTemperature/100) // (newEnergy >= oldEnergy) 
    {
    temperature = exp(-counter/temperatureIncrement)*initialTemperature;
    if (myParameterReader.verbose) cout<<"andomizeRNACoordinates.h] temperature = "<<temperature<<endl;
    counter++;
    //temperature -= temperatureIncrement;
    double myRand = rand();///RAND_MAX;//
    myRand = myRand/RAND_MAX;
    //int i = myRand*myChain.getNumAtoms();//for (int i = 0; i < myChain.getNumAtoms(); i++) 
    //if ((myChain.getAtomName(Compound::AtomIndex(i)).compare("O2*") != 0 ) &&
    //    (i>0) )
    int numMoves = 2;
    int numFixedAtoms = 0;
    //if (myParameterReader.weldToGround)	numFixedAtoms = (myChain.getResidue(ResidueInfo::Index(0))).getNumAtoms();	
    for (int k = 0; k<numMoves; k++)
    {
        int i = myRand*(matter.getNBodies  ());
        //int i = myRand*(myChain.getNumAtoms()-numFixedAtoms) + numFixedAtoms; // make sure we don't get an i = 0.
        //if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<"myChain.getAtomName(Compound::AtomIndex(i)"<< myChain.getAtomName(Compound::AtomIndex(i))<<endl;
        

        MobilizedBodyIndex bodyId = matter.updMobilizedBody(SimTK::MobilizedBodyIndex(i));
        //MobilizedBodyIndex bodyId = matter.updMobilizedBody((myChain.getAtomMobilizedBodyIndex(Compound::AtomIndex(i))));//SimTK::MobilizedBodyIndex(i));

/*
            stringstream ss3;
            ss3<<0<<"/N9" ;
            MobilizedBodyIndex bodyId = matter.updMobilizedBody(myChain.getAtomMobilizedBodyIndex(Compound::AtomIndex(myChain.getAtomIndex(ss3.str()))));
*/
        if (myParameterReader.verbose) cout<<"[RandomizeRNACoordinates.h] bodyId :"<<bodyId<<endl; 
        MobilizedBody & myBody (matter.updMobilizedBody(bodyId));
        Vector myQVector(myBody.getQAsVector(state));
        if (myParameterReader.verbose) cout<<"[RandomizeRNACoordinates.h] myBody.getQ()"<<myBody.getQAsVector(state)<<endl;//ctor.size()<<endl; 
        if (myParameterReader.verbose) cout<<"[RandomizeRNACoordinates.h]  AtomIndex, myQVector,myQVector.size() :"<<i<<","<<myQVector<< ","<<myQVector.size()<<endl; 
        //cout<<"[RandomizeRNACoordinates.h]  AtomName :"<<myChain.getAtomName(Compound::AtomIndex(i))<<endl;//<<","<<myQVector<< ","<<myQVector.size()<<endl; 
        if (myQVector.size  () == 1) {
            //myBody.setQFromVector(state,myQVector);//Vector(rand()*360/Rad2Deg)); // from main
            //myBody.setQFromVector(state,Vector(rand()*360/Rad2Deg)); // from main
            myRand = rand();///RAND_MAX;//
            myRand = myRand/RAND_MAX;
            if (myParameterReader.verbose) cout<<"[RandomizeRNACoordinates.h] rand(),RAND_MAX,myRand :"<<rand()<<","<<RAND_MAX<<","<<myRand<<endl;
            myQVector[0]=myRand*360/Rad2Deg;
            myBody.setQFromVector(state,myQVector);//(1.4439));              // from main
            //myBody.setQFromVector(state,Vector(1.4439));              // from main
            if (myParameterReader.verbose) cout<<"[RandomizeRNACoordinates.h] myBody.getQ(), after setQ"<<myBody.getQAsVector(state)<<endl;//ctor.size()<<endl; 
	    //myBody.setQFromVector(state,Vector(0.0));//rand()*360/Rad2Deg));
	    //myBody.setQFromVector(state,Vector(rand()*360/Rad2Deg));
	} else if (myQVector.size  () == 7) {
            Transform newTransform;
            Transform oldTransform = myBody.getBodyTransform(state)	;
            myRandx = rand();
            myRandx = myRandx/RAND_MAX -.5;
            myRandy = rand();
            myRandy = myRandy/RAND_MAX - .5;
            myRandz = rand();
            myRandz = myRandz/RAND_MAX - .5;
            newTranslation =                   Vec3(myRandx,myRandy,myRandz);
            if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" myRandx, myRandy,myRandz :"<< myRandx<<","<< myRandy<<","<<myRandz   <<endl;
            myRandx = rand();
            myRandx = myRandx/RAND_MAX -.5;
            myRandy = rand();
            myRandy = myRandy/RAND_MAX - .5;
            myRandz = rand();
            myRandz = myRandz/RAND_MAX - .5;
            myRandA = rand();
            myRandA = (myRandA/RAND_MAX - .5)*90 /Rad2Deg;
            if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" myRandx, myRandy,myRandz,myRandA :" << myRandx<<","<< myRandy<<","<<myRandz<<","<<myRandA   <<endl;
            newTransform.set( 
                      Rotation(myRandA,Vec3(myRandx,myRandy,myRandz))*oldTransform.R(), 
                      oldTransform.T() +  newTranslation //Vec3(myRandx,myRandy,myRandz)
                ); 
            if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" body's old Transform:" << myBody.getBodyTransform(state)    <<endl;
            if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" body's old Q:" << myBody.getQAsVector(state)                <<endl;
            myBody.setQToFitTransform(state,newTransform);
            if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" newTransform :"<<newTransform<<endl; 
            if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" setting a FREE body's rotation and translation"<<endl;
            if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" body's new Transform:" << myBody.getBodyTransform(state)    <<endl;
            if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" body's new Q:" << myBody.getQAsVector(state)                <<endl;
        }
    }
    system.realize(state,Stage::Dynamics);
    float newEnergy = system.calcPotentialEnergy(state);
    double myRand2 = rand();///RAND_MAX;//
    assert(myRand2 > 0);
    myRand2 = myRand2/RAND_MAX;
    if (myParameterReader.verbose) cout<< "[RandomizeRNACoordinates.h] oldEnergy ="<<oldEnergy<<endl;
    if (myParameterReader.verbose) cout<< "[RandomizeRNACoordinates.h] newEnergy ="<<newEnergy<<endl;
    if (myParameterReader.verbose) cout <<"[MonteCarlo.cpp] exp(-(newEnergy - oldEnergy)/temperature), myRand2 :"<<(exp(-(newEnergy - oldEnergy)/temperature))<<","<< myRand2<<endl;
    if (!(exp(-(newEnergy - oldEnergy)/temperature) > myRand2)  ) 
    {
    state = oldState;
    if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" rejecting new state"<<endl;
    } 
    else 
    {
        oldState = state;  
        if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" accepting new state"<<endl;
        oldEnergy = newEnergy; 
        myoutfile << "MODEL "<<temperature<<endl;
        for (SimTK::CompoundSystem::CompoundIndex c(0); c < system.getNumCompounds(); ++c)
            (system.getCompound(c)).writePdb(state,myoutfile    ,Transform(Vec3(0)));
        //myChain.writePdb(state,myoutfile);
        myoutfile << "REMARK Potential Energy = "<<setprecision(8)<<newEnergy<<endl;
        //stringstream ss1;
        //ss1.clear();
        //sprintf ("REMARK Potential Energy = %f ", newEnergy);
        myoutfile << "ENDMDL"<<endl;
    }
//if (myParameterReader.verbose) cout<< "[RandomizeRNACoordinates.h] calcEnergy"<<system.calcEnergy(state)<<endl;//  ialEnergy(state);
    }
    myoutfile.close();
	return(0);
}
#endif
