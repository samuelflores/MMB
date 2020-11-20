/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Molmodel                            *   
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *   
 * Simbios, the NIH National Center for Physics-Based Simulation of           *   
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *   
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *   
 *                                                                            *   
 * Portions copyright (c) 2008 Stanford University and the Authors.           *   
 * Authors: Samuel Flores                                                     *   
 * Contributors: Christopher Bruns, Peter Eastman                             *   
 *                                                                            *   
 * Permission is hereby granted, free of charge, to any person obtaining a    *   
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *   
 * and/or sell copies of the Software, and to permit persons to whom the      *   
 * Software is furnished to do so, subject to the following conditions:       *   
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

/* This is based on WaterDroplet.h, but instead of generating a droplet of    *
 * water, it generates a "droplet" of P12 ligands.  Peter recommends using    *
*  a factory class -- I leave that to future generations to consider.         *
*/


#include "molmodel/internal/common.h"
#include "molmodel/internal/Compound.h"
#include "MMBLogger.h"
#include <iostream>
#include <fstream>
//#include "/Users/samuelflores/rna-dynamics/Water.h"
//#include <ios>

namespace SimTK 
{

class LigandDroplet : public Compound { public:
        LigandDroplet(CompoundSystem &system 
		, TinkerDuMMForceFieldSubsystem &dumm
		, GeneralForceSubsystem& forces
		, Vec3 dropletCenter = Vec3(0,0,0 )
		, float dropletRadius = 2.50  //in nm, of course
		, float arista = 1.6  // this is enough to separate out the ligands sufficiently     
		)
        {   
	Vec3 atomLocation;
	AtomPathName myAtomName;
	float shellRadius =10.0; //nanometers
	int cellsAcross = dropletRadius*2/arista+1; // number of cells across.  add 1 b.c. water molecules are at corners of cells.
	std::cout<<"cellsAcross = "<<cellsAcross<<std::endl;
	int totalAtoms =0;
	int      myWaterVecOccupation[cellsAcross][cellsAcross][cellsAcross];
	int numWater =0;//lsAcross*cellsAcross*cellsAcross;

        VanderWallSphere myVanderWallSphere( forces,  dumm,dropletCenter,shellRadius, .1            , 1.0           );

	float myRadius=0.;
        for (int xi = 0; xi < cellsAcross; xi++)
        	for (int yi = 0; yi < cellsAcross; yi++)
        		for (int zi = 0; zi < cellsAcross; zi++)
                	{
				cout<<"for indices "<<xi<<","<<yi<<","<<zi<<endl;
				
				cout<<" checking condition: "<<dropletCenter[0]<<","<<xi<<" , = "<<pow((dropletCenter[0]-(xi*arista-dropletRadius)),2)<<endl;
				myRadius = pow((pow((xi*arista-dropletRadius),float(2))+ pow((yi*arista-dropletRadius),float(2))+ pow((zi*arista-dropletRadius),float(2))),float(.5) );
				cout<<"radius = "<<myRadius<<endl;
				if (myRadius <= dropletRadius) //(pow((pow((dropletCenter[0]-(xi*arista+dropletRadius)),2)+ pow((dropletCenter[1]-(yi*arista+dropletRadius)),2)+ pow((dropletCenter[2]-(zi*arista+dropletRadius)),2)),.5 ) <= dropletRadius)  
					{
						
					numWater++;myWaterVecOccupation[xi][yi][zi] = 1;
					}
				else
		  			{cout<<"outside of water droplet radius; no water to be placed here."<<endl; myWaterVecOccupation[xi][yi][zi] = 0;}       
			}	
	//Count the atoms in the CompoundSystem
	
        for (int i =0 ; i< system.getNumCompounds(); i++)
	{
		cout<<"counting atoms in compound # "<<i<<endl;	
		totalAtoms += system.updCompound(Compound::Index(i)).getNAtoms();
	}	
	cout<<"counted "<<totalAtoms<<" total atoms."<<endl;
	//Create an array to hold their locations
	Vec3 atomLocationArray[totalAtoms];
	int myAtomIndex = 0;
	//iState state = system.realizeTopology();
	for (int i =0 ; i< system.getNumCompounds(); i++)
	{
		for (int j = 0; j< system.updCompound(Compound::Index(i)).getNAtoms(); j++)
		{
			Compound & myCompound = system.updCompound(Compound::Index(i));
			myAtomName = myCompound.getAtomName(Compound::AtomIndex(j));
			atomLocationArray[myAtomIndex] = myCompound.calcDefaultAtomLocationInGroundFrame(myAtomName);
			//atomLocationArray[myAtomIndex] = myCompound.getDefaultAtomLocation(myAtomName);
			//atomLocationArray[myAtomIndex] = myCompound.getTopLevelTransform() * atomLocationArray[myAtomIndex] ;  
			//atomLocationArray[myAtomIndex] = system.updCompound(Compound::Index(i)).calcAtomLocationInGroundFrame(state,Compound::AtomIndex(j));
			cout << "myAtomName= "<<myAtomName<<" atomLocation= "<<atomLocationArray[myAtomIndex]<<endl;
			int myxi = (atomLocationArray[myAtomIndex][0]-dropletCenter[0]+dropletRadius)/arista+.5;
			int myyi = (atomLocationArray[myAtomIndex][1]-dropletCenter[1]+dropletRadius)/arista+.5;
			int myzi = (atomLocationArray[myAtomIndex][2]-dropletCenter[2]+dropletRadius)/arista+.5;
			cout<<" indices = "<<myxi<<","<<myyi<<","<<myzi<<endl;
			if (((myxi < cellsAcross) && (myxi >=0))
				&& ((myyi < cellsAcross) && (myyi >=0))
				&& ((myzi < cellsAcross) && (myzi >=0))
				&& (myWaterVecOccupation[myxi][myyi][myzi] == 1))
			{
					numWater--;
					cout<<"setting an occupation to 0"<<endl;					
					myWaterVecOccupation[myxi][myyi][myzi] = 0;			
			}
			myAtomIndex++;
		}
	}

	cout<<"about to create array of "<<numWater<<" water molecules."<<endl;
	P12   * myWaterVec[numWater];
        for (int k =0;k<numWater;k++)	{
		myWaterVec[k]=new   P12(dumm);
		(*myWaterVec[k]).setPdbResidueName("P12");
		(*myWaterVec[k]).setPdbResidueNumber(k+1  );
		(*myWaterVec[k]).setPdbChainId('D');
	}

	int myWaterIndex =-1;
        for (int xi = 0; xi < cellsAcross; xi++)
        	for (int yi = 0; yi < cellsAcross; yi++)
        		for (int zi = 0; zi < cellsAcross; zi++)
			{
				cout<<"checking occupation for "<<xi<<","<<yi<<","<<zi<<endl;
				if (myWaterVecOccupation[xi][yi][zi]) 
				{
					myWaterIndex++;
					if (myWaterIndex > numWater) {
						MMBBLOG_FILE_FUNC_LINE(CRITICAL, "Attempted to adopt more water than exists in array"<<endl);
					}
				  	cout<<"adopting water at "<< xi*arista-dropletRadius+dropletCenter[0]<<","<<yi*arista-dropletRadius+dropletCenter[1]<<"=,"<<zi*arista-dropletRadius+dropletCenter[2]<<endl;
					Rotation myRotation(45.*Deg2Rad, YAxis);
					myRotation =  myRotation * Rotation(45.*Deg2Rad, XAxis);
					Transform myTransform(myRotation,Vec3(xi*arista-dropletRadius+dropletCenter[0],yi*arista-dropletRadius+dropletCenter[1],zi*arista-dropletRadius+dropletCenter[2]));
				        system.adoptCompound(*myWaterVec[myWaterIndex],myTransform);//Vec3(xi*arista-dropletRadius+dropletCenter[0],yi*arista-dropletRadius+dropletCenter[1],zi*arista-dropletRadius+dropletCenter[2]));
				}	
			}
		
		
	

	}
};
}
