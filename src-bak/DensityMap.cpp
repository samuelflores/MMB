/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#include "Utils.h"
//#include "float.h"
#include <fstream>
#include "DensityMap.h"
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

using namespace std;
using namespace SimTK;

/*
GridPoint::GridPoint() {}
void GridPoint::initializeGradient(){
		ddyPositiveXGradient = 0;
		ddzPositiveXGradient = 0;
		ddxPositiveYGradient = 0;
		ddzPositiveYGradient = 0;
		ddxPositiveZGradient = 0;
		ddyPositiveZGradient = 0;
                firstQuadrantGradient = fVec3(0,0,0);
	}

void GridPoint::initialize(){
		initializeGradient();
		density = 0; 
		position = fVec3(0);	
	}

void GridPoint::validatePosition(fVec3 myPosition) const{
                //ValidateVec3(myPosition);
	}

void GridPoint::validateDensity(double myDensity) const{
		ValidateReal(myDensity);
	
}

void GridPoint::validate() const{
		validatePosition(position);
		validateDensity(density);
		// write this code later.
	}

void GridPoint::setDensity(Real myDensity)	{validateDensity(myDensity); density = (float)myDensity;	}

void GridPoint::setPosition(fVec3 myPosition)	{
 		validatePosition(myPosition);
		position = myPosition;	
	}

Quadrant GridPoint::calcQuadrant(fVec3 queryPosition) const			{
		Quadrant tempQuadrant;    
		if (queryPosition[0] < position[0]) {tempQuadrant.positiveX = false;} else { tempQuadrant.positiveX = true;}		
		if (queryPosition[1] < position[1]) {tempQuadrant.positiveY = false;} else { tempQuadrant.positiveY = true;}
		if (queryPosition[2] < position[2]) {tempQuadrant.positiveZ = false;} else { tempQuadrant.positiveZ = true;}
		return tempQuadrant;
	}

fVec3 GridPoint::fetchGradient(fVec3 queryPosition) const {
		Quadrant tempQuadrant = calcQuadrant(queryPosition);
		fVec3 myGradient;
          	ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" this function is obsolete, fetchGradient doesn't work anymore"<<endl;
        	ErrorManager::instance.treatError();
		return myGradient;

	}

fVec3 GridPoint::fetchFirstQuadrantGradient() const {
                return firstQuadrantGradient;
}

double GridPoint::getDensity() const	{return density ;}
	
double GridPoint::getDensity( fVec3 queryPosition) const	{
                fVec3 myVectorToGridPoint = queryPosition - position;
                fVec3 myGradient = calcInterpolatedFirstQuadrantGradient(  queryPosition); 
                return (density + (myGradient[0]* myVectorToGridPoint[0])  + (myGradient[1]* myVectorToGridPoint[1]) + (myGradient[2]* myVectorToGridPoint[2]));// DotProduct(myGradient, myVectorToGridPoint));
}
void GridPoint::setPositiveXGradient(Real myPositiveXGradient) { firstQuadrantGradient[0] = myPositiveXGradient; }
void GridPoint::setddyPositiveXGradient(Real value) { ddyPositiveXGradient = (float)value; }
void GridPoint::setddzPositiveXGradient(Real value) { ddzPositiveXGradient = (float)value; }
void GridPoint::setPositiveYGradient(Real myPositiveYGradient) { firstQuadrantGradient[1] = myPositiveYGradient; }
void GridPoint::setddxPositiveYGradient(Real value) { ddxPositiveYGradient = (float)value; }
void GridPoint::setddzPositiveYGradient(Real value) { ddzPositiveYGradient = (float)value; }
void GridPoint::setPositiveZGradient(Real myPositiveZGradient) { firstQuadrantGradient[2] = myPositiveZGradient; }
void GridPoint::setddxPositiveZGradient(Real value) { ddxPositiveZGradient = (float)value; }
void GridPoint::setddyPositiveZGradient(Real value) { ddyPositiveZGradient = (float)value; }
void GridPoint::setNegativeXGradient(Real myNegativeXGradient) {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" this function is only for calculating the gradient in the first quadrant"<<endl;
        ErrorManager::instance.treatError();}
void GridPoint::setNegativeYGradient(Real myNegativeYGradient) { 
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" this function is only for calculating the gradient in the first quadrant"<<endl;
        ErrorManager::instance.treatError();}
void GridPoint::setNegativeZGradient(Real myNegativeZGradient) { 
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" this function is only for calculating the gradient in the first quadrant"<<endl;
        ErrorManager::instance.treatError();}

fVec3 GridPoint::calcInterpolatedFirstQuadrantGradient(fVec3 queryPosition) const {
    fVec3 dxdydz = queryPosition - position; // the first term is the query position, the second term is the grid point position in cartesian space
    if ((dxdydz[0] < 0) || (dxdydz[1] <0 ) || (dxdydz[2] < 0)) { // see if we can make this trap unnecessary implicitly

        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" this function is only for calculating the gradient in the first quadrant"<<endl;
        ErrorManager::instance.treatError();
    }
    fVec3 myGradient;
    myGradient[0] = firstQuadrantGradient[0]  + ddyPositiveXGradient*dxdydz[1] + ddzPositiveXGradient*dxdydz[2];
    myGradient[1] = firstQuadrantGradient[1] + ddxPositiveYGradient*dxdydz[0] + ddzPositiveYGradient*dxdydz[2];
    myGradient[2] = firstQuadrantGradient[2] + ddxPositiveZGradient*dxdydz[0] + ddyPositiveZGradient*dxdydz[1];
    return myGradient;
}

fVec3 GridPoint::getFirstQuadrantGradient(){
    return firstQuadrantGradient;
}
*/

GridIndices::GridIndices( int myXIndex,  int myYIndex,  int myZIndex) {
		xGridPoint = ValidateInt(myXIndex);
		yGridPoint = ValidateInt(myYIndex);
		zGridPoint = ValidateInt(myZIndex);
	}
int GridIndices::getXGridIndex () const {return xGridPoint;}
int GridIndices::getYGridIndex () const {return yGridPoint;}
int GridIndices::getZGridIndex () const {return zGridPoint;}

DensityMap::DensityMap(){
    initializeMap();
};

DensityMap::~DensityMap(){
    ArrayOfGridPoints.clear();
};
void DensityMap::initializeMap() {
    unitCellNumGridX = 0;
    unitCellNumGridY = 0;
    unitCellNumGridZ = 0;
    totalNumGridX = 0;
    totalNumGridY = 0;
    totalNumGridZ = 0;
    minX = 0;
    minY = 0;
    minZ = 0;
    maxX = 0;
    maxY = 0;
    maxY = 0;
    gridXSpacing = 0;
    gridYSpacing = 0;
    gridZSpacing = 0;
};

void DensityMap::validateGridParameters() {
    cout<<__FILE__<<":"<<__LINE__<< " Checking maxX minX  gridXSpacing totalNumGridX (in nm, nm, nm, unitless): "<<maxX<<" , "<<minX<<" , "<<gridXSpacing<<" , "<<totalNumGridX<<std::endl;
    cout<<__FILE__<<":"<<__LINE__<< " Checking maxY minY  gridYSpacing totalNumGridY (in nm, nm, nm, unitless): "<<maxY<<" , "<<minY<<" , "<<gridYSpacing<<" , "<<totalNumGridY<<std::endl;
    cout<<__FILE__<<":"<<__LINE__<< " Checking maxZ minZ  gridZSpacing totalNumGridZ (in nm, nm, nm, unitless): "<<maxZ<<" , "<<minZ<<" , "<<gridZSpacing<<" , "<<totalNumGridZ<<std::endl;
    /*cout<<__FILE__<<":"<<__LINE__<< " int((maxX-minX)/gridXSpacing +1) "<<int((maxX-minX)/gridXSpacing +1)<<endl;
    cout<<__FILE__<<":"<<__LINE__<< " int((maxY-minY)/gridYSpacing +1) "<<int((maxY-minY)/gridYSpacing +1)<<endl;
    cout<<__FILE__<<":"<<__LINE__<< " int((maxZ-minZ)/gridZSpacing +1) "<<int((maxZ-minZ)/gridZSpacing +1)<<endl;
    cout<<__FILE__<<":"<<__LINE__<< " ((maxX-minX)/gridXSpacing +1 "<<0.0+((maxX-minX)/gridXSpacing +1.0)<<endl;
    cout<<__FILE__<<":"<<__LINE__<< " ((maxY-minY)/gridYSpacing +1 "<<0.0+((maxY-minY)/gridYSpacing +1.0)<<endl;
    cout<<__FILE__<<":"<<__LINE__<< " ((maxZ-minZ)/gridZSpacing +1 "<<0.0+((maxZ-minZ)/gridZSpacing +1.0)<<endl;*/
    //cout<<__FILE__<<":"<<__LINE__<< " Checking maxX minX  gridXSpacing : "<<maxX<<" , "<<minX<<" , "<<gridXSpacing<<std::endl;
    if (
	/*(int((maxX-minX)/gridXSpacing +1) != totalNumGridX  ) ||
        (int((maxY-minY)/gridYSpacing +1) != totalNumGridY  ) ||
        (int((maxZ-minZ)/gridZSpacing +1) != totalNumGridZ  ) */
	( abs(((maxX-minX)/gridXSpacing +1) - totalNumGridX) > 1E-7  ) ||
        ( abs(((maxY-minY)/gridYSpacing +1) - totalNumGridY) > 1E-7  ) ||
        ( abs(((maxZ-minZ)/gridZSpacing +1) - totalNumGridZ) > 1E-7  ) 
       ) {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<< " There is a problem with the max- , min- (XYZ) or totalNumGrid (XYZ) or grid (XYZ) Spacing parameters "<< ((maxX-minX)/gridXSpacing +1) - totalNumGridX <<endl;	
        ErrorManager::instance.treatError();
    }
}
		
bool        DensityMap::hasGridPoint(GridIndices myGridIndices){

        //cout<<__FILE__<<":"<<__LINE__<< " ArrayOfGridPoints[0].size() = "  <<ArrayOfGridPoints[0].size() <<endl;
        //cout<<__FILE__<<":"<<__LINE__<< " ArrayOfGridPoints[0][0].size() = "  <<ArrayOfGridPoints[0][0].size() <<endl;
	bool myHasGridPoint;
	if (myGridIndices.getZGridIndex() >= totalNumGridZ) 
        {
		//cout<<__FILE__<<":"<<__LINE__<< " The Z-index of "<<myGridIndices.getZGridIndex()<<" exceeds the Z-dimension, "<< totalNumGridZ <<" of the grid map. No force applied."<<endl;	
		myHasGridPoint = false;
	}
	else if (myGridIndices.getZGridIndex() <  0 ) {
		//cout<<__FILE__<<":"<<__LINE__<< " The Z-index of "<<myGridIndices.getZGridIndex()<<" is less than zero. No force applied."<<endl; 
		myHasGridPoint = false;
	}
	else if (myGridIndices.getYGridIndex() >= totalNumGridY) {
		//cout<<__FILE__<<":"<<__LINE__<< " The Y-index of "<<myGridIndices.getYGridIndex()<<" exceeds the Y-dimension, "<<ArrayOfGridPoints[myGridIndices.getZGridIndex()].size() <<" of the grid map. No force applied."<<endl;	
		myHasGridPoint = false;
	}
	else if (myGridIndices.getYGridIndex() <  0 ) {
		//cout<<__FILE__<<":"<<__LINE__<< " The Y-index of "<<myGridIndices.getYGridIndex()<<" is less than zero. No force applied."<<endl; 
		myHasGridPoint = false;
	}
	else if (myGridIndices.getXGridIndex() >= totalNumGridX) {
		//cout<<__FILE__<<":"<<__LINE__<< " The X-index of "<<myGridIndices.getXGridIndex()<<" exceeds the X-dimension, "<<totalNumGridX <<" of the grid map. No force applied."<<endl;	
		myHasGridPoint = false;
	}
	else if (myGridIndices.getXGridIndex() <  0 ) {
		//cout<<__FILE__<<":"<<__LINE__<< " The X-index of "<<myGridIndices.getXGridIndex()<<" is less than zero. No force applied."<<endl; 
		myHasGridPoint = false;
	}
	else { myHasGridPoint = true;}
	return myHasGridPoint;

}


GridPoint & DensityMap::updGridPoint(GridIndices myGridIndices){
		if (! hasGridPoint(myGridIndices)) {
			ErrorManager::instance <<__FILE__<<":"<<__LINE__<<	" No nearby grid point.  The point you requested "<< myGridIndices.getZGridIndex()<<" , "<<  myGridIndices.getYGridIndex()<<" , "<< myGridIndices.getXGridIndex() << " is off the density map.";
			ErrorManager::instance.treatError();
		} else {
			GridPoint & myGridPoint = ArrayOfGridPoints[myGridIndices.getZGridIndex()][myGridIndices.getYGridIndex()][myGridIndices.getXGridIndex()] ;
			//myGridPoint.validate();
			return myGridPoint;
		}
  		}
   
const GridPoint DensityMap::getGridPoint(GridIndices myGridIndices)  {
		//return & ArrayOfGridPoints[0][0][0];
		return ArrayOfGridPoints[myGridIndices.getZGridIndex()][myGridIndices.getYGridIndex()][myGridIndices.getXGridIndex()] ;
  		}
   
void	    DensityMap::validateGridPoint(GridIndices myGridIndices){
			if (! hasGridPoint(myGridIndices)) {
				ErrorManager::instance <<__FILE__<<":"<<__LINE__<<	" No nearby grid point.  The point you requested is off the density map.";
				ErrorManager::instance.treatError();
			} else validate(updGridPoint(myGridIndices));

}


const bool DensityMap::hasNearbyGridIndices(fVec3 position)
		{
			int tempXIndex 	= 	int((position[0] + gridXSpacing/2)/gridXSpacing);
			int tempYIndex 	= 	int((position[1] + gridYSpacing/2)/gridYSpacing);
			int tempZIndex 	= 	int((position[2] + gridZSpacing/2)/gridZSpacing);
			//cout<<__FILE__<<":"<<__LINE__<<	" in No nearby grid point.  The point you requested is off the density map.";
			GridIndices tempGridIndices(tempXIndex, tempYIndex, tempZIndex);
			if (hasGridPoint (tempGridIndices)) {
				return true;            
			} else {
				return false;
			}
		}
		

GridIndices DensityMap::calcNearestGridIndices(fVec3 position)
		{
			int tempXIndex 	= int (((position[0]-minX) )/gridXSpacing);
			int tempYIndex 	= int (((position[1]-minY) )/gridYSpacing);
			int tempZIndex 	= int (((position[2]-minZ) )/gridZSpacing);
			//int tempYIndex 	= 	int((position[1] + gridYSpacing/2)/gridYSpacing);
			//int tempZIndex 	= 	int((position[2] + gridZSpacing/2)/gridZSpacing);
			GridIndices tempGridIndices(tempXIndex, tempYIndex, tempZIndex);
			return tempGridIndices;
		}
		
GridIndices DensityMap::calcLowerLeftGridIndices(fVec3 position)
                {
                        int tempXIndex  = int (floor(((position[0]-minX) )/gridXSpacing));
                        int tempYIndex  = int (floor(((position[1]-minY) )/gridYSpacing));
                        int tempZIndex  = int (floor(((position[2]-minZ) )/gridZSpacing));

                        return GridIndices(tempXIndex, tempYIndex, tempZIndex);

                }
		
const GridPoint DensityMap::getGridPoint(fVec3 myPosition)  {
	return getGridPoint(calcNearestGridIndices( myPosition));				
}

GridPoint & DensityMap::updGridPoint(fVec3 myPosition)  {
	return updGridPoint(calcNearestGridIndices( myPosition));				
}

const double DensityMap::getDensity(fVec3 myPosition) {
	if ((myPosition[0] < minX) || (myPosition[0] > maxX) ||	
	   (myPosition[1] < minY) || (myPosition[1] > maxY)  ||	
	   (myPosition[2] < minZ) || (myPosition[2] > maxZ)) // if outside bounds of density map
	{
		return 0.0; //return zero density
	} else {
        //if (hasNearbyGridIndices(myPosition)) 
                GridIndices myLowerLeftGridIndices = calcLowerLeftGridIndices(myPosition);
                
                if (hasGridPoint(myLowerLeftGridIndices)) {
                    //fVec3 myDensity = getGridPoint(myLowerLeftGridIndices).getDensity(myPosition);
                    //ValidatefVec3(myDensity); 
		    //return myDensity;
		    return getDensity(updGridPoint(myLowerLeftGridIndices),myPosition);
                } else {
                    return 0.0;
                }
		//return getGridPoint(myPosition).getDensity(myPosition);
	} 
}

const double DensityMap::getDensity(SimTK::Vec3 myPosition) {
    return getDensity(fVec3(myPosition[0], myPosition[1],myPosition[2] ));
}
void DensityMap::initializeArrayOfGridPoints(){
        validateGridParameters();
        fVec3 tempPosition(0,0,0);
        ArrayOfGridPoints.resize(totalNumGridZ);
        for ( int zIndex = 0; zIndex < totalNumGridZ; zIndex++) {
        	ArrayOfGridPoints[zIndex].resize(totalNumGridY);
        	for ( int yIndex = 0; yIndex < totalNumGridY; yIndex++) {
        		ArrayOfGridPoints[zIndex][yIndex].resize(totalNumGridX);
        	}}
        for ( int xIndex = 0; xIndex < totalNumGridX; xIndex++) {
        	for ( int yIndex = 0; yIndex < totalNumGridY; yIndex++) {
        		for ( int zIndex = 0; zIndex < totalNumGridZ; zIndex++) {
        			GridPoint tempGridPoint;// = updGridPoint(GridIndices(xIndex, yIndex, zIndex)) ;
        			initialize(tempGridPoint);
        			tempPosition = fVec3(xIndex * gridXSpacing + minX, yIndex * gridYSpacing + minY,zIndex * gridZSpacing + minZ);
        			setPosition(tempGridPoint,tempPosition);
        			validate(tempGridPoint);
        			updGridPoint(GridIndices(xIndex, yIndex, zIndex)) = tempGridPoint;
        			validate(updGridPoint(GridIndices(xIndex, yIndex, zIndex))); // just being paranoid
        			//ArrayOfGridPoints[zIndex, yIndex, xIndex] = tempGridPoint;			
        }}}
        if (ArrayOfGridPoints.size() != totalNumGridZ) {
            ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Wrong number of grid points in Z direction! Found :"<< ArrayOfGridPoints.size()<<" expected : " << totalNumGridZ<<","<< totalNumGridZ<<endl;
            ErrorManager::instance.treatError();
        }
        for ( int zIndex = 0; zIndex < totalNumGridZ; zIndex++) {
        	if (ArrayOfGridPoints[zIndex].size() != totalNumGridY) {
                   ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Wrong number of grid points in Y direction! Found :"<< ArrayOfGridPoints[zIndex].size()<<" expected : " << totalNumGridY<<","<< totalNumGridZ<<endl;
                    ErrorManager::instance.treatError();
                }
        	for ( int yIndex = 0; yIndex < totalNumGridY; yIndex++) {
        	    if (ArrayOfGridPoints[zIndex][yIndex].size() != (totalNumGridX) ) {
        		             ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Wrong number of grid points in X direction! Found :"<< ArrayOfGridPoints[zIndex][yIndex].size()<<" expected : " << totalNumGridX<<endl;
                             ErrorManager::instance.treatError();
                }
        	}
        } // of for zIndex
}

void DensityMap::loadParametersAndDensity(const String densityFileName)
{
    unsigned int extIndex = densityFileName.rfind(".");
    String extension = densityFileName.substr(extIndex);
    if(extension == ".xplor")
        loadParametersAndDensity_XPLOR(densityFileName);
    else if(extension == ".dx")
        loadParametersAndDensity_OpenDX(densityFileName);
    else if(extension == ".situs" || extension == ".sit")
        loadParametersAndDensity_Situs(densityFileName);
    else
    {
        ErrorManager::instance << "DensityMap: Extension unknown for " << densityFileName << endl;
        ErrorManager::instance.treatError();
    }
}

// Expects .xplor density maps to be in the following format:
// http://psb11.snv.jussieu.fr/doc-logiciels/msi/xplor981/formats.html -- not accessible anymore
// http://chem5.nchc.org.tw/software/document2007/insightII/doc/xplor/formats.html
void DensityMap::loadParametersAndDensity_XPLOR(const String densityFileName) {

        ifstream inFile(densityFileName.c_str(),ifstream::in);
        int densitiesPerLine = 6;	
        	//const int numFields = 10 ;
        vector<String> mystring;
        //String mystring2[numFields];
        	stringstream u;
        	String tempString;
        if (! (inFile.is_open())) {

            ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Unable to open density map : "<<densityFileName<<endl;
            ErrorManager::instance.treatError();
        }
        readAndParseLine  (inFile);
        readAndParseLine  (inFile);
        readAndParseLine  (inFile);
        readAndParseLine  (inFile);
        mystring = readAndParseLine  (inFile);
        unitCellNumGridX = ValidateNonNegativeInt (atoi (((mystring[0].c_str()))) );
        unitCellNumGridY = ValidateNonNegativeInt (atoi (((mystring[3].c_str()))) );
        unitCellNumGridZ = ValidateNonNegativeInt (atoi (((mystring[6].c_str()))) );
        // are these being read in the right order?
        totalNumGridX = ValidateNonNegativeInt (atoi (((mystring[2].c_str()))) - atoi (((mystring[1].c_str()))) + 1);
        totalNumGridY = ValidateNonNegativeInt (atoi (((mystring[5].c_str()))) - atoi (((mystring[4].c_str()))) + 1);
        totalNumGridZ = ValidateNonNegativeInt (atoi (((mystring[8].c_str()))) - atoi (((mystring[7].c_str()))) + 1);
        cout<<__FILE__<<":"<<__LINE__<<" set totalNumGridZ,Y,X: "<<totalNumGridX<<","<< totalNumGridY<<","<< totalNumGridZ<<endl;

        cout<<__FILE__<<":"<<__LINE__<<" "<< mystring[1]<<endl;
        minX = atoi(mystring[1].c_str());
        minY = atoi(mystring[4].c_str());
        minZ = atoi(mystring[7].c_str());

        cout<<__FILE__<<":"<<__LINE__<<" "<< mystring[2]<<endl;
        maxX = atoi((mystring[2].c_str()));	//in grid points
        maxY = atoi((mystring[5].c_str()));	//in grid units
        maxZ = atoi((mystring[8].c_str()));	//in grid units

        //initializeArrayOfGridPoints();
        mystring = readAndParseLine  (inFile);
        cout<<__FILE__<<":"<<__LINE__<<" "<< mystring[0]<<endl;
        gridXSpacing = (atof (mystring[0].c_str()))/unitCellNumGridX/10 ;  //dividing by 10 to convert from Å(the XPLOR format, per Alwyn Jones) to nm (molmodel units).   
        gridYSpacing = (atof (mystring[1].c_str()))/unitCellNumGridY/10 ;   
        gridZSpacing = (atof (mystring[2].c_str()))/unitCellNumGridZ/10 ;   

        cout<<__FILE__<<":"<<__LINE__<<" set grid Spacing in Z,Y,X to: "<<gridZSpacing*10<<","<<gridYSpacing*10<<","<<gridXSpacing*10<<" (Å,Å,Å)"<<endl;

        minX = minX * gridXSpacing;
        minY = minY * gridYSpacing;
        minZ = minZ * gridZSpacing;

        maxX = maxX * gridXSpacing;
        maxY = maxY * gridYSpacing;
        maxZ = maxZ * gridZSpacing;
        initializeArrayOfGridPoints();

        cout<<__FILE__<<":"<<__LINE__<<" minimum extents in  Z,Y,X : "<<minZ*10<<","<< minY*10<<","<< minX*10<<endl;
        cout<<__FILE__<<":"<<__LINE__<<" maximum extents in  Z,Y,X : "<<maxZ*10<<","<< maxY*10<<","<< maxX*10<<endl;
        if (atof (     (mystring[3].c_str())) != 90) {
        			ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Found alpha = "<< mystring[3]<<".  At the moment we can only handle cubic lattices (alpha = 90)." <<endl;
        		ErrorManager::instance.treatError();
        }
        if (atof (     (mystring[4].c_str())) != 90) {
        			ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Found beta = "<< mystring[4]<<".  At the moment we can only handle cubic lattices (beta = 90)." <<endl;
        		ErrorManager::instance.treatError();
        }
        if (atof (     (mystring[5].c_str())) != 90) {
        			ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Found gamma = "<< mystring[5]<<".  At the moment we can only handle cubic lattices (gamma = 90)." <<endl;
        		ErrorManager::instance.treatError();
        }
        mystring = readAndParseLine  (inFile);
        if ((mystring[0]).compare("ZYX") != 0    ) {ErrorManager::instance <<__FILE__<<":"<<__LINE__<<"expected ZYX, got :"<<mystring[0]<<endl; ErrorManager::instance.treatError();}


        for ( int zIndex = 0; zIndex < totalNumGridZ; zIndex ++) {
        	// read z-index
                mystring = readAndParseLine  (inFile);
        	if (atoi ((mystring[0].c_str())) != zIndex) {
        			ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Something is wrong with the input file reading.  expected to read zIndex "<<zIndex<< " and instead read : "<<mystring[0]<<endl;
        		ErrorManager::instance.treatError();
        	}	
        	

        	for ( int yIndex = 0; yIndex < totalNumGridY; yIndex ++) 
        	for ( int xIndex = 0; xIndex < totalNumGridX; xIndex = xIndex + 0) 

        	{
        		mystring = readAndParseOnColWidth (inFile,12);
        		if ((int)mystring.size() > densitiesPerLine) {
        				ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Too many densities on this line!  expected "<<densitiesPerLine<< "and found "<<mystring.size()<<endl;
        			ErrorManager::instance.treatError();
        		}
        		//xIndex --;
        		for (int i = 0; i<(int)mystring.size() ;i++) {
        			if (xIndex == totalNumGridX){xIndex = 0; yIndex++; }
        			if ((xIndex+1 + (yIndex*totalNumGridX) ) <= (totalNumGridX* totalNumGridY))
        			{
        				GridIndices centralGridIndices((xIndex ),  yIndex, zIndex);
        				GridPoint & centralGridPoint = updGridPoint(centralGridIndices);
        				setDensity(centralGridPoint, atof(mystring[i].c_str()));	
        			} else{ 
        					ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Something went wrong .. too many densities !  "<<endl;
        				ErrorManager::instance.treatError();
        		        }
        			xIndex++;
        		}
        		
        	}
        } // of zIndex
        inFile.close();
}

// Expects Situs density maps
// http://situs.biomachina.org/fmap.pdf
void DensityMap::loadParametersAndDensity_Situs(const String densityFileName) {

        ifstream inFile(densityFileName.c_str(),ifstream::in);
        int densitiesPerLine = 10;   
            //const int numFields = 10 ;
        vector<String> mystring;
        //String mystring2[numFields];
            stringstream u;
            String tempString;
        if (! (inFile.is_open())) {

            ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Unable to open density map : "<<densityFileName<<endl;
            ErrorManager::instance.treatError();
        }
        mystring = readAndParseLine  (inFile);
        gridXSpacing = atof(mystring[0].c_str())/10.0;
        gridYSpacing = gridXSpacing;
        gridZSpacing = gridXSpacing;
        
        unitCellNumGridX = ValidateNonNegativeInt (atoi (((mystring[4].c_str()))) );
        unitCellNumGridY = ValidateNonNegativeInt (atoi (((mystring[5].c_str()))) );
        unitCellNumGridZ = ValidateNonNegativeInt (atoi (((mystring[6].c_str()))) );
        // are these being read in the right order?
        totalNumGridX = unitCellNumGridX;
        totalNumGridY = unitCellNumGridY;
        totalNumGridZ = unitCellNumGridZ;
        cout<<__FILE__<<":"<<__LINE__<<" set totalNumGridZ,Y,X: "<<totalNumGridX<<","<< totalNumGridY<<","<< totalNumGridZ<<endl;

        cout<<__FILE__<<":"<<__LINE__<<" "<< mystring[1]<<endl;
        minX = atoi(mystring[1].c_str());
        minY = atoi(mystring[2].c_str());
        minZ = atoi(mystring[3].c_str());

        cout<<__FILE__<<":"<<__LINE__<<" set grid Spacing in Z,Y,X to: "<<gridZSpacing*10<<","<<gridYSpacing*10<<","<<gridXSpacing*10<<" (Å,Å,Å)"<<endl;

        minX = minX/10.0;
        minY = minY/10.0;
        minZ = minZ/10.0;

        maxX = minX + ((totalNumGridX-1) * gridXSpacing);
        maxY = minY + ((totalNumGridY-1) * gridYSpacing);
        maxZ = minZ + ((totalNumGridZ-1) * gridZSpacing);
        initializeArrayOfGridPoints();

        cout<<__FILE__<<":"<<__LINE__<<" minimum extents in  Z,Y,X : "<<minZ*10<<","<< minY*10<<","<< minX*10<<endl;
        cout<<__FILE__<<":"<<__LINE__<<" maximum extents in  Z,Y,X : "<<maxZ*10<<","<< maxY*10<<","<< maxX*10<<endl;

        readAndParseLine  (inFile);



        for ( int zIndex = 0; zIndex < totalNumGridZ; zIndex ++) 
        {
            for ( int yIndex = 0; yIndex < totalNumGridY; yIndex ++) 
            {   
                for ( int xIndex = 0; xIndex < totalNumGridX; xIndex = xIndex + 0) 
                {
                    mystring = readAndParseLine(inFile);
                    if ((int)mystring.size()-1 > densitiesPerLine) {
                        cout << __FILE__<<":"<<__LINE__<<": " << mystring[0]<< " " << mystring[10] << endl;
                        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Too many densities on this line!  expected "<<densitiesPerLine<< " and found "<<mystring.size()<<endl;
                        ErrorManager::instance.treatError();
                    }
                    //xIndex --;
                    for (int i = 0; i<(int)mystring.size()-1 ;i++) {
                        if (xIndex == totalNumGridX){
                            xIndex = 0; yIndex++; 
                        }
                        if ((xIndex+1 + (yIndex*totalNumGridX) ) <= (totalNumGridX* totalNumGridY))
                        {
                            GridIndices centralGridIndices((xIndex ),  yIndex, zIndex);
                            GridPoint & centralGridPoint = updGridPoint(centralGridIndices);
                            setDensity(centralGridPoint, atof(mystring[i].c_str()));    
                        } else{
                            cout <<  xIndex+1 + (yIndex*totalNumGridX) << " " << totalNumGridX* totalNumGridY << endl;
                            ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Something went wrong .. too many densities !  "<<endl;
                            ErrorManager::instance.treatError();
                            }
                        xIndex++;
                    }
                
                }
            }
        } // of zIndex
        inFile.close();
}

void DensityMap::precomputeGradient() {
			for ( int xIndex = 0; xIndex < totalNumGridX; xIndex ++) 
			for ( int yIndex = 0; yIndex < totalNumGridY; yIndex ++) 
			for ( int zIndex = 0; zIndex < totalNumGridZ; zIndex ++) 

			{
				GridIndices centralGridIndices(xIndex,  yIndex, zIndex);
				GridPoint & centralGridPoint = updGridPoint(centralGridIndices);
				initializeGradient(centralGridPoint); // sets gradient to zero; if the if statements below are false that component of the gradient remains zero.

				if (xIndex < (totalNumGridX-1)) {// if grid point is not at the +X boundary
					setPositiveXGradient (centralGridPoint,(getDensity(updGridPoint(GridIndices(xIndex+1,  yIndex, zIndex))) -getDensity(centralGridPoint))/ gridXSpacing);
				} else if (xIndex == (totalNumGridX-1)) { // if grid point is at the +X boundary
					setPositiveXGradient (centralGridPoint,(0. -getDensity(centralGridPoint)) / gridXSpacing );
                                }

				if (yIndex < (totalNumGridY-1)) {// if grid point is not at the +Y boundary
					setPositiveYGradient (centralGridPoint,(getDensity(updGridPoint(GridIndices(xIndex,  yIndex+1, zIndex))) -getDensity(centralGridPoint))/ gridYSpacing);
				} else if (yIndex == (totalNumGridY-1)) {
					setPositiveYGradient (centralGridPoint,(0. -getDensity(centralGridPoint)) / gridYSpacing );
                                }



				if (zIndex < (totalNumGridZ-1)) {// if grid point is not at the +Z boundary:
					setPositiveZGradient (centralGridPoint,(getDensity(updGridPoint(GridIndices(xIndex,  yIndex, zIndex+1))) -getDensity(centralGridPoint))/ gridZSpacing);
				
				} else if (zIndex == (totalNumGridZ-1)) {
					setPositiveZGradient (centralGridPoint,(0. -getDensity(centralGridPoint)) / gridZSpacing );
                                }
				//ValidateVec3((fetchFirstQuadrantGradient(centralGridPoint)));



			}
		}

void DensityMap::precomputeGradientDerivatives() {
			for ( int xIndex = 0; xIndex < totalNumGridX; xIndex ++) 
			for ( int yIndex = 0; yIndex < totalNumGridY; yIndex ++) 
			for ( int zIndex = 0; zIndex < totalNumGridZ; zIndex ++) 

			{
				GridIndices centralGridIndices(xIndex,  yIndex, zIndex);
				GridPoint & centralGridPoint = updGridPoint(centralGridIndices);

				if (xIndex < (totalNumGridX-1)) {// if grid point is not at the +X boundary
					setddxPositiveYGradient (centralGridPoint,(fetchFirstQuadrantGradient(updGridPoint(GridIndices(xIndex+1,  yIndex, zIndex)))[1] -fetchFirstQuadrantGradient(centralGridPoint)[1])/ gridXSpacing);
					setddxPositiveZGradient (centralGridPoint,(fetchFirstQuadrantGradient(updGridPoint(GridIndices(xIndex+1,  yIndex, zIndex)))[2] -fetchFirstQuadrantGradient(centralGridPoint)[2])/ gridXSpacing);
				} else if (xIndex == (totalNumGridX-1)) { // if grid point is at the + X boundary
					setddxPositiveYGradient (centralGridPoint,(0 - fetchFirstQuadrantGradient(centralGridPoint)[1])/ gridXSpacing);
					setddxPositiveZGradient (centralGridPoint,(0 -fetchFirstQuadrantGradient(centralGridPoint)[2])/ gridXSpacing);
                                }


				if (yIndex < (totalNumGridY-1)) {// if grid point is not at the +Y boundary
					setddyPositiveXGradient (centralGridPoint,(fetchFirstQuadrantGradient(updGridPoint(GridIndices(xIndex,  yIndex+1, zIndex)))[0] -fetchFirstQuadrantGradient(centralGridPoint)[0])/ gridYSpacing);
					setddyPositiveZGradient (centralGridPoint,(fetchFirstQuadrantGradient(updGridPoint(GridIndices(xIndex,  yIndex+1, zIndex)))[2] -fetchFirstQuadrantGradient(centralGridPoint)[2])/ gridYSpacing);
				
				} else if (yIndex == (totalNumGridY-1)) {
					setddyPositiveXGradient (centralGridPoint,(0 - fetchFirstQuadrantGradient(centralGridPoint)[0])/ gridYSpacing);
					setddyPositiveZGradient (centralGridPoint,(0 -fetchFirstQuadrantGradient(centralGridPoint)[2])/ gridYSpacing);
                                }

				if (zIndex < (totalNumGridZ-1)) {// if grid point is not at the +Z boundary:
					setddzPositiveXGradient ( centralGridPoint,(fetchFirstQuadrantGradient(updGridPoint(GridIndices(xIndex,  yIndex, zIndex+1)))[0] -fetchFirstQuadrantGradient(centralGridPoint)[0])/ gridZSpacing);
					setddzPositiveYGradient ( centralGridPoint, ( fetchFirstQuadrantGradient(updGridPoint(GridIndices(xIndex,  yIndex, zIndex+1)))[1] -fetchFirstQuadrantGradient(centralGridPoint)[1])/ gridZSpacing);
				
				} else if (zIndex == (totalNumGridZ-1)) {
					setddzPositiveXGradient (centralGridPoint, (0 - fetchFirstQuadrantGradient(centralGridPoint)[0])/ gridZSpacing);
					setddzPositiveYGradient (centralGridPoint, (0 -fetchFirstQuadrantGradient(centralGridPoint)[1])/ gridZSpacing);
                                }

			}
		}


/*
fVec3 DensityMap::fetchGradient(fVec3 position)  {
			GridIndices myNearestGridIndices = calcNearestGridIndices(   position);
                        if (hasGridPoint(myNearestGridIndices)) {
				GridPoint myGridPoint =  updGridPoint(myNearestGridIndices);
				return myGridPoint.fetchGradient (position) ;	
			} else {
				return fVec3(0);
			}
		}
		
*/


fVec3 DensityMap::fetchGradient(fVec3 position)  {
                     	GridIndices myNearestGridIndices = calcNearestGridIndices(   position);
                        //GridPoint myGridPoint; 
     			//myGridPoint.initialize();
                         if (hasGridPoint(myNearestGridIndices)) {
                                 //myGridPoint =  updGridPoint(myNearestGridIndices);
                                 return fetchGradient(updGridPoint(myNearestGridIndices), position) ;
                         } else if (hasGridPoint(GridIndices (myNearestGridIndices.getXGridIndex()+1, myNearestGridIndices.getYGridIndex(),myNearestGridIndices.getZGridIndex()))) {
                                GridPoint myGridPoint; 
     			        initialize(myGridPoint);
				setPositiveXGradient(myGridPoint,(getDensity(updGridPoint(GridIndices(myNearestGridIndices.getXGridIndex() +1,  myNearestGridIndices.getYGridIndex() , myNearestGridIndices.getZGridIndex() ))) - 0.) / gridXSpacing) ;	
         			
          			ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" this function is obsolete, fetchGradient doesn't work anymore"<<endl;
        			ErrorManager::instance.treatError();
                         } else if (hasGridPoint(GridIndices (myNearestGridIndices.getXGridIndex()-1, myNearestGridIndices.getYGridIndex(),myNearestGridIndices.getZGridIndex()))) {
                                GridPoint myGridPoint; 
     			        initialize(myGridPoint);
				setNegativeXGradient(myGridPoint,(0. - getDensity(updGridPoint(GridIndices(myNearestGridIndices.getXGridIndex() -1,  myNearestGridIndices.getYGridIndex() , myNearestGridIndices.getZGridIndex() )))) / gridXSpacing) ;	 
          			ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" this function is obsolete, fetchGradient doesn't work anymore"<<endl;
        			ErrorManager::instance.treatError();
                         } else if (hasGridPoint(GridIndices (myNearestGridIndices.getXGridIndex(), myNearestGridIndices.getYGridIndex()+1,myNearestGridIndices.getZGridIndex()))) {
                                GridPoint myGridPoint; 
     			        initialize(myGridPoint);
				setPositiveYGradient(myGridPoint,(getDensity(updGridPoint(GridIndices(myNearestGridIndices.getXGridIndex() ,  myNearestGridIndices.getYGridIndex() +1, myNearestGridIndices.getZGridIndex() ))) - 0.) / gridYSpacing) ;	
                                return fetchGradient (myGridPoint,position) ;
                         } else if (hasGridPoint(GridIndices (myNearestGridIndices.getXGridIndex(), myNearestGridIndices.getYGridIndex()-1,myNearestGridIndices.getZGridIndex()))) {
                                GridPoint myGridPoint; 
     			        initialize(myGridPoint);
				setNegativeYGradient(myGridPoint ,(0. - getDensity(updGridPoint(GridIndices(myNearestGridIndices.getXGridIndex() -1,  myNearestGridIndices.getYGridIndex() -1 , myNearestGridIndices.getZGridIndex() )))) / gridYSpacing) ;
                                return fetchGradient(myGridPoint,position) ;
                         } else if (hasGridPoint(GridIndices (myNearestGridIndices.getXGridIndex(), myNearestGridIndices.getYGridIndex(),myNearestGridIndices.getZGridIndex()+1))) {
                                GridPoint myGridPoint; 
     			        initialize(myGridPoint );
				setPositiveZGradient(myGridPoint,(getDensity(updGridPoint(GridIndices(myNearestGridIndices.getXGridIndex() ,  myNearestGridIndices.getYGridIndex() , myNearestGridIndices.getZGridIndex() +1 ))) - 0.) / gridXSpacing) ;	
                                 return fetchGradient (myGridPoint,position) ;
                         } else if (hasGridPoint(GridIndices (myNearestGridIndices.getXGridIndex(), myNearestGridIndices.getYGridIndex(),myNearestGridIndices.getZGridIndex()-1))) {
                                GridPoint myGridPoint; 
     			        initialize(myGridPoint);
				setNegativeZGradient(myGridPoint,(0. - getDensity(updGridPoint(GridIndices(myNearestGridIndices.getXGridIndex() ,  myNearestGridIndices.getYGridIndex() , myNearestGridIndices.getZGridIndex() - 1 )))) / gridXSpacing) ;
                                 return fetchGradient (myGridPoint,position) ;
                         } else {
                                 return fVec3(0);
                         }
                 }



fVec3 DensityMap::fetchFirstQuadrantGradient(fVec3 position)  {

                        GridIndices myLowerLeftGridIndex = calcLowerLeftGridIndices(   position);
                         if (hasGridPoint(myLowerLeftGridIndex)) {
                                 return fetchFirstQuadrantGradient(ArrayOfGridPoints[myLowerLeftGridIndex.getZGridIndex()][myLowerLeftGridIndex.getYGridIndex()][myLowerLeftGridIndex.getXGridIndex()]);

                         } else if (hasGridPoint(GridIndices (myLowerLeftGridIndex.getXGridIndex()+1, myLowerLeftGridIndex.getYGridIndex(),myLowerLeftGridIndex.getZGridIndex()))) {
                                GridPoint myGridPoint;
                                initialize(myGridPoint);
                                setPositiveXGradient(myGridPoint,(getDensity(updGridPoint(GridIndices(myLowerLeftGridIndex.getXGridIndex() +1,  myLowerLeftGridIndex.getYGridIndex() , myLowerLeftGridIndex.getZGridIndex() ))) - 0.) / gridXSpacing) ;
                                return fetchFirstQuadrantGradient(myGridPoint) ;

                         } else if (hasGridPoint(GridIndices (myLowerLeftGridIndex.getXGridIndex(), myLowerLeftGridIndex.getYGridIndex()+1,myLowerLeftGridIndex.getZGridIndex()))) {
                                GridPoint myGridPoint;
                                initialize(myGridPoint);
                                setPositiveYGradient(myGridPoint,(getDensity(updGridPoint(GridIndices(myLowerLeftGridIndex.getXGridIndex() ,  myLowerLeftGridIndex.getYGridIndex() +1, myLowerLeftGridIndex.getZGridIndex() ))) - 0.) / gridYSpacing) ;
                                return fetchFirstQuadrantGradient(myGridPoint) ;

                         } else if (hasGridPoint(GridIndices (myLowerLeftGridIndex.getXGridIndex(), myLowerLeftGridIndex.getYGridIndex(),myLowerLeftGridIndex.getZGridIndex()+1))) {
                                GridPoint myGridPoint;
                                initialize(myGridPoint);
                                setPositiveZGradient(myGridPoint,(getDensity(updGridPoint(GridIndices(myLowerLeftGridIndex.getXGridIndex() ,  myLowerLeftGridIndex.getYGridIndex() , myLowerLeftGridIndex.getZGridIndex() +1 ))) - 0.) / gridXSpacing) ;
                                // return myGridPoint.fetchGradient (position) ;
                                return fetchFirstQuadrantGradient(myGridPoint) ;

                         } else {
                                 fVec3 tempVec3(0);
                                 return tempVec3;
                         }

}

fVec3 DensityMap::calcInterpolatedFirstQuadrantGradient(fVec3 position)  {

                        GridIndices myLowerLeftGridIndex = calcLowerLeftGridIndices(   position);
                         if (hasGridPoint(myLowerLeftGridIndex)) {
                                 return calcInterpolatedFirstQuadrantGradient(ArrayOfGridPoints[myLowerLeftGridIndex.getZGridIndex()][myLowerLeftGridIndex.getYGridIndex()][myLowerLeftGridIndex.getXGridIndex()],position);

                         } 
                         else { // might want to trap the conditions at the boundaries of the map, to get the minimizer to work
                                 fVec3 tempVec3(0);
                                 return tempVec3;
                         }
}

SimTK::Vec3 DensityMap::calcInterpolatedFirstQuadrantGradient(SimTK::Vec3 position)  {
    fVec3 myVec3Float = calcInterpolatedFirstQuadrantGradient(fVec3(  position[0],position[1],position[2]));
    return Vec3(myVec3Float[0], myVec3Float[1], myVec3Float[2]);
}

// Functions which were moved from GridPoint to DensityMap for memory savings

void DensityMap::initializeGradient(GridPoint & gridPoint){
		setddyPositiveXGradient(gridPoint,  0);
		setddzPositiveXGradient(gridPoint,  0);
		setddxPositiveYGradient(gridPoint,  0);
		setddzPositiveYGradient(gridPoint,  0);
		setddxPositiveZGradient(gridPoint,  0);
		setddyPositiveZGradient(gridPoint,  0);
                setFirstQuadrantGradient(gridPoint,  fVec3(0,0,0));
	}

void DensityMap::initialize(GridPoint & gridPoint){
		initializeGradient(gridPoint);
                setFirstQuadrantGradient(gridPoint,fVec3(0));
		gridPoint.density = 0; 
		gridPoint.position = fVec3(0);	
                //cout<<__FILE__<<":"<<__LINE__<<" Set firstQuadrantGradient =  "<<fetchFirstQuadrantGradient(gridPoint)<<"  position = "<<gridPoint.position <<" density = "<<gridPoint.density <<endl;
	}

void DensityMap::validatePosition(GridPoint & gridPoint, fVec3 myPosition) const{
               
          	//cout<<__FILE__<<":"<<__LINE__<<" about to validate position = "<<myPosition<<endl;
                //ValidateVec3( myPosition);
	}

void DensityMap::validateDensity(GridPoint & gridPoint, double myDensity) const{
		ValidateReal(myDensity);
	
}

void DensityMap::validate(GridPoint & gridPoint) const{
		validatePosition(gridPoint,gridPoint.position);
		validateDensity(gridPoint,gridPoint.density);
		// write this code later.
	}

void DensityMap::setDensity(GridPoint & gridPoint,Real myDensity)	{
    validateDensity(gridPoint,myDensity); 
    gridPoint.density = (float)myDensity;	
    //cout<<__FILE__<<":"<<__LINE__<<" Just set density to "<<gridPoint.density<<" for grid point at position "<<gridPoint.position<<endl;
}

void DensityMap::setPosition(GridPoint & gridPoint, fVec3 myPosition)	{
 		validatePosition(gridPoint,myPosition);
		gridPoint.position = myPosition;	
	}

Quadrant DensityMap::calcQuadrant(GridPoint & gridPoint, fVec3 queryPosition) const			{
		Quadrant tempQuadrant;    
		if (queryPosition[0] < gridPoint.position[0]) {tempQuadrant.positiveX = false;} else { tempQuadrant.positiveX = true;}		
		if (queryPosition[1] < gridPoint.position[1]) {tempQuadrant.positiveY = false;} else { tempQuadrant.positiveY = true;}
		if (queryPosition[2] < gridPoint.position[2]) {tempQuadrant.positiveZ = false;} else { tempQuadrant.positiveZ = true;}
		return tempQuadrant;
	}

fVec3 DensityMap::fetchGradient(GridPoint & gridPoint,fVec3 queryPosition) const {
		Quadrant tempQuadrant = calcQuadrant(gridPoint,queryPosition);
		fVec3 myGradient;
          	ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" this function is obsolete, fetchGradient doesn't work anymore"<<endl;
        	ErrorManager::instance.treatError();
		return myGradient;

	}

fVec3 DensityMap::fetchFirstQuadrantGradient(GridPoint & gridPoint) const {
                return gridPoint.firstQuadrantGradient;
}



double DensityMap::getDensity(GridPoint & gridPoint) const	{return gridPoint.density ;	}
	
double DensityMap::getDensity( GridPoint & gridPoint, fVec3 queryPosition) const	{
                fVec3 myVectorToGridMap = queryPosition - gridPoint.position;
                fVec3 myGradient = calcInterpolatedFirstQuadrantGradient(gridPoint,  queryPosition); 
                return (gridPoint.density + (myGradient[0]* myVectorToGridMap[0])  + (myGradient[1]* myVectorToGridMap[1]) + (myGradient[2]* myVectorToGridMap[2]));
}




void DensityMap::setPositiveXGradient(GridPoint & gridPoint,Real myPositiveXGradient) {  gridPoint.firstQuadrantGradient[0] = myPositiveXGradient; }
void DensityMap::setddyPositiveXGradient(GridPoint & gridPoint,Real value) { gridPoint.ddyPositiveXGradient = (float)value; }
void DensityMap::setddzPositiveXGradient(GridPoint & gridPoint,Real value) { gridPoint.ddzPositiveXGradient = (float)value; }
void DensityMap::setPositiveYGradient(GridPoint & gridPoint, Real myPositiveYGradient) { gridPoint.firstQuadrantGradient[1] = myPositiveYGradient; }
void DensityMap::setddxPositiveYGradient(GridPoint & gridPoint,Real value) { gridPoint.ddxPositiveYGradient = (float)value; }
void DensityMap::setddzPositiveYGradient(GridPoint & gridPoint,Real value) { gridPoint.ddzPositiveYGradient = (float)value; }
void DensityMap::setPositiveZGradient(GridPoint & gridPoint,Real myPositiveZGradient) {  gridPoint.firstQuadrantGradient[2] = myPositiveZGradient; }
void DensityMap::setddxPositiveZGradient(GridPoint & gridPoint,Real value) { gridPoint.ddxPositiveZGradient = (float)value; }
void DensityMap::setddyPositiveZGradient(GridPoint & gridPoint,Real value) { gridPoint.ddyPositiveZGradient = (float)value; }
void DensityMap::setNegativeXGradient(GridPoint & gridPoint,Real myNegativeXGradient) { 
        //ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" this function is only for calculating the gradient in the first quadrant"<<endl;
        ErrorManager::instance.treatError();}
void DensityMap::setNegativeYGradient(GridPoint & gridPoint,Real myNegativeYGradient) { 
        //ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" this function is only for calculating the gradient in the first quadrant"<<endl;
        ErrorManager::instance.treatError();}
void DensityMap::setNegativeZGradient(GridPoint & gridPoint,Real myNegativeZGradient) { 
        //ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" this function is only for calculating the gradient in the first quadrant"<<endl;
        ErrorManager::instance.treatError();}
void DensityMap::setFirstQuadrantGradient(GridPoint & gridPoint,fVec3 gradient){
        //cout<<__FILE__<<":"<<__LINE__<<" Setting firstQuadrantGradient to "<<gradient<<" for grid point at position : "<<gridPoint.position <<endl;
        gridPoint.firstQuadrantGradient = gradient;
    }
fVec3 DensityMap::calcInterpolatedFirstQuadrantGradient(GridPoint & gridPoint,fVec3 queryPosition) const {
    fVec3 dxdydz = queryPosition - gridPoint.position; // the first term is the query position, the second term is the grid point position in cartesian space
    if ((dxdydz[0] < 0) || (dxdydz[1] <0 ) || (dxdydz[2] < 0)) { // see if we can make this trap unnecessary implicitly

        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" this function is only for calculating the gradient in the first quadrant"<<endl;
        ErrorManager::instance.treatError();
    }
    fVec3 myGradient;
    myGradient[0] = fetchFirstQuadrantGradient(gridPoint)[0] + gridPoint.ddyPositiveXGradient*dxdydz[1] + gridPoint.ddzPositiveXGradient*dxdydz[2];
    myGradient[1] = fetchFirstQuadrantGradient(gridPoint)[1] + gridPoint.ddxPositiveYGradient*dxdydz[0] + gridPoint.ddzPositiveYGradient*dxdydz[2];
    myGradient[2] = fetchFirstQuadrantGradient(gridPoint)[2] + gridPoint.ddxPositiveZGradient*dxdydz[0] + gridPoint.ddyPositiveZGradient*dxdydz[1];
    return myGradient;
}

/*fVec3 DensityMap::getFirstQuadrantGradient(GridPoint & gridPoint){
    return fetchFirstQuadrantGradient(gridPoint);
}*/

// Expects OpenDX dx density maps with Angstrom as spatial unit.
// Adapted from BioSpring http://sourceforge.net/projects/biospring/
void DensityMap::loadParametersAndDensity_OpenDX(const String densityFileName)
{   
    FILE * fdx = NULL;
    int totalsize = 0;
    char inbuf[LINESIZE];
    unsigned  gx = 0;
    unsigned gy = 0;
    unsigned gz = 0;

    unsigned sizei=0,sizej=0,sizek=0;
    double offsetx=0.0,offsety=0.0,offsetz=0.0;


    double grid[3]={0.0};
    double delta[3][3];
    
    if ((fdx=fopen(densityFileName,"r"))==0)
    {
        ErrorManager::instance << "Can't open file " <<  densityFileName << endl;
        ErrorManager::instance.treatError();
    }

    // Jumping comments in the .dx file
    do
    {
        if(dxGets(inbuf, LINESIZE, fdx) == NULL)
        {
            ErrorManager::instance << "Error while jumping comments in " <<  densityFileName << endl;
            ErrorManager::instance.treatError();
        }
    }
    while (inbuf[0]=='#');


    //Get number of grid points
    if (sscanf(inbuf, "object 1 class gridpositions counts %d %d %d", &sizei,&sizej,&sizek)!=3)
    {
            ErrorManager::instance << "Error while reading dimensions in " <<  densityFileName << endl;
            ErrorManager::instance.treatError();
    }
    else
    {
        totalNumGridX = sizei;
        totalNumGridY = sizej;
        totalNumGridZ = sizek;
    }
    //Grid origin
    if (dxGets(inbuf, LINESIZE, fdx) == NULL) 
    {
        fprintf(stderr, "dxGets -> Erreur de lecture origine de la grille!\n");
    }
    else if (sscanf(inbuf, "origin %lf %lf %lf", &(offsetx), &(offsety), &(offsetz))!=3)
    {   
        ErrorManager::instance << "Error while reading origin coordinates in " <<  densityFileName << endl;
        ErrorManager::instance.treatError();
    }
    else
    {
        minX = offsetx/10.0;
        minY = offsety/10.0;
        minZ = offsetz/10.0;
        // cout.precision(5);
        cout << "Grid Origin " << minX << " " << minY << " " << minZ << endl;
    }

    // Grid size
    // dimension X
    if (dxGets(inbuf, LINESIZE, fdx) == NULL)
    {
        ErrorManager::instance << "Error while reading grid size X in " <<  densityFileName << endl;
        ErrorManager::instance.treatError();
    }
    //cout<<"inbuff:"<<inbuf<<endl;
    else if (sscanf(inbuf, "delta %lf %lf %lf", &(delta[0][0]),&(delta[0][1]),&(delta[0][2]))!=3)
    {
        ErrorManager::instance << "Error while reading grid size X in " <<  densityFileName << endl;
        ErrorManager::instance.treatError();
    }
    // dimension Y
    else if (dxGets(inbuf, LINESIZE, fdx) == NULL)
    {
        ErrorManager::instance << "Error while reading grid size Y in " <<  densityFileName << endl;
        ErrorManager::instance.treatError();
    }
    else if (sscanf(inbuf, "delta %lf %lf %lf", &(delta[1][0]), &(delta[1][1]), &(delta[1][2]))!=3)
    {
        ErrorManager::instance << "Error while reading grid size Y in " <<  densityFileName << endl;
        ErrorManager::instance.treatError();
    }
    // dimension Z
    else if (dxGets(inbuf, LINESIZE, fdx) == NULL)
    {
        ErrorManager::instance << "Error while reading grid size Z in " <<  densityFileName << endl;
        ErrorManager::instance.treatError();
    }

    else if (sscanf(inbuf, "delta %lf %lf %lf", &(delta[2][0]), &(delta[2][1]), &(delta[2][2]))!=3)
    {
        ErrorManager::instance << "Error while reading grid size Z in " <<  densityFileName << endl;
        ErrorManager::instance.treatError();
    }
    else 
    {
        gridXSpacing = delta[0][0]/10.0;
        gridYSpacing = delta[1][1]/10.0;
        gridZSpacing = delta[2][2]/10.0;
    }
    //cout<<"delta "<<_deltax<<" "<<_deltay<<" "<<_deltaz<<endl;
    // skipping "object 2 class gridconnections counts xn yn zn"
    if (dxGets(inbuf, LINESIZE, fdx) == NULL)
    {
        ErrorManager::instance << "Error while reading " <<  densityFileName << endl;
        ErrorManager::instance.treatError();
    }
    // skipping "object 3 class array type double rank 0 items xn*yn*zn [binary] data follows
    if (dxGets(inbuf, LINESIZE, fdx) == NULL)
    {   
        ErrorManager::instance << "Error while reading " <<  densityFileName << endl;
        ErrorManager::instance.treatError();
    }


    totalsize = sizei*sizej * sizek;
    cout << "Grid Spacing " << gridXSpacing << " " << gridYSpacing << " " << gridZSpacing << endl;

    maxX = minX + ( (totalNumGridX-1) * gridXSpacing);
    maxY = minY + ( (totalNumGridY-1) * gridYSpacing);
    maxZ = minZ + ( (totalNumGridZ-1) * gridZSpacing);
    initializeArrayOfGridPoints();
    
    float unityconvert=1.0;
    for (int count =0; count<(totalsize/3); count++)
    {
        if (dxGets(inbuf,LINESIZE,fdx)==NULL)
        {
            ErrorManager::instance << "Error while reading value in " <<  densityFileName << endl;
            ErrorManager::instance.treatError();
        }
        if (sscanf(inbuf,"%lf %lf %lf", &grid[0], &grid[1], &grid[2]) != 3)
        {
            ErrorManager::instance << "Error while reading value in " <<  densityFileName << endl;
            ErrorManager::instance.treatError();
        }
        
        for (int i=0;i<3;i++)
        {
            GridIndices centralGridIndices(gx, gy, gz);
            GridPoint & centralGridPoint = updGridPoint(centralGridIndices);
            setDensity(centralGridPoint, grid[i]*unityconvert);
            // cout<<"Potential "<<gx<<" "<<gy<<" "<<gz<<" "<<grid[i]<<endl;
            gz++;
            if(gz>=sizek)
            {
                gz=0;
                gy++;
                if(gy>=sizej)
                {
                    gy=0;
                    gx++;
                }
            }
        }
    }

    if((totalsize%3)!=0)
    {   
        //printf("dans total\n");
        if(dxGets(inbuf,LINESIZE,fdx)==NULL)
        {
            ErrorManager::instance << "Error while reading values in " <<  densityFileName << endl;
            ErrorManager::instance.treatError();
        }
        int count = sscanf(inbuf,"%lf %lf %lf", &grid[0], &grid[1], &grid[2]);
        if(count!=(totalsize%3))
        {
            ErrorManager::instance << "Error, too many data point in " <<  densityFileName << endl;
            ErrorManager::instance.treatError();
        }
        for(int i=0;i<count;i++)
        {
            GridIndices centralGridIndices(gx, gy, gz);
            GridPoint & centralGridPoint = updGridPoint(centralGridIndices);
            setDensity(centralGridPoint, grid[i]*unityconvert);
            gz++;
        }
    }
    fclose(fdx);  
}

char * DensityMap::dxGets( char *s, int n, FILE *stream)
{
    char *resdxGets;

    if (feof(stream)) 
    {
        ErrorManager::instance << __FILE__<<":"<<__FILE__<< " Error: unexpected End of File" << endl;
        ErrorManager::instance.treatError();
    }
    else if (ferror(stream)) 
    {
        ErrorManager::instance << __FILE__<<":"<<__FILE__<<  " Error while reading dx file" << endl;
        ErrorManager::instance.treatError();
    }
    else    
    {
        resdxGets = fgets(s,n,stream);
        if (resdxGets == NULL)
        {
            ErrorManager::instance << __FILE__<<":"<<__FILE__<< "Error while reading dx file " <<  endl;
            ErrorManager::instance.treatError();
        }
    }
    return resdxGets;
}