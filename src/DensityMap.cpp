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
#include "UnitCellParameters.h"
#include <time.h>       /* time */

//================================================ Include Gemmi headers
#include <gemmi/ccp4.hpp>
#include <gemmi/gz.hpp>

using namespace std;
using namespace SimTK;


GridIndices::GridIndices( int myXIndex,  int myYIndex,  int myZIndex) {
		xGridPoint = ValidateInt(myXIndex);
		yGridPoint = ValidateInt(myYIndex);
		zGridPoint = ValidateInt(myZIndex);
	}
int GridIndices::getXGridIndex () const {return xGridPoint;}
int GridIndices::getYGridIndex () const {return yGridPoint;}
int GridIndices::getZGridIndex () const {return zGridPoint;}

DensityMap::DensityMap(){
    unitCellParameters.setDefaultParameters(); 
    //initializeMap();
};

DensityMap::~DensityMap(){
    ArrayOfGridPoints.clear();
};

void DensityMap::initializeMap() {
    unitCellParameters.setDefaultParameters(); 
    setNoiseTemperature(00.);
    setNoiseScale(.0);
    /*
    unitCellNumGridX = 0;
    unitCellNumGridY = 0;
    unitCellNumGridZ = 0;
    totalNumGridX = 0;
    totalNumGridY = 0;
    totalNumGridZ = 0;
    unitCellParameters.geta()*unitCellParameters.getaMin() = 0;
    minY = 0;
    minZ = 0;
    maxX = 0;
    maxY = 0;
    maxY = 0;
    unitCellParameters.geta() = 0;
    gridYSpacing = 0;
    gridZSpacing = 0;*/	
};

void DensityMap::validateGridParameters() {
    unitCellParameters.validate();
    //cout<<__FILE__<<":"<<__LINE__<< " Checking maxX minX  gridXSpacing totalNumGridX (in nm, nm, nm, unitless): "<<maxX<<" , "<<minX<<" , "<<gridXSpacing<<" , "<<totalNumGridX<<std::endl;
    //cout<<__FILE__<<":"<<__LINE__<< " Checking maxY minY  gridYSpacing totalNumGridY (in nm, nm, nm, unitless): "<<maxY<<" , "<<minY<<" , "<<gridYSpacing<<" , "<<totalNumGridY<<std::endl;
    //cout<<__FILE__<<":"<<__LINE__<< " Checking maxZ minZ  gridZSpacing totalNumGridZ (in nm, nm, nm, unitless): "<<maxZ<<" , "<<minZ<<" , "<<gridZSpacing<<" , "<<totalNumGridZ<<std::endl;
    /*
    if (
	( abs(((maxX-minX)/gridXSpacing +1) - totalNumGridX) > 1E-7  ) ||
        ( abs(((maxY-minY)/gridYSpacing +1) - totalNumGridY) > 1E-7  ) ||
        ( abs(((maxZ-minZ)/gridZSpacing +1) - totalNumGridZ) > 1E-7  ) 
       ) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL,  " There is a problem with the max- , min- (XYZ) or totalNumGrid (XYZ) or grid (XYZ) Spacing parameters "<< ((maxX-minX)/gridXSpacing +1) - totalNumGridX <<endl;	
    }*/
}
		
bool        DensityMap::hasGridPoint(GridIndices myGridIndices){

        //cout<<__FILE__<<":"<<__LINE__<< " ArrayOfGridPoints[0].size() = "  <<ArrayOfGridPoints[0].size() <<endl;
        //cout<<__FILE__<<":"<<__LINE__<< " ArrayOfGridPoints[0][0].size() = "  <<ArrayOfGridPoints[0][0].size() <<endl;
	bool myHasGridPoint;
	if (myGridIndices.getZGridIndex() >= unitCellParameters.getNc()) 
        {
		//cout<<__FILE__<<":"<<__LINE__<< " The Z-index of "<<myGridIndices.getZGridIndex()<<" exceeds the Z-dimension, "<< unitCellParameters.getNc() <<" of the grid map. No force applied."<<endl;	
		myHasGridPoint = false;
	}
	else if (myGridIndices.getZGridIndex() <  0 ) {
		//cout<<__FILE__<<":"<<__LINE__<< " The Z-index of "<<myGridIndices.getZGridIndex()<<" is less than zero. No force applied."<<endl; 
		myHasGridPoint = false;
	}
	else if (myGridIndices.getYGridIndex() >= unitCellParameters.getNb()) {
		//cout<<__FILE__<<":"<<__LINE__<< " The Y-index of "<<myGridIndices.getYGridIndex()<<" exceeds the Y-dimension, "<<ArrayOfGridPoints[myGridIndices.getZGridIndex()].size() <<" of the grid map. No force applied."<<endl;	
		myHasGridPoint = false;
	}
	else if (myGridIndices.getYGridIndex() <  0 ) {
		//cout<<__FILE__<<":"<<__LINE__<< " The Y-index of "<<myGridIndices.getYGridIndex()<<" is less than zero. No force applied."<<endl; 
		myHasGridPoint = false;
	}
	else if (myGridIndices.getXGridIndex() >= unitCellParameters.getNa()) {
		//cout<<__FILE__<<":"<<__LINE__<< " The X-index of "<<myGridIndices.getXGridIndex()<<" exceeds the X-dimension, "<<unitCellParameters.getNa() <<" of the grid map. No force applied."<<endl;	
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
			MMBLOG_FILE_FUNC_LINE(CRITICAL, "No nearby grid point.  The point you requested "<< myGridIndices.getZGridIndex()<<" , "<<  myGridIndices.getYGridIndex()<<" , "<< myGridIndices.getXGridIndex() << " is off the density map."<<endl);
		} else {
			GridPoint & myGridPoint = ArrayOfGridPoints[myGridIndices.getZGridIndex()][myGridIndices.getYGridIndex()][myGridIndices.getXGridIndex()] ;
			//myGridPoint.validate();
			return myGridPoint;
		}
  		}
   
GridPoint DensityMap::getGridPoint(GridIndices myGridIndices) const  {
		//return & ArrayOfGridPoints[0][0][0];
		return ArrayOfGridPoints[myGridIndices.getZGridIndex()][myGridIndices.getYGridIndex()][myGridIndices.getXGridIndex()] ;
  		}
   
void	    DensityMap::validateGridPoint(GridIndices myGridIndices){
			if (! hasGridPoint(myGridIndices)) {
				MMBLOG_FILE_FUNC_LINE(CRITICAL, "No nearby grid point.  The point you requested is off the density map."<<endl);
			} else validate(updGridPoint(myGridIndices));

}


/*const bool DensityMap::hasNearbyGridIndices(Vec3 position)
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
		
*/
GridIndices DensityMap::calcNearestGridIndices(const Vec3 &position)
		{
                        iVec3 tempIndexVector = {0,0,0};
                        tempIndexVector = unitCellParameters.convertCartesianVectorToNearestIndexVector(position);
                        
			//int tempXIndex 	= int (((position[0]-minX) )/gridXSpacing);
			//int tempYIndex 	= int (((position[1]-minY) )/gridYSpacing);
			//int tempZIndex 	= int (((position[2]-minZ) )/gridZSpacing);
			//int tempYIndex 	= 	int((position[1] + gridYSpacing/2)/gridYSpacing);
			//int tempZIndex 	= 	int((position[2] + gridZSpacing/2)/gridZSpacing);
			GridIndices tempGridIndices(tempIndexVector[0], tempIndexVector[1], tempIndexVector[2]); //tempXIndex, tempYIndex, tempZIndex);
			return tempGridIndices;
		}
		
GridIndices DensityMap::calcLowerLeftGridIndices(const Vec3 &position)
                {
                        iVec3 tempIndexVector = unitCellParameters.convertCartesianVectorToLowerIndexVector(position); //convertFractionalVectorToLowerIndexVector(position);
                        //int tempXIndex  = int (floor(((position[0]-minX) )/gridXSpacing));
                        //int tempYIndex  = int (floor(((position[1]-minY) )/gridYSpacing));
                        //int tempZIndex  = int (floor(((position[2]-minZ) )/gridZSpacing));

                        return GridIndices(tempIndexVector[0], tempIndexVector[1], tempIndexVector[2]); /// tempXIndex, tempYIndex, tempZIndex);

                }
		
GridPoint DensityMap::getGridPoint(const Vec3 &myPosition)  {
	return getGridPoint(calcNearestGridIndices( myPosition));				
}

GridPoint & DensityMap::updGridPoint(const Vec3 &myPosition)  {
	return updGridPoint(calcNearestGridIndices( myPosition));				
}

double DensityMap::getDensity(const Vec3 &myPosition) {
        //MMBLOG_FILE_FUNC_LINE("Taking myPosition = "<<myPosition<<std::endl;
        Vec3 myFractionalVector = unitCellParameters.convertCartesianVectorToFractionalVector(myPosition);
        //MMBLOG_FILE_FUNC_LINE( myFractionalVector[0] <<", "<<  myFractionalVector[1]   <<", "<< myFractionalVector[2]   <<" ..preceding should be unitCellParameters.convertCartesianVectorToFractionalVector(myPosition)"<<std::endl;
        if (!(unitCellParameters.fractionalVectorIsInsideMapBoundaries(myFractionalVector)))
	{
		return 0.0; //return zero density
	} else {
                GridIndices myLowerLeftGridIndices = calcLowerLeftGridIndices(myPosition);
                    MMBLOG_FILE_FUNC_LINE(DEBUG, endl);
		    return getDensity(  updGridPoint(myLowerLeftGridIndices), myPosition); // need to verify that this interpolated density is reasonable
	} 
}

//const double DensityMap::getDensity(SimTK::Vec3 myPosition) {
//    return getDensity(Vec3(myPosition[0], myPosition[1],myPosition[2] ));
//}
void DensityMap::initializeArrayOfGridPoints(){
        //validateGridParameters();
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        unitCellParameters.validate();
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        initializeVectorOfAmplitudeAndRandomPhases();
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        Vec3 tempPosition(0,0,0);
        ArrayOfGridPoints.resize(unitCellParameters.getNc());
        for ( int zIndex = 0; zIndex < unitCellParameters.getNc(); zIndex++) {
        	ArrayOfGridPoints[zIndex].resize(unitCellParameters.getNb());
        	for ( int yIndex = 0; yIndex < unitCellParameters.getNb(); yIndex++) {
        		ArrayOfGridPoints[zIndex][yIndex].resize(unitCellParameters.getNa());
        	}}
        for ( int xIndex = 0; xIndex < unitCellParameters.getNa(); xIndex++) {
        	for ( int yIndex = 0; yIndex < unitCellParameters.getNb(); yIndex++) {
        		for ( int zIndex = 0; zIndex < unitCellParameters.getNc(); zIndex++) {
        			GridPoint tempGridPoint;// = updGridPoint(GridIndices(xIndex, yIndex, zIndex)) ;
        			initialize(tempGridPoint);
                                //std::cout<<" We will need a function here called convertFractionalVectorToCartesianVector .. I think? or do we?";
                                //exit(1);
        			tempPosition = Vec3(xIndex * unitCellParameters.geta() + unitCellParameters.geta()*unitCellParameters.getaMin(), yIndex * unitCellParameters.getb() + unitCellParameters.getb()*unitCellParameters.getbMin(),zIndex * unitCellParameters.getc() + unitCellParameters.getc()*unitCellParameters.getcMin());
        			setPosition(tempGridPoint,tempPosition);
        			validate(tempGridPoint);
        			updGridPoint(GridIndices(xIndex, yIndex, zIndex)) = tempGridPoint;
        			validate(updGridPoint(GridIndices(xIndex, yIndex, zIndex))); // just being paranoid
        			//ArrayOfGridPoints[zIndex, yIndex, xIndex] = tempGridPoint;			
        }}}
        if (ArrayOfGridPoints.size() != unitCellParameters.getNc()) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Wrong number of grid points in Z direction! Found :"<< ArrayOfGridPoints.size()<<" expected : " << unitCellParameters.getNc()<<","<< unitCellParameters.getNc()<<endl);
        }
        for ( int zIndex = 0; zIndex < unitCellParameters.getNc(); zIndex++) {
        	if (ArrayOfGridPoints[zIndex].size() != unitCellParameters.getNb()) {
                   MMBLOG_FILE_FUNC_LINE(CRITICAL, "Wrong number of grid points in Y direction! Found :"<< ArrayOfGridPoints[zIndex].size()<<" expected : " << unitCellParameters.getNb()<<","<< unitCellParameters.getNc()<<endl);
                }
        	for ( int yIndex = 0; yIndex < unitCellParameters.getNb(); yIndex++) {
        	    if (ArrayOfGridPoints[zIndex][yIndex].size() != (unitCellParameters.getNa()) ) {
        		             MMBLOG_FILE_FUNC_LINE(CRITICAL, "Wrong number of grid points in X direction! Found :"<< ArrayOfGridPoints[zIndex][yIndex].size()<<" expected : " << unitCellParameters.getNa()<<endl);
                }
        	}
        } // of for zIndex
}

double plancksLaw(double temperature, double frequency){
    // This computes number density
    // frequency is "nu" .. not that it matters so much here
    // lambda is wavelength.
    // pi and c are both taken to be unity
    // nu = 1/lambda
    //if (temperature >= 1000000){ 
    //    MMBLOG_FILE_FUNC_LINE(" You have specified a very high temperature : "<<temperature<<" .. returning white noise spectrum"<<std::endl;  
    //    return 1.;
    //}

    return (1 / (exp(frequency/temperature) - 1));
}



void DensityMap::resizeVectorOfAmplitudeAndRandomPhases(){ 
        MMBLOG_FILE_FUNC_LINE(INFO, "resizing Z component vectorOfAmplitudeFrequencyAndRandomPhases to "<<(unitCellParameters.calcMaxFrequencyDoublingsZ () +1)<<endl);
        vectorOfAmplitudeFrequencyAndRandomPhases.resize(unitCellParameters.calcMaxFrequencyDoublingsZ () +1);                                         
        MMBLOG_FILE_FUNC_LINE(DEBUG, endl);
        for ( int zIndex = 0; zIndex <= unitCellParameters.calcMaxFrequencyDoublingsZ(); zIndex++) {
                MMBLOG_FILE_FUNC_LINE(DEBUG, endl);
                MMBLOG_FILE_FUNC_LINE(DEBUG, "resizing  vectorOfAmplitudeFrequencyAndRandomPhases["<<zIndex <<"] to "<<(unitCellParameters.calcMaxFrequencyDoublingsZ () +1)<<endl);
        	vectorOfAmplitudeFrequencyAndRandomPhases[zIndex].resize(unitCellParameters.calcMaxFrequencyDoublingsY () +1);
        	for ( int yIndex = 0; yIndex <=  unitCellParameters.calcMaxFrequencyDoublingsY(); yIndex++) {
                    //MMBLOG_FILE_FUNC_LINE(std::endl;
                    //MMBLOG_FILE_FUNC_LINE(" resizing  vectorOfAmplitudeFrequencyAndRandomPhases["<<zIndex <<"]["<<yIndex<< "] to "<<(unitCellParameters.calcMaxFrequencyDoublingsX () +1)<<std::endl;
                    vectorOfAmplitudeFrequencyAndRandomPhases[zIndex][yIndex].resize(unitCellParameters.calcMaxFrequencyDoublingsX () +1);
        	}}
        MMBLOG_FILE_FUNC_LINE(INFO, endl);

        // Now make sure size was set correctly:
        if (vectorOfAmplitudeFrequencyAndRandomPhases.size() != (unitCellParameters.calcMaxFrequencyDoublingsZ()+1)) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Wrong number of grid points in Z direction! Found :"<< vectorOfAmplitudeFrequencyAndRandomPhases.size()<<" expected : " << unitCellParameters.calcMaxFrequencyDoublingsZ()+1             <<endl);
        }
        for ( int zIndex = 0; zIndex < unitCellParameters.calcMaxFrequencyDoublingsZ(); zIndex++) {
        	if (vectorOfAmplitudeFrequencyAndRandomPhases[zIndex].size() != (unitCellParameters.calcMaxFrequencyDoublingsY()+1) ) {
                   MMBLOG_FILE_FUNC_LINE(CRITICAL, "Wrong number of grid points in Y direction! Found :"<< vectorOfAmplitudeFrequencyAndRandomPhases[zIndex].size()<<" expected : " << unitCellParameters.calcMaxFrequencyDoublingsY() +1            <<endl);
                }
	        for ( int yIndex = 0; yIndex < unitCellParameters.calcMaxFrequencyDoublingsY(); yIndex++) {
        	    if (vectorOfAmplitudeFrequencyAndRandomPhases[zIndex][yIndex].size() != (unitCellParameters.calcMaxFrequencyDoublingsX()+1 )) {
                             MMBLOG_FILE_FUNC_LINE(DEBUG, endl);
        		     MMBLOG_FILE_FUNC_LINE(CRITICAL, "Check 1. Wrong number of grid points in X direction! Found :"<< vectorOfAmplitudeFrequencyAndRandomPhases[zIndex][yIndex].size()<<" expected : " <<unitCellParameters.calcMaxFrequencyDoublingsX()+1 <<endl);
                    }
        	}
        } // of for zIndex

}

// Compute autocorrelation functions for density maps
void DensityMap::densityAutocorrelation(const bool computeNoiseAutocorrelation, const bool computeDensityAutocorrelation ) const {
    struct autoCorrelationStruct{ 
        //double distanceSquared;
        long double sumOfObservations;
        int    numOfObservations;
    };
    autoCorrelationStruct myAutoCorrelationStruct ;
    std::map <double ,autoCorrelationStruct> autoCorrelationMap;
    std::map <double ,autoCorrelationStruct>::iterator autoCorrelationMapIterator;
    if (!(computeNoiseAutocorrelation ^computeDensityAutocorrelation)){ // ^is XOR. So we demand that the user specify exactly one of computeNoiseAutocorrelation or computeDensityAutocorrelation, but not both.
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have specified computeNoiseAutocorrelation = "<<computeNoiseAutocorrelation<<" and computeDensityAutocorrelation = "<< computeDensityAutocorrelation << " .. you must specify exactly one of these." <<endl);
    }
    double distanceSquared = 0.;
    double xSquared = 0.;
    double ySquared = 0.;
    double zSquared = 0.;
    double ySquaredplusXsquared = 0.;
    
    double correlation = 0.;
    int xIndex2Corrected = 0;
    int yIndex2Corrected = 0;
    int zIndex2Corrected = 0;

    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    for ( int xIndex1 = 0; xIndex1 < unitCellParameters.getNa();  xIndex1++){
        for ( int yIndex1 = 0; yIndex1 < unitCellParameters.getNb();  yIndex1++){
            for ( int zIndex1 = 0; zIndex1 < unitCellParameters.getNc();  zIndex1++)
            {
                //for ( int xIndex2 = 0; xIndex2 < unitCellParameters.getNa();  xIndex2++)
                int xIndex2 = xIndex1;
                {
                    if (xIndex2 < xIndex1) { // this is the case that Index2 is in the next unit cell, in the positive direction.
                        xIndex2Corrected = xIndex2 + unitCellParameters.getNa() - 1;
                    } else xIndex2Corrected = xIndex2;

                    xSquared        = (xIndex2Corrected-xIndex1)*(xIndex2Corrected-xIndex1); // x^2
                    xSquared       *= unitCellParameters.geta()*unitCellParameters.geta(); // convert from index to distance. 
                    //for ( int yIndex2 = 0; yIndex2 < unitCellParameters.getNb();  yIndex2++)
                    int yIndex2 = yIndex1;
                    {
                        if (yIndex2 < yIndex1) { // this is the case that Index2 is in the next unit cell, in the positive direction.
                            yIndex2Corrected = yIndex2 + unitCellParameters.getNb() - 1;
                        } else yIndex2Corrected = yIndex2;
                        ySquaredplusXsquared    = xSquared + (yIndex2Corrected-yIndex1)*(yIndex2Corrected-yIndex1)*unitCellParameters.getb()*unitCellParameters.getb(); // add y^2
                        //ySquaredplusXsquared   *= unitCellParameters.getb()*unitCellParameters.getb();
                        for ( int zIndex2 = 0; zIndex2 < unitCellParameters.getNc();  zIndex2++){
                            if (zIndex2 < zIndex1) { // this is the case that Index2 is in the next unit cell, in the positive direction.
                                zIndex2Corrected = zIndex2 + unitCellParameters.getNc() - 1;
                            } else zIndex2Corrected = zIndex2;
                            distanceSquared = ySquaredplusXsquared + (zIndex2Corrected-zIndex1)*(zIndex2Corrected-zIndex1)*unitCellParameters.getc()*unitCellParameters.getc(); // add z^2
                            if (computeDensityAutocorrelation){
                                correlation = getGridPoint(GridIndices(xIndex1,yIndex1,zIndex1)).noiseFreeDensity * getGridPoint(GridIndices(xIndex2,yIndex2,zIndex2)).noiseFreeDensity ; 

                            } else if (computeNoiseAutocorrelation){
                                correlation = getGridPoint(GridIndices(xIndex1,yIndex1,zIndex1)).noise * getGridPoint(GridIndices(xIndex2,yIndex2,zIndex2)).noise ; 
                            }
                            autoCorrelationMapIterator = autoCorrelationMap.find(distanceSquared); 
                            if (autoCorrelationMapIterator != autoCorrelationMap.end()){autoCorrelationMapIterator->second.sumOfObservations += correlation;
                                autoCorrelationMapIterator->second.numOfObservations ++;  
                                //MMBLOG_FILE_FUNC_LINE("distanceSquared, sumCorrelation, numCorrelationObservations, correlation, "<<autoCorrelationMapIterator->first<<","<<autoCorrelationMapIterator->second.sumOfObservations <<","<<autoCorrelationMapIterator->second.numOfObservations<<","<<autoCorrelationMapIterator->second.sumOfObservations / autoCorrelationMapIterator->second.numOfObservations<<std::endl;
             
                            }
                            else {
                                myAutoCorrelationStruct.sumOfObservations = correlation;
                                myAutoCorrelationStruct.numOfObservations = 1;
                                autoCorrelationMap.insert(std::make_pair(distanceSquared,myAutoCorrelationStruct));
                                //MMBLOG_FILE_FUNC_LINE("distanceSquared, sumCorrelation, numCorrelationObservations, correlation, "<<distanceSquared<<","<<myAutoCorrelationStruct.sumOfObservations  <<","<<myAutoCorrelationStruct.numOfObservations<<std::endl;
                            }

                            autoCorrelationMap.insert(std::make_pair(distanceSquared,myAutoCorrelationStruct));
                            //MMBLOG_FILE_FUNC_LINE( " computeDensityAutocorrelation = " << computeDensityAutocorrelation <<" computeNoiseAutocorrelation = "<<computeNoiseAutocorrelation <<",  distanceSquared, correlation = ,"<<distanceSquared<<","<<correlation<<std::endl;
    }}}}}}
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    MMBLOG_FILE_FUNC_LINE(INFO, std::string(40,'*')<<endl);
    MMBLOG_FILE_FUNC_LINE(INFO, "Starting autocorrelation section. "<<endl);
    MMBLOG_FILE_FUNC_LINE(INFO, std::string(40,'*')<<endl);
    MMBLOG_FILE_FUNC_LINE(INFO, "distanceSquared, sumCorrelation, numCorrelationObservations, correlation "<<endl);
    for(autoCorrelationMapIterator = autoCorrelationMap.begin(); autoCorrelationMapIterator != autoCorrelationMap.end(); autoCorrelationMapIterator++){
        MMBLOG_FILE_FUNC_LINE(INFO, autoCorrelationMapIterator->first<<","<<autoCorrelationMapIterator->second.sumOfObservations <<","<<autoCorrelationMapIterator->second.numOfObservations<<","<<autoCorrelationMapIterator->second.sumOfObservations / autoCorrelationMapIterator->second.numOfObservations<<endl);
    }
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
}

void DensityMap::normalizeNoiseMap(const double totalNoiseEverywhere){
    // this part would set the average noise per grid point to unity:
    double normalization = ( unitCellParameters.getNa() *  unitCellParameters.getNb() * unitCellParameters.getNc()) / totalNoiseEverywhere ;
    // Then the user's noiseScale could vary that linearly:
    normalization  *= noiseScale;
    double totalSignal = 0.;
    double totalNoise = 0.;
    double newTotalNoiseEverywhere=0.;
    for ( int xIndex = 0; xIndex < unitCellParameters.getNa();  xIndex++){
        for ( int yIndex = 0; yIndex < unitCellParameters.getNb();  yIndex++){
            for ( int zIndex = 0; zIndex < unitCellParameters.getNc();  zIndex++){
                updGridPoint(GridIndices(xIndex,yIndex,zIndex)).noise *= normalization   ;
                updGridPoint(GridIndices(xIndex,yIndex,zIndex)).density = updGridPoint(GridIndices(xIndex,yIndex,zIndex)).noise + updGridPoint(GridIndices(xIndex,yIndex,zIndex)).noiseFreeDensity;
  
                if (updGridPoint(GridIndices(xIndex,yIndex,zIndex)).noiseFreeDensity > 0.0000001) { // We don't bother averaging in regions of zero density. That would give us an unnaturally low SNR.
                    totalSignal += updGridPoint(GridIndices(xIndex,yIndex,zIndex)).noiseFreeDensity ; // Don't bother dividing through by the number of grid points. That will cancel as we just want the ratio.
                    totalNoise +=  updGridPoint(GridIndices(xIndex,yIndex,zIndex)).noise;
                }
                newTotalNoiseEverywhere += updGridPoint(GridIndices(xIndex,yIndex,zIndex)).noise;
    }}}
    double signalToNoiseRatio = totalSignal / totalNoise;
    MMBLOG_FILE_FUNC_LINE(INFO, "total signal = "<<totalSignal<<", total noise in dense regions = "<<totalNoise<<", total noise everywhere = "<<newTotalNoiseEverywhere<<",  signalToNoiseRatio (masked, dense regions only) = "<< signalToNoiseRatio << endl);
    
}


// SCF
void DensityMap::populateNoiseMap(){
    // The ground state wavefunction has wavelength equal to twice the box dimension in real space
    double groundFrequencyX = 1/unitCellParameters.getNa()/unitCellParameters.geta()/2;
    double groundFrequencyY = 1/unitCellParameters.getNb()/unitCellParameters.getb()/2;
    double groundFrequencyZ = 1/unitCellParameters.getNc()/unitCellParameters.getc()/2;
    // Max frequency is the Nyquist frequency
    double maxFrequencyX = 1/unitCellParameters.geta()/2;
    double maxFrequencyY = 1/unitCellParameters.getb()/2;
    double maxFrequencyZ = 1/unitCellParameters.getc()/2;
    double myNoise = .0;
    double averageSignal = 0.;
    double averageNoise = 0.;
    double totalNoiseEverywhere = 0.;
    double signalToNoiseRatio = 0.;
    double inverseNoiseTemperatureToPower4 = 1/pow(noiseTemperature,4); //Stefan–Boltzmann law says total radiance goes like T^4. So we normalize by this number to keep noise sort of constant with temperatuere. Of course our oven has a maximum wavenumber, so this won't be perfect.
    AmplitudeFrequencyAndRandomPhases myAmpFreqRandPhase;
    double myTotalVolume = unitCellParameters.totalVolume(); // compute this only once, cuz there are many operations (perhaps surprisingly)
    MMBLOG_FILE_FUNC_LINE(INFO, "unitCellParameters.totalVolume() = "<<  unitCellParameters.totalVolume() << endl);
    for ( int xIndex = 0; xIndex < unitCellParameters.getNa();  xIndex++){
        
        //MMBLOG_FILE_FUNC_LINE( " xIndex = "<<xIndex <<  std::endl;
        for ( int yIndex = 0; yIndex < unitCellParameters.getNb();  yIndex++){
            for ( int zIndex = 0; zIndex < unitCellParameters.getNc();  zIndex++){
                updGridPoint(GridIndices(xIndex,yIndex,zIndex)).noise=0.;
                for ( int xPhaseIndex = 1; xPhaseIndex <= unitCellParameters.calcMaxFrequencyDoublingsX ();  xPhaseIndex++){
                    for ( int yPhaseIndex = 1; yPhaseIndex <= unitCellParameters.calcMaxFrequencyDoublingsY () ;  yPhaseIndex++){
                        for ( int zPhaseIndex = 1; zPhaseIndex <= unitCellParameters.calcMaxFrequencyDoublingsZ () ;  zPhaseIndex++){
                            myAmpFreqRandPhase = vectorOfAmplitudeFrequencyAndRandomPhases[zPhaseIndex][yPhaseIndex][xPhaseIndex];
                            //MMBLOG_FILE_FUNC_LINE( "  myAmpFreqRandPhase.phaseX = "<< myAmpFreqRandPhase.phaseX<<std::endl;
                            //MMBLOG_FILE_FUNC_LINE( "  myAmpFreqRandPhase.phaseY = "<< myAmpFreqRandPhase.phaseY<<std::endl;
                            //MMBLOG_FILE_FUNC_LINE( "  myAmpFreqRandPhase.phaseZ = "<< myAmpFreqRandPhase.phaseZ<<std::endl;
                            myNoise = //noiseScale  // this is the global noise scale. Moved this down so it is out of the inner loop.      
                                myAmpFreqRandPhase.amplitude                                                                        //  amplitude for this wavenumber vector, taken from planck's law
                                * sin(xIndex * unitCellParameters.geta() * myAmpFreqRandPhase.frequencyX + myAmpFreqRandPhase.phaseX) // position in X, times frequency in X, with random phase in X.
                                * sin(yIndex * unitCellParameters.getb() * myAmpFreqRandPhase.frequencyY + myAmpFreqRandPhase.phaseY)
                                * sin(zIndex * unitCellParameters.getc() * myAmpFreqRandPhase.frequencyZ + myAmpFreqRandPhase.phaseZ) ;
                            updGridPoint(GridIndices(xIndex,yIndex,zIndex)).noise += myNoise; // Save the noise, for later debugging.
                            //noiseMap[zIndex][yIndex][xIndex] = myNoise;
                }}} // of for zPhaseIndex, yPhaseIndex, xPhaseIndex
                //updGridPoint(GridIndices(xIndex,yIndex,zIndex)).noise *= noiseScale;                                             // We scale the noise BEFORE squaring
                updGridPoint(GridIndices(xIndex,yIndex,zIndex)).noise *= updGridPoint(GridIndices(xIndex,yIndex,zIndex)).noise ; // We want intensity, not amplitude. So take the square
                updGridPoint(GridIndices(xIndex,yIndex,zIndex)).noise /= myTotalVolume;                                          // Now we normalize by total volume.                   
                //updGridPoint(GridIndices(xIndex,yIndex,zIndex)).noise *= noiseScale;                                             // Scale AFTER squaring.                                 
                //  we keep noise separate from signal here:
                updGridPoint(GridIndices(xIndex,yIndex,zIndex)).noiseFreeDensity = updGridPoint(GridIndices(xIndex,yIndex,zIndex)).density ;
                if (updGridPoint(GridIndices(xIndex,yIndex,zIndex)).noiseFreeDensity > 0.0000001) { // We don't bother averaging in regions of zero density. That would give us an unnaturally low SNR.
                    averageSignal += updGridPoint(GridIndices(xIndex,yIndex,zIndex)).noiseFreeDensity ; // Don't bother dividing through by the number of grid points. That will cancel as we just want the ratio.
                    averageNoise +=  updGridPoint(GridIndices(xIndex,yIndex,zIndex)).noise;
                }
                totalNoiseEverywhere +=  updGridPoint(GridIndices(xIndex,yIndex,zIndex)).noise; // Unlike the averageNoise which is masked, this is added everywhere.
                // inverseNoiseTemperatureToPower4 is a normalization based on the Stefan–Boltzmann law
                // might want to remove now. then again, may make for better numerical behavior
                // This is now added in normalizeNoiseMap:
                // updGridPoint(GridIndices(xIndex,yIndex,zIndex)).density += inverseNoiseTemperatureToPower4 * updGridPoint(GridIndices(xIndex,yIndex,zIndex)).noise; // Go ahead and add the noise to the density.
    }}} // of for zIndex, yIndex, xIndex
    if (averageNoise>0.) signalToNoiseRatio = averageSignal / averageNoise ; //signalToNoiseRatio / (unitCellParameters.getNa() * unitCellParameters.getNb() * unitCellParameters.getNc()); // Divide through by the total number of map points.
    //MMBLOG_FILE_FUNC_LINE( " total signal = "<<averageSignal<<", total noise in dense regions = "<<averageNoise<<", total noise everywhere = "<<totalNoiseEverywhere<<",  signalToNoiseRatio (masked, dense regions only) = "<< signalToNoiseRatio << std::endl;
    normalizeNoiseMap(totalNoiseEverywhere);
    // moved to Repel.cpp:
    //densityAutocorrelation(1,0); // Arguments are : calculate noise correlation = 1, calculate density corrleation = 0
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    //densityAutocorrelation(0,1);
    MMBLOG_FILE_FUNC_LINE(INFO, endl);

}

/*void DensityMap::addNoiseToDensity(){I
    for ( int xIndex = 0; xIndex < unitCellParameters.getNa(); xIndex++){
        for ( int yIndex = 0; yIndex < unitCellParameters.getNb();  yIndex++){
            for ( int zIndex = 0; zIndex < unitCellParameters.getNc();  zIndex++){
                updGridPoint(GridIndices(xIndex,yIndex,zIndex)).density +=  noiseMap[zIndex][yIndex][xIndex] ;
    }}}
}*/


void DensityMap::initializeVectorOfAmplitudeAndRandomPhases(){
        //validateGridParameters();
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        resizeVectorOfAmplitudeAndRandomPhases();
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        //validateVectorOfAmplitudeFrequencyAndRandomPhasesSize();
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        unitCellParameters.validate();
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        Vec3 tempPosition(0,0,0);	
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        AmplitudeFrequencyAndRandomPhases myAmpFreqPhase ;
        double myFrequency = 0.0;
        srand (time(NULL)); // initialize rand()
        // Just never access index 0 in x,y,z. Makes the math slightly easier.
        MMBLOG_FILE_FUNC_LINE(INFO, "unitCellParameters.geta() = "<<unitCellParameters.geta()<<endl);
        MMBLOG_FILE_FUNC_LINE(INFO, "unitCellParameters.getNa() = "<<unitCellParameters.getNa()<<endl);
        for ( int xIndex = 1; xIndex <= unitCellParameters.calcMaxFrequencyDoublingsX ();xIndex++) {
                for ( int yIndex = 1; yIndex <= unitCellParameters.calcMaxFrequencyDoublingsY ();yIndex++) {
                        for ( int zIndex = 1; zIndex <= unitCellParameters.calcMaxFrequencyDoublingsZ ();zIndex++) {
                            // Recall we take c = 1: So omega = 2 * pi / lambda. (getNa()-1)*geta() is lambda/2. xIndex is the  frequency multiplier :
                            // Recall we take c = 1: So omega = 2 * pi / lambda. (getNa()-1)*geta() is lambda/2. xIndex is the number of frequency doublings. :
                            //                                             <- ranges from 1 to number of grid spacings->     <------------max lambda = cell width * 2----------------->   
                            myAmpFreqPhase.frequencyX = 2*(double)SimTK::Pi* xIndex                                      /(  ((unitCellParameters.getNa()-1)*unitCellParameters.geta()*2)   );  
                            myAmpFreqPhase.frequencyY = 2*(double)SimTK::Pi* yIndex                                      /(  ((unitCellParameters.getNb()-1)*unitCellParameters.getb()*2)   );  
                            myAmpFreqPhase.frequencyZ = 2*(double)SimTK::Pi* zIndex                                      /(  ((unitCellParameters.getNc()-1)*unitCellParameters.getc()*2)   );  
                            //myAmpFreqPhase.frequencyX = 2*(double)SimTK::Pi/((unitCellParameters.getNa()-1)*unitCellParameters.geta()*2*std::pow(2,(-xIndex)));  
                            //myAmpFreqPhase.frequencyY = 2*(double)SimTK::Pi/((unitCellParameters.getNb()-1)*unitCellParameters.getb()*2*std::pow(2,(-yIndex)));  
                            //myAmpFreqPhase.frequencyZ = 2*(double)SimTK::Pi/((unitCellParameters.getNc()-1)*unitCellParameters.getc()*2*std::pow(2,(-zIndex)));  
                            // There is actually a much faster way to take powers of 2: https://stackoverflow.com/questions/39693509/fast-integer-power-of-two
                            //myAmpFreqPhase.frequencyY = 2*(double)SimTK::Pi/(yIndex*unitCellParameters.getb());  
                            //myAmpFreqPhase.frequencyZ = 2*(double)SimTK::Pi/(zIndex*unitCellParameters.getc());  
                            //MMBLOG_FILE_FUNC_LINE( "  myAmpFreqPhase.frequencyX = "<< myAmpFreqPhase.frequencyX <<std::endl;
                            //MMBLOG_FILE_FUNC_LINE( "  myAmpFreqPhase.frequencyY = "<< myAmpFreqPhase.frequencyX <<std::endl;
                            //MMBLOG_FILE_FUNC_LINE( "  myAmpFreqPhase.frequencyZ = "<< myAmpFreqPhase.frequencyX <<std::endl;
                            myFrequency=  (sqrt(myAmpFreqPhase.frequencyX*myAmpFreqPhase.frequencyX + myAmpFreqPhase.frequencyY*myAmpFreqPhase.frequencyY + myAmpFreqPhase.frequencyZ*myAmpFreqPhase.frequencyZ )) ;
                            myAmpFreqPhase.amplitude = plancksLaw(noiseTemperature, myFrequency);
                            //MMBLOG_FILE_FUNC_LINE( "wavenumber doublings x,y,z, noiseTemperature, myFrequency,freqX,freqY,freqZ, plancksLaw(noiseTemperature, myFrequency) =              "<<","<< xIndex <<","<< yIndex <<","<<zIndex <<","  <<noiseTemperature<<","<<myFrequency<<",>"<<  myAmpFreqPhase.frequencyX<<","<<   myAmpFreqPhase.frequencyY  <<","<< myAmpFreqPhase.frequencyZ <<"<,"<<    myAmpFreqPhase.amplitude<<std::endl  ;
                            // random phases in radians:
                            myAmpFreqPhase.phaseX = rand() / (double)RAND_MAX * SimTK::Pi * 2.0; // Have to cast (double)RAND_MAX to prevent the result of the division from being cast as int.
                            myAmpFreqPhase.phaseY = rand() / (double)RAND_MAX * SimTK::Pi * 2. ;
                            myAmpFreqPhase.phaseZ = rand() / (double)RAND_MAX * SimTK::Pi * 2. ;
                            //myAmpFreqPhase.phaseX = 0.;
                            //myAmpFreqPhase.phaseY = 0.;
                            //myAmpFreqPhase.phaseZ = 0.;
                            //MMBLOG_FILE_FUNC_LINE(" setting vectorOfAmplitudeFrequencyAndRandomPhases"<<zIndex<<","<<   yIndex <<","<<  xIndex<<std::endl;
                            vectorOfAmplitudeFrequencyAndRandomPhases[zIndex][yIndex][xIndex]=myAmpFreqPhase;  
                            //MMBLOG_FILE_FUNC_LINE(std::endl;
        }}}
        //validateVectorOfAmplitudeFrequencyAndRandomPhasesSize()
}


void DensityMap::loadParametersAndDensity(const String &densityFileName)
{
    unsigned int extIndex = densityFileName.rfind(".");
    String extension = densityFileName.substr(extIndex);
    if(extension == ".xplor")
        loadParametersAndDensity_XPLOR(densityFileName);
    else if(extension == ".dx")
        { //loadParametersAndDensity_OpenDX(densityFileName);
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unable to open density map : "<<densityFileName<<" .. DX is not a supported format at the moment."<<endl);
        }
    else if(extension == ".situs" || extension == ".sit")
        {//loadParametersAndDensity_Situs(densityFileName);
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unable to open density map : "<<densityFileName<<" .. Situs is nota  supported format at the moment"<<endl);
        }
    else if (extension == ".map"    || extension == ".ccp4"    || extension == ".mrc"    ||
             extension == ".map.gz" || extension == ".ccp4.gz" || extension == ".mrc.gz" )
    {
        loadParametersAndDensity_CCP4MAP              ( densityFileName );
    }
    else
    {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "DensityMap: Extension unknown for " << densityFileName << endl);
    }
}

/*! \brief Function responsible for reading in the CCP4's MAP format.

    This function uses the Gemmi library to read in the CCP4 formatted map data. It makes use of the setup function from
    Gemmi, so that the symmetry operations and other unpleasantries should be already deal with when the map is read.
 
    \warning This function requires the Gemmi library and will not work unless this is installed.
    \warning This function does only limitted reporting to the log, please add appropriate log outputs for the procedure.

    \param[in] densityFileName String containing the path to the density file which should read in.
*/
void DensityMap::loadParametersAndDensity_CCP4MAP(const String &densityFileName)
{
    //================================================ Open the file using Gammi
    gemmi::Ccp4<float> map;
    map.read_ccp4                                     ( gemmi::MaybeGzipped ( densityFileName.c_str() ) );
    
    //================================================ Read in the axes starting points before it is modified by the gemmi set-up
    int xFrom                                         = static_cast<int> ( map.header_i32   ( 5  ) );
    int yFrom                                         = static_cast<int> ( map.header_i32   ( 6  ) );
    int zFrom                                         = static_cast<int> ( map.header_i32   ( 7  ) );
    
    //================================================ Convert to XYZ and create complete map, if need be
    map.setup                                         ( gemmi::GridSetup::ReorderOnly, NAN );
    
    //================================================ Parse the map header for cell info
    int xDimInds                                      = static_cast<int> ( map.header_i32   ( 1  ) );
    int yDimInds                                      = static_cast<int> ( map.header_i32   ( 2  ) );
    int zDimInds                                      = static_cast<int> ( map.header_i32   ( 3  ) );
    
    float xDim                                        = static_cast<float> ( map.header_float ( 11 ) );
    float yDim                                        = static_cast<float> ( map.header_float ( 12 ) );
    float zDim                                        = static_cast<float> ( map.header_float ( 13 ) );
    
    float aAng                                        = static_cast<float> ( map.header_float ( 14 ) );
    float bAng                                        = static_cast<float> ( map.header_float ( 15 ) );
    float cAng                                        = static_cast<float> ( map.header_float ( 16 ) );
    
    int xAxOrigin                                     = static_cast<int> ( map.header_i32   ( 50 ) ) + xFrom;
    int yAxOrigin                                     = static_cast<int> ( map.header_i32   ( 51 ) ) + yFrom;
    int zAxOrigin                                     = static_cast<int> ( map.header_i32   ( 52 ) ) + zFrom;
    
    int xAxOrder                                      = static_cast<int> ( map.header_i32   ( 17 ) );
    int yAxOrder                                      = static_cast<int> ( map.header_i32   ( 18 ) );
    int zAxOrder                                      = static_cast<int> ( map.header_i32   ( 19 ) );
    
    int xGridInds                                     = xDimInds;
    int yGridInds                                     = yDimInds;
    int zGridInds                                     = zDimInds;
    
    //================================================ Set the N for unitCellParameters object
    unitCellParameters.setN                           ( xDimInds,                   // Number of indices along the x-axis
                                                        xFrom,                      // X-axis from index
                                                        xFrom + xDimInds - 1,       // X-axis to index
                                                        yDimInds,                   // Number of indices along the y-axis
                                                        yFrom,                      // Y-axis from index
                                                        yFrom + yDimInds - 1,       // Y-axis to index
                                                        zDimInds,                   // Number of indices along the z-axis
                                                        zFrom,                      // Z-axis from index
                                                        zFrom + zDimInds - 1 );     // Z-axis to index
    
    //================================================ Set dimensions and angles for unitCellParameters object
    unitCellParameters.setabc                         ( ( xDim / ( xDimInds - 1 ) ) / 10.0 ,  // Distance between two indices along the x-axis in nm ( divide by 10 to get nm from A )
                                                        ( yDim / ( yDimInds - 1 ) ) / 10.0 ,  // Distance between two indices along the y-axis in nm ( divide by 10 to get nm from A )
                                                        ( zDim / ( zDimInds - 1 ) ) / 10.0 ); // Distance between two indices along the z-axis in nm ( divide by 10 to get nm from A )
    unitCellParameters.setAlphaUsingDegrees           ( aAng );
    unitCellParameters.setBetaUsingDegrees            ( bAng );
    unitCellParameters.setGammaUsingDegrees           ( cAng );
    
    //================================================ De-orthoginalisation matrix
    unitCellParameters.setDeOrthogonalizationMatrix   ( );
    
    //================================================ Prepate grid
    initializeArrayOfGridPoints                       ( );
    
    //================================================ Read in the map
    {
        //============================================ Initialise internal variables
        int axOrdArr[3];
        int axDimArr[3];
        int newU, newV, newW;

        //============================================ Fill in values
        axDimArr[0]                                   = xDimInds;
        axDimArr[1]                                   = yDimInds;
        axDimArr[2]                                   = zDimInds;

        //============================================ Copy read in data to internal map variable
        for ( axOrdArr[0] = 0; axOrdArr[0] < axDimArr[xAxOrder-1]; axOrdArr[0]++ )
        {
            for ( axOrdArr[1] = 0; axOrdArr[1] < axDimArr[yAxOrder-1]; axOrdArr[1]++ )
            {
                for ( axOrdArr[2] = 0; axOrdArr[2] < axDimArr[zAxOrder-1]; axOrdArr[2]++ )
                {
                    GridIndices centralGridIndices    ( axOrdArr[0],  axOrdArr[1], axOrdArr[2] );
                    GridPoint & centralGridPoint      = updGridPoint ( centralGridIndices );
                    setDensity                        ( centralGridPoint, static_cast < float > ( map.grid.get_value_q( axOrdArr[xAxOrder-1], axOrdArr[yAxOrder-1], axOrdArr[zAxOrder-1] ) ) );
                }
            }
        }
    }

    //======================================== DONE!
    return ;
}

// Expects .xplor density maps to be in the following format:
// http://psb11.snv.jussieu.fr/doc-logiciels/msi/xplor981/formats.html -- not accessible anymore
// http://chem5.nchc.org.tw/software/document2007/insightII/doc/xplor/formats.html
void DensityMap::loadParametersAndDensity_XPLOR(const String &densityFileName) {

        ifstream inFile(densityFileName.c_str(),ifstream::in);
        int densitiesPerLine = 6;	
        	//const int numFields = 10 ;
        vector<String> mystring;
        //String mystring2[numFields];
        	stringstream u;
        	String tempString;
        if (! (inFile.is_open())) {

            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unable to open density map : "<<densityFileName<<endl);
        }
        readAndParseLine  (inFile); // line 1. this is expected to be a blank line
        mystring =readAndParseLine  (inFile); // line 2. This should have number of following remark/title lines. then the title itself. We will need the number of blank lines.
        int skipLines = atoi (mystring[0].c_str());
        for (int i = 0; i < skipLines; i++) {
            readAndParseLine  (inFile);}
        //readAndParseLine  (inFile);
        mystring = readAndParseLine  (inFile); // Now time to read the number of grid points
        unitCellParameters.setN(
            atoi (mystring[0].c_str()),
            atoi (mystring[1].c_str()),
            atoi (mystring[2].c_str()),
            atoi (mystring[3].c_str()),
            atoi (mystring[4].c_str()),
            atoi (mystring[5].c_str()),
            atoi (mystring[6].c_str()),
            atoi (mystring[7].c_str()),
            atoi (mystring[8].c_str())
            );

        //initializeArrayOfGridPoints();
        mystring = readAndParseLine  (inFile);
        unitCellParameters.setabc(
            atof (mystring[0].c_str())/(unitCellParameters.getNa()-1) /10. , //dividing by 10 to convert from Å(the XPLOR format, per Alwyn Jones) to nm (molmodel units). the "-1" is because we want the space between grid points, not the grid points themselves.
            atof (mystring[1].c_str())/(unitCellParameters.getNb()-1) /10. ,
            atof (mystring[2].c_str())/(unitCellParameters.getNc()-1) /10. );

        unitCellParameters.setAlphaUsingDegrees(atof(mystring[3].c_str()));
        unitCellParameters.setBetaUsingDegrees (atof(mystring[4].c_str()));
        unitCellParameters.setGammaUsingDegrees(atof(mystring[5].c_str()));
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        unitCellParameters.setDeOrthogonalizationMatrix ();
        //unitCellParameters.setOrthogonalizationMatrix ();
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        initializeArrayOfGridPoints();
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        mystring = readAndParseLine  (inFile);
        if ((mystring[0]).compare("ZYX") != 0    ) {MMBLOG_FILE_FUNC_LINE(CRITICAL, "expected ZYX, got :"<<mystring[0]<<endl);}


        for ( int zIndex = 0; zIndex < unitCellParameters.getNc(); zIndex ++) {
        	// read z-index
                //MMBLOG_FILE_FUNC_LINE(" unitCellParameters.getNc() " <<unitCellParameters.getNc()<<std::endl;
                mystring = readAndParseLine  (inFile);
        	if (atoi ((mystring[0].c_str())) != zIndex) {
        			MMBLOG_FILE_FUNC_LINE(CRITICAL, "Something is wrong with the input file reading.  expected to read zIndex "<<zIndex<< " and instead read : "<<mystring[0]<<endl);
        	}	
        	

        	for ( int yIndex = 0; yIndex < unitCellParameters.getNb(); yIndex ++) 
        	for ( int xIndex = 0; xIndex < unitCellParameters.getNa(); xIndex = xIndex + 0) 

        	{
                        //MMBLOG_FILE_FUNC_LINE(std::endl;
        		mystring = readAndParseOnColWidth (inFile,12);
        		if ((int)mystring.size() > densitiesPerLine) {
        			MMBLOG_FILE_FUNC_LINE(CRITICAL, "Too many densities on this line!  expected "<<densitiesPerLine<< "and found "<<mystring.size()<<endl);
        		}
        		//xIndex --;
        		for (int i = 0; i<(int)mystring.size() ;i++) {
        			if (xIndex == unitCellParameters.getNa()){xIndex = 0; yIndex++; }
        			if ((xIndex+1 + (yIndex*unitCellParameters.getNa()) ) <= (unitCellParameters.getNa()* unitCellParameters.getNb()))
        			{
        				GridIndices centralGridIndices((xIndex ),  yIndex, zIndex);
        				GridPoint & centralGridPoint = updGridPoint(centralGridIndices);
                                        //MMBLOG_FILE_FUNC_LINE(" For point xIndex, yIndex, zIndex : "<<  xIndex<<", "<< yIndex<<", "<< zIndex<<" setting density to "<<mystring[i].c_str()<<std::endl;

        				setDensity(centralGridPoint, atof(mystring[i].c_str()));	
        			} else{ 
        				MMBLOG_FILE_FUNC_LINE(CRITICAL, " Something went wrong .. too many densities !  "<<endl);
        		        }
        			xIndex++;
        		}
        		
        	}
        } // of zIndex
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        inFile.close();
}

void incrementer(int & myIndex, const int maxIndex) {
    if (myIndex < maxIndex) myIndex ++;
    if (myIndex == maxIndex) {myIndex = 0;}
    if (myIndex > maxIndex) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unexplained error! index of  "<<myIndex<<" exceeds maximum of "<<maxIndex<<endl);
    }
}

void incrementer(int & xIndex, const int maxX, int & yIndex, const int maxY, int & zIndex, const int maxZ, int & counter, const int maxCounter){
    //MMBLOG_FILE_FUNC_LINE(" xIndex, maxX : "<<xIndex<<", "<<maxX <<" yIndex, maxY : "<<yIndex<<", "<<maxY <<" zIndex, maxZ : "<<zIndex<<", "<<maxZ<<std::endl;
    incrementer(xIndex, maxX);
    if (xIndex == 0) {incrementer(yIndex,maxY);} // if x rolls over, then increment y
    if ((yIndex == 0) && (xIndex == 0)) {
        incrementer(zIndex, maxZ);
    } // if x and y both roll over, increment z
    incrementer(counter, maxCounter); // Always increment counter
}

void DensityMap::writeDensityMapXplor(const String &densityFileName, const bool writeDensity, const bool writeNoise ) {
        ofstream outFile(densityFileName.c_str(),ofstream::out);
        int densitiesPerLine = 6;	
        vector<String> mystring;
       	stringstream u;
       	String tempString;
        if (! (outFile.is_open())) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unable to open density map : "<<densityFileName<<endl);
        }
        outFile <<std::endl; // write blank line out.
        outFile << "0 !NTITLE"<< std::endl; // Number of lines to skip
        
        outFile <<unitCellParameters.getNa()<<" "<< unitCellParameters. getaMin()<<" "<< unitCellParameters.getaMax()<<" " <<unitCellParameters.getNb()<<" "<< unitCellParameters.getbMin()<<" "<< unitCellParameters.getbMax()<<" " <<unitCellParameters.getNc()<<" "<< unitCellParameters.getcMin()<<" "<< unitCellParameters.getcMax()<<std::endl;
        outFile << 10*unitCellParameters.geta()*(unitCellParameters.getNa() - 1) <<" " << 10*unitCellParameters.getb()*(unitCellParameters.getNb() - 1)<<" " <<10*unitCellParameters.getc()*(unitCellParameters.getNc() - 1) << " " << unitCellParameters.getAlpha() / SimTK::Pi * 180.0  << " " << unitCellParameters.getBeta() / SimTK::Pi * 180.0  << " " << unitCellParameters.getGamma() / SimTK::Pi * 180.0  <<std::endl;
        
        MMBLOG_FILE_FUNC_LINE(INFO, std::endl);
        outFile << "ZYX"; //<<std::endl;

        int xIndex = 0; 
        int yIndex = 0; 
        int zIndex = 0;
        int counter = 0;
                 bool keepGoing = 1;
                 while(keepGoing)
                 {
                     //MMBLOG_FILE_FUNC_LINE(" xIndex = "<<xIndex<<" yIndex = "<<yIndex<<" zIndex = "<<zIndex<<std::endl;
                     if ((yIndex == 0) && (xIndex == 0)) {
                         if (counter !=0) {outFile <<std::endl; // Carriage return. To terminate the incomplete density line.     .
                             counter = 0;
                         }
                         MMBLOG_FILE_FUNC_LINE(DEBUG, "xIndex = "<<xIndex<<" yIndex = "<<yIndex<<" zIndex = "<<zIndex<<endl);
                         outFile <<std::endl; // Carriage return.
                         outFile <<zIndex<<std::endl;
                     } else if (counter == 0 ){
                         //MMBLOG_FILE_FUNC_LINE(" xIndex = "<<xIndex<<" yIndex = "<<yIndex<<" zIndex = "<<zIndex<<std::endl;
                         outFile <<std::endl; // Carriage return. but only if we did not just print zIndex above.
                     }
                     //MMBLOG_FILE_FUNC_LINE(" xIndex = "<<xIndex<<" yIndex = "<<yIndex<<" zIndex = "<<zIndex<<std::endl;
                     // Write a density:
                     GridIndices centralGridIndices((xIndex ),  yIndex, zIndex);
                     GridPoint & centralGridPoint = updGridPoint(centralGridIndices);
                     outFile <<" "<<setw(11)<< (centralGridPoint.noiseFreeDensity * writeDensity) + (centralGridPoint.noise * writeNoise);	// One density for each 12-characgter column. Depending on parameters passed, one can write the original density, the noise, or the sum of the two.
                     // now increment indices:
                     incrementer(xIndex, unitCellParameters.getNa(), yIndex, unitCellParameters.getNb(), zIndex , unitCellParameters.getNc(), counter, 6);
                      

                     if ((zIndex == 0) && (yIndex == 0) && (xIndex == 0)) {
                         //MMBLOG_FILE_FUNC_LINE(" xIndex = "<<xIndex<<" yIndex = "<<yIndex<<" zIndex = "<<zIndex<<std::endl;
                         keepGoing = 0;  // All numbers have rolled over.  exit while
                         //zIndex =  unitCellParameters.getNc(); // exit while
                     }
                 } // of while(keepGoing)
                 outFile <<std::endl; // Carriage return.
                 outFile <<-9999<<std::endl;
                 outFile <<"1 1"<<std::endl; // This is supposed to be average and std dev of density. But right now just left at a constant value.
        //  // of zIndex
        outFile.close();
}

// Expects Situs density maps
// http://situs.biomachina.org/fmap.pdf
/*
void DensityMap::loadParametersAndDensity_Situs(const String densityFileName) {

        ifstream inFile(densityFileName.c_str(),ifstream::in);
        int densitiesPerLine = 10;   
            //const int numFields = 10 ;
        vector<String> mystring;
        //String mystring2[numFields];
            stringstream u;
            String tempString;
        if (! (inFile.is_open())) {

            MMBLOG_FILE_FUNC_LINE(CRITICAL, " Unable to open density map : "<<densityFileName<<endl;
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
        //unitCellParameters.getNc() = unitCellNumGridZ;
        cout<<__FILE__<<":"<<__LINE__<<" set unitCellParameters.getNc(),Y,X: "<<totalNumGridX<<","<< totalNumGridY<<","<< unitCellParameters.getNc()<<endl;

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
        maxZ = minZ + ((unitCellParameters.getNc()-1) * gridZSpacing);
        initializeArrayOfGridPoints();

        cout<<__FILE__<<":"<<__LINE__<<" minimum extents in  Z,Y,X : "<<minZ*10<<","<< minY*10<<","<< minX*10<<endl;
        cout<<__FILE__<<":"<<__LINE__<<" maximum extents in  Z,Y,X : "<<maxZ*10<<","<< maxY*10<<","<< maxX*10<<endl;

        readAndParseLine  (inFile);



        for ( int zIndex = 0; zIndex < unitCellParameters.getNc(); zIndex ++) 
        {
            for ( int yIndex = 0; yIndex < totalNumGridY; yIndex ++) 
            {   
                for ( int xIndex = 0; xIndex < totalNumGridX; xIndex = xIndex + 0) 
                {
                    mystring = readAndParseLine(inFile);
                    if ((int)mystring.size()-1 > densitiesPerLine) {
                        cout << __FILE__<<":"<<__LINE__<<": " << mystring[0]<< " " << mystring[10] << endl;
                        MMBLOG_FILE_FUNC_LINE(CRITICAL, " Too many densities on this line!  expected "<<densitiesPerLine<< " and found "<<mystring.size()<<endl;
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
                            MMBLOG_FILE_FUNC_LINE(CRITICAL, " Something went wrong .. too many densities !  "<<endl;
                            }
                        xIndex++;
                    }
                
                }
            }
        } // of zIndex
        inFile.close();
} */

void DensityMap::precomputeGradient() {
			for ( int xIndex = 0; xIndex < unitCellParameters.getNa(); xIndex ++) 
			for ( int yIndex = 0; yIndex < unitCellParameters.getNb(); yIndex ++) 
			for ( int zIndex = 0; zIndex < unitCellParameters.getNc(); zIndex ++) 

			{
				GridIndices centralGridIndices(xIndex,  yIndex, zIndex);
				GridPoint & centralGridPoint = updGridPoint(centralGridIndices);
				initializeGradient(centralGridPoint); // sets gradient to zero; if the if statements below are false that component of the gradient remains zero.

				if (xIndex < (unitCellParameters.getNa()-1)) {// if grid point is not at the +X boundary
                                        MMBLOG_FILE_FUNC_LINE(DEBUG, endl);
					setPositiveXGradient (centralGridPoint,(getDensity(updGridPoint(GridIndices(xIndex+1,  yIndex, zIndex))) -getDensity(centralGridPoint))/ unitCellParameters.geta());
				} else if (xIndex == (unitCellParameters.getNa()-1)) { // if grid point is at the +X boundary
                                        MMBLOG_FILE_FUNC_LINE(DEBUG, endl);
					setPositiveXGradient (centralGridPoint,(0. -getDensity(centralGridPoint)) / unitCellParameters.geta() );
                                }

				if (yIndex < (unitCellParameters.getNb()-1)) {// if grid point is not at the +Y boundary
                                        MMBLOG_FILE_FUNC_LINE(DEBUG, endl);
					setPositiveYGradient (centralGridPoint,(getDensity(updGridPoint(GridIndices(xIndex,  yIndex+1, zIndex))) -getDensity(centralGridPoint))/ unitCellParameters.getb());
				} else if (yIndex == (unitCellParameters.getNb()-1)) {
                                        MMBLOG_FILE_FUNC_LINE(DEBUG, endl);
					setPositiveYGradient (centralGridPoint,(0. -getDensity(centralGridPoint)) / unitCellParameters.getb() );
                                }



				if (zIndex < (unitCellParameters.getNc()-1)) {// if grid point is not at the +Z boundary:
                                        MMBLOG_FILE_FUNC_LINE(DEBUG, endl);
					setPositiveZGradient (centralGridPoint,(getDensity(updGridPoint(GridIndices(xIndex,  yIndex, zIndex+1))) -getDensity(centralGridPoint))/ unitCellParameters.getc());
				
				} else if (zIndex == (unitCellParameters.getNc()-1)) {
                                        MMBLOG_FILE_FUNC_LINE(DEBUG, endl);
					setPositiveZGradient (centralGridPoint,(0. -getDensity(centralGridPoint)) / unitCellParameters.getc() );
                                }
				//ValidateVec3((fetchFirstQuadrantGradient(centralGridPoint)));



			}
		}

void DensityMap::precomputeGradientDerivatives() {
			for ( int xIndex = 0; xIndex < unitCellParameters.getNa(); xIndex ++) 
			for ( int yIndex = 0; yIndex < unitCellParameters.getNb(); yIndex ++) 
			for ( int zIndex = 0; zIndex < unitCellParameters.getNc(); zIndex ++) 

			{
				GridIndices centralGridIndices(xIndex,  yIndex, zIndex);
				GridPoint & centralGridPoint = updGridPoint(centralGridIndices);

				if (xIndex < (unitCellParameters.getNa()-1)) {// if grid point is not at the +X boundary
					setddxPositiveXGradient (centralGridPoint,(fetchFirstQuadrantGradient(updGridPoint(GridIndices(xIndex+1,  yIndex, zIndex)))[0] -fetchFirstQuadrantGradient(centralGridPoint)[0]) );//  / unitCellParameters.geta());
					setddxPositiveYGradient (centralGridPoint,(fetchFirstQuadrantGradient(updGridPoint(GridIndices(xIndex+1,  yIndex, zIndex)))[1] -fetchFirstQuadrantGradient(centralGridPoint)[1]));// / unitCellParameters.geta());
					setddxPositiveZGradient (centralGridPoint,(fetchFirstQuadrantGradient(updGridPoint(GridIndices(xIndex+1,  yIndex, zIndex)))[2] -fetchFirstQuadrantGradient(centralGridPoint)[2]) );//  / unitCellParameters.geta());
				} else if (xIndex == (unitCellParameters.getNa()-1)) { // if grid point is at the + X boundary
					setddxPositiveXGradient (centralGridPoint,(0 - fetchFirstQuadrantGradient(centralGridPoint)[0]) );// / unitCellParameters.geta());
					setddxPositiveYGradient (centralGridPoint,(0 - fetchFirstQuadrantGradient(centralGridPoint)[1]) );// / unitCellParameters.geta());
					setddxPositiveZGradient (centralGridPoint,(0 -fetchFirstQuadrantGradient(centralGridPoint)[2])  );// / unitCellParameters.geta());
                                }


				if (yIndex < (unitCellParameters.getNb()-1)) {// if grid point is not at the +Y boundary
                                        // Actually since Hessians are usually symmetric, we should just set this equal to ddxPositiveYGradient .
					setddyPositiveXGradient (centralGridPoint,(fetchFirstQuadrantGradient(updGridPoint(GridIndices(xIndex,  yIndex+1, zIndex)))[0] -fetchFirstQuadrantGradient(centralGridPoint)[0]));/// unitCellParameters.getb());
					setddyPositiveYGradient (centralGridPoint,(fetchFirstQuadrantGradient(updGridPoint(GridIndices(xIndex,  yIndex+1, zIndex)))[1] -fetchFirstQuadrantGradient(centralGridPoint)[1]));/// unitCellParameters.getb());
					setddyPositiveZGradient (centralGridPoint,(fetchFirstQuadrantGradient(updGridPoint(GridIndices(xIndex,  yIndex+1, zIndex)))[2] -fetchFirstQuadrantGradient(centralGridPoint)[2]));/// unitCellParameters.getb());
				
				} else if (yIndex == (unitCellParameters.getNb()-1)) {
					setddyPositiveXGradient (centralGridPoint,(0 - fetchFirstQuadrantGradient(centralGridPoint)[0]));/// unitCellParameters.getb());
					setddyPositiveYGradient (centralGridPoint,(0 -fetchFirstQuadrantGradient(centralGridPoint)[1]));/// unitCellParameters.getb());
					setddyPositiveZGradient (centralGridPoint,(0 -fetchFirstQuadrantGradient(centralGridPoint)[2]));/// unitCellParameters.getb());
                                }

				if (zIndex < (unitCellParameters.getNc()-1)) {// if grid point is not at the +Z boundary:
					setddzPositiveXGradient ( centralGridPoint,(fetchFirstQuadrantGradient(updGridPoint(GridIndices(xIndex,  yIndex, zIndex+1)))[0] -fetchFirstQuadrantGradient(centralGridPoint)[0]));/// unitCellParameters.getc());
					setddzPositiveYGradient ( centralGridPoint, ( fetchFirstQuadrantGradient(updGridPoint(GridIndices(xIndex,  yIndex, zIndex+1)))[1] -fetchFirstQuadrantGradient(centralGridPoint)[1]));/// unitCellParameters.getc());
					setddzPositiveZGradient ( centralGridPoint, ( fetchFirstQuadrantGradient(updGridPoint(GridIndices(xIndex,  yIndex, zIndex+1)))[2] -fetchFirstQuadrantGradient(centralGridPoint)[2]));/// unitCellParameters.getc());
				
				} else if (zIndex == (unitCellParameters.getNc()-1)) {
					setddzPositiveXGradient (centralGridPoint, (0 - fetchFirstQuadrantGradient(centralGridPoint)[0]));/// unitCellParameters.getc());
					setddzPositiveYGradient (centralGridPoint, (0 -fetchFirstQuadrantGradient(centralGridPoint)[1]));/// unitCellParameters.getc());
					setddzPositiveZGradient (centralGridPoint, (0 -fetchFirstQuadrantGradient(centralGridPoint)[2]));/// unitCellParameters.getc());
                                }

			}
		}


/*
Vec3 DensityMap::fetchGradient(Vec3 position)  {
			GridIndices myNearestGridIndices = calcNearestGridIndices(   position);
                        if (hasGridPoint(myNearestGridIndices)) {
				GridPoint myGridPoint =  updGridPoint(myNearestGridIndices);
				return myGridPoint.fetchGradient (position) ;	
			} else {
				return Vec3(0);
			}
		}
		



Vec3 DensityMap::fetchGradient(Vec3 position)  {
                     	GridIndices myNearestGridIndices = calcNearestGridIndices(   position);
                        //GridPoint myGridPoint; 
     			//myGridPoint.initialize();
                         if (hasGridPoint(myNearestGridIndices)) {
                                 //myGridPoint =  updGridPoint(myNearestGridIndices);
                                 return fetchGradient(updGridPoint(myNearestGridIndices), position) ;
                         } else if (hasGridPoint(GridIndices (myNearestGridIndices.getXGridIndex()+1, myNearestGridIndices.getYGridIndex(),myNearestGridIndices.getZGridIndex()))) {
                                GridPoint myGridPoint; 
     			        initialize(myGridPoint);
				setPositiveXGradient(myGridPoint,(getDensity(updGridPoint(GridIndices(myNearestGridIndices.getXGridIndex() +1,  myNearestGridIndices.getYGridIndex() , myNearestGridIndices.getZGridIndex() ))) - 0.) / unitCellParameters.geta()) ;	
         			
          			MMBLOG_FILE_FUNC_LINE(CRITICAL, " this function is obsolete, fetchGradient doesn't work anymore"<<endl;
                         } else if (hasGridPoint(GridIndices (myNearestGridIndices.getXGridIndex()-1, myNearestGridIndices.getYGridIndex(),myNearestGridIndices.getZGridIndex()))) {
                                GridPoint myGridPoint; 
     			        initialize(myGridPoint);
				setNegativeXGradient(myGridPoint,(0. - getDensity(updGridPoint(GridIndices(myNearestGridIndices.getXGridIndex() -1,  myNearestGridIndices.getYGridIndex() , myNearestGridIndices.getZGridIndex() )))) / gridXSpacing) ;	 
          			MMBLOG_FILE_FUNC_LINE(CRITICAL, " this function is obsolete, fetchGradient doesn't work anymore"<<endl;
                         } else if (hasGridPoint(GridIndices (myNearestGridIndices.getXGridIndex(), myNearestGridIndices.getYGridIndex()+1,myNearestGridIndices.getZGridIndex()))) {
                                GridPoint myGridPoint; 
     			        initialize(myGridPoint);
				setPositiveYGradient(myGridPoint,(getDensity(updGridPoint(GridIndices(myNearestGridIndices.getXGridIndex() ,  myNearestGridIndices.getYGridIndex() +1, myNearestGridIndices.getZGridIndex() ))) - 0.) / unitCellParameters.getb()) ;	
                                return fetchGradient (myGridPoint,position) ;
                         } else if (hasGridPoint(GridIndices (myNearestGridIndices.getXGridIndex(), myNearestGridIndices.getYGridIndex()-1,myNearestGridIndices.getZGridIndex()))) {
                                GridPoint myGridPoint; 
     			        initialize(myGridPoint);
				setNegativeYGradient(myGridPoint ,(0. - getDensity(updGridPoint(GridIndices(myNearestGridIndices.getXGridIndex() -1,  myNearestGridIndices.getYGridIndex() -1 , myNearestGridIndices.getZGridIndex() )))) / unitCellParameters.getb()) ;
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
                                 return Vec3(0);
                         }
                 }

*/

Vec3 DensityMap::fetchFirstQuadrantGradient(const Vec3 &position)  {

                        GridIndices myLowerLeftGridIndex = calcLowerLeftGridIndices(   position);
                         if (hasGridPoint(myLowerLeftGridIndex)) {
                                 return fetchFirstQuadrantGradient(ArrayOfGridPoints[myLowerLeftGridIndex.getZGridIndex()][myLowerLeftGridIndex.getYGridIndex()][myLowerLeftGridIndex.getXGridIndex()]);

                         } else if (hasGridPoint(GridIndices (myLowerLeftGridIndex.getXGridIndex()+1, myLowerLeftGridIndex.getYGridIndex(),myLowerLeftGridIndex.getZGridIndex()))) {
                                GridPoint myGridPoint;
                                initialize(myGridPoint);
                                MMBLOG_FILE_FUNC_LINE(DEBUG, endl);
                                setPositiveXGradient(myGridPoint,(getDensity(updGridPoint(GridIndices(myLowerLeftGridIndex.getXGridIndex() +1,  myLowerLeftGridIndex.getYGridIndex() , myLowerLeftGridIndex.getZGridIndex() ))) - 0.) / unitCellParameters.geta()) ;
                                return fetchFirstQuadrantGradient(myGridPoint) ;

                         } else if (hasGridPoint(GridIndices (myLowerLeftGridIndex.getXGridIndex(), myLowerLeftGridIndex.getYGridIndex()+1,myLowerLeftGridIndex.getZGridIndex()))) {
                                GridPoint myGridPoint;
                                initialize(myGridPoint);  
                                MMBLOG_FILE_FUNC_LINE(DEBUG, endl);
                                setPositiveYGradient(myGridPoint,(getDensity(updGridPoint(GridIndices(myLowerLeftGridIndex.getXGridIndex() ,  myLowerLeftGridIndex.getYGridIndex() +1, myLowerLeftGridIndex.getZGridIndex() ))) - 0.) / unitCellParameters.getb()) ;
                                return fetchFirstQuadrantGradient(myGridPoint) ;

                         } else if (hasGridPoint(GridIndices (myLowerLeftGridIndex.getXGridIndex(), myLowerLeftGridIndex.getYGridIndex(),myLowerLeftGridIndex.getZGridIndex()+1))) {
                                GridPoint myGridPoint;
                                initialize(myGridPoint);
                                MMBLOG_FILE_FUNC_LINE(DEBUG, endl);
                                setPositiveZGradient(myGridPoint,(getDensity(updGridPoint(GridIndices(myLowerLeftGridIndex.getXGridIndex() ,  myLowerLeftGridIndex.getYGridIndex() , myLowerLeftGridIndex.getZGridIndex() +1 ))) - 0.) / unitCellParameters.geta()) ;
                                // return myGridPoint.fetchGradient (position) ;
                                return fetchFirstQuadrantGradient(myGridPoint) ;

                         } else {
                                 Vec3 tempVec3(0);
                                 return tempVec3;
                         }

}

Vec3 DensityMap::calcInterpolatedFirstQuadrantGradient(const Vec3 &position)  {

                        GridIndices myLowerLeftGridIndex = calcLowerLeftGridIndices(   position);
                         if (hasGridPoint(myLowerLeftGridIndex)) {
                                 MMBLOG_FILE_FUNC_LINE(DEBUG, endl);
                                 Vec3 tempVec3 = calcInterpolatedFirstQuadrantGradient(ArrayOfGridPoints[myLowerLeftGridIndex.getZGridIndex()][myLowerLeftGridIndex.getYGridIndex()][myLowerLeftGridIndex.getXGridIndex()],position);
                                 //cout<<__FILE__<<":"<<__LINE__<<":"<<__FUNCTION__<<" Returning NON-ZERO force "<< tempVec3 <<" for grid point at position "<<position<<", lower left indices "<<myLowerLeftGridIndex.getXGridIndex() <<", "<< myLowerLeftGridIndex.getYGridIndex()   <<", "<< myLowerLeftGridIndex.getZGridIndex()  <<  endl;
                                 return tempVec3; //calcInterpolatedFirstQuadrantGradient(ArrayOfGridPoints[myLowerLeftGridIndex.getZGridIndex()][myLowerLeftGridIndex.getYGridIndex()][myLowerLeftGridIndex.getXGridIndex()],position);

                         } 
                         else { // might want to trap the conditions at the boundaries of the map, to get the minimizer to work
                                 Vec3 tempVec3(0);
                                 //cout<<__FILE__<<":"<<__LINE__<<":"<<__FUNCTION__<<" Returning ZERO force "<< tempVec3 <<" for grid point at position "<<position<<endl;
                                 return tempVec3;
                         }
}

//SimTK::Vec3 DensityMap::calcInterpolatedFirstQuadrantGradient(SimTK::Vec3 position)  {
//    Vec3 myVec3Float = calcInterpolatedFirstQuadrantGradient(Vec3(  position[0],position[1],position[2]));
//    return Vec3(myVec3Float[0], myVec3Float[1], myVec3Float[2]);
//}

// Functions which were moved from GridPoint to DensityMap for memory savings

void DensityMap::initializeGradient(GridPoint & gridPoint){
		setddxPositiveXGradient(gridPoint,  0);
		setddyPositiveXGradient(gridPoint,  0);
		setddzPositiveXGradient(gridPoint,  0);
		setddxPositiveYGradient(gridPoint,  0);
		setddyPositiveYGradient(gridPoint,  0);
		setddzPositiveYGradient(gridPoint,  0);
		setddxPositiveZGradient(gridPoint,  0);
		setddyPositiveZGradient(gridPoint,  0);
		setddzPositiveZGradient(gridPoint,  0);
                //setFirstQuadrantGradient(gridPoint,  Vec3(0,0,0));
	}

void DensityMap::initialize(GridPoint & gridPoint){
		initializeGradient(gridPoint);
                //setFirstQuadrantGradient(gridPoint,Vec3(0));
		gridPoint.density = 0; 
		gridPoint.noise = 0; 
		gridPoint.position = Vec3(0);	
                //cout<<__FILE__<<":"<<__LINE__<<" Set firstQuadrantGradient =  "<<fetchFirstQuadrantGradient(gridPoint)<<"  position = "<<gridPoint.position <<" density = "<<gridPoint.density <<endl;
	}

void DensityMap::validatePosition(GridPoint & gridPoint, const Vec3 &myPosition) const{
               
          	//cout<<__FILE__<<":"<<__LINE__<<" about to validate position = "<<myPosition<<endl;
                //ValidateVec3( myPosition);
	}

void DensityMap::validateDensity(GridPoint & gridPoint, double myDensity) const{
		ValidateDouble(myDensity);
	
}

void DensityMap::validate(GridPoint & gridPoint) const{
		validatePosition(gridPoint,gridPoint.position);
		validateDensity(gridPoint,gridPoint.density);
		// write this code later.
	}

void DensityMap::setDensity(GridPoint & gridPoint,double myDensity)	{
    validateDensity(gridPoint,myDensity); 
    gridPoint.density = myDensity;	
    //cout<<__FILE__<<":"<<__LINE__<<":"<<__FUNCTION__<<" Just set density to "<<gridPoint.density<<" for grid point at position "<<gridPoint.position<<endl;
}

void DensityMap::setPosition(GridPoint & gridPoint, const Vec3 &myPosition)	{
 		validatePosition(gridPoint,myPosition);
		gridPoint.position = myPosition;	
	}

/*
Quadrant DensityMap::calcQuadrant(GridPoint & gridPoint, Vec3 queryPosition) const			{
		Quadrant tempQuadrant;    
		if (queryPosition[0] < gridPoint.position[0]) {tempQuadrant.positiveX = true;} else { tempQuadrant.positiveX = true;}		
		if (queryPosition[1] < gridPoint.position[1]) {tempQuadrant.positiveY = true;} else { tempQuadrant.positiveY = true;}
		if (queryPosition[2] < gridPoint.position[2]) {tempQuadrant.positiveZ = true;} else { tempQuadrant.positiveZ = true;}
		return tempQuadrant;
	}
Vec3 DensityMap::fetchGradient(GridPoint & gridPoint,Vec3 queryPosition) const {
		Quadrant tempQuadrant = calcQuadrant(gridPoint,queryPosition);
		Vec3 myGradient;
          	MMBLOG_FILE_FUNC_LINE(CRITICAL, " this function is obsolete, fetchGradient doesn't work anymore"<<endl;
		return myGradient;

	}
*/
Vec3 DensityMap::fetchFirstQuadrantGradient(GridPoint & gridPoint) const {
                return gridPoint.firstQuadrantGradient;
}



double DensityMap::getDensity(GridPoint & gridPoint) const	{return gridPoint.density ;	}
//double DensityMap::getDensity(GridPoint & gridPoint) const	{return gridPoint.density ;	}
	
double DensityMap::getDensity( GridPoint & gridPoint, Vec3 queryPosition) const	{
                Vec3 myVectorToGridMap = queryPosition - gridPoint.position;
                MMBLOG_FILE_FUNC_LINE(DEBUG, endl);
                Vec3 myGradient = calcInterpolatedFirstQuadrantGradient(gridPoint,  queryPosition); 
                double myDensity =  (gridPoint.density + (myGradient[0]* myVectorToGridMap[0])  + (myGradient[1]* myVectorToGridMap[1]) + (myGradient[2]* myVectorToGridMap[2]));
                //MMBLOG_FILE_FUNC_LINE(" Returning interpolated density of : "<<myDensity<<std::endl;
                return myDensity;
}




void DensityMap::setPositiveXGradient(GridPoint & gridPoint,Real myPositiveXGradient) {  gridPoint.firstQuadrantGradient[0] = myPositiveXGradient; }
void DensityMap::setddxPositiveXGradient(GridPoint & gridPoint,Real value) { gridPoint.ddxPositiveXGradient = (float)value; }
void DensityMap::setddyPositiveXGradient(GridPoint & gridPoint,Real value) { gridPoint.ddyPositiveXGradient = (float)value; }
void DensityMap::setddzPositiveXGradient(GridPoint & gridPoint,Real value) { gridPoint.ddzPositiveXGradient = (float)value; }
void DensityMap::setPositiveYGradient(GridPoint & gridPoint, Real myPositiveYGradient) { gridPoint.firstQuadrantGradient[1] = myPositiveYGradient; }
void DensityMap::setddxPositiveYGradient(GridPoint & gridPoint,Real value) { gridPoint.ddxPositiveYGradient = (float)value; }
void DensityMap::setddyPositiveYGradient(GridPoint & gridPoint,Real value) { gridPoint.ddyPositiveYGradient = (float)value; }
void DensityMap::setddzPositiveYGradient(GridPoint & gridPoint,Real value) { gridPoint.ddzPositiveYGradient = (float)value; }
void DensityMap::setPositiveZGradient(GridPoint & gridPoint,Real myPositiveZGradient) {  gridPoint.firstQuadrantGradient[2] = myPositiveZGradient; }
void DensityMap::setddxPositiveZGradient(GridPoint & gridPoint,Real value) { gridPoint.ddxPositiveZGradient = (float)value; }
void DensityMap::setddyPositiveZGradient(GridPoint & gridPoint,Real value) { gridPoint.ddyPositiveZGradient = (float)value; }
void DensityMap::setddzPositiveZGradient(GridPoint & gridPoint,Real value) { gridPoint.ddzPositiveZGradient = (float)value; }

void DensityMap::printSecondDerivatives(GridPoint & gridPoint) const {
    MMBLOG_FILE_FUNC_LINE(INFO, "Printing gridPoint.ddyPositiveXGradient, gridPoint.ddzPositiveXGradient : "
        <<gridPoint.ddyPositiveXGradient <<" "<<gridPoint.ddzPositiveXGradient<<endl);
    MMBLOG_FILE_FUNC_LINE(INFO, "Printing gridPoint.ddxPositiveYGradient, gridPoint.ddzPositiveYGradient : "
        <<gridPoint.ddxPositiveYGradient <<" "<<gridPoint.ddzPositiveYGradient<<endl);
    MMBLOG_FILE_FUNC_LINE(INFO, "Printing gridPoint.ddxPositiveZGradient, gridPoint.ddyPositiveZGradient : "
        <<gridPoint.ddxPositiveZGradient <<" "<<gridPoint.ddyPositiveZGradient<<endl);
}

void DensityMap::setNegativeXGradient(GridPoint & gridPoint,Real myNegativeXGradient) { 
        //MMBLOG_FILE_FUNC_LINE(CRITICAL, " this function is only for calculating the gradient in the first quadrant"<<endl;
        MMBLOG_FILE_FUNC_LINE(CRITICAL, ""); //Seriously????
}
void DensityMap::setNegativeYGradient(GridPoint & gridPoint,Real myNegativeYGradient) { 
        //MMBLOG_FILE_FUNC_LINE(CRITICAL, " this function is only for calculating the gradient in the first quadrant"<<endl;
        MMBLOG_FILE_FUNC_LINE(CRITICAL, ""); //Seriously????
}
void DensityMap::setNegativeZGradient(GridPoint & gridPoint,Real myNegativeZGradient) { 
        //MMBLOG_FILE_FUNC_LINE(CRITICAL, " this function is only for calculating the gradient in the first quadrant"<<endl;
        MMBLOG_FILE_FUNC_LINE(CRITICAL, ""); //Seriously????
}
/*void DensityMap::setFirstQuadrantGradient(GridPoint & gridPoint,Vec3 gradient){
        cout<<__FILE__<<":"<<__LINE__<<" Setting firstQuadrantGradient to "<<gradient<<" for grid point at position : "<<gridPoint.position <<endl;
        gridPoint.firstQuadrantGradient = gradient;
        MMBLOG_FILE_FUNC_LINE(CRITICAL, " this function is obsolete!"<<endl;

    }*/
Vec3 DensityMap::calcInterpolatedFirstQuadrantGradient(GridPoint & gridPoint, const Vec3 &queryPosition) const {
    //MMBLOG_FILE_FUNC_LINE(std::endl;
    MMBLOG_FILE_FUNC_LINE(DEBUG, endl);
    Vec3 dxdydz = unitCellParameters.convertFractionalVectorToFractionFromLowerLeft(unitCellParameters.convertCartesianVectorToFractionalVector(queryPosition)); //queryPosition - gridPoint.position; // the first term is the query position, the second term is the grid point position in cartesian space
    if ((dxdydz[0] < 0) || (dxdydz[1] <0 ) || (dxdydz[2] < 0)) { // see if we can make this trap unnecessary implicitly
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "this function is only for calculating the gradient in the first quadrant"<<endl);
    }
    if ((dxdydz[0] > 1.0) || (dxdydz[1] >1.0 ) || (dxdydz[2] > 1.0)) { // see if we can make this trap unnecessary implicitly
        MMBLOG_FILE_FUNC_LINE(INFO, " dxdydz = "<<dxdydz<<endl);
        MMBLOG_FILE_FUNC_LINE(
            CRITICAL,
            "this function is only for calculating the gradient in the first quadrant, and only within the cell in question. "
            "This value goes outside the voxel."<<endl
        );
    }
    //MMBLOG_FILE_FUNC_LINE(" dxdydz = "<<dxdydz<<std::endl;
    //printSecondDerivatives(gridPoint);
    Vec3 myGradient;
    myGradient[0] = fetchFirstQuadrantGradient(gridPoint)[0] ; //                                           + gridPoint.ddyPositiveXGradient*dxdydz[1] + gridPoint.ddzPositiveXGradient*dxdydz[2];
    myGradient[1] = fetchFirstQuadrantGradient(gridPoint)[1] ; // + gridPoint.ddxPositiveYGradient*dxdydz[0]                                           + gridPoint.ddzPositiveYGradient*dxdydz[2];
    myGradient[2] = fetchFirstQuadrantGradient(gridPoint)[2] ; // + gridPoint.ddxPositiveZGradient*dxdydz[0] + gridPoint.ddyPositiveZGradient*dxdydz[1];
    MMBLOG_FILE_FUNC_LINE(DEBUG, " myGradient = "<<myGradient <<endl);

    // Separated out ther second derivatives for debugging:
    /*if (0) {
        myGradient[0] +=  gridPoint.ddxPositiveXGradient*dxdydz[0] + gridPoint.ddyPositiveXGradient*dxdydz[1] + gridPoint.ddzPositiveXGradient*dxdydz[2];
        myGradient[1] +=  gridPoint.ddxPositiveYGradient*dxdydz[0] + gridPoint.ddyPositiveYGradient*dxdydz[1] + gridPoint.ddzPositiveYGradient*dxdydz[2];
        myGradient[2] +=  gridPoint.ddxPositiveZGradient*dxdydz[0] + gridPoint.ddyPositiveZGradient*dxdydz[1] + gridPoint.ddzPositiveZGradient*dxdydz[2];
        //return myGradient;
    }*/
    return myGradient;
}
