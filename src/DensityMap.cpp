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
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<< " There is a problem with the max- , min- (XYZ) or totalNumGrid (XYZ) or grid (XYZ) Spacing parameters "<< ((maxX-minX)/gridXSpacing +1) - totalNumGridX <<endl;	
        ErrorManager::instance.treatError();
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
			ErrorManager::instance <<__FILE__<<":"<<__LINE__<<	" No nearby grid point.  The point you requested "<< myGridIndices.getZGridIndex()<<" , "<<  myGridIndices.getYGridIndex()<<" , "<< myGridIndices.getXGridIndex() << " is off the density map.";
			ErrorManager::instance.treatError();
		} else {
			GridPoint & myGridPoint = ArrayOfGridPoints[myGridIndices.getZGridIndex()][myGridIndices.getYGridIndex()][myGridIndices.getXGridIndex()] ;
			//myGridPoint.validate();
			return myGridPoint;
		}
  		}
   
const GridPoint DensityMap::getGridPoint(GridIndices myGridIndices) const  {
		//return & ArrayOfGridPoints[0][0][0];
		return ArrayOfGridPoints[myGridIndices.getZGridIndex()][myGridIndices.getYGridIndex()][myGridIndices.getXGridIndex()] ;
  		}
   
void	    DensityMap::validateGridPoint(GridIndices myGridIndices){
			if (! hasGridPoint(myGridIndices)) {
				ErrorManager::instance <<__FILE__<<":"<<__LINE__<<	" No nearby grid point.  The point you requested is off the density map.";
				ErrorManager::instance.treatError();
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
GridIndices DensityMap::calcNearestGridIndices(Vec3 position)
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
		
GridIndices DensityMap::calcLowerLeftGridIndices(Vec3 position)
                {
                        iVec3 tempIndexVector = unitCellParameters.convertCartesianVectorToLowerIndexVector(position); //convertFractionalVectorToLowerIndexVector(position);
                        //int tempXIndex  = int (floor(((position[0]-minX) )/gridXSpacing));
                        //int tempYIndex  = int (floor(((position[1]-minY) )/gridYSpacing));
                        //int tempZIndex  = int (floor(((position[2]-minZ) )/gridZSpacing));

                        return GridIndices(tempIndexVector[0], tempIndexVector[1], tempIndexVector[2]); /// tempXIndex, tempYIndex, tempZIndex);

                }
		
const GridPoint DensityMap::getGridPoint(Vec3 myPosition)  {
	return getGridPoint(calcNearestGridIndices( myPosition));				
}

GridPoint & DensityMap::updGridPoint(Vec3 myPosition)  {
	return updGridPoint(calcNearestGridIndices( myPosition));				
}

const double DensityMap::getDensity(Vec3 myPosition) {
        //std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<"Taking myPosition = "<<myPosition<<std::endl;
        Vec3 myFractionalVector = unitCellParameters.convertCartesianVectorToFractionalVector(myPosition);
        //std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<< myFractionalVector[0] <<", "<<  myFractionalVector[1]   <<", "<< myFractionalVector[2]   <<" ..preceding should be unitCellParameters.convertCartesianVectorToFractionalVector(myPosition)"<<std::endl;
        if (!(unitCellParameters.fractionalVectorIsInsideMapBoundaries(myFractionalVector)))
	{
		return 0.0; //return zero density
	} else {
                GridIndices myLowerLeftGridIndices = calcLowerLeftGridIndices(myPosition);
		    return getDensity(  updGridPoint(myLowerLeftGridIndices), myPosition); // need to verify that this interpolated density is reasonable
	} 
}

//const double DensityMap::getDensity(SimTK::Vec3 myPosition) {
//    return getDensity(Vec3(myPosition[0], myPosition[1],myPosition[2] ));
//}
void DensityMap::initializeArrayOfGridPoints(){
        //validateGridParameters();
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
        unitCellParameters.validate();
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
        initializeVectorOfAmplitudeAndRandomPhases();
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
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
            ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Wrong number of grid points in Z direction! Found :"<< ArrayOfGridPoints.size()<<" expected : " << unitCellParameters.getNc()<<","<< unitCellParameters.getNc()<<endl;
            ErrorManager::instance.treatError();
        }
        for ( int zIndex = 0; zIndex < unitCellParameters.getNc(); zIndex++) {
        	if (ArrayOfGridPoints[zIndex].size() != unitCellParameters.getNb()) {
                   ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Wrong number of grid points in Y direction! Found :"<< ArrayOfGridPoints[zIndex].size()<<" expected : " << unitCellParameters.getNb()<<","<< unitCellParameters.getNc()<<endl;
                    ErrorManager::instance.treatError();
                }
        	for ( int yIndex = 0; yIndex < unitCellParameters.getNb(); yIndex++) {
        	    if (ArrayOfGridPoints[zIndex][yIndex].size() != (unitCellParameters.getNa()) ) {
        		             ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Wrong number of grid points in X direction! Found :"<< ArrayOfGridPoints[zIndex][yIndex].size()<<" expected : " << unitCellParameters.getNa()<<endl;
                             ErrorManager::instance.treatError();
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
    //    std::cout<<__FILE__<<":"<<__LINE__<<" You have specified a very high temperature : "<<temperature<<" .. returning white noise spectrum"<<std::endl;  
    //    return 1.;
    //}

    return (1 / (exp(frequency/temperature) - 1));
}



void DensityMap::resizeVectorOfAmplitudeAndRandomPhases(){ 
        std::cout<<__FILE__<<":"<<__LINE__<<" resizing Z component vectorOfAmplitudeFrequencyAndRandomPhases to "<<(unitCellParameters.calcMaxFrequencyDoublingsZ () +1)<<std::endl;
        vectorOfAmplitudeFrequencyAndRandomPhases.resize(unitCellParameters.calcMaxFrequencyDoublingsZ () +1);                                         
        std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
        for ( int zIndex = 0; zIndex <= unitCellParameters.calcMaxFrequencyDoublingsZ(); zIndex++) {
                std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
                std::cout<<__FILE__<<":"<<__LINE__<<" resizing  vectorOfAmplitudeFrequencyAndRandomPhases["<<zIndex <<"] to "<<(unitCellParameters.calcMaxFrequencyDoublingsZ () +1)<<std::endl;
        	vectorOfAmplitudeFrequencyAndRandomPhases[zIndex].resize(unitCellParameters.calcMaxFrequencyDoublingsY () +1);
        	for ( int yIndex = 0; yIndex <=  unitCellParameters.calcMaxFrequencyDoublingsY(); yIndex++) {
                    //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
                    //std::cout<<__FILE__<<":"<<__LINE__<<" resizing  vectorOfAmplitudeFrequencyAndRandomPhases["<<zIndex <<"]["<<yIndex<< "] to "<<(unitCellParameters.calcMaxFrequencyDoublingsX () +1)<<std::endl;
                    vectorOfAmplitudeFrequencyAndRandomPhases[zIndex][yIndex].resize(unitCellParameters.calcMaxFrequencyDoublingsX () +1);
        	}}
        std::cout<<__FILE__<<":"<<__LINE__<<std::endl;

        // Now make sure size was set correctly:
        if (vectorOfAmplitudeFrequencyAndRandomPhases.size() != (unitCellParameters.calcMaxFrequencyDoublingsZ()+1)) {
            ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Wrong number of grid points in Z direction! Found :"<< vectorOfAmplitudeFrequencyAndRandomPhases.size()<<" expected : " << unitCellParameters.calcMaxFrequencyDoublingsZ()+1             <<endl;
            ErrorManager::instance.treatError();
        }
        for ( int zIndex = 0; zIndex < unitCellParameters.calcMaxFrequencyDoublingsZ(); zIndex++) {
        	if (vectorOfAmplitudeFrequencyAndRandomPhases[zIndex].size() != (unitCellParameters.calcMaxFrequencyDoublingsY()+1) ) {
                   ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Wrong number of grid points in Y direction! Found :"<< vectorOfAmplitudeFrequencyAndRandomPhases[zIndex].size()<<" expected : " << unitCellParameters.calcMaxFrequencyDoublingsY() +1            <<endl;
                    ErrorManager::instance.treatError();
                }
        	for ( int yIndex = 0; yIndex < unitCellParameters.calcMaxFrequencyDoublingsY(); yIndex++) {
        	    if (vectorOfAmplitudeFrequencyAndRandomPhases[zIndex][yIndex].size() != (unitCellParameters.calcMaxFrequencyDoublingsX()+1 )) {
                             std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
        		     ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Check 1. Wrong number of grid points in X direction! Found :"<< vectorOfAmplitudeFrequencyAndRandomPhases[zIndex][yIndex].size()<<" expected : " <<unitCellParameters.calcMaxFrequencyDoublingsX()+1 <<endl;
                             ErrorManager::instance.treatError();
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
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" You have specified computeNoiseAutocorrelation = "<<computeNoiseAutocorrelation<<" and computeDensityAutocorrelation = "<< computeDensityAutocorrelation << " .. you must specify exactly one of these." <<endl;
        ErrorManager::instance.treatError();
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

    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
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
                                //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"distanceSquared, sumCorrelation, numCorrelationObservations, correlation, "<<autoCorrelationMapIterator->first<<","<<autoCorrelationMapIterator->second.sumOfObservations <<","<<autoCorrelationMapIterator->second.numOfObservations<<","<<autoCorrelationMapIterator->second.sumOfObservations / autoCorrelationMapIterator->second.numOfObservations<<std::endl;
             
                            }
                            else {
                                myAutoCorrelationStruct.sumOfObservations = correlation;
                                myAutoCorrelationStruct.numOfObservations = 1;
                                autoCorrelationMap.insert(std::make_pair(distanceSquared,myAutoCorrelationStruct));
                                //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"distanceSquared, sumCorrelation, numCorrelationObservations, correlation, "<<distanceSquared<<","<<myAutoCorrelationStruct.sumOfObservations  <<","<<myAutoCorrelationStruct.numOfObservations<<std::endl;
                            }

                            autoCorrelationMap.insert(std::make_pair(distanceSquared,myAutoCorrelationStruct));
                            //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<< " computeDensityAutocorrelation = " << computeDensityAutocorrelation <<" computeNoiseAutocorrelation = "<<computeNoiseAutocorrelation <<",  distanceSquared, correlation = ,"<<distanceSquared<<","<<correlation<<std::endl;
    }}}}}}
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
    std::cout<<std::string(40,'*')<<std::endl;
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"   Starting autocorrelation section. "<<std::endl;
    std::cout<<std::string(40,'*')<<std::endl;
    std::cout<<"distanceSquared, sumCorrelation, numCorrelationObservations, correlation "<<std::endl;
    for(autoCorrelationMapIterator = autoCorrelationMap.begin(); autoCorrelationMapIterator != autoCorrelationMap.end(); autoCorrelationMapIterator++){
        std::cout<<autoCorrelationMapIterator->first<<","<<autoCorrelationMapIterator->second.sumOfObservations <<","<<autoCorrelationMapIterator->second.numOfObservations<<","<<autoCorrelationMapIterator->second.sumOfObservations / autoCorrelationMapIterator->second.numOfObservations<<std::endl;
    }
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
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
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<< " total signal = "<<totalSignal<<", total noise in dense regions = "<<totalNoise<<", total noise everywhere = "<<newTotalNoiseEverywhere<<",  signalToNoiseRatio (masked, dense regions only) = "<< signalToNoiseRatio << std::endl;
    
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
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<< "  unitCellParameters.totalVolume() = "<<  unitCellParameters.totalVolume() << std::endl;
    for ( int xIndex = 0; xIndex < unitCellParameters.getNa();  xIndex++){
        
        //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<< " xIndex = "<<xIndex <<  std::endl;
        for ( int yIndex = 0; yIndex < unitCellParameters.getNb();  yIndex++){
            for ( int zIndex = 0; zIndex < unitCellParameters.getNc();  zIndex++){
                updGridPoint(GridIndices(xIndex,yIndex,zIndex)).noise=0.;
                for ( int xPhaseIndex = 1; xPhaseIndex <= unitCellParameters.calcMaxFrequencyDoublingsX ();  xPhaseIndex++){
                    for ( int yPhaseIndex = 1; yPhaseIndex <= unitCellParameters.calcMaxFrequencyDoublingsY () ;  yPhaseIndex++){
                        for ( int zPhaseIndex = 1; zPhaseIndex <= unitCellParameters.calcMaxFrequencyDoublingsZ () ;  zPhaseIndex++){
                            myAmpFreqRandPhase = vectorOfAmplitudeFrequencyAndRandomPhases[zPhaseIndex][yPhaseIndex][xPhaseIndex];
                            //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<< "  myAmpFreqRandPhase.phaseX = "<< myAmpFreqRandPhase.phaseX<<std::endl;
                            //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<< "  myAmpFreqRandPhase.phaseY = "<< myAmpFreqRandPhase.phaseY<<std::endl;
                            //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<< "  myAmpFreqRandPhase.phaseZ = "<< myAmpFreqRandPhase.phaseZ<<std::endl;
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
    //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<< " total signal = "<<averageSignal<<", total noise in dense regions = "<<averageNoise<<", total noise everywhere = "<<totalNoiseEverywhere<<",  signalToNoiseRatio (masked, dense regions only) = "<< signalToNoiseRatio << std::endl;
    normalizeNoiseMap(totalNoiseEverywhere);
    // moved to Repel.cpp:
    //densityAutocorrelation(1,0); // Arguments are : calculate noise correlation = 1, calculate density corrleation = 0
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<< std::endl;
    //densityAutocorrelation(0,1);
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<< std::endl;

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
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
        resizeVectorOfAmplitudeAndRandomPhases();
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
        //validateVectorOfAmplitudeFrequencyAndRandomPhasesSize();
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
        unitCellParameters.validate();
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
        Vec3 tempPosition(0,0,0);	
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
        AmplitudeFrequencyAndRandomPhases myAmpFreqPhase ;
        double myFrequency = 0.0;
        srand (time(NULL)); // initialize rand()
        // Just never access index 0 in x,y,z. Makes the math slightly easier.
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<< "  unitCellParameters.geta() = "<<unitCellParameters.geta()<<std::endl;
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<< "  unitCellParameters.getNa() = "<<unitCellParameters.getNa()<<std::endl;
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
                            //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<< "  myAmpFreqPhase.frequencyX = "<< myAmpFreqPhase.frequencyX <<std::endl;
                            //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<< "  myAmpFreqPhase.frequencyY = "<< myAmpFreqPhase.frequencyX <<std::endl;
                            //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<< "  myAmpFreqPhase.frequencyZ = "<< myAmpFreqPhase.frequencyX <<std::endl;
                            myFrequency=  (sqrt(myAmpFreqPhase.frequencyX*myAmpFreqPhase.frequencyX + myAmpFreqPhase.frequencyY*myAmpFreqPhase.frequencyY + myAmpFreqPhase.frequencyZ*myAmpFreqPhase.frequencyZ )) ;
                            myAmpFreqPhase.amplitude = plancksLaw(noiseTemperature, myFrequency);
                            //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<< "wavenumber doublings x,y,z, noiseTemperature, myFrequency,freqX,freqY,freqZ, plancksLaw(noiseTemperature, myFrequency) =              "<<","<< xIndex <<","<< yIndex <<","<<zIndex <<","  <<noiseTemperature<<","<<myFrequency<<",>"<<  myAmpFreqPhase.frequencyX<<","<<   myAmpFreqPhase.frequencyY  <<","<< myAmpFreqPhase.frequencyZ <<"<,"<<    myAmpFreqPhase.amplitude<<std::endl  ;
                            // random phases in radians:
                            myAmpFreqPhase.phaseX = rand() / (double)RAND_MAX * SimTK::Pi * 2.0; // Have to cast (double)RAND_MAX to prevent the result of the division from being cast as int.
                            myAmpFreqPhase.phaseY = rand() / (double)RAND_MAX * SimTK::Pi * 2. ;
                            myAmpFreqPhase.phaseZ = rand() / (double)RAND_MAX * SimTK::Pi * 2. ;
                            //myAmpFreqPhase.phaseX = 0.;
                            //myAmpFreqPhase.phaseY = 0.;
                            //myAmpFreqPhase.phaseZ = 0.;
                            //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" setting vectorOfAmplitudeFrequencyAndRandomPhases"<<zIndex<<","<<   yIndex <<","<<  xIndex<<std::endl;
                            vectorOfAmplitudeFrequencyAndRandomPhases[zIndex][yIndex][xIndex]=myAmpFreqPhase;  
                            //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
        }}}
        //validateVectorOfAmplitudeFrequencyAndRandomPhasesSize()
}


void DensityMap::loadParametersAndDensity(const String densityFileName)
{
    unsigned int extIndex = densityFileName.rfind(".");
    String extension = densityFileName.substr(extIndex);
    if(extension == ".xplor")
        loadParametersAndDensity_XPLOR(densityFileName);
    else if(extension == ".dx")
        { //loadParametersAndDensity_OpenDX(densityFileName);
            ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Unable to open density map : "<<densityFileName<<" .. DX is not a supported format at the moment."<<endl;
            ErrorManager::instance.treatError();
        }
    else if(extension == ".situs" || extension == ".sit")
        {//loadParametersAndDensity_Situs(densityFileName);
            ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Unable to open density map : "<<densityFileName<<" .. Situs is nota  supported format at the moment"<<endl;
            ErrorManager::instance.treatError();
        }
    else if (extension == ".map" || extension == ".ccp4" || extension == ".mrc" )
    {
#ifdef CPP4_MAPS_USAGE
        loadParametersAndDensity_CCP4MAP  ( densityFileName );
#else
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Unable to open density map : "<<densityFileName<<" .. CCP4 maps loading feature was not allowed when MMB installation was done on this machine."<<endl;
        ErrorManager::instance.treatError();
#endif
    }
    else
    {
        ErrorManager::instance << "DensityMap: Extension unknown for " << densityFileName << endl;
        ErrorManager::instance.treatError();
    }
}

#ifdef CPP4_MAPS_USAGE
/*! \brief Function responsible for reading in the CCP4's MAP format.

This function uses the CMAPLIB ccp4's library (the C version of the original fortran library) to read both, the
header and the data portions of the file format. It should be able to deal with both, the CCP4 MAP and the MRC
formats.
 
More specifically, the function firstly opens the supplied MAP file, proceeding to read in the cell information (float
array of 6 numbers, x, y and z dimensions in Angstroms and angles a, b and c) from the file header. It then continues
to read the dims (int array of 3 with the x, y and z axes sizes in indices) from the file header. Then, it reads the origin
(int array of 3 with x, y and z axes origin positions) and the axis order (int array of 3 with the x, y and z axes order given
by values 0, 1 and 2) from the file header.
 
Next, this function fills in the unitCellParameters structure with the correct indexing information, index distances and
angles. Subsequently, the function calls the DeOrthogonalisation matrix computation function and the grid preparing
functions. Finally, this function reads in the data portion of the file, taking into account the axis orders are detected
from the file header and saving the densities to the GridIndices object. It then closes the MAP file and terminates the
function.
 
\warning This function requires the ccp4's CMAPLIB library and will not work unless this is installed and properly linked.
\warning This function does only limitted reporting to the log, please add appropriate log outputs for the procedure.
\warning This function needs more testing - it was tested only for a single case, please do more testing and then remore this warning.

\param[in] densityFileName String containing the path to the density file which should read in.
*/
void DensityMap::loadParametersAndDensity_CCP4MAP(const String densityFileName) {
    //======================================== Create CCP4MAPfile object
    CMap_io::CMMFile *mapFile                 = NULL;
    
    //======================================== Set local variables
    int myMapMode                             = O_RDONLY;
    float xDim, yDim, zDim, aAng, bAng, cAng; // The cell dimensions in Angstroms and the angles of the cell axes in degrees
    int xInds, yInds, zInds;                  // The number of indices along the three axes.
    int xAxisOrder, yAxisOrder, zAxisOrder;   // The order of the axes (XYZ, YXZ, ...)
    int xOrigin, yOrigin, zOrigin;            // The origin of the indexing system.
    
    //======================================== Open file for reading and check if it can be opened
    mapFile                                   = reinterpret_cast<CMap_io::CMMFile*> ( CMap_io::ccp4_cmap_open ( densityFileName.c_str() , myMapMode ) );
    if ( mapFile == NULL )
    {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Unable to open density map : "<<densityFileName<<endl;
        ErrorManager::instance.treatError();
    }
    
    //======================================== Read in the cell information
    {
         //=================================== Initialise variables
         float *cell                          = NULL;
         
         //=================================== Allocate memory
         cell                                 = (float*) malloc ( 6 * sizeof ( float ) );
         
         //=================================== Check memory allocation
         if ( cell == NULL ) { ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Failed to allocate memory." << endl; ErrorManager::instance.treatError(); }
         
         //=================================== Read the data
         CMap_io::ccp4_cmap_get_cell          ( mapFile, cell );
         
         //=================================== Save the data
         xDim                                 = cell[0];
         yDim                                 = cell[1];
         zDim                                 = cell[2];
         aAng                                 = cell[3];
         bAng                                 = cell[4];
         cAng                                 = cell[5];
         
         //=================================== Release memory
         free                                 ( cell );
    }
    
    //======================================== Read in the dimensions information
    {
         //=================================== Initialise variables
         int *dim                             = NULL;
         
         //=================================== Allocate memory
         dim                                  = (int*) malloc (3 * sizeof ( int ) );
         
         //=================================== Check memory allocation
         if ( dim == NULL ) { ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Failed to allocate memory." << endl; ErrorManager::instance.treatError(); }
         
         //=================================== Read the data
         CMap_io::ccp4_cmap_get_dim           ( mapFile, dim );
         
         //=================================== Save the data
         xInds                                = dim[0];
         yInds                                = dim[1];
         zInds                                = dim[2];
         
         //=================================== Release memory
         free                                 ( dim );
    }
    
    //======================================== Read in the system origin information
    {
         //=================================== Initialise variables
         int *origin                          = NULL;
         
         //=================================== Allocate memory
         origin                               = (int*) malloc (3 * sizeof ( int ) );
         
         //=================================== Check memory allocation
         if ( origin == NULL ) { ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Failed to allocate memory." << endl; ErrorManager::instance.treatError(); }
         
         //=================================== Read the data
         CMap_io::ccp4_cmap_get_origin        ( mapFile, origin );
         
         //=================================== Save the data
         xOrigin                              = origin[0];
         yOrigin                              = origin[1];
         zOrigin                              = origin[2];
         
         //=================================== Release memory
         free                                 ( origin );
    }
    
    //======================================== Read in the axes order information
    {
         //=================================== Initialise variables
         int *order                           = NULL;
         
         //=================================== Allocate memory
         order                                = (int*) malloc (3 * sizeof ( int ) );
         
         //=================================== Check memory allocation
         if ( order == NULL ) { ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Failed to allocate memory." << endl; ErrorManager::instance.treatError(); }
         
         //=================================== Read the data
         CMap_io::ccp4_cmap_get_order         ( mapFile, order );
         
         //=================================== Save the data
         xAxisOrder                           = order[0];
         yAxisOrder                           = order[1];
         zAxisOrder                           = order[2];
         
         //=================================== Release memory
         free                                 ( order );
    }
    
    //======================================== Set the N for unitCellParameters object
    unitCellParameters.setN                   ( xInds,                   // Number of atoms along the x-axis
                                                xOrigin,                 // X-axis from index
                                                xOrigin + xInds - 1,     // X-axis to index
                                                yInds,                   // Number of atoms along the y-axis
                                                yOrigin,                 // Y-axis from index
                                                yOrigin + yInds - 1,     // Y-axis to index
                                                zInds,                   // Number of atoms along the z-axis
                                                zOrigin,                 // Z-axis from index
                                                zOrigin + zInds - 1 );   // Z-axis to index
    
    //======================================== Set dimensions and angles for unitCellParameters object
    unitCellParameters.setabc                 ( ( xDim / ( xInds - 1 ) ) / 10.0 ,  // Distance between two indices along the x-axis in nm ( divide by 10 to get nm from A )
                                                ( yDim / ( yInds - 1 ) ) / 10.0 ,  // Distance between two indices along the y-axis in nm ( divide by 10 to get nm from A )
                                                ( zDim / ( zInds - 1 ) ) / 10.0 ); // Distance between two indices along the z-axis in nm ( divide by 10 to get nm from A )
    unitCellParameters.setAlphaUsingDegrees   ( aAng );
    unitCellParameters.setBetaUsingDegrees    ( bAng );
    unitCellParameters.setGammaUsingDegrees   ( cAng );
    
    //======================================== De-orthoginalisation matrix
    unitCellParameters.setDeOrthogonalizationMatrix ( );
    
    //======================================== Prepate grid
    initializeArrayOfGridPoints               ( );
    
    //======================================== Read in the densities
    {
        //==================================== Check the map mode
        int mapMode                           = CMap_io::ccp4_cmap_get_datamode ( mapFile );
        if ( ( mapMode != 0 ) && ( mapMode != 2 ) )
        {
            ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Allowed map modes are 0 or 2, supplied map mode is: " << mapMode << endl; ErrorManager::instance.treatError();
        }
        
        //==================================== Initialise local variables
        int index;
        int iters[3];
        int maxLim[3];
        int XYZOrder[3];
        int newU, newV, newW;
        int arrPos;
        int maxMapV = 0, maxMapW = 0;
        int order[3]; order[0] = xAxisOrder;    order[1] = yAxisOrder;    order[2] = zAxisOrder;
        int orig[3];  orig[0]  = xOrigin;       orig[1]  = yOrigin;       orig[2]  = zOrigin;
        int dim[3];   dim[0]   = xInds;         dim[1]   = yInds;         dim[2]   = zInds;
        
        //==================================== Set the dimensions for XYZ indexing
        if ( order[1] == 1 ) { maxMapV = dim[0] - 1; }
        if ( order[1] == 2 ) { maxMapV = dim[1] - 1; }
        if ( order[1] == 3 ) { maxMapV = dim[2] - 1; }
        
        if ( order[2] == 1 ) { maxMapW = dim[0] - 1; }
        if ( order[2] == 2 ) { maxMapW = dim[1] - 1; }
        if ( order[2] == 3 ) { maxMapW = dim[2] - 1; }
        
        //==================================== Set the XYZ indexing order and indices
        for ( int iter = 0; iter < 3; iter++ )
        {
            maxLim[iter]                      = orig[iter] + dim[iter] - 1;
            XYZOrder[order[iter]-1]           = iter;
        }
        
        //==================================== Solve the dimensions and sizes for reading
        int fastDimSize                       = ( maxLim[0] - orig[0] + 1 );
        int midDimSize                        = ( maxLim[1] - orig[1] + 1 ) * fastDimSize;
        std::vector < float > section         ( midDimSize );
        
        //==================================== Read in the map data from the correct axis order, which can be any.
        for ( iters[2] = orig[2]; iters[2] <= maxLim[2]; iters[2]++ )
        {
            index                             = 0;
            CMap_io::ccp4_cmap_read_section( mapFile, &section[0] );
            
            //================================ Deal with mode 0
            if ( mapMode == 0 )
            {
                for ( int iter = ( midDimSize - 1 ); iter >= 0; iter-- )
                {
                    section[iter]             = static_cast < float > ( ( reinterpret_cast<unsigned char*> (&section[0]) )[iter] );
                }
            }
            
            for ( iters[1] = orig[1]; iters[1] <= maxLim[1]; iters[1]++ )
            {
                for ( iters[0] = orig[0]; iters[0] <= maxLim[0]; iters[0]++ )
                {
                    //======================== Compute the correct XYZ positions
                    newU                      = iters[XYZOrder[0]] - orig[XYZOrder[0]];
                    newV                      = iters[XYZOrder[1]] - orig[XYZOrder[1]];
                    newW                      = iters[XYZOrder[2]] - orig[XYZOrder[2]];
                    
                    //======================== Save the read in density into the correct MMB grid position
                    GridIndices centralGridIndices ( newU,  newV, newW );
                    GridPoint & centralGridPoint = updGridPoint ( centralGridIndices );
                    setDensity                ( centralGridPoint, static_cast < float > ( section[ index++ ] ) );
                }
            }
        }
    }
    
    //======================================== Close the map file
    CMap_io::ccp4_cmap_close                  ( mapFile );

    //======================================== DONE!
    return ;
}
#endif

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
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
        unitCellParameters.setDeOrthogonalizationMatrix ();
        //unitCellParameters.setOrthogonalizationMatrix ();
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
        initializeArrayOfGridPoints();
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
        mystring = readAndParseLine  (inFile);
        if ((mystring[0]).compare("ZYX") != 0    ) {ErrorManager::instance <<__FILE__<<":"<<__LINE__<<"expected ZYX, got :"<<mystring[0]<<endl; ErrorManager::instance.treatError();}


        for ( int zIndex = 0; zIndex < unitCellParameters.getNc(); zIndex ++) {
        	// read z-index
                //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" unitCellParameters.getNc() " <<unitCellParameters.getNc()<<std::endl;
                mystring = readAndParseLine  (inFile);
        	if (atoi ((mystring[0].c_str())) != zIndex) {
        			ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Something is wrong with the input file reading.  expected to read zIndex "<<zIndex<< " and instead read : "<<mystring[0]<<endl;
        		ErrorManager::instance.treatError();
        	}	
        	

        	for ( int yIndex = 0; yIndex < unitCellParameters.getNb(); yIndex ++) 
        	for ( int xIndex = 0; xIndex < unitCellParameters.getNa(); xIndex = xIndex + 0) 

        	{
                        //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
        		mystring = readAndParseOnColWidth (inFile,12);
        		if ((int)mystring.size() > densitiesPerLine) {
        			ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Too many densities on this line!  expected "<<densitiesPerLine<< "and found "<<mystring.size()<<endl;
        			ErrorManager::instance.treatError();
        		}
        		//xIndex --;
        		for (int i = 0; i<(int)mystring.size() ;i++) {
        			if (xIndex == unitCellParameters.getNa()){xIndex = 0; yIndex++; }
        			if ((xIndex+1 + (yIndex*unitCellParameters.getNa()) ) <= (unitCellParameters.getNa()* unitCellParameters.getNb()))
        			{
        				GridIndices centralGridIndices((xIndex ),  yIndex, zIndex);
        				GridPoint & centralGridPoint = updGridPoint(centralGridIndices);
                                        //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" For point xIndex, yIndex, zIndex : "<<  xIndex<<", "<< yIndex<<", "<< zIndex<<" setting density to "<<mystring[i].c_str()<<std::endl;

        				setDensity(centralGridPoint, atof(mystring[i].c_str()));	
        			} else{ 
        				ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Something went wrong .. too many densities !  "<<endl;
        				ErrorManager::instance.treatError();
        		        }
        			xIndex++;
        		}
        		
        	}
        } // of zIndex
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
        inFile.close();
}

void incrementer(int & myIndex, const int maxIndex) {
    if (myIndex < maxIndex) myIndex ++;
    if (myIndex == maxIndex) {myIndex = 0;}
    if (myIndex > maxIndex) {
            ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Unexplained error! index of  "<<myIndex<<" exceeds maximum of "<<maxIndex<<endl;
            ErrorManager::instance.treatError();
    }
}

void incrementer(int & xIndex, const int maxX, int & yIndex, const int maxY, int & zIndex, const int maxZ, int & counter, const int maxCounter){
    //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" xIndex, maxX : "<<xIndex<<", "<<maxX <<" yIndex, maxY : "<<yIndex<<", "<<maxY <<" zIndex, maxZ : "<<zIndex<<", "<<maxZ<<std::endl;
    incrementer(xIndex, maxX);
    if (xIndex == 0) {incrementer(yIndex,maxY);} // if x rolls over, then increment y
    if ((yIndex == 0) && (xIndex == 0)) {
        incrementer(zIndex, maxZ);
    } // if x and y both roll over, increment z
    incrementer(counter, maxCounter); // Always increment counter
}

void DensityMap::writeDensityMapXplor(const String densityFileName, const bool writeDensity, const bool writeNoise ) {
        ofstream outFile(densityFileName.c_str(),ofstream::out);
        int densitiesPerLine = 6;	
        vector<String> mystring;
       	stringstream u;
       	String tempString;
        if (! (outFile.is_open())) {
            ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Unable to open density map : "<<densityFileName<<endl;
            ErrorManager::instance.treatError();
        }
        outFile <<std::endl; // write blank line out.
        outFile << "0 !NTITLE"<< std::endl; // Number of lines to skip
        
        outFile <<unitCellParameters.getNa()<<" "<< unitCellParameters. getaMin()<<" "<< unitCellParameters.getaMax()<<" " <<unitCellParameters.getNb()<<" "<< unitCellParameters.getbMin()<<" "<< unitCellParameters.getbMax()<<" " <<unitCellParameters.getNc()<<" "<< unitCellParameters.getcMin()<<" "<< unitCellParameters.getcMax()<<std::endl;
        outFile << 10*unitCellParameters.geta()*(unitCellParameters.getNa() - 1) <<" " << 10*unitCellParameters.getb()*(unitCellParameters.getNb() - 1)<<" " <<10*unitCellParameters.getc()*(unitCellParameters.getNc() - 1) << " " << unitCellParameters.getAlpha() / SimTK::Pi * 180.0  << " " << unitCellParameters.getBeta() / SimTK::Pi * 180.0  << " " << unitCellParameters.getGamma() / SimTK::Pi * 180.0  <<std::endl;
        
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
        outFile << "ZYX"; //<<std::endl;

        int xIndex = 0; 
        int yIndex = 0; 
        int zIndex = 0;
        int counter = 0;
                 bool keepGoing = 1;
                 while(keepGoing)
                 {
                     //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" xIndex = "<<xIndex<<" yIndex = "<<yIndex<<" zIndex = "<<zIndex<<std::endl;
                     if ((yIndex == 0) && (xIndex == 0)) {
                         if (counter !=0) {outFile <<std::endl; // Carriage return. To terminate the incomplete density line.     .
                             counter = 0;
                         }
                         std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" xIndex = "<<xIndex<<" yIndex = "<<yIndex<<" zIndex = "<<zIndex<<std::endl;
                         outFile <<std::endl; // Carriage return.
                         outFile <<zIndex<<std::endl;
                     } else if (counter == 0 ){
                         //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" xIndex = "<<xIndex<<" yIndex = "<<yIndex<<" zIndex = "<<zIndex<<std::endl;
                         outFile <<std::endl; // Carriage return. but only if we did not just print zIndex above.
                     }
                     //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" xIndex = "<<xIndex<<" yIndex = "<<yIndex<<" zIndex = "<<zIndex<<std::endl;
                     // Write a density:
                     GridIndices centralGridIndices((xIndex ),  yIndex, zIndex);
                     GridPoint & centralGridPoint = updGridPoint(centralGridIndices);
                     outFile <<" "<<setw(11)<< (centralGridPoint.noiseFreeDensity * writeDensity) + (centralGridPoint.noise * writeNoise);	// One density for each 12-characgter column. Depending on parameters passed, one can write the original density, the noise, or the sum of the two.
                     // now increment indices:
                     incrementer(xIndex, unitCellParameters.getNa(), yIndex, unitCellParameters.getNb(), zIndex , unitCellParameters.getNc(), counter, 6);
                      

                     if ((zIndex == 0) && (yIndex == 0) && (xIndex == 0)) {
                         //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" xIndex = "<<xIndex<<" yIndex = "<<yIndex<<" zIndex = "<<zIndex<<std::endl;
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
					setPositiveXGradient (centralGridPoint,(getDensity(updGridPoint(GridIndices(xIndex+1,  yIndex, zIndex))) -getDensity(centralGridPoint))/ unitCellParameters.geta());
				} else if (xIndex == (unitCellParameters.getNa()-1)) { // if grid point is at the +X boundary
					setPositiveXGradient (centralGridPoint,(0. -getDensity(centralGridPoint)) / unitCellParameters.geta() );
                                }

				if (yIndex < (unitCellParameters.getNb()-1)) {// if grid point is not at the +Y boundary
					setPositiveYGradient (centralGridPoint,(getDensity(updGridPoint(GridIndices(xIndex,  yIndex+1, zIndex))) -getDensity(centralGridPoint))/ unitCellParameters.getb());
				} else if (yIndex == (unitCellParameters.getNb()-1)) {
					setPositiveYGradient (centralGridPoint,(0. -getDensity(centralGridPoint)) / unitCellParameters.getb() );
                                }



				if (zIndex < (unitCellParameters.getNc()-1)) {// if grid point is not at the +Z boundary:
					setPositiveZGradient (centralGridPoint,(getDensity(updGridPoint(GridIndices(xIndex,  yIndex, zIndex+1))) -getDensity(centralGridPoint))/ unitCellParameters.getc());
				
				} else if (zIndex == (unitCellParameters.getNc()-1)) {
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

Vec3 DensityMap::fetchFirstQuadrantGradient(Vec3 position)  {

                        GridIndices myLowerLeftGridIndex = calcLowerLeftGridIndices(   position);
                         if (hasGridPoint(myLowerLeftGridIndex)) {
                                 return fetchFirstQuadrantGradient(ArrayOfGridPoints[myLowerLeftGridIndex.getZGridIndex()][myLowerLeftGridIndex.getYGridIndex()][myLowerLeftGridIndex.getXGridIndex()]);

                         } else if (hasGridPoint(GridIndices (myLowerLeftGridIndex.getXGridIndex()+1, myLowerLeftGridIndex.getYGridIndex(),myLowerLeftGridIndex.getZGridIndex()))) {
                                GridPoint myGridPoint;
                                initialize(myGridPoint);
                                setPositiveXGradient(myGridPoint,(getDensity(updGridPoint(GridIndices(myLowerLeftGridIndex.getXGridIndex() +1,  myLowerLeftGridIndex.getYGridIndex() , myLowerLeftGridIndex.getZGridIndex() ))) - 0.) / unitCellParameters.geta()) ;
                                return fetchFirstQuadrantGradient(myGridPoint) ;

                         } else if (hasGridPoint(GridIndices (myLowerLeftGridIndex.getXGridIndex(), myLowerLeftGridIndex.getYGridIndex()+1,myLowerLeftGridIndex.getZGridIndex()))) {
                                GridPoint myGridPoint;
                                initialize(myGridPoint);
                                setPositiveYGradient(myGridPoint,(getDensity(updGridPoint(GridIndices(myLowerLeftGridIndex.getXGridIndex() ,  myLowerLeftGridIndex.getYGridIndex() +1, myLowerLeftGridIndex.getZGridIndex() ))) - 0.) / unitCellParameters.getb()) ;
                                return fetchFirstQuadrantGradient(myGridPoint) ;

                         } else if (hasGridPoint(GridIndices (myLowerLeftGridIndex.getXGridIndex(), myLowerLeftGridIndex.getYGridIndex(),myLowerLeftGridIndex.getZGridIndex()+1))) {
                                GridPoint myGridPoint;
                                initialize(myGridPoint);
                                setPositiveZGradient(myGridPoint,(getDensity(updGridPoint(GridIndices(myLowerLeftGridIndex.getXGridIndex() ,  myLowerLeftGridIndex.getYGridIndex() , myLowerLeftGridIndex.getZGridIndex() +1 ))) - 0.) / unitCellParameters.geta()) ;
                                // return myGridPoint.fetchGradient (position) ;
                                return fetchFirstQuadrantGradient(myGridPoint) ;

                         } else {
                                 Vec3 tempVec3(0);
                                 return tempVec3;
                         }

}

Vec3 DensityMap::calcInterpolatedFirstQuadrantGradient(Vec3 position)  {

                        GridIndices myLowerLeftGridIndex = calcLowerLeftGridIndices(   position);
                         if (hasGridPoint(myLowerLeftGridIndex)) {
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

void DensityMap::validatePosition(GridPoint & gridPoint, Vec3 myPosition) const{
               
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

void DensityMap::setPosition(GridPoint & gridPoint, Vec3 myPosition)	{
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
          	ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" this function is obsolete, fetchGradient doesn't work anymore"<<endl;
        	ErrorManager::instance.treatError();
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
                Vec3 myGradient = calcInterpolatedFirstQuadrantGradient(gridPoint,  queryPosition); 
                double myDensity =  (gridPoint.density + (myGradient[0]* myVectorToGridMap[0])  + (myGradient[1]* myVectorToGridMap[1]) + (myGradient[2]* myVectorToGridMap[2]));
                //std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<" Returning interpolated density of : "<<myDensity<<std::endl;
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
    std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<" Printing gridPoint.ddyPositiveXGradient, gridPoint.ddzPositiveXGradient : "
        <<gridPoint.ddyPositiveXGradient <<" "<<gridPoint.ddzPositiveXGradient<<std::endl;
    std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<" Printing gridPoint.ddxPositiveYGradient, gridPoint.ddzPositiveYGradient : "
        <<gridPoint.ddxPositiveYGradient <<" "<<gridPoint.ddzPositiveYGradient<<std::endl;
    std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<" Printing gridPoint.ddxPositiveZGradient, gridPoint.ddyPositiveZGradient : "
        <<gridPoint.ddxPositiveZGradient <<" "<<gridPoint.ddyPositiveZGradient<<std::endl;
}

void DensityMap::setNegativeXGradient(GridPoint & gridPoint,Real myNegativeXGradient) { 
        //ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" this function is only for calculating the gradient in the first quadrant"<<endl;
        ErrorManager::instance.treatError();}
void DensityMap::setNegativeYGradient(GridPoint & gridPoint,Real myNegativeYGradient) { 
        //ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" this function is only for calculating the gradient in the first quadrant"<<endl;
        ErrorManager::instance.treatError();}
void DensityMap::setNegativeZGradient(GridPoint & gridPoint,Real myNegativeZGradient) { 
        //ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" this function is only for calculating the gradient in the first quadrant"<<endl;
        ErrorManager::instance.treatError();}
/*void DensityMap::setFirstQuadrantGradient(GridPoint & gridPoint,Vec3 gradient){
        cout<<__FILE__<<":"<<__LINE__<<" Setting firstQuadrantGradient to "<<gradient<<" for grid point at position : "<<gridPoint.position <<endl;
        gridPoint.firstQuadrantGradient = gradient;
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" this function is obsolete!"<<endl;
        ErrorManager::instance.treatError();

    }*/
Vec3 DensityMap::calcInterpolatedFirstQuadrantGradient(GridPoint & gridPoint,Vec3 queryPosition) const {
    //std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<std::endl;
    Vec3 dxdydz = unitCellParameters.convertFractionalVectorToFractionFromLowerLeft(unitCellParameters.convertCartesianVectorToFractionalVector(queryPosition)); //queryPosition - gridPoint.position; // the first term is the query position, the second term is the grid point position in cartesian space
    if ((dxdydz[0] < 0) || (dxdydz[1] <0 ) || (dxdydz[2] < 0)) { // see if we can make this trap unnecessary implicitly

        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" this function is only for calculating the gradient in the first quadrant"<<endl;
        ErrorManager::instance.treatError();
    }
    if ((dxdydz[0] > 1.0) || (dxdydz[1] >1.0 ) || (dxdydz[2] > 1.0)) { // see if we can make this trap unnecessary implicitly
        std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<" dxdydz = "<<dxdydz<<std::endl;
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" this function is only for calculating the gradient in the first quadrant, and only within the cell in question. This value goes outside the voxel."<<endl;
        ErrorManager::instance.treatError();
    }
    //std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<" dxdydz = "<<dxdydz<<std::endl;
    //printSecondDerivatives(gridPoint);
    Vec3 myGradient;
    myGradient[0] = fetchFirstQuadrantGradient(gridPoint)[0] ; //                                           + gridPoint.ddyPositiveXGradient*dxdydz[1] + gridPoint.ddzPositiveXGradient*dxdydz[2];
    myGradient[1] = fetchFirstQuadrantGradient(gridPoint)[1] ; // + gridPoint.ddxPositiveYGradient*dxdydz[0]                                           + gridPoint.ddzPositiveYGradient*dxdydz[2];
    myGradient[2] = fetchFirstQuadrantGradient(gridPoint)[2] ; // + gridPoint.ddxPositiveZGradient*dxdydz[0] + gridPoint.ddyPositiveZGradient*dxdydz[1];

    // Separated out ther second derivatives for debugging:
    if (0) {
    myGradient[0] +=  gridPoint.ddxPositiveXGradient*dxdydz[0] + gridPoint.ddyPositiveXGradient*dxdydz[1] + gridPoint.ddzPositiveXGradient*dxdydz[2];
    myGradient[1] +=  gridPoint.ddxPositiveYGradient*dxdydz[0] + gridPoint.ddyPositiveYGradient*dxdydz[1] + gridPoint.ddzPositiveYGradient*dxdydz[2];
    myGradient[2] +=  gridPoint.ddxPositiveZGradient*dxdydz[0] + gridPoint.ddyPositiveZGradient*dxdydz[1] + gridPoint.ddzPositiveZGradient*dxdydz[2];
    return myGradient;
    }
}

/*Vec3 DensityMap::getFirstQuadrantGradient(GridPoint & gridPoint){
    return fetchFirstQuadrantGradient(gridPoint);
}*/

// Expects OpenDX dx density maps with Angstrom as spatial unit.
// Adapted from BioSpring http://sourceforge.net/projects/biospring/
/*
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
        unitCellParameters.getNb() = sizej;
        unitCellParameters.getNc() = sizek;
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
        unitCellParameters.geta() = delta[0][0]/10.0;
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
    cout << "Grid Spacing " << unitCellParameters.geta() << " " << gridYSpacing << " " << gridZSpacing << endl;

    maxX = minX + ( (totalNumGridX-1) * unitCellParameters.geta());
    maxY = minY + ( (unitCellParameters.getNb()-1) * gridYSpacing);
    maxZ = minZ + ( (unitCellParameters.getNc()-1) * gridZSpacing);
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
}*/
