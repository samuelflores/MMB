/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#ifndef DensityMap_H_
#define DensityMap_H_

#include "SimTKmolmodel.h"

#define LINESIZE 1024

using namespace std;  
using namespace SimTK;

class Quadrant {
public:
	bool positiveX, positiveY, positiveZ;
};

struct GridPoint   {
//private:
//public:
	float density; 
	float ddyPositiveXGradient;
	float ddzPositiveXGradient;
        float ddxPositiveYGradient;
        float ddzPositiveYGradient;
	float ddxPositiveZGradient;
	float ddyPositiveZGradient;
	fVec3 position;
        fVec3 firstQuadrantGradient; 
        //float dummy1;
        //float dummy1;
//public:  // obsolete, moved to DensityMap to save memory
        /*
        GridPoint();
	void initializeGradient();
	void initialize();
	void validatePosition(fVec3 myPosition)const ;
	void validateDensity (double          ) const;
	void validate() const;
	void setDensity(Real myDensity);
	void setPosition(fVec3 myPosition);
	Quadrant calcQuadrant(fVec3 queryPosition) const;
	fVec3  fetchGradient(fVec3 queryPosition) const;
        fVec3 fetchFirstQuadrantGradient() const ;
	double getDensity() const;
	double getDensity(fVec3 myPosition) const;   
	void setPositiveXGradient(Real);
	void setddyPositiveXGradient(Real);
	void setddzPositiveXGradient(Real);
	void setPositiveYGradient(Real);
	void setddxPositiveYGradient(Real);
	void setddzPositiveYGradient(Real);
	void setPositiveZGradient(Real);
	void setddxPositiveZGradient(Real);
	void setddyPositiveZGradient(Real);
	void setNegativeXGradient(Real);
	void setNegativeYGradient(Real);
	void setNegativeZGradient(Real);
        fVec3 firstQuadrantGradient; // might be obsolete
	fVec3 calcInterpolatedFirstQuadrantGradient(fVec3 queryPosition) const;
	fVec3 getFirstQuadrantGradient();
        */	
};

class GridIndices {
private:
	int  xGridPoint; 
	int  yGridPoint;
	int  zGridPoint;
public:
	GridIndices( int myXIndex,  int myYIndex,  int myZIndex) ;
	int getXGridIndex () const;
	int getYGridIndex () const;
	int getZGridIndex () const;
};


class MMB_EXPORT DensityMap {
	protected:
	 	double minX, minY, minZ, maxX, maxY, maxZ, gridXSpacing, gridYSpacing, gridZSpacing;
		int unitCellNumGridX;
		int unitCellNumGridY;
		int unitCellNumGridZ;
		int totalNumGridX;
		int totalNumGridY;
		int totalNumGridZ;

        char * dxGets( char *s, int n, FILE *stream);
		
	public:
        DensityMap();
        ~DensityMap();
        void initializeMap();
        void validateGridParameters();
        std::vector<std::vector<std::vector<GridPoint> > > ArrayOfGridPoints;        
		bool hasGridPoint(GridIndices);
		GridPoint     & updGridPoint(GridIndices);
		const GridPoint getGridPoint(GridIndices) ;
		void validateGridPoint(GridIndices myGridIndices);
		const bool hasNearbyGridIndices(fVec3 position);
		GridIndices calcNearestGridIndices(fVec3 position);
                GridIndices calcLowerLeftGridIndices(fVec3 position);
		const GridPoint getGridPoint(fVec3);
		GridPoint     & updGridPoint(fVec3);
		const double getDensity(fVec3);
		const double getDensity(SimTK::Vec3);
		void initializeArrayOfGridPoints();
        void loadParametersAndDensity(const String densityFileName) ;
        void loadParametersAndDensity_XPLOR(const String densityFileName) ;
        void loadParametersAndDensity_OpenDX(const String densityFileName) ;
		void loadParametersAndDensity_Situs(const String densityFileName) ;
		void precomputeGradient();
		void precomputeGradientDerivatives();
		fVec3 fetchGradient(fVec3 position);
                fVec3 fetchFirstQuadrantGradient(fVec3 position);
                fVec3 calcInterpolatedFirstQuadrantGradient(fVec3 position);
                SimTK::Vec3 calcInterpolatedFirstQuadrantGradient(SimTK::Vec3 position) ;
                // Functions which were moved from GridPoint to DensityMap for memory savings
		void initializeGradient(GridPoint & gridPoint );
		void initialize(GridPoint & gridPoint );
		void validatePosition(GridPoint & gridPoint, fVec3 myPosition)const ;
		void validateDensity (GridPoint & gridPoint, double          ) const;
		void validate(GridPoint & gridPoint) const;
		void setDensity(GridPoint & gridPoint, Real myDensity);
		void setPosition(GridPoint & gridPoint, fVec3 myPosition);
		Quadrant calcQuadrant(GridPoint & gridPoint, fVec3 queryPosition) const;
		fVec3  fetchGradient(GridPoint & gridPoint, fVec3 queryPosition) const;
		fVec3 fetchFirstQuadrantGradient(GridPoint & gridPoint) const ;
		double getDensity(GridPoint & gridPoint) const;
		double getDensity(GridPoint & gridPoint, fVec3 myPosition) const;   
		void setPositiveXGradient(GridPoint & gridPoint, Real);
		void setddyPositiveXGradient(GridPoint & gridPoint, Real);
		void setddzPositiveXGradient(GridPoint & gridPoint, Real);
		void setPositiveYGradient(GridPoint & gridPoint,Real);
		void setddxPositiveYGradient(GridPoint & gridPoint, Real);
		void setddzPositiveYGradient(GridPoint & gridPoint, Real);
		void setPositiveZGradient(GridPoint & gridPoint, Real);
		void setddxPositiveZGradient(GridPoint & gridPoint, Real);
		void setddyPositiveZGradient(GridPoint & gridPoint, Real);
		void setNegativeXGradient(GridPoint & gridPoint, Real);
		void setNegativeYGradient(GridPoint & gridPoint, Real);
		void setNegativeZGradient(GridPoint & gridPoint, Real);
		void setFirstQuadrantGradient(GridPoint & gridPoint,fVec3 gradient); // might be obsolete
		//fVec3 getFirstQuadrantGradient(GridPoint & gridPoint); 
		//fVec3 firstQuadrantGradient(GridPoint & gridPoint); // might be obsolete
		fVec3 calcInterpolatedFirstQuadrantGradient(GridPoint & gridPoint, fVec3 queryPosition) const;
		//fVec3 getFirstQuadrantGradient(GridPoint & gridPoint);

                
};

// #define LINESIZE 1024
// class MMB_EXPORT DensityMap_OpenDX : DensityMap
// {
// protected:
//     char * dxGets( char *s, int n, FILE *stream);
// public:
//     DensityMap_OpenDX();
//     ~DensityMap_OpenDX();
//     void loadParametersAndDensity(const String densityFileName) ;
// };



#endif

