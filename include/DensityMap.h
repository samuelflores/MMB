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
#include "UnitCellParameters.h"
#include <math.h>       /* exp */

#define LINESIZE 1024

using namespace std;  
using namespace SimTK;

class Quadrant {
public:
	bool positiveX, positiveY, positiveZ;
};
struct AmplitudeFrequencyAndRandomPhases {
    double amplitude;
    double frequencyX;
    double frequencyY;
    double frequencyZ;
    double phaseX;
    double phaseY;
    double phaseZ;
};

struct GridPoint   {
	double noiseFreeDensity; 
	double density; 
        double noise;/*
	float ddxPositiveXGradient;
	float ddyPositiveXGradient;
	float ddzPositiveXGradient;
        float ddxPositiveYGradient;
        float ddyPositiveYGradient;
        float ddzPositiveYGradient;
	float ddxPositiveZGradient;
	float ddyPositiveZGradient;
	float ddzPositiveZGradient;*/
	Vec3 position;
        Vec3 firstQuadrantGradient; 
};

class GridIndices {
public:
	GridIndices(int x, int y, int z) :
        z{z},
        y{y},
        x{x}
    {}

	const int z;
	const int y;
	const int x;
};


class MMB_EXPORT DensityMap {
    protected:
        //double minX, minY, minZ, maxX, maxY, maxZ, gridXSpacing, gridYSpacing, gridZSpacing;
        int unitCellNumGridX;
        int unitCellNumGridY;
        int unitCellNumGridZ;
        int totalNumGridX;
        int totalNumGridY;
        int totalNumGridZ;
        UnitCellParameters unitCellParameters;
        char * dxGets( char *s, int n, FILE *stream);
    private:
        double noiseTemperature;    		
        double noiseScale;
	double forceConstant;
	String densityFileName;

        std::vector<GridPoint> m_gridPoints;
        const GridPoint & gridPoint(int z, int y, int x) const;
        GridPoint & gridPoint(int z, int y, int x);

    public:
        DensityMap();
        void initializeMap();
        void validateGridParameters();
        size_t getSizeOfArrayOfGridPoints() const ;//
        // 	{return ArrayOfGridPoints.size()*ArrayOfGridPoints[0].size()*ArrayOfGridPoints[0][0].size();};
        std::vector<std::vector<std::vector<AmplitudeFrequencyAndRandomPhases> > > vectorOfAmplitudeFrequencyAndRandomPhases;     
        bool hasGridPoint(const GridIndices &indices) const;
        const GridPoint & getGridPoint(const GridIndices &indices) const;
        const GridPoint & getGridPoint(const Vec3 &pos) const;
        GridPoint & updGridPoint(const GridIndices &indices);
        GridPoint & updGridPoint(const Vec3 &pos);
        void validateGridPoint(const GridIndices &indices) const;
        //const bool hasNearbyGridIndices(Vec3 position);
        GridIndices calcNearestGridIndices(const Vec3 &position) const;
        GridIndices calcLowerLeftGridIndices(const Vec3 &position) const;
        //const double getDensity(Vec3);
        double getDensity(const SimTK::Vec3 &);
        void initializeArrayOfGridPoints();
        void setNoiseTemperature(double myTemperature){noiseTemperature=myTemperature;};
        void setNoiseScale(double myNoiseScale){noiseScale = myNoiseScale;};
        double getNoiseScale(){return noiseScale ;};
        //void addNoiseToMap(double temperature, double amplitude);
        void initializeVectorOfAmplitudeAndRandomPhases();
        //void resizeNoiseMap();
        void resizeVectorOfAmplitudeAndRandomPhases();
        void normalizeNoiseMap(const double totalNoiseEverywhere);
        void densityAutocorrelation(const bool computeNoiseAutocorrelation, const bool computeDensityAutocorrelation ) const;
        void populateNoiseMap();
        double getForceConstant  () const {return forceConstant; };
        void setForceConstant  (double myForceConstant  ) {forceConstant   = myForceConstant  ;};
        String getDensityFileName() const {return densityFileName;};
        void setDensityFileName(String myDensityFileName) {densityFileName = myDensityFileName;};
        void loadParametersAndDensity() ;
        void loadParametersAndDensity_XPLOR() ;
        void loadParametersAndDensity_CCP4MAP() ;
        void writeDensityMapXplor(const String outDensityFileName, const bool writeDensity = 1, const bool writeNoise =1);
        //void loadParametersAndDensity_OpenDX(const String densityFileName) ;
        //void loadParametersAndDensity_Situs(const String densityFileName) ;
        void precomputeGradient();
        //void precomputeGradientDerivatives();
        Vec3 fetchGradient(const Vec3 &position);
        const Vec3 & fetchFirstQuadrantGradient(const Vec3 &position);
        //Vec3 calcInterpolatedFirstQuadrantGradient(Vec3 position);
        const Vec3 & calcInterpolatedFirstQuadrantGradient(const Vec3 &position) const;
        // Functions which were moved from GridPoint to DensityMap for memory savings
        void initializeGradient(GridPoint & gridPoint );
        void initialize(GridPoint & gridPoint );
        void validatePosition(const GridPoint & gridPoint, const Vec3 &myPosition)const ;
        void validateDensity (const GridPoint & gridPoint, double          ) const;
        void validate(const GridPoint & gridPoint) const;
        void setDensity(GridPoint & gridPoint, Real myDensity);
        void setPosition(GridPoint & gridPoint, const Vec3 &myPosition);
        //Quadrant calcQuadrant(GridPoint & gridPoint, Vec3 queryPosition) const;
        //Vec3  fetchGradient(GridPoint & gridPoint, Vec3 queryPosition) const;
        const Vec3 & fetchFirstQuadrantGradient(const GridPoint & gridPoint) const ;
        double getDensity(const GridPoint & gridPoint) const;
        double getDensity(const GridPoint & gridPoint, const Vec3 &myPosition) const;   
        void setPositiveXGradient(GridPoint & gridPoint, Real);
        //void setddxPositiveXGradient(GridPoint & gridPoint, Real);
        //void setddyPositiveXGradient(GridPoint & gridPoint, Real);
        //void setddzPositiveXGradient(GridPoint & gridPoint, Real);
        void setPositiveYGradient(GridPoint & gridPoint,Real);
        //void setddxPositiveYGradient(GridPoint & gridPoint, Real);
        //void setddyPositiveYGradient(GridPoint & gridPoint, Real);
        //void setddzPositiveYGradient(GridPoint & gridPoint, Real);
        void setPositiveZGradient(GridPoint & gridPoint, Real);
        //void setddxPositiveZGradient(GridPoint & gridPoint, Real);
        //void setddyPositiveZGradient(GridPoint & gridPoint, Real);
        //void setddzPositiveZGradient(GridPoint & gridPoint, Real);
        void setNegativeXGradient(GridPoint & gridPoint, Real);
        void setNegativeYGradient(GridPoint & gridPoint, Real);
        void setNegativeZGradient(GridPoint & gridPoint, Real);
        //void printSecondDerivatives(GridPoint & gridPoint) const;
        const Vec3 & calcInterpolatedFirstQuadrantGradient(const GridPoint & gridPoint, const Vec3 &queryPosition) const;
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

