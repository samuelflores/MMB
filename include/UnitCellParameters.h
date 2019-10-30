#ifndef UnitCellParameters_H_
#define UnitCellParameters_H_
#include  "SimTKmolmodel.h" // Actually this could just need simtkcommon or some such..
#include  <vector>
//typedef SimTK::Vec<3,int> iVec3;
typedef std::array<int, 3> iVec3; // This is not really a vector, it is an array. But I think it should work fine for us.
class UnitCellParameters{
	private:
		bool valid;
		double a,b,c;
		double alpha, beta, gamma;
		int  na, aMin, aMax, nb, bMin, bMax, nc, cMin, cMax;
		SimTK::Mat33 orthogonalizationMatrix;
		SimTK::Mat33 deOrthogonalizationMatrix;
	public:
		void setDefaultParameters();

		void setDeOrthogonalizationMatrix ();
                const SimTK::Mat33 getDeOrthogonalizationMatrix () const ;
		//SimTK::Vec3 multiplyMat33TimesVec3  (const SimTK::Mat33 myMat33 , const SimTK::Vec3 myVec3 ) ;
		// takes a cartesian vector, returns a vector in fractional coordinates
		SimTK::Vec3 convertCartesianVectorToFractionalVector  (const SimTK::Vec3 cartesianVector) const ;
                // takes a fractional vector, drops the integer part and returns only the fractional part:
		SimTK::Vec3 convertFractionalVectorToFractionFromLowerLeft  (const SimTK::Vec3 & fractionalVector) const ;
                // Rounds fractional vector to nearest integers
		iVec3 convertCartesianVectorToNearestIndexVector  (const SimTK::Vec3 & cartesianVector) const;
                // Truncates fractional vector components to find integer vector.
		iVec3 convertFractionalVectorToLowerIndexVector  (const SimTK::Vec3 & fractionalVector) const;
                // first calls convertCartesianVectorToFractionalVector, then calls convertFractionalVectorToLowerIndexVector:
		iVec3 convertCartesianVectorToLowerIndexVector  (const SimTK::Vec3 & cartesianVector) const;
                // Checks whether the given fractional vector is inside the unit cell.
		const bool fractionalVectorIsInsideMapBoundaries(const SimTK::Vec3 & fractionalVector);
                // Default constructor. Mainly, just sets valid to false.
		UnitCellParameters();
		const double angleInRange (const double angleInDegrees) ;
                // makes sure the n is more unity, and that nMax-nMin+1 = n.  Usable on a,b,c dimensions
		void validateNnMinnMax(const int n, const int nMin, const int nMax);
                // validates inputAngle, dies if out of reasonable range.
                const int getNa() const {return na;};
                const int getNb() const {return nb;};
                const int getNc() const {return nc;};
                const double geta() const {return a;};
                const double getb() const {return b;};
                const double getc() const {return c;};
                const double getaMin() const {return aMin;};
                const double getbMin() const {return bMin;};
                const double getcMin() const {return cMin;};
                const double getaMax() const {return aMax;};
                const double getbMax() const {return bMax;};
                const double getcMax() const {return cMax;};
                const int calcMaxFrequencyDoublingsX ();
                const int calcMaxFrequencyDoublingsY ();
                const int calcMaxFrequencyDoublingsZ ();
                const double getMinHalfWavelengthX () ;                   
                const double getMinHalfWavelengthY () ;                   
                const double getMinHalfWavelengthZ () ;                   
                const double getMaxHalfWavelengthX () {return (geta()*(getNa()-1));};// The longest half-wavelength is the width of the unit cell.
                const double getMaxHalfWavelengthY () {return (getb()*(getNb()-1));};// The longest half-wavelength is the width of the unit cell.
                const double getMaxHalfWavelengthZ () {return (getc()*(getNc()-1));};// The longest half-wavelength is the width of the unit cell.
		void setAlphaUsingDegrees(double inputAngle);
		const double getAlpha();
		void setBetaUsingDegrees(double inputAngle);
		const double getBeta();
		void setGammaUsingDegrees(double inputAngle);
		void setabc(const double mya, const double myb, const double myc);
		const double getGamma();
                // calls validateNnMinnMax three times for validation, then sets the corresponding unit cell parameters:
		void setN (const int myna, const int myaMin,const int myaMax,const int mynb,const int mybMin,const int mybMax,const int mync,const int mycMin,const int mycMax);
                // Dies if valid = 0.  This can be called liberally whenever we want to use any parameters.
		const bool exitIfNotValid() const;
                // Runs checks on the parameters:
		void validate();
		// This returns the nonorthogonal basis vector va, in orthogonal (cartesian) coordinates. 
		// fractionalMagnitude counts number of a's from origin. Can be negative. For now, we allow this to even access points outside the provided density map.
		// if fractionalMagnitude argument is unity (default), magnitude of va is a.
		// if fractionalMagnitude argument is anything else, magnitude of va is fractionalMagnitude * a
		SimTK::Vec3 va (const double fractionalMagnitude = 1);
		SimTK::Vec3 vb (const double fractionalMagnitude = 1);
		// From http://www.webmineral.com/help/CellDimensions.shtml#.XWjGX5MzYmJ 
                double totalVolume ();
		double volume ();
		SimTK::Vec3 vc (const double fractionalMagnitude = 1);
};
#endif

