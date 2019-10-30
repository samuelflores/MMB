#include<cmath>
#include "UnitCellParameters.h"
#include "ErrorManager.h"
//#define RadToDeg 3.14159265/180
//#define sinDegrees(x) sin((x) * M_PI / 180.0)
void UnitCellParameters::setDefaultParameters(){
    valid = 0;
    a = 0; b = 0; c=0;
    alpha  = 0; beta  = 0; gamma  = 0;
    na = 0; aMin= 0; aMax= 0; nb= 0; bMin= 0; bMax= 0; nc= 0; cMin= 0; cMax= 0;
}

void UnitCellParameters::setDeOrthogonalizationMatrix (){ 
    //exitIfNotValid();
    /*float alpha = myUnitCellParameters.alpha;
    float beta = myUnitCellParameters.beta;
    float gamma = myUnitCellParameters.gamma;
    Float cg = cos(gamma)
    float a = myUnitCellParameters.a;
    float b = myUnitCellParameters.b;
    float c = myUnitCellParameters.c;*/
    std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<std::endl;
    //float returnMatrix[3][3] = {{0,0,0}{0,0,0}{0,0,0}};
    deOrthogonalizationMatrix[0][0] = 1 / a;
    deOrthogonalizationMatrix[1][0] = 0;
    deOrthogonalizationMatrix[2][0] = 0;
    deOrthogonalizationMatrix[0][1] = - cos(gamma)/ a / sin (gamma) ;
    deOrthogonalizationMatrix[1][1] =  1 / b / sin(gamma);
    deOrthogonalizationMatrix[2][1] = 0;
    deOrthogonalizationMatrix[0][2] = b * cos(gamma) * c*( (cos(alpha) - cos(beta)*cos(gamma))/sin(gamma)  -b*c*cos(beta)*sin(gamma)  ) / volume();
    deOrthogonalizationMatrix[1][2] =  -a*c*(cos(alpha) - cos(beta)*cos(gamma))/volume()/sin(gamma);
    deOrthogonalizationMatrix[2][2] = a*b*sin(gamma)/volume();
}
const SimTK::Mat33 UnitCellParameters::getDeOrthogonalizationMatrix () const{
    exitIfNotValid();
    return deOrthogonalizationMatrix;
}
SimTK::Vec3 multiplyMat33TimesVec3  (const SimTK::Mat33 myMat33 , const SimTK::Vec3 myVec3 ){
    SimTK::Vec3 returnVec3 = {0,0,0};
    //std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<" provided SimTK::Mat33 contains : "<<std::endl;
    //for (int i = 0; i < 3; i++) for (int j = 0 ; j < 3; j++) std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<" myMat33["<<i<<"]["<<j<<"] = "<<myMat33[i][j]<<std::endl; 
    returnVec3[0] =    myMat33[0][0] * myVec3[0] + myMat33[0][1] * myVec3[1]  + myMat33[0][2] * myVec3[2];
    returnVec3[1] =    myMat33[1][0] * myVec3[0] + myMat33[1][1] * myVec3[1]  + myMat33[1][2] * myVec3[2];
    returnVec3[2] =    myMat33[2][0] * myVec3[0] + myMat33[2][1] * myVec3[1]  + myMat33[2][2] * myVec3[2];
    return returnVec3;
}

// takes a cartesian vector, returns a vector in fractional coordinates
SimTK::Vec3 UnitCellParameters::convertCartesianVectorToFractionalVector  (const SimTK::Vec3 cartesianVector) const {
    //std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<" Converting cartesian vector "<<cartesianVector<<std::endl;
    exitIfNotValid();

    SimTK::Vec3 fractionalVector = multiplyMat33TimesVec3(getDeOrthogonalizationMatrix(), cartesianVector);
    //std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<" To fractional vector "<< fractionalVector <<std::endl;
    return fractionalVector;
}
SimTK::Vec3 UnitCellParameters::convertFractionalVectorToFractionFromLowerLeft  (const SimTK::Vec3 & fractionalVector) const{
    //std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<" Converting fractional vector "<<fractionalVector<<std::endl;
    //Vec3 fractionalVector = multiplyMat33TimesVec3(deOrthogonalizationMatrix, cartesianVector);
    SimTK::Vec3 fractionalPartOfFractionalVector = {0,0,0};
    SimTK::Vec3 roundingBuffer = {0.0000001,0.0000001,0.0000001};
    iVec3 truncatedFractionalVector = convertFractionalVectorToLowerIndexVector(fractionalVector+roundingBuffer); // The roundingBuffer is so slight negative fractions will round up to the nearest grid point. That way fractions can be on the open to closed interval (0,1] .
    truncatedFractionalVector[0] += getaMin();
    truncatedFractionalVector[1] += getbMin();
    truncatedFractionalVector[2] += getcMin();
    //std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<" To truncated (non-index) fractional vector "<<truncatedFractionalVector[0] <<", "<<  truncatedFractionalVector[1]  <<", "<<  truncatedFractionalVector[2] <<std::endl;
    for (int i = 0; i < 3; i++){
        fractionalPartOfFractionalVector[i] = fractionalVector[i] - truncatedFractionalVector[i];
        if ((fractionalPartOfFractionalVector[i] < 0.0000001) && (fractionalPartOfFractionalVector[i] > -0.0000001)) {fractionalPartOfFractionalVector[i] = 0.0;} // Get rid of tiny numbers.
    }
    //fractionalPartOfFractionalVector[0] = fractionalVector[0] - truncatedFractionalVector[0];
    //fractionalPartOfFractionalVector[1] = fractionalVector[1] - truncatedFractionalVector[1];
    //fractionalPartOfFractionalVector[2] = fractionalVector[2] - truncatedFractionalVector[2];
    //if (
    /*fractionalPartOfFractionalVector[0] = fractionalVector[0] - trunc(fractionalVector[0]);
    if (fractionalVector[0] < 0.) {
        fractionalPartOfFractionalVector[0] = 1 - fractionalPartOfFractionalVector[0];
    }
    fractionalPartOfFractionalVector[1] = fractionalVector[1] - trunc(fractionalVector[1]);
    if (fractionalVector[1] < 0.) {fractionalPartOfFractionalVector[1] = 1 - fractionalPartOfFractionalVector[1];}
    fractionalPartOfFractionalVector[2] = fractionalVector[2] - trunc(fractionalVector[2]);
    if (fractionalVector[2] < 0.) {fractionalPartOfFractionalVector[2] = 1 - fractionalPartOfFractionalVector[2];}*/
    //std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<" To fraction of unit cell from lower left: "<< fractionalPartOfFractionalVector <<std::endl;
    if ((fractionalPartOfFractionalVector[0] >= 1) || (fractionalPartOfFractionalVector[0] < 0) ||
        (fractionalPartOfFractionalVector[1] >= 1) || (fractionalPartOfFractionalVector[1] < 0) ||
        (fractionalPartOfFractionalVector[2] >= 1) || (fractionalPartOfFractionalVector[2] < 0) ) {
            ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Unexplained error! fractionalPartOfFractionalVector = " << fractionalPartOfFractionalVector <<std::endl; 
            ErrorManager::instance.treatError();
        }
    return fractionalPartOfFractionalVector;
}
iVec3 UnitCellParameters::convertCartesianVectorToNearestIndexVector  (const SimTK::Vec3 & cartesianVector) const {
    //std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<" Converting cartesian vector "<<cartesianVector<<std::endl;
    SimTK::Vec3 fractionalVector = convertCartesianVectorToFractionalVector( cartesianVector);
    iVec3 indexVector = {0,0,0};
    indexVector[0] = round(fractionalVector[0]) - aMin;
    indexVector[1] = round(fractionalVector[1]) - bMin;
    indexVector[2] = round(fractionalVector[2]) - cMin;
    //std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<" To rounded index vector "<< indexVector[0] <<","<< indexVector[1]<<","<<  indexVector[2] <<std::endl;
    return indexVector;}
iVec3 UnitCellParameters::convertFractionalVectorToLowerIndexVector  (const SimTK::Vec3 & fractionalVector) const {
    #ifdef _DEBUG_FLAGS_ON_
    //std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<" Converting fractional coordinate vector "<<fractionalVector<<std::endl;
    #endif
    //Vec3 fractionalVector = convertCartesianVectorToFractionalVector( fractionalVector);
    iVec3 indexVector = {0,0,0};
    indexVector[0] = trunc(fractionalVector[0]) - aMin;
    if ((fractionalVector[0] - trunc(fractionalVector[0])) < 0 ) {indexVector[0] -= 1;} // If this component of the fractional vector is negative, the lower left index is computed a bit differently
    indexVector[1] = trunc(fractionalVector[1]) - bMin;
    if ((fractionalVector[1] - trunc(fractionalVector[1])) < 0 ) {indexVector[1] -= 1;}
    indexVector[2] = trunc(fractionalVector[2]) - cMin;
    if ((fractionalVector[2] - trunc(fractionalVector[2])) < 0 ) {indexVector[2] -= 1;}
    //std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<" To truncated (lower left) index vector "<< indexVector[0] <<","<< indexVector[1]<<","<<  indexVector[2]  <<std::endl;
    return indexVector;
}
iVec3 UnitCellParameters::convertCartesianVectorToLowerIndexVector  (const SimTK::Vec3 & cartesianVector) const {
    //std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<" Converting cartesian vector "<<cartesianVector<<std::endl;
    SimTK::Vec3 fractionalVector = convertCartesianVectorToFractionalVector( cartesianVector);
    //std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<" .. to fractional coordinate vector "<<fractionalVector<<std::endl;
    iVec3 indexVector = {0,0,0};
    indexVector = convertFractionalVectorToLowerIndexVector (fractionalVector);
    //indexVector[0] = trunc(fractionalVector[0]) - aMin;
    //indexVector[1] = trunc(fractionalVector[1]) - bMin;
    //indexVector[2] = trunc(fractionalVector[2]) - cMin;
    //std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<" .. and finally to truncated index vector "<< indexVector[0] <<","<< indexVector[1]<<","<<  indexVector[2]  <<std::endl;
    return indexVector;
}
const bool UnitCellParameters::fractionalVectorIsInsideMapBoundaries(const SimTK::Vec3 & fractionalVector){
    //std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<std::endl;
    exitIfNotValid();
    //std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<std::endl;
    if ((fractionalVector[0] < aMin) || (fractionalVector[0] > aMax) ||
        (fractionalVector[1] < bMin) || (fractionalVector[1] > bMax) ||
        (fractionalVector[2] < cMin) || (fractionalVector[2] > cMax) ) {
        //std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<std::endl;
        return 0;} 
    else {
        //std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<std::endl;
        return 1;}
}

        
// default constructor. mainly here to set valid to false.
UnitCellParameters::UnitCellParameters(){  setDefaultParameters();valid = 0;}
        
const double UnitCellParameters::angleInRange  (const double angleInDegrees){
    const double minAngle = 4.0*SimTK::Pi/180; // we switched to rads // by setting the minimum angle to a little over pi, I intend to prevent people from providing input in rads. 4 degrees is quite acute anyway.
    const double maxAngle = 170.0*SimTK::Pi/180; // we switched to rads // Technically this could be 180 degrees, but why would anyone use such an obtuse angle?
    if (angleInDegrees >= maxAngle) {
        valid = 0;
        ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" You have provided an angle of : "<<angleInDegrees<< " degrees, which matches or exceeds our maximum of "<< maxAngle <<" degrees. ";
        ErrorManager::instance.treatError();
    }
    else if (angleInDegrees <= minAngle) {
        valid = 0;
        ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" You have provided an angle of : "<<angleInDegrees<< " degrees, which matches or is under our minimum of "<< minAngle <<" degrees. ";
        ErrorManager::instance.treatError();
    } else {
        std::cout << __FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" You have provided an angle of : "<<angleInDegrees<< " degrees, which is between our allowed extremes of " << minAngle << " and " << maxAngle << " degrees. Validation passed."<<std::endl;
        return angleInDegrees;
    }
    // no idea how we would ever get to this point, so trip error.
    ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Unexplained error! ";
    ErrorManager::instance.treatError();
 }

    
void UnitCellParameters::validateNnMinnMax(const int n, const int nMin, const int nMax){
        if ( n < 2 ){
            valid = 0;
            ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" You have provided a number of grid points : "<< n << " which is under our minimum of 2.";
            ErrorManager::instance.treatError();
        }
        if ( n > 10000 ){
            valid = 0;
            ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" You have provided a number of grid points : "<< n << " which exceeds our maximum.";
            ErrorManager::instance.treatError();
        }
        if ( n != (nMax - nMin +1 ) ){
            valid = 0;
            ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" You have provided a number of grid points : "<< n << ". Your grid starts at point " << nMin << " and ends at point " << nMax << " These numbers are not consistent.";
            ErrorManager::instance.treatError();
        }
        // If we get this far, everything is kosher. Do nothing.

    }
const int calcMaxFrequencyDoublings (const int numGridPoints) {  
    /* 
    double maxFrequencyDoublings = log2 (numGridPoints -1);
    //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" The Nyquist half-wavelength, set by a grid spacing of "<<gridSpacing<<" is 2^"<<maxFrequencyDoublings<<" times smaller than than the max half-wavelength, set by the unit cell width of "<<unitCellWidth;//<<std::endl;
    maxFrequencyDoublings = floor (maxFrequencyDoublings);
    //std::cout<<" we round down to "<<maxFrequencyDoublings<< " frequency doublings."<<std::endl;
    return maxFrequencyDoublings; // The rounded-down number should be cast as int            
    */
    return (numGridPoints-1);
}
const int calcMaxFrequencyDoublings (const double unitCellWidth,     const double gridSpacing) {
    double maxFrequencyDoublings = log2 (unitCellWidth / gridSpacing);
    //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" The Nyquist half-wavelength, set by a grid spacing of "<<gridSpacing<<" is 2^"<<maxFrequencyDoublings<<" times smaller than than the max half-wavelength, set by the unit cell width of "<<unitCellWidth;//<<std::endl;
    maxFrequencyDoublings = floor (maxFrequencyDoublings);
    //std::cout<<" we round down to "<<maxFrequencyDoublings<< " frequency doublings."<<std::endl;
    return maxFrequencyDoublings; // The rounded-down number should be cast as int            
}

const int UnitCellParameters::calcMaxFrequencyDoublingsX (){
    //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
    return calcMaxFrequencyDoublings((getNa()));
}
const int UnitCellParameters::calcMaxFrequencyDoublingsY (){
    //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
    //return calcMaxFrequencyDoublings(((getNb()-1)*getb()), getb());
    return calcMaxFrequencyDoublings((getNb()));
}
const int UnitCellParameters::calcMaxFrequencyDoublingsZ (){
    //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
    //return calcMaxFrequencyDoublings(((getNc()-1)*getc()), getc());
    return calcMaxFrequencyDoublings((getNc()));
}

void UnitCellParameters::setAlphaUsingDegrees(double inputAngle){
    alpha = angleInRange (inputAngle*M_PI/180);
}
const double UnitCellParameters::getAlpha(){ return alpha;}
void UnitCellParameters::setBetaUsingDegrees(double inputAngle){
    beta = angleInRange (inputAngle*M_PI/180);
}
const double UnitCellParameters::getBeta(){ return beta;}
void UnitCellParameters::setGammaUsingDegrees(double inputAngle){
    gamma = angleInRange (inputAngle*M_PI/180);
}
void UnitCellParameters::setabc(const double mya, const double myb, const double  myc){
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" About to set a, b, c to : "<<mya <<", "<<  myb <<", "<< myc<<std::endl;
    if (mya <= 0.0){
        ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" You have provided an invalid value for a   : "<< mya <<std::endl;
        ErrorManager::instance.treatError();
    }else {a = mya;}
    if (myb <= 0.0){
        ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" You have provided an  invalid value for b : "<< myb <<std::endl;
        ErrorManager::instance.treatError();
    }else {b = myb;}
    if (myc <= 0.0){
        ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" You have provided an  invalid value for c : "<< myc <<std::endl;
        ErrorManager::instance.treatError();
    }else {c = myc;}
}
const double UnitCellParameters::getGamma(){ return gamma;}

void UnitCellParameters::setN ( const int myna, const int myaMin,const int myaMax,const int mynb,const int mybMin,const int mybMax,const int mync,const int mycMin,const int mycMax ){
   std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" About to set na, aMin, aMax to : "<<myna <<", "<<  myaMin <<", "<< myaMax<<" nm "<<std::endl;
   validateNnMinnMax(myna,  myaMin, myaMax); // validateNnMinnMax will kill the program if any of its arguments are invalid.
   na = myna; aMin = myaMin; aMax = myaMax;
   std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Confirming na, aMin, aMax are now : "<<getNa()<<", "<<getaMin()<<", "<<getaMax()<<std::endl;
   std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" About to set nb, bMin, bMax to : "<<mynb <<", "<<  mybMin <<", "<< mybMax<<std::endl;
   validateNnMinnMax(mynb,  mybMin, mybMax);
   nb = mynb; bMin = mybMin; bMax = mybMax;
   std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Confirming nb, bMin, bMax are now : "<<getNb()<<", "<<getbMin()<<", "<<getbMax()<<std::endl;
   std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" About to set nc, cMin, cMax to : "<<mync <<", "<<  mycMin <<", "<< mycMax<<std::endl;
   validateNnMinnMax(mync,  mycMin, mycMax);
   nc = mync; cMin = mycMin; cMax = mycMax;
   std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Confirming nc, cMin, cMax are now : "<<getNc()<<", "<<getcMin()<<", "<<getcMax()<<std::endl;
}
const bool UnitCellParameters::exitIfNotValid() const{
    if (valid) { 
        return true; // All is fine.  
    } else {
    ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" The variable named \'valid\' is set to " << valid << ". This means your parameters have not been checked. Exiting now. ";
    ErrorManager::instance.treatError();
    }
}
 void UnitCellParameters::validate (){
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" About to check your parameters. Currently the variable named \'valid\' is set to " << valid <<std::endl;
    valid = 0;
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" About to check your parameters. The variable named \'valid\' is now deliberately set to " << valid <<std::endl;
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Confirming na, aMin, aMax are now : "<<getNa()<<", "<<getaMin()<<", "<<getaMax()<<std::endl;
    validateNnMinnMax(getNa(),getaMin(),getaMax());
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Confirming nb, bMin, bMax are now : "<<getNb()<<", "<<getbMin()<<", "<<getbMax()<<std::endl;
    validateNnMinnMax(getNb(),getbMin(),getbMax());
    validateNnMinnMax(getNc(),getcMin(),getcMax());
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" getAlpha() = "<<getAlpha()<<std::endl;
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" angleInRange(getAlpha() = "<<angleInRange(getAlpha())<<std::endl;
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" getBeta() = "<<getBeta()<<std::endl;
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" getGamma() = "<<getGamma()<<std::endl;
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;

    angleInRange(getAlpha());
    if (angleInRange(getAlpha()) &&
        angleInRange(getBeta()) &&
        angleInRange(getGamma()) ) 
    {
        //valid = true;
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Your parameters have been checked and found to be kosher. Pending further checks, we still set \'valid\' to "<<valid<<std::endl;
    } 
    else {
        valid = false;
        ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" The variable named \'valid\' is set to " << valid << ". This means your parameters have not been checked. Exiting now. ";
        ErrorManager::instance.treatError();
    }
    if (deOrthogonalizationMatrix[0][0] != 1 / geta() ) {valid = 0; ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Found an error in deOrthogonalizationMatrix " ; ErrorManager::instance.treatError();}
    if (deOrthogonalizationMatrix[1][0] != 0 ) {valid = 0; ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Found an error in deOrthogonalizationMatrix " ; ErrorManager::instance.treatError();}
    if (deOrthogonalizationMatrix[2][0] != 0) {valid = 0; ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Found an error in deOrthogonalizationMatrix " ; ErrorManager::instance.treatError();}
    if (deOrthogonalizationMatrix[0][1] != - cos(gamma)/ a / sin (gamma) ) {valid = 0; ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Found an error in deOrthogonalizationMatrix " ; ErrorManager::instance.treatError();}
    if (deOrthogonalizationMatrix[1][1] !=  1 / b / sin(gamma)) {valid = 0; ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Found an error in deOrthogonalizationMatrix " ; ErrorManager::instance.treatError();}
    if (deOrthogonalizationMatrix[2][1] != 0) {valid = 0; ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Found an error in deOrthogonalizationMatrix " ; ErrorManager::instance.treatError();}
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"deOrthogonalizationMatrix[0][2] = "<< deOrthogonalizationMatrix[0][2]<<std::endl;
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"volume() = "<< volume()<<std::endl;
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" b * cos(gamma) * c*( (cos(alpha) - cos(beta)*cos(gamma))/sin(gamma)  -b*c*cos(beta)*sin(gamma)  ) / volume() = " <<b * cos(gamma) * c*( (cos(alpha) - cos(beta)*cos(gamma))/sin(gamma)  -b*c*cos(beta)*sin(gamma)  ) / volume() <<std::endl;
    if (deOrthogonalizationMatrix[0][2] != b * cos(gamma) * c*( (cos(alpha) - cos(beta)*cos(gamma))/sin(gamma)  -b*c*cos(beta)*sin(gamma)  ) / volume()) {valid = 0; ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Found an error in deOrthogonalizationMatrix " ; ErrorManager::instance.treatError();}
    if (deOrthogonalizationMatrix[1][2] !=  -a*c*(cos(alpha) - cos(beta)*cos(gamma))/volume()/sin(gamma)) {valid = 0; ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Found an error in deOrthogonalizationMatrix " ; ErrorManager::instance.treatError();}
    if (deOrthogonalizationMatrix[2][2] != a*b*sin(gamma)/volume()) {valid = 0; ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Found an error in deOrthogonalizationMatrix " ; ErrorManager::instance.treatError();}
    valid = 1; // If we got this far, all is kosher
}

// This returns the nonorthogonal basis vector va, in orthogonal (cartesian) coordinates. 
// fractionalMagnitude counts number of a's from origin. Can be negative. For now, we allow this to even access points outside the provided density map.
// if fractionalMagnitude argument is unity (default), magnitude of va is a.
// if fractionalMagnitude argument is anything else, magnitude of va is fractionalMagnitude * a
SimTK::Vec3 UnitCellParameters::va  (const double fractionalMagnitude ){
    SimTK::Vec3 returnVec = {0,0,0};
    returnVec[0] = a * fractionalMagnitude;
    returnVec[1] = 0;
    returnVec[2] = 0;
    return returnVec;
}

SimTK::Vec3 UnitCellParameters::vb  (const double fractionalMagnitude ){ 
    SimTK::Vec3 returnVec = {0,0,0};
    returnVec[0] = b * fractionalMagnitude * cos(gamma);
    returnVec[0] = b * fractionalMagnitude * sin(gamma);
    returnVec[0] = 0;
    return returnVec;
}

double UnitCellParameters::totalVolume (){ // volume() is for one voxel. This is for the whole unit cell
    return (volume() * (getNa() - 1) * (getNb() - 1) * (getNc() - 1));    
}

// From http://www.webmineral.com/help/CellDimensions.shtml#.XWjGX5MzYmJ 
double UnitCellParameters::volume ()
{   
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"  ((cos(alpha) * cos(beta) * cos(gamma)) ) = "<<  ((cos(alpha) * cos(beta) * cos(gamma)) )<<std::endl;
    double arg1 = cos(alpha) * cos(beta) * cos(gamma);
    if ((arg1 < 0.0000001) && (arg1  > -0.0000001)){arg1 = 0.0;} // had trouble with extremely tiny negative values of arg1
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" pow ((cos(alpha) * cos(beta) * cos(gamma)) ,0.5) = "<< pow ((cos(alpha) * cos(beta) * cos(gamma)) ,0.5)<<std::endl;
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" a*b*c = "<<a*b*c<<std::endl; 
    double returnVolume = a*b*c* (1- pow(cos (alpha),2) - pow(cos (beta) , 2) - pow( cos (gamma) , 2)) + 2 * pow ((arg1) ,0.5);
    return returnVolume; 
}


SimTK::Vec3 UnitCellParameters::vc  (const double fractionalMagnitude){ 
    SimTK::Vec3 returnVec = {0,0,0};
    returnVec[0] = c * fractionalMagnitude * cos(beta);
    returnVec[1] = c * fractionalMagnitude * (cos(alpha) - cos(beta) * cos(gamma)  ) / sin ( gamma );
    returnVec[2] = fractionalMagnitude * volume () /  (a*b * sin (gamma)) ;
    return returnVec;
}


