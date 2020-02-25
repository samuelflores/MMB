
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "../../../wrappers/OpenMMCWrapper.h"
#include "AmoebaOpenMMCWrapper.h"
#include <cstring>
#include <vector>

using namespace OpenMM;
using namespace std;

extern "C" {

/* OpenMM_2D_IntArray */
OPENMM_EXPORT_AMOEBA OpenMM_2D_IntArray* OpenMM_2D_IntArray_create(int size) {
    return reinterpret_cast<OpenMM_2D_IntArray*>(new vector<vector<int> >(size));
}
OPENMM_EXPORT_AMOEBA void OpenMM_2D_IntArray_destroy(OpenMM_2D_IntArray* array) {
    delete reinterpret_cast<vector<vector<int> >*>(array);
}
OPENMM_EXPORT_AMOEBA int OpenMM_2D_IntArray_getSize(const OpenMM_2D_IntArray* array) {
    return reinterpret_cast<const vector<vector<int> >*>(array)->size();
}
OPENMM_EXPORT_AMOEBA void OpenMM_2D_IntArray_resize(OpenMM_2D_IntArray* array, int size) {
    reinterpret_cast<vector<vector<int> >*>(array)->resize(size);
}
OPENMM_EXPORT_AMOEBA void OpenMM_2D_IntArray_append(OpenMM_2D_IntArray* array, int index1, int value) {
    vector<vector<int> >* array2DInt = reinterpret_cast<vector<vector<int> >*>(array);
    if (array2DInt->size() <= index1) {
        array2DInt->resize(index1+1);
    }
    (*array2DInt)[index1].push_back(value);
}
OPENMM_EXPORT_AMOEBA void OpenMM_2D_IntArray_set(OpenMM_2D_IntArray* array, int index1, int index2, int value) {
    vector<vector<int> >* array2DInt = reinterpret_cast<vector<vector<int> >*>(array);
    if (array2DInt->size() <= index1) {
        array2DInt->resize(index1+1);
    }
    if (array2DInt[index1].size() <= index2) {
        array2DInt[index1].resize(index2+1);
    }
    (*array2DInt)[index1][index2] = value;
}
OPENMM_EXPORT_AMOEBA void OpenMM_2D_IntArray_get(const OpenMM_2D_IntArray* array, int index1, int index2, int* value) {
    const vector<vector<int> >* array2DInt = reinterpret_cast<const vector<vector<int> >*>(array);
    if (array2DInt->size() <= index1)
        throw OpenMMException("OpenMM_2D_IntArray_get: first index out of range.");

    if ((*array2DInt)[index1].size() <= index2)
        throw OpenMMException("OpenMM_2D_IntArray_get: second index out of range.");
    *value = (*array2DInt)[index1][index2];
}

/* OpenMM_3D_DoubleArray */
OPENMM_EXPORT_AMOEBA OpenMM_3D_DoubleArray* OpenMM_3D_DoubleArray_create(int size1, int size2, int size3) {
    int ii, jj;  
    std::vector< std::vector< std::vector<double> > >* v3D_Array = new std::vector<std::vector<std::vector<double> > >(size1);

    for (ii = 0; ii < size1; ii++) {
        (*v3D_Array)[ii].resize(size2);
        for (jj = 0; jj < size2; jj++) {
           (*v3D_Array)[ii][jj].resize(size3);
        }    
    }    
    return reinterpret_cast<OpenMM_3D_DoubleArray*>(v3D_Array);
}

OPENMM_EXPORT_AMOEBA void OpenMM_3D_DoubleArray_set(OpenMM_3D_DoubleArray* array, int index1, int index2, OpenMM_DoubleArray* values) {
    unsigned int ii;
    std::vector< std::vector< std::vector<double> > >* v3D_Array = reinterpret_cast<std::vector<std::vector<std::vector<double> > >*>(array);
    std::vector<double> * value_array                            = reinterpret_cast<std::vector<double> *>(values);
    for (ii = 0; ii < (*value_array).size(); ii++) {
        (*v3D_Array)[index1][index2][ii] = (*value_array)[ii];
    }    
}

OPENMM_EXPORT_AMOEBA void OpenMM_3D_DoubleArray_destroy(OpenMM_3D_DoubleArray* array) {
    delete reinterpret_cast<std::vector<std::vector<std::vector<double> > >*>(array);
}

/* OpenMM::AmoebaMultipoleForce */
OPENMM_EXPORT_AMOEBA OpenMM_AmoebaMultipoleForce* OpenMM_AmoebaMultipoleForce_create() {
    return reinterpret_cast<OpenMM_AmoebaMultipoleForce*>(new OpenMM::AmoebaMultipoleForce());
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_destroy(OpenMM_AmoebaMultipoleForce* target) {
    delete reinterpret_cast<OpenMM::AmoebaMultipoleForce*>(target);
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaMultipoleForce_getNumMultipoles(const OpenMM_AmoebaMultipoleForce* target) {
    int result = reinterpret_cast<const OpenMM::AmoebaMultipoleForce*>(target)->getNumMultipoles();
    return result;
}
OPENMM_EXPORT_AMOEBA OpenMM_AmoebaMultipoleForce_NonbondedMethod OpenMM_AmoebaMultipoleForce_getNonbondedMethod(const OpenMM_AmoebaMultipoleForce* target) {
    OpenMM::AmoebaMultipoleForce::NonbondedMethod result = reinterpret_cast<const OpenMM::AmoebaMultipoleForce*>(target)->getNonbondedMethod();
    return static_cast<OpenMM_AmoebaMultipoleForce_NonbondedMethod>(result);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_setNonbondedMethod(OpenMM_AmoebaMultipoleForce* target, OpenMM_AmoebaMultipoleForce_NonbondedMethod method) {
    reinterpret_cast<OpenMM::AmoebaMultipoleForce*>(target)->setNonbondedMethod(static_cast<OpenMM::AmoebaMultipoleForce::NonbondedMethod>(method));
}
OPENMM_EXPORT_AMOEBA OpenMM_AmoebaMultipoleForce_PolarizationType OpenMM_AmoebaMultipoleForce_getPolarizationType(const OpenMM_AmoebaMultipoleForce* target) {
    OpenMM::AmoebaMultipoleForce::PolarizationType result = reinterpret_cast<const OpenMM::AmoebaMultipoleForce*>(target)->getPolarizationType();
    return static_cast<OpenMM_AmoebaMultipoleForce_PolarizationType>(result);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_setPolarizationType(OpenMM_AmoebaMultipoleForce* target, OpenMM_AmoebaMultipoleForce_PolarizationType type) {
    reinterpret_cast<OpenMM::AmoebaMultipoleForce*>(target)->setPolarizationType(static_cast<OpenMM::AmoebaMultipoleForce::PolarizationType>(type));
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaMultipoleForce_getCutoffDistance(const OpenMM_AmoebaMultipoleForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaMultipoleForce*>(target)->getCutoffDistance();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_setCutoffDistance(OpenMM_AmoebaMultipoleForce* target, double distance) {
    reinterpret_cast<OpenMM::AmoebaMultipoleForce*>(target)->setCutoffDistance(distance);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_getPMEParameters(const OpenMM_AmoebaMultipoleForce* target, double* alpha, int* nx, int* ny, int* nz) {
    reinterpret_cast<const OpenMM::AmoebaMultipoleForce*>(target)->getPMEParameters(*reinterpret_cast<double*>(alpha), *reinterpret_cast<int*>(nx), *reinterpret_cast<int*>(ny), *reinterpret_cast<int*>(nz));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_setPMEParameters(OpenMM_AmoebaMultipoleForce* target, double alpha, int nx, int ny, int nz) {
    reinterpret_cast<OpenMM::AmoebaMultipoleForce*>(target)->setPMEParameters(alpha, nx, ny, nz);
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaMultipoleForce_getAEwald(const OpenMM_AmoebaMultipoleForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaMultipoleForce*>(target)->getAEwald();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_setAEwald(OpenMM_AmoebaMultipoleForce* target, double aewald) {
    reinterpret_cast<OpenMM::AmoebaMultipoleForce*>(target)->setAEwald(aewald);
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaMultipoleForce_getPmeBSplineOrder(const OpenMM_AmoebaMultipoleForce* target) {
    int result = reinterpret_cast<const OpenMM::AmoebaMultipoleForce*>(target)->getPmeBSplineOrder();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_getPmeGridDimensions(const OpenMM_AmoebaMultipoleForce* target, OpenMM_IntArray* gridDimension) {
    reinterpret_cast<const OpenMM::AmoebaMultipoleForce*>(target)->getPmeGridDimensions(*reinterpret_cast<std::vector< int >*>(gridDimension));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_setPmeGridDimensions(OpenMM_AmoebaMultipoleForce* target, const OpenMM_IntArray* gridDimension) {
    reinterpret_cast<OpenMM::AmoebaMultipoleForce*>(target)->setPmeGridDimensions(*reinterpret_cast<const std::vector< int >*>(gridDimension));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_getPMEParametersInContext(const OpenMM_AmoebaMultipoleForce* target, const OpenMM_Context* context, double* alpha, int* nx, int* ny, int* nz) {
    reinterpret_cast<const OpenMM::AmoebaMultipoleForce*>(target)->getPMEParametersInContext(*reinterpret_cast<const Context*>(context), *reinterpret_cast<double*>(alpha), *reinterpret_cast<int*>(nx), *reinterpret_cast<int*>(ny), *reinterpret_cast<int*>(nz));
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaMultipoleForce_addMultipole(OpenMM_AmoebaMultipoleForce* target, double charge, const OpenMM_DoubleArray* molecularDipole, const OpenMM_DoubleArray* molecularQuadrupole, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY, double thole, double dampingFactor, double polarity) {
    int result = reinterpret_cast<OpenMM::AmoebaMultipoleForce*>(target)->addMultipole(charge, *reinterpret_cast<const std::vector< double >*>(molecularDipole), *reinterpret_cast<const std::vector< double >*>(molecularQuadrupole), axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY, thole, dampingFactor, polarity);
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_getMultipoleParameters(const OpenMM_AmoebaMultipoleForce* target, int index, double* charge, OpenMM_DoubleArray* molecularDipole, OpenMM_DoubleArray* molecularQuadrupole, int* axisType, int* multipoleAtomZ, int* multipoleAtomX, int* multipoleAtomY, double* thole, double* dampingFactor, double* polarity) {
    reinterpret_cast<const OpenMM::AmoebaMultipoleForce*>(target)->getMultipoleParameters(index, *reinterpret_cast<double*>(charge), *reinterpret_cast<std::vector< double >*>(molecularDipole), *reinterpret_cast<std::vector< double >*>(molecularQuadrupole), *reinterpret_cast<int*>(axisType), *reinterpret_cast<int*>(multipoleAtomZ), *reinterpret_cast<int*>(multipoleAtomX), *reinterpret_cast<int*>(multipoleAtomY), *reinterpret_cast<double*>(thole), *reinterpret_cast<double*>(dampingFactor), *reinterpret_cast<double*>(polarity));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_setMultipoleParameters(OpenMM_AmoebaMultipoleForce* target, int index, double charge, const OpenMM_DoubleArray* molecularDipole, const OpenMM_DoubleArray* molecularQuadrupole, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY, double thole, double dampingFactor, double polarity) {
    reinterpret_cast<OpenMM::AmoebaMultipoleForce*>(target)->setMultipoleParameters(index, charge, *reinterpret_cast<const std::vector< double >*>(molecularDipole), *reinterpret_cast<const std::vector< double >*>(molecularQuadrupole), axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY, thole, dampingFactor, polarity);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_setCovalentMap(OpenMM_AmoebaMultipoleForce* target, int index, OpenMM_AmoebaMultipoleForce_CovalentType typeId, const OpenMM_IntArray* covalentAtoms) {
    reinterpret_cast<OpenMM::AmoebaMultipoleForce*>(target)->setCovalentMap(index, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(typeId), *reinterpret_cast<const std::vector< int >*>(covalentAtoms));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_getCovalentMap(const OpenMM_AmoebaMultipoleForce* target, int index, OpenMM_AmoebaMultipoleForce_CovalentType typeId, OpenMM_IntArray* covalentAtoms) {
    reinterpret_cast<const OpenMM::AmoebaMultipoleForce*>(target)->getCovalentMap(index, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(typeId), *reinterpret_cast<std::vector< int >*>(covalentAtoms));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_getCovalentMaps(const OpenMM_AmoebaMultipoleForce* target, int index, OpenMM_2D_IntArray* covalentLists) {
    reinterpret_cast<const OpenMM::AmoebaMultipoleForce*>(target)->getCovalentMaps(index, *reinterpret_cast<std::vector< std::vector< int > >*>(covalentLists));
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaMultipoleForce_getMutualInducedMaxIterations(const OpenMM_AmoebaMultipoleForce* target) {
    int result = reinterpret_cast<const OpenMM::AmoebaMultipoleForce*>(target)->getMutualInducedMaxIterations();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations(OpenMM_AmoebaMultipoleForce* target, int inputMutualInducedMaxIterations) {
    reinterpret_cast<OpenMM::AmoebaMultipoleForce*>(target)->setMutualInducedMaxIterations(inputMutualInducedMaxIterations);
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaMultipoleForce_getMutualInducedTargetEpsilon(const OpenMM_AmoebaMultipoleForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaMultipoleForce*>(target)->getMutualInducedTargetEpsilon();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon(OpenMM_AmoebaMultipoleForce* target, double inputMutualInducedTargetEpsilon) {
    reinterpret_cast<OpenMM::AmoebaMultipoleForce*>(target)->setMutualInducedTargetEpsilon(inputMutualInducedTargetEpsilon);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients(OpenMM_AmoebaMultipoleForce* target, const OpenMM_DoubleArray* coefficients) {
    reinterpret_cast<OpenMM::AmoebaMultipoleForce*>(target)->setExtrapolationCoefficients(*reinterpret_cast<const std::vector< double >*>(coefficients));
}
OPENMM_EXPORT_AMOEBA const OpenMM_DoubleArray* OpenMM_AmoebaMultipoleForce_getExtrapolationCoefficients(const OpenMM_AmoebaMultipoleForce* target) {
    const std::vector< double >* result = &reinterpret_cast<const OpenMM::AmoebaMultipoleForce*>(target)->getExtrapolationCoefficients();
    return reinterpret_cast<const OpenMM_DoubleArray*>(result);
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaMultipoleForce_getEwaldErrorTolerance(const OpenMM_AmoebaMultipoleForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaMultipoleForce*>(target)->getEwaldErrorTolerance();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance(OpenMM_AmoebaMultipoleForce* target, double tol) {
    reinterpret_cast<OpenMM::AmoebaMultipoleForce*>(target)->setEwaldErrorTolerance(tol);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_getLabFramePermanentDipoles(OpenMM_AmoebaMultipoleForce* target, OpenMM_Context* context, OpenMM_Vec3Array* dipoles) {
    reinterpret_cast<OpenMM::AmoebaMultipoleForce*>(target)->getLabFramePermanentDipoles(*reinterpret_cast<Context*>(context), *reinterpret_cast<std::vector< Vec3 >*>(dipoles));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_getInducedDipoles(OpenMM_AmoebaMultipoleForce* target, OpenMM_Context* context, OpenMM_Vec3Array* dipoles) {
    reinterpret_cast<OpenMM::AmoebaMultipoleForce*>(target)->getInducedDipoles(*reinterpret_cast<Context*>(context), *reinterpret_cast<std::vector< Vec3 >*>(dipoles));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_getTotalDipoles(OpenMM_AmoebaMultipoleForce* target, OpenMM_Context* context, OpenMM_Vec3Array* dipoles) {
    reinterpret_cast<OpenMM::AmoebaMultipoleForce*>(target)->getTotalDipoles(*reinterpret_cast<Context*>(context), *reinterpret_cast<std::vector< Vec3 >*>(dipoles));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_getElectrostaticPotential(OpenMM_AmoebaMultipoleForce* target, const OpenMM_Vec3Array* inputGrid, OpenMM_Context* context, OpenMM_DoubleArray* outputElectrostaticPotential) {
    reinterpret_cast<OpenMM::AmoebaMultipoleForce*>(target)->getElectrostaticPotential(*reinterpret_cast<const std::vector< Vec3 >*>(inputGrid), *reinterpret_cast<Context*>(context), *reinterpret_cast<std::vector< double >*>(outputElectrostaticPotential));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_getSystemMultipoleMoments(OpenMM_AmoebaMultipoleForce* target, OpenMM_Context* context, OpenMM_DoubleArray* outputMultipoleMoments) {
    reinterpret_cast<OpenMM::AmoebaMultipoleForce*>(target)->getSystemMultipoleMoments(*reinterpret_cast<Context*>(context), *reinterpret_cast<std::vector< double >*>(outputMultipoleMoments));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_updateParametersInContext(OpenMM_AmoebaMultipoleForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::AmoebaMultipoleForce*>(target)->updateParametersInContext(*reinterpret_cast<Context*>(context));
}
OPENMM_EXPORT_AMOEBA OpenMM_Boolean OpenMM_AmoebaMultipoleForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaMultipoleForce* target) {
    bool result = reinterpret_cast<const OpenMM::AmoebaMultipoleForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::AmoebaStretchBendForce */
OPENMM_EXPORT_AMOEBA OpenMM_AmoebaStretchBendForce* OpenMM_AmoebaStretchBendForce_create() {
    return reinterpret_cast<OpenMM_AmoebaStretchBendForce*>(new OpenMM::AmoebaStretchBendForce());
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaStretchBendForce_destroy(OpenMM_AmoebaStretchBendForce* target) {
    delete reinterpret_cast<OpenMM::AmoebaStretchBendForce*>(target);
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaStretchBendForce_getNumStretchBends(const OpenMM_AmoebaStretchBendForce* target) {
    int result = reinterpret_cast<const OpenMM::AmoebaStretchBendForce*>(target)->getNumStretchBends();
    return result;
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaStretchBendForce_addStretchBend(OpenMM_AmoebaStretchBendForce* target, int particle1, int particle2, int particle3, double lengthAB, double lengthCB, double angle, double k1, double k2) {
    int result = reinterpret_cast<OpenMM::AmoebaStretchBendForce*>(target)->addStretchBend(particle1, particle2, particle3, lengthAB, lengthCB, angle, k1, k2);
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaStretchBendForce_getStretchBendParameters(const OpenMM_AmoebaStretchBendForce* target, int index, int* particle1, int* particle2, int* particle3, double* lengthAB, double* lengthCB, double* angle, double* k1, double* k2) {
    reinterpret_cast<const OpenMM::AmoebaStretchBendForce*>(target)->getStretchBendParameters(index, *reinterpret_cast<int*>(particle1), *reinterpret_cast<int*>(particle2), *reinterpret_cast<int*>(particle3), *reinterpret_cast<double*>(lengthAB), *reinterpret_cast<double*>(lengthCB), *reinterpret_cast<double*>(angle), *reinterpret_cast<double*>(k1), *reinterpret_cast<double*>(k2));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaStretchBendForce_setStretchBendParameters(OpenMM_AmoebaStretchBendForce* target, int index, int particle1, int particle2, int particle3, double lengthAB, double lengthCB, double angle, double k1, double k2) {
    reinterpret_cast<OpenMM::AmoebaStretchBendForce*>(target)->setStretchBendParameters(index, particle1, particle2, particle3, lengthAB, lengthCB, angle, k1, k2);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaStretchBendForce_updateParametersInContext(OpenMM_AmoebaStretchBendForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::AmoebaStretchBendForce*>(target)->updateParametersInContext(*reinterpret_cast<Context*>(context));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaStretchBendForce_setUsesPeriodicBoundaryConditions(OpenMM_AmoebaStretchBendForce* target, OpenMM_Boolean periodic) {
    reinterpret_cast<OpenMM::AmoebaStretchBendForce*>(target)->setUsesPeriodicBoundaryConditions(periodic);
}
OPENMM_EXPORT_AMOEBA OpenMM_Boolean OpenMM_AmoebaStretchBendForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaStretchBendForce* target) {
    bool result = reinterpret_cast<const OpenMM::AmoebaStretchBendForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::AmoebaVdwForce */
OPENMM_EXPORT_AMOEBA OpenMM_AmoebaVdwForce* OpenMM_AmoebaVdwForce_create() {
    return reinterpret_cast<OpenMM_AmoebaVdwForce*>(new OpenMM::AmoebaVdwForce());
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_destroy(OpenMM_AmoebaVdwForce* target) {
    delete reinterpret_cast<OpenMM::AmoebaVdwForce*>(target);
}
OPENMM_EXPORT_AMOEBA const char* OpenMM_AmoebaVdwForce_Lambda() {
    const std::string* result = &OpenMM::AmoebaVdwForce::Lambda();
    return result->c_str();
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaVdwForce_getNumParticles(const OpenMM_AmoebaVdwForce* target) {
    int result = reinterpret_cast<const OpenMM::AmoebaVdwForce*>(target)->getNumParticles();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_setParticleParameters(OpenMM_AmoebaVdwForce* target, int particleIndex, int parentIndex, double sigma, double epsilon, double reductionFactor, OpenMM_Boolean isAlchemical) {
    reinterpret_cast<OpenMM::AmoebaVdwForce*>(target)->setParticleParameters(particleIndex, parentIndex, sigma, epsilon, reductionFactor, isAlchemical);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_getParticleParameters(const OpenMM_AmoebaVdwForce* target, int particleIndex, int* parentIndex, double* sigma, double* epsilon, double* reductionFactor, OpenMM_Boolean* isAlchemical) {
    reinterpret_cast<const OpenMM::AmoebaVdwForce*>(target)->getParticleParameters(particleIndex, *reinterpret_cast<int*>(parentIndex), *reinterpret_cast<double*>(sigma), *reinterpret_cast<double*>(epsilon), *reinterpret_cast<double*>(reductionFactor), *reinterpret_cast<bool*>(isAlchemical));
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaVdwForce_addParticle(OpenMM_AmoebaVdwForce* target, int parentIndex, double sigma, double epsilon, double reductionFactor, OpenMM_Boolean isAlchemical) {
    int result = reinterpret_cast<OpenMM::AmoebaVdwForce*>(target)->addParticle(parentIndex, sigma, epsilon, reductionFactor, isAlchemical);
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_setSigmaCombiningRule(OpenMM_AmoebaVdwForce* target, const char* sigmaCombiningRule) {
    reinterpret_cast<OpenMM::AmoebaVdwForce*>(target)->setSigmaCombiningRule(std::string(sigmaCombiningRule));
}
OPENMM_EXPORT_AMOEBA const char* OpenMM_AmoebaVdwForce_getSigmaCombiningRule(const OpenMM_AmoebaVdwForce* target) {
    const std::string* result = &reinterpret_cast<const OpenMM::AmoebaVdwForce*>(target)->getSigmaCombiningRule();
    return result->c_str();
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_setEpsilonCombiningRule(OpenMM_AmoebaVdwForce* target, const char* epsilonCombiningRule) {
    reinterpret_cast<OpenMM::AmoebaVdwForce*>(target)->setEpsilonCombiningRule(std::string(epsilonCombiningRule));
}
OPENMM_EXPORT_AMOEBA const char* OpenMM_AmoebaVdwForce_getEpsilonCombiningRule(const OpenMM_AmoebaVdwForce* target) {
    const std::string* result = &reinterpret_cast<const OpenMM::AmoebaVdwForce*>(target)->getEpsilonCombiningRule();
    return result->c_str();
}
OPENMM_EXPORT_AMOEBA OpenMM_Boolean OpenMM_AmoebaVdwForce_getUseDispersionCorrection(const OpenMM_AmoebaVdwForce* target) {
    bool result = reinterpret_cast<const OpenMM::AmoebaVdwForce*>(target)->getUseDispersionCorrection();
    return (result ? OpenMM_True : OpenMM_False);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_setUseDispersionCorrection(OpenMM_AmoebaVdwForce* target, OpenMM_Boolean useCorrection) {
    reinterpret_cast<OpenMM::AmoebaVdwForce*>(target)->setUseDispersionCorrection(useCorrection);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_setParticleExclusions(OpenMM_AmoebaVdwForce* target, int particleIndex, const OpenMM_IntArray* exclusions) {
    reinterpret_cast<OpenMM::AmoebaVdwForce*>(target)->setParticleExclusions(particleIndex, *reinterpret_cast<const std::vector< int >*>(exclusions));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_getParticleExclusions(const OpenMM_AmoebaVdwForce* target, int particleIndex, OpenMM_IntArray* exclusions) {
    reinterpret_cast<const OpenMM::AmoebaVdwForce*>(target)->getParticleExclusions(particleIndex, *reinterpret_cast<std::vector< int >*>(exclusions));
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaVdwForce_getCutoffDistance(const OpenMM_AmoebaVdwForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaVdwForce*>(target)->getCutoffDistance();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_setCutoffDistance(OpenMM_AmoebaVdwForce* target, double distance) {
    reinterpret_cast<OpenMM::AmoebaVdwForce*>(target)->setCutoffDistance(distance);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_setCutoff(OpenMM_AmoebaVdwForce* target, double cutoff) {
    reinterpret_cast<OpenMM::AmoebaVdwForce*>(target)->setCutoff(cutoff);
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaVdwForce_getCutoff(const OpenMM_AmoebaVdwForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaVdwForce*>(target)->getCutoff();
    return result;
}
OPENMM_EXPORT_AMOEBA OpenMM_AmoebaVdwForce_NonbondedMethod OpenMM_AmoebaVdwForce_getNonbondedMethod(const OpenMM_AmoebaVdwForce* target) {
    OpenMM::AmoebaVdwForce::NonbondedMethod result = reinterpret_cast<const OpenMM::AmoebaVdwForce*>(target)->getNonbondedMethod();
    return static_cast<OpenMM_AmoebaVdwForce_NonbondedMethod>(result);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_setNonbondedMethod(OpenMM_AmoebaVdwForce* target, OpenMM_AmoebaVdwForce_NonbondedMethod method) {
    reinterpret_cast<OpenMM::AmoebaVdwForce*>(target)->setNonbondedMethod(static_cast<OpenMM::AmoebaVdwForce::NonbondedMethod>(method));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_setSoftcorePower(OpenMM_AmoebaVdwForce* target, int n) {
    reinterpret_cast<OpenMM::AmoebaVdwForce*>(target)->setSoftcorePower(n);
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaVdwForce_getSoftcorePower(const OpenMM_AmoebaVdwForce* target) {
    int result = reinterpret_cast<const OpenMM::AmoebaVdwForce*>(target)->getSoftcorePower();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_setSoftcoreAlpha(OpenMM_AmoebaVdwForce* target, double alpha) {
    reinterpret_cast<OpenMM::AmoebaVdwForce*>(target)->setSoftcoreAlpha(alpha);
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaVdwForce_getSoftcoreAlpha(const OpenMM_AmoebaVdwForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaVdwForce*>(target)->getSoftcoreAlpha();
    return result;
}
OPENMM_EXPORT_AMOEBA OpenMM_AmoebaVdwForce_AlchemicalMethod OpenMM_AmoebaVdwForce_getAlchemicalMethod(const OpenMM_AmoebaVdwForce* target) {
    OpenMM::AmoebaVdwForce::AlchemicalMethod result = reinterpret_cast<const OpenMM::AmoebaVdwForce*>(target)->getAlchemicalMethod();
    return static_cast<OpenMM_AmoebaVdwForce_AlchemicalMethod>(result);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_setAlchemicalMethod(OpenMM_AmoebaVdwForce* target, OpenMM_AmoebaVdwForce_AlchemicalMethod method) {
    reinterpret_cast<OpenMM::AmoebaVdwForce*>(target)->setAlchemicalMethod(static_cast<OpenMM::AmoebaVdwForce::AlchemicalMethod>(method));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_updateParametersInContext(OpenMM_AmoebaVdwForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::AmoebaVdwForce*>(target)->updateParametersInContext(*reinterpret_cast<Context*>(context));
}
OPENMM_EXPORT_AMOEBA OpenMM_Boolean OpenMM_AmoebaVdwForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaVdwForce* target) {
    bool result = reinterpret_cast<const OpenMM::AmoebaVdwForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::AmoebaGeneralizedKirkwoodForce */
OPENMM_EXPORT_AMOEBA OpenMM_AmoebaGeneralizedKirkwoodForce* OpenMM_AmoebaGeneralizedKirkwoodForce_create() {
    return reinterpret_cast<OpenMM_AmoebaGeneralizedKirkwoodForce*>(new OpenMM::AmoebaGeneralizedKirkwoodForce());
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaGeneralizedKirkwoodForce_destroy(OpenMM_AmoebaGeneralizedKirkwoodForce* target) {
    delete reinterpret_cast<OpenMM::AmoebaGeneralizedKirkwoodForce*>(target);
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaGeneralizedKirkwoodForce_getNumParticles(const OpenMM_AmoebaGeneralizedKirkwoodForce* target) {
    int result = reinterpret_cast<const OpenMM::AmoebaGeneralizedKirkwoodForce*>(target)->getNumParticles();
    return result;
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle(OpenMM_AmoebaGeneralizedKirkwoodForce* target, double charge, double radius, double scalingFactor) {
    int result = reinterpret_cast<OpenMM::AmoebaGeneralizedKirkwoodForce*>(target)->addParticle(charge, radius, scalingFactor);
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaGeneralizedKirkwoodForce_getParticleParameters(const OpenMM_AmoebaGeneralizedKirkwoodForce* target, int index, double* charge, double* radius, double* scalingFactor) {
    reinterpret_cast<const OpenMM::AmoebaGeneralizedKirkwoodForce*>(target)->getParticleParameters(index, *reinterpret_cast<double*>(charge), *reinterpret_cast<double*>(radius), *reinterpret_cast<double*>(scalingFactor));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters(OpenMM_AmoebaGeneralizedKirkwoodForce* target, int index, double charge, double radius, double scalingFactor) {
    reinterpret_cast<OpenMM::AmoebaGeneralizedKirkwoodForce*>(target)->setParticleParameters(index, charge, radius, scalingFactor);
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaGeneralizedKirkwoodForce_getSolventDielectric(const OpenMM_AmoebaGeneralizedKirkwoodForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaGeneralizedKirkwoodForce*>(target)->getSolventDielectric();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric(OpenMM_AmoebaGeneralizedKirkwoodForce* target, double dielectric) {
    reinterpret_cast<OpenMM::AmoebaGeneralizedKirkwoodForce*>(target)->setSolventDielectric(dielectric);
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaGeneralizedKirkwoodForce_getSoluteDielectric(const OpenMM_AmoebaGeneralizedKirkwoodForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaGeneralizedKirkwoodForce*>(target)->getSoluteDielectric();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric(OpenMM_AmoebaGeneralizedKirkwoodForce* target, double dielectric) {
    reinterpret_cast<OpenMM::AmoebaGeneralizedKirkwoodForce*>(target)->setSoluteDielectric(dielectric);
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaGeneralizedKirkwoodForce_getIncludeCavityTerm(const OpenMM_AmoebaGeneralizedKirkwoodForce* target) {
    int result = reinterpret_cast<const OpenMM::AmoebaGeneralizedKirkwoodForce*>(target)->getIncludeCavityTerm();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm(OpenMM_AmoebaGeneralizedKirkwoodForce* target, int includeCavityTerm) {
    reinterpret_cast<OpenMM::AmoebaGeneralizedKirkwoodForce*>(target)->setIncludeCavityTerm(includeCavityTerm);
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaGeneralizedKirkwoodForce_getProbeRadius(const OpenMM_AmoebaGeneralizedKirkwoodForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaGeneralizedKirkwoodForce*>(target)->getProbeRadius();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius(OpenMM_AmoebaGeneralizedKirkwoodForce* target, double probeRadius) {
    reinterpret_cast<OpenMM::AmoebaGeneralizedKirkwoodForce*>(target)->setProbeRadius(probeRadius);
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaGeneralizedKirkwoodForce_getSurfaceAreaFactor(const OpenMM_AmoebaGeneralizedKirkwoodForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaGeneralizedKirkwoodForce*>(target)->getSurfaceAreaFactor();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor(OpenMM_AmoebaGeneralizedKirkwoodForce* target, double surfaceAreaFactor) {
    reinterpret_cast<OpenMM::AmoebaGeneralizedKirkwoodForce*>(target)->setSurfaceAreaFactor(surfaceAreaFactor);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext(OpenMM_AmoebaGeneralizedKirkwoodForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::AmoebaGeneralizedKirkwoodForce*>(target)->updateParametersInContext(*reinterpret_cast<Context*>(context));
}
OPENMM_EXPORT_AMOEBA OpenMM_Boolean OpenMM_AmoebaGeneralizedKirkwoodForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaGeneralizedKirkwoodForce* target) {
    bool result = reinterpret_cast<const OpenMM::AmoebaGeneralizedKirkwoodForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::AmoebaOutOfPlaneBendForce */
OPENMM_EXPORT_AMOEBA OpenMM_AmoebaOutOfPlaneBendForce* OpenMM_AmoebaOutOfPlaneBendForce_create() {
    return reinterpret_cast<OpenMM_AmoebaOutOfPlaneBendForce*>(new OpenMM::AmoebaOutOfPlaneBendForce());
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaOutOfPlaneBendForce_destroy(OpenMM_AmoebaOutOfPlaneBendForce* target) {
    delete reinterpret_cast<OpenMM::AmoebaOutOfPlaneBendForce*>(target);
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaOutOfPlaneBendForce_getNumOutOfPlaneBends(const OpenMM_AmoebaOutOfPlaneBendForce* target) {
    int result = reinterpret_cast<const OpenMM::AmoebaOutOfPlaneBendForce*>(target)->getNumOutOfPlaneBends();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendCubic(OpenMM_AmoebaOutOfPlaneBendForce* target, double cubicK) {
    reinterpret_cast<OpenMM::AmoebaOutOfPlaneBendForce*>(target)->setAmoebaGlobalOutOfPlaneBendCubic(cubicK);
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaOutOfPlaneBendForce_getAmoebaGlobalOutOfPlaneBendCubic(const OpenMM_AmoebaOutOfPlaneBendForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaOutOfPlaneBendForce*>(target)->getAmoebaGlobalOutOfPlaneBendCubic();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendQuartic(OpenMM_AmoebaOutOfPlaneBendForce* target, double quarticK) {
    reinterpret_cast<OpenMM::AmoebaOutOfPlaneBendForce*>(target)->setAmoebaGlobalOutOfPlaneBendQuartic(quarticK);
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaOutOfPlaneBendForce_getAmoebaGlobalOutOfPlaneBendQuartic(const OpenMM_AmoebaOutOfPlaneBendForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaOutOfPlaneBendForce*>(target)->getAmoebaGlobalOutOfPlaneBendQuartic();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendPentic(OpenMM_AmoebaOutOfPlaneBendForce* target, double penticK) {
    reinterpret_cast<OpenMM::AmoebaOutOfPlaneBendForce*>(target)->setAmoebaGlobalOutOfPlaneBendPentic(penticK);
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaOutOfPlaneBendForce_getAmoebaGlobalOutOfPlaneBendPentic(const OpenMM_AmoebaOutOfPlaneBendForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaOutOfPlaneBendForce*>(target)->getAmoebaGlobalOutOfPlaneBendPentic();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendSextic(OpenMM_AmoebaOutOfPlaneBendForce* target, double sexticK) {
    reinterpret_cast<OpenMM::AmoebaOutOfPlaneBendForce*>(target)->setAmoebaGlobalOutOfPlaneBendSextic(sexticK);
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaOutOfPlaneBendForce_getAmoebaGlobalOutOfPlaneBendSextic(const OpenMM_AmoebaOutOfPlaneBendForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaOutOfPlaneBendForce*>(target)->getAmoebaGlobalOutOfPlaneBendSextic();
    return result;
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaOutOfPlaneBendForce_addOutOfPlaneBend(OpenMM_AmoebaOutOfPlaneBendForce* target, int particle1, int particle2, int particle3, int particle4, double k) {
    int result = reinterpret_cast<OpenMM::AmoebaOutOfPlaneBendForce*>(target)->addOutOfPlaneBend(particle1, particle2, particle3, particle4, k);
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaOutOfPlaneBendForce_getOutOfPlaneBendParameters(const OpenMM_AmoebaOutOfPlaneBendForce* target, int index, int* particle1, int* particle2, int* particle3, int* particle4, double* k) {
    reinterpret_cast<const OpenMM::AmoebaOutOfPlaneBendForce*>(target)->getOutOfPlaneBendParameters(index, *reinterpret_cast<int*>(particle1), *reinterpret_cast<int*>(particle2), *reinterpret_cast<int*>(particle3), *reinterpret_cast<int*>(particle4), *reinterpret_cast<double*>(k));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaOutOfPlaneBendForce_setOutOfPlaneBendParameters(OpenMM_AmoebaOutOfPlaneBendForce* target, int index, int particle1, int particle2, int particle3, int particle4, double k) {
    reinterpret_cast<OpenMM::AmoebaOutOfPlaneBendForce*>(target)->setOutOfPlaneBendParameters(index, particle1, particle2, particle3, particle4, k);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaOutOfPlaneBendForce_updateParametersInContext(OpenMM_AmoebaOutOfPlaneBendForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::AmoebaOutOfPlaneBendForce*>(target)->updateParametersInContext(*reinterpret_cast<Context*>(context));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaOutOfPlaneBendForce_setUsesPeriodicBoundaryConditions(OpenMM_AmoebaOutOfPlaneBendForce* target, OpenMM_Boolean periodic) {
    reinterpret_cast<OpenMM::AmoebaOutOfPlaneBendForce*>(target)->setUsesPeriodicBoundaryConditions(periodic);
}
OPENMM_EXPORT_AMOEBA OpenMM_Boolean OpenMM_AmoebaOutOfPlaneBendForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaOutOfPlaneBendForce* target) {
    bool result = reinterpret_cast<const OpenMM::AmoebaOutOfPlaneBendForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::HippoNonbondedForce */
OPENMM_EXPORT_AMOEBA OpenMM_HippoNonbondedForce* OpenMM_HippoNonbondedForce_create() {
    return reinterpret_cast<OpenMM_HippoNonbondedForce*>(new OpenMM::HippoNonbondedForce());
}
OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_destroy(OpenMM_HippoNonbondedForce* target) {
    delete reinterpret_cast<OpenMM::HippoNonbondedForce*>(target);
}
OPENMM_EXPORT_AMOEBA int OpenMM_HippoNonbondedForce_getNumParticles(const OpenMM_HippoNonbondedForce* target) {
    int result = reinterpret_cast<const OpenMM::HippoNonbondedForce*>(target)->getNumParticles();
    return result;
}
OPENMM_EXPORT_AMOEBA int OpenMM_HippoNonbondedForce_getNumExceptions(const OpenMM_HippoNonbondedForce* target) {
    int result = reinterpret_cast<const OpenMM::HippoNonbondedForce*>(target)->getNumExceptions();
    return result;
}
OPENMM_EXPORT_AMOEBA OpenMM_HippoNonbondedForce_NonbondedMethod OpenMM_HippoNonbondedForce_getNonbondedMethod(const OpenMM_HippoNonbondedForce* target) {
    OpenMM::HippoNonbondedForce::NonbondedMethod result = reinterpret_cast<const OpenMM::HippoNonbondedForce*>(target)->getNonbondedMethod();
    return static_cast<OpenMM_HippoNonbondedForce_NonbondedMethod>(result);
}
OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_setNonbondedMethod(OpenMM_HippoNonbondedForce* target, OpenMM_HippoNonbondedForce_NonbondedMethod method) {
    reinterpret_cast<OpenMM::HippoNonbondedForce*>(target)->setNonbondedMethod(static_cast<OpenMM::HippoNonbondedForce::NonbondedMethod>(method));
}
OPENMM_EXPORT_AMOEBA double OpenMM_HippoNonbondedForce_getCutoffDistance(const OpenMM_HippoNonbondedForce* target) {
    double result = reinterpret_cast<const OpenMM::HippoNonbondedForce*>(target)->getCutoffDistance();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_setCutoffDistance(OpenMM_HippoNonbondedForce* target, double distance) {
    reinterpret_cast<OpenMM::HippoNonbondedForce*>(target)->setCutoffDistance(distance);
}
OPENMM_EXPORT_AMOEBA double OpenMM_HippoNonbondedForce_getSwitchingDistance(const OpenMM_HippoNonbondedForce* target) {
    double result = reinterpret_cast<const OpenMM::HippoNonbondedForce*>(target)->getSwitchingDistance();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_setSwitchingDistance(OpenMM_HippoNonbondedForce* target, double distance) {
    reinterpret_cast<OpenMM::HippoNonbondedForce*>(target)->setSwitchingDistance(distance);
}
OPENMM_EXPORT_AMOEBA const OpenMM_DoubleArray* OpenMM_HippoNonbondedForce_getExtrapolationCoefficients(const OpenMM_HippoNonbondedForce* target) {
    const std::vector< double >* result = &reinterpret_cast<const OpenMM::HippoNonbondedForce*>(target)->getExtrapolationCoefficients();
    return reinterpret_cast<const OpenMM_DoubleArray*>(result);
}
OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_setExtrapolationCoefficients(OpenMM_HippoNonbondedForce* target, const OpenMM_DoubleArray* coefficients) {
    reinterpret_cast<OpenMM::HippoNonbondedForce*>(target)->setExtrapolationCoefficients(*reinterpret_cast<const std::vector< double >*>(coefficients));
}
OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_getPMEParameters(const OpenMM_HippoNonbondedForce* target, double* alpha, int* nx, int* ny, int* nz) {
    reinterpret_cast<const OpenMM::HippoNonbondedForce*>(target)->getPMEParameters(*reinterpret_cast<double*>(alpha), *reinterpret_cast<int*>(nx), *reinterpret_cast<int*>(ny), *reinterpret_cast<int*>(nz));
}
OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_getDPMEParameters(const OpenMM_HippoNonbondedForce* target, double* alpha, int* nx, int* ny, int* nz) {
    reinterpret_cast<const OpenMM::HippoNonbondedForce*>(target)->getDPMEParameters(*reinterpret_cast<double*>(alpha), *reinterpret_cast<int*>(nx), *reinterpret_cast<int*>(ny), *reinterpret_cast<int*>(nz));
}
OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_setPMEParameters(OpenMM_HippoNonbondedForce* target, double alpha, int nx, int ny, int nz) {
    reinterpret_cast<OpenMM::HippoNonbondedForce*>(target)->setPMEParameters(alpha, nx, ny, nz);
}
OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_setDPMEParameters(OpenMM_HippoNonbondedForce* target, double alpha, int nx, int ny, int nz) {
    reinterpret_cast<OpenMM::HippoNonbondedForce*>(target)->setDPMEParameters(alpha, nx, ny, nz);
}
OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_getPMEParametersInContext(const OpenMM_HippoNonbondedForce* target, const OpenMM_Context* context, double* alpha, int* nx, int* ny, int* nz) {
    reinterpret_cast<const OpenMM::HippoNonbondedForce*>(target)->getPMEParametersInContext(*reinterpret_cast<const Context*>(context), *reinterpret_cast<double*>(alpha), *reinterpret_cast<int*>(nx), *reinterpret_cast<int*>(ny), *reinterpret_cast<int*>(nz));
}
OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_getDPMEParametersInContext(const OpenMM_HippoNonbondedForce* target, const OpenMM_Context* context, double* alpha, int* nx, int* ny, int* nz) {
    reinterpret_cast<const OpenMM::HippoNonbondedForce*>(target)->getDPMEParametersInContext(*reinterpret_cast<const Context*>(context), *reinterpret_cast<double*>(alpha), *reinterpret_cast<int*>(nx), *reinterpret_cast<int*>(ny), *reinterpret_cast<int*>(nz));
}
OPENMM_EXPORT_AMOEBA int OpenMM_HippoNonbondedForce_addParticle(OpenMM_HippoNonbondedForce* target, double charge, const OpenMM_DoubleArray* dipole, const OpenMM_DoubleArray* quadrupole, double coreCharge, double alpha, double epsilon, double damping, double c6, double pauliK, double pauliQ, double pauliAlpha, double polarizability, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY) {
    int result = reinterpret_cast<OpenMM::HippoNonbondedForce*>(target)->addParticle(charge, *reinterpret_cast<const std::vector< double >*>(dipole), *reinterpret_cast<const std::vector< double >*>(quadrupole), coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha, polarizability, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY);
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_getParticleParameters(const OpenMM_HippoNonbondedForce* target, int index, double* charge, OpenMM_DoubleArray* dipole, OpenMM_DoubleArray* quadrupole, double* coreCharge, double* alpha, double* epsilon, double* damping, double* c6, double* pauliK, double* pauliQ, double* pauliAlpha, double* polarizability, int* axisType, int* multipoleAtomZ, int* multipoleAtomX, int* multipoleAtomY) {
    reinterpret_cast<const OpenMM::HippoNonbondedForce*>(target)->getParticleParameters(index, *reinterpret_cast<double*>(charge), *reinterpret_cast<std::vector< double >*>(dipole), *reinterpret_cast<std::vector< double >*>(quadrupole), *reinterpret_cast<double*>(coreCharge), *reinterpret_cast<double*>(alpha), *reinterpret_cast<double*>(epsilon), *reinterpret_cast<double*>(damping), *reinterpret_cast<double*>(c6), *reinterpret_cast<double*>(pauliK), *reinterpret_cast<double*>(pauliQ), *reinterpret_cast<double*>(pauliAlpha), *reinterpret_cast<double*>(polarizability), *reinterpret_cast<int*>(axisType), *reinterpret_cast<int*>(multipoleAtomZ), *reinterpret_cast<int*>(multipoleAtomX), *reinterpret_cast<int*>(multipoleAtomY));
}
OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_setParticleParameters(OpenMM_HippoNonbondedForce* target, int index, double charge, const OpenMM_DoubleArray* dipole, const OpenMM_DoubleArray* quadrupole, double coreCharge, double alpha, double epsilon, double damping, double c6, double pauliK, double pauliQ, double pauliAlpha, double polarizability, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY) {
    reinterpret_cast<OpenMM::HippoNonbondedForce*>(target)->setParticleParameters(index, charge, *reinterpret_cast<const std::vector< double >*>(dipole), *reinterpret_cast<const std::vector< double >*>(quadrupole), coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha, polarizability, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY);
}
OPENMM_EXPORT_AMOEBA int OpenMM_HippoNonbondedForce_addException(OpenMM_HippoNonbondedForce* target, int particle1, int particle2, double multipoleMultipoleScale, double dipoleMultipoleScale, double dipoleDipoleScale, double dispersionScale, double repulsionScale, double chargeTransferScale, OpenMM_Boolean replace) {
    int result = reinterpret_cast<OpenMM::HippoNonbondedForce*>(target)->addException(particle1, particle2, multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale, replace);
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_getExceptionParameters(const OpenMM_HippoNonbondedForce* target, int index, int* particle1, int* particle2, double* multipoleMultipoleScale, double* dipoleMultipoleScale, double* dipoleDipoleScale, double* dispersionScale, double* repulsionScale, double* chargeTransferScale) {
    reinterpret_cast<const OpenMM::HippoNonbondedForce*>(target)->getExceptionParameters(index, *reinterpret_cast<int*>(particle1), *reinterpret_cast<int*>(particle2), *reinterpret_cast<double*>(multipoleMultipoleScale), *reinterpret_cast<double*>(dipoleMultipoleScale), *reinterpret_cast<double*>(dipoleDipoleScale), *reinterpret_cast<double*>(dispersionScale), *reinterpret_cast<double*>(repulsionScale), *reinterpret_cast<double*>(chargeTransferScale));
}
OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_setExceptionParameters(OpenMM_HippoNonbondedForce* target, int index, int particle1, int particle2, double multipoleMultipoleScale, double dipoleMultipoleScale, double dipoleDipoleScale, double dispersionScale, double repulsionScale, double chargeTransferScale) {
    reinterpret_cast<OpenMM::HippoNonbondedForce*>(target)->setExceptionParameters(index, particle1, particle2, multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale);
}
OPENMM_EXPORT_AMOEBA double OpenMM_HippoNonbondedForce_getEwaldErrorTolerance(const OpenMM_HippoNonbondedForce* target) {
    double result = reinterpret_cast<const OpenMM::HippoNonbondedForce*>(target)->getEwaldErrorTolerance();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_setEwaldErrorTolerance(OpenMM_HippoNonbondedForce* target, double tol) {
    reinterpret_cast<OpenMM::HippoNonbondedForce*>(target)->setEwaldErrorTolerance(tol);
}
OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_getLabFramePermanentDipoles(OpenMM_HippoNonbondedForce* target, OpenMM_Context* context, OpenMM_Vec3Array* dipoles) {
    reinterpret_cast<OpenMM::HippoNonbondedForce*>(target)->getLabFramePermanentDipoles(*reinterpret_cast<Context*>(context), *reinterpret_cast<std::vector< Vec3 >*>(dipoles));
}
OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_getInducedDipoles(OpenMM_HippoNonbondedForce* target, OpenMM_Context* context, OpenMM_Vec3Array* dipoles) {
    reinterpret_cast<OpenMM::HippoNonbondedForce*>(target)->getInducedDipoles(*reinterpret_cast<Context*>(context), *reinterpret_cast<std::vector< Vec3 >*>(dipoles));
}
OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_updateParametersInContext(OpenMM_HippoNonbondedForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::HippoNonbondedForce*>(target)->updateParametersInContext(*reinterpret_cast<Context*>(context));
}
OPENMM_EXPORT_AMOEBA OpenMM_Boolean OpenMM_HippoNonbondedForce_usesPeriodicBoundaryConditions(const OpenMM_HippoNonbondedForce* target) {
    bool result = reinterpret_cast<const OpenMM::HippoNonbondedForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::AmoebaAngleForce */
OPENMM_EXPORT_AMOEBA OpenMM_AmoebaAngleForce* OpenMM_AmoebaAngleForce_create() {
    return reinterpret_cast<OpenMM_AmoebaAngleForce*>(new OpenMM::AmoebaAngleForce());
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaAngleForce_destroy(OpenMM_AmoebaAngleForce* target) {
    delete reinterpret_cast<OpenMM::AmoebaAngleForce*>(target);
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaAngleForce_getNumAngles(const OpenMM_AmoebaAngleForce* target) {
    int result = reinterpret_cast<const OpenMM::AmoebaAngleForce*>(target)->getNumAngles();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleCubic(OpenMM_AmoebaAngleForce* target, double cubicK) {
    reinterpret_cast<OpenMM::AmoebaAngleForce*>(target)->setAmoebaGlobalAngleCubic(cubicK);
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaAngleForce_getAmoebaGlobalAngleCubic(const OpenMM_AmoebaAngleForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaAngleForce*>(target)->getAmoebaGlobalAngleCubic();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleQuartic(OpenMM_AmoebaAngleForce* target, double quarticK) {
    reinterpret_cast<OpenMM::AmoebaAngleForce*>(target)->setAmoebaGlobalAngleQuartic(quarticK);
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaAngleForce_getAmoebaGlobalAngleQuartic(const OpenMM_AmoebaAngleForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaAngleForce*>(target)->getAmoebaGlobalAngleQuartic();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaAngleForce_setAmoebaGlobalAnglePentic(OpenMM_AmoebaAngleForce* target, double penticK) {
    reinterpret_cast<OpenMM::AmoebaAngleForce*>(target)->setAmoebaGlobalAnglePentic(penticK);
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaAngleForce_getAmoebaGlobalAnglePentic(const OpenMM_AmoebaAngleForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaAngleForce*>(target)->getAmoebaGlobalAnglePentic();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleSextic(OpenMM_AmoebaAngleForce* target, double sexticK) {
    reinterpret_cast<OpenMM::AmoebaAngleForce*>(target)->setAmoebaGlobalAngleSextic(sexticK);
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaAngleForce_getAmoebaGlobalAngleSextic(const OpenMM_AmoebaAngleForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaAngleForce*>(target)->getAmoebaGlobalAngleSextic();
    return result;
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaAngleForce_addAngle(OpenMM_AmoebaAngleForce* target, int particle1, int particle2, int particle3, double length, double quadraticK) {
    int result = reinterpret_cast<OpenMM::AmoebaAngleForce*>(target)->addAngle(particle1, particle2, particle3, length, quadraticK);
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaAngleForce_getAngleParameters(const OpenMM_AmoebaAngleForce* target, int index, int* particle1, int* particle2, int* particle3, double* length, double* quadraticK) {
    reinterpret_cast<const OpenMM::AmoebaAngleForce*>(target)->getAngleParameters(index, *reinterpret_cast<int*>(particle1), *reinterpret_cast<int*>(particle2), *reinterpret_cast<int*>(particle3), *reinterpret_cast<double*>(length), *reinterpret_cast<double*>(quadraticK));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaAngleForce_setAngleParameters(OpenMM_AmoebaAngleForce* target, int index, int particle1, int particle2, int particle3, double length, double quadraticK) {
    reinterpret_cast<OpenMM::AmoebaAngleForce*>(target)->setAngleParameters(index, particle1, particle2, particle3, length, quadraticK);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaAngleForce_updateParametersInContext(OpenMM_AmoebaAngleForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::AmoebaAngleForce*>(target)->updateParametersInContext(*reinterpret_cast<Context*>(context));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaAngleForce_setUsesPeriodicBoundaryConditions(OpenMM_AmoebaAngleForce* target, OpenMM_Boolean periodic) {
    reinterpret_cast<OpenMM::AmoebaAngleForce*>(target)->setUsesPeriodicBoundaryConditions(periodic);
}
OPENMM_EXPORT_AMOEBA OpenMM_Boolean OpenMM_AmoebaAngleForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaAngleForce* target) {
    bool result = reinterpret_cast<const OpenMM::AmoebaAngleForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::AmoebaInPlaneAngleForce */
OPENMM_EXPORT_AMOEBA OpenMM_AmoebaInPlaneAngleForce* OpenMM_AmoebaInPlaneAngleForce_create() {
    return reinterpret_cast<OpenMM_AmoebaInPlaneAngleForce*>(new OpenMM::AmoebaInPlaneAngleForce());
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaInPlaneAngleForce_destroy(OpenMM_AmoebaInPlaneAngleForce* target) {
    delete reinterpret_cast<OpenMM::AmoebaInPlaneAngleForce*>(target);
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaInPlaneAngleForce_getNumAngles(const OpenMM_AmoebaInPlaneAngleForce* target) {
    int result = reinterpret_cast<const OpenMM::AmoebaInPlaneAngleForce*>(target)->getNumAngles();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleCubic(OpenMM_AmoebaInPlaneAngleForce* target, double cubicK) {
    reinterpret_cast<OpenMM::AmoebaInPlaneAngleForce*>(target)->setAmoebaGlobalInPlaneAngleCubic(cubicK);
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaInPlaneAngleForce_getAmoebaGlobalInPlaneAngleCubic(const OpenMM_AmoebaInPlaneAngleForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaInPlaneAngleForce*>(target)->getAmoebaGlobalInPlaneAngleCubic();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleQuartic(OpenMM_AmoebaInPlaneAngleForce* target, double quarticK) {
    reinterpret_cast<OpenMM::AmoebaInPlaneAngleForce*>(target)->setAmoebaGlobalInPlaneAngleQuartic(quarticK);
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaInPlaneAngleForce_getAmoebaGlobalInPlaneAngleQuartic(const OpenMM_AmoebaInPlaneAngleForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaInPlaneAngleForce*>(target)->getAmoebaGlobalInPlaneAngleQuartic();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAnglePentic(OpenMM_AmoebaInPlaneAngleForce* target, double penticK) {
    reinterpret_cast<OpenMM::AmoebaInPlaneAngleForce*>(target)->setAmoebaGlobalInPlaneAnglePentic(penticK);
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaInPlaneAngleForce_getAmoebaGlobalInPlaneAnglePentic(const OpenMM_AmoebaInPlaneAngleForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaInPlaneAngleForce*>(target)->getAmoebaGlobalInPlaneAnglePentic();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleSextic(OpenMM_AmoebaInPlaneAngleForce* target, double sexticK) {
    reinterpret_cast<OpenMM::AmoebaInPlaneAngleForce*>(target)->setAmoebaGlobalInPlaneAngleSextic(sexticK);
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaInPlaneAngleForce_getAmoebaGlobalInPlaneAngleSextic(const OpenMM_AmoebaInPlaneAngleForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaInPlaneAngleForce*>(target)->getAmoebaGlobalInPlaneAngleSextic();
    return result;
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaInPlaneAngleForce_addAngle(OpenMM_AmoebaInPlaneAngleForce* target, int particle1, int particle2, int particle3, int particle4, double length, double quadraticK) {
    int result = reinterpret_cast<OpenMM::AmoebaInPlaneAngleForce*>(target)->addAngle(particle1, particle2, particle3, particle4, length, quadraticK);
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaInPlaneAngleForce_getAngleParameters(const OpenMM_AmoebaInPlaneAngleForce* target, int index, int* particle1, int* particle2, int* particle3, int* particle4, double* length, double* quadraticK) {
    reinterpret_cast<const OpenMM::AmoebaInPlaneAngleForce*>(target)->getAngleParameters(index, *reinterpret_cast<int*>(particle1), *reinterpret_cast<int*>(particle2), *reinterpret_cast<int*>(particle3), *reinterpret_cast<int*>(particle4), *reinterpret_cast<double*>(length), *reinterpret_cast<double*>(quadraticK));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaInPlaneAngleForce_setAngleParameters(OpenMM_AmoebaInPlaneAngleForce* target, int index, int particle1, int particle2, int particle3, int particle4, double length, double quadraticK) {
    reinterpret_cast<OpenMM::AmoebaInPlaneAngleForce*>(target)->setAngleParameters(index, particle1, particle2, particle3, particle4, length, quadraticK);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaInPlaneAngleForce_updateParametersInContext(OpenMM_AmoebaInPlaneAngleForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::AmoebaInPlaneAngleForce*>(target)->updateParametersInContext(*reinterpret_cast<Context*>(context));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaInPlaneAngleForce_setUsesPeriodicBoundaryConditions(OpenMM_AmoebaInPlaneAngleForce* target, OpenMM_Boolean periodic) {
    reinterpret_cast<OpenMM::AmoebaInPlaneAngleForce*>(target)->setUsesPeriodicBoundaryConditions(periodic);
}
OPENMM_EXPORT_AMOEBA OpenMM_Boolean OpenMM_AmoebaInPlaneAngleForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaInPlaneAngleForce* target) {
    bool result = reinterpret_cast<const OpenMM::AmoebaInPlaneAngleForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::AmoebaBondForce */
OPENMM_EXPORT_AMOEBA OpenMM_AmoebaBondForce* OpenMM_AmoebaBondForce_create() {
    return reinterpret_cast<OpenMM_AmoebaBondForce*>(new OpenMM::AmoebaBondForce());
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaBondForce_destroy(OpenMM_AmoebaBondForce* target) {
    delete reinterpret_cast<OpenMM::AmoebaBondForce*>(target);
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaBondForce_getNumBonds(const OpenMM_AmoebaBondForce* target) {
    int result = reinterpret_cast<const OpenMM::AmoebaBondForce*>(target)->getNumBonds();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaBondForce_setAmoebaGlobalBondCubic(OpenMM_AmoebaBondForce* target, double cubicK) {
    reinterpret_cast<OpenMM::AmoebaBondForce*>(target)->setAmoebaGlobalBondCubic(cubicK);
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaBondForce_getAmoebaGlobalBondCubic(const OpenMM_AmoebaBondForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaBondForce*>(target)->getAmoebaGlobalBondCubic();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaBondForce_setAmoebaGlobalBondQuartic(OpenMM_AmoebaBondForce* target, double quarticK) {
    reinterpret_cast<OpenMM::AmoebaBondForce*>(target)->setAmoebaGlobalBondQuartic(quarticK);
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaBondForce_getAmoebaGlobalBondQuartic(const OpenMM_AmoebaBondForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaBondForce*>(target)->getAmoebaGlobalBondQuartic();
    return result;
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaBondForce_addBond(OpenMM_AmoebaBondForce* target, int particle1, int particle2, double length, double quadraticK) {
    int result = reinterpret_cast<OpenMM::AmoebaBondForce*>(target)->addBond(particle1, particle2, length, quadraticK);
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaBondForce_getBondParameters(const OpenMM_AmoebaBondForce* target, int index, int* particle1, int* particle2, double* length, double* quadraticK) {
    reinterpret_cast<const OpenMM::AmoebaBondForce*>(target)->getBondParameters(index, *reinterpret_cast<int*>(particle1), *reinterpret_cast<int*>(particle2), *reinterpret_cast<double*>(length), *reinterpret_cast<double*>(quadraticK));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaBondForce_setBondParameters(OpenMM_AmoebaBondForce* target, int index, int particle1, int particle2, double length, double quadraticK) {
    reinterpret_cast<OpenMM::AmoebaBondForce*>(target)->setBondParameters(index, particle1, particle2, length, quadraticK);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaBondForce_updateParametersInContext(OpenMM_AmoebaBondForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::AmoebaBondForce*>(target)->updateParametersInContext(*reinterpret_cast<Context*>(context));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaBondForce_setUsesPeriodicBoundaryConditions(OpenMM_AmoebaBondForce* target, OpenMM_Boolean periodic) {
    reinterpret_cast<OpenMM::AmoebaBondForce*>(target)->setUsesPeriodicBoundaryConditions(periodic);
}
OPENMM_EXPORT_AMOEBA OpenMM_Boolean OpenMM_AmoebaBondForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaBondForce* target) {
    bool result = reinterpret_cast<const OpenMM::AmoebaBondForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::AmoebaWcaDispersionForce */
OPENMM_EXPORT_AMOEBA OpenMM_AmoebaWcaDispersionForce* OpenMM_AmoebaWcaDispersionForce_create() {
    return reinterpret_cast<OpenMM_AmoebaWcaDispersionForce*>(new OpenMM::AmoebaWcaDispersionForce());
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaWcaDispersionForce_destroy(OpenMM_AmoebaWcaDispersionForce* target) {
    delete reinterpret_cast<OpenMM::AmoebaWcaDispersionForce*>(target);
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaWcaDispersionForce_getNumParticles(const OpenMM_AmoebaWcaDispersionForce* target) {
    int result = reinterpret_cast<const OpenMM::AmoebaWcaDispersionForce*>(target)->getNumParticles();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaWcaDispersionForce_setParticleParameters(OpenMM_AmoebaWcaDispersionForce* target, int particleIndex, double radius, double epsilon) {
    reinterpret_cast<OpenMM::AmoebaWcaDispersionForce*>(target)->setParticleParameters(particleIndex, radius, epsilon);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaWcaDispersionForce_getParticleParameters(const OpenMM_AmoebaWcaDispersionForce* target, int particleIndex, double* radius, double* epsilon) {
    reinterpret_cast<const OpenMM::AmoebaWcaDispersionForce*>(target)->getParticleParameters(particleIndex, *reinterpret_cast<double*>(radius), *reinterpret_cast<double*>(epsilon));
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaWcaDispersionForce_addParticle(OpenMM_AmoebaWcaDispersionForce* target, double radius, double epsilon) {
    int result = reinterpret_cast<OpenMM::AmoebaWcaDispersionForce*>(target)->addParticle(radius, epsilon);
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaWcaDispersionForce_updateParametersInContext(OpenMM_AmoebaWcaDispersionForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::AmoebaWcaDispersionForce*>(target)->updateParametersInContext(*reinterpret_cast<Context*>(context));
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaWcaDispersionForce_getEpso(const OpenMM_AmoebaWcaDispersionForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaWcaDispersionForce*>(target)->getEpso();
    return result;
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaWcaDispersionForce_getEpsh(const OpenMM_AmoebaWcaDispersionForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaWcaDispersionForce*>(target)->getEpsh();
    return result;
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaWcaDispersionForce_getRmino(const OpenMM_AmoebaWcaDispersionForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaWcaDispersionForce*>(target)->getRmino();
    return result;
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaWcaDispersionForce_getRminh(const OpenMM_AmoebaWcaDispersionForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaWcaDispersionForce*>(target)->getRminh();
    return result;
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaWcaDispersionForce_getAwater(const OpenMM_AmoebaWcaDispersionForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaWcaDispersionForce*>(target)->getAwater();
    return result;
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaWcaDispersionForce_getShctd(const OpenMM_AmoebaWcaDispersionForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaWcaDispersionForce*>(target)->getShctd();
    return result;
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaWcaDispersionForce_getDispoff(const OpenMM_AmoebaWcaDispersionForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaWcaDispersionForce*>(target)->getDispoff();
    return result;
}
OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaWcaDispersionForce_getSlevy(const OpenMM_AmoebaWcaDispersionForce* target) {
    double result = reinterpret_cast<const OpenMM::AmoebaWcaDispersionForce*>(target)->getSlevy();
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaWcaDispersionForce_setEpso(OpenMM_AmoebaWcaDispersionForce* target, double inputValue) {
    reinterpret_cast<OpenMM::AmoebaWcaDispersionForce*>(target)->setEpso(inputValue);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaWcaDispersionForce_setEpsh(OpenMM_AmoebaWcaDispersionForce* target, double inputValue) {
    reinterpret_cast<OpenMM::AmoebaWcaDispersionForce*>(target)->setEpsh(inputValue);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaWcaDispersionForce_setRmino(OpenMM_AmoebaWcaDispersionForce* target, double inputValue) {
    reinterpret_cast<OpenMM::AmoebaWcaDispersionForce*>(target)->setRmino(inputValue);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaWcaDispersionForce_setRminh(OpenMM_AmoebaWcaDispersionForce* target, double inputValue) {
    reinterpret_cast<OpenMM::AmoebaWcaDispersionForce*>(target)->setRminh(inputValue);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaWcaDispersionForce_setAwater(OpenMM_AmoebaWcaDispersionForce* target, double inputValue) {
    reinterpret_cast<OpenMM::AmoebaWcaDispersionForce*>(target)->setAwater(inputValue);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaWcaDispersionForce_setShctd(OpenMM_AmoebaWcaDispersionForce* target, double inputValue) {
    reinterpret_cast<OpenMM::AmoebaWcaDispersionForce*>(target)->setShctd(inputValue);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaWcaDispersionForce_setDispoff(OpenMM_AmoebaWcaDispersionForce* target, double inputValue) {
    reinterpret_cast<OpenMM::AmoebaWcaDispersionForce*>(target)->setDispoff(inputValue);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaWcaDispersionForce_setSlevy(OpenMM_AmoebaWcaDispersionForce* target, double inputValue) {
    reinterpret_cast<OpenMM::AmoebaWcaDispersionForce*>(target)->setSlevy(inputValue);
}
OPENMM_EXPORT_AMOEBA OpenMM_Boolean OpenMM_AmoebaWcaDispersionForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaWcaDispersionForce* target) {
    bool result = reinterpret_cast<const OpenMM::AmoebaWcaDispersionForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::AmoebaPiTorsionForce */
OPENMM_EXPORT_AMOEBA OpenMM_AmoebaPiTorsionForce* OpenMM_AmoebaPiTorsionForce_create() {
    return reinterpret_cast<OpenMM_AmoebaPiTorsionForce*>(new OpenMM::AmoebaPiTorsionForce());
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaPiTorsionForce_destroy(OpenMM_AmoebaPiTorsionForce* target) {
    delete reinterpret_cast<OpenMM::AmoebaPiTorsionForce*>(target);
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaPiTorsionForce_getNumPiTorsions(const OpenMM_AmoebaPiTorsionForce* target) {
    int result = reinterpret_cast<const OpenMM::AmoebaPiTorsionForce*>(target)->getNumPiTorsions();
    return result;
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaPiTorsionForce_addPiTorsion(OpenMM_AmoebaPiTorsionForce* target, int particle1, int particle2, int particle3, int particle4, int particle5, int particle6, double k) {
    int result = reinterpret_cast<OpenMM::AmoebaPiTorsionForce*>(target)->addPiTorsion(particle1, particle2, particle3, particle4, particle5, particle6, k);
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaPiTorsionForce_getPiTorsionParameters(const OpenMM_AmoebaPiTorsionForce* target, int index, int* particle1, int* particle2, int* particle3, int* particle4, int* particle5, int* particle6, double* k) {
    reinterpret_cast<const OpenMM::AmoebaPiTorsionForce*>(target)->getPiTorsionParameters(index, *reinterpret_cast<int*>(particle1), *reinterpret_cast<int*>(particle2), *reinterpret_cast<int*>(particle3), *reinterpret_cast<int*>(particle4), *reinterpret_cast<int*>(particle5), *reinterpret_cast<int*>(particle6), *reinterpret_cast<double*>(k));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaPiTorsionForce_setPiTorsionParameters(OpenMM_AmoebaPiTorsionForce* target, int index, int particle1, int particle2, int particle3, int particle4, int particle5, int particle6, double k) {
    reinterpret_cast<OpenMM::AmoebaPiTorsionForce*>(target)->setPiTorsionParameters(index, particle1, particle2, particle3, particle4, particle5, particle6, k);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaPiTorsionForce_updateParametersInContext(OpenMM_AmoebaPiTorsionForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::AmoebaPiTorsionForce*>(target)->updateParametersInContext(*reinterpret_cast<Context*>(context));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaPiTorsionForce_setUsesPeriodicBoundaryConditions(OpenMM_AmoebaPiTorsionForce* target, OpenMM_Boolean periodic) {
    reinterpret_cast<OpenMM::AmoebaPiTorsionForce*>(target)->setUsesPeriodicBoundaryConditions(periodic);
}
OPENMM_EXPORT_AMOEBA OpenMM_Boolean OpenMM_AmoebaPiTorsionForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaPiTorsionForce* target) {
    bool result = reinterpret_cast<const OpenMM::AmoebaPiTorsionForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::AmoebaTorsionTorsionForce */
OPENMM_EXPORT_AMOEBA OpenMM_AmoebaTorsionTorsionForce* OpenMM_AmoebaTorsionTorsionForce_create() {
    return reinterpret_cast<OpenMM_AmoebaTorsionTorsionForce*>(new OpenMM::AmoebaTorsionTorsionForce());
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaTorsionTorsionForce_destroy(OpenMM_AmoebaTorsionTorsionForce* target) {
    delete reinterpret_cast<OpenMM::AmoebaTorsionTorsionForce*>(target);
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsions(const OpenMM_AmoebaTorsionTorsionForce* target) {
    int result = reinterpret_cast<const OpenMM::AmoebaTorsionTorsionForce*>(target)->getNumTorsionTorsions();
    return result;
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsionGrids(const OpenMM_AmoebaTorsionTorsionForce* target) {
    int result = reinterpret_cast<const OpenMM::AmoebaTorsionTorsionForce*>(target)->getNumTorsionTorsionGrids();
    return result;
}
OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion(OpenMM_AmoebaTorsionTorsionForce* target, int particle1, int particle2, int particle3, int particle4, int particle5, int chiralCheckAtomIndex, int gridIndex) {
    int result = reinterpret_cast<OpenMM::AmoebaTorsionTorsionForce*>(target)->addTorsionTorsion(particle1, particle2, particle3, particle4, particle5, chiralCheckAtomIndex, gridIndex);
    return result;
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionParameters(const OpenMM_AmoebaTorsionTorsionForce* target, int index, int* particle1, int* particle2, int* particle3, int* particle4, int* particle5, int* chiralCheckAtomIndex, int* gridIndex) {
    reinterpret_cast<const OpenMM::AmoebaTorsionTorsionForce*>(target)->getTorsionTorsionParameters(index, *reinterpret_cast<int*>(particle1), *reinterpret_cast<int*>(particle2), *reinterpret_cast<int*>(particle3), *reinterpret_cast<int*>(particle4), *reinterpret_cast<int*>(particle5), *reinterpret_cast<int*>(chiralCheckAtomIndex), *reinterpret_cast<int*>(gridIndex));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionParameters(OpenMM_AmoebaTorsionTorsionForce* target, int index, int particle1, int particle2, int particle3, int particle4, int particle5, int chiralCheckAtomIndex, int gridIndex) {
    reinterpret_cast<OpenMM::AmoebaTorsionTorsionForce*>(target)->setTorsionTorsionParameters(index, particle1, particle2, particle3, particle4, particle5, chiralCheckAtomIndex, gridIndex);
}
OPENMM_EXPORT_AMOEBA const OpenMM_3D_DoubleArray* OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionGrid(const OpenMM_AmoebaTorsionTorsionForce* target, int index) {
    const std::vector< std::vector< std::vector< double > > >* result = &reinterpret_cast<const OpenMM::AmoebaTorsionTorsionForce*>(target)->getTorsionTorsionGrid(index);
    return reinterpret_cast<const OpenMM_3D_DoubleArray*>(result);
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid(OpenMM_AmoebaTorsionTorsionForce* target, int index, const OpenMM_3D_DoubleArray* grid) {
    reinterpret_cast<OpenMM::AmoebaTorsionTorsionForce*>(target)->setTorsionTorsionGrid(index, *reinterpret_cast<const std::vector< std::vector< std::vector< double > > >*>(grid));
}
OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaTorsionTorsionForce_setUsesPeriodicBoundaryConditions(OpenMM_AmoebaTorsionTorsionForce* target, OpenMM_Boolean periodic) {
    reinterpret_cast<OpenMM::AmoebaTorsionTorsionForce*>(target)->setUsesPeriodicBoundaryConditions(periodic);
}
OPENMM_EXPORT_AMOEBA OpenMM_Boolean OpenMM_AmoebaTorsionTorsionForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaTorsionTorsionForce* target) {
    bool result = reinterpret_cast<const OpenMM::AmoebaTorsionTorsionForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

}

