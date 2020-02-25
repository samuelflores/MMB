
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "../../../wrappers/OpenMMCWrapper.h"
#include "AmoebaOpenMMCWrapper.h"
#include <cstring>
#include <vector>

using namespace OpenMM;
using namespace std;

/* Utilities for dealing with Fortran's blank-padded strings. */
static void copyAndPadString(char* dest, const char* source, int length) {
    bool reachedEnd = false;
    for (int i = 0; i < length; i++) {
        if (source[i] == 0)
            reachedEnd = true;
        dest[i] = (reachedEnd ? ' ' : source[i]);
    }
}

static string makeString(const char* fsrc, int length) {
    while (length && fsrc[length-1]==' ')
        --length;
    return string(fsrc, length);
}

extern "C" {


/* OpenMM::AmoebaMultipoleForce */
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_create_(OpenMM_AmoebaMultipoleForce*& result) {
    result = OpenMM_AmoebaMultipoleForce_create();
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_CREATE(OpenMM_AmoebaMultipoleForce*& result) {
    result = OpenMM_AmoebaMultipoleForce_create();
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_destroy_(OpenMM_AmoebaMultipoleForce*& destroy) {
    OpenMM_AmoebaMultipoleForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_DESTROY(OpenMM_AmoebaMultipoleForce*& destroy) {
    OpenMM_AmoebaMultipoleForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT_AMOEBA int openmm_amoebamultipoleforce_getnummultipoles_(const OpenMM_AmoebaMultipoleForce*& target) {
    return OpenMM_AmoebaMultipoleForce_getNumMultipoles(target);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBAMULTIPOLEFORCE_GETNUMMULTIPOLES(const OpenMM_AmoebaMultipoleForce*& target) {
    return OpenMM_AmoebaMultipoleForce_getNumMultipoles(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_getnonbondedmethod_(const OpenMM_AmoebaMultipoleForce*& target, OpenMM_AmoebaMultipoleForce_NonbondedMethod& result) {
    result = OpenMM_AmoebaMultipoleForce_getNonbondedMethod(target);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_GETNONBONDEDMETHOD(const OpenMM_AmoebaMultipoleForce*& target, OpenMM_AmoebaMultipoleForce_NonbondedMethod& result) {
    result = OpenMM_AmoebaMultipoleForce_getNonbondedMethod(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_setnonbondedmethod_(OpenMM_AmoebaMultipoleForce*& target, OpenMM_AmoebaMultipoleForce_NonbondedMethod& method) {
    OpenMM_AmoebaMultipoleForce_setNonbondedMethod(target, method);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_SETNONBONDEDMETHOD(OpenMM_AmoebaMultipoleForce*& target, OpenMM_AmoebaMultipoleForce_NonbondedMethod& method) {
    OpenMM_AmoebaMultipoleForce_setNonbondedMethod(target, method);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_getpolarizationtype_(const OpenMM_AmoebaMultipoleForce*& target, OpenMM_AmoebaMultipoleForce_PolarizationType& result) {
    result = OpenMM_AmoebaMultipoleForce_getPolarizationType(target);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_GETPOLARIZATIONTYPE(const OpenMM_AmoebaMultipoleForce*& target, OpenMM_AmoebaMultipoleForce_PolarizationType& result) {
    result = OpenMM_AmoebaMultipoleForce_getPolarizationType(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_setpolarizationtype_(OpenMM_AmoebaMultipoleForce*& target, OpenMM_AmoebaMultipoleForce_PolarizationType& type) {
    OpenMM_AmoebaMultipoleForce_setPolarizationType(target, type);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_SETPOLARIZATIONTYPE(OpenMM_AmoebaMultipoleForce*& target, OpenMM_AmoebaMultipoleForce_PolarizationType& type) {
    OpenMM_AmoebaMultipoleForce_setPolarizationType(target, type);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebamultipoleforce_getcutoffdistance_(const OpenMM_AmoebaMultipoleForce*& target) {
    return OpenMM_AmoebaMultipoleForce_getCutoffDistance(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAMULTIPOLEFORCE_GETCUTOFFDISTANCE(const OpenMM_AmoebaMultipoleForce*& target) {
    return OpenMM_AmoebaMultipoleForce_getCutoffDistance(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_setcutoffdistance_(OpenMM_AmoebaMultipoleForce*& target, double const& distance) {
    OpenMM_AmoebaMultipoleForce_setCutoffDistance(target, distance);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_SETCUTOFFDISTANCE(OpenMM_AmoebaMultipoleForce*& target, double const& distance) {
    OpenMM_AmoebaMultipoleForce_setCutoffDistance(target, distance);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_getpmeparameters_(const OpenMM_AmoebaMultipoleForce*& target, double* alpha, int* nx, int* ny, int* nz) {
    OpenMM_AmoebaMultipoleForce_getPMEParameters(target, alpha, nx, ny, nz);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_GETPMEPARAMETERS(const OpenMM_AmoebaMultipoleForce*& target, double* alpha, int* nx, int* ny, int* nz) {
    OpenMM_AmoebaMultipoleForce_getPMEParameters(target, alpha, nx, ny, nz);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_setpmeparameters_(OpenMM_AmoebaMultipoleForce*& target, double const& alpha, int const& nx, int const& ny, int const& nz) {
    OpenMM_AmoebaMultipoleForce_setPMEParameters(target, alpha, nx, ny, nz);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_SETPMEPARAMETERS(OpenMM_AmoebaMultipoleForce*& target, double const& alpha, int const& nx, int const& ny, int const& nz) {
    OpenMM_AmoebaMultipoleForce_setPMEParameters(target, alpha, nx, ny, nz);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebamultipoleforce_getaewald_(const OpenMM_AmoebaMultipoleForce*& target) {
    return OpenMM_AmoebaMultipoleForce_getAEwald(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAMULTIPOLEFORCE_GETAEWALD(const OpenMM_AmoebaMultipoleForce*& target) {
    return OpenMM_AmoebaMultipoleForce_getAEwald(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_setaewald_(OpenMM_AmoebaMultipoleForce*& target, double const& aewald) {
    OpenMM_AmoebaMultipoleForce_setAEwald(target, aewald);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_SETAEWALD(OpenMM_AmoebaMultipoleForce*& target, double const& aewald) {
    OpenMM_AmoebaMultipoleForce_setAEwald(target, aewald);
}
OPENMM_EXPORT_AMOEBA int openmm_amoebamultipoleforce_getpmebsplineorder_(const OpenMM_AmoebaMultipoleForce*& target) {
    return OpenMM_AmoebaMultipoleForce_getPmeBSplineOrder(target);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBAMULTIPOLEFORCE_GETPMEBSPLINEORDER(const OpenMM_AmoebaMultipoleForce*& target) {
    return OpenMM_AmoebaMultipoleForce_getPmeBSplineOrder(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_getpmegriddimensions_(const OpenMM_AmoebaMultipoleForce*& target, OpenMM_IntArray*& gridDimension) {
    OpenMM_AmoebaMultipoleForce_getPmeGridDimensions(target, gridDimension);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_GETPMEGRIDDIMENSIONS(const OpenMM_AmoebaMultipoleForce*& target, OpenMM_IntArray*& gridDimension) {
    OpenMM_AmoebaMultipoleForce_getPmeGridDimensions(target, gridDimension);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_setpmegriddimensions_(OpenMM_AmoebaMultipoleForce*& target, const OpenMM_IntArray*& gridDimension) {
    OpenMM_AmoebaMultipoleForce_setPmeGridDimensions(target, gridDimension);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_SETPMEGRIDDIMENSIONS(OpenMM_AmoebaMultipoleForce*& target, const OpenMM_IntArray*& gridDimension) {
    OpenMM_AmoebaMultipoleForce_setPmeGridDimensions(target, gridDimension);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_getpmeparametersincontext_(const OpenMM_AmoebaMultipoleForce*& target, const OpenMM_Context*& context, double* alpha, int* nx, int* ny, int* nz) {
    OpenMM_AmoebaMultipoleForce_getPMEParametersInContext(target, context, alpha, nx, ny, nz);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_GETPMEPARAMETERSINCONTEXT(const OpenMM_AmoebaMultipoleForce*& target, const OpenMM_Context*& context, double* alpha, int* nx, int* ny, int* nz) {
    OpenMM_AmoebaMultipoleForce_getPMEParametersInContext(target, context, alpha, nx, ny, nz);
}
OPENMM_EXPORT_AMOEBA int openmm_amoebamultipoleforce_addmultipole_(OpenMM_AmoebaMultipoleForce*& target, double const& charge, const OpenMM_DoubleArray*& molecularDipole, const OpenMM_DoubleArray*& molecularQuadrupole, int const& axisType, int const& multipoleAtomZ, int const& multipoleAtomX, int const& multipoleAtomY, double const& thole, double const& dampingFactor, double const& polarity) {
    return OpenMM_AmoebaMultipoleForce_addMultipole(target, charge, molecularDipole, molecularQuadrupole, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY, thole, dampingFactor, polarity);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBAMULTIPOLEFORCE_ADDMULTIPOLE(OpenMM_AmoebaMultipoleForce*& target, double const& charge, const OpenMM_DoubleArray*& molecularDipole, const OpenMM_DoubleArray*& molecularQuadrupole, int const& axisType, int const& multipoleAtomZ, int const& multipoleAtomX, int const& multipoleAtomY, double const& thole, double const& dampingFactor, double const& polarity) {
    return OpenMM_AmoebaMultipoleForce_addMultipole(target, charge, molecularDipole, molecularQuadrupole, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY, thole, dampingFactor, polarity);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_getmultipoleparameters_(const OpenMM_AmoebaMultipoleForce*& target, int const& index, double* charge, OpenMM_DoubleArray*& molecularDipole, OpenMM_DoubleArray*& molecularQuadrupole, int* axisType, int* multipoleAtomZ, int* multipoleAtomX, int* multipoleAtomY, double* thole, double* dampingFactor, double* polarity) {
    OpenMM_AmoebaMultipoleForce_getMultipoleParameters(target, index, charge, molecularDipole, molecularQuadrupole, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY, thole, dampingFactor, polarity);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_GETMULTIPOLEPARAMETERS(const OpenMM_AmoebaMultipoleForce*& target, int const& index, double* charge, OpenMM_DoubleArray*& molecularDipole, OpenMM_DoubleArray*& molecularQuadrupole, int* axisType, int* multipoleAtomZ, int* multipoleAtomX, int* multipoleAtomY, double* thole, double* dampingFactor, double* polarity) {
    OpenMM_AmoebaMultipoleForce_getMultipoleParameters(target, index, charge, molecularDipole, molecularQuadrupole, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY, thole, dampingFactor, polarity);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_setmultipoleparameters_(OpenMM_AmoebaMultipoleForce*& target, int const& index, double const& charge, const OpenMM_DoubleArray*& molecularDipole, const OpenMM_DoubleArray*& molecularQuadrupole, int const& axisType, int const& multipoleAtomZ, int const& multipoleAtomX, int const& multipoleAtomY, double const& thole, double const& dampingFactor, double const& polarity) {
    OpenMM_AmoebaMultipoleForce_setMultipoleParameters(target, index, charge, molecularDipole, molecularQuadrupole, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY, thole, dampingFactor, polarity);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_SETMULTIPOLEPARAMETERS(OpenMM_AmoebaMultipoleForce*& target, int const& index, double const& charge, const OpenMM_DoubleArray*& molecularDipole, const OpenMM_DoubleArray*& molecularQuadrupole, int const& axisType, int const& multipoleAtomZ, int const& multipoleAtomX, int const& multipoleAtomY, double const& thole, double const& dampingFactor, double const& polarity) {
    OpenMM_AmoebaMultipoleForce_setMultipoleParameters(target, index, charge, molecularDipole, molecularQuadrupole, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY, thole, dampingFactor, polarity);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_setcovalentmap_(OpenMM_AmoebaMultipoleForce*& target, int const& index, OpenMM_AmoebaMultipoleForce_CovalentType& typeId, const OpenMM_IntArray*& covalentAtoms) {
    OpenMM_AmoebaMultipoleForce_setCovalentMap(target, index, typeId, covalentAtoms);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_SETCOVALENTMAP(OpenMM_AmoebaMultipoleForce*& target, int const& index, OpenMM_AmoebaMultipoleForce_CovalentType& typeId, const OpenMM_IntArray*& covalentAtoms) {
    OpenMM_AmoebaMultipoleForce_setCovalentMap(target, index, typeId, covalentAtoms);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_getcovalentmap_(const OpenMM_AmoebaMultipoleForce*& target, int const& index, OpenMM_AmoebaMultipoleForce_CovalentType& typeId, OpenMM_IntArray*& covalentAtoms) {
    OpenMM_AmoebaMultipoleForce_getCovalentMap(target, index, typeId, covalentAtoms);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_GETCOVALENTMAP(const OpenMM_AmoebaMultipoleForce*& target, int const& index, OpenMM_AmoebaMultipoleForce_CovalentType& typeId, OpenMM_IntArray*& covalentAtoms) {
    OpenMM_AmoebaMultipoleForce_getCovalentMap(target, index, typeId, covalentAtoms);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_getcovalentmaps_(const OpenMM_AmoebaMultipoleForce*& target, int const& index, OpenMM_2D_IntArray*& covalentLists) {
    OpenMM_AmoebaMultipoleForce_getCovalentMaps(target, index, covalentLists);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_GETCOVALENTMAPS(const OpenMM_AmoebaMultipoleForce*& target, int const& index, OpenMM_2D_IntArray*& covalentLists) {
    OpenMM_AmoebaMultipoleForce_getCovalentMaps(target, index, covalentLists);
}
OPENMM_EXPORT_AMOEBA int openmm_amoebamultipoleforce_getmutualinducedmaxiterations_(const OpenMM_AmoebaMultipoleForce*& target) {
    return OpenMM_AmoebaMultipoleForce_getMutualInducedMaxIterations(target);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBAMULTIPOLEFORCE_GETMUTUALINDUCEDMAXITERATIONS(const OpenMM_AmoebaMultipoleForce*& target) {
    return OpenMM_AmoebaMultipoleForce_getMutualInducedMaxIterations(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_setmutualinducedmaxiterations_(OpenMM_AmoebaMultipoleForce*& target, int const& inputMutualInducedMaxIterations) {
    OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations(target, inputMutualInducedMaxIterations);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_SETMUTUALINDUCEDMAXITERATIONS(OpenMM_AmoebaMultipoleForce*& target, int const& inputMutualInducedMaxIterations) {
    OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations(target, inputMutualInducedMaxIterations);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebamultipoleforce_getmutualinducedtargetepsilon_(const OpenMM_AmoebaMultipoleForce*& target) {
    return OpenMM_AmoebaMultipoleForce_getMutualInducedTargetEpsilon(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAMULTIPOLEFORCE_GETMUTUALINDUCEDTARGETEPSILON(const OpenMM_AmoebaMultipoleForce*& target) {
    return OpenMM_AmoebaMultipoleForce_getMutualInducedTargetEpsilon(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_setmutualinducedtargetepsilon_(OpenMM_AmoebaMultipoleForce*& target, double const& inputMutualInducedTargetEpsilon) {
    OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon(target, inputMutualInducedTargetEpsilon);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_SETMUTUALINDUCEDTARGETEPSILON(OpenMM_AmoebaMultipoleForce*& target, double const& inputMutualInducedTargetEpsilon) {
    OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon(target, inputMutualInducedTargetEpsilon);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_setextrapolationcoefficients_(OpenMM_AmoebaMultipoleForce*& target, const OpenMM_DoubleArray*& coefficients) {
    OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients(target, coefficients);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_SETEXTRAPOLATIONCOEFFICIENTS(OpenMM_AmoebaMultipoleForce*& target, const OpenMM_DoubleArray*& coefficients) {
    OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients(target, coefficients);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_getextrapolationcoefficients_(const OpenMM_AmoebaMultipoleForce*& target, const OpenMM_DoubleArray*& result) {
    result = OpenMM_AmoebaMultipoleForce_getExtrapolationCoefficients(target);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_GETEXTRAPOLATIONCOEFFICIENTS(const OpenMM_AmoebaMultipoleForce*& target, const OpenMM_DoubleArray*& result) {
    result = OpenMM_AmoebaMultipoleForce_getExtrapolationCoefficients(target);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebamultipoleforce_getewalderrortolerance_(const OpenMM_AmoebaMultipoleForce*& target) {
    return OpenMM_AmoebaMultipoleForce_getEwaldErrorTolerance(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAMULTIPOLEFORCE_GETEWALDERRORTOLERANCE(const OpenMM_AmoebaMultipoleForce*& target) {
    return OpenMM_AmoebaMultipoleForce_getEwaldErrorTolerance(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_setewalderrortolerance_(OpenMM_AmoebaMultipoleForce*& target, double const& tol) {
    OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance(target, tol);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_SETEWALDERRORTOLERANCE(OpenMM_AmoebaMultipoleForce*& target, double const& tol) {
    OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance(target, tol);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_getlabframepermanentdipoles_(OpenMM_AmoebaMultipoleForce*& target, OpenMM_Context*& context, OpenMM_Vec3Array*& dipoles) {
    OpenMM_AmoebaMultipoleForce_getLabFramePermanentDipoles(target, context, dipoles);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_GETLABFRAMEPERMANENTDIPOLES(OpenMM_AmoebaMultipoleForce*& target, OpenMM_Context*& context, OpenMM_Vec3Array*& dipoles) {
    OpenMM_AmoebaMultipoleForce_getLabFramePermanentDipoles(target, context, dipoles);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_getinduceddipoles_(OpenMM_AmoebaMultipoleForce*& target, OpenMM_Context*& context, OpenMM_Vec3Array*& dipoles) {
    OpenMM_AmoebaMultipoleForce_getInducedDipoles(target, context, dipoles);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_GETINDUCEDDIPOLES(OpenMM_AmoebaMultipoleForce*& target, OpenMM_Context*& context, OpenMM_Vec3Array*& dipoles) {
    OpenMM_AmoebaMultipoleForce_getInducedDipoles(target, context, dipoles);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_gettotaldipoles_(OpenMM_AmoebaMultipoleForce*& target, OpenMM_Context*& context, OpenMM_Vec3Array*& dipoles) {
    OpenMM_AmoebaMultipoleForce_getTotalDipoles(target, context, dipoles);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_GETTOTALDIPOLES(OpenMM_AmoebaMultipoleForce*& target, OpenMM_Context*& context, OpenMM_Vec3Array*& dipoles) {
    OpenMM_AmoebaMultipoleForce_getTotalDipoles(target, context, dipoles);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_getelectrostaticpotential_(OpenMM_AmoebaMultipoleForce*& target, const OpenMM_Vec3Array*& inputGrid, OpenMM_Context*& context, OpenMM_DoubleArray*& outputElectrostaticPotential) {
    OpenMM_AmoebaMultipoleForce_getElectrostaticPotential(target, inputGrid, context, outputElectrostaticPotential);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_GETELECTROSTATICPOTENTIAL(OpenMM_AmoebaMultipoleForce*& target, const OpenMM_Vec3Array*& inputGrid, OpenMM_Context*& context, OpenMM_DoubleArray*& outputElectrostaticPotential) {
    OpenMM_AmoebaMultipoleForce_getElectrostaticPotential(target, inputGrid, context, outputElectrostaticPotential);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_getsystemmultipolemoments_(OpenMM_AmoebaMultipoleForce*& target, OpenMM_Context*& context, OpenMM_DoubleArray*& outputMultipoleMoments) {
    OpenMM_AmoebaMultipoleForce_getSystemMultipoleMoments(target, context, outputMultipoleMoments);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_GETSYSTEMMULTIPOLEMOMENTS(OpenMM_AmoebaMultipoleForce*& target, OpenMM_Context*& context, OpenMM_DoubleArray*& outputMultipoleMoments) {
    OpenMM_AmoebaMultipoleForce_getSystemMultipoleMoments(target, context, outputMultipoleMoments);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_updateparametersincontext_(OpenMM_AmoebaMultipoleForce*& target, OpenMM_Context*& context) {
    OpenMM_AmoebaMultipoleForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_AmoebaMultipoleForce*& target, OpenMM_Context*& context) {
    OpenMM_AmoebaMultipoleForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebamultipoleforce_usesperiodicboundaryconditions_(const OpenMM_AmoebaMultipoleForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_AmoebaMultipoleForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAMULTIPOLEFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_AmoebaMultipoleForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_AmoebaMultipoleForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::AmoebaStretchBendForce */
OPENMM_EXPORT_AMOEBA void openmm_amoebastretchbendforce_create_(OpenMM_AmoebaStretchBendForce*& result) {
    result = OpenMM_AmoebaStretchBendForce_create();
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBASTRETCHBENDFORCE_CREATE(OpenMM_AmoebaStretchBendForce*& result) {
    result = OpenMM_AmoebaStretchBendForce_create();
}
OPENMM_EXPORT_AMOEBA void openmm_amoebastretchbendforce_destroy_(OpenMM_AmoebaStretchBendForce*& destroy) {
    OpenMM_AmoebaStretchBendForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBASTRETCHBENDFORCE_DESTROY(OpenMM_AmoebaStretchBendForce*& destroy) {
    OpenMM_AmoebaStretchBendForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT_AMOEBA int openmm_amoebastretchbendforce_getnumstretchbends_(const OpenMM_AmoebaStretchBendForce*& target) {
    return OpenMM_AmoebaStretchBendForce_getNumStretchBends(target);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBASTRETCHBENDFORCE_GETNUMSTRETCHBENDS(const OpenMM_AmoebaStretchBendForce*& target) {
    return OpenMM_AmoebaStretchBendForce_getNumStretchBends(target);
}
OPENMM_EXPORT_AMOEBA int openmm_amoebastretchbendforce_addstretchbend_(OpenMM_AmoebaStretchBendForce*& target, int const& particle1, int const& particle2, int const& particle3, double const& lengthAB, double const& lengthCB, double const& angle, double const& k1, double const& k2) {
    return OpenMM_AmoebaStretchBendForce_addStretchBend(target, particle1, particle2, particle3, lengthAB, lengthCB, angle, k1, k2);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBASTRETCHBENDFORCE_ADDSTRETCHBEND(OpenMM_AmoebaStretchBendForce*& target, int const& particle1, int const& particle2, int const& particle3, double const& lengthAB, double const& lengthCB, double const& angle, double const& k1, double const& k2) {
    return OpenMM_AmoebaStretchBendForce_addStretchBend(target, particle1, particle2, particle3, lengthAB, lengthCB, angle, k1, k2);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebastretchbendforce_getstretchbendparameters_(const OpenMM_AmoebaStretchBendForce*& target, int const& index, int* particle1, int* particle2, int* particle3, double* lengthAB, double* lengthCB, double* angle, double* k1, double* k2) {
    OpenMM_AmoebaStretchBendForce_getStretchBendParameters(target, index, particle1, particle2, particle3, lengthAB, lengthCB, angle, k1, k2);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBASTRETCHBENDFORCE_GETSTRETCHBENDPARAMETERS(const OpenMM_AmoebaStretchBendForce*& target, int const& index, int* particle1, int* particle2, int* particle3, double* lengthAB, double* lengthCB, double* angle, double* k1, double* k2) {
    OpenMM_AmoebaStretchBendForce_getStretchBendParameters(target, index, particle1, particle2, particle3, lengthAB, lengthCB, angle, k1, k2);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebastretchbendforce_setstretchbendparameters_(OpenMM_AmoebaStretchBendForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, double const& lengthAB, double const& lengthCB, double const& angle, double const& k1, double const& k2) {
    OpenMM_AmoebaStretchBendForce_setStretchBendParameters(target, index, particle1, particle2, particle3, lengthAB, lengthCB, angle, k1, k2);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBASTRETCHBENDFORCE_SETSTRETCHBENDPARAMETERS(OpenMM_AmoebaStretchBendForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, double const& lengthAB, double const& lengthCB, double const& angle, double const& k1, double const& k2) {
    OpenMM_AmoebaStretchBendForce_setStretchBendParameters(target, index, particle1, particle2, particle3, lengthAB, lengthCB, angle, k1, k2);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebastretchbendforce_updateparametersincontext_(OpenMM_AmoebaStretchBendForce*& target, OpenMM_Context*& context) {
    OpenMM_AmoebaStretchBendForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBASTRETCHBENDFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_AmoebaStretchBendForce*& target, OpenMM_Context*& context) {
    OpenMM_AmoebaStretchBendForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebastretchbendforce_setusesperiodicboundaryconditions_(OpenMM_AmoebaStretchBendForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_AmoebaStretchBendForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBASTRETCHBENDFORCE_SETUSESPERIODICBOUNDARYCONDITIONS(OpenMM_AmoebaStretchBendForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_AmoebaStretchBendForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebastretchbendforce_usesperiodicboundaryconditions_(const OpenMM_AmoebaStretchBendForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_AmoebaStretchBendForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBASTRETCHBENDFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_AmoebaStretchBendForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_AmoebaStretchBendForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::AmoebaVdwForce */
OPENMM_EXPORT_AMOEBA void openmm_amoebavdwforce_create_(OpenMM_AmoebaVdwForce*& result) {
    result = OpenMM_AmoebaVdwForce_create();
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAVDWFORCE_CREATE(OpenMM_AmoebaVdwForce*& result) {
    result = OpenMM_AmoebaVdwForce_create();
}
OPENMM_EXPORT_AMOEBA void openmm_amoebavdwforce_destroy_(OpenMM_AmoebaVdwForce*& destroy) {
    OpenMM_AmoebaVdwForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAVDWFORCE_DESTROY(OpenMM_AmoebaVdwForce*& destroy) {
    OpenMM_AmoebaVdwForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT_AMOEBA void openmm_amoebavdwforce_lambda_(char* result, int result_length) {
    const char* result_chars = OpenMM_AmoebaVdwForce_Lambda();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAVDWFORCE_LAMBDA(char* result, int result_length) {
    const char* result_chars = OpenMM_AmoebaVdwForce_Lambda();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT_AMOEBA int openmm_amoebavdwforce_getnumparticles_(const OpenMM_AmoebaVdwForce*& target) {
    return OpenMM_AmoebaVdwForce_getNumParticles(target);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBAVDWFORCE_GETNUMPARTICLES(const OpenMM_AmoebaVdwForce*& target) {
    return OpenMM_AmoebaVdwForce_getNumParticles(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebavdwforce_setparticleparameters_(OpenMM_AmoebaVdwForce*& target, int const& particleIndex, int const& parentIndex, double const& sigma, double const& epsilon, double const& reductionFactor, OpenMM_Boolean& isAlchemical) {
    OpenMM_AmoebaVdwForce_setParticleParameters(target, particleIndex, parentIndex, sigma, epsilon, reductionFactor, isAlchemical);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAVDWFORCE_SETPARTICLEPARAMETERS(OpenMM_AmoebaVdwForce*& target, int const& particleIndex, int const& parentIndex, double const& sigma, double const& epsilon, double const& reductionFactor, OpenMM_Boolean& isAlchemical) {
    OpenMM_AmoebaVdwForce_setParticleParameters(target, particleIndex, parentIndex, sigma, epsilon, reductionFactor, isAlchemical);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebavdwforce_getparticleparameters_(const OpenMM_AmoebaVdwForce*& target, int const& particleIndex, int* parentIndex, double* sigma, double* epsilon, double* reductionFactor, OpenMM_Boolean*& isAlchemical) {
    OpenMM_AmoebaVdwForce_getParticleParameters(target, particleIndex, parentIndex, sigma, epsilon, reductionFactor, isAlchemical);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAVDWFORCE_GETPARTICLEPARAMETERS(const OpenMM_AmoebaVdwForce*& target, int const& particleIndex, int* parentIndex, double* sigma, double* epsilon, double* reductionFactor, OpenMM_Boolean*& isAlchemical) {
    OpenMM_AmoebaVdwForce_getParticleParameters(target, particleIndex, parentIndex, sigma, epsilon, reductionFactor, isAlchemical);
}
OPENMM_EXPORT_AMOEBA int openmm_amoebavdwforce_addparticle_(OpenMM_AmoebaVdwForce*& target, int const& parentIndex, double const& sigma, double const& epsilon, double const& reductionFactor, OpenMM_Boolean& isAlchemical) {
    return OpenMM_AmoebaVdwForce_addParticle(target, parentIndex, sigma, epsilon, reductionFactor, isAlchemical);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBAVDWFORCE_ADDPARTICLE(OpenMM_AmoebaVdwForce*& target, int const& parentIndex, double const& sigma, double const& epsilon, double const& reductionFactor, OpenMM_Boolean& isAlchemical) {
    return OpenMM_AmoebaVdwForce_addParticle(target, parentIndex, sigma, epsilon, reductionFactor, isAlchemical);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebavdwforce_setsigmacombiningrule_(OpenMM_AmoebaVdwForce*& target, const char* sigmaCombiningRule, int sigmaCombiningRule_length) {
    OpenMM_AmoebaVdwForce_setSigmaCombiningRule(target, makeString(sigmaCombiningRule, sigmaCombiningRule_length).c_str());
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAVDWFORCE_SETSIGMACOMBININGRULE(OpenMM_AmoebaVdwForce*& target, const char* sigmaCombiningRule, int sigmaCombiningRule_length) {
    OpenMM_AmoebaVdwForce_setSigmaCombiningRule(target, makeString(sigmaCombiningRule, sigmaCombiningRule_length).c_str());
}
OPENMM_EXPORT_AMOEBA void openmm_amoebavdwforce_getsigmacombiningrule_(const OpenMM_AmoebaVdwForce*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_AmoebaVdwForce_getSigmaCombiningRule(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAVDWFORCE_GETSIGMACOMBININGRULE(const OpenMM_AmoebaVdwForce*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_AmoebaVdwForce_getSigmaCombiningRule(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebavdwforce_setepsiloncombiningrule_(OpenMM_AmoebaVdwForce*& target, const char* epsilonCombiningRule, int epsilonCombiningRule_length) {
    OpenMM_AmoebaVdwForce_setEpsilonCombiningRule(target, makeString(epsilonCombiningRule, epsilonCombiningRule_length).c_str());
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAVDWFORCE_SETEPSILONCOMBININGRULE(OpenMM_AmoebaVdwForce*& target, const char* epsilonCombiningRule, int epsilonCombiningRule_length) {
    OpenMM_AmoebaVdwForce_setEpsilonCombiningRule(target, makeString(epsilonCombiningRule, epsilonCombiningRule_length).c_str());
}
OPENMM_EXPORT_AMOEBA void openmm_amoebavdwforce_getepsiloncombiningrule_(const OpenMM_AmoebaVdwForce*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_AmoebaVdwForce_getEpsilonCombiningRule(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAVDWFORCE_GETEPSILONCOMBININGRULE(const OpenMM_AmoebaVdwForce*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_AmoebaVdwForce_getEpsilonCombiningRule(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebavdwforce_getusedispersioncorrection_(const OpenMM_AmoebaVdwForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_AmoebaVdwForce_getUseDispersionCorrection(target);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAVDWFORCE_GETUSEDISPERSIONCORRECTION(const OpenMM_AmoebaVdwForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_AmoebaVdwForce_getUseDispersionCorrection(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebavdwforce_setusedispersioncorrection_(OpenMM_AmoebaVdwForce*& target, OpenMM_Boolean& useCorrection) {
    OpenMM_AmoebaVdwForce_setUseDispersionCorrection(target, useCorrection);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAVDWFORCE_SETUSEDISPERSIONCORRECTION(OpenMM_AmoebaVdwForce*& target, OpenMM_Boolean& useCorrection) {
    OpenMM_AmoebaVdwForce_setUseDispersionCorrection(target, useCorrection);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebavdwforce_setparticleexclusions_(OpenMM_AmoebaVdwForce*& target, int const& particleIndex, const OpenMM_IntArray*& exclusions) {
    OpenMM_AmoebaVdwForce_setParticleExclusions(target, particleIndex, exclusions);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAVDWFORCE_SETPARTICLEEXCLUSIONS(OpenMM_AmoebaVdwForce*& target, int const& particleIndex, const OpenMM_IntArray*& exclusions) {
    OpenMM_AmoebaVdwForce_setParticleExclusions(target, particleIndex, exclusions);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebavdwforce_getparticleexclusions_(const OpenMM_AmoebaVdwForce*& target, int const& particleIndex, OpenMM_IntArray*& exclusions) {
    OpenMM_AmoebaVdwForce_getParticleExclusions(target, particleIndex, exclusions);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAVDWFORCE_GETPARTICLEEXCLUSIONS(const OpenMM_AmoebaVdwForce*& target, int const& particleIndex, OpenMM_IntArray*& exclusions) {
    OpenMM_AmoebaVdwForce_getParticleExclusions(target, particleIndex, exclusions);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebavdwforce_getcutoffdistance_(const OpenMM_AmoebaVdwForce*& target) {
    return OpenMM_AmoebaVdwForce_getCutoffDistance(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAVDWFORCE_GETCUTOFFDISTANCE(const OpenMM_AmoebaVdwForce*& target) {
    return OpenMM_AmoebaVdwForce_getCutoffDistance(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebavdwforce_setcutoffdistance_(OpenMM_AmoebaVdwForce*& target, double const& distance) {
    OpenMM_AmoebaVdwForce_setCutoffDistance(target, distance);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAVDWFORCE_SETCUTOFFDISTANCE(OpenMM_AmoebaVdwForce*& target, double const& distance) {
    OpenMM_AmoebaVdwForce_setCutoffDistance(target, distance);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebavdwforce_setcutoff_(OpenMM_AmoebaVdwForce*& target, double const& cutoff) {
    OpenMM_AmoebaVdwForce_setCutoff(target, cutoff);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAVDWFORCE_SETCUTOFF(OpenMM_AmoebaVdwForce*& target, double const& cutoff) {
    OpenMM_AmoebaVdwForce_setCutoff(target, cutoff);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebavdwforce_getcutoff_(const OpenMM_AmoebaVdwForce*& target) {
    return OpenMM_AmoebaVdwForce_getCutoff(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAVDWFORCE_GETCUTOFF(const OpenMM_AmoebaVdwForce*& target) {
    return OpenMM_AmoebaVdwForce_getCutoff(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebavdwforce_getnonbondedmethod_(const OpenMM_AmoebaVdwForce*& target, OpenMM_AmoebaVdwForce_NonbondedMethod& result) {
    result = OpenMM_AmoebaVdwForce_getNonbondedMethod(target);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAVDWFORCE_GETNONBONDEDMETHOD(const OpenMM_AmoebaVdwForce*& target, OpenMM_AmoebaVdwForce_NonbondedMethod& result) {
    result = OpenMM_AmoebaVdwForce_getNonbondedMethod(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebavdwforce_setnonbondedmethod_(OpenMM_AmoebaVdwForce*& target, OpenMM_AmoebaVdwForce_NonbondedMethod& method) {
    OpenMM_AmoebaVdwForce_setNonbondedMethod(target, method);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAVDWFORCE_SETNONBONDEDMETHOD(OpenMM_AmoebaVdwForce*& target, OpenMM_AmoebaVdwForce_NonbondedMethod& method) {
    OpenMM_AmoebaVdwForce_setNonbondedMethod(target, method);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebavdwforce_setsoftcorepower_(OpenMM_AmoebaVdwForce*& target, int const& n) {
    OpenMM_AmoebaVdwForce_setSoftcorePower(target, n);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAVDWFORCE_SETSOFTCOREPOWER(OpenMM_AmoebaVdwForce*& target, int const& n) {
    OpenMM_AmoebaVdwForce_setSoftcorePower(target, n);
}
OPENMM_EXPORT_AMOEBA int openmm_amoebavdwforce_getsoftcorepower_(const OpenMM_AmoebaVdwForce*& target) {
    return OpenMM_AmoebaVdwForce_getSoftcorePower(target);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBAVDWFORCE_GETSOFTCOREPOWER(const OpenMM_AmoebaVdwForce*& target) {
    return OpenMM_AmoebaVdwForce_getSoftcorePower(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebavdwforce_setsoftcorealpha_(OpenMM_AmoebaVdwForce*& target, double const& alpha) {
    OpenMM_AmoebaVdwForce_setSoftcoreAlpha(target, alpha);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAVDWFORCE_SETSOFTCOREALPHA(OpenMM_AmoebaVdwForce*& target, double const& alpha) {
    OpenMM_AmoebaVdwForce_setSoftcoreAlpha(target, alpha);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebavdwforce_getsoftcorealpha_(const OpenMM_AmoebaVdwForce*& target) {
    return OpenMM_AmoebaVdwForce_getSoftcoreAlpha(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAVDWFORCE_GETSOFTCOREALPHA(const OpenMM_AmoebaVdwForce*& target) {
    return OpenMM_AmoebaVdwForce_getSoftcoreAlpha(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebavdwforce_getalchemicalmethod_(const OpenMM_AmoebaVdwForce*& target, OpenMM_AmoebaVdwForce_AlchemicalMethod& result) {
    result = OpenMM_AmoebaVdwForce_getAlchemicalMethod(target);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAVDWFORCE_GETALCHEMICALMETHOD(const OpenMM_AmoebaVdwForce*& target, OpenMM_AmoebaVdwForce_AlchemicalMethod& result) {
    result = OpenMM_AmoebaVdwForce_getAlchemicalMethod(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebavdwforce_setalchemicalmethod_(OpenMM_AmoebaVdwForce*& target, OpenMM_AmoebaVdwForce_AlchemicalMethod& method) {
    OpenMM_AmoebaVdwForce_setAlchemicalMethod(target, method);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAVDWFORCE_SETALCHEMICALMETHOD(OpenMM_AmoebaVdwForce*& target, OpenMM_AmoebaVdwForce_AlchemicalMethod& method) {
    OpenMM_AmoebaVdwForce_setAlchemicalMethod(target, method);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebavdwforce_updateparametersincontext_(OpenMM_AmoebaVdwForce*& target, OpenMM_Context*& context) {
    OpenMM_AmoebaVdwForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAVDWFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_AmoebaVdwForce*& target, OpenMM_Context*& context) {
    OpenMM_AmoebaVdwForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebavdwforce_usesperiodicboundaryconditions_(const OpenMM_AmoebaVdwForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_AmoebaVdwForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAVDWFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_AmoebaVdwForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_AmoebaVdwForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::AmoebaGeneralizedKirkwoodForce */
OPENMM_EXPORT_AMOEBA void openmm_amoebageneralizedkirkwoodforce_create_(OpenMM_AmoebaGeneralizedKirkwoodForce*& result) {
    result = OpenMM_AmoebaGeneralizedKirkwoodForce_create();
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAGENERALIZEDKIRKWOODFORCE_CREATE(OpenMM_AmoebaGeneralizedKirkwoodForce*& result) {
    result = OpenMM_AmoebaGeneralizedKirkwoodForce_create();
}
OPENMM_EXPORT_AMOEBA void openmm_amoebageneralizedkirkwoodforce_destroy_(OpenMM_AmoebaGeneralizedKirkwoodForce*& destroy) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAGENERALIZEDKIRKWOODFORCE_DESTROY(OpenMM_AmoebaGeneralizedKirkwoodForce*& destroy) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT_AMOEBA int openmm_amoebageneralizedkirkwoodforce_getnumparticles_(const OpenMM_AmoebaGeneralizedKirkwoodForce*& target) {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getNumParticles(target);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBAGENERALIZEDKIRKWOODFORCE_GETNUMPARTICLES(const OpenMM_AmoebaGeneralizedKirkwoodForce*& target) {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getNumParticles(target);
}
OPENMM_EXPORT_AMOEBA int openmm_amoebageneralizedkirkwoodforce_addparticle_(OpenMM_AmoebaGeneralizedKirkwoodForce*& target, double const& charge, double const& radius, double const& scalingFactor) {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle(target, charge, radius, scalingFactor);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBAGENERALIZEDKIRKWOODFORCE_ADDPARTICLE(OpenMM_AmoebaGeneralizedKirkwoodForce*& target, double const& charge, double const& radius, double const& scalingFactor) {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle(target, charge, radius, scalingFactor);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebageneralizedkirkwoodforce_getparticleparameters_(const OpenMM_AmoebaGeneralizedKirkwoodForce*& target, int const& index, double* charge, double* radius, double* scalingFactor) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_getParticleParameters(target, index, charge, radius, scalingFactor);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAGENERALIZEDKIRKWOODFORCE_GETPARTICLEPARAMETERS(const OpenMM_AmoebaGeneralizedKirkwoodForce*& target, int const& index, double* charge, double* radius, double* scalingFactor) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_getParticleParameters(target, index, charge, radius, scalingFactor);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebageneralizedkirkwoodforce_setparticleparameters_(OpenMM_AmoebaGeneralizedKirkwoodForce*& target, int const& index, double const& charge, double const& radius, double const& scalingFactor) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters(target, index, charge, radius, scalingFactor);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAGENERALIZEDKIRKWOODFORCE_SETPARTICLEPARAMETERS(OpenMM_AmoebaGeneralizedKirkwoodForce*& target, int const& index, double const& charge, double const& radius, double const& scalingFactor) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters(target, index, charge, radius, scalingFactor);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebageneralizedkirkwoodforce_getsolventdielectric_(const OpenMM_AmoebaGeneralizedKirkwoodForce*& target) {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getSolventDielectric(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAGENERALIZEDKIRKWOODFORCE_GETSOLVENTDIELECTRIC(const OpenMM_AmoebaGeneralizedKirkwoodForce*& target) {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getSolventDielectric(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebageneralizedkirkwoodforce_setsolventdielectric_(OpenMM_AmoebaGeneralizedKirkwoodForce*& target, double const& dielectric) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric(target, dielectric);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAGENERALIZEDKIRKWOODFORCE_SETSOLVENTDIELECTRIC(OpenMM_AmoebaGeneralizedKirkwoodForce*& target, double const& dielectric) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric(target, dielectric);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebageneralizedkirkwoodforce_getsolutedielectric_(const OpenMM_AmoebaGeneralizedKirkwoodForce*& target) {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getSoluteDielectric(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAGENERALIZEDKIRKWOODFORCE_GETSOLUTEDIELECTRIC(const OpenMM_AmoebaGeneralizedKirkwoodForce*& target) {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getSoluteDielectric(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebageneralizedkirkwoodforce_setsolutedielectric_(OpenMM_AmoebaGeneralizedKirkwoodForce*& target, double const& dielectric) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric(target, dielectric);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAGENERALIZEDKIRKWOODFORCE_SETSOLUTEDIELECTRIC(OpenMM_AmoebaGeneralizedKirkwoodForce*& target, double const& dielectric) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric(target, dielectric);
}
OPENMM_EXPORT_AMOEBA int openmm_amoebageneralizedkirkwoodforce_getincludecavityterm_(const OpenMM_AmoebaGeneralizedKirkwoodForce*& target) {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getIncludeCavityTerm(target);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBAGENERALIZEDKIRKWOODFORCE_GETINCLUDECAVITYTERM(const OpenMM_AmoebaGeneralizedKirkwoodForce*& target) {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getIncludeCavityTerm(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebageneralizedkirkwoodforce_setincludecavityterm_(OpenMM_AmoebaGeneralizedKirkwoodForce*& target, int const& includeCavityTerm) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm(target, includeCavityTerm);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAGENERALIZEDKIRKWOODFORCE_SETINCLUDECAVITYTERM(OpenMM_AmoebaGeneralizedKirkwoodForce*& target, int const& includeCavityTerm) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm(target, includeCavityTerm);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebageneralizedkirkwoodforce_getproberadius_(const OpenMM_AmoebaGeneralizedKirkwoodForce*& target) {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getProbeRadius(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAGENERALIZEDKIRKWOODFORCE_GETPROBERADIUS(const OpenMM_AmoebaGeneralizedKirkwoodForce*& target) {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getProbeRadius(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebageneralizedkirkwoodforce_setproberadius_(OpenMM_AmoebaGeneralizedKirkwoodForce*& target, double const& probeRadius) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius(target, probeRadius);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAGENERALIZEDKIRKWOODFORCE_SETPROBERADIUS(OpenMM_AmoebaGeneralizedKirkwoodForce*& target, double const& probeRadius) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius(target, probeRadius);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebageneralizedkirkwoodforce_getsurfaceareafactor_(const OpenMM_AmoebaGeneralizedKirkwoodForce*& target) {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getSurfaceAreaFactor(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAGENERALIZEDKIRKWOODFORCE_GETSURFACEAREAFACTOR(const OpenMM_AmoebaGeneralizedKirkwoodForce*& target) {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getSurfaceAreaFactor(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebageneralizedkirkwoodforce_setsurfaceareafactor_(OpenMM_AmoebaGeneralizedKirkwoodForce*& target, double const& surfaceAreaFactor) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor(target, surfaceAreaFactor);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAGENERALIZEDKIRKWOODFORCE_SETSURFACEAREAFACTOR(OpenMM_AmoebaGeneralizedKirkwoodForce*& target, double const& surfaceAreaFactor) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor(target, surfaceAreaFactor);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebageneralizedkirkwoodforce_updateparametersincontext_(OpenMM_AmoebaGeneralizedKirkwoodForce*& target, OpenMM_Context*& context) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAGENERALIZEDKIRKWOODFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_AmoebaGeneralizedKirkwoodForce*& target, OpenMM_Context*& context) {
    OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebageneralizedkirkwoodforce_usesperiodicboundaryconditions_(const OpenMM_AmoebaGeneralizedKirkwoodForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_AmoebaGeneralizedKirkwoodForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAGENERALIZEDKIRKWOODFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_AmoebaGeneralizedKirkwoodForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_AmoebaGeneralizedKirkwoodForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::AmoebaOutOfPlaneBendForce */
OPENMM_EXPORT_AMOEBA void openmm_amoebaoutofplanebendforce_create_(OpenMM_AmoebaOutOfPlaneBendForce*& result) {
    result = OpenMM_AmoebaOutOfPlaneBendForce_create();
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAOUTOFPLANEBENDFORCE_CREATE(OpenMM_AmoebaOutOfPlaneBendForce*& result) {
    result = OpenMM_AmoebaOutOfPlaneBendForce_create();
}
OPENMM_EXPORT_AMOEBA void openmm_amoebaoutofplanebendforce_destroy_(OpenMM_AmoebaOutOfPlaneBendForce*& destroy) {
    OpenMM_AmoebaOutOfPlaneBendForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAOUTOFPLANEBENDFORCE_DESTROY(OpenMM_AmoebaOutOfPlaneBendForce*& destroy) {
    OpenMM_AmoebaOutOfPlaneBendForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT_AMOEBA int openmm_amoebaoutofplanebendforce_getnumoutofplanebends_(const OpenMM_AmoebaOutOfPlaneBendForce*& target) {
    return OpenMM_AmoebaOutOfPlaneBendForce_getNumOutOfPlaneBends(target);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBAOUTOFPLANEBENDFORCE_GETNUMOUTOFPLANEBENDS(const OpenMM_AmoebaOutOfPlaneBendForce*& target) {
    return OpenMM_AmoebaOutOfPlaneBendForce_getNumOutOfPlaneBends(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebaoutofplanebendforce_setamoebaglobaloutofplanebendcubic_(OpenMM_AmoebaOutOfPlaneBendForce*& target, double const& cubicK) {
    OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendCubic(target, cubicK);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAOUTOFPLANEBENDFORCE_SETAMOEBAGLOBALOUTOFPLANEBENDCUBIC(OpenMM_AmoebaOutOfPlaneBendForce*& target, double const& cubicK) {
    OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendCubic(target, cubicK);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebaoutofplanebendforce_getamoebaglobaloutofplanebendcubic_(const OpenMM_AmoebaOutOfPlaneBendForce*& target) {
    return OpenMM_AmoebaOutOfPlaneBendForce_getAmoebaGlobalOutOfPlaneBendCubic(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAOUTOFPLANEBENDFORCE_GETAMOEBAGLOBALOUTOFPLANEBENDCUBIC(const OpenMM_AmoebaOutOfPlaneBendForce*& target) {
    return OpenMM_AmoebaOutOfPlaneBendForce_getAmoebaGlobalOutOfPlaneBendCubic(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebaoutofplanebendforce_setamoebaglobaloutofplanebendquartic_(OpenMM_AmoebaOutOfPlaneBendForce*& target, double const& quarticK) {
    OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendQuartic(target, quarticK);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAOUTOFPLANEBENDFORCE_SETAMOEBAGLOBALOUTOFPLANEBENDQUARTIC(OpenMM_AmoebaOutOfPlaneBendForce*& target, double const& quarticK) {
    OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendQuartic(target, quarticK);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebaoutofplanebendforce_getamoebaglobaloutofplanebendquartic_(const OpenMM_AmoebaOutOfPlaneBendForce*& target) {
    return OpenMM_AmoebaOutOfPlaneBendForce_getAmoebaGlobalOutOfPlaneBendQuartic(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAOUTOFPLANEBENDFORCE_GETAMOEBAGLOBALOUTOFPLANEBENDQUARTIC(const OpenMM_AmoebaOutOfPlaneBendForce*& target) {
    return OpenMM_AmoebaOutOfPlaneBendForce_getAmoebaGlobalOutOfPlaneBendQuartic(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebaoutofplanebendforce_setamoebaglobaloutofplanebendpentic_(OpenMM_AmoebaOutOfPlaneBendForce*& target, double const& penticK) {
    OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendPentic(target, penticK);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAOUTOFPLANEBENDFORCE_SETAMOEBAGLOBALOUTOFPLANEBENDPENTIC(OpenMM_AmoebaOutOfPlaneBendForce*& target, double const& penticK) {
    OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendPentic(target, penticK);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebaoutofplanebendforce_getamoebaglobaloutofplanebendpentic_(const OpenMM_AmoebaOutOfPlaneBendForce*& target) {
    return OpenMM_AmoebaOutOfPlaneBendForce_getAmoebaGlobalOutOfPlaneBendPentic(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAOUTOFPLANEBENDFORCE_GETAMOEBAGLOBALOUTOFPLANEBENDPENTIC(const OpenMM_AmoebaOutOfPlaneBendForce*& target) {
    return OpenMM_AmoebaOutOfPlaneBendForce_getAmoebaGlobalOutOfPlaneBendPentic(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebaoutofplanebendforce_setamoebaglobaloutofplanebendsextic_(OpenMM_AmoebaOutOfPlaneBendForce*& target, double const& sexticK) {
    OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendSextic(target, sexticK);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAOUTOFPLANEBENDFORCE_SETAMOEBAGLOBALOUTOFPLANEBENDSEXTIC(OpenMM_AmoebaOutOfPlaneBendForce*& target, double const& sexticK) {
    OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendSextic(target, sexticK);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebaoutofplanebendforce_getamoebaglobaloutofplanebendsextic_(const OpenMM_AmoebaOutOfPlaneBendForce*& target) {
    return OpenMM_AmoebaOutOfPlaneBendForce_getAmoebaGlobalOutOfPlaneBendSextic(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAOUTOFPLANEBENDFORCE_GETAMOEBAGLOBALOUTOFPLANEBENDSEXTIC(const OpenMM_AmoebaOutOfPlaneBendForce*& target) {
    return OpenMM_AmoebaOutOfPlaneBendForce_getAmoebaGlobalOutOfPlaneBendSextic(target);
}
OPENMM_EXPORT_AMOEBA int openmm_amoebaoutofplanebendforce_addoutofplanebend_(OpenMM_AmoebaOutOfPlaneBendForce*& target, int const& particle1, int const& particle2, int const& particle3, int const& particle4, double const& k) {
    return OpenMM_AmoebaOutOfPlaneBendForce_addOutOfPlaneBend(target, particle1, particle2, particle3, particle4, k);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBAOUTOFPLANEBENDFORCE_ADDOUTOFPLANEBEND(OpenMM_AmoebaOutOfPlaneBendForce*& target, int const& particle1, int const& particle2, int const& particle3, int const& particle4, double const& k) {
    return OpenMM_AmoebaOutOfPlaneBendForce_addOutOfPlaneBend(target, particle1, particle2, particle3, particle4, k);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebaoutofplanebendforce_getoutofplanebendparameters_(const OpenMM_AmoebaOutOfPlaneBendForce*& target, int const& index, int* particle1, int* particle2, int* particle3, int* particle4, double* k) {
    OpenMM_AmoebaOutOfPlaneBendForce_getOutOfPlaneBendParameters(target, index, particle1, particle2, particle3, particle4, k);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAOUTOFPLANEBENDFORCE_GETOUTOFPLANEBENDPARAMETERS(const OpenMM_AmoebaOutOfPlaneBendForce*& target, int const& index, int* particle1, int* particle2, int* particle3, int* particle4, double* k) {
    OpenMM_AmoebaOutOfPlaneBendForce_getOutOfPlaneBendParameters(target, index, particle1, particle2, particle3, particle4, k);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebaoutofplanebendforce_setoutofplanebendparameters_(OpenMM_AmoebaOutOfPlaneBendForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, int const& particle4, double const& k) {
    OpenMM_AmoebaOutOfPlaneBendForce_setOutOfPlaneBendParameters(target, index, particle1, particle2, particle3, particle4, k);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAOUTOFPLANEBENDFORCE_SETOUTOFPLANEBENDPARAMETERS(OpenMM_AmoebaOutOfPlaneBendForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, int const& particle4, double const& k) {
    OpenMM_AmoebaOutOfPlaneBendForce_setOutOfPlaneBendParameters(target, index, particle1, particle2, particle3, particle4, k);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebaoutofplanebendforce_updateparametersincontext_(OpenMM_AmoebaOutOfPlaneBendForce*& target, OpenMM_Context*& context) {
    OpenMM_AmoebaOutOfPlaneBendForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAOUTOFPLANEBENDFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_AmoebaOutOfPlaneBendForce*& target, OpenMM_Context*& context) {
    OpenMM_AmoebaOutOfPlaneBendForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebaoutofplanebendforce_setusesperiodicboundaryconditions_(OpenMM_AmoebaOutOfPlaneBendForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_AmoebaOutOfPlaneBendForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAOUTOFPLANEBENDFORCE_SETUSESPERIODICBOUNDARYCONDITIONS(OpenMM_AmoebaOutOfPlaneBendForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_AmoebaOutOfPlaneBendForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebaoutofplanebendforce_usesperiodicboundaryconditions_(const OpenMM_AmoebaOutOfPlaneBendForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_AmoebaOutOfPlaneBendForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAOUTOFPLANEBENDFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_AmoebaOutOfPlaneBendForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_AmoebaOutOfPlaneBendForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::HippoNonbondedForce */
OPENMM_EXPORT_AMOEBA void openmm_hippononbondedforce_create_(OpenMM_HippoNonbondedForce*& result) {
    result = OpenMM_HippoNonbondedForce_create();
}
OPENMM_EXPORT_AMOEBA void OPENMM_HIPPONONBONDEDFORCE_CREATE(OpenMM_HippoNonbondedForce*& result) {
    result = OpenMM_HippoNonbondedForce_create();
}
OPENMM_EXPORT_AMOEBA void openmm_hippononbondedforce_destroy_(OpenMM_HippoNonbondedForce*& destroy) {
    OpenMM_HippoNonbondedForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT_AMOEBA void OPENMM_HIPPONONBONDEDFORCE_DESTROY(OpenMM_HippoNonbondedForce*& destroy) {
    OpenMM_HippoNonbondedForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT_AMOEBA int openmm_hippononbondedforce_getnumparticles_(const OpenMM_HippoNonbondedForce*& target) {
    return OpenMM_HippoNonbondedForce_getNumParticles(target);
}
OPENMM_EXPORT_AMOEBA int OPENMM_HIPPONONBONDEDFORCE_GETNUMPARTICLES(const OpenMM_HippoNonbondedForce*& target) {
    return OpenMM_HippoNonbondedForce_getNumParticles(target);
}
OPENMM_EXPORT_AMOEBA int openmm_hippononbondedforce_getnumexceptions_(const OpenMM_HippoNonbondedForce*& target) {
    return OpenMM_HippoNonbondedForce_getNumExceptions(target);
}
OPENMM_EXPORT_AMOEBA int OPENMM_HIPPONONBONDEDFORCE_GETNUMEXCEPTIONS(const OpenMM_HippoNonbondedForce*& target) {
    return OpenMM_HippoNonbondedForce_getNumExceptions(target);
}
OPENMM_EXPORT_AMOEBA void openmm_hippononbondedforce_getnonbondedmethod_(const OpenMM_HippoNonbondedForce*& target, OpenMM_HippoNonbondedForce_NonbondedMethod& result) {
    result = OpenMM_HippoNonbondedForce_getNonbondedMethod(target);
}
OPENMM_EXPORT_AMOEBA void OPENMM_HIPPONONBONDEDFORCE_GETNONBONDEDMETHOD(const OpenMM_HippoNonbondedForce*& target, OpenMM_HippoNonbondedForce_NonbondedMethod& result) {
    result = OpenMM_HippoNonbondedForce_getNonbondedMethod(target);
}
OPENMM_EXPORT_AMOEBA void openmm_hippononbondedforce_setnonbondedmethod_(OpenMM_HippoNonbondedForce*& target, OpenMM_HippoNonbondedForce_NonbondedMethod& method) {
    OpenMM_HippoNonbondedForce_setNonbondedMethod(target, method);
}
OPENMM_EXPORT_AMOEBA void OPENMM_HIPPONONBONDEDFORCE_SETNONBONDEDMETHOD(OpenMM_HippoNonbondedForce*& target, OpenMM_HippoNonbondedForce_NonbondedMethod& method) {
    OpenMM_HippoNonbondedForce_setNonbondedMethod(target, method);
}
OPENMM_EXPORT_AMOEBA double openmm_hippononbondedforce_getcutoffdistance_(const OpenMM_HippoNonbondedForce*& target) {
    return OpenMM_HippoNonbondedForce_getCutoffDistance(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_HIPPONONBONDEDFORCE_GETCUTOFFDISTANCE(const OpenMM_HippoNonbondedForce*& target) {
    return OpenMM_HippoNonbondedForce_getCutoffDistance(target);
}
OPENMM_EXPORT_AMOEBA void openmm_hippononbondedforce_setcutoffdistance_(OpenMM_HippoNonbondedForce*& target, double const& distance) {
    OpenMM_HippoNonbondedForce_setCutoffDistance(target, distance);
}
OPENMM_EXPORT_AMOEBA void OPENMM_HIPPONONBONDEDFORCE_SETCUTOFFDISTANCE(OpenMM_HippoNonbondedForce*& target, double const& distance) {
    OpenMM_HippoNonbondedForce_setCutoffDistance(target, distance);
}
OPENMM_EXPORT_AMOEBA double openmm_hippononbondedforce_getswitchingdistance_(const OpenMM_HippoNonbondedForce*& target) {
    return OpenMM_HippoNonbondedForce_getSwitchingDistance(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_HIPPONONBONDEDFORCE_GETSWITCHINGDISTANCE(const OpenMM_HippoNonbondedForce*& target) {
    return OpenMM_HippoNonbondedForce_getSwitchingDistance(target);
}
OPENMM_EXPORT_AMOEBA void openmm_hippononbondedforce_setswitchingdistance_(OpenMM_HippoNonbondedForce*& target, double const& distance) {
    OpenMM_HippoNonbondedForce_setSwitchingDistance(target, distance);
}
OPENMM_EXPORT_AMOEBA void OPENMM_HIPPONONBONDEDFORCE_SETSWITCHINGDISTANCE(OpenMM_HippoNonbondedForce*& target, double const& distance) {
    OpenMM_HippoNonbondedForce_setSwitchingDistance(target, distance);
}
OPENMM_EXPORT_AMOEBA void openmm_hippononbondedforce_getextrapolationcoefficients_(const OpenMM_HippoNonbondedForce*& target, const OpenMM_DoubleArray*& result) {
    result = OpenMM_HippoNonbondedForce_getExtrapolationCoefficients(target);
}
OPENMM_EXPORT_AMOEBA void OPENMM_HIPPONONBONDEDFORCE_GETEXTRAPOLATIONCOEFFICIENTS(const OpenMM_HippoNonbondedForce*& target, const OpenMM_DoubleArray*& result) {
    result = OpenMM_HippoNonbondedForce_getExtrapolationCoefficients(target);
}
OPENMM_EXPORT_AMOEBA void openmm_hippononbondedforce_setextrapolationcoefficients_(OpenMM_HippoNonbondedForce*& target, const OpenMM_DoubleArray*& coefficients) {
    OpenMM_HippoNonbondedForce_setExtrapolationCoefficients(target, coefficients);
}
OPENMM_EXPORT_AMOEBA void OPENMM_HIPPONONBONDEDFORCE_SETEXTRAPOLATIONCOEFFICIENTS(OpenMM_HippoNonbondedForce*& target, const OpenMM_DoubleArray*& coefficients) {
    OpenMM_HippoNonbondedForce_setExtrapolationCoefficients(target, coefficients);
}
OPENMM_EXPORT_AMOEBA void openmm_hippononbondedforce_getpmeparameters_(const OpenMM_HippoNonbondedForce*& target, double* alpha, int* nx, int* ny, int* nz) {
    OpenMM_HippoNonbondedForce_getPMEParameters(target, alpha, nx, ny, nz);
}
OPENMM_EXPORT_AMOEBA void OPENMM_HIPPONONBONDEDFORCE_GETPMEPARAMETERS(const OpenMM_HippoNonbondedForce*& target, double* alpha, int* nx, int* ny, int* nz) {
    OpenMM_HippoNonbondedForce_getPMEParameters(target, alpha, nx, ny, nz);
}
OPENMM_EXPORT_AMOEBA void openmm_hippononbondedforce_getdpmeparameters_(const OpenMM_HippoNonbondedForce*& target, double* alpha, int* nx, int* ny, int* nz) {
    OpenMM_HippoNonbondedForce_getDPMEParameters(target, alpha, nx, ny, nz);
}
OPENMM_EXPORT_AMOEBA void OPENMM_HIPPONONBONDEDFORCE_GETDPMEPARAMETERS(const OpenMM_HippoNonbondedForce*& target, double* alpha, int* nx, int* ny, int* nz) {
    OpenMM_HippoNonbondedForce_getDPMEParameters(target, alpha, nx, ny, nz);
}
OPENMM_EXPORT_AMOEBA void openmm_hippononbondedforce_setpmeparameters_(OpenMM_HippoNonbondedForce*& target, double const& alpha, int const& nx, int const& ny, int const& nz) {
    OpenMM_HippoNonbondedForce_setPMEParameters(target, alpha, nx, ny, nz);
}
OPENMM_EXPORT_AMOEBA void OPENMM_HIPPONONBONDEDFORCE_SETPMEPARAMETERS(OpenMM_HippoNonbondedForce*& target, double const& alpha, int const& nx, int const& ny, int const& nz) {
    OpenMM_HippoNonbondedForce_setPMEParameters(target, alpha, nx, ny, nz);
}
OPENMM_EXPORT_AMOEBA void openmm_hippononbondedforce_setdpmeparameters_(OpenMM_HippoNonbondedForce*& target, double const& alpha, int const& nx, int const& ny, int const& nz) {
    OpenMM_HippoNonbondedForce_setDPMEParameters(target, alpha, nx, ny, nz);
}
OPENMM_EXPORT_AMOEBA void OPENMM_HIPPONONBONDEDFORCE_SETDPMEPARAMETERS(OpenMM_HippoNonbondedForce*& target, double const& alpha, int const& nx, int const& ny, int const& nz) {
    OpenMM_HippoNonbondedForce_setDPMEParameters(target, alpha, nx, ny, nz);
}
OPENMM_EXPORT_AMOEBA void openmm_hippononbondedforce_getpmeparametersincontext_(const OpenMM_HippoNonbondedForce*& target, const OpenMM_Context*& context, double* alpha, int* nx, int* ny, int* nz) {
    OpenMM_HippoNonbondedForce_getPMEParametersInContext(target, context, alpha, nx, ny, nz);
}
OPENMM_EXPORT_AMOEBA void OPENMM_HIPPONONBONDEDFORCE_GETPMEPARAMETERSINCONTEXT(const OpenMM_HippoNonbondedForce*& target, const OpenMM_Context*& context, double* alpha, int* nx, int* ny, int* nz) {
    OpenMM_HippoNonbondedForce_getPMEParametersInContext(target, context, alpha, nx, ny, nz);
}
OPENMM_EXPORT_AMOEBA void openmm_hippononbondedforce_getdpmeparametersincontext_(const OpenMM_HippoNonbondedForce*& target, const OpenMM_Context*& context, double* alpha, int* nx, int* ny, int* nz) {
    OpenMM_HippoNonbondedForce_getDPMEParametersInContext(target, context, alpha, nx, ny, nz);
}
OPENMM_EXPORT_AMOEBA void OPENMM_HIPPONONBONDEDFORCE_GETDPMEPARAMETERSINCONTEXT(const OpenMM_HippoNonbondedForce*& target, const OpenMM_Context*& context, double* alpha, int* nx, int* ny, int* nz) {
    OpenMM_HippoNonbondedForce_getDPMEParametersInContext(target, context, alpha, nx, ny, nz);
}
OPENMM_EXPORT_AMOEBA int openmm_hippononbondedforce_addparticle_(OpenMM_HippoNonbondedForce*& target, double const& charge, const OpenMM_DoubleArray*& dipole, const OpenMM_DoubleArray*& quadrupole, double const& coreCharge, double const& alpha, double const& epsilon, double const& damping, double const& c6, double const& pauliK, double const& pauliQ, double const& pauliAlpha, double const& polarizability, int const& axisType, int const& multipoleAtomZ, int const& multipoleAtomX, int const& multipoleAtomY) {
    return OpenMM_HippoNonbondedForce_addParticle(target, charge, dipole, quadrupole, coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha, polarizability, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY);
}
OPENMM_EXPORT_AMOEBA int OPENMM_HIPPONONBONDEDFORCE_ADDPARTICLE(OpenMM_HippoNonbondedForce*& target, double const& charge, const OpenMM_DoubleArray*& dipole, const OpenMM_DoubleArray*& quadrupole, double const& coreCharge, double const& alpha, double const& epsilon, double const& damping, double const& c6, double const& pauliK, double const& pauliQ, double const& pauliAlpha, double const& polarizability, int const& axisType, int const& multipoleAtomZ, int const& multipoleAtomX, int const& multipoleAtomY) {
    return OpenMM_HippoNonbondedForce_addParticle(target, charge, dipole, quadrupole, coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha, polarizability, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY);
}
OPENMM_EXPORT_AMOEBA void openmm_hippononbondedforce_getparticleparameters_(const OpenMM_HippoNonbondedForce*& target, int const& index, double* charge, OpenMM_DoubleArray*& dipole, OpenMM_DoubleArray*& quadrupole, double* coreCharge, double* alpha, double* epsilon, double* damping, double* c6, double* pauliK, double* pauliQ, double* pauliAlpha, double* polarizability, int* axisType, int* multipoleAtomZ, int* multipoleAtomX, int* multipoleAtomY) {
    OpenMM_HippoNonbondedForce_getParticleParameters(target, index, charge, dipole, quadrupole, coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha, polarizability, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY);
}
OPENMM_EXPORT_AMOEBA void OPENMM_HIPPONONBONDEDFORCE_GETPARTICLEPARAMETERS(const OpenMM_HippoNonbondedForce*& target, int const& index, double* charge, OpenMM_DoubleArray*& dipole, OpenMM_DoubleArray*& quadrupole, double* coreCharge, double* alpha, double* epsilon, double* damping, double* c6, double* pauliK, double* pauliQ, double* pauliAlpha, double* polarizability, int* axisType, int* multipoleAtomZ, int* multipoleAtomX, int* multipoleAtomY) {
    OpenMM_HippoNonbondedForce_getParticleParameters(target, index, charge, dipole, quadrupole, coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha, polarizability, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY);
}
OPENMM_EXPORT_AMOEBA void openmm_hippononbondedforce_setparticleparameters_(OpenMM_HippoNonbondedForce*& target, int const& index, double const& charge, const OpenMM_DoubleArray*& dipole, const OpenMM_DoubleArray*& quadrupole, double const& coreCharge, double const& alpha, double const& epsilon, double const& damping, double const& c6, double const& pauliK, double const& pauliQ, double const& pauliAlpha, double const& polarizability, int const& axisType, int const& multipoleAtomZ, int const& multipoleAtomX, int const& multipoleAtomY) {
    OpenMM_HippoNonbondedForce_setParticleParameters(target, index, charge, dipole, quadrupole, coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha, polarizability, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY);
}
OPENMM_EXPORT_AMOEBA void OPENMM_HIPPONONBONDEDFORCE_SETPARTICLEPARAMETERS(OpenMM_HippoNonbondedForce*& target, int const& index, double const& charge, const OpenMM_DoubleArray*& dipole, const OpenMM_DoubleArray*& quadrupole, double const& coreCharge, double const& alpha, double const& epsilon, double const& damping, double const& c6, double const& pauliK, double const& pauliQ, double const& pauliAlpha, double const& polarizability, int const& axisType, int const& multipoleAtomZ, int const& multipoleAtomX, int const& multipoleAtomY) {
    OpenMM_HippoNonbondedForce_setParticleParameters(target, index, charge, dipole, quadrupole, coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha, polarizability, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY);
}
OPENMM_EXPORT_AMOEBA int openmm_hippononbondedforce_addexception_(OpenMM_HippoNonbondedForce*& target, int const& particle1, int const& particle2, double const& multipoleMultipoleScale, double const& dipoleMultipoleScale, double const& dipoleDipoleScale, double const& dispersionScale, double const& repulsionScale, double const& chargeTransferScale, OpenMM_Boolean& replace) {
    return OpenMM_HippoNonbondedForce_addException(target, particle1, particle2, multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale, replace);
}
OPENMM_EXPORT_AMOEBA int OPENMM_HIPPONONBONDEDFORCE_ADDEXCEPTION(OpenMM_HippoNonbondedForce*& target, int const& particle1, int const& particle2, double const& multipoleMultipoleScale, double const& dipoleMultipoleScale, double const& dipoleDipoleScale, double const& dispersionScale, double const& repulsionScale, double const& chargeTransferScale, OpenMM_Boolean& replace) {
    return OpenMM_HippoNonbondedForce_addException(target, particle1, particle2, multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale, replace);
}
OPENMM_EXPORT_AMOEBA void openmm_hippononbondedforce_getexceptionparameters_(const OpenMM_HippoNonbondedForce*& target, int const& index, int* particle1, int* particle2, double* multipoleMultipoleScale, double* dipoleMultipoleScale, double* dipoleDipoleScale, double* dispersionScale, double* repulsionScale, double* chargeTransferScale) {
    OpenMM_HippoNonbondedForce_getExceptionParameters(target, index, particle1, particle2, multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale);
}
OPENMM_EXPORT_AMOEBA void OPENMM_HIPPONONBONDEDFORCE_GETEXCEPTIONPARAMETERS(const OpenMM_HippoNonbondedForce*& target, int const& index, int* particle1, int* particle2, double* multipoleMultipoleScale, double* dipoleMultipoleScale, double* dipoleDipoleScale, double* dispersionScale, double* repulsionScale, double* chargeTransferScale) {
    OpenMM_HippoNonbondedForce_getExceptionParameters(target, index, particle1, particle2, multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale);
}
OPENMM_EXPORT_AMOEBA void openmm_hippononbondedforce_setexceptionparameters_(OpenMM_HippoNonbondedForce*& target, int const& index, int const& particle1, int const& particle2, double const& multipoleMultipoleScale, double const& dipoleMultipoleScale, double const& dipoleDipoleScale, double const& dispersionScale, double const& repulsionScale, double const& chargeTransferScale) {
    OpenMM_HippoNonbondedForce_setExceptionParameters(target, index, particle1, particle2, multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale);
}
OPENMM_EXPORT_AMOEBA void OPENMM_HIPPONONBONDEDFORCE_SETEXCEPTIONPARAMETERS(OpenMM_HippoNonbondedForce*& target, int const& index, int const& particle1, int const& particle2, double const& multipoleMultipoleScale, double const& dipoleMultipoleScale, double const& dipoleDipoleScale, double const& dispersionScale, double const& repulsionScale, double const& chargeTransferScale) {
    OpenMM_HippoNonbondedForce_setExceptionParameters(target, index, particle1, particle2, multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale);
}
OPENMM_EXPORT_AMOEBA double openmm_hippononbondedforce_getewalderrortolerance_(const OpenMM_HippoNonbondedForce*& target) {
    return OpenMM_HippoNonbondedForce_getEwaldErrorTolerance(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_HIPPONONBONDEDFORCE_GETEWALDERRORTOLERANCE(const OpenMM_HippoNonbondedForce*& target) {
    return OpenMM_HippoNonbondedForce_getEwaldErrorTolerance(target);
}
OPENMM_EXPORT_AMOEBA void openmm_hippononbondedforce_setewalderrortolerance_(OpenMM_HippoNonbondedForce*& target, double const& tol) {
    OpenMM_HippoNonbondedForce_setEwaldErrorTolerance(target, tol);
}
OPENMM_EXPORT_AMOEBA void OPENMM_HIPPONONBONDEDFORCE_SETEWALDERRORTOLERANCE(OpenMM_HippoNonbondedForce*& target, double const& tol) {
    OpenMM_HippoNonbondedForce_setEwaldErrorTolerance(target, tol);
}
OPENMM_EXPORT_AMOEBA void openmm_hippononbondedforce_getlabframepermanentdipoles_(OpenMM_HippoNonbondedForce*& target, OpenMM_Context*& context, OpenMM_Vec3Array*& dipoles) {
    OpenMM_HippoNonbondedForce_getLabFramePermanentDipoles(target, context, dipoles);
}
OPENMM_EXPORT_AMOEBA void OPENMM_HIPPONONBONDEDFORCE_GETLABFRAMEPERMANENTDIPOLES(OpenMM_HippoNonbondedForce*& target, OpenMM_Context*& context, OpenMM_Vec3Array*& dipoles) {
    OpenMM_HippoNonbondedForce_getLabFramePermanentDipoles(target, context, dipoles);
}
OPENMM_EXPORT_AMOEBA void openmm_hippononbondedforce_getinduceddipoles_(OpenMM_HippoNonbondedForce*& target, OpenMM_Context*& context, OpenMM_Vec3Array*& dipoles) {
    OpenMM_HippoNonbondedForce_getInducedDipoles(target, context, dipoles);
}
OPENMM_EXPORT_AMOEBA void OPENMM_HIPPONONBONDEDFORCE_GETINDUCEDDIPOLES(OpenMM_HippoNonbondedForce*& target, OpenMM_Context*& context, OpenMM_Vec3Array*& dipoles) {
    OpenMM_HippoNonbondedForce_getInducedDipoles(target, context, dipoles);
}
OPENMM_EXPORT_AMOEBA void openmm_hippononbondedforce_updateparametersincontext_(OpenMM_HippoNonbondedForce*& target, OpenMM_Context*& context) {
    OpenMM_HippoNonbondedForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT_AMOEBA void OPENMM_HIPPONONBONDEDFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_HippoNonbondedForce*& target, OpenMM_Context*& context) {
    OpenMM_HippoNonbondedForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT_AMOEBA void openmm_hippononbondedforce_usesperiodicboundaryconditions_(const OpenMM_HippoNonbondedForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_HippoNonbondedForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT_AMOEBA void OPENMM_HIPPONONBONDEDFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_HippoNonbondedForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_HippoNonbondedForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::AmoebaAngleForce */
OPENMM_EXPORT_AMOEBA void openmm_amoebaangleforce_create_(OpenMM_AmoebaAngleForce*& result) {
    result = OpenMM_AmoebaAngleForce_create();
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAANGLEFORCE_CREATE(OpenMM_AmoebaAngleForce*& result) {
    result = OpenMM_AmoebaAngleForce_create();
}
OPENMM_EXPORT_AMOEBA void openmm_amoebaangleforce_destroy_(OpenMM_AmoebaAngleForce*& destroy) {
    OpenMM_AmoebaAngleForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAANGLEFORCE_DESTROY(OpenMM_AmoebaAngleForce*& destroy) {
    OpenMM_AmoebaAngleForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT_AMOEBA int openmm_amoebaangleforce_getnumangles_(const OpenMM_AmoebaAngleForce*& target) {
    return OpenMM_AmoebaAngleForce_getNumAngles(target);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBAANGLEFORCE_GETNUMANGLES(const OpenMM_AmoebaAngleForce*& target) {
    return OpenMM_AmoebaAngleForce_getNumAngles(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebaangleforce_setamoebaglobalanglecubic_(OpenMM_AmoebaAngleForce*& target, double const& cubicK) {
    OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleCubic(target, cubicK);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAANGLEFORCE_SETAMOEBAGLOBALANGLECUBIC(OpenMM_AmoebaAngleForce*& target, double const& cubicK) {
    OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleCubic(target, cubicK);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebaangleforce_getamoebaglobalanglecubic_(const OpenMM_AmoebaAngleForce*& target) {
    return OpenMM_AmoebaAngleForce_getAmoebaGlobalAngleCubic(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAANGLEFORCE_GETAMOEBAGLOBALANGLECUBIC(const OpenMM_AmoebaAngleForce*& target) {
    return OpenMM_AmoebaAngleForce_getAmoebaGlobalAngleCubic(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebaangleforce_setamoebaglobalanglequartic_(OpenMM_AmoebaAngleForce*& target, double const& quarticK) {
    OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleQuartic(target, quarticK);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAANGLEFORCE_SETAMOEBAGLOBALANGLEQUARTIC(OpenMM_AmoebaAngleForce*& target, double const& quarticK) {
    OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleQuartic(target, quarticK);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebaangleforce_getamoebaglobalanglequartic_(const OpenMM_AmoebaAngleForce*& target) {
    return OpenMM_AmoebaAngleForce_getAmoebaGlobalAngleQuartic(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAANGLEFORCE_GETAMOEBAGLOBALANGLEQUARTIC(const OpenMM_AmoebaAngleForce*& target) {
    return OpenMM_AmoebaAngleForce_getAmoebaGlobalAngleQuartic(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebaangleforce_setamoebaglobalanglepentic_(OpenMM_AmoebaAngleForce*& target, double const& penticK) {
    OpenMM_AmoebaAngleForce_setAmoebaGlobalAnglePentic(target, penticK);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAANGLEFORCE_SETAMOEBAGLOBALANGLEPENTIC(OpenMM_AmoebaAngleForce*& target, double const& penticK) {
    OpenMM_AmoebaAngleForce_setAmoebaGlobalAnglePentic(target, penticK);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebaangleforce_getamoebaglobalanglepentic_(const OpenMM_AmoebaAngleForce*& target) {
    return OpenMM_AmoebaAngleForce_getAmoebaGlobalAnglePentic(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAANGLEFORCE_GETAMOEBAGLOBALANGLEPENTIC(const OpenMM_AmoebaAngleForce*& target) {
    return OpenMM_AmoebaAngleForce_getAmoebaGlobalAnglePentic(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebaangleforce_setamoebaglobalanglesextic_(OpenMM_AmoebaAngleForce*& target, double const& sexticK) {
    OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleSextic(target, sexticK);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAANGLEFORCE_SETAMOEBAGLOBALANGLESEXTIC(OpenMM_AmoebaAngleForce*& target, double const& sexticK) {
    OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleSextic(target, sexticK);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebaangleforce_getamoebaglobalanglesextic_(const OpenMM_AmoebaAngleForce*& target) {
    return OpenMM_AmoebaAngleForce_getAmoebaGlobalAngleSextic(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAANGLEFORCE_GETAMOEBAGLOBALANGLESEXTIC(const OpenMM_AmoebaAngleForce*& target) {
    return OpenMM_AmoebaAngleForce_getAmoebaGlobalAngleSextic(target);
}
OPENMM_EXPORT_AMOEBA int openmm_amoebaangleforce_addangle_(OpenMM_AmoebaAngleForce*& target, int const& particle1, int const& particle2, int const& particle3, double const& length, double const& quadraticK) {
    return OpenMM_AmoebaAngleForce_addAngle(target, particle1, particle2, particle3, length, quadraticK);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBAANGLEFORCE_ADDANGLE(OpenMM_AmoebaAngleForce*& target, int const& particle1, int const& particle2, int const& particle3, double const& length, double const& quadraticK) {
    return OpenMM_AmoebaAngleForce_addAngle(target, particle1, particle2, particle3, length, quadraticK);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebaangleforce_getangleparameters_(const OpenMM_AmoebaAngleForce*& target, int const& index, int* particle1, int* particle2, int* particle3, double* length, double* quadraticK) {
    OpenMM_AmoebaAngleForce_getAngleParameters(target, index, particle1, particle2, particle3, length, quadraticK);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAANGLEFORCE_GETANGLEPARAMETERS(const OpenMM_AmoebaAngleForce*& target, int const& index, int* particle1, int* particle2, int* particle3, double* length, double* quadraticK) {
    OpenMM_AmoebaAngleForce_getAngleParameters(target, index, particle1, particle2, particle3, length, quadraticK);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebaangleforce_setangleparameters_(OpenMM_AmoebaAngleForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, double const& length, double const& quadraticK) {
    OpenMM_AmoebaAngleForce_setAngleParameters(target, index, particle1, particle2, particle3, length, quadraticK);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAANGLEFORCE_SETANGLEPARAMETERS(OpenMM_AmoebaAngleForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, double const& length, double const& quadraticK) {
    OpenMM_AmoebaAngleForce_setAngleParameters(target, index, particle1, particle2, particle3, length, quadraticK);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebaangleforce_updateparametersincontext_(OpenMM_AmoebaAngleForce*& target, OpenMM_Context*& context) {
    OpenMM_AmoebaAngleForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAANGLEFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_AmoebaAngleForce*& target, OpenMM_Context*& context) {
    OpenMM_AmoebaAngleForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebaangleforce_setusesperiodicboundaryconditions_(OpenMM_AmoebaAngleForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_AmoebaAngleForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAANGLEFORCE_SETUSESPERIODICBOUNDARYCONDITIONS(OpenMM_AmoebaAngleForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_AmoebaAngleForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebaangleforce_usesperiodicboundaryconditions_(const OpenMM_AmoebaAngleForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_AmoebaAngleForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAANGLEFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_AmoebaAngleForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_AmoebaAngleForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::AmoebaInPlaneAngleForce */
OPENMM_EXPORT_AMOEBA void openmm_amoebainplaneangleforce_create_(OpenMM_AmoebaInPlaneAngleForce*& result) {
    result = OpenMM_AmoebaInPlaneAngleForce_create();
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAINPLANEANGLEFORCE_CREATE(OpenMM_AmoebaInPlaneAngleForce*& result) {
    result = OpenMM_AmoebaInPlaneAngleForce_create();
}
OPENMM_EXPORT_AMOEBA void openmm_amoebainplaneangleforce_destroy_(OpenMM_AmoebaInPlaneAngleForce*& destroy) {
    OpenMM_AmoebaInPlaneAngleForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAINPLANEANGLEFORCE_DESTROY(OpenMM_AmoebaInPlaneAngleForce*& destroy) {
    OpenMM_AmoebaInPlaneAngleForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT_AMOEBA int openmm_amoebainplaneangleforce_getnumangles_(const OpenMM_AmoebaInPlaneAngleForce*& target) {
    return OpenMM_AmoebaInPlaneAngleForce_getNumAngles(target);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBAINPLANEANGLEFORCE_GETNUMANGLES(const OpenMM_AmoebaInPlaneAngleForce*& target) {
    return OpenMM_AmoebaInPlaneAngleForce_getNumAngles(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebainplaneangleforce_setamoebaglobalinplaneanglecubic_(OpenMM_AmoebaInPlaneAngleForce*& target, double const& cubicK) {
    OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleCubic(target, cubicK);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAINPLANEANGLEFORCE_SETAMOEBAGLOBALINPLANEANGLECUBIC(OpenMM_AmoebaInPlaneAngleForce*& target, double const& cubicK) {
    OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleCubic(target, cubicK);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebainplaneangleforce_getamoebaglobalinplaneanglecubic_(const OpenMM_AmoebaInPlaneAngleForce*& target) {
    return OpenMM_AmoebaInPlaneAngleForce_getAmoebaGlobalInPlaneAngleCubic(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAINPLANEANGLEFORCE_GETAMOEBAGLOBALINPLANEANGLECUBIC(const OpenMM_AmoebaInPlaneAngleForce*& target) {
    return OpenMM_AmoebaInPlaneAngleForce_getAmoebaGlobalInPlaneAngleCubic(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebainplaneangleforce_setamoebaglobalinplaneanglequartic_(OpenMM_AmoebaInPlaneAngleForce*& target, double const& quarticK) {
    OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleQuartic(target, quarticK);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAINPLANEANGLEFORCE_SETAMOEBAGLOBALINPLANEANGLEQUARTIC(OpenMM_AmoebaInPlaneAngleForce*& target, double const& quarticK) {
    OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleQuartic(target, quarticK);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebainplaneangleforce_getamoebaglobalinplaneanglequartic_(const OpenMM_AmoebaInPlaneAngleForce*& target) {
    return OpenMM_AmoebaInPlaneAngleForce_getAmoebaGlobalInPlaneAngleQuartic(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAINPLANEANGLEFORCE_GETAMOEBAGLOBALINPLANEANGLEQUARTIC(const OpenMM_AmoebaInPlaneAngleForce*& target) {
    return OpenMM_AmoebaInPlaneAngleForce_getAmoebaGlobalInPlaneAngleQuartic(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebainplaneangleforce_setamoebaglobalinplaneanglepentic_(OpenMM_AmoebaInPlaneAngleForce*& target, double const& penticK) {
    OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAnglePentic(target, penticK);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAINPLANEANGLEFORCE_SETAMOEBAGLOBALINPLANEANGLEPENTIC(OpenMM_AmoebaInPlaneAngleForce*& target, double const& penticK) {
    OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAnglePentic(target, penticK);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebainplaneangleforce_getamoebaglobalinplaneanglepentic_(const OpenMM_AmoebaInPlaneAngleForce*& target) {
    return OpenMM_AmoebaInPlaneAngleForce_getAmoebaGlobalInPlaneAnglePentic(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAINPLANEANGLEFORCE_GETAMOEBAGLOBALINPLANEANGLEPENTIC(const OpenMM_AmoebaInPlaneAngleForce*& target) {
    return OpenMM_AmoebaInPlaneAngleForce_getAmoebaGlobalInPlaneAnglePentic(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebainplaneangleforce_setamoebaglobalinplaneanglesextic_(OpenMM_AmoebaInPlaneAngleForce*& target, double const& sexticK) {
    OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleSextic(target, sexticK);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAINPLANEANGLEFORCE_SETAMOEBAGLOBALINPLANEANGLESEXTIC(OpenMM_AmoebaInPlaneAngleForce*& target, double const& sexticK) {
    OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleSextic(target, sexticK);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebainplaneangleforce_getamoebaglobalinplaneanglesextic_(const OpenMM_AmoebaInPlaneAngleForce*& target) {
    return OpenMM_AmoebaInPlaneAngleForce_getAmoebaGlobalInPlaneAngleSextic(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAINPLANEANGLEFORCE_GETAMOEBAGLOBALINPLANEANGLESEXTIC(const OpenMM_AmoebaInPlaneAngleForce*& target) {
    return OpenMM_AmoebaInPlaneAngleForce_getAmoebaGlobalInPlaneAngleSextic(target);
}
OPENMM_EXPORT_AMOEBA int openmm_amoebainplaneangleforce_addangle_(OpenMM_AmoebaInPlaneAngleForce*& target, int const& particle1, int const& particle2, int const& particle3, int const& particle4, double const& length, double const& quadraticK) {
    return OpenMM_AmoebaInPlaneAngleForce_addAngle(target, particle1, particle2, particle3, particle4, length, quadraticK);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBAINPLANEANGLEFORCE_ADDANGLE(OpenMM_AmoebaInPlaneAngleForce*& target, int const& particle1, int const& particle2, int const& particle3, int const& particle4, double const& length, double const& quadraticK) {
    return OpenMM_AmoebaInPlaneAngleForce_addAngle(target, particle1, particle2, particle3, particle4, length, quadraticK);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebainplaneangleforce_getangleparameters_(const OpenMM_AmoebaInPlaneAngleForce*& target, int const& index, int* particle1, int* particle2, int* particle3, int* particle4, double* length, double* quadraticK) {
    OpenMM_AmoebaInPlaneAngleForce_getAngleParameters(target, index, particle1, particle2, particle3, particle4, length, quadraticK);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAINPLANEANGLEFORCE_GETANGLEPARAMETERS(const OpenMM_AmoebaInPlaneAngleForce*& target, int const& index, int* particle1, int* particle2, int* particle3, int* particle4, double* length, double* quadraticK) {
    OpenMM_AmoebaInPlaneAngleForce_getAngleParameters(target, index, particle1, particle2, particle3, particle4, length, quadraticK);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebainplaneangleforce_setangleparameters_(OpenMM_AmoebaInPlaneAngleForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, int const& particle4, double const& length, double const& quadraticK) {
    OpenMM_AmoebaInPlaneAngleForce_setAngleParameters(target, index, particle1, particle2, particle3, particle4, length, quadraticK);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAINPLANEANGLEFORCE_SETANGLEPARAMETERS(OpenMM_AmoebaInPlaneAngleForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, int const& particle4, double const& length, double const& quadraticK) {
    OpenMM_AmoebaInPlaneAngleForce_setAngleParameters(target, index, particle1, particle2, particle3, particle4, length, quadraticK);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebainplaneangleforce_updateparametersincontext_(OpenMM_AmoebaInPlaneAngleForce*& target, OpenMM_Context*& context) {
    OpenMM_AmoebaInPlaneAngleForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAINPLANEANGLEFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_AmoebaInPlaneAngleForce*& target, OpenMM_Context*& context) {
    OpenMM_AmoebaInPlaneAngleForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebainplaneangleforce_setusesperiodicboundaryconditions_(OpenMM_AmoebaInPlaneAngleForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_AmoebaInPlaneAngleForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAINPLANEANGLEFORCE_SETUSESPERIODICBOUNDARYCONDITIONS(OpenMM_AmoebaInPlaneAngleForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_AmoebaInPlaneAngleForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebainplaneangleforce_usesperiodicboundaryconditions_(const OpenMM_AmoebaInPlaneAngleForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_AmoebaInPlaneAngleForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAINPLANEANGLEFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_AmoebaInPlaneAngleForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_AmoebaInPlaneAngleForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::AmoebaBondForce */
OPENMM_EXPORT_AMOEBA void openmm_amoebabondforce_create_(OpenMM_AmoebaBondForce*& result) {
    result = OpenMM_AmoebaBondForce_create();
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBABONDFORCE_CREATE(OpenMM_AmoebaBondForce*& result) {
    result = OpenMM_AmoebaBondForce_create();
}
OPENMM_EXPORT_AMOEBA void openmm_amoebabondforce_destroy_(OpenMM_AmoebaBondForce*& destroy) {
    OpenMM_AmoebaBondForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBABONDFORCE_DESTROY(OpenMM_AmoebaBondForce*& destroy) {
    OpenMM_AmoebaBondForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT_AMOEBA int openmm_amoebabondforce_getnumbonds_(const OpenMM_AmoebaBondForce*& target) {
    return OpenMM_AmoebaBondForce_getNumBonds(target);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBABONDFORCE_GETNUMBONDS(const OpenMM_AmoebaBondForce*& target) {
    return OpenMM_AmoebaBondForce_getNumBonds(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebabondforce_setamoebaglobalbondcubic_(OpenMM_AmoebaBondForce*& target, double const& cubicK) {
    OpenMM_AmoebaBondForce_setAmoebaGlobalBondCubic(target, cubicK);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBABONDFORCE_SETAMOEBAGLOBALBONDCUBIC(OpenMM_AmoebaBondForce*& target, double const& cubicK) {
    OpenMM_AmoebaBondForce_setAmoebaGlobalBondCubic(target, cubicK);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebabondforce_getamoebaglobalbondcubic_(const OpenMM_AmoebaBondForce*& target) {
    return OpenMM_AmoebaBondForce_getAmoebaGlobalBondCubic(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBABONDFORCE_GETAMOEBAGLOBALBONDCUBIC(const OpenMM_AmoebaBondForce*& target) {
    return OpenMM_AmoebaBondForce_getAmoebaGlobalBondCubic(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebabondforce_setamoebaglobalbondquartic_(OpenMM_AmoebaBondForce*& target, double const& quarticK) {
    OpenMM_AmoebaBondForce_setAmoebaGlobalBondQuartic(target, quarticK);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBABONDFORCE_SETAMOEBAGLOBALBONDQUARTIC(OpenMM_AmoebaBondForce*& target, double const& quarticK) {
    OpenMM_AmoebaBondForce_setAmoebaGlobalBondQuartic(target, quarticK);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebabondforce_getamoebaglobalbondquartic_(const OpenMM_AmoebaBondForce*& target) {
    return OpenMM_AmoebaBondForce_getAmoebaGlobalBondQuartic(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBABONDFORCE_GETAMOEBAGLOBALBONDQUARTIC(const OpenMM_AmoebaBondForce*& target) {
    return OpenMM_AmoebaBondForce_getAmoebaGlobalBondQuartic(target);
}
OPENMM_EXPORT_AMOEBA int openmm_amoebabondforce_addbond_(OpenMM_AmoebaBondForce*& target, int const& particle1, int const& particle2, double const& length, double const& quadraticK) {
    return OpenMM_AmoebaBondForce_addBond(target, particle1, particle2, length, quadraticK);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBABONDFORCE_ADDBOND(OpenMM_AmoebaBondForce*& target, int const& particle1, int const& particle2, double const& length, double const& quadraticK) {
    return OpenMM_AmoebaBondForce_addBond(target, particle1, particle2, length, quadraticK);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebabondforce_getbondparameters_(const OpenMM_AmoebaBondForce*& target, int const& index, int* particle1, int* particle2, double* length, double* quadraticK) {
    OpenMM_AmoebaBondForce_getBondParameters(target, index, particle1, particle2, length, quadraticK);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBABONDFORCE_GETBONDPARAMETERS(const OpenMM_AmoebaBondForce*& target, int const& index, int* particle1, int* particle2, double* length, double* quadraticK) {
    OpenMM_AmoebaBondForce_getBondParameters(target, index, particle1, particle2, length, quadraticK);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebabondforce_setbondparameters_(OpenMM_AmoebaBondForce*& target, int const& index, int const& particle1, int const& particle2, double const& length, double const& quadraticK) {
    OpenMM_AmoebaBondForce_setBondParameters(target, index, particle1, particle2, length, quadraticK);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBABONDFORCE_SETBONDPARAMETERS(OpenMM_AmoebaBondForce*& target, int const& index, int const& particle1, int const& particle2, double const& length, double const& quadraticK) {
    OpenMM_AmoebaBondForce_setBondParameters(target, index, particle1, particle2, length, quadraticK);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebabondforce_updateparametersincontext_(OpenMM_AmoebaBondForce*& target, OpenMM_Context*& context) {
    OpenMM_AmoebaBondForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBABONDFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_AmoebaBondForce*& target, OpenMM_Context*& context) {
    OpenMM_AmoebaBondForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebabondforce_setusesperiodicboundaryconditions_(OpenMM_AmoebaBondForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_AmoebaBondForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBABONDFORCE_SETUSESPERIODICBOUNDARYCONDITIONS(OpenMM_AmoebaBondForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_AmoebaBondForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebabondforce_usesperiodicboundaryconditions_(const OpenMM_AmoebaBondForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_AmoebaBondForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBABONDFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_AmoebaBondForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_AmoebaBondForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::AmoebaWcaDispersionForce */
OPENMM_EXPORT_AMOEBA void openmm_amoebawcadispersionforce_create_(OpenMM_AmoebaWcaDispersionForce*& result) {
    result = OpenMM_AmoebaWcaDispersionForce_create();
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAWCADISPERSIONFORCE_CREATE(OpenMM_AmoebaWcaDispersionForce*& result) {
    result = OpenMM_AmoebaWcaDispersionForce_create();
}
OPENMM_EXPORT_AMOEBA void openmm_amoebawcadispersionforce_destroy_(OpenMM_AmoebaWcaDispersionForce*& destroy) {
    OpenMM_AmoebaWcaDispersionForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAWCADISPERSIONFORCE_DESTROY(OpenMM_AmoebaWcaDispersionForce*& destroy) {
    OpenMM_AmoebaWcaDispersionForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT_AMOEBA int openmm_amoebawcadispersionforce_getnumparticles_(const OpenMM_AmoebaWcaDispersionForce*& target) {
    return OpenMM_AmoebaWcaDispersionForce_getNumParticles(target);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBAWCADISPERSIONFORCE_GETNUMPARTICLES(const OpenMM_AmoebaWcaDispersionForce*& target) {
    return OpenMM_AmoebaWcaDispersionForce_getNumParticles(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebawcadispersionforce_setparticleparameters_(OpenMM_AmoebaWcaDispersionForce*& target, int const& particleIndex, double const& radius, double const& epsilon) {
    OpenMM_AmoebaWcaDispersionForce_setParticleParameters(target, particleIndex, radius, epsilon);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAWCADISPERSIONFORCE_SETPARTICLEPARAMETERS(OpenMM_AmoebaWcaDispersionForce*& target, int const& particleIndex, double const& radius, double const& epsilon) {
    OpenMM_AmoebaWcaDispersionForce_setParticleParameters(target, particleIndex, radius, epsilon);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebawcadispersionforce_getparticleparameters_(const OpenMM_AmoebaWcaDispersionForce*& target, int const& particleIndex, double* radius, double* epsilon) {
    OpenMM_AmoebaWcaDispersionForce_getParticleParameters(target, particleIndex, radius, epsilon);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAWCADISPERSIONFORCE_GETPARTICLEPARAMETERS(const OpenMM_AmoebaWcaDispersionForce*& target, int const& particleIndex, double* radius, double* epsilon) {
    OpenMM_AmoebaWcaDispersionForce_getParticleParameters(target, particleIndex, radius, epsilon);
}
OPENMM_EXPORT_AMOEBA int openmm_amoebawcadispersionforce_addparticle_(OpenMM_AmoebaWcaDispersionForce*& target, double const& radius, double const& epsilon) {
    return OpenMM_AmoebaWcaDispersionForce_addParticle(target, radius, epsilon);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBAWCADISPERSIONFORCE_ADDPARTICLE(OpenMM_AmoebaWcaDispersionForce*& target, double const& radius, double const& epsilon) {
    return OpenMM_AmoebaWcaDispersionForce_addParticle(target, radius, epsilon);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebawcadispersionforce_updateparametersincontext_(OpenMM_AmoebaWcaDispersionForce*& target, OpenMM_Context*& context) {
    OpenMM_AmoebaWcaDispersionForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAWCADISPERSIONFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_AmoebaWcaDispersionForce*& target, OpenMM_Context*& context) {
    OpenMM_AmoebaWcaDispersionForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebawcadispersionforce_getepso_(const OpenMM_AmoebaWcaDispersionForce*& target) {
    return OpenMM_AmoebaWcaDispersionForce_getEpso(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAWCADISPERSIONFORCE_GETEPSO(const OpenMM_AmoebaWcaDispersionForce*& target) {
    return OpenMM_AmoebaWcaDispersionForce_getEpso(target);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebawcadispersionforce_getepsh_(const OpenMM_AmoebaWcaDispersionForce*& target) {
    return OpenMM_AmoebaWcaDispersionForce_getEpsh(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAWCADISPERSIONFORCE_GETEPSH(const OpenMM_AmoebaWcaDispersionForce*& target) {
    return OpenMM_AmoebaWcaDispersionForce_getEpsh(target);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebawcadispersionforce_getrmino_(const OpenMM_AmoebaWcaDispersionForce*& target) {
    return OpenMM_AmoebaWcaDispersionForce_getRmino(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAWCADISPERSIONFORCE_GETRMINO(const OpenMM_AmoebaWcaDispersionForce*& target) {
    return OpenMM_AmoebaWcaDispersionForce_getRmino(target);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebawcadispersionforce_getrminh_(const OpenMM_AmoebaWcaDispersionForce*& target) {
    return OpenMM_AmoebaWcaDispersionForce_getRminh(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAWCADISPERSIONFORCE_GETRMINH(const OpenMM_AmoebaWcaDispersionForce*& target) {
    return OpenMM_AmoebaWcaDispersionForce_getRminh(target);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebawcadispersionforce_getawater_(const OpenMM_AmoebaWcaDispersionForce*& target) {
    return OpenMM_AmoebaWcaDispersionForce_getAwater(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAWCADISPERSIONFORCE_GETAWATER(const OpenMM_AmoebaWcaDispersionForce*& target) {
    return OpenMM_AmoebaWcaDispersionForce_getAwater(target);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebawcadispersionforce_getshctd_(const OpenMM_AmoebaWcaDispersionForce*& target) {
    return OpenMM_AmoebaWcaDispersionForce_getShctd(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAWCADISPERSIONFORCE_GETSHCTD(const OpenMM_AmoebaWcaDispersionForce*& target) {
    return OpenMM_AmoebaWcaDispersionForce_getShctd(target);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebawcadispersionforce_getdispoff_(const OpenMM_AmoebaWcaDispersionForce*& target) {
    return OpenMM_AmoebaWcaDispersionForce_getDispoff(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAWCADISPERSIONFORCE_GETDISPOFF(const OpenMM_AmoebaWcaDispersionForce*& target) {
    return OpenMM_AmoebaWcaDispersionForce_getDispoff(target);
}
OPENMM_EXPORT_AMOEBA double openmm_amoebawcadispersionforce_getslevy_(const OpenMM_AmoebaWcaDispersionForce*& target) {
    return OpenMM_AmoebaWcaDispersionForce_getSlevy(target);
}
OPENMM_EXPORT_AMOEBA double OPENMM_AMOEBAWCADISPERSIONFORCE_GETSLEVY(const OpenMM_AmoebaWcaDispersionForce*& target) {
    return OpenMM_AmoebaWcaDispersionForce_getSlevy(target);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebawcadispersionforce_setepso_(OpenMM_AmoebaWcaDispersionForce*& target, double const& inputValue) {
    OpenMM_AmoebaWcaDispersionForce_setEpso(target, inputValue);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAWCADISPERSIONFORCE_SETEPSO(OpenMM_AmoebaWcaDispersionForce*& target, double const& inputValue) {
    OpenMM_AmoebaWcaDispersionForce_setEpso(target, inputValue);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebawcadispersionforce_setepsh_(OpenMM_AmoebaWcaDispersionForce*& target, double const& inputValue) {
    OpenMM_AmoebaWcaDispersionForce_setEpsh(target, inputValue);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAWCADISPERSIONFORCE_SETEPSH(OpenMM_AmoebaWcaDispersionForce*& target, double const& inputValue) {
    OpenMM_AmoebaWcaDispersionForce_setEpsh(target, inputValue);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebawcadispersionforce_setrmino_(OpenMM_AmoebaWcaDispersionForce*& target, double const& inputValue) {
    OpenMM_AmoebaWcaDispersionForce_setRmino(target, inputValue);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAWCADISPERSIONFORCE_SETRMINO(OpenMM_AmoebaWcaDispersionForce*& target, double const& inputValue) {
    OpenMM_AmoebaWcaDispersionForce_setRmino(target, inputValue);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebawcadispersionforce_setrminh_(OpenMM_AmoebaWcaDispersionForce*& target, double const& inputValue) {
    OpenMM_AmoebaWcaDispersionForce_setRminh(target, inputValue);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAWCADISPERSIONFORCE_SETRMINH(OpenMM_AmoebaWcaDispersionForce*& target, double const& inputValue) {
    OpenMM_AmoebaWcaDispersionForce_setRminh(target, inputValue);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebawcadispersionforce_setawater_(OpenMM_AmoebaWcaDispersionForce*& target, double const& inputValue) {
    OpenMM_AmoebaWcaDispersionForce_setAwater(target, inputValue);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAWCADISPERSIONFORCE_SETAWATER(OpenMM_AmoebaWcaDispersionForce*& target, double const& inputValue) {
    OpenMM_AmoebaWcaDispersionForce_setAwater(target, inputValue);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebawcadispersionforce_setshctd_(OpenMM_AmoebaWcaDispersionForce*& target, double const& inputValue) {
    OpenMM_AmoebaWcaDispersionForce_setShctd(target, inputValue);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAWCADISPERSIONFORCE_SETSHCTD(OpenMM_AmoebaWcaDispersionForce*& target, double const& inputValue) {
    OpenMM_AmoebaWcaDispersionForce_setShctd(target, inputValue);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebawcadispersionforce_setdispoff_(OpenMM_AmoebaWcaDispersionForce*& target, double const& inputValue) {
    OpenMM_AmoebaWcaDispersionForce_setDispoff(target, inputValue);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAWCADISPERSIONFORCE_SETDISPOFF(OpenMM_AmoebaWcaDispersionForce*& target, double const& inputValue) {
    OpenMM_AmoebaWcaDispersionForce_setDispoff(target, inputValue);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebawcadispersionforce_setslevy_(OpenMM_AmoebaWcaDispersionForce*& target, double const& inputValue) {
    OpenMM_AmoebaWcaDispersionForce_setSlevy(target, inputValue);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAWCADISPERSIONFORCE_SETSLEVY(OpenMM_AmoebaWcaDispersionForce*& target, double const& inputValue) {
    OpenMM_AmoebaWcaDispersionForce_setSlevy(target, inputValue);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebawcadispersionforce_usesperiodicboundaryconditions_(const OpenMM_AmoebaWcaDispersionForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_AmoebaWcaDispersionForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAWCADISPERSIONFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_AmoebaWcaDispersionForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_AmoebaWcaDispersionForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::AmoebaPiTorsionForce */
OPENMM_EXPORT_AMOEBA void openmm_amoebapitorsionforce_create_(OpenMM_AmoebaPiTorsionForce*& result) {
    result = OpenMM_AmoebaPiTorsionForce_create();
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAPITORSIONFORCE_CREATE(OpenMM_AmoebaPiTorsionForce*& result) {
    result = OpenMM_AmoebaPiTorsionForce_create();
}
OPENMM_EXPORT_AMOEBA void openmm_amoebapitorsionforce_destroy_(OpenMM_AmoebaPiTorsionForce*& destroy) {
    OpenMM_AmoebaPiTorsionForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAPITORSIONFORCE_DESTROY(OpenMM_AmoebaPiTorsionForce*& destroy) {
    OpenMM_AmoebaPiTorsionForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT_AMOEBA int openmm_amoebapitorsionforce_getnumpitorsions_(const OpenMM_AmoebaPiTorsionForce*& target) {
    return OpenMM_AmoebaPiTorsionForce_getNumPiTorsions(target);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBAPITORSIONFORCE_GETNUMPITORSIONS(const OpenMM_AmoebaPiTorsionForce*& target) {
    return OpenMM_AmoebaPiTorsionForce_getNumPiTorsions(target);
}
OPENMM_EXPORT_AMOEBA int openmm_amoebapitorsionforce_addpitorsion_(OpenMM_AmoebaPiTorsionForce*& target, int const& particle1, int const& particle2, int const& particle3, int const& particle4, int const& particle5, int const& particle6, double const& k) {
    return OpenMM_AmoebaPiTorsionForce_addPiTorsion(target, particle1, particle2, particle3, particle4, particle5, particle6, k);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBAPITORSIONFORCE_ADDPITORSION(OpenMM_AmoebaPiTorsionForce*& target, int const& particle1, int const& particle2, int const& particle3, int const& particle4, int const& particle5, int const& particle6, double const& k) {
    return OpenMM_AmoebaPiTorsionForce_addPiTorsion(target, particle1, particle2, particle3, particle4, particle5, particle6, k);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebapitorsionforce_getpitorsionparameters_(const OpenMM_AmoebaPiTorsionForce*& target, int const& index, int* particle1, int* particle2, int* particle3, int* particle4, int* particle5, int* particle6, double* k) {
    OpenMM_AmoebaPiTorsionForce_getPiTorsionParameters(target, index, particle1, particle2, particle3, particle4, particle5, particle6, k);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAPITORSIONFORCE_GETPITORSIONPARAMETERS(const OpenMM_AmoebaPiTorsionForce*& target, int const& index, int* particle1, int* particle2, int* particle3, int* particle4, int* particle5, int* particle6, double* k) {
    OpenMM_AmoebaPiTorsionForce_getPiTorsionParameters(target, index, particle1, particle2, particle3, particle4, particle5, particle6, k);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebapitorsionforce_setpitorsionparameters_(OpenMM_AmoebaPiTorsionForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, int const& particle4, int const& particle5, int const& particle6, double const& k) {
    OpenMM_AmoebaPiTorsionForce_setPiTorsionParameters(target, index, particle1, particle2, particle3, particle4, particle5, particle6, k);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAPITORSIONFORCE_SETPITORSIONPARAMETERS(OpenMM_AmoebaPiTorsionForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, int const& particle4, int const& particle5, int const& particle6, double const& k) {
    OpenMM_AmoebaPiTorsionForce_setPiTorsionParameters(target, index, particle1, particle2, particle3, particle4, particle5, particle6, k);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebapitorsionforce_updateparametersincontext_(OpenMM_AmoebaPiTorsionForce*& target, OpenMM_Context*& context) {
    OpenMM_AmoebaPiTorsionForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAPITORSIONFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_AmoebaPiTorsionForce*& target, OpenMM_Context*& context) {
    OpenMM_AmoebaPiTorsionForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebapitorsionforce_setusesperiodicboundaryconditions_(OpenMM_AmoebaPiTorsionForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_AmoebaPiTorsionForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAPITORSIONFORCE_SETUSESPERIODICBOUNDARYCONDITIONS(OpenMM_AmoebaPiTorsionForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_AmoebaPiTorsionForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebapitorsionforce_usesperiodicboundaryconditions_(const OpenMM_AmoebaPiTorsionForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_AmoebaPiTorsionForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBAPITORSIONFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_AmoebaPiTorsionForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_AmoebaPiTorsionForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::AmoebaTorsionTorsionForce */
OPENMM_EXPORT_AMOEBA void openmm_amoebatorsiontorsionforce_create_(OpenMM_AmoebaTorsionTorsionForce*& result) {
    result = OpenMM_AmoebaTorsionTorsionForce_create();
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBATORSIONTORSIONFORCE_CREATE(OpenMM_AmoebaTorsionTorsionForce*& result) {
    result = OpenMM_AmoebaTorsionTorsionForce_create();
}
OPENMM_EXPORT_AMOEBA void openmm_amoebatorsiontorsionforce_destroy_(OpenMM_AmoebaTorsionTorsionForce*& destroy) {
    OpenMM_AmoebaTorsionTorsionForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBATORSIONTORSIONFORCE_DESTROY(OpenMM_AmoebaTorsionTorsionForce*& destroy) {
    OpenMM_AmoebaTorsionTorsionForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT_AMOEBA int openmm_amoebatorsiontorsionforce_getnumtorsiontorsions_(const OpenMM_AmoebaTorsionTorsionForce*& target) {
    return OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsions(target);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBATORSIONTORSIONFORCE_GETNUMTORSIONTORSIONS(const OpenMM_AmoebaTorsionTorsionForce*& target) {
    return OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsions(target);
}
OPENMM_EXPORT_AMOEBA int openmm_amoebatorsiontorsionforce_getnumtorsiontorsiongrids_(const OpenMM_AmoebaTorsionTorsionForce*& target) {
    return OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsionGrids(target);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBATORSIONTORSIONFORCE_GETNUMTORSIONTORSIONGRIDS(const OpenMM_AmoebaTorsionTorsionForce*& target) {
    return OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsionGrids(target);
}
OPENMM_EXPORT_AMOEBA int openmm_amoebatorsiontorsionforce_addtorsiontorsion_(OpenMM_AmoebaTorsionTorsionForce*& target, int const& particle1, int const& particle2, int const& particle3, int const& particle4, int const& particle5, int const& chiralCheckAtomIndex, int const& gridIndex) {
    return OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion(target, particle1, particle2, particle3, particle4, particle5, chiralCheckAtomIndex, gridIndex);
}
OPENMM_EXPORT_AMOEBA int OPENMM_AMOEBATORSIONTORSIONFORCE_ADDTORSIONTORSION(OpenMM_AmoebaTorsionTorsionForce*& target, int const& particle1, int const& particle2, int const& particle3, int const& particle4, int const& particle5, int const& chiralCheckAtomIndex, int const& gridIndex) {
    return OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion(target, particle1, particle2, particle3, particle4, particle5, chiralCheckAtomIndex, gridIndex);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebatorsiontorsionforce_gettorsiontorsionparameters_(const OpenMM_AmoebaTorsionTorsionForce*& target, int const& index, int* particle1, int* particle2, int* particle3, int* particle4, int* particle5, int* chiralCheckAtomIndex, int* gridIndex) {
    OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionParameters(target, index, particle1, particle2, particle3, particle4, particle5, chiralCheckAtomIndex, gridIndex);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBATORSIONTORSIONFORCE_GETTORSIONTORSIONPARAMETERS(const OpenMM_AmoebaTorsionTorsionForce*& target, int const& index, int* particle1, int* particle2, int* particle3, int* particle4, int* particle5, int* chiralCheckAtomIndex, int* gridIndex) {
    OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionParameters(target, index, particle1, particle2, particle3, particle4, particle5, chiralCheckAtomIndex, gridIndex);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebatorsiontorsionforce_settorsiontorsionparameters_(OpenMM_AmoebaTorsionTorsionForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, int const& particle4, int const& particle5, int const& chiralCheckAtomIndex, int const& gridIndex) {
    OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionParameters(target, index, particle1, particle2, particle3, particle4, particle5, chiralCheckAtomIndex, gridIndex);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBATORSIONTORSIONFORCE_SETTORSIONTORSIONPARAMETERS(OpenMM_AmoebaTorsionTorsionForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, int const& particle4, int const& particle5, int const& chiralCheckAtomIndex, int const& gridIndex) {
    OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionParameters(target, index, particle1, particle2, particle3, particle4, particle5, chiralCheckAtomIndex, gridIndex);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebatorsiontorsionforce_gettorsiontorsiongrid_(const OpenMM_AmoebaTorsionTorsionForce*& target, int const& index, const OpenMM_3D_DoubleArray*& result) {
    result = OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionGrid(target, index);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBATORSIONTORSIONFORCE_GETTORSIONTORSIONGRID(const OpenMM_AmoebaTorsionTorsionForce*& target, int const& index, const OpenMM_3D_DoubleArray*& result) {
    result = OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionGrid(target, index);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebatorsiontorsionforce_settorsiontorsiongrid_(OpenMM_AmoebaTorsionTorsionForce*& target, int const& index, const OpenMM_3D_DoubleArray*& grid) {
    OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid(target, index, grid);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBATORSIONTORSIONFORCE_SETTORSIONTORSIONGRID(OpenMM_AmoebaTorsionTorsionForce*& target, int const& index, const OpenMM_3D_DoubleArray*& grid) {
    OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid(target, index, grid);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebatorsiontorsionforce_setusesperiodicboundaryconditions_(OpenMM_AmoebaTorsionTorsionForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_AmoebaTorsionTorsionForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBATORSIONTORSIONFORCE_SETUSESPERIODICBOUNDARYCONDITIONS(OpenMM_AmoebaTorsionTorsionForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_AmoebaTorsionTorsionForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT_AMOEBA void openmm_amoebatorsiontorsionforce_usesperiodicboundaryconditions_(const OpenMM_AmoebaTorsionTorsionForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_AmoebaTorsionTorsionForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT_AMOEBA void OPENMM_AMOEBATORSIONTORSIONFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_AmoebaTorsionTorsionForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_AmoebaTorsionTorsionForce_usesPeriodicBoundaryConditions(target);
}

}
