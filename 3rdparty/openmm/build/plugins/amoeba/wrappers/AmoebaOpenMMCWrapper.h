
#ifndef AMOEBA_OPENMM_CWRAPPER_H_
#define AMOEBA_OPENMM_CWRAPPER_H_

#ifndef OPENMM_EXPORT_AMOEBA
#define OPENMM_EXPORT_AMOEBA
#endif
/* Global Constants */


/* Type Declarations */

typedef struct OpenMM_AmoebaMultipoleForce_struct OpenMM_AmoebaMultipoleForce;
typedef struct OpenMM_AmoebaStretchBendForce_struct OpenMM_AmoebaStretchBendForce;
typedef struct OpenMM_AmoebaVdwForce_struct OpenMM_AmoebaVdwForce;
typedef struct OpenMM_AmoebaGeneralizedKirkwoodForce_struct OpenMM_AmoebaGeneralizedKirkwoodForce;
typedef struct OpenMM_AmoebaOutOfPlaneBendForce_struct OpenMM_AmoebaOutOfPlaneBendForce;
typedef struct OpenMM_HippoNonbondedForce_struct OpenMM_HippoNonbondedForce;
typedef struct OpenMM_AmoebaAngleForce_struct OpenMM_AmoebaAngleForce;
typedef struct OpenMM_AmoebaInPlaneAngleForce_struct OpenMM_AmoebaInPlaneAngleForce;
typedef struct OpenMM_AmoebaBondForce_struct OpenMM_AmoebaBondForce;
typedef struct OpenMM_AmoebaWcaDispersionForce_struct OpenMM_AmoebaWcaDispersionForce;
typedef struct OpenMM_AmoebaPiTorsionForce_struct OpenMM_AmoebaPiTorsionForce;
typedef struct OpenMM_AmoebaTorsionTorsionForce_struct OpenMM_AmoebaTorsionTorsionForce;

typedef struct OpenMM_2D_IntArray_struct OpenMM_2D_IntArray;
typedef struct OpenMM_3D_DoubleArray_struct OpenMM_3D_DoubleArray;

#if defined(__cplusplus)
extern "C" {
#endif

/* OpenMM_3D_DoubleArray */
OPENMM_EXPORT_AMOEBA OpenMM_3D_DoubleArray* OpenMM_3D_DoubleArray_create(int size1, int size2, int size3);
OPENMM_EXPORT_AMOEBA void OpenMM_3D_DoubleArray_set(OpenMM_3D_DoubleArray* array, int index1, int index2, OpenMM_DoubleArray* values);
OPENMM_EXPORT_AMOEBA void OpenMM_3D_DoubleArray_destroy(OpenMM_3D_DoubleArray* array);

/* AmoebaMultipoleForce */
typedef enum {
  OpenMM_AmoebaMultipoleForce_NoCutoff = 0, OpenMM_AmoebaMultipoleForce_PME = 1
} OpenMM_AmoebaMultipoleForce_NonbondedMethod;
typedef enum {
  OpenMM_AmoebaMultipoleForce_Mutual = 0, OpenMM_AmoebaMultipoleForce_Direct = 1, OpenMM_AmoebaMultipoleForce_Extrapolated = 2
} OpenMM_AmoebaMultipoleForce_PolarizationType;
typedef enum {
  OpenMM_AmoebaMultipoleForce_ZThenX = 0, OpenMM_AmoebaMultipoleForce_Bisector = 1, OpenMM_AmoebaMultipoleForce_ZBisect = 2, OpenMM_AmoebaMultipoleForce_ThreeFold = 3, OpenMM_AmoebaMultipoleForce_ZOnly = 4, OpenMM_AmoebaMultipoleForce_NoAxisType = 5, OpenMM_AmoebaMultipoleForce_LastAxisTypeIndex = 6
} OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes;
typedef enum {
  OpenMM_AmoebaMultipoleForce_Covalent12 = 0, OpenMM_AmoebaMultipoleForce_Covalent13 = 1, OpenMM_AmoebaMultipoleForce_Covalent14 = 2, OpenMM_AmoebaMultipoleForce_Covalent15 = 3, OpenMM_AmoebaMultipoleForce_PolarizationCovalent11 = 4, OpenMM_AmoebaMultipoleForce_PolarizationCovalent12 = 5, OpenMM_AmoebaMultipoleForce_PolarizationCovalent13 = 6, OpenMM_AmoebaMultipoleForce_PolarizationCovalent14 = 7, OpenMM_AmoebaMultipoleForce_CovalentEnd = 8
} OpenMM_AmoebaMultipoleForce_CovalentType;

extern OPENMM_EXPORT_AMOEBA OpenMM_AmoebaMultipoleForce* OpenMM_AmoebaMultipoleForce_create();
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_destroy(OpenMM_AmoebaMultipoleForce* target);
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaMultipoleForce_getNumMultipoles(const OpenMM_AmoebaMultipoleForce* target);
extern OPENMM_EXPORT_AMOEBA OpenMM_AmoebaMultipoleForce_NonbondedMethod OpenMM_AmoebaMultipoleForce_getNonbondedMethod(const OpenMM_AmoebaMultipoleForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_setNonbondedMethod(OpenMM_AmoebaMultipoleForce* target, OpenMM_AmoebaMultipoleForce_NonbondedMethod method);
extern OPENMM_EXPORT_AMOEBA OpenMM_AmoebaMultipoleForce_PolarizationType OpenMM_AmoebaMultipoleForce_getPolarizationType(const OpenMM_AmoebaMultipoleForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_setPolarizationType(OpenMM_AmoebaMultipoleForce* target, OpenMM_AmoebaMultipoleForce_PolarizationType type);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaMultipoleForce_getCutoffDistance(const OpenMM_AmoebaMultipoleForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_setCutoffDistance(OpenMM_AmoebaMultipoleForce* target, double distance);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_getPMEParameters(const OpenMM_AmoebaMultipoleForce* target, double* alpha, int* nx, int* ny, int* nz);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_setPMEParameters(OpenMM_AmoebaMultipoleForce* target, double alpha, int nx, int ny, int nz);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaMultipoleForce_getAEwald(const OpenMM_AmoebaMultipoleForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_setAEwald(OpenMM_AmoebaMultipoleForce* target, double aewald);
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaMultipoleForce_getPmeBSplineOrder(const OpenMM_AmoebaMultipoleForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_getPmeGridDimensions(const OpenMM_AmoebaMultipoleForce* target, OpenMM_IntArray* gridDimension);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_setPmeGridDimensions(OpenMM_AmoebaMultipoleForce* target, const OpenMM_IntArray* gridDimension);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_getPMEParametersInContext(const OpenMM_AmoebaMultipoleForce* target, const OpenMM_Context* context, double* alpha, int* nx, int* ny, int* nz);
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaMultipoleForce_addMultipole(OpenMM_AmoebaMultipoleForce* target, double charge, const OpenMM_DoubleArray* molecularDipole, const OpenMM_DoubleArray* molecularQuadrupole, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY, double thole, double dampingFactor, double polarity);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_getMultipoleParameters(const OpenMM_AmoebaMultipoleForce* target, int index, double* charge, OpenMM_DoubleArray* molecularDipole, OpenMM_DoubleArray* molecularQuadrupole, int* axisType, int* multipoleAtomZ, int* multipoleAtomX, int* multipoleAtomY, double* thole, double* dampingFactor, double* polarity);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_setMultipoleParameters(OpenMM_AmoebaMultipoleForce* target, int index, double charge, const OpenMM_DoubleArray* molecularDipole, const OpenMM_DoubleArray* molecularQuadrupole, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY, double thole, double dampingFactor, double polarity);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_setCovalentMap(OpenMM_AmoebaMultipoleForce* target, int index, OpenMM_AmoebaMultipoleForce_CovalentType typeId, const OpenMM_IntArray* covalentAtoms);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_getCovalentMap(const OpenMM_AmoebaMultipoleForce* target, int index, OpenMM_AmoebaMultipoleForce_CovalentType typeId, OpenMM_IntArray* covalentAtoms);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_getCovalentMaps(const OpenMM_AmoebaMultipoleForce* target, int index, OpenMM_2D_IntArray* covalentLists);
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaMultipoleForce_getMutualInducedMaxIterations(const OpenMM_AmoebaMultipoleForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations(OpenMM_AmoebaMultipoleForce* target, int inputMutualInducedMaxIterations);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaMultipoleForce_getMutualInducedTargetEpsilon(const OpenMM_AmoebaMultipoleForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon(OpenMM_AmoebaMultipoleForce* target, double inputMutualInducedTargetEpsilon);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients(OpenMM_AmoebaMultipoleForce* target, const OpenMM_DoubleArray* coefficients);
extern OPENMM_EXPORT_AMOEBA const OpenMM_DoubleArray* OpenMM_AmoebaMultipoleForce_getExtrapolationCoefficients(const OpenMM_AmoebaMultipoleForce* target);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaMultipoleForce_getEwaldErrorTolerance(const OpenMM_AmoebaMultipoleForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance(OpenMM_AmoebaMultipoleForce* target, double tol);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_getLabFramePermanentDipoles(OpenMM_AmoebaMultipoleForce* target, OpenMM_Context* context, OpenMM_Vec3Array* dipoles);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_getInducedDipoles(OpenMM_AmoebaMultipoleForce* target, OpenMM_Context* context, OpenMM_Vec3Array* dipoles);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_getTotalDipoles(OpenMM_AmoebaMultipoleForce* target, OpenMM_Context* context, OpenMM_Vec3Array* dipoles);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_getElectrostaticPotential(OpenMM_AmoebaMultipoleForce* target, const OpenMM_Vec3Array* inputGrid, OpenMM_Context* context, OpenMM_DoubleArray* outputElectrostaticPotential);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_getSystemMultipoleMoments(OpenMM_AmoebaMultipoleForce* target, OpenMM_Context* context, OpenMM_DoubleArray* outputMultipoleMoments);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaMultipoleForce_updateParametersInContext(OpenMM_AmoebaMultipoleForce* target, OpenMM_Context* context);
extern OPENMM_EXPORT_AMOEBA OpenMM_Boolean OpenMM_AmoebaMultipoleForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaMultipoleForce* target);

/* AmoebaStretchBendForce */
extern OPENMM_EXPORT_AMOEBA OpenMM_AmoebaStretchBendForce* OpenMM_AmoebaStretchBendForce_create();
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaStretchBendForce_destroy(OpenMM_AmoebaStretchBendForce* target);
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaStretchBendForce_getNumStretchBends(const OpenMM_AmoebaStretchBendForce* target);
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaStretchBendForce_addStretchBend(OpenMM_AmoebaStretchBendForce* target, int particle1, int particle2, int particle3, double lengthAB, double lengthCB, double angle, double k1, double k2);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaStretchBendForce_getStretchBendParameters(const OpenMM_AmoebaStretchBendForce* target, int index, int* particle1, int* particle2, int* particle3, double* lengthAB, double* lengthCB, double* angle, double* k1, double* k2);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaStretchBendForce_setStretchBendParameters(OpenMM_AmoebaStretchBendForce* target, int index, int particle1, int particle2, int particle3, double lengthAB, double lengthCB, double angle, double k1, double k2);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaStretchBendForce_updateParametersInContext(OpenMM_AmoebaStretchBendForce* target, OpenMM_Context* context);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaStretchBendForce_setUsesPeriodicBoundaryConditions(OpenMM_AmoebaStretchBendForce* target, OpenMM_Boolean periodic);
extern OPENMM_EXPORT_AMOEBA OpenMM_Boolean OpenMM_AmoebaStretchBendForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaStretchBendForce* target);

/* AmoebaVdwForce */
typedef enum {
  OpenMM_AmoebaVdwForce_NoCutoff = 0, OpenMM_AmoebaVdwForce_CutoffPeriodic = 1
} OpenMM_AmoebaVdwForce_NonbondedMethod;
typedef enum {
  OpenMM_AmoebaVdwForce_None = 0, OpenMM_AmoebaVdwForce_Decouple = 1, OpenMM_AmoebaVdwForce_Annihilate = 2
} OpenMM_AmoebaVdwForce_AlchemicalMethod;

extern OPENMM_EXPORT_AMOEBA OpenMM_AmoebaVdwForce* OpenMM_AmoebaVdwForce_create();
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_destroy(OpenMM_AmoebaVdwForce* target);
extern OPENMM_EXPORT_AMOEBA const char* OpenMM_AmoebaVdwForce_Lambda();
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaVdwForce_getNumParticles(const OpenMM_AmoebaVdwForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_setParticleParameters(OpenMM_AmoebaVdwForce* target, int particleIndex, int parentIndex, double sigma, double epsilon, double reductionFactor, OpenMM_Boolean isAlchemical);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_getParticleParameters(const OpenMM_AmoebaVdwForce* target, int particleIndex, int* parentIndex, double* sigma, double* epsilon, double* reductionFactor, OpenMM_Boolean* isAlchemical);
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaVdwForce_addParticle(OpenMM_AmoebaVdwForce* target, int parentIndex, double sigma, double epsilon, double reductionFactor, OpenMM_Boolean isAlchemical);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_setSigmaCombiningRule(OpenMM_AmoebaVdwForce* target, const char* sigmaCombiningRule);
extern OPENMM_EXPORT_AMOEBA const char* OpenMM_AmoebaVdwForce_getSigmaCombiningRule(const OpenMM_AmoebaVdwForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_setEpsilonCombiningRule(OpenMM_AmoebaVdwForce* target, const char* epsilonCombiningRule);
extern OPENMM_EXPORT_AMOEBA const char* OpenMM_AmoebaVdwForce_getEpsilonCombiningRule(const OpenMM_AmoebaVdwForce* target);
extern OPENMM_EXPORT_AMOEBA OpenMM_Boolean OpenMM_AmoebaVdwForce_getUseDispersionCorrection(const OpenMM_AmoebaVdwForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_setUseDispersionCorrection(OpenMM_AmoebaVdwForce* target, OpenMM_Boolean useCorrection);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_setParticleExclusions(OpenMM_AmoebaVdwForce* target, int particleIndex, const OpenMM_IntArray* exclusions);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_getParticleExclusions(const OpenMM_AmoebaVdwForce* target, int particleIndex, OpenMM_IntArray* exclusions);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaVdwForce_getCutoffDistance(const OpenMM_AmoebaVdwForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_setCutoffDistance(OpenMM_AmoebaVdwForce* target, double distance);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_setCutoff(OpenMM_AmoebaVdwForce* target, double cutoff);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaVdwForce_getCutoff(const OpenMM_AmoebaVdwForce* target);
extern OPENMM_EXPORT_AMOEBA OpenMM_AmoebaVdwForce_NonbondedMethod OpenMM_AmoebaVdwForce_getNonbondedMethod(const OpenMM_AmoebaVdwForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_setNonbondedMethod(OpenMM_AmoebaVdwForce* target, OpenMM_AmoebaVdwForce_NonbondedMethod method);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_setSoftcorePower(OpenMM_AmoebaVdwForce* target, int n);
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaVdwForce_getSoftcorePower(const OpenMM_AmoebaVdwForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_setSoftcoreAlpha(OpenMM_AmoebaVdwForce* target, double alpha);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaVdwForce_getSoftcoreAlpha(const OpenMM_AmoebaVdwForce* target);
extern OPENMM_EXPORT_AMOEBA OpenMM_AmoebaVdwForce_AlchemicalMethod OpenMM_AmoebaVdwForce_getAlchemicalMethod(const OpenMM_AmoebaVdwForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_setAlchemicalMethod(OpenMM_AmoebaVdwForce* target, OpenMM_AmoebaVdwForce_AlchemicalMethod method);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaVdwForce_updateParametersInContext(OpenMM_AmoebaVdwForce* target, OpenMM_Context* context);
extern OPENMM_EXPORT_AMOEBA OpenMM_Boolean OpenMM_AmoebaVdwForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaVdwForce* target);

/* AmoebaGeneralizedKirkwoodForce */
extern OPENMM_EXPORT_AMOEBA OpenMM_AmoebaGeneralizedKirkwoodForce* OpenMM_AmoebaGeneralizedKirkwoodForce_create();
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaGeneralizedKirkwoodForce_destroy(OpenMM_AmoebaGeneralizedKirkwoodForce* target);
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaGeneralizedKirkwoodForce_getNumParticles(const OpenMM_AmoebaGeneralizedKirkwoodForce* target);
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle(OpenMM_AmoebaGeneralizedKirkwoodForce* target, double charge, double radius, double scalingFactor);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaGeneralizedKirkwoodForce_getParticleParameters(const OpenMM_AmoebaGeneralizedKirkwoodForce* target, int index, double* charge, double* radius, double* scalingFactor);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters(OpenMM_AmoebaGeneralizedKirkwoodForce* target, int index, double charge, double radius, double scalingFactor);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaGeneralizedKirkwoodForce_getSolventDielectric(const OpenMM_AmoebaGeneralizedKirkwoodForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric(OpenMM_AmoebaGeneralizedKirkwoodForce* target, double dielectric);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaGeneralizedKirkwoodForce_getSoluteDielectric(const OpenMM_AmoebaGeneralizedKirkwoodForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric(OpenMM_AmoebaGeneralizedKirkwoodForce* target, double dielectric);
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaGeneralizedKirkwoodForce_getIncludeCavityTerm(const OpenMM_AmoebaGeneralizedKirkwoodForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm(OpenMM_AmoebaGeneralizedKirkwoodForce* target, int includeCavityTerm);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaGeneralizedKirkwoodForce_getProbeRadius(const OpenMM_AmoebaGeneralizedKirkwoodForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius(OpenMM_AmoebaGeneralizedKirkwoodForce* target, double probeRadius);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaGeneralizedKirkwoodForce_getSurfaceAreaFactor(const OpenMM_AmoebaGeneralizedKirkwoodForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor(OpenMM_AmoebaGeneralizedKirkwoodForce* target, double surfaceAreaFactor);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext(OpenMM_AmoebaGeneralizedKirkwoodForce* target, OpenMM_Context* context);
extern OPENMM_EXPORT_AMOEBA OpenMM_Boolean OpenMM_AmoebaGeneralizedKirkwoodForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaGeneralizedKirkwoodForce* target);

/* AmoebaOutOfPlaneBendForce */
extern OPENMM_EXPORT_AMOEBA OpenMM_AmoebaOutOfPlaneBendForce* OpenMM_AmoebaOutOfPlaneBendForce_create();
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaOutOfPlaneBendForce_destroy(OpenMM_AmoebaOutOfPlaneBendForce* target);
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaOutOfPlaneBendForce_getNumOutOfPlaneBends(const OpenMM_AmoebaOutOfPlaneBendForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendCubic(OpenMM_AmoebaOutOfPlaneBendForce* target, double cubicK);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaOutOfPlaneBendForce_getAmoebaGlobalOutOfPlaneBendCubic(const OpenMM_AmoebaOutOfPlaneBendForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendQuartic(OpenMM_AmoebaOutOfPlaneBendForce* target, double quarticK);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaOutOfPlaneBendForce_getAmoebaGlobalOutOfPlaneBendQuartic(const OpenMM_AmoebaOutOfPlaneBendForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendPentic(OpenMM_AmoebaOutOfPlaneBendForce* target, double penticK);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaOutOfPlaneBendForce_getAmoebaGlobalOutOfPlaneBendPentic(const OpenMM_AmoebaOutOfPlaneBendForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendSextic(OpenMM_AmoebaOutOfPlaneBendForce* target, double sexticK);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaOutOfPlaneBendForce_getAmoebaGlobalOutOfPlaneBendSextic(const OpenMM_AmoebaOutOfPlaneBendForce* target);
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaOutOfPlaneBendForce_addOutOfPlaneBend(OpenMM_AmoebaOutOfPlaneBendForce* target, int particle1, int particle2, int particle3, int particle4, double k);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaOutOfPlaneBendForce_getOutOfPlaneBendParameters(const OpenMM_AmoebaOutOfPlaneBendForce* target, int index, int* particle1, int* particle2, int* particle3, int* particle4, double* k);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaOutOfPlaneBendForce_setOutOfPlaneBendParameters(OpenMM_AmoebaOutOfPlaneBendForce* target, int index, int particle1, int particle2, int particle3, int particle4, double k);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaOutOfPlaneBendForce_updateParametersInContext(OpenMM_AmoebaOutOfPlaneBendForce* target, OpenMM_Context* context);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaOutOfPlaneBendForce_setUsesPeriodicBoundaryConditions(OpenMM_AmoebaOutOfPlaneBendForce* target, OpenMM_Boolean periodic);
extern OPENMM_EXPORT_AMOEBA OpenMM_Boolean OpenMM_AmoebaOutOfPlaneBendForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaOutOfPlaneBendForce* target);

/* HippoNonbondedForce */
typedef enum {
  OpenMM_HippoNonbondedForce_NoCutoff = 0, OpenMM_HippoNonbondedForce_PME = 1
} OpenMM_HippoNonbondedForce_NonbondedMethod;
typedef enum {
  OpenMM_HippoNonbondedForce_ZThenX = 0, OpenMM_HippoNonbondedForce_Bisector = 1, OpenMM_HippoNonbondedForce_ZBisect = 2, OpenMM_HippoNonbondedForce_ThreeFold = 3, OpenMM_HippoNonbondedForce_ZOnly = 4, OpenMM_HippoNonbondedForce_NoAxisType = 5
} OpenMM_HippoNonbondedForce_ParticleAxisTypes;

extern OPENMM_EXPORT_AMOEBA OpenMM_HippoNonbondedForce* OpenMM_HippoNonbondedForce_create();
extern OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_destroy(OpenMM_HippoNonbondedForce* target);
extern OPENMM_EXPORT_AMOEBA int OpenMM_HippoNonbondedForce_getNumParticles(const OpenMM_HippoNonbondedForce* target);
extern OPENMM_EXPORT_AMOEBA int OpenMM_HippoNonbondedForce_getNumExceptions(const OpenMM_HippoNonbondedForce* target);
extern OPENMM_EXPORT_AMOEBA OpenMM_HippoNonbondedForce_NonbondedMethod OpenMM_HippoNonbondedForce_getNonbondedMethod(const OpenMM_HippoNonbondedForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_setNonbondedMethod(OpenMM_HippoNonbondedForce* target, OpenMM_HippoNonbondedForce_NonbondedMethod method);
extern OPENMM_EXPORT_AMOEBA double OpenMM_HippoNonbondedForce_getCutoffDistance(const OpenMM_HippoNonbondedForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_setCutoffDistance(OpenMM_HippoNonbondedForce* target, double distance);
extern OPENMM_EXPORT_AMOEBA double OpenMM_HippoNonbondedForce_getSwitchingDistance(const OpenMM_HippoNonbondedForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_setSwitchingDistance(OpenMM_HippoNonbondedForce* target, double distance);
extern OPENMM_EXPORT_AMOEBA const OpenMM_DoubleArray* OpenMM_HippoNonbondedForce_getExtrapolationCoefficients(const OpenMM_HippoNonbondedForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_setExtrapolationCoefficients(OpenMM_HippoNonbondedForce* target, const OpenMM_DoubleArray* coefficients);
extern OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_getPMEParameters(const OpenMM_HippoNonbondedForce* target, double* alpha, int* nx, int* ny, int* nz);
extern OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_getDPMEParameters(const OpenMM_HippoNonbondedForce* target, double* alpha, int* nx, int* ny, int* nz);
extern OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_setPMEParameters(OpenMM_HippoNonbondedForce* target, double alpha, int nx, int ny, int nz);
extern OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_setDPMEParameters(OpenMM_HippoNonbondedForce* target, double alpha, int nx, int ny, int nz);
extern OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_getPMEParametersInContext(const OpenMM_HippoNonbondedForce* target, const OpenMM_Context* context, double* alpha, int* nx, int* ny, int* nz);
extern OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_getDPMEParametersInContext(const OpenMM_HippoNonbondedForce* target, const OpenMM_Context* context, double* alpha, int* nx, int* ny, int* nz);
extern OPENMM_EXPORT_AMOEBA int OpenMM_HippoNonbondedForce_addParticle(OpenMM_HippoNonbondedForce* target, double charge, const OpenMM_DoubleArray* dipole, const OpenMM_DoubleArray* quadrupole, double coreCharge, double alpha, double epsilon, double damping, double c6, double pauliK, double pauliQ, double pauliAlpha, double polarizability, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY);
extern OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_getParticleParameters(const OpenMM_HippoNonbondedForce* target, int index, double* charge, OpenMM_DoubleArray* dipole, OpenMM_DoubleArray* quadrupole, double* coreCharge, double* alpha, double* epsilon, double* damping, double* c6, double* pauliK, double* pauliQ, double* pauliAlpha, double* polarizability, int* axisType, int* multipoleAtomZ, int* multipoleAtomX, int* multipoleAtomY);
extern OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_setParticleParameters(OpenMM_HippoNonbondedForce* target, int index, double charge, const OpenMM_DoubleArray* dipole, const OpenMM_DoubleArray* quadrupole, double coreCharge, double alpha, double epsilon, double damping, double c6, double pauliK, double pauliQ, double pauliAlpha, double polarizability, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY);
extern OPENMM_EXPORT_AMOEBA int OpenMM_HippoNonbondedForce_addException(OpenMM_HippoNonbondedForce* target, int particle1, int particle2, double multipoleMultipoleScale, double dipoleMultipoleScale, double dipoleDipoleScale, double dispersionScale, double repulsionScale, double chargeTransferScale, OpenMM_Boolean replace);
extern OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_getExceptionParameters(const OpenMM_HippoNonbondedForce* target, int index, int* particle1, int* particle2, double* multipoleMultipoleScale, double* dipoleMultipoleScale, double* dipoleDipoleScale, double* dispersionScale, double* repulsionScale, double* chargeTransferScale);
extern OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_setExceptionParameters(OpenMM_HippoNonbondedForce* target, int index, int particle1, int particle2, double multipoleMultipoleScale, double dipoleMultipoleScale, double dipoleDipoleScale, double dispersionScale, double repulsionScale, double chargeTransferScale);
extern OPENMM_EXPORT_AMOEBA double OpenMM_HippoNonbondedForce_getEwaldErrorTolerance(const OpenMM_HippoNonbondedForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_setEwaldErrorTolerance(OpenMM_HippoNonbondedForce* target, double tol);
extern OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_getLabFramePermanentDipoles(OpenMM_HippoNonbondedForce* target, OpenMM_Context* context, OpenMM_Vec3Array* dipoles);
extern OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_getInducedDipoles(OpenMM_HippoNonbondedForce* target, OpenMM_Context* context, OpenMM_Vec3Array* dipoles);
extern OPENMM_EXPORT_AMOEBA void OpenMM_HippoNonbondedForce_updateParametersInContext(OpenMM_HippoNonbondedForce* target, OpenMM_Context* context);
extern OPENMM_EXPORT_AMOEBA OpenMM_Boolean OpenMM_HippoNonbondedForce_usesPeriodicBoundaryConditions(const OpenMM_HippoNonbondedForce* target);

/* AmoebaAngleForce */
extern OPENMM_EXPORT_AMOEBA OpenMM_AmoebaAngleForce* OpenMM_AmoebaAngleForce_create();
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaAngleForce_destroy(OpenMM_AmoebaAngleForce* target);
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaAngleForce_getNumAngles(const OpenMM_AmoebaAngleForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleCubic(OpenMM_AmoebaAngleForce* target, double cubicK);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaAngleForce_getAmoebaGlobalAngleCubic(const OpenMM_AmoebaAngleForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleQuartic(OpenMM_AmoebaAngleForce* target, double quarticK);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaAngleForce_getAmoebaGlobalAngleQuartic(const OpenMM_AmoebaAngleForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaAngleForce_setAmoebaGlobalAnglePentic(OpenMM_AmoebaAngleForce* target, double penticK);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaAngleForce_getAmoebaGlobalAnglePentic(const OpenMM_AmoebaAngleForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleSextic(OpenMM_AmoebaAngleForce* target, double sexticK);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaAngleForce_getAmoebaGlobalAngleSextic(const OpenMM_AmoebaAngleForce* target);
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaAngleForce_addAngle(OpenMM_AmoebaAngleForce* target, int particle1, int particle2, int particle3, double length, double quadraticK);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaAngleForce_getAngleParameters(const OpenMM_AmoebaAngleForce* target, int index, int* particle1, int* particle2, int* particle3, double* length, double* quadraticK);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaAngleForce_setAngleParameters(OpenMM_AmoebaAngleForce* target, int index, int particle1, int particle2, int particle3, double length, double quadraticK);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaAngleForce_updateParametersInContext(OpenMM_AmoebaAngleForce* target, OpenMM_Context* context);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaAngleForce_setUsesPeriodicBoundaryConditions(OpenMM_AmoebaAngleForce* target, OpenMM_Boolean periodic);
extern OPENMM_EXPORT_AMOEBA OpenMM_Boolean OpenMM_AmoebaAngleForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaAngleForce* target);

/* AmoebaInPlaneAngleForce */
extern OPENMM_EXPORT_AMOEBA OpenMM_AmoebaInPlaneAngleForce* OpenMM_AmoebaInPlaneAngleForce_create();
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaInPlaneAngleForce_destroy(OpenMM_AmoebaInPlaneAngleForce* target);
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaInPlaneAngleForce_getNumAngles(const OpenMM_AmoebaInPlaneAngleForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleCubic(OpenMM_AmoebaInPlaneAngleForce* target, double cubicK);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaInPlaneAngleForce_getAmoebaGlobalInPlaneAngleCubic(const OpenMM_AmoebaInPlaneAngleForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleQuartic(OpenMM_AmoebaInPlaneAngleForce* target, double quarticK);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaInPlaneAngleForce_getAmoebaGlobalInPlaneAngleQuartic(const OpenMM_AmoebaInPlaneAngleForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAnglePentic(OpenMM_AmoebaInPlaneAngleForce* target, double penticK);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaInPlaneAngleForce_getAmoebaGlobalInPlaneAnglePentic(const OpenMM_AmoebaInPlaneAngleForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleSextic(OpenMM_AmoebaInPlaneAngleForce* target, double sexticK);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaInPlaneAngleForce_getAmoebaGlobalInPlaneAngleSextic(const OpenMM_AmoebaInPlaneAngleForce* target);
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaInPlaneAngleForce_addAngle(OpenMM_AmoebaInPlaneAngleForce* target, int particle1, int particle2, int particle3, int particle4, double length, double quadraticK);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaInPlaneAngleForce_getAngleParameters(const OpenMM_AmoebaInPlaneAngleForce* target, int index, int* particle1, int* particle2, int* particle3, int* particle4, double* length, double* quadraticK);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaInPlaneAngleForce_setAngleParameters(OpenMM_AmoebaInPlaneAngleForce* target, int index, int particle1, int particle2, int particle3, int particle4, double length, double quadraticK);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaInPlaneAngleForce_updateParametersInContext(OpenMM_AmoebaInPlaneAngleForce* target, OpenMM_Context* context);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaInPlaneAngleForce_setUsesPeriodicBoundaryConditions(OpenMM_AmoebaInPlaneAngleForce* target, OpenMM_Boolean periodic);
extern OPENMM_EXPORT_AMOEBA OpenMM_Boolean OpenMM_AmoebaInPlaneAngleForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaInPlaneAngleForce* target);

/* AmoebaBondForce */
extern OPENMM_EXPORT_AMOEBA OpenMM_AmoebaBondForce* OpenMM_AmoebaBondForce_create();
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaBondForce_destroy(OpenMM_AmoebaBondForce* target);
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaBondForce_getNumBonds(const OpenMM_AmoebaBondForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaBondForce_setAmoebaGlobalBondCubic(OpenMM_AmoebaBondForce* target, double cubicK);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaBondForce_getAmoebaGlobalBondCubic(const OpenMM_AmoebaBondForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaBondForce_setAmoebaGlobalBondQuartic(OpenMM_AmoebaBondForce* target, double quarticK);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaBondForce_getAmoebaGlobalBondQuartic(const OpenMM_AmoebaBondForce* target);
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaBondForce_addBond(OpenMM_AmoebaBondForce* target, int particle1, int particle2, double length, double quadraticK);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaBondForce_getBondParameters(const OpenMM_AmoebaBondForce* target, int index, int* particle1, int* particle2, double* length, double* quadraticK);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaBondForce_setBondParameters(OpenMM_AmoebaBondForce* target, int index, int particle1, int particle2, double length, double quadraticK);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaBondForce_updateParametersInContext(OpenMM_AmoebaBondForce* target, OpenMM_Context* context);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaBondForce_setUsesPeriodicBoundaryConditions(OpenMM_AmoebaBondForce* target, OpenMM_Boolean periodic);
extern OPENMM_EXPORT_AMOEBA OpenMM_Boolean OpenMM_AmoebaBondForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaBondForce* target);

/* AmoebaWcaDispersionForce */
extern OPENMM_EXPORT_AMOEBA OpenMM_AmoebaWcaDispersionForce* OpenMM_AmoebaWcaDispersionForce_create();
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaWcaDispersionForce_destroy(OpenMM_AmoebaWcaDispersionForce* target);
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaWcaDispersionForce_getNumParticles(const OpenMM_AmoebaWcaDispersionForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaWcaDispersionForce_setParticleParameters(OpenMM_AmoebaWcaDispersionForce* target, int particleIndex, double radius, double epsilon);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaWcaDispersionForce_getParticleParameters(const OpenMM_AmoebaWcaDispersionForce* target, int particleIndex, double* radius, double* epsilon);
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaWcaDispersionForce_addParticle(OpenMM_AmoebaWcaDispersionForce* target, double radius, double epsilon);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaWcaDispersionForce_updateParametersInContext(OpenMM_AmoebaWcaDispersionForce* target, OpenMM_Context* context);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaWcaDispersionForce_getEpso(const OpenMM_AmoebaWcaDispersionForce* target);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaWcaDispersionForce_getEpsh(const OpenMM_AmoebaWcaDispersionForce* target);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaWcaDispersionForce_getRmino(const OpenMM_AmoebaWcaDispersionForce* target);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaWcaDispersionForce_getRminh(const OpenMM_AmoebaWcaDispersionForce* target);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaWcaDispersionForce_getAwater(const OpenMM_AmoebaWcaDispersionForce* target);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaWcaDispersionForce_getShctd(const OpenMM_AmoebaWcaDispersionForce* target);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaWcaDispersionForce_getDispoff(const OpenMM_AmoebaWcaDispersionForce* target);
extern OPENMM_EXPORT_AMOEBA double OpenMM_AmoebaWcaDispersionForce_getSlevy(const OpenMM_AmoebaWcaDispersionForce* target);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaWcaDispersionForce_setEpso(OpenMM_AmoebaWcaDispersionForce* target, double inputValue);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaWcaDispersionForce_setEpsh(OpenMM_AmoebaWcaDispersionForce* target, double inputValue);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaWcaDispersionForce_setRmino(OpenMM_AmoebaWcaDispersionForce* target, double inputValue);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaWcaDispersionForce_setRminh(OpenMM_AmoebaWcaDispersionForce* target, double inputValue);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaWcaDispersionForce_setAwater(OpenMM_AmoebaWcaDispersionForce* target, double inputValue);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaWcaDispersionForce_setShctd(OpenMM_AmoebaWcaDispersionForce* target, double inputValue);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaWcaDispersionForce_setDispoff(OpenMM_AmoebaWcaDispersionForce* target, double inputValue);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaWcaDispersionForce_setSlevy(OpenMM_AmoebaWcaDispersionForce* target, double inputValue);
extern OPENMM_EXPORT_AMOEBA OpenMM_Boolean OpenMM_AmoebaWcaDispersionForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaWcaDispersionForce* target);

/* AmoebaPiTorsionForce */
extern OPENMM_EXPORT_AMOEBA OpenMM_AmoebaPiTorsionForce* OpenMM_AmoebaPiTorsionForce_create();
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaPiTorsionForce_destroy(OpenMM_AmoebaPiTorsionForce* target);
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaPiTorsionForce_getNumPiTorsions(const OpenMM_AmoebaPiTorsionForce* target);
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaPiTorsionForce_addPiTorsion(OpenMM_AmoebaPiTorsionForce* target, int particle1, int particle2, int particle3, int particle4, int particle5, int particle6, double k);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaPiTorsionForce_getPiTorsionParameters(const OpenMM_AmoebaPiTorsionForce* target, int index, int* particle1, int* particle2, int* particle3, int* particle4, int* particle5, int* particle6, double* k);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaPiTorsionForce_setPiTorsionParameters(OpenMM_AmoebaPiTorsionForce* target, int index, int particle1, int particle2, int particle3, int particle4, int particle5, int particle6, double k);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaPiTorsionForce_updateParametersInContext(OpenMM_AmoebaPiTorsionForce* target, OpenMM_Context* context);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaPiTorsionForce_setUsesPeriodicBoundaryConditions(OpenMM_AmoebaPiTorsionForce* target, OpenMM_Boolean periodic);
extern OPENMM_EXPORT_AMOEBA OpenMM_Boolean OpenMM_AmoebaPiTorsionForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaPiTorsionForce* target);

/* AmoebaTorsionTorsionForce */
extern OPENMM_EXPORT_AMOEBA OpenMM_AmoebaTorsionTorsionForce* OpenMM_AmoebaTorsionTorsionForce_create();
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaTorsionTorsionForce_destroy(OpenMM_AmoebaTorsionTorsionForce* target);
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsions(const OpenMM_AmoebaTorsionTorsionForce* target);
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsionGrids(const OpenMM_AmoebaTorsionTorsionForce* target);
extern OPENMM_EXPORT_AMOEBA int OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion(OpenMM_AmoebaTorsionTorsionForce* target, int particle1, int particle2, int particle3, int particle4, int particle5, int chiralCheckAtomIndex, int gridIndex);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionParameters(const OpenMM_AmoebaTorsionTorsionForce* target, int index, int* particle1, int* particle2, int* particle3, int* particle4, int* particle5, int* chiralCheckAtomIndex, int* gridIndex);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionParameters(OpenMM_AmoebaTorsionTorsionForce* target, int index, int particle1, int particle2, int particle3, int particle4, int particle5, int chiralCheckAtomIndex, int gridIndex);
extern OPENMM_EXPORT_AMOEBA const OpenMM_3D_DoubleArray* OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionGrid(const OpenMM_AmoebaTorsionTorsionForce* target, int index);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid(OpenMM_AmoebaTorsionTorsionForce* target, int index, const OpenMM_3D_DoubleArray* grid);
extern OPENMM_EXPORT_AMOEBA void OpenMM_AmoebaTorsionTorsionForce_setUsesPeriodicBoundaryConditions(OpenMM_AmoebaTorsionTorsionForce* target, OpenMM_Boolean periodic);
extern OPENMM_EXPORT_AMOEBA OpenMM_Boolean OpenMM_AmoebaTorsionTorsionForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaTorsionTorsionForce* target);


#if defined(__cplusplus)
}
#endif

#endif /*AMOEBA_OPENMM_CWRAPPER_H_*/
