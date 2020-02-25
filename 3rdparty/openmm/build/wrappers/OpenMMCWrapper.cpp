
#include "OpenMM.h"
#include "OpenMMCWrapper.h"
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <vector>

using namespace OpenMM;
using namespace std;

extern "C" {

/* OpenMM_Vec3 */
OPENMM_EXPORT OpenMM_Vec3 OpenMM_Vec3_scale(const OpenMM_Vec3 vec, double scale) {
    OpenMM_Vec3 result = {vec.x*scale, vec.y*scale, vec.z*scale};
    return result;
}

/* OpenMM_Vec3Array */
OPENMM_EXPORT OpenMM_Vec3Array* OpenMM_Vec3Array_create(int size) {
    return reinterpret_cast<OpenMM_Vec3Array*>(new vector<Vec3>(size));
}
OPENMM_EXPORT void OpenMM_Vec3Array_destroy(OpenMM_Vec3Array* array) {
    delete reinterpret_cast<vector<Vec3>*>(array);
}
OPENMM_EXPORT int OpenMM_Vec3Array_getSize(const OpenMM_Vec3Array* array) {
    return reinterpret_cast<const vector<Vec3>*>(array)->size();
}
OPENMM_EXPORT void OpenMM_Vec3Array_resize(OpenMM_Vec3Array* array, int size) {
    reinterpret_cast<vector<Vec3>*>(array)->resize(size);
}
OPENMM_EXPORT void OpenMM_Vec3Array_append(OpenMM_Vec3Array* array, const OpenMM_Vec3 vec) {
    reinterpret_cast<vector<Vec3>*>(array)->push_back(Vec3(vec.x, vec.y, vec.z));
}
OPENMM_EXPORT void OpenMM_Vec3Array_set(OpenMM_Vec3Array* array, int index, const OpenMM_Vec3 vec) {
    (*reinterpret_cast<vector<Vec3>*>(array))[index] = Vec3(vec.x, vec.y, vec.z);
}
OPENMM_EXPORT const OpenMM_Vec3* OpenMM_Vec3Array_get(const OpenMM_Vec3Array* array, int index) {
    return reinterpret_cast<const OpenMM_Vec3*>((&(*reinterpret_cast<const vector<Vec3>*>(array))[index]));
}

/* OpenMM_StringArray */
OPENMM_EXPORT OpenMM_StringArray* OpenMM_StringArray_create(int size) {
    return reinterpret_cast<OpenMM_StringArray*>(new vector<string>(size));
}
OPENMM_EXPORT void OpenMM_StringArray_destroy(OpenMM_StringArray* array) {
    delete reinterpret_cast<vector<string>*>(array);
}
OPENMM_EXPORT int OpenMM_StringArray_getSize(const OpenMM_StringArray* array) {
    return reinterpret_cast<const vector<string>*>(array)->size();
}
OPENMM_EXPORT void OpenMM_StringArray_resize(OpenMM_StringArray* array, int size) {
    reinterpret_cast<vector<string>*>(array)->resize(size);
}
OPENMM_EXPORT void OpenMM_StringArray_append(OpenMM_StringArray* array, const char* str) {
    reinterpret_cast<vector<string>*>(array)->push_back(string(str));
}
OPENMM_EXPORT void OpenMM_StringArray_set(OpenMM_StringArray* array, int index, const char* str) {
    (*reinterpret_cast<vector<string>*>(array))[index] = string(str);
}
OPENMM_EXPORT const char* OpenMM_StringArray_get(const OpenMM_StringArray* array, int index) {
    return (*reinterpret_cast<const vector<string>*>(array))[index].c_str();
}

/* OpenMM_BondArray */
OPENMM_EXPORT OpenMM_BondArray* OpenMM_BondArray_create(int size) {
    return reinterpret_cast<OpenMM_BondArray*>(new vector<pair<int, int> >(size));
}
OPENMM_EXPORT void OpenMM_BondArray_destroy(OpenMM_BondArray* array) {
    delete reinterpret_cast<vector<pair<int, int> >*>(array);
}
OPENMM_EXPORT int OpenMM_BondArray_getSize(const OpenMM_BondArray* array) {
    return reinterpret_cast<const vector<pair<int, int> >*>(array)->size();
}
OPENMM_EXPORT void OpenMM_BondArray_resize(OpenMM_BondArray* array, int size) {
    reinterpret_cast<vector<pair<int, int> >*>(array)->resize(size);
}
OPENMM_EXPORT void OpenMM_BondArray_append(OpenMM_BondArray* array, int particle1, int particle2) {
    reinterpret_cast<vector<pair<int, int> >*>(array)->push_back(pair<int, int>(particle1, particle2));
}
OPENMM_EXPORT void OpenMM_BondArray_set(OpenMM_BondArray* array, int index, int particle1, int particle2) {
    (*reinterpret_cast<vector<pair<int, int> >*>(array))[index] = pair<int, int>(particle1, particle2);
}
OPENMM_EXPORT void OpenMM_BondArray_get(const OpenMM_BondArray* array, int index, int* particle1, int* particle2) {
    pair<int, int> particles = (*reinterpret_cast<const vector<pair<int, int> >*>(array))[index];
    *particle1 = particles.first;
    *particle2 = particles.second;
}

/* OpenMM_ParameterArray */
OPENMM_EXPORT int OpenMM_ParameterArray_getSize(const OpenMM_ParameterArray* array) {
    return reinterpret_cast<const map<string, double>*>(array)->size();
}
OPENMM_EXPORT double OpenMM_ParameterArray_get(const OpenMM_ParameterArray* array, const char* name) {
    const map<string, double>* params = reinterpret_cast<const map<string, double>*>(array);
    const map<string, double>::const_iterator iter = params->find(string(name));
    if (iter == params->end())
        throw OpenMMException("OpenMM_ParameterArray_get: No such parameter");
    return iter->second;
}

/* OpenMM_PropertyArray */
OPENMM_EXPORT int OpenMM_PropertyArray_getSize(const OpenMM_PropertyArray* array) {
    return reinterpret_cast<const map<string, double>*>(array)->size();
}
OPENMM_EXPORT const char* OpenMM_PropertyArray_get(const OpenMM_PropertyArray* array, const char* name) {
    const map<string, string>* params = reinterpret_cast<const map<string, string>*>(array);
    const map<string, string>::const_iterator iter = params->find(string(name));
    if (iter == params->end())
        throw OpenMMException("OpenMM_PropertyArray_get: No such property");
    return iter->second.c_str();
}

/* OpenMM_DoubleArray */
OPENMM_EXPORT OpenMM_DoubleArray* OpenMM_DoubleArray_create(int size) {
    return reinterpret_cast<OpenMM_DoubleArray*>(new vector<double>(size));
}
OPENMM_EXPORT void OpenMM_DoubleArray_destroy(OpenMM_DoubleArray* array) {
    delete reinterpret_cast<vector<double>*>(array);
}
OPENMM_EXPORT int OpenMM_DoubleArray_getSize(const OpenMM_DoubleArray* array) {
    return reinterpret_cast<const vector<double>*>(array)->size();
}
OPENMM_EXPORT void OpenMM_DoubleArray_resize(OpenMM_DoubleArray* array, int size) {
    reinterpret_cast<vector<double>*>(array)->resize(size);
}
OPENMM_EXPORT void OpenMM_DoubleArray_append(OpenMM_DoubleArray* array, double value) {
    reinterpret_cast<vector<double>*>(array)->push_back(value);
}
OPENMM_EXPORT void OpenMM_DoubleArray_set(OpenMM_DoubleArray* array, int index, double value) {
    (*reinterpret_cast<vector<double>*>(array))[index] = value;
}
OPENMM_EXPORT double OpenMM_DoubleArray_get(const OpenMM_DoubleArray* array, int index) {
    return (*reinterpret_cast<const vector<double>*>(array))[index];
}

/* OpenMM_IntArray */
OPENMM_EXPORT OpenMM_IntArray* OpenMM_IntArray_create(int size) {
    return reinterpret_cast<OpenMM_IntArray*>(new vector<int>(size));
}
OPENMM_EXPORT void OpenMM_IntArray_destroy(OpenMM_IntArray* array) {
    delete reinterpret_cast<vector<int>*>(array);
}
OPENMM_EXPORT int OpenMM_IntArray_getSize(const OpenMM_IntArray* array) {
    return reinterpret_cast<const vector<int>*>(array)->size();
}
OPENMM_EXPORT void OpenMM_IntArray_resize(OpenMM_IntArray* array, int size) {
    reinterpret_cast<vector<int>*>(array)->resize(size);
}
OPENMM_EXPORT void OpenMM_IntArray_append(OpenMM_IntArray* array, int value) {
    reinterpret_cast<vector<int>*>(array)->push_back(value);
}
OPENMM_EXPORT void OpenMM_IntArray_set(OpenMM_IntArray* array, int index, int value) {
    (*reinterpret_cast<vector<int>*>(array))[index] = value;
}
OPENMM_EXPORT int OpenMM_IntArray_get(const OpenMM_IntArray* array, int index) {
    return (*reinterpret_cast<const vector<int>*>(array))[index];
}

/* OpenMM_IntSet */
OPENMM_EXPORT OpenMM_IntSet* OpenMM_IntSet_create() {
    return reinterpret_cast<OpenMM_IntSet*>(new set<int>());
}
OPENMM_EXPORT void OpenMM_IntSet_destroy(OpenMM_IntSet* s) {
    delete reinterpret_cast<set<int>*>(s);
}
OPENMM_EXPORT int OpenMM_IntSet_getSize(const OpenMM_IntSet* s) {
    return reinterpret_cast<const set<int>*>(s)->size();
}
OPENMM_EXPORT void OpenMM_IntSet_insert(OpenMM_IntSet* s, int value) {
    reinterpret_cast<set<int>*>(s)->insert(value);
}

/* These methods need to be handled specially, since their C++ APIs cannot be directly translated to C.
   Unlike the C++ versions, the return value is allocated on the heap, and you must delete it yourself. */
OPENMM_EXPORT OpenMM_State* OpenMM_Context_getState(const OpenMM_Context* target, int types, int enforcePeriodicBox) {
    State result = reinterpret_cast<const Context*>(target)->getState(types, enforcePeriodicBox);
    return reinterpret_cast<OpenMM_State*>(new State(result));
}
OPENMM_EXPORT OpenMM_State* OpenMM_Context_getState_2(const OpenMM_Context* target, int types, int enforcePeriodicBox, int groups) {
    State result = reinterpret_cast<const Context*>(target)->getState(types, enforcePeriodicBox, groups);
    return reinterpret_cast<OpenMM_State*>(new State(result));
}
OPENMM_EXPORT OpenMM_StringArray* OpenMM_Platform_loadPluginsFromDirectory(const char* directory) {
    vector<string> result = Platform::loadPluginsFromDirectory(string(directory));
    return reinterpret_cast<OpenMM_StringArray*>(new vector<string>(result));
}
OPENMM_EXPORT OpenMM_StringArray* OpenMM_Platform_getPluginLoadFailures() {
    vector<string> result = Platform::getPluginLoadFailures();
    return reinterpret_cast<OpenMM_StringArray*>(new vector<string>(result));
}
static char* createStringFromStream(stringstream& stream) {
    int length = stream.str().size();
    char* result = (char*) malloc(length+1);
    stream.str().copy(result, length);
    result[length] = 0;
    return result;
}
OPENMM_EXPORT char* OpenMM_XmlSerializer_serializeSystem(const OpenMM_System* system) {
    stringstream stream;
    OpenMM::XmlSerializer::serialize<OpenMM::System>(reinterpret_cast<const OpenMM::System*>(system), "System", stream);
    return createStringFromStream(stream);
}
OPENMM_EXPORT char* OpenMM_XmlSerializer_serializeState(const OpenMM_State* state) {
    stringstream stream;
    OpenMM::XmlSerializer::serialize<OpenMM::State>(reinterpret_cast<const OpenMM::State*>(state), "State", stream);
    return createStringFromStream(stream);
}
OPENMM_EXPORT char* OpenMM_XmlSerializer_serializeIntegrator(const OpenMM_Integrator* integrator) {
    stringstream stream;
    OpenMM::XmlSerializer::serialize<OpenMM::Integrator>(reinterpret_cast<const OpenMM::Integrator*>(integrator), "Integrator", stream);
    return createStringFromStream(stream);
}
OPENMM_EXPORT OpenMM_System* OpenMM_XmlSerializer_deserializeSystem(const char* xml) {
    string input(xml);
    stringstream stream(input);
    return reinterpret_cast<OpenMM_System*>(OpenMM::XmlSerializer::deserialize<OpenMM::System>(stream));
}
OPENMM_EXPORT OpenMM_State* OpenMM_XmlSerializer_deserializeState(const char* xml) {
    string input(xml);
    stringstream stream(input);
    return reinterpret_cast<OpenMM_State*>(OpenMM::XmlSerializer::deserialize<OpenMM::State>(stream));
}
OPENMM_EXPORT OpenMM_Integrator* OpenMM_XmlSerializer_deserializeIntegrator(const char* xml) {
    string input(xml);
    stringstream stream(input);
    return reinterpret_cast<OpenMM_Integrator*>(OpenMM::XmlSerializer::deserialize<OpenMM::Integrator>(stream));
}

/* OpenMM::Force */
OPENMM_EXPORT void OpenMM_Force_destroy(OpenMM_Force* target) {
    delete reinterpret_cast<OpenMM::Force*>(target);
}
OPENMM_EXPORT int OpenMM_Force_getForceGroup(const OpenMM_Force* target) {
    int result = reinterpret_cast<const OpenMM::Force*>(target)->getForceGroup();
    return result;
}
OPENMM_EXPORT void OpenMM_Force_setForceGroup(OpenMM_Force* target, int group) {
    reinterpret_cast<OpenMM::Force*>(target)->setForceGroup(group);
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_Force_usesPeriodicBoundaryConditions(const OpenMM_Force* target) {
    bool result = reinterpret_cast<const OpenMM::Force*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::CustomAngleForce */
OPENMM_EXPORT OpenMM_CustomAngleForce* OpenMM_CustomAngleForce_create(const char* energy) {
    return reinterpret_cast<OpenMM_CustomAngleForce*>(new OpenMM::CustomAngleForce(std::string(energy)));
}
OPENMM_EXPORT void OpenMM_CustomAngleForce_destroy(OpenMM_CustomAngleForce* target) {
    delete reinterpret_cast<OpenMM::CustomAngleForce*>(target);
}
OPENMM_EXPORT int OpenMM_CustomAngleForce_getNumAngles(const OpenMM_CustomAngleForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomAngleForce*>(target)->getNumAngles();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomAngleForce_getNumPerAngleParameters(const OpenMM_CustomAngleForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomAngleForce*>(target)->getNumPerAngleParameters();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomAngleForce_getNumGlobalParameters(const OpenMM_CustomAngleForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomAngleForce*>(target)->getNumGlobalParameters();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomAngleForce_getNumEnergyParameterDerivatives(const OpenMM_CustomAngleForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomAngleForce*>(target)->getNumEnergyParameterDerivatives();
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomAngleForce_getEnergyFunction(const OpenMM_CustomAngleForce* target) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomAngleForce*>(target)->getEnergyFunction();
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomAngleForce_setEnergyFunction(OpenMM_CustomAngleForce* target, const char* energy) {
    reinterpret_cast<OpenMM::CustomAngleForce*>(target)->setEnergyFunction(std::string(energy));
}
OPENMM_EXPORT int OpenMM_CustomAngleForce_addPerAngleParameter(OpenMM_CustomAngleForce* target, const char* name) {
    int result = reinterpret_cast<OpenMM::CustomAngleForce*>(target)->addPerAngleParameter(std::string(name));
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomAngleForce_getPerAngleParameterName(const OpenMM_CustomAngleForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomAngleForce*>(target)->getPerAngleParameterName(index);
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomAngleForce_setPerAngleParameterName(OpenMM_CustomAngleForce* target, int index, const char* name) {
    reinterpret_cast<OpenMM::CustomAngleForce*>(target)->setPerAngleParameterName(index, std::string(name));
}
OPENMM_EXPORT int OpenMM_CustomAngleForce_addGlobalParameter(OpenMM_CustomAngleForce* target, const char* name, double defaultValue) {
    int result = reinterpret_cast<OpenMM::CustomAngleForce*>(target)->addGlobalParameter(std::string(name), defaultValue);
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomAngleForce_getGlobalParameterName(const OpenMM_CustomAngleForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomAngleForce*>(target)->getGlobalParameterName(index);
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomAngleForce_setGlobalParameterName(OpenMM_CustomAngleForce* target, int index, const char* name) {
    reinterpret_cast<OpenMM::CustomAngleForce*>(target)->setGlobalParameterName(index, std::string(name));
}
OPENMM_EXPORT double OpenMM_CustomAngleForce_getGlobalParameterDefaultValue(const OpenMM_CustomAngleForce* target, int index) {
    double result = reinterpret_cast<const OpenMM::CustomAngleForce*>(target)->getGlobalParameterDefaultValue(index);
    return result;
}
OPENMM_EXPORT void OpenMM_CustomAngleForce_setGlobalParameterDefaultValue(OpenMM_CustomAngleForce* target, int index, double defaultValue) {
    reinterpret_cast<OpenMM::CustomAngleForce*>(target)->setGlobalParameterDefaultValue(index, defaultValue);
}
OPENMM_EXPORT void OpenMM_CustomAngleForce_addEnergyParameterDerivative(OpenMM_CustomAngleForce* target, const char* name) {
    reinterpret_cast<OpenMM::CustomAngleForce*>(target)->addEnergyParameterDerivative(std::string(name));
}
OPENMM_EXPORT const char* OpenMM_CustomAngleForce_getEnergyParameterDerivativeName(const OpenMM_CustomAngleForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomAngleForce*>(target)->getEnergyParameterDerivativeName(index);
    return result->c_str();
}
OPENMM_EXPORT int OpenMM_CustomAngleForce_addAngle(OpenMM_CustomAngleForce* target, int particle1, int particle2, int particle3, const OpenMM_DoubleArray* parameters) {
    int result = reinterpret_cast<OpenMM::CustomAngleForce*>(target)->addAngle(particle1, particle2, particle3, *reinterpret_cast<const std::vector< double >*>(parameters));
    return result;
}
OPENMM_EXPORT void OpenMM_CustomAngleForce_getAngleParameters(const OpenMM_CustomAngleForce* target, int index, int* particle1, int* particle2, int* particle3, OpenMM_DoubleArray* parameters) {
    reinterpret_cast<const OpenMM::CustomAngleForce*>(target)->getAngleParameters(index, *reinterpret_cast<int*>(particle1), *reinterpret_cast<int*>(particle2), *reinterpret_cast<int*>(particle3), *reinterpret_cast<std::vector< double >*>(parameters));
}
OPENMM_EXPORT void OpenMM_CustomAngleForce_setAngleParameters(OpenMM_CustomAngleForce* target, int index, int particle1, int particle2, int particle3, const OpenMM_DoubleArray* parameters) {
    reinterpret_cast<OpenMM::CustomAngleForce*>(target)->setAngleParameters(index, particle1, particle2, particle3, *reinterpret_cast<const std::vector< double >*>(parameters));
}
OPENMM_EXPORT void OpenMM_CustomAngleForce_updateParametersInContext(OpenMM_CustomAngleForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::CustomAngleForce*>(target)->updateParametersInContext(*reinterpret_cast<OpenMM::Context*>(context));
}
OPENMM_EXPORT void OpenMM_CustomAngleForce_setUsesPeriodicBoundaryConditions(OpenMM_CustomAngleForce* target, OpenMM_Boolean periodic) {
    reinterpret_cast<OpenMM::CustomAngleForce*>(target)->setUsesPeriodicBoundaryConditions(periodic);
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_CustomAngleForce_usesPeriodicBoundaryConditions(const OpenMM_CustomAngleForce* target) {
    bool result = reinterpret_cast<const OpenMM::CustomAngleForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::MonteCarloBarostat */
OPENMM_EXPORT OpenMM_MonteCarloBarostat* OpenMM_MonteCarloBarostat_create(double defaultPressure, double defaultTemperature, int frequency) {
    return reinterpret_cast<OpenMM_MonteCarloBarostat*>(new OpenMM::MonteCarloBarostat(defaultPressure, defaultTemperature, frequency));
}
OPENMM_EXPORT void OpenMM_MonteCarloBarostat_destroy(OpenMM_MonteCarloBarostat* target) {
    delete reinterpret_cast<OpenMM::MonteCarloBarostat*>(target);
}
OPENMM_EXPORT const char* OpenMM_MonteCarloBarostat_Pressure() {
    const std::string* result = &OpenMM::MonteCarloBarostat::Pressure();
    return result->c_str();
}
OPENMM_EXPORT const char* OpenMM_MonteCarloBarostat_Temperature() {
    const std::string* result = &OpenMM::MonteCarloBarostat::Temperature();
    return result->c_str();
}
OPENMM_EXPORT double OpenMM_MonteCarloBarostat_getDefaultPressure(const OpenMM_MonteCarloBarostat* target) {
    double result = reinterpret_cast<const OpenMM::MonteCarloBarostat*>(target)->getDefaultPressure();
    return result;
}
OPENMM_EXPORT void OpenMM_MonteCarloBarostat_setDefaultPressure(OpenMM_MonteCarloBarostat* target, double pressure) {
    reinterpret_cast<OpenMM::MonteCarloBarostat*>(target)->setDefaultPressure(pressure);
}
OPENMM_EXPORT int OpenMM_MonteCarloBarostat_getFrequency(const OpenMM_MonteCarloBarostat* target) {
    int result = reinterpret_cast<const OpenMM::MonteCarloBarostat*>(target)->getFrequency();
    return result;
}
OPENMM_EXPORT void OpenMM_MonteCarloBarostat_setFrequency(OpenMM_MonteCarloBarostat* target, int freq) {
    reinterpret_cast<OpenMM::MonteCarloBarostat*>(target)->setFrequency(freq);
}
OPENMM_EXPORT double OpenMM_MonteCarloBarostat_getDefaultTemperature(const OpenMM_MonteCarloBarostat* target) {
    double result = reinterpret_cast<const OpenMM::MonteCarloBarostat*>(target)->getDefaultTemperature();
    return result;
}
OPENMM_EXPORT void OpenMM_MonteCarloBarostat_setDefaultTemperature(OpenMM_MonteCarloBarostat* target, double temp) {
    reinterpret_cast<OpenMM::MonteCarloBarostat*>(target)->setDefaultTemperature(temp);
}
OPENMM_EXPORT int OpenMM_MonteCarloBarostat_getRandomNumberSeed(const OpenMM_MonteCarloBarostat* target) {
    int result = reinterpret_cast<const OpenMM::MonteCarloBarostat*>(target)->getRandomNumberSeed();
    return result;
}
OPENMM_EXPORT void OpenMM_MonteCarloBarostat_setRandomNumberSeed(OpenMM_MonteCarloBarostat* target, int seed) {
    reinterpret_cast<OpenMM::MonteCarloBarostat*>(target)->setRandomNumberSeed(seed);
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_MonteCarloBarostat_usesPeriodicBoundaryConditions(const OpenMM_MonteCarloBarostat* target) {
    bool result = reinterpret_cast<const OpenMM::MonteCarloBarostat*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::RBTorsionForce */
OPENMM_EXPORT OpenMM_RBTorsionForce* OpenMM_RBTorsionForce_create() {
    return reinterpret_cast<OpenMM_RBTorsionForce*>(new OpenMM::RBTorsionForce());
}
OPENMM_EXPORT void OpenMM_RBTorsionForce_destroy(OpenMM_RBTorsionForce* target) {
    delete reinterpret_cast<OpenMM::RBTorsionForce*>(target);
}
OPENMM_EXPORT int OpenMM_RBTorsionForce_getNumTorsions(const OpenMM_RBTorsionForce* target) {
    int result = reinterpret_cast<const OpenMM::RBTorsionForce*>(target)->getNumTorsions();
    return result;
}
OPENMM_EXPORT int OpenMM_RBTorsionForce_addTorsion(OpenMM_RBTorsionForce* target, int particle1, int particle2, int particle3, int particle4, double c0, double c1, double c2, double c3, double c4, double c5) {
    int result = reinterpret_cast<OpenMM::RBTorsionForce*>(target)->addTorsion(particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
    return result;
}
OPENMM_EXPORT void OpenMM_RBTorsionForce_getTorsionParameters(const OpenMM_RBTorsionForce* target, int index, int* particle1, int* particle2, int* particle3, int* particle4, double* c0, double* c1, double* c2, double* c3, double* c4, double* c5) {
    reinterpret_cast<const OpenMM::RBTorsionForce*>(target)->getTorsionParameters(index, *reinterpret_cast<int*>(particle1), *reinterpret_cast<int*>(particle2), *reinterpret_cast<int*>(particle3), *reinterpret_cast<int*>(particle4), *reinterpret_cast<double*>(c0), *reinterpret_cast<double*>(c1), *reinterpret_cast<double*>(c2), *reinterpret_cast<double*>(c3), *reinterpret_cast<double*>(c4), *reinterpret_cast<double*>(c5));
}
OPENMM_EXPORT void OpenMM_RBTorsionForce_setTorsionParameters(OpenMM_RBTorsionForce* target, int index, int particle1, int particle2, int particle3, int particle4, double c0, double c1, double c2, double c3, double c4, double c5) {
    reinterpret_cast<OpenMM::RBTorsionForce*>(target)->setTorsionParameters(index, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
}
OPENMM_EXPORT void OpenMM_RBTorsionForce_updateParametersInContext(OpenMM_RBTorsionForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::RBTorsionForce*>(target)->updateParametersInContext(*reinterpret_cast<OpenMM::Context*>(context));
}
OPENMM_EXPORT void OpenMM_RBTorsionForce_setUsesPeriodicBoundaryConditions(OpenMM_RBTorsionForce* target, OpenMM_Boolean periodic) {
    reinterpret_cast<OpenMM::RBTorsionForce*>(target)->setUsesPeriodicBoundaryConditions(periodic);
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_RBTorsionForce_usesPeriodicBoundaryConditions(const OpenMM_RBTorsionForce* target) {
    bool result = reinterpret_cast<const OpenMM::RBTorsionForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::Integrator */
OPENMM_EXPORT void OpenMM_Integrator_destroy(OpenMM_Integrator* target) {
    delete reinterpret_cast<OpenMM::Integrator*>(target);
}
OPENMM_EXPORT double OpenMM_Integrator_getStepSize(const OpenMM_Integrator* target) {
    double result = reinterpret_cast<const OpenMM::Integrator*>(target)->getStepSize();
    return result;
}
OPENMM_EXPORT void OpenMM_Integrator_setStepSize(OpenMM_Integrator* target, double size) {
    reinterpret_cast<OpenMM::Integrator*>(target)->setStepSize(size);
}
OPENMM_EXPORT double OpenMM_Integrator_getConstraintTolerance(const OpenMM_Integrator* target) {
    double result = reinterpret_cast<const OpenMM::Integrator*>(target)->getConstraintTolerance();
    return result;
}
OPENMM_EXPORT void OpenMM_Integrator_setConstraintTolerance(OpenMM_Integrator* target, double tol) {
    reinterpret_cast<OpenMM::Integrator*>(target)->setConstraintTolerance(tol);
}
OPENMM_EXPORT void OpenMM_Integrator_step(OpenMM_Integrator* target, int steps) {
    reinterpret_cast<OpenMM::Integrator*>(target)->step(steps);
}

/* OpenMM::VerletIntegrator */
OPENMM_EXPORT OpenMM_VerletIntegrator* OpenMM_VerletIntegrator_create(double stepSize) {
    return reinterpret_cast<OpenMM_VerletIntegrator*>(new OpenMM::VerletIntegrator(stepSize));
}
OPENMM_EXPORT void OpenMM_VerletIntegrator_destroy(OpenMM_VerletIntegrator* target) {
    delete reinterpret_cast<OpenMM::VerletIntegrator*>(target);
}
OPENMM_EXPORT void OpenMM_VerletIntegrator_step(OpenMM_VerletIntegrator* target, int steps) {
    reinterpret_cast<OpenMM::VerletIntegrator*>(target)->step(steps);
}

/* OpenMM::LangevinIntegrator */
OPENMM_EXPORT OpenMM_LangevinIntegrator* OpenMM_LangevinIntegrator_create(double temperature, double frictionCoeff, double stepSize) {
    return reinterpret_cast<OpenMM_LangevinIntegrator*>(new OpenMM::LangevinIntegrator(temperature, frictionCoeff, stepSize));
}
OPENMM_EXPORT void OpenMM_LangevinIntegrator_destroy(OpenMM_LangevinIntegrator* target) {
    delete reinterpret_cast<OpenMM::LangevinIntegrator*>(target);
}
OPENMM_EXPORT double OpenMM_LangevinIntegrator_getTemperature(const OpenMM_LangevinIntegrator* target) {
    double result = reinterpret_cast<const OpenMM::LangevinIntegrator*>(target)->getTemperature();
    return result;
}
OPENMM_EXPORT void OpenMM_LangevinIntegrator_setTemperature(OpenMM_LangevinIntegrator* target, double temp) {
    reinterpret_cast<OpenMM::LangevinIntegrator*>(target)->setTemperature(temp);
}
OPENMM_EXPORT double OpenMM_LangevinIntegrator_getFriction(const OpenMM_LangevinIntegrator* target) {
    double result = reinterpret_cast<const OpenMM::LangevinIntegrator*>(target)->getFriction();
    return result;
}
OPENMM_EXPORT void OpenMM_LangevinIntegrator_setFriction(OpenMM_LangevinIntegrator* target, double coeff) {
    reinterpret_cast<OpenMM::LangevinIntegrator*>(target)->setFriction(coeff);
}
OPENMM_EXPORT int OpenMM_LangevinIntegrator_getRandomNumberSeed(const OpenMM_LangevinIntegrator* target) {
    int result = reinterpret_cast<const OpenMM::LangevinIntegrator*>(target)->getRandomNumberSeed();
    return result;
}
OPENMM_EXPORT void OpenMM_LangevinIntegrator_setRandomNumberSeed(OpenMM_LangevinIntegrator* target, int seed) {
    reinterpret_cast<OpenMM::LangevinIntegrator*>(target)->setRandomNumberSeed(seed);
}
OPENMM_EXPORT void OpenMM_LangevinIntegrator_step(OpenMM_LangevinIntegrator* target, int steps) {
    reinterpret_cast<OpenMM::LangevinIntegrator*>(target)->step(steps);
}

/* OpenMM::BAOABLangevinIntegrator */
OPENMM_EXPORT OpenMM_BAOABLangevinIntegrator* OpenMM_BAOABLangevinIntegrator_create(double temperature, double frictionCoeff, double stepSize) {
    return reinterpret_cast<OpenMM_BAOABLangevinIntegrator*>(new OpenMM::BAOABLangevinIntegrator(temperature, frictionCoeff, stepSize));
}
OPENMM_EXPORT void OpenMM_BAOABLangevinIntegrator_destroy(OpenMM_BAOABLangevinIntegrator* target) {
    delete reinterpret_cast<OpenMM::BAOABLangevinIntegrator*>(target);
}
OPENMM_EXPORT double OpenMM_BAOABLangevinIntegrator_getTemperature(const OpenMM_BAOABLangevinIntegrator* target) {
    double result = reinterpret_cast<const OpenMM::BAOABLangevinIntegrator*>(target)->getTemperature();
    return result;
}
OPENMM_EXPORT void OpenMM_BAOABLangevinIntegrator_setTemperature(OpenMM_BAOABLangevinIntegrator* target, double temp) {
    reinterpret_cast<OpenMM::BAOABLangevinIntegrator*>(target)->setTemperature(temp);
}
OPENMM_EXPORT double OpenMM_BAOABLangevinIntegrator_getFriction(const OpenMM_BAOABLangevinIntegrator* target) {
    double result = reinterpret_cast<const OpenMM::BAOABLangevinIntegrator*>(target)->getFriction();
    return result;
}
OPENMM_EXPORT void OpenMM_BAOABLangevinIntegrator_setFriction(OpenMM_BAOABLangevinIntegrator* target, double coeff) {
    reinterpret_cast<OpenMM::BAOABLangevinIntegrator*>(target)->setFriction(coeff);
}
OPENMM_EXPORT int OpenMM_BAOABLangevinIntegrator_getRandomNumberSeed(const OpenMM_BAOABLangevinIntegrator* target) {
    int result = reinterpret_cast<const OpenMM::BAOABLangevinIntegrator*>(target)->getRandomNumberSeed();
    return result;
}
OPENMM_EXPORT void OpenMM_BAOABLangevinIntegrator_setRandomNumberSeed(OpenMM_BAOABLangevinIntegrator* target, int seed) {
    reinterpret_cast<OpenMM::BAOABLangevinIntegrator*>(target)->setRandomNumberSeed(seed);
}
OPENMM_EXPORT void OpenMM_BAOABLangevinIntegrator_step(OpenMM_BAOABLangevinIntegrator* target, int steps) {
    reinterpret_cast<OpenMM::BAOABLangevinIntegrator*>(target)->step(steps);
}

/* OpenMM::TabulatedFunction */
OPENMM_EXPORT void OpenMM_TabulatedFunction_destroy(OpenMM_TabulatedFunction* target) {
    delete reinterpret_cast<OpenMM::TabulatedFunction*>(target);
}
OPENMM_EXPORT OpenMM_TabulatedFunction* OpenMM_TabulatedFunction_Copy(const OpenMM_TabulatedFunction* target) {
    TabulatedFunction * result = reinterpret_cast<const OpenMM::TabulatedFunction*>(target)->Copy();
    return reinterpret_cast<OpenMM_TabulatedFunction*>(result);
}

/* OpenMM::Context */
OPENMM_EXPORT OpenMM_Context* OpenMM_Context_create(const OpenMM_System* system, OpenMM_Integrator* integrator) {
    return reinterpret_cast<OpenMM_Context*>(new OpenMM::Context(*reinterpret_cast<const System*>(system), *reinterpret_cast<OpenMM::Integrator*>(integrator)));
}
OPENMM_EXPORT OpenMM_Context* OpenMM_Context_create_2(const OpenMM_System* system, OpenMM_Integrator* integrator, OpenMM_Platform* platform) {
    return reinterpret_cast<OpenMM_Context*>(new OpenMM::Context(*reinterpret_cast<const System*>(system), *reinterpret_cast<OpenMM::Integrator*>(integrator), *reinterpret_cast<OpenMM::Platform*>(platform)));
}
OPENMM_EXPORT OpenMM_Context* OpenMM_Context_create_3(const OpenMM_System* system, OpenMM_Integrator* integrator, OpenMM_Platform* platform, const OpenMM_PropertyArray* properties) {
    return reinterpret_cast<OpenMM_Context*>(new OpenMM::Context(*reinterpret_cast<const System*>(system), *reinterpret_cast<OpenMM::Integrator*>(integrator), *reinterpret_cast<OpenMM::Platform*>(platform), *reinterpret_cast<const std::map< std::string, std::string >*>(properties)));
}
OPENMM_EXPORT void OpenMM_Context_destroy(OpenMM_Context* target) {
    delete reinterpret_cast<OpenMM::Context*>(target);
}
OPENMM_EXPORT const OpenMM_System* OpenMM_Context_getSystem(const OpenMM_Context* target) {
    const System* result = &reinterpret_cast<const OpenMM::Context*>(target)->getSystem();
    return reinterpret_cast<const OpenMM_System*>(result);
}
OPENMM_EXPORT OpenMM_Integrator* OpenMM_Context_getIntegrator(OpenMM_Context* target) {
    Integrator* result = &reinterpret_cast<OpenMM::Context*>(target)->getIntegrator();
    return reinterpret_cast<OpenMM_Integrator*>(result);
}
OPENMM_EXPORT OpenMM_Platform* OpenMM_Context_getPlatform(OpenMM_Context* target) {
    Platform* result = &reinterpret_cast<OpenMM::Context*>(target)->getPlatform();
    return reinterpret_cast<OpenMM_Platform*>(result);
}
OPENMM_EXPORT void OpenMM_Context_setState(OpenMM_Context* target, const OpenMM_State* state) {
    reinterpret_cast<OpenMM::Context*>(target)->setState(*reinterpret_cast<const State*>(state));
}
OPENMM_EXPORT void OpenMM_Context_setTime(OpenMM_Context* target, double time) {
    reinterpret_cast<OpenMM::Context*>(target)->setTime(time);
}
OPENMM_EXPORT void OpenMM_Context_setPositions(OpenMM_Context* target, const OpenMM_Vec3Array* positions) {
    reinterpret_cast<OpenMM::Context*>(target)->setPositions(*reinterpret_cast<const std::vector< Vec3 >*>(positions));
}
OPENMM_EXPORT void OpenMM_Context_setVelocities(OpenMM_Context* target, const OpenMM_Vec3Array* velocities) {
    reinterpret_cast<OpenMM::Context*>(target)->setVelocities(*reinterpret_cast<const std::vector< Vec3 >*>(velocities));
}
OPENMM_EXPORT void OpenMM_Context_setVelocitiesToTemperature(OpenMM_Context* target, double temperature, int randomSeed) {
    reinterpret_cast<OpenMM::Context*>(target)->setVelocitiesToTemperature(temperature, randomSeed);
}
OPENMM_EXPORT const OpenMM_ParameterArray* OpenMM_Context_getParameters(const OpenMM_Context* target) {
    const std::map< std::string, double >* result = &reinterpret_cast<const OpenMM::Context*>(target)->getParameters();
    return reinterpret_cast<const OpenMM_ParameterArray*>(result);
}
OPENMM_EXPORT double OpenMM_Context_getParameter(const OpenMM_Context* target, const char* name) {
    double result = reinterpret_cast<const OpenMM::Context*>(target)->getParameter(std::string(name));
    return result;
}
OPENMM_EXPORT void OpenMM_Context_setParameter(OpenMM_Context* target, const char* name, double value) {
    reinterpret_cast<OpenMM::Context*>(target)->setParameter(std::string(name), value);
}
OPENMM_EXPORT void OpenMM_Context_setPeriodicBoxVectors(OpenMM_Context* target, const OpenMM_Vec3* a, const OpenMM_Vec3* b, const OpenMM_Vec3* c) {
    reinterpret_cast<OpenMM::Context*>(target)->setPeriodicBoxVectors(*reinterpret_cast<const Vec3*>(a), *reinterpret_cast<const Vec3*>(b), *reinterpret_cast<const Vec3*>(c));
}
OPENMM_EXPORT void OpenMM_Context_applyConstraints(OpenMM_Context* target, double tol) {
    reinterpret_cast<OpenMM::Context*>(target)->applyConstraints(tol);
}
OPENMM_EXPORT void OpenMM_Context_applyVelocityConstraints(OpenMM_Context* target, double tol) {
    reinterpret_cast<OpenMM::Context*>(target)->applyVelocityConstraints(tol);
}
OPENMM_EXPORT void OpenMM_Context_computeVirtualSites(OpenMM_Context* target) {
    reinterpret_cast<OpenMM::Context*>(target)->computeVirtualSites();
}
OPENMM_EXPORT void OpenMM_Context_reinitialize(OpenMM_Context* target, OpenMM_Boolean preserveState) {
    reinterpret_cast<OpenMM::Context*>(target)->reinitialize(preserveState);
}

/* OpenMM::CMMotionRemover */
OPENMM_EXPORT OpenMM_CMMotionRemover* OpenMM_CMMotionRemover_create(int frequency) {
    return reinterpret_cast<OpenMM_CMMotionRemover*>(new OpenMM::CMMotionRemover(frequency));
}
OPENMM_EXPORT void OpenMM_CMMotionRemover_destroy(OpenMM_CMMotionRemover* target) {
    delete reinterpret_cast<OpenMM::CMMotionRemover*>(target);
}
OPENMM_EXPORT int OpenMM_CMMotionRemover_getFrequency(const OpenMM_CMMotionRemover* target) {
    int result = reinterpret_cast<const OpenMM::CMMotionRemover*>(target)->getFrequency();
    return result;
}
OPENMM_EXPORT void OpenMM_CMMotionRemover_setFrequency(OpenMM_CMMotionRemover* target, int freq) {
    reinterpret_cast<OpenMM::CMMotionRemover*>(target)->setFrequency(freq);
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_CMMotionRemover_usesPeriodicBoundaryConditions(const OpenMM_CMMotionRemover* target) {
    bool result = reinterpret_cast<const OpenMM::CMMotionRemover*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::CustomBondForce */
OPENMM_EXPORT OpenMM_CustomBondForce* OpenMM_CustomBondForce_create(const char* energy) {
    return reinterpret_cast<OpenMM_CustomBondForce*>(new OpenMM::CustomBondForce(std::string(energy)));
}
OPENMM_EXPORT void OpenMM_CustomBondForce_destroy(OpenMM_CustomBondForce* target) {
    delete reinterpret_cast<OpenMM::CustomBondForce*>(target);
}
OPENMM_EXPORT int OpenMM_CustomBondForce_getNumBonds(const OpenMM_CustomBondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomBondForce*>(target)->getNumBonds();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomBondForce_getNumPerBondParameters(const OpenMM_CustomBondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomBondForce*>(target)->getNumPerBondParameters();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomBondForce_getNumGlobalParameters(const OpenMM_CustomBondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomBondForce*>(target)->getNumGlobalParameters();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomBondForce_getNumEnergyParameterDerivatives(const OpenMM_CustomBondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomBondForce*>(target)->getNumEnergyParameterDerivatives();
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomBondForce_getEnergyFunction(const OpenMM_CustomBondForce* target) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomBondForce*>(target)->getEnergyFunction();
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomBondForce_setEnergyFunction(OpenMM_CustomBondForce* target, const char* energy) {
    reinterpret_cast<OpenMM::CustomBondForce*>(target)->setEnergyFunction(std::string(energy));
}
OPENMM_EXPORT int OpenMM_CustomBondForce_addPerBondParameter(OpenMM_CustomBondForce* target, const char* name) {
    int result = reinterpret_cast<OpenMM::CustomBondForce*>(target)->addPerBondParameter(std::string(name));
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomBondForce_getPerBondParameterName(const OpenMM_CustomBondForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomBondForce*>(target)->getPerBondParameterName(index);
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomBondForce_setPerBondParameterName(OpenMM_CustomBondForce* target, int index, const char* name) {
    reinterpret_cast<OpenMM::CustomBondForce*>(target)->setPerBondParameterName(index, std::string(name));
}
OPENMM_EXPORT int OpenMM_CustomBondForce_addGlobalParameter(OpenMM_CustomBondForce* target, const char* name, double defaultValue) {
    int result = reinterpret_cast<OpenMM::CustomBondForce*>(target)->addGlobalParameter(std::string(name), defaultValue);
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomBondForce_getGlobalParameterName(const OpenMM_CustomBondForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomBondForce*>(target)->getGlobalParameterName(index);
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomBondForce_setGlobalParameterName(OpenMM_CustomBondForce* target, int index, const char* name) {
    reinterpret_cast<OpenMM::CustomBondForce*>(target)->setGlobalParameterName(index, std::string(name));
}
OPENMM_EXPORT double OpenMM_CustomBondForce_getGlobalParameterDefaultValue(const OpenMM_CustomBondForce* target, int index) {
    double result = reinterpret_cast<const OpenMM::CustomBondForce*>(target)->getGlobalParameterDefaultValue(index);
    return result;
}
OPENMM_EXPORT void OpenMM_CustomBondForce_setGlobalParameterDefaultValue(OpenMM_CustomBondForce* target, int index, double defaultValue) {
    reinterpret_cast<OpenMM::CustomBondForce*>(target)->setGlobalParameterDefaultValue(index, defaultValue);
}
OPENMM_EXPORT void OpenMM_CustomBondForce_addEnergyParameterDerivative(OpenMM_CustomBondForce* target, const char* name) {
    reinterpret_cast<OpenMM::CustomBondForce*>(target)->addEnergyParameterDerivative(std::string(name));
}
OPENMM_EXPORT const char* OpenMM_CustomBondForce_getEnergyParameterDerivativeName(const OpenMM_CustomBondForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomBondForce*>(target)->getEnergyParameterDerivativeName(index);
    return result->c_str();
}
OPENMM_EXPORT int OpenMM_CustomBondForce_addBond(OpenMM_CustomBondForce* target, int particle1, int particle2, const OpenMM_DoubleArray* parameters) {
    int result = reinterpret_cast<OpenMM::CustomBondForce*>(target)->addBond(particle1, particle2, *reinterpret_cast<const std::vector< double >*>(parameters));
    return result;
}
OPENMM_EXPORT void OpenMM_CustomBondForce_getBondParameters(const OpenMM_CustomBondForce* target, int index, int* particle1, int* particle2, OpenMM_DoubleArray* parameters) {
    reinterpret_cast<const OpenMM::CustomBondForce*>(target)->getBondParameters(index, *reinterpret_cast<int*>(particle1), *reinterpret_cast<int*>(particle2), *reinterpret_cast<std::vector< double >*>(parameters));
}
OPENMM_EXPORT void OpenMM_CustomBondForce_setBondParameters(OpenMM_CustomBondForce* target, int index, int particle1, int particle2, const OpenMM_DoubleArray* parameters) {
    reinterpret_cast<OpenMM::CustomBondForce*>(target)->setBondParameters(index, particle1, particle2, *reinterpret_cast<const std::vector< double >*>(parameters));
}
OPENMM_EXPORT void OpenMM_CustomBondForce_updateParametersInContext(OpenMM_CustomBondForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::CustomBondForce*>(target)->updateParametersInContext(*reinterpret_cast<OpenMM::Context*>(context));
}
OPENMM_EXPORT void OpenMM_CustomBondForce_setUsesPeriodicBoundaryConditions(OpenMM_CustomBondForce* target, OpenMM_Boolean periodic) {
    reinterpret_cast<OpenMM::CustomBondForce*>(target)->setUsesPeriodicBoundaryConditions(periodic);
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_CustomBondForce_usesPeriodicBoundaryConditions(const OpenMM_CustomBondForce* target) {
    bool result = reinterpret_cast<const OpenMM::CustomBondForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::MonteCarloMembraneBarostat */
OPENMM_EXPORT OpenMM_MonteCarloMembraneBarostat* OpenMM_MonteCarloMembraneBarostat_create(double defaultPressure, double defaultSurfaceTension, double defaultTemperature, OpenMM_MonteCarloMembraneBarostat_XYMode xymode, OpenMM_MonteCarloMembraneBarostat_ZMode zmode, int frequency) {
    return reinterpret_cast<OpenMM_MonteCarloMembraneBarostat*>(new OpenMM::MonteCarloMembraneBarostat(defaultPressure, defaultSurfaceTension, defaultTemperature, static_cast<OpenMM::MonteCarloMembraneBarostat::XYMode>(xymode), static_cast<OpenMM::MonteCarloMembraneBarostat::ZMode>(zmode), frequency));
}
OPENMM_EXPORT void OpenMM_MonteCarloMembraneBarostat_destroy(OpenMM_MonteCarloMembraneBarostat* target) {
    delete reinterpret_cast<OpenMM::MonteCarloMembraneBarostat*>(target);
}
OPENMM_EXPORT const char* OpenMM_MonteCarloMembraneBarostat_Pressure() {
    const std::string* result = &OpenMM::MonteCarloMembraneBarostat::Pressure();
    return result->c_str();
}
OPENMM_EXPORT const char* OpenMM_MonteCarloMembraneBarostat_SurfaceTension() {
    const std::string* result = &OpenMM::MonteCarloMembraneBarostat::SurfaceTension();
    return result->c_str();
}
OPENMM_EXPORT const char* OpenMM_MonteCarloMembraneBarostat_Temperature() {
    const std::string* result = &OpenMM::MonteCarloMembraneBarostat::Temperature();
    return result->c_str();
}
OPENMM_EXPORT double OpenMM_MonteCarloMembraneBarostat_getDefaultPressure(const OpenMM_MonteCarloMembraneBarostat* target) {
    double result = reinterpret_cast<const OpenMM::MonteCarloMembraneBarostat*>(target)->getDefaultPressure();
    return result;
}
OPENMM_EXPORT void OpenMM_MonteCarloMembraneBarostat_setDefaultPressure(OpenMM_MonteCarloMembraneBarostat* target, double pressure) {
    reinterpret_cast<OpenMM::MonteCarloMembraneBarostat*>(target)->setDefaultPressure(pressure);
}
OPENMM_EXPORT double OpenMM_MonteCarloMembraneBarostat_getDefaultSurfaceTension(const OpenMM_MonteCarloMembraneBarostat* target) {
    double result = reinterpret_cast<const OpenMM::MonteCarloMembraneBarostat*>(target)->getDefaultSurfaceTension();
    return result;
}
OPENMM_EXPORT void OpenMM_MonteCarloMembraneBarostat_setDefaultSurfaceTension(OpenMM_MonteCarloMembraneBarostat* target, double surfaceTension) {
    reinterpret_cast<OpenMM::MonteCarloMembraneBarostat*>(target)->setDefaultSurfaceTension(surfaceTension);
}
OPENMM_EXPORT int OpenMM_MonteCarloMembraneBarostat_getFrequency(const OpenMM_MonteCarloMembraneBarostat* target) {
    int result = reinterpret_cast<const OpenMM::MonteCarloMembraneBarostat*>(target)->getFrequency();
    return result;
}
OPENMM_EXPORT void OpenMM_MonteCarloMembraneBarostat_setFrequency(OpenMM_MonteCarloMembraneBarostat* target, int freq) {
    reinterpret_cast<OpenMM::MonteCarloMembraneBarostat*>(target)->setFrequency(freq);
}
OPENMM_EXPORT double OpenMM_MonteCarloMembraneBarostat_getDefaultTemperature(const OpenMM_MonteCarloMembraneBarostat* target) {
    double result = reinterpret_cast<const OpenMM::MonteCarloMembraneBarostat*>(target)->getDefaultTemperature();
    return result;
}
OPENMM_EXPORT void OpenMM_MonteCarloMembraneBarostat_setDefaultTemperature(OpenMM_MonteCarloMembraneBarostat* target, double temp) {
    reinterpret_cast<OpenMM::MonteCarloMembraneBarostat*>(target)->setDefaultTemperature(temp);
}
OPENMM_EXPORT OpenMM_MonteCarloMembraneBarostat_XYMode OpenMM_MonteCarloMembraneBarostat_getXYMode(const OpenMM_MonteCarloMembraneBarostat* target) {
    OpenMM::MonteCarloMembraneBarostat::XYMode result = reinterpret_cast<const OpenMM::MonteCarloMembraneBarostat*>(target)->getXYMode();
    return static_cast<OpenMM_MonteCarloMembraneBarostat_XYMode>(result);
}
OPENMM_EXPORT void OpenMM_MonteCarloMembraneBarostat_setXYMode(OpenMM_MonteCarloMembraneBarostat* target, OpenMM_MonteCarloMembraneBarostat_XYMode mode) {
    reinterpret_cast<OpenMM::MonteCarloMembraneBarostat*>(target)->setXYMode(static_cast<OpenMM::MonteCarloMembraneBarostat::XYMode>(mode));
}
OPENMM_EXPORT OpenMM_MonteCarloMembraneBarostat_ZMode OpenMM_MonteCarloMembraneBarostat_getZMode(const OpenMM_MonteCarloMembraneBarostat* target) {
    OpenMM::MonteCarloMembraneBarostat::ZMode result = reinterpret_cast<const OpenMM::MonteCarloMembraneBarostat*>(target)->getZMode();
    return static_cast<OpenMM_MonteCarloMembraneBarostat_ZMode>(result);
}
OPENMM_EXPORT void OpenMM_MonteCarloMembraneBarostat_setZMode(OpenMM_MonteCarloMembraneBarostat* target, OpenMM_MonteCarloMembraneBarostat_ZMode mode) {
    reinterpret_cast<OpenMM::MonteCarloMembraneBarostat*>(target)->setZMode(static_cast<OpenMM::MonteCarloMembraneBarostat::ZMode>(mode));
}
OPENMM_EXPORT int OpenMM_MonteCarloMembraneBarostat_getRandomNumberSeed(const OpenMM_MonteCarloMembraneBarostat* target) {
    int result = reinterpret_cast<const OpenMM::MonteCarloMembraneBarostat*>(target)->getRandomNumberSeed();
    return result;
}
OPENMM_EXPORT void OpenMM_MonteCarloMembraneBarostat_setRandomNumberSeed(OpenMM_MonteCarloMembraneBarostat* target, int seed) {
    reinterpret_cast<OpenMM::MonteCarloMembraneBarostat*>(target)->setRandomNumberSeed(seed);
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_MonteCarloMembraneBarostat_usesPeriodicBoundaryConditions(const OpenMM_MonteCarloMembraneBarostat* target) {
    bool result = reinterpret_cast<const OpenMM::MonteCarloMembraneBarostat*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::VirtualSite */
OPENMM_EXPORT void OpenMM_VirtualSite_destroy(OpenMM_VirtualSite* target) {
    delete reinterpret_cast<OpenMM::VirtualSite*>(target);
}
OPENMM_EXPORT int OpenMM_VirtualSite_getNumParticles(const OpenMM_VirtualSite* target) {
    int result = reinterpret_cast<const OpenMM::VirtualSite*>(target)->getNumParticles();
    return result;
}
OPENMM_EXPORT int OpenMM_VirtualSite_getParticle(const OpenMM_VirtualSite* target, int particle) {
    int result = reinterpret_cast<const OpenMM::VirtualSite*>(target)->getParticle(particle);
    return result;
}

/* OpenMM::OutOfPlaneSite */
OPENMM_EXPORT OpenMM_OutOfPlaneSite* OpenMM_OutOfPlaneSite_create(int particle1, int particle2, int particle3, double weight12, double weight13, double weightCross) {
    return reinterpret_cast<OpenMM_OutOfPlaneSite*>(new OpenMM::OutOfPlaneSite(particle1, particle2, particle3, weight12, weight13, weightCross));
}
OPENMM_EXPORT void OpenMM_OutOfPlaneSite_destroy(OpenMM_OutOfPlaneSite* target) {
    delete reinterpret_cast<OpenMM::OutOfPlaneSite*>(target);
}
OPENMM_EXPORT double OpenMM_OutOfPlaneSite_getWeight12(const OpenMM_OutOfPlaneSite* target) {
    double result = reinterpret_cast<const OpenMM::OutOfPlaneSite*>(target)->getWeight12();
    return result;
}
OPENMM_EXPORT double OpenMM_OutOfPlaneSite_getWeight13(const OpenMM_OutOfPlaneSite* target) {
    double result = reinterpret_cast<const OpenMM::OutOfPlaneSite*>(target)->getWeight13();
    return result;
}
OPENMM_EXPORT double OpenMM_OutOfPlaneSite_getWeightCross(const OpenMM_OutOfPlaneSite* target) {
    double result = reinterpret_cast<const OpenMM::OutOfPlaneSite*>(target)->getWeightCross();
    return result;
}

/* OpenMM::CustomTorsionForce */
OPENMM_EXPORT OpenMM_CustomTorsionForce* OpenMM_CustomTorsionForce_create(const char* energy) {
    return reinterpret_cast<OpenMM_CustomTorsionForce*>(new OpenMM::CustomTorsionForce(std::string(energy)));
}
OPENMM_EXPORT void OpenMM_CustomTorsionForce_destroy(OpenMM_CustomTorsionForce* target) {
    delete reinterpret_cast<OpenMM::CustomTorsionForce*>(target);
}
OPENMM_EXPORT int OpenMM_CustomTorsionForce_getNumTorsions(const OpenMM_CustomTorsionForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomTorsionForce*>(target)->getNumTorsions();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomTorsionForce_getNumPerTorsionParameters(const OpenMM_CustomTorsionForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomTorsionForce*>(target)->getNumPerTorsionParameters();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomTorsionForce_getNumGlobalParameters(const OpenMM_CustomTorsionForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomTorsionForce*>(target)->getNumGlobalParameters();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomTorsionForce_getNumEnergyParameterDerivatives(const OpenMM_CustomTorsionForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomTorsionForce*>(target)->getNumEnergyParameterDerivatives();
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomTorsionForce_getEnergyFunction(const OpenMM_CustomTorsionForce* target) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomTorsionForce*>(target)->getEnergyFunction();
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomTorsionForce_setEnergyFunction(OpenMM_CustomTorsionForce* target, const char* energy) {
    reinterpret_cast<OpenMM::CustomTorsionForce*>(target)->setEnergyFunction(std::string(energy));
}
OPENMM_EXPORT int OpenMM_CustomTorsionForce_addPerTorsionParameter(OpenMM_CustomTorsionForce* target, const char* name) {
    int result = reinterpret_cast<OpenMM::CustomTorsionForce*>(target)->addPerTorsionParameter(std::string(name));
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomTorsionForce_getPerTorsionParameterName(const OpenMM_CustomTorsionForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomTorsionForce*>(target)->getPerTorsionParameterName(index);
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomTorsionForce_setPerTorsionParameterName(OpenMM_CustomTorsionForce* target, int index, const char* name) {
    reinterpret_cast<OpenMM::CustomTorsionForce*>(target)->setPerTorsionParameterName(index, std::string(name));
}
OPENMM_EXPORT int OpenMM_CustomTorsionForce_addGlobalParameter(OpenMM_CustomTorsionForce* target, const char* name, double defaultValue) {
    int result = reinterpret_cast<OpenMM::CustomTorsionForce*>(target)->addGlobalParameter(std::string(name), defaultValue);
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomTorsionForce_getGlobalParameterName(const OpenMM_CustomTorsionForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomTorsionForce*>(target)->getGlobalParameterName(index);
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomTorsionForce_setGlobalParameterName(OpenMM_CustomTorsionForce* target, int index, const char* name) {
    reinterpret_cast<OpenMM::CustomTorsionForce*>(target)->setGlobalParameterName(index, std::string(name));
}
OPENMM_EXPORT double OpenMM_CustomTorsionForce_getGlobalParameterDefaultValue(const OpenMM_CustomTorsionForce* target, int index) {
    double result = reinterpret_cast<const OpenMM::CustomTorsionForce*>(target)->getGlobalParameterDefaultValue(index);
    return result;
}
OPENMM_EXPORT void OpenMM_CustomTorsionForce_setGlobalParameterDefaultValue(OpenMM_CustomTorsionForce* target, int index, double defaultValue) {
    reinterpret_cast<OpenMM::CustomTorsionForce*>(target)->setGlobalParameterDefaultValue(index, defaultValue);
}
OPENMM_EXPORT void OpenMM_CustomTorsionForce_addEnergyParameterDerivative(OpenMM_CustomTorsionForce* target, const char* name) {
    reinterpret_cast<OpenMM::CustomTorsionForce*>(target)->addEnergyParameterDerivative(std::string(name));
}
OPENMM_EXPORT const char* OpenMM_CustomTorsionForce_getEnergyParameterDerivativeName(const OpenMM_CustomTorsionForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomTorsionForce*>(target)->getEnergyParameterDerivativeName(index);
    return result->c_str();
}
OPENMM_EXPORT int OpenMM_CustomTorsionForce_addTorsion(OpenMM_CustomTorsionForce* target, int particle1, int particle2, int particle3, int particle4, const OpenMM_DoubleArray* parameters) {
    int result = reinterpret_cast<OpenMM::CustomTorsionForce*>(target)->addTorsion(particle1, particle2, particle3, particle4, *reinterpret_cast<const std::vector< double >*>(parameters));
    return result;
}
OPENMM_EXPORT void OpenMM_CustomTorsionForce_getTorsionParameters(const OpenMM_CustomTorsionForce* target, int index, int* particle1, int* particle2, int* particle3, int* particle4, OpenMM_DoubleArray* parameters) {
    reinterpret_cast<const OpenMM::CustomTorsionForce*>(target)->getTorsionParameters(index, *reinterpret_cast<int*>(particle1), *reinterpret_cast<int*>(particle2), *reinterpret_cast<int*>(particle3), *reinterpret_cast<int*>(particle4), *reinterpret_cast<std::vector< double >*>(parameters));
}
OPENMM_EXPORT void OpenMM_CustomTorsionForce_setTorsionParameters(OpenMM_CustomTorsionForce* target, int index, int particle1, int particle2, int particle3, int particle4, const OpenMM_DoubleArray* parameters) {
    reinterpret_cast<OpenMM::CustomTorsionForce*>(target)->setTorsionParameters(index, particle1, particle2, particle3, particle4, *reinterpret_cast<const std::vector< double >*>(parameters));
}
OPENMM_EXPORT void OpenMM_CustomTorsionForce_updateParametersInContext(OpenMM_CustomTorsionForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::CustomTorsionForce*>(target)->updateParametersInContext(*reinterpret_cast<OpenMM::Context*>(context));
}
OPENMM_EXPORT void OpenMM_CustomTorsionForce_setUsesPeriodicBoundaryConditions(OpenMM_CustomTorsionForce* target, OpenMM_Boolean periodic) {
    reinterpret_cast<OpenMM::CustomTorsionForce*>(target)->setUsesPeriodicBoundaryConditions(periodic);
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_CustomTorsionForce_usesPeriodicBoundaryConditions(const OpenMM_CustomTorsionForce* target) {
    bool result = reinterpret_cast<const OpenMM::CustomTorsionForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::Discrete2DFunction */
OPENMM_EXPORT OpenMM_Discrete2DFunction* OpenMM_Discrete2DFunction_create(int xsize, int ysize, const OpenMM_DoubleArray* values) {
    return reinterpret_cast<OpenMM_Discrete2DFunction*>(new OpenMM::Discrete2DFunction(xsize, ysize, *reinterpret_cast<const std::vector< double >*>(values)));
}
OPENMM_EXPORT void OpenMM_Discrete2DFunction_destroy(OpenMM_Discrete2DFunction* target) {
    delete reinterpret_cast<OpenMM::Discrete2DFunction*>(target);
}
OPENMM_EXPORT void OpenMM_Discrete2DFunction_getFunctionParameters(const OpenMM_Discrete2DFunction* target, int* xsize, int* ysize, OpenMM_DoubleArray* values) {
    reinterpret_cast<const OpenMM::Discrete2DFunction*>(target)->getFunctionParameters(*reinterpret_cast<int*>(xsize), *reinterpret_cast<int*>(ysize), *reinterpret_cast<std::vector< double >*>(values));
}
OPENMM_EXPORT void OpenMM_Discrete2DFunction_setFunctionParameters(OpenMM_Discrete2DFunction* target, int xsize, int ysize, const OpenMM_DoubleArray* values) {
    reinterpret_cast<OpenMM::Discrete2DFunction*>(target)->setFunctionParameters(xsize, ysize, *reinterpret_cast<const std::vector< double >*>(values));
}
OPENMM_EXPORT OpenMM_Discrete2DFunction* OpenMM_Discrete2DFunction_Copy(const OpenMM_Discrete2DFunction* target) {
    Discrete2DFunction * result = reinterpret_cast<const OpenMM::Discrete2DFunction*>(target)->Copy();
    return reinterpret_cast<OpenMM_Discrete2DFunction*>(result);
}

/* OpenMM::AndersenThermostat */
OPENMM_EXPORT OpenMM_AndersenThermostat* OpenMM_AndersenThermostat_create(double defaultTemperature, double defaultCollisionFrequency) {
    return reinterpret_cast<OpenMM_AndersenThermostat*>(new OpenMM::AndersenThermostat(defaultTemperature, defaultCollisionFrequency));
}
OPENMM_EXPORT void OpenMM_AndersenThermostat_destroy(OpenMM_AndersenThermostat* target) {
    delete reinterpret_cast<OpenMM::AndersenThermostat*>(target);
}
OPENMM_EXPORT const char* OpenMM_AndersenThermostat_Temperature() {
    const std::string* result = &OpenMM::AndersenThermostat::Temperature();
    return result->c_str();
}
OPENMM_EXPORT const char* OpenMM_AndersenThermostat_CollisionFrequency() {
    const std::string* result = &OpenMM::AndersenThermostat::CollisionFrequency();
    return result->c_str();
}
OPENMM_EXPORT double OpenMM_AndersenThermostat_getDefaultTemperature(const OpenMM_AndersenThermostat* target) {
    double result = reinterpret_cast<const OpenMM::AndersenThermostat*>(target)->getDefaultTemperature();
    return result;
}
OPENMM_EXPORT void OpenMM_AndersenThermostat_setDefaultTemperature(OpenMM_AndersenThermostat* target, double temperature) {
    reinterpret_cast<OpenMM::AndersenThermostat*>(target)->setDefaultTemperature(temperature);
}
OPENMM_EXPORT double OpenMM_AndersenThermostat_getDefaultCollisionFrequency(const OpenMM_AndersenThermostat* target) {
    double result = reinterpret_cast<const OpenMM::AndersenThermostat*>(target)->getDefaultCollisionFrequency();
    return result;
}
OPENMM_EXPORT void OpenMM_AndersenThermostat_setDefaultCollisionFrequency(OpenMM_AndersenThermostat* target, double frequency) {
    reinterpret_cast<OpenMM::AndersenThermostat*>(target)->setDefaultCollisionFrequency(frequency);
}
OPENMM_EXPORT int OpenMM_AndersenThermostat_getRandomNumberSeed(const OpenMM_AndersenThermostat* target) {
    int result = reinterpret_cast<const OpenMM::AndersenThermostat*>(target)->getRandomNumberSeed();
    return result;
}
OPENMM_EXPORT void OpenMM_AndersenThermostat_setRandomNumberSeed(OpenMM_AndersenThermostat* target, int seed) {
    reinterpret_cast<OpenMM::AndersenThermostat*>(target)->setRandomNumberSeed(seed);
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_AndersenThermostat_usesPeriodicBoundaryConditions(const OpenMM_AndersenThermostat* target) {
    bool result = reinterpret_cast<const OpenMM::AndersenThermostat*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::ThreeParticleAverageSite */
OPENMM_EXPORT OpenMM_ThreeParticleAverageSite* OpenMM_ThreeParticleAverageSite_create(int particle1, int particle2, int particle3, double weight1, double weight2, double weight3) {
    return reinterpret_cast<OpenMM_ThreeParticleAverageSite*>(new OpenMM::ThreeParticleAverageSite(particle1, particle2, particle3, weight1, weight2, weight3));
}
OPENMM_EXPORT void OpenMM_ThreeParticleAverageSite_destroy(OpenMM_ThreeParticleAverageSite* target) {
    delete reinterpret_cast<OpenMM::ThreeParticleAverageSite*>(target);
}
OPENMM_EXPORT double OpenMM_ThreeParticleAverageSite_getWeight(const OpenMM_ThreeParticleAverageSite* target, int particle) {
    double result = reinterpret_cast<const OpenMM::ThreeParticleAverageSite*>(target)->getWeight(particle);
    return result;
}

/* OpenMM::CustomCVForce */
OPENMM_EXPORT OpenMM_CustomCVForce* OpenMM_CustomCVForce_create(const char* energy) {
    return reinterpret_cast<OpenMM_CustomCVForce*>(new OpenMM::CustomCVForce(std::string(energy)));
}
OPENMM_EXPORT void OpenMM_CustomCVForce_destroy(OpenMM_CustomCVForce* target) {
    delete reinterpret_cast<OpenMM::CustomCVForce*>(target);
}
OPENMM_EXPORT int OpenMM_CustomCVForce_getNumCollectiveVariables(const OpenMM_CustomCVForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomCVForce*>(target)->getNumCollectiveVariables();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomCVForce_getNumGlobalParameters(const OpenMM_CustomCVForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomCVForce*>(target)->getNumGlobalParameters();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomCVForce_getNumEnergyParameterDerivatives(const OpenMM_CustomCVForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomCVForce*>(target)->getNumEnergyParameterDerivatives();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomCVForce_getNumTabulatedFunctions(const OpenMM_CustomCVForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomCVForce*>(target)->getNumTabulatedFunctions();
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomCVForce_getEnergyFunction(const OpenMM_CustomCVForce* target) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomCVForce*>(target)->getEnergyFunction();
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomCVForce_setEnergyFunction(OpenMM_CustomCVForce* target, const char* energy) {
    reinterpret_cast<OpenMM::CustomCVForce*>(target)->setEnergyFunction(std::string(energy));
}
OPENMM_EXPORT int OpenMM_CustomCVForce_addCollectiveVariable(OpenMM_CustomCVForce* target, const char* name, OpenMM_Force* variable) {
    int result = reinterpret_cast<OpenMM::CustomCVForce*>(target)->addCollectiveVariable(std::string(name), reinterpret_cast<Force *>(variable));
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomCVForce_getCollectiveVariableName(const OpenMM_CustomCVForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomCVForce*>(target)->getCollectiveVariableName(index);
    return result->c_str();
}
OPENMM_EXPORT OpenMM_Force* OpenMM_CustomCVForce_getCollectiveVariable(OpenMM_CustomCVForce* target, int index) {
    Force* result = &reinterpret_cast<OpenMM::CustomCVForce*>(target)->getCollectiveVariable(index);
    return reinterpret_cast<OpenMM_Force*>(result);
}
OPENMM_EXPORT int OpenMM_CustomCVForce_addGlobalParameter(OpenMM_CustomCVForce* target, const char* name, double defaultValue) {
    int result = reinterpret_cast<OpenMM::CustomCVForce*>(target)->addGlobalParameter(std::string(name), defaultValue);
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomCVForce_getGlobalParameterName(const OpenMM_CustomCVForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomCVForce*>(target)->getGlobalParameterName(index);
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomCVForce_setGlobalParameterName(OpenMM_CustomCVForce* target, int index, const char* name) {
    reinterpret_cast<OpenMM::CustomCVForce*>(target)->setGlobalParameterName(index, std::string(name));
}
OPENMM_EXPORT double OpenMM_CustomCVForce_getGlobalParameterDefaultValue(const OpenMM_CustomCVForce* target, int index) {
    double result = reinterpret_cast<const OpenMM::CustomCVForce*>(target)->getGlobalParameterDefaultValue(index);
    return result;
}
OPENMM_EXPORT void OpenMM_CustomCVForce_setGlobalParameterDefaultValue(OpenMM_CustomCVForce* target, int index, double defaultValue) {
    reinterpret_cast<OpenMM::CustomCVForce*>(target)->setGlobalParameterDefaultValue(index, defaultValue);
}
OPENMM_EXPORT void OpenMM_CustomCVForce_addEnergyParameterDerivative(OpenMM_CustomCVForce* target, const char* name) {
    reinterpret_cast<OpenMM::CustomCVForce*>(target)->addEnergyParameterDerivative(std::string(name));
}
OPENMM_EXPORT const char* OpenMM_CustomCVForce_getEnergyParameterDerivativeName(const OpenMM_CustomCVForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomCVForce*>(target)->getEnergyParameterDerivativeName(index);
    return result->c_str();
}
OPENMM_EXPORT int OpenMM_CustomCVForce_addTabulatedFunction(OpenMM_CustomCVForce* target, const char* name, OpenMM_TabulatedFunction* function) {
    int result = reinterpret_cast<OpenMM::CustomCVForce*>(target)->addTabulatedFunction(std::string(name), reinterpret_cast<TabulatedFunction *>(function));
    return result;
}
OPENMM_EXPORT OpenMM_TabulatedFunction* OpenMM_CustomCVForce_getTabulatedFunction(OpenMM_CustomCVForce* target, int index) {
    TabulatedFunction* result = &reinterpret_cast<OpenMM::CustomCVForce*>(target)->getTabulatedFunction(index);
    return reinterpret_cast<OpenMM_TabulatedFunction*>(result);
}
OPENMM_EXPORT const char* OpenMM_CustomCVForce_getTabulatedFunctionName(const OpenMM_CustomCVForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomCVForce*>(target)->getTabulatedFunctionName(index);
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomCVForce_getCollectiveVariableValues(OpenMM_CustomCVForce* target, OpenMM_Context* context, OpenMM_DoubleArray* values) {
    reinterpret_cast<OpenMM::CustomCVForce*>(target)->getCollectiveVariableValues(*reinterpret_cast<OpenMM::Context*>(context), *reinterpret_cast<std::vector< double >*>(values));
}
OPENMM_EXPORT OpenMM_Context* OpenMM_CustomCVForce_getInnerContext(OpenMM_CustomCVForce* target, OpenMM_Context* context) {
    Context* result = &reinterpret_cast<OpenMM::CustomCVForce*>(target)->getInnerContext(*reinterpret_cast<OpenMM::Context*>(context));
    return reinterpret_cast<OpenMM_Context*>(result);
}
OPENMM_EXPORT void OpenMM_CustomCVForce_updateParametersInContext(OpenMM_CustomCVForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::CustomCVForce*>(target)->updateParametersInContext(*reinterpret_cast<OpenMM::Context*>(context));
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_CustomCVForce_usesPeriodicBoundaryConditions(const OpenMM_CustomCVForce* target) {
    bool result = reinterpret_cast<const OpenMM::CustomCVForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::CustomGBForce */
OPENMM_EXPORT OpenMM_CustomGBForce* OpenMM_CustomGBForce_create() {
    return reinterpret_cast<OpenMM_CustomGBForce*>(new OpenMM::CustomGBForce());
}
OPENMM_EXPORT void OpenMM_CustomGBForce_destroy(OpenMM_CustomGBForce* target) {
    delete reinterpret_cast<OpenMM::CustomGBForce*>(target);
}
OPENMM_EXPORT int OpenMM_CustomGBForce_getNumParticles(const OpenMM_CustomGBForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomGBForce*>(target)->getNumParticles();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomGBForce_getNumExclusions(const OpenMM_CustomGBForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomGBForce*>(target)->getNumExclusions();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomGBForce_getNumPerParticleParameters(const OpenMM_CustomGBForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomGBForce*>(target)->getNumPerParticleParameters();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomGBForce_getNumGlobalParameters(const OpenMM_CustomGBForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomGBForce*>(target)->getNumGlobalParameters();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomGBForce_getNumEnergyParameterDerivatives(const OpenMM_CustomGBForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomGBForce*>(target)->getNumEnergyParameterDerivatives();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomGBForce_getNumTabulatedFunctions(const OpenMM_CustomGBForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomGBForce*>(target)->getNumTabulatedFunctions();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomGBForce_getNumFunctions(const OpenMM_CustomGBForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomGBForce*>(target)->getNumFunctions();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomGBForce_getNumComputedValues(const OpenMM_CustomGBForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomGBForce*>(target)->getNumComputedValues();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomGBForce_getNumEnergyTerms(const OpenMM_CustomGBForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomGBForce*>(target)->getNumEnergyTerms();
    return result;
}
OPENMM_EXPORT OpenMM_CustomGBForce_NonbondedMethod OpenMM_CustomGBForce_getNonbondedMethod(const OpenMM_CustomGBForce* target) {
    OpenMM::CustomGBForce::NonbondedMethod result = reinterpret_cast<const OpenMM::CustomGBForce*>(target)->getNonbondedMethod();
    return static_cast<OpenMM_CustomGBForce_NonbondedMethod>(result);
}
OPENMM_EXPORT void OpenMM_CustomGBForce_setNonbondedMethod(OpenMM_CustomGBForce* target, OpenMM_CustomGBForce_NonbondedMethod method) {
    reinterpret_cast<OpenMM::CustomGBForce*>(target)->setNonbondedMethod(static_cast<OpenMM::CustomGBForce::NonbondedMethod>(method));
}
OPENMM_EXPORT double OpenMM_CustomGBForce_getCutoffDistance(const OpenMM_CustomGBForce* target) {
    double result = reinterpret_cast<const OpenMM::CustomGBForce*>(target)->getCutoffDistance();
    return result;
}
OPENMM_EXPORT void OpenMM_CustomGBForce_setCutoffDistance(OpenMM_CustomGBForce* target, double distance) {
    reinterpret_cast<OpenMM::CustomGBForce*>(target)->setCutoffDistance(distance);
}
OPENMM_EXPORT int OpenMM_CustomGBForce_addPerParticleParameter(OpenMM_CustomGBForce* target, const char* name) {
    int result = reinterpret_cast<OpenMM::CustomGBForce*>(target)->addPerParticleParameter(std::string(name));
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomGBForce_getPerParticleParameterName(const OpenMM_CustomGBForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomGBForce*>(target)->getPerParticleParameterName(index);
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomGBForce_setPerParticleParameterName(OpenMM_CustomGBForce* target, int index, const char* name) {
    reinterpret_cast<OpenMM::CustomGBForce*>(target)->setPerParticleParameterName(index, std::string(name));
}
OPENMM_EXPORT int OpenMM_CustomGBForce_addGlobalParameter(OpenMM_CustomGBForce* target, const char* name, double defaultValue) {
    int result = reinterpret_cast<OpenMM::CustomGBForce*>(target)->addGlobalParameter(std::string(name), defaultValue);
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomGBForce_getGlobalParameterName(const OpenMM_CustomGBForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomGBForce*>(target)->getGlobalParameterName(index);
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomGBForce_setGlobalParameterName(OpenMM_CustomGBForce* target, int index, const char* name) {
    reinterpret_cast<OpenMM::CustomGBForce*>(target)->setGlobalParameterName(index, std::string(name));
}
OPENMM_EXPORT double OpenMM_CustomGBForce_getGlobalParameterDefaultValue(const OpenMM_CustomGBForce* target, int index) {
    double result = reinterpret_cast<const OpenMM::CustomGBForce*>(target)->getGlobalParameterDefaultValue(index);
    return result;
}
OPENMM_EXPORT void OpenMM_CustomGBForce_setGlobalParameterDefaultValue(OpenMM_CustomGBForce* target, int index, double defaultValue) {
    reinterpret_cast<OpenMM::CustomGBForce*>(target)->setGlobalParameterDefaultValue(index, defaultValue);
}
OPENMM_EXPORT void OpenMM_CustomGBForce_addEnergyParameterDerivative(OpenMM_CustomGBForce* target, const char* name) {
    reinterpret_cast<OpenMM::CustomGBForce*>(target)->addEnergyParameterDerivative(std::string(name));
}
OPENMM_EXPORT const char* OpenMM_CustomGBForce_getEnergyParameterDerivativeName(const OpenMM_CustomGBForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomGBForce*>(target)->getEnergyParameterDerivativeName(index);
    return result->c_str();
}
OPENMM_EXPORT int OpenMM_CustomGBForce_addParticle(OpenMM_CustomGBForce* target, const OpenMM_DoubleArray* parameters) {
    int result = reinterpret_cast<OpenMM::CustomGBForce*>(target)->addParticle(*reinterpret_cast<const std::vector< double >*>(parameters));
    return result;
}
OPENMM_EXPORT void OpenMM_CustomGBForce_getParticleParameters(const OpenMM_CustomGBForce* target, int index, OpenMM_DoubleArray* parameters) {
    reinterpret_cast<const OpenMM::CustomGBForce*>(target)->getParticleParameters(index, *reinterpret_cast<std::vector< double >*>(parameters));
}
OPENMM_EXPORT void OpenMM_CustomGBForce_setParticleParameters(OpenMM_CustomGBForce* target, int index, const OpenMM_DoubleArray* parameters) {
    reinterpret_cast<OpenMM::CustomGBForce*>(target)->setParticleParameters(index, *reinterpret_cast<const std::vector< double >*>(parameters));
}
OPENMM_EXPORT int OpenMM_CustomGBForce_addComputedValue(OpenMM_CustomGBForce* target, const char* name, const char* expression, OpenMM_CustomGBForce_ComputationType type) {
    int result = reinterpret_cast<OpenMM::CustomGBForce*>(target)->addComputedValue(std::string(name), std::string(expression), static_cast<OpenMM::CustomGBForce::ComputationType>(type));
    return result;
}
OPENMM_EXPORT void OpenMM_CustomGBForce_getComputedValueParameters(const OpenMM_CustomGBForce* target, int index, char** name, char** expression, OpenMM_CustomGBForce_ComputationType* type) {
    reinterpret_cast<const OpenMM::CustomGBForce*>(target)->getComputedValueParameters(index, *reinterpret_cast<std::string*>(name), *reinterpret_cast<std::string*>(expression), *reinterpret_cast<OpenMM::CustomGBForce::ComputationType*>(type));
}
OPENMM_EXPORT void OpenMM_CustomGBForce_setComputedValueParameters(OpenMM_CustomGBForce* target, int index, const char* name, const char* expression, OpenMM_CustomGBForce_ComputationType type) {
    reinterpret_cast<OpenMM::CustomGBForce*>(target)->setComputedValueParameters(index, std::string(name), std::string(expression), static_cast<OpenMM::CustomGBForce::ComputationType>(type));
}
OPENMM_EXPORT int OpenMM_CustomGBForce_addEnergyTerm(OpenMM_CustomGBForce* target, const char* expression, OpenMM_CustomGBForce_ComputationType type) {
    int result = reinterpret_cast<OpenMM::CustomGBForce*>(target)->addEnergyTerm(std::string(expression), static_cast<OpenMM::CustomGBForce::ComputationType>(type));
    return result;
}
OPENMM_EXPORT void OpenMM_CustomGBForce_getEnergyTermParameters(const OpenMM_CustomGBForce* target, int index, char** expression, OpenMM_CustomGBForce_ComputationType* type) {
    reinterpret_cast<const OpenMM::CustomGBForce*>(target)->getEnergyTermParameters(index, *reinterpret_cast<std::string*>(expression), *reinterpret_cast<OpenMM::CustomGBForce::ComputationType*>(type));
}
OPENMM_EXPORT void OpenMM_CustomGBForce_setEnergyTermParameters(OpenMM_CustomGBForce* target, int index, const char* expression, OpenMM_CustomGBForce_ComputationType type) {
    reinterpret_cast<OpenMM::CustomGBForce*>(target)->setEnergyTermParameters(index, std::string(expression), static_cast<OpenMM::CustomGBForce::ComputationType>(type));
}
OPENMM_EXPORT int OpenMM_CustomGBForce_addExclusion(OpenMM_CustomGBForce* target, int particle1, int particle2) {
    int result = reinterpret_cast<OpenMM::CustomGBForce*>(target)->addExclusion(particle1, particle2);
    return result;
}
OPENMM_EXPORT void OpenMM_CustomGBForce_getExclusionParticles(const OpenMM_CustomGBForce* target, int index, int* particle1, int* particle2) {
    reinterpret_cast<const OpenMM::CustomGBForce*>(target)->getExclusionParticles(index, *reinterpret_cast<int*>(particle1), *reinterpret_cast<int*>(particle2));
}
OPENMM_EXPORT void OpenMM_CustomGBForce_setExclusionParticles(OpenMM_CustomGBForce* target, int index, int particle1, int particle2) {
    reinterpret_cast<OpenMM::CustomGBForce*>(target)->setExclusionParticles(index, particle1, particle2);
}
OPENMM_EXPORT int OpenMM_CustomGBForce_addTabulatedFunction(OpenMM_CustomGBForce* target, const char* name, OpenMM_TabulatedFunction* function) {
    int result = reinterpret_cast<OpenMM::CustomGBForce*>(target)->addTabulatedFunction(std::string(name), reinterpret_cast<TabulatedFunction *>(function));
    return result;
}
OPENMM_EXPORT OpenMM_TabulatedFunction* OpenMM_CustomGBForce_getTabulatedFunction(OpenMM_CustomGBForce* target, int index) {
    TabulatedFunction* result = &reinterpret_cast<OpenMM::CustomGBForce*>(target)->getTabulatedFunction(index);
    return reinterpret_cast<OpenMM_TabulatedFunction*>(result);
}
OPENMM_EXPORT const char* OpenMM_CustomGBForce_getTabulatedFunctionName(const OpenMM_CustomGBForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomGBForce*>(target)->getTabulatedFunctionName(index);
    return result->c_str();
}
OPENMM_EXPORT int OpenMM_CustomGBForce_addFunction(OpenMM_CustomGBForce* target, const char* name, const OpenMM_DoubleArray* values, double min, double max) {
    int result = reinterpret_cast<OpenMM::CustomGBForce*>(target)->addFunction(std::string(name), *reinterpret_cast<const std::vector< double >*>(values), min, max);
    return result;
}
OPENMM_EXPORT void OpenMM_CustomGBForce_getFunctionParameters(const OpenMM_CustomGBForce* target, int index, char** name, OpenMM_DoubleArray* values, double* min, double* max) {
    reinterpret_cast<const OpenMM::CustomGBForce*>(target)->getFunctionParameters(index, *reinterpret_cast<std::string*>(name), *reinterpret_cast<std::vector< double >*>(values), *reinterpret_cast<double*>(min), *reinterpret_cast<double*>(max));
}
OPENMM_EXPORT void OpenMM_CustomGBForce_setFunctionParameters(OpenMM_CustomGBForce* target, int index, const char* name, const OpenMM_DoubleArray* values, double min, double max) {
    reinterpret_cast<OpenMM::CustomGBForce*>(target)->setFunctionParameters(index, std::string(name), *reinterpret_cast<const std::vector< double >*>(values), min, max);
}
OPENMM_EXPORT void OpenMM_CustomGBForce_updateParametersInContext(OpenMM_CustomGBForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::CustomGBForce*>(target)->updateParametersInContext(*reinterpret_cast<OpenMM::Context*>(context));
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_CustomGBForce_usesPeriodicBoundaryConditions(const OpenMM_CustomGBForce* target) {
    bool result = reinterpret_cast<const OpenMM::CustomGBForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::Discrete3DFunction */
OPENMM_EXPORT OpenMM_Discrete3DFunction* OpenMM_Discrete3DFunction_create(int xsize, int ysize, int zsize, const OpenMM_DoubleArray* values) {
    return reinterpret_cast<OpenMM_Discrete3DFunction*>(new OpenMM::Discrete3DFunction(xsize, ysize, zsize, *reinterpret_cast<const std::vector< double >*>(values)));
}
OPENMM_EXPORT void OpenMM_Discrete3DFunction_destroy(OpenMM_Discrete3DFunction* target) {
    delete reinterpret_cast<OpenMM::Discrete3DFunction*>(target);
}
OPENMM_EXPORT void OpenMM_Discrete3DFunction_getFunctionParameters(const OpenMM_Discrete3DFunction* target, int* xsize, int* ysize, int* zsize, OpenMM_DoubleArray* values) {
    reinterpret_cast<const OpenMM::Discrete3DFunction*>(target)->getFunctionParameters(*reinterpret_cast<int*>(xsize), *reinterpret_cast<int*>(ysize), *reinterpret_cast<int*>(zsize), *reinterpret_cast<std::vector< double >*>(values));
}
OPENMM_EXPORT void OpenMM_Discrete3DFunction_setFunctionParameters(OpenMM_Discrete3DFunction* target, int xsize, int ysize, int zsize, const OpenMM_DoubleArray* values) {
    reinterpret_cast<OpenMM::Discrete3DFunction*>(target)->setFunctionParameters(xsize, ysize, zsize, *reinterpret_cast<const std::vector< double >*>(values));
}
OPENMM_EXPORT OpenMM_Discrete3DFunction* OpenMM_Discrete3DFunction_Copy(const OpenMM_Discrete3DFunction* target) {
    Discrete3DFunction * result = reinterpret_cast<const OpenMM::Discrete3DFunction*>(target)->Copy();
    return reinterpret_cast<OpenMM_Discrete3DFunction*>(result);
}

/* OpenMM::PeriodicTorsionForce */
OPENMM_EXPORT OpenMM_PeriodicTorsionForce* OpenMM_PeriodicTorsionForce_create() {
    return reinterpret_cast<OpenMM_PeriodicTorsionForce*>(new OpenMM::PeriodicTorsionForce());
}
OPENMM_EXPORT void OpenMM_PeriodicTorsionForce_destroy(OpenMM_PeriodicTorsionForce* target) {
    delete reinterpret_cast<OpenMM::PeriodicTorsionForce*>(target);
}
OPENMM_EXPORT int OpenMM_PeriodicTorsionForce_getNumTorsions(const OpenMM_PeriodicTorsionForce* target) {
    int result = reinterpret_cast<const OpenMM::PeriodicTorsionForce*>(target)->getNumTorsions();
    return result;
}
OPENMM_EXPORT int OpenMM_PeriodicTorsionForce_addTorsion(OpenMM_PeriodicTorsionForce* target, int particle1, int particle2, int particle3, int particle4, int periodicity, double phase, double k) {
    int result = reinterpret_cast<OpenMM::PeriodicTorsionForce*>(target)->addTorsion(particle1, particle2, particle3, particle4, periodicity, phase, k);
    return result;
}
OPENMM_EXPORT void OpenMM_PeriodicTorsionForce_getTorsionParameters(const OpenMM_PeriodicTorsionForce* target, int index, int* particle1, int* particle2, int* particle3, int* particle4, int* periodicity, double* phase, double* k) {
    reinterpret_cast<const OpenMM::PeriodicTorsionForce*>(target)->getTorsionParameters(index, *reinterpret_cast<int*>(particle1), *reinterpret_cast<int*>(particle2), *reinterpret_cast<int*>(particle3), *reinterpret_cast<int*>(particle4), *reinterpret_cast<int*>(periodicity), *reinterpret_cast<double*>(phase), *reinterpret_cast<double*>(k));
}
OPENMM_EXPORT void OpenMM_PeriodicTorsionForce_setTorsionParameters(OpenMM_PeriodicTorsionForce* target, int index, int particle1, int particle2, int particle3, int particle4, int periodicity, double phase, double k) {
    reinterpret_cast<OpenMM::PeriodicTorsionForce*>(target)->setTorsionParameters(index, particle1, particle2, particle3, particle4, periodicity, phase, k);
}
OPENMM_EXPORT void OpenMM_PeriodicTorsionForce_updateParametersInContext(OpenMM_PeriodicTorsionForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::PeriodicTorsionForce*>(target)->updateParametersInContext(*reinterpret_cast<OpenMM::Context*>(context));
}
OPENMM_EXPORT void OpenMM_PeriodicTorsionForce_setUsesPeriodicBoundaryConditions(OpenMM_PeriodicTorsionForce* target, OpenMM_Boolean periodic) {
    reinterpret_cast<OpenMM::PeriodicTorsionForce*>(target)->setUsesPeriodicBoundaryConditions(periodic);
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_PeriodicTorsionForce_usesPeriodicBoundaryConditions(const OpenMM_PeriodicTorsionForce* target) {
    bool result = reinterpret_cast<const OpenMM::PeriodicTorsionForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::LocalCoordinatesSite */
OPENMM_EXPORT OpenMM_LocalCoordinatesSite* OpenMM_LocalCoordinatesSite_create(const OpenMM_IntArray* particles, const OpenMM_DoubleArray* originWeights, const OpenMM_DoubleArray* xWeights, const OpenMM_DoubleArray* yWeights, const OpenMM_Vec3* localPosition) {
    return reinterpret_cast<OpenMM_LocalCoordinatesSite*>(new OpenMM::LocalCoordinatesSite(*reinterpret_cast<const std::vector< int >*>(particles), *reinterpret_cast<const std::vector< double >*>(originWeights), *reinterpret_cast<const std::vector< double >*>(xWeights), *reinterpret_cast<const std::vector< double >*>(yWeights), *reinterpret_cast<const Vec3*>(localPosition)));
}
OPENMM_EXPORT OpenMM_LocalCoordinatesSite* OpenMM_LocalCoordinatesSite_create_2(int particle1, int particle2, int particle3, const OpenMM_Vec3* originWeights, const OpenMM_Vec3* xWeights, const OpenMM_Vec3* yWeights, const OpenMM_Vec3* localPosition) {
    return reinterpret_cast<OpenMM_LocalCoordinatesSite*>(new OpenMM::LocalCoordinatesSite(particle1, particle2, particle3, *reinterpret_cast<const Vec3*>(originWeights), *reinterpret_cast<const Vec3*>(xWeights), *reinterpret_cast<const Vec3*>(yWeights), *reinterpret_cast<const Vec3*>(localPosition)));
}
OPENMM_EXPORT void OpenMM_LocalCoordinatesSite_destroy(OpenMM_LocalCoordinatesSite* target) {
    delete reinterpret_cast<OpenMM::LocalCoordinatesSite*>(target);
}
OPENMM_EXPORT void OpenMM_LocalCoordinatesSite_getOriginWeights(const OpenMM_LocalCoordinatesSite* target, OpenMM_DoubleArray* weights) {
    reinterpret_cast<const OpenMM::LocalCoordinatesSite*>(target)->getOriginWeights(*reinterpret_cast<std::vector< double >*>(weights));
}
OPENMM_EXPORT void OpenMM_LocalCoordinatesSite_getXWeights(const OpenMM_LocalCoordinatesSite* target, OpenMM_DoubleArray* weights) {
    reinterpret_cast<const OpenMM::LocalCoordinatesSite*>(target)->getXWeights(*reinterpret_cast<std::vector< double >*>(weights));
}
OPENMM_EXPORT void OpenMM_LocalCoordinatesSite_getYWeights(const OpenMM_LocalCoordinatesSite* target, OpenMM_DoubleArray* weights) {
    reinterpret_cast<const OpenMM::LocalCoordinatesSite*>(target)->getYWeights(*reinterpret_cast<std::vector< double >*>(weights));
}
OPENMM_EXPORT const OpenMM_Vec3* OpenMM_LocalCoordinatesSite_getLocalPosition(const OpenMM_LocalCoordinatesSite* target) {
    const Vec3* result = &reinterpret_cast<const OpenMM::LocalCoordinatesSite*>(target)->getLocalPosition();
    return reinterpret_cast<const OpenMM_Vec3*>(result);
}

/* OpenMM::State */
OPENMM_EXPORT OpenMM_State* OpenMM_State_create() {
    return reinterpret_cast<OpenMM_State*>(new OpenMM::State());
}
OPENMM_EXPORT void OpenMM_State_destroy(OpenMM_State* target) {
    delete reinterpret_cast<OpenMM::State*>(target);
}
OPENMM_EXPORT double OpenMM_State_getTime(const OpenMM_State* target) {
    double result = reinterpret_cast<const OpenMM::State*>(target)->getTime();
    return result;
}
OPENMM_EXPORT const OpenMM_Vec3Array* OpenMM_State_getPositions(const OpenMM_State* target) {
    const std::vector< Vec3 >* result = &reinterpret_cast<const OpenMM::State*>(target)->getPositions();
    return reinterpret_cast<const OpenMM_Vec3Array*>(result);
}
OPENMM_EXPORT const OpenMM_Vec3Array* OpenMM_State_getVelocities(const OpenMM_State* target) {
    const std::vector< Vec3 >* result = &reinterpret_cast<const OpenMM::State*>(target)->getVelocities();
    return reinterpret_cast<const OpenMM_Vec3Array*>(result);
}
OPENMM_EXPORT const OpenMM_Vec3Array* OpenMM_State_getForces(const OpenMM_State* target) {
    const std::vector< Vec3 >* result = &reinterpret_cast<const OpenMM::State*>(target)->getForces();
    return reinterpret_cast<const OpenMM_Vec3Array*>(result);
}
OPENMM_EXPORT double OpenMM_State_getKineticEnergy(const OpenMM_State* target) {
    double result = reinterpret_cast<const OpenMM::State*>(target)->getKineticEnergy();
    return result;
}
OPENMM_EXPORT double OpenMM_State_getPotentialEnergy(const OpenMM_State* target) {
    double result = reinterpret_cast<const OpenMM::State*>(target)->getPotentialEnergy();
    return result;
}
OPENMM_EXPORT void OpenMM_State_getPeriodicBoxVectors(const OpenMM_State* target, OpenMM_Vec3* a, OpenMM_Vec3* b, OpenMM_Vec3* c) {
    reinterpret_cast<const OpenMM::State*>(target)->getPeriodicBoxVectors(*reinterpret_cast<Vec3*>(a), *reinterpret_cast<Vec3*>(b), *reinterpret_cast<Vec3*>(c));
}
OPENMM_EXPORT double OpenMM_State_getPeriodicBoxVolume(const OpenMM_State* target) {
    double result = reinterpret_cast<const OpenMM::State*>(target)->getPeriodicBoxVolume();
    return result;
}
OPENMM_EXPORT const OpenMM_ParameterArray* OpenMM_State_getParameters(const OpenMM_State* target) {
    const std::map< std::string, double >* result = &reinterpret_cast<const OpenMM::State*>(target)->getParameters();
    return reinterpret_cast<const OpenMM_ParameterArray*>(result);
}
OPENMM_EXPORT const OpenMM_ParameterArray* OpenMM_State_getEnergyParameterDerivatives(const OpenMM_State* target) {
    const std::map< std::string, double >* result = &reinterpret_cast<const OpenMM::State*>(target)->getEnergyParameterDerivatives();
    return reinterpret_cast<const OpenMM_ParameterArray*>(result);
}
OPENMM_EXPORT int OpenMM_State_getDataTypes(const OpenMM_State* target) {
    int result = reinterpret_cast<const OpenMM::State*>(target)->getDataTypes();
    return result;
}

/* OpenMM::OpenMMException */
OPENMM_EXPORT OpenMM_OpenMMException* OpenMM_OpenMMException_create(const char* message) {
    return reinterpret_cast<OpenMM_OpenMMException*>(new OpenMM::OpenMMException(std::string(message)));
}
OPENMM_EXPORT void OpenMM_OpenMMException_destroy(OpenMM_OpenMMException* target) {
    delete reinterpret_cast<OpenMM::OpenMMException*>(target);
}
OPENMM_EXPORT const char* OpenMM_OpenMMException_what(const OpenMM_OpenMMException* target) {
    const char * result = reinterpret_cast<const OpenMM::OpenMMException*>(target)->what();
    return reinterpret_cast<const char*>(result);
}

/* OpenMM::CMAPTorsionForce */
OPENMM_EXPORT OpenMM_CMAPTorsionForce* OpenMM_CMAPTorsionForce_create() {
    return reinterpret_cast<OpenMM_CMAPTorsionForce*>(new OpenMM::CMAPTorsionForce());
}
OPENMM_EXPORT void OpenMM_CMAPTorsionForce_destroy(OpenMM_CMAPTorsionForce* target) {
    delete reinterpret_cast<OpenMM::CMAPTorsionForce*>(target);
}
OPENMM_EXPORT int OpenMM_CMAPTorsionForce_getNumMaps(const OpenMM_CMAPTorsionForce* target) {
    int result = reinterpret_cast<const OpenMM::CMAPTorsionForce*>(target)->getNumMaps();
    return result;
}
OPENMM_EXPORT int OpenMM_CMAPTorsionForce_getNumTorsions(const OpenMM_CMAPTorsionForce* target) {
    int result = reinterpret_cast<const OpenMM::CMAPTorsionForce*>(target)->getNumTorsions();
    return result;
}
OPENMM_EXPORT int OpenMM_CMAPTorsionForce_addMap(OpenMM_CMAPTorsionForce* target, int size, const OpenMM_DoubleArray* energy) {
    int result = reinterpret_cast<OpenMM::CMAPTorsionForce*>(target)->addMap(size, *reinterpret_cast<const std::vector< double >*>(energy));
    return result;
}
OPENMM_EXPORT void OpenMM_CMAPTorsionForce_getMapParameters(const OpenMM_CMAPTorsionForce* target, int index, int* size, OpenMM_DoubleArray* energy) {
    reinterpret_cast<const OpenMM::CMAPTorsionForce*>(target)->getMapParameters(index, *reinterpret_cast<int*>(size), *reinterpret_cast<std::vector< double >*>(energy));
}
OPENMM_EXPORT void OpenMM_CMAPTorsionForce_setMapParameters(OpenMM_CMAPTorsionForce* target, int index, int size, const OpenMM_DoubleArray* energy) {
    reinterpret_cast<OpenMM::CMAPTorsionForce*>(target)->setMapParameters(index, size, *reinterpret_cast<const std::vector< double >*>(energy));
}
OPENMM_EXPORT int OpenMM_CMAPTorsionForce_addTorsion(OpenMM_CMAPTorsionForce* target, int map, int a1, int a2, int a3, int a4, int b1, int b2, int b3, int b4) {
    int result = reinterpret_cast<OpenMM::CMAPTorsionForce*>(target)->addTorsion(map, a1, a2, a3, a4, b1, b2, b3, b4);
    return result;
}
OPENMM_EXPORT void OpenMM_CMAPTorsionForce_getTorsionParameters(const OpenMM_CMAPTorsionForce* target, int index, int* map, int* a1, int* a2, int* a3, int* a4, int* b1, int* b2, int* b3, int* b4) {
    reinterpret_cast<const OpenMM::CMAPTorsionForce*>(target)->getTorsionParameters(index, *reinterpret_cast<int*>(map), *reinterpret_cast<int*>(a1), *reinterpret_cast<int*>(a2), *reinterpret_cast<int*>(a3), *reinterpret_cast<int*>(a4), *reinterpret_cast<int*>(b1), *reinterpret_cast<int*>(b2), *reinterpret_cast<int*>(b3), *reinterpret_cast<int*>(b4));
}
OPENMM_EXPORT void OpenMM_CMAPTorsionForce_setTorsionParameters(OpenMM_CMAPTorsionForce* target, int index, int map, int a1, int a2, int a3, int a4, int b1, int b2, int b3, int b4) {
    reinterpret_cast<OpenMM::CMAPTorsionForce*>(target)->setTorsionParameters(index, map, a1, a2, a3, a4, b1, b2, b3, b4);
}
OPENMM_EXPORT void OpenMM_CMAPTorsionForce_updateParametersInContext(OpenMM_CMAPTorsionForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::CMAPTorsionForce*>(target)->updateParametersInContext(*reinterpret_cast<OpenMM::Context*>(context));
}
OPENMM_EXPORT void OpenMM_CMAPTorsionForce_setUsesPeriodicBoundaryConditions(OpenMM_CMAPTorsionForce* target, OpenMM_Boolean periodic) {
    reinterpret_cast<OpenMM::CMAPTorsionForce*>(target)->setUsesPeriodicBoundaryConditions(periodic);
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_CMAPTorsionForce_usesPeriodicBoundaryConditions(const OpenMM_CMAPTorsionForce* target) {
    bool result = reinterpret_cast<const OpenMM::CMAPTorsionForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::Discrete1DFunction */
OPENMM_EXPORT OpenMM_Discrete1DFunction* OpenMM_Discrete1DFunction_create(const OpenMM_DoubleArray* values) {
    return reinterpret_cast<OpenMM_Discrete1DFunction*>(new OpenMM::Discrete1DFunction(*reinterpret_cast<const std::vector< double >*>(values)));
}
OPENMM_EXPORT void OpenMM_Discrete1DFunction_destroy(OpenMM_Discrete1DFunction* target) {
    delete reinterpret_cast<OpenMM::Discrete1DFunction*>(target);
}
OPENMM_EXPORT void OpenMM_Discrete1DFunction_getFunctionParameters(const OpenMM_Discrete1DFunction* target, OpenMM_DoubleArray* values) {
    reinterpret_cast<const OpenMM::Discrete1DFunction*>(target)->getFunctionParameters(*reinterpret_cast<std::vector< double >*>(values));
}
OPENMM_EXPORT void OpenMM_Discrete1DFunction_setFunctionParameters(OpenMM_Discrete1DFunction* target, const OpenMM_DoubleArray* values) {
    reinterpret_cast<OpenMM::Discrete1DFunction*>(target)->setFunctionParameters(*reinterpret_cast<const std::vector< double >*>(values));
}
OPENMM_EXPORT OpenMM_Discrete1DFunction* OpenMM_Discrete1DFunction_Copy(const OpenMM_Discrete1DFunction* target) {
    Discrete1DFunction * result = reinterpret_cast<const OpenMM::Discrete1DFunction*>(target)->Copy();
    return reinterpret_cast<OpenMM_Discrete1DFunction*>(result);
}

/* OpenMM::CustomManyParticleForce */
OPENMM_EXPORT OpenMM_CustomManyParticleForce* OpenMM_CustomManyParticleForce_create(int particlesPerSet, const char* energy) {
    return reinterpret_cast<OpenMM_CustomManyParticleForce*>(new OpenMM::CustomManyParticleForce(particlesPerSet, std::string(energy)));
}
OPENMM_EXPORT void OpenMM_CustomManyParticleForce_destroy(OpenMM_CustomManyParticleForce* target) {
    delete reinterpret_cast<OpenMM::CustomManyParticleForce*>(target);
}
OPENMM_EXPORT int OpenMM_CustomManyParticleForce_getNumParticlesPerSet(const OpenMM_CustomManyParticleForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomManyParticleForce*>(target)->getNumParticlesPerSet();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomManyParticleForce_getNumParticles(const OpenMM_CustomManyParticleForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomManyParticleForce*>(target)->getNumParticles();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomManyParticleForce_getNumExclusions(const OpenMM_CustomManyParticleForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomManyParticleForce*>(target)->getNumExclusions();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomManyParticleForce_getNumPerParticleParameters(const OpenMM_CustomManyParticleForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomManyParticleForce*>(target)->getNumPerParticleParameters();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomManyParticleForce_getNumGlobalParameters(const OpenMM_CustomManyParticleForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomManyParticleForce*>(target)->getNumGlobalParameters();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomManyParticleForce_getNumTabulatedFunctions(const OpenMM_CustomManyParticleForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomManyParticleForce*>(target)->getNumTabulatedFunctions();
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomManyParticleForce_getEnergyFunction(const OpenMM_CustomManyParticleForce* target) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomManyParticleForce*>(target)->getEnergyFunction();
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomManyParticleForce_setEnergyFunction(OpenMM_CustomManyParticleForce* target, const char* energy) {
    reinterpret_cast<OpenMM::CustomManyParticleForce*>(target)->setEnergyFunction(std::string(energy));
}
OPENMM_EXPORT OpenMM_CustomManyParticleForce_NonbondedMethod OpenMM_CustomManyParticleForce_getNonbondedMethod(const OpenMM_CustomManyParticleForce* target) {
    OpenMM::CustomManyParticleForce::NonbondedMethod result = reinterpret_cast<const OpenMM::CustomManyParticleForce*>(target)->getNonbondedMethod();
    return static_cast<OpenMM_CustomManyParticleForce_NonbondedMethod>(result);
}
OPENMM_EXPORT void OpenMM_CustomManyParticleForce_setNonbondedMethod(OpenMM_CustomManyParticleForce* target, OpenMM_CustomManyParticleForce_NonbondedMethod method) {
    reinterpret_cast<OpenMM::CustomManyParticleForce*>(target)->setNonbondedMethod(static_cast<OpenMM::CustomManyParticleForce::NonbondedMethod>(method));
}
OPENMM_EXPORT OpenMM_CustomManyParticleForce_PermutationMode OpenMM_CustomManyParticleForce_getPermutationMode(const OpenMM_CustomManyParticleForce* target) {
    OpenMM::CustomManyParticleForce::PermutationMode result = reinterpret_cast<const OpenMM::CustomManyParticleForce*>(target)->getPermutationMode();
    return static_cast<OpenMM_CustomManyParticleForce_PermutationMode>(result);
}
OPENMM_EXPORT void OpenMM_CustomManyParticleForce_setPermutationMode(OpenMM_CustomManyParticleForce* target, OpenMM_CustomManyParticleForce_PermutationMode mode) {
    reinterpret_cast<OpenMM::CustomManyParticleForce*>(target)->setPermutationMode(static_cast<OpenMM::CustomManyParticleForce::PermutationMode>(mode));
}
OPENMM_EXPORT double OpenMM_CustomManyParticleForce_getCutoffDistance(const OpenMM_CustomManyParticleForce* target) {
    double result = reinterpret_cast<const OpenMM::CustomManyParticleForce*>(target)->getCutoffDistance();
    return result;
}
OPENMM_EXPORT void OpenMM_CustomManyParticleForce_setCutoffDistance(OpenMM_CustomManyParticleForce* target, double distance) {
    reinterpret_cast<OpenMM::CustomManyParticleForce*>(target)->setCutoffDistance(distance);
}
OPENMM_EXPORT int OpenMM_CustomManyParticleForce_addPerParticleParameter(OpenMM_CustomManyParticleForce* target, const char* name) {
    int result = reinterpret_cast<OpenMM::CustomManyParticleForce*>(target)->addPerParticleParameter(std::string(name));
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomManyParticleForce_getPerParticleParameterName(const OpenMM_CustomManyParticleForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomManyParticleForce*>(target)->getPerParticleParameterName(index);
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomManyParticleForce_setPerParticleParameterName(OpenMM_CustomManyParticleForce* target, int index, const char* name) {
    reinterpret_cast<OpenMM::CustomManyParticleForce*>(target)->setPerParticleParameterName(index, std::string(name));
}
OPENMM_EXPORT int OpenMM_CustomManyParticleForce_addGlobalParameter(OpenMM_CustomManyParticleForce* target, const char* name, double defaultValue) {
    int result = reinterpret_cast<OpenMM::CustomManyParticleForce*>(target)->addGlobalParameter(std::string(name), defaultValue);
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomManyParticleForce_getGlobalParameterName(const OpenMM_CustomManyParticleForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomManyParticleForce*>(target)->getGlobalParameterName(index);
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomManyParticleForce_setGlobalParameterName(OpenMM_CustomManyParticleForce* target, int index, const char* name) {
    reinterpret_cast<OpenMM::CustomManyParticleForce*>(target)->setGlobalParameterName(index, std::string(name));
}
OPENMM_EXPORT double OpenMM_CustomManyParticleForce_getGlobalParameterDefaultValue(const OpenMM_CustomManyParticleForce* target, int index) {
    double result = reinterpret_cast<const OpenMM::CustomManyParticleForce*>(target)->getGlobalParameterDefaultValue(index);
    return result;
}
OPENMM_EXPORT void OpenMM_CustomManyParticleForce_setGlobalParameterDefaultValue(OpenMM_CustomManyParticleForce* target, int index, double defaultValue) {
    reinterpret_cast<OpenMM::CustomManyParticleForce*>(target)->setGlobalParameterDefaultValue(index, defaultValue);
}
OPENMM_EXPORT int OpenMM_CustomManyParticleForce_addParticle(OpenMM_CustomManyParticleForce* target, const OpenMM_DoubleArray* parameters, int type) {
    int result = reinterpret_cast<OpenMM::CustomManyParticleForce*>(target)->addParticle(*reinterpret_cast<const std::vector< double >*>(parameters), type);
    return result;
}
OPENMM_EXPORT void OpenMM_CustomManyParticleForce_getParticleParameters(const OpenMM_CustomManyParticleForce* target, int index, OpenMM_DoubleArray* parameters, int* type) {
    reinterpret_cast<const OpenMM::CustomManyParticleForce*>(target)->getParticleParameters(index, *reinterpret_cast<std::vector< double >*>(parameters), *reinterpret_cast<int*>(type));
}
OPENMM_EXPORT void OpenMM_CustomManyParticleForce_setParticleParameters(OpenMM_CustomManyParticleForce* target, int index, const OpenMM_DoubleArray* parameters, int type) {
    reinterpret_cast<OpenMM::CustomManyParticleForce*>(target)->setParticleParameters(index, *reinterpret_cast<const std::vector< double >*>(parameters), type);
}
OPENMM_EXPORT int OpenMM_CustomManyParticleForce_addExclusion(OpenMM_CustomManyParticleForce* target, int particle1, int particle2) {
    int result = reinterpret_cast<OpenMM::CustomManyParticleForce*>(target)->addExclusion(particle1, particle2);
    return result;
}
OPENMM_EXPORT void OpenMM_CustomManyParticleForce_getExclusionParticles(const OpenMM_CustomManyParticleForce* target, int index, int* particle1, int* particle2) {
    reinterpret_cast<const OpenMM::CustomManyParticleForce*>(target)->getExclusionParticles(index, *reinterpret_cast<int*>(particle1), *reinterpret_cast<int*>(particle2));
}
OPENMM_EXPORT void OpenMM_CustomManyParticleForce_setExclusionParticles(OpenMM_CustomManyParticleForce* target, int index, int particle1, int particle2) {
    reinterpret_cast<OpenMM::CustomManyParticleForce*>(target)->setExclusionParticles(index, particle1, particle2);
}
OPENMM_EXPORT void OpenMM_CustomManyParticleForce_createExclusionsFromBonds(OpenMM_CustomManyParticleForce* target, const OpenMM_BondArray* bonds, int bondCutoff) {
    reinterpret_cast<OpenMM::CustomManyParticleForce*>(target)->createExclusionsFromBonds(*reinterpret_cast<const std::vector< std::pair< int, int > >*>(bonds), bondCutoff);
}
OPENMM_EXPORT void OpenMM_CustomManyParticleForce_getTypeFilter(const OpenMM_CustomManyParticleForce* target, int index, OpenMM_IntSet* types) {
    reinterpret_cast<const OpenMM::CustomManyParticleForce*>(target)->getTypeFilter(index, *reinterpret_cast<std::set< int >*>(types));
}
OPENMM_EXPORT void OpenMM_CustomManyParticleForce_setTypeFilter(OpenMM_CustomManyParticleForce* target, int index, const OpenMM_IntSet* types) {
    reinterpret_cast<OpenMM::CustomManyParticleForce*>(target)->setTypeFilter(index, *reinterpret_cast<const std::set< int >*>(types));
}
OPENMM_EXPORT int OpenMM_CustomManyParticleForce_addTabulatedFunction(OpenMM_CustomManyParticleForce* target, const char* name, OpenMM_TabulatedFunction* function) {
    int result = reinterpret_cast<OpenMM::CustomManyParticleForce*>(target)->addTabulatedFunction(std::string(name), reinterpret_cast<TabulatedFunction *>(function));
    return result;
}
OPENMM_EXPORT OpenMM_TabulatedFunction* OpenMM_CustomManyParticleForce_getTabulatedFunction(OpenMM_CustomManyParticleForce* target, int index) {
    TabulatedFunction* result = &reinterpret_cast<OpenMM::CustomManyParticleForce*>(target)->getTabulatedFunction(index);
    return reinterpret_cast<OpenMM_TabulatedFunction*>(result);
}
OPENMM_EXPORT const char* OpenMM_CustomManyParticleForce_getTabulatedFunctionName(const OpenMM_CustomManyParticleForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomManyParticleForce*>(target)->getTabulatedFunctionName(index);
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomManyParticleForce_updateParametersInContext(OpenMM_CustomManyParticleForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::CustomManyParticleForce*>(target)->updateParametersInContext(*reinterpret_cast<OpenMM::Context*>(context));
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_CustomManyParticleForce_usesPeriodicBoundaryConditions(const OpenMM_CustomManyParticleForce* target) {
    bool result = reinterpret_cast<const OpenMM::CustomManyParticleForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::HarmonicAngleForce */
OPENMM_EXPORT OpenMM_HarmonicAngleForce* OpenMM_HarmonicAngleForce_create() {
    return reinterpret_cast<OpenMM_HarmonicAngleForce*>(new OpenMM::HarmonicAngleForce());
}
OPENMM_EXPORT void OpenMM_HarmonicAngleForce_destroy(OpenMM_HarmonicAngleForce* target) {
    delete reinterpret_cast<OpenMM::HarmonicAngleForce*>(target);
}
OPENMM_EXPORT int OpenMM_HarmonicAngleForce_getNumAngles(const OpenMM_HarmonicAngleForce* target) {
    int result = reinterpret_cast<const OpenMM::HarmonicAngleForce*>(target)->getNumAngles();
    return result;
}
OPENMM_EXPORT int OpenMM_HarmonicAngleForce_addAngle(OpenMM_HarmonicAngleForce* target, int particle1, int particle2, int particle3, double angle, double k) {
    int result = reinterpret_cast<OpenMM::HarmonicAngleForce*>(target)->addAngle(particle1, particle2, particle3, angle, k);
    return result;
}
OPENMM_EXPORT void OpenMM_HarmonicAngleForce_getAngleParameters(const OpenMM_HarmonicAngleForce* target, int index, int* particle1, int* particle2, int* particle3, double* angle, double* k) {
    reinterpret_cast<const OpenMM::HarmonicAngleForce*>(target)->getAngleParameters(index, *reinterpret_cast<int*>(particle1), *reinterpret_cast<int*>(particle2), *reinterpret_cast<int*>(particle3), *reinterpret_cast<double*>(angle), *reinterpret_cast<double*>(k));
}
OPENMM_EXPORT void OpenMM_HarmonicAngleForce_setAngleParameters(OpenMM_HarmonicAngleForce* target, int index, int particle1, int particle2, int particle3, double angle, double k) {
    reinterpret_cast<OpenMM::HarmonicAngleForce*>(target)->setAngleParameters(index, particle1, particle2, particle3, angle, k);
}
OPENMM_EXPORT void OpenMM_HarmonicAngleForce_updateParametersInContext(OpenMM_HarmonicAngleForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::HarmonicAngleForce*>(target)->updateParametersInContext(*reinterpret_cast<OpenMM::Context*>(context));
}
OPENMM_EXPORT void OpenMM_HarmonicAngleForce_setUsesPeriodicBoundaryConditions(OpenMM_HarmonicAngleForce* target, OpenMM_Boolean periodic) {
    reinterpret_cast<OpenMM::HarmonicAngleForce*>(target)->setUsesPeriodicBoundaryConditions(periodic);
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_HarmonicAngleForce_usesPeriodicBoundaryConditions(const OpenMM_HarmonicAngleForce* target) {
    bool result = reinterpret_cast<const OpenMM::HarmonicAngleForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::CompoundIntegrator */
OPENMM_EXPORT OpenMM_CompoundIntegrator* OpenMM_CompoundIntegrator_create() {
    return reinterpret_cast<OpenMM_CompoundIntegrator*>(new OpenMM::CompoundIntegrator());
}
OPENMM_EXPORT void OpenMM_CompoundIntegrator_destroy(OpenMM_CompoundIntegrator* target) {
    delete reinterpret_cast<OpenMM::CompoundIntegrator*>(target);
}
OPENMM_EXPORT int OpenMM_CompoundIntegrator_getNumIntegrators(const OpenMM_CompoundIntegrator* target) {
    int result = reinterpret_cast<const OpenMM::CompoundIntegrator*>(target)->getNumIntegrators();
    return result;
}
OPENMM_EXPORT int OpenMM_CompoundIntegrator_addIntegrator(OpenMM_CompoundIntegrator* target, OpenMM_Integrator* integrator) {
    int result = reinterpret_cast<OpenMM::CompoundIntegrator*>(target)->addIntegrator(reinterpret_cast<Integrator *>(integrator));
    return result;
}
OPENMM_EXPORT OpenMM_Integrator* OpenMM_CompoundIntegrator_getIntegrator(OpenMM_CompoundIntegrator* target, int index) {
    Integrator* result = &reinterpret_cast<OpenMM::CompoundIntegrator*>(target)->getIntegrator(index);
    return reinterpret_cast<OpenMM_Integrator*>(result);
}
OPENMM_EXPORT int OpenMM_CompoundIntegrator_getCurrentIntegrator(const OpenMM_CompoundIntegrator* target) {
    int result = reinterpret_cast<const OpenMM::CompoundIntegrator*>(target)->getCurrentIntegrator();
    return result;
}
OPENMM_EXPORT void OpenMM_CompoundIntegrator_setCurrentIntegrator(OpenMM_CompoundIntegrator* target, int index) {
    reinterpret_cast<OpenMM::CompoundIntegrator*>(target)->setCurrentIntegrator(index);
}
OPENMM_EXPORT double OpenMM_CompoundIntegrator_getStepSize(const OpenMM_CompoundIntegrator* target) {
    double result = reinterpret_cast<const OpenMM::CompoundIntegrator*>(target)->getStepSize();
    return result;
}
OPENMM_EXPORT void OpenMM_CompoundIntegrator_setStepSize(OpenMM_CompoundIntegrator* target, double size) {
    reinterpret_cast<OpenMM::CompoundIntegrator*>(target)->setStepSize(size);
}
OPENMM_EXPORT double OpenMM_CompoundIntegrator_getConstraintTolerance(const OpenMM_CompoundIntegrator* target) {
    double result = reinterpret_cast<const OpenMM::CompoundIntegrator*>(target)->getConstraintTolerance();
    return result;
}
OPENMM_EXPORT void OpenMM_CompoundIntegrator_setConstraintTolerance(OpenMM_CompoundIntegrator* target, double tol) {
    reinterpret_cast<OpenMM::CompoundIntegrator*>(target)->setConstraintTolerance(tol);
}
OPENMM_EXPORT void OpenMM_CompoundIntegrator_step(OpenMM_CompoundIntegrator* target, int steps) {
    reinterpret_cast<OpenMM::CompoundIntegrator*>(target)->step(steps);
}

/* OpenMM::MonteCarloAnisotropicBarostat */
OPENMM_EXPORT OpenMM_MonteCarloAnisotropicBarostat* OpenMM_MonteCarloAnisotropicBarostat_create(const OpenMM_Vec3* defaultPressure, double defaultTemperature, OpenMM_Boolean scaleX, OpenMM_Boolean scaleY, OpenMM_Boolean scaleZ, int frequency) {
    return reinterpret_cast<OpenMM_MonteCarloAnisotropicBarostat*>(new OpenMM::MonteCarloAnisotropicBarostat(*reinterpret_cast<const Vec3*>(defaultPressure), defaultTemperature, scaleX, scaleY, scaleZ, frequency));
}
OPENMM_EXPORT void OpenMM_MonteCarloAnisotropicBarostat_destroy(OpenMM_MonteCarloAnisotropicBarostat* target) {
    delete reinterpret_cast<OpenMM::MonteCarloAnisotropicBarostat*>(target);
}
OPENMM_EXPORT const char* OpenMM_MonteCarloAnisotropicBarostat_PressureX() {
    const std::string* result = &OpenMM::MonteCarloAnisotropicBarostat::PressureX();
    return result->c_str();
}
OPENMM_EXPORT const char* OpenMM_MonteCarloAnisotropicBarostat_PressureY() {
    const std::string* result = &OpenMM::MonteCarloAnisotropicBarostat::PressureY();
    return result->c_str();
}
OPENMM_EXPORT const char* OpenMM_MonteCarloAnisotropicBarostat_PressureZ() {
    const std::string* result = &OpenMM::MonteCarloAnisotropicBarostat::PressureZ();
    return result->c_str();
}
OPENMM_EXPORT const char* OpenMM_MonteCarloAnisotropicBarostat_Temperature() {
    const std::string* result = &OpenMM::MonteCarloAnisotropicBarostat::Temperature();
    return result->c_str();
}
OPENMM_EXPORT const OpenMM_Vec3* OpenMM_MonteCarloAnisotropicBarostat_getDefaultPressure(const OpenMM_MonteCarloAnisotropicBarostat* target) {
    const Vec3* result = &reinterpret_cast<const OpenMM::MonteCarloAnisotropicBarostat*>(target)->getDefaultPressure();
    return reinterpret_cast<const OpenMM_Vec3*>(result);
}
OPENMM_EXPORT void OpenMM_MonteCarloAnisotropicBarostat_setDefaultPressure(OpenMM_MonteCarloAnisotropicBarostat* target, const OpenMM_Vec3* pressure) {
    reinterpret_cast<OpenMM::MonteCarloAnisotropicBarostat*>(target)->setDefaultPressure(*reinterpret_cast<const Vec3*>(pressure));
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_MonteCarloAnisotropicBarostat_getScaleX(const OpenMM_MonteCarloAnisotropicBarostat* target) {
    bool result = reinterpret_cast<const OpenMM::MonteCarloAnisotropicBarostat*>(target)->getScaleX();
    return (result ? OpenMM_True : OpenMM_False);
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_MonteCarloAnisotropicBarostat_getScaleY(const OpenMM_MonteCarloAnisotropicBarostat* target) {
    bool result = reinterpret_cast<const OpenMM::MonteCarloAnisotropicBarostat*>(target)->getScaleY();
    return (result ? OpenMM_True : OpenMM_False);
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_MonteCarloAnisotropicBarostat_getScaleZ(const OpenMM_MonteCarloAnisotropicBarostat* target) {
    bool result = reinterpret_cast<const OpenMM::MonteCarloAnisotropicBarostat*>(target)->getScaleZ();
    return (result ? OpenMM_True : OpenMM_False);
}
OPENMM_EXPORT int OpenMM_MonteCarloAnisotropicBarostat_getFrequency(const OpenMM_MonteCarloAnisotropicBarostat* target) {
    int result = reinterpret_cast<const OpenMM::MonteCarloAnisotropicBarostat*>(target)->getFrequency();
    return result;
}
OPENMM_EXPORT void OpenMM_MonteCarloAnisotropicBarostat_setFrequency(OpenMM_MonteCarloAnisotropicBarostat* target, int freq) {
    reinterpret_cast<OpenMM::MonteCarloAnisotropicBarostat*>(target)->setFrequency(freq);
}
OPENMM_EXPORT double OpenMM_MonteCarloAnisotropicBarostat_getDefaultTemperature(const OpenMM_MonteCarloAnisotropicBarostat* target) {
    double result = reinterpret_cast<const OpenMM::MonteCarloAnisotropicBarostat*>(target)->getDefaultTemperature();
    return result;
}
OPENMM_EXPORT void OpenMM_MonteCarloAnisotropicBarostat_setDefaultTemperature(OpenMM_MonteCarloAnisotropicBarostat* target, double temp) {
    reinterpret_cast<OpenMM::MonteCarloAnisotropicBarostat*>(target)->setDefaultTemperature(temp);
}
OPENMM_EXPORT int OpenMM_MonteCarloAnisotropicBarostat_getRandomNumberSeed(const OpenMM_MonteCarloAnisotropicBarostat* target) {
    int result = reinterpret_cast<const OpenMM::MonteCarloAnisotropicBarostat*>(target)->getRandomNumberSeed();
    return result;
}
OPENMM_EXPORT void OpenMM_MonteCarloAnisotropicBarostat_setRandomNumberSeed(OpenMM_MonteCarloAnisotropicBarostat* target, int seed) {
    reinterpret_cast<OpenMM::MonteCarloAnisotropicBarostat*>(target)->setRandomNumberSeed(seed);
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_MonteCarloAnisotropicBarostat_usesPeriodicBoundaryConditions(const OpenMM_MonteCarloAnisotropicBarostat* target) {
    bool result = reinterpret_cast<const OpenMM::MonteCarloAnisotropicBarostat*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::NoseHooverIntegrator */
OPENMM_EXPORT OpenMM_NoseHooverIntegrator* OpenMM_NoseHooverIntegrator_create(double stepSize) {
    return reinterpret_cast<OpenMM_NoseHooverIntegrator*>(new OpenMM::NoseHooverIntegrator(stepSize));
}
OPENMM_EXPORT OpenMM_NoseHooverIntegrator* OpenMM_NoseHooverIntegrator_create_2(double temperature, double collisionFrequency, double stepSize, int chainLength, int numMTS, int numYoshidaSuzuki) {
    return reinterpret_cast<OpenMM_NoseHooverIntegrator*>(new OpenMM::NoseHooverIntegrator(temperature, collisionFrequency, stepSize, chainLength, numMTS, numYoshidaSuzuki));
}
OPENMM_EXPORT void OpenMM_NoseHooverIntegrator_destroy(OpenMM_NoseHooverIntegrator* target) {
    delete reinterpret_cast<OpenMM::NoseHooverIntegrator*>(target);
}
OPENMM_EXPORT void OpenMM_NoseHooverIntegrator_step(OpenMM_NoseHooverIntegrator* target, int steps) {
    reinterpret_cast<OpenMM::NoseHooverIntegrator*>(target)->step(steps);
}
OPENMM_EXPORT int OpenMM_NoseHooverIntegrator_addThermostat(OpenMM_NoseHooverIntegrator* target, double temperature, double collisionFrequency, int chainLength, int numMTS, int numYoshidaSuzuki) {
    int result = reinterpret_cast<OpenMM::NoseHooverIntegrator*>(target)->addThermostat(temperature, collisionFrequency, chainLength, numMTS, numYoshidaSuzuki);
    return result;
}
OPENMM_EXPORT int OpenMM_NoseHooverIntegrator_addSubsystemThermostat(OpenMM_NoseHooverIntegrator* target, const OpenMM_IntArray* thermostatedParticles, const OpenMM_BondArray* thermostatedPairs, double temperature, double collisionFrequency, double relativeTemperature, double relativeCollisionFrequency, int chainLength, int numMTS, int numYoshidaSuzuki) {
    int result = reinterpret_cast<OpenMM::NoseHooverIntegrator*>(target)->addSubsystemThermostat(*reinterpret_cast<const std::vector< int >*>(thermostatedParticles), *reinterpret_cast<const std::vector< std::pair< int, int > >*>(thermostatedPairs), temperature, collisionFrequency, relativeTemperature, relativeCollisionFrequency, chainLength, numMTS, numYoshidaSuzuki);
    return result;
}
OPENMM_EXPORT double OpenMM_NoseHooverIntegrator_getTemperature(const OpenMM_NoseHooverIntegrator* target, int chainID) {
    double result = reinterpret_cast<const OpenMM::NoseHooverIntegrator*>(target)->getTemperature(chainID);
    return result;
}
OPENMM_EXPORT void OpenMM_NoseHooverIntegrator_setTemperature(OpenMM_NoseHooverIntegrator* target, double temperature, int chainID) {
    reinterpret_cast<OpenMM::NoseHooverIntegrator*>(target)->setTemperature(temperature, chainID);
}
OPENMM_EXPORT double OpenMM_NoseHooverIntegrator_getRelativeTemperature(const OpenMM_NoseHooverIntegrator* target, int chainID) {
    double result = reinterpret_cast<const OpenMM::NoseHooverIntegrator*>(target)->getRelativeTemperature(chainID);
    return result;
}
OPENMM_EXPORT void OpenMM_NoseHooverIntegrator_setRelativeTemperature(OpenMM_NoseHooverIntegrator* target, double temperature, int chainID) {
    reinterpret_cast<OpenMM::NoseHooverIntegrator*>(target)->setRelativeTemperature(temperature, chainID);
}
OPENMM_EXPORT double OpenMM_NoseHooverIntegrator_getCollisionFrequency(const OpenMM_NoseHooverIntegrator* target, int chainID) {
    double result = reinterpret_cast<const OpenMM::NoseHooverIntegrator*>(target)->getCollisionFrequency(chainID);
    return result;
}
OPENMM_EXPORT void OpenMM_NoseHooverIntegrator_setCollisionFrequency(OpenMM_NoseHooverIntegrator* target, double frequency, int chainID) {
    reinterpret_cast<OpenMM::NoseHooverIntegrator*>(target)->setCollisionFrequency(frequency, chainID);
}
OPENMM_EXPORT double OpenMM_NoseHooverIntegrator_getRelativeCollisionFrequency(const OpenMM_NoseHooverIntegrator* target, int chainID) {
    double result = reinterpret_cast<const OpenMM::NoseHooverIntegrator*>(target)->getRelativeCollisionFrequency(chainID);
    return result;
}
OPENMM_EXPORT void OpenMM_NoseHooverIntegrator_setRelativeCollisionFrequency(OpenMM_NoseHooverIntegrator* target, double frequency, int chainID) {
    reinterpret_cast<OpenMM::NoseHooverIntegrator*>(target)->setRelativeCollisionFrequency(frequency, chainID);
}
OPENMM_EXPORT double OpenMM_NoseHooverIntegrator_computeHeatBathEnergy(OpenMM_NoseHooverIntegrator* target) {
    double result = reinterpret_cast<OpenMM::NoseHooverIntegrator*>(target)->computeHeatBathEnergy();
    return result;
}
OPENMM_EXPORT int OpenMM_NoseHooverIntegrator_getNumThermostats(const OpenMM_NoseHooverIntegrator* target) {
    int result = reinterpret_cast<const OpenMM::NoseHooverIntegrator*>(target)->getNumThermostats();
    return result;
}
OPENMM_EXPORT const OpenMM_NoseHooverChain* OpenMM_NoseHooverIntegrator_getThermostat(const OpenMM_NoseHooverIntegrator* target, int chainID) {
    const NoseHooverChain* result = &reinterpret_cast<const OpenMM::NoseHooverIntegrator*>(target)->getThermostat(chainID);
    return reinterpret_cast<const OpenMM_NoseHooverChain*>(result);
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_NoseHooverIntegrator_hasSubsystemThermostats(const OpenMM_NoseHooverIntegrator* target) {
    bool result = reinterpret_cast<const OpenMM::NoseHooverIntegrator*>(target)->hasSubsystemThermostats();
    return (result ? OpenMM_True : OpenMM_False);
}
OPENMM_EXPORT double OpenMM_NoseHooverIntegrator_getMaximumPairDistance(const OpenMM_NoseHooverIntegrator* target) {
    double result = reinterpret_cast<const OpenMM::NoseHooverIntegrator*>(target)->getMaximumPairDistance();
    return result;
}
OPENMM_EXPORT void OpenMM_NoseHooverIntegrator_setMaximumPairDistance(OpenMM_NoseHooverIntegrator* target, double distance) {
    reinterpret_cast<OpenMM::NoseHooverIntegrator*>(target)->setMaximumPairDistance(distance);
}

/* OpenMM::VariableVerletIntegrator */
OPENMM_EXPORT OpenMM_VariableVerletIntegrator* OpenMM_VariableVerletIntegrator_create(double errorTol) {
    return reinterpret_cast<OpenMM_VariableVerletIntegrator*>(new OpenMM::VariableVerletIntegrator(errorTol));
}
OPENMM_EXPORT void OpenMM_VariableVerletIntegrator_destroy(OpenMM_VariableVerletIntegrator* target) {
    delete reinterpret_cast<OpenMM::VariableVerletIntegrator*>(target);
}
OPENMM_EXPORT double OpenMM_VariableVerletIntegrator_getErrorTolerance(const OpenMM_VariableVerletIntegrator* target) {
    double result = reinterpret_cast<const OpenMM::VariableVerletIntegrator*>(target)->getErrorTolerance();
    return result;
}
OPENMM_EXPORT void OpenMM_VariableVerletIntegrator_setErrorTolerance(OpenMM_VariableVerletIntegrator* target, double tol) {
    reinterpret_cast<OpenMM::VariableVerletIntegrator*>(target)->setErrorTolerance(tol);
}
OPENMM_EXPORT double OpenMM_VariableVerletIntegrator_getMaximumStepSize(const OpenMM_VariableVerletIntegrator* target) {
    double result = reinterpret_cast<const OpenMM::VariableVerletIntegrator*>(target)->getMaximumStepSize();
    return result;
}
OPENMM_EXPORT void OpenMM_VariableVerletIntegrator_setMaximumStepSize(OpenMM_VariableVerletIntegrator* target, double size) {
    reinterpret_cast<OpenMM::VariableVerletIntegrator*>(target)->setMaximumStepSize(size);
}
OPENMM_EXPORT void OpenMM_VariableVerletIntegrator_step(OpenMM_VariableVerletIntegrator* target, int steps) {
    reinterpret_cast<OpenMM::VariableVerletIntegrator*>(target)->step(steps);
}
OPENMM_EXPORT void OpenMM_VariableVerletIntegrator_stepTo(OpenMM_VariableVerletIntegrator* target, double time) {
    reinterpret_cast<OpenMM::VariableVerletIntegrator*>(target)->stepTo(time);
}

/* OpenMM::CustomHbondForce */
OPENMM_EXPORT OpenMM_CustomHbondForce* OpenMM_CustomHbondForce_create(const char* energy) {
    return reinterpret_cast<OpenMM_CustomHbondForce*>(new OpenMM::CustomHbondForce(std::string(energy)));
}
OPENMM_EXPORT void OpenMM_CustomHbondForce_destroy(OpenMM_CustomHbondForce* target) {
    delete reinterpret_cast<OpenMM::CustomHbondForce*>(target);
}
OPENMM_EXPORT int OpenMM_CustomHbondForce_getNumDonors(const OpenMM_CustomHbondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomHbondForce*>(target)->getNumDonors();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomHbondForce_getNumAcceptors(const OpenMM_CustomHbondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomHbondForce*>(target)->getNumAcceptors();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomHbondForce_getNumExclusions(const OpenMM_CustomHbondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomHbondForce*>(target)->getNumExclusions();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomHbondForce_getNumPerDonorParameters(const OpenMM_CustomHbondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomHbondForce*>(target)->getNumPerDonorParameters();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomHbondForce_getNumPerAcceptorParameters(const OpenMM_CustomHbondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomHbondForce*>(target)->getNumPerAcceptorParameters();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomHbondForce_getNumGlobalParameters(const OpenMM_CustomHbondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomHbondForce*>(target)->getNumGlobalParameters();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomHbondForce_getNumTabulatedFunctions(const OpenMM_CustomHbondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomHbondForce*>(target)->getNumTabulatedFunctions();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomHbondForce_getNumFunctions(const OpenMM_CustomHbondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomHbondForce*>(target)->getNumFunctions();
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomHbondForce_getEnergyFunction(const OpenMM_CustomHbondForce* target) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomHbondForce*>(target)->getEnergyFunction();
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomHbondForce_setEnergyFunction(OpenMM_CustomHbondForce* target, const char* energy) {
    reinterpret_cast<OpenMM::CustomHbondForce*>(target)->setEnergyFunction(std::string(energy));
}
OPENMM_EXPORT OpenMM_CustomHbondForce_NonbondedMethod OpenMM_CustomHbondForce_getNonbondedMethod(const OpenMM_CustomHbondForce* target) {
    OpenMM::CustomHbondForce::NonbondedMethod result = reinterpret_cast<const OpenMM::CustomHbondForce*>(target)->getNonbondedMethod();
    return static_cast<OpenMM_CustomHbondForce_NonbondedMethod>(result);
}
OPENMM_EXPORT void OpenMM_CustomHbondForce_setNonbondedMethod(OpenMM_CustomHbondForce* target, OpenMM_CustomHbondForce_NonbondedMethod method) {
    reinterpret_cast<OpenMM::CustomHbondForce*>(target)->setNonbondedMethod(static_cast<OpenMM::CustomHbondForce::NonbondedMethod>(method));
}
OPENMM_EXPORT double OpenMM_CustomHbondForce_getCutoffDistance(const OpenMM_CustomHbondForce* target) {
    double result = reinterpret_cast<const OpenMM::CustomHbondForce*>(target)->getCutoffDistance();
    return result;
}
OPENMM_EXPORT void OpenMM_CustomHbondForce_setCutoffDistance(OpenMM_CustomHbondForce* target, double distance) {
    reinterpret_cast<OpenMM::CustomHbondForce*>(target)->setCutoffDistance(distance);
}
OPENMM_EXPORT int OpenMM_CustomHbondForce_addPerDonorParameter(OpenMM_CustomHbondForce* target, const char* name) {
    int result = reinterpret_cast<OpenMM::CustomHbondForce*>(target)->addPerDonorParameter(std::string(name));
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomHbondForce_getPerDonorParameterName(const OpenMM_CustomHbondForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomHbondForce*>(target)->getPerDonorParameterName(index);
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomHbondForce_setPerDonorParameterName(OpenMM_CustomHbondForce* target, int index, const char* name) {
    reinterpret_cast<OpenMM::CustomHbondForce*>(target)->setPerDonorParameterName(index, std::string(name));
}
OPENMM_EXPORT int OpenMM_CustomHbondForce_addPerAcceptorParameter(OpenMM_CustomHbondForce* target, const char* name) {
    int result = reinterpret_cast<OpenMM::CustomHbondForce*>(target)->addPerAcceptorParameter(std::string(name));
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomHbondForce_getPerAcceptorParameterName(const OpenMM_CustomHbondForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomHbondForce*>(target)->getPerAcceptorParameterName(index);
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomHbondForce_setPerAcceptorParameterName(OpenMM_CustomHbondForce* target, int index, const char* name) {
    reinterpret_cast<OpenMM::CustomHbondForce*>(target)->setPerAcceptorParameterName(index, std::string(name));
}
OPENMM_EXPORT int OpenMM_CustomHbondForce_addGlobalParameter(OpenMM_CustomHbondForce* target, const char* name, double defaultValue) {
    int result = reinterpret_cast<OpenMM::CustomHbondForce*>(target)->addGlobalParameter(std::string(name), defaultValue);
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomHbondForce_getGlobalParameterName(const OpenMM_CustomHbondForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomHbondForce*>(target)->getGlobalParameterName(index);
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomHbondForce_setGlobalParameterName(OpenMM_CustomHbondForce* target, int index, const char* name) {
    reinterpret_cast<OpenMM::CustomHbondForce*>(target)->setGlobalParameterName(index, std::string(name));
}
OPENMM_EXPORT double OpenMM_CustomHbondForce_getGlobalParameterDefaultValue(const OpenMM_CustomHbondForce* target, int index) {
    double result = reinterpret_cast<const OpenMM::CustomHbondForce*>(target)->getGlobalParameterDefaultValue(index);
    return result;
}
OPENMM_EXPORT void OpenMM_CustomHbondForce_setGlobalParameterDefaultValue(OpenMM_CustomHbondForce* target, int index, double defaultValue) {
    reinterpret_cast<OpenMM::CustomHbondForce*>(target)->setGlobalParameterDefaultValue(index, defaultValue);
}
OPENMM_EXPORT int OpenMM_CustomHbondForce_addDonor(OpenMM_CustomHbondForce* target, int d1, int d2, int d3, const OpenMM_DoubleArray* parameters) {
    int result = reinterpret_cast<OpenMM::CustomHbondForce*>(target)->addDonor(d1, d2, d3, *reinterpret_cast<const std::vector< double >*>(parameters));
    return result;
}
OPENMM_EXPORT void OpenMM_CustomHbondForce_getDonorParameters(const OpenMM_CustomHbondForce* target, int index, int* d1, int* d2, int* d3, OpenMM_DoubleArray* parameters) {
    reinterpret_cast<const OpenMM::CustomHbondForce*>(target)->getDonorParameters(index, *reinterpret_cast<int*>(d1), *reinterpret_cast<int*>(d2), *reinterpret_cast<int*>(d3), *reinterpret_cast<std::vector< double >*>(parameters));
}
OPENMM_EXPORT void OpenMM_CustomHbondForce_setDonorParameters(OpenMM_CustomHbondForce* target, int index, int d1, int d2, int d3, const OpenMM_DoubleArray* parameters) {
    reinterpret_cast<OpenMM::CustomHbondForce*>(target)->setDonorParameters(index, d1, d2, d3, *reinterpret_cast<const std::vector< double >*>(parameters));
}
OPENMM_EXPORT int OpenMM_CustomHbondForce_addAcceptor(OpenMM_CustomHbondForce* target, int a1, int a2, int a3, const OpenMM_DoubleArray* parameters) {
    int result = reinterpret_cast<OpenMM::CustomHbondForce*>(target)->addAcceptor(a1, a2, a3, *reinterpret_cast<const std::vector< double >*>(parameters));
    return result;
}
OPENMM_EXPORT void OpenMM_CustomHbondForce_getAcceptorParameters(const OpenMM_CustomHbondForce* target, int index, int* a1, int* a2, int* a3, OpenMM_DoubleArray* parameters) {
    reinterpret_cast<const OpenMM::CustomHbondForce*>(target)->getAcceptorParameters(index, *reinterpret_cast<int*>(a1), *reinterpret_cast<int*>(a2), *reinterpret_cast<int*>(a3), *reinterpret_cast<std::vector< double >*>(parameters));
}
OPENMM_EXPORT void OpenMM_CustomHbondForce_setAcceptorParameters(OpenMM_CustomHbondForce* target, int index, int a1, int a2, int a3, const OpenMM_DoubleArray* parameters) {
    reinterpret_cast<OpenMM::CustomHbondForce*>(target)->setAcceptorParameters(index, a1, a2, a3, *reinterpret_cast<const std::vector< double >*>(parameters));
}
OPENMM_EXPORT int OpenMM_CustomHbondForce_addExclusion(OpenMM_CustomHbondForce* target, int donor, int acceptor) {
    int result = reinterpret_cast<OpenMM::CustomHbondForce*>(target)->addExclusion(donor, acceptor);
    return result;
}
OPENMM_EXPORT void OpenMM_CustomHbondForce_getExclusionParticles(const OpenMM_CustomHbondForce* target, int index, int* donor, int* acceptor) {
    reinterpret_cast<const OpenMM::CustomHbondForce*>(target)->getExclusionParticles(index, *reinterpret_cast<int*>(donor), *reinterpret_cast<int*>(acceptor));
}
OPENMM_EXPORT void OpenMM_CustomHbondForce_setExclusionParticles(OpenMM_CustomHbondForce* target, int index, int donor, int acceptor) {
    reinterpret_cast<OpenMM::CustomHbondForce*>(target)->setExclusionParticles(index, donor, acceptor);
}
OPENMM_EXPORT int OpenMM_CustomHbondForce_addTabulatedFunction(OpenMM_CustomHbondForce* target, const char* name, OpenMM_TabulatedFunction* function) {
    int result = reinterpret_cast<OpenMM::CustomHbondForce*>(target)->addTabulatedFunction(std::string(name), reinterpret_cast<TabulatedFunction *>(function));
    return result;
}
OPENMM_EXPORT OpenMM_TabulatedFunction* OpenMM_CustomHbondForce_getTabulatedFunction(OpenMM_CustomHbondForce* target, int index) {
    TabulatedFunction* result = &reinterpret_cast<OpenMM::CustomHbondForce*>(target)->getTabulatedFunction(index);
    return reinterpret_cast<OpenMM_TabulatedFunction*>(result);
}
OPENMM_EXPORT const char* OpenMM_CustomHbondForce_getTabulatedFunctionName(const OpenMM_CustomHbondForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomHbondForce*>(target)->getTabulatedFunctionName(index);
    return result->c_str();
}
OPENMM_EXPORT int OpenMM_CustomHbondForce_addFunction(OpenMM_CustomHbondForce* target, const char* name, const OpenMM_DoubleArray* values, double min, double max) {
    int result = reinterpret_cast<OpenMM::CustomHbondForce*>(target)->addFunction(std::string(name), *reinterpret_cast<const std::vector< double >*>(values), min, max);
    return result;
}
OPENMM_EXPORT void OpenMM_CustomHbondForce_getFunctionParameters(const OpenMM_CustomHbondForce* target, int index, char** name, OpenMM_DoubleArray* values, double* min, double* max) {
    reinterpret_cast<const OpenMM::CustomHbondForce*>(target)->getFunctionParameters(index, *reinterpret_cast<std::string*>(name), *reinterpret_cast<std::vector< double >*>(values), *reinterpret_cast<double*>(min), *reinterpret_cast<double*>(max));
}
OPENMM_EXPORT void OpenMM_CustomHbondForce_setFunctionParameters(OpenMM_CustomHbondForce* target, int index, const char* name, const OpenMM_DoubleArray* values, double min, double max) {
    reinterpret_cast<OpenMM::CustomHbondForce*>(target)->setFunctionParameters(index, std::string(name), *reinterpret_cast<const std::vector< double >*>(values), min, max);
}
OPENMM_EXPORT void OpenMM_CustomHbondForce_updateParametersInContext(OpenMM_CustomHbondForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::CustomHbondForce*>(target)->updateParametersInContext(*reinterpret_cast<OpenMM::Context*>(context));
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_CustomHbondForce_usesPeriodicBoundaryConditions(const OpenMM_CustomHbondForce* target) {
    bool result = reinterpret_cast<const OpenMM::CustomHbondForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::CustomExternalForce */
OPENMM_EXPORT OpenMM_CustomExternalForce* OpenMM_CustomExternalForce_create(const char* energy) {
    return reinterpret_cast<OpenMM_CustomExternalForce*>(new OpenMM::CustomExternalForce(std::string(energy)));
}
OPENMM_EXPORT void OpenMM_CustomExternalForce_destroy(OpenMM_CustomExternalForce* target) {
    delete reinterpret_cast<OpenMM::CustomExternalForce*>(target);
}
OPENMM_EXPORT int OpenMM_CustomExternalForce_getNumParticles(const OpenMM_CustomExternalForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomExternalForce*>(target)->getNumParticles();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomExternalForce_getNumPerParticleParameters(const OpenMM_CustomExternalForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomExternalForce*>(target)->getNumPerParticleParameters();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomExternalForce_getNumGlobalParameters(const OpenMM_CustomExternalForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomExternalForce*>(target)->getNumGlobalParameters();
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomExternalForce_getEnergyFunction(const OpenMM_CustomExternalForce* target) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomExternalForce*>(target)->getEnergyFunction();
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomExternalForce_setEnergyFunction(OpenMM_CustomExternalForce* target, const char* energy) {
    reinterpret_cast<OpenMM::CustomExternalForce*>(target)->setEnergyFunction(std::string(energy));
}
OPENMM_EXPORT int OpenMM_CustomExternalForce_addPerParticleParameter(OpenMM_CustomExternalForce* target, const char* name) {
    int result = reinterpret_cast<OpenMM::CustomExternalForce*>(target)->addPerParticleParameter(std::string(name));
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomExternalForce_getPerParticleParameterName(const OpenMM_CustomExternalForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomExternalForce*>(target)->getPerParticleParameterName(index);
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomExternalForce_setPerParticleParameterName(OpenMM_CustomExternalForce* target, int index, const char* name) {
    reinterpret_cast<OpenMM::CustomExternalForce*>(target)->setPerParticleParameterName(index, std::string(name));
}
OPENMM_EXPORT int OpenMM_CustomExternalForce_addGlobalParameter(OpenMM_CustomExternalForce* target, const char* name, double defaultValue) {
    int result = reinterpret_cast<OpenMM::CustomExternalForce*>(target)->addGlobalParameter(std::string(name), defaultValue);
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomExternalForce_getGlobalParameterName(const OpenMM_CustomExternalForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomExternalForce*>(target)->getGlobalParameterName(index);
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomExternalForce_setGlobalParameterName(OpenMM_CustomExternalForce* target, int index, const char* name) {
    reinterpret_cast<OpenMM::CustomExternalForce*>(target)->setGlobalParameterName(index, std::string(name));
}
OPENMM_EXPORT double OpenMM_CustomExternalForce_getGlobalParameterDefaultValue(const OpenMM_CustomExternalForce* target, int index) {
    double result = reinterpret_cast<const OpenMM::CustomExternalForce*>(target)->getGlobalParameterDefaultValue(index);
    return result;
}
OPENMM_EXPORT void OpenMM_CustomExternalForce_setGlobalParameterDefaultValue(OpenMM_CustomExternalForce* target, int index, double defaultValue) {
    reinterpret_cast<OpenMM::CustomExternalForce*>(target)->setGlobalParameterDefaultValue(index, defaultValue);
}
OPENMM_EXPORT int OpenMM_CustomExternalForce_addParticle(OpenMM_CustomExternalForce* target, int particle, const OpenMM_DoubleArray* parameters) {
    int result = reinterpret_cast<OpenMM::CustomExternalForce*>(target)->addParticle(particle, *reinterpret_cast<const std::vector< double >*>(parameters));
    return result;
}
OPENMM_EXPORT void OpenMM_CustomExternalForce_getParticleParameters(const OpenMM_CustomExternalForce* target, int index, int* particle, OpenMM_DoubleArray* parameters) {
    reinterpret_cast<const OpenMM::CustomExternalForce*>(target)->getParticleParameters(index, *reinterpret_cast<int*>(particle), *reinterpret_cast<std::vector< double >*>(parameters));
}
OPENMM_EXPORT void OpenMM_CustomExternalForce_setParticleParameters(OpenMM_CustomExternalForce* target, int index, int particle, const OpenMM_DoubleArray* parameters) {
    reinterpret_cast<OpenMM::CustomExternalForce*>(target)->setParticleParameters(index, particle, *reinterpret_cast<const std::vector< double >*>(parameters));
}
OPENMM_EXPORT void OpenMM_CustomExternalForce_updateParametersInContext(OpenMM_CustomExternalForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::CustomExternalForce*>(target)->updateParametersInContext(*reinterpret_cast<OpenMM::Context*>(context));
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_CustomExternalForce_usesPeriodicBoundaryConditions(const OpenMM_CustomExternalForce* target) {
    bool result = reinterpret_cast<const OpenMM::CustomExternalForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::CustomIntegrator */
OPENMM_EXPORT OpenMM_CustomIntegrator* OpenMM_CustomIntegrator_create(double stepSize) {
    return reinterpret_cast<OpenMM_CustomIntegrator*>(new OpenMM::CustomIntegrator(stepSize));
}
OPENMM_EXPORT void OpenMM_CustomIntegrator_destroy(OpenMM_CustomIntegrator* target) {
    delete reinterpret_cast<OpenMM::CustomIntegrator*>(target);
}
OPENMM_EXPORT int OpenMM_CustomIntegrator_getNumGlobalVariables(const OpenMM_CustomIntegrator* target) {
    int result = reinterpret_cast<const OpenMM::CustomIntegrator*>(target)->getNumGlobalVariables();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomIntegrator_getNumPerDofVariables(const OpenMM_CustomIntegrator* target) {
    int result = reinterpret_cast<const OpenMM::CustomIntegrator*>(target)->getNumPerDofVariables();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomIntegrator_getNumComputations(const OpenMM_CustomIntegrator* target) {
    int result = reinterpret_cast<const OpenMM::CustomIntegrator*>(target)->getNumComputations();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomIntegrator_getNumTabulatedFunctions(const OpenMM_CustomIntegrator* target) {
    int result = reinterpret_cast<const OpenMM::CustomIntegrator*>(target)->getNumTabulatedFunctions();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomIntegrator_addGlobalVariable(OpenMM_CustomIntegrator* target, const char* name, double initialValue) {
    int result = reinterpret_cast<OpenMM::CustomIntegrator*>(target)->addGlobalVariable(std::string(name), initialValue);
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomIntegrator_getGlobalVariableName(const OpenMM_CustomIntegrator* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomIntegrator*>(target)->getGlobalVariableName(index);
    return result->c_str();
}
OPENMM_EXPORT int OpenMM_CustomIntegrator_addPerDofVariable(OpenMM_CustomIntegrator* target, const char* name, double initialValue) {
    int result = reinterpret_cast<OpenMM::CustomIntegrator*>(target)->addPerDofVariable(std::string(name), initialValue);
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomIntegrator_getPerDofVariableName(const OpenMM_CustomIntegrator* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomIntegrator*>(target)->getPerDofVariableName(index);
    return result->c_str();
}
OPENMM_EXPORT double OpenMM_CustomIntegrator_getGlobalVariable(const OpenMM_CustomIntegrator* target, int index) {
    double result = reinterpret_cast<const OpenMM::CustomIntegrator*>(target)->getGlobalVariable(index);
    return result;
}
OPENMM_EXPORT double OpenMM_CustomIntegrator_getGlobalVariableByName(const OpenMM_CustomIntegrator* target, const char* name) {
    double result = reinterpret_cast<const OpenMM::CustomIntegrator*>(target)->getGlobalVariableByName(std::string(name));
    return result;
}
OPENMM_EXPORT void OpenMM_CustomIntegrator_setGlobalVariable(OpenMM_CustomIntegrator* target, int index, double value) {
    reinterpret_cast<OpenMM::CustomIntegrator*>(target)->setGlobalVariable(index, value);
}
OPENMM_EXPORT void OpenMM_CustomIntegrator_setGlobalVariableByName(OpenMM_CustomIntegrator* target, const char* name, double value) {
    reinterpret_cast<OpenMM::CustomIntegrator*>(target)->setGlobalVariableByName(std::string(name), value);
}
OPENMM_EXPORT void OpenMM_CustomIntegrator_getPerDofVariable(const OpenMM_CustomIntegrator* target, int index, OpenMM_Vec3Array* values) {
    reinterpret_cast<const OpenMM::CustomIntegrator*>(target)->getPerDofVariable(index, *reinterpret_cast<std::vector< Vec3 >*>(values));
}
OPENMM_EXPORT void OpenMM_CustomIntegrator_getPerDofVariableByName(const OpenMM_CustomIntegrator* target, const char* name, OpenMM_Vec3Array* values) {
    reinterpret_cast<const OpenMM::CustomIntegrator*>(target)->getPerDofVariableByName(std::string(name), *reinterpret_cast<std::vector< Vec3 >*>(values));
}
OPENMM_EXPORT void OpenMM_CustomIntegrator_setPerDofVariable(OpenMM_CustomIntegrator* target, int index, const OpenMM_Vec3Array* values) {
    reinterpret_cast<OpenMM::CustomIntegrator*>(target)->setPerDofVariable(index, *reinterpret_cast<const std::vector< Vec3 >*>(values));
}
OPENMM_EXPORT void OpenMM_CustomIntegrator_setPerDofVariableByName(OpenMM_CustomIntegrator* target, const char* name, const OpenMM_Vec3Array* values) {
    reinterpret_cast<OpenMM::CustomIntegrator*>(target)->setPerDofVariableByName(std::string(name), *reinterpret_cast<const std::vector< Vec3 >*>(values));
}
OPENMM_EXPORT int OpenMM_CustomIntegrator_addComputeGlobal(OpenMM_CustomIntegrator* target, const char* variable, const char* expression) {
    int result = reinterpret_cast<OpenMM::CustomIntegrator*>(target)->addComputeGlobal(std::string(variable), std::string(expression));
    return result;
}
OPENMM_EXPORT int OpenMM_CustomIntegrator_addComputePerDof(OpenMM_CustomIntegrator* target, const char* variable, const char* expression) {
    int result = reinterpret_cast<OpenMM::CustomIntegrator*>(target)->addComputePerDof(std::string(variable), std::string(expression));
    return result;
}
OPENMM_EXPORT int OpenMM_CustomIntegrator_addComputeSum(OpenMM_CustomIntegrator* target, const char* variable, const char* expression) {
    int result = reinterpret_cast<OpenMM::CustomIntegrator*>(target)->addComputeSum(std::string(variable), std::string(expression));
    return result;
}
OPENMM_EXPORT int OpenMM_CustomIntegrator_addConstrainPositions(OpenMM_CustomIntegrator* target) {
    int result = reinterpret_cast<OpenMM::CustomIntegrator*>(target)->addConstrainPositions();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomIntegrator_addConstrainVelocities(OpenMM_CustomIntegrator* target) {
    int result = reinterpret_cast<OpenMM::CustomIntegrator*>(target)->addConstrainVelocities();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomIntegrator_addUpdateContextState(OpenMM_CustomIntegrator* target) {
    int result = reinterpret_cast<OpenMM::CustomIntegrator*>(target)->addUpdateContextState();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomIntegrator_beginIfBlock(OpenMM_CustomIntegrator* target, const char* condition) {
    int result = reinterpret_cast<OpenMM::CustomIntegrator*>(target)->beginIfBlock(std::string(condition));
    return result;
}
OPENMM_EXPORT int OpenMM_CustomIntegrator_beginWhileBlock(OpenMM_CustomIntegrator* target, const char* condition) {
    int result = reinterpret_cast<OpenMM::CustomIntegrator*>(target)->beginWhileBlock(std::string(condition));
    return result;
}
OPENMM_EXPORT int OpenMM_CustomIntegrator_endBlock(OpenMM_CustomIntegrator* target) {
    int result = reinterpret_cast<OpenMM::CustomIntegrator*>(target)->endBlock();
    return result;
}
OPENMM_EXPORT void OpenMM_CustomIntegrator_getComputationStep(const OpenMM_CustomIntegrator* target, int index, OpenMM_CustomIntegrator_ComputationType* type, char** variable, char** expression) {
    reinterpret_cast<const OpenMM::CustomIntegrator*>(target)->getComputationStep(index, *reinterpret_cast<OpenMM::CustomIntegrator::ComputationType*>(type), *reinterpret_cast<std::string*>(variable), *reinterpret_cast<std::string*>(expression));
}
OPENMM_EXPORT int OpenMM_CustomIntegrator_addTabulatedFunction(OpenMM_CustomIntegrator* target, const char* name, OpenMM_TabulatedFunction* function) {
    int result = reinterpret_cast<OpenMM::CustomIntegrator*>(target)->addTabulatedFunction(std::string(name), reinterpret_cast<TabulatedFunction *>(function));
    return result;
}
OPENMM_EXPORT OpenMM_TabulatedFunction* OpenMM_CustomIntegrator_getTabulatedFunction(OpenMM_CustomIntegrator* target, int index) {
    TabulatedFunction* result = &reinterpret_cast<OpenMM::CustomIntegrator*>(target)->getTabulatedFunction(index);
    return reinterpret_cast<OpenMM_TabulatedFunction*>(result);
}
OPENMM_EXPORT const char* OpenMM_CustomIntegrator_getTabulatedFunctionName(const OpenMM_CustomIntegrator* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomIntegrator*>(target)->getTabulatedFunctionName(index);
    return result->c_str();
}
OPENMM_EXPORT const char* OpenMM_CustomIntegrator_getKineticEnergyExpression(const OpenMM_CustomIntegrator* target) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomIntegrator*>(target)->getKineticEnergyExpression();
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomIntegrator_setKineticEnergyExpression(OpenMM_CustomIntegrator* target, const char* expression) {
    reinterpret_cast<OpenMM::CustomIntegrator*>(target)->setKineticEnergyExpression(std::string(expression));
}
OPENMM_EXPORT int OpenMM_CustomIntegrator_getRandomNumberSeed(const OpenMM_CustomIntegrator* target) {
    int result = reinterpret_cast<const OpenMM::CustomIntegrator*>(target)->getRandomNumberSeed();
    return result;
}
OPENMM_EXPORT void OpenMM_CustomIntegrator_setRandomNumberSeed(OpenMM_CustomIntegrator* target, int seed) {
    reinterpret_cast<OpenMM::CustomIntegrator*>(target)->setRandomNumberSeed(seed);
}
OPENMM_EXPORT void OpenMM_CustomIntegrator_step(OpenMM_CustomIntegrator* target, int steps) {
    reinterpret_cast<OpenMM::CustomIntegrator*>(target)->step(steps);
}

/* OpenMM::Continuous2DFunction */
OPENMM_EXPORT OpenMM_Continuous2DFunction* OpenMM_Continuous2DFunction_create(int xsize, int ysize, const OpenMM_DoubleArray* values, double xmin, double xmax, double ymin, double ymax) {
    return reinterpret_cast<OpenMM_Continuous2DFunction*>(new OpenMM::Continuous2DFunction(xsize, ysize, *reinterpret_cast<const std::vector< double >*>(values), xmin, xmax, ymin, ymax));
}
OPENMM_EXPORT void OpenMM_Continuous2DFunction_destroy(OpenMM_Continuous2DFunction* target) {
    delete reinterpret_cast<OpenMM::Continuous2DFunction*>(target);
}
OPENMM_EXPORT void OpenMM_Continuous2DFunction_getFunctionParameters(const OpenMM_Continuous2DFunction* target, int* xsize, int* ysize, OpenMM_DoubleArray* values, double* xmin, double* xmax, double* ymin, double* ymax) {
    reinterpret_cast<const OpenMM::Continuous2DFunction*>(target)->getFunctionParameters(*reinterpret_cast<int*>(xsize), *reinterpret_cast<int*>(ysize), *reinterpret_cast<std::vector< double >*>(values), *reinterpret_cast<double*>(xmin), *reinterpret_cast<double*>(xmax), *reinterpret_cast<double*>(ymin), *reinterpret_cast<double*>(ymax));
}
OPENMM_EXPORT void OpenMM_Continuous2DFunction_setFunctionParameters(OpenMM_Continuous2DFunction* target, int xsize, int ysize, const OpenMM_DoubleArray* values, double xmin, double xmax, double ymin, double ymax) {
    reinterpret_cast<OpenMM::Continuous2DFunction*>(target)->setFunctionParameters(xsize, ysize, *reinterpret_cast<const std::vector< double >*>(values), xmin, xmax, ymin, ymax);
}
OPENMM_EXPORT OpenMM_Continuous2DFunction* OpenMM_Continuous2DFunction_Copy(const OpenMM_Continuous2DFunction* target) {
    Continuous2DFunction * result = reinterpret_cast<const OpenMM::Continuous2DFunction*>(target)->Copy();
    return reinterpret_cast<OpenMM_Continuous2DFunction*>(result);
}

/* OpenMM::GayBerneForce */
OPENMM_EXPORT OpenMM_GayBerneForce* OpenMM_GayBerneForce_create() {
    return reinterpret_cast<OpenMM_GayBerneForce*>(new OpenMM::GayBerneForce());
}
OPENMM_EXPORT void OpenMM_GayBerneForce_destroy(OpenMM_GayBerneForce* target) {
    delete reinterpret_cast<OpenMM::GayBerneForce*>(target);
}
OPENMM_EXPORT int OpenMM_GayBerneForce_getNumParticles(const OpenMM_GayBerneForce* target) {
    int result = reinterpret_cast<const OpenMM::GayBerneForce*>(target)->getNumParticles();
    return result;
}
OPENMM_EXPORT int OpenMM_GayBerneForce_getNumExceptions(const OpenMM_GayBerneForce* target) {
    int result = reinterpret_cast<const OpenMM::GayBerneForce*>(target)->getNumExceptions();
    return result;
}
OPENMM_EXPORT OpenMM_GayBerneForce_NonbondedMethod OpenMM_GayBerneForce_getNonbondedMethod(const OpenMM_GayBerneForce* target) {
    OpenMM::GayBerneForce::NonbondedMethod result = reinterpret_cast<const OpenMM::GayBerneForce*>(target)->getNonbondedMethod();
    return static_cast<OpenMM_GayBerneForce_NonbondedMethod>(result);
}
OPENMM_EXPORT void OpenMM_GayBerneForce_setNonbondedMethod(OpenMM_GayBerneForce* target, OpenMM_GayBerneForce_NonbondedMethod method) {
    reinterpret_cast<OpenMM::GayBerneForce*>(target)->setNonbondedMethod(static_cast<OpenMM::GayBerneForce::NonbondedMethod>(method));
}
OPENMM_EXPORT double OpenMM_GayBerneForce_getCutoffDistance(const OpenMM_GayBerneForce* target) {
    double result = reinterpret_cast<const OpenMM::GayBerneForce*>(target)->getCutoffDistance();
    return result;
}
OPENMM_EXPORT void OpenMM_GayBerneForce_setCutoffDistance(OpenMM_GayBerneForce* target, double distance) {
    reinterpret_cast<OpenMM::GayBerneForce*>(target)->setCutoffDistance(distance);
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_GayBerneForce_getUseSwitchingFunction(const OpenMM_GayBerneForce* target) {
    bool result = reinterpret_cast<const OpenMM::GayBerneForce*>(target)->getUseSwitchingFunction();
    return (result ? OpenMM_True : OpenMM_False);
}
OPENMM_EXPORT void OpenMM_GayBerneForce_setUseSwitchingFunction(OpenMM_GayBerneForce* target, OpenMM_Boolean use) {
    reinterpret_cast<OpenMM::GayBerneForce*>(target)->setUseSwitchingFunction(use);
}
OPENMM_EXPORT double OpenMM_GayBerneForce_getSwitchingDistance(const OpenMM_GayBerneForce* target) {
    double result = reinterpret_cast<const OpenMM::GayBerneForce*>(target)->getSwitchingDistance();
    return result;
}
OPENMM_EXPORT void OpenMM_GayBerneForce_setSwitchingDistance(OpenMM_GayBerneForce* target, double distance) {
    reinterpret_cast<OpenMM::GayBerneForce*>(target)->setSwitchingDistance(distance);
}
OPENMM_EXPORT int OpenMM_GayBerneForce_addParticle(OpenMM_GayBerneForce* target, double sigma, double epsilon, int xparticle, int yparticle, double sx, double sy, double sz, double ex, double ey, double ez) {
    int result = reinterpret_cast<OpenMM::GayBerneForce*>(target)->addParticle(sigma, epsilon, xparticle, yparticle, sx, sy, sz, ex, ey, ez);
    return result;
}
OPENMM_EXPORT void OpenMM_GayBerneForce_getParticleParameters(const OpenMM_GayBerneForce* target, int index, double* sigma, double* epsilon, int* xparticle, int* yparticle, double* sx, double* sy, double* sz, double* ex, double* ey, double* ez) {
    reinterpret_cast<const OpenMM::GayBerneForce*>(target)->getParticleParameters(index, *reinterpret_cast<double*>(sigma), *reinterpret_cast<double*>(epsilon), *reinterpret_cast<int*>(xparticle), *reinterpret_cast<int*>(yparticle), *reinterpret_cast<double*>(sx), *reinterpret_cast<double*>(sy), *reinterpret_cast<double*>(sz), *reinterpret_cast<double*>(ex), *reinterpret_cast<double*>(ey), *reinterpret_cast<double*>(ez));
}
OPENMM_EXPORT void OpenMM_GayBerneForce_setParticleParameters(OpenMM_GayBerneForce* target, int index, double sigma, double epsilon, int xparticle, int yparticle, double sx, double sy, double sz, double ex, double ey, double ez) {
    reinterpret_cast<OpenMM::GayBerneForce*>(target)->setParticleParameters(index, sigma, epsilon, xparticle, yparticle, sx, sy, sz, ex, ey, ez);
}
OPENMM_EXPORT int OpenMM_GayBerneForce_addException(OpenMM_GayBerneForce* target, int particle1, int particle2, double sigma, double epsilon, OpenMM_Boolean replace) {
    int result = reinterpret_cast<OpenMM::GayBerneForce*>(target)->addException(particle1, particle2, sigma, epsilon, replace);
    return result;
}
OPENMM_EXPORT void OpenMM_GayBerneForce_getExceptionParameters(const OpenMM_GayBerneForce* target, int index, int* particle1, int* particle2, double* sigma, double* epsilon) {
    reinterpret_cast<const OpenMM::GayBerneForce*>(target)->getExceptionParameters(index, *reinterpret_cast<int*>(particle1), *reinterpret_cast<int*>(particle2), *reinterpret_cast<double*>(sigma), *reinterpret_cast<double*>(epsilon));
}
OPENMM_EXPORT void OpenMM_GayBerneForce_setExceptionParameters(OpenMM_GayBerneForce* target, int index, int particle1, int particle2, double sigma, double epsilon) {
    reinterpret_cast<OpenMM::GayBerneForce*>(target)->setExceptionParameters(index, particle1, particle2, sigma, epsilon);
}
OPENMM_EXPORT void OpenMM_GayBerneForce_updateParametersInContext(OpenMM_GayBerneForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::GayBerneForce*>(target)->updateParametersInContext(*reinterpret_cast<OpenMM::Context*>(context));
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_GayBerneForce_usesPeriodicBoundaryConditions(const OpenMM_GayBerneForce* target) {
    bool result = reinterpret_cast<const OpenMM::GayBerneForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::NonbondedForce */
OPENMM_EXPORT OpenMM_NonbondedForce* OpenMM_NonbondedForce_create() {
    return reinterpret_cast<OpenMM_NonbondedForce*>(new OpenMM::NonbondedForce());
}
OPENMM_EXPORT void OpenMM_NonbondedForce_destroy(OpenMM_NonbondedForce* target) {
    delete reinterpret_cast<OpenMM::NonbondedForce*>(target);
}
OPENMM_EXPORT int OpenMM_NonbondedForce_getNumParticles(const OpenMM_NonbondedForce* target) {
    int result = reinterpret_cast<const OpenMM::NonbondedForce*>(target)->getNumParticles();
    return result;
}
OPENMM_EXPORT int OpenMM_NonbondedForce_getNumExceptions(const OpenMM_NonbondedForce* target) {
    int result = reinterpret_cast<const OpenMM::NonbondedForce*>(target)->getNumExceptions();
    return result;
}
OPENMM_EXPORT int OpenMM_NonbondedForce_getNumGlobalParameters(const OpenMM_NonbondedForce* target) {
    int result = reinterpret_cast<const OpenMM::NonbondedForce*>(target)->getNumGlobalParameters();
    return result;
}
OPENMM_EXPORT int OpenMM_NonbondedForce_getNumParticleParameterOffsets(const OpenMM_NonbondedForce* target) {
    int result = reinterpret_cast<const OpenMM::NonbondedForce*>(target)->getNumParticleParameterOffsets();
    return result;
}
OPENMM_EXPORT int OpenMM_NonbondedForce_getNumExceptionParameterOffsets(const OpenMM_NonbondedForce* target) {
    int result = reinterpret_cast<const OpenMM::NonbondedForce*>(target)->getNumExceptionParameterOffsets();
    return result;
}
OPENMM_EXPORT OpenMM_NonbondedForce_NonbondedMethod OpenMM_NonbondedForce_getNonbondedMethod(const OpenMM_NonbondedForce* target) {
    OpenMM::NonbondedForce::NonbondedMethod result = reinterpret_cast<const OpenMM::NonbondedForce*>(target)->getNonbondedMethod();
    return static_cast<OpenMM_NonbondedForce_NonbondedMethod>(result);
}
OPENMM_EXPORT void OpenMM_NonbondedForce_setNonbondedMethod(OpenMM_NonbondedForce* target, OpenMM_NonbondedForce_NonbondedMethod method) {
    reinterpret_cast<OpenMM::NonbondedForce*>(target)->setNonbondedMethod(static_cast<OpenMM::NonbondedForce::NonbondedMethod>(method));
}
OPENMM_EXPORT double OpenMM_NonbondedForce_getCutoffDistance(const OpenMM_NonbondedForce* target) {
    double result = reinterpret_cast<const OpenMM::NonbondedForce*>(target)->getCutoffDistance();
    return result;
}
OPENMM_EXPORT void OpenMM_NonbondedForce_setCutoffDistance(OpenMM_NonbondedForce* target, double distance) {
    reinterpret_cast<OpenMM::NonbondedForce*>(target)->setCutoffDistance(distance);
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_NonbondedForce_getUseSwitchingFunction(const OpenMM_NonbondedForce* target) {
    bool result = reinterpret_cast<const OpenMM::NonbondedForce*>(target)->getUseSwitchingFunction();
    return (result ? OpenMM_True : OpenMM_False);
}
OPENMM_EXPORT void OpenMM_NonbondedForce_setUseSwitchingFunction(OpenMM_NonbondedForce* target, OpenMM_Boolean use) {
    reinterpret_cast<OpenMM::NonbondedForce*>(target)->setUseSwitchingFunction(use);
}
OPENMM_EXPORT double OpenMM_NonbondedForce_getSwitchingDistance(const OpenMM_NonbondedForce* target) {
    double result = reinterpret_cast<const OpenMM::NonbondedForce*>(target)->getSwitchingDistance();
    return result;
}
OPENMM_EXPORT void OpenMM_NonbondedForce_setSwitchingDistance(OpenMM_NonbondedForce* target, double distance) {
    reinterpret_cast<OpenMM::NonbondedForce*>(target)->setSwitchingDistance(distance);
}
OPENMM_EXPORT double OpenMM_NonbondedForce_getReactionFieldDielectric(const OpenMM_NonbondedForce* target) {
    double result = reinterpret_cast<const OpenMM::NonbondedForce*>(target)->getReactionFieldDielectric();
    return result;
}
OPENMM_EXPORT void OpenMM_NonbondedForce_setReactionFieldDielectric(OpenMM_NonbondedForce* target, double dielectric) {
    reinterpret_cast<OpenMM::NonbondedForce*>(target)->setReactionFieldDielectric(dielectric);
}
OPENMM_EXPORT double OpenMM_NonbondedForce_getEwaldErrorTolerance(const OpenMM_NonbondedForce* target) {
    double result = reinterpret_cast<const OpenMM::NonbondedForce*>(target)->getEwaldErrorTolerance();
    return result;
}
OPENMM_EXPORT void OpenMM_NonbondedForce_setEwaldErrorTolerance(OpenMM_NonbondedForce* target, double tol) {
    reinterpret_cast<OpenMM::NonbondedForce*>(target)->setEwaldErrorTolerance(tol);
}
OPENMM_EXPORT void OpenMM_NonbondedForce_getPMEParameters(const OpenMM_NonbondedForce* target, double* alpha, int* nx, int* ny, int* nz) {
    reinterpret_cast<const OpenMM::NonbondedForce*>(target)->getPMEParameters(*reinterpret_cast<double*>(alpha), *reinterpret_cast<int*>(nx), *reinterpret_cast<int*>(ny), *reinterpret_cast<int*>(nz));
}
OPENMM_EXPORT void OpenMM_NonbondedForce_getLJPMEParameters(const OpenMM_NonbondedForce* target, double* alpha, int* nx, int* ny, int* nz) {
    reinterpret_cast<const OpenMM::NonbondedForce*>(target)->getLJPMEParameters(*reinterpret_cast<double*>(alpha), *reinterpret_cast<int*>(nx), *reinterpret_cast<int*>(ny), *reinterpret_cast<int*>(nz));
}
OPENMM_EXPORT void OpenMM_NonbondedForce_setPMEParameters(OpenMM_NonbondedForce* target, double alpha, int nx, int ny, int nz) {
    reinterpret_cast<OpenMM::NonbondedForce*>(target)->setPMEParameters(alpha, nx, ny, nz);
}
OPENMM_EXPORT void OpenMM_NonbondedForce_setLJPMEParameters(OpenMM_NonbondedForce* target, double alpha, int nx, int ny, int nz) {
    reinterpret_cast<OpenMM::NonbondedForce*>(target)->setLJPMEParameters(alpha, nx, ny, nz);
}
OPENMM_EXPORT void OpenMM_NonbondedForce_getPMEParametersInContext(const OpenMM_NonbondedForce* target, const OpenMM_Context* context, double* alpha, int* nx, int* ny, int* nz) {
    reinterpret_cast<const OpenMM::NonbondedForce*>(target)->getPMEParametersInContext(*reinterpret_cast<const Context*>(context), *reinterpret_cast<double*>(alpha), *reinterpret_cast<int*>(nx), *reinterpret_cast<int*>(ny), *reinterpret_cast<int*>(nz));
}
OPENMM_EXPORT void OpenMM_NonbondedForce_getLJPMEParametersInContext(const OpenMM_NonbondedForce* target, const OpenMM_Context* context, double* alpha, int* nx, int* ny, int* nz) {
    reinterpret_cast<const OpenMM::NonbondedForce*>(target)->getLJPMEParametersInContext(*reinterpret_cast<const Context*>(context), *reinterpret_cast<double*>(alpha), *reinterpret_cast<int*>(nx), *reinterpret_cast<int*>(ny), *reinterpret_cast<int*>(nz));
}
OPENMM_EXPORT int OpenMM_NonbondedForce_addParticle(OpenMM_NonbondedForce* target, double charge, double sigma, double epsilon) {
    int result = reinterpret_cast<OpenMM::NonbondedForce*>(target)->addParticle(charge, sigma, epsilon);
    return result;
}
OPENMM_EXPORT void OpenMM_NonbondedForce_getParticleParameters(const OpenMM_NonbondedForce* target, int index, double* charge, double* sigma, double* epsilon) {
    reinterpret_cast<const OpenMM::NonbondedForce*>(target)->getParticleParameters(index, *reinterpret_cast<double*>(charge), *reinterpret_cast<double*>(sigma), *reinterpret_cast<double*>(epsilon));
}
OPENMM_EXPORT void OpenMM_NonbondedForce_setParticleParameters(OpenMM_NonbondedForce* target, int index, double charge, double sigma, double epsilon) {
    reinterpret_cast<OpenMM::NonbondedForce*>(target)->setParticleParameters(index, charge, sigma, epsilon);
}
OPENMM_EXPORT int OpenMM_NonbondedForce_addException(OpenMM_NonbondedForce* target, int particle1, int particle2, double chargeProd, double sigma, double epsilon, OpenMM_Boolean replace) {
    int result = reinterpret_cast<OpenMM::NonbondedForce*>(target)->addException(particle1, particle2, chargeProd, sigma, epsilon, replace);
    return result;
}
OPENMM_EXPORT void OpenMM_NonbondedForce_getExceptionParameters(const OpenMM_NonbondedForce* target, int index, int* particle1, int* particle2, double* chargeProd, double* sigma, double* epsilon) {
    reinterpret_cast<const OpenMM::NonbondedForce*>(target)->getExceptionParameters(index, *reinterpret_cast<int*>(particle1), *reinterpret_cast<int*>(particle2), *reinterpret_cast<double*>(chargeProd), *reinterpret_cast<double*>(sigma), *reinterpret_cast<double*>(epsilon));
}
OPENMM_EXPORT void OpenMM_NonbondedForce_setExceptionParameters(OpenMM_NonbondedForce* target, int index, int particle1, int particle2, double chargeProd, double sigma, double epsilon) {
    reinterpret_cast<OpenMM::NonbondedForce*>(target)->setExceptionParameters(index, particle1, particle2, chargeProd, sigma, epsilon);
}
OPENMM_EXPORT void OpenMM_NonbondedForce_createExceptionsFromBonds(OpenMM_NonbondedForce* target, const OpenMM_BondArray* bonds, double coulomb14Scale, double lj14Scale) {
    reinterpret_cast<OpenMM::NonbondedForce*>(target)->createExceptionsFromBonds(*reinterpret_cast<const std::vector< std::pair< int, int > >*>(bonds), coulomb14Scale, lj14Scale);
}
OPENMM_EXPORT int OpenMM_NonbondedForce_addGlobalParameter(OpenMM_NonbondedForce* target, const char* name, double defaultValue) {
    int result = reinterpret_cast<OpenMM::NonbondedForce*>(target)->addGlobalParameter(std::string(name), defaultValue);
    return result;
}
OPENMM_EXPORT const char* OpenMM_NonbondedForce_getGlobalParameterName(const OpenMM_NonbondedForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::NonbondedForce*>(target)->getGlobalParameterName(index);
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_NonbondedForce_setGlobalParameterName(OpenMM_NonbondedForce* target, int index, const char* name) {
    reinterpret_cast<OpenMM::NonbondedForce*>(target)->setGlobalParameterName(index, std::string(name));
}
OPENMM_EXPORT double OpenMM_NonbondedForce_getGlobalParameterDefaultValue(const OpenMM_NonbondedForce* target, int index) {
    double result = reinterpret_cast<const OpenMM::NonbondedForce*>(target)->getGlobalParameterDefaultValue(index);
    return result;
}
OPENMM_EXPORT void OpenMM_NonbondedForce_setGlobalParameterDefaultValue(OpenMM_NonbondedForce* target, int index, double defaultValue) {
    reinterpret_cast<OpenMM::NonbondedForce*>(target)->setGlobalParameterDefaultValue(index, defaultValue);
}
OPENMM_EXPORT int OpenMM_NonbondedForce_addParticleParameterOffset(OpenMM_NonbondedForce* target, const char* parameter, int particleIndex, double chargeScale, double sigmaScale, double epsilonScale) {
    int result = reinterpret_cast<OpenMM::NonbondedForce*>(target)->addParticleParameterOffset(std::string(parameter), particleIndex, chargeScale, sigmaScale, epsilonScale);
    return result;
}
OPENMM_EXPORT void OpenMM_NonbondedForce_getParticleParameterOffset(const OpenMM_NonbondedForce* target, int index, char** parameter, int* particleIndex, double* chargeScale, double* sigmaScale, double* epsilonScale) {
    reinterpret_cast<const OpenMM::NonbondedForce*>(target)->getParticleParameterOffset(index, *reinterpret_cast<std::string*>(parameter), *reinterpret_cast<int*>(particleIndex), *reinterpret_cast<double*>(chargeScale), *reinterpret_cast<double*>(sigmaScale), *reinterpret_cast<double*>(epsilonScale));
}
OPENMM_EXPORT void OpenMM_NonbondedForce_setParticleParameterOffset(OpenMM_NonbondedForce* target, int index, const char* parameter, int particleIndex, double chargeScale, double sigmaScale, double epsilonScale) {
    reinterpret_cast<OpenMM::NonbondedForce*>(target)->setParticleParameterOffset(index, std::string(parameter), particleIndex, chargeScale, sigmaScale, epsilonScale);
}
OPENMM_EXPORT int OpenMM_NonbondedForce_addExceptionParameterOffset(OpenMM_NonbondedForce* target, const char* parameter, int exceptionIndex, double chargeProdScale, double sigmaScale, double epsilonScale) {
    int result = reinterpret_cast<OpenMM::NonbondedForce*>(target)->addExceptionParameterOffset(std::string(parameter), exceptionIndex, chargeProdScale, sigmaScale, epsilonScale);
    return result;
}
OPENMM_EXPORT void OpenMM_NonbondedForce_getExceptionParameterOffset(const OpenMM_NonbondedForce* target, int index, char** parameter, int* exceptionIndex, double* chargeProdScale, double* sigmaScale, double* epsilonScale) {
    reinterpret_cast<const OpenMM::NonbondedForce*>(target)->getExceptionParameterOffset(index, *reinterpret_cast<std::string*>(parameter), *reinterpret_cast<int*>(exceptionIndex), *reinterpret_cast<double*>(chargeProdScale), *reinterpret_cast<double*>(sigmaScale), *reinterpret_cast<double*>(epsilonScale));
}
OPENMM_EXPORT void OpenMM_NonbondedForce_setExceptionParameterOffset(OpenMM_NonbondedForce* target, int index, const char* parameter, int exceptionIndex, double chargeProdScale, double sigmaScale, double epsilonScale) {
    reinterpret_cast<OpenMM::NonbondedForce*>(target)->setExceptionParameterOffset(index, std::string(parameter), exceptionIndex, chargeProdScale, sigmaScale, epsilonScale);
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_NonbondedForce_getUseDispersionCorrection(const OpenMM_NonbondedForce* target) {
    bool result = reinterpret_cast<const OpenMM::NonbondedForce*>(target)->getUseDispersionCorrection();
    return (result ? OpenMM_True : OpenMM_False);
}
OPENMM_EXPORT void OpenMM_NonbondedForce_setUseDispersionCorrection(OpenMM_NonbondedForce* target, OpenMM_Boolean useCorrection) {
    reinterpret_cast<OpenMM::NonbondedForce*>(target)->setUseDispersionCorrection(useCorrection);
}
OPENMM_EXPORT int OpenMM_NonbondedForce_getReciprocalSpaceForceGroup(const OpenMM_NonbondedForce* target) {
    int result = reinterpret_cast<const OpenMM::NonbondedForce*>(target)->getReciprocalSpaceForceGroup();
    return result;
}
OPENMM_EXPORT void OpenMM_NonbondedForce_setReciprocalSpaceForceGroup(OpenMM_NonbondedForce* target, int group) {
    reinterpret_cast<OpenMM::NonbondedForce*>(target)->setReciprocalSpaceForceGroup(group);
}
OPENMM_EXPORT void OpenMM_NonbondedForce_updateParametersInContext(OpenMM_NonbondedForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::NonbondedForce*>(target)->updateParametersInContext(*reinterpret_cast<OpenMM::Context*>(context));
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_NonbondedForce_usesPeriodicBoundaryConditions(const OpenMM_NonbondedForce* target) {
    bool result = reinterpret_cast<const OpenMM::NonbondedForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_NonbondedForce_getExceptionsUsePeriodicBoundaryConditions(const OpenMM_NonbondedForce* target) {
    bool result = reinterpret_cast<const OpenMM::NonbondedForce*>(target)->getExceptionsUsePeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}
OPENMM_EXPORT void OpenMM_NonbondedForce_setExceptionsUsePeriodicBoundaryConditions(OpenMM_NonbondedForce* target, OpenMM_Boolean periodic) {
    reinterpret_cast<OpenMM::NonbondedForce*>(target)->setExceptionsUsePeriodicBoundaryConditions(periodic);
}

/* OpenMM::NoseHooverChain */
OPENMM_EXPORT OpenMM_NoseHooverChain* OpenMM_NoseHooverChain_create(double temperature, double relativeTemperature, double collisionFrequency, double relativeCollisionFrequency, int numDOFs, int chainLength, int numMTS, int numYoshidaSuzuki, int chainID, const OpenMM_IntArray* thermostatedAtoms, const OpenMM_BondArray* thermostatedPairs) {
    return reinterpret_cast<OpenMM_NoseHooverChain*>(new OpenMM::NoseHooverChain(temperature, relativeTemperature, collisionFrequency, relativeCollisionFrequency, numDOFs, chainLength, numMTS, numYoshidaSuzuki, chainID, *reinterpret_cast<const std::vector< int >*>(thermostatedAtoms), *reinterpret_cast<const std::vector< std::pair< int, int > >*>(thermostatedPairs)));
}
OPENMM_EXPORT void OpenMM_NoseHooverChain_destroy(OpenMM_NoseHooverChain* target) {
    delete reinterpret_cast<OpenMM::NoseHooverChain*>(target);
}
OPENMM_EXPORT double OpenMM_NoseHooverChain_getTemperature(const OpenMM_NoseHooverChain* target) {
    double result = reinterpret_cast<const OpenMM::NoseHooverChain*>(target)->getTemperature();
    return result;
}
OPENMM_EXPORT void OpenMM_NoseHooverChain_setTemperature(OpenMM_NoseHooverChain* target, double temperature) {
    reinterpret_cast<OpenMM::NoseHooverChain*>(target)->setTemperature(temperature);
}
OPENMM_EXPORT double OpenMM_NoseHooverChain_getRelativeTemperature(const OpenMM_NoseHooverChain* target) {
    double result = reinterpret_cast<const OpenMM::NoseHooverChain*>(target)->getRelativeTemperature();
    return result;
}
OPENMM_EXPORT void OpenMM_NoseHooverChain_setRelativeTemperature(OpenMM_NoseHooverChain* target, double temperature) {
    reinterpret_cast<OpenMM::NoseHooverChain*>(target)->setRelativeTemperature(temperature);
}
OPENMM_EXPORT double OpenMM_NoseHooverChain_getCollisionFrequency(const OpenMM_NoseHooverChain* target) {
    double result = reinterpret_cast<const OpenMM::NoseHooverChain*>(target)->getCollisionFrequency();
    return result;
}
OPENMM_EXPORT void OpenMM_NoseHooverChain_setCollisionFrequency(OpenMM_NoseHooverChain* target, double frequency) {
    reinterpret_cast<OpenMM::NoseHooverChain*>(target)->setCollisionFrequency(frequency);
}
OPENMM_EXPORT double OpenMM_NoseHooverChain_getRelativeCollisionFrequency(const OpenMM_NoseHooverChain* target) {
    double result = reinterpret_cast<const OpenMM::NoseHooverChain*>(target)->getRelativeCollisionFrequency();
    return result;
}
OPENMM_EXPORT void OpenMM_NoseHooverChain_setRelativeCollisionFrequency(OpenMM_NoseHooverChain* target, double frequency) {
    reinterpret_cast<OpenMM::NoseHooverChain*>(target)->setRelativeCollisionFrequency(frequency);
}
OPENMM_EXPORT int OpenMM_NoseHooverChain_getNumDegreesOfFreedom(const OpenMM_NoseHooverChain* target) {
    int result = reinterpret_cast<const OpenMM::NoseHooverChain*>(target)->getNumDegreesOfFreedom();
    return result;
}
OPENMM_EXPORT void OpenMM_NoseHooverChain_setNumDegreesOfFreedom(OpenMM_NoseHooverChain* target, int numDOF) {
    reinterpret_cast<OpenMM::NoseHooverChain*>(target)->setNumDegreesOfFreedom(numDOF);
}
OPENMM_EXPORT int OpenMM_NoseHooverChain_getChainLength(const OpenMM_NoseHooverChain* target) {
    int result = reinterpret_cast<const OpenMM::NoseHooverChain*>(target)->getChainLength();
    return result;
}
OPENMM_EXPORT int OpenMM_NoseHooverChain_getNumMultiTimeSteps(const OpenMM_NoseHooverChain* target) {
    int result = reinterpret_cast<const OpenMM::NoseHooverChain*>(target)->getNumMultiTimeSteps();
    return result;
}
OPENMM_EXPORT int OpenMM_NoseHooverChain_getNumYoshidaSuzukiTimeSteps(const OpenMM_NoseHooverChain* target) {
    int result = reinterpret_cast<const OpenMM::NoseHooverChain*>(target)->getNumYoshidaSuzukiTimeSteps();
    return result;
}
OPENMM_EXPORT int OpenMM_NoseHooverChain_getChainID(const OpenMM_NoseHooverChain* target) {
    int result = reinterpret_cast<const OpenMM::NoseHooverChain*>(target)->getChainID();
    return result;
}
OPENMM_EXPORT const OpenMM_IntArray* OpenMM_NoseHooverChain_getThermostatedAtoms(const OpenMM_NoseHooverChain* target) {
    const std::vector< int >* result = &reinterpret_cast<const OpenMM::NoseHooverChain*>(target)->getThermostatedAtoms();
    return reinterpret_cast<const OpenMM_IntArray*>(result);
}
OPENMM_EXPORT void OpenMM_NoseHooverChain_setThermostatedAtoms(OpenMM_NoseHooverChain* target, const OpenMM_IntArray* atomIDs) {
    reinterpret_cast<OpenMM::NoseHooverChain*>(target)->setThermostatedAtoms(*reinterpret_cast<const std::vector< int >*>(atomIDs));
}
OPENMM_EXPORT const OpenMM_BondArray* OpenMM_NoseHooverChain_getThermostatedPairs(const OpenMM_NoseHooverChain* target) {
    const std::vector< std::pair< int, int > >* result = &reinterpret_cast<const OpenMM::NoseHooverChain*>(target)->getThermostatedPairs();
    return reinterpret_cast<const OpenMM_BondArray*>(result);
}
OPENMM_EXPORT void OpenMM_NoseHooverChain_setThermostatedPairs(OpenMM_NoseHooverChain* target, const OpenMM_BondArray* pairIDs) {
    reinterpret_cast<OpenMM::NoseHooverChain*>(target)->setThermostatedPairs(*reinterpret_cast<const std::vector< std::pair< int, int > >*>(pairIDs));
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_NoseHooverChain_usesPeriodicBoundaryConditions(const OpenMM_NoseHooverChain* target) {
    bool result = reinterpret_cast<const OpenMM::NoseHooverChain*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::HarmonicBondForce */
OPENMM_EXPORT OpenMM_HarmonicBondForce* OpenMM_HarmonicBondForce_create() {
    return reinterpret_cast<OpenMM_HarmonicBondForce*>(new OpenMM::HarmonicBondForce());
}
OPENMM_EXPORT void OpenMM_HarmonicBondForce_destroy(OpenMM_HarmonicBondForce* target) {
    delete reinterpret_cast<OpenMM::HarmonicBondForce*>(target);
}
OPENMM_EXPORT int OpenMM_HarmonicBondForce_getNumBonds(const OpenMM_HarmonicBondForce* target) {
    int result = reinterpret_cast<const OpenMM::HarmonicBondForce*>(target)->getNumBonds();
    return result;
}
OPENMM_EXPORT int OpenMM_HarmonicBondForce_addBond(OpenMM_HarmonicBondForce* target, int particle1, int particle2, double length, double k) {
    int result = reinterpret_cast<OpenMM::HarmonicBondForce*>(target)->addBond(particle1, particle2, length, k);
    return result;
}
OPENMM_EXPORT void OpenMM_HarmonicBondForce_getBondParameters(const OpenMM_HarmonicBondForce* target, int index, int* particle1, int* particle2, double* length, double* k) {
    reinterpret_cast<const OpenMM::HarmonicBondForce*>(target)->getBondParameters(index, *reinterpret_cast<int*>(particle1), *reinterpret_cast<int*>(particle2), *reinterpret_cast<double*>(length), *reinterpret_cast<double*>(k));
}
OPENMM_EXPORT void OpenMM_HarmonicBondForce_setBondParameters(OpenMM_HarmonicBondForce* target, int index, int particle1, int particle2, double length, double k) {
    reinterpret_cast<OpenMM::HarmonicBondForce*>(target)->setBondParameters(index, particle1, particle2, length, k);
}
OPENMM_EXPORT void OpenMM_HarmonicBondForce_updateParametersInContext(OpenMM_HarmonicBondForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::HarmonicBondForce*>(target)->updateParametersInContext(*reinterpret_cast<OpenMM::Context*>(context));
}
OPENMM_EXPORT void OpenMM_HarmonicBondForce_setUsesPeriodicBoundaryConditions(OpenMM_HarmonicBondForce* target, OpenMM_Boolean periodic) {
    reinterpret_cast<OpenMM::HarmonicBondForce*>(target)->setUsesPeriodicBoundaryConditions(periodic);
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_HarmonicBondForce_usesPeriodicBoundaryConditions(const OpenMM_HarmonicBondForce* target) {
    bool result = reinterpret_cast<const OpenMM::HarmonicBondForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::CustomNonbondedForce */
OPENMM_EXPORT OpenMM_CustomNonbondedForce* OpenMM_CustomNonbondedForce_create(const char* energy) {
    return reinterpret_cast<OpenMM_CustomNonbondedForce*>(new OpenMM::CustomNonbondedForce(std::string(energy)));
}
OPENMM_EXPORT OpenMM_CustomNonbondedForce* OpenMM_CustomNonbondedForce_create_2(const OpenMM_CustomNonbondedForce* rhs) {
    return reinterpret_cast<OpenMM_CustomNonbondedForce*>(new OpenMM::CustomNonbondedForce(*reinterpret_cast<const CustomNonbondedForce*>(rhs)));
}
OPENMM_EXPORT void OpenMM_CustomNonbondedForce_destroy(OpenMM_CustomNonbondedForce* target) {
    delete reinterpret_cast<OpenMM::CustomNonbondedForce*>(target);
}
OPENMM_EXPORT int OpenMM_CustomNonbondedForce_getNumParticles(const OpenMM_CustomNonbondedForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomNonbondedForce*>(target)->getNumParticles();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomNonbondedForce_getNumExclusions(const OpenMM_CustomNonbondedForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomNonbondedForce*>(target)->getNumExclusions();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomNonbondedForce_getNumPerParticleParameters(const OpenMM_CustomNonbondedForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomNonbondedForce*>(target)->getNumPerParticleParameters();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomNonbondedForce_getNumGlobalParameters(const OpenMM_CustomNonbondedForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomNonbondedForce*>(target)->getNumGlobalParameters();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomNonbondedForce_getNumTabulatedFunctions(const OpenMM_CustomNonbondedForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomNonbondedForce*>(target)->getNumTabulatedFunctions();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomNonbondedForce_getNumFunctions(const OpenMM_CustomNonbondedForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomNonbondedForce*>(target)->getNumFunctions();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomNonbondedForce_getNumInteractionGroups(const OpenMM_CustomNonbondedForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomNonbondedForce*>(target)->getNumInteractionGroups();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomNonbondedForce_getNumEnergyParameterDerivatives(const OpenMM_CustomNonbondedForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomNonbondedForce*>(target)->getNumEnergyParameterDerivatives();
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomNonbondedForce_getEnergyFunction(const OpenMM_CustomNonbondedForce* target) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomNonbondedForce*>(target)->getEnergyFunction();
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomNonbondedForce_setEnergyFunction(OpenMM_CustomNonbondedForce* target, const char* energy) {
    reinterpret_cast<OpenMM::CustomNonbondedForce*>(target)->setEnergyFunction(std::string(energy));
}
OPENMM_EXPORT OpenMM_CustomNonbondedForce_NonbondedMethod OpenMM_CustomNonbondedForce_getNonbondedMethod(const OpenMM_CustomNonbondedForce* target) {
    OpenMM::CustomNonbondedForce::NonbondedMethod result = reinterpret_cast<const OpenMM::CustomNonbondedForce*>(target)->getNonbondedMethod();
    return static_cast<OpenMM_CustomNonbondedForce_NonbondedMethod>(result);
}
OPENMM_EXPORT void OpenMM_CustomNonbondedForce_setNonbondedMethod(OpenMM_CustomNonbondedForce* target, OpenMM_CustomNonbondedForce_NonbondedMethod method) {
    reinterpret_cast<OpenMM::CustomNonbondedForce*>(target)->setNonbondedMethod(static_cast<OpenMM::CustomNonbondedForce::NonbondedMethod>(method));
}
OPENMM_EXPORT double OpenMM_CustomNonbondedForce_getCutoffDistance(const OpenMM_CustomNonbondedForce* target) {
    double result = reinterpret_cast<const OpenMM::CustomNonbondedForce*>(target)->getCutoffDistance();
    return result;
}
OPENMM_EXPORT void OpenMM_CustomNonbondedForce_setCutoffDistance(OpenMM_CustomNonbondedForce* target, double distance) {
    reinterpret_cast<OpenMM::CustomNonbondedForce*>(target)->setCutoffDistance(distance);
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_CustomNonbondedForce_getUseSwitchingFunction(const OpenMM_CustomNonbondedForce* target) {
    bool result = reinterpret_cast<const OpenMM::CustomNonbondedForce*>(target)->getUseSwitchingFunction();
    return (result ? OpenMM_True : OpenMM_False);
}
OPENMM_EXPORT void OpenMM_CustomNonbondedForce_setUseSwitchingFunction(OpenMM_CustomNonbondedForce* target, OpenMM_Boolean use) {
    reinterpret_cast<OpenMM::CustomNonbondedForce*>(target)->setUseSwitchingFunction(use);
}
OPENMM_EXPORT double OpenMM_CustomNonbondedForce_getSwitchingDistance(const OpenMM_CustomNonbondedForce* target) {
    double result = reinterpret_cast<const OpenMM::CustomNonbondedForce*>(target)->getSwitchingDistance();
    return result;
}
OPENMM_EXPORT void OpenMM_CustomNonbondedForce_setSwitchingDistance(OpenMM_CustomNonbondedForce* target, double distance) {
    reinterpret_cast<OpenMM::CustomNonbondedForce*>(target)->setSwitchingDistance(distance);
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_CustomNonbondedForce_getUseLongRangeCorrection(const OpenMM_CustomNonbondedForce* target) {
    bool result = reinterpret_cast<const OpenMM::CustomNonbondedForce*>(target)->getUseLongRangeCorrection();
    return (result ? OpenMM_True : OpenMM_False);
}
OPENMM_EXPORT void OpenMM_CustomNonbondedForce_setUseLongRangeCorrection(OpenMM_CustomNonbondedForce* target, OpenMM_Boolean use) {
    reinterpret_cast<OpenMM::CustomNonbondedForce*>(target)->setUseLongRangeCorrection(use);
}
OPENMM_EXPORT int OpenMM_CustomNonbondedForce_addPerParticleParameter(OpenMM_CustomNonbondedForce* target, const char* name) {
    int result = reinterpret_cast<OpenMM::CustomNonbondedForce*>(target)->addPerParticleParameter(std::string(name));
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomNonbondedForce_getPerParticleParameterName(const OpenMM_CustomNonbondedForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomNonbondedForce*>(target)->getPerParticleParameterName(index);
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomNonbondedForce_setPerParticleParameterName(OpenMM_CustomNonbondedForce* target, int index, const char* name) {
    reinterpret_cast<OpenMM::CustomNonbondedForce*>(target)->setPerParticleParameterName(index, std::string(name));
}
OPENMM_EXPORT int OpenMM_CustomNonbondedForce_addGlobalParameter(OpenMM_CustomNonbondedForce* target, const char* name, double defaultValue) {
    int result = reinterpret_cast<OpenMM::CustomNonbondedForce*>(target)->addGlobalParameter(std::string(name), defaultValue);
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomNonbondedForce_getGlobalParameterName(const OpenMM_CustomNonbondedForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomNonbondedForce*>(target)->getGlobalParameterName(index);
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomNonbondedForce_setGlobalParameterName(OpenMM_CustomNonbondedForce* target, int index, const char* name) {
    reinterpret_cast<OpenMM::CustomNonbondedForce*>(target)->setGlobalParameterName(index, std::string(name));
}
OPENMM_EXPORT double OpenMM_CustomNonbondedForce_getGlobalParameterDefaultValue(const OpenMM_CustomNonbondedForce* target, int index) {
    double result = reinterpret_cast<const OpenMM::CustomNonbondedForce*>(target)->getGlobalParameterDefaultValue(index);
    return result;
}
OPENMM_EXPORT void OpenMM_CustomNonbondedForce_setGlobalParameterDefaultValue(OpenMM_CustomNonbondedForce* target, int index, double defaultValue) {
    reinterpret_cast<OpenMM::CustomNonbondedForce*>(target)->setGlobalParameterDefaultValue(index, defaultValue);
}
OPENMM_EXPORT void OpenMM_CustomNonbondedForce_addEnergyParameterDerivative(OpenMM_CustomNonbondedForce* target, const char* name) {
    reinterpret_cast<OpenMM::CustomNonbondedForce*>(target)->addEnergyParameterDerivative(std::string(name));
}
OPENMM_EXPORT const char* OpenMM_CustomNonbondedForce_getEnergyParameterDerivativeName(const OpenMM_CustomNonbondedForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomNonbondedForce*>(target)->getEnergyParameterDerivativeName(index);
    return result->c_str();
}
OPENMM_EXPORT int OpenMM_CustomNonbondedForce_addParticle(OpenMM_CustomNonbondedForce* target, const OpenMM_DoubleArray* parameters) {
    int result = reinterpret_cast<OpenMM::CustomNonbondedForce*>(target)->addParticle(*reinterpret_cast<const std::vector< double >*>(parameters));
    return result;
}
OPENMM_EXPORT void OpenMM_CustomNonbondedForce_getParticleParameters(const OpenMM_CustomNonbondedForce* target, int index, OpenMM_DoubleArray* parameters) {
    reinterpret_cast<const OpenMM::CustomNonbondedForce*>(target)->getParticleParameters(index, *reinterpret_cast<std::vector< double >*>(parameters));
}
OPENMM_EXPORT void OpenMM_CustomNonbondedForce_setParticleParameters(OpenMM_CustomNonbondedForce* target, int index, const OpenMM_DoubleArray* parameters) {
    reinterpret_cast<OpenMM::CustomNonbondedForce*>(target)->setParticleParameters(index, *reinterpret_cast<const std::vector< double >*>(parameters));
}
OPENMM_EXPORT int OpenMM_CustomNonbondedForce_addExclusion(OpenMM_CustomNonbondedForce* target, int particle1, int particle2) {
    int result = reinterpret_cast<OpenMM::CustomNonbondedForce*>(target)->addExclusion(particle1, particle2);
    return result;
}
OPENMM_EXPORT void OpenMM_CustomNonbondedForce_getExclusionParticles(const OpenMM_CustomNonbondedForce* target, int index, int* particle1, int* particle2) {
    reinterpret_cast<const OpenMM::CustomNonbondedForce*>(target)->getExclusionParticles(index, *reinterpret_cast<int*>(particle1), *reinterpret_cast<int*>(particle2));
}
OPENMM_EXPORT void OpenMM_CustomNonbondedForce_setExclusionParticles(OpenMM_CustomNonbondedForce* target, int index, int particle1, int particle2) {
    reinterpret_cast<OpenMM::CustomNonbondedForce*>(target)->setExclusionParticles(index, particle1, particle2);
}
OPENMM_EXPORT void OpenMM_CustomNonbondedForce_createExclusionsFromBonds(OpenMM_CustomNonbondedForce* target, const OpenMM_BondArray* bonds, int bondCutoff) {
    reinterpret_cast<OpenMM::CustomNonbondedForce*>(target)->createExclusionsFromBonds(*reinterpret_cast<const std::vector< std::pair< int, int > >*>(bonds), bondCutoff);
}
OPENMM_EXPORT int OpenMM_CustomNonbondedForce_addTabulatedFunction(OpenMM_CustomNonbondedForce* target, const char* name, OpenMM_TabulatedFunction* function) {
    int result = reinterpret_cast<OpenMM::CustomNonbondedForce*>(target)->addTabulatedFunction(std::string(name), reinterpret_cast<TabulatedFunction *>(function));
    return result;
}
OPENMM_EXPORT OpenMM_TabulatedFunction* OpenMM_CustomNonbondedForce_getTabulatedFunction(OpenMM_CustomNonbondedForce* target, int index) {
    TabulatedFunction* result = &reinterpret_cast<OpenMM::CustomNonbondedForce*>(target)->getTabulatedFunction(index);
    return reinterpret_cast<OpenMM_TabulatedFunction*>(result);
}
OPENMM_EXPORT const char* OpenMM_CustomNonbondedForce_getTabulatedFunctionName(const OpenMM_CustomNonbondedForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomNonbondedForce*>(target)->getTabulatedFunctionName(index);
    return result->c_str();
}
OPENMM_EXPORT int OpenMM_CustomNonbondedForce_addFunction(OpenMM_CustomNonbondedForce* target, const char* name, const OpenMM_DoubleArray* values, double min, double max) {
    int result = reinterpret_cast<OpenMM::CustomNonbondedForce*>(target)->addFunction(std::string(name), *reinterpret_cast<const std::vector< double >*>(values), min, max);
    return result;
}
OPENMM_EXPORT void OpenMM_CustomNonbondedForce_getFunctionParameters(const OpenMM_CustomNonbondedForce* target, int index, char** name, OpenMM_DoubleArray* values, double* min, double* max) {
    reinterpret_cast<const OpenMM::CustomNonbondedForce*>(target)->getFunctionParameters(index, *reinterpret_cast<std::string*>(name), *reinterpret_cast<std::vector< double >*>(values), *reinterpret_cast<double*>(min), *reinterpret_cast<double*>(max));
}
OPENMM_EXPORT void OpenMM_CustomNonbondedForce_setFunctionParameters(OpenMM_CustomNonbondedForce* target, int index, const char* name, const OpenMM_DoubleArray* values, double min, double max) {
    reinterpret_cast<OpenMM::CustomNonbondedForce*>(target)->setFunctionParameters(index, std::string(name), *reinterpret_cast<const std::vector< double >*>(values), min, max);
}
OPENMM_EXPORT int OpenMM_CustomNonbondedForce_addInteractionGroup(OpenMM_CustomNonbondedForce* target, const OpenMM_IntSet* set1, const OpenMM_IntSet* set2) {
    int result = reinterpret_cast<OpenMM::CustomNonbondedForce*>(target)->addInteractionGroup(*reinterpret_cast<const std::set< int >*>(set1), *reinterpret_cast<const std::set< int >*>(set2));
    return result;
}
OPENMM_EXPORT void OpenMM_CustomNonbondedForce_getInteractionGroupParameters(const OpenMM_CustomNonbondedForce* target, int index, OpenMM_IntSet* set1, OpenMM_IntSet* set2) {
    reinterpret_cast<const OpenMM::CustomNonbondedForce*>(target)->getInteractionGroupParameters(index, *reinterpret_cast<std::set< int >*>(set1), *reinterpret_cast<std::set< int >*>(set2));
}
OPENMM_EXPORT void OpenMM_CustomNonbondedForce_setInteractionGroupParameters(OpenMM_CustomNonbondedForce* target, int index, const OpenMM_IntSet* set1, const OpenMM_IntSet* set2) {
    reinterpret_cast<OpenMM::CustomNonbondedForce*>(target)->setInteractionGroupParameters(index, *reinterpret_cast<const std::set< int >*>(set1), *reinterpret_cast<const std::set< int >*>(set2));
}
OPENMM_EXPORT void OpenMM_CustomNonbondedForce_updateParametersInContext(OpenMM_CustomNonbondedForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::CustomNonbondedForce*>(target)->updateParametersInContext(*reinterpret_cast<OpenMM::Context*>(context));
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_CustomNonbondedForce_usesPeriodicBoundaryConditions(const OpenMM_CustomNonbondedForce* target) {
    bool result = reinterpret_cast<const OpenMM::CustomNonbondedForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::System */
OPENMM_EXPORT OpenMM_System* OpenMM_System_create() {
    return reinterpret_cast<OpenMM_System*>(new OpenMM::System());
}
OPENMM_EXPORT void OpenMM_System_destroy(OpenMM_System* target) {
    delete reinterpret_cast<OpenMM::System*>(target);
}
OPENMM_EXPORT int OpenMM_System_getNumParticles(const OpenMM_System* target) {
    int result = reinterpret_cast<const OpenMM::System*>(target)->getNumParticles();
    return result;
}
OPENMM_EXPORT int OpenMM_System_addParticle(OpenMM_System* target, double mass) {
    int result = reinterpret_cast<OpenMM::System*>(target)->addParticle(mass);
    return result;
}
OPENMM_EXPORT double OpenMM_System_getParticleMass(const OpenMM_System* target, int index) {
    double result = reinterpret_cast<const OpenMM::System*>(target)->getParticleMass(index);
    return result;
}
OPENMM_EXPORT void OpenMM_System_setParticleMass(OpenMM_System* target, int index, double mass) {
    reinterpret_cast<OpenMM::System*>(target)->setParticleMass(index, mass);
}
OPENMM_EXPORT void OpenMM_System_setVirtualSite(OpenMM_System* target, int index, OpenMM_VirtualSite* virtualSite) {
    reinterpret_cast<OpenMM::System*>(target)->setVirtualSite(index, reinterpret_cast<VirtualSite *>(virtualSite));
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_System_isVirtualSite(const OpenMM_System* target, int index) {
    bool result = reinterpret_cast<const OpenMM::System*>(target)->isVirtualSite(index);
    return (result ? OpenMM_True : OpenMM_False);
}
OPENMM_EXPORT const OpenMM_VirtualSite* OpenMM_System_getVirtualSite(const OpenMM_System* target, int index) {
    const VirtualSite* result = &reinterpret_cast<const OpenMM::System*>(target)->getVirtualSite(index);
    return reinterpret_cast<const OpenMM_VirtualSite*>(result);
}
OPENMM_EXPORT int OpenMM_System_getNumConstraints(const OpenMM_System* target) {
    int result = reinterpret_cast<const OpenMM::System*>(target)->getNumConstraints();
    return result;
}
OPENMM_EXPORT int OpenMM_System_addConstraint(OpenMM_System* target, int particle1, int particle2, double distance) {
    int result = reinterpret_cast<OpenMM::System*>(target)->addConstraint(particle1, particle2, distance);
    return result;
}
OPENMM_EXPORT void OpenMM_System_getConstraintParameters(const OpenMM_System* target, int index, int* particle1, int* particle2, double* distance) {
    reinterpret_cast<const OpenMM::System*>(target)->getConstraintParameters(index, *reinterpret_cast<int*>(particle1), *reinterpret_cast<int*>(particle2), *reinterpret_cast<double*>(distance));
}
OPENMM_EXPORT void OpenMM_System_setConstraintParameters(OpenMM_System* target, int index, int particle1, int particle2, double distance) {
    reinterpret_cast<OpenMM::System*>(target)->setConstraintParameters(index, particle1, particle2, distance);
}
OPENMM_EXPORT void OpenMM_System_removeConstraint(OpenMM_System* target, int index) {
    reinterpret_cast<OpenMM::System*>(target)->removeConstraint(index);
}
OPENMM_EXPORT int OpenMM_System_addForce(OpenMM_System* target, OpenMM_Force* force) {
    int result = reinterpret_cast<OpenMM::System*>(target)->addForce(reinterpret_cast<Force *>(force));
    return result;
}
OPENMM_EXPORT int OpenMM_System_getNumForces(const OpenMM_System* target) {
    int result = reinterpret_cast<const OpenMM::System*>(target)->getNumForces();
    return result;
}
OPENMM_EXPORT OpenMM_Force* OpenMM_System_getForce(OpenMM_System* target, int index) {
    Force* result = &reinterpret_cast<OpenMM::System*>(target)->getForce(index);
    return reinterpret_cast<OpenMM_Force*>(result);
}
OPENMM_EXPORT void OpenMM_System_removeForce(OpenMM_System* target, int index) {
    reinterpret_cast<OpenMM::System*>(target)->removeForce(index);
}
OPENMM_EXPORT void OpenMM_System_getDefaultPeriodicBoxVectors(const OpenMM_System* target, OpenMM_Vec3* a, OpenMM_Vec3* b, OpenMM_Vec3* c) {
    reinterpret_cast<const OpenMM::System*>(target)->getDefaultPeriodicBoxVectors(*reinterpret_cast<Vec3*>(a), *reinterpret_cast<Vec3*>(b), *reinterpret_cast<Vec3*>(c));
}
OPENMM_EXPORT void OpenMM_System_setDefaultPeriodicBoxVectors(OpenMM_System* target, const OpenMM_Vec3* a, const OpenMM_Vec3* b, const OpenMM_Vec3* c) {
    reinterpret_cast<OpenMM::System*>(target)->setDefaultPeriodicBoxVectors(*reinterpret_cast<const Vec3*>(a), *reinterpret_cast<const Vec3*>(b), *reinterpret_cast<const Vec3*>(c));
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_System_usesPeriodicBoundaryConditions(const OpenMM_System* target) {
    bool result = reinterpret_cast<const OpenMM::System*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::Continuous1DFunction */
OPENMM_EXPORT OpenMM_Continuous1DFunction* OpenMM_Continuous1DFunction_create(const OpenMM_DoubleArray* values, double min, double max) {
    return reinterpret_cast<OpenMM_Continuous1DFunction*>(new OpenMM::Continuous1DFunction(*reinterpret_cast<const std::vector< double >*>(values), min, max));
}
OPENMM_EXPORT void OpenMM_Continuous1DFunction_destroy(OpenMM_Continuous1DFunction* target) {
    delete reinterpret_cast<OpenMM::Continuous1DFunction*>(target);
}
OPENMM_EXPORT void OpenMM_Continuous1DFunction_getFunctionParameters(const OpenMM_Continuous1DFunction* target, OpenMM_DoubleArray* values, double* min, double* max) {
    reinterpret_cast<const OpenMM::Continuous1DFunction*>(target)->getFunctionParameters(*reinterpret_cast<std::vector< double >*>(values), *reinterpret_cast<double*>(min), *reinterpret_cast<double*>(max));
}
OPENMM_EXPORT void OpenMM_Continuous1DFunction_setFunctionParameters(OpenMM_Continuous1DFunction* target, const OpenMM_DoubleArray* values, double min, double max) {
    reinterpret_cast<OpenMM::Continuous1DFunction*>(target)->setFunctionParameters(*reinterpret_cast<const std::vector< double >*>(values), min, max);
}
OPENMM_EXPORT OpenMM_Continuous1DFunction* OpenMM_Continuous1DFunction_Copy(const OpenMM_Continuous1DFunction* target) {
    Continuous1DFunction * result = reinterpret_cast<const OpenMM::Continuous1DFunction*>(target)->Copy();
    return reinterpret_cast<OpenMM_Continuous1DFunction*>(result);
}

/* OpenMM::Platform */
OPENMM_EXPORT void OpenMM_Platform_destroy(OpenMM_Platform* target) {
    delete reinterpret_cast<OpenMM::Platform*>(target);
}
OPENMM_EXPORT void OpenMM_Platform_registerPlatform(OpenMM_Platform* platform) {
    OpenMM::Platform::registerPlatform(reinterpret_cast<Platform *>(platform));
}
OPENMM_EXPORT int OpenMM_Platform_getNumPlatforms() {
    int result = OpenMM::Platform::getNumPlatforms();
    return result;
}
OPENMM_EXPORT OpenMM_Platform* OpenMM_Platform_getPlatform(int index) {
    Platform* result = &OpenMM::Platform::getPlatform(index);
    return reinterpret_cast<OpenMM_Platform*>(result);
}
OPENMM_EXPORT OpenMM_Platform* OpenMM_Platform_getPlatformByName(const char* name) {
    Platform* result = &OpenMM::Platform::getPlatformByName(std::string(name));
    return reinterpret_cast<OpenMM_Platform*>(result);
}
OPENMM_EXPORT OpenMM_Platform* OpenMM_Platform_findPlatform(const OpenMM_StringArray* kernelNames) {
    Platform* result = &OpenMM::Platform::findPlatform(*reinterpret_cast<const std::vector< std::string >*>(kernelNames));
    return reinterpret_cast<OpenMM_Platform*>(result);
}
OPENMM_EXPORT void OpenMM_Platform_loadPluginLibrary(const char* file) {
    OpenMM::Platform::loadPluginLibrary(std::string(file));
}
OPENMM_EXPORT const char* OpenMM_Platform_getDefaultPluginsDirectory() {
    const std::string* result = &OpenMM::Platform::getDefaultPluginsDirectory();
    return result->c_str();
}
OPENMM_EXPORT const char* OpenMM_Platform_getOpenMMVersion() {
    const std::string* result = &OpenMM::Platform::getOpenMMVersion();
    return result->c_str();
}
OPENMM_EXPORT const char* OpenMM_Platform_getName(const OpenMM_Platform* target) {
    const std::string* result = &reinterpret_cast<const OpenMM::Platform*>(target)->getName();
    return result->c_str();
}
OPENMM_EXPORT double OpenMM_Platform_getSpeed(const OpenMM_Platform* target) {
    double result = reinterpret_cast<const OpenMM::Platform*>(target)->getSpeed();
    return result;
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_Platform_supportsDoublePrecision(const OpenMM_Platform* target) {
    bool result = reinterpret_cast<const OpenMM::Platform*>(target)->supportsDoublePrecision();
    return (result ? OpenMM_True : OpenMM_False);
}
OPENMM_EXPORT const OpenMM_StringArray* OpenMM_Platform_getPropertyNames(const OpenMM_Platform* target) {
    const std::vector< std::string >* result = &reinterpret_cast<const OpenMM::Platform*>(target)->getPropertyNames();
    return reinterpret_cast<const OpenMM_StringArray*>(result);
}
OPENMM_EXPORT const char* OpenMM_Platform_getPropertyValue(const OpenMM_Platform* target, const OpenMM_Context* context, const char* property) {
    const std::string* result = &reinterpret_cast<const OpenMM::Platform*>(target)->getPropertyValue(*reinterpret_cast<const Context*>(context), std::string(property));
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_Platform_setPropertyValue(const OpenMM_Platform* target, OpenMM_Context* context, const char* property, const char* value) {
    reinterpret_cast<const OpenMM::Platform*>(target)->setPropertyValue(*reinterpret_cast<OpenMM::Context*>(context), std::string(property), std::string(value));
}
OPENMM_EXPORT const char* OpenMM_Platform_getPropertyDefaultValue(const OpenMM_Platform* target, const char* property) {
    const std::string* result = &reinterpret_cast<const OpenMM::Platform*>(target)->getPropertyDefaultValue(std::string(property));
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_Platform_setPropertyDefaultValue(OpenMM_Platform* target, const char* property, const char* value) {
    reinterpret_cast<OpenMM::Platform*>(target)->setPropertyDefaultValue(std::string(property), std::string(value));
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_Platform_supportsKernels(const OpenMM_Platform* target, const OpenMM_StringArray* kernelNames) {
    bool result = reinterpret_cast<const OpenMM::Platform*>(target)->supportsKernels(*reinterpret_cast<const std::vector< std::string >*>(kernelNames));
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::Continuous3DFunction */
OPENMM_EXPORT OpenMM_Continuous3DFunction* OpenMM_Continuous3DFunction_create(int xsize, int ysize, int zsize, const OpenMM_DoubleArray* values, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) {
    return reinterpret_cast<OpenMM_Continuous3DFunction*>(new OpenMM::Continuous3DFunction(xsize, ysize, zsize, *reinterpret_cast<const std::vector< double >*>(values), xmin, xmax, ymin, ymax, zmin, zmax));
}
OPENMM_EXPORT void OpenMM_Continuous3DFunction_destroy(OpenMM_Continuous3DFunction* target) {
    delete reinterpret_cast<OpenMM::Continuous3DFunction*>(target);
}
OPENMM_EXPORT void OpenMM_Continuous3DFunction_getFunctionParameters(const OpenMM_Continuous3DFunction* target, int* xsize, int* ysize, int* zsize, OpenMM_DoubleArray* values, double* xmin, double* xmax, double* ymin, double* ymax, double* zmin, double* zmax) {
    reinterpret_cast<const OpenMM::Continuous3DFunction*>(target)->getFunctionParameters(*reinterpret_cast<int*>(xsize), *reinterpret_cast<int*>(ysize), *reinterpret_cast<int*>(zsize), *reinterpret_cast<std::vector< double >*>(values), *reinterpret_cast<double*>(xmin), *reinterpret_cast<double*>(xmax), *reinterpret_cast<double*>(ymin), *reinterpret_cast<double*>(ymax), *reinterpret_cast<double*>(zmin), *reinterpret_cast<double*>(zmax));
}
OPENMM_EXPORT void OpenMM_Continuous3DFunction_setFunctionParameters(OpenMM_Continuous3DFunction* target, int xsize, int ysize, int zsize, const OpenMM_DoubleArray* values, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) {
    reinterpret_cast<OpenMM::Continuous3DFunction*>(target)->setFunctionParameters(xsize, ysize, zsize, *reinterpret_cast<const std::vector< double >*>(values), xmin, xmax, ymin, ymax, zmin, zmax);
}
OPENMM_EXPORT OpenMM_Continuous3DFunction* OpenMM_Continuous3DFunction_Copy(const OpenMM_Continuous3DFunction* target) {
    Continuous3DFunction * result = reinterpret_cast<const OpenMM::Continuous3DFunction*>(target)->Copy();
    return reinterpret_cast<OpenMM_Continuous3DFunction*>(result);
}

/* OpenMM::GBSAOBCForce */
OPENMM_EXPORT OpenMM_GBSAOBCForce* OpenMM_GBSAOBCForce_create() {
    return reinterpret_cast<OpenMM_GBSAOBCForce*>(new OpenMM::GBSAOBCForce());
}
OPENMM_EXPORT void OpenMM_GBSAOBCForce_destroy(OpenMM_GBSAOBCForce* target) {
    delete reinterpret_cast<OpenMM::GBSAOBCForce*>(target);
}
OPENMM_EXPORT int OpenMM_GBSAOBCForce_getNumParticles(const OpenMM_GBSAOBCForce* target) {
    int result = reinterpret_cast<const OpenMM::GBSAOBCForce*>(target)->getNumParticles();
    return result;
}
OPENMM_EXPORT int OpenMM_GBSAOBCForce_addParticle(OpenMM_GBSAOBCForce* target, double charge, double radius, double scalingFactor) {
    int result = reinterpret_cast<OpenMM::GBSAOBCForce*>(target)->addParticle(charge, radius, scalingFactor);
    return result;
}
OPENMM_EXPORT void OpenMM_GBSAOBCForce_getParticleParameters(const OpenMM_GBSAOBCForce* target, int index, double* charge, double* radius, double* scalingFactor) {
    reinterpret_cast<const OpenMM::GBSAOBCForce*>(target)->getParticleParameters(index, *reinterpret_cast<double*>(charge), *reinterpret_cast<double*>(radius), *reinterpret_cast<double*>(scalingFactor));
}
OPENMM_EXPORT void OpenMM_GBSAOBCForce_setParticleParameters(OpenMM_GBSAOBCForce* target, int index, double charge, double radius, double scalingFactor) {
    reinterpret_cast<OpenMM::GBSAOBCForce*>(target)->setParticleParameters(index, charge, radius, scalingFactor);
}
OPENMM_EXPORT double OpenMM_GBSAOBCForce_getSolventDielectric(const OpenMM_GBSAOBCForce* target) {
    double result = reinterpret_cast<const OpenMM::GBSAOBCForce*>(target)->getSolventDielectric();
    return result;
}
OPENMM_EXPORT void OpenMM_GBSAOBCForce_setSolventDielectric(OpenMM_GBSAOBCForce* target, double dielectric) {
    reinterpret_cast<OpenMM::GBSAOBCForce*>(target)->setSolventDielectric(dielectric);
}
OPENMM_EXPORT double OpenMM_GBSAOBCForce_getSoluteDielectric(const OpenMM_GBSAOBCForce* target) {
    double result = reinterpret_cast<const OpenMM::GBSAOBCForce*>(target)->getSoluteDielectric();
    return result;
}
OPENMM_EXPORT void OpenMM_GBSAOBCForce_setSoluteDielectric(OpenMM_GBSAOBCForce* target, double dielectric) {
    reinterpret_cast<OpenMM::GBSAOBCForce*>(target)->setSoluteDielectric(dielectric);
}
OPENMM_EXPORT double OpenMM_GBSAOBCForce_getSurfaceAreaEnergy(const OpenMM_GBSAOBCForce* target) {
    double result = reinterpret_cast<const OpenMM::GBSAOBCForce*>(target)->getSurfaceAreaEnergy();
    return result;
}
OPENMM_EXPORT void OpenMM_GBSAOBCForce_setSurfaceAreaEnergy(OpenMM_GBSAOBCForce* target, double energy) {
    reinterpret_cast<OpenMM::GBSAOBCForce*>(target)->setSurfaceAreaEnergy(energy);
}
OPENMM_EXPORT OpenMM_GBSAOBCForce_NonbondedMethod OpenMM_GBSAOBCForce_getNonbondedMethod(const OpenMM_GBSAOBCForce* target) {
    OpenMM::GBSAOBCForce::NonbondedMethod result = reinterpret_cast<const OpenMM::GBSAOBCForce*>(target)->getNonbondedMethod();
    return static_cast<OpenMM_GBSAOBCForce_NonbondedMethod>(result);
}
OPENMM_EXPORT void OpenMM_GBSAOBCForce_setNonbondedMethod(OpenMM_GBSAOBCForce* target, OpenMM_GBSAOBCForce_NonbondedMethod method) {
    reinterpret_cast<OpenMM::GBSAOBCForce*>(target)->setNonbondedMethod(static_cast<OpenMM::GBSAOBCForce::NonbondedMethod>(method));
}
OPENMM_EXPORT double OpenMM_GBSAOBCForce_getCutoffDistance(const OpenMM_GBSAOBCForce* target) {
    double result = reinterpret_cast<const OpenMM::GBSAOBCForce*>(target)->getCutoffDistance();
    return result;
}
OPENMM_EXPORT void OpenMM_GBSAOBCForce_setCutoffDistance(OpenMM_GBSAOBCForce* target, double distance) {
    reinterpret_cast<OpenMM::GBSAOBCForce*>(target)->setCutoffDistance(distance);
}
OPENMM_EXPORT void OpenMM_GBSAOBCForce_updateParametersInContext(OpenMM_GBSAOBCForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::GBSAOBCForce*>(target)->updateParametersInContext(*reinterpret_cast<OpenMM::Context*>(context));
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_GBSAOBCForce_usesPeriodicBoundaryConditions(const OpenMM_GBSAOBCForce* target) {
    bool result = reinterpret_cast<const OpenMM::GBSAOBCForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::CustomCentroidBondForce */
OPENMM_EXPORT OpenMM_CustomCentroidBondForce* OpenMM_CustomCentroidBondForce_create(int numGroups, const char* energy) {
    return reinterpret_cast<OpenMM_CustomCentroidBondForce*>(new OpenMM::CustomCentroidBondForce(numGroups, std::string(energy)));
}
OPENMM_EXPORT void OpenMM_CustomCentroidBondForce_destroy(OpenMM_CustomCentroidBondForce* target) {
    delete reinterpret_cast<OpenMM::CustomCentroidBondForce*>(target);
}
OPENMM_EXPORT int OpenMM_CustomCentroidBondForce_getNumGroupsPerBond(const OpenMM_CustomCentroidBondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomCentroidBondForce*>(target)->getNumGroupsPerBond();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomCentroidBondForce_getNumGroups(const OpenMM_CustomCentroidBondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomCentroidBondForce*>(target)->getNumGroups();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomCentroidBondForce_getNumBonds(const OpenMM_CustomCentroidBondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomCentroidBondForce*>(target)->getNumBonds();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomCentroidBondForce_getNumPerBondParameters(const OpenMM_CustomCentroidBondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomCentroidBondForce*>(target)->getNumPerBondParameters();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomCentroidBondForce_getNumGlobalParameters(const OpenMM_CustomCentroidBondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomCentroidBondForce*>(target)->getNumGlobalParameters();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomCentroidBondForce_getNumEnergyParameterDerivatives(const OpenMM_CustomCentroidBondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomCentroidBondForce*>(target)->getNumEnergyParameterDerivatives();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomCentroidBondForce_getNumTabulatedFunctions(const OpenMM_CustomCentroidBondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomCentroidBondForce*>(target)->getNumTabulatedFunctions();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomCentroidBondForce_getNumFunctions(const OpenMM_CustomCentroidBondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomCentroidBondForce*>(target)->getNumFunctions();
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomCentroidBondForce_getEnergyFunction(const OpenMM_CustomCentroidBondForce* target) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomCentroidBondForce*>(target)->getEnergyFunction();
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomCentroidBondForce_setEnergyFunction(OpenMM_CustomCentroidBondForce* target, const char* energy) {
    reinterpret_cast<OpenMM::CustomCentroidBondForce*>(target)->setEnergyFunction(std::string(energy));
}
OPENMM_EXPORT int OpenMM_CustomCentroidBondForce_addPerBondParameter(OpenMM_CustomCentroidBondForce* target, const char* name) {
    int result = reinterpret_cast<OpenMM::CustomCentroidBondForce*>(target)->addPerBondParameter(std::string(name));
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomCentroidBondForce_getPerBondParameterName(const OpenMM_CustomCentroidBondForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomCentroidBondForce*>(target)->getPerBondParameterName(index);
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomCentroidBondForce_setPerBondParameterName(OpenMM_CustomCentroidBondForce* target, int index, const char* name) {
    reinterpret_cast<OpenMM::CustomCentroidBondForce*>(target)->setPerBondParameterName(index, std::string(name));
}
OPENMM_EXPORT int OpenMM_CustomCentroidBondForce_addGlobalParameter(OpenMM_CustomCentroidBondForce* target, const char* name, double defaultValue) {
    int result = reinterpret_cast<OpenMM::CustomCentroidBondForce*>(target)->addGlobalParameter(std::string(name), defaultValue);
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomCentroidBondForce_getGlobalParameterName(const OpenMM_CustomCentroidBondForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomCentroidBondForce*>(target)->getGlobalParameterName(index);
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomCentroidBondForce_setGlobalParameterName(OpenMM_CustomCentroidBondForce* target, int index, const char* name) {
    reinterpret_cast<OpenMM::CustomCentroidBondForce*>(target)->setGlobalParameterName(index, std::string(name));
}
OPENMM_EXPORT double OpenMM_CustomCentroidBondForce_getGlobalParameterDefaultValue(const OpenMM_CustomCentroidBondForce* target, int index) {
    double result = reinterpret_cast<const OpenMM::CustomCentroidBondForce*>(target)->getGlobalParameterDefaultValue(index);
    return result;
}
OPENMM_EXPORT void OpenMM_CustomCentroidBondForce_setGlobalParameterDefaultValue(OpenMM_CustomCentroidBondForce* target, int index, double defaultValue) {
    reinterpret_cast<OpenMM::CustomCentroidBondForce*>(target)->setGlobalParameterDefaultValue(index, defaultValue);
}
OPENMM_EXPORT void OpenMM_CustomCentroidBondForce_addEnergyParameterDerivative(OpenMM_CustomCentroidBondForce* target, const char* name) {
    reinterpret_cast<OpenMM::CustomCentroidBondForce*>(target)->addEnergyParameterDerivative(std::string(name));
}
OPENMM_EXPORT const char* OpenMM_CustomCentroidBondForce_getEnergyParameterDerivativeName(const OpenMM_CustomCentroidBondForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomCentroidBondForce*>(target)->getEnergyParameterDerivativeName(index);
    return result->c_str();
}
OPENMM_EXPORT int OpenMM_CustomCentroidBondForce_addGroup(OpenMM_CustomCentroidBondForce* target, const OpenMM_IntArray* particles, const OpenMM_DoubleArray* weights) {
    int result = reinterpret_cast<OpenMM::CustomCentroidBondForce*>(target)->addGroup(*reinterpret_cast<const std::vector< int >*>(particles), *reinterpret_cast<const std::vector< double >*>(weights));
    return result;
}
OPENMM_EXPORT void OpenMM_CustomCentroidBondForce_getGroupParameters(const OpenMM_CustomCentroidBondForce* target, int index, OpenMM_IntArray* particles, OpenMM_DoubleArray* weights) {
    reinterpret_cast<const OpenMM::CustomCentroidBondForce*>(target)->getGroupParameters(index, *reinterpret_cast<std::vector< int >*>(particles), *reinterpret_cast<std::vector< double >*>(weights));
}
OPENMM_EXPORT void OpenMM_CustomCentroidBondForce_setGroupParameters(OpenMM_CustomCentroidBondForce* target, int index, const OpenMM_IntArray* particles, const OpenMM_DoubleArray* weights) {
    reinterpret_cast<OpenMM::CustomCentroidBondForce*>(target)->setGroupParameters(index, *reinterpret_cast<const std::vector< int >*>(particles), *reinterpret_cast<const std::vector< double >*>(weights));
}
OPENMM_EXPORT int OpenMM_CustomCentroidBondForce_addBond(OpenMM_CustomCentroidBondForce* target, const OpenMM_IntArray* groups, const OpenMM_DoubleArray* parameters) {
    int result = reinterpret_cast<OpenMM::CustomCentroidBondForce*>(target)->addBond(*reinterpret_cast<const std::vector< int >*>(groups), *reinterpret_cast<const std::vector< double >*>(parameters));
    return result;
}
OPENMM_EXPORT void OpenMM_CustomCentroidBondForce_getBondParameters(const OpenMM_CustomCentroidBondForce* target, int index, OpenMM_IntArray* groups, OpenMM_DoubleArray* parameters) {
    reinterpret_cast<const OpenMM::CustomCentroidBondForce*>(target)->getBondParameters(index, *reinterpret_cast<std::vector< int >*>(groups), *reinterpret_cast<std::vector< double >*>(parameters));
}
OPENMM_EXPORT void OpenMM_CustomCentroidBondForce_setBondParameters(OpenMM_CustomCentroidBondForce* target, int index, const OpenMM_IntArray* groups, const OpenMM_DoubleArray* parameters) {
    reinterpret_cast<OpenMM::CustomCentroidBondForce*>(target)->setBondParameters(index, *reinterpret_cast<const std::vector< int >*>(groups), *reinterpret_cast<const std::vector< double >*>(parameters));
}
OPENMM_EXPORT int OpenMM_CustomCentroidBondForce_addTabulatedFunction(OpenMM_CustomCentroidBondForce* target, const char* name, OpenMM_TabulatedFunction* function) {
    int result = reinterpret_cast<OpenMM::CustomCentroidBondForce*>(target)->addTabulatedFunction(std::string(name), reinterpret_cast<TabulatedFunction *>(function));
    return result;
}
OPENMM_EXPORT OpenMM_TabulatedFunction* OpenMM_CustomCentroidBondForce_getTabulatedFunction(OpenMM_CustomCentroidBondForce* target, int index) {
    TabulatedFunction* result = &reinterpret_cast<OpenMM::CustomCentroidBondForce*>(target)->getTabulatedFunction(index);
    return reinterpret_cast<OpenMM_TabulatedFunction*>(result);
}
OPENMM_EXPORT const char* OpenMM_CustomCentroidBondForce_getTabulatedFunctionName(const OpenMM_CustomCentroidBondForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomCentroidBondForce*>(target)->getTabulatedFunctionName(index);
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomCentroidBondForce_updateParametersInContext(OpenMM_CustomCentroidBondForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::CustomCentroidBondForce*>(target)->updateParametersInContext(*reinterpret_cast<OpenMM::Context*>(context));
}
OPENMM_EXPORT void OpenMM_CustomCentroidBondForce_setUsesPeriodicBoundaryConditions(OpenMM_CustomCentroidBondForce* target, OpenMM_Boolean periodic) {
    reinterpret_cast<OpenMM::CustomCentroidBondForce*>(target)->setUsesPeriodicBoundaryConditions(periodic);
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_CustomCentroidBondForce_usesPeriodicBoundaryConditions(const OpenMM_CustomCentroidBondForce* target) {
    bool result = reinterpret_cast<const OpenMM::CustomCentroidBondForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::LocalEnergyMinimizer */
OPENMM_EXPORT void OpenMM_LocalEnergyMinimizer_destroy(OpenMM_LocalEnergyMinimizer* target) {
    delete reinterpret_cast<OpenMM::LocalEnergyMinimizer*>(target);
}
OPENMM_EXPORT void OpenMM_LocalEnergyMinimizer_minimize(OpenMM_Context* context, double tolerance, int maxIterations) {
    OpenMM::LocalEnergyMinimizer::minimize(*reinterpret_cast<OpenMM::Context*>(context), tolerance, maxIterations);
}

/* OpenMM::TwoParticleAverageSite */
OPENMM_EXPORT OpenMM_TwoParticleAverageSite* OpenMM_TwoParticleAverageSite_create(int particle1, int particle2, double weight1, double weight2) {
    return reinterpret_cast<OpenMM_TwoParticleAverageSite*>(new OpenMM::TwoParticleAverageSite(particle1, particle2, weight1, weight2));
}
OPENMM_EXPORT void OpenMM_TwoParticleAverageSite_destroy(OpenMM_TwoParticleAverageSite* target) {
    delete reinterpret_cast<OpenMM::TwoParticleAverageSite*>(target);
}
OPENMM_EXPORT double OpenMM_TwoParticleAverageSite_getWeight(const OpenMM_TwoParticleAverageSite* target, int particle) {
    double result = reinterpret_cast<const OpenMM::TwoParticleAverageSite*>(target)->getWeight(particle);
    return result;
}

/* OpenMM::CustomCompoundBondForce */
OPENMM_EXPORT OpenMM_CustomCompoundBondForce* OpenMM_CustomCompoundBondForce_create(int numParticles, const char* energy) {
    return reinterpret_cast<OpenMM_CustomCompoundBondForce*>(new OpenMM::CustomCompoundBondForce(numParticles, std::string(energy)));
}
OPENMM_EXPORT void OpenMM_CustomCompoundBondForce_destroy(OpenMM_CustomCompoundBondForce* target) {
    delete reinterpret_cast<OpenMM::CustomCompoundBondForce*>(target);
}
OPENMM_EXPORT int OpenMM_CustomCompoundBondForce_getNumParticlesPerBond(const OpenMM_CustomCompoundBondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomCompoundBondForce*>(target)->getNumParticlesPerBond();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomCompoundBondForce_getNumBonds(const OpenMM_CustomCompoundBondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomCompoundBondForce*>(target)->getNumBonds();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomCompoundBondForce_getNumPerBondParameters(const OpenMM_CustomCompoundBondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomCompoundBondForce*>(target)->getNumPerBondParameters();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomCompoundBondForce_getNumGlobalParameters(const OpenMM_CustomCompoundBondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomCompoundBondForce*>(target)->getNumGlobalParameters();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomCompoundBondForce_getNumEnergyParameterDerivatives(const OpenMM_CustomCompoundBondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomCompoundBondForce*>(target)->getNumEnergyParameterDerivatives();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomCompoundBondForce_getNumTabulatedFunctions(const OpenMM_CustomCompoundBondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomCompoundBondForce*>(target)->getNumTabulatedFunctions();
    return result;
}
OPENMM_EXPORT int OpenMM_CustomCompoundBondForce_getNumFunctions(const OpenMM_CustomCompoundBondForce* target) {
    int result = reinterpret_cast<const OpenMM::CustomCompoundBondForce*>(target)->getNumFunctions();
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomCompoundBondForce_getEnergyFunction(const OpenMM_CustomCompoundBondForce* target) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomCompoundBondForce*>(target)->getEnergyFunction();
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomCompoundBondForce_setEnergyFunction(OpenMM_CustomCompoundBondForce* target, const char* energy) {
    reinterpret_cast<OpenMM::CustomCompoundBondForce*>(target)->setEnergyFunction(std::string(energy));
}
OPENMM_EXPORT int OpenMM_CustomCompoundBondForce_addPerBondParameter(OpenMM_CustomCompoundBondForce* target, const char* name) {
    int result = reinterpret_cast<OpenMM::CustomCompoundBondForce*>(target)->addPerBondParameter(std::string(name));
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomCompoundBondForce_getPerBondParameterName(const OpenMM_CustomCompoundBondForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomCompoundBondForce*>(target)->getPerBondParameterName(index);
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomCompoundBondForce_setPerBondParameterName(OpenMM_CustomCompoundBondForce* target, int index, const char* name) {
    reinterpret_cast<OpenMM::CustomCompoundBondForce*>(target)->setPerBondParameterName(index, std::string(name));
}
OPENMM_EXPORT int OpenMM_CustomCompoundBondForce_addGlobalParameter(OpenMM_CustomCompoundBondForce* target, const char* name, double defaultValue) {
    int result = reinterpret_cast<OpenMM::CustomCompoundBondForce*>(target)->addGlobalParameter(std::string(name), defaultValue);
    return result;
}
OPENMM_EXPORT const char* OpenMM_CustomCompoundBondForce_getGlobalParameterName(const OpenMM_CustomCompoundBondForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomCompoundBondForce*>(target)->getGlobalParameterName(index);
    return result->c_str();
}
OPENMM_EXPORT void OpenMM_CustomCompoundBondForce_setGlobalParameterName(OpenMM_CustomCompoundBondForce* target, int index, const char* name) {
    reinterpret_cast<OpenMM::CustomCompoundBondForce*>(target)->setGlobalParameterName(index, std::string(name));
}
OPENMM_EXPORT double OpenMM_CustomCompoundBondForce_getGlobalParameterDefaultValue(const OpenMM_CustomCompoundBondForce* target, int index) {
    double result = reinterpret_cast<const OpenMM::CustomCompoundBondForce*>(target)->getGlobalParameterDefaultValue(index);
    return result;
}
OPENMM_EXPORT void OpenMM_CustomCompoundBondForce_setGlobalParameterDefaultValue(OpenMM_CustomCompoundBondForce* target, int index, double defaultValue) {
    reinterpret_cast<OpenMM::CustomCompoundBondForce*>(target)->setGlobalParameterDefaultValue(index, defaultValue);
}
OPENMM_EXPORT void OpenMM_CustomCompoundBondForce_addEnergyParameterDerivative(OpenMM_CustomCompoundBondForce* target, const char* name) {
    reinterpret_cast<OpenMM::CustomCompoundBondForce*>(target)->addEnergyParameterDerivative(std::string(name));
}
OPENMM_EXPORT const char* OpenMM_CustomCompoundBondForce_getEnergyParameterDerivativeName(const OpenMM_CustomCompoundBondForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomCompoundBondForce*>(target)->getEnergyParameterDerivativeName(index);
    return result->c_str();
}
OPENMM_EXPORT int OpenMM_CustomCompoundBondForce_addBond(OpenMM_CustomCompoundBondForce* target, const OpenMM_IntArray* particles, const OpenMM_DoubleArray* parameters) {
    int result = reinterpret_cast<OpenMM::CustomCompoundBondForce*>(target)->addBond(*reinterpret_cast<const std::vector< int >*>(particles), *reinterpret_cast<const std::vector< double >*>(parameters));
    return result;
}
OPENMM_EXPORT void OpenMM_CustomCompoundBondForce_getBondParameters(const OpenMM_CustomCompoundBondForce* target, int index, OpenMM_IntArray* particles, OpenMM_DoubleArray* parameters) {
    reinterpret_cast<const OpenMM::CustomCompoundBondForce*>(target)->getBondParameters(index, *reinterpret_cast<std::vector< int >*>(particles), *reinterpret_cast<std::vector< double >*>(parameters));
}
OPENMM_EXPORT void OpenMM_CustomCompoundBondForce_setBondParameters(OpenMM_CustomCompoundBondForce* target, int index, const OpenMM_IntArray* particles, const OpenMM_DoubleArray* parameters) {
    reinterpret_cast<OpenMM::CustomCompoundBondForce*>(target)->setBondParameters(index, *reinterpret_cast<const std::vector< int >*>(particles), *reinterpret_cast<const std::vector< double >*>(parameters));
}
OPENMM_EXPORT int OpenMM_CustomCompoundBondForce_addTabulatedFunction(OpenMM_CustomCompoundBondForce* target, const char* name, OpenMM_TabulatedFunction* function) {
    int result = reinterpret_cast<OpenMM::CustomCompoundBondForce*>(target)->addTabulatedFunction(std::string(name), reinterpret_cast<TabulatedFunction *>(function));
    return result;
}
OPENMM_EXPORT OpenMM_TabulatedFunction* OpenMM_CustomCompoundBondForce_getTabulatedFunction(OpenMM_CustomCompoundBondForce* target, int index) {
    TabulatedFunction* result = &reinterpret_cast<OpenMM::CustomCompoundBondForce*>(target)->getTabulatedFunction(index);
    return reinterpret_cast<OpenMM_TabulatedFunction*>(result);
}
OPENMM_EXPORT const char* OpenMM_CustomCompoundBondForce_getTabulatedFunctionName(const OpenMM_CustomCompoundBondForce* target, int index) {
    const std::string* result = &reinterpret_cast<const OpenMM::CustomCompoundBondForce*>(target)->getTabulatedFunctionName(index);
    return result->c_str();
}
OPENMM_EXPORT int OpenMM_CustomCompoundBondForce_addFunction(OpenMM_CustomCompoundBondForce* target, const char* name, const OpenMM_DoubleArray* values, double min, double max) {
    int result = reinterpret_cast<OpenMM::CustomCompoundBondForce*>(target)->addFunction(std::string(name), *reinterpret_cast<const std::vector< double >*>(values), min, max);
    return result;
}
OPENMM_EXPORT void OpenMM_CustomCompoundBondForce_getFunctionParameters(const OpenMM_CustomCompoundBondForce* target, int index, char** name, OpenMM_DoubleArray* values, double* min, double* max) {
    reinterpret_cast<const OpenMM::CustomCompoundBondForce*>(target)->getFunctionParameters(index, *reinterpret_cast<std::string*>(name), *reinterpret_cast<std::vector< double >*>(values), *reinterpret_cast<double*>(min), *reinterpret_cast<double*>(max));
}
OPENMM_EXPORT void OpenMM_CustomCompoundBondForce_setFunctionParameters(OpenMM_CustomCompoundBondForce* target, int index, const char* name, const OpenMM_DoubleArray* values, double min, double max) {
    reinterpret_cast<OpenMM::CustomCompoundBondForce*>(target)->setFunctionParameters(index, std::string(name), *reinterpret_cast<const std::vector< double >*>(values), min, max);
}
OPENMM_EXPORT void OpenMM_CustomCompoundBondForce_updateParametersInContext(OpenMM_CustomCompoundBondForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::CustomCompoundBondForce*>(target)->updateParametersInContext(*reinterpret_cast<OpenMM::Context*>(context));
}
OPENMM_EXPORT void OpenMM_CustomCompoundBondForce_setUsesPeriodicBoundaryConditions(OpenMM_CustomCompoundBondForce* target, OpenMM_Boolean periodic) {
    reinterpret_cast<OpenMM::CustomCompoundBondForce*>(target)->setUsesPeriodicBoundaryConditions(periodic);
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_CustomCompoundBondForce_usesPeriodicBoundaryConditions(const OpenMM_CustomCompoundBondForce* target) {
    bool result = reinterpret_cast<const OpenMM::CustomCompoundBondForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::RMSDForce */
OPENMM_EXPORT OpenMM_RMSDForce* OpenMM_RMSDForce_create(const OpenMM_Vec3Array* referencePositions, const OpenMM_IntArray* particles) {
    return reinterpret_cast<OpenMM_RMSDForce*>(new OpenMM::RMSDForce(*reinterpret_cast<const std::vector< Vec3 >*>(referencePositions), *reinterpret_cast<const std::vector< int >*>(particles)));
}
OPENMM_EXPORT void OpenMM_RMSDForce_destroy(OpenMM_RMSDForce* target) {
    delete reinterpret_cast<OpenMM::RMSDForce*>(target);
}
OPENMM_EXPORT const OpenMM_Vec3Array* OpenMM_RMSDForce_getReferencePositions(const OpenMM_RMSDForce* target) {
    const std::vector< Vec3 >* result = &reinterpret_cast<const OpenMM::RMSDForce*>(target)->getReferencePositions();
    return reinterpret_cast<const OpenMM_Vec3Array*>(result);
}
OPENMM_EXPORT void OpenMM_RMSDForce_setReferencePositions(OpenMM_RMSDForce* target, const OpenMM_Vec3Array* positions) {
    reinterpret_cast<OpenMM::RMSDForce*>(target)->setReferencePositions(*reinterpret_cast<const std::vector< Vec3 >*>(positions));
}
OPENMM_EXPORT const OpenMM_IntArray* OpenMM_RMSDForce_getParticles(const OpenMM_RMSDForce* target) {
    const std::vector< int >* result = &reinterpret_cast<const OpenMM::RMSDForce*>(target)->getParticles();
    return reinterpret_cast<const OpenMM_IntArray*>(result);
}
OPENMM_EXPORT void OpenMM_RMSDForce_setParticles(OpenMM_RMSDForce* target, const OpenMM_IntArray* particles) {
    reinterpret_cast<OpenMM::RMSDForce*>(target)->setParticles(*reinterpret_cast<const std::vector< int >*>(particles));
}
OPENMM_EXPORT void OpenMM_RMSDForce_updateParametersInContext(OpenMM_RMSDForce* target, OpenMM_Context* context) {
    reinterpret_cast<OpenMM::RMSDForce*>(target)->updateParametersInContext(*reinterpret_cast<OpenMM::Context*>(context));
}
OPENMM_EXPORT OpenMM_Boolean OpenMM_RMSDForce_usesPeriodicBoundaryConditions(const OpenMM_RMSDForce* target) {
    bool result = reinterpret_cast<const OpenMM::RMSDForce*>(target)->usesPeriodicBoundaryConditions();
    return (result ? OpenMM_True : OpenMM_False);
}

/* OpenMM::BrownianIntegrator */
OPENMM_EXPORT OpenMM_BrownianIntegrator* OpenMM_BrownianIntegrator_create(double temperature, double frictionCoeff, double stepSize) {
    return reinterpret_cast<OpenMM_BrownianIntegrator*>(new OpenMM::BrownianIntegrator(temperature, frictionCoeff, stepSize));
}
OPENMM_EXPORT void OpenMM_BrownianIntegrator_destroy(OpenMM_BrownianIntegrator* target) {
    delete reinterpret_cast<OpenMM::BrownianIntegrator*>(target);
}
OPENMM_EXPORT double OpenMM_BrownianIntegrator_getTemperature(const OpenMM_BrownianIntegrator* target) {
    double result = reinterpret_cast<const OpenMM::BrownianIntegrator*>(target)->getTemperature();
    return result;
}
OPENMM_EXPORT void OpenMM_BrownianIntegrator_setTemperature(OpenMM_BrownianIntegrator* target, double temp) {
    reinterpret_cast<OpenMM::BrownianIntegrator*>(target)->setTemperature(temp);
}
OPENMM_EXPORT double OpenMM_BrownianIntegrator_getFriction(const OpenMM_BrownianIntegrator* target) {
    double result = reinterpret_cast<const OpenMM::BrownianIntegrator*>(target)->getFriction();
    return result;
}
OPENMM_EXPORT void OpenMM_BrownianIntegrator_setFriction(OpenMM_BrownianIntegrator* target, double coeff) {
    reinterpret_cast<OpenMM::BrownianIntegrator*>(target)->setFriction(coeff);
}
OPENMM_EXPORT int OpenMM_BrownianIntegrator_getRandomNumberSeed(const OpenMM_BrownianIntegrator* target) {
    int result = reinterpret_cast<const OpenMM::BrownianIntegrator*>(target)->getRandomNumberSeed();
    return result;
}
OPENMM_EXPORT void OpenMM_BrownianIntegrator_setRandomNumberSeed(OpenMM_BrownianIntegrator* target, int seed) {
    reinterpret_cast<OpenMM::BrownianIntegrator*>(target)->setRandomNumberSeed(seed);
}
OPENMM_EXPORT void OpenMM_BrownianIntegrator_step(OpenMM_BrownianIntegrator* target, int steps) {
    reinterpret_cast<OpenMM::BrownianIntegrator*>(target)->step(steps);
}

/* OpenMM::VariableLangevinIntegrator */
OPENMM_EXPORT OpenMM_VariableLangevinIntegrator* OpenMM_VariableLangevinIntegrator_create(double temperature, double frictionCoeff, double errorTol) {
    return reinterpret_cast<OpenMM_VariableLangevinIntegrator*>(new OpenMM::VariableLangevinIntegrator(temperature, frictionCoeff, errorTol));
}
OPENMM_EXPORT void OpenMM_VariableLangevinIntegrator_destroy(OpenMM_VariableLangevinIntegrator* target) {
    delete reinterpret_cast<OpenMM::VariableLangevinIntegrator*>(target);
}
OPENMM_EXPORT double OpenMM_VariableLangevinIntegrator_getTemperature(const OpenMM_VariableLangevinIntegrator* target) {
    double result = reinterpret_cast<const OpenMM::VariableLangevinIntegrator*>(target)->getTemperature();
    return result;
}
OPENMM_EXPORT void OpenMM_VariableLangevinIntegrator_setTemperature(OpenMM_VariableLangevinIntegrator* target, double temp) {
    reinterpret_cast<OpenMM::VariableLangevinIntegrator*>(target)->setTemperature(temp);
}
OPENMM_EXPORT double OpenMM_VariableLangevinIntegrator_getFriction(const OpenMM_VariableLangevinIntegrator* target) {
    double result = reinterpret_cast<const OpenMM::VariableLangevinIntegrator*>(target)->getFriction();
    return result;
}
OPENMM_EXPORT void OpenMM_VariableLangevinIntegrator_setFriction(OpenMM_VariableLangevinIntegrator* target, double coeff) {
    reinterpret_cast<OpenMM::VariableLangevinIntegrator*>(target)->setFriction(coeff);
}
OPENMM_EXPORT double OpenMM_VariableLangevinIntegrator_getErrorTolerance(const OpenMM_VariableLangevinIntegrator* target) {
    double result = reinterpret_cast<const OpenMM::VariableLangevinIntegrator*>(target)->getErrorTolerance();
    return result;
}
OPENMM_EXPORT void OpenMM_VariableLangevinIntegrator_setErrorTolerance(OpenMM_VariableLangevinIntegrator* target, double tol) {
    reinterpret_cast<OpenMM::VariableLangevinIntegrator*>(target)->setErrorTolerance(tol);
}
OPENMM_EXPORT double OpenMM_VariableLangevinIntegrator_getMaximumStepSize(const OpenMM_VariableLangevinIntegrator* target) {
    double result = reinterpret_cast<const OpenMM::VariableLangevinIntegrator*>(target)->getMaximumStepSize();
    return result;
}
OPENMM_EXPORT void OpenMM_VariableLangevinIntegrator_setMaximumStepSize(OpenMM_VariableLangevinIntegrator* target, double size) {
    reinterpret_cast<OpenMM::VariableLangevinIntegrator*>(target)->setMaximumStepSize(size);
}
OPENMM_EXPORT int OpenMM_VariableLangevinIntegrator_getRandomNumberSeed(const OpenMM_VariableLangevinIntegrator* target) {
    int result = reinterpret_cast<const OpenMM::VariableLangevinIntegrator*>(target)->getRandomNumberSeed();
    return result;
}
OPENMM_EXPORT void OpenMM_VariableLangevinIntegrator_setRandomNumberSeed(OpenMM_VariableLangevinIntegrator* target, int seed) {
    reinterpret_cast<OpenMM::VariableLangevinIntegrator*>(target)->setRandomNumberSeed(seed);
}
OPENMM_EXPORT void OpenMM_VariableLangevinIntegrator_step(OpenMM_VariableLangevinIntegrator* target, int steps) {
    reinterpret_cast<OpenMM::VariableLangevinIntegrator*>(target)->step(steps);
}
OPENMM_EXPORT void OpenMM_VariableLangevinIntegrator_stepTo(OpenMM_VariableLangevinIntegrator* target, double time) {
    reinterpret_cast<OpenMM::VariableLangevinIntegrator*>(target)->stepTo(time);
}

}

