
#include "OpenMM.h"
#include "OpenMMCWrapper.h"
#include <cstring>
#include <vector>
#include <cstdlib>

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

static void convertStringToChars(char* source, char*& cstr, int& length) {
	length = strlen(source);
	cstr = new char[length+1];
	strcpy(cstr, source);
    free(source);
}

extern "C" {

/* OpenMM_Vec3 */
OPENMM_EXPORT void openmm_vec3_scale_(const OpenMM_Vec3& vec, double const& scale, OpenMM_Vec3& result) {
    result = OpenMM_Vec3_scale(vec, scale);
}
OPENMM_EXPORT void OPENMM_VEC3_SCALE(const OpenMM_Vec3& vec, double const& scale, OpenMM_Vec3& result) {
    result = OpenMM_Vec3_scale(vec, scale);
}

/* OpenMM_Vec3Array */
OPENMM_EXPORT void openmm_vec3array_create_(OpenMM_Vec3Array*& result, const int& size) {
    result = OpenMM_Vec3Array_create(size);
}
OPENMM_EXPORT void OPENMM_VEC3ARRAY_CREATE(OpenMM_Vec3Array*& result, const int& size) {
    result = OpenMM_Vec3Array_create(size);
}
OPENMM_EXPORT void openmm_vec3array_destroy_(OpenMM_Vec3Array*& array) {
    OpenMM_Vec3Array_destroy(array);
    array = 0;
}
OPENMM_EXPORT void OPENMM_VEC3ARRAY_DESTROY(OpenMM_Vec3Array*& array) {
    OpenMM_Vec3Array_destroy(array);
    array = 0;
}
OPENMM_EXPORT int openmm_vec3array_getsize_(const OpenMM_Vec3Array* const& array) {
    return OpenMM_Vec3Array_getSize(array);
}
OPENMM_EXPORT int OPENMM_VEC3ARRAY_GETSIZE(const OpenMM_Vec3Array* const& array) {
    return OpenMM_Vec3Array_getSize(array);
}
OPENMM_EXPORT void openmm_vec3array_resize_(OpenMM_Vec3Array* const& array, const int& size) {
    OpenMM_Vec3Array_resize(array, size);
}
OPENMM_EXPORT void OPENMM_VEC3ARRAY_RESIZE(OpenMM_Vec3Array* const& array, const int& size) {
    OpenMM_Vec3Array_resize(array, size);
}
OPENMM_EXPORT void openmm_vec3array_append_(OpenMM_Vec3Array* const& array, const OpenMM_Vec3& vec) {
    OpenMM_Vec3Array_append(array, vec);
}
OPENMM_EXPORT void OPENMM_VEC3ARRAY_APPEND(OpenMM_Vec3Array* const& array, const OpenMM_Vec3& vec) {
    OpenMM_Vec3Array_append(array, vec);
}
OPENMM_EXPORT void openmm_vec3array_set_(OpenMM_Vec3Array* const& array, const int& index, const OpenMM_Vec3& vec) {
    OpenMM_Vec3Array_set(array, index-1, vec);
}
OPENMM_EXPORT void OPENMM_VEC3ARRAY_SET(OpenMM_Vec3Array* const& array, const int& index, const OpenMM_Vec3& vec) {
    OpenMM_Vec3Array_set(array, index-1, vec);
}
OPENMM_EXPORT void openmm_vec3array_get_(const OpenMM_Vec3Array* const& array, const int& index, OpenMM_Vec3& result) {
    result = *OpenMM_Vec3Array_get(array, index-1);
}
OPENMM_EXPORT void OPENMM_VEC3ARRAY_GET(const OpenMM_Vec3Array* const& array, const int& index, OpenMM_Vec3& result) {
    result = *OpenMM_Vec3Array_get(array, index-1);
}

/* OpenMM_StringArray */
OPENMM_EXPORT void openmm_stringarray_create_(OpenMM_StringArray*& result, const int& size) {
    result = OpenMM_StringArray_create(size);
}
OPENMM_EXPORT void OPENMM_STRINGARRAY_CREATE(OpenMM_StringArray*& result, const int& size) {
    result = OpenMM_StringArray_create(size);
}
OPENMM_EXPORT void openmm_stringarray_destroy_(OpenMM_StringArray*& array) {
    OpenMM_StringArray_destroy(array);
    array = 0;
}
OPENMM_EXPORT void OPENMM_STRINGARRAY_DESTROY(OpenMM_StringArray*& array) {
    OpenMM_StringArray_destroy(array);
    array = 0;
}
OPENMM_EXPORT int openmm_stringarray_getsize_(const OpenMM_StringArray* const& array) {
    return OpenMM_StringArray_getSize(array);
}
OPENMM_EXPORT int OPENMM_STRINGARRAY_GETSIZE(const OpenMM_StringArray* const& array) {
    return OpenMM_StringArray_getSize(array);
}
OPENMM_EXPORT void openmm_stringarray_resize_(OpenMM_StringArray* const& array, const int& size) {
    OpenMM_StringArray_resize(array, size);
}
OPENMM_EXPORT void OPENMM_STRINGARRAY_RESIZE(OpenMM_StringArray* const& array, const int& size) {
    OpenMM_StringArray_resize(array, size);
}
OPENMM_EXPORT void openmm_stringarray_append_(OpenMM_StringArray* const& array, const char* str, int length) {
    OpenMM_StringArray_append(array, makeString(str, length).c_str());
}
OPENMM_EXPORT void OPENMM_STRINGARRAY_APPEND(OpenMM_StringArray* const& array, const char* str, int length) {
    OpenMM_StringArray_append(array, makeString(str, length).c_str());
}
OPENMM_EXPORT void openmm_stringarray_set_(OpenMM_StringArray* const& array, const int& index, const char* str, int length) {
  OpenMM_StringArray_set(array, index-1, makeString(str, length).c_str());
  }
OPENMM_EXPORT void OPENMM_STRINGARRAY_SET(OpenMM_StringArray* const& array, const int& index, const char* str, int length) {
  OpenMM_StringArray_set(array, index-1, makeString(str, length).c_str());
}
OPENMM_EXPORT void openmm_stringarray_get_(const OpenMM_StringArray* const& array, const int& index, char* result, int length) {
    const char* str = OpenMM_StringArray_get(array, index-1);
    copyAndPadString(result, str, length);
}
OPENMM_EXPORT void OPENMM_STRINGARRAY_GET(const OpenMM_StringArray* const& array, const int& index, char* result, int length) {
    const char* str = OpenMM_StringArray_get(array, index-1);
    copyAndPadString(result, str, length);
}

/* OpenMM_BondArray */
OPENMM_EXPORT void openmm_bondarray_create_(OpenMM_BondArray*& result, const int& size) {
    result = OpenMM_BondArray_create(size);
}
OPENMM_EXPORT void OPENMM_BONDARRAY_CREATE(OpenMM_BondArray*& result, const int& size) {
    result = OpenMM_BondArray_create(size);
}
OPENMM_EXPORT void openmm_bondarray_destroy_(OpenMM_BondArray*& array) {
    OpenMM_BondArray_destroy(array);
    array = 0;
}
OPENMM_EXPORT void OPENMM_BONDARRAY_DESTROY(OpenMM_BondArray*& array) {
    OpenMM_BondArray_destroy(array);
    array = 0;
}
OPENMM_EXPORT int openmm_bondarray_getsize_(const OpenMM_BondArray* const& array) {
    return OpenMM_BondArray_getSize(array);
}
OPENMM_EXPORT int OPENMM_BONDARRAY_GETSIZE(const OpenMM_BondArray* const& array) {
    return OpenMM_BondArray_getSize(array);
}
OPENMM_EXPORT void openmm_bondarray_resize_(OpenMM_BondArray* const& array, const int& size) {
    OpenMM_BondArray_resize(array, size);
}
OPENMM_EXPORT void OPENMM_BONDARRAY_RESIZE(OpenMM_BondArray* const& array, const int& size) {
    OpenMM_BondArray_resize(array, size);
}
OPENMM_EXPORT void openmm_bondarray_append_(OpenMM_BondArray* const& array, const int& particle1, const int& particle2) {
    OpenMM_BondArray_append(array, particle1, particle2);
}
OPENMM_EXPORT void OPENMM_BONDARRAY_APPEND(OpenMM_BondArray* const& array, const int& particle1, const int& particle2) {
    OpenMM_BondArray_append(array, particle1, particle2);
}
OPENMM_EXPORT void openmm_bondarray_set_(OpenMM_BondArray* const& array, const int& index, const int& particle1, const int& particle2) {
    OpenMM_BondArray_set(array, index-1, particle1, particle2);
}
OPENMM_EXPORT void OPENMM_BONDARRAY_SET(OpenMM_BondArray* const& array, const int& index, const int& particle1, const int& particle2) {
    OpenMM_BondArray_set(array, index-1, particle1, particle2);
}
OPENMM_EXPORT void openmm_bondarray_get_(const OpenMM_BondArray* const& array, const int& index, int* particle1, int* particle2) {
    OpenMM_BondArray_get(array, index-1, particle1, particle2);
}
OPENMM_EXPORT void OPENMM_BONDARRAY_GET(const OpenMM_BondArray* const& array, const int& index, int* particle1, int* particle2) {
    OpenMM_BondArray_get(array, index-1, particle1, particle2);
}

/* OpenMM_ParameterArray */
OPENMM_EXPORT int openmm_parameterarray_getsize_(const OpenMM_ParameterArray* const& array) {
    return OpenMM_ParameterArray_getSize(array);
}
OPENMM_EXPORT int OPENMM_PARAMETERARRAY_GETSIZE(const OpenMM_ParameterArray* const& array) {
    return OpenMM_ParameterArray_getSize(array);
}
OPENMM_EXPORT double openmm_parameterarray_get_(const OpenMM_ParameterArray* const& array, const char* name, int length) {
    return OpenMM_ParameterArray_get(array, makeString(name, length).c_str());
}
OPENMM_EXPORT double OPENMM_PARAMETERARRAY_GET(const OpenMM_ParameterArray* const& array, const char* name, int length) {
    return OpenMM_ParameterArray_get(array, makeString(name, length).c_str());
}

/* OpenMM_PropertyArray */
OPENMM_EXPORT int openmm_propertyarray_getsize_(const OpenMM_PropertyArray* const& array) {
    return OpenMM_PropertyArray_getSize(array);
}
OPENMM_EXPORT int OPENMM_PROPERTYARRAY_GETSIZE(const OpenMM_PropertyArray* const& array) {
    return OpenMM_PropertyArray_getSize(array);
}
OPENMM_EXPORT const char* openmm_propertyarray_get_(const OpenMM_PropertyArray* const& array, const char* name, int length) {
    return OpenMM_PropertyArray_get(array, makeString(name, length).c_str());
}
OPENMM_EXPORT const char* OPENMM_PROPERTYARRAY_GET(const OpenMM_PropertyArray* const& array, const char* name, int length) {
    return OpenMM_PropertyArray_get(array, makeString(name, length).c_str());
}

/* OpenMM_DoubleArray */
OPENMM_EXPORT void openmm_doublearray_create_(OpenMM_DoubleArray*& result, const int& size) {
    result = OpenMM_DoubleArray_create(size);
}
OPENMM_EXPORT void OPENMM_DOUBLEARRAY_CREATE(OpenMM_DoubleArray*& result, const int& size) {
    result = OpenMM_DoubleArray_create(size);
}
OPENMM_EXPORT void openmm_doublearray_destroy_(OpenMM_DoubleArray*& array) {
    OpenMM_DoubleArray_destroy(array);
    array = 0;
}
OPENMM_EXPORT void OPENMM_DOUBLEARRAY_DESTROY(OpenMM_DoubleArray*& array) {
    OpenMM_DoubleArray_destroy(array);
    array = 0;
}
OPENMM_EXPORT int openmm_doublearray_getsize_(const OpenMM_DoubleArray* const& array) {
    return OpenMM_DoubleArray_getSize(array);
}
OPENMM_EXPORT int OPENMM_DOUBLEARRAY_GETSIZE(const OpenMM_DoubleArray* const& array) {
    return OpenMM_DoubleArray_getSize(array);
}
OPENMM_EXPORT void openmm_doublearray_resize_(OpenMM_DoubleArray* const& array, const int& size) {
    OpenMM_DoubleArray_resize(array, size);
}
OPENMM_EXPORT void OPENMM_DOUBLEARRAY_RESIZE(OpenMM_DoubleArray* const& array, const int& size) {
    OpenMM_DoubleArray_resize(array, size);
}
OPENMM_EXPORT void openmm_doublearray_append_(OpenMM_DoubleArray* const& array, const double& value) {
    OpenMM_DoubleArray_append(array, value);
}
OPENMM_EXPORT void OPENMM_DOUBLEARRAY_APPEND(OpenMM_DoubleArray* const& array, const double& value) {
    OpenMM_DoubleArray_append(array, value);
}
OPENMM_EXPORT void openmm_doublearray_set_(OpenMM_DoubleArray* const& array, const int& index, const double& value) {
    OpenMM_DoubleArray_set(array, index-1, value);
}
OPENMM_EXPORT void OPENMM_DOUBLEARRAY_SET(OpenMM_DoubleArray* const& array, const int& index, const double& value) {
    OpenMM_DoubleArray_set(array, index-1, value);
}
OPENMM_EXPORT void openmm_doublearray_get_(const OpenMM_DoubleArray* const& array, const int& index, double& result) {
    result = OpenMM_DoubleArray_get(array, index-1);
}
OPENMM_EXPORT void OPENMM_DOUBLEARRAY_GET(const OpenMM_DoubleArray* const& array, const int& index, double& result) {
    result = OpenMM_DoubleArray_get(array, index-1);
}

/* OpenMM_IntArray */
OPENMM_EXPORT void openmm_intarray_create_(OpenMM_IntArray*& result, const int& size) {
    result = OpenMM_IntArray_create(size);
}
OPENMM_EXPORT void OPENMM_INTARRAY_CREATE(OpenMM_IntArray*& result, const int& size) {
    result = OpenMM_IntArray_create(size);
}
OPENMM_EXPORT void openmm_intarray_destroy_(OpenMM_IntArray*& array) {
    OpenMM_IntArray_destroy(array);
    array = 0;
}
OPENMM_EXPORT void OPENMM_INTARRAY_DESTROY(OpenMM_IntArray*& array) {
    OpenMM_IntArray_destroy(array);
    array = 0;
}
OPENMM_EXPORT int openmm_intarray_getsize_(const OpenMM_IntArray* const& array) {
    return OpenMM_IntArray_getSize(array);
}
OPENMM_EXPORT int OPENMM_INTARRAY_GETSIZE(const OpenMM_IntArray* const& array) {
    return OpenMM_IntArray_getSize(array);
}
OPENMM_EXPORT void openmm_intarray_resize_(OpenMM_IntArray* const& array, const int& size) {
    OpenMM_IntArray_resize(array, size);
}
OPENMM_EXPORT void OPENMM_INTARRAY_RESIZE(OpenMM_IntArray* const& array, const int& size) {
    OpenMM_IntArray_resize(array, size);
}
OPENMM_EXPORT void openmm_intarray_append_(OpenMM_IntArray* const& array, const int& value) {
    OpenMM_IntArray_append(array, value);
}
OPENMM_EXPORT void OPENMM_INTARRAY_APPEND(OpenMM_IntArray* const& array, const int& value) {
    OpenMM_IntArray_append(array, value);
}
OPENMM_EXPORT void openmm_intarray_set_(OpenMM_IntArray* const& array, const int& index, const int& value) {
    OpenMM_IntArray_set(array, index-1, value);
}
OPENMM_EXPORT void OPENMM_INTARRAY_SET(OpenMM_IntArray* const& array, const int& index, const int& value) {
    OpenMM_IntArray_set(array, index-1, value);
}
OPENMM_EXPORT void openmm_intarray_get_(const OpenMM_IntArray* const& array, const int& index, int& result) {
    result = OpenMM_IntArray_get(array, index-1);
}
OPENMM_EXPORT void OPENMM_INTARRAY_GET(const OpenMM_IntArray* const& array, const int& index, int& result) {
    result = OpenMM_IntArray_get(array, index-1);
}

/* OpenMM_IntSet */
OPENMM_EXPORT void openmm_intset_create_(OpenMM_IntSet*& result) {
    result = OpenMM_IntSet_create();
}
OPENMM_EXPORT void OPENMM_INTSET_CREATE(OpenMM_IntSet*& result) {
    result = OpenMM_IntSet_create();
}
OPENMM_EXPORT void openmm_intset_destroy_(OpenMM_IntSet*& array) {
    OpenMM_IntSet_destroy(array);
    array = 0;
}
OPENMM_EXPORT void OPENMM_INTSET_DESTROY(OpenMM_IntSet*& array) {
    OpenMM_IntSet_destroy(array);
    array = 0;
}
OPENMM_EXPORT int openmm_intset_getsize_(const OpenMM_IntSet* const& array) {
    return OpenMM_IntSet_getSize(array);
}
OPENMM_EXPORT int OPENMM_INTSET_GETSIZE(const OpenMM_IntSet* const& array) {
    return OpenMM_IntSet_getSize(array);
}
OPENMM_EXPORT void openmm_intset_insert_(OpenMM_IntSet* const& array, const int& value) {
    OpenMM_IntSet_insert(array, value);
}
OPENMM_EXPORT void OPENMM_INTSET_INSERT(OpenMM_IntSet* const& array, const int& value) {
    OpenMM_IntSet_insert(array, value);
}

/* These methods need to be handled specially, since their C++ APIs cannot be directly translated to C.
   Unlike the C++ versions, the return value is allocated on the heap, and you must delete it yourself. */
OPENMM_EXPORT void openmm_context_getstate_(const OpenMM_Context*& target, int const& types, int const& enforcePeriodicBox, OpenMM_State*& result) {
    result = OpenMM_Context_getState(target, types, enforcePeriodicBox);
}
OPENMM_EXPORT void OPENMM_CONTEXT_GETSTATE(const OpenMM_Context*& target, int const& types, int const& enforcePeriodicBox, OpenMM_State*& result) {
    result = OpenMM_Context_getState(target, types, enforcePeriodicBox);
}
OPENMM_EXPORT void openmm_context_getstate_2_(const OpenMM_Context*& target, int const& types, int const& enforcePeriodicBox, int const& groups, OpenMM_State*& result) {
    result = OpenMM_Context_getState_2(target, types, enforcePeriodicBox, groups);
}
OPENMM_EXPORT void OPENMM_CONTEXT_GETSTATE_2(const OpenMM_Context*& target, int const& types, int const& enforcePeriodicBox, int const& groups, OpenMM_State*& result) {
    result = OpenMM_Context_getState_2(target, types, enforcePeriodicBox, groups);
}
OPENMM_EXPORT void openmm_platform_loadpluginsfromdirectory_(const char* directory, OpenMM_StringArray*& result, int length) {
    result = OpenMM_Platform_loadPluginsFromDirectory(makeString(directory, length).c_str());
}
OPENMM_EXPORT void OPENMM_PLATFORM_LOADPLUGINSFROMDIRECTORY(const char* directory, OpenMM_StringArray*& result, int length) {
    result = OpenMM_Platform_loadPluginsFromDirectory(makeString(directory, length).c_str());
}
OPENMM_EXPORT void openmm_platform_getpluginloadfailures_(OpenMM_StringArray*& result) {
    result = OpenMM_Platform_getPluginLoadFailures();
}
OPENMM_EXPORT void OPENMM_PLATFORM_GETPLUGINLOADFAILURES(OpenMM_StringArray*& result) {
    result = OpenMM_Platform_getPluginLoadFailures();
}
OPENMM_EXPORT void openmm_xmlserializer_serializesystemtoc_(OpenMM_System*& system, char*& result, int& result_length) {
    convertStringToChars(OpenMM_XmlSerializer_serializeSystem(system), result, result_length);
}
OPENMM_EXPORT void OPENMM_XMLSERIALIZER_SERIALIZESYSTEMTOC(OpenMM_System*& system, char*& result, int& result_length) {
    convertStringToChars(OpenMM_XmlSerializer_serializeSystem(system), result, result_length);
}
OPENMM_EXPORT void openmm_xmlserializer_serializestatetoc_(OpenMM_State*& state, char*& result, int& result_length) {
    convertStringToChars(OpenMM_XmlSerializer_serializeState(state), result, result_length);
}
OPENMM_EXPORT void OPENMM_XMLSERIALIZER_SERIALIZESTATETOC(OpenMM_State*& state, char*& result, int& result_length) {
    convertStringToChars(OpenMM_XmlSerializer_serializeState(state), result, result_length);
}
OPENMM_EXPORT void openmm_xmlserializer_serializeintegratortoc_(OpenMM_Integrator*& integrator, char*& result, int& result_length) {
    convertStringToChars(OpenMM_XmlSerializer_serializeIntegrator(integrator), result, result_length);
}
OPENMM_EXPORT void OPENMM_XMLSERIALIZER_SERIALIZEINTEGRATORTOC(OpenMM_Integrator*& integrator, char*& result, int& result_length) {
    convertStringToChars(OpenMM_XmlSerializer_serializeIntegrator(integrator), result, result_length);
}
OPENMM_EXPORT void openmm_xmlserializer_deserializesystem_(const char* xml, OpenMM_System*& result, int length) {
    result = OpenMM_XmlSerializer_deserializeSystem(makeString(xml, length).c_str());
}
OPENMM_EXPORT void OPENMM_XMLSERIALIZER_DESERIALIZESYSTEM(const char* xml, OpenMM_System*& result, int length) {
    result = OpenMM_XmlSerializer_deserializeSystem(makeString(xml, length).c_str());
}
OPENMM_EXPORT void openmm_xmlserializer_deserializestate_(const char* xml, OpenMM_State*& result, int length) {
    result = OpenMM_XmlSerializer_deserializeState(makeString(xml, length).c_str());
}
OPENMM_EXPORT void OPENMM_XMLSERIALIZER_DESERIALIZESTATE(const char* xml, OpenMM_State*& result, int length) {
    result = OpenMM_XmlSerializer_deserializeState(makeString(xml, length).c_str());
}
OPENMM_EXPORT void openmm_xmlserializer_deserializeintegrator_(const char* xml, OpenMM_Integrator*& result, int length) {
    result = OpenMM_XmlSerializer_deserializeIntegrator(makeString(xml, length).c_str());
}
OPENMM_EXPORT void OPENMM_XMLSERIALIZER_DESERIALIZEINTEGRATOR(const char* xml, OpenMM_Integrator*& result, int length) {
    result = OpenMM_XmlSerializer_deserializeIntegrator(makeString(xml, length).c_str());
}

/* OpenMM::Force */
OPENMM_EXPORT void openmm_force_destroy_(OpenMM_Force*& destroy) {
    OpenMM_Force_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_FORCE_DESTROY(OpenMM_Force*& destroy) {
    OpenMM_Force_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT int openmm_force_getforcegroup_(const OpenMM_Force*& target) {
    return OpenMM_Force_getForceGroup(target);
}
OPENMM_EXPORT int OPENMM_FORCE_GETFORCEGROUP(const OpenMM_Force*& target) {
    return OpenMM_Force_getForceGroup(target);
}
OPENMM_EXPORT void openmm_force_setforcegroup_(OpenMM_Force*& target, int const& group) {
    OpenMM_Force_setForceGroup(target, group);
}
OPENMM_EXPORT void OPENMM_FORCE_SETFORCEGROUP(OpenMM_Force*& target, int const& group) {
    OpenMM_Force_setForceGroup(target, group);
}
OPENMM_EXPORT void openmm_force_usesperiodicboundaryconditions_(const OpenMM_Force*& target, OpenMM_Boolean& result) {
    result = OpenMM_Force_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_FORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_Force*& target, OpenMM_Boolean& result) {
    result = OpenMM_Force_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::CustomAngleForce */
OPENMM_EXPORT void openmm_customangleforce_create_(OpenMM_CustomAngleForce*& result, const char* energy, int energy_length) {
    result = OpenMM_CustomAngleForce_create(makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMANGLEFORCE_CREATE(OpenMM_CustomAngleForce*& result, const char* energy, int energy_length) {
    result = OpenMM_CustomAngleForce_create(makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void openmm_customangleforce_destroy_(OpenMM_CustomAngleForce*& destroy) {
    OpenMM_CustomAngleForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_CUSTOMANGLEFORCE_DESTROY(OpenMM_CustomAngleForce*& destroy) {
    OpenMM_CustomAngleForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT int openmm_customangleforce_getnumangles_(const OpenMM_CustomAngleForce*& target) {
    return OpenMM_CustomAngleForce_getNumAngles(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMANGLEFORCE_GETNUMANGLES(const OpenMM_CustomAngleForce*& target) {
    return OpenMM_CustomAngleForce_getNumAngles(target);
}
OPENMM_EXPORT int openmm_customangleforce_getnumperangleparameters_(const OpenMM_CustomAngleForce*& target) {
    return OpenMM_CustomAngleForce_getNumPerAngleParameters(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMANGLEFORCE_GETNUMPERANGLEPARAMETERS(const OpenMM_CustomAngleForce*& target) {
    return OpenMM_CustomAngleForce_getNumPerAngleParameters(target);
}
OPENMM_EXPORT int openmm_customangleforce_getnumglobalparameters_(const OpenMM_CustomAngleForce*& target) {
    return OpenMM_CustomAngleForce_getNumGlobalParameters(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMANGLEFORCE_GETNUMGLOBALPARAMETERS(const OpenMM_CustomAngleForce*& target) {
    return OpenMM_CustomAngleForce_getNumGlobalParameters(target);
}
OPENMM_EXPORT int openmm_customangleforce_getnumenergyparameterderivatives_(const OpenMM_CustomAngleForce*& target) {
    return OpenMM_CustomAngleForce_getNumEnergyParameterDerivatives(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMANGLEFORCE_GETNUMENERGYPARAMETERDERIVATIVES(const OpenMM_CustomAngleForce*& target) {
    return OpenMM_CustomAngleForce_getNumEnergyParameterDerivatives(target);
}
OPENMM_EXPORT void openmm_customangleforce_getenergyfunction_(const OpenMM_CustomAngleForce*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomAngleForce_getEnergyFunction(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMANGLEFORCE_GETENERGYFUNCTION(const OpenMM_CustomAngleForce*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomAngleForce_getEnergyFunction(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customangleforce_setenergyfunction_(OpenMM_CustomAngleForce*& target, const char* energy, int energy_length) {
    OpenMM_CustomAngleForce_setEnergyFunction(target, makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMANGLEFORCE_SETENERGYFUNCTION(OpenMM_CustomAngleForce*& target, const char* energy, int energy_length) {
    OpenMM_CustomAngleForce_setEnergyFunction(target, makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT int openmm_customangleforce_addperangleparameter_(OpenMM_CustomAngleForce*& target, const char* name, int name_length) {
    return OpenMM_CustomAngleForce_addPerAngleParameter(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT int OPENMM_CUSTOMANGLEFORCE_ADDPERANGLEPARAMETER(OpenMM_CustomAngleForce*& target, const char* name, int name_length) {
    return OpenMM_CustomAngleForce_addPerAngleParameter(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void openmm_customangleforce_getperangleparametername_(const OpenMM_CustomAngleForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomAngleForce_getPerAngleParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMANGLEFORCE_GETPERANGLEPARAMETERNAME(const OpenMM_CustomAngleForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomAngleForce_getPerAngleParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customangleforce_setperangleparametername_(OpenMM_CustomAngleForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomAngleForce_setPerAngleParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMANGLEFORCE_SETPERANGLEPARAMETERNAME(OpenMM_CustomAngleForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomAngleForce_setPerAngleParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT int openmm_customangleforce_addglobalparameter_(OpenMM_CustomAngleForce*& target, const char* name, double const& defaultValue, int name_length) {
    return OpenMM_CustomAngleForce_addGlobalParameter(target, makeString(name, name_length).c_str(), defaultValue);
}
OPENMM_EXPORT int OPENMM_CUSTOMANGLEFORCE_ADDGLOBALPARAMETER(OpenMM_CustomAngleForce*& target, const char* name, double const& defaultValue, int name_length) {
    return OpenMM_CustomAngleForce_addGlobalParameter(target, makeString(name, name_length).c_str(), defaultValue);
}
OPENMM_EXPORT void openmm_customangleforce_getglobalparametername_(const OpenMM_CustomAngleForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomAngleForce_getGlobalParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMANGLEFORCE_GETGLOBALPARAMETERNAME(const OpenMM_CustomAngleForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomAngleForce_getGlobalParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customangleforce_setglobalparametername_(OpenMM_CustomAngleForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomAngleForce_setGlobalParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMANGLEFORCE_SETGLOBALPARAMETERNAME(OpenMM_CustomAngleForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomAngleForce_setGlobalParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT double openmm_customangleforce_getglobalparameterdefaultvalue_(const OpenMM_CustomAngleForce*& target, int const& index) {
    return OpenMM_CustomAngleForce_getGlobalParameterDefaultValue(target, index);
}
OPENMM_EXPORT double OPENMM_CUSTOMANGLEFORCE_GETGLOBALPARAMETERDEFAULTVALUE(const OpenMM_CustomAngleForce*& target, int const& index) {
    return OpenMM_CustomAngleForce_getGlobalParameterDefaultValue(target, index);
}
OPENMM_EXPORT void openmm_customangleforce_setglobalparameterdefaultvalue_(OpenMM_CustomAngleForce*& target, int const& index, double const& defaultValue) {
    OpenMM_CustomAngleForce_setGlobalParameterDefaultValue(target, index, defaultValue);
}
OPENMM_EXPORT void OPENMM_CUSTOMANGLEFORCE_SETGLOBALPARAMETERDEFAULTVALUE(OpenMM_CustomAngleForce*& target, int const& index, double const& defaultValue) {
    OpenMM_CustomAngleForce_setGlobalParameterDefaultValue(target, index, defaultValue);
}
OPENMM_EXPORT void openmm_customangleforce_addenergyparameterderivative_(OpenMM_CustomAngleForce*& target, const char* name, int name_length) {
    OpenMM_CustomAngleForce_addEnergyParameterDerivative(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMANGLEFORCE_ADDENERGYPARAMETERDERIVATIVE(OpenMM_CustomAngleForce*& target, const char* name, int name_length) {
    OpenMM_CustomAngleForce_addEnergyParameterDerivative(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void openmm_customangleforce_getenergyparameterderivativename_(const OpenMM_CustomAngleForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomAngleForce_getEnergyParameterDerivativeName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMANGLEFORCE_GETENERGYPARAMETERDERIVATIVENAME(const OpenMM_CustomAngleForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomAngleForce_getEnergyParameterDerivativeName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT int openmm_customangleforce_addangle_(OpenMM_CustomAngleForce*& target, int const& particle1, int const& particle2, int const& particle3, const OpenMM_DoubleArray*& parameters) {
    return OpenMM_CustomAngleForce_addAngle(target, particle1, particle2, particle3, parameters);
}
OPENMM_EXPORT int OPENMM_CUSTOMANGLEFORCE_ADDANGLE(OpenMM_CustomAngleForce*& target, int const& particle1, int const& particle2, int const& particle3, const OpenMM_DoubleArray*& parameters) {
    return OpenMM_CustomAngleForce_addAngle(target, particle1, particle2, particle3, parameters);
}
OPENMM_EXPORT void openmm_customangleforce_getangleparameters_(const OpenMM_CustomAngleForce*& target, int const& index, int* particle1, int* particle2, int* particle3, OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomAngleForce_getAngleParameters(target, index, particle1, particle2, particle3, parameters);
}
OPENMM_EXPORT void OPENMM_CUSTOMANGLEFORCE_GETANGLEPARAMETERS(const OpenMM_CustomAngleForce*& target, int const& index, int* particle1, int* particle2, int* particle3, OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomAngleForce_getAngleParameters(target, index, particle1, particle2, particle3, parameters);
}
OPENMM_EXPORT void openmm_customangleforce_setangleparameters_(OpenMM_CustomAngleForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, const OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomAngleForce_setAngleParameters(target, index, particle1, particle2, particle3, parameters);
}
OPENMM_EXPORT void OPENMM_CUSTOMANGLEFORCE_SETANGLEPARAMETERS(OpenMM_CustomAngleForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, const OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomAngleForce_setAngleParameters(target, index, particle1, particle2, particle3, parameters);
}
OPENMM_EXPORT void openmm_customangleforce_updateparametersincontext_(OpenMM_CustomAngleForce*& target, OpenMM_Context*& context) {
    OpenMM_CustomAngleForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void OPENMM_CUSTOMANGLEFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_CustomAngleForce*& target, OpenMM_Context*& context) {
    OpenMM_CustomAngleForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void openmm_customangleforce_setusesperiodicboundaryconditions_(OpenMM_CustomAngleForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_CustomAngleForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT void OPENMM_CUSTOMANGLEFORCE_SETUSESPERIODICBOUNDARYCONDITIONS(OpenMM_CustomAngleForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_CustomAngleForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT void openmm_customangleforce_usesperiodicboundaryconditions_(const OpenMM_CustomAngleForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CustomAngleForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_CUSTOMANGLEFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_CustomAngleForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CustomAngleForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::MonteCarloBarostat */
OPENMM_EXPORT void openmm_montecarlobarostat_create_(OpenMM_MonteCarloBarostat*& result, double const& defaultPressure, double const& defaultTemperature, int const& frequency) {
    result = OpenMM_MonteCarloBarostat_create(defaultPressure, defaultTemperature, frequency);
}
OPENMM_EXPORT void OPENMM_MONTECARLOBAROSTAT_CREATE(OpenMM_MonteCarloBarostat*& result, double const& defaultPressure, double const& defaultTemperature, int const& frequency) {
    result = OpenMM_MonteCarloBarostat_create(defaultPressure, defaultTemperature, frequency);
}
OPENMM_EXPORT void openmm_montecarlobarostat_destroy_(OpenMM_MonteCarloBarostat*& destroy) {
    OpenMM_MonteCarloBarostat_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_MONTECARLOBAROSTAT_DESTROY(OpenMM_MonteCarloBarostat*& destroy) {
    OpenMM_MonteCarloBarostat_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void openmm_montecarlobarostat_pressure_(char* result, int result_length) {
    const char* result_chars = OpenMM_MonteCarloBarostat_Pressure();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_MONTECARLOBAROSTAT_PRESSURE(char* result, int result_length) {
    const char* result_chars = OpenMM_MonteCarloBarostat_Pressure();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_montecarlobarostat_temperature_(char* result, int result_length) {
    const char* result_chars = OpenMM_MonteCarloBarostat_Temperature();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_MONTECARLOBAROSTAT_TEMPERATURE(char* result, int result_length) {
    const char* result_chars = OpenMM_MonteCarloBarostat_Temperature();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT double openmm_montecarlobarostat_getdefaultpressure_(const OpenMM_MonteCarloBarostat*& target) {
    return OpenMM_MonteCarloBarostat_getDefaultPressure(target);
}
OPENMM_EXPORT double OPENMM_MONTECARLOBAROSTAT_GETDEFAULTPRESSURE(const OpenMM_MonteCarloBarostat*& target) {
    return OpenMM_MonteCarloBarostat_getDefaultPressure(target);
}
OPENMM_EXPORT void openmm_montecarlobarostat_setdefaultpressure_(OpenMM_MonteCarloBarostat*& target, double const& pressure) {
    OpenMM_MonteCarloBarostat_setDefaultPressure(target, pressure);
}
OPENMM_EXPORT void OPENMM_MONTECARLOBAROSTAT_SETDEFAULTPRESSURE(OpenMM_MonteCarloBarostat*& target, double const& pressure) {
    OpenMM_MonteCarloBarostat_setDefaultPressure(target, pressure);
}
OPENMM_EXPORT int openmm_montecarlobarostat_getfrequency_(const OpenMM_MonteCarloBarostat*& target) {
    return OpenMM_MonteCarloBarostat_getFrequency(target);
}
OPENMM_EXPORT int OPENMM_MONTECARLOBAROSTAT_GETFREQUENCY(const OpenMM_MonteCarloBarostat*& target) {
    return OpenMM_MonteCarloBarostat_getFrequency(target);
}
OPENMM_EXPORT void openmm_montecarlobarostat_setfrequency_(OpenMM_MonteCarloBarostat*& target, int const& freq) {
    OpenMM_MonteCarloBarostat_setFrequency(target, freq);
}
OPENMM_EXPORT void OPENMM_MONTECARLOBAROSTAT_SETFREQUENCY(OpenMM_MonteCarloBarostat*& target, int const& freq) {
    OpenMM_MonteCarloBarostat_setFrequency(target, freq);
}
OPENMM_EXPORT double openmm_montecarlobarostat_getdefaulttemperature_(const OpenMM_MonteCarloBarostat*& target) {
    return OpenMM_MonteCarloBarostat_getDefaultTemperature(target);
}
OPENMM_EXPORT double OPENMM_MONTECARLOBAROSTAT_GETDEFAULTTEMPERATURE(const OpenMM_MonteCarloBarostat*& target) {
    return OpenMM_MonteCarloBarostat_getDefaultTemperature(target);
}
OPENMM_EXPORT void openmm_montecarlobarostat_setdefaulttemperature_(OpenMM_MonteCarloBarostat*& target, double const& temp) {
    OpenMM_MonteCarloBarostat_setDefaultTemperature(target, temp);
}
OPENMM_EXPORT void OPENMM_MONTECARLOBAROSTAT_SETDEFAULTTEMPERATURE(OpenMM_MonteCarloBarostat*& target, double const& temp) {
    OpenMM_MonteCarloBarostat_setDefaultTemperature(target, temp);
}
OPENMM_EXPORT int openmm_montecarlobarostat_getrandomnumberseed_(const OpenMM_MonteCarloBarostat*& target) {
    return OpenMM_MonteCarloBarostat_getRandomNumberSeed(target);
}
OPENMM_EXPORT int OPENMM_MONTECARLOBAROSTAT_GETRANDOMNUMBERSEED(const OpenMM_MonteCarloBarostat*& target) {
    return OpenMM_MonteCarloBarostat_getRandomNumberSeed(target);
}
OPENMM_EXPORT void openmm_montecarlobarostat_setrandomnumberseed_(OpenMM_MonteCarloBarostat*& target, int const& seed) {
    OpenMM_MonteCarloBarostat_setRandomNumberSeed(target, seed);
}
OPENMM_EXPORT void OPENMM_MONTECARLOBAROSTAT_SETRANDOMNUMBERSEED(OpenMM_MonteCarloBarostat*& target, int const& seed) {
    OpenMM_MonteCarloBarostat_setRandomNumberSeed(target, seed);
}
OPENMM_EXPORT void openmm_montecarlobarostat_usesperiodicboundaryconditions_(const OpenMM_MonteCarloBarostat*& target, OpenMM_Boolean& result) {
    result = OpenMM_MonteCarloBarostat_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_MONTECARLOBAROSTAT_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_MonteCarloBarostat*& target, OpenMM_Boolean& result) {
    result = OpenMM_MonteCarloBarostat_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::RBTorsionForce */
OPENMM_EXPORT void openmm_rbtorsionforce_create_(OpenMM_RBTorsionForce*& result) {
    result = OpenMM_RBTorsionForce_create();
}
OPENMM_EXPORT void OPENMM_RBTORSIONFORCE_CREATE(OpenMM_RBTorsionForce*& result) {
    result = OpenMM_RBTorsionForce_create();
}
OPENMM_EXPORT void openmm_rbtorsionforce_destroy_(OpenMM_RBTorsionForce*& destroy) {
    OpenMM_RBTorsionForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_RBTORSIONFORCE_DESTROY(OpenMM_RBTorsionForce*& destroy) {
    OpenMM_RBTorsionForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT int openmm_rbtorsionforce_getnumtorsions_(const OpenMM_RBTorsionForce*& target) {
    return OpenMM_RBTorsionForce_getNumTorsions(target);
}
OPENMM_EXPORT int OPENMM_RBTORSIONFORCE_GETNUMTORSIONS(const OpenMM_RBTorsionForce*& target) {
    return OpenMM_RBTorsionForce_getNumTorsions(target);
}
OPENMM_EXPORT int openmm_rbtorsionforce_addtorsion_(OpenMM_RBTorsionForce*& target, int const& particle1, int const& particle2, int const& particle3, int const& particle4, double const& c0, double const& c1, double const& c2, double const& c3, double const& c4, double const& c5) {
    return OpenMM_RBTorsionForce_addTorsion(target, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
}
OPENMM_EXPORT int OPENMM_RBTORSIONFORCE_ADDTORSION(OpenMM_RBTorsionForce*& target, int const& particle1, int const& particle2, int const& particle3, int const& particle4, double const& c0, double const& c1, double const& c2, double const& c3, double const& c4, double const& c5) {
    return OpenMM_RBTorsionForce_addTorsion(target, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
}
OPENMM_EXPORT void openmm_rbtorsionforce_gettorsionparameters_(const OpenMM_RBTorsionForce*& target, int const& index, int* particle1, int* particle2, int* particle3, int* particle4, double* c0, double* c1, double* c2, double* c3, double* c4, double* c5) {
    OpenMM_RBTorsionForce_getTorsionParameters(target, index, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
}
OPENMM_EXPORT void OPENMM_RBTORSIONFORCE_GETTORSIONPARAMETERS(const OpenMM_RBTorsionForce*& target, int const& index, int* particle1, int* particle2, int* particle3, int* particle4, double* c0, double* c1, double* c2, double* c3, double* c4, double* c5) {
    OpenMM_RBTorsionForce_getTorsionParameters(target, index, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
}
OPENMM_EXPORT void openmm_rbtorsionforce_settorsionparameters_(OpenMM_RBTorsionForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, int const& particle4, double const& c0, double const& c1, double const& c2, double const& c3, double const& c4, double const& c5) {
    OpenMM_RBTorsionForce_setTorsionParameters(target, index, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
}
OPENMM_EXPORT void OPENMM_RBTORSIONFORCE_SETTORSIONPARAMETERS(OpenMM_RBTorsionForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, int const& particle4, double const& c0, double const& c1, double const& c2, double const& c3, double const& c4, double const& c5) {
    OpenMM_RBTorsionForce_setTorsionParameters(target, index, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
}
OPENMM_EXPORT void openmm_rbtorsionforce_updateparametersincontext_(OpenMM_RBTorsionForce*& target, OpenMM_Context*& context) {
    OpenMM_RBTorsionForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void OPENMM_RBTORSIONFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_RBTorsionForce*& target, OpenMM_Context*& context) {
    OpenMM_RBTorsionForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void openmm_rbtorsionforce_setusesperiodicboundaryconditions_(OpenMM_RBTorsionForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_RBTorsionForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT void OPENMM_RBTORSIONFORCE_SETUSESPERIODICBOUNDARYCONDITIONS(OpenMM_RBTorsionForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_RBTorsionForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT void openmm_rbtorsionforce_usesperiodicboundaryconditions_(const OpenMM_RBTorsionForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_RBTorsionForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_RBTORSIONFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_RBTorsionForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_RBTorsionForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::Integrator */
OPENMM_EXPORT void openmm_integrator_destroy_(OpenMM_Integrator*& destroy) {
    OpenMM_Integrator_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_INTEGRATOR_DESTROY(OpenMM_Integrator*& destroy) {
    OpenMM_Integrator_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT double openmm_integrator_getstepsize_(const OpenMM_Integrator*& target) {
    return OpenMM_Integrator_getStepSize(target);
}
OPENMM_EXPORT double OPENMM_INTEGRATOR_GETSTEPSIZE(const OpenMM_Integrator*& target) {
    return OpenMM_Integrator_getStepSize(target);
}
OPENMM_EXPORT void openmm_integrator_setstepsize_(OpenMM_Integrator*& target, double const& size) {
    OpenMM_Integrator_setStepSize(target, size);
}
OPENMM_EXPORT void OPENMM_INTEGRATOR_SETSTEPSIZE(OpenMM_Integrator*& target, double const& size) {
    OpenMM_Integrator_setStepSize(target, size);
}
OPENMM_EXPORT double openmm_integrator_getconstrainttolerance_(const OpenMM_Integrator*& target) {
    return OpenMM_Integrator_getConstraintTolerance(target);
}
OPENMM_EXPORT double OPENMM_INTEGRATOR_GETCONSTRAINTTOLERANCE(const OpenMM_Integrator*& target) {
    return OpenMM_Integrator_getConstraintTolerance(target);
}
OPENMM_EXPORT void openmm_integrator_setconstrainttolerance_(OpenMM_Integrator*& target, double const& tol) {
    OpenMM_Integrator_setConstraintTolerance(target, tol);
}
OPENMM_EXPORT void OPENMM_INTEGRATOR_SETCONSTRAINTTOLERANCE(OpenMM_Integrator*& target, double const& tol) {
    OpenMM_Integrator_setConstraintTolerance(target, tol);
}
OPENMM_EXPORT void openmm_integrator_step_(OpenMM_Integrator*& target, int const& steps) {
    OpenMM_Integrator_step(target, steps);
}
OPENMM_EXPORT void OPENMM_INTEGRATOR_STEP(OpenMM_Integrator*& target, int const& steps) {
    OpenMM_Integrator_step(target, steps);
}

/* OpenMM::VerletIntegrator */
OPENMM_EXPORT void openmm_verletintegrator_create_(OpenMM_VerletIntegrator*& result, double const& stepSize) {
    result = OpenMM_VerletIntegrator_create(stepSize);
}
OPENMM_EXPORT void OPENMM_VERLETINTEGRATOR_CREATE(OpenMM_VerletIntegrator*& result, double const& stepSize) {
    result = OpenMM_VerletIntegrator_create(stepSize);
}
OPENMM_EXPORT void openmm_verletintegrator_destroy_(OpenMM_VerletIntegrator*& destroy) {
    OpenMM_VerletIntegrator_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_VERLETINTEGRATOR_DESTROY(OpenMM_VerletIntegrator*& destroy) {
    OpenMM_VerletIntegrator_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void openmm_verletintegrator_step_(OpenMM_VerletIntegrator*& target, int const& steps) {
    OpenMM_VerletIntegrator_step(target, steps);
}
OPENMM_EXPORT void OPENMM_VERLETINTEGRATOR_STEP(OpenMM_VerletIntegrator*& target, int const& steps) {
    OpenMM_VerletIntegrator_step(target, steps);
}

/* OpenMM::LangevinIntegrator */
OPENMM_EXPORT void openmm_langevinintegrator_create_(OpenMM_LangevinIntegrator*& result, double const& temperature, double const& frictionCoeff, double const& stepSize) {
    result = OpenMM_LangevinIntegrator_create(temperature, frictionCoeff, stepSize);
}
OPENMM_EXPORT void OPENMM_LANGEVININTEGRATOR_CREATE(OpenMM_LangevinIntegrator*& result, double const& temperature, double const& frictionCoeff, double const& stepSize) {
    result = OpenMM_LangevinIntegrator_create(temperature, frictionCoeff, stepSize);
}
OPENMM_EXPORT void openmm_langevinintegrator_destroy_(OpenMM_LangevinIntegrator*& destroy) {
    OpenMM_LangevinIntegrator_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_LANGEVININTEGRATOR_DESTROY(OpenMM_LangevinIntegrator*& destroy) {
    OpenMM_LangevinIntegrator_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT double openmm_langevinintegrator_gettemperature_(const OpenMM_LangevinIntegrator*& target) {
    return OpenMM_LangevinIntegrator_getTemperature(target);
}
OPENMM_EXPORT double OPENMM_LANGEVININTEGRATOR_GETTEMPERATURE(const OpenMM_LangevinIntegrator*& target) {
    return OpenMM_LangevinIntegrator_getTemperature(target);
}
OPENMM_EXPORT void openmm_langevinintegrator_settemperature_(OpenMM_LangevinIntegrator*& target, double const& temp) {
    OpenMM_LangevinIntegrator_setTemperature(target, temp);
}
OPENMM_EXPORT void OPENMM_LANGEVININTEGRATOR_SETTEMPERATURE(OpenMM_LangevinIntegrator*& target, double const& temp) {
    OpenMM_LangevinIntegrator_setTemperature(target, temp);
}
OPENMM_EXPORT double openmm_langevinintegrator_getfriction_(const OpenMM_LangevinIntegrator*& target) {
    return OpenMM_LangevinIntegrator_getFriction(target);
}
OPENMM_EXPORT double OPENMM_LANGEVININTEGRATOR_GETFRICTION(const OpenMM_LangevinIntegrator*& target) {
    return OpenMM_LangevinIntegrator_getFriction(target);
}
OPENMM_EXPORT void openmm_langevinintegrator_setfriction_(OpenMM_LangevinIntegrator*& target, double const& coeff) {
    OpenMM_LangevinIntegrator_setFriction(target, coeff);
}
OPENMM_EXPORT void OPENMM_LANGEVININTEGRATOR_SETFRICTION(OpenMM_LangevinIntegrator*& target, double const& coeff) {
    OpenMM_LangevinIntegrator_setFriction(target, coeff);
}
OPENMM_EXPORT int openmm_langevinintegrator_getrandomnumberseed_(const OpenMM_LangevinIntegrator*& target) {
    return OpenMM_LangevinIntegrator_getRandomNumberSeed(target);
}
OPENMM_EXPORT int OPENMM_LANGEVININTEGRATOR_GETRANDOMNUMBERSEED(const OpenMM_LangevinIntegrator*& target) {
    return OpenMM_LangevinIntegrator_getRandomNumberSeed(target);
}
OPENMM_EXPORT void openmm_langevinintegrator_setrandomnumberseed_(OpenMM_LangevinIntegrator*& target, int const& seed) {
    OpenMM_LangevinIntegrator_setRandomNumberSeed(target, seed);
}
OPENMM_EXPORT void OPENMM_LANGEVININTEGRATOR_SETRANDOMNUMBERSEED(OpenMM_LangevinIntegrator*& target, int const& seed) {
    OpenMM_LangevinIntegrator_setRandomNumberSeed(target, seed);
}
OPENMM_EXPORT void openmm_langevinintegrator_step_(OpenMM_LangevinIntegrator*& target, int const& steps) {
    OpenMM_LangevinIntegrator_step(target, steps);
}
OPENMM_EXPORT void OPENMM_LANGEVININTEGRATOR_STEP(OpenMM_LangevinIntegrator*& target, int const& steps) {
    OpenMM_LangevinIntegrator_step(target, steps);
}

/* OpenMM::BAOABLangevinIntegrator */
OPENMM_EXPORT void openmm_baoablangevinintegrator_create_(OpenMM_BAOABLangevinIntegrator*& result, double const& temperature, double const& frictionCoeff, double const& stepSize) {
    result = OpenMM_BAOABLangevinIntegrator_create(temperature, frictionCoeff, stepSize);
}
OPENMM_EXPORT void OPENMM_BAOABLANGEVININTEGRATOR_CREATE(OpenMM_BAOABLangevinIntegrator*& result, double const& temperature, double const& frictionCoeff, double const& stepSize) {
    result = OpenMM_BAOABLangevinIntegrator_create(temperature, frictionCoeff, stepSize);
}
OPENMM_EXPORT void openmm_baoablangevinintegrator_destroy_(OpenMM_BAOABLangevinIntegrator*& destroy) {
    OpenMM_BAOABLangevinIntegrator_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_BAOABLANGEVININTEGRATOR_DESTROY(OpenMM_BAOABLangevinIntegrator*& destroy) {
    OpenMM_BAOABLangevinIntegrator_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT double openmm_baoablangevinintegrator_gettemperature_(const OpenMM_BAOABLangevinIntegrator*& target) {
    return OpenMM_BAOABLangevinIntegrator_getTemperature(target);
}
OPENMM_EXPORT double OPENMM_BAOABLANGEVININTEGRATOR_GETTEMPERATURE(const OpenMM_BAOABLangevinIntegrator*& target) {
    return OpenMM_BAOABLangevinIntegrator_getTemperature(target);
}
OPENMM_EXPORT void openmm_baoablangevinintegrator_settemperature_(OpenMM_BAOABLangevinIntegrator*& target, double const& temp) {
    OpenMM_BAOABLangevinIntegrator_setTemperature(target, temp);
}
OPENMM_EXPORT void OPENMM_BAOABLANGEVININTEGRATOR_SETTEMPERATURE(OpenMM_BAOABLangevinIntegrator*& target, double const& temp) {
    OpenMM_BAOABLangevinIntegrator_setTemperature(target, temp);
}
OPENMM_EXPORT double openmm_baoablangevinintegrator_getfriction_(const OpenMM_BAOABLangevinIntegrator*& target) {
    return OpenMM_BAOABLangevinIntegrator_getFriction(target);
}
OPENMM_EXPORT double OPENMM_BAOABLANGEVININTEGRATOR_GETFRICTION(const OpenMM_BAOABLangevinIntegrator*& target) {
    return OpenMM_BAOABLangevinIntegrator_getFriction(target);
}
OPENMM_EXPORT void openmm_baoablangevinintegrator_setfriction_(OpenMM_BAOABLangevinIntegrator*& target, double const& coeff) {
    OpenMM_BAOABLangevinIntegrator_setFriction(target, coeff);
}
OPENMM_EXPORT void OPENMM_BAOABLANGEVININTEGRATOR_SETFRICTION(OpenMM_BAOABLangevinIntegrator*& target, double const& coeff) {
    OpenMM_BAOABLangevinIntegrator_setFriction(target, coeff);
}
OPENMM_EXPORT int openmm_baoablangevinintegrator_getrandomnumberseed_(const OpenMM_BAOABLangevinIntegrator*& target) {
    return OpenMM_BAOABLangevinIntegrator_getRandomNumberSeed(target);
}
OPENMM_EXPORT int OPENMM_BAOABLANGEVININTEGRATOR_GETRANDOMNUMBERSEED(const OpenMM_BAOABLangevinIntegrator*& target) {
    return OpenMM_BAOABLangevinIntegrator_getRandomNumberSeed(target);
}
OPENMM_EXPORT void openmm_baoablangevinintegrator_setrandomnumberseed_(OpenMM_BAOABLangevinIntegrator*& target, int const& seed) {
    OpenMM_BAOABLangevinIntegrator_setRandomNumberSeed(target, seed);
}
OPENMM_EXPORT void OPENMM_BAOABLANGEVININTEGRATOR_SETRANDOMNUMBERSEED(OpenMM_BAOABLangevinIntegrator*& target, int const& seed) {
    OpenMM_BAOABLangevinIntegrator_setRandomNumberSeed(target, seed);
}
OPENMM_EXPORT void openmm_baoablangevinintegrator_step_(OpenMM_BAOABLangevinIntegrator*& target, int const& steps) {
    OpenMM_BAOABLangevinIntegrator_step(target, steps);
}
OPENMM_EXPORT void OPENMM_BAOABLANGEVININTEGRATOR_STEP(OpenMM_BAOABLangevinIntegrator*& target, int const& steps) {
    OpenMM_BAOABLangevinIntegrator_step(target, steps);
}

/* OpenMM::TabulatedFunction */
OPENMM_EXPORT void openmm_tabulatedfunction_destroy_(OpenMM_TabulatedFunction*& destroy) {
    OpenMM_TabulatedFunction_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_TABULATEDFUNCTION_DESTROY(OpenMM_TabulatedFunction*& destroy) {
    OpenMM_TabulatedFunction_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void openmm_tabulatedfunction_copy_(const OpenMM_TabulatedFunction*& target, OpenMM_TabulatedFunction*& result) {
    result = OpenMM_TabulatedFunction_Copy(target);
}
OPENMM_EXPORT void OPENMM_TABULATEDFUNCTION_COPY(const OpenMM_TabulatedFunction*& target, OpenMM_TabulatedFunction*& result) {
    result = OpenMM_TabulatedFunction_Copy(target);
}

/* OpenMM::Context */
OPENMM_EXPORT void openmm_context_create_(OpenMM_Context*& result, const OpenMM_System*& system, OpenMM_Integrator*& integrator) {
    result = OpenMM_Context_create(system, integrator);
}
OPENMM_EXPORT void OPENMM_CONTEXT_CREATE(OpenMM_Context*& result, const OpenMM_System*& system, OpenMM_Integrator*& integrator) {
    result = OpenMM_Context_create(system, integrator);
}
OPENMM_EXPORT void openmm_context_create_2_(OpenMM_Context*& result, const OpenMM_System*& system, OpenMM_Integrator*& integrator, OpenMM_Platform*& platform) {
    result = OpenMM_Context_create_2(system, integrator, platform);
}
OPENMM_EXPORT void OPENMM_CONTEXT_CREATE_2(OpenMM_Context*& result, const OpenMM_System*& system, OpenMM_Integrator*& integrator, OpenMM_Platform*& platform) {
    result = OpenMM_Context_create_2(system, integrator, platform);
}
OPENMM_EXPORT void openmm_context_create_3_(OpenMM_Context*& result, const OpenMM_System*& system, OpenMM_Integrator*& integrator, OpenMM_Platform*& platform, const OpenMM_PropertyArray*& properties) {
    result = OpenMM_Context_create_3(system, integrator, platform, properties);
}
OPENMM_EXPORT void OPENMM_CONTEXT_CREATE_3(OpenMM_Context*& result, const OpenMM_System*& system, OpenMM_Integrator*& integrator, OpenMM_Platform*& platform, const OpenMM_PropertyArray*& properties) {
    result = OpenMM_Context_create_3(system, integrator, platform, properties);
}
OPENMM_EXPORT void openmm_context_destroy_(OpenMM_Context*& destroy) {
    OpenMM_Context_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_CONTEXT_DESTROY(OpenMM_Context*& destroy) {
    OpenMM_Context_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void openmm_context_getsystem_(const OpenMM_Context*& target, const OpenMM_System*& result) {
    result = OpenMM_Context_getSystem(target);
}
OPENMM_EXPORT void OPENMM_CONTEXT_GETSYSTEM(const OpenMM_Context*& target, const OpenMM_System*& result) {
    result = OpenMM_Context_getSystem(target);
}
OPENMM_EXPORT void openmm_context_getintegrator_(OpenMM_Context*& target, OpenMM_Integrator*& result) {
    result = OpenMM_Context_getIntegrator(target);
}
OPENMM_EXPORT void OPENMM_CONTEXT_GETINTEGRATOR(OpenMM_Context*& target, OpenMM_Integrator*& result) {
    result = OpenMM_Context_getIntegrator(target);
}
OPENMM_EXPORT void openmm_context_getplatform_(OpenMM_Context*& target, OpenMM_Platform*& result) {
    result = OpenMM_Context_getPlatform(target);
}
OPENMM_EXPORT void OPENMM_CONTEXT_GETPLATFORM(OpenMM_Context*& target, OpenMM_Platform*& result) {
    result = OpenMM_Context_getPlatform(target);
}
OPENMM_EXPORT void openmm_context_setstate_(OpenMM_Context*& target, const OpenMM_State*& state) {
    OpenMM_Context_setState(target, state);
}
OPENMM_EXPORT void OPENMM_CONTEXT_SETSTATE(OpenMM_Context*& target, const OpenMM_State*& state) {
    OpenMM_Context_setState(target, state);
}
OPENMM_EXPORT void openmm_context_settime_(OpenMM_Context*& target, double const& time) {
    OpenMM_Context_setTime(target, time);
}
OPENMM_EXPORT void OPENMM_CONTEXT_SETTIME(OpenMM_Context*& target, double const& time) {
    OpenMM_Context_setTime(target, time);
}
OPENMM_EXPORT void openmm_context_setpositions_(OpenMM_Context*& target, const OpenMM_Vec3Array*& positions) {
    OpenMM_Context_setPositions(target, positions);
}
OPENMM_EXPORT void OPENMM_CONTEXT_SETPOSITIONS(OpenMM_Context*& target, const OpenMM_Vec3Array*& positions) {
    OpenMM_Context_setPositions(target, positions);
}
OPENMM_EXPORT void openmm_context_setvelocities_(OpenMM_Context*& target, const OpenMM_Vec3Array*& velocities) {
    OpenMM_Context_setVelocities(target, velocities);
}
OPENMM_EXPORT void OPENMM_CONTEXT_SETVELOCITIES(OpenMM_Context*& target, const OpenMM_Vec3Array*& velocities) {
    OpenMM_Context_setVelocities(target, velocities);
}
OPENMM_EXPORT void openmm_context_setvelocitiestotemperature_(OpenMM_Context*& target, double const& temperature, int const& randomSeed) {
    OpenMM_Context_setVelocitiesToTemperature(target, temperature, randomSeed);
}
OPENMM_EXPORT void OPENMM_CONTEXT_SETVELOCITIESTOTEMPERATURE(OpenMM_Context*& target, double const& temperature, int const& randomSeed) {
    OpenMM_Context_setVelocitiesToTemperature(target, temperature, randomSeed);
}
OPENMM_EXPORT void openmm_context_getparameters_(const OpenMM_Context*& target, const OpenMM_ParameterArray*& result) {
    result = OpenMM_Context_getParameters(target);
}
OPENMM_EXPORT void OPENMM_CONTEXT_GETPARAMETERS(const OpenMM_Context*& target, const OpenMM_ParameterArray*& result) {
    result = OpenMM_Context_getParameters(target);
}
OPENMM_EXPORT double openmm_context_getparameter_(const OpenMM_Context*& target, const char* name, int name_length) {
    return OpenMM_Context_getParameter(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT double OPENMM_CONTEXT_GETPARAMETER(const OpenMM_Context*& target, const char* name, int name_length) {
    return OpenMM_Context_getParameter(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void openmm_context_setparameter_(OpenMM_Context*& target, const char* name, double const& value, int name_length) {
    OpenMM_Context_setParameter(target, makeString(name, name_length).c_str(), value);
}
OPENMM_EXPORT void OPENMM_CONTEXT_SETPARAMETER(OpenMM_Context*& target, const char* name, double const& value, int name_length) {
    OpenMM_Context_setParameter(target, makeString(name, name_length).c_str(), value);
}
OPENMM_EXPORT void openmm_context_setperiodicboxvectors_(OpenMM_Context*& target, const OpenMM_Vec3* a, const OpenMM_Vec3* b, const OpenMM_Vec3* c) {
    OpenMM_Context_setPeriodicBoxVectors(target, a, b, c);
}
OPENMM_EXPORT void OPENMM_CONTEXT_SETPERIODICBOXVECTORS(OpenMM_Context*& target, const OpenMM_Vec3* a, const OpenMM_Vec3* b, const OpenMM_Vec3* c) {
    OpenMM_Context_setPeriodicBoxVectors(target, a, b, c);
}
OPENMM_EXPORT void openmm_context_applyconstraints_(OpenMM_Context*& target, double const& tol) {
    OpenMM_Context_applyConstraints(target, tol);
}
OPENMM_EXPORT void OPENMM_CONTEXT_APPLYCONSTRAINTS(OpenMM_Context*& target, double const& tol) {
    OpenMM_Context_applyConstraints(target, tol);
}
OPENMM_EXPORT void openmm_context_applyvelocityconstraints_(OpenMM_Context*& target, double const& tol) {
    OpenMM_Context_applyVelocityConstraints(target, tol);
}
OPENMM_EXPORT void OPENMM_CONTEXT_APPLYVELOCITYCONSTRAINTS(OpenMM_Context*& target, double const& tol) {
    OpenMM_Context_applyVelocityConstraints(target, tol);
}
OPENMM_EXPORT void openmm_context_computevirtualsites_(OpenMM_Context*& target) {
    OpenMM_Context_computeVirtualSites(target);
}
OPENMM_EXPORT void OPENMM_CONTEXT_COMPUTEVIRTUALSITES(OpenMM_Context*& target) {
    OpenMM_Context_computeVirtualSites(target);
}
OPENMM_EXPORT void openmm_context_reinitialize_(OpenMM_Context*& target, OpenMM_Boolean& preserveState) {
    OpenMM_Context_reinitialize(target, preserveState);
}
OPENMM_EXPORT void OPENMM_CONTEXT_REINITIALIZE(OpenMM_Context*& target, OpenMM_Boolean& preserveState) {
    OpenMM_Context_reinitialize(target, preserveState);
}

/* OpenMM::CMMotionRemover */
OPENMM_EXPORT void openmm_cmmotionremover_create_(OpenMM_CMMotionRemover*& result, int const& frequency) {
    result = OpenMM_CMMotionRemover_create(frequency);
}
OPENMM_EXPORT void OPENMM_CMMOTIONREMOVER_CREATE(OpenMM_CMMotionRemover*& result, int const& frequency) {
    result = OpenMM_CMMotionRemover_create(frequency);
}
OPENMM_EXPORT void openmm_cmmotionremover_destroy_(OpenMM_CMMotionRemover*& destroy) {
    OpenMM_CMMotionRemover_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_CMMOTIONREMOVER_DESTROY(OpenMM_CMMotionRemover*& destroy) {
    OpenMM_CMMotionRemover_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT int openmm_cmmotionremover_getfrequency_(const OpenMM_CMMotionRemover*& target) {
    return OpenMM_CMMotionRemover_getFrequency(target);
}
OPENMM_EXPORT int OPENMM_CMMOTIONREMOVER_GETFREQUENCY(const OpenMM_CMMotionRemover*& target) {
    return OpenMM_CMMotionRemover_getFrequency(target);
}
OPENMM_EXPORT void openmm_cmmotionremover_setfrequency_(OpenMM_CMMotionRemover*& target, int const& freq) {
    OpenMM_CMMotionRemover_setFrequency(target, freq);
}
OPENMM_EXPORT void OPENMM_CMMOTIONREMOVER_SETFREQUENCY(OpenMM_CMMotionRemover*& target, int const& freq) {
    OpenMM_CMMotionRemover_setFrequency(target, freq);
}
OPENMM_EXPORT void openmm_cmmotionremover_usesperiodicboundaryconditions_(const OpenMM_CMMotionRemover*& target, OpenMM_Boolean& result) {
    result = OpenMM_CMMotionRemover_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_CMMOTIONREMOVER_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_CMMotionRemover*& target, OpenMM_Boolean& result) {
    result = OpenMM_CMMotionRemover_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::CustomBondForce */
OPENMM_EXPORT void openmm_custombondforce_create_(OpenMM_CustomBondForce*& result, const char* energy, int energy_length) {
    result = OpenMM_CustomBondForce_create(makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMBONDFORCE_CREATE(OpenMM_CustomBondForce*& result, const char* energy, int energy_length) {
    result = OpenMM_CustomBondForce_create(makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void openmm_custombondforce_destroy_(OpenMM_CustomBondForce*& destroy) {
    OpenMM_CustomBondForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_CUSTOMBONDFORCE_DESTROY(OpenMM_CustomBondForce*& destroy) {
    OpenMM_CustomBondForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT int openmm_custombondforce_getnumbonds_(const OpenMM_CustomBondForce*& target) {
    return OpenMM_CustomBondForce_getNumBonds(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMBONDFORCE_GETNUMBONDS(const OpenMM_CustomBondForce*& target) {
    return OpenMM_CustomBondForce_getNumBonds(target);
}
OPENMM_EXPORT int openmm_custombondforce_getnumperbondparameters_(const OpenMM_CustomBondForce*& target) {
    return OpenMM_CustomBondForce_getNumPerBondParameters(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMBONDFORCE_GETNUMPERBONDPARAMETERS(const OpenMM_CustomBondForce*& target) {
    return OpenMM_CustomBondForce_getNumPerBondParameters(target);
}
OPENMM_EXPORT int openmm_custombondforce_getnumglobalparameters_(const OpenMM_CustomBondForce*& target) {
    return OpenMM_CustomBondForce_getNumGlobalParameters(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMBONDFORCE_GETNUMGLOBALPARAMETERS(const OpenMM_CustomBondForce*& target) {
    return OpenMM_CustomBondForce_getNumGlobalParameters(target);
}
OPENMM_EXPORT int openmm_custombondforce_getnumenergyparameterderivatives_(const OpenMM_CustomBondForce*& target) {
    return OpenMM_CustomBondForce_getNumEnergyParameterDerivatives(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMBONDFORCE_GETNUMENERGYPARAMETERDERIVATIVES(const OpenMM_CustomBondForce*& target) {
    return OpenMM_CustomBondForce_getNumEnergyParameterDerivatives(target);
}
OPENMM_EXPORT void openmm_custombondforce_getenergyfunction_(const OpenMM_CustomBondForce*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomBondForce_getEnergyFunction(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMBONDFORCE_GETENERGYFUNCTION(const OpenMM_CustomBondForce*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomBondForce_getEnergyFunction(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_custombondforce_setenergyfunction_(OpenMM_CustomBondForce*& target, const char* energy, int energy_length) {
    OpenMM_CustomBondForce_setEnergyFunction(target, makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMBONDFORCE_SETENERGYFUNCTION(OpenMM_CustomBondForce*& target, const char* energy, int energy_length) {
    OpenMM_CustomBondForce_setEnergyFunction(target, makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT int openmm_custombondforce_addperbondparameter_(OpenMM_CustomBondForce*& target, const char* name, int name_length) {
    return OpenMM_CustomBondForce_addPerBondParameter(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT int OPENMM_CUSTOMBONDFORCE_ADDPERBONDPARAMETER(OpenMM_CustomBondForce*& target, const char* name, int name_length) {
    return OpenMM_CustomBondForce_addPerBondParameter(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void openmm_custombondforce_getperbondparametername_(const OpenMM_CustomBondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomBondForce_getPerBondParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMBONDFORCE_GETPERBONDPARAMETERNAME(const OpenMM_CustomBondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomBondForce_getPerBondParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_custombondforce_setperbondparametername_(OpenMM_CustomBondForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomBondForce_setPerBondParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMBONDFORCE_SETPERBONDPARAMETERNAME(OpenMM_CustomBondForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomBondForce_setPerBondParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT int openmm_custombondforce_addglobalparameter_(OpenMM_CustomBondForce*& target, const char* name, double const& defaultValue, int name_length) {
    return OpenMM_CustomBondForce_addGlobalParameter(target, makeString(name, name_length).c_str(), defaultValue);
}
OPENMM_EXPORT int OPENMM_CUSTOMBONDFORCE_ADDGLOBALPARAMETER(OpenMM_CustomBondForce*& target, const char* name, double const& defaultValue, int name_length) {
    return OpenMM_CustomBondForce_addGlobalParameter(target, makeString(name, name_length).c_str(), defaultValue);
}
OPENMM_EXPORT void openmm_custombondforce_getglobalparametername_(const OpenMM_CustomBondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomBondForce_getGlobalParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMBONDFORCE_GETGLOBALPARAMETERNAME(const OpenMM_CustomBondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomBondForce_getGlobalParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_custombondforce_setglobalparametername_(OpenMM_CustomBondForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomBondForce_setGlobalParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMBONDFORCE_SETGLOBALPARAMETERNAME(OpenMM_CustomBondForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomBondForce_setGlobalParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT double openmm_custombondforce_getglobalparameterdefaultvalue_(const OpenMM_CustomBondForce*& target, int const& index) {
    return OpenMM_CustomBondForce_getGlobalParameterDefaultValue(target, index);
}
OPENMM_EXPORT double OPENMM_CUSTOMBONDFORCE_GETGLOBALPARAMETERDEFAULTVALUE(const OpenMM_CustomBondForce*& target, int const& index) {
    return OpenMM_CustomBondForce_getGlobalParameterDefaultValue(target, index);
}
OPENMM_EXPORT void openmm_custombondforce_setglobalparameterdefaultvalue_(OpenMM_CustomBondForce*& target, int const& index, double const& defaultValue) {
    OpenMM_CustomBondForce_setGlobalParameterDefaultValue(target, index, defaultValue);
}
OPENMM_EXPORT void OPENMM_CUSTOMBONDFORCE_SETGLOBALPARAMETERDEFAULTVALUE(OpenMM_CustomBondForce*& target, int const& index, double const& defaultValue) {
    OpenMM_CustomBondForce_setGlobalParameterDefaultValue(target, index, defaultValue);
}
OPENMM_EXPORT void openmm_custombondforce_addenergyparameterderivative_(OpenMM_CustomBondForce*& target, const char* name, int name_length) {
    OpenMM_CustomBondForce_addEnergyParameterDerivative(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMBONDFORCE_ADDENERGYPARAMETERDERIVATIVE(OpenMM_CustomBondForce*& target, const char* name, int name_length) {
    OpenMM_CustomBondForce_addEnergyParameterDerivative(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void openmm_custombondforce_getenergyparameterderivativename_(const OpenMM_CustomBondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomBondForce_getEnergyParameterDerivativeName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMBONDFORCE_GETENERGYPARAMETERDERIVATIVENAME(const OpenMM_CustomBondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomBondForce_getEnergyParameterDerivativeName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT int openmm_custombondforce_addbond_(OpenMM_CustomBondForce*& target, int const& particle1, int const& particle2, const OpenMM_DoubleArray*& parameters) {
    return OpenMM_CustomBondForce_addBond(target, particle1, particle2, parameters);
}
OPENMM_EXPORT int OPENMM_CUSTOMBONDFORCE_ADDBOND(OpenMM_CustomBondForce*& target, int const& particle1, int const& particle2, const OpenMM_DoubleArray*& parameters) {
    return OpenMM_CustomBondForce_addBond(target, particle1, particle2, parameters);
}
OPENMM_EXPORT void openmm_custombondforce_getbondparameters_(const OpenMM_CustomBondForce*& target, int const& index, int* particle1, int* particle2, OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomBondForce_getBondParameters(target, index, particle1, particle2, parameters);
}
OPENMM_EXPORT void OPENMM_CUSTOMBONDFORCE_GETBONDPARAMETERS(const OpenMM_CustomBondForce*& target, int const& index, int* particle1, int* particle2, OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomBondForce_getBondParameters(target, index, particle1, particle2, parameters);
}
OPENMM_EXPORT void openmm_custombondforce_setbondparameters_(OpenMM_CustomBondForce*& target, int const& index, int const& particle1, int const& particle2, const OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomBondForce_setBondParameters(target, index, particle1, particle2, parameters);
}
OPENMM_EXPORT void OPENMM_CUSTOMBONDFORCE_SETBONDPARAMETERS(OpenMM_CustomBondForce*& target, int const& index, int const& particle1, int const& particle2, const OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomBondForce_setBondParameters(target, index, particle1, particle2, parameters);
}
OPENMM_EXPORT void openmm_custombondforce_updateparametersincontext_(OpenMM_CustomBondForce*& target, OpenMM_Context*& context) {
    OpenMM_CustomBondForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void OPENMM_CUSTOMBONDFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_CustomBondForce*& target, OpenMM_Context*& context) {
    OpenMM_CustomBondForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void openmm_custombondforce_setusesperiodicboundaryconditions_(OpenMM_CustomBondForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_CustomBondForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT void OPENMM_CUSTOMBONDFORCE_SETUSESPERIODICBOUNDARYCONDITIONS(OpenMM_CustomBondForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_CustomBondForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT void openmm_custombondforce_usesperiodicboundaryconditions_(const OpenMM_CustomBondForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CustomBondForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_CUSTOMBONDFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_CustomBondForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CustomBondForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::MonteCarloMembraneBarostat */
OPENMM_EXPORT void openmm_montecarlomembranebarostat_create_(OpenMM_MonteCarloMembraneBarostat*& result, double const& defaultPressure, double const& defaultSurfaceTension, double const& defaultTemperature, OpenMM_MonteCarloMembraneBarostat_XYMode& xymode, OpenMM_MonteCarloMembraneBarostat_ZMode& zmode, int const& frequency) {
    result = OpenMM_MonteCarloMembraneBarostat_create(defaultPressure, defaultSurfaceTension, defaultTemperature, xymode, zmode, frequency);
}
OPENMM_EXPORT void OPENMM_MONTECARLOMEMBRANEBAROSTAT_CREATE(OpenMM_MonteCarloMembraneBarostat*& result, double const& defaultPressure, double const& defaultSurfaceTension, double const& defaultTemperature, OpenMM_MonteCarloMembraneBarostat_XYMode& xymode, OpenMM_MonteCarloMembraneBarostat_ZMode& zmode, int const& frequency) {
    result = OpenMM_MonteCarloMembraneBarostat_create(defaultPressure, defaultSurfaceTension, defaultTemperature, xymode, zmode, frequency);
}
OPENMM_EXPORT void openmm_montecarlomembranebarostat_destroy_(OpenMM_MonteCarloMembraneBarostat*& destroy) {
    OpenMM_MonteCarloMembraneBarostat_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_MONTECARLOMEMBRANEBAROSTAT_DESTROY(OpenMM_MonteCarloMembraneBarostat*& destroy) {
    OpenMM_MonteCarloMembraneBarostat_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void openmm_montecarlomembranebarostat_pressure_(char* result, int result_length) {
    const char* result_chars = OpenMM_MonteCarloMembraneBarostat_Pressure();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_MONTECARLOMEMBRANEBAROSTAT_PRESSURE(char* result, int result_length) {
    const char* result_chars = OpenMM_MonteCarloMembraneBarostat_Pressure();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_montecarlomembranebarostat_surfacetension_(char* result, int result_length) {
    const char* result_chars = OpenMM_MonteCarloMembraneBarostat_SurfaceTension();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_MONTECARLOMEMBRANEBAROSTAT_SURFACETENSION(char* result, int result_length) {
    const char* result_chars = OpenMM_MonteCarloMembraneBarostat_SurfaceTension();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_montecarlomembranebarostat_temperature_(char* result, int result_length) {
    const char* result_chars = OpenMM_MonteCarloMembraneBarostat_Temperature();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_MONTECARLOMEMBRANEBAROSTAT_TEMPERATURE(char* result, int result_length) {
    const char* result_chars = OpenMM_MonteCarloMembraneBarostat_Temperature();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT double openmm_montecarlomembranebarostat_getdefaultpressure_(const OpenMM_MonteCarloMembraneBarostat*& target) {
    return OpenMM_MonteCarloMembraneBarostat_getDefaultPressure(target);
}
OPENMM_EXPORT double OPENMM_MONTECARLOMEMBRANEBAROSTAT_GETDEFAULTPRESSURE(const OpenMM_MonteCarloMembraneBarostat*& target) {
    return OpenMM_MonteCarloMembraneBarostat_getDefaultPressure(target);
}
OPENMM_EXPORT void openmm_montecarlomembranebarostat_setdefaultpressure_(OpenMM_MonteCarloMembraneBarostat*& target, double const& pressure) {
    OpenMM_MonteCarloMembraneBarostat_setDefaultPressure(target, pressure);
}
OPENMM_EXPORT void OPENMM_MONTECARLOMEMBRANEBAROSTAT_SETDEFAULTPRESSURE(OpenMM_MonteCarloMembraneBarostat*& target, double const& pressure) {
    OpenMM_MonteCarloMembraneBarostat_setDefaultPressure(target, pressure);
}
OPENMM_EXPORT double openmm_montecarlomembranebarostat_getdefaultsurfacetension_(const OpenMM_MonteCarloMembraneBarostat*& target) {
    return OpenMM_MonteCarloMembraneBarostat_getDefaultSurfaceTension(target);
}
OPENMM_EXPORT double OPENMM_MONTECARLOMEMBRANEBAROSTAT_GETDEFAULTSURFACETENSION(const OpenMM_MonteCarloMembraneBarostat*& target) {
    return OpenMM_MonteCarloMembraneBarostat_getDefaultSurfaceTension(target);
}
OPENMM_EXPORT void openmm_montecarlomembranebarostat_setdefaultsurfacetension_(OpenMM_MonteCarloMembraneBarostat*& target, double const& surfaceTension) {
    OpenMM_MonteCarloMembraneBarostat_setDefaultSurfaceTension(target, surfaceTension);
}
OPENMM_EXPORT void OPENMM_MONTECARLOMEMBRANEBAROSTAT_SETDEFAULTSURFACETENSION(OpenMM_MonteCarloMembraneBarostat*& target, double const& surfaceTension) {
    OpenMM_MonteCarloMembraneBarostat_setDefaultSurfaceTension(target, surfaceTension);
}
OPENMM_EXPORT int openmm_montecarlomembranebarostat_getfrequency_(const OpenMM_MonteCarloMembraneBarostat*& target) {
    return OpenMM_MonteCarloMembraneBarostat_getFrequency(target);
}
OPENMM_EXPORT int OPENMM_MONTECARLOMEMBRANEBAROSTAT_GETFREQUENCY(const OpenMM_MonteCarloMembraneBarostat*& target) {
    return OpenMM_MonteCarloMembraneBarostat_getFrequency(target);
}
OPENMM_EXPORT void openmm_montecarlomembranebarostat_setfrequency_(OpenMM_MonteCarloMembraneBarostat*& target, int const& freq) {
    OpenMM_MonteCarloMembraneBarostat_setFrequency(target, freq);
}
OPENMM_EXPORT void OPENMM_MONTECARLOMEMBRANEBAROSTAT_SETFREQUENCY(OpenMM_MonteCarloMembraneBarostat*& target, int const& freq) {
    OpenMM_MonteCarloMembraneBarostat_setFrequency(target, freq);
}
OPENMM_EXPORT double openmm_montecarlomembranebarostat_getdefaulttemperature_(const OpenMM_MonteCarloMembraneBarostat*& target) {
    return OpenMM_MonteCarloMembraneBarostat_getDefaultTemperature(target);
}
OPENMM_EXPORT double OPENMM_MONTECARLOMEMBRANEBAROSTAT_GETDEFAULTTEMPERATURE(const OpenMM_MonteCarloMembraneBarostat*& target) {
    return OpenMM_MonteCarloMembraneBarostat_getDefaultTemperature(target);
}
OPENMM_EXPORT void openmm_montecarlomembranebarostat_setdefaulttemperature_(OpenMM_MonteCarloMembraneBarostat*& target, double const& temp) {
    OpenMM_MonteCarloMembraneBarostat_setDefaultTemperature(target, temp);
}
OPENMM_EXPORT void OPENMM_MONTECARLOMEMBRANEBAROSTAT_SETDEFAULTTEMPERATURE(OpenMM_MonteCarloMembraneBarostat*& target, double const& temp) {
    OpenMM_MonteCarloMembraneBarostat_setDefaultTemperature(target, temp);
}
OPENMM_EXPORT void openmm_montecarlomembranebarostat_getxymode_(const OpenMM_MonteCarloMembraneBarostat*& target, OpenMM_MonteCarloMembraneBarostat_XYMode& result) {
    result = OpenMM_MonteCarloMembraneBarostat_getXYMode(target);
}
OPENMM_EXPORT void OPENMM_MONTECARLOMEMBRANEBAROSTAT_GETXYMODE(const OpenMM_MonteCarloMembraneBarostat*& target, OpenMM_MonteCarloMembraneBarostat_XYMode& result) {
    result = OpenMM_MonteCarloMembraneBarostat_getXYMode(target);
}
OPENMM_EXPORT void openmm_montecarlomembranebarostat_setxymode_(OpenMM_MonteCarloMembraneBarostat*& target, OpenMM_MonteCarloMembraneBarostat_XYMode& mode) {
    OpenMM_MonteCarloMembraneBarostat_setXYMode(target, mode);
}
OPENMM_EXPORT void OPENMM_MONTECARLOMEMBRANEBAROSTAT_SETXYMODE(OpenMM_MonteCarloMembraneBarostat*& target, OpenMM_MonteCarloMembraneBarostat_XYMode& mode) {
    OpenMM_MonteCarloMembraneBarostat_setXYMode(target, mode);
}
OPENMM_EXPORT void openmm_montecarlomembranebarostat_getzmode_(const OpenMM_MonteCarloMembraneBarostat*& target, OpenMM_MonteCarloMembraneBarostat_ZMode& result) {
    result = OpenMM_MonteCarloMembraneBarostat_getZMode(target);
}
OPENMM_EXPORT void OPENMM_MONTECARLOMEMBRANEBAROSTAT_GETZMODE(const OpenMM_MonteCarloMembraneBarostat*& target, OpenMM_MonteCarloMembraneBarostat_ZMode& result) {
    result = OpenMM_MonteCarloMembraneBarostat_getZMode(target);
}
OPENMM_EXPORT void openmm_montecarlomembranebarostat_setzmode_(OpenMM_MonteCarloMembraneBarostat*& target, OpenMM_MonteCarloMembraneBarostat_ZMode& mode) {
    OpenMM_MonteCarloMembraneBarostat_setZMode(target, mode);
}
OPENMM_EXPORT void OPENMM_MONTECARLOMEMBRANEBAROSTAT_SETZMODE(OpenMM_MonteCarloMembraneBarostat*& target, OpenMM_MonteCarloMembraneBarostat_ZMode& mode) {
    OpenMM_MonteCarloMembraneBarostat_setZMode(target, mode);
}
OPENMM_EXPORT int openmm_montecarlomembranebarostat_getrandomnumberseed_(const OpenMM_MonteCarloMembraneBarostat*& target) {
    return OpenMM_MonteCarloMembraneBarostat_getRandomNumberSeed(target);
}
OPENMM_EXPORT int OPENMM_MONTECARLOMEMBRANEBAROSTAT_GETRANDOMNUMBERSEED(const OpenMM_MonteCarloMembraneBarostat*& target) {
    return OpenMM_MonteCarloMembraneBarostat_getRandomNumberSeed(target);
}
OPENMM_EXPORT void openmm_montecarlomembranebarostat_setrandomnumberseed_(OpenMM_MonteCarloMembraneBarostat*& target, int const& seed) {
    OpenMM_MonteCarloMembraneBarostat_setRandomNumberSeed(target, seed);
}
OPENMM_EXPORT void OPENMM_MONTECARLOMEMBRANEBAROSTAT_SETRANDOMNUMBERSEED(OpenMM_MonteCarloMembraneBarostat*& target, int const& seed) {
    OpenMM_MonteCarloMembraneBarostat_setRandomNumberSeed(target, seed);
}
OPENMM_EXPORT void openmm_montecarlomembranebarostat_usesperiodicboundarycondition_(const OpenMM_MonteCarloMembraneBarostat*& target, OpenMM_Boolean& result) {
    result = OpenMM_MonteCarloMembraneBarostat_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_MONTECARLOMEMBRANEBAROSTAT_USESPERIODICBOUNDARYCONDITION(const OpenMM_MonteCarloMembraneBarostat*& target, OpenMM_Boolean& result) {
    result = OpenMM_MonteCarloMembraneBarostat_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::VirtualSite */
OPENMM_EXPORT void openmm_virtualsite_destroy_(OpenMM_VirtualSite*& destroy) {
    OpenMM_VirtualSite_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_VIRTUALSITE_DESTROY(OpenMM_VirtualSite*& destroy) {
    OpenMM_VirtualSite_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT int openmm_virtualsite_getnumparticles_(const OpenMM_VirtualSite*& target) {
    return OpenMM_VirtualSite_getNumParticles(target);
}
OPENMM_EXPORT int OPENMM_VIRTUALSITE_GETNUMPARTICLES(const OpenMM_VirtualSite*& target) {
    return OpenMM_VirtualSite_getNumParticles(target);
}
OPENMM_EXPORT int openmm_virtualsite_getparticle_(const OpenMM_VirtualSite*& target, int const& particle) {
    return OpenMM_VirtualSite_getParticle(target, particle);
}
OPENMM_EXPORT int OPENMM_VIRTUALSITE_GETPARTICLE(const OpenMM_VirtualSite*& target, int const& particle) {
    return OpenMM_VirtualSite_getParticle(target, particle);
}

/* OpenMM::OutOfPlaneSite */
OPENMM_EXPORT void openmm_outofplanesite_create_(OpenMM_OutOfPlaneSite*& result, int const& particle1, int const& particle2, int const& particle3, double const& weight12, double const& weight13, double const& weightCross) {
    result = OpenMM_OutOfPlaneSite_create(particle1, particle2, particle3, weight12, weight13, weightCross);
}
OPENMM_EXPORT void OPENMM_OUTOFPLANESITE_CREATE(OpenMM_OutOfPlaneSite*& result, int const& particle1, int const& particle2, int const& particle3, double const& weight12, double const& weight13, double const& weightCross) {
    result = OpenMM_OutOfPlaneSite_create(particle1, particle2, particle3, weight12, weight13, weightCross);
}
OPENMM_EXPORT void openmm_outofplanesite_destroy_(OpenMM_OutOfPlaneSite*& destroy) {
    OpenMM_OutOfPlaneSite_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_OUTOFPLANESITE_DESTROY(OpenMM_OutOfPlaneSite*& destroy) {
    OpenMM_OutOfPlaneSite_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT double openmm_outofplanesite_getweight12_(const OpenMM_OutOfPlaneSite*& target) {
    return OpenMM_OutOfPlaneSite_getWeight12(target);
}
OPENMM_EXPORT double OPENMM_OUTOFPLANESITE_GETWEIGHT12(const OpenMM_OutOfPlaneSite*& target) {
    return OpenMM_OutOfPlaneSite_getWeight12(target);
}
OPENMM_EXPORT double openmm_outofplanesite_getweight13_(const OpenMM_OutOfPlaneSite*& target) {
    return OpenMM_OutOfPlaneSite_getWeight13(target);
}
OPENMM_EXPORT double OPENMM_OUTOFPLANESITE_GETWEIGHT13(const OpenMM_OutOfPlaneSite*& target) {
    return OpenMM_OutOfPlaneSite_getWeight13(target);
}
OPENMM_EXPORT double openmm_outofplanesite_getweightcross_(const OpenMM_OutOfPlaneSite*& target) {
    return OpenMM_OutOfPlaneSite_getWeightCross(target);
}
OPENMM_EXPORT double OPENMM_OUTOFPLANESITE_GETWEIGHTCROSS(const OpenMM_OutOfPlaneSite*& target) {
    return OpenMM_OutOfPlaneSite_getWeightCross(target);
}

/* OpenMM::CustomTorsionForce */
OPENMM_EXPORT void openmm_customtorsionforce_create_(OpenMM_CustomTorsionForce*& result, const char* energy, int energy_length) {
    result = OpenMM_CustomTorsionForce_create(makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMTORSIONFORCE_CREATE(OpenMM_CustomTorsionForce*& result, const char* energy, int energy_length) {
    result = OpenMM_CustomTorsionForce_create(makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void openmm_customtorsionforce_destroy_(OpenMM_CustomTorsionForce*& destroy) {
    OpenMM_CustomTorsionForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_CUSTOMTORSIONFORCE_DESTROY(OpenMM_CustomTorsionForce*& destroy) {
    OpenMM_CustomTorsionForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT int openmm_customtorsionforce_getnumtorsions_(const OpenMM_CustomTorsionForce*& target) {
    return OpenMM_CustomTorsionForce_getNumTorsions(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMTORSIONFORCE_GETNUMTORSIONS(const OpenMM_CustomTorsionForce*& target) {
    return OpenMM_CustomTorsionForce_getNumTorsions(target);
}
OPENMM_EXPORT int openmm_customtorsionforce_getnumpertorsionparameters_(const OpenMM_CustomTorsionForce*& target) {
    return OpenMM_CustomTorsionForce_getNumPerTorsionParameters(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMTORSIONFORCE_GETNUMPERTORSIONPARAMETERS(const OpenMM_CustomTorsionForce*& target) {
    return OpenMM_CustomTorsionForce_getNumPerTorsionParameters(target);
}
OPENMM_EXPORT int openmm_customtorsionforce_getnumglobalparameters_(const OpenMM_CustomTorsionForce*& target) {
    return OpenMM_CustomTorsionForce_getNumGlobalParameters(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMTORSIONFORCE_GETNUMGLOBALPARAMETERS(const OpenMM_CustomTorsionForce*& target) {
    return OpenMM_CustomTorsionForce_getNumGlobalParameters(target);
}
OPENMM_EXPORT int openmm_customtorsionforce_getnumenergyparameterderivatives_(const OpenMM_CustomTorsionForce*& target) {
    return OpenMM_CustomTorsionForce_getNumEnergyParameterDerivatives(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMTORSIONFORCE_GETNUMENERGYPARAMETERDERIVATIVES(const OpenMM_CustomTorsionForce*& target) {
    return OpenMM_CustomTorsionForce_getNumEnergyParameterDerivatives(target);
}
OPENMM_EXPORT void openmm_customtorsionforce_getenergyfunction_(const OpenMM_CustomTorsionForce*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomTorsionForce_getEnergyFunction(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMTORSIONFORCE_GETENERGYFUNCTION(const OpenMM_CustomTorsionForce*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomTorsionForce_getEnergyFunction(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customtorsionforce_setenergyfunction_(OpenMM_CustomTorsionForce*& target, const char* energy, int energy_length) {
    OpenMM_CustomTorsionForce_setEnergyFunction(target, makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMTORSIONFORCE_SETENERGYFUNCTION(OpenMM_CustomTorsionForce*& target, const char* energy, int energy_length) {
    OpenMM_CustomTorsionForce_setEnergyFunction(target, makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT int openmm_customtorsionforce_addpertorsionparameter_(OpenMM_CustomTorsionForce*& target, const char* name, int name_length) {
    return OpenMM_CustomTorsionForce_addPerTorsionParameter(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT int OPENMM_CUSTOMTORSIONFORCE_ADDPERTORSIONPARAMETER(OpenMM_CustomTorsionForce*& target, const char* name, int name_length) {
    return OpenMM_CustomTorsionForce_addPerTorsionParameter(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void openmm_customtorsionforce_getpertorsionparametername_(const OpenMM_CustomTorsionForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomTorsionForce_getPerTorsionParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMTORSIONFORCE_GETPERTORSIONPARAMETERNAME(const OpenMM_CustomTorsionForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomTorsionForce_getPerTorsionParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customtorsionforce_setpertorsionparametername_(OpenMM_CustomTorsionForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomTorsionForce_setPerTorsionParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMTORSIONFORCE_SETPERTORSIONPARAMETERNAME(OpenMM_CustomTorsionForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomTorsionForce_setPerTorsionParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT int openmm_customtorsionforce_addglobalparameter_(OpenMM_CustomTorsionForce*& target, const char* name, double const& defaultValue, int name_length) {
    return OpenMM_CustomTorsionForce_addGlobalParameter(target, makeString(name, name_length).c_str(), defaultValue);
}
OPENMM_EXPORT int OPENMM_CUSTOMTORSIONFORCE_ADDGLOBALPARAMETER(OpenMM_CustomTorsionForce*& target, const char* name, double const& defaultValue, int name_length) {
    return OpenMM_CustomTorsionForce_addGlobalParameter(target, makeString(name, name_length).c_str(), defaultValue);
}
OPENMM_EXPORT void openmm_customtorsionforce_getglobalparametername_(const OpenMM_CustomTorsionForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomTorsionForce_getGlobalParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMTORSIONFORCE_GETGLOBALPARAMETERNAME(const OpenMM_CustomTorsionForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomTorsionForce_getGlobalParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customtorsionforce_setglobalparametername_(OpenMM_CustomTorsionForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomTorsionForce_setGlobalParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMTORSIONFORCE_SETGLOBALPARAMETERNAME(OpenMM_CustomTorsionForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomTorsionForce_setGlobalParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT double openmm_customtorsionforce_getglobalparameterdefaultvalue_(const OpenMM_CustomTorsionForce*& target, int const& index) {
    return OpenMM_CustomTorsionForce_getGlobalParameterDefaultValue(target, index);
}
OPENMM_EXPORT double OPENMM_CUSTOMTORSIONFORCE_GETGLOBALPARAMETERDEFAULTVALUE(const OpenMM_CustomTorsionForce*& target, int const& index) {
    return OpenMM_CustomTorsionForce_getGlobalParameterDefaultValue(target, index);
}
OPENMM_EXPORT void openmm_customtorsionforce_setglobalparameterdefaultvalue_(OpenMM_CustomTorsionForce*& target, int const& index, double const& defaultValue) {
    OpenMM_CustomTorsionForce_setGlobalParameterDefaultValue(target, index, defaultValue);
}
OPENMM_EXPORT void OPENMM_CUSTOMTORSIONFORCE_SETGLOBALPARAMETERDEFAULTVALUE(OpenMM_CustomTorsionForce*& target, int const& index, double const& defaultValue) {
    OpenMM_CustomTorsionForce_setGlobalParameterDefaultValue(target, index, defaultValue);
}
OPENMM_EXPORT void openmm_customtorsionforce_addenergyparameterderivative_(OpenMM_CustomTorsionForce*& target, const char* name, int name_length) {
    OpenMM_CustomTorsionForce_addEnergyParameterDerivative(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMTORSIONFORCE_ADDENERGYPARAMETERDERIVATIVE(OpenMM_CustomTorsionForce*& target, const char* name, int name_length) {
    OpenMM_CustomTorsionForce_addEnergyParameterDerivative(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void openmm_customtorsionforce_getenergyparameterderivativename_(const OpenMM_CustomTorsionForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomTorsionForce_getEnergyParameterDerivativeName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMTORSIONFORCE_GETENERGYPARAMETERDERIVATIVENAME(const OpenMM_CustomTorsionForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomTorsionForce_getEnergyParameterDerivativeName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT int openmm_customtorsionforce_addtorsion_(OpenMM_CustomTorsionForce*& target, int const& particle1, int const& particle2, int const& particle3, int const& particle4, const OpenMM_DoubleArray*& parameters) {
    return OpenMM_CustomTorsionForce_addTorsion(target, particle1, particle2, particle3, particle4, parameters);
}
OPENMM_EXPORT int OPENMM_CUSTOMTORSIONFORCE_ADDTORSION(OpenMM_CustomTorsionForce*& target, int const& particle1, int const& particle2, int const& particle3, int const& particle4, const OpenMM_DoubleArray*& parameters) {
    return OpenMM_CustomTorsionForce_addTorsion(target, particle1, particle2, particle3, particle4, parameters);
}
OPENMM_EXPORT void openmm_customtorsionforce_gettorsionparameters_(const OpenMM_CustomTorsionForce*& target, int const& index, int* particle1, int* particle2, int* particle3, int* particle4, OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomTorsionForce_getTorsionParameters(target, index, particle1, particle2, particle3, particle4, parameters);
}
OPENMM_EXPORT void OPENMM_CUSTOMTORSIONFORCE_GETTORSIONPARAMETERS(const OpenMM_CustomTorsionForce*& target, int const& index, int* particle1, int* particle2, int* particle3, int* particle4, OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomTorsionForce_getTorsionParameters(target, index, particle1, particle2, particle3, particle4, parameters);
}
OPENMM_EXPORT void openmm_customtorsionforce_settorsionparameters_(OpenMM_CustomTorsionForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, int const& particle4, const OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomTorsionForce_setTorsionParameters(target, index, particle1, particle2, particle3, particle4, parameters);
}
OPENMM_EXPORT void OPENMM_CUSTOMTORSIONFORCE_SETTORSIONPARAMETERS(OpenMM_CustomTorsionForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, int const& particle4, const OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomTorsionForce_setTorsionParameters(target, index, particle1, particle2, particle3, particle4, parameters);
}
OPENMM_EXPORT void openmm_customtorsionforce_updateparametersincontext_(OpenMM_CustomTorsionForce*& target, OpenMM_Context*& context) {
    OpenMM_CustomTorsionForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void OPENMM_CUSTOMTORSIONFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_CustomTorsionForce*& target, OpenMM_Context*& context) {
    OpenMM_CustomTorsionForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void openmm_customtorsionforce_setusesperiodicboundaryconditions_(OpenMM_CustomTorsionForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_CustomTorsionForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT void OPENMM_CUSTOMTORSIONFORCE_SETUSESPERIODICBOUNDARYCONDITIONS(OpenMM_CustomTorsionForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_CustomTorsionForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT void openmm_customtorsionforce_usesperiodicboundaryconditions_(const OpenMM_CustomTorsionForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CustomTorsionForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_CUSTOMTORSIONFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_CustomTorsionForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CustomTorsionForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::Discrete2DFunction */
OPENMM_EXPORT void openmm_discrete2dfunction_create_(OpenMM_Discrete2DFunction*& result, int const& xsize, int const& ysize, const OpenMM_DoubleArray*& values) {
    result = OpenMM_Discrete2DFunction_create(xsize, ysize, values);
}
OPENMM_EXPORT void OPENMM_DISCRETE2DFUNCTION_CREATE(OpenMM_Discrete2DFunction*& result, int const& xsize, int const& ysize, const OpenMM_DoubleArray*& values) {
    result = OpenMM_Discrete2DFunction_create(xsize, ysize, values);
}
OPENMM_EXPORT void openmm_discrete2dfunction_destroy_(OpenMM_Discrete2DFunction*& destroy) {
    OpenMM_Discrete2DFunction_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_DISCRETE2DFUNCTION_DESTROY(OpenMM_Discrete2DFunction*& destroy) {
    OpenMM_Discrete2DFunction_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void openmm_discrete2dfunction_getfunctionparameters_(const OpenMM_Discrete2DFunction*& target, int* xsize, int* ysize, OpenMM_DoubleArray*& values) {
    OpenMM_Discrete2DFunction_getFunctionParameters(target, xsize, ysize, values);
}
OPENMM_EXPORT void OPENMM_DISCRETE2DFUNCTION_GETFUNCTIONPARAMETERS(const OpenMM_Discrete2DFunction*& target, int* xsize, int* ysize, OpenMM_DoubleArray*& values) {
    OpenMM_Discrete2DFunction_getFunctionParameters(target, xsize, ysize, values);
}
OPENMM_EXPORT void openmm_discrete2dfunction_setfunctionparameters_(OpenMM_Discrete2DFunction*& target, int const& xsize, int const& ysize, const OpenMM_DoubleArray*& values) {
    OpenMM_Discrete2DFunction_setFunctionParameters(target, xsize, ysize, values);
}
OPENMM_EXPORT void OPENMM_DISCRETE2DFUNCTION_SETFUNCTIONPARAMETERS(OpenMM_Discrete2DFunction*& target, int const& xsize, int const& ysize, const OpenMM_DoubleArray*& values) {
    OpenMM_Discrete2DFunction_setFunctionParameters(target, xsize, ysize, values);
}
OPENMM_EXPORT void openmm_discrete2dfunction_copy_(const OpenMM_Discrete2DFunction*& target, OpenMM_Discrete2DFunction*& result) {
    result = OpenMM_Discrete2DFunction_Copy(target);
}
OPENMM_EXPORT void OPENMM_DISCRETE2DFUNCTION_COPY(const OpenMM_Discrete2DFunction*& target, OpenMM_Discrete2DFunction*& result) {
    result = OpenMM_Discrete2DFunction_Copy(target);
}

/* OpenMM::AndersenThermostat */
OPENMM_EXPORT void openmm_andersenthermostat_create_(OpenMM_AndersenThermostat*& result, double const& defaultTemperature, double const& defaultCollisionFrequency) {
    result = OpenMM_AndersenThermostat_create(defaultTemperature, defaultCollisionFrequency);
}
OPENMM_EXPORT void OPENMM_ANDERSENTHERMOSTAT_CREATE(OpenMM_AndersenThermostat*& result, double const& defaultTemperature, double const& defaultCollisionFrequency) {
    result = OpenMM_AndersenThermostat_create(defaultTemperature, defaultCollisionFrequency);
}
OPENMM_EXPORT void openmm_andersenthermostat_destroy_(OpenMM_AndersenThermostat*& destroy) {
    OpenMM_AndersenThermostat_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_ANDERSENTHERMOSTAT_DESTROY(OpenMM_AndersenThermostat*& destroy) {
    OpenMM_AndersenThermostat_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void openmm_andersenthermostat_temperature_(char* result, int result_length) {
    const char* result_chars = OpenMM_AndersenThermostat_Temperature();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_ANDERSENTHERMOSTAT_TEMPERATURE(char* result, int result_length) {
    const char* result_chars = OpenMM_AndersenThermostat_Temperature();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_andersenthermostat_collisionfrequency_(char* result, int result_length) {
    const char* result_chars = OpenMM_AndersenThermostat_CollisionFrequency();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_ANDERSENTHERMOSTAT_COLLISIONFREQUENCY(char* result, int result_length) {
    const char* result_chars = OpenMM_AndersenThermostat_CollisionFrequency();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT double openmm_andersenthermostat_getdefaulttemperature_(const OpenMM_AndersenThermostat*& target) {
    return OpenMM_AndersenThermostat_getDefaultTemperature(target);
}
OPENMM_EXPORT double OPENMM_ANDERSENTHERMOSTAT_GETDEFAULTTEMPERATURE(const OpenMM_AndersenThermostat*& target) {
    return OpenMM_AndersenThermostat_getDefaultTemperature(target);
}
OPENMM_EXPORT void openmm_andersenthermostat_setdefaulttemperature_(OpenMM_AndersenThermostat*& target, double const& temperature) {
    OpenMM_AndersenThermostat_setDefaultTemperature(target, temperature);
}
OPENMM_EXPORT void OPENMM_ANDERSENTHERMOSTAT_SETDEFAULTTEMPERATURE(OpenMM_AndersenThermostat*& target, double const& temperature) {
    OpenMM_AndersenThermostat_setDefaultTemperature(target, temperature);
}
OPENMM_EXPORT double openmm_andersenthermostat_getdefaultcollisionfrequency_(const OpenMM_AndersenThermostat*& target) {
    return OpenMM_AndersenThermostat_getDefaultCollisionFrequency(target);
}
OPENMM_EXPORT double OPENMM_ANDERSENTHERMOSTAT_GETDEFAULTCOLLISIONFREQUENCY(const OpenMM_AndersenThermostat*& target) {
    return OpenMM_AndersenThermostat_getDefaultCollisionFrequency(target);
}
OPENMM_EXPORT void openmm_andersenthermostat_setdefaultcollisionfrequency_(OpenMM_AndersenThermostat*& target, double const& frequency) {
    OpenMM_AndersenThermostat_setDefaultCollisionFrequency(target, frequency);
}
OPENMM_EXPORT void OPENMM_ANDERSENTHERMOSTAT_SETDEFAULTCOLLISIONFREQUENCY(OpenMM_AndersenThermostat*& target, double const& frequency) {
    OpenMM_AndersenThermostat_setDefaultCollisionFrequency(target, frequency);
}
OPENMM_EXPORT int openmm_andersenthermostat_getrandomnumberseed_(const OpenMM_AndersenThermostat*& target) {
    return OpenMM_AndersenThermostat_getRandomNumberSeed(target);
}
OPENMM_EXPORT int OPENMM_ANDERSENTHERMOSTAT_GETRANDOMNUMBERSEED(const OpenMM_AndersenThermostat*& target) {
    return OpenMM_AndersenThermostat_getRandomNumberSeed(target);
}
OPENMM_EXPORT void openmm_andersenthermostat_setrandomnumberseed_(OpenMM_AndersenThermostat*& target, int const& seed) {
    OpenMM_AndersenThermostat_setRandomNumberSeed(target, seed);
}
OPENMM_EXPORT void OPENMM_ANDERSENTHERMOSTAT_SETRANDOMNUMBERSEED(OpenMM_AndersenThermostat*& target, int const& seed) {
    OpenMM_AndersenThermostat_setRandomNumberSeed(target, seed);
}
OPENMM_EXPORT void openmm_andersenthermostat_usesperiodicboundaryconditions_(const OpenMM_AndersenThermostat*& target, OpenMM_Boolean& result) {
    result = OpenMM_AndersenThermostat_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_ANDERSENTHERMOSTAT_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_AndersenThermostat*& target, OpenMM_Boolean& result) {
    result = OpenMM_AndersenThermostat_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::ThreeParticleAverageSite */
OPENMM_EXPORT void openmm_threeparticleaveragesite_create_(OpenMM_ThreeParticleAverageSite*& result, int const& particle1, int const& particle2, int const& particle3, double const& weight1, double const& weight2, double const& weight3) {
    result = OpenMM_ThreeParticleAverageSite_create(particle1, particle2, particle3, weight1, weight2, weight3);
}
OPENMM_EXPORT void OPENMM_THREEPARTICLEAVERAGESITE_CREATE(OpenMM_ThreeParticleAverageSite*& result, int const& particle1, int const& particle2, int const& particle3, double const& weight1, double const& weight2, double const& weight3) {
    result = OpenMM_ThreeParticleAverageSite_create(particle1, particle2, particle3, weight1, weight2, weight3);
}
OPENMM_EXPORT void openmm_threeparticleaveragesite_destroy_(OpenMM_ThreeParticleAverageSite*& destroy) {
    OpenMM_ThreeParticleAverageSite_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_THREEPARTICLEAVERAGESITE_DESTROY(OpenMM_ThreeParticleAverageSite*& destroy) {
    OpenMM_ThreeParticleAverageSite_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT double openmm_threeparticleaveragesite_getweight_(const OpenMM_ThreeParticleAverageSite*& target, int const& particle) {
    return OpenMM_ThreeParticleAverageSite_getWeight(target, particle);
}
OPENMM_EXPORT double OPENMM_THREEPARTICLEAVERAGESITE_GETWEIGHT(const OpenMM_ThreeParticleAverageSite*& target, int const& particle) {
    return OpenMM_ThreeParticleAverageSite_getWeight(target, particle);
}

/* OpenMM::CustomCVForce */
OPENMM_EXPORT void openmm_customcvforce_create_(OpenMM_CustomCVForce*& result, const char* energy, int energy_length) {
    result = OpenMM_CustomCVForce_create(makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMCVFORCE_CREATE(OpenMM_CustomCVForce*& result, const char* energy, int energy_length) {
    result = OpenMM_CustomCVForce_create(makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void openmm_customcvforce_destroy_(OpenMM_CustomCVForce*& destroy) {
    OpenMM_CustomCVForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_CUSTOMCVFORCE_DESTROY(OpenMM_CustomCVForce*& destroy) {
    OpenMM_CustomCVForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT int openmm_customcvforce_getnumcollectivevariables_(const OpenMM_CustomCVForce*& target) {
    return OpenMM_CustomCVForce_getNumCollectiveVariables(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMCVFORCE_GETNUMCOLLECTIVEVARIABLES(const OpenMM_CustomCVForce*& target) {
    return OpenMM_CustomCVForce_getNumCollectiveVariables(target);
}
OPENMM_EXPORT int openmm_customcvforce_getnumglobalparameters_(const OpenMM_CustomCVForce*& target) {
    return OpenMM_CustomCVForce_getNumGlobalParameters(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMCVFORCE_GETNUMGLOBALPARAMETERS(const OpenMM_CustomCVForce*& target) {
    return OpenMM_CustomCVForce_getNumGlobalParameters(target);
}
OPENMM_EXPORT int openmm_customcvforce_getnumenergyparameterderivatives_(const OpenMM_CustomCVForce*& target) {
    return OpenMM_CustomCVForce_getNumEnergyParameterDerivatives(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMCVFORCE_GETNUMENERGYPARAMETERDERIVATIVES(const OpenMM_CustomCVForce*& target) {
    return OpenMM_CustomCVForce_getNumEnergyParameterDerivatives(target);
}
OPENMM_EXPORT int openmm_customcvforce_getnumtabulatedfunctions_(const OpenMM_CustomCVForce*& target) {
    return OpenMM_CustomCVForce_getNumTabulatedFunctions(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMCVFORCE_GETNUMTABULATEDFUNCTIONS(const OpenMM_CustomCVForce*& target) {
    return OpenMM_CustomCVForce_getNumTabulatedFunctions(target);
}
OPENMM_EXPORT void openmm_customcvforce_getenergyfunction_(const OpenMM_CustomCVForce*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCVForce_getEnergyFunction(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMCVFORCE_GETENERGYFUNCTION(const OpenMM_CustomCVForce*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCVForce_getEnergyFunction(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customcvforce_setenergyfunction_(OpenMM_CustomCVForce*& target, const char* energy, int energy_length) {
    OpenMM_CustomCVForce_setEnergyFunction(target, makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMCVFORCE_SETENERGYFUNCTION(OpenMM_CustomCVForce*& target, const char* energy, int energy_length) {
    OpenMM_CustomCVForce_setEnergyFunction(target, makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT int openmm_customcvforce_addcollectivevariable_(OpenMM_CustomCVForce*& target, const char* name, OpenMM_Force*& variable, int name_length) {
    return OpenMM_CustomCVForce_addCollectiveVariable(target, makeString(name, name_length).c_str(), variable);
}
OPENMM_EXPORT int OPENMM_CUSTOMCVFORCE_ADDCOLLECTIVEVARIABLE(OpenMM_CustomCVForce*& target, const char* name, OpenMM_Force*& variable, int name_length) {
    return OpenMM_CustomCVForce_addCollectiveVariable(target, makeString(name, name_length).c_str(), variable);
}
OPENMM_EXPORT void openmm_customcvforce_getcollectivevariablename_(const OpenMM_CustomCVForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCVForce_getCollectiveVariableName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMCVFORCE_GETCOLLECTIVEVARIABLENAME(const OpenMM_CustomCVForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCVForce_getCollectiveVariableName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customcvforce_getcollectivevariable_(OpenMM_CustomCVForce*& target, int const& index, OpenMM_Force*& result) {
    result = OpenMM_CustomCVForce_getCollectiveVariable(target, index);
}
OPENMM_EXPORT void OPENMM_CUSTOMCVFORCE_GETCOLLECTIVEVARIABLE(OpenMM_CustomCVForce*& target, int const& index, OpenMM_Force*& result) {
    result = OpenMM_CustomCVForce_getCollectiveVariable(target, index);
}
OPENMM_EXPORT int openmm_customcvforce_addglobalparameter_(OpenMM_CustomCVForce*& target, const char* name, double const& defaultValue, int name_length) {
    return OpenMM_CustomCVForce_addGlobalParameter(target, makeString(name, name_length).c_str(), defaultValue);
}
OPENMM_EXPORT int OPENMM_CUSTOMCVFORCE_ADDGLOBALPARAMETER(OpenMM_CustomCVForce*& target, const char* name, double const& defaultValue, int name_length) {
    return OpenMM_CustomCVForce_addGlobalParameter(target, makeString(name, name_length).c_str(), defaultValue);
}
OPENMM_EXPORT void openmm_customcvforce_getglobalparametername_(const OpenMM_CustomCVForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCVForce_getGlobalParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMCVFORCE_GETGLOBALPARAMETERNAME(const OpenMM_CustomCVForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCVForce_getGlobalParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customcvforce_setglobalparametername_(OpenMM_CustomCVForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomCVForce_setGlobalParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMCVFORCE_SETGLOBALPARAMETERNAME(OpenMM_CustomCVForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomCVForce_setGlobalParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT double openmm_customcvforce_getglobalparameterdefaultvalue_(const OpenMM_CustomCVForce*& target, int const& index) {
    return OpenMM_CustomCVForce_getGlobalParameterDefaultValue(target, index);
}
OPENMM_EXPORT double OPENMM_CUSTOMCVFORCE_GETGLOBALPARAMETERDEFAULTVALUE(const OpenMM_CustomCVForce*& target, int const& index) {
    return OpenMM_CustomCVForce_getGlobalParameterDefaultValue(target, index);
}
OPENMM_EXPORT void openmm_customcvforce_setglobalparameterdefaultvalue_(OpenMM_CustomCVForce*& target, int const& index, double const& defaultValue) {
    OpenMM_CustomCVForce_setGlobalParameterDefaultValue(target, index, defaultValue);
}
OPENMM_EXPORT void OPENMM_CUSTOMCVFORCE_SETGLOBALPARAMETERDEFAULTVALUE(OpenMM_CustomCVForce*& target, int const& index, double const& defaultValue) {
    OpenMM_CustomCVForce_setGlobalParameterDefaultValue(target, index, defaultValue);
}
OPENMM_EXPORT void openmm_customcvforce_addenergyparameterderivative_(OpenMM_CustomCVForce*& target, const char* name, int name_length) {
    OpenMM_CustomCVForce_addEnergyParameterDerivative(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMCVFORCE_ADDENERGYPARAMETERDERIVATIVE(OpenMM_CustomCVForce*& target, const char* name, int name_length) {
    OpenMM_CustomCVForce_addEnergyParameterDerivative(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void openmm_customcvforce_getenergyparameterderivativename_(const OpenMM_CustomCVForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCVForce_getEnergyParameterDerivativeName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMCVFORCE_GETENERGYPARAMETERDERIVATIVENAME(const OpenMM_CustomCVForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCVForce_getEnergyParameterDerivativeName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT int openmm_customcvforce_addtabulatedfunction_(OpenMM_CustomCVForce*& target, const char* name, OpenMM_TabulatedFunction*& function, int name_length) {
    return OpenMM_CustomCVForce_addTabulatedFunction(target, makeString(name, name_length).c_str(), function);
}
OPENMM_EXPORT int OPENMM_CUSTOMCVFORCE_ADDTABULATEDFUNCTION(OpenMM_CustomCVForce*& target, const char* name, OpenMM_TabulatedFunction*& function, int name_length) {
    return OpenMM_CustomCVForce_addTabulatedFunction(target, makeString(name, name_length).c_str(), function);
}
OPENMM_EXPORT void openmm_customcvforce_gettabulatedfunction_(OpenMM_CustomCVForce*& target, int const& index, OpenMM_TabulatedFunction*& result) {
    result = OpenMM_CustomCVForce_getTabulatedFunction(target, index);
}
OPENMM_EXPORT void OPENMM_CUSTOMCVFORCE_GETTABULATEDFUNCTION(OpenMM_CustomCVForce*& target, int const& index, OpenMM_TabulatedFunction*& result) {
    result = OpenMM_CustomCVForce_getTabulatedFunction(target, index);
}
OPENMM_EXPORT void openmm_customcvforce_gettabulatedfunctionname_(const OpenMM_CustomCVForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCVForce_getTabulatedFunctionName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMCVFORCE_GETTABULATEDFUNCTIONNAME(const OpenMM_CustomCVForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCVForce_getTabulatedFunctionName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customcvforce_getcollectivevariablevalues_(OpenMM_CustomCVForce*& target, OpenMM_Context*& context, OpenMM_DoubleArray*& values) {
    OpenMM_CustomCVForce_getCollectiveVariableValues(target, context, values);
}
OPENMM_EXPORT void OPENMM_CUSTOMCVFORCE_GETCOLLECTIVEVARIABLEVALUES(OpenMM_CustomCVForce*& target, OpenMM_Context*& context, OpenMM_DoubleArray*& values) {
    OpenMM_CustomCVForce_getCollectiveVariableValues(target, context, values);
}
OPENMM_EXPORT void openmm_customcvforce_getinnercontext_(OpenMM_CustomCVForce*& target, OpenMM_Context*& context, OpenMM_Context*& result) {
    result = OpenMM_CustomCVForce_getInnerContext(target, context);
}
OPENMM_EXPORT void OPENMM_CUSTOMCVFORCE_GETINNERCONTEXT(OpenMM_CustomCVForce*& target, OpenMM_Context*& context, OpenMM_Context*& result) {
    result = OpenMM_CustomCVForce_getInnerContext(target, context);
}
OPENMM_EXPORT void openmm_customcvforce_updateparametersincontext_(OpenMM_CustomCVForce*& target, OpenMM_Context*& context) {
    OpenMM_CustomCVForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void OPENMM_CUSTOMCVFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_CustomCVForce*& target, OpenMM_Context*& context) {
    OpenMM_CustomCVForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void openmm_customcvforce_usesperiodicboundaryconditions_(const OpenMM_CustomCVForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CustomCVForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_CUSTOMCVFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_CustomCVForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CustomCVForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::CustomGBForce */
OPENMM_EXPORT void openmm_customgbforce_create_(OpenMM_CustomGBForce*& result) {
    result = OpenMM_CustomGBForce_create();
}
OPENMM_EXPORT void OPENMM_CUSTOMGBFORCE_CREATE(OpenMM_CustomGBForce*& result) {
    result = OpenMM_CustomGBForce_create();
}
OPENMM_EXPORT void openmm_customgbforce_destroy_(OpenMM_CustomGBForce*& destroy) {
    OpenMM_CustomGBForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_CUSTOMGBFORCE_DESTROY(OpenMM_CustomGBForce*& destroy) {
    OpenMM_CustomGBForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT int openmm_customgbforce_getnumparticles_(const OpenMM_CustomGBForce*& target) {
    return OpenMM_CustomGBForce_getNumParticles(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMGBFORCE_GETNUMPARTICLES(const OpenMM_CustomGBForce*& target) {
    return OpenMM_CustomGBForce_getNumParticles(target);
}
OPENMM_EXPORT int openmm_customgbforce_getnumexclusions_(const OpenMM_CustomGBForce*& target) {
    return OpenMM_CustomGBForce_getNumExclusions(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMGBFORCE_GETNUMEXCLUSIONS(const OpenMM_CustomGBForce*& target) {
    return OpenMM_CustomGBForce_getNumExclusions(target);
}
OPENMM_EXPORT int openmm_customgbforce_getnumperparticleparameters_(const OpenMM_CustomGBForce*& target) {
    return OpenMM_CustomGBForce_getNumPerParticleParameters(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMGBFORCE_GETNUMPERPARTICLEPARAMETERS(const OpenMM_CustomGBForce*& target) {
    return OpenMM_CustomGBForce_getNumPerParticleParameters(target);
}
OPENMM_EXPORT int openmm_customgbforce_getnumglobalparameters_(const OpenMM_CustomGBForce*& target) {
    return OpenMM_CustomGBForce_getNumGlobalParameters(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMGBFORCE_GETNUMGLOBALPARAMETERS(const OpenMM_CustomGBForce*& target) {
    return OpenMM_CustomGBForce_getNumGlobalParameters(target);
}
OPENMM_EXPORT int openmm_customgbforce_getnumenergyparameterderivatives_(const OpenMM_CustomGBForce*& target) {
    return OpenMM_CustomGBForce_getNumEnergyParameterDerivatives(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMGBFORCE_GETNUMENERGYPARAMETERDERIVATIVES(const OpenMM_CustomGBForce*& target) {
    return OpenMM_CustomGBForce_getNumEnergyParameterDerivatives(target);
}
OPENMM_EXPORT int openmm_customgbforce_getnumtabulatedfunctions_(const OpenMM_CustomGBForce*& target) {
    return OpenMM_CustomGBForce_getNumTabulatedFunctions(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMGBFORCE_GETNUMTABULATEDFUNCTIONS(const OpenMM_CustomGBForce*& target) {
    return OpenMM_CustomGBForce_getNumTabulatedFunctions(target);
}
OPENMM_EXPORT int openmm_customgbforce_getnumfunctions_(const OpenMM_CustomGBForce*& target) {
    return OpenMM_CustomGBForce_getNumFunctions(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMGBFORCE_GETNUMFUNCTIONS(const OpenMM_CustomGBForce*& target) {
    return OpenMM_CustomGBForce_getNumFunctions(target);
}
OPENMM_EXPORT int openmm_customgbforce_getnumcomputedvalues_(const OpenMM_CustomGBForce*& target) {
    return OpenMM_CustomGBForce_getNumComputedValues(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMGBFORCE_GETNUMCOMPUTEDVALUES(const OpenMM_CustomGBForce*& target) {
    return OpenMM_CustomGBForce_getNumComputedValues(target);
}
OPENMM_EXPORT int openmm_customgbforce_getnumenergyterms_(const OpenMM_CustomGBForce*& target) {
    return OpenMM_CustomGBForce_getNumEnergyTerms(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMGBFORCE_GETNUMENERGYTERMS(const OpenMM_CustomGBForce*& target) {
    return OpenMM_CustomGBForce_getNumEnergyTerms(target);
}
OPENMM_EXPORT void openmm_customgbforce_getnonbondedmethod_(const OpenMM_CustomGBForce*& target, OpenMM_CustomGBForce_NonbondedMethod& result) {
    result = OpenMM_CustomGBForce_getNonbondedMethod(target);
}
OPENMM_EXPORT void OPENMM_CUSTOMGBFORCE_GETNONBONDEDMETHOD(const OpenMM_CustomGBForce*& target, OpenMM_CustomGBForce_NonbondedMethod& result) {
    result = OpenMM_CustomGBForce_getNonbondedMethod(target);
}
OPENMM_EXPORT void openmm_customgbforce_setnonbondedmethod_(OpenMM_CustomGBForce*& target, OpenMM_CustomGBForce_NonbondedMethod& method) {
    OpenMM_CustomGBForce_setNonbondedMethod(target, method);
}
OPENMM_EXPORT void OPENMM_CUSTOMGBFORCE_SETNONBONDEDMETHOD(OpenMM_CustomGBForce*& target, OpenMM_CustomGBForce_NonbondedMethod& method) {
    OpenMM_CustomGBForce_setNonbondedMethod(target, method);
}
OPENMM_EXPORT double openmm_customgbforce_getcutoffdistance_(const OpenMM_CustomGBForce*& target) {
    return OpenMM_CustomGBForce_getCutoffDistance(target);
}
OPENMM_EXPORT double OPENMM_CUSTOMGBFORCE_GETCUTOFFDISTANCE(const OpenMM_CustomGBForce*& target) {
    return OpenMM_CustomGBForce_getCutoffDistance(target);
}
OPENMM_EXPORT void openmm_customgbforce_setcutoffdistance_(OpenMM_CustomGBForce*& target, double const& distance) {
    OpenMM_CustomGBForce_setCutoffDistance(target, distance);
}
OPENMM_EXPORT void OPENMM_CUSTOMGBFORCE_SETCUTOFFDISTANCE(OpenMM_CustomGBForce*& target, double const& distance) {
    OpenMM_CustomGBForce_setCutoffDistance(target, distance);
}
OPENMM_EXPORT int openmm_customgbforce_addperparticleparameter_(OpenMM_CustomGBForce*& target, const char* name, int name_length) {
    return OpenMM_CustomGBForce_addPerParticleParameter(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT int OPENMM_CUSTOMGBFORCE_ADDPERPARTICLEPARAMETER(OpenMM_CustomGBForce*& target, const char* name, int name_length) {
    return OpenMM_CustomGBForce_addPerParticleParameter(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void openmm_customgbforce_getperparticleparametername_(const OpenMM_CustomGBForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomGBForce_getPerParticleParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMGBFORCE_GETPERPARTICLEPARAMETERNAME(const OpenMM_CustomGBForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomGBForce_getPerParticleParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customgbforce_setperparticleparametername_(OpenMM_CustomGBForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomGBForce_setPerParticleParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMGBFORCE_SETPERPARTICLEPARAMETERNAME(OpenMM_CustomGBForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomGBForce_setPerParticleParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT int openmm_customgbforce_addglobalparameter_(OpenMM_CustomGBForce*& target, const char* name, double const& defaultValue, int name_length) {
    return OpenMM_CustomGBForce_addGlobalParameter(target, makeString(name, name_length).c_str(), defaultValue);
}
OPENMM_EXPORT int OPENMM_CUSTOMGBFORCE_ADDGLOBALPARAMETER(OpenMM_CustomGBForce*& target, const char* name, double const& defaultValue, int name_length) {
    return OpenMM_CustomGBForce_addGlobalParameter(target, makeString(name, name_length).c_str(), defaultValue);
}
OPENMM_EXPORT void openmm_customgbforce_getglobalparametername_(const OpenMM_CustomGBForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomGBForce_getGlobalParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMGBFORCE_GETGLOBALPARAMETERNAME(const OpenMM_CustomGBForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomGBForce_getGlobalParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customgbforce_setglobalparametername_(OpenMM_CustomGBForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomGBForce_setGlobalParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMGBFORCE_SETGLOBALPARAMETERNAME(OpenMM_CustomGBForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomGBForce_setGlobalParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT double openmm_customgbforce_getglobalparameterdefaultvalue_(const OpenMM_CustomGBForce*& target, int const& index) {
    return OpenMM_CustomGBForce_getGlobalParameterDefaultValue(target, index);
}
OPENMM_EXPORT double OPENMM_CUSTOMGBFORCE_GETGLOBALPARAMETERDEFAULTVALUE(const OpenMM_CustomGBForce*& target, int const& index) {
    return OpenMM_CustomGBForce_getGlobalParameterDefaultValue(target, index);
}
OPENMM_EXPORT void openmm_customgbforce_setglobalparameterdefaultvalue_(OpenMM_CustomGBForce*& target, int const& index, double const& defaultValue) {
    OpenMM_CustomGBForce_setGlobalParameterDefaultValue(target, index, defaultValue);
}
OPENMM_EXPORT void OPENMM_CUSTOMGBFORCE_SETGLOBALPARAMETERDEFAULTVALUE(OpenMM_CustomGBForce*& target, int const& index, double const& defaultValue) {
    OpenMM_CustomGBForce_setGlobalParameterDefaultValue(target, index, defaultValue);
}
OPENMM_EXPORT void openmm_customgbforce_addenergyparameterderivative_(OpenMM_CustomGBForce*& target, const char* name, int name_length) {
    OpenMM_CustomGBForce_addEnergyParameterDerivative(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMGBFORCE_ADDENERGYPARAMETERDERIVATIVE(OpenMM_CustomGBForce*& target, const char* name, int name_length) {
    OpenMM_CustomGBForce_addEnergyParameterDerivative(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void openmm_customgbforce_getenergyparameterderivativename_(const OpenMM_CustomGBForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomGBForce_getEnergyParameterDerivativeName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMGBFORCE_GETENERGYPARAMETERDERIVATIVENAME(const OpenMM_CustomGBForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomGBForce_getEnergyParameterDerivativeName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT int openmm_customgbforce_addparticle_(OpenMM_CustomGBForce*& target, const OpenMM_DoubleArray*& parameters) {
    return OpenMM_CustomGBForce_addParticle(target, parameters);
}
OPENMM_EXPORT int OPENMM_CUSTOMGBFORCE_ADDPARTICLE(OpenMM_CustomGBForce*& target, const OpenMM_DoubleArray*& parameters) {
    return OpenMM_CustomGBForce_addParticle(target, parameters);
}
OPENMM_EXPORT void openmm_customgbforce_getparticleparameters_(const OpenMM_CustomGBForce*& target, int const& index, OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomGBForce_getParticleParameters(target, index, parameters);
}
OPENMM_EXPORT void OPENMM_CUSTOMGBFORCE_GETPARTICLEPARAMETERS(const OpenMM_CustomGBForce*& target, int const& index, OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomGBForce_getParticleParameters(target, index, parameters);
}
OPENMM_EXPORT void openmm_customgbforce_setparticleparameters_(OpenMM_CustomGBForce*& target, int const& index, const OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomGBForce_setParticleParameters(target, index, parameters);
}
OPENMM_EXPORT void OPENMM_CUSTOMGBFORCE_SETPARTICLEPARAMETERS(OpenMM_CustomGBForce*& target, int const& index, const OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomGBForce_setParticleParameters(target, index, parameters);
}
OPENMM_EXPORT int openmm_customgbforce_addcomputedvalue_(OpenMM_CustomGBForce*& target, const char* name, const char* expression, OpenMM_CustomGBForce_ComputationType& type, int name_length, int expression_length) {
    return OpenMM_CustomGBForce_addComputedValue(target, makeString(name, name_length).c_str(), makeString(expression, expression_length).c_str(), type);
}
OPENMM_EXPORT int OPENMM_CUSTOMGBFORCE_ADDCOMPUTEDVALUE(OpenMM_CustomGBForce*& target, const char* name, const char* expression, OpenMM_CustomGBForce_ComputationType& type, int name_length, int expression_length) {
    return OpenMM_CustomGBForce_addComputedValue(target, makeString(name, name_length).c_str(), makeString(expression, expression_length).c_str(), type);
}
OPENMM_EXPORT void openmm_customgbforce_getcomputedvalueparameters_(const OpenMM_CustomGBForce*& target, int const& index, char** name, char** expression, OpenMM_CustomGBForce_ComputationType*& type) {
    OpenMM_CustomGBForce_getComputedValueParameters(target, index, name, expression, type);
}
OPENMM_EXPORT void OPENMM_CUSTOMGBFORCE_GETCOMPUTEDVALUEPARAMETERS(const OpenMM_CustomGBForce*& target, int const& index, char** name, char** expression, OpenMM_CustomGBForce_ComputationType*& type) {
    OpenMM_CustomGBForce_getComputedValueParameters(target, index, name, expression, type);
}
OPENMM_EXPORT void openmm_customgbforce_setcomputedvalueparameters_(OpenMM_CustomGBForce*& target, int const& index, const char* name, const char* expression, OpenMM_CustomGBForce_ComputationType& type, int name_length, int expression_length) {
    OpenMM_CustomGBForce_setComputedValueParameters(target, index, makeString(name, name_length).c_str(), makeString(expression, expression_length).c_str(), type);
}
OPENMM_EXPORT void OPENMM_CUSTOMGBFORCE_SETCOMPUTEDVALUEPARAMETERS(OpenMM_CustomGBForce*& target, int const& index, const char* name, const char* expression, OpenMM_CustomGBForce_ComputationType& type, int name_length, int expression_length) {
    OpenMM_CustomGBForce_setComputedValueParameters(target, index, makeString(name, name_length).c_str(), makeString(expression, expression_length).c_str(), type);
}
OPENMM_EXPORT int openmm_customgbforce_addenergyterm_(OpenMM_CustomGBForce*& target, const char* expression, OpenMM_CustomGBForce_ComputationType& type, int expression_length) {
    return OpenMM_CustomGBForce_addEnergyTerm(target, makeString(expression, expression_length).c_str(), type);
}
OPENMM_EXPORT int OPENMM_CUSTOMGBFORCE_ADDENERGYTERM(OpenMM_CustomGBForce*& target, const char* expression, OpenMM_CustomGBForce_ComputationType& type, int expression_length) {
    return OpenMM_CustomGBForce_addEnergyTerm(target, makeString(expression, expression_length).c_str(), type);
}
OPENMM_EXPORT void openmm_customgbforce_getenergytermparameters_(const OpenMM_CustomGBForce*& target, int const& index, char** expression, OpenMM_CustomGBForce_ComputationType*& type) {
    OpenMM_CustomGBForce_getEnergyTermParameters(target, index, expression, type);
}
OPENMM_EXPORT void OPENMM_CUSTOMGBFORCE_GETENERGYTERMPARAMETERS(const OpenMM_CustomGBForce*& target, int const& index, char** expression, OpenMM_CustomGBForce_ComputationType*& type) {
    OpenMM_CustomGBForce_getEnergyTermParameters(target, index, expression, type);
}
OPENMM_EXPORT void openmm_customgbforce_setenergytermparameters_(OpenMM_CustomGBForce*& target, int const& index, const char* expression, OpenMM_CustomGBForce_ComputationType& type, int expression_length) {
    OpenMM_CustomGBForce_setEnergyTermParameters(target, index, makeString(expression, expression_length).c_str(), type);
}
OPENMM_EXPORT void OPENMM_CUSTOMGBFORCE_SETENERGYTERMPARAMETERS(OpenMM_CustomGBForce*& target, int const& index, const char* expression, OpenMM_CustomGBForce_ComputationType& type, int expression_length) {
    OpenMM_CustomGBForce_setEnergyTermParameters(target, index, makeString(expression, expression_length).c_str(), type);
}
OPENMM_EXPORT int openmm_customgbforce_addexclusion_(OpenMM_CustomGBForce*& target, int const& particle1, int const& particle2) {
    return OpenMM_CustomGBForce_addExclusion(target, particle1, particle2);
}
OPENMM_EXPORT int OPENMM_CUSTOMGBFORCE_ADDEXCLUSION(OpenMM_CustomGBForce*& target, int const& particle1, int const& particle2) {
    return OpenMM_CustomGBForce_addExclusion(target, particle1, particle2);
}
OPENMM_EXPORT void openmm_customgbforce_getexclusionparticles_(const OpenMM_CustomGBForce*& target, int const& index, int* particle1, int* particle2) {
    OpenMM_CustomGBForce_getExclusionParticles(target, index, particle1, particle2);
}
OPENMM_EXPORT void OPENMM_CUSTOMGBFORCE_GETEXCLUSIONPARTICLES(const OpenMM_CustomGBForce*& target, int const& index, int* particle1, int* particle2) {
    OpenMM_CustomGBForce_getExclusionParticles(target, index, particle1, particle2);
}
OPENMM_EXPORT void openmm_customgbforce_setexclusionparticles_(OpenMM_CustomGBForce*& target, int const& index, int const& particle1, int const& particle2) {
    OpenMM_CustomGBForce_setExclusionParticles(target, index, particle1, particle2);
}
OPENMM_EXPORT void OPENMM_CUSTOMGBFORCE_SETEXCLUSIONPARTICLES(OpenMM_CustomGBForce*& target, int const& index, int const& particle1, int const& particle2) {
    OpenMM_CustomGBForce_setExclusionParticles(target, index, particle1, particle2);
}
OPENMM_EXPORT int openmm_customgbforce_addtabulatedfunction_(OpenMM_CustomGBForce*& target, const char* name, OpenMM_TabulatedFunction*& function, int name_length) {
    return OpenMM_CustomGBForce_addTabulatedFunction(target, makeString(name, name_length).c_str(), function);
}
OPENMM_EXPORT int OPENMM_CUSTOMGBFORCE_ADDTABULATEDFUNCTION(OpenMM_CustomGBForce*& target, const char* name, OpenMM_TabulatedFunction*& function, int name_length) {
    return OpenMM_CustomGBForce_addTabulatedFunction(target, makeString(name, name_length).c_str(), function);
}
OPENMM_EXPORT void openmm_customgbforce_gettabulatedfunction_(OpenMM_CustomGBForce*& target, int const& index, OpenMM_TabulatedFunction*& result) {
    result = OpenMM_CustomGBForce_getTabulatedFunction(target, index);
}
OPENMM_EXPORT void OPENMM_CUSTOMGBFORCE_GETTABULATEDFUNCTION(OpenMM_CustomGBForce*& target, int const& index, OpenMM_TabulatedFunction*& result) {
    result = OpenMM_CustomGBForce_getTabulatedFunction(target, index);
}
OPENMM_EXPORT void openmm_customgbforce_gettabulatedfunctionname_(const OpenMM_CustomGBForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomGBForce_getTabulatedFunctionName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMGBFORCE_GETTABULATEDFUNCTIONNAME(const OpenMM_CustomGBForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomGBForce_getTabulatedFunctionName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT int openmm_customgbforce_addfunction_(OpenMM_CustomGBForce*& target, const char* name, const OpenMM_DoubleArray*& values, double const& min, double const& max, int name_length) {
    return OpenMM_CustomGBForce_addFunction(target, makeString(name, name_length).c_str(), values, min, max);
}
OPENMM_EXPORT int OPENMM_CUSTOMGBFORCE_ADDFUNCTION(OpenMM_CustomGBForce*& target, const char* name, const OpenMM_DoubleArray*& values, double const& min, double const& max, int name_length) {
    return OpenMM_CustomGBForce_addFunction(target, makeString(name, name_length).c_str(), values, min, max);
}
OPENMM_EXPORT void openmm_customgbforce_getfunctionparameters_(const OpenMM_CustomGBForce*& target, int const& index, char** name, OpenMM_DoubleArray*& values, double* min, double* max) {
    OpenMM_CustomGBForce_getFunctionParameters(target, index, name, values, min, max);
}
OPENMM_EXPORT void OPENMM_CUSTOMGBFORCE_GETFUNCTIONPARAMETERS(const OpenMM_CustomGBForce*& target, int const& index, char** name, OpenMM_DoubleArray*& values, double* min, double* max) {
    OpenMM_CustomGBForce_getFunctionParameters(target, index, name, values, min, max);
}
OPENMM_EXPORT void openmm_customgbforce_setfunctionparameters_(OpenMM_CustomGBForce*& target, int const& index, const char* name, const OpenMM_DoubleArray*& values, double const& min, double const& max, int name_length) {
    OpenMM_CustomGBForce_setFunctionParameters(target, index, makeString(name, name_length).c_str(), values, min, max);
}
OPENMM_EXPORT void OPENMM_CUSTOMGBFORCE_SETFUNCTIONPARAMETERS(OpenMM_CustomGBForce*& target, int const& index, const char* name, const OpenMM_DoubleArray*& values, double const& min, double const& max, int name_length) {
    OpenMM_CustomGBForce_setFunctionParameters(target, index, makeString(name, name_length).c_str(), values, min, max);
}
OPENMM_EXPORT void openmm_customgbforce_updateparametersincontext_(OpenMM_CustomGBForce*& target, OpenMM_Context*& context) {
    OpenMM_CustomGBForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void OPENMM_CUSTOMGBFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_CustomGBForce*& target, OpenMM_Context*& context) {
    OpenMM_CustomGBForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void openmm_customgbforce_usesperiodicboundaryconditions_(const OpenMM_CustomGBForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CustomGBForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_CUSTOMGBFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_CustomGBForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CustomGBForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::Discrete3DFunction */
OPENMM_EXPORT void openmm_discrete3dfunction_create_(OpenMM_Discrete3DFunction*& result, int const& xsize, int const& ysize, int const& zsize, const OpenMM_DoubleArray*& values) {
    result = OpenMM_Discrete3DFunction_create(xsize, ysize, zsize, values);
}
OPENMM_EXPORT void OPENMM_DISCRETE3DFUNCTION_CREATE(OpenMM_Discrete3DFunction*& result, int const& xsize, int const& ysize, int const& zsize, const OpenMM_DoubleArray*& values) {
    result = OpenMM_Discrete3DFunction_create(xsize, ysize, zsize, values);
}
OPENMM_EXPORT void openmm_discrete3dfunction_destroy_(OpenMM_Discrete3DFunction*& destroy) {
    OpenMM_Discrete3DFunction_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_DISCRETE3DFUNCTION_DESTROY(OpenMM_Discrete3DFunction*& destroy) {
    OpenMM_Discrete3DFunction_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void openmm_discrete3dfunction_getfunctionparameters_(const OpenMM_Discrete3DFunction*& target, int* xsize, int* ysize, int* zsize, OpenMM_DoubleArray*& values) {
    OpenMM_Discrete3DFunction_getFunctionParameters(target, xsize, ysize, zsize, values);
}
OPENMM_EXPORT void OPENMM_DISCRETE3DFUNCTION_GETFUNCTIONPARAMETERS(const OpenMM_Discrete3DFunction*& target, int* xsize, int* ysize, int* zsize, OpenMM_DoubleArray*& values) {
    OpenMM_Discrete3DFunction_getFunctionParameters(target, xsize, ysize, zsize, values);
}
OPENMM_EXPORT void openmm_discrete3dfunction_setfunctionparameters_(OpenMM_Discrete3DFunction*& target, int const& xsize, int const& ysize, int const& zsize, const OpenMM_DoubleArray*& values) {
    OpenMM_Discrete3DFunction_setFunctionParameters(target, xsize, ysize, zsize, values);
}
OPENMM_EXPORT void OPENMM_DISCRETE3DFUNCTION_SETFUNCTIONPARAMETERS(OpenMM_Discrete3DFunction*& target, int const& xsize, int const& ysize, int const& zsize, const OpenMM_DoubleArray*& values) {
    OpenMM_Discrete3DFunction_setFunctionParameters(target, xsize, ysize, zsize, values);
}
OPENMM_EXPORT void openmm_discrete3dfunction_copy_(const OpenMM_Discrete3DFunction*& target, OpenMM_Discrete3DFunction*& result) {
    result = OpenMM_Discrete3DFunction_Copy(target);
}
OPENMM_EXPORT void OPENMM_DISCRETE3DFUNCTION_COPY(const OpenMM_Discrete3DFunction*& target, OpenMM_Discrete3DFunction*& result) {
    result = OpenMM_Discrete3DFunction_Copy(target);
}

/* OpenMM::PeriodicTorsionForce */
OPENMM_EXPORT void openmm_periodictorsionforce_create_(OpenMM_PeriodicTorsionForce*& result) {
    result = OpenMM_PeriodicTorsionForce_create();
}
OPENMM_EXPORT void OPENMM_PERIODICTORSIONFORCE_CREATE(OpenMM_PeriodicTorsionForce*& result) {
    result = OpenMM_PeriodicTorsionForce_create();
}
OPENMM_EXPORT void openmm_periodictorsionforce_destroy_(OpenMM_PeriodicTorsionForce*& destroy) {
    OpenMM_PeriodicTorsionForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_PERIODICTORSIONFORCE_DESTROY(OpenMM_PeriodicTorsionForce*& destroy) {
    OpenMM_PeriodicTorsionForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT int openmm_periodictorsionforce_getnumtorsions_(const OpenMM_PeriodicTorsionForce*& target) {
    return OpenMM_PeriodicTorsionForce_getNumTorsions(target);
}
OPENMM_EXPORT int OPENMM_PERIODICTORSIONFORCE_GETNUMTORSIONS(const OpenMM_PeriodicTorsionForce*& target) {
    return OpenMM_PeriodicTorsionForce_getNumTorsions(target);
}
OPENMM_EXPORT int openmm_periodictorsionforce_addtorsion_(OpenMM_PeriodicTorsionForce*& target, int const& particle1, int const& particle2, int const& particle3, int const& particle4, int const& periodicity, double const& phase, double const& k) {
    return OpenMM_PeriodicTorsionForce_addTorsion(target, particle1, particle2, particle3, particle4, periodicity, phase, k);
}
OPENMM_EXPORT int OPENMM_PERIODICTORSIONFORCE_ADDTORSION(OpenMM_PeriodicTorsionForce*& target, int const& particle1, int const& particle2, int const& particle3, int const& particle4, int const& periodicity, double const& phase, double const& k) {
    return OpenMM_PeriodicTorsionForce_addTorsion(target, particle1, particle2, particle3, particle4, periodicity, phase, k);
}
OPENMM_EXPORT void openmm_periodictorsionforce_gettorsionparameters_(const OpenMM_PeriodicTorsionForce*& target, int const& index, int* particle1, int* particle2, int* particle3, int* particle4, int* periodicity, double* phase, double* k) {
    OpenMM_PeriodicTorsionForce_getTorsionParameters(target, index, particle1, particle2, particle3, particle4, periodicity, phase, k);
}
OPENMM_EXPORT void OPENMM_PERIODICTORSIONFORCE_GETTORSIONPARAMETERS(const OpenMM_PeriodicTorsionForce*& target, int const& index, int* particle1, int* particle2, int* particle3, int* particle4, int* periodicity, double* phase, double* k) {
    OpenMM_PeriodicTorsionForce_getTorsionParameters(target, index, particle1, particle2, particle3, particle4, periodicity, phase, k);
}
OPENMM_EXPORT void openmm_periodictorsionforce_settorsionparameters_(OpenMM_PeriodicTorsionForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, int const& particle4, int const& periodicity, double const& phase, double const& k) {
    OpenMM_PeriodicTorsionForce_setTorsionParameters(target, index, particle1, particle2, particle3, particle4, periodicity, phase, k);
}
OPENMM_EXPORT void OPENMM_PERIODICTORSIONFORCE_SETTORSIONPARAMETERS(OpenMM_PeriodicTorsionForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, int const& particle4, int const& periodicity, double const& phase, double const& k) {
    OpenMM_PeriodicTorsionForce_setTorsionParameters(target, index, particle1, particle2, particle3, particle4, periodicity, phase, k);
}
OPENMM_EXPORT void openmm_periodictorsionforce_updateparametersincontext_(OpenMM_PeriodicTorsionForce*& target, OpenMM_Context*& context) {
    OpenMM_PeriodicTorsionForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void OPENMM_PERIODICTORSIONFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_PeriodicTorsionForce*& target, OpenMM_Context*& context) {
    OpenMM_PeriodicTorsionForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void openmm_periodictorsionforce_setusesperiodicboundaryconditions_(OpenMM_PeriodicTorsionForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_PeriodicTorsionForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT void OPENMM_PERIODICTORSIONFORCE_SETUSESPERIODICBOUNDARYCONDITIONS(OpenMM_PeriodicTorsionForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_PeriodicTorsionForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT void openmm_periodictorsionforce_usesperiodicboundaryconditions_(const OpenMM_PeriodicTorsionForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_PeriodicTorsionForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_PERIODICTORSIONFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_PeriodicTorsionForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_PeriodicTorsionForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::LocalCoordinatesSite */
OPENMM_EXPORT void openmm_localcoordinatessite_create_(OpenMM_LocalCoordinatesSite*& result, const OpenMM_IntArray*& particles, const OpenMM_DoubleArray*& originWeights, const OpenMM_DoubleArray*& xWeights, const OpenMM_DoubleArray*& yWeights, const OpenMM_Vec3* localPosition) {
    result = OpenMM_LocalCoordinatesSite_create(particles, originWeights, xWeights, yWeights, localPosition);
}
OPENMM_EXPORT void OPENMM_LOCALCOORDINATESSITE_CREATE(OpenMM_LocalCoordinatesSite*& result, const OpenMM_IntArray*& particles, const OpenMM_DoubleArray*& originWeights, const OpenMM_DoubleArray*& xWeights, const OpenMM_DoubleArray*& yWeights, const OpenMM_Vec3* localPosition) {
    result = OpenMM_LocalCoordinatesSite_create(particles, originWeights, xWeights, yWeights, localPosition);
}
OPENMM_EXPORT void openmm_localcoordinatessite_create_2_(OpenMM_LocalCoordinatesSite*& result, int const& particle1, int const& particle2, int const& particle3, const OpenMM_Vec3* originWeights, const OpenMM_Vec3* xWeights, const OpenMM_Vec3* yWeights, const OpenMM_Vec3* localPosition) {
    result = OpenMM_LocalCoordinatesSite_create_2(particle1, particle2, particle3, originWeights, xWeights, yWeights, localPosition);
}
OPENMM_EXPORT void OPENMM_LOCALCOORDINATESSITE_CREATE_2(OpenMM_LocalCoordinatesSite*& result, int const& particle1, int const& particle2, int const& particle3, const OpenMM_Vec3* originWeights, const OpenMM_Vec3* xWeights, const OpenMM_Vec3* yWeights, const OpenMM_Vec3* localPosition) {
    result = OpenMM_LocalCoordinatesSite_create_2(particle1, particle2, particle3, originWeights, xWeights, yWeights, localPosition);
}
OPENMM_EXPORT void openmm_localcoordinatessite_destroy_(OpenMM_LocalCoordinatesSite*& destroy) {
    OpenMM_LocalCoordinatesSite_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_LOCALCOORDINATESSITE_DESTROY(OpenMM_LocalCoordinatesSite*& destroy) {
    OpenMM_LocalCoordinatesSite_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void openmm_localcoordinatessite_getoriginweights_(const OpenMM_LocalCoordinatesSite*& target, OpenMM_DoubleArray*& weights) {
    OpenMM_LocalCoordinatesSite_getOriginWeights(target, weights);
}
OPENMM_EXPORT void OPENMM_LOCALCOORDINATESSITE_GETORIGINWEIGHTS(const OpenMM_LocalCoordinatesSite*& target, OpenMM_DoubleArray*& weights) {
    OpenMM_LocalCoordinatesSite_getOriginWeights(target, weights);
}
OPENMM_EXPORT void openmm_localcoordinatessite_getxweights_(const OpenMM_LocalCoordinatesSite*& target, OpenMM_DoubleArray*& weights) {
    OpenMM_LocalCoordinatesSite_getXWeights(target, weights);
}
OPENMM_EXPORT void OPENMM_LOCALCOORDINATESSITE_GETXWEIGHTS(const OpenMM_LocalCoordinatesSite*& target, OpenMM_DoubleArray*& weights) {
    OpenMM_LocalCoordinatesSite_getXWeights(target, weights);
}
OPENMM_EXPORT void openmm_localcoordinatessite_getyweights_(const OpenMM_LocalCoordinatesSite*& target, OpenMM_DoubleArray*& weights) {
    OpenMM_LocalCoordinatesSite_getYWeights(target, weights);
}
OPENMM_EXPORT void OPENMM_LOCALCOORDINATESSITE_GETYWEIGHTS(const OpenMM_LocalCoordinatesSite*& target, OpenMM_DoubleArray*& weights) {
    OpenMM_LocalCoordinatesSite_getYWeights(target, weights);
}
OPENMM_EXPORT void openmm_localcoordinatessite_getlocalposition_(const OpenMM_LocalCoordinatesSite*& target, const OpenMM_Vec3*& result) {
    result = OpenMM_LocalCoordinatesSite_getLocalPosition(target);
}
OPENMM_EXPORT void OPENMM_LOCALCOORDINATESSITE_GETLOCALPOSITION(const OpenMM_LocalCoordinatesSite*& target, const OpenMM_Vec3*& result) {
    result = OpenMM_LocalCoordinatesSite_getLocalPosition(target);
}

/* OpenMM::State */
OPENMM_EXPORT void openmm_state_create_(OpenMM_State*& result) {
    result = OpenMM_State_create();
}
OPENMM_EXPORT void OPENMM_STATE_CREATE(OpenMM_State*& result) {
    result = OpenMM_State_create();
}
OPENMM_EXPORT void openmm_state_destroy_(OpenMM_State*& destroy) {
    OpenMM_State_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_STATE_DESTROY(OpenMM_State*& destroy) {
    OpenMM_State_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT double openmm_state_gettime_(const OpenMM_State*& target) {
    return OpenMM_State_getTime(target);
}
OPENMM_EXPORT double OPENMM_STATE_GETTIME(const OpenMM_State*& target) {
    return OpenMM_State_getTime(target);
}
OPENMM_EXPORT void openmm_state_getpositions_(const OpenMM_State*& target, const OpenMM_Vec3Array*& result) {
    result = OpenMM_State_getPositions(target);
}
OPENMM_EXPORT void OPENMM_STATE_GETPOSITIONS(const OpenMM_State*& target, const OpenMM_Vec3Array*& result) {
    result = OpenMM_State_getPositions(target);
}
OPENMM_EXPORT void openmm_state_getvelocities_(const OpenMM_State*& target, const OpenMM_Vec3Array*& result) {
    result = OpenMM_State_getVelocities(target);
}
OPENMM_EXPORT void OPENMM_STATE_GETVELOCITIES(const OpenMM_State*& target, const OpenMM_Vec3Array*& result) {
    result = OpenMM_State_getVelocities(target);
}
OPENMM_EXPORT void openmm_state_getforces_(const OpenMM_State*& target, const OpenMM_Vec3Array*& result) {
    result = OpenMM_State_getForces(target);
}
OPENMM_EXPORT void OPENMM_STATE_GETFORCES(const OpenMM_State*& target, const OpenMM_Vec3Array*& result) {
    result = OpenMM_State_getForces(target);
}
OPENMM_EXPORT double openmm_state_getkineticenergy_(const OpenMM_State*& target) {
    return OpenMM_State_getKineticEnergy(target);
}
OPENMM_EXPORT double OPENMM_STATE_GETKINETICENERGY(const OpenMM_State*& target) {
    return OpenMM_State_getKineticEnergy(target);
}
OPENMM_EXPORT double openmm_state_getpotentialenergy_(const OpenMM_State*& target) {
    return OpenMM_State_getPotentialEnergy(target);
}
OPENMM_EXPORT double OPENMM_STATE_GETPOTENTIALENERGY(const OpenMM_State*& target) {
    return OpenMM_State_getPotentialEnergy(target);
}
OPENMM_EXPORT void openmm_state_getperiodicboxvectors_(const OpenMM_State*& target, OpenMM_Vec3* a, OpenMM_Vec3* b, OpenMM_Vec3* c) {
    OpenMM_State_getPeriodicBoxVectors(target, a, b, c);
}
OPENMM_EXPORT void OPENMM_STATE_GETPERIODICBOXVECTORS(const OpenMM_State*& target, OpenMM_Vec3* a, OpenMM_Vec3* b, OpenMM_Vec3* c) {
    OpenMM_State_getPeriodicBoxVectors(target, a, b, c);
}
OPENMM_EXPORT double openmm_state_getperiodicboxvolume_(const OpenMM_State*& target) {
    return OpenMM_State_getPeriodicBoxVolume(target);
}
OPENMM_EXPORT double OPENMM_STATE_GETPERIODICBOXVOLUME(const OpenMM_State*& target) {
    return OpenMM_State_getPeriodicBoxVolume(target);
}
OPENMM_EXPORT void openmm_state_getparameters_(const OpenMM_State*& target, const OpenMM_ParameterArray*& result) {
    result = OpenMM_State_getParameters(target);
}
OPENMM_EXPORT void OPENMM_STATE_GETPARAMETERS(const OpenMM_State*& target, const OpenMM_ParameterArray*& result) {
    result = OpenMM_State_getParameters(target);
}
OPENMM_EXPORT void openmm_state_getenergyparameterderivatives_(const OpenMM_State*& target, const OpenMM_ParameterArray*& result) {
    result = OpenMM_State_getEnergyParameterDerivatives(target);
}
OPENMM_EXPORT void OPENMM_STATE_GETENERGYPARAMETERDERIVATIVES(const OpenMM_State*& target, const OpenMM_ParameterArray*& result) {
    result = OpenMM_State_getEnergyParameterDerivatives(target);
}
OPENMM_EXPORT int openmm_state_getdatatypes_(const OpenMM_State*& target) {
    return OpenMM_State_getDataTypes(target);
}
OPENMM_EXPORT int OPENMM_STATE_GETDATATYPES(const OpenMM_State*& target) {
    return OpenMM_State_getDataTypes(target);
}

/* OpenMM::OpenMMException */
OPENMM_EXPORT void openmm_openmmexception_create_(OpenMM_OpenMMException*& result, const char* message, int message_length) {
    result = OpenMM_OpenMMException_create(makeString(message, message_length).c_str());
}
OPENMM_EXPORT void OPENMM_OPENMMEXCEPTION_CREATE(OpenMM_OpenMMException*& result, const char* message, int message_length) {
    result = OpenMM_OpenMMException_create(makeString(message, message_length).c_str());
}
OPENMM_EXPORT void openmm_openmmexception_destroy_(OpenMM_OpenMMException*& destroy) {
    OpenMM_OpenMMException_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_OPENMMEXCEPTION_DESTROY(OpenMM_OpenMMException*& destroy) {
    OpenMM_OpenMMException_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void openmm_openmmexception_what_(const OpenMM_OpenMMException*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_OpenMMException_what(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_OPENMMEXCEPTION_WHAT(const OpenMM_OpenMMException*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_OpenMMException_what(target);
    copyAndPadString(result, result_chars, result_length);
}

/* OpenMM::CMAPTorsionForce */
OPENMM_EXPORT void openmm_cmaptorsionforce_create_(OpenMM_CMAPTorsionForce*& result) {
    result = OpenMM_CMAPTorsionForce_create();
}
OPENMM_EXPORT void OPENMM_CMAPTORSIONFORCE_CREATE(OpenMM_CMAPTorsionForce*& result) {
    result = OpenMM_CMAPTorsionForce_create();
}
OPENMM_EXPORT void openmm_cmaptorsionforce_destroy_(OpenMM_CMAPTorsionForce*& destroy) {
    OpenMM_CMAPTorsionForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_CMAPTORSIONFORCE_DESTROY(OpenMM_CMAPTorsionForce*& destroy) {
    OpenMM_CMAPTorsionForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT int openmm_cmaptorsionforce_getnummaps_(const OpenMM_CMAPTorsionForce*& target) {
    return OpenMM_CMAPTorsionForce_getNumMaps(target);
}
OPENMM_EXPORT int OPENMM_CMAPTORSIONFORCE_GETNUMMAPS(const OpenMM_CMAPTorsionForce*& target) {
    return OpenMM_CMAPTorsionForce_getNumMaps(target);
}
OPENMM_EXPORT int openmm_cmaptorsionforce_getnumtorsions_(const OpenMM_CMAPTorsionForce*& target) {
    return OpenMM_CMAPTorsionForce_getNumTorsions(target);
}
OPENMM_EXPORT int OPENMM_CMAPTORSIONFORCE_GETNUMTORSIONS(const OpenMM_CMAPTorsionForce*& target) {
    return OpenMM_CMAPTorsionForce_getNumTorsions(target);
}
OPENMM_EXPORT int openmm_cmaptorsionforce_addmap_(OpenMM_CMAPTorsionForce*& target, int const& size, const OpenMM_DoubleArray*& energy) {
    return OpenMM_CMAPTorsionForce_addMap(target, size, energy);
}
OPENMM_EXPORT int OPENMM_CMAPTORSIONFORCE_ADDMAP(OpenMM_CMAPTorsionForce*& target, int const& size, const OpenMM_DoubleArray*& energy) {
    return OpenMM_CMAPTorsionForce_addMap(target, size, energy);
}
OPENMM_EXPORT void openmm_cmaptorsionforce_getmapparameters_(const OpenMM_CMAPTorsionForce*& target, int const& index, int* size, OpenMM_DoubleArray*& energy) {
    OpenMM_CMAPTorsionForce_getMapParameters(target, index, size, energy);
}
OPENMM_EXPORT void OPENMM_CMAPTORSIONFORCE_GETMAPPARAMETERS(const OpenMM_CMAPTorsionForce*& target, int const& index, int* size, OpenMM_DoubleArray*& energy) {
    OpenMM_CMAPTorsionForce_getMapParameters(target, index, size, energy);
}
OPENMM_EXPORT void openmm_cmaptorsionforce_setmapparameters_(OpenMM_CMAPTorsionForce*& target, int const& index, int const& size, const OpenMM_DoubleArray*& energy) {
    OpenMM_CMAPTorsionForce_setMapParameters(target, index, size, energy);
}
OPENMM_EXPORT void OPENMM_CMAPTORSIONFORCE_SETMAPPARAMETERS(OpenMM_CMAPTorsionForce*& target, int const& index, int const& size, const OpenMM_DoubleArray*& energy) {
    OpenMM_CMAPTorsionForce_setMapParameters(target, index, size, energy);
}
OPENMM_EXPORT int openmm_cmaptorsionforce_addtorsion_(OpenMM_CMAPTorsionForce*& target, int const& map, int const& a1, int const& a2, int const& a3, int const& a4, int const& b1, int const& b2, int const& b3, int const& b4) {
    return OpenMM_CMAPTorsionForce_addTorsion(target, map, a1, a2, a3, a4, b1, b2, b3, b4);
}
OPENMM_EXPORT int OPENMM_CMAPTORSIONFORCE_ADDTORSION(OpenMM_CMAPTorsionForce*& target, int const& map, int const& a1, int const& a2, int const& a3, int const& a4, int const& b1, int const& b2, int const& b3, int const& b4) {
    return OpenMM_CMAPTorsionForce_addTorsion(target, map, a1, a2, a3, a4, b1, b2, b3, b4);
}
OPENMM_EXPORT void openmm_cmaptorsionforce_gettorsionparameters_(const OpenMM_CMAPTorsionForce*& target, int const& index, int* map, int* a1, int* a2, int* a3, int* a4, int* b1, int* b2, int* b3, int* b4) {
    OpenMM_CMAPTorsionForce_getTorsionParameters(target, index, map, a1, a2, a3, a4, b1, b2, b3, b4);
}
OPENMM_EXPORT void OPENMM_CMAPTORSIONFORCE_GETTORSIONPARAMETERS(const OpenMM_CMAPTorsionForce*& target, int const& index, int* map, int* a1, int* a2, int* a3, int* a4, int* b1, int* b2, int* b3, int* b4) {
    OpenMM_CMAPTorsionForce_getTorsionParameters(target, index, map, a1, a2, a3, a4, b1, b2, b3, b4);
}
OPENMM_EXPORT void openmm_cmaptorsionforce_settorsionparameters_(OpenMM_CMAPTorsionForce*& target, int const& index, int const& map, int const& a1, int const& a2, int const& a3, int const& a4, int const& b1, int const& b2, int const& b3, int const& b4) {
    OpenMM_CMAPTorsionForce_setTorsionParameters(target, index, map, a1, a2, a3, a4, b1, b2, b3, b4);
}
OPENMM_EXPORT void OPENMM_CMAPTORSIONFORCE_SETTORSIONPARAMETERS(OpenMM_CMAPTorsionForce*& target, int const& index, int const& map, int const& a1, int const& a2, int const& a3, int const& a4, int const& b1, int const& b2, int const& b3, int const& b4) {
    OpenMM_CMAPTorsionForce_setTorsionParameters(target, index, map, a1, a2, a3, a4, b1, b2, b3, b4);
}
OPENMM_EXPORT void openmm_cmaptorsionforce_updateparametersincontext_(OpenMM_CMAPTorsionForce*& target, OpenMM_Context*& context) {
    OpenMM_CMAPTorsionForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void OPENMM_CMAPTORSIONFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_CMAPTorsionForce*& target, OpenMM_Context*& context) {
    OpenMM_CMAPTorsionForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void openmm_cmaptorsionforce_setusesperiodicboundaryconditions_(OpenMM_CMAPTorsionForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_CMAPTorsionForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT void OPENMM_CMAPTORSIONFORCE_SETUSESPERIODICBOUNDARYCONDITIONS(OpenMM_CMAPTorsionForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_CMAPTorsionForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT void openmm_cmaptorsionforce_usesperiodicboundaryconditions_(const OpenMM_CMAPTorsionForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CMAPTorsionForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_CMAPTORSIONFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_CMAPTorsionForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CMAPTorsionForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::Discrete1DFunction */
OPENMM_EXPORT void openmm_discrete1dfunction_create_(OpenMM_Discrete1DFunction*& result, const OpenMM_DoubleArray*& values) {
    result = OpenMM_Discrete1DFunction_create(values);
}
OPENMM_EXPORT void OPENMM_DISCRETE1DFUNCTION_CREATE(OpenMM_Discrete1DFunction*& result, const OpenMM_DoubleArray*& values) {
    result = OpenMM_Discrete1DFunction_create(values);
}
OPENMM_EXPORT void openmm_discrete1dfunction_destroy_(OpenMM_Discrete1DFunction*& destroy) {
    OpenMM_Discrete1DFunction_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_DISCRETE1DFUNCTION_DESTROY(OpenMM_Discrete1DFunction*& destroy) {
    OpenMM_Discrete1DFunction_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void openmm_discrete1dfunction_getfunctionparameters_(const OpenMM_Discrete1DFunction*& target, OpenMM_DoubleArray*& values) {
    OpenMM_Discrete1DFunction_getFunctionParameters(target, values);
}
OPENMM_EXPORT void OPENMM_DISCRETE1DFUNCTION_GETFUNCTIONPARAMETERS(const OpenMM_Discrete1DFunction*& target, OpenMM_DoubleArray*& values) {
    OpenMM_Discrete1DFunction_getFunctionParameters(target, values);
}
OPENMM_EXPORT void openmm_discrete1dfunction_setfunctionparameters_(OpenMM_Discrete1DFunction*& target, const OpenMM_DoubleArray*& values) {
    OpenMM_Discrete1DFunction_setFunctionParameters(target, values);
}
OPENMM_EXPORT void OPENMM_DISCRETE1DFUNCTION_SETFUNCTIONPARAMETERS(OpenMM_Discrete1DFunction*& target, const OpenMM_DoubleArray*& values) {
    OpenMM_Discrete1DFunction_setFunctionParameters(target, values);
}
OPENMM_EXPORT void openmm_discrete1dfunction_copy_(const OpenMM_Discrete1DFunction*& target, OpenMM_Discrete1DFunction*& result) {
    result = OpenMM_Discrete1DFunction_Copy(target);
}
OPENMM_EXPORT void OPENMM_DISCRETE1DFUNCTION_COPY(const OpenMM_Discrete1DFunction*& target, OpenMM_Discrete1DFunction*& result) {
    result = OpenMM_Discrete1DFunction_Copy(target);
}

/* OpenMM::CustomManyParticleForce */
OPENMM_EXPORT void openmm_custommanyparticleforce_create_(OpenMM_CustomManyParticleForce*& result, int const& particlesPerSet, const char* energy, int energy_length) {
    result = OpenMM_CustomManyParticleForce_create(particlesPerSet, makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMMANYPARTICLEFORCE_CREATE(OpenMM_CustomManyParticleForce*& result, int const& particlesPerSet, const char* energy, int energy_length) {
    result = OpenMM_CustomManyParticleForce_create(particlesPerSet, makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void openmm_custommanyparticleforce_destroy_(OpenMM_CustomManyParticleForce*& destroy) {
    OpenMM_CustomManyParticleForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_CUSTOMMANYPARTICLEFORCE_DESTROY(OpenMM_CustomManyParticleForce*& destroy) {
    OpenMM_CustomManyParticleForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT int openmm_custommanyparticleforce_getnumparticlesperset_(const OpenMM_CustomManyParticleForce*& target) {
    return OpenMM_CustomManyParticleForce_getNumParticlesPerSet(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMMANYPARTICLEFORCE_GETNUMPARTICLESPERSET(const OpenMM_CustomManyParticleForce*& target) {
    return OpenMM_CustomManyParticleForce_getNumParticlesPerSet(target);
}
OPENMM_EXPORT int openmm_custommanyparticleforce_getnumparticles_(const OpenMM_CustomManyParticleForce*& target) {
    return OpenMM_CustomManyParticleForce_getNumParticles(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMMANYPARTICLEFORCE_GETNUMPARTICLES(const OpenMM_CustomManyParticleForce*& target) {
    return OpenMM_CustomManyParticleForce_getNumParticles(target);
}
OPENMM_EXPORT int openmm_custommanyparticleforce_getnumexclusions_(const OpenMM_CustomManyParticleForce*& target) {
    return OpenMM_CustomManyParticleForce_getNumExclusions(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMMANYPARTICLEFORCE_GETNUMEXCLUSIONS(const OpenMM_CustomManyParticleForce*& target) {
    return OpenMM_CustomManyParticleForce_getNumExclusions(target);
}
OPENMM_EXPORT int openmm_custommanyparticleforce_getnumperparticleparameters_(const OpenMM_CustomManyParticleForce*& target) {
    return OpenMM_CustomManyParticleForce_getNumPerParticleParameters(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMMANYPARTICLEFORCE_GETNUMPERPARTICLEPARAMETERS(const OpenMM_CustomManyParticleForce*& target) {
    return OpenMM_CustomManyParticleForce_getNumPerParticleParameters(target);
}
OPENMM_EXPORT int openmm_custommanyparticleforce_getnumglobalparameters_(const OpenMM_CustomManyParticleForce*& target) {
    return OpenMM_CustomManyParticleForce_getNumGlobalParameters(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMMANYPARTICLEFORCE_GETNUMGLOBALPARAMETERS(const OpenMM_CustomManyParticleForce*& target) {
    return OpenMM_CustomManyParticleForce_getNumGlobalParameters(target);
}
OPENMM_EXPORT int openmm_custommanyparticleforce_getnumtabulatedfunctions_(const OpenMM_CustomManyParticleForce*& target) {
    return OpenMM_CustomManyParticleForce_getNumTabulatedFunctions(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMMANYPARTICLEFORCE_GETNUMTABULATEDFUNCTIONS(const OpenMM_CustomManyParticleForce*& target) {
    return OpenMM_CustomManyParticleForce_getNumTabulatedFunctions(target);
}
OPENMM_EXPORT void openmm_custommanyparticleforce_getenergyfunction_(const OpenMM_CustomManyParticleForce*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomManyParticleForce_getEnergyFunction(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMMANYPARTICLEFORCE_GETENERGYFUNCTION(const OpenMM_CustomManyParticleForce*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomManyParticleForce_getEnergyFunction(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_custommanyparticleforce_setenergyfunction_(OpenMM_CustomManyParticleForce*& target, const char* energy, int energy_length) {
    OpenMM_CustomManyParticleForce_setEnergyFunction(target, makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMMANYPARTICLEFORCE_SETENERGYFUNCTION(OpenMM_CustomManyParticleForce*& target, const char* energy, int energy_length) {
    OpenMM_CustomManyParticleForce_setEnergyFunction(target, makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void openmm_custommanyparticleforce_getnonbondedmethod_(const OpenMM_CustomManyParticleForce*& target, OpenMM_CustomManyParticleForce_NonbondedMethod& result) {
    result = OpenMM_CustomManyParticleForce_getNonbondedMethod(target);
}
OPENMM_EXPORT void OPENMM_CUSTOMMANYPARTICLEFORCE_GETNONBONDEDMETHOD(const OpenMM_CustomManyParticleForce*& target, OpenMM_CustomManyParticleForce_NonbondedMethod& result) {
    result = OpenMM_CustomManyParticleForce_getNonbondedMethod(target);
}
OPENMM_EXPORT void openmm_custommanyparticleforce_setnonbondedmethod_(OpenMM_CustomManyParticleForce*& target, OpenMM_CustomManyParticleForce_NonbondedMethod& method) {
    OpenMM_CustomManyParticleForce_setNonbondedMethod(target, method);
}
OPENMM_EXPORT void OPENMM_CUSTOMMANYPARTICLEFORCE_SETNONBONDEDMETHOD(OpenMM_CustomManyParticleForce*& target, OpenMM_CustomManyParticleForce_NonbondedMethod& method) {
    OpenMM_CustomManyParticleForce_setNonbondedMethod(target, method);
}
OPENMM_EXPORT void openmm_custommanyparticleforce_getpermutationmode_(const OpenMM_CustomManyParticleForce*& target, OpenMM_CustomManyParticleForce_PermutationMode& result) {
    result = OpenMM_CustomManyParticleForce_getPermutationMode(target);
}
OPENMM_EXPORT void OPENMM_CUSTOMMANYPARTICLEFORCE_GETPERMUTATIONMODE(const OpenMM_CustomManyParticleForce*& target, OpenMM_CustomManyParticleForce_PermutationMode& result) {
    result = OpenMM_CustomManyParticleForce_getPermutationMode(target);
}
OPENMM_EXPORT void openmm_custommanyparticleforce_setpermutationmode_(OpenMM_CustomManyParticleForce*& target, OpenMM_CustomManyParticleForce_PermutationMode& mode) {
    OpenMM_CustomManyParticleForce_setPermutationMode(target, mode);
}
OPENMM_EXPORT void OPENMM_CUSTOMMANYPARTICLEFORCE_SETPERMUTATIONMODE(OpenMM_CustomManyParticleForce*& target, OpenMM_CustomManyParticleForce_PermutationMode& mode) {
    OpenMM_CustomManyParticleForce_setPermutationMode(target, mode);
}
OPENMM_EXPORT double openmm_custommanyparticleforce_getcutoffdistance_(const OpenMM_CustomManyParticleForce*& target) {
    return OpenMM_CustomManyParticleForce_getCutoffDistance(target);
}
OPENMM_EXPORT double OPENMM_CUSTOMMANYPARTICLEFORCE_GETCUTOFFDISTANCE(const OpenMM_CustomManyParticleForce*& target) {
    return OpenMM_CustomManyParticleForce_getCutoffDistance(target);
}
OPENMM_EXPORT void openmm_custommanyparticleforce_setcutoffdistance_(OpenMM_CustomManyParticleForce*& target, double const& distance) {
    OpenMM_CustomManyParticleForce_setCutoffDistance(target, distance);
}
OPENMM_EXPORT void OPENMM_CUSTOMMANYPARTICLEFORCE_SETCUTOFFDISTANCE(OpenMM_CustomManyParticleForce*& target, double const& distance) {
    OpenMM_CustomManyParticleForce_setCutoffDistance(target, distance);
}
OPENMM_EXPORT int openmm_custommanyparticleforce_addperparticleparameter_(OpenMM_CustomManyParticleForce*& target, const char* name, int name_length) {
    return OpenMM_CustomManyParticleForce_addPerParticleParameter(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT int OPENMM_CUSTOMMANYPARTICLEFORCE_ADDPERPARTICLEPARAMETER(OpenMM_CustomManyParticleForce*& target, const char* name, int name_length) {
    return OpenMM_CustomManyParticleForce_addPerParticleParameter(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void openmm_custommanyparticleforce_getperparticleparametername_(const OpenMM_CustomManyParticleForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomManyParticleForce_getPerParticleParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMMANYPARTICLEFORCE_GETPERPARTICLEPARAMETERNAME(const OpenMM_CustomManyParticleForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomManyParticleForce_getPerParticleParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_custommanyparticleforce_setperparticleparametername_(OpenMM_CustomManyParticleForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomManyParticleForce_setPerParticleParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMMANYPARTICLEFORCE_SETPERPARTICLEPARAMETERNAME(OpenMM_CustomManyParticleForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomManyParticleForce_setPerParticleParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT int openmm_custommanyparticleforce_addglobalparameter_(OpenMM_CustomManyParticleForce*& target, const char* name, double const& defaultValue, int name_length) {
    return OpenMM_CustomManyParticleForce_addGlobalParameter(target, makeString(name, name_length).c_str(), defaultValue);
}
OPENMM_EXPORT int OPENMM_CUSTOMMANYPARTICLEFORCE_ADDGLOBALPARAMETER(OpenMM_CustomManyParticleForce*& target, const char* name, double const& defaultValue, int name_length) {
    return OpenMM_CustomManyParticleForce_addGlobalParameter(target, makeString(name, name_length).c_str(), defaultValue);
}
OPENMM_EXPORT void openmm_custommanyparticleforce_getglobalparametername_(const OpenMM_CustomManyParticleForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomManyParticleForce_getGlobalParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMMANYPARTICLEFORCE_GETGLOBALPARAMETERNAME(const OpenMM_CustomManyParticleForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomManyParticleForce_getGlobalParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_custommanyparticleforce_setglobalparametername_(OpenMM_CustomManyParticleForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomManyParticleForce_setGlobalParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMMANYPARTICLEFORCE_SETGLOBALPARAMETERNAME(OpenMM_CustomManyParticleForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomManyParticleForce_setGlobalParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT double openmm_custommanyparticleforce_getglobalparameterdefaultvalue_(const OpenMM_CustomManyParticleForce*& target, int const& index) {
    return OpenMM_CustomManyParticleForce_getGlobalParameterDefaultValue(target, index);
}
OPENMM_EXPORT double OPENMM_CUSTOMMANYPARTICLEFORCE_GETGLOBALPARAMETERDEFAULTVALUE(const OpenMM_CustomManyParticleForce*& target, int const& index) {
    return OpenMM_CustomManyParticleForce_getGlobalParameterDefaultValue(target, index);
}
OPENMM_EXPORT void openmm_custommanyparticleforce_setglobalparameterdefaultvalue_(OpenMM_CustomManyParticleForce*& target, int const& index, double const& defaultValue) {
    OpenMM_CustomManyParticleForce_setGlobalParameterDefaultValue(target, index, defaultValue);
}
OPENMM_EXPORT void OPENMM_CUSTOMMANYPARTICLEFORCE_SETGLOBALPARAMETERDEFAULTVALUE(OpenMM_CustomManyParticleForce*& target, int const& index, double const& defaultValue) {
    OpenMM_CustomManyParticleForce_setGlobalParameterDefaultValue(target, index, defaultValue);
}
OPENMM_EXPORT int openmm_custommanyparticleforce_addparticle_(OpenMM_CustomManyParticleForce*& target, const OpenMM_DoubleArray*& parameters, int const& type) {
    return OpenMM_CustomManyParticleForce_addParticle(target, parameters, type);
}
OPENMM_EXPORT int OPENMM_CUSTOMMANYPARTICLEFORCE_ADDPARTICLE(OpenMM_CustomManyParticleForce*& target, const OpenMM_DoubleArray*& parameters, int const& type) {
    return OpenMM_CustomManyParticleForce_addParticle(target, parameters, type);
}
OPENMM_EXPORT void openmm_custommanyparticleforce_getparticleparameters_(const OpenMM_CustomManyParticleForce*& target, int const& index, OpenMM_DoubleArray*& parameters, int* type) {
    OpenMM_CustomManyParticleForce_getParticleParameters(target, index, parameters, type);
}
OPENMM_EXPORT void OPENMM_CUSTOMMANYPARTICLEFORCE_GETPARTICLEPARAMETERS(const OpenMM_CustomManyParticleForce*& target, int const& index, OpenMM_DoubleArray*& parameters, int* type) {
    OpenMM_CustomManyParticleForce_getParticleParameters(target, index, parameters, type);
}
OPENMM_EXPORT void openmm_custommanyparticleforce_setparticleparameters_(OpenMM_CustomManyParticleForce*& target, int const& index, const OpenMM_DoubleArray*& parameters, int const& type) {
    OpenMM_CustomManyParticleForce_setParticleParameters(target, index, parameters, type);
}
OPENMM_EXPORT void OPENMM_CUSTOMMANYPARTICLEFORCE_SETPARTICLEPARAMETERS(OpenMM_CustomManyParticleForce*& target, int const& index, const OpenMM_DoubleArray*& parameters, int const& type) {
    OpenMM_CustomManyParticleForce_setParticleParameters(target, index, parameters, type);
}
OPENMM_EXPORT int openmm_custommanyparticleforce_addexclusion_(OpenMM_CustomManyParticleForce*& target, int const& particle1, int const& particle2) {
    return OpenMM_CustomManyParticleForce_addExclusion(target, particle1, particle2);
}
OPENMM_EXPORT int OPENMM_CUSTOMMANYPARTICLEFORCE_ADDEXCLUSION(OpenMM_CustomManyParticleForce*& target, int const& particle1, int const& particle2) {
    return OpenMM_CustomManyParticleForce_addExclusion(target, particle1, particle2);
}
OPENMM_EXPORT void openmm_custommanyparticleforce_getexclusionparticles_(const OpenMM_CustomManyParticleForce*& target, int const& index, int* particle1, int* particle2) {
    OpenMM_CustomManyParticleForce_getExclusionParticles(target, index, particle1, particle2);
}
OPENMM_EXPORT void OPENMM_CUSTOMMANYPARTICLEFORCE_GETEXCLUSIONPARTICLES(const OpenMM_CustomManyParticleForce*& target, int const& index, int* particle1, int* particle2) {
    OpenMM_CustomManyParticleForce_getExclusionParticles(target, index, particle1, particle2);
}
OPENMM_EXPORT void openmm_custommanyparticleforce_setexclusionparticles_(OpenMM_CustomManyParticleForce*& target, int const& index, int const& particle1, int const& particle2) {
    OpenMM_CustomManyParticleForce_setExclusionParticles(target, index, particle1, particle2);
}
OPENMM_EXPORT void OPENMM_CUSTOMMANYPARTICLEFORCE_SETEXCLUSIONPARTICLES(OpenMM_CustomManyParticleForce*& target, int const& index, int const& particle1, int const& particle2) {
    OpenMM_CustomManyParticleForce_setExclusionParticles(target, index, particle1, particle2);
}
OPENMM_EXPORT void openmm_custommanyparticleforce_createexclusionsfrombonds_(OpenMM_CustomManyParticleForce*& target, const OpenMM_BondArray*& bonds, int const& bondCutoff) {
    OpenMM_CustomManyParticleForce_createExclusionsFromBonds(target, bonds, bondCutoff);
}
OPENMM_EXPORT void OPENMM_CUSTOMMANYPARTICLEFORCE_CREATEEXCLUSIONSFROMBONDS(OpenMM_CustomManyParticleForce*& target, const OpenMM_BondArray*& bonds, int const& bondCutoff) {
    OpenMM_CustomManyParticleForce_createExclusionsFromBonds(target, bonds, bondCutoff);
}
OPENMM_EXPORT void openmm_custommanyparticleforce_gettypefilter_(const OpenMM_CustomManyParticleForce*& target, int const& index, OpenMM_IntSet*& types) {
    OpenMM_CustomManyParticleForce_getTypeFilter(target, index, types);
}
OPENMM_EXPORT void OPENMM_CUSTOMMANYPARTICLEFORCE_GETTYPEFILTER(const OpenMM_CustomManyParticleForce*& target, int const& index, OpenMM_IntSet*& types) {
    OpenMM_CustomManyParticleForce_getTypeFilter(target, index, types);
}
OPENMM_EXPORT void openmm_custommanyparticleforce_settypefilter_(OpenMM_CustomManyParticleForce*& target, int const& index, const OpenMM_IntSet*& types) {
    OpenMM_CustomManyParticleForce_setTypeFilter(target, index, types);
}
OPENMM_EXPORT void OPENMM_CUSTOMMANYPARTICLEFORCE_SETTYPEFILTER(OpenMM_CustomManyParticleForce*& target, int const& index, const OpenMM_IntSet*& types) {
    OpenMM_CustomManyParticleForce_setTypeFilter(target, index, types);
}
OPENMM_EXPORT int openmm_custommanyparticleforce_addtabulatedfunction_(OpenMM_CustomManyParticleForce*& target, const char* name, OpenMM_TabulatedFunction*& function, int name_length) {
    return OpenMM_CustomManyParticleForce_addTabulatedFunction(target, makeString(name, name_length).c_str(), function);
}
OPENMM_EXPORT int OPENMM_CUSTOMMANYPARTICLEFORCE_ADDTABULATEDFUNCTION(OpenMM_CustomManyParticleForce*& target, const char* name, OpenMM_TabulatedFunction*& function, int name_length) {
    return OpenMM_CustomManyParticleForce_addTabulatedFunction(target, makeString(name, name_length).c_str(), function);
}
OPENMM_EXPORT void openmm_custommanyparticleforce_gettabulatedfunction_(OpenMM_CustomManyParticleForce*& target, int const& index, OpenMM_TabulatedFunction*& result) {
    result = OpenMM_CustomManyParticleForce_getTabulatedFunction(target, index);
}
OPENMM_EXPORT void OPENMM_CUSTOMMANYPARTICLEFORCE_GETTABULATEDFUNCTION(OpenMM_CustomManyParticleForce*& target, int const& index, OpenMM_TabulatedFunction*& result) {
    result = OpenMM_CustomManyParticleForce_getTabulatedFunction(target, index);
}
OPENMM_EXPORT void openmm_custommanyparticleforce_gettabulatedfunctionname_(const OpenMM_CustomManyParticleForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomManyParticleForce_getTabulatedFunctionName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMMANYPARTICLEFORCE_GETTABULATEDFUNCTIONNAME(const OpenMM_CustomManyParticleForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomManyParticleForce_getTabulatedFunctionName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_custommanyparticleforce_updateparametersincontext_(OpenMM_CustomManyParticleForce*& target, OpenMM_Context*& context) {
    OpenMM_CustomManyParticleForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void OPENMM_CUSTOMMANYPARTICLEFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_CustomManyParticleForce*& target, OpenMM_Context*& context) {
    OpenMM_CustomManyParticleForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void openmm_custommanyparticleforce_usesperiodicboundaryconditions_(const OpenMM_CustomManyParticleForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CustomManyParticleForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_CUSTOMMANYPARTICLEFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_CustomManyParticleForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CustomManyParticleForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::HarmonicAngleForce */
OPENMM_EXPORT void openmm_harmonicangleforce_create_(OpenMM_HarmonicAngleForce*& result) {
    result = OpenMM_HarmonicAngleForce_create();
}
OPENMM_EXPORT void OPENMM_HARMONICANGLEFORCE_CREATE(OpenMM_HarmonicAngleForce*& result) {
    result = OpenMM_HarmonicAngleForce_create();
}
OPENMM_EXPORT void openmm_harmonicangleforce_destroy_(OpenMM_HarmonicAngleForce*& destroy) {
    OpenMM_HarmonicAngleForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_HARMONICANGLEFORCE_DESTROY(OpenMM_HarmonicAngleForce*& destroy) {
    OpenMM_HarmonicAngleForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT int openmm_harmonicangleforce_getnumangles_(const OpenMM_HarmonicAngleForce*& target) {
    return OpenMM_HarmonicAngleForce_getNumAngles(target);
}
OPENMM_EXPORT int OPENMM_HARMONICANGLEFORCE_GETNUMANGLES(const OpenMM_HarmonicAngleForce*& target) {
    return OpenMM_HarmonicAngleForce_getNumAngles(target);
}
OPENMM_EXPORT int openmm_harmonicangleforce_addangle_(OpenMM_HarmonicAngleForce*& target, int const& particle1, int const& particle2, int const& particle3, double const& angle, double const& k) {
    return OpenMM_HarmonicAngleForce_addAngle(target, particle1, particle2, particle3, angle, k);
}
OPENMM_EXPORT int OPENMM_HARMONICANGLEFORCE_ADDANGLE(OpenMM_HarmonicAngleForce*& target, int const& particle1, int const& particle2, int const& particle3, double const& angle, double const& k) {
    return OpenMM_HarmonicAngleForce_addAngle(target, particle1, particle2, particle3, angle, k);
}
OPENMM_EXPORT void openmm_harmonicangleforce_getangleparameters_(const OpenMM_HarmonicAngleForce*& target, int const& index, int* particle1, int* particle2, int* particle3, double* angle, double* k) {
    OpenMM_HarmonicAngleForce_getAngleParameters(target, index, particle1, particle2, particle3, angle, k);
}
OPENMM_EXPORT void OPENMM_HARMONICANGLEFORCE_GETANGLEPARAMETERS(const OpenMM_HarmonicAngleForce*& target, int const& index, int* particle1, int* particle2, int* particle3, double* angle, double* k) {
    OpenMM_HarmonicAngleForce_getAngleParameters(target, index, particle1, particle2, particle3, angle, k);
}
OPENMM_EXPORT void openmm_harmonicangleforce_setangleparameters_(OpenMM_HarmonicAngleForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, double const& angle, double const& k) {
    OpenMM_HarmonicAngleForce_setAngleParameters(target, index, particle1, particle2, particle3, angle, k);
}
OPENMM_EXPORT void OPENMM_HARMONICANGLEFORCE_SETANGLEPARAMETERS(OpenMM_HarmonicAngleForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, double const& angle, double const& k) {
    OpenMM_HarmonicAngleForce_setAngleParameters(target, index, particle1, particle2, particle3, angle, k);
}
OPENMM_EXPORT void openmm_harmonicangleforce_updateparametersincontext_(OpenMM_HarmonicAngleForce*& target, OpenMM_Context*& context) {
    OpenMM_HarmonicAngleForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void OPENMM_HARMONICANGLEFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_HarmonicAngleForce*& target, OpenMM_Context*& context) {
    OpenMM_HarmonicAngleForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void openmm_harmonicangleforce_setusesperiodicboundaryconditions_(OpenMM_HarmonicAngleForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_HarmonicAngleForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT void OPENMM_HARMONICANGLEFORCE_SETUSESPERIODICBOUNDARYCONDITIONS(OpenMM_HarmonicAngleForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_HarmonicAngleForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT void openmm_harmonicangleforce_usesperiodicboundaryconditions_(const OpenMM_HarmonicAngleForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_HarmonicAngleForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_HARMONICANGLEFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_HarmonicAngleForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_HarmonicAngleForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::CompoundIntegrator */
OPENMM_EXPORT void openmm_compoundintegrator_create_(OpenMM_CompoundIntegrator*& result) {
    result = OpenMM_CompoundIntegrator_create();
}
OPENMM_EXPORT void OPENMM_COMPOUNDINTEGRATOR_CREATE(OpenMM_CompoundIntegrator*& result) {
    result = OpenMM_CompoundIntegrator_create();
}
OPENMM_EXPORT void openmm_compoundintegrator_destroy_(OpenMM_CompoundIntegrator*& destroy) {
    OpenMM_CompoundIntegrator_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_COMPOUNDINTEGRATOR_DESTROY(OpenMM_CompoundIntegrator*& destroy) {
    OpenMM_CompoundIntegrator_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT int openmm_compoundintegrator_getnumintegrators_(const OpenMM_CompoundIntegrator*& target) {
    return OpenMM_CompoundIntegrator_getNumIntegrators(target);
}
OPENMM_EXPORT int OPENMM_COMPOUNDINTEGRATOR_GETNUMINTEGRATORS(const OpenMM_CompoundIntegrator*& target) {
    return OpenMM_CompoundIntegrator_getNumIntegrators(target);
}
OPENMM_EXPORT int openmm_compoundintegrator_addintegrator_(OpenMM_CompoundIntegrator*& target, OpenMM_Integrator*& integrator) {
    return OpenMM_CompoundIntegrator_addIntegrator(target, integrator);
}
OPENMM_EXPORT int OPENMM_COMPOUNDINTEGRATOR_ADDINTEGRATOR(OpenMM_CompoundIntegrator*& target, OpenMM_Integrator*& integrator) {
    return OpenMM_CompoundIntegrator_addIntegrator(target, integrator);
}
OPENMM_EXPORT void openmm_compoundintegrator_getintegrator_(OpenMM_CompoundIntegrator*& target, int const& index, OpenMM_Integrator*& result) {
    result = OpenMM_CompoundIntegrator_getIntegrator(target, index);
}
OPENMM_EXPORT void OPENMM_COMPOUNDINTEGRATOR_GETINTEGRATOR(OpenMM_CompoundIntegrator*& target, int const& index, OpenMM_Integrator*& result) {
    result = OpenMM_CompoundIntegrator_getIntegrator(target, index);
}
OPENMM_EXPORT int openmm_compoundintegrator_getcurrentintegrator_(const OpenMM_CompoundIntegrator*& target) {
    return OpenMM_CompoundIntegrator_getCurrentIntegrator(target);
}
OPENMM_EXPORT int OPENMM_COMPOUNDINTEGRATOR_GETCURRENTINTEGRATOR(const OpenMM_CompoundIntegrator*& target) {
    return OpenMM_CompoundIntegrator_getCurrentIntegrator(target);
}
OPENMM_EXPORT void openmm_compoundintegrator_setcurrentintegrator_(OpenMM_CompoundIntegrator*& target, int const& index) {
    OpenMM_CompoundIntegrator_setCurrentIntegrator(target, index);
}
OPENMM_EXPORT void OPENMM_COMPOUNDINTEGRATOR_SETCURRENTINTEGRATOR(OpenMM_CompoundIntegrator*& target, int const& index) {
    OpenMM_CompoundIntegrator_setCurrentIntegrator(target, index);
}
OPENMM_EXPORT double openmm_compoundintegrator_getstepsize_(const OpenMM_CompoundIntegrator*& target) {
    return OpenMM_CompoundIntegrator_getStepSize(target);
}
OPENMM_EXPORT double OPENMM_COMPOUNDINTEGRATOR_GETSTEPSIZE(const OpenMM_CompoundIntegrator*& target) {
    return OpenMM_CompoundIntegrator_getStepSize(target);
}
OPENMM_EXPORT void openmm_compoundintegrator_setstepsize_(OpenMM_CompoundIntegrator*& target, double const& size) {
    OpenMM_CompoundIntegrator_setStepSize(target, size);
}
OPENMM_EXPORT void OPENMM_COMPOUNDINTEGRATOR_SETSTEPSIZE(OpenMM_CompoundIntegrator*& target, double const& size) {
    OpenMM_CompoundIntegrator_setStepSize(target, size);
}
OPENMM_EXPORT double openmm_compoundintegrator_getconstrainttolerance_(const OpenMM_CompoundIntegrator*& target) {
    return OpenMM_CompoundIntegrator_getConstraintTolerance(target);
}
OPENMM_EXPORT double OPENMM_COMPOUNDINTEGRATOR_GETCONSTRAINTTOLERANCE(const OpenMM_CompoundIntegrator*& target) {
    return OpenMM_CompoundIntegrator_getConstraintTolerance(target);
}
OPENMM_EXPORT void openmm_compoundintegrator_setconstrainttolerance_(OpenMM_CompoundIntegrator*& target, double const& tol) {
    OpenMM_CompoundIntegrator_setConstraintTolerance(target, tol);
}
OPENMM_EXPORT void OPENMM_COMPOUNDINTEGRATOR_SETCONSTRAINTTOLERANCE(OpenMM_CompoundIntegrator*& target, double const& tol) {
    OpenMM_CompoundIntegrator_setConstraintTolerance(target, tol);
}
OPENMM_EXPORT void openmm_compoundintegrator_step_(OpenMM_CompoundIntegrator*& target, int const& steps) {
    OpenMM_CompoundIntegrator_step(target, steps);
}
OPENMM_EXPORT void OPENMM_COMPOUNDINTEGRATOR_STEP(OpenMM_CompoundIntegrator*& target, int const& steps) {
    OpenMM_CompoundIntegrator_step(target, steps);
}

/* OpenMM::MonteCarloAnisotropicBarostat */
OPENMM_EXPORT void openmm_montecarloanisotropicbarostat_create_(OpenMM_MonteCarloAnisotropicBarostat*& result, const OpenMM_Vec3* defaultPressure, double const& defaultTemperature, OpenMM_Boolean& scaleX, OpenMM_Boolean& scaleY, OpenMM_Boolean& scaleZ, int const& frequency) {
    result = OpenMM_MonteCarloAnisotropicBarostat_create(defaultPressure, defaultTemperature, scaleX, scaleY, scaleZ, frequency);
}
OPENMM_EXPORT void OPENMM_MONTECARLOANISOTROPICBAROSTAT_CREATE(OpenMM_MonteCarloAnisotropicBarostat*& result, const OpenMM_Vec3* defaultPressure, double const& defaultTemperature, OpenMM_Boolean& scaleX, OpenMM_Boolean& scaleY, OpenMM_Boolean& scaleZ, int const& frequency) {
    result = OpenMM_MonteCarloAnisotropicBarostat_create(defaultPressure, defaultTemperature, scaleX, scaleY, scaleZ, frequency);
}
OPENMM_EXPORT void openmm_montecarloanisotropicbarostat_destroy_(OpenMM_MonteCarloAnisotropicBarostat*& destroy) {
    OpenMM_MonteCarloAnisotropicBarostat_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_MONTECARLOANISOTROPICBAROSTAT_DESTROY(OpenMM_MonteCarloAnisotropicBarostat*& destroy) {
    OpenMM_MonteCarloAnisotropicBarostat_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void openmm_montecarloanisotropicbarostat_pressurex_(char* result, int result_length) {
    const char* result_chars = OpenMM_MonteCarloAnisotropicBarostat_PressureX();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_MONTECARLOANISOTROPICBAROSTAT_PRESSUREX(char* result, int result_length) {
    const char* result_chars = OpenMM_MonteCarloAnisotropicBarostat_PressureX();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_montecarloanisotropicbarostat_pressurey_(char* result, int result_length) {
    const char* result_chars = OpenMM_MonteCarloAnisotropicBarostat_PressureY();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_MONTECARLOANISOTROPICBAROSTAT_PRESSUREY(char* result, int result_length) {
    const char* result_chars = OpenMM_MonteCarloAnisotropicBarostat_PressureY();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_montecarloanisotropicbarostat_pressurez_(char* result, int result_length) {
    const char* result_chars = OpenMM_MonteCarloAnisotropicBarostat_PressureZ();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_MONTECARLOANISOTROPICBAROSTAT_PRESSUREZ(char* result, int result_length) {
    const char* result_chars = OpenMM_MonteCarloAnisotropicBarostat_PressureZ();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_montecarloanisotropicbarostat_temperature_(char* result, int result_length) {
    const char* result_chars = OpenMM_MonteCarloAnisotropicBarostat_Temperature();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_MONTECARLOANISOTROPICBAROSTAT_TEMPERATURE(char* result, int result_length) {
    const char* result_chars = OpenMM_MonteCarloAnisotropicBarostat_Temperature();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_montecarloanisotropicbarostat_getdefaultpressure_(const OpenMM_MonteCarloAnisotropicBarostat*& target, const OpenMM_Vec3*& result) {
    result = OpenMM_MonteCarloAnisotropicBarostat_getDefaultPressure(target);
}
OPENMM_EXPORT void OPENMM_MONTECARLOANISOTROPICBAROSTAT_GETDEFAULTPRESSURE(const OpenMM_MonteCarloAnisotropicBarostat*& target, const OpenMM_Vec3*& result) {
    result = OpenMM_MonteCarloAnisotropicBarostat_getDefaultPressure(target);
}
OPENMM_EXPORT void openmm_montecarloanisotropicbarostat_setdefaultpressure_(OpenMM_MonteCarloAnisotropicBarostat*& target, const OpenMM_Vec3* pressure) {
    OpenMM_MonteCarloAnisotropicBarostat_setDefaultPressure(target, pressure);
}
OPENMM_EXPORT void OPENMM_MONTECARLOANISOTROPICBAROSTAT_SETDEFAULTPRESSURE(OpenMM_MonteCarloAnisotropicBarostat*& target, const OpenMM_Vec3* pressure) {
    OpenMM_MonteCarloAnisotropicBarostat_setDefaultPressure(target, pressure);
}
OPENMM_EXPORT void openmm_montecarloanisotropicbarostat_getscalex_(const OpenMM_MonteCarloAnisotropicBarostat*& target, OpenMM_Boolean& result) {
    result = OpenMM_MonteCarloAnisotropicBarostat_getScaleX(target);
}
OPENMM_EXPORT void OPENMM_MONTECARLOANISOTROPICBAROSTAT_GETSCALEX(const OpenMM_MonteCarloAnisotropicBarostat*& target, OpenMM_Boolean& result) {
    result = OpenMM_MonteCarloAnisotropicBarostat_getScaleX(target);
}
OPENMM_EXPORT void openmm_montecarloanisotropicbarostat_getscaley_(const OpenMM_MonteCarloAnisotropicBarostat*& target, OpenMM_Boolean& result) {
    result = OpenMM_MonteCarloAnisotropicBarostat_getScaleY(target);
}
OPENMM_EXPORT void OPENMM_MONTECARLOANISOTROPICBAROSTAT_GETSCALEY(const OpenMM_MonteCarloAnisotropicBarostat*& target, OpenMM_Boolean& result) {
    result = OpenMM_MonteCarloAnisotropicBarostat_getScaleY(target);
}
OPENMM_EXPORT void openmm_montecarloanisotropicbarostat_getscalez_(const OpenMM_MonteCarloAnisotropicBarostat*& target, OpenMM_Boolean& result) {
    result = OpenMM_MonteCarloAnisotropicBarostat_getScaleZ(target);
}
OPENMM_EXPORT void OPENMM_MONTECARLOANISOTROPICBAROSTAT_GETSCALEZ(const OpenMM_MonteCarloAnisotropicBarostat*& target, OpenMM_Boolean& result) {
    result = OpenMM_MonteCarloAnisotropicBarostat_getScaleZ(target);
}
OPENMM_EXPORT int openmm_montecarloanisotropicbarostat_getfrequency_(const OpenMM_MonteCarloAnisotropicBarostat*& target) {
    return OpenMM_MonteCarloAnisotropicBarostat_getFrequency(target);
}
OPENMM_EXPORT int OPENMM_MONTECARLOANISOTROPICBAROSTAT_GETFREQUENCY(const OpenMM_MonteCarloAnisotropicBarostat*& target) {
    return OpenMM_MonteCarloAnisotropicBarostat_getFrequency(target);
}
OPENMM_EXPORT void openmm_montecarloanisotropicbarostat_setfrequency_(OpenMM_MonteCarloAnisotropicBarostat*& target, int const& freq) {
    OpenMM_MonteCarloAnisotropicBarostat_setFrequency(target, freq);
}
OPENMM_EXPORT void OPENMM_MONTECARLOANISOTROPICBAROSTAT_SETFREQUENCY(OpenMM_MonteCarloAnisotropicBarostat*& target, int const& freq) {
    OpenMM_MonteCarloAnisotropicBarostat_setFrequency(target, freq);
}
OPENMM_EXPORT double openmm_montecarloanisotropicbarostat_getdefaulttemperature_(const OpenMM_MonteCarloAnisotropicBarostat*& target) {
    return OpenMM_MonteCarloAnisotropicBarostat_getDefaultTemperature(target);
}
OPENMM_EXPORT double OPENMM_MONTECARLOANISOTROPICBAROSTAT_GETDEFAULTTEMPERATURE(const OpenMM_MonteCarloAnisotropicBarostat*& target) {
    return OpenMM_MonteCarloAnisotropicBarostat_getDefaultTemperature(target);
}
OPENMM_EXPORT void openmm_montecarloanisotropicbarostat_setdefaulttemperature_(OpenMM_MonteCarloAnisotropicBarostat*& target, double const& temp) {
    OpenMM_MonteCarloAnisotropicBarostat_setDefaultTemperature(target, temp);
}
OPENMM_EXPORT void OPENMM_MONTECARLOANISOTROPICBAROSTAT_SETDEFAULTTEMPERATURE(OpenMM_MonteCarloAnisotropicBarostat*& target, double const& temp) {
    OpenMM_MonteCarloAnisotropicBarostat_setDefaultTemperature(target, temp);
}
OPENMM_EXPORT int openmm_montecarloanisotropicbarostat_getrandomnumberseed_(const OpenMM_MonteCarloAnisotropicBarostat*& target) {
    return OpenMM_MonteCarloAnisotropicBarostat_getRandomNumberSeed(target);
}
OPENMM_EXPORT int OPENMM_MONTECARLOANISOTROPICBAROSTAT_GETRANDOMNUMBERSEED(const OpenMM_MonteCarloAnisotropicBarostat*& target) {
    return OpenMM_MonteCarloAnisotropicBarostat_getRandomNumberSeed(target);
}
OPENMM_EXPORT void openmm_montecarloanisotropicbarostat_setrandomnumberseed_(OpenMM_MonteCarloAnisotropicBarostat*& target, int const& seed) {
    OpenMM_MonteCarloAnisotropicBarostat_setRandomNumberSeed(target, seed);
}
OPENMM_EXPORT void OPENMM_MONTECARLOANISOTROPICBAROSTAT_SETRANDOMNUMBERSEED(OpenMM_MonteCarloAnisotropicBarostat*& target, int const& seed) {
    OpenMM_MonteCarloAnisotropicBarostat_setRandomNumberSeed(target, seed);
}
OPENMM_EXPORT void openmm_montecarloanisotropicbarostat_usesperiodicboundarycondit_(const OpenMM_MonteCarloAnisotropicBarostat*& target, OpenMM_Boolean& result) {
    result = OpenMM_MonteCarloAnisotropicBarostat_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_MONTECARLOANISOTROPICBAROSTAT_USESPERIODICBOUNDARYCONDIT(const OpenMM_MonteCarloAnisotropicBarostat*& target, OpenMM_Boolean& result) {
    result = OpenMM_MonteCarloAnisotropicBarostat_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::NoseHooverIntegrator */
OPENMM_EXPORT void openmm_nosehooverintegrator_create_(OpenMM_NoseHooverIntegrator*& result, double const& stepSize) {
    result = OpenMM_NoseHooverIntegrator_create(stepSize);
}
OPENMM_EXPORT void OPENMM_NOSEHOOVERINTEGRATOR_CREATE(OpenMM_NoseHooverIntegrator*& result, double const& stepSize) {
    result = OpenMM_NoseHooverIntegrator_create(stepSize);
}
OPENMM_EXPORT void openmm_nosehooverintegrator_create_2_(OpenMM_NoseHooverIntegrator*& result, double const& temperature, double const& collisionFrequency, double const& stepSize, int const& chainLength, int const& numMTS, int const& numYoshidaSuzuki) {
    result = OpenMM_NoseHooverIntegrator_create_2(temperature, collisionFrequency, stepSize, chainLength, numMTS, numYoshidaSuzuki);
}
OPENMM_EXPORT void OPENMM_NOSEHOOVERINTEGRATOR_CREATE_2(OpenMM_NoseHooverIntegrator*& result, double const& temperature, double const& collisionFrequency, double const& stepSize, int const& chainLength, int const& numMTS, int const& numYoshidaSuzuki) {
    result = OpenMM_NoseHooverIntegrator_create_2(temperature, collisionFrequency, stepSize, chainLength, numMTS, numYoshidaSuzuki);
}
OPENMM_EXPORT void openmm_nosehooverintegrator_destroy_(OpenMM_NoseHooverIntegrator*& destroy) {
    OpenMM_NoseHooverIntegrator_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_NOSEHOOVERINTEGRATOR_DESTROY(OpenMM_NoseHooverIntegrator*& destroy) {
    OpenMM_NoseHooverIntegrator_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void openmm_nosehooverintegrator_step_(OpenMM_NoseHooverIntegrator*& target, int const& steps) {
    OpenMM_NoseHooverIntegrator_step(target, steps);
}
OPENMM_EXPORT void OPENMM_NOSEHOOVERINTEGRATOR_STEP(OpenMM_NoseHooverIntegrator*& target, int const& steps) {
    OpenMM_NoseHooverIntegrator_step(target, steps);
}
OPENMM_EXPORT int openmm_nosehooverintegrator_addthermostat_(OpenMM_NoseHooverIntegrator*& target, double const& temperature, double const& collisionFrequency, int const& chainLength, int const& numMTS, int const& numYoshidaSuzuki) {
    return OpenMM_NoseHooverIntegrator_addThermostat(target, temperature, collisionFrequency, chainLength, numMTS, numYoshidaSuzuki);
}
OPENMM_EXPORT int OPENMM_NOSEHOOVERINTEGRATOR_ADDTHERMOSTAT(OpenMM_NoseHooverIntegrator*& target, double const& temperature, double const& collisionFrequency, int const& chainLength, int const& numMTS, int const& numYoshidaSuzuki) {
    return OpenMM_NoseHooverIntegrator_addThermostat(target, temperature, collisionFrequency, chainLength, numMTS, numYoshidaSuzuki);
}
OPENMM_EXPORT int openmm_nosehooverintegrator_addsubsystemthermostat_(OpenMM_NoseHooverIntegrator*& target, const OpenMM_IntArray*& thermostatedParticles, const OpenMM_BondArray*& thermostatedPairs, double const& temperature, double const& collisionFrequency, double const& relativeTemperature, double const& relativeCollisionFrequency, int const& chainLength, int const& numMTS, int const& numYoshidaSuzuki) {
    return OpenMM_NoseHooverIntegrator_addSubsystemThermostat(target, thermostatedParticles, thermostatedPairs, temperature, collisionFrequency, relativeTemperature, relativeCollisionFrequency, chainLength, numMTS, numYoshidaSuzuki);
}
OPENMM_EXPORT int OPENMM_NOSEHOOVERINTEGRATOR_ADDSUBSYSTEMTHERMOSTAT(OpenMM_NoseHooverIntegrator*& target, const OpenMM_IntArray*& thermostatedParticles, const OpenMM_BondArray*& thermostatedPairs, double const& temperature, double const& collisionFrequency, double const& relativeTemperature, double const& relativeCollisionFrequency, int const& chainLength, int const& numMTS, int const& numYoshidaSuzuki) {
    return OpenMM_NoseHooverIntegrator_addSubsystemThermostat(target, thermostatedParticles, thermostatedPairs, temperature, collisionFrequency, relativeTemperature, relativeCollisionFrequency, chainLength, numMTS, numYoshidaSuzuki);
}
OPENMM_EXPORT double openmm_nosehooverintegrator_gettemperature_(const OpenMM_NoseHooverIntegrator*& target, int const& chainID) {
    return OpenMM_NoseHooverIntegrator_getTemperature(target, chainID);
}
OPENMM_EXPORT double OPENMM_NOSEHOOVERINTEGRATOR_GETTEMPERATURE(const OpenMM_NoseHooverIntegrator*& target, int const& chainID) {
    return OpenMM_NoseHooverIntegrator_getTemperature(target, chainID);
}
OPENMM_EXPORT void openmm_nosehooverintegrator_settemperature_(OpenMM_NoseHooverIntegrator*& target, double const& temperature, int const& chainID) {
    OpenMM_NoseHooverIntegrator_setTemperature(target, temperature, chainID);
}
OPENMM_EXPORT void OPENMM_NOSEHOOVERINTEGRATOR_SETTEMPERATURE(OpenMM_NoseHooverIntegrator*& target, double const& temperature, int const& chainID) {
    OpenMM_NoseHooverIntegrator_setTemperature(target, temperature, chainID);
}
OPENMM_EXPORT double openmm_nosehooverintegrator_getrelativetemperature_(const OpenMM_NoseHooverIntegrator*& target, int const& chainID) {
    return OpenMM_NoseHooverIntegrator_getRelativeTemperature(target, chainID);
}
OPENMM_EXPORT double OPENMM_NOSEHOOVERINTEGRATOR_GETRELATIVETEMPERATURE(const OpenMM_NoseHooverIntegrator*& target, int const& chainID) {
    return OpenMM_NoseHooverIntegrator_getRelativeTemperature(target, chainID);
}
OPENMM_EXPORT void openmm_nosehooverintegrator_setrelativetemperature_(OpenMM_NoseHooverIntegrator*& target, double const& temperature, int const& chainID) {
    OpenMM_NoseHooverIntegrator_setRelativeTemperature(target, temperature, chainID);
}
OPENMM_EXPORT void OPENMM_NOSEHOOVERINTEGRATOR_SETRELATIVETEMPERATURE(OpenMM_NoseHooverIntegrator*& target, double const& temperature, int const& chainID) {
    OpenMM_NoseHooverIntegrator_setRelativeTemperature(target, temperature, chainID);
}
OPENMM_EXPORT double openmm_nosehooverintegrator_getcollisionfrequency_(const OpenMM_NoseHooverIntegrator*& target, int const& chainID) {
    return OpenMM_NoseHooverIntegrator_getCollisionFrequency(target, chainID);
}
OPENMM_EXPORT double OPENMM_NOSEHOOVERINTEGRATOR_GETCOLLISIONFREQUENCY(const OpenMM_NoseHooverIntegrator*& target, int const& chainID) {
    return OpenMM_NoseHooverIntegrator_getCollisionFrequency(target, chainID);
}
OPENMM_EXPORT void openmm_nosehooverintegrator_setcollisionfrequency_(OpenMM_NoseHooverIntegrator*& target, double const& frequency, int const& chainID) {
    OpenMM_NoseHooverIntegrator_setCollisionFrequency(target, frequency, chainID);
}
OPENMM_EXPORT void OPENMM_NOSEHOOVERINTEGRATOR_SETCOLLISIONFREQUENCY(OpenMM_NoseHooverIntegrator*& target, double const& frequency, int const& chainID) {
    OpenMM_NoseHooverIntegrator_setCollisionFrequency(target, frequency, chainID);
}
OPENMM_EXPORT double openmm_nosehooverintegrator_getrelativecollisionfrequency_(const OpenMM_NoseHooverIntegrator*& target, int const& chainID) {
    return OpenMM_NoseHooverIntegrator_getRelativeCollisionFrequency(target, chainID);
}
OPENMM_EXPORT double OPENMM_NOSEHOOVERINTEGRATOR_GETRELATIVECOLLISIONFREQUENCY(const OpenMM_NoseHooverIntegrator*& target, int const& chainID) {
    return OpenMM_NoseHooverIntegrator_getRelativeCollisionFrequency(target, chainID);
}
OPENMM_EXPORT void openmm_nosehooverintegrator_setrelativecollisionfrequency_(OpenMM_NoseHooverIntegrator*& target, double const& frequency, int const& chainID) {
    OpenMM_NoseHooverIntegrator_setRelativeCollisionFrequency(target, frequency, chainID);
}
OPENMM_EXPORT void OPENMM_NOSEHOOVERINTEGRATOR_SETRELATIVECOLLISIONFREQUENCY(OpenMM_NoseHooverIntegrator*& target, double const& frequency, int const& chainID) {
    OpenMM_NoseHooverIntegrator_setRelativeCollisionFrequency(target, frequency, chainID);
}
OPENMM_EXPORT double openmm_nosehooverintegrator_computeheatbathenergy_(OpenMM_NoseHooverIntegrator*& target) {
    return OpenMM_NoseHooverIntegrator_computeHeatBathEnergy(target);
}
OPENMM_EXPORT double OPENMM_NOSEHOOVERINTEGRATOR_COMPUTEHEATBATHENERGY(OpenMM_NoseHooverIntegrator*& target) {
    return OpenMM_NoseHooverIntegrator_computeHeatBathEnergy(target);
}
OPENMM_EXPORT int openmm_nosehooverintegrator_getnumthermostats_(const OpenMM_NoseHooverIntegrator*& target) {
    return OpenMM_NoseHooverIntegrator_getNumThermostats(target);
}
OPENMM_EXPORT int OPENMM_NOSEHOOVERINTEGRATOR_GETNUMTHERMOSTATS(const OpenMM_NoseHooverIntegrator*& target) {
    return OpenMM_NoseHooverIntegrator_getNumThermostats(target);
}
OPENMM_EXPORT void openmm_nosehooverintegrator_getthermostat_(const OpenMM_NoseHooverIntegrator*& target, int const& chainID, const OpenMM_NoseHooverChain*& result) {
    result = OpenMM_NoseHooverIntegrator_getThermostat(target, chainID);
}
OPENMM_EXPORT void OPENMM_NOSEHOOVERINTEGRATOR_GETTHERMOSTAT(const OpenMM_NoseHooverIntegrator*& target, int const& chainID, const OpenMM_NoseHooverChain*& result) {
    result = OpenMM_NoseHooverIntegrator_getThermostat(target, chainID);
}
OPENMM_EXPORT void openmm_nosehooverintegrator_hassubsystemthermostats_(const OpenMM_NoseHooverIntegrator*& target, OpenMM_Boolean& result) {
    result = OpenMM_NoseHooverIntegrator_hasSubsystemThermostats(target);
}
OPENMM_EXPORT void OPENMM_NOSEHOOVERINTEGRATOR_HASSUBSYSTEMTHERMOSTATS(const OpenMM_NoseHooverIntegrator*& target, OpenMM_Boolean& result) {
    result = OpenMM_NoseHooverIntegrator_hasSubsystemThermostats(target);
}
OPENMM_EXPORT double openmm_nosehooverintegrator_getmaximumpairdistance_(const OpenMM_NoseHooverIntegrator*& target) {
    return OpenMM_NoseHooverIntegrator_getMaximumPairDistance(target);
}
OPENMM_EXPORT double OPENMM_NOSEHOOVERINTEGRATOR_GETMAXIMUMPAIRDISTANCE(const OpenMM_NoseHooverIntegrator*& target) {
    return OpenMM_NoseHooverIntegrator_getMaximumPairDistance(target);
}
OPENMM_EXPORT void openmm_nosehooverintegrator_setmaximumpairdistance_(OpenMM_NoseHooverIntegrator*& target, double const& distance) {
    OpenMM_NoseHooverIntegrator_setMaximumPairDistance(target, distance);
}
OPENMM_EXPORT void OPENMM_NOSEHOOVERINTEGRATOR_SETMAXIMUMPAIRDISTANCE(OpenMM_NoseHooverIntegrator*& target, double const& distance) {
    OpenMM_NoseHooverIntegrator_setMaximumPairDistance(target, distance);
}

/* OpenMM::VariableVerletIntegrator */
OPENMM_EXPORT void openmm_variableverletintegrator_create_(OpenMM_VariableVerletIntegrator*& result, double const& errorTol) {
    result = OpenMM_VariableVerletIntegrator_create(errorTol);
}
OPENMM_EXPORT void OPENMM_VARIABLEVERLETINTEGRATOR_CREATE(OpenMM_VariableVerletIntegrator*& result, double const& errorTol) {
    result = OpenMM_VariableVerletIntegrator_create(errorTol);
}
OPENMM_EXPORT void openmm_variableverletintegrator_destroy_(OpenMM_VariableVerletIntegrator*& destroy) {
    OpenMM_VariableVerletIntegrator_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_VARIABLEVERLETINTEGRATOR_DESTROY(OpenMM_VariableVerletIntegrator*& destroy) {
    OpenMM_VariableVerletIntegrator_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT double openmm_variableverletintegrator_geterrortolerance_(const OpenMM_VariableVerletIntegrator*& target) {
    return OpenMM_VariableVerletIntegrator_getErrorTolerance(target);
}
OPENMM_EXPORT double OPENMM_VARIABLEVERLETINTEGRATOR_GETERRORTOLERANCE(const OpenMM_VariableVerletIntegrator*& target) {
    return OpenMM_VariableVerletIntegrator_getErrorTolerance(target);
}
OPENMM_EXPORT void openmm_variableverletintegrator_seterrortolerance_(OpenMM_VariableVerletIntegrator*& target, double const& tol) {
    OpenMM_VariableVerletIntegrator_setErrorTolerance(target, tol);
}
OPENMM_EXPORT void OPENMM_VARIABLEVERLETINTEGRATOR_SETERRORTOLERANCE(OpenMM_VariableVerletIntegrator*& target, double const& tol) {
    OpenMM_VariableVerletIntegrator_setErrorTolerance(target, tol);
}
OPENMM_EXPORT double openmm_variableverletintegrator_getmaximumstepsize_(const OpenMM_VariableVerletIntegrator*& target) {
    return OpenMM_VariableVerletIntegrator_getMaximumStepSize(target);
}
OPENMM_EXPORT double OPENMM_VARIABLEVERLETINTEGRATOR_GETMAXIMUMSTEPSIZE(const OpenMM_VariableVerletIntegrator*& target) {
    return OpenMM_VariableVerletIntegrator_getMaximumStepSize(target);
}
OPENMM_EXPORT void openmm_variableverletintegrator_setmaximumstepsize_(OpenMM_VariableVerletIntegrator*& target, double const& size) {
    OpenMM_VariableVerletIntegrator_setMaximumStepSize(target, size);
}
OPENMM_EXPORT void OPENMM_VARIABLEVERLETINTEGRATOR_SETMAXIMUMSTEPSIZE(OpenMM_VariableVerletIntegrator*& target, double const& size) {
    OpenMM_VariableVerletIntegrator_setMaximumStepSize(target, size);
}
OPENMM_EXPORT void openmm_variableverletintegrator_step_(OpenMM_VariableVerletIntegrator*& target, int const& steps) {
    OpenMM_VariableVerletIntegrator_step(target, steps);
}
OPENMM_EXPORT void OPENMM_VARIABLEVERLETINTEGRATOR_STEP(OpenMM_VariableVerletIntegrator*& target, int const& steps) {
    OpenMM_VariableVerletIntegrator_step(target, steps);
}
OPENMM_EXPORT void openmm_variableverletintegrator_stepto_(OpenMM_VariableVerletIntegrator*& target, double const& time) {
    OpenMM_VariableVerletIntegrator_stepTo(target, time);
}
OPENMM_EXPORT void OPENMM_VARIABLEVERLETINTEGRATOR_STEPTO(OpenMM_VariableVerletIntegrator*& target, double const& time) {
    OpenMM_VariableVerletIntegrator_stepTo(target, time);
}

/* OpenMM::CustomHbondForce */
OPENMM_EXPORT void openmm_customhbondforce_create_(OpenMM_CustomHbondForce*& result, const char* energy, int energy_length) {
    result = OpenMM_CustomHbondForce_create(makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMHBONDFORCE_CREATE(OpenMM_CustomHbondForce*& result, const char* energy, int energy_length) {
    result = OpenMM_CustomHbondForce_create(makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void openmm_customhbondforce_destroy_(OpenMM_CustomHbondForce*& destroy) {
    OpenMM_CustomHbondForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_CUSTOMHBONDFORCE_DESTROY(OpenMM_CustomHbondForce*& destroy) {
    OpenMM_CustomHbondForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT int openmm_customhbondforce_getnumdonors_(const OpenMM_CustomHbondForce*& target) {
    return OpenMM_CustomHbondForce_getNumDonors(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMHBONDFORCE_GETNUMDONORS(const OpenMM_CustomHbondForce*& target) {
    return OpenMM_CustomHbondForce_getNumDonors(target);
}
OPENMM_EXPORT int openmm_customhbondforce_getnumacceptors_(const OpenMM_CustomHbondForce*& target) {
    return OpenMM_CustomHbondForce_getNumAcceptors(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMHBONDFORCE_GETNUMACCEPTORS(const OpenMM_CustomHbondForce*& target) {
    return OpenMM_CustomHbondForce_getNumAcceptors(target);
}
OPENMM_EXPORT int openmm_customhbondforce_getnumexclusions_(const OpenMM_CustomHbondForce*& target) {
    return OpenMM_CustomHbondForce_getNumExclusions(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMHBONDFORCE_GETNUMEXCLUSIONS(const OpenMM_CustomHbondForce*& target) {
    return OpenMM_CustomHbondForce_getNumExclusions(target);
}
OPENMM_EXPORT int openmm_customhbondforce_getnumperdonorparameters_(const OpenMM_CustomHbondForce*& target) {
    return OpenMM_CustomHbondForce_getNumPerDonorParameters(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMHBONDFORCE_GETNUMPERDONORPARAMETERS(const OpenMM_CustomHbondForce*& target) {
    return OpenMM_CustomHbondForce_getNumPerDonorParameters(target);
}
OPENMM_EXPORT int openmm_customhbondforce_getnumperacceptorparameters_(const OpenMM_CustomHbondForce*& target) {
    return OpenMM_CustomHbondForce_getNumPerAcceptorParameters(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMHBONDFORCE_GETNUMPERACCEPTORPARAMETERS(const OpenMM_CustomHbondForce*& target) {
    return OpenMM_CustomHbondForce_getNumPerAcceptorParameters(target);
}
OPENMM_EXPORT int openmm_customhbondforce_getnumglobalparameters_(const OpenMM_CustomHbondForce*& target) {
    return OpenMM_CustomHbondForce_getNumGlobalParameters(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMHBONDFORCE_GETNUMGLOBALPARAMETERS(const OpenMM_CustomHbondForce*& target) {
    return OpenMM_CustomHbondForce_getNumGlobalParameters(target);
}
OPENMM_EXPORT int openmm_customhbondforce_getnumtabulatedfunctions_(const OpenMM_CustomHbondForce*& target) {
    return OpenMM_CustomHbondForce_getNumTabulatedFunctions(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMHBONDFORCE_GETNUMTABULATEDFUNCTIONS(const OpenMM_CustomHbondForce*& target) {
    return OpenMM_CustomHbondForce_getNumTabulatedFunctions(target);
}
OPENMM_EXPORT int openmm_customhbondforce_getnumfunctions_(const OpenMM_CustomHbondForce*& target) {
    return OpenMM_CustomHbondForce_getNumFunctions(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMHBONDFORCE_GETNUMFUNCTIONS(const OpenMM_CustomHbondForce*& target) {
    return OpenMM_CustomHbondForce_getNumFunctions(target);
}
OPENMM_EXPORT void openmm_customhbondforce_getenergyfunction_(const OpenMM_CustomHbondForce*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomHbondForce_getEnergyFunction(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMHBONDFORCE_GETENERGYFUNCTION(const OpenMM_CustomHbondForce*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomHbondForce_getEnergyFunction(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customhbondforce_setenergyfunction_(OpenMM_CustomHbondForce*& target, const char* energy, int energy_length) {
    OpenMM_CustomHbondForce_setEnergyFunction(target, makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMHBONDFORCE_SETENERGYFUNCTION(OpenMM_CustomHbondForce*& target, const char* energy, int energy_length) {
    OpenMM_CustomHbondForce_setEnergyFunction(target, makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void openmm_customhbondforce_getnonbondedmethod_(const OpenMM_CustomHbondForce*& target, OpenMM_CustomHbondForce_NonbondedMethod& result) {
    result = OpenMM_CustomHbondForce_getNonbondedMethod(target);
}
OPENMM_EXPORT void OPENMM_CUSTOMHBONDFORCE_GETNONBONDEDMETHOD(const OpenMM_CustomHbondForce*& target, OpenMM_CustomHbondForce_NonbondedMethod& result) {
    result = OpenMM_CustomHbondForce_getNonbondedMethod(target);
}
OPENMM_EXPORT void openmm_customhbondforce_setnonbondedmethod_(OpenMM_CustomHbondForce*& target, OpenMM_CustomHbondForce_NonbondedMethod& method) {
    OpenMM_CustomHbondForce_setNonbondedMethod(target, method);
}
OPENMM_EXPORT void OPENMM_CUSTOMHBONDFORCE_SETNONBONDEDMETHOD(OpenMM_CustomHbondForce*& target, OpenMM_CustomHbondForce_NonbondedMethod& method) {
    OpenMM_CustomHbondForce_setNonbondedMethod(target, method);
}
OPENMM_EXPORT double openmm_customhbondforce_getcutoffdistance_(const OpenMM_CustomHbondForce*& target) {
    return OpenMM_CustomHbondForce_getCutoffDistance(target);
}
OPENMM_EXPORT double OPENMM_CUSTOMHBONDFORCE_GETCUTOFFDISTANCE(const OpenMM_CustomHbondForce*& target) {
    return OpenMM_CustomHbondForce_getCutoffDistance(target);
}
OPENMM_EXPORT void openmm_customhbondforce_setcutoffdistance_(OpenMM_CustomHbondForce*& target, double const& distance) {
    OpenMM_CustomHbondForce_setCutoffDistance(target, distance);
}
OPENMM_EXPORT void OPENMM_CUSTOMHBONDFORCE_SETCUTOFFDISTANCE(OpenMM_CustomHbondForce*& target, double const& distance) {
    OpenMM_CustomHbondForce_setCutoffDistance(target, distance);
}
OPENMM_EXPORT int openmm_customhbondforce_addperdonorparameter_(OpenMM_CustomHbondForce*& target, const char* name, int name_length) {
    return OpenMM_CustomHbondForce_addPerDonorParameter(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT int OPENMM_CUSTOMHBONDFORCE_ADDPERDONORPARAMETER(OpenMM_CustomHbondForce*& target, const char* name, int name_length) {
    return OpenMM_CustomHbondForce_addPerDonorParameter(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void openmm_customhbondforce_getperdonorparametername_(const OpenMM_CustomHbondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomHbondForce_getPerDonorParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMHBONDFORCE_GETPERDONORPARAMETERNAME(const OpenMM_CustomHbondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomHbondForce_getPerDonorParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customhbondforce_setperdonorparametername_(OpenMM_CustomHbondForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomHbondForce_setPerDonorParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMHBONDFORCE_SETPERDONORPARAMETERNAME(OpenMM_CustomHbondForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomHbondForce_setPerDonorParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT int openmm_customhbondforce_addperacceptorparameter_(OpenMM_CustomHbondForce*& target, const char* name, int name_length) {
    return OpenMM_CustomHbondForce_addPerAcceptorParameter(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT int OPENMM_CUSTOMHBONDFORCE_ADDPERACCEPTORPARAMETER(OpenMM_CustomHbondForce*& target, const char* name, int name_length) {
    return OpenMM_CustomHbondForce_addPerAcceptorParameter(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void openmm_customhbondforce_getperacceptorparametername_(const OpenMM_CustomHbondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomHbondForce_getPerAcceptorParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMHBONDFORCE_GETPERACCEPTORPARAMETERNAME(const OpenMM_CustomHbondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomHbondForce_getPerAcceptorParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customhbondforce_setperacceptorparametername_(OpenMM_CustomHbondForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomHbondForce_setPerAcceptorParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMHBONDFORCE_SETPERACCEPTORPARAMETERNAME(OpenMM_CustomHbondForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomHbondForce_setPerAcceptorParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT int openmm_customhbondforce_addglobalparameter_(OpenMM_CustomHbondForce*& target, const char* name, double const& defaultValue, int name_length) {
    return OpenMM_CustomHbondForce_addGlobalParameter(target, makeString(name, name_length).c_str(), defaultValue);
}
OPENMM_EXPORT int OPENMM_CUSTOMHBONDFORCE_ADDGLOBALPARAMETER(OpenMM_CustomHbondForce*& target, const char* name, double const& defaultValue, int name_length) {
    return OpenMM_CustomHbondForce_addGlobalParameter(target, makeString(name, name_length).c_str(), defaultValue);
}
OPENMM_EXPORT void openmm_customhbondforce_getglobalparametername_(const OpenMM_CustomHbondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomHbondForce_getGlobalParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMHBONDFORCE_GETGLOBALPARAMETERNAME(const OpenMM_CustomHbondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomHbondForce_getGlobalParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customhbondforce_setglobalparametername_(OpenMM_CustomHbondForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomHbondForce_setGlobalParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMHBONDFORCE_SETGLOBALPARAMETERNAME(OpenMM_CustomHbondForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomHbondForce_setGlobalParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT double openmm_customhbondforce_getglobalparameterdefaultvalue_(const OpenMM_CustomHbondForce*& target, int const& index) {
    return OpenMM_CustomHbondForce_getGlobalParameterDefaultValue(target, index);
}
OPENMM_EXPORT double OPENMM_CUSTOMHBONDFORCE_GETGLOBALPARAMETERDEFAULTVALUE(const OpenMM_CustomHbondForce*& target, int const& index) {
    return OpenMM_CustomHbondForce_getGlobalParameterDefaultValue(target, index);
}
OPENMM_EXPORT void openmm_customhbondforce_setglobalparameterdefaultvalue_(OpenMM_CustomHbondForce*& target, int const& index, double const& defaultValue) {
    OpenMM_CustomHbondForce_setGlobalParameterDefaultValue(target, index, defaultValue);
}
OPENMM_EXPORT void OPENMM_CUSTOMHBONDFORCE_SETGLOBALPARAMETERDEFAULTVALUE(OpenMM_CustomHbondForce*& target, int const& index, double const& defaultValue) {
    OpenMM_CustomHbondForce_setGlobalParameterDefaultValue(target, index, defaultValue);
}
OPENMM_EXPORT int openmm_customhbondforce_adddonor_(OpenMM_CustomHbondForce*& target, int const& d1, int const& d2, int const& d3, const OpenMM_DoubleArray*& parameters) {
    return OpenMM_CustomHbondForce_addDonor(target, d1, d2, d3, parameters);
}
OPENMM_EXPORT int OPENMM_CUSTOMHBONDFORCE_ADDDONOR(OpenMM_CustomHbondForce*& target, int const& d1, int const& d2, int const& d3, const OpenMM_DoubleArray*& parameters) {
    return OpenMM_CustomHbondForce_addDonor(target, d1, d2, d3, parameters);
}
OPENMM_EXPORT void openmm_customhbondforce_getdonorparameters_(const OpenMM_CustomHbondForce*& target, int const& index, int* d1, int* d2, int* d3, OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomHbondForce_getDonorParameters(target, index, d1, d2, d3, parameters);
}
OPENMM_EXPORT void OPENMM_CUSTOMHBONDFORCE_GETDONORPARAMETERS(const OpenMM_CustomHbondForce*& target, int const& index, int* d1, int* d2, int* d3, OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomHbondForce_getDonorParameters(target, index, d1, d2, d3, parameters);
}
OPENMM_EXPORT void openmm_customhbondforce_setdonorparameters_(OpenMM_CustomHbondForce*& target, int const& index, int const& d1, int const& d2, int const& d3, const OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomHbondForce_setDonorParameters(target, index, d1, d2, d3, parameters);
}
OPENMM_EXPORT void OPENMM_CUSTOMHBONDFORCE_SETDONORPARAMETERS(OpenMM_CustomHbondForce*& target, int const& index, int const& d1, int const& d2, int const& d3, const OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomHbondForce_setDonorParameters(target, index, d1, d2, d3, parameters);
}
OPENMM_EXPORT int openmm_customhbondforce_addacceptor_(OpenMM_CustomHbondForce*& target, int const& a1, int const& a2, int const& a3, const OpenMM_DoubleArray*& parameters) {
    return OpenMM_CustomHbondForce_addAcceptor(target, a1, a2, a3, parameters);
}
OPENMM_EXPORT int OPENMM_CUSTOMHBONDFORCE_ADDACCEPTOR(OpenMM_CustomHbondForce*& target, int const& a1, int const& a2, int const& a3, const OpenMM_DoubleArray*& parameters) {
    return OpenMM_CustomHbondForce_addAcceptor(target, a1, a2, a3, parameters);
}
OPENMM_EXPORT void openmm_customhbondforce_getacceptorparameters_(const OpenMM_CustomHbondForce*& target, int const& index, int* a1, int* a2, int* a3, OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomHbondForce_getAcceptorParameters(target, index, a1, a2, a3, parameters);
}
OPENMM_EXPORT void OPENMM_CUSTOMHBONDFORCE_GETACCEPTORPARAMETERS(const OpenMM_CustomHbondForce*& target, int const& index, int* a1, int* a2, int* a3, OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomHbondForce_getAcceptorParameters(target, index, a1, a2, a3, parameters);
}
OPENMM_EXPORT void openmm_customhbondforce_setacceptorparameters_(OpenMM_CustomHbondForce*& target, int const& index, int const& a1, int const& a2, int const& a3, const OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomHbondForce_setAcceptorParameters(target, index, a1, a2, a3, parameters);
}
OPENMM_EXPORT void OPENMM_CUSTOMHBONDFORCE_SETACCEPTORPARAMETERS(OpenMM_CustomHbondForce*& target, int const& index, int const& a1, int const& a2, int const& a3, const OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomHbondForce_setAcceptorParameters(target, index, a1, a2, a3, parameters);
}
OPENMM_EXPORT int openmm_customhbondforce_addexclusion_(OpenMM_CustomHbondForce*& target, int const& donor, int const& acceptor) {
    return OpenMM_CustomHbondForce_addExclusion(target, donor, acceptor);
}
OPENMM_EXPORT int OPENMM_CUSTOMHBONDFORCE_ADDEXCLUSION(OpenMM_CustomHbondForce*& target, int const& donor, int const& acceptor) {
    return OpenMM_CustomHbondForce_addExclusion(target, donor, acceptor);
}
OPENMM_EXPORT void openmm_customhbondforce_getexclusionparticles_(const OpenMM_CustomHbondForce*& target, int const& index, int* donor, int* acceptor) {
    OpenMM_CustomHbondForce_getExclusionParticles(target, index, donor, acceptor);
}
OPENMM_EXPORT void OPENMM_CUSTOMHBONDFORCE_GETEXCLUSIONPARTICLES(const OpenMM_CustomHbondForce*& target, int const& index, int* donor, int* acceptor) {
    OpenMM_CustomHbondForce_getExclusionParticles(target, index, donor, acceptor);
}
OPENMM_EXPORT void openmm_customhbondforce_setexclusionparticles_(OpenMM_CustomHbondForce*& target, int const& index, int const& donor, int const& acceptor) {
    OpenMM_CustomHbondForce_setExclusionParticles(target, index, donor, acceptor);
}
OPENMM_EXPORT void OPENMM_CUSTOMHBONDFORCE_SETEXCLUSIONPARTICLES(OpenMM_CustomHbondForce*& target, int const& index, int const& donor, int const& acceptor) {
    OpenMM_CustomHbondForce_setExclusionParticles(target, index, donor, acceptor);
}
OPENMM_EXPORT int openmm_customhbondforce_addtabulatedfunction_(OpenMM_CustomHbondForce*& target, const char* name, OpenMM_TabulatedFunction*& function, int name_length) {
    return OpenMM_CustomHbondForce_addTabulatedFunction(target, makeString(name, name_length).c_str(), function);
}
OPENMM_EXPORT int OPENMM_CUSTOMHBONDFORCE_ADDTABULATEDFUNCTION(OpenMM_CustomHbondForce*& target, const char* name, OpenMM_TabulatedFunction*& function, int name_length) {
    return OpenMM_CustomHbondForce_addTabulatedFunction(target, makeString(name, name_length).c_str(), function);
}
OPENMM_EXPORT void openmm_customhbondforce_gettabulatedfunction_(OpenMM_CustomHbondForce*& target, int const& index, OpenMM_TabulatedFunction*& result) {
    result = OpenMM_CustomHbondForce_getTabulatedFunction(target, index);
}
OPENMM_EXPORT void OPENMM_CUSTOMHBONDFORCE_GETTABULATEDFUNCTION(OpenMM_CustomHbondForce*& target, int const& index, OpenMM_TabulatedFunction*& result) {
    result = OpenMM_CustomHbondForce_getTabulatedFunction(target, index);
}
OPENMM_EXPORT void openmm_customhbondforce_gettabulatedfunctionname_(const OpenMM_CustomHbondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomHbondForce_getTabulatedFunctionName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMHBONDFORCE_GETTABULATEDFUNCTIONNAME(const OpenMM_CustomHbondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomHbondForce_getTabulatedFunctionName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT int openmm_customhbondforce_addfunction_(OpenMM_CustomHbondForce*& target, const char* name, const OpenMM_DoubleArray*& values, double const& min, double const& max, int name_length) {
    return OpenMM_CustomHbondForce_addFunction(target, makeString(name, name_length).c_str(), values, min, max);
}
OPENMM_EXPORT int OPENMM_CUSTOMHBONDFORCE_ADDFUNCTION(OpenMM_CustomHbondForce*& target, const char* name, const OpenMM_DoubleArray*& values, double const& min, double const& max, int name_length) {
    return OpenMM_CustomHbondForce_addFunction(target, makeString(name, name_length).c_str(), values, min, max);
}
OPENMM_EXPORT void openmm_customhbondforce_getfunctionparameters_(const OpenMM_CustomHbondForce*& target, int const& index, char** name, OpenMM_DoubleArray*& values, double* min, double* max) {
    OpenMM_CustomHbondForce_getFunctionParameters(target, index, name, values, min, max);
}
OPENMM_EXPORT void OPENMM_CUSTOMHBONDFORCE_GETFUNCTIONPARAMETERS(const OpenMM_CustomHbondForce*& target, int const& index, char** name, OpenMM_DoubleArray*& values, double* min, double* max) {
    OpenMM_CustomHbondForce_getFunctionParameters(target, index, name, values, min, max);
}
OPENMM_EXPORT void openmm_customhbondforce_setfunctionparameters_(OpenMM_CustomHbondForce*& target, int const& index, const char* name, const OpenMM_DoubleArray*& values, double const& min, double const& max, int name_length) {
    OpenMM_CustomHbondForce_setFunctionParameters(target, index, makeString(name, name_length).c_str(), values, min, max);
}
OPENMM_EXPORT void OPENMM_CUSTOMHBONDFORCE_SETFUNCTIONPARAMETERS(OpenMM_CustomHbondForce*& target, int const& index, const char* name, const OpenMM_DoubleArray*& values, double const& min, double const& max, int name_length) {
    OpenMM_CustomHbondForce_setFunctionParameters(target, index, makeString(name, name_length).c_str(), values, min, max);
}
OPENMM_EXPORT void openmm_customhbondforce_updateparametersincontext_(OpenMM_CustomHbondForce*& target, OpenMM_Context*& context) {
    OpenMM_CustomHbondForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void OPENMM_CUSTOMHBONDFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_CustomHbondForce*& target, OpenMM_Context*& context) {
    OpenMM_CustomHbondForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void openmm_customhbondforce_usesperiodicboundaryconditions_(const OpenMM_CustomHbondForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CustomHbondForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_CUSTOMHBONDFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_CustomHbondForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CustomHbondForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::CustomExternalForce */
OPENMM_EXPORT void openmm_customexternalforce_create_(OpenMM_CustomExternalForce*& result, const char* energy, int energy_length) {
    result = OpenMM_CustomExternalForce_create(makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMEXTERNALFORCE_CREATE(OpenMM_CustomExternalForce*& result, const char* energy, int energy_length) {
    result = OpenMM_CustomExternalForce_create(makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void openmm_customexternalforce_destroy_(OpenMM_CustomExternalForce*& destroy) {
    OpenMM_CustomExternalForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_CUSTOMEXTERNALFORCE_DESTROY(OpenMM_CustomExternalForce*& destroy) {
    OpenMM_CustomExternalForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT int openmm_customexternalforce_getnumparticles_(const OpenMM_CustomExternalForce*& target) {
    return OpenMM_CustomExternalForce_getNumParticles(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMEXTERNALFORCE_GETNUMPARTICLES(const OpenMM_CustomExternalForce*& target) {
    return OpenMM_CustomExternalForce_getNumParticles(target);
}
OPENMM_EXPORT int openmm_customexternalforce_getnumperparticleparameters_(const OpenMM_CustomExternalForce*& target) {
    return OpenMM_CustomExternalForce_getNumPerParticleParameters(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMEXTERNALFORCE_GETNUMPERPARTICLEPARAMETERS(const OpenMM_CustomExternalForce*& target) {
    return OpenMM_CustomExternalForce_getNumPerParticleParameters(target);
}
OPENMM_EXPORT int openmm_customexternalforce_getnumglobalparameters_(const OpenMM_CustomExternalForce*& target) {
    return OpenMM_CustomExternalForce_getNumGlobalParameters(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMEXTERNALFORCE_GETNUMGLOBALPARAMETERS(const OpenMM_CustomExternalForce*& target) {
    return OpenMM_CustomExternalForce_getNumGlobalParameters(target);
}
OPENMM_EXPORT void openmm_customexternalforce_getenergyfunction_(const OpenMM_CustomExternalForce*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomExternalForce_getEnergyFunction(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMEXTERNALFORCE_GETENERGYFUNCTION(const OpenMM_CustomExternalForce*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomExternalForce_getEnergyFunction(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customexternalforce_setenergyfunction_(OpenMM_CustomExternalForce*& target, const char* energy, int energy_length) {
    OpenMM_CustomExternalForce_setEnergyFunction(target, makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMEXTERNALFORCE_SETENERGYFUNCTION(OpenMM_CustomExternalForce*& target, const char* energy, int energy_length) {
    OpenMM_CustomExternalForce_setEnergyFunction(target, makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT int openmm_customexternalforce_addperparticleparameter_(OpenMM_CustomExternalForce*& target, const char* name, int name_length) {
    return OpenMM_CustomExternalForce_addPerParticleParameter(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT int OPENMM_CUSTOMEXTERNALFORCE_ADDPERPARTICLEPARAMETER(OpenMM_CustomExternalForce*& target, const char* name, int name_length) {
    return OpenMM_CustomExternalForce_addPerParticleParameter(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void openmm_customexternalforce_getperparticleparametername_(const OpenMM_CustomExternalForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomExternalForce_getPerParticleParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMEXTERNALFORCE_GETPERPARTICLEPARAMETERNAME(const OpenMM_CustomExternalForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomExternalForce_getPerParticleParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customexternalforce_setperparticleparametername_(OpenMM_CustomExternalForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomExternalForce_setPerParticleParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMEXTERNALFORCE_SETPERPARTICLEPARAMETERNAME(OpenMM_CustomExternalForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomExternalForce_setPerParticleParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT int openmm_customexternalforce_addglobalparameter_(OpenMM_CustomExternalForce*& target, const char* name, double const& defaultValue, int name_length) {
    return OpenMM_CustomExternalForce_addGlobalParameter(target, makeString(name, name_length).c_str(), defaultValue);
}
OPENMM_EXPORT int OPENMM_CUSTOMEXTERNALFORCE_ADDGLOBALPARAMETER(OpenMM_CustomExternalForce*& target, const char* name, double const& defaultValue, int name_length) {
    return OpenMM_CustomExternalForce_addGlobalParameter(target, makeString(name, name_length).c_str(), defaultValue);
}
OPENMM_EXPORT void openmm_customexternalforce_getglobalparametername_(const OpenMM_CustomExternalForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomExternalForce_getGlobalParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMEXTERNALFORCE_GETGLOBALPARAMETERNAME(const OpenMM_CustomExternalForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomExternalForce_getGlobalParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customexternalforce_setglobalparametername_(OpenMM_CustomExternalForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomExternalForce_setGlobalParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMEXTERNALFORCE_SETGLOBALPARAMETERNAME(OpenMM_CustomExternalForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomExternalForce_setGlobalParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT double openmm_customexternalforce_getglobalparameterdefaultvalue_(const OpenMM_CustomExternalForce*& target, int const& index) {
    return OpenMM_CustomExternalForce_getGlobalParameterDefaultValue(target, index);
}
OPENMM_EXPORT double OPENMM_CUSTOMEXTERNALFORCE_GETGLOBALPARAMETERDEFAULTVALUE(const OpenMM_CustomExternalForce*& target, int const& index) {
    return OpenMM_CustomExternalForce_getGlobalParameterDefaultValue(target, index);
}
OPENMM_EXPORT void openmm_customexternalforce_setglobalparameterdefaultvalue_(OpenMM_CustomExternalForce*& target, int const& index, double const& defaultValue) {
    OpenMM_CustomExternalForce_setGlobalParameterDefaultValue(target, index, defaultValue);
}
OPENMM_EXPORT void OPENMM_CUSTOMEXTERNALFORCE_SETGLOBALPARAMETERDEFAULTVALUE(OpenMM_CustomExternalForce*& target, int const& index, double const& defaultValue) {
    OpenMM_CustomExternalForce_setGlobalParameterDefaultValue(target, index, defaultValue);
}
OPENMM_EXPORT int openmm_customexternalforce_addparticle_(OpenMM_CustomExternalForce*& target, int const& particle, const OpenMM_DoubleArray*& parameters) {
    return OpenMM_CustomExternalForce_addParticle(target, particle, parameters);
}
OPENMM_EXPORT int OPENMM_CUSTOMEXTERNALFORCE_ADDPARTICLE(OpenMM_CustomExternalForce*& target, int const& particle, const OpenMM_DoubleArray*& parameters) {
    return OpenMM_CustomExternalForce_addParticle(target, particle, parameters);
}
OPENMM_EXPORT void openmm_customexternalforce_getparticleparameters_(const OpenMM_CustomExternalForce*& target, int const& index, int* particle, OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomExternalForce_getParticleParameters(target, index, particle, parameters);
}
OPENMM_EXPORT void OPENMM_CUSTOMEXTERNALFORCE_GETPARTICLEPARAMETERS(const OpenMM_CustomExternalForce*& target, int const& index, int* particle, OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomExternalForce_getParticleParameters(target, index, particle, parameters);
}
OPENMM_EXPORT void openmm_customexternalforce_setparticleparameters_(OpenMM_CustomExternalForce*& target, int const& index, int const& particle, const OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomExternalForce_setParticleParameters(target, index, particle, parameters);
}
OPENMM_EXPORT void OPENMM_CUSTOMEXTERNALFORCE_SETPARTICLEPARAMETERS(OpenMM_CustomExternalForce*& target, int const& index, int const& particle, const OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomExternalForce_setParticleParameters(target, index, particle, parameters);
}
OPENMM_EXPORT void openmm_customexternalforce_updateparametersincontext_(OpenMM_CustomExternalForce*& target, OpenMM_Context*& context) {
    OpenMM_CustomExternalForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void OPENMM_CUSTOMEXTERNALFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_CustomExternalForce*& target, OpenMM_Context*& context) {
    OpenMM_CustomExternalForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void openmm_customexternalforce_usesperiodicboundaryconditions_(const OpenMM_CustomExternalForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CustomExternalForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_CUSTOMEXTERNALFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_CustomExternalForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CustomExternalForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::CustomIntegrator */
OPENMM_EXPORT void openmm_customintegrator_create_(OpenMM_CustomIntegrator*& result, double const& stepSize) {
    result = OpenMM_CustomIntegrator_create(stepSize);
}
OPENMM_EXPORT void OPENMM_CUSTOMINTEGRATOR_CREATE(OpenMM_CustomIntegrator*& result, double const& stepSize) {
    result = OpenMM_CustomIntegrator_create(stepSize);
}
OPENMM_EXPORT void openmm_customintegrator_destroy_(OpenMM_CustomIntegrator*& destroy) {
    OpenMM_CustomIntegrator_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_CUSTOMINTEGRATOR_DESTROY(OpenMM_CustomIntegrator*& destroy) {
    OpenMM_CustomIntegrator_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT int openmm_customintegrator_getnumglobalvariables_(const OpenMM_CustomIntegrator*& target) {
    return OpenMM_CustomIntegrator_getNumGlobalVariables(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMINTEGRATOR_GETNUMGLOBALVARIABLES(const OpenMM_CustomIntegrator*& target) {
    return OpenMM_CustomIntegrator_getNumGlobalVariables(target);
}
OPENMM_EXPORT int openmm_customintegrator_getnumperdofvariables_(const OpenMM_CustomIntegrator*& target) {
    return OpenMM_CustomIntegrator_getNumPerDofVariables(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMINTEGRATOR_GETNUMPERDOFVARIABLES(const OpenMM_CustomIntegrator*& target) {
    return OpenMM_CustomIntegrator_getNumPerDofVariables(target);
}
OPENMM_EXPORT int openmm_customintegrator_getnumcomputations_(const OpenMM_CustomIntegrator*& target) {
    return OpenMM_CustomIntegrator_getNumComputations(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMINTEGRATOR_GETNUMCOMPUTATIONS(const OpenMM_CustomIntegrator*& target) {
    return OpenMM_CustomIntegrator_getNumComputations(target);
}
OPENMM_EXPORT int openmm_customintegrator_getnumtabulatedfunctions_(const OpenMM_CustomIntegrator*& target) {
    return OpenMM_CustomIntegrator_getNumTabulatedFunctions(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMINTEGRATOR_GETNUMTABULATEDFUNCTIONS(const OpenMM_CustomIntegrator*& target) {
    return OpenMM_CustomIntegrator_getNumTabulatedFunctions(target);
}
OPENMM_EXPORT int openmm_customintegrator_addglobalvariable_(OpenMM_CustomIntegrator*& target, const char* name, double const& initialValue, int name_length) {
    return OpenMM_CustomIntegrator_addGlobalVariable(target, makeString(name, name_length).c_str(), initialValue);
}
OPENMM_EXPORT int OPENMM_CUSTOMINTEGRATOR_ADDGLOBALVARIABLE(OpenMM_CustomIntegrator*& target, const char* name, double const& initialValue, int name_length) {
    return OpenMM_CustomIntegrator_addGlobalVariable(target, makeString(name, name_length).c_str(), initialValue);
}
OPENMM_EXPORT void openmm_customintegrator_getglobalvariablename_(const OpenMM_CustomIntegrator*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomIntegrator_getGlobalVariableName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMINTEGRATOR_GETGLOBALVARIABLENAME(const OpenMM_CustomIntegrator*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomIntegrator_getGlobalVariableName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT int openmm_customintegrator_addperdofvariable_(OpenMM_CustomIntegrator*& target, const char* name, double const& initialValue, int name_length) {
    return OpenMM_CustomIntegrator_addPerDofVariable(target, makeString(name, name_length).c_str(), initialValue);
}
OPENMM_EXPORT int OPENMM_CUSTOMINTEGRATOR_ADDPERDOFVARIABLE(OpenMM_CustomIntegrator*& target, const char* name, double const& initialValue, int name_length) {
    return OpenMM_CustomIntegrator_addPerDofVariable(target, makeString(name, name_length).c_str(), initialValue);
}
OPENMM_EXPORT void openmm_customintegrator_getperdofvariablename_(const OpenMM_CustomIntegrator*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomIntegrator_getPerDofVariableName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMINTEGRATOR_GETPERDOFVARIABLENAME(const OpenMM_CustomIntegrator*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomIntegrator_getPerDofVariableName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT double openmm_customintegrator_getglobalvariable_(const OpenMM_CustomIntegrator*& target, int const& index) {
    return OpenMM_CustomIntegrator_getGlobalVariable(target, index);
}
OPENMM_EXPORT double OPENMM_CUSTOMINTEGRATOR_GETGLOBALVARIABLE(const OpenMM_CustomIntegrator*& target, int const& index) {
    return OpenMM_CustomIntegrator_getGlobalVariable(target, index);
}
OPENMM_EXPORT double openmm_customintegrator_getglobalvariablebyname_(const OpenMM_CustomIntegrator*& target, const char* name, int name_length) {
    return OpenMM_CustomIntegrator_getGlobalVariableByName(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT double OPENMM_CUSTOMINTEGRATOR_GETGLOBALVARIABLEBYNAME(const OpenMM_CustomIntegrator*& target, const char* name, int name_length) {
    return OpenMM_CustomIntegrator_getGlobalVariableByName(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void openmm_customintegrator_setglobalvariable_(OpenMM_CustomIntegrator*& target, int const& index, double const& value) {
    OpenMM_CustomIntegrator_setGlobalVariable(target, index, value);
}
OPENMM_EXPORT void OPENMM_CUSTOMINTEGRATOR_SETGLOBALVARIABLE(OpenMM_CustomIntegrator*& target, int const& index, double const& value) {
    OpenMM_CustomIntegrator_setGlobalVariable(target, index, value);
}
OPENMM_EXPORT void openmm_customintegrator_setglobalvariablebyname_(OpenMM_CustomIntegrator*& target, const char* name, double const& value, int name_length) {
    OpenMM_CustomIntegrator_setGlobalVariableByName(target, makeString(name, name_length).c_str(), value);
}
OPENMM_EXPORT void OPENMM_CUSTOMINTEGRATOR_SETGLOBALVARIABLEBYNAME(OpenMM_CustomIntegrator*& target, const char* name, double const& value, int name_length) {
    OpenMM_CustomIntegrator_setGlobalVariableByName(target, makeString(name, name_length).c_str(), value);
}
OPENMM_EXPORT void openmm_customintegrator_getperdofvariable_(const OpenMM_CustomIntegrator*& target, int const& index, OpenMM_Vec3Array*& values) {
    OpenMM_CustomIntegrator_getPerDofVariable(target, index, values);
}
OPENMM_EXPORT void OPENMM_CUSTOMINTEGRATOR_GETPERDOFVARIABLE(const OpenMM_CustomIntegrator*& target, int const& index, OpenMM_Vec3Array*& values) {
    OpenMM_CustomIntegrator_getPerDofVariable(target, index, values);
}
OPENMM_EXPORT void openmm_customintegrator_getperdofvariablebyname_(const OpenMM_CustomIntegrator*& target, const char* name, OpenMM_Vec3Array*& values, int name_length) {
    OpenMM_CustomIntegrator_getPerDofVariableByName(target, makeString(name, name_length).c_str(), values);
}
OPENMM_EXPORT void OPENMM_CUSTOMINTEGRATOR_GETPERDOFVARIABLEBYNAME(const OpenMM_CustomIntegrator*& target, const char* name, OpenMM_Vec3Array*& values, int name_length) {
    OpenMM_CustomIntegrator_getPerDofVariableByName(target, makeString(name, name_length).c_str(), values);
}
OPENMM_EXPORT void openmm_customintegrator_setperdofvariable_(OpenMM_CustomIntegrator*& target, int const& index, const OpenMM_Vec3Array*& values) {
    OpenMM_CustomIntegrator_setPerDofVariable(target, index, values);
}
OPENMM_EXPORT void OPENMM_CUSTOMINTEGRATOR_SETPERDOFVARIABLE(OpenMM_CustomIntegrator*& target, int const& index, const OpenMM_Vec3Array*& values) {
    OpenMM_CustomIntegrator_setPerDofVariable(target, index, values);
}
OPENMM_EXPORT void openmm_customintegrator_setperdofvariablebyname_(OpenMM_CustomIntegrator*& target, const char* name, const OpenMM_Vec3Array*& values, int name_length) {
    OpenMM_CustomIntegrator_setPerDofVariableByName(target, makeString(name, name_length).c_str(), values);
}
OPENMM_EXPORT void OPENMM_CUSTOMINTEGRATOR_SETPERDOFVARIABLEBYNAME(OpenMM_CustomIntegrator*& target, const char* name, const OpenMM_Vec3Array*& values, int name_length) {
    OpenMM_CustomIntegrator_setPerDofVariableByName(target, makeString(name, name_length).c_str(), values);
}
OPENMM_EXPORT int openmm_customintegrator_addcomputeglobal_(OpenMM_CustomIntegrator*& target, const char* variable, const char* expression, int variable_length, int expression_length) {
    return OpenMM_CustomIntegrator_addComputeGlobal(target, makeString(variable, variable_length).c_str(), makeString(expression, expression_length).c_str());
}
OPENMM_EXPORT int OPENMM_CUSTOMINTEGRATOR_ADDCOMPUTEGLOBAL(OpenMM_CustomIntegrator*& target, const char* variable, const char* expression, int variable_length, int expression_length) {
    return OpenMM_CustomIntegrator_addComputeGlobal(target, makeString(variable, variable_length).c_str(), makeString(expression, expression_length).c_str());
}
OPENMM_EXPORT int openmm_customintegrator_addcomputeperdof_(OpenMM_CustomIntegrator*& target, const char* variable, const char* expression, int variable_length, int expression_length) {
    return OpenMM_CustomIntegrator_addComputePerDof(target, makeString(variable, variable_length).c_str(), makeString(expression, expression_length).c_str());
}
OPENMM_EXPORT int OPENMM_CUSTOMINTEGRATOR_ADDCOMPUTEPERDOF(OpenMM_CustomIntegrator*& target, const char* variable, const char* expression, int variable_length, int expression_length) {
    return OpenMM_CustomIntegrator_addComputePerDof(target, makeString(variable, variable_length).c_str(), makeString(expression, expression_length).c_str());
}
OPENMM_EXPORT int openmm_customintegrator_addcomputesum_(OpenMM_CustomIntegrator*& target, const char* variable, const char* expression, int variable_length, int expression_length) {
    return OpenMM_CustomIntegrator_addComputeSum(target, makeString(variable, variable_length).c_str(), makeString(expression, expression_length).c_str());
}
OPENMM_EXPORT int OPENMM_CUSTOMINTEGRATOR_ADDCOMPUTESUM(OpenMM_CustomIntegrator*& target, const char* variable, const char* expression, int variable_length, int expression_length) {
    return OpenMM_CustomIntegrator_addComputeSum(target, makeString(variable, variable_length).c_str(), makeString(expression, expression_length).c_str());
}
OPENMM_EXPORT int openmm_customintegrator_addconstrainpositions_(OpenMM_CustomIntegrator*& target) {
    return OpenMM_CustomIntegrator_addConstrainPositions(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMINTEGRATOR_ADDCONSTRAINPOSITIONS(OpenMM_CustomIntegrator*& target) {
    return OpenMM_CustomIntegrator_addConstrainPositions(target);
}
OPENMM_EXPORT int openmm_customintegrator_addconstrainvelocities_(OpenMM_CustomIntegrator*& target) {
    return OpenMM_CustomIntegrator_addConstrainVelocities(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMINTEGRATOR_ADDCONSTRAINVELOCITIES(OpenMM_CustomIntegrator*& target) {
    return OpenMM_CustomIntegrator_addConstrainVelocities(target);
}
OPENMM_EXPORT int openmm_customintegrator_addupdatecontextstate_(OpenMM_CustomIntegrator*& target) {
    return OpenMM_CustomIntegrator_addUpdateContextState(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMINTEGRATOR_ADDUPDATECONTEXTSTATE(OpenMM_CustomIntegrator*& target) {
    return OpenMM_CustomIntegrator_addUpdateContextState(target);
}
OPENMM_EXPORT int openmm_customintegrator_beginifblock_(OpenMM_CustomIntegrator*& target, const char* condition, int condition_length) {
    return OpenMM_CustomIntegrator_beginIfBlock(target, makeString(condition, condition_length).c_str());
}
OPENMM_EXPORT int OPENMM_CUSTOMINTEGRATOR_BEGINIFBLOCK(OpenMM_CustomIntegrator*& target, const char* condition, int condition_length) {
    return OpenMM_CustomIntegrator_beginIfBlock(target, makeString(condition, condition_length).c_str());
}
OPENMM_EXPORT int openmm_customintegrator_beginwhileblock_(OpenMM_CustomIntegrator*& target, const char* condition, int condition_length) {
    return OpenMM_CustomIntegrator_beginWhileBlock(target, makeString(condition, condition_length).c_str());
}
OPENMM_EXPORT int OPENMM_CUSTOMINTEGRATOR_BEGINWHILEBLOCK(OpenMM_CustomIntegrator*& target, const char* condition, int condition_length) {
    return OpenMM_CustomIntegrator_beginWhileBlock(target, makeString(condition, condition_length).c_str());
}
OPENMM_EXPORT int openmm_customintegrator_endblock_(OpenMM_CustomIntegrator*& target) {
    return OpenMM_CustomIntegrator_endBlock(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMINTEGRATOR_ENDBLOCK(OpenMM_CustomIntegrator*& target) {
    return OpenMM_CustomIntegrator_endBlock(target);
}
OPENMM_EXPORT void openmm_customintegrator_getcomputationstep_(const OpenMM_CustomIntegrator*& target, int const& index, OpenMM_CustomIntegrator_ComputationType*& type, char** variable, char** expression) {
    OpenMM_CustomIntegrator_getComputationStep(target, index, type, variable, expression);
}
OPENMM_EXPORT void OPENMM_CUSTOMINTEGRATOR_GETCOMPUTATIONSTEP(const OpenMM_CustomIntegrator*& target, int const& index, OpenMM_CustomIntegrator_ComputationType*& type, char** variable, char** expression) {
    OpenMM_CustomIntegrator_getComputationStep(target, index, type, variable, expression);
}
OPENMM_EXPORT int openmm_customintegrator_addtabulatedfunction_(OpenMM_CustomIntegrator*& target, const char* name, OpenMM_TabulatedFunction*& function, int name_length) {
    return OpenMM_CustomIntegrator_addTabulatedFunction(target, makeString(name, name_length).c_str(), function);
}
OPENMM_EXPORT int OPENMM_CUSTOMINTEGRATOR_ADDTABULATEDFUNCTION(OpenMM_CustomIntegrator*& target, const char* name, OpenMM_TabulatedFunction*& function, int name_length) {
    return OpenMM_CustomIntegrator_addTabulatedFunction(target, makeString(name, name_length).c_str(), function);
}
OPENMM_EXPORT void openmm_customintegrator_gettabulatedfunction_(OpenMM_CustomIntegrator*& target, int const& index, OpenMM_TabulatedFunction*& result) {
    result = OpenMM_CustomIntegrator_getTabulatedFunction(target, index);
}
OPENMM_EXPORT void OPENMM_CUSTOMINTEGRATOR_GETTABULATEDFUNCTION(OpenMM_CustomIntegrator*& target, int const& index, OpenMM_TabulatedFunction*& result) {
    result = OpenMM_CustomIntegrator_getTabulatedFunction(target, index);
}
OPENMM_EXPORT void openmm_customintegrator_gettabulatedfunctionname_(const OpenMM_CustomIntegrator*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomIntegrator_getTabulatedFunctionName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMINTEGRATOR_GETTABULATEDFUNCTIONNAME(const OpenMM_CustomIntegrator*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomIntegrator_getTabulatedFunctionName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customintegrator_getkineticenergyexpression_(const OpenMM_CustomIntegrator*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomIntegrator_getKineticEnergyExpression(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMINTEGRATOR_GETKINETICENERGYEXPRESSION(const OpenMM_CustomIntegrator*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomIntegrator_getKineticEnergyExpression(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customintegrator_setkineticenergyexpression_(OpenMM_CustomIntegrator*& target, const char* expression, int expression_length) {
    OpenMM_CustomIntegrator_setKineticEnergyExpression(target, makeString(expression, expression_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMINTEGRATOR_SETKINETICENERGYEXPRESSION(OpenMM_CustomIntegrator*& target, const char* expression, int expression_length) {
    OpenMM_CustomIntegrator_setKineticEnergyExpression(target, makeString(expression, expression_length).c_str());
}
OPENMM_EXPORT int openmm_customintegrator_getrandomnumberseed_(const OpenMM_CustomIntegrator*& target) {
    return OpenMM_CustomIntegrator_getRandomNumberSeed(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMINTEGRATOR_GETRANDOMNUMBERSEED(const OpenMM_CustomIntegrator*& target) {
    return OpenMM_CustomIntegrator_getRandomNumberSeed(target);
}
OPENMM_EXPORT void openmm_customintegrator_setrandomnumberseed_(OpenMM_CustomIntegrator*& target, int const& seed) {
    OpenMM_CustomIntegrator_setRandomNumberSeed(target, seed);
}
OPENMM_EXPORT void OPENMM_CUSTOMINTEGRATOR_SETRANDOMNUMBERSEED(OpenMM_CustomIntegrator*& target, int const& seed) {
    OpenMM_CustomIntegrator_setRandomNumberSeed(target, seed);
}
OPENMM_EXPORT void openmm_customintegrator_step_(OpenMM_CustomIntegrator*& target, int const& steps) {
    OpenMM_CustomIntegrator_step(target, steps);
}
OPENMM_EXPORT void OPENMM_CUSTOMINTEGRATOR_STEP(OpenMM_CustomIntegrator*& target, int const& steps) {
    OpenMM_CustomIntegrator_step(target, steps);
}

/* OpenMM::Continuous2DFunction */
OPENMM_EXPORT void openmm_continuous2dfunction_create_(OpenMM_Continuous2DFunction*& result, int const& xsize, int const& ysize, const OpenMM_DoubleArray*& values, double const& xmin, double const& xmax, double const& ymin, double const& ymax) {
    result = OpenMM_Continuous2DFunction_create(xsize, ysize, values, xmin, xmax, ymin, ymax);
}
OPENMM_EXPORT void OPENMM_CONTINUOUS2DFUNCTION_CREATE(OpenMM_Continuous2DFunction*& result, int const& xsize, int const& ysize, const OpenMM_DoubleArray*& values, double const& xmin, double const& xmax, double const& ymin, double const& ymax) {
    result = OpenMM_Continuous2DFunction_create(xsize, ysize, values, xmin, xmax, ymin, ymax);
}
OPENMM_EXPORT void openmm_continuous2dfunction_destroy_(OpenMM_Continuous2DFunction*& destroy) {
    OpenMM_Continuous2DFunction_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_CONTINUOUS2DFUNCTION_DESTROY(OpenMM_Continuous2DFunction*& destroy) {
    OpenMM_Continuous2DFunction_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void openmm_continuous2dfunction_getfunctionparameters_(const OpenMM_Continuous2DFunction*& target, int* xsize, int* ysize, OpenMM_DoubleArray*& values, double* xmin, double* xmax, double* ymin, double* ymax) {
    OpenMM_Continuous2DFunction_getFunctionParameters(target, xsize, ysize, values, xmin, xmax, ymin, ymax);
}
OPENMM_EXPORT void OPENMM_CONTINUOUS2DFUNCTION_GETFUNCTIONPARAMETERS(const OpenMM_Continuous2DFunction*& target, int* xsize, int* ysize, OpenMM_DoubleArray*& values, double* xmin, double* xmax, double* ymin, double* ymax) {
    OpenMM_Continuous2DFunction_getFunctionParameters(target, xsize, ysize, values, xmin, xmax, ymin, ymax);
}
OPENMM_EXPORT void openmm_continuous2dfunction_setfunctionparameters_(OpenMM_Continuous2DFunction*& target, int const& xsize, int const& ysize, const OpenMM_DoubleArray*& values, double const& xmin, double const& xmax, double const& ymin, double const& ymax) {
    OpenMM_Continuous2DFunction_setFunctionParameters(target, xsize, ysize, values, xmin, xmax, ymin, ymax);
}
OPENMM_EXPORT void OPENMM_CONTINUOUS2DFUNCTION_SETFUNCTIONPARAMETERS(OpenMM_Continuous2DFunction*& target, int const& xsize, int const& ysize, const OpenMM_DoubleArray*& values, double const& xmin, double const& xmax, double const& ymin, double const& ymax) {
    OpenMM_Continuous2DFunction_setFunctionParameters(target, xsize, ysize, values, xmin, xmax, ymin, ymax);
}
OPENMM_EXPORT void openmm_continuous2dfunction_copy_(const OpenMM_Continuous2DFunction*& target, OpenMM_Continuous2DFunction*& result) {
    result = OpenMM_Continuous2DFunction_Copy(target);
}
OPENMM_EXPORT void OPENMM_CONTINUOUS2DFUNCTION_COPY(const OpenMM_Continuous2DFunction*& target, OpenMM_Continuous2DFunction*& result) {
    result = OpenMM_Continuous2DFunction_Copy(target);
}

/* OpenMM::GayBerneForce */
OPENMM_EXPORT void openmm_gayberneforce_create_(OpenMM_GayBerneForce*& result) {
    result = OpenMM_GayBerneForce_create();
}
OPENMM_EXPORT void OPENMM_GAYBERNEFORCE_CREATE(OpenMM_GayBerneForce*& result) {
    result = OpenMM_GayBerneForce_create();
}
OPENMM_EXPORT void openmm_gayberneforce_destroy_(OpenMM_GayBerneForce*& destroy) {
    OpenMM_GayBerneForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_GAYBERNEFORCE_DESTROY(OpenMM_GayBerneForce*& destroy) {
    OpenMM_GayBerneForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT int openmm_gayberneforce_getnumparticles_(const OpenMM_GayBerneForce*& target) {
    return OpenMM_GayBerneForce_getNumParticles(target);
}
OPENMM_EXPORT int OPENMM_GAYBERNEFORCE_GETNUMPARTICLES(const OpenMM_GayBerneForce*& target) {
    return OpenMM_GayBerneForce_getNumParticles(target);
}
OPENMM_EXPORT int openmm_gayberneforce_getnumexceptions_(const OpenMM_GayBerneForce*& target) {
    return OpenMM_GayBerneForce_getNumExceptions(target);
}
OPENMM_EXPORT int OPENMM_GAYBERNEFORCE_GETNUMEXCEPTIONS(const OpenMM_GayBerneForce*& target) {
    return OpenMM_GayBerneForce_getNumExceptions(target);
}
OPENMM_EXPORT void openmm_gayberneforce_getnonbondedmethod_(const OpenMM_GayBerneForce*& target, OpenMM_GayBerneForce_NonbondedMethod& result) {
    result = OpenMM_GayBerneForce_getNonbondedMethod(target);
}
OPENMM_EXPORT void OPENMM_GAYBERNEFORCE_GETNONBONDEDMETHOD(const OpenMM_GayBerneForce*& target, OpenMM_GayBerneForce_NonbondedMethod& result) {
    result = OpenMM_GayBerneForce_getNonbondedMethod(target);
}
OPENMM_EXPORT void openmm_gayberneforce_setnonbondedmethod_(OpenMM_GayBerneForce*& target, OpenMM_GayBerneForce_NonbondedMethod& method) {
    OpenMM_GayBerneForce_setNonbondedMethod(target, method);
}
OPENMM_EXPORT void OPENMM_GAYBERNEFORCE_SETNONBONDEDMETHOD(OpenMM_GayBerneForce*& target, OpenMM_GayBerneForce_NonbondedMethod& method) {
    OpenMM_GayBerneForce_setNonbondedMethod(target, method);
}
OPENMM_EXPORT double openmm_gayberneforce_getcutoffdistance_(const OpenMM_GayBerneForce*& target) {
    return OpenMM_GayBerneForce_getCutoffDistance(target);
}
OPENMM_EXPORT double OPENMM_GAYBERNEFORCE_GETCUTOFFDISTANCE(const OpenMM_GayBerneForce*& target) {
    return OpenMM_GayBerneForce_getCutoffDistance(target);
}
OPENMM_EXPORT void openmm_gayberneforce_setcutoffdistance_(OpenMM_GayBerneForce*& target, double const& distance) {
    OpenMM_GayBerneForce_setCutoffDistance(target, distance);
}
OPENMM_EXPORT void OPENMM_GAYBERNEFORCE_SETCUTOFFDISTANCE(OpenMM_GayBerneForce*& target, double const& distance) {
    OpenMM_GayBerneForce_setCutoffDistance(target, distance);
}
OPENMM_EXPORT void openmm_gayberneforce_getuseswitchingfunction_(const OpenMM_GayBerneForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_GayBerneForce_getUseSwitchingFunction(target);
}
OPENMM_EXPORT void OPENMM_GAYBERNEFORCE_GETUSESWITCHINGFUNCTION(const OpenMM_GayBerneForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_GayBerneForce_getUseSwitchingFunction(target);
}
OPENMM_EXPORT void openmm_gayberneforce_setuseswitchingfunction_(OpenMM_GayBerneForce*& target, OpenMM_Boolean& use) {
    OpenMM_GayBerneForce_setUseSwitchingFunction(target, use);
}
OPENMM_EXPORT void OPENMM_GAYBERNEFORCE_SETUSESWITCHINGFUNCTION(OpenMM_GayBerneForce*& target, OpenMM_Boolean& use) {
    OpenMM_GayBerneForce_setUseSwitchingFunction(target, use);
}
OPENMM_EXPORT double openmm_gayberneforce_getswitchingdistance_(const OpenMM_GayBerneForce*& target) {
    return OpenMM_GayBerneForce_getSwitchingDistance(target);
}
OPENMM_EXPORT double OPENMM_GAYBERNEFORCE_GETSWITCHINGDISTANCE(const OpenMM_GayBerneForce*& target) {
    return OpenMM_GayBerneForce_getSwitchingDistance(target);
}
OPENMM_EXPORT void openmm_gayberneforce_setswitchingdistance_(OpenMM_GayBerneForce*& target, double const& distance) {
    OpenMM_GayBerneForce_setSwitchingDistance(target, distance);
}
OPENMM_EXPORT void OPENMM_GAYBERNEFORCE_SETSWITCHINGDISTANCE(OpenMM_GayBerneForce*& target, double const& distance) {
    OpenMM_GayBerneForce_setSwitchingDistance(target, distance);
}
OPENMM_EXPORT int openmm_gayberneforce_addparticle_(OpenMM_GayBerneForce*& target, double const& sigma, double const& epsilon, int const& xparticle, int const& yparticle, double const& sx, double const& sy, double const& sz, double const& ex, double const& ey, double const& ez) {
    return OpenMM_GayBerneForce_addParticle(target, sigma, epsilon, xparticle, yparticle, sx, sy, sz, ex, ey, ez);
}
OPENMM_EXPORT int OPENMM_GAYBERNEFORCE_ADDPARTICLE(OpenMM_GayBerneForce*& target, double const& sigma, double const& epsilon, int const& xparticle, int const& yparticle, double const& sx, double const& sy, double const& sz, double const& ex, double const& ey, double const& ez) {
    return OpenMM_GayBerneForce_addParticle(target, sigma, epsilon, xparticle, yparticle, sx, sy, sz, ex, ey, ez);
}
OPENMM_EXPORT void openmm_gayberneforce_getparticleparameters_(const OpenMM_GayBerneForce*& target, int const& index, double* sigma, double* epsilon, int* xparticle, int* yparticle, double* sx, double* sy, double* sz, double* ex, double* ey, double* ez) {
    OpenMM_GayBerneForce_getParticleParameters(target, index, sigma, epsilon, xparticle, yparticle, sx, sy, sz, ex, ey, ez);
}
OPENMM_EXPORT void OPENMM_GAYBERNEFORCE_GETPARTICLEPARAMETERS(const OpenMM_GayBerneForce*& target, int const& index, double* sigma, double* epsilon, int* xparticle, int* yparticle, double* sx, double* sy, double* sz, double* ex, double* ey, double* ez) {
    OpenMM_GayBerneForce_getParticleParameters(target, index, sigma, epsilon, xparticle, yparticle, sx, sy, sz, ex, ey, ez);
}
OPENMM_EXPORT void openmm_gayberneforce_setparticleparameters_(OpenMM_GayBerneForce*& target, int const& index, double const& sigma, double const& epsilon, int const& xparticle, int const& yparticle, double const& sx, double const& sy, double const& sz, double const& ex, double const& ey, double const& ez) {
    OpenMM_GayBerneForce_setParticleParameters(target, index, sigma, epsilon, xparticle, yparticle, sx, sy, sz, ex, ey, ez);
}
OPENMM_EXPORT void OPENMM_GAYBERNEFORCE_SETPARTICLEPARAMETERS(OpenMM_GayBerneForce*& target, int const& index, double const& sigma, double const& epsilon, int const& xparticle, int const& yparticle, double const& sx, double const& sy, double const& sz, double const& ex, double const& ey, double const& ez) {
    OpenMM_GayBerneForce_setParticleParameters(target, index, sigma, epsilon, xparticle, yparticle, sx, sy, sz, ex, ey, ez);
}
OPENMM_EXPORT int openmm_gayberneforce_addexception_(OpenMM_GayBerneForce*& target, int const& particle1, int const& particle2, double const& sigma, double const& epsilon, OpenMM_Boolean& replace) {
    return OpenMM_GayBerneForce_addException(target, particle1, particle2, sigma, epsilon, replace);
}
OPENMM_EXPORT int OPENMM_GAYBERNEFORCE_ADDEXCEPTION(OpenMM_GayBerneForce*& target, int const& particle1, int const& particle2, double const& sigma, double const& epsilon, OpenMM_Boolean& replace) {
    return OpenMM_GayBerneForce_addException(target, particle1, particle2, sigma, epsilon, replace);
}
OPENMM_EXPORT void openmm_gayberneforce_getexceptionparameters_(const OpenMM_GayBerneForce*& target, int const& index, int* particle1, int* particle2, double* sigma, double* epsilon) {
    OpenMM_GayBerneForce_getExceptionParameters(target, index, particle1, particle2, sigma, epsilon);
}
OPENMM_EXPORT void OPENMM_GAYBERNEFORCE_GETEXCEPTIONPARAMETERS(const OpenMM_GayBerneForce*& target, int const& index, int* particle1, int* particle2, double* sigma, double* epsilon) {
    OpenMM_GayBerneForce_getExceptionParameters(target, index, particle1, particle2, sigma, epsilon);
}
OPENMM_EXPORT void openmm_gayberneforce_setexceptionparameters_(OpenMM_GayBerneForce*& target, int const& index, int const& particle1, int const& particle2, double const& sigma, double const& epsilon) {
    OpenMM_GayBerneForce_setExceptionParameters(target, index, particle1, particle2, sigma, epsilon);
}
OPENMM_EXPORT void OPENMM_GAYBERNEFORCE_SETEXCEPTIONPARAMETERS(OpenMM_GayBerneForce*& target, int const& index, int const& particle1, int const& particle2, double const& sigma, double const& epsilon) {
    OpenMM_GayBerneForce_setExceptionParameters(target, index, particle1, particle2, sigma, epsilon);
}
OPENMM_EXPORT void openmm_gayberneforce_updateparametersincontext_(OpenMM_GayBerneForce*& target, OpenMM_Context*& context) {
    OpenMM_GayBerneForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void OPENMM_GAYBERNEFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_GayBerneForce*& target, OpenMM_Context*& context) {
    OpenMM_GayBerneForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void openmm_gayberneforce_usesperiodicboundaryconditions_(const OpenMM_GayBerneForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_GayBerneForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_GAYBERNEFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_GayBerneForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_GayBerneForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::NonbondedForce */
OPENMM_EXPORT void openmm_nonbondedforce_create_(OpenMM_NonbondedForce*& result) {
    result = OpenMM_NonbondedForce_create();
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_CREATE(OpenMM_NonbondedForce*& result) {
    result = OpenMM_NonbondedForce_create();
}
OPENMM_EXPORT void openmm_nonbondedforce_destroy_(OpenMM_NonbondedForce*& destroy) {
    OpenMM_NonbondedForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_DESTROY(OpenMM_NonbondedForce*& destroy) {
    OpenMM_NonbondedForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT int openmm_nonbondedforce_getnumparticles_(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getNumParticles(target);
}
OPENMM_EXPORT int OPENMM_NONBONDEDFORCE_GETNUMPARTICLES(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getNumParticles(target);
}
OPENMM_EXPORT int openmm_nonbondedforce_getnumexceptions_(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getNumExceptions(target);
}
OPENMM_EXPORT int OPENMM_NONBONDEDFORCE_GETNUMEXCEPTIONS(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getNumExceptions(target);
}
OPENMM_EXPORT int openmm_nonbondedforce_getnumglobalparameters_(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getNumGlobalParameters(target);
}
OPENMM_EXPORT int OPENMM_NONBONDEDFORCE_GETNUMGLOBALPARAMETERS(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getNumGlobalParameters(target);
}
OPENMM_EXPORT int openmm_nonbondedforce_getnumparticleparameteroffsets_(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getNumParticleParameterOffsets(target);
}
OPENMM_EXPORT int OPENMM_NONBONDEDFORCE_GETNUMPARTICLEPARAMETEROFFSETS(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getNumParticleParameterOffsets(target);
}
OPENMM_EXPORT int openmm_nonbondedforce_getnumexceptionparameteroffsets_(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getNumExceptionParameterOffsets(target);
}
OPENMM_EXPORT int OPENMM_NONBONDEDFORCE_GETNUMEXCEPTIONPARAMETEROFFSETS(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getNumExceptionParameterOffsets(target);
}
OPENMM_EXPORT void openmm_nonbondedforce_getnonbondedmethod_(const OpenMM_NonbondedForce*& target, OpenMM_NonbondedForce_NonbondedMethod& result) {
    result = OpenMM_NonbondedForce_getNonbondedMethod(target);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_GETNONBONDEDMETHOD(const OpenMM_NonbondedForce*& target, OpenMM_NonbondedForce_NonbondedMethod& result) {
    result = OpenMM_NonbondedForce_getNonbondedMethod(target);
}
OPENMM_EXPORT void openmm_nonbondedforce_setnonbondedmethod_(OpenMM_NonbondedForce*& target, OpenMM_NonbondedForce_NonbondedMethod& method) {
    OpenMM_NonbondedForce_setNonbondedMethod(target, method);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_SETNONBONDEDMETHOD(OpenMM_NonbondedForce*& target, OpenMM_NonbondedForce_NonbondedMethod& method) {
    OpenMM_NonbondedForce_setNonbondedMethod(target, method);
}
OPENMM_EXPORT double openmm_nonbondedforce_getcutoffdistance_(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getCutoffDistance(target);
}
OPENMM_EXPORT double OPENMM_NONBONDEDFORCE_GETCUTOFFDISTANCE(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getCutoffDistance(target);
}
OPENMM_EXPORT void openmm_nonbondedforce_setcutoffdistance_(OpenMM_NonbondedForce*& target, double const& distance) {
    OpenMM_NonbondedForce_setCutoffDistance(target, distance);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_SETCUTOFFDISTANCE(OpenMM_NonbondedForce*& target, double const& distance) {
    OpenMM_NonbondedForce_setCutoffDistance(target, distance);
}
OPENMM_EXPORT void openmm_nonbondedforce_getuseswitchingfunction_(const OpenMM_NonbondedForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_NonbondedForce_getUseSwitchingFunction(target);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_GETUSESWITCHINGFUNCTION(const OpenMM_NonbondedForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_NonbondedForce_getUseSwitchingFunction(target);
}
OPENMM_EXPORT void openmm_nonbondedforce_setuseswitchingfunction_(OpenMM_NonbondedForce*& target, OpenMM_Boolean& use) {
    OpenMM_NonbondedForce_setUseSwitchingFunction(target, use);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_SETUSESWITCHINGFUNCTION(OpenMM_NonbondedForce*& target, OpenMM_Boolean& use) {
    OpenMM_NonbondedForce_setUseSwitchingFunction(target, use);
}
OPENMM_EXPORT double openmm_nonbondedforce_getswitchingdistance_(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getSwitchingDistance(target);
}
OPENMM_EXPORT double OPENMM_NONBONDEDFORCE_GETSWITCHINGDISTANCE(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getSwitchingDistance(target);
}
OPENMM_EXPORT void openmm_nonbondedforce_setswitchingdistance_(OpenMM_NonbondedForce*& target, double const& distance) {
    OpenMM_NonbondedForce_setSwitchingDistance(target, distance);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_SETSWITCHINGDISTANCE(OpenMM_NonbondedForce*& target, double const& distance) {
    OpenMM_NonbondedForce_setSwitchingDistance(target, distance);
}
OPENMM_EXPORT double openmm_nonbondedforce_getreactionfielddielectric_(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getReactionFieldDielectric(target);
}
OPENMM_EXPORT double OPENMM_NONBONDEDFORCE_GETREACTIONFIELDDIELECTRIC(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getReactionFieldDielectric(target);
}
OPENMM_EXPORT void openmm_nonbondedforce_setreactionfielddielectric_(OpenMM_NonbondedForce*& target, double const& dielectric) {
    OpenMM_NonbondedForce_setReactionFieldDielectric(target, dielectric);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_SETREACTIONFIELDDIELECTRIC(OpenMM_NonbondedForce*& target, double const& dielectric) {
    OpenMM_NonbondedForce_setReactionFieldDielectric(target, dielectric);
}
OPENMM_EXPORT double openmm_nonbondedforce_getewalderrortolerance_(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getEwaldErrorTolerance(target);
}
OPENMM_EXPORT double OPENMM_NONBONDEDFORCE_GETEWALDERRORTOLERANCE(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getEwaldErrorTolerance(target);
}
OPENMM_EXPORT void openmm_nonbondedforce_setewalderrortolerance_(OpenMM_NonbondedForce*& target, double const& tol) {
    OpenMM_NonbondedForce_setEwaldErrorTolerance(target, tol);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_SETEWALDERRORTOLERANCE(OpenMM_NonbondedForce*& target, double const& tol) {
    OpenMM_NonbondedForce_setEwaldErrorTolerance(target, tol);
}
OPENMM_EXPORT void openmm_nonbondedforce_getpmeparameters_(const OpenMM_NonbondedForce*& target, double* alpha, int* nx, int* ny, int* nz) {
    OpenMM_NonbondedForce_getPMEParameters(target, alpha, nx, ny, nz);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_GETPMEPARAMETERS(const OpenMM_NonbondedForce*& target, double* alpha, int* nx, int* ny, int* nz) {
    OpenMM_NonbondedForce_getPMEParameters(target, alpha, nx, ny, nz);
}
OPENMM_EXPORT void openmm_nonbondedforce_getljpmeparameters_(const OpenMM_NonbondedForce*& target, double* alpha, int* nx, int* ny, int* nz) {
    OpenMM_NonbondedForce_getLJPMEParameters(target, alpha, nx, ny, nz);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_GETLJPMEPARAMETERS(const OpenMM_NonbondedForce*& target, double* alpha, int* nx, int* ny, int* nz) {
    OpenMM_NonbondedForce_getLJPMEParameters(target, alpha, nx, ny, nz);
}
OPENMM_EXPORT void openmm_nonbondedforce_setpmeparameters_(OpenMM_NonbondedForce*& target, double const& alpha, int const& nx, int const& ny, int const& nz) {
    OpenMM_NonbondedForce_setPMEParameters(target, alpha, nx, ny, nz);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_SETPMEPARAMETERS(OpenMM_NonbondedForce*& target, double const& alpha, int const& nx, int const& ny, int const& nz) {
    OpenMM_NonbondedForce_setPMEParameters(target, alpha, nx, ny, nz);
}
OPENMM_EXPORT void openmm_nonbondedforce_setljpmeparameters_(OpenMM_NonbondedForce*& target, double const& alpha, int const& nx, int const& ny, int const& nz) {
    OpenMM_NonbondedForce_setLJPMEParameters(target, alpha, nx, ny, nz);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_SETLJPMEPARAMETERS(OpenMM_NonbondedForce*& target, double const& alpha, int const& nx, int const& ny, int const& nz) {
    OpenMM_NonbondedForce_setLJPMEParameters(target, alpha, nx, ny, nz);
}
OPENMM_EXPORT void openmm_nonbondedforce_getpmeparametersincontext_(const OpenMM_NonbondedForce*& target, const OpenMM_Context*& context, double* alpha, int* nx, int* ny, int* nz) {
    OpenMM_NonbondedForce_getPMEParametersInContext(target, context, alpha, nx, ny, nz);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_GETPMEPARAMETERSINCONTEXT(const OpenMM_NonbondedForce*& target, const OpenMM_Context*& context, double* alpha, int* nx, int* ny, int* nz) {
    OpenMM_NonbondedForce_getPMEParametersInContext(target, context, alpha, nx, ny, nz);
}
OPENMM_EXPORT void openmm_nonbondedforce_getljpmeparametersincontext_(const OpenMM_NonbondedForce*& target, const OpenMM_Context*& context, double* alpha, int* nx, int* ny, int* nz) {
    OpenMM_NonbondedForce_getLJPMEParametersInContext(target, context, alpha, nx, ny, nz);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_GETLJPMEPARAMETERSINCONTEXT(const OpenMM_NonbondedForce*& target, const OpenMM_Context*& context, double* alpha, int* nx, int* ny, int* nz) {
    OpenMM_NonbondedForce_getLJPMEParametersInContext(target, context, alpha, nx, ny, nz);
}
OPENMM_EXPORT int openmm_nonbondedforce_addparticle_(OpenMM_NonbondedForce*& target, double const& charge, double const& sigma, double const& epsilon) {
    return OpenMM_NonbondedForce_addParticle(target, charge, sigma, epsilon);
}
OPENMM_EXPORT int OPENMM_NONBONDEDFORCE_ADDPARTICLE(OpenMM_NonbondedForce*& target, double const& charge, double const& sigma, double const& epsilon) {
    return OpenMM_NonbondedForce_addParticle(target, charge, sigma, epsilon);
}
OPENMM_EXPORT void openmm_nonbondedforce_getparticleparameters_(const OpenMM_NonbondedForce*& target, int const& index, double* charge, double* sigma, double* epsilon) {
    OpenMM_NonbondedForce_getParticleParameters(target, index, charge, sigma, epsilon);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_GETPARTICLEPARAMETERS(const OpenMM_NonbondedForce*& target, int const& index, double* charge, double* sigma, double* epsilon) {
    OpenMM_NonbondedForce_getParticleParameters(target, index, charge, sigma, epsilon);
}
OPENMM_EXPORT void openmm_nonbondedforce_setparticleparameters_(OpenMM_NonbondedForce*& target, int const& index, double const& charge, double const& sigma, double const& epsilon) {
    OpenMM_NonbondedForce_setParticleParameters(target, index, charge, sigma, epsilon);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_SETPARTICLEPARAMETERS(OpenMM_NonbondedForce*& target, int const& index, double const& charge, double const& sigma, double const& epsilon) {
    OpenMM_NonbondedForce_setParticleParameters(target, index, charge, sigma, epsilon);
}
OPENMM_EXPORT int openmm_nonbondedforce_addexception_(OpenMM_NonbondedForce*& target, int const& particle1, int const& particle2, double const& chargeProd, double const& sigma, double const& epsilon, OpenMM_Boolean& replace) {
    return OpenMM_NonbondedForce_addException(target, particle1, particle2, chargeProd, sigma, epsilon, replace);
}
OPENMM_EXPORT int OPENMM_NONBONDEDFORCE_ADDEXCEPTION(OpenMM_NonbondedForce*& target, int const& particle1, int const& particle2, double const& chargeProd, double const& sigma, double const& epsilon, OpenMM_Boolean& replace) {
    return OpenMM_NonbondedForce_addException(target, particle1, particle2, chargeProd, sigma, epsilon, replace);
}
OPENMM_EXPORT void openmm_nonbondedforce_getexceptionparameters_(const OpenMM_NonbondedForce*& target, int const& index, int* particle1, int* particle2, double* chargeProd, double* sigma, double* epsilon) {
    OpenMM_NonbondedForce_getExceptionParameters(target, index, particle1, particle2, chargeProd, sigma, epsilon);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_GETEXCEPTIONPARAMETERS(const OpenMM_NonbondedForce*& target, int const& index, int* particle1, int* particle2, double* chargeProd, double* sigma, double* epsilon) {
    OpenMM_NonbondedForce_getExceptionParameters(target, index, particle1, particle2, chargeProd, sigma, epsilon);
}
OPENMM_EXPORT void openmm_nonbondedforce_setexceptionparameters_(OpenMM_NonbondedForce*& target, int const& index, int const& particle1, int const& particle2, double const& chargeProd, double const& sigma, double const& epsilon) {
    OpenMM_NonbondedForce_setExceptionParameters(target, index, particle1, particle2, chargeProd, sigma, epsilon);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_SETEXCEPTIONPARAMETERS(OpenMM_NonbondedForce*& target, int const& index, int const& particle1, int const& particle2, double const& chargeProd, double const& sigma, double const& epsilon) {
    OpenMM_NonbondedForce_setExceptionParameters(target, index, particle1, particle2, chargeProd, sigma, epsilon);
}
OPENMM_EXPORT void openmm_nonbondedforce_createexceptionsfrombonds_(OpenMM_NonbondedForce*& target, const OpenMM_BondArray*& bonds, double const& coulomb14Scale, double const& lj14Scale) {
    OpenMM_NonbondedForce_createExceptionsFromBonds(target, bonds, coulomb14Scale, lj14Scale);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_CREATEEXCEPTIONSFROMBONDS(OpenMM_NonbondedForce*& target, const OpenMM_BondArray*& bonds, double const& coulomb14Scale, double const& lj14Scale) {
    OpenMM_NonbondedForce_createExceptionsFromBonds(target, bonds, coulomb14Scale, lj14Scale);
}
OPENMM_EXPORT int openmm_nonbondedforce_addglobalparameter_(OpenMM_NonbondedForce*& target, const char* name, double const& defaultValue, int name_length) {
    return OpenMM_NonbondedForce_addGlobalParameter(target, makeString(name, name_length).c_str(), defaultValue);
}
OPENMM_EXPORT int OPENMM_NONBONDEDFORCE_ADDGLOBALPARAMETER(OpenMM_NonbondedForce*& target, const char* name, double const& defaultValue, int name_length) {
    return OpenMM_NonbondedForce_addGlobalParameter(target, makeString(name, name_length).c_str(), defaultValue);
}
OPENMM_EXPORT void openmm_nonbondedforce_getglobalparametername_(const OpenMM_NonbondedForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_NonbondedForce_getGlobalParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_GETGLOBALPARAMETERNAME(const OpenMM_NonbondedForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_NonbondedForce_getGlobalParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_nonbondedforce_setglobalparametername_(OpenMM_NonbondedForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_NonbondedForce_setGlobalParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_SETGLOBALPARAMETERNAME(OpenMM_NonbondedForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_NonbondedForce_setGlobalParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT double openmm_nonbondedforce_getglobalparameterdefaultvalue_(const OpenMM_NonbondedForce*& target, int const& index) {
    return OpenMM_NonbondedForce_getGlobalParameterDefaultValue(target, index);
}
OPENMM_EXPORT double OPENMM_NONBONDEDFORCE_GETGLOBALPARAMETERDEFAULTVALUE(const OpenMM_NonbondedForce*& target, int const& index) {
    return OpenMM_NonbondedForce_getGlobalParameterDefaultValue(target, index);
}
OPENMM_EXPORT void openmm_nonbondedforce_setglobalparameterdefaultvalue_(OpenMM_NonbondedForce*& target, int const& index, double const& defaultValue) {
    OpenMM_NonbondedForce_setGlobalParameterDefaultValue(target, index, defaultValue);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_SETGLOBALPARAMETERDEFAULTVALUE(OpenMM_NonbondedForce*& target, int const& index, double const& defaultValue) {
    OpenMM_NonbondedForce_setGlobalParameterDefaultValue(target, index, defaultValue);
}
OPENMM_EXPORT int openmm_nonbondedforce_addparticleparameteroffset_(OpenMM_NonbondedForce*& target, const char* parameter, int const& particleIndex, double const& chargeScale, double const& sigmaScale, double const& epsilonScale, int parameter_length) {
    return OpenMM_NonbondedForce_addParticleParameterOffset(target, makeString(parameter, parameter_length).c_str(), particleIndex, chargeScale, sigmaScale, epsilonScale);
}
OPENMM_EXPORT int OPENMM_NONBONDEDFORCE_ADDPARTICLEPARAMETEROFFSET(OpenMM_NonbondedForce*& target, const char* parameter, int const& particleIndex, double const& chargeScale, double const& sigmaScale, double const& epsilonScale, int parameter_length) {
    return OpenMM_NonbondedForce_addParticleParameterOffset(target, makeString(parameter, parameter_length).c_str(), particleIndex, chargeScale, sigmaScale, epsilonScale);
}
OPENMM_EXPORT void openmm_nonbondedforce_getparticleparameteroffset_(const OpenMM_NonbondedForce*& target, int const& index, char** parameter, int* particleIndex, double* chargeScale, double* sigmaScale, double* epsilonScale) {
    OpenMM_NonbondedForce_getParticleParameterOffset(target, index, parameter, particleIndex, chargeScale, sigmaScale, epsilonScale);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_GETPARTICLEPARAMETEROFFSET(const OpenMM_NonbondedForce*& target, int const& index, char** parameter, int* particleIndex, double* chargeScale, double* sigmaScale, double* epsilonScale) {
    OpenMM_NonbondedForce_getParticleParameterOffset(target, index, parameter, particleIndex, chargeScale, sigmaScale, epsilonScale);
}
OPENMM_EXPORT void openmm_nonbondedforce_setparticleparameteroffset_(OpenMM_NonbondedForce*& target, int const& index, const char* parameter, int const& particleIndex, double const& chargeScale, double const& sigmaScale, double const& epsilonScale, int parameter_length) {
    OpenMM_NonbondedForce_setParticleParameterOffset(target, index, makeString(parameter, parameter_length).c_str(), particleIndex, chargeScale, sigmaScale, epsilonScale);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_SETPARTICLEPARAMETEROFFSET(OpenMM_NonbondedForce*& target, int const& index, const char* parameter, int const& particleIndex, double const& chargeScale, double const& sigmaScale, double const& epsilonScale, int parameter_length) {
    OpenMM_NonbondedForce_setParticleParameterOffset(target, index, makeString(parameter, parameter_length).c_str(), particleIndex, chargeScale, sigmaScale, epsilonScale);
}
OPENMM_EXPORT int openmm_nonbondedforce_addexceptionparameteroffset_(OpenMM_NonbondedForce*& target, const char* parameter, int const& exceptionIndex, double const& chargeProdScale, double const& sigmaScale, double const& epsilonScale, int parameter_length) {
    return OpenMM_NonbondedForce_addExceptionParameterOffset(target, makeString(parameter, parameter_length).c_str(), exceptionIndex, chargeProdScale, sigmaScale, epsilonScale);
}
OPENMM_EXPORT int OPENMM_NONBONDEDFORCE_ADDEXCEPTIONPARAMETEROFFSET(OpenMM_NonbondedForce*& target, const char* parameter, int const& exceptionIndex, double const& chargeProdScale, double const& sigmaScale, double const& epsilonScale, int parameter_length) {
    return OpenMM_NonbondedForce_addExceptionParameterOffset(target, makeString(parameter, parameter_length).c_str(), exceptionIndex, chargeProdScale, sigmaScale, epsilonScale);
}
OPENMM_EXPORT void openmm_nonbondedforce_getexceptionparameteroffset_(const OpenMM_NonbondedForce*& target, int const& index, char** parameter, int* exceptionIndex, double* chargeProdScale, double* sigmaScale, double* epsilonScale) {
    OpenMM_NonbondedForce_getExceptionParameterOffset(target, index, parameter, exceptionIndex, chargeProdScale, sigmaScale, epsilonScale);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_GETEXCEPTIONPARAMETEROFFSET(const OpenMM_NonbondedForce*& target, int const& index, char** parameter, int* exceptionIndex, double* chargeProdScale, double* sigmaScale, double* epsilonScale) {
    OpenMM_NonbondedForce_getExceptionParameterOffset(target, index, parameter, exceptionIndex, chargeProdScale, sigmaScale, epsilonScale);
}
OPENMM_EXPORT void openmm_nonbondedforce_setexceptionparameteroffset_(OpenMM_NonbondedForce*& target, int const& index, const char* parameter, int const& exceptionIndex, double const& chargeProdScale, double const& sigmaScale, double const& epsilonScale, int parameter_length) {
    OpenMM_NonbondedForce_setExceptionParameterOffset(target, index, makeString(parameter, parameter_length).c_str(), exceptionIndex, chargeProdScale, sigmaScale, epsilonScale);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_SETEXCEPTIONPARAMETEROFFSET(OpenMM_NonbondedForce*& target, int const& index, const char* parameter, int const& exceptionIndex, double const& chargeProdScale, double const& sigmaScale, double const& epsilonScale, int parameter_length) {
    OpenMM_NonbondedForce_setExceptionParameterOffset(target, index, makeString(parameter, parameter_length).c_str(), exceptionIndex, chargeProdScale, sigmaScale, epsilonScale);
}
OPENMM_EXPORT void openmm_nonbondedforce_getusedispersioncorrection_(const OpenMM_NonbondedForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_NonbondedForce_getUseDispersionCorrection(target);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_GETUSEDISPERSIONCORRECTION(const OpenMM_NonbondedForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_NonbondedForce_getUseDispersionCorrection(target);
}
OPENMM_EXPORT void openmm_nonbondedforce_setusedispersioncorrection_(OpenMM_NonbondedForce*& target, OpenMM_Boolean& useCorrection) {
    OpenMM_NonbondedForce_setUseDispersionCorrection(target, useCorrection);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_SETUSEDISPERSIONCORRECTION(OpenMM_NonbondedForce*& target, OpenMM_Boolean& useCorrection) {
    OpenMM_NonbondedForce_setUseDispersionCorrection(target, useCorrection);
}
OPENMM_EXPORT int openmm_nonbondedforce_getreciprocalspaceforcegroup_(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getReciprocalSpaceForceGroup(target);
}
OPENMM_EXPORT int OPENMM_NONBONDEDFORCE_GETRECIPROCALSPACEFORCEGROUP(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getReciprocalSpaceForceGroup(target);
}
OPENMM_EXPORT void openmm_nonbondedforce_setreciprocalspaceforcegroup_(OpenMM_NonbondedForce*& target, int const& group) {
    OpenMM_NonbondedForce_setReciprocalSpaceForceGroup(target, group);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_SETRECIPROCALSPACEFORCEGROUP(OpenMM_NonbondedForce*& target, int const& group) {
    OpenMM_NonbondedForce_setReciprocalSpaceForceGroup(target, group);
}
OPENMM_EXPORT void openmm_nonbondedforce_updateparametersincontext_(OpenMM_NonbondedForce*& target, OpenMM_Context*& context) {
    OpenMM_NonbondedForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_NonbondedForce*& target, OpenMM_Context*& context) {
    OpenMM_NonbondedForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void openmm_nonbondedforce_usesperiodicboundaryconditions_(const OpenMM_NonbondedForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_NonbondedForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_NonbondedForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_NonbondedForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void openmm_nonbondedforce_getexceptionsuseperiodicboundarycondition_(const OpenMM_NonbondedForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_NonbondedForce_getExceptionsUsePeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_GETEXCEPTIONSUSEPERIODICBOUNDARYCONDITION(const OpenMM_NonbondedForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_NonbondedForce_getExceptionsUsePeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void openmm_nonbondedforce_setexceptionsuseperiodicboundarycondition_(OpenMM_NonbondedForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_NonbondedForce_setExceptionsUsePeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT void OPENMM_NONBONDEDFORCE_SETEXCEPTIONSUSEPERIODICBOUNDARYCONDITION(OpenMM_NonbondedForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_NonbondedForce_setExceptionsUsePeriodicBoundaryConditions(target, periodic);
}

/* OpenMM::NoseHooverChain */
OPENMM_EXPORT void openmm_nosehooverchain_create_(OpenMM_NoseHooverChain*& result, double const& temperature, double const& relativeTemperature, double const& collisionFrequency, double const& relativeCollisionFrequency, int const& numDOFs, int const& chainLength, int const& numMTS, int const& numYoshidaSuzuki, int const& chainID, const OpenMM_IntArray*& thermostatedAtoms, const OpenMM_BondArray*& thermostatedPairs) {
    result = OpenMM_NoseHooverChain_create(temperature, relativeTemperature, collisionFrequency, relativeCollisionFrequency, numDOFs, chainLength, numMTS, numYoshidaSuzuki, chainID, thermostatedAtoms, thermostatedPairs);
}
OPENMM_EXPORT void OPENMM_NOSEHOOVERCHAIN_CREATE(OpenMM_NoseHooverChain*& result, double const& temperature, double const& relativeTemperature, double const& collisionFrequency, double const& relativeCollisionFrequency, int const& numDOFs, int const& chainLength, int const& numMTS, int const& numYoshidaSuzuki, int const& chainID, const OpenMM_IntArray*& thermostatedAtoms, const OpenMM_BondArray*& thermostatedPairs) {
    result = OpenMM_NoseHooverChain_create(temperature, relativeTemperature, collisionFrequency, relativeCollisionFrequency, numDOFs, chainLength, numMTS, numYoshidaSuzuki, chainID, thermostatedAtoms, thermostatedPairs);
}
OPENMM_EXPORT void openmm_nosehooverchain_destroy_(OpenMM_NoseHooverChain*& destroy) {
    OpenMM_NoseHooverChain_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_NOSEHOOVERCHAIN_DESTROY(OpenMM_NoseHooverChain*& destroy) {
    OpenMM_NoseHooverChain_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT double openmm_nosehooverchain_gettemperature_(const OpenMM_NoseHooverChain*& target) {
    return OpenMM_NoseHooverChain_getTemperature(target);
}
OPENMM_EXPORT double OPENMM_NOSEHOOVERCHAIN_GETTEMPERATURE(const OpenMM_NoseHooverChain*& target) {
    return OpenMM_NoseHooverChain_getTemperature(target);
}
OPENMM_EXPORT void openmm_nosehooverchain_settemperature_(OpenMM_NoseHooverChain*& target, double const& temperature) {
    OpenMM_NoseHooverChain_setTemperature(target, temperature);
}
OPENMM_EXPORT void OPENMM_NOSEHOOVERCHAIN_SETTEMPERATURE(OpenMM_NoseHooverChain*& target, double const& temperature) {
    OpenMM_NoseHooverChain_setTemperature(target, temperature);
}
OPENMM_EXPORT double openmm_nosehooverchain_getrelativetemperature_(const OpenMM_NoseHooverChain*& target) {
    return OpenMM_NoseHooverChain_getRelativeTemperature(target);
}
OPENMM_EXPORT double OPENMM_NOSEHOOVERCHAIN_GETRELATIVETEMPERATURE(const OpenMM_NoseHooverChain*& target) {
    return OpenMM_NoseHooverChain_getRelativeTemperature(target);
}
OPENMM_EXPORT void openmm_nosehooverchain_setrelativetemperature_(OpenMM_NoseHooverChain*& target, double const& temperature) {
    OpenMM_NoseHooverChain_setRelativeTemperature(target, temperature);
}
OPENMM_EXPORT void OPENMM_NOSEHOOVERCHAIN_SETRELATIVETEMPERATURE(OpenMM_NoseHooverChain*& target, double const& temperature) {
    OpenMM_NoseHooverChain_setRelativeTemperature(target, temperature);
}
OPENMM_EXPORT double openmm_nosehooverchain_getcollisionfrequency_(const OpenMM_NoseHooverChain*& target) {
    return OpenMM_NoseHooverChain_getCollisionFrequency(target);
}
OPENMM_EXPORT double OPENMM_NOSEHOOVERCHAIN_GETCOLLISIONFREQUENCY(const OpenMM_NoseHooverChain*& target) {
    return OpenMM_NoseHooverChain_getCollisionFrequency(target);
}
OPENMM_EXPORT void openmm_nosehooverchain_setcollisionfrequency_(OpenMM_NoseHooverChain*& target, double const& frequency) {
    OpenMM_NoseHooverChain_setCollisionFrequency(target, frequency);
}
OPENMM_EXPORT void OPENMM_NOSEHOOVERCHAIN_SETCOLLISIONFREQUENCY(OpenMM_NoseHooverChain*& target, double const& frequency) {
    OpenMM_NoseHooverChain_setCollisionFrequency(target, frequency);
}
OPENMM_EXPORT double openmm_nosehooverchain_getrelativecollisionfrequency_(const OpenMM_NoseHooverChain*& target) {
    return OpenMM_NoseHooverChain_getRelativeCollisionFrequency(target);
}
OPENMM_EXPORT double OPENMM_NOSEHOOVERCHAIN_GETRELATIVECOLLISIONFREQUENCY(const OpenMM_NoseHooverChain*& target) {
    return OpenMM_NoseHooverChain_getRelativeCollisionFrequency(target);
}
OPENMM_EXPORT void openmm_nosehooverchain_setrelativecollisionfrequency_(OpenMM_NoseHooverChain*& target, double const& frequency) {
    OpenMM_NoseHooverChain_setRelativeCollisionFrequency(target, frequency);
}
OPENMM_EXPORT void OPENMM_NOSEHOOVERCHAIN_SETRELATIVECOLLISIONFREQUENCY(OpenMM_NoseHooverChain*& target, double const& frequency) {
    OpenMM_NoseHooverChain_setRelativeCollisionFrequency(target, frequency);
}
OPENMM_EXPORT int openmm_nosehooverchain_getnumdegreesoffreedom_(const OpenMM_NoseHooverChain*& target) {
    return OpenMM_NoseHooverChain_getNumDegreesOfFreedom(target);
}
OPENMM_EXPORT int OPENMM_NOSEHOOVERCHAIN_GETNUMDEGREESOFFREEDOM(const OpenMM_NoseHooverChain*& target) {
    return OpenMM_NoseHooverChain_getNumDegreesOfFreedom(target);
}
OPENMM_EXPORT void openmm_nosehooverchain_setnumdegreesoffreedom_(OpenMM_NoseHooverChain*& target, int const& numDOF) {
    OpenMM_NoseHooverChain_setNumDegreesOfFreedom(target, numDOF);
}
OPENMM_EXPORT void OPENMM_NOSEHOOVERCHAIN_SETNUMDEGREESOFFREEDOM(OpenMM_NoseHooverChain*& target, int const& numDOF) {
    OpenMM_NoseHooverChain_setNumDegreesOfFreedom(target, numDOF);
}
OPENMM_EXPORT int openmm_nosehooverchain_getchainlength_(const OpenMM_NoseHooverChain*& target) {
    return OpenMM_NoseHooverChain_getChainLength(target);
}
OPENMM_EXPORT int OPENMM_NOSEHOOVERCHAIN_GETCHAINLENGTH(const OpenMM_NoseHooverChain*& target) {
    return OpenMM_NoseHooverChain_getChainLength(target);
}
OPENMM_EXPORT int openmm_nosehooverchain_getnummultitimesteps_(const OpenMM_NoseHooverChain*& target) {
    return OpenMM_NoseHooverChain_getNumMultiTimeSteps(target);
}
OPENMM_EXPORT int OPENMM_NOSEHOOVERCHAIN_GETNUMMULTITIMESTEPS(const OpenMM_NoseHooverChain*& target) {
    return OpenMM_NoseHooverChain_getNumMultiTimeSteps(target);
}
OPENMM_EXPORT int openmm_nosehooverchain_getnumyoshidasuzukitimesteps_(const OpenMM_NoseHooverChain*& target) {
    return OpenMM_NoseHooverChain_getNumYoshidaSuzukiTimeSteps(target);
}
OPENMM_EXPORT int OPENMM_NOSEHOOVERCHAIN_GETNUMYOSHIDASUZUKITIMESTEPS(const OpenMM_NoseHooverChain*& target) {
    return OpenMM_NoseHooverChain_getNumYoshidaSuzukiTimeSteps(target);
}
OPENMM_EXPORT int openmm_nosehooverchain_getchainid_(const OpenMM_NoseHooverChain*& target) {
    return OpenMM_NoseHooverChain_getChainID(target);
}
OPENMM_EXPORT int OPENMM_NOSEHOOVERCHAIN_GETCHAINID(const OpenMM_NoseHooverChain*& target) {
    return OpenMM_NoseHooverChain_getChainID(target);
}
OPENMM_EXPORT void openmm_nosehooverchain_getthermostatedatoms_(const OpenMM_NoseHooverChain*& target, const OpenMM_IntArray*& result) {
    result = OpenMM_NoseHooverChain_getThermostatedAtoms(target);
}
OPENMM_EXPORT void OPENMM_NOSEHOOVERCHAIN_GETTHERMOSTATEDATOMS(const OpenMM_NoseHooverChain*& target, const OpenMM_IntArray*& result) {
    result = OpenMM_NoseHooverChain_getThermostatedAtoms(target);
}
OPENMM_EXPORT void openmm_nosehooverchain_setthermostatedatoms_(OpenMM_NoseHooverChain*& target, const OpenMM_IntArray*& atomIDs) {
    OpenMM_NoseHooverChain_setThermostatedAtoms(target, atomIDs);
}
OPENMM_EXPORT void OPENMM_NOSEHOOVERCHAIN_SETTHERMOSTATEDATOMS(OpenMM_NoseHooverChain*& target, const OpenMM_IntArray*& atomIDs) {
    OpenMM_NoseHooverChain_setThermostatedAtoms(target, atomIDs);
}
OPENMM_EXPORT void openmm_nosehooverchain_getthermostatedpairs_(const OpenMM_NoseHooverChain*& target, const OpenMM_BondArray*& result) {
    result = OpenMM_NoseHooverChain_getThermostatedPairs(target);
}
OPENMM_EXPORT void OPENMM_NOSEHOOVERCHAIN_GETTHERMOSTATEDPAIRS(const OpenMM_NoseHooverChain*& target, const OpenMM_BondArray*& result) {
    result = OpenMM_NoseHooverChain_getThermostatedPairs(target);
}
OPENMM_EXPORT void openmm_nosehooverchain_setthermostatedpairs_(OpenMM_NoseHooverChain*& target, const OpenMM_BondArray*& pairIDs) {
    OpenMM_NoseHooverChain_setThermostatedPairs(target, pairIDs);
}
OPENMM_EXPORT void OPENMM_NOSEHOOVERCHAIN_SETTHERMOSTATEDPAIRS(OpenMM_NoseHooverChain*& target, const OpenMM_BondArray*& pairIDs) {
    OpenMM_NoseHooverChain_setThermostatedPairs(target, pairIDs);
}
OPENMM_EXPORT void openmm_nosehooverchain_usesperiodicboundaryconditions_(const OpenMM_NoseHooverChain*& target, OpenMM_Boolean& result) {
    result = OpenMM_NoseHooverChain_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_NOSEHOOVERCHAIN_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_NoseHooverChain*& target, OpenMM_Boolean& result) {
    result = OpenMM_NoseHooverChain_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::HarmonicBondForce */
OPENMM_EXPORT void openmm_harmonicbondforce_create_(OpenMM_HarmonicBondForce*& result) {
    result = OpenMM_HarmonicBondForce_create();
}
OPENMM_EXPORT void OPENMM_HARMONICBONDFORCE_CREATE(OpenMM_HarmonicBondForce*& result) {
    result = OpenMM_HarmonicBondForce_create();
}
OPENMM_EXPORT void openmm_harmonicbondforce_destroy_(OpenMM_HarmonicBondForce*& destroy) {
    OpenMM_HarmonicBondForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_HARMONICBONDFORCE_DESTROY(OpenMM_HarmonicBondForce*& destroy) {
    OpenMM_HarmonicBondForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT int openmm_harmonicbondforce_getnumbonds_(const OpenMM_HarmonicBondForce*& target) {
    return OpenMM_HarmonicBondForce_getNumBonds(target);
}
OPENMM_EXPORT int OPENMM_HARMONICBONDFORCE_GETNUMBONDS(const OpenMM_HarmonicBondForce*& target) {
    return OpenMM_HarmonicBondForce_getNumBonds(target);
}
OPENMM_EXPORT int openmm_harmonicbondforce_addbond_(OpenMM_HarmonicBondForce*& target, int const& particle1, int const& particle2, double const& length, double const& k) {
    return OpenMM_HarmonicBondForce_addBond(target, particle1, particle2, length, k);
}
OPENMM_EXPORT int OPENMM_HARMONICBONDFORCE_ADDBOND(OpenMM_HarmonicBondForce*& target, int const& particle1, int const& particle2, double const& length, double const& k) {
    return OpenMM_HarmonicBondForce_addBond(target, particle1, particle2, length, k);
}
OPENMM_EXPORT void openmm_harmonicbondforce_getbondparameters_(const OpenMM_HarmonicBondForce*& target, int const& index, int* particle1, int* particle2, double* length, double* k) {
    OpenMM_HarmonicBondForce_getBondParameters(target, index, particle1, particle2, length, k);
}
OPENMM_EXPORT void OPENMM_HARMONICBONDFORCE_GETBONDPARAMETERS(const OpenMM_HarmonicBondForce*& target, int const& index, int* particle1, int* particle2, double* length, double* k) {
    OpenMM_HarmonicBondForce_getBondParameters(target, index, particle1, particle2, length, k);
}
OPENMM_EXPORT void openmm_harmonicbondforce_setbondparameters_(OpenMM_HarmonicBondForce*& target, int const& index, int const& particle1, int const& particle2, double const& length, double const& k) {
    OpenMM_HarmonicBondForce_setBondParameters(target, index, particle1, particle2, length, k);
}
OPENMM_EXPORT void OPENMM_HARMONICBONDFORCE_SETBONDPARAMETERS(OpenMM_HarmonicBondForce*& target, int const& index, int const& particle1, int const& particle2, double const& length, double const& k) {
    OpenMM_HarmonicBondForce_setBondParameters(target, index, particle1, particle2, length, k);
}
OPENMM_EXPORT void openmm_harmonicbondforce_updateparametersincontext_(OpenMM_HarmonicBondForce*& target, OpenMM_Context*& context) {
    OpenMM_HarmonicBondForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void OPENMM_HARMONICBONDFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_HarmonicBondForce*& target, OpenMM_Context*& context) {
    OpenMM_HarmonicBondForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void openmm_harmonicbondforce_setusesperiodicboundaryconditions_(OpenMM_HarmonicBondForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_HarmonicBondForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT void OPENMM_HARMONICBONDFORCE_SETUSESPERIODICBOUNDARYCONDITIONS(OpenMM_HarmonicBondForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_HarmonicBondForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT void openmm_harmonicbondforce_usesperiodicboundaryconditions_(const OpenMM_HarmonicBondForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_HarmonicBondForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_HARMONICBONDFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_HarmonicBondForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_HarmonicBondForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::CustomNonbondedForce */
OPENMM_EXPORT void openmm_customnonbondedforce_create_(OpenMM_CustomNonbondedForce*& result, const char* energy, int energy_length) {
    result = OpenMM_CustomNonbondedForce_create(makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_CREATE(OpenMM_CustomNonbondedForce*& result, const char* energy, int energy_length) {
    result = OpenMM_CustomNonbondedForce_create(makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void openmm_customnonbondedforce_create_2_(OpenMM_CustomNonbondedForce*& result, const OpenMM_CustomNonbondedForce*& rhs) {
    result = OpenMM_CustomNonbondedForce_create_2(rhs);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_CREATE_2(OpenMM_CustomNonbondedForce*& result, const OpenMM_CustomNonbondedForce*& rhs) {
    result = OpenMM_CustomNonbondedForce_create_2(rhs);
}
OPENMM_EXPORT void openmm_customnonbondedforce_destroy_(OpenMM_CustomNonbondedForce*& destroy) {
    OpenMM_CustomNonbondedForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_DESTROY(OpenMM_CustomNonbondedForce*& destroy) {
    OpenMM_CustomNonbondedForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT int openmm_customnonbondedforce_getnumparticles_(const OpenMM_CustomNonbondedForce*& target) {
    return OpenMM_CustomNonbondedForce_getNumParticles(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMNONBONDEDFORCE_GETNUMPARTICLES(const OpenMM_CustomNonbondedForce*& target) {
    return OpenMM_CustomNonbondedForce_getNumParticles(target);
}
OPENMM_EXPORT int openmm_customnonbondedforce_getnumexclusions_(const OpenMM_CustomNonbondedForce*& target) {
    return OpenMM_CustomNonbondedForce_getNumExclusions(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMNONBONDEDFORCE_GETNUMEXCLUSIONS(const OpenMM_CustomNonbondedForce*& target) {
    return OpenMM_CustomNonbondedForce_getNumExclusions(target);
}
OPENMM_EXPORT int openmm_customnonbondedforce_getnumperparticleparameters_(const OpenMM_CustomNonbondedForce*& target) {
    return OpenMM_CustomNonbondedForce_getNumPerParticleParameters(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMNONBONDEDFORCE_GETNUMPERPARTICLEPARAMETERS(const OpenMM_CustomNonbondedForce*& target) {
    return OpenMM_CustomNonbondedForce_getNumPerParticleParameters(target);
}
OPENMM_EXPORT int openmm_customnonbondedforce_getnumglobalparameters_(const OpenMM_CustomNonbondedForce*& target) {
    return OpenMM_CustomNonbondedForce_getNumGlobalParameters(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMNONBONDEDFORCE_GETNUMGLOBALPARAMETERS(const OpenMM_CustomNonbondedForce*& target) {
    return OpenMM_CustomNonbondedForce_getNumGlobalParameters(target);
}
OPENMM_EXPORT int openmm_customnonbondedforce_getnumtabulatedfunctions_(const OpenMM_CustomNonbondedForce*& target) {
    return OpenMM_CustomNonbondedForce_getNumTabulatedFunctions(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMNONBONDEDFORCE_GETNUMTABULATEDFUNCTIONS(const OpenMM_CustomNonbondedForce*& target) {
    return OpenMM_CustomNonbondedForce_getNumTabulatedFunctions(target);
}
OPENMM_EXPORT int openmm_customnonbondedforce_getnumfunctions_(const OpenMM_CustomNonbondedForce*& target) {
    return OpenMM_CustomNonbondedForce_getNumFunctions(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMNONBONDEDFORCE_GETNUMFUNCTIONS(const OpenMM_CustomNonbondedForce*& target) {
    return OpenMM_CustomNonbondedForce_getNumFunctions(target);
}
OPENMM_EXPORT int openmm_customnonbondedforce_getnuminteractiongroups_(const OpenMM_CustomNonbondedForce*& target) {
    return OpenMM_CustomNonbondedForce_getNumInteractionGroups(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMNONBONDEDFORCE_GETNUMINTERACTIONGROUPS(const OpenMM_CustomNonbondedForce*& target) {
    return OpenMM_CustomNonbondedForce_getNumInteractionGroups(target);
}
OPENMM_EXPORT int openmm_customnonbondedforce_getnumenergyparameterderivatives_(const OpenMM_CustomNonbondedForce*& target) {
    return OpenMM_CustomNonbondedForce_getNumEnergyParameterDerivatives(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMNONBONDEDFORCE_GETNUMENERGYPARAMETERDERIVATIVES(const OpenMM_CustomNonbondedForce*& target) {
    return OpenMM_CustomNonbondedForce_getNumEnergyParameterDerivatives(target);
}
OPENMM_EXPORT void openmm_customnonbondedforce_getenergyfunction_(const OpenMM_CustomNonbondedForce*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomNonbondedForce_getEnergyFunction(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_GETENERGYFUNCTION(const OpenMM_CustomNonbondedForce*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomNonbondedForce_getEnergyFunction(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customnonbondedforce_setenergyfunction_(OpenMM_CustomNonbondedForce*& target, const char* energy, int energy_length) {
    OpenMM_CustomNonbondedForce_setEnergyFunction(target, makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_SETENERGYFUNCTION(OpenMM_CustomNonbondedForce*& target, const char* energy, int energy_length) {
    OpenMM_CustomNonbondedForce_setEnergyFunction(target, makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void openmm_customnonbondedforce_getnonbondedmethod_(const OpenMM_CustomNonbondedForce*& target, OpenMM_CustomNonbondedForce_NonbondedMethod& result) {
    result = OpenMM_CustomNonbondedForce_getNonbondedMethod(target);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_GETNONBONDEDMETHOD(const OpenMM_CustomNonbondedForce*& target, OpenMM_CustomNonbondedForce_NonbondedMethod& result) {
    result = OpenMM_CustomNonbondedForce_getNonbondedMethod(target);
}
OPENMM_EXPORT void openmm_customnonbondedforce_setnonbondedmethod_(OpenMM_CustomNonbondedForce*& target, OpenMM_CustomNonbondedForce_NonbondedMethod& method) {
    OpenMM_CustomNonbondedForce_setNonbondedMethod(target, method);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_SETNONBONDEDMETHOD(OpenMM_CustomNonbondedForce*& target, OpenMM_CustomNonbondedForce_NonbondedMethod& method) {
    OpenMM_CustomNonbondedForce_setNonbondedMethod(target, method);
}
OPENMM_EXPORT double openmm_customnonbondedforce_getcutoffdistance_(const OpenMM_CustomNonbondedForce*& target) {
    return OpenMM_CustomNonbondedForce_getCutoffDistance(target);
}
OPENMM_EXPORT double OPENMM_CUSTOMNONBONDEDFORCE_GETCUTOFFDISTANCE(const OpenMM_CustomNonbondedForce*& target) {
    return OpenMM_CustomNonbondedForce_getCutoffDistance(target);
}
OPENMM_EXPORT void openmm_customnonbondedforce_setcutoffdistance_(OpenMM_CustomNonbondedForce*& target, double const& distance) {
    OpenMM_CustomNonbondedForce_setCutoffDistance(target, distance);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_SETCUTOFFDISTANCE(OpenMM_CustomNonbondedForce*& target, double const& distance) {
    OpenMM_CustomNonbondedForce_setCutoffDistance(target, distance);
}
OPENMM_EXPORT void openmm_customnonbondedforce_getuseswitchingfunction_(const OpenMM_CustomNonbondedForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CustomNonbondedForce_getUseSwitchingFunction(target);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_GETUSESWITCHINGFUNCTION(const OpenMM_CustomNonbondedForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CustomNonbondedForce_getUseSwitchingFunction(target);
}
OPENMM_EXPORT void openmm_customnonbondedforce_setuseswitchingfunction_(OpenMM_CustomNonbondedForce*& target, OpenMM_Boolean& use) {
    OpenMM_CustomNonbondedForce_setUseSwitchingFunction(target, use);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_SETUSESWITCHINGFUNCTION(OpenMM_CustomNonbondedForce*& target, OpenMM_Boolean& use) {
    OpenMM_CustomNonbondedForce_setUseSwitchingFunction(target, use);
}
OPENMM_EXPORT double openmm_customnonbondedforce_getswitchingdistance_(const OpenMM_CustomNonbondedForce*& target) {
    return OpenMM_CustomNonbondedForce_getSwitchingDistance(target);
}
OPENMM_EXPORT double OPENMM_CUSTOMNONBONDEDFORCE_GETSWITCHINGDISTANCE(const OpenMM_CustomNonbondedForce*& target) {
    return OpenMM_CustomNonbondedForce_getSwitchingDistance(target);
}
OPENMM_EXPORT void openmm_customnonbondedforce_setswitchingdistance_(OpenMM_CustomNonbondedForce*& target, double const& distance) {
    OpenMM_CustomNonbondedForce_setSwitchingDistance(target, distance);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_SETSWITCHINGDISTANCE(OpenMM_CustomNonbondedForce*& target, double const& distance) {
    OpenMM_CustomNonbondedForce_setSwitchingDistance(target, distance);
}
OPENMM_EXPORT void openmm_customnonbondedforce_getuselongrangecorrection_(const OpenMM_CustomNonbondedForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CustomNonbondedForce_getUseLongRangeCorrection(target);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_GETUSELONGRANGECORRECTION(const OpenMM_CustomNonbondedForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CustomNonbondedForce_getUseLongRangeCorrection(target);
}
OPENMM_EXPORT void openmm_customnonbondedforce_setuselongrangecorrection_(OpenMM_CustomNonbondedForce*& target, OpenMM_Boolean& use) {
    OpenMM_CustomNonbondedForce_setUseLongRangeCorrection(target, use);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_SETUSELONGRANGECORRECTION(OpenMM_CustomNonbondedForce*& target, OpenMM_Boolean& use) {
    OpenMM_CustomNonbondedForce_setUseLongRangeCorrection(target, use);
}
OPENMM_EXPORT int openmm_customnonbondedforce_addperparticleparameter_(OpenMM_CustomNonbondedForce*& target, const char* name, int name_length) {
    return OpenMM_CustomNonbondedForce_addPerParticleParameter(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT int OPENMM_CUSTOMNONBONDEDFORCE_ADDPERPARTICLEPARAMETER(OpenMM_CustomNonbondedForce*& target, const char* name, int name_length) {
    return OpenMM_CustomNonbondedForce_addPerParticleParameter(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void openmm_customnonbondedforce_getperparticleparametername_(const OpenMM_CustomNonbondedForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomNonbondedForce_getPerParticleParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_GETPERPARTICLEPARAMETERNAME(const OpenMM_CustomNonbondedForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomNonbondedForce_getPerParticleParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customnonbondedforce_setperparticleparametername_(OpenMM_CustomNonbondedForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomNonbondedForce_setPerParticleParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_SETPERPARTICLEPARAMETERNAME(OpenMM_CustomNonbondedForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomNonbondedForce_setPerParticleParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT int openmm_customnonbondedforce_addglobalparameter_(OpenMM_CustomNonbondedForce*& target, const char* name, double const& defaultValue, int name_length) {
    return OpenMM_CustomNonbondedForce_addGlobalParameter(target, makeString(name, name_length).c_str(), defaultValue);
}
OPENMM_EXPORT int OPENMM_CUSTOMNONBONDEDFORCE_ADDGLOBALPARAMETER(OpenMM_CustomNonbondedForce*& target, const char* name, double const& defaultValue, int name_length) {
    return OpenMM_CustomNonbondedForce_addGlobalParameter(target, makeString(name, name_length).c_str(), defaultValue);
}
OPENMM_EXPORT void openmm_customnonbondedforce_getglobalparametername_(const OpenMM_CustomNonbondedForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomNonbondedForce_getGlobalParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_GETGLOBALPARAMETERNAME(const OpenMM_CustomNonbondedForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomNonbondedForce_getGlobalParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customnonbondedforce_setglobalparametername_(OpenMM_CustomNonbondedForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomNonbondedForce_setGlobalParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_SETGLOBALPARAMETERNAME(OpenMM_CustomNonbondedForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomNonbondedForce_setGlobalParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT double openmm_customnonbondedforce_getglobalparameterdefaultvalue_(const OpenMM_CustomNonbondedForce*& target, int const& index) {
    return OpenMM_CustomNonbondedForce_getGlobalParameterDefaultValue(target, index);
}
OPENMM_EXPORT double OPENMM_CUSTOMNONBONDEDFORCE_GETGLOBALPARAMETERDEFAULTVALUE(const OpenMM_CustomNonbondedForce*& target, int const& index) {
    return OpenMM_CustomNonbondedForce_getGlobalParameterDefaultValue(target, index);
}
OPENMM_EXPORT void openmm_customnonbondedforce_setglobalparameterdefaultvalue_(OpenMM_CustomNonbondedForce*& target, int const& index, double const& defaultValue) {
    OpenMM_CustomNonbondedForce_setGlobalParameterDefaultValue(target, index, defaultValue);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_SETGLOBALPARAMETERDEFAULTVALUE(OpenMM_CustomNonbondedForce*& target, int const& index, double const& defaultValue) {
    OpenMM_CustomNonbondedForce_setGlobalParameterDefaultValue(target, index, defaultValue);
}
OPENMM_EXPORT void openmm_customnonbondedforce_addenergyparameterderivative_(OpenMM_CustomNonbondedForce*& target, const char* name, int name_length) {
    OpenMM_CustomNonbondedForce_addEnergyParameterDerivative(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_ADDENERGYPARAMETERDERIVATIVE(OpenMM_CustomNonbondedForce*& target, const char* name, int name_length) {
    OpenMM_CustomNonbondedForce_addEnergyParameterDerivative(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void openmm_customnonbondedforce_getenergyparameterderivativename_(const OpenMM_CustomNonbondedForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomNonbondedForce_getEnergyParameterDerivativeName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_GETENERGYPARAMETERDERIVATIVENAME(const OpenMM_CustomNonbondedForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomNonbondedForce_getEnergyParameterDerivativeName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT int openmm_customnonbondedforce_addparticle_(OpenMM_CustomNonbondedForce*& target, const OpenMM_DoubleArray*& parameters) {
    return OpenMM_CustomNonbondedForce_addParticle(target, parameters);
}
OPENMM_EXPORT int OPENMM_CUSTOMNONBONDEDFORCE_ADDPARTICLE(OpenMM_CustomNonbondedForce*& target, const OpenMM_DoubleArray*& parameters) {
    return OpenMM_CustomNonbondedForce_addParticle(target, parameters);
}
OPENMM_EXPORT void openmm_customnonbondedforce_getparticleparameters_(const OpenMM_CustomNonbondedForce*& target, int const& index, OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomNonbondedForce_getParticleParameters(target, index, parameters);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_GETPARTICLEPARAMETERS(const OpenMM_CustomNonbondedForce*& target, int const& index, OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomNonbondedForce_getParticleParameters(target, index, parameters);
}
OPENMM_EXPORT void openmm_customnonbondedforce_setparticleparameters_(OpenMM_CustomNonbondedForce*& target, int const& index, const OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomNonbondedForce_setParticleParameters(target, index, parameters);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_SETPARTICLEPARAMETERS(OpenMM_CustomNonbondedForce*& target, int const& index, const OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomNonbondedForce_setParticleParameters(target, index, parameters);
}
OPENMM_EXPORT int openmm_customnonbondedforce_addexclusion_(OpenMM_CustomNonbondedForce*& target, int const& particle1, int const& particle2) {
    return OpenMM_CustomNonbondedForce_addExclusion(target, particle1, particle2);
}
OPENMM_EXPORT int OPENMM_CUSTOMNONBONDEDFORCE_ADDEXCLUSION(OpenMM_CustomNonbondedForce*& target, int const& particle1, int const& particle2) {
    return OpenMM_CustomNonbondedForce_addExclusion(target, particle1, particle2);
}
OPENMM_EXPORT void openmm_customnonbondedforce_getexclusionparticles_(const OpenMM_CustomNonbondedForce*& target, int const& index, int* particle1, int* particle2) {
    OpenMM_CustomNonbondedForce_getExclusionParticles(target, index, particle1, particle2);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_GETEXCLUSIONPARTICLES(const OpenMM_CustomNonbondedForce*& target, int const& index, int* particle1, int* particle2) {
    OpenMM_CustomNonbondedForce_getExclusionParticles(target, index, particle1, particle2);
}
OPENMM_EXPORT void openmm_customnonbondedforce_setexclusionparticles_(OpenMM_CustomNonbondedForce*& target, int const& index, int const& particle1, int const& particle2) {
    OpenMM_CustomNonbondedForce_setExclusionParticles(target, index, particle1, particle2);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_SETEXCLUSIONPARTICLES(OpenMM_CustomNonbondedForce*& target, int const& index, int const& particle1, int const& particle2) {
    OpenMM_CustomNonbondedForce_setExclusionParticles(target, index, particle1, particle2);
}
OPENMM_EXPORT void openmm_customnonbondedforce_createexclusionsfrombonds_(OpenMM_CustomNonbondedForce*& target, const OpenMM_BondArray*& bonds, int const& bondCutoff) {
    OpenMM_CustomNonbondedForce_createExclusionsFromBonds(target, bonds, bondCutoff);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_CREATEEXCLUSIONSFROMBONDS(OpenMM_CustomNonbondedForce*& target, const OpenMM_BondArray*& bonds, int const& bondCutoff) {
    OpenMM_CustomNonbondedForce_createExclusionsFromBonds(target, bonds, bondCutoff);
}
OPENMM_EXPORT int openmm_customnonbondedforce_addtabulatedfunction_(OpenMM_CustomNonbondedForce*& target, const char* name, OpenMM_TabulatedFunction*& function, int name_length) {
    return OpenMM_CustomNonbondedForce_addTabulatedFunction(target, makeString(name, name_length).c_str(), function);
}
OPENMM_EXPORT int OPENMM_CUSTOMNONBONDEDFORCE_ADDTABULATEDFUNCTION(OpenMM_CustomNonbondedForce*& target, const char* name, OpenMM_TabulatedFunction*& function, int name_length) {
    return OpenMM_CustomNonbondedForce_addTabulatedFunction(target, makeString(name, name_length).c_str(), function);
}
OPENMM_EXPORT void openmm_customnonbondedforce_gettabulatedfunction_(OpenMM_CustomNonbondedForce*& target, int const& index, OpenMM_TabulatedFunction*& result) {
    result = OpenMM_CustomNonbondedForce_getTabulatedFunction(target, index);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_GETTABULATEDFUNCTION(OpenMM_CustomNonbondedForce*& target, int const& index, OpenMM_TabulatedFunction*& result) {
    result = OpenMM_CustomNonbondedForce_getTabulatedFunction(target, index);
}
OPENMM_EXPORT void openmm_customnonbondedforce_gettabulatedfunctionname_(const OpenMM_CustomNonbondedForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomNonbondedForce_getTabulatedFunctionName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_GETTABULATEDFUNCTIONNAME(const OpenMM_CustomNonbondedForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomNonbondedForce_getTabulatedFunctionName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT int openmm_customnonbondedforce_addfunction_(OpenMM_CustomNonbondedForce*& target, const char* name, const OpenMM_DoubleArray*& values, double const& min, double const& max, int name_length) {
    return OpenMM_CustomNonbondedForce_addFunction(target, makeString(name, name_length).c_str(), values, min, max);
}
OPENMM_EXPORT int OPENMM_CUSTOMNONBONDEDFORCE_ADDFUNCTION(OpenMM_CustomNonbondedForce*& target, const char* name, const OpenMM_DoubleArray*& values, double const& min, double const& max, int name_length) {
    return OpenMM_CustomNonbondedForce_addFunction(target, makeString(name, name_length).c_str(), values, min, max);
}
OPENMM_EXPORT void openmm_customnonbondedforce_getfunctionparameters_(const OpenMM_CustomNonbondedForce*& target, int const& index, char** name, OpenMM_DoubleArray*& values, double* min, double* max) {
    OpenMM_CustomNonbondedForce_getFunctionParameters(target, index, name, values, min, max);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_GETFUNCTIONPARAMETERS(const OpenMM_CustomNonbondedForce*& target, int const& index, char** name, OpenMM_DoubleArray*& values, double* min, double* max) {
    OpenMM_CustomNonbondedForce_getFunctionParameters(target, index, name, values, min, max);
}
OPENMM_EXPORT void openmm_customnonbondedforce_setfunctionparameters_(OpenMM_CustomNonbondedForce*& target, int const& index, const char* name, const OpenMM_DoubleArray*& values, double const& min, double const& max, int name_length) {
    OpenMM_CustomNonbondedForce_setFunctionParameters(target, index, makeString(name, name_length).c_str(), values, min, max);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_SETFUNCTIONPARAMETERS(OpenMM_CustomNonbondedForce*& target, int const& index, const char* name, const OpenMM_DoubleArray*& values, double const& min, double const& max, int name_length) {
    OpenMM_CustomNonbondedForce_setFunctionParameters(target, index, makeString(name, name_length).c_str(), values, min, max);
}
OPENMM_EXPORT int openmm_customnonbondedforce_addinteractiongroup_(OpenMM_CustomNonbondedForce*& target, const OpenMM_IntSet*& set1, const OpenMM_IntSet*& set2) {
    return OpenMM_CustomNonbondedForce_addInteractionGroup(target, set1, set2);
}
OPENMM_EXPORT int OPENMM_CUSTOMNONBONDEDFORCE_ADDINTERACTIONGROUP(OpenMM_CustomNonbondedForce*& target, const OpenMM_IntSet*& set1, const OpenMM_IntSet*& set2) {
    return OpenMM_CustomNonbondedForce_addInteractionGroup(target, set1, set2);
}
OPENMM_EXPORT void openmm_customnonbondedforce_getinteractiongroupparameters_(const OpenMM_CustomNonbondedForce*& target, int const& index, OpenMM_IntSet*& set1, OpenMM_IntSet*& set2) {
    OpenMM_CustomNonbondedForce_getInteractionGroupParameters(target, index, set1, set2);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_GETINTERACTIONGROUPPARAMETERS(const OpenMM_CustomNonbondedForce*& target, int const& index, OpenMM_IntSet*& set1, OpenMM_IntSet*& set2) {
    OpenMM_CustomNonbondedForce_getInteractionGroupParameters(target, index, set1, set2);
}
OPENMM_EXPORT void openmm_customnonbondedforce_setinteractiongroupparameters_(OpenMM_CustomNonbondedForce*& target, int const& index, const OpenMM_IntSet*& set1, const OpenMM_IntSet*& set2) {
    OpenMM_CustomNonbondedForce_setInteractionGroupParameters(target, index, set1, set2);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_SETINTERACTIONGROUPPARAMETERS(OpenMM_CustomNonbondedForce*& target, int const& index, const OpenMM_IntSet*& set1, const OpenMM_IntSet*& set2) {
    OpenMM_CustomNonbondedForce_setInteractionGroupParameters(target, index, set1, set2);
}
OPENMM_EXPORT void openmm_customnonbondedforce_updateparametersincontext_(OpenMM_CustomNonbondedForce*& target, OpenMM_Context*& context) {
    OpenMM_CustomNonbondedForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_CustomNonbondedForce*& target, OpenMM_Context*& context) {
    OpenMM_CustomNonbondedForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void openmm_customnonbondedforce_usesperiodicboundaryconditions_(const OpenMM_CustomNonbondedForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CustomNonbondedForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_CUSTOMNONBONDEDFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_CustomNonbondedForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CustomNonbondedForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::System */
OPENMM_EXPORT void openmm_system_create_(OpenMM_System*& result) {
    result = OpenMM_System_create();
}
OPENMM_EXPORT void OPENMM_SYSTEM_CREATE(OpenMM_System*& result) {
    result = OpenMM_System_create();
}
OPENMM_EXPORT void openmm_system_destroy_(OpenMM_System*& destroy) {
    OpenMM_System_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_SYSTEM_DESTROY(OpenMM_System*& destroy) {
    OpenMM_System_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT int openmm_system_getnumparticles_(const OpenMM_System*& target) {
    return OpenMM_System_getNumParticles(target);
}
OPENMM_EXPORT int OPENMM_SYSTEM_GETNUMPARTICLES(const OpenMM_System*& target) {
    return OpenMM_System_getNumParticles(target);
}
OPENMM_EXPORT int openmm_system_addparticle_(OpenMM_System*& target, double const& mass) {
    return OpenMM_System_addParticle(target, mass);
}
OPENMM_EXPORT int OPENMM_SYSTEM_ADDPARTICLE(OpenMM_System*& target, double const& mass) {
    return OpenMM_System_addParticle(target, mass);
}
OPENMM_EXPORT double openmm_system_getparticlemass_(const OpenMM_System*& target, int const& index) {
    return OpenMM_System_getParticleMass(target, index);
}
OPENMM_EXPORT double OPENMM_SYSTEM_GETPARTICLEMASS(const OpenMM_System*& target, int const& index) {
    return OpenMM_System_getParticleMass(target, index);
}
OPENMM_EXPORT void openmm_system_setparticlemass_(OpenMM_System*& target, int const& index, double const& mass) {
    OpenMM_System_setParticleMass(target, index, mass);
}
OPENMM_EXPORT void OPENMM_SYSTEM_SETPARTICLEMASS(OpenMM_System*& target, int const& index, double const& mass) {
    OpenMM_System_setParticleMass(target, index, mass);
}
OPENMM_EXPORT void openmm_system_setvirtualsite_(OpenMM_System*& target, int const& index, OpenMM_VirtualSite*& virtualSite) {
    OpenMM_System_setVirtualSite(target, index, virtualSite);
}
OPENMM_EXPORT void OPENMM_SYSTEM_SETVIRTUALSITE(OpenMM_System*& target, int const& index, OpenMM_VirtualSite*& virtualSite) {
    OpenMM_System_setVirtualSite(target, index, virtualSite);
}
OPENMM_EXPORT void openmm_system_isvirtualsite_(const OpenMM_System*& target, int const& index, OpenMM_Boolean& result) {
    result = OpenMM_System_isVirtualSite(target, index);
}
OPENMM_EXPORT void OPENMM_SYSTEM_ISVIRTUALSITE(const OpenMM_System*& target, int const& index, OpenMM_Boolean& result) {
    result = OpenMM_System_isVirtualSite(target, index);
}
OPENMM_EXPORT void openmm_system_getvirtualsite_(const OpenMM_System*& target, int const& index, const OpenMM_VirtualSite*& result) {
    result = OpenMM_System_getVirtualSite(target, index);
}
OPENMM_EXPORT void OPENMM_SYSTEM_GETVIRTUALSITE(const OpenMM_System*& target, int const& index, const OpenMM_VirtualSite*& result) {
    result = OpenMM_System_getVirtualSite(target, index);
}
OPENMM_EXPORT int openmm_system_getnumconstraints_(const OpenMM_System*& target) {
    return OpenMM_System_getNumConstraints(target);
}
OPENMM_EXPORT int OPENMM_SYSTEM_GETNUMCONSTRAINTS(const OpenMM_System*& target) {
    return OpenMM_System_getNumConstraints(target);
}
OPENMM_EXPORT int openmm_system_addconstraint_(OpenMM_System*& target, int const& particle1, int const& particle2, double const& distance) {
    return OpenMM_System_addConstraint(target, particle1, particle2, distance);
}
OPENMM_EXPORT int OPENMM_SYSTEM_ADDCONSTRAINT(OpenMM_System*& target, int const& particle1, int const& particle2, double const& distance) {
    return OpenMM_System_addConstraint(target, particle1, particle2, distance);
}
OPENMM_EXPORT void openmm_system_getconstraintparameters_(const OpenMM_System*& target, int const& index, int* particle1, int* particle2, double* distance) {
    OpenMM_System_getConstraintParameters(target, index, particle1, particle2, distance);
}
OPENMM_EXPORT void OPENMM_SYSTEM_GETCONSTRAINTPARAMETERS(const OpenMM_System*& target, int const& index, int* particle1, int* particle2, double* distance) {
    OpenMM_System_getConstraintParameters(target, index, particle1, particle2, distance);
}
OPENMM_EXPORT void openmm_system_setconstraintparameters_(OpenMM_System*& target, int const& index, int const& particle1, int const& particle2, double const& distance) {
    OpenMM_System_setConstraintParameters(target, index, particle1, particle2, distance);
}
OPENMM_EXPORT void OPENMM_SYSTEM_SETCONSTRAINTPARAMETERS(OpenMM_System*& target, int const& index, int const& particle1, int const& particle2, double const& distance) {
    OpenMM_System_setConstraintParameters(target, index, particle1, particle2, distance);
}
OPENMM_EXPORT void openmm_system_removeconstraint_(OpenMM_System*& target, int const& index) {
    OpenMM_System_removeConstraint(target, index);
}
OPENMM_EXPORT void OPENMM_SYSTEM_REMOVECONSTRAINT(OpenMM_System*& target, int const& index) {
    OpenMM_System_removeConstraint(target, index);
}
OPENMM_EXPORT int openmm_system_addforce_(OpenMM_System*& target, OpenMM_Force*& force) {
    return OpenMM_System_addForce(target, force);
}
OPENMM_EXPORT int OPENMM_SYSTEM_ADDFORCE(OpenMM_System*& target, OpenMM_Force*& force) {
    return OpenMM_System_addForce(target, force);
}
OPENMM_EXPORT int openmm_system_getnumforces_(const OpenMM_System*& target) {
    return OpenMM_System_getNumForces(target);
}
OPENMM_EXPORT int OPENMM_SYSTEM_GETNUMFORCES(const OpenMM_System*& target) {
    return OpenMM_System_getNumForces(target);
}
OPENMM_EXPORT void openmm_system_getforce_(OpenMM_System*& target, int const& index, OpenMM_Force*& result) {
    result = OpenMM_System_getForce(target, index);
}
OPENMM_EXPORT void OPENMM_SYSTEM_GETFORCE(OpenMM_System*& target, int const& index, OpenMM_Force*& result) {
    result = OpenMM_System_getForce(target, index);
}
OPENMM_EXPORT void openmm_system_removeforce_(OpenMM_System*& target, int const& index) {
    OpenMM_System_removeForce(target, index);
}
OPENMM_EXPORT void OPENMM_SYSTEM_REMOVEFORCE(OpenMM_System*& target, int const& index) {
    OpenMM_System_removeForce(target, index);
}
OPENMM_EXPORT void openmm_system_getdefaultperiodicboxvectors_(const OpenMM_System*& target, OpenMM_Vec3* a, OpenMM_Vec3* b, OpenMM_Vec3* c) {
    OpenMM_System_getDefaultPeriodicBoxVectors(target, a, b, c);
}
OPENMM_EXPORT void OPENMM_SYSTEM_GETDEFAULTPERIODICBOXVECTORS(const OpenMM_System*& target, OpenMM_Vec3* a, OpenMM_Vec3* b, OpenMM_Vec3* c) {
    OpenMM_System_getDefaultPeriodicBoxVectors(target, a, b, c);
}
OPENMM_EXPORT void openmm_system_setdefaultperiodicboxvectors_(OpenMM_System*& target, const OpenMM_Vec3* a, const OpenMM_Vec3* b, const OpenMM_Vec3* c) {
    OpenMM_System_setDefaultPeriodicBoxVectors(target, a, b, c);
}
OPENMM_EXPORT void OPENMM_SYSTEM_SETDEFAULTPERIODICBOXVECTORS(OpenMM_System*& target, const OpenMM_Vec3* a, const OpenMM_Vec3* b, const OpenMM_Vec3* c) {
    OpenMM_System_setDefaultPeriodicBoxVectors(target, a, b, c);
}
OPENMM_EXPORT void openmm_system_usesperiodicboundaryconditions_(const OpenMM_System*& target, OpenMM_Boolean& result) {
    result = OpenMM_System_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_SYSTEM_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_System*& target, OpenMM_Boolean& result) {
    result = OpenMM_System_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::Continuous1DFunction */
OPENMM_EXPORT void openmm_continuous1dfunction_create_(OpenMM_Continuous1DFunction*& result, const OpenMM_DoubleArray*& values, double const& min, double const& max) {
    result = OpenMM_Continuous1DFunction_create(values, min, max);
}
OPENMM_EXPORT void OPENMM_CONTINUOUS1DFUNCTION_CREATE(OpenMM_Continuous1DFunction*& result, const OpenMM_DoubleArray*& values, double const& min, double const& max) {
    result = OpenMM_Continuous1DFunction_create(values, min, max);
}
OPENMM_EXPORT void openmm_continuous1dfunction_destroy_(OpenMM_Continuous1DFunction*& destroy) {
    OpenMM_Continuous1DFunction_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_CONTINUOUS1DFUNCTION_DESTROY(OpenMM_Continuous1DFunction*& destroy) {
    OpenMM_Continuous1DFunction_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void openmm_continuous1dfunction_getfunctionparameters_(const OpenMM_Continuous1DFunction*& target, OpenMM_DoubleArray*& values, double* min, double* max) {
    OpenMM_Continuous1DFunction_getFunctionParameters(target, values, min, max);
}
OPENMM_EXPORT void OPENMM_CONTINUOUS1DFUNCTION_GETFUNCTIONPARAMETERS(const OpenMM_Continuous1DFunction*& target, OpenMM_DoubleArray*& values, double* min, double* max) {
    OpenMM_Continuous1DFunction_getFunctionParameters(target, values, min, max);
}
OPENMM_EXPORT void openmm_continuous1dfunction_setfunctionparameters_(OpenMM_Continuous1DFunction*& target, const OpenMM_DoubleArray*& values, double const& min, double const& max) {
    OpenMM_Continuous1DFunction_setFunctionParameters(target, values, min, max);
}
OPENMM_EXPORT void OPENMM_CONTINUOUS1DFUNCTION_SETFUNCTIONPARAMETERS(OpenMM_Continuous1DFunction*& target, const OpenMM_DoubleArray*& values, double const& min, double const& max) {
    OpenMM_Continuous1DFunction_setFunctionParameters(target, values, min, max);
}
OPENMM_EXPORT void openmm_continuous1dfunction_copy_(const OpenMM_Continuous1DFunction*& target, OpenMM_Continuous1DFunction*& result) {
    result = OpenMM_Continuous1DFunction_Copy(target);
}
OPENMM_EXPORT void OPENMM_CONTINUOUS1DFUNCTION_COPY(const OpenMM_Continuous1DFunction*& target, OpenMM_Continuous1DFunction*& result) {
    result = OpenMM_Continuous1DFunction_Copy(target);
}

/* OpenMM::Platform */
OPENMM_EXPORT void openmm_platform_destroy_(OpenMM_Platform*& destroy) {
    OpenMM_Platform_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_PLATFORM_DESTROY(OpenMM_Platform*& destroy) {
    OpenMM_Platform_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void openmm_platform_registerplatform_(OpenMM_Platform*& platform) {
    OpenMM_Platform_registerPlatform(platform);
}
OPENMM_EXPORT void OPENMM_PLATFORM_REGISTERPLATFORM(OpenMM_Platform*& platform) {
    OpenMM_Platform_registerPlatform(platform);
}
OPENMM_EXPORT int openmm_platform_getnumplatforms_() {
    return OpenMM_Platform_getNumPlatforms();
}
OPENMM_EXPORT int OPENMM_PLATFORM_GETNUMPLATFORMS() {
    return OpenMM_Platform_getNumPlatforms();
}
OPENMM_EXPORT void openmm_platform_getplatform_(int const& index, OpenMM_Platform*& result) {
    result = OpenMM_Platform_getPlatform(index);
}
OPENMM_EXPORT void OPENMM_PLATFORM_GETPLATFORM(int const& index, OpenMM_Platform*& result) {
    result = OpenMM_Platform_getPlatform(index);
}
OPENMM_EXPORT void openmm_platform_getplatformbyname_(const char* name, OpenMM_Platform*& result, int name_length) {
    result = OpenMM_Platform_getPlatformByName(makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_PLATFORM_GETPLATFORMBYNAME(const char* name, OpenMM_Platform*& result, int name_length) {
    result = OpenMM_Platform_getPlatformByName(makeString(name, name_length).c_str());
}
OPENMM_EXPORT void openmm_platform_findplatform_(const OpenMM_StringArray*& kernelNames, OpenMM_Platform*& result) {
    result = OpenMM_Platform_findPlatform(kernelNames);
}
OPENMM_EXPORT void OPENMM_PLATFORM_FINDPLATFORM(const OpenMM_StringArray*& kernelNames, OpenMM_Platform*& result) {
    result = OpenMM_Platform_findPlatform(kernelNames);
}
OPENMM_EXPORT void openmm_platform_loadpluginlibrary_(const char* file, int file_length) {
    OpenMM_Platform_loadPluginLibrary(makeString(file, file_length).c_str());
}
OPENMM_EXPORT void OPENMM_PLATFORM_LOADPLUGINLIBRARY(const char* file, int file_length) {
    OpenMM_Platform_loadPluginLibrary(makeString(file, file_length).c_str());
}
OPENMM_EXPORT void openmm_platform_getdefaultpluginsdirectory_(char* result, int result_length) {
    const char* result_chars = OpenMM_Platform_getDefaultPluginsDirectory();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_PLATFORM_GETDEFAULTPLUGINSDIRECTORY(char* result, int result_length) {
    const char* result_chars = OpenMM_Platform_getDefaultPluginsDirectory();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_platform_getopenmmversion_(char* result, int result_length) {
    const char* result_chars = OpenMM_Platform_getOpenMMVersion();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_PLATFORM_GETOPENMMVERSION(char* result, int result_length) {
    const char* result_chars = OpenMM_Platform_getOpenMMVersion();
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_platform_getname_(const OpenMM_Platform*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_Platform_getName(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_PLATFORM_GETNAME(const OpenMM_Platform*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_Platform_getName(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT double openmm_platform_getspeed_(const OpenMM_Platform*& target) {
    return OpenMM_Platform_getSpeed(target);
}
OPENMM_EXPORT double OPENMM_PLATFORM_GETSPEED(const OpenMM_Platform*& target) {
    return OpenMM_Platform_getSpeed(target);
}
OPENMM_EXPORT void openmm_platform_supportsdoubleprecision_(const OpenMM_Platform*& target, OpenMM_Boolean& result) {
    result = OpenMM_Platform_supportsDoublePrecision(target);
}
OPENMM_EXPORT void OPENMM_PLATFORM_SUPPORTSDOUBLEPRECISION(const OpenMM_Platform*& target, OpenMM_Boolean& result) {
    result = OpenMM_Platform_supportsDoublePrecision(target);
}
OPENMM_EXPORT void openmm_platform_getpropertynames_(const OpenMM_Platform*& target, const OpenMM_StringArray*& result) {
    result = OpenMM_Platform_getPropertyNames(target);
}
OPENMM_EXPORT void OPENMM_PLATFORM_GETPROPERTYNAMES(const OpenMM_Platform*& target, const OpenMM_StringArray*& result) {
    result = OpenMM_Platform_getPropertyNames(target);
}
OPENMM_EXPORT void openmm_platform_getpropertyvalue_(const OpenMM_Platform*& target, const OpenMM_Context*& context, const char* property, char* result, int property_length, int result_length) {
    const char* result_chars = OpenMM_Platform_getPropertyValue(target, context, makeString(property, property_length).c_str());
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_PLATFORM_GETPROPERTYVALUE(const OpenMM_Platform*& target, const OpenMM_Context*& context, const char* property, char* result, int property_length, int result_length) {
    const char* result_chars = OpenMM_Platform_getPropertyValue(target, context, makeString(property, property_length).c_str());
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_platform_setpropertyvalue_(const OpenMM_Platform*& target, OpenMM_Context*& context, const char* property, const char* value, int property_length, int value_length) {
    OpenMM_Platform_setPropertyValue(target, context, makeString(property, property_length).c_str(), makeString(value, value_length).c_str());
}
OPENMM_EXPORT void OPENMM_PLATFORM_SETPROPERTYVALUE(const OpenMM_Platform*& target, OpenMM_Context*& context, const char* property, const char* value, int property_length, int value_length) {
    OpenMM_Platform_setPropertyValue(target, context, makeString(property, property_length).c_str(), makeString(value, value_length).c_str());
}
OPENMM_EXPORT void openmm_platform_getpropertydefaultvalue_(const OpenMM_Platform*& target, const char* property, char* result, int property_length, int result_length) {
    const char* result_chars = OpenMM_Platform_getPropertyDefaultValue(target, makeString(property, property_length).c_str());
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_PLATFORM_GETPROPERTYDEFAULTVALUE(const OpenMM_Platform*& target, const char* property, char* result, int property_length, int result_length) {
    const char* result_chars = OpenMM_Platform_getPropertyDefaultValue(target, makeString(property, property_length).c_str());
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_platform_setpropertydefaultvalue_(OpenMM_Platform*& target, const char* property, const char* value, int property_length, int value_length) {
    OpenMM_Platform_setPropertyDefaultValue(target, makeString(property, property_length).c_str(), makeString(value, value_length).c_str());
}
OPENMM_EXPORT void OPENMM_PLATFORM_SETPROPERTYDEFAULTVALUE(OpenMM_Platform*& target, const char* property, const char* value, int property_length, int value_length) {
    OpenMM_Platform_setPropertyDefaultValue(target, makeString(property, property_length).c_str(), makeString(value, value_length).c_str());
}
OPENMM_EXPORT void openmm_platform_supportskernels_(const OpenMM_Platform*& target, const OpenMM_StringArray*& kernelNames, OpenMM_Boolean& result) {
    result = OpenMM_Platform_supportsKernels(target, kernelNames);
}
OPENMM_EXPORT void OPENMM_PLATFORM_SUPPORTSKERNELS(const OpenMM_Platform*& target, const OpenMM_StringArray*& kernelNames, OpenMM_Boolean& result) {
    result = OpenMM_Platform_supportsKernels(target, kernelNames);
}

/* OpenMM::Continuous3DFunction */
OPENMM_EXPORT void openmm_continuous3dfunction_create_(OpenMM_Continuous3DFunction*& result, int const& xsize, int const& ysize, int const& zsize, const OpenMM_DoubleArray*& values, double const& xmin, double const& xmax, double const& ymin, double const& ymax, double const& zmin, double const& zmax) {
    result = OpenMM_Continuous3DFunction_create(xsize, ysize, zsize, values, xmin, xmax, ymin, ymax, zmin, zmax);
}
OPENMM_EXPORT void OPENMM_CONTINUOUS3DFUNCTION_CREATE(OpenMM_Continuous3DFunction*& result, int const& xsize, int const& ysize, int const& zsize, const OpenMM_DoubleArray*& values, double const& xmin, double const& xmax, double const& ymin, double const& ymax, double const& zmin, double const& zmax) {
    result = OpenMM_Continuous3DFunction_create(xsize, ysize, zsize, values, xmin, xmax, ymin, ymax, zmin, zmax);
}
OPENMM_EXPORT void openmm_continuous3dfunction_destroy_(OpenMM_Continuous3DFunction*& destroy) {
    OpenMM_Continuous3DFunction_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_CONTINUOUS3DFUNCTION_DESTROY(OpenMM_Continuous3DFunction*& destroy) {
    OpenMM_Continuous3DFunction_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void openmm_continuous3dfunction_getfunctionparameters_(const OpenMM_Continuous3DFunction*& target, int* xsize, int* ysize, int* zsize, OpenMM_DoubleArray*& values, double* xmin, double* xmax, double* ymin, double* ymax, double* zmin, double* zmax) {
    OpenMM_Continuous3DFunction_getFunctionParameters(target, xsize, ysize, zsize, values, xmin, xmax, ymin, ymax, zmin, zmax);
}
OPENMM_EXPORT void OPENMM_CONTINUOUS3DFUNCTION_GETFUNCTIONPARAMETERS(const OpenMM_Continuous3DFunction*& target, int* xsize, int* ysize, int* zsize, OpenMM_DoubleArray*& values, double* xmin, double* xmax, double* ymin, double* ymax, double* zmin, double* zmax) {
    OpenMM_Continuous3DFunction_getFunctionParameters(target, xsize, ysize, zsize, values, xmin, xmax, ymin, ymax, zmin, zmax);
}
OPENMM_EXPORT void openmm_continuous3dfunction_setfunctionparameters_(OpenMM_Continuous3DFunction*& target, int const& xsize, int const& ysize, int const& zsize, const OpenMM_DoubleArray*& values, double const& xmin, double const& xmax, double const& ymin, double const& ymax, double const& zmin, double const& zmax) {
    OpenMM_Continuous3DFunction_setFunctionParameters(target, xsize, ysize, zsize, values, xmin, xmax, ymin, ymax, zmin, zmax);
}
OPENMM_EXPORT void OPENMM_CONTINUOUS3DFUNCTION_SETFUNCTIONPARAMETERS(OpenMM_Continuous3DFunction*& target, int const& xsize, int const& ysize, int const& zsize, const OpenMM_DoubleArray*& values, double const& xmin, double const& xmax, double const& ymin, double const& ymax, double const& zmin, double const& zmax) {
    OpenMM_Continuous3DFunction_setFunctionParameters(target, xsize, ysize, zsize, values, xmin, xmax, ymin, ymax, zmin, zmax);
}
OPENMM_EXPORT void openmm_continuous3dfunction_copy_(const OpenMM_Continuous3DFunction*& target, OpenMM_Continuous3DFunction*& result) {
    result = OpenMM_Continuous3DFunction_Copy(target);
}
OPENMM_EXPORT void OPENMM_CONTINUOUS3DFUNCTION_COPY(const OpenMM_Continuous3DFunction*& target, OpenMM_Continuous3DFunction*& result) {
    result = OpenMM_Continuous3DFunction_Copy(target);
}

/* OpenMM::GBSAOBCForce */
OPENMM_EXPORT void openmm_gbsaobcforce_create_(OpenMM_GBSAOBCForce*& result) {
    result = OpenMM_GBSAOBCForce_create();
}
OPENMM_EXPORT void OPENMM_GBSAOBCFORCE_CREATE(OpenMM_GBSAOBCForce*& result) {
    result = OpenMM_GBSAOBCForce_create();
}
OPENMM_EXPORT void openmm_gbsaobcforce_destroy_(OpenMM_GBSAOBCForce*& destroy) {
    OpenMM_GBSAOBCForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_GBSAOBCFORCE_DESTROY(OpenMM_GBSAOBCForce*& destroy) {
    OpenMM_GBSAOBCForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT int openmm_gbsaobcforce_getnumparticles_(const OpenMM_GBSAOBCForce*& target) {
    return OpenMM_GBSAOBCForce_getNumParticles(target);
}
OPENMM_EXPORT int OPENMM_GBSAOBCFORCE_GETNUMPARTICLES(const OpenMM_GBSAOBCForce*& target) {
    return OpenMM_GBSAOBCForce_getNumParticles(target);
}
OPENMM_EXPORT int openmm_gbsaobcforce_addparticle_(OpenMM_GBSAOBCForce*& target, double const& charge, double const& radius, double const& scalingFactor) {
    return OpenMM_GBSAOBCForce_addParticle(target, charge, radius, scalingFactor);
}
OPENMM_EXPORT int OPENMM_GBSAOBCFORCE_ADDPARTICLE(OpenMM_GBSAOBCForce*& target, double const& charge, double const& radius, double const& scalingFactor) {
    return OpenMM_GBSAOBCForce_addParticle(target, charge, radius, scalingFactor);
}
OPENMM_EXPORT void openmm_gbsaobcforce_getparticleparameters_(const OpenMM_GBSAOBCForce*& target, int const& index, double* charge, double* radius, double* scalingFactor) {
    OpenMM_GBSAOBCForce_getParticleParameters(target, index, charge, radius, scalingFactor);
}
OPENMM_EXPORT void OPENMM_GBSAOBCFORCE_GETPARTICLEPARAMETERS(const OpenMM_GBSAOBCForce*& target, int const& index, double* charge, double* radius, double* scalingFactor) {
    OpenMM_GBSAOBCForce_getParticleParameters(target, index, charge, radius, scalingFactor);
}
OPENMM_EXPORT void openmm_gbsaobcforce_setparticleparameters_(OpenMM_GBSAOBCForce*& target, int const& index, double const& charge, double const& radius, double const& scalingFactor) {
    OpenMM_GBSAOBCForce_setParticleParameters(target, index, charge, radius, scalingFactor);
}
OPENMM_EXPORT void OPENMM_GBSAOBCFORCE_SETPARTICLEPARAMETERS(OpenMM_GBSAOBCForce*& target, int const& index, double const& charge, double const& radius, double const& scalingFactor) {
    OpenMM_GBSAOBCForce_setParticleParameters(target, index, charge, radius, scalingFactor);
}
OPENMM_EXPORT double openmm_gbsaobcforce_getsolventdielectric_(const OpenMM_GBSAOBCForce*& target) {
    return OpenMM_GBSAOBCForce_getSolventDielectric(target);
}
OPENMM_EXPORT double OPENMM_GBSAOBCFORCE_GETSOLVENTDIELECTRIC(const OpenMM_GBSAOBCForce*& target) {
    return OpenMM_GBSAOBCForce_getSolventDielectric(target);
}
OPENMM_EXPORT void openmm_gbsaobcforce_setsolventdielectric_(OpenMM_GBSAOBCForce*& target, double const& dielectric) {
    OpenMM_GBSAOBCForce_setSolventDielectric(target, dielectric);
}
OPENMM_EXPORT void OPENMM_GBSAOBCFORCE_SETSOLVENTDIELECTRIC(OpenMM_GBSAOBCForce*& target, double const& dielectric) {
    OpenMM_GBSAOBCForce_setSolventDielectric(target, dielectric);
}
OPENMM_EXPORT double openmm_gbsaobcforce_getsolutedielectric_(const OpenMM_GBSAOBCForce*& target) {
    return OpenMM_GBSAOBCForce_getSoluteDielectric(target);
}
OPENMM_EXPORT double OPENMM_GBSAOBCFORCE_GETSOLUTEDIELECTRIC(const OpenMM_GBSAOBCForce*& target) {
    return OpenMM_GBSAOBCForce_getSoluteDielectric(target);
}
OPENMM_EXPORT void openmm_gbsaobcforce_setsolutedielectric_(OpenMM_GBSAOBCForce*& target, double const& dielectric) {
    OpenMM_GBSAOBCForce_setSoluteDielectric(target, dielectric);
}
OPENMM_EXPORT void OPENMM_GBSAOBCFORCE_SETSOLUTEDIELECTRIC(OpenMM_GBSAOBCForce*& target, double const& dielectric) {
    OpenMM_GBSAOBCForce_setSoluteDielectric(target, dielectric);
}
OPENMM_EXPORT double openmm_gbsaobcforce_getsurfaceareaenergy_(const OpenMM_GBSAOBCForce*& target) {
    return OpenMM_GBSAOBCForce_getSurfaceAreaEnergy(target);
}
OPENMM_EXPORT double OPENMM_GBSAOBCFORCE_GETSURFACEAREAENERGY(const OpenMM_GBSAOBCForce*& target) {
    return OpenMM_GBSAOBCForce_getSurfaceAreaEnergy(target);
}
OPENMM_EXPORT void openmm_gbsaobcforce_setsurfaceareaenergy_(OpenMM_GBSAOBCForce*& target, double const& energy) {
    OpenMM_GBSAOBCForce_setSurfaceAreaEnergy(target, energy);
}
OPENMM_EXPORT void OPENMM_GBSAOBCFORCE_SETSURFACEAREAENERGY(OpenMM_GBSAOBCForce*& target, double const& energy) {
    OpenMM_GBSAOBCForce_setSurfaceAreaEnergy(target, energy);
}
OPENMM_EXPORT void openmm_gbsaobcforce_getnonbondedmethod_(const OpenMM_GBSAOBCForce*& target, OpenMM_GBSAOBCForce_NonbondedMethod& result) {
    result = OpenMM_GBSAOBCForce_getNonbondedMethod(target);
}
OPENMM_EXPORT void OPENMM_GBSAOBCFORCE_GETNONBONDEDMETHOD(const OpenMM_GBSAOBCForce*& target, OpenMM_GBSAOBCForce_NonbondedMethod& result) {
    result = OpenMM_GBSAOBCForce_getNonbondedMethod(target);
}
OPENMM_EXPORT void openmm_gbsaobcforce_setnonbondedmethod_(OpenMM_GBSAOBCForce*& target, OpenMM_GBSAOBCForce_NonbondedMethod& method) {
    OpenMM_GBSAOBCForce_setNonbondedMethod(target, method);
}
OPENMM_EXPORT void OPENMM_GBSAOBCFORCE_SETNONBONDEDMETHOD(OpenMM_GBSAOBCForce*& target, OpenMM_GBSAOBCForce_NonbondedMethod& method) {
    OpenMM_GBSAOBCForce_setNonbondedMethod(target, method);
}
OPENMM_EXPORT double openmm_gbsaobcforce_getcutoffdistance_(const OpenMM_GBSAOBCForce*& target) {
    return OpenMM_GBSAOBCForce_getCutoffDistance(target);
}
OPENMM_EXPORT double OPENMM_GBSAOBCFORCE_GETCUTOFFDISTANCE(const OpenMM_GBSAOBCForce*& target) {
    return OpenMM_GBSAOBCForce_getCutoffDistance(target);
}
OPENMM_EXPORT void openmm_gbsaobcforce_setcutoffdistance_(OpenMM_GBSAOBCForce*& target, double const& distance) {
    OpenMM_GBSAOBCForce_setCutoffDistance(target, distance);
}
OPENMM_EXPORT void OPENMM_GBSAOBCFORCE_SETCUTOFFDISTANCE(OpenMM_GBSAOBCForce*& target, double const& distance) {
    OpenMM_GBSAOBCForce_setCutoffDistance(target, distance);
}
OPENMM_EXPORT void openmm_gbsaobcforce_updateparametersincontext_(OpenMM_GBSAOBCForce*& target, OpenMM_Context*& context) {
    OpenMM_GBSAOBCForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void OPENMM_GBSAOBCFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_GBSAOBCForce*& target, OpenMM_Context*& context) {
    OpenMM_GBSAOBCForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void openmm_gbsaobcforce_usesperiodicboundaryconditions_(const OpenMM_GBSAOBCForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_GBSAOBCForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_GBSAOBCFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_GBSAOBCForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_GBSAOBCForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::CustomCentroidBondForce */
OPENMM_EXPORT void openmm_customcentroidbondforce_create_(OpenMM_CustomCentroidBondForce*& result, int const& numGroups, const char* energy, int energy_length) {
    result = OpenMM_CustomCentroidBondForce_create(numGroups, makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMCENTROIDBONDFORCE_CREATE(OpenMM_CustomCentroidBondForce*& result, int const& numGroups, const char* energy, int energy_length) {
    result = OpenMM_CustomCentroidBondForce_create(numGroups, makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void openmm_customcentroidbondforce_destroy_(OpenMM_CustomCentroidBondForce*& destroy) {
    OpenMM_CustomCentroidBondForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_CUSTOMCENTROIDBONDFORCE_DESTROY(OpenMM_CustomCentroidBondForce*& destroy) {
    OpenMM_CustomCentroidBondForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT int openmm_customcentroidbondforce_getnumgroupsperbond_(const OpenMM_CustomCentroidBondForce*& target) {
    return OpenMM_CustomCentroidBondForce_getNumGroupsPerBond(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMCENTROIDBONDFORCE_GETNUMGROUPSPERBOND(const OpenMM_CustomCentroidBondForce*& target) {
    return OpenMM_CustomCentroidBondForce_getNumGroupsPerBond(target);
}
OPENMM_EXPORT int openmm_customcentroidbondforce_getnumgroups_(const OpenMM_CustomCentroidBondForce*& target) {
    return OpenMM_CustomCentroidBondForce_getNumGroups(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMCENTROIDBONDFORCE_GETNUMGROUPS(const OpenMM_CustomCentroidBondForce*& target) {
    return OpenMM_CustomCentroidBondForce_getNumGroups(target);
}
OPENMM_EXPORT int openmm_customcentroidbondforce_getnumbonds_(const OpenMM_CustomCentroidBondForce*& target) {
    return OpenMM_CustomCentroidBondForce_getNumBonds(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMCENTROIDBONDFORCE_GETNUMBONDS(const OpenMM_CustomCentroidBondForce*& target) {
    return OpenMM_CustomCentroidBondForce_getNumBonds(target);
}
OPENMM_EXPORT int openmm_customcentroidbondforce_getnumperbondparameters_(const OpenMM_CustomCentroidBondForce*& target) {
    return OpenMM_CustomCentroidBondForce_getNumPerBondParameters(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMCENTROIDBONDFORCE_GETNUMPERBONDPARAMETERS(const OpenMM_CustomCentroidBondForce*& target) {
    return OpenMM_CustomCentroidBondForce_getNumPerBondParameters(target);
}
OPENMM_EXPORT int openmm_customcentroidbondforce_getnumglobalparameters_(const OpenMM_CustomCentroidBondForce*& target) {
    return OpenMM_CustomCentroidBondForce_getNumGlobalParameters(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMCENTROIDBONDFORCE_GETNUMGLOBALPARAMETERS(const OpenMM_CustomCentroidBondForce*& target) {
    return OpenMM_CustomCentroidBondForce_getNumGlobalParameters(target);
}
OPENMM_EXPORT int openmm_customcentroidbondforce_getnumenergyparameterderivatives_(const OpenMM_CustomCentroidBondForce*& target) {
    return OpenMM_CustomCentroidBondForce_getNumEnergyParameterDerivatives(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMCENTROIDBONDFORCE_GETNUMENERGYPARAMETERDERIVATIVES(const OpenMM_CustomCentroidBondForce*& target) {
    return OpenMM_CustomCentroidBondForce_getNumEnergyParameterDerivatives(target);
}
OPENMM_EXPORT int openmm_customcentroidbondforce_getnumtabulatedfunctions_(const OpenMM_CustomCentroidBondForce*& target) {
    return OpenMM_CustomCentroidBondForce_getNumTabulatedFunctions(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMCENTROIDBONDFORCE_GETNUMTABULATEDFUNCTIONS(const OpenMM_CustomCentroidBondForce*& target) {
    return OpenMM_CustomCentroidBondForce_getNumTabulatedFunctions(target);
}
OPENMM_EXPORT int openmm_customcentroidbondforce_getnumfunctions_(const OpenMM_CustomCentroidBondForce*& target) {
    return OpenMM_CustomCentroidBondForce_getNumFunctions(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMCENTROIDBONDFORCE_GETNUMFUNCTIONS(const OpenMM_CustomCentroidBondForce*& target) {
    return OpenMM_CustomCentroidBondForce_getNumFunctions(target);
}
OPENMM_EXPORT void openmm_customcentroidbondforce_getenergyfunction_(const OpenMM_CustomCentroidBondForce*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCentroidBondForce_getEnergyFunction(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMCENTROIDBONDFORCE_GETENERGYFUNCTION(const OpenMM_CustomCentroidBondForce*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCentroidBondForce_getEnergyFunction(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customcentroidbondforce_setenergyfunction_(OpenMM_CustomCentroidBondForce*& target, const char* energy, int energy_length) {
    OpenMM_CustomCentroidBondForce_setEnergyFunction(target, makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMCENTROIDBONDFORCE_SETENERGYFUNCTION(OpenMM_CustomCentroidBondForce*& target, const char* energy, int energy_length) {
    OpenMM_CustomCentroidBondForce_setEnergyFunction(target, makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT int openmm_customcentroidbondforce_addperbondparameter_(OpenMM_CustomCentroidBondForce*& target, const char* name, int name_length) {
    return OpenMM_CustomCentroidBondForce_addPerBondParameter(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT int OPENMM_CUSTOMCENTROIDBONDFORCE_ADDPERBONDPARAMETER(OpenMM_CustomCentroidBondForce*& target, const char* name, int name_length) {
    return OpenMM_CustomCentroidBondForce_addPerBondParameter(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void openmm_customcentroidbondforce_getperbondparametername_(const OpenMM_CustomCentroidBondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCentroidBondForce_getPerBondParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMCENTROIDBONDFORCE_GETPERBONDPARAMETERNAME(const OpenMM_CustomCentroidBondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCentroidBondForce_getPerBondParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customcentroidbondforce_setperbondparametername_(OpenMM_CustomCentroidBondForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomCentroidBondForce_setPerBondParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMCENTROIDBONDFORCE_SETPERBONDPARAMETERNAME(OpenMM_CustomCentroidBondForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomCentroidBondForce_setPerBondParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT int openmm_customcentroidbondforce_addglobalparameter_(OpenMM_CustomCentroidBondForce*& target, const char* name, double const& defaultValue, int name_length) {
    return OpenMM_CustomCentroidBondForce_addGlobalParameter(target, makeString(name, name_length).c_str(), defaultValue);
}
OPENMM_EXPORT int OPENMM_CUSTOMCENTROIDBONDFORCE_ADDGLOBALPARAMETER(OpenMM_CustomCentroidBondForce*& target, const char* name, double const& defaultValue, int name_length) {
    return OpenMM_CustomCentroidBondForce_addGlobalParameter(target, makeString(name, name_length).c_str(), defaultValue);
}
OPENMM_EXPORT void openmm_customcentroidbondforce_getglobalparametername_(const OpenMM_CustomCentroidBondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCentroidBondForce_getGlobalParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMCENTROIDBONDFORCE_GETGLOBALPARAMETERNAME(const OpenMM_CustomCentroidBondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCentroidBondForce_getGlobalParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customcentroidbondforce_setglobalparametername_(OpenMM_CustomCentroidBondForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomCentroidBondForce_setGlobalParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMCENTROIDBONDFORCE_SETGLOBALPARAMETERNAME(OpenMM_CustomCentroidBondForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomCentroidBondForce_setGlobalParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT double openmm_customcentroidbondforce_getglobalparameterdefaultvalue_(const OpenMM_CustomCentroidBondForce*& target, int const& index) {
    return OpenMM_CustomCentroidBondForce_getGlobalParameterDefaultValue(target, index);
}
OPENMM_EXPORT double OPENMM_CUSTOMCENTROIDBONDFORCE_GETGLOBALPARAMETERDEFAULTVALUE(const OpenMM_CustomCentroidBondForce*& target, int const& index) {
    return OpenMM_CustomCentroidBondForce_getGlobalParameterDefaultValue(target, index);
}
OPENMM_EXPORT void openmm_customcentroidbondforce_setglobalparameterdefaultvalue_(OpenMM_CustomCentroidBondForce*& target, int const& index, double const& defaultValue) {
    OpenMM_CustomCentroidBondForce_setGlobalParameterDefaultValue(target, index, defaultValue);
}
OPENMM_EXPORT void OPENMM_CUSTOMCENTROIDBONDFORCE_SETGLOBALPARAMETERDEFAULTVALUE(OpenMM_CustomCentroidBondForce*& target, int const& index, double const& defaultValue) {
    OpenMM_CustomCentroidBondForce_setGlobalParameterDefaultValue(target, index, defaultValue);
}
OPENMM_EXPORT void openmm_customcentroidbondforce_addenergyparameterderivative_(OpenMM_CustomCentroidBondForce*& target, const char* name, int name_length) {
    OpenMM_CustomCentroidBondForce_addEnergyParameterDerivative(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMCENTROIDBONDFORCE_ADDENERGYPARAMETERDERIVATIVE(OpenMM_CustomCentroidBondForce*& target, const char* name, int name_length) {
    OpenMM_CustomCentroidBondForce_addEnergyParameterDerivative(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void openmm_customcentroidbondforce_getenergyparameterderivativename_(const OpenMM_CustomCentroidBondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCentroidBondForce_getEnergyParameterDerivativeName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMCENTROIDBONDFORCE_GETENERGYPARAMETERDERIVATIVENAME(const OpenMM_CustomCentroidBondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCentroidBondForce_getEnergyParameterDerivativeName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT int openmm_customcentroidbondforce_addgroup_(OpenMM_CustomCentroidBondForce*& target, const OpenMM_IntArray*& particles, const OpenMM_DoubleArray*& weights) {
    return OpenMM_CustomCentroidBondForce_addGroup(target, particles, weights);
}
OPENMM_EXPORT int OPENMM_CUSTOMCENTROIDBONDFORCE_ADDGROUP(OpenMM_CustomCentroidBondForce*& target, const OpenMM_IntArray*& particles, const OpenMM_DoubleArray*& weights) {
    return OpenMM_CustomCentroidBondForce_addGroup(target, particles, weights);
}
OPENMM_EXPORT void openmm_customcentroidbondforce_getgroupparameters_(const OpenMM_CustomCentroidBondForce*& target, int const& index, OpenMM_IntArray*& particles, OpenMM_DoubleArray*& weights) {
    OpenMM_CustomCentroidBondForce_getGroupParameters(target, index, particles, weights);
}
OPENMM_EXPORT void OPENMM_CUSTOMCENTROIDBONDFORCE_GETGROUPPARAMETERS(const OpenMM_CustomCentroidBondForce*& target, int const& index, OpenMM_IntArray*& particles, OpenMM_DoubleArray*& weights) {
    OpenMM_CustomCentroidBondForce_getGroupParameters(target, index, particles, weights);
}
OPENMM_EXPORT void openmm_customcentroidbondforce_setgroupparameters_(OpenMM_CustomCentroidBondForce*& target, int const& index, const OpenMM_IntArray*& particles, const OpenMM_DoubleArray*& weights) {
    OpenMM_CustomCentroidBondForce_setGroupParameters(target, index, particles, weights);
}
OPENMM_EXPORT void OPENMM_CUSTOMCENTROIDBONDFORCE_SETGROUPPARAMETERS(OpenMM_CustomCentroidBondForce*& target, int const& index, const OpenMM_IntArray*& particles, const OpenMM_DoubleArray*& weights) {
    OpenMM_CustomCentroidBondForce_setGroupParameters(target, index, particles, weights);
}
OPENMM_EXPORT int openmm_customcentroidbondforce_addbond_(OpenMM_CustomCentroidBondForce*& target, const OpenMM_IntArray*& groups, const OpenMM_DoubleArray*& parameters) {
    return OpenMM_CustomCentroidBondForce_addBond(target, groups, parameters);
}
OPENMM_EXPORT int OPENMM_CUSTOMCENTROIDBONDFORCE_ADDBOND(OpenMM_CustomCentroidBondForce*& target, const OpenMM_IntArray*& groups, const OpenMM_DoubleArray*& parameters) {
    return OpenMM_CustomCentroidBondForce_addBond(target, groups, parameters);
}
OPENMM_EXPORT void openmm_customcentroidbondforce_getbondparameters_(const OpenMM_CustomCentroidBondForce*& target, int const& index, OpenMM_IntArray*& groups, OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomCentroidBondForce_getBondParameters(target, index, groups, parameters);
}
OPENMM_EXPORT void OPENMM_CUSTOMCENTROIDBONDFORCE_GETBONDPARAMETERS(const OpenMM_CustomCentroidBondForce*& target, int const& index, OpenMM_IntArray*& groups, OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomCentroidBondForce_getBondParameters(target, index, groups, parameters);
}
OPENMM_EXPORT void openmm_customcentroidbondforce_setbondparameters_(OpenMM_CustomCentroidBondForce*& target, int const& index, const OpenMM_IntArray*& groups, const OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomCentroidBondForce_setBondParameters(target, index, groups, parameters);
}
OPENMM_EXPORT void OPENMM_CUSTOMCENTROIDBONDFORCE_SETBONDPARAMETERS(OpenMM_CustomCentroidBondForce*& target, int const& index, const OpenMM_IntArray*& groups, const OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomCentroidBondForce_setBondParameters(target, index, groups, parameters);
}
OPENMM_EXPORT int openmm_customcentroidbondforce_addtabulatedfunction_(OpenMM_CustomCentroidBondForce*& target, const char* name, OpenMM_TabulatedFunction*& function, int name_length) {
    return OpenMM_CustomCentroidBondForce_addTabulatedFunction(target, makeString(name, name_length).c_str(), function);
}
OPENMM_EXPORT int OPENMM_CUSTOMCENTROIDBONDFORCE_ADDTABULATEDFUNCTION(OpenMM_CustomCentroidBondForce*& target, const char* name, OpenMM_TabulatedFunction*& function, int name_length) {
    return OpenMM_CustomCentroidBondForce_addTabulatedFunction(target, makeString(name, name_length).c_str(), function);
}
OPENMM_EXPORT void openmm_customcentroidbondforce_gettabulatedfunction_(OpenMM_CustomCentroidBondForce*& target, int const& index, OpenMM_TabulatedFunction*& result) {
    result = OpenMM_CustomCentroidBondForce_getTabulatedFunction(target, index);
}
OPENMM_EXPORT void OPENMM_CUSTOMCENTROIDBONDFORCE_GETTABULATEDFUNCTION(OpenMM_CustomCentroidBondForce*& target, int const& index, OpenMM_TabulatedFunction*& result) {
    result = OpenMM_CustomCentroidBondForce_getTabulatedFunction(target, index);
}
OPENMM_EXPORT void openmm_customcentroidbondforce_gettabulatedfunctionname_(const OpenMM_CustomCentroidBondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCentroidBondForce_getTabulatedFunctionName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMCENTROIDBONDFORCE_GETTABULATEDFUNCTIONNAME(const OpenMM_CustomCentroidBondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCentroidBondForce_getTabulatedFunctionName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customcentroidbondforce_updateparametersincontext_(OpenMM_CustomCentroidBondForce*& target, OpenMM_Context*& context) {
    OpenMM_CustomCentroidBondForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void OPENMM_CUSTOMCENTROIDBONDFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_CustomCentroidBondForce*& target, OpenMM_Context*& context) {
    OpenMM_CustomCentroidBondForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void openmm_customcentroidbondforce_setusesperiodicboundarycondition_(OpenMM_CustomCentroidBondForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_CustomCentroidBondForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT void OPENMM_CUSTOMCENTROIDBONDFORCE_SETUSESPERIODICBOUNDARYCONDITION(OpenMM_CustomCentroidBondForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_CustomCentroidBondForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT void openmm_customcentroidbondforce_usesperiodicboundaryconditions_(const OpenMM_CustomCentroidBondForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CustomCentroidBondForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_CUSTOMCENTROIDBONDFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_CustomCentroidBondForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CustomCentroidBondForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::LocalEnergyMinimizer */
OPENMM_EXPORT void openmm_localenergyminimizer_destroy_(OpenMM_LocalEnergyMinimizer*& destroy) {
    OpenMM_LocalEnergyMinimizer_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_LOCALENERGYMINIMIZER_DESTROY(OpenMM_LocalEnergyMinimizer*& destroy) {
    OpenMM_LocalEnergyMinimizer_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void openmm_localenergyminimizer_minimize_(OpenMM_Context*& context, double const& tolerance, int const& maxIterations) {
    OpenMM_LocalEnergyMinimizer_minimize(context, tolerance, maxIterations);
}
OPENMM_EXPORT void OPENMM_LOCALENERGYMINIMIZER_MINIMIZE(OpenMM_Context*& context, double const& tolerance, int const& maxIterations) {
    OpenMM_LocalEnergyMinimizer_minimize(context, tolerance, maxIterations);
}

/* OpenMM::TwoParticleAverageSite */
OPENMM_EXPORT void openmm_twoparticleaveragesite_create_(OpenMM_TwoParticleAverageSite*& result, int const& particle1, int const& particle2, double const& weight1, double const& weight2) {
    result = OpenMM_TwoParticleAverageSite_create(particle1, particle2, weight1, weight2);
}
OPENMM_EXPORT void OPENMM_TWOPARTICLEAVERAGESITE_CREATE(OpenMM_TwoParticleAverageSite*& result, int const& particle1, int const& particle2, double const& weight1, double const& weight2) {
    result = OpenMM_TwoParticleAverageSite_create(particle1, particle2, weight1, weight2);
}
OPENMM_EXPORT void openmm_twoparticleaveragesite_destroy_(OpenMM_TwoParticleAverageSite*& destroy) {
    OpenMM_TwoParticleAverageSite_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_TWOPARTICLEAVERAGESITE_DESTROY(OpenMM_TwoParticleAverageSite*& destroy) {
    OpenMM_TwoParticleAverageSite_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT double openmm_twoparticleaveragesite_getweight_(const OpenMM_TwoParticleAverageSite*& target, int const& particle) {
    return OpenMM_TwoParticleAverageSite_getWeight(target, particle);
}
OPENMM_EXPORT double OPENMM_TWOPARTICLEAVERAGESITE_GETWEIGHT(const OpenMM_TwoParticleAverageSite*& target, int const& particle) {
    return OpenMM_TwoParticleAverageSite_getWeight(target, particle);
}

/* OpenMM::CustomCompoundBondForce */
OPENMM_EXPORT void openmm_customcompoundbondforce_create_(OpenMM_CustomCompoundBondForce*& result, int const& numParticles, const char* energy, int energy_length) {
    result = OpenMM_CustomCompoundBondForce_create(numParticles, makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMCOMPOUNDBONDFORCE_CREATE(OpenMM_CustomCompoundBondForce*& result, int const& numParticles, const char* energy, int energy_length) {
    result = OpenMM_CustomCompoundBondForce_create(numParticles, makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void openmm_customcompoundbondforce_destroy_(OpenMM_CustomCompoundBondForce*& destroy) {
    OpenMM_CustomCompoundBondForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_CUSTOMCOMPOUNDBONDFORCE_DESTROY(OpenMM_CustomCompoundBondForce*& destroy) {
    OpenMM_CustomCompoundBondForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT int openmm_customcompoundbondforce_getnumparticlesperbond_(const OpenMM_CustomCompoundBondForce*& target) {
    return OpenMM_CustomCompoundBondForce_getNumParticlesPerBond(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMCOMPOUNDBONDFORCE_GETNUMPARTICLESPERBOND(const OpenMM_CustomCompoundBondForce*& target) {
    return OpenMM_CustomCompoundBondForce_getNumParticlesPerBond(target);
}
OPENMM_EXPORT int openmm_customcompoundbondforce_getnumbonds_(const OpenMM_CustomCompoundBondForce*& target) {
    return OpenMM_CustomCompoundBondForce_getNumBonds(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMCOMPOUNDBONDFORCE_GETNUMBONDS(const OpenMM_CustomCompoundBondForce*& target) {
    return OpenMM_CustomCompoundBondForce_getNumBonds(target);
}
OPENMM_EXPORT int openmm_customcompoundbondforce_getnumperbondparameters_(const OpenMM_CustomCompoundBondForce*& target) {
    return OpenMM_CustomCompoundBondForce_getNumPerBondParameters(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMCOMPOUNDBONDFORCE_GETNUMPERBONDPARAMETERS(const OpenMM_CustomCompoundBondForce*& target) {
    return OpenMM_CustomCompoundBondForce_getNumPerBondParameters(target);
}
OPENMM_EXPORT int openmm_customcompoundbondforce_getnumglobalparameters_(const OpenMM_CustomCompoundBondForce*& target) {
    return OpenMM_CustomCompoundBondForce_getNumGlobalParameters(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMCOMPOUNDBONDFORCE_GETNUMGLOBALPARAMETERS(const OpenMM_CustomCompoundBondForce*& target) {
    return OpenMM_CustomCompoundBondForce_getNumGlobalParameters(target);
}
OPENMM_EXPORT int openmm_customcompoundbondforce_getnumenergyparameterderivatives_(const OpenMM_CustomCompoundBondForce*& target) {
    return OpenMM_CustomCompoundBondForce_getNumEnergyParameterDerivatives(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMCOMPOUNDBONDFORCE_GETNUMENERGYPARAMETERDERIVATIVES(const OpenMM_CustomCompoundBondForce*& target) {
    return OpenMM_CustomCompoundBondForce_getNumEnergyParameterDerivatives(target);
}
OPENMM_EXPORT int openmm_customcompoundbondforce_getnumtabulatedfunctions_(const OpenMM_CustomCompoundBondForce*& target) {
    return OpenMM_CustomCompoundBondForce_getNumTabulatedFunctions(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMCOMPOUNDBONDFORCE_GETNUMTABULATEDFUNCTIONS(const OpenMM_CustomCompoundBondForce*& target) {
    return OpenMM_CustomCompoundBondForce_getNumTabulatedFunctions(target);
}
OPENMM_EXPORT int openmm_customcompoundbondforce_getnumfunctions_(const OpenMM_CustomCompoundBondForce*& target) {
    return OpenMM_CustomCompoundBondForce_getNumFunctions(target);
}
OPENMM_EXPORT int OPENMM_CUSTOMCOMPOUNDBONDFORCE_GETNUMFUNCTIONS(const OpenMM_CustomCompoundBondForce*& target) {
    return OpenMM_CustomCompoundBondForce_getNumFunctions(target);
}
OPENMM_EXPORT void openmm_customcompoundbondforce_getenergyfunction_(const OpenMM_CustomCompoundBondForce*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCompoundBondForce_getEnergyFunction(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMCOMPOUNDBONDFORCE_GETENERGYFUNCTION(const OpenMM_CustomCompoundBondForce*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCompoundBondForce_getEnergyFunction(target);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customcompoundbondforce_setenergyfunction_(OpenMM_CustomCompoundBondForce*& target, const char* energy, int energy_length) {
    OpenMM_CustomCompoundBondForce_setEnergyFunction(target, makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMCOMPOUNDBONDFORCE_SETENERGYFUNCTION(OpenMM_CustomCompoundBondForce*& target, const char* energy, int energy_length) {
    OpenMM_CustomCompoundBondForce_setEnergyFunction(target, makeString(energy, energy_length).c_str());
}
OPENMM_EXPORT int openmm_customcompoundbondforce_addperbondparameter_(OpenMM_CustomCompoundBondForce*& target, const char* name, int name_length) {
    return OpenMM_CustomCompoundBondForce_addPerBondParameter(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT int OPENMM_CUSTOMCOMPOUNDBONDFORCE_ADDPERBONDPARAMETER(OpenMM_CustomCompoundBondForce*& target, const char* name, int name_length) {
    return OpenMM_CustomCompoundBondForce_addPerBondParameter(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void openmm_customcompoundbondforce_getperbondparametername_(const OpenMM_CustomCompoundBondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCompoundBondForce_getPerBondParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMCOMPOUNDBONDFORCE_GETPERBONDPARAMETERNAME(const OpenMM_CustomCompoundBondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCompoundBondForce_getPerBondParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customcompoundbondforce_setperbondparametername_(OpenMM_CustomCompoundBondForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomCompoundBondForce_setPerBondParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMCOMPOUNDBONDFORCE_SETPERBONDPARAMETERNAME(OpenMM_CustomCompoundBondForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomCompoundBondForce_setPerBondParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT int openmm_customcompoundbondforce_addglobalparameter_(OpenMM_CustomCompoundBondForce*& target, const char* name, double const& defaultValue, int name_length) {
    return OpenMM_CustomCompoundBondForce_addGlobalParameter(target, makeString(name, name_length).c_str(), defaultValue);
}
OPENMM_EXPORT int OPENMM_CUSTOMCOMPOUNDBONDFORCE_ADDGLOBALPARAMETER(OpenMM_CustomCompoundBondForce*& target, const char* name, double const& defaultValue, int name_length) {
    return OpenMM_CustomCompoundBondForce_addGlobalParameter(target, makeString(name, name_length).c_str(), defaultValue);
}
OPENMM_EXPORT void openmm_customcompoundbondforce_getglobalparametername_(const OpenMM_CustomCompoundBondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCompoundBondForce_getGlobalParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMCOMPOUNDBONDFORCE_GETGLOBALPARAMETERNAME(const OpenMM_CustomCompoundBondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCompoundBondForce_getGlobalParameterName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void openmm_customcompoundbondforce_setglobalparametername_(OpenMM_CustomCompoundBondForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomCompoundBondForce_setGlobalParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMCOMPOUNDBONDFORCE_SETGLOBALPARAMETERNAME(OpenMM_CustomCompoundBondForce*& target, int const& index, const char* name, int name_length) {
    OpenMM_CustomCompoundBondForce_setGlobalParameterName(target, index, makeString(name, name_length).c_str());
}
OPENMM_EXPORT double openmm_customcompoundbondforce_getglobalparameterdefaultvalue_(const OpenMM_CustomCompoundBondForce*& target, int const& index) {
    return OpenMM_CustomCompoundBondForce_getGlobalParameterDefaultValue(target, index);
}
OPENMM_EXPORT double OPENMM_CUSTOMCOMPOUNDBONDFORCE_GETGLOBALPARAMETERDEFAULTVALUE(const OpenMM_CustomCompoundBondForce*& target, int const& index) {
    return OpenMM_CustomCompoundBondForce_getGlobalParameterDefaultValue(target, index);
}
OPENMM_EXPORT void openmm_customcompoundbondforce_setglobalparameterdefaultvalue_(OpenMM_CustomCompoundBondForce*& target, int const& index, double const& defaultValue) {
    OpenMM_CustomCompoundBondForce_setGlobalParameterDefaultValue(target, index, defaultValue);
}
OPENMM_EXPORT void OPENMM_CUSTOMCOMPOUNDBONDFORCE_SETGLOBALPARAMETERDEFAULTVALUE(OpenMM_CustomCompoundBondForce*& target, int const& index, double const& defaultValue) {
    OpenMM_CustomCompoundBondForce_setGlobalParameterDefaultValue(target, index, defaultValue);
}
OPENMM_EXPORT void openmm_customcompoundbondforce_addenergyparameterderivative_(OpenMM_CustomCompoundBondForce*& target, const char* name, int name_length) {
    OpenMM_CustomCompoundBondForce_addEnergyParameterDerivative(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void OPENMM_CUSTOMCOMPOUNDBONDFORCE_ADDENERGYPARAMETERDERIVATIVE(OpenMM_CustomCompoundBondForce*& target, const char* name, int name_length) {
    OpenMM_CustomCompoundBondForce_addEnergyParameterDerivative(target, makeString(name, name_length).c_str());
}
OPENMM_EXPORT void openmm_customcompoundbondforce_getenergyparameterderivativename_(const OpenMM_CustomCompoundBondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCompoundBondForce_getEnergyParameterDerivativeName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMCOMPOUNDBONDFORCE_GETENERGYPARAMETERDERIVATIVENAME(const OpenMM_CustomCompoundBondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCompoundBondForce_getEnergyParameterDerivativeName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT int openmm_customcompoundbondforce_addbond_(OpenMM_CustomCompoundBondForce*& target, const OpenMM_IntArray*& particles, const OpenMM_DoubleArray*& parameters) {
    return OpenMM_CustomCompoundBondForce_addBond(target, particles, parameters);
}
OPENMM_EXPORT int OPENMM_CUSTOMCOMPOUNDBONDFORCE_ADDBOND(OpenMM_CustomCompoundBondForce*& target, const OpenMM_IntArray*& particles, const OpenMM_DoubleArray*& parameters) {
    return OpenMM_CustomCompoundBondForce_addBond(target, particles, parameters);
}
OPENMM_EXPORT void openmm_customcompoundbondforce_getbondparameters_(const OpenMM_CustomCompoundBondForce*& target, int const& index, OpenMM_IntArray*& particles, OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomCompoundBondForce_getBondParameters(target, index, particles, parameters);
}
OPENMM_EXPORT void OPENMM_CUSTOMCOMPOUNDBONDFORCE_GETBONDPARAMETERS(const OpenMM_CustomCompoundBondForce*& target, int const& index, OpenMM_IntArray*& particles, OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomCompoundBondForce_getBondParameters(target, index, particles, parameters);
}
OPENMM_EXPORT void openmm_customcompoundbondforce_setbondparameters_(OpenMM_CustomCompoundBondForce*& target, int const& index, const OpenMM_IntArray*& particles, const OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomCompoundBondForce_setBondParameters(target, index, particles, parameters);
}
OPENMM_EXPORT void OPENMM_CUSTOMCOMPOUNDBONDFORCE_SETBONDPARAMETERS(OpenMM_CustomCompoundBondForce*& target, int const& index, const OpenMM_IntArray*& particles, const OpenMM_DoubleArray*& parameters) {
    OpenMM_CustomCompoundBondForce_setBondParameters(target, index, particles, parameters);
}
OPENMM_EXPORT int openmm_customcompoundbondforce_addtabulatedfunction_(OpenMM_CustomCompoundBondForce*& target, const char* name, OpenMM_TabulatedFunction*& function, int name_length) {
    return OpenMM_CustomCompoundBondForce_addTabulatedFunction(target, makeString(name, name_length).c_str(), function);
}
OPENMM_EXPORT int OPENMM_CUSTOMCOMPOUNDBONDFORCE_ADDTABULATEDFUNCTION(OpenMM_CustomCompoundBondForce*& target, const char* name, OpenMM_TabulatedFunction*& function, int name_length) {
    return OpenMM_CustomCompoundBondForce_addTabulatedFunction(target, makeString(name, name_length).c_str(), function);
}
OPENMM_EXPORT void openmm_customcompoundbondforce_gettabulatedfunction_(OpenMM_CustomCompoundBondForce*& target, int const& index, OpenMM_TabulatedFunction*& result) {
    result = OpenMM_CustomCompoundBondForce_getTabulatedFunction(target, index);
}
OPENMM_EXPORT void OPENMM_CUSTOMCOMPOUNDBONDFORCE_GETTABULATEDFUNCTION(OpenMM_CustomCompoundBondForce*& target, int const& index, OpenMM_TabulatedFunction*& result) {
    result = OpenMM_CustomCompoundBondForce_getTabulatedFunction(target, index);
}
OPENMM_EXPORT void openmm_customcompoundbondforce_gettabulatedfunctionname_(const OpenMM_CustomCompoundBondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCompoundBondForce_getTabulatedFunctionName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT void OPENMM_CUSTOMCOMPOUNDBONDFORCE_GETTABULATEDFUNCTIONNAME(const OpenMM_CustomCompoundBondForce*& target, int const& index, char* result, int result_length) {
    const char* result_chars = OpenMM_CustomCompoundBondForce_getTabulatedFunctionName(target, index);
    copyAndPadString(result, result_chars, result_length);
}
OPENMM_EXPORT int openmm_customcompoundbondforce_addfunction_(OpenMM_CustomCompoundBondForce*& target, const char* name, const OpenMM_DoubleArray*& values, double const& min, double const& max, int name_length) {
    return OpenMM_CustomCompoundBondForce_addFunction(target, makeString(name, name_length).c_str(), values, min, max);
}
OPENMM_EXPORT int OPENMM_CUSTOMCOMPOUNDBONDFORCE_ADDFUNCTION(OpenMM_CustomCompoundBondForce*& target, const char* name, const OpenMM_DoubleArray*& values, double const& min, double const& max, int name_length) {
    return OpenMM_CustomCompoundBondForce_addFunction(target, makeString(name, name_length).c_str(), values, min, max);
}
OPENMM_EXPORT void openmm_customcompoundbondforce_getfunctionparameters_(const OpenMM_CustomCompoundBondForce*& target, int const& index, char** name, OpenMM_DoubleArray*& values, double* min, double* max) {
    OpenMM_CustomCompoundBondForce_getFunctionParameters(target, index, name, values, min, max);
}
OPENMM_EXPORT void OPENMM_CUSTOMCOMPOUNDBONDFORCE_GETFUNCTIONPARAMETERS(const OpenMM_CustomCompoundBondForce*& target, int const& index, char** name, OpenMM_DoubleArray*& values, double* min, double* max) {
    OpenMM_CustomCompoundBondForce_getFunctionParameters(target, index, name, values, min, max);
}
OPENMM_EXPORT void openmm_customcompoundbondforce_setfunctionparameters_(OpenMM_CustomCompoundBondForce*& target, int const& index, const char* name, const OpenMM_DoubleArray*& values, double const& min, double const& max, int name_length) {
    OpenMM_CustomCompoundBondForce_setFunctionParameters(target, index, makeString(name, name_length).c_str(), values, min, max);
}
OPENMM_EXPORT void OPENMM_CUSTOMCOMPOUNDBONDFORCE_SETFUNCTIONPARAMETERS(OpenMM_CustomCompoundBondForce*& target, int const& index, const char* name, const OpenMM_DoubleArray*& values, double const& min, double const& max, int name_length) {
    OpenMM_CustomCompoundBondForce_setFunctionParameters(target, index, makeString(name, name_length).c_str(), values, min, max);
}
OPENMM_EXPORT void openmm_customcompoundbondforce_updateparametersincontext_(OpenMM_CustomCompoundBondForce*& target, OpenMM_Context*& context) {
    OpenMM_CustomCompoundBondForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void OPENMM_CUSTOMCOMPOUNDBONDFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_CustomCompoundBondForce*& target, OpenMM_Context*& context) {
    OpenMM_CustomCompoundBondForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void openmm_customcompoundbondforce_setusesperiodicboundarycondition_(OpenMM_CustomCompoundBondForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_CustomCompoundBondForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT void OPENMM_CUSTOMCOMPOUNDBONDFORCE_SETUSESPERIODICBOUNDARYCONDITION(OpenMM_CustomCompoundBondForce*& target, OpenMM_Boolean& periodic) {
    OpenMM_CustomCompoundBondForce_setUsesPeriodicBoundaryConditions(target, periodic);
}
OPENMM_EXPORT void openmm_customcompoundbondforce_usesperiodicboundaryconditions_(const OpenMM_CustomCompoundBondForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CustomCompoundBondForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_CUSTOMCOMPOUNDBONDFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_CustomCompoundBondForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_CustomCompoundBondForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::RMSDForce */
OPENMM_EXPORT void openmm_rmsdforce_create_(OpenMM_RMSDForce*& result, const OpenMM_Vec3Array*& referencePositions, const OpenMM_IntArray*& particles) {
    result = OpenMM_RMSDForce_create(referencePositions, particles);
}
OPENMM_EXPORT void OPENMM_RMSDFORCE_CREATE(OpenMM_RMSDForce*& result, const OpenMM_Vec3Array*& referencePositions, const OpenMM_IntArray*& particles) {
    result = OpenMM_RMSDForce_create(referencePositions, particles);
}
OPENMM_EXPORT void openmm_rmsdforce_destroy_(OpenMM_RMSDForce*& destroy) {
    OpenMM_RMSDForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_RMSDFORCE_DESTROY(OpenMM_RMSDForce*& destroy) {
    OpenMM_RMSDForce_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void openmm_rmsdforce_getreferencepositions_(const OpenMM_RMSDForce*& target, const OpenMM_Vec3Array*& result) {
    result = OpenMM_RMSDForce_getReferencePositions(target);
}
OPENMM_EXPORT void OPENMM_RMSDFORCE_GETREFERENCEPOSITIONS(const OpenMM_RMSDForce*& target, const OpenMM_Vec3Array*& result) {
    result = OpenMM_RMSDForce_getReferencePositions(target);
}
OPENMM_EXPORT void openmm_rmsdforce_setreferencepositions_(OpenMM_RMSDForce*& target, const OpenMM_Vec3Array*& positions) {
    OpenMM_RMSDForce_setReferencePositions(target, positions);
}
OPENMM_EXPORT void OPENMM_RMSDFORCE_SETREFERENCEPOSITIONS(OpenMM_RMSDForce*& target, const OpenMM_Vec3Array*& positions) {
    OpenMM_RMSDForce_setReferencePositions(target, positions);
}
OPENMM_EXPORT void openmm_rmsdforce_getparticles_(const OpenMM_RMSDForce*& target, const OpenMM_IntArray*& result) {
    result = OpenMM_RMSDForce_getParticles(target);
}
OPENMM_EXPORT void OPENMM_RMSDFORCE_GETPARTICLES(const OpenMM_RMSDForce*& target, const OpenMM_IntArray*& result) {
    result = OpenMM_RMSDForce_getParticles(target);
}
OPENMM_EXPORT void openmm_rmsdforce_setparticles_(OpenMM_RMSDForce*& target, const OpenMM_IntArray*& particles) {
    OpenMM_RMSDForce_setParticles(target, particles);
}
OPENMM_EXPORT void OPENMM_RMSDFORCE_SETPARTICLES(OpenMM_RMSDForce*& target, const OpenMM_IntArray*& particles) {
    OpenMM_RMSDForce_setParticles(target, particles);
}
OPENMM_EXPORT void openmm_rmsdforce_updateparametersincontext_(OpenMM_RMSDForce*& target, OpenMM_Context*& context) {
    OpenMM_RMSDForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void OPENMM_RMSDFORCE_UPDATEPARAMETERSINCONTEXT(OpenMM_RMSDForce*& target, OpenMM_Context*& context) {
    OpenMM_RMSDForce_updateParametersInContext(target, context);
}
OPENMM_EXPORT void openmm_rmsdforce_usesperiodicboundaryconditions_(const OpenMM_RMSDForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_RMSDForce_usesPeriodicBoundaryConditions(target);
}
OPENMM_EXPORT void OPENMM_RMSDFORCE_USESPERIODICBOUNDARYCONDITIONS(const OpenMM_RMSDForce*& target, OpenMM_Boolean& result) {
    result = OpenMM_RMSDForce_usesPeriodicBoundaryConditions(target);
}

/* OpenMM::BrownianIntegrator */
OPENMM_EXPORT void openmm_brownianintegrator_create_(OpenMM_BrownianIntegrator*& result, double const& temperature, double const& frictionCoeff, double const& stepSize) {
    result = OpenMM_BrownianIntegrator_create(temperature, frictionCoeff, stepSize);
}
OPENMM_EXPORT void OPENMM_BROWNIANINTEGRATOR_CREATE(OpenMM_BrownianIntegrator*& result, double const& temperature, double const& frictionCoeff, double const& stepSize) {
    result = OpenMM_BrownianIntegrator_create(temperature, frictionCoeff, stepSize);
}
OPENMM_EXPORT void openmm_brownianintegrator_destroy_(OpenMM_BrownianIntegrator*& destroy) {
    OpenMM_BrownianIntegrator_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_BROWNIANINTEGRATOR_DESTROY(OpenMM_BrownianIntegrator*& destroy) {
    OpenMM_BrownianIntegrator_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT double openmm_brownianintegrator_gettemperature_(const OpenMM_BrownianIntegrator*& target) {
    return OpenMM_BrownianIntegrator_getTemperature(target);
}
OPENMM_EXPORT double OPENMM_BROWNIANINTEGRATOR_GETTEMPERATURE(const OpenMM_BrownianIntegrator*& target) {
    return OpenMM_BrownianIntegrator_getTemperature(target);
}
OPENMM_EXPORT void openmm_brownianintegrator_settemperature_(OpenMM_BrownianIntegrator*& target, double const& temp) {
    OpenMM_BrownianIntegrator_setTemperature(target, temp);
}
OPENMM_EXPORT void OPENMM_BROWNIANINTEGRATOR_SETTEMPERATURE(OpenMM_BrownianIntegrator*& target, double const& temp) {
    OpenMM_BrownianIntegrator_setTemperature(target, temp);
}
OPENMM_EXPORT double openmm_brownianintegrator_getfriction_(const OpenMM_BrownianIntegrator*& target) {
    return OpenMM_BrownianIntegrator_getFriction(target);
}
OPENMM_EXPORT double OPENMM_BROWNIANINTEGRATOR_GETFRICTION(const OpenMM_BrownianIntegrator*& target) {
    return OpenMM_BrownianIntegrator_getFriction(target);
}
OPENMM_EXPORT void openmm_brownianintegrator_setfriction_(OpenMM_BrownianIntegrator*& target, double const& coeff) {
    OpenMM_BrownianIntegrator_setFriction(target, coeff);
}
OPENMM_EXPORT void OPENMM_BROWNIANINTEGRATOR_SETFRICTION(OpenMM_BrownianIntegrator*& target, double const& coeff) {
    OpenMM_BrownianIntegrator_setFriction(target, coeff);
}
OPENMM_EXPORT int openmm_brownianintegrator_getrandomnumberseed_(const OpenMM_BrownianIntegrator*& target) {
    return OpenMM_BrownianIntegrator_getRandomNumberSeed(target);
}
OPENMM_EXPORT int OPENMM_BROWNIANINTEGRATOR_GETRANDOMNUMBERSEED(const OpenMM_BrownianIntegrator*& target) {
    return OpenMM_BrownianIntegrator_getRandomNumberSeed(target);
}
OPENMM_EXPORT void openmm_brownianintegrator_setrandomnumberseed_(OpenMM_BrownianIntegrator*& target, int const& seed) {
    OpenMM_BrownianIntegrator_setRandomNumberSeed(target, seed);
}
OPENMM_EXPORT void OPENMM_BROWNIANINTEGRATOR_SETRANDOMNUMBERSEED(OpenMM_BrownianIntegrator*& target, int const& seed) {
    OpenMM_BrownianIntegrator_setRandomNumberSeed(target, seed);
}
OPENMM_EXPORT void openmm_brownianintegrator_step_(OpenMM_BrownianIntegrator*& target, int const& steps) {
    OpenMM_BrownianIntegrator_step(target, steps);
}
OPENMM_EXPORT void OPENMM_BROWNIANINTEGRATOR_STEP(OpenMM_BrownianIntegrator*& target, int const& steps) {
    OpenMM_BrownianIntegrator_step(target, steps);
}

/* OpenMM::VariableLangevinIntegrator */
OPENMM_EXPORT void openmm_variablelangevinintegrator_create_(OpenMM_VariableLangevinIntegrator*& result, double const& temperature, double const& frictionCoeff, double const& errorTol) {
    result = OpenMM_VariableLangevinIntegrator_create(temperature, frictionCoeff, errorTol);
}
OPENMM_EXPORT void OPENMM_VARIABLELANGEVININTEGRATOR_CREATE(OpenMM_VariableLangevinIntegrator*& result, double const& temperature, double const& frictionCoeff, double const& errorTol) {
    result = OpenMM_VariableLangevinIntegrator_create(temperature, frictionCoeff, errorTol);
}
OPENMM_EXPORT void openmm_variablelangevinintegrator_destroy_(OpenMM_VariableLangevinIntegrator*& destroy) {
    OpenMM_VariableLangevinIntegrator_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT void OPENMM_VARIABLELANGEVININTEGRATOR_DESTROY(OpenMM_VariableLangevinIntegrator*& destroy) {
    OpenMM_VariableLangevinIntegrator_destroy(destroy);
    destroy = 0;
}
OPENMM_EXPORT double openmm_variablelangevinintegrator_gettemperature_(const OpenMM_VariableLangevinIntegrator*& target) {
    return OpenMM_VariableLangevinIntegrator_getTemperature(target);
}
OPENMM_EXPORT double OPENMM_VARIABLELANGEVININTEGRATOR_GETTEMPERATURE(const OpenMM_VariableLangevinIntegrator*& target) {
    return OpenMM_VariableLangevinIntegrator_getTemperature(target);
}
OPENMM_EXPORT void openmm_variablelangevinintegrator_settemperature_(OpenMM_VariableLangevinIntegrator*& target, double const& temp) {
    OpenMM_VariableLangevinIntegrator_setTemperature(target, temp);
}
OPENMM_EXPORT void OPENMM_VARIABLELANGEVININTEGRATOR_SETTEMPERATURE(OpenMM_VariableLangevinIntegrator*& target, double const& temp) {
    OpenMM_VariableLangevinIntegrator_setTemperature(target, temp);
}
OPENMM_EXPORT double openmm_variablelangevinintegrator_getfriction_(const OpenMM_VariableLangevinIntegrator*& target) {
    return OpenMM_VariableLangevinIntegrator_getFriction(target);
}
OPENMM_EXPORT double OPENMM_VARIABLELANGEVININTEGRATOR_GETFRICTION(const OpenMM_VariableLangevinIntegrator*& target) {
    return OpenMM_VariableLangevinIntegrator_getFriction(target);
}
OPENMM_EXPORT void openmm_variablelangevinintegrator_setfriction_(OpenMM_VariableLangevinIntegrator*& target, double const& coeff) {
    OpenMM_VariableLangevinIntegrator_setFriction(target, coeff);
}
OPENMM_EXPORT void OPENMM_VARIABLELANGEVININTEGRATOR_SETFRICTION(OpenMM_VariableLangevinIntegrator*& target, double const& coeff) {
    OpenMM_VariableLangevinIntegrator_setFriction(target, coeff);
}
OPENMM_EXPORT double openmm_variablelangevinintegrator_geterrortolerance_(const OpenMM_VariableLangevinIntegrator*& target) {
    return OpenMM_VariableLangevinIntegrator_getErrorTolerance(target);
}
OPENMM_EXPORT double OPENMM_VARIABLELANGEVININTEGRATOR_GETERRORTOLERANCE(const OpenMM_VariableLangevinIntegrator*& target) {
    return OpenMM_VariableLangevinIntegrator_getErrorTolerance(target);
}
OPENMM_EXPORT void openmm_variablelangevinintegrator_seterrortolerance_(OpenMM_VariableLangevinIntegrator*& target, double const& tol) {
    OpenMM_VariableLangevinIntegrator_setErrorTolerance(target, tol);
}
OPENMM_EXPORT void OPENMM_VARIABLELANGEVININTEGRATOR_SETERRORTOLERANCE(OpenMM_VariableLangevinIntegrator*& target, double const& tol) {
    OpenMM_VariableLangevinIntegrator_setErrorTolerance(target, tol);
}
OPENMM_EXPORT double openmm_variablelangevinintegrator_getmaximumstepsize_(const OpenMM_VariableLangevinIntegrator*& target) {
    return OpenMM_VariableLangevinIntegrator_getMaximumStepSize(target);
}
OPENMM_EXPORT double OPENMM_VARIABLELANGEVININTEGRATOR_GETMAXIMUMSTEPSIZE(const OpenMM_VariableLangevinIntegrator*& target) {
    return OpenMM_VariableLangevinIntegrator_getMaximumStepSize(target);
}
OPENMM_EXPORT void openmm_variablelangevinintegrator_setmaximumstepsize_(OpenMM_VariableLangevinIntegrator*& target, double const& size) {
    OpenMM_VariableLangevinIntegrator_setMaximumStepSize(target, size);
}
OPENMM_EXPORT void OPENMM_VARIABLELANGEVININTEGRATOR_SETMAXIMUMSTEPSIZE(OpenMM_VariableLangevinIntegrator*& target, double const& size) {
    OpenMM_VariableLangevinIntegrator_setMaximumStepSize(target, size);
}
OPENMM_EXPORT int openmm_variablelangevinintegrator_getrandomnumberseed_(const OpenMM_VariableLangevinIntegrator*& target) {
    return OpenMM_VariableLangevinIntegrator_getRandomNumberSeed(target);
}
OPENMM_EXPORT int OPENMM_VARIABLELANGEVININTEGRATOR_GETRANDOMNUMBERSEED(const OpenMM_VariableLangevinIntegrator*& target) {
    return OpenMM_VariableLangevinIntegrator_getRandomNumberSeed(target);
}
OPENMM_EXPORT void openmm_variablelangevinintegrator_setrandomnumberseed_(OpenMM_VariableLangevinIntegrator*& target, int const& seed) {
    OpenMM_VariableLangevinIntegrator_setRandomNumberSeed(target, seed);
}
OPENMM_EXPORT void OPENMM_VARIABLELANGEVININTEGRATOR_SETRANDOMNUMBERSEED(OpenMM_VariableLangevinIntegrator*& target, int const& seed) {
    OpenMM_VariableLangevinIntegrator_setRandomNumberSeed(target, seed);
}
OPENMM_EXPORT void openmm_variablelangevinintegrator_step_(OpenMM_VariableLangevinIntegrator*& target, int const& steps) {
    OpenMM_VariableLangevinIntegrator_step(target, steps);
}
OPENMM_EXPORT void OPENMM_VARIABLELANGEVININTEGRATOR_STEP(OpenMM_VariableLangevinIntegrator*& target, int const& steps) {
    OpenMM_VariableLangevinIntegrator_step(target, steps);
}
OPENMM_EXPORT void openmm_variablelangevinintegrator_stepto_(OpenMM_VariableLangevinIntegrator*& target, double const& time) {
    OpenMM_VariableLangevinIntegrator_stepTo(target, time);
}
OPENMM_EXPORT void OPENMM_VARIABLELANGEVININTEGRATOR_STEPTO(OpenMM_VariableLangevinIntegrator*& target, double const& time) {
    OpenMM_VariableLangevinIntegrator_stepTo(target, time);
}

}
