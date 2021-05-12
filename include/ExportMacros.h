/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-13 by the Authors.                                      *
 * Authors: Samuel Flores, Alex Tek                                           *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#ifndef ExportMacros_H_
#define ExportMacros_H_  

/*
 * Shared libraries are messy in Visual Studio. We have to distinguish three
 * cases:
 *   (1) this header is being used to build the MMB shared library (dllexport)
 *   (2) this header is being used by a *client* of the MMB shared
 *       library (dllimport)
 *   (3) we are building the MMB static library, or the client is
 *       being compiled with the expectation of linking with the
 *       MMB static library (nothing special needed)
 * In the CMake script for building this library, we define one of the symbols
 *     MMB_BUILDING_{SHARED|STATIC}_LIBRARY
 * Client code normally has no special symbol defined, in which case we'll
 * assume it wants to use the shared library. However, if the client defines
 * the symbol MMB_USE_STATIC_LIBRARIES we'll suppress the dllimport so
 * that the client code can be linked with static libraries. If this flag
 * is set, we'll also define SimTK_USE_STATIC_LIBRARIES.
 */

#ifdef _WIN32
    #ifdef _MSC_VER
    #pragma warning(disable:4231) /*need to use 'extern' template explicit instantiation*/
    #pragma warning(disable:4251) /*no DLL interface for type of member of exported class*/
    #pragma warning(disable:4275) /*no DLL interface for base class of exported class*/
    #pragma warning(disable:4345) /*warning about PODs being default-initialized*/
    #endif
    #if defined(MMB_BUILDING_SHARED_LIBRARY)
        #define MMB_EXPORT __declspec(dllexport)
        /* Keep MS VC++ quiet when it tries to instantiate incomplete template classes in a DLL. */
        #ifdef _MSC_VER
        #pragma warning(disable:4661)
        #endif
    #elif defined(MMB_BUILDING_STATIC_LIBRARY) || defined(MMB_USE_STATIC_LIBRARIES)
        #define MMB_EXPORT
        #ifndef SimTK_USE_STATIC_LIBRARIES
            #define SimTK_USE_STATIC_LIBRARIES
        #endif
    #else
        #define MMB_EXPORT __declspec(dllimport) /*i.e., a client of a shared library*/
    #endif
#else
    #if defined(MMB_BUILDING_SHARED_LIBRARY)
	#define MMB_EXPORT __attribute__ ((visibility ("default")))
    #else
        #define MMB_EXPORT // Linux, Mac
    #endif
#endif

#endif
