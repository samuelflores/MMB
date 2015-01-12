/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#ifndef DisplacementContainer_H_
#define DisplacementContainer_H_

#include "BiopolymerClass.h"
#include "Utils.h"

class DisplacementContainer {
private:
        vector <Displacement> displacementVector;
public:
        void clear() {displacementVector.clear();};
        DisplacementContainer() {clear(); };
        void validateDisplacement(const Displacement displacement ,   BiopolymerClassContainer & myBiopolymerContainer);
        void initializeDisplacement(Displacement & displacement) {displacement.chain = ""; displacement.displacement = Vec3(0);};
        void add(Displacement displacement,BiopolymerClassContainer & myBiopolymerContainer ){validateDisplacement(displacement,myBiopolymerContainer); displacementVector.push_back(displacement); };
        Displacement & updDisplacement(int displacementIndex) {return displacementVector[displacementIndex]; };
        Displacement & updDisplacement(String chain) ; 
        int numDisplacements() {return displacementVector.size(); };
        bool hasChain(String chain);
        Vec3 getInitialDisplacementVec3(String chain) ;
        vector <Displacement> getDisplacementVector( ){return displacementVector; };
};




#endif

