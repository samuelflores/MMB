/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#ifndef DensityContainer_H_
#define DensityContainer_H_

#include "BiopolymerClass.h"
#include "ResidueStretchContainer.h"

class DensityContainer : public ResidueStretchContainer <DensityStretch>{
private :
    std::vector <DensityStretch> & densityStretchVector = ResidueStretchContainer::residueStretchVector; // alias -- use parent's residueStretchVector member.Since this is private, alias may be unnecessary..
public:
    void clear();
    //void validate(const DensityStretch & myDensityStretch,BiopolymerClassContainer & myBiopolymerClassContainer); // Already defined in parent.
    //void	add(const DensityStretch & myDensityStretch,BiopolymerClassContainer & myBiopolymerClassContainer){addStretch(}; // just use parent's addStretch(ResidueStretchType newStretch) 
    //void	add(const DensityStretch & myDensityStretch,BiopolymerClassContainer & myBiopolymerClassContainer){addStretch(}; // just use parent's addStretch(ResidueStretchType newStretch) 
    void addStretch(const ResidueStretch myResidueStretch){addStretch(DensityStretch(myResidueStretch));}; // If user provides a ResidueStretch as argument, upcast ResidueStretch as a DensityStretch, to use the addStretch member already defined in the parent.
    DensityStretch getDensityStretch(int densityStretchIndex){return getResidueStretch(densityStretchIndex);}; // just use parent's getResidueStretch
    int numDensityStretches(){return getNumResidueStretches();}; // use parent member
    void stuffDensityStretchVector(BiopolymerClassContainer & myBiopolymerClassContainer );
    void printDensityStretches(){printResidueStretchVector();}; // just use printResidueStretchVector()
    std::vector <DensityStretch> getDensityStretchVector() { return getResidueStretchVector();} //residueStretchVector;}//densityStretchVector; } // just use parent's getResidueStretchVector
    void updateDensityStretch(int id, const DensityStretch & stretch, BiopolymerClassContainer & myBiopolymerClassContainer);
    void deleteDensityStretch(int id){deleteResidueStretch(id);}; // just use parent's deleteResidueStretch(int id)
};
#endif

