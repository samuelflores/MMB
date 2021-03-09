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

//class MMB_EXPORT MobilizerContainer : public ResidueStretchContainer <MobilizerStretch> 
class MMB_EXPORT DensityContainer : public ResidueStretchContainer <DensityStretch> {
private :
	//std::vector <DensityStretch> residueStretchVector; // should override 
public:
	void clear();
	void validate(const DensityStretch & myDensityStretch,BiopolymerClassContainer & myBiopolymerClassContainer);
	void	add(const DensityStretch & myDensityStretch,BiopolymerClassContainer & myBiopolymerClassContainer);
	DensityStretch getDensityStretch(int densityStretchIndex);
	int numDensityStretches();
	void stuffDensityStretchVector(BiopolymerClassContainer & myBiopolymerClassContainer );
    void printDensityStretches();

    std::vector <DensityStretch> getDensityStretchVector() { return residueStretchVector; }

    void updateDensityStretch(int id, const DensityStretch & stretch, BiopolymerClassContainer & myBiopolymerClassContainer);
    void deleteDensityStretch(int id);
};
#endif

