/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Molmodel                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "molmodel/internal/PDBReader.h"
#include "molmodel/internal/Protein.h"
#include "molmodel/internal/RNA.h"
#include "mol.h"
#include "SimTKsimbody.h"
#include <cctype>
#include <map>
#include <vector>
#include <algorithm>

using namespace SimTK;
using std::map;
using std::string;
using std::vector;

class PDBReader::PDBReaderImpl {
public:
    PDBReaderImpl(string filename) : hasBuiltSystem(false) {
        mol_DbRead("mol", filename.c_str(), MOL_DB_PDB, &model);
        MolStructure *structure;
        int numStructures;
        mol_MolModelStructuresGet(model, &numStructures, &structure);
        SimTK_APIARGCHECK_ALWAYS(numStructures == 1, "PDBReaderImpl", "PDBReaderImpl", "The PDB file must contain exactly one structure");
    }

    ~PDBReaderImpl() {
        mol_MemFreeModel(model);
    }

    void createCompounds(CompoundSystem& system) {
        SimTK_APIARGCHECK_ALWAYS(!hasBuiltSystem, "PDBReaderImpl", "createSystem", "createSystem() has already been called");
        MolStructure *structure;
        int numStructures;
        mol_MolModelStructuresGet(model, &numStructures, &structure);
        MolAtom* atoms;
        int numAtoms;
        mol_StructureAtomsGet(structure, &numAtoms, &atoms);
        MolChain* chains;
        int numChains;
        mol_StructureChainsGet(structure, &numChains, &chains);
        
        // Loop over chains and create a Biopolymer from each one.
        
        while (chains) {
            MolResidue** chainResidues;
            int numResidues;
            mol_ChainResiduesGet(chains, &numResidues, &chainResidues);
            
            // Create a string of the sequence.
            
            string sequence;
            for (int i = 0; i < numResidues; ++i) {
                sequence += mol_res_names[chainResidues[i]->type][2];
            }
            std::transform(sequence.begin(), sequence.end(), sequence.begin(), (int(*)(int)) std::toupper);
            if (chainResidues[0]->type > 24) {
                
                // Create an RNA.
                
                RNA rna(sequence);
                for (int i = 0; i < rna.getNumResidues(); i ++) {
                    rna.updResidue(ResidueInfo::Index(i)).setPdbResidueNumber(chainResidues[i]->id);
                    rna.updResidue(ResidueInfo::Index(i)).setPdbInsertionCode(chainResidues[i]->insertion_code);//    ('101A');
                }
                rna.assignBiotypes();
                compounds.push_back(rna);
            }
            else if (chainResidues[0]->type < 21) {
                
                // Create a Protein.
                // scf changed this to use one more parameter in the Protein constructor, set to "Torsion".  Default is "Rigid".  Now the peptide bond will not be rigid.
                Protein protein(sequence,BondMobility::Torsion);

                    protein.updResidue(ResidueInfo::Index(0)).setPdbResidueNumber(chainResidues[0]->id-1);
                    protein.updResidue(ResidueInfo::Index(0)).setPdbInsertionCode(' ');
                    std::cout<<__FILE__<<":"<<__LINE__<<" i, residue number, insertion code, residue type : "<< protein.updResidue(ResidueInfo::Index(0)).getPdbResidueNumber()  <<", "<< protein.updResidue(ResidueInfo::Index(0)).getPdbInsertionCode()<<", "<<protein.updResidue(ResidueInfo::Index(0)).getOneLetterCode()<<std::endl;
                    //std::cout<<__FILE__<<":"<<__LINE__<<" i, residue number, insertion code, residue type : "<<chainResidues[0]->id<<", "<<chainResidues[0]->insertion_code<<", "<< protein.updResidue(ResidueInfo::Index(0)).getOneLetterCode()<<std::endl;
                for (int i = 1; i < (protein.getNumResidues() - 1); i ++) { // assume protein capping is ON.  this means first and last residues are end caps, we ignore at this stage.
                    protein.updResidue(ResidueInfo::Index(i)).setPdbResidueNumber(chainResidues[i-1]->id);
                    protein.updResidue(ResidueInfo::Index(i)).setPdbInsertionCode(chainResidues[i-1]->insertion_code);
                    std::cout<<__FILE__<<":"<<__LINE__<<" i, residue number, insertion code, residue type : "<< protein.updResidue(ResidueInfo::Index(i)).getPdbResidueNumber()  <<", "<< protein.updResidue(ResidueInfo::Index(i)).getPdbInsertionCode()<<", "<<protein.updResidue(ResidueInfo::Index(i)).getOneLetterCode()<<std::endl;
                }
	            // was unable to retrieve C-terminal end cap for some reason:	 	
                    //protein.updResidue(ResidueInfo::Index((protein.getNumResidues() - 1) )).setPdbResidueNumber(chainResidues[(protein.getNumResidues() - 1) ]->id+1);
                    //protein.updResidue(ResidueInfo::Index((protein.getNumResidues() - 1) )).setPdbInsertionCode(' ');
                    std::cout<<__FILE__<<":"<<__LINE__<<" i, residue number, insertion code, residue type : "<< protein.updResidue(ResidueInfo::Index((protein.getNumResidues() - 1) )).getPdbResidueNumber()  <<", "<< protein.updResidue(ResidueInfo::Index( (protein.getNumResidues() - 1) )).getPdbInsertionCode()<<", "<<protein.updResidue(ResidueInfo::Index( (protein.getNumResidues() - 1) )).getOneLetterCode()<<std::endl;
                protein.assignBiotypes();
                compounds.push_back(protein);
            }
            compounds[compounds.size()-1].setPdbChainId((*chains).id );
            chains = chains->next;
        }
        
        // Add them to the system.
        
        for (int i = 0; i < (int)compounds.size(); ++i)
            system.adoptCompound(compounds[i]);
        hasBuiltSystem = true;
    }
    
    Real createState(const CompoundSystem& system, State& state) const {
        SimTK_APIARGCHECK_ALWAYS(hasBuiltSystem, "PDBReaderImpl", "createState", "createSystem() has not yet been called");

        MolStructure *structure;
        int numStructures;
        mol_MolModelStructuresGet(model, &numStructures, &structure);
        MolAtom* atoms;
        int numAtoms;
        mol_StructureAtomsGet(structure, &numAtoms, &atoms);
        MolChain* chains;
        int numChains;
        mol_StructureChainsGet(structure, &numChains, &chains);

        // Loop over atoms, match each one to the appropriate mobilized body, and
        // create a list of stations that will be used for fitting the State.

        map<MobilizedBodyIndex, vector<Vec3> > stations;
        map<MobilizedBodyIndex, vector<Vec3> > targetLocations;
        int proteinIndex = 0;
        while (chains) {
            MolResidue** chainResidues;
            int numResidues;
            mol_ChainResiduesGet(chains, &numResidues, &chainResidues);
            for (int res = 0; res < numResidues; ++res) {
                MolAtomList resAtoms = chainResidues[res]->atoms;
                char resName[32];
                sprintf(resName, "%d", res);
                const Biopolymer& compound = compounds[proteinIndex];
                const ResidueInfo::Index resIx = compound.getResidue(resName).getIndex();
                const ResidueInfo& residue = compound.getResidue(resName);
// residue.setBondMobility(
                for (int i = 0; i < resAtoms.num; ++i) {
                    MolAtom atom = atoms[resAtoms.list[i]];
                    string atomName = atom.name;
                    atomName = atomName.erase(atomName.find_last_not_of(' ')+1);
                    atomName.erase(0, atomName.find_first_not_of(' '));
                    for (ResidueInfo::AtomIndex atomId(0); atomId < residue.getNumAtoms(); ++atomId) {
                        if (atomName == residue.getAtomName(atomId)) {
                            MobilizedBodyIndex bodyId = compound.getResidueAtomMobilizedBodyIndex(resIx, atomId);
                            stations[bodyId].push_back(compound.getResidueAtomLocationInMobilizedBodyFrame(resIx, atomId));
                            targetLocations[bodyId].push_back(Vec3(atom.pos[0], atom.pos[1], atom.pos[2]));
                        }
                    }
                }
            }
            proteinIndex++;
            chains = chains->next;
        }

        // Now perform the fitting.
        
        vector<MobilizedBodyIndex> bodyList;
        vector<vector<Vec3> > stationList;
        vector<vector<Vec3> > targetList;
        for (map<MobilizedBodyIndex, vector<Vec3> >::const_iterator iter = stations.begin(); iter != stations.end(); iter++) {
            bodyList.push_back(iter->first);
            stationList.push_back(iter->second);
            targetList.push_back(targetLocations.find(iter->first)->second);
        }
        // sherm 100307: Optimizers now use relative tolerance.
        Real tolerance = .001; // 0.1%
        return ObservedPointFitter::findBestFit(system, state, bodyList, stationList, targetList, tolerance);
    }
        
private:
    MolModel *model;
    vector<Biopolymer> compounds;
    bool hasBuiltSystem;
};

PDBReader::PDBReader(string filename) {
    impl = new PDBReaderImpl(filename);
}

PDBReader::~PDBReader() {
    delete impl;
}

void PDBReader::createCompounds(CompoundSystem& system) {
    impl->createCompounds(system);
}

Real PDBReader::createState(const CompoundSystem& system, State& state) const {
    return impl->createState(system, state);
}

