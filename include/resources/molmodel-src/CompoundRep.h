#ifndef SimTK_COMPOUNDREP_H_
#define SimTK_COMPOUNDREP_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Molmodel                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
 * Authors: Christopher Bruns                                                 *
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

#include "SimTKcommon.h"

#include "molmodel/internal/bondGeometry.h"
#include "molmodel/internal/Compound.h"
#include "molmodel/internal/CompoundSystem.h"
#include "molmodel/internal/Superpose.h"

#include "CompoundAtom.h"

#include <vector>
#include <map>
#include <string>
#include <set>
using std::string;

namespace SimTK {

static const String InboardBondName = "inboard bond";


void invalidateAtomFrameCache(std::vector<Transform>& atomFrameCache, int numAtoms);


// Manager class for uniform handling of bonded, locally placed, and deeper subcompounds.
// These are indexed by Compound::Index.
//class CompoundInfo {
//public:
//    // Constructor for directly bonded subcompounds
//    // TODO - two other constructors for nonbonded and distant descendants
//    CompoundInfo(Compound::Index compoundId, const Compound::Name& scName, const Compound& subcompound, Compound::BondIndex bondIndex)
//        : local(false), bonded(true), id(compoundId), name(scName), compoundBondIndex(bondIndex), compound(subcompound)
//    {}
//
//    // Constructor for local non-bonded compounds
//    CompoundInfo(Compound::Index compoundId, const Compound::Name& scName, const Compound& subcompound, const Transform& transform)
//        : local(true), bonded(false), id(compoundId), name(scName), 
//        // localSubcompoundId(localId), 
//        localFrameInParentFrame(transform),
//        compound(subcompound)
//    {}
//
//    // Constructor for indirectly contained subsubcompounds
//    CompoundInfo(Compound::Index compoundId, Compound::Index topSubcompoundId, Compound::Index childSubcompoundId)
//        : local(false), bonded(false), id(compoundId)
//        , intermediateSubcompoundId(topSubcompoundId)
//        , intermediateSubcompoundSubcompoundId(childSubcompoundId)
//    {
//    }
//
//
//    Compound&       updCompound()       {return compound;}
//    const Compound& getCompound() const {return compound;}
//
//    bool isLocal()  const {return local;}
//    bool isBonded() const {return bonded;}
//
//    CompoundInfo& setBondIndex(Compound::BondIndex bId) {
//        assert(isBonded());
//        assert(!compoundBondIndex.isValid());
//
//        compoundBondIndex = bId;
//
//        assert(compoundBondIndex.isValid());
//
//        return *this;
//    }
//
//    Compound::Index getIndex() const {
//        return id;
//    }
//
//    //int getLocalSubcompoundId() const {
//    //    assert( isLocal() );
//    //    return localSubcompoundId;
//    //}
//
//    const Compound::Name& getName() const {return name;}
//
//    const Transform& getLocalFrameInParentFrame() const {
//        assert( isLocal() );
//        assert( ! isBonded() );
//        return localFrameInParentFrame;
//    }
//
//    Compound::BondIndex getCompoundBondIndex() const {
//        assert( ! isLocal() );
//        assert( isBonded() );
//        return compoundBondIndex;
//    }
//
//    Compound::Index getIntermediateSubcompoundId() const {
//        assert( ! isLocal() );
//        assert( ! isBonded() );
//        return intermediateSubcompoundId;
//    }
//    Compound::Index getIntermediateSubcompoundSubcompoundId() const {
//        assert( ! isLocal() );
//        assert( ! isBonded() );
//        return intermediateSubcompoundSubcompoundId;
//    }
//
//private:
//    bool local;  // locally placed in this compound
//    bool bonded; // attached to an atom via a bond
//
//    Compound::Index id;
//
//    Compound compound; // shared by cases local and bonded
//
//    // 1) if is local
//    // int localSubcompoundId;
//    Transform localFrameInParentFrame;
//    Compound::Name name;
//
//    // 2) if is bonded
//    Compound::BondIndex compoundBondIndex;
//
//    // 3) if belongs to a subsubcompound
//    Compound::Index intermediateSubcompoundId;
//    Compound::Index intermediateSubcompoundSubcompoundId;
//};


class DihedralAngle {
public:
    DihedralAngle() 
        : nomenclatureOffset(0*Deg2Rad)
    {}

    DihedralAngle(Compound::BondCenterIndex bc1, Compound::BondCenterIndex bc2, Angle o = 0*Deg2Rad) 
        : bondCenter1(bc1), bondCenter2(bc2), nomenclatureOffset(o)
        // , internalOffset(0)
    {}

    Compound::BondCenterIndex getBondCenter1Id() const {return bondCenter1;}
    Compound::BondCenterIndex getBondCenter2Id() const {return bondCenter2;}
    Angle getNomenclatureOffset() const {return nomenclatureOffset;}

    // Angle getInternalOffset() const {return internalOffset;}
    // void setInternalOffset(Angle a) {internalOffset = a;}

private:
    // Angle internalOffset; // this dihedral may be offset from canonical dihedral for this bond
    Angle nomenclatureOffset; // internal + offset = nominal
    Compound::BondCenterIndex bondCenter1;
    Compound::BondCenterIndex bondCenter2;
};

    //////////////////
    // COMPOUND REP //
    //////////////////

class CompoundRep : public PIMPLImplementation<Compound,CompoundRep> 
{
public:
    friend class CompoundSystem;
    friend class ResidueInfo;

    explicit CompoundRep(const String& n="UnknownCompoundType", const Transform& transform = Transform()) 
      : ownerSystem(0), 
      name(n), 
      pdbChainId(' '), 
      pdbResidueName("UNK"), 
      pdbResidueNumber(-111111), 
      // haveParentCompound(false),
      topLevelTransform(transform)
    {
    	addCompoundSynonym(n);
    }

    CompoundRep&  setCompoundName(const Compound::Name& n) {
        name=n; 
        addCompoundSynonym(name);
        return *this;
    }
    CompoundRep& addCompoundSynonym(Compound::Name synonym) {
        synonyms.insert(synonym);
        return *this;
    }
    const Compound::Name& getCompoundName() const          {return name;}

    // This is concrete, but can be extended.
    virtual ~CompoundRep() { }
    virtual CompoundRep* clone() const {return new CompoundRep(*this);}

    // int getNumSubcompounds() const {return allSubcompounds.size();}

    void setMultibodySystem(MultibodySystem& system) 
    {
        ownerSystem = &system;

        // Also give subcompounds access to the system, but lacking a valid index
        //for (Compound::Index c(0); c < getNumSubcompounds(); ++c)
        //{
        //    CompoundRep& scRep = updSubcompound(c).updImpl();
        //    scRep.setMultibodySystem(system);
        //}
    }
    //void setCompoundSystem(CompoundSystem& system, Compound::Index id) 
    //{
    //    ownerSystem = &system;
    //    ixWithinOwnerSystem = id;

    //    // Also give subcompounds access to the system, but lacking a valid index
    //    for (Compound::Index c(0); c < getNumSubcompounds(); ++c)
    //    {
    //        CompoundRep& scRep = updSubcompound(c).updImpl();

    //        // I am confused about how to pass an invalid compound id
    //        Compound::Index invalidCompoundIx;
    //        scRep.setCompoundSystem(system, invalidCompoundIx);
    //    }
    //}

    bool isOwnedBySystem() const {return ownerSystem != 0;}
    // const CompoundSystem& getOwnerCompoundSystem() const {assert(ownerSystem); return *ownerSystem;}
    const MultibodySystem& getOwnerMultibodySystem() const {assert(ownerSystem); return *ownerSystem;}
    // Compound::Index getIdWithinOwnerCompoundSystem() const {assert(ownerSystem); return ixWithinOwnerSystem;}

    // Add one simple atom unconnected to anything else
    CompoundRep& setBaseAtom(
        const Compound::AtomName& name, 
        const Element& element,
        const Transform& location);

    CompoundRep& setBaseAtom(
        const Compound::AtomName& name, 
        const Biotype& biotype,
        const Transform& location);

    // Add a subcompound containing exactly one atom, so the Compound::AtomName can be reused for the Compound::Name
    // This atom is not connected to anything else
    CompoundRep& setBaseAtom(
        const Compound::SingleAtom& compound,
        const Transform& location);

    // Add a subcompound without attaching it to anything
    Compound::BondCenterIndex setBaseCompound(
        const Compound::Name& n, 
        const Compound& c,
        const Transform& location);

    // Add a subcompound containing exactly one atom, so the Compound::AtomName can be reused for the Compound::Name
    // This atom is connected to existing material
    CompoundRep& bondAtom(
        const Compound::SingleAtom&   compound, 
        const Compound::BondCenterPathName& parentBondName, 
        mdunits::Length                      distance,
        Angle                         dihedral = 0,
        BondMobility::Mobility        mobility = BondMobility::Default
        );

    // Bond atom using default bond length and dihedral angle
    CompoundRep& bondAtom(
        const Compound::SingleAtom& compound, 
        const Compound::BondCenterPathName& parentBondName) 
    {
        // assert(! hasParentCompound());
        // assert(! compound.getImpl().hasParentCompound());

        // There are two choice for how to delegate this method.
        // 1) deduce the bond length and dihedral and call the other 
        //    bondAtom() that takes those parameters
        // 2) deduce the compound name and call the bondCompound()
        //    that does not take geometry parameters.
        // 
        // It is better to choose course 2), because it is less complex
        // (at the time of this writing)

        const CompoundRep& scRep = compound.getImpl();
        Compound::AtomName atomName = scRep.atomIdsByName.begin()->first;

        bondCompound(atomName, compound, parentBondName);
        inheritAtomNames(atomName);

        // assert(getSubcompound(atomName).getImpl().hasParentCompound());

        return *this;
    }

    // Add a subcompound attached by a bond to an existing atom
    // bondCompound("H1", MonovalentAtom(Element::Hydrogen()), "bond", "C/bond2", C_Hdistance );
    CompoundRep& bondCompound(
        const Compound::Name& name, 
        const Compound& subcompound, 
        const Compound::BondCenterPathName& parentBondName, 
        mdunits::Length distance,
        Angle dihedral = 0,
        BondMobility::Mobility mobility = BondMobility::Default
        );

    // Shorter version uses default bond length and dihedral angle
    CompoundRep& bondCompound(
        const Compound::Name& n, 
        const Compound& c, 
        const Compound::BondCenterPathName& parentBondName);
    // sam added polymorphism
  
    CompoundRep& bondCompound(
        const Compound::Name& n,
        const Compound& c,
        const Compound::BondCenterPathName& parentBondName,
        BondMobility::Mobility mobility
        );
    // deprecate removeSubcompound for now -- I'm not using it
    // CompoundRep& removeSubcompound(const Compound::Name& name);

    CompoundRep& setInboardBondCenter(
        const Compound::BondCenterName& centerName, 
        const Compound::AtomName& atomName, 
        Angle zRotation,
        Angle oldXRotation);

    CompoundRep& setDefaultInboardBondLength(mdunits::Length d) {
        updInboardBondCenter().setDefaultBondLength(d);
        return *this;
    }

    CompoundRep& setDefaultInboardDihedralAngle(Angle a) {
        updInboardBondCenter().setDefaultDihedralAngle(a);
        return *this;
    }

    //CompoundRep& addBondCenter(
    //    const Compound::BondCenterName& centerName, 
    //    const Compound::AtomName& atomName, 
    //    Angle zRotation,
    //    Angle oldXrotation);

    CompoundRep& addFirstBondCenter(
        const Compound::BondCenterName& centerName, 
        const Compound::AtomName& atomName);

    CompoundRep& addSecondBondCenter(
        const Compound::BondCenterName& centerName, 
        const Compound::AtomName& atomName,
        Angle bondAngle1
        );

    CompoundRep& addPlanarBondCenter(
        const Compound::BondCenterName& centerName, 
        const Compound::AtomName& atomName,
        Angle bondAngle1,
        Angle bondAngle2);

    CompoundRep& addRightHandedBondCenter(
        const Compound::BondCenterName& centerName, 
        const Compound::AtomName& atomName,
        Angle bondAngle1,
        Angle bondAngle2
        );

    CompoundRep& addLeftHandedBondCenter(
        const Compound::BondCenterName& centerName, 
        const Compound::AtomName& atomName,
        Angle bondAngle1,
        Angle bondAngle2
        );


    CompoundRep& addBondCenterInfo(
        const Compound::AtomIndex   atomId,
        const CompoundAtom::BondCenterIndex atomCenterIndex);

    CompoundRep& addRingClosingBond(
        const Compound::BondCenterName& centerName1, 
        const Compound::BondCenterName& centerName2 
        );
    CompoundRep& addRingClosingBond(
        const Compound::BondCenterName& centerName1, 
        const Compound::BondCenterName& centerName2,
        mdunits::Length bondLength,
        Angle dihedral,
        BondMobility::Mobility mobility
        );

    int getNumAtoms() const;

    const Compound::AtomName getAtomName(Compound::AtomIndex) const;

    const Element& getAtomElement(Compound::AtomIndex atomIndex) const {
        return getAtom(atomIndex).getElement();
    }

    const Element& getAtomElement(const Compound::AtomName& atomName) const {
        return getAtom(atomName).getElement();
    }


    BiotypeIndex getAtomBiotypeIndex(Compound::AtomIndex) const;
    void setDuMMAtomIndex(Compound::AtomIndex aid, DuMM::AtomIndex dummId) {
        updAtom(aid).setDuMMAtomIndex(dummId);
    }
    DuMM::AtomIndex getDuMMAtomIndex(Compound::AtomIndex aid) const {
        return getAtom(aid).getDuMMAtomIndex();
    }

    size_t getNumBondCenters() const;
	size_t getNumBondCenters(Compound::AtomIndex atomIndex) const;

    // const Compound::BondCenterName& getBondCenterName(Compound::BondCenterIndex bondCenterIndex) const;

    CompoundRep& nameAtom(const Compound::AtomName& newName, Compound::AtomIndex atomId);
    CompoundRep& nameAtom(const Compound::AtomName& newName, const Compound::AtomPathName& oldName);

    CompoundRep& nameAtom(
        const Compound::AtomName& newName, 
        const Compound::AtomPathName& oldName, 
        BiotypeIndex biotype);

    // setBiotype("C", Biotype::MethaneC);
    CompoundRep& setBiotypeIndex(const Compound::AtomName& atomName, BiotypeIndex biotype);
        
    CompoundRep& nameBondCenter(Compound::BondCenterName newName, Compound::BondCenterPathName oldName);

    // Use atoms names as found in subcompound
    CompoundRep& inheritAtomNames(const Compound::Name& scName);
    CompoundRep& inheritBondCenterNames(const Compound::Name& scName);

    bool hasDihedral(const Compound::DihedralName& angleName) const {
        return ( dihedralAnglesByName.find(angleName) != dihedralAnglesByName.end() );
    }

    bool atomsAreBonded(const AtomInfo& atom1, const AtomInfo& atom2) const 
    {
        std::pair<Compound::AtomIndex, Compound::AtomIndex> key(atom1.getIndex(), atom2.getIndex());
        return (bondIndicesByAtomIndexPair.find(key) != bondIndicesByAtomIndexPair.end());
    }


    CompoundRep& defineDihedralAngle(
        const Compound::DihedralName& angleName,
        const Compound::AtomName& atom1,
        const Compound::AtomName& atom2,
        const Compound::AtomName& atom3,
        const Compound::AtomName& atom4,
        Angle nomenclatureOffset
        ) 
    {
        assert( ! hasDihedral(angleName) );
        assert( atomsAreBonded(getAtomInfo(atom1), getAtomInfo(atom2)) );
        assert( atomsAreBonded(getAtomInfo(atom2), getAtomInfo(atom3)) );
        assert( atomsAreBonded(getAtomInfo(atom3), getAtomInfo(atom4)) );

        const BondCenterInfo& bond1 = getBondCenterInfo(atom2, atom1);
        const BondCenterInfo& bond2 = getBondCenterInfo(atom3, atom4);

        defineDihedralAngle( angleName, bond1, bond2, nomenclatureOffset );

        assert( hasDihedral(angleName) );

        return *this;
    }

    CompoundRep& defineDihedralAngle(
        const Compound::DihedralName& angleName,
        const Compound::BondCenterName& bondName1,
        const Compound::BondCenterName& bondName2,
        Angle nomenclatureOffset
        ) 
    {
        // assert( ! hasDihedral(angleName) );

        const BondCenterInfo& bond1 = getBondCenterInfo(bondName1);
        const BondCenterInfo& bond2 = getBondCenterInfo(bondName2);

        defineDihedralAngle(angleName, bond1, bond2, nomenclatureOffset);

        assert( hasDihedral(angleName) ); 

        return *this;
    }

    CompoundRep& defineDihedralAngle(
        const Compound::DihedralName& angleName,
        const BondCenterInfo& bond1,
        const BondCenterInfo& bond2,
        Angle nomenclatureOffset
        )
    {
        assert(dihedralAnglesByName.find(angleName) == dihedralAnglesByName.end());
    
        dihedralAnglesByName[angleName] = DihedralAngle(bond1.getIndex(), bond2.getIndex(), nomenclatureOffset);

        assert(dihedralAnglesByName.find(angleName) != dihedralAnglesByName.end());

        //// Define internal offset
        //DihedralAngle& dihedral = dihedralAnglesByName.find(angleName)->second;
        //const BondCenterInfo& bc21 = getBondCenterInfo(dihedral.getBondCenter1Id());
        //const BondCenterInfo& bc34 = getBondCenterInfo(dihedral.getBondCenter2Id());

        //// Find bond axis to project onto
        //const AtomInfo& atom2 = getAtomInfo(bc21.getAtomIndex());
        //const AtomInfo& atom3 = getAtomInfo(bc34.getAtomIndex());
        //const BondCenterInfo& bondBondCenter = getBondCenterInfo(atom2, atom3);

        //assert(bondBondCenter.isBonded());
        //assert(bondBondCenter.getIndex() != bc21.getIndex());
        //assert(bondBondCenter.getIndex() != bc34.getIndex());
        //assert(bc21.getIndex() != bc34.getIndex());

        //UnitVec3 xAxis(1,0,0);

        //// vector v1: from atom 1 to atom 2
        //Transform C_X_A2 = calcDefaultAtomFrameInCompoundFrame(atom2.getIndex());
        //Transform A2_X_BC21 = calcDefaultBondCenterFrameInAtomFrame(bc21);
        //Transform C_X_BC21 = C_X_A2 * A2_X_BC21;
        //UnitVec3 v1(C_X_BC21 * -xAxis); // negative x-axis because want 1->2, not 2->1 vector

        //// vector v2: from atom 2 to atom 3
        //Transform A2_X_BCB = calcDefaultBondCenterFrameInAtomFrame(bondBondCenter);
        //Transform C_X_BCB = C_X_A2 * A2_X_BCB;
        //UnitVec3 v2(C_X_BCB * xAxis);

        //// vector v3: from atom 3 to atom 4
        //Transform C_X_A3 = calcDefaultAtomFrameInCompoundFrame(atom3.getIndex());
        //Transform A3_X_BC34 = calcDefaultBondCenterFrameInAtomFrame(bc34);
        //Transform C_X_BC34 = C_X_A3 * A3_X_BC34;
        //UnitVec3 v3(C_X_BC34 * xAxis);

        //Angle nominalDihedralAngle = calcDihedralAngle(v1, v2, v3);

        //const Bond& bond = getBond(getBondInfo(bondBondCenter.getBondIndex()));
        //Angle internalDihedralAngle = bond.getDefaultDihedralAngle();

        //// internal + offset = nominal
        //Angle offset = nominalDihedralAngle - internalDihedralAngle;
        //dihedral.setInternalOffset(offset);

        return *this;
    }

    Bond& updBondByDihedral(DihedralAngle& dihedral) 
    {
        const BondCenterInfo& bc1 = getBondCenterInfo(dihedral.getBondCenter1Id());
        const BondCenterInfo& bc2 = getBondCenterInfo(dihedral.getBondCenter2Id());

        const AtomInfo& atom1 = getAtomInfo(bc1.getAtomIndex());
        const AtomInfo& atom2 = getAtomInfo(bc2.getAtomIndex());
        assert( atomsAreBonded(atom1, atom2) );

        BondInfo& bondInfo = updBondInfo(atom1, atom2);
        Bond& bond = updBond(bondInfo);

        return bond;
    }

    Bond& updBondByDihedralName(const String& bondName) 
    {
        assert( dihedralAnglesByName.find(bondName) != dihedralAnglesByName.end() );
        DihedralAngle& dihedral = dihedralAnglesByName.find(bondName)->second;
        return  updBondByDihedral(dihedral);
    }

    const Bond& getBondByDihedral(const DihedralAngle& dihedral) const 
    {
        const BondCenterInfo& bc1 = getBondCenterInfo(dihedral.getBondCenter1Id());
        const BondCenterInfo& bc2 = getBondCenterInfo(dihedral.getBondCenter2Id());

        const AtomInfo& atom1 = getAtomInfo(bc1.getAtomIndex());
        const AtomInfo& atom2 = getAtomInfo(bc2.getAtomIndex());
        assert( atomsAreBonded(atom1, atom2) );

        const BondInfo& bondInfo = getBondInfo(atom1, atom2);
        const Bond& bond = getBond(bondInfo);

        return bond;
    }

    const Bond& getBondByDihedralName(const String& dihedralName) const {
        assert( dihedralAnglesByName.find(dihedralName) != dihedralAnglesByName.end() );

        const DihedralAngle& dihedral = dihedralAnglesByName.find(dihedralName)->second;

        return getBondByDihedral(dihedral);
    }


    /**
     * \brief Sets dihedral angle without modifying bond-length or bond-angles.
     *
     * Modifying bond-angles, on the other hand, can modify those dihedral angles that involve
     * BondCenters other than the first two BondCenters on each atom.
     */
    CompoundRep& setDefaultDihedralAngle( 
            Angle angle, 
            Compound::AtomIndex atomIndex1, 
            Compound::AtomIndex atomIndex2, 
            Compound::AtomIndex atomIndex3, 
            Compound::AtomIndex atomIndex4)
    {
        const BondCenterInfo& bondCenterInfo21 = getBondCenterInfo( getAtomInfo(atomIndex2), getAtomInfo(atomIndex1) );
        const BondCenterInfo& bondCenterInfo34 = getBondCenterInfo( getAtomInfo(atomIndex3), getAtomInfo(atomIndex4) );

        // for debugging
        //String atom1Name = getAtomName(atomIndex1);
        //String atom2Name = getAtomName(atomIndex2);
        //String atom3Name = getAtomName(atomIndex3);
        //String atom4Name = getAtomName(atomIndex4);
        //std::cout << atom1Name << "->" << atom2Name << "->" << atom3Name << "->" << atom4Name << std::endl;

        return setDefaultDihedralAngle( angle, bondCenterInfo21.getIndex(), bondCenterInfo34.getIndex() );
    }


    CompoundRep& setDefaultDihedralAngle( 
            Angle angle, 
            Compound::AtomName atom1, 
            Compound::AtomName atom2, 
            Compound::AtomName atom3, 
            Compound::AtomName atom4)
    {
    	return setDefaultDihedralAngle(angle, 
    			getAtomInfo(atom1).getIndex(),
    			getAtomInfo(atom2).getIndex(),
    			getAtomInfo(atom3).getIndex(),
    			getAtomInfo(atom4).getIndex() );
    }
    
    // determine difference, in radians, between dihedral defined by these bond centers (nominal),
    // and dihedral defined by "canonical" bond centers (internal).
    // nominal = internal + offset => offset = nominal - internal
    Angle calcDefaultInternalDihedralOffsetAngle(
            Compound::BondCenterIndex bondCenterIndex21, 
            Compound::BondCenterIndex bondCenterIndex34) const
    {
        Compound::AtomIndex atomIndex2 = getBondCenterInfo(bondCenterIndex21).getAtomIndex();
        Compound::AtomIndex atomIndex3 = getBondCenterInfo(bondCenterIndex34).getAtomIndex();

        const AtomInfo& atomInfo2 = getAtomInfo(atomIndex2);
        const AtomInfo& atomInfo3 = getAtomInfo(atomIndex3);

        // Sanity check topology
        assert( atomsAreBonded(atomInfo2, atomInfo3) ); // absolutely required

        // Find central bond
        const BondInfo& bondInfo23 = getBondInfo(atomInfo2, atomInfo3);
        const Bond& bond23 = getBond(bondInfo23);

        // Identify the bond centers associated with the atom2-atom3 bond
        const BondCenterInfo& bondCenterInfo23 = getBondCenterInfo(atomInfo2, atomInfo3);
        const BondCenterInfo& bondCenterInfo32 = getBondCenterInfo(atomInfo3, atomInfo2);
        // sanity check those central bond centers
        assert(bondCenterInfo23.getAtomIndex() == atomIndex2);
        assert(bondCenterInfo32.getAtomIndex() == atomIndex3);

        // 1) Identify canonical bond centers for internal dihedral angle
        // Usually bond-center number zero(0), unless zero participates in the atom2-atom3 bond
        CompoundAtom::BondCenterIndex canonicalCenterIndex2(0); // default to zero
        if (bondCenterInfo23.getAtomBondCenterIndex() == 0) // unless zero is used for 2->3 bond
            canonicalCenterIndex2 = CompoundAtom::BondCenterIndex(1);

        CompoundAtom::BondCenterIndex canonicalCenterIndex3(0); // default to zero
        if (bondCenterInfo32.getAtomBondCenterIndex() == 0) // unless zero is used for 2->3 bond
            canonicalCenterIndex3 = CompoundAtom::BondCenterIndex(1);

        // debug
        // Compound::AtomName n2 = getAtomName(atomIndex2);
        // Compound::AtomName n3 = getAtomName(atomIndex3);

        // 2) Compute offsets for actual bond centers
        // * offsetAngle1 is counter-clockwise angle from canonical bond center on atom2 to atom1, viewed
        // down the atom3-atom2 axis.
        const BondCenterInfo& bondCenterInfo21 = getBondCenterInfo(bondCenterIndex21);
        Angle offsetAngle1;
        if (canonicalCenterIndex2 == bondCenterInfo21.getAtomBondCenterIndex())
            offsetAngle1 = 0.0;
        else
        {
            // trick the bond-vector version of calcDihedralAngle into giving the offset angle at the atom
            const CompoundAtom& atom2 = getAtom(atomIndex2);
            UnitVec3 dirAtom1    = -atom2.getBondCenterDirectionInAtomFrame(bondCenterInfo21.getAtomBondCenterIndex());
            UnitVec3 dirBond     = atom2.getBondCenterDirectionInAtomFrame(bondCenterInfo23.getAtomBondCenterIndex());
            UnitVec3 dirRefAtom1 = atom2.getBondCenterDirectionInAtomFrame(canonicalCenterIndex2);

            // Sometimes bond direction is colinear with atom direction, if chirality is hosed
            double problemCheck = std::abs(dot(dirBond, dirAtom1));
            if (problemCheck > 0.999)
                offsetAngle1 = 0.0;
            else
                offsetAngle1 = SimTK::calcDihedralAngle(dirRefAtom1, dirBond, dirAtom1);

            // assert(offsetAngle1 != 0);
        }

        // * offsetAngle4 is counter-clockwise angle from canonical bond center on atom3 to atom4, viewed
        // down the atom3-atom2 axis.
        const BondCenterInfo& bondCenterInfo34 = getBondCenterInfo(bondCenterIndex34);
        Angle offsetAngle4;
        if (canonicalCenterIndex3 == bondCenterInfo34.getAtomBondCenterIndex())
            offsetAngle4 = 0.0;
        else
        {
            // trick the bond-vector version of calcDihedralAngle into giving the offset angle at the atom
            const CompoundAtom& atom3 = getAtom(atomIndex3);
            UnitVec3 dirAtom4    = -atom3.getBondCenterDirectionInAtomFrame(bondCenterInfo34.getAtomBondCenterIndex());
            UnitVec3 dirBond     = -atom3.getBondCenterDirectionInAtomFrame(bondCenterInfo32.getAtomBondCenterIndex());
            UnitVec3 dirRefAtom4 = atom3.getBondCenterDirectionInAtomFrame(canonicalCenterIndex3);

            // Sometimes bond direction is colinear with atom direction, if chirality is hosed
            double problemCheck = std::abs(dot(dirBond, dirAtom4));
            if (problemCheck > 0.999)
                offsetAngle1 = 0.0;
            else
                offsetAngle4 = SimTK::calcDihedralAngle(dirRefAtom4, dirBond, dirAtom4);

            // assert(offsetAngle4 != 0);
        }

        // nominal = internal + offset
        // offset = nominal - internal
        // internal = nominal - offset
        // Angle internalDihedralAngle = angle + offsetAngle1 - offsetAngle4;

        Angle offsetAngle = offsetAngle4 - offsetAngle1;

        while ( -SimTK::Pi >= offsetAngle ) offsetAngle += 2 * SimTK::Pi;
        while ( SimTK::Pi < offsetAngle ) offsetAngle -= 2 * SimTK::Pi;

        // debugging
        //std::cout << "  total offset = " << offsetAngle * DuMM::Rad2Deg;
        //std::cout << "; offset1 = " << offsetAngle1 * DuMM::Rad2Deg;
        //std::cout << "; offset4 = " << offsetAngle4 * DuMM::Rad2Deg << std::endl;

        return offsetAngle;
    }


    /**
     * \brief Sets dihedral angle without modifying bond-length or bond-angles.
     *
     * Modifying bond-angles, on the other hand, can modify those dihedral angles that involve
     * BondCenters other than the first two BondCenters on each atom.
     */
    CompoundRep& setDefaultDihedralAngle( 
            Angle angle, 
            Compound::BondCenterIndex bondCenterIndex21, 
            Compound::BondCenterIndex bondCenterIndex34)
    {
        Compound::AtomIndex atomIndex2 = getBondCenterInfo(bondCenterIndex21).getAtomIndex();
        Compound::AtomIndex atomIndex3 = getBondCenterInfo(bondCenterIndex34).getAtomIndex();

        const AtomInfo& atomInfo2 = getAtomInfo(atomIndex2);
        const AtomInfo& atomInfo3 = getAtomInfo(atomIndex3);

        // Sanity check topology
        assert( atomsAreBonded(atomInfo2, atomInfo3) ); // absolutely required

        // Find central bond
        BondInfo& bondInfo23 = updBondInfo(atomInfo2, atomInfo3);
        Bond& bond23 = updBond(bondInfo23);

        Angle internalDihedralAngle = angle - calcDefaultInternalDihedralOffsetAngle(bondCenterIndex21, bondCenterIndex34);

        while ( -SimTK::Pi >= internalDihedralAngle ) internalDihedralAngle += 2 * SimTK::Pi;
        while ( SimTK::Pi < internalDihedralAngle ) internalDihedralAngle -= 2 * SimTK::Pi;

        //std::cout << "old internal angle = " << bond23.getDefaultDihedralAngle() * DuMM::Rad2Deg << std::endl;
        //std::cout << "new internal angle = " << internalDihedralAngle * DuMM::Rad2Deg << std::endl;

		// debug - notice when angle changes
		//Real diff = internalDihedralAngle - bond23.getDefaultDihedral();
  //      while ( -SimTK::Pi >= diff ) diff += 2 * SimTK::Pi;
  //      while ( SimTK::Pi < diff ) diff -= 2 * SimTK::Pi;
		//diff = diff < 0 ? -diff : diff;
		//if (diff > 0.005) 
		//{
		//	int x = 1;
		//}

        bond23.setDefaultDihedralAngle(internalDihedralAngle);

        return *this;
    }

    // setDefaultDihedral changes no bond lengths or bond angles
    CompoundRep& setDefaultDihedralAngle(const String& dihedralName, Angle finalNominalAngle) 
    {
        Bond& bond = updBondByDihedralName(dihedralName);
        DihedralAngle& dihedral = dihedralAnglesByName.find(dihedralName)->second;

        // internal = nominal - offset
        Angle angle = finalNominalAngle - dihedral.getNomenclatureOffset();

        setDefaultDihedralAngle(angle, dihedral.getBondCenter1Id(), dihedral.getBondCenter2Id());
        // Angle internalAngle = finalNominalAngle - dihedral.getInternalOffset() - dihedral.getNomenclatureOffset();

        // bond.setDefaultDihedralAngle(internalAngle);

        return *this;
    }


    Angle calcDefaultDihedralAngle(const String& dihedralName) const 
    {
        assert( dihedralAnglesByName.find(dihedralName) != dihedralAnglesByName.end() );

        const DihedralAngle& dihedral = dihedralAnglesByName.find(dihedralName)->second;

        return calcDefaultDihedralAngle(dihedral);
    }

    Angle calcDefaultDihedralAngle(            
            Compound::BondCenterIndex bondCenterIndex21, 
            Compound::BondCenterIndex bondCenterIndex34)
    {
        Compound::AtomIndex atomIndex2 = getBondCenterInfo(bondCenterIndex21).getAtomIndex();
        Compound::AtomIndex atomIndex3 = getBondCenterInfo(bondCenterIndex34).getAtomIndex();

        const AtomInfo& atomInfo2 = getAtomInfo(atomIndex2);
        const AtomInfo& atomInfo3 = getAtomInfo(atomIndex3);

        // Sanity check topology
        assert( atomsAreBonded(atomInfo2, atomInfo3) ); // absolutely required

        // Find central bond
        const BondInfo& bondInfo23 = getBondInfo(atomInfo2, atomInfo3);
        const Bond& bond23 = getBond(bondInfo23);

        Angle internalDihedralAngle = bond23.getDefaultDihedralAngle();

        Angle nominalDihedralAngle = internalDihedralAngle + 
            calcDefaultInternalDihedralOffsetAngle(bondCenterIndex21, bondCenterIndex34);

        return nominalDihedralAngle;
    }

    Angle calcDefaultDihedralAngle( 
            Compound::AtomIndex atomIndex1, 
            Compound::AtomIndex atomIndex2, 
            Compound::AtomIndex atomIndex3, 
            Compound::AtomIndex atomIndex4)
    {
        const BondCenterInfo& bondCenterInfo21 = getBondCenterInfo( getAtomInfo(atomIndex2), getAtomInfo(atomIndex1) );
        const BondCenterInfo& bondCenterInfo34 = getBondCenterInfo( getAtomInfo(atomIndex3), getAtomInfo(atomIndex4) );

        return calcDefaultDihedralAngle( bondCenterInfo21.getIndex(), bondCenterInfo34.getIndex() );
    }

    CompoundRep& setDihedralAngle(State& state, const String& dihedralName, Angle angleInRadians) 
    {
        assert(ownerSystem != NULL);
        Bond& bond = updBondByDihedralName(dihedralName);
        DihedralAngle& dihedral = dihedralAnglesByName.find(dihedralName)->second;

        // case1 : Pin dihedral
        if (bond.getPinJointId().isValid()) {
            SimbodyMatterSubsystem& matter = ownerSystem->updMatterSubsystem();
            MobilizedBody::Pin& body = (MobilizedBody::Pin&) matter.updMobilizedBody(bond.getPinJointId());

            // TODO - create calcDihedralOffset(State&...) method and use it here, instead of default
            Angle internalOffset = calcDefaultInternalDihedralOffsetAngle(dihedral.getBondCenter1Id(), dihedral.getBondCenter2Id());

            // nominal = internal + offset
            Angle internalAngle = angleInRadians - internalOffset - dihedral.getNomenclatureOffset();
            body.setAngle(state, internalAngle);
        }

        else  // TODO
        {
            assert(false);

            // Dihedral may be offset from "standard" dihedral for bond
            assert(ownerSystem);
            SimbodyMatterSubsystem& matter = ownerSystem->updMatterSubsystem();

            Angle previousInternalDihedral = bond.getDihedralAngle(state, matter);

            Angle previousNominalDihedral = calcDihedralAngle(state, dihedralName); // requires realizePosition
            // Nominal = internal + offset
            Angle offsetAngle = previousNominalDihedral - previousInternalDihedral;
            Angle newInternalDihedral = angleInRadians - offsetAngle;
            // Restrict to range +-Pi
            while (newInternalDihedral <= -Pi) newInternalDihedral += 2*Pi;
            while (newInternalDihedral > Pi) newInternalDihedral -= 2*Pi;
            bond.setDihedralAngle(state, matter, newInternalDihedral); // clears realizePosition

            Angle testAngle = calcDihedralAngle(state, dihedralName); // requires realizePosition
            Angle error = calcDihedralAngle(state, dihedralName) - testAngle;
            while (error <= -Pi) error += 2*Pi;
            while (error > Pi) error -= 2*Pi;
            error *= error;
            assert(error < 1e-6);
        }

        return *this;
    }

    CompoundRep& setDefaultBond1Angle(const String& bondName, Angle angle) {
        updBondCenter(bondName).setDefaultBond1Angle(angle);
        return *this;
    }
    CompoundRep& setDefaultBond2Angle(const String& bondName, Angle angle) {
        updBondCenter(bondName).setDefaultBond2Angle(angle);
        return *this;
    }

    CompoundRep& setDefaultBondAngle(
        Angle angle, 
        const Compound::AtomName& atom1Name, 
        const Compound::AtomName& atom2Name, 
        const Compound::AtomName& atom3Name) 
    {
        const Compound::AtomIndex atom1Id   = getAtomInfo(atom1Name).getIndex();
        const Compound::AtomIndex atom2Id   = getAtomInfo(atom2Name).getIndex();
        const Compound::AtomIndex atom3Id   = getAtomInfo(atom3Name).getIndex();
        
        return setDefaultBondAngle(angle, atom1Id, atom2Id, atom3Id);
    }
    
    // Version that takes IDs, to reduce string lookups
    CompoundRep& setDefaultBondAngle(
        Angle angle, 
        const Compound::AtomIndex atom1Id, 
        const Compound::AtomIndex atom2Id, 
        const Compound::AtomIndex atom3Id) 
    {
        CompoundAtom& atom2 = updAtom(atom2Id);
        const AtomInfo&        atom2Info = getAtomInfo(atom2Id);

        // go through bond centers on atom2
        CompoundAtom::BondCenterIndex center1Id;
        CompoundAtom::BondCenterIndex center3Id;
        for (CompoundAtom::BondCenterIndex b(0); b < atom2.getNumBonds(); ++b) {
            const BondCenterInfo& bondCenterInfo = getBondCenterInfo(atom2Info.getIndex(), b);
            if (bondCenterInfo.isBonded()) {
                const BondCenterInfo& partnerBondCenterInfo = 
                    getBondCenterInfo(bondCenterInfo.getBondPartnerBondCenterIndex());
                if (partnerBondCenterInfo.getAtomIndex() == atom1Id)
                    center1Id = b;
                else if (partnerBondCenterInfo.getAtomIndex() == atom3Id)
                    center3Id = b;
            }
        }

        assert(center1Id.isValid());
        assert(center3Id.isValid());
        assert(center1Id != center3Id);

        CompoundAtom::BondCenterIndex largerId, smallerId;
        if (center1Id > center3Id) {
            largerId = center1Id;
            smallerId = center3Id;
        } else {
            largerId = center3Id;
            smallerId = center1Id;
        }

        // one of the bond centers must be bond1 or bond2
        // assert(smallerId < 2);

        if (smallerId == 0) {
            atom2.updBondCenter(largerId).setDefaultBond1Angle(angle);
        }
        else if (smallerId == 1) {
            atom2.updBondCenter(largerId).setDefaultBond2Angle(angle);
        }
        else {
            // assert(false); // TODO
        }

        return *this;
    }

    CompoundRep& setDefaultBondLength(mdunits::Length length, const Compound::AtomName& atom1Name, const Compound::AtomName& atom2Name) 
    {
        const CompoundAtom&             atom2       = getAtom(atom2Name);
        const AtomInfo&         atom2Info   = getAtomInfo(atom2Name);
        const Compound::AtomIndex  atom1Id     = getAtomInfo(atom1Name).getIndex();

        // go through bond centers on atom2
        CompoundAtom::BondCenterIndex center1Id;
        for (CompoundAtom::BondCenterIndex b(0); b < atom2.getNumBonds(); ++b) {
            const BondCenterInfo& bondCenterInfo = getBondCenterInfo(atom2Info.getIndex(), b);
            if (getBondCenter(bondCenterInfo).isBonded()) {
                const BondCenterInfo& partnerBondCenterInfo = getBondCenterInfo(bondCenterInfo.getBondPartnerBondCenterIndex());
                if (partnerBondCenterInfo.getAtomIndex() == atom1Id) 
                {
                    center1Id = b;

                    // set length from atom1 direction
                    updBondCenter(bondCenterInfo).setDefaultBondLength(length);

                    // for good measure, set length from atom2 direction
                    updBondCenter(partnerBondCenterInfo).setDefaultBondLength(length);

                    break;
                }
            }
        }

        assert(center1Id.isValid());

        return *this;
    }

    // More efficient getting of all atoms at once
    // Compound::AtomTargetLocations calcDefaultAtomLocationsInCompoundFrame1() const;
    void calcDefaultAtomFramesInCompoundFrame(std::vector<Transform>& atomFrameCache) const;
//    Compound::AtomTargetLocations calcDefaultAtomLocationsInCompoundFrame() const;
//    void calcDefaultAtomLocationsInParentFrame(
//            Compound::AtomTargetLocations& locations,
//            std::vector<Compound::AtomIndex>& indexes,  // map atom index to location index
//            const Transform& X_ParentChild
//            ) const;  
    
    
    // Compute atom location in local compound frame
    Transform calcDefaultAtomFrameInCompoundFrame(const Compound::AtomName& name) const;
    Transform calcDefaultAtomFrameInGroundFrame(const Compound::AtomName& name) const;
    Vec3 calcDefaultAtomLocationInCompoundFrame(const Compound::AtomName& name) const;
    Vec3 calcDefaultAtomLocationInGroundFrame(const Compound::AtomName& name) const;

    MobilizedBodyIndex getAtomMobilizedBodyIndex(Compound::AtomIndex atomId) const 
    {
        const CompoundAtom& atom = getAtom(atomId);
        return atom.getMobilizedBodyIndex();
    }
    Vec3 getAtomLocationInMobilizedBodyFrame(Compound::AtomIndex atomId) const {
        const CompoundAtom& atom = getAtom(atomId);
        return atom.getLocationInMobilizedBodyFrame();
    }
    Vec3 calcAtomLocationInGroundFrame(const State& state, Compound::AtomIndex atomId) const {
        ownerSystem->realize(state, Stage::Position);
        const CompoundAtom& atom = getAtom(atomId);
        Vec3 loc = atom.getLocationInMobilizedBodyFrame();
        const SimbodyMatterSubsystem& matter = ownerSystem->getMatterSubsystem();
        const MobilizedBody& body = matter.getMobilizedBody(getAtomMobilizedBodyIndex(atomId));
        return body.getBodyTransform(state)*loc;
    }
    Vec3 calcAtomVelocityInGroundFrame(const State& state, Compound::AtomIndex atomId) const {
        ownerSystem->realize(state, Stage::Velocity);
        const CompoundAtom& atom = getAtom(atomId);
        Vec3 loc = atom.getLocationInMobilizedBodyFrame();
        const SimbodyMatterSubsystem& matter = ownerSystem->getMatterSubsystem();
        const MobilizedBody& body = matter.getMobilizedBody(getAtomMobilizedBodyIndex(atomId));
        return body.findStationVelocityInGround(state, loc);
    }
    Vec3 calcAtomAccelerationInGroundFrame(const State& state, Compound::AtomIndex atomId) const {
        ownerSystem->realize(state, Stage::Acceleration);
        const CompoundAtom& atom = getAtom(atomId);
        Vec3 loc = atom.getLocationInMobilizedBodyFrame();
        const SimbodyMatterSubsystem& matter = ownerSystem->getMatterSubsystem();
        const MobilizedBody& body = matter.getMobilizedBody(getAtomMobilizedBodyIndex(atomId));
        return body.findStationAccelerationInGround(state, loc);
    }

    
    Transform calcDefaultBondCenterFrameInCompoundFrame(const BondCenterInfo& info) const;

    Transform calcDefaultAtomFrameInCompoundFrame(Compound::AtomIndex atomId) const;
    // Version with caching for O(n) performance
    const Transform& calcDefaultAtomFrameInCompoundFrame(Compound::AtomIndex atomId, std::vector<Transform>& atomFrameCache) const;

    Transform calcDefaultAtomFrameInGroundFrame(Compound::AtomIndex atomId) const;

    typedef std::vector<Compound::AtomIndex> AtomIndexList;
    // get list of all runs of consecutive bonded atoms of run-length n from the atoms mentions in an AtomTargetLocations structure
    // for example, to get a list of all bonded pairs, set run-length to 2.
    std::vector< AtomIndexList > getBondedAtomRuns(int atomRunCount, const Compound::AtomTargetLocations& atomTargets) const 
    {
        std::vector< AtomIndexList >  answer;

        typedef std::map<Compound::AtomIndex, AtomIndexList > BondMap;
        BondMap bondMap;

        // 1) Hash bonding data
        for (Compound::BondIndex b(0); b < getNumBonds(); ++b) 
        {
            const BondInfo& bondInfo = getBondInfo(b);

            // ignore bonds without known atom positions at both ends
            // i.e. keep the previous default bond lengths for those
            Compound::AtomIndex atomIndex1 = getBondCenterInfo(bondInfo.getParentBondCenterIndex()).getAtomIndex();
            Compound::AtomIndex atomIndex2 = getBondCenterInfo(bondInfo.getChildBondCenterIndex()).getAtomIndex();
            if (atomTargets.find(atomIndex1) == atomTargets.end()) continue;
            if (atomTargets.find(atomIndex2) == atomTargets.end()) continue;

            assert(atomIndex1 != atomIndex2);

            if (bondMap.find(atomIndex1) == bondMap.end()) bondMap[atomIndex1] = AtomIndexList();
            bondMap[atomIndex1].push_back(atomIndex2);

            if (bondMap.find(atomIndex2) == bondMap.end()) bondMap[atomIndex2] = AtomIndexList();
            bondMap[atomIndex2].push_back(atomIndex1);

            // Update n==2 version of answer
            answer.push_back(AtomIndexList());
            answer.back().push_back(atomIndex1);
            answer.back().push_back(atomIndex2);

            // and in reverse order
            answer.push_back(AtomIndexList());
            answer.back().push_back(atomIndex2);
            answer.back().push_back(atomIndex1);
        }

        // 2) Create bonded atom run lists
        for (int n = 3; n <= atomRunCount; ++n)
        {
            std::vector< AtomIndexList >  newAnswer;
            std::vector< AtomIndexList >::const_iterator oldRunI;
            for (oldRunI = answer.begin(); oldRunI != answer.end(); ++oldRunI)
            {
                // remember which atoms are already in this run
                std::set<Compound::AtomIndex> oldAtoms;
                AtomIndexList::const_iterator oldAtomI;
                for (oldAtomI = oldRunI->begin(); oldAtomI != oldRunI->end(); ++oldAtomI)
                    oldAtoms.insert(*oldAtomI);

                // look for new atoms bonded to end of old run
                Compound::AtomIndex oldTailIndex = oldRunI->back(); // final atom of shorter run
                const AtomIndexList& newTailCandidates = bondMap[oldTailIndex];
                AtomIndexList::const_iterator newTailI;
                for (newTailI = newTailCandidates.begin(); newTailI != newTailCandidates.end(); ++newTailI)
                {
                    if (oldAtoms.find(*newTailI) != oldAtoms.end()) continue; // ignore all atoms already in this run

                    // first place a copy of the old run in the new set
                    newAnswer.push_back(*oldRunI);
                    // then add the new atom to the end of the run
                    newAnswer.back().push_back(*newTailI);
                }
            }

            // overwrite the old set of shorter runs with the new one containing longer runs.
            answer = newAnswer;
        }

        return answer;
    }




    CompoundRep& matchDefaultAtomChirality(const Compound::AtomTargetLocations& atomTargets, Angle breakPlanarityThreshold, bool flipAll=true)
    {
        std::vector< AtomIndexList > atomPairs = getBondedAtomRuns(2, atomTargets);

        // 1 Ignore atoms with less than three known atoms bonded
        // 1a count bonds to each atom from other target atoms
    
        typedef std::set<Compound::AtomIndex> AtomIndexSet;
        typedef std::map< Compound::AtomIndex, AtomIndexSet > AtomPartnerList;
        AtomPartnerList atomPartnerList;
        // Initialize counts to zero
        Compound::AtomTargetLocations::const_iterator atomI;
        for (atomI = atomTargets.begin(); atomI != atomTargets.end(); ++atomI)
            atomPartnerList[atomI->first] = AtomIndexSet();

        // Enumerate the bonds for each atom
        std::vector< AtomIndexList >::const_iterator bondI;
        for (bondI = atomPairs.begin(); bondI != atomPairs.end(); ++bondI)
        {
            Compound::AtomIndex atomIndex1 = (*bondI)[0];
            Compound::AtomIndex atomIndex2 = (*bondI)[1];

            atomPartnerList[atomIndex1].insert(atomIndex2);
            atomPartnerList[atomIndex2].insert(atomIndex1);
        }

        // Check the chirality of each atom
        // Loop over atoms
        for (atomI = atomTargets.begin(); atomI != atomTargets.end(); ++atomI) {
            Compound::AtomIndex atomIndex = atomI->first;
            const AtomIndexSet& neighborAtomIndexes = atomPartnerList[atomIndex];

            int numberOfBonds = neighborAtomIndexes.size();

            // std::cout << "number of bonds = " << numberOfBonds << std::endl;

            if (numberOfBonds < 3)
                continue; // cannot be a chiral disagreement

            Vec3 targetCenter = atomI->second;
            CompoundAtom& atom = updAtom(atomIndex);

            // Note direction and chirality of each bond center in target structure
            // Compute (unit) bond vectors for source and target
            std::vector< std::pair<UnitVec3, UnitVec3> > bondVectors;
            std::vector< Compound::BondCenterIndex > bondCenterIndexes;
            AtomIndexSet::const_iterator neighborI;
            for (neighborI = neighborAtomIndexes.begin(); neighborI != neighborAtomIndexes.end(); ++neighborI)
            {
                // 1) Target bond direction
                Compound::AtomIndex neighborIndex = atomTargets.find(*neighborI)->first;
                UnitVec3 targetDirection(atomTargets.find(*neighborI)->second - targetCenter);

                // 2) Compound bond direction
                // 2a identify bond center
                const BondCenterInfo& bondCenterInfo = 
                    getBondCenterInfo( getAtomInfo(atomIndex), getAtomInfo(neighborIndex) );
                CompoundAtom::BondCenterIndex atomBondCenterIndex = bondCenterInfo.getAtomBondCenterIndex();
                // 2b get bond center frame with respect to atom
                Transform bondCenterFrame = atom.calcDefaultBondCenterFrameInAtomFrame(atomBondCenterIndex);
                // 2c apply frame to x-axis (abitrary axis) to get relative direction
                UnitVec3 sourceDirection(bondCenterFrame * UnitVec3(1, 0, 0));

                // store molecule and target structure bond directions...
                bondVectors.push_back(std::pair<Vec3, Vec3>(sourceDirection, targetDirection));

                // ... and store bond center index
                bondCenterIndexes.push_back( bondCenterInfo.getIndex() );
            }
            assert(numberOfBonds == bondVectors.size());

            // Because we need to distinguish left-handed from right-handed geometry in the atom frame,
            // we should use bondcenters number 0 and 1 to define the plane, so that the target structure
            // geometry matches that of the internal atom geometry
            // So identify a mapping between internal atom BondCenter indices and the recently constructed
            // neighbor atom indices
            int zeroBondCenterIndex = 0;
            int oneBondCenterIndex = 1;
            int twoBondCenterIndex = 2;
            for (int bondIx = 0; bondIx < (int) bondVectors.size(); ++bondIx) {
                const BondCenterInfo& bondCenterInfo = getBondCenterInfo(bondCenterIndexes[bondIx]);
                if (bondCenterInfo.getAtomBondCenterIndex() == 0) {
                    // permute indices if necessary
                    if ( oneBondCenterIndex == bondIx ) oneBondCenterIndex = zeroBondCenterIndex;
                    else if ( twoBondCenterIndex == bondIx ) twoBondCenterIndex = zeroBondCenterIndex;

                    zeroBondCenterIndex = bondIx;
                }
                else if (bondCenterInfo.getAtomBondCenterIndex() == 1) {
                    // permute indices if necessary
                    if ( zeroBondCenterIndex == bondIx ) zeroBondCenterIndex = oneBondCenterIndex;
                    else if ( twoBondCenterIndex == bondIx ) twoBondCenterIndex = oneBondCenterIndex;

                    oneBondCenterIndex = bondIx;
                }
            }

            // Use the first three atoms to detect chirality
            // This should work well for 3 and four atom case
            // With more than four bonded atoms, well... that's tricky.
            // std::cout << "Chiral center found: " << atomI->first << std::endl;

            // Three ordered source vectors
            UnitVec3 s1 = bondVectors[zeroBondCenterIndex].first;
            UnitVec3 s2 = bondVectors[oneBondCenterIndex].first;
            UnitVec3 s3 = bondVectors[twoBondCenterIndex].first;

            // And three ordered target vectors
            UnitVec3 t1 = bondVectors[zeroBondCenterIndex].second;
            UnitVec3 t2 = bondVectors[oneBondCenterIndex].second;
            UnitVec3 t3 = bondVectors[twoBondCenterIndex].second;

            // Break planar groups if target structure is farther than <tolerance> from planar
            // Determine if any supposedly planar bonds from this atom are 
            // significantly out of plane in the target structure
            UnitVec3 targetPlaneNormal(cross(t1, t2));
            bool doBreakPlane = false;
            for (int bondIx = 0; bondIx < (int) bondVectors.size(); ++bondIx) {

                const BondCenter& bondCenter = getBondCenter(bondCenterIndexes[bondIx]);

                // Can't break planarity if it's not planar to begin with
                if (bondCenter.getChirality() != BondCenter::Planar) 
                    continue;

                // Don't break planarity if atom is close enough to planar in target structure
                Angle sinePlaneDeviation = dot(bondVectors[bondIx].second, targetPlaneNormal);
                if ( std::abs(sinePlaneDeviation) < std::sin(breakPlanarityThreshold) )
                    continue;

                doBreakPlane = true;

                break;
            }
            // If one bond is out of plane, put them all out of plane
            if (doBreakPlane) {
                for (int bondIx = 0; bondIx < (int) bondVectors.size(); ++bondIx) {

                    BondCenter& bondCenter = updBondCenter(bondCenterIndexes[bondIx]);

                    // Can't break planarity if it's not planar to begin with
                    if (bondCenter.getChirality() != BondCenter::Planar) 
                        continue;

                    // Don't try to break planarity of first and second bond centers
                    // because they define the plane
                    const BondCenterInfo& bondCenterInfo = getBondCenterInfo(bondCenterIndexes[bondIx]);
                    if (bondCenterInfo.getAtomBondCenterIndex() == 0) continue;
                    if (bondCenterInfo.getAtomBondCenterIndex() == 1) continue;

                    // OK, if we got this far, we must break planarity
                    std::cerr << "WARNING: matching out-of-plane atoms about atom ";
                    std::cerr << getAtomName(atomIndex);
                    std::cerr << std::endl;

                    Angle sinePlaneDeviation = dot(bondVectors[bondIx].second, targetPlaneNormal);

                    if (sinePlaneDeviation < 0) 
                        bondCenter.setChirality(BondCenter::LeftHanded);
                    else
                        bondCenter.setChirality(BondCenter::RightHanded);
                }
            }

            if (flipAll) { // flip entire atom or nothing
                Real sourceChirality = dot(cross(s1, s2),s3);
                Real targetChirality = dot(cross(t1, t2),t3);

                // std::cout << "Source chirality = " << sourceChirality << std::endl;
                // std::cout << "Target chirality = " << targetChirality << std::endl;


                // Reverse chirality of bond centers that differ from those in target structure
                // same sign means same chirality
                if (sourceChirality * targetChirality < 0)
                { // mismatch
                    std::cerr << "WARNING: Using unexpected chirality about atom ";
                    std::cerr << getAtomName(atomIndex);
                    std::cerr << std::endl;

                    // flip the chirality of every handed bond center in the target atom
                    for (CompoundAtom::BondCenterIndex bcIx(0); bcIx < atom.getNumBonds(); ++bcIx)
                    {
                        BondCenter& bondCenter = updBondCenter(getBondCenterInfo(atomIndex, bcIx));
                        BondCenter::Chirality newChirality = bondCenter.getChirality();
                        switch(bondCenter.getChirality()) {
                            case BondCenter::RightHanded:
                                newChirality = BondCenter::LeftHanded;
                                break;
                            case BondCenter::LeftHanded:
                                newChirality = BondCenter::RightHanded;
                                break;
                            default:
                                break;
                        }
                        bondCenter.setChirality(newChirality);
                    }
                } // end if chirality differs
            }
            else { // flip on a bondcenter by bondcenter basis
                // Flip chirality on a BondCenter by BondCenter basis
                for (int bondIx = 0; bondIx < (int) bondVectors.size(); ++bondIx) 
                {
                    // First two bond centers cannot be chiral
                    if (bondIx == zeroBondCenterIndex) continue;
                    if (bondIx == oneBondCenterIndex) continue;

                    // Measure source and target chiralities
                    UnitVec3 sourceBondVec = bondVectors[bondIx].first;
                    UnitVec3 targetBondVec = bondVectors[bondIx].second;
                    Real sourceChirality = dot(cross(s1, s2),sourceBondVec);
                    Real targetChirality = dot(cross(t1, t2),targetBondVec);

                    if (sourceChirality * targetChirality < 0) // mismatched chirality
                    {
                        const BondCenterInfo& bondCenterInfo = getBondCenterInfo(bondCenterIndexes[bondIx]);
                        Compound::AtomIndex partnerAtomIndex = 
                            getBondCenterInfo(bondCenterInfo.getBondPartnerBondCenterIndex())
                            .getAtomIndex();
                        std::cerr << "WARNING: Flipping chirality of bond from atom ";
                        std::cerr << getAtomName(atomIndex);
                        std::cerr << " to atom ";
                        std::cerr << getAtomName(partnerAtomIndex);
                        std::cerr << std::endl;
                        BondCenter& bondCenter = updBondCenter(bondCenterInfo);
                        switch(bondCenter.getChirality()) {
                            case BondCenter::RightHanded:
                                bondCenter.setChirality(BondCenter::LeftHanded);
                                break;
                            case BondCenter::LeftHanded:
                                bondCenter.setChirality(BondCenter::RightHanded);
                                break;
                        }
                    }
                }
            }
        } // end for atoms

        return *this;
    }

    CompoundRep& matchDefaultBondLengths(const Compound::AtomTargetLocations& atomTargets) 
    {
        std::vector< AtomIndexList > atomPairs = getBondedAtomRuns(2, atomTargets);

        // Loop over those pairs of atoms and set the bond length default to the target distances
        // This method is broken into two parts like this to serve as an example for the more
        // complex methods to follow.
        std::vector< AtomIndexList >::const_iterator bonds12Ix;
        for (bonds12Ix = atomPairs.begin(); bonds12Ix != atomPairs.end(); ++bonds12Ix) 
        {
            Compound::AtomIndex atomIndex1 = (*bonds12Ix)[0];
            Compound::AtomIndex atomIndex2 = (*bonds12Ix)[1];

            // For efficiency, only set bonds lengths once per bond
            if (atomIndex2 > atomIndex1) continue;

            // compute distance
            Vec3 d3 = atomTargets.find(atomIndex1)->second - atomTargets.find(atomIndex2)->second;
            Real distance = std::sqrt(dot(d3, d3));

            // set bond length
            updBond(updBondInfo(updAtomInfo(atomIndex1), updAtomInfo(atomIndex2))).setDefaultBondLength(distance);
        }
        return *this;
    }

    CompoundRep& matchDefaultBondAngles(const Compound::AtomTargetLocations& atomTargets) 
    {
        std::vector< AtomIndexList > atomTriples = getBondedAtomRuns(3, atomTargets);

        std::vector< AtomIndexList >::const_iterator bonds13Ix;
        for (bonds13Ix = atomTriples.begin(); bonds13Ix != atomTriples.end(); ++bonds13Ix) 
        {
            Compound::AtomIndex atomIndex1 = (*bonds13Ix)[0];
            Compound::AtomIndex atomIndex2 = (*bonds13Ix)[1];
            Compound::AtomIndex atomIndex3 = (*bonds13Ix)[2];

            // for efficiency, set each angle only once, not both 3->2->1 and 1->2->3
            if (atomIndex3 < atomIndex1) continue;

            UnitVec3 v1(atomTargets.find(atomIndex1)->second - atomTargets.find(atomIndex2)->second);
            UnitVec3 v2(atomTargets.find(atomIndex3)->second - atomTargets.find(atomIndex2)->second);

            Real dotProduct = dot(v1, v2);
            assert(dotProduct < 1.1);
            assert(dotProduct > -1.1);
            if (dotProduct > 1.0) dotProduct = 1.0;
            if (dotProduct < -1.0) dotProduct = -1.0;
            Real angle = std::acos(dotProduct);

            // std::cerr << angle / SimTK::Deg2Rad << std::endl;

            setDefaultBondAngle(angle, atomIndex1, atomIndex2, atomIndex3);
        }

        return *this;
    }

    // Helper method for matchDefaultDihedralAngles
    bool isPlanarBond(
            Compound::AtomIndex atomIndex2,
            Compound::AtomIndex atomIndex3) 
    {
        const CompoundAtom& atom2 = getAtom(atomIndex2);
        const CompoundAtom& atom3 = getAtom(atomIndex3);

        // Three criteria for whether bond is planar

        // 1) both central atoms have three bonds
        if (atom2.getNumBondCenters() != 3) return false;
        if (atom3.getNumBondCenters() != 3) return false;

        // 2) third bond center on each of those atoms is planar
        if (atom2.getBondCenter(CompoundAtom::BondCenterIndex(2)).getChirality() != BondCenter::Planar)
            return false;
        if (atom3.getBondCenter(CompoundAtom::BondCenterIndex(2)).getChirality() != BondCenter::Planar)
            return false;

        // 3) initial dihedral angle is near 0 or 180 degrees
        const BondCenter& bondCenter23 = 
            getBondCenter(getBondCenterInfo(getAtomInfo(atomIndex2), getAtomInfo(atomIndex3)).getIndex());
        Angle initialAngle = bondCenter23.getDefaultDihedralAngle();
        // normalize to be near zero
        while (initialAngle < -90.0 * Deg2Rad) initialAngle += 180.0 * Deg2Rad;
        while (initialAngle > 90.0 * Deg2Rad) initialAngle -= 180.0 * Deg2Rad;
        if (std::abs(initialAngle) > 0.001) return false;

        // If we got this far, it must be planar
        return true;
    }

    CompoundRep& matchDefaultDihedralAngles(
            const Compound::AtomTargetLocations& atomTargets, 
            Compound::PlanarBondMatchingPolicy policy) 
    {
        std::vector< AtomIndexList > atomQuads = getBondedAtomRuns(4, atomTargets);

        std::vector< AtomIndexList >::const_iterator bonds14Ix;
        for (bonds14Ix = atomQuads.begin(); bonds14Ix != atomQuads.end(); ++bonds14Ix) 
        {
            Compound::AtomIndex atomIndex1 = (*bonds14Ix)[0];
            Compound::AtomIndex atomIndex2 = (*bonds14Ix)[1];
            Compound::AtomIndex atomIndex3 = (*bonds14Ix)[2];
            Compound::AtomIndex atomIndex4 = (*bonds14Ix)[3];
            //std::cout<<__LINE__<<" "<<getAtom(atomIndex1)<<","<<getAtom(atomIndex2)<<","<<getAtom(atomIndex3)<<","<<getAtom(atomIndex4)<<std::endl; 
            // for efficiency, set each dihedral only once
            if (atomIndex4 < atomIndex1) continue;

			// Don't set dihedrals involving ring-closing bonds, as these can damage "real" dihedrals
			if ( getBond(atomIndex2, atomIndex1).isRingClosingBond() ) continue;
			if ( getBond(atomIndex3, atomIndex4).isRingClosingBond() ) continue;

            // Compute and set dihedral angle
            UnitVec3 bond12(atomTargets.find(atomIndex2)->second - atomTargets.find(atomIndex1)->second);
            UnitVec3 bond23(atomTargets.find(atomIndex3)->second - atomTargets.find(atomIndex2)->second);
            UnitVec3 bond34(atomTargets.find(atomIndex4)->second - atomTargets.find(atomIndex3)->second);
            //std::cout<<__FILE__<<":"<<__LINE__<<" "<<atomTargets.find(atomIndex1)->second << ","<< atomTargets.find(atomIndex2)->second<< ","<< atomTargets.find(atomIndex3)->second<< ","<< atomTargets.find(atomIndex4)->second<<std::endl;
            Angle angle = SimTK::calcDihedralAngle(bond12, bond23, bond34);

            // assert(false);  // need to implement general setDefaultDihedralAngle method

            // Don't set torsion for planar bonds, except maybe to Flip them
            if (policy == Compound::KeepPlanarBonds)
            {
                if ( isPlanarBond(atomIndex2, atomIndex3) )
                    continue;
            }

            if (policy == Compound::FlipPlanarBonds)
                if ( isPlanarBond(atomIndex2, atomIndex3) ) {
                    // TODO - decide whether to flip the dihedral angle 180 degrees
                    Angle initialAngle = calcDefaultDihedralAngle(                
                        atomIndex1, 
                        atomIndex2, 
                        atomIndex3, 
                        atomIndex4);
                    Angle diffAngle = angle - initialAngle;
                    // normalize to range (-180 degrees, 180 degrees)
                    while ( -SimTK::Pi >= diffAngle ) diffAngle += 2 * SimTK::Pi;
                    while ( SimTK::Pi < diffAngle )   diffAngle -= 2 * SimTK::Pi;
                    // Either flip the dihedral 180 degrees...
                    if (std::abs(diffAngle) > 0.5 * SimTK::Pi)
                        angle = initialAngle + SimTK::Pi;          
                    else // ... or do nothing.
                        continue;
                }

            setDefaultDihedralAngle( 
                angle, 
                atomIndex1, 
                atomIndex2, 
                atomIndex3, 
                atomIndex4);
        }

        return *this;
    }

    CompoundRep& matchDefaultTopLevelTransform(const Compound::AtomTargetLocations& atomTargets) 
    {
        Transform adjustment = getTransformAndResidual(atomTargets).transform;
        setTopLevelTransform( adjustment * getTopLevelTransform() );

        return *this;
    }

    TransformAndResidual getTransformAndResidual(const Compound::AtomTargetLocations& atomTargets) const
    {
        Kabsch78::VectorSet vecPairs;

        // Try for more efficient calculation of atom starting locations
        std::vector<Transform> atomSourceFrames(getNumAtoms());
        invalidateAtomFrameCache(atomSourceFrames, getNumAtoms());
        calcDefaultAtomFramesInCompoundFrame(atomSourceFrames);
        Real weight = 1.0;
        
        Compound::AtomTargetLocations::const_iterator tI;
        for (tI = atomTargets.begin(); tI != atomTargets.end(); ++tI) 
        {
            Compound::AtomIndex atomIndex = tI->first;
            const Vec3& target = tI->second;
            
            // slow
            // Vec3 source = calcDefaultAtomLocationInGroundFrame(getAtomName(atomIndex));
            
            // faster
            const Vec3 source = 
                    getTopLevelTransform() * atomSourceFrames[atomIndex].T();

            vecPairs.push_back(Vec3Pair(source, target, weight));
        }

        return Kabsch78::superpose(vecPairs);
    }

    const std::set<Compound::AtomName>& getAtomSynonyms(Compound::AtomIndex a) const
    {
        const AtomInfo& atomInfo = getAtomInfo(a);
        return atomInfo.getNames();
    }
    // scf added a new parameter .. when guessCoordinates is true, default atom positions from the biopolymer are pushed into the returned AtomTargetLocations.
    virtual Compound::AtomTargetLocations createAtomTargets(const PdbStructure& targetStructure, const bool guessCoordinates = false) const 
    {
        Compound::AtomTargetLocations answer;

        for (Compound::AtomIndex a(0); a < getNumAtoms(); ++a) 
        {
            int residueNumber = getPdbResidueNumber();
            char insertionCode = ' '; // TODO - make this an attribute of the residue?
            char chainId = getPdbChainId();

            String atomName = getAtomName(a);
            //std::cout<<__FILE__<<":"<<__LINE__<<" "<<a<<" "<<atomName<<" "<<residueNumber<<std::endl;
            // search synonyms if we cannot find this atom in the structure
            if (! targetStructure.hasAtom(atomName, PdbResidueId(residueNumber, insertionCode), chainId) )
            {
                const std::set<Compound::AtomName>& atomNames = getAtomSynonyms(a);
                std::set<Compound::AtomName>::const_iterator nameIx;
                for (nameIx = atomNames.begin(); nameIx != atomNames.end(); ++nameIx)
                {
                    atomName = *nameIx;
                    if ( targetStructure.hasAtom(atomName, PdbResidueId(residueNumber, insertionCode), chainId) )
                        break;
                }
            }

            if ( targetStructure.hasAtom(atomName, PdbResidueId(residueNumber, insertionCode), chainId) ) {
                const PdbAtom& pdbAtom = targetStructure.getAtom( atomName, PdbResidueId(residueNumber, insertionCode), chainId );
                if (pdbAtom.hasLocation())
                    answer[a] = pdbAtom.getLocation();
            }
            else {
                // std::cerr << atomName << std::endl;
            }
        }

        return answer;
    }

    virtual Compound::AtomTargetLocations createAtomTargets(const PdbChain& targetChain, const bool guessCoordinates = false) const 
    {
        Compound::AtomTargetLocations answer;

        for (Compound::AtomIndex a(0); a < getNumAtoms(); ++a) 
        {
            int residueNumber = getPdbResidueNumber();
            char insertionCode = ' '; // TODO - make this an attribute of the residue?
            char chainId = getPdbChainId();

            String atomName = getAtomName(a);
            // search synonyms if we cannot find this atom in the structure
            if (! targetChain.hasAtom(atomName, PdbResidueId(residueNumber, insertionCode)) )
            {
                const std::set<Compound::AtomName>& atomNames = getAtomSynonyms(a);
                std::set<Compound::AtomName>::const_iterator nameIx;
                for (nameIx = atomNames.begin(); nameIx != atomNames.end(); ++nameIx)
                {
                    atomName = *nameIx;
                    if ( targetChain.hasAtom(atomName, PdbResidueId(residueNumber, insertionCode)) )
                        break;
                }
            }

            if ( targetChain.hasAtom(atomName, PdbResidueId(residueNumber, insertionCode)) ) {
                const PdbAtom& pdbAtom = targetChain.getAtom( atomName, PdbResidueId(residueNumber, insertionCode) );
                if (pdbAtom.hasLocation())
                    answer[a] = pdbAtom.getLocation();
            }
            else {
                // std::cerr << atomName << std::endl;
            }
        }

        return answer;
    }

	// TODO - too much copy paste in these createAtomTargets methods
    virtual Compound::AtomTargetLocations createAtomTargets(const PdbResidue& targetResidue, const bool guessCoordinates = false) const 
    {
        Compound::AtomTargetLocations answer;
		if (getPdbResidueNumber() != targetResidue.getPdbResidueNumber()) return answer;

        for (Compound::AtomIndex a(0); a < getNumAtoms(); ++a) 
        {
            int residueNumber = getPdbResidueNumber();
            char insertionCode = ' '; // TODO - make this an attribute of the residue?
            char chainId = getPdbChainId();

            String atomName = getAtomName(a);
            // search synonyms if we cannot find this atom in the structure
            if (! targetResidue.hasAtom(atomName) )
            {
                const std::set<Compound::AtomName>& atomNames = getAtomSynonyms(a);
                std::set<Compound::AtomName>::const_iterator nameIx;
                for (nameIx = atomNames.begin(); nameIx != atomNames.end(); ++nameIx)
                {
                    atomName = *nameIx;
                    if ( targetResidue.hasAtom(atomName) )
                        break;
                }
            }

            if ( targetResidue.hasAtom(atomName) ) {
                const PdbAtom& pdbAtom = targetResidue.getAtom( atomName );
                if (pdbAtom.hasLocation())
                    answer[a] = pdbAtom.getLocation();
            }
            else { // skip atoms not found in pdbResidue
                // std::cerr << atomName << std::endl;
            }
        }

        return answer;
    }

    /// New way to do PDB writing: create intermediate PdbChain object
    /// Write current default(initial) Compound configuration into a PdbChain object
    virtual const CompoundRep& populateDefaultPdbChain(
        class PdbChain& pdbChain, 
        int& defaultNextResidueNumber,
        const Transform& transform) const 
    {
		Transform myTransform = transform;
            // if (!hasParentCompound()) {
                myTransform = myTransform * getTopLevelTransform();
            // }

        int residueNumber = getPdbResidueNumber();

        // try to guess when to use internal PdbResidueNumber vs. defaultNextResidueNumber
        if (-9999 > getPdbResidueNumber()) 
            residueNumber = defaultNextResidueNumber;

        // In case of residue number conflicts, find a new number
        if (pdbChain.hasResidue(PdbResidueId(residueNumber)))
            residueNumber = defaultNextResidueNumber;
        while (pdbChain.hasResidue(PdbResidueId(residueNumber)))
        {
            ++defaultNextResidueNumber;
            residueNumber = defaultNextResidueNumber;        
        }

        // In case of residue number conflicts, find a new number
        if (pdbChain.hasResidue(PdbResidueId(residueNumber)))
            residueNumber = defaultNextResidueNumber;
        while (pdbChain.hasResidue(PdbResidueId(residueNumber)))
        {
            ++defaultNextResidueNumber;
            residueNumber = defaultNextResidueNumber;        
        }

        pdbChain.appendResidue( PdbResidue(getOwnerHandle(), residueNumber, myTransform) );

        defaultNextResidueNumber = residueNumber + 1;

        return *this;
    }

    /// New way to do PDB writing: create intermediate PdbChain object
    /// Write current default(initial) Compound configuration into a PdbChain object
    virtual const CompoundRep& populatePdbChain(
        const State& state, 
        class PdbChain& pdbChain, 
        int& defaultNextResidueNumber,
        const Transform& transform) const 
    {

		Transform myTransform = transform;
		// Don't apply top-level transform for state-taking methods!!!
        //    if (!hasParentCompound()) {
        //        myTransform = myTransform * getTopLevelTransform();
        //    }

        int residueNumber = getPdbResidueNumber();

        // try to guess when to use internal PdbResidueNumber vs. defaultNextResidueNumber
        if (-9999 > getPdbResidueNumber()) 
            residueNumber = defaultNextResidueNumber;

        // In case of residue number conflicts, find a new number
        if (pdbChain.hasResidue(PdbResidueId(residueNumber)))
            residueNumber = defaultNextResidueNumber;
        while (pdbChain.hasResidue(PdbResidueId(residueNumber)))
        {
            ++defaultNextResidueNumber;
            residueNumber = defaultNextResidueNumber;        
        }

        pdbChain.appendResidue( PdbResidue(state, getOwnerHandle(), residueNumber, myTransform) );

        defaultNextResidueNumber = residueNumber + 1;

        return *this;
    }

    // One argument version of writeDefaultPdb begins numbering atoms at 1
    std::ostream& writeDefaultPdb(std::ostream& os, const Transform& transform) const;
    std::ostream& writeDefaultPdb(std::ostream& os, int& nextSerialNumber, const Transform& transform) const;

    //std::ostream& writeDefaultAtomPdb(
    //    const Compound::AtomName& name, 
    //    std::ostream& os, 
    //    int& nextSerialNumber,
    //    const Transform& transform
    //    ) const;

    std::ostream& writePdb(
        const State& state, 
        std::ostream& os, 
        const Transform& transform) const;

    std::ostream& writePdb(
        const State& state, 
        std::ostream& os, 
        int& nextSerialNumber,
        const Transform& transform) const;

    //std::ostream& writeAtomPdb(
    //    const State& state, 
    //    const Compound::AtomName&   name, 
    //    std::ostream& os, 
    //    int& nextSerialNumber, 
    //    const Transform& transform) const;

    //std::ostream& writeAtomPdb(
    //    const Compound::AtomName&   name, 
    //    std::ostream&               os, 
    //    int&                        nextSerialNumber,
    //    const Vec3&                 location
    //    ) const;
// protected:

    bool hasInboardBondCenter() const;

    CompoundRep& convertInboardBondCenterToOutboard();

    const BondCenter& getInboardBondCenter() const;
    BondCenter& updInboardBondCenter();
    const BondCenterInfo& getInboardBondCenterInfo() const;
    BondCenterInfo& updInboardBondCenterInfo();

    CompoundRep& setInboardBondCenter(const Compound::BondCenterName& n);
    CompoundRep& setInboardBondCenter(Compound::BondCenterIndex id);

    Compound::BondCenterIndex addLocalCompound(
        const Compound::Name& scName, 
        const Compound& subcompound,
        const Transform& location = Transform());

    // Copy atoms etc.
    // Returns new bond center index of absorbed inboard bond center
    Compound::BondCenterIndex absorbSubcompound(const Compound::Name& scName, const Compound& subcompound, bool isBase);

    const BondInfo& getBondInfo(Compound::BondIndex bi) const {
        return allBonds[bi];
    }
    BondInfo& updBondInfo(Compound::BondIndex bi) {
        return allBonds[bi];
    }

    const BondInfo& getBondInfo(const AtomInfo& a1, const AtomInfo& a2) const {
        std::pair<Compound::AtomIndex, Compound::AtomIndex> key(a1.getIndex(), a2.getIndex());
        Compound::BondIndex bi = bondIndicesByAtomIndexPair.find(key)->second;

        return allBonds[bi];
    }
    BondInfo& updBondInfo(const AtomInfo& a1, const AtomInfo& a2) {
        std::pair<Compound::AtomIndex, Compound::AtomIndex> key(a1.getIndex(), a2.getIndex());
        Compound::BondIndex bi = bondIndicesByAtomIndexPair.find(key)->second;

        return allBonds[bi];
    }

	const Bond& getBond(Compound::AtomIndex atom1, Compound::AtomIndex atom2) {
		return getBond( getBondInfo(getAtomInfo(atom1), getAtomInfo(atom2)) );
	}

    const Bond& getBond(const BondInfo& bondInfo) const {
        return bondInfo.getBond();
        //if ( bondInfo.isLocalBond() || bondInfo.isRingClosingBond() );
        //else {
        //    assert(false);
        //}
        //    assert(bondInfo.isSubcompoundBond());
        //    Compound::Index subcompoundId = bondInfo.getSubcompoundId();
        //    // const CompoundRep& scRep = bondInfo.getSubcompound().getImpl();
        //    const CompoundRep& scRep = getSubcompound(subcompoundId).getImpl();
        //    return scRep.getBond(scRep.getBondInfo(bondInfo.getSubcompoundBondIndex()));
        //}
    }
    Bond& updBond(BondInfo& bondInfo) {
        return bondInfo.updBond();
        //if ( bondInfo.isLocalBond() || bondInfo.isRingClosingBond() );
        //else {
        //    assert(false);
        //}
        //    assert(bondInfo.isSubcompoundBond());
        //    CompoundRep& scRep = updSubcompound(bondInfo.getSubcompoundId()).updImpl();
        //    // CompoundRep& scRep = bondInfo.updSubcompound().updImpl();
        //    return scRep.updBond(scRep.updBondInfo(bondInfo.getSubcompoundBondIndex()));
        //}
    }

    Transform calcDefaultBondCenterFrameInAtomFrame(const BondCenterInfo& info) const;
    const Transform calcDefaultBondCenterFrameInCompoundFrame(const Compound::BondCenterName name) const;

    const Transform calcDefaultBondCenterFrameInCompoundFrame(const BondCenterInfo& info, std::vector<Transform>& atomFrameCache) const;

    Compound::BondCenterIndex getBondCenterIndex(const Compound::BondCenterName& name) const;

    // return the bond center on atom1 that is attached to atom2
    const BondCenterInfo& getBondCenterInfo(const Compound::AtomName& atom1, const Compound::AtomName& atom2) const;
    BondCenterInfo& updBondCenterInfo(const Compound::AtomName& atom1, const Compound::AtomName& atom2);
    const BondCenterInfo& getBondCenterInfo(const AtomInfo& atom1, const AtomInfo& atom2) const;
    BondCenterInfo& updBondCenterInfo(const AtomInfo& atom1, const AtomInfo& atom2) ;
    BondCenterInfo& updBondCenterInfo(const Compound::BondCenterName&);
    const BondCenterInfo& getBondCenterInfo(const Compound::BondCenterName&) const;
    BondCenterInfo& updBondCenterInfo(Compound::BondCenterIndex);
    const BondCenterInfo& getBondCenterInfo(Compound::BondCenterIndex) const;
    BondCenterInfo& updBondCenterInfo(Compound::AtomIndex atomId, CompoundAtom::BondCenterIndex atomBondCenterIndex);
    const BondCenterInfo& getBondCenterInfo(Compound::AtomIndex atomId, CompoundAtom::BondCenterIndex atomBondCenterIndex) const;
    BondCenterInfo& updBondCenterInfo(BondCenterInfo::AtomKey key);
    const BondCenterInfo& getBondCenterInfo(BondCenterInfo::AtomKey key) const;

    bool hasBondCenter(const Compound::BondCenterName&) const;
    bool hasBondCenter(Compound::AtomIndex atomId, CompoundAtom::BondCenterIndex atomBondCenterIndex) const {
        return hasBondCenter(BondCenterInfo::AtomKey(atomId, atomBondCenterIndex));
    }
    bool hasBondCenter(Compound::BondCenterIndex id) const;
    bool hasBondCenter(const BondCenterInfo::AtomKey& key) const {
        return bondCenterIndicesByAtomKey.find(key) != bondCenterIndicesByAtomKey.end();
    }

    BondCenter& updBondCenter(const Compound::BondCenterName& name);
    const BondCenter& getBondCenter(const Compound::BondCenterName& name) const;
    BondCenter& updBondCenter(Compound::BondCenterIndex id);

    const BondCenter& getBondCenter(Compound::BondCenterIndex id) const;

    const BondCenter& getBondCenter(const BondCenterInfo& info) const;

    BondCenter& updBondCenter(const BondCenterInfo& info);

    Compound::AtomIndex getAtomIndex(const Compound::AtomName& atomName) const;

    AtomInfo& updAtomInfo(const Compound::AtomName&);

    const AtomInfo& getAtomInfo(const Compound::AtomName& name) const {
        // assert(CompoundPathName::isValidAtomName(name));
        assert(hasAtom(name));
        return getAtomInfo(atomIdsByName.find(name)->second);
    }

    AtomInfo& updAtomInfo(Compound::AtomIndex);
    const AtomInfo& getAtomInfo(Compound::AtomIndex id) const {
        assert((Compound::AtomIndex)allAtoms.size() > id);
        assert(0 <= id);
        return allAtoms[id];
    }
    //AtomInfo& updAtomInfo(Compound::Index subcompoundId, Compound::AtomIndex subAtomIndex) {
    //    const CompoundRep&      scRep           = getSubcompound(subcompoundId).getImpl();
    //    const AtomInfo&         scAtomInfo      = scRep.getAtomInfo(subAtomIndex);
    //    const Compound::AtomIndex  parentAtomIndex    = scAtomInfo.getParentCompoundAtomIndex();
    //    return updAtomInfo(parentAtomIndex);
    //}
    //const AtomInfo& getAtomInfo(Compound::Index subcompoundId, Compound::AtomIndex subAtomIndex) const {
    //    const CompoundRep&      scRep           = getSubcompound(subcompoundId).getImpl();
    //    const AtomInfo&         scAtomInfo      = scRep.getAtomInfo(subAtomIndex);
    //    const Compound::AtomIndex  parentAtomIndex    = scAtomInfo.getParentCompoundAtomIndex();
    //    return getAtomInfo(parentAtomIndex);
    //}

    bool hasAtom(const Compound::AtomName& name) const;

    bool hasAtom(Compound::AtomIndex atomId) const {
        if (atomId < 0) return false;
        if (atomId >= (Compound::AtomIndex) allAtoms.size()) return false;

        return true;
    }

    const CompoundAtom& getAtom(const Compound::AtomName& name) const {
        return getAtom(getAtomInfo(name));
    }
    CompoundAtom& updAtom(const Compound::AtomName& name) {
        return updAtom(updAtomInfo(name));
    }
    const CompoundAtom& getAtom(Compound::AtomIndex id) const {
        return getAtom(getAtomInfo(id));
    }
    CompoundAtom& updAtom(Compound::AtomIndex id) {
        return updAtom(updAtomInfo(id));
    }
    CompoundAtom& updAtom(AtomInfo& info) {
        return info.updAtom();
    }
    const CompoundAtom& getAtom(const AtomInfo& info) const {
        return info.getAtom();
    }

    Compound::BondIndex getNumBonds() const {return Compound::BondIndex(allBonds.size());}
    Compound::AtomIndex getBondAtomIndex(Compound::BondIndex bid, int which) const;

    //const CompoundInfo& getSubcompoundInfo(const Compound::Name& name) const 
    //{
    //    assert( hasSubcompound(name) );

    //    // TODO - parse "X/Y" indirect subcompound identifiers
    //    Compound::Index subcompoundId;

    //    // First check for simple subcompound name without any "/" separators
    //    if (CompoundPathName::isValidSubcompoundName(name)) 
    //    {
    //        subcompoundId = subcompoundIdsByName.find(name)->second;
    //    }
    //    else // parse "X/Y" path type subcompound names
    //    {
    //        std::vector<String> tokens;
    //        if (CompoundPathName::isValidSubcompoundPathName(name, &tokens)) 
    //        {
    //            String subcompoundName = tokens[0];
    //            if (! hasSubcompound(subcompoundName)) 
    //            {
    //                assert(false); // TODO - raise exception - hasSubcompound() check should have caught this
    //            }

    //            Compound::Index topSubcompoundId = getSubcompoundInfo(subcompoundName).getIndex();

    //            const CompoundRep& scRep = getSubcompound(subcompoundName).getImpl();
    //            Compound::Index childSubcompoundId = scRep.getSubcompoundInfo(CompoundPathName::shiftLeftPathName(name)).getIndex();

    //            // TODO find CompoundInfo that matches topSubcompoundId and childSubcompoundId
    //            // TODO this is not efficient, checking every subcompound
    //            std::vector<CompoundInfo>::const_iterator scI;
    //            for (scI = allSubcompounds.begin(); scI != allSubcompounds.end(); ++scI) 
    //            {
    //                if (scI->isLocal()) continue;
    //                if (scI->isBonded()) continue;
    //                if (scI->getIntermediateSubcompoundId() != topSubcompoundId) continue;
    //                if (scI->getIntermediateSubcompoundSubcompoundId() != childSubcompoundId) continue;

    //                // If we get this far, we have found the correct subcompound
    //                subcompoundId = scI->getIndex();
    //                break;
    //            }

    //            Compound::Index invalidCompoundId;
    //            assert(subcompoundId != invalidCompoundId);
    //        }
    //        else { // string not well formed
    //            assert(false);
    //            // TODO raise exception - hasSubcompound() check should have caught this
    //        }
    //    }

    //    return getSubcompoundInfo(subcompoundId);
    //}
    //CompoundInfo& updSubcompoundInfo(const Compound::Name& name) {
    //    assert( hasSubcompound(name) );

    //    const Compound::Index id = subcompoundIdsByName.find(name)->second;
    //    return updSubcompoundInfo(id);
    //}
    //const CompoundInfo& getSubcompoundInfo(Compound::Index id) const {
    //    assert (0 <= id);
    //    assert ((Compound::Index)allSubcompounds.size() > id);

    //    return allSubcompounds[id];
    //}
    //CompoundInfo& updSubcompoundInfo(Compound::Index id) {
    //    assert (0 <= id);
    //    assert ((Compound::Index)allSubcompounds.size() > id);

    //    return allSubcompounds[id];
    //}

    //Compound& updSubcompound(const Compound::Name& name) {
    //    CompoundInfo& info = updSubcompoundInfo(name);
    //    return updSubcompound(info);
    //}
    //const Compound& getSubcompound(const Compound::Name& name) const {
    //    const CompoundInfo& info = getSubcompoundInfo(name);
    //    return getSubcompound(info);
    //}
    //Compound& updSubcompound(Compound::Index id) {
    //    CompoundInfo& info = updSubcompoundInfo(id);
    //    return updSubcompound(info);
    //}
    //const Compound& getSubcompound(Compound::Index id) const {
    //    const CompoundInfo& info = getSubcompoundInfo(id);
    //    return getSubcompound(info);
    //}

    //const Compound& getSubcompound(const CompoundInfo& info) const {
    //    if (info.isLocal()) {
    //        return info.getCompound();
    //    }
    //    else if (info.isBonded()) {
    //        return info.getCompound();
    //    }
    //    else { // subcompound of subcompound
    //        const CompoundRep& sc1 = getSubcompound(info.getIntermediateSubcompoundId()).getImpl();
    //        return sc1.getSubcompound(info.getIntermediateSubcompoundSubcompoundId());
    //    }
    //}

    //Compound& updSubcompound(CompoundInfo& info) {
    //    if (info.isLocal()) {
    //        return info.updCompound();
    //    }
    //    else if (info.isBonded()) {
    //        return info.updCompound();
    //    }
    //    else { // subcompound of subcompound
    //        CompoundRep& sc1 = updSubcompound(info.getIntermediateSubcompoundId()).updImpl();
    //        return sc1.updSubcompound(info.getIntermediateSubcompoundSubcompoundId());
    //    }
    //}


    //// const Compound& getSubcompound(int) const;
    //// Compound& getSubcompound(int);

    //CompoundRep& nameSubcompound(Compound::Name newName, Compound::Name olderName) 
    //{
    //    assert(hasSubcompound(olderName));
    //    assert(!hasSubcompound(newName));
    //    Compound::Index subcompoundId = getSubcompoundInfo(olderName).getIndex();
    //    nameSubcompound(newName, subcompoundId);
    //    assert(hasSubcompound(newName));
    //    return *this;
    //}

    //CompoundRep& nameSubcompound(Compound::Name newName, Compound::Index subcompoundId) 
    //{
    //    assert(hasSubcompound(subcompoundId));
    //    assert(! hasSubcompound(newName));
    //    subcompoundIdsByName[newName] = subcompoundId;
    //    assert(hasSubcompound(newName));
    //    return *this;
    //}


    //bool hasSubcompound(const Compound::Name& name) const 
    //{
    //    // First check for simple subcompound name without any "/" separators
    //    if (CompoundPathName::isValidSubcompoundName(name)) {
    //        return subcompoundIdsByName.find(name) != subcompoundIdsByName.end();
    //    }
    //    else // parse "X/Y" path type subcompound names
    //    {
    //        std::vector<String> tokens;
    //        if (CompoundPathName::isValidSubcompoundPathName(name, &tokens)) 
    //        {
    //            String subcompoundName = tokens[0];
    //            if (! hasSubcompound(subcompoundName)) 
    //                return false;

    //            const CompoundRep& scRep = getSubcompound(subcompoundName).getImpl();
    //            return scRep.hasSubcompound(CompoundPathName::shiftLeftPathName(name));
    //        }
    //        else { // string not well formed
    //            return false;
    //        }
    //    }
    //}


    //bool hasSubcompound(const Compound::Index cId) const {
    //    if (cId < 0) return false;
    //    if (cId >= (Compound::Index)allSubcompounds.size()) return false;
    //    return true;
    //}

    static mdunits::Length getConsensusBondLength(const BondCenter& c1, const BondCenter& c2) 
    {
        mdunits::Length d1 = c1.getDefaultBondLength();
        mdunits::Length d2 = c2.getDefaultBondLength();

        mdunits::Length answer = d1;

        // Need to address the following cases:
        // 1) neither center has distance defined -> raise error
        // 2) one or other has distance defined -> OK, use that distance
        // 3) both have distance defined and agree -> OK, use that distance
        // 4) both have distance defined and disagree -> raise error

        // No information available to determine bond length
        if (isNaN(d1) && isNaN(d2)) assert(false); // case 1

        else if (isNaN(d1)) answer = d2; // case 2a
        else if (isNaN(d2)) answer = d1; // case 2b
        else if (d1 == d2) answer = d1; // case 3
        else assert(false); // case 4

        return answer;
    }

    static Angle getConsensusDihedralAngle(const BondCenter& c1, const BondCenter& c2) 
    {
        Angle a1 = c1.getDefaultDihedralAngle();
        Angle a2 = c2.getDefaultDihedralAngle();

        Angle answer = a1;

        // The same rules for dihedral angle as for bond length,
        // except that if both are NaN, default to 180 degrees
        if (isNaN(a1) && isNaN(a2)) answer = 180*Deg2Rad; // case 1
        else if (isNaN(a1)) answer = a2; // case 2a
        else if (isNaN(a2)) answer = a1; // case 2b
        else if (a1 == a2) answer = a1; // case 3
        else assert(false); // case 4

        return answer;
    }

    //Compound::BondIndex bondBondCenters(
    //    Compound::BondCenterIndex outboardId, 
    //    Compound::BondCenterIndex inboardId
    //    ) 
    //{
    //    mdunits::Length bondLength = getConsensusBondLength   (getBondCenter(outboardId), getBondCenter(inboardId));
    //    Angle    dihedral   = getConsensusDihedralAngle(getBondCenter(outboardId), getBondCenter(inboardId));

    //    return bondBondCenters(outboardId, inboardId, bondLength, dihedral);
    //}

    // TODO When a new bond is created, always call indexNewBond() to establish cross references
    //Compound::BondIndex bondBondCenters(
    //    Compound::BondCenterIndex outboardId, 
    //    Compound::BondCenterIndex inboardId,
    //    mdunits::Length               distance,
    //    Angle                  dihedral
    //    ) 

    void indexNewBond(const BondInfo& newBondInfo)
    {
        const Compound::BondCenterIndex outboardId = newBondInfo.getChildBondCenterIndex();
        const Compound::BondCenterIndex inboardId = newBondInfo.getParentBondCenterIndex();
        const Compound::BondIndex bondIndex = newBondInfo.getIndex();

        BondCenter&     outboardBondCenter      = updBondCenter(outboardId);
        BondCenter&     inboardBondCenter       = updBondCenter(inboardId);
        BondCenterInfo& outboardBondCenterInfo  = updBondCenterInfo(outboardId);
        BondCenterInfo& inboardBondCenterInfo   = updBondCenterInfo(inboardId);

        assert(! outboardBondCenter.isBonded() );
        assert(! inboardBondCenter.isBonded() );
        assert(! outboardBondCenterInfo.isBonded() );
        assert(! inboardBondCenterInfo.isBonded() );

        // Add a new BondInfo to this compound.
        // const Compound::BondIndex bondIndex = Compound::BondIndex(allBonds.size());

        // TODO - BondInfo constructor should depend on type of bond
        // allBonds.push_back( BondInfo(bondIndex, outboardId, inboardId, distance, dihedral) );

        // Mark the two newly-connected BondCenterInfos so they know they're connected.
        outboardBondCenterInfo.setBondPartnerBondCenterIndex(inboardBondCenterInfo.getIndex());
        inboardBondCenterInfo.setBondPartnerBondCenterIndex(outboardBondCenterInfo.getIndex());
        outboardBondCenterInfo.setBondIndex( bondIndex );
        inboardBondCenterInfo.setBondIndex( bondIndex );

        // Reach down to the Atoms and mark the physical bonds as in use, although they
        // don't know to whom they are connected.
        outboardBondCenter.setBonded(true);
        inboardBondCenter.setBonded(true);

        // Build a map so that we can find the connecting BondCenterInfo given the
        // AtomInfo indexes in either order.
        std::pair<Compound::AtomIndex, Compound::AtomIndex> 
            key1(outboardBondCenterInfo.getAtomIndex(), inboardBondCenterInfo.getAtomIndex());
        std::pair<Compound::AtomIndex, Compound::AtomIndex> 
            key2(inboardBondCenterInfo.getAtomIndex(), outboardBondCenterInfo.getAtomIndex());
        bondIndicesByAtomIndexPair[key1] = bondIndex;
        bondIndicesByAtomIndexPair[key2] = bondIndex;

        assert( outboardBondCenter.isBonded() );
        assert( inboardBondCenter.isBonded() );
        assert( outboardBondCenterInfo.isBonded() );
        assert( inboardBondCenterInfo.isBonded() );

        // return bondIndex;
    }

    //bool hasLocalSubcompound(const Compound::Name& name) const {
    //    return localSubcompoundIdsByName.find(name) != localSubcompoundIdsByName.end();
    //}

    //bool hasBondedSubcompound(const Compound::Name& name) const {
    //    return bondCenterIndexesByCompoundName.find(name) != bondCenterIndexesByCompoundName.end();
    //}

    CompoundRep& setBondMobility(BondMobility::Mobility mobility, const Compound::AtomName& atom1, const Compound::AtomName& atom2) 
    {
        AtomInfo& atomInfo1 = updAtomInfo(atom1);
        AtomInfo& atomInfo2 = updAtomInfo(atom2);
        BondInfo& bondInfo = updBondInfo(atomInfo1, atomInfo2);
        Bond& bond = updBond(bondInfo);

        bond.setMobility(mobility);
        // bond.setRotatable(isRotatable);
        
        return *this;
    }

    CompoundRep& setBondMobility(BondMobility::Mobility mobility, const Compound::BondIndex bondIndex) 
    {
        BondInfo& bondInfo = updBondInfo(bondIndex);
        Bond& bond = updBond(bondInfo);

        bond.setMobility(mobility);
        // bond.setRotatable(isRotatable);
        
        return *this;
    }


    CompoundRep& setPdbResidueNumber(int);
    int getPdbResidueNumber() const;

    CompoundRep& setPdbResidueName(const String&);
    const String& getPdbResidueName() const;

    CompoundRep& setPdbChainId(char);
    char getPdbChainId() const;


    // const Bond& getBond(const String& name) const;
    // Bond& getBond(const String& name);
    // const Bond& getBond(int id) const;
    // Bond& getBond(int id);
    // bool hasBond(const String& name) const;

     
    std::ostream& dumpCompoundRepToStream(std::ostream& o, int level=0) const;

    // bool hasParentCompound() const {return haveParentCompound;}


    CompoundRep& setTopLevelTransform(const Transform& transform) {
        // assert(!hasParentCompound());
        topLevelTransform = transform;

        return *this;
    }

    const Transform& getTopLevelTransform() const {
        // assert(!hasParentCompound());
        return topLevelTransform;
    }

protected:

    // offset + internal = nominal; offset + Bond = Dihedral
    Angle calcDefaultDihedralAngle(const DihedralAngle& dihedral) const
    {
        const BondCenterInfo& bc1 = getBondCenterInfo(dihedral.getBondCenter1Id());
        const BondCenterInfo& bc2 = getBondCenterInfo(dihedral.getBondCenter2Id());
        const AtomInfo& atom1 = getAtomInfo(bc1.getAtomIndex());
        const AtomInfo& atom2 = getAtomInfo(bc2.getAtomIndex());
        assert( atomsAreBonded(atom1, atom2) );
        const BondInfo& bondInfo = getBondInfo(atom1, atom2);
        const Bond& bond = getBond(bondInfo);

        Angle internalAngle = bond.getDefaultDihedralAngle();

        Angle internalOffset = calcDefaultInternalDihedralOffsetAngle(bc1.getIndex(), bc2.getIndex());

        Angle nominalAngle = internalAngle + internalOffset + dihedral.getNomenclatureOffset();

        return nominalAngle;
    }


    Angle calcDihedralAngle(const State& state, const String& dihedralName) const
    {
        assert( dihedralAnglesByName.find(dihedralName) != dihedralAnglesByName.end() );

        const DihedralAngle& dihedral = dihedralAnglesByName.find(dihedralName)->second;

        return calcDihedralAngle(state, dihedral);
    }

    Transform calcAtomFrameInGroundFrame(const State& state, Compound::AtomIndex atomId) const
    {
        const CompoundAtom& atom = getAtom(atomId);

        // Frame of parent body
        //DuMM::AtomIndex dummAtomIndex = atom.getDuMMAtomIndex();
        //const DuMMForceFieldSubsystem& dumm = ownerSystem->getMolecularMechanicsForceSubsystem();

        MobilizedBodyIndex bodyId = atom.getMobilizedBodyIndex();

        const SimbodyMatterSubsystem& matter = ownerSystem->getMatterSubsystem();
        const MobilizedBody& body = matter.getMobilizedBody(bodyId);
        const Transform& G_X_B = body.getBodyTransform(state);

        // Frame of atom in body;
        Transform B_X_A = atom.getFrameInMobilizedBodyFrame();

        return G_X_B * B_X_A;
    }
    
    Angle calcDihedralAngle(const State& state, const DihedralAngle& dihedral) const
    {
        const Bond& bond = getBondByDihedral(dihedral);
        MobilizedBodyIndex bodyId = bond.getPinJointId();
        if (bodyId.isValid())
        {
            assert(ownerSystem != NULL);
            const SimbodyMatterSubsystem& matter = ownerSystem->getMatterSubsystem();
            const MobilizedBody::Pin& body = (const MobilizedBody::Pin&)matter.getMobilizedBody(bodyId);
            Angle internalAngle = body.getAngle(state);

            // TODO - use simtime offset, not default

            Angle internalOffset = calcDefaultInternalDihedralOffsetAngle(dihedral.getBondCenter1Id(), dihedral.getBondCenter2Id());

            Angle nominalAngle = internalAngle + internalOffset + dihedral.getNomenclatureOffset();

            return nominalAngle;
        }
        else
        {
            assert(false); // TODO

            assert(ownerSystem != NULL);
            // Requires that system be realized to Position stage
            ownerSystem->realize(state, Stage::Position);

            const BondCenterInfo& bc21 = getBondCenterInfo(dihedral.getBondCenter1Id());
            const BondCenterInfo& bc34 = getBondCenterInfo(dihedral.getBondCenter2Id());

            // Find bond axis to project onto
            const AtomInfo& atom2 = getAtomInfo(bc21.getAtomIndex());
            const AtomInfo& atom3 = getAtomInfo(bc34.getAtomIndex());
            const BondCenterInfo& bondBondCenter = getBondCenterInfo(atom2, atom3);

            // Debug TODO - remove this stanza
            const Bond& bond = getBond(getBondInfo(bondBondCenter.getBondIndex()));
            MobilizedBodyIndex bodyId = bond.getPinJointId();
            std::cout << "Pin joint id = " << bodyId;
            std::cout << "\t";
            const SimbodyMatterSubsystem& matter = ownerSystem->getMatterSubsystem();
            const MobilizedBody::Pin& body = (const MobilizedBody::Pin&)matter.getMobilizedBody(bodyId);
            std::cout << "angle = " << body.getAngle(state) * DuMM::Rad2Deg << " degrees";
            std::cout << std::endl;

            assert(bondBondCenter.isBonded());
            assert(bondBondCenter.getIndex() != bc21.getIndex());
            assert(bondBondCenter.getIndex() != bc34.getIndex());
            assert(bc21.getIndex() != bc34.getIndex());

            UnitVec3 xAxis(1,0,0);

            // vector v1: from atom 1 to atom 2
            Transform G_X_A2 = calcDefaultAtomFrameInCompoundFrame(atom2.getIndex());
            Transform A2_X_BC21 = calcDefaultBondCenterFrameInAtomFrame(bc21);
            Transform G_X_BC21 = G_X_A2 * A2_X_BC21;
            UnitVec3 v1(G_X_BC21 * -xAxis); // negative x-axis because want 1->2, not 2->1 vector

            // vector v2: from atom 2 to atom 3
            Transform A2_X_BCB = calcDefaultBondCenterFrameInAtomFrame(bondBondCenter);
            Transform G_X_BCB = G_X_A2 * A2_X_BCB;
            UnitVec3 v2(G_X_BCB * xAxis);

            // vector v3: from atom 3 to atom 4
            Transform G_X_A3 = calcDefaultAtomFrameInCompoundFrame(atom3.getIndex());
            Transform A3_X_BC34 = calcDefaultBondCenterFrameInAtomFrame(bc34);
            Transform G_X_BC34 = G_X_A3 * A3_X_BC34;
            UnitVec3 v3(G_X_BC34 * xAxis);

            Angle nominalDihedralAngle = SimTK::calcDihedralAngle(v1, v2, v3) + dihedral.getNomenclatureOffset();

            return nominalDihedralAngle;
        }
    }

    // synonyms are alternate names for this compound type
    // synonyms should include the primary name of the compound
    // one use of synonyms is to help resolve biotypes atoms found in tinker parameter files
    std::set<Compound::Name> synonyms;

    // AtomInfo references
    std::vector<AtomInfo>                          allAtoms;    // [Compound::AtomIndex]
    // std::vector<CompoundAtom>                              localAtoms;  // [Compound::LocalAtomIndex]
    std::map<Compound::AtomName, Compound::AtomIndex> atomIdsByName;

    // bool haveParentCompound;

private:
    friend class Compound;
    friend class Bond;

    // ownerSystem is being used in two ways:
    // 1) ownerSystem plus ixWithinOwnerSystem represent handle for compounds directly owned by a CompoundSystem
    // 2) ownerSystem with invalid ixWithinOwnerSystem represents a handle to the system for subcompounds of
    //    a compound that is in turn directly owned by a CompoundSystem
    MultibodySystem* ownerSystem;
    // CompoundSystem* ownerSystem;
    // Compound::Index          ixWithinOwnerSystem;

    // Global transform that applies only to top level compound
    Transform topLevelTransform;

    Compound::Name name; // set on construction; means whatever you like

    // The following comment may be wrong cmb Feb 2009
    // local subcompounds placed directly - NOT those placed by bonds
    // std::vector<CompoundInfo> allSubcompounds; // [Compound::Index]
    // std::vector<Compound> localSubcompounds;
    // std::map<Compound::Name, Compound::Index> subcompoundIdsByName;
    // std::map<int, Transform> localSubcompoundTransformsById;

    // BondCenters
    std::vector<BondCenterInfo>              allBondCenters; // [Compound::BondCenterIndex]
    std::map<String, Compound::BondCenterIndex> bondCenterIndicesByName;
    BondCenterInfo::AtomKeyMap               bondCenterIndicesByAtomKey;

    // Bonds
    // bonds do not contain subcompounds
    std::vector<BondInfo> allBonds; // [Compound::BondIndex]
    std::map< std::pair<Compound::AtomIndex, Compound::AtomIndex>, Compound::BondIndex > 
                          bondIndicesByAtomIndexPair;

    // Dihedral Angles
    std::map<String, DihedralAngle> dihedralAnglesByName;

    int    pdbResidueNumber;
    String pdbResidueName;
    char   pdbChainId;
    
    class MemberForDebuggingCopyCtor {
    public:
        MemberForDebuggingCopyCtor() {}
        MemberForDebuggingCopyCtor(const MemberForDebuggingCopyCtor&) {
            // Put a breakpoint here to notice when CompoundRep copyCtor is called
            int x = 5;
        }
    };
    // MemberForDebuggingCopyCtor testMember;
};


BiotypeIndex SimTK_MOLMODEL_EXPORT getBiotypeIndex(
                        const Compound::Name& resName, 
                        const Compound::AtomName& atomName, 
                        Ordinality::Residue ordinality = Ordinality::Any);

// Return a biotype index matching any of the residue/atom names supplied
// Returns invalid index if not found
BiotypeIndex SimTK_MOLMODEL_EXPORT getBiotypeIndex(
                        const std::set<Compound::Name>& resNames, 
                        const std::set<Compound::AtomName>& atomNames, 
                        Ordinality::Residue ordinality = Ordinality::Any);

class BiopolymerResidueRep : public CompoundRep {
public:
    /*virtual*/ ~BiopolymerResidueRep() { }
    /*virtual*/ BiopolymerResidueRep* clone() const {return new BiopolymerResidueRep(*this);}

    BiopolymerResidueRep(String name, String tlc = "Unk", char olc = '?')
        : residueName(name), threeLetterCode(tlc), oneLetterCode(olc)
        {}

    BiopolymerResidueRep& setOneLetterCode(char olc) {
        oneLetterCode = olc;
        return *this;
    }
    BiopolymerResidueRep& setThreeLetterCode(const String& tlc) {
        threeLetterCode = tlc;
        return *this;
    }
    BiopolymerResidueRep& setResidueTypeName(const String& name) {
        residueName = name;
        return *this;
    }

    const String& getResidueTypeName() const {return residueName;}
    const String& getThreeLetterCode() const {return threeLetterCode;}
    char getOneLetterCode() const {return oneLetterCode;}

    // Attempt to deduce correct biotypes from global biotypes database
    /// @return true if all atoms have a valid biotype assigned, false otherwise
    bool assignBiotypes(Ordinality::Residue ordinality = Ordinality::Any) 
    {
        bool answer = true; // start optimistic

        // for each named atom, look up resname, atomname, ordinality
        // if biotype is still undefined, raise exception
        std::vector<AtomInfo>::iterator atomI;
        for (atomI = allAtoms.begin(); atomI != allAtoms.end(); ++atomI) 
        {
            // debugging
            // int atomIndex = atomI->getIndex();
            // std::cout << "biotype for atom " << atomI->getIndex() << std::endl;
            
            CompoundAtom& atom = updAtom(*atomI);
            
			// if the atom already has a valid biotype, keep it.
			if (atom.getBiotypeIndex().isValid()) continue;

            // Examine all possible residue names
            const std::set<Compound::Name>& residueNames = synonyms;

            // Create a container to hold variations of atom name
            const std::set<Compound::AtomName>& atomNames = getAtomSynonyms(atomI->getIndex());

            // Loop over residue names and atom names until a match is found
            bool foundBiotype = false;
            BiotypeIndex index = getBiotypeIndex(residueNames, atomNames, ordinality);
            if (index.isValid()) {
                atom.setBiotypeIndex(index);
                foundBiotype = true;
            }
            else {
                foundBiotype = false;
            }

            if (!foundBiotype) {
                answer = false;
                //std::cout<<__FILE__<<":"<<__LINE__<<" residueNames = "<<residueNames<<std::endl;
                std::set<Compound::AtomName>::iterator atomNamesIterator; 
                for (atomNamesIterator = atomNames.begin(); atomNamesIterator != atomNames.end(); atomNamesIterator++){
                    std::cout<<__FILE__<<":"<<__LINE__<<" atomNames = "<<string(*atomNamesIterator)<<std::endl;}
            }
            // Perhaps the atom already had a usable biotype...
            /// assert(atom.getBiotypeIndex().isValid());
            assert(foundBiotype);

        }

        return answer;
    }

private:
    String residueName;
    String threeLetterCode;
    char   oneLetterCode;
};



class BiopolymerRep : public CompoundRep {
public:
    friend class Biopolymer;

    BiopolymerRep* clone() const {return new BiopolymerRep(*this);}
    //const std::vector<String>& getResidueNames() const {
    //    return residueNames;
    //}

    //std::vector<String>& updResidueNames() {
    //    return residueNames;
    //}


    int getNumResidues() const;
    const String& getResidueName(int residueIndex) const;

    const ResidueInfo& getResidue(ResidueInfo::Index residueIndex) const;
    ResidueInfo& updResidue(ResidueInfo::Index residueIndex);

    const ResidueInfo& getResidue(const Compound::Name& residueName) const;
    ResidueInfo& updResidue(const Compound::Name& residueName);

    Transform calcDefaultResidueFrameInBiopolymerFrame(ResidueInfo::Index r, const std::vector<Transform>& atomFrameCache) const {
        // frame of inboard atom
        Compound::AtomIndex inboardAtomIndex = getResidue(r).getAtomIndex(ResidueInfo::AtomIndex(0));
        return atomFrameCache[inboardAtomIndex];
        // return calcDefaultAtomFrameInCompoundFrame(inboardAtomIndex);
    }

    virtual const CompoundRep& populateResidueDefaultPdbChain(
        ResidueInfo::Index r,
        class PdbChain& pdbChain, 
        int& defaultNextResidueNumber,
        const Transform& transform,
        const std::vector<Transform>& atomFrameCache) const 
    {
        Transform myTransform = transform * calcDefaultResidueFrameInBiopolymerFrame(r, atomFrameCache);
        const ResidueInfo& residue = getResidue(r);

        int residueNumber = residue.getPdbResidueNumber();
        if (-9999 > residueNumber) residueNumber = defaultNextResidueNumber;
        // In case of residue number conflicts, find a new number
        if (pdbChain.hasResidue(PdbResidueId(residueNumber)))
            residueNumber = defaultNextResidueNumber;
        while (pdbChain.hasResidue(PdbResidueId(residueNumber)))
            ++residueNumber;

        PdbResidue pdbResidue(residue.getPdbResidueName(), PdbResidueId(residueNumber));
        for (ResidueInfo::AtomIndex a(0); a < residue.getNumAtoms(); ++a) 
        {
            Compound::AtomIndex atomIx = residue.getAtomIndex(a);
            const CompoundAtom& atom = getAtom(atomIx);

            PdbAtom pdbAtom(residue.getAtomName(a), getAtomElement(atomIx));
            pdbAtom.setLocation(PdbAtomLocation(transform * atomFrameCache[atomIx].p()));
            pdbResidue.addAtom(pdbAtom);
        }
        pdbChain.appendResidue(pdbResidue);

        defaultNextResidueNumber = residueNumber + 1;

        return *this;
    }

    /// New way to do PDB writing: create intermediate PdbChain object
    /// Write current default(initial) Compound configuration into a PdbChain object
    virtual const CompoundRep& populateDefaultPdbChain(
        class PdbChain& pdbChain, 
        int& defaultNextResidueNumber,
        const Transform& transform) const 
    {
        std::vector<Transform> atomFrameCache(getNumAtoms());
        invalidateAtomFrameCache(atomFrameCache, getNumAtoms());
        calcDefaultAtomFramesInCompoundFrame(atomFrameCache);

        for (ResidueInfo::Index r(0); r < getNumResidues(); ++r) 
            populateResidueDefaultPdbChain(r, pdbChain, defaultNextResidueNumber, transform, atomFrameCache);

        return *this;
    }

    /// Write current default(initial) Compound configuration into a PdbChain object
    virtual const CompoundRep& populatePdbChain(
        const State& state, 
        class PdbChain& pdbChain, 
        int& defaultNextResidueNumber,
        const Transform& transform) const 
    {
        for (ResidueInfo::Index r(0); r < getNumResidues(); ++r) 
            populateResiduePdbChain(state, r, pdbChain, defaultNextResidueNumber, transform);

        return *this;
    }

    const CompoundRep& populateResiduePdbChain(
        const State& state,
        ResidueInfo::Index r,
        class PdbChain& pdbChain, 
        int& defaultNextResidueNumber,
        const Transform& transform) const 
    {
        // Transform myTransform = transform * calcResidueFrameInBiopolymerFrame(r);
        const ResidueInfo& residue = getResidue(r);

        int residueNumber = residue.getPdbResidueNumber();
        int insertionCode = residue.getPdbInsertionCode();
        if (-9999 > residueNumber) residueNumber = defaultNextResidueNumber;
        // In case of residue number conflicts, find a new number
        if (pdbChain.hasResidue(PdbResidueId(residueNumber,insertionCode)))
            residueNumber = defaultNextResidueNumber;
        while (pdbChain.hasResidue(PdbResidueId(residueNumber,insertionCode)))
            ++residueNumber;

        PdbResidue pdbResidue(residue.getPdbResidueName(), PdbResidueId(residueNumber,insertionCode ));
        for (ResidueInfo::AtomIndex a(0); a < residue.getNumAtoms(); ++a) 
        {
            Compound::AtomIndex atomIx = residue.getAtomIndex(a);
            const CompoundAtom& atom = getAtom(atomIx);

            PdbAtom pdbAtom(residue.getAtomName(a), getAtomElement(atomIx));
            pdbAtom.setLocation(
                PdbAtomLocation(transform * calcAtomLocationInGroundFrame(state, atomIx))
            );
            pdbResidue.addAtom(pdbAtom);
        }
        pdbChain.appendResidue(pdbResidue);

        defaultNextResidueNumber = residueNumber + 1;

        return *this;
    }


    //virtual std::ostream& writeDefaultPdb(
    //    std::ostream& os, 
    //    int& nextSerialNumber, 
    //    const Transform& transform) const;

    //virtual std::ostream& writePdb(
    //    const State& state, 
    //    std::ostream& os, 
    //    int& nextSerialNumber,
    //    const Transform& transform) const;

private:
    // std::vector<String> residueNames;
    std::vector<ResidueInfo> residues;
    std::map<Compound::Name, ResidueInfo::Index> residueIdsByName;
};

} // namespace SimTK

#endif // SimTK_COMPOUNDREP_H_
