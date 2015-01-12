/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Molmodel(tm)                         *
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

#include "SimTKsimbody.h"
#include "molmodel/internal/common.h"

#define DO_INSTANTIATE_COMPOUND_PIMPL_HANDLE
#include "molmodel/internal/Compound.h"
#include "molmodel/internal/Protein.h"
#include "molmodel/internal/RNA.h"
#include "molmodel/internal/DNA.h"
#undef DO_INSTANTIATE_COMPOUND_PIMPL_HANDLE

#include "CompoundRep.h"
#include "SimTKcommon/internal/PrivateImplementation_Defs.h"

#include "molmodel/internal/CompoundSystem.h"

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>
#include <cctype> // std::toupper
#include <stdexcept>
#include <set>

using namespace std;


namespace SimTK {


BiotypeIndex SimTK_MOLMODEL_EXPORT getBiotypeIndex(
                        const Compound::Name& resName, 
                        const Compound::AtomName& atomName, 
                        Ordinality::Residue ordinality) 
{
    if ( Biotype::exists(resName, atomName, ordinality) ) 
        return Biotype::get(resName, atomName, ordinality).getIndex();
    else
        return BiotypeIndex(); // invalid index
}

// Return a biotype index matching any of the residue/atom names supplied
// Returns invalid index if not found
BiotypeIndex SimTK_MOLMODEL_EXPORT getBiotypeIndex(
                        const std::set<Compound::Name>& resNames, 
                        const std::set<Compound::AtomName>& atomNames, 
                        Ordinality::Residue ordinality) 
{
    std::set<Compound::Name>::const_iterator resI;
    std::set<Compound::AtomName>::const_iterator atomI;

    // Create a modified list of atom names, with leading and trailing digits removed
    // Sometimes atom name has trailing digit
    // e.g. phenylalanine "CD1" should be biotype "CD"
    // but leucine "CG1" should be biotype "CG1"
    // so for generality try both ways
    // i.e. create a version of the atom name with trailing digits removed
    // also loop over atom synonyms
    std::set<Compound::AtomName> moreAtomNames;
    for (atomI = atomNames.begin(); atomI != atomNames.end(); ++atomI)
    {
        String atomName = *atomI;

        // Atom type never has leading digit found on hydrogen instances
        int startPos = atomName.find_first_not_of("0123456789");
        if (startPos != 0)
            atomName = atomName.substr(startPos);

        moreAtomNames.insert(atomName); // First atom name is the standard one

        Compound::Name shortAtomName = atomName;
        int endPos = atomName.find_last_not_of("0123456789");
        if ( endPos != (atomName.length() - 1) ) {
            shortAtomName = atomName.substr(0, endPos + 1);
            moreAtomNames.insert(shortAtomName);
        }
    }

    for (resI = resNames.begin(); resI != resNames.end(); ++resI)
        for (atomI = moreAtomNames.begin(); atomI != moreAtomNames.end(); ++atomI)
        {
            BiotypeIndex index = getBiotypeIndex(*resI, *atomI, ordinality);
            if (index.isValid()) 
                return index;
        }

    return BiotypeIndex(); // invalid index
}


// 1 Ring closing bond
BondInfo::BondInfo(
    Compound::BondIndex bId,
    Compound::BondCenterIndex bondCenter1,
    Compound::BondCenterIndex bondCenter2,
    const Bond& b
    )
    : id(bId), 
        // amLocalBond(false), 
        // amRingClosingBond(b.isRingClosingBond()), 
        parentBondCenterIndex(bondCenter1), 
        childBondCenterIndex(bondCenter2),
        bond(b)
{
    // bond.setRingClosingBond(true);
}

//// 2 Local subcompound bond
//// (but not to a "local subcompound")
//// sorry that's confusing...
//BondInfo::BondInfo(
//    Compound::BondIndex bId,
//    const Compound& c, 
//    Compound::BondCenterIndex bondCenter1,
//    Compound::BondCenterIndex bondCenter2,
//    mdunits::Length l, 
//    Angle t
//     )
//    : id(bId), 
//        // amLocalBond(true), 
//        amRingClosingBond(false), 
//        localBondSubcompound(c),
//        parentBondCenterIndex(bondCenter1), childBondCenterIndex(bondCenter2),
//        bond(l, t)
//{}

// 3 Bond inherited from a subcompound
//BondInfo::BondInfo(
//    Compound::BondIndex bId,
//    Compound::Index cId,
//    Compound::BondIndex scBId,
//    Compound::BondCenterIndex bondCenter1,
//    Compound::BondCenterIndex bondCenter2
//    //mdunits::Length l, 
//    //Angle t
//    )
//    : id(bId), 
//        //defaultLength(l), defaultDihedral(t), isRotatable(true),
//        //amRingClosingBond(false), 
//        amSubcompoundBond(true), amLocalBond(false), amRingClosingBond(false), 
//        subcompoundId(cId), subcompoundBondIndex(scBId),
//        parentBondCenterIndex(bondCenter1), childBondCenterIndex(bondCenter2)
//{}


BondCenter::BondCenter() 
    : 
    bonded(false), 
    inboard(false), 
    chirality(BondCenter::Planar),
    defaultBondLength(NaN), 
    defaultDihedralAngle(NaN)
{}


BondCenter::BondCenter(Angle angle1, Angle angle2, int yCenter, BondCenter::Chirality c)
    : 
    defaultBond1Angle(angle1), 
    defaultBond2Angle(angle2),
    defaultDihedralReferenceCenter(yCenter), 
    chirality(c), 
    defaultBondLength(NaN), 
    defaultDihedralAngle(NaN), 
    bonded(false),
    inboard(false)
{}

bool BondCenter::isBonded() const { 
    return (bonded);
}

BondCenter& BondCenter::setBonded(bool b) {
    bonded = b;
    return *this;
}


//BondCenter& BondCenter::setBonded(bool b) {
//    bonded = b;
//
//    return *this;
//}



//////////////////////
/// BondCenterInfo ///
//////////////////////

class CompoundRep;

// const Compound::BondCenterName& BondCenterInfo::getName() const {return name;}



////////////
/// Atom ///
////////////
    
CompoundAtom::CompoundAtom(const Element& element, Transform location)
  : 
    element(element), 
    biotypeIx(InvalidBiotypeIndex), 
    localTransform(location)
{
}

const Element& CompoundAtom::getElement() const {
    return element;
}



////////////////
/// AtomInfo ///
////////////////

// One constructor for local atoms
AtomInfo::AtomInfo(Compound::AtomIndex i, const CompoundAtom& atom, bool isBase /*, const Transform& xform */ )
  : id(i), 
    atom(atom), 
    bIsBaseAtom(isBase) // ,
    // frameInCompound(xform)
{}

///////////////////
/// CompoundRep ///
///////////////////

template class PIMPLImplementation<Compound,CompoundRep>;

// Add one simple atom unconnected to anything else
CompoundRep& CompoundRep::setBaseAtom(
    const Compound::AtomName& name, 
    const Element& element,
    const Transform& location) 
{
    // Add an AtomInfo reference in the list of all atoms
    const Compound::AtomIndex compoundAtomIndex = Compound::AtomIndex(allAtoms.size());
    allAtoms.push_back(AtomInfo(compoundAtomIndex, CompoundAtom(element, location), true));

    // Set name
    nameAtom(name, compoundAtomIndex);

    return *this;
}

CompoundRep& CompoundRep::setBaseAtom(
    const Compound::AtomName& name, 
    const Biotype& biotype,
    const Transform& location) 
{
    const Element& element = biotype.getElement();

    setBaseAtom(name, element, location);

    setBiotypeIndex(name, biotype.getIndex());

    return *this;
}

// Add a subcompound containing exactly one atom, so the Compound::AtomName can be reused for the Compound::Name
// This atom is not connected to anything else
CompoundRep& CompoundRep::setBaseAtom(
    const Compound::SingleAtom& compound,
    const Transform& location)
{
    // Only top level compounds can construct new topology
    // assert(! hasParentCompound() );
    // assert(! compound.getImpl().hasParentCompound() );

    Compound::AtomName atomName = compound.getAtomName(Compound::AtomIndex(0));
    setBaseCompound(atomName, compound, location);
    inheritAtomNames(atomName);

    // assert(getSubcompound(atomName).getImpl().hasParentCompound());

    return *this;
}


// Add a subcompound without attaching it to anything
Compound::BondCenterIndex CompoundRep::setBaseCompound(
    const Compound::Name& n, 
    const Compound& c,
    const Transform& location) 
{
    Compound::BondCenterIndex inboardIndex = addLocalCompound(n, c, location);

    // Inherit inboard bond, if available
    if (inboardIndex.isValid()) setInboardBondCenter(inboardIndex);

    return inboardIndex;
}

// Add a subcompound containing exactly one atom, so the Compound::AtomName can be reused for the Compound::Name
// This atom is connected to existing material
CompoundRep& CompoundRep::bondAtom(
    const Compound::SingleAtom&   compound, 
    const Compound::BondCenterPathName& parentBondName, 
    mdunits::Length                      distance,
    Angle                         dihedral,
    BondMobility::Mobility        mobility
    ) 
{
    // Only top level compounds can construct new topology
    // assert(! hasParentCompound() );
    // assert(! compound.getImpl().hasParentCompound() );

    const Compound::AtomName atomName = compound.getAtomName(Compound::AtomIndex(0));
    bondCompound(atomName, compound, parentBondName, distance, dihedral, mobility);
    inheritAtomNames(atomName);

    return *this;
}

// Compute atom location in local compound frame
Transform CompoundRep::calcDefaultAtomFrameInCompoundFrame(const Compound::AtomName& name) const {
    assert(hasAtom(name));

    const Compound::AtomIndex atomId = getAtomInfo(name).getIndex();
    return calcDefaultAtomFrameInCompoundFrame(atomId);
}

// Compute atom location in local compound frame
Vec3 CompoundRep::calcDefaultAtomLocationInCompoundFrame(const Compound::AtomName& name) const 
{
    return calcDefaultAtomFrameInCompoundFrame(name).p();
}

// Compute atom location in local compound frame
Transform CompoundRep::calcDefaultAtomFrameInGroundFrame(const Compound::AtomName& name) const 
{
    // TODO - this only works for top level compounds
    // assert(! hasParentCompound() );

    return getTopLevelTransform() * calcDefaultAtomFrameInCompoundFrame(name);
}

// Compute atom location in local compound frame
Vec3 CompoundRep::calcDefaultAtomLocationInGroundFrame(const Compound::AtomName& name) const 
{
    return calcDefaultAtomFrameInGroundFrame(name).p();
}

// Add a subcompound attached by a bond to an existing atom
// bondCompound("H1", MonovalentAtom(Element::Hydrogen()), "bond", "C/bond2", C_Hdistance );
CompoundRep& CompoundRep::bondCompound(
    const Compound::Name&           name, 
    const Compound&                 subcompoundArg, 
    const Compound::BondCenterPathName&   parentBondName, 
    mdunits::Length                        distance,
    Angle                           dihedral,
    BondMobility::Mobility          mobility
    ) 
{
    assert(! isNaN(dihedral) );

    const Compound::BondCenterIndex inboardBondCenterIndex = 
        absorbSubcompound(name, subcompoundArg, false);

    // Get atoms to bond
    const Compound::BondCenterIndex outboardBondCenterIndex = getBondCenterInfo(parentBondName).getIndex();

    // Don't bond this compound's official inboard bond center as outboard
    if (hasInboardBondCenter()) {
        const Compound::BondCenterIndex primaryInboardId = getInboardBondCenterInfo().getIndex();
        if (primaryInboardId == outboardBondCenterIndex) {
            convertInboardBondCenterToOutboard();
            // either that, or raise an exception...            
        }
    }
    assert( ! getBondCenter(outboardBondCenterIndex).isInboard() );

    // update bond info using subcompound
    Compound::BondIndex bondIndex(allBonds.size());
    allBonds.push_back(BondInfo(bondIndex, outboardBondCenterIndex, inboardBondCenterIndex, Bond(distance, dihedral, false)));
    indexNewBond(updBondInfo(bondIndex));

    //const Compound::BondIndex bondIndex = 
    //    bondBondCenters(outboardBondCenterIndex, inboardBondCenterIndex, distance, dihedral);

    const BondInfo& bondInfo = getBondInfo(bondIndex);

    // Set bond mobility
    Bond& bond = updBond(updBondInfo(bondIndex));
    bond.setMobility(mobility);

    return *this;
}


// Shorter version uses default bond length and dihedral angle
CompoundRep& CompoundRep::bondCompound(
    const Compound::Name&         n, 
    const Compound&               c, 
    const Compound::BondCenterPathName& parentBondName)
{
    // Only top level compounds can construct new topology
    // assert(! hasParentCompound() );
    // assert(! c.getImpl().hasParentCompound() );

    // Get distance and dihedral from child inboard bond center
    const CompoundRep& scRep = c.getImpl();
    const BondCenter& childCenter = scRep.getInboardBondCenter();
    const BondCenter& parentCenter = getBondCenter(parentBondName);

    mdunits::Length d = getConsensusBondLength(parentCenter, childCenter);
    Angle    a = getConsensusDihedralAngle(parentCenter, childCenter);

    bondCompound(n, c, parentBondName, d, a);

    return *this;
}


//Sam added bond mobility parameter
CompoundRep& CompoundRep::bondCompound(
    const Compound::Name&         n,
    const Compound&               c,
    const Compound::BondCenterPathName& parentBondName,
    BondMobility::Mobility mobility)
{
    // Only top level compounds can construct new topology
    // assert(! hasParentCompound() );
    // assert(! c.getImpl().hasParentCompound() );

    // Get distance and dihedral from child inboard bond center
    const CompoundRep& scRep = c.getImpl();
    const BondCenter& childCenter = scRep.getInboardBondCenter();
    const BondCenter& parentCenter = getBondCenter(parentBondName);

    mdunits::Length d = getConsensusBondLength(parentCenter, childCenter);
    Angle    a = getConsensusDihedralAngle(parentCenter, childCenter);
    bondCompound(n, c, parentBondName, d, a,mobility);

    return *this;
}



//CompoundRep& CompoundRep::removeAtom(const Compound::AtomName& name) {
//    removeAtom(getAtomInfo(name).getIndex());
//}
//
//
//CompoundRep& CompoundRep::removeAtom(const Compound::AtomIndex& atomId) {
//    AtomInfo& atomInfo = getAtomInfo(atomId);
//    Atom& atom = getAtom(atomId);
//
//    // TODO remove bonds and bond centers
//    // remove bond centers
//    std::vector<BondCenterInfo>::iterator bc;
//    for (bc = bondCenters.begin(); bc != bondCenters.end(); ++bc) {
//        if (bc->getAtomIndex() == atomId)
//            removeBondCenter(bc->getIndex());
//    }
//
//    // TODO remove bonds
//
//    // TODO remove name table entries
//    // TODO what if our parent has a reference to this atom?
//
//    // Finally, clear the data structures referencing the atoms
//    if (atomInfo.isLocal()) {
//        localAtoms[atomInfo.getLocalAtomIndex()] = Atom();
//    }
//    else {
//        Compound& subcompound = getSubcompound(atomInfo.getSubcompoundName());
//        subcompound.removeAtom(atomInfo.getSubcompoundAtomName());
//    }
//}


//CompoundRep& CompoundRep::setInboardBondCenter(
//    const Compound::BondCenterName& centerName, 
//    const Compound::AtomName& atomName, 
//    Angle zRotation,
//    Angle oldXRotation)
//{
//    addBondCenter(centerName, atomName, zRotation, oldXRotation);
//    setInboardBondCenter(centerName);
//
//    return *this;
//}

CompoundRep& CompoundRep::addFirstBondCenter(
    const Compound::BondCenterName& centerName, 
    const Compound::AtomName& atomName
    ) 
{
    // Only top level compounds can construct new topology
    // assert(! hasParentCompound() );

    assert (!hasBondCenter(centerName));

    const Compound::AtomIndex atomId = getAtomInfo(atomName).getIndex();
    const CompoundAtom::BondCenterIndex centerIndex = updAtom(atomId).addFirstBondCenter();

    addBondCenterInfo(atomId, centerIndex);
    bondCenterIndicesByName[centerName] = getBondCenterInfo(atomId, centerIndex).getIndex();

    assert (hasBondCenter(centerName));

    return *this;
}

CompoundRep& CompoundRep::addSecondBondCenter(
    const Compound::BondCenterName& centerName, 
    const Compound::AtomName& atomName,
    Angle bondAngle1
    ) 
{
    // Only top level compounds can construct new topology
    // assert(! hasParentCompound() );

    assert (!hasBondCenter(centerName));

    const Compound::AtomIndex atomId = getAtomInfo(atomName).getIndex();
    const CompoundAtom::BondCenterIndex centerIndex = updAtom(atomId).addSecondBondCenter(bondAngle1);

    addBondCenterInfo(atomId, centerIndex);
    bondCenterIndicesByName[centerName] = getBondCenterInfo(atomId, centerIndex).getIndex();

    assert (hasBondCenter(centerName));

    return *this;
}

CompoundRep& CompoundRep::addPlanarBondCenter(
    const Compound::BondCenterName& centerName, 
    const Compound::AtomName& atomName,
    Angle bondAngle1,
    Angle bondAngle2
    ) 
{
    // Only top level compounds can construct new topology
    // assert(! hasParentCompound() );

    assert (!hasBondCenter(centerName));

    const Compound::AtomIndex atomId = getAtomInfo(atomName).getIndex();
    const CompoundAtom::BondCenterIndex centerIndex = updAtom(atomId).addPlanarBondCenter(bondAngle1, bondAngle2);

    addBondCenterInfo(atomId, centerIndex);
    bondCenterIndicesByName[centerName] = getBondCenterInfo(atomId, centerIndex).getIndex();

    assert (hasBondCenter(centerName));

    return *this;
}

CompoundRep& CompoundRep::addRightHandedBondCenter(
    const Compound::BondCenterName& centerName, 
    const Compound::AtomName& atomName,
    Angle bondAngle1,
    Angle bondAngle2
    ) 
{
    // Only top level compounds can construct new topology
    // assert(! hasParentCompound() );

    assert (!hasBondCenter(centerName));

    const Compound::AtomIndex atomId = getAtomInfo(atomName).getIndex();
    const CompoundAtom::BondCenterIndex centerIndex = updAtom(atomId).addRightHandedBondCenter(bondAngle1, bondAngle2);

    addBondCenterInfo(atomId, centerIndex);
    bondCenterIndicesByName[centerName] = getBondCenterInfo(atomId, centerIndex).getIndex();

    assert (hasBondCenter(centerName));

    return *this;
}

CompoundRep& CompoundRep::addLeftHandedBondCenter(
    const Compound::BondCenterName& centerName, 
    const Compound::AtomName& atomName,
    Angle bondAngle1,
    Angle bondAngle2
    ) 
{
    // Only top level compounds can construct new topology
    // assert(! hasParentCompound() );

    assert (!hasBondCenter(centerName));

    const Compound::AtomIndex atomId = getAtomInfo(atomName).getIndex();
    const CompoundAtom::BondCenterIndex centerIndex = updAtom(atomId).addLeftHandedBondCenter(bondAngle1, bondAngle2);

    addBondCenterInfo(atomId, centerIndex);
    bondCenterIndicesByName[centerName] = getBondCenterInfo(atomId, centerIndex).getIndex();

    assert (hasBondCenter(centerName));

    return *this;
}

CompoundRep& CompoundRep::addBondCenterInfo(
    Compound::AtomIndex atomId,
    CompoundAtom::BondCenterIndex atomCenterIndex)
{
    // Only top level compounds can construct new topology
    // assert(! hasParentCompound() );

    assert( ! hasBondCenter( atomId, atomCenterIndex ) );

    // assert( ! getAtom(atomId).getBondCenter(atomCenterIndex).isBonded() );

    const Compound::BondCenterIndex infoId = Compound::BondCenterIndex(allBondCenters.size());
    allBondCenters.push_back( BondCenterInfo(infoId, atomId, atomCenterIndex) );
    bondCenterIndicesByAtomKey[BondCenterInfo::AtomKey(atomId, atomCenterIndex)] = infoId;
    // bondCenterIndicesByName[centerName] = infoId;

    // Set bonded bit
    BondCenterInfo bondCenterInfo = getBondCenterInfo(infoId);
    const BondCenter& bondCenter = getBondCenter(infoId);
    // bondCenterInfo.setBonded(bondCenter.isBonded());

    assert( hasBondCenter( atomId, atomCenterIndex ) );

    return *this;
}


CompoundRep& CompoundRep::addRingClosingBond(
    const Compound::BondCenterName& centerName1, 
    const Compound::BondCenterName& centerName2,
    mdunits::Length bondLength,
    Angle dihedral,
    BondMobility::Mobility mobility 
    ) 
{
    SimTK_ERRCHK1_ALWAYS(hasBondCenter(centerName1), "Compound::addRingClosingBond()",
        "Couldn't find BondCenter '%s'.\n", centerName1.c_str());
    SimTK_ERRCHK1_ALWAYS(hasBondCenter(centerName2), "Compound::addRingClosingBond()",
        "Couldn't find BondCenter '%s'.\n", centerName2.c_str());

    const Compound::BondCenterIndex id1 = getBondCenterInfo(centerName1).getIndex();
    const Compound::BondCenterIndex id2 = getBondCenterInfo(centerName2).getIndex();

    Compound::BondIndex bondIndex(allBonds.size());
    allBonds.push_back(BondInfo(bondIndex, id1, id2, Bond(bondLength, dihedral, true)));
    indexNewBond(updBondInfo(bondIndex));

    Bond& bond = updBond(updBondInfo(bondIndex));
    bond.setMobility(mobility);

    return *this;
}


CompoundRep& CompoundRep::addRingClosingBond(
    const Compound::BondCenterName& centerName1, 
    const Compound::BondCenterName& centerName2
    ) 
{
    SimTK_ERRCHK1_ALWAYS(hasBondCenter(centerName1), "Compound::addRingClosingBond()",
        "Couldn't find BondCenter '%s'.\n", centerName1.c_str());
    SimTK_ERRCHK1_ALWAYS(hasBondCenter(centerName2), "Compound::addRingClosingBond()",
        "Couldn't find BondCenter '%s'.\n", centerName2.c_str());

    const Compound::BondCenterIndex id1 = getBondCenterInfo(centerName1).getIndex();
    const Compound::BondCenterIndex id2 = getBondCenterInfo(centerName2).getIndex();

    mdunits::Length bondLength = getConsensusBondLength   (getBondCenter(id1), getBondCenter(id2));
    Angle    dihedral   = getConsensusDihedralAngle(getBondCenter(id1), getBondCenter(id2));

    Compound::BondIndex bondIndex(allBonds.size());
    allBonds.push_back(BondInfo(bondIndex, id1, id2, Bond(bondLength, dihedral, true)));
    indexNewBond(updBondInfo(bondIndex));

    return *this;
}

int CompoundRep::getNumAtoms() const {
    return allAtoms.size();
}

const Compound::AtomName CompoundRep::getAtomName(Compound::AtomIndex aid) const {
    const AtomInfo& atomInfo = getAtomInfo(aid);
    // assert (atomInfo.hasValidAtomName()); 
    return atomInfo.getName(); // local name exists
}

BiotypeIndex CompoundRep::getAtomBiotypeIndex(const Compound::AtomIndex aid) const {
    return getAtom(aid).getBiotypeIndex();
}

Compound::AtomIndex CompoundRep::getBondAtomIndex(Compound::BondIndex bid, int which) const {
    assert(which==0 || which==1);
    const BondInfo& binfo = allBonds[bid];
    const Compound::BondCenterIndex parentB = binfo.getParentBondCenterIndex();
    const Compound::BondCenterIndex childB = binfo.getChildBondCenterIndex();
    const BondCenterInfo& pbc = allBondCenters[parentB];
    const BondCenterInfo& cbc = allBondCenters[childB];

    assert(pbc.getAtomIndex() >= 0 && cbc.getAtomIndex() >= 0);

    return which==0 ? pbc.getAtomIndex() : cbc.getAtomIndex();
}

size_t CompoundRep::getNumBondCenters() const {
    return allBondCenters.size();
}

size_t CompoundRep::getNumBondCenters(Compound::AtomIndex atomIndex) const {
    return getAtom(atomIndex).getNumBondCenters();
}

//const Compound::BondCenterName& CompoundRep::getBondCenterName(Compound::BondCenterIndex bondCenterIndex) const {
//    return allBondCenters[bondCenterIndex].getName();
//}

// nameAtom("methyl1/C", "C1", Biotype::EthaneC());
CompoundRep& CompoundRep::nameAtom(const Compound::AtomName& newName, const Compound::AtomPathName& oldName)
{
    const Compound::AtomIndex atomId = getAtomInfo(oldName).getIndex();
    nameAtom(newName, atomId);

    return *this;
}

// nameAtom("methyl1/C", "C1", Biotype::EthaneC());
CompoundRep& CompoundRep::nameAtom(const Compound::AtomName& newName, Compound::AtomIndex atomId)
{
    assert( !hasAtom(newName) );
    assert( CompoundPathName::isValidAtomName(newName) );

    // notice that that all atom path names are placed into index when adding subcompounds
    AtomInfo& atomInfo = updAtomInfo(atomId);

    // NOTE - most recent name is the official name
    atomInfo.addName(newName);

    atomIdsByName[newName] = atomId;

    assert( hasAtom(atomId) );
    assert( hasAtom(newName) );

    return *this;
}

// nameAtom("methyl1/C", "C1", Biotype::EthaneC());
CompoundRep& CompoundRep::nameAtom(const Compound::AtomName& newName, const Compound::AtomPathName& oldName, BiotypeIndex biotypeIx)
{
    nameAtom(newName, oldName);
    setBiotypeIndex(newName, biotypeIx);

    return *this;
}

// setBiotype("C", Biotype::MethaneC);
CompoundRep& CompoundRep::setBiotypeIndex(const Compound::AtomName& atomName, BiotypeIndex biotype) {
    // AtomInfo& atom = atoms[atomIdsByName[atomName]];
    CompoundAtom& atom = updAtom(atomName);
    atom.setBiotypeIndex(biotype);

    return *this;
}
    
CompoundRep& CompoundRep::nameBondCenter(Compound::BondCenterName newName, Compound::BondCenterPathName oldName) {
    assert( hasBondCenter(oldName) );
    // assert( ! hasBondCenter(newName) );

    const BondCenterInfo& bondCenterInfo = getBondCenterInfo(oldName);
    bondCenterIndicesByName[newName] = bondCenterInfo.getIndex();
    // bondCenterInfo.setName(newName);

    assert( hasBondCenter(newName) );

    return *this;
}

Compound::BondCenterIndex CompoundRep::getBondCenterIndex(const Compound::BondCenterName& name) const
{
    assert(hasBondCenter(name));
    return bondCenterIndicesByName.find(name)->second;
}

BondCenterInfo& CompoundRep::updBondCenterInfo(const Compound::BondCenterName& name) 
{
    return updBondCenterInfo(getBondCenterIndex(name));
}

const BondCenterInfo& CompoundRep::getBondCenterInfo(const Compound::BondCenterName& name) const
{
    return getBondCenterInfo(getBondCenterIndex(name));
}

BondCenterInfo& CompoundRep::updBondCenterInfo(Compound::BondCenterIndex id) {
    assert(0 <= id);
    assert((Compound::BondCenterIndex)allBondCenters.size() > id);
    return allBondCenters[id];
}  
const BondCenterInfo& CompoundRep::getBondCenterInfo(Compound::BondCenterIndex id) const {
    assert(0 <= id);
    assert((Compound::BondCenterIndex)allBondCenters.size() > id);
    return allBondCenters[id];
}


// return the bond center on atom1 that is attached to atom2
const BondCenterInfo& CompoundRep::getBondCenterInfo(const Compound::AtomName& atom1, const Compound::AtomName& atom2) const 
{
    return getBondCenterInfo(getAtomInfo(atom1), getAtomInfo(atom2));
}
BondCenterInfo& CompoundRep::updBondCenterInfo(const Compound::AtomName& atom1, const Compound::AtomName& atom2) 
{
    return updBondCenterInfo(getAtomInfo(atom1), getAtomInfo(atom2));
}

const BondCenterInfo& CompoundRep::getBondCenterInfo(const AtomInfo& atom1, const AtomInfo& atom2) const 
{
    const std::pair<Compound::AtomIndex, Compound::AtomIndex> key( atom1.getIndex(), atom2.getIndex() );
    assert( bondIndicesByAtomIndexPair.find(key)!= bondIndicesByAtomIndexPair.end() );
    Compound::BondIndex bondIndex = bondIndicesByAtomIndexPair.find(key)->second;

    const BondInfo& bondInfo = getBondInfo(bondIndex);
    const BondCenterInfo& bc1 = getBondCenterInfo( bondInfo.getParentBondCenterIndex() );
    const BondCenterInfo& bc2 = getBondCenterInfo( bondInfo.getChildBondCenterIndex() );

    assert(bc1.isBonded());
    assert(bc2.isBonded());

    if (bc1.getAtomIndex() == atom1.getIndex()) return bc1;
    else if (bc2.getAtomIndex() == atom1.getIndex()) return bc2;
    else {
        assert(!"neither bond center belonged to atom1");
        return bc1; // to make compiler happy
    }
}

BondCenterInfo& CompoundRep::updBondCenterInfo(const AtomInfo& atom1, const AtomInfo& atom2) 
{
    const std::pair<Compound::AtomIndex, Compound::AtomIndex> key( atom1.getIndex(), atom2.getIndex() );
    assert( bondIndicesByAtomIndexPair.find(key)!= bondIndicesByAtomIndexPair.end() );
    const Compound::BondIndex bondIndex = bondIndicesByAtomIndexPair.find(key)->second;

    const BondInfo& bondInfo = getBondInfo(bondIndex);
    BondCenterInfo& bc1 = updBondCenterInfo( bondInfo.getParentBondCenterIndex() );
    BondCenterInfo& bc2 = updBondCenterInfo( bondInfo.getChildBondCenterIndex() );

    assert(bc1.isBonded());
    assert(bc2.isBonded());

    if (bc1.getAtomIndex() == atom1.getIndex()) return bc1;
    else if (bc2.getAtomIndex() == atom1.getIndex()) return bc2;
    else {
        assert(false);
        return bc1; // to make compiler happy
    }
}

BondCenterInfo& CompoundRep::updBondCenterInfo(Compound::AtomIndex atomId, CompoundAtom::BondCenterIndex atomBondCenterIndex) {
    BondCenterInfo::AtomKey key(atomId, atomBondCenterIndex);
    return updBondCenterInfo(key);
}

const BondCenterInfo& CompoundRep::getBondCenterInfo(Compound::AtomIndex atomId, CompoundAtom::BondCenterIndex atomBondCenterIndex) const {
    BondCenterInfo::AtomKey key(atomId, atomBondCenterIndex);
    return getBondCenterInfo(key);
}

BondCenterInfo& CompoundRep::updBondCenterInfo(BondCenterInfo::AtomKey key) {
    assert(hasBondCenter(key));
    const Compound::BondCenterIndex bcId = bondCenterIndicesByAtomKey.find(key)->second;
    return updBondCenterInfo(bcId);
}
const BondCenterInfo& CompoundRep::getBondCenterInfo(BondCenterInfo::AtomKey key) const {
    assert(hasBondCenter(key));
    const Compound::BondCenterIndex bcId = bondCenterIndicesByAtomKey.find(key)->second;
    return getBondCenterInfo(bcId);
}

bool CompoundRep::hasBondCenter(Compound::BondCenterIndex id) const {
    if ( id < 0 ) return false;
    if ( id >= (Compound::BondCenterIndex)allBondCenters.size() ) return false;
    return true;
}

bool CompoundRep::hasBondCenter(const Compound::BondCenterName& name) const {
    return ( bondCenterIndicesByName.find(name) != bondCenterIndicesByName.end() );
}

// Use atoms names as found in subcompound
// But only "simple" names
CompoundRep& CompoundRep::inheritAtomNames(const Compound::Name& scName) 
{
    std::string prefix(scName + "/");
    std::set<Compound::AtomName> newNames; // for evaluating bond center names
    std::map<Compound::AtomName, Compound::AtomIndex>::const_iterator anI;
    for (anI = atomIdsByName.begin(); anI != atomIdsByName.end(); ++anI) {
        if (anI->first.find(prefix) != 0) continue;
        Compound::AtomName newName = anI->first.substr(prefix.length());
        // But don't inherit sub-sub names
        if (newName.find_first_of("/") != std::string::npos) continue;
        nameAtom(newName, anI->second);
        newNames.insert(newName);
    }

    // TODO - also inherit bond center names with atoms in them
    std::map<String, Compound::BondCenterIndex>::const_iterator scBc;
    for (scBc = bondCenterIndicesByName.begin(); scBc != bondCenterIndicesByName.end(); ++scBc) 
    {
        // Ensure bond center name begins with this subcompound name
        if (scBc->first.find(prefix) != 0) continue;

        const String centerName(scBc->first.substr(prefix.length()));

        // Ensure truncated name begins with an atom name
        if (centerName.find("/") == std::string::npos) continue;
        int slashPos = centerName.find("/");
        String possibleAtomName = centerName.substr(0, slashPos);
        if (newNames.find(possibleAtomName) == newNames.end()) continue;

        nameBondCenter(centerName, scBc->first);
    }
    // assert(false);

    return *this;
}

CompoundRep& CompoundRep::inheritBondCenterNames(const Compound::Name& scName) 
{
    std::string prefix(scName + "/");
    std::map<String, Compound::BondCenterIndex>::const_iterator scBc;
    for (scBc = bondCenterIndicesByName.begin(); scBc != bondCenterIndicesByName.end(); ++scBc) 
    {
        if (scBc->first.find(prefix) != 0) continue;

        const String centerName(scBc->first.substr(prefix.length()));
        // but don't inherit special "inboard bond" name
        if (centerName == InboardBondName) continue;

        nameBondCenter(centerName, scBc->first);
    }

    return *this;
}

// One argument version of writeDefaultPdb begins numbering atoms at 1
std::ostream& CompoundRep::writeDefaultPdb(std::ostream& os, const Transform& transform) const 
{
    PdbChain chain(getOwnerHandle(), transform);

    int nextSerialNumber = 1; // so we can pass a reference
    chain.write(os, nextSerialNumber);

    // writeDefaultPdb(os, nextSerialNumber, transform);

    return os;
}

std::ostream& CompoundRep::writeDefaultPdb(std::ostream& os, int& nextSerialNumber, const Transform& t) const 
{
    PdbChain chain(getOwnerHandle(), t);

    chain.write(os, nextSerialNumber);

    return os;
}


ostream& CompoundRep::writePdb(const State& state, ostream& os, const Transform& transform) const  
{
    // writePdb(state, os, nextSerialNumber, transform);

    PdbChain chain(state, getOwnerHandle(), transform);

    int nextSerialNumber = 1; // so we can pass a reference
    chain.write(os, nextSerialNumber);

    return os;
}

std::ostream& CompoundRep::writePdb(const State& state, std::ostream& os, int& nextSerialNumber, const Transform& transform) const 
{
    PdbChain chain(state, getOwnerHandle(), transform);

    chain.write(os, nextSerialNumber);

    return os;
}

// protected:

bool CompoundRep::hasInboardBondCenter() const {
    return bondCenterIndicesByName.find(InboardBondName) != bondCenterIndicesByName.end();
}

CompoundRep& CompoundRep::convertInboardBondCenterToOutboard() 
{
    // Only top level compounds can construct new topology
    // assert(! hasParentCompound() );

    if (hasInboardBondCenter()) {
        assert(!getInboardBondCenter().isBonded());
        updInboardBondCenter().setInboard(false);
        bondCenterIndicesByName.erase(InboardBondName);
    }
    return *this;
};

const BondCenter& CompoundRep::getInboardBondCenter() const {
    return getBondCenter(InboardBondName);
}
BondCenter& CompoundRep::updInboardBondCenter() {
    return updBondCenter(InboardBondName);
}


const BondCenterInfo& CompoundRep::getInboardBondCenterInfo() const {
    return getBondCenterInfo(InboardBondName);
}
BondCenterInfo& CompoundRep::updInboardBondCenterInfo() {
    return updBondCenterInfo(InboardBondName);
}


CompoundRep& CompoundRep::setInboardBondCenter(const Compound::BondCenterName& n) 
{
    // Only top level compounds can construct new topology
    // assert(! hasParentCompound() );

    // Order operations so that final internal bond name is centerName, rather
    // than InboardBondName
    assert(! hasInboardBondCenter() );
    assert( hasBondCenter(n) );

    const BondCenterInfo& bondCenterInfo = getBondCenterInfo(n);
    BondCenter& bondCenter = updBondCenter(n);

    assert(! bondCenter.isBonded() );

    bondCenter.setInboard(true);

    bondCenterIndicesByName[InboardBondName] = bondCenterInfo.getIndex();

    assert(hasInboardBondCenter());

    return *this;
}

CompoundRep& CompoundRep::setInboardBondCenter(Compound::BondCenterIndex id) 
{
    // Only top level compounds can construct new topology
    // assert(! hasParentCompound() );

    // Order operations so that final internal bond name is centerName, rather
    // than InboardBondName
    assert(! hasInboardBondCenter() );
    assert( hasBondCenter(id) );

    const BondCenterInfo& bondCenterInfo = getBondCenterInfo(id);
    BondCenter& bondCenter = updBondCenter(id);

    assert(! bondCenter.isBonded() );

    bondCenterIndicesByName[InboardBondName] = bondCenterInfo.getIndex();

    assert(hasInboardBondCenter());

    return *this;
}

// Returns new bond center index for inboard center of new subcompound
Compound::BondCenterIndex CompoundRep::addLocalCompound(
    const Compound::Name&   scName, 
    const Compound&         subcompound,
    const Transform& location) 
{
    Compound::BondCenterIndex inboardIndex = absorbSubcompound(scName, subcompound, true);
    if (inboardIndex.isValid()) {
        const BondCenterInfo& bc = getBondCenterInfo(inboardIndex);
        CompoundAtom& atom = updAtom(bc.getAtomIndex());
        const AtomInfo& atomInfo = getAtomInfo(bc.getAtomIndex());
        assert(atomInfo.isBaseAtom());
        atom.setDefaultFrameInCompoundFrame(location * atom.getDefaultFrameInCompoundFrame());
    }
    else {
        CompoundAtom& atom = updAtom(Compound::AtomIndex(0));
        const AtomInfo& atomInfo = getAtomInfo(Compound::AtomIndex(0));
        assert(atomInfo.isBaseAtom());
        atom.setDefaultFrameInCompoundFrame(location * atom.getDefaultFrameInCompoundFrame());
    }
    return inboardIndex;
}


Compound::BondCenterIndex CompoundRep::absorbSubcompound(
    const Compound::Name& scName,
    const Compound& subcompound,
    bool isBaseCompound)
{
    const CompoundRep& subcompoundRep  = subcompound.getImpl();

    // copy atoms
    std::map<Compound::AtomIndex, Compound::AtomIndex> parentAtomIndicesBySubcompoundAtomIndex;
    for (Compound::AtomIndex a(0); a < subcompound.getNumAtoms(); ++a) {
        const CompoundAtom& childAtom = subcompoundRep.getAtom(a);
        const AtomInfo& childAtomInfo = subcompoundRep.getAtomInfo(a);

        const Compound::AtomIndex parentAtomIndex = Compound::AtomIndex(allAtoms.size());
        bool isBase = childAtomInfo.isBaseAtom();
        if (!isBaseCompound) isBase = false;
        const AtomInfo atomInfo(parentAtomIndex, childAtom, isBase);

        allAtoms.push_back(atomInfo);
        parentAtomIndicesBySubcompoundAtomIndex[a] = parentAtomIndex;
    }
    // Copy atom names
    std::map<Compound::AtomName, Compound::AtomIndex>::const_iterator an;
    for (an = subcompoundRep.atomIdsByName.begin(); an != subcompoundRep.atomIdsByName.end(); ++an)
    {
        Compound::AtomIndex parentIx = parentAtomIndicesBySubcompoundAtomIndex[an->second];
        Compound::AtomName parentName = scName + "/" + an->first;
        AtomInfo& parentAtomInfo = updAtomInfo(parentIx);
        parentAtomInfo.addName(parentName);
        atomIdsByName[parentName] = parentIx;
    }
    // Set "main" atom name last, to make it stick
    for (Compound::AtomIndex a(0); a < subcompound.getNumAtoms(); ++a) {
        const AtomInfo& childAtomInfo = subcompoundRep.getAtomInfo(a);
        Compound::AtomIndex parentIx = parentAtomIndicesBySubcompoundAtomIndex[a];
        AtomInfo& parentAtomInfo = updAtomInfo(parentIx);
        Compound::AtomName parentName = scName + "/" + childAtomInfo.getName();
        parentAtomInfo.addName(parentName);
        atomIdsByName[parentName] = parentIx;
    }

    // copy bondCenters
    // copy even bonded centers, for use in dihedral nomenclature
    std::map<Compound::BondCenterIndex, Compound::BondCenterIndex> parentBondCenterIndicesBySubcompoundBondCenterIndex;
    for (Compound::BondCenterIndex bond(0); bond < subcompound.getNumBondCenters(); ++bond) 
    {
        const BondCenter&   scBc     = subcompoundRep.getBondCenter(bond);
        const BondCenterInfo&     scBcInfo = subcompoundRep.getBondCenterInfo(bond);

        const AtomInfo& atomInfo = getAtomInfo(parentAtomIndicesBySubcompoundAtomIndex[scBcInfo.getAtomIndex()]);

        addBondCenterInfo( atomInfo.getIndex(), scBcInfo.getAtomBondCenterIndex() );
        const BondCenterInfo& parentBondCenterInfo = getBondCenterInfo( atomInfo.getIndex(), scBcInfo.getAtomBondCenterIndex() );

        parentBondCenterIndicesBySubcompoundBondCenterIndex[bond] = parentBondCenterInfo.getIndex();
    }
    // copy BondCenter names
    std::map<String, Compound::BondCenterIndex>::const_iterator bcn;
    for (bcn = subcompoundRep.bondCenterIndicesByName.begin(); bcn != subcompoundRep.bondCenterIndicesByName.end(); ++bcn)
    {
        Compound::BondCenterIndex parentIx = parentBondCenterIndicesBySubcompoundBondCenterIndex[bcn->second];
        Compound::BondCenterName parentName = scName + "/" + bcn->first;
        bondCenterIndicesByName[parentName] = parentIx;
    }

    // copy bonds
    std::map<Compound::BondIndex, Compound::BondIndex> parentBondIndicesBySubcompoundBondIndex;
    for (Compound::BondIndex bond(0); bond < subcompoundRep.getNumBonds(); ++bond)
    {
        // 1) Index Info relative to subcompound
        const BondInfo&       scBondInfo = subcompoundRep.getBondInfo(bond);
        const BondCenterInfo& scBc1      = subcompoundRep.getBondCenterInfo( scBondInfo.getParentBondCenterIndex() );
        const BondCenterInfo& scBc2      = subcompoundRep.getBondCenterInfo( scBondInfo.getChildBondCenterIndex() );
        const AtomInfo&       scAtom1    = subcompoundRep.getAtomInfo(scBc1.getAtomIndex());
        const AtomInfo&       scAtom2    = subcompoundRep.getAtomInfo(scBc2.getAtomIndex());

        // 2) Index Info relative to parent compound
        const AtomInfo& atom1Info = getAtomInfo(parentAtomIndicesBySubcompoundAtomIndex[scBc1.getAtomIndex()]);
        const AtomInfo& atom2Info = getAtomInfo(parentAtomIndicesBySubcompoundAtomIndex[scBc2.getAtomIndex()]);
        BondCenterInfo& bc1Info   = updBondCenterInfo(atom1Info.getIndex(), scBc1.getAtomBondCenterIndex());
        BondCenterInfo& bc2Info   = updBondCenterInfo(atom2Info.getIndex(), scBc2.getAtomBondCenterIndex());

        const Compound::BondIndex bondIndex = Compound::BondIndex(allBonds.size());
        allBonds.push_back( BondInfo(
            bondIndex, 
            bc1Info.getIndex(), 
            bc2Info.getIndex(),
            scBondInfo.getBond()
            ) );
        const BondInfo& bondInfo = getBondInfo(bondIndex);

        parentBondIndicesBySubcompoundBondIndex[bond] = bondIndex;

        // Update cross references in BondCenters
        bc1Info.setBondIndex(bondIndex);
        bc2Info.setBondIndex(bondIndex);
        bc1Info.setBondPartnerBondCenterIndex(bc2Info.getIndex());
        bc2Info.setBondPartnerBondCenterIndex(bc1Info.getIndex());

        // map bonds to atomId pairs
        const std::pair<Compound::AtomIndex, Compound::AtomIndex> key1(atom1Info.getIndex(), atom2Info.getIndex());
        const std::pair<Compound::AtomIndex, Compound::AtomIndex> key2(atom2Info.getIndex(), atom1Info.getIndex());
        bondIndicesByAtomIndexPair[key1] = bondIndex;
        bondIndicesByAtomIndexPair[key2] = bondIndex;
    }

    // copy dihedral angles
    std::map<String, DihedralAngle>::const_iterator dn;
    for (dn = subcompoundRep.dihedralAnglesByName.begin(); dn != subcompoundRep.dihedralAnglesByName.end(); ++dn)
    {
        const DihedralAngle& scAngle = dn->second;
        Compound::BondCenterIndex bc1 = parentBondCenterIndicesBySubcompoundBondCenterIndex[scAngle.getBondCenter1Id()];
        Compound::BondCenterIndex bc2 = parentBondCenterIndicesBySubcompoundBondCenterIndex[scAngle.getBondCenter2Id()];

        Compound::AtomName parentName = scName + "/" + dn->first;
        dihedralAnglesByName[parentName] = DihedralAngle(bc1, bc2, scAngle.getNomenclatureOffset());
    }

    if (!subcompoundRep.hasInboardBondCenter())
        return Compound::BondCenterIndex(); // invalid index
    else
        return parentBondCenterIndicesBySubcompoundBondCenterIndex[subcompoundRep.getInboardBondCenterInfo().getIndex()];
}


bool CompoundRep::hasAtom(const Compound::AtomName& atomName) const 
{
    return atomIdsByName.find(atomName) != atomIdsByName.end();
}


Compound::AtomIndex CompoundRep::getAtomIndex(const Compound::AtomName& atomName) const {
    assert(hasAtom(atomName));
    return atomIdsByName.find(atomName)->second;
}

AtomInfo& CompoundRep::updAtomInfo(const Compound::AtomName& atomName) {
    assert(hasAtom(atomName));

    const Compound::AtomIndex compoundAtomIndex = getAtomIndex(atomName);
    return updAtomInfo(compoundAtomIndex);
}

AtomInfo& CompoundRep::updAtomInfo(Compound::AtomIndex i) {
    assert(hasAtom(i));

    return allAtoms[i];
}


BondCenter& CompoundRep::updBondCenter(const Compound::BondCenterName& name) {
    return updBondCenter(getBondCenterInfo(name));
}

const BondCenter& CompoundRep::getBondCenter(const Compound::BondCenterName& name) const {
    return getBondCenter(getBondCenterInfo(name));
}

BondCenter& CompoundRep::updBondCenter(Compound::BondCenterIndex id) {
    return updBondCenter(getBondCenterInfo(id));
}

const BondCenter& CompoundRep::getBondCenter(Compound::BondCenterIndex id) const {
    return getBondCenter(getBondCenterInfo(id));
}

const BondCenter& CompoundRep::getBondCenter(const BondCenterInfo& info) const {
    const CompoundAtom& atom = getAtom(info.getAtomIndex());
    return atom.getBondCenter(info.getAtomBondCenterIndex());
}

BondCenter& CompoundRep::updBondCenter(const BondCenterInfo& info) {
    CompoundAtom& atom = updAtom(info.getAtomIndex());
    return atom.updBondCenter(info.getAtomBondCenterIndex());
}

Transform CompoundRep::calcDefaultBondCenterFrameInAtomFrame(const BondCenterInfo& info) const 
{
    return getAtom(info.getAtomIndex()).calcDefaultBondCenterFrameInAtomFrame(info.getAtomBondCenterIndex());
}

// Cache method used in O(n) all atom Frame computation
const Transform CompoundRep::calcDefaultBondCenterFrameInCompoundFrame(const BondCenterInfo& info, std::vector<Transform>& atomFrameCache) const
{
    Transform X_compound_center;
    const BondCenter& bondCenter = getBondCenter(info);

    if (bondCenter.isInboard()) {
        assert(bondCenter.isBonded());
        const BondInfo& bondInfo = getBondInfo(info.getBondIndex());
        const Bond& bond = getBond(bondInfo);
        Transform X_parentBC_childBC = bond.getDefaultBondCenterFrameInOtherBondCenterFrame();
        assert(info.getIndex() == bondInfo.getChildBondCenterIndex());
        const BondCenterInfo& parentBondCenterInfo = getBondCenterInfo(bondInfo.getParentBondCenterIndex());
        Transform X_compound_parentBC = calcDefaultBondCenterFrameInCompoundFrame(parentBondCenterInfo, atomFrameCache);
        X_compound_center = X_compound_parentBC * X_parentBC_childBC;
    } 
    else {
        Transform X_compound_atom = calcDefaultAtomFrameInCompoundFrame(info.getAtomIndex(), atomFrameCache);
        Transform X_atom_center = getAtom(info.getAtomIndex()).calcDefaultBondCenterFrameInAtomFrame(info.getAtomBondCenterIndex());
        X_compound_center = X_compound_atom * X_atom_center;
    }

    return X_compound_center;
}

Transform CompoundRep::calcDefaultBondCenterFrameInCompoundFrame(const BondCenterInfo& info) const 
{
    Transform X_compound_center;
    const BondCenter& bondCenter = getBondCenter(info);

    if (bondCenter.isInboard()) {
        assert(bondCenter.isBonded());
        const BondInfo& bondInfo = getBondInfo(info.getBondIndex());
        const Bond& bond = getBond(bondInfo);
        Transform X_parentBC_childBC = bond.getDefaultBondCenterFrameInOtherBondCenterFrame();
        assert(info.getIndex() == bondInfo.getChildBondCenterIndex());
        const BondCenterInfo& parentBondCenterInfo = getBondCenterInfo(bondInfo.getParentBondCenterIndex());
        Transform X_compound_parentBC = calcDefaultBondCenterFrameInCompoundFrame(parentBondCenterInfo);
        X_compound_center = X_compound_parentBC * X_parentBC_childBC;
    } 
    else {
        Transform X_compound_atom = calcDefaultAtomFrameInCompoundFrame(info.getAtomIndex());
        Transform X_atom_center = getAtom(info.getAtomIndex()).calcDefaultBondCenterFrameInAtomFrame(info.getAtomBondCenterIndex());
        X_compound_center = X_compound_atom * X_atom_center;
    }

    return X_compound_center;
}

// for O(n) version of all atom Frame computation
const Transform& CompoundRep::calcDefaultAtomFrameInCompoundFrame(Compound::AtomIndex atomId, std::vector<Transform>& atomFrameCache) const
{
    // Is it already cached?
    Transform& candidate = atomFrameCache[atomId];
    if (! isNaN(candidate.p()[0])) 
        return candidate;

    const AtomInfo& atomInfo = getAtomInfo(atomId);
    const CompoundAtom& atom = atomInfo.getAtom();

    Transform parent_X_atom; // atom frame in this compound's frame
    if (atomInfo.isBaseAtom()) {
        parent_X_atom = atom.getDefaultFrameInCompoundFrame();
    }
    else {
        // Delegate frame to inboard bond center
        Transform inboard_X_atom = atom.calcDefaultFrameInInboardCenterFrame();
        const BondCenterInfo& center = getBondCenterInfo(atomId, atom.getInboardBondCenterIndex());
        Transform parent_X_inboard = calcDefaultBondCenterFrameInCompoundFrame(center, atomFrameCache);
        parent_X_atom = parent_X_inboard * inboard_X_atom;
    }

    candidate = parent_X_atom;
    return candidate;
}

Transform CompoundRep::calcDefaultAtomFrameInCompoundFrame(Compound::AtomIndex atomId) const {
    const AtomInfo& atomInfo = getAtomInfo(atomId);
    const CompoundAtom& atom = atomInfo.getAtom();

    Transform parent_X_atom; // atom frame in this compound's frame
    if (atomInfo.isBaseAtom()) {
        parent_X_atom = atom.getDefaultFrameInCompoundFrame();
    }
    else {
        // Delegate frame to inboard bond center
        Transform inboard_X_atom = atom.calcDefaultFrameInInboardCenterFrame();
        const BondCenterInfo& center = getBondCenterInfo(atomId, atom.getInboardBondCenterIndex());
        Transform parent_X_inboard = calcDefaultBondCenterFrameInCompoundFrame(center);
        parent_X_atom = parent_X_inboard * inboard_X_atom;
    }

    return parent_X_atom;
}

// Get all atom locations at once, for greater efficiency
void CompoundRep::calcDefaultAtomFramesInCompoundFrame(std::vector<Transform>& atomFrameCache) const 
{
    std::vector<AtomInfo>::const_iterator aI;
    for (aI = allAtoms.begin(); aI != allAtoms.end(); ++aI) {
        Transform& transform = atomFrameCache[aI->getIndex()];
        if (isNaN(transform.p()[0]))
            transform = calcDefaultAtomFrameInCompoundFrame(aI->getIndex(), atomFrameCache);
    }
    for (int a = 0; a < getNumAtoms(); ++a) {
        assert(! isNaN(atomFrameCache[a].p()[0]));
    }
}

CompoundRep& CompoundRep::setPdbResidueNumber(int n) {
    pdbResidueNumber = n;
    return *this;
}
int CompoundRep::getPdbResidueNumber() const {return pdbResidueNumber;}

CompoundRep& CompoundRep::setPdbResidueName(const String& n) {
    pdbResidueName = n;
    return *this;
}
const String& CompoundRep::getPdbResidueName() const {return pdbResidueName;}

CompoundRep& CompoundRep::setPdbChainId(char c) {
    pdbChainId = c;

    return *this;
}
char CompoundRep::getPdbChainId() const {return pdbChainId;}



//////////////
// COMPOUND //
//////////////

static String indent(int level) {
    String padding;
    for (int i=0; i<level; ++i)
        padding += "  ";
    return padding;
}


Compound& Compound::setCompoundBondMobility(BondMobility::Mobility mobility) 
{
    for (Compound::BondIndex r(0); r <  getNumBonds(); r++)
        setBondMobility(mobility, r);

    return *this;
}

std::ostream& CompoundRep::dumpCompoundRepToStream(std::ostream& o, int level) const {
    using std::endl;

    o << indent(level) << getCompoundName() << " CompoundRep@" << this;
    if (ownerSystem) o << " is owned by System@" << ownerSystem;
    // << " with Compound Index=" << ixWithinOwnerSystem;
    else o << " is not owned by a System";
    o << " {" << endl;


    o << indent(level) << "all atoms:";
    if (allAtoms.size()) {
        o << endl;
        for (Compound::AtomIndex i(0); i < (int)allAtoms.size(); ++i) {
            o << indent(level+1) << allAtoms[i] << endl;
            o << indent(level+1) << "--> " << calcDefaultAtomFrameInCompoundFrame(i).p() << endl;

            const CompoundAtom& atom = getAtom(i);
            o << indent(level+1) << atom << endl;
        }
    }
    else o << " NONE\n";

    o << indent(level) << "all bonds:";
    if (allBonds.size()) {
        o << endl;
        for (Compound::BondIndex i(0); i < (int)allBonds.size(); ++i) {
            o << indent(level+1) << allBonds[i] << endl;
        }
    }
    else o << " NONE\n";

    return o << "}" << endl;
}

template class PIMPLHandle<Compound,CompoundRep>; // instantiate everything

Compound::Compound()
  : HandleBase( new CompoundRep() ) // ensure one-to-one Compound/Rep pairing
{
}

/*explicit*/
Compound::Compound(const String& type)
  : HandleBase( new CompoundRep(type) ) // ensure one-to-one Compound/Rep pairing
{
}

//void Compound::setCompoundSystem(CompoundSystem& system, Compound::Index id) {
//    updImpl().setCompoundSystem(system,id);
//}
void Compound::setMultibodySystem(MultibodySystem& system) {
    updImpl().setMultibodySystem(system);
}

Compound& Compound::setCompoundName(const String& name) {
    updImpl().setCompoundName(name);
    return *this;
}
const String& Compound::getCompoundName() const {
    return getImpl().getCompoundName();
}
Compound& Compound::addCompoundSynonym(const Compound::Name& synonym) {
    updImpl().addCompoundSynonym(synonym);
    return *this;
}

Compound& Compound::setBaseAtom(const Compound::AtomName& name, const Element& element, const Transform& location) 
{
    updImpl().setBaseAtom(name, element, location);
    return *this;
}
Compound& Compound::setBaseAtom(const Compound::SingleAtom& compound, const Transform& location) 
{
    updImpl().setBaseAtom(compound, location);
    return *this;
}
Compound& Compound::setBaseAtom(
    const Compound::AtomName& name, 
    const Biotype& biotype, 
    const Transform& location) 
{
    updImpl().setBaseAtom(name, biotype, location);
    return *this;
}
Compound& Compound::setBaseCompound(
    const Compound::Name& n, 
    const Compound& c,
    const Transform& location) 
{
    updImpl().setBaseCompound(n, c, location);
    return *this;
}
Compound& Compound::bondAtom(
    const Compound::SingleAtom& c, 
    const BondCenterPathName& parentBondName, 
    mdunits::Length distance,
    Angle dihedral,
    BondMobility::Mobility mobility) 
{
    updImpl().bondAtom(c, parentBondName, distance, dihedral, mobility);
    return *this;
}
Compound& Compound::bondAtom(
    const Compound::SingleAtom& c, 
    const BondCenterPathName& parentBondName) 
{
    updImpl().bondAtom(c, parentBondName);
    return *this;
}
Compound& Compound::setInboardBondCenter(const Compound::BondCenterName& centerName)
{
    updImpl().setInboardBondCenter(centerName);
    return *this;
}


//Compound& Compound::addBondCenter(
//        const Compound::BondCenterName& centerName, 
//        const Compound::AtomName& atomName, 
//        Angle zRotation,
//        Angle oldXRotation)
//{
//    updImpl().addBondCenter(centerName, atomName, zRotation, oldXRotation);
//    return *this;
//}


Compound& Compound::addFirstBondCenter(
    const Compound::BondCenterName& centerName, 
    const Compound::AtomName& atomName)
{
    updImpl().addFirstBondCenter(centerName, atomName);
    return *this;
}


Compound& Compound::addSecondBondCenter(
    const Compound::BondCenterName& centerName, 
    const Compound::AtomName& atomName,
    Angle bondAngle1
    )
{
    updImpl().addSecondBondCenter(centerName, atomName, bondAngle1);
    return *this;
}


Compound& Compound::addPlanarBondCenter(
    const Compound::BondCenterName& centerName, 
    const Compound::AtomName& atomName,
    Angle bondAngle1,
    Angle bondAngle2
    )
{
    updImpl().addPlanarBondCenter(centerName, atomName, bondAngle1, bondAngle2);
    return *this;
}


Compound& Compound::addRightHandedBondCenter(
    const Compound::BondCenterName& centerName, 
    const Compound::AtomName& atomName,
    Angle bondAngle1,
    Angle bondAngle2
    )
{
    updImpl().addRightHandedBondCenter(centerName, atomName, bondAngle1, bondAngle2);
    return *this;
}


Compound& Compound::addLeftHandedBondCenter(
    const Compound::BondCenterName& centerName, 
    const Compound::AtomName& atomName,
    Angle bondAngle1,
    Angle bondAngle2
    )
{
    updImpl().addLeftHandedBondCenter(centerName, atomName, bondAngle1, bondAngle2);
    return *this;
}


Compound& Compound::bondCompound(
    const Compound::Name& n, 
    const Compound& c, 
    const BondCenterPathName& parentBondName, 
    mdunits::Length distance,
    Angle dihedral,
    BondMobility::Mobility mobility
    ) 
{
    updImpl().bondCompound(n, c, parentBondName, distance, dihedral, mobility);
    return *this;
}

Compound& Compound::bondCompound(
    const Compound::Name& n, 
    const Compound& c, 
    const BondCenterPathName& parentBondName) 
{
    updImpl().bondCompound(n, c, parentBondName);
    return *this;
}

    //sam added this polymorphism
    /// This is to make it easy to add a residue to a chain and have the bond connecting the residue to the existing chain
    /// have a non-default mobility, e.g. Rigid.
Compound&    Compound::bondCompound(
        const Compound::Name& n,
        const Compound& c,
        const BondCenterPathName& parentBondName,
        BondMobility::Mobility mobility )
{
    updImpl().bondCompound(n, c, parentBondName,mobility);
    return *this;

}

//Compound& removeSubcompound(const Compound::Name& name) {
//    ((CompoundRep*)rep)->removeSubcompound(name);
//    return *this;
//}

Compound& Compound::addRingClosingBond(
    const Compound::BondCenterName& centerName1, 
    const Compound::BondCenterName& centerName2,
    mdunits::Length bondLength,
    Angle dihedral,
    BondMobility::Mobility mobility
    ) 
{
    updImpl().addRingClosingBond(centerName1, centerName2, bondLength, dihedral, mobility);
    return *this;
}

Compound& Compound::addRingClosingBond(
    const Compound::BondCenterName& centerName1, 
    const Compound::BondCenterName& centerName2
    ) 
{
    updImpl().addRingClosingBond(centerName1, centerName2);
    return *this;
}

Compound& Compound::setDefaultInboardBondLength(mdunits::Length distance) {
    updImpl().setDefaultInboardBondLength(distance);
    return *this;
}

Compound& Compound::setDefaultInboardDihedralAngle(Angle angle) {
    updImpl().setDefaultInboardDihedralAngle(angle);
    return *this;
}

Compound& Compound::defineDihedralAngle(
    const Compound::DihedralName& angleName,
    const Compound::AtomName& atom1,
    const Compound::AtomName& atom2,
    const Compound::AtomName& atom3,
    const Compound::AtomName& atom4,
    Angle offset
    ) 
{
    updImpl().defineDihedralAngle(angleName, atom1, atom2, atom3, atom4, offset);
    return *this;
}

Compound& Compound::defineDihedralAngle(
    const Compound::DihedralName& angleName,
    const Compound::BondCenterName& bond1,
    const Compound::BondCenterName& bond2,
    Angle offset
    )
{
    updImpl().defineDihedralAngle(angleName, bond1, bond2, offset);
    return *this;
}


//Compound& Compound::nameBond(const String& newBondName, const String& oldBondName) {
//    ((CompoundRep*)rep)->nameBond(newBondName, oldBondName);
//    return *this;
//}
Compound& Compound::setDefaultDihedralAngle(const String& bondName, Angle angle) {
    updImpl().setDefaultDihedralAngle(bondName, angle);
    return *this;
}

Compound& Compound::setDefaultDihedralAngle(
			Angle angle, 
			Compound::AtomIndex atom1,
			Compound::AtomIndex atom2,
			Compound::AtomIndex atom3,
			Compound::AtomIndex atom4
			) {
    updImpl().setDefaultDihedralAngle(angle, atom1, atom2, atom3, atom4);
    return *this;
}

Compound& Compound::setDefaultDihedralAngle(
			Angle angle, 
			const Compound::AtomName& atom1,
			const Compound::AtomName& atom2,
			const Compound::AtomName& atom3,
			const Compound::AtomName& atom4
			) {
    updImpl().setDefaultDihedralAngle(angle, atom1, atom2, atom3, atom4);
    return *this;
}

Angle Compound::calcDefaultDihedralAngle(const String& bondName) const {
    return getImpl().calcDefaultDihedralAngle(bondName);
}


Compound& Compound::setDihedralAngle(State& state, const String& bondName, Angle angleInRadians)
{
    updImpl().setDihedralAngle(state, bondName, angleInRadians);
    return *this;
}

Angle Compound::calcDihedralAngle(const State& state, const String& bondName) const
{
    return getImpl().calcDihedralAngle(state, bondName);
}

Compound& Compound::setDefaultBondAngle(Angle angle, const AtomName& atom1, const AtomName& atom2, const AtomName& atom3) {
    updImpl().setDefaultBondAngle(angle, atom1, atom2, atom3);
    return *this;
}

Compound& Compound::setDefaultBondLength(mdunits::Length length, const AtomName& atom1, const AtomName& atom2) {
    updImpl().setDefaultBondLength(length, atom1, atom2);
    return *this;
}

int Compound::getNumAtoms() const {
    return getImpl().getNumAtoms();
}

const Compound::AtomName Compound::getAtomName(Compound::AtomIndex aid) const {
    return getImpl().getAtomName(aid);
}

Compound::AtomIndex Compound::getAtomIndex(const Compound::AtomName& name) const {
    return getImpl().getAtomInfo(name).getIndex();
}

const Element& Compound::getAtomElement(Compound::AtomIndex atomIndex) const {
    return getImpl().getAtomElement(atomIndex);
}

const Element& Compound::getAtomElement(const Compound::AtomName& atomName) const {
    return getImpl().getAtomElement(atomName);
}

BiotypeIndex Compound::getAtomBiotypeIndex(Compound::AtomIndex aid) const {
    return getImpl().getAtom(aid).getBiotypeIndex();
}

Transform Compound::calcDefaultAtomFrameInCompoundFrame(Compound::AtomIndex aid) const {
    return getImpl().calcDefaultAtomFrameInCompoundFrame(aid);
}

MobilizedBodyIndex Compound::getAtomMobilizedBodyIndex(Compound::AtomIndex atomId) const {
    return getImpl().getAtomMobilizedBodyIndex(atomId);
}

Vec3 Compound::getAtomLocationInMobilizedBodyFrame(Compound::AtomIndex atomId) const {
    return getImpl().getAtomLocationInMobilizedBodyFrame(atomId);
}

Vec3 Compound::calcAtomLocationInGroundFrame(const State& state, Compound::AtomIndex atomId) const {
    return getImpl().calcAtomLocationInGroundFrame(state, atomId);
}

Vec3 Compound::calcAtomVelocityInGroundFrame(const State& state, Compound::AtomIndex atomId) const {
    return getImpl().calcAtomVelocityInGroundFrame(state, atomId);
}

Vec3 Compound::calcAtomAccelerationInGroundFrame(const State& state, Compound::AtomIndex atomId) const {
    return getImpl().calcAtomAccelerationInGroundFrame(state, atomId);
}

size_t Compound::getNumBonds() const {
    return getImpl().getNumBonds();
}

Compound::AtomIndex Compound::getBondAtomIndex(Compound::BondIndex bid, int which) const {
    return getImpl().getBondAtomIndex(bid,which);
}


size_t Compound::getNumBondCenters() const {
    return getImpl().getNumBondCenters();
}

size_t Compound::getNumBondCenters(Compound::AtomIndex atomIndex) const {
    return getImpl().getNumBondCenters(atomIndex);
}

//const Compound::BondCenterName& Compound::getBondCenterName(Compound::BondCenterIndex bondCenterIndex) const {
//    return getImpl().getBondCenterName(bondCenterIndex);
//}

Compound& Compound::nameAtom(const Compound::AtomName& newName, const AtomPathName& oldName) 
{
    updImpl().nameAtom(newName, oldName);
    return *this;
}
Compound& Compound::nameAtom(const Compound::AtomName& newName, const AtomPathName& oldName, BiotypeIndex biotype) 
{
    updImpl().nameAtom(newName, oldName, biotype);
    return *this;
}
Compound& Compound::setBiotypeIndex(const Compound::AtomName& atomName, BiotypeIndex biotypeIx) 
{
    updImpl().setBiotypeIndex(atomName, biotypeIx);
    return *this;
}
Compound& Compound::nameBondCenter(const Compound::BondCenterName& centerName, const BondCenterPathName& pathName) 
{
    updImpl().nameBondCenter(centerName, pathName);
    return *this;
}

Compound& Compound::inheritAtomNames(const Compound::Name& compoundName) 
{
    updImpl().inheritAtomNames(compoundName);
    return *this;
}

Compound& Compound::inheritBondCenterNames(const Compound::Name& compoundName) 
{
    updImpl().inheritBondCenterNames(compoundName);
    return *this;
}

ostream& Compound::writeDefaultPdb(ostream& os, int& firstSerialNumber, const Transform& transform) const  
{
    return getImpl().writeDefaultPdb(os, firstSerialNumber, transform);
}
ostream& Compound::writeDefaultPdb(ostream& os, const Transform& transform) const  
{
    return getImpl().writeDefaultPdb(os, transform);
}

/**
 * /brief This polymorphism of writeDefaultPdb takes a char* file name rather than an ostream, to save the user a couple of lines of code.
 *
 */
 
void Compound::writeDefaultPdb(const char* outFileName, const Transform& transform) const  
{
    std::ofstream os(outFileName);
    getImpl().writeDefaultPdb(os, transform);
    os.close();
}


ostream& Compound::writePdb(const SimTK::State& state, ostream& os, const Transform& transform) const  
{
    getImpl().writePdb(state, os, transform);
    return os;
}

ostream& Compound::writePdb(const SimTK::State& state, ostream& os, int& nextAtomSerialNumber, const Transform& transform) const  
{
    getImpl().writePdb(state, os, nextAtomSerialNumber, transform);
    return os;
}

Compound& Compound::convertInboardBondCenterToOutboard() {
    updImpl().convertInboardBondCenterToOutboard();
    return *this;
};

Vec3 Compound::calcDefaultAtomLocationInGroundFrame(const AtomName& name) const {
    return getImpl().calcDefaultAtomLocationInGroundFrame(name);
}

// Compute atom location in local compound frame
Vec3 Compound::calcDefaultAtomLocationInCompoundFrame(const AtomName& name) const {
    return getImpl().calcDefaultAtomLocationInCompoundFrame(name);
}

Compound& Compound::setPdbResidueNumber(int resNum) {
    updImpl().setPdbResidueNumber(resNum);
    return *this;
}
int Compound::getPdbResidueNumber() const {
    return getImpl().getPdbResidueNumber();
}

Compound& Compound::setPdbResidueName(const String& resName) {
    updImpl().setPdbResidueName(resName);
    return *this;
}
const String& Compound::getPdbResidueName() const {
    return getImpl().getPdbResidueName();
}

Compound& Compound::setPdbChainId(char chainId) {
    updImpl().setPdbChainId(chainId);
    return *this;
}
char Compound::getPdbChainId() const {
    return getImpl().getPdbChainId();
}

bool Compound::hasBondCenter(const BondCenterName& n) const {
    return getImpl().hasBondCenter(n);
}
bool Compound::hasAtom(const AtomName& n) const {
    return getImpl().hasAtom(n);
}

void Compound::setDuMMAtomIndex(Compound::AtomIndex aid, DuMM::AtomIndex dummId) {
    updImpl().setDuMMAtomIndex(aid, dummId);
}
DuMM::AtomIndex Compound::getDuMMAtomIndex(Compound::AtomIndex aid) const {
    return getImpl().getDuMMAtomIndex(aid);
}

Compound& Compound::setBondMobility(BondMobility::Mobility mobility, const AtomName& atom1, const AtomName& atom2) 
{
    updImpl().setBondMobility(mobility, atom1, atom2);
    return *this;
}
Compound& Compound::setBondMobility(BondMobility::Mobility mobility, const Compound::BondIndex bondIndex) 
{
    updImpl().setBondMobility(mobility, bondIndex);
    return *this;
}
/*
Biopolymer& Biopolymer::setBondMobility(BondMobility::Mobility mobility) 
{
    for (int q = 0; q <  (..getNumResidues()); q++)
    {
    	for (int r = 0; r <  (getResidue(Compound::Index(q)).getNumBonds(); r++)
	{
		updResidue(Compound::Index(q)).setBondMobility(    mobility,Compound::BondIndex(r));
		//updImpl().setBondMobility(mobility, Compound::BondIndex(r));

	} 
    }	

//updImpl().setBondMobility(mobility, bondIndex);
    return *this;
}
*/

Compound::AtomTargetLocations Compound::createAtomTargets(const PdbStructure& targetStructure, const bool guessCoordinates  ) const 
{
    return getImpl().createAtomTargets(targetStructure,guessCoordinates);
}

Compound::AtomTargetLocations Compound::createAtomTargets(const PdbChain& targetChain,const bool guessCoordinates ) const 
{
    return getImpl().createAtomTargets(targetChain);
}

Compound::AtomTargetLocations Compound::createAtomTargets(const PdbResidue& targetResidue,const bool guessCoordinates ) const 
{
    return getImpl().createAtomTargets(targetResidue, guessCoordinates);
}

Compound& Compound::setTopLevelTransform(const Transform& transform)
{
	updImpl().setTopLevelTransform(transform);
	return *this;
}

const Transform& Compound::getTopLevelTransform() const
{
	return getImpl().getTopLevelTransform();
}

Compound& Compound::matchDefaultAtomChirality(const AtomTargetLocations& atomTargets, Angle planarityTolerance, bool flipAll)
{
    updImpl().matchDefaultAtomChirality(atomTargets, planarityTolerance, flipAll);
    return *this;
}

Compound& Compound::matchDefaultBondLengths(const AtomTargetLocations& atomTargets)
{
    updImpl().matchDefaultBondLengths(atomTargets);
    return *this;
}
Compound& Compound::matchDefaultBondAngles(const AtomTargetLocations& atomTargets)
{
    updImpl().matchDefaultBondAngles(atomTargets);
    return *this;
}
Compound& Compound::matchDefaultDihedralAngles(
        const AtomTargetLocations& atomTargets, 
        Compound::PlanarBondMatchingPolicy policy)
{
    updImpl().matchDefaultDihedralAngles(atomTargets, policy);
    return *this;
}
Compound& Compound::matchDefaultTopLevelTransform(const AtomTargetLocations& atomTargets)
{
    updImpl().matchDefaultTopLevelTransform(atomTargets);
    return *this;
}
TransformAndResidual Compound::getTransformAndResidual(const Compound::AtomTargetLocations& atomTargets) const
{
    return getImpl().getTransformAndResidual(atomTargets);
}

//Compound& Compound::matchDefaultConfiguration(const  AtomTargetLocations& atomTargets, Angle planarityTolerance)
//{
//    
//    matchDefaultAtomChirality(atomTargets,planarityTolerance);
//    matchDefaultBondLengths(atomTargets);
//    matchDefaultBondAngles(atomTargets);
//    matchDefaultDihedralAngles(atomTargets);
//    matchDefaultTopLevelTransform(atomTargets);
//
//    return *this;
//}

Compound& Compound::fitDefaultConfiguration(
        const Compound::AtomTargetLocations& atomTargets,
        SimTK::Real targetRms,
        const bool useObservedPointFitter,
        const Real minimizerTolerance
        ) 
{
    Compound compoundCopy = *this;
    // TODO - DuMM should not be required
    CompoundSystem matchingSystem;
    SimbodyMatterSubsystem matchingMatter(matchingSystem);
    DuMMForceFieldSubsystem dumm(matchingSystem);
    dumm.loadAmber99Parameters();
    dumm.setAllGlobalScaleFactors(0);
    GeneralForceSubsystem forces(matchingSystem);
    matchingSystem.adoptCompound(compoundCopy);
    matchingSystem.modelCompounds();
    matchingSystem.realizeTopology();
    State& state = matchingSystem.updDefaultState();
    matchingSystem.realize(state, Stage::Position);
    // cout << "Number of atom matches(2) = " << optimizationAtomTargets.size() << endl;
    std::map<MobilizedBodyIndex, std::vector<Vec3> > stations;
    std::map<MobilizedBodyIndex, std::vector<Vec3> > targetLocations;
    for (Compound::AtomTargetLocations::const_iterator targetIx = atomTargets.begin(); 
         targetIx != atomTargets.end(); 
         ++targetIx)
    {
        Compound::AtomIndex atomId = targetIx->first;
        MobilizedBodyIndex bodyId = compoundCopy.getAtomMobilizedBodyIndex(atomId);
        stations[bodyId].push_back(compoundCopy.getAtomLocationInMobilizedBodyFrame(atomId));
        targetLocations[bodyId].push_back(targetIx->second);            
    }
    
    // Use ObservedPointFitter to optimize geometry
    std::vector<MobilizedBodyIndex> bodyList;
    std::vector<std::vector<Vec3> > stationList;
    std::vector<std::vector<Vec3> > targetList;
    for (std::map<MobilizedBodyIndex, std::vector<Vec3> >::const_iterator iter = stations.begin(); iter != stations.end(); iter++) {
        bodyList.push_back(iter->first);
        stationList.push_back(iter->second);
        targetList.push_back(targetLocations.find(iter->first)->second);
    }


    // ObservedPointFitter takes a while, and occasionally aborts with line search trouble,
    // So lets try a minimization using custom forces
    //bool useObservedPointFitter = true;
    if (useObservedPointFitter) {
        // sherm 100307: Optimizers now use relative tolerance.
        Real tolerance = .001; // 0.1%
        ObservedPointFitter::findBestFit(matchingSystem, state, bodyList, stationList, targetList, tolerance);
    }
    else {
            const MobilizedBody& groundBody = matchingMatter.getGround();
            for (int b = 0; b < (int)bodyList.size(); ++b)
            {
                    const MobilizedBody& atomBody = matchingMatter.getMobilizedBody(bodyList[b]);
                    for (int s = 0; s < (int)stationList[b].size(); ++s)
                    {
                            const Vec3& atomLocation = stationList[b][s];
                            const Vec3& targetLocation = targetList[b][s];
                            Force::TwoPointLinearSpring(forces, atomBody, atomLocation, groundBody, targetLocation, 1000000.0, 0.0);
                    }
            }


            state = matchingSystem.realizeTopology();
            matchingSystem.realize(state, Stage::Position);
        //Real tolerance = .001; // 0.1%
        //ObservedPointFitter::findBestFit(matchingSystem, state, bodyList, stationList, targetList, tolerance);
            matchingSystem.realize(state, Stage::Position);
		LocalEnergyMinimizer::minimizeEnergy(matchingSystem, state,  minimizerTolerance);

	// Stuff optimized coordinates into a string
	std::ostringstream optimizedPdbStringOut;
	matchingSystem.realize(state, Stage::Position);
	compoundCopy.writePdb(state, optimizedPdbStringOut);
	// Create another PdbStructure, and match the dihedral angles to that
	std::istringstream optimizedPdbStringIn(optimizedPdbStringOut.str());
	PdbStructure optimizedStructure(optimizedPdbStringIn);
	//Compound::AtomTargetLocations optimizedAtomTargets = 
	//        createAtomTargets(optimizedStructure,false); // scf set guessCoordinates to false here to make sure it's done exactly as before
    std::ofstream myofstream("match1e.pdb");
    optimizedStructure.write(   myofstream, SimTK::Transform(Vec3(0)));


    }
    // Stuff optimized coordinates into a string
    std::ostringstream optimizedPdbStringOut;
    matchingSystem.realize(state, Stage::Position);
    compoundCopy.writePdb(state, optimizedPdbStringOut);
    // Create another PdbStructure, and match the dihedral angles to that
    std::istringstream optimizedPdbStringIn(optimizedPdbStringOut.str());
    PdbStructure optimizedStructure(optimizedPdbStringIn);
    Compound::AtomTargetLocations optimizedAtomTargets = 
            createAtomTargets(optimizedStructure,false); // scf set guessCoordinates to false here to make sure it's done exactly as before
    bool matchHydrogenAtomLocations = false;
    if (! matchHydrogenAtomLocations)
        {
           map<Compound::AtomIndex, Vec3>::iterator it;
           map<Compound::AtomIndex, Vec3>::iterator next;
           next = optimizedAtomTargets.begin();
           while (next != optimizedAtomTargets.end())
           {
               it = next;
               Compound::AtomIndex m = (*it).first;
               Element myAtomElement = getAtomElement(m);
               next++;
               if  ((myAtomElement.getName()).compare("hydrogen") == 0)
                    {
                    //cout<<__FILE__<<":"<<__LINE__<<" erasing "<<m<<endl;
                    optimizedAtomTargets.erase(it);
                    }
               //cout<<__FILE__<<":"<<__LINE__<<" "<<m<<","<<(getAtomName(m))<<endl;
           }
        }
    matchDefaultBondLengths(optimizedAtomTargets);
    matchDefaultBondAngles(optimizedAtomTargets);
    matchDefaultDihedralAngles(optimizedAtomTargets);
    // Use original atom locations for top level transform
    matchDefaultTopLevelTransform(optimizedAtomTargets); 
    return *this;
}

/// Write current default(initial) Compound configuration into a PdbChain object
const Compound& Compound::populateDefaultPdbChain(
    class PdbChain& pdbChain, 
    int& defaultNextResidueNumber,
    const Transform& transform) const
{
    getImpl().populateDefaultPdbChain(pdbChain, defaultNextResidueNumber, transform);
    return *this;
}

/// Write dynamic Compound configuration into a PdbChain object
const Compound& Compound::populatePdbChain(
    const State& state, 
    class PdbChain& pdbChain, 
    int& defaultNextResidueNumber,
    const Transform& transform) const
{
    getImpl().populatePdbChain(state, pdbChain, defaultNextResidueNumber, transform);
    return *this;
}

Compound& Compound::inheritCompoundSynonyms(const Compound& otherCompound) {
    const std::set<Compound::Name>& n = otherCompound.getImpl().synonyms;
    updImpl().synonyms.insert(n.begin(), n.end());
    return *this;
}

Compound::Compound(CompoundRep* rep) 
  : HandleBase(rep)
{}

std::ostream& operator<<(std::ostream& o, const Compound& c) {
    c.getImpl().dumpCompoundRepToStream(o);
    return o;
}

////////////////
/// Molecule ///
////////////////

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Molecule, CompoundRep, Compound);

Molecule::Molecule() : Compound(new CompoundRep()) {
}

Molecule::Molecule(CompoundRep* rep) : Compound(rep) {
}

/////////////////////
/// BiopolymerRep ///
/////////////////////

int BiopolymerRep::getNumResidues() const {
    return residues.size();
}

const String& BiopolymerRep::getResidueName(int residueIndex) const {
    return residues[residueIndex].getName();
}

const ResidueInfo& BiopolymerRep::getResidue(ResidueInfo::Index residueIndex) const {
    return residues[residueIndex];
}

ResidueInfo& BiopolymerRep::updResidue(ResidueInfo::Index residueIndex) {
    return residues[residueIndex];
}

const ResidueInfo& BiopolymerRep::getResidue(const Compound::Name& residueName) const {
    return residues[residueIdsByName.find(residueName)->second];
}

ResidueInfo& BiopolymerRep::updResidue(const Compound::Name& residueName) {
    return residues[residueIdsByName.find(residueName)->second];
}


///////////////////
/// ResidueInfo ///
///////////////////

ResidueInfo::ResidueInfo(ResidueInfo::Index ix, 
                         const Compound::Name& name, 
                         const BiopolymerResidue& res, 
                         Compound::AtomIndex atomOffset,
                         char insertionCode) 
    : index(ix),
      nameInCompound(name), 
      pdbResidueName(res.getPdbResidueName()), 
      pdbResidueNumber(res.getPdbResidueNumber()), 
      pdbInsertionCode(insertionCode),
      synonyms(res.getImpl().synonyms)
{
    setOneLetterCode(res.getOneLetterCode());
    for (Compound::AtomIndex a(0); a < res.getNumAtoms(); ++a) {
        const SimTK::AtomInfo& compoundAtom = res.getImpl().getAtomInfo(a);
        ResidueInfo::AtomIndex raIx = addAtom(Compound::AtomIndex(a + atomOffset), res.getAtomName(a));
        ResidueInfo::AtomInfo& residueAtom = updAtomInfo(raIx);
        std::set<Compound::AtomName>::const_iterator nameI;
        for (nameI = compoundAtom.getNames().begin(); nameI != compoundAtom.getNames().end(); ++nameI) {
            residueAtom.synonyms.insert(*nameI);
            atomIdsByName[*nameI] = raIx;
        }
    }
}


//////////////////
/// Biopolymer ///
//////////////////

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Biopolymer, BiopolymerRep, Molecule);

Biopolymer::Biopolymer() : Molecule(new BiopolymerRep()) {
}

int Biopolymer::getNumResidues() const {
    return getImpl().getNumResidues();
}

const ResidueInfo& Biopolymer::getResidue(ResidueInfo::Index residueIndex) const {
    return getImpl().getResidue(residueIndex);
}

ResidueInfo& Biopolymer::updResidue(ResidueInfo::Index residueIndex) {
    return updImpl().updResidue(residueIndex);
}

const ResidueInfo& Biopolymer::getResidue(Compound::Name residueName) const {
    return getImpl().getResidue(residueName);
}

ResidueInfo& Biopolymer::updResidue(Compound::Name residueName) {
    return updImpl().updResidue(residueName);
}

const Compound::Name& Biopolymer::getResidueName(ResidueInfo::Index residueIndex) const {
    return getImpl().getResidueName(residueIndex);
}

// This just sets the residue numbering to start at the provided integer and proceed consecutively.
// Added by Samuel Flores
Biopolymer& Biopolymer::renumberPdbResidues(int firstPdbResidueNumber) 
{
    for (ResidueInfo::Index q(0); q < getNumResidues(); q++)
	    updResidue(q).setPdbResidueNumber(firstPdbResidueNumber+q);

    return *this;
}



void Biopolymer::assignBiotypes() {
    for (ResidueInfo::Index r(0); r < getNumResidues(); ++r) {

        // Assign biotypes depending upon position in chain
        Ordinality::Residue residueOrdinality = Ordinality::Any; // default
        if (r == 0) residueOrdinality = Ordinality::Initial;
        else if (r == (getNumResidues() - 1)) residueOrdinality = Ordinality::Final;

        ResidueInfo& residue = updResidue(r);

        assignResidueBiotypes(r, residueOrdinality);
    }
}

ResidueInfo::Index Biopolymer::appendResidue(const String& resName, const BiopolymerResidue& r) 
{
    ResidueInfo::Index resId(getImpl().residues.size());
    updImpl().residues.push_back(ResidueInfo(resId, resName, r, Compound::AtomIndex(getNumAtoms())));
    updImpl().residueIdsByName[resName] = resId;
    ResidueInfo& residue = updResidue(resId);

    // Note that new (empty) residue has already been added to residue count
    if (getNumResidues() <= 1) {
        setBaseCompound(resName, r);
    }
    else {
        const String& previousResidueName = getResidue(ResidueInfo::Index(getNumResidues()-2)).getName();
        // attach residue directly to the biopolymer instead of to the residue
        bondCompound(resName, r, previousResidueName + "/bondNext");
    }

    return resId;
}

//Sam added polymorphism with mobility parameter
ResidueInfo::Index Biopolymer::appendResidue(const String& resName, const BiopolymerResidue& residue, BondMobility::Mobility mobility) {
    ResidueInfo::Index r = appendResidue(resName, residue);
    setResidueBondMobility(r, mobility);
    return r;
}

Compound::AtomTargetLocations Biopolymer::createAtomTargets(const PdbStructure& targetStructure, const bool guessCoordinates ) const 
{
    if (targetStructure.getNumModels() < 1)  // no models? ...
        return Compound::AtomTargetLocations(); // ... no atom targets
    // Use first model only
    const PdbModel& targetModel = targetStructure.getModel(Pdb::ModelIndex(0));
    char chainId = getPdbChainId();
    if (!targetModel.hasChain(chainId)) // no such chain? ...
        return Compound::AtomTargetLocations(); // ... no atom targets
    // Delegate to method that takes a PdbChain& argument
    const PdbChain& targetChain = targetModel.getChain(chainId);
    return createAtomTargets(targetChain,guessCoordinates);
}

Compound::AtomTargetLocations Biopolymer::createAtomTargets(const PdbChain& targetChain, 
    const bool guessCoordinates // optional parameter, defaults to FALSE
    ) const 
{
    Compound::AtomTargetLocations answer;

    // If chain id does not match, we won't find any matches
    char chainId = getPdbChainId();
    if (targetChain.getChainId() != chainId)
        return answer;

    // Compare residues one at a time
    for (ResidueInfo::Index r(0); r < getNumResidues(); ++r) {
        const ResidueInfo& residue = getResidue(r);
        int residueNumber = residue.getPdbResidueNumber();
        char insertionCode = residue.getPdbInsertionCode(); // TODO - make this an attribute of the residue?
        PdbResidueId residueId(residueNumber, insertionCode);
        // scf added guessCoords exception
        if ((!targetChain.hasResidue(residueId)) && (! guessCoordinates)) 
            continue;

        for(ResidueInfo::AtomIndex a(0); a < residue.getNumAtoms(); ++a) {
            String atomName = residue.getAtomName(a);            
            String oldAtomName = atomName; // atomName can potentially be stuffed with odd values in the loop below.  Let's keep a copy with the original name.
            const std::set<Compound::AtomName>& atomNames = residue.getAtomSynonyms(a);
            std::set<Compound::AtomName>::const_iterator nameIx;
            int myCount =0;
            for (nameIx = atomNames.begin(); nameIx != atomNames.end(); ++nameIx)
            {
                //cout <<__FILE__<<":"<<__LINE__<<" count, *nameIx = "<<myCount<<", "<<*nameIx<<endl;
                myCount++;
            }
            for (nameIx = atomNames.begin(); nameIx != atomNames.end(); ++nameIx)
            {
                //cout <<__FILE__<<":"<<__LINE__<<" *nameIx = "<<*nameIx<<endl;
                atomName = *nameIx;
                if ( targetChain.hasAtom(atomName, residueId) )
                    break;
            }
            if ( targetChain.hasAtom(atomName, residueId) ) {
                const PdbAtom& pdbAtom = targetChain.getAtom(atomName, residueId);
                if (pdbAtom.hasLocation()) {
                    // atom index in Biopolymer, as opposed to in residue
                    Compound::AtomIndex atomIndex = residue.getAtomIndex(a);
                    answer[atomIndex] = pdbAtom.getLocation();

                    std::stringstream sstm;
                    sstm << int(r) <<"/"<< oldAtomName;
                    string atomPath = sstm.str();
                    //cout <<__FILE__<<":"<<__LINE__<<"pdb, default coords: "<<pdbAtom.getLocation()<<" "<<calcDefaultAtomLocationInGroundFrame(atomPath)<<endl;
                }
            }
            else if (guessCoordinates)  //&&
                    //(getAtomElement(residue.getAtomIndex(a) ).getName().compare("hydrogen") != 0)) // scf added
            {
                    std::stringstream sstm;
                    sstm << int(r) <<"/"<< oldAtomName;
                    string atomPath = sstm.str();
                    Compound::AtomIndex atomIndex = residue.getAtomIndex(a);
                    answer[atomIndex] = calcDefaultAtomLocationInGroundFrame(atomPath );
            }



        }
    }

    return answer;
}


// Attempt to deduce correct biotypes from global biotypes database
/// @return true if all atoms have a valid biotype assigned, false otherwise
bool Biopolymer::assignResidueBiotypes(ResidueInfo::Index r, Ordinality::Residue ordinality = Ordinality::Any) 
{
    bool answer = true; // start optimistic

    const ResidueInfo& residue = getResidue(r);

    const std::set<Compound::Name>& resNames = residue.getNames();
    for(ResidueInfo::AtomIndex a(0); a < residue.getNumAtoms(); ++a) {
        CompoundAtom& atom = updImpl().updAtom(residue.getAtomIndex(a));
        // Retain biotypes previously assigned
        if (atom.getBiotypeIndex().isValid()) continue;
        const std::set<Compound::AtomName> atomNames = residue.getAtomInfo(a).getNames();
        BiotypeIndex index = getBiotypeIndex(resNames, atomNames, ordinality);
        if (index.isValid())
            atom.setBiotypeIndex(index);
        else
            answer = false;
    }

    return answer;
}


Biopolymer& Biopolymer::setResidueBondMobility(ResidueInfo::Index r, BondMobility::Mobility mobility) 
{
    ResidueInfo& residue = updResidue(r);

    // Create a list of all atoms in the residue
    std::set<Compound::AtomIndex> residueAtoms;
    for (ResidueInfo::AtomIndex a(0); a < residue.getNumAtoms(); ++a)
        residueAtoms.insert(residue.getAtomIndex(a));

    // Set mobility on bonds that are within residue
    // atoms
    std::set<Compound::AtomIndex>::const_iterator aI;
    for (aI = residueAtoms.begin(); aI !=residueAtoms.end(); ++aI) {
        const CompoundAtom& atom = getImpl().getAtom(*aI);
        // atoms->bondCenters
        for (CompoundAtom::BondCenterIndex bc(0); bc < atom.getNumBondCenters(); ++bc) {
            BondCenterInfo::AtomKey key(*aI, bc);
            const BondCenterInfo& bondCenter = getImpl().getBondCenterInfo(key);
            if (! bondCenter.isBonded()) continue; // skip unbonded centers
            const BondCenterInfo& otherBondCenter = 
                getImpl().getBondCenterInfo(bondCenter.getBondPartnerBondCenterIndex());
            Compound::AtomIndex otherAtomIndex = otherBondCenter.getAtomIndex();
            // skip bonds that go outside of this residue
            if (residueAtoms.find(otherAtomIndex) == residueAtoms.end()) continue;
            BondInfo& bondInfo = updImpl().updBondInfo(bondCenter.getBondIndex());
            updImpl().updBond(bondInfo).setMobility(mobility);
        }
    }

    return *this;
}


////////////////////////////
/// BiopolymerResidueRep /// 
////////////////////////////
 


////////////////////////
/// BipolymerResidue ///
////////////////////////

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(BiopolymerResidue,BiopolymerResidueRep,Compound);

// Attempt to deduce correct biotypes from global biotypes database
bool BiopolymerResidue::assignBiotypes(Ordinality::Residue ordinality) {
    return updImpl().assignBiotypes(ordinality);
}

BiopolymerResidue::BiopolymerResidue(const String& name, const String& threeLetterCode, char oneLetterCode)
    : Compound( new BiopolymerResidueRep(name, threeLetterCode, oneLetterCode) )
{
    // set pdb residue name
    String pdbRes = threeLetterCode;
    // convert to upper case
    std::transform(pdbRes.begin(), pdbRes.end(), pdbRes.begin(), (int(*)(int)) toupper);
    setPdbResidueName(pdbRes);

    setCompoundName(name);
    if ( (threeLetterCode != "UNK") && (threeLetterCode.length() > 1) )
        addCompoundSynonym(threeLetterCode);
}

BiopolymerResidue& BiopolymerResidue::setOneLetterCode(char olc) {
    updImpl().setOneLetterCode(olc);
    return *this;
}
BiopolymerResidue& BiopolymerResidue::setThreeLetterCode(const String& tlc) {
    updImpl().setThreeLetterCode(tlc);
    return *this;
}
BiopolymerResidue& BiopolymerResidue::setResidueTypeName(const String& name) {
    updImpl().setResidueTypeName(name);
    return *this;
}

char BiopolymerResidue::getOneLetterCode() const {
    return getImpl().getOneLetterCode();
}
const String& BiopolymerResidue::getThreeLetterCode() const {
    return getImpl().getThreeLetterCode();
}
const String& BiopolymerResidue::getResidueTypeName() const {
    return getImpl().getResidueTypeName();
}

////////////////////////
/// AminoAcidResidue ///
////////////////////////

// AminoAcidResidue has three unsatisfied BondCenters:
//  1) at the amino nitrogen, for the preceding amino acid residue
//  2) at the carbonyl carbon, for the next amino acid residue
//  3) at the alpha carbon, for the side chain

AminoAcidResidue AminoAcidResidue::create(const PdbResidue& pdbResidue) 
{
    const String& residueName(pdbResidue.getName());
    
    AminoAcidResidue answer = Alanine();
    
    if      (residueName == "ALA") answer = Alanine();
    else if (residueName == "CYS") answer = Cysteine();
    else if (residueName == "ASP") answer = Aspartate();
    else if (residueName == "GLU") answer = Glutamate();
    else if (residueName == "PHE") answer = Phenylalanine();
    else if (residueName == "GLY") answer = Glycine();
    else if (residueName == "HIS") answer = Histidine();
    else if (residueName == "ILE") answer = Isoleucine();
    else if (residueName == "LYS") answer = Lysine();
    else if (residueName == "LEU") answer = Leucine();
    else if (residueName == "MET") answer = Methionine();
    else if (residueName == "ASN") answer = Asparagine();
    else if (residueName == "PRO") answer = Proline();
    else if (residueName == "GLN") answer = Glutamine();
    else if (residueName == "ARG") answer = Arginine();
    else if (residueName == "SER") answer = Serine();
    else if (residueName == "THR") answer = Threonine();
    else if (residueName == "VAL") answer = Valine();
    else if (residueName == "TRP") answer = Tryptophan();
    else if (residueName == "TYR") answer = Tyrosine();
    
    else 
        SimTK_THROW1(Exception::UndefinedAminoAcidResidue, residueName);
    
    answer.setPdbResidueNumber( pdbResidue.getPdbResidueNumber() );
    
    return answer;
}

AminoAcidResidue AminoAcidResidue::create(char oneLetterCode) 
{
    switch(oneLetterCode) {
    case 'A':
        return Alanine();
    case 'C':
        return Cysteine();
    case 'D':
        return Aspartate();
    case 'E':
        return Glutamate();
    case 'F':
        return Phenylalanine();
    case 'G':
        return Glycine();
    case 'H':
        return Histidine();
    case 'I':
        return Isoleucine();
    case 'K':
        return Lysine();
    case 'L':
        return Leucine();
    case 'M':
        return Methionine();
    case 'N':
        return Asparagine();
    case 'P':
        return Proline();
    case 'Q':
        return Glutamine();
    case 'R':
        return Arginine();
    case 'S':
        return Serine();
    case 'T':
        return Threonine();
    case 'V':
        return Valine();
    case 'W':
        return Tryptophan();
    case 'Y':
        return Tyrosine();
    }

    assert(false);

    return Alanine();
}

std::ostream& operator<<(std::ostream& o, const BondInfo& binfo) {
    o << "BONDINFO: BondIndex=" << binfo.id << "inb,outb BondCenterIndices=" 
      << binfo.parentBondCenterIndex << ", " << binfo.childBondCenterIndex ;
    //   << " rotatable=" << binfo.isRotatable << " ringclosing=" << binfo.amRingClosingBond;
    return o;
}


std::ostream& operator<<(std::ostream& o, const CompoundAtom& a) {
    o << "ATOM " << a.element
      << " biotype=" << a.getBiotypeIndex() << " NbondCenters=" << a.bondCenters.size()
      << " localPos=" << a.localTransform.p() << std::endl;
    o << "REMARK-SIMTK-COORDS"<<std::endl;
    return o;
}


std::ostream& operator<<(std::ostream& o, const AtomInfo& ai) {
    o << "ATOMINFO '" << ai.name << "': id=" << ai.id;
      // << " parentCompoundAtomIndex=" << ai.parentCompoundAtomIndex;
    return o << std::endl;
}


/////////////////////////////
/// RibonucleotideResidue ///
/////////////////////////////

RibonucleotideResidue RibonucleotideResidue::create(char oneLetterCode) {
    switch(oneLetterCode) {
    case 'A':
        return Adenylate();
        break;
    case 'C':
        return Cytidylate();
        break;
    case 'G':
        return Guanylate();
        break;
    case 'U':
        return Uridylate();
        break;
    default:
        assert(false);
    }

    assert(false);

    return Adenylate();
}

RibonucleotideResidue RibonucleotideResidue::create(const PdbResidue& pdbResidue, char chainId) 
{
	RibonucleotideResidue answer("Unknown residue");

	String pdbResidueName(pdbResidue.getName());
    // The residue name should be "A  ", but we'll be relaxed about where the
    // spaces are.
    pdbResidueName.trimWhiteSpace();
	if      (pdbResidueName == "A")
		answer = Adenylate();
	else if (pdbResidueName == "C")
		answer = Cytidylate();
	else if (pdbResidueName == "G")
		answer = Guanylate();
	else if (pdbResidueName == "U")
		answer = Uridylate();
	else 
	{
		throw std::range_error(String("Unrecognized Ribonucleotide residue name: '") + pdbResidueName + "'");
	}

	answer.setPdbResidueNumber(pdbResidue.getPdbResidueNumber());
        answer.setPdbChainId(chainId);
	answer.matchDefaultDihedralAngles(answer.createAtomTargets(pdbResidue,false));

	return answer;
}

//////////////////////////////////
/// DeoxyribonucleotideResidue ///
//////////////////////////////////

DeoxyribonucleotideResidue DeoxyribonucleotideResidue::create(char oneLetterCode) {
    switch(oneLetterCode) {
    case 'A':
        return Deoxyadenosine();
        break;
    case 'C':
        return Deoxycytidine();
        break;
    case 'G':
        return Deoxyguanosine();
        break;
    case 'T':
        return Deoxythymidine();
        break;
    default:
        assert(false);
    }

    assert(false);

    return Deoxyadenosine();
}

DeoxyribonucleotideResidue DeoxyribonucleotideResidue::create(const PdbResidue& pdbResidue, char chainId) 
{
	DeoxyribonucleotideResidue answer("Unknown residue");

	String pdbResidueName(pdbResidue.getName());
    // The residue name should be "A  ", but we'll be relaxed about where the
    // spaces are.
    pdbResidueName.trimWhiteSpace();
	if      (pdbResidueName == "A")
		answer = Deoxyadenosine();
	else if (pdbResidueName == "C")
		answer = Deoxycytidine();
	else if (pdbResidueName == "G")
		answer = Deoxyguanosine();
	else if (pdbResidueName == "T")
		answer = Deoxythymidine();
	else 
	{
		throw std::range_error(String("Unrecognized Deoxyribonucleotide residue name: '") + pdbResidueName + "'");
	}

	answer.setPdbResidueNumber(pdbResidue.getPdbResidueNumber());
        answer.setPdbChainId(chainId);
	answer.matchDefaultDihedralAngles(answer.createAtomTargets(pdbResidue,false));

	return answer;
}

void Protein::loadFromPdbChain(const PdbChain& pdbChain, SimTK::Real targetRms)
{
	setPdbChainId( pdbChain.getChainId() );

    // 1) Link canonical residue types into a chain
    // (the conformation will not yet match the pdbChain though)
    for (Pdb::ResidueIndex resIx(0); resIx < pdbChain.getNumResidues();  ++resIx)
    {
        const PdbResidue& pdbResidue = pdbChain.getResidue(resIx);
        
        AminoAcidResidue newResidue = AminoAcidResidue::create(pdbResidue);
        newResidue.assignBiotypes();
        newResidue.setPdbChainId( pdbChain.getChainId() );
        
        const String& residueName = String(pdbResidue.getName()) + String(pdbResidue.getPdbResidueNumber());
        appendResidue( residueName, newResidue );
    }
    
    // Set the dihedrals to those observed, especially on rigid bonds, so ObservedPointFitter will have
    // an easier time.
    // matchDefaultAtomChirality() will prevent optimization trouble for inverted chiral centers
    // matchDefaultDihedralAngles will prevent trouble for cis-peptide bonds
    Compound::AtomTargetLocations atomTargets = createAtomTargets(pdbChain,false);

    matchDefaultTopLevelTransform(atomTargets);
    //cout << "Initial residual = " << getTransformAndResidual(atomTargets).residual << endl;
    matchDefaultAtomChirality(atomTargets); // Chirality must agree for a good match
    matchDefaultDihedralAngles(atomTargets); // e.g. rigid cis-peptide bonds must agree for good match
	matchDefaultTopLevelTransform(atomTargets);
    // Make a copy of the current structure, so we can attach it to a temporary CompoundSystem for optimization        
    Protein proteinCopy = *this;
    proteinCopy.assignBiotypes();
    CompoundSystem matchingSystem;
    DuMMForceFieldSubsystem dumm(matchingSystem);
    SimbodyMatterSubsystem matchingMatter(matchingSystem);
    dumm.loadAmber99Parameters();
    dumm.setAllGlobalScaleFactors(0);
    GeneralForceSubsystem forces(matchingSystem);
    matchingSystem.adoptCompound(proteinCopy);
    matchingSystem.modelCompounds();
    matchingSystem.realizeTopology();
    State& state = matchingSystem.updDefaultState();
    matchingSystem.realize(state, Stage::Position);
    Compound::AtomTargetLocations optimizationAtomTargets = proteinCopy.createAtomTargets(pdbChain,false);
    std::map<MobilizedBodyIndex, std::vector<Vec3> > stations;
    std::map<MobilizedBodyIndex, std::vector<Vec3> > targetLocations;
    for (Compound::AtomTargetLocations::const_iterator targetIx = optimizationAtomTargets.begin(); targetIx != optimizationAtomTargets.end(); ++targetIx)
    {
        Compound::AtomIndex atomId = targetIx->first;
        MobilizedBodyIndex bodyId = proteinCopy.getAtomMobilizedBodyIndex(atomId);
        stations[bodyId].push_back(proteinCopy.getAtomLocationInMobilizedBodyFrame(atomId));
        targetLocations[bodyId].push_back(targetIx->second);            
    }
    
    // Use ObservedPointFitter to optimize geometry
    std::vector<MobilizedBodyIndex> bodyList;
    std::vector<std::vector<Vec3> > stationList;
    std::vector<std::vector<Vec3> > targetList;
    for (std::map<MobilizedBodyIndex, std::vector<Vec3> >::const_iterator iter = stations.begin(); iter != stations.end(); iter++) {
        bodyList.push_back(iter->first);
        stationList.push_back(iter->second);
        targetList.push_back(targetLocations.find(iter->first)->second);
    }


	// ObservedPointFitter takes a while, and occasionally aborts with line search trouble,
	// So lets try a minimization using custom forces
	bool useObservedPointFitter = true;
	if (useObservedPointFitter) {
        // sherm 100307: Optimizers now use relative tolerance.
        Real tolerance = .001; // 0.1%
	Real rmsd = ObservedPointFitter::findBestFit(matchingSystem, state, bodyList, stationList, targetList, tolerance);
        //std::cout << "Rmsd is " << rmsd << std::endl;
	}
	else {
		const MobilizedBody& groundBody = matchingMatter.getGround();
		for (int b = 0; b < (int)bodyList.size(); ++b)
		{
			const MobilizedBody& atomBody = matchingMatter.getMobilizedBody(bodyList[b]);
			for (int s = 0; s < (int)stationList[b].size(); ++s)
			{
				const Vec3& atomLocation = stationList[b][s];
				const Vec3& targetLocation = targetList[b][s];
				Force::TwoPointLinearSpring(forces, atomBody, atomLocation, groundBody, targetLocation, 1000.0, 0.0);
			}
		}
		state = matchingSystem.realizeTopology();
		matchingSystem.realize(state, Stage::Position);
		LocalEnergyMinimizer::minimizeEnergy(matchingSystem, state, 15.0);
	}

  /*  
    // Stuff optimized coordinates into a string
    std::ostringstream optimizedPdbStringOut;
    matchingSystem.realize(state, Stage::Position);
    proteinCopy.writePdb(state, optimizedPdbStringOut);
    
	// std::ofstream match1dOut("match1d.pdb"); // 
	// proteinCopy.writePdb(state, match1dOut);

    // Create another PdbStructure, and match the dihedral angles to that
    std::istringstream optimizedPdbStringIn(optimizedPdbStringOut.str());
    PdbStructure optimizedStructure(optimizedPdbStringIn);
    Compound::AtomTargetLocations optimizedAtomTargets = createAtomTargets(optimizedStructure,false);

	// optimizedStructure.write(std::ofstream("match1e.pdb"));

	// Here is the trouble...
	matchDefaultDihedralAngles(optimizedAtomTargets);

	// Use original atom locations for top level transform
    matchDefaultTopLevelTransform(optimizedAtomTargets);        
*/
	// cout << "Residual after optimization = " << getTransformAndResidual(atomTargets).residual << endl;

	// std::ofstream match2Out("match2.pdb");
	// writeDefaultPdb(match2Out);
}

void invalidateAtomFrameCache(std::vector<Transform>& atomFrameCache, int numAtoms) {
    for (int a = 0; a < numAtoms; ++a)
        atomFrameCache[a].updP() = Vec3(NaN);
}

} // namespace SimTK
