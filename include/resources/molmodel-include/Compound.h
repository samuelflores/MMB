#ifndef SimTK_MOLMODEL_COMPOUND_H_
#define SimTK_MOLMODEL_COMPOUND_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Molmodel                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007-2008 Stanford University and the Authors.      *
 * Authors: Michael Sherman, Christopher Bruns                                *
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
#include "DuMMForceFieldSubsystem.h"
#include "molmodel/internal/Pdb.h"
#include "molmodel/internal/Superpose.h"
#include "molmodel/internal/units.h"
#include <map>

#include <iosfwd> // declare ostream without all the definitions

namespace SimTK {
class CompoundSystem;

// Because of private data, put Compound implementation behind the curtain
class Compound;
class CompoundRep;

// Make sure class is instantiated just once in the library
// #ifndef DO_INSTANTIATE_COMPOUND_PIMPL_HANDLE
//    template class SimTK_MOLMODEL_EXPORT PIMPLHandle<Compound, CompoundRep>;
// #endif

/**
 * \brief Namespace for description of allowed bond motions.
 */
namespace BondMobility {
    /**
     * \brief Which motions are allowed for a particular covalent bond.
     *
     * An enum applied to each covalent bond in a Compound.
     * Used to specify the degrees of freedom of that
     * bond during subsequent construction of a multibody model in simbody.
     */
    enum Mobility {
        Free = 1, ///< Unrestricted bond, permitting changes in stretch, bend, and torsion modes
        Torsion = 2, ///< Bond has fixed length and angles, but permits rotation about the bond axis
        Rigid = 3 ///< Bond links both atoms to the same rigid unit
    };
    static Mobility Default = Torsion;
}

/**
 * \brief The base class for atoms, molecules, and chemical groups.
 *
 * The Compound class is the base for all molecular entities in the SimTK
 * Molmodel API.
 */
class SimTK_MOLMODEL_EXPORT Compound : public PIMPLHandle<Compound,CompoundRep> {
public:

    enum MatchStratagem {
        Match_Exact, //< Try to match input atom positions precisely
        Match_Idealized, //< Try to match atom positions with idealized geometry
        Match_TopologyOnly //< Don't match atom positions, only topology
    };
    
    /// \name Data types that identify subcomponents of molecular compounds
    /// @{

    /**
     * \brief Type for name of a particular subcompound within a Compound.
     *
     * Compound::Name is used as an index to subcompounds of a parent compound, e.g. "methyl group gamma"
     */
    typedef String Name;

    /**
     * \brief Type for name of a particular atom within a Compound.
     *
     * Compound::AtomName is not intrinsic to the atom itself, but rather the relationship
     * between an atom and a particular Compound.  A particular atom may have different
     * names in a compound and its subcompounds.  For example, a particular atom might
     * have the name "methyl H1" in a Compound, and have the name "H1" within a
     * subcompound of that Compound.
     * Unlike AtomPathNames, AtomNames must not be qualified by subcompound indirection.
     * e.g. "H1" is a valid AtomName
     * and a valid AtomPathName, while "methyl/H1" (because of the "/" character) is not a valid AtomName,
     * but is a valid AtomPathName.
     */
    typedef String AtomName;

    /**
     * \brief Type for name of a particular bond center within a Compound or atom.
     *
     * Compound::BondCenterName is not intrinsic to the bond center itself, but rather the 
     * relationship between a bond center and a particular atom or Compound.
     * BondCenterNames are local to a particular compound.  In contrast, BondCenterPathNames
     * may be qualified by subcompound indirection.  e.g. "bond1" is a valid BondCenterName
     * and a valid BondCenterPathName, while "methyl gamma/bond1" is not a valid BondCenterName,
     * but is a valid BondCenterPathName.
     */
    typedef String BondCenterName;

    /**
     * \brief Type for name of a particular named dihedral angle within a Compound.
     *
     * Compound::DihedralName is not intrinsic to the dihedral angle itself, but rather the 
     * relationship between a dihedral angle and a particular Compound.
     */
    typedef String DihedralName;

    /**
     * \brief Type for name of a particular bond center within a Compound, possibly including subcompound indirection.
     *
     * Compound::BondCenterPathName is not intrinsic to the bond center itself, but rather the 
     * relationship between a bond center and a particular Compound.
     * BondCenterPathNames may be qualified by subcompound indirection, in contrast to BondCenterNames.
     * e.g. "bond1" is a valid BondCenterName
     * and a valid BondCenterPathName, while "methyl gamma/bond1" is not a valid BondCenterName,
     * but is a valid BondCenterPathName.
     */
    typedef String BondCenterPathName;

    /**
     * \brief Type for name of a particular atom within a Compound, possibly including subcompound path.
     *
     * Compound::AtomPathName is not intrinsic to the atom itself, but rather the relationship
     * between an atom and a particular Compound.  A particular atom may have different
     * names in a compound and its subcompounds.  For example, a particular atom might
     * have the name "methyl H1" in a Compound, and have the name "H1" within a
     * subcompound of that Compound.
     * AtomPathNames may be qualified by subcompound indirection, in contrast to AtomNames.
     * e.g. "H1" is a valid AtomName
     * and a valid AtomPathName, while "methyl/H1" (because of the "/" character) is not a valid AtomName,
     * but is a valid AtomPathName.
     */
    typedef String AtomPathName;

    // Unique classes which behave like type-safe ints. These are local
    // to Compound.

    /**
     * Compound::Index type is an integer index into subcompounds of a Compound.  It is NOT
     * instrinsic to the subcompound, but represents the relationship between a subcompound
     * and precisely one of its parent compounds.
     */
    // SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(Compound,Index);


    /**
     * Compound::AtomIndex type is an integer index into atoms of a Compound.  It is NOT
     * instrinsic to the atom, but represents the relationship between an atom
     * and precisely one of its parent compounds.
     */
    SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(Compound,AtomIndex);

    /**
     * Compound::LocalAtomIndex type is an integer index into atoms directly attached to a Compound,
     * that is to say that the atom does not belong to any subcompounds of the Compound.  This index 
     * type should ordinarily not be used by clients of the Molmodel API
     */
    SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(Compound,LocalAtomIndex);

    /**
     * Compound::BondCenterIndex type is an integer index into BondCenters of a Compound.  It is NOT
     * instrinsic to the BondCenter, but represents the relationship between a BondCenter
     * and precisely one of its parent compounds.
     */
    SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(Compound,BondCenterIndex);

    /**
     * Compound::BondIndex type is an integer index into Bonds of a Compound.  It is NOT
     * instrinsic to the Bond, but represents the relationship between a Bond
     * and precisely one of its parent compounds.
     */
    SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(Compound,BondIndex);

    /// Type for set of target atom locations to be used for structure matching
    typedef std::map<AtomIndex, Vec3> AtomTargetLocations;

    class SingleAtom;

    /// @}


    /// \name Molecular topology methods
    /// @{

    /**
     * \brief default Compound constructor.
     *
     * Create an empty compound object representing a simulatable molecular structure.
     */
    Compound();

    /**
     * \brief Construct a Compound with a type name.
     *
     * Create an empty compound object representing a simulatable molecular structure.
     * The name is intended to represent the class of molecule being represented.
     * e.g. "ethane", not "ethane number 3"
     */
    explicit Compound(const Name& name ///< name of the compound type, not the compound instance
        );

    /**
     * \return the number of atoms in the Compound, including those within subcompounds.
     */
    int getNumAtoms() const;

    /// \return total number of BondCenters in a Compound, including its subcompounds
    size_t getNumBondCenters() const;

    /// \return total number of BondCenters in a  particular atom of a Compound
	size_t getNumBondCenters(Compound::AtomIndex atomIndex) const;

    /// \return total number of Bonds in a Compound, including its subcompounds
    size_t getNumBonds() const;

    /// \return total number of subcompounds in a compoud, including sub-subcompounds etc.
    // int getNumSubcompounds() const;

    /// \return true if Compound contains an atom by that name (relative to this Compound)
    bool hasAtom(const AtomPathName& name ///< name of the atom relative to this compound.
        ) const;

    /// \return true if Compound contains a BondCenter (half-bond) by that name (relative to this Compound)
    bool hasBondCenter(const BondCenterPathName& bondCenter ///< name of the BondCenter relative to this compound.
        ) const;

    /// \return true if Compound contains a subcompound by that name (relative to this Compound)
    // bool hasSubcompound(const Name& name ///< name of the subcompound relative to this compound.
    //     ) const;

    /// \return integer index of one of the two atoms involved in a particular covalent bond
    Compound::AtomIndex getBondAtomIndex(
        Compound::BondIndex bondIndex, ///< integer index of a bond in this Compound
        int which                      ///< zero(0) for parent-side (rootward) atom of bond, one(1) for child-side (leafward) atom
        ) const;


    /** 
     *  Add the first atom unconnected to anything else (yet).
     *
     * \return a reference to this compound object
     */
    Compound& setBaseAtom(
        const Compound::AtomName& name,   ///< name of the new atom being created
        const Element& element,           ///< chemical element of the new atom being created
        const Transform& location = Vec3(0)  ///< default location of the new atom being created in orthogonal nanometers.  Defaults to (0,0,0)
        );
    /** 
     *  Add the first atom unconnected to anything else (yet).
     *
     * \return a reference to this compound object
     */
    Compound& setBaseAtom(
        const Compound::AtomName& name,   ///< name of the new atom being created
        const Biotype& biotype,           ///< Biotype of the new atom being created
        const Transform& location = Vec3(0)    ///< default location of the new atom being created in orthogonal nanometers.  Defaults to (0,0,0)
        );  

    /** 
     *  Add a first subcompound containing exactly one atom, so the Compound::AtomName can be reused for the Compound::Name.
     *  This atom is not connected to anything else (yet).
     * 
     * \return a reference to this compound object
     */
    Compound& setBaseAtom(
        const Compound::SingleAtom& c, ///< single-atom Compound to add as first subcompound of this Compound
        const Transform& = Transform() ///< default location and orientation of the new atom being created
        ); 

    /** 
     *  Add a first subcompound without attaching it to anything.
     * 
     * \return a reference to this compound object
     */
    Compound& setBaseCompound(
        const Compound::Name& n,         ///< name for new subcompound, from viewpoint of parent 
        const Compound& c,               ///< new subcompound to be copied and attached to parent Compound
        const Transform& location = Transform()///< default location of the new atom being created in orthogonal nanometers.  Defaults to (0,0,0)
        );

    /** 
     *  Add a subcompound containing exactly one atom, so the Compound::AtomName can be reused for the Compound::Name.
     *  This atom is connected to existing material.
     * 
     * \return a reference to this compound object
     */
    Compound& bondAtom(
        const Compound::SingleAtom& c, ///< the new subcompound to attach to this compound (a copy of the subcompound will be attached) 
        const BondCenterPathName& parentBondName, ///< name of the bond center on the parent Compound to which the new subcompound will be attached
        mdunits::Length distance, ///< default bond length in nanometers of new bond connecting new subcompound to parent Compound
        Angle dihedral = 0, ///< default dihedral angle about new bond, with respect to the first bond center on each atom
        BondMobility::Mobility = BondMobility::Default ///< allowed motion in new bond
        );

    /** 
     *  Bond atom using default bond length and dihedral angle.
     *
     *  Bond length and dihedral angle must have already been predefined in the 
     *  parent *or* the child BondCenter, but not both
     * 
     * \return a reference to this compound object
     */
    Compound& bondAtom(
        const Compound::SingleAtom& c, ///< the new subcompound to attach to this compound (a copy of the subcompound will be attached) 
        const BondCenterPathName& parentBondName ///< name of the bond center on the parent Compound to which the new subcompound will be attached
        );

    /** 
     *  \brief Add a subcompound attached by its inboard bond to an existing bond center
     * 
     * \return a reference to this compound object
     */
    Compound& bondCompound(
        const Compound::Name& n, ///< name for the new subcompound from the viewpoint of the parent Compound
        const Compound& c, ///< the new subcompound to attach to this compound (a copy of the subcompound will be attached)
        const BondCenterPathName& parentBondName, ///< name of the bond center on the parent Compound to which the new subcompound will be attached
        mdunits::Length distance, ///< default bond length in nanometers of new bond connecting new subcompound to parent Compound
        Angle dihedral = 180*Deg2Rad, ///< default dihedral angle about new bond, with respect to the first other bond center on each atom
        BondMobility::Mobility mobility = BondMobility::Default ///< allowed motion in new bond
        );

    /** 
     *  Add a subcompound attached by its inboard bond to an existing bond center
     * 
     *  Shorter version uses default bond length and dihedral angle,
     *  which must have been predefined in the parent *or* the child BondCenter,
     *  but not both.
     * 
     * \return a reference to this compound object
     */
    Compound& bondCompound(
        const Compound::Name& n, ///< name for the new subcompound from the viewpoint of the parent Compound
        const Compound& c,  ///< the new subcompound to attach to this compound (a copy of the subcompound will be attached)
        const BondCenterPathName& parentBondName ///< name of the bond center on the parent Compound to which the new subcompound will be attached
        );

    /** 
     *  \brief Add a subcompound attached by its inboard bond to an existing bond center
     *  
     *  Shorter version uses default bond length and dihedral angle,
     *  which must have been predefined in the parent *or* the child BondCenter,
     *  but not both.
     *
     * \return a reference to this compound object
     */
    Compound& bondCompound(
        const Compound::Name& n, ///< name for the new subcompound from the viewpoint of the parent Compound
        const Compound& c, ///< the new subcompound to attach to this compound (a copy of the subcompound will be attached)
        const BondCenterPathName& parentBondName,  ///< name of the bond center on the parent Compound to which the new subcompound will be attached
        BondMobility::Mobility mobility ///< type of motion permitted in the bond connecting parent to new subcompound
        );

    /** 
     *  setInboardBondCenter assigns special status to a bond center.
     *  There can be at most one inboard bond center in a Compound.
     *  Only an inboard bond center can be used to bond to a parent compound
     * 
     * \return a reference to this compound object
     */
    Compound& setInboardBondCenter(
        const Compound::BondCenterName& centerName ///< name of existing BondCenter to use as inboard BondCenter
        );

    /** Make so that this compound can no longer be a child to the geometry of another compound
     *  Raises an error if the inboard bond center is already bonded.
     * 
     * \return a reference to this compound object
     */
    Compound& convertInboardBondCenterToOutboard();

    /** Assign the very first bond center on a particular atom
     *  \return a reference to this Compound object
     *
     * \return a reference to this compound object
     */
    Compound& addFirstBondCenter(
        const Compound::BondCenterName& centerName, ///< name for the new BondCenter, must be unique within the Compound
        const Compound::AtomPathName& atomName ///< name of the Atom to attach the bond center to
        );

    /**
     * Place a second bond center on an atom, placed a particular
     * angle away from the first bond center.
     *
     * \return a reference to this compound object
     */
    Compound& addSecondBondCenter(
        const Compound::BondCenterName& centerName, ///< name for the new BondCenter, must be unique within the Compound
        const Compound::AtomName& atomName, ///< name of the Atom to attach the bond center to
        Angle bondAngle1 ///< default bond bend angle relating the first two bond centers of the atom
        );

    /**
     * Places a third or later bond center on an atom, in the same
     * plane as the first two bond centers.
     * The angle represents the angular distance from the first
     * bond center on the atom.
     * The second angle is used to determine whether the new bond center is on the same
     * side of the plane as the second bond center.
     *
     * \return a reference to this compound object
     */
    Compound& addPlanarBondCenter(
        const Compound::BondCenterName& centerName, ///< name for the new BondCenter, must be unique within the Compound
        const Compound::AtomName& atomName, ///< name of the Atom to attach the bond center to
        Angle bondAngle1, ///< default bond bend angle relating the first bond center to the new bond center
        Angle bondAngle2 ///< APPROXIMATE default bond bend angle relating the second bond center to the new bond center
        );

    /**
     * Place a third or later bond center on an atom, defined
     * by two angular displacements from the first and second
     * bond centers, respectively.  Placed the new bond center
     * such that the directions of the first, second, and new
     * bond centers form a right-handed geometry.
     *
     * \return a reference to this compound object
     */
    Compound& addRightHandedBondCenter(
        const Compound::BondCenterName& centerName, ///< name for the new BondCenter, must be unique within the Compound
        const Compound::AtomName& atomName, ///< name of the Atom to attach the bond center to
        Angle bondAngle1, ///< default bond bend angle relating the first bond center to the new bond center
        Angle bondAngle2 ///< default bond bend angle relating the second bond center to the new bond center
        );

    /**
     * Place a third or later bond center on an atom, defined
     * by two angular displacements from the first and second
     * bond centers, respectively.  Placed the new bond center
     * such that the directions of the first, second, and new
     * bond centers form a left-handed geometry.
     *
     * \return a reference to this compound object
     */
    Compound& addLeftHandedBondCenter(
        const Compound::BondCenterName& centerName, ///< name for the new BondCenter, must be unique within the Compound
        const Compound::AtomName& atomName, ///< name of the Atom to attach the bond center to
        Angle bondAngle1, ///< default bond bend angle relating the first bond center to the new bond center
        Angle bondAngle2 ///< default bond bend angle relating the second bond center to the new bond center
        );

    /**
     * Adds a covalent bond that is not part of the main
     * tree-topology of bonds in a Compound.  Such bonds are
     * necessary when there are ring structures.  These bonds
     * will not correspond to MobilizedBodies when the molecule
     * is realized as a simbody multibody model, and thus the
     * bond lengths, bond angles, dihedral angles, and mobility for a ring
     * closing bond might not be enforced in the simulation.
     *
     * \return a reference to this compound object
     */
    Compound& addRingClosingBond(
        const Compound::BondCenterPathName& centerName1, ///< name of the first existing bond center in the new bond
        const Compound::BondCenterPathName& centerName2, ///< name of the other existing bond center in the new bond
        mdunits::Length bondLength, ///< default bond length in nanometers of the new bond (might have no effect)
        Angle dihedral = 180*Deg2Rad, ///< default dihedral angle about new bond, with respect to the first other bond center on each atom
        BondMobility::Mobility mobility = BondMobility::Default ///< type of motion permitted in the bond connecting parent to new subcompound
        );

    /**
     * Adds a covalent bond that is not part of the main
     * tree-topology of bonds in a molecule.  Such bonds are
     * necessary when there are ring structures.  These bonds
     * will not correspond to MobilizedBodies when the molecule
     * is realized as a simbody multibody model, and thus the
     * bond lengths, bond angles, and dihedral angles for a ring
     * closing bond might not be enforced in the simulation.
     *
     * \return a reference to this compound object
     */
    Compound& addRingClosingBond(
        const Compound::BondCenterPathName& centerName1, ///< name of the first existing bond center in the new bond
        const Compound::BondCenterPathName& centerName2 ///< name of the other existing bond center in the new bond
        );

    /// \return read-only reference to a subcompound of this Compound
    //const Compound& getSubcompound(Compound::Index ///< integer index of subcompound
    //    ) const;

    /// \return mutable reference to a subcompound of this Compound
    //Compound& updSubcompound(Compound::Index ///< integer index of subcompound
    //    );

    /// \return read-only reference to a subcompound of this Compound
    const Compound& getSubcompound(const Compound::Name& subcompoundName ///< name of subcompound
        ) const;

    /// \return mutable reference to a subcompound of this Compound
    Compound& updSubcompound(const Compound::Name& subcompoundName ///< name of subcompound
        );

    /// @}
    // end topology methods


    /// \name Molecular geometry methods
    /// @{

    /**
     * Compute atom location in local Compound frame
     * 
     * \return default (initial) location of atom in orthogonal nanometers with respect to Ground frame
     */
    Vec3 calcDefaultAtomLocationInGroundFrame(const AtomPathName& atomName ///< name of the atom from viewpoint of Compound
        ) const;

    Vec3 calcDefaultAtomLocationInCompoundFrame(const AtomPathName& atomName ///< name of the atom from viewpoint of Compound
        ) const;

    /** 
     *  \brief Sets default (initial) bond length of current or future Bond using this Compound's inboard BondCenter
     *
     *  Default inboard bond length and dihedral should only be set
     *  when the future bonding partners of this compound are 
     *  restricted to a particular atom type
     * 
     * \return a reference to this compound object
     */
    Compound& setDefaultInboardBondLength(mdunits::Length ///< bond length in nanometers
        );

    /**
     * \brief Stores a default (initial) dihedral angle in the inboard BondCenter of this Compound.
     *
     * If exactly one of the two BondCenters that are joined to form a Bond has default geometry
     * set, then the bondCompound() lacking explicit geometry arguments may be used to create the Bond.
     *
     * \return a reference to this compound object
     */
    Compound& setDefaultInboardDihedralAngle(Angle ///< dihedral angle in radians
        );

    /** 
     *  \brief Sets a default(initial) bond angle defined by three atoms
     *
     *  \warning setDefaultBondAngle only works if one of the two bond centers is bond center 1 or 2
     *  \todo - remove this restriction
     * 
     * \return a reference to this compound object
     */
    Compound& setDefaultBondAngle(
        Angle angle, ///< new default bond angle in radians
        const AtomPathName& atom1, ///< name of first atom defining bond angle
        const AtomPathName& atom2, ///< name of middle atom defining bond angle
        const AtomPathName& atom3 ///< name of third atom defining bond angle
        );

    /**
     * \brief Sets a default (inital) bond length defined by two atoms
     *
     * \return a reference to this compound object
     */
    Compound& setDefaultBondLength(
        mdunits::Length length, ///< default bond length in nanometers
        const AtomPathName& atom1, ///< name of first atom in bond
        const AtomPathName& atom2 ///< name of second atom in bond
        );

    /**
     * \brief Sets a default (initial) dihedral angle of a previously named dihedral
     *
     * \return a reference to this compound object
     */
    Compound& setDefaultDihedralAngle(
        const DihedralName& dihedralName, ///< name of predefined dihedral angle with respect to this Compound
        Angle angleInRadians ///< new dihedral angle in radians
        );

    Compound& setDefaultDihedralAngle(
    			Angle angle, 
    			Compound::AtomIndex atom1,
    			Compound::AtomIndex atom2,
    			Compound::AtomIndex atom3,
    			Compound::AtomIndex atom4
    			);
    
    Compound& setDefaultDihedralAngle(
    			Angle angle, 
    			const Compound::AtomName& atom1,
    			const Compound::AtomName& atom2,
    			const Compound::AtomName& atom3,
    			const Compound::AtomName& atom4
    			);
    
   /**
     * \brief Computes default (initial) dihedral angle of a previously named dihedral
     *
     * \return dihedral angle in radians with respect to first other bond centers of bonded atoms
     */
    Angle calcDefaultDihedralAngle(
        const DihedralName& dihedralName ///< name of predefined dihedral angle with respect to this Compound
        ) const;

    /**
     *  \brief Sets dynamic dihedral angle of a previously named dihedral
     *
     *  The dihedral angle must be an adjustable parameter of the dynamic model.
     *
     *  \return a reference to this compound object
     */
    Compound& setDihedralAngle(
        State& state, ///< simbody State object representing the current configuration
        const DihedralName& dihedralName, ///< name of predefined dihedral angle with respect to this Compound
        Angle ///< dihedral angle in radians
        );

    /**
     * \brief Computes dynamic dihedral angle of a previously named dihedral
     *
     * \return dihedral angle in radians with respect to first other bond centers of bonded atoms
     */
    Angle calcDihedralAngle(
        const State& state, ///< simbody State object representing the current configuration
        const DihedralName& dihedralName ///< name of predefined dihedral angle with respect to this Compound
        ) const;

    /// \return location and orientation of Atom with respect to this Compound
    Transform calcDefaultAtomFrameInCompoundFrame(Compound::AtomIndex ///< integer index of Atom with respect to this Compound
        ) const;

    /**
     * \return location in orthogonal nanometers
     */
    Vec3 calcAtomLocationInGroundFrame(
        const State& state, ///< simbody State representing current configuration
        Compound::AtomIndex atomId ///< integer index of Atom with respect to this Compound
        ) const;

    Vec3 calcAtomLocationInCompoundFrame(
        const State& state, ///< simbody State representing current configuration
        Compound::AtomIndex atomId ///< integer index of Atom with respect to this Compound
        ) const;

    /**
     * \return vector velocity in nanometers per picosecond
     */
    Vec3 calcAtomVelocityInGroundFrame(
        const State& state, ///< simbody State representing current configuration
        Compound::AtomIndex atomId  ///< integer index of Atom with respect to this Compound
        ) const;

    /**
     * \return vector acceleration in nanometers per picosecond squared
     */
    Vec3 calcAtomAccelerationInGroundFrame(
        const State& state, ///< simbody State representing current configuration
        Compound::AtomIndex atomId   ///< integer index of Atom with respect to this Compound
        ) const;

    /**
     * \return default (initial) location and orientation of a subcompound with respect to this Compound
     */
    Transform getSubcompoundFrameInParentFrame(
        const Compound::Name& subcompoundName ///< name of subcompound with respect to this Compound
        ) const;

    /**
     * \brief Create a mapping between this Compound and atom locations in a PdbStructure
     *
     * For use in the process of importing a configuration from a PDB file
     */
    virtual AtomTargetLocations createAtomTargets(const class PdbStructure& targetStructure,const bool guessCoordinates = false ) const;
    virtual AtomTargetLocations createAtomTargets(const class PdbChain& targetChain,const bool guessCoordinates  = false) const;
    virtual AtomTargetLocations createAtomTargets(const class PdbResidue& targetResidue,const bool guessCoordinates  = false) const;

    /**
     * \brief Adjust stereochemistry about chiral atoms to match that seen in a set of atomic locations
     *
     * Choose a small value of planarityTolerance parameter to break the planarity of out-of-plane atoms from the target set 
     *
     * \return a reference to this compound
     */
    Compound& matchDefaultAtomChirality(
            const AtomTargetLocations& atomTargets, ///< another set of atom locations to match chirality of
            Angle planarityTolerance = 90.0 * Deg2Rad, ///< break planarity of BondCenters more than this angle out of plane (radians)
            bool flipAll = true ///< whether to invert each mismatched BondCenter, or all at once
            );
    Compound& matchDefaultBondLengths(const AtomTargetLocations& atomTargets);
    Compound& matchDefaultBondAngles(const AtomTargetLocations& atomTargets);

    enum PlanarBondMatchingPolicy {KeepPlanarBonds, DistortPlanarBonds, FlipPlanarBonds};
    Compound& matchDefaultDihedralAngles(
        const AtomTargetLocations& atomTargets,///< set of atom locations to match dihedrals of
        PlanarBondMatchingPolicy policy = FlipPlanarBonds ///< whether to keep ideal torsions on bonds between planar atoms
        );

    Compound& matchDefaultTopLevelTransform(const AtomTargetLocations& atomTargets);
    TransformAndResidual getTransformAndResidual(const Compound::AtomTargetLocations& atomTargets) const;
    
    /// Adjust internal coordinates to match a collection of atom targets.
    Compound& matchDefaultConfiguration(
            const AtomTargetLocations& atomTargets, 
            MatchStratagem matchStratagem = Match_Exact,
            const bool useObservedPointFitter = true,
            const Real minimizerTolerance = 150.0
            ) 
    {
        if (matchStratagem == Compound::Match_TopologyOnly)
            // do nothing, topology is already there
            return *this;
        
        else if (matchStratagem == Compound::Match_Exact) 
        {
            // low tolerance breaks planarity just about everywhere
            matchDefaultAtomChirality(atomTargets, 0.01, false);

            matchDefaultBondLengths(atomTargets);
            matchDefaultBondAngles(atomTargets);
            
            // Set dihedral angles even when bonded atoms are planar
            matchDefaultDihedralAngles(atomTargets, Compound::DistortPlanarBonds);

            matchDefaultTopLevelTransform(atomTargets);
            
            // No further optimization should be needed
            
            return *this;
        }
        
        else if (matchStratagem == Compound::Match_Idealized) 
        {
            // Break planarity constraint only when input atoms are 90 degrees out of plane
            matchDefaultAtomChirality(atomTargets, 90*Deg2Rad, true);
            
            // Avoid setting non zero/180 dihedral angles on bonded planar atoms
            matchDefaultDihedralAngles(atomTargets, Compound::FlipPlanarBonds);
            
            matchDefaultTopLevelTransform(atomTargets);

            // At this point the structure match is approximate and may have
            // many bad contacts.  Optimization is needed
            // TODO - ObservedPointFitter on internal coordinates
            fitDefaultConfiguration(atomTargets, 0.005,useObservedPointFitter,minimizerTolerance);
            
            return *this;
        }
            
        else assert(false); // unknown match stratagem

        return *this;
    }

    /// Optimize adjustable degrees of freedom to best match atom targets
    Compound& fitDefaultConfiguration(
            const AtomTargetLocations& atomTargets,
            SimTK::Real targetRms,
            const bool useObservedPointFitter = true,
            const Real minimizerTolerance = 150.0
            );
    
    /// Write current default(initial) Compound configuration into a PdbChain object
    const Compound& populateDefaultPdbChain(
        class PdbChain&, 
        int& defaultNextResidueNumber,
        const Transform& transform = Transform()) const;

    /// Write dynamic Compound configuration into a PdbChain object
    const Compound& populatePdbChain(
        const State& state, 
        class PdbChain&, 
        int& defaultNextResidueNumber,
        const Transform& transform = Transform()) const;

    /**
     * \brief Write the default (initial) configuration in Protein Data Bank (PDB) format.
     *
     * Starting with the first atom's serial number as one(1).
     */
    std::ostream& writeDefaultPdb(
        std::ostream& os, ///< output stream to write PDB coordinates to
        const Transform& transform = Transform()  ///< optional change to location and orientation of molecule
        ) const;

   /**
    * /brief This polymorphism takes a char* file name rather than ostream, to save the user a couple of lines of code.
    *
    */

    void writeDefaultPdb(
        const char* outFileName, 
        const Transform& transform
        ) const;

    /** 
     * \brief Write the default (initial) configuration in Protein Data Bank (PDB) format.
     *
     * integer nextAtomSerialNumber reference is incremented within writeDefaultPdb method, so that subsequent
     * calls to this method using the same integer variable will continue numbering
     * atoms where the previous call left off.
     */
    std::ostream& writeDefaultPdb(
        std::ostream& os, ///< output stream to write PDB coordinates to
        int& nextAtomSerialNumber, ///< mutable integer reference containing the next desired atom serial number
        const Transform& transform = Transform() ///< optional change to location and orientation of molecule
        ) const;

    /**
     * \brief Write the dynamic Compound configuration in Protein Data Bank (PDB) format.
     *
     * Starting with the first atom's serial number as one(1).
     */
    std::ostream& writePdb(
        const State& state, ///< simbody state representing the current configuration of the molecule
        std::ostream& os,  ///< output stream to write PDB coordinates to
        const Transform& transform = Transform() ///< optional change to location and orientation of molecule
        ) const;


    /** 
     * \brief Write the dynamic Compound configuration in Protein Data Bank (PDB) format.
     *
     * integer nextAtomSerialNumber reference is incremented within writePdb method, so that subsequent
     * calls to this method using the same integer variable will continue numbering
     * atoms where the previous call left off.
     */
    std::ostream& writePdb(
        const State& state, ///< simbody state representing the current configuration of the molecule
        std::ostream& os,  ///< output stream to write PDB coordinates to
        int& nextAtomSerialNumber,  ///< mutable integer reference containing the next desired atom serial number
        const Transform& transform = Transform() ///< optional change to location and orientation of molecule
        ) const;

    /// @} 
    // end configuration section


    /// \name Compound nomenclature methods
    /// @{

    /**
     * \brief Name a type of Compound
     *
     * Sets the name for this class of compound, e.g. "ethane", not "ethane number 3".
     * The compound name can be important for automatic resolution of atom types for
     * during force field resolution.
     *
     * \return a reference to this compound object
     */
    Compound& setCompoundName(const Name& ///< new name for Compound type
        );


    /// \return name of Compound type
    const Name& getCompoundName() const;

    /**
     * \brief Add an additional name for this class of compound.
     *
     * Adding synonyms can be useful for automatically resolving biotypes, if the compound name in
     * the force field parameter file differs from that used to define the compound.
     * See also setCompoundName()
     * \return a reference to this compound object
     */
    Compound& addCompoundSynonym(const Name& ///< new synonym for Compound type
        );

    /**
     * Returns the most recently assigned name, if any, given to
     * an atom in this compound.  
     */
    const AtomPathName getAtomName(Compound::AtomIndex ///< integer index for Atom with respect to this Compound.
        ) const;

    /// \return chemical element of atom
    const Element& getAtomElement(Compound::AtomIndex ///< integer index for Atom with respect to this Compound.
        ) const;

    const Element& getAtomElement(const Compound::AtomName& ///< name for Atom in the context of this Compound.
        ) const;

    /**
     * \return  reference to this compound object
     */
    Compound& nameAtom(
        const Compound::AtomName& newName, ///< new name for this atom; name is local to this Compound
        const AtomPathName& oldName ///< previous name; can be a subcompound-qualified Atom name
        );

    /**
     * \return a reference to this compound object
     */
    Compound& nameAtom(
        const Compound::AtomName& newName,  ///< new name for this atom; name is local to this Compound
        const AtomPathName& oldName,  ///< previous name; can be a subcompound-qualified Atom name
        BiotypeIndex biotype ///< new Biotype for the atom, with respect to this Compound
        );

    /**
     * \brief Define a named dihedral angle using four atoms
     *
     * The offset parameter is used in cases where the dihedral angle definition differs from 
     * the IUPAC definition of dihedral angles using four atoms.
     * For example the definition of protein phi and psi angles is 180 degrees offset from the 
     * standard four-atom dihedral angle definition.
     * 
     * \return a reference to this compound object
     */
    Compound& defineDihedralAngle(
        const Compound::DihedralName& angleName, ///< unique name for new dihedral angle definition
        const Compound::AtomPathName& atom1, ///< first atom name
        const Compound::AtomPathName& atom2, ///< second atom name
        const Compound::AtomPathName& atom3, ///< third atom name
        const Compound::AtomPathName& atom4, ///< fourth atom name
        Angle offset = 0*Deg2Rad ///< nomenclature offset
        );

    /**
     * \brief Define a named dihedral in terms of two bond centers.
     *
     * Useful in cases where not all four atoms have
     * been placed yet.  The bond centers must be from two atoms that are bonded (atoms 2 and 3 of 4)
     * But are NOT the bond centers that connect the two middle atoms.
     *
     * The offset parameter is used in cases where the dihedral angle definition differs from 
     * the IUPAC definition of dihedral angles using four atoms.
     * For example the definition of protein phi and psi angles is 180 degrees offset from the 
     * standard four-atom dihedral angle definition.
     * 
     * \return a reference to this compound object
     */
    Compound& defineDihedralAngle(
        const Compound::DihedralName& angleName, ///< unique name for new dihedral angle definition
        const Compound::BondCenterPathName& bond1, ///< first bond center, connecting atom 2 to atom 1
        const Compound::BondCenterPathName& bond2, ///< second bond center, connecting atom 3 to atom 4
        Angle offset = 0*Deg2Rad ///< nomenclature offset
        );

    /**
     * \brief Set value to populate "residue number" field for PDB file output
     *
     * \return a reference to this compound object
     */
    Compound& setPdbResidueNumber(int resNum);

    /// \return integer that would appear in "residue number" field for this Compound in a PDB file
    int getPdbResidueNumber() const;

    /**
     * \brief Set string to populate "residue type" field for PDB file output
     *
     * \return a reference to this compound object
     */
    Compound& setPdbResidueName(const String& resName);

    /// \return string that would be used to populate "residue type" field for this Compound in a PDB file
    const String& getPdbResidueName() const;

    /**
     * \brief Set character to populate "chain id" field for PDB file output
     *
     * \return a reference to this compound object
     */
    Compound& setPdbChainId(char chainId ///< Protein Data Bank "chain Id" for this Compound
        );

    /// \return character that would populate "chain id" field for PDB file output for this Compound
    char getPdbChainId() const;

    /**
     * \brief Define a local name for a particular BondCenter in this Compound
     *
     * \return a reference to this compound object
     */
    Compound& nameBondCenter(
        const Compound::BondCenterName& newName, ///< new local name for bond center
        const BondCenterPathName& oldName ///< previous name, can be a subcompound-qualified name
        );

    /**
     * \brief Convenience method to locally import all of the atom names that a particular subcompound uses.
     *
     * \return a reference to this compound object
     */
    Compound& inheritAtomNames(const Compound::Name& ///< name of subcompound with having desired atom names
        );
    
    Compound& inheritCompoundSynonyms(const Compound& otherCompound);
    
    /**
     * \brief Convenience method to locally import all of the BondCenter names that a particular subcompound uses.
     *
     * \return a reference to this compound object
     */
    Compound& inheritBondCenterNames(const Compound::Name& ///< name of subcompound with having desired BondCenter names
        );

    /// \return integer index of named atom, with repect to this Compound
    Compound::AtomIndex getAtomIndex(const Compound::AtomPathName& ///< atom name with respect to this Compound
        ) const;

    /// @}
    // end nomenclature methods


    /// \name Compound simulation methods
    /// @{

    /** 
     * \brief Add this Compound, including all of its subcompounds, to a CompoundSystem, for simulation etc.
     *
     * This Compound will become a top-level object with the CompoundSystem.
     */
    //void setCompoundSystem(CompoundSystem& system, ///< The CompoundSystem that will take ownership of this Compound.
    //                       Compound::Index);       ///< The index within the new owner CompoundSystem that this Compound will have.

    void setMultibodySystem(MultibodySystem& system);

    /** 
     *  Override the default rotatability of a bond
     * 
     * \return a reference to this compound object
     */
    Compound& setBondMobility(
        BondMobility::Mobility mobility, ///< the new allowed motion of the bond
        const AtomPathName& atom1, ///< the name of the first atom defining the Bond
        const AtomPathName& atom2 ///< the name of the second atom defining the Bond
        );

    /** 
     *  Override the default rotatability of a bond
     * 
     * \return a reference to this compound object
     */
    Compound& setBondMobility(
        BondMobility::Mobility mobility,  ///< the new allowed motion of the bond
        Compound::BondIndex bondIndex ///< the integer index of the Bond in the context of this Compound
        );

    /**
     * \brief Set BondMobility for every bond in the Compound.
     * 
     *
     * \return A reference to this Compound  
     *
     */	
    Compound  & setCompoundBondMobility(BondMobility::Mobility mobility);

    /**
     * \brief get the simbody MobilizedBody to which a particular atom is attached
     *
     * Requires that this Compound has already been modeled in a CompoundSystem
     *
     * \return the integer index of the MobilizedBody in the CompoundSystem
     */
    MobilizedBodyIndex getAtomMobilizedBodyIndex(Compound::AtomIndex ///< integer index of the Atom in this Compound context.
        ) const;

    /**
     * \brief get the location of an Atom in the frame of the simbody MobilizedBody to which the Atom is attached
     *
     * Requires that this Compound has already been modeled in a CompoundSystem
     *
     * \return location in orthogonal nanometers with respect to the rigid body to which the atom is attached
     */
    Vec3 getAtomLocationInMobilizedBodyFrame(Compound::AtomIndex ///< integer index of the Atom in this Compound context.
        ) const;

    /**
     * \brief define the Biotype for an Atom in this Compound
     *
     * \return a reference to this compound object
     */
    Compound& setBiotypeIndex(
        const Compound::AtomPathName& atomName, ///< name of the atom in the context of this Compound
        BiotypeIndex biotype ///< integer index of an existing Biotype known to the Biotype class
        );
    
    Compound& setAtomBiotype(
            const Compound::AtomPathName& atomName,
            const String& biotypeResidueName,
            const String& biotypeAtomName,
            SimTK::Ordinality::Residue ordinality = SimTK::Ordinality::Any
            )
    {
        Compound::AtomIndex atomIndex = getAtomIndex(atomName);
        
        if ( ! Biotype::exists(biotypeResidueName, biotypeAtomName, ordinality) )
			Biotype::defineBiotype(
                getAtomElement(atomIndex),
                getNumBondCenters(atomIndex), 
                biotypeResidueName, 
                biotypeAtomName, 
                ordinality
            );
        
        const Biotype& biotype = Biotype::get(biotypeResidueName, biotypeAtomName, ordinality);
        
        assert( biotype.getElement() == getAtomElement(atomIndex) );
        assert( biotype.getValence() == getNumBondCenters(atomIndex) );
        
        BiotypeIndex biotypeIndex = biotype.getIndex();
        setBiotypeIndex(atomName, biotypeIndex);

        return *this;
    }

    /// \return the biotype assigned to an Atom in this Compound
    BiotypeIndex getAtomBiotypeIndex(Compound::AtomIndex ///< integer index of an Atom in this Compound
        ) const;

    // RUNTIME INTERFACE

    /**
     * \brief get DuMMForceFieldSubsystem atom index corresponding to an atom in this Compound
     *
     * Requires that this Compound has already been modeled in a CompoundSystem
     *
     * \return an AtomIndex in the DuMMForceFieldSubystem
     */
    DuMM::AtomIndex getDuMMAtomIndex(Compound::AtomIndex) const;

    /// @}
    // end simulation methods

    /// Set the orientation and location of this Compound
    Compound& setTopLevelTransform(const Transform& transform);

    /// \return the orientation and location of this Compound
    const Transform& getTopLevelTransform() const;



protected:

    /**
     * \brief Stores relationship between a Compound Atom and an Atom defined in a DuMMForcefieldSubsystem
     */
    void setDuMMAtomIndex(
        Compound::AtomIndex, ///< integer index of an existing Atom in this Compound
        DuMM::AtomIndex ///< integer index of an Atom in a DuMMForceFieldSubsystem
        );

    explicit Compound(CompoundRep* ip);
    friend class CompoundSystem;

private:
    // OBSOLETE: use getNumAtoms()
    int getNAtoms() const {return getNumAtoms();}
    // OBSOLETE: use getNumBondCenters()
    size_t getNBondCenters() const {return getNumBondCenters();}
    // OBSOLETE: use getNumBondCenters()
    size_t getNBondCenters(Compound::AtomIndex atomIndex) const {return getNumBondCenters(atomIndex);}
    // OBSOLETE: use getNumBonds()
    size_t getNBonds() const {return getNumBonds();}
};


/**
 * \brief Dump debugging information about compound structure to a stream.
 *
 *  This method does NOT produce PDB files or anything like it.  It is used for 
 *  debugging the internal structure of instances of the Compound class and is thus
 *  not intended for typical API client use.
 */
SimTK_MOLMODEL_EXPORT std::ostream& operator<<(std::ostream& o, const Compound& c);

/**
 * Base class for single-atom Compound building blocks
 *
 * There is not an explicit Atom class in the public Molmodel API.  The 
 * Compound::SingleAtom class is intended to provide single atom Compounds
 * that can be linked together during the construction of Molecules.  Many 
 * of the predefined Molecule types in the Molmodel API
 * use SingleAtoms in their construction.
 */
class  Compound::SingleAtom : public Compound {
public:
    SingleAtom(
        const Compound::AtomName& atomName, ///< name for new atom
        const Element& element ///< chemical element of new atom
        ) 
    {
        setBaseAtom( atomName, element );

        setCompoundName("SingleAtom"); // should be overridden by derived class constructors
    }
};


/**
 * \brief Base class for atoms with exaclty one covalent bond partner.
 *
 * Lone bond center is called "bond".
 *
 * About the use of "Univalent" instead of "Monovalent":
 * "valence" is derived from Latin, and thus should properly be combined
 * with other Latin roots, though both Greek and Latin forms are observed
 * for numbered valences.
 * Latin naming convention is univalent, bivalent, trivalent, quadrivalent,
 * quinqivalent, sexivalent, septivalent, octivalent.
 * Greek-Latin mongrel words are also sometimes used: monovalent, divalent,
 * trivalent, tetravalent, pentavalent, hexavalent, heptavalent, octavalent.
 */
class  UnivalentAtom : public Compound::SingleAtom {
public:
    UnivalentAtom(
        const Compound::AtomName& atomName, ///< name for new atom
        const Element& element ///< chemical element for new atom
        ) 
        : Compound::SingleAtom(atomName, element)
    {
        addFirstBondCenter( "bond", atomName);
        setInboardBondCenter("bond");

        setCompoundName("UnivalentAtom"); // should be overridden by derived class constructors
    }
};


/**
 *  \brief Base class for atoms having exactly two covalent bonds.
 *
 *  Bond centers are named "bond1" and "bond2"
 */
class  BivalentAtom : public Compound::SingleAtom {
public:
    BivalentAtom(
        const Compound::AtomName& atomName,  ///< name for new atom
        const Element& element,  ///< chemical element for new atom
        Angle angle = 180*Deg2Rad ///< default (initial) bond angle between new atom's two BondCenters
        ) 
        : Compound::SingleAtom(atomName, element)
    {        
        // BondCenter1 dihedral will be relative to BondCenter2
        addFirstBondCenter( "bond1", atomName);
        // conversely, bond center 2 dihedral is relative to bond center 1
        addSecondBondCenter( "bond2", atomName, angle);
        setInboardBondCenter("bond1"); // without loss of generality

        setCompoundName("BivalentAtom"); // should be overridden by derived class constructors
    }
};


/**
 *  \brief Base class for atoms having exactly three covalent bonds.
 *
 *  Initial default configuration is planar.
 *  Bond centers are named "bond1", "bond2", and "bond3"
 */
class  TrivalentAtom : public Compound::SingleAtom {
public:
    TrivalentAtom(
        const Compound::AtomName& atomName,   ///< name for new atom
        const Element& element,   ///< chemical element for new atom
        Angle angle1 = 120*Deg2Rad, ///< angle between first and second BondCenters
        Angle angle2 = 120*Deg2Rad ///< angle between first and third BondCenters
        ) 
        : Compound::SingleAtom(atomName, element)
    {
        // BondCenter1 dihedral will be relative to BondCenter2
        addFirstBondCenter( "bond1", atomName );
        // bond centers 2 and 3 dihedrals relative to bond center 1
        addSecondBondCenter( "bond2", atomName,  angle1);
        addPlanarBondCenter( "bond3", atomName, angle2, 360*Deg2Rad - angle1 - angle2);

        // Choice of inboard bond may differ from bond priority - user may change this
        setInboardBondCenter("bond1");

        setCompoundName("TrivalentAtom"); // should be overridden by derived class constructors
    }
};

/**
 *  \brief Base class for atoms having exactly four covalent bonds.
 *
 *  Bond centers are named "bond1", "bond2", "bond3", and "bond4"
 */
class  QuadrivalentAtom : public Compound::SingleAtom {
public:
    QuadrivalentAtom(
        const Compound::AtomName& atomName, ///< name for new atom
        const Element& element ///< chemical element for new atom
        ) 
        : Compound::SingleAtom(atomName, element)
    {
        static const Angle TetrahedralAngle = 109.47 * Deg2Rad;

        // BondCenter1 dihedral will be relative to BondCenter2
        // Using Rotation() constructor that takes x axis and approximate y axis
        addFirstBondCenter( "bond1", atomName );
        // bond centers 2, 3, and 4 dihedrals will be relative to bond1
        addSecondBondCenter( "bond2", atomName, TetrahedralAngle );
        addLeftHandedBondCenter( "bond3", atomName, TetrahedralAngle, TetrahedralAngle );
        addRightHandedBondCenter( "bond4", atomName, TetrahedralAngle, TetrahedralAngle );

        // Choice of inboard bond may differ from bond priority - user may change this
        setInboardBondCenter("bond1");

        // default dihedral should be left blank, so it could be
        // defined from the other side.
        // the ultimate default will be 180 degrees anyway.
        // setDefaultInboardDihedralAngle(180*Deg2Rad);

        setCompoundName("QuadrivalentAtom"); // should be overridden by derived class constructors
    }
        
    QuadrivalentAtom(
        const Compound::AtomName& atomName, ///< name for new atom
        const Element& element, ///< chemical element for new atom
        Angle bond12Angle, ///< bond angle in radians between bond center 1 and bond center 2
        Angle bond13Angle, ///< bond angle in radians between bond center 1 and bond center 3
        Angle bond14Angle, ///< bond angle in radians between bond center 1 and bond center 4
        Angle dihedral3, ///< dihedral angle of bond center 3, relative to dihedral of bond center 2
        Angle dihedral4 ///< dihedral angle of bond center 4, relative to dihedral of bond center 2
        ) 
        : Compound::SingleAtom(atomName, element)
    {
        // BondCenter1 dihedral will be relative to BondCenter2
        // Using Rotation() constructor that takes x axis and approximate y axis
        addFirstBondCenter( "bond1", atomName );
        // bond centers 2, 3, and 4 dihedrals will be relative to bond1
        addSecondBondCenter( "bond2", atomName, bond12Angle );
        
        // Compute bond23Angle and bond24Angle
        SimTK::Vec3 v1(0,0,1);
        SimTK::Vec3 v2 = SimTK::Rotation(bond12Angle, ZAxis) * SimTK::Vec3(0,0,1);
        SimTK::Vec3 v3 = SimTK::Rotation(dihedral3, YAxis) * SimTK::Rotation(bond12Angle, ZAxis) * SimTK::Vec3(0,0,1);
        SimTK::Vec3 v4 = SimTK::Rotation(dihedral4, YAxis) * SimTK::Rotation(bond12Angle, ZAxis) * SimTK::Vec3(0,0,1);
        
        Angle bond23Angle = std::acos(SimTK::dot(v2, v3));
        Angle bond24Angle = std::acos(SimTK::dot(v2, v4));
        
        // restrict angle range to (-Pi,Pi)
        while (dihedral3 >   Pi) dihedral3 -= 2.0*Pi;
        while (dihedral3 <= -Pi) dihedral3 += 2.0*Pi;
        while (dihedral4 >   Pi) dihedral4 -= 2.0*Pi;
        while (dihedral4 <= -Pi) dihedral4 += 2.0*Pi;
        
        if (dihedral3 > 0)
        	addLeftHandedBondCenter( "bond3", atomName, bond13Angle, bond23Angle );
        else
        	addRightHandedBondCenter( "bond3", atomName, bond13Angle, bond23Angle );
        
        if (dihedral4 > 0)
        	addLeftHandedBondCenter( "bond4", atomName, bond14Angle, bond24Angle );
        else 
        	addRightHandedBondCenter( "bond4", atomName, bond14Angle, bond24Angle );

        // Choice of inboard bond may differ from bond priority - user may change this
        setInboardBondCenter("bond1");

        // default dihedral should be left blank, so it could be
        // defined from the other side.
        // the ultimate default will be 180 degrees anyway.
        // setDefaultInboardDihedralAngle(180*Deg2Rad);

        setCompoundName("QuadrivalentAtom"); // should be overridden by derived class constructors
    }
};


/**
 *  AliphaticHydrogen is a hydrogen atom for bonding to AliphaticCarbon atoms (see below)
 */
class  AliphaticHydrogen : public UnivalentAtom {
public:
    explicit AliphaticHydrogen(const AtomName& atomName = "H" ///< name for new atom, defaults to "H"
        ) 
        : UnivalentAtom(atomName, Element::Hydrogen())
    {
        setDefaultInboardBondLength(0.1112); // for bonding to aliphatic carbon

        setCompoundName("AliphaticHydrogenAtom");
    }
};


/**
 * \brief AliphaticCarbon is a tetrahedral sp3 carbon atom for bonding to four other things
 *
 * The default geometry is perfectly terahedral, which is not bad for many uses.
 */
class  AliphaticCarbon : public QuadrivalentAtom {
public:
    explicit AliphaticCarbon(const AtomName& atomName = "C" ///< name for new atom, defaults to "C"
        ) 
        : QuadrivalentAtom(atomName, Element::Carbon()) 
    {
        // In case this bonds to another aliphatic carbon
        setDefaultInboardBondLength(0.15620); // for bonding to another aliphatic carbon

        setCompoundName("AliphaticCarbon");
    }
};


/**
 * MethyleneGroup is -CH2- group for bonding to aliphatic carbon and a second something
 */
class  MethyleneGroup : public Compound {
public:
    MethyleneGroup() {
        static const mdunits::Length C_Hdistance = 0.1112; // nanometers

        setBaseAtom(AliphaticCarbon("C"));
        bondAtom(AliphaticHydrogen("H1"), "C/bond3");
        bondAtom(AliphaticHydrogen("H2"), "C/bond4");
        nameBondCenter("bond1", "C/bond1");
        nameBondCenter("bond2", "C/bond2");

        setDefaultInboardBondLength(0.15620); // for bonding to another aliphatic carbon

        setCompoundName("MethyleneGroup");
    }
};


/**
 *  MethylGroup is CH3 for attaching to aliphatic carbon
 */
class  MethylGroup : public Compound {
public:
    MethylGroup() {
        setBaseAtom(AliphaticCarbon("C"));
        bondAtom(AliphaticHydrogen("H1"), "C/bond2");
        bondAtom(AliphaticHydrogen("H2"), "C/bond3");
        bondAtom(AliphaticHydrogen("H3"), "C/bond4");
        nameBondCenter("bond", "C/bond1");

        setDefaultInboardBondLength(0.15620); // for bonding to another aliphatic carbon
        setCompoundName("MethylGroup");
    }
};

/**
 * \brief Two atom C-H group in six membered aromatic rings.
 *
 * Was intended to be used as a building block for aromatic rings like benzene or phenylalanine,
 * but it turns out to be just as easy to build those using one atom at a time
 */
class  AromaticSixMemberedCHGroup : public Compound {
public:
    AromaticSixMemberedCHGroup() {
        static const mdunits::Length C_Hdistance = 0.1080; // in nanometers, from tinker amber99.dat
        static const mdunits::Length C_Cdistance = 0.1400; // in nanometers, from tinker amber99.dat

        setBaseAtom( TrivalentAtom("C", Element::Carbon(), 120*Deg2Rad, 120*Deg2Rad) );
        bondAtom( UnivalentAtom("H", Element::Hydrogen()), "C/bond3", C_Hdistance);

        // remember to change default length for bonding other other things
        setDefaultInboardBondLength(C_Cdistance); // for bonding to other aromatic CH

        setCompoundName("AromaticSixMemberedCHGroup");
    }
};


/**
 * AlcoholOHGroup is OH group for attachment to aliphatic carbon
 */
class  AlcoholOHGroup : public Compound {
public:
    AlcoholOHGroup() {
        static const mdunits::Length O_Hdistance = 0.0960; // nanometers

        setBaseAtom( BivalentAtom("O", Element::Oxygen(), 108.50 * Deg2Rad) );

        bondAtom( UnivalentAtom("H", Element::Hydrogen()), "O/bond2", O_Hdistance );
        nameBondCenter("bond", "O/bond1");

        // In case of bonding to aliphatic carbon
        setDefaultInboardBondLength(0.1410); // for bonding to aliphatic carbon

        setCompoundName("AlcoholOHGroup");
    }
};


/**
 * PrimaryAmineGroup is NH3<sup>+</sup> for attachment to tetrahedral carbon
 */
class  PrimaryAmineGroup : public Compound {
public:
    PrimaryAmineGroup() {
        static const mdunits::Length H_Ndistance = 0.1010; // nanometers
        static const mdunits::Length N_Cdistance = 0.1471; // in nanometers, for bonding to aliphatic carbon

        setBaseAtom( QuadrivalentAtom("N", Element::Nitrogen()) );
        bondAtom( UnivalentAtom("H1", Element::Hydrogen()), "N/bond2", H_Ndistance);
        bondAtom( UnivalentAtom("H2", Element::Hydrogen()), "N/bond3", H_Ndistance);
        bondAtom( UnivalentAtom("H3", Element::Hydrogen()), "N/bond4", H_Ndistance);

        setDefaultInboardBondLength(N_Cdistance); // for bonding to aliphatic carbon
    }
};


/**
 * CaboxylateGroup is COO<sup>-</sup> for attachment to tetrahedral carbon
 */
class  CarboxylateGroup : public Compound {
public:
    CarboxylateGroup() {
        // parameters from Tinker amber99.param
        static const mdunits::Length O_Cdistance = 0.1250; // nanometers
        static const mdunits::Length C_CtDistance = 0.15220; // nanometers
        static const Angle O_C_Oangle = 126*Deg2Rad;

        // we actually use the inboard-C-O angle during construction
        static const Angle O_C_CtAngle = 180*Deg2Rad - 0.5*O_C_Oangle;

        setBaseAtom( TrivalentAtom("C", Element::Carbon(), O_C_CtAngle, O_C_CtAngle) );
        bondAtom( UnivalentAtom("O1", Element::Oxygen()), "C/bond2", O_Cdistance);
        bondAtom( UnivalentAtom("O2", Element::Oxygen()), "C/bond3", O_Cdistance);

        setDefaultInboardBondLength(C_CtDistance);

        nameBondCenter("bond", "C/bond1");
    }
};

/**
 * \brief Base class for complete covalently connected molecules.
 *
 * Molecule class is parent class for complete covalently connected molecules
 * Presently it is used to document the intention that derived classes
 * are single covalently connected entities
 *
 * Derives from Compound.
 * (though that may not be obvious in the automatically
 * generated API documentation).
 */
class  SimTK_MOLMODEL_EXPORT Molecule : public Compound {
public:
    Molecule();
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Molecule, CompoundRep, Compound);
protected:
    explicit Molecule(CompoundRep* rep);
};


/**
 * \brief The noble gas argon, which does not bond with other atoms
 */
class  Argon : public Molecule {
public:
    Argon() {
        setPdbResidueName("AR ");

        setBaseAtom( "Ar", Biotype::Argon() );

        setCompoundName("Argon");
    }
};


/**
 * \brief The simplest hydrocarbon methane, CH<sub>4</sub>
 */
class  Methane : public Molecule {
public:
    Methane() 
    {
        setBaseCompound("methyl", MethylGroup());
        inheritAtomNames("methyl");

        // Ordinarily, methyl group bonds to aliphatic carbon,
        // and has a default bond length to match.
        // Here we turn off the methyl inboard bond, so the
        // AliphaticHydrogen will, as desired, dictate the bond length
        convertInboardBondCenterToOutboard();

        bondAtom(AliphaticHydrogen("H4"), "methyl/bond", 0.1112);

        setBiotypeIndex( "C", Biotype::MethaneC().getIndex() );
        setBiotypeIndex( "H1", Biotype::MethaneH().getIndex() );
        setBiotypeIndex( "H2", Biotype::MethaneH().getIndex() );
        setBiotypeIndex( "H3", Biotype::MethaneH().getIndex() );
        setBiotypeIndex( "H4", Biotype::MethaneH().getIndex() );

        setCompoundName("Methane");


    }
};

/**
 * \brief The small hydrocarbon ethane, C<sub>2</sub>H<sub>6</sub>, which has a single torsion degree of freedom.
 */
class  Ethane : public Molecule {
public:
    Ethane() 
    {
        setPdbResidueName("EHN");

        setBaseCompound("methyl1", MethylGroup());
        nameAtom("C1", "methyl1/C", Biotype::EthaneC().getIndex() );
        nameAtom("H1", "methyl1/H1", Biotype::EthaneH().getIndex() );
        nameAtom("H2", "methyl1/H2", Biotype::EthaneH().getIndex() );
        nameAtom("H3", "methyl1/H3", Biotype::EthaneH().getIndex() );

        // This first methyl is a base, not a decoration
        convertInboardBondCenterToOutboard();

        bondCompound("methyl2", MethylGroup(), "methyl1/bond");

        nameAtom("C2", "methyl2/C", Biotype::EthaneC().getIndex() );
        nameAtom("H4", "methyl2/H1", Biotype::EthaneH().getIndex() );
        nameAtom("H5", "methyl2/H2", Biotype::EthaneH().getIndex() );
        nameAtom("H6", "methyl2/H3", Biotype::EthaneH().getIndex() );

        defineDihedralAngle("torsion", "H1", "C1", "C2", "H4");

        setDefaultTorsionAngle(180*Deg2Rad);
        // setBondRotatable(false, "C1", "C2");

        setCompoundName("Ethane");

    }

    /**
     * \brief Sets dihedral angle about C-C bond 
     *
     * \return a reference to this Ethane object
     */
    Ethane& setDefaultTorsionAngle(Angle angle ///< dihedral angle about C-C bond in radians
        ) {
        setDefaultDihedralAngle("torsion", angle);

        return *this;
    }

    /**
     * \return the default (initial) dihedral angle about the C-C bond
     */
    Angle calcDefaultTorsionAngle() const {
        return calcDefaultDihedralAngle("torsion");
    }
};

/**
 * Base class for individual residue building blocks that comprise a Biopolymer chain
 *
 * Derives from Compound.
 * (though that may not be obvious in the automatically generated API documentation).
 */
class BiopolymerResidueRep;
// class SimTK_MOLMODEL_EXPORT BiopolymerResidue : public Compound {
class SimTK_MOLMODEL_EXPORT BiopolymerResidue : public Compound {
public:   
    /**
     * \brief Constructor for BiopolymerResidue
     */
    BiopolymerResidue(
        const Compound::Name& residueTypeName, ///< name for the type of BiopolymerResidue, e.g. "glycine"
        const String& threeLetterCode, ///< three letter code for the type of residue, e.g. "GLY".  This name will be used in the ResidueType field when PDB files are created.
        char oneLetterCode ///< one letter code for the type of residue.  e.g 'G'.  Use 'X' if the residue type is non-canonical.
        );

    /**
     * \return a reference to this BiopolymerResidue
     */
    BiopolymerResidue& setOneLetterCode(char olc ///< The one-letter-code for this type of residue.  Currently not used for anything.
        );

    /**
     * \return a reference to this BiopolymerResidue
     */
    BiopolymerResidue& setThreeLetterCode(const String& tlc ///< three letter code for the type of residue, e.g. "GLY".  This name will be used in the ResidueType field when PDB files are created.
        );

    /**
     * \return a reference to this BiopolymerResidue
     */
    BiopolymerResidue& setResidueTypeName(const String& name ///< name for the type of BiopolymerResidue, e.g. "glycine"
        );
    
    char getOneLetterCode() const;
    const String& getThreeLetterCode() const;
    const String& getResidueTypeName() const;

    /**
     * \brief Attempt to automatically assign a Biotype to each atom.
     *
     * Biotype assignments for each atom in a compound are required to match atoms to force field parameters.
     * Residue name and atom name, including synonyms, are used to identify Biotypes.  Previous loading of relevant Biotypes from 
     * a force field definition may be required for the assignBiotypes() method to work.
     *
     * \return true if Biotypes were successfully assigned
     */
    bool assignBiotypes(
        Ordinality::Residue ordinality = Ordinality::Any ///< whether this residue is at the beginning, middle, or end of the Biopolymer chain.
        );

    Compound::AtomIndex getParentCompoundAtomIndex(AtomIndex residueAtomIndex) const;

    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(BiopolymerResidue,BiopolymerResidueRep,Compound);
};

class BiopolymerRep;

// This class to point to residue atoms in a Biopolymer
class ResidueInfo {
public:
    /**
     * BiopolymerResidue::Index type is an integer index into residues of a Biopolymer.
     */
    SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(ResidueInfo,Index);
    SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(ResidueInfo,AtomIndex);

    class AtomInfo {
    public:
        friend class ResidueInfo;

        AtomInfo(Compound::AtomIndex index, const Compound::AtomName& name) 
            : pdbAtomName(name), biopolymerAtomIndex(index)
        {
            synonyms.insert(name);
            //std::cout<<__FILE__<<":"<<__LINE__<<" "<<name<<std::endl;
        }

        const std::set<Compound::AtomName>& getNames() const {
            return synonyms;
        }

    private:
        Compound::AtomIndex biopolymerAtomIndex;
        String pdbAtomName;
        std::set<Compound::AtomName> synonyms;
    };


    ResidueInfo(
        ResidueInfo::Index ix, 
        const Compound::Name& name, 
        const BiopolymerResidue& res,
        Compound::AtomIndex atomOffset,
        char insertionCode = ' ');

    AtomIndex addAtom(Compound::AtomIndex index, const Compound::AtomName& name) {
        ResidueInfo::AtomIndex answer(atoms.size());
        atoms.push_back(AtomInfo(index, name));
        atomIdsByName[name] = answer;
        return answer;
    }
    char getOneLetterCode() const {
        return oneLetterCode;
    }

    char setOneLetterCode(char myCode) {
        oneLetterCode = myCode	 ;
    }
    size_t getNumAtoms() const {
        return atoms.size();
    }

    const AtomInfo& getAtomInfo(AtomIndex a) const {
        return atoms[a];
    }
    AtomInfo& updAtomInfo(AtomIndex a) {
        return atoms[a];
    }

    const std::set<Compound::Name>& getNames() const {return synonyms;}

    const Compound::Name& getPdbResidueName() const {
        return pdbResidueName;
    }

    const int getPdbResidueNumber() const {
        return pdbResidueNumber;
    }

    const char getPdbInsertionCode() const {
        return pdbInsertionCode;
    }

    Compound::AtomIndex getAtomIndex(ResidueInfo::AtomIndex a) const {
        return atoms[a].biopolymerAtomIndex;
    }

    Compound::AtomIndex getAtomIndex(Compound::AtomName a) const {
        return atoms[atomIdsByName.find(a)->second].biopolymerAtomIndex;
    }

    const String& getAtomName(ResidueInfo::AtomIndex a) const {
        return atoms[a].pdbAtomName;
    }

    const std::set<Compound::AtomName>& getAtomSynonyms(ResidueInfo::AtomIndex a) const {
        return atoms[a].synonyms;
    }

    const Compound::Name& getName() const {return nameInCompound;}

    ResidueInfo& setPdbResidueNumber(int num) {
        pdbResidueNumber = num;
        return *this;
    }


    ResidueInfo& setPdbInsertionCode(char insertionCode) {
        pdbInsertionCode = insertionCode;
        return *this;
    }


    ResidueInfo::Index getIndex() const {return index;}

private:
    ResidueInfo::Index index;
    Compound::Name nameInCompound;
    Compound::Name pdbResidueName;
    char   oneLetterCode;
    std::set<Compound::Name> synonyms;
    int pdbResidueNumber;
    char pdbInsertionCode;
    std::vector<AtomInfo> atoms;
    std::map<const Compound::Name, AtomIndex> atomIdsByName;
    // TODO - put private data in ResidueInfoRep
};

/**
 * \brief The base class for DNA, RNA, and Protein molecules.
 *
 * Contains an ordered list of BiopolymerResidue subcompounds.
 * Derives from Molecule.
 * (though that may not be obvious in the automatically
 * generated API documentation).
 */
class SimTK_MOLMODEL_EXPORT Biopolymer : public Molecule
{
public:
    /**
     * String representing the sequence of a Protein or nucleic acid (DNA or RNA).  For use in Biopolymer constructor.
     */
    typedef String Sequence;
    
    /**
     * \brief Default constructore for Biopolymer.  Produces a Biopolymer with no atoms nor residues.
     */
    Biopolymer();

    /**
     * \return The number of residues in the polymer chain.
     */
    int getNumResidues() const;

    /**
     * \warning The residueIndex is ordinarily NOT the same as the PdbResidueNumber of the BiopolymerResidue.
     *
     * \return A read-only reference to a BiopolymerResidue subcompound of this Biopolymer molecule
     */
    const ResidueInfo& getResidue(
        ResidueInfo::Index residueIndex ///< integer index of residue in context of this Biopolymer, in the range zero (0) to (getNumResidues() - 1).
        ) const;

    /**
     * \warning The residueIndex is ordinarily NOT the same as the PdbResidueNumber of the BiopolymerResidue.
     *
     * \return A mutable reference to a BiopolymerResidue subcompound of this Biopolymer molecule
     */
    ResidueInfo& updResidue(
        ResidueInfo::Index residueIndex ///< integer index of residue in context of this Biopolymer, in the range zero (0) to (getNumResidues() - 1).
        ); 

    /**
     * \return A read-only reference to a BiopolymerResidue subcompound of this Biopolymer molecule
     */
    const ResidueInfo& getResidue(
        Compound::Name residueName ///< a Compound::Name for the BiopolymerResidue from the viewpoint of this Biopolymer.
        ) const;

    /**
     * \return A mutable reference to a BiopolymerResidue subcompound of this Biopolymer molecule
     */
    ResidueInfo& updResidue(
        Compound::Name residueName ///< a Compound::Name for the BiopolymerResidue from the viewpoint of this Biopolymer.
        );

    /**
     * \warning The residueIndex is ordinarily NOT the same as the PdbResidueNumber of the BiopolymerResidue.
     *
     * \return A Compound::Name for a BiopolymerResidue from the viewpoint of this containing Biopolymer.
     */
    const Compound::Name& getResidueName(
        ResidueInfo::Index residueIndex ///< integer index of residue in context of this Biopolymer, in the range zero (0) to (getNumResidues() - 1).
        ) const;

    // BiopolymerResidue createResidueCompound(ResidueInfo::Index r) const;

    /**
     * \brief This method renumbers the PDB Residue Numbers of the biopolymer chain to start at the supplied integer and increase consecutively.
     *
     * \return A reference to this Biopolymer.
     *
     */ 
    Biopolymer& renumberPdbResidues(int firstPdbResidueNumber); 
    
    bool assignResidueBiotypes(ResidueInfo::Index, Ordinality::Residue);

    /**
     * \brief Attempt to automatically assign a Biotype to each atom.
     *
     * Biotype assignments for each atom in a compound are required to match atoms to force field parameters.
     * Residue name and atom name, including synonyms, are used to identify Biotypes.  Previous loading of relevant Biotypes from 
     * a force field definition may be required for the assignBiotypes() method to work.
     */
    void assignBiotypes();

    /**
     * \brief Attach a new residue onto the end of the current Biopolymer chain
     */
    ResidueInfo::Index appendResidue(
        const Compound::Name& resName, ///< new name for the new residue, local to this parent Biopolymer.
        const BiopolymerResidue& residue ///< template residue to copy onto the end of the chain.  The residue will be copied, not incorporated.
        );

    /**
     * \brief Attach a new residue onto the end of the current Biopolymer chain
     */
    ResidueInfo::Index appendResidue(
        const Compound::Name& resName, ///< new residue name, from the viewpoint of this parent Biopolymer
        const BiopolymerResidue& residue, ///< template residue to copy onto the end of the chain.  The residue will be copied, not incorporated.
        BondMobility::Mobility mobility ///< allowed motion of the new bond connecting the new residue to the rest of the chain
        );

    Biopolymer& setResidueBondMobility(ResidueInfo::Index, BondMobility::Mobility);
    MobilizedBodyIndex getResidueAtomMobilizedBodyIndex(ResidueInfo::Index res, ResidueInfo::AtomIndex a) const {
        return getAtomMobilizedBodyIndex(getResidue(res).getAtomIndex(a));
    }
    Vec3 getResidueAtomLocationInMobilizedBodyFrame(ResidueInfo::Index res, ResidueInfo::AtomIndex a) const {
        return getAtomLocationInMobilizedBodyFrame(getResidue(res).getAtomIndex(a));
    }

    virtual AtomTargetLocations createAtomTargets(const class PdbStructure& targetStructure,const bool guessCoordinates  = false) const;
    virtual AtomTargetLocations createAtomTargets(const class PdbChain& targetChain,const bool guessCoordinates  = false) const;

    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Biopolymer, BiopolymerRep, Molecule);

    // OBSOLETE; TODO: remove in SimTK 2.0
    int getNResidues() const {return getNumResidues();}

private:
};


} // namespace SimTK

#endif // SimTK_MOLMODEL_COMPOUND_H_
