#include "SimTKsimbody.h" 
#include "/usr/local/SimTK/core/include/simbody/internal/ForceSubsystem.h" 
 
using namespace SimTK; 



class ExampleSubsystemImpl : public ForceSubsystem::Guts { 
public: 
    ExampleSubsystemImpl() : ForceSubsystem::Guts("Example", "1.0") { 
    } 
    Subsystem::Guts* cloneImpl() const { 
        return new ExampleSubsystemImpl(); 
    } 
    int realizeSubsystemDynamicsImpl(const State& state) const { 
        const MultibodySystem& system =  
                static_cast<const MultibodySystem&>(getSystem()); 
        const SimbodyMatterSubsystem& matter = system.getMatterSubsystem(); 
        Vector_<SpatialVec>& forces = system.updRigidBodyForces(state,  
                Stage::Dynamics); 
        for (MobilizedBodyIndex i(0); i < matter.getNBodies(); i++) { 
            const MobilizedBody& body1 = matter.getMobilizedBody(i); 
            for (MobilizedBodyIndex j(0); j < i; j++) { 
                const MobilizedBody& body2 = matter.getMobilizedBody(j); 
                Vec3 r = body1.getBodyOriginLocation(state)- 
                         body2.getBodyOriginLocation(state); 
                Real distance = r.norm(); 
                Vec3 force = r/(distance*distance*distance); 
                forces[i][1] += force; 
                forces[j][1] -= force; 
            } 
        } 

        return 0; 
    } 
    Real calcPotentialEnergy(const State& state) const { 
        const MultibodySystem& system = static_cast<const  
                MultibodySystem&>(getSystem()); 
        const SimbodyMatterSubsystem& matter = system.getMatterSubsystem(); 
        double energy = 0.0; 
        for (MobilizedBodyIndex i(0); i < matter.getNBodies(); i++) { 
            const MobilizedBody& body1 = matter.getMobilizedBody(i); 
            for (MobilizedBodyIndex j(0); j < i; j++) { 
                const MobilizedBody& body2 = matter.getMobilizedBody(j); 
                Vec3 r = body1.getBodyOriginLocation(state)- 
                         body2.getBodyOriginLocation(state); 
                energy -= 1.0/r.norm(); 
            } 
        } 
        return energy; 
    } 
}; 
 
class ExampleSubsystem : public ForceSubsystem { 
public: 
    ExampleSubsystem(MultibodySystem& system) { 
        adoptSubsystemGuts(new ExampleSubsystemImpl()); 
        system.addForceSubsystem(*this); 
    } 
}; 

