/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#include "TetherForce.h"  
    TetherForce::TetherForce (
                SimbodyMatterSubsystem &    matter,
                GeneralForceSubsystem  &    forces,                
                const MobilizedBodyIndex    body1Index,
                const Vec3 &    station1,
                const MobilizedBodyIndex    body2Index,
                const Vec3 &    station2,
                Real    k,
                Real    x0
                ):matter(matter),forces(forces), body1Index(body1Index), station1(station1),  body2Index(body2Index), station2(station2),k(k),x0(x0) 

        { 
    };    
    void TetherForce::calcForce(const State& state, Vector_<SpatialVec>& bodyForces,  
            Vector_<Vec3>& particleForces, Vector& mobilityForces) const 
        {
          const MobilizedBody body1 = matter.getMobilizedBody(body1Index);
          const MobilizedBody body2 = matter.getMobilizedBody(body2Index);
          
          const Transform& X_GB1 = /*matter.getMobilizedBody*/(body1).getBodyTransform(state);
          const Transform& X_GB2 = /*matter.getMobilizedBody*/(body2).getBodyTransform(state);

          const Vec3 s1_G = X_GB1.R() * station1;
          const Vec3 s2_G = X_GB2.R() * station2;
	  //cout<<__FILE__<<":"<<__LINE__<<" : station 1, 2: "<<station1<<" , "<<station2<<endl;
	  //cout<<__FILE__<<":"<<__LINE__<<" : body index 1, 2: "<<body1Index<<" , "<<body2Index <<endl;
          const Vec3 p1_G = X_GB1.p() + s1_G; // station measured from ground origin
          const Vec3 p2_G = X_GB2.p() + s2_G;

          const Vec3 r_G               = p2_G - p1_G; // vector from point1 to point2
          const Real d                 = r_G.norm();  // distance between the points
          const Real stretch   = d - x0;              // + -> tension, - -> compression
          Real frcScalar =0 ;
          if (d <= x0) {  // correct is <= ..  for quique did >= .
              frcScalar = 0;
          } else {
              frcScalar = k*stretch;
          }
          const Vec3 f1_G = (frcScalar/d) * r_G;
          //cout<<"scf5 "<<d<<" , "<<stretch<<" , "<< f1_G<< endl;
          bodyForces[body1.getMobilizedBodyIndex()] +=  SpatialVec(s1_G % f1_G, f1_G);
          bodyForces[body2.getMobilizedBodyIndex()] -=  SpatialVec(s2_G % f1_G, f1_G);
        };
    Real TetherForce::calcPotentialEnergy(const State& state) const { 

        const MobilizedBody body1 = matter.getMobilizedBody(body1Index);
        const MobilizedBody body2 = matter.getMobilizedBody(body2Index);
        const Transform& X_GB1 = (body1).getBodyTransform(state);
        const Transform& X_GB2 = (body2).getBodyTransform(state);

        const Vec3 s1_G = X_GB1.R() * station1;
        const Vec3 s2_G = X_GB2.R() * station2;

        const Vec3 p1_G = X_GB1.p() + s1_G; // station measured from ground origin
        const Vec3 p2_G = X_GB2.p() + s2_G;

        const Vec3 r_G = p2_G - p1_G; // vector from point1 to point2
        const Real d   = r_G.norm();  // distance between the points
        const Real stretch   = d - x0; // + -> tension, - -> compression
        MMBLOG_FILE_FUNC_LINE(DEBUG , "tether extension : "<<d<<" maximum : "<<x0<<" difference: "<<stretch<<" nm"); // used to convert to Å, now using nm directly
        //cout<<__FILE__<<":"<<__LINE__<<" tether extension : "<<d<<" maximum : "<<x0<<" difference: "<<stretch<<" nm"<<endl; // used to convert to Å, now using nm directly
        if (d <= x0) {
            return 0;
        } else {
            return 0.5*k*stretch*stretch;
        }
    };
    bool TetherForce::dependsOnlyOnPositions() const  { 
        return true; 
    };    
