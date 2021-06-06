/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */
#include "NtCForces.h"  
#include <string.h>
#include <sstream>
#include <Utils.h>
#include "MMBLogger.h"
#include "ParameterReader.h"

    NTC_Torque::NTC_Torque (SimbodyMatterSubsystem& matter,ParameterReader& myParameterReader,  NTC_PAR_Class& myNTC_PAR_Class, BiopolymerClassContainer & myBiopolymerClassContainer, std::ostream& outputStream ) : matter(matter),myParameterReader(myParameterReader), myNTC_PAR_Class (myNTC_PAR_Class), myBiopolymerClassContainer(myBiopolymerClassContainer), outputStream(outputStream)
        { 
    };    
    
    void NTC_Torque::calcForce(const State& state, Vector_<SpatialVec>& bodyForces,  
            Vector_<Vec3>& particleForces, Vector& mobilityForces) const 
        {  
        MobilizedBody body1;
        MobilizedBody body2;
        MobilizedBody body3;
        MobilizedBody body4;
        Transform transform1;
        Transform transform2;
        double torqueConstant;
        double pot_angle;
        double angle;
        double dih,bias;
        int    value,i;
        double prob[361];
        
        AtomSpringContainer atomSpringContainer;
        
        for (int r=0;r<myParameterReader.ntc_class_container.numNTC_Torsions();r++) 
        { 
        
            String chainId1=(myParameterReader.ntc_class_container.getNTC_Class(r)).NtC_FirstBPChain;            
            NTC_PAR_BondRow myNTC_PAR_BondRow = myNTC_PAR_Class.myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[(myParameterReader.ntc_class_container.getNTC_Class(r)).NTC_PAR_BondRowIndex];
            String basePairIsTwoTransformForce="ntcstep";
            ResidueID residueNumber1=(myParameterReader.ntc_class_container.getNTC_Class(r)).FirstBPResidue;
            ResidueID residueNumber2=(myParameterReader.ntc_class_container.getNTC_Class(r)).SecondBPResidue;
            ResidueID myResidueNumber,myResidueNumber2;

	    if ( myNTC_PAR_BondRow.bondLength[0] == 0.0 ) {
//		cout << "torsion " << r << " is real torsion" << endl;
            
         //   myResidueNumber = residueNumber1;            
         //   myResidueNumber.ResidueNumber += stoi(myNTC_PAR_BondRow.atom_shift[0]);             
             Vec3  state_1,state_2,state_3,state_4;
        
            if(stoi(myNTC_PAR_BondRow.atom_shift[0]) == 0) myResidueNumber = residueNumber1;
            if(stoi(myNTC_PAR_BondRow.atom_shift[0]) == 1) myResidueNumber = residueNumber2;
            
            body1 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[0]);
            state_1 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[0]); 
            
            if(stoi(myNTC_PAR_BondRow.atom_shift[1]) == 0) myResidueNumber = residueNumber1;
            if(stoi(myNTC_PAR_BondRow.atom_shift[1]) == 1) myResidueNumber = residueNumber2;
            
            body2 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[1]);
            state_2 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[1]); 
            
            if(stoi(myNTC_PAR_BondRow.atom_shift[2]) == 0) myResidueNumber = residueNumber1;
            if(stoi(myNTC_PAR_BondRow.atom_shift[2]) == 1) myResidueNumber = residueNumber2;
            
            body3 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[2]);
            state_3 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[2]); 

            if(stoi(myNTC_PAR_BondRow.atom_shift[3]) == 0) myResidueNumber = residueNumber1;
            if(stoi(myNTC_PAR_BondRow.atom_shift[3]) == 1) myResidueNumber = residueNumber2;

            body4 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[3]);
            state_4 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[3]);               
            
            torqueConstant = myNTC_PAR_BondRow.torqueConstant;            
            
            Vec3 d_d1,d_d2,d_d3;
            
            d_d1 = state_2 - state_1;
            d_d2 = state_3 - state_2;
            d_d3 = state_4 - state_3;
            
            Vec3 cross_1, cross_2;
            
            cross_1 = d_d1 % d_d2;
            cross_2 = d_d2 % d_d3;
            
            cross_1 = cross_1 / cross_1.norm();
            cross_2 = cross_2 / cross_2.norm();
        
            Vec3 cross_3;
            
            cross_3 = cross_1 % cross_2;
            
            
     angle = return_angle(cross_1,cross_2,cross_3,d_d2);
     
     double dist_ang = return_dist_ang(angle,myNTC_PAR_BondRow.rotationAngle);

    // pot_angle = torqueConstant*(dist_ang/57.295779513);//(exp(-(pow(dist_ang,2)/(2.0*pow(l_param,2)))))*(dist_ang/57.295779513)/(2.0*pow(l_param/57.295779513,2))*(sin(dist_ang/57.295779513));       
     
       if((myParameterReader.ntc_class_container.getNTC_Class(r)).meta == 0) pot_angle = torqueConstant*(myParameterReader.ntc_class_container.getNTC_Class(r)).weight*(-sin((dist_ang + 180.0)/57.295779513))*(360.0/57.295779513 + 1.0)/(1.0 + myNTC_PAR_BondRow.CONFALVALUE)/(360.0/57.295779513);
     
     Vec3 torque;    
     
       if((myParameterReader.ntc_class_container.getNTC_Class(r)).meta == 0) 
          
       { 
           
      //   cout << pot_angle << "  -- pot_angle " << endl;  
           
         torque = d_d2/d_d2.norm()*pot_angle;
   
         bodyForces[body1.getMobilizedBodyIndex()] += SpatialVec(torque,Vec3(0));
         bodyForces[body4.getMobilizedBodyIndex()] -= SpatialVec(torque,Vec3(0));
         bodyForces[body2.getMobilizedBodyIndex()] += SpatialVec(torque,Vec3(0));
         bodyForces[body3.getMobilizedBodyIndex()] -= SpatialVec(torque,Vec3(0));  
       
        };
       
       if((myParameterReader.ntc_class_container.getNTC_Class(r)).meta == 1) {
           
        dih = 0;
        value = -1;
        bias = 0;
           
        angle *= 57.295779513;
         
       if(isfinite(angle) == 1) { 
        
        i = (int) round(angle);
        
        myBiopolymerClassContainer.hist[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][i] += 1.0;
        value = i;
                  
        myBiopolymerClassContainer.counter[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += 1.0;
             
        if(myBiopolymerClassContainer.counter[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] > 0.0) myBiopolymerClassContainer.prob[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][i] = myBiopolymerClassContainer.hist[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][i]/myBiopolymerClassContainer.counter[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]; 
               
        if(value > 0 && myBiopolymerClassContainer.prob[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][value] > 1e-3 && myBiopolymerClassContainer.prob[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][value+1] > 1e-3 && value < 360) bias = 2.479*(log(myBiopolymerClassContainer.prob[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][value]/(1e-3))-log(myBiopolymerClassContainer.prob[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][value+1]/(1e-3)));
        if(prob[value] > 1e-3 && myBiopolymerClassContainer.prob[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][value+1] > 1e-3 && value == 360) bias = 2.479*(log(myBiopolymerClassContainer.prob[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][value]/(1e-3))-log(myBiopolymerClassContainer.prob[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][1]/(1e-3))); 
           
       /* if( isfinite(log(prob[value]/(1e-3))) == 1 && isfinite(log(prob[value+1]/(1e-3))) == 1) {   
           
         cout << bias << " value " << value << " bias "<< endl;
         cout << pot_angle << " pot angle " << endl;  
         cout << log(prob[value]/(1e-3)) << log(prob[value+1]/(1e-3)) << " log prob " << endl; 
         
        };*/
         
         pot_angle = torqueConstant*(-sin((dist_ang + 180.0)/57.295779513))*(360.0/57.295779513 + 1.0)/(1.0 + myNTC_PAR_BondRow.CONFALVALUE)/(360.0/57.295779513);
       
         torque = d_d2/d_d2.norm()*(pot_angle)/(1.0+(myParameterReader.ntc_class_container.getNTC_Class(r)).weight2)*(myParameterReader.ntc_class_container.getNTC_Class(r)).weight;
   
         bodyForces[body1.getMobilizedBodyIndex()] += SpatialVec(torque,Vec3(0));
         bodyForces[body4.getMobilizedBodyIndex()] -= SpatialVec(torque,Vec3(0));
         bodyForces[body2.getMobilizedBodyIndex()] += SpatialVec(torque,Vec3(0));
         bodyForces[body3.getMobilizedBodyIndex()] -= SpatialVec(torque,Vec3(0));          
         
        // cout << (pot_angle)/(1.0+(myParameterReader.ntc_class_container.getNTC_Class(r)).weight) << " pot angle " << endl;
        // cout << (myParameterReader.ntc_class_container.getNTC_Class(r)).weight << " weight " << endl;
        // cout << pot_angle << " pot angle " << endl;
         
       if( isfinite(log(myBiopolymerClassContainer.prob[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][value]/(1e-3))) == 1 && isfinite(log(myBiopolymerClassContainer.prob[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][value+1]/(1e-3))) == 1) {       
         
       if( isfinite(bias) == 1 && sqrt(pow(bias,2)) > 0.0) { 
           
        torque = -d_d2/d_d2.norm()*(bias)*sqrt(pow(pot_angle,2))/sqrt(pow(bias,2))*(myParameterReader.ntc_class_container.getNTC_Class(r)).weight2*(myParameterReader.ntc_class_container.getNTC_Class(r)).weight;
         
         bodyForces[body1.getMobilizedBodyIndex()] += SpatialVec(torque,Vec3(0));
         bodyForces[body4.getMobilizedBodyIndex()] -= SpatialVec(torque,Vec3(0));
         bodyForces[body2.getMobilizedBodyIndex()] += SpatialVec(torque,Vec3(0));
         bodyForces[body3.getMobilizedBodyIndex()] -= SpatialVec(torque,Vec3(0));  
  
        };
       };
     };
         
   }; 
       
// end real torsions
	}
// bonds
	else {

            Vec3  state_1;
            Vec3  state_2;
        
            body1 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,residueNumber1,myNTC_PAR_BondRow.residue1Atom[0]);
            state_1 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,residueNumber1,myNTC_PAR_BondRow.residue1Atom[0]); 
            
            body2 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,residueNumber2,myNTC_PAR_BondRow.residue1Atom[1]);
            state_2 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,residueNumber2,myNTC_PAR_BondRow.residue1Atom[1]); 
            
            Vec3 ptp = state_2 - state_1;
            double d = ptp.norm(); //sqrt(pow(x_d2-x_d1,2) + pow(y_d2-y_d1,2) + pow(z_d2-z_d1,2));
            double frc;
            Vec3   frcVec;
            
    if((myParameterReader.ntc_class_container.getNTC_Class(r)).meta == 0) {
      
      frc = (1.0-exp(-(2.0*myNTC_PAR_BondRow.CONFALVALUE)*(d-myNTC_PAR_BondRow.bondLength[0])))*(-exp(-(2.0*myNTC_PAR_BondRow.CONFALVALUE)*(d-myNTC_PAR_BondRow.bondLength[0])))*myNTC_PAR_BondRow.springConstant[0]*(myParameterReader.ntc_class_container.getNTC_Class(r)).weight;
      frcVec = (frc)*ptp/d;
      
       bodyForces[body1.getMobilizedBodyIndex()] -=  SpatialVec(frcVec,Vec3(1));     
	   bodyForces[body2.getMobilizedBodyIndex()] +=  SpatialVec(frcVec,Vec3(1));       
      
    };
  
    
  if((myParameterReader.ntc_class_container.getNTC_Class(r)).meta == 1) {  
  
      
      bias = 0.0;
      
       if(d < 3.0) {
           
       i = (int)round((d)*10.0);
       
       myBiopolymerClassContainer.hist_d[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][i] += 1.0;
       value = i;
                   
           
       myBiopolymerClassContainer.counter_d[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += 1.0;
           
       if(myBiopolymerClassContainer.counter_d[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] > 0.0) myBiopolymerClassContainer.prob_d[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][i] = myBiopolymerClassContainer.hist_d[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][i]/myBiopolymerClassContainer.counter_d[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]; 
           
        if(value > 0 && myBiopolymerClassContainer.prob_d[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][value] > 1e-3 && myBiopolymerClassContainer.prob_d[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][value+1] > 1e-3 && value < 31) bias = 2.479*(log(myBiopolymerClassContainer.prob_d[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][value]/(1e-3))-log(myBiopolymerClassContainer.prob_d[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][value+1]/(1e-3)));

        frc =  (1.0-exp(-(2.0*myNTC_PAR_BondRow.CONFALVALUE)*(d-myNTC_PAR_BondRow.bondLength[0])))*(-exp(-(2.0*myNTC_PAR_BondRow.CONFALVALUE)*(d-myNTC_PAR_BondRow.bondLength[0])))*myNTC_PAR_BondRow.springConstant[0]; 
        frc =  frc/(1.0 + (myParameterReader.ntc_class_container.getNTC_Class(r)).weight2);
        
        frcVec = (frc)*ptp/d*(myParameterReader.ntc_class_container.getNTC_Class(r)).weight;
        
	    bodyForces[body1.getMobilizedBodyIndex()] -=  SpatialVec(frcVec,Vec3(1));     
	    bodyForces[body2.getMobilizedBodyIndex()] +=  SpatialVec(frcVec,Vec3(1));        
       
     if( isfinite(bias) == 1 && sqrt(pow(bias,2)) > 0.0) {   
        
        bias = bias*(sqrt(pow(frc,2))/(sqrt(pow(bias,2))))*(myParameterReader.ntc_class_container.getNTC_Class(r)).weight2; 
        frcVec = (bias)*ptp/d*(myParameterReader.ntc_class_container.getNTC_Class(r)).weight;
        
	    bodyForces[body1.getMobilizedBodyIndex()] +=  SpatialVec(frcVec,Vec3(1));     
	    bodyForces[body2.getMobilizedBodyIndex()] -=  SpatialVec(frcVec,Vec3(1));          
        
      };
     };
     
   }; 
  };
    
  };
    
  };
        
        
    Real NTC_Torque::calcPotentialEnergy(const State& state) const { 

        double energy = 0.0;  
        double rms    = 0.0;
        double rmsTorsionAngleForThisNtCAndDinucleotide    = 0.0;
	String oldNtCClassString = "ZZZZZZ";
        MobilizedBody body1;
        MobilizedBody body2;
        MobilizedBody body3;
        MobilizedBody body4;
        Transform transform1;
        Transform transform2;
        double torqueConstant;
        double angle;

        
        for (int r=0;r<myParameterReader.ntc_class_container.numNTC_Torsions();r++) 
        { 
           // If we have changed our NtC class type, meaning we are computing a new NtC Class, not just an additional torsion in the same class:		
	   if (!((myParameterReader.ntc_class_container.getNTC_Class(r)).NtC_Class_String == oldNtCClassString)) {
               if (r > 0) {
                   MMBLOG_FILE_FUNC_LINE(
                       INFO,
		       "RMSD Angle sum for NtC of type "<<(myParameterReader.ntc_class_container.getNTC_Class(r)).NtC_step_ID <<" NtC Index: "<<  (myParameterReader.ntc_class_container.getNTC_Class(r)).NtC_INDEX <<" "  <<   oldNtCClassString <<" is "  <<rmsTorsionAngleForThisNtCAndDinucleotide << endl
                   );
               }

	       rmsTorsionAngleForThisNtCAndDinucleotide = 0.; 

	   }
            //rmsTorsionAngleForThisNtCAndDinucleotide = 0;		
            MMBLOG_FILE_FUNC_LINE(DEBUG, "doing base pair #"<<r<<endl);
        
            String chainId1=(myParameterReader.ntc_class_container.getNTC_Class(r)).NtC_FirstBPChain;            
            NTC_PAR_BondRow myNTC_PAR_BondRow = myNTC_PAR_Class.myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[(myParameterReader.ntc_class_container.getNTC_Class(r)).NTC_PAR_BondRowIndex];
            String basePairIsTwoTransformForce="ntcstep";
            ResidueID residueNumber1=(myParameterReader.ntc_class_container.getNTC_Class(r)).FirstBPResidue;
            ResidueID residueNumber2=(myParameterReader.ntc_class_container.getNTC_Class(r)).SecondBPResidue;
            ResidueID myResidueNumber,myResidueNumber2;

	    if ( myNTC_PAR_BondRow.bondLength[0] == 0.0 ) {
//		cout << "torsion " << r << " is real torsion" << endl;
            
             Vec3  state_1,state_2,state_3,state_4;
        
            if(stoi(myNTC_PAR_BondRow.atom_shift[0]) == 0) myResidueNumber = residueNumber1;
            if(stoi(myNTC_PAR_BondRow.atom_shift[0]) == 1) myResidueNumber = residueNumber2;
            
            body1 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[0]);
            state_1 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[0]); 
            
            if(stoi(myNTC_PAR_BondRow.atom_shift[1]) == 0) myResidueNumber = residueNumber1;
            if(stoi(myNTC_PAR_BondRow.atom_shift[1]) == 1) myResidueNumber = residueNumber2;
            
            body2 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[1]);
            state_2 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[1]); 
            
            if(stoi(myNTC_PAR_BondRow.atom_shift[2]) == 0) myResidueNumber = residueNumber1;
            if(stoi(myNTC_PAR_BondRow.atom_shift[2]) == 1) myResidueNumber = residueNumber2;
            
            body3 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[2]);
            state_3 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[2]); 

            if(stoi(myNTC_PAR_BondRow.atom_shift[3]) == 0) myResidueNumber = residueNumber1;
            if(stoi(myNTC_PAR_BondRow.atom_shift[3]) == 1) myResidueNumber = residueNumber2;

            body4 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[3]);
            state_4 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[3]);               
            
            torqueConstant = myNTC_PAR_BondRow.torqueConstant;            

            Vec3 d_d1,d_d2,d_d3;
            
            d_d1 = state_2 - state_1;
            d_d2 = state_3 - state_2;
            d_d3 = state_4 - state_3;
            
            Vec3 cross_1, cross_2;
            
            cross_1 = d_d1 % d_d2;
            cross_2 = d_d2 % d_d3;
            
            cross_1 = cross_1 / cross_1.norm();
            cross_2 = cross_2 / cross_2.norm();
        
            Vec3 cross_3;
            
            cross_3 = cross_1 % cross_2;
        
     angle = return_angle(cross_1,cross_2,cross_3,d_d2);
     
     double dist_ang = return_dist_ang(angle,myNTC_PAR_BondRow.rotationAngle);  

    energy += torqueConstant*cos((dist_ang + 180.0)/57.295779513)*(360.0/57.295779513+1.0)/(1.0+myNTC_PAR_BondRow.CONFALVALUE)/(360.0/57.295779513);//torqueConstant*pow(dist_ang/57.295779513,2);//-torqueConstant*(-cos(dist_ang/57.295779513)+(exp(-(pow(dist_ang,2)/(2.0*pow(l_param,2))))));
    MMBLOG_FILE_FUNC_LINE(DEBUG, " NTC sampling - CHAIN ID = " << chainId1 << ", residuenumber " << myResidueNumber.ResidueNumber  << " difference-angle = "<< dist_ang << " , CONFALVALUE = " << myNTC_PAR_BondRow.CONFALVALUE << " , " << angle*57.295779513 << " = angle at time t for atoms  = " << myNTC_PAR_BondRow.residue1Atom[0] << " , " << myNTC_PAR_BondRow.residue1Atom[1] << " , " << myNTC_PAR_BondRow.residue1Atom[2] << " , " << myNTC_PAR_BondRow.residue1Atom[3] << " , "<< myNTC_PAR_BondRow.rotationAngle*57.295779513 << " = angle_0 from  input , " << "energy = " << energy << endl);    
    //cout << " NTC sampling - CHAIN ID = " << chainId1 << ", residuenumber " << myResidueNumber.ResidueNumber  << " difference-angle = "<< dist_ang << " , CONFALVALUE = " << myNTC_PAR_BondRow.CONFALVALUE << " , " << angle*57.295779513 << " = angle at time t for atoms  = " << myNTC_PAR_BondRow.residue1Atom[0] << " , " << myNTC_PAR_BondRow.residue1Atom[1] << " , " << myNTC_PAR_BondRow.residue1Atom[2] << " , " << myNTC_PAR_BondRow.residue1Atom[3] << " , "<< myNTC_PAR_BondRow.rotationAngle*57.295779513 << " = angle_0 from  input , " << "energy = " << energy << endl;
        
    rms   += sqrt(pow(dist_ang,2));    
    rmsTorsionAngleForThisNtCAndDinucleotide += sqrt(pow(dist_ang,2));
// end real torsions
	}
// bonds
	else {
	
	    Vec3  state_1;
	    Vec3  state_2;
                    
            body1 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,residueNumber1,myNTC_PAR_BondRow.residue1Atom[0]);
            state_1 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,residueNumber1,myNTC_PAR_BondRow.residue1Atom[0]); 
            
            
            body2 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,residueNumber2,myNTC_PAR_BondRow.residue1Atom[1]);
            state_2 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,residueNumber2,myNTC_PAR_BondRow.residue1Atom[1]); 
            
            Vec3 diff = state_2 - state_1;
            
            double d = diff.norm();

              MMBLOG_FILE_FUNC_LINE(DEBUG, "NN|CC difference: " << (d - myNTC_PAR_BondRow.bondLength[0]) << ", current value: " << d << ", equilibrium value: " << myNTC_PAR_BondRow.bondLength[0] << endl);
	      //cout << "NN|CC difference: " << (d - myNTC_PAR_BondRow.bondLength[0]) << ", current value: " << d << ", equilibrium value: " << myNTC_PAR_BondRow.bondLength[0] << endl;

          energy += myNTC_PAR_BondRow.springConstant[0]*pow(1.0-exp(-(2.0*myNTC_PAR_BondRow.CONFALVALUE)*(d - myNTC_PAR_BondRow.bondLength[0])),2);
          
	     };
           //cout << __LINE__<< "  RMSD Angle sum for NtC of type "<<(myParameterReader.ntc_class_container.getNTC_Class(r)).NtC_step_ID <<" NtC Index: "<<  (myParameterReader.ntc_class_container.getNTC_Class(r)).NtC_INDEX <<" "  <<   (myParameterReader.ntc_class_container.getNTC_Class(r)).NtC_Class_String<<" is "  <<rmsTorsionAngleForThisNtCAndDinucleotide << endl;
	   // If this is the last torsion:
	   if (r == (myParameterReader.ntc_class_container.numNTC_Torsions()-1)){
               MMBLOG_FILE_FUNC_LINE(INFO, "RMSD Angle sum for NtC of type "<<(myParameterReader.ntc_class_container.getNTC_Class(r)).NtC_step_ID <<" NtC Index: "<<  (myParameterReader.ntc_class_container.getNTC_Class(r)).NtC_INDEX <<" "  <<   (myParameterReader.ntc_class_container.getNTC_Class(r)).NtC_Class_String<<" is "  <<rmsTorsionAngleForThisNtCAndDinucleotide << endl);
           }
	   oldNtCClassString = (myParameterReader.ntc_class_container.getNTC_Class(r)).NtC_Class_String;
       }

       MMBLOG_PLAIN_NOSEV(INFO, rms << " = RMSD Angle sum " << endl);
       
       return energy;   
          

         
    };
    bool NTC_Torque::dependsOnlyOnPositions() const  { 
        return true; 
    };  
    
    
Real NTC_Torque::return_dist_ang(double angle,double rotationAngle) const
{
       
     double ang_diff = (angle - rotationAngle)*57.295779513; // Deg
     double dist_ang = 180.0 - abs(180.0 - abs(ang_diff));
     int angle_1 = int(round(angle*57.295779513));
     int angle_2 = int(round(rotationAngle*57.295779513));

     int interval_begin = angle_2;
     int interval_end   = (interval_begin + 180) % 360;

     if (interval_end > interval_begin) {
         if (angle_1 < interval_begin || angle_1 > interval_end ) {
             dist_ang = -dist_ang;
        }
     } else {
         if (angle_1 < interval_begin && angle_1 > interval_end ) {
             dist_ang = -dist_ang;
        };
    };    
    
  return dist_ang;  
} 

Real NTC_Torque::return_angle(Vec3 cross_1,Vec3 cross_2,Vec3 cross_3,Vec3 d_d2) const 
{
   
   double angle; 
   double PI = 3.14159265359;
           
   Vec3 direction;
            
   direction[0] = cross_3[0]*d_d2[0];
   direction[1] = cross_3[1]*d_d2[1];  
   direction[2] = cross_3[2]*d_d2[2];
            
   double scalar_product = dot(cross_1,cross_2);
        
   if(scalar_product > 1.0) scalar_product = 1.0;
   if(scalar_product < -1.0) scalar_product = -1.0;       
      
   angle = acos(scalar_product)*180.0/PI;
        
   if(direction[0] < 0.0 && direction[1] < 0.0 && direction[2] < 0.0) {
            
      angle = -angle;
            
   };          
            
   if(angle < 0.0) angle = angle + 360.0;
        
   angle /= 57.295779513;    
    
  return angle;  
}
    
