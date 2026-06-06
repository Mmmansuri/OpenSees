/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.1 $
// $Date: 2020-03-09 22:27:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/SelfCenteringMaterial.cpp,v $

// Written: Mansuri
// Created: 2019, checked:2026
// Revision: A
//
// Description: This file contains the class implementation for 
// SelfCenteringMaterial. 

#include <SelfCenteringMaterial.h>
#include <Vector.h>
#include <Channel.h>
#include <Matrix.h>
#include <Information.h>
#include <Parameter.h>

#include <math.h>
#include <float.h>
#include <elementAPI.h>

void* OPS_SelfCenteringMaterial()
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 16) {
	opserr << "WARNING: Insufficient arguments\n";
	opserr << "Want: uniaxialMaterial SelfCentering tag? k1? k2? k3? SecondFuse? epsA? epsB? ForcB? a3? a4? a5? a6? a7? a8? ";
	opserr << "ActF? beta? <SlipDef? BearDef? rBear?>" << endln;
	return 0;
    }

    int tag;
    numdata = 1;
    if (OPS_GetIntInput(&numdata,&tag) < 0) {
	opserr << "WARNING invalid tag\n";
	return 0;
    }

    double data[18] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 18) {
	numdata = 18;
    }
    if (OPS_GetDoubleInput(&numdata,data)) {
	opserr << "WARNING invalid double inputs\n";
	return 0;
    }

    UniaxialMaterial* mat = new SelfCenteringMaterial(tag,data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8],data[9],data[10],data[11],data[12],data[13],data[14],data[15],data[16],data[17]);
    if (mat == 0) {
	opserr << "WARNING: failed to create SelfCenteringMaterial material\n";
	return 0;
    }

    return mat;
}

SelfCenteringMaterial::SelfCenteringMaterial(int tag, double K1, double K2, double K3,double SecondFusee, double epsAA, double epsBB, double ForcBB, double a33, double a44, double a55, double a66, double a77, double a88,
					     double fa, double b, double sd, 
					     double bd, double rb)
  : UniaxialMaterial(tag,MAT_TAG_SelfCentering), 
    k1(K1), k2(K2), k3(K3), SecondFuse(SecondFusee) , epsA(epsAA) , epsB(epsBB) , ForcB(ForcBB) , a3(a33) , a4(a44) , a5(a55) , a6(a66) , a7(a77) , a8(a88) , ActF(fa), beta(b), SlipDef(sd), BearDef(bd), rBear(rb)
{
  // Find Equivalent Slip Force
  ActDef = ActF / k1;
  SlipF = ActF + (SlipDef - ActDef) * k2;
  if (BearDef != 0) {
    if (SlipDef != 0 && SlipDef < BearDef) {
      BearF = SlipF;
    }
    else {
      BearF = ActF + (BearDef - ActDef) * k2;
    }
  }
  
  // Initialize variables
  this->revertToStart();
  
}

SelfCenteringMaterial::SelfCenteringMaterial()
  : UniaxialMaterial(0,MAT_TAG_SelfCentering),
    k1(0.0), k2(0.0), k3(0.0), SecondFuse(0.0) , epsA(0.0) , epsB(0.0) , ForcB(0.0) , a3(0.0) , a4(0.0) , a5(0.0) , a6(0.0) , a7(0.0) , a8(0.0) , ActF(0.0), beta(0.0), 
    SlipDef(0.0), BearDef(0.0), rBear(0.0)
{
  // Initialize variables
  SlipF = 0;
  ActDef = 0;
  BearF = 0;
  
  epsMaxnP = 0; 
  epsMaxnN = 0; 
  SecondFuseAP = 0; 
  SecondFuseAN = 0;
	
  x1p = 0; 
  y1p = 0; 
  m1p = 0; 
  x2p = 0; 
  y2p = 0; 
  m2p = 0;
  x3p = 0; 
  y3p = 0; 
  m3p = 0;
  x4p = 0; 
  y4p = 0; 
  kupp = 0;
	
  x1n = 0; 
  y1n = 0; 
  m1n = 0; 
  x2n = 0; 
  y2n = 0; 
  m2n = 0;
  x3n = 0; 
  y3n = 0; 
  m3n = 0;
  x4n = 0; 
  y4n = 0; 
  kupn = 0;
	
  xp = 0; 
  yp = 0;
  xn = 0; 
  yn = 0;

  x1tp = 0; 
  y1tp = 0;
  x3tp = 0; 
  y3tp = 0;
  xtp = 0;
	
  x1tn = 0; 
  y1tn = 0;
  x3tn = 0; 
  y3tn = 0;
  xtn = 0;
	
  
  this->revertToStart();
  
}

SelfCenteringMaterial::~SelfCenteringMaterial()
{
  
}

int 
SelfCenteringMaterial::setTrialStrain (double strain, double strainRate)
{
  
  diffStrain = strain - Cstrain;
  
  if (fabs(diffStrain) < DBL_EPSILON)
    return 0;
  
  // Set total strain
  Tstrain = strain;
  noSlipStrain = Tstrain - CslipStrain;
  

	// Positive Quadrant (Top Right) where strain >= 0
	if (noSlipStrain >= 0) {
	// Start of Positive General Parameters
	if (epsMaxnP < Cstrain) {    
    epsMaxnP = Cstrain;   
    }
	// First Cycle of Second Fuse
	if (SecondFuse != 0 ) {
	if (epsMaxnP >= epsA && SecondFuseAP == 0) {  // Positive    
    SecondFuseAP = 1;  
	x1p= ActF / k1 ;
	y1p= ActF ;
	x3p= epsA ;
	y3p= ActF + (epsA-ActF / k1)*k2 ;
	x2p= (1*x1p+99*x3p)/100 ;
	y2p= (1*y1p+99*y3p)/100 ;
	m1p=k2;
	m2p=k2;
	m3p=(a6-y3p)/(a5-epsA);
	x4p= epsMaxnP;
	y4p=m3p*(x4p-x3p)+y3p;
    }
	} // End of Positive General Parameters
	if (SecondFuse != 0 && SecondFuseAP != 0) { //Start of Positive FUSE 2
	
	if ((epsMaxnP) >= epsB)  {//Start of Positive FUSE 2 Collapse
	// if ((epsMaxnP) <= (ForcB/(2*100*k1)+epsB)){
	// Tstress = ForcB + (noSlipStrain - epsB) * (-100) * k1;
	// Ttangent = -100 * k1;			 
	// }
	 //else {
		     Tstress = 0.000001*k1 * noSlipStrain;
			 Ttangent = 0.000001*k1;
	 //}
	 }//End of Positive FUSE 2 Collapse
	 
	else { //Start of Positive FUSE 2 No Collapse
		xp=x4p-epsA;
		yp=xp/a3;
		// Positive Fuse 2 Linear Part
	if ( (fabs(noSlipStrain) <= ( ((1-beta) * ActF-.3226*yp) / k1 )) && (noSlipStrain >= 0) ) {
    Tstress = k1 * noSlipStrain;
    Ttangent = k1;
  }
  
  else if ( (noSlipStrain) > ( ((1-beta) * ActF-0.3226*yp) / k1 ) ) { // Positive Fuse 2 NonLinear Part
    
    // Positive Quadrant (Top Right) where strain >= 0
      // Positive Fuse 2 NonLinear Part
      // External Fuse Bolt Bearing
      if (BearDef != 0 && Tstrain > epsB) {
	Tstress = ForcB + (Tstrain - epsB) * rBear * k1;
	Ttangent = rBear * k1;
      }
      // Positive Fuse 2 NonLinear Part
      // External Fuse Slip
      else if (SlipDef != 0 && noSlipStrain > SlipDef) {
	Tstress = SlipF;
	TslipStrain = CslipStrain + diffStrain;
      }
      // Positive Fuse 2 NonLinear Part
      // Linear range movement (no upper or lower activation)
      else if ( (noSlipStrain >= ClowerStrainPos) && 
		(noSlipStrain <= CupperStrainPos) ) {
	
	Tstress = (noSlipStrain - CactivStrainPos) * k1;
	Ttangent = k1;
	
      }
	  // Positive Fuse 2 NonLinear Part
      // Upper Activation
      else if ( noSlipStrain > CupperStrainPos ) {
		  if ( (noSlipStrain >= x1p) && (noSlipStrain < x2p) ) {
			  kupp=m1p;
		  }
		  else if ( (noSlipStrain >= x2p) && (noSlipStrain < x3p) ) {
			  kupp=m2p;
		  }
		  else {
			  
			//  
		  if ( (epsMaxnP) < a5 ) {
			  m3p=(a6-(ActF + (epsA-ActF / k1)*k2))/(a5-epsA);
		  }
		  else if ( (epsMaxnP) < a7 ) {
			  m3p=(a8-a6)/(a7-a5);
		  }
		  else {
			  m3p=(ForcB-a8)/(epsB-a7);
		  }
			  
			//  
			  kupp=m3p;
		  }
		  
		  
	
	TupperStressPos = CupperStressPos + 
	  (noSlipStrain - CupperStrainPos) * kupp;
	TupperStrainPos = noSlipStrain;
	x1tp=TupperStrainPos;
	y1tp=TupperStressPos;
	x3tp=( ((1-beta) * ActF-0.3226*yp) / k1 );
	y3tp=( ((1-beta) * ActF-0.3226*yp) );
	xtp=((y3tp-y1tp)+k1*x1tp-k3*x3tp)/(k1-k3);

	TlowerStrainPos = xtp;

	TlowerStressPos = TupperStressPos - ((noSlipStrain-TlowerStrainPos)*k1);
	Tstress = TupperStressPos;
	TactivStrainPos = TupperStrainPos - Tstress / k1;
	
	Ttangent = kupp;
      }
	  // Positive Fuse 2 NonLinear Part
      // Lower Activation
      else { // Tstrain < ClowerStrainPos
	
	TlowerStressPos = ClowerStressPos +
	  (noSlipStrain - ClowerStrainPos) * k3;
	TlowerStrainPos = noSlipStrain;
	x1tp=TlowerStrainPos;
	y1tp=TlowerStressPos;
	kupp=m1p;
	x3tp=x1p;
	y3tp=y1p;
	xtp=((y3tp-y1tp)+k1*x1tp-kupp*x3tp)/(k1-kupp);
	
		if  (xtp >= x2p)  {
			  kupp=m2p;
			  x3tp=x2p;
			  y3tp=y2p;
			  xtp=((y3tp-y1tp)+k1*x1tp-kupp*x3tp)/(k1-kupp);
		  }
	    if  (xtp >= x3p) {
			//  
		  if ( (epsMaxnP) < a5 ) {
			  m3p=(a6-(ActF + (epsA-ActF / k1)*k2))/(a5-epsA);
		  }
		  else if ( (epsMaxnP) < a7 ) {
			  m3p=(a8-a6)/(a7-a5);
		  }
		  else {
			  m3p=(ForcB-a8)/(epsB-a7);
		  }
			  
			//  
			  kupp=m3p;
			  x3tp=x3p;
			  y3tp=y3p;
			  xtp=((y3tp-y1tp)+k1*x1tp-kupp*x3tp)/(k1-kupp);
		  }

	//TupperStrainPos = noSlipStrain + ((beta * ActF)+(k2-k3)*(noSlipStrain-( (1-beta) * ActF / k1 ))) / k1;
	TupperStrainPos = xtp;
	TupperStressPos = TlowerStressPos +((TupperStrainPos-noSlipStrain)*k1) ;
	Tstress = TlowerStressPos;
	TactivStrainPos = TlowerStrainPos - Tstress / k1;
	
	Ttangent = k3;
      }// Positive Fuse 2  End of 5 state
  } // End of Positive Fuse 1 NonLinear Part
  
  // Updating parameters
  	// if (epsMaxnP < noSlipStrain) {    
    // epsMaxnP = noSlipStrain;   
    // }

    xp=epsMaxnP-epsA;
	yp=xp/a3;
    x1p= (ActF-yp) / k1 ;
	y1p= ActF-yp ;
	x3p= epsMaxnP ;
	//  
		  if ( (epsMaxnP) < a5 ) {
			  m3p=(a6-(ActF + (epsA-ActF / k1)*k2))/(a5-epsA);
			  y3p= (ActF + (epsA-ActF / k1)*k2)+ m3p*(x3p-epsA);
		  }
		  else if ( (epsMaxnP) < a7 ) {
			  m3p=(a8-a6)/(a7-a5);
			  y3p= (ActF + (epsA-ActF / k1)*k2)+ ((a6-(ActF + (epsA-ActF / k1)*k2))/(a5-epsA))*(a5-epsA)+ m3p*(x3p-a5);
		  }
		  else {
			  m3p=(ForcB-a8)/(epsB-a7);
			  y3p= (ActF + (epsA-ActF / k1)*k2)+ ((a6-(ActF + (epsA-ActF / k1)*k2))/(a5-epsA))*(a5-epsA)+ ((a8-a6)/(a7-a5))*(a7-a5)+ m3p*(x3p-a7);
		  }
			  
	// 
	
	//y3p= ActF + (epsA-ActF / k1)*k2 ;
	//m3p=(ForcB-y3p)/(epsB-epsA);
	//y3p= (ActF + (epsA-ActF / k1)*k2)+ m3p*(x3p-epsA);
	
	x2p=a4*x3p;
	y2p=y1p+k2*(x2p-x1p);

	m1p=k2;
	m2p=(y3p-y2p)/(x3p-x2p);
	
	x4p= epsMaxnP ;
	y4p=y3p;
  // End of Updating parameters
  // Fuse 2 end
  // Fuse 2 end
	}    // End of IF Positive Fuse 2 No Collapse
	}	// End of IF Positive Fuse 2 
	else { //Start of Fuse 1
	// Positive Fuse 1 Linear Part
	if ( (fabs(noSlipStrain) <= ( (1-beta) * ActF / k1 )) && (noSlipStrain >= 0) ) {
    Tstress = k1 * noSlipStrain;
    Ttangent = k1; 
  }
  // Positive Fuse 1 NonLinear Part
  else if ( (noSlipStrain) > ( (1-beta) * ActF / k1 ) ) {
    
      // Positive Fuse 1 NonLinear Part
      // External Fuse Bolt Bearing
      if (BearDef != 0 && Tstrain > BearDef) {
	Tstress = BearF + (Tstrain - BearDef) * rBear * k1;
	Ttangent = rBear * k1;
      }
      // Positive Fuse 1 NonLinear Part
      // External Fuse Slip
      else if (SlipDef != 0 && noSlipStrain > SlipDef) {
	Tstress = SlipF;
	TslipStrain = CslipStrain + diffStrain;
      }
      // Positive Fuse 1 NonLinear Part
      // Linear range movement (no upper or
      //     lower activation)
      else if ( (noSlipStrain >= ClowerStrainPos) && 
		(noSlipStrain <= CupperStrainPos) ) {
	
	Tstress = (noSlipStrain - CactivStrainPos) * k1;
	Ttangent = k1;
	
      }
	  // Positive Fuse 1 NonLinear Part
      // Upper Activation
      else if ( noSlipStrain > CupperStrainPos ) {
	
	TupperStressPos = CupperStressPos + 
	  (noSlipStrain - CupperStrainPos) * k2;
	TupperStrainPos = noSlipStrain;
	// TlowerStrainPos = noSlipStrain - ((beta * ActF)+(k2-k3)*(noSlipStrain-(ActF / k1))) / k1;
	TlowerStrainPos = ((TupperStressPos)-(k1*noSlipStrain)+(k3/k1*(1-beta) * ActF)-((1-beta) * ActF))/(k3-k1);
	//TlowerStressPos = TupperStressPos - ((beta * ActF)+(k2-k3)*(noSlipStrain-(ActF / k1)));
	TlowerStressPos = TupperStressPos - ((noSlipStrain-TlowerStrainPos)*k1);
	Tstress = TupperStressPos;
	TactivStrainPos = TupperStrainPos - Tstress / k1;
	
	Ttangent = k2;
      }
	  // Positive Fuse 1 NonLinear Part
      // Lower Activation
      else { // Tstrain < ClowerStrainPos
	
	TlowerStressPos = ClowerStressPos +
	  (noSlipStrain - ClowerStrainPos) * k3;
	TlowerStrainPos = noSlipStrain;
	//TupperStrainPos = noSlipStrain + ((beta * ActF)+(k2-k3)*(noSlipStrain-( (1-beta) * ActF / k1 ))) / k1;
	TupperStrainPos = ((TlowerStressPos)+(k2*ActF/k1)-(ActF)-(k1*noSlipStrain))/(k2-k1);
	TupperStressPos = TlowerStressPos +((TupperStrainPos-noSlipStrain)*k1) ;
	Tstress = TlowerStressPos;
	TactivStrainPos = TlowerStrainPos - Tstress / k1;
	
	Ttangent = k3;
      } // End of 5 state

  } // End of Positive Fuse 1 NonLinear Part
	}   // End of Else Fuse 1
	}   // End of if (noSlipStrain >= 0) {
		
    // Negative Quadrant (Bottom Left) where strain < 0

	if (noSlipStrain < 0) {
	 // Start of Negative General Parameters
	 if (epsMaxnN > Cstrain) {    
     epsMaxnN = Cstrain;   
     }
	 if (SecondFuse != 0 ) {
	 if ((-1*epsMaxnN) >= epsA && SecondFuseAN == 0) {    // Negative  
     SecondFuseAN = 1;  
	 x1n= -1*ActF / k1 ;
	 y1n= -1*ActF ;
	 x3n= -1*epsA ;
	 y3n= -1*(ActF + (epsA-ActF / k1)*k2) ;
	 x2n= (1*x1n+99*x3n)/100 ;
	 y2n= (1*y1n+99*y3n)/100 ;
	 m1n=k2;
	 m2n=k2;
	 m3n=(a6-(ActF + (epsA-ActF / k1)*k2))/(a5-epsA);
	 x4n= epsMaxnN ;
	 y4n=m3n*(x4n-x3n)+y3n;
     }
	 } // End of Negative General Parameters
	 if (SecondFuse != 0 && SecondFuseAN != 0) {//Start of Negative FUSE 2 
	 if ((-1*epsMaxnN) >= epsB)  {//Start of Negative FUSE 2 Collapse
	 //if ((-1*epsMaxnN) <= (ForcB/(2*100*k1)+epsB)){
	//Tstress = -1*ForcB + (noSlipStrain + epsB) * (-100) * k1;
	//Ttangent = -100 * k1;			 
	 //}
	 //else {
		     Tstress = 0.000001*k1 * noSlipStrain;
			 Ttangent = 0.000001*k1;
	 //}
 
	 
	 }//End of Negative FUSE 2 Collapse
	 else {  //Start of Negative FUSE 2 No Collapse
		xn=x4n+epsA;
		yn=xn/a3;
		// Negative Fuse 2 Linear Part
	if ( (fabs(noSlipStrain) <= ( ((1-beta) * ActF+0.3226*yn) / k1 )) && (noSlipStrain < 0) ) {
    
    Tstress = k1 * noSlipStrain;
    Ttangent = k1;
    
  }
  
  else if ( (-1*noSlipStrain) > ( ((1-beta) * ActF+0.3226*yn) / k1 ) ) {// Negative Fuse 2 NonLinear Part
    
    // Negative Quadrant (Bottom Left) where strain < 0

      // Negative Fuse 2 NonLinear Part
      // External Fuse Bolt Bearing
      if (BearDef != 0 && Tstrain < (-1*epsB)) {
	Tstress = -1*ForcB + (Tstrain + epsB) * rBear * k1;
	Ttangent = rBear * k1;
      }
      // Negative Fuse 2 NonLinear Part
      // External Fuse Slip
      else if (SlipDef != 0 && noSlipStrain < (-1*SlipDef)) {
	Tstress = -1*SlipF;
	TslipStrain = CslipStrain + diffStrain;
      }
      // Negative Fuse 2 NonLinear Part
      // Linear range movement (no upper or
      //     lower activation)
      else if ( (noSlipStrain <= ClowerStrainNeg) && 
		(noSlipStrain >= CupperStrainNeg) ) {
	
	Tstress = (noSlipStrain - CactivStrainNeg) * k1;
	Ttangent = k1;
	
      }
	  // Negative Fuse 2 NonLinear Part
      // Upper Activation
      else if ( noSlipStrain < CupperStrainNeg ) {
		  if ( (noSlipStrain <= x1n) && (noSlipStrain > x2n) ) {
			  kupn=m1n;
		  }
		  else if ( (noSlipStrain <= x2n) && (noSlipStrain > x3n) ) {
			  kupn=m2n;
		  }
		  else {
			  //  
		  if ( (-1*epsMaxnN) < a5 ) {
			  m3n=(a6-(ActF + (epsA-ActF / k1)*k2))/(a5-epsA);
		  }
		  else if ( (-1*epsMaxnN) < a7 ) {
			  m3n=(a8-a6)/(a7-a5);
		  }
		  else {
			  m3n=(ForcB-a8)/(epsB-a7);
		  }
			  
			//
			  kupn=m3n;
		  }
		  
		  
	
	TupperStressNeg = CupperStressNeg + 
	  (noSlipStrain - CupperStrainNeg) * kupn;
	TupperStrainNeg = noSlipStrain;
	x1tn=TupperStrainNeg;
	y1tn=TupperStressNeg;
	x3tn=-1*(((1-beta) * ActF+0.3226*yn) / k1 );
	y3tn=-1*(((1-beta) * ActF+0.3226*yn) );
	xtn=((y3tn-y1tn)+k1*x1tn-k3*x3tn)/(k1-k3);

	TlowerStrainNeg = xtn;

	TlowerStressNeg = TupperStressNeg - ((noSlipStrain-TlowerStrainNeg)*k1);
	Tstress = TupperStressNeg;
	TactivStrainNeg = TupperStrainNeg - Tstress / k1;
	
	Ttangent = kupn;
      }
	  // Negative Fuse 2 NonLinear Part
      // Lower Activation
      else { // Tstrain > ClowerStrainNeg
	
	TlowerStressNeg = ClowerStressNeg +
	  (noSlipStrain - ClowerStrainNeg) * k3;
	TlowerStrainNeg = noSlipStrain;
	x1tn=TlowerStrainNeg;
	y1tn=TlowerStressNeg;
	kupn=m1n;
	x3tn=x1n;
	y3tn=y1n;
	xtn=((y3tn-y1tn)+k1*x1tn-kupn*x3tn)/(k1-kupn);
	
		if  (xtn <= x2n)  {
			  kupn=m2n;
			  x3tn=x2n;
			  y3tn=y2n;
			  xtn=((y3tn-y1tn)+k1*x1tn-kupn*x3tn)/(k1-kupn);
		  }
	    if  (xtn <= x3n) {
			//  
		  if ( (-1*epsMaxnN) < a5 ) {
			  m3n=(a6-(ActF + (epsA-ActF / k1)*k2))/(a5-epsA);
		  }
		  else if ( (-1*epsMaxnN) < a7 ) {
			  m3n=(a8-a6)/(a7-a5);
		  }
		  else {
			  m3n=(ForcB-a8)/(epsB-a7);
		  }
			  
			//
			  kupn=m3n;
			  x3tn=x3n;
			  y3tn=y3n;
			  xtn=((y3tn-y1tn)+k1*x1tn-kupn*x3tn)/(k1-kupn);
		  }

	//TupperStrainPos = noSlipStrain + ((beta * ActF)+(k2-k3)*(noSlipStrain-( (1-beta) * ActF / k1 ))) / k1;
	TupperStrainNeg = xtn;
	TupperStressNeg = TlowerStressNeg +((TupperStrainNeg-noSlipStrain)*k1) ;
	Tstress = TlowerStressNeg;
	TactivStrainNeg = TlowerStrainNeg - Tstress / k1;
	
	Ttangent = k3;
      }   // Negative Fuse 2  End of 5 state
   
  }    // End of Negative Fuse 1 NonLinear Part
  
  // Updating parameters
  	// if (epsMaxnN > noSlipStrain) {    
    // epsMaxnN = noSlipStrain;   
    // }

    xn=epsMaxnN+epsA;
	yn=xn/a3;
    x1n= -1*(ActF+yn) / k1 ;
	y1n= -1*(ActF+yn) ;
	x3n= epsMaxnN ;
	//  
		  if ( (-1*epsMaxnN) < a5 ) {
			  m3n=(a6-(ActF + (epsA-ActF / k1)*k2))/(a5-epsA);
			  y3n= -1*(ActF + (epsA-ActF / k1)*k2)+ m3n*(x3n+epsA);
		  }
		  else if ( (-1*epsMaxnN) < a7 ) {
			  m3n=(a8-a6)/(a7-a5);
			  y3n= -1*(ActF + (epsA-ActF / k1)*k2)+ ((a6-(ActF + (epsA-ActF / k1)*k2))/(a5-epsA))*(-1*a5+epsA)+ m3n*(x3n+a5);
		  }
		  else {
			  m3n=(ForcB-a8)/(epsB-a7);
			  y3n= -1*(ActF + (epsA-ActF / k1)*k2)+ ((a6-(ActF + (epsA-ActF / k1)*k2))/(a5-epsA))*(-1*a5+epsA)+ ((a8-a6)/(a7-a5))*(-1*a7+a5)+ m3n*(x3n+a7);
		  }
			  
	//
	//y3n= -1*(ActF + (epsA-ActF / k1)*k2) ;
	//m3n=(-1*ForcB-y3n)/(-1*epsB+epsA);
	//y3n= (-1*ActF - (epsA-ActF / k1)*k2)+ m3n*(x3n+epsA);
	x2n=a4*x3n;
	y2n=y1n+k2*(x2n-x1n);

	m1n=k2;
	m2n=(y3n-y2n)/(x3n-x2n);
	
	x4n= epsMaxnN ;
	y4n=y3n;
  // End of Updating parameters
  // Fuse 2 end
    
 
		// FUSE 2
	}     // End of IF Negative Fuse 2 No Collapse
	 }   // End of IF Negative Fuse 2 
	else { //Start of Fuse 1
	// Negative Fuse 1 Linear Part
	  if ( (fabs(noSlipStrain) <= ( (1-beta) * ActF / k1 )) && (noSlipStrain < 0) ) {
    
    Tstress = k1 * noSlipStrain;
    Ttangent = k1;
    
  }
	// Negative Fuse 1 NonLinear Part
	else if ( (-1*noSlipStrain) > ( (1-beta) * ActF / k1 ) ) {   
	   // Tstrain < 0)
      // Negative Fuse 1 NonLinear Part
      // External Fuse Bolt Bearing
      if (BearDef != 0 && Tstrain < -BearDef) {
	Tstress = -BearF + (Tstrain + BearDef) * rBear * k1;
	Ttangent = rBear * k1;
      }
      // Negative Fuse 1 NonLinear Part
      // External Fuse Slip
      else if (SlipDef != 0 && noSlipStrain < -SlipDef) {
	Tstress = -SlipF;
	TslipStrain = CslipStrain + diffStrain;
      }
      // Negative Fuse 1 NonLinear Part
      // Linear range movement (no upper or
      //     lower activation)
      else if ( (noSlipStrain <= ClowerStrainNeg) && 
		(noSlipStrain >= CupperStrainNeg) ) {
	
	Tstress = (noSlipStrain - CactivStrainNeg) * k1;
	Ttangent = k1;
	
      }
	  // Negative Fuse 1 NonLinear Part
      // Upper Activation
      else if ( noSlipStrain < CupperStrainNeg ) {
	
	TupperStressNeg = CupperStressNeg + 
	  (noSlipStrain - CupperStrainNeg) * k2;
	TupperStrainNeg = noSlipStrain;
	TlowerStrainNeg = ((TupperStressNeg)-(k1*noSlipStrain)-(k3/k1*(1-beta) * ActF)+((1-beta) * ActF))/(k3-k1);
	TlowerStressNeg = TupperStressNeg - ((noSlipStrain-TlowerStrainNeg)*k1);
	Tstress = TupperStressNeg;
	TactivStrainNeg = TupperStrainNeg - Tstress / k1;
	
	Ttangent = k2;
      }
	  // Negative Fuse 1 NonLinear Part
      // Lower Activation
      else { // Tstrain > ClowerStrainNeg
	
	TlowerStressNeg = ClowerStressNeg +
	  (noSlipStrain - ClowerStrainNeg) * k3;
	TlowerStrainNeg = noSlipStrain;
	TupperStrainNeg = ((TlowerStressNeg)-(k2*ActF/k1)+(ActF)-(k1*noSlipStrain))/(k2-k1);
	TupperStressNeg = TlowerStressNeg +((TupperStrainNeg-noSlipStrain)*k1) ;
	Tstress = TlowerStressNeg;
	TactivStrainNeg = TlowerStrainNeg - Tstress / k1;
	
	Ttangent = k3;
      } // End of 5 state
    
    
  }  // End of Negative Fuse 1 NonLinear Part
	}   // End of Else Fuse 1
	}   // End of if (noSlipStrain < 0) {
  return 0;
}      // End of setTrialStrain

double 
SelfCenteringMaterial::getStress(void)
{
  return Tstress;
}

double 
SelfCenteringMaterial::getTangent(void)
{
  return Ttangent;
}

double 
SelfCenteringMaterial::getStrain(void)
{
  return Tstrain;
}

int 
SelfCenteringMaterial::commitState(void)
{
  // Commit trial history variables
  CactivStrainPos = TactivStrainPos;
  CactivStrainNeg = TactivStrainNeg;
  CslipStrain = TslipStrain;
  CupperStrainPos = TupperStrainPos;
  ClowerStrainPos = TlowerStrainPos;	
  CupperStressPos = TupperStressPos;
  ClowerStressPos = TlowerStressPos;
  CupperStrainNeg = TupperStrainNeg;
  ClowerStrainNeg = TlowerStrainNeg;
  CupperStressNeg = TupperStressNeg;
  ClowerStressNeg = TlowerStressNeg;
  Cstrain = Tstrain;
  Cstress = Tstress;
  Ctangent = Ttangent;
  
  return 0;
}

int 
SelfCenteringMaterial::revertToLastCommit(void)
{
  Tstrain = Cstrain;
  Tstress = Cstress;
  Ttangent = Ctangent;

  return 0;
}

int 
SelfCenteringMaterial::revertToStart(void)
{
  // Reset committed history variables
  CactivStrainPos = 0.0;
  CactivStrainNeg = 0.0;
  CslipStrain = 0.0;
  CupperStrainPos = ActDef;
  ClowerStrainPos = (1-beta) * ActDef;	
  CupperStressPos = ActF;
  ClowerStressPos = (1-beta) * ActF;
  CupperStrainNeg = -CupperStrainPos;
  ClowerStrainNeg = -ClowerStrainPos;
  CupperStressNeg = -CupperStressPos;
  ClowerStressNeg = -ClowerStressPos;
  
  // Reset trial history variables
  TactivStrainPos = 0.0;
  TactivStrainNeg = 0.0;
  TslipStrain = 0.0;
  TupperStrainPos = ActDef;
  TlowerStrainPos = (1-beta) * ActDef;	
  TupperStressPos = ActF;
  TlowerStressPos = (1-beta) * ActF;
  TupperStrainNeg = -CupperStrainPos;
  TlowerStrainNeg = -ClowerStrainPos;
  TupperStressNeg = -CupperStressPos;
  TlowerStressNeg = -ClowerStressPos;
  
  // Initialize state variables
  Tstrain = 0.0;
  Tstress = 0.0;
  Ttangent = k1;
  
  Cstrain = 0.0;
  
  return 0;
}

UniaxialMaterial *
SelfCenteringMaterial::getCopy(void)
{
  SelfCenteringMaterial *theCopy =
    new SelfCenteringMaterial(this->getTag(), k1, k2, k3, SecondFuse, epsA, epsB, ForcB, a3, a4, a5, a6, a7, a8, ActF, beta,
			      SlipDef, BearDef, rBear);
  
  // Copy committed history variables
  theCopy->CactivStrainPos = CactivStrainPos;
  theCopy->CactivStrainNeg = CactivStrainNeg;
  theCopy->CslipStrain = CslipStrain;
  theCopy->CupperStrainPos = CupperStrainPos;
  theCopy->ClowerStrainPos = ClowerStrainPos;	
  theCopy->CupperStressPos = CupperStressPos;
  theCopy->ClowerStressPos = ClowerStressPos;
  theCopy->CupperStrainNeg = CupperStrainNeg;
  theCopy->ClowerStrainNeg = ClowerStrainNeg;
  theCopy->CupperStressNeg = CupperStressNeg;
  theCopy->ClowerStressNeg = ClowerStressNeg;
  
  // Copy trial history variables
  theCopy->TactivStrainPos = TactivStrainPos;
  theCopy->TactivStrainNeg = TactivStrainNeg;
  theCopy->TslipStrain = TslipStrain;
  theCopy->TupperStrainPos = TupperStrainPos;
  theCopy->TlowerStrainPos = TlowerStrainPos;	
  theCopy->TupperStressPos = TupperStressPos;
  theCopy->TlowerStressPos = TlowerStressPos;
  theCopy->TupperStrainNeg = TupperStrainNeg;
  theCopy->TlowerStrainNeg = TlowerStrainNeg;
  theCopy->TupperStressNeg = TupperStressNeg;
  theCopy->TlowerStressNeg = TlowerStressNeg;
  
  // Copy trial state variables
  theCopy->Tstrain = Tstrain;
  theCopy->Tstress = Tstress;
  theCopy->Ttangent = Ttangent;
  
  theCopy->Cstrain = Cstrain;
  
  return theCopy;
}

int 
SelfCenteringMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(79);
  
  data(0) = this->getTag();
  data(1) = k1;
  data(2) = k2;
  data(3) = k3;
  data(4) = SecondFuse;
  data(5) = epsA;
  data(6) = epsB;
  data(7) = ForcB;
  data(8) = a3;
  data(9) = a4;
  data(10) = a5;
  data(11) = a6;
  data(12) = a7;
  data(13) = a8;
  data(14) = ActF;
  data(15) = beta;
  data(16) = rBear;
  data(17) = SlipDef;
  data(18) = BearDef;
  data(19) = SlipF;
  data(20) = ActDef;
  data(21) = BearF;
  data(22) = epsMaxnP;
  data(23) = epsMaxnN;
  data(24) = SecondFuseAP;
  data(25) = SecondFuseAN;	
  data(26) = x1p;
  data(27) = y1p;
  data(28) = m1p;
  data(29) = x2p;
  data(30) = y2p;
  data(31) = m2p;
  data(32) = x3p;
  data(33) = y3p;
  data(34) = m3p;
  data(35) = x4p;
  data(36) = y4p;
  data(37) = kupp;	
  data(38) = x1n;
  data(39) = y1n;
  data(40) = m1n;
  data(41) = x2n;
  data(42) = y2n;
  data(43) = m2n;
  data(44) = x3n;
  data(45) = y3n;
  data(46) = m3n;
  data(47) = x4n;
  data(48) = y4n;
  data(49) = kupn;	
  data(50) = xp;
  data(51) = yp;
  data(52) = xn;
  data(53) = yn;
  data(54) = x1tp;
  data(55) = y1tp;
  data(56) = x3tp;
  data(57) = y3tp;
  data(58) = xtp;	
  data(59) = x1tn;
  data(60) = y1tn;
  data(61) = x3tn;
  data(62) = y3tn;
  data(63) = xtn;
  data(64) = CactivStrainPos;
  data(65) = CactivStrainNeg;
  data(66) = CslipStrain;
  data(67) = CupperStrainPos;
  data(68) = ClowerStrainPos;
  data(69) = CupperStressPos;
  data(70) = ClowerStressPos;
  data(71) = CupperStrainNeg;
  data(72) = ClowerStrainNeg;
  data(73) = CupperStressNeg;
  data(74) = ClowerStressNeg;
  data(75) = Tstrain;
  data(76) = Tstress;
  data(77) = Ttangent;
  data(78) = Cstrain;
  
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "SelfCenteringMaterial::sendSelf() - failed to send data\n";

  return res;
}

int 
SelfCenteringMaterial::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(79);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "SelfCenteringMaterial::recvSelf() - failed to receive data\n";
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
    k1 = data(1);
	k2 = data(2);
	k3 = data(3);
	SecondFuse = data(4);
	epsA = data(5);
	epsB = data(6);
	ForcB = data(7);
	a3 = data(8);
	a4 = data(9);
	a5 = data(10);
	a6 = data(11);
	a7 = data(12);
	a8 = data(13);
	ActF = data(14);
	beta = data(15);
	rBear = data(16);
	SlipDef = data(17);
	BearDef = data(18);
	SlipF = data(19);
	ActDef = data(20);
	BearF = data(21);
	epsMaxnP = data(22);
	epsMaxnN = data(23);
	SecondFuseAP = data(24);
	SecondFuseAN = data(25);	
	x1p = data(26);
	y1p = data(27);
	m1p = data(28);
	x2p = data(29);
	y2p = data(30);
	m2p = data(31);
	x3p = data(32);
	y3p = data(33);
	m3p = data(34);
	x4p = data(35);
	y4p = data(36);
	kupp = data(37);	
	x1n = data(38);
	y1n = data(39);
	m1n = data(40);
	x2n = data(41);
	y2n = data(42);
	m2n = data(43);
	x3n = data(44);
	y3n = data(45);
	m3n = data(46);
	x4n = data(47);
	y4n = data(48);
	kupn = data(49);	
	xp = data(50);
	yp = data(51);
	xn = data(52);
	yn = data(53);
	x1tp = data(54);
	y1tp = data(55);
	x3tp = data(56);
	y3tp = data(57);
	xtp = data(58);	
	x1tn = data(59);
	y1tn = data(60);
	x3tn = data(61);
	y3tn = data(62);
	xtn = data(63);
	CactivStrainPos = data(64);
	CactivStrainNeg = data(65);
	CslipStrain = data(66);
	CupperStrainPos = data(67);
	ClowerStrainPos = data(68);
	CupperStressPos = data(69);
	ClowerStressPos = data(70);
	CupperStrainNeg = data(71);
	ClowerStrainNeg = data(72);
	CupperStressNeg = data(73);
	ClowerStressNeg = data(74);
	Tstrain = data(75);
	Tstress = data(76);
	Ttangent = data(77);
	Cstrain = data(78);
  }
    
  return res;
}

void 
SelfCenteringMaterial::Print(OPS_Stream &s, int flag)
{
    s << "SelfCenteringMaterial, tag: " << this->getTag() << endln;
    s << "  k1: " << k1 << endln;
    s << "  k2: " << k2 << endln;
	s << "  k3: " << k3 << endln;
	s << "  SecondFuse: " << SecondFuse << endln;
	s << "  epsA: " << epsA << endln;
	s << "  epsB: " << epsB << endln;
	s << "  ForcB: " << ForcB << endln;
	s << "  a3: " << a3 << endln;
	s << "  a4: " << a4 << endln;
	s << "  a5: " << a5 << endln;
	s << "  a6: " << a6 << endln;
	s << "  a7: " << a7 << endln;
	s << "  a8: " << a8 << endln;
    s << "  ActF: " << ActF << endln;
    s << "  beta: " << beta << endln;
    s << "  rBear: " << rBear << endln;
    s << "  SlipDef: " << SlipDef << endln;
    s << "  BearDef: " << BearDef << endln;
}

