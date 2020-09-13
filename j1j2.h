#ifndef __J1J2_H
#define __J1J2_H

#include "itensor/all.h"

namespace itensor {

MPS inline
computeGroundState(SpinHalf const& sites, 
                   Real Jp)
    {
    auto ampo = AutoMPO(sites);
    auto N = sites.N();
    Real J1 = 1.0;
    Real Ju = 0.76;
        
    for (int j = 1; j < N-2; j += 2)
    {
        ampo += 0.5*J1,"S+",j,"S-",j+2;
        ampo += 0.5*J1,"S-",j,"S+",j+2;
        ampo +=     J1,"Sz",j,"Sz",j+2;
            
        ampo += 0.5*J1,"S+",j+1,"S-",j+3;
        ampo += 0.5*J1,"S-",j+1,"S+",j+3;
        ampo +=     J1,"Sz",j+1,"Sz",j+3;
    }
    for (int j = 1; j < N; j += 2)
    {
        ampo += 0.5*Jp,"S+",j,"S-",j+1;
        ampo += 0.5*Jp,"S-",j,"S+",j+1;
        ampo +=     Jp,"Sz",j,"Sz",j+1;
    }
        
    for (int j = 1; j < N-2; j += 2)
    {
        //Jx = -0.5Jp
        ampo += -0.25*Jp,"S+",j,"S-",j+3;
        ampo += -0.25*Jp,"S-",j,"S+",j+3;
        ampo += -0.5*Jp,"Sz",j,"Sz",j+3;
            
        ampo += -0.25*Jp,"S+",j+1,"S-",j+2;
        ampo += -0.25*Jp,"S-",j+1,"S+",j+2;
        ampo += -0.5*Jp,"Sz",j+1,"Sz",j+2;
    }
        
    for (int j = 1; j < N-2; j += 2)
    {
        ampo += -0.25*Ju,"S+",j,"S-",j+2,"S+",j+1,"S-",j+3;
        ampo += -0.25*Ju,"S+",j,"S-",j+2,"S-",j+1,"S+",j+3;
        ampo += -0.25*Ju,"S-",j,"S+",j+2,"S+",j+1,"S-",j+3;
        ampo += -0.25*Ju,"S-",j,"S+",j+2,"S-",j+1,"S+",j+3;
            
        ampo += -0.5*Ju,"S+",j,"S-",j+2,"Sz",j+1,"Sz",j+3;
        ampo += -0.5*Ju,"S-",j,"S+",j+2,"Sz",j+1,"Sz",j+3;
        ampo += -0.5*Ju,"Sz",j,"Sz",j+2,"S+",j+1,"S-",j+3;
        ampo += -0.5*Ju,"Sz",j,"Sz",j+2,"S-",j+1,"S+",j+3;
            
        ampo += -1.0*Ju,"Sz",j,"Sz",j+2,"Sz",j+1,"Sz",j+3;
    }
    auto H = MPO(ampo);

    auto psi = MPS(sites);

    auto sweeps = Sweeps(8);
    sweeps.maxm() = 50,100,200,300,400,500,600,700;
       //,1400,1500;
      //  50,100,200,300,400;
    sweeps.cutoff() = 1E-9;

    println("Starting ground state calculation for Jp = ",Jp);

    dmrg(psi,H,sweeps,{"Quiet",true});

    println("Done with ground state calculation.");

    return psi;
    }

} //namespace itensor

#endif
