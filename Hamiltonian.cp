#include "itensor/all.h"
#include <fstream>
#include <iomanip>

using namespace std;
using namespace itensor;

int main()
{
    int Nx = 24;
    //24,28,32,40,48,64,80,96; 88; 72,56
    auto N = 2*Nx;
    auto J1 = 1.0;
    auto Jp = 0.1;
    //auto Jp = 0.1;
    //auto Ju = -0.55;
    auto Jx = -0.05;
    //auto Ju = 0.1;
    
for(Real Ju = 0.0; Ju <= 2.0; Ju += 0.1)
 {
    auto sites = SpinHalf(N);
    
    auto ampo = AutoMPO(sites);
    for (int j = 1; j <= N-3; j += 2)
    {
        ampo += J1/2,"S+",j,"S-",j+2;
        ampo += J1/2,"S-",j,"S+",j+2;
        ampo +=   J1,"Sz",j,"Sz",j+2;
        
        ampo += J1/2,"S+",j+1,"S-",j+3;
        ampo += J1/2,"S-",j+1,"S+",j+3;
        ampo +=   J1,"Sz",j+1,"Sz",j+3;
    }
    for (int j = 1; j <= N-1; j += 2)
    {
        ampo += Jp/2,"S+",j,"S-",j+1;
        ampo += Jp/2,"S-",j,"S+",j+1;
        ampo +=   Jp,"Sz",j,"Sz",j+1;
    }
    for (int j = 1; j <= N-3; j += 2)
    {
        ampo += Jx/2,"S+",j,"S-",j+3;
        ampo += Jx/2,"S-",j,"S+",j+3;
        ampo +=   Jx,"Sz",j,"Sz",j+3;
        
        ampo += Jx/2,"S+",j+1,"S-",j+2;
        ampo += Jx/2,"S-",j+1,"S+",j+2;
        ampo +=   Jx,"Sz",j+1,"Sz",j+2;
    }
    for (int j = 1; j <= N-3; j += 2)
    {
        ampo += -Ju/4,"S+",j,"S-",j+2,"S+",j+1,"S-",j+3;
        ampo += -Ju/4,"S+",j,"S-",j+2,"S-",j+1,"S+",j+3;
        ampo += -Ju/4,"S-",j,"S+",j+2,"S+",j+1,"S-",j+3;
        ampo += -Ju/4,"S-",j,"S+",j+2,"S-",j+1,"S+",j+3;
        
        ampo += -Ju/2,"S+",j,"S-",j+2,"Sz",j+1,"Sz",j+3;
        ampo += -Ju/2,"S-",j,"S+",j+2,"Sz",j+1,"Sz",j+3;
        ampo += -Ju/2,"Sz",j,"Sz",j+2,"S+",j+1,"S-",j+3;
        ampo += -Ju/2,"Sz",j,"Sz",j+2,"S-",j+1,"S+",j+3;
        
        ampo +=   -Ju,"Sz",j,"Sz",j+2,"Sz",j+1,"Sz",j+3;
    }
    auto H = MPO(ampo);
   
    printfln("Ju = %.10f",Ju);   
 
    auto sweeps = Sweeps(5);
    sweeps.maxm() = 50,100,200,300,400;
    sweeps.cutoff() = 1E-10;
    println(sweeps);
    
    //we can also set the initial sate as a ramdom initial state
    auto psi0 = MPS(sites);
    
    //auto init = InitState(sites);
    //for (int i = 1; i <= N; ++i)
    //{
      //  if(i%2 == 1) init.set(i,"Up");
      //  else         init.set(i,"Dn");
    //}
    //auto psi = MPS(init);
    
    auto en0 = dmrg(psi0,H,sweeps,{"Quiet",true});
    
    println("\n----------------------\n");
    
    //
    // Make a vector of previous wavefunctions;
    // code will penalize future wavefunctions
    // for having any overlap with these
    //
    auto wfs = std::vector<MPS>(1);
    wfs.at(0) = psi0;
    
    auto psi1 = MPS(sites);
    
    //
    // Here the Weight option sets the energy penalty for
    // psi1 having any overlap with psi0
    //
    auto en1 = dmrg(psi1,H,wfs,sweeps,{"Quiet=",true,"Weight=",20.0});
    
    //
    // Print the final energies reported by DMRG
    //
    printfln("\nGround State Energy = %.10f",en0);
    printfln("\nExcited State Energy = %.10f",en1);
    
    //
    //
    // (The DMRG gap will have finite-size corrections.)
    //
    printfln("\nDMRG energy gap = %.10f",en1-en0);
    //printfln("\nTheoretical gap = %.10f",2*std::fabs(h-1));
    
    //
    // The overlap <psi0|psi1> should be very close to zero
    //
    printfln("\nOverlap <psi0|psi1> = %.2E",overlap(psi0,psi1));
    
    //writeToFile(format("GS_energy_%d",N),en0);
    
    ofstream FileOutput;
    //FileOutput.open (format("GS_energy_%d.txt",N),ios::out|ios::app);
    //FileOutput.open(format("Energy_Jp=%d_Jx=%d_Ju=%d.txt",Jp,Jx,Ju),ios::out|ios::app);
    //FileOutput << fixed << setprecision(10) << Nx << " " << en0 << " " << en1 << " " << en1-en0 << endl;
    FileOutput.open(format("Energy_Jp=%d_Jx=%d_JuChange.txt",Jp,Jx),ios::out|ios::app);
    FileOutput << fixed << setprecision(10) << Nx << " " << -Ju << " " << en0 << " " <<  en1 << " " << en1-en0 << endl;
    FileOutput.close();
 }
    
    
    return 0;
}
