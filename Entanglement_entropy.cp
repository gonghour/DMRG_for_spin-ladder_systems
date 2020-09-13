#include "j1j2.h"
#include <fstream>
#include <iomanip>

using namespace std;

using std::vector;
using namespace itensor;


int main()
{
   // int Nxset[7] = {256,128,64,32,24,400,512}; 
   int Nxset[1] = {512}; 
    Real Jp = -0.2;

    for(int a = 0 ; a < 1; a += 1)
    {
    
        auto Nx = Nxset[a];
        auto N = 2*Nx;
        auto sites = SpinHalf(N);
        //
        // Compute ground state using
        // "black box" routine (or have a look at j1j2.h)
        //
        MPS psi = computeGroundState(sites,Jp);
        
        //Given an MPS or IQMPS called "psi",
        //and some particular bond "b" (1 <= b < psi.N())
        //across which we want to compute the von Neumann entanglement
            
        //"Gauge" the MPS to site b
        auto i1 = N/2;
        psi.position(i1);
            
        //Here assuming an MPS of ITensors, but same code works
        //for IQMPS by replacing ITensor -> IQTensor
            
        //Compute t:wo-site wavefunction for sites (b,b+1)
        ITensor wf = psi.A(i1)*psi.A(i1+1);
            
        //SVD this wavefunction to get the spectrum
        //of density-matrix eigenvalues
        auto U = psi.A(i1);
        ITensor S,V;
        auto spectrum = svd(wf,U,S,V);
            
        //Apply von Neumann formula
        //spectrum.eigs() is a Vector containing
        //the density matrix eigenvalues
        //(squares of the singular values)
        Real SvN = 0.;
        for(auto p : spectrum.eigs())
        {
            if(p > 1E-12) SvN += -p*log(p);
        }
        printfln("Across bond b=%d, SvN = %.10f",i1,SvN);
        
        
        ofstream FileOutput;
        FileOutput.open(format("EE_BD=1500_Jp=%d.txt",Jp),ios::out|ios::app);
        FileOutput << fixed << setprecision(10) << Nx  << " " << SvN << endl;
        FileOutput.close();
        //FileOutput.open(format("Correlator_SS_N=%d.txt",N),ios::out|ios::app);
        //FileOutput << fixed << setprecision(10)  << i << " " << G << endl;
        //FileOutput.close();
            
            
   }


  return 0;
}

