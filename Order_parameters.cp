#include "j1j2.h"
#include <fstream>
#include <iomanip>

using namespace std;

using std::vector;
using namespace itensor;

int main()
{
      //for(int Nx = 16; Nx< 100; Nx += 16)
       //{   
  int Nxset[1] = {48}; 
  //int Nxset[1] = {200};  
  // int Nx = 96;
    //16,24,32,48,64,80,96
 for(int i = 0; i < 1; i += 1)
  { 
    auto Nx = Nxset[i];
    auto N = 2*Nx;    
        
    auto sites = SpinHalf(N);

    //auto Jps = vector<Real>();
    //auto dimer = vector<Real>();

    for(Real Jp = -0.5; Jp > -2.6; Jp -= 0.1)
        {
        //
        // Compute ground state using
        // "black box" routine (or have a look at j1j2.h)
        //
        MPS psi = computeGroundState(sites,Jp);

        Real D = 0.;
        //
        // We will add code below to
        // measure the dimer order parameter:
        //
        // D = B_(N/2) - B_(N/2-2)
        //

        //
        // Compute:  <B_(Nx/2)>
        //
        auto i1 = N/2;
        auto i2 = N/2 + 2;
        ITensor Szi1 = sites.op("Sz",i1);
        ITensor Szi2 = sites.op("Sz",i2);
        ITensor Spi1 = sites.op("Sp",i1);
        ITensor Spi2 = sites.op("Sp",i2);
        ITensor Smi1 = sites.op("Sm",i1);
        ITensor Smi2 = sites.op("Sm",i2);

        psi.position(i1);
            
        //index linking i to i+1:
        auto ira = commonIndex(psi.A(i1),psi.A(i1+1),Link);
            
        auto Cazz = psi.A(i1)*Szi1*dag(prime(psi.A(i1),Site,ira));
        auto Capm = 0.5*psi.A(i1)*Spi1*dag(prime(psi.A(i1),Site,ira));
        auto Camp = 0.5*psi.A(i1)*Smi1*dag(prime(psi.A(i1),Site,ira));
            
        for(int k = i1+1; k < i2; ++k)
        {
            Cazz *= psi.A(k);
            Cazz *= dag(prime(psi.A(k),Link));
            
            Capm *= psi.A(k);
            Capm *= dag(prime(psi.A(k),Link));
            
            Camp *= psi.A(k);
            Camp *= dag(prime(psi.A(k),Link));
        }
        Cazz *= psi.A(i2);
        Cazz *= Szi2;
        Capm *= psi.A(i2);
        Capm *= Smi2;
        Camp *= psi.A(i2);
        Camp *= Spi2;
            
        //index linking j to j-1:
        auto jla = commonIndex(psi.A(i2),psi.A(i2-1),Link);
        Cazz *= dag(prime(psi.A(i2),jla,Site));
        Capm *= dag(prime(psi.A(i2),jla,Site));
        Camp *= dag(prime(psi.A(i2),jla,Site));
        auto Ca = Cazz + Capm + Camp;
        D += Ca.real();

        //
        // Compute: - <B_(Nx/2-2)>
        //
        auto i3 = N/2 - 2;
        ITensor Szi3 = sites.op("Sz",i3);
        ITensor Spi3 = sites.op("Sp",i3);
        ITensor Smi3 = sites.op("Sm",i3);
            
        psi.position(i3);
    
        //index linking i to i+1:
        auto irb = commonIndex(psi.A(i3),psi.A(i3+1),Link);
            
        auto Cbzz = psi.A(i3)*Szi3*dag(prime(psi.A(i3),Site,irb));
        auto Cbpm = 0.5*psi.A(i3)*Spi3*dag(prime(psi.A(i3),Site,irb));
        auto Cbmp = 0.5*psi.A(i3)*Smi3*dag(prime(psi.A(i3),Site,irb));
            
        for(int k = i3+1; k < i1; ++k)
        {
            Cbzz *= psi.A(k);
            Cbzz *= dag(prime(psi.A(k),Link));
                
            Cbpm *= psi.A(k);
            Cbpm *= dag(prime(psi.A(k),Link));
                
            Cbmp *= psi.A(k);
            Cbmp *= dag(prime(psi.A(k),Link));
        }
        Cbzz *= psi.A(i1);
        Cbzz *= Szi1;
        Cbpm *= psi.A(i1);
        Cbpm *= Smi1;
        Cbmp *= psi.A(i1);
        Cbmp *= Spi1;
            
        //index linking j to j-1:
        auto jlb = commonIndex(psi.A(i1),psi.A(i1-1),Link);
        Cbzz *= dag(prime(psi.A(i1),jlb,Site));
        Cbpm *= dag(prime(psi.A(i1),jlb,Site));
        Cbmp *= dag(prime(psi.A(i1),jlb,Site));
        auto Cb = Cbzz + Cbpm + Cbmp;
        D += ((-1.0)*Cb).real();

            
        //
        // Compute: - <B_(Nx/2+2)>
        //
        //auto i3 = Nx+2;
        //psi.position(i3);
        //auto B3 = makeB(sites,i3);
        //auto wf3 = psi.A(i3)*psi.A(i3+2);
        //TODO: ADD CODE to compute <wf2|B2|wf2>
        //D +=  ((-0.5)*dag(prime(wf3,Site)) * B3 * wf3).real();
            
        printfln("Jp = %.5f, Dimer order parameter = %.10f",Jp,D);
        printfln("---------------------------------------------");
            
        ofstream FileOutput;
        FileOutput.open(format("OP_rough_N=%d_Ju=-1.txt",Nx),ios::out|ios::app);
        FileOutput << fixed << setprecision(10) << Nx << " " << Jp << " " << D << endl;
        FileOutput.close();

            
        }
  }


 return 0;
 }

