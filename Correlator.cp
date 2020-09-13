#include "j1j2.h"
#include <fstream>
#include <iomanip>

using namespace std;

using std::vector;
using namespace itensor;

int main()
    {
    int Nx = 128;
    //24,32,48,64,80,96
    auto N = 2*Nx;
    Real Jp = 0.041;
        

    auto sites = SpinHalf(N);

    //auto Jps = vector<Real>();
    //auto dimer = vector<Real>();

    for(int i1 = 52; i1 < 160; i1 += 6)
        {
        //
        // Compute ground state using
        // "black box" routine (or have a look at j1j2.h)
        //
        MPS psi = computeGroundState(sites,Jp);
        
        Real G = 0.;
        Real Gzz = 0.;
        //
        // G = B_(i)
        //
        auto i2 = i1 + 2;
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
        G = Ca.real();
        Gzz = Cazz.real();
            
        printfln("b = %.5f, S.S = %.10f, Sz.Sz = %.10f",i1,G,Gzz);
        
        ofstream FileOutput;
        FileOutput.open(format("SS_N=%d_Jp=%d_c.txt",Nx,Jp),ios::out|ios::app);
        FileOutput << fixed << setprecision(10) << Nx << " " << i1 << " " << G << " " << Gzz << endl;
        FileOutput.close();
        //FileOutput.open(format("Correlator_SS_N=%d.txt",N),ios::out|ios::app);
        //FileOutput << fixed << setprecision(10)  << i << " " << G << endl;
        //FileOutput.close();
            
            
        }



    return 0;
    }

