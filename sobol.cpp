//+______________________________________________________________________________
// Example of generating quasi-random numbers

#include "TApplication.h"
#include "Math/QuasiRandom.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TGClient.h"
#include "TStyle.h"
#include <cstdlib>
#include <iostream>

using namespace std;
using namespace ROOT::Math;


// This helper class is used to instantiate the sobel sequence generator for a given dimension
// skip may be used to start at some offset from the beginning of the sequence

class SobolN {
public:
  SobolN(int ndim, unsigned long skip){
    r1=new QuasiRandomSobol(ndim);
    r1->Skip(skip);
  }
  void RndmArray(double *array) {
    r1->Next(array);
  }
private:
  QuasiRandomSobol *r1;
};



// run as
// ./sobol
// or to adjust the numebr of points and offset in the sobel sequence
// eg ./sobol 2000 10000

int main(int argc, char **argv) {

  unsigned long npoints=5000;
  unsigned long nskip=0;
  
  // This allows you to view ROOT-based graphics in your C++ program
  // If you don't want view graphics (eg just read/process/write data files), 
  // this can be ignored
  TApplication theApp("App", &argc, argv);

  if (argc>1) npoints=atol(argv[1]);
  if (argc>2) nskip=atol(argv[2]);
  
  gStyle->SetOptStat(0);

  TGraph *gr = new TGraph();
  TGraph *gs = new TGraph();

  TRandom3 *tr3 = new TRandom3(0);
  
  // quasi random numbers need to be created giving the dimension of the sequence
  // in this case we generate a 2-d sequence
  const int ndim=2;
  SobolN rsobel(ndim,nskip);

  TH2F *h = new TH2F("h","",10,0,1,10,-1,1);
   
  // generate n random number pairs
  double *x = new double[ndim];

  for (int i = 0; i < npoints; ++i){
    tr3->RndmArray(ndim,x);
    gr->SetPoint(i,x[0],x[1]);
    rsobel.RndmArray(x);
    gs->SetPoint(i,x[0],x[1]);
  }
   
  UInt_t dw = gClient->GetDisplayWidth();
  UInt_t dh = gClient->GetDisplayHeight();

  TCanvas * c1 = new TCanvas("c1","Random sequence",dw/2,dh/3);
  c1->Divide(2,1);


  auto hr=new TH2F("hr","PseudoRandom;x1;x2",10,0,1,10,0,1);
  auto hs=new TH2F("hs","Sobol;x1;x2",10,0,1,10,0,1);
  c1->cd(1);
  hr->Draw();
  gr->Draw("p");
  
  c1->cd(2);
  hs->Draw();
  gs->Draw("p");

  c1->SaveAs("random.pdf");
  
  // view graphics in ROOT if we are in an interactive session
  if (!gROOT->IsBatch()) {
    cout << "To exit, quit ROOT from the File menu of the plot (or use control-C)" << endl;
    theApp.SetIdleTimer(30,".q"); // set up a failsafe timer to end the program
    theApp.Run(true);
  }

  delete [] x;
  return 0;
}
