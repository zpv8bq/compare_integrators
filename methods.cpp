// This code may be compiled to make a stand alone exe
// or it can be run from the ROOT command line as:

#include "TApplication.h"
#include "Math/QuasiRandom.h"
#include "TMath.h"
#include "TDatime.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TString.h"
#include "TROOT.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TVectorD.h"
#include "TLegend.h"
#include "TRandom3.h"
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;
using namespace ROOT::Math;


class SobolN {
public:
  SobolN(int ndim, long skip){
    r1=new QuasiRandomSobol(ndim);
    r1->Skip(skip);
  }
  void RndmArray(double *array) {
    r1->Next(array);
  }
  QuasiRandomSobol *r1;
};


// test for point in hypersphere of radius =1
int TestPoint(Double_t *x){
  double val=0;
  for (int i=0; i<5; i++) val+=x[i]*x[i];
  if (val<1.0) return 1;
  else return 0;
}

// chop up 5d space
long TestGrid(long ndiv){
  static TRandom3 tr3(0);
  long hits=0;
  double step=1.0/ndiv;
  double offset=1.0/(ndiv+1);
  double x[5];
  for (int i=0; i<ndiv; i++) {
    x[0]=offset+step*i;
    for (int j=0; j<ndiv; j++) {
      x[1]=offset+step*j;
      for (int k=0; k<ndiv; k++) {
	x[2]=offset+step*k;
	for (int l=0; l<ndiv; l++) {
	  x[3]=offset+step*l;
	  for (int m=0; m<ndiv; m++) {
	    x[4]=offset+step*m;
	    double val=x[4]*x[4];
	    for (int n=0; n<4; n++) val+=x[n]*x[n];
	    if (val<1) hits++;
	  }
	}
      }
    }
  }
  return hits;
}

// chop up 5d space
long TestSubR(long ndiv){
  static TRandom3 tr3(0);
  long hits=0;
  double step=1.0/ndiv;
  double x[5];
  for (int i=0; i<ndiv; i++) {
    x[0]=step*i;
    x[0]+=tr3.Uniform(step);
    for (int j=0; j<ndiv; j++) {
      x[1]=step*j;
      x[1]+=tr3.Uniform(step);
      for (int k=0; k<ndiv; k++) {
	x[2]=step*k;
	x[2]+=tr3.Uniform(step);
	for (int l=0; l<ndiv; l++) {
	  x[3]=step*l;
	  x[3]+=tr3.Uniform(step);
	  for (int m=0; m<ndiv; m++) {
	    x[4]=step*m;
	    x[4]+=tr3.Uniform(step);
	    double val=x[4]*x[4];
	    for (int n=0; n<4; n++) val+=x[n]*x[n];
	    if (val<1) hits++;
	  }
	}
      }
    }
  }
  return hits;
}


// integrations of 5dimensional sphere
void Integral5d(long ne6=20, long skip=0){
  gStyle->SetOptStat(0);
  // setup generators
  TRandom3 tr3(0);  
  SobolN sob(5,skip);
  long hits[4];   // keep track of hits for different sets of samples
  TStopwatch ts;
  TDatime td;
  TFile *tf=new TFile(TString::Format("run_%d_%ld.root",td.GetTime(),skip),"recreate");
  TGraph *graphs[4];  // convergence graphs
  for (int i=0; i<4; i++) {
    graphs[i]=new TGraph();
    graphs[i]->SetName(TString::Format("tg%d",i));
    graphs[i]->SetLineWidth(3);
    graphs[i]->SetLineColor(i+1);
    hits[i]=0;
  }
  long NTRIALS = 1000*1000*ne6;
  const long MINTRIALS=100;
  long showat=100;
  const double SHOWMULT=TMath::Power(10.,1./30);  // 10^1/n, plot ~ n pts/decade
  const double PISQ=9.86960440108936;      // pi^2
  const double V5D = 8*PISQ/15/32;         // Volume of positive hyperquadrant
  const double VBOX = 1;

  cout << "Study Integrals using " << NTRIALS << " points" << endl;

  
  Double_t x[5];
  for (long n=MINTRIALS; n<=NTRIALS; n++){

    tr3.RndmArray(5,x);
    hits[0]+=TestPoint(x);
    sob.RndmArray(x);
    hits[1]+=TestPoint(x);
    long nGrid;
    if (n==showat){
      cout << "Trial: " << n << " of: " << NTRIALS << endl;
      // calc the fixed and subrandom grids here
      if (n<1000*1000*100){
	long nDiv=TMath::Power(n,0.2);   // divisions per dimension
	nGrid=nDiv*nDiv*nDiv*nDiv*nDiv;  // total points sampled
	hits[2]=TestGrid(nDiv);
	hits[3]=TestSubR(nDiv);
      }
      showat*=SHOWMULT;

      for (int j=0;j<4;j++) {
	double p;
	if (j<2) p = (double) hits[j]/n;
	else p = (double) hits[j]/nGrid;
	double V = VBOX * p;
	if (j==0) cout << V << endl;
	double err = TMath::Abs(V-V5D)/V5D;
	if (n<1000*1000*100 || j<2)
	  graphs[j]->SetPoint(graphs[j]->GetN(),sqrt(n),err);
      }
    }
  }
  

  TCanvas *tc=new TCanvas("c1","Convergence tests");
  tc->SetLogy();
  //tc->SetLogx();
  TH2F *hax=new TH2F("hax",";sqrt{Number of samples};Relative error",
		     10,MINTRIALS-1,sqrt(NTRIALS),10,1e-8,2);
  hax->Draw();
  for (int i=0; i<4; i++) {
    graphs[i]->Draw("LP");
  }

  TLegend *tl=new TLegend(0.1,0.1,0.5,0.4);
  tl->AddEntry(graphs[2],"Fixed grid","l");
  tl->AddEntry(graphs[0],"Pseudorandom","l");
  tl->AddEntry(graphs[3],"Subrandom","l");
  tl->AddEntry(graphs[1],"Sobol","l");
  tl->Draw();
  tc->Update();
  
  tc->Print("methods.pdf");
  tc->Print("methods.png");
  for (int i=0; i<4; i++) graphs[i]->Write();
  tf->Close();
  ts.Stop();
  cout << "Wall time: " << ts.RealTime() << endl;
  cout << "CPU time: " << ts.CpuTime() << endl;
}

int main(int argc, char **argv) {
  int ne6=100;  // 100 million samples
  long skip=0;

  // This allows you to view ROOT-based graphics in your C++ program
  // If you don't want view graphics (eg just read/process/write data files), 
  // this can be ignored
  TApplication theApp("App", &argc, argv);

  if (argc>1) ne6=atoi(argv[1]);
  if (argc>2) skip=atol(argv[1]);
  
  Integral5d(ne6,skip);

  // view graphics in ROOT if we are in an interactive session
  if (!gROOT->IsBatch()) {
      cout << "To exit, quit ROOT from the File menu of the plot (or use control-C)" << endl;
    theApp.SetIdleTimer(30,".q"); // set up a failsafe timer to end the program
    theApp.Run(true);
  }
  return 0;
}
  
