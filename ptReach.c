//
//  ptReach.c
//  
//
//  Created by Javier Martin Blanco on 20/06/2018.
//

#include <stdio.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TMathText.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TVectorD.h>
#include "TLine.h"

// Constants
const double maxLumi = 100000; //mub-1
const double lumipbpb_ABCD = 351*1.049; //mub-1
const double NMB = 2.366003e9*1.049;
const double TAA = 5.607; //mb-1
const double deltaLumi = 500.0;
const double deltapT = 0.005;
const double refLumi = 10000; //mub-1

double normfactor = (NMB*TAA*1e-6); // the 1e-3 factor is because taa is in mb-1 while lumis are in mub-1
double maxPT = 130.0;
double pTbinWidth = 15.0;
double pTlow = 35.0;
double pTup = 50.0;
double deltaRap = 4.8;

TH1D* histoJpsi();
TH1D* histoUpsilon1S();
TH1D* histoUpsilon2S();
TH1D* histoUpsilon3S();

void ptReach(const char* partName = "Jpsi" // partName can be "Jpsi", "Upsilon1S", "Upsilon2S" or "Upsilon3S"
)
{
  TH1D* hxSec(0x0);
  const char* title = "";
  bool isNorm = false;
  bool noLowPTfit = false;
  
  if (!strcmp(partName,"Jpsi")) {hxSec = histoJpsi(); title = "Prompt J/#psi"; isNorm = true; noLowPTfit = true;}
  else if (!strcmp(partName,"Upsilon1S")) {hxSec = histoUpsilon1S(); title = "#Upsilon(1S)";}
  else if (!strcmp(partName,"Upsilon2S")) {hxSec = histoUpsilon2S(); title = "#Upsilon(2S)";}
  else if (!strcmp(partName,"Upsilon3S")) {hxSec = histoUpsilon3S(); title = "#Upsilon(3S)"; /*noLowPTfit = true;*/}
  else {std::cout << "[ERROR] valid options for partName are: Jpsi, Upsilon1S, Upsilon2S or Upsilon3S" << std::endl; return;}
  
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  gStyle->SetOptStat(0);
  gStyle->SetPadTickY(1);
  c1->Range(-16.07516,-7.082789,108.7683,1.671024);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetLogy();
  c1->SetLeftMargin(0.1287625);
  c1->SetRightMargin(0.07023411);
  c1->SetTopMargin(0.07665505);
  c1->SetBottomMargin(0.1236934);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderMode(0);
  
  TH1D* haxis = new TH1D("haxis","",1,0.0,maxPT);
  haxis->GetXaxis()->SetTitle("p_{T} (GeV/#it{c})");
  haxis->GetYaxis()->SetTitle(isNorm?"#bf{#it{#Beta}} #times (1/T_{AA} N_{MB}) dN/dp_{T} (nb / GeV/#it{c})":"dN/dp_{T} (events / GeV/#it{c})");
  haxis->GetYaxis()->SetTitleOffset(1.55);
  haxis->GetXaxis()->SetTitleOffset(1.25);
  haxis->GetXaxis()->CenterTitle(true);
  
  TF1* f1(0x0);
  if (noLowPTfit) f1 = new TF1("f1","[0]/x^[1]",isNorm?6.5:0.0,maxPT);
  else {
    f1 = new TF1("f1","[0]*x/(1+(x*x/([1]*[1])))^[2]",isNorm?6.5:0.0,maxPT);
    f1->SetParameter(0,isNorm?3E6:2000.0);
    f1->SetParameter(1,isNorm?0.5:10.0);
    f1->SetParameter(2,3.0);
  }
    
  hxSec->Fit(f1,"I","",isNorm?6.5:0.0,maxPT);
  //  gxSec->Fit(f1,"I","",6.5,maxPT);
  
  haxis->GetYaxis()->SetRangeUser(f1->Eval(maxPT),isNorm?100.0:100000.0);
  haxis->Draw();
  //  gxSec->Draw("P");
  hxSec->Draw("sameP");
  
  TLatex *tex = new TLatex(0.2,0.77,"#splitline{|y| < 2.4}{0-100 %}");
  tex->SetNDC();
  tex->SetTextSize(0.037);
  tex->SetTextFont(42);
  tex->Draw();
  
  TLatex *tl = new TLatex(0.2,0.85,title);
  tl->SetNDC();
  tl->SetTextFont(42);
  tl->SetTextSize(0.055);
  tl->Draw();
  
  TLatex *t2 = new TLatex(0.70,0.83,"CMS");
  t2->SetNDC();
  t2->SetTextFont(61);
  t2->SetTextSize(0.075);
  t2->Draw();
  
  TLatex *t22 = new TLatex(0.64,0.78,"#it{Projection}");
  t22->SetNDC();
  t22->SetTextFont(42);
  t22->SetTextSize(0.050);
  t22->Draw();
  
  c1->SaveAs(Form("./%s_dNdpT.pdf",partName));
  
  
  ///////////////////////
  // pT bin vs luminosity
  //////////////////////
  
  double Nevents = f1->Integral(pTlow,pTup)*normfactor*deltaRap;
  //  double Nevents = y[n-1]*normfactor*2*ex[n-1]*deltaRap;
  std::cout << Nevents << std::endl;
  
  double Nevents_newLumi = (Nevents/lumipbpb_ABCD)*maxLumi;
  //  std::cout << "Events last bin current lumi = " << Nevents << std::endl;
  //  std::cout << "Events last bin max lumi = " << Nevents_newLumi << std::endl;
  
  
  double newLumi = lumipbpb_ABCD + deltaLumi;
  const int nPoints = int((maxLumi - lumipbpb_ABCD)/deltaLumi);
  double vx[nPoints];
  double vy_min[nPoints];
  double vy_max[nPoints];
  double vex[nPoints];
  double vey[nPoints];
  
  vx[0] = lumipbpb_ABCD;
  vy_max[0] = pTup;
  vy_min[0] = pTlow;
  vex[0] = 0.0;
  vey[0] = 0.0;
  
  for (int i = 1 ; i < nPoints ; i++)
  {
    Nevents_newLumi = 0.0;
    int j = 0;
    double pt_min = maxPT - pTbinWidth;
    double pt_max = maxPT;
    while (Nevents_newLumi < Nevents)
    {
      pt_min -= deltapT*j;
      pt_max -= deltapT*j;
      Nevents_newLumi = (f1->Integral(pt_min,pt_max)*normfactor/lumipbpb_ABCD)*newLumi*deltaRap;
      j++;
      if (pt_min <= pTlow) break;
    }
    
    vx[i] = newLumi;
    vy_max[i] = pt_max;
    vy_min[i] = pt_min;
    vex[i] = 0.0;
    vey[i] = 0.0;
    
    //    std::cout << "Events last bin new lumi 2 = " << Nevents_newLumi << std::endl;
    //    std::cout << "New pT bin = " << pt_min << " - " << pt_max << std::endl;
    
    //    newLumi = newLumi + deltaLumi;
    newLumi = newLumi*2.0;
  }
  
  TGraphErrors* gptVSlumi_max = new TGraphErrors(nPoints,vx,vy_max,vex,vey);
  gptVSlumi_max->SetLineColor(2);
  
  TGraphErrors* gptVSlumi_min = new TGraphErrors(nPoints,vx,vy_min,vex,vey);
  gptVSlumi_min->SetLineColor(2);
  
  TGraphErrors* gptVSlumi_fill = new TGraphErrors(2*nPoints);
  for(int i = 0 ; i < nPoints ; i++)
  {
    gptVSlumi_fill->SetPoint(i,vx[i],vy_max[i]);
    gptVSlumi_fill->SetPoint(nPoints+i,vx[nPoints-i-1],vy_min[nPoints-i-1]);
  }
  
  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  gStyle->SetOptStat(0);
  gStyle->SetPadTickY(1);
  c2->Range(-16.07516,-7.082789,108.7683,1.671024);
//  c2->SetGridy();
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetBorderSize(2);
  //  c2->SetLogy();
  c2->SetLogx();
  c2->SetLeftMargin(0.1287625);
  c2->SetRightMargin(0.07023411);
  c2->SetTopMargin(0.07665505);
  c2->SetBottomMargin(0.1236934);
  c2->SetFrameBorderMode(0);
  c2->SetFrameBorderMode(0);
  
  TH1D* hptVSlumi_axis = new TH1D("hptVSlumi",Form("%0.1f GeV/#it{c} wide bins with constant N_{%s}",pTup-pTlow,title),1,200.0,maxLumi);
  hptVSlumi_axis->GetXaxis()->SetTitle("Lumi (#mub^{-1})");
  hptVSlumi_axis->GetYaxis()->SetTitle("p_{T} (GeV/#it{c})");
  hptVSlumi_axis->GetYaxis()->SetTitleOffset(1.55);
  hptVSlumi_axis->GetXaxis()->SetTitleOffset(1.25);
  hptVSlumi_axis->GetXaxis()->CenterTitle(true);
  hptVSlumi_axis->GetYaxis()->SetRangeUser(pTlow-3,maxPT);
  hptVSlumi_axis->Draw();
  
  gptVSlumi_fill->SetFillStyle(1001);
  gptVSlumi_fill->SetFillColorAlpha(2, 0.2);
  //  gptVSlumi_fill->SetFillColor(2);
  gptVSlumi_fill->Draw("f");
  gptVSlumi_max->Draw("l");
  gptVSlumi_min->Draw("l");
  
  TLatex *tex2 = new TLatex(0.2,0.77,"#splitline{|y| < 2.4}{0-100 %}");
  tex2->SetNDC();
  tex2->SetTextSize(0.037);
  tex2->SetTextFont(42);
  tex2->Draw();
  
  TLatex *tl2 = new TLatex(0.2,0.85,title);
  tl2->SetNDC();
  tl2->SetTextFont(42);
  tl2->SetTextSize(0.055);
  tl2->Draw();
  
  TLatex *t3 = new TLatex(0.70,0.83,"CMS");
  t3->SetNDC();
  t3->SetTextFont(61);
  t3->SetTextSize(0.075);
  t3->Draw();
  
  TLatex *t33 = new TLatex(0.64,0.78,"#it{Projection}");
  t33->SetNDC();
  t33->SetTextFont(42);
  t33->SetTextSize(0.050);
  t33->Draw();
  
  TLine* line = new TLine(refLumi, pTlow-3, refLumi, gptVSlumi_max->Eval(refLumi));
  line->SetLineStyle(2);
  line->SetLineWidth(3);
  line->Draw();
  
  TLine* lineHmax = new TLine(200.0, gptVSlumi_max->Eval(refLumi,0,"S"), refLumi, gptVSlumi_max->Eval(refLumi,0,"S"));
  lineHmax->SetLineStyle(2);
  lineHmax->SetLineWidth(3);
  lineHmax->Draw();
  
  TLine* lineHmin = new TLine(200.0, gptVSlumi_min->Eval(refLumi,0,"S"), refLumi, gptVSlumi_min->Eval(refLumi,0,"S"));
  lineHmin->SetLineStyle(2);
  lineHmin->SetLineWidth(3);
  lineHmin->Draw();
  
  c2->SaveAs(Form("./%s_pTvsLumi.pdf",partName));
}

TH1D* histoJpsi()
{
  normfactor = (NMB*TAA*1e-6); // the 1e-3 factor is because taa is in mb-1 while lumis are in mub-1
  maxPT = 130.0;
  pTbinWidth = 15.0;
  pTlow = 35.0;
  pTup = 50.0;
  deltaRap = 4.8;
  
  const int n = 12;
  const double x[n] = {7.0, 8.0, 9.0, 10.25, 12.0, 14.0, 16.25, 18.75, 22.5, 27.5, 32.5 , 42.5};
  const double ex[n] = {0.5, 0.5, 0.5, 0.75, 1.0, 1.0, 1.25, 1.25, 2.5, 2.5, 2.5, 7.5};
  const double xlow[n+1] = {6.5, 7.5, 8.5, 9.5, 11.0, 13.0, 15.0, 17.5, 20.0, 25.0, 30.0 , 35.0, 50.0};
  
  const double y[n] = {3.624, 1.818, 0.945, 0.476, 0.2021, 0.0889, 0.03721, 0.01616, 0.00688, 0.00258, 0.001017, 0.000241}; // \mathcal{B} \times (1/(N_{MB}T_{AA})) \times dN/dp_{T}, units: nb/GeV
  const double ey_sta[n] = {0.118, 0.050, 0.024, 0.011, 0.0047, 0.0025, 0.0012, 0.00069, 0.00031, 0.00018, 0.000098, 0.000029};
  const double ey_sys[n] = {0.281, 0.111, 0.049, 0.022, 0.0090, 0.0040, 0.0017, 0.00095, 0.00034, 0.00018, 0.000127, 0.000023};
  double ey[n] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (int i = 0 ; i < n ; i++)
  {
    ey[i] = sqrt(ey_sta[i]*ey_sta[i] + ey_sys[i]*ey_sys[i]);
  }
  
//  TGraphErrors* gxSec = new TGraphErrors(n,x,y,ex,ey);
  
  TH1D* hxSec = new TH1D("hxSec","",n,xlow);
  for (int i = 1 ; i <= n ; i++)
  {
    hxSec->SetBinContent(i,y[i-1]);
    hxSec->SetBinError(i,ey[i-1]);
  }
  
  return hxSec;
}

TH1D* histoUpsilon1S()
{
  normfactor = 1.0;
  maxPT = 100.0;
  pTbinWidth = 18.0;
  pTlow = 12.0;
  pTup = 30.0;
  deltaRap = 1.0;
  
  const int n = 6;
  const double x[n] = {1.0, 3.0, 5.0, 7.0, 11.0, 21.0};
  const double ex[n] = {1.0, 1.0, 1.0, 2.0, 2.0, 9.0};
  const double xlow[n+1] = {0.0, 2.0, 4.0, 6.0, 9.0, 12.0, 30.0};
  
  double y[n] = {4057.33, 9233.16, 9029.76, 7410.23, 3220.13, 2315.88}; // N_{Upsilon}
  double ey_sta[n] = {346.697, 686.329, 733.311, 493.686, 253.713, 131.588};
  double ey_sys[n] = {1535.06, 693.973, 618.72, 417.505, 330.711, 170.358};
  double ey[n] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (int i = 0 ; i < n ; i++)
  {
    ey[i] = sqrt((ey_sta[i]*ey_sta[i] + ey_sys[i]*ey_sys[i])/(4*ex[i]*ex[i]));
    y[i] = y[i]/(2*ex[i]); // To convert it to dN/dpT
  }
  
//  TGraphErrors* gxSec = new TGraphErrors(n,x,y,ex,ey);
  
  TH1D* hxSec = new TH1D("hxSec","",n,xlow);
  for (int i = 1 ; i <= n ; i++)
  {
    hxSec->SetBinContent(i,y[i-1]);
    hxSec->SetBinError(i,ey[i-1]);
  }
  
  return hxSec;
}

TH1D* histoUpsilon2S()
{
  normfactor = 1.0;
  maxPT = 60.0;
  pTbinWidth = 21.0;
  pTlow = 9.0;
  pTup = 30.0;
  deltaRap = 1.0;
  
  const int n = 3;
  const double x[n] = {2.0, 6.5, 19.5};
  const double ex[n] = {2.0, 2.5, 10.5};
  const double xlow[n+1] = {0.0, 4.0, 9.0, 30.0};
  
  double y[n] = {702.497, 1368.4, 555.391}; // N_{Upsilon}
  double ey_sta[n] = {323.249, 399.333, 316.797};
  double ey_sys[n] = {770.364, 142.122, 68.6075};
  double ey[n] = {0.0, 0.0, 0.0};
  for (int i = 0 ; i < n ; i++)
  {
    ey[i] = sqrt((ey_sta[i]*ey_sta[i] + ey_sys[i]*ey_sys[i])/(4*ex[i]*ex[i]));
    y[i] = y[i]/(2*ex[i]); // To convert it to dN/dpT
  }
  
  //  TGraphErrors* gxSec = new TGraphErrors(n,x,y,ex,ey);
  
  TH1D* hxSec = new TH1D("hxSec","",n,xlow);
  for (int i = 1 ; i <= n ; i++)
  {
    hxSec->SetBinContent(i,y[i-1]);
    hxSec->SetBinError(i,ey[i-1]);
  }
  
  return hxSec;
}

TH1D* histoUpsilon3S()
{
  normfactor = 1.0;
  maxPT = 60.0;
  pTbinWidth = 24.0;
  pTlow = 6.0;
  pTup = 30.0;
  deltaRap = 1.0;
  
  const int n = 2;
  const double x[n] = {3.0, 18.0};
  const double ex[n] = {3.0, 12.0};
  const double xlow[n+1] = {0.0, 6.0, 30.0};
  
  double y[n] = {144.191, 41.3866}; // N_{Upsilon}
  double ey_sta[n] = {374.275, 371.695};
  double ey_sys[n] = {185.691, 226.309};
  double ey[n] = {0.0, 0.0};
  for (int i = 0 ; i < n ; i++)
  {
    ey[i] = sqrt((ey_sta[i]*ey_sta[i] + ey_sys[i]*ey_sys[i])/(4*ex[i]*ex[i]));
    y[i] = y[i]/(2*ex[i]); // To convert it to dN/dpT
  }
  
  //  TGraphErrors* gxSec = new TGraphErrors(n,x,y,ex,ey);
  
  TH1D* hxSec = new TH1D("hxSec","",n,xlow);
  for (int i = 1 ; i <= n ; i++)
  {
    hxSec->SetBinContent(i,y[i-1]);
    hxSec->SetBinError(i,ey[i-1]);
  }
  
  return hxSec;
}
