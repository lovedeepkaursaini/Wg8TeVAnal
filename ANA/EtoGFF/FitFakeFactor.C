#include "iostream"
#include <cmath>
#include "sstream"
#include "TGraph.h"
#include "RooFormulaVar.h"
double GetErr(double A, double errA, double B, double errB)
{
  return err = (A/B)*sqrt(((errA/A)*(errA/A))+((errB/B)*(errB/B)));
}

RooRealVar FitMeg(TH1D* Z_mass,TString saveas)
{
  
  gSystem->Load("libRooFit");
  
  using namespace RooFit;
  
  //TH1F* Z_mass =hName;// (TH1F*)f.Get(hName);
  //Z_mass->Rebin(1.5);
  // Z_mass=(TH1F*)hist->Clone();
  double hmin = 65.;//Z_mass->GetXaxis()->GetXmin();
  double hmax = 200.;//Z_mass->GetXaxis()->GetXmax();
  double nMax = Z_mass->Integral();
  
  // Declare observable x
  RooRealVar xtmp("xtmp","xtmp",40,200) ;
  RooRealVar x("x","x",hmin,hmax) ;
  RooDataHist dh("dh","dh",x,Import(*Z_mass)) ;
  RooDataHist x1("x1","x1",xtmp,Import(*Z_mass)) ;


  //Construct the signal P.D.F., a gaussian function
  RooRealVar mean_bw("mean_bw","mean of bw",92,75.,100.);
  RooRealVar mean_cb("mean_cb","mean of cb",0,-5.0,10);
  RooRealVar sigma_bw("sigma_bw","width of bw",2.5,0,5.);
  RooRealVar sigma_cb("sigma_cb","width of cb",2.0, 0.00001, 8.0);//this range worked for overall pt range
  // RooRealVar sigma_cb("sigma_cb","width of cb",5.0, -20,20.);
  RooRealVar n("n","", 6.7489,1.0,15.0);
  // RooRealVar alpha("alpha","", 1.0,0.0002,5.0);//this range worked for overall pt range
  RooRealVar alpha("alpha","", 1.0,-50.,50);

  //RooGaussian SigPDF("gauss","Signal PDF",x,mean,sigma);  
  RooBreitWigner bwPDF("BreitWigner","BreitWigner",x,mean_bw,sigma_bw);
  // RooBreitWigner SigPDF("BreitWigner","BreitWigner",x,mean,sigma);
  //  RooCBShape SigPDF("cball", "crystal ball", x, mean, sigma, alpha, n);
  RooCBShape cbPDF("cball", "crystal ball", x, mean_cb, sigma_cb, alpha, n);
  RooFFTConvPdf SigPDF("bwxCryBall", "FFT Conv CryBall and BW", x, bwPDF, cbPDF); 
  //Now define the background P.D.F, a simple exponential
  RooRealVar tau("tau","exponential function parameter",0,-10.,10.);//this range worked for overall pt range
  //RooRealVar tau("tau","exponential function parameter",-0.001,-5.0,5.0);
  RooExponential exp("exponential","Background PDF",x,tau);
   RooFormulaVar eff("eff","0.5*(TMath::Erf((x-1)/0.5)+1)",x) ;
   // RooEffProd BkgPDF("BkgPDF","exp. * Err.fn",exp,eff) ;
  RooAddPdf BkgPDF("BkgPDF","exp. * Err.fn",exp,eff) ;

  RooRealVar a0("a0","a0",0.,-5.,5.);
  RooRealVar a1("a1","a1",0.,-1.,1.);
  RooRealVar a2("a2","a2",0.,-1.,1.);
  //RooChebychev BkgPDF("bkg1","bkg1",x,RooArgList(a0,a1));

  //Now construct the total PDF. We need to define the number of signal and background events in the model
  RooRealVar Nsig("Nsig","Number of signal events",0.5*nMax,0.,nMax);
  RooRealVar Nbkg("Nbkg","Number of background events",0.2*nMax,0.,nMax);

  RooAddPdf PDFtot("PDFtot","PDFtot",RooArgList(SigPDF,BkgPDF),RooArgList(Nsig,Nbkg));
  //Now generate a sample with the total PDF
  // RooDataSet *data = PDFtot.generate(RooArgSet(mass),NumEvents(Nsig.getVal()+Nbkg.getVal()),Extended(1));

  //Now fit the PDF to the data
  //For low statistics, fix the mean and the width
  // mean.setConstant(kTRUE);
  // sigma.setConstant(kTRUE);
  PDFtot.fitTo(dh);//ML fit is default

  // Print values of mean and sigma (that now reflect fitted values and errors, unless you fixed them)
  mean_bw.Print();
  // sigma.Print();
  //Now plot the data and the fitted PDF
  RooPlot* frame = x.frame() ;
  x1.plotOn(frame);//MarkerColor(2),MarkerSize(0.9),MarkerStyle(21));  //this will show histogram data points on canvas 
  dh.plotOn(frame);//MarkerColor(2),MarkerSize(0.9),MarkerStyle(21));  //this will show histogram data points on canvas 
 //dh.statOn(frame);  //this will display hist stat on canvas

  PDFtot.plotOn(frame,LineColor(4));//this will show fit overlay on canvas 
  PDFtot.paramOn(frame,Layout(0.70,0.98,0.8)); //this will display the fit parameters on canvas
  //One can also plot the single components of the total PDF, like the background component
  PDFtot.plotOn(frame,Components(BkgPDF),LineStyle(kDashed),LineColor(kRed));

  //Actually plot the result
  TCanvas c1;
  c1.cd();
  c1.SetLogy();
  frame->Draw();
  //double chi2 = RooChi2Var("chi2","chi2",PDFtot,dh);
  double chi2 = frame->chiSquare('PDFtot');
  cout<<chi2<<endl;
  c1.SaveAs("png/"+saveas+".png");
  //c1.SaveAs(Z_mass);

  cout<<Nsig.getVal()<<endl;
  return Nsig;
}

void ffRate(TString File) {
  gStyle->SetOptFit();
  TCanvas *c1 = new TCanvas("c1","multigraph",200,10,700,500);
  c1->SetGrid();
  
  // draw a frame to define the range
  TMultiGraph *mg = new TMultiGraph();
  
  // create first graph
  const Int_t n1 = 8;//
  
  Double_t x1[]  = {0.0,0.1,0.5,1.0,1.5,2.1,2.2,2.4,2.5};
  Double_t y1[] = {0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t ex1[] = {0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t ey1[] = {0,0,0,0,0,0,0,0,0,0,0,0};

  Double_t y2[]  = {0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t ex2[] = {0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t ey2[] = {0,0,0,0,0,0,0,0,0,0,0,0};

  Double_t y3[]  = {0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t ex3[] = {0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t ey3[] = {0,0,0,0,0,0,0,0,0,0,0,0};

  Double_t y4[]  = {0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t ex4[] = {0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t ey4[] = {0,0,0,0,0,0,0,0,0,0,0,0};

  TFile *file = new TFile(File);
  for(int i=0;i<n1;i++){
    std::stringstream ii;
    ii<<i+1;
    TH1D* hpass1 = file->Get(("passPSV_Meg_pt1_eta"+ii.str()).c_str());
    double passPSV1=FitMeg(hpass1,("passPSV_Meg_pt1_eta"+ii.str()).c_str()).getVal();
    double ErrpassPSV1=FitMeg(hpass1,("passPSV_Meg_pt1_eta"+ii.str()).c_str()).getError();
    TH1D* hfail1 = file->Get(("failPSV_Meg_pt1_eta"+ii.str()).c_str());
    double failPSV1   =FitMeg(hfail1,("failPSV_Meg_pt1_eta"+ii.str()).c_str()).getVal();
    double ErrfailPSV1=FitMeg(hfail1,("failPSV_Meg_pt1_eta"+ii.str()).c_str()).getError();
    y1[i] = passPSV1/failPSV1;
    ey1[i] = GetErr(passPSV1,ErrpassPSV1,failPSV1,ErrfailPSV1);
    // cout<<passPSV1<<'\t'<<ErrpassPSV1<<'\t'<<failPSV1<<'\t'<<ErrfailPSV1<<endl;
    // cout<<y1[i]<<": "<<ey1[i]<<endl;
    TH1D* hpass2 = file->Get(("passPSV_Meg_pt2_eta"+ii.str()).c_str());
    double passPSV2=FitMeg(hpass2,("passPSV_Meg_pt2_eta"+ii.str()).c_str()).getVal();
    double ErrpassPSV2=FitMeg(hpass2,("passPSV_Meg_pt2_eta"+ii.str()).c_str()).getError();
    TH1D* hfail2 = file->Get(("failPSV_Meg_pt2_eta"+ii.str()).c_str());
    double failPSV2=FitMeg(hfail2,("failPSV_Meg_pt2_eta"+ii.str()).c_str()).getVal();
    double ErrfailPSV2=FitMeg(hfail2,("failPSV_Meg_pt2_eta"+ii.str()).c_str()).getError();
    y2[i] = passPSV2/failPSV2;
    ey2[i] = GetErr(passPSV2,ErrpassPSV2,failPSV2,ErrfailPSV2);
    cout<<y2[i]<<endl;

    TH1D* hpass3 = file->Get(("passPSV_Meg_pt3_eta"+ii.str()).c_str());
    double passPSV3=FitMeg(hpass3,("passPSV_Meg_pt3_eta"+ii.str()).c_str()).getVal();
    double ErrpassPSV3=FitMeg(hpass3,("passPSV_Meg_pt3_eta"+ii.str()).c_str()).getError();
    TH1D* hfail3 = file->Get(("failPSV_Meg_pt3_eta"+ii.str()).c_str());
    double failPSV3=FitMeg(hfail3,("failPSV_Meg_pt3_eta"+ii.str()).c_str()).getVal();
    double ErrfailPSV3=FitMeg(hfail3,("failPSV_Meg_pt3_eta"+ii.str()).c_str()).getError();
    y3[i] = passPSV3/failPSV3;
    ey3[i] = GetErr(passPSV3,ErrpassPSV3,failPSV3,ErrfailPSV3);
    cout<<y3[i]<<endl;

    TH1D* hpass4 = file->Get(("passPSV_Meg_pt4_eta"+ii.str()).c_str());
    double passPSV4=FitMeg(hpass4,("passPSV_Meg_pt4_eta"+ii.str()).c_str()).getVal();
    double ErrpassPSV4=FitMeg(hpass4,("passPSV_Meg_pt4_eta"+ii.str()).c_str()).getError();
    TH1D* hfail4 = file->Get(("failPSV_Meg_pt4_eta"+ii.str()).c_str());
    double failPSV4=FitMeg(hfail4,("failPSV_Meg_pt4_eta"+ii.str()).c_str()).getVal();
    double ErrfailPSV4=FitMeg(hfail4,("failPSV_Meg_pt4_eta"+ii.str()).c_str()).getError();
    y4[i] = passPSV4/failPSV4;
    ey4[i] = GetErr(passPSV4,ErrpassPSV4,failPSV4,ErrfailPSV4);
    cout<<y4[i]<<endl;


  }
  // create second graph
  TGraphErrors *gr1 = new TGraphErrors(n1,x1,y1,ex1,ey1);
  gr1->SetMarkerColor(kBlack);
  gr1->SetMarkerStyle(21);
  // gr1->Fit("pol6","q");
  
  
  
  TGraphErrors *gr2 = new TGraphErrors(n1,x1,y2,ex2,ey2);
  gr2->SetMarkerColor(kRed);
  gr2->SetMarkerStyle(20);

  TGraphErrors *gr3 = new TGraphErrors(n1,x1,y3,ex3,ey3);
  gr3->SetMarkerColor(kBlue);
  gr3->SetMarkerStyle(22);

  TGraphErrors *gr4 = new TGraphErrors(n1,x1,y4,ex4,ey4);
  gr4->SetMarkerColor(kMagenta);
  gr4->SetMarkerStyle(22);


  mg->Add(gr4);
  mg->Add(gr3);
  mg->Add(gr2);  
  mg->Add(gr1);
  mg->Draw("ap");
  
   //force drawing of canvas to generate the fit TPaveStats
  c1->Update();
   mg->GetYaxis()->SetTitle("fake factor");
  mg->GetXaxis()->SetTitle("|#eta(#gamma)|");
  //mg->GetXaxis()->SetRangeUser(38,200);
  //mg->GetYaxis()->SetRangeUser(0.5,1.2);
  //force drawing of canvas to generate the fit TPaveStats
  c1->Update();
   TLegend *leg = new TLegend(0.5,0.67,0.88,0.88,NULL,"brNDC");
   leg->AddEntry(gr1,"15 < pT(#gamma) < 25","p");
   leg->AddEntry(gr2,"25 < pT(#gamma) < 40","p");
   leg->AddEntry(gr3,"40 < pT(#gamma) < 70","p");
   leg->AddEntry(gr4,"70 < pT(#gamma) < 1000","p");
   // leg->AddEntry(gr3,"Intrinsic","p");
   leg->Draw();
   c1->Modified();
}





