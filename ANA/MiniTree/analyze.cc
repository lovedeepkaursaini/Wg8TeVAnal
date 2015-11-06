#define analyze_cxx
#include "analyze.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include <TMath.h>

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include "photon2012.cc"
#include <fstream>
//ofstream fout("output.txt");

using namespace std;

int   pho_ID = 1;
float dr_pho_ele = 0.7;

// Tree variables
int trigger;
int nvtx;
float IDWeight;

//electron
float ele_pt;
float ele_phi;
float ele_eta;

// photon
int iPho1=-999;
float genPhoPt;
float pho_pt;
float pho_phi;
float pho_eta;
float pho_SCEta;
float pho_Sihih;
float pho_ChHadIso;
float pho_SCRChIso;
float pho_PFPhoIso;
float pho_SCRPhoIso;
float pho_RandConeChIso;
int   mc_truth;
int   PhoSihih_pass;
int   PhoIso_pass;
int   PhoChHadIso_pass;
float PFMET=-9;
float PFMETPhi=-9;
float MatchRecoEleRecoJet  =-9;
float MatchRecoEle2RecoJet =-9;
float MatchRecoPhoRecoJet  =-9;

float MatchRecoEleGenEle=-9;
float MatchRecoPhoGenEle=-9;
float RecoJetPt_drRjRp0p5=-9;

int genPartPdgId_MatchRecoEleGenEle = -999;
int genPartPdgId_MatchRecoPhoGenEle = -999;
int genPartMomPdgId_MatchRecoEleGenEle = -999;
int genPartMomPdgId_MatchRecoPhoGenEle = -999;
float genPartPt_MatchRecoPhoGenPho = -999;
float genPartEta_MatchRecoPhoGenPho = -999;
float genPartPt_MatchRecoEleGenPho = -999;
float genPartEta_MatchRecoEleGenPho = -999;
float genPartPt_MatchRecoPhoGenEle = -999;
float genPartEta_MatchRecoPhoGenEle = -999;
float genPartPt_MatchRecoEleGenEle = -999;
float genPartEta_MatchRecoEleGenEle = -999;

float MatchRecoPhoGenPho=-9;
float MatchRecoEleGenPho=-9;

float InvMass = -9;
float W_MT=-9;
float evt_weight = 1;
int PhohasPixelSeed = -9.;
int event=0;
int ilep1=-999;

float lepton1Pt=-99;
float lepton1Eta=-99;
float lepton1Phi=-99;
float lepton1SCEta=-99;
float lepton1SCPhi=-99;
float lepton1Sieie=-99;
float lepton1HoE = -99;
float lepton1RelIso=-99;

int lepton2VetoID = -99;
float lepton2Pt=-99;
float lepton2Eta=-99;
float lepton2Phi=-99;
float lepton2Sieie=-99;
float lepton2HoE = -99;
float lepton2RelIso=-99;


int run=-999 ;

bool analyze::HasMatchingGSFelectron(vector<int> ielectrons, int ipho)
{
  for (int iele=0; iele<ielectrons.size(); iele++){
    if (elePt->at(iele)<2) break;
    float delR = dR(phoSCEta->at(ipho),phoSCPhi->at(ipho),eleSCEta->at(iele),eleSCPhi->at(iele));
    if (delR<0.02) return 1;
  }
  return 0;
}
//https://twiki.cern.ch/twiki/bin/view/Main/EGammaScaleFactors2012#2012_8_TeV_Jan22_Re_recoed_data
float analyze::SF_TightEID(float pt, float eta) {
  if ( pt <= 40){
    if (fabs(eta) <= 0.8)                       { return 0.978;}
    if (fabs(eta) > 0.8 && fabs(eta) <= 1.4442)  { return 0.958;}
    if (fabs(eta) > 1.566 && fabs(eta) <= 2.0)  { return 0.909;}
    if (fabs(eta) > 2.0 )    { return 0.987;}
  }
  if ( pt > 40 && pt <= 50){
    if (fabs(eta) <= 0.8)                       { return 0.981;}
    if (fabs(eta) > 0.8 && fabs(eta) <= 1.4442)  { return 0.969;}
    if (fabs(eta) > 1.566 && fabs(eta) <= 2.0)  { return 0.942;}
    if (fabs(eta) > 2.0 )    { return 0.991;}
  }
  if ( pt > 50 ) {
    if (fabs(eta) <= 0.8)                       { return 0.982;}
    if (fabs(eta) > 0.8 && fabs(eta) <= 1.4442)  { return 0.969;}
    if (fabs(eta) > 1.566 && fabs(eta) <= 2.0)  { return 0.957;}
    if (fabs(eta) > 2.0 )    { return 0.999;}
  }
}

float analyze::SF_MediumEID(float pt, float eta) {
  if ( pt <= 40){
    if (fabs(eta) <= 0.8)                       { return 1.002;}
    if (fabs(eta) > 0.8 && fabs(eta) <= 1.4442)  { return 0.980;}
    if (fabs(eta) > 1.566 && fabs(eta) <= 2.0)  { return 0.967;}
    if (fabs(eta) > 2.0 )    { return 1.021;}
  }
  if ( pt > 40 && pt <= 50){
    if (fabs(eta) <= 0.8)                       { return 1.005;}
    if (fabs(eta) > 0.8 && fabs(eta) <= 1.4442)  { return 0.988;}
    if (fabs(eta) > 1.566 && fabs(eta) <= 2.0)  { return 0.992;}
    if (fabs(eta) > 2.0 )    { return 1.019;}
  }
  if ( pt > 50 ) {
    if (fabs(eta) <= 0.8)                       { return 1.004;}
    if (fabs(eta) > 0.8 && fabs(eta) <= 1.4442)  { return 0.988;}
    if (fabs(eta) > 1.566 && fabs(eta) <= 2.0)  { return 1.000;}
    if (fabs(eta) > 2.0 )    { return 1.022;}
  }
}

float analyze::SF_LooseEID(float pt, float eta) {
  if ( pt <= 40){
    if (fabs(eta) <= 0.8)                       { return 1.002;}
    if (fabs(eta) > 0.8 && fabs(eta) <= 1.4442)  { return 0.985;}
    if (fabs(eta) > 1.566 && fabs(eta) <= 2.0)  { return 0.980;}
    if (fabs(eta) > 2.0 )    { return 1.019;}
  }
  if ( pt > 40 && pt <= 50){
    if (fabs(eta) <= 0.8)                       { return 1.005;}
    if (fabs(eta) > 0.8 && fabs(eta) <= 1.4442)  { return 0.989;}
    if (fabs(eta) > 1.566 && fabs(eta) <= 2.0)  { return 0.999;}
    if (fabs(eta) > 2.0 )    { return 1.019;}
  }
  if ( pt > 50 ) {
    if (fabs(eta) <= 0.8)                       { return 1.005;}
    if (fabs(eta) > 0.8 && fabs(eta) <= 1.4442)  { return 0.989;}
    if (fabs(eta) > 1.566 && fabs(eta) <= 2.0)  { return 1.004;}
    if (fabs(eta) > 2.0 )    { return 1.023;}
  }
}

// photon medium ID SF

// https://indico.cern.ch/event/305105/contribution/3/material/slides/0.pdf (p-5 * p-8-Medium-colmn)

float analyze::SF_mediumPhoID_EleVeto(float pt, float eta){
  if( pt > 15 && pt <= 20){
    if (fabs(eta) <= 0.8)                     { return 0.9424;}
    if (fabs(eta) > 0.8 && fabs(eta) <= 1.4442){ return 0.9879;}
    if (fabs(eta) > 1.566 && fabs(eta) <= 2.0){ return 0.9929;}
    if (fabs(eta) > 2.0 && fabs(eta) <= 2.5)  { return 1.0084;}
  }
  if( pt > 20 && pt <= 30 ){
    if (fabs(eta) <= 0.8)                     { return 0.9592;}
    if (fabs(eta) > 0.8 && fabs(eta) <= 1.4442){ return 0.9682;}
    if (fabs(eta) > 1.566 && fabs(eta) <= 2.0){ return 0.9841;}
    if (fabs(eta) > 2.0 && fabs(eta) <= 2.5)  { return 1.0052;}
  }
  if( pt > 30 && pt <= 40 ){
    if (fabs(eta) <= 0.8)                     { return 0.9709;}
    if (fabs(eta) > 0.8 && fabs(eta) <= 1.4442){ return 0.9722;}
    if (fabs(eta) > 1.566 && fabs(eta) <= 2.0){ return 0.9851;}
    if (fabs(eta) > 2.0 && fabs(eta) <= 2.5)  { return 0.9974;}
  }
  if( pt > 40 && pt <= 50 ){
    if (fabs(eta) <= 0.8)                     { return 0.9727;}
    if (fabs(eta) > 0.8 && fabs(eta) <= 1.4442){ return 0.9762;}
    if (fabs(eta) > 1.566 && fabs(eta) <= 2.0){ return 0.9788;}
    if (fabs(eta) > 2.0 && fabs(eta) <= 2.5)  { return 0.9893;}
  }
  if( pt > 50 ) {
    if (fabs(eta) <= 0.8)                     { return 0.9732;}
    if (fabs(eta) > 0.8 && fabs(eta) <= 1.4442){ return 0.9767;}
    if (fabs(eta) > 1.566 && fabs(eta) <= 2.0){ return 1.0051;}
    if (fabs(eta) > 2.0 && fabs(eta) <= 2.5)  { return 1.0154;}
  }
}

// https://indico.cern.ch/event/305105/contribution/3/material/slides/0.pdf (p-5) medium photon ID only
float analyze::SF_mediumPhoID(float pt, float eta){
  if( pt > 15 && pt <= 20){
    if (fabs(eta) <= 0.8)                     { return 0.9462;}
    if (fabs(eta) > 0.8 && fabs(eta) <= 1.4442){ return 0.9919;}
    if (fabs(eta) > 1.566 && fabs(eta) <= 2.0){ return 1.0013;}
    if (fabs(eta) > 2.0 && fabs(eta) <= 2.5)  { return 1.0169;}
  }
  if( pt > 20 && pt <= 30 ){
    if (fabs(eta) <= 0.8)                     { return 0.9639;}
    if (fabs(eta) > 0.8 && fabs(eta) <= 1.4442){ return 0.9730;}
    if (fabs(eta) > 1.566 && fabs(eta) <= 2.0){ return 0.9835;}
    if (fabs(eta) > 2.0 && fabs(eta) <= 2.5)  { return 1.0046;}
  }
  if( pt > 30 && pt <= 40 ){
    if (fabs(eta) <= 0.8)                     { return 0.9764;}
    if (fabs(eta) > 0.8 && fabs(eta) <= 1.4442){ return 0.9777;}
    if (fabs(eta) > 1.566 && fabs(eta) <= 2.0){ return 0.9919;}
    if (fabs(eta) > 2.0 && fabs(eta) <= 2.5)  { return 1.0043;}
  }
  if( pt > 40 && pt <= 50 ){
    if (fabs(eta) <= 0.8)                     { return 0.9804;}
    if (fabs(eta) > 0.8 && fabs(eta) <= 1.4442){ return 0.9840;}
    if (fabs(eta) > 1.566 && fabs(eta) <= 2.0){ return 0.9959;}
    if (fabs(eta) > 2.0 && fabs(eta) <= 2.5)  { return 1.0066;}
  }
  if( pt > 50 ) {
    if (fabs(eta) <= 0.8)                     { return 0.9787;}
    if (fabs(eta) > 0.8 && fabs(eta) <= 1.4442){ return 0.9822;}
    if (fabs(eta) > 1.566 && fabs(eta) <= 2.0){ return 0.9973;}
    if (fabs(eta) > 2.0 && fabs(eta) <= 2.5)  { return 1.0075;}
  }
}

float analyze::SF_tightPhoID(float pt, float eta){
  if( pt > 15 && pt <= 20){
    if (fabs(eta) <= 0.8)                     { return 0.9496;}
    if (fabs(eta) > 0.8 && fabs(eta) <= 1.4442){ return 0.9803;}
    if (fabs(eta) > 1.566 && fabs(eta) <= 2.0){ return 1.0005;}
    if (fabs(eta) > 2.0 && fabs(eta) <= 2.5)  { return 1.0171;}
  }
  if( pt > 20 && pt <= 30 ){
    if (fabs(eta) <= 0.8)                     { return 0.9672;}
    if (fabs(eta) > 0.8 && fabs(eta) <= 1.4442){ return 0.9724;}
    if (fabs(eta) > 1.566 && fabs(eta) <= 2.0){ return 0.9867;}
    if (fabs(eta) > 2.0 && fabs(eta) <= 2.5)  { return 1.0130;}
  }
  if( pt > 30 && pt <= 40 ){
    if (fabs(eta) <= 0.8)                     { return 0.9711;}
    if (fabs(eta) > 0.8 && fabs(eta) <= 1.4442){ return 0.9688;}
    if (fabs(eta) > 1.566 && fabs(eta) <= 2.0){ return 0.9971;}
    if (fabs(eta) > 2.0 && fabs(eta) <= 2.5)  { return 1.0143;}
  }
  if( pt > 40 && pt <= 50 ){
    if (fabs(eta) <= 0.8)                     { return 0.9766;}
    if (fabs(eta) > 0.8 && fabs(eta) <= 1.4442){ return 0.9805;}
    if (fabs(eta) > 1.566 && fabs(eta) <= 2.0){ return 0.9996;}
    if (fabs(eta) > 2.0 && fabs(eta) <= 2.5)  { return 1.0129;}
  }
  if( pt > 50 ) {
    if (fabs(eta) <= 0.8)                     { return 0.9893;}
    if (fabs(eta) > 0.8 && fabs(eta) <= 1.4442){ return 0.9924;}
    if (fabs(eta) > 1.566 && fabs(eta) <= 2.0){ return 1.0004;}
    if (fabs(eta) > 2.0 && fabs(eta) <= 2.5)  { return 1.0019;}
  }
}

float analyze::PixelSeedVeto(float pt, float eta){

  if(fabs(eta) <= 1.4442){
    if(pt > 15 && pt <= 20) {return 0.996;}
    if(pt > 20 && pt <= 25) {return 0.994;}
    if(pt > 25 && pt <= 30) {return 0.996;}
    if(pt > 30 && pt <= 40) {return 0.999;}
    if(pt > 40 && pt <= 50) {return 1.009;}
    if(pt > 50 && pt <= 70) {return 0.993;}
    if(pt > 70 )            {return 1.047;}
  } else if(fabs(eta) > 1.566 && fabs(eta) <= 2.5){
    if(pt > 15 && pt <= 20) {return 0.960;}
    if(pt > 20 && pt <= 25) {return 0.977;}
    if(pt > 25 && pt <= 30) {return 0.951;}
    if(pt > 30 && pt <= 40) {return 1.029;}
    if(pt > 40 && pt <= 50) {return 0.971;}
    if(pt > 50 && pt <= 70) {return 0.965;}
    if(pt > 70 )           {return 1.145;}
  } 
}

void analyze::clearVariables(){
  trigger = 0;
  nvtx = 0;
  IDWeight = 1;

  ele_pt = 0;
  ele_phi = 0;
  ele_eta = 0;
  genPhoPt=-9;
  pho_pt = 0;
  pho_phi = 0;
  pho_eta = 0;
  pho_SCEta = 0;
  pho_Sihih = 0;
  pho_ChHadIso = 0;
  pho_SCRChIso = 0;
  pho_PFPhoIso = 0;
  pho_SCRPhoIso = 0;
  pho_RandConeChIso = 0;
  mc_truth = 0;
  PhoSihih_pass = 0;
  PhoIso_pass = 0;
  PhoChHadIso_pass = 0;
  W_MT=0;
  InvMass=0;
  PFMET=-9;
  PFMETPhi=-9;

  MatchRecoEleRecoJet  =-9;
  MatchRecoEle2RecoJet =-9;
  MatchRecoPhoRecoJet  =-9;
  RecoJetPt_drRjRp0p5 = -9;
  MatchRecoPhoGenPho = -9;
  MatchRecoEleGenPho = -9;
  MatchRecoEleGenEle = -9;
  MatchRecoPhoGenEle = -9;
  genPartPdgId_MatchRecoPhoGenEle = -999;
  genPartPdgId_MatchRecoEleGenEle = -999;
  genPartMomPdgId_MatchRecoPhoGenEle = -999;
  genPartMomPdgId_MatchRecoEleGenEle = -999;
  genPartPt_MatchRecoPhoGenPho = -999;
  genPartEta_MatchRecoPhoGenPho = -999;
  genPartPt_MatchRecoEleGenPho = -999;
  genPartEta_MatchRecoEleGenPho = -999;
  genPartPt_MatchRecoPhoGenEle = -999;
  genPartEta_MatchRecoPhoGenEle = -999;
  genPartPt_MatchRecoEleGenEle = -999;
  genPartEta_MatchRecoEleGenEle =-999;
  evt_weight=1;
  PhohasPixelSeed = 0;
  iPho1=-999;
  event=0;
  ilep1=-999;
  lepton1Pt=-99;
  lepton1Eta=-99;
  lepton1Phi=-99;
  lepton1SCEta=-99;
  lepton1SCPhi=-99;
  lepton1Sieie=-99;
  lepton1HoE = -99;
  lepton1RelIso=-99;
  
  lepton2VetoID = -99;
  lepton2Pt=-99;
  lepton2Eta=-99;
  lepton2Phi=-99;
  lepton2Sieie=-99;
  lepton2HoE = -99;
  lepton2RelIso=-99;
  run=-999 ;
}

float analyze::geteIDEA(int i, float eta){
//https://github.com/cms-analysis/EgammaAnalysis-ElectronTools/blob/master/interface/ElectronEffectiveArea.h#L114
  double EffectiveArea=0.0;
  // effective area for isolation 
  if (fabs(eta) >= 0.0 && fabs(eta) < 1.0 ) EffectiveArea = 0.130;
  if (fabs(eta) >= 1.0 && fabs(eta) < 1.479 ) EffectiveArea = 0.137;
  if (fabs(eta) >= 1.479 && fabs(eta) < 2.0 ) EffectiveArea = 0.067;
  if (fabs(eta) >= 2.0 && fabs(eta) < 2.2 ) EffectiveArea = 0.089;
  if (fabs(eta) >= 2.2 && fabs(eta) < 2.3 ) EffectiveArea = 0.107; 
  if (fabs(eta) >= 2.3 && fabs(eta) < 2.4 ) EffectiveArea = 0.110;
  if (fabs(eta) >= 2.4) EffectiveArea = 0.138;
  
  return EffectiveArea;
}

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification#Barrel_Cuts_eta_supercluster_1_4
bool analyze::eIDTight2012(int i){ 
  double eleEt = (*elePt)[i]; 
  double eleEta = (*eleSCEta)[i];
  float EffectiveArea = geteIDEA(i, (*eleSCEta)[i]);
  if (eleEt < 30.0 ) return false;
  if(fabs(eleEta)>1.4442 &&fabs(eleEta)<1.566 ) return false;
  if(fabs(eleEta)>2.5) return false;
  bool pass = true;

  if (fabs(eleEta) < 1.4442) {
    if (fabs((*eledEtaAtVtx)[i]) > 0.004) pass = false;
    if (fabs((*eledPhiAtVtx)[i]) > 0.03) pass = false;
    if ((*eleSigmaIEtaIEta)[i] > 0.01) pass = false;
    if ((*eleHoverE)[i] > 0.12) pass = false;
    if (!( fabs((*eleD0GV)[i]) < 0.02 )) pass = false;
    if (!( fabs((*eleDzGV)[i]) < 0.1 )) pass = false;
    if (!( fabs(1./(*eleSCEn)[i] - (*eleEoverP)[i]/(*eleSCEn)[i]) < 0.05 )) pass = false;
    if ((*eleConvVtxFit)[i] == 1) pass = false;
    if ((*eleMissHits)[i] > 0) pass = false;
    if ( ((*elePFChIso03)[i] + TMath::Max((Double_t) 0.,(Double_t) (*elePFPhoIso03)[i] + (*elePFNeuIso03)[i] - EffectiveArea*rho2012))/(*elePt)[i] > 0.10) pass = false;
}
  // //endcap
  else if(fabs(eleEta) > 1.566 && fabs(eleEta) <2.5){
    if (fabs((*eledEtaAtVtx)[i]) > 0.005) pass = false;
    if (fabs((*eledPhiAtVtx)[i]) > 0.02) pass = false;
    if ((*eleSigmaIEtaIEta)[i] > 0.03) pass = false;
    if ((*eleHoverE)[i] > 0.10) pass = false;
    if (!( fabs((*eleD0GV)[i]) < 0.02 )) pass = false;
    if (!( fabs((*eleDzGV)[i]) < 0.1 )) pass = false;
    if (!( fabs(1./(*eleSCEn)[i] - (*eleEoverP)[i]/(*eleSCEn)[i]) < 0.05 )) pass = false; 
    if ((*eleConvVtxFit)[i] == 1) pass = false;
    if ((*eleMissHits)[i] > 0) pass = false;
    if ((*elePt)[i] < 20){
     if ( ((*elePFChIso03)[i] + TMath::Max((Double_t) 0., (Double_t) (*elePFPhoIso03)[i] + (*elePFNeuIso03)[i] - EffectiveArea*rho2012))/(*elePt)[i] > 0.07) pass = false;
    } else {
      if ( ((*elePFChIso03)[i] + TMath::Max((Double_t) 0., (Double_t) (*elePFPhoIso03)[i] + (*elePFNeuIso03)[i] - EffectiveArea*rho2012))/(*elePt)[i] > 0.10) pass = false;
    }
  }
  return pass;
}
bool analyze::eIDLoose2012(int i){
  double eleEt = (*elePt)[i];
  double eleEta = (*eleSCEta)[i];
  float EffectiveArea = geteIDEA(i, (*eleSCEta)[i]);
  if (eleEt < 30.0 ) return false;
  if(fabs(eleEta)>1.4442 &&fabs(eleEta)<1.566 ) return false;
  if(fabs(eleEta)>2.5) return false;
  bool pass = true;

  if (fabs(eleEta) < 1.4442) {
    if (fabs((*eledEtaAtVtx)[i]) > 0.007) pass = false;
    if (fabs((*eledPhiAtVtx)[i]) > 0.15) pass = false;
    if ((*eleSigmaIEtaIEta)[i] > 0.01) pass = false;
    if ((*eleHoverE)[i] > 0.12) pass = false;
    if (!( fabs((*eleD0GV)[i]) < 0.02 )) pass = false;
    if (!( fabs((*eleDzGV)[i]) < 0.2 )) pass = false;
    if (!( fabs(1./(*eleSCEn)[i] - (*eleEoverP)[i]/(*eleSCEn)[i]) < 0.05 )) pass = false;
    if ((*eleConvVtxFit)[i] == 1) pass = false;
    if ((*eleMissHits)[i] > 1) pass = false;
    if ( ((*elePFChIso03)[i] + TMath::Max((Double_t) 0.,(Double_t) (*elePFPhoIso03)[i] + (*elePFNeuIso03)[i] - EffectiveArea*rho2012))/(*elePt)[i] > 0.15) pass = false;
}
  // //endcap
  else if(fabs(eleEta) > 1.566 && fabs(eleEta) <2.5){
    if (fabs((*eledEtaAtVtx)[i]) > 0.009) pass = false;
    if (fabs((*eledPhiAtVtx)[i]) > 0.10) pass = false;
    if ((*eleSigmaIEtaIEta)[i] > 0.03) pass = false;
    if ((*eleHoverE)[i] > 0.10) pass = false;
    if (!( fabs((*eleD0GV)[i]) < 0.02 )) pass = false;
    if (!( fabs((*eleDzGV)[i]) < 0.2 )) pass = false;
    if (!( fabs(1./(*eleSCEn)[i] - (*eleEoverP)[i]/(*eleSCEn)[i]) < 0.05 )) pass = false;
    if ((*eleConvVtxFit)[i] == 1) pass = false;
    if ((*eleMissHits)[i] > 1) pass = false;
    if ((*elePt)[i] < 20){
     if ( ((*elePFChIso03)[i] + TMath::Max((Double_t) 0., (Double_t) (*elePFPhoIso03)[i] + (*elePFNeuIso03)[i] - EffectiveArea*rho2012))/(*elePt)[i] > 0.10) pass = false;
    } else {
      if ( ((*elePFChIso03)[i] + TMath::Max((Double_t) 0., (Double_t) (*elePFPhoIso03)[i] + (*elePFNeuIso03)[i] - EffectiveArea*rho2012))/(*elePt)[i] > 0.15) pass = false;
    }
  }
  return pass;
}

bool analyze::eIDMedium2012(int i){
  double eleEt = (*elePt)[i];
  double eleEta = (*eleSCEta)[i];
  float EffectiveArea = geteIDEA(i, (*eleSCEta)[i]);
  //  cout<<eleEta<<'\t'<<EffectiveArea<<endl;                                                                                         
  if (eleEt < 30.0 ) return false;
  if(fabs(eleEta)>1.4442 &&fabs(eleEta)<1.566 ) return false;
  if(fabs(eleEta)>2.5) return false;
  bool pass = true;

  if (fabs(eleEta) < 1.4442) {
    if (fabs((*eledEtaAtVtx)[i]) > 0.004) pass = false;
    if (fabs((*eledPhiAtVtx)[i]) > 0.06) pass = false;
    if ((*eleSigmaIEtaIEta)[i] > 0.01) pass = false;
    if ((*eleHoverE)[i] > 0.12) pass = false;
    //    if (fabs((*eleD0Vtx)[i][0]) > 0.02) pass = false;                                                                         
    //    if (fabs((*eleDzVtx)[i][0]) > 0.1) pass = false;                                                                          
    //  if (fabs(1./(*eleEcalEn)[i] - 1./(*elePin)[i]) > 0.05) pass = false;                                                        
    //Sachiko changed due to Wgg. Sep 8, 2015                                                                                       
    if (!( fabs((*eleD0GV)[i]) < 0.02 )) pass = false;
    if (!( fabs((*eleDzGV)[i]) < 0.1 )) pass = false;
    if (!( fabs(1./(*eleSCEn)[i] - (*eleEoverP)[i]/(*eleSCEn)[i]) < 0.05 )) pass = false;
    if ((*eleConvVtxFit)[i] == 1) pass = false;
    if ((*eleMissHits)[i] > 1) pass = false;
    if ( ((*elePFChIso03)[i] + TMath::Max((Double_t) 0.,(Double_t) (*elePFPhoIso03)[i] + (*elePFNeuIso03)[i] - EffectiveArea*rho2012))\
	 /(*elePt)[i] > 0.15) pass = false;

  }
  // //endcap                                                                                                                          
  else if(fabs(eleEta) > 1.566 && fabs(eleEta) <2.5){
    if (fabs((*eledEtaAtVtx)[i]) > 0.007) pass = false;
    if (fabs((*eledPhiAtVtx)[i]) > 0.03) pass = false;
    if ((*eleSigmaIEtaIEta)[i] > 0.03) pass = false;
    if ((*eleHoverE)[i] > 0.10) pass = false;
    //    if (fabs((*eleD0Vtx)[i][0]) > 0.02) pass = false;
    //    if (fabs((*eleDzVtx)[i][0]) > 0.1) pass = false;
    //    if (fabs(1./(*eleEcalEn)[i] - 1./(*elePin)[i]) > 0.05) pass = false;
    if (!( fabs((*eleD0GV)[i]) < 0.02 )) pass = false;
    if (!( fabs((*eleDzGV)[i]) < 0.1 )) pass = false;
    if (!( fabs(1./(*eleSCEn)[i] - (*eleEoverP)[i]/(*eleSCEn)[i]) < 0.05 )) pass = false;

    if ((*eleConvVtxFit)[i] == 1) pass = false;
    if ((*eleMissHits)[i] > 1) pass = false;
    if ((*elePt)[i] < 20){
      if ( ((*elePFChIso03)[i] + TMath::Max((Double_t) 0., (Double_t) (*elePFPhoIso03)[i] + (*elePFNeuIso03)[i] - EffectiveArea*rho2012))/(*elePt)[i] > 0.10) pass = false;
    } else {
      if ( ((*elePFChIso03)[i] + TMath::Max((Double_t) 0., (Double_t) (*elePFPhoIso03)[i] + (*elePFNeuIso03)[i] - EffectiveArea*rho2012))/(*elePt)[i] > 0.15) pass = false;
    }
  }
  return pass;
}

bool analyze::eIDVeto2012(int i){
  double eleEt = (*elePt)[i];
  double eleEta = (*eleSCEta)[i];
  float EffectiveArea = geteIDEA(i, (*eleSCEta)[i]);
  //cout<<eleEta<<'\t'<<EffectiveArea<<endl;  
  if (eleEt < 10.0 ) return false;
  if(fabs(eleEta)>1.4442 &&fabs(eleEta)<1.566 ) return false;
  if(fabs(eleEta)>2.5) return false;
  bool pass = true;
  
  if (fabs(eleEta) < 1.4442) {
    if (fabs((*eledEtaAtVtx)[i]) > 0.007) pass = false;
    if (fabs((*eledPhiAtVtx)[i]) > 0.8) pass = false;
    if ((*eleSigmaIEtaIEta)[i] > 0.01) pass = false;
//    if ((*eleSigmaIEtaIEta)[i] > 0.015) pass = false;
    if ((*eleHoverE)[i] > 0.15) pass = false;
    if (fabs((*eleD0GV)[i]) > 0.04) pass = false;
    if (fabs((*eleDzGV)[i]) > 0.2) pass = false;
//    if ( ((*elePFChIso03)[i] + TMath::Max((Double_t) 0.,(Double_t) (*elePFPhoIso03)[i] + (*elePFNeuIso03)[i] - EffectiveArea*rho2012))/(*elePt)[i] > 0.15) pass = false;
}
  // //endcap
  else {
     if (fabs((*eledEtaAtVtx)[i]) > 0.01) pass = false;
//cout<<fabs((*eledEtaAtVtx)[i]) << '\t'<<pass<<endl;
     if (fabs((*eledPhiAtVtx)[i]) > 0.7) pass = false;
//cout<<pass<<endl;
     if ((*eleSigmaIEtaIEta)[i] > 0.03) pass = false;
//cout<<pass<<endl;
     if (fabs((*eleD0GV)[i]) > 0.04) pass = false;
//cout<<pass<<endl;
     if (fabs((*eleDzGV)[i]) > 0.2) pass = false;
//cout<<pass<<endl;
  //   if ( ((*elePFChIso03)[i] + TMath::Max((Double_t) 0., (Double_t) (*elePFPhoIso03)[i] + (*elePFNeuIso03)[i] - EffectiveArea*rho2012))/(*elePt)[i] > 0.15) pass = false;
//cout<<pass<<endl;

  }
  return pass;
}


// MCTruthMatch is for closure test
int analyze::MCTruthMatch(int jpho){

  int phoInd = -1;
  for(int imc = 0; imc < nMC; ++imc){

    if((*mcPID)[imc] != 22) continue;
    
    bool match_gen = dR((*mcEta)[imc], (*mcPhi)[imc], (*phoEta)[jpho], (*phoPhi)[jpho]) < 0.1 &&
      (fabs((*phoEt)[jpho] - (*mcPt)[imc]) / (*mcPt)[imc] < 1.0);
    
    if(match_gen && phoInd < 0) phoInd = imc;
  }

  if(phoInd >= 0){
    if(((*mcParentage)[phoInd]& 4)==0) return 1; // true signal
    else 
      return 2; // true background
  } else {
    return 3; // fake background
  }
}

void analyze::Loop(TString name)
{
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  TFile* tmp = TFile::Open(name, "RECREATE");
  
  TTree* miniTree = new TTree("miniTree", "miniTree");
  /*  miniTree->Branch("trigger", &trigger, "trigger/I");
      miniTree->Branch("nvtx", &nvtx, "nvtx/I");
      // electrons
      miniTree->Branch("ele_pt", &ele_pt, "ele_pt/F");
      miniTree->Branch("ele_phi", &ele_phi, "ele_phi/F");
      miniTree->Branch("ele_eta", &ele_eta, "ele_eta/F");
      // photons
      miniTree->Branch("pho_PFPhoIso", &pho_PFPhoIso, "pho_PFPhoIso/F");
      miniTree->Branch("pho_SCRPhoIso", &pho_SCRPhoIso, "pho_SCRPhoIso/F");
      
      miniTree->Branch("mc_truth", &mc_truth, "mc_truth/I");
      miniTree->Branch("PhoSihih_pass", &PhoSihih_pass, "PhoSihih_pass/I");
      miniTree->Branch("PhoIso_pass", &PhoIso_pass, "PhoIso_pass/I");
      miniTree->Branch("PhoChHadIso_pass", &PhoChHadIso_pass, "PhoChHadIso_pass/I");
      
      miniTree->Branch("totalWeight", &totalWeight, "totalWeight/F");
      miniTree->Branch("totalWeight_m5", &totalWeight_m5, "totalWeight_m5/F");
      miniTree->Branch("totalWeight_p5", &totalWeight_p5, "totalWeight_p5/F");
      miniTree->Branch("processCode", &processCode, "processCode/I");
      miniTree->Branch("xsWeight", &xsWeight, "xsWeight/F");
      miniTree->Branch("puWeight", &puWeight, "puWeight/F");
      miniTree->Branch("puWeight_m5", &puWeight_m5, "puWeight_m5/F");
      miniTree->Branch("puWeight_p5", &puWeight_p5, "puWeight_p5/F");
      miniTree->Branch("IDWeight", &IDWeight, "IDWeight/F");
  */
  miniTree->Branch("event", &event);
  miniTree->Branch("run", &run);
  miniTree->Branch("ilep1", &ilep1);
  miniTree->Branch("iPho1", &iPho1);


  miniTree->Branch("lepton1Pt", &lepton1Pt, "lepton1Pt/F");
  miniTree->Branch("lepton1Eta", &lepton1Eta, "lepton1Eta/F");
  miniTree->Branch("lepton1Phi", &lepton1Phi, "lepton1Phi/F");
  miniTree->Branch("lepton1SCEta", &lepton1SCEta, "lepton1SCEta/F");
  miniTree->Branch("lepton1SCPhi", &lepton1SCPhi, "lepton1SCPhi/F");
  miniTree->Branch("lepton1Sieie", &lepton1Sieie, "lepton1Sieie/F");
  miniTree->Branch("lepton1HoE", &lepton1HoE, "lepton1HoE/F");
  miniTree->Branch("lepton1RelIso", &lepton1RelIso, "lepton1RelIso/F");

  miniTree->Branch("lepton2VetoID",&lepton2VetoID,"lepton2VetoID/I");
  miniTree->Branch("lepton2Pt", &lepton2Pt, "lepton2Pt/F");
  miniTree->Branch("lepton2Eta", &lepton2Eta, "lepton2Eta/F");
  miniTree->Branch("lepton2Phi", &lepton2Phi, "lepton2Phi/F");
  miniTree->Branch("lepton2Sieie", &lepton2Sieie, "lepton2Sieie/F");
  miniTree->Branch("lepton2HoE", &lepton2HoE, "lepton2HoE/F");
  miniTree->Branch("lepton2RelIso", &lepton2RelIso, "lepton2RelIso/F");

  miniTree->Branch("evt_weight", &evt_weight, "evt_weight/F");
  miniTree->Branch("genPhoPt",&genPhoPt,"genPhoPt/F");
  miniTree->Branch("pho_pt", &pho_pt, "pho_pt/F");
  miniTree->Branch("pho_eta", &pho_eta, "pho_eta/F");
  miniTree->Branch("pho_phi", &pho_phi, "pho_phi/F");
  miniTree->Branch("pho_SCEta", &pho_SCEta, "pho_SCEta/F");
  miniTree->Branch("pho_Sihih", &pho_Sihih, "pho_Sihih/F");
  miniTree->Branch("pho_ChHadIso", &pho_ChHadIso, "pho_ChHadIso/F");
  miniTree->Branch("pho_SCRChIso", &pho_SCRChIso, "pho_SCRChIso/F");
  miniTree->Branch("pho_RandConeChIso", &pho_RandConeChIso, "pho_RandConeChIso/F");
  miniTree->Branch("W_MT", &W_MT, "W_MT/F");
  miniTree->Branch("PFMET", &PFMET, "PFMET/F");
  miniTree->Branch("PFMETPhi", &PFMETPhi, "PFMETPhi/F");

  miniTree->Branch("MatchRecoEleRecoJet", &MatchRecoEleRecoJet, "MatchRecoEleRecoJet/F");
  miniTree->Branch("MatchRecoEle2RecoJet", &MatchRecoEle2RecoJet, "MatchRecoEle2RecoJet/F");
  miniTree->Branch("MatchRecoPhoRecoJet", &MatchRecoPhoRecoJet, "MatchRecoPhoRecoJet/F");
  miniTree->Branch("RecoJetPt_drRjRp0p5", &RecoJetPt_drRjRp0p5, "RecoJetPt_drRjRp0p5/F");

  miniTree->Branch("MatchRecoPhoGenPho", &MatchRecoPhoGenPho, "MatchRecoPhoGenPho/F");
  miniTree->Branch("genPartPt_MatchRecoPhoGenPho", &genPartPt_MatchRecoPhoGenPho, "genPartPt_MatchRecoPhoGenPho/F");
  miniTree->Branch("genPartEta_MatchRecoPhoGenPho", &genPartEta_MatchRecoPhoGenPho, "genPartEta_MatchRecoPhoGenPho/F");
  miniTree->Branch("MatchRecoEleGenPho", &MatchRecoEleGenPho, "MatchRecoEleGenPho/F");
  miniTree->Branch("genPartPt_MatchRecoEleGenPho", &genPartPt_MatchRecoEleGenPho, "genPartPt_MatchRecoEleGenPho/F");
  miniTree->Branch("genPartEta_MatchRecoEleGenPho", &genPartEta_MatchRecoEleGenPho, "genPartEta_MatchRecoEleGenPho/F");
  miniTree->Branch("MatchRecoEleGenEle", &MatchRecoEleGenEle, "MatchRecoEleGenEle/F");
  miniTree->Branch("MatchRecoPhoGenEle", &MatchRecoPhoGenEle, "MatchRecoPhoGenEle/F");
  miniTree->Branch("genPartPdgId_MatchRecoPhoGenEle", &genPartPdgId_MatchRecoPhoGenEle, "genPartPdgId_MatchRecoPhoGenEle/I");
  miniTree->Branch("genPartMomPdgId_MatchRecoPhoGenEle", &genPartMomPdgId_MatchRecoPhoGenEle, "genPartMomPdgId_MatchRecoPhoGenEle/I");
  miniTree->Branch("genPartPt_MatchRecoPhoGenEle", &genPartPt_MatchRecoPhoGenEle, "genPartPt_MatchRecoPhoGenEle/F");
  miniTree->Branch("genPartEta_MatchRecoPhoGenEle", &genPartEta_MatchRecoPhoGenEle, "genPartEta_MatchRecoPhoGenEle/F");

  miniTree->Branch("genPartPdgId_MatchRecoEleGenEle", &genPartPdgId_MatchRecoEleGenEle, "genPartPdgId_MatchRecoEleGenEle/I");
  miniTree->Branch("genPartMomPdgId_MatchRecoEleGenEle", &genPartMomPdgId_MatchRecoEleGenEle, "genPartMomPdgId_MatchRecoEleGenEle/I");
  miniTree->Branch("genPartPt_MatchRecoEleGenEle", &genPartPt_MatchRecoEleGenEle, "genPartPt_MatchRecoEleGenEle/F");
  miniTree->Branch("genPartEta_MatchRecoEleGenEle", &genPartEta_MatchRecoEleGenEle, "genPartEta_MatchRecoEleGenEle/F");
  
  miniTree->Branch("PhohasPixelSeed",&PhohasPixelSeed,"PhohasPixelSeed/I");
  miniTree->Branch("InvMass", &InvMass, "InvMass/F");

  TH1F *WMT=new TH1F("WMT","Transverse mass of W",60,0.,300);
  TH1F *MET=new TH1F("MET","MET",30,0.,300);
  TH1F *pT_gamma=new TH1F("pT_gamma","pT_gamma",30,0.,300);
  const Int_t NBINS = 9;
  Double_t edges[NBINS + 1] = {15,20,25,30,35,40,55,75,95,500};
  TH1F *pT_photon=new TH1F("pT_photon","pT_photon",NBINS,edges);
  
  
  TH1F *h_elec_pho_dR=new TH1F("h_elec_pho_dR","h_elec_pho_dR",50,0.,5.);
  int countWenu = 0, countWmunu=0, countWtaunu=0, leftEvents=0, pfmet30=0,WMt50=0,WMt70=0;
  TH1F *h_left=new TH1F("h_left","h_left",10,0,10);
  TH1F* h_nEle=new TH1F("h_nEle","h_nEle",10,0,10);
  TH1F* h_nSelElec=new TH1F("h_nSelElec","h_nSelElec",10,0,10);
  TH1F* h_SelleadElecPt=new TH1F("h_SelleadElecPt","h_SelleadElecPt",30,0.,300.);
  TH1F* h_SelsubleadElecPt=new TH1F("h_SelsubleadElecPt","h_SelsubleadElecPt",30,0.,300.);

  TH1F* h_pfMET=new TH1F("h_pfMET","h_pfMET",30,0.,300);
  TH1F* h_WMT=new TH1F("h_WMT","h_WMT",30,0.,300);
  TH1F* h_pfMET_EB=new TH1F("h_pfMET_EB","h_pfMET_EB",30,0.,300);
  TH1F* h_WMT_EB=new TH1F("h_WMT_EB","h_WMT_EB",30,0.,300);
  TH1F* h_pfMET_EE=new TH1F("h_pfMET_EE","h_pfMET_EE",30,0.,300);
  TH1F* h_WMT_EE=new TH1F("h_WMT_EE","h_WMT_EE",30,0.,300);
  TH1F* h_WMT_EBmet35=new TH1F("h_WMT_EBmet35","h_pfMET_EBmet35",30,0.,300);
  TH1F* h_WMT_EEmet35=new TH1F("h_WMT_EEmet35","h_pfMET_EEmet35",30,0.,300);

  TH1F* h_InvMass=new TH1F("h_InvMass","h_InvMass",30,0.,300);
  TH1F* h_HasMatchingGSFelectron=new TH1F("h_HasMatchingGSFelectron","h_HasMatchingGSFelectron",100.,0.,5);
  TH1F* hCounter=new TH1F("hCounter","hCounter",20,0.,20);
  TH1F *h_genPartPdgId=new TH1F("h_genPartPdgId","h_genPartPdgId",100,0,100); 
  TH1F* h_checkOverlap=new TH1F("h_checkOverlap","h_checkOverlap", 100,0.,4.);
  TH1F* h_dr_recogenPho=new TH1F("h_dr_recogenPho","h_dr_recogenPho", 100,0.,4.);
  
  Long64_t nbytes = 0, nb = 0;
  
  int overlapping=0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    
    clearVariables();
    
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    //for blinding
//      if(processCode==0 && !(event % 20 == 0))continue;
  //want to separate W->enu and W->taunu
    bool isWenu=0;
    bool isWmunu = 0;
    bool isWtaunu=0;
    if(processCode==1) 
      {
	for (int iMC=0; iMC<nMC; iMC++)
	  {
	    if (abs(mcPID->at(iMC))==24) //W PDGID
	      {
		if (mcDecayType->at(iMC)==2) isWenu = 1; 
		else if (mcDecayType->at(iMC)==3) isWmunu = 1;
		else if (mcDecayType->at(iMC)==4)isWtaunu = 1;
	      }    
	  }
      }
    // if (Cut(ientry) < 0) continue;
    if(jentry % 10000 == 0) cout << "Processed " << jentry
				 << " events out of " <<nentries<< endl;
    //if(!(run == 194644  && event == 408284122))continue;
    // if( run==193336 && event == 115039380) cout<<run<<'\t'<<event<<endl;
    //if(!(event== 1075101220 ))continue;
//if(!((event==11127540 || event == 14128141) && run ==190706))continue;
    //cout<<event<<endl;
    //if(processCode==0 )
    //met filters as Yutaro do
    hCounter->Fill(0);
    if (!metFilters[0]) continue; //[0]CSC Beam Halo
    hCounter->Fill(1);
    if (!metFilters[1]) continue; //[1]HBHE Noise
    hCounter->Fill(2);
    if (!metFilters[2]) continue; //[2]HCAL laser
    hCounter->Fill(3);
    if (!metFilters[3]) continue; //[3]ECAL dead cell
    hCounter->Fill(4);
    if (!metFilters[4]) continue; //[4]Tracking failure
    hCounter->Fill(5);
    if (!metFilters[5]) continue; //[5]Bad EE SC
    hCounter->Fill(6);
    if(!(IsVtxGood>=0)) continue;
    nvtx = nVtx;
    //Single Electron Trigger
    trigger = 0;
    hCounter->Fill(7);
    //cout<<HLT[HLTIndex[17]]<<endl;
    if(HLT[HLTIndex[17]] != -1 && HLT[HLTIndex[17]] != 0 ) trigger = 1;
    if(!trigger){
      continue;
    }
    hCounter->Fill(8);

    //katya do so
    //    if(nPho<1 || nEle<1) continue; 
    if( nEle<1) continue; 
    hCounter->Fill(9);
  //  cout<<nPho<<'\t'<<nEle<<endl;
    // LOOP OVER EM OBJECTS
    
    //h_nEle->Fill(nEle,totalWeight);
    vector <int> ielectrons;
    int ecount=0;
    for (int iele = 0; iele < nEle; ++iele){
      //cout<<iele<<'\t'<<(*elePt)[iele]<<'\t'<<eleTrg->at(iele)%2<<'\t'<<(*eleSCEta)[iele]<<'\t'<<eIDVeto2012(iele)<<endl;
      //      if(!(eleTrg->at(iele)%2))continue;// {cout<<nEle<<'\t'<<iele<<" matched"<<endl;}
      
      if ((*elePt)[iele] < 10) continue;
      //      if ((*elePt)[iele] < 30) continue;
      if ( fabs((*eleSCEta)[iele]) > 1.4442 && fabs((*eleSCEta)[iele]) < 1.566 ) continue;
      if ( fabs((*eleSCEta)[iele]) > 2.5 ) continue;
      //   if( !eIDTight2012(iele) )continue;
      if( !eIDVeto2012(iele) )continue;
      ielectrons.push_back(iele);
    //  cout<<"aftr cuts: "<<iele<<'\t'<<(*elePt)[iele]<<'\t'<<eleTrg->at(iele)%2<<'\t'<<(*eleSCEta)[iele]<<'\t'<<eIDVeto2012(iele)<<'\t'<<eIDTight2012(iele)<<endl;
      
    }
    // if ( ielectrons.size()<1 ) continue;
  //  cout<<ielectrons.size()<<endl;
    h_nSelElec->Fill(ielectrons.size(),totalWeight);
    if(ielectrons.size()>0)h_SelleadElecPt->Fill(elePt->at(ielectrons[0]),totalWeight);
    if(ielectrons.size()>1)h_SelsubleadElecPt->Fill(elePt->at(ielectrons[1]),totalWeight);
    if ( ielectrons.size()<1 ) continue;
    hCounter->Fill(10);
//cout<<eIDTight2012(ielectrons[0])<<endl;

    if( !eIDTight2012(ielectrons[0])) continue;
    hCounter->Fill(11);

if(elePt->at(ielectrons[0])<30.)continue;    
    if(!(eleTrg->at(ielectrons[0])%2))continue;// {cout<<nEle<<'\t'<<iele<<" matched"<<endl;}
    hCounter->Fill(12);
    
    if(ielectrons.size()>1 && elePt->at(ielectrons[1])>10 && eIDVeto2012(ielectrons[1])==1 ) continue;
//    if(ielectrons.size()>1 && elePt->at(ielectrons[1])>10) continue;
    hCounter->Fill(13);
    
    
    //cout<<"all: "<<ielectrons.size()<<'\t'<<elePt->at(ielectrons[0])<<endl;
    if(ielectrons.size()>1 && ( fabs((*eleSCEta)[ielectrons[1]]) > 1.4442 && fabs((*eleSCEta)[ielectrons[1]]) < 1.566 ) ) continue;
    hCounter->Fill(14);
    
    TLorentzVector ele;
    int i=-9;
    if( eIDTight2012(ielectrons[0])) i = ielectrons[0];
    ele.SetPtEtaPhiM((*elePt)[i], (*eleEta)[i], (*elePhi)[i], 0);

    TLorentzVector ele2;
    int i2 = -9;
    if(ielectrons.size()>1)
	{
	  ele2.SetPtEtaPhiM((*elePt)[ielectrons[1]], (*eleEta)[ielectrons[1]], (*elePhi)[ielectrons[1]], 0);
	 i2 = ielectrons[1]; 	
	}
    //else ele2.SetPtEtaPhiM(0,-9,-9,0);
    //cout<<(*elePt)[i]<<'\t'<<(*eleEta)[i]<<'\t'<<(*elePhi)[i]<<endl;   
    ilep1=i;
    lepton1Pt = (*elePt)[i];
    lepton1SCEta = (*eleSCEta)[i];
    lepton1Eta = (*eleEta)[i];
    lepton1Phi = (*elePhi)[i];
    lepton1SCPhi = (*eleSCPhi)[i];
    float EffectiveArea1 = geteIDEA(i, (*eleSCEta)[i]);
    float EffectiveArea2 = geteIDEA(i2, (*eleSCEta)[i2]);

    lepton1RelIso = ((*elePFChIso03)[i] + TMath::Max((Double_t) 0., (Double_t) (*elePFPhoIso03)[i] + (*elePFNeuIso03)[i] - EffectiveArea1*rho2012))/(*elePt)[i];
    lepton1Sieie = (*eleSigmaIEtaIEta)[i];
    lepton1HoE = (*eleHoverE)[i];

    if(ielectrons.size()>1) {
	lepton2RelIso = ((*elePFChIso03)[i2] + TMath::Max((Double_t) 0., (Double_t) (*elePFPhoIso03)[i2] + (*elePFNeuIso03)[i2] - EffectiveArea2*rho2012))/(*elePt)[i2];
     lepton2Pt = elePt->at(i2);
     lepton2Eta = eleEta->at(i2);
     lepton2VetoID = eIDVeto2012(i2);
         lepton2Sieie = (*eleSigmaIEtaIEta)[i2];
    lepton2HoE = (*eleHoverE)[i2];
}
    //for (int m=0;m<10;m++)if(metFilters[m])continue;
    //h_pfMET->Fill(pfMET,totalWeight);
    //if(pfMET<40.)continue;
    
    //overlapping events??
    bool overlapped=0;
    if((processCode==2) || (processCode==4) || (processCode==5) || (processCode==6)){
      float dr_min = 99999;
      float dr;
      for (int genPho = 0; genPho < nMC; ++genPho){
	if( (*mcPID)[genPho] != 22 ) continue;
	if( ((*mcParentage)[genPho]& 4) == 4) continue;
	if( (*mcPt)[genPho] < 10 ) continue;
	for (int genEle = 0; genEle < nMC; ++genEle){
	   if( fabs((*mcPID)[genEle]) != 11 ) continue;
	   dr = dR((*mcEta)[genEle],(*mcPhi)[genEle],(*mcEta)[genPho],(*mcPhi)[genPho]);
            if(dr < dr_min) dr_min = dr;
	}//end loop over genEle
      }//end loop over for gen-photon

      h_checkOverlap->Fill(dr_min);
      if(dr_min > 0.4 && dr_min != 99999 ) {
	overlapped = true;
      }	  
    }
    
    if(overlapped){overlapping++; continue;} 
    
    //if(pfMET<70.) continue;
    double mass_w = sqrt( 2 * (*elePt)[i] * pfMET * (1 - cos(dR(0.0, (*elePhi)[i], 0.0, pfMETPhi))) );
    h_WMT->Fill(mass_w,totalWeight);
    h_pfMET->Fill(pfMET,totalWeight);
    if((fabs((*eleEta)[i]))<1.4442)
      {
	h_WMT_EB->Fill(mass_w,totalWeight);
	h_pfMET_EB->Fill(pfMET,totalWeight);
	if(pfMET>35)h_WMT_EBmet35->Fill(mass_w,totalWeight);
      }
    else 
      {
	h_WMT_EE->Fill(mass_w,totalWeight);
	h_pfMET_EE->Fill(pfMET,totalWeight);
	if(pfMET>35)h_WMT_EEmet35->Fill(mass_w,totalWeight);
      }
    //cout<<mass_w<<endl; 
    //    if(mass_w < 40.) continue;
    hCounter->Fill(15);
    
    // LOOP OVER FOR PHOTON
    vector <int> iphotons;
    vector <int> pre_photons;
    vector <int> sihih_pass;
    vector <int> phoIso_pass;
    vector <int> chHadIso_pass;
    if(nPho<1) continue;
    for (int ipho = 0; ipho < nPho; ++ipho){
      
      double Pho03NeuHadIso =    (*phoPFNeuIso)[ipho]  - rho2012 * phoEffArea03NeuHad((*phoSCEta)[ipho]);
      double Pho03PhoIso =       (*phoPFPhoIso)[ipho]  - rho2012 * phoEffArea03Pho((*phoSCEta)[ipho]);
      
      float Pho03ChHadIso = (*phoPFChIso)[ipho] - rho2012 * phoEffArea03ChHad((*phoSCEta)[ipho]);
      float Pho03ChHadSCRIso = (*phoSCRChIso)[ipho] - rho2012 * phoEffArea03ChHad((*phoSCEta)[ipho]);
      //cout<<ipho<<", pho: pt "<<(*phoEt)[ipho]<<", eta: "<<(*phoSCEta)[ipho]
      /*<<", phi: "<<(*phoPhi)[ipho]<<", ID: "<<passPhotonID(ipho, 1)
	<<"phoEleVeto: "<<(*phoEleVeto)[ipho]
	<<", phoHoverE12: "<<(*phoHoverE12)[ipho]
	<<", phoPFNeuIso: "<<(*phoPFNeuIso)[ipho]
	<<", rho2012: "<< rho2012
	<<", phoEffArea03NeuHad: "<< phoEffArea03NeuHad((*phoSCEta)[ipho])
	<<", Pho03NeuHadIso: "<<Pho03NeuHadIso
	<<", phoPFPhoIso: "<<(*phoPFPhoIso)[ipho] 
	<<", phoEffArea03Pho: "<<phoEffArea03Pho((*phoSCEta)[ipho])
	<<", Pho03PhoIso: "<<Pho03PhoIso
	<<endl;
      */
      
      // PRE-PHOTON SELECTION
      if((*phoEt)[ipho] < 15) continue;
      if( fabs((*phoSCEta)[ipho])>2.5) continue;
      if( fabs((*phoSCEta)[ipho])<1.566 && fabs((*phoSCEta)[ipho])>1.4442) continue;
      if (!passPhotonID(ipho, 1)) continue;
      TLorentzVector pho;
      pho.SetPtEtaPhiM((*phoEt)[ipho], (*phoEta)[ipho], (*phoPhi)[ipho], 0);
      // REJECT PHOTONS THAT ARE TOO CLOSE TO ELECTRON
      if( pho.DeltaR(ele) < 0.7 ) continue;
      TLorentzVector ElecPhopair;
      ElecPhopair = pho+ele;
      InvMass=ElecPhopair.M();
      //cout<<"Zmass: "<<InvMass<<endl;
      //    if(InvMass>70. && InvMass<110)continue;
      //hCounter->Fill(17);
      int region = 1;
      if(fabs((*phoSCEta)[ipho]) < 1.4442) region = 0;
      
      //    bool pass_all = (*phoSigmaIEtaIEta)[ipho] < photonID_SigmaIEtaIEta[region][pho_ID];
      // &&      	Pho03ChHadIso < photonID_RhoCorrR03ChHadIso[region][pho_ID];
      //         if(!pass_all)continue;
      
      //      if(phohasPixelSeed)continue;
      //check overlapping photons
      /*    if((processCode==2) || (processCode==4) || (processCode==5) || (processCode==6)){
	    bool pho_matched_dy = false;
	    float dr_min = 1000.;
	    float dr;
	    for (int genPho = 0; genPho < nMC; ++genPho){
	    if( (*mcPID)[genPho] != 22 ) continue;
	    if( ((*mcParentage)[genPho]& 4) == 4) continue;
	    if( (*mcPt)[genPho] < 10 ) continue;
	    // CALCULATE DR(gen-pho,reco-pho)
	    float deta = (*phoEta)[ipho] - (*mcEta)[genPho];
	    float dphi = acos(cos((*phoPhi)[ipho] - (*mcPhi)[genPho]));
	    dr = sqrt(deta*deta + dphi*dphi);
	    if(dr < dr_min) dr_min = dr;
	    if(dr < 0.2) {
	    h_dr_recogenPho->Fill(dr,totalWeight);
	    pho_matched_dy = true;
	    break;
	    }
	    } // loop over for gen-photon
	    if(pho_matched_dy) continue;
	    } // loop over for processCode
	    //end check overlapping photons*/
      iphotons.push_back(ipho);
    } // loop over for photon
    
    if(iphotons.size() < 1 ) continue;
    hCounter->Fill(16);
    
    TLorentzVector pho;
    int j = iphotons[0];
    pho.SetPtEtaPhiM((*phoEt)[j], (*phoEta)[j], (*phoPhi)[j], 0);
    
    iPho1= j; 
    
    TLorentzVector ElecPhopair;
    ElecPhopair = pho+ele;
    InvMass=ElecPhopair.M();
    h_InvMass->Fill(InvMass,totalWeight);
    //cout<<"Zmass: "<<InvMass<<endl;    
    // if(InvMass>70. && InvMass<110)continue;
    hCounter->Fill(17);
    
    //    if(mass_w < 70.) continue;
    hCounter->Fill(18);
    
    //  if(phohasPixelSeed->at(j)!=0)continue;
    //   hCounter->Fill(18);
    //miniTree->Fill();
    //if(iphotons.size() < 1 ) continue;
    /*	bool HasMatchingGSFelectron=0;
      for(int ipho=0;ipho<iphotons.size();ipho++)
      {
      cout<<iphotons.size()<<" "<<ielectrons.size()<<endl;
      for (int iele=0; iele<ielectrons.size(); iele++){
      if (elePt->at(iele)<2) break;
      float delR = dR(phoSCEta->at(ipho),phoSCPhi->at(ipho),eleEta->at(iele),elePhi->at(iele));
      h_HasMatchingGSFelectron->Fill(delR,totalWeight);
      cout<<"???????????????????????????????????? dR"<<delR<<endl;
      if (delR<0.2) {HasMatchingGSFelectron= 1;break;}
      }
      //if(HasMatchingGSFelectron(ielectrons,ipho))continue;
      }
      
      if(HasMatchingGSFelectron==1)continue;
    */
    // OBTAIN WEIGHTS 
    double PhoIdweight = SF_mediumPhoID( (*phoEt)[j], (*phoSCEta)[j] )*PixelSeedVeto((*phoEt)[j], (*phoSCEta)[j]);
    //cout<<"pho weight: "<<PhoIdweight<<", eID weight: "<<SF_TightEID( (*elePt)[i], (*eleSCEta)[i] )<<endl;
    float scale = 1.0;
    scale *= SF_TightEID( (*elePt)[i], (*eleSCEta)[i] );
    scale *= PhoIdweight;
    
    IDWeight = scale;
    
    if ( processCode != 0 ) {
      totalWeight *= IDWeight;
      totalWeight_m5 *= IDWeight;
      totalWeight_p5 *= IDWeight;
      puWeight    *= IDWeight;
      puWeight_m5 *= IDWeight;
      puWeight_p5 *= IDWeight;
    }
   else {
     puWeight = 1.0;
     totalWeight = 1.0;
     puWeight_p5 = 1.0;
     totalWeight_p5 = 1.0;
     puWeight_m5 = 1.0;
     totalWeight_m5 = 1.0;
   }
    
    if(pfMET>30){pfmet30+=totalWeight;}
    
    if(mass_w > 50.) {WMt50+=totalWeight;}
    if(mass_w > 70.) {WMt70+=totalWeight;}
    
    
    //for blinding 
    //if ( processCode != 0 ) totalWeight  = totalWeight /20;
    
    pho_pt = pho.Pt();
    pho_eta = pho.Eta();
    pho_phi = pho.Phi();
    pho_SCEta = (*phoSCEta)[j];
    PFMET = pfMET;
    PFMETPhi = pfMETPhi;
    evt_weight = totalWeight;
    W_MT=mass_w;
    event = event;
    run = run;
    
    PhohasPixelSeed = phohasPixelSeed->at(j);
    pho_RandConeChIso = (*phoRandConeChIso)[j] - rho2012 * phoEffArea03ChHad((*phoSCEta)[j]);
    pho_ChHadIso = (*phoPFChIso)[j] - rho2012 * phoEffArea03ChHad((*phoSCEta)[j]);
    pho_SCRChIso = (*phoSCRChIso)[j] - rho2012 * phoEffArea03ChHad((*phoSCEta)[j]);
    pho_Sihih = (*phoSigmaIEtaIEta)[j];
    
    //check recopho match to reco jet
    float mindrRjRp=1111;
    float mindrRjRe=1111;
    float mindrRjRe2=1111;
    int jInd = -9;
    for(int j=0;j<nJet;j++)
      {
	float drRjRe = dR(ele.Eta(),ele.Phi(),(*jetEta)[j],(*jetPhi)[j]);
	if(drRjRe < mindrRjRe) mindrRjRe = drRjRe;
	//if(drRjRe<0.5)break;
	if(ele2.Pt()>0)
	  {float drRjRe2 = dR(ele2.Eta(),ele2.Phi(),(*jetEta)[j],(*jetPhi)[j]);
	    if(drRjRe2 < mindrRjRe2) mindrRjRe2 = drRjRe2;
	  }
	//if(drRjRe2<0.5)break;
	float drRjRp = dR(pho_eta,pho.Phi(),(*jetEta)[j],(*jetPhi)[j]);
	if(drRjRp < mindrRjRp) mindrRjRp = drRjRp;
	
	if(drRjRp<0.5)jInd = j;
	//cout<<jentry<<" #Jet: "<<nJet<<", pt: "<<(*jetPt)[j]<<", eta: "<<(*jetEta)[j]<<endl;
      }      
    
    MatchRecoEleRecoJet = mindrRjRe;
    MatchRecoEle2RecoJet = mindrRjRe2;
    MatchRecoPhoRecoJet = mindrRjRp;
    if(MatchRecoPhoRecoJet<0.5){RecoJetPt_drRjRp0p5 = (*jetPt)[jInd];}
    else RecoJetPt_drRjRp0p5 = -9;
    
    if(processCode!=0 ){
      
      //check recopho match for genpho
      double dr_min0 = 1000;
      double dr0 = 1111;
      double dr_min1 = 1000;
      double dr1 = 1111;
      int matchingGenPhoIndtoRecoPho = -9; int matchingGenPhoIndtoRecoEle = -9;
      for (int genPho = 0; genPho < nMC; ++genPho){
        if( abs((*mcPID)[genPho]) != 22  ) continue;
	//        if( (*mcStatus)[genPho] != 1 ) continue;
        //if( (*mcPt)[genPho] < 10 ) continue;
        dr0 = dR(ele.Eta(),ele.Phi(),(*mcEta)[genPho],(*mcPhi)[genPho]);
        if(dr0< dr_min0) {dr_min0=dr0; matchingGenPhoIndtoRecoEle = genPho;}

        dr1 = dR(pho_eta,pho.Phi(),(*mcEta)[genPho],(*mcPhi)[genPho]);
	if(dr1< dr_min1) {dr_min1=dr1; matchingGenPhoIndtoRecoPho = genPho;}
      }
      MatchRecoEleGenPho = dr_min0;
      MatchRecoPhoGenPho = dr_min1;
        if(matchingGenPhoIndtoRecoPho<0)
	  {	
	    genPartPt_MatchRecoPhoGenPho = -999;
	    genPartEta_MatchRecoPhoGenPho = -999;
	  }
        else
	  {
	    genPartPt_MatchRecoPhoGenPho = (*mcPt)[matchingGenPhoIndtoRecoPho];
	    genPartEta_MatchRecoPhoGenPho = (*mcEta)[matchingGenPhoIndtoRecoPho];
	  }
	
	if(matchingGenPhoIndtoRecoEle<0)
	  {
	    genPartPt_MatchRecoEleGenPho = -999;
	    genPartEta_MatchRecoEleGenPho = -999;
	  }
	else
	  {
	    genPartPt_MatchRecoEleGenPho = (*mcPt)[matchingGenPhoIndtoRecoEle];
	    genPartEta_MatchRecoEleGenPho = (*mcEta)[matchingGenPhoIndtoRecoEle];
	  }
	
	double dr_min2 = 1000;
	double dr2 = 1111;      
	
	double dr_min3 = 1000;
	double dr3 = 1111;
	int matchingGenIndtoRecoPho = -9;
	int matchingGenIndtoRecoEle = -9;
	for (int genPho = 0; genPho < nMC; ++genPho){
	  if( abs((*mcPID)[genPho]) != 11  ) continue;
	  if( (*mcMomPID)[genPho] !=23)continue;
	  //   if( (*mcStatus)[genPho] != 1 ) continue;
	  //      if( (*mcPt)[genPho] < 10 ) continue;
	  //check recoele match for genele
	  dr2 = dR(ele.Eta(),ele.Phi(),(*mcEta)[genPho],(*mcPhi)[genPho]);
	  if(dr2< dr_min2) {dr_min2=dr2; matchingGenIndtoRecoEle = genPho;}
	  //check recopho match for genele
	  dr3 = dR(pho_eta,pho.Phi(),(*mcEta)[genPho],(*mcPhi)[genPho]);
	  //	    cout<<dr3<<endl;
	  if(dr3 < dr_min3){dr_min3 = dr3; matchingGenIndtoRecoPho = genPho;}
	}
	
	
	//cout<<dr_min3<<endl; 
	MatchRecoPhoGenEle = dr_min3;
	MatchRecoEleGenEle = dr_min2;
	// if(MatchRecoPhoGenEle<0.2){
	if(matchingGenIndtoRecoPho<0){
	  genPartPt_MatchRecoPhoGenEle = -999;
	  genPartEta_MatchRecoPhoGenEle =-999;
	  genPartPdgId_MatchRecoPhoGenEle = -999;
	  genPartMomPdgId_MatchRecoPhoGenEle = -999;
	}	
	else {
	  genPartPt_MatchRecoPhoGenEle = (*mcPt)[matchingGenIndtoRecoPho];
	  genPartEta_MatchRecoPhoGenEle = (*mcEta)[matchingGenIndtoRecoPho];
	  genPartPdgId_MatchRecoPhoGenEle = (*mcPID)[matchingGenIndtoRecoPho];
	  genPartMomPdgId_MatchRecoPhoGenEle = (*mcMomPID)[matchingGenIndtoRecoPho];
	}
	if(matchingGenIndtoRecoEle<0){
	  genPartPt_MatchRecoEleGenEle = -999;
	  genPartEta_MatchRecoEleGenEle = -999;
	  genPartPdgId_MatchRecoEleGenEle = -999;
	  genPartMomPdgId_MatchRecoEleGenEle =-999;
	}
	else {
	  genPartPt_MatchRecoEleGenEle = (*mcPt)[matchingGenIndtoRecoEle];
	  genPartEta_MatchRecoEleGenEle = (*mcEta)[matchingGenIndtoRecoEle];
	  
	  genPartPdgId_MatchRecoEleGenEle = (*mcPID)[matchingGenIndtoRecoEle];
	  genPartMomPdgId_MatchRecoEleGenEle = (*mcMomPID)[matchingGenIndtoRecoEle];
	}
    	
	
    }
    WMT->Fill(mass_w,totalWeight);
    MET->Fill(pfMET,totalWeight);
    pT_gamma->Fill(pho.Pt(),totalWeight);
    pT_photon->Fill(pho.Pt(),totalWeight);
    h_elec_pho_dR->Fill(pho.DeltaR(ele),totalWeight);
    
    //    leftEvents +=totalWeight;
    leftEvents++;
    h_left->Fill(leftEvents);
    if(processCode==1){
      if(isWenu==1) countWenu++;
      else if(isWmunu==1) countWmunu++;  
      else if(isWtaunu==1) countWtaunu++;  
    }
    //if(isWenu==1)
    miniTree->Fill();
    
    //if(isWenu==1 && fabs((*phoSCEta)[j])>1.442 )
    //fout<<run<<'\t'<<event<<endl;//<<" isWenu: "<<isWenu<<" or isWtaunu: "<<isWtaunu
    //	<<" or isWmunu: "<<isWmunu<<endl;
    
  } // loop over for events
  miniTree->Write();
  cout<<"overlapping: "<<overlapping<<",  left: "<<leftEvents<<"( Wenu = "<<countWenu<<
    ", Wmunu = "<<countWmunu<<", Wtaunu = "<<countWtaunu<<" )"<<", pfMET30: "<<pfmet30<<", WMt50: "<<WMt50<<", WMt70: "<<WMt70<<endl;
  //cout<<leftEvents<<endl;
  tmp->Write();
  tmp->Close();
}

int main(int argc, char* argv[]) {
  TString fName = argv[1];
  
  analyze t(fName);
  //t.Loop("OnlyWenuUsingMCDecayType_miniTree_"+fName+".root");
  t.Loop("2ndEVetoNoIsominiTreewSFsNoSieieIsoMegPSVwMt_"+fName+".root");
  //t.Loop("miniTreewSFs_"+fName+".root");
  
  //t.Loop("tmp_scaleDown.root", 2);
  //t.Loop("tmp_res.root", 3);
  
  return 0;
}



			  
