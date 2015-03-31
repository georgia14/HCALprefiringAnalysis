#define l1ntuple_cxx
#include "L1Ntuple.h"

double L1Ntuple::deltaPhi(double phi1, double phi2) {

  if (phi1<0) phi1+=2.*TMath::Pi();
  if (phi2<0) phi2+=2.*TMath::Pi();

  double result = phi1 - phi2;
  if(fabs(result) > 9999) return result;
  while (result > TMath::Pi()) result -= 2*TMath::Pi();
  while (result <= -TMath::Pi()) result += 2*TMath::Pi();
  return result;
}

double L1Ntuple::deltaR(double eta1, double phi1, double eta2, double phi2) {
  double deta = eta1 - eta2;
  double dphi = deltaPhi(phi1, phi2);
  return sqrt(deta*deta + dphi*dphi);
}

int L1Ntuple::ietaBin(double eta) const {

  assert (eta > -5.191 && eta < 5.191);

 int etaBin = 10000;

 if( eta < -4.889 )      etaBin = -41;
 else if( eta < -4.716 ) etaBin = -40;
 else if( eta < -4.538 ) etaBin = -39;
 else if( eta < -4.363 ) etaBin = -38;
 else if( eta < -4.191 ) etaBin = -37;
 else if( eta < -4.013 ) etaBin = -36;
 else if( eta < -3.839 ) etaBin = -35;
 else if( eta < -3.664 ) etaBin = -34;
 else if( eta < -3.489 ) etaBin = -33;
 else if( eta < -3.314 ) etaBin = -32;
 else if( eta < -3.139 ) etaBin = -31;
 else if( eta < -2.964 ) etaBin = -30;
 else if( eta < -2.853 ) etaBin = -29; 
 else if( eta <  -2.65 ) etaBin = -28;
 else if( eta <   -2.5 ) etaBin = -27;
 else if( eta < -2.322 ) etaBin = -26;
 else if( eta < -2.172 ) etaBin = -25;
 else if( eta < -2.043 ) etaBin = -24;
 else if( eta <  -1.93 ) etaBin = -23;
 else if( eta <  -1.83 ) etaBin = -22;
 else if( eta <  -1.74 ) etaBin = -21;
 else if( eta < -1.653 ) etaBin = -20;
 else if( eta < -1.566 ) etaBin = -19;
 else if( eta < -1.479 ) etaBin = -18;
 else if( eta < -1.392 ) etaBin = -17;
 else if( eta < -1.305 ) etaBin = -16;
 else if( eta < -1.218 ) etaBin = -15;
 else if( eta < -1.131 ) etaBin = -14;
 else if( eta < -1.044 ) etaBin = -13;
 else if( eta < -0.957 ) etaBin = -12;
 else if( eta < -0.879 ) etaBin = -11;
 else if( eta < -0.783 ) etaBin = -10;
 else if( eta < -0.696 ) etaBin = -9;
 else if( eta < -0.609 ) etaBin = -8;
 else if( eta < -0.522 ) etaBin = -7;
 else if( eta < -0.435 ) etaBin = -6;
 else if( eta < -0.348 ) etaBin = -5;
 else if( eta < -0.261 ) etaBin = -4;
 else if( eta < -0.174 ) etaBin = -3;
 else if( eta < -0.087 ) etaBin = -2;
 else if( eta <      0 ) etaBin = -1;
 else if( eta <  0.087 ) etaBin =  1;
 else if( eta <  0.174 ) etaBin =  2;
 else if( eta <  0.261 ) etaBin =  3;
 else if( eta <  0.348 ) etaBin =  4;
 else if( eta <  0.435 ) etaBin =  5;
 else if( eta <  0.522 ) etaBin =  6;
 else if( eta <  0.609 ) etaBin =  7;
 else if( eta <  0.696 ) etaBin =  8;
 else if( eta <  0.783 ) etaBin =  9;
 else if( eta <  0.879 ) etaBin =  10;
 else if( eta <  0.957 ) etaBin =  11;
 else if( eta <  1.044 ) etaBin =  12;
 else if( eta <  1.131 ) etaBin =  13;
 else if( eta <  1.218 ) etaBin =  14;
 else if( eta <  1.305 ) etaBin =  15;
 else if( eta <  1.392 ) etaBin =  16;
 else if( eta <  1.479 ) etaBin =  17;
 else if( eta <  1.566 ) etaBin =  18;
 else if( eta <  1.653 ) etaBin =  19;
 else if( eta <   1.74 ) etaBin =  20;
 else if( eta <   1.83 ) etaBin =  21;
 else if( eta <   1.93 ) etaBin =  22;
 else if( eta <  2.043 ) etaBin =  23;
 else if( eta <  2.172 ) etaBin =  24;
 else if( eta <  2.322 ) etaBin =  25;
 else if( eta <    2.5 ) etaBin =  26;
 else if( eta <   2.65 ) etaBin =  27;
 else if( eta <  2.853 ) etaBin =  28;
 else if( eta <  2.964 ) etaBin =  29;
 else if( eta <  3.139 ) etaBin =  30;
 else if( eta <  3.314 ) etaBin =  31;
 else if( eta <  3.489 ) etaBin =  32;
 else if( eta <  3.664 ) etaBin =  33;
 else if( eta <  3.839 ) etaBin =  34;
 else if( eta <  4.013 ) etaBin =  35;
 else if( eta <  4.191 ) etaBin =  36;
 else if( eta <  4.363 ) etaBin =  37;
 else if( eta <  4.538 ) etaBin =  38;
 else if( eta <  4.716 ) etaBin =  39;
 else if( eta <  4.889 ) etaBin =  40;
 else if( eta <  5.191 ) etaBin =  41;
 
 return etaBin;

}

int L1Ntuple::iphiBin(double &phi) const
 {

   int iphi=10000;
  
   const static double dPhi =  2 * TMath::Pi() / 72;
   if(phi < 0) phi += 2 * TMath::Pi();
   iphi = (int)(phi / dPhi);
   iphi--; 
   if (iphi==0) iphi = 72;
   if (iphi==-1) iphi = 71;
   //   if(iphi > 72) iphi = 72; 

   //   if(useTowerCenterEtaPhi_) {
   phi = (iphi-0.5) * dPhi;
   // }
   assert(phi < 2 * TMath::Pi());
   assert(iphi > 0);
   assert(iphi <= 72);

   return iphi;

 }


double L1Ntuple::etaVal(int ieta) {

  double ietaBins=58.;
  double ietaBinsHF=3.;

  double etavl;

  if (abs(ieta)>=29) { // HF Jet
    etavl=(double)ieta*(0.2/ietaBinsHF);
  } else { // central Jet
    etavl=(double)ieta*(0.6/ietaBins);
  }

  return etavl;

}

double L1Ntuple::phiVal(int iphi) {

  double phiBins=72.;

  double phivl;
  phivl=double(iphi)*(2.*TMath::Pi()/phiBins);

  return phivl;

}

int L1Ntuple::matchUpgradedJet(double etaj, double phij) {

  int ipass=-1;

  // New PMTs DetId --> convert to eta/phi of L1jet
  int ieta, iphi;

  ieta=ietaBin(etaj);
  iphi=iphiBin(phij);

  /*
  if (iphi==43) { // Normally is iphi==11
    std::cout << "FwdJet (eta, phi)= (" << etaj << ", " << phij << 
       ") is (ieta, iphi)= (" << ieta << ", " << iphi << ") " << std::endl;
   }
  */
  //  if (ieta<-28 && iphi==39) { ipass=999; }
  if (ieta<-28 && iphi==43) { ipass=999; }

  return ipass;
}



void L1Ntuple::Test()
{ 

  if (fChain == 0)  return;
 
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  unsigned int nevents =0;

  std::cout << nentries << " events to process"<<std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    //fChain->GetEvent(jentry);
  
    nevents++;
    if (nevents<9)  //eight first events
      { 
    std::cout << "--------------------- Event "<<jentry<<" ---------------------"<<std::endl;

    //event_
    std::cout << "L1Tree         : event_->run = "<<event_->run<<std::endl;

    }
  }
   
}

void L1Ntuple::prefiringHCAL()
{ 

  if (fChain == 0)  return;
 
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  //unsigned int nevents =0;

  std::cout << nentries << " events to process"<<std::endl;

  TFile *total=new TFile("total.root","RECREATE");
  TDirectory *dir = total->mkdir("l1ExtraJets","Histograms for l1Extra Jet candidates; timing studies");

  // Central Jets
  TH1D *ptJet_early=new TH1D("ptJet_early","",200,0.,400.);
  TH1D *ptJet_central=new TH1D("ptJet_central","",200,0.,400.);
  TH1D *ptJet_late=new TH1D("ptJet_late","",200,0.,400.);

  TH1D *etaJet_early=new TH1D("etaJet_early","",30,-3.,3.);
  TH1D *etaJet_central=new TH1D("etaJet_central","",30,-3.,3.);
  TH1D *etaJet_late=new TH1D("etaJet_late","",30,-3.,3.);
 //  TH1D *etaJet_early=new TH1D("etaJet_early","",3000,-3.,3.);
//   TH1D *etaJet_central=new TH1D("etaJet_central","",3000,-3.,3.);
//   TH1D *etaJet_late=new TH1D("etaJet_late","",3000,-3.,3.);

  TH1D *phiJet_early=new TH1D("phiJet_early","",30,-TMath::Pi(),TMath::Pi());
  TH1D *phiJet_central=new TH1D("phiJet_central","",30,-TMath::Pi(),TMath::Pi());
  TH1D *phiJet_late=new TH1D("phiJet_late","",30,-TMath::Pi(),TMath::Pi());
 //  TH1D *phiJet_early=new TH1D("phiJet_early","",3000,-TMath::Pi(),TMath::Pi());
//   TH1D *phiJet_central=new TH1D("phiJet_central","",3000,-TMath::Pi(),TMath::Pi());
//   TH1D *phiJet_late=new TH1D("phiJet_late","",3000,-TMath::Pi(),TMath::Pi());

  // Debugging plots
  TH1D *etaJet_early_H=new TH1D("etaJet_early_H","",30,-3.,3.);
  TH1D *etaJet_central_H=new TH1D("etaJet_central_H","",30,-3.,3.);
  TH1D *phiJet_early_H=new TH1D("phiJet_early_H","",30,-TMath::Pi(),TMath::Pi());
  TH1D *phiJet_central_H=new TH1D("phiJet_central_H","",30,-TMath::Pi(),TMath::Pi());

  // HCAL TPs
  TH1D *ietaTP_H=new TH1D("ietaTP_H","",69,-34.5,34.5);
  TH1D *iphiTP_H=new TH1D("iphiTP_H","",77,-0.5,76.5);
  TH2D *ietaTP_iphiTP_H=new TH2D("ietaTP_iphiTP_H","",69,-34.5,34.5,77,-0.5,76.5);
  TH2D *ietaJet_iphiJet_H=new TH2D("ietaJet_iphiJet_H","",69,-34.5,34.5,77,-0.5,76.5);
  TH2D *ietaJet_iphiJet_central_H=new TH2D("ietaJet_iphiJet_central_H","",69,-34.5,34.5,77,-0.5,76.5);
  TH2D *ietaTP_iphiTP_H_2=new TH2D("ietaTP_iphiTP_H_2","",69,-34.5,34.5,77,-0.5,76.5);
  TH3D *ietaTP_iphiTP_etTP_H=new TH3D("ietaTP_iphiTP_etTP_H","",69,-34.5,34.5,77,-0.5,76.5,500,0.,50.);

  TH1D *etTP_H=new TH1D("etTP_H","",500,0.,50.);
  TH1D *compEtTP_H=new TH1D("compEtTP_H","",1500,0.,150.);

  // ECAL TPs
  TH1D *e_ietaTP_H=new TH1D("e_ietaTP_H","",69,-34.5,34.5);
  TH1D *e_iphiTP_H=new TH1D("e_iphiTP_H","",77,-0.5,76.5);
  TH2D *e_ietaTP_iphiTP_H=new TH2D("e_ietaTP_iphiTP_H","",69,-34.5,34.5,77,-0.5,76.5);
  TH2D *e_ietaJet_iphiJet_H=new TH2D("e_ietaJet_iphiJet_H","",69,-34.5,34.5,77,-0.5,76.5);
  TH2D *e_ietaJet_iphiJet_central_H=new TH2D("e_ietaJet_iphiJet_central_H","",69,-34.5,34.5,77,-0.5,76.5);
  TH2D *e_ietaTP_iphiTP_H_2=new TH2D("e_ietaTP_iphiTP_H_2","",69,-34.5,34.5,77,-0.5,76.5);
  TH3D *e_ietaTP_iphiTP_etTP_H=new TH3D("e_ietaTP_iphiTP_etTP_H","",69,-34.5,34.5,77,-0.5,76.5,500,0.,50.);

  TH1D *e_etTP_H=new TH1D("e_etTP_H","",500,0.,50.);
  TH1D *e_compEtTP_H=new TH1D("e_compEtTP_H","",1500,0.,150.);

  TH2D *ietaJet_iphiJet=new TH2D("ietaJet_iphiJet","",69,-34.5,34.5,77,-0.5,76.5);
  TH2D *ietaJet_iphiJet_central=new TH2D("ietaJet_iphiJet_central","",69,-34.5,34.5,77,-0.5,76.5);

  // HCAL TPs
  TH1D *ietaTP_Hs=new TH1D("ietaTP_Hs","",69,-34.5,34.5);
  TH1D *iphiTP_Hs=new TH1D("iphiTP_Hs","",77,-0.5,76.5);
  TH2D *ietaTP_iphiTP_Hs=new TH2D("ietaTP_iphiTP_Hs","",69,-34.5,34.5,77,-0.5,76.5);
  TH2D *ietaJet_iphiJet_Hs=new TH2D("ietaJet_iphiJet_Hs","",69,-34.5,34.5,77,-0.5,76.5);
  TH2D *ietaJet_iphiJet_central_Hs=new TH2D("ietaJet_iphiJet_central_Hs","",69,-34.5,34.5,77,-0.5,76.5);
  TH2D *ietaTP_iphiTP_Hs_2=new TH2D("ietaTP_iphiTP_Hs_2","",69,-34.5,34.5,77,-0.5,76.5);
  TH3D *ietaTP_iphiTP_etTP_Hs=new TH3D("ietaTP_iphiTP_etTP_Hs","",69,-34.5,34.5,77,-0.5,76.5,500,0.,50.);

  TH1D *etTP_Hs=new TH1D("etTP_Hs","",500,0.,50.);
  TH1D *compEtTP_Hs=new TH1D("compEtTP_Hs","",1500,0.,150.);

  // ECAL TPs
  TH1D *e_ietaTP_Hs=new TH1D("e_ietaTP_Hs","",69,-34.5,34.5);
  TH1D *e_iphiTP_Hs=new TH1D("e_iphiTP_Hs","",77,-0.5,76.5);
  TH2D *e_ietaTP_iphiTP_Hs=new TH2D("e_ietaTP_iphiTP_Hs","",69,-34.5,34.5,77,-0.5,76.5);
  TH2D *e_ietaJet_iphiJet_Hs=new TH2D("e_ietaJet_iphiJet_Hs","",69,-34.5,34.5,77,-0.5,76.5);
  TH2D *e_ietaJet_iphiJet_central_Hs=new TH2D("e_ietaJet_iphiJet_central_Hs","",69,-34.5,34.5,77,-0.5,76.5);
  TH2D *e_ietaTP_iphiTP_Hs_2=new TH2D("e_ietaTP_iphiTP_Hs_2","",69,-34.5,34.5,77,-0.5,76.5);
  TH3D *e_ietaTP_iphiTP_etTP_Hs=new TH3D("e_ietaTP_iphiTP_etTP_Hs","",69,-34.5,34.5,77,-0.5,76.5,500,0.,50.);

  TH1D *e_etTP_Hs=new TH1D("e_etTP_Hs","",500,0.,50.);
  TH1D *e_compEtTP_Hs=new TH1D("e_compEtTP_Hs","",1500,0.,150.);

  // HCAL TPs
  TH1D *ietaTP=new TH1D("ietaTP","",69,-34.5,34.5);
  TH1D *iphiTP=new TH1D("iphiTP","",77,-0.5,76.5);
  TH2D *ietaTP_iphiTP=new TH2D("ietaTP_iphiTP","",69,-34.5,34.5,77,-0.5,76.5);
  TH2D *ietaTP_iphiTP_2=new TH2D("ietaTP_iphiTP_2","",69,-34.5,34.5,77,-0.5,76.5);
  TH3D *ietaTP_iphiTP_etTP=new TH3D("ietaTP_iphiTP_etTP","",69,-34.5,34.5,77,-0.5,76.5,500,0.,50.);

  TH1D *etTP=new TH1D("etTP","",500,0.,50.);
  TH1D *compEtTP=new TH1D("compEtTP","",1500,0.,150.);

  // ECAL TPs
  TH1D *e_ietaTP=new TH1D("e_ietaTP","",69,-34.5,34.5);
  TH1D *e_iphiTP=new TH1D("e_iphiTP","",77,-0.5,76.5);
  TH2D *e_ietaTP_iphiTP=new TH2D("e_ietaTP_iphiTP","",69,-34.5,34.5,77,-0.5,76.5);
  TH2D *e_ietaTP_iphiTP_2=new TH2D("e_ietaTP_iphiTP_2","",69,-34.5,34.5,77,-0.5,76.5);
  TH3D *e_ietaTP_iphiTP_etTP=new TH3D("e_ietaTP_iphiTP_etTP","",69,-34.5,34.5,77,-0.5,76.5,500,0.,50.);

  TH1D *e_etTP=new TH1D("e_etTP","",500,0.,50.);
  TH1D *e_compEtTP=new TH1D("e_compEtTP","",1500,0.,150.);

  // HFJets
  TH1D *ptJetHF_early=new TH1D("ptJetHF_early","",200,0.,400.);
  TH1D *ptJetHF_central=new TH1D("ptJetHF_central","",200,0.,400.);
  TH1D *ptJetHF_late=new TH1D("ptJetHF_late","",200,0.,400.);

  TH1D *etaJetHF_early=new TH1D("etaJetHF_early","",50,-5.,5.);
  TH1D *etaJetHF_central=new TH1D("etaJetHF_central","",50,-5.,5.);
  TH1D *etaJetHF_late=new TH1D("etaJetHF_late","",50,-5.,5.);

  TH1D *phiJetHF_early=new TH1D("phiJetHF_early","",30,-TMath::Pi(),TMath::Pi());
  TH1D *phiJetHF_central=new TH1D("phiJetHF_central","",30,-TMath::Pi(),TMath::Pi());
  TH1D *phiJetHF_late=new TH1D("phiJetHF_late","",30,-TMath::Pi(),TMath::Pi());

  // Debugging plots
  TH1D *etaJetHF_early_H=new TH1D("etaJetHF_early_H","",50,-5.,5.);
  TH1D *etaJetHF_central_H=new TH1D("etaJetHF_central_H","",50,-5.,5.);
  TH1D *phiJetHF_early_H=new TH1D("phiJetHF_early_H","",30,-TMath::Pi(),TMath::Pi());
  TH1D *phiJetHF_central_H=new TH1D("phiJetHF_central_H","",30,-TMath::Pi(),TMath::Pi());

  // HCAL TPs
  TH1D *ietaTP_hf_H=new TH1D("ietaTP_hf_H","",69,-34.5,34.5);
  TH1D *iphiTP_hf_H=new TH1D("iphiTP_hf_H","",77,-0.5,76.5);
  TH2D *ietaTP_iphiTP_hf_H=new TH2D("ietaTP_iphiTP_hf_H","",69,-34.5,34.5,77,-0.5,76.5);
  TH2D *ietaJet_iphiJet_hf_H=new TH2D("ietaJet_iphiJet_hf_H","",69,-34.5,34.5,77,-0.5,76.5);
  TH2D *ietaJet_iphiJet_hf_central_H=new TH2D("ietaJet_iphiJet_hf_central_H","",69,-34.5,34.5,77,-0.5,76.5);
  TH2D *ietaTP_iphiTP_hf_H_2=new TH2D("ietaTP_iphiTP_hf_H_2","",69,-34.5,34.5,77,-0.5,76.5);
  TH3D *ietaTP_iphiTP_etTP_hf_H=new TH3D("ietaTP_iphiTP_etTP_hf_H","",69,-34.5,34.5,77,-0.5,76.5,500,0.,50.);

  TH1D *etTP_hf_H=new TH1D("etTP_hf_H","",500,0.,50.);
  TH1D *compEtTP_hf_H=new TH1D("compEtTP_hf_H","",1500,0.,150.);

  TH2D *ietaJet_iphiJet_hf=new TH2D("ietaJet_iphiJet_hf","",69,-34.5,34.5,77,-0.5,76.5);
  TH2D *ietaJet_iphiJet_hf_central=new TH2D("ietaJet_iphiJet_hf_central","",69,-34.5,34.5,77,-0.5,76.5);

  TH1D *ietaTP_hf_Hs=new TH1D("ietaTP_hf_Hs","",69,-34.5,34.5);
  TH1D *iphiTP_hf_Hs=new TH1D("iphiTP_hf_Hs","",77,-0.5,76.5);
  TH2D *ietaTP_iphiTP_hf_Hs=new TH2D("ietaTP_iphiTP_hf_Hs","",69,-34.5,34.5,77,-0.5,76.5);
  TH2D *ietaJet_iphiJet_hf_Hs=new TH2D("ietaJet_iphiJet_hf_Hs","",69,-34.5,34.5,77,-0.5,76.5);
  TH2D *ietaJet_iphiJet_hf_central_Hs=new TH2D("ietaJet_iphiJet_hf_central_Hs","",69,-34.5,34.5,77,-0.5,76.5);
  TH2D *ietaTP_iphiTP_hf_Hs_2=new TH2D("ietaTP_iphiTP_hf_Hs_2","",69,-34.5,34.5,77,-0.5,76.5);
  TH3D *ietaTP_iphiTP_etTP_hf_Hs=new TH3D("ietaTP_iphiTP_etTP_hf_Hs","",69,-34.5,34.5,77,-0.5,76.5,500,0.,50.);

  TH1D *etTP_hf_Hs=new TH1D("etTP_hf_Hs","",500,0.,50.);
  TH1D *compEtTP_hf_Hs=new TH1D("compEtTP_hf_Hs","",1500,0.,150.);

  Int_t runPeak=0;

  for (Long64_t jentry=0; jentry<nentries;jentry++)
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if(!(ientry%100000) && ientry) { std::cout << "processing event " << ientry << "\r" << std::endl; }

    //    if (jentry!=2492259 && event_->run!=205666) continue;

    Int_t ie,ic,il;
    ie=ic=il=0;

    Int_t ie1,ic1,il1;
    ie1=ic1=il1=-1;

    // Write your code here
    
    // Loop over central Jets in l1extra
    for (Int_t i=0; i<l1extra_->nCenJets; i++) {
      // eaarly candidates 
      if (l1extra_->cenJetBx[i]==-1 ) {
	ie++; if(ie>1) continue;
	ie1=i;
      }
      // central cadndidates
      if (l1extra_->cenJetBx[i]==0) {
	ic++; if (ic>1) continue;
	ic1=i;
      }
      // eaarly candidates
      if (l1extra_->cenJetBx[i]==1 ) {
	il++; if (il>1) continue;
	il1=i;
      }
    }

    // Run Analysis for L1Jets only that fall into the new PMTs detIds.
    //    Double_t isUpgraded=-1.;
    Double_t yHBHE_early, yHBHE_central;
    yHBHE_early=yHBHE_central=-1.;

    // For debugging
    Double_t idH_early, idH_central;
    idH_early=idH_central=-1.;

    Double_t etaj_early, phij_early;
    Double_t etaj_central, phij_central;

    etaj_early=etaj_central=9999.;
    phij_early=phij_central=9999.;

      // Loop over tau Jets in l1extra
    ie=ic=il=0;

    for (Int_t i=0; i<l1extra_->nTauJets; i++) {
	// early candidates
      if (l1extra_->tauJetBx[i]==-1 ) {
	ie++; if (ie>1) continue;

	yHBHE_early=999.;

	if (ie1>=0) {
	  if (l1extra_->tauJetEt[i]>=l1extra_->cenJetEt[ie1]) {
	    ptJet_early->Fill(l1extra_->tauJetEt[i]);
	    etaJet_early->Fill(l1extra_->tauJetEta[i]);
	    phiJet_early->Fill(l1extra_->tauJetPhi[i]);

	    etaj_early=l1extra_->tauJetEta[i];
	    phij_early=l1extra_->tauJetPhi[i];
	    
	    // Debugging
	    if (l1extra_->tauJetEt[i]>=240.) {
	      etaJet_early_H->Fill(l1extra_->tauJetEta[i]);
	      phiJet_early_H->Fill(l1extra_->tauJetPhi[i]);
	      idH_early=999.;
	    }

	  } else {
	    ptJet_early->Fill(l1extra_->cenJetEt[ie1]);
	    etaJet_early->Fill(l1extra_->cenJetEta[ie1]);
	    phiJet_early->Fill(l1extra_->cenJetPhi[ie1]);

	    etaj_early=l1extra_->cenJetEta[ie1];
	    phij_early=l1extra_->cenJetPhi[ie1];

	    if (l1extra_->cenJetEt[ie1]>=240.) {
	      etaJet_early_H->Fill(l1extra_->cenJetEta[ie1]);
	      phiJet_early_H->Fill(l1extra_->cenJetPhi[ie1]);
	      idH_early=999.;
	    }

	  }
	} else {
	  ptJet_early->Fill(l1extra_->tauJetEt[i]);
	  etaJet_early->Fill(l1extra_->tauJetEta[i]);
	  phiJet_early->Fill(l1extra_->tauJetPhi[i]);
	  
	  etaj_early=l1extra_->tauJetEta[i];
	  phij_early=l1extra_->tauJetPhi[i];

	  if (l1extra_->tauJetEt[i]>=240.) {
	    etaJet_early_H->Fill(l1extra_->tauJetEta[i]);
	    phiJet_early_H->Fill(l1extra_->tauJetPhi[i]);
	    idH_early=999.;
	  }

	}
      }
      // central candidates
      if (l1extra_->tauJetBx[i]==0) {
	ic++; if (ic>1) continue;

	yHBHE_central=999.;

	if (ic1>=0) {
	  if (l1extra_->tauJetEt[i]>=l1extra_->cenJetEt[ic1]) {
	    ptJet_central->Fill(l1extra_->tauJetEt[i]);
	    etaJet_central->Fill(l1extra_->tauJetEta[i]);
	    phiJet_central->Fill(l1extra_->tauJetPhi[i]);

	    etaj_central=l1extra_->tauJetEta[i];
	    phij_central=l1extra_->tauJetPhi[i];

	    if (l1extra_->tauJetEt[i]>=240.) {
              etaJet_central_H->Fill(l1extra_->tauJetEta[i]);
              phiJet_central_H->Fill(l1extra_->tauJetPhi[i]);
	      idH_central=999.;
            }

	  } else {
	    ptJet_central->Fill(l1extra_->cenJetEt[ic1]);
	    etaJet_central->Fill(l1extra_->cenJetEta[ic1]);
	    phiJet_central->Fill(l1extra_->cenJetPhi[ic1]);

	    etaj_central=l1extra_->cenJetEta[ic1];
	    phij_central=l1extra_->cenJetPhi[ic1];

	    if (l1extra_->cenJetEt[ie1]>=240.) {
              etaJet_central_H->Fill(l1extra_->cenJetEta[ie1]);
              phiJet_central_H->Fill(l1extra_->cenJetPhi[ie1]);
	      idH_central=999.;
            }


	  }
	} else {
	  ptJet_central->Fill(l1extra_->tauJetEt[i]);
	  etaJet_central->Fill(l1extra_->tauJetEta[i]);
	  phiJet_central->Fill(l1extra_->tauJetPhi[i]);

	  etaj_central=l1extra_->tauJetEta[i];
	  phij_central=l1extra_->tauJetPhi[i];

	  if (l1extra_->tauJetEt[i]>=240.) {
	    etaJet_central_H->Fill(l1extra_->tauJetEta[i]);
	    phiJet_central_H->Fill(l1extra_->tauJetPhi[i]);
	    idH_central=999.;
	  }

	}
      }
      // late candidates
      if (l1extra_->tauJetBx[i]==1 ) {
	il++; if (il>1) continue;
	if (il1>=0) {
	  if (l1extra_->tauJetEt[i]>=l1extra_->cenJetEt[il1]) {
	    ptJet_late->Fill(l1extra_->tauJetEt[i]);
	    etaJet_late->Fill(l1extra_->tauJetEta[i]);
	    phiJet_late->Fill(l1extra_->tauJetPhi[i]);
	  } else {
	    ptJet_late->Fill(l1extra_->cenJetEt[il1]);
	    etaJet_late->Fill(l1extra_->cenJetEta[il1]);
	    phiJet_late->Fill(l1extra_->cenJetPhi[il1]);
	  }
	} else {
	  ptJet_late->Fill(l1extra_->tauJetEt[i]);
	  etaJet_late->Fill(l1extra_->tauJetEta[i]);
	  phiJet_late->Fill(l1extra_->tauJetPhi[i]);
	}
      }
   
    }

      // All TP info
    for (Int_t j=0; j<caloTP_->nHCALTP; j++) {
      
      ietaTP->Fill(caloTP_->hcalTPieta[j]);
      iphiTP->Fill(caloTP_->hcalTPCaliphi[j]);
      ietaTP_iphiTP->Fill(caloTP_->hcalTPieta[j],caloTP_->hcalTPCaliphi[j]);
      ietaTP_iphiTP_2->Fill(caloTP_->hcalTPieta[j],caloTP_->hcalTPCaliphi[j],caloTP_->hcalTPet[j]);
      ietaTP_iphiTP_etTP->Fill(caloTP_->hcalTPieta[j],caloTP_->hcalTPCaliphi[j],caloTP_->hcalTPet[j]);
      
      etTP->Fill(caloTP_->hcalTPet[j]);
      compEtTP->Fill(caloTP_->hcalTPcompEt[j]);
    }

    for (Int_t j=0; j<caloTP_->nECALTP; j++) {
      
      e_ietaTP->Fill(caloTP_->ecalTPieta[j]);
      e_iphiTP->Fill(caloTP_->ecalTPCaliphi[j]);
      e_ietaTP_iphiTP->Fill(caloTP_->ecalTPieta[j],caloTP_->ecalTPCaliphi[j]);
      e_ietaTP_iphiTP_2->Fill(caloTP_->ecalTPieta[j],caloTP_->ecalTPCaliphi[j],caloTP_->ecalTPet[j]);
      e_ietaTP_iphiTP_etTP->Fill(caloTP_->ecalTPieta[j],caloTP_->ecalTPCaliphi[j],caloTP_->ecalTPet[j]);
      
      e_etTP->Fill(caloTP_->ecalTPet[j]);
      e_compEtTP->Fill(caloTP_->ecalTPcompEt[j]);
    }
    
      // ieta/iphi of BX=-1 and BX=0 L1Jets
    if (etaj_early<9999.) { ietaJet_iphiJet->Fill(ietaBin(etaj_early),iphiBin(phij_early)); }
    if (etaj_central<9999.) { ietaJet_iphiJet_central->Fill(ietaBin(etaj_central),iphiBin(phij_central)); }
      
      
      // If L1Jet BX=-1 and saturated ET
    if (idH_early>0) {
      /*	
      std::cout << " " << std::endl;
      std::cout << "Found event with early Jet of saturated Et!!! " << std::endl;
      std::cout << "--------------------- Event "<<jentry<<" ---------------------"<<std::endl;
      std::cout << "L1Tree         : event_->run = "<<event_->run<<std::endl;
      */
       if (etaj_early<9999.) ietaJet_iphiJet_Hs->Fill(ietaBin(etaj_early),iphiBin(phij_early));
       if (etaj_central<9999.) ietaJet_iphiJet_central_Hs->Fill(ietaBin(etaj_central),iphiBin(phij_central));
	  
       for (Int_t j=0; j<caloTP_->nHCALTP; j++) {
	 
	 ietaTP_Hs->Fill(caloTP_->hcalTPieta[j]);
	 iphiTP_Hs->Fill(caloTP_->hcalTPCaliphi[j]);
	 ietaTP_iphiTP_Hs->Fill(caloTP_->hcalTPieta[j],caloTP_->hcalTPCaliphi[j]);
	  
	 ietaTP_iphiTP_Hs_2->Fill(caloTP_->hcalTPieta[j],caloTP_->hcalTPCaliphi[j],caloTP_->hcalTPet[j]);
	 ietaTP_iphiTP_etTP_Hs->Fill(caloTP_->hcalTPieta[j],caloTP_->hcalTPCaliphi[j],caloTP_->hcalTPet[j]);
	 
	 etTP_Hs->Fill(caloTP_->hcalTPet[j]);
	 compEtTP_Hs->Fill(caloTP_->hcalTPcompEt[j]);
       }

       for (Int_t j=0; j<caloTP_->nECALTP; j++) {
	 
	 e_ietaTP_Hs->Fill(caloTP_->ecalTPieta[j]);
	 e_iphiTP_Hs->Fill(caloTP_->ecalTPCaliphi[j]);
	 e_ietaTP_iphiTP_Hs->Fill(caloTP_->ecalTPieta[j],caloTP_->ecalTPCaliphi[j]);
	  
	 e_ietaTP_iphiTP_Hs_2->Fill(caloTP_->ecalTPieta[j],caloTP_->ecalTPCaliphi[j],caloTP_->ecalTPet[j]);
	 e_ietaTP_iphiTP_etTP_Hs->Fill(caloTP_->ecalTPieta[j],caloTP_->ecalTPCaliphi[j],caloTP_->ecalTPet[j]);
	 
	 e_etTP_Hs->Fill(caloTP_->ecalTPet[j]);
	 e_compEtTP_Hs->Fill(caloTP_->ecalTPcompEt[j]);
       }

    }
		
    if (idH_central>0) { // ET of central Jet is saturated
      /*
      std::cout << " " << std::endl;
      std::cout << "Found event with early Jet of saturated Et!!! " << std::endl;
      std::cout << "--------------------- Event "<<jentry<<" ---------------------"<<std::endl;
      std::cout << "L1Tree : event_->run = "<<event_->run<<std::endl;
      */
      runPeak++; //=999;

      // ieta, iphi of the Jet
      if (etaj_early<9999.) ietaJet_iphiJet_H->Fill(ietaBin(etaj_early),iphiBin(phij_early));
      if (etaj_central<9999.) ietaJet_iphiJet_central_H->Fill(ietaBin(etaj_central),iphiBin(phij_central));
      
      // TP info
      for (Int_t j=0; j<caloTP_->nHCALTP; j++) {

	ietaTP_H->Fill(caloTP_->hcalTPieta[j]);
	iphiTP_H->Fill(caloTP_->hcalTPCaliphi[j]);
	ietaTP_iphiTP_H->Fill(caloTP_->hcalTPieta[j],caloTP_->hcalTPCaliphi[j]);
	
	ietaTP_iphiTP_H_2->Fill(caloTP_->hcalTPieta[j],caloTP_->hcalTPCaliphi[j],caloTP_->hcalTPet[j]);
	ietaTP_iphiTP_etTP_H->Fill(caloTP_->hcalTPieta[j],caloTP_->hcalTPCaliphi[j],caloTP_->hcalTPet[j]);
	
	etTP_H->Fill(caloTP_->hcalTPet[j]);
	compEtTP_H->Fill(caloTP_->hcalTPcompEt[j]);
	
      }

      for (Int_t j=0; j<caloTP_->nECALTP; j++) {

	e_ietaTP_H->Fill(caloTP_->ecalTPieta[j]);
	e_iphiTP_H->Fill(caloTP_->ecalTPCaliphi[j]);
	e_ietaTP_iphiTP_H->Fill(caloTP_->ecalTPieta[j],caloTP_->ecalTPCaliphi[j]);
	
	e_ietaTP_iphiTP_H_2->Fill(caloTP_->ecalTPieta[j],caloTP_->ecalTPCaliphi[j],caloTP_->ecalTPet[j]);
	e_ietaTP_iphiTP_etTP_H->Fill(caloTP_->ecalTPieta[j],caloTP_->ecalTPCaliphi[j],caloTP_->ecalTPet[j]);
	
	e_etTP_H->Fill(caloTP_->ecalTPet[j]);
	e_compEtTP_H->Fill(caloTP_->ecalTPcompEt[j]);
	
      }
    }

    // HFjets
    Double_t etaj_hf_early, phij_hf_early;
    Double_t etaj_hf_central, phij_hf_central;

    etaj_hf_early=etaj_hf_central=9999.;
    phij_hf_early=phij_hf_central=9999.;

    // For debuggig
    Double_t idHf_early, idHf_central;  // if HF jet at ET->250 GeV
    idHf_early=idHf_central=-1.;

    Double_t yHF_early, yHF_central; // if HF jet
    yHF_early=yHF_central=-1.;

    //    Int_t ie,ic,il;
    ie=ic=il=0;
    for (Int_t i=0; i<l1extra_->nFwdJets; i++) {
      // eaarly candidates
      if (l1extra_->fwdJetBx[i]==-1 ) {

	//	if (matchUpgradedJet(l1extra_->fwdJetEta[i],l1extra_->fwdJetPhi[i])<0) continue;
	ie++; if (ie>1) continue;

	yHF_early=999.;

        ptJetHF_early->Fill(l1extra_->fwdJetEt[i]);
        etaJetHF_early->Fill(l1extra_->fwdJetEta[i]);
        phiJetHF_early->Fill(l1extra_->fwdJetPhi[i]);

	etaj_hf_early=l1extra_->fwdJetEta[i];
	phij_hf_early=l1extra_->fwdJetPhi[i];

	if (l1extra_->fwdJetEt[i]>=240.) {
	  etaJetHF_early_H->Fill(l1extra_->fwdJetEta[i]);
	  phiJetHF_early_H->Fill(l1extra_->fwdJetPhi[i]);
	  idHf_early=999.;
	}

      }
      // fwdtral cadndidates
      if (l1extra_->fwdJetBx[i]==0) {

	//	if (matchUpgradedJet(l1extra_->fwdJetEta[i],l1extra_->fwdJetPhi[i])<0) continue;
	ic++; if (ic>1) continue;

	yHF_central=999.;

        ptJetHF_central->Fill(l1extra_->fwdJetEt[i]);
        etaJetHF_central->Fill(l1extra_->fwdJetEta[i]);
        phiJetHF_central->Fill(l1extra_->fwdJetPhi[i]);

	etaj_hf_central=l1extra_->fwdJetEta[i];
	phij_hf_central=l1extra_->fwdJetPhi[i];

	if (l1extra_->fwdJetEt[i]>=240.) {
          etaJetHF_central_H->Fill(l1extra_->fwdJetEta[i]);
          phiJetHF_central_H->Fill(l1extra_->fwdJetPhi[i]);
	  idHf_central=999.;
        }

      }
      // late candidates
      if (l1extra_->fwdJetBx[i]==1 ) {
	il++; if (il>1) continue;
        ptJetHF_late->Fill(l1extra_->fwdJetEt[i]);
        etaJetHF_late->Fill(l1extra_->fwdJetEta[i]);
        phiJetHF_late->Fill(l1extra_->fwdJetPhi[i]);
      }
    }

     // ieta/iphi of BX=-1 and BX=0 L1Jets
    if (etaj_hf_early<9999.) { ietaJet_iphiJet_hf->Fill(ietaBin(etaj_hf_early),iphiBin(phij_hf_early)); }
    if (etaj_hf_central<9999.) { ietaJet_iphiJet_hf_central->Fill(ietaBin(etaj_hf_central),iphiBin(phij_hf_central)); }


    //If found a central forward Jet
    //    if (yHF_central<0) continue;

    if (idHf_early>0) {

       if (etaj_hf_early<9999.) ietaJet_iphiJet_hf_Hs->Fill(ietaBin(etaj_hf_early),iphiBin(phij_hf_early));
       if (etaj_hf_central<9999.) ietaJet_iphiJet_hf_central_Hs->Fill(ietaBin(etaj_hf_central),iphiBin(phij_hf_central));
      
      for (Int_t j=0; j<caloTP_->nHCALTP; j++) {
	
	ietaTP_hf_Hs->Fill(caloTP_->hcalTPieta[j]);
	iphiTP_hf_Hs->Fill(caloTP_->hcalTPCaliphi[j]);
	ietaTP_iphiTP_hf_Hs->Fill(caloTP_->hcalTPieta[j],caloTP_->hcalTPCaliphi[j]);
	
	ietaTP_iphiTP_hf_Hs_2->Fill(caloTP_->hcalTPieta[j],caloTP_->hcalTPCaliphi[j],caloTP_->hcalTPet[j]);
	ietaTP_iphiTP_etTP_hf_Hs->Fill(caloTP_->hcalTPieta[j],caloTP_->hcalTPCaliphi[j],caloTP_->hcalTPet[j]);
	
	etTP_hf_Hs->Fill(caloTP_->hcalTPet[j]);
	compEtTP_hf_Hs->Fill(caloTP_->hcalTPcompEt[j]);
	
      }
    }

    if (idHf_central>0) {
	// ieta, iphi of the Jet
      if (etaj_hf_early<9999.) ietaJet_iphiJet_hf_H->Fill(ietaBin(etaj_hf_early),iphiBin(phij_hf_early));
      if (etaj_hf_central<9999.) ietaJet_iphiJet_hf_central_H->Fill(ietaBin(etaj_hf_central),iphiBin(phij_hf_central));
      
	// TP info
      for (Int_t j=0; j<caloTP_->nHCALTP; j++) {
	
	ietaTP_hf_H->Fill(caloTP_->hcalTPieta[j]);
	iphiTP_hf_H->Fill(caloTP_->hcalTPCaliphi[j]);
	ietaTP_iphiTP_hf_H->Fill(caloTP_->hcalTPieta[j],caloTP_->hcalTPCaliphi[j]);
	
	ietaTP_iphiTP_hf_H_2->Fill(caloTP_->hcalTPieta[j],caloTP_->hcalTPCaliphi[j],caloTP_->hcalTPet[j]);
	ietaTP_iphiTP_etTP_hf_H->Fill(caloTP_->hcalTPieta[j],caloTP_->hcalTPCaliphi[j],caloTP_->hcalTPet[j]);
	
	etTP_hf_H->Fill(caloTP_->hcalTPet[j]);
	compEtTP_hf_H->Fill(caloTP_->hcalTPcompEt[j]);
      }
      
    }

    //    if (runPeak>0) break;
  } // End loop over all Events

  ptJet_early->Sumw2(); ptJet_early->SetDirectory(dir);
  ptJet_central->Sumw2(); ptJet_central->SetDirectory(dir);
  ptJet_late->Sumw2(); ptJet_late->SetDirectory(dir);
  
  etaJet_early->Sumw2(); etaJet_early->SetDirectory(dir);
  etaJet_central->Sumw2(); etaJet_central->SetDirectory(dir);
  etaJet_late->Sumw2(); etaJet_late->SetDirectory(dir);

  phiJet_early->Sumw2(); phiJet_early->SetDirectory(dir);
  phiJet_central->Sumw2(); phiJet_central->SetDirectory(dir);
  phiJet_late->Sumw2(); phiJet_late->SetDirectory(dir);

  etaJet_early_H->Sumw2(); etaJet_early_H->SetDirectory(dir);
  etaJet_central_H->Sumw2(); etaJet_central_H->SetDirectory(dir);
  phiJet_early_H->Sumw2(); phiJet_early_H->SetDirectory(dir);
  phiJet_central_H->Sumw2(); phiJet_central_H->SetDirectory(dir);


  // HF jets
  ptJetHF_early->Sumw2(); ptJetHF_early->SetDirectory(dir);
  ptJetHF_central->Sumw2(); ptJetHF_central->SetDirectory(dir);
  ptJetHF_late->Sumw2(); ptJetHF_late->SetDirectory(dir);

  etaJetHF_early->Sumw2(); etaJetHF_early->SetDirectory(dir);
  etaJetHF_central->Sumw2(); etaJetHF_central->SetDirectory(dir);
  etaJetHF_late->Sumw2(); etaJetHF_late->SetDirectory(dir);

  phiJetHF_early->Sumw2(); phiJetHF_early->SetDirectory(dir);
  phiJetHF_central->Sumw2(); phiJetHF_central->SetDirectory(dir);
  phiJetHF_late->Sumw2(); phiJetHF_late->SetDirectory(dir);
  
  etaJetHF_early_H->Sumw2(); etaJetHF_early_H->SetDirectory(dir);
  etaJetHF_central_H->Sumw2(); etaJetHF_central_H->SetDirectory(dir);
  phiJetHF_early_H->Sumw2(); phiJetHF_early_H->SetDirectory(dir);
  phiJetHF_central_H->Sumw2(); phiJetHF_central_H->SetDirectory(dir);

  // TP info
  ietaTP_H->Sumw2(); ietaTP_H->SetDirectory(dir);
  iphiTP_H->Sumw2(); iphiTP_H->SetDirectory(dir);
  ietaTP_iphiTP_H->Sumw2(); ietaTP_iphiTP_H->SetDirectory(dir);
  ietaJet_iphiJet_H->Sumw2(); ietaJet_iphiJet_H->SetDirectory(dir);
  ietaJet_iphiJet_central_H->Sumw2(); ietaJet_iphiJet_central_H->SetDirectory(dir);

  e_ietaTP_H->Sumw2(); e_ietaTP_H->SetDirectory(dir);
  e_iphiTP_H->Sumw2(); e_iphiTP_H->SetDirectory(dir);
  e_ietaTP_iphiTP_H->Sumw2(); e_ietaTP_iphiTP_H->SetDirectory(dir);
  e_ietaJet_iphiJet_H->Sumw2(); e_ietaJet_iphiJet_H->SetDirectory(dir);
  e_ietaJet_iphiJet_central_H->Sumw2(); e_ietaJet_iphiJet_central_H->SetDirectory(dir);

  ietaJet_iphiJet->Sumw2(); ietaJet_iphiJet->SetDirectory(dir);
  ietaJet_iphiJet_central->Sumw2(); ietaJet_iphiJet_central->SetDirectory(dir);

  ietaTP_iphiTP_H_2->Sumw2(); ietaTP_iphiTP_H_2->SetDirectory(dir);
  ietaTP_iphiTP_etTP_H->Sumw2(); ietaTP_iphiTP_etTP_H->SetDirectory(dir);
  etTP_H->Sumw2(); etTP_H->SetDirectory(dir);
  compEtTP_H->Sumw2(); compEtTP_H->SetDirectory(dir);

  ietaTP_Hs->Sumw2(); ietaTP_Hs->SetDirectory(dir);
  iphiTP_Hs->Sumw2(); iphiTP_Hs->SetDirectory(dir);
  ietaTP_iphiTP_Hs->Sumw2(); ietaTP_iphiTP_Hs->SetDirectory(dir);

  e_ietaTP_iphiTP_H_2->Sumw2(); e_ietaTP_iphiTP_H_2->SetDirectory(dir);
  e_ietaTP_iphiTP_etTP_H->Sumw2(); e_ietaTP_iphiTP_etTP_H->SetDirectory(dir);
  e_etTP_H->Sumw2(); e_etTP_H->SetDirectory(dir);
  e_compEtTP_H->Sumw2(); e_compEtTP_H->SetDirectory(dir);

  e_ietaTP_Hs->Sumw2(); e_ietaTP_Hs->SetDirectory(dir);
  e_iphiTP_Hs->Sumw2(); e_iphiTP_Hs->SetDirectory(dir);
  e_ietaTP_iphiTP_Hs->Sumw2(); e_ietaTP_iphiTP_Hs->SetDirectory(dir);

  ietaJet_iphiJet_Hs->Sumw2(); ietaJet_iphiJet_Hs->SetDirectory(dir);
  ietaJet_iphiJet_central_Hs->Sumw2(); ietaJet_iphiJet_central_Hs->SetDirectory(dir);

  ietaTP_iphiTP_Hs_2->Sumw2(); ietaTP_iphiTP_Hs_2->SetDirectory(dir);
  ietaTP_iphiTP_etTP_Hs->Sumw2(); ietaTP_iphiTP_etTP_Hs->SetDirectory(dir);
  etTP_Hs->Sumw2(); etTP_Hs->SetDirectory(dir);
  compEtTP_Hs->Sumw2(); compEtTP_Hs->SetDirectory(dir);

  ietaTP->Sumw2(); ietaTP->SetDirectory(dir);
  iphiTP->Sumw2(); iphiTP->SetDirectory(dir);
  ietaTP_iphiTP->Sumw2(); ietaTP_iphiTP->SetDirectory(dir);
  ietaTP_iphiTP_2->Sumw2(); ietaTP_iphiTP_2->SetDirectory(dir);
  ietaTP_iphiTP_etTP->Sumw2(); ietaTP_iphiTP_etTP->SetDirectory(dir);
  etTP->Sumw2(); etTP->SetDirectory(dir);
  compEtTP->Sumw2(); compEtTP->SetDirectory(dir);

  e_ietaTP_iphiTP_Hs_2->Sumw2(); e_ietaTP_iphiTP_Hs_2->SetDirectory(dir);
  e_ietaTP_iphiTP_etTP_Hs->Sumw2(); e_ietaTP_iphiTP_etTP_Hs->SetDirectory(dir);
  e_etTP_Hs->Sumw2(); e_etTP_Hs->SetDirectory(dir);
  e_compEtTP_Hs->Sumw2(); e_compEtTP_Hs->SetDirectory(dir);

  e_ietaTP->Sumw2(); e_ietaTP->SetDirectory(dir);
  e_iphiTP->Sumw2(); e_iphiTP->SetDirectory(dir);
  e_ietaTP_iphiTP->Sumw2(); e_ietaTP_iphiTP->SetDirectory(dir);
  e_ietaTP_iphiTP_2->Sumw2(); e_ietaTP_iphiTP_2->SetDirectory(dir);
  e_ietaTP_iphiTP_etTP->Sumw2(); e_ietaTP_iphiTP_etTP->SetDirectory(dir);
  e_etTP->Sumw2(); e_etTP->SetDirectory(dir);
  e_compEtTP->Sumw2(); e_compEtTP->SetDirectory(dir);

  ietaTP_hf_H->Sumw2(); ietaTP_hf_H->SetDirectory(dir);
  iphiTP_hf_H->Sumw2(); iphiTP_hf_H->SetDirectory(dir);
  ietaTP_iphiTP_hf_H->Sumw2(); ietaTP_iphiTP_hf_H->SetDirectory(dir);
  ietaJet_iphiJet_hf_H->Sumw2(); ietaJet_iphiJet_hf_H->SetDirectory(dir);
  ietaJet_iphiJet_hf_central_H->Sumw2(); ietaJet_iphiJet_hf_central_H->SetDirectory(dir);

  ietaTP_iphiTP_hf_H_2->Sumw2(); ietaTP_iphiTP_hf_H_2->SetDirectory(dir);
  ietaTP_iphiTP_etTP_hf_H->Sumw2(); ietaTP_iphiTP_etTP_hf_H->SetDirectory(dir);
  etTP_hf_H->Sumw2(); etTP_hf_H->SetDirectory(dir);
  compEtTP_hf_H->Sumw2(); compEtTP_hf_H->SetDirectory(dir);

  ietaTP_hf_Hs->Sumw2(); ietaTP_hf_Hs->SetDirectory(dir);
  iphiTP_hf_Hs->Sumw2(); iphiTP_hf_Hs->SetDirectory(dir);
  ietaTP_iphiTP_hf_Hs->Sumw2(); ietaTP_iphiTP_hf_Hs->SetDirectory(dir);
  ietaJet_iphiJet_hf_Hs->Sumw2(); ietaJet_iphiJet_hf_Hs->SetDirectory(dir);
  ietaJet_iphiJet_hf_central_Hs->Sumw2(); ietaJet_iphiJet_hf_central_Hs->SetDirectory(dir);

  ietaTP_iphiTP_hf_Hs_2->Sumw2(); ietaTP_iphiTP_hf_Hs_2->SetDirectory(dir);
  ietaTP_iphiTP_etTP_hf_Hs->Sumw2(); ietaTP_iphiTP_etTP_hf_Hs->SetDirectory(dir);
  etTP_hf_Hs->Sumw2(); etTP_hf_Hs->SetDirectory(dir);
  compEtTP_hf_Hs->Sumw2(); compEtTP_hf_Hs->SetDirectory(dir);

  ietaJet_iphiJet_hf->Sumw2(); ietaJet_iphiJet_hf->SetDirectory(dir);
  ietaJet_iphiJet_hf_central->Sumw2(); ietaJet_iphiJet_hf_central->SetDirectory(dir);

  total->Write(); total->Close();

}

