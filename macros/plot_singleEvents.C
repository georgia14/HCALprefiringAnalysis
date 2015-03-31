#include "setTDRStyle.C"

TH1D* readHist(TString nameHist,TString nameFile, int rebin)
{
 TFile* file = new TFile(nameFile);

 TH1D* hist = (TH1D*)file->Get(nameHist);
 hist->GetSumw2();
 // hist->SetLineWidth(2);
 if(rebin>0) hist->Rebin(rebin);
 hist->GetXaxis()->SetTitleSize(.055);
 hist->GetYaxis()->SetTitleSize(.055);
 hist->GetXaxis()->SetLabelSize(.05);
 hist->GetYaxis()->SetLabelSize(.05);
 hist->SetStats(kFALSE);
 return hist;
}

TH2D* readHist2D(TString nameHist,TString nameFile, int rebin)
{
 TFile* file = new TFile(nameFile);

 TH2D* hist = (TH2D*)file->Get(nameHist);
 hist->GetSumw2();
 // hist->SetLineWidth(2);
 // if(rebin>0) hist->RebinX(rebin); hist->RebinY(rebin);
 hist->GetXaxis()->SetTitleSize(.055);
 hist->GetYaxis()->SetTitleSize(.055);
 hist->GetXaxis()->SetLabelSize(.05);
 hist->GetYaxis()->SetLabelSize(.05);
 hist->SetStats(kFALSE);
 return hist;
}

TH3D* readHist3D(TString nameHist,TString nameFile, int rebin)
{
 TFile* file = new TFile(nameFile);

 TH3D* hist = (TH3D*)file->Get(nameHist);
 hist->GetSumw2();
 // hist->SetLineWidth(2);
 // if(rebin>0) hist->RebinX(rebin); hist->RebinY(rebin);
 hist->GetXaxis()->SetTitleSize(.055);
 hist->GetYaxis()->SetTitleSize(.055);
 hist->GetXaxis()->SetLabelSize(.05);
 hist->GetYaxis()->SetLabelSize(.05);
 hist->SetStats(kFALSE);
 return hist;
}

TCanvas* getaCanvas(TString name)
{

  TCanvas* aCanvas = new TCanvas(name,"",292,55,500,700);//,"",181,237,1575,492);

  aCanvas->SetFillColor(0);
  aCanvas->SetBottomMargin(0.125);
  aCanvas->SetLeftMargin(0.125);
  aCanvas->SetFrameFillColor(0);
  aCanvas->SetFrameBorderMode(0);
  aCanvas->SetFrameLineWidth(2);
  return aCanvas;
}

TLegend *legend() {

 TLegend *leg2 = new TLegend(0.52,0.67,0.92,0.90);
 leg2->SetFillStyle(0);
 leg2->SetBorderSize(0);
 leg2->SetTextSize(0.05);
 leg2->SetTextFont(42);

 return leg2;

}

TPad *getaPad_up(TString name){

  TPad *pad1 = new TPad("name", "The pad with the function",0.05,0.4,0.95,0.1);
  //   pad1->Draw();

   //   pad1->Range(-112.6742,-73.17708,1143.438,551.3021);
   pad1->SetFillColor(0);
   pad1->SetBorderMode(0);
   pad1->SetBorderSize(2);
   pad1->SetGridx();
   pad1->SetGridy();
   pad1->SetLeftMargin(0.1271439);
   pad1->SetRightMargin(0.07307979);
   pad1->SetTopMargin(0.08215179);
   pad1->SetBottomMargin(0.117181);
   pad1->SetFrameBorderMode(0);
   pad1->SetFrameBorderMode(0);
  
   return pad1;

}

TPad *getaPad_dn(TString name){

  
 TPad *pad2 = new TPad("pad2", "The pad with the histogram",0.05,0.95,0.95,0.5);
 //   pad2->Draw();

   //   pad2->Range(-153.3652,-2.142584,1185.427,5.367464);
   pad2->SetFillColor(0);
   pad2->SetBorderMode(0);
   pad2->SetBorderSize(2);
   pad2->SetGridx();
   pad2->SetGridy();
   pad2->SetLeftMargin(0.1215511);
   pad2->SetRightMargin(0.07867263);
   pad2->SetTopMargin(0.04892967);
   pad2->SetBottomMargin(0.1521407);
   pad2->SetFrameBorderMode(0);
   pad2->SetFrameBorderMode(0);

   return pad2;

}

TH1D *hTemp(TString ifile, TString dirname, TString hname) {

  TH1D *h=readHist(dirname+hname,ifile,2);

  TH1D *h_temp=h->Clone();
  h_temp->Reset();

  Double_t con, err;

  Int_t nB=h->GetNbinsX();
  for (Int_t i=0; i<nB; i++) {

    con=h->IntegralAndError(i,nB+1,err);

    h_temp->SetBinContent(i,con);
    h_temp->SetBinError(i,err);
  }

  return h_temp;

}

void plot_singleEvents(TString sample="Minbias", TString mode="rate", Bool_t debug=true) {

  TString ifile, data;

  if (sample=="Minbias") {
    //    ifile="total_evt12609423_run207515.root";
    ifile="total_evt2492259_run205666.root";
    //    ifile="total_evt13345698_run207905.root";
    //    ifile="total_evt583532_run203909.root";
    //total_evt3688542_run205833.root";
    //  ifile="total_evt3349097_run205718.root";
    data=sample+" Run2012D";
  } else if (sample=="Singlemu") {
    ifile="l1Ntuple_IsoMu_RunC_1509.root";
    data="singleMu 2012C";
  }

  TString dirname="l1ExtraJets/";
  
  TH1D *ptJet_e, *ptJet_c, *ptJet_l;
  TH1D *ptJetHF_e, *ptJetHF_c, *ptJetHF_l;
  
 // Compute cummulative histograms

  if (mode=="rate") {
 //  // Central + Tau jets
    ptJet_e=hTemp(ifile,dirname,"ptJet_early");
    ptJet_c=hTemp(ifile,dirname,"ptJet_central");
    ptJet_l=hTemp(ifile,dirname,"ptJet_late");
    // HFjets
    ptJetHF_e=hTemp(ifile,dirname,"ptJetHF_early");
    ptJetHF_c=hTemp(ifile,dirname,"ptJetHF_central");
    ptJetHF_l=hTemp(ifile,dirname,"ptJetHF_late");


  } else if (mode=="norate") {
    ptJet_e=readHist(dirname+"ptJet_early",ifile,4);
    ptJet_c=readHist(dirname+"ptJet_central",ifile,4);
    ptJet_l=readHist(dirname+"ptJet_late",ifile,4);
    
    ptJetHF_e=readHist(dirname+"ptJetHF_early",ifile,4);
    ptJetHF_c=readHist(dirname+"ptJetHF_central",ifile,4);
    ptJetHF_l=readHist(dirname+"ptJetHF_late",ifile,4);
  }

  TH1D *etaJet_e=readHist(dirname+"etaJet_early",ifile,2);
  TH1D *etaJet_c=readHist(dirname+"etaJet_central",ifile,2);
  TH1D *etaJet_l=readHist(dirname+"etaJet_late",ifile,2);

  TH1D *etaJetHF_e=readHist(dirname+"etaJetHF_early",ifile,2);
  TH1D *etaJetHF_c=readHist(dirname+"etaJetHF_central",ifile,2);
  TH1D *etaJetHF_l=readHist(dirname+"etaJetHF_late",ifile,2);
  
  TH1D *phiJet_e=readHist(dirname+"phiJet_early",ifile,0);
  TH1D *phiJet_c=readHist(dirname+"phiJet_central",ifile,0);
  TH1D *phiJet_l=readHist(dirname+"phiJet_late",ifile,0);
 

  setTDRStyle();

  if (debug) {

    TH1D *h_matchDR_early=readHist(dirname+"h_matchDR_early",ifile,10);
    TH1D *h_matchDR_central=readHist(dirname+"h_matchDR_central",ifile,10);
    TH1D *h_matchDR_early_H=readHist(dirname+"h_matchDR_early_H",ifile,10);
    TH1D *h_matchDR_central_H=readHist(dirname+"h_matchDR_central_H",ifile,10);

    TH1D *h_matchDR_hf_early=readHist(dirname+"h_matchDR_hf_early",ifile,10);
    TH1D *h_matchDR_hf_central=readHist(dirname+"h_matchDR_hf_central",ifile,10);
    TH1D *h_matchDR_hf_early_H=readHist(dirname+"h_matchDR_hf_early_H",ifile,10);
    TH1D *h_matchDR_hf_central_H=readHist(dirname+"h_matchDR_hf_central_H",ifile,10);

    TH1D *etaJet_e_H=readHist(dirname+"etaJet_early_H",ifile,2);
    TH1D *etaJet_c_H=readHist(dirname+"etaJet_central_H",ifile,2);
    TH1D *phiJet_e_H=readHist(dirname+"phiJet_early_H",ifile,2);
    TH1D *phiJet_c_H=readHist(dirname+"phiJet_central_H",ifile,2);
    
    TH1D *etaJetHF_e_H=readHist(dirname+"etaJetHF_early_H",ifile,2);
    TH1D *etaJetHF_c_H=readHist(dirname+"etaJetHF_central_H",ifile,2);
    TH1D *phiJetHF_e_H=readHist(dirname+"phiJetHF_early_H",ifile,2);
    TH1D *phiJetHF_c_H=readHist(dirname+"phiJetHF_central_H",ifile,2);

    // HCAL Trigger Primitives info
    TH2D *ietaJet_iphiJet_central_H=readHist2D(dirname+"ietaJet_iphiJet_central_H",ifile,0);
    TH2D *ietaJet_iphiJet_hf_central_H=readHist2D(dirname+"ietaJet_iphiJet_hf_central_H",ifile,0);
    TH2D *ietaJet_iphiJet_central_Hs=readHist2D(dirname+"ietaJet_iphiJet_central_Hs",ifile,0);
    TH2D *ietaJet_iphiJet_hf_central_Hs=readHist2D(dirname+"ietaJet_iphiJet_hf_central_Hs",ifile,0);
    
 // For L1Jet ET>240 (BX=0)
    TH1D *ietaTP_H=readHist(dirname+"ietaTP_H",ifile,0);
    TH1D *iphiTP_H=readHist(dirname+"iphiTP_H",ifile,0);
    TH2D *ietaTP_iphiTP_H=readHist2D(dirname+"ietaTP_iphiTP_H",ifile,0);
    TH2D *ietaJet_iphiJet_H=readHist2D(dirname+"ietaJet_iphiJet_H",ifile,0);
    TH2D *ietaTP_iphiTP_H_2=readHist2D(dirname+"ietaTP_iphiTP_H_2",ifile,0);
    TH1D *etTP_H=readHist(dirname+"etTP_H",ifile,0);
    TH1D *compEtTP_H=readHist(dirname+"compEtTP_H",ifile,0);  
    //HF
    TH1D *ietaTP_hf_H=readHist(dirname+"ietaTP_hf_H",ifile,0);
    TH1D *iphiTP_hf_H=readHist(dirname+"iphiTP_hf_H",ifile,0);
    TH2D *ietaTP_iphiTP_hf_H=readHist2D(dirname+"ietaTP_iphiTP_hf_H",ifile,0);
    TH2D *ietaJet_iphiJet_hf_H=readHist2D(dirname+"ietaJet_iphiJet_hf_H",ifile,0);
    TH2D *ietaTP_iphiTP_hf_H_2=readHist2D(dirname+"ietaTP_iphiTP_hf_H_2",ifile,0);

    // For L1Jet ET>240 (BX=-1)
    TH1D *ietaTP_Hs=readHist(dirname+"ietaTP_Hs",ifile,0);
    TH1D *iphiTP_Hs=readHist(dirname+"iphiTP_Hs",ifile,0);
    TH2D *ietaTP_iphiTP_Hs=readHist2D(dirname+"ietaTP_iphiTP_Hs",ifile,0);
    TH2D *ietaJet_iphiJet_Hs=readHist2D(dirname+"ietaJet_iphiJet_Hs",ifile,0);
    TH2D *ietaTP_iphiTP_Hs_2=readHist2D(dirname+"ietaTP_iphiTP_Hs_2",ifile,0);
    TH1D *etTP_Hs=readHist(dirname+"etTP_Hs",ifile,0);
    TH1D *compEtTP_Hs=readHist(dirname+"compEtTP_Hs",ifile,0);  
    //HF
    TH1D *ietaTP_hf_Hs=readHist(dirname+"ietaTP_hf_Hs",ifile,0);
    TH1D *iphiTP_hf_Hs=readHist(dirname+"iphiTP_hf_Hs",ifile,0);
    TH2D *ietaTP_iphiTP_hf_Hs=readHist2D(dirname+"ietaTP_iphiTP_hf_Hs",ifile,0);
    TH2D *ietaJet_iphiJet_hf_Hs=readHist2D(dirname+"ietaJet_iphiJet_hf_Hs",ifile,0);
    TH2D *ietaTP_iphiTP_hf_Hs_2=readHist2D(dirname+"ietaTP_iphiTP_hf_Hs_2",ifile,0);

    // For all Events
    TH1D *ietaTP=readHist(dirname+"ietaTP",ifile,0);
    TH1D *iphiTP=readHist(dirname+"iphiTP",ifile,0);
    TH2D *ietaTP_iphiTP=readHist2D(dirname+"ietaTP_iphiTP",ifile,0);
    TH2D *ietaTP_iphiTP_2=readHist2D(dirname+"ietaTP_iphiTP_2",ifile,0);

    gStyle->SetPalette(1);

    TCanvas *tp1=getaCanvas("tp1");
    tp1->Divide(2,1);

    tp1->cd(1); 
    TLegend *leg = legend();
    leg->SetHeader("HCAL TPs");
    leg->AddEntry(ietaTP,"All events","LF");
    leg->AddEntry(ietaTP_H,"Events with BX=0 L1jet of ET>240","LF");
    // leg->AddEntry(ietaTP_hf_H,"Events with BX=-1 fwd L1jet ET>240","LF");

    ietaTP->Scale(1./ietaTP->Integral());
    ietaTP_H->Scale(1./ietaTP_H->Integral());
    ietaTP_hf_H->Scale(1./ietaTP_hf_H->Integral());

    ietaTP->Draw("EHIST"); ietaTP->SetLineColor(1);
    ietaTP_H->Draw("EHISTSAME");  ietaTP_H->SetLineWidth(2);
    ietaTP->SetTitle(";hcalTP i#eta;a.u.");
    //  ietaTP_hf_H->Draw("EHISTSAME");  ietaTP_hf_H->SetLineWidth(2);
    //  ietaTP_hf_H->SetLineColor(3);

    leg->Draw("SAME");
    
    tp1->cd(2); 
    TLegend *leg = legend();
    leg->SetHeader("HCAL TPs");
    leg->AddEntry(iphiTP,"All events","LF");
    leg->AddEntry(iphiTP_H,"Events with BX=0 L1jet ET>240","LF");
    //  leg->AddEntry(iphiTP_hf_H,"Events with BX=-1 fwd L1jet ET>240","LF");

    iphiTP->Scale(1./iphiTP->Integral());
    iphiTP_H->Scale(1./iphiTP_H->Integral());
    iphiTP_hf_H->Scale(1./iphiTP_hf_H->Integral());

    //  iphiTP->Draw("EHIST"); iphiTP->SetLineColor(1);
    iphiTP->Draw("EHIST");  iphiTP_H->SetLineWidth(2);
    iphiTP_H->Draw("EHISTSAME"); iphiTP->SetLineColor(1);
    iphiTP_H->SetTitle(";hcalTP i#phi;a.u.");
    //  iphiTP_hf_H->Draw("EHISTSAME");  iphiTP_hf_H->SetLineWidth(2);
    //  iphiTP_hf_H->SetLineColor(3);

    leg->Draw("SAME");

    //return;


    TCanvas *tp3=getaCanvas("tp3");
    tp3->Divide(2,1);

    tp3->cd(1); 
    TLegend *leg = legend();
    leg->SetHeader("hcalTPs for events with BX=0 Jets of saturated E_{T}");

    ietaTP_iphiTP->Draw("ZCOL");
    ietaTP_iphiTP->SetTitle(";hcalTP i#eta;hcalTP i#phi");

    leg->Draw("SAME");

    tp3->cd(2);
    TLegend *leg = legend();
    leg->SetHeader("hcalTPs for events with BX=0 Jets of saturated E_{T}");

    ietaTP_iphiTP_2->Draw("ZCOL");
    ietaTP_iphiTP_2->SetTitle(";hcalTP i#eta;hcalTP i#phi");
    leg->Draw("SAME");


    TCanvas *tp33=getaCanvas("tp33");
    tp33->Divide(3,1);


    tp33->cd(1); 
    TLegend *leg = legend();
    leg->SetHeader("hcalTPs for events with a BX=-1 Jet of saturated E_{T}");

    ietaTP_iphiTP_H->Draw("ZCOL");
    ietaTP_iphiTP_H->SetTitle(";hcalTP i#eta;hcalTP i#phi");

    leg->Draw("SAME");

    tp33->cd(2);
    TLegend *leg = legend();
    leg->SetHeader("HB/HE: BX=0 Jets with saturated E_{T}");

    ietaJet_iphiJet_central_H->Draw("ZCOL");
    ietaJet_iphiJet_central_H->SetTitle(";BX=0 L1Jet i#eta;BX=0 L1Jet i#phi");

    leg->Draw("SAME");

    tp33->cd(3);
    TLegend *leg = legend();
    leg->SetHeader("HB/HE: BX=0 Jets with saturated E_{T}");

    ietaJet_iphiJet_H->Draw("ZCOL");
    ietaJet_iphiJet_H->SetTitle(";BX=-1 L1Jet i#eta;BX=-1 L1Jet i#phi");

    leg->Draw("SAME");

    return;
    
    TCanvas *tp334=getaCanvas("tp334");
    tp334->Divide(3,1);

    tp334->cd(1); 
    TLegend *leg = legend();
    leg->SetHeader("hcalTPs for events with BX=0 HFJets of saturated E_{T}");

    ietaTP_iphiTP_hf_H->Draw("ZCOL");
    ietaTP_iphiTP_hf_H->SetTitle(";hcalTP i#eta;hcalTP i#phi");

    // leg->Draw("SAME");

    tp334->cd(2);
    TLegend *leg = legend();
    leg->SetHeader("HF: BX=0 Jets with saturated E_{T}");

    ietaJet_iphiJet_hf_central_H->Draw("ZCOL");
    ietaJet_iphiJet_hf_central_H->SetTitle(";L1 fwdJet i#eta;L1 fwdJet i#phi");

    leg->Draw("SAME");

    tp334->cd(3);
    TLegend *leg = legend();
    leg->SetHeader("HF: BX=-1 Jets with saturated E_{T}");

    ietaJet_iphiJet_hf_H->Draw("ZCOL");
    ietaJet_iphiJet_hf_H->SetTitle(";L1 fwdJet i#eta;L1 fwdJet i#phi");

    leg->Draw("SAME");
    
    // return;

    TCanvas *tp2=getaCanvas("tp2");
    tp2->Divide(2,1);

    tp2->cd(1);
    etTP_H->Draw("EHIST");
    etTP_H->SetLineWidth(2); etTP_H->SetTitle(";hcalTP E_{T};a.u");
    etTP_H->GetXaxis()->SetRangeUser(0.,20.);

    tp2->cd(2);
    compEtTP_H->Draw("EHIST");
    compEtTP_H->SetLineWidth(2); compEtTP_H->SetTitle(";hcalTP comp E_{T};a.u");
    compEtTP_H->GetXaxis()->SetRangeUser(0.,20.);

    // HBHE eta/phi
    // ----------> Early triggers
    TCanvas *d1=getaCanvas("d1");
    d1->Divide(2,1);

    d1->cd(1); 
    etaJet_e_H->Draw("EHIST"); etaJet_e_H->SetLineWidth(2);
    etaJet_e_H->SetTitle(";(BX=-1) L1 Jet #eta;a.u.");

    d1->cd(2);
    phiJet_e_H->Draw("EHIST"); phiJet_e_H->SetLineWidth(2);
    phiJet_e_H->SetTitle(";(BX=-1) L1 Jet #phi;a.u.");

     // ----------> Central triggers
    TCanvas *d2=getaCanvas("d2");
    d2->Divide(2,1);

    d2->cd(1); 
    etaJet_c_H->Draw("EHIST"); etaJet_c_H->SetLineWidth(2);
    etaJet_c_H->SetTitle(";(BX=0) L1 Jet #eta;a.u.");

    d2->cd(2);
    phiJet_c_H->Draw("EHIST"); phiJet_c_H->SetLineWidth(2);
    phiJet_c_H->SetTitle(";(BX=0) L1 Jet #phi;a.u.");

    // HF
     // ----------> Early triggers
    TCanvas *d3=getaCanvas("d3");
    d3->Divide(2,1);

    d3->cd(1); 
    etaJetHF_e_H->Draw("EHIST"); etaJetHF_e_H->SetLineWidth(2);
    etaJetHF_e_H->SetTitle(";(BX=-1) L1 fwdJet #eta;a.u.");

    d3->cd(2);
    phiJetHF_e_H->Draw("EHIST"); phiJetHF_e_H->SetLineWidth(2);
    phiJetHF_e_H->SetTitle(";(BX=-1) L1 fwdJet #phi;a.u.");

 // ----------> Central triggers
    TCanvas *d4=getaCanvas("d4");
    d4->Divide(2,1);

    d4->cd(1); 
    etaJetHF_c_H->Draw("EHIST"); etaJetHF_c_H->SetLineWidth(2);
    etaJetHF_c_H->SetTitle(";(BX=0) L1 fwdJet #eta;a.u.");

    d4->cd(2);
    phiJetHF_c_H->Draw("EHIST"); phiJetHF_c_H->SetLineWidth(2);
    phiJetHF_c_H->SetTitle(";(BX=0) L1 fwdJet #phi;a.u.");


  }

  //  return;

  // Central+Tau Jets
  TCanvas *c1=getaCanvas("c1");
  c1->Divide(1,2);

  TPad *pad1 = new TPad("pad1", "The pad with the function",0.05,0.4,0.95,0.95);
  pad1->Draw();
  
   //   pad1->Range(-112.6742,-73.17708,1143.438,551.3021);
  pad1->SetFillColor(0);
  pad1->SetBorderMode(0);
  pad1->SetBorderSize(2);
 //  pad1->SetGridx();
//   pad1->SetGridy();
  pad1->SetLeftMargin(0.1271439);
  pad1->SetRightMargin(0.07307979);
  pad1->SetTopMargin(0.08215179);
  pad1->SetBottomMargin(0.117181);
  pad1->SetFrameBorderMode(0);
  pad1->SetFrameBorderMode(0);
  TPad *pad2 = new TPad("pad2", "The pad with the histogram",0.05,0.05,0.95,0.4);
  pad2->Draw();
  
  //   pad2->Range(-153.3652,-2.142584,1185.427,5.367464);
  pad2->SetFillColor(0);
  pad2->SetBorderMode(0);
  pad2->SetBorderSize(2);
  pad2->SetGridx();
  pad2->SetGridy();
  pad2->SetLeftMargin(0.1215511);
  pad2->SetRightMargin(0.07867263);
  pad2->SetTopMargin(0.04892967);
  pad2->SetBottomMargin(0.3521407);
  pad2->SetFrameBorderMode(0);
  pad2->SetFrameBorderMode(0);
  
  pad1->cd();
  gPad->SetLogy();
  //  gPad->SetGridx(); gPad->SetGridy();
 
  TLegend *leg = legend();
  leg->SetHeader(data);
  //  leg->AddEntry(ptJet_e,"BX=-1","LF");
  leg->AddEntry(ptJet_c,"BX=0","LF");
  leg->AddEntry(ptJet_e,"BX=-1","LF");
  //  leg->AddEntry(ptJet_l,"BX=+1","LF");

  ptJet_c->Draw("EHIST"); 
  ptJet_c->SetLineColor(1); ptJet_c->SetLineWidth(2);
  ptJet_c->GetXaxis()->SetRangeUser(0.,260.);
  ptJet_c->GetYaxis()->SetRangeUser(1.01,100000000.);

  if (mode=="rate") {
    ptJet_c->SetTitle(";L1Jet E_{T} Threshold (GeV);a.u.");
  } else {
    ptJet_c->SetTitle(";L1Jet E_{T} (GeV);a.u.");
  }
  ptJet_e->Draw("EHISTSAMES"); 
  ptJet_e->SetLineColor(1); ptJet_e->SetLineColor(2);
  ptJet_e->SetLineStyle(1);ptJet_e->SetLineWidth(2);

//   ptJet_l->Draw("EHISTSAMES"); 
//   ptJet_l->SetLineColor(1); ptJet_l->SetLineColor(4);
//   ptJet_l->SetLineStyle(3);

  leg->Draw("SAME");

 //  TPad *p2=getaPad_dn("p2"); p2->Draw();
//   p2->cd();
  //c1->cd(2);
  pad2->cd();
  gPad->SetGridx(); gPad->SetGridy();

  TH1D *r1_e = ptJet_e->Clone();
  r1_e->Divide(ptJet_c);
  TH1D *r1_l = ptJet_l->Clone();
  r1_l->Divide(ptJet_c);

  TLegend *leg = legend();
  leg->SetHeader(data);
  leg->AddEntry(r1_e,"BX=-1","LF");
  //leg->AddEntry(r1_l,"BX=+1","LF");

  r1_e->Draw("EHIST");r1_e->SetLineColor(4);
  //  r1_l->Draw("EHISTSAME"); r1_l->SetLineColor(4);

  r1_e->GetYaxis()->SetRangeUser(-0.01,1.0);

  r1_e->GetXaxis()->SetRangeUser(0.,260.);
  r1_e->SetLineWidth(2); r1_l->SetLineWidth(2);

  if (mode=="rate") {
    r1_e->SetTitle(";L1Jet E_{T} Threshold (GeV);Ratio");
  } else {
    r1_e->SetTitle(";L1Jet E_{T} (GeV);Ratio");
  }
  //  r1_e->SetTitle(";L1Jet E_{T} (GeV);ratio");

  r1_e->GetXaxis()->SetLabelSize(.1);
  r1_e->GetYaxis()->SetLabelSize(.1);
  r1_e->GetXaxis()->SetTitleSize(.09);
  r1_e->GetYaxis()->SetTitleSize(.095);
  r1_e->GetXaxis()->SetTitleOffset(1.0);
  r1_e->GetYaxis()->SetTitleOffset(0.6);

  //leg->Draw("SAME");

  //return;

  // Eta Jets
  TCanvas *c2=getaCanvas("c2");
  c2->Divide(1,2);
  // TPad *p1=getaPad_up("p1"); p1->Draw();
//   p1->cd();
  TPad *pad3 = new TPad("pad3", "The pad with the function",0.05,0.4,0.95,0.95);
  pad3->Draw();
  
   //   pad3->Range(-112.6742,-73.17708,1143.438,551.3021);
  pad3->SetFillColor(0);
  pad3->SetBorderMode(0);
  pad3->SetBorderSize(2);
//   pad3->SetGridx();
//   pad3->SetGridy();
  pad3->SetLeftMargin(0.1271439);
  pad3->SetRightMargin(0.07307979);
  pad3->SetTopMargin(0.08215179);
  pad3->SetBottomMargin(0.117181);
  pad3->SetFrameBorderMode(0);
  pad3->SetFrameBorderMode(0);
  TPad *pad4 = new TPad("pad4", "The pad with the histogram",0.05,0.1,0.95,0.4);
  pad4->Draw();
  
  //   pad4->Range(-153.3652,-2.142584,1185.427,5.367464);
  pad4->SetFillColor(0);
  pad4->SetBorderMode(0);
  pad4->SetBorderSize(2);
  pad4->SetGridx();
  pad4->SetGridy();
  pad4->SetLeftMargin(0.1215511);
  pad4->SetRightMargin(0.07867263);
  pad4->SetTopMargin(0.04892967);
  pad4->SetBottomMargin(0.3521407);
  pad4->SetFrameBorderMode(0);
  pad4->SetFrameBorderMode(0);



  pad3->cd();
  // c2->cd(1);
  //  gPad->SetGridx(); gPad->SetGridy();

  TLegend *leg = legend();
  leg->SetHeader(data);
  //  leg->AddEntry(etaJet_e,"BX=-1","LF");
  leg->AddEntry(etaJet_c,"BX=0","LF");
  // leg->AddEntry(etaJet_l,"BX=+1","LF");
  leg->AddEntry(etaJet_e,"BX=-1","LF");

  etaJet_c->Draw("EHIST"); 
  etaJet_c->SetLineColor(1); etaJet_c->SetLineWidth(2);
  etaJet_c->SetTitle(";L1Jet #eta;a.u.");
  //  etaJet_c->GetYaxis()->SetRangeUser(1.01,100000000.);
  etaJet_e->Draw("EHISTSAMES"); 
  etaJet_e->SetLineColor(1); etaJet_e->SetLineColor(2);  etaJet_e->SetLineStyle(1);etaJet_e->SetLineWidth(2);
 //  etaJet_l->Draw("EHISTSAMES"); 
//   etaJet_l->SetLineColor(1); etaJet_l->SetLineColor(4);   etaJet_l->SetLineStyle(3);

  leg->Draw("SAME");
 
  pad4->cd();
  //  c2->cd(2);
  gPad->SetGridx(); gPad->SetGridy();

  TH1D *r2_e = etaJet_e->Clone();
  r2_e->Divide(etaJet_c);
  TH1D *r2_l = etaJet_l->Clone();
  r2_l->Divide(etaJet_c);

  TLegend *leg = legend();
  leg->SetHeader(data);
  leg->AddEntry(r2_e,"BX=-1","LF");
  // leg->AddEntry(etaJet_c,"BX=0","LF");
  // leg->AddEntry(r2_l,"BX=+1","LF");

  r2_e->Draw("HIST");r2_e->SetLineColor(4);r2_e->SetLineWidth(2);
  // r2_l->Draw("HISTSAME"); r2_l->SetLineColor(4);

  r2_e->GetYaxis()->SetRangeUser(-0.01,1.0);

  //r2_e->GetXaxis()->SetRangeUser(0.,300.);
  r2_e->SetLineWidth(2); r2_l->SetLineWidth(2);
  r2_e->SetTitle(";L1Jet #eta;Ratio ");

  r2_e->GetXaxis()->SetTitleSize(.09);
  r2_e->GetYaxis()->SetTitleSize(.095);
  r2_e->GetXaxis()->SetTitleOffset(1.0);
  r2_e->GetYaxis()->SetTitleOffset(0.6);
  r2_e->GetXaxis()->SetLabelSize(.1);
  r2_e->GetYaxis()->SetLabelSize(.1);

 
  // leg->Draw("SAME");


  // HF jets


 TCanvas *c3=getaCanvas("c3");
  c3->Divide(1,2);

 TPad *pad5 = new TPad("pad5", "The pad with the function",0.05,0.4,0.95,0.95);
  pad5->Draw();
  
   //   pad5->Range(-112.6742,-73.17708,1143.438,551.3021);
  pad5->SetFillColor(0);
  pad5->SetBorderMode(0);
  pad5->SetBorderSize(2);
  //  pad5->SetGridx();
  // pad5->SetGridy();
  pad5->SetLeftMargin(0.1271439);
  pad5->SetRightMargin(0.07307979);
  pad5->SetTopMargin(0.08215179);
  pad5->SetBottomMargin(0.117181);
  pad5->SetFrameBorderMode(0);
  pad5->SetFrameBorderMode(0);
  TPad *pad6 = new TPad("pad6", "The pad with the histogram",0.05,0.1,0.95,0.4);
  pad6->Draw();
  
  //   pad6->Range(-153.3652,-2.142584,1185.427,5.367464);
  pad6->SetFillColor(0);
  pad6->SetBorderMode(0);
  pad6->SetBorderSize(2);
  pad6->SetGridx();
  pad6->SetGridy();
  pad6->SetLeftMargin(0.1215511);
  pad6->SetRightMargin(0.07867263);
  pad6->SetTopMargin(0.04892967);
  pad6->SetBottomMargin(0.3521407);
  pad6->SetFrameBorderMode(0);
  pad6->SetFrameBorderMode(0);

  pad5->cd();
  //  c3->cd(1);
  gPad->SetLogy();
  //  gPad->SetGridx(); gPad->SetGridy();

  TLegend *leg = legend();
  leg->SetHeader(data);
  //  leg->AddEntry(ptJetHF_e,"BX=-1","LF");
  leg->AddEntry(ptJetHF_c,"BX=0","LF");
  //  leg->AddEntry(ptJetHF_l,"BX=+1","LF");
  leg->AddEntry(ptJetHF_e,"BX=-1","LF");

  ptJetHF_c->Draw("EHIST"); 
  ptJetHF_c->SetLineColor(1); ptJetHF_c->SetLineWidth(2);
  ptJetHF_c->GetXaxis()->SetRangeUser(0.,260.);
  ptJetHF_c->GetYaxis()->SetRangeUser(1.01,100000000.);
  //  ptJetHF_c->SetTitle(";L1 FwdJet E_{T} (GeV);a.u.");
 if (mode=="rate") {
    ptJetHF_c->SetTitle(";L1 FwdJet E_{T} Threshold (GeV);a.u.");
  } else {
    ptJetHF_c->SetTitle(";L1 FwdJet E_{T} (GeV);a.u.");
  }

  ptJetHF_e->Draw("EHISTSAMES"); 
  ptJetHF_e->SetLineColor(1); ptJetHF_e->SetLineColor(2);
  ptJetHF_e->SetLineStyle(1); ptJetHF_e->SetLineWidth(2);
 //  ptJetHF_l->Draw("EHISTSAMES"); 
//   ptJetHF_l->SetLineColor(1); ptJetHF_l->SetLineColor(4);
//   ptJetHF_l->SetLineStyle(3);

  leg->Draw("SAME");

 //  TPad *p2=getaPad_dn("p2"); p2->Draw();
//   p2->cd();
  // c3->cd(2);
  pad6->cd();
  gPad->SetGridx(); gPad->SetGridy();

  TH1D *r1HF_e = ptJetHF_e->Clone();
  r1HF_e->Divide(ptJetHF_c);
  TH1D *r1HF_l = ptJetHF_l->Clone();
  r1HF_l->Divide(ptJetHF_c);

  TLegend *leg = legend();
  leg->AddEntry(r1HF_e,"BX=-1","LF");
  // leg->AddEntry(ptJetHF_c,"BX=0","LF");
  //  leg->AddEntry(r1HF_l,"BX=+1","LF");

  r1HF_e->Draw("EHIST");r1HF_e->SetLineColor(4);
  // r1HF_l->Draw("EHISTSAME"); r1HF_l->SetLineColor(4);

  r1HF_e->GetYaxis()->SetRangeUser(-0.01,1.0);

  r1HF_e->GetXaxis()->SetRangeUser(0.,260.);
  r1HF_e->SetLineWidth(2); r1HF_l->SetLineWidth(2);
  //  r1HF_e->SetTitle(";L1 FwdJet E_{T} (GeV);ratio");
  if (mode=="rate") {
    r1HF_e->SetTitle(";L1 FwdJet E_{T} Threshold (GeV);Ratio");
  } else {
    r1HF_e->SetTitle(";L1 FwdJet E_{T} (GeV);Ratio");
  }

  r1HF_e->GetXaxis()->SetTitleSize(.09);
  r1HF_e->GetYaxis()->SetTitleSize(.095);
  r1HF_e->GetXaxis()->SetTitleOffset(1.0);
  r1HF_e->GetYaxis()->SetTitleOffset(0.6);
  r1HF_e->GetXaxis()->SetLabelSize(.1);
  r1HF_e->GetYaxis()->SetLabelSize(.1);

  //  leg->Draw("SAME");

}
