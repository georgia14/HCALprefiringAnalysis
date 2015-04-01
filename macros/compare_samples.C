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
TCanvas* getanotherCanvas(TString name)
{

  TCanvas* aCanvas = new TCanvas(name,"",292,55,500,465);//,"",181,237,1575,492);

  aCanvas->SetFillColor(0);
  aCanvas->SetBottomMargin(0.125);
  aCanvas->SetLeftMargin(0.125);
  aCanvas->SetFrameFillColor(0);
  aCanvas->SetFrameBorderMode(0);
  aCanvas->SetFrameLineWidth(2);
  return aCanvas;
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

void compare_samples(TString sample="Singlemu", TString mode="rate") {

  TString ifile1, ifile2, data;

  if (sample=="Minbias") {
    ifile1="l1Ntuple_ZeroBias_RunD_oldPMTs_HF-.root";
    ifile2="l1Ntuple_ZeroBias_RunD_upgradedPMTs_HF-.root";
    data="MinBias Run2012D";
  } else if (sample=="Singlemu") {
    //    ifile1="l1Ntuple_IsoMu_RunC_1509.root";
    ifile1="l1Ntuple_SingleMu_RunD_oldPMTs.root";
    ifile2="l1Ntuple_SingleMu_RunD_upgradedPMTs.root";
    data="singleMu 2012D";
  }

  TString dirname="l1ExtraJets/";
  
  TH1D *ptJet_e, *ptJet_c, *ptJet_l;
  TH1D *ptJetHF_e, *ptJetHF_c, *ptJetHF_l;
  
 // Compute cummulative histograms

  if (mode=="rate") {
 //  // Central + Tau jets
    ptJet_e=hTemp(ifile1,dirname,"ptJet_early");
    ptJet_c=hTemp(ifile1,dirname,"ptJet_central");
    ptJet_l=hTemp(ifile1,dirname,"ptJet_late");
    // HFjets
    ptJetHF_e=hTemp(ifile1,dirname,"ptJetHF_early");
    ptJetHF_c=hTemp(ifile1,dirname,"ptJetHF_central");
    ptJetHF_l=hTemp(ifile1,dirname,"ptJetHF_late");

    ptJetHF2_e=hTemp(ifile2,dirname,"ptJetHF_early");
    ptJetHF2_c=hTemp(ifile2,dirname,"ptJetHF_central");
    ptJetHF2_l=hTemp(ifile2,dirname,"ptJetHF_late");

  } else if (mode=="norate") {
    ptJet_e=readHist(dirname+"ptJet_early",ifile1,4);
    ptJet_c=readHist(dirname+"ptJet_central",ifile1,4);
    ptJet_l=readHist(dirname+"ptJet_late",ifile1,4);
    
    ptJetHF_e=readHist(dirname+"ptJetHF_early",ifile1,4);
    ptJetHF_c=readHist(dirname+"ptJetHF_central",ifile1,4);
    ptJetHF_l=readHist(dirname+"ptJetHF_late",ifile1,4);

    ptJetHF2_e=readHist(dirname+"ptJetHF_early",ifile2,4);
    ptJetHF2_c=readHist(dirname+"ptJetHF_central",ifile2,4);
    ptJetHF2_l=readHist(dirname+"ptJetHF_late",ifile2,4);

  }

  TH1D *etaJet_e=readHist(dirname+"etaJet_early",ifile1,2);
  TH1D *etaJet_c=readHist(dirname+"etaJet_central",ifile1,2);
  TH1D *etaJet_l=readHist(dirname+"etaJet_late",ifile1,2);

  TH1D *etaJetHF_e=readHist(dirname+"etaJetHF_early",ifile1,2);
  TH1D *etaJetHF_c=readHist(dirname+"etaJetHF_central",ifile1,2);
  TH1D *etaJetHF_l=readHist(dirname+"etaJetHF_late",ifile1,2);
  
  TH1D *phiJet_e=readHist(dirname+"phiJet_early",ifile1,0);
  TH1D *phiJet_c=readHist(dirname+"phiJet_central",ifile1,0);
  TH1D *phiJet_l=readHist(dirname+"phiJet_late",ifile1,0);
 

  setTDRStyle();

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
  leg->AddEntry(ptJetHF_c,"old PMTs region: BX=0 Jet","LF");
  leg->AddEntry(ptJetHF2_c,"upgraded PMTs region: BX=0 Jet","LF");
  leg->AddEntry(ptJetHF_e,"old PMTs region: BX=-1 Jet","LF");
  leg->AddEntry(ptJetHF2_e,"upgraded PMTs region: BX=-1 Jet","LF");

  ptJetHF_c->Rebin(2); ptJetHF2_c->Rebin(2);
  ptJetHF_e->Rebin(2); ptJetHF2_e->Rebin(2);

  ptJetHF_c->Draw("EHIST"); 
  ptJetHF2_c->Draw("EHISTSAME");
  ptJetHF2_c->SetLineColor(1); //ptJetHF2_c->SetLineWidth(2);  
  //  ptJetHF2_c->SetFillStyle(3002); ptJetHF2_c->SetFillColor(1);
  

  ptJetHF_c->SetLineColor(1); ptJetHF_c->SetLineWidth(2);
  ptJetHF_c->GetXaxis()->SetRangeUser(0.,260.);
  ptJetHF_c->GetYaxis()->SetRangeUser(1.01,100000000.);
  //  ptJetHF_c->SetTitle(";L1 FwdJet E_{T} (GeV);a.u.");
 if (mode=="rate") {
    ptJetHF_c->SetTitle(";L1 FwdJet E_{T} Threshold (GeV);a.u.");
  } else {
    ptJetHF_c->SetTitle(";L1 FwdJet E_{T} (GeV);a.u.");
  }

  ptJetHF_e->Draw("EHISTSAME"); 
  ptJetHF2_e->Draw("EHISTSAME");
  ptJetHF2_e->SetLineColor(2);// ptJetHF2_e->SetLineWidth(2); 
  //  ptJetHF2_e->SetFillStyle(3002);ptJetHF2_e->SetFillColor(2);

  ptJetHF_e->SetLineColor(2);
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

  TH1D *r1HF2_e = ptJetHF2_e->Clone();
  r1HF2_e->Divide(ptJetHF2_c);

  TLegend *leg = legend();
  leg->AddEntry(r1HF_e,"old PMTs region","LF");
  leg->AddEntry(r1HF2_e,"upgraded PMTs region","LF");
  // leg->AddEntry(ptJetHF_c,"BX=0","LF");
  //  leg->AddEntry(r1HF_l,"BX=+1","LF");

  r1HF_e->Draw("EHIST");r1HF_e->SetLineColor(4);
  // r1HF_l->Draw("EHISTSAME"); r1HF_l->SetLineColor(4);
  r1HF2_e->Draw("EHISTSAME");  //r1HF2_e->SetLineColor(kBlack-6);
  //r1HF2_e->SetFillColor(4); 
  r1HF2_e->SetLineColor(4); r1HF2_e->SetLineWidth(0.5);

  r1HF_e->GetYaxis()->SetRangeUser(-0.01,1.0);

  r1HF_e->GetXaxis()->SetRangeUser(0.,260.);
  r1HF_e->SetLineWidth(2); //r1HF2_e->SetLineWidth(1);
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

  leg->Draw("SAME");

  TCanvas *cf=getanotherCanvas("cf");
  gPad->SetGridx(); gPad->SetGridy();

   TLegend *leg = legend();
  leg->AddEntry(r1HF_e,"old PMTs region","LF");
  leg->AddEntry(r1HF2_e,"upgraded PMTs region","LF");

  r1HF_e->Draw("EHIST");r1HF_e->SetLineColor(4);
  //  r1HF2_e->SetLineColor(1);
  r1HF2_e->Draw("EHISTSAME");  //r1HF2_e->SetLineColor(kBlack-6);
  r1HF2_e->SetLineWidth(2); r1HF2_e->SetLineStyle(2);


  r1HF_e->GetXaxis()->SetRangeUser(0.,170.);
  r1HF_e->GetYaxis()->SetRangeUser(0.001,0.3);

  r1HF_e->GetXaxis()->SetTitleSize(.05);
  r1HF_e->GetYaxis()->SetTitleSize(.05);
  r1HF_e->GetXaxis()->SetTitleOffset(1.0);
  r1HF_e->GetYaxis()->SetTitleOffset(0.6);
  r1HF_e->GetXaxis()->SetLabelSize(.05);
  r1HF_e->GetYaxis()->SetLabelSize(.05);

leg->Draw("SAME");
}
