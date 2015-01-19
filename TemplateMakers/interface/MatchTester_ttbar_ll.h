#ifndef _MatchTester_ttbar_ll_h
#define _MatchTester_ttbar_ll_h

#include "ttH-13TeVMultiLeptons/TemplateMakers/interface/KinematicVariable.h"
#include "ttH-13TeVMultiLeptons/TemplateMakers/interface/BranchInfo.h"
//#include <typeinfo>

class MatchTester_ttbar_ll: public KinematicVariable<double> {

public:

  //Store branch values so they are accessible to other classes
  vector<BranchInfo<double>> myVars;

  BNleptonCollection **leptons;
  BNjetCollection **jets;

  TFile * weight_file;
  TH1 * ratio_top_2jet_CSV;
  TH1 * ratio_top_3jet_CSV;
  TH1 * ratio_top_4jet_CSV;
  TH1 * ratio_top_jet_charge;
  TH1 * ratio_antiTop_jet_charge;
  TH1 * ratio_top_mass_lep_b;
  TH1 * ratio_ttbar_MT_mass_ratio_B_b;

  MatchTester_ttbar_ll(BNleptonCollection **_leptons, BNjetCollection **_jets);
  ~MatchTester_ttbar_ll();
  void evaluate();
};

MatchTester_ttbar_ll::MatchTester_ttbar_ll(BNleptonCollection **_leptons, BNjetCollection **_jets):
  leptons(_leptons), jets(_jets) {

  //std::cout << "Initializing MatchTester_ttbar_ll" << std::endl;

  //this->resetVal = 0.0; //What does this do?
  this->resetVal = KinematicVariableConstants::FLOAT_INIT; //What does this do?

  branches["Match_ttbar_ll_B"] = BranchInfo<double>("Match_ttbar_ll_B");
  branches["Match_ttbar_ll_b"] = BranchInfo<double>("Match_ttbar_ll_b");
  branches["Match_ttbar_ll_Bb"] = BranchInfo<double>("Match_ttbar_ll_Bb");

  branches["Full_match_ttbar_ll_B"] = BranchInfo<double>("Full_match_ttbar_ll_B");
  branches["Full_match_ttbar_ll_b"] = BranchInfo<double>("Full_match_ttbar_ll_b");
  branches["Full_match_ttbar_ll_Bb"] = BranchInfo<double>("Full_match_ttbar_ll_Bb");

  branches["ttbar_ll_B_CSV"] = BranchInfo<double>("ttbar_ll_B_CSV");
  branches["ttbar_ll_b_CSV"] = BranchInfo<double>("ttbar_ll_b_CSV");
  branches["ttbar_ll_B_charge"] = BranchInfo<double>("ttbar_ll_B_charge");
  branches["ttbar_ll_b_charge"] = BranchInfo<double>("ttbar_ll_b_charge");
  branches["ttbar_ll_top_mass_lep_B"] = BranchInfo<double>("ttbar_ll_top_mass_lep_B");
  branches["ttbar_ll_antiTop_mass_lep_b"] = BranchInfo<double>("ttbar_ll_antiTop_mass_lep_b");
  branches["ttbar_ll_ttbar_MT_mass_ratio_B_b"] = BranchInfo<double>("ttbar_ll_ttbar_MT_mass_ratio_B_b");

  branches["Match_ttbar_ll_B"].branchVal = KinematicVariableConstants::FLOAT_INIT; 
  branches["Match_ttbar_ll_b"].branchVal = KinematicVariableConstants::FLOAT_INIT; 
  branches["Match_ttbar_ll_Bb"].branchVal = KinematicVariableConstants::FLOAT_INIT; 

  branches["Full_match_ttbar_ll_B"].branchVal = 0.0; 
  branches["Full_match_ttbar_ll_b"].branchVal = 0.0; 
  branches["Full_match_ttbar_ll_Bb"].branchVal = 0.0; 

  branches["ttbar_ll_B_CSV"].branchVal = KinematicVariableConstants::FLOAT_INIT; 
  branches["ttbar_ll_b_CSV"].branchVal = KinematicVariableConstants::FLOAT_INIT;
  branches["ttbar_ll_B_charge"].branchVal = KinematicVariableConstants::FLOAT_INIT;
  branches["ttbar_ll_b_charge"].branchVal = KinematicVariableConstants::FLOAT_INIT;
  branches["ttbar_ll_top_mass_lep_B"].branchVal = KinematicVariableConstants::FLOAT_INIT;
  branches["ttbar_ll_antiTop_mass_lep_b"].branchVal = KinematicVariableConstants::FLOAT_INIT;
  branches["ttbar_ll_ttbar_MT_mass_ratio_B_b"].branchVal = KinematicVariableConstants::FLOAT_INIT;

  //std::cout << "Getting weight file" << std::endl;
  string directory = (string(getenv("CMSSW_BASE"))+"/src/ttH-13TeVMultiLeptons/TemplateMakers/data/NOVa/matchbox/").c_str();
//   TString weight_file_name = Form("%smatch_ttbarW_3l.root", directory.c_str());
  TString weight_file_name = Form("%smatch_ttbar_ll_3l.root", directory.c_str());
  weight_file = TFile::Open(weight_file_name);
  //std::cout << weight_file_name << std::endl;

  //std::cout << "Cloning histograms" << std::endl;  
  ratio_top_2jet_CSV = (TH1*)weight_file->Get("ratio_top_2jet_CSV")->Clone();
  ratio_top_3jet_CSV = (TH1*)weight_file->Get("ratio_top_3jet_CSV")->Clone();
  ratio_top_4jet_CSV = (TH1*)weight_file->Get("ratio_top_4jet_CSV")->Clone();
  ratio_top_jet_charge = (TH1*)weight_file->Get("ratio_top_jet_charge")->Clone();
  ratio_antiTop_jet_charge = (TH1*)weight_file->Get("ratio_antiTop_jet_charge")->Clone();
  ratio_top_mass_lep_b = (TH1*)weight_file->Get("ratio_top_mass_lep_b")->Clone();
  ratio_ttbar_MT_mass_ratio_B_b = (TH1*)weight_file->Get("ratio_ttbar_MT_mass_ratio_B_b")->Clone();

  //std::cout << "Finished initializing MatchTester_ttbar_ll" << std::endl;
}

void MatchTester_ttbar_ll::evaluate() {
  if (this->evaluatedThisEvent) return;
  evaluatedThisEvent = true;

  //std::cout << "Evaluating MatchTester_ttbar_ll" << std::endl;
  
  if ( (*leptons)->size() < 2) return;
  if ( (*jets)->size() < 1) return;

  //std::cout << "There's at least two leptons" << std::endl;
  
  TLorentzVector lep1_vect;
  TLorentzVector lep2_vect;
  TLorentzVector jet1_vect;
  TLorentzVector jet2_vect;
  
  TLorentzVector jet1_vect_trans;
  TLorentzVector jet2_vect_trans;

  TLorentzVector lep1_B_vect;
  TLorentzVector lep2_b_vect;
  TLorentzVector B_b_vect;
  TLorentzVector B_b_vect_trans;

  double ratio_B = 0.0000000001;
  double ratio_b = 0.0000000001;
  double ratio_Bb = 0.0000000001;

  int bin = 0;
  int lep1_charge = 0;
  int lep2_charge = 0;

  //lep1 is the lepton from the top
  //lep2 is the lepton from the anti-top
  //B is the b from the top
  //b is the b from the anti-top
  
  //std::cout << "About to enter jet loop" << std::endl;

  for (unsigned int iLep1 = 0; iLep1 < (*leptons)->size(); iLep1++) {
    for (unsigned int iLep2 = 0; iLep2 < (*leptons)->size(); iLep2++) {
        if (iLep1 == iLep2) continue;
        
        lep1_vect.SetPtEtaPhiE((*leptons)->at(iLep1)->pt, (*leptons)->at(iLep1)->eta, (*leptons)->at(iLep1)->phi, (*leptons)->at(iLep1)->energy);
        lep2_vect.SetPtEtaPhiE((*leptons)->at(iLep2)->pt, (*leptons)->at(iLep2)->eta, (*leptons)->at(iLep2)->phi, (*leptons)->at(iLep2)->energy);
        lep1_charge = std::max((*leptons)->at(iLep1)->tkCharge, -1);
        lep2_charge = std::max((*leptons)->at(iLep2)->tkCharge, -1);
        if (lep1_charge != 1 || lep2_charge != -1) continue;
        
        for (unsigned int iJet1 = 0; iJet1 < (*jets)->size(); iJet1++) {    
          //std::cout << "iJet1 = " << iJet1 << std::endl;    
          jet1_vect.SetPtEtaPhiE((*jets)->at(iJet1).pt, (*jets)->at(iJet1).eta, (*jets)->at(iJet1).phi, (*jets)->at(iJet1).energy);
          jet1_vect_trans.SetPtEtaPhiE((*jets)->at(iJet1).pt, 0.0, (*jets)->at(iJet1).phi, (*jets)->at(iJet1).pt);

          //ratio_B
          lep1_B_vect = lep1_vect+jet1_vect;
          
          bin = std::max(1, std::min(ratio_top_2jet_CSV->GetNbinsX(), ratio_top_2jet_CSV->GetXaxis()->FindBin((*jets)->at(iJet1).btagCombinedSecVertex)) );
          ratio_B = ratio_top_2jet_CSV->GetBinContent(bin);
          if ((*jets)->size() == 3) {
            bin = std::max(1, std::min(ratio_top_3jet_CSV->GetNbinsX(), ratio_top_3jet_CSV->GetXaxis()->FindBin((*jets)->at(iJet1).btagCombinedSecVertex)) );
            ratio_B = ratio_top_3jet_CSV->GetBinContent(bin); }
          if ((*jets)->size() >= 4) {
            bin = std::max(1, std::min(ratio_top_4jet_CSV->GetNbinsX(), ratio_top_4jet_CSV->GetXaxis()->FindBin((*jets)->at(iJet1).btagCombinedSecVertex)) );
            ratio_B = ratio_top_4jet_CSV->GetBinContent(bin); }
          bin = std::max(1, std::min(ratio_top_jet_charge->GetNbinsX(), ratio_top_jet_charge->GetXaxis()->FindBin((*jets)->at(iJet1).charge)) );
          ratio_B *= ratio_top_jet_charge->GetBinContent(bin);
          bin = std::max(1, std::min(ratio_top_mass_lep_b->GetNbinsX(), ratio_top_mass_lep_b->GetXaxis()->FindBin(lep1_B_vect.M())) );
          ratio_B *= ratio_top_mass_lep_b->GetBinContent(bin);
          
          if (log(ratio_B) > branches["Match_ttbar_ll_B"].branchVal) {
            branches["Match_ttbar_ll_B"].branchVal = log(ratio_B);
            branches["Full_match_ttbar_ll_B"].branchVal = ((*jets)->at(iJet1).genPartonMotherId == 6);
            branches["ttbar_ll_B_CSV"].branchVal = (*jets)->at(iJet1).btagCombinedSecVertex;
            branches["ttbar_ll_B_charge"].branchVal = (*jets)->at(iJet1).charge;
            branches["ttbar_ll_top_mass_lep_B"].branchVal = lep1_B_vect.M();
          }
          
          //ratio_b
          lep2_b_vect = lep2_vect+jet1_vect;
          
          bin = std::max(1, std::min(ratio_top_2jet_CSV->GetNbinsX(), ratio_top_2jet_CSV->GetXaxis()->FindBin((*jets)->at(iJet1).btagCombinedSecVertex)) );
          ratio_b = ratio_top_2jet_CSV->GetBinContent(bin);
          if ((*jets)->size() == 3) {
            bin = std::max(1, std::min(ratio_top_3jet_CSV->GetNbinsX(), ratio_top_3jet_CSV->GetXaxis()->FindBin((*jets)->at(iJet1).btagCombinedSecVertex)) );
            ratio_b = ratio_top_3jet_CSV->GetBinContent(bin); }
          if ((*jets)->size() >= 4) {
            bin = std::max(1, std::min(ratio_top_4jet_CSV->GetNbinsX(), ratio_top_4jet_CSV->GetXaxis()->FindBin((*jets)->at(iJet1).btagCombinedSecVertex)) );
            ratio_b = ratio_top_4jet_CSV->GetBinContent(bin); }
          bin = std::max(1, std::min(ratio_antiTop_jet_charge->GetNbinsX(), ratio_antiTop_jet_charge->GetXaxis()->FindBin((*jets)->at(iJet1).charge)) );
          ratio_b *= ratio_antiTop_jet_charge->GetBinContent(bin);
          bin = std::max(1, std::min(ratio_top_mass_lep_b->GetNbinsX(), ratio_top_mass_lep_b->GetXaxis()->FindBin(lep2_b_vect.M())) );
          ratio_b *= ratio_top_mass_lep_b->GetBinContent(bin);
          
          if (log(ratio_b) > branches["Match_ttbar_ll_b"].branchVal) {
            branches["Match_ttbar_ll_b"].branchVal = log(ratio_b);
            branches["Full_match_ttbar_ll_b"].branchVal = ((*jets)->at(iJet1).genPartonMotherId == -6);
            branches["ttbar_ll_b_CSV"].branchVal = (*jets)->at(iJet1).btagCombinedSecVertex;
            branches["ttbar_ll_b_charge"].branchVal = (*jets)->at(iJet1).charge;
            branches["ttbar_ll_antiTop_mass_lep_b"].branchVal = lep2_b_vect.M();
          }
          
          for (unsigned int iJet2 = 0; iJet2 < (*jets)->size(); iJet2++) {
            if (iJet2 == iJet1) continue;
            //std::cout << "iJet2 = " << iJet1 << std::endl;      
            jet2_vect.SetPtEtaPhiE((*jets)->at(iJet2).pt, (*jets)->at(iJet2).eta, (*jets)->at(iJet2).phi, (*jets)->at(iJet2).energy);
            jet2_vect_trans.SetPtEtaPhiE((*jets)->at(iJet2).pt, 0, (*jets)->at(iJet2).phi, (*jets)->at(iJet2).pt);

            //ratio_Bb
            lep1_B_vect = lep1_vect+jet1_vect;
            lep2_b_vect = lep2_vect+jet2_vect;
            B_b_vect = jet1_vect+jet2_vect;
            B_b_vect_trans = jet1_vect_trans+jet2_vect_trans;
            
            bin = std::max(1, std::min(ratio_top_2jet_CSV->GetNbinsX(), ratio_top_2jet_CSV->GetXaxis()->FindBin((*jets)->at(iJet1).btagCombinedSecVertex)) );
            ratio_Bb = ratio_top_2jet_CSV->GetBinContent(bin);
            bin = std::max(1, std::min(ratio_top_2jet_CSV->GetNbinsX(), ratio_top_2jet_CSV->GetXaxis()->FindBin((*jets)->at(iJet2).btagCombinedSecVertex)) );
            ratio_Bb *= ratio_top_2jet_CSV->GetBinContent(bin);
            if ((*jets)->size() == 3) {
              bin = std::max(1, std::min(ratio_top_3jet_CSV->GetNbinsX(), ratio_top_3jet_CSV->GetXaxis()->FindBin((*jets)->at(iJet1).btagCombinedSecVertex)) );
              ratio_Bb = ratio_top_3jet_CSV->GetBinContent(bin);
              bin = std::max(1, std::min(ratio_top_3jet_CSV->GetNbinsX(), ratio_top_3jet_CSV->GetXaxis()->FindBin((*jets)->at(iJet2).btagCombinedSecVertex)) );
              ratio_Bb *= ratio_top_3jet_CSV->GetBinContent(bin); }
            if ((*jets)->size() >= 4) {
              bin = std::max(1, std::min(ratio_top_4jet_CSV->GetNbinsX(), ratio_top_4jet_CSV->GetXaxis()->FindBin((*jets)->at(iJet1).btagCombinedSecVertex)) );
              ratio_Bb = ratio_top_4jet_CSV->GetBinContent(bin);
              bin = std::max(1, std::min(ratio_top_4jet_CSV->GetNbinsX(), ratio_top_4jet_CSV->GetXaxis()->FindBin((*jets)->at(iJet2).btagCombinedSecVertex)) );
              ratio_Bb *= ratio_top_4jet_CSV->GetBinContent(bin); }
            bin = std::max(1, std::min(ratio_top_jet_charge->GetNbinsX(), ratio_top_jet_charge->GetXaxis()->FindBin((*jets)->at(iJet1).charge)) );
            ratio_Bb *= ratio_top_jet_charge->GetBinContent(bin);
            bin = std::max(1, std::min(ratio_antiTop_jet_charge->GetNbinsX(), ratio_antiTop_jet_charge->GetXaxis()->FindBin((*jets)->at(iJet2).charge)) );
            ratio_Bb *= ratio_antiTop_jet_charge->GetBinContent(bin);
            bin = std::max(1, std::min(ratio_top_mass_lep_b->GetNbinsX(), ratio_top_mass_lep_b->GetXaxis()->FindBin(lep1_B_vect.M())) );
            ratio_Bb *= ratio_top_mass_lep_b->GetBinContent(bin);
            bin = std::max(1, std::min(ratio_top_mass_lep_b->GetNbinsX(), ratio_top_mass_lep_b->GetXaxis()->FindBin(lep2_b_vect.M())) );
            ratio_Bb *= ratio_top_mass_lep_b->GetBinContent(bin);
            bin = std::max(1, std::min(ratio_ttbar_MT_mass_ratio_B_b->GetNbinsX(), ratio_ttbar_MT_mass_ratio_B_b->GetXaxis()->FindBin(B_b_vect_trans.M()/B_b_vect.M())) );
            ratio_Bb *= ratio_ttbar_MT_mass_ratio_B_b->GetBinContent(bin);
              
            if (log(ratio_Bb) > branches["Match_ttbar_ll_Bb"].branchVal) {
              branches["Match_ttbar_ll_Bb"].branchVal = log(ratio_Bb);
              branches["Full_match_ttbar_ll_Bb"].branchVal = ((*jets)->at(iJet1).genPartonMotherId == 6 && (*jets)->at(iJet2).genPartonMotherId == -6);
              branches["ttbar_ll_B_CSV"].branchVal = (*jets)->at(iJet1).btagCombinedSecVertex;
              branches["ttbar_ll_b_CSV"].branchVal = (*jets)->at(iJet2).btagCombinedSecVertex;
              branches["ttbar_ll_B_charge"].branchVal = (*jets)->at(iJet1).charge;
              branches["ttbar_ll_b_charge"].branchVal = (*jets)->at(iJet2).charge;
              branches["ttbar_ll_top_mass_lep_B"].branchVal = lep1_B_vect.M();
              branches["ttbar_ll_antiTop_mass_lep_b"].branchVal = lep2_b_vect.M();
              branches["ttbar_ll_ttbar_MT_mass_ratio_B_b"].branchVal = B_b_vect_trans.M()/B_b_vect.M();
            }
          } //end for iJet2
        } //end for iJet1
    } //end for iLep2
  } //end for iLep1

  //std::cout << "Finished jet loop" << std::endl;
  
  myVars.clear();
  
  for (typename map<TString, BranchInfo<double>>::iterator iBranch = branches.begin();
       iBranch != branches.end(); iBranch++) {
    myVars.push_back(iBranch->second);
  }
  
}

MatchTester_ttbar_ll::~MatchTester_ttbar_ll() {

  //Delete histograms BEFORE closing file
  
  delete ratio_top_2jet_CSV;
  delete ratio_top_3jet_CSV;
  delete ratio_top_4jet_CSV;
  delete ratio_top_jet_charge;
  delete ratio_antiTop_jet_charge;
  delete ratio_top_mass_lep_b;
  delete ratio_ttbar_MT_mass_ratio_B_b;
  weight_file->Close();

}


#endif
