#ifndef _DataDrivenFR_h
#define _DataDrivenFR_h

#include "ttHMultileptonAnalysis/TemplateMakers/interface/KinematicVariable.h"
#include "ttHMultileptonAnalysis/TemplateMakers/interface/BranchInfo.h"
#include <typeinfo>

template <class collectionType>
class DataDrivenFR: public KinematicVariable<double> {

public:

  //Store branch values so they are accessible to other classes
  vector<BranchInfo<double>> myVars;

  HelperLeptonCore * myHelper;
  collectionType **selCollection;
  unsigned int number_of_leptons;
  double working_point;
  string file_name_NP;
  string file_name_QF;

  TFile * weight_file_NP;
  TFile * weight_file_QF;
  TH2 * FR_NP_tight_mu; //FR for < 2 b-jets
  TH2 * FR_NP_tight2_mu; //FR for >= 2 b-jets
  TH2 * FR_NP_tight_el; //FR for < 2 b-jets
  TH2 * FR_NP_tight2_el; //FR for >= 2 b-jets
  TH2 * FR_QF_el; //Charge flip FR
  double FR_NP[6];
  double FR_QF[6];
  int QF_charge[6];
  
  DataDrivenFR(HelperLeptonCore *input_myHelper, collectionType **input_selCollection, int input_number_of_leptons,
               double input_working_point, string input_file_name_NP, string input_file_name_QF);

  ~DataDrivenFR();
  
  void evaluate();

};

template <class collectionType> 
DataDrivenFR<collectionType>::DataDrivenFR(HelperLeptonCore *input_myHelper, collectionType **input_selCollection, int input_number_of_leptons,
                                           double input_working_point, string input_file_name_NP, string input_file_name_QF):
  myHelper(input_myHelper), selCollection(input_selCollection), number_of_leptons(input_number_of_leptons),
  working_point(input_working_point), file_name_NP(input_file_name_NP), file_name_QF(input_file_name_QF) {

  this->resetVal = 1.0;
  
  branches["DataDrivenFR"] = BranchInfo<double>("DataDrivenFR");

  string directory = "../data/fakerate/";
  TString weight_file_name_NP = Form("%s%s.root", directory.c_str(), file_name_NP.c_str());
  weight_file_NP = TFile::Open(weight_file_name_NP);
  TString weight_file_name_QF = Form("%s%s.root", directory.c_str(), file_name_QF.c_str());
  weight_file_QF = TFile::Open(weight_file_name_QF);

  FR_NP_tight_mu = (TH2*)weight_file_NP->Get("FR_tight_mu")->Clone();
  FR_NP_tight2_mu = (TH2*)weight_file_NP->Get("FR_tight2_mu")->Clone();
  FR_NP_tight_el = (TH2*)weight_file_NP->Get("FR_tight_el")->Clone();
  FR_NP_tight2_el = (TH2*)weight_file_NP->Get("FR_tight2_el")->Clone();
  FR_QF_el = (TH2*)weight_file_QF->Get("QF_el_data")->Clone();

}

template <class collectionType> 
void DataDrivenFR<collectionType>::evaluate() {

  if (this->evaluatedThisEvent) return;
  evaluatedThisEvent = true;
  
  //--------

  BEANhelper * beanHelper = &(myHelper->bHelp);
  
  FR_NP[0] = 1.0; FR_NP[1] = 1.0; FR_NP[2] = 1.0; FR_NP[3] = 1.0; FR_NP[4] = 1.0; FR_NP[5] = 1.0; 
  FR_QF[0] = 1.0; FR_QF[1] = 1.0; FR_QF[2] = 1.0; FR_QF[3] = 1.0; FR_QF[4] = 1.0; FR_QF[5] = 1.0; 
  QF_charge[0] = 0; QF_charge[1] = 0; QF_charge[2] = 0; QF_charge[3] = 0; QF_charge[4] = 0; QF_charge[5] = 0; 
  int num_fail = 0;
  int num_leptons = 0;
  int num_electrons = 0;

  for (unsigned int iObj = 0; iObj < (*selCollection)->size(); iObj++) {
    if ( iObj < number_of_leptons ) {

      double lep_pt = ptr((*selCollection)->at(iObj))->pt;
      double lep_eta = abs(ptr((*selCollection)->at(iObj))->eta);

      if ( ptr((*selCollection)->at(iObj))->isMuon ) {
        if ( beanHelper->GetMuonLepMVA(*(BNmuon*)ptr((*selCollection)->at(iObj)),
                                           this->blocks->jetsForLepMVACollection) < working_point) {

          if (this->blocks->jetCollectionMediumCSV->size() < 2) {
            int pt_bin  = std::max(1, std::min(FR_NP_tight_mu->GetNbinsX(), FR_NP_tight_mu->GetXaxis()->FindBin(lep_pt)));
            int eta_bin  = std::max(1, std::min(FR_NP_tight_mu->GetNbinsY(), FR_NP_tight_mu->GetYaxis()->FindBin(lep_eta)));
            FR_NP[num_fail] = FR_NP_tight_mu->GetBinContent(pt_bin,eta_bin);            
          }
          else {
            int pt_bin  = std::max(1, std::min(FR_NP_tight2_mu->GetNbinsX(), FR_NP_tight2_mu->GetXaxis()->FindBin(lep_pt)));
            int eta_bin  = std::max(1, std::min(FR_NP_tight2_mu->GetNbinsY(), FR_NP_tight2_mu->GetYaxis()->FindBin(lep_eta)));
            FR_NP[num_fail] = FR_NP_tight_mu->GetBinContent(pt_bin,eta_bin);
          } 
          num_fail += 1;
        }
        else {
          QF_charge[num_leptons] = ((BNmuon*)ptr((*selCollection)->at(iObj)))->tkCharge;
          num_leptons += 1;
        }
      }
      else if ( ptr((*selCollection)->at(iObj))->isElectron ) {
        if ( beanHelper->GetElectronLepMVA(*(BNelectron*)ptr((*selCollection)->at(iObj)),
                                           this->blocks->jetsForLepMVACollection) < working_point) {

          if (this->blocks->jetCollectionMediumCSV->size() < 2) {
            int pt_bin  = std::max(1, std::min(FR_NP_tight_el->GetNbinsX(), FR_NP_tight_el->GetXaxis()->FindBin(lep_pt)));
            int eta_bin  = std::max(1, std::min(FR_NP_tight_el->GetNbinsY(), FR_NP_tight_el->GetYaxis()->FindBin(lep_eta)));
            FR_NP[num_fail] = FR_NP_tight_el->GetBinContent(pt_bin,eta_bin);
          }
          else {
            int pt_bin  = std::max(1, std::min(FR_NP_tight2_el->GetNbinsX(), FR_NP_tight2_el->GetXaxis()->FindBin(lep_pt)));
            int eta_bin  = std::max(1, std::min(FR_NP_tight2_el->GetNbinsY(), FR_NP_tight2_el->GetYaxis()->FindBin(lep_eta)));
            FR_NP[num_fail] = FR_NP_tight_el->GetBinContent(pt_bin,eta_bin);
          }
          num_fail += 1;
        }
        else {
          int pt_bin  = std::max(1, std::min(FR_QF_el->GetNbinsX(), FR_QF_el->GetXaxis()->FindBin(lep_pt)));
          int eta_bin  = std::max(1, std::min(FR_QF_el->GetNbinsY(), FR_QF_el->GetYaxis()->FindBin(lep_eta)));
          FR_QF[num_electrons] = FR_QF_el->GetBinContent(pt_bin,eta_bin);
          QF_charge[num_leptons] = ((BNelectron*)ptr((*selCollection)->at(iObj)))->tkCharge;
          num_electrons += 1;
          num_leptons += 1;
        }
      }
      else std::cout << "Lepton is neither muon nor electron" << std::endl;

    } //end if ( iObj < number_of_leptons )
  } //end loop over iObj

  if (number_of_leptons == 2) {
    //FR for events where both leptons fail
    if (FR_NP[1] != 1.0) {
      branches["DataDrivenFR"].branchVal = -FR_NP[0]*FR_NP[1]/((1-FR_NP[0])*(1-FR_NP[1]));
    }
    //FR for events where one lepton fails
    else if (FR_NP[0] != 1.0 ) {
      branches["DataDrivenFR"].branchVal = FR_NP[0]/(1-FR_NP[0]);
    }
    //FR for events where no leptons fail
    else branches["DataDrivenFR"].branchVal = 1.0;

    //FR for opposite-sign events with electrons
    if (QF_charge[0] * QF_charge[1] == -1 && FR_NP[0] == 1.0) {
      if (num_electrons == 0) branches["DataDrivenFR"].branchVal = 1.0;
      if (num_electrons == 1) branches["DataDrivenFR"].branchVal = FR_QF[0]; 
      if (num_electrons == 2) branches["DataDrivenFR"].branchVal = FR_QF[0] + FR_QF[1]; 
    }
  }

  //Clean out values from last event
  myVars.clear();
  
  for ( typename map<TString, BranchInfo<double>>::iterator iBranch = branches.begin();
        iBranch != branches.end(); iBranch++ ) {
    myVars.push_back(iBranch->second);
  }
  
}

template <class collectionType> 
DataDrivenFR<collectionType>::~DataDrivenFR() {

  //Delete histograms BEFORE closing file
  delete FR_NP_tight_mu;
  delete FR_NP_tight2_mu;
  delete FR_NP_tight_el;
  delete FR_NP_tight2_el;
  delete FR_QF_el;
  weight_file_NP->Close();
  weight_file_QF->Close();

}


#endif 