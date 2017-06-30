#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <map>
#include "TSystem.h"
#include "TH1.h"
#include "TChain.h"
#include "TString.h"
#include "TH1D.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include "helperToolsEFT.h"
#include "loadEFTSamples.h"
#include "selectionEFT.h"

/////////////////////////////////////////
///
/// usage: root -l makeEFTSelectionTree.C+
///
/////////////////////////////////////////

void printProgress(int current_index, int total_entries) {
    if (current_index % max(int(total_entries/100.),1) == 0) {
        float fraction = 100.*current_index/total_entries;
        cout << int(fraction) << " % processed" << endl;
    }
}

void run_it(TString sample_name, vector<TString> bin_names, TString output_file_name, int job_no) {
    FileLoader file_loader(sample_name,job_no);
    TChain* chain = file_loader.chain;

    int total_count = file_loader.total_count;
    int chainentries = chain->GetEntries();   
    int last_entry = chainentries;
    int first_entry = 0;

    cout << "job_no: " << job_no << endl;
    cout << "chainentries: " << chainentries << endl;
    cout << "first entry: " << first_entry << endl;
    cout << "last entry: " << last_entry << endl;

    vector<ttH::GenParticle> *pruned_genParticles_intree = 0;
    vector<ttH::GenParticle> *genJets_intree = 0;

    vector<ttH::Electron> *tight_electrons_intree = 0;
    vector<ttH::Muon>     *tight_muons_intree = 0;
    vector<ttH::Electron> *preselected_electrons_intree = 0;
    vector<ttH::Muon>     *preselected_muons_intree = 0;
    vector<ttH::Jet>      *preselected_jets_intree = 0;

    double mcwgt_intree = -999.;


    chain->SetBranchAddress("mcwgt", &mcwgt_intree);
    chain->SetBranchAddress("pruned_genParticles",  &pruned_genParticles_intree);
    chain->SetBranchAddress("genJets", &genJets_intree);

    chain->SetBranchAddress("tight_electrons",       &tight_electrons_intree);
    chain->SetBranchAddress("tight_muons",           &tight_muons_intree);
    chain->SetBranchAddress("preselected_electrons", &preselected_electrons_intree);
    chain->SetBranchAddress("preselected_muons",     &preselected_muons_intree);
    chain->SetBranchAddress("preselected_jets",      &preselected_jets_intree);

    Int_t cachesize = 250000000;   //500 MBytes
    //  Int_t cachesize = 1024000000;   //1 GBytes
    chain->SetCacheSize(cachesize);

    TH1D* inv_W_mass_hist = new TH1D("W Mass Distribution","W Mass Distribution",150,-1,100);
    TH1D* WZ_njets_hist = new TH1D("N Jet Distribution","N Jet Distribution",16,-0.5,15.5);
    TH1D* ZZ_njets_hist = new TH1D("N Jet Distribution","N Jet Distribution",16,-0.5,15.5);
    TH1D* WW_njets_hist = new TH1D("N Jet Distribution","N Jet Distribution",16,-0.5,15.5);

    inv_W_mass_hist->GetXaxis()->SetTitle("Mass (GeV)");
    inv_W_mass_hist->GetYaxis()->SetTitle("Counts");

    TString cutflow_var = "Pt";
    //vector<double> cutflow_points {10,12,14,16,18,20,22,24,25,26,28,30};
    vector<double> cutflow_points {20};
    std::map<TString,std::map<int,double> > cutflow_counts;
    for (auto &bin_name: bin_names) {
        for (auto &cutflow_point: cutflow_points) {
            cutflow_counts[bin_name][cutflow_point] = 0;
        }
    }
    
    double lep_eta_cut = 2.8;
    double jet_pt_cut = 25;
    double jet_eta_cut = 2.8;

    //last_entry = 5;
    for (int i = first_entry; i <= last_entry; i++) {
        printProgress(i,last_entry);
        chain->GetEntry(i);

        // Sample specific checks
        if (sample_name == "ttW") {
            // Plot the invariant W mass distribution
            double invariant_W_mass = getInvWMass(*pruned_genParticles_intree);
            inv_W_mass_hist->Fill(invariant_W_mass);
        } else if (sample_name == "WZ_to3lnu") {
            // Plot the N jet distributions
            double njets = getNJets(*genJets_intree,jet_pt_cut,jet_eta_cut);
            WZ_njets_hist->Fill(njets,1);
        } else if (sample_name == "ZZ_to4l") {
            // Plot the N jet distributions
            double njets = getNJets(*genJets_intree,jet_pt_cut,jet_eta_cut);
            ZZ_njets_hist->Fill(njets,1);
        } else if (sample_name == "WW_2l2nu") {
            // Plot the N jet distributions
            double njets = getNJets(*genJets_intree,jet_pt_cut,jet_eta_cut);
            WW_njets_hist->Fill(njets,1);
        }

        //vector<ttH::GenParticle> input_leptons = getGenLeptons(*pruned_genParticles_intree);
        //vector<ttH::Lepton> input_leptons = getRecoLeptons(*tight_electrons_intree,*tight_muons_intree);
        vector<ttH::Lepton> input_leptons = getRecoLeptons(*preselected_electrons_intree,*preselected_muons_intree);

        //vector<ttH::GenParticle> input_jets = *genJets_intree;
        vector<ttH::Jet> input_jets = *preselected_jets_intree;


        double lep_pt_veto = 10.0
        int lep_req = 2;
        int jet_req = 6;
        int b_jet_req = 1;
        int sign = 0;                   //NOTE: SS+ --> +1, SS- --> -1, OS --> 0
        bool req_exact_lep = true;
        bool req_exact_jet = true;

        // Fill acceptance tables for each bin_name
        for (auto &bin_name: bin_names) {
            for (auto &cutflow_point: cutflow_points) {
                double lep_pt_cut = cutflow_point;
                bool passes = false;
                if (bin_name == "2los_6jets") {
                    lep_req = 2;
                    jet_req = 6;
                    sign = 0;
                } else if (bin_name == "2lss_p_6jets") {
                    lep_req = 2;
                    jet_req = 6;
                    sign = 1;
                } else if (bin_name == "2lss_m_6jets") {
                    lep_req = 2;
                    jet_req = 6;
                    sign = -1;
                } else if (bin_name == "3l_ppm_4jets") {
                    lep_req = 3;
                    jet_req = 4;
                    sign = 1;
                } else if (bin_name == "3l_mmp_4jets") {
                    lep_req = 3;
                    jet_req = 4;
                    sign = -1;
                } else if (bin_name == "4l_2jets") {
                    lep_req = 4;
                    jet_req = 2;
                    sign = 0;
                } else if (bin_name == "full_WW_2l2nu") {
                    lep_req = 0;
                    jet_req = 0;
                    sign = 0;
                    req_exact_lep = false;
                    req_exact_jet = false;
                    b_jet_req = 0;
                } else if (bin_name == "1l_0jets_0bjets_WW_2l2nu") {
                    lep_req = 1;
                    jet_req = 0;
                    sign = 0;
                    b_jet_req = 0;
                } else if (bin_name == "2l_0jets_0bjets_WW_2l2nu") {
                    lep_req = 2;
                    jet_req = 0;
                    sign = 0;
                    b_jet_req = 0;
                } else if (bin_name == "2l_1jets_0bjets_WW_2l2nu") {
                    lep_req = 2;
                    jet_req = 1;
                    sign = 0;
                    b_jet_req = 0;
                } else if (bin_name == "2l_2jets_0bjets_WW_2l2nu") {
                    lep_req = 2;
                    jet_req = 2;
                    sign = 0;
                    b_jet_req = 0;
                } else if (bin_name == "2l_3jets_0bjets_WW_2l2nu") {
                    lep_req = 2;
                    jet_req = 3;
                    sign = 0;
                    b_jet_req = 0;
                } else if (bin_name == "2l_4jets_0bjets_WW_2l2nu") {
                    lep_req = 2;
                    jet_req = 4;
                    sign = 0;
                    b_jet_req = 0;
                } else if (bin_name == "2l_5jets_0bjets_WW_2l2nu") {
                    lep_req = 2;
                    jet_req = 5;
                    sign = 0;
                    b_jet_req = 0;
                } else if (bin_name == "2l_6jets_0bjets_WW_2l2nu") {
                    lep_req = 2;
                    jet_req = 6;
                    sign = 0;
                    b_jet_req = 0;
                } else if (bin_name == "2l_6jets_1bjets_WW_2l2nu") {
                    lep_req = 2;
                    jet_req = 6;
                    sign = 0;
                    b_jet_req = 1;
                }

                passes = pass_selection(
                    input_leptons,
                    input_jets,
                    *pruned_genParticles_intree,
                    lep_pt_cut,
                    lep_eta_cut,
                    jet_pt_cut,
                    jet_eta_cut,
                    lep_pt_veto,
                    sign,
                    lep_req,
                    jet_req,
                    b_jet_req,
                    req_exact_lep,
                    req_exact_jet
                );

                if (passes) {
                    cutflow_counts[bin_name][cutflow_point] += 1*mcwgt_intree;
                }
            }
        }
    }

    cout << "Total Counts: " << total_count << endl;

    TString public_path = "/afs/crc.nd.edu/user/a/awightma/Public/";
    //if (sample_name == "ttW") {
    //    TString output_root_file = "W_mass_dist.root";
    //    cout << "Output File: " << public_path + output_root_file << endl;
    //    TFile *out_f = new TFile(public_path + output_root_file, "RECREATE");
    //    inv_W_mass_hist->Write();
    //    out_f->Close();
    //}

    //if (sample_name == "WZ_to3lnu") {
    //    TString output_root_file = "WZ_Njets_dist.root";
    //    cout << "Output File: " << public_path + output_root_file << endl;
    //    TFile *out_f = new TFile(public_path + output_root_file,"RECREATE");
    //    WZ_njets_hist->Write();
    //    out_f->Close();
    //} else if (sample_name == "ZZ_to4l") {
    //    TString output_root_file = "ZZ_Njets_dist.root";
    //    cout << "Output File: " << public_path + output_root_file << endl;
    //    TFile *out_f = new TFile(public_path + output_root_file,"RECREATE");
    //    ZZ_njets_hist->Write();
    //    out_f->Close();
    //} else if (sample_name == "WW_2l2nu") {
    //    TString output_root_file = "WW_Njets_dist.root";
    //    cout << "Output File: " << public_path + output_root_file << endl;
    //    TFile *out_f = new TFile(public_path + output_root_file,"RECREATE");
    //    WW_njets_hist->Write();
    //    out_f->Close();
    //}

    cout << "Saving file... " << output_file_name << endl;
    std::ofstream ofile(output_file_name);
    ofile << "Total: " << total_count << endl;
    ofile << "Sample: " << sample_name << endl;
    for (auto &bin_name: bin_names) {
        ofile << "Bin Name: " << bin_name << endl;
        for (auto &cutflow_point: cutflow_points) {
            ofile << "\t" << cutflow_var << " " << cutflow_point << ": " << cutflow_counts[bin_name][cutflow_point] << endl;
        }
    }
    ofile.close();
}
/*
    2los_6jets
    2lss_p_6jets
    2lss_m_6jets
    3l_ppm_4jets
    3l_mmp_4jets
    4l_2jets

    5l_p_0bjets
    5l_m_0bjets
    6l_0bjets
*/

void makeEFTSelectionTree(TString sample_name="ttW",int job_no=-1) {
    TString output_dir = "/afs/crc.nd.edu/user/a/awightma/Public/";
    TString output_file_name = output_dir + sample_name + "_acceptance_table.txt";

    //vector<TString> bin_names {};
    vector<TString> bin_names {"2los_6jets","2lss_p_6jets","2lss_m_6jets","3l_ppm_4jets","3l_mmp_4jets","4l_2jets"};
    //vector<TString> bin_names {
    //    "full_WW_2l2nu",
    //    "1l_0jets_0bjets_WW_2l2nu",
    //    "2l_0jets_0bjets_WW_2l2nu",
    //    "2l_1jets_0bjets_WW_2l2nu",
    //    "2l_2jets_0bjets_WW_2l2nu",
    //    "2l_3jets_0bjets_WW_2l2nu",
    //    "2l_4jets_0bjets_WW_2l2nu",
    //    "2l_5jets_0bjets_WW_2l2nu",
    //    "2l_6jets_0bjets_WW_2l2nu",
    //    "2l_6jets_1bjets_WW_2l2nu"
    //};

    run_it(sample_name,bin_names,output_file_name,job_no);
}