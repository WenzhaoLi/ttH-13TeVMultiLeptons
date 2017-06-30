#include "helperToolsEFT.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename T1, typename T2>bool pass_selection(
    vector<T1> input_leptons,
    vector<T2> input_jets,
    vector<ttH::GenParticle> gen_particles,
    double lep_pt_cut,
    double lep_eta_cut,
    double jet_pt_cut,
    double jet_eta_cut,
    double lep_pt_veto,
    int sign,
    uint lep_req,
    uint jet_req,
    uint b_jet_req,
    bool req_exact_lep,
    bool req_exact_jet
) {
    vector<T1> selected_leptons = applyPtCut(input_leptons,lep_pt_cut);
    selected_leptons = applyEtaCut(selected_leptons,lep_eta_cut);

    vector<T2> selected_jets = applyPtCut(input_jets,jet_pt_cut);
    selected_jets = applyEtaCut(selected_jets,jet_eta_cut);

    vector<ttH::GenParticle> cleaned_particles = applyPtCut(gen_particles,1.0);
    vector<T2> matched_b_jets = getBJets(cleaned_particles,selected_jets);

    if (req_exact_lep) {
        if (selected_leptons.size() != lep_req) {
            // Exact number of leptons
            return false;
        }
    } else {
        if (selected_leptons.size() < lep_req) {
            return false;
        }
    }

    if (req_exact_jet) {
        if (selected_jets.size() != jet_req) {
            // Exact number of jets
            return false;
        }
    } else {
        if (selected_jets.size() < jet_req) {
            return false;
        }
    }

    if (lep_req == 2) {
        // We are in the 2l category
        int sum_sign = getSign(selected_leptons[0].charge + selected_leptons[1].charge);
        if (sum_sign != sign) {
            return false;
        }
    } else if (lep_req == 3) {
        // We are in the 3l category
        int sum_sign = getSign(selected_leptons[0].charge + selected_leptons[1].charge + selected_leptons[2].charge);
        if (sum_sign != sign) {
            return false;
        }
    }

    if (matched_b_jets.size() < b_jet_req) {
        return false;
    }

    return true;
}