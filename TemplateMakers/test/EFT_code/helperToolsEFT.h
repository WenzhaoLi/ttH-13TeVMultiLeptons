#ifndef HELPERTOOLSEFT_H_
#define HELPERTOOLSEFT_H_

#include <vector>
#include <typeinfo>

#include "TLorentzVector.h"
#include "ttH-13TeVMultiLeptons/TemplateMakers/src/classes.h"
#include "DataFormats/Math/interface/deltaR.h"

vector<ttH::GenParticle> sortGenParticles(vector<ttH::GenParticle> gen_particles) {
    vector<ttH::GenParticle> gen_list(gen_particles.begin(),gen_particles.end());
    std::sort(gen_list.begin(), gen_list.end(), [] (ttH::GenParticle a, ttH::GenParticle b) {return a.obj.Pt() > b.obj.Pt();});
    return gen_list;
}

template <typename T> int getSign(T val) {
    return (T(0) < val) - (val < T(0));
}

template <typename inObj> TLorentzVector setTlv(const inObj inputObj ) {
  TLorentzVector tlv; 
  tlv.SetPxPyPzE( inputObj.obj.px(), inputObj.obj.py(), inputObj.obj.pz(), inputObj.obj.E() );
  return tlv;
}

template <typename inObj1, typename inObj2> double getDeltaR(const inObj1 obj1, const inObj2 obj2) {
    TLorentzVector obj1_tlv = setTlv(obj1);
    TLorentzVector obj2_tlv = setTlv(obj2);
    if (obj1_tlv.Pt()*obj2_tlv.Pt() <= 1e-10) {
        return -1.;
    }

    return obj1_tlv.DeltaR( obj2_tlv );
}

vector<ttH::GenParticle> getGenLeptons(vector<ttH::GenParticle> gen_particles) {
    vector<ttH::GenParticle> leptons;
    for (auto &gen_particle: gen_particles) {
        // Must be prompt lepton
        if (gen_particle.isPromptFinalState || gen_particle.isDirectPromptTauDecayProductFinalState) {        
            int pdg_ID = fabs(gen_particle.pdgID);

            // Must be electron or muon
            if (pdg_ID != 11 && pdg_ID != 13) {
                continue;
            }

            leptons.push_back(gen_particle);
        }
    }
    return leptons;
}

template <typename inObj1, typename inObj2> vector<ttH::Lepton> getRecoLeptons(const vector<inObj1> collection_1, const vector<inObj2> collection_2) {
    vector<ttH::Lepton> lep_collection(collection_1.begin(),collection_1.end());
    lep_collection.insert(lep_collection.end(),collection_2.begin(),collection_2.end());
    std::sort(lep_collection.begin(), lep_collection.end(), [] (ttH::Lepton a, ttH::Lepton b) {
        return a.obj.Pt() > b.obj.Pt();
    });
    return lep_collection;
}

vector<ttH::GenParticle> getJets(vector<ttH::GenParticle> gen_particles) {
    vector<ttH::GenParticle> jets;
    for (auto &gen_particle: gen_particles) {
        jets.push_back(gen_particle);
    }
    return jets;
}

template <typename inObj> vector<inObj> applyPtCut(vector<inObj> particles, double pt_cut) {
    vector<inObj> selected_particles;
    for (auto &particle: particles) {
        if (particle.obj.Pt() < pt_cut) {
            continue;
        }
        selected_particles.push_back(particle);
    }
    return selected_particles;
}

template <typename inObj> vector<inObj> applyEtaCut(vector<inObj> particles, double eta_cut) {
    vector<inObj> selected_particles;
    for (auto &particle: particles) {
        if (particle.obj.Eta() > eta_cut) {
            continue;
        }
        selected_particles.push_back(particle);
    }
    return selected_particles;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
vector<ttH::GenParticle> getPromptParticles(vector<ttH::GenParticle> gen_particles) {
    vector<ttH::GenParticle> prompt_particles;
    for (auto &gen_particle: gen_particles) {
        if (!gen_particle.isPromptFinalState) {
            continue;
        }
        prompt_particles.push_back(gen_particle);
    }
    return prompt_particles;
}


TString getParticleName(int pdg_id) {
    pdg_id = abs(pdg_id);
    if (pdg_id == 11) {
        return TString("Electron");
    } else if (pdg_id == 12) {
        return TString("Neutrino (e)");
    } else if (pdg_id == 13) {
        return TString("Muon");
    } else if (pdg_id == 14) {
        return TString("Neutrino (m)");
    } else if (pdg_id == 15) {
        return TString("Tau");
    } else if (pdg_id == 16) {
        return TString("Neutrino (t)");
    } else if (pdg_id == 21) {
        return TString("Gluon");
    } else if (pdg_id == 22) {
        return TString("Photon");
    } else if (pdg_id == 24) {
        return TString("W");
    } else if (pdg_id == 1) {
        return TString("d Quark");
    } else if (pdg_id == 2) {
        return TString("u Quark");
    } else if (pdg_id == 3) {
        return TString("s Quark");
    } else if (pdg_id == 4) {
        return TString("c Quark");
    } else if (pdg_id == 5) {
        return TString("b Quark");
    } else if (pdg_id == 6) {
        return TString("t Quark");
    } else if (pdg_id == 111) {
        return TString("Pi0 Meson");
    } else if (pdg_id == 211) {
        return TString("Pi+ Meson");
    } else if (pdg_id == 310) {
        return TString("K0_S Meson");
    } else if (pdg_id == 413) {
        return TString("D*(2010)+ Meson");
    } else if (pdg_id == 421) {
        return TString("D0 Meson");
    } else if (pdg_id == 423) {
        return TString("D*(2007)0 Meson");
    } else if (pdg_id == 511) {
        return TString("B0 Meson");
    } else if (pdg_id == 523) {
        return TString("B+ Meson");
    } else if (pdg_id == 2212) {
        return TString("Proton");
    } else {
        return TString(std::to_string(pdg_id));
    }
}

vector<ttH::GenParticle> getParticlesByID(int id,vector<ttH::GenParticle> gen_particles) {
    vector<ttH::GenParticle> particles;
    for (auto &gen_particle: gen_particles) {
        int pdg_ID = abs(gen_particle.pdgID);
        if (pdg_ID != id) {
            continue;
        }
        particles.push_back(gen_particle);
    }
    return particles;
}

ttH::GenParticle getMotherParticle(ttH::GenParticle gen_particle, vector<ttH::GenParticle> gen_particles) {
    ttH::GenParticle mother_particle = gen_particles.at(gen_particle.mother);
    return mother_particle;
}

ttH::GenParticle getChildParticle(ttH::GenParticle gen_particle, vector<ttH::GenParticle> gen_particles, int child_choice = 0) {
    ttH::GenParticle child_particle;
    if (child_choice == 0) {
        child_particle = gen_particles.at(gen_particle.child0);
    } else {
        child_particle = gen_particles.at(gen_particle.child1);
    }
    return child_particle;
}

TString getIndentString(int depth = 0) {
    TString indent = "";
    for (int i = 0; i < depth; i++) {
        indent += "\t";
    }
    return indent;
}

void readParticleInfo(ttH::GenParticle gen_particle, int particle_index, int depth = 0) {
    TString indent = getIndentString(depth);

    TString particle_name = getParticleName(gen_particle.pdgID);

    cout << indent << "Particle: " << particle_name << endl;
    cout << indent << "\tIndex:   " << particle_index << endl;
    cout << indent << "\tpdgID:   " << gen_particle.pdgID << endl;
    cout << indent << "\tCharge:  " << gen_particle.charge << endl;
    cout << indent << "\tStatus:  " << gen_particle.status << endl;
    cout << indent << "\tPrompt:  " << gen_particle.isPromptFinalState << endl;
    cout << indent << "\tObj:     " << gen_particle.obj << endl;
    cout << indent << "\tM():     " << gen_particle.obj.M() << endl;
    cout << indent << "\tPt():    " << gen_particle.obj.Pt() << endl;
    cout << indent << "\tP():     " << gen_particle.obj.P() << endl;
    cout << indent << "\tE():     " << gen_particle.obj.E() << endl;
    cout << indent << "\tMother:  " << gen_particle.mother << endl;
    cout << indent << "\tGMother: " << gen_particle.grandmother << endl;
    cout << indent << "\tChild0:  " << gen_particle.child0 << endl;
    cout << indent << "\tChild1:  " << gen_particle.child1 << endl;
}

template <typename inObj> void readCollectionInfo(vector<inObj> collection, int depth = 0) {
    TString indent = getIndentString(depth);

    for (auto &particle: collection) {

        cout << indent << "Particle: " << getParticleName(particle.genPdgID) << endl;
        cout << indent << "\tMother:   " << getParticleName(particle.genMotherPdgID) << endl;
        cout << indent << "\tGMother:  " << getParticleName(particle.genGrandMotherPdgID) << endl;
    }
}

// Returns b-jets matched to their gen_particle counterparts
template <typename T> vector<T> getBJets(vector<ttH::GenParticle> gen_particles, vector<T> jets) {
    double deltaR_cut = 0.4;
    vector<T> b_jets;
    for (auto &gen_jet: jets) {
        int index = 0;
        for (auto &gen_particle: gen_particles) {
            double delta_r = getDeltaR(gen_jet,gen_particle);
            if (delta_r > deltaR_cut || delta_r < 0) {
                index += 1;
                continue;
            }

            if (fabs(gen_particle.pdgID) == 5) {
                b_jets.push_back(gen_jet);
                break;
            }
            index += 1;
        }
    }
    return b_jets;
}

// Determines if the given particle (index) corresponds directly to the mother's children
bool isOriginal(uint particle_index, ttH::GenParticle mother_particle) {
    if (particle_index == mother_particle.child0 || particle_index == mother_particle.child1) {
        return true;
    } else {
        return false;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

double getInvWMass(vector<ttH::GenParticle> gen_particles) {
    double inv_w_mass = -999.;
    vector<ttH::GenParticle> prompt_particles = getPromptParticles(gen_particles);
    for (uint i = 0; i < gen_particles.size(); i++) {
        uint particle_index = i;
        ttH::GenParticle gen_particle = gen_particles.at(i);

        if (fabs(gen_particle.pdgID) != 24) {
            // Only look at W's
            continue;
        }


        if (gen_particle.mother == 9999) {
            // The particle has no mother
            continue;
        }

        ttH::GenParticle mother_particle = getMotherParticle(gen_particle,gen_particles);
        
        if (fabs(mother_particle.pdgID) == 6) {
            // Ignore W's from decaying top
            continue;
        }

        if (!isOriginal(particle_index,mother_particle)) {
            // This particle wasn't the original one spawned from the mother
            continue;
        }

        inv_w_mass = gen_particle.obj.M();
        break;
    }

    return inv_w_mass;
}

double getNJets(vector<ttH::GenParticle> gen_jets, double jet_pt_cut, double jet_eta_cut) {
    vector<ttH::GenParticle> selected_jets = applyPtCut(gen_jets,jet_pt_cut);
    selected_jets = applyEtaCut(selected_jets,jet_eta_cut);

    return selected_jets.size();
}

#endif
/* HELPERTOOLSEFT */