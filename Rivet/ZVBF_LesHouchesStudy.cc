/*
 * ZVBF_LesHouchesStudy.cc
 *
 *  Created on: Dec. 8, 2017
 *      Author: Les Houches MC working group
 */

// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

#include <iostream>
#include <algorithm>
#include <string> 

static const float ZMASS = 91.1876;

namespace Rivet {

class ZVBF_LesHouchesStudy: public Analysis {
public:
    enum CutFlow {
        inSample,
        passingLepAcceptance,
        passingLepVeto,
        passingZconstraint,
        passing2jselection,
        passingAll
    };
    float totalEvents = 0;
    float sumSelectedWeights = 0;
    float sumWeightsLeps = 0;
    float sumWeights2j = 0;
    float sumWeightsZ = 0;
    float sumWeightsVeto = 0;
    float selectedEvents = 0;

    ZVBF_LesHouchesStudy() : Analysis("ZVBF_LesHouchesStudy"){};

    void init() {
        double lepConeSize = 0.1;
        double lepMaxEta = 2.5;
  
        Cut lepton_cut   = (Cuts::abseta < lepMaxEta);
  
        // Initialise and register projections
        FinalState fs(-2.5,2.5,0.0*GeV);
        FinalState fsm(-5,5,0.0*GeV);
        addProjection(fs, "FS");
        addProjection(fsm, "FSM");
  
        ChargedLeptons charged_leptons(fs);
        IdentifiedFinalState photons(fs);
        photons.acceptIdPair(PID::PHOTON);
  
        PromptFinalState prompt_leptons(charged_leptons);
        prompt_leptons.acceptMuonDecays(true);
        prompt_leptons.acceptTauDecays(false);
  
        PromptFinalState prompt_photons(photons);
        prompt_photons.acceptMuonDecays(true);
        prompt_photons.acceptTauDecays(false);
  
        DressedLeptons dressed_leptons = DressedLeptons(prompt_photons, prompt_leptons, lepConeSize, lepton_cut, /*cluster*/ true, /*useDecayPhotons*/ true);
        addProjection(dressed_leptons, "DressedLeptons");

        addProjection(FastJets(FinalState(-4.7, 4.7), FastJets::ANTIKT, 0.4),"jets");

        // Bookkeeping variables
        hists1D_["yieldByChannel"] = bookHisto1D("yieldByChannel", 80, -40, 40);
        bookChannelHist("final_xsec", 1, 0, 10);
        bookChannelHist("cut_flow", 7, 0, 7);

        bookChannelHist("lep1_Pt", 100, 0, 1000);
        bookChannelHist("lep1_Eta", 20, -2.5, 2.5);
        bookChannelHist("lep2_Pt", 50, 0, 250);
        bookChannelHist("lep2_Eta", 20, -2.5, 2.5);
        bookChannelHist("nJets", 10, 2, 10);

        bookChannelHist("jet1_Pt", 100, 0, 1000);
        bookChannelHist("jet1_Eta", 20, -4.7, 4.7);
        bookChannelHist("jet2_Pt", 100, 0, 1000);
        bookChannelHist("jet2_Eta", 20, -4.7, 4.7);
        bookChannelHist("jet3_Pt", 60, 0, 300);
        bookChannelHist("jet3_Eta", 20, -4.7, 4.7);

        // Composite variables
        bookChannelHist("ZMass", 32, 75, 107);
        bookChannelHist("ZPt", 100, 0, 1000);
        bookChannelHist("ZEta", 20, -3, 3);

        // VBS variables
        bookChannelHist("mjj", 100, 0, 4000);
        bookChannelHist("dEtajj", 64, -8, 8);
        bookChannelHist("dPhijj", 64, -8, 8);
        bookChannelHist("dRjj", 60, 0, 15);
        bookChannelHist("etasj3", 40, -5, 5);
        bookChannelHist("etasll", 40, -5, 5);
        hists1D_["veto_j3"] = bookHisto1D("veto_j3", 60, 0, 300);
    }

    // Change to true for CMS tight fiducial definition
    bool fullFiducial = false;
    void analyze(const Event& event) {
        double weight = event.weight();
        totalEvents++;
        channelHists_["cut_flow"].fill(CutFlow::inSample, weight, 0);
    
        float leppt_cut = 25;
        Particles leptons = applyProjection<DressedLeptons>(event, "DressedLeptons").particlesByPt(leppt_cut*GeV);

        if (leptons.size() < 2) {
            vetoEvent;
        }

        int sumPids = leptons.at(0).pdgId()+leptons.at(1).pdgId();
        int sumAbsPids = std::abs(leptons.at(0).pdgId())+std::abs(leptons.at(1).pdgId());

        // Unique identifier = sum(abs(pdgid)) * 1 if W-, -1 if W+
        int chanId = sumAbsPids;

        sumWeightsLeps += weight;

        // Channel division doesn't work here because of > 3 lepton events
        channelHists_["cut_flow"].fill(CutFlow::passingLepAcceptance, weight, 0);

        if (leptons.size() != 2) {
            vetoEvent;
        }
        sumWeightsVeto += weight;
        channelHists_["cut_flow"].fill(CutFlow::passingLepVeto, weight, 0);

        if (sumPids != 0) {
            vetoEvent;
            std::cerr << "WARNING: Found expected event without sum(pdgID) = " 
                      << sumPids << std::endl;
        }

        FourMomentum zcand = leptons.at(0).momentum() + leptons.at(1).momentum();

        if (std::abs(zcand.mass() - ZMASS) > 15)
            vetoEvent;

        sumWeightsZ += weight;
        channelHists_["cut_flow"].fill(CutFlow::passingZconstraint, weight, chanId);

        Jets jets;
        foreach (const Jet& jet, applyProjection<FastJets>(event, "jets").jetsByPt(30*GeV) ) {
            bool isolated = true;
            foreach (const Particle& lepton, leptons) {
                if (deltaR(lepton, jet) < 0.4) {
                    isolated = false;
                    break;
                }
            }
            if (isolated)
                jets.push_back(jet);
        }

        if (jets.size() < 2) {  
            vetoEvent;
        }
        channelHists_["cut_flow"].fill(CutFlow::passing2jselection, weight, chanId);
     
        FourMomentum dijet_system = jets.at(0).momentum() + jets.at(1).momentum();
        
        float mjj = dijet_system.mass();
        float dEtajj = jets.at(0).momentum().eta() - jets.at(1).momentum().eta();
        float dPhijj = jets.at(0).momentum().phi() - jets.at(1).momentum().phi();
        float dRjj = std::sqrt(dEtajj*dEtajj+dPhijj*dPhijj);
        float etasll = zcand.eta() - 0.5*(jets.at(0).momentum().eta()  + jets.at(1).momentum().eta());
        float etasj3 = jets.size() > 2 ? jets.at(2).eta() - 0.5*(jets.at(0).momentum().eta()  + jets.at(1).momentum().eta()) : -999;
        sumWeights2j += weight;

        if (mjj < 500 || std::abs(dEtajj) < 2.5)
            vetoEvent;

        sumSelectedWeights += weight;
        selectedEvents++;
        channelHists_["cut_flow"].fill(CutFlow::passingAll, weight, chanId);
        channelHists_["final_xsec"].fill(1, weight, chanId);
        hists1D_["yieldByChannel"]->fill(chanId, weight);

        // Primitive variables
        channelHists_["lep1_Pt"].fill(leptons.at(0).pt(), weight, chanId);
        channelHists_["lep1_Eta"].fill(leptons.at(0).eta(), weight, chanId);
        channelHists_["lep2_Pt"].fill(leptons.at(1).pt(), weight, chanId);
        channelHists_["lep2_Eta"].fill(leptons.at(1).eta(), weight, chanId);
        channelHists_["nJets"].fill(jets.size(), weight, chanId);

        channelHists_["jet1_Pt"].fill(jets.at(0).pt(), weight, chanId);
        channelHists_["jet1_Eta"].fill(jets.at(0).eta(), weight, chanId);
        channelHists_["jet2_Pt"].fill(jets.at(1).pt(), weight, chanId);
        channelHists_["jet2_Eta"].fill(jets.at(1).eta(), weight, chanId);

        // Composite variables
        channelHists_["ZMass"].fill(zcand.mass(), weight, chanId);
        channelHists_["ZPt"].fill(zcand.pt(), weight, chanId);
        channelHists_["ZEta"].fill(zcand.eta(), weight, chanId);

        // VBS variables
        channelHists_["mjj"].fill(mjj, weight, chanId);
        channelHists_["dEtajj"].fill(dEtajj, weight, chanId);
        channelHists_["dRjj"].fill(dRjj, weight, chanId);
        channelHists_["dPhijj"].fill(dPhijj, weight, chanId);
        channelHists_["etasll"].fill(etasll, weight, chanId);
        if (etasj3 != -999) {
            channelHists_["etasj3"].fill(etasj3, weight, chanId);
            channelHists_["jet3_Pt"].fill(jets.at(2).pt(), weight, chanId);
            channelHists_["jet3_Eta"].fill(jets.at(2).eta(), weight, chanId);
        }
    }

    void finalize() {
        std::cout << "Finalizing..." << endl;  

        float efficiency= selectedEvents/totalEvents; 
        double xsec = crossSection();

        std::cout << "SumOfWeights() processed = " << sumOfWeights() << endl;
        std::cout << "Selected sumWeights = " << sumSelectedWeights << endl;
        std::cout << "Sum weights passing lepton selection = " << sumWeightsLeps << endl;
        std::cout << "    cross section = " << sumWeightsLeps*xsec/sumOfWeights() << endl;
        std::cout << "Sum weights passing lepton veto = " << sumWeightsVeto << endl;
        std::cout << "    cross section = " << sumWeightsVeto*xsec/sumOfWeights() << endl;
        std::cout << "Sum weights passing Z selection = " << sumWeightsZ << endl;
        std::cout << "    cross section = " << sumWeightsZ*xsec/sumOfWeights() << endl;
        std::cout << "Sum weights passing 2j selection = " << sumWeights2j << endl;
        std::cout << "    cross section = " << sumWeights2j*xsec/sumOfWeights() << endl;
        std::cout << std::endl << "Selected Events = " << selectedEvents << ", Total= " << totalEvents << endl;
        std::cout << "Efficiency = " << efficiency << endl;    
        
        std::cout << "Initial cross section was: " << xsec << std::endl;
        std::cout << "Fiducial cross section is: " << xsec*sumSelectedWeights/sumOfWeights() << std::endl;
        
        // If you're using scale weights, using the sum of events is correct since the
        // cross section is associated with the central value (always 1) and the weights 
        // should not be unitary
        auto j3hist = channelHists_["jet3_Pt"].GetHist(0);
        for (size_t i = 0; i < j3hist->numBins(); i++) {
            double totalInt = j3hist->integral();
            double vetoEff = totalInt != 0 ? j3hist->integralTo(i)/totalInt : 0;
            hists1D_["veto_j3"]->fillBin(i, vetoEff);
        }

        removeEmptyHists();
        float sumWeights = sumOfWeights();
        for (const auto& hist : hists1D_) {
            if (hist.first.find("veto") != std::string::npos)
                continue;
            scale(hist.second, xsec / sumWeights);
        }
        for (auto& chanHist : channelHists_) {
            chanHist.second.Scale(xsec, sumWeights);
        }

    }        

private:
    // Label = sum{abs(pdgid)}*sign(sum{pdgid})
    // For now we don't distiguish between e/m states
    const std::map <int, std::string> channels_ = {
        { 26 , "Z_ee"},
        { 22 , "Z_mm"},
    };
    class ByChannelHist {
        private:
            // Needed to access functions in Analysis (e.g. Analysis::scale)
            Analysis* analysis_;
            std::map<int, Histo1DPtr> hists_;

        public:
            ByChannelHist() : analysis_(NULL) {};
            ByChannelHist(Analysis* analysis) : analysis_(analysis) {};

            void SetHist(Histo1DPtr hist, int channel) {
                hists_[channel] = hist;
            };

            Histo1DPtr GetHist(int channel) {
                if (hists_.find(channel) != hists_.end())
                    return hists_[channel];
                return nullptr;
            };

            void Scale(double xsec, double sumWeights) {
                std::set<Histo1DPtr> unique_hists;
                for (const auto& hist : hists_)
                    unique_hists.insert(hist.second);
                for (const auto& hist : unique_hists)
                    analysis_->scale(hist, xsec / sumWeights);
            };

            std::vector<Histo1DPtr> GetEmpty() {
                std::vector<Histo1DPtr> emptyHists = {};
                for (const auto& hist : hists_) {
                    if (hist.second->integral() == 0)
                        emptyHists.push_back(hist.second);
                }
                return emptyHists;
            };

            void fill(double value, double weight, int channel) {
                if (hists_.find(channel) == hists_.end()) {
                    throw std::runtime_error("Attempt to fill hist for invalid channel ID " 
                            + std::to_string(channel));
                }

                hists_[channel]->fill(value, weight);
                if (channel != 0)
                    hists_[0]->fill(value, weight);
            };

            void fillBin(int bin, double value, int channel) {
                if (hists_.find(channel) == hists_.end()) {
                    throw std::runtime_error("Attempt to fill hist for invalid channel ID " 
                            + std::to_string(channel));
                }

                hists_[channel]->fillBin(bin, value);
                std::cout << "Bin is " << bin << " Value is " << value << std::endl;
                if (channel != 0)
                    hists_[0]->fillBin(bin, value);
            };
    };
    std::map<std::string, Histo1DPtr> hists1D_;
    //std::map<std::string, Scatter2DPtr> scatters_;
    std::map<std::string, ByChannelHist> channelHists_;

    void bookChannelHist (const std::string histname, int bins, float binlow, float binhigh) {
        channelHists_[histname] = ByChannelHist(this);
        // Central hist
        channelHists_[histname].SetHist(bookHisto1D(histname, bins, binlow, binhigh), 0);

        std::map<std::string, int> uniqueNames;
        for (auto& channel : channels_) {
            const std::string chanHistName = channel.second + "_" + histname;
            if (uniqueNames.find(chanHistName) == uniqueNames.end()) {
                channelHists_[histname].SetHist(bookHisto1D(chanHistName, bins, binlow, binhigh), channel.first);
                uniqueNames[chanHistName] = channel.first;
            }
            // Get the existing histogram if channels share a histogram (i.e., they're combined into one)
            else
                channelHists_[histname].SetHist(
                    channelHists_[histname].GetHist(uniqueNames[chanHistName]), channel.first);
        }
    };
    void removeEmptyHists () {
        for (auto& hist : hists1D_) {
            if (hist.second->integral() == 0)
                removeAnalysisObject(hist.second);
        }
        for (auto& chanHist : channelHists_) {
            for (auto& hist : chanHist.second.GetEmpty()) {
                removeAnalysisObject(hist);
            }
        }
    };
};

DECLARE_RIVET_PLUGIN (ZVBF_LesHouchesStudy);
}
