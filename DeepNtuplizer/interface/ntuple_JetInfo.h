/*
 * ntuple_JetInfo.h
 *
 *  Created on: 13 Feb 2017
 *      Author: jkiesele
 */


#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_JETINFO_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_JETINFO_H_


#include "ntuple_content.h"
#include "TRandom3.h"
#include <map>
#include <string>


/*
 * For global jet info such as eta, pt, gen info
 */
class ntuple_JetInfo : public ntuple_content {

    public:

        ntuple_JetInfo();
        ~ntuple_JetInfo();

        void getInput(const edm::ParameterSet& iConfig);
        void initBranches(TTree*);
        void readEvent(const edm::Event& iEvent);
        void initContainers();
        void clearContainers();
        void deleteContainers();

        // Use either of these functions
        // bool fillBranches(const pat::Jet&, const size_t& jetidx, const edm::View<pat::Jet>* coll = 0);
        bool fillBranches(const pat::Jet&, const size_t& jetidx, const edm::View<pat::Jet>* coll = 0, float EventTime = -1);
        void fillBranches(bool applySelection);

        void setAxis2Token(edm::EDGetTokenT<edm::ValueMap<float>> axis2Token) {
            axis2Token_ = axis2Token;
        }

        void setMultToken(edm::EDGetTokenT<edm::ValueMap<int>> multToken) {
            multToken_ = multToken;
        }

        void setPtDToken(edm::EDGetTokenT<edm::ValueMap<float>> ptDToken) {
            ptDToken_ = ptDToken;
        }

        void setQglToken(edm::EDGetTokenT<edm::ValueMap<float>> qglToken) {
            qglToken_ = qglToken;
        }

        void setGenJetMatchReclusterToken(edm::EDGetTokenT<edm::Association<reco::GenJetCollection>> genJetMatchReclusterToken) {
            genJetMatchReclusterToken_ = genJetMatchReclusterToken;
        }

        void setGenJetMatchWithNuToken(edm::EDGetTokenT<edm::Association<reco::GenJetCollection>> genJetMatchWithNuToken) {
            genJetMatchWithNuToken_ = genJetMatchWithNuToken;
        }

        void setGenParticlesToken(edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken) {
            genParticlesToken_ = genParticlesToken;
        }

        void setPUInfoToken(edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puInfoToken) {
            puInfoToken_ = puInfoToken;
        }

        void setMuonsToken(edm::EDGetTokenT<pat::MuonCollection> muonsToken) {
            muonsToken_ = muonsToken;
        }

        void setElectronsToken(edm::EDGetTokenT<pat::ElectronCollection> electronsToken) {
            electronsToken_ = electronsToken;
        }

        void setUseHerwigCompatibleMatching(const bool use) {
            useherwcompat_matching_ = use;
        }

        void setIsHerwig(const bool use) {
            isherwig_ = use;
        }

    // private:

        double jetPtMin_;
        double jetPtMax_;
        double jetAbsEtaMin_;
        double jetAbsEtaMax_;

        // Quark gluon likelihood
        edm::EDGetTokenT<edm::ValueMap<float>> qglToken_;
        edm::EDGetTokenT<edm::ValueMap<float>> ptDToken_;
        edm::EDGetTokenT<edm::ValueMap<float>> axis2Token_;
        edm::EDGetTokenT<edm::ValueMap<int>> multToken_;

        edm::Handle<edm::ValueMap<float>> qglHandle_;
        edm::Handle<edm::ValueMap<float>> ptDHandle_;
        edm::Handle<edm::ValueMap<float>> axis2Handle_;
        edm::Handle<edm::ValueMap<int>> multHandle_;

        edm::EDGetTokenT<edm::Association<reco::GenJetCollection>> genJetMatchReclusterToken_;
        edm::EDGetTokenT<edm::Association<reco::GenJetCollection>> genJetMatchWithNuToken_;

        edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
        edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puInfoToken_;

        edm::EDGetTokenT<pat::MuonCollection> muonsToken_;
        edm::EDGetTokenT<pat::ElectronCollection> electronsToken_;

        edm::Handle<edm::Association<reco::GenJetCollection>> genJetMatchRecluster_;
        edm::Handle<edm::Association<reco::GenJetCollection>> genJetMatchWithNu_;

        edm::Handle<reco::GenParticleCollection> genParticlesHandle_;
        edm::Handle<std::vector <PileupSummaryInfo>> PUInfo_;

        edm::Handle<pat::MuonCollection> muonsHandle_;
        edm::Handle<pat::ElectronCollection> electronsHandle_;

        TRandom3 TRandom_;
        float gluonReduction_;

        std::vector <reco::GenParticle> neutrinosLepB_;
        std::vector <reco::GenParticle> neutrinosLepB_C_;

        std::vector<reco::GenParticle> gToBB_;
        std::vector<reco::GenParticle> gToCC_;
        std::vector<reco::GenParticle> alltaus_;

        std::vector<reco::GenParticle> Bhadron_;
        std::vector<reco::GenParticle> Bhadron_daughter_;

        bool useherwcompat_matching_;
        bool isherwig_;

        // Branches

        // Labels (MC truth)

        // Classification
        // int isB_;
        // int isGBB_;
        // int isBB_;
        // int isC_;
        // int isGCC_;
        // int isCC_;
        // int isUD_;
        // int isS_;
        // int isG_;
        // int isUndefined_;
        // float genDecay_;
        // int isLeptonicB_;
        // int isLeptonicB_C_;
        // int isTau_;
        std::vector<int>* isB_;
        std::vector<int>* isGBB_;
        std::vector<int>* isBB_;
        std::vector<int>* isC_;
        std::vector<int>* isGCC_;
        std::vector<int>* isCC_;
        std::vector<int>* isUD_;
        std::vector<int>* isS_;
        std::vector<int>* isG_;
        std::vector<int>* isUndefined_;
        std::vector<int>* isLeptonicB_;
        std::vector<int>* isLeptonicB_C_;
        std::vector<int>* isTau_;
        std::vector<float>* genDecay_;

        // Truth labeling with fallback to physics definition for light/gluon/undefined of standard flavor definition
        // int isPhysB_;
        // int isPhysGBB_;
        // int isPhysBB_;
        // int isPhysC_;
        // int isPhysGCC_;
        // int isPhysCC_;
        // int isPhysUD_;
        // int isPhysS_;
        // int isPhysG_;
        // int isPhysUndefined_;
        // int isPhysLeptonicB_;
        // int isPhysLeptonicB_C_;
        // int isPhysTau_;
        std::vector<int>* isPhysB_;
        std::vector<int>* isPhysGBB_;
        std::vector<int>* isPhysBB_;
        std::vector<int>* isPhysC_;
        std::vector<int>* isPhysGCC_;
        std::vector<int>* isPhysCC_;
        std::vector<int>* isPhysUD_;
        std::vector<int>* isPhysS_;
        std::vector<int>* isPhysG_;
        std::vector<int>* isPhysUndefined_;
        std::vector<int>* isPhysLeptonicB_;
        std::vector<int>* isPhysLeptonicB_C_;
        std::vector<int>* isPhysTau_;

        // Global variables
        int nJets_;
        float npv_;
        float npv_0_z_;
        float PU_rho_;
        float ntrueInt_;
        float rho_;
        unsigned int event_no_;
        // unsigned int jet_no_;

        // Jet variables
        // float jet_pt_;
        // float jet_corr_pt_;
        // float jet_eta_;
        // float jet_phi_;
        // float jet_mass_;
        // float jet_energy_;
        // float jet_looseId_;
        std::vector<float>* jet_pt_;
        std::vector<float>* jet_corr_pt_;
        std::vector<float>* jet_eta_;
        std::vector<float>* jet_phi_;
        std::vector<float>* jet_mass_;
        std::vector<float>* jet_energy_;
        std::vector<float>* jet_looseId_;

        // Quark/gluon
        // float jet_qgl_;
        // float QG_ptD_;
        // float QG_axis2_;
        // float QG_mult_;
        std::vector<float>* jet_qgl_;
        std::vector<float>* QG_ptD_;
        std::vector<float>* QG_axis2_;
        std::vector<float>* QG_mult_;

        // float y_multiplicity_;
        // float y_charged_multiplicity_;
        // float y_neutral_multiplicity_;
        // float y_ptD_;
        // float y_axis1_;
        // float y_axis2_;
        // float y_pt_dr_log_;
        std::vector<float>* y_multiplicity_;
        std::vector<float>* y_charged_multiplicity_;
        std::vector<float>* y_neutral_multiplicity_;
        std::vector<float>* y_ptD_;
        std::vector<float>* y_axis1_;
        std::vector<float>* y_axis2_;
        std::vector<float>* y_pt_dr_log_;

        static constexpr std::size_t max_num_lept_ = 5;
        // int muons_isLooseMuon_[max_num_lept_];
        // int muons_isTightMuon_[max_num_lept_];
        // int muons_isSoftMuon_[max_num_lept_];
        // int muons_isHighPtMuon_[max_num_lept_];
        // float muons_pt_[max_num_lept_];
        // float muons_relEta_[max_num_lept_];
        // float muons_relPhi_[max_num_lept_];
        // float muons_energy_[max_num_lept_];
        // float electrons_pt_[max_num_lept_];
        // float electrons_relEta_[max_num_lept_];
        // float electrons_relPhi_[max_num_lept_];
        // float electrons_energy_[max_num_lept_];

        // int muons_number_ = 0;
        // int electrons_number_ = 0;

        // These vectors have the same number of entries per event as the number of jets
        std::vector<int>* muons_number_;
        std::vector<int>* electrons_number_;

        // These vectors will probably have more entries per event than the number of jets
        std::vector<int>* muons_jetIdx_; // The index of the jet (that the muon at entry i is matched to)
        std::vector<int>* muons_isLooseMuon_;
        std::vector<int>* muons_isTightMuon_;
        std::vector<int>* muons_isSoftMuon_;
        std::vector<int>* muons_isHighPtMuon_;
        std::vector<float>* muons_pt_;
        std::vector<float>* muons_relEta_;
        std::vector<float>* muons_relPhi_;
        std::vector<float>* muons_energy_;

        std::vector<int>* electrons_jetIdx_; // The index of the jet (that the electron at entry i is matched to)
        std::vector<float>* electrons_pt_;
        std::vector<float>* electrons_relEta_;
        std::vector<float>* electrons_relPhi_;
        std::vector<float>* electrons_energy_;

        // Regressions pt, deta, dphi (truth)
        // float gen_pt_;
        // float gen_pt_Recluster_;
        // float gen_pt_WithNu_;
        // float Delta_gen_pt_;
        // float Delta_gen_pt_Recluster_;
        // float Delta_gen_pt_WithNu_;
        std::vector<float>* gen_pt_;
        std::vector<float>* gen_pt_Recluster_;
        std::vector<float>* gen_pt_WithNu_;
        std::vector<float>* Delta_gen_pt_;
        std::vector<float>* Delta_gen_pt_Recluster_;
        std::vector<float>* Delta_gen_pt_WithNu_;

        std::map<std::string, float> discriminators_;
};


#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_JETINFO_H_ */
