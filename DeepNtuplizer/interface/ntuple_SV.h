/*
 * ntuple_SV.h
 *
 *  Created on: 13 Feb 2017
 *      Author: jkiesele
 */


#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_SV_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_SV_H_


#include "ntuple_content.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
// #include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"


class ntuple_SV : public ntuple_content {

    public:

        ntuple_SV(std::string prefix = "", double jetR = 0.4);
        ~ntuple_SV();

        void getInput(const edm::ParameterSet& iConfig);
        void initBranches(TTree*);
        void readEvent(const edm::Event& iEvent);
        void readSetup(const edm::EventSetup& iSetup);
        void initContainers();
        void clearContainers();
        void deleteContainers();

        // Use either of these functions
        // bool fillBranches(const pat::Jet&, const size_t& jetidx, const edm::View<pat::Jet>* coll = 0);
        bool fillBranches(const pat::Jet&, const size_t& jetidx, const edm::View<pat::Jet>* coll = 0, float EventTime = -1) { return true; }
        int fillBranches(bool applySelection, float EventTime = -1);

        void setTrackBuilderToken(const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord>& track_builder_token) {
            track_builder_token_ = track_builder_token;
        }

        void setGenParticlesToken(const edm::EDGetTokenT<pat::PackedGenParticleCollection>& genParticles_token) {
            genParticles_token_ = genParticles_token;
        }

        void setGenJetFlavourInfoToken(const edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection>& genJetFlavourInfo_token) {
            genJetFlavourInfo_token_ = genJetFlavourInfo_token;
        }

        void setSimTracksToken(const edm::EDGetTokenT<edm::SimTrackContainer>& simTracks_token) {
            simTracks_token_ = simTracks_token;
        }

        // void setPFCandToken(const edm::EDGetTokenT<pat::PackedCandidateCollection>& pf_cand_token) {
        //     pf_cand_token_ = pf_cand_token;
        // }

        // void setPFMCMatchToken(const edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>>& pf_mcmatch_token) {
        //     pf_mcmatch_token_ = pf_mcmatch_token;
        // }

        void setRecoTracksToken(const edm::EDGetTokenT<reco::TrackCollection>& recoTracks_token) {
            recoTracks_token_ = recoTracks_token;
        }

        void setTimeValueMapToken(const edm::EDGetTokenT<edm::ValueMap<float>> timeValueMap_token) {
            timeValueMap_token_ = timeValueMap_token;
        }

        void setTimeErrorMapToken(const edm::EDGetTokenT<edm::ValueMap<float>> timeErrorMap_token) {
            timeErrorMap_token_ = timeErrorMap_token;
        }

        void setTimeQualityMapToken(const edm::EDGetTokenT<edm::ValueMap<float>> timeQualityMap_token) {
            timeQualityMap_token_ = timeQualityMap_token;
        }

        // void setGenVertexToken(const edm::EDGetTokenT<TrackingVertexCollection>& genVertices_token) {
        //     genVertices_token_ = genVertices_token;
        // }

        void setPVsToken(const edm::EDGetTokenT<reco::VertexCollection>& PVs_token) {
            PVs_token_ = PVs_token;
        }

        void setInclusiveSVsToken(const edm::EDGetTokenT<reco::VertexCollection>& inclusiveSVs_token) {
            inclusiveSVs_token_ = inclusiveSVs_token;
        }

        void setInclusiveSVsMTDTimingToken(const edm::EDGetTokenT<reco::VertexCollection>& inclusiveSVsMTDTiming_token) {
            inclusiveSVsMTDTiming_token_ = inclusiveSVsMTDTiming_token;
        }

    private:

        double timeQualityCut_;
        double matchGVdR_;
        double jetPtMin_;
        double jetPtMax_;
        double jetAbsEtaMin_;
        double jetAbsEtaMax_;
        double genJetMatchdR_;

        // SV candidates
        int sv_num_;
        // float nsv_;
        std::string prefix_;

        edm::ESHandle<TransientTrackBuilder> builder_;
        edm::Handle<pat::PackedGenParticleCollection> genParticles_;
        edm::Handle<reco::JetFlavourInfoMatchingCollection> genJetFlavourInfo_;
        edm::Handle<edm::SimTrackContainer> simTracks_;
        // edm::Handle<pat::PackedCandidateCollection> pf_cand_;
        // edm::Handle<edm::Association<reco::GenParticleCollection>> pf_mcmatch_;
        edm::Handle<reco::TrackCollection> recoTracks_;
        edm::Handle<edm::ValueMap<float>> timeValueMap_;
        edm::Handle<edm::ValueMap<float>> timeErrorMap_;
        edm::Handle<edm::ValueMap<float>> timeQualityMap_;
        // edm::Handle<TrackingVertexCollection> genVertices_;
        edm::Handle<reco::VertexCollection> PVs_;
        edm::Handle<reco::VertexCollection> inclusiveSVs_;
        edm::Handle<reco::VertexCollection> inclusiveSVsMTDTiming_;

        edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> track_builder_token_;
        edm::EDGetTokenT<pat::PackedGenParticleCollection> genParticles_token_;
        edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> genJetFlavourInfo_token_;
        edm::EDGetTokenT<edm::SimTrackContainer> simTracks_token_;
        // edm::EDGetTokenT<pat::PackedCandidateCollection> pf_cand_token_;
        // edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>> pf_mcmatch_token_;
        edm::EDGetTokenT<reco::TrackCollection> recoTracks_token_;
        edm::EDGetTokenT<edm::ValueMap<float>> timeValueMap_token_;
        edm::EDGetTokenT<edm::ValueMap<float>> timeErrorMap_token_;
        edm::EDGetTokenT<edm::ValueMap<float>> timeQualityMap_token_;
        // edm::EDGetTokenT<TrackingVertexCollection> genVertices_token_;
        edm::EDGetTokenT<reco::VertexCollection> PVs_token_;
        edm::EDGetTokenT<reco::VertexCollection> inclusiveSVs_token_;
        edm::EDGetTokenT<reco::VertexCollection> inclusiveSVsMTDTiming_token_;

        static constexpr size_t max_sv_ = 1000;

        std::vector<int>* gp_strange_status_;
        std::vector<int>* gp_strange_pdgId_;
        std::vector<int>* gp_strange_nDaughters_;
        std::vector<int>* gp_strange_daughterPdgIds_;
        std::vector<float>* gp_strange_pt_;
        std::vector<float>* gp_strange_eta_;

        std::vector<float>* st_pt_;
        std::vector<float>* st_eta_;
        std::vector<float>* st_phi_;
        std::vector<int>* st_charge_;

        std::vector<float>* rt_pt_;
        std::vector<float>* rt_eta_;
        std::vector<float>* rt_phi_;
        std::vector<int>* rt_charge_;
        std::vector<float>* rt_chi2_;
        std::vector<float>* rt_ndof_;
        std::vector<float>* rt_chi2dof_;
        std::vector<float>* rt_dxy_;
        std::vector<float>* rt_dxyerr_;
        std::vector<float>* rt_d0_;
        std::vector<float>* rt_d0err_;
        std::vector<float>* rt_dz_;
        std::vector<float>* rt_dzerr_;

        std::vector<float>* gv_pt_;
        std::vector<float>* gv_x_;
        std::vector<float>* gv_y_;
        std::vector<float>* gv_z_;
        std::vector<float>* gv_motherEta_;
        std::vector<int>* gv_motherPdgId_;
        std::vector<int>* gv_nDaughters_;
        std::vector<float>* gv_daughterPt_;
        std::vector<float>* gv_daughterEta_;
        std::vector<int>* gv_daughterPdgId_;
        std::vector<float>* gv_dRtoClosestSV_;
        std::vector<float>* gv_dRtoAllSV_;
        std::vector<int>* gv_nTimesMatchedToSV_;
        std::vector<int>* gv_nTimesMatchedToSVinJet_;

        std::vector<float>* gv_matchSV_pt_;
        std::vector<float>* gv_matchSV_x_;
        std::vector<float>* gv_matchSV_y_;
        std::vector<float>* gv_matchSV_z_;
        std::vector<float>* gv_matchSV_motherEta_;
        std::vector<int>* gv_matchSV_motherPdgId_;
        std::vector<int>* gv_matchSV_nDaughters_;
        std::vector<float>* gv_matchSV_daughterPt_;
        std::vector<float>* gv_matchSV_daughterEta_;
        std::vector<int>* gv_matchSV_daughterPdgId_;
        std::vector<float>* gv_matchSV_dRtoClosestSV_;

        std::vector<float>* sv_x_;
        std::vector<float>* sv_y_;
        std::vector<float>* sv_z_;
        std::vector<float>* sv_pt_;
        std::vector<float>* sv_eta_;
        std::vector<float>* sv_phi_;
        std::vector<float>* sv_mass_;
        std::vector<float>* sv_energy_;
        std::vector<float>* sv_chi2_;
        std::vector<float>* sv_ndof_;
        std::vector<float>* sv_chi2dof_;
        std::vector<float>* sv_dxy_; // from PV
        std::vector<float>* sv_dxyerr_;
        std::vector<float>* sv_dxysig_;
        std::vector<float>* sv_dz_; // from PV
        std::vector<float>* sv_d3D_; // from GV
        std::vector<float>* sv_d3Derr_;
        std::vector<float>* sv_d3Dsig_;
        std::vector<int>* sv_nDaughters_;

        std::vector<float>* sv_matchGV_x_;
        std::vector<float>* sv_matchGV_y_;
        std::vector<float>* sv_matchGV_z_;
        std::vector<float>* sv_matchGV_pt_;
        std::vector<float>* sv_matchGV_eta_;
        std::vector<float>* sv_matchGV_phi_;
        std::vector<float>* sv_matchGV_mass_;
        std::vector<float>* sv_matchGV_energy_;
        std::vector<float>* sv_matchGV_chi2_;
        std::vector<float>* sv_matchGV_ndof_;
        std::vector<float>* sv_matchGV_chi2dof_;
        std::vector<float>* sv_matchGV_dxy_;
        std::vector<float>* sv_matchGV_dxyerr_;
        std::vector<float>* sv_matchGV_dxysig_;
        std::vector<float>* sv_matchGV_dz_;
        std::vector<float>* sv_matchGV_d3D_;
        std::vector<float>* sv_matchGV_d3Derr_;
        std::vector<float>* sv_matchGV_d3Dsig_;
        std::vector<float>* sv_matchGV_SVdRtoGV_;
        std::vector<int>* sv_matchGV_nDaughters_;
        std::vector<int>* sv_matchGV_GVmotherPdgId_;

        std::vector<float>* sv_matchGV_matchJet_x_;
        std::vector<float>* sv_matchGV_matchJet_y_;
        std::vector<float>* sv_matchGV_matchJet_z_;
        std::vector<float>* sv_matchGV_matchJet_pt_;
        std::vector<float>* sv_matchGV_matchJet_eta_;
        std::vector<float>* sv_matchGV_matchJet_phi_;
        std::vector<float>* sv_matchGV_matchJet_mass_;
        std::vector<float>* sv_matchGV_matchJet_energy_;
        std::vector<float>* sv_matchGV_matchJet_chi2_;
        std::vector<float>* sv_matchGV_matchJet_ndof_;
        std::vector<float>* sv_matchGV_matchJet_chi2dof_;
        std::vector<float>* sv_matchGV_matchJet_dxy_;
        std::vector<float>* sv_matchGV_matchJet_dxyerr_;
        std::vector<float>* sv_matchGV_matchJet_dxysig_;
        std::vector<float>* sv_matchGV_matchJet_dz_;
        std::vector<float>* sv_matchGV_matchJet_d3D_;
        std::vector<float>* sv_matchGV_matchJet_d3Derr_;
        std::vector<float>* sv_matchGV_matchJet_d3Dsig_;
        std::vector<float>* sv_matchGV_matchJet_SVdRtoGV_;
        std::vector<int>* sv_matchGV_matchJet_nJets_;
        std::vector<int>* sv_matchGV_matchJet_genJetHadFlav_;
        std::vector<int>* sv_matchGV_matchJet_nDaughters_;
        std::vector<int>* sv_matchGV_matchJet_GVmotherPdgId_;

        int nPV_;
        int nSV_;
        int nSVMTDTiming_;

        // Track timing information
        std::vector<float>* pv_trk_timeVal_;
        std::vector<float>* pv_trk_timeErr_;
        std::vector<float>* pv_trk_timeQual_;
        std::vector<float>* sv_trk_timeVal_;
        std::vector<float>* sv_trk_timeErr_;
        std::vector<float>* sv_trk_timeQual_;

        // Jets
        std::vector<float>* jet_pt_;
        std::vector<float>* jet_eta_;
        std::vector<float>* jet_phi_;
        std::vector<float>* jet_radius_;
        std::vector<int>* jet_hadFlav_;
        std::vector<int>* jet_partFlav_;
        std::vector<int>* jet_genHadFlav_;
        std::vector<int>* jet_genPartFlav_;

        std::vector<float>* jet_matchSV_pt_;
        std::vector<float>* jet_matchSV_eta_;
        std::vector<float>* jet_matchSV_phi_;
        std::vector<float>* jet_matchSV_radius_;
        std::vector<int>* jet_matchSV_nSV_;
        std::vector<int>* jet_matchSV_hadFlav_;
        std::vector<int>* jet_matchSV_partFlav_;
        std::vector<int>* jet_matchSV_genHadFlav_;
        std::vector<int>* jet_matchSV_genPartFlav_;

        std::vector<float>* jet_matchSV_matchGV_pt_;
        std::vector<float>* jet_matchSV_matchGV_eta_;
        std::vector<float>* jet_matchSV_matchGV_phi_;
        std::vector<float>* jet_matchSV_matchGV_radius_;
        std::vector<float>* jet_matchSV_matchGV_SVdRtoGV_;
        std::vector<int>* jet_matchSV_matchGV_nGV_;
        std::vector<int>* jet_matchSV_matchGV_hadFlav_;
        std::vector<int>* jet_matchSV_matchGV_partFlav_;
        std::vector<int>* jet_matchSV_matchGV_genHadFlav_;
        std::vector<int>* jet_matchSV_matchGV_genPartFlav_;

        // float sv_pt_[max_sv_];
        // float sv_eta_[max_sv_];
        // float sv_phi_[max_sv_];
        // float sv_e_[max_sv_];

        std::vector<float>* sv_etarel_;
        std::vector<float>* sv_phirel_;
        std::vector<float>* sv_deltaR_;

        // float sv_mass_[max_sv_];
        // float sv_phirel_[max_sv_];
        // float sv_etarel_[max_sv_];
        float sv_ntracks_[max_sv_];
        float sv_nMatchPFCand_[max_sv_];
        float sv_nUnmatchedTrk_[max_sv_];
        // float sv_chi2_[max_sv_];
        // float sv_ndf_[max_sv_];
        // float sv_normchi2_[max_sv_];
        // float sv_dxy_[max_sv_];
        // float sv_dxyerr_[max_sv_];
        // float sv_dxysig_[max_sv_];
        // float sv_d3d_[max_sv_];
        // float sv_d3derr_[max_sv_];
        // float sv_d3dsig_[max_sv_];
        float sv_costhetasvpv_[max_sv_];

        std::vector<float>* sv_enratio_;

        float sv_calo_frac_[max_sv_];
        float sv_hcal_frac_[max_sv_];
        float sv_puppiw_[max_sv_];
        // float sv_dz_[max_sv_];
        float sv_charge_sum_[max_sv_];

        std::vector<float>* sv_pfd2dval_;
        std::vector<float>* sv_pfd2dsig_;
        std::vector<float>* sv_pfd3dval_;
        std::vector<float>* sv_pfd3dsig_;

        std::vector<float>* sv_time_;

        static const reco::Vertex* spvp_;

        // Helper functions
        static Measurement1D vertexDxy(const reco::VertexCompositePtrCandidate& svcand, const reco::Vertex& pv);
        static Measurement1D vertexD3d(const reco::VertexCompositePtrCandidate& svcand, const reco::Vertex& pv);
        static float vertexDdotP(const reco::VertexCompositePtrCandidate& sv, const reco::Vertex& pv);
        static bool compareDxyDxyErr(const reco::VertexCompositePtrCandidate& sva, const reco::VertexCompositePtrCandidate& svb);
        // int findPFCandIdx(const pat::PackedCandidate& trk, const pat::PackedCandidateCollection& pcands);
        static float deltaR3D(float x1, float x2, float y1, float y2, float z1, float z2);
        static int genPartID(int pdgId);
};


#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_SV_H_ */
