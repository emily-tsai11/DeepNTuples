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


class ntuple_SV : public ntuple_content {

    public:

        ntuple_SV(std::string prefix = "", double jetR = 0.4);
        ~ntuple_SV();

        void getInput(const edm::ParameterSet& iConfig);
        void initBranches(TTree*);
        void readEvent(const edm::Event& iEvent);
        void readSetup(const edm::EventSetup& iSetup);

        // Use either of these functions
        // bool fillBranches(const pat::Jet&, const size_t& jetidx, const edm::View<pat::Jet>* coll = 0);
        bool fillBranches(const pat::Jet&, const size_t& jetidx, const edm::View<pat::Jet>* coll = 0, float EventTime = -1);

        void setTrackBuilderToken(const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord>& track_builder_token) {
            track_builder_token_ = track_builder_token;
        }

        void setPFCandToken(const edm::EDGetTokenT<pat::PackedCandidateCollection>& pf_cand_token) {
            pf_cand_token_ = pf_cand_token;
        }

        void setLostTracksToken(const edm::EDGetTokenT<pat::PackedCandidateCollection>& lost_tracks_token) {
            lost_tracks_token_ = lost_tracks_token;
        }

        void setPFMCMatchToken(const edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>>& pf_mcmatch_token) {
            pf_mcmatch_token_ = pf_mcmatch_token;
        }

        void setLTMCMatchToken(const edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>>& lt_mcmatch_token) {
            lt_mcmatch_token_ = lt_mcmatch_token;
        }

    private:

        // SV candidates
        int sv_num_;
        float nsv_;
        std::string prefix_;

        edm::ESHandle<TransientTrackBuilder> builder_;
        edm::Handle<pat::PackedCandidateCollection> pf_cand_;
        edm::Handle<pat::PackedCandidateCollection> lost_tracks_;
        edm::Handle<edm::Association<reco::GenParticleCollection>> pf_mcmatch_;
        edm::Handle<edm::Association<reco::GenParticleCollection>> lt_mcmatch_;

        edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> track_builder_token_;
        edm::EDGetTokenT<pat::PackedCandidateCollection> pf_cand_token_;
        edm::EDGetTokenT<pat::PackedCandidateCollection> lost_tracks_token_;
        edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>> pf_mcmatch_token_;
        edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>> lt_mcmatch_token_;

        static constexpr size_t max_sv = 10;

        float sv_pt_[max_sv];
        float sv_eta_[max_sv];
        float sv_phi_[max_sv];
        float sv_e_[max_sv];
        float sv_etarel_[max_sv];
        float sv_phirel_[max_sv];
        float sv_deltaR_[max_sv];
        float sv_mass_[max_sv];
        // float sv_phirel_[max_sv];
        // float sv_etarel_[max_sv];
        float sv_ntracks_[max_sv];
        float sv_chi2_[max_sv];
        float sv_ndf_[max_sv];
        float sv_normchi2_[max_sv];
        float sv_dxy_[max_sv];
        float sv_dxyerr_[max_sv];
        float sv_dxysig_[max_sv];
        float sv_d3d_[max_sv];
        float sv_d3derr_[max_sv];
        float sv_d3dsig_[max_sv];
        float sv_costhetasvpv_[max_sv];
        float sv_enratio_[max_sv];
        float sv_hcal_frac_[max_sv];
        float sv_calo_frac_[max_sv];
        float sv_dz_[max_sv];
        float sv_pfd2dval_[max_sv];
        float sv_pfd2dsig_[max_sv];
        float sv_pfd3dval_[max_sv];
        float sv_pfd3dsig_[max_sv];
        float sv_puppiw_[max_sv];
        float sv_charge_sum_[max_sv];
        float sv_time_[max_sv];

        static const reco::Vertex* spvp_;

        // Helper functions
        static bool compareDxyDxyErr(const reco::VertexCompositePtrCandidate& sva, const reco::VertexCompositePtrCandidate& svb);
        static Measurement1D vertexDxy(const reco::VertexCompositePtrCandidate& svcand, const reco::Vertex& pv);
        static Measurement1D vertexD3d(const reco::VertexCompositePtrCandidate& svcand, const reco::Vertex& pv);
        static float vertexDdotP(const reco::VertexCompositePtrCandidate& sv, const reco::Vertex& pv);
        static bool compareVtxEta(TrackingVertex& gva, TrackingVertex& gvb);
        int findPFCandIdx(const pat::PackedCandidate& trk, const pat::PackedCandidateCollection& pcands);
        int findLostTrackIdx(const pat::PackedCandidate& trk, const pat::PackedCandidateCollection& lts);
};


#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_SV_H_ */
