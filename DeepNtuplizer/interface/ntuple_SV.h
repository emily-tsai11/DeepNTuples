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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include <map>


class GenVertex;

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

        void setGenParticlesToken(const edm::EDGetTokenT<reco::GenParticleCollection>& genParticles_token) {
            genParticles_token_ = genParticles_token;
        }

        void setGenJetFlavourInfoToken(const edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection>& genJetFlavourInfo_token) {
            genJetFlavourInfo_token_ = genJetFlavourInfo_token;
        }

        void setSimTracksToken(const edm::EDGetTokenT<edm::SimTrackContainer>& simTracks_token) {
            simTracks_token_ = simTracks_token;
        }

        // void setTPsToken(const edm::EDGetTokenT<TrackingParticleCollection>& TPs_token) {
        //     TPs_token_ = TPs_token;
        // }

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

        void setTrackMCMatchToken(const edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>>& trackMCMatch_token) {
            trackMCMatch_token_ = trackMCMatch_token;
        }

        // void setTVsToken(const edm::EDGetTokenT<TrackingVertexCollection>& TVs_token) {
        //     TVs_token_ = TVs_token;
        // }

        void setPVsToken(const edm::EDGetTokenT<reco::VertexCollection>& PVs_token) {
            PVs_token_ = PVs_token;
        }

        void setInclusiveSVsToken(const edm::EDGetTokenT<reco::VertexCollection>& inclusiveSVs_token) {
            inclusiveSVs_token_ = inclusiveSVs_token;
        }

        void setIVFClustersToken(const edm::EDGetTokenT<unsigned int> IVFclusters_token) {
            IVFclusters_token_ = IVFclusters_token;
        }

        void setInclusiveSVsMTDTimingToken(const edm::EDGetTokenT<reco::VertexCollection>& inclusiveSVsMTDTiming_token) {
            inclusiveSVsMTDTiming_token_ = inclusiveSVsMTDTiming_token;
        }

        void setIVFClustersMTDTimingToken(const edm::EDGetTokenT<unsigned int> IVFclustersMTDTiming_token) {
            IVFclustersMTDTiming_token_ = IVFclustersMTDTiming_token;
        }

    private:

        double absEtaMin_;
        double absEtaMax_;
        double genPartPtCut_;
        double genDauPtCut_;
        double genToSimTrackdR_;
        double genToSimTrackPt_;
        double trackPtCut_;
        double timeQualityCut_;
        double recoToGenTrackdR_;
        double recoToGenTrackPt_;
        double recoToSimTrackdR_;
        double recoToSimTrackPt_;
        double GVSimTrackMatchFrac_;
        double SVtoGVdR_;
        double SVGenPartMatchFrac_;
        double jetPtMin_;
        double jetPtMax_;
        double genJetMatchdR_;

        bool debug_;

        // SV candidates
        int sv_num_;
        // float nsv_;
        std::string prefix_;

        edm::ESHandle<TransientTrackBuilder> builder_;
        edm::Handle<reco::GenParticleCollection> genParticles_;
        edm::Handle<edm::SimTrackContainer> simTracks_;
        // edm::Handle<TrackingParticleCollection> TPs_;
        // edm::Handle<pat::PackedCandidateCollection> pf_cand_;
        // edm::Handle<edm::Association<reco::GenParticleCollection>> pf_mcmatch_;
        edm::Handle<reco::TrackCollection> recoTracks_;
        edm::Handle<edm::ValueMap<float>> timeValueMap_;
        edm::Handle<edm::ValueMap<float>> timeErrorMap_;
        edm::Handle<edm::ValueMap<float>> timeQualityMap_;
        edm::Handle<edm::Association<reco::GenParticleCollection>> trackMCMatch_;
        // edm::Handle<TrackingVertexCollection> TVs_;
        edm::Handle<reco::VertexCollection> PVs_;
        edm::Handle<reco::VertexCollection> inclusiveSVs_;
        edm::Handle<unsigned int> IVFclusters_;
        edm::Handle<reco::VertexCollection> inclusiveSVsMTDTiming_;
        edm::Handle<unsigned int> IVFclustersMTDTiming_;
        edm::Handle<reco::JetFlavourInfoMatchingCollection> genJetFlavourInfo_;

        edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> track_builder_token_;
        edm::EDGetTokenT<reco::GenParticleCollection> genParticles_token_;
        edm::EDGetTokenT<edm::SimTrackContainer> simTracks_token_;
        // edm::EDGetTokenT<TrackingParticleCollection> TPs_token_;
        // edm::EDGetTokenT<pat::PackedCandidateCollection> pf_cand_token_;
        // edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>> pf_mcmatch_token_;
        edm::EDGetTokenT<reco::TrackCollection> recoTracks_token_;
        edm::EDGetTokenT<edm::ValueMap<float>> timeValueMap_token_;
        edm::EDGetTokenT<edm::ValueMap<float>> timeErrorMap_token_;
        edm::EDGetTokenT<edm::ValueMap<float>> timeQualityMap_token_;
        edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>> trackMCMatch_token_;
        // edm::EDGetTokenT<TrackingVertexCollection> TVs_token_;
        edm::EDGetTokenT<reco::VertexCollection> PVs_token_;
        edm::EDGetTokenT<reco::VertexCollection> inclusiveSVs_token_;
        edm::EDGetTokenT<unsigned int> IVFclusters_token_;
        edm::EDGetTokenT<reco::VertexCollection> inclusiveSVsMTDTiming_token_;
        edm::EDGetTokenT<unsigned int> IVFclustersMTDTiming_token_;
        edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> genJetFlavourInfo_token_;

        std::vector<TString> n_;
        std::map<TString, std::vector<float>*> b_;

        // std::vector<TString> gp_collections_ = {
        //     "strange",
        // };

        // std::vector<TString> gp_branches_ = {
        //     "pt",
        //     "eta",
        //     "phi"
        //     "pdgId",
        //     "status",
        //     "nDaughters",
        //     "daughterPdgIds", // Special case
        // };

        std::vector<TString> trk_collections_ = {
            "st",           // SimTrack
            "all",          // All reco tracks
            "match_gp",     // Via trackMCMatch association map
            "pu",           // Pileup
            "gv",           // GenVertex (the daughters)
            "gvs",          // GenVertex (the daughters)
            "pv",           // PrimaryVertex
            "sv",           // SecondaryVertex
            "matchsv_gv",   // SV matched to GV
            "matchsv_gvs",  // SV matched to GV w/Sim
            "svt",          // SecondaryVertex w/MTD timing
            "matchsvt_gv",  // SV w/MTD timing matched to GV
            "matchsvt_gvs", // SV w/MTD timing matched to GV w/Sim
        };

        std::vector<TString> trk_branches_ = {
            "tval",
            "terr",
            "tsig",
            "tqual",
            "x",
            "y",
            "z",
            "pt",
            "pterr",
            "eta",
            "etaerr",
            "phi",
            "phierr",
            "dxy",
            "dxyerr",
            "dxysig",
            "dz",
            "dzerr",
            "dzsig",
            "d3d",
            "d3derr",
            "d3dsig",
            "d0",
            "d0err",
            "d0sig",
            "charge",
            "pdgId",
            "chi2",
            "ndof",
            "chi2dof",
        };

        std::vector<TString> vtx_collections_ {
            "gv",           // GenVertex
            "matchgv_sv",   // GV matched to SV
            "matchgv_svt",  // GV matched to SV w/MTD timing
            "gvs",          // GenVertex w/Sim
            "matchgvs_sv",  // GV w/Sim matched to SV
            "matchgvs_svt", // GV w/Sim matched to SV w/MTD timing
            "pv",           // PrimaryVertex
            "sv",           // SecondaryVertex
            "svt",          // SecondaryVertex w/MTD timing
        };

        std::vector<TString> vtx_branches_ {
            "timeavg",
            "timerange",
            "x",
            "y",
            "z",
            "pt",
            "eta",
            "phi",
            "dxy",
            "dxyerr",
            "dxysig",
            "dz",
            "dzerr",
            "dzsig",
            "d3d",
            "d3derr",
            "d3dsig",
            "chi2",
            "ndof",
            "chi2dof",
            "ntracks", // (daughters for GV)
            "motherPdgId",
            "matchdR",
            "nmatch",
        };

        std::vector<TString> jet_collections_ = {
            "all",
            "match_gen",
            // "match_gv",
            "match_gensv",
            "match_gensvgv",
            "match_gensvt",
            "match_gensvtgv",
        };

        std::vector<TString> jet_branches_ = {
            "pt",
            "eta",
            "phi",
            "hadflav",
            // "partflav",
            "nmatch",
        };

        std::vector<TString> evt_branches_ = {
            "nGV",
            "nGVs",
            "nPV",
            "nSV",
            "nSVt",
            "nClusters",
            "nClusterst",
        };

        static constexpr size_t max_sv_ = 10;

        // float sv_pt_[max_sv_];
        // float sv_eta_[max_sv_];
        // float sv_phi_[max_sv_];
        float sv_e_[max_sv_];

        std::vector<float>* sv_etarel_;
        std::vector<float>* sv_phirel_;
        std::vector<float>* sv_deltaR_;

        float sv_mass_[max_sv_];
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

        // std::vector<float>* sv_time_;

        static const reco::Vertex* spvp_;

        // Quality cuts
        template <class P> static bool goodGenParticle(const P* gp, float ptCut, float etaCut);
        // static bool goodGenVertex(const GenVertex& gv, float motherPtCut, float dauPtCut, float etaCut);
        template <class T> static bool goodTrack(const T& trkRef, const edm::ValueMap<float>& timeValueMap, 
                const edm::ValueMap<float>& timeErrorMap, const edm::ValueMap<float>& timeQualityMap,
                float trackPtCut, float trackEtaCut, float timeQualityCut);
        static bool goodRecoVertex(const reco::Vertex& rv, const edm::ValueMap<float>& timeValueMap,
                const edm::ValueMap<float>& timeErrorMap, const edm::ValueMap<float>& timeQualityMap,
                float trackPtCut, float trackEtaCut, float timeQualityCut);
        static bool goodJet(const pat::Jet& jet, float ptMin, float ptMax, float etaCut);
        static int genPartID(int pdgId);

        // Matching
        static bool isPileupTrack(const SimTrack& st);
        template <class P> static bool matchGenToSimTrack(const P* gt, const SimTrack& st,
                float drCut, float ptCut);
        template <class T, class P> static bool matchRecoToGenTrack(const T& trkRef, const P& gt,
                float drCut, float ptCut);
        template <class T> static int matchRecoToSimTrack(const T& trkRef,
                const edm::SimTrackContainer& simTracks, float drCut, float ptCut);
        static bool matchGenToSimVertex(const GenVertex& gv, const edm::SimTrackContainer& simTracks,
                float matchFrac, float trkDrCut, float trkPtCut);
        static bool matchRecoToGenVertex(const reco::Vertex& v, const GenVertex& gv,
                float matchFrac, float trkDrCut, float trkPtCut,
                const edm::ValueMap<float>& timeValueMap,
                const edm::ValueMap<float>& timeErrorMap,
                const edm::ValueMap<float>& timeQualityMap,
                float trackPtCut, float trackEtaCut, float timeQualityCut);

        // Helper functions
        // int findPFCandIdx(const pat::PackedCandidate& trk, const pat::PackedCandidateCollection& pcands);
        template <class T> static float trackD3d(const T& trkRef);
        template <class T> static float trackD3dErr(const T& trkRef);
        static float vertexPt(const reco::Vertex& sv);
        static float vertexEta(const reco::Vertex& sv);
        static float vertexPhi(const reco::Vertex& sv);
        static float vertexDxy(const GenVertex& gv, const reco::Vertex& v);
        static float vertexDxyErr(const GenVertex& gv, const reco::Vertex& v);
        static float vertexD3d(const GenVertex& gv, const reco::Vertex& v);
        static float vertexD3dErr(const GenVertex& gv, const reco::Vertex& v);
        static Measurement1D vertexDxy(const reco::Vertex& sv, const reco::Vertex& pv);
        static Measurement1D vertexD3d(const reco::Vertex& sv, const reco::Vertex& pv);
        static Measurement1D vertexDxy(
                const reco::VertexCompositePtrCandidate& cand, const reco::Vertex& pv);
        static Measurement1D vertexD3d(
                const reco::VertexCompositePtrCandidate& cand, const reco::Vertex& pv);
        static bool gvCompareDxyDxySig(const GenVertex& gva, const GenVertex& gvb);
        static bool svCompareDxyDxySig(const reco::Vertex& sva, const reco::Vertex& svb);
        static bool candCompareDxyDxySig(
                const reco::VertexCompositePtrCandidate& sva, const reco::VertexCompositePtrCandidate& svb);
        static float vertexDdotP(const reco::VertexCompositePtrCandidate& sv, const reco::Vertex& pv);
        static float jetRadius(float jet_radius, const pat::Jet& jet);
};


#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_SV_H_ */
