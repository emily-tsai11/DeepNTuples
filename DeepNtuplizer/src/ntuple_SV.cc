/*
 * ntuple_SV.cc
 *
 *  Created on: 13 Feb 2017
 *      Author: jkiesele
 */


#include "../interface/ntuple_SV.h"

// For IVF
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TLorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "TVector3.h"

#include <algorithm>


class SVTrackInfoBuilder {

    public:

        SVTrackInfoBuilder(edm::ESHandle<TransientTrackBuilder>& build) :
            builder_(build),
            trackMomentum_(0),
            trackEta_(0),
            trackEtaRel_(0),
            trackPtRel_(0),
            trackPPar_(0),
            trackDeltaR_(0),
            trackPtRatio_(0),
            trackPParRatio_(0),
            trackSip2dVal_(0),
            trackSip2dSig_(0),
            trackSip3dVal_(0),
            trackSip3dSig_(0),
            trackJetDecayLen_(0),
            trackJetDistVal_(0),
            trackJetDistSig_(0),
            ttrack_(0) {}

        void buildTrackInfo(const pat::PackedCandidate* PackedCandidate_, const math::XYZVector& jetDir, GlobalVector refjetdirection, const reco::Vertex& pv) {

            TVector3 jetDir3(jetDir.x(), jetDir.y(), jetDir.z());
            if(!PackedCandidate_->hasTrackDetails()) {
                TVector3 trackMom3(
                    PackedCandidate_->momentum().x(),
                    PackedCandidate_->momentum().y(),
                    PackedCandidate_->momentum().z()
                );
                trackMomentum_ = PackedCandidate_->p();
                trackEta_ = PackedCandidate_->eta();
                trackEtaRel_ = reco::btau::etaRel(jetDir, PackedCandidate_->momentum());
                trackPtRel_ = trackMom3.Perp(jetDir3);
                trackPPar_ = jetDir.Dot(PackedCandidate_->momentum());
                trackDeltaR_ = reco::deltaR(PackedCandidate_->momentum(), jetDir);
                trackPtRatio_ = trackMom3.Perp(jetDir3) / PackedCandidate_->p();
                trackPParRatio_ = jetDir.Dot(PackedCandidate_->momentum()) / PackedCandidate_->p();
                trackSip2dVal_ = 0.0;
                trackSip2dSig_ = 0.0;
                trackSip3dVal_ = 0.0;
                trackSip3dSig_ = 0.0;
                trackJetDecayLen_ = 0.0;
                trackJetDistVal_ = 0.0;
                trackJetDistSig_ = 0.0;
                return;
            }

            const reco::Track& PseudoTrack = PackedCandidate_->pseudoTrack();

            reco::TransientTrack transientTrack;
            transientTrack = builder_->build(PseudoTrack);
            Measurement1D meas_ip2d = IPTools::signedTransverseImpactParameter(transientTrack, refjetdirection, pv).second;
            Measurement1D meas_ip3d = IPTools::signedImpactParameter3D(transientTrack, refjetdirection, pv).second;
            Measurement1D jetdist = IPTools::jetTrackDistance(transientTrack, refjetdirection, pv).second;
            Measurement1D decayl = IPTools::signedDecayLength3D(transientTrack, refjetdirection, pv).second;
            math::XYZVector trackMom = PseudoTrack.momentum();
            double trackMag = std::sqrt(trackMom.Mag2());
            TVector3 trackMom3(trackMom.x(), trackMom.y(), trackMom.z());

            trackMomentum_ = std::sqrt(trackMom.Mag2());
            trackEta_ = trackMom.Eta();
            trackEtaRel_ = reco::btau::etaRel(jetDir, trackMom);
            trackPtRel_ = trackMom3.Perp(jetDir3);
            trackPPar_ = jetDir.Dot(trackMom);
            trackDeltaR_ = reco::deltaR(trackMom, jetDir);
            trackPtRatio_ = trackMom3.Perp(jetDir3) / trackMag;
            trackPParRatio_ = jetDir.Dot(trackMom) / trackMag;

            trackSip2dVal_ = meas_ip2d.value();
            trackSip2dSig_ = meas_ip2d.significance();
            trackSip3dVal_ = meas_ip3d.value();
            trackSip3dSig_ = meas_ip3d.significance();

            trackJetDecayLen_ = decayl.value();
            trackJetDistVal_ = jetdist.value();
            trackJetDistSig_ = jetdist.significance();

            ttrack_ = transientTrack;
        }

        const float& getTrackDeltaR() const { return trackDeltaR_; }
        const float& getTrackEta() const { return trackEta_; }
        const float& getTrackEtaRel() const { return trackEtaRel_; }
        const float& getTrackJetDecayLen() const { return trackJetDecayLen_; }
        const float& getTrackJetDistSig() const { return trackJetDistSig_; }
        const float& getTrackJetDistVal() const { return trackJetDistVal_; }
        const float& getTrackMomentum() const { return trackMomentum_; }
        const float& getTrackPPar() const { return trackPPar_; }
        const float& getTrackPParRatio() const { return trackPParRatio_; }
        const float& getTrackPtRatio() const { return trackPtRatio_; }
        const float& getTrackPtRel() const { return trackPtRel_; }
        const float& getTrackSip2dSig() const { return trackSip2dSig_; }
        const float& getTrackSip2dVal() const { return trackSip2dVal_; }
        const float& getTrackSip3dSig() const { return trackSip3dSig_; }
        const float& getTrackSip3dVal() const { return trackSip3dVal_; }
        const reco::TransientTrack getTTrack() const { return ttrack_; }

    private:

        edm::ESHandle<TransientTrackBuilder>& builder_;
        edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> track_builder_token_;

        float trackMomentum_;
        float trackEta_;
        float trackEtaRel_;
        float trackPtRel_;
        float trackPPar_;
        float trackDeltaR_;
        float trackPtRatio_;
        float trackPParRatio_;
        float trackSip2dVal_;
        float trackSip2dSig_;
        float trackSip3dVal_;
        float trackSip3dSig_;

        float trackJetDecayLen_;
        float trackJetDistVal_;
        float trackJetDistSig_;
        reco::TransientTrack ttrack_;
};


class GenVertex {

    public:

        GenVertex() : pt_(-1.0), x_(-1.0), y_(-1.0), z_(-1.0), daughters_(0) {}

        const float pt() const { return pt_; }
        const float x() const { return x_; }
        const float y() const { return y_; }
        const float z() const { return z_; }
        const float motherEta() const { return mother_->eta(); }
        const int motherPdgId() const { return mother_->pdgId(); }
        const unsigned int nDaughters() const { return nDaughters_; }
        const std::vector<const reco::Candidate*> daughters() const { return daughters_; }

        void setMother(const pat::PackedGenParticle* mother) { mother_ = mother; }
        void addDaughter(const reco::Candidate* daughter) { daughters_.push_back(daughter); }
        void setGenVertexAttributes() {
            if (!mother_) {
                std::cout << "WARNING in GenVertex::setGenVertexAttributes() in ntuple_SV.cc: add mother particle first!" << std::endl;
                return;
            }
            nDaughters_ = daughters_.size();
            if (nDaughters() <= 1) {
                std::cout << "WARNING in GenVertex::setGenVertexAttributes() in ntuple_SV.cc: add daughter particles first!" << std::endl;
                return;
            }
            pt_ = mother_->pt();
            x_ = daughters_.at(0)->vx();
            y_ = daughters_.at(0)->vy();
            z_ = daughters_.at(0)->vz();
        }

        void print() {
            std::cout << "GenVertex:" << std::endl;
            std::cout << "    vertex pt       = " << pt() << std::endl;
            std::cout << "    vertex x        = " << x() << std::endl;
            std::cout << "    vertex y        = " << y() << std::endl;
            std::cout << "    vertex z        = " << z() << std::endl;
            std::cout << "    mother pdg id   = " << motherPdgId() << std::endl;
            std::cout << "    mother pt       = " << mother_->pt() << std::endl;
            std::cout << "    mother eta      = " << motherEta() << std::endl;
            std::cout << "    mother phi      = " << mother_->phi() << std::endl;
            std::cout << "    nDaughters      = " << nDaughters() << std::endl;
            std::cout << "    daughter pdgIds = ";
            for (unsigned int iDau = 0; iDau < nDaughters(); iDau++) {
                std::cout << daughters().at(iDau)->pdgId();
                if (iDau < nDaughters() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
            std::cout << "    daughter etas   = ";
            for (unsigned int iDau = 0; iDau < nDaughters(); iDau++) {
                std::cout << daughters().at(iDau)->eta();
                if (iDau < nDaughters() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
            std::cout << "    daughter phis   = ";
            for (unsigned int iDau = 0; iDau < nDaughters(); iDau++) {
                std::cout << daughters().at(iDau)->phi();
                if (iDau < nDaughters() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
        }

    private:

        const pat::PackedGenParticle* mother_;
        float pt_; // the pt of the incoming particle
        float x_; // take vertex position of daughters as position of GenVertex
        float y_;
        float z_;
        unsigned int nDaughters_;
        std::vector<const reco::Candidate*> daughters_;
};


const reco::Vertex* ntuple_SV::spvp_;


ntuple_SV::ntuple_SV(std::string prefix, double jetR) : ntuple_content(jetR), sv_num_(0) {

    prefix_ = prefix;
}


ntuple_SV::~ntuple_SV() {}


void ntuple_SV::getInput(const edm::ParameterSet& iConfig) {

    timeQualityCut_ = iConfig.getParameter<double>("timeQualityCut");
    matchGVdR_ = iConfig.getParameter<double>("matchGVdR");
    jetPtMin_ = iConfig.getParameter<double>("jetPtMin");
    jetPtMax_ = iConfig.getParameter<double>("jetPtMax");
    jetAbsEtaMin_ = iConfig.getParameter<double>("jetAbsEtaMin");
    jetAbsEtaMax_ = iConfig.getParameter<double>("jetAbsEtaMax");
    genJetMatchdR_ = iConfig.getParameter<double>("genJetMatchdR");
}


void ntuple_SV::initBranches(TTree* tree) {

    // SV candidates
    addBranch(tree, (prefix_ + "n_sv").c_str(), &sv_num_, (prefix_ + "sv_num_/I").c_str());
    // addBranch(tree, (prefix_ + "nsv").c_str(), &nsv_, (prefix_ + "nsv_/F").c_str());
    // addBranch(tree, (prefix_ + "sv_pt").c_str(), &sv_pt_, (prefix_ + "sv_pt_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_eta").c_str(), &sv_eta_, (prefix_ + "sv_eta_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_phi").c_str(), &sv_phi_, (prefix_ + "sv_phi_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_e").c_str(), &sv_e_, (prefix_ + "sv_e_[" + prefix_ + "sv_num_]/F").c_str());
    addBranch(tree, (prefix_ + "sv_etarel").c_str(), &sv_etarel_);
    addBranch(tree, (prefix_ + "sv_phirel").c_str(), &sv_phirel_);
    addBranch(tree, (prefix_ + "sv_deltaR").c_str(), &sv_deltaR_);
    // addBranch(tree, (prefix_ + "sv_mass").c_str(), &sv_mass_, (prefix_ + "sv_mass_[" + prefix_ + "sv_num_]/F").c_str());
    addBranch(tree, (prefix_ + "sv_ntracks").c_str(), &sv_ntracks_, (prefix_ + "sv_ntracks_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_nMatchPFCand").c_str(), &sv_nMatchPFCand_, (prefix_ + "sv_nMatchPFCand_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_nUnmatchedTrk").c_str(), &sv_nUnmatchedTrk_, (prefix_ + "sv_nUnmatchedTrk_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_chi2").c_str(), &sv_chi2_, (prefix_ + "sv_chi2_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_ndf").c_str(), &sv_ndf_, (prefix_ + "sv_ndf_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_normchi2").c_str(), &sv_normchi2_, (prefix_ + "sv_normchi2_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_dxy").c_str(), &sv_dxy_, (prefix_ + "sv_dxy_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_dxyerr").c_str(), &sv_dxyerr_, (prefix_ + "sv_dxyerr_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_dxysig").c_str(), &sv_dxysig_, (prefix_ + "sv_dxysig_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_d3d").c_str(), &sv_d3d_, (prefix_ + "sv_d3d_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_d3derr").c_str(), &sv_d3derr_, (prefix_ + "sv_d3err_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_d3dsig").c_str(), &sv_d3dsig_, (prefix_ + "sv_d3dsig_[" + prefix_ + "sv_num_]/F").c_str());
    addBranch(tree, (prefix_ + "sv_costhetasvpv").c_str(), &sv_costhetasvpv_, (prefix_ + "sv_costhetasvpv_[" + prefix_ + "sv_num_]/F").c_str());
    addBranch(tree, (prefix_ + "sv_enratio").c_str(), &sv_enratio_);
    addBranch(tree, (prefix_ + "sv_calo_frac").c_str(), &sv_calo_frac_, (prefix_ + "sv_calo_frac_[" + prefix_ + "sv_num_]/F").c_str());
    addBranch(tree, (prefix_ + "sv_hcal_frac").c_str(), &sv_hcal_frac_, (prefix_ + "sv_hcal_frac_[" + prefix_ + "sv_num_]/F").c_str());
    addBranch(tree, (prefix_ + "sv_puppiw").c_str(), &sv_puppiw_, (prefix_ + "sv_puppiw_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_dz").c_str(), &sv_dz_, (prefix_ + "sv_dz_[" + prefix_ + "sv_num_]/F").c_str());
    addBranch(tree, (prefix_ + "sv_charge_sum").c_str(), &sv_charge_sum_, (prefix_ + "sv_charge_sum_[" + prefix_ + "sv_num_]/F").c_str());
    addBranch(tree, (prefix_ + "sv_pfd2dval").c_str(), &sv_pfd2dval_);
    addBranch(tree, (prefix_ + "sv_pfd2dsig").c_str(), &sv_pfd2dsig_);
    addBranch(tree, (prefix_ + "sv_pfd3dval").c_str(), &sv_pfd3dval_);
    addBranch(tree, (prefix_ + "sv_pfd3dsig").c_str(), &sv_pfd3dsig_);
    addBranch(tree, (prefix_ + "sv_time").c_str(), &sv_time_);

    addBranch(tree, (prefix_ + "gp_strange_status").c_str(), &gp_strange_status_);
    addBranch(tree, (prefix_ + "gp_strange_pdgId").c_str(), &gp_strange_pdgId_);
    addBranch(tree, (prefix_ + "gp_strange_nDaughters").c_str(), &gp_strange_nDaughters_);
    addBranch(tree, (prefix_ + "gp_strange_daughterPdgIds").c_str(), &gp_strange_daughterPdgIds_);
    addBranch(tree, (prefix_ + "gp_strange_pt").c_str(), &gp_strange_pt_);
    addBranch(tree, (prefix_ + "gp_strange_eta").c_str(), &gp_strange_eta_);

    addBranch(tree, (prefix_ + "st_pt").c_str(), &st_pt_);
    addBranch(tree, (prefix_ + "st_eta").c_str(), &st_eta_);
    addBranch(tree, (prefix_ + "st_phi").c_str(), &st_phi_);
    addBranch(tree, (prefix_ + "st_charge").c_str(), &st_charge_);

    addBranch(tree, (prefix_ + "rt_pt").c_str(), &rt_pt_);
    addBranch(tree, (prefix_ + "rt_eta").c_str(), &rt_eta_);
    addBranch(tree, (prefix_ + "rt_phi").c_str(), &rt_phi_);
    addBranch(tree, (prefix_ + "rt_charge").c_str(), &rt_charge_);
    addBranch(tree, (prefix_ + "rt_chi2").c_str(), &rt_chi2_);
    addBranch(tree, (prefix_ + "rt_ndof").c_str(), &rt_ndof_);
    addBranch(tree, (prefix_ + "rt_chi2dof").c_str(), &rt_chi2dof_);
    addBranch(tree, (prefix_ + "rt_dxy").c_str(), &rt_dxy_);
    addBranch(tree, (prefix_ + "rt_dxyerr").c_str(), &rt_dxyerr_);
    addBranch(tree, (prefix_ + "rt_d0").c_str(), &rt_d0_);
    addBranch(tree, (prefix_ + "rt_d0err").c_str(), &rt_d0err_);
    addBranch(tree, (prefix_ + "rt_dz").c_str(), &rt_dz_);
    addBranch(tree, (prefix_ + "rt_dzerr").c_str(), &rt_dzerr_);

    addBranch(tree, (prefix_ + "gv_pt").c_str(), &gv_pt_);
    addBranch(tree, (prefix_ + "gv_x").c_str(), &gv_x_);
    addBranch(tree, (prefix_ + "gv_y").c_str(), &gv_y_);
    addBranch(tree, (prefix_ + "gv_z").c_str(), &gv_z_);
    addBranch(tree, (prefix_ + "gv_motherEta").c_str(), &gv_motherEta_);
    addBranch(tree, (prefix_ + "gv_motherPdgId").c_str(), &gv_motherPdgId_);
    addBranch(tree, (prefix_ + "gv_nDaughters").c_str(), &gv_nDaughters_);
    addBranch(tree, (prefix_ + "gv_daughterPt").c_str(), &gv_daughterPt_);
    addBranch(tree, (prefix_ + "gv_daughterEta").c_str(), &gv_daughterEta_);
    addBranch(tree, (prefix_ + "gv_daughterPdgId").c_str(), &gv_daughterPdgId_);
    addBranch(tree, (prefix_ + "gv_dRtoClosestSV").c_str(), &gv_dRtoClosestSV_);
    addBranch(tree, (prefix_ + "gv_dRtoAllSV").c_str(), &gv_dRtoAllSV_);
    addBranch(tree, (prefix_ + "gv_nTimesMatchedToSV").c_str(), &gv_nTimesMatchedToSV_);
    addBranch(tree, (prefix_ + "gv_nTimesMatchedToSVinJet").c_str(), &gv_nTimesMatchedToSVinJet_);

    addBranch(tree, (prefix_ + "gv_matchSV_pt").c_str(), &gv_matchSV_pt_);
    addBranch(tree, (prefix_ + "gv_matchSV_x").c_str(), &gv_matchSV_x_);
    addBranch(tree, (prefix_ + "gv_matchSV_y").c_str(), &gv_matchSV_y_);
    addBranch(tree, (prefix_ + "gv_matchSV_z").c_str(), &gv_matchSV_z_);
    addBranch(tree, (prefix_ + "gv_matchSV_motherEta").c_str(), &gv_matchSV_motherEta_);
    addBranch(tree, (prefix_ + "gv_matchSV_motherPdgId").c_str(), &gv_matchSV_motherPdgId_);
    addBranch(tree, (prefix_ + "gv_matchSV_nDaughters").c_str(), &gv_matchSV_nDaughters_);
    addBranch(tree, (prefix_ + "gv_matchSV_daughterPt").c_str(), &gv_matchSV_daughterPt_);
    addBranch(tree, (prefix_ + "gv_matchSV_daughterEta").c_str(), &gv_matchSV_daughterEta_);
    addBranch(tree, (prefix_ + "gv_matchSV_daughterPdgId").c_str(), &gv_matchSV_daughterPdgId_);
    addBranch(tree, (prefix_ + "gv_matchSV_dRtoClosestSV").c_str(), &gv_matchSV_dRtoClosestSV_);

    addBranch(tree, (prefix_ + "sv_x").c_str(), &sv_x_);
    addBranch(tree, (prefix_ + "sv_y").c_str(), &sv_y_);
    addBranch(tree, (prefix_ + "sv_z").c_str(), &sv_z_);
    addBranch(tree, (prefix_ + "sv_pt").c_str(), &sv_pt_);
    addBranch(tree, (prefix_ + "sv_eta").c_str(), &sv_eta_);
    addBranch(tree, (prefix_ + "sv_phi").c_str(), &sv_phi_);
    addBranch(tree, (prefix_ + "sv_mass").c_str(), &sv_mass_);
    addBranch(tree, (prefix_ + "sv_energy").c_str(), &sv_energy_);
    addBranch(tree, (prefix_ + "sv_chi2").c_str(), &sv_chi2_);
    addBranch(tree, (prefix_ + "sv_ndof").c_str(), &sv_ndof_);
    addBranch(tree, (prefix_ + "sv_chi2dof").c_str(), &sv_chi2dof_);
    addBranch(tree, (prefix_ + "sv_dxy").c_str(), &sv_dxy_);
    addBranch(tree, (prefix_ + "sv_dxyerr").c_str(), &sv_dxyerr_);
    addBranch(tree, (prefix_ + "sv_dxysig").c_str(), &sv_dxysig_);
    addBranch(tree, (prefix_ + "sv_dz").c_str(), &sv_dz_);
    addBranch(tree, (prefix_ + "sv_d3D").c_str(), &sv_d3D_);
    addBranch(tree, (prefix_ + "sv_d3Derr").c_str(), &sv_d3Derr_);
    addBranch(tree, (prefix_ + "sv_d3Dsig").c_str(), &sv_d3Dsig_);
    addBranch(tree, (prefix_ + "sv_nDaughters").c_str(), &sv_nDaughters_);

    addBranch(tree, (prefix_ + "sv_matchGV_x").c_str(), &sv_matchGV_x_);
    addBranch(tree, (prefix_ + "sv_matchGV_y").c_str(), &sv_matchGV_y_);
    addBranch(tree, (prefix_ + "sv_matchGV_z").c_str(), &sv_matchGV_z_);
    addBranch(tree, (prefix_ + "sv_matchGV_pt").c_str(), &sv_matchGV_pt_);
    addBranch(tree, (prefix_ + "sv_matchGV_eta").c_str(), &sv_matchGV_eta_);
    addBranch(tree, (prefix_ + "sv_matchGV_phi").c_str(), &sv_matchGV_phi_);
    addBranch(tree, (prefix_ + "sv_matchGV_mass").c_str(), &sv_matchGV_mass_);
    addBranch(tree, (prefix_ + "sv_matchGV_energy").c_str(), &sv_matchGV_energy_);
    addBranch(tree, (prefix_ + "sv_matchGV_chi2").c_str(), &sv_matchGV_chi2_);
    addBranch(tree, (prefix_ + "sv_matchGV_ndof").c_str(), &sv_matchGV_ndof_);
    addBranch(tree, (prefix_ + "sv_matchGV_chi2dof").c_str(), &sv_matchGV_chi2dof_);
    addBranch(tree, (prefix_ + "sv_matchGV_dxy").c_str(), &sv_matchGV_dxy_);
    addBranch(tree, (prefix_ + "sv_matchGV_dxyerr").c_str(), &sv_matchGV_dxyerr_);
    addBranch(tree, (prefix_ + "sv_matchGV_dxysig").c_str(), &sv_matchGV_dxysig_);
    addBranch(tree, (prefix_ + "sv_matchGV_dz").c_str(), &sv_matchGV_dz_);
    addBranch(tree, (prefix_ + "sv_matchGV_d3D").c_str(), &sv_matchGV_d3D_);
    addBranch(tree, (prefix_ + "sv_matchGV_d3Derr").c_str(), &sv_matchGV_d3Derr_);
    addBranch(tree, (prefix_ + "sv_matchGV_d3Dsig").c_str(), &sv_matchGV_d3Dsig_);
    addBranch(tree, (prefix_ + "sv_matchGV_SVdRtoGV").c_str(), &sv_matchGV_SVdRtoGV_);
    addBranch(tree, (prefix_ + "sv_matchGV_nDaughters").c_str(), &sv_matchGV_nDaughters_);
    addBranch(tree, (prefix_ + "sv_matchGV_GVmotherPdgId").c_str(), &sv_matchGV_GVmotherPdgId_);

    addBranch(tree, (prefix_ + "sv_matchGV_matchJet_x").c_str(), &sv_matchGV_matchJet_x_);
    addBranch(tree, (prefix_ + "sv_matchGV_matchJet_y").c_str(), &sv_matchGV_matchJet_y_);
    addBranch(tree, (prefix_ + "sv_matchGV_matchJet_z").c_str(), &sv_matchGV_matchJet_z_);
    addBranch(tree, (prefix_ + "sv_matchGV_matchJet_pt").c_str(), &sv_matchGV_matchJet_pt_);
    addBranch(tree, (prefix_ + "sv_matchGV_matchJet_eta").c_str(), &sv_matchGV_matchJet_eta_);
    addBranch(tree, (prefix_ + "sv_matchGV_matchJet_phi").c_str(), &sv_matchGV_matchJet_phi_);
    addBranch(tree, (prefix_ + "sv_matchGV_matchJet_mass").c_str(), &sv_matchGV_matchJet_mass_);
    addBranch(tree, (prefix_ + "sv_matchGV_matchJet_energy").c_str(), &sv_matchGV_matchJet_energy_);
    addBranch(tree, (prefix_ + "sv_matchGV_matchJet_chi2").c_str(), &sv_matchGV_matchJet_chi2_);
    addBranch(tree, (prefix_ + "sv_matchGV_matchJet_ndof").c_str(), &sv_matchGV_matchJet_ndof_);
    addBranch(tree, (prefix_ + "sv_matchGV_matchJet_chi2dof").c_str(), &sv_matchGV_matchJet_chi2dof_);
    addBranch(tree, (prefix_ + "sv_matchGV_matchJet_dxy").c_str(), &sv_matchGV_matchJet_dxy_);
    addBranch(tree, (prefix_ + "sv_matchGV_matchJet_dxyerr").c_str(), &sv_matchGV_matchJet_dxyerr_);
    addBranch(tree, (prefix_ + "sv_matchGV_matchJet_dxysig").c_str(), &sv_matchGV_matchJet_dxysig_);
    addBranch(tree, (prefix_ + "sv_matchGV_matchJet_dz").c_str(), &sv_matchGV_matchJet_dz_);
    addBranch(tree, (prefix_ + "sv_matchGV_matchJet_d3D").c_str(), &sv_matchGV_matchJet_d3D_);
    addBranch(tree, (prefix_ + "sv_matchGV_matchJet_d3Derr").c_str(), &sv_matchGV_matchJet_d3Derr_);
    addBranch(tree, (prefix_ + "sv_matchGV_matchJet_d3Dsig").c_str(), &sv_matchGV_matchJet_d3Dsig_);
    addBranch(tree, (prefix_ + "sv_matchGV_matchJet_SVdRtoGV").c_str(), &sv_matchGV_matchJet_SVdRtoGV_);
    addBranch(tree, (prefix_ + "sv_matchGV_matchJet_nJets").c_str(), &sv_matchGV_matchJet_nJets_);
    addBranch(tree, (prefix_ + "sv_matchGV_matchJet_genJetFlav").c_str(), &sv_matchGV_matchJet_genJetHadFlav_);
    addBranch(tree, (prefix_ + "sv_matchGV_matchJet_nDaughters").c_str(), &sv_matchGV_matchJet_nDaughters_);
    addBranch(tree, (prefix_ + "sv_matchGV_matchJet_GVmotherPdgId").c_str(), &sv_matchGV_matchJet_GVmotherPdgId_);

    addBranch(tree, (prefix_ + "nPV").c_str(), &nPV_);
    addBranch(tree, (prefix_ + "nSV").c_str(), &nSV_);
    addBranch(tree, (prefix_ + "nSVMTDTiming").c_str(), &nSVMTDTiming_);

    addBranch(tree, (prefix_ + "pv_trk_timeVal").c_str(), &pv_trk_timeVal_);
    addBranch(tree, (prefix_ + "pv_trk_timeErr").c_str(), &pv_trk_timeErr_);
    addBranch(tree, (prefix_ + "pv_trk_timeQual").c_str(), &pv_trk_timeQual_);
    addBranch(tree, (prefix_ + "sv_trk_timeVal").c_str(), &sv_trk_timeVal_);
    addBranch(tree, (prefix_ + "sv_trk_timeErr").c_str(), &sv_trk_timeErr_);
    addBranch(tree, (prefix_ + "sv_trk_timeQual").c_str(), &sv_trk_timeQual_);

    addBranch(tree, (prefix_ + "jet_pt").c_str(), &jet_pt_);
    addBranch(tree, (prefix_ + "jet_eta").c_str(), &jet_eta_);
    addBranch(tree, (prefix_ + "jet_phi").c_str(), &jet_phi_);
    addBranch(tree, (prefix_ + "jet_radius").c_str(), &jet_radius_);
    addBranch(tree, (prefix_ + "jet_hadFlav").c_str(), &jet_hadFlav_);
    addBranch(tree, (prefix_ + "jet_partFlav").c_str(), &jet_partFlav_);
    addBranch(tree, (prefix_ + "jet_genHadFlav").c_str(), &jet_genHadFlav_);
    addBranch(tree, (prefix_ + "jet_genPartFlav").c_str(), &jet_genPartFlav_);

    addBranch(tree, (prefix_ + "jet_matchSV_pt").c_str(), &jet_matchSV_pt_);
    addBranch(tree, (prefix_ + "jet_matchSV_eta").c_str(), &jet_matchSV_eta_);
    addBranch(tree, (prefix_ + "jet_matchSV_phi").c_str(), &jet_matchSV_phi_);
    addBranch(tree, (prefix_ + "jet_matchSV_radius").c_str(), &jet_matchSV_radius_);
    addBranch(tree, (prefix_ + "jet_matchSV_nSV").c_str(), &jet_matchSV_nSV_);
    addBranch(tree, (prefix_ + "jet_matchSV_hadFlav").c_str(), &jet_matchSV_hadFlav_);
    addBranch(tree, (prefix_ + "jet_matchSV_partFlav").c_str(), &jet_matchSV_partFlav_);
    addBranch(tree, (prefix_ + "jet_matchSV_genHadFlav").c_str(), &jet_matchSV_genHadFlav_);
    addBranch(tree, (prefix_ + "jet_matchSV_genPartFlav").c_str(), &jet_matchSV_genPartFlav_);

    addBranch(tree, (prefix_ + "jet_matchSV_matchGV_pt").c_str(), &jet_matchSV_matchGV_pt_);
    addBranch(tree, (prefix_ + "jet_matchSV_matchGV_eta").c_str(), &jet_matchSV_matchGV_eta_);
    addBranch(tree, (prefix_ + "jet_matchSV_matchGV_phi").c_str(), &jet_matchSV_matchGV_phi_);
    addBranch(tree, (prefix_ + "jet_matchSV_matchGV_radius").c_str(), &jet_matchSV_matchGV_radius_);
    addBranch(tree, (prefix_ + "jet_matchSV_matchGV_SVdRtoGV").c_str(), &jet_matchSV_matchGV_SVdRtoGV_);
    addBranch(tree, (prefix_ + "jet_matchSV_matchGV_nGV").c_str(), &jet_matchSV_matchGV_nGV_);
    addBranch(tree, (prefix_ + "jet_matchSV_matchGV_hadFlav").c_str(), &jet_matchSV_matchGV_hadFlav_);
    addBranch(tree, (prefix_ + "jet_matchSV_matchGV_partFlav").c_str(), &jet_matchSV_matchGV_partFlav_);
    addBranch(tree, (prefix_ + "jet_matchSV_matchGV_genHadFlav").c_str(), &jet_matchSV_matchGV_genHadFlav_);
    addBranch(tree, (prefix_ + "jet_matchSV_matchGV_genPartFlav").c_str(), &jet_matchSV_matchGV_genPartFlav_);
}


void ntuple_SV::readSetup(const edm::EventSetup& iSetup) {

    builder_ = iSetup.getHandle(track_builder_token_);
}


void ntuple_SV::readEvent(const edm::Event& iEvent) {

    iEvent.getByToken(genParticles_token_, genParticles_);
    iEvent.getByToken(genJetFlavourInfo_token_, genJetFlavourInfo_);
    iEvent.getByToken(simTracks_token_, simTracks_);
    // iEvent.getByToken(pf_cand_token_, pf_cand_);
    // iEvent.getByToken(pf_mcmatch_token_, pf_mcmatch_);
    iEvent.getByToken(recoTracks_token_, recoTracks_);
    iEvent.getByToken(timeValueMap_token_, timeValueMap_);
    iEvent.getByToken(timeErrorMap_token_, timeErrorMap_);
    iEvent.getByToken(timeQualityMap_token_, timeQualityMap_);
    // iEvent.getByToken(genVertices_token_, genVertices_);
    iEvent.getByToken(PVs_token_, PVs_);
    iEvent.getByToken(inclusiveSVs_token_, inclusiveSVs_);
    iEvent.getByToken(inclusiveSVsMTDTiming_token_, inclusiveSVsMTDTiming_);
}


void ntuple_SV::initContainers() {

    sv_etarel_ = new std::vector<float>;
    sv_phirel_ = new std::vector<float>;
    sv_deltaR_ = new std::vector<float>;
    sv_enratio_ = new std::vector<float>;
    sv_pfd2dval_ = new std::vector<float>;
    sv_pfd2dsig_ = new std::vector<float>;
    sv_pfd3dval_ = new std::vector<float>;
    sv_pfd3dsig_ = new std::vector<float>;

    sv_time_ = new std::vector<float>;

    gp_strange_status_ = new std::vector<int>;
    gp_strange_pdgId_ = new std::vector<int>;
    gp_strange_nDaughters_ = new std::vector<int>;
    gp_strange_daughterPdgIds_ = new std::vector<int>;
    gp_strange_pt_ = new std::vector<float>;
    gp_strange_eta_ = new std::vector<float>;

    st_pt_ = new std::vector<float>;
    st_eta_ = new std::vector<float>;
    st_phi_ = new std::vector<float>;
    st_charge_ = new std::vector<int>;

    rt_pt_ = new std::vector<float>;
    rt_eta_ = new std::vector<float>;
    rt_phi_ = new std::vector<float>;
    rt_charge_ = new std::vector<int>;
    rt_chi2_ = new std::vector<float>;
    rt_ndof_ = new std::vector<float>;
    rt_chi2dof_ = new std::vector<float>;
    rt_dxy_ = new std::vector<float>;
    rt_dxyerr_ = new std::vector<float>;
    rt_d0_ = new std::vector<float>;
    rt_d0err_ = new std::vector<float>;
    rt_dz_ = new std::vector<float>;
    rt_dzerr_ = new std::vector<float>;

    gv_pt_ = new std::vector<float>;
    gv_x_ = new std::vector<float>;
    gv_y_ = new std::vector<float>;
    gv_z_ = new std::vector<float>;
    gv_motherEta_ = new std::vector<float>;
    gv_motherPdgId_ = new std::vector<int>;
    gv_nDaughters_ = new std::vector<int>;
    gv_daughterPt_ = new std::vector<float>;
    gv_daughterEta_ = new std::vector<float>;
    gv_daughterPdgId_ = new std::vector<int>;
    gv_dRtoClosestSV_ = new std::vector<float>;
    gv_dRtoAllSV_ = new std::vector<float>;
    gv_nTimesMatchedToSV_ = new std::vector<int>;
    gv_nTimesMatchedToSVinJet_ = new std::vector<int>;

    gv_matchSV_pt_ = new std::vector<float>;
    gv_matchSV_x_ = new std::vector<float>;
    gv_matchSV_y_ = new std::vector<float>;
    gv_matchSV_z_ = new std::vector<float>;
    gv_matchSV_motherEta_ = new std::vector<float>;
    gv_matchSV_motherPdgId_ = new std::vector<int>;
    gv_matchSV_nDaughters_ = new std::vector<int>;
    gv_matchSV_daughterPt_ = new std::vector<float>;
    gv_matchSV_daughterEta_ = new std::vector<float>;
    gv_matchSV_daughterPdgId_ = new std::vector<int>;
    gv_matchSV_dRtoClosestSV_ = new std::vector<float>;

    sv_x_ = new std::vector<float>;
    sv_y_ = new std::vector<float>;
    sv_z_ = new std::vector<float>;
    sv_pt_ = new std::vector<float>;
    sv_eta_ = new std::vector<float>;
    sv_phi_ = new std::vector<float>;
    sv_mass_ = new std::vector<float>;
    sv_energy_ = new std::vector<float>;
    sv_chi2_ = new std::vector<float>;
    sv_ndof_ = new std::vector<float>;
    sv_chi2dof_ = new std::vector<float>;
    sv_dxy_ = new std::vector<float>;
    sv_dxyerr_ = new std::vector<float>;
    sv_dxysig_ = new std::vector<float>;
    sv_dz_ = new std::vector<float>;
    sv_d3D_ = new std::vector<float>;
    sv_d3Derr_ = new std::vector<float>;
    sv_d3Dsig_ = new std::vector<float>;
    sv_nDaughters_ = new std::vector<int>;

    sv_matchGV_x_ = new std::vector<float>;
    sv_matchGV_y_ = new std::vector<float>;
    sv_matchGV_z_ = new std::vector<float>;
    sv_matchGV_pt_ = new std::vector<float>;
    sv_matchGV_eta_ = new std::vector<float>;
    sv_matchGV_phi_ = new std::vector<float>;
    sv_matchGV_mass_ = new std::vector<float>;
    sv_matchGV_energy_ = new std::vector<float>;
    sv_matchGV_chi2_ = new std::vector<float>;
    sv_matchGV_ndof_ = new std::vector<float>;
    sv_matchGV_chi2dof_ = new std::vector<float>;
    sv_matchGV_dxy_ = new std::vector<float>;
    sv_matchGV_dxyerr_ = new std::vector<float>;
    sv_matchGV_dxysig_ = new std::vector<float>;
    sv_matchGV_dz_ = new std::vector<float>;
    sv_matchGV_d3D_ = new std::vector<float>;
    sv_matchGV_d3Derr_ = new std::vector<float>;
    sv_matchGV_d3Dsig_ = new std::vector<float>;
    sv_matchGV_SVdRtoGV_ = new std::vector<float>;
    sv_matchGV_nDaughters_ = new std::vector<int>;
    sv_matchGV_GVmotherPdgId_ = new std::vector<int>;

    sv_matchGV_matchJet_x_ = new std::vector<float>;
    sv_matchGV_matchJet_y_ = new std::vector<float>;
    sv_matchGV_matchJet_z_ = new std::vector<float>;
    sv_matchGV_matchJet_pt_ = new std::vector<float>;
    sv_matchGV_matchJet_eta_ = new std::vector<float>;
    sv_matchGV_matchJet_phi_ = new std::vector<float>;
    sv_matchGV_matchJet_mass_ = new std::vector<float>;
    sv_matchGV_matchJet_energy_ = new std::vector<float>;
    sv_matchGV_matchJet_chi2_ = new std::vector<float>;
    sv_matchGV_matchJet_ndof_ = new std::vector<float>;
    sv_matchGV_matchJet_chi2dof_ = new std::vector<float>;
    sv_matchGV_matchJet_dxy_ = new std::vector<float>;
    sv_matchGV_matchJet_dxyerr_ = new std::vector<float>;
    sv_matchGV_matchJet_dxysig_ = new std::vector<float>;
    sv_matchGV_matchJet_dz_ = new std::vector<float>;
    sv_matchGV_matchJet_d3D_ = new std::vector<float>;
    sv_matchGV_matchJet_d3Derr_ = new std::vector<float>;
    sv_matchGV_matchJet_d3Dsig_ = new std::vector<float>;
    sv_matchGV_matchJet_SVdRtoGV_ = new std::vector<float>;
    sv_matchGV_matchJet_nJets_ = new std::vector<int>;
    sv_matchGV_matchJet_genJetHadFlav_ = new std::vector<int>;
    sv_matchGV_matchJet_nDaughters_ = new std::vector<int>;
    sv_matchGV_matchJet_GVmotherPdgId_ = new std::vector<int>;

    pv_trk_timeVal_ = new std::vector<float>;
    pv_trk_timeErr_ = new std::vector<float>;
    pv_trk_timeQual_ = new std::vector<float>;
    sv_trk_timeVal_ = new std::vector<float>;
    sv_trk_timeErr_ = new std::vector<float>;
    sv_trk_timeQual_ = new std::vector<float>;

    jet_pt_ = new std::vector<float>;
    jet_eta_ = new std::vector<float>;
    jet_phi_ = new std::vector<float>;
    jet_radius_ = new std::vector<float>;
    jet_hadFlav_ = new std::vector<int>;
    jet_partFlav_ = new std::vector<int>;
    jet_genHadFlav_ = new std::vector<int>;
    jet_genPartFlav_ = new std::vector<int>;

    jet_matchSV_pt_ = new std::vector<float>;
    jet_matchSV_eta_ = new std::vector<float>;
    jet_matchSV_phi_ = new std::vector<float>;
    jet_matchSV_radius_ = new std::vector<float>;
    jet_matchSV_nSV_ = new std::vector<int>;
    jet_matchSV_hadFlav_ = new std::vector<int>;
    jet_matchSV_partFlav_ = new std::vector<int>;
    jet_matchSV_genHadFlav_ = new std::vector<int>;
    jet_matchSV_genPartFlav_ = new std::vector<int>;

    jet_matchSV_matchGV_pt_ = new std::vector<float>;
    jet_matchSV_matchGV_eta_ = new std::vector<float>;
    jet_matchSV_matchGV_phi_ = new std::vector<float>;
    jet_matchSV_matchGV_radius_ = new std::vector<float>;
    jet_matchSV_matchGV_SVdRtoGV_ = new std::vector<float>;
    jet_matchSV_matchGV_nGV_ = new std::vector<int>;
    jet_matchSV_matchGV_hadFlav_ = new std::vector<int>;
    jet_matchSV_matchGV_partFlav_ = new std::vector<int>;
    jet_matchSV_matchGV_genHadFlav_ = new std::vector<int>;
    jet_matchSV_matchGV_genPartFlav_ = new std::vector<int>;
}


void ntuple_SV::clearContainers() {

    sv_etarel_->clear();
    sv_phirel_->clear();
    sv_deltaR_->clear();
    sv_enratio_->clear();
    sv_pfd2dval_->clear();
    sv_pfd2dsig_->clear();
    sv_pfd3dval_->clear();
    sv_pfd3dsig_->clear();

    sv_time_->clear();

    gp_strange_status_->clear();
    gp_strange_pdgId_->clear();
    gp_strange_nDaughters_->clear();
    gp_strange_daughterPdgIds_->clear();
    gp_strange_pt_->clear();
    gp_strange_eta_->clear();

    st_pt_->clear();
    st_eta_->clear();
    st_phi_->clear();
    st_charge_->clear();

    rt_pt_->clear();
    rt_eta_->clear();
    rt_phi_->clear();
    rt_charge_->clear();
    rt_chi2_->clear();
    rt_ndof_->clear();
    rt_chi2dof_->clear();
    rt_dxy_->clear();
    rt_dxyerr_->clear();
    rt_d0_->clear();
    rt_d0err_->clear();
    rt_dz_->clear();
    rt_dzerr_->clear();

    gv_pt_->clear();
    gv_x_->clear();
    gv_y_->clear();
    gv_z_->clear();
    gv_motherEta_->clear();
    gv_motherPdgId_->clear();
    gv_nDaughters_->clear();
    gv_daughterPt_->clear();
    gv_daughterEta_->clear();
    gv_daughterPdgId_->clear();
    gv_dRtoClosestSV_->clear();
    gv_dRtoAllSV_->clear();
    gv_nTimesMatchedToSV_->clear();
    gv_nTimesMatchedToSVinJet_->clear();

    gv_matchSV_pt_->clear();
    gv_matchSV_x_->clear();
    gv_matchSV_y_->clear();
    gv_matchSV_z_->clear();
    gv_matchSV_motherEta_->clear();
    gv_matchSV_motherPdgId_->clear();
    gv_matchSV_nDaughters_->clear();
    gv_matchSV_daughterPt_->clear();
    gv_matchSV_daughterEta_->clear();
    gv_matchSV_daughterPdgId_->clear();
    gv_matchSV_dRtoClosestSV_->clear();

    sv_x_->clear();
    sv_y_->clear();
    sv_z_->clear();
    sv_pt_->clear();
    sv_eta_->clear();
    sv_phi_->clear();
    sv_mass_->clear();
    sv_energy_->clear();
    sv_chi2_->clear();
    sv_ndof_->clear();
    sv_chi2dof_->clear();
    sv_dxy_->clear();
    sv_dxyerr_->clear();
    sv_dxysig_->clear();
    sv_dz_->clear();
    sv_d3D_->clear();
    sv_d3Derr_->clear();
    sv_d3Dsig_->clear();
    sv_nDaughters_->clear();

    sv_matchGV_x_->clear();
    sv_matchGV_y_->clear();
    sv_matchGV_z_->clear();
    sv_matchGV_pt_->clear();
    sv_matchGV_eta_->clear();
    sv_matchGV_phi_->clear();
    sv_matchGV_mass_->clear();
    sv_matchGV_energy_->clear();
    sv_matchGV_chi2_->clear();
    sv_matchGV_ndof_->clear();
    sv_matchGV_chi2dof_->clear();
    sv_matchGV_dxy_->clear();
    sv_matchGV_dxyerr_->clear();
    sv_matchGV_dxysig_->clear();
    sv_matchGV_dz_->clear();
    sv_matchGV_d3D_->clear();
    sv_matchGV_d3Derr_->clear();
    sv_matchGV_d3Dsig_->clear();
    sv_matchGV_SVdRtoGV_->clear();
    sv_matchGV_nDaughters_->clear();
    sv_matchGV_GVmotherPdgId_->clear();

    sv_matchGV_matchJet_x_->clear();
    sv_matchGV_matchJet_y_->clear();
    sv_matchGV_matchJet_z_->clear();
    sv_matchGV_matchJet_pt_->clear();
    sv_matchGV_matchJet_eta_->clear();
    sv_matchGV_matchJet_phi_->clear();
    sv_matchGV_matchJet_mass_->clear();
    sv_matchGV_matchJet_energy_->clear();
    sv_matchGV_matchJet_chi2_->clear();
    sv_matchGV_matchJet_ndof_->clear();
    sv_matchGV_matchJet_chi2dof_->clear();
    sv_matchGV_matchJet_dxy_->clear();
    sv_matchGV_matchJet_dxyerr_->clear();
    sv_matchGV_matchJet_dxysig_->clear();
    sv_matchGV_matchJet_dz_->clear();
    sv_matchGV_matchJet_d3D_->clear();
    sv_matchGV_matchJet_d3Derr_->clear();
    sv_matchGV_matchJet_d3Dsig_->clear();
    sv_matchGV_matchJet_SVdRtoGV_->clear();
    sv_matchGV_matchJet_nJets_->clear();
    sv_matchGV_matchJet_genJetHadFlav_->clear();
    sv_matchGV_matchJet_nDaughters_->clear();
    sv_matchGV_matchJet_GVmotherPdgId_->clear();

    pv_trk_timeVal_->clear();
    pv_trk_timeErr_->clear();
    pv_trk_timeQual_->clear();
    sv_trk_timeVal_->clear();
    sv_trk_timeErr_->clear();
    sv_trk_timeQual_->clear();

    jet_pt_->clear();
    jet_eta_->clear();
    jet_phi_->clear();
    jet_radius_->clear();
    jet_hadFlav_->clear();
    jet_partFlav_->clear();
    jet_genHadFlav_->clear();
    jet_genPartFlav_->clear();

    jet_matchSV_pt_->clear();
    jet_matchSV_eta_->clear();
    jet_matchSV_phi_->clear();
    jet_matchSV_radius_->clear();
    jet_matchSV_nSV_->clear();
    jet_matchSV_hadFlav_->clear();
    jet_matchSV_partFlav_->clear();
    jet_matchSV_genHadFlav_->clear();
    jet_matchSV_genPartFlav_->clear();

    jet_matchSV_matchGV_pt_->clear();
    jet_matchSV_matchGV_eta_->clear();
    jet_matchSV_matchGV_phi_->clear();
    jet_matchSV_matchGV_radius_->clear();
    jet_matchSV_matchGV_SVdRtoGV_->clear();
    jet_matchSV_matchGV_nGV_->clear();
    jet_matchSV_matchGV_hadFlav_->clear();
    jet_matchSV_matchGV_partFlav_->clear();
    jet_matchSV_matchGV_genHadFlav_->clear();
    jet_matchSV_matchGV_genPartFlav_->clear();
}


void ntuple_SV::deleteContainers() {

    delete sv_etarel_;
    delete sv_phirel_;
    delete sv_deltaR_;
    delete sv_enratio_;
    delete sv_pfd2dval_;
    delete sv_pfd2dsig_;
    delete sv_pfd3dval_;
    delete sv_pfd3dsig_;

    delete sv_time_;

    delete gp_strange_status_;
    delete gp_strange_pdgId_;
    delete gp_strange_nDaughters_;
    delete gp_strange_daughterPdgIds_;
    delete gp_strange_pt_;
    delete gp_strange_eta_;

    delete st_pt_;
    delete st_eta_;
    delete st_phi_;
    delete st_charge_;

    delete rt_pt_;
    delete rt_eta_;
    delete rt_phi_;
    delete rt_charge_;
    delete rt_chi2_;
    delete rt_ndof_;
    delete rt_chi2dof_;
    delete rt_dxy_;
    delete rt_dxyerr_;
    delete rt_d0_;
    delete rt_d0err_;
    delete rt_dz_;
    delete rt_dzerr_;

    delete gv_pt_;
    delete gv_x_;
    delete gv_y_;
    delete gv_z_;
    delete gv_motherEta_;
    delete gv_motherPdgId_;
    delete gv_nDaughters_;
    delete gv_daughterPt_;
    delete gv_daughterEta_;
    delete gv_daughterPdgId_;
    delete gv_dRtoClosestSV_;
    delete gv_dRtoAllSV_;
    delete gv_nTimesMatchedToSV_;
    delete gv_nTimesMatchedToSVinJet_;

    delete gv_matchSV_pt_;
    delete gv_matchSV_x_;
    delete gv_matchSV_y_;
    delete gv_matchSV_z_;
    delete gv_matchSV_motherEta_;
    delete gv_matchSV_motherPdgId_;
    delete gv_matchSV_nDaughters_;
    delete gv_matchSV_daughterPt_;
    delete gv_matchSV_daughterEta_;
    delete gv_matchSV_daughterPdgId_;
    delete gv_matchSV_dRtoClosestSV_;

    delete sv_x_;
    delete sv_y_;
    delete sv_z_;
    delete sv_pt_;
    delete sv_eta_;
    delete sv_phi_;
    delete sv_mass_;
    delete sv_energy_;
    delete sv_chi2_;
    delete sv_ndof_;
    delete sv_chi2dof_;
    delete sv_dxy_;
    delete sv_dxyerr_;
    delete sv_dxysig_;
    delete sv_dz_;
    delete sv_d3D_;
    delete sv_d3Derr_;
    delete sv_d3Dsig_;
    delete sv_nDaughters_;

    delete sv_matchGV_x_;
    delete sv_matchGV_y_;
    delete sv_matchGV_z_;
    delete sv_matchGV_pt_;
    delete sv_matchGV_eta_;
    delete sv_matchGV_phi_;
    delete sv_matchGV_mass_;
    delete sv_matchGV_energy_;
    delete sv_matchGV_chi2_;
    delete sv_matchGV_ndof_;
    delete sv_matchGV_chi2dof_;
    delete sv_matchGV_dxy_;
    delete sv_matchGV_dxyerr_;
    delete sv_matchGV_dxysig_;
    delete sv_matchGV_dz_;
    delete sv_matchGV_d3D_;
    delete sv_matchGV_d3Derr_;
    delete sv_matchGV_d3Dsig_;
    delete sv_matchGV_SVdRtoGV_;
    delete sv_matchGV_nDaughters_;
    delete sv_matchGV_GVmotherPdgId_;

    delete sv_matchGV_matchJet_x_;
    delete sv_matchGV_matchJet_y_;
    delete sv_matchGV_matchJet_z_;
    delete sv_matchGV_matchJet_pt_;
    delete sv_matchGV_matchJet_eta_;
    delete sv_matchGV_matchJet_phi_;
    delete sv_matchGV_matchJet_mass_;
    delete sv_matchGV_matchJet_energy_;
    delete sv_matchGV_matchJet_chi2_;
    delete sv_matchGV_matchJet_ndof_;
    delete sv_matchGV_matchJet_chi2dof_;
    delete sv_matchGV_matchJet_dxy_;
    delete sv_matchGV_matchJet_dxyerr_;
    delete sv_matchGV_matchJet_dxysig_;
    delete sv_matchGV_matchJet_dz_;
    delete sv_matchGV_matchJet_d3D_;
    delete sv_matchGV_matchJet_d3Derr_;
    delete sv_matchGV_matchJet_d3Dsig_;
    delete sv_matchGV_matchJet_SVdRtoGV_;
    delete sv_matchGV_matchJet_nJets_;
    delete sv_matchGV_matchJet_genJetHadFlav_;
    delete sv_matchGV_matchJet_nDaughters_;
    delete sv_matchGV_matchJet_GVmotherPdgId_;

    delete pv_trk_timeVal_;
    delete pv_trk_timeErr_;
    delete pv_trk_timeQual_;
    delete sv_trk_timeVal_;
    delete sv_trk_timeErr_;
    delete sv_trk_timeQual_;

    delete jet_pt_;
    delete jet_eta_;
    delete jet_phi_;
    delete jet_radius_;
    delete jet_hadFlav_;
    delete jet_partFlav_;
    delete jet_genHadFlav_;
    delete jet_genPartFlav_;

    delete jet_matchSV_pt_;
    delete jet_matchSV_eta_;
    delete jet_matchSV_phi_;
    delete jet_matchSV_radius_;
    delete jet_matchSV_nSV_;
    delete jet_matchSV_hadFlav_;
    delete jet_matchSV_partFlav_;
    delete jet_matchSV_genHadFlav_;
    delete jet_matchSV_genPartFlav_;

    delete jet_matchSV_matchGV_pt_;
    delete jet_matchSV_matchGV_eta_;
    delete jet_matchSV_matchGV_phi_;
    delete jet_matchSV_matchGV_radius_;
    delete jet_matchSV_matchGV_SVdRtoGV_;
    delete jet_matchSV_matchGV_nGV_;
    delete jet_matchSV_matchGV_hadFlav_;
    delete jet_matchSV_matchGV_partFlav_;
    delete jet_matchSV_matchGV_genHadFlav_;
    delete jet_matchSV_matchGV_genPartFlav_;
}


int ntuple_SV::fillBranches(bool applySelection, float EventTime) {

    clearContainers();

    // Sorting described here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideOfflinePrimaryVertexProduction
    const reco::Vertex& pv = vertices()->at(0); // Most likely the signal vertex
    spvp_ = &vertices()->at(0);

    reco::VertexCompositePtrCandidateCollection cpvtx = *secVertices();
    std::sort(cpvtx.begin(), cpvtx.end(), ntuple_SV::compareDxyDxyErr);

    const edm::View<pat::Jet> jetCollection = *jets();

    edm::SimTrackContainer simTracks = *(simTracks_.product());
    reco::TrackCollection recoTracks = *(recoTracks_.product());
    const edm::ValueMap<float>& timeValueMap = *(timeValueMap_.product());
    const edm::ValueMap<float>& timeErrorMap = *(timeErrorMap_.product());
    const edm::ValueMap<float>& timeQualityMap = *(timeQualityMap_.product());
    reco::VertexCollection PVs = *(PVs_.product());
    reco::VertexCollection inclusiveSVs = *(inclusiveSVs_.product());
    reco::VertexCollection inclusiveSVsMTDTimings = *(inclusiveSVsMTDTiming_.product());

    nPV_ = PVs.size();
    nSV_ = inclusiveSVs.size();
    nSVMTDTiming_ = inclusiveSVsMTDTimings.size();

    // Save SimTrack and Track information
    for (unsigned int iST = 0; iST < simTracks.size(); iST++) {
        SimTrack st = simTracks.at(iST);
        st_pt_->push_back(st.momentum().Pt());
        st_eta_->push_back(st.momentum().Eta());
        st_phi_->push_back(st.momentum().Phi());
        st_charge_->push_back(st.charge());
    }
    for (unsigned int iRT = 0; iRT < recoTracks.size(); iRT++) {
        reco::Track rt = recoTracks.at(iRT);
        rt_pt_->push_back(rt.pt());
        rt_eta_->push_back(rt.eta());
        rt_phi_->push_back(rt.phi());
        rt_charge_->push_back(rt.charge());
        rt_chi2_->push_back(rt.chi2());
        rt_ndof_->push_back(rt.ndof());
        rt_chi2dof_->push_back(rt.normalizedChi2());
        rt_dxy_->push_back(rt.dxy());
        rt_dxyerr_->push_back(rt.dxyError());
        rt_d0_->push_back(rt.d0());
        rt_d0err_->push_back(rt.d0Error());
        rt_dz_->push_back(rt.dz());
        rt_dzerr_->push_back(rt.dzError());
    }

    // Construct gen vertices by looping through list of gen particles
    std::vector<GenVertex*> genVertices(0);
    for (unsigned int iGenPart = 0; iGenPart < genParticles_->size(); iGenPart++) {
        const pat::PackedGenParticle gp = genParticles_->at(iGenPart);
        int motherPartID = genPartID(gp.pdgId());

        // APPLY THESE CUTS IN PYTHON SCRIPT
        // if (gp.pt() < 10.0) continue;
        // if (abs(gp.eta()) > 2.5) continue;
        if (motherPartID < 0) continue; // Mother is not interesting hadron

        if (motherPartID == 4) {
            GenVertex* newGenVtx = new GenVertex();
            newGenVtx->setMother(gp.clone());
            // std::cout << "NUMBER OF S BARYON DAUGHTERS = " << gp.numberOfDaughters() << std::endl;
            gp_strange_status_->push_back((int) (gp.statusFlags().flags_).to_ulong());
            gp_strange_pdgId_->push_back(gp.pdgId());
            gp_strange_nDaughters_->push_back(gp.numberOfDaughters());
            gp_strange_pt_->push_back(gp.pt());
            gp_strange_eta_->push_back(gp.eta());
            for (unsigned int iDau = 0; iDau < gp.numberOfDaughters(); iDau++) {
                const reco::Candidate* dau = gp.daughter(iDau)->clone();
                gp_strange_daughterPdgIds_->push_back(dau->pdgId());
                newGenVtx->addDaughter(dau);
            }
            newGenVtx->setGenVertexAttributes();
            // if (gp.numberOfDaughters() > 0)
            //     newGenVtx->print();
        }

        // Check for last instance of interesting hadron
        bool lastInstance = true;
        unsigned int nGoodDaughters = 0;
        for (unsigned int iDau = 0; iDau < gp.numberOfDaughters(); iDau++) {
            const reco::Candidate* dau = gp.daughter(iDau);
            int daughterPartID = genPartID(dau->pdgId());

            // if (dau->pt() < 10.0) continue;
            // if (abs(dau->eta()) > 2.5) continue;

            nGoodDaughters++;
            if (daughterPartID == motherPartID) {
                lastInstance = false; // Not last instance of interesting hadron
                break;
            }
        }

        if (lastInstance) {

            if (nGoodDaughters < 2) continue; // Make sure a vertex exists

            // Make new GenVertex
            GenVertex* newGenVtx = new GenVertex();
            newGenVtx->setMother(gp.clone());
            for (unsigned int iDau = 0; iDau < gp.numberOfDaughters(); iDau++) {
                const reco::Candidate* dau = gp.daughter(iDau)->clone();

                // Only add "good" daughters
                // if (dau->pt() < 10.0) continue;
                // if (abs(dau->eta()) > 2.5) continue;

                newGenVtx->addDaughter(dau);
            }
            newGenVtx->setGenVertexAttributes();
            genVertices.push_back(newGenVtx);
        }
    }

    // GenVertex stuff
    for (unsigned int iGV = 0; iGV < genVertices.size(); iGV++) {
        const GenVertex* genVtx = genVertices.at(iGV);
        gv_pt_->push_back(genVtx->pt());
        gv_x_->push_back(genVtx->x());
        gv_y_->push_back(genVtx->y());
        gv_z_->push_back(genVtx->z());
        gv_motherEta_->push_back(genVtx->motherEta());
        gv_motherPdgId_->push_back(genVtx->motherPdgId());
        gv_nDaughters_->push_back(genVtx->nDaughters());
        for (const reco::Candidate* dau : genVtx->daughters()) {
            gv_daughterPt_->push_back(dau->pt());
            gv_daughterEta_->push_back(dau->eta());
            gv_daughterPdgId_->push_back(dau->pdgId());
        }

        // Match GV to SV
        int nMatch = 0;
        int closestSVIdx = -1;
        float mindR = matchGVdR_; // some maximum value
        for (unsigned int iSV = 0; iSV < cpvtx.size(); iSV++) {
            if (iSV >= max_sv_) break;
            const reco::VertexCompositePtrCandidate& sv = cpvtx.at(iSV);
            float dR3D = deltaR3D(genVtx->x(), sv.vertex().x(), genVtx->y(), sv.vertex().y(), genVtx->z(), sv.vertex().z());
            if (dR3D < matchGVdR_) {
                gv_dRtoAllSV_->push_back(mindR);
                nMatch++;
                if (dR3D < mindR) {
                    closestSVIdx = iSV;
                    mindR = dR3D;
                }
            }
        }
        gv_dRtoClosestSV_->push_back(mindR);
        gv_nTimesMatchedToSV_->push_back(nMatch);
        if (closestSVIdx >= 0) {
            gv_matchSV_pt_->push_back(genVtx->pt());
            gv_matchSV_x_->push_back(genVtx->x());
            gv_matchSV_y_->push_back(genVtx->y());
            gv_matchSV_z_->push_back(genVtx->z());
            gv_matchSV_motherEta_->push_back(genVtx->motherEta());
            gv_matchSV_motherPdgId_->push_back(genVtx->motherPdgId());
            gv_matchSV_nDaughters_->push_back(genVtx->nDaughters());
            for (const reco::Candidate* dau : genVtx->daughters()) {
                gv_matchSV_daughterPt_->push_back(dau->pt());
                gv_matchSV_daughterEta_->push_back(dau->eta());
                gv_matchSV_daughterPdgId_->push_back(dau->pdgId());
            }
            gv_matchSV_dRtoClosestSV_->push_back(mindR);
        }
    }

    // Secondary Vertex stuff
    for (unsigned int iGV = 0; iGV < genVertices.size(); iGV++) {
        gv_nTimesMatchedToSVinJet_->push_back(0);
    }
    for (unsigned int iSV = 0; iSV < cpvtx.size(); iSV++) {
        if (iSV >= max_sv_) break;
        const reco::VertexCompositePtrCandidate& sv = cpvtx.at(iSV);
        sv_x_->push_back(sv.vertex().x());
        sv_y_->push_back(sv.vertex().y());
        sv_z_->push_back(sv.vertex().z());
        sv_pt_->push_back(sv.pt());
        sv_eta_->push_back(sv.eta());
        sv_phi_->push_back(sv.phi());
        sv_mass_->push_back(sv.mass());
        sv_energy_->push_back(sv.energy());
        sv_chi2_->push_back(sv.vertexChi2());
        sv_ndof_->push_back(sv.vertexNdof());
        sv_chi2dof_->push_back(catchInfsAndBound(sv.vertexChi2() / sv.vertexNdof(), 1000, -1000, 1000));
        sv_dxy_->push_back(vertexDxy(sv, pv).value());
        sv_dxyerr_->push_back(catchInfsAndBound(vertexDxy(sv, pv).error() - 2, 0, -2, 0));
        sv_dxysig_->push_back(catchInfsAndBound(vertexDxy(sv, pv).value() / vertexDxy(sv, pv).error(), 0, -1, 800));
        sv_dz_->push_back(TMath::Abs(sv.vertex().z() - pv.z()));
        sv_d3D_->push_back(vertexD3d(sv, pv).value());
        sv_d3Derr_->push_back(catchInfsAndBound(vertexD3d(sv, pv).error() - 2, 0, -2, 0));
        sv_d3Dsig_->push_back(catchInfsAndBound(vertexD3d(sv, pv).value() / vertexD3d(sv, pv).error(), 0, -1, 800));
        sv_nDaughters_->push_back(sv.numberOfDaughters());

        // Match SV to GV
        int closestGVIdx = -1;
        float mindR = 99999999.99; // some maximum value
        for (unsigned int iGV = 0; iGV < genVertices.size(); iGV++) {
            const GenVertex* genVtx = genVertices.at(iGV);
            float dR3D = deltaR3D(genVtx->x(), sv.vertex().x(), genVtx->y(), sv.vertex().y(), genVtx->z(), sv.vertex().z());
            if ((dR3D < matchGVdR_) && (dR3D < mindR)) {
                closestGVIdx = iGV;
                mindR = dR3D;
            }
        }
        if (closestGVIdx >= 0) {
            sv_matchGV_x_->push_back(sv.vertex().x());
            sv_matchGV_y_->push_back(sv.vertex().y());
            sv_matchGV_z_->push_back(sv.vertex().z());
            sv_matchGV_pt_->push_back(sv.pt());
            sv_matchGV_eta_->push_back(sv.eta());
            sv_matchGV_phi_->push_back(sv.phi());
            sv_matchGV_mass_->push_back(sv.mass());
            sv_matchGV_energy_->push_back(sv.energy());
            sv_matchGV_chi2_->push_back(sv.vertexChi2());
            sv_matchGV_ndof_->push_back(sv.vertexNdof());
            sv_matchGV_chi2dof_->push_back(catchInfsAndBound(sv.vertexChi2() / sv.vertexNdof(), 1000, -1000, 1000));
            sv_matchGV_dxy_->push_back(vertexDxy(sv, pv).value());
            sv_matchGV_dxyerr_->push_back(catchInfsAndBound(vertexDxy(sv, pv).error() - 2, 0, -2, 0));
            sv_matchGV_dxysig_->push_back(catchInfsAndBound(vertexDxy(sv, pv).value() / vertexDxy(sv, pv).error(), 0, -1, 800));
            sv_matchGV_dz_->push_back(TMath::Abs(sv.vertex().z() - pv.z()));
            sv_matchGV_d3D_->push_back(vertexD3d(sv, pv).value());
            sv_matchGV_d3Derr_->push_back(catchInfsAndBound(vertexD3d(sv, pv).error() - 2, 0, -2, 0));
            sv_matchGV_d3Dsig_->push_back(catchInfsAndBound(vertexD3d(sv, pv).value() / vertexD3d(sv, pv).error(), 0, -1, 800));
            sv_matchGV_SVdRtoGV_->push_back(mindR);
            sv_matchGV_nDaughters_->push_back(sv.numberOfDaughters());
            sv_matchGV_GVmotherPdgId_->push_back(genVertices.at(closestGVIdx)->motherPdgId());

            // Match SV to Jet
            unsigned int nJets = 0;
            for (unsigned int iJet = 0; iJet < jetCollection.size(); iJet++) {
                const pat::Jet& jet = jetCollection.at(iJet);

                // Some cuts
                // if (jet.pt() < jetPtMin_) continue;
                // if (jet.pt() > jetPtMax_) continue;
                // if (TMath::Abs(jet.eta()) < jetAbsEtaMin_) continue;
                // if (TMath::Abs(jet.eta()) > jetAbsEtaMax_) continue;

                // Calculate Jet radius
                double jet_radius = jetR();
                if (jet_radius < 0) {
                    for (unsigned idau = 0; idau < jet.numberOfDaughters(); ++idau) {
                        double dR = reco::deltaR(*jet.daughter(idau), jet);
                        if (dR > jet_radius)
                            jet_radius = dR;
                    }
                }

                // Get GenJet flavours (might be the same as reco?)
                int genJetHadFlav = -1000;
                // int genJetPartFlav = -1000;
                const reco::GenJet* genJet = jet.genJet();
                if (genJet) {
                    for (const reco::JetFlavourInfoMatching& genJetFlavInfo : *(genJetFlavourInfo_.product())) {
                        if (reco::deltaR(genJet->p4(), genJetFlavInfo.first->p4()) < genJetMatchdR_) {
                            genJetHadFlav = genJetFlavInfo.second.getHadronFlavour();
                            // genJetPartFlav = genJetFlavInfo.second.getPartonFlavour();
                            break;
                        }
                    }
                }

                if (reco::deltaR(sv, jet) < jet_radius) {
                    nJets++;
                    sv_matchGV_matchJet_genJetHadFlav_->push_back(genJetHadFlav);
                }
            }

            if (nJets > 0) {
                (gv_nTimesMatchedToSVinJet_->at(closestGVIdx))++;
                sv_matchGV_matchJet_x_->push_back(sv.vertex().x());
                sv_matchGV_matchJet_y_->push_back(sv.vertex().y());
                sv_matchGV_matchJet_z_->push_back(sv.vertex().z());
                sv_matchGV_matchJet_pt_->push_back(sv.pt());
                sv_matchGV_matchJet_eta_->push_back(sv.eta());
                sv_matchGV_matchJet_phi_->push_back(sv.phi());
                sv_matchGV_matchJet_mass_->push_back(sv.mass());
                sv_matchGV_matchJet_energy_->push_back(sv.energy());
                sv_matchGV_matchJet_chi2_->push_back(sv.vertexChi2());
                sv_matchGV_matchJet_ndof_->push_back(sv.vertexNdof());
                sv_matchGV_matchJet_chi2dof_->push_back(catchInfsAndBound(sv.vertexChi2() / sv.vertexNdof(), 1000, -1000, 1000));
                sv_matchGV_matchJet_dxy_->push_back(vertexDxy(sv, pv).value());
                sv_matchGV_matchJet_dxyerr_->push_back(catchInfsAndBound(vertexDxy(sv, pv).error() - 2, 0, -2, 0));
                sv_matchGV_matchJet_dxysig_->push_back(catchInfsAndBound(vertexDxy(sv, pv).value() / vertexDxy(sv, pv).error(), 0, -1, 800));
                sv_matchGV_matchJet_dz_->push_back(TMath::Abs(sv.vertex().z() - pv.z()));
                sv_matchGV_matchJet_d3D_->push_back(vertexD3d(sv, pv).value());
                sv_matchGV_matchJet_d3Derr_->push_back(catchInfsAndBound(vertexD3d(sv, pv).error() - 2, 0, -2, 0));
                sv_matchGV_matchJet_d3Dsig_->push_back(catchInfsAndBound(vertexD3d(sv, pv).value() / vertexD3d(sv, pv).error(), 0, -1, 800));
                sv_matchGV_matchJet_SVdRtoGV_->push_back(mindR);
                sv_matchGV_matchJet_nJets_->push_back(nJets);
                sv_matchGV_matchJet_nDaughters_->push_back(sv.numberOfDaughters());
                sv_matchGV_matchJet_GVmotherPdgId_->push_back(genVertices.at(closestGVIdx)->motherPdgId());
            }
        }
    }

    // Jet stuff
    for (unsigned int iJet = 0; iJet < jetCollection.size(); iJet++) {
        const pat::Jet& jet = jetCollection.at(iJet);

        // Some cuts
        // if (jet.pt() < jetPtMin_) continue;
        // if (jet.pt() > jetPtMax_) continue;
        // if (TMath::Abs(jet.eta()) < jetAbsEtaMin_) continue;
        // if (TMath::Abs(jet.eta()) > jetAbsEtaMax_) continue;

        // Calculate Jet radius
        double jet_radius = jetR();
        if (jet_radius < 0) {
            for (unsigned idau = 0; idau < jet.numberOfDaughters(); ++idau) {
                double dR = reco::deltaR(*jet.daughter(idau), jet);
                if (dR > jet_radius)
                    jet_radius = dR;
            }
        }

        // Get Gen Jet flavours (might be the same as reco?)
        int genJetHadFlav = -1000;
        int genJetPartFlav = -1000;
        const reco::GenJet* genJet = jet.genJet();
        if (genJet) {
            for (const reco::JetFlavourInfoMatching& genJetFlavInfo : *(genJetFlavourInfo_.product())) {
                if (reco::deltaR(genJet->p4(), genJetFlavInfo.first->p4()) < genJetMatchdR_) {
                    genJetHadFlav = genJetFlavInfo.second.getHadronFlavour();
                    genJetPartFlav = genJetFlavInfo.second.getPartonFlavour();
                    break;
                }
            }
        }

        jet_pt_->push_back(jet.pt());
        jet_eta_->push_back(jet.eta());
        jet_phi_->push_back(jet.phi());
        jet_radius_->push_back(jet_radius);
        jet_hadFlav_->push_back(jet.hadronFlavour());
        jet_partFlav_->push_back(jet.partonFlavour());
        jet_genHadFlav_->push_back(genJetHadFlav);
        jet_genPartFlav_->push_back(genJetPartFlav);

        // Match Jet to SV
        bool matchedToSV = false;
        bool matchedToGV = true;
        unsigned int nSV = 0;
        unsigned int nGV = 0;
        for (unsigned int iSV = 0; iSV < cpvtx.size(); iSV++) {
            if (iSV >= max_sv_) break;
            const reco::VertexCompositePtrCandidate& sv = cpvtx.at(iSV);
            if (reco::deltaR(sv, jet) < jet_radius) {
                matchedToSV = true;
                nSV++;

                // Match (jet-matched) SV to GV
                int closestGVIdx = -1;
                float mindR = 99999999.99; // some maximum value
                for (unsigned int iGV = 0; iGV < genVertices.size(); iGV++) {
                    const GenVertex* genVtx = genVertices.at(iGV);
                    float dR3D = deltaR3D(genVtx->x(), sv.vertex().x(), genVtx->y(), sv.vertex().y(), genVtx->z(), sv.vertex().z());
                    if ((dR3D < matchGVdR_) && (dR3D < mindR)) {
                        closestGVIdx = iGV;
                        mindR = dR3D;
                    }
                }
                if (closestGVIdx >= 0) {
                    matchedToGV = true;
                    nGV++;
                    jet_matchSV_matchGV_SVdRtoGV_->push_back(mindR);
                }
            }
        }

        if (matchedToSV) {
            jet_matchSV_pt_->push_back(jet.pt());
            jet_matchSV_eta_->push_back(jet.eta());
            jet_matchSV_phi_->push_back(jet.phi());
            jet_matchSV_radius_->push_back(jet_radius);
            jet_matchSV_nSV_->push_back(nSV);
            jet_matchSV_hadFlav_->push_back(jet.hadronFlavour());
            jet_matchSV_partFlav_->push_back(jet.partonFlavour());
            jet_matchSV_genHadFlav_->push_back(genJetHadFlav);
            jet_matchSV_genPartFlav_->push_back(genJetPartFlav);
        }
        if (matchedToSV && matchedToGV) {
            jet_matchSV_matchGV_pt_->push_back(jet.pt());
            jet_matchSV_matchGV_eta_->push_back(jet.eta());
            jet_matchSV_matchGV_phi_->push_back(jet.phi());
            jet_matchSV_matchGV_radius_->push_back(jet_radius);
            jet_matchSV_matchGV_nGV_->push_back(nGV);
            jet_matchSV_matchGV_hadFlav_->push_back(jet.hadronFlavour());
            jet_matchSV_matchGV_partFlav_->push_back(jet.partonFlavour());
            jet_matchSV_matchGV_genHadFlav_->push_back(genJetHadFlav);
            jet_matchSV_matchGV_genPartFlav_->push_back(genJetPartFlav);
        }
    }

    // Fill track timing information
    for (unsigned int iPV = 0; iPV < PVs.size(); iPV++) {
        const reco::Vertex& prv = PVs.at(iPV);
        if (prv.isFake()) std::cout << "ntuple_SV::filBranches: This is a fake PV." << std::endl;
        for (reco::Vertex::trackRef_iterator trk_it = prv.tracks_begin(); trk_it != prv.tracks_end(); trk_it++) {
            reco::TrackBaseRef trkRef = *trk_it;
            if (timeValueMap.contains(trkRef.id())) {
                // if (timeQualityMap[trkRef] < timeQualityCut_) continue;
                pv_trk_timeVal_->push_back(timeValueMap[trkRef]);
                pv_trk_timeErr_->push_back(timeErrorMap[trkRef]);
                pv_trk_timeQual_->push_back(timeQualityMap[trkRef]);
            }
            else std::cout << "ntuple_SV::filBranches: PV track is not in the track maps." << std::endl;
        }
    }
    for (unsigned int iSV = 0; iSV < inclusiveSVs.size(); iSV++) {
        const reco::Vertex& sv = inclusiveSVs.at(iSV);
        for (reco::Vertex::trackRef_iterator trk_it = sv.tracks_begin(); trk_it != sv.tracks_end(); trk_it++) {
            reco::TrackBaseRef trkRef = *trk_it;
            if (timeValueMap.contains(trkRef.id())) {
                // if (timeQualityMap[trkRef] < timeQualityCut_) continue;
                sv_trk_timeVal_->push_back(timeValueMap[trkRef]);
                sv_trk_timeErr_->push_back(timeErrorMap[trkRef]);
                sv_trk_timeQual_->push_back(timeQualityMap[trkRef]);
            }
            else std::cout << "ntuple_SV::filBranches: SV track is not in the track maps." << std::endl;
        }
    }

    sv_num_ = 0;
    for (const reco::VertexCompositePtrCandidate& sv : cpvtx) {
        if (sv_num_ < (int) max_sv_) { // Limit number of SVs

            // sv_pt_[sv_num_] = sv.pt();
            // sv_eta_[sv_num_] = sv.eta();
            // sv_phi_[sv_num_] = sv.phi();
            // sv_mass_[sv_num_] = sv.mass();
            // sv_e_[sv_num_] = sv.energy();
            sv_ntracks_[sv_num_] = sv.numberOfDaughters();
            // sv_chi2_[sv_num_] = sv.vertexChi2();
            // sv_ndf_[sv_num_] = sv.vertexNdof();
            // sv_normchi2_[sv_num_] = catchInfsAndBound(sv_chi2_[sv_num_] / sv_ndf_[sv_num_], 1000, -1000, 1000);
            // sv_dxy_[sv_num_] = vertexDxy(sv, pv).value();
            // sv_dxyerr_[sv_num_] = catchInfsAndBound(vertexDxy(sv, pv).error() - 2, 0, -2, 0);
            // sv_dxysig_[sv_num_] = catchInfsAndBound(sv_dxy_[sv_num_] / vertexDxy(sv, pv).error(), 0, -1, 800);
            // sv_d3d_[sv_num_] = vertexD3d(sv, pv).value();
            // sv_d3derr_[sv_num_] = catchInfsAndBound(vertexD3d(sv, pv).error() - 2, 0, -2, 0);
            // sv_d3dsig_[sv_num_] = catchInfsAndBound(vertexD3d(sv, pv).value() / vertexD3d(sv, pv).error(), 0, -1, 800);
            sv_costhetasvpv_[sv_num_] = vertexDdotP(sv, pv); // The pointing angle (i.e. the angle between the sum of the momentum of the tracks in the SV and the flight direction betwen PV and SV)

            SVTrackInfoBuilder trackinfo(builder_);
            float calo_frac = 0.0;
            float hcal_frac = 0.0;
            float puppiw = 0.0;
            float dz = 0.0;
            float charge = 0.0;
            for (unsigned idx = 0; idx < sv.numberOfDaughters(); ++idx) {
                const pat::PackedCandidate* PackedCandidate_ = dynamic_cast<const pat::PackedCandidate*>(sv.daughter(idx));
                calo_frac = calo_frac + PackedCandidate_->caloFraction();
                hcal_frac = hcal_frac + PackedCandidate_->hcalFraction();
                puppiw = puppiw + PackedCandidate_->puppiWeight();
                dz = dz + PackedCandidate_->dz();
                charge = charge + PackedCandidate_->charge();
            }
            sv_calo_frac_[sv_num_] = calo_frac / sv.numberOfDaughters();
            sv_hcal_frac_[sv_num_] = hcal_frac / sv.numberOfDaughters();
            sv_puppiw_[sv_num_] = puppiw / sv.numberOfDaughters();
            // sv_dz_[sv_num_] = dz / sv.numberOfDaughters();
            sv_charge_sum_[sv_num_] = charge;

            for (unsigned int j = 0; j < jetCollection.size(); j++) {

                const pat::Jet& jet = jetCollection.at(j);

                // Match to a jet based on deltaR
                double jet_radius = jetR();
                if (jet_radius < 0) {
                    // Subjets: use maxDR(subjet, pfcand)
                    for (unsigned idau = 0; idau < jet.numberOfDaughters(); ++idau) {
                        double dR = reco::deltaR(*jet.daughter(idau), jet);
                        if (dR > jet_radius)
                            jet_radius = dR;
                    }
                }
                if (reco::deltaR(sv, jet) > jet_radius) continue;

                sv_etarel_->push_back(catchInfsAndBound(fabs(sv.eta() - jet.eta()) - 0.5, 0, -2, 0));
                sv_phirel_->push_back(catchInfsAndBound(fabs(reco::deltaPhi(sv.phi(), jet.phi())) - 0.5, 0, -2, 0));
                sv_deltaR_->push_back(catchInfsAndBound(fabs(reco::deltaR(sv, jet)) - 0.5, 0, -2, 0));
                sv_enratio_->push_back(sv.energy() / jet.correctedJet("Uncorrected").energy());

                math::XYZVector jetDir = jet.momentum().Unit();
                GlobalVector jetRefTrackDir(jet.px(), jet.py(), jet.pz());
                float pfd3dval = 0.0;
                float pfd3dsig = 0.0;
                float pfd2dval = 0.0;
                float pfd2dsig = 0.0;
                float pfcount = 0.0;
                for (unsigned idx = 0; idx < sv.numberOfDaughters(); ++idx) {
                    const pat::PackedCandidate* PackedCandidate_ = dynamic_cast<const pat::PackedCandidate*>(sv.daughter(idx));
                    if (PackedCandidate_->charge() != 0 and PackedCandidate_->pt() > 0.95) { // TODO: understand these "track" cuts
                        trackinfo.buildTrackInfo(PackedCandidate_, jetDir, jetRefTrackDir, pv);
                        pfd3dval = pfd3dval + catchInfsAndBound(trackinfo.getTrackSip3dVal(), 0, -1, 1e5);
                        pfd3dsig = pfd3dsig + catchInfsAndBound(trackinfo.getTrackSip3dSig(), 0, -1, 4e4);
                        pfd2dval = pfd2dval + catchInfsAndBound(trackinfo.getTrackSip2dVal(), 0, -1, 70);
                        pfd2dsig = pfd2dsig + catchInfsAndBound(trackinfo.getTrackSip2dSig(), 0, -1, 4e4);
                        pfcount = pfcount + 1.0;
                    }
                }
                sv_pfd3dval_->push_back(pfd3dval / pfcount);
                sv_pfd3dsig_->push_back(pfd3dsig / pfcount);
                sv_pfd2dval_->push_back(pfd2dval / pfcount);
                sv_pfd2dsig_->push_back(pfd2dsig / pfcount);

                // Get SV tag info
                int nSV = -2;
                const reco::CandSecondaryVertexTagInfo* candSVTagInfo = jet.tagInfoCandSecondaryVertex("pfInclusiveSecondaryVertexFinder");
                if (candSVTagInfo != nullptr) nSV = candSVTagInfo->nVertices();
                if (nSV > 0 && candSVTagInfo->vertexTracks().size() == 0) nSV = -1;
                // Loop over SV tag info
                float vertex_time = 0;
                float vertex_timeWeight = 0;
                float vertex_timeNtk = 0;
                if (nSV > 0 && sv.pt() > 0.0) {
                    for (unsigned int isv = 0; isv < candSVTagInfo->nVertices(); ++isv) {
                        float dSVpt = TMath::Abs(candSVTagInfo->secondaryVertex(isv).pt() / sv.pt() - 1.0);
                        float dSVeta = TMath::Abs(candSVTagInfo->secondaryVertex(isv).eta() - sv.eta());
                        float dSVphi = TMath::Abs(candSVTagInfo->secondaryVertex(isv).phi() - sv.phi());
                        if (dSVphi > 3.141593) dSVphi -= 2.0 * 3.141593; // Make sure phi is within -pi to pi
                        // Match reco SV to SV tag info (to combine information from both collections?)
                        if (dSVpt > 0.01 || dSVeta > 0.01 || dSVphi > 0.01) continue;

                        for (unsigned int it = 0; it < candSVTagInfo->nVertexTracks(isv); ++it) {
                            for (unsigned int i = 0; i < jet.numberOfDaughters(); i++) {
                                // Get best track in jet track
                                const pat::PackedCandidate* PackedCandidate = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(i));
                                if (!PackedCandidate) continue;
                                if (PackedCandidate->charge() == 0) continue; // Take out photons
                                auto track = PackedCandidate->bestTrack();
                                if (!track) continue;

                                // Check charge of SV tag info track and jet track match
                                if (candSVTagInfo->vertexTracks(isv)[it]->charge() != track->charge()) continue;

                                float track_time = track->t0();
                                float track_timeError = track->covt0t0();
                                float track_pt = track->pt();
                                float time_weight = track_pt * track_pt;
                                if (track_timeError < 0.0 || abs(track_time) > 1) continue;

                                float dpt = TMath::Abs(candSVTagInfo->vertexTracks(isv)[it]->pt() / track->pt() - 1.0);
                                float deta = TMath::Abs(candSVTagInfo->vertexTracks(isv)[it]->eta() - track->eta());
                                float dphi = TMath::Abs(candSVTagInfo->vertexTracks(isv)[it]->phi() - track->phi());
                                if (dphi > 3.141593) dphi -= 2.0 * 3.141593; // Make sure phi is within -pi to pi
                                // Match SV track to jet track
                                if (dpt < 0.01 && deta < 0.01 && dphi < 0.01) {
                                    vertex_timeNtk += 1;
                                    vertex_timeWeight += time_weight;
                                    vertex_time += track_time * time_weight;
                                    // std::cout << "  => matched track " << it << " to " << i << " time " << track_time << std::endl;
                                }
                            } // End loop on all tracks in jet
                        } // End loop on tracks from SV in jet
                    } // End loop on SVs in jet
                    if (vertex_timeNtk > 0 && EventTime > -1) {
                        vertex_time = vertex_time / vertex_timeWeight - EventTime;
                        vertex_time = TMath::Abs(vertex_time);
                    } // Time of flight?
                    else vertex_time = -1.0; 
                } // End if nSV > 0 && sv.pt() > 0.0
                else vertex_time = -1.0;
                // std::cout << " NTuple sv " << sv_num_ << " pt eta phi " << sv.pt() << " " << sv.eta() << " " << sv.phi() << " time " << vertex_time << std::endl;
                sv_time_->push_back(vertex_time);
            } // End loop through jets
        } // End if sv_num_ < max_sv_
        sv_num_++;
    } // End loop through SVs
    // nsv_ = sv_num_;

    return 0;
}


// Helper functions


Measurement1D ntuple_SV::vertexDxy(const reco::VertexCompositePtrCandidate& svcand, const reco::Vertex& pv) {

    VertexDistanceXY dist;
    reco::Vertex::CovarianceMatrix csv;
    svcand.fillVertexCovariance(csv);
    reco::Vertex svtx(svcand.vertex(), csv);
    return dist.distance(svtx, pv);
}


Measurement1D ntuple_SV::vertexD3d(const reco::VertexCompositePtrCandidate& svcand, const reco::Vertex& pv) {

    VertexDistance3D dist;
    reco::Vertex::CovarianceMatrix csv;
    svcand.fillVertexCovariance(csv);
    reco::Vertex svtx(svcand.vertex(), csv);
    return dist.distance(svtx, pv);
}


bool ntuple_SV::compareDxyDxyErr(const reco::VertexCompositePtrCandidate& sva, const reco::VertexCompositePtrCandidate& svb) {

    reco::Vertex pv = *spvp_;
    float adxy = ntuple_SV::vertexDxy(sva, pv).value();
    float bdxy = ntuple_SV::vertexDxy(svb, pv).value();
    float aerr = ntuple_SV::vertexDxy(sva, pv).error();
    float berr = ntuple_SV::vertexDxy(svb, pv).error();

    float asig = ntuple_SV::catchInfs(adxy / aerr, 0.0);
    float bsig = ntuple_SV::catchInfs(bdxy / berr, 0.0);

    return bsig < asig;
}


float ntuple_SV::vertexDdotP(const reco::VertexCompositePtrCandidate& sv, const reco::Vertex& pv) {

    reco::Candidate::Vector p = sv.momentum();
    reco::Candidate::Vector d(sv.vx() - pv.x(), sv.vy() - pv.y(), sv.vz() - pv.z());
    return p.Unit().Dot(d.Unit());
}


float ntuple_SV::deltaR3D(float x1, float x2, float y1, float y2, float z1, float z2) {

    return TMath::Sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
}


int ntuple_SV::genPartID(int pdgId) {

    int checkPdgId = abs(pdgId);
    // B meson
    if ((checkPdgId == 521) ||
        (checkPdgId == 511) ||
        (checkPdgId == 531) ||
        (checkPdgId == 541)) return 0;
    // B baryon
    else if (
        (checkPdgId == 5122) ||
        (checkPdgId == 5112) ||
        (checkPdgId == 5212) ||
        (checkPdgId == 5222) ||
        (checkPdgId == 5132) ||
        (checkPdgId == 5232) ||
        (checkPdgId == 5332) ||
        (checkPdgId == 5142) ||
        (checkPdgId == 5242) ||
        (checkPdgId == 5342) ||
        (checkPdgId == 5512) ||
        (checkPdgId == 5532) ||
        (checkPdgId == 5542) ||
        (checkPdgId == 5554)) return 1;
    // C meson
    else if (
        (checkPdgId == 411) ||
        (checkPdgId == 421) ||
        (checkPdgId == 431)) return 2;
    // C baryon
    else if (
        (checkPdgId == 4122) ||
        (checkPdgId == 4222) ||
        (checkPdgId == 4212) ||
        (checkPdgId == 4112) ||
        (checkPdgId == 4232) ||
        (checkPdgId == 4132) ||
        (checkPdgId == 4332) ||
        (checkPdgId == 4412) ||
        (checkPdgId == 4422) ||
        (checkPdgId == 4432) ||
        (checkPdgId == 4444)) return 3;
    // S baryon
    else if (
        (checkPdgId == 3122) ||
        (checkPdgId == 3222) ||
        (checkPdgId == 3212) ||
        (checkPdgId == 3112) ||
        (checkPdgId == 3322) ||
        (checkPdgId == 3312) ||
        (checkPdgId == 3334)) return 4;
    return -1;
}


// int ntuple_SV::findPFCandIdx(const pat::PackedCandidate& trk, const pat::PackedCandidateCollection& pcands) {

//     int nmatches = 0;
//     int matchIdx = 0;
//     for (unsigned int trkIdx = 0; trkIdx < pcands.size(); trkIdx++) {
//         if (&trk == &pcands.at(trkIdx)) {
//             matchIdx = trkIdx;
//             nmatches += 1;
//         }
//     }
//     if (nmatches != 1) {
//         // std::cout << "ntuple_SV.cc: 0 or more than one PF Candidate match found! Returning -1." << std::endl;
//         return -1;
//     }
//     return matchIdx;
// }
