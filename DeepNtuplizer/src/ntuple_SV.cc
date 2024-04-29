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

        GenVertex() {}

        const float x() const { return daughters_.at(0)->vx(); }
        const float y() const { return daughters_.at(0)->vy(); }
        const float z() const { return daughters_.at(0)->vz(); }
        const float pt() const { return mother_->pt(); }
        const float eta() const { return mother_->eta(); }
        const float phi() const { return mother_->phi(); }
        const unsigned int nDaughters() const { return daughters_.size(); }
        const int motherPdgId() const { return mother_->pdgId(); }

        const reco::GenParticle* mother() const { return mother_; }
        const std::vector<const reco::Candidate*> daughters() const { return daughters_; }

        void setMother(const reco::GenParticle* mother) { mother_ = mother; }
        void addDaughter(const reco::Candidate* daughter) { daughters_.push_back(daughter); }
        bool isValid() {
            bool valid = true;
            if (!mother_) valid = false;
            if (daughters_.size() < 2) valid = false;
            return valid;
        }

        void print() {
            std::cout << "GenVertex:" << std::endl;
            std::cout << "    vertex x        = " << x() << std::endl;
            std::cout << "    vertex y        = " << y() << std::endl;
            std::cout << "    vertex z        = " << z() << std::endl;
            std::cout << "    vertex pt       = " << pt() << std::endl;
            std::cout << "    vertex eta      = " << eta() << std::endl;
            std::cout << "    vertex phi      = " << phi() << std::endl;
            std::cout << "    nDaughters      = " << nDaughters() << std::endl;
            std::cout << "    mother pdg id   = " << motherPdgId() << std::endl;

            std::cout << "    daughter pdgIds = ";
            for (unsigned int iDau = 0; iDau < nDaughters(); iDau++) {
                std::cout << daughters().at(iDau)->pdgId();
                if (iDau < nDaughters() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
            std::cout << "    daughter pts   = ";
            for (unsigned int iDau = 0; iDau < nDaughters(); iDau++) {
                std::cout << daughters().at(iDau)->pt();
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

        const reco::GenParticle* mother_;
        std::vector<const reco::Candidate*> daughters_;
};


const reco::Vertex* ntuple_SV::spvp_;


ntuple_SV::ntuple_SV(std::string prefix, double jetR) : ntuple_content(jetR), sv_num_(0) {

    prefix_ = prefix;
}


ntuple_SV::~ntuple_SV() {}


void ntuple_SV::getInput(const edm::ParameterSet& iConfig) {

    absEtaMin_ = iConfig.getParameter<double>("absEtaMin");
    absEtaMax_ = iConfig.getParameter<double>("absEtaMax");
    genPartPtCut_ = iConfig.getParameter<double>("genPartPtCut");
    genDauPtCut_ = iConfig.getParameter<double>("genDauPtCut");
    trackPtCut_ = iConfig.getParameter<double>("trackPtCut");
    timeQualityCut_ = iConfig.getParameter<double>("timeQualityCut");
    matchGVdR_ = iConfig.getParameter<double>("matchGVdR");
    jetPtMin_ = iConfig.getParameter<double>("jetPtMin");
    jetPtMax_ = iConfig.getParameter<double>("jetPtMax");
    genJetMatchdR_ = iConfig.getParameter<double>("genJetMatchdR");
}


void ntuple_SV::initBranches(TTree* tree) {

    // for (TString collection : gp_collections_) {
    //     for (TString branch : gp_branches_) {
    //         TString name = prefix_ + "gp_" + collection + "_" + branch;
    //         n_.push_back(name);
    //     }
    // }
    for (TString collection : trk_collections_) {
        for (TString branch : trk_branches_) {
            TString name = prefix_ + "trk_" + collection + "_" + branch;
            n_.push_back(name);
        }
    }
    for (TString collection : vtx_collections_) {
        for (TString branch : vtx_branches_) {
            TString name = prefix_ + "vtx_" + collection + "_" + branch;
            n_.push_back(name);
        }
    }
    for (TString collection : jet_collections_) {
        for (TString branch : jet_branches_) {
            TString name = prefix_ + "jet_" + collection + "_" + branch;
            n_.push_back(name);
        }
    }
    for (TString branch : evt_branches_) {
        TString name = prefix_ + "evt_" + branch;
        n_.push_back(name);
    }

    for (TString name : n_) {
        b_[name] = new std::vector<float>;
        addBranch(tree, name, &b_[name]);
    }

    // SV candidates
    addBranch(tree, (prefix_ + "n_sv").c_str(), &sv_num_, (prefix_ + "sv_num_/I").c_str());
    // addBranch(tree, (prefix_ + "nsv").c_str(), &nsv_, (prefix_ + "nsv_/F").c_str());
    // addBranch(tree, (prefix_ + "sv_pt").c_str(), &sv_pt_, (prefix_ + "sv_pt_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_eta").c_str(), &sv_eta_, (prefix_ + "sv_eta_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_phi").c_str(), &sv_phi_, (prefix_ + "sv_phi_[" + prefix_ + "sv_num_]/F").c_str());
    addBranch(tree, (prefix_ + "sv_e").c_str(), &sv_e_, (prefix_ + "sv_e_[" + prefix_ + "sv_num_]/F").c_str());
    addBranch(tree, (prefix_ + "sv_etarel").c_str(), &sv_etarel_);
    addBranch(tree, (prefix_ + "sv_phirel").c_str(), &sv_phirel_);
    addBranch(tree, (prefix_ + "sv_deltaR").c_str(), &sv_deltaR_);
    addBranch(tree, (prefix_ + "sv_mass").c_str(), &sv_mass_, (prefix_ + "sv_mass_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_ntracks").c_str(), &sv_ntracks_, (prefix_ + "sv_ntracks_[" + prefix_ + "sv_num_]/F").c_str());
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
    // addBranch(tree, (prefix_ + "sv_time").c_str(), &sv_time_);
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
    iEvent.getByToken(trackMCMatch_token_, trackMCMatch_);
    iEvent.getByToken(timeValueMap_token_, timeValueMap_);
    iEvent.getByToken(timeErrorMap_token_, timeErrorMap_);
    iEvent.getByToken(timeQualityMap_token_, timeQualityMap_);
    // iEvent.getByToken(genVertices_token_, genVertices_);
    iEvent.getByToken(PVs_token_, PVs_);
    iEvent.getByToken(inclusiveSVs_token_, inclusiveSVs_);
    iEvent.getByToken(inclusiveSVsMTDTiming_token_, inclusiveSVsMTDTiming_);
}


void ntuple_SV::initContainers() {

    for (TString name : n_) b_[name] = new std::vector<float>;

    sv_etarel_ = new std::vector<float>;
    sv_phirel_ = new std::vector<float>;
    sv_deltaR_ = new std::vector<float>;
    sv_enratio_ = new std::vector<float>;
    sv_pfd2dval_ = new std::vector<float>;
    sv_pfd2dsig_ = new std::vector<float>;
    sv_pfd3dval_ = new std::vector<float>;
    sv_pfd3dsig_ = new std::vector<float>;
    // sv_time_ = new std::vector<float>;
}


void ntuple_SV::clearContainers() {

    for (TString name : n_) b_[name]->clear();

    sv_etarel_->clear();
    sv_phirel_->clear();
    sv_deltaR_->clear();
    sv_enratio_->clear();
    sv_pfd2dval_->clear();
    sv_pfd2dsig_->clear();
    sv_pfd3dval_->clear();
    sv_pfd3dsig_->clear();
    // sv_time_->clear();
}


void ntuple_SV::deleteContainers() {

    for (TString name : n_) delete b_[name];

    delete sv_etarel_;
    delete sv_phirel_;
    delete sv_deltaR_;
    delete sv_enratio_;
    delete sv_pfd2dval_;
    delete sv_pfd2dsig_;
    delete sv_pfd3dval_;
    delete sv_pfd3dsig_;
    // delete sv_time_;
}


int ntuple_SV::fillBranches(bool applySelection, float EventTime) {

    clearContainers();

    // Sorting described here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideOfflinePrimaryVertexProduction
    const reco::Vertex& pv = vertices()->at(0); // Most likely the signal vertex
    spvp_ = &vertices()->at(0);

    reco::VertexCompositePtrCandidateCollection cpvtx = *secVertices();
    const edm::View<pat::Jet> jetCollection = *jets();
    edm::SimTrackContainer simTracks = *(simTracks_.product());
    reco::TrackCollection recoTracks = *(recoTracks_.product());
    const edm::Association<reco::GenParticleCollection> trackMCMatch = *(trackMCMatch_.product());
    const edm::ValueMap<float>& timeValueMap = *(timeValueMap_.product());
    const edm::ValueMap<float>& timeErrorMap = *(timeErrorMap_.product());
    const edm::ValueMap<float>& timeQualityMap = *(timeQualityMap_.product());
    reco::VertexCollection PVs = *(PVs_.product()); // Not the slimmed collection
    reco::VertexCollection inclusiveSVs = *(inclusiveSVs_.product());
    reco::VertexCollection inclusiveSVsMTDTiming = *(inclusiveSVsMTDTiming_.product());

    // Sort SVs by dxy significance to PV
    std::sort(cpvtx.begin(), cpvtx.end(), ntuple_SV::candCompareDxyDxyErr);
    std::sort(inclusiveSVs.begin(), inclusiveSVs.end(), ntuple_SV::vertexCompareDxyDxyErr);
    std::sort(inclusiveSVsMTDTiming.begin(), inclusiveSVsMTDTiming.end(), ntuple_SV::vertexCompareDxyDxyErr);

    // Construct GenVertex collection
    std::vector<GenVertex*> genVertices(0);
    for (unsigned int iGP = 0; iGP < genParticles_->size(); iGP++) {
        const reco::GenParticle gp = genParticles_->at(iGP);
        int motherPartID = genPartID(gp.pdgId());
        if (motherPartID < 0) continue; // Mother is not interesting hadron
        if (gp.numberOfDaughters() < 2) continue; // Not a vertex
        // Check for last instance of interesting hadron
        bool lastInstance = true;
        for (unsigned int iDau = 0; iDau < gp.numberOfDaughters(); iDau++) {
            const reco::Candidate* dau = gp.daughter(iDau);
            int daughterPartID = genPartID(dau->pdgId());
            if (daughterPartID == motherPartID) {
                lastInstance = false; // Not last instance of interesting hadron
                break;
            }
        }
        if (lastInstance) {
            GenVertex* newGV = new GenVertex();
            newGV->setMother(gp.clone());
            for (unsigned int iDau = 0; iDau < gp.numberOfDaughters(); iDau++) {
                const reco::Candidate* dau = gp.daughter(iDau)->clone();
                newGV->addDaughter(dau);
            }
            if (newGV->isValid()) genVertices.push_back(newGV);
        }

        // if (motherPartID == 4) {
        //     GenVertex* newGenVtx = new GenVertex();
        //     newGenVtx->setMother(gp.clone());
        //     // std::cout << "NUMBER OF S BARYON DAUGHTERS = " << gp.numberOfDaughters() << std::endl;
        //     // gp_strange_status_->push_back((int) (gp.statusFlags().flags_).to_ulong());
        //     // gp_strange_pdgId_->push_back(gp.pdgId());
        //     // gp_strange_nDaughters_->push_back(gp.numberOfDaughters());
        //     // gp_strange_pt_->push_back(gp.pt());
        //     // gp_strange_eta_->push_back(gp.eta());
        //     for (unsigned int iDau = 0; iDau < gp.numberOfDaughters(); iDau++) {
        //         const reco::Candidate* dau = gp.daughter(iDau)->clone();
        //         // gp_strange_daughterPdgIds_->push_back(dau->pdgId());
        //         newGenVtx->addDaughter(dau);
        //     }
        //     newGenVtx->setGenVertexAttributes();
        //     // if (gp.numberOfDaughters() > 0)
        //     //     newGenVtx->print();
        // }
    }

    // All the matching stuff
    std::vector<int> GV_matchtoSV;
    std::vector<int> GV_matchtoSVt;
    for (unsigned int iGV = 0; iGV < genVertices.size(); iGV++) {
        GenVertex& gv = *(genVertices.at(iGV));

        int closestSVIdx = -1;
        float mindR = matchGVdR_;
        for (unsigned int iSV = 0; iSV < inclusiveSVs.size(); iSV++) {
            const reco::Vertex& sv = inclusiveSVs.at(iSV);
            float dR3D = deltaR3D(gv.x(), sv.x(), gv.y(), sv.y(), gv.z(), sv.z());
            if (dR3D < mindR) {
                closestSVIdx = iSV;
                mindR = dR3D;
            }
        }
        GV_matchtoSV.push_back(closestSVIdx);

        int closestSVtIdx = -1;
        mindR = matchGVdR_;
        for (unsigned int iSV = 0; iSV < inclusiveSVsMTDTiming.size(); iSV++) {
            const reco::Vertex& sv = inclusiveSVsMTDTiming.at(iSV);
            float dR3D = deltaR3D(gv.x(), sv.x(), gv.y(), sv.y(), gv.z(), sv.z());
            if (dR3D < mindR) {
                closestSVtIdx = iSV;
                mindR = dR3D;
            }
        }
        GV_matchtoSVt.push_back(closestSVtIdx);
    }

    for (unsigned int iST = 0; iST < simTracks.size(); iST++) {
        SimTrack st = simTracks.at(iST);
        b_["trk_st_pt"]->push_back(st.momentum().Pt());
        b_["trk_st_eta"]->push_back(st.momentum().Eta());
        b_["trk_st_phi"]->push_back(st.momentum().Phi());
        b_["trk_st_charge"]->push_back(st.charge());
    }

    for (unsigned int iRT = 0; iRT < recoTracks.size(); iRT++) {
        reco::TrackRef trkRef(recoTracks_, iRT);

        if (!goodTrack(trkRef, timeValueMap, timeErrorMap, timeQualityMap, trackPtCut_, timeQualityCut_)) continue;

        float tval = timeValueMap[trkRef];
        float terr = timeErrorMap[trkRef];
        float tsig = tval / terr;
        float tqual = timeQualityMap[trkRef];
        float d3d = TMath::Sqrt(trkRef->dxy()*trkRef->dxy() + trkRef->dz()*trkRef->dz());
        float dxy2err = 2*trkRef->dxyError()*trkRef->dxy();
        float dz2err = 2*trkRef->dzError()*trkRef->dz();
        float sum = d3d*d3d;
        float sumerr = TMath::Sqrt(dxy2err*dxy2err + dz2err*dz2err);
        float d3derr = 0.5*sumerr/TMath::Sqrt(sum);
        float d3dsig = d3d / d3derr;

        b_["trk_all_tval"]->push_back(tval);
        b_["trk_all_terr"]->push_back(terr);
        b_["trk_all_tsig"]->push_back(tsig);
        b_["trk_all_tqual"]->push_back(tqual);
        b_["trk_all_x"]->push_back(trkRef->vx());
        b_["trk_all_y"]->push_back(trkRef->vy());
        b_["trk_all_z"]->push_back(trkRef->vz());
        b_["trk_all_pt"]->push_back(trkRef->pt());
        b_["trk_all_pterr"]->push_back(trkRef->ptError());
        b_["trk_all_eta"]->push_back(trkRef->eta());
        b_["trk_all_etaerr"]->push_back(trkRef->etaError());
        b_["trk_all_phi"]->push_back(trkRef->phi());
        b_["trk_all_phierr"]->push_back(trkRef->phiError());
        b_["trk_all_dxy"]->push_back(trkRef->dxy());
        b_["trk_all_dxyerr"]->push_back(trkRef->dxyError());
        b_["trk_all_dxysig"]->push_back(trkRef->dxy() / trkRef->dxyError());
        b_["trk_all_dz"]->push_back(trkRef->dz());
        b_["trk_all_dzerr"]->push_back(trkRef->dzError());
        b_["trk_all_dzsig"]->push_back(trkRef->dz() / trkRef->dzError());
        b_["trk_all_d3d"]->push_back(d3d);
        b_["trk_all_d3derr"]->push_back(d3derr);
        b_["trk_all_d3dsig"]->push_back(d3dsig);
        b_["trk_all_d0"]->push_back(trkRef->d0());
        b_["trk_all_d0err"]->push_back(trkRef->d0Error());
        b_["trk_all_d0sig"]->push_back(trkRef->d0() / trkRef->d0Error());
        b_["trk_all_charge"]->push_back(trkRef->charge());
        b_["trk_all_chi2"]->push_back(trkRef->chi2());
        b_["trk_all_ndof"]->push_back(trkRef->ndof());
        b_["trk_all_chi2dof"]->push_back(trkRef->normalizedChi2());

        reco::GenParticleRef trkMCRef = trackMCMatch[trkRef];
        if (trkMCRef.id().isValid()) {
            // Tracked matched to GenParticle
            b_["trk_match_gp_tval"]->push_back(tval);
            b_["trk_match_gp_terr"]->push_back(terr);
            b_["trk_match_gp_tsig"]->push_back(tsig);
            b_["trk_match_gp_tqual"]->push_back(tqual);
            b_["trk_match_gp_x"]->push_back(trkRef->vx());
            b_["trk_match_gp_y"]->push_back(trkRef->vy());
            b_["trk_match_gp_z"]->push_back(trkRef->vz());
            b_["trk_match_gp_pt"]->push_back(trkRef->pt());
            b_["trk_match_gp_pterr"]->push_back(trkRef->ptError());
            b_["trk_match_gp_eta"]->push_back(trkRef->eta());
            b_["trk_match_gp_etaerr"]->push_back(trkRef->etaError());
            b_["trk_match_gp_phi"]->push_back(trkRef->phi());
            b_["trk_match_gp_phierr"]->push_back(trkRef->phiError());
            b_["trk_match_gp_dxy"]->push_back(trkRef->dxy());
            b_["trk_match_gp_dxyerr"]->push_back(trkRef->dxyError());
            b_["trk_match_gp_dxysig"]->push_back(trkRef->dxy() / trkRef->dxyError());
            b_["trk_match_gp_dz"]->push_back(trkRef->dz());
            b_["trk_match_gp_dzerr"]->push_back(trkRef->dzError());
            b_["trk_match_gp_dzsig"]->push_back(trkRef->dz() / trkRef->dzError());
            b_["trk_match_gp_d3d"]->push_back(d3d);
            b_["trk_match_gp_d3derr"]->push_back(d3derr);
            b_["trk_match_gp_d3dsig"]->push_back(d3dsig);
            b_["trk_match_gp_d0"]->push_back(trkRef->d0());
            b_["trk_match_gp_d0err"]->push_back(trkRef->d0Error());
            b_["trk_match_gp_d0sig"]->push_back(trkRef->d0() / trkRef->d0Error());
            b_["trk_match_gp_charge"]->push_back(trkRef->charge());
            b_["trk_match_gp_chi2"]->push_back(trkRef->chi2());
            b_["trk_match_gp_ndof"]->push_back(trkRef->ndof());
            b_["trk_match_gp_chi2dof"]->push_back(trkRef->normalizedChi2());

            // Check for track from pileup using GEN information
            if (!trkMCRef->isHardProcess() &&
                !trkMCRef->fromHardProcessDecayed() &&
                !trkMCRef->fromHardProcessFinalState() &&
                !trkMCRef->isDirectHardProcessTauDecayProductFinalState()) {

                b_["trk_pu_tval"]->push_back(tval);
                b_["trk_pu_terr"]->push_back(terr);
                b_["trk_pu_tsig"]->push_back(tsig);
                b_["trk_pu_tqual"]->push_back(tqual);
                b_["trk_pu_x"]->push_back(trkRef->vx());
                b_["trk_pu_y"]->push_back(trkRef->vy());
                b_["trk_pu_z"]->push_back(trkRef->vz());
                b_["trk_pu_pt"]->push_back(trkRef->pt());
                b_["trk_pu_pterr"]->push_back(trkRef->ptError());
                b_["trk_pu_eta"]->push_back(trkRef->eta());
                b_["trk_pu_etaerr"]->push_back(trkRef->etaError());
                b_["trk_pu_phi"]->push_back(trkRef->phi());
                b_["trk_pu_phierr"]->push_back(trkRef->phiError());
                b_["trk_pu_dxy"]->push_back(trkRef->dxy());
                b_["trk_pu_dxyerr"]->push_back(trkRef->dxyError());
                b_["trk_pu_dxysig"]->push_back(trkRef->dxy() / trkRef->dxyError());
                b_["trk_pu_dz"]->push_back(trkRef->dz());
                b_["trk_pu_dzerr"]->push_back(trkRef->dzError());
                b_["trk_pu_dzsig"]->push_back(trkRef->dz() / trkRef->dzError());
                b_["trk_pu_d3d"]->push_back(d3d);
                b_["trk_pu_d3derr"]->push_back(d3derr);
                b_["trk_pu_d3dsig"]->push_back(d3dsig);
                b_["trk_pu_d0"]->push_back(trkRef->d0());
                b_["trk_pu_d0err"]->push_back(trkRef->d0Error());
                b_["trk_pu_d0sig"]->push_back(trkRef->d0() / trkRef->d0Error());
                b_["trk_pu_charge"]->push_back(trkRef->charge());
                b_["trk_pu_chi2"]->push_back(trkRef->chi2());
                b_["trk_pu_ndof"]->push_back(trkRef->ndof());
                b_["trk_pu_chi2dof"]->push_back(trkRef->normalizedChi2());
            } // End pileup
        } // End gen matched
    } // End all reco tracks

    unsigned int evt_nGV = 0;
    for (unsigned int iGV = 0; iGV < genVertices.size(); iGV++) {
        const GenVertex* gv = genVertices.at(iGV);

        // Skip vertices not passing cuts
        if (!goodGenVertex(*gv, genPartPtCut_, genDauPtCut_, absEtaMax_)) continue;

        evt_nGV++;
        for (const reco::Candidate* dau : gv->daughters()) {
            if (!goodGenParticle(*dau, genDauPtCut_, absEtaMax_)) continue;

            b_["trk_gv_x"]->push_back(dau->vx());
            b_["trk_gv_y"]->push_back(dau->vy());
            b_["trk_gv_z"]->push_back(dau->vz());
            b_["trk_gv_pt"]->push_back(dau->pt());
            b_["trk_gv_eta"]->push_back(dau->eta());
            b_["trk_gv_phi"]->push_back(dau->phi());
            b_["trk_gv_charge"]->push_back(dau->charge());
            b_["trk_gv_pdgId"]->push_back(dau->pdgId());
        }

        b_["vtx_gv_x"]->push_back(gv->x());
        b_["vtx_gv_y"]->push_back(gv->y());
        b_["vtx_gv_z"]->push_back(gv->z());
        b_["vtx_gv_pt"]->push_back(gv->pt());
        b_["vtx_gv_eta"]->push_back(gv->eta());
        b_["vtx_gv_phi"]->push_back(gv->phi());
        b_["vtx_gv_ntracks"]->push_back(gv->nDaughters());
        b_["vtx_gv_motherPdgId"]->push_back(gv->motherPdgId());
    }

    // Fill PrimaryVertex information
    unsigned int evt_nPV = 0;
    for (unsigned int iPV = 0; iPV < PVs.size(); iPV++) {
        const reco::Vertex& prv = PVs.at(iPV);

        if (prv.isFake()) continue;
        // Skip vertices with less than 2 tracks passing cuts
        if (!goodRecoVertex(prv, timeValueMap, timeErrorMap, timeQualityMap, trackPtCut_, timeQualityCut_)) continue;

        evt_nPV++;
        unsigned int pv_ntrks = 0;
        for (reco::Vertex::trackRef_iterator trk_it = prv.tracks_begin(); trk_it != prv.tracks_end(); trk_it++) {
            reco::TrackBaseRef trkRef = *trk_it;

            if (!goodTrack(trkRef, timeValueMap, timeErrorMap, timeQualityMap, trackPtCut_, timeQualityCut_)) continue;

            pv_ntrks++;
            float tval = timeValueMap[trkRef];
            float terr = timeErrorMap[trkRef];
            float tsig = tval / terr;
            float tqual = timeQualityMap[trkRef];
            float d3d = TMath::Sqrt(trkRef->dxy()*trkRef->dxy() + trkRef->dz()*trkRef->dz());
            float dxy2err = 2*trkRef->dxyError()*trkRef->dxy();
            float dz2err = 2*trkRef->dzError()*trkRef->dz();
            float sum = d3d*d3d;
            float sumerr = TMath::Sqrt(dxy2err*dxy2err + dz2err*dz2err);
            float d3derr = 0.5*sumerr/TMath::Sqrt(sum);
            float d3dsig = d3d / d3derr;

            b_["trk_pv_tval"]->push_back(tval);
            b_["trk_pv_terr"]->push_back(terr);
            b_["trk_pv_tsig"]->push_back(tsig);
            b_["trk_pv_tqual"]->push_back(tqual);
            b_["trk_pv_x"]->push_back(trkRef->vx());
            b_["trk_pv_y"]->push_back(trkRef->vy());
            b_["trk_pv_z"]->push_back(trkRef->vz());
            b_["trk_pv_pt"]->push_back(trkRef->pt());
            b_["trk_pv_pterr"]->push_back(trkRef->ptError());
            b_["trk_pv_eta"]->push_back(trkRef->eta());
            b_["trk_pv_etaerr"]->push_back(trkRef->etaError());
            b_["trk_pv_phi"]->push_back(trkRef->phi());
            b_["trk_pv_phierr"]->push_back(trkRef->phiError());
            b_["trk_pv_dxy"]->push_back(trkRef->dxy());
            b_["trk_pv_dxyerr"]->push_back(trkRef->dxyError());
            b_["trk_pv_dxysig"]->push_back(trkRef->dxy() / trkRef->dxyError());
            b_["trk_pv_dz"]->push_back(trkRef->dz());
            b_["trk_pv_dzerr"]->push_back(trkRef->dzError());
            b_["trk_pv_dzsig"]->push_back(trkRef->dz() / trkRef->dzError());
            b_["trk_pv_d3d"]->push_back(d3d);
            b_["trk_pv_d3derr"]->push_back(d3derr);
            b_["trk_pv_d3dsig"]->push_back(d3dsig);
            b_["trk_pv_d0"]->push_back(trkRef->d0());
            b_["trk_pv_d0err"]->push_back(trkRef->d0Error());
            b_["trk_pv_d0sig"]->push_back(trkRef->d0() / trkRef->d0Error());
            b_["trk_pv_charge"]->push_back(trkRef->charge());
            b_["trk_pv_chi2"]->push_back(trkRef->chi2());
            b_["trk_pv_ndof"]->push_back(trkRef->ndof());
            b_["trk_pv_chi2dof"]->push_back(trkRef->normalizedChi2());
        } // End PV tracks

        b_["vtx_pv_ntracks"]->push_back(pv_ntrks);
    }

    unsigned int evt_nSV = 0;
    for (unsigned int iSV = 0; iSV < inclusiveSVs.size(); iSV++) {
        const reco::Vertex& sv = inclusiveSVs.at(iSV);

        if (!goodRecoVertex(sv, timeValueMap, timeErrorMap, timeQualityMap, trackPtCut_, timeQualityCut_)) continue;

        evt_nSV++;
        float timeMin = 99999.9;
        float timeMax = -99999.9;
        float timeAvg = 0.0;
        unsigned int sv_ntrks = 0;
        for (reco::Vertex::trackRef_iterator trk_it = sv.tracks_begin(); trk_it != sv.tracks_end(); trk_it++) {
            reco::TrackBaseRef trkRef = *trk_it;

            if (!goodTrack(trkRef, timeValueMap, timeErrorMap, timeQualityMap, trackPtCut_, timeQualityCut_)) continue;

            sv_ntrks++;
            float tval = timeValueMap[trkRef];
            float terr = timeErrorMap[trkRef];
            float tsig = tval / terr;
            float tqual = timeQualityMap[trkRef];
            float d3d = TMath::Sqrt(trkRef->dxy()*trkRef->dxy() + trkRef->dz()*trkRef->dz());
            float dxy2err = 2*trkRef->dxyError()*trkRef->dxy();
            float dz2err = 2*trkRef->dzError()*trkRef->dz();
            float sum = d3d*d3d;
            float sumerr = TMath::Sqrt(dxy2err*dxy2err + dz2err*dz2err);
            float d3derr = 0.5*sumerr/TMath::Sqrt(sum);
            float d3dsig = d3d / d3derr;

            timeMin = std::min(timeMin, tval);
            timeMax = std::max(timeMax, tval);
            timeAvg += tval;

            b_["trk_sv_tval"]->push_back(tval);
            b_["trk_sv_terr"]->push_back(terr);
            b_["trk_sv_tsig"]->push_back(tsig);
            b_["trk_sv_tqual"]->push_back(tqual);
            b_["trk_sv_x"]->push_back(trkRef->vx());
            b_["trk_sv_y"]->push_back(trkRef->vy());
            b_["trk_sv_z"]->push_back(trkRef->vz());
            b_["trk_sv_pt"]->push_back(trkRef->pt());
            b_["trk_sv_pterr"]->push_back(trkRef->ptError());
            b_["trk_sv_eta"]->push_back(trkRef->eta());
            b_["trk_sv_etaerr"]->push_back(trkRef->etaError());
            b_["trk_sv_phi"]->push_back(trkRef->phi());
            b_["trk_sv_phierr"]->push_back(trkRef->phiError());
            b_["trk_sv_dxy"]->push_back(trkRef->dxy());
            b_["trk_sv_dxyerr"]->push_back(trkRef->dxyError());
            b_["trk_sv_dxysig"]->push_back(trkRef->dxy() / trkRef->dxyError());
            b_["trk_sv_dz"]->push_back(trkRef->dz());
            b_["trk_sv_dzerr"]->push_back(trkRef->dzError());
            b_["trk_sv_dzsig"]->push_back(trkRef->dz() / trkRef->dzError());
            b_["trk_sv_d3d"]->push_back(d3d);
            b_["trk_sv_d3derr"]->push_back(d3derr);
            b_["trk_sv_d3dsig"]->push_back(d3dsig);
            b_["trk_sv_d0"]->push_back(trkRef->d0());
            b_["trk_sv_d0err"]->push_back(trkRef->d0Error());
            b_["trk_sv_d0sig"]->push_back(trkRef->d0() / trkRef->d0Error());
            b_["trk_sv_charge"]->push_back(trkRef->charge());
            b_["trk_sv_chi2"]->push_back(trkRef->chi2());
            b_["trk_sv_ndof"]->push_back(trkRef->ndof());
            b_["trk_sv_chi2dof"]->push_back(trkRef->normalizedChi2());
        }

        timeAvg /= sv_ntrks;
        float timeRange = timeMax - timeMin;

        b_["vtx_sv_timeavg"]->push_back(timeAvg);
        b_["vtx_sv_timerange"]->push_back(timeRange);
        b_["vtx_sv_x"]->push_back(sv.x());
        b_["vtx_sv_y"]->push_back(sv.y());
        b_["vtx_sv_z"]->push_back(sv.z());
        b_["vtx_sv_pt"]->push_back(vertexPt(inclusiveSVs.at(iSV)));
        b_["vtx_sv_eta"]->push_back(vertexEta(inclusiveSVs.at(iSV)));
        b_["vtx_sv_phi"]->push_back(vertexPhi(inclusiveSVs.at(iSV)));
        b_["vtx_sv_dxy"]->push_back(vertexDxy(sv, pv).value());
        b_["vtx_sv_dxyerr"]->push_back(catchInfsAndBound(vertexDxy(sv, pv).error() - 2, 0, -2, 0));
        b_["vtx_sv_dxysig"]->push_back(catchInfsAndBound(vertexDxy(sv, pv).value() / vertexDxy(sv, pv).error(), 0, -1, 800));
        b_["vtx_sv_dz"]->push_back(TMath::Abs(sv.z() - pv.z()));
        b_["vtx_sv_d3d"]->push_back(vertexD3d(sv, pv).value());
        b_["vtx_sv_d3derr"]->push_back(catchInfsAndBound(vertexD3d(sv, pv).error() - 2, 0, -2, 0));
        b_["vtx_sv_d3dsig"]->push_back(catchInfsAndBound(vertexD3d(sv, pv).value() / vertexD3d(sv, pv).error(), 0, -1, 800));
        b_["vtx_sv_chi2"]->push_back(sv.chi2());
        b_["vtx_sv_ndof"]->push_back(sv.ndof());
        b_["vtx_sv_chi2dof"]->push_back(catchInfsAndBound(sv.chi2() / sv.ndof(), 1000, -1000, 1000));
        b_["vtx_sv_ntracks"]->push_back(sv_ntrks);
    }

    unsigned int evt_nSVt = 0;
    for (unsigned int iSVt = 0; iSVt < inclusiveSVsMTDTiming.size(); iSVt++) {
        const reco::Vertex& sv = inclusiveSVsMTDTiming.at(iSVt);

        // Skip vertices with less than 2 tracks passing cuts
        if (!goodRecoVertex(sv, timeValueMap, timeErrorMap, timeQualityMap, trackPtCut_, timeQualityCut_)) continue;

        evt_nSVt++;
        float timeMin = 99999.9;
        float timeMax = -99999.9;
        float timeAvg = 0.0;
        unsigned int svt_ntrks = 0;
        for (reco::Vertex::trackRef_iterator trk_it = sv.tracks_begin(); trk_it != sv.tracks_end(); trk_it++) {
            reco::TrackBaseRef trkRef = *trk_it;

            if (!goodTrack(trkRef, timeValueMap, timeErrorMap, timeQualityMap, trackPtCut_, timeQualityCut_)) continue;

            svt_ntrks++;
            float tval = timeValueMap[trkRef];
            float terr = timeErrorMap[trkRef];
            float tsig = tval / terr;
            float tqual = timeQualityMap[trkRef];
            float d3d = TMath::Sqrt(trkRef->dxy()*trkRef->dxy() + trkRef->dz()*trkRef->dz());
            float dxy2err = 2*trkRef->dxyError()*trkRef->dxy();
            float dz2err = 2*trkRef->dzError()*trkRef->dz();
            float sum = d3d*d3d;
            float sumerr = TMath::Sqrt(dxy2err*dxy2err + dz2err*dz2err);
            float d3derr = 0.5*sumerr/TMath::Sqrt(sum);
            float d3dsig = d3d / d3derr;

            timeMin = std::min(timeMin, tval);
            timeMax = std::max(timeMax, tval);
            timeAvg += tval;

            b_["trk_svt_tval"]->push_back(tval);
            b_["trk_svt_terr"]->push_back(terr);
            b_["trk_svt_tsig"]->push_back(tsig);
            b_["trk_svt_tqual"]->push_back(tqual);
            b_["trk_svt_x"]->push_back(trkRef->vx());
            b_["trk_svt_y"]->push_back(trkRef->vy());
            b_["trk_svt_z"]->push_back(trkRef->vz());
            b_["trk_svt_pt"]->push_back(trkRef->pt());
            b_["trk_svt_pterr"]->push_back(trkRef->ptError());
            b_["trk_svt_eta"]->push_back(trkRef->eta());
            b_["trk_svt_etaerr"]->push_back(trkRef->etaError());
            b_["trk_svt_phi"]->push_back(trkRef->phi());
            b_["trk_svt_phierr"]->push_back(trkRef->phiError());
            b_["trk_svt_dxy"]->push_back(trkRef->dxy());
            b_["trk_svt_dxyerr"]->push_back(trkRef->dxyError());
            b_["trk_svt_dxysig"]->push_back(trkRef->dxy() / trkRef->dxyError());
            b_["trk_svt_dz"]->push_back(trkRef->dz());
            b_["trk_svt_dzerr"]->push_back(trkRef->dzError());
            b_["trk_svt_dzsig"]->push_back(trkRef->dz() / trkRef->dzError());
            b_["trk_svt_d3d"]->push_back(d3d);
            b_["trk_svt_d3derr"]->push_back(d3derr);
            b_["trk_svt_d3dsig"]->push_back(d3dsig);
            b_["trk_svt_d0"]->push_back(trkRef->d0());
            b_["trk_svt_d0err"]->push_back(trkRef->d0Error());
            b_["trk_svt_d0sig"]->push_back(trkRef->d0() / trkRef->d0Error());
            b_["trk_svt_charge"]->push_back(trkRef->charge());
            b_["trk_svt_chi2"]->push_back(trkRef->chi2());
            b_["trk_svt_ndof"]->push_back(trkRef->ndof());
            b_["trk_svt_chi2dof"]->push_back(trkRef->normalizedChi2());
        }

        timeAvg /= svt_ntrks;
        float timeRange = timeMax - timeMin;

        b_["vtx_svt_timeavg"]->push_back(timeAvg);
        b_["vtx_svt_timerange"]->push_back(timeRange);
        b_["vtx_svt_x"]->push_back(sv.x());
        b_["vtx_svt_y"]->push_back(sv.y());
        b_["vtx_svt_z"]->push_back(sv.z());
        b_["vtx_svt_pt"]->push_back(vertexPt(inclusiveSVs.at(iSVt)));
        b_["vtx_svt_eta"]->push_back(vertexEta(inclusiveSVs.at(iSVt)));
        b_["vtx_svt_phi"]->push_back(vertexPhi(inclusiveSVs.at(iSVt)));
        b_["vtx_svt_dxy"]->push_back(vertexDxy(sv, pv).value());
        b_["vtx_svt_dxyerr"]->push_back(catchInfsAndBound(vertexDxy(sv, pv).error() - 2, 0, -2, 0));
        b_["vtx_svt_dxysig"]->push_back(catchInfsAndBound(vertexDxy(sv, pv).value() / vertexDxy(sv, pv).error(), 0, -1, 800));
        b_["vtx_svt_dz"]->push_back(TMath::Abs(sv.z() - pv.z()));
        b_["vtx_svt_d3d"]->push_back(vertexD3d(sv, pv).value());
        b_["vtx_svt_d3derr"]->push_back(catchInfsAndBound(vertexD3d(sv, pv).error() - 2, 0, -2, 0));
        b_["vtx_svt_d3dsig"]->push_back(catchInfsAndBound(vertexD3d(sv, pv).value() / vertexD3d(sv, pv).error(), 0, -1, 800));
        b_["vtx_svt_chi2"]->push_back(sv.chi2());
        b_["vtx_svt_ndof"]->push_back(sv.ndof());
        b_["vtx_svt_chi2dof"]->push_back(catchInfsAndBound(sv.chi2() / sv.ndof(), 1000, -1000, 1000));
        b_["vtx_svt_ntracks"]->push_back(svt_ntrks);
    }

    // Jet stuff
    for (unsigned int iJet = 0; iJet < jetCollection.size(); iJet++) {
        const pat::Jet& jet = jetCollection.at(iJet);

        // Some cuts
        // if (jet.pt() < jetPtMin_) continue;
        // if (jet.pt() > jetPtMax_) continue;
        // if (TMath::Abs(jet.eta()) < absEtaMin_) continue;
        // if (TMath::Abs(jet.eta()) > absEtaMax_) continue;

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
        // int genJetHadFlav = -1000;
        // int genJetPartFlav = -1000;
        // const reco::GenJet* genJet = jet.genJet();
        // if (genJet) {
        //     for (const reco::JetFlavourInfoMatching& genJetFlavInfo : *(genJetFlavourInfo_.product())) {
        //         if (reco::deltaR(genJet->p4(), genJetFlavInfo.first->p4()) < genJetMatchdR_) {
        //             genJetHadFlav = genJetFlavInfo.second.getHadronFlavour();
        //             genJetPartFlav = genJetFlavInfo.second.getPartonFlavour();
        //             break;
        //         }
        //     }
        // }

        // jet_pt_->push_back(jet.pt());
        // jet_eta_->push_back(jet.eta());
        // jet_phi_->push_back(jet.phi());
        // jet_radius_->push_back(jet_radius);
        // jet_hadFlav_->push_back(jet.hadronFlavour());
        // jet_partFlav_->push_back(jet.partonFlavour());
        // jet_genHadFlav_->push_back(genJetHadFlav);
        // jet_genPartFlav_->push_back(genJetPartFlav);

        // Match Jet to SV
        bool matchedToSV = false;
        bool matchedToGV = true;
        unsigned int nSV = 0;
        unsigned int nGV = 0;
        for (unsigned int iSV = 0; iSV < inclusiveSVs.size(); iSV++) {
            if (iSV >= max_sv_) break;
            const reco::Vertex& sv = inclusiveSVs.at(iSV);
            if (reco::deltaR(vertexEta(inclusiveSVs.at(iSV)), vertexPhi(inclusiveSVs.at(iSV)), jet.eta(), jet.phi()) < jet_radius) {
                matchedToSV = true;
                nSV++;

                // Match (jet-matched) SV to GV
                int closestGVIdx = -1;
                float mindR = 99999999.99; // some maximum value
                for (unsigned int iGV = 0; iGV < genVertices.size(); iGV++) {
                    const GenVertex* genVtx = genVertices.at(iGV);
                    float dR3D = deltaR3D(genVtx->x(), sv.x(), genVtx->y(), sv.y(), genVtx->z(), sv.z());
                    if ((dR3D < matchGVdR_) && (dR3D < mindR)) {
                        closestGVIdx = iGV;
                        mindR = dR3D;
                    }
                }
                if (closestGVIdx >= 0) {
                    matchedToGV = true;
                    nGV++;
                    // jet_matchSV_matchGV_SVdRtoGV_->push_back(mindR);
                }
            }
        }

        if (matchedToSV) {
            // jet_matchSV_pt_->push_back(jet.pt());
            // jet_matchSV_eta_->push_back(jet.eta());
            // jet_matchSV_phi_->push_back(jet.phi());
            // jet_matchSV_radius_->push_back(jet_radius);
            // jet_matchSV_nSV_->push_back(nSV);
            // jet_matchSV_hadFlav_->push_back(jet.hadronFlavour());
            // jet_matchSV_partFlav_->push_back(jet.partonFlavour());
            // jet_matchSV_genHadFlav_->push_back(genJetHadFlav);
            // jet_matchSV_genPartFlav_->push_back(genJetPartFlav);
        }
        if (matchedToSV && matchedToGV) {
            // jet_matchSV_matchGV_pt_->push_back(jet.pt());
            // jet_matchSV_matchGV_eta_->push_back(jet.eta());
            // jet_matchSV_matchGV_phi_->push_back(jet.phi());
            // jet_matchSV_matchGV_radius_->push_back(jet_radius);
            // jet_matchSV_matchGV_nGV_->push_back(nGV);
            // jet_matchSV_matchGV_hadFlav_->push_back(jet.hadronFlavour());
            // jet_matchSV_matchGV_partFlav_->push_back(jet.partonFlavour());
            // jet_matchSV_matchGV_genHadFlav_->push_back(genJetHadFlav);
            // jet_matchSV_matchGV_genPartFlav_->push_back(genJetPartFlav);
        }
    }

    b_["evt_nPV"]->push_back(evt_nPV);
    b_["evt_nSV"]->push_back(evt_nSV);
    b_["evt_nSVt"]->push_back(evt_nSVt);

    sv_num_ = 0;
    for (const reco::VertexCompositePtrCandidate& sv : cpvtx) {
        if (sv_num_ < (int) max_sv_) { // Limit number of SVs

            // sv_pt_[sv_num_] = sv.pt();
            // sv_eta_[sv_num_] = sv.eta();
            // sv_phi_[sv_num_] = sv.phi();
            sv_mass_[sv_num_] = sv.mass();
            sv_e_[sv_num_] = sv.energy();
            // sv_ntracks_[sv_num_] = sv.numberOfDaughters();
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
                // sv_time_->push_back(vertex_time);
            } // End loop through jets
        } // End if sv_num_ < max_sv_
        sv_num_++;
    } // End loop through SVs
    // nsv_ = sv_num_;

    return 0;
}


// Helper functions


template <class P>
bool ntuple_SV::goodGenParticle(const P& gp, float ptCut, float etaCut) {

    bool pass = true;
    if (gp.pt() < ptCut) pass = false;
    if (abs(gp.eta()) > etaCut) pass = false;
    return pass;
}


bool ntuple_SV::goodGenVertex(const GenVertex& gv,
        float motherPtCut, float dauPtCut, float etaCut) {

    bool pass = true;
    if (!goodGenParticle(*(gv.mother()), motherPtCut, etaCut)) pass = false;
    int nDaughters = gv.nDaughters();
    for (const reco::Candidate* dau : gv.daughters()) {
        if (!goodGenParticle(*dau, dauPtCut, etaCut)) nDaughters--;
    }
    if (nDaughters < 2) pass = false;
    return pass;
}


template <class T>
bool ntuple_SV::goodTrack(const T& trkRef, const edm::ValueMap<float>& timeValueMap,
        const edm::ValueMap<float>& timeErrorMap, const edm::ValueMap<float>& timeQualityMap,
        float trackPtCut, float timeQualityCut) {

    bool trkPass = true;
    if (trkRef->pt() < trackPtCut) trkPass = false;
    if (!timeValueMap.contains(trkRef.id())) trkPass = false;
    if (timeValueMap.contains(trkRef.id()) && timeErrorMap[trkRef] == -1.0) trkPass = false;
    // if (timeValueMap.contains(trkRef.id()) && timeQualityMap[trkRef] < timeQualityCut) trkPass = false;
    return trkPass;
}


bool ntuple_SV::goodRecoVertex(const reco::Vertex& rv, const edm::ValueMap<float>& timeValueMap,
        const edm::ValueMap<float>& timeErrorMap, const edm::ValueMap<float>& timeQualityMap,
        float trackPtCut, float timeQualityCut) {

    bool pass = true;
    int ntracks = rv.tracksSize();
    for (reco::Vertex::trackRef_iterator trk_it = rv.tracks_begin(); trk_it != rv.tracks_end(); trk_it++) {
        reco::TrackBaseRef trkRef = *trk_it;
        if (!goodTrack(trkRef, timeValueMap, timeErrorMap, timeQualityMap, trackPtCut, timeQualityCut)) ntracks--;
    }
    if (ntracks < 2) pass = false;
    return pass;
}


float ntuple_SV::vertexPt(const reco::Vertex& sv) {

    float pt = 0.0;
    for (reco::Vertex::trackRef_iterator trk_it = sv.tracks_begin(); trk_it != sv.tracks_end(); trk_it++) {
        reco::TrackBaseRef trkRef = *trk_it;
        pt += trkRef->pt();
    }
    return pt;
}


float ntuple_SV::vertexEta(const reco::Vertex& sv) {

    float eta = 0.0;
    for (reco::Vertex::trackRef_iterator trk_it = sv.tracks_begin(); trk_it != sv.tracks_end(); trk_it++) {
        reco::TrackBaseRef trkRef = *trk_it;
        eta += trkRef->eta();
    }
    return eta;
}


float ntuple_SV::vertexPhi(const reco::Vertex& sv) {

    float phi = 0.0;
    for (reco::Vertex::trackRef_iterator trk_it = sv.tracks_begin(); trk_it != sv.tracks_end(); trk_it++) {
        reco::TrackBaseRef trkRef = *trk_it;
        phi += trkRef->phi();
    }
    return phi;
}


Measurement1D ntuple_SV::vertexDxy(const reco::VertexCompositePtrCandidate& svcand, const reco::Vertex& pv) {

    VertexDistanceXY dist;
    reco::Vertex::CovarianceMatrix csv;
    svcand.fillVertexCovariance(csv);
    reco::Vertex svtx(svcand.vertex(), csv);
    return dist.distance(svtx, pv);
}


Measurement1D ntuple_SV::vertexDxy(const reco::Vertex& sv, const reco::Vertex& pv) {

    VertexDistanceXY dist;
    reco::Vertex::CovarianceMatrix csv;
    sv.fill(csv);
    reco::Vertex svtx(sv.position(), csv);
    return dist.distance(svtx, pv);
}


Measurement1D ntuple_SV::vertexD3d(const reco::VertexCompositePtrCandidate& svcand, const reco::Vertex& pv) {

    VertexDistance3D dist;
    reco::Vertex::CovarianceMatrix csv;
    svcand.fillVertexCovariance(csv);
    reco::Vertex svtx(svcand.vertex(), csv);
    return dist.distance(svtx, pv);
}


Measurement1D ntuple_SV::vertexD3d(const reco::Vertex& sv, const reco::Vertex& pv) {

    VertexDistance3D dist;
    reco::Vertex::CovarianceMatrix csv;
    sv.fill(csv);
    reco::Vertex svtx(sv.position(), csv);
    return dist.distance(svtx, pv);
}


bool ntuple_SV::candCompareDxyDxyErr(const reco::VertexCompositePtrCandidate& sva, const reco::VertexCompositePtrCandidate& svb) {

    reco::Vertex pv = *spvp_;
    float adxy = ntuple_SV::vertexDxy(sva, pv).value();
    float bdxy = ntuple_SV::vertexDxy(svb, pv).value();
    float aerr = ntuple_SV::vertexDxy(sva, pv).error();
    float berr = ntuple_SV::vertexDxy(svb, pv).error();

    float asig = ntuple_SV::catchInfs(adxy / aerr, 0.0);
    float bsig = ntuple_SV::catchInfs(bdxy / berr, 0.0);

    return bsig < asig;
}


bool ntuple_SV::vertexCompareDxyDxyErr(const reco::Vertex& sva, const reco::Vertex& svb) {

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
