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

        GenVertex(const reco::GenParticle* mother, std::vector<const reco::Candidate*>* daughters) : 
                mother_(mother), daughters_(daughters) {}

        const float x() const { return daughters_->at(0)->vx(); }
        const float y() const { return daughters_->at(0)->vy(); }
        const float z() const { return daughters_->at(0)->vz(); }
        const float pt() const { return mother_->pt(); }
        const float eta() const { return mother_->eta(); }
        const float phi() const { return mother_->phi(); }
        const unsigned int nDaughters() const { return daughters_->size(); }
        const int motherPdgId() const { return mother_->pdgId(); }

        const reco::GenParticle* mother() const { return mother_; }
        const std::vector<const reco::Candidate*>* daughters() const { return daughters_; }

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
                std::cout << daughters_->at(iDau)->pdgId();
                if (iDau < nDaughters() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
            std::cout << "    daughter pts   = ";
            for (unsigned int iDau = 0; iDau < nDaughters(); iDau++) {
                std::cout << daughters_->at(iDau)->pt();
                if (iDau < nDaughters() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
            std::cout << "    daughter etas   = ";
            for (unsigned int iDau = 0; iDau < nDaughters(); iDau++) {
                std::cout << daughters_->at(iDau)->eta();
                if (iDau < nDaughters() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
            std::cout << "    daughter phis   = ";
            for (unsigned int iDau = 0; iDau < nDaughters(); iDau++) {
                std::cout << daughters_->at(iDau)->phi();
                if (iDau < nDaughters() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
        }

    private:

        const reco::GenParticle* mother_;
        std::vector<const reco::Candidate*>* daughters_;
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
    genToSimTrackdR_ = iConfig.getParameter<double>("genToSimTrackdR");
    genToSimTrackPt_ = iConfig.getParameter<double>("genToSimTrackPt");
    trackPtCut_ = iConfig.getParameter<double>("trackPtCut");
    timeQualityCut_ = iConfig.getParameter<double>("timeQualityCut");
    recoToGenTrackdR_ = iConfig.getParameter<double>("recoToGenTrackdR");
    recoToGenTrackPt_ = iConfig.getParameter<double>("recoToGenTrackPt");
    recoToSimTrackdR_ = iConfig.getParameter<double>("recoToSimTrackdR");
    recoToSimTrackPt_ = iConfig.getParameter<double>("recoToSimTrackPt");
    GVSimTrackMatchFrac_ = iConfig.getParameter<double>("GVSimTrackMatchFrac");
    SVtoGVdR_ = iConfig.getParameter<double>("SVtoGVdR");
    SVGenPartMatchFrac_ = iConfig.getParameter<double>("SVGenPartMatchFrac");
    jetPtMin_ = iConfig.getParameter<double>("jetPtMin");
    jetPtMax_ = iConfig.getParameter<double>("jetPtMax");
    genJetMatchdR_ = iConfig.getParameter<double>("genJetMatchdR");

    debug_ = iConfig.getParameter<bool>("debug");
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
    iEvent.getByToken(simTracks_token_, simTracks_);
    // iEvent.getByToken(TPs_token_, TPs_);
    // iEvent.getByToken(pf_cand_token_, pf_cand_);
    // iEvent.getByToken(pf_mcmatch_token_, pf_mcmatch_);
    iEvent.getByToken(recoTracks_token_, recoTracks_);
    iEvent.getByToken(timeValueMap_token_, timeValueMap_);
    iEvent.getByToken(timeErrorMap_token_, timeErrorMap_);
    iEvent.getByToken(timeQualityMap_token_, timeQualityMap_);
    iEvent.getByToken(trackMCMatch_token_, trackMCMatch_);
    // iEvent.getByToken(TVs_token_, TVs_);
    iEvent.getByToken(PVs_token_, PVs_);
    iEvent.getByToken(inclusiveSVs_token_, inclusiveSVs_);
    iEvent.getByToken(IVFclusters_token_, IVFclusters_);
    iEvent.getByToken(inclusiveSVsMTDTiming_token_, inclusiveSVsMTDTiming_);
    iEvent.getByToken(IVFclustersMTDTiming_token_, IVFclustersMTDTiming_);
    iEvent.getByToken(genJetFlavourInfo_token_, genJetFlavourInfo_);
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
    const reco::Vertex& spv = vertices()->at(0); // Most likely the signal vertex
    spvp_ = &vertices()->at(0);

    edm::SimTrackContainer simTracks = *(simTracks_.product());
    reco::TrackCollection recoTracks = *(recoTracks_.product());
    const edm::ValueMap<float>& timeValueMap = *(timeValueMap_.product());
    const edm::ValueMap<float>& timeErrorMap = *(timeErrorMap_.product());
    const edm::ValueMap<float>& timeQualityMap = *(timeQualityMap_.product());
    const edm::Association<reco::GenParticleCollection> trackMCMatch = *(trackMCMatch_.product());
    std::vector<GenVertex> GVs;
    std::vector<GenVertex> simGVs;
    reco::VertexCollection PVs = *(PVs_.product()); // Not the slimmed collection
    reco::VertexCollection inclusiveSVs = *(inclusiveSVs_.product());
    unsigned int nClusters = *(IVFclusters_.product());
    unsigned int nClusterst = *(IVFclustersMTDTiming_.product());
    reco::VertexCollection inclusiveSVsMTDTiming = *(inclusiveSVsMTDTiming_.product());
    reco::VertexCompositePtrCandidateCollection cpvtx = *secVertices();
    const edm::View<pat::Jet> jetCollection = *jets();

    // Construct GenVertex collection ("good" GVs only!!)
    for (unsigned int iGP = 0; iGP < genParticles_->size(); iGP++) {
        const reco::GenParticle* gp = (genParticles_->at(iGP)).clone();
        int motherPartID = genPartID(gp->pdgId());

        if (motherPartID < 0) continue; // Mother is not interesting hadron
        if (!goodGenParticle(gp, genPartPtCut_, absEtaMax_)) continue; // Cut on mother
        if (gp->numberOfDaughters() < 2) continue; // Not a vertex

        bool lastInstance = true; // Check for last instance of interesting hadron
        for (unsigned int iDau = 0; iDau < gp->numberOfDaughters(); iDau++) {
            if (genPartID((gp->daughter(iDau))->pdgId()) == motherPartID) { // Cut on daughters
                lastInstance = false; // Not last instance of interesting hadron
                break;
            }
        }

        if (lastInstance) {
            std::vector<const reco::Candidate*>* goodDaughters = new std::vector<const reco::Candidate*>;
            for (unsigned int iDau = 0; iDau < gp->numberOfDaughters(); iDau++) {
                const reco::Candidate* dau = gp->daughter(iDau);
                if (goodGenParticle(dau, genDauPtCut_, absEtaMax_))
                    goodDaughters->push_back(dau);
            }
            if (goodDaughters->size() < 2) continue;

            GenVertex newGV(gp, goodDaughters);
            GVs.push_back(newGV);
            if (matchGenToSimVertex(newGV, simTracks, GVSimTrackMatchFrac_, genToSimTrackdR_, genToSimTrackPt_))
                simGVs.push_back(newGV);
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

    // Sort vertices by dxy significance to PV
    std::sort(GVs.begin(), GVs.end(), ntuple_SV::gvCompareDxyDxySig);
    std::sort(simGVs.begin(), simGVs.end(), ntuple_SV::gvCompareDxyDxySig);
    std::sort(inclusiveSVs.begin(), inclusiveSVs.end(), ntuple_SV::svCompareDxyDxySig);
    std::sort(inclusiveSVsMTDTiming.begin(), inclusiveSVsMTDTiming.end(), ntuple_SV::svCompareDxyDxySig);
    std::sort(cpvtx.begin(), cpvtx.end(), ntuple_SV::candCompareDxyDxySig);

    // All matching
    std::vector<int> GV_matchtoSV(GVs.size(), -1); // Save matched indices to other collection
    std::vector<int> SV_matchtoGV(inclusiveSVs.size(), -1);
    std::vector<int> GV_matchtoSVt(GVs.size(), -1);
    std::vector<int> SVt_matchtoGV(inclusiveSVsMTDTiming.size(), -1);
    std::vector<int> simGV_matchtoSV(simGVs.size(), -1);
    std::vector<int> SV_matchtoSimGV(inclusiveSVs.size(), -1);
    std::vector<int> simGV_matchtoSVt(simGVs.size(), -1);
    std::vector<int> SVt_matchtoSimGV(inclusiveSVsMTDTiming.size(), -1);
    std::vector<int> jet_matchGen(jetCollection.size(), -1000); // Save gen jet flavours
    std::vector<int> jet_matchtoSV(jetCollection.size(), 0); // Save number of matches
    std::vector<int> jet_matchtoSVGV(jetCollection.size(), 0);
    std::vector<int> jet_matchtoSVt(jetCollection.size(), 0);
    std::vector<int> jet_matchtoSVtGV(jetCollection.size(), 0);

    if (debug_) std::cout << "Matching vertices" << std::endl;
    for (unsigned int iGV = 0; iGV < GVs.size(); iGV++) {
        const GenVertex& gv = GVs.at(iGV);

        int matchSVIdx = -1;
        float mindR = SVtoGVdR_;
        for (unsigned int iSV = 0; iSV < inclusiveSVs.size(); iSV++) {
            const reco::Vertex& sv = inclusiveSVs.at(iSV);
            if (!goodRecoVertex(sv, timeValueMap, timeErrorMap, timeQualityMap, trackPtCut_, absEtaMax_, timeQualityCut_)) continue;
            float dR3D = vertexD3d(gv, sv);
            if (dR3D < mindR) {
                if (!matchRecoToGenVertex(sv, gv, SVGenPartMatchFrac_, recoToGenTrackdR_, recoToGenTrackPt_,
                        timeValueMap, timeErrorMap, timeQualityMap, trackPtCut_, absEtaMax_, timeQualityCut_)) continue;
                matchSVIdx = iSV;
                mindR = dR3D;
            }
        }
        if (matchSVIdx >= 0) {
            GV_matchtoSV.at(iGV) = matchSVIdx;
            SV_matchtoGV.at(matchSVIdx) = iGV;
        }

        int matchSVtIdx = -1;
        mindR = SVtoGVdR_;
        for (unsigned int iSVt = 0; iSVt < inclusiveSVsMTDTiming.size(); iSVt++) {
            const reco::Vertex& svt = inclusiveSVsMTDTiming.at(iSVt);
            if (!goodRecoVertex(svt, timeValueMap, timeErrorMap, timeQualityMap, trackPtCut_, absEtaMax_, timeQualityCut_)) continue;
            float dR3D = vertexD3d(gv, svt);
            if (dR3D < mindR) {
                if (!matchRecoToGenVertex(svt, gv, SVGenPartMatchFrac_, recoToGenTrackdR_, recoToGenTrackPt_,
                        timeValueMap, timeErrorMap, timeQualityMap, trackPtCut_, absEtaMax_, timeQualityCut_)) continue;
                matchSVtIdx = iSVt;
                mindR = dR3D;
            }
        }
        if (matchSVtIdx >= 0) {
            GV_matchtoSVt.at(iGV) = matchSVtIdx;
            SVt_matchtoGV.at(matchSVtIdx) = iGV;
        }
    }

    if (debug_) std::cout << "Matching vertices w/sim" << std::endl;
    for (unsigned int iGVs = 0; iGVs < simGVs.size(); iGVs++) {
        const GenVertex& gvs = simGVs.at(iGVs);

        int matchSVIdx = -1;
        float mindR = SVtoGVdR_;
        for (unsigned int iSV = 0; iSV < inclusiveSVs.size(); iSV++) {
            const reco::Vertex& sv = inclusiveSVs.at(iSV);
            if (!goodRecoVertex(sv, timeValueMap, timeErrorMap, timeQualityMap, trackPtCut_, absEtaMax_, timeQualityCut_)) continue;
            float dR3D = vertexD3d(gvs, sv);
            if (dR3D < mindR) {
                if (!matchRecoToGenVertex(sv, gvs, SVGenPartMatchFrac_, recoToGenTrackdR_, recoToGenTrackPt_,
                        timeValueMap, timeErrorMap, timeQualityMap, trackPtCut_, absEtaMax_, timeQualityCut_)) continue;
                matchSVIdx = iSV;
                mindR = dR3D;
            }
        }
        if (matchSVIdx >= 0) {
            simGV_matchtoSV.at(iGVs) = matchSVIdx;
            SV_matchtoSimGV.at(matchSVIdx) = iGVs;
        }

        int matchSVtIdx = -1;
        mindR = SVtoGVdR_;
        for (unsigned int iSVt = 0; iSVt < inclusiveSVsMTDTiming.size(); iSVt++) {
            const reco::Vertex& svt = inclusiveSVsMTDTiming.at(iSVt);
            if (!goodRecoVertex(svt, timeValueMap, timeErrorMap, timeQualityMap, trackPtCut_, absEtaMax_, timeQualityCut_)) continue;
            float dR3D = vertexD3d(gvs, svt);
            if (dR3D < mindR) {
                if (!matchRecoToGenVertex(svt, gvs, SVGenPartMatchFrac_, recoToGenTrackdR_, recoToGenTrackPt_,
                        timeValueMap, timeErrorMap, timeQualityMap, trackPtCut_, absEtaMax_, timeQualityCut_)) continue;
                matchSVtIdx = iSVt;
                mindR = dR3D;
            }
        }
        if (matchSVtIdx >= 0) {
            simGV_matchtoSVt.at(iGVs) = matchSVtIdx;
            SVt_matchtoSimGV.at(matchSVtIdx) = iGVs;
        }
    }

    if (debug_) std::cout << "Matching jets" << std::endl;
    for (unsigned int iJet = 0; iJet < jetCollection.size(); iJet++) {
        const pat::Jet& jet = jetCollection.at(iJet);

        // Get Gen Jet flavours
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
        jet_matchGen.at(iJet) = genJetHadFlav;

        float jet_radius = jetRadius(jetR(), jetCollection.at(iJet));

        // Match to SV and GV
        for (unsigned int iSV = 0; iSV < inclusiveSVs.size(); iSV++) {
            const reco::Vertex& sv = inclusiveSVs.at(iSV);
            if (reco::deltaR(sv.p4().Eta(), sv.p4().Phi(), jet.eta(), jet.phi()) < jet_radius) {
                jet_matchtoSV.at(iJet)++;
                if (SV_matchtoGV.at(iSV) >= 0) jet_matchtoSVGV.at(iJet)++;
            }
        }

        // Match to SV w/MTD timing and GV
        for (unsigned int iSV = 0; iSV < inclusiveSVsMTDTiming.size(); iSV++) {
            const reco::Vertex& sv = inclusiveSVsMTDTiming.at(iSV);
            if (reco::deltaR(sv.p4().Eta(), sv.p4().Phi(), jet.eta(), jet.phi()) < jet_radius) {
                jet_matchtoSVt.at(iJet)++;
                if (SVt_matchtoGV.at(iSV) >= 0) jet_matchtoSVtGV.at(iJet)++;
            }
        }
    }

    if (debug_) std::cout << "Sim tracks" << std::endl;
    for (unsigned int iST = 0; iST < simTracks.size(); iST++) {
        SimTrack st = simTracks.at(iST);
        b_["trk_st_pt"]->push_back(st.momentum().Pt());
        b_["trk_st_eta"]->push_back(st.momentum().Eta());
        b_["trk_st_phi"]->push_back(st.momentum().Phi());
        b_["trk_st_charge"]->push_back(st.charge());
    }

    if (debug_) std::cout << "Reco tracks" << std::endl;
    for (unsigned int iRT = 0; iRT < recoTracks.size(); iRT++) {
        reco::TrackRef trkRef(recoTracks_, iRT);

        if (!goodTrack(trkRef, timeValueMap, timeErrorMap, timeQualityMap, trackPtCut_, absEtaMax_, timeQualityCut_)) continue;

        b_["trk_all_tval"]->push_back(timeValueMap[trkRef]);
        b_["trk_all_terr"]->push_back(timeErrorMap[trkRef]);
        b_["trk_all_tsig"]->push_back(timeValueMap[trkRef] / timeErrorMap[trkRef]);
        b_["trk_all_tqual"]->push_back(timeQualityMap[trkRef]);
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
        b_["trk_all_d3d"]->push_back(trackD3d(trkRef));
        b_["trk_all_d3derr"]->push_back(trackD3dErr(trkRef));
        b_["trk_all_d3dsig"]->push_back(trackD3d(trkRef) / trackD3dErr(trkRef));
        b_["trk_all_d0"]->push_back(trkRef->d0());
        b_["trk_all_d0err"]->push_back(trkRef->d0Error());
        b_["trk_all_d0sig"]->push_back(trkRef->d0() / trkRef->d0Error());
        b_["trk_all_charge"]->push_back(trkRef->charge());
        b_["trk_all_chi2"]->push_back(trkRef->chi2());
        b_["trk_all_ndof"]->push_back(trkRef->ndof());
        b_["trk_all_chi2dof"]->push_back(trkRef->normalizedChi2());

        reco::GenParticleRef trkMCRef = trackMCMatch[trkRef];
        const int iST = matchRecoToSimTrack(trkRef, simTracks, recoToSimTrackdR_, recoToSimTrackPt_);
        if (trkMCRef.id().isValid() && iST >= 0) {
            // Tracked matched to GenParticle
            b_["trk_match_gp_tval"]->push_back(timeValueMap[trkRef]);
            b_["trk_match_gp_terr"]->push_back(timeErrorMap[trkRef]);
            b_["trk_match_gp_tsig"]->push_back(timeValueMap[trkRef] / timeErrorMap[trkRef]);
            b_["trk_match_gp_tqual"]->push_back(timeQualityMap[trkRef]);
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
            b_["trk_match_gp_d3d"]->push_back(trackD3d(trkRef));
            b_["trk_match_gp_d3derr"]->push_back(trackD3dErr(trkRef));
            b_["trk_match_gp_d3dsig"]->push_back(trackD3d(trkRef) / trackD3dErr(trkRef));
            b_["trk_match_gp_d0"]->push_back(trkRef->d0());
            b_["trk_match_gp_d0err"]->push_back(trkRef->d0Error());
            b_["trk_match_gp_d0sig"]->push_back(trkRef->d0() / trkRef->d0Error());
            b_["trk_match_gp_charge"]->push_back(trkRef->charge());
            b_["trk_match_gp_chi2"]->push_back(trkRef->chi2());
            b_["trk_match_gp_ndof"]->push_back(trkRef->ndof());
            b_["trk_match_gp_chi2dof"]->push_back(trkRef->normalizedChi2());

            // Check for track from pileup using SIM information
            if (isPileupTrack(simTracks.at(iST))) {
                b_["trk_pu_tval"]->push_back(timeValueMap[trkRef]);
                b_["trk_pu_terr"]->push_back(timeErrorMap[trkRef]);
                b_["trk_pu_tsig"]->push_back(timeValueMap[trkRef] / timeErrorMap[trkRef]);
                b_["trk_pu_tqual"]->push_back(timeQualityMap[trkRef]);
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
                b_["trk_pu_d3d"]->push_back(trackD3d(trkRef));
                b_["trk_pu_d3derr"]->push_back(trackD3dErr(trkRef));
                b_["trk_pu_d3dsig"]->push_back(trackD3d(trkRef) / trackD3dErr(trkRef));
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

    if (debug_) std::cout << "Gen vertices" << std::endl;
    unsigned int evt_nGV = 0;
    for (unsigned int iGV = 0; iGV < GVs.size(); iGV++) {
        const GenVertex& gv = GVs.at(iGV);
        // By construction, all GVs are "good"

        for (const reco::Candidate* dau : *(gv.daughters())) {
            b_["trk_gv_x"]->push_back(dau->vx());
            b_["trk_gv_y"]->push_back(dau->vy());
            b_["trk_gv_z"]->push_back(dau->vz());
            b_["trk_gv_pt"]->push_back(dau->pt());
            b_["trk_gv_eta"]->push_back(dau->eta());
            b_["trk_gv_phi"]->push_back(dau->phi());
            b_["trk_gv_charge"]->push_back(dau->charge());
            b_["trk_gv_pdgId"]->push_back(dau->pdgId());
        }

        evt_nGV++;
        b_["vtx_gv_x"]->push_back(gv.x());
        b_["vtx_gv_y"]->push_back(gv.y());
        b_["vtx_gv_z"]->push_back(gv.z());
        b_["vtx_gv_pt"]->push_back(gv.pt());
        b_["vtx_gv_eta"]->push_back(gv.eta());
        b_["vtx_gv_phi"]->push_back(gv.phi());
        b_["vtx_gv_dxy"]->push_back(vertexDxy(gv, spv));
        b_["vtx_gv_dxyerr"]->push_back(vertexDxyErr(gv, spv));
        b_["vtx_gv_dxysig"]->push_back(vertexDxy(gv, spv) / vertexDxyErr(gv, spv));
        b_["vtx_gv_dz"]->push_back(TMath::Abs(spv.z() - gv.z()));
        b_["vtx_gv_d3d"]->push_back(vertexD3d(gv, spv));
        b_["vtx_gv_d3derr"]->push_back(vertexD3dErr(gv, spv));
        b_["vtx_gv_d3dsig"]->push_back(vertexD3d(gv, spv) / vertexD3dErr(gv, spv));
        b_["vtx_gv_ntracks"]->push_back(gv.nDaughters());
        b_["vtx_gv_motherPdgId"]->push_back(gv.motherPdgId());

        // Match to SV
        if (debug_) std::cout << "Gen vertices matched to secondary vertices" << std::endl;
        int iSV = GV_matchtoSV.at(iGV);
        if (iSV >= 0) {
            const reco::Vertex& sv = inclusiveSVs.at(iSV);
            b_["vtx_matchgv_sv_x"]->push_back(gv.x());
            b_["vtx_matchgv_sv_y"]->push_back(gv.y());
            b_["vtx_matchgv_sv_z"]->push_back(gv.z());
            b_["vtx_matchgv_sv_pt"]->push_back(gv.pt());
            b_["vtx_matchgv_sv_eta"]->push_back(gv.eta());
            b_["vtx_matchgv_sv_phi"]->push_back(gv.phi());
            b_["vtx_matchgv_sv_dxy"]->push_back(vertexDxy(gv, spv));
            b_["vtx_matchgv_sv_dxyerr"]->push_back(vertexDxyErr(gv, spv));
            b_["vtx_matchgv_sv_dxysig"]->push_back(vertexDxy(gv, spv) / vertexDxyErr(gv, spv));
            b_["vtx_matchgv_sv_dz"]->push_back(TMath::Abs(spv.z() - gv.z()));
            b_["vtx_matchgv_sv_d3d"]->push_back(vertexD3d(gv, spv));
            b_["vtx_matchgv_sv_d3derr"]->push_back(vertexD3dErr(gv, spv));
            b_["vtx_matchgv_sv_d3dsig"]->push_back(vertexD3d(gv, spv) / vertexD3dErr(gv, spv));
            b_["vtx_matchgv_sv_ntracks"]->push_back(gv.nDaughters());
            b_["vtx_matchgv_sv_motherPdgId"]->push_back(gv.motherPdgId());
            b_["vtx_matchgv_sv_matchdR"]->push_back(vertexD3d(gv, sv));
        }

        // Match to SV with MTD timing
        if (debug_) std::cout << "Gen vertices matched to secondary vertices w/timing" << std::endl;
        int iSVt = GV_matchtoSVt.at(iGV);
        if (iSVt >= 0) {
            const reco::Vertex& svt = inclusiveSVsMTDTiming.at(iSVt);
            b_["vtx_matchgv_svt_x"]->push_back(gv.x());
            b_["vtx_matchgv_svt_y"]->push_back(gv.y());
            b_["vtx_matchgv_svt_z"]->push_back(gv.z());
            b_["vtx_matchgv_svt_pt"]->push_back(gv.pt());
            b_["vtx_matchgv_svt_eta"]->push_back(gv.eta());
            b_["vtx_matchgv_svt_phi"]->push_back(gv.phi());
            b_["vtx_matchgv_svt_dxy"]->push_back(vertexDxy(gv, spv));
            b_["vtx_matchgv_svt_dxyerr"]->push_back(vertexDxyErr(gv, spv));
            b_["vtx_matchgv_svt_dxysig"]->push_back(vertexDxy(gv, spv) / vertexDxyErr(gv, spv));
            b_["vtx_matchgv_svt_dz"]->push_back(TMath::Abs(spv.z() - gv.z()));
            b_["vtx_matchgv_svt_d3d"]->push_back(vertexD3d(gv, spv));
            b_["vtx_matchgv_svt_d3derr"]->push_back(vertexD3dErr(gv, spv));
            b_["vtx_matchgv_svt_d3dsig"]->push_back(vertexD3d(gv, spv) / vertexD3dErr(gv, spv));
            b_["vtx_matchgv_svt_ntracks"]->push_back(gv.nDaughters());
            b_["vtx_matchgv_svt_motherPdgId"]->push_back(gv.motherPdgId());
            b_["vtx_matchgv_svt_matchdR"]->push_back(vertexD3d(gv, svt));
        }
    }

    if (debug_) std::cout << "Gen vertices w/Sim" << std::endl;
    unsigned int evt_nGVs = 0;
    for (unsigned int iGVs = 0; iGVs < simGVs.size(); iGVs++) {
        const GenVertex& gvs = simGVs.at(iGVs);
        // By construction, all GVs are "good"

        for (const reco::Candidate* dau : *(gvs.daughters())) {
            b_["trk_gvs_x"]->push_back(dau->vx());
            b_["trk_gvs_y"]->push_back(dau->vy());
            b_["trk_gvs_z"]->push_back(dau->vz());
            b_["trk_gvs_pt"]->push_back(dau->pt());
            b_["trk_gvs_eta"]->push_back(dau->eta());
            b_["trk_gvs_phi"]->push_back(dau->phi());
            b_["trk_gvs_charge"]->push_back(dau->charge());
            b_["trk_gvs_pdgId"]->push_back(dau->pdgId());
        }

        evt_nGVs++;
        b_["vtx_gvs_x"]->push_back(gvs.x());
        b_["vtx_gvs_y"]->push_back(gvs.y());
        b_["vtx_gvs_z"]->push_back(gvs.z());
        b_["vtx_gvs_pt"]->push_back(gvs.pt());
        b_["vtx_gvs_eta"]->push_back(gvs.eta());
        b_["vtx_gvs_phi"]->push_back(gvs.phi());
        b_["vtx_gvs_dxy"]->push_back(vertexDxy(gvs, spv));
        b_["vtx_gvs_dxyerr"]->push_back(vertexDxyErr(gvs, spv));
        b_["vtx_gvs_dxysig"]->push_back(vertexDxy(gvs, spv) / vertexDxyErr(gvs, spv));
        b_["vtx_gvs_dz"]->push_back(TMath::Abs(spv.z() - gvs.z()));
        b_["vtx_gvs_d3d"]->push_back(vertexD3d(gvs, spv));
        b_["vtx_gvs_d3derr"]->push_back(vertexD3dErr(gvs, spv));
        b_["vtx_gvs_d3dsig"]->push_back(vertexD3d(gvs, spv) / vertexD3dErr(gvs, spv));
        b_["vtx_gvs_ntracks"]->push_back(gvs.nDaughters());
        b_["vtx_gvs_motherPdgId"]->push_back(gvs.motherPdgId());

        // Match to SV
        if (debug_) std::cout << "Gen vertices w/sim matched to secondary vertices" << std::endl;
        int iSV = simGV_matchtoSV.at(iGVs);
        if (iSV >= 0) {
            const reco::Vertex& sv = inclusiveSVs.at(iSV);
            b_["vtx_matchgvs_sv_x"]->push_back(gvs.x());
            b_["vtx_matchgvs_sv_y"]->push_back(gvs.y());
            b_["vtx_matchgvs_sv_z"]->push_back(gvs.z());
            b_["vtx_matchgvs_sv_pt"]->push_back(gvs.pt());
            b_["vtx_matchgvs_sv_eta"]->push_back(gvs.eta());
            b_["vtx_matchgvs_sv_phi"]->push_back(gvs.phi());
            b_["vtx_matchgvs_sv_dxy"]->push_back(vertexDxy(gvs, spv));
            b_["vtx_matchgvs_sv_dxyerr"]->push_back(vertexDxyErr(gvs, spv));
            b_["vtx_matchgvs_sv_dxysig"]->push_back(vertexDxy(gvs, spv) / vertexDxyErr(gvs, spv));
            b_["vtx_matchgvs_sv_dz"]->push_back(TMath::Abs(spv.z() - gvs.z()));
            b_["vtx_matchgvs_sv_d3d"]->push_back(vertexD3d(gvs, spv));
            b_["vtx_matchgvs_sv_d3derr"]->push_back(vertexD3dErr(gvs, spv));
            b_["vtx_matchgvs_sv_d3dsig"]->push_back(vertexD3d(gvs, spv) / vertexD3dErr(gvs, spv));
            b_["vtx_matchgvs_sv_ntracks"]->push_back(gvs.nDaughters());
            b_["vtx_matchgvs_sv_motherPdgId"]->push_back(gvs.motherPdgId());
            b_["vtx_matchgvs_sv_matchdR"]->push_back(vertexD3d(gvs, sv));
        }

        // Match to SV with MTD timing
        if (debug_) std::cout << "Gen vertices w/sim matched to secondary vertices w/timing" << std::endl;
        int iSVt = simGV_matchtoSVt.at(iGVs);
        if (iSVt >= 0) {
            const reco::Vertex& svt = inclusiveSVsMTDTiming.at(iSVt);
            b_["vtx_matchgvs_svt_x"]->push_back(gvs.x());
            b_["vtx_matchgvs_svt_y"]->push_back(gvs.y());
            b_["vtx_matchgvs_svt_z"]->push_back(gvs.z());
            b_["vtx_matchgvs_svt_pt"]->push_back(gvs.pt());
            b_["vtx_matchgvs_svt_eta"]->push_back(gvs.eta());
            b_["vtx_matchgvs_svt_phi"]->push_back(gvs.phi());
            b_["vtx_matchgvs_svt_dxy"]->push_back(vertexDxy(gvs, spv));
            b_["vtx_matchgvs_svt_dxyerr"]->push_back(vertexDxyErr(gvs, spv));
            b_["vtx_matchgvs_svt_dxysig"]->push_back(vertexDxy(gvs, spv) / vertexDxyErr(gvs, spv));
            b_["vtx_matchgvs_svt_dz"]->push_back(TMath::Abs(spv.z() - gvs.z()));
            b_["vtx_matchgvs_svt_d3d"]->push_back(vertexD3d(gvs, spv));
            b_["vtx_matchgvs_svt_d3derr"]->push_back(vertexD3dErr(gvs, spv));
            b_["vtx_matchgvs_svt_d3dsig"]->push_back(vertexD3d(gvs, spv) / vertexD3dErr(gvs, spv));
            b_["vtx_matchgvs_svt_ntracks"]->push_back(gvs.nDaughters());
            b_["vtx_matchgvs_svt_motherPdgId"]->push_back(gvs.motherPdgId());
            b_["vtx_matchgvs_svt_matchdR"]->push_back(vertexD3d(gvs, svt));
        }
    }

    // Fill PrimaryVertex information
    if (debug_) std::cout << "Primary vertices" << std::endl;
    unsigned int evt_nPV = 0;
    for (unsigned int iPV = 0; iPV < PVs.size(); iPV++) {
        const reco::Vertex& pv = PVs.at(iPV);

        if (pv.isFake()) continue;
        // Skip vertices with less than 2 tracks passing cuts
        if (!goodRecoVertex(pv, timeValueMap, timeErrorMap, timeQualityMap, trackPtCut_, absEtaMax_, timeQualityCut_)) continue;

        evt_nPV++;
        unsigned int pv_ntrks = 0;
        for (reco::Vertex::trackRef_iterator trk_it = pv.tracks_begin(); trk_it != pv.tracks_end(); trk_it++) {
            reco::TrackBaseRef trkRef = *trk_it;

            if (!goodTrack(trkRef, timeValueMap, timeErrorMap, timeQualityMap, trackPtCut_, absEtaMax_, timeQualityCut_)) continue;

            pv_ntrks++;
            b_["trk_pv_tval"]->push_back(timeValueMap[trkRef]);
            b_["trk_pv_terr"]->push_back(timeErrorMap[trkRef]);
            b_["trk_pv_tsig"]->push_back(timeValueMap[trkRef] / timeErrorMap[trkRef]);
            b_["trk_pv_tqual"]->push_back(timeQualityMap[trkRef]);
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
            b_["trk_pv_d3d"]->push_back(trackD3d(trkRef));
            b_["trk_pv_d3derr"]->push_back(trackD3dErr(trkRef));
            b_["trk_pv_d3dsig"]->push_back(trackD3d(trkRef) / trackD3dErr(trkRef));
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

    if (debug_) std::cout << "Secondary vertices" << std::endl;
    unsigned int evt_nSV = 0;
    for (unsigned int iSV = 0; iSV < inclusiveSVs.size(); iSV++) {
        const reco::Vertex& sv = inclusiveSVs.at(iSV);

        if (!goodRecoVertex(sv, timeValueMap, timeErrorMap, timeQualityMap, trackPtCut_, absEtaMax_, timeQualityCut_)) continue;

        evt_nSV++;
        float timeMin = 99999.9;
        float timeMax = -99999.9;
        float timeAvg = 0.0;
        unsigned int sv_ntrks = 0;
        for (reco::Vertex::trackRef_iterator trk_it = sv.tracks_begin(); trk_it != sv.tracks_end(); trk_it++) {
            reco::TrackBaseRef trkRef = *trk_it;

            if (!goodTrack(trkRef, timeValueMap, timeErrorMap, timeQualityMap, trackPtCut_, absEtaMax_, timeQualityCut_)) continue;

            sv_ntrks++;
            timeMin = std::min(timeMin, timeValueMap[trkRef]);
            timeMax = std::max(timeMax, timeValueMap[trkRef]);
            timeAvg += timeValueMap[trkRef];
            b_["trk_sv_tval"]->push_back(timeValueMap[trkRef]);
            b_["trk_sv_terr"]->push_back(timeErrorMap[trkRef]);
            b_["trk_sv_tsig"]->push_back(timeValueMap[trkRef] / timeErrorMap[trkRef]);
            b_["trk_sv_tqual"]->push_back(timeQualityMap[trkRef]);
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
            b_["trk_sv_d3d"]->push_back(trackD3d(trkRef));
            b_["trk_sv_d3derr"]->push_back(trackD3dErr(trkRef));
            b_["trk_sv_d3dsig"]->push_back(trackD3d(trkRef) / trackD3dErr(trkRef));
            b_["trk_sv_d0"]->push_back(trkRef->d0());
            b_["trk_sv_d0err"]->push_back(trkRef->d0Error());
            b_["trk_sv_d0sig"]->push_back(trkRef->d0() / trkRef->d0Error());
            b_["trk_sv_charge"]->push_back(trkRef->charge());
            b_["trk_sv_chi2"]->push_back(trkRef->chi2());
            b_["trk_sv_ndof"]->push_back(trkRef->ndof());
            b_["trk_sv_chi2dof"]->push_back(trkRef->normalizedChi2());

            // Match to GV
            if (debug_) std::cout << "Secondary vertex tracks matched to gen vertices" << std::endl;
            int iGV = SV_matchtoGV.at(iSV);
            if (iGV >= 0) {
                b_["trk_matchsv_gv_tval"]->push_back(timeValueMap[trkRef]);
                b_["trk_matchsv_gv_terr"]->push_back(timeErrorMap[trkRef]);
                b_["trk_matchsv_gv_tsig"]->push_back(timeValueMap[trkRef] / timeErrorMap[trkRef]);
                b_["trk_matchsv_gv_tqual"]->push_back(timeQualityMap[trkRef]);
                b_["trk_matchsv_gv_x"]->push_back(trkRef->vx());
                b_["trk_matchsv_gv_y"]->push_back(trkRef->vy());
                b_["trk_matchsv_gv_z"]->push_back(trkRef->vz());
                b_["trk_matchsv_gv_pt"]->push_back(trkRef->pt());
                b_["trk_matchsv_gv_pterr"]->push_back(trkRef->ptError());
                b_["trk_matchsv_gv_eta"]->push_back(trkRef->eta());
                b_["trk_matchsv_gv_etaerr"]->push_back(trkRef->etaError());
                b_["trk_matchsv_gv_phi"]->push_back(trkRef->phi());
                b_["trk_matchsv_gv_phierr"]->push_back(trkRef->phiError());
                b_["trk_matchsv_gv_dxy"]->push_back(trkRef->dxy());
                b_["trk_matchsv_gv_dxyerr"]->push_back(trkRef->dxyError());
                b_["trk_matchsv_gv_dxysig"]->push_back(trkRef->dxy() / trkRef->dxyError());
                b_["trk_matchsv_gv_dz"]->push_back(trkRef->dz());
                b_["trk_matchsv_gv_dzerr"]->push_back(trkRef->dzError());
                b_["trk_matchsv_gv_dzsig"]->push_back(trkRef->dz() / trkRef->dzError());
                b_["trk_matchsv_gv_d3d"]->push_back(trackD3d(trkRef));
                b_["trk_matchsv_gv_d3derr"]->push_back(trackD3dErr(trkRef));
                b_["trk_matchsv_gv_d3dsig"]->push_back(trackD3d(trkRef) / trackD3dErr(trkRef));
                b_["trk_matchsv_gv_d0"]->push_back(trkRef->d0());
                b_["trk_matchsv_gv_d0err"]->push_back(trkRef->d0Error());
                b_["trk_matchsv_gv_d0sig"]->push_back(trkRef->d0() / trkRef->d0Error());
                b_["trk_matchsv_gv_charge"]->push_back(trkRef->charge());
                b_["trk_matchsv_gv_chi2"]->push_back(trkRef->chi2());
                b_["trk_matchsv_gv_ndof"]->push_back(trkRef->ndof());
                b_["trk_matchsv_gv_chi2dof"]->push_back(trkRef->normalizedChi2());
            }

            // Match to GV w/sim
            if (debug_) std::cout << "Secondary vertex tracks matched to GV w/sim" << std::endl;
            int iGVs = SV_matchtoSimGV.at(iSV);
            if (iGVs >= 0) {
                b_["trk_matchsv_gvs_tval"]->push_back(timeValueMap[trkRef]);
                b_["trk_matchsv_gvs_terr"]->push_back(timeErrorMap[trkRef]);
                b_["trk_matchsv_gvs_tsig"]->push_back(timeValueMap[trkRef] / timeErrorMap[trkRef]);
                b_["trk_matchsv_gvs_tqual"]->push_back(timeQualityMap[trkRef]);
                b_["trk_matchsv_gvs_x"]->push_back(trkRef->vx());
                b_["trk_matchsv_gvs_y"]->push_back(trkRef->vy());
                b_["trk_matchsv_gvs_z"]->push_back(trkRef->vz());
                b_["trk_matchsv_gvs_pt"]->push_back(trkRef->pt());
                b_["trk_matchsv_gvs_pterr"]->push_back(trkRef->ptError());
                b_["trk_matchsv_gvs_eta"]->push_back(trkRef->eta());
                b_["trk_matchsv_gvs_etaerr"]->push_back(trkRef->etaError());
                b_["trk_matchsv_gvs_phi"]->push_back(trkRef->phi());
                b_["trk_matchsv_gvs_phierr"]->push_back(trkRef->phiError());
                b_["trk_matchsv_gvs_dxy"]->push_back(trkRef->dxy());
                b_["trk_matchsv_gvs_dxyerr"]->push_back(trkRef->dxyError());
                b_["trk_matchsv_gvs_dxysig"]->push_back(trkRef->dxy() / trkRef->dxyError());
                b_["trk_matchsv_gvs_dz"]->push_back(trkRef->dz());
                b_["trk_matchsv_gvs_dzerr"]->push_back(trkRef->dzError());
                b_["trk_matchsv_gvs_dzsig"]->push_back(trkRef->dz() / trkRef->dzError());
                b_["trk_matchsv_gvs_d3d"]->push_back(trackD3d(trkRef));
                b_["trk_matchsv_gvs_d3derr"]->push_back(trackD3dErr(trkRef));
                b_["trk_matchsv_gvs_d3dsig"]->push_back(trackD3d(trkRef) / trackD3dErr(trkRef));
                b_["trk_matchsv_gvs_d0"]->push_back(trkRef->d0());
                b_["trk_matchsv_gvs_d0err"]->push_back(trkRef->d0Error());
                b_["trk_matchsv_gvs_d0sig"]->push_back(trkRef->d0() / trkRef->d0Error());
                b_["trk_matchsv_gvs_charge"]->push_back(trkRef->charge());
                b_["trk_matchsv_gvs_chi2"]->push_back(trkRef->chi2());
                b_["trk_matchsv_gvs_ndof"]->push_back(trkRef->ndof());
                b_["trk_matchsv_gvs_chi2dof"]->push_back(trkRef->normalizedChi2());
            }
        }

        timeAvg /= sv_ntrks;
        float timeRange = timeMax - timeMin;

        b_["vtx_sv_timeavg"]->push_back(timeAvg);
        b_["vtx_sv_timerange"]->push_back(timeRange);
        b_["vtx_sv_x"]->push_back(sv.x());
        b_["vtx_sv_y"]->push_back(sv.y());
        b_["vtx_sv_z"]->push_back(sv.z());
        b_["vtx_sv_pt"]->push_back(sv.p4().Pt());
        b_["vtx_sv_eta"]->push_back(sv.p4().Eta());
        b_["vtx_sv_phi"]->push_back(sv.p4().Phi());
        b_["vtx_sv_dxy"]->push_back(vertexDxy(sv, spv).value());
        b_["vtx_sv_dxyerr"]->push_back(vertexDxy(sv, spv).error());
        b_["vtx_sv_dxysig"]->push_back(vertexDxy(sv, spv).value() / vertexDxy(sv, spv).error());
        b_["vtx_sv_dz"]->push_back(TMath::Abs(sv.z() - spv.z()));
        b_["vtx_sv_d3d"]->push_back(vertexD3d(sv, spv).value());
        b_["vtx_sv_d3derr"]->push_back(vertexD3d(sv, spv).error() );
        b_["vtx_sv_d3dsig"]->push_back(vertexD3d(sv, spv).value() / vertexD3d(sv, spv).error());
        b_["vtx_sv_chi2"]->push_back(sv.chi2());
        b_["vtx_sv_ndof"]->push_back(sv.ndof());
        b_["vtx_sv_chi2dof"]->push_back(sv.chi2() / sv.ndof());
        b_["vtx_sv_ntracks"]->push_back(sv_ntrks);
    }

    if (debug_) std::cout << "Secondary vertices w/timing" << std::endl;
    unsigned int evt_nSVt = 0;
    for (unsigned int iSVt = 0; iSVt < inclusiveSVsMTDTiming.size(); iSVt++) {
        const reco::Vertex& sv = inclusiveSVsMTDTiming.at(iSVt);

        // Skip vertices with less than 2 tracks passing cuts
        if (!goodRecoVertex(sv, timeValueMap, timeErrorMap, timeQualityMap, trackPtCut_, absEtaMax_, timeQualityCut_)) continue;

        evt_nSVt++;
        float timeMin = 99999.9;
        float timeMax = -99999.9;
        float timeAvg = 0.0;
        unsigned int svt_ntrks = 0;
        for (reco::Vertex::trackRef_iterator trk_it = sv.tracks_begin(); trk_it != sv.tracks_end(); trk_it++) {
            reco::TrackBaseRef trkRef = *trk_it;

            if (!goodTrack(trkRef, timeValueMap, timeErrorMap, timeQualityMap, trackPtCut_, absEtaMax_, timeQualityCut_)) continue;

            svt_ntrks++;
            timeMin = std::min(timeMin, timeValueMap[trkRef]);
            timeMax = std::max(timeMax, timeValueMap[trkRef]);
            timeAvg += timeValueMap[trkRef];
            b_["trk_svt_tval"]->push_back(timeValueMap[trkRef]);
            b_["trk_svt_terr"]->push_back(timeErrorMap[trkRef]);
            b_["trk_svt_tsig"]->push_back(timeValueMap[trkRef] / timeErrorMap[trkRef]);
            b_["trk_svt_tqual"]->push_back(timeQualityMap[trkRef]);
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
            b_["trk_svt_d3d"]->push_back(trackD3d(trkRef));
            b_["trk_svt_d3derr"]->push_back(trackD3dErr(trkRef));
            b_["trk_svt_d3dsig"]->push_back(trackD3d(trkRef) / trackD3dErr(trkRef));
            b_["trk_svt_d0"]->push_back(trkRef->d0());
            b_["trk_svt_d0err"]->push_back(trkRef->d0Error());
            b_["trk_svt_d0sig"]->push_back(trkRef->d0() / trkRef->d0Error());
            b_["trk_svt_charge"]->push_back(trkRef->charge());
            b_["trk_svt_chi2"]->push_back(trkRef->chi2());
            b_["trk_svt_ndof"]->push_back(trkRef->ndof());
            b_["trk_svt_chi2dof"]->push_back(trkRef->normalizedChi2());

            // Match to GV
            if (debug_) std::cout << "Secondary vertex tracks w/timing matched to gen vertices" << std::endl;
            int iGV = SVt_matchtoGV.at(iSVt);
            if (iGV >= 0) {
                b_["trk_matchsvt_gv_tval"]->push_back(timeValueMap[trkRef]);
                b_["trk_matchsvt_gv_terr"]->push_back(timeErrorMap[trkRef]);
                b_["trk_matchsvt_gv_tsig"]->push_back(timeValueMap[trkRef] / timeErrorMap[trkRef]);
                b_["trk_matchsvt_gv_tqual"]->push_back(timeQualityMap[trkRef]);
                b_["trk_matchsvt_gv_x"]->push_back(trkRef->vx());
                b_["trk_matchsvt_gv_y"]->push_back(trkRef->vy());
                b_["trk_matchsvt_gv_z"]->push_back(trkRef->vz());
                b_["trk_matchsvt_gv_pt"]->push_back(trkRef->pt());
                b_["trk_matchsvt_gv_pterr"]->push_back(trkRef->ptError());
                b_["trk_matchsvt_gv_eta"]->push_back(trkRef->eta());
                b_["trk_matchsvt_gv_etaerr"]->push_back(trkRef->etaError());
                b_["trk_matchsvt_gv_phi"]->push_back(trkRef->phi());
                b_["trk_matchsvt_gv_phierr"]->push_back(trkRef->phiError());
                b_["trk_matchsvt_gv_dxy"]->push_back(trkRef->dxy());
                b_["trk_matchsvt_gv_dxyerr"]->push_back(trkRef->dxyError());
                b_["trk_matchsvt_gv_dxysig"]->push_back(trkRef->dxy() / trkRef->dxyError());
                b_["trk_matchsvt_gv_dz"]->push_back(trkRef->dz());
                b_["trk_matchsvt_gv_dzerr"]->push_back(trkRef->dzError());
                b_["trk_matchsvt_gv_dzsig"]->push_back(trkRef->dz() / trkRef->dzError());
                b_["trk_matchsvt_gv_d3d"]->push_back(trackD3d(trkRef));
                b_["trk_matchsvt_gv_d3derr"]->push_back(trackD3dErr(trkRef));
                b_["trk_matchsvt_gv_d3dsig"]->push_back(trackD3d(trkRef) / trackD3dErr(trkRef));
                b_["trk_matchsvt_gv_d0"]->push_back(trkRef->d0());
                b_["trk_matchsvt_gv_d0err"]->push_back(trkRef->d0Error());
                b_["trk_matchsvt_gv_d0sig"]->push_back(trkRef->d0() / trkRef->d0Error());
                b_["trk_matchsvt_gv_charge"]->push_back(trkRef->charge());
                b_["trk_matchsvt_gv_chi2"]->push_back(trkRef->chi2());
                b_["trk_matchsvt_gv_ndof"]->push_back(trkRef->ndof());
                b_["trk_matchsvt_gv_chi2dof"]->push_back(trkRef->normalizedChi2());
            }

            // Match to GV w/sim
            if (debug_) std::cout << "Secondary vertex tracks w/timing matched to GV w/sim" << std::endl;
            int iGVs = SVt_matchtoSimGV.at(iSVt);
            if (iGVs >= 0) {
                b_["trk_matchsvt_gvs_tval"]->push_back(timeValueMap[trkRef]);
                b_["trk_matchsvt_gvs_terr"]->push_back(timeErrorMap[trkRef]);
                b_["trk_matchsvt_gvs_tsig"]->push_back(timeValueMap[trkRef] / timeErrorMap[trkRef]);
                b_["trk_matchsvt_gvs_tqual"]->push_back(timeQualityMap[trkRef]);
                b_["trk_matchsvt_gvs_x"]->push_back(trkRef->vx());
                b_["trk_matchsvt_gvs_y"]->push_back(trkRef->vy());
                b_["trk_matchsvt_gvs_z"]->push_back(trkRef->vz());
                b_["trk_matchsvt_gvs_pt"]->push_back(trkRef->pt());
                b_["trk_matchsvt_gvs_pterr"]->push_back(trkRef->ptError());
                b_["trk_matchsvt_gvs_eta"]->push_back(trkRef->eta());
                b_["trk_matchsvt_gvs_etaerr"]->push_back(trkRef->etaError());
                b_["trk_matchsvt_gvs_phi"]->push_back(trkRef->phi());
                b_["trk_matchsvt_gvs_phierr"]->push_back(trkRef->phiError());
                b_["trk_matchsvt_gvs_dxy"]->push_back(trkRef->dxy());
                b_["trk_matchsvt_gvs_dxyerr"]->push_back(trkRef->dxyError());
                b_["trk_matchsvt_gvs_dxysig"]->push_back(trkRef->dxy() / trkRef->dxyError());
                b_["trk_matchsvt_gvs_dz"]->push_back(trkRef->dz());
                b_["trk_matchsvt_gvs_dzerr"]->push_back(trkRef->dzError());
                b_["trk_matchsvt_gvs_dzsig"]->push_back(trkRef->dz() / trkRef->dzError());
                b_["trk_matchsvt_gvs_d3d"]->push_back(trackD3d(trkRef));
                b_["trk_matchsvt_gvs_d3derr"]->push_back(trackD3dErr(trkRef));
                b_["trk_matchsvt_gvs_d3dsig"]->push_back(trackD3d(trkRef) / trackD3dErr(trkRef));
                b_["trk_matchsvt_gvs_d0"]->push_back(trkRef->d0());
                b_["trk_matchsvt_gvs_d0err"]->push_back(trkRef->d0Error());
                b_["trk_matchsvt_gvs_d0sig"]->push_back(trkRef->d0() / trkRef->d0Error());
                b_["trk_matchsvt_gvs_charge"]->push_back(trkRef->charge());
                b_["trk_matchsvt_gvs_chi2"]->push_back(trkRef->chi2());
                b_["trk_matchsvt_gvs_ndof"]->push_back(trkRef->ndof());
                b_["trk_matchsvt_gvs_chi2dof"]->push_back(trkRef->normalizedChi2());
            }
        }

        timeAvg /= svt_ntrks;
        float timeRange = timeMax - timeMin;

        b_["vtx_svt_timeavg"]->push_back(timeAvg);
        b_["vtx_svt_timerange"]->push_back(timeRange);
        b_["vtx_svt_x"]->push_back(sv.x());
        b_["vtx_svt_y"]->push_back(sv.y());
        b_["vtx_svt_z"]->push_back(sv.z());
        b_["vtx_svt_pt"]->push_back(sv.p4().Pt());
        b_["vtx_svt_eta"]->push_back(sv.p4().Eta());
        b_["vtx_svt_phi"]->push_back(sv.p4().Phi());
        b_["vtx_svt_dxy"]->push_back(vertexDxy(sv, spv).value());
        b_["vtx_svt_dxyerr"]->push_back(vertexDxy(sv, spv).error());
        b_["vtx_svt_dxysig"]->push_back(vertexDxy(sv, spv).value() / vertexDxy(sv, spv).error());
        b_["vtx_svt_dz"]->push_back(TMath::Abs(sv.z() - spv.z()));
        b_["vtx_svt_d3d"]->push_back(vertexD3d(sv, spv).value());
        b_["vtx_svt_d3derr"]->push_back(vertexD3d(sv, spv).error());
        b_["vtx_svt_d3dsig"]->push_back(vertexD3d(sv, spv).value() / vertexD3d(sv, spv).error());
        b_["vtx_svt_chi2"]->push_back(sv.chi2());
        b_["vtx_svt_ndof"]->push_back(sv.ndof());
        b_["vtx_svt_chi2dof"]->push_back(sv.chi2() / sv.ndof());
        b_["vtx_svt_ntracks"]->push_back(svt_ntrks);
    }

    if (debug_) std::cout << "Jets" << std::endl;
    for (unsigned int iJet = 0; iJet < jetCollection.size(); iJet++) {
        const pat::Jet& jet = jetCollection.at(iJet);

        if (!goodJet(jet, jetPtMin_, jetPtMax_, absEtaMax_)) continue;

        b_["jet_all_pt"]->push_back(jet.pt());
        b_["jet_all_eta"]->push_back(jet.eta());
        b_["jet_all_phi"]->push_back(jet.phi());

        // Match to Gen Jet
        if (debug_) std::cout << "Jets matched to gen jet" << std::endl;
        int hadFlav = jet_matchGen.at(iJet);
        if (hadFlav >= 0) {
            b_["jet_match_gen_pt"]->push_back(jet.pt());
            b_["jet_match_gen_eta"]->push_back(jet.eta());
            b_["jet_match_gen_phi"]->push_back(jet.phi());
            b_["jet_match_gen_hadflav"]->push_back(hadFlav);

            // Match to SV
            if (debug_) std::cout << "Jets matched to gen jet with SV" << std::endl;
            int nSV = jet_matchtoSV.at(iJet);
            if (nSV > 0) {
                b_["jet_match_gensv_pt"]->push_back(jet.pt());
                b_["jet_match_gensv_eta"]->push_back(jet.eta());
                b_["jet_match_gensv_phi"]->push_back(jet.phi());
                b_["jet_match_gensv_hadflav"]->push_back(hadFlav);
                b_["jet_match_gensv_nmatch"]->push_back(nSV);

                // Match to GV
                if (debug_) std::cout << "Jets matched to gen jet with SV matched to GV" << std::endl;
                int nSVGV = jet_matchtoSVGV.at(iJet);
                if (nSVGV > 0) {
                    b_["jet_match_gensvgv_pt"]->push_back(jet.pt());
                    b_["jet_match_gensvgv_eta"]->push_back(jet.eta());
                    b_["jet_match_gensvgv_phi"]->push_back(jet.phi());
                    b_["jet_match_gensvgv_hadflav"]->push_back(hadFlav);
                    b_["jet_match_gensvgv_nmatch"]->push_back(nSVGV);
                }
            }

            // Match to SV w/timing
            if (debug_) std::cout << "Jets matched to gen jet with SV w/timing" << std::endl;
            int nSVt = jet_matchtoSVt.at(iJet);
            if (nSVt > 0) {
                b_["jet_match_gensvt_pt"]->push_back(jet.pt());
                b_["jet_match_gensvt_eta"]->push_back(jet.eta());
                b_["jet_match_gensvt_phi"]->push_back(jet.phi());
                b_["jet_match_gensvt_hadflav"]->push_back(hadFlav);
                b_["jet_match_gensvt_nmatch"]->push_back(nSVt);

                // Match to GV
                if (debug_) std::cout << "Jets matched to gen jet with SV w/timing matched to GV" << std::endl;
                int nSVtGV = jet_matchtoSVtGV.at(iJet);
                if (nSVtGV > 0) {
                    b_["jet_match_gensvtgv_pt"]->push_back(jet.pt());
                    b_["jet_match_gensvtgv_eta"]->push_back(jet.eta());
                    b_["jet_match_gensvtgv_phi"]->push_back(jet.phi());
                    b_["jet_match_gensvtgv_hadflav"]->push_back(hadFlav);
                    b_["jet_match_gensvtgv_nmatch"]->push_back(nSVtGV);
                }
            }
        }
    }

    b_["evt_nGV"]->push_back(evt_nGV);
    b_["evt_nGVs"]->push_back(evt_nGVs);
    b_["evt_nPV"]->push_back(evt_nPV);
    b_["evt_nSV"]->push_back(evt_nSV);
    b_["evt_nSVt"]->push_back(evt_nSVt);
    b_["evt_nClusters"]->push_back(nClusters);
    b_["evt_nClusterst"]->push_back(nClusterst);

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
            // sv_dxy_[sv_num_] = vertexDxy(sv, spv).value();
            // sv_dxyerr_[sv_num_] = catchInfsAndBound(vertexDxy(sv, spv).error() - 2, 0, -2, 0);
            // sv_dxysig_[sv_num_] = catchInfsAndBound(sv_dxy_[sv_num_] / vertexDxy(sv, spv).error(), 0, -1, 800);
            // sv_d3d_[sv_num_] = vertexD3d(sv, spv).value();
            // sv_d3derr_[sv_num_] = catchInfsAndBound(vertexD3d(sv, spv).error() - 2, 0, -2, 0);
            // sv_d3dsig_[sv_num_] = catchInfsAndBound(vertexD3d(sv, spv).value() / vertexD3d(sv, spv).error(), 0, -1, 800);
            sv_costhetasvpv_[sv_num_] = vertexDdotP(sv, spv); // The pointing angle (i.e. the angle between the sum of the momentum of the tracks in the SV and the flight direction betwen PV and SV)

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
                        trackinfo.buildTrackInfo(PackedCandidate_, jetDir, jetRefTrackDir, spv);
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


// Quality cuts


template <class P>
bool ntuple_SV::goodGenParticle(const P* gp, float ptCut, float etaCut) {

    bool pass = true;
    if (gp->pt() < ptCut) pass = false;
    if (abs(gp->eta()) > etaCut) pass = false;
    return pass;
}


// bool ntuple_SV::goodGenVertex(const GenVertex& gv,
//         float motherPtCut, float dauPtCut, float etaCut) {

//     bool pass = true;
//     if (!goodGenParticle(gv.mother(), motherPtCut, etaCut)) pass = false;
//     int nDaughters = gv.nDaughters();
//     for (const reco::Candidate* dau : *(gv.daughters())) {
//         if (!goodGenParticle(dau, dauPtCut, etaCut)) nDaughters--;
//     }
//     if (nDaughters < 2) pass = false;
//     return pass;
// }


template <class T>
bool ntuple_SV::goodTrack(const T& trkRef, const edm::ValueMap<float>& timeValueMap,
        const edm::ValueMap<float>& timeErrorMap, const edm::ValueMap<float>& timeQualityMap,
        float trackPtCut, float trackEtaCut, float timeQualityCut) {

    bool trkPass = true;
    if (trkRef->pt() < trackPtCut) trkPass = false;
    if (abs(trkRef->eta()) > trackEtaCut) trkPass = false;
    if (!timeValueMap.contains(trkRef.id())) trkPass = false;
    if (timeValueMap.contains(trkRef.id()) && timeErrorMap[trkRef] == -1.0) trkPass = false;
    // if (timeValueMap.contains(trkRef.id()) && timeQualityMap[trkRef] < timeQualityCut) trkPass = false;
    return trkPass;
}


bool ntuple_SV::goodRecoVertex(const reco::Vertex& rv, const edm::ValueMap<float>& timeValueMap,
        const edm::ValueMap<float>& timeErrorMap, const edm::ValueMap<float>& timeQualityMap,
        float trackPtCut, float trackEtaCut, float timeQualityCut) {

    bool pass = true;
    int ntracks = rv.tracksSize();
    for (reco::Vertex::trackRef_iterator trk_it = rv.tracks_begin(); trk_it != rv.tracks_end(); trk_it++) {
        reco::TrackBaseRef trkRef = *trk_it;
        if (!goodTrack(trkRef, timeValueMap, timeErrorMap, timeQualityMap, trackPtCut, trackEtaCut, timeQualityCut)) ntracks--;
    }
    if (ntracks < 2) pass = false;
    return pass;
}


bool ntuple_SV::goodJet(const pat::Jet& jet, float ptMin, float ptMax, float etaCut) {

    bool pass = true;
    if (jet.pt() < ptMin) pass = false;
    if (jet.pt() > ptMax) pass = false;
    if (abs(jet.eta()) > etaCut) pass = false;
    return pass;
}


// Matching


bool ntuple_SV::isPileupTrack(const SimTrack& st) {

    return st.eventId().event() > 0;
}


template <class P>
bool ntuple_SV::matchGenToSimTrack(const P* gt, const SimTrack& st, float drCut, float ptCut) {

    bool match = true;
    if (reco::deltaR(gt->eta(), gt->phi(), st.momentum().Eta(), st.momentum().Phi()) > drCut) match = false;
    float dpt = abs(gt->pt() - st.momentum().Pt()) / (gt->pt() + st.momentum().Pt());
    if (dpt > ptCut) match = false;
    return match;
}


template <class T, class P>
bool ntuple_SV::matchRecoToGenTrack(const T& trkRef, const P* gt, float drCut, float ptCut) {

    bool match = true;
    if (reco::deltaR(trkRef->eta(), trkRef->phi(), gt->eta(), gt->phi()) > drCut) match = false;
    float dpt = abs(gt->pt() - trkRef->pt()) / (gt->pt() + trkRef->pt());
    if (dpt > ptCut) match = false;
    return match;
}


template <class T>
int ntuple_SV::matchRecoToSimTrack(const T& trkRef, const edm::SimTrackContainer& simTracks,
        float drCut, float ptCut) {

    int stIdx = -1;
    for (const SimTrack& st : simTracks) {
        stIdx++;
        if (reco::deltaR(trkRef->eta(), trkRef->phi(), st.momentum().Eta(), st.momentum().Phi()) > drCut) continue;
        float dpt = abs(trkRef->pt() - st.momentum().Pt()) / (trkRef->pt() + st.momentum().Pt());
        if (dpt > ptCut) continue;
        return stIdx;
    }
    return -1;
}


bool ntuple_SV::matchGenToSimVertex(const GenVertex& gv, const edm::SimTrackContainer& simTracks,
        float matchFrac, float trkDrCut, float trkPtCut) {

    float nmatch = 0;
    for (const reco::Candidate* dau : *(gv.daughters())) {
        for (const SimTrack& st : simTracks) {
            if (matchGenToSimTrack(dau, st, trkDrCut, trkPtCut)) {
                nmatch++;
                break;
            }
        }
    }
    return (nmatch/(float)gv.nDaughters()) >= matchFrac;
}


bool ntuple_SV::matchRecoToGenVertex(const reco::Vertex& v, const GenVertex& gv,
        float matchFrac, float trkDrCut, float trkPtCut,
        const edm::ValueMap<float>& timeValueMap,
        const edm::ValueMap<float>& timeErrorMap,
        const edm::ValueMap<float>& timeQualityMap,
        float trackPtCut, float trackEtaCut, float timeQualityCut) {

    float nmatch = 0;
    for (const reco::Candidate* dau : *(gv.daughters())) {
        for (reco::Vertex::trackRef_iterator trk_it = v.tracks_begin(); trk_it != v.tracks_end(); trk_it++) {
            reco::TrackBaseRef trkRef = *trk_it;
            if (!goodTrack(trkRef, timeValueMap, timeErrorMap, timeQualityMap, trackPtCut, trackEtaCut, timeQualityCut)) continue;
            if (matchRecoToGenTrack(trkRef, dau, trkDrCut, trkPtCut)) {
                nmatch++;
                break;
            }
        }
    }
    return (nmatch/(float)gv.nDaughters()) >= matchFrac;
}


// Helper functions


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


template <class T>
float ntuple_SV::trackD3d(const T& trkRef) {

    return TMath::Sqrt(trkRef->dxy()*trkRef->dxy() + trkRef->dz()*trkRef->dz());
}


template <class T>
float ntuple_SV::trackD3dErr(const T& trkRef) {

    float dxy2 = trkRef->dxy()*trkRef->dxy();
    float dz2 = trkRef->dz()*trkRef->dz();
    float dxy2err = 2*trkRef->dxyError()*trkRef->dxy();
    float dz2err = 2*trkRef->dzError()*trkRef->dz();
    float d3d2 = dxy2 + dz2;
    float d3d2err = TMath::Sqrt(dxy2err*dxy2err + dz2err*dz2err);
    return 0.5*d3d2err/TMath::Sqrt(d3d2);
}


float ntuple_SV::vertexDxy(const GenVertex& gv, const reco::Vertex& v) {

    float dx = v.x() - gv.x();
    float dy = v.y() - gv.y();
    return TMath::Sqrt(dx*dx + dy*dy);
}


float ntuple_SV::vertexDxyErr(const GenVertex& gv, const reco::Vertex& v) {

    float dx = v.x() - gv.x();
    float dy = v.y() - gv.y();
    float dxerr = v.xError();
    float dyerr = v.yError();
    float dx2 = dx*dx;
    float dy2 = dy*dy;
    float dx2err = 2*dx*dxerr;
    float dy2err = 2*dy*dyerr;
    float dxy2 = dx2 + dy2;
    float dxy2err = TMath::Sqrt(dx2err*dx2err + dy2err*dy2err);
    return 0.5*dxy2err/TMath::Sqrt(dxy2);
}


float ntuple_SV::vertexD3d(const GenVertex& gv, const reco::Vertex& v) {

    float dx = v.x() - gv.x();
    float dy = v.y() - gv.y();
    float dz = v.z() - gv.z();
    return TMath::Sqrt(dx*dx + dy*dy * dz*dz);
}


float ntuple_SV::vertexD3dErr(const GenVertex& gv, const reco::Vertex& v) {

    float dxy = ntuple_SV::vertexDxy(gv, v);
    float dz = v.z() - gv.z();
    float dxyerr = ntuple_SV::vertexDxyErr(gv, v);
    float dzerr = v.zError();
    float dxy2 = dxy*dxy;
    float dz2 = dz*dz;
    float dxy2err = 2*dxy*dxyerr;
    float dz2err = 2*dz*dzerr;
    float d3d2 = dxy2 + dz2;
    float d3d2err = TMath::Sqrt(dxy2err*dxy2err + dz2err*dz2err);
    return 0.5*d3d2err/TMath::Sqrt(d3d2);
}


Measurement1D ntuple_SV::vertexDxy(const reco::Vertex& sv, const reco::Vertex& pv) {

    VertexDistanceXY dist;
    reco::Vertex::CovarianceMatrix csv;
    sv.fill(csv);
    reco::Vertex svtx(sv.position(), csv);
    return dist.distance(svtx, pv);
}


Measurement1D ntuple_SV::vertexD3d(const reco::Vertex& sv, const reco::Vertex& pv) {

    VertexDistance3D dist;
    reco::Vertex::CovarianceMatrix csv;
    sv.fill(csv);
    reco::Vertex svtx(sv.position(), csv);
    return dist.distance(svtx, pv);
}


Measurement1D ntuple_SV::vertexDxy(const reco::VertexCompositePtrCandidate& cand, const reco::Vertex& pv) {

    VertexDistanceXY dist;
    reco::Vertex::CovarianceMatrix csv;
    cand.fillVertexCovariance(csv);
    reco::Vertex svtx(cand.vertex(), csv);
    return dist.distance(svtx, pv);
}


Measurement1D ntuple_SV::vertexD3d(const reco::VertexCompositePtrCandidate& cand, const reco::Vertex& pv) {

    VertexDistance3D dist;
    reco::Vertex::CovarianceMatrix csv;
    cand.fillVertexCovariance(csv);
    reco::Vertex svtx(cand.vertex(), csv);
    return dist.distance(svtx, pv);
}


bool ntuple_SV::gvCompareDxyDxySig(const GenVertex& gva, const GenVertex& gvb) {

    reco::Vertex pv = *spvp_;
    float adxy = ntuple_SV::vertexDxy(gva, pv);
    float bdxy = ntuple_SV::vertexDxy(gvb, pv);
    float aerr = ntuple_SV::vertexDxyErr(gva, pv);
    float berr = ntuple_SV::vertexDxyErr(gvb, pv);
    float asig = ntuple_SV::catchInfs(adxy / aerr, 0.0);
    float bsig = ntuple_SV::catchInfs(bdxy / berr, 0.0);
    return bsig < asig;
}


bool ntuple_SV::svCompareDxyDxySig(const reco::Vertex& sva, const reco::Vertex& svb) {

    reco::Vertex pv = *spvp_;
    float adxy = ntuple_SV::vertexDxy(sva, pv).value();
    float bdxy = ntuple_SV::vertexDxy(svb, pv).value();
    float aerr = ntuple_SV::vertexDxy(sva, pv).error();
    float berr = ntuple_SV::vertexDxy(svb, pv).error();
    float asig = ntuple_SV::catchInfs(adxy / aerr, 0.0);
    float bsig = ntuple_SV::catchInfs(bdxy / berr, 0.0);
    return bsig < asig;
}


bool ntuple_SV::candCompareDxyDxySig(const reco::VertexCompositePtrCandidate& sva, const reco::VertexCompositePtrCandidate& svb) {

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


float ntuple_SV::jetRadius(float jet_rad, const pat::Jet& jet) {

    float jet_radius = jet_rad;
    if (jet_radius < 0.0) {
        for (unsigned idau = 0; idau < jet.numberOfDaughters(); ++idau) {
            double dR = reco::deltaR(*jet.daughter(idau), jet);
            if (dR > jet_radius)
                jet_radius = dR;
        }
    }
    return jet_radius;
}
