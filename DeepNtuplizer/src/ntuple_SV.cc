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

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "TVector3.h"


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


const reco::Vertex* ntuple_SV::spvp_;


ntuple_SV::ntuple_SV(std::string prefix, double jetR) : ntuple_content(jetR), sv_num_(0) {

    prefix_ = prefix;
}


ntuple_SV::~ntuple_SV() {}


void ntuple_SV::getInput(const edm::ParameterSet& iConfig) {}


void ntuple_SV::initBranches(TTree* tree) {

    // SV candidates
    addBranch(tree, (prefix_ + "n_sv").c_str(), &sv_num_, (prefix_ + "sv_num_/I").c_str());
    addBranch(tree, (prefix_ + "nsv").c_str(), &nsv_, (prefix_ + "nsv_/F").c_str());
    addBranch(tree, (prefix_ + "sv_pt").c_str(), &sv_pt_, (prefix_ + "sv_pt_[" + prefix_ + "sv_num_]/F").c_str());
    addBranch(tree, (prefix_ + "sv_eta").c_str(), &sv_eta_, (prefix_ + "sv_eta_[" + prefix_ + "sv_num_]/F").c_str());
    addBranch(tree, (prefix_ + "sv_phi").c_str(), &sv_phi_, (prefix_ + "sv_phi_[" + prefix_ + "sv_num_]/F").c_str());
    addBranch(tree, (prefix_ + "sv_e").c_str(), &sv_e_, (prefix_ + "sv_e_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_etarel").c_str(), &sv_etarel_, (prefix_ + "sv_etarel_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_phirel").c_str(), &sv_phirel_, (prefix_ + "sv_phirel_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_deltaR").c_str(), &sv_deltaR_, (prefix_ + "sv_deltaR_[" + prefix_ + "sv_num_]/F").c_str());
    addBranch(tree, (prefix_ + "sv_mass").c_str(), &sv_mass_, (prefix_ + "sv_mass_[" + prefix_ + "sv_num_]/F").c_str());
    addBranch(tree, (prefix_ + "sv_ntracks").c_str(), &sv_ntracks_, (prefix_ + "sv_ntracks_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_nMatchPFCand").c_str(), &sv_nMatchPFCand_, (prefix_ + "sv_nMatchPFCand_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_nMatchLostTrk").c_str(), &sv_nMatchLostTrk_, (prefix_ + "sv_nMatchLostTrk_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_nUnmatchedTrk").c_str(), &sv_nUnmatchedTrk_, (prefix_ + "sv_nUnmatchedTrk_[" + prefix_ + "sv_num_]/F").c_str());
    addBranch(tree, (prefix_ + "sv_chi2").c_str(), &sv_chi2_, (prefix_ + "sv_chi2_[" + prefix_ + "sv_num_]/F").c_str());
    addBranch(tree, (prefix_ + "sv_ndf").c_str(), &sv_ndf_, (prefix_ + "sv_ndf_[" + prefix_ + "sv_num_]/F").c_str());
    addBranch(tree, (prefix_ + "sv_normchi2").c_str(), &sv_normchi2_, (prefix_ + "sv_normchi2_[" + prefix_ + "sv_num_]/F").c_str());
    addBranch(tree, (prefix_ + "sv_dxy").c_str(), &sv_dxy_, (prefix_ + "sv_dxy_[" + prefix_ + "sv_num_]/F").c_str());
    addBranch(tree, (prefix_ + "sv_dxyerr").c_str(), &sv_dxyerr_, (prefix_ + "sv_dxyerr_[" + prefix_ + "sv_num_]/F").c_str());
    addBranch(tree, (prefix_ + "sv_dxysig").c_str(), &sv_dxysig_, (prefix_ + "sv_dxysig_[" + prefix_ + "sv_num_]/F").c_str());
    addBranch(tree, (prefix_ + "sv_d3d").c_str(), &sv_d3d_, (prefix_ + "sv_d3d_[" + prefix_ + "sv_num_]/F").c_str());
    addBranch(tree, (prefix_ + "sv_d3derr").c_str(), &sv_d3derr_, (prefix_ + "sv_d3err_[" + prefix_ + "sv_num_]/F").c_str());
    addBranch(tree, (prefix_ + "sv_d3dsig").c_str(), &sv_d3dsig_, (prefix_ + "sv_d3dsig_[" + prefix_ + "sv_num_]/F").c_str());
    addBranch(tree, (prefix_ + "sv_costhetasvpv").c_str(), &sv_costhetasvpv_, (prefix_ + "sv_costhetasvpv_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_enratio").c_str(), &sv_enratio_, (prefix_ + "sv_enratio_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_hcal_frac").c_str(), &sv_hcal_frac_, (prefix_ + "sv_hcal_frac_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_calo_frac").c_str(), &sv_calo_frac_, (prefix_ + "sv_calo_frac_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_dz").c_str(), &sv_dz_, (prefix_ + "sv_dz_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_pfd2dval").c_str(), &sv_pfd2dval_, (prefix_ + "sv_pfd2dval_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_pfd2dsig").c_str(), &sv_pfd2dsig_, (prefix_ + "sv_pfd2dsig_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_pfd3dval").c_str(), &sv_pfd3dval_, (prefix_ + "sv_pfd3dval_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_pfd3dsig").c_str(), &sv_pfd3dsig_, (prefix_ + "sv_pfd3dsig_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_puppiw").c_str(), &sv_puppiw_, (prefix_ + "sv_puppiw_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_charge_sum").c_str(), &sv_charge_sum_, (prefix_ + "sv_charge_sum_[" + prefix_ + "sv_num_]/F").c_str());
    // addBranch(tree, (prefix_ + "sv_time").c_str(), &sv_time_, (prefix_ + "sv_time_[" + prefix_ + "sv_num_]/F").c_str());

    // No leaf type exists for STL containers
    tree->Branch((prefix_ + "sv_jetIdx").c_str(), &sv_jetIdx_);
    tree->Branch((prefix_ + "sv_matchJetIdx").c_str(), &sv_matchJetIdx_);

    // add2DBranch(tree, (prefix_ + "sv_ptGenVsPtReco").c_str(), &sv_ptGenVsPtReco_, (prefix_ + "sv_ptGenVsPtReco_[" + prefix_ + "sv_num_]/F").c_str());
}


void ntuple_SV::readSetup(const edm::EventSetup& iSetup) {

    builder_ = iSetup.getHandle(track_builder_token_);
}


void ntuple_SV::readEvent(const edm::Event& iEvent) {

    iEvent.getByToken(pf_cand_token_, pf_cand_);
    iEvent.getByToken(lost_tracks_token_, lost_tracks_);
    iEvent.getByToken(pf_mcmatch_token_, pf_mcmatch_);
    iEvent.getByToken(lt_mcmatch_token_, lt_mcmatch_);
}


void ntuple_SV::initContainers() {

    sv_jetIdx_ = new std::vector<int>;
    sv_matchJetIdx_ = new std::vector<int>;
}


void ntuple_SV::clearContainers() {

    sv_jetIdx_->clear();
    sv_matchJetIdx_->clear();
}


void ntuple_SV::deleteContainers() {

    delete sv_jetIdx_;
    delete sv_matchJetIdx_;
}


// Use either of these functions
// bool ntuple_SV::fillBranches(const pat::Jet& jet, const size_t& jetidx, const edm::View<pat::Jet>* coll) {
bool ntuple_SV::fillBranches(const pat::Jet& jet, const size_t& jetidx, const edm::View<pat::Jet>* coll, float EventTime) {

    const float jet_uncorr_e = jet.correctedJet("Uncorrected").energy();

    // Sorting described here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideOfflinePrimaryVertexProduction
    const reco::Vertex& pv = vertices()->at(0); // Most likely the signal vertex

    math::XYZVector jetDir = jet.momentum().Unit();
    GlobalVector jetRefTrackDir(jet.px(), jet.py(), jet.pz());

    sv_num_ = 0;

    reco::VertexCompositePtrCandidateCollection cpvtx = *secVertices();

    spvp_ = &vertices()->at(0);
    std::sort(cpvtx.begin(), cpvtx.end(), ntuple_SV::compareDxyDxyErr);

    SVTrackInfoBuilder trackinfo(builder_);

    const pat::PackedCandidateCollection& pfCands(*(pf_cand_.product()));
    const pat::PackedCandidateCollection& lostTracks(*(lost_tracks_.product()));

    const edm::Association<reco::GenParticleCollection>& PFCandMCTruth(*(pf_mcmatch_.product()));
    const edm::Association<reco::GenParticleCollection>& LostTrackMCTruth(*(lt_mcmatch_.product()));

    float etasign = 1;
    etasign++; // avoid unused warning
    if (jet.eta() < 0) etasign =- 1;

    double jet_radius = jetR();
    if (jet_radius < 0) {
        // subjets: use maxDR(subjet, pfcand)
        for (unsigned idau = 0; idau < jet.numberOfDaughters(); ++idau) {
            double dR = reco::deltaR(*jet.daughter(idau), jet);
            if (dR > jet_radius)
                jet_radius = dR;
        }
    }

    int nSV = -2;
    const reco::CandSecondaryVertexTagInfo* candSVTagInfo = jet.tagInfoCandSecondaryVertex("pfInclusiveSecondaryVertexFinder");
    if (candSVTagInfo != nullptr) nSV = candSVTagInfo->nVertices();
    if (nSV > 0 && candSVTagInfo->vertexTracks().size() == 0) nSV = -1;

    // std::cout << " NTuple jet   pt eta phi " << jet.pt() << " " << jet.eta() << " " << jet.phi() << "   nSV " << nSV << std::endl;

    for (const reco::VertexCompositePtrCandidate& sv : cpvtx) {
        if (reco::deltaR(sv, jet) > jet_radius) continue;

        if ((int) max_sv > sv_num_) { // limit number of SVs
            std::cout << "new vertex ---------------------------------------------------------" << std::endl;
            // sv_pt_[sv_num_] = sv.pt();
            // sv_eta_[sv_num_] = sv.eta();
            // sv_phi_[sv_num_] = sv.phi();
            sv_etarel_[sv_num_] = catchInfsAndBound(fabs(sv.eta() - jet.eta()) - 0.5, 0, -2, 0);
            sv_phirel_[sv_num_] = catchInfsAndBound(fabs(reco::deltaPhi(sv.phi(), jet.phi())) - 0.5, 0, -2, 0);
            sv_deltaR_[sv_num_] = catchInfsAndBound(fabs(reco::deltaR(sv, jet)) - 0.5, 0, -2, 0);
            // sv_mass_[sv_num_] = sv.mass();
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
            // sv_costhetasvpv_[sv_num_] = vertexDdotP(sv, pv); // the pointing angle (i.e. the angle between the sum of the momentum of the tracks in the SV and the flight direction betwen PV and SV)
            sv_enratio_[sv_num_] = sv.energy() / jet_uncorr_e;
            // sv_e_[sv_num_] = sv.energy();

            float calo_frac = 0.0;
            float hcal_frac = 0.0;
            float puppiw = 0.0;
            float charge = 0.0;
            float dz = 0.0;

            float pfd3dval = 0.0;
            float pfd3dsig = 0.0;
            float pfd2dval = 0.0;
            float pfd2dsig = 0.0;
            float pfcount = 0.0;

            float pfCandMatchCount = 0.0;
            float lostTrackMatchCount = 0.0;
            for (unsigned idx = 0; idx < sv.numberOfDaughters(); ++idx) {
                const pat::PackedCandidate* PackedCandidate_ = dynamic_cast<const pat::PackedCandidate*>(sv.daughter(idx));

                // int pfCandMatchIdx = findPFCandIdx(*PackedCandidate_, pfCands);
                // int ltMatchIdx = findLostTrackIdx(*PackedCandidate_, lostTracks);
                // if (pfCandMatchIdx >= 0) {
                //     pat::PackedCandidateRef tempTrkRef = pat::PackedCandidateRef(pf_cand_, (unsigned int) pfCandMatchIdx);
                //     reco::GenParticleRef trkTruthRef = PFCandMCTruth[tempTrkRef];
                //     if (trkTruthRef.isNonnull()) {
                //         // std::cout << "matched pf candidate" << std::endl;
                //         // std::cout << "pt: " << tempTrkRef->pt() << ", " << trkTruthRef->pt() << std::endl;
                //         std::cout << "number of mothers = " << trkTruthRef->numberOfMothers() << std::endl;
                //         if (trkTruthRef->numberOfMothers() == 1) {
                //             const reco::Candidate* truthMother = trkTruthRef->mother();
                //             std::cout << "mother pt = " << truthMother->pt() << std::endl;
                //         }
                //         else {
                //             std::cout << "wrong number of mothers!!" << std::endl;
                //         }
                //         pfCandMatchCount++;
                //     }
                //     else {
                //         // std::cout << "no matched pf candidates found" << std::endl;
                //     }
                // }
                // else if (ltMatchIdx >= 0) {
                //     pat::PackedCandidateRef tempTrkRef = pat::PackedCandidateRef(lost_tracks_, (unsigned int) ltMatchIdx);
                //     reco::GenParticleRef trkTruthRef = LostTrackMCTruth[tempTrkRef];
                //     if (trkTruthRef.isNonnull()) {
                //         std::cout << "matched lost track" << std::endl;
                //         std::cout << "pt: " << tempTrkRef->pt() << ", " << trkTruthRef->pt() << std::endl;
                //         lostTrackMatchCount++;
                //     }
                //     else {
                //         // std::cout << "no matched lost tracks found" << std::endl;
                //     }
                // }
                // else {
                //     // std::cout << "no matched pf candidates or lost tracks found" << std::endl;
                // }

                calo_frac = calo_frac + PackedCandidate_->caloFraction();
                hcal_frac = hcal_frac + PackedCandidate_->hcalFraction();
                puppiw = puppiw + PackedCandidate_->puppiWeight();
                charge = charge + PackedCandidate_->charge();
                dz = dz + PackedCandidate_->dz();
                if (PackedCandidate_->charge() != 0 and PackedCandidate_->pt() > 0.95) { // TODO: understand these "track" cuts
                    trackinfo.buildTrackInfo(PackedCandidate_, jetDir, jetRefTrackDir, pv);
                    pfd3dval = pfd3dval + catchInfsAndBound(trackinfo.getTrackSip3dVal(), 0, -1, 1e5);
                    pfd3dsig = pfd3dsig + catchInfsAndBound(trackinfo.getTrackSip3dSig(), 0, -1, 4e4);
                    pfd2dval = pfd2dval + catchInfsAndBound(trackinfo.getTrackSip2dVal(), 0, -1, 70);
                    pfd2dsig = pfd2dsig + catchInfsAndBound(trackinfo.getTrackSip2dSig(), 0, -1, 4e4);
                    pfcount = pfcount + 1.0;
                }
            }

            sv_nMatchPFCand_[sv_num_] = pfCandMatchCount;
            sv_nMatchLostTrk_[sv_num_] = lostTrackMatchCount;
            sv_nUnmatchedTrk_[sv_num_] = sv.numberOfDaughters() - (pfCandMatchCount + lostTrackMatchCount);

            sv_calo_frac_[sv_num_] = calo_frac / sv.numberOfDaughters();
            sv_hcal_frac_[sv_num_] = hcal_frac / sv.numberOfDaughters();
            sv_puppiw_[sv_num_] = puppiw / sv.numberOfDaughters();
            sv_dz_[sv_num_] = dz / sv.numberOfDaughters();
            sv_charge_sum_[sv_num_] = charge;

            sv_pfd3dval_[sv_num_] = pfd3dval / pfcount;
            sv_pfd3dsig_[sv_num_] = pfd3dsig / pfcount;
            sv_pfd2dval_[sv_num_] = pfd2dval / pfcount;
            sv_pfd2dsig_[sv_num_] = pfd2dsig / pfcount;

            // Get the vertex time
            // Matching VertexCompositePtrCandidate (reconstructed) and tagInfoCandSecondaryVertex (information used to compute b (&c?) tag discriminator)...
            // Maybe a trick to get the SVTagInfo features for each reconstructed vertex, as we loop over the latter and not all features are calculated -- by matching in eta, phi, you can get the same vertex and combine information from different collections
            float vertex_time = 0;
            float vertex_timeWeight = 0;
            float vertex_timeNtk = 0;

            if (nSV > 0 && sv.pt() > 0.0) {
                for (unsigned int isv = 0; isv < candSVTagInfo->nVertices(); ++isv) {
                    float dSVpt = TMath::Abs(candSVTagInfo->secondaryVertex(isv).pt() / sv.pt() - 1.0);
                    float dSVeta = TMath::Abs(candSVTagInfo->secondaryVertex(isv).eta() - sv.eta());
                    float dSVphi = TMath::Abs(candSVTagInfo->secondaryVertex(isv).phi() - sv.phi());
                    if (dSVphi > 3.141593) dSVphi -= 2.0 * 3.141593;
                    if (!(dSVpt < 0.01 && dSVeta < 0.01 && dSVphi < 0.01)) continue;
                    // std::cout << "  => matched sv " << sv_num_ << " to " << isv << std::endl;

                    for (unsigned int it = 0; it < candSVTagInfo->nVertexTracks(isv); ++it) {
                        for (unsigned int i = 0; i < jet.numberOfDaughters(); i++) {
                            const pat::PackedCandidate* PackedCandidate = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(i));
                            if (!PackedCandidate) continue;
                            if (PackedCandidate->charge() == 0) continue;
                            auto track = PackedCandidate->bestTrack();
                            if (!track) continue;
                            if (candSVTagInfo->vertexTracks(isv)[it]->charge() != track->charge()) continue;
                            float track_time = track->t0();
                            float track_timeError = track->covt0t0();
                            float track_pt = track->pt();
                            float time_weight = track_pt * track_pt;
                            if (!(track_timeError > 0.0 && abs(track_time) < 1)) continue;

                            float dpt = TMath::Abs(candSVTagInfo->vertexTracks(isv)[it]->pt() / track->pt() - 1.0);
                            float deta = TMath::Abs(candSVTagInfo->vertexTracks(isv)[it]->eta() - track->eta());
                            float dphi = TMath::Abs(candSVTagInfo->vertexTracks(isv)[it]->phi() - track->phi());
                            if (dphi > 3.141593) dphi -= 2.0 * 3.141593;
                            if (dpt < 0.01 && deta < 0.01 && dphi < 0.01) {
                                vertex_timeNtk += 1;
                                vertex_timeWeight += time_weight;
                                vertex_time += track_time * time_weight;
                                // std::cout << "  => matched track " << it << " to " << i << " time " << track_time << std::endl;
                            }
                        } // end loop on all tracks in jet
                    } // end loop on tracks from SV in jet
                } // end loop on SVs in jet
                if (vertex_timeNtk > 0 && EventTime > -1) {
                    vertex_time = vertex_time / vertex_timeWeight - EventTime;
                    vertex_time = TMath::Abs(vertex_time);
                } // time of flight?
                else vertex_time = -1; 
            }
            else vertex_time = -1.0;

            // std::cout << " NTuple sv " << sv_num_ << " pt eta phi " << sv.pt() << " " << sv.eta() << " " << sv.phi() << " time " << vertex_time << std::endl;

            sv_time_[sv_num_] = vertex_time;

            sv_num_++;
        } // end if max_sv > sv_num_
    } // end of looping over the reconstructed secondary vertices
    nsv_ = sv_num_;

    return true;
}


bool ntuple_SV::fillBranches() {

    std::cout << "the other fillBranches() in ntuple_SV" << std::endl;

    clearContainers();

    // Sorting described here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideOfflinePrimaryVertexProduction
    const reco::Vertex& pv = vertices()->at(0); // Most likely the signal vertex
    spvp_ = &vertices()->at(0);

    reco::VertexCompositePtrCandidateCollection cpvtx = *secVertices();
    std::sort(cpvtx.begin(), cpvtx.end(), ntuple_SV::compareDxyDxyErr);

    const edm::View<pat::Jet> jetCollection = *jets();

    sv_num_ = 0;
    for (const reco::VertexCompositePtrCandidate& sv : cpvtx) {
        if (sv_num_ < (int) max_sv) {
            sv_pt_[sv_num_] = sv.pt();
            sv_eta_[sv_num_] = sv.eta();
            sv_phi_[sv_num_] = sv.phi();
            sv_mass_[sv_num_] = sv.mass();
            sv_e_[sv_num_] = sv.energy();
            sv_ntracks_[sv_num_] = sv.numberOfDaughters();
            sv_chi2_[sv_num_] = sv.vertexChi2();
            sv_ndf_[sv_num_] = sv.vertexNdof();
            sv_normchi2_[sv_num_] = catchInfsAndBound(sv_chi2_[sv_num_] / sv_ndf_[sv_num_], 1000, -1000, 1000);
            sv_dxy_[sv_num_] = vertexDxy(sv, pv).value();
            sv_dxyerr_[sv_num_] = catchInfsAndBound(vertexDxy(sv, pv).error() - 2, 0, -2, 0);
            sv_dxysig_[sv_num_] = catchInfsAndBound(sv_dxy_[sv_num_] / vertexDxy(sv, pv).error(), 0, -1, 800);
            sv_d3d_[sv_num_] = vertexD3d(sv, pv).value();
            sv_d3derr_[sv_num_] = catchInfsAndBound(vertexD3d(sv, pv).error() - 2, 0, -2, 0);
            sv_d3dsig_[sv_num_] = catchInfsAndBound(vertexD3d(sv, pv).value() / vertexD3d(sv, pv).error(), 0, -1, 800);
            sv_costhetasvpv_[sv_num_] = vertexDdotP(sv, pv); // the pointing angle (i.e. the angle between the sum of the momentum of the tracks in the SV and the flight direction betwen PV and SV)

            // Match to a jet based on deltaR
            for (unsigned int j = 0; j < jetCollection.size(); j++) {

                const pat::Jet& jet = jetCollection.at(j);

                double jet_radius = jetR();
                if (jet_radius < 0) {
                    // subjets: use maxDR(subjet, pfcand)
                    for (unsigned idau = 0; idau < jet.numberOfDaughters(); ++idau) {
                        double dR = reco::deltaR(*jet.daughter(idau), jet);
                        if (dR > jet_radius)
                            jet_radius = dR;
                    }
                }

                if (reco::deltaR(sv, jet) > jet_radius) continue;

                sv_jetIdx_->push_back((int) j);
                sv_matchJetIdx_->push_back(sv_num_);
            }
        }
        sv_num_++;
    }
    nsv_ = sv_num_;

    return false;
}


// Helpers functions


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


int ntuple_SV::findPFCandIdx(const pat::PackedCandidate& trk, const pat::PackedCandidateCollection& pcands) {

    int nmatches = 0;
    int matchIdx = 0;
    for (unsigned int trkIdx = 0; trkIdx < pcands.size(); trkIdx++) {
        if (&trk == &pcands.at(trkIdx)) {
            matchIdx = trkIdx;
            nmatches += 1;
        }
    }
    if (nmatches != 1) {
        // std::cout << "ntuple_SV.cc: 0 or more than one PF Candidate match found! Returning -1." << std::endl;
        return -1;
    }
    return matchIdx;
}


int ntuple_SV::findLostTrackIdx(const pat::PackedCandidate& trk, const pat::PackedCandidateCollection& lts) {

    int nmatches = 0;
    int matchIdx = 0;
    for (unsigned int trkIdx = 0; trkIdx < lts.size(); trkIdx++) {
        if (&trk == &lts.at(trkIdx)) {
            matchIdx = trkIdx;
            nmatches += 1;
        }
    }
    if (nmatches != 1) {
        // std::cout << "ntuple_SV.cc: 0 or more than one Lost Track match found! Returning -1." << std::endl;
        return -1;
    }
    return matchIdx;
}
