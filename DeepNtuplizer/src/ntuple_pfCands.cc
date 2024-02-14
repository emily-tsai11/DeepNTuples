#include "../interface/ntuple_pfCands.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "../interface/sorting_modules.h"


#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "TVector3.h"

#include "DataFormats/GeometrySurface/interface/Line.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalTrajectoryExtrapolatorToLine.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"

class TrackInfoBuilder{
public:
    TrackInfoBuilder(edm::ESHandle<TransientTrackBuilder> & build):
        builder(build),
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

        trackJetDistVal_(0),
        trackJetDistSig_(0)
{


}

    void buildTrackInfo(const pat::PackedCandidate* PackedCandidate_ ,const math::XYZVector&  jetDir, GlobalVector refjetdirection, const reco::Vertex & pv){
			TVector3 jetDir3(jetDir.x(),jetDir.y(),jetDir.z());
			if(!PackedCandidate_->hasTrackDetails()) {
				TVector3 trackMom3(
					PackedCandidate_->momentum().x(),
					PackedCandidate_->momentum().y(),
					PackedCandidate_->momentum().z()
					);
				trackMomentum_=PackedCandidate_->p();
				trackEta_= PackedCandidate_->eta();
				trackEtaRel_=reco::btau::etaRel(jetDir, PackedCandidate_->momentum());
				trackPtRel_=trackMom3.Perp(jetDir3);
				trackPPar_=jetDir.Dot(PackedCandidate_->momentum());
				trackDeltaR_=reco::deltaR(PackedCandidate_->momentum(), jetDir);
				trackPtRatio_=trackMom3.Perp(jetDir3) / PackedCandidate_->p();
				trackPParRatio_=jetDir.Dot(PackedCandidate_->momentum()) / PackedCandidate_->p();
				trackSip2dVal_=0.;
				trackSip2dSig_=0.;
				trackSip3dVal_=0.;
				trackSip3dSig_=0.;
				trackJetDistVal_=0.;
				trackJetDistSig_=0.;
				return;
			}

        const reco::Track & PseudoTrack =  PackedCandidate_->pseudoTrack();

        reco::TransientTrack transientTrack;
        transientTrack=builder->build(PseudoTrack);
        Measurement1D meas_ip2d=IPTools::signedTransverseImpactParameter(transientTrack, refjetdirection, pv).second;
        Measurement1D meas_ip3d=IPTools::signedImpactParameter3D(transientTrack, refjetdirection, pv).second;
        Measurement1D jetdist=IPTools::jetTrackDistance(transientTrack, refjetdirection, pv).second;
        math::XYZVector trackMom = PseudoTrack.momentum();
        double trackMag = std::sqrt(trackMom.Mag2());
        TVector3 trackMom3(trackMom.x(),trackMom.y(),trackMom.z());


        trackMomentum_=std::sqrt(trackMom.Mag2());
        trackEta_= trackMom.Eta();
        trackEtaRel_=reco::btau::etaRel(jetDir, trackMom);
        trackPtRel_=trackMom3.Perp(jetDir3);
        trackPPar_=jetDir.Dot(trackMom);
        trackDeltaR_=reco::deltaR(trackMom, jetDir);
        trackPtRatio_=trackMom3.Perp(jetDir3) / trackMag;
        trackPParRatio_=jetDir.Dot(trackMom) / trackMag;
        trackSip2dVal_=(meas_ip2d.value());

        trackSip2dSig_=(meas_ip2d.significance());
        trackSip3dVal_=(meas_ip3d.value());


        trackSip3dSig_=meas_ip3d.significance();
        trackJetDistVal_= jetdist.value();
        trackJetDistSig_= jetdist.significance();

    }

    const float& getTrackDeltaR() const {return trackDeltaR_;}
    const float& getTrackEta() const {return trackEta_;}
    const float& getTrackEtaRel() const {return trackEtaRel_;}
    const float& getTrackJetDistSig() const {return trackJetDistSig_;}
    const float& getTrackJetDistVal() const {return trackJetDistVal_;}
    const float& getTrackMomentum() const {return trackMomentum_;}
    const float& getTrackPPar() const {return trackPPar_;}
    const float& getTrackPParRatio() const {return trackPParRatio_;}
    const float& getTrackPtRatio() const {return trackPtRatio_;}
    const float& getTrackPtRel() const {return trackPtRel_;}
    const float& getTrackSip2dSig() const {return trackSip2dSig_;}
    const float& getTrackSip2dVal() const {return trackSip2dVal_;}
    const float& getTrackSip3dSig() const {return trackSip3dSig_;}
    const float& getTrackSip3dVal() const {return trackSip3dVal_;}

private:

    edm::ESHandle<TransientTrackBuilder>& builder;
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

    float trackJetDistVal_;
    float trackJetDistSig_;

};




void ntuple_pfCands::readSetup(const edm::EventSetup& iSetup){

	builder = iSetup.getHandle(track_builder_token_);
//    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);

}

void ntuple_pfCands::getInput(const edm::ParameterSet& iConfig){
	min_candidate_pt_ = (iConfig.getParameter<double>("minCandidatePt"));
}

void ntuple_pfCands::initBranches(TTree* tree){

    addBranch(tree,"n_Cpfcand", &n_Cpfcand_,"n_Cpfcand_/I");

    addBranch(tree,"nCpfcand", &nCpfcand_,"nCpfcand_/F");

    addBranch(tree,"Cpfcan_pt", &Cpfcan_pt_,"Cpfcan_pt_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_eta", &Cpfcan_eta_,"Cpfcan_eta_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_phi", &Cpfcan_phi_,"Cpfcan_phi_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_ptrel", &Cpfcan_ptrel_,"Cpfcan_ptrel_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_e", &Cpfcan_e_,"Cpfcan_e_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_erel", &Cpfcan_erel_,"Cpfcan_erel_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_phirel",&Cpfcan_phirel_,"Cpfcan_phirel_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_etarel",&Cpfcan_etarel_,"Cpfcan_etarel_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_deltaR",&Cpfcan_deltaR_,"Cpfcan_deltaR_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_puppiw",&Cpfcan_puppiw_,"Cpfcan_puppiw_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_dxy",&Cpfcan_dxy_,"Cpfcan_dxy_[n_Cpfcand_]/F");

    addBranch(tree,"Cpfcan_dxyerrinv",&Cpfcan_dxyerrinv_,"Cpfcan_dxyerrinv_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_dxysig",&Cpfcan_dxysig_,"Cpfcan_dxysig_[n_Cpfcand_]/F");

    addBranch(tree,"Cpfcan_dz",&Cpfcan_dz_,"Cpfcan_dz_[n_Cpfcand_]/F");

    addBranch(tree,"Cpfcan_VTX_ass",&Cpfcan_VTX_ass_,"Cpfcan_VTX_ass_[n_Cpfcand_]/F");

    addBranch(tree,"Cpfcan_fromPV",&Cpfcan_fromPV_,"Cpfcan_fromPV_[n_Cpfcand_]/F");

    addBranch(tree,"Cpfcan_drminsv",&Cpfcan_drminsv_,"Cpfcan_drminsv_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_distminsv",&Cpfcan_distminsv_,"Cpfcan_distminsv_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_distminsv2",&Cpfcan_distminsv2_,"Cpfcan_distminsv2_[n_Cpfcand_]/F");    

    //commented ones don't work
    /**///addBranch(tree,"Cpfcan_vertexChi2",&Cpfcan_vertexChi2_,"Cpfcan_vertexChi2_[n_Cpfcand_]/F");
    /**///addBranch(tree,"Cpfcan_vertexNdof",&Cpfcan_vertexNdof_,"Cpfcan_vertexNdof_[n_Cpfcand_]/F");
    /**///addBranch(tree,"Cpfcan_vertexNormalizedChi2",&Cpfcan_vertexNormalizedChi2_,"Cpfcan_vertexNormalizedChi2_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_vertex_rho",&Cpfcan_vertex_rho_,"Cpfcan_vertex_rho_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_vertex_phirel",&Cpfcan_vertex_phirel_,"Cpfcan_vertex_phirel_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_vertex_etarel",&Cpfcan_vertex_etarel_,"Cpfcan_vertex_etarel_[n_Cpfcand_]/F");
    /**///addBranch(tree,"Cpfcan_vertexRef_mass",&Cpfcan_vertexRef_mass_,"Cpfcan_vertexRef_mass_[n_Cpfcand_]/F");

    /*
    addBranch(tree,"Cpfcan_dptdpt",&Cpfcan_dptdpt_,"Cpfcan_dptdpt_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_detadeta",&Cpfcan_detadeta_,"Cpfcan_detadeta_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_dphidphi",&Cpfcan_dphidphi_,"Cpfcan_dphidphi_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_dxydxy",&Cpfcan_dxydxy_,"Cpfcan_dxydxy_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_dzdz",&Cpfcan_dzdz_,"Cpfcan_dzdz_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_dxydz",&Cpfcan_dxydz_,"Cpfcan_dxydz_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_dphidxy",&Cpfcan_dphidxy_,"Cpfcan_dphidxy_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_dlambdadz",&Cpfcan_dlambdadz_,"Cpfcan_dlambdadz_[n_Cpfcand_]/F");
     */

//$$
    addBranch(tree,"Cpfcan_time", &Cpfcan_time_,"Cpfcan_time_[n_Cpfcand_]/F");
//$$

    addBranch(tree,"Cpfcan_BtagPf_trackMomentum",&Cpfcan_BtagPf_trackMomentum_,"Cpfcan_BtagPf_trackMomentum_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_BtagPf_trackEta",&Cpfcan_BtagPf_trackEta_,"Cpfcan_BtagPf_trackEta_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_BtagPf_trackEtaRel",&Cpfcan_BtagPf_trackEtaRel_,"Cpfcan_BtagPf_trackEtaRel_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_BtagPf_trackPtRel",&Cpfcan_BtagPf_trackPtRel_,"Cpfcan_BtagPf_trackPtRel_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_BtagPf_trackPPar",&Cpfcan_BtagPf_trackPPar_,"Cpfcan_BtagPf_trackPPar_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_BtagPf_trackDeltaR",&Cpfcan_BtagPf_trackDeltaR_,"Cpfcan_BtagPf_trackDeltaR_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_BtagPf_trackPtRatio",&Cpfcan_BtagPf_trackPtRatio_,"Cpfcan_BtagPf_trackPtRatio_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_BtagPf_trackPParRatio",&Cpfcan_BtagPf_trackPParRatio_,"Cpfcan_BtagPf_trackPParRatio[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_BtagPf_trackSip3dVal",&Cpfcan_BtagPf_trackSip3dVal_,"Cpfcan_BtagPf_trackSip3dVal_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_BtagPf_trackSip3dSig",&Cpfcan_BtagPf_trackSip3dSig_,"Cpfcan_BtagPf_trackSip3dSig_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_BtagPf_trackSip2dVal",&Cpfcan_BtagPf_trackSip2dVal_,"Cpfcan_BtagPf_trackSip2dVal_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_BtagPf_trackSip2dSig",&Cpfcan_BtagPf_trackSip2dSig_,"Cpfcan_BtagPf_trackSip2dSig_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_BtagPf_trackDecayLen",&Cpfcan_BtagPf_trackDecayLen_,"Cpfcan_BtagPf_trackDecayLen_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_BtagPf_trackJetDistVal",&Cpfcan_BtagPf_trackJetDistVal_,"Cpfcan_BtagPf_trackJetDistVal_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_BtagPf_trackJetDistSig",&Cpfcan_BtagPf_trackJetDistSig_,"Cpfcan_BtagPf_trackJetDistSig_[n_Cpfcand_]/F");




    addBranch(tree,"Cpfcan_isMu",&Cpfcan_isMu_,"Cpfcan_isMu_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_isEl",&Cpfcan_isEl_,"Cpfcan_isEl_[n_Cpfcand_]/F");

    //in16 conversion broken
    addBranch(tree,"Cpfcan_lostInnerHits",&Cpfcan_lostInnerHits_,"Cpfcan_lostInnerHits_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_numberOfPixelHits",&Cpfcan_numberOfPixelHits_,"Cpfcan_numberOfPixelHits_[n_Cpfcand_]/F");

    addBranch(tree,"Cpfcan_chi2",&Cpfcan_chi2_,"Cpfcan_chi2_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_quality",&Cpfcan_quality_,"Cpfcan_quality_[n_Cpfcand_]/F");

    // did not give integers !!
    //  addBranch(tree,"Cpfcan_charge",&Cpfcan_charge_,"Cpfcan_charge_[n_Cpfcand_]/I");
	addBranch(tree,"Cpfcan_charge",&Cpfcan_charge_,"Cpfcan_charge_[n_Cpfcand_]/F");

    addBranch(tree,"Cpfcan_pdgID",&Cpfcan_pdgID_,"Cpfcan_pdgID_[n_Cpfcand_]/F");

    addBranch(tree,"Cpfcan_HadFrac",&Cpfcan_HadFrac_,"Cpfcan_HadFrac_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_CaloFrac",&Cpfcan_CaloFrac_,"Cpfcan_CaloFrac_[n_Cpfcand_]/F");
	

//hit pattern variables, as defined here https://github.com/cms-sw/cmssw/blob/master/DataFormats/TrackReco/interface/HitPattern.h
    //Tracker per layer
    //Pixel barrel 
    addBranch(tree,"Cpfcan_nhitpixelBarrelLayer1",&Cpfcan_nhitpixelBarrelLayer1_,"Cpfcan_nhitpixelBarrelLayer1_[n_Cpfcand_]/F"); 
    addBranch(tree,"Cpfcan_nhitpixelBarrelLayer2",&Cpfcan_nhitpixelBarrelLayer2_,"Cpfcan_nhitpixelBarrelLayer2_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_nhitpixelBarrelLayer3",&Cpfcan_nhitpixelBarrelLayer3_,"Cpfcan_nhitpixelBarrelLayer3_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_nhitpixelBarrelLayer4",&Cpfcan_nhitpixelBarrelLayer4_,"Cpfcan_nhitpixelBarrelLayer4_[n_Cpfcand_]/F");
    //Pixel Endcap 
    addBranch(tree,"Cpfcan_nhitpixelEndcapLayer1",&Cpfcan_nhitpixelEndcapLayer1_,"Cpfcan_nhitpixelEndcapLayer1_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_nhitpixelEndcapLayer2",&Cpfcan_nhitpixelEndcapLayer2_,"Cpfcan_nhitpixelEndcapLayer2_[n_Cpfcand_]/F");
    //Strip TIB
    addBranch(tree,"Cpfcan_nhitstripTIBLayer1",&Cpfcan_nhitstripTIBLayer1_,"Cpfcan_nhitstripTIBLayer1_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_nhitstripTIBLayer2",&Cpfcan_nhitstripTIBLayer2_,"Cpfcan_nhitstripTIBLayer2_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_nhitstripTIBLayer3",&Cpfcan_nhitstripTIBLayer3_,"Cpfcan_nhitstripTIBLayer3_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_nhitstripTIBLayer4",&Cpfcan_nhitstripTIBLayer4_,"Cpfcan_nhitstripTIBLayer4_[n_Cpfcand_]/F");
    //Strip TID
    addBranch(tree,"Cpfcan_nhitstripTIDLayer1",&Cpfcan_nhitstripTIDLayer1_,"Cpfcan_nhitstripTIDLayer1_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_nhitstripTIDLayer2",&Cpfcan_nhitstripTIDLayer2_,"Cpfcan_nhitstripTIDLayer2_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_nhitstripTIDLayer3",&Cpfcan_nhitstripTIDLayer3_,"Cpfcan_nhitstripTIDLayer3_[n_Cpfcand_]/F");
    //Strip TOB
    addBranch(tree,"Cpfcan_nhitstripTOBLayer1",&Cpfcan_nhitstripTOBLayer1_,"Cpfcan_nhitstripTOBLayer1_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_nhitstripTOBLayer2",&Cpfcan_nhitstripTOBLayer2_,"Cpfcan_nhitstripTOBLayer2_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_nhitstripTOBLayer3",&Cpfcan_nhitstripTOBLayer3_,"Cpfcan_nhitstripTOBLayer3_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_nhitstripTOBLayer4",&Cpfcan_nhitstripTOBLayer4_,"Cpfcan_nhitstripTOBLayer4_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_nhitstripTOBLayer5",&Cpfcan_nhitstripTOBLayer5_,"Cpfcan_nhitstripTOBLayer5_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_nhitstripTOBLayer6",&Cpfcan_nhitstripTOBLayer6_,"Cpfcan_nhitstripTOBLayer6_[n_Cpfcand_]/F");
    //Strip TEC
    addBranch(tree,"Cpfcan_nhitstripTECLayer1",&Cpfcan_nhitstripTECLayer1_,"Cpfcan_nhitstripTECLayer1_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_nhitstripTECLayer2",&Cpfcan_nhitstripTECLayer2_,"Cpfcan_nhitstripTECLayer2_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_nhitstripTECLayer3",&Cpfcan_nhitstripTECLayer3_,"Cpfcan_nhitstripTECLayer3_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_nhitstripTECLayer4",&Cpfcan_nhitstripTECLayer4_,"Cpfcan_nhitstripTECLayer4_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_nhitstripTECLayer5",&Cpfcan_nhitstripTECLayer5_,"Cpfcan_nhitstripTECLayer5_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_nhitstripTECLayer6",&Cpfcan_nhitstripTECLayer6_,"Cpfcan_nhitstripTECLayer6_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_nhitstripTECLayer7",&Cpfcan_nhitstripTECLayer7_,"Cpfcan_nhitstripTECLayer7_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_nhitstripTECLayer8",&Cpfcan_nhitstripTECLayer8_,"Cpfcan_nhitstripTECLayer8_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_nhitstripTECLayer9",&Cpfcan_nhitstripTECLayer9_,"Cpfcan_nhitstripTECLayer9_[n_Cpfcand_]/F");
    //Tracker all layers together
    //Valid hits
    addBranch(tree,"Cpfcan_numberOfValidHits",&Cpfcan_numberOfValidHits_,"Cpfcan_numberOfValidHits_[n_Cpfcand_]/F"); 
    addBranch(tree,"Cpfcan_numberOfValidTrackerHits",&Cpfcan_numberOfValidTrackerHits_,"Cpfcan_numberOfValidTrackerHits_[n_Cpfcand_]/F"); 
    addBranch(tree,"Cpfcan_numberOfValidPixelHits",&Cpfcan_numberOfValidPixelHits_,"Cpfcan_numberOfValidPixelHits_[n_Cpfcand_]/F"); 
    addBranch(tree,"Cpfcan_numberOfValidPixelBarrelHits",&Cpfcan_numberOfValidPixelBarrelHits_,"Cpfcan_numberOfValidPixelBarrelHits_[n_Cpfcand_]/F"); 
    addBranch(tree,"Cpfcan_numberOfValidPixelEndcapHits",&Cpfcan_numberOfValidPixelEndcapHits_,"Cpfcan_numberOfValidPixelEndcapHits_[n_Cpfcand_]/F"); 
    addBranch(tree,"Cpfcan_numberOfValidStripHits",&Cpfcan_numberOfValidStripHits_,"Cpfcan_numberOfValidStripHits_[n_Cpfcand_]/F"); 
    addBranch(tree,"Cpfcan_numberOfValidStripTIBHits",&Cpfcan_numberOfValidStripTIBHits_,"Cpfcan_numberOfValidStripTIBHits_[n_Cpfcand_]/F"); 
    addBranch(tree,"Cpfcan_numberOfValidStripTIDHits",&Cpfcan_numberOfValidStripTIDHits_,"Cpfcan_numberOfValidStripTIDHits_[n_Cpfcand_]/F"); 
    addBranch(tree,"Cpfcan_numberOfValidStripTOBHits",&Cpfcan_numberOfValidStripTOBHits_,"Cpfcan_numberOfValidStripTOBHits_[n_Cpfcand_]/F"); 
    addBranch(tree,"Cpfcan_numberOfValidStripTECHits",&Cpfcan_numberOfValidStripTECHits_,"Cpfcan_numberOfValidStripTECHits_[n_Cpfcand_]/F"); 
    //LayersWithMeasuremen
    addBranch(tree,"Cpfcan_trackerLayersWithMeasurementOld",&Cpfcan_trackerLayersWithMeasurementOld_,"Cpfcan_trackerLayersWithMeasurementOld_[n_Cpfcand_]/F"); 
    addBranch(tree,"Cpfcan_trackerLayersWithMeasurement",&Cpfcan_trackerLayersWithMeasurement_,"Cpfcan_trackerLayersWithMeasurement_[n_Cpfcand_]/F"); 
    addBranch(tree,"Cpfcan_pixelLayersWithMeasurementOld",&Cpfcan_pixelLayersWithMeasurementOld_,"Cpfcan_pixelLayersWithMeasurementOld_[n_Cpfcand_]/F"); 
    addBranch(tree,"Cpfcan_pixelLayersWithMeasurement",&Cpfcan_pixelLayersWithMeasurement_,"Cpfcan_pixelLayersWithMeasurement_[n_Cpfcand_]/F"); 
    addBranch(tree,"Cpfcan_stripLayersWithMeasurement",&Cpfcan_stripLayersWithMeasurement_,"Cpfcan_stripLayersWithMeasurement_[n_Cpfcand_]/F"); 
    addBranch(tree,"Cpfcan_pixelBarrelLayersWithMeasurement",&Cpfcan_pixelBarrelLayersWithMeasurement_,"Cpfcan_pixelBarrelLayersWithMeasurement_[n_Cpfcand_]/F"); 
    addBranch(tree,"Cpfcan_pixelEndcapLayersWithMeasurement",&Cpfcan_pixelEndcapLayersWithMeasurement_,"Cpfcan_pixelEndcapLayersWithMeasurement_[n_Cpfcand_]/F"); 
    addBranch(tree,"Cpfcan_stripTIBLayersWithMeasurement",&Cpfcan_stripTIBLayersWithMeasurement_,"Cpfcan_stripTIBLayersWithMeasurement_[n_Cpfcand_]/F"); 
    addBranch(tree,"Cpfcan_stripTIDLayersWithMeasurement",&Cpfcan_stripTIDLayersWithMeasurement_,"Cpfcan_stripTIDLayersWithMeasurement_[n_Cpfcand_]/F"); 
    addBranch(tree,"Cpfcan_stripTOBLayersWithMeasurement",&Cpfcan_stripTOBLayersWithMeasurement_,"Cpfcan_stripTOBLayersWithMeasurement_[n_Cpfcand_]/F"); 
    addBranch(tree,"Cpfcan_stripTECLayersWithMeasurement",&Cpfcan_stripTECLayersWithMeasurement_,"Cpfcan_stripTECLayersWithMeasurement_[n_Cpfcand_]/F"); 
    //Null
    addBranch(tree,"Cpfcan_trackerLayersNull",&Cpfcan_trackerLayersNull_,"Cpfcan_trackerLayersNull_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_pixelLayersNull",&Cpfcan_pixelLayersNull_,"Cpfcan_pixelLayersNull_[n_Cpfcand_]/F"); 
    addBranch(tree,"Cpfcan_stripLayersNull",&Cpfcan_stripLayersNull_,"Cpfcan_stripLayersNull_[n_Cpfcand_]/F"); 
    addBranch(tree,"Cpfcan_pixelBarrelLayersNull",&Cpfcan_pixelBarrelLayersNull_,"Cpfcan_pixelBarrelLayersNull_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_pixelEndcapLayersNull",&Cpfcan_pixelEndcapLayersNull_,"Cpfcan_pixelEndcapLayersNull_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_stripTIBLayersNull",&Cpfcan_stripTIBLayersNull_,"Cpfcan_stripTIBLayersNull_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_stripTIDLayersNull",&Cpfcan_stripTIDLayersNull_,"Cpfcan_stripTIDLayersNull_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_stripTOBLayersNull",&Cpfcan_stripTOBLayersNull_,"Cpfcan_stripTOBLayersNull_[n_Cpfcand_]/F");
    addBranch(tree,"Cpfcan_stripTECLayersNull",&Cpfcan_stripTECLayersNull_,"Cpfcan_stripTECLayersNull_[n_Cpfcand_]/F");

    //Neutral Pf candidates
    addBranch(tree,"n_Npfcand", &n_Npfcand_,"n_Npfcand_/I");
    addBranch(tree,"nNpfcand", &nNpfcand_,"nNpfcand/F");

    addBranch(tree,"Npfcan_pt", &Npfcan_pt_,"Npfcan_pt_[n_Npfcand_]/F");
    addBranch(tree,"Npfcan_eta", &Npfcan_eta_,"Npfcan_eta_[n_Npfcand_]/F");
    addBranch(tree,"Npfcan_phi", &Npfcan_phi_,"Npfcan_phi_[n_Npfcand_]/F");
    addBranch(tree,"Npfcan_ptrel", &Npfcan_ptrel_,"Npfcan_ptrel_[n_Npfcand_]/F");
    addBranch(tree,"Npfcan_e", &Npfcan_e_,"Npfcan_e_[n_Npfcand_]/F");
    addBranch(tree,"Npfcan_erel", &Npfcan_erel_,"Npfcan_erel_[n_Npfcand_]/F");

    addBranch(tree,"Npfcan_puppiw", &Npfcan_puppiw_,"Npfcan_puppiw_[n_Npfcand_]/F");


    addBranch(tree,"Npfcan_phirel",&Npfcan_phirel_,"Npfcan_phirel_[n_Npfcand_]/F");
    addBranch(tree,"Npfcan_etarel",&Npfcan_etarel_,"Npfcan_etarel_[n_Npfcand_]/F");
    addBranch(tree,"Npfcan_deltaR",&Npfcan_deltaR_,"Npfcan_deltaR_[n_Npfcand_]/F");
    addBranch(tree,"Npfcan_isGamma",&Npfcan_isGamma_,"Npfcan_isGamma_[n_Npfcand_]/F");
    addBranch(tree,"Npfcan_HadFrac",&Npfcan_HadFrac_,"Npfcan_HadFrac_[n_Npfcand_]/F");
    addBranch(tree,"Npfcan_CaloFrac",&Npfcan_CaloFrac_,"Npfcan_CaloFrac_[n_Npfcand_]/F");
    addBranch(tree,"Npfcan_drminsv",&Npfcan_drminsv_,"Npfcan_drminsv_[n_Npfcand_]/F");

	addBranch(tree,"Npfcan_pdgID",&Npfcan_pdgID_,"Npfcan_pdgID_[n_Npfcand_]/F");

}

void ntuple_pfCands::readEvent(const edm::Event& iEvent){


    n_Npfcand_=0;
    n_Cpfcand_=0;

}



//use either of these functions

//$$ bool ntuple_pfCands::fillBranches(const pat::Jet & jet, const size_t& jetidx, const  edm::View<pat::Jet> * coll){
bool ntuple_pfCands::fillBranches(const pat::Jet & jet, const size_t& jetidx, const  edm::View<pat::Jet> * coll, float EventTime){
    float etasign = 1.;
    if (jet.eta()<0) etasign =-1.;
    math::XYZVector jetDir = jet.momentum().Unit();
    GlobalVector jetRefTrackDir(jet.px(),jet.py(),jet.pz());
    const reco::Vertex & pv = vertices()->at(0);


    std::vector<sorting::sortingClass<size_t> > sortedcharged, sortedneutrals;

    const float jet_uncorr_pt=jet.correctedJet("Uncorrected").pt();
    const float jet_uncorr_e=jet.correctedJet("Uncorrected").energy();

    TrackInfoBuilder trackinfo(builder);
    //create collection first, to be able to do some sorting
    for (unsigned int i = 0; i <  jet.numberOfDaughters(); i++){
        const pat::PackedCandidate* PackedCandidate = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(i));
        if(PackedCandidate){
            if(PackedCandidate->pt() < min_candidate_pt_) continue; 
            if(PackedCandidate->charge()!=0){
                trackinfo.buildTrackInfo(PackedCandidate,jetDir,jetRefTrackDir,pv);
                sortedcharged.push_back(sorting::sortingClass<size_t>
                (i, trackinfo.getTrackSip2dSig(),
                        -mindrsvpfcand(PackedCandidate), PackedCandidate->pt()/jet_uncorr_pt));
            }
            else{
                sortedneutrals.push_back(sorting::sortingClass<size_t>
                (i, -1, -mindrsvpfcand(PackedCandidate), PackedCandidate->pt()/jet_uncorr_pt));
            }
        }
    }
	std::sort(sortedcharged.begin(),sortedcharged.end(),sorting::sortingClass<size_t>::compareByABCInv);
    n_Cpfcand_ = std::min(sortedcharged.size(),max_pfcand_);

//$$
//   std::cout << " Ntuple jet pt eta phi " << jet.pt() << " " << jet.eta() << " " << jet.phi() 
//             << " ntrks " << n_Cpfcand_ << " EventTime " << EventTime << std::endl;

    float track_time = -1;
//$$

    std::sort(sortedneutrals.begin(),sortedneutrals.end(),sorting::sortingClass<size_t>::compareByABCInv);
    std::vector<size_t> sortedchargedindices,sortedneutralsindices;
    n_Npfcand_ = std::min(sortedneutrals.size(),max_pfcand_);
	sortedchargedindices=sorting::invertSortingVector(sortedcharged);
	sortedneutralsindices=sorting::invertSortingVector(sortedneutrals);

    for (unsigned int i = 0; i <  jet.numberOfDaughters(); i++){
        const pat::PackedCandidate* PackedCandidate_ = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(i));
        //const auto& PackedCandidate_=s.get();
        if(!PackedCandidate_) continue;
        if(PackedCandidate_->pt() < min_candidate_pt_) continue; 

        // get the dr with the closest sv
        float drminpfcandsv_ = mindrsvpfcand(PackedCandidate_);


        /// This might include more than PF candidates, e.g. Reco muons and could
        /// be double counting. Needs to be checked.!!!!
        ///
        /// Split to charged and neutral candidates
        if(PackedCandidate_->charge()!=0 ){

            size_t fillntupleentry= sortedchargedindices.at(i);
            if(fillntupleentry>=max_pfcand_) continue;


            Cpfcan_pt_[fillntupleentry] = PackedCandidate_->pt();
            Cpfcan_eta_[fillntupleentry] = PackedCandidate_->eta();
            Cpfcan_phi_[fillntupleentry] = PackedCandidate_->phi();
            Cpfcan_ptrel_[fillntupleentry] = catchInfsAndBound(PackedCandidate_->pt()/jet_uncorr_pt,0,-1,0,-1);
            Cpfcan_erel_[fillntupleentry] = catchInfsAndBound(PackedCandidate_->energy()/jet_uncorr_e,0,-1,0,-1);
            Cpfcan_e_[fillntupleentry] = PackedCandidate_->energy();
            Cpfcan_phirel_[fillntupleentry] = catchInfsAndBound(fabs(reco::deltaPhi(PackedCandidate_->phi(),jet.phi())),0,-2,0,-0.5);
            Cpfcan_etarel_[fillntupleentry] = catchInfsAndBound(fabs(PackedCandidate_->eta()-jet.eta()),0,-2,0,-0.5);
            Cpfcan_deltaR_[fillntupleentry] =catchInfsAndBound(reco::deltaR(*PackedCandidate_,jet),0,-0.6,0,-0.6);
            Cpfcan_dxy_[fillntupleentry] = catchInfsAndBound(fabs(PackedCandidate_->dxy()),0,-50,50);

            Cpfcan_dxyerrinv_[fillntupleentry]= PackedCandidate_->hasTrackDetails() ? catchInfsAndBound(1/PackedCandidate_->dxyError(),0,-1, 10000.) : -1;

            Cpfcan_dxysig_[fillntupleentry]= PackedCandidate_->hasTrackDetails() ? catchInfsAndBound(fabs(PackedCandidate_->dxy()/PackedCandidate_->dxyError()),0.,-2000,2000) : 0.;


            Cpfcan_dz_[fillntupleentry] = PackedCandidate_->dz();
            Cpfcan_VTX_ass_[fillntupleentry] = PackedCandidate_->pvAssociationQuality();

            Cpfcan_fromPV_[fillntupleentry] = PackedCandidate_->fromPV();

            float tempdontopt=PackedCandidate_->vx();
            tempdontopt++;

            Cpfcan_vertexChi2_[fillntupleentry]=PackedCandidate_->vertexChi2();
            Cpfcan_vertexNdof_[fillntupleentry]=PackedCandidate_->vertexNdof();
            
			Cpfcan_CaloFrac_[fillntupleentry] = PackedCandidate_->caloFraction();
            Cpfcan_HadFrac_[fillntupleentry] = PackedCandidate_->hcalFraction();

            //divided
            Cpfcan_vertexNormalizedChi2_[fillntupleentry]=PackedCandidate_->vertexNormalizedChi2();
            Cpfcan_vertex_rho_[fillntupleentry]=catchInfsAndBound(PackedCandidate_->vertex().rho(),0,-1,50);
            Cpfcan_vertex_phirel_[fillntupleentry]=reco::deltaPhi(PackedCandidate_->vertex().phi(),jet.phi());
            Cpfcan_vertex_etarel_[fillntupleentry]=etasign*(PackedCandidate_->vertex().eta()-jet.eta());
            Cpfcan_vertexRef_mass_[fillntupleentry]=PackedCandidate_->vertexRef()->p4().M();


            Cpfcan_puppiw_[fillntupleentry] = PackedCandidate_->puppiWeight();


            /*
            reco::Track::CovarianceMatrix myCov = PseudoTrack.covariance ();
            //https://github.com/cms-sw/cmssw/blob/CMSSW_9_0_X/DataFormats/PatCandidates/interface/PackedCandidate.h#L394
            Cpfcan_dptdpt_[fillntupleentry] =    catchInfsAndBound(myCov[0][0],0,-1,1);
            Cpfcan_detadeta_[fillntupleentry]=   catchInfsAndBound(myCov[1][1],0,-1,0.01);
            Cpfcan_dphidphi_[fillntupleentry]=   catchInfsAndBound(myCov[2][2],0,-1,0.1);
            Cpfcan_dxydxy_[fillntupleentry] =    catchInfsAndBound(myCov[3][3],7.,-1,7); //zero if pvAssociationQuality ==7 ?
            Cpfcan_dzdz_[fillntupleentry] =      catchInfsAndBound(myCov[4][4],6.5,-1,6.5); //zero if pvAssociationQuality ==7 ?
            Cpfcan_dxydz_[fillntupleentry] =     catchInfsAndBound(myCov[3][4],6.,-6,6); //zero if pvAssociationQuality ==7 ?
            Cpfcan_dphidxy_[fillntupleentry] =   catchInfs(myCov[2][3],-0.03); //zero if pvAssociationQuality ==7 ?
            Cpfcan_dlambdadz_[fillntupleentry]=  catchInfs(myCov[1][4],-0.03); //zero if pvAssociationQuality ==7 ?
             */

            trackinfo.buildTrackInfo(PackedCandidate_,jetDir,jetRefTrackDir,pv);

            Cpfcan_BtagPf_trackMomentum_[fillntupleentry]   =catchInfsAndBound(trackinfo.getTrackMomentum(),0,0 ,1000);
            Cpfcan_BtagPf_trackEta_[fillntupleentry]        =catchInfsAndBound(trackinfo.getTrackEta()   ,  0,-5,5);
            Cpfcan_BtagPf_trackEtaRel_[fillntupleentry]     =catchInfsAndBound(trackinfo.getTrackEtaRel(),  0,-5,15);
            Cpfcan_BtagPf_trackPtRel_[fillntupleentry]      =catchInfsAndBound(trackinfo.getTrackPtRel(),   0,-1,4);
            Cpfcan_BtagPf_trackPPar_[fillntupleentry]       =catchInfsAndBound(trackinfo.getTrackPPar(),    0,-1e5,1e5 );
            Cpfcan_BtagPf_trackDeltaR_[fillntupleentry]     =catchInfsAndBound(trackinfo.getTrackDeltaR(),  0,-5,5 );
            Cpfcan_BtagPf_trackPtRatio_[fillntupleentry]    =catchInfsAndBound(trackinfo.getTrackPtRatio(), 0,-1,10 );
            Cpfcan_BtagPf_trackPParRatio_[fillntupleentry]  =catchInfsAndBound(trackinfo.getTrackPParRatio(),0,-10,100);
            Cpfcan_BtagPf_trackSip3dVal_[fillntupleentry]   =catchInfsAndBound(trackinfo.getTrackSip3dVal(), 0, -1,1e5 );
            Cpfcan_BtagPf_trackSip3dSig_[fillntupleentry]   =catchInfsAndBound(trackinfo.getTrackSip3dSig(), 0, -1,4e4 );
            Cpfcan_BtagPf_trackSip2dVal_[fillntupleentry]   =catchInfsAndBound(trackinfo.getTrackSip2dVal(), 0, -1,70 );
            Cpfcan_BtagPf_trackSip2dSig_[fillntupleentry]   =catchInfsAndBound(trackinfo.getTrackSip2dSig(), 0, -1,4e4 );
            Cpfcan_BtagPf_trackDecayLen_[fillntupleentry]   =0;
            Cpfcan_BtagPf_trackJetDistVal_[fillntupleentry] =catchInfsAndBound(trackinfo.getTrackJetDistVal(),0,-20,1 );
            Cpfcan_BtagPf_trackJetDistSig_[fillntupleentry] =catchInfsAndBound(trackinfo.getTrackJetDistSig(),0,-1,1e5 );

//$$
	    track_time = -1.;
            auto track = PackedCandidate_->bestTrack();
            if ( track && EventTime != -1 ) {
              if ( track->covt0t0() > 0. && abs(track->t0()) < 1 ) {
	        track_time = track->t0() - EventTime;
	        track_time = TMath::Abs(track_time);
	      }
	    }
            Cpfcan_time_[fillntupleentry] = track_time;

//   std::cout << " Ntuple Cpfcan " << fillntupleentry << " pt " << PackedCandidate_->pt()
//             << " IPsig " << trackinfo.getTrackSip2dSig()   
//             << " time " << track_time 
// 	       << std::endl;     
//$$

            // TO DO: we can do better than that by including reco::muon informations
            Cpfcan_isMu_[fillntupleentry] = 0;
            if(abs(PackedCandidate_->pdgId())==13) {
                Cpfcan_isMu_[fillntupleentry] = 1;
            }
            // TO DO: we can do better than that by including reco::electron informations
            Cpfcan_isEl_[fillntupleentry] = 0;
            if(abs(PackedCandidate_->pdgId())==11) {
                Cpfcan_isEl_[fillntupleentry] = 1;

            }

            Cpfcan_charge_[fillntupleentry] = PackedCandidate_->charge();
            Cpfcan_lostInnerHits_[fillntupleentry] = catchInfs(PackedCandidate_->lostInnerHits(),2);
			Cpfcan_numberOfPixelHits_[fillntupleentry] = catchInfs(PackedCandidate_->numberOfPixelHits(),-1);

	    //std::cout << PackedCandidate_->lostInnerHits()<< " inner hits " <<std::endl;
	    //std::cout << PackedCandidate_->numberOfPixelHits()<< " Pixel hits + masked " <<std::endl;
	    //std::cout <<PackedCandidate_->pixelLayersWithMeasurement()<< " Pixel hits " <<std::endl;

			Cpfcan_chi2_[fillntupleentry] = PackedCandidate_->hasTrackDetails() ? \
			catchInfsAndBound(PackedCandidate_->pseudoTrack().normalizedChi2(),300,-1,300) : -1;
			//for some reason this returns the quality enum not a mask.
			Cpfcan_quality_[fillntupleentry] = PackedCandidate_->hasTrackDetails() ? 
			PackedCandidate_->pseudoTrack().qualityMask() : (1 << reco::TrackBase::loose);

            Cpfcan_drminsv_[fillntupleentry] = catchInfsAndBound(drminpfcandsv_,0,-0.4,0,-0.4);

          //hit pattern variables, as defined here https://github.com/cms-sw/cmssw/blob/master/DataFormats/TrackReco/interface/HitPattern.h
            //get track associated to a jet constituent
            const reco::Track *track_ptr = nullptr;
            auto pf_candidate = dynamic_cast<const reco::PFCandidate *>(PackedCandidate_);
            auto packed_candidate = dynamic_cast<const pat::PackedCandidate *>(PackedCandidate_);
            if(pf_candidate){
             track_ptr = pf_candidate->bestTrack(); //trackRef was sometimes null
            }else if(packed_candidate && packed_candidate->hasTrackDetails()){//if PackedCandidate does not have TrackDetails this gives an Exception because unpackCovariance might be called for pseudoTrack/bestTrack
             track_ptr = &(packed_candidate->pseudoTrack());
            }
            //get hit pattern information
            if(track_ptr){
             const reco::HitPattern &p = track_ptr->hitPattern();
             //Tracker per layer
             //Pixel barrel 
             int Cpfcan_nhitpixelBarrelLayer1 = 0;
             int Cpfcan_nhitpixelBarrelLayer2 = 0;
             int Cpfcan_nhitpixelBarrelLayer3 = 0;
             int Cpfcan_nhitpixelBarrelLayer4 = 0;
             //Pixel Endcap 
             int Cpfcan_nhitpixelEndcapLayer1 = 0;
             int Cpfcan_nhitpixelEndcapLayer2 = 0;
             //Strip TIB
             int Cpfcan_nhitstripTIBLayer1 = 0;
             int Cpfcan_nhitstripTIBLayer2 = 0;
             int Cpfcan_nhitstripTIBLayer3 = 0;
             int Cpfcan_nhitstripTIBLayer4 = 0;
             //Strip TID
             int Cpfcan_nhitstripTIDLayer1 = 0;
             int Cpfcan_nhitstripTIDLayer2 = 0;
             int Cpfcan_nhitstripTIDLayer3 = 0;
             //Strip TOB
             int Cpfcan_nhitstripTOBLayer1 = 0;
             int Cpfcan_nhitstripTOBLayer2 = 0;
             int Cpfcan_nhitstripTOBLayer3 = 0;
             int Cpfcan_nhitstripTOBLayer4 = 0;
             int Cpfcan_nhitstripTOBLayer5 = 0;
             int Cpfcan_nhitstripTOBLayer6 = 0;
             //Strip TEC
             int Cpfcan_nhitstripTECLayer1 = 0;
             int Cpfcan_nhitstripTECLayer2 = 0;
             int Cpfcan_nhitstripTECLayer3 = 0;
             int Cpfcan_nhitstripTECLayer4 = 0;
             int Cpfcan_nhitstripTECLayer5 = 0;
             int Cpfcan_nhitstripTECLayer6 = 0;
             int Cpfcan_nhitstripTECLayer7 = 0;
             int Cpfcan_nhitstripTECLayer8 = 0;
             int Cpfcan_nhitstripTECLayer9 = 0;
             // loop over the hits of the track.
             //const static unsigned short LayerOffset = 3;
             //const static unsigned short LayerMask = 0xF;
             for(int nh = 0; nh < p.numberOfAllHits(reco::HitPattern::TRACK_HITS); nh++){
              uint32_t hit = p.getHitPattern(reco::HitPattern::TRACK_HITS, nh);
              if(p.validHitFilter(hit)){// if the hit is valid
               //Pixel Barrel // it is in pixel barrel
               if(p.pixelBarrelHitFilter(hit)){
                //std::cout << "valid hit found in pixel Barrel layer " << p.getLayer(hit) << std::endl;
                //if(p.getLayer(hit)==1){
                // std::cout<< (hit >> LayerOffset) << " " << ((hit >> LayerOffset) & LayerMask) << std::endl;
                //}  
                if(p.getLayer(hit)==1) Cpfcan_nhitpixelBarrelLayer1 = Cpfcan_nhitpixelBarrelLayer1+1;
                if(p.getLayer(hit)==2) Cpfcan_nhitpixelBarrelLayer2 = Cpfcan_nhitpixelBarrelLayer2+1;
                if(p.getLayer(hit)==3) Cpfcan_nhitpixelBarrelLayer3 = Cpfcan_nhitpixelBarrelLayer3+1;
                if(p.getLayer(hit)==4) Cpfcan_nhitpixelBarrelLayer4 = Cpfcan_nhitpixelBarrelLayer4+1;
               } 
               //Pixel Endcap
               if(p.pixelEndcapHitFilter(hit)){
                //std::cout << "valid hit found in pixel Endcap layer " << p.getLayer(hit) << std::endl;
                if(p.getLayer(hit)==1) Cpfcan_nhitpixelEndcapLayer1 = Cpfcan_nhitpixelEndcapLayer1+1;
                if(p.getLayer(hit)==2) Cpfcan_nhitpixelEndcapLayer2 = Cpfcan_nhitpixelEndcapLayer2+1;
               } 
               //Strip TIB
               if(p.stripTIBHitFilter(hit)){
                //std::cout << "valid hit found in TIB layer " << p.getLayer(hit) << std::endl;
                if(p.getLayer(hit)==1) Cpfcan_nhitstripTIBLayer1 = Cpfcan_nhitstripTIBLayer1+1;
                if(p.getLayer(hit)==2) Cpfcan_nhitstripTIBLayer2 = Cpfcan_nhitstripTIBLayer2+1;
                if(p.getLayer(hit)==3) Cpfcan_nhitstripTIBLayer3 = Cpfcan_nhitstripTIBLayer3+1;
                if(p.getLayer(hit)==4) Cpfcan_nhitstripTIBLayer4 = Cpfcan_nhitstripTIBLayer4+1;
               } 
               //Strip TID
               if(p.stripTIDHitFilter(hit)){
                //std::cout << "valid hit found in TID layer " << p.getLayer(hit) << std::endl;
                if(p.getLayer(hit)==1) Cpfcan_nhitstripTIDLayer1 = Cpfcan_nhitstripTIDLayer1+1;
                if(p.getLayer(hit)==2) Cpfcan_nhitstripTIDLayer2 = Cpfcan_nhitstripTIDLayer2+1;
                if(p.getLayer(hit)==3) Cpfcan_nhitstripTIDLayer3 = Cpfcan_nhitstripTIDLayer3+1;
               } 
               //Strip TOB
               if(p.stripTOBHitFilter(hit)){
                //std::cout << "valid hit found in TOB layer " << p.getLayer(hit) << std::endl;
                if(p.getLayer(hit)==1) Cpfcan_nhitstripTOBLayer1 = Cpfcan_nhitstripTOBLayer1+1;
                if(p.getLayer(hit)==2) Cpfcan_nhitstripTOBLayer2 = Cpfcan_nhitstripTOBLayer2+1;
                if(p.getLayer(hit)==3) Cpfcan_nhitstripTOBLayer3 = Cpfcan_nhitstripTOBLayer3+1;
                if(p.getLayer(hit)==4) Cpfcan_nhitstripTOBLayer4 = Cpfcan_nhitstripTOBLayer4+1;
                if(p.getLayer(hit)==5) Cpfcan_nhitstripTOBLayer5 = Cpfcan_nhitstripTOBLayer5+1;
                if(p.getLayer(hit)==6) Cpfcan_nhitstripTOBLayer6 = Cpfcan_nhitstripTOBLayer6+1;
               } 
               //Strip TEC
               if(p.stripTECHitFilter(hit)){
                //std::cout << "valid hit found in TEC layer " << p.getLayer(hit) << std::endl;
                if(p.getLayer(hit)==1) Cpfcan_nhitstripTECLayer1 = Cpfcan_nhitstripTECLayer1+1;
                if(p.getLayer(hit)==2) Cpfcan_nhitstripTECLayer2 = Cpfcan_nhitstripTECLayer2+1;
                if(p.getLayer(hit)==3) Cpfcan_nhitstripTECLayer3 = Cpfcan_nhitstripTECLayer3+1;
                if(p.getLayer(hit)==4) Cpfcan_nhitstripTECLayer4 = Cpfcan_nhitstripTECLayer4+1;
                if(p.getLayer(hit)==5) Cpfcan_nhitstripTECLayer5 = Cpfcan_nhitstripTECLayer5+1;
                if(p.getLayer(hit)==6) Cpfcan_nhitstripTECLayer6 = Cpfcan_nhitstripTECLayer6+1;
                if(p.getLayer(hit)==7) Cpfcan_nhitstripTECLayer7 = Cpfcan_nhitstripTECLayer7+1;
                if(p.getLayer(hit)==8) Cpfcan_nhitstripTECLayer8 = Cpfcan_nhitstripTECLayer8+1;
                if(p.getLayer(hit)==9) Cpfcan_nhitstripTECLayer9 = Cpfcan_nhitstripTECLayer9+1;
               } 
              }
              //// expert level: printout the hit in 11-bit binary format
              //std::cout << "hit in 11-bit binary format = ";
              //for (int j = 10; j >= 0; j--){
              //    int bit = (hit >> j) & 0x1;
              //    std::cout << bit;
              //}
              //std::cout << std::endl;
             }
             //Pixel Barrel 
             Cpfcan_nhitpixelBarrelLayer1_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidPixelBarrelHits()) ? catchInfsAndBound(Cpfcan_nhitpixelBarrelLayer1,-1,0,100,0) : 0;
             Cpfcan_nhitpixelBarrelLayer2_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidPixelBarrelHits()) ? catchInfsAndBound(Cpfcan_nhitpixelBarrelLayer2,-1,0,100,0) : 0;
             Cpfcan_nhitpixelBarrelLayer3_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidPixelBarrelHits()) ? catchInfsAndBound(Cpfcan_nhitpixelBarrelLayer3,-1,0,100,0) : 0;
             Cpfcan_nhitpixelBarrelLayer4_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidPixelBarrelHits()) ? catchInfsAndBound(Cpfcan_nhitpixelBarrelLayer4,-1,0,100,0) : 0;
             //Pixel Endcap 
             Cpfcan_nhitpixelEndcapLayer1_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidPixelEndcapHits()) ? catchInfsAndBound(Cpfcan_nhitpixelEndcapLayer1,-1,0,100,0) : 0;
             Cpfcan_nhitpixelEndcapLayer2_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidPixelEndcapHits()) ? catchInfsAndBound(Cpfcan_nhitpixelEndcapLayer2,-1,0,100,0) : 0;
             //Strip TIB
             Cpfcan_nhitstripTIBLayer1_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripTIBHits()) ? catchInfsAndBound(Cpfcan_nhitstripTIBLayer1,-1,0,100,0) : 0;
             Cpfcan_nhitstripTIBLayer2_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripTIBHits()) ? catchInfsAndBound(Cpfcan_nhitstripTIBLayer2,-1,0,100,0) : 0;
             Cpfcan_nhitstripTIBLayer3_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripTIBHits()) ? catchInfsAndBound(Cpfcan_nhitstripTIBLayer3,-1,0,100,0) : 0;
             Cpfcan_nhitstripTIBLayer4_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripTIBHits()) ? catchInfsAndBound(Cpfcan_nhitstripTIBLayer4,-1,0,100,0) : 0;
             //Strip TID
             Cpfcan_nhitstripTIDLayer1_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripTIDHits()) ? catchInfsAndBound(Cpfcan_nhitstripTIDLayer1,-1,0,100,0) : 0;
             Cpfcan_nhitstripTIDLayer2_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripTIDHits()) ? catchInfsAndBound(Cpfcan_nhitstripTIDLayer2,-1,0,100,0) : 0;
             Cpfcan_nhitstripTIDLayer3_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripTIDHits()) ? catchInfsAndBound(Cpfcan_nhitstripTIDLayer3,-1,0,100,0) : 0;
             //Strip TOB
             Cpfcan_nhitstripTOBLayer1_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripTOBHits()) ? catchInfsAndBound(Cpfcan_nhitstripTOBLayer1,-1,0,100,0) : 0;
             Cpfcan_nhitstripTOBLayer2_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripTOBHits()) ? catchInfsAndBound(Cpfcan_nhitstripTOBLayer2,-1,0,100,0) : 0;
             Cpfcan_nhitstripTOBLayer3_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripTOBHits()) ? catchInfsAndBound(Cpfcan_nhitstripTOBLayer3,-1,0,100,0) : 0;
             Cpfcan_nhitstripTOBLayer4_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripTOBHits()) ? catchInfsAndBound(Cpfcan_nhitstripTOBLayer4,-1,0,100,0) : 0;
             Cpfcan_nhitstripTOBLayer5_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripTOBHits()) ? catchInfsAndBound(Cpfcan_nhitstripTOBLayer5,-1,0,100,0) : 0;
             Cpfcan_nhitstripTOBLayer6_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripTOBHits()) ? catchInfsAndBound(Cpfcan_nhitstripTOBLayer6,-1,0,100,0) : 0;
             //Strip TEC
             Cpfcan_nhitstripTECLayer1_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripTECHits()) ? catchInfsAndBound(Cpfcan_nhitstripTECLayer1,-1,0,100,0) : 0;
             Cpfcan_nhitstripTECLayer2_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripTECHits()) ? catchInfsAndBound(Cpfcan_nhitstripTECLayer2,-1,0,100,0) : 0;
             Cpfcan_nhitstripTECLayer3_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripTECHits()) ? catchInfsAndBound(Cpfcan_nhitstripTECLayer3,-1,0,100,0) : 0;
             Cpfcan_nhitstripTECLayer4_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripTECHits()) ? catchInfsAndBound(Cpfcan_nhitstripTECLayer4,-1,0,100,0) : 0;
             Cpfcan_nhitstripTECLayer5_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripTECHits()) ? catchInfsAndBound(Cpfcan_nhitstripTECLayer5,-1,0,100,0) : 0;
             Cpfcan_nhitstripTECLayer6_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripTECHits()) ? catchInfsAndBound(Cpfcan_nhitstripTECLayer6,-1,0,100,0) : 0;
             Cpfcan_nhitstripTECLayer7_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripTECHits()) ? catchInfsAndBound(Cpfcan_nhitstripTECLayer7,-1,0,100,0) : 0;
             Cpfcan_nhitstripTECLayer8_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripTECHits()) ? catchInfsAndBound(Cpfcan_nhitstripTECLayer8,-1,0,100,0) : 0;
             Cpfcan_nhitstripTECLayer9_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripTECHits()) ? catchInfsAndBound(Cpfcan_nhitstripTECLayer9,-1,0,100,0) : 0;
             //Tracker all layers together   
             //Valid hits
             Cpfcan_numberOfValidHits_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidHits()) ? catchInfsAndBound(track_ptr->hitPattern().numberOfValidHits(),-1,0,100,0) : 0;
             Cpfcan_numberOfValidTrackerHits_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidTrackerHits()) ? catchInfsAndBound(track_ptr->hitPattern().numberOfValidTrackerHits(),-1,0,100,0) : 0;
             Cpfcan_numberOfValidPixelHits_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidPixelHits()) ? catchInfsAndBound(track_ptr->hitPattern().numberOfValidPixelHits(),-1,0,100,0) : 0;
             Cpfcan_numberOfValidPixelBarrelHits_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidPixelBarrelHits()) ? catchInfsAndBound(track_ptr->hitPattern().numberOfValidPixelBarrelHits(),-1,0,100,0) : 0;
             Cpfcan_numberOfValidPixelEndcapHits_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidPixelEndcapHits()) ? catchInfsAndBound(track_ptr->hitPattern().numberOfValidPixelEndcapHits(),-1,0,100,0) : 0;
             Cpfcan_numberOfValidStripHits_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripHits()) ? catchInfsAndBound(track_ptr->hitPattern().numberOfValidStripHits(),-1,0,100,0) : -1;
             Cpfcan_numberOfValidStripTIBHits_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripTIBHits()) ? catchInfsAndBound(track_ptr->hitPattern().numberOfValidStripTIBHits(),-1,0,100,0) : 0;
             Cpfcan_numberOfValidStripTIDHits_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripTIDHits()) ? catchInfsAndBound(track_ptr->hitPattern().numberOfValidStripTIDHits(),-1,0,100,0) : 0;
             Cpfcan_numberOfValidStripTOBHits_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripTOBHits()) ? catchInfsAndBound(track_ptr->hitPattern().numberOfValidStripTOBHits(),-1,0,100,0) : 0;
             Cpfcan_numberOfValidStripTECHits_[fillntupleentry] = (track_ptr->hitPattern().numberOfValidStripTECHits()) ? catchInfsAndBound(track_ptr->hitPattern().numberOfValidStripTECHits(),-1,0,100,0) : 0;
             //LayersWithMeasurement
             Cpfcan_trackerLayersWithMeasurementOld_[fillntupleentry] = (track_ptr->hitPattern().trackerLayersWithMeasurementOld()) ? catchInfsAndBound(track_ptr->hitPattern().trackerLayersWithMeasurementOld(),-1,0,100,0) : 0;
             Cpfcan_trackerLayersWithMeasurement_[fillntupleentry] = (track_ptr->hitPattern().trackerLayersWithMeasurement()) ? catchInfsAndBound(track_ptr->hitPattern().trackerLayersWithMeasurement(),-1,0,100,0) : 0;
             Cpfcan_pixelLayersWithMeasurementOld_[fillntupleentry] = (track_ptr->hitPattern().pixelLayersWithMeasurementOld()) ? catchInfsAndBound(track_ptr->hitPattern().pixelLayersWithMeasurementOld(),-1,0,100,0) : 0;
             Cpfcan_pixelLayersWithMeasurement_[fillntupleentry] = (track_ptr->hitPattern().pixelLayersWithMeasurement()) ? catchInfsAndBound(track_ptr->hitPattern().pixelLayersWithMeasurement(),-1,0,100,0) : 0;
             Cpfcan_stripLayersWithMeasurement_[fillntupleentry] = (track_ptr->hitPattern().stripLayersWithMeasurement()) ? catchInfsAndBound(track_ptr->hitPattern().stripLayersWithMeasurement(),-1,0,100,0) : 0;
             Cpfcan_pixelBarrelLayersWithMeasurement_[fillntupleentry] = (track_ptr->hitPattern().pixelBarrelLayersWithMeasurement()) ? catchInfsAndBound(track_ptr->hitPattern().pixelBarrelLayersWithMeasurement(),-1,0,100,0) : 0;
             Cpfcan_pixelEndcapLayersWithMeasurement_[fillntupleentry] = (track_ptr->hitPattern().pixelEndcapLayersWithMeasurement()) ? catchInfsAndBound(track_ptr->hitPattern().pixelEndcapLayersWithMeasurement(),-1,0,100,0) : 0;
             Cpfcan_stripTIBLayersWithMeasurement_[fillntupleentry] = (track_ptr->hitPattern().stripTIBLayersWithMeasurement()) ? catchInfsAndBound(track_ptr->hitPattern().stripTIBLayersWithMeasurement(),-1,0,100,0) : 0;
             Cpfcan_stripTIDLayersWithMeasurement_[fillntupleentry] = (track_ptr->hitPattern().stripTIDLayersWithMeasurement()) ? catchInfsAndBound(track_ptr->hitPattern().stripTIDLayersWithMeasurement(),-1,0,100,0) : 0;
             Cpfcan_stripTOBLayersWithMeasurement_[fillntupleentry] = (track_ptr->hitPattern().stripTOBLayersWithMeasurement()) ? catchInfsAndBound(track_ptr->hitPattern().stripTOBLayersWithMeasurement(),-1,0,100,0) : 0;
             Cpfcan_stripTECLayersWithMeasurement_[fillntupleentry] = (track_ptr->hitPattern().stripTECLayersWithMeasurement()) ? catchInfsAndBound(track_ptr->hitPattern().stripTECLayersWithMeasurement(),-1,0,100,0) : 0;
             //Null
             Cpfcan_trackerLayersNull_[fillntupleentry] = (track_ptr->hitPattern().trackerLayersNull()) ? catchInfsAndBound(track_ptr->hitPattern().trackerLayersNull(),-1,0,100,0) : 0;
             Cpfcan_pixelLayersNull_[fillntupleentry] = (track_ptr->hitPattern().pixelLayersNull()) ? catchInfsAndBound(track_ptr->hitPattern().pixelLayersNull(),-1,0,100,0) : 0;
             Cpfcan_stripLayersNull_[fillntupleentry] = (track_ptr->hitPattern().stripLayersNull()) ? catchInfsAndBound(track_ptr->hitPattern().stripLayersNull(),-1,0,100,0) : 0;
             Cpfcan_pixelBarrelLayersNull_[fillntupleentry] = (track_ptr->hitPattern().pixelBarrelLayersNull()) ? catchInfsAndBound(track_ptr->hitPattern().pixelBarrelLayersNull(),-1,0,100,0) : 0;
             Cpfcan_pixelEndcapLayersNull_[fillntupleentry] = (track_ptr->hitPattern().pixelEndcapLayersNull()) ? catchInfsAndBound(track_ptr->hitPattern().pixelEndcapLayersNull(),-1,0,100,0) : 0;
             Cpfcan_stripTIBLayersNull_[fillntupleentry] = (track_ptr->hitPattern().stripTIBLayersNull()) ? catchInfsAndBound(track_ptr->hitPattern().stripTIBLayersNull(),-1,0,100,0) : 0;
             Cpfcan_stripTIDLayersNull_[fillntupleentry] = (track_ptr->hitPattern().stripTIDLayersNull()) ? catchInfsAndBound(track_ptr->hitPattern().stripTIDLayersNull(),-1,0,100,0) : 0;
             Cpfcan_stripTOBLayersNull_[fillntupleentry] = (track_ptr->hitPattern().stripTOBLayersNull()) ? catchInfsAndBound(track_ptr->hitPattern().stripTOBLayersNull(),-1,0,100,0) : 0;
             Cpfcan_stripTECLayersNull_[fillntupleentry] = (track_ptr->hitPattern().stripTECLayersNull()) ? catchInfsAndBound(track_ptr->hitPattern().stripTECLayersNull(),-1,0,100,0) : 0;
            }else{
             //Tracker per layer
             //Pixel barrel 
             Cpfcan_nhitpixelBarrelLayer1_[fillntupleentry] = 0;
             Cpfcan_nhitpixelBarrelLayer2_[fillntupleentry] = 0;
             Cpfcan_nhitpixelBarrelLayer3_[fillntupleentry] = 0;
             Cpfcan_nhitpixelBarrelLayer4_[fillntupleentry] = 0;
             //Pixel Endcap 
             Cpfcan_nhitpixelEndcapLayer1_[fillntupleentry] = 0;
             Cpfcan_nhitpixelEndcapLayer2_[fillntupleentry] = 0;
             //Strip TIB
             Cpfcan_nhitstripTIBLayer1_[fillntupleentry] = 0;
             Cpfcan_nhitstripTIBLayer2_[fillntupleentry] = 0;
             Cpfcan_nhitstripTIBLayer3_[fillntupleentry] = 0;
             Cpfcan_nhitstripTIBLayer4_[fillntupleentry] = 0;
             //Strip TID
             Cpfcan_nhitstripTIDLayer1_[fillntupleentry] = 0;
             Cpfcan_nhitstripTIDLayer2_[fillntupleentry] = 0;
             Cpfcan_nhitstripTIDLayer3_[fillntupleentry] = 0;
             //Strip TOB
             Cpfcan_nhitstripTOBLayer1_[fillntupleentry] = 0;
             Cpfcan_nhitstripTOBLayer2_[fillntupleentry] = 0;
             Cpfcan_nhitstripTOBLayer3_[fillntupleentry] = 0;
             Cpfcan_nhitstripTOBLayer4_[fillntupleentry] = 0;
             Cpfcan_nhitstripTOBLayer5_[fillntupleentry] = 0;
             Cpfcan_nhitstripTOBLayer6_[fillntupleentry] = 0;
             //Strip TEC
             Cpfcan_nhitstripTECLayer1_[fillntupleentry] = 0;
             Cpfcan_nhitstripTECLayer2_[fillntupleentry] = 0;
             Cpfcan_nhitstripTECLayer3_[fillntupleentry] = 0;
             Cpfcan_nhitstripTECLayer4_[fillntupleentry] = 0;
             Cpfcan_nhitstripTECLayer5_[fillntupleentry] = 0;
             Cpfcan_nhitstripTECLayer6_[fillntupleentry] = 0;
             Cpfcan_nhitstripTECLayer7_[fillntupleentry] = 0;
             Cpfcan_nhitstripTECLayer8_[fillntupleentry] = 0;
             Cpfcan_nhitstripTECLayer9_[fillntupleentry] = 0;
             //Tracker all layers together
             //Valid hits
             Cpfcan_numberOfValidHits_[fillntupleentry] = 0; 
             Cpfcan_numberOfValidTrackerHits_[fillntupleentry] = 0; 
             Cpfcan_numberOfValidPixelHits_[fillntupleentry] = 0; 
             Cpfcan_numberOfValidPixelBarrelHits_[fillntupleentry] = 0; 
             Cpfcan_numberOfValidPixelEndcapHits_[fillntupleentry] = 0; 
             Cpfcan_numberOfValidStripHits_[fillntupleentry] = 0; 
             Cpfcan_numberOfValidStripTIBHits_[fillntupleentry] = 0; 
             Cpfcan_numberOfValidStripTIDHits_[fillntupleentry] = 0; 
             Cpfcan_numberOfValidStripTOBHits_[fillntupleentry] = 0; 
             Cpfcan_numberOfValidStripTECHits_[fillntupleentry] = 0; 
             //LayersWithMeasurement
             Cpfcan_trackerLayersWithMeasurementOld_[fillntupleentry] = 0; 
             Cpfcan_trackerLayersWithMeasurement_[fillntupleentry] = 0; 
             Cpfcan_pixelLayersWithMeasurementOld_[fillntupleentry] = 0; 
             Cpfcan_pixelLayersWithMeasurement_[fillntupleentry] = 0; 
             Cpfcan_stripLayersWithMeasurement_[fillntupleentry] = 0; 
             Cpfcan_pixelBarrelLayersWithMeasurement_[fillntupleentry] = 0; 
             Cpfcan_pixelEndcapLayersWithMeasurement_[fillntupleentry] = 0; 
             Cpfcan_stripTIBLayersWithMeasurement_[fillntupleentry] = 0; 
             Cpfcan_stripTIDLayersWithMeasurement_[fillntupleentry] = 0; 
             Cpfcan_stripTOBLayersWithMeasurement_[fillntupleentry] = 0; 
             Cpfcan_stripTECLayersWithMeasurement_[fillntupleentry] = 0; 
             //Null
             Cpfcan_trackerLayersNull_[fillntupleentry] = 0;
             Cpfcan_pixelLayersNull_[fillntupleentry] = 0; 
             Cpfcan_stripLayersNull_[fillntupleentry] = 0; 
             Cpfcan_pixelBarrelLayersNull_[fillntupleentry] = 0;
             Cpfcan_pixelEndcapLayersNull_[fillntupleentry] = 0;
             Cpfcan_stripTIBLayersNull_[fillntupleentry] = 0;
             Cpfcan_stripTIDLayersNull_[fillntupleentry] = 0;
             Cpfcan_stripTOBLayersNull_[fillntupleentry] = 0;
             Cpfcan_stripTECLayersNull_[fillntupleentry] = 0;
			}
        }
        else{// neutral candidates


            size_t fillntupleentry= sortedneutralsindices.at(i);
            if(fillntupleentry>=max_pfcand_) continue;

            Npfcan_pt_[fillntupleentry] = PackedCandidate_->pt();
            Npfcan_eta_[fillntupleentry] = PackedCandidate_->eta();
            Npfcan_phi_[fillntupleentry] = PackedCandidate_->phi();
            Npfcan_ptrel_[fillntupleentry] = catchInfsAndBound(PackedCandidate_->pt()/jet_uncorr_pt,0,-1,0,-1);
            Npfcan_erel_[fillntupleentry] = catchInfsAndBound(PackedCandidate_->energy()/jet_uncorr_e,0,-1,0,-1);
            Npfcan_e_[fillntupleentry] = PackedCandidate_->energy();
            Npfcan_puppiw_[fillntupleentry] = PackedCandidate_->puppiWeight();
            Npfcan_phirel_[fillntupleentry] = catchInfsAndBound(fabs(reco::deltaPhi(PackedCandidate_->phi(),jet.phi())),0,-2,0,-0.5);
            Npfcan_etarel_[fillntupleentry] = catchInfsAndBound(fabs(PackedCandidate_->eta()-jet.eta()),0,-2,0,-0.5);
            Npfcan_deltaR_[fillntupleentry] = catchInfsAndBound(reco::deltaR(*PackedCandidate_,jet),0,-0.6,0,-0.6);
            Npfcan_isGamma_[fillntupleentry] = 0;
            if(fabs(PackedCandidate_->pdgId())==22)  Npfcan_isGamma_[fillntupleentry] = 1;
            Npfcan_HadFrac_[fillntupleentry] = PackedCandidate_->hcalFraction();
			Npfcan_CaloFrac_[fillntupleentry] = PackedCandidate_->caloFraction();
            Npfcan_drminsv_[fillntupleentry] = catchInfsAndBound(drminpfcandsv_,0,-0.4,0,-0.4);

        }

    } // end loop over jet.numberOfDaughters()

    /*
    std::cout <<"numbers charged/neutrals"<<std::endl;
    std::cout << n_Cpfcand_ << std::endl;
    std::cout << n_Npfcand_ << std::endl;
    std::cout <<"charged IPs"<<std::endl;
    for(size_t i=0;i<n_Cpfcand_;i++){
        std::cout << Cpfcan_BtagPf_trackSip2dSig_[i] << " " << Npfcan_drminsv_[i] << " " << Npfcan_ptrel_[i]<<std::endl;
    }
    std::cout <<"neutrals minDR"<<std::endl;
    for(size_t i=0;i<n_Npfcand_;i++){
        std::cout << Npfcan_drminsv_[i] << " " << Npfcan_ptrel_[i]<<std::endl;
    }
     */

    nCpfcand_=n_Cpfcand_;
    nNpfcand_=n_Npfcand_;

    return true; //for making cuts
}


float ntuple_pfCands::mindrsvpfcand(const pat::PackedCandidate* pfcand) {

    float mindr_ = jetradius_;
    for (unsigned int i=0; i<secVertices()->size(); ++i) {
        if(!pfcand) continue;
        //if(!svs.at(i)) continue;
        float tempdr_ = reco::deltaR(secVertices()->at(i),*pfcand);
        if (tempdr_<mindr_) { mindr_ = tempdr_; }

    }
    return mindr_;
}

float ntuple_pfCands::mindistsvpfcand(const reco::TransientTrack track) {

  float mindist_ = 999.999;
  float out_dist = 0.0;
  for (unsigned int i=0; i<secVertices()->size(); ++i) {
    if (!track.isValid()) {continue;}
    reco::Vertex::CovarianceMatrix csv; secVertices()->at(i).fillVertexCovariance(csv);
    reco::Vertex vertex(secVertices()->at(i).vertex(), csv);
    if (!vertex.isValid()) {continue;}

    GlobalVector direction(secVertices()->at(i).px(),secVertices()->at(i).py(),secVertices()->at(i).pz());


    AnalyticalImpactPointExtrapolator extrapolator(track.field());
    TrajectoryStateOnSurface tsos =  extrapolator.extrapolate(track.impactPointState(), RecoVertex::convertPos(vertex.position()));
    
    VertexDistance3D dist;

    if (!tsos.isValid()) {continue;}
    GlobalPoint refPoint = tsos.globalPosition();
    GlobalError refPointErr = tsos.cartesianError().position();
    GlobalPoint vertexPosition = RecoVertex::convertPos(vertex.position());
    GlobalError vertexPositionErr = RecoVertex::convertError(vertex.error());

    std::pair<bool, Measurement1D> result(true, dist.distance(VertexState(vertexPosition, vertexPositionErr), VertexState(refPoint, refPointErr)));
    if (!result.first) {continue;}

    GlobalPoint impactPoint = tsos.globalPosition();
    GlobalVector IPVec(impactPoint.x() - vertex.x(), impactPoint.y() - vertex.y(), impactPoint.z() - vertex.z());
    double prod = IPVec.dot(direction);
    double sign = (prod >= 0) ? 1. : -1.;

    if(result.second.value() < mindist_){
      float inv_dist = 1.0 / (sign * result.second.value() + 0.001);
      out_dist = inv_dist;
      mindist_ = result.second.value();
    } 
  }
  return out_dist;
}
