/*
 * ntuple_pfcands.h
 *
 *  Created on: 13 Feb 2017
 *      Author: jkiesele
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_PFCANDS_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_PFCANDS_H_

#include "ntuple_content.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

class ntuple_pfCands: public ntuple_content{
public:

    ntuple_pfCands():ntuple_content(),jetradius_(0.4),
    n_Cpfcand_(0),n_Npfcand_(0){}

    void getInput(const edm::ParameterSet& iConfig);
    void initBranches(TTree* );
    void readEvent(const edm::Event& iEvent);
    void readSetup(const edm::EventSetup& iSetup);
    void initContainers() {}
    void clearContainers() {}
    void deleteContainers() {}

    //use either of these functions

//$$    bool fillBranches(const pat::Jet &, const size_t& jetidx, const  edm::View<pat::Jet> * coll=0);
    bool fillBranches(const pat::Jet &, const size_t& jetidx, const  edm::View<pat::Jet> * coll=0, float EventTime = -1);
    bool fillBranches() { return false; }

    void setJetRadius(const float& radius){jetradius_=radius;}
    void setTrackBuilderToken(const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord>& track_builder_token) {
		track_builder_token_ = track_builder_token;
	}

private:

    float jetradius_;
		float min_candidate_pt_ = -1;

    edm::ESHandle<TransientTrackBuilder> builder;
	edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> track_builder_token_;

    unsigned int n_Cpfcand_;
    float nCpfcand_;

    static constexpr size_t max_pfcand_=50;

    float  Cpfcan_pt_[max_pfcand_];
    float  Cpfcan_eta_[max_pfcand_];
    float  Cpfcan_phi_[max_pfcand_];
    float  Cpfcan_ptrel_[max_pfcand_];
    float  Cpfcan_e_[max_pfcand_];
    float  Cpfcan_erel_[max_pfcand_];
    float  Cpfcan_phirel_[max_pfcand_];
    float  Cpfcan_etarel_[max_pfcand_];
    float  Cpfcan_deltaR_[max_pfcand_];
    float  Cpfcan_puppiw_[max_pfcand_];
    float   Cpfcan_VTX_ass_[max_pfcand_];

    float   Cpfcan_fromPV_[max_pfcand_];

    float Cpfcan_vertexChi2_[max_pfcand_];
    float Cpfcan_vertexNdof_[max_pfcand_];
    float Cpfcan_vertexNormalizedChi2_[max_pfcand_];
    float Cpfcan_vertex_rho_[max_pfcand_];
    float Cpfcan_vertex_phirel_[max_pfcand_];
    float Cpfcan_vertex_etarel_[max_pfcand_];
    float Cpfcan_vertexRef_mass_[max_pfcand_];

    // covariance
    float  Cpfcan_dz_[max_pfcand_];
    float  Cpfcan_dxy_[max_pfcand_];

    float  Cpfcan_dxyerrinv_[max_pfcand_];
    float  Cpfcan_dxysig_[max_pfcand_];

    float  Cpfcan_dptdpt_[max_pfcand_];
    float  Cpfcan_detadeta_[max_pfcand_];
    float  Cpfcan_dphidphi_[max_pfcand_];
    float  Cpfcan_dxydxy_[max_pfcand_];
    float  Cpfcan_dzdz_[max_pfcand_];
    float  Cpfcan_dxydz_[max_pfcand_];
    float  Cpfcan_dphidxy_[max_pfcand_];
    float  Cpfcan_dlambdadz_[max_pfcand_];

//$$
    float  Cpfcan_time_[max_pfcand_];
//$$

    float Cpfcan_BtagPf_trackMomentum_[max_pfcand_];
    float Cpfcan_BtagPf_trackEta_[max_pfcand_];
    float Cpfcan_BtagPf_trackEtaRel_[max_pfcand_];
    float Cpfcan_BtagPf_trackPtRel_[max_pfcand_];
    float Cpfcan_BtagPf_trackPPar_[max_pfcand_];
    float Cpfcan_BtagPf_trackDeltaR_[max_pfcand_];
    float Cpfcan_BtagPf_trackPtRatio_[max_pfcand_];
    float Cpfcan_BtagPf_trackPParRatio_[max_pfcand_];
    float Cpfcan_BtagPf_trackSip3dVal_[max_pfcand_];
    float Cpfcan_BtagPf_trackSip3dSig_[max_pfcand_];
    float Cpfcan_BtagPf_trackSip2dVal_[max_pfcand_];
    float Cpfcan_BtagPf_trackSip2dSig_[max_pfcand_];

    float Cpfcan_BtagPf_trackDecayLen_[max_pfcand_];

    float Cpfcan_BtagPf_trackJetDistVal_[max_pfcand_];
    float Cpfcan_BtagPf_trackJetDistSig_[max_pfcand_];

    // ID, skipped "charged hadron" as that is true if now the other
    // TODO (comment of Markus Stoye) add reco information
    float Cpfcan_isMu_[max_pfcand_]; // pitty that the quality is missing
    float Cpfcan_isEl_[max_pfcand_]; // pitty that the quality is missing
    float Cpfcan_charge_[max_pfcand_];

    // track quality
    float Cpfcan_lostInnerHits_[max_pfcand_];
    float Cpfcan_numberOfPixelHits_[max_pfcand_];

    float Cpfcan_chi2_[max_pfcand_];
    float Cpfcan_quality_[max_pfcand_];

    float Cpfcan_drminsv_[max_pfcand_];
	float Cpfcan_distminsv_[max_pfcand_];
    float Cpfcan_distminsv2_[max_pfcand_];
    
    float  Cpfcan_pdgID_[max_pfcand_];

    float  Cpfcan_HadFrac_[max_pfcand_];
    float  Cpfcan_CaloFrac_[max_pfcand_];

 //hit pattern variables, as defined here https://github.com/cms-sw/cmssw/blob/master/DataFormats/TrackReco/interface/HitPattern.h
    //Tracker per layer 
    //Pixel barrel 
    float Cpfcan_nhitpixelBarrelLayer1_[max_pfcand_]; 
    float Cpfcan_nhitpixelBarrelLayer2_[max_pfcand_];
    float Cpfcan_nhitpixelBarrelLayer3_[max_pfcand_];
    float Cpfcan_nhitpixelBarrelLayer4_[max_pfcand_];
    //Pixel Endcap 
    float Cpfcan_nhitpixelEndcapLayer1_[max_pfcand_];
    float Cpfcan_nhitpixelEndcapLayer2_[max_pfcand_];
    //Strip TIB
    float Cpfcan_nhitstripTIBLayer1_[max_pfcand_];
    float Cpfcan_nhitstripTIBLayer2_[max_pfcand_];
    float Cpfcan_nhitstripTIBLayer3_[max_pfcand_];
    float Cpfcan_nhitstripTIBLayer4_[max_pfcand_];
    //Strip TID
    float Cpfcan_nhitstripTIDLayer1_[max_pfcand_];
    float Cpfcan_nhitstripTIDLayer2_[max_pfcand_];
    float Cpfcan_nhitstripTIDLayer3_[max_pfcand_];
    //Strip TOB
    float Cpfcan_nhitstripTOBLayer1_[max_pfcand_];
    float Cpfcan_nhitstripTOBLayer2_[max_pfcand_];
    float Cpfcan_nhitstripTOBLayer3_[max_pfcand_];
    float Cpfcan_nhitstripTOBLayer4_[max_pfcand_];
    float Cpfcan_nhitstripTOBLayer5_[max_pfcand_];
    float Cpfcan_nhitstripTOBLayer6_[max_pfcand_];
    //Strip TEC
    float Cpfcan_nhitstripTECLayer1_[max_pfcand_];
    float Cpfcan_nhitstripTECLayer2_[max_pfcand_];
    float Cpfcan_nhitstripTECLayer3_[max_pfcand_];
    float Cpfcan_nhitstripTECLayer4_[max_pfcand_];
    float Cpfcan_nhitstripTECLayer5_[max_pfcand_];
    float Cpfcan_nhitstripTECLayer6_[max_pfcand_];
    float Cpfcan_nhitstripTECLayer7_[max_pfcand_];
    float Cpfcan_nhitstripTECLayer8_[max_pfcand_];
    float Cpfcan_nhitstripTECLayer9_[max_pfcand_];
    //Tracker all layers together  
    //Valid hits
    float Cpfcan_numberOfValidHits_[max_pfcand_]; 
    float Cpfcan_numberOfValidTrackerHits_[max_pfcand_]; 
    float Cpfcan_numberOfValidPixelHits_[max_pfcand_]; 
    float Cpfcan_numberOfValidPixelBarrelHits_[max_pfcand_]; 
    float Cpfcan_numberOfValidPixelEndcapHits_[max_pfcand_]; 
    float Cpfcan_numberOfValidStripHits_[max_pfcand_]; 
    float Cpfcan_numberOfValidStripTIBHits_[max_pfcand_]; 
    float Cpfcan_numberOfValidStripTIDHits_[max_pfcand_]; 
    float Cpfcan_numberOfValidStripTOBHits_[max_pfcand_]; 
    float Cpfcan_numberOfValidStripTECHits_[max_pfcand_]; 
    //LayersWithMeasurement
    float Cpfcan_trackerLayersWithMeasurementOld_[max_pfcand_]; 
    float Cpfcan_trackerLayersWithMeasurement_[max_pfcand_]; 
    float Cpfcan_pixelLayersWithMeasurementOld_[max_pfcand_]; 
    float Cpfcan_pixelLayersWithMeasurement_[max_pfcand_]; 
    float Cpfcan_stripLayersWithMeasurement_[max_pfcand_]; 
    float Cpfcan_pixelBarrelLayersWithMeasurement_[max_pfcand_]; 
    float Cpfcan_pixelEndcapLayersWithMeasurement_[max_pfcand_]; 
    float Cpfcan_stripTIBLayersWithMeasurement_[max_pfcand_]; 
    float Cpfcan_stripTIDLayersWithMeasurement_[max_pfcand_]; 
    float Cpfcan_stripTOBLayersWithMeasurement_[max_pfcand_]; 
    float Cpfcan_stripTECLayersWithMeasurement_[max_pfcand_]; 
    //Null
    float Cpfcan_trackerLayersNull_[max_pfcand_];
    float Cpfcan_pixelLayersNull_[max_pfcand_]; 
    float Cpfcan_stripLayersNull_[max_pfcand_]; 
    float Cpfcan_pixelBarrelLayersNull_[max_pfcand_];
    float Cpfcan_pixelEndcapLayersNull_[max_pfcand_];
    float Cpfcan_stripTIBLayersNull_[max_pfcand_];
    float Cpfcan_stripTIDLayersNull_[max_pfcand_];
    float Cpfcan_stripTOBLayersNull_[max_pfcand_];
    float Cpfcan_stripTECLayersNull_[max_pfcand_];

    //Neutral Pf candidates
    unsigned int n_Npfcand_;
    float nNpfcand_;
    float  Npfcan_pt_[max_pfcand_];
    float  Npfcan_eta_[max_pfcand_];
    float  Npfcan_phi_[max_pfcand_];
    float  Npfcan_ptrel_[max_pfcand_];
    float  Npfcan_e_[max_pfcand_];
    float  Npfcan_erel_[max_pfcand_];
    float  Npfcan_puppiw_[max_pfcand_];
    float  Npfcan_phirel_[max_pfcand_];
    float  Npfcan_etarel_[max_pfcand_];
    float  Npfcan_deltaR_[max_pfcand_];
    float  Npfcan_isGamma_[max_pfcand_];
    float  Npfcan_HadFrac_[max_pfcand_];
    float  Npfcan_CaloFrac_[max_pfcand_];
    float  Npfcan_drminsv_[max_pfcand_];
	float  Npfcan_pdgID_[max_pfcand_];

    float mindrsvpfcand(const pat::PackedCandidate* pfcand);
	float mindistsvpfcand(const reco::TransientTrack track);
};


#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_PFCANDS_H_ */
