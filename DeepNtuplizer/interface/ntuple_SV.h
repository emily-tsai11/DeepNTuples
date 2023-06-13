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
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

class ntuple_SV: public ntuple_content{
public:

    ntuple_SV(std::string prefix = "", double jetR = 0.4);
    ~ntuple_SV();

    void getInput(const edm::ParameterSet& iConfig);
    void initBranches(TTree* );
    void readEvent(const edm::Event& iEvent);
	void readSetup(const edm::EventSetup& iSetup);
    //use either of these functions

//$$    bool fillBranches(const pat::Jet &, const size_t& jetidx, const  edm::View<pat::Jet> * coll=0);
    bool fillBranches(const pat::Jet &, const size_t& jetidx, const  edm::View<pat::Jet> * coll=0, float EventTime = -1);
	
	void setTrackBuilderToken(const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord>& track_builder_token) {
		track_builder_token_ = track_builder_token;
    }

private:

    // SV candidates
    int   sv_num_;
    float nsv_;
    std::string prefix_;

    edm::ESHandle<TransientTrackBuilder> builder;
    edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> track_builder_token_;

    static constexpr size_t max_sv=10;

    float sv_pt_[max_sv];
    float sv_eta_[max_sv];
    float sv_phi_[max_sv];

    float sv_e_[max_sv];

    float sv_etarel_[max_sv];
    float sv_phirel_[max_sv];
    float sv_deltaR_[max_sv];
    float sv_mass_[max_sv];
    //  float sv_phirel_[max_sv];
    //  float sv_etarel_[max_sv];
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
//$$
    float sv_time_[max_sv];
//$$


    static const reco::Vertex * spvp_;

    static bool compareDxyDxyErr(const reco::VertexCompositePtrCandidate &sva,const reco::VertexCompositePtrCandidate &svb);

    //helper functions:
    static Measurement1D vertexDxy(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv)  ;
    static Measurement1D vertexD3d(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv)  ;
    static float vertexDdotP(const reco::VertexCompositePtrCandidate &sv, const reco::Vertex &pv)  ;

};



#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_SV_H_ */
