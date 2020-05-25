// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "../interface/ntuple_content.h"
#include "../interface/ntuple_SV.h"
#include "../interface/ntuple_JetInfo.h"
#include "../interface/ntuple_pfCands.h"
#include "../interface/ntuple_bTagVars.h"
#include "../interface/ntuple_FatJetInfo.h"
#include "../interface/ntuple_DeepVertex.h"

//ROOT includes
#include "TTree.h"
#include <TFile.h>
#include <TROOT.h>
#include "TBranch.h"
#include <string>
#include <vector>
#include "TSystem.h"
#include <TRandom.h>

//CMSSW includes
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

// for ivf
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TLorentzVector.h"

#include <algorithm>

#include "DataFormats/Math/interface/deltaR.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"

#if defined( __GXX_EXPERIMENTAL_CXX0X__)
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#endif

struct MagneticField;
const reco::TrackBaseRef toTrackRef(const reco::PFCandidate * pfcand);

class DeepNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
    explicit DeepNtuplizer(const edm::ParameterSet&);
    ~DeepNtuplizer();
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    Measurement1D vertexDxy(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv) const;
    Measurement1D vertexD3d(const reco::VertexCompositePtrCandidate &sv, const reco::Vertex &pv) const ;
    float vertexDdotP(const reco::VertexCompositePtrCandidate &sv, const reco::Vertex &pv) const ;

    // ----------member data ---------------------------
    edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
    edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> svToken_;
    edm::EDGetTokenT<edm::View<pat::Jet> >      jetToken_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puToken_;
    edm::EDGetTokenT<double> rhoToken_;

    std::string t_qgtagger;

    edm::Service<TFileService> fs;
    TTree *tree_;

//$$
    float event_time_ = 0;
    float jet_time_   = 0;

    const reco::Vertex  *pv;
    float jet_vertex_time_   = 0;
//$$

    size_t njetstotal_;
    size_t njetswithgenjet_;
    size_t njetsselected_;

	ntuple_content * addModule(ntuple_content *m, std::string name = ""){
        modules_.push_back(m);
	module_names_.push_back(name);
        return m;
    }
    std::vector<ntuple_content* > modules_;
    std::vector<std::string> module_names_;

    bool applySelection_;
};

DeepNtuplizer::DeepNtuplizer(const edm::ParameterSet& iConfig):
                            vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
                            svToken_(consumes<reco::VertexCompositePtrCandidateCollection>(
                                    iConfig.getParameter<edm::InputTag>("secVertices"))),
                                    jetToken_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"))),
                                    puToken_(consumes<std::vector<PileupSummaryInfo >>(iConfig.getParameter<edm::InputTag>("pupInfo"))),
                                    rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoInfo"))),
                                    t_qgtagger(iConfig.getParameter<std::string>("qgtagger"))
{
    /*
     *  Initialise the modules here
     *  Everything else does not need to be changed if
     *  modules don't interact.
     */

    // read configuration parameters
    const double jetR = iConfig.getParameter<double>("jetR");
    const bool  runFatJets_ = iConfig.getParameter<bool>("runFatJet");

    //not implemented yet
    const bool useHerwigCompatibleMatching=iConfig.getParameter<bool>("useHerwigCompatible");
    const bool isHerwig=iConfig.getParameter<bool>("isHerwig");

    ntuple_content::useoffsets = iConfig.getParameter<bool>("useOffsets");

    applySelection_=iConfig.getParameter<bool>("applySelection");

    ntuple_SV* svmodule=new ntuple_SV("", jetR);
    addModule(svmodule, "SVNtuple");

    //Loose IVF vertices
    //ntuple_SV* svmodule_LooseIVF=new ntuple_SV("LooseIVF_", jetR);
    //svmodule_LooseIVF->setSVToken(
    //        consumes<reco::VertexCompositePtrCandidateCollection>(
    //                iConfig.getParameter<edm::InputTag>("LooseSVs")));
    //removed LooseIVF module
    //addModule(svmodule_LooseIVF);

    // DeepVertex info
    // ntuple_DeepVertex* deepvertexmodule=new ntuple_DeepVertex(jetR);
    // deepvertexmodule->setCandidatesToken(consumes<edm::View<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("candidates")));
    // addModule(deepvertexmodule);

    ntuple_JetInfo* jetinfo=new ntuple_JetInfo();
    jetinfo->setQglToken(consumes<edm::ValueMap<float>>(edm::InputTag(t_qgtagger, "qgLikelihood")));
    jetinfo->setPtDToken(consumes<edm::ValueMap<float>>(edm::InputTag(t_qgtagger, "ptD")));
    jetinfo->setAxis2Token(consumes<edm::ValueMap<float>>(edm::InputTag(t_qgtagger, "axis2")));
    jetinfo->setMultToken(consumes<edm::ValueMap<int>>(edm::InputTag(t_qgtagger, "mult")));

    jetinfo->setUseHerwigCompatibleMatching(useHerwigCompatibleMatching);
    jetinfo->setIsHerwig(isHerwig);

    jetinfo->setGenJetMatchReclusterToken(
            consumes<edm::Association<reco::GenJetCollection> >(
                    iConfig.getParameter<edm::InputTag>( "genJetMatchRecluster" )));
    jetinfo->setGenJetMatchWithNuToken(
            consumes<edm::Association<reco::GenJetCollection> >(
                    iConfig.getParameter<edm::InputTag>( "genJetMatchWithNu" )));

    jetinfo->setGenParticlesToken(
            consumes<reco::GenParticleCollection>(
                    iConfig.getParameter<edm::InputTag>("pruned")));

    jetinfo->setMuonsToken(
            consumes<pat::MuonCollection>(
                    iConfig.getParameter<edm::InputTag>("muons")));

    jetinfo->setElectronsToken(
            consumes<pat::ElectronCollection>(
                    iConfig.getParameter<edm::InputTag>("electrons")));

    jetinfo->setPUInfoToken(
        consumes<std::vector<PileupSummaryInfo>>(
                    edm::InputTag("slimmedAddPileupInfo")));

    addModule(jetinfo, "jetinfo");

    ntuple_pfCands * pfcands = new ntuple_pfCands();
    pfcands->setJetRadius(jetR);

    addModule(pfcands, "pfcands");

    addModule(new ntuple_bTagVars(), "bTagVars");

    if(runFatJets_){
        auto *fatjetinfo = new ntuple_FatJetInfo(jetR);
        fatjetinfo->setGenParticleToken(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("pruned")));
        fatjetinfo->setFatJetToken(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets")));
        addModule(fatjetinfo, "fatJets");
    }
    /*
     *
     * Modules initialized
     *
     * parse the input parameters (if any)
     */
    for(auto& m: modules_)
        m->getInput(iConfig);

}


DeepNtuplizer::~DeepNtuplizer()
{
    return;
    for(auto& m:modules_)
        delete m;
}


// ------------ method called for each event  ------------
void
DeepNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    //global info

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);
    if (vertices->empty()) return; // skip the event if no PV found
//$$
    pv = &(*vertices->begin());
//$$

    edm::Handle<std::vector<reco::VertexCompositePtrCandidate> > secvertices;
    iEvent.getByToken(svToken_, secvertices);

    edm::Handle<std::vector<PileupSummaryInfo> > pupInfo;
    iEvent.getByToken(puToken_, pupInfo);

    edm::Handle<double> rhoInfo;
    iEvent.getByToken(rhoToken_,rhoInfo);

    edm::Handle<edm::View<pat::Jet> > jets;
    iEvent.getByToken(jetToken_, jets);

    for(auto& m:modules_){
        m->setPrimaryVertices(vertices.product());
        m->setSecVertices(secvertices.product());
        m->setPuInfo(pupInfo.product());
        m->setRhoInfo(rhoInfo.product());
        m->readSetup(iSetup);
        m->readEvent(iEvent);
    }

    std::vector<size_t> indices(jets->size());
    for(size_t i=0;i<jets->size();i++)
        indices.at(i)=i;

    if(applySelection_)
        std::random_shuffle (indices.begin(),indices.end());
 
    edm::View<pat::Jet>::const_iterator jetIter;
    // loop over the jets
    //for (edm::View<pat::Jet>::const_iterator jetIter = jets->begin(); jetIter != jets->end(); ++jetIter) {

//$$
    float event_time       = 0;
    float event_timeWeight = 0;
    float event_timeNtk    = 0;

    for (size_t j=0; j<indices.size(); j++) {
      size_t jetidx=indices.at(j);
      jetIter = jets->begin()+jetidx;
      const pat::Jet& jet = *jetIter;
      for (unsigned int i = 0; i <  jet.numberOfDaughters(); i++) {
        const pat::PackedCandidate* PackedCandidate = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(i));
      if ( !PackedCandidate ) continue;
      if ( PackedCandidate->charge() == 0 ) continue;
     	auto track = PackedCandidate->bestTrack();
      if ( !track ) continue;
        float track_dxy = track->dxy(pv->position());
        float track_dz  = track->dz(pv->position());
        float track_time    = track->t0();
        float track_timeError = track->covt0t0();
        float track_pt    = track->pt();
        float time_weight = track_pt * track_pt;
        if ( track_timeError > 0. && abs(track_time) < 1 
	     && abs(track_dxy) < 0.05 && abs(track_dz) < 0.10 ) {
          event_timeNtk    += 1;
          event_timeWeight += time_weight;
          event_time	   += track_time * time_weight;
        }
      }
    }
    if ( event_timeNtk > 0 ) event_time /= event_timeWeight;
    else                     event_time = -1;

//     std::cout << std::endl;
//     std::cout << " event PVz time " << (*vertices)[0].z() << " " << event_time << std::endl;

    event_time_ = event_time;

    int njet = 0;
    for (size_t j=0; j<indices.size(); j++) {
      njetstotal_++;
      size_t jetidx=indices.at(j);
      jetIter = jets->begin()+jetidx;
      const pat::Jet& jet = *jetIter;

      if ( jet.genJet() ) njetswithgenjet_++;

      bool writejet=true;
      size_t idx = 0;
      for ( auto& m:modules_ ) {
        //std::cout << module_names_[idx] << std::endl;
//$$        if(! m->fillBranches(jet, jetidx, jets.product())){
        if(! m->fillBranches(jet, jetidx, jets.product(), event_time)){
            writejet=false;
            if(applySelection_) break;
        }
        idx++;
      }

//$$
      float jet_time	   = 0;
      float jet_timeWeight = 0;
      float jet_timeNtk    = 0;

      float jet_vertex_time	  = 0;
      float jet_vertex_timeWeight = 0;
      float jet_vertex_timeNtk    = 0;

      int  nSV = -2;
      const reco::CandSecondaryVertexTagInfo *candSVTagInfo = jet.tagInfoCandSecondaryVertex("pfInclusiveSecondaryVertexFinder");
      if ( candSVTagInfo != nullptr ) nSV = candSVTagInfo->nVertices();
      if ( nSV > 0 && candSVTagInfo->vertexTracks().size() == 0 ) nSV = -1;

//   std::cout << "        jet " << j << " pt eta phi " 
// 	 << jet.pt() << " " << jet.eta() << " " << jet.phi() << "   nSV " << nSV << std::endl;

      for (unsigned int i = 0; i <  jet.numberOfDaughters(); i++) {
        const pat::PackedCandidate* PackedCandidate = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(i));
      if ( !PackedCandidate ) continue;
      if ( PackedCandidate->charge() == 0 ) continue;
        auto track = PackedCandidate->bestTrack();
      if ( !track ) continue;
        float track_time      = track->t0();
        float track_timeError = track->covt0t0();
        float track_pt    = track->pt();
        float time_weight = track_pt * track_pt;
        if ( track_timeError > 0. && abs(track_time) < 1 ) {
          jet_timeNtk += 1;
          jet_timeWeight += time_weight;
          jet_time += track_time * time_weight;
        }

  	bool SVtrack = false;
	if ( nSV > 0  && track_timeError > 0. && abs(track_time) < 1 ) {
	  for (unsigned int isv=0; isv<candSVTagInfo->nVertices(); ++isv) {
	    for (unsigned int it=0; it<candSVTagInfo->nVertexTracks(isv); ++it) {
	    if ( candSVTagInfo->vertexTracks(isv)[it]->charge() != track->charge() ) continue;
	      float dpt  = TMath::Abs(candSVTagInfo->vertexTracks(isv)[it]->pt()  / track->pt() - 1.);
              float deta = TMath::Abs(candSVTagInfo->vertexTracks(isv)[it]->eta() - track->eta());
              float dphi = TMath::Abs(candSVTagInfo->vertexTracks(isv)[it]->phi() - track->phi());
	      if (dphi > 3.141593 ) dphi -= 2.*3.141593;
	      if (dpt < 0.01 && deta < 0.01 && dphi < 0.01) SVtrack = true;
	    }
	  }
	}
	if ( SVtrack ) {
          jet_vertex_timeNtk     += 1;
          jet_vertex_timeWeight  += time_weight;
          jet_vertex_time        += track_time * time_weight;
	}
      }  // end loop on tracks in jets

      if ( nSV > 0 ) {
        if ( jet_vertex_timeNtk > 0 && event_timeNtk > 0 ) {
          jet_vertex_time = jet_vertex_time / jet_vertex_timeWeight - event_time;
	  jet_vertex_time = TMath::Abs(jet_vertex_time);
	}
        else jet_vertex_time = -1; 
      }
      else jet_vertex_time = -1;

      if ( jet_timeNtk > 0 && event_timeNtk > 0 ) {
        jet_time = jet_time / jet_timeWeight - event_time;
	jet_time = TMath::Abs(jet_time);
      }
      else jet_time = -1; 

      jet_time_ = jet_time;
      jet_vertex_time_ = jet_vertex_time;
//$$

      if( (writejet&&applySelection_) || !applySelection_ ){
          tree_->Fill();
          njetsselected_++;
      }

      njet++;
    } // end of looping over the jets

}


// ------------ method called once each job just before starting event loop  ------------
    void
DeepNtuplizer::beginJob()
{
    if( !fs ){
        throw edm::Exception( edm::errors::Configuration,
                "TFile Service is not registered in cfg file" );
    }
    tree_=(fs->make<TTree>("tree" ,"tree" ));
//$$
    tree_->Branch("Event_time", &event_time_, "Event_time/F");
    tree_->Branch("Jet_time", &jet_time_, "Jet_time/F");
    tree_->Branch("Jet_vertex_time", &jet_vertex_time_, "Jet_vertex_time/F");
//$$


    for(auto& m:modules_)
        m->initBranches(tree_);

    njetstotal_=0;
    njetswithgenjet_=0;
    njetsselected_=0;
}

// ------------ method called once each job just after ending the event loop  ------------
    void
DeepNtuplizer::endJob()
{

    std::cout << "total number of processed jets: " << njetstotal_<<std::endl;
    std::cout << "total number of jets with gen:  " << njetswithgenjet_<<std::endl;
    std::cout << "total number of selected jets:  "<< njetsselected_<<std::endl;
    std::cout << "fraction of selected jets:      "<< (float)njetsselected_/(float)njetstotal_<<std::endl;
    std::cout << "fraction of selected jets with gen: "<< (float)njetsselected_/(float)njetswithgenjet_<<std::endl;

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DeepNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

const reco::TrackBaseRef toTrackRef(const reco::PFCandidate * pfcand){
    if ( (std::abs(pfcand->pdgId()) == 11 || pfcand->pdgId() == 22) && pfcand->gsfTrackRef().isNonnull() && pfcand->gsfTrackRef().isAvailable() )
      return reco::TrackBaseRef(pfcand->gsfTrackRef());
    else if ( pfcand->trackRef().isNonnull() && pfcand->trackRef().isAvailable() )
          return reco::TrackBaseRef(pfcand->trackRef());
    else
        return reco::TrackBaseRef();
}


//define this as a plug-in
DEFINE_FWK_MODULE(DeepNtuplizer);
