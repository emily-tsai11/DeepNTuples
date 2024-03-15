#include <string>
#include <vector>
#include <algorithm>

// CMSSW includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
// #include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#if defined( __GXX_EXPERIMENTAL_CXX0X__)
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#endif

// For IVF
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"

// User includes
#include "../interface/ntuple_content.h"
#include "../interface/ntuple_SV.h"
#include "../interface/ntuple_JetInfo.h"
#include "../interface/ntuple_pfCands.h"
#include "../interface/ntuple_bTagVars.h"
#include "../interface/ntuple_FatJetInfo.h"
#include "../interface/ntuple_DeepVertex.h"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TROOT.h"
#include "TBranch.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TLorentzVector.h"


// struct MagneticField;
// const reco::TrackBaseRef toTrackRef(const reco::PFCandidate* pfcand);

class DeepNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {

    public:

        explicit DeepNtuplizer(const edm::ParameterSet&);
        ~DeepNtuplizer();
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:

        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        Measurement1D vertexDxy(const reco::VertexCompositePtrCandidate& svcand, const reco::Vertex& pv) const;
        Measurement1D vertexD3d(const reco::VertexCompositePtrCandidate& sv, const reco::Vertex& pv) const;
        float vertexDdotP(const reco::VertexCompositePtrCandidate& sv, const reco::Vertex& pv) const;

        // Member data
        edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
        edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> svToken_;
        // edm::EDGetTokenT<TrackingVertexCollection> gvToken_;
        edm::EDGetTokenT<edm::View<pat::Jet>> jetToken_;
        edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puToken_;
        edm::EDGetTokenT<edm::View<reco::BaseTagInfo>> pixHitsToken_;
        edm::EDGetTokenT<double> rhoToken_;
        edm::EDGetTokenT<float> genT0Token_;

        std::string t_qgtagger;

        edm::Service<TFileService> fs;
        TTree* tree_;

        float event_time_ = 0;
        const reco::Vertex* pv;

        size_t njetstotal_;
        size_t njetswithgenjet_;
        size_t njetsselected_;

        size_t njetstotal_crossCheck_;
        size_t njetswithgenjet_crossCheck_;
        size_t njetsselected_crossCheck_;

        ntuple_content* addModule(ntuple_content* m, std::string name = "") {
            modules_.push_back(m);
            module_names_.push_back(name);
            return m;
        }

        std::vector<ntuple_content*> modules_;
        std::vector<std::string> module_names_;

        bool applySelection_;
};


DeepNtuplizer::DeepNtuplizer(const edm::ParameterSet& iConfig) :
    vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    svToken_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secVertices"))),
    // gvToken_(consumes<TrackingVertexCollection>(iConfig.getParameter<edm::InputTag>("genVertices"))),
    jetToken_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"))),
    puToken_(consumes<std::vector<PileupSummaryInfo >>(iConfig.getParameter<edm::InputTag>("pupInfo"))),
    pixHitsToken_(consumes< edm::View<reco::BaseTagInfo> > (iConfig.getParameter<edm::InputTag>("pixelhit"))),
    rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoInfo"))),
    genT0Token_(consumes<float>(iConfig.getParameter<edm::InputTag>("genT0Tag"))),
    t_qgtagger(iConfig.getParameter<std::string>("qgtagger")) {

    // Initialise the modules here -- everything else does not need to be
    // changed if modules don't interact

    // Read configuration parameters
    const double jetR = iConfig.getParameter<double>("jetR");
    const bool runFatJets_ = iConfig.getParameter<bool>("runFatJet");

    // Not implemented yet
    const bool useHerwigCompatibleMatching = iConfig.getParameter<bool>("useHerwigCompatible");
    const bool isHerwig = iConfig.getParameter<bool>("isHerwig");

    ntuple_content::useoffsets = iConfig.getParameter<bool>("useOffsets");

    applySelection_ = iConfig.getParameter<bool>("applySelection");

    ntuple_SV* svmodule = new ntuple_SV("", jetR);
    svmodule->setTrackBuilderToken(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder")));

    svmodule->setPFCandToken(consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidates")));
    svmodule->setLostTracksToken(consumes<pat::PackedCandidateCollection>(edm::InputTag("lostTracks")));

    svmodule->setPFMCMatchToken(consumes<edm::Association<reco::GenParticleCollection>>(edm::InputTag("packedPFCandidateToGenAssociation")));
    svmodule->setLTMCMatchToken(consumes<edm::Association<reco::GenParticleCollection>>(edm::InputTag("lostTracksToGenAssociation")));

    addModule(svmodule, "SVNtuple");

    // ntuple_DeepVertex* dvmodule = new ntuple_DeepVertex(jetR);
    // dvmodule->setTrackBuilderToken(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder")));
    // addModule(dvmodule, "DVNtuple");
 
    // Loose IVF vertices
    // ntuple_SV* svmodule_LooseIVF = new ntuple_SV("LooseIVF_", jetR);
    // svmodule_LooseIVF->setSVToken(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("LooseSVs")));
    // Removed LooseIVF module
    // addModule(svmodule_LooseIVF);

    // DeepVertex info
    // ntuple_DeepVertex* deepvertexmodule = new ntuple_DeepVertex(jetR);
    // deepvertexmodule->setCandidatesToken(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("candidates")));
    // addModule(deepvertexmodule);

    ntuple_JetInfo* jetinfo = new ntuple_JetInfo();
    jetinfo->setQglToken(consumes<edm::ValueMap<float>>(edm::InputTag(t_qgtagger, "qgLikelihood")));
    jetinfo->setPtDToken(consumes<edm::ValueMap<float>>(edm::InputTag(t_qgtagger, "ptD")));
    jetinfo->setAxis2Token(consumes<edm::ValueMap<float>>(edm::InputTag(t_qgtagger, "axis2")));
    jetinfo->setMultToken(consumes<edm::ValueMap<int>>(edm::InputTag(t_qgtagger, "mult")));

    jetinfo->setUseHerwigCompatibleMatching(useHerwigCompatibleMatching);
    jetinfo->setIsHerwig(isHerwig);

    jetinfo->setGenJetMatchReclusterToken(consumes<edm::Association<reco::GenJetCollection>>(iConfig.getParameter<edm::InputTag>("genJetMatchRecluster")));
    jetinfo->setGenJetMatchWithNuToken(consumes<edm::Association<reco::GenJetCollection>>(iConfig.getParameter<edm::InputTag>("genJetMatchWithNu")));

    jetinfo->setGenParticlesToken(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("pruned")));

    jetinfo->setMuonsToken(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons")));
    jetinfo->setElectronsToken(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons")));

    jetinfo->setPUInfoToken(consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("slimmedAddPileupInfo")));

    addModule(jetinfo, "jetinfo");

    ntuple_pfCands* pfcands = new ntuple_pfCands();
    pfcands->setJetRadius(jetR);
	pfcands->setTrackBuilderToken(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder")));

    addModule(pfcands, "pfcands");

    // std::cout << "PFcands Check" << std::endl;

    addModule(new ntuple_bTagVars(), "bTagVars");

    if (runFatJets_) {
        auto *fatjetinfo = new ntuple_FatJetInfo(jetR);
        fatjetinfo->setGenParticleToken(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("pruned")));
        fatjetinfo->setFatJetToken(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets")));
        addModule(fatjetinfo, "fatJets");
    }

    // Initialize modules and parse input parameters (if any)
    for (auto& m : modules_)
        m->getInput(iConfig);
}


DeepNtuplizer::~DeepNtuplizer() {

    return;
    for (auto& m : modules_)
        delete m;
}


// Method called for each event
void DeepNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    // Global info
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);
    if (vertices->empty()) return; // Skip the event if no PV found
    pv = &(*vertices->begin());

    edm::Handle<std::vector<reco::VertexCompositePtrCandidate>> secvertices;
    iEvent.getByToken(svToken_, secvertices);

    // edm::Handle<TrackingVertexCollection> genvertices;
    // iEvent.getByToken(gvToken_, genvertices);

    edm::Handle<std::vector<PileupSummaryInfo>> pupInfo;
    iEvent.getByToken(puToken_, pupInfo);

    edm::Handle<double> rhoInfo;
    iEvent.getByToken(rhoToken_, rhoInfo);

    edm::Handle<edm::View<pat::Jet>> jets;
    iEvent.getByToken(jetToken_, jets);

    edm::Handle<edm::View<reco::BaseTagInfo>> pixHits;
    iEvent.getByToken(pixHitsToken_, pixHits);

    for (auto& m : modules_) {
        m->setPrimaryVertices(vertices.product());
        m->setSecVertices(secvertices.product());
        // m->setGenVertices(genvertices.product());
        m->setJets(jets.product());
        m->setPuInfo(pupInfo.product());
        m->setRhoInfo(rhoInfo.product());
        m->readSetup(iSetup);
        m->readEvent(iEvent);
        m->initContainers();
    }

    // reco::VertexCompositePtrCandidateCollection cpvtx = *(secvertices.product());
    // for (const reco::VertexCompositePtrCandidate &sv : cpvtx) {
    //     cout << "SV_PT: " << sv.pt() << endl;
    // }

    std::vector<size_t> indices(jets->size());
    for (size_t i = 0; i < jets->size(); i++)
        indices.at(i) = i;

    if (applySelection_)
        std::random_shuffle(indices.begin(), indices.end());
 
    edm::View<pat::Jet>::const_iterator jetIter;
    // Loop over the jets
    // for (edm::View<pat::Jet>::const_iterator jetIter = jets->begin(); jetIter != jets->end(); ++jetIter) {

    float event_time = 0;
    // float event_timeWeight = 0;
    // float event_timeNtk = 0;

    // for (size_t j = 0; j < indices.size(); j++) {
    //     size_t jetidx = indices.at(j);
    //     jetIter = jets->begin() + jetidx;
    //     const pat::Jet& jet = *jetIter;
    //     for (unsigned int i = 0; i < jet.numberOfDaughters(); i++) {
    //         const pat::PackedCandidate* PackedCandidate = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(i));
    //         if (!PackedCandidate) continue;
    //         if (PackedCandidate->charge() == 0) continue;
    //         auto track = PackedCandidate->bestTrack();
    //         if (!track) continue;
    //         float track_dxy = track->dxy(pv->position());
    //         float track_dz = track->dz(pv->position());
    //         float track_time = track->t0();
    //         float track_timeError = track->covt0t0();
    //         float track_pt = track->pt();
    //         float time_weight = track_pt * track_pt;
    //         if (track_timeError > 0.0 && abs(track_time) < 1 && abs(track_dxy) < 0.05 && abs(track_dz) < 0.10 ) {
    //             event_timeNtk += 1;
    //             event_timeWeight += time_weight;
    //             event_time += track_time * time_weight;

    //             // std::cout << "    Jet  pt " << jet.pt
    //             //           << "   Track pt dxy dz time " << track_pt
    //             //           << " " << track_dxy << " " << track_dz << " " << track_time
    //             //           << std::endl;
    //         }
    //     }
    // }
    // if (event_timeNtk > 0) event_time /= event_timeWeight;
    // else event_time = -1;

    // Event time from generated t0:
    edm::Handle<float> genT0Handle;
    iEvent.getByToken(genT0Token_, genT0Handle);
    event_time = *genT0Handle.product();

    // std::cout << std::endl;
    // std::cout << " in DeepNtuplizer: Event time " << event_time << std::endl;

    event_time_ = event_time;

    int njet = 0;
    for (size_t j = 0; j < indices.size(); j++) {
        njetstotal_++;
        size_t jetidx = indices.at(j);
        jetIter = jets->begin() + jetidx;
        const pat::Jet& jet = *jetIter;

        if (jet.genJet()) njetswithgenjet_++;

        bool writejet = true;
        size_t idx = 0;
        for (auto& m : modules_) {
            // std::cout << module_names_[idx] << std::endl;
            // if (module_names_[idx] == "SVNtuple") writejet = false;
            // if (!m->fillBranches(jet, jetidx, jets.product())) {
            if (!m->fillBranches(jet, jetidx, jets.product(), event_time)) {
                writejet = false;
                if (applySelection_) break;
            }
            idx++;
        }

        if ((writejet && applySelection_) || !applySelection_) {
            // tree_->Fill();
            njetsselected_++;
        }

        njet++;
    } // End of looping over the jets

    njetstotal_crossCheck_ += jets->size();
    for (size_t j = 0; j < jets->size(); j++) {
        if (jets->at(j).genJet()) njetswithgenjet_crossCheck_++;
    }
    size_t index = 0;
    for (auto& m : modules_) {
        if (module_names_[index] == "SVNtuple" || module_names_[index] == "jetinfo") {
            m->fillBranches();
        }
        index++;
    }
    tree_->Fill();

    for (auto& m : modules_) {
        m->deleteContainers();
    }
}


// Method called once each job just before starting event loop
void DeepNtuplizer::beginJob() {

    if (!fs)
        throw edm::Exception(edm::errors::Configuration, "TFile Service is not registered in cfg file");

    tree_ = (fs->make<TTree>("tree", "tree"));
    tree_->Branch("Event_time", &event_time_, "Event_time/F");
    for (auto& m : modules_)
        m->initBranches(tree_);

    njetstotal_ = 0;
    njetswithgenjet_ = 0;
    njetsselected_ = 0;

    njetstotal_crossCheck_ = 0;
    njetswithgenjet_crossCheck_ = 0;
    njetsselected_crossCheck_ = 0;
}


// Method called once each job just after ending the event loop
void DeepNtuplizer::endJob() {

    std::cout << "total number of processed jets:     " << njetstotal_ << std::endl;
    std::cout << "total number of jets with gen:      " << njetswithgenjet_ << std::endl;
    std::cout << "total number of selected jets:      " << njetsselected_ << std::endl;
    std::cout << "fraction of selected jets:          " << (float) njetsselected_ / (float) njetstotal_ << std::endl;
    std::cout << "fraction of selected jets with gen: " << (float) njetsselected_ / (float) njetswithgenjet_ << std::endl;

    std::cout << "and now checking with new code" << std::endl;
    std::cout << "total number of processed jets (check):     " << njetstotal_crossCheck_ << std::endl;
    std::cout << "total number of jets with gen (check):      " << njetswithgenjet_crossCheck_ << std::endl;
    std::cout << "total number of selected jets (check):      " << njetsselected_crossCheck_ << std::endl;
    std::cout << "fraction of selected jets (check):          " << (float) njetsselected_crossCheck_ / (float) njetstotal_crossCheck_ << std::endl;
    std::cout << "fraction of selected jets with gen (check): " << (float) njetsselected_crossCheck_ / (float) njetswithgenjet_crossCheck_ << std::endl;
}

// Method fills 'descriptions' with the allowed parameters for the module
void DeepNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    // The following says we do not know what parameters are allowed so do no
    // validation -- please change this to state exactly what you do use, even
    // if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}


// TODO: Comment
// const reco::TrackBaseRef toTrackRef(const reco::PFCandidate* pfcand) {

//     if ((std::abs(pfcand->pdgId()) == 11 || pfcand->pdgId() == 22) && pfcand->gsfTrackRef().isNonnull() && pfcand->gsfTrackRef().isAvailable())
//         return reco::TrackBaseRef(pfcand->gsfTrackRef());
//     else if (pfcand->trackRef().isNonnull() && pfcand->trackRef().isAvailable())
//         return reco::TrackBaseRef(pfcand->trackRef());
//     else
//         return reco::TrackBaseRef();
// }


// Define this as a plug-in
DEFINE_FWK_MODULE(DeepNtuplizer);
