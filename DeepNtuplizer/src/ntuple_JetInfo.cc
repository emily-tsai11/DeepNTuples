/*
 * ntuple_JetInfo.cc
 *
 *  Created on: 13 Feb 2017
 *      Author: jkiesele
 */


#include "../interface/ntuple_JetInfo.h"
#include "../interface/helpers.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <TF1.h>
#include <vector>
#include <algorithm>
#include <memory>

// User include files
// #include "FWCore/Framework/interface/Frameworkfwd.h"
// #include "FWCore/Framework/interface/EDAnalyzer.h"
// #include "FWCore/Framework/interface/Event.h"
// #include "FWCore/Framework/interface/MakerMacros.h"
// #include "FWCore/ParameterSet/interface/ParameterSet.h"
// #include "FWCore/ServiceRegistry/interface/Service.h"
// #include "CommonTools/UtilAlgos/interface/TFileService.h"
// #include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
// #include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"


using namespace std;


ntuple_JetInfo::ntuple_JetInfo() : ntuple_content(), gluonReduction_(0), useherwcompat_matching_(false), isherwig_(false) {}


ntuple_JetInfo::~ntuple_JetInfo() {}


void ntuple_JetInfo::getInput(const edm::ParameterSet& iConfig) {

    gluonReduction_ = (iConfig.getParameter<double>("gluonReduction"));
    jetPtMin_ = (iConfig.getParameter<double>("jetPtMin"));
    jetPtMax_ = (iConfig.getParameter<double>("jetPtMax"));
    jetAbsEtaMin_ = (iConfig.getParameter<double>("jetAbsEtaMin"));
    jetAbsEtaMax_ = (iConfig.getParameter<double>("jetAbsEtaMax"));

    vector<string> disc_names = iConfig.getParameter<vector<string>>("bDiscriminators");
    for (auto& name : disc_names) {
        discriminators_[name] = 0;
    }
}


void ntuple_JetInfo::initBranches(TTree* tree) {

    // More general event info, here applied per jet
    addBranch(tree, "nJets", &nJets_, "nJets/I");
    addBranch(tree, "npv", &npv_, "npv/F");
    addBranch(tree, "npv_0_z", &npv_0_z_, "npv_0_z/F");
    addBranch(tree, "PU_rho", &PU_rho_, "PU_rho/F");
    addBranch(tree, "ntrueInt", &ntrueInt_, "ntrueInt/F");
    addBranch(tree, "rho", &rho_, "rho/F");
    addBranch(tree, "event_no", &event_no_, "event_no/I");
    // addBranch(tree, "jet_no", &jet_no_, "jet_no/I");

    // Truth labels
    addBranch(tree, "isB", &isB_);
    addBranch(tree, "isGBB", &isGBB_);
    addBranch(tree, "isBB", &isBB_);
    addBranch(tree, "isLeptonicB", &isLeptonicB_);
    addBranch(tree, "isLeptonicB_C", &isLeptonicB_C_);
    addBranch(tree, "isC", &isC_);
    addBranch(tree, "isGCC", &isGCC_);
    addBranch(tree, "isCC", &isCC_);
    // addBranch(tree, "isTau", &isTau_);
    addBranch(tree, "isUD", &isUD_);
    addBranch(tree, "isS", &isS_);
    addBranch(tree, "isG", &isG_);
    addBranch(tree, "isUndefined", &isUndefined_);
    addBranch(tree, "genDecay", &genDecay_); // dxy corresponds to the distance the Bhadron traveled

    // Truth labeling with fallback to physics definition for light/gluon/undefined of standard flavor definition
    addBranch(tree, "isPhysB", &isPhysB_);
    addBranch(tree, "isPhysGBB", &isPhysGBB_);
    addBranch(tree, "isPhysBB", &isPhysBB_);
    addBranch(tree, "isPhysLeptonicB", &isPhysLeptonicB_);
    addBranch(tree, "isPhysLeptonicB_C", &isPhysLeptonicB_C_);
    addBranch(tree, "isPhysC", &isPhysC_);
    addBranch(tree, "isPhysGCC", &isPhysGCC_);
    addBranch(tree, "isPhysCC", &isPhysCC_);
    // addBranch(tree, "isPhysTau", &isPhysTau_);
    addBranch(tree, "isPhysUD", &isPhysUD_);
    addBranch(tree, "isPhysS", &isPhysS_);
    addBranch(tree, "isPhysG", &isPhysG_);
    addBranch(tree, "isPhysUndefined", &isPhysUndefined_);

    // Jet variables
    // b = tree->Branch("jet_pt", &jet_pt_);
    addBranch(tree, "jet_pt", &jet_pt_);
    addBranch(tree, "jet_corr_pt", &jet_corr_pt_);
    addBranch(tree, "jet_eta", &jet_eta_);
    addBranch(tree, "jet_phi", &jet_phi_);
    addBranch(tree, "jet_mass", &jet_mass_);
    addBranch(tree, "jet_energy", &jet_energy_);

    // Jet id
    addBranch(tree, "jet_looseId", &jet_looseId_);

    // Quark gluon
    addBranch(tree, "jet_qgl", &jet_qgl_); // qg tagger from jmar
    addBranch(tree, "QG_ptD", &QG_ptD_); // momentum fraction per jet constituent
    addBranch(tree, "QG_axis2", &QG_axis2_); // jet shape i.e. gluon are wider than quarks
    addBranch(tree, "QG_mult", &QG_mult_); // multiplicity i.e. total num of PFcands reconstructed

    // Yutas quark-gluon info
    addBranch(tree, "y_multiplicity", &y_multiplicity_);
    addBranch(tree, "y_charged_multiplicity", &y_charged_multiplicity_);
    addBranch(tree, "y_neutral_multiplicity", &y_neutral_multiplicity_);
    addBranch(tree, "y_ptD", &y_ptD_);
    addBranch(tree, "y_axis1", &y_axis1_);
    addBranch(tree, "y_axis2", &y_axis2_);
    addBranch(tree, "y_pt_dr_log", &y_pt_dr_log_);

    // In the jet
    addBranch(tree, "muons_number", &muons_number_);
    addBranch(tree, "muons_isLooseMuon", &muons_isLooseMuon_);
    addBranch(tree, "muons_isTightMuon", &muons_isTightMuon_);
    addBranch(tree, "muons_isSoftMuon", &muons_isSoftMuon_);
    addBranch(tree, "muons_isHighPtMuon", &muons_isHighPtMuon_);
    addBranch(tree, "muons_pt", &muons_pt_);
    addBranch(tree, "muons_relEta", &muons_relEta_);
    addBranch(tree, "muons_relPhi", &muons_relPhi_);
    addBranch(tree, "muons_energy", &muons_energy_);
    addBranch(tree, "electrons_number", &electrons_number_);
    addBranch(tree, "electrons_pt", &electrons_pt_);
    addBranch(tree, "electrons_relEta", &electrons_relEta_);
    addBranch(tree, "electrons_relPhi", &electrons_relPhi_);
    addBranch(tree, "electrons_energy", &electrons_energy_);

    // More truth labels
    addBranch(tree, "gen_pt", &gen_pt_);
    addBranch(tree, "gen_pt_Recluster", &gen_pt_Recluster_);
    addBranch(tree, "gen_pt_WithNu", &gen_pt_WithNu_);
    addBranch(tree, "Delta_gen_pt", &Delta_gen_pt_);
    addBranch(tree, "Delta_gen_pt_Recluster", &Delta_gen_pt_Recluster_);
    addBranch(tree, "Delta_gen_pt_WithNu", &Delta_gen_pt_WithNu_);

    if (1) { // discriminators might need to be filled differently -- FIXME
        for (auto& entry : discriminators_) {
            string better_name(entry.first);
            std::replace(better_name.begin(), better_name.end(), ':', '_');
            addBranch(tree, better_name.c_str(), &entry.second);
        }
    }
}


void ntuple_JetInfo::readEvent(const edm::Event& iEvent) {

    iEvent.getByToken(qglToken_, qglHandle_);
    iEvent.getByToken(ptDToken_, ptDHandle_);
    iEvent.getByToken(axis2Token_, axis2Handle_);
    iEvent.getByToken(multToken_, multHandle_);

    iEvent.getByToken(genJetMatchReclusterToken_, genJetMatchRecluster_);
    iEvent.getByToken(genJetMatchWithNuToken_, genJetMatchWithNu_);

    iEvent.getByToken(genParticlesToken_, genParticlesHandle_);

    iEvent.getByToken(muonsToken_, muonsHandle_);
    iEvent.getByToken(electronsToken_, electronsHandle_);

    iEvent.getByToken(puInfoToken_, PUInfo_);

    event_no_ = iEvent.id().event();

    // Presumably this whole part can be removed!
    neutrinosLepB_.clear();
    neutrinosLepB_C_.clear();
    gToBB_.clear();
    gToCC_.clear();
    alltaus_.clear();
    Bhadron_.clear();
    Bhadron_daughter_.clear();

    // std::cout << "start search for a b in this event " << std::endl;
    for (const reco::Candidate& genC : *genParticlesHandle_) {
        const reco::GenParticle& gen = static_cast<const reco::GenParticle&>(genC);

        if ((abs(gen.pdgId()) > 500 && abs(gen.pdgId()) < 600) || (abs(gen.pdgId()) > 5000 && abs(gen.pdgId()) < 6000)) {
            // std::cout << gen.end_vertex() << endl;

            Bhadron_.push_back(gen);
            if (gen.numberOfDaughters() > 0) {
                if ((abs(gen.daughter(0)->pdgId()) > 500 && abs(gen.daughter(0)->pdgId()) < 600)
                        || (abs(gen.daughter(0)->pdgId()) > 5000 && abs(gen.daughter(0)->pdgId()) < 6000)) {

                    if (gen.daughter(0)->numberOfDaughters() > 0) {
                        const reco::GenParticle& daughter_ = static_cast<const reco::GenParticle&>(*(gen.daughter(0)->daughter(0)));
                        if (daughter_.vx() != gen.vx()) { 
                            Bhadron_daughter_.push_back(daughter_);
                        }
                        else Bhadron_daughter_.push_back(gen);
                        //     std::cout << "only b daughters " << endl;
                        // }
                    }
                    else Bhadron_daughter_.push_back(gen);
                }
                else {
                    // std::cout << gen.daughter(0)->vx() << " oh " << gen.vx() << " " << gen.pt() << " " << gen.daughter(0)->pdgId() << std::endl; 
                    const reco::GenParticle& daughter_ = static_cast<const reco::GenParticle&>(*gen.daughter(0));
                    Bhadron_daughter_.push_back(daughter_);
                }
            } // End if daughter is there
            else {
                // std::cout << "lonely B hadron, has NO daughter???" << std::endl;
                Bhadron_daughter_.push_back(gen);
            }
        }
    }

    for (const reco::Candidate& genC : *genParticlesHandle_) {
        const reco::GenParticle& gen = static_cast<const reco::GenParticle&>(genC);
        if (abs(gen.pdgId()) == 12 || abs(gen.pdgId()) == 14 || abs(gen.pdgId()) == 16) {
            const reco::GenParticle* mother = static_cast<const reco::GenParticle*>(gen.mother());
            if (mother != NULL) {
                if ((abs(mother->pdgId()) > 500 && abs(mother->pdgId()) < 600) || (abs(mother->pdgId()) > 5000 && abs(mother->pdgId()) < 6000)) {
                    neutrinosLepB_.emplace_back(gen);
                }
                if ((abs(mother->pdgId()) > 400 && abs(mother->pdgId()) < 500) || (abs(mother->pdgId()) > 4000 && abs(mother->pdgId()) < 5000)) {
                    neutrinosLepB_C_.emplace_back(gen);
                }
            }
            else {
                std::cout << "No mother" << std::endl;
            }
        }

        int id(std::abs(gen.pdgId())); 
        int status(gen.status());

        if (id == 21 && status >= 21 && status <= 59) { // Pythia8 hard scatter, ISR, or FSR
            if (gen.numberOfDaughters() == 2) {
                const reco::Candidate* d0 = gen.daughter(0);
                const reco::Candidate* d1 = gen.daughter(1);
                if (std::abs(d0->pdgId()) == 5 && std::abs(d1->pdgId()) == 5
                        && d0->pdgId() * d1->pdgId() < 0 && reco::deltaR(*d0, *d1) < 0.4) gToBB_.push_back(gen);
                if (std::abs(d0->pdgId()) == 4 && std::abs(d1->pdgId()) == 4
                        && d0->pdgId() * d1->pdgId() < 0 && reco::deltaR(*d0, *d1) < 0.4) gToCC_.push_back(gen);
            }
        }

        if (id == 15 && false) {
            alltaus_.push_back(gen);
        }
    }
    // Technically a branch fill but per event, therefore here
}


void ntuple_JetInfo::initContainers() {

    isB_ = new std::vector<int>;
    isGBB_ = new std::vector<int>;
    isBB_ = new std::vector<int>;
    isC_ = new std::vector<int>;
    isGCC_ = new std::vector<int>;
    isCC_ = new std::vector<int>;
    isUD_ = new std::vector<int>;
    isS_ = new std::vector<int>;
    isG_ = new std::vector<int>;
    isUndefined_ = new std::vector<int>;
    isLeptonicB_ = new std::vector<int>;
    isLeptonicB_C_ = new std::vector<int>;
    isTau_ = new std::vector<int>;

    genDecay_ = new std::vector<float>;

    isPhysB_ = new std::vector<int>;
    isPhysGBB_ = new std::vector<int>;
    isPhysBB_ = new std::vector<int>;
    isPhysC_ = new std::vector<int>;
    isPhysGCC_ = new std::vector<int>;
    isPhysCC_ = new std::vector<int>;
    isPhysUD_ = new std::vector<int>;
    isPhysS_ = new std::vector<int>;
    isPhysG_ = new std::vector<int>;
    isPhysUndefined_ = new std::vector<int>;
    isPhysLeptonicB_ = new std::vector<int>;
    isPhysLeptonicB_C_ = new std::vector<int>;
    isPhysTau_ = new std::vector<int>;

    jet_qgl_ = new std::vector<float>;
    QG_ptD_ = new std::vector<float>;
    QG_axis2_ = new std::vector<float>;
    QG_mult_ = new std::vector<float>;

    muons_number_ = new std::vector<int>;
    electrons_number_ = new std::vector<int>;

    muons_jetIdx_ = new std::vector<int>;
    muons_isLooseMuon_ = new std::vector<int>;
    muons_isTightMuon_ = new std::vector<int>;
    muons_isSoftMuon_ = new std::vector<int>;
    muons_isHighPtMuon_ = new std::vector<int>;
    muons_pt_ = new std::vector<float>;
    muons_relEta_ = new std::vector<float>;
    muons_relPhi_ = new std::vector<float>;
    muons_energy_ = new std::vector<float>;

    electrons_jetIdx_ = new std::vector<int>;
    electrons_pt_ = new std::vector<float>;
    electrons_relEta_ = new std::vector<float>;
    electrons_relPhi_ = new std::vector<float>;
    electrons_energy_ = new std::vector<float>;

    jet_pt_ = new std::vector<float>;
    jet_corr_pt_ = new std::vector<float>;
    jet_eta_ = new std::vector<float>;
    jet_phi_ = new std::vector<float>;
    jet_mass_ = new std::vector<float>;
    jet_energy_ = new std::vector<float>;
    jet_looseId_ = new std::vector<float>;

    gen_pt_ = new std::vector<float>;
    gen_pt_Recluster_ = new std::vector<float>;
    gen_pt_WithNu_ = new std::vector<float>;
    Delta_gen_pt_ = new std::vector<float>;
    Delta_gen_pt_Recluster_ = new std::vector<float>;
    Delta_gen_pt_WithNu_ = new std::vector<float>;

    y_multiplicity_ = new std::vector<float>;
    y_charged_multiplicity_ = new std::vector<float>;
    y_neutral_multiplicity_ = new std::vector<float>;
    y_ptD_ = new std::vector<float>;
    y_axis1_ = new std::vector<float>;
    y_axis2_ = new std::vector<float>;
    y_pt_dr_log_ = new std::vector<float>;

    // Add if statement later, like above?
    for (auto& entry : discriminators_) {
        entry.second = new std::vector<float>;
    }
}


void ntuple_JetInfo::clearContainers() {

    isB_->clear();
    isGBB_->clear();
    isBB_->clear();
    isC_->clear();
    isGCC_->clear();
    isCC_->clear();
    isUD_->clear();
    isS_->clear();
    isG_->clear();
    isUndefined_->clear();
    isLeptonicB_->clear();
    isLeptonicB_C_->clear();
    isTau_->clear();

    genDecay_->clear();

    isPhysB_->clear();
    isPhysGBB_->clear();
    isPhysBB_->clear();
    isPhysC_->clear();
    isPhysGCC_->clear();
    isPhysCC_->clear();
    isPhysUD_->clear();
    isPhysS_->clear();
    isPhysG_->clear();
    isPhysUndefined_->clear();
    isPhysLeptonicB_->clear();
    isPhysLeptonicB_C_->clear();
    isPhysTau_->clear();

    jet_qgl_->clear();
    QG_ptD_->clear();
    QG_axis2_->clear();
    QG_mult_->clear();

    muons_number_->clear();
    electrons_number_->clear();

    muons_jetIdx_->clear();
    muons_isLooseMuon_->clear();
    muons_isTightMuon_->clear();
    muons_isSoftMuon_->clear();
    muons_isHighPtMuon_->clear();
    muons_pt_->clear();
    muons_relEta_->clear();
    muons_relPhi_->clear();
    muons_energy_->clear();

    electrons_jetIdx_->clear();
    electrons_pt_->clear();
    electrons_relEta_->clear();
    electrons_relPhi_->clear();
    electrons_energy_->clear();

    jet_pt_->clear();
    jet_corr_pt_->clear();
    jet_eta_->clear();
    jet_phi_->clear();
    jet_mass_->clear();
    jet_energy_->clear();
    jet_looseId_->clear();

    gen_pt_->clear();
    gen_pt_Recluster_->clear();
    gen_pt_WithNu_->clear();
    Delta_gen_pt_->clear();
    Delta_gen_pt_Recluster_->clear();
    Delta_gen_pt_WithNu_->clear();

    y_multiplicity_->clear();
    y_charged_multiplicity_->clear();
    y_neutral_multiplicity_->clear();
    y_ptD_->clear();
    y_axis1_->clear();
    y_axis2_->clear();
    y_pt_dr_log_->clear();

    for (auto& entry : discriminators_) {
        entry.second->clear();
    }
}


void ntuple_JetInfo::deleteContainers() {

    delete isB_;
    delete isGBB_;
    delete isBB_;
    delete isC_;
    delete isGCC_;
    delete isCC_;
    delete isUD_;
    delete isS_;
    delete isG_;
    delete isUndefined_;
    delete isLeptonicB_;
    delete isLeptonicB_C_;
    delete isTau_;

    delete genDecay_;

    delete isPhysB_;
    delete isPhysGBB_;
    delete isPhysBB_;
    delete isPhysC_;
    delete isPhysGCC_;
    delete isPhysCC_;
    delete isPhysUD_;
    delete isPhysS_;
    delete isPhysG_;
    delete isPhysUndefined_;
    delete isPhysLeptonicB_;
    delete isPhysLeptonicB_C_;
    delete isPhysTau_;
    
    delete jet_qgl_;
    delete QG_ptD_;
    delete QG_axis2_;
    delete QG_mult_;

    delete muons_number_;
    delete electrons_number_;

    delete muons_jetIdx_;
    delete muons_isLooseMuon_;
    delete muons_isTightMuon_;
    delete muons_isSoftMuon_;
    delete muons_isHighPtMuon_;
    delete muons_pt_;
    delete muons_relEta_;
    delete muons_relPhi_;
    delete muons_energy_;

    delete electrons_jetIdx_;
    delete electrons_pt_;
    delete electrons_relEta_;
    delete electrons_relPhi_;
    delete electrons_energy_;

    delete jet_pt_;
    delete jet_corr_pt_;
    delete jet_eta_;
    delete jet_phi_;
    delete jet_mass_;
    delete jet_energy_;
    delete jet_looseId_;

    delete gen_pt_;
    delete gen_pt_Recluster_;
    delete gen_pt_WithNu_;
    delete Delta_gen_pt_;
    delete Delta_gen_pt_Recluster_;
    delete Delta_gen_pt_WithNu_;

    delete y_multiplicity_;
    delete y_charged_multiplicity_;
    delete y_neutral_multiplicity_;
    delete y_ptD_;
    delete y_axis1_;
    delete y_axis2_;
    delete y_pt_dr_log_;

    for (auto& entry : discriminators_) {
        delete entry.second;
    }
}


/*
// Use either of these functions
// bool ntuple_JetInfo::fillBranches(const pat::Jet& jet, const size_t& jetidx, const edm::View<pat::Jet>* coll) {
bool ntuple_JetInfo::fillBranches(const pat::Jet& jet, const size_t& jetidx, const edm::View<pat::Jet>* coll, float EventTime) {

    if (!coll) throw std::runtime_error("ntuple_JetInfo::fillBranches: no jet collection");

    // Cuts
    bool returnval = true;

    // Some cuts to constrain training region
    if (jet.pt() < jetPtMin_ || jet.pt() > jetPtMax_) returnval = false; // Apply jet pT cut
    if (fabs(jet.eta()) < jetAbsEtaMin_ || fabs(jet.eta()) > jetAbsEtaMax_) returnval = false; // Apply jet eta cut

    // Often we have way too many gluons that we do not need -- this randomly reduces the gluons
    if (gluonReduction_ > 0 && jet.partonFlavour() == 21)
        if (TRandom_.Uniform() > gluonReduction_) returnval = false;

    if (jet.genJet() == NULL) returnval = false;

    // Branch fills
    for (auto& entry : discriminators_) {
        entry.second = catchInfs(jet.bDiscriminator(entry.first), -0.1);
    }

    // npv_ = vertices()->size();

    // npv_0_z_ = vertices()->at(0).z();

    float PUrho = 0;
    std::vector<PileupSummaryInfo>::const_iterator ipu;
    for (ipu = PUInfo_->begin(); ipu != PUInfo_->end(); ++ipu) {
        if (ipu->getBunchCrossing() != 0) continue; // storing detailed PU info only for BX=0

        for (unsigned int i = 0; i < ipu->getPU_zpositions().size(); ++i) {
            auto PU_z = (ipu->getPU_zpositions())[i];
            // std::cout << i << " " << PU_z << std::endl;
            if (std::abs(PU_z - npv_0_z_) < 1) PUrho++;
        }
    }
    PUrho /= 20.0;

    // std::cout << "====================" << std::endl; 
    // std::cout << "PUrho: " << PUrho << std::endl;
    // std::cout << "====================" << std::endl; 

    PU_rho_ = PUrho;

    // float nPU = 200;
    // float cons = 0.09319;
    // float mean = -0.0045;
    // float sig = 4.28;
    // float numPU = nPU;

    // TF1 fungaus("fungaus", "gaus", -25.0, 25.0);
    // fungaus.SetParameters(cons, mean, sig);

    // float PUrho = fungaus.Integral(npv_0_z_ - 0.05, npv_0_z_ + 0.05);
    // PU_rho_ = PUrho * numPU;

    for (auto const& v : *pupInfo()) {
        int bx = v.getBunchCrossing();
        if (bx == 0) {
            ntrueInt_ = v.getTrueNumInteractions();
        }
    }
    rho_ = rhoInfo()[0];

    jet_no_ = jetidx;

    const auto jetRef = reco::CandidatePtr(coll->ptrs().at(jetidx));

    jet_qgl_ = (*qglHandle_)[jetRef];
    QG_ptD_ = (*ptDHandle_)[jetRef];
    QG_axis2_ = (*axis2Handle_)[jetRef];
    QG_mult_ = (*multHandle_)[jetRef];

    // std::vector<Ptr<pat::Jet>> p = coll->ptrs();

    isB_ = 0; isGBB_ = 0; isBB_ = 0; isC_ = 0; isGCC_ = 0; isCC_ = 0; isUD_ = 0; isTau_ = 0;
    isS_ = 0; isG_ = 0; isLeptonicB_ = 0; isLeptonicB_C_ = 0; isUndefined_ = 0;
    auto muIds = deep_ntuples::jet_muonsIds(jet, *muonsHandle_);
    auto elecIds = deep_ntuples::jet_electronsIds(jet, *electronsHandle_);

    muons_number_ = muIds.size();
    electrons_number_ = elecIds.size();

    float etasign = 1.0;
    if (jet.eta() < 0) etasign = -1.0;

    for (std::size_t i = 0; i < max_num_lept_; i++) {
        if (i < muIds.size()) {
            const auto& muon = (*muonsHandle_).at(muIds.at(i));
            muons_isLooseMuon_[i] = muon.isLooseMuon();
            muons_isTightMuon_[i] = muon.isTightMuon(vertices()->at(0));
            muons_isSoftMuon_[i] = muon.isSoftMuon(vertices()->at(0));
            muons_isHighPtMuon_[i] = muon.isHighPtMuon(vertices()->at(0));
            muons_pt_[i] = muon.pt();
            muons_relEta_[i] = etasign * (muon.eta() - jet.eta());
            muons_relPhi_[i] = reco::deltaPhi(muon.phi(), jet.phi());
            muons_energy_[i] = muon.energy() / jet.energy();
        }
        if (i < elecIds.size()) {
            const auto& electron = (*electronsHandle_).at(elecIds.at(i));
            electrons_pt_[i] = electron.pt();
            electrons_relEta_[i] = etasign * (electron.eta() - jet.eta());
            electrons_relPhi_[i] = reco::deltaPhi(electron.phi(), jet.phi());
            electrons_energy_[i] = electron.energy() / jet.energy();
        }
    }

    // Note that jets with gluon->bb (cc) and x->bb (cc) are in the same categories
    if (jet.genJet() != NULL) {
        switch(deep_ntuples::jet_flavour(jet, gToBB_, gToCC_, neutrinosLepB_, neutrinosLepB_C_, alltaus_)) {
            case deep_ntuples::JetFlavor::B:           isB_ = 1;           break;
            case deep_ntuples::JetFlavor::LeptonicB:   isLeptonicB_ = 1;   break;
            case deep_ntuples::JetFlavor::LeptonicB_C: isLeptonicB_C_ = 1; break;
            case deep_ntuples::JetFlavor::GBB:         isGBB_ = 1;         break;
            case deep_ntuples::JetFlavor::BB:          isBB_ = 1;          break;
            case deep_ntuples::JetFlavor::C:           isC_ = 1;           break;
            case deep_ntuples::JetFlavor::GCC:         isGCC_ = 1;         break;
            case deep_ntuples::JetFlavor::CC:          isCC_ = 1;          break;
            case deep_ntuples::JetFlavor::TAU:         isTau_ = 1;         break;
            case deep_ntuples::JetFlavor::G:           isG_ = 1;           break;
            case deep_ntuples::JetFlavor::UD:          isUD_ = 1;          break;
            case deep_ntuples::JetFlavor::S:           isS_ = 1;           break;
            default:                                   isUndefined_ = 1;   break;
        }
    }

    // Truth labeling with fallback to physics definition for light/gluon/undefined of standard flavor definition
    // Note that jets with gluon->bb (cc) and x->bb (cc) are in the same categories
    isPhysB_ = 0; isPhysBB_ = 0; isPhysGBB_ = 0; isPhysC_ = 0; isPhysCC_ = 0;
    isPhysGCC_ = 0; isPhysUD_ = 0; isPhysS_ = 0; isPhysG_ = 0, isPhysLeptonicB_ = 0, isPhysLeptonicB_C_ = 0, isPhysUndefined_ = 0;
    isPhysTau_ = 0;
    if (jet.genJet() != NULL) {
        switch(deep_ntuples::jet_flavour(jet, gToBB_, gToCC_, neutrinosLepB_, neutrinosLepB_C_, alltaus_, true)) {
            case deep_ntuples::JetFlavor::UD:          isPhysUD_ = 1;          break;
            case deep_ntuples::JetFlavor::S:           isPhysS_ = 1;           break;
            case deep_ntuples::JetFlavor::B:           isPhysB_ = 1;           break;
            case deep_ntuples::JetFlavor::BB:          isPhysBB_ = 1;          break;
            case deep_ntuples::JetFlavor::GBB:         isPhysGBB_ = 1;         break;
            case deep_ntuples::JetFlavor::C:           isPhysC_ = 1;           break;
            case deep_ntuples::JetFlavor::CC:          isPhysCC_ = 1;          break;
            case deep_ntuples::JetFlavor::GCC:         isPhysGCC_ = 1;         break;
            case deep_ntuples::JetFlavor::TAU:         isPhysTau_ = 1;         break;
            case deep_ntuples::JetFlavor::G:           isPhysG_ = 1;           break;
            case deep_ntuples::JetFlavor::LeptonicB:   isPhysLeptonicB_ = 1;   break;
            case deep_ntuples::JetFlavor::LeptonicB_C: isPhysLeptonicB_C_ = 1; break;
            default:                                   isPhysUndefined_ = 1;   break;
        }
    }

    if (!jet.genJet()) { // for data
        isUndefined_ = 1; isPhysUndefined_ = 1;
    }

    // Skip event, if neither standard flavor definition nor physics definition fallback define a "proper flavor"
    if (isUndefined_ && isPhysUndefined_) returnval = false;

    pat::JetCollection h;

    jet_pt_ = jet.correctedJet("Uncorrected").pt();
    jet_eta_ = jet.eta();
    jet_phi_ = jet.phi();
    jet_corr_pt_ = jet.pt();
    jet_mass_ = jet.mass();
    jet_energy_ = jet.energy();

    genDecay_ = -1.0;

    try {
        reco::GenParticleRefVector Bhadrons_in_jet = jet.jetFlavourInfo().getbHadrons();

        if (Bhadrons_in_jet.size() > 0) { 
            for (unsigned int idx = 0; idx < Bhadron_.size(); ++idx) {
                reco::GenParticle bhad = Bhadron_[idx];
                bool bhad_is_in_jet = false;

                for (reco::GenParticleRefVector::const_iterator bhad_in_jet = Bhadrons_in_jet.begin(); bhad_in_jet != Bhadrons_in_jet.end(); ++bhad_in_jet) {
                    // Check if bhad is identical to bhad_in_jet
                    if ((*bhad_in_jet)->pt() == bhad.pt() && (*bhad_in_jet)->eta() == bhad.eta()
                            && (*bhad_in_jet)->phi() == bhad.phi() && (*bhad_in_jet)->pdgId() == bhad.pdgId())
                        bhad_is_in_jet = true;
                }
                if (bhad_is_in_jet) {
                    if (Bhadron_daughter_[idx].vx() != bhad.vx()) {
                        float vx = Bhadron_daughter_[idx].vx() - bhad.vx();
                        float vy = Bhadron_daughter_[idx].vy() - bhad.vy();

                        float dxy = sqrt(vx * vx + vy * vy);
                        if (dxy > genDecay_)
                            genDecay_ = dxy;
                    }
                    else if (genDecay_ < 0) 
                        genDecay_ = -0.1;
                }
            }
        }
    }
    catch (const cms::Exception& e) {
        genDecay_ = -1.0;
    }

    // https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
    try {
        float NHF = jet.neutralHadronEnergyFraction();
        float NEMF = jet.neutralEmEnergyFraction();
        float CHF = jet.chargedHadronEnergyFraction();
        // float MUF = jet.muonEnergyFraction();
        float CEMF = jet.chargedEmEnergyFraction();
        float NumConst = jet.chargedMultiplicity() + jet.neutralMultiplicity();
        float NumNeutralParticles = jet.neutralMultiplicity();
        float CHM = jet.chargedMultiplicity();

        jet_looseId_ = ((NHF < 0.99 && NEMF < 0.99 && NumConst > 1)
                        && ((abs(jet_eta_) <= 2.4 && CHF > 0 && CHM > 0 && CEMF < 0.99) || abs(jet_eta_) > 2.4)
                        && abs(jet_eta_) <= 2.7)
                || (NHF < 0.98 && NEMF > 0.01 && NumNeutralParticles > 2 && abs(jet_eta_) > 2.7 && abs(jet_eta_) <= 3.0)
                || (NEMF < 0.90 && NumNeutralParticles > 10 && abs(jet_eta_) > 3.0);
    }
    catch(const cms::Exception& e) {
        jet_looseId_ = 1;
    }

    gen_pt_ = 0;
    Delta_gen_pt_ = 0;
    gen_pt_Recluster_ = 0;
    gen_pt_WithNu_ = 0;
    Delta_gen_pt_Recluster_ = 0;
    Delta_gen_pt_WithNu_ = 0;

    if (jet.genJet()) {
        gen_pt_ = jet.genJet()->pt();
        Delta_gen_pt_ = jet.genJet()->pt() - jet_pt_;

        const edm::RefToBase<pat::Jet> patJetRef = coll->refAt(jetidx);
        reco::GenJetRef genjetRecluster = (*genJetMatchRecluster_)[patJetRef];

        gen_pt_Recluster_ = 0.0;
        if (genjetRecluster.isNonnull() && genjetRecluster.isAvailable()) {
            gen_pt_Recluster_ = genjetRecluster->pt();
        }
        reco::GenJetRef genjetWithNu = (*genJetMatchWithNu_)[patJetRef];

        gen_pt_WithNu_ = 0.0;
        if (genjetWithNu.isNonnull() && genjetWithNu.isAvailable()) {
            gen_pt_WithNu_ = genjetWithNu->pt();
        }

        Delta_gen_pt_Recluster_ = gen_pt_Recluster_ - jet.pt();
        Delta_gen_pt_WithNu_ = gen_pt_WithNu_ - jet.pt();
    }

    auto qgtuple = yuta::calcVariables(&jet);
    // (multiplicity, charged_multiplicity, neutral_multiplicity, ptD, axis1, axis2, pt_dr_log);

    y_multiplicity_ = std::get<0>(qgtuple);
    y_charged_multiplicity_ = std::get<1>(qgtuple);
    y_neutral_multiplicity_ = std::get<2>(qgtuple);
    y_ptD_ = std::get<3>(qgtuple);
    y_axis1_ = std::get<4>(qgtuple);
    y_axis2_ = std::get<5>(qgtuple);
    y_pt_dr_log_ = std::get<6>(qgtuple);

    return returnval;
}
*/


// Returns number of selected jets after cuts
int ntuple_JetInfo::fillBranches(bool applySelection, float EventTime) {

    std::cout << "the other fillBranches() in ntuple_JetInfo" << std::endl;

    clearContainers();

    const edm::View<pat::Jet>* jetCollection = jets();

    if (jetCollection->size() == 0) throw std::runtime_error("ntuple_JetInfo::fillBranches (other): no jet collection");

    nJets_ = (int) jetCollection->size();
    npv_ = vertices()->size();
    npv_0_z_ = vertices()->at(0).z();

    float PUrho = 0;
    std::vector<PileupSummaryInfo>::const_iterator ipu;
    for (ipu = PUInfo_->begin(); ipu != PUInfo_->end(); ++ipu) {
        if (ipu->getBunchCrossing() != 0) continue; // storing detailed PU info only for BX=0

        for (unsigned int i = 0; i < ipu->getPU_zpositions().size(); ++i) {
            auto PU_z = (ipu->getPU_zpositions())[i];
            if (std::abs(PU_z - npv_0_z_) < 1) PUrho++;
        }
    }
    PUrho /= 20.0;
    PU_rho_ = PUrho;

    for (auto const& v : *pupInfo()) {
        int bx = v.getBunchCrossing();
        if (bx == 0) {
            ntrueInt_ = v.getTrueNumInteractions();
        }
    }
    rho_ = rhoInfo()[0];

    std::vector<size_t> indices(nJets_);
    for (size_t i = 0; i < (unsigned int) nJets_; i++)
        indices.at(i) = i;
    if (applySelection)
        std::random_shuffle(indices.begin(), indices.end());

    int nJetsSelected = 0;
    edm::View<pat::Jet>::const_iterator jetIter;
    for (size_t jetIdx : indices) {

        std::cout << "jetIdx = " << jetIdx << std::endl;

        jetIter = jetCollection->begin() + jetIdx;
        const pat::Jet& jet = *jetIter;

        // Some cuts to constrain training region
        if (jet.pt() < jetPtMin_ || jet.pt() > jetPtMax_) continue; // Apply jet pT cut
        if (fabs(jet.eta()) < jetAbsEtaMin_ || fabs(jet.eta()) > jetAbsEtaMax_) continue; // Apply jet eta cut

        // Often we have way too many gluons that we do not need -- this randomly reduces the gluons
        if (gluonReduction_ > 0 && jet.partonFlavour() == 21) {
            if (TRandom_.Uniform() > gluonReduction_) continue;
        }

        if (jet.genJet() == NULL) continue;

        // Note that jets with gluon->bb (cc) and x->bb (cc) are in the same categories
        int isB = 0;
        int isGBB = 0;
        int isBB = 0;
        int isC = 0;
        int isGCC = 0;
        int isCC = 0;
        int isUD = 0;
        int isTau = 0;
        int isS = 0;
        int isG = 0;
        int isLeptonicB = 0;
        int isLeptonicB_C = 0;
        int isUndefined = 0;
        if (jet.genJet() != NULL) { // Maybe this if condition is redundant?
            switch(deep_ntuples::jet_flavour(jet, gToBB_, gToCC_, neutrinosLepB_, neutrinosLepB_C_, alltaus_)) {
                case deep_ntuples::JetFlavor::B:           isB = 1;           break;
                case deep_ntuples::JetFlavor::LeptonicB:   isLeptonicB = 1;   break;
                case deep_ntuples::JetFlavor::LeptonicB_C: isLeptonicB_C = 1; break;
                case deep_ntuples::JetFlavor::GBB:         isGBB = 1;         break;
                case deep_ntuples::JetFlavor::BB:          isBB = 1;          break;
                case deep_ntuples::JetFlavor::C:           isC = 1;           break;
                case deep_ntuples::JetFlavor::GCC:         isGCC = 1;         break;
                case deep_ntuples::JetFlavor::CC:          isCC = 1;          break;
                case deep_ntuples::JetFlavor::TAU:         isTau = 1;         break;
                case deep_ntuples::JetFlavor::G:           isG = 1;           break;
                case deep_ntuples::JetFlavor::UD:          isUD = 1;          break;
                case deep_ntuples::JetFlavor::S:           isS = 1;           break;
                default:                                   isUndefined = 1;   break;
            }
        }

        // Truth labeling with fallback to physics definition for light/gluon/undefined of standard flavor definition
        // Note that jets with gluon->bb (cc) and x->bb (cc) are in the same categories
        int isPhysB = 0;
        int isPhysBB = 0;
        int isPhysGBB = 0;
        int isPhysC = 0;
        int isPhysCC = 0;
        int isPhysGCC = 0;
        int isPhysUD = 0;
        int isPhysS = 0;
        int isPhysG = 0;
        int isPhysLeptonicB = 0;
        int isPhysLeptonicB_C = 0;
        int isPhysUndefined = 0;
        int isPhysTau = 0;
        if (jet.genJet() != NULL) { // Maybe this if condition is redundant?
            switch(deep_ntuples::jet_flavour(jet, gToBB_, gToCC_, neutrinosLepB_, neutrinosLepB_C_, alltaus_, true)) {
                case deep_ntuples::JetFlavor::UD:          isPhysUD = 1;          break;
                case deep_ntuples::JetFlavor::S:           isPhysS = 1;           break;
                case deep_ntuples::JetFlavor::B:           isPhysB = 1;           break;
                case deep_ntuples::JetFlavor::BB:          isPhysBB = 1;          break;
                case deep_ntuples::JetFlavor::GBB:         isPhysGBB = 1;         break;
                case deep_ntuples::JetFlavor::C:           isPhysC = 1;           break;
                case deep_ntuples::JetFlavor::CC:          isPhysCC = 1;          break;
                case deep_ntuples::JetFlavor::GCC:         isPhysGCC = 1;         break;
                case deep_ntuples::JetFlavor::TAU:         isPhysTau = 1;         break;
                case deep_ntuples::JetFlavor::G:           isPhysG = 1;           break;
                case deep_ntuples::JetFlavor::LeptonicB:   isPhysLeptonicB = 1;   break;
                case deep_ntuples::JetFlavor::LeptonicB_C: isPhysLeptonicB_C = 1; break;
                default:                                   isPhysUndefined = 1;   break;
            }
        }

        if (!jet.genJet()) { // for data
            isUndefined = 1;
            isPhysUndefined = 1;
        }

        // Skip event, if neither standard flavor definition nor physics definition fallback define a "proper flavor"
        if (isUndefined && isPhysUndefined) continue;

        isB_->push_back(isB);
        isGBB_->push_back(isGBB);
        isBB_->push_back(isBB);
        isC_->push_back(isC);
        isGCC_->push_back(isGCC);
        isCC_->push_back(isCC);
        isUD_->push_back(isUD);
        isTau_->push_back(isTau);
        isS_->push_back(isS);
        isG_->push_back(isG);
        isLeptonicB_->push_back(isLeptonicB);
        isLeptonicB_C_->push_back(isLeptonicB_C);
        isUndefined_->push_back(isUndefined);

        isPhysB_->push_back(isPhysB);
        isPhysBB_->push_back(isPhysBB);
        isPhysGBB_->push_back(isPhysGBB);
        isPhysC_->push_back(isPhysC);
        isPhysCC_->push_back(isPhysCC);
        isPhysGCC_->push_back(isPhysGCC);
        isPhysUD_->push_back(isPhysUD);
        isPhysS_->push_back(isPhysS);
        isPhysG_->push_back(isPhysG);
        isPhysLeptonicB_->push_back(isPhysLeptonicB);
        isPhysLeptonicB_C_->push_back(isPhysLeptonicB_C);
        isPhysTau_->push_back(isPhysTau);
        isPhysUndefined_->push_back(isPhysUndefined);

        // Branch fills
        for (auto& entry : discriminators_) {
            entry.second->push_back(catchInfs(jet.bDiscriminator(entry.first), -0.1));
        }

        const auto jetRef = reco::CandidatePtr(jetCollection->ptrs().at(jetIdx));
        jet_qgl_->push_back((*qglHandle_)[jetRef]);
        QG_ptD_->push_back((*ptDHandle_)[jetRef]);
        QG_axis2_->push_back((*axis2Handle_)[jetRef]);
        QG_mult_->push_back((*multHandle_)[jetRef]);

        auto muIds = deep_ntuples::jet_muonsIds(jet, *muonsHandle_);
        auto elecIds = deep_ntuples::jet_electronsIds(jet, *electronsHandle_);

        muons_number_->push_back(muIds.size());
        electrons_number_->push_back(elecIds.size());

        float etasign = 1.0;
        if (jet.eta() < 0) etasign = -1.0;

        for (std::size_t i = 0; i < max_num_lept_; i++) {
            if (i < muIds.size()) {
                const auto& muon = (*muonsHandle_).at(muIds.at(i));
                muons_jetIdx_->push_back(jetIdx);
                muons_isLooseMuon_->push_back(muon.isLooseMuon());
                muons_isTightMuon_->push_back(muon.isTightMuon(vertices()->at(0)));
                muons_isSoftMuon_->push_back(muon.isSoftMuon(vertices()->at(0)));
                muons_isHighPtMuon_->push_back(muon.isHighPtMuon(vertices()->at(0)));
                muons_pt_->push_back(muon.pt());
                muons_relEta_->push_back(etasign * (muon.eta() - jet.eta()));
                muons_relPhi_->push_back(reco::deltaPhi(muon.phi(), jet.phi()));
                muons_energy_->push_back(muon.energy() / jet.energy());
            }
            if (i < elecIds.size()) {
                const auto& electron = (*electronsHandle_).at(elecIds.at(i));
                electrons_jetIdx_->push_back(jetIdx);
                electrons_pt_->push_back(electron.pt());
                electrons_relEta_->push_back(etasign * (electron.eta() - jet.eta()));
                electrons_relPhi_->push_back(reco::deltaPhi(electron.phi(), jet.phi()));
                electrons_energy_->push_back(electron.energy() / jet.energy());
            }
        }

        float jet_pt = jet.correctedJet("Uncorrected").pt();
        float jet_eta = jet.eta();
        float jet_corr_pt = jet.pt();
        jet_pt_->push_back(jet_pt);
        jet_eta_->push_back(jet_eta);
        jet_phi_->push_back(jet.phi());
        jet_corr_pt_->push_back(jet_corr_pt);
        jet_mass_->push_back(jet.mass());
        jet_energy_->push_back(jet.energy());

        float genDecay = -1.0;
        try {
            reco::GenParticleRefVector Bhadrons_in_jet = jet.jetFlavourInfo().getbHadrons();

            if (Bhadrons_in_jet.size() > 0) { 
                for (unsigned int idx = 0; idx < Bhadron_.size(); ++idx) {
                    reco::GenParticle bhad = Bhadron_[idx];
                    bool bhad_is_in_jet = false;

                    for (reco::GenParticleRefVector::const_iterator bhad_in_jet = Bhadrons_in_jet.begin(); bhad_in_jet != Bhadrons_in_jet.end(); ++bhad_in_jet) {
                        // Check if bhad is identical to bhad_in_jet
                        if ((*bhad_in_jet)->pt() == bhad.pt() && (*bhad_in_jet)->eta() == bhad.eta()
                                && (*bhad_in_jet)->phi() == bhad.phi() && (*bhad_in_jet)->pdgId() == bhad.pdgId())
                            bhad_is_in_jet = true;
                    }
                    if (bhad_is_in_jet) {
                        if (Bhadron_daughter_[idx].vx() != bhad.vx()) {
                            float vx = Bhadron_daughter_[idx].vx() - bhad.vx();
                            float vy = Bhadron_daughter_[idx].vy() - bhad.vy();

                            float dxy = sqrt(vx * vx + vy * vy);
                            if (dxy > genDecay)
                                genDecay = dxy;
                        }
                        else if (genDecay < 0) 
                            genDecay = -0.1;
                    }
                }
            }
        }
        catch (const cms::Exception& e) {
            genDecay = -1.0;
        }
        genDecay_->push_back(genDecay);

        // https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
        float jet_looseId;
        try {
            float NHF = jet.neutralHadronEnergyFraction();
            float NEMF = jet.neutralEmEnergyFraction();
            float CHF = jet.chargedHadronEnergyFraction();
            // float MUF = jet.muonEnergyFraction();
            float CEMF = jet.chargedEmEnergyFraction();
            float NumConst = jet.chargedMultiplicity() + jet.neutralMultiplicity();
            float NumNeutralParticles = jet.neutralMultiplicity();
            float CHM = jet.chargedMultiplicity();

            jet_looseId = ((NHF < 0.99 && NEMF < 0.99 && NumConst > 1)
                            && ((abs(jet_eta) <= 2.4 && CHF > 0 && CHM > 0 && CEMF < 0.99) || abs(jet_eta) > 2.4)
                            && abs(jet_eta) <= 2.7)
                    || (NHF < 0.98 && NEMF > 0.01 && NumNeutralParticles > 2 && abs(jet_eta) > 2.7 && abs(jet_eta) <= 3.0)
                    || (NEMF < 0.90 && NumNeutralParticles > 10 && abs(jet_eta) > 3.0);
        }
        catch(const cms::Exception& e) {
            jet_looseId = 1;
        }
        jet_looseId_->push_back(jet_looseId);

        float gen_pt_Recluster = 0.0;
        float gen_pt_WithNu = 0.0;
        if (jet.genJet()) { // Maybe this if condition is redundant?
            gen_pt_->push_back(jet.genJet()->pt());
            Delta_gen_pt_->push_back(jet.genJet()->pt() - jet_pt);

            const edm::RefToBase<pat::Jet> patJetRef = jetCollection->refAt(jetIdx);
            reco::GenJetRef genjetRecluster = (*genJetMatchRecluster_)[patJetRef];

            gen_pt_Recluster = 0.0;
            if (genjetRecluster.isNonnull() && genjetRecluster.isAvailable()) {
                gen_pt_Recluster = genjetRecluster->pt();
            }
            reco::GenJetRef genjetWithNu = (*genJetMatchWithNu_)[patJetRef];

            gen_pt_WithNu = 0.0;
            if (genjetWithNu.isNonnull() && genjetWithNu.isAvailable()) {
                gen_pt_WithNu = genjetWithNu->pt();
            }

            gen_pt_Recluster_->push_back(gen_pt_Recluster);
            gen_pt_WithNu_->push_back(gen_pt_WithNu);
            Delta_gen_pt_Recluster_->push_back(gen_pt_Recluster - jet_corr_pt);
            Delta_gen_pt_WithNu_->push_back(gen_pt_WithNu - jet_corr_pt);
        }

        auto qgtuple = yuta::calcVariables(&jet);
        // (multiplicity, charged_multiplicity, neutral_multiplicity, ptD, axis1, axis2, pt_dr_log);
        y_multiplicity_->push_back(std::get<0>(qgtuple));
        y_charged_multiplicity_->push_back(std::get<1>(qgtuple));
        y_neutral_multiplicity_->push_back(std::get<2>(qgtuple));
        y_ptD_->push_back(std::get<3>(qgtuple));
        y_axis1_->push_back(std::get<4>(qgtuple));
        y_axis2_->push_back(std::get<5>(qgtuple));
        y_pt_dr_log_->push_back(std::get<6>(qgtuple));

        nJetsSelected++;
    } // End of loop over jets

    return nJetsSelected;
}
