/*
 * ntuple_content.h
 *
 *  Created on: 13 Feb 2017
 *      Author: jkiesele
 */


#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_CONTENT_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_CONTENT_H_


#include <math.h>
#include <iostream>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TTree.h"
#include "TString.h"


// Base class for modules to inherit from
class ntuple_content {

    public:

        ntuple_content() : ntuple_content(0.4) {}
        // ntuple_content(double jetR) : vertices_(0), secvertices_(0), genvertices_(0), pupInfo_(0), rhoInfo_(0), jetR_(jetR), read_(false) {}
        ntuple_content(double jetR) : vertices_(0), secvertices_(0), pupInfo_(0), rhoInfo_(0), jetR_(jetR), read_(false) {}
        virtual ~ntuple_content();

        virtual void getInput(const edm::ParameterSet& iConfig) {}
        virtual void initBranches(TTree*) = 0;
        virtual void readEvent(const edm::Event& iEvent) = 0;
        virtual void readSetup(const edm::EventSetup& iSetup) {}
        virtual void initContainers() = 0;
        virtual void clearContainers() = 0;
        virtual void deleteContainers() = 0;

        // Use either of these functions
        // virtual bool fillBranches(const pat::Jet&, const size_t& jetidx, const edm::View<pat::Jet>* coll = 0) = 0;
        virtual bool fillBranches(const pat::Jet&, const size_t& jetidx, const edm::View<pat::Jet>* coll = 0, float EventTime = -1) = 0;
        virtual int fillBranches(bool applySelection, float EventTime = -1) = 0;

        void setPrimaryVertices(const reco::VertexCollection* v) { vertices_ = v; }
        void setSecVertices(const std::vector<reco::VertexCompositePtrCandidate>* v) { secvertices_ = v; }
        // void setGenVertices(const TrackingVertexCollection* v) { genvertices_ = v; }
        void setJets(const edm::View<pat::Jet>* v) { jets_ = v; }
        void setPuInfo(const std::vector<PileupSummaryInfo>* v) { pupInfo_ = v; }
        void setRhoInfo(const double* v) { rhoInfo_ = v; }

        void setIsRead(bool isread) { read_ = isread; }

        std::vector<TString> getListOfBranches() {

            if (allbranches_.size()) return allbranches_;
            else {
                TTree *t = new TTree();
                initBranches(t);
                return allbranches_;
            }
        }

        // std::vector<TString> getListOf2DBranches() {

        //     if (all2Dbranches_.size()) return all2Dbranches_;
        //     else {
        //         TTree *t = new TTree();
        //         initBranches(t);
        //         return all2Dbranches_;
        //     }
        // }

        static bool useoffsets;

    protected:

        const reco::VertexCollection* vertices() const;
        const std::vector<reco::VertexCompositePtrCandidate>* secVertices() const;
        // const TrackingVertexCollection* genVertices() const;
        const edm::View<pat::Jet>* jets() const;
        const std::vector<PileupSummaryInfo>* pupInfo() const;
        const double* rhoInfo() const;

        template <class T>
        void addBranch(TTree* t, const char* name, T*, const char* leaflist = 0);
        // template <class T>
        // void add2DBranch(TTree* t, const char* name, T*, const char* leaflist = 0);

        double jetR() const { return jetR_; }

        static inline const float& catchInfs(const float& in, const float& replace_value) {

            if (in == in) {
                if (std::isinf(in))
                    return replace_value;
                else if (in < -1e32 || in > 1e32)
                    return replace_value;
                return in;
            }
            return replace_value;
        }

        static inline float catchInfsAndBound(const float& in, const float& replace_value,
                const float& lowerbound, const float& upperbound, const float offset = 0) {

            float withoutinfs = catchInfs(in, replace_value);
            if (withoutinfs + offset < lowerbound) return lowerbound;
            if (withoutinfs + offset > upperbound) return upperbound;
            if (useoffsets)
                withoutinfs += offset;
            return withoutinfs;
        }

    private:

        const reco::VertexCollection* vertices_;
        const std::vector<reco::VertexCompositePtrCandidate>* secvertices_;
        // const TrackingVertexCollection* genvertices_;
        const edm::View<pat::Jet>* jets_;
        const std::vector<PileupSummaryInfo>* pupInfo_;
        const double* rhoInfo_;
        double jetR_;
        bool read_;

        std::vector<TString> allbranches_;
        // std::vector<TString> all2Dbranches_;
};


template <class T>
void ntuple_content::addBranch(TTree* t, const char* name, T* address, const char* leaflist) {

    if (read_) {
        t->SetBranchAddress(name, address);
    }
    else {
        if (leaflist)
            t->Branch(name, address, leaflist);
        else
            t->Branch(name, address);
    }
    allbranches_.push_back((TString) name);
}


// template <class T>
// void ntuple_content::add2DBranch(TTree* t, const char* name, T* address, const char* leaflist) {

//     if (read_) {
//         t->SetBranchAddress(name, address);
//     }
//     else {
//         if (leaflist)
//             t->Branch(name, address, leaflist);
//         else
//             t->Branch(name, address);
//     }
//     all2Dbranches_.push_back((TString) name);
// }


#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_CONTENT_H_ */
