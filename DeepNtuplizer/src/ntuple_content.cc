/*
 * ntuple_content.cc
 *
 *  Created on: 13 Feb 2017
 *      Author: jkiesele
 */


#include "../interface/ntuple_content.h"
#include <stdexcept>


bool ntuple_content::useoffsets = true;


ntuple_content::~ntuple_content() {}


const reco::VertexCollection* ntuple_content::vertices() const {

    if (vertices_) return vertices_;
    throw std::runtime_error("ntuple_content: vertices not assigned");
}


const std::vector<reco::VertexCompositePtrCandidate>* ntuple_content::secVertices() const {

    if (secvertices_) return secvertices_;
    throw std::runtime_error("ntuple_content: secvertices not assigned");
}


const edm::View<pat::Jet>* ntuple_content::jets() const {

    return jets_;
}


const std::vector<PileupSummaryInfo>* ntuple_content::pupInfo() const {

    return pupInfo_;
}


const double* ntuple_content::rhoInfo() const {

    return rhoInfo_;
}
