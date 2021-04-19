#include <pybind11.h>
#include <stl.h>
#include <iostream.h>
#include <stdio.h>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"



//#include "interface/HelperFunctions.h"
//#include "interface/THStatus.h"
// #include "TTreeReader.h"
// #include "TTreeReaderValue.h"


namespace helper {
  template <typename T>
  void SetBranchAddress(TTree *tree, const char *bname, T **add, TBranch **ptr = 0);

  template <typename T>
  void SetBranchAddress(TTree *tree, const char *bname, T *add, TBranch **ptr = 0);
}

template <typename T>
void helper::SetBranchAddress(TTree *tree, const char *bname, T **add,
                              TBranch **ptr) {

  tree->SetBranchStatus(bname, 1);
  Int_t result = tree->SetBranchAddress(bname, add, ptr);
//  if (result != 0) {
//    std::stringstream sserr;
//    sserr << "Failed to set branch address for branch " << bname;
//    throw std::runtime_error(sserr.str());
//  }
}

template <typename T>
void helper::SetBranchAddress(TTree *tree, const char *bname, T *add,
                              TBranch **ptr) {

  tree->SetBranchStatus(bname, 1);
  Int_t result = tree->SetBranchAddress(bname, add, ptr);
//  if (result != 0) {
//    std::stringstream sserr;
//    sserr << "Failed to set branch address for branch " << bname;
//    throw std::runtime_error(sserr.str());
//  }
}

void rootf(std::string input_file_name) {

  std::cout << input_file_name << std::endl;
  TFile *infile = TFile::Open(input_file_name.c_str());

  TTree *intree = (TTree *)infile->Get("deepntuplizer/tree");

  auto njets = intree->GetEntries();
  std::cout << "njets: " << njets << std::endl;

  TFile *ofile = TFile::Open(("vertexCategory1Fix/" + input_file_name + ".tmp").c_str(), "RECREATE");
  TDirectory *odir = ofile->mkdir("deepntuplizer");
  odir->cd();

  TTree *otree = (TTree *)intree->CloneTree(0);

  Float_t TagVarCSV_flightDistance2dVal = -99;
  Float_t TagVarCSV_flightDistance2dSig = -99.;
  Float_t TagVarCSV_flightDistance3dVal = -99;
  Float_t TagVarCSV_flightDistance3dSig = -99;
  Float_t vertexCategory = -99;

  Float_t TagVarCSV_flightDistance2dVal_new = -99;
  Float_t TagVarCSV_flightDistance2dSig_new = -99;
  Float_t TagVarCSV_flightDistance3dVal_new = -99;
  Float_t TagVarCSV_flightDistance3dSig_new = -99;

  otree->Branch("TagVarCSV_flightDistance2dVal_new",
                &TagVarCSV_flightDistance2dVal_new,
                "TagVarCSV_flightDistance2dVal_new/F");
  otree->Branch("TagVarCSV_flightDistance2dSig_new",
                &TagVarCSV_flightDistance2dSig_new,
                "TagVarCSV_flightDistance2dSig_new/F");
  otree->Branch("TagVarCSV_flightDistance3dSig_new",
                &TagVarCSV_flightDistance3dSig_new,
                "TagVarCSV_flightDistance3dSig_new/F");
  otree->Branch("TagVarCSV_flightDistance3dVal_new",
                &TagVarCSV_flightDistance3dVal_new,
                "TagVarCSV_flightDistance3dVal_new/F");

  helper::SetBranchAddress(intree, "TagVarCSV_vertexCategory", &vertexCategory);
  helper::SetBranchAddress(intree, "TagVarCSV_flightDistance2dVal", &TagVarCSV_flightDistance2dVal);
  helper::SetBranchAddress(intree, "TagVarCSV_flightDistance2dSig", &TagVarCSV_flightDistance2dSig);
  helper::SetBranchAddress(intree, "TagVarCSV_flightDistance3dVal", &TagVarCSV_flightDistance3dVal);
  helper::SetBranchAddress(intree, "TagVarCSV_flightDistance3dSig", &TagVarCSV_flightDistance3dSig);


  for (int ijet = 0; ijet < njets; ++ijet) {

    intree->GetEntry(ijet);

    if (vertexCategory == 1 || vertexCategory == 2) {
      TagVarCSV_flightDistance2dVal_new = 0;
      TagVarCSV_flightDistance2dSig_new = 0;
      TagVarCSV_flightDistance3dVal_new = 0;
      TagVarCSV_flightDistance3dSig_new = 0;
    } else {
      TagVarCSV_flightDistance2dVal_new = TagVarCSV_flightDistance2dVal;
      TagVarCSV_flightDistance2dSig_new = TagVarCSV_flightDistance2dSig;
      TagVarCSV_flightDistance3dVal_new = TagVarCSV_flightDistance3dVal;
      TagVarCSV_flightDistance3dSig_new = TagVarCSV_flightDistance3dSig;
    }

    otree->Fill();
  }

  otree->Write("", TObject::kOverwrite);

  infile->Close();
  ofile->Close();

  std::string tmp_out = "vertexCategory1Fix/" + input_file_name + ".tmp";
  std::string new_out = "vertexCategory1Fix/" + input_file_name;
  char* tmp_out_char = new char [tmp_out.length()+1];
  char* new_out_char = new char [new_out.length()+1];
  std::strcpy (tmp_out_char, tmp_out.c_str());
  std::strcpy (new_out_char, new_out.c_str());

  rename(tmp_out_char, new_out_char);
}

void rootf_2(std::string input_file_name) {

  std::cout << input_file_name << std::endl;

  TFile *infile =
      TFile::Open(("vertexCategory1Fix/" + input_file_name).c_str());

  TTree *intree = (TTree *)infile->Get("deepntuplizer/tree");
  auto njets = intree->GetEntries();

  TFile *ofile = TFile::Open(
      ("vertexCategory1Fix/" + input_file_name + ".tmp").c_str(), "RECREATE");
  TDirectory *odir = ofile->mkdir("deepntuplizer");
  odir->cd();

  intree->SetBranchStatus("TagVarCSV_flightDistance2dVal", 0);
  intree->SetBranchStatus("TagVarCSV_flightDistance2dSig", 0);
  intree->SetBranchStatus("TagVarCSV_flightDistance3dVal", 0);
  intree->SetBranchStatus("TagVarCSV_flightDistance3dSig", 0);

  TTree *otree = (TTree *)intree->CloneTree(0);

  Float_t TagVarCSV_flightDistance2dVal = -99;
  Float_t TagVarCSV_flightDistance2dSig = -99.;
  Float_t TagVarCSV_flightDistance3dVal = -99;
  Float_t TagVarCSV_flightDistance3dSig = -99;

  otree->Branch("TagVarCSV_flightDistance2dVal", &TagVarCSV_flightDistance2dVal,
                "TagVarCSV_flightDistance2dVal/F");
  otree->Branch("TagVarCSV_flightDistance2dSig", &TagVarCSV_flightDistance2dSig,
                "TagVarCSV_flightDistance2dSig/F");
  otree->Branch("TagVarCSV_flightDistance3dSig", &TagVarCSV_flightDistance3dSig,
                "TagVarCSV_flightDistance3dSig/F");
  otree->Branch("TagVarCSV_flightDistance3dVal", &TagVarCSV_flightDistance3dVal,
                "TagVarCSV_flightDistance3dVal/F");

  helper::SetBranchAddress(intree, "TagVarCSV_flightDistance2dVal_new",
                           &TagVarCSV_flightDistance2dVal);
  helper::SetBranchAddress(intree, "TagVarCSV_flightDistance2dSig_new",
                           &TagVarCSV_flightDistance2dSig);
  helper::SetBranchAddress(intree, "TagVarCSV_flightDistance3dVal_new",
                           &TagVarCSV_flightDistance3dVal);
  helper::SetBranchAddress(intree, "TagVarCSV_flightDistance3dSig_new",
                           &TagVarCSV_flightDistance3dSig);


  for (int ijet = 0; ijet < njets; ++ijet) {

    intree->GetEntry(ijet);

    otree->Fill();
  }

  otree->Write("", TObject::kOverwrite);

  infile->Close();
  ofile->Close();

  std::string tmp_out = "vertexCategory1Fix/" + input_file_name + ".tmp";
  std::string new_out = "vertexCategory1Fix/" + input_file_name;
  char* tmp_out_char = new char [tmp_out.length()+1];
  char* new_out_char = new char [new_out.length()+1];
  std::strcpy (tmp_out_char, tmp_out.c_str());
  std::strcpy (new_out_char, new_out.c_str());

  rename(tmp_out_char, new_out_char);
}
namespace py = pybind11;

PYBIND11_MODULE(btv_vc1Fix, m) {
  m.def("rootf", &rootf);
  m.def("rootf_2", &rootf_2);
}
