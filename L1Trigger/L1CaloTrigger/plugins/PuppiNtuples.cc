// -*- C++ -*-
//
// Package:    L1Trigger/PuppiNtuples
// Class:      PuppiNtuples
//
/**\class PuppiNtuples PuppiNtuples.cc L1Trigger/PuppiNtuples/plugins/PuppiNtuples.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Emyr Clement
//         Created:  Thu, 26 Jan 2023 12:57:51 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


class PuppiNtuples : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit PuppiNtuples(const edm::ParameterSet&);
  ~PuppiNtuples() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  void clear();

  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<l1t::PFCandidate>> l1PFToken;

  // output file
  edm::Service<TFileService> fs_;

  // tree
  TTree* tree_;

  int puppi_n_;
  std::vector<float> puppi_pt_;
  std::vector<float> puppi_eta_;
  std::vector<float> puppi_phi_;
  std::vector<float> puppi_z0_;
  std::vector<float> puppi_dxy_;
  std::vector<int> puppi_hw_z0_;
  std::vector<int> puppi_hw_dxy_;
  std::vector<float> puppi_weight_;
  std::vector<int> puppi_pid_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PuppiNtuples::PuppiNtuples(const edm::ParameterSet& iConfig)
    :       l1PFToken(consumes<std::vector<l1t::PFCandidate>>(iConfig.getUntrackedParameter<edm::InputTag>("inputPFCands"))){
  //now do what ever initialization is needed
  // set up output
  tree_ = fs_->make<TTree>("Puppi", "Puppi");
  tree_->Branch("N", &puppi_n_);
  tree_->Branch("pt", &puppi_pt_);
  tree_->Branch("eta", &puppi_eta_);
  tree_->Branch("phi", &puppi_phi_);
  tree_->Branch("z0", &puppi_z0_);
  tree_->Branch("dxy", &puppi_dxy_);
  tree_->Branch("hwZ0", &puppi_hw_z0_);
  tree_->Branch("hwDxy", &puppi_hw_dxy_);
  tree_->Branch("puppiWeight", &puppi_weight_);
  tree_->Branch("pid", &puppi_pid_);
}

PuppiNtuples::~PuppiNtuples() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void PuppiNtuples::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  edm::Handle<l1t::PFCandidateCollection> l1PFCandidates;
  iEvent.getByToken(l1PFToken, l1PFCandidates);

  clear();

  puppi_n_ = l1PFCandidates->size();
  for (const auto& l1PFCand : *l1PFCandidates) {
    puppi_pt_.push_back( l1PFCand.pt() );
    puppi_eta_.push_back( l1PFCand.eta() );
    puppi_phi_.push_back( l1PFCand.phi() );
    puppi_z0_.push_back( l1PFCand.z0() );

    double d0 = 0;
    if (l1PFCand.pfTrack().isNonnull()) {
      d0 = std::hypot( l1PFCand.pfTrack()->vx(), l1PFCand.pfTrack()->vy() );
    }

    puppi_dxy_.push_back( d0 );
    puppi_hw_z0_.push_back( l1PFCand.hwZ0() );
    puppi_hw_dxy_.push_back( l1PFCand.hwDxy() );
    puppi_weight_.push_back( l1PFCand.puppiWeight() );
    puppi_pid_.push_back( l1PFCand.id() );
    // std::cout << "Puppi cand : " << l1PFCand.pt() << " " << l1PFCand.eta() << " " << l1PFCand.phi() << std::endl;
  }

  tree_->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void PuppiNtuples::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void PuppiNtuples::endJob() {
  // please remove this method if not needed
}

void PuppiNtuples::clear() {
  puppi_n_ = 0;
  puppi_pt_.clear();
  puppi_eta_.clear();
  puppi_phi_.clear();
  puppi_z0_.clear();
  puppi_dxy_.clear();
  puppi_hw_z0_.clear();
  puppi_hw_dxy_.clear();
  puppi_weight_.clear();
  puppi_pid_.clear();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void PuppiNtuples::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.addUntracked<edm::InputTag>("inputPFCands", edm::InputTag("l1ctLayer1","Puppi"));
  descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PuppiNtuples);

