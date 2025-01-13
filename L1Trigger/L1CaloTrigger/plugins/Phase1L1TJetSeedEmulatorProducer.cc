// -*- C++ -*-
//
// Package: L1CaloTrigger
// Class: Phase1L1TJetSeedEmulatorProducer
//
/**\class Phase1L1TJetSeedEmulatorProducer Phase1L1TJetSeedEmulatorProducer.cc L1Trigger/L1CaloTrigger/plugin/Phase1L1TJetSeedEmulatorProducer.cc
*/

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1TParticleFlow/interface/PFCluster.h"
#include "DataFormats/L1Trigger/interface/L1Candidate.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/L1TParticleFlow/interface/puppi.h"
#include "DataFormats/L1TParticleFlow/interface/gt_datatypes.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/common/bitonic_hybrid_sort_ref.h"
#include "TH2F.h"

#include <cmath>

#include <algorithm>
#include "L1Trigger/L1CaloTrigger/interface/Phase1L1TJetSeedEmulator.h"

class Phase1L1TJetSeedEmulatorProducer : public edm::one::EDProducer<> {
public:
  explicit Phase1L1TJetSeedEmulatorProducer(const edm::ParameterSet&);
  ~Phase1L1TJetSeedEmulatorProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  

  edm::EDGetTokenT<edm::View<reco::Candidate>> inputCollectionTag_;
  // histogram containing our clustered inputs
  
  bool debug;
  std::vector<double> etaBinning;
  size_t nBinsEta;
  unsigned int nBinsPhi;
  unsigned int jetIEtaSize;
  unsigned int jetIPhiSize;
  bool trimmedGrid;
  double seedPtThreshold;
  double ptlsb;
  double philsb;
  double etalsb;
  // Eta and phi edges of input PF regions
  std::vector<double> etaRegionEdges;
  std::vector<double> phiRegionEdges;
  // Maximum number of candidates per input PF region
  unsigned int maxInputsPerRegion;
  std::unique_ptr<TH2F> caloGrid_;
  Phase1L1TJetSeedEmulator emulator;
  std::string outputCollectionName;

};

Phase1L1TJetSeedEmulatorProducer::Phase1L1TJetSeedEmulatorProducer(const edm::ParameterSet& iConfig)
    : inputCollectionTag_{
          consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("inputCollectionTag"))},
      debug(iConfig.getParameter<bool>("debug")),
      etaBinning(iConfig.getParameter<std::vector<double>>("etaBinning")),
      nBinsEta(etaBinning.size() - 1),
      nBinsPhi(iConfig.getParameter<unsigned int>("nBinsPhi")),
      jetIEtaSize(iConfig.getParameter<unsigned int>("jetIEtaSize")),
      jetIPhiSize(iConfig.getParameter<unsigned int>("jetIPhiSize")),
      trimmedGrid(iConfig.getParameter<bool>("trimmedGrid")),
      seedPtThreshold(iConfig.getParameter<double>("seedPtThreshold")),
      ptlsb(iConfig.getParameter<double>("ptlsb")),
      philsb(iConfig.getParameter<double>("philsb")),
      etalsb(iConfig.getParameter<double>("etalsb")),
      etaRegionEdges(iConfig.getParameter<std::vector<double>>("etaRegions")),
      phiRegionEdges(iConfig.getParameter<std::vector<double>>("phiRegions")),
      maxInputsPerRegion(iConfig.getParameter<unsigned int>("maxInputsPerRegion")),
      caloGrid_(std::make_unique<TH2F>("caloGrid", "Calorimeter grid", nBinsEta, etaBinning.data(), nBinsPhi, phiRegionEdges.front(), phiRegionEdges.back())),
      emulator(debug, caloGrid_.get(), etaBinning, nBinsPhi, jetIEtaSize, jetIPhiSize, trimmedGrid, 
               seedPtThreshold, ptlsb, philsb, etalsb, etaRegionEdges, phiRegionEdges, maxInputsPerRegion),  // pass references and pointers
      outputCollectionName(iConfig.getParameter<std::string>("outputCollectionName")) {
    
    // Set axis titles (optional, but useful for visualization)
    caloGrid_->GetXaxis()->SetTitle("#eta");
    caloGrid_->GetYaxis()->SetTitle("#phi");

    // Register output collection
    produce<l1t::PFCandidateCollection>(outputCollectionName);
}

Phase1L1TJetSeedEmulatorProducer::~Phase1L1TJetSeedEmulatorProducer() {}


void Phase1L1TJetSeedEmulatorProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<reco::Candidate>> inputCollectionHandle;
  iEvent.getByToken(inputCollectionTag_, inputCollectionHandle);

  // sort inputs into PF regions
  std::vector<std::vector<reco::CandidatePtr>> inputsInRegions = prepareInputsIntoRegions<reco::Candidate>(inputCollectionHandle);

  // histogramming the data
  caloGrid_->Reset();
  for (unsigned int iInputRegion = 0; iInputRegion < inputsInRegions.size(); ++iInputRegion) {
    fillCaloGrid<reco::Candidate>(*(caloGrid_), inputsInRegions[iInputRegion], iInputRegion);
  }

  // find the seeds
  const auto& seedsVector = emulator.findSeeds(seedPtThreshold);  // seedPtThreshold = 5

  // sort by pt
  l1t::PFCandidateCollection sortedSeeds;
  emulator.sortSeeds( seedsVector, sortedSeeds );

  auto seedsVectorPtr = std::make_unique<l1t::PFCandidateCollection>(sortedSeeds);
  iEvent.put(std::move(seedsVectorPtr), outputCollectionName );

  return;
}



void Phase1L1TJetSeedEmulatorProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<bool>("debug", false);
  desc.add<edm::InputTag>("inputCollectionTag", edm::InputTag("l1pfCandidates", "Puppi"));
  desc.add<std::vector<double>>("etaBinning");
  desc.add<unsigned int>("nBinsPhi", 72);
  desc.add<unsigned int>("jetIEtaSize", 7);
  desc.add<unsigned int>("jetIPhiSize", 7);
  desc.add<bool>("trimmedGrid", false);
  desc.add<double>("seedPtThreshold", 5);
  desc.add<double>("ptlsb", 0.25), desc.add<double>("philsb", 0.0043633231), desc.add<double>("etalsb", 0.0043633231),
  desc.add<string>("outputCollectionName", "UncalibratedPhase1L1TJetFromPfCandidates");
  desc.add<std::vector<double>>("etaRegions");
  desc.add<std::vector<double>>("phiRegions");
  desc.add<unsigned int>("maxInputsPerRegion", 18);
  descriptions.add("Phase1L1TJetSeedEmulatorProducer", desc);
}

DEFINE_FWK_MODULE(Phase1L1TJetSeedEmulatorProducer);
