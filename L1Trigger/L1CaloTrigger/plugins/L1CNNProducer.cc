// -*- C++ -*-
//
// Package: L1CaloTrigger
// Class: Phase1L1TJetProducer
//
/**\class Phase1L1TJetProducer Phase1L1TJetProducer.cc L1Trigger/L1CaloTrigger/plugin/Phase1L1TJetProducer.cc

Description: Produces jets with a phase-1 like sliding window algorithm using a collection of reco::Candidates in input.  Also calculates MET from the histogram used to find the jets.

*** INPUT PARAMETERS ***
  * nBinsPhi: uint32, number of bins in phi
  * phiLow: double, min phi (typically -pi)
  * phiUp: double, max phi (typically +pi)
  * jetIEtaSize: uint32, jet cluster size in ieta
  * jetIPhiSize: uint32, jet cluster size in iphi
  * trimmedGrid: Flag (bool) to remove three bins in each corner of grid in jet finding
  * seedPtThreshold: double, threshold of the seed tower
  * pt/eta/philsb : lsb of quantities used in firmware implementation
  * puSubtraction: bool, runs chunky doughnut pile-up subtraction, 9x9 jet only
  * eta/phiRegionEdges: Boundaries of the input (PF) regions
  * maxInputsPerRegion: Truncate number of inputes per input (PF) region
  * sin/cosPhi: Value of sin/cos phi in the middle of each bin of the grid.
  * met{HF}AbsETaCut: Eta selection of input candidates for calculation of MET
  * outputCollectionName: string, tag for the output collection
  * vetoZeroPt: bool, controls whether jets with 0 pt should be save. 
    It matters if PU is ON, as you can get negative or zero pt jets after it.
  * inputCollectionTag: inputtag, collection of reco::candidates used as input to the algo

*/
//
// Original Simone Bologna
// Created: Mon Jul 02 2018
// Modified 2020 Emyr Clement
//

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1TParticleFlow/interface/PFCluster.h"
#include "DataFormats/L1Trigger/interface/L1Candidate.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <iostream>
#include "TH2F.h"
#include "TFile.h"

#include <cmath>

#include <algorithm>
constexpr int x_scroll_min = -4;
constexpr int x_scroll_max = 4;
constexpr int y_scroll_min = 0;
constexpr int y_scroll_max = 3;

class CNNProducer : public edm::one::EDProducer<> {
public:
  explicit CNNProducer(const edm::ParameterSet&);
  ~CNNProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    void produce(edm::Event&, const edm::EventSetup&) override;

  /// Get the energy of a certain tower while correctly handling phi periodicity in case of overflow
  // float getTowerEnergy(int iEta, int iPhi) const;

  // <3 handy method to fill the calogrid with whatever type
  template <class Container>
  void fillCaloGrid(TH2F& caloGrid, const Container& triggerPrimitives, const unsigned int regionIndex);

  // Digitise the eta and phi coordinates of input candidates
  // This converts the quantities to integers to reduce precision
  // And takes account of bin edge effects i.e. makes sure the
  // candidate ends up in the correct (i.e. same behaviour as the firmware) bin of caloGrid_
  std::pair<float, float> getCandidateDigiEtaPhi(const float eta,
                                                 const float phi,
                                                 const unsigned int regionIndex) const;

  // Sorts the input candidates into the PF regions they arrive in
  // Truncates the inputs.  Takes the first N candidates as they are provided, without any sorting (this may be needed in the future and/or provided in this way from emulation of layer 1)
  template <class Handle>
  std::vector<std::vector<reco::CandidatePtr>> prepareInputsIntoRegions(const Handle triggerPrimitives);

  // Converts phi and eta (PF) region indices to a single index
  unsigned int getRegionIndex(const unsigned int phiRegion, const unsigned int etaRegion) const;
  // From the single index, calculated by getRegionIndex, provides the lower eta and phi boundaries of the input (PF) region index
  std::pair<double, double> regionEtaPhiLowEdges(const unsigned int regionIndex) const;
  // From the single index, calculated by getRegionIndex, provides the upper eta and phi boundaries of the input (PF) region index
  std::pair<double, double> regionEtaPhiUpEdges(const unsigned int regionIndex) const;

  void printHistogram(std::unique_ptr<TH2F>& histogram) const;

  edm::EDGetTokenT<edm::View<reco::Candidate>> inputCollectionTag_;
  edm::EDGetTokenT<std::vector<l1t::EtSum>> metTag_;
  edm::EDGetTokenT<std::vector<l1t::EtSum>> jetSumsTag_;

  // histogram containing our clustered inputs
  std::unique_ptr<TH2F> caloGrid_;
  std::unique_ptr<TH2F> caloGridWriteToFile_;
  std::unique_ptr<TH2F> caloGridTest_;

  double phiLow_;
  double phiUp_;
  double etaLow_;
  double etaUp_;
  double ptlsb_;
  double philsb_;
  double etalsb_;
  

  // Eta and phi edges of input PF regions
  std::vector<double> etaRegionEdges_;
  std::vector<double> phiRegionEdges_;
  // Maximum number of candidates per input PF region
  unsigned int maxInputsPerRegion_;
  std::string outputCollectionName_;

  bool runInference_;
  edm::FileInPath cnnGraph_;
  unsigned int cnnImageSizePhi_;
  unsigned int cnnImageSizeEta_;
  unsigned int cnnNPadPhi_;
  unsigned int cnnNPadEta_;
  float maxPixelValue_;

  bool writeImages_;
  unsigned int etaNBins_writeToFile_;
  double etaMin_writeToFile_;
  double etaMax_writeToFile_;
  unsigned int phiNBins_writeToFile_;
  double phiMin_writeToFile_;
  double phiMax_writeToFile_;
  std::string imagesFileName_;
  unsigned int imageCounter_;

  TFile* imageOutputFile_;

};

CNNProducer::CNNProducer(const edm::ParameterSet& iConfig)
    :  // getting configuration settings
      inputCollectionTag_{
          consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("inputCollectionTag"))},
      metTag_{consumes<std::vector<l1t::EtSum>>(iConfig.getParameter<edm::InputTag>("metTag"))},
      jetSumsTag_{consumes<std::vector<l1t::EtSum>>(iConfig.getParameter<edm::InputTag>("jetSumsTag"))},
      phiLow_(iConfig.getParameter<double>("phiLow")),
      phiUp_(iConfig.getParameter<double>("phiUp")),
      etaLow_(iConfig.getParameter<double>("etaLow")),
      etaUp_(iConfig.getParameter<double>("etaUp")),
      ptlsb_(iConfig.getParameter<double>("ptlsb")),
      philsb_(iConfig.getParameter<double>("philsb")),
      etalsb_(iConfig.getParameter<double>("etalsb")),
      etaRegionEdges_(iConfig.getParameter<std::vector<double>>("etaRegions")),
      phiRegionEdges_(iConfig.getParameter<std::vector<double>>("phiRegions")),
      maxInputsPerRegion_(iConfig.getParameter<unsigned int>("maxInputsPerRegion")),
      outputCollectionName_(iConfig.getParameter<std::string>("outputCollectionName")),
      runInference_(iConfig.getParameter<bool>("runInference")),
      cnnGraph_(iConfig.getParameter<edm::FileInPath>("cnnGraph")),
      cnnImageSizePhi_(iConfig.getParameter<unsigned int>("cnnImageSizePhi")),
      cnnImageSizeEta_(iConfig.getParameter<unsigned int>("cnnImageSizeEta")),
      cnnNPadPhi_(iConfig.getParameter<unsigned int>("cnnNPadPhi")),
      cnnNPadEta_(iConfig.getParameter<unsigned int>("cnnNPadEta")),
      maxPixelValue_(iConfig.getParameter<double>("maxPixelValue")),
      writeImages_(iConfig.getParameter<bool>("writeImages")),
      etaNBins_writeToFile_(iConfig.getParameter<unsigned int >("etaNBins_writeToFile")),
      etaMin_writeToFile_(iConfig.getParameter<double>("etaMin_writeToFile")),
      etaMax_writeToFile_(iConfig.getParameter<double>("etaMax_writeToFile")),
      phiNBins_writeToFile_(iConfig.getParameter<unsigned int >("phiNBins_writeToFile")),
      phiMin_writeToFile_(iConfig.getParameter<double>("phiMin_writeToFile")),
      phiMax_writeToFile_(iConfig.getParameter<double>("phiMax_writeToFile")),
      imagesFileName_(iConfig.getParameter<std::string>("imagesFileName")),
      imageCounter_(0) {

  caloGrid_ =
      std::make_unique<TH2F>("caloGrid", "Calorimeter grid", cnnImageSizeEta_, etaLow_, etaUp_, cnnImageSizePhi_, phiLow_, phiUp_);
  caloGrid_->GetXaxis()->SetTitle("#eta");
  caloGrid_->GetYaxis()->SetTitle("#phi");
  if ( writeImages_ ) {
    imageOutputFile_ = new TFile(imagesFileName_.c_str(), "RECREATE");
  caloGridWriteToFile_ =
      std::make_unique<TH2F>("caloGrid", "Calorimeter grid", etaNBins_writeToFile_, etaMin_writeToFile_, etaMax_writeToFile_, phiNBins_writeToFile_, phiMin_writeToFile_, phiMax_writeToFile_);
  caloGridWriteToFile_->GetXaxis()->SetTitle("#eta");
  caloGridWriteToFile_->GetYaxis()->SetTitle("#phi");
  }

  caloGridTest_ =
      std::make_unique<TH2F>("caloGridTest", "Calorimeter grid", etaNBins_writeToFile_, etaMin_writeToFile_, etaMax_writeToFile_, phiNBins_writeToFile_, phiMin_writeToFile_, phiMax_writeToFile_);
  caloGridTest_->GetXaxis()->SetTitle("#eta");
  caloGridTest_->GetYaxis()->SetTitle("#phi");

  produces<double>(outputCollectionName_ + "CNN").setBranchAlias(outputCollectionName_ + "CNN");
  
}
CNNProducer::~CNNProducer() {}

void CNNProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<reco::Candidate>> inputCollectionHandle;
  iEvent.getByToken(inputCollectionTag_, inputCollectionHandle);

  edm::Handle<std::vector<l1t::EtSum>> met;
  iEvent.getByToken(metTag_, met);

  edm::Handle<std::vector<l1t::EtSum>> jetSums;
  iEvent.getByToken(jetSumsTag_, jetSums);

  // std::cout << "MET : " << met->at(0).et() << std::endl;
  // std::cout << "HT : " << jetSums->at(0).et() << std::endl;

  // sort inputs into PF regions
  std::vector<std::vector<reco::CandidatePtr>> inputsInRegions_1 = prepareInputsIntoRegions<>(inputCollectionHandle);

  // histogramming the data
  caloGrid_->Reset();
  for (unsigned int iInputRegion = 0; iInputRegion < inputsInRegions_1.size(); ++iInputRegion) {
    fillCaloGrid<>(*(caloGrid_), inputsInRegions_1[iInputRegion], iInputRegion);
  }

  // std::cout << "Histogram for inference" << std::endl;
  // printHistgoram(caloGrid_);

  if ( writeImages_ ) {
    caloGridWriteToFile_->Reset();
    for (unsigned int iInputRegion = 0; iInputRegion < inputsInRegions_1.size(); ++iInputRegion) {
      fillCaloGrid<>(*(caloGridWriteToFile_), inputsInRegions_1[iInputRegion], iInputRegion);
    }

    caloGridTest_->Reset();
    for (unsigned int i = 0; i < inputCollectionHandle->size(); ++i) {
      reco::CandidatePtr tp(inputCollectionHandle, i);
      if ( abs(tp->eta()) > 3 ) continue;
      caloGridTest_->Fill(tp->eta(), tp->phi(), tp->pt());
    }

    imageOutputFile_->cd();
    // saving histo and closing file
    caloGridTest_->Write(std::to_string(imageCounter_).c_str());
    ++imageCounter_;
  }



  double cnnScore = 0;
  if ( runInference_ ) {
    //load in the network

    tensorflow::GraphDef* graphDef = tensorflow::loadGraphDef(cnnGraph_.fullPath());

    //create a session
    tensorflow::Session* session = tensorflow::createSession(graphDef);

    //create the input tensor
    std::vector<float> energyTowers;
    unsigned int fullImageSizePhi = cnnImageSizePhi_ + 2*cnnNPadPhi_;
    unsigned int fullImageSizeEta = cnnImageSizeEta_ + 2*cnnNPadEta_;
    tensorflow::Tensor input_1(tensorflow::DT_FLOAT, {1,fullImageSizeEta,fullImageSizePhi,1}); // single batch of dimension 16x18x1

    for (uint eta_ind=0; eta_ind < fullImageSizeEta; eta_ind++){
      for (uint phi_ind=0; phi_ind < fullImageSizePhi; phi_ind++) {
        input_1.tensor<float,4>()(0,eta_ind,phi_ind,0) = 0;
      }
    }

    for (uint eta_ind=0; eta_ind < cnnImageSizeEta_; eta_ind++){
        for (uint phi_ind=0; phi_ind < cnnImageSizePhi_; phi_ind++) {
          // float towerEnergy = getTowerEnergy(eta_ind+1, phi_ind+1);
          float towerEnergy = caloGrid_->GetBinContent(eta_ind+1, phi_ind+1);
          if ( towerEnergy > maxPixelValue_ ) towerEnergy = maxPixelValue_;
          float scaledTowerEnergy = towerEnergy / maxPixelValue_;
          // std::cout << "Setting index : " << eta_ind+cnnNPadEta_ << " " << phi_ind+cnnNPadPhi_ << " " << scaledTowerEnergy << std::endl;
          input_1.tensor<float,4>()(0,eta_ind+cnnNPadEta_,phi_ind+cnnNPadPhi_,0) = scaledTowerEnergy;
        }
    }

    // Handle padding
    // Padding in eta + corners
    for (uint eta_ind=0; eta_ind < cnnNPadEta_; eta_ind++){
      for (uint phi_ind=0; phi_ind < fullImageSizePhi; phi_ind++) {
        input_1.tensor<float,4>()(0,eta_ind,phi_ind,0) = 0;
        input_1.tensor<float,4>()(0,fullImageSizeEta-eta_ind-1,phi_ind,0) = 0;
      }
    }

    // Padding in phi
    for (uint eta_ind=cnnNPadEta_; eta_ind < cnnImageSizeEta_+cnnNPadEta_; eta_ind++){
      for (uint phi_ind=0; phi_ind < cnnNPadPhi_; phi_ind++) {
        unsigned int phiIndexPadLow = cnnImageSizePhi_ + phi_ind;
        input_1.tensor<float,4>()(0,eta_ind,phi_ind,0) = input_1.tensor<float,4>()(0,eta_ind,phiIndexPadLow,0);
        unsigned int phiIndexPadHigh = cnnNPadPhi_ + phi_ind;
        input_1.tensor<float,4>()(0,eta_ind,fullImageSizePhi-cnnNPadPhi_+phi_ind,0) = input_1.tensor<float,4>()(0,eta_ind,phiIndexPadHigh,0);
      }
    }

    std::vector<tensorflow::Tensor> outputs;
    tensorflow::run(session, { { "input_1", input_1 } }, { "Identity" }, &outputs);
    
    cnnScore = outputs[0].flat<float>()(0);
    // std::cout << cnnScore << std::endl;
    // if ( cnnScore != cnnScore ) {
    //   std::cout << "Nan score : " << cnnScore << std::endl;
    //   printHistogram(caloGrid_);
    // }

    //clean-up
    tensorflow::closeSession(session);
    delete graphDef;
  }

  std::unique_ptr<double> cnnScorePtr = std::make_unique<double>(cnnScore);
  iEvent.put(std::move(cnnScorePtr), this->outputCollectionName_ + "CNN");

}

template <class Container>
void CNNProducer::fillCaloGrid(TH2F& caloGrid,
                                        const Container& triggerPrimitives,
                                        const unsigned int regionIndex) {
  //Filling the calo grid with the primitives
  for (const auto& primitiveIterator : triggerPrimitives) {
    // Get digitised (floating point with reduced precision) eta and phi
    std::pair<float, float> digi_EtaPhi =
        getCandidateDigiEtaPhi(primitiveIterator->eta(), primitiveIterator->phi(), regionIndex);

    caloGrid.Fill(static_cast<float>(digi_EtaPhi.first),
                  static_cast<float>(digi_EtaPhi.second),
                  static_cast<float>(primitiveIterator->pt()));
  }
}

std::pair<float, float> CNNProducer::getCandidateDigiEtaPhi(const float eta,
                                                                     const float phi,
                                                                     const unsigned int regionIndex) const {
  std::pair<double, double> regionLowEdges = regionEtaPhiLowEdges(regionIndex);

  int digitisedEta = floor((eta - regionLowEdges.second) / etalsb_);
  int digitisedPhi = floor((phi - regionLowEdges.first) / philsb_);

  // If eta or phi is on a bin edge
  // Put in bin above, to match behaviour of HLS
  // Unless it's on the last bin of this pf region
  // Then it is placed in the last bin, not the overflow
  TAxis* etaAxis = caloGrid_->GetXaxis();
  std::pair<double, double> regionUpEdges = regionEtaPhiUpEdges(regionIndex);
  int digiEtaEdgeLastBinUp = floor((regionUpEdges.second - regionLowEdges.second) / etalsb_);
  // If the digi eta is outside the last bin of this pf region
  // Set the digitised quantity so it would be in the last bin
  // These cases could be avoided by sorting input candidates based on digitised eta/phi
  if (digitisedEta >= digiEtaEdgeLastBinUp) {
    digitisedEta = digiEtaEdgeLastBinUp - 1;
  } else {
    for (int i = 0; i < etaAxis->GetNbins(); ++i) {
      if (etaAxis->GetBinUpEdge(i) < regionLowEdges.second)
        continue;
      int digiEdgeBinUp = floor((etaAxis->GetBinUpEdge(i) - regionLowEdges.second) / etalsb_);
      if (digiEdgeBinUp == digitisedEta) {
        digitisedEta += 1;
      }
    }
  }

  // Similar for phi
  TAxis* phiAxis = caloGrid_->GetYaxis();
  int digiPhiEdgeLastBinUp = floor((regionUpEdges.first - regionLowEdges.first) / philsb_);
  if (digitisedPhi >= digiPhiEdgeLastBinUp) {
    digitisedPhi = digiPhiEdgeLastBinUp - 1;
  } else {
    for (int i = 0; i < phiAxis->GetNbins(); ++i) {
      if (phiAxis->GetBinUpEdge(i) < regionLowEdges.first)
        continue;
      int digiEdgeBinUp = floor((phiAxis->GetBinUpEdge(i) - regionLowEdges.first) / philsb_);
      if (digiEdgeBinUp == digitisedPhi) {
        digitisedPhi += 1;
      }
    }
  }

  // Convert digitised eta and phi back to floating point quantities with reduced precision
  float floatDigitisedEta = (digitisedEta + 0.5) * etalsb_ + regionLowEdges.second;
  float floatDigitisedPhi = (digitisedPhi + 0.5) * philsb_ + regionLowEdges.first;

  return std::pair<float, float>{floatDigitisedEta, floatDigitisedPhi};
}

void CNNProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("inputCollectionTag", edm::InputTag("l1pfCandidates", "Puppi"));
  desc.add<edm::InputTag>("metTag", edm::InputTag("Phase1L1TJetProducer9x9trimmed" ,   "UncalibratedPhase1L1TJetFromPfCandidatesMET"));
  desc.add<edm::InputTag>("jetSumsTag", edm::InputTag("Phase1L1TJetSumsProducer9x9trimmed" ,   "Sums"));
  desc.add<double>("phiLow", -M_PI);
  desc.add<double>("phiUp", M_PI);
  desc.add<double>("etaLow", -3);
  desc.add<double>("etaUp", 3);
  desc.add<double>("ptlsb", 0.25), desc.add<double>("philsb", 0.0043633231), desc.add<double>("etalsb", 0.0043633231),
  desc.add<string>("outputCollectionName", "UncalibratedPhase1L1TJetFromPfCandidates");
  desc.add<std::vector<double>>("etaRegions");
  desc.add<std::vector<double>>("phiRegions");
  desc.add<unsigned int>("maxInputsPerRegion", 18);
  desc.add<bool>("runInference", false);
  desc.add<edm::FileInPath>("cnnGraph",edm::FileInPath("L1Trigger/L1CaloTrigger/data/Qcnn_graph_bigNetwork.pb"));
  desc.add<unsigned int>("cnnImageSizeEta", 12);
  desc.add<unsigned int>("cnnImageSizePhi", 12);
  desc.add<unsigned int>("cnnNPadEta", 3);
  desc.add<unsigned int>("cnnNPadPhi", 3);
  desc.add<double>("maxPixelValue", 63);

  desc.add<bool>("writeImages", false);
  desc.add<unsigned int >("etaNBins_writeToFile", 120 );
  desc.add<double>("etaMin_writeToFile", -5);
  desc.add<double>("etaMax_writeToFile", 5);
  desc.add<unsigned int >("phiNBins_writeToFile", 72);
  desc.add<double>("phiMin_writeToFile", -M_PI);
  desc.add<double>("phiMax_writeToFile", M_PI);
  desc.add<string>("imagesFileName", "images.root");
  descriptions.add("CNNProducer", desc);
}

template <class Handle>
std::vector<std::vector<edm::Ptr<reco::Candidate>>> CNNProducer::prepareInputsIntoRegions(
    const Handle triggerPrimitives) {
  std::vector<std::vector<reco::CandidatePtr>> inputsInRegions{etaRegionEdges_.size() * (phiRegionEdges_.size() - 1)};

  for (unsigned int i = 0; i < triggerPrimitives->size(); ++i) {
    reco::CandidatePtr tp(triggerPrimitives, i);

    if (tp->phi() < phiRegionEdges_.front() || tp->phi() >= phiRegionEdges_.back() ||
        tp->eta() < etaRegionEdges_.front() || tp->eta() >= etaRegionEdges_.back())
      continue;

    // Which phi region does this tp belong to
    auto it_phi = phiRegionEdges_.begin();
    it_phi = std::upper_bound(phiRegionEdges_.begin(), phiRegionEdges_.end(), tp->phi()) - 1;

    // Which eta region does this tp belong to
    auto it_eta = etaRegionEdges_.begin();
    it_eta = std::upper_bound(etaRegionEdges_.begin(), etaRegionEdges_.end(), tp->eta()) - 1;

    if (it_phi != phiRegionEdges_.end() && it_eta != etaRegionEdges_.end()) {
      auto phiRegion = it_phi - phiRegionEdges_.begin();
      auto etaRegion = it_eta - etaRegionEdges_.begin();
      inputsInRegions[getRegionIndex(phiRegion, etaRegion)].emplace_back(tp);
    }
  }

  // Truncate number of inputs in each pf region
  for (auto& inputs : inputsInRegions) {

    

    if (inputs.size() > maxInputsPerRegion_) {

      std::sort(inputs.begin(), inputs.end(), [](const reco::CandidatePtr& input1, const reco::CandidatePtr& input2) {
        return input1->pt() > input2->pt();
      });

      inputs.resize(maxInputsPerRegion_);
    }
  }

  return inputsInRegions;
}

unsigned int CNNProducer::getRegionIndex(const unsigned int phiRegion, const unsigned int etaRegion) const {
  return etaRegion * (phiRegionEdges_.size() - 1) + phiRegion;
}

std::pair<double, double> CNNProducer::regionEtaPhiLowEdges(const unsigned int regionIndex) const {
  unsigned int phiRegion = regionIndex % (phiRegionEdges_.size() - 1);
  unsigned int etaRegion = (regionIndex - phiRegion) / (phiRegionEdges_.size() - 1);
  return std::pair<double, double>{phiRegionEdges_.at(phiRegion), etaRegionEdges_.at(etaRegion)};
}

std::pair<double, double> CNNProducer::regionEtaPhiUpEdges(const unsigned int regionIndex) const {
  unsigned int phiRegion = regionIndex % (phiRegionEdges_.size() - 1);
  unsigned int etaRegion = (regionIndex - phiRegion) / (phiRegionEdges_.size() - 1);
  if (phiRegion == phiRegionEdges_.size() - 1) {
    return std::pair<double, double>{phiRegionEdges_.at(phiRegion), etaRegionEdges_.at(etaRegion + 1)};
  } else if (etaRegion == etaRegionEdges_.size() - 1) {
    return std::pair<double, double>{phiRegionEdges_.at(phiRegion + 1), etaRegionEdges_.at(etaRegion)};
  }

  return std::pair<double, double>{phiRegionEdges_.at(phiRegion + 1), etaRegionEdges_.at(etaRegion + 1)};
}

  void CNNProducer::printHistogram(std::unique_ptr<TH2F>& histogram) const {
    unsigned int nBinsEta = histogram->GetNbinsX();
    unsigned int nBinsPhi = histogram->GetNbinsY();
    for (unsigned int iPhi = 1; iPhi < nBinsPhi+1; ++iPhi) {
      for (unsigned int iEta = 1; iEta < nBinsEta+1; ++iEta) {
        std::cout << histogram->GetBinContent(iEta,iPhi) << " ";
      }
      std::cout << std::endl;
    }
  }

    
DEFINE_FWK_MODULE(CNNProducer);