// -*- C++ -*-
//
// Package: L1CaloTrigger
// Class: Phase1L1TJetSeedProducer
//
/**\class Phase1L1TJetSeedProducer Phase1L1TJetSeedProducer.cc L1Trigger/L1CaloTrigger/plugin/Phase1L1TJetSeedProducer.cc
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
#include "L1Trigger/Phase2L1ParticleFlow/src/newfirmware/dataformats/puppi.h"
#include "L1Trigger/Phase2L1ParticleFlow/src/newfirmware/common/bitonic_hybrid_sort_ref.h"

#include "TH2F.h"

#include <cmath>

#include <algorithm>

class Phase1L1TJetSeedProducer : public edm::one::EDProducer<> {
public:
  explicit Phase1L1TJetSeedProducer(const edm::ParameterSet&);
  ~Phase1L1TJetSeedProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  /// Finds the seeds in the caloGrid, seeds are saved in a vector that contain the index in the TH2F of each seed
  l1t::PFCandidateCollection findSeeds(float seedThreshold) const;

  /// Get the energy of a certain tower while correctly handling phi periodicity in case of overflow
  float getTowerEnergy(int iEta, int iPhi) const;

  // Sort the seeds as in firmware
  void sortSeeds(const l1t::PFCandidateCollection unsortedJets, l1t::PFCandidateCollection& sortedJets );
  template <typename T>
  void hybridBitonicMergeRef(std::vector<T>& a, int N, int low, bool dir);
  template <typename T>
  void hybridBitonicSortRef(std::vector<T>& a, int N, int low, bool dir);
  template <typename T>
  void hybrid_bitonic_sort_and_crop_ref(unsigned int nIn, unsigned int nOut, std::vector<T> in, std::vector<T>& out);
  template <typename T>
  void compAndSwap(std::vector<T>& a, unsigned int i, unsigned int j, bool dir=false);
  template <typename T>
  void swap(T& a, T& b);

  // <3 handy method to fill the calogrid with whatever type
  template <class Container>
  void fillCaloGrid(TH2F& caloGrid, const Container& triggerPrimitives, const unsigned int regionIndex);

  // Digitise the eta and phi coordinates of input candidates
  // This converts the quantities to integers to reduce precision
  // And takes account of bin edge effects i.e. makes sure the
  // candidate ends up in the correct (i.e. same behaviour as the firmware) bin of caloGrid_
  std::pair<unsigned, unsigned> getCandidateBin(const float eta,
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
  std::pair<unsigned, unsigned> regionEtaPhiBinOffset(const unsigned int regionIndex) const;

  // Determine if this tower should be trimmed or not
  // Used only when trimmedGrid_ option is set to true
  // Trim means removing 3 towers in each corner of the square grid
  // giving a cross shaped grid, which is a bit more circular in shape than a square
  bool trimTower(const int etaIndex, const int phiIndex) const;

  edm::EDGetTokenT<edm::View<reco::Candidate>> inputCollectionTag_;
  // histogram containing our clustered inputs
  std::unique_ptr<TH2F> caloGrid_;

  std::vector<double> etaBinning_;
  size_t nBinsEta_;
  unsigned int nBinsPhi_;
  double phiLow_;
  double phiUp_;
  unsigned int jetIEtaSize_;
  unsigned int jetIPhiSize_;
  bool trimmedGrid_;
  double seedPtThreshold_;
  double ptlsb_;
  double philsb_;
  double etalsb_;
  // Eta and phi edges of input PF regions
  std::vector<double> etaRegionEdges_;
  std::vector<double> phiRegionEdges_;
  // Maximum number of candidates per input PF region
  unsigned int maxInputsPerRegion_;
  std::string outputCollectionName_;

};

Phase1L1TJetSeedProducer::Phase1L1TJetSeedProducer(const edm::ParameterSet& iConfig)
    :  // getting configuration settings
      inputCollectionTag_{
          consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("inputCollectionTag"))},
      etaBinning_(iConfig.getParameter<std::vector<double>>("etaBinning")),
      nBinsEta_(etaBinning_.size() - 1),
      nBinsPhi_(iConfig.getParameter<unsigned int>("nBinsPhi")),
      phiLow_(iConfig.getParameter<double>("phiLow")),
      phiUp_(iConfig.getParameter<double>("phiUp")),
      jetIEtaSize_(iConfig.getParameter<unsigned int>("jetIEtaSize")),
      jetIPhiSize_(iConfig.getParameter<unsigned int>("jetIPhiSize")),
      trimmedGrid_(iConfig.getParameter<bool>("trimmedGrid")),
      seedPtThreshold_(iConfig.getParameter<double>("seedPtThreshold")),
      ptlsb_(iConfig.getParameter<double>("ptlsb")),
      philsb_(iConfig.getParameter<double>("philsb")),
      etalsb_(iConfig.getParameter<double>("etalsb")),
      etaRegionEdges_(iConfig.getParameter<std::vector<double>>("etaRegions")),
      phiRegionEdges_(iConfig.getParameter<std::vector<double>>("phiRegions")),
      maxInputsPerRegion_(iConfig.getParameter<unsigned int>("maxInputsPerRegion")),
      outputCollectionName_(iConfig.getParameter<std::string>("outputCollectionName")) {
  caloGrid_ =
      std::make_unique<TH2F>("caloGrid", "Calorimeter grid", nBinsEta_, etaBinning_.data(), nBinsPhi_, phiLow_, phiUp_);
  caloGrid_->GetXaxis()->SetTitle("#eta");
  caloGrid_->GetYaxis()->SetTitle("#phi");
  // produces<l1t::PFCandidateCollection>(outputCollectionName_).setBranchAlias(outputCollectionName_);
  produces<l1t::PFCandidateCollection>(outputCollectionName_);

}

Phase1L1TJetSeedProducer::~Phase1L1TJetSeedProducer() {}

float Phase1L1TJetSeedProducer::getTowerEnergy(int iEta, int iPhi) const {
  // We return the pt of a certain bin in the calo grid, taking account of the phi periodicity when overflowing (e.g. phi > phiSize), and returning 0 for the eta out of bounds

  int nBinsEta = caloGrid_->GetNbinsX();
  int nBinsPhi = caloGrid_->GetNbinsY();
  while (iPhi < 1) {
    iPhi += nBinsPhi;
  }
  while (iPhi > nBinsPhi) {
    iPhi -= nBinsPhi;
  }
  if (iEta < 1) {
    return 0;
  }
  if (iEta > nBinsEta) {
    return 0;
  }
  return caloGrid_->GetBinContent(iEta, iPhi);
}

void Phase1L1TJetSeedProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<reco::Candidate>> inputCollectionHandle;
  iEvent.getByToken(inputCollectionTag_, inputCollectionHandle);

  // sort inputs into PF regions
  std::vector<std::vector<reco::CandidatePtr>> inputsInRegions = prepareInputsIntoRegions<>(inputCollectionHandle);

  // histogramming the data
  caloGrid_->Reset();
  for (unsigned int iInputRegion = 0; iInputRegion < inputsInRegions.size(); ++iInputRegion) {
    fillCaloGrid<>(*(caloGrid_), inputsInRegions[iInputRegion], iInputRegion);
  }

  // find the seeds
  const auto& seedsVector = findSeeds(seedPtThreshold_);  // seedPtThreshold = 5

  // sort by pt
  l1t::PFCandidateCollection sortedSeeds;
  sortSeeds( seedsVector, sortedSeeds );

  auto seedsVectorPtr = std::make_unique<l1t::PFCandidateCollection>(sortedSeeds);
  iEvent.put(std::move(seedsVectorPtr), outputCollectionName_ );

  return;
}

l1t::PFCandidateCollection Phase1L1TJetSeedProducer::findSeeds(float seedThreshold) const {
  int nBinsX = caloGrid_->GetNbinsX();
  int nBinsY = caloGrid_->GetNbinsY();

  l1t::PFCandidateCollection seeds;

  int etaHalfSize = (int)jetIEtaSize_ / 2;
  int phiHalfSize = (int)jetIPhiSize_ / 2;

  // for each point of the grid check if it is a local maximum
  // to do so I take a point, and look if is greater than the points around it (in the 9x9 neighborhood)
  // to prevent mutual exclusion, I check greater or equal for points above and right to the one I am considering (including the top-left point)
  // to prevent mutual exclusion, I check greater for points below and left to the one I am considering (including the bottom-right point)

  for (int iPhi = 1; iPhi <= nBinsY; iPhi++) {
    for (int iEta = 1; iEta <= nBinsX; iEta++) {
      float centralPt = caloGrid_->GetBinContent(iEta, iPhi);
      if (centralPt < seedThreshold)
        continue;
      bool isLocalMaximum = true;
      // Scanning through the grid centered on the seed
      for (int etaIndex = -etaHalfSize; etaIndex <= etaHalfSize; etaIndex++) {
        for (int phiIndex = -phiHalfSize; phiIndex <= phiHalfSize; phiIndex++) {
          if (trimmedGrid_) {
            if (trimTower(etaIndex, phiIndex))
              continue;
          }

          if ((etaIndex == 0) && (phiIndex == 0))
            continue;
          if (etaIndex > 0) {
            isLocalMaximum = ((isLocalMaximum) && (centralPt > getTowerEnergy(iEta + etaIndex, iPhi + phiIndex)));
          } else if ( etaIndex < 0 ) {
            isLocalMaximum = ((isLocalMaximum) && (centralPt >= getTowerEnergy(iEta + etaIndex, iPhi + phiIndex)));
          }
          else {
            if ( phiIndex > 0 ) {
              isLocalMaximum = ((isLocalMaximum) && (centralPt > getTowerEnergy(iEta + etaIndex, iPhi + phiIndex)));
            }
            else {
              isLocalMaximum = ((isLocalMaximum) && (centralPt >= getTowerEnergy(iEta + etaIndex, iPhi + phiIndex)));
            }
          }
        }
      }

      if (isLocalMaximum) {
        l1t::PFCandidate p;
        reco::Candidate::PolarLorentzVector pfVector;

        const float etaLSB = 1.5 / 18;
        double etaBinCentre = -3 + (iEta-1+0.5)*etaLSB;

        const float phiLSB = 2. * M_PI / 72;
        double phiBinCentre = -3.15 + ( iPhi-1+0.5 ) * phiLSB;

        pfVector.SetPt(centralPt);
        pfVector.SetPhi(phiBinCentre);
        pfVector.SetEta(etaBinCentre);
        p.setP4( pfVector );

        seeds.emplace_back(p);
      }
    }
  }

  return seeds;
}

void Phase1L1TJetSeedProducer::sortSeeds(const l1t::PFCandidateCollection unsortedSeeds, l1t::PFCandidateCollection& sortedSeeds ) {

  const unsigned int nEtaRegions = 4;
  const unsigned int nInputsPerSortModule = 18;
  const unsigned int nOutputSeedsPerEtaRegion = 4;
  const unsigned int nOutputSeedsToGT = 12;

  unsigned int nUnsortedSeeds = unsortedSeeds.size();
  // Get seeds into the regions and time ordering seen in firmware
  std::vector< std::vector< std::vector< l1t::PFCandidateCollection > > > seedsPerEtaPhiRegions( 
    nEtaRegions, std::vector< std::vector< l1t::PFCandidateCollection > > (
      2, std::vector< l1t::PFCandidateCollection > (
        nInputsPerSortModule, l1t::PFCandidateCollection() ) ) );


  for ( const auto& seed : unsortedSeeds ) { 
    unsigned int etaRegion = (seed.eta()+3)/1.5;
    unsigned int seedEtaBin = floor( ( seed.eta() + (2 - 1.0*etaRegion) * 1.5 ) / 0.0833 );
    unsigned int seedPhiBin = floor( ( seed.phi() + M_PI ) / 0.0875 );
    unsigned int phiRegion = ( ( seedPhiBin ) % 4 ) / 2;
    seedsPerEtaPhiRegions[etaRegion][phiRegion][seedPhiBin/4].push_back(seed);
  }

  // Rotate to first phi region found in firmware
  for ( unsigned iEtaRegion = 0; iEtaRegion < nEtaRegions; ++iEtaRegion ) {
    for ( unsigned iPhiRegion = 0; iPhiRegion < 2; ++ iPhiRegion ) {
      std::rotate( seedsPerEtaPhiRegions[iEtaRegion][iPhiRegion].begin(), seedsPerEtaPhiRegions[iEtaRegion][iPhiRegion].begin()+8, seedsPerEtaPhiRegions[iEtaRegion][iPhiRegion].end() );
    }
  }

  // Push seeds in first phi bin to back, as these are found last after receiving all bins (i.e. handling of phi wrap-around)
  for ( unsigned iEtaRegion = 0; iEtaRegion < nEtaRegions; ++iEtaRegion ) {
    for ( unsigned iPhiRegion = 0; iPhiRegion < 2; ++ iPhiRegion ) {
      std::rotate( seedsPerEtaPhiRegions[iEtaRegion][iPhiRegion].begin(), seedsPerEtaPhiRegions[iEtaRegion][iPhiRegion].begin()+1, seedsPerEtaPhiRegions[iEtaRegion][iPhiRegion].end() );
    }
  }
  std::vector< l1t::PFCandidate > sortedSeedsAllEta;
  for ( unsigned iEtaRegion = 0; iEtaRegion < nEtaRegions; ++iEtaRegion ) {
    std::vector<l1t::PFCandidate > sortedSeedsInEtaRegion;
    for ( unsigned iPhiRegion = 0; iPhiRegion < 2; ++ iPhiRegion ) {
      std::vector<l1t::PFCandidate > sortedSeeds( 4, l1t::PFCandidate() );
      for ( unsigned int iInputClock = 0; iInputClock < nInputsPerSortModule; ++iInputClock ) {

        // Sort input seeds
        l1t::PFCandidateCollection inputSeeds = seedsPerEtaPhiRegions[iEtaRegion][iPhiRegion][iInputClock];
        // First by eta
        std::sort(inputSeeds.begin(), inputSeeds.end(), [](l1t::PFCandidate seed1, l1t::PFCandidate seed2) {
          return seed1.eta() < seed2.eta();
        });
        inputSeeds.resize(nOutputSeedsPerEtaRegion);
        hybrid_bitonic_sort_and_crop_ref(4,4,inputSeeds,inputSeeds);

        // Add to list of top 4 seeds so far
        // Merge with top 4 seeds so far, and sort
        sortedSeeds.insert( sortedSeeds.end(), inputSeeds.begin(), inputSeeds.end() );
        std::reverse(sortedSeeds.begin(),sortedSeeds.begin()+nOutputSeedsPerEtaRegion);
        for (int i = 0; i < 4; i++) {
            compAndSwap(sortedSeeds, i, i + 4, 0);
        }

        sortedSeeds.resize(nOutputSeedsPerEtaRegion);
        std::reverse(sortedSeeds.begin(),sortedSeeds.end());
        compAndSwap(sortedSeeds, 0, 2); 
        compAndSwap(sortedSeeds, 1, 3); 
        //---
        compAndSwap(sortedSeeds, 0, 1); 
        compAndSwap(sortedSeeds, 2, 3); 
      }

      if ( iPhiRegion % 2 == 0 ) {
        sortedSeedsInEtaRegion.insert( sortedSeedsInEtaRegion.end(), sortedSeeds.rbegin(), sortedSeeds.rend() );
      }
      else {
        sortedSeedsInEtaRegion.insert( sortedSeedsInEtaRegion.end(), sortedSeeds.begin(), sortedSeeds.end() );
      }
    }
    // Sort 8 seeds in each eta region
    // std::cout << "8 seeds in one of the regions, before merge" << std::endl;
    std::reverse(sortedSeedsInEtaRegion.begin(),sortedSeedsInEtaRegion.end());
    hybridBitonicMergeRef(sortedSeedsInEtaRegion,nOutputSeedsPerEtaRegion*2,0,false);

    if ( iEtaRegion % 2 == 0 ) {
      sortedSeedsAllEta.insert(sortedSeedsAllEta.end(), sortedSeedsInEtaRegion.rbegin(), sortedSeedsInEtaRegion.rend() );
    }
    else {
      sortedSeedsAllEta.insert(sortedSeedsAllEta.end(), sortedSeedsInEtaRegion.begin(), sortedSeedsInEtaRegion.end() );
    }
  }
  hybridBitonicMergeRef(sortedSeedsAllEta,nOutputSeedsPerEtaRegion*2*2,0,false);
  hybridBitonicMergeRef(sortedSeedsAllEta,nOutputSeedsPerEtaRegion*2*2,nOutputSeedsPerEtaRegion*2*2,false);
  std::reverse(sortedSeedsAllEta.begin(),sortedSeedsAllEta.begin()+nOutputSeedsPerEtaRegion*2*2);

  for ( unsigned int iJet = 0; iJet < nOutputSeedsPerEtaRegion*2*2 - nOutputSeedsToGT; ++iJet ) {
    sortedSeedsAllEta.erase(sortedSeedsAllEta.begin());
    sortedSeedsAllEta.erase(sortedSeedsAllEta.end()-1);
  }

  hybridBitonicMergeRef(sortedSeedsAllEta,nOutputSeedsToGT*2,0,false);
  sortedSeedsAllEta.resize(nOutputSeedsToGT);
  unsigned int nSeedsGT0=0;
  for ( const auto& iJet : sortedSeedsAllEta ) {
    if ( iJet.pt() > 0 ) {
      sortedSeeds.push_back( iJet );
      ++nSeedsGT0;
    }
  }
}

template <typename T>
void Phase1L1TJetSeedProducer::hybridBitonicMergeRef(std::vector<T>& a, int N, int low, bool dir) {
  int k = hybridBitonicSortUtils::PowerOf2LessThan(N);
  int k2 = N - k;
  if (N > 1) {
    for (int i = low; i < low + k; i++) {
      if (i + k < low + N)
        compAndSwap(a, i, i + k, dir);
    }
    if (N > 2) {
      hybridBitonicMergeRef(a, k, low, dir);
      hybridBitonicMergeRef(a, k2, low + k, dir);
    }
  }
}

template <typename T>
void Phase1L1TJetSeedProducer::hybridBitonicSortRef(std::vector<T>& a, int N, int low, bool dir) {
  // general case
  if (N > 1) {
    int lowerSize = N / 2;
    int upperSize = N - N / 2;
    bool notDir = not dir;
    hybridBitonicSortRef(a, lowerSize, low, notDir);
    hybridBitonicSortRef(a, upperSize, low + lowerSize, dir);
    hybridBitonicMergeRef(a, N, low, dir);
  }
}

template <typename T>
void Phase1L1TJetSeedProducer::hybrid_bitonic_sort_and_crop_ref(
    unsigned int nIn, unsigned int nOut, std::vector<T> in, std::vector<T>& out) {  // just an interface
  std::vector<T> work(nIn, T());
  for (unsigned int i = 0; i < nIn; ++i) {
    work[i] = in[i];
  }
  hybridBitonicSortRef(work, nIn, 0, false);
  for (unsigned int i = 0; i < nOut; ++i) {
    out[i] = work[i];
  }
}

template <typename T>
void Phase1L1TJetSeedProducer::compAndSwap(std::vector<T>& a, unsigned int i, unsigned int j, bool dir) {
  if ( i >= a.size() || j >= a.size() || i == j ) return;

  if (dir) {
    if (a[j].pt()<a[i].pt()) {
      std::swap(a[i],a[j]);
    }
  }
  else {
    if (a[i].pt()<a[j].pt()) {
      std::swap(a[i],a[j]);
    }
  }
}

template <typename T>
void Phase1L1TJetSeedProducer::swap(T& a, T& b) {
  T c = a;
  a=b;
  b=c;
  return;
}

template <class Container>
void Phase1L1TJetSeedProducer::fillCaloGrid(TH2F& caloGrid,
                                        const Container& triggerPrimitives,
                                        const unsigned int regionIndex) {
  //Filling the calo grid with the primitives
  for (const auto& primitiveIterator : triggerPrimitives) {
    // if ( primitiveIterator->pt() != 12.25 ) continue;
    // Get digitised (floating point with reduced precision) eta and phi
    std::pair<unsigned, unsigned> bin_EtaPhi =
        getCandidateBin(primitiveIterator->eta(), primitiveIterator->phi(), regionIndex);
    unsigned int globalBin = caloGrid.GetBin( bin_EtaPhi.second, bin_EtaPhi.first );
    caloGrid.AddBinContent(globalBin,
                 float( l1ct::pt_t(primitiveIterator->pt()) ) );


  }
}

std::pair<unsigned, unsigned> Phase1L1TJetSeedProducer::getCandidateBin(const float eta,
                                                                     const float phi,
                                                                     const unsigned int regionIndex) const {

  l1ct::glbeta_t glbEta = l1ct::Scales::makeGlbEta( eta );
  l1ct::glbphi_t glbPhi = l1ct::Scales::makeGlbPhi( phi );

  std::pair<double, double> regionLowEdges = regionEtaPhiLowEdges(regionIndex);
  l1ct::glbeta_t etaOffset = l1ct::Scales::makeGlbEta( regionLowEdges.second );
  l1ct::glbphi_t phiOffset = l1ct::Scales::makeGlbPhi( regionLowEdges.first );

  int etaBin = ( glbEta - etaOffset ) / 19 + 1;
  int phiBin = ( glbPhi - phiOffset ) / 20 + 1;

  if ( regionLowEdges.second == -2.5 || regionLowEdges.second == 1.5 ) {
    if ( etaBin >= 12 ) etaBin = 12;
  }
  else if ( etaBin >= 6 ) etaBin = 6;
  if ( phiBin >= 8 ) phiBin = 8;

  std::pair<unsigned, unsigned> binOffsets = regionEtaPhiBinOffset(regionIndex);

  return std::pair<unsigned, unsigned>{phiBin + binOffsets.first, etaBin + binOffsets.second };
}

void Phase1L1TJetSeedProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("inputCollectionTag", edm::InputTag("l1pfCandidates", "Puppi"));
  desc.add<std::vector<double>>("etaBinning");
  desc.add<unsigned int>("nBinsPhi", 72);
  desc.add<double>("phiLow", -M_PI);
  desc.add<double>("phiUp", M_PI);
  desc.add<unsigned int>("jetIEtaSize", 7);
  desc.add<unsigned int>("jetIPhiSize", 7);
  desc.add<bool>("trimmedGrid", false);
  desc.add<double>("seedPtThreshold", 5);
  desc.add<double>("ptlsb", 0.25), desc.add<double>("philsb", 0.0043633231), desc.add<double>("etalsb", 0.0043633231),
  desc.add<string>("outputCollectionName", "UncalibratedPhase1L1TJetFromPfCandidates");
  desc.add<std::vector<double>>("etaRegions");
  desc.add<std::vector<double>>("phiRegions");
  desc.add<unsigned int>("maxInputsPerRegion", 18);
  descriptions.add("Phase1L1TJetSeedProducer", desc);
}

template <class Handle>
std::vector<std::vector<edm::Ptr<reco::Candidate>>> Phase1L1TJetSeedProducer::prepareInputsIntoRegions(
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
    if ( l1ct::Scales::makeGlbPhi( *(it_phi+1) ) == l1ct::Scales::makeGlbPhi( tp->phi() ) ) {
      it_phi += 1;
    }

    // Which eta region does this tp belong to
    auto it_eta = etaRegionEdges_.begin();
    it_eta = std::upper_bound(etaRegionEdges_.begin(), etaRegionEdges_.end(), tp->eta()) - 1;
    if ( l1ct::Scales::makeGlbEta( *(it_eta+1) ) == l1ct::Scales::makeGlbEta( tp->eta() ) ) {
      it_eta += 1;
    }


    if (it_phi != phiRegionEdges_.end() && it_eta != etaRegionEdges_.end()) {
      auto phiRegion = it_phi - phiRegionEdges_.begin();
      auto etaRegion = it_eta - etaRegionEdges_.begin();
      inputsInRegions[getRegionIndex(phiRegion, etaRegion)].emplace_back(tp);
    }
  }

  // Truncate number of inputs in each pf region
  for (auto& inputs : inputsInRegions) {
    if (inputs.size() > maxInputsPerRegion_) {
      inputs.resize(maxInputsPerRegion_);
    }
  }

  return inputsInRegions;
}

unsigned int Phase1L1TJetSeedProducer::getRegionIndex(const unsigned int phiRegion, const unsigned int etaRegion) const {
  return etaRegion * (phiRegionEdges_.size() - 1) + phiRegion;
}

std::pair<double, double> Phase1L1TJetSeedProducer::regionEtaPhiLowEdges(const unsigned int regionIndex) const {
  unsigned int phiRegion = regionIndex % (phiRegionEdges_.size() - 1);
  unsigned int etaRegion = (regionIndex - phiRegion) / (phiRegionEdges_.size() - 1);
  return std::pair<double, double>{phiRegionEdges_.at(phiRegion), etaRegionEdges_.at(etaRegion)};
}


std::pair<unsigned, unsigned> Phase1L1TJetSeedProducer::regionEtaPhiBinOffset(const unsigned int regionIndex) const {
  unsigned int phiRegion = regionIndex % (phiRegionEdges_.size() - 1);
  unsigned int etaRegion = (regionIndex - phiRegion) / (phiRegionEdges_.size() - 1);

  float etaBinOffset = ( 3 + etaRegionEdges_.at(etaRegion) ) / 0.5 * 6;
  float phiBinOffset = ( 3.15 + phiRegionEdges_.at(phiRegion) ) / 0.7 * 8;
  return std::pair<unsigned, unsigned>{phiBinOffset, etaBinOffset};
}

std::pair<double, double> Phase1L1TJetSeedProducer::regionEtaPhiUpEdges(const unsigned int regionIndex) const {
  unsigned int phiRegion = regionIndex % (phiRegionEdges_.size() - 1);
  unsigned int etaRegion = (regionIndex - phiRegion) / (phiRegionEdges_.size() - 1);
  if (phiRegion == phiRegionEdges_.size() - 1) {
    return std::pair<double, double>{phiRegionEdges_.at(phiRegion), etaRegionEdges_.at(etaRegion + 1)};
  } else if (etaRegion == etaRegionEdges_.size() - 1) {
    return std::pair<double, double>{phiRegionEdges_.at(phiRegion + 1), etaRegionEdges_.at(etaRegion)};
  }

  return std::pair<double, double>{phiRegionEdges_.at(phiRegion + 1), etaRegionEdges_.at(etaRegion + 1)};
}

bool Phase1L1TJetSeedProducer::trimTower(const int etaIndex, const int phiIndex) const {
  int etaHalfSize = jetIEtaSize_ / 2;
  int phiHalfSize = jetIPhiSize_ / 2;

  if (etaIndex == -etaHalfSize || etaIndex == etaHalfSize) {
    if (phiIndex <= -phiHalfSize + 1 || phiIndex >= phiHalfSize - 1) {
      return true;
    }
  } else if (etaIndex == -etaHalfSize + 1 || etaIndex == etaHalfSize - 1) {
    if (phiIndex == -phiHalfSize || phiIndex == phiHalfSize) {
      return true;
    }
  }

  return false;
}
DEFINE_FWK_MODULE(Phase1L1TJetSeedProducer);
