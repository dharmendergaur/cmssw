// -*- C++ -*-
//
// Package:    L1Trigger/L1THGCalUtilities
// Class:      Stage2FileReader
//
/**\class Stage2FileReader Stage2FileReader.cc L1Trigger/L1THGCalUtilities/plugins/patternFiles/Stage2FileReader.cc

 Description: EDAnalyzer for reading I/O buffer files for hardware/firmware tests of HGC stage 2

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Emyr Clement
//         Created:  Tue, 13 Apr 2022
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "L1Trigger/DemonstratorTools/interface/BoardDataReader.h"
#include "L1Trigger/DemonstratorTools/interface/utilities.h"
#include "L1Trigger/L1THGCalUtilities/interface/patternFiles/codecs_clusters.h"

//
// class declaration
//

class Stage2FileReader : public edm::stream::EDProducer<> {
public:
  explicit Stage2FileReader(const edm::ParameterSet&);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // ----------constants, enums and typedefs ---------
  static constexpr size_t kFramesPerTMUXPeriod = 9;
  static constexpr size_t kGapLengthOutput = 0;
  static constexpr size_t kS2BoardTMUX = 18;
  static constexpr size_t kEmptyFrames = 0;

  const std::map<l1t::demo::LinkId, std::pair<l1t::demo::ChannelSpec, std::vector<size_t>>>
      kChannelSpecsOutputToL1T = {
          /* logical channel within time slice -> {{link TMUX, inter-packet gap}, vector of channel indices} */
          {{"towersAndClusters", 0}, {{kS2BoardTMUX, kGapLengthOutput}, {0}}},
          {{"towersAndClusters", 1}, {{kS2BoardTMUX, kGapLengthOutput}, {1}}},
          {{"towersAndClusters", 2}, {{kS2BoardTMUX, kGapLengthOutput}, {2}}},
          {{"towersAndClusters", 3}, {{kS2BoardTMUX, kGapLengthOutput}, {3}}}
        };

  // ----------member functions ----------------------
  void produce(edm::Event&, const edm::EventSetup&) override;

  // ----------member data ---------------------------
  l1t::demo::BoardDataReader fileReader_;
};

//
// class implementation
//

Stage2FileReader::Stage2FileReader(const edm::ParameterSet& iConfig)
    : fileReader_(l1t::demo::parseFileFormat(iConfig.getUntrackedParameter<std::string>("format")),
                  iConfig.getParameter<std::vector<std::string>>("files"),
                  kFramesPerTMUXPeriod,
                  kS2BoardTMUX,
                  kEmptyFrames,
                  kChannelSpecsOutputToL1T) {
  produces<l1t::HGCalMulticlusterBxCollection>();
}

// ------------ method called to produce the data  ------------
void Stage2FileReader::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
//   using namespace l1t::demo::codecs;

  l1t::demo::EventData eventData(fileReader_.getNextEvent());

  const unsigned bitsPerWord = 64;
  const unsigned wordsPerCluster = 4;
  std::array<std::vector<ap_uint<bitsPerWord>>, wordsPerCluster> clusterWords;
  for ( unsigned iWord = 0; iWord < wordsPerCluster; ++iWord ) {
      clusterWords.at(iWord) = eventData.at({"towersAndClusters", iWord});
  }


  l1thgcfirmware::decodeClusters(clusterWords); 

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Stage2FileReader::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // Stage2FileReader
  edm::ParameterSetDescription desc;
  desc.add<std::vector<std::string>>("files",
                                     {
                                         "HGCS2OutputToL1TFile_Sector0_0.txt",
                                     });
  desc.addUntracked<std::string>("format", "EMPv2");
  descriptions.add("Stage2FileReader", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Stage2FileReader);