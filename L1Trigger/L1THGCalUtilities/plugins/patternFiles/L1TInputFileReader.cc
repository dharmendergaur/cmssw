// -*- C++ -*-
//
// Package:    L1Trigger/L1THGCalUtilities
// Class:      Stage2FileWriter
//
/**\class Stage2FileWriter Stage2FileWriter.cc L1Trigger/L1THGCalUtilities/plugins/patternFiles/Stage2FileWriter.cc

 Description: EDAnalyzer for read I/O buffer files for hardware/firmware tests of HGC stage 2
    Implemented as inputs appear at CTL1, not at output of HGC Stage2

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Emyr Clement
//         Created:  Wed, 20 Apr 2022
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

#include "L1Trigger/DemonstratorTools/interface/BoardDataReader.h"
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "L1Trigger/DemonstratorTools/interface/utilities.h"

//
// class declaration
//

class L1TInputFileReader : public edm::stream::EDProducer<> {
public:
  explicit L1TInputFileReader(const edm::ParameterSet&);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // ----------constants, enums and typedefs ---------
  static constexpr size_t kFramesPerTMUXPeriod = 9;
  static constexpr size_t kGapLength = 0;
  static constexpr size_t kS2BoardTMUX = 18;
  static constexpr size_t kCTL1BoardTMUX = 6;
  static constexpr size_t kMaxLinesPerFile = 1024;
  static constexpr size_t kEmptyFrames = 0;

//   const std::map<std::string, l1t::demo::ChannelSpec> kChannelSpecs = {
//       /* interface name -> {link TMUX, inter-packet gap} */
//       {"towersAndClusters", {kS2BoardTMUX, kGapLength}}};

//   const std::map<l1t::demo::LinkId, std::vector<size_t>> kChannelIds = {
  const std::map<l1t::demo::LinkId, std::pair<l1t::demo::ChannelSpec, std::vector<size_t>>> kChannelIds = {
      /* logical channel within time slice -> vector of channel indices (one entry per time slice) */
      {{"towersAndClusters", 0}, {{kS2BoardTMUX, kGapLength}, {44, 48, 52}}},
      {{"towersAndClusters", 1}, {{kS2BoardTMUX, kGapLength}, {45, 49, 53}}},
      {{"towersAndClusters", 2}, {{kS2BoardTMUX, kGapLength}, {46, 50, 54}}},
      {{"towersAndClusters", 3}, {{kS2BoardTMUX, kGapLength}, {47, 51, 55}}},
      {{"towersAndClusters", 4}, {{kS2BoardTMUX, kGapLength}, {56, 60, 64}}},
      {{"towersAndClusters", 5}, {{kS2BoardTMUX, kGapLength}, {57, 61, 65}}},
      {{"towersAndClusters", 6}, {{kS2BoardTMUX, kGapLength}, {58, 62, 66}}},
      {{"towersAndClusters", 7}, {{kS2BoardTMUX, kGapLength}, {59, 63, 67}}},
      {{"towersAndClusters", 8}, {{kS2BoardTMUX, kGapLength}, {68, 72, 76}}},
      {{"towersAndClusters", 9}, {{kS2BoardTMUX, kGapLength}, {69, 73, 77}}},
      {{"towersAndClusters", 10}, {{kS2BoardTMUX, kGapLength}, {70, 74, 78}}},
      {{"towersAndClusters", 11}, {{kS2BoardTMUX, kGapLength}, {71, 75, 79}}}};
  // ----------member functions ----------------------
  void produce(edm::Event&, const edm::EventSetup&) override;

//   void decodePacket();

  // ----------member data ---------------------------
  l1t::demo::BoardDataReader fileReader_;
};

//
// class implementation
//

L1TInputFileReader::L1TInputFileReader(const edm::ParameterSet& iConfig)
    : fileReader_(l1t::demo::parseFileFormat(iConfig.getUntrackedParameter<std::string>("format")),
                  iConfig.getParameter<std::vector<std::string>>("files"),
                  kFramesPerTMUXPeriod,
                  kCTL1BoardTMUX,
                  kEmptyFrames,
                  kChannelIds) {
    std::cout << "Constructor" << std::endl;
  produces<l1t::HGCalMulticlusterBxCollection>();
}

// ------------ method called to produce the data  ------------
void L1TInputFileReader::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
//   using namespace l1t::demo::codecs;
    std::cout << "In produce" << std::endl;
  l1t::demo::EventData eventData(fileReader_.getNextEvent());
    std::cout << "Got event data" << std::endl;

    const auto& link0 = eventData.at({"towersAndClusters", 0});
    const auto& link1 = eventData.at({"towersAndClusters", 1});
    for ( unsigned int i = 0; i< link0.size(); ++i) {
        if ( i < 31 ) continue;
        const auto& firstWord = link0.at(i);
        const auto& secondWord = link1.at(i);
        if ( firstWord == 0 ) continue;
        std::cout << "Got data : " << std::hex << " " << firstWord << " " << secondWord << std::endl;
        std::cout << "Got data : "  << firstWord.to_string() << " " << secondWord.to_string() << std::endl;

        // Extract pt's
        ap_uint<14> pt = firstWord(13,0);
        std::cout << "Got pt : " << pt << std::endl;

        ap_uint<9> pt_em = firstWord(27,14);
        std::cout << "Got pt em : " << pt_em << std::endl;

        // Extract eta
        constexpr float ETAPHI_LSB = M_PI / 720;
        ap_int<9> eta = secondWord(8,0);
        std::cout << "Got eta : " << eta << " " << eta.to_string() << " " << secondWord(8,0).to_string() << " " << eta * ETAPHI_LSB << std::endl;
        std::cout << ap_int<9>(secondWord(8,0)) << " " << ap_uint<9>(secondWord(8,0)) << std::endl;
        std::cout << ap_int<9>(1) << " " << ap_uint<9>(1) << std::endl;
        // Extract phi
        ap_int<9> phi = secondWord(17,9);
        std::cout << "Got phi : " << phi << " " << phi * ETAPHI_LSB << std::endl;

    }
//   l1t::VertexWordCollection vertices(decodeVertices(eventData.at({"vertices", 0})));

//   edm::LogInfo("L1TInputFileReader") << vertices.size() << " vertices found";

  iEvent.put(std::make_unique<l1t::HGCalMulticlusterBxCollection>());
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void L1TInputFileReader::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // L1TInputFileReader
  edm::ParameterSetDescription desc;
  desc.add<std::vector<std::string>>("files",
                                     {
                                         "input-emp-vu9p.001.txt",
                                     });
  desc.addUntracked<std::string>("format", "EMP");
  descriptions.add("L1TInputFileReader", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TInputFileReader);