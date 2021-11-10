// -*- C++ -*-
//
// Package:    L1Trigger/DemonstratorTools
// Class:      CTL2FileWriter
//
/**\class CTL2FileWriter CTL2FileWriter.cc L1Trigger/DemonstratorTools/plugins/CTL2FileWriter.cc

 Description: Example EDAnalyzer class, illustrating how BoardDataWriter can be used to
   write I/O buffer files for hardware/firmware tests

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Thomas Williams <thomas.williams@stfc.ac.uk>
//         Created:  Mon, 15 Feb 2021 00:39:44 GMT
//
//

// system include files
#include <memory>

#include "ap_int.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1TParticleFlow/interface/PFCluster.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"
#include "DataFormats/Common/interface/View.h"

#include "L1Trigger/DemonstratorTools/interface/BoardDataWriter.h"
#include "L1Trigger/DemonstratorTools/interface/codecs/puppiCands.h"
#include "L1Trigger/DemonstratorTools/interface/codecs/jetMet.h"
#include "L1Trigger/DemonstratorTools/interface/utilities.h"

//
// class declaration
//

class CTL2FileWriter : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit CTL2FileWriter(const edm::ParameterSet&);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // ----------constants, enums and typedefs ---------
  // NOTE: At least some of the info from these constants will eventually come from config files
  static constexpr size_t kFramesPerTMUXPeriod = 9;
  static constexpr size_t kGapLength = 0;
  static constexpr size_t kCTL2TMUX = 6;
  static constexpr size_t kMaxLinesPerFile = 1024;

  const std::map<l1t::demo::LinkId, std::vector<size_t>> kChannelIdsInput = {
      /* logical channel within time slice -> vector of channel indices (one entry per time slice) */
      {{"puppiCands", 0}, {24}},
      {{"puppiCands", 1}, {25}},
      {{"puppiCands", 2}, {26}},
      {{"puppiCands", 3}, {27}},
      {{"puppiCands", 4}, {28}},
      {{"puppiCands", 5}, {29}},
      {{"puppiCands", 6}, {30}},
      {{"puppiCands", 7}, {31}},
      {{"puppiCands", 8}, {32}},
      {{"puppiCands", 9}, {33}},
      {{"puppiCands", 10}, {34}},
      {{"puppiCands", 11}, {35}},
      {{"puppiCands", 12}, {36}},
      {{"puppiCands", 13}, {37}},
      {{"puppiCands", 14}, {38}},
      {{"puppiCands", 15}, {39}},
      {{"puppiCands", 16}, {80}},
      {{"puppiCands", 17}, {81}},

      {{"puppiCands", 18}, {82}},
      {{"puppiCands", 19}, {83}},
      {{"puppiCands", 20}, {84}},
      {{"puppiCands", 21}, {85}},
      {{"puppiCands", 22}, {86}},
      {{"puppiCands", 23}, {87}},

      {{"puppiCands", 24}, {88}},
      {{"puppiCands", 25}, {89}},
      {{"puppiCands", 26}, {90}}

      
      };

  const std::map<std::string, l1t::demo::ChannelSpec> kChannelSpecsInput = {
      /* interface name -> {link TMUX, inter-packet gap} */
      {"puppiCands", {kCTL2TMUX, kGapLength}}};

  const std::map<l1t::demo::LinkId, std::pair<l1t::demo::ChannelSpec, std::vector<size_t>>>
      kChannelSpecsOutputToGT = {
          /* logical channel within time slice -> {{link TMUX, inter-packet gap}, vector of channel indices} */
          {{"jetMet", 0}, {{kCTL2TMUX, kGapLength}, {0}}}};

  // ----------member functions ----------------------
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::View<reco::Candidate>> puppiToken_;
  edm::EDGetTokenT<edm::View<reco::CaloJet>> jetsToken_;
  edm::EDGetTokenT<edm::View<l1t::EtSum>> metToken_;
  edm::EDGetTokenT<edm::View<l1t::EtSum>> htToken_;

  l1t::demo::BoardDataWriter fileWriterInputPuppiCands_;

  l1t::demo::BoardDataWriter fileWriterOutputJetsAndSums_;

};

//
// class implementation
//

CTL2FileWriter::CTL2FileWriter(const edm::ParameterSet& iConfig)
    : puppiToken_(consumes<edm::View<reco::Candidate>>(iConfig.getUntrackedParameter<edm::InputTag>("puppiCands"))),
      jetsToken_(consumes<edm::View<reco::CaloJet>>(iConfig.getUntrackedParameter<edm::InputTag>("jets"))),
      metToken_(consumes<edm::View<l1t::EtSum>>(iConfig.getUntrackedParameter<edm::InputTag>("met"))),
      htToken_(consumes<edm::View<l1t::EtSum>>(iConfig.getUntrackedParameter<edm::InputTag>("ht"))),
      fileWriterInputPuppiCands_(l1t::demo::parseFileFormat(iConfig.getUntrackedParameter<std::string>("format")),
                             iConfig.getUntrackedParameter<std::string>("inputFilename"),
                             kFramesPerTMUXPeriod,
                             kCTL2TMUX,
                             kMaxLinesPerFile,
                             kChannelIdsInput,
                             kChannelSpecsInput),
        fileWriterOutputJetsAndSums_(l1t::demo::parseFileFormat(iConfig.getUntrackedParameter<std::string>("format")),
                             iConfig.getUntrackedParameter<std::string>("outputFilename"),
                             kFramesPerTMUXPeriod,
                             kCTL2TMUX,
                             kMaxLinesPerFile,
                             kChannelSpecsOutputToGT ) {}

void CTL2FileWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace l1t::demo::codecs;

  // // 1) Encode track information onto vectors containing link data
  const auto puppiData(encodePuppiCands(iEvent.get(puppiToken_)));
  const auto jetMetData(encodeJetMet(iEvent.get(jetsToken_),iEvent.get(metToken_),iEvent.get(htToken_)));

  // // 2) Pack track information into 'event data' object, and pass that to file writer
  l1t::demo::EventData eventDataPuppiCands;
  for (size_t i = 0; i < 27; i++) {
    eventDataPuppiCands.add({"puppiCands", i}, puppiData.at(i));
  }

  l1t::demo::EventData eventDataJetMet;
  eventDataJetMet.add({"jetMet", 0}, jetMetData.at(0));


  fileWriterInputPuppiCands_.addEvent(eventDataPuppiCands);
  fileWriterOutputJetsAndSums_.addEvent(eventDataJetMet);

}

// ------------ method called once each job just after ending the event loop  ------------
void CTL2FileWriter::endJob() {
  // Writing pending events to file before exiting
  fileWriterInputPuppiCands_.flush();
  fileWriterOutputJetsAndSums_.flush();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void CTL2FileWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CTL2FileWriter);
