// -*- C++ -*-
//
// Package:    L1Trigger/DemonstratorTools
// Class:      CTL2JetSeedWriter
//
/**\class CTL2JetSeedWriter CTL2JetSeedWriter.cc L1Trigger/DemonstratorTools/plugins/CTL2JetSeedWriter.cc

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
#include "DataFormats/Common/interface/View.h"

#include "L1Trigger/DemonstratorTools/interface/BoardDataWriter.h"
#include "L1Trigger/DemonstratorTools/interface/codecs/puppiCands.h"
#include "L1Trigger/DemonstratorTools/interface/codecs/jetSeeds.h"
#include "L1Trigger/DemonstratorTools/interface/utilities.h"

//
// class declaration
//

class CTL2JetSeedWriter : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit CTL2JetSeedWriter(const edm::ParameterSet&);

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

      // CTL1->CTL2 test
      // Barrel, 6*3
      {{"puppiCands", 0}, {62}},
      {{"puppiCands", 1}, {63}},
      {{"puppiCands", 2}, {64}},
      {{"puppiCands", 3}, {65}},
      {{"puppiCands", 4}, {66}},
      {{"puppiCands", 5}, {67}},
      {{"puppiCands", 6}, {68}},
      {{"puppiCands", 7}, {69}},
      {{"puppiCands", 8}, {70}},
      {{"puppiCands", 9}, {71}},
      {{"puppiCands", 10}, {72}},
      {{"puppiCands", 11}, {73}},
      {{"puppiCands", 12}, {74}},
      {{"puppiCands", 13}, {75}},
      {{"puppiCands", 14}, {76}},
      {{"puppiCands", 15}, {77}},
      {{"puppiCands", 16}, {78}},
      {{"puppiCands", 17}, {79}},

      // Endcap with tracks, 3*2
      {{"puppiCands", 18}, {43}},
      {{"puppiCands", 19}, {44}},
      {{"puppiCands", 20}, {45}},
      {{"puppiCands", 21}, {56}},
      {{"puppiCands", 22}, {57}},
      {{"puppiCands", 23}, {58}},

      // Endcap no tracks, 3*1
      {{"puppiCands", 24}, {40}},
      {{"puppiCands", 25}, {41}},
      {{"puppiCands", 26}, {42}}

      
      };

  const std::map<std::string, l1t::demo::ChannelSpec> kChannelSpecsInput = {
      /* interface name -> {link TMUX, inter-packet gap} */
      {"puppiCands", {kCTL2TMUX, kGapLength}}};

  const std::map<l1t::demo::LinkId, std::pair<l1t::demo::ChannelSpec, std::vector<size_t>>>
      kChannelSpecsJetSeedsOutput = {
          /* logical channel within time slice -> {{link TMUX, inter-packet gap}, vector of channel indices} */
          {{"jetSeeds", 0}, {{kCTL2TMUX, kGapLength}, {58}}},
          {{"jetSeeds", 1}, {{kCTL2TMUX, kGapLength}, {59}}},
          {{"jetSeeds", 2}, {{kCTL2TMUX, kGapLength}, {60}}},
          {{"jetSeeds", 3}, {{kCTL2TMUX, kGapLength}, {61}}},
          {{"jetSeeds", 4}, {{kCTL2TMUX, kGapLength}, {62}}},
          {{"jetSeeds", 5}, {{kCTL2TMUX, kGapLength}, {63}}},
          {{"jetSeeds", 6}, {{kCTL2TMUX, kGapLength}, {64}}},
          {{"jetSeeds", 7}, {{kCTL2TMUX, kGapLength}, {65}}},
          {{"jetSeeds", 8}, {{kCTL2TMUX, kGapLength}, {66}}},
          {{"jetSeeds", 9}, {{kCTL2TMUX, kGapLength}, {67}}},
          {{"jetSeeds", 10}, {{kCTL2TMUX, kGapLength}, {68}}},
          {{"jetSeeds", 11}, {{kCTL2TMUX, kGapLength}, {69}}}          
      };

  // ----------member functions ----------------------
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::View<reco::Candidate>> puppiToken_;
  edm::EDGetTokenT<edm::View<reco::Candidate>> jetSeedsToken_;

  l1t::demo::BoardDataWriter fileWriterInputPuppiCands_;

  l1t::demo::BoardDataWriter fileWriterOutputJetSeeds_;

};

//
// class implementation
//

CTL2JetSeedWriter::CTL2JetSeedWriter(const edm::ParameterSet& iConfig)
    : puppiToken_(consumes<edm::View<reco::Candidate>>(iConfig.getUntrackedParameter<edm::InputTag>("puppiCands"))),
      jetSeedsToken_(consumes<edm::View<reco::Candidate>>(iConfig.getUntrackedParameter<edm::InputTag>("jetSeeds"))),
      fileWriterInputPuppiCands_(l1t::demo::parseFileFormat(iConfig.getUntrackedParameter<std::string>("format")),
                             iConfig.getUntrackedParameter<std::string>("inputFilename"),
                             kFramesPerTMUXPeriod,
                             kCTL2TMUX,
                             kMaxLinesPerFile,
                             kChannelIdsInput,
                             kChannelSpecsInput),
        fileWriterOutputJetSeeds_(l1t::demo::parseFileFormat(iConfig.getUntrackedParameter<std::string>("format")),
                             iConfig.getUntrackedParameter<std::string>("outputFilename"),
                             kFramesPerTMUXPeriod,
                             kCTL2TMUX,
                             kMaxLinesPerFile,
                             kChannelSpecsJetSeedsOutput ) {}

void CTL2JetSeedWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace l1t::demo::codecs;

  // // 1) Encode track information onto vectors containing link data
  const auto puppiData(encodePuppiCands(iEvent.get(puppiToken_)));
  const auto jetSeedsData(encodeJetSeeds(iEvent.get(jetSeedsToken_)));

  // // 2) Pack track information into 'event data' object, and pass that to file writer
  l1t::demo::EventData eventDataPuppiCands;
  for (size_t i = 0; i < 27; i++) {
    eventDataPuppiCands.add({"puppiCands", i}, puppiData.at(i));
  }

  l1t::demo::EventData eventDataJetSeeds;
  for (size_t i = 0; i < 12; i++) {
    eventDataJetSeeds.add({"jetSeeds", i}, jetSeedsData.at(i));
  }

  fileWriterInputPuppiCands_.addEvent(eventDataPuppiCands);
  fileWriterOutputJetSeeds_.addEvent(eventDataJetSeeds);

}

// ------------ method called once each job just after ending the event loop  ------------
void CTL2JetSeedWriter::endJob() {
  // Writing pending events to file before exiting
  fileWriterInputPuppiCands_.flush();
  fileWriterOutputJetSeeds_.flush();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void CTL2JetSeedWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(CTL2JetSeedWriter);
