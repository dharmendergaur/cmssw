
#ifndef L1Trigger_DemonstratorTools_codecs_jetSeeds_h
#define L1Trigger_DemonstratorTools_codecs_jetSeeds_h

#include <array>
#include <vector>

#include "ap_int.h"
#include "ap_fixed.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1TParticleFlow/interface/PFJet.h"
#include "DataFormats/L1TParticleFlow/interface/puppi.h"
#include "L1Trigger/DemonstratorTools/interface/BoardData.h"

namespace l1t::demo::codecs {

  ap_uint<64> encodeJetSeeds(const reco::Candidate& c);

  std::array<std::vector<ap_uint<64>>, 12> encodeJetSeeds(const edm::View<reco::Candidate>&);

  std::array<std::vector<ap_uint<64>>, 12> encodeJets(const edm::View<l1t::PFJet>&);

}  // namespace l1t::demo::codecs

#endif