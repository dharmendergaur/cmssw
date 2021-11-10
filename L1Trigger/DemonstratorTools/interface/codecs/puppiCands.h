
#ifndef L1Trigger_DemonstratorTools_codecs_puppiCands_h
#define L1Trigger_DemonstratorTools_codecs_puppiCands_h

#include <array>
#include <vector>

#include "ap_int.h"
#include "ap_fixed.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"

#include "L1Trigger/DemonstratorTools/interface/BoardData.h"

namespace l1t::demo::codecs {

  ap_uint<64> encodePuppiCand(const reco::Candidate& c, unsigned phiEdge, unsigned etaEdge);

  std::array<std::vector<ap_uint<64>>, 27> encodePuppiCands(const edm::View<reco::Candidate>&);

}  // namespace l1t::demo::codecs

#endif