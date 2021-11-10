
#ifndef L1Trigger_DemonstratorTools_codecs_jetMet_h
#define L1Trigger_DemonstratorTools_codecs_jetMet_h

#include <array>
#include <vector>

#include "ap_int.h"
#include "ap_fixed.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"
#include "L1Trigger/Phase2L1ParticleFlow/src/newfirmware/dataformats/jets.h"
#include "L1Trigger/Phase2L1ParticleFlow/src/newfirmware/dataformats/gt_datatypes.h"

#include "L1Trigger/DemonstratorTools/interface/BoardData.h"

namespace l1t::demo::codecs {

  ap_uint<64> encodeJet(const reco::CaloJet& c);
  ap_uint<64> encodeMet(const edm::View<l1t::EtSum>& met );

  std::array<std::vector<ap_uint<64>>, 1> encodeJetMet(const edm::View<reco::CaloJet>&, const edm::View<l1t::EtSum>&, const edm::View<l1t::EtSum>&);

}  // namespace l1t::demo::codecs

#endif