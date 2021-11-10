
#include "L1Trigger/DemonstratorTools/interface/codecs/jetSeeds.h"
#include <numeric>

namespace l1t::demo::codecs {

  ap_uint<64> encodeJetSeed(const reco::Candidate& j ) { 

    l1ct::pt_t pt = j.pt();
    l1ct::glbeta_t eta = l1ct::Scales::makeGlbEta( j.eta() );
    l1ct::glbphi_t phi = l1ct::Scales::makeGlbPhi( j.phi() );

    l1ct::PuppiObj puppiCand;
    puppiCand.clear();
    puppiCand.hwPt = pt;
    puppiCand.hwEta = eta;
    puppiCand.hwPhi = phi;
    ap_uint<64> candWord = puppiCand.pack();

    return candWord;
  }

  std::array<std::vector<ap_uint<64>>, 12> encodeJetSeeds(const edm::View<reco::Candidate>& jetSeeds) {
    std::array<std::vector<ap_uint<64>>, 12> linkData;
    for ( unsigned iLink=0; iLink<12; ++iLink ) {
      linkData.at(iLink).resize(54, {0});
    }

    unsigned int link = 0;
    for ( const auto& jetSeed : jetSeeds ) {
      linkData.at(link).at(0) = encodeJetSeed(jetSeed);
      link += 1;
    }

    return linkData;
  }

std::array<std::vector<ap_uint<64>>, 12> encodeJets(const edm::View<l1t::PFJet>& jets) {
    std::array<std::vector<ap_uint<64>>, 12> linkData;

    for ( const auto& jet : jets ) {
      const auto gtJet = jet.getHWJetGT();
    }

    return linkData;
  }
}  // namespace l1t::demo::codecs