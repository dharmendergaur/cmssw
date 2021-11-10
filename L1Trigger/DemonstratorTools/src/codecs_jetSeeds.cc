
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



    // l1gt::pt_t pt = j.pt();

    // // Get eta bin number and eta region
    // double etaFloat = j.eta();
    // unsigned int etaRegion = 2;

    // if ( etaFloat > 1.5 ) {
    //   etaFloat -= 1.5;
    //   etaRegion = 3;
    // }
    // else if ( etaFloat < 0 ) {
    //   while (etaFloat < 0 ) {
    //     etaFloat += 1.5;
    //     etaRegion -= 1;
    //   }
    // }

    // // Get eta of bin centre, convert to GT eta
    // double etaLSB = 6. / 72;
    // ap_uint<10> etaBin = etaFloat / etaLSB;
    // double etaBinCentre = -3 +1.5*etaRegion+(etaBin+0.5)*(etaLSB);
    // l1gt::eta_t gtEta = etaBinCentre / l1gt::Scales::ETAPHI_LSB;
    // // std::cout << "Jet : " << j.pt() << " " << j.eta() << " " << etaBin << " " << etaRegion << " " << etaBinCentre << " " << gtEta << " " << gtEta.to_string() << std::endl;

    // double phiLSB = 2 * M_PI / 72;
    // ap_uint<10> phiBin = (j.phi() + 3.14) / phiLSB;
    // double phiBinCentre = -3.14 + ( phiBin+0.5 ) * phiLSB;
    // l1gt::phi_t gtPhi = phiBinCentre / l1gt::Scales::ETAPHI_LSB;
    // // std::cout << "Phi : " << j.phi() << " " << phiBin << " " << phiBinCentre << " " << gtPhi << " " << gtPhi.to_string() << std::endl;

    std::cout << "Jet seed, float : " << j.pt() << " " << j.eta() << " " << j.phi() << std::endl;
    std::cout << "Jet seed : " << pt << " " << eta << " " << phi << std::endl;
    // l1gt::Jet jet;
    // jet.valid = 1;
    // jet.v3.pt = pt;
    // jet.v3.eta = gtEta;
    // jet.v3.phi = gtPhi;
    // jet.z0 = 0;
    // ap_uint<64> candWord = jet.pack()[0];

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

}  // namespace l1t::demo::codecs