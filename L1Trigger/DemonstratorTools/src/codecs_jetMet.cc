
#include "L1Trigger/DemonstratorTools/interface/codecs/jetMet.h"
#include <numeric>

namespace l1t::demo::codecs {

  ap_uint<64> encodeJet(const reco::CaloJet& j ) { 

    l1gt::pt_t pt = j.pt();

    // Get eta bin number and eta region
    double etaFloat = j.eta();
    unsigned int etaRegion = 2;

    if ( etaFloat > 1.5 ) {
      etaFloat -= 1.5;
      etaRegion = 3;
    }
    else if ( etaFloat < 0 ) {
      while (etaFloat < 0 ) {
        etaFloat += 1.5;
        etaRegion -= 1;
      }
    }

    // Get eta of bin centre, convert to GT eta
    double etaLSB = 6. / 72;
    ap_uint<10> etaBin = etaFloat / etaLSB;
    double etaBinCentre = -3 +1.5*etaRegion+(etaBin+0.5)*(etaLSB);
    l1gt::eta_t gtEta = etaBinCentre / l1gt::Scales::ETAPHI_LSB;

    double phiLSB = 2 * M_PI / 72;
    ap_uint<10> phiBin = (j.phi() + 3.14) / phiLSB;
    double phiBinCentre = -3.14 + ( phiBin+0.5 ) * phiLSB;
    l1gt::phi_t gtPhi = phiBinCentre / l1gt::Scales::ETAPHI_LSB;

    std::cout << "Jet : " << pt << " " << etaBinCentre << " " << phiBinCentre << std::endl;
    std::cout << pt.to_string() << " " << gtEta << " " << gtPhi << std::endl;
    l1gt::Jet jet;
    jet.valid = 1;
    jet.v3.pt = pt;
    jet.v3.eta = gtEta;
    jet.v3.phi = gtPhi;
    jet.z0 = 0;
    ap_uint<64> candWord = jet.pack()[0];

    return candWord;
  }

  ap_uint<64> encodeMet(const edm::View<l1t::EtSum>& met ) {
    ap_uint<64> candWord = 0;
    l1gt::Sum gtMet;
    gtMet.valid = 1;
    gtMet.vector_phi = 0;
    for ( const auto& sum : met ) {
      if ( sum.getType() == l1t::EtSum::EtSumType::kMissingEt ) {
        gtMet.vector_pt = sum.pt();
      }
      else if ( sum.getType() == l1t::EtSum::EtSumType::kTotalEt ) {
        gtMet.scalar_pt = sum.pt();
      }
    }
    candWord = gtMet.pack();
    candWord(63,46) = 0;
    std::cout << "MET word : " << candWord << std::endl;
    return candWord;
  }

  ap_uint<64> encodeHt(const edm::View<l1t::EtSum>& ht ) {
    ap_uint<64> candWord = 0;
    l1gt::Sum gtHt;
    gtHt.valid = 1;
    gtHt.vector_phi = 0;
    for ( const auto& sum : ht ) {
      if ( sum.getType() == l1t::EtSum::EtSumType::kTotalHt ) {
        gtHt.scalar_pt = sum.pt();
      }
      else if ( sum.getType() == l1t::EtSum::EtSumType::kMissingHt ) {
        gtHt.vector_pt = sum.pt();
      }
    }
    candWord = gtHt.pack();
    candWord(63,46) = 0;

    return candWord;
  }

  std::array<std::vector<ap_uint<64>>, 1> encodeJetMet(const edm::View<reco::CaloJet>& jets, const edm::View<l1t::EtSum>& met, const edm::View<l1t::EtSum>& ht) {
    std::array<std::vector<ap_uint<64>>, 1> linkData;
    linkData.at(0).resize(54, {0});

    unsigned int frame = 0;
    for ( const auto& jet : jets ) {
      linkData.at(0).at(frame) = encodeJet(jet);
      frame += 2;
    }

    ap_uint<64> htWord = encodeHt( ht );
    linkData.at(0).at(24) = htWord;

    ap_uint<64> metWord = encodeMet( met );
    linkData.at(0).at(25) = metWord;

    return linkData;
  }

}  // namespace l1t::demo::codecs