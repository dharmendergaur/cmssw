
#include "L1Trigger/DemonstratorTools/interface/codecs/jetMet.h"
#include <numeric>

namespace l1t::demo::codecs {

  ap_uint<64> encodeJet(const reco::CaloJet& j ) { 
    // std::cout << "Jet : " << j.pt() << " " << j.eta() << " " << j.phi() << std::endl;

    ap_ufixed<14, 12, AP_TRN, AP_SAT> pt = j.pt();

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
    ap_uint<10> etaBin = etaFloat / 0.0833;

    // Get eta of bin centre, convert to GT eta
    double etaBinCentre = -3 +1.5*etaRegion+(etaBin+0.5)*0.0833;
    auto gtEta = l1ct::Scales::makeGlbEta( etaBinCentre );
    // std::cout << "Jet : " << j.pt() << " " << j.eta() << " " << etaBin << " " << etaRegion << " " << etaBinCentre << " " << gtEta << std::endl;

    ap_uint<10> phiBin = (j.phi() + 3.14) / 0.0875;
    double phiBinCentre = -3.14 + ( phiBin+0.5 ) * 0.0875;
    auto gtPhi = l1ct::Scales::makeGlbPhi( phiBinCentre );
    std::cout << "Phi : " << j.phi() << " " << phiBin << " " << phiBinCentre << " " << gtPhi << std::endl;
    l1ct::Jet jet;
    jet.hwPt = j.pt();
    jet.hwEta = gtEta;
    jet.hwPhi = gtPhi;

    ap_uint<64> candWord = jet.pack_ap();

    // candWord(14-1, 0) = pt(13, 0);
    // candWord(14+8-1, 14) = eta(7,0);
    // candWord(14+8+8-1, 14+8) = phi(7,0);
    // std::cout << "Digi jet : " << pt << " " << eta << " " << phi << " " << candWord << std::endl;
    // std::cout << "Float eta : " << j.eta() << std::endl;
    // std::cout << "Eta bin : " << eta << std::endl;
    // std::cout << "GT eta : " << l1ct::Scales::makeGlbEta( j.eta() ) << " " << l1ct::Scales::makeGlbEta( -3 ) << " " << l1ct::Scales::makeGlbEta( 0 ) << " " << l1ct::Scales::makeGlbEta( 3 ) << std::endl;
    // std::cout << "Region edges : " << l1ct::Scales::makeGlbEta( -3 ) << " " << l1ct::Scales::makeGlbEta( -1.5 ) << " " << l1ct::Scales::makeGlbEta( 0 ) << " " << l1ct::Scales::makeGlbEta( 1.5 ) << " " << l1ct::Scales::makeGlbEta( 3 ) << std::endl;
    return candWord;
  }

  ap_uint<64> encodeMet(const edm::View<l1t::EtSum>& met ) {
    ap_uint<64> candWord = 0;
    for ( const auto& sum : met ) {
      if ( sum.getType() != l1t::EtSum::EtSumType::kMissingEt ) continue;
      ap_ufixed<14, 12, AP_TRN, AP_SAT> met = sum.pt();
      candWord(14-1, 0) = met(13, 0);
      std::cout << "MET : " << met << std::endl;
    }
    return candWord;
  }

  ap_uint<64> encodeHt(const edm::View<l1t::EtSum>& ht ) {
    ap_uint<64> candWord = 0;
    for ( const auto& sum : ht ) {
      if ( sum.getType() == l1t::EtSum::EtSumType::kTotalHt ) {
        ap_ufixed<14, 12, AP_TRN, AP_SAT> ht = sum.pt();
        candWord(38, 25) = ht(13, 0);
        std::cout << "HT : " << ht << std::endl;
      }
      else if ( sum.getType() == l1t::EtSum::EtSumType::kMissingHt ) {
        ap_ufixed<14, 12, AP_TRN, AP_SAT> mht = sum.pt();
        candWord(14-1, 0) = mht(13, 0);
        std::cout << "MHT : " << mht << std::endl;
      }
    }
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