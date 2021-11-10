
#include "L1Trigger/DemonstratorTools/interface/codecs/puppiCands.h"
#include <numeric>

namespace l1t::demo::codecs {

  ap_uint<64> encodePuppiCand(const reco::Candidate& c, float phiEdge, float etaEdge ) { 
    ap_ufixed<14, 12, AP_TRN, AP_SAT> pt = c.pt();
    ap_uint<10> phi = (c.phi() - phiEdge ) / 0.0043633231;
    ap_uint<10> eta = (c.eta() - etaEdge ) / 0.0043633231;
    ap_uint<64> candWord = 0;
    // if ( abs(c.pt() - 8.75) < 0.1 ) {
      // std::cout << "Encoding puppi cand : " << c.pt() << " " << c.phi() << " " << c.eta() << std::endl;
      // std::cout << phiEdge << " " << etaEdge << std::endl;
      // std::cout << pt << " " << phi << " " << eta << std::endl;
    // }

    candWord(14-1, 0) = pt(13, 0);
    candWord(14+10-1, 14) = phi(9,0);
    candWord(14+10+10-1, 14+10) = eta(9,0);
    // std::cout << candWord << std::endl;
    return candWord;
  }

  std::array<std::vector<ap_uint<64>>, 27> encodePuppiCands(const edm::View<reco::Candidate>& puppiCands) {
    std::array<std::vector<ap_uint<64>>, 27> linkData;
    for (size_t i = 0; i < linkData.size(); i++) {
      linkData.at(i).resize(54, {0});
      // linkData.at(i).at(0) = 5839227683695833246;
    }

    // Sort puppi candidates into input regions
    // What's currently done by Phase1L1TJetProducer::prepareInputsIntoRegions in emulator
    std::vector<float> etaRegionEdges{-3, -2.5, -1.5, -1.0, -0.5, 0, 0.5, 1, 1.5, 2.5, 3};
    // std::vector<float> phiRegionEdges{-3.5, -2.8, -2.1, -1.4, -0.7, 0, 0.7, 1.4, 2.1, 2.8, 3.5};
    std::vector<float> phiRegionEdges{-3.15, -2.45, -1.75, -1.05, -0.35, 0.35, 1.05, 1.75, 2.45, 3.15};
    std::vector<std::vector<ap_uint<64>>> inputsInRegions{ (etaRegionEdges.size()-1) * (phiRegionEdges.size() - 1)};
    for (const auto& puppiCand : puppiCands ) {

      if (puppiCand.phi() < phiRegionEdges.front() || puppiCand.phi() >= phiRegionEdges.back() ||
          puppiCand.eta() < etaRegionEdges.front() || puppiCand.eta() >= etaRegionEdges.back())
        continue;
      // std::cout << "Input puppi cand : " << puppiCand.pt() << " " << puppiCand.eta() << " " << puppiCand.phi() << std::endl;
      // Which phi region does this cand belong to
      auto it_phi = phiRegionEdges.begin();
      it_phi = std::upper_bound(phiRegionEdges.begin(), phiRegionEdges.end(), puppiCand.phi()) - 1;

      // Which eta region does this cand belong to
      auto it_eta = etaRegionEdges.begin();
      it_eta = std::upper_bound(etaRegionEdges.begin(), etaRegionEdges.end(), puppiCand.eta()) - 1;

      if (it_phi != phiRegionEdges.end() && it_eta != etaRegionEdges.end()) {
        auto phiRegion = it_phi - phiRegionEdges.begin();
        auto etaRegion = it_eta - etaRegionEdges.begin();
        unsigned regionIndex = etaRegion + phiRegion * (etaRegionEdges.size() - 1);
        // std::cout << "Packed cand : " << encodePuppiCand(puppiCand, *it_phi, *it_eta) << std::endl;
        // std::cout << *it_phi << " " << *it_eta << std::endl;
        // std::cout << "Region index : " << regionIndex << " " << etaRegion << " " << phiRegion << std::endl;
        inputsInRegions[regionIndex].push_back(encodePuppiCand(puppiCand, *it_phi, *it_eta));
      }
    }

    // Truncate number of inputs in each pf region
    unsigned maxInputsPerRegion = 18;
    for (auto& inputs : inputsInRegions) {
      if (inputs.size() > maxInputsPerRegion) {
        inputs.resize(maxInputsPerRegion);
      }
    }

    for ( unsigned iRegion=0; iRegion<inputsInRegions.size(); ++iRegion ) {
      // Get links and frames for this region
      unsigned etaRegionIndex = iRegion % (etaRegionEdges.size()-1);
      unsigned phiRegionIndex = iRegion / (etaRegionEdges.size()-1);
      unsigned firstLink=0;
      unsigned nLinks=0;
      unsigned firstFrame=0;
      unsigned nFrames=0;
      if ( etaRegionIndex == 0 || etaRegionIndex == 9 ) { // Endcap no tracks
        firstLink = 24;
        nLinks = 3;
        firstFrame = (etaRegionIndex == 0) ? phiRegionIndex*6 : phiRegionIndex*6 + 3;
        nFrames = 3;
      }
      else if ( etaRegionIndex == 1 || etaRegionIndex == 8 ) { // Endcap with tracks
        firstLink = (etaRegionIndex == 1) ? 18 : 21;
        nLinks = 3;
        firstFrame=phiRegionIndex*6;
        nFrames = 6;
      }
      else { // Barrel
        unsigned bigRegion = ( etaRegionIndex - 2 ) / 2;
        unsigned smallRegion = ( etaRegionIndex - 2 ) % 2;
        firstLink = bigRegion * 6;
        nLinks = 6;
        firstFrame = (smallRegion == 0) ? phiRegionIndex*6 : phiRegionIndex*6 + 3;
        nFrames = 3;
      }

      std::vector<unsigned> linksForRegion( nLinks, 0);
      std::vector<unsigned> framesForRegion( nFrames, 0);
      std::iota( linksForRegion.begin(), linksForRegion.end(), firstLink );
      std::iota( framesForRegion.begin(), framesForRegion.end(), firstFrame );

      // std::cout << "Regions : " << iRegion << " " << etaRegionIndex << " " << phiRegionIndex << std::endl;
      // std::cout << "N links, frames : " << nLinks << " " << nFrames << std::endl;
      // std::cout << "Links : ";
      // for ( const auto& link : linksForRegion ) std::cout << link << " ";
      // std::cout << std::endl;
      // std::cout << "Frames : ";
      // for ( const auto& frame : framesForRegion ) std::cout << frame << " ";
      // std::cout << std::endl;
      // std::cout << std::endl;

      const auto& cands = inputsInRegions.at(iRegion);
      unsigned iCand=0;
      for ( const auto& cand : cands ) {
        unsigned link = iCand % nLinks + firstLink;
        unsigned frame = iCand / nLinks + firstFrame;
        // std::cout << "Input cand : " << iCand << " " << cand << " " << link << " " << frame << std::endl;
        // std::cout << std::hex << cand << std::dec << std::endl;
        ++iCand;
        linkData.at(link).at(frame) = cand;
      }
      // std::cout << std::endl;
      // std::cout << std::endl;
    }

    return linkData;
  }

}  // namespace l1t::demo::codecs