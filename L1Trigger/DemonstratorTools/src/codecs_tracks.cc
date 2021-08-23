
#include "L1Trigger/DemonstratorTools/interface/codecs/tracks.h"

namespace l1t::demo::codecs {

  ap_uint<96> encodeTrack(const TTTrack_TrackWord& t) { return t.getTrackWord(); }

  // Encodes track collection onto 18 output links (2x9 eta-phi sectors; first 9 negative eta)
  std::array<std::vector<ap_uint<64>>, 18> encodeTracks(const edm::View<TTTrack<Ref_Phase2TrackerDigi_>>& tracks) {
    std::array<std::vector<ap_uint<96>>, 18> trackWords;

    for (const auto& track : tracks)
      trackWords.at((track.eta() >= 0 ? 9 : 0) + track.phiSector()).push_back(encodeTrack(track));

    std::array<std::vector<ap_uint<64>>, 18> linkData;

    for (size_t i = 0; i < linkData.size(); i++) {
      // Pad track vectors -> full packet length (156 frames = 104 tracks)
      trackWords.at(i).resize(104, 0);
      linkData.at(i).resize(156, {0});

      for (size_t j = 0; (j < trackWords.size()); j += 2) {
        linkData.at(i).at(3 * j / 2) = trackWords.at(i).at(j)(63, 0);
        linkData.at(i).at(3 * j / 2 + 1) =
            (ap_uint<32>(trackWords.at(i).at(j + 1)(31, 0)), ap_uint<32>(trackWords.at(i).at(j)(95, 64)));
        linkData.at(i).at(3 * j / 2 + 2) = trackWords.at(i).at(j + 1)(95, 32);
      }
    }

    return linkData;
  }

}  // namespace l1t::demo::codecs