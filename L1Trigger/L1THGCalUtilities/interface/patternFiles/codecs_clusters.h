#ifndef L1Trigger_L1TriggerUtilities_patternFiles_codecs_clusters_h
#define L1Trigger_L1TriggerUtilities_patternFiles_codecs_clusters_h

#include <array>
#include <vector>

#include "ap_int.h"

#include "DataFormats/Common/interface/View.h"
#include "L1Trigger/DemonstratorTools/interface/BoardData.h"

namespace l1thgcfirmware {
    void decodeClusters(const std::array<std::vector<ap_uint<64>>, 4> &);

    void decodeFirstWord( const ap_uint<64> );
    void decodeSecondWord( const ap_uint<64> );
    void decodeThirdWord( const ap_uint<64> );
    void decodeFourthWord( const ap_uint<64> );

}

#endif