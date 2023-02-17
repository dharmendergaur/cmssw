#include "L1Trigger/L1THGCalUtilities/interface/patternFiles/codecs_clusters.h"
#include <iostream>

namespace l1thgcfirmware {
    std::vector<HGCalCluster_HW> decodeClusters(const std::array<std::vector<ap_uint<64>>, 4>& frames) {

        std::vector<HGCalCluster_HW> clusters;
        unsigned nFrames = frames.at(0).size();
        for ( unsigned int iFrame=0; iFrame < nFrames; ++iFrame ) {
            if ( iFrame < 31 ) continue; // Skip header and towers

            ClusterWords clusterWords;
            for ( unsigned int i=0; i<nWordsPerCluster; ++i ) {
                clusterWords[i] = frames.at(i).at(iFrame);
            }
            HGCalCluster_HW cluster = HGCalCluster_HW::unpack( clusterWords );
            // std::cout << "Got a cluster on frame : " << iFrame << " " << cluster.e << std::endl;
            if (cluster.e > 0 ) clusters.emplace_back( cluster );
        }
        return clusters;
    }
}

//     void decodeFirstWord( const ap_uint<64> word) {
//         // First word
//         const unsigned nb_e = 14;
//         const unsigned nb_e_em = 14;
//         const unsigned nb_gctEGSelectBits = 4;
//         const unsigned nb_fractionInCE_E = 8;
//         const unsigned nb_fractionInCoreCE_E = 8;
//         const unsigned nb_fractionInEarlyCE_H = 8;
//         const unsigned nb_firstLayer = 6;
//         const unsigned nb_spare = 2;

//         unsigned int lsb = 0;
//         unsigned int msb = nb_e-1;
//         ap_uint<nb_e> hw_e = word.range(msb,lsb);
//         // std::cout << "Cluster E : " << lsb << " " << msb << " " << word.to_string() << " " << hw_e.to_string() << std::endl;

//         lsb += nb_e;
//         msb += nb_e_em;
//         ap_uint<nb_e_em> hw_e_em = word.range(msb,lsb);

//         lsb += nb_e_em;
//         msb += nb_gctEGSelectBits;
//         ap_uint<nb_gctEGSelectBits> hw_gctEGSelectBits = word.range(msb,lsb);

//         lsb += nb_gctEGSelectBits;
//         msb += nb_fractionInCE_E;
//         ap_uint<nb_fractionInCE_E> hw_fractionInCE_E = word.range(msb,lsb);

//         lsb += nb_fractionInCE_E;
//         msb += nb_fractionInCoreCE_E;
//         ap_uint<nb_fractionInCoreCE_E> hw_fractionInCoreCE_E = word.range(msb,lsb);

//         lsb += nb_fractionInCoreCE_E;
//         msb += nb_fractionInEarlyCE_H;
//         ap_uint<nb_fractionInEarlyCE_H> hw_fractionInEarlyCE_H = word.range(msb,lsb);

//         lsb += nb_fractionInEarlyCE_H;
//         msb += nb_firstLayer;
//         ap_uint<nb_firstLayer> hw_firstLayer = word.range(msb,lsb);

//         if ( word != 0 ) {
//             // std::cout << "First cluster word : " << hw_firstLayer.to_string() << " " << hw_fractionInEarlyCE_H.to_string() << " " << hw_fractionInCoreCE_E.to_string() << " " << hw_fractionInCE_E.to_string() << " " << hw_gctEGSelectBits.to_string() << " " << hw_e_em.to_string() << " " << hw_e.to_string() << std::endl;
//             std::cout << "First cluster word : " << hw_fractionInEarlyCE_H.to_string() << " " << hw_fractionInCoreCE_E.to_string() << " " << hw_fractionInCE_E.to_string() << " " << hw_gctEGSelectBits.to_string() << " " << hw_e_em.to_string() << " " << hw_e.to_string() << std::endl;
//         }

//     }

//     void decodeSecondWord( const ap_uint<64> word) {
//         const unsigned nb_eta = 10;
//         const unsigned nb_phi = 9;
//         const unsigned nb_z = 12;
//         const unsigned nb_spare_0 = 1;
//         const unsigned nb_nTC = 10;
//         const unsigned nb_satTC = 1;
//         const unsigned nb_qualFracCE_E = 1;
//         const unsigned nb_qualFracCoreCE_E = 1;
//         const unsigned nb_qualFracEarlyCE_H = 1;
//         const unsigned nb_qualSigmasMeans = 1;
//         const unsigned nb_satPhi = 1;
//         const unsigned nb_nominalPhi = 1;
//         const unsigned nb_spare_1 = 15;

//         unsigned int lsb = 0;
//         unsigned int msb = nb_eta-1;
//         ap_uint<nb_eta> hw_eta = word.range(msb,lsb);

//         lsb += nb_eta;
//         msb += nb_phi;
//         ap_int<nb_phi> hw_phi = word.range(msb,lsb);

//         lsb += nb_phi;
//         msb += nb_z;
//         ap_uint<nb_z> hw_z = word.range(msb,lsb);

//         lsb += nb_z;
//         msb += nb_spare_0;
//         ap_uint<nb_spare_0> hw_spare_0 = word.range(msb,lsb);

//         lsb += nb_spare_0;
//         msb += nb_nTC;
//         ap_uint<nb_nTC> hw_nTC = word.range(msb,lsb);

//         lsb += nb_nTC;
//         msb += nb_satTC;
//         ap_uint<nb_satTC> hw_satTC = word.range(msb,lsb);

//         lsb += nb_satTC;
//         msb += nb_qualFracCE_E;
//         ap_uint<nb_qualFracCE_E> hw_qualFracCE_E = word.range(msb,lsb);

//         lsb += nb_qualFracCE_E;
//         msb += nb_qualFracCoreCE_E;
//         ap_uint<nb_qualFracCoreCE_E> hw_qualFracCoreCE_E = word.range(msb,lsb);

//         lsb += nb_qualFracCoreCE_E;
//         msb += nb_qualFracEarlyCE_H;
//         ap_uint<nb_qualFracEarlyCE_H> hw_qualFracEarlyCE_H = word.range(msb,lsb);

//         lsb += nb_qualFracEarlyCE_H;
//         msb += nb_qualSigmasMeans;
//         ap_uint<nb_qualSigmasMeans> hw_qualSigmasMeans = word.range(msb,lsb);

//         lsb += nb_qualSigmasMeans;
//         msb += nb_satPhi;
//         ap_uint<nb_satPhi> hw_satPhi = word.range(msb,lsb);

//         lsb += nb_satPhi;
//         msb += nb_nominalPhi;
//         ap_uint<nb_nominalPhi> hw_nominalPhi = word.range(msb,lsb);

//         lsb += nb_nominalPhi;
//         msb += nb_spare_1;
//         ap_uint<nb_spare_1> hw_spare_1 = word.range(msb,lsb);

//         if ( word != 0 ) {
//             std::cout << "Second cluster word : " << hw_spare_1.to_string() << " " << hw_nominalPhi.to_string() << " " << hw_satPhi.to_string() << " " << hw_qualSigmasMeans.to_string() << " " << hw_qualFracEarlyCE_H.to_string() << " " << hw_qualFracCoreCE_E.to_string() << " " << hw_qualFracCE_E.to_string() << " " << hw_satTC.to_string() << " " << hw_nTC.to_string() << " " << hw_spare_0.to_string() << " " << hw_z.to_string() << " " << hw_phi.to_string() << " " << hw_eta.to_string() << std::endl;
//         }
//     }

//     void decodeThirdWord( const ap_uint<64> word) {
//         const unsigned nb_sigmaE = 7;
//         const unsigned nb_lastLayer = 6;
//         const unsigned nb_showerLength = 6;
//         const unsigned nb_spare = 13;
//         const unsigned nb_sigmaZZ = 7;
//         const unsigned nb_sigmaPhiPhi = 7;
//         const unsigned nb_coreShowerLength = 6;
//         const unsigned nb_sigmaEtaEta = 5;
//         const unsigned nb_sigmaRozRoz = 7;

//         unsigned int lsb = 0;
//         unsigned int msb = nb_sigmaE-1;
//         ap_uint<nb_sigmaE> hw_sigmaE = word.range(msb,lsb);
    
//         lsb += nb_sigmaE;
//         msb += nb_lastLayer;
//         ap_int<nb_lastLayer> hw_lastLayer = word.range(msb,lsb);
    
//         lsb += nb_lastLayer;
//         msb += nb_showerLength;
//         ap_uint<nb_showerLength> hw_showerLength = word.range(msb,lsb);
    
//         lsb += nb_showerLength;
//         msb += nb_spare;
//         ap_uint<nb_spare> hw_spare = word.range(msb,lsb);
    
//         lsb += nb_spare;
//         msb += nb_sigmaZZ;
//         ap_uint<nb_sigmaZZ> hw_sigmaZZ = word.range(msb,lsb);
    
//         lsb += nb_sigmaZZ;
//         msb += nb_sigmaPhiPhi;
//         ap_uint<nb_sigmaPhiPhi> hw_sigmaPhiPhi = word.range(msb,lsb);
    
//         lsb += nb_sigmaPhiPhi;
//         msb += nb_coreShowerLength;
//         ap_uint<nb_coreShowerLength> hw_coreShowerLength = word.range(msb,lsb);

//         lsb += nb_coreShowerLength;
//         msb += nb_sigmaEtaEta;
//         ap_int<nb_sigmaEtaEta> hw_sigmaEtaEta = word.range(msb,lsb);

//         lsb += nb_sigmaEtaEta;
//         msb += nb_sigmaRozRoz;
//         ap_int<nb_sigmaRozRoz> hw_sigmaRozRoz = word.range(msb,lsb);

//         if ( word != 0 ) {
//             // std::cout << "Third cluster word : " << hw_sigmaRozRoz.to_string() << " " << hw_sigmaEtaEta.to_string() << " " << hw_coreShowerLength.to_string() << " " << hw_sigmaPhiPhi.to_string() << " " << hw_sigmaZZ.to_string() << " " << hw_spare.to_string() << " " << hw_showerLength.to_string() << " " << hw_lastLayer.to_string() << " " << hw_sigmaE.to_string() << std::endl;
//             std::cout << "Third cluster word : " << hw_sigmaRozRoz.to_string() << " " << hw_sigmaEtaEta.to_string() << " " << hw_sigmaPhiPhi.to_string() << " " << hw_sigmaZZ.to_string() << " " << hw_sigmaE.to_string() << std::endl;

//         }
//     }

// void decodeFourthWord( const ap_uint<64> word) {
//         // const unsigned nb_sigmaEtaEta = 9;
//         // const unsigned nb_sigmaRozRoz = 13;
//         // const unsigned nb_spare10Bits = 10;
//         // const unsigned nb_spare32Bits = 32;

//         unsigned int lsb = 0;
//         unsigned int msb = 63;
//         ap_uint<64> hw_spare = word.range(msb,lsb);

//         if ( word != 0 ) {
//             std::cout << "Fourth cluster word : " << hw_spare.to_string() << std::endl;
//         }
//     }
// }