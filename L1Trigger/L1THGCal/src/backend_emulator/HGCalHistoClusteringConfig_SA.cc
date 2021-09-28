#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringConfig_SA.h"
#include <math.h>
#include <iostream>
using namespace std;
using namespace l1thgcfirmware;
ClusterAlgoConfig::ClusterAlgoConfig() :
  histogramOffset_(216), // All constants to be read from config file
  clusterizerOffset_(276),
  cClocks_(250),
  cInputs_(72),
  cInputs2_(75),
  cInt_(90),
  cColumns_(108),
  cRows_(44),
  kernelWidths_(),
  areaNormalizations_(),
  thresholdMaximaConstants_(),
  cosLUT_(),
  clusterizerMagicTime_(434),
  stepLatency_({
    { Input  ,  0 },
    { Dist0  ,  1 },
    { Dist1  ,  1 },
    { Dist2  ,  6 },
    { Dist3  ,  0 },
    { Dist4  ,  1 },
    { Dist5  ,  7 },
    { TcToHc , 1 },
    { Hist   , histogramOffset_ + 1 },
    { Smearing1D , 5 },
    { NormArea   , 3 },
    { Smearing2D , 6 },
    { Maxima1D   , 0 },
    { Maxima2D   , 2 },
    { CalcAverage , 6 },
    { Clusterizer , 0 },
    { TriggerCellToCluster , 18 },
    { ClusterSum , 0 }
  }),
  depths_({0 , // No zero layer
          0 , 30 , 59 , 89 , 118 , 148 , 178 , 208 , 237 , 267 , 297 , 327 , 356 , 386 , // CE-E
          415 , 445 , 475 , 505 , 534 , 564 , 594 , 624 , 653 , 683 , 712 , 742 , 772 , 802 , // CE-E
          911 , 1020 , 1129 , 1238 , 1347 , 1456 , 1565 , 1674 , // CE-FH
          1783 , 1892 , 2001 , 2110 , 2281 , 2452 , 2623 , 2794 , 2965 , 3136 , 3307 , 3478 , 3649 , 3820 , // CE-BH
          0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0}), // Padding
  triggerLayers_({0 , // No zero layer
                  1 , 0 , 2 , 0 , 3 , 0 , // CE-E (early)
                  4 , 0 , 5 , 0 , 6 , 0 , 7 , 0 , 8 , 0 , // CE-E (core)
                  9 , 0 , 10 , 0 , 11 , 0 , 12 , 0 , 13 , 0 , 14 , 0 , // CE-E (back)
                  15 , 16 , 17 , 18 , // CE-H (early)
                  19 , 20 , 21 , 22 , 23 , 24 , 25 , 26 , 27 , 28 , 29 , 30 , 31 , 32 , 33 , 34 , 35 , 36 , 0 , 0 , // CE-H (back)
                  0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0  }), // Padding
  layerWeights_E_({0 , // No zero layer
                 1 , 1 , 1 , // CE-E (early)
                 1 , 1 , 1 , 1 , 1 , // CE-E (core)
                 1 , 1 , 1 , 1 , 1 , 1 , // CE-E (back)
                 1 , 1 , 1 , 1 , // CE-H (early)
                 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1}), // CE-H (back)
  layerWeights_E_EM_({0 , // No zero layer
                    252969 , 254280 , 255590 , // CE-E (early)
                    256901 , 258212 , 259523 , 260833 , 262144 , // CE-E (core)
                    263455 , 264765 , 266076 , 267387 , 268698 , 270008 , // CE-E (back)
                    0 , 0 , 0 , 0 , // CE-H (early)
                    0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }), // CE-H (back)
  layerWeights_E_EM_core_({0 , // No zero layer
                         0 , 0 , 0 , // CE-E (early)
                         256901 , 258212 , 259523 , 260833 , 262144 , // CE-E (core)
                         0 , 0 , 0 , 0 , 0 , 0 , // CE-E (back)
                         0 , 0 , 0 , 0 , // CE-H (early)
                         0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }), // CE-H (back)
 layerWeights_E_H_early_({ 0 , // No zero layer
                         0 , 0 , 0 , // CE-E (early)
                         0 , 0 , 0 , 0 , 0 , // CE-E (core)
                         0 , 0 , 0 , 0 , 0 , 0 , // CE-E (back)
                         1 , 1 , 1 , 1 , // CE-H (early)
                         0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0}), // CE-H 
  correction_(131071), // 0b011111111111111111
  saturation_(524287) // (2 ** 19) - 1
{
  initializeSmearingKernelConstants( cRows_, 440, 64 ); // Magic numbers
  initializeThresholdMaximaConstants( cRows_ );
  initializeCosLUT();
}

unsigned int ClusterAlgoConfig::getLatencyUpToAndIncluding( const Step step ) {
  unsigned int latency = 0;
  for ( int iStep = 0; iStep <= step; ++iStep ) latency += getStepLatency(Step(iStep));
  return latency;
}

void ClusterAlgoConfig::initializeSmearingKernelConstants( unsigned int bins, unsigned int offset, unsigned int height ) {
  const unsigned int lWidth0 = offset + (0.5*height);
  const unsigned int lTarget = int( 6.5*lWidth0 - 0.5 );

  for ( unsigned int iBin = 0; iBin < bins; ++iBin ) {
    unsigned int lCentre = lWidth0 + ( height * iBin );
    const unsigned int lBins = int( round(1.0 * lTarget / lCentre) );

    kernelWidths_.push_back( lBins );

    lCentre *= lBins;

    const unsigned int lRatio = int( round(1.0*lTarget/lCentre * pow(2,17) ) ); // Magic numbers

    areaNormalizations_.push_back( lRatio );
  }
}

void ClusterAlgoConfig::initializeThresholdMaximaConstants( unsigned int bins ) {
  unsigned int a = 18000; // Magic numbers
  unsigned int b = 800; // Magic numbers
  int c = -20; // Magic numbers

  for ( unsigned int iBin = 0; iBin < bins; ++iBin ) {
    int threshold = a + b*iBin + c*iBin*iBin;
    thresholdMaximaConstants_.push_back( threshold );
  }
}

void ClusterAlgoConfig::initializeCosLUT() {
  for ( unsigned int iBin = 0; iBin < 78; ++iBin ) { // Magic numbers
    unsigned int cosBin = round( pow(2,18) * ( 1 - cos(iBin*M_PI/1944) ) ); // Magic numbers
    cosLUT_.push_back( cosBin );
  }
}


