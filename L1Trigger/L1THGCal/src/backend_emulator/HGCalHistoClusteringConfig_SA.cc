#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringConfig_SA.h"
#include <iostream>
#include <fstream>

using namespace std;
using namespace l1thgcfirmware;
ClusterAlgoConfig::ClusterAlgoConfig() :
  // clusterizerOffset_(0),
  TriggerCellDistributionLUT_( OpenMif("S2.CombinedTD.Balanced60.MixedTypes.NoSplit.mif") ),
  TriggerCellAddressLUT_( OpenMif("S2.TCaddr.CombinedTD.Balanced60.MixedTypes.NoSplit.mif") ),  
  cClocks_(0),
  cInputs_(0),
  cInputs2_(0),
  cInt_(0),
  cColumns_(124),
  cRows_(54),
  rOverZHistOffset_(440),
  rOverZBinSize_(64),
  kernelWidths_(),
  areaNormalizations_(),
  thresholdMaximaParam_a_(400000),
  thresholdMaximaParam_b_(2000),
  thresholdMaximaParam_c_(-80),   
  thresholdMaximaConstants_(),
  maximaWidths_{ 6 , 5 , 5 , 5 , 4 , 4 , 4 , 4 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 },
  fanoutWidths_{ 6 , 5 , 5 , 5 , 4 , 4 , 4 , 4 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 },
  cosLUT_(),
  clusterizerMagicTime_(0),
  stepLatency_({
    { UnpackLinks , 3 },        
    { TriggerCellDistribution , 2 },
    { UnpackTriggerCells , 4 },
    { TcToHc , 2 },
    { Hist   , 231 },
    { Smearing1D , 6 },
    { Interleaving , 1 },
    { NormArea   , 3 },
    { Smearing2D , 5 },
    { Deinterleaved , 4 },    
    { Maxima1D   , 5 },
    { Interleaving2 , 0 },  // Unused
    { Maxima2D   , 6 },
    { ThresholdMaxima , 3 },  // Unused?
    { CalcAverage , 4 },
    { Deinterleaving2, 0 }, // Unused
    { MaximaFanout , 9 },
    { Clusterizer , 0 },
    { TriggerCellToCluster , 8 }
    // { ClusterSum , 0 }
  }),
  depths_{ (unsigned)(-1) , 0 , 30 , 59 , 89 , 118 , 148 , 178 , 208 , 237 , 267 , 297 , 327 , 356 , 386 , 415 , 445 , 475 , 505 , 534 , 564 , 594 , 624 , 653 , 683 , 712 , 742 , 772 , 802 , 911 , 1020 , 1129 , 1238 , 1347 , 1456 , 1565 , 1674 , 1783 , 1892 , 2001 , 2110 , 2281 , 2452 , 2623 , 2794 , 2965 , 3136 , 3307 , 3478 , 3649 , 3820 },
  triggerLayers_{ (unsigned)(-1) , 1 , 0 , 2 , 0 , 3 , 0 , 4 , 0 , 5 , 0 , 6 , 0 , 7 , 0 , 8 , 0 , 9 , 0 , 10 , 0 , 11 , 0 , 12 , 0 , 13 , 0 , 14 , 0 , 15 , 16 , 17 , 18 , 19 , 20 , 21 , 22 , 23 , 24 , 25 , 26 , 27 , 28 , 29 , 30 , 31 , 32 , 33 , 34 , 35 , 36 },
  layerWeights_E_{ (unsigned)(-1) , 1 , 0 , 1 , 0 , 1 , 0 , 1 , 0 , 1 , 0 , 1 , 0 , 1 , 0 , 1 , 0 , 1 , 0 , 1 , 0 , 1 , 0 , 1 , 0 , 1 , 0 , 1 , 0 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 },
  layerWeights_E_EM_{ (unsigned)(-1) , 252969 , 0 , 254280 , 0 , 255590 , 0 , 256901 , 0 , 258212 , 0 , 259523 , 0 , 260833 , 0 , 262144 , 0 , 263455 , 0 , 264765 , 0 , 266076 , 0 , 267387 , 0 , 268698 , 0 , 270008 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
  layerWeights_E_EM_core_{ (unsigned)(-1) , 0 , 0 , 0 , 0 , 0 , 0 , 256901 , 0 , 258212 , 0 , 259523 , 0 , 260833 , 0 , 262144 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0},
  layerWeights_E_H_early_{ (unsigned)(-1) , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 1 , 1 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
  correction_( 131071 ), // 0b011111111111111111
  saturation_( (2<<19) - 1 )
{}

void ClusterAlgoConfig::setStepLatencies( const std::vector<unsigned int> latencies ) {
  // Add check that stepLatency is at least same size as latencies
  // But not as cms.exception
  for ( unsigned int iStep = 0; iStep < latencies.size(); ++iStep ) {
    stepLatency_.at(Step(iStep)) = latencies.at(iStep);
  }
}

unsigned int ClusterAlgoConfig::getLatencyUpToAndIncluding( const Step step ) const {
  unsigned int latency = 0;
  for ( int iStep = 0; iStep <= step; ++iStep ) latency += getStepLatency(Step(iStep));
  return latency;
}




std::vector<uint64_t> ClusterAlgoConfig::OpenMif( const std::string& aFilename )
{
  std::vector<uint64_t> lRet;
  std::ifstream file( aFilename );
  std::string str; 
  while ( std::getline( file , str ) ) lRet.push_back( std::stoull( str, nullptr, 16 ) );
  return lRet;
}





void ClusterAlgoConfig::initializeLUTs() {
  initializeSmearingKernelConstants( cRows_, rOverZHistOffset_, rOverZBinSize_ );
  initializeThresholdMaximaConstants( cRows_, thresholdMaximaParam_a_, thresholdMaximaParam_b_, thresholdMaximaParam_c_  );
  initializeCosLUT();
  initializeRoZToEtaLUT();
  initializeSigmaRoZToSigmaEtaLUT();
}

void ClusterAlgoConfig::initializeSmearingKernelConstants( unsigned int bins, unsigned int offset, unsigned int height ) {
  const unsigned int lWidth0 = offset + (0.5*height);
  const unsigned int lTarget = int( 6.5*lWidth0 - 0.5 ); // Magic numbers

  for ( unsigned int iBin = 0; iBin < bins; ++iBin ) {
    unsigned int lCentre = lWidth0 + ( height * iBin );
    const unsigned int lBins = int( round(1.0 * lTarget / lCentre) );

    kernelWidths_.push_back( lBins );

    lCentre *= lBins;

    const unsigned int lRatio = int( round(1.0*lTarget/lCentre * pow(2,16) ) ); // Magic numbers

    areaNormalizations_.push_back( lRatio );
  }
}

void ClusterAlgoConfig::initializeThresholdMaximaConstants( unsigned int bins, unsigned int a, unsigned int b, int c ) {
  for ( unsigned int iBin = 0; iBin < bins; ++iBin ) {
    int threshold = a + b*iBin + c*iBin*iBin;
    thresholdMaximaConstants_.push_back( threshold );
  }
}

void ClusterAlgoConfig::initializeCosLUT() {
  constexpr double Phiat120deg = 4092;
  constexpr double _120degInrad = 2*M_PI/3;
  constexpr double AngleScale = _120degInrad / Phiat120deg;
  constexpr unsigned int bitscale = pow(2,19);
  
  for ( unsigned int iBin = 0; iBin < 173; ++iBin ) { 
    unsigned int cosBin = round( bitscale * ( 1 - cos( iBin * AngleScale ) ) );
    cosLUT_.push_back( cosBin );
  }
  cosLUT_.push_back( 2047 ); // 2056 by the above formula
}

void ClusterAlgoConfig::initializeRoZToEtaLUT() {
  const float eta_min = 320 * M_PI/720; // ~1.4
  const float eta_max = 687 * M_PI/720; // ~3.0
  const float eta_LSB = M_PI/720;

  const float roz_min = 1/(sinh(eta_max)); // coming from L1T eta_max
  const float roz_max = 1/(sinh(eta_min)); // coming from L1T eta_min
  const float roz_LSB = (roz_max-roz_min)/pow(2,10); // coming from L1T eta_LSB; CP block needs to calculate mean(r/z) with at least 10b if we are to meet L1T eta LSB

  unsigned nEntries = 0x1<<10;
  rozToEtaLUT_.reserve(nEntries);
  for ( unsigned i = 0; i<nEntries; ++i ) {
    float roz = i * roz_LSB + roz_min;
    float eta = asinh(1./roz);
    unsigned eta_int_local = int( (eta - eta_min) / eta_LSB );
    unsigned eta_int_global = int( eta_int_local + eta_min / eta_LSB );
    rozToEtaLUT_.push_back(eta_int_global);
  }
}

void ClusterAlgoConfig::initializeSigmaRoZToSigmaEtaLUT() {

  // const float eta_min = 320 * M_PI/720; // ~1.4
  // const float eta_max = 687 * M_PI/720; // ~3.0
  // const float eta_LSB = M_PI/720;

  const float roz_min = 0.078953360909371;
  const float roz_max = 0.4869046060590794;
  const unsigned roz_nbits = 6;
  const float roz_LSB = (roz_max-roz_min)/pow(2,roz_nbits);

  // dEta/dMu LUT
  unsigned nEntries_dEta_dMu_LUT = 0x1<<roz_nbits;
  vector<double> dEta_dMu_LUT;
  dEta_dMu_LUT.reserve(nEntries_dEta_dMu_LUT);
  // std::cout << "dEta/dMu LUT" << std::endl;
  for ( unsigned i = 0; i<nEntries_dEta_dMu_LUT; ++i ) {
    float roz = i * roz_LSB + roz_min;
    float dEta_dMu = 1. / ( sqrt( pow(roz,-2)+1 ) * roz * roz );
    // std::cout << i << " " << roz_min << " " << roz_LSB << " " << roz << " " << dEta_dMu << std::endl;
    dEta_dMu_LUT.push_back(dEta_dMu);
  }

  // Sigma RoZ LUT
  const float sigmaRoz_min = 0;
  const float sigmaRoz_max = 0.0282874636266928;
  const unsigned sigmaRoz_nbits = 6;
  const float sigma_roz_LSB = (sigmaRoz_max-sigmaRoz_min)/pow(2,sigmaRoz_nbits);
  unsigned nEntries_sigmaRoz_LUT = 0x1<<sigmaRoz_nbits;
  vector<double> sigmaRoz_LUT;
  sigmaRoz_LUT.reserve(nEntries_sigmaRoz_LUT);
  for ( unsigned i = 0; i<nEntries_sigmaRoz_LUT; ++i ) {
    float sigma_roz = i * sigma_roz_LSB + sigmaRoz_min;
    sigmaRoz_LUT.push_back(sigma_roz);
  }

  // Final roz, sigma roz LUT
  const float sigmaEta_min = 0;
  const float sigmaEta_max = 0.1459855741275447;
  const unsigned sigmaEta_nbits = 5;
  const float sigma_eta_LSB = (sigmaEta_max-sigmaEta_min)/pow(2,sigmaEta_nbits);
  unsigned nEntries_roz_sigmaRoz_LUT = 0x1<<(roz_nbits+sigmaRoz_nbits);
  // vector<unsigned> roz_SigmaRoz_LUT(nEntries_roz_sigmaRoz_LUT);
  sigmaRozToSigmaEtaLUT_.reserve(nEntries_roz_sigmaRoz_LUT);
  for ( unsigned mu=0; mu<nEntries_dEta_dMu_LUT; ++mu) {
    double dEta_dMu = dEta_dMu_LUT.at(mu);
    for ( unsigned sigma=0; sigma<nEntries_sigmaRoz_LUT; ++sigma) {
      double sigmaRoz = sigmaRoz_LUT.at(sigma);

      float sigmaEta = dEta_dMu * sigmaRoz / sigma_eta_LSB;
      if ( sigmaEta > ( (0x1<<sigmaEta_nbits) - 1 ) ) sigmaEta = (0x1<<sigmaEta_nbits) - 1;
      sigmaRozToSigmaEtaLUT_.push_back(unsigned(sigmaEta));
      // std::cout << mu << " " << sigma << " " << dEta_dMu << " " << sigmaRoz << " " << sigma_eta_LSB << " " << dEta_dMu * sigmaRoz / sigma_eta_LSB << " " << unsigned(sigmaEta) << std::endl;
    }
  }
}

void ClusterAlgoConfig::printConfiguration() const {
  cout << "Config params" << endl;
  // cout << "Clusterizer offset : " << clusterizerOffset() << endl;
  cout << "Latencies : ";
  for ( const auto& latency : stepLatency_ ) cout << latency.first << ":" << latency.second << ", ";
  cout << endl;
  cout << "cClocks : " << cClocks() << endl;
  cout << "cInputs : " << cInputs() << endl;
  cout << "cInputs2 : " << cInputs2() << endl;
  cout << "cColumns : " << cColumns() << endl;
  cout << "cRows : " << cRows() << endl;
  cout << "rOverZHistOffset : " << rOverZHistOffset() << endl;
  cout << "rOverZBinSize : " << rOverZBinSize() << endl;
  cout << "nBinsCosLUT : " << nBinsCosLUT() << endl;
  cout << "CosLUT : ";
  for ( const auto& value : cosLUT_ ) cout << value << ", ";
  cout << endl;
  
  cout << "clusterizerMagicTime : " << clusterizerMagicTime() << endl;
  cout << "depths : ";
  for ( const auto& depth : depths() ) cout << depth << " ";
  cout << endl;
  cout << "triggerLayers : ";
  for ( const auto& triggerLayer : triggerLayers() ) cout << triggerLayer << " ";
  cout << endl;
  cout << "layerWeights_E : ";
  for ( const auto& weight : layerWeights_E() ) cout << weight << " ";
  cout << endl;
  cout << "layerWeights_E_EM : ";
  for ( const auto& weight : layerWeights_E_EM() ) cout << weight << " ";
  cout << endl;
  cout << "layerWeights_E_EM_core : ";
  for ( const auto& weight : layerWeights_E_EM_core() ) cout << weight << " ";
  cout << endl;
  cout << "layerWeights_E_H_early : ";
  for ( const auto& weight : layerWeights_E_H_early() ) cout << weight << " ";
  cout << endl;
  cout << "correction : " << correction() << endl;
  cout << "saturation : " << saturation() << endl;

  cout << "Trigger cell LUT size : " << TriggerCellDistributionLUT_.size() << endl;
}


