#ifndef __L1Trigger_L1THGCal_HGCalHistoCluteringConfig_SA_h__
#define __L1Trigger_L1THGCal_HGCalHistoCluteringConfig_SA_h__

#include <map>
#include <vector>
#include <cstdint>
#include <string>
#include <math.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif



namespace l1thgcfirmware {

  enum Step { Uninitialized = -1, 
              UnpackLinks = 0,              
              TriggerCellDistribution = 1,
              UnpackTriggerCells = 2,
              TcToHc = 3,
              Hist = 4,
              Smearing1D = 5,
              Interleaving = 6,
              NormArea = 7,
              Smearing2D = 8,
              Deinterleaved = 9,
              Maxima1D = 10,
              Interleaving2 = 11,              
              Maxima2D = 12,
              ThresholdMaxima = 13,
              CalcAverage = 14,
              Deinterleaving2 = 15,
              MaximaFanout = 16,
              Clusterizer = 17,
              TriggerCellToCluster = 18
            };

  class ClusterAlgoConfig {
  public:
    ClusterAlgoConfig();
    ~ClusterAlgoConfig() {};

    void setParameters() {}

    void setSector( const unsigned int sector ) { sector_ = sector; }
    unsigned int sector() const { return sector_; }

    void setZSide( const int zSide ) { zSide_ = zSide; }
    int zSide() const { return zSide_; }

    void setStepLatencies( const std::vector<unsigned int> latencies );
    unsigned int getStepLatency( const Step step ) const { return stepLatency_.at(step); }
    unsigned int getLatencyUpToAndIncluding( const Step step ) const;

    // void setClusterizerOffset( const unsigned clusterizerOffset ) { clusterizerOffset_ = clusterizerOffset ;}
    // unsigned int clusterizerOffset() const { return clusterizerOffset_; }

    void setCClocks( const unsigned cClocks ) { cClocks_ = cClocks;}
    unsigned int cClocks() const { return cClocks_; }

    void setCInputs( const unsigned cInputs ) { cInputs_ = cInputs;}
    unsigned int cInputs() const { return cInputs_; }

    void setCInputs2( const unsigned cInputs2 ) { cInputs2_ = cInputs2;}
    unsigned int cInputs2() const { return cInputs2_; }

    void setCInt( const unsigned cInt ) { cInt_ = cInt;}
    unsigned int cInt() const { return cInt_; }

    void setCColumns( const unsigned cColumns ) { cColumns_ = cColumns;}
    unsigned int cColumns() const { return cColumns_; }

    void setCRows( const unsigned cRows ) { cRows_ = cRows;}
    unsigned int cRows() const { return cRows_; }

    void setROverZHistOffset( const unsigned rOverZHistOffset ) { rOverZHistOffset_ = rOverZHistOffset;}
    unsigned int rOverZHistOffset() const { return rOverZHistOffset_; }

    void setROverZBinSize( const unsigned rOverZBinSize ) { rOverZBinSize_ = rOverZBinSize;}
    unsigned int rOverZBinSize() const { return rOverZBinSize_; }

    unsigned int kernelWidth( unsigned int iBin ) const { return kernelWidths_.at(iBin); }
    unsigned int areaNormalization( unsigned int iBin ) const { return areaNormalizations_.at(iBin); }
    
    void setROverZRange( const float rOverZRange ) { rOverZRange_ = rOverZRange_; }
    unsigned int rOverZRange() const { return rOverZRange_; }
    void setROverZNValues( const float rOverZNValues ) { rOverZNValues_ = rOverZNValues_; }
    unsigned int rOverZNValues() const { return rOverZNValues_; }
    void setPhiRange( const float phiRange ) { phiRange_ = phiRange_; }
    unsigned int phiRange() const { return phiRange_; }
    void setPhiNValues( const float phiNValues ) { phiNValues_ = phiNValues_; }
    unsigned int phiNValues() const { return phiNValues_; }
    void setPtDigiFactor(const float ptDigiFactor) { ptDigiFactor_ = ptDigiFactor; }
    float ptDigiFactor() const { return ptDigiFactor_; }
    void setMinClusterPtOut(const float pt) { minClusterPtOut_ = pt; }
    float minClusterPtOut() const { return minClusterPtOut_; }

    void setMaxClustersPerLink(const unsigned maxClustersPerLink) { maxClustersPerLink_ = maxClustersPerLink; }
    unsigned maxClustersPerLink() const { return maxClustersPerLink_; }
    void setNInputLinks(const unsigned nInputLinks) { nInputLinks_ = nInputLinks; }
    unsigned nInputLinks() const { return nInputLinks_; }

    void setThresholdParams( const unsigned int a, const unsigned int b, const int c ) { 
      thresholdMaximaParam_a_ = a;
      thresholdMaximaParam_b_ = b;
      thresholdMaximaParam_c_ = c;
    }
    unsigned int thresholdMaxima( unsigned int iBin ) const { return thresholdMaximaConstants_.at(iBin); }

    unsigned int maximaWidths( unsigned int iBin ) const { return maximaWidths_.at(iBin); }

    unsigned int fanoutWidths( unsigned int iBin ) const { return fanoutWidths_.at(iBin); }


    unsigned int nBinsCosLUT() const { return cosLUT_.size(); }

    unsigned int cosLUT( unsigned int iBin ) const { return cosLUT_.at(iBin); }

    void setClusterizerMagicTime( const unsigned clusterizerMagicTime ) { clusterizerMagicTime_ = clusterizerMagicTime;}
    unsigned int clusterizerMagicTime() const { return clusterizerMagicTime_; }

    void setDepths( const std::vector<unsigned int> depths ) {
      depths_.clear();
      for ( const auto& depth : depths ) depths_.push_back( depth );
    }
    std::vector<unsigned int> depths() const { return depths_; }
    unsigned int depth( unsigned int iLayer ) const { return depths_.at(iLayer); }

    void setTriggerLayers( const std::vector<unsigned int> triggerLayers ) {
      triggerLayers_.clear();
      for ( const auto& triggerLayer : triggerLayers ) triggerLayers_.push_back( triggerLayer );
    }
    std::vector<unsigned int> triggerLayers() const { return triggerLayers_; }
    unsigned int triggerLayer( unsigned int iLayer ) const { return triggerLayers_.at(iLayer); }
    
    void setNTriggerLayers(const unsigned n) { nTriggerLayers_ = n; }
    unsigned int nTriggerLayers() const { return nTriggerLayers_; }

    void setLayerWeights_E( const std::vector<unsigned int> layerWeights_E ) {
      layerWeights_E_.clear();
      for ( const auto& weight : layerWeights_E ) layerWeights_E_.push_back( weight );
    }
    std::vector<unsigned int> layerWeights_E() const { return layerWeights_E_; }
    unsigned int layerWeight_E( unsigned int iTriggerLayer ) const { return layerWeights_E_.at(iTriggerLayer); }

    void setLayerWeights_E_EM( const std::vector<unsigned int> layerWeights_E_EM ) {
      layerWeights_E_EM_.clear();
      for ( const auto& weight : layerWeights_E_EM ) layerWeights_E_EM_.push_back( weight );
    }
    std::vector<unsigned int> layerWeights_E_EM() const { return layerWeights_E_EM_; }
    unsigned int layerWeight_E_EM( unsigned int iTriggerLayer ) const { return layerWeights_E_EM_.at(iTriggerLayer); }
    
    void setLayerWeights_E_EM_core( const std::vector<unsigned int> layerWeights_E_EM_core ) {
      layerWeights_E_EM_core_.clear();
      for ( const auto& weight : layerWeights_E_EM_core ) layerWeights_E_EM_core_.push_back( weight );
    }
    std::vector<unsigned int> layerWeights_E_EM_core() const { return layerWeights_E_EM_core_; }
    unsigned int layerWeight_E_EM_core( unsigned int iTriggerLayer ) const { return layerWeights_E_EM_core_.at(iTriggerLayer); }


    void setLayerWeights_E_H_early( const std::vector<unsigned int> layerWeights_E_H_early ) {
      layerWeights_E_H_early_.clear();
      for ( const auto& weight : layerWeights_E_H_early ) layerWeights_E_H_early_.push_back( weight );
    }
    std::vector<unsigned int> layerWeights_E_H_early() const { return layerWeights_E_H_early_; }
    unsigned int layerWeight_E_H_early( unsigned int iTriggerLayer ) const { return layerWeights_E_H_early_.at(iTriggerLayer); }

    void setCorrection( const unsigned correction ) { correction_ = correction;}
    unsigned int correction() const { return correction_; }

    void setSaturation( const unsigned saturation ) { saturation_ = saturation;}
    unsigned int saturation() const { return saturation_; }

    unsigned int rozToEtaLUT(unsigned int iBin) const { return rozToEtaLUT_.at(iBin); }
    unsigned int sigmaRozToSigmaEtaLUT(unsigned int iBin) const { return sigmaRozToSigmaEtaLUT_.at(iBin); }

    const uint64_t& TriggerCellDistributionLUT( unsigned int addr ) const { return TriggerCellDistributionLUT_.at(addr); }
    const uint64_t& TriggerCellAddressLUT( unsigned int addr ) const { return TriggerCellAddressLUT_.at(addr); }

    void initializeLUTs();

    void printConfiguration() const; // For debugging

  private:
    void initializeSmearingKernelConstants( unsigned int bins, unsigned int offset, unsigned int height );
    void initializeThresholdMaximaConstants( unsigned int bins, unsigned int a, unsigned int b, int c );
    void initializeCosLUT();
    void initializeRoZToEtaLUT();
    void initializeSigmaRoZToSigmaEtaLUT();
    static std::vector<uint64_t> OpenMif( const std::string& aFilename );

    std::vector<uint64_t> TriggerCellDistributionLUT_;
    std::vector<uint64_t> TriggerCellAddressLUT_;
    
    
    // unsigned int histogramOffset_;
    // unsigned int clusterizerOffset_;
    unsigned int cClocks_;
    unsigned int cInputs_;
    unsigned int cInputs2_; // Better name for variable?
    unsigned int cInt_;
    unsigned int cColumns_;
    unsigned int cRows_;
    unsigned int rOverZHistOffset_;
    unsigned int rOverZBinSize_;

    float rOverZRange_;
    unsigned rOverZNValues_;
    float phiRange_;
    unsigned phiNValues_;
    float ptDigiFactor_;

    // Selection on output clusters
    float minClusterPtOut_;

    // Input link params
    unsigned int maxClustersPerLink_;
    unsigned int nInputLinks_;

    std::vector<unsigned int> kernelWidths_;
    std::vector<unsigned int> areaNormalizations_;

    unsigned int thresholdMaximaParam_a_;
    unsigned int thresholdMaximaParam_b_;
    int thresholdMaximaParam_c_;
    std::vector<int> thresholdMaximaConstants_;
    std::vector<int> maximaWidths_;
    std::vector<int> fanoutWidths_;

    std::vector<unsigned int> cosLUT_;

    unsigned int clusterizerMagicTime_;

    std::map<Step,unsigned int> stepLatency_;

    // Parameters for triggerCellToCluster
    std::vector<unsigned int> depths_;
    std::vector<unsigned int> triggerLayers_;
    std::vector<unsigned int> layerWeights_E_;
    std::vector<unsigned int> layerWeights_E_EM_;
    std::vector<unsigned int> layerWeights_E_EM_core_;
    std::vector<unsigned int> layerWeights_E_H_early_;
    unsigned int correction_;
    unsigned int saturation_;

    std::vector<unsigned int> rozToEtaLUT_;
    std::vector<unsigned int> sigmaRozToSigmaEtaLUT_;

    unsigned int nTriggerLayers_;

    unsigned int sector_;
    int zSide_;

  };

}  // namespace l1thgcfirmware

#endif
