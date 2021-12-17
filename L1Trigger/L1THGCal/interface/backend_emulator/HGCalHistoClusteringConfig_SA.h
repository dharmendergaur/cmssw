#ifndef __L1Trigger_L1THGCal_HGCalHistoCluteringConfig_SA_h__
#define __L1Trigger_L1THGCal_HGCalHistoCluteringConfig_SA_h__

#include <map>
#include <vector>

namespace l1thgcfirmware {

  enum Step { Uninitialized = -1, 
              Input = 0,
              Dist0 = 1,
              Dist1 = 2,
              Dist2 = 3,
              Dist3 = 4,
              Dist4 = 5,
              Dist5 = 6,
              TcToHc = 7,
              Hist = 8,
              Smearing1D = 9,
              NormArea = 10,
              Smearing2D = 11,
              Maxima1D = 12, // Not actually used currently
              Maxima2D = 13,
              CalcAverage = 14,
              Clusterizer = 15,
              TriggerCellToCluster = 16,
              ClusterSum = 17
            };

  class ClusterAlgoConfig {
  public:
    ClusterAlgoConfig();

    void setParameters() {}

    void setSector( const unsigned int sector ) { sector_ = sector; }
    unsigned int sector() const { return sector_; }

    void setZSide( const int zSide ) { zSide_ = zSide; }
    int zSide() const { return zSide_; }

    unsigned int getStepLatency( const Step step ) const { return stepLatency_.at(step); }
    unsigned int getLatencyUpToAndIncluding( const Step step );
    unsigned int clusterizerOffset() const { return clusterizerOffset_; }
    unsigned int cClocks() const { return cClocks_; }
    unsigned int cInputs() const { return cInputs_; }
    unsigned int cInputs2() const { return cInputs2_; }
    unsigned int cInt() const { return cInt_; }
    unsigned int cColumns() const { return cColumns_; }
    unsigned int cRows() const { return cRows_; }
    unsigned int kernelWidth( unsigned int iBin ) const { return kernelWidths_.at(iBin); }
    unsigned int areaNormalization( unsigned int iBin ) const { return areaNormalizations_.at(iBin); }
    unsigned int thresholdMaxima( unsigned int iBin ) const { return thresholdMaximaConstants_.at(iBin); }
    unsigned int cosLUT( unsigned int iBin ) const { return cosLUT_.at(iBin); }
    unsigned int clusterizerMagicTime() const { return clusterizerMagicTime_; }
    unsigned int depth( unsigned int iLayer ) const { return depths_.at(iLayer); }
    unsigned int triggerLayer( unsigned int iLayer ) const { return triggerLayers_.at(iLayer); }
    unsigned int layerWeight_E( unsigned int iTriggerLayer ) const { return layerWeights_E_.at(iTriggerLayer); }
    unsigned int layerWeight_E_EM( unsigned int iTriggerLayer ) const { return layerWeights_E_EM_.at(iTriggerLayer); }
    unsigned int layerWeight_E_EM_core( unsigned int iTriggerLayer ) const { return layerWeights_E_EM_core_.at(iTriggerLayer); }
    unsigned int layerWeight_E_H_early( unsigned int iTriggerLayer ) const { return layerWeights_E_H_early_.at(iTriggerLayer); }
    unsigned int correction() const { return correction_; }
    unsigned int saturation() const { return saturation_; }

  private:
    void initializeSmearingKernelConstants( unsigned int bins, unsigned int offset, unsigned int height );
    void initializeThresholdMaximaConstants( unsigned int bins );
    void initializeCosLUT();

    unsigned int histogramOffset_;
    unsigned int clusterizerOffset_;
    unsigned int cClocks_;
    unsigned int cInputs_;
    unsigned int cInputs2_; // Better name for variable?
    unsigned int cInt_;
    unsigned int cColumns_;
    unsigned int cRows_;

    std::vector<unsigned int> kernelWidths_;
    std::vector<unsigned int> areaNormalizations_;
    std::vector<int> thresholdMaximaConstants_;

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

    unsigned int sector_;
    int zSide_;

  };

}  // namespace l1thgcfirmware

#endif
