#ifndef __L1Trigger_L1THGCal_HGCalHistoClusteringDebugSA_h__
#define __L1Trigger_L1THGCal_HGCalHistoClusteringDebugSA_h__

#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringImpl_SA.h"

namespace l1thgcfirmware {

class HGCalHistoClusteringDebugSA : public HGCalHistoClusteringImplSA {
public:
  HGCalHistoClusteringDebugSA( ClusterAlgoConfig& config );
  virtual ~HGCalHistoClusteringDebugSA() = default;

  void runAlgorithm( const l1thgcfirmware::HGCalLinkTriggerCellSAPtrCollection& linksData );


private:

  std::size_t cHistSize;
  std::size_t cOffset;
  const bool DebugIntermediate;
  
};

}

#endif
