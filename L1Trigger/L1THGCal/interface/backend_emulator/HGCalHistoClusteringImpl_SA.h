#ifndef __L1Trigger_L1THGCal_HGCalHistoClusteringImplSA_h__
#define __L1Trigger_L1THGCal_HGCalHistoClusteringImplSA_h__

#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalLinkTriggerCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalTriggerCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistogramCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalCluster_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringConfig_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalLinkUnpacking_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoSeeding_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClustering_SA.h"

#include <vector>
#include <memory>

namespace l1thgcfirmware {

class HGCalHistoClusteringImplSA {
public:
  HGCalHistoClusteringImplSA(ClusterAlgoConfig& config);
  virtual ~HGCalHistoClusteringImplSA(){};

  virtual void runAlgorithm( const l1thgcfirmware::HGCalLinkTriggerCellSAPtrCollection& linksData,                  
                              l1thgcfirmware::HGCalClusterSAPtrCollection& clusterSums );

protected:
  const ClusterAlgoConfig& config_;
  l1thgcfirmware::HGCalLinkUnpacking linkUnpacking_;
  l1thgcfirmware::HGCalHistoSeeding seeding_;
  l1thgcfirmware::HGCalHistoClustering clustering_;

};

}

#endif
