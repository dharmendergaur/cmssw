#ifndef __L1Trigger_L1THGCal_HGCalHistoClustering_h__
#define __L1Trigger_L1THGCal_HGCalHistoClustering_h__

#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalTriggerCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistogramCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/CentroidHelper.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringConfig_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalCluster_SA.h"

#include <list>

namespace l1thgcfirmware {

  class HGCalHistoClustering {
  public:
    HGCalHistoClustering(const l1thgcfirmware::ClusterAlgoConfig& config);
    ~HGCalHistoClustering() {}

    void runClustering(l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCellsIn,
                       l1thgcfirmware::HGCalHistogramCellSAPtrCollection& histogramIn,
                       HGCalClusterSAPtrCollection& clustersOut) const;

    // Clustering
    void clusterizer( l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCells,
                      l1thgcfirmware::HGCalHistogramCellSAPtrCollection& histogram,
                      l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCellsRamOut,
                      l1thgcfirmware::HGCalHistogramCellSAPtrCollection& maximaFifoOut ) const;
    void triggerCellToCluster(const l1thgcfirmware::HGCalTriggerCellSAPtrCollection& clusteredTriggerCells,
                              const l1thgcfirmware::HGCalHistogramCellSAPtrCollection& histogram,
                              l1thgcfirmware::HGCalClusterSAPtrCollection& clustersOut) const;
    void clusterAccumulator( l1thgcfirmware::HGCalClusterSAPtrCollection& clusters, const l1thgcfirmware::HGCalHistogramCellSAPtrCollection& histogram ) const;
    void clusterTree( l1thgcfirmware::HGCalClusterSAPtrCollection& clusters ) const;

  private:
    const l1thgcfirmware::ClusterAlgoConfig& config_;
  };

}  // namespace l1thgcfirmware

#endif