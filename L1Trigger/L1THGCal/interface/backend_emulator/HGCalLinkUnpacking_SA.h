#ifndef __L1Trigger_L1THGCal_HGCalLinkUnpacking_h__
#define __L1Trigger_L1THGCal_HGCalLinkUnpacking_h__

#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalTriggerCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalLinkTriggerCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringConfig_SA.h"

#include <vector>

namespace l1thgcfirmware {

  class HGCalLinkUnpacking {
  public:
    HGCalLinkUnpacking(const l1thgcfirmware::ClusterAlgoConfig& config);
    ~HGCalLinkUnpacking() {}

    void runLinkUnpacking(const l1thgcfirmware::HGCalLinkTriggerCellSAPtrCollection& linksIn,
                                    l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCellsOut) const;

    void unpackLinks( const l1thgcfirmware::HGCalLinkTriggerCellSAPtrCollection& LinksIn , l1thgcfirmware::HGCalLinkTriggerCellSAPtrCollection& TriggerCells , l1thgcfirmware::HGCalLinkTriggerCellSAPtrCollection& Towers ) const;
    l1thgcfirmware::HGCalTriggerCellSAPtrCollection triggerCellDistribution( const l1thgcfirmware::HGCalLinkTriggerCellSAPtrCollection& TriggerCells ) const;
    void unpackTriggerCells( const l1thgcfirmware::HGCalTriggerCellSAPtrCollection& TriggerCells ) const;

  private:

    const l1thgcfirmware::ClusterAlgoConfig& config_;
  };
}  // namespace l1thgcfirmware

#endif