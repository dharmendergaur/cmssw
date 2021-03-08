#ifndef __L1Trigger_L1THGCal_HGCalAlgoWrapperBase_h__
#define __L1Trigger_L1THGCal_HGCalAlgoWrapperBase_h__

#include "L1Trigger/L1THGCal/interface/backend/HGCalAlgoWrapperBaseT.h"

#include "L1Trigger/L1THGCal/interface/backend/HGCalCluster_SA.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalSeed_SA.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalMulticluster_SA.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalHistoCluteringConfig_SA.h"

typedef HGCalAlgoWrapperBaseT<std::pair< const l1t::HGCalClusterSACollection, const l1t::HGCalSeedSACollection >, std::pair< l1t::HGCalMulticlusterSACollection, l1t::HGCalClusterSACollection>, l1t::clusterAlgoConfig_SA> HGCalHistoClusteringWrapperBase;

#endif
