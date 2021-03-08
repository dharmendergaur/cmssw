#ifndef __L1Trigger_L1THGCal_HGCalHistoClusteringImplSA_h__
#define __L1Trigger_L1THGCal_HGCalHistoClusteringImplSA_h__

#include "L1Trigger/L1THGCal/interface/backend/HGCalCluster_SA.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalSeed_SA.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalMulticluster_SA.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalHistoCluteringConfig_SA.h"

#include <string>
#include <map>
#include <vector>
#include <utility>
#include <memory>

class HGCalHistoClusteringImplSA {
public:
  HGCalHistoClusteringImplSA();
  ~HGCalHistoClusteringImplSA() {}

  void runAlgorithm() const;

  std::vector<l1t::HGCalMulticluster_SA> clusterSeedMulticluster_SA(const std::vector<l1t::HGCalCluster_SA>& clusters,
                                                                    const std::vector<l1t::HGCalSeed_SA>& seeds,
                                                                    std::vector<l1t::HGCalCluster_SA>& rejected_clusters,
                                                                    const std::unique_ptr<l1t::clusterAlgoConfig_SA>& configuration) const;

  void finalizeClusters_SA(std::vector<l1t::HGCalMulticluster_SA>&,
                           const std::vector<l1t::HGCalCluster_SA>&,
                           std::vector<l1t::HGCalMulticluster_SA>&,
                           std::vector<l1t::HGCalCluster_SA>&,
                           const std::unique_ptr<l1t::clusterAlgoConfig_SA>& configuration) const;

private:
  enum ClusterAssociationStrategy { NearestNeighbour };

  std::string cluster_association_input_;
  ClusterAssociationStrategy cluster_association_strategy_;
};

#endif
