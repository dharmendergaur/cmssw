#ifndef __L1Trigger_L1THGCal_HGCalHistoClusteringImplSA_h__
#define __L1Trigger_L1THGCal_HGCalHistoClusteringImplSA_h__

#include "DataFormats/L1THGCal/interface/HGCalCluster_SA.h"
#include "DataFormats/L1THGCal/interface/HGCalSeed_SA.h"
#include "DataFormats/L1THGCal/interface/HGCalMulticluster_SA.h"

#include <string>
#include <map>
#include <vector>
#include <utility>

namespace l1t {

struct clusterAlgoConfig_SA {
  clusterAlgoConfig_SA( double kMidRadius,
                        double dr, 
                        std::vector<double> dr_byLayer_coefficientA, 
                        std::vector<double> dr_byLayer_coefficientB,
                        float ptC3dThreshold) :
                            kMidRadius_(kMidRadius),
                            dr_(dr),
                            dr_byLayer_coefficientA_(dr_byLayer_coefficientA), 
                            dr_byLayer_coefficientB_(dr_byLayer_coefficientB),
                            ptC3dThreshold_(ptC3dThreshold)
                             {}
  const double kMidRadius_;
  const double dr_;
  const std::vector<double> dr_byLayer_coefficientA_;
  const std::vector<double> dr_byLayer_coefficientB_;
  const float ptC3dThreshold_;
};

}

class HGCalHistoClusteringImplSA {
public:
  HGCalHistoClusteringImplSA();
  ~HGCalHistoClusteringImplSA() {}

  void runAlgorithm() const;

  std::vector<l1t::HGCalMulticluster_SA> clusterSeedMulticluster_SA(std::vector<l1t::HGCalCluster_SA>& clusters,
                                  std::vector<l1t::HGCalSeed_SA>& seeds,
                                  std::vector<l1t::HGCalCluster_SA>& rejected_clusters,
                                  l1t::clusterAlgoConfig_SA& configuration) const;

  void finalizeClusters_SA(std::vector<l1t::HGCalMulticluster_SA>&,
                        const std::vector<l1t::HGCalCluster_SA>&,
                        std::vector<l1t::HGCalMulticluster_SA>&,
                        std::vector<l1t::HGCalCluster_SA>&,
                        l1t::clusterAlgoConfig_SA& configuration) const;
                        
private:
  enum ClusterAssociationStrategy { NearestNeighbour };



  std::string cluster_association_input_;
  ClusterAssociationStrategy cluster_association_strategy_;

};

#endif
