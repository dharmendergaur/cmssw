#ifndef __L1Trigger_L1THGCal_HGCalHistoClusteringWrapper_h__
#define __L1Trigger_L1THGCal_HGCalHistoClusteringWrapper_h__

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "L1Trigger/L1THGCal/interface/backend/HGCalAlgoWrapperBase.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalHistoClusteringImpl_SA.h"

#include "DataFormats/L1THGCal/interface/HGCalCluster.h"
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"

#include "L1Trigger/L1THGCal/interface/backend/HGCalCluster_SA.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalSeed_SA.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalMulticluster_SA.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalShowerShape.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerTools.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalTriggerClusterIdentificationBase.h"

class HGCalHistoClusteringWrapper : public HGCalHistoClusteringWrapperBase {
public:
  HGCalHistoClusteringWrapper(const edm::ParameterSet& conf);
  ~HGCalHistoClusteringWrapper() {}

  void configure(const l1t::clusterAlgoConfig_SA& parameters) override;

  void process(const std::pair< const l1t::HGCalClusterSACollection, const l1t::HGCalSeedSACollection >& inputClustersAndSeeds, std::pair<  l1t::HGCalMulticlusterSACollection, l1t::HGCalClusterSACollection>& outputMulticlustersAndRejectedClusters) const override;

  void convertCMSSWInputs(const std::vector<edm::Ptr<l1t::HGCalCluster>>& clustersPtrs,
                          l1t::HGCalClusterSACollection& clusters_SA,
                          const std::vector<std::pair<GlobalPoint, double>>& seeds,
                          l1t::HGCalSeedSACollection& seeds_SA) const;
  void convertAlgorithmOutputs(const l1t::HGCalMulticlusterSACollection& multiclusters_out,
                               const l1t::HGCalClusterSACollection& rejected_clusters_out,
                               const std::vector<edm::Ptr<l1t::HGCalCluster>>& clustersPtrs,
                               l1t::HGCalMulticlusterBxCollection& multiclusters,
                               l1t::HGCalClusterBxCollection& rejected_clusters) const;

  void eventSetup(const edm::EventSetup& es) {
    triggerTools_.eventSetup(es);
    shape_.eventSetup(es);
    if ((!dr_byLayer_coefficientA_.empty() && (dr_byLayer_coefficientA_.size() - 1) < triggerTools_.lastLayerBH()) ||
        (!dr_byLayer_coefficientB_.empty() && (dr_byLayer_coefficientB_.size() - 1) < triggerTools_.lastLayerBH())) {
      throw cms::Exception("Configuration")
          << "The per-layer dR values go up to " << (dr_byLayer_coefficientA_.size() - 1) << "(A) and "
          << (dr_byLayer_coefficientB_.size() - 1) << "(B), while layers go up to " << triggerTools_.lastLayerBH()
          << "\n";
    }
  }

  void clusterizeHisto(const std::vector<edm::Ptr<l1t::HGCalCluster>>& clustersPtr,
                       const std::vector<std::pair<GlobalPoint, double>>& seedPositionsEnergy,
                       const HGCalTriggerGeometryBase& triggerGeometry,
                       l1t::HGCalMulticlusterBxCollection& multiclusters,
                       l1t::HGCalClusterBxCollection& rejected_clusters) const;

private:
  enum ClusterAssociationStrategy { NearestNeighbour };

  double dr_;
  std::vector<double> dr_byLayer_coefficientA_;
  std::vector<double> dr_byLayer_coefficientB_;
  double ptC3dThreshold_;

  std::string cluster_association_input_;
  ClusterAssociationStrategy cluster_association_strategy_;

  HGCalShowerShape shape_;
  HGCalTriggerTools triggerTools_;
  std::unique_ptr<HGCalTriggerClusterIdentificationBase> id_;

  HGCalHistoClusteringImplSA theAlgo_;

  std::unique_ptr<l1t::clusterAlgoConfig_SA> theConfiguration_;

  static constexpr double kMidRadius_ = 2.3;
};

#endif
