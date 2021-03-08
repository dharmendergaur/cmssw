#ifndef __L1Trigger_L1THGCal_HGCalHistoClusteringWrapper_h__
#define __L1Trigger_L1THGCal_HGCalHistoClusteringWrapper_h__

#include "FWCore/ParameterSet/interface/ParameterSet.h"

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

// EJC Could this class inherit from a generic wrapper class?  Would we need "convert" functions to have same number of input collections for each algorithm?
class HGCalHistoClusteringWrapper {
public:
  HGCalHistoClusteringWrapper(const edm::ParameterSet& conf);
  ~HGCalHistoClusteringWrapper() {}

  void convertCMSSWInputs(const std::vector<edm::Ptr<l1t::HGCalCluster>>& clustersPtrs,
                          std::vector<l1t::HGCalCluster_SA>& clusters_SA,
                          const std::vector<std::pair<GlobalPoint, double>>& seeds,
                          std::vector<l1t::HGCalSeed_SA>& seeds_SA) const;
  void convertAlgorithmOutputs(const std::vector<l1t::HGCalMulticluster_SA>& multiclusters_out,
                               const std::vector<l1t::HGCalCluster_SA>& rejected_clusters_out,
                               const std::vector<edm::Ptr<l1t::HGCalCluster>>& clustersPtrs,
                               l1t::HGCalMulticlusterBxCollection& multiclusters,
                               l1t::HGCalClusterBxCollection& rejected_clusters) const;

  void runAlgorithm() const;

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

  static constexpr double kMidRadius_ = 2.3;
};

#endif
