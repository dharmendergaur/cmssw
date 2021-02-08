#include "L1Trigger/L1THGCal/interface/backend/HGCalHistoClusteringImpl_SA.h"

#include <cmath>

#include <iostream>

HGCalHistoClusteringImplSA::HGCalHistoClusteringImplSA() {}

// void HGCalHistoClusteringImplSA::runAlgorithm() const {

// }

std::vector<l1t::HGCalMulticluster_SA> HGCalHistoClusteringImplSA::clusterSeedMulticluster_SA(
    std::vector<l1t::HGCalCluster_SA>& clusters,
    std::vector<l1t::HGCalSeed_SA>& seeds,
    std::vector<l1t::HGCalCluster_SA>& rejected_clusters,
    l1t::clusterAlgoConfig_SA& configuration) const {
  std::map<int, l1t::HGCalMulticluster_SA> mapSeedMulticluster;
  std::vector<l1t::HGCalMulticluster_SA> multiclustersOut;

  for (const auto& clu : clusters) {
    int z_side = clu.zside();

    double radiusCoefficientA = configuration.dr_byLayer_coefficientA_.empty()
                                    ? configuration.dr_
                                    : configuration.dr_byLayer_coefficientA_[clu.layer()];
    double radiusCoefficientB =
        configuration.dr_byLayer_coefficientB_.empty() ? 0 : configuration.dr_byLayer_coefficientB_[clu.layer()];

    double minDist = radiusCoefficientA + radiusCoefficientB * (configuration.kMidRadius_ - std::abs(clu.eta()));

    std::vector<std::pair<int, double>> targetSeedsEnergy;

    unsigned int iseed = 0;
    for (const auto& seed : seeds) {
      if (z_side * seed.z() < 0) {
        ++iseed;
        continue;
      }

      double seedEnergy = seed.energy();

      double d = sqrt((clu.x() - seed.x()) * (clu.x() - seed.x()) + (clu.y() - seed.y()) * (clu.y() - seed.y()));

      if (d < minDist) {
        // NearestNeighbour
        minDist = d;

        if (targetSeedsEnergy.empty()) {
          targetSeedsEnergy.emplace_back(iseed, seedEnergy);
        } else {
          targetSeedsEnergy.at(0).first = iseed;
          targetSeedsEnergy.at(0).second = seedEnergy;
        }
      }
      ++iseed;
    }

    if (targetSeedsEnergy.empty()) {
      rejected_clusters.emplace_back(clu);
      continue;
    }

    // N.B. as I have only implemented NearestNeighbour option
    // then targetSeedsEnergy has at most 1 seed for this cluster
    // Leaving in some redundant functionality in case we need
    // EnergySplit option

    //Loop over target seeds and divide up the clusters energy
    double totalTargetSeedEnergy = 0;
    for (const auto& energy : targetSeedsEnergy) {
      totalTargetSeedEnergy += energy.second;
    }

    for (const auto& energy : targetSeedsEnergy) {
      double seedWeight = 1;
      if (mapSeedMulticluster[energy.first].size() == 0) {
        mapSeedMulticluster[energy.first] = l1t::HGCalMulticluster_SA(clu, 1);
      } else {
        mapSeedMulticluster[energy.first].addConstituent(clu, true, seedWeight);
      }
    }
  }

  multiclustersOut.reserve(mapSeedMulticluster.size());
  for (const auto& mclu : mapSeedMulticluster)
    multiclustersOut.emplace_back(mclu.second);

  return multiclustersOut;
}

void HGCalHistoClusteringImplSA::finalizeClusters_SA(std::vector<l1t::HGCalMulticluster_SA>& multiclusters_in,
                                                     const std::vector<l1t::HGCalCluster_SA>& rejected_clusters_in,
                                                     std::vector<l1t::HGCalMulticluster_SA>& multiclusters_out,
                                                     std::vector<l1t::HGCalCluster_SA>& rejected_clusters_out,
                                                     l1t::clusterAlgoConfig_SA& configuration) const {
  for (const auto& tc : rejected_clusters_in) {
    rejected_clusters_out.push_back(tc);
  }

  for (auto& multicluster : multiclusters_in) {
    if (multicluster.sumPt() > configuration.ptC3dThreshold_) {
      multiclusters_out.push_back(multicluster);
    } else {
      for (const auto& tc : multicluster.constituents()) {
        rejected_clusters_out.push_back(tc);
      }
    }
  }
}