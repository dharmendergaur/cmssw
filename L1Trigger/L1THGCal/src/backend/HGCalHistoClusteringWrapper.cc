#include "L1Trigger/L1THGCal/interface/backend/HGCalHistoClusteringWrapper.h"

#include <cmath>

#include <iostream>

HGCalHistoClusteringWrapper::HGCalHistoClusteringWrapper(const edm::ParameterSet& conf)
    : dr_(conf.getParameter<double>("dR_multicluster")),
      dr_byLayer_coefficientA_(conf.existsAs<std::vector<double>>("dR_multicluster_byLayer_coefficientA")
                                   ? conf.getParameter<std::vector<double>>("dR_multicluster_byLayer_coefficientA")
                                   : std::vector<double>()),
      dr_byLayer_coefficientB_(conf.existsAs<std::vector<double>>("dR_multicluster_byLayer_coefficientB")
                                   ? conf.getParameter<std::vector<double>>("dR_multicluster_byLayer_coefficientB")
                                   : std::vector<double>()),
      ptC3dThreshold_(conf.getParameter<double>("minPt_multicluster")),
      cluster_association_input_(conf.getParameter<string>("cluster_association")),
      shape_(conf),
      theAlgo_(), 
      theConfiguration_( std::make_unique<l1t::clusterAlgoConfig_SA>( kMidRadius_, dr_, dr_byLayer_coefficientA_, dr_byLayer_coefficientB_, ptC3dThreshold_) ) {
  if (cluster_association_input_ == "NearestNeighbour") {
    cluster_association_strategy_ = NearestNeighbour;
  } else {
    throw cms::Exception("HGCTriggerParameterError")
        << "Unknown cluster association strategy'" << cluster_association_strategy_;
  }

  edm::LogInfo("HGCalMulticlusterParameters")
      << "Multicluster dR: " << dr_ << "\nMulticluster minimum transverse-momentum: " << ptC3dThreshold_;

  id_ = std::unique_ptr<HGCalTriggerClusterIdentificationBase>{
      HGCalTriggerClusterIdentificationFactory::get()->create("HGCalTriggerClusterIdentificationBDT")};
  id_->initialize(conf.getParameter<edm::ParameterSet>("EGIdentification"));
}

void HGCalHistoClusteringWrapper::convertCMSSWInputs(const std::vector<edm::Ptr<l1t::HGCalCluster>>& clustersPtrs,
                                                     std::vector<l1t::HGCalCluster_SA>& clusters_SA,
                                                     const std::vector<std::pair<GlobalPoint, double>>& seeds,
                                                     std::vector<l1t::HGCalSeed_SA>& seeds_SA) const {
  clusters_SA.clear();
  clusters_SA.reserve(clustersPtrs.size());
  unsigned int clusterIndex = 0;
  for (const auto& cluster : clustersPtrs) {
    clusters_SA.emplace_back(l1t::HGCalCluster_SA(cluster->centreProj().x(),
                                                  cluster->centreProj().y(),
                                                  triggerTools_.zside(cluster->detId()),
                                                  triggerTools_.layerWithOffset(cluster->detId()),
                                                  cluster->eta(),
                                                  cluster->phi(),
                                                  cluster->pt(),
                                                  cluster->mipPt(),
                                                  clusterIndex));
    ++clusterIndex;
  }

  seeds_SA.clear();
  seeds_SA.reserve(seeds.size());
  for (const auto& seed : seeds) {
    seeds_SA.emplace_back(l1t::HGCalSeed_SA(seed.first.x(), seed.first.y(), seed.first.z(), seed.second));
  }
}

void HGCalHistoClusteringWrapper::convertAlgorithmOutputs(
    const std::vector<l1t::HGCalMulticluster_SA>& multiclusters_out,
    const std::vector<l1t::HGCalCluster_SA>& rejected_clusters_out,
    const std::vector<edm::Ptr<l1t::HGCalCluster>>& clustersPtrs,
    l1t::HGCalMulticlusterBxCollection& multiclustersBXCollection,
    l1t::HGCalClusterBxCollection& rejected_clusters) const {
  // Not doing completely the correct thing here
  // Taking the multiclusters from the stand alone emulation
  // Getting their consistuent clusters (stand alone objects)
  // Linking back to the original CMSSW-type cluster
  // And creating a CMSSW-type multicluster based from these clusters
  // So the output multiclusters will not be storing bit accurate quantities (or whatever was derived by the stand along emulation)
  // As these inherit from L1Candidate, could set their HW quantities to the bit accurate ones
  for (const auto& rejected_cluster : rejected_clusters_out) {
    rejected_clusters.push_back(0, *clustersPtrs.at(rejected_cluster.index_cmssw()));
  }

  std::vector<l1t::HGCalMulticluster> multiclusters;
  multiclusters.reserve(multiclusters_out.size());
  for (unsigned int imulticluster = 0; imulticluster < multiclusters_out.size(); ++imulticluster) {
    bool firstConstituent = true;
    for (const auto& constituent : multiclusters_out[imulticluster].constituents()) {
      if (firstConstituent) {
        multiclusters.emplace_back(l1t::HGCalMulticluster(clustersPtrs.at(constituent.index_cmssw()), 1.));
      } else {
        multiclusters.at(imulticluster).addConstituent(clustersPtrs.at(constituent.index_cmssw()), 1.);
      }
      firstConstituent = false;
    }
  }

  for (const auto& multicluster : multiclusters) {
    multiclustersBXCollection.push_back(0, multicluster);
  }
}

void HGCalHistoClusteringWrapper::process(const std::pair< const std::vector<edm::Ptr<l1t::HGCalCluster>>, const std::vector<std::pair<GlobalPoint, double>> >& inputClustersAndSeeds, std::pair<  l1t::HGCalMulticlusterBxCollection, l1t::HGCalClusterBxCollection>& outputMulticlustersAndRejectedClusters) const {

  l1t::HGCalClusterSACollection clusters_SA;
  l1t::HGCalSeedSACollection seeds_SA;
  convertCMSSWInputs(inputClustersAndSeeds.first, clusters_SA, inputClustersAndSeeds.second, seeds_SA);

  l1t::HGCalClusterSACollection rejected_clusters_finalized_SA;
  l1t::HGCalMulticlusterSACollection multiclusters_finalized_SA;
  clusterizeHisto( clusters_SA, seeds_SA, multiclusters_finalized_SA, rejected_clusters_finalized_SA );

  convertAlgorithmOutputs(
      multiclusters_finalized_SA, rejected_clusters_finalized_SA, inputClustersAndSeeds.first, outputMulticlustersAndRejectedClusters.first, outputMulticlustersAndRejectedClusters.second);
}

void HGCalHistoClusteringWrapper::clusterizeHisto(const l1t::HGCalClusterSACollection& inputClusters,
                       const l1t::HGCalSeedSACollection& inputSeeds,
                       l1t::HGCalMulticlusterSACollection& outputMulticlusters,
                       l1t::HGCalClusterSACollection& outputRejectedClusters) const {
  // Call SA clustering
  std::vector<l1t::HGCalCluster_SA> rejected_clusters_vec_SA;
  std::vector<l1t::HGCalMulticluster_SA> multiclusters_vec_SA =
      theAlgo_.clusterSeedMulticluster_SA(inputClusters, inputSeeds, rejected_clusters_vec_SA, theConfiguration_);

  theAlgo_.finalizeClusters_SA(multiclusters_vec_SA,
                               rejected_clusters_vec_SA,
                               outputMulticlusters,
                               outputRejectedClusters,
                               theConfiguration_);
}

void HGCalHistoClusteringWrapper::configure( const std::pair<const edm::EventSetup&, const edm::ParameterSet& >& configuration ) {
  eventSetup(configuration.first);
  // theConfiguration_->setParameters( ... );
};
