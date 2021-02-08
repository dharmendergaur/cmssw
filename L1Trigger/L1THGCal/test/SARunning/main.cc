#include "L1Trigger/L1THGCal/interface/backend/HGCalHistoClusteringImpl_SA.h"
#include <iostream>
#include <cmath>

int main() {

	HGCalHistoClusteringImplSA theAlgo;

	// Fill input clusters and seeds
	std::vector<l1t::HGCalCluster_SA> clusters_SA = {
		#include "input_clusters.icc"
	};
	std::vector<l1t::HGCalSeed_SA> seeds_SA = {
		#include "input_seeds.icc"
	};

	// Set configuration
	std::vector<double> dr_byLayer_coefficientA_{
		#include "drCoeffA.icc"
	};

	std::vector<double> dr_byLayer_coefficientB_( dr_byLayer_coefficientA_.size(), 0  );
	l1t::clusterAlgoConfig_SA configuration( 2.3, 0.0, dr_byLayer_coefficientA_, dr_byLayer_coefficientB_, 0.5 );

	// Run the clustering step
	std::vector<l1t::HGCalCluster_SA> rejected_clusters_vec_SA;
	std::vector<l1t::HGCalMulticluster_SA> multiclusters_vec_SA = theAlgo.clusterSeedMulticluster_SA(clusters_SA, seeds_SA, rejected_clusters_vec_SA, configuration);

	// Run the finalize step
	std::vector<l1t::HGCalCluster_SA> rejected_clusters_finalized_SA;
	std::vector<l1t::HGCalMulticluster_SA> multiclusters_finalized_SA;
	theAlgo.finalizeClusters_SA(multiclusters_vec_SA, rejected_clusters_vec_SA, multiclusters_finalized_SA, rejected_clusters_finalized_SA, configuration);

	std::cout << "N multiclusters, rejected clusters : " << multiclusters_finalized_SA.size() << " " << rejected_clusters_finalized_SA.size() << std::endl;
}
