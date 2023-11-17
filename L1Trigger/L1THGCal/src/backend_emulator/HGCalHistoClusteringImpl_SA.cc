#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringImpl_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalLinkTriggerCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalTriggerCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistogramCell_SA.h"

// #include <math.h>
// #include <map>
#include <vector>
#include <array>
#include <list>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits.h>

using namespace std;
using namespace l1thgcfirmware;

// =======================================================================================================================================================
HGCalHistoClusteringImplSA::HGCalHistoClusteringImplSA( ClusterAlgoConfig& config ) : config_(config), linkUnpacking_(config), seeding_(config), clustering_(config)
{}
// =======================================================================================================================================================

// =======================================================================================================================================================
void HGCalHistoClusteringImplSA::runAlgorithm( const l1thgcfirmware::HGCalLinkTriggerCellSAPtrCollection& linkData,
                                                l1thgcfirmware::HGCalClusterSAPtrCollection& clusterSums )
{
  std::cout << "---------------------------------------" << std::endl;
  config_.printConfiguration();
  std::cout << "---------------------------------------" << std::endl;

  // // Hack for now, feed in some simple data
  // HGCalLinkTriggerCellSAPtrCollection LinksInData;
  // for ( unsigned int i = 0; i < 54432; ++ i ) {
  //   LinksInData.emplace_back(make_unique<HGCalLinkTriggerCell>() );
  //   if ( i == 264 ) LinksInData.back()->data_ = 99;
  // }

  HGCalTriggerCellSAPtrCollection unpackedTCs;
  linkUnpacking_.runLinkUnpacking( linkData, unpackedTCs);

  for ( const auto& tc : unpackedTCs) {
    if (tc->energy() > 0 ) {
      std::cout << "TC : " << tc->energy() << std::endl;
    }
  }

  HGCalHistogramCellSAPtrCollection histogram;
  seeding_.runSeeding(unpackedTCs, histogram);

  // for ( const auto& bin : histogram ) {
  //   if ( bin->S() > 0 )
  //     std::cout << "Histo bin : " << bin->S() << std::endl;
  // }

  clustering_.runClustering(unpackedTCs, histogram, clusterSums);

  for ( const auto& cluster : clusterSums ) {
    if ( cluster->n_tc().value_ > 0 )
      std::cout << "Cluster : " << cluster->n_tc() << std::endl;
  }

  // Cluster properties  
  // clusterProperties( clusters );
}
