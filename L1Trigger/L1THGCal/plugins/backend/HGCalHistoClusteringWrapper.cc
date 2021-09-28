#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "L1Trigger/L1THGCal/interface/HGCalAlgoWrapperBase.h"

#include "DataFormats/L1THGCal/interface/HGCalCluster.h"
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"

#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringImpl_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringConfig_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalTriggerCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalCluster_SA.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerBackendDetId.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerTools.h"

class HGCalHistoClusteringWrapper : public HGCalHistoClusteringWrapperBase {
public:
  HGCalHistoClusteringWrapper(const edm::ParameterSet& conf);
  ~HGCalHistoClusteringWrapper() override {}

  void configure(const std::pair<const edm::EventSetup&, const edm::ParameterSet&>& configuration) override;

  void process(const std::vector<std::vector<edm::Ptr<l1t::HGCalCluster>>>& inputClusters,
               std::pair<l1t::HGCalMulticlusterBxCollection&, l1t::HGCalClusterBxCollection&>&
                   outputMulticlustersAndRejectedClusters) const override;

private:
  void convertCMSSWInputs(const std::vector<std::vector<edm::Ptr<l1t::HGCalCluster>>>& clustersPtrs,
                          l1thgcfirmware::HGCalTriggerCellSAPtrCollections& clusters_SA) const;
  void convertAlgorithmOutputs() const {}

  void clusterizeHisto( l1thgcfirmware::HGCalTriggerCellSAPtrCollections& triggerCells_in_SA, l1thgcfirmware::HGCalTriggerCellSAPtrCollection& clusteredTCs, l1thgcfirmware::HGCalTriggerCellSAPtrCollection& unclusteredTCs, l1thgcfirmware::CentroidHelperPtrCollection& prioritizedMaxima, l1thgcfirmware::CentroidHelperPtrCollection& readoutFlags ) const;

  void eventSetup(const edm::EventSetup& es) { triggerTools_.eventSetup(es);
                                               es.get<CaloGeometryRecord>().get("", triggerGeometry_);
                                             }

  HGCalTriggerTools triggerTools_;

  l1thgcfirmware::ClusterAlgoConfig theConfiguration_;

  l1thgcfirmware::HGCalHistoClusteringImplSA theAlgo_;

  edm::ESHandle<HGCalTriggerGeometryBase> triggerGeometry_;
};

HGCalHistoClusteringWrapper::HGCalHistoClusteringWrapper(const edm::ParameterSet& conf)
    : HGCalHistoClusteringWrapperBase(conf),
      theConfiguration_(),
      theAlgo_(theConfiguration_) {}

void HGCalHistoClusteringWrapper::convertCMSSWInputs(const std::vector<std::vector<edm::Ptr<l1t::HGCalCluster>>>& clustersPtrs, l1thgcfirmware::HGCalTriggerCellSAPtrCollections& clusters_SA ) const {

  // Convert trigger cells to format required by emulator
  l1thgcfirmware::HGCalTriggerCellSAPtrCollections clusters_SA_perSector60(3, l1thgcfirmware::HGCalTriggerCellSAPtrCollection() );
  unsigned iSector60 = 0;
  for (const auto& sector60 : clustersPtrs) {
    for (const auto& cluster : sector60) {
      const GlobalPoint& position = cluster->position();
      double x = position.x();
      double y = position.y();
      double z = position.z();
      unsigned int digi_rOverZ = ( std::sqrt(x * x + y * y) / std::abs(z) ) * 4096/ 0.7; // Magic numbers
      double phi = cluster->phi();
      phi += ( phi < 0 ) ? 2 * M_PI : 0;
      unsigned int digi_phi = ( phi ) * 1944 / M_PI; // Magic numbers
      unsigned int digi_energy = ( cluster->energy() ) * 10000; // Magic numbers
      clusters_SA_perSector60[iSector60].emplace_back( 
                                            std::make_shared<l1thgcfirmware::HGCalTriggerCell>(
                                              true,
                                              true,
                                              digi_rOverZ,
                                              digi_phi,
                                              triggerTools_.layerWithOffset(cluster->detId()),
                                              digi_energy
                                          ) );
    }
      ++iSector60;
  }

  // Distribute trigger cells to links and frames
  // Current firmware expects trigger cells from each S1 FPGA are ordered by r/z in time (r/z increase with frame number), and links from same 60 degree sector are grouped together
  // As first (optimistic) step, all trigger cells within a 60 degree sector are combined and sorted in r/z
  // Ultimately, links/ordering in time should come from S1 emulation
  // Sort by r/z in each 60 degree sector
  for (auto& clusters : clusters_SA_perSector60) {
    std::sort(clusters.begin(), clusters.end(), [](l1thgcfirmware::HGCalTriggerCellSAPtr a, l1thgcfirmware::HGCalTriggerCellSAPtr b) { return a->rOverZ()<b->rOverZ(); });
  }

  // Distribute to links
  clusters_SA.clear();
  clusters_SA.resize( 212, l1thgcfirmware::HGCalTriggerCellSAPtrCollection() ); // Magic numbers
  for ( auto& clusters : clusters_SA ) {
    clusters.resize(72, l1thgcfirmware::HGCalTriggerCellSAPtr( std::make_shared<l1thgcfirmware::HGCalTriggerCell>() ) ); // Magic numbers
  }
  iSector60 = 0;
  for (const auto& sector60 : clusters_SA_perSector60) {
    unsigned iCluster = 0;
    for ( const auto& cluster : sector60 ) {
      // Leave first two frames empty
      unsigned frame = 2 + iCluster / 24; // Magic numbers
      unsigned link = iCluster % 24 + iSector60 * 24; // Magic numbers
      if ( frame >= 212 ) break;
      clusters_SA[frame][link] = cluster;
      ++iCluster;
    }
    ++iSector60;
  }
}

void HGCalHistoClusteringWrapper::process(
    const std::vector<std::vector<edm::Ptr<l1t::HGCalCluster>>>&
        inputClusters,
    std::pair<l1t::HGCalMulticlusterBxCollection&, l1t::HGCalClusterBxCollection&>&
        outputMulticlustersAndRejectedClusters) const {

  l1thgcfirmware::HGCalTriggerCellSAPtrCollections triggerCells_in_SA;
  convertCMSSWInputs(inputClusters, triggerCells_in_SA);

  l1thgcfirmware::HGCalTriggerCellSAPtrCollection clusteredTCs_out_SA;
  l1thgcfirmware::HGCalTriggerCellSAPtrCollection unclusteredTCs_out_SA;
  l1thgcfirmware::CentroidHelperPtrCollection prioritizedMaxima_out_SA;
  l1thgcfirmware::CentroidHelperPtrCollection readoutFlags_out_SA;
  clusterizeHisto(triggerCells_in_SA, clusteredTCs_out_SA, unclusteredTCs_out_SA, prioritizedMaxima_out_SA, readoutFlags_out_SA);

  // convertAlgorithmOutputs();
}

void HGCalHistoClusteringWrapper::clusterizeHisto( l1thgcfirmware::HGCalTriggerCellSAPtrCollections& triggerCells_in_SA, l1thgcfirmware::HGCalTriggerCellSAPtrCollection& clusteredTCs, l1thgcfirmware::HGCalTriggerCellSAPtrCollection& unclusteredTCs, l1thgcfirmware::CentroidHelperPtrCollection& prioritizedMaxima, l1thgcfirmware::CentroidHelperPtrCollection& readoutFlags ) const {

  theAlgo_.runAlgorithm( triggerCells_in_SA, clusteredTCs, unclusteredTCs, prioritizedMaxima, readoutFlags );
}

void HGCalHistoClusteringWrapper::configure(
    const std::pair<const edm::EventSetup&, const edm::ParameterSet&>& configuration) {
  eventSetup(configuration.first);

  // theConfiguration_.setParameters( ... );
};

DEFINE_EDM_PLUGIN(HGCalHistoClusteringWrapperBaseFactory, HGCalHistoClusteringWrapper, "HGCalHistoClusteringWrapper");