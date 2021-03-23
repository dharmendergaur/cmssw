#include "L1Trigger/L1THGCal/interface/HGCalProcessorBase.h"

#include "DataFormats/L1THGCal/interface/HGCalTriggerCell.h"
#include "DataFormats/L1THGCal/interface/HGCalCluster.h"
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalMulticlusteringImpl.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalHistoSeedingImpl.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalHistoClusteringImpl.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalHistoClusteringWrapper.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalTriggerClusterInterpreterBase.h"

#include <utility>

// Can never get MessageLogger to do what I want...
#include <iostream>

class HGCalBackendLayer2Processor3DClusteringSA : public HGCalBackendLayer2ProcessorBase {
public:
  HGCalBackendLayer2Processor3DClusteringSA(const edm::ParameterSet& conf) : HGCalBackendLayer2ProcessorBase(conf), conf_(conf) {
    std::cout << "In stand alone Clustering Processor" << std::endl;

    std::string typeMulticluster(conf.getParameterSet("C3d_parameters").getParameter<std::string>("type_multicluster"));
    if (typeMulticluster == "Histo") {
      multiclusteringAlgoType_ = HistoC3d;
      multiclusteringHistoSeeding_ = std::make_unique<HGCalHistoSeedingImpl>(
          conf.getParameterSet("C3d_parameters").getParameterSet("histoMax_C3d_seeding_parameters"));
      multiclusteringHistoClustering_ = std::make_unique<HGCalHistoClusteringImpl>(
          conf.getParameterSet("C3d_parameters").getParameterSet("histoMax_C3d_clustering_parameters"));
    } else if (typeMulticluster == "SAHisto") {
      multiclusteringAlgoType_ = SAHistoC3d;
      multiclusteringHistoSeeding_ = std::make_unique<HGCalHistoSeedingImpl>(
          conf.getParameterSet("C3d_parameters").getParameterSet("histoMax_C3d_seeding_parameters"));
      multiclusteringHistoClusteringWrapper_ = std::make_unique<HGCalHistoClusteringWrapper>(
          conf.getParameterSet("C3d_parameters").getParameterSet("histoMax_C3d_clustering_parameters"));
    } else {
      throw cms::Exception("HGCTriggerParameterError") << "Unknown Multiclustering type '" << typeMulticluster << "'";
    }

    for (const auto& interpretationPset : conf.getParameter<std::vector<edm::ParameterSet>>("energy_interpretations")) {
      std::unique_ptr<HGCalTriggerClusterInterpreterBase> interpreter{
          HGCalTriggerClusterInterpreterFactory::get()->create(interpretationPset.getParameter<std::string>("type"))};
      interpreter->initialize(interpretationPset);
      energy_interpreters_.push_back(std::move(interpreter));
    }
  }

  void run(const edm::Handle<l1t::HGCalClusterBxCollection>& collHandle,
           std::pair<l1t::HGCalMulticlusterBxCollection, l1t::HGCalClusterBxCollection>& be_output,
           const edm::EventSetup& es) override {
    es.get<CaloGeometryRecord>().get("", triggerGeometry_);
    if (multiclusteringHistoSeeding_)
      multiclusteringHistoSeeding_->eventSetup(es);
    if (multiclusteringHistoClustering_)
      multiclusteringHistoClustering_->eventSetup(es);

    auto& collCluster3D = be_output.first;
    auto& rejectedClusters = be_output.second;

    /* create a persistent vector of pointers to the trigger-cells */
    std::vector<edm::Ptr<l1t::HGCalCluster>> clustersPtrs;
    for (unsigned i = 0; i < collHandle->size(); ++i) {
      edm::Ptr<l1t::HGCalCluster> ptr(collHandle, i);
      clustersPtrs.push_back(ptr);
    }

    /* create a vector of seed positions and their energy*/
    std::vector<std::pair<GlobalPoint, double>> seedPositionsEnergy;

    /* call to multiclustering and compute shower shape*/
    switch (multiclusteringAlgoType_) {
      case HistoC3d:
        multiclusteringHistoSeeding_->findHistoSeeds(clustersPtrs, seedPositionsEnergy);
        multiclusteringHistoClustering_->clusterizeHisto(
            clustersPtrs, seedPositionsEnergy, *triggerGeometry_, collCluster3D, rejectedClusters);
        break;
      case SAHistoC3d:
        {
          multiclusteringHistoSeeding_->findHistoSeeds(clustersPtrs, seedPositionsEnergy);
          
          // Inputs
          std::pair< const std::vector<edm::Ptr<l1t::HGCalCluster>>, const std::vector<std::pair<GlobalPoint, double > > > inputClustersAndSeeds{ clustersPtrs, seedPositionsEnergy };
          // Outputs
          std::pair< l1t::HGCalMulticlusterBxCollection, l1t::HGCalClusterBxCollection > outputMulticlustersAndRejectedClusters{ collCluster3D, rejectedClusters };
          // Configuration
          const std::pair<const edm::EventSetup&, const edm::ParameterSet& > configuration{ es, conf_ };

          // Configure and process
          multiclusteringHistoClusteringWrapper_->configure( configuration );
          multiclusteringHistoClusteringWrapper_->process( inputClustersAndSeeds, outputMulticlustersAndRejectedClusters );
        }
        break;
      default:
        // Should not happen, clustering type checked in constructor
        break;
    }

    // Call all the energy interpretation modules on the cluster collection
    for (const auto& interpreter : energy_interpreters_) {
      interpreter->eventSetup(es);
      interpreter->interpret(collCluster3D);
    }
  }

private:
  enum MulticlusterType { HistoC3d, SAHistoC3d };

  edm::ESHandle<HGCalTriggerGeometryBase> triggerGeometry_;

  /* algorithms instances */
  std::unique_ptr<HGCalHistoSeedingImpl> multiclusteringHistoSeeding_;
  std::unique_ptr<HGCalHistoClusteringImpl> multiclusteringHistoClustering_;
  std::unique_ptr<HGCalHistoClusteringWrapper> multiclusteringHistoClusteringWrapper_;

  /* algorithm type */
  MulticlusterType multiclusteringAlgoType_;

  std::vector<std::unique_ptr<HGCalTriggerClusterInterpreterBase>> energy_interpreters_;

  const edm::ParameterSet conf_;
};

DEFINE_EDM_PLUGIN(HGCalBackendLayer2Factory,
                  HGCalBackendLayer2Processor3DClusteringSA,
                  "HGCalBackendLayer2Processor3DClusteringSA");
