#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "L1Trigger/L1THGCal/interface/HGCalAlgoWrapperBase.h"

#include "DataFormats/L1THGCal/interface/HGCalCluster.h"
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"

#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringImpl_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringConfig_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalLinkTriggerCell_SA.h"
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

  void configure(
      const std::tuple<const HGCalTriggerGeometryBase* const, const edm::ParameterSet&, const unsigned int, const int>&
          configuration) override;

  void process(const std::map<unsigned int, std::vector<edm::Ptr<l1t::HGCalCluster>>>& inputClusters,
               std::pair<l1t::HGCalMulticlusterBxCollection&, l1t::HGCalClusterBxCollection&>&
                   outputMulticlustersAndRejectedClusters) override;

private:
<<<<<<< HEAD
  void convertCMSSWInputs(const std::vector<std::vector<edm::Ptr<l1t::HGCalCluster>>>& clustersPtrs,
                          l1thgcfirmware::HGCalLinkTriggerCellSAPtrCollection& linkData_SA) const;
=======
  void convertCMSSWInputs(const std::map<unsigned int, std::vector<edm::Ptr<l1t::HGCalCluster>>>& clustersPtrs,
                          l1thgcfirmware::HGCalTriggerCellSAPtrCollections& clusters_SA) const;
>>>>>>> 2a520b2a914e (Group TCs per module, sort by phi)
  void convertAlgorithmOutputs(l1thgcfirmware::HGCalClusterSAPtrCollection& clusterSums,
                               l1t::HGCalMulticlusterBxCollection& multiClusters_out,
                               const std::map<unsigned int, std::vector<edm::Ptr<l1t::HGCalCluster>>>& inputClustersPtrs) const;

  void clusterizeHisto( const l1thgcfirmware::HGCalLinkTriggerCellSAPtrCollection& linkData_in_SA,
                       l1thgcfirmware::HGCalClusterSAPtrCollection& clusterSums);

  void setGeometry(const HGCalTriggerGeometryBase* const geom) { triggerTools_.setGeometry(geom); }

  double rotatePhiToSectorZero(const double phi, const unsigned sector) const;
  double rotatePhiFromSectorZero(const double phi, const unsigned sector) const;

  HGCalTriggerTools triggerTools_;

  l1thgcfirmware::ClusterAlgoConfig theConfiguration_;

  l1thgcfirmware::HGCalHistoClusteringImplSA theAlgo_;

  edm::ESHandle<HGCalTriggerGeometryBase> triggerGeometry_;
};

HGCalHistoClusteringWrapper::HGCalHistoClusteringWrapper(const edm::ParameterSet& conf)
    : HGCalHistoClusteringWrapperBase(conf), theConfiguration_(), theAlgo_(theConfiguration_) {}

void HGCalHistoClusteringWrapper::convertCMSSWInputs(
    const std::vector<std::vector<edm::Ptr<l1t::HGCalCluster>>>& clustersPtrs,
    l1thgcfirmware::HGCalLinkTriggerCellSAPtrCollection& linkData_in) const {
  // std::cout << "N TCs : " << clustersPtrs.size() << " " << clustersPtrs.at(0).size() << std::endl;
  // for (const auto& sector60 : clustersPtrs) {
  //   for (const auto& cluster : sector60) {
  //     linkData_in.emplace_back(std::make_unique<l1thgcfirmware::HGCalLinkTriggerCell>());
  //     linkData_in.back()->data_.value_ = 99;
  //     std::cout << "Added a link data : " << linkData_in.size() << " " << linkData_in.back()->data_.value_ << std::endl;
  //     break;
  //   }
  //   break;
  // }

  // for (const auto& linkData : linkData_in ) {
  //   if (linkData->data_.value_ != 0 ) std::cout << "Got input link data : " << linkData->data_ << std::endl;
  // }
    // const std::map<unsigned int, std::vector<edm::Ptr<l1t::HGCalCluster>>>& clustersPtrs,
    // l1thgcfirmware::HGCalTriggerCellSAPtrCollections& clusters_SA) const {
  // Convert trigger cells to format required by emulator

  const int numLinks = 84;
  const int numFrames = 162;

  std::vector<std::vector<int64_t>> LinkData(numLinks, std::vector<int64_t>(numFrames, 0));

  l1thgcfirmware::HGCalTriggerCellSAPtrCollections clusters_SA_permodule(clustersPtrs.size());
  for (const auto& module : clustersPtrs) {
    unsigned iCluster = 0;
    std::cout << module.first << std::endl;
    for (const auto& cluster : module.second) {
      const GlobalPoint& position = cluster->position();
      double x = position.x();
      double y = position.y();
      double z = position.z();
      unsigned int digi_rOverZ = (std::sqrt(x * x + y * y) / std::abs(z)) * theConfiguration_.rOverZNValues() /
                                 theConfiguration_.rOverZRange();

      if (z > 0)
        x = -x;
      double phi = std::atan2(y, x);
      // Rotate phi to sector 0
      auto sector = theConfiguration_.sector();
      phi = rotatePhiToSectorZero(phi, sector);

      // Ignore TCs that are outside of the nominal 180 degree S2 sector
      // Assume these cannot be part of a cluster found within the central 120 degrees of the S2 sector?
      if (phi < 0 || phi > M_PI) {
        continue;
      }

      unsigned int digi_phi = (phi)*theConfiguration_.phiNValues() / theConfiguration_.phiRange();
      unsigned int digi_energy = (cluster->pt()) * theConfiguration_.ptDigiFactor();

      clusters_SA_permodule[iCluster].emplace_back(std::make_unique<l1thgcfirmware::HGCalTriggerCell>(
          true, true, digi_rOverZ, digi_phi, triggerTools_.layerWithOffset(cluster->detId()), digi_energy));
      clusters_SA_permodule[iCluster].back()->setCmsswIndex(std::pair<int, int>{iCluster, module.first});
      ++iCluster;
    }
  }

  // Distribute to links
  clusters_SA.clear();
  const unsigned empty_frames = 2;
  clusters_SA.resize(theConfiguration_.maxClustersPerLink() + empty_frames);
  for (auto& clusters : clusters_SA) {
    for (unsigned int iCluster = 0; iCluster < theConfiguration_.nInputLinks(); ++iCluster) {
      clusters.push_back(std::make_unique<l1thgcfirmware::HGCalTriggerCell>());
    }
  }
  // iSector60 = 0;
  // unsigned int nLinksPerSector60 = theConfiguration_.nInputLinks() / 3;
  // for (auto& sector60 : clusters_SA_permodule) {
  //   unsigned iCluster = 0;
  //   for (auto& cluster : sector60) {
  //     const unsigned empty_frames = 2;
  //     unsigned frame = empty_frames + iCluster / nLinksPerSector60;
  //     unsigned link = iCluster % nLinksPerSector60 + iSector60 * nLinksPerSector60;
  //     if (frame >= theConfiguration_.maxClustersPerLink() + empty_frames)
  //       break;
  //     clusters_SA[frame][link] = std::move(cluster);
  //     ++iCluster;
  //   }
  //   ++iSector60;
  // }
}

void HGCalHistoClusteringWrapper::convertAlgorithmOutputs(
    l1thgcfirmware::HGCalClusterSAPtrCollection& clusterSums,
    l1t::HGCalMulticlusterBxCollection& multiClusters_out,
    const std::map<unsigned int, std::vector<edm::Ptr<l1t::HGCalCluster>>>& inputClustersPtrs) const {
  for (const auto& cluster : clusterSums) {

    if (cluster->w() == 0 || cluster->e() == 0)
      continue;

    // Convert to format expected by L1T
    const auto hwData = cluster->hwCluster();

    // Remove if below pt threshold
    double pt = l1thgcfirmware::Scales::floatEt(hwData.e);
    if ( pt <= theConfiguration_.minClusterPtOut())
      continue;

    double eta = l1thgcfirmware::Scales::floatEta(hwData.w_eta);
    eta *= theConfiguration_.zSide();

    auto sector = theConfiguration_.sector();

    double phi = l1thgcfirmware::Scales::floatPhi(hwData.w_phi);
    phi = rotatePhiFromSectorZero(phi, sector);
    if (theConfiguration_.zSide() == 1) {
      phi = M_PI - phi;
    }
    phi -= (phi > M_PI) ? 2 * M_PI : 0;
    math::PtEtaPhiMLorentzVector clusterP4(pt, eta, phi, 0.);
    l1t::HGCalMulticluster multicluster;
    multicluster.setP4(clusterP4);

    // for (const auto& tc : cluster->constituents()) {
    //   const auto& tc_cmssw = inputClustersPtrs.at(tc->cmsswIndex().first).at(tc->cmsswIndex().second);
    //   // Add tc as constituent, but don't update any other properties of the multicluster i.e. leave them unchanged from those calculated by the emulator
    //   multicluster.addConstituent(tc_cmssw, false, 0.);
    // }

    double emIntfraction = float(hwData.e_em) / hwData.e * multicluster.energy(); // Hack to get multicluster.iPt to return hwData.e_em
    multicluster.saveEnergyInterpretation(l1t::HGCalMulticluster::EnergyInterpretation::EM,
                                          emIntfraction);

    // double emCoreIntfraction = float(cluster->e_em_core()) / cluster->e();
    // multicluster.saveEnergyInterpretation(l1t::HGCalMulticluster::EnergyInterpretation::EM_CORE,
    //                                       emCoreIntfraction * multicluster.energy());

    // double emHEarlyIntfraction = float(cluster->e_h_early()) / cluster->e();
    // multicluster.saveEnergyInterpretation(l1t::HGCalMulticluster::EnergyInterpretation::H_EARLY,
    //                                       emHEarlyIntfraction * multicluster.energy());

    multicluster.zBarycenter(l1thgcfirmware::Scales::floatZ(hwData.w_z));

    multicluster.sigmaRRTot( l1thgcfirmware::Scales::floatSigmaRozRoz(hwData.sigma_roz) );

    // Set cluster shower shape properties
    multicluster.showerLength(cluster->showerLen());
    multicluster.coreShowerLength(cluster->coreShowerLen());
    multicluster.firstLayer(cluster->firstLayer());

    multicluster.setHwData( hwData.pack() );
    multicluster.setHwSector( sector );
    multicluster.setHwZSide( theConfiguration_.zSide() );

    const auto unpacked = l1thgcfirmware::HGCalCluster_HW::unpack(multicluster.getHwData());
    if ( hwData != unpacked ) {
      std::cout << "Problem in unpacking words!" << std::endl;
    }

    const auto hwClusterSumData = cluster->formatClusterSumWords( theConfiguration_ );
    multicluster.setHwClusterSumData( hwClusterSumData );

    multiClusters_out.push_back(0, multicluster);
  }
}

void HGCalHistoClusteringWrapper::process(const std::map<unsigned int, std::vector<edm::Ptr<l1t::HGCalCluster>>>& inputClusters,
                                          std::pair<l1t::HGCalMulticlusterBxCollection&, l1t::HGCalClusterBxCollection&>&
                                              outputMulticlustersAndRejectedClusters) {
  l1thgcfirmware::HGCalLinkTriggerCellSAPtrCollection linkData_in_SA;
  convertCMSSWInputs(inputClusters, linkData_in_SA);

  l1thgcfirmware::HGCalClusterSAPtrCollection clusterSums_out_SA;
  clusterizeHisto(linkData_in_SA, clusterSums_out_SA);

  convertAlgorithmOutputs(clusterSums_out_SA, outputMulticlustersAndRejectedClusters.first, inputClusters);
}

void HGCalHistoClusteringWrapper::clusterizeHisto(
    const l1thgcfirmware::HGCalLinkTriggerCellSAPtrCollection& linkData_in_SA,
    l1thgcfirmware::HGCalClusterSAPtrCollection& clusterSums) {
  // theAlgo_.runAlgorithm(triggerCells_in_SA, clusteredTCs, clusterSums);
  theAlgo_.runAlgorithm( linkData_in_SA, clusterSums );
}

void HGCalHistoClusteringWrapper::configure(
    const std::tuple<const HGCalTriggerGeometryBase* const, const edm::ParameterSet&, const unsigned int, const int>&
        configuration) {
  setGeometry(std::get<0>(configuration));

  theConfiguration_.setNTriggerLayers(std::get<0>(configuration)->lastTriggerLayer());
  theConfiguration_.setTriggerLayers(std::get<0>(configuration)->triggerLayers());

  theConfiguration_.setSector(std::get<2>(configuration));
  theConfiguration_.setZSide(std::get<3>(configuration));

  const edm::ParameterSet pset = std::get<1>(configuration)
                                     .getParameterSet("C3d_parameters")
                                     .getParameterSet("histoMax_C3d_clustering_parameters")
                                     .getParameterSet("layer2FwClusteringParameters");

  // theConfiguration_.setClusterizerOffset(pset.getParameter<unsigned int>("clusterizerOffset"));
  theConfiguration_.setStepLatencies(pset.getParameter<std::vector<unsigned int>>("stepLatencies"));
  theConfiguration_.setCClocks(pset.getParameter<unsigned int>("cClocks"));
  theConfiguration_.setCInputs(pset.getParameter<unsigned int>("cInputs"));
  theConfiguration_.setCInputs2(pset.getParameter<unsigned int>("cInputs2"));
  theConfiguration_.setCInt(pset.getParameter<unsigned int>("cInt"));
  theConfiguration_.setCColumns(pset.getParameter<unsigned int>("cColumns"));
  theConfiguration_.setCRows(pset.getParameter<unsigned int>("cRows"));
  theConfiguration_.setROverZHistOffset(pset.getParameter<unsigned int>("rOverZHistOffset"));
  theConfiguration_.setROverZBinSize(pset.getParameter<unsigned int>("rOverZBinSize"));
  theConfiguration_.setDepths(pset.getParameter<std::vector<unsigned int>>("depths"));
  theConfiguration_.setLayerWeights_E(pset.getParameter<std::vector<unsigned int>>("layerWeights_E"));
  theConfiguration_.setLayerWeights_E_EM(pset.getParameter<std::vector<unsigned int>>("layerWeights_E_EM"));
  theConfiguration_.setLayerWeights_E_EM_core(pset.getParameter<std::vector<unsigned int>>("layerWeights_E_EM_core"));
  theConfiguration_.setLayerWeights_E_H_early(pset.getParameter<std::vector<unsigned int>>("layerWeights_E_H_early"));
  theConfiguration_.setCorrection(pset.getParameter<unsigned int>("correction"));
  theConfiguration_.setSaturation(pset.getParameter<unsigned int>("saturation"));
  const edm::ParameterSet& thresholdParams = pset.getParameterSet("thresholdMaximaParams");
  theConfiguration_.setThresholdParams(thresholdParams.getParameter<unsigned int>("a"),
                                       thresholdParams.getParameter<unsigned int>("b"),
                                       thresholdParams.getParameter<int>("c"));

  // Digitization parameters
  const edm::ParameterSet& digitizationPset = pset.getParameterSet("digiParams");
  theConfiguration_.setROverZRange(digitizationPset.getParameter<double>("rOverZRange"));
  theConfiguration_.setROverZNValues(digitizationPset.getParameter<double>("rOverZNValues"));
  theConfiguration_.setPhiRange(digitizationPset.getParameter<double>("phiRange"));
  theConfiguration_.setPhiNValues(digitizationPset.getParameter<double>("phiNValues"));
  theConfiguration_.setPtDigiFactor(digitizationPset.getParameter<double>("ptDigiFactor"));

  // Parameters for selecting output clusters
  theConfiguration_.setMinClusterPtOut(pset.getParameter<double>("minClusterPtOut"));

  // Input links parameters
  const edm::ParameterSet& inputLinksPset = pset.getParameterSet("inputLinkParams");
  theConfiguration_.setMaxClustersPerLink(inputLinksPset.getParameter<unsigned int>("maxClustersPerLink"));
  theConfiguration_.setNInputLinks(inputLinksPset.getParameter<unsigned int>("nInputLinks"));

  // // TC distribution parameters
  // const edm::ParameterSet& tcDistPset = pset.getParameterSet("tcDistParams");
  // theConfiguration_.setN60Sectors(tcDistPset.getParameter<unsigned int>("n60Sectors"));
  // theConfiguration_.setNCoarsePhiDist1(tcDistPset.getParameter<unsigned int>("nCoarsePhiRegionsDist1"));
  // theConfiguration_.setNDistServers1(tcDistPset.getParameter<unsigned int>("nDistServers1"));
  // theConfiguration_.setDistServer1_nIn(tcDistPset.getParameter<unsigned int>("distServer1_nIn"));
  // theConfiguration_.setDistServer1_nOut(tcDistPset.getParameter<unsigned int>("distServer1_nOut"));
  // theConfiguration_.setDistServer1_nInterleave(tcDistPset.getParameter<unsigned int>("distServer1_nInterleave"));
  // theConfiguration_.setNCoarsePhiDist2(tcDistPset.getParameter<unsigned int>("nCoarsePhiRegionsDist2"));
  // theConfiguration_.setNDistServers2(tcDistPset.getParameter<unsigned int>("nDistServers2"));
  // theConfiguration_.setDistServer2_nIn(tcDistPset.getParameter<unsigned int>("distServer2_nIn"));
  // theConfiguration_.setDistServer2_nOut(tcDistPset.getParameter<unsigned int>("distServer2_nOut"));
  // theConfiguration_.setDistServer2_nInterleave(tcDistPset.getParameter<unsigned int>("distServer2_nInterleave"));

  // // Smearing parameters
  // const edm::ParameterSet& smearingPset = pset.getParameterSet("smearingParams");
  // theConfiguration_.setMaxBinsSmearing1D(smearingPset.getParameter<unsigned int>("maxBinsSmearing1D"));
  // theConfiguration_.setNBitsAreaNormLUT(smearingPset.getParameter<unsigned int>("nBitsAreaNormLUT"));

  // // Clusterizer parameters
  // const edm::ParameterSet& clusterizerPset = pset.getParameterSet("clusterizerParams");
  // theConfiguration_.setNBinsCosLUT(clusterizerPset.getParameter<unsigned int>("nBinsCosLUT"));
  // theConfiguration_.setNBitsCosLUT(clusterizerPset.getParameter<unsigned int>("nBitsCosLUT"));
  // theConfiguration_.setNFifos(clusterizerPset.getParameter<unsigned int>("nFifos"));
  // theConfiguration_.setNColumnsPerFifo(clusterizerPset.getParameter<unsigned int>("nColumnsPerFifo"));
  // theConfiguration_.setClusterizerMagicTime(clusterizerPset.getParameter<unsigned int>("clusterizerMagicTime"));
  // theConfiguration_.setFirstSeedBin(clusterizerPset.getParameter<unsigned int>("firstSeedBin"));
  // theConfiguration_.setNColumnFifoVeto(clusterizerPset.getParameter<unsigned int>("nColumnsFifoVeto"));
  // theConfiguration_.setDeltaR2Cut(clusterizerPset.getParameter<unsigned int>("deltaR2Cut"));
  // theConfiguration_.setNColumnsForClustering(clusterizerPset.getParameter<unsigned int>("nColumnsForClustering"));
  // theConfiguration_.setNRowsForClustering(clusterizerPset.getParameter<unsigned int>("nRowsForClustering"));

  theConfiguration_.initializeLUTs();
};

double HGCalHistoClusteringWrapper::rotatePhiToSectorZero(const double phi, const unsigned sector) const {
  double rotatedPhi = phi;
  if (sector == 1) {
    if (rotatedPhi < M_PI and rotatedPhi > 0)
      rotatedPhi = rotatedPhi - (2. * M_PI / 3.);
    else
      rotatedPhi = rotatedPhi + (4. * M_PI / 3.);
  } else if (sector == 2) {
    rotatedPhi = rotatedPhi + (2. * M_PI / 3.);
  }
  return rotatedPhi;
}

double HGCalHistoClusteringWrapper::rotatePhiFromSectorZero(const double phi, const unsigned sector) const {
  double rotatedPhi = phi;
  if (sector == 1) {
    rotatedPhi += (2. * M_PI / 3.);
  } else if (sector == 2) {
    rotatedPhi += (4. * M_PI / 3.);
  }
  return rotatedPhi;
}

DEFINE_EDM_PLUGIN(HGCalHistoClusteringWrapperBaseFactory, HGCalHistoClusteringWrapper, "HGCalHistoClusteringWrapper");
