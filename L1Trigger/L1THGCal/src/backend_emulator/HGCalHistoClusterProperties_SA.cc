#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusterProperties_SA.h"

#include <cmath>
#include <algorithm>
#include <bitset>

using namespace std;
using namespace l1thgcfirmware;

HGCalHistoClusterProperties::HGCalHistoClusterProperties(const ClusterAlgoConfig& config) : config_(config) {}

void HGCalHistoClusterProperties::runClusterProperties( HGCalClusterSAPtrCollection& clusters) const {

  // Cluster properties
  clusterProperties(clusters);
}

// Calculates properties of clusters from accumulated quantities
void HGCalHistoClusterProperties::clusterProperties(HGCalClusterSAPtrCollection& clusters) const {

  for (auto& c : clusters) {

    HGCalCluster_HW& hwCluster = c->hwCluster();
    hwCluster.e = Scales::HGCaltoL1_et(c->e().value_);
    hwCluster.e_em = Scales::HGCaltoL1_et(c->e_em().value_);
    hwCluster.fractionInCE_E = Scales::makeL1EFraction(c->e_em().value_, c->e().value_);
    hwCluster.fractionInCoreCE_E = Scales::makeL1EFraction(c->e_em_core().value_, c->e_em().value_);
    hwCluster.fractionInEarlyCE_E = Scales::makeL1EFraction(c->e_h_early().value_, c->e().value_);
    hwCluster.setGCTBits();
    std::vector<int> layeroutput = showerLengthProperties(c->layerbits().value_);
    c->set_firstLayer(layeroutput[0]);
    c->set_lastLayer(layeroutput[1]);
    c->set_showerLen(layeroutput[2]);
    c->set_coreShowerLen(layeroutput[3]);
    hwCluster.firstLayer = c->firstLayer();
    hwCluster.lastLayer = c->lastLayer();
    hwCluster.showerLength = c->showerLen();
    hwCluster.coreShowerLength = c->coreShowerLen();
    hwCluster.nTC = c->n_tc().value_;

    if (c->n_tc_w() == 0)
      continue;

    hwCluster.w_eta = convertRozToEta( c );
    bool saturatedPhi = false;
    bool nominalPhi = false;
    hwCluster.w_phi = Scales::HGCaltoL1_phi(float(c->wphi().value_)/c->w().value_, saturatedPhi, nominalPhi);
    hwCluster.w_z = Scales::HGCaltoL1_z( float(c->wz().value_) / c->w().value_ );

    // Quality flags are placeholders at the moment
    hwCluster.setQualityFlags(Scales::HGCaltoL1_et(c->e_em_core().value_), Scales::HGCaltoL1_et(c->e_h_early().value_), c->sat_tc().value_, c->shapeq().value_, saturatedPhi, nominalPhi);

    const double sigma_E_scale = 0.008982944302260876;
    hwCluster.sigma_E = sigma_coordinate(c->n_tc_w().value_, c->w2().value_, c->w().value_, sigma_E_scale);

    const double sigma_z_scale = 0.08225179463624954;
    hwCluster.sigma_z = sigma_coordinate(c->w().value_, c->wz2().value_, c->wz().value_, sigma_z_scale);

    const double sigma_phi_scale = 0.907465934753418;
    hwCluster.sigma_phi = sigma_coordinate(c->w().value_, c->wphi2().value_, c->wphi().value_, sigma_phi_scale);

    hwCluster.sigma_eta = convertSigmaRozRozToSigmaEtaEta(c);

    const double sigma_roz_scale = 0.5073223114013672;
    unsigned int sigma_roz = sigma_coordinate(c->w().value_, c->wroz2().value_, c->wroz().value_, sigma_roz_scale);
    // Emulation of a bug in firmware
    // if ( sigma_roz >=256 ) sigma_roz -= 256;
    while (sigma_roz >= 256) sigma_roz -= 256;
    if ( sigma_roz > 127 ) sigma_roz = 127;
    hwCluster.sigma_roz = sigma_roz;
  }
}

unsigned int HGCalHistoClusterProperties::sigma_coordinate(unsigned int w,
                                                            unsigned long int wc2,
                                                            unsigned int wc,
                                                            double scale ) const {
  if ( w == 0 ) return 0;
  unsigned int sigma = round(sqrt( (float(w)*float(wc2) - float(wc) * float(wc))  / ( float(w) * float(w) ) ) * scale);
  return sigma;
}

std::vector<int> HGCalHistoClusterProperties::showerLengthProperties(unsigned long int layerBits) const {
  int counter = 0;
  int firstLayer = 0;
  bool firstLayerFound = false;
  int lastLayer = 0;
  std::vector<int> layerBits_array;

  bitset<34> layerBitsBitset(layerBits);
  for (size_t i = 0; i < layerBitsBitset.size(); ++i) {
      bool bit = layerBitsBitset[34-1-i];
      if ( bit ) {
        if ( !firstLayerFound ) {
          firstLayer = i + 1;
          firstLayerFound = true;
        }
        lastLayer = i+1;
        counter += 1;
      } else {
        layerBits_array.push_back(counter);
        counter = 0;
      }
  }

  int showerLen = lastLayer - firstLayer + 1;
  int coreShowerLen = config_.nTriggerLayers();
  if (!layerBits_array.empty()) {
    coreShowerLen = *std::max_element(layerBits_array.begin(), layerBits_array.end());
  }
  return {firstLayer, lastLayer, showerLen, coreShowerLen};
}

double HGCalHistoClusterProperties::convertRozToEta( HGCalClusterSAPtr& cluster ) const {
  // TODO : named constants for magic numbers
  double roz = double(cluster->wroz().value_)/cluster->w().value_;
  if ( roz < 1026.9376220703125 ) roz = 1026.9376220703125;
  else if ( roz > 5412.17138671875 ) roz = 5412.17138671875;
  roz -= 1026.9376220703125;
  roz *= 0.233510936;
  roz = int(round(roz));
  if ( roz > 1023 ) roz = 1023;
  return config_.rozToEtaLUT(roz);
}

double HGCalHistoClusterProperties::convertSigmaRozRozToSigmaEtaEta( HGCalClusterSAPtr& cluster ) const {
  // TODO : named constants for magic numbers
  // Sigma eta eta calculation
  double roz = cluster->wroz().value_/cluster->w().value_;
  const double min_roz = 809.9324340820312;
  const double max_roz = 4996.79833984375;
  if ( roz < min_roz ) roz = min_roz;
  else if ( roz > max_roz ) roz = max_roz;
  roz -= min_roz;
  const double scale = 0.015286154113709927;
  roz *= scale;
  roz = int(round(roz));
  if ( roz > 63 ) roz = 63;

  const double sigma_roz_scale = 0.220451220870018;
  double sigmaRoz = sigma_coordinate(cluster->w().value_, cluster->wroz2().value_, cluster->wroz().value_, sigma_roz_scale);

  sigmaRoz = int(round(sigmaRoz));
  if ( sigmaRoz > 63 ) sigmaRoz = 63;
  unsigned int lutAddress = roz * 64 + sigmaRoz;
  if ( lutAddress >= 4096 ) lutAddress = 4095;
  return config_.sigmaRozToSigmaEtaLUT(lutAddress);
}
