#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalCluster_SA.h"
#include <algorithm>
#include <cmath>

// Many useful functions for handling and manipulation of bitsets
// From L1 Track Finder DTC
// But could perhaps still just use ap_uint?
#include "DataFormats/L1TrackTrigger/interface/TTBV.h"

using namespace l1thgcfirmware;

const HGCalCluster& HGCalCluster::operator+=(HGCalCluster& c) {
  // Not handling field widths
  HGCalCluster original(*this);
  this->set_n_tc(this->n_tc() + c.n_tc());
  this->set_e(this->e() + c.e());
  this->set_e_em(this->e_em() + c.e_em());
  this->set_e_em_core(this->e_em_core() + c.e_em_core());
  this->set_e_h_early(this->e_h_early() + c.e_h_early());
  this->set_w(this->w() + c.w());
  this->set_n_tc_w(this->n_tc_w() + c.n_tc_w());
  this->set_w2(this->w2() + c.w2());
  this->set_wz(this->wz() + c.wz());
  this->set_weta(this->weta() + c.weta());
  this->set_wphi(this->wphi() + c.wphi());
  this->set_wroz(this->wroz() + c.wroz());
  this->set_wz2(this->wz2() + c.wz2());
  this->set_weta2(this->weta2() + c.weta2());
  this->set_wphi2(this->wphi2() + c.wphi2());
  this->set_wroz2(this->wroz2() + c.wroz2());

  this->set_layerbits(this->layerbits() | c.layerbits());
  this->set_sat_tc(this->sat_tc() | c.sat_tc());
  this->set_shapeq(this->shapeq() | c.shapeq());

  const unsigned clusterWeightSat80 = ((1 << 16) - 1) * 0.8;  // 52428
  if (w_ <= clusterWeightSat80 && original.shapeq() == 1 && c.shapeq() == 1) {
    this->set_shapeq(1);
  } else {
    this->set_shapeq(0);

    if (this->w() > c.w()) {
      this->set_w(original.w());
      this->set_w2(original.w2());
      this->set_wz(original.wz());
      this->set_weta(original.weta());
      this->set_wphi(original.wphi());
      this->set_wroz(original.wroz());
      this->set_wz2(original.wz2());
      this->set_weta2(original.weta2());
      this->set_wphi2(original.wphi2());
      this->set_wroz2(original.wroz2());
      this->set_n_tc_w(original.n_tc_w());
    } else {
      this->set_w(c.w());
      this->set_w2(c.w2());
      this->set_wz(c.wz());
      this->set_weta(c.weta());
      this->set_wphi(c.wphi());
      this->set_wroz(c.wroz());
      this->set_wz2(c.wz2());
      this->set_weta2(c.weta2());
      this->set_wphi2(c.wphi2());
      this->set_wroz2(c.wroz2());
      this->set_n_tc_w(c.n_tc_w());
    }
  }

  for (const auto& constituent : c.constituents()) {
    this->add_constituent(constituent);
  }

  return *this;
}

void HGCalCluster::clearClusterWords() {
  for ( auto& clusterWord : packedData_ ) {
    clusterWord = 0;
  }
}

HGCalCluster::ClusterWords HGCalCluster::formatClusterWords( const ClusterAlgoConfig& config ) {
  clearClusterWords();
  packedData_[0] = formatFirstWord( config );
  packedData_[1] = formatSecondWord( config );
  packedData_[2] = formatThirdWord( config );
  packedData_[3] = formatFourthWord( config );

  return packedData_;
}

HGCalCluster::ClusterWord HGCalCluster::formatFirstWord( const ClusterAlgoConfig& config ) {

  // const TTBV hw_Et(round(e()/config.ptDigiFactor()), 0.25, 14);
  int pt = round(e()/config.ptDigiFactor()/0.25);
  const TTBV hw_Et(pt, 14, false);
  std::cout << "Converting pt : " << e() << " " << e()/config.ptDigiFactor() << " " << round(e()/config.ptDigiFactor())/0.25 << " " << hw_Et.str() << std::endl;
  // const TTBV hw_EtEgamma(e_em()/config.ptDigiFactor(), 0.25, 14);
  const TTBV hw_EtEgamma(0, 14);
  const TTBV hw_GCTSelectBits(0, 4);

  // const TTBV hw_fracCEE(int(E_EM_over_E_Fraction()), 8);
  // const TTBV hw_fracCEECore(int(E_EM_core_over_E_EM_Fraction()), 8);
  // const TTBV hw_fracCEHFront(int(E_H_early_over_E_Fraction()), 8);
  const TTBV hw_fracCEE(0, 8);
  const TTBV hw_fracCEECore(0, 8);
  const TTBV hw_fracCEHFront(0, 8);
  const TTBV hw_firstLayerZ(0, 6);
  const TTBV hw_spare(0, 2);

  // ClusterWord clusterWord = TTBV(hw_Et + hw_EtEgamma + hw_GCTSelectBits + hw_fracCEE + hw_fracCEECore + w_fracCEHFront + hw_firstLayerZ + hw_spare).bs();
  ClusterWord clusterWord = TTBV(
    hw_spare +
    hw_firstLayerZ + 
    hw_fracCEHFront + 
    hw_fracCEECore + 
    hw_fracCEECore + 
    hw_GCTSelectBits + 
    hw_EtEgamma + 
    hw_Et
    ).bs();

  return clusterWord;
}
HGCalCluster::ClusterWord HGCalCluster::formatSecondWord( const ClusterAlgoConfig& config ) {

  // Not being careful with round/floor here
  constexpr float ETAPHI_LSB = M_PI / 720;
  double rOverZ = ( 1.0 * wroz() / w() ) * config.rOverZRange() / config.rOverZNValues();
  double eta = -1.0 * std::log( tan( atan( rOverZ ) / 2 ) );
  // eta *= theConfiguration_.zSide();

  // Interface document definition
  // eta -= 1.5; 
  // int l1Eta = eta / ETAPHI_LSB;
  // const TTBV hw_Eta(l1Eta, 9, false);
  // double etaTemp = eta - 1.5;
  // const TTBV hw_Eta_interface(int(etaTemp / ETAPHI_LSB), 9, false);

  // Current correlator definition
  eta -= 2.25;
  int l1Eta = round( eta / ETAPHI_LSB );
  const TTBV hw_Eta(l1Eta, 9, true);
  std::cout << "Global eta, l1 eta : " << -1.0 * std::log( tan( atan( rOverZ ) / 2 ) ) << " " << l1Eta << " " << hw_Eta.str() << std::endl;
  // int l1Phi = ( ( 1.0 * wphi() / w() ) * config.phiRange() / config.phiNValues() - M_PI/3 )  / ETAPHI_LSB;

  double phi = ( 1.0 * wphi() / w() ) * config.phiRange() / config.phiNValues();
  if ( config.zSide() == 1 ) {
    phi = M_PI - phi;
  }
  // phi -= ( phi > M_PI ) ? 2 * M_PI : 0;

  int l1Phi = round( ( phi - M_PI/2  )  / ETAPHI_LSB );


  std::cout << "Global phi, l1 phi : " << ( 1.0 * wphi() / w() ) * config.phiRange() / config.phiNValues() << " " << l1Phi << std::endl;

  // std::cout << "Phi : " << 1.0 * wphi() / w() << " " << Mean_phi_Quotient() << " " << Mean_phi_Fraction() << td::endl;
  // std::cout << l1Phi << std::endl;
  const TTBV hw_Phi(l1Phi, 9, true);
  const TTBV hw_z(0, 12);
  const TTBV hw_spare1(0, 2);

  // const TTBV hw_nTCs(int(n_tc()), 10);
  const TTBV hw_nTCs(0, 10);
  const TTBV hw_qualityFlags(0, 10);
  const TTBV hw_spare2(0, 12);

  // ClusterWord clusterWord = TTBV(hw_Eta + hw_Phi + hw_z + hw_spare1 + hw_nTCs + hw_qualityFlags + hw_spare2).bs();
  ClusterWord clusterWord = TTBV(    
    hw_spare2 +
    hw_qualityFlags +
    hw_nTCs + 
    hw_spare1 + 
    hw_z + 
    hw_Phi + 
    hw_Eta
    ).bs();

  return clusterWord;
}
HGCalCluster::ClusterWord HGCalCluster::formatThirdWord( const ClusterAlgoConfig& config ) {
  ClusterWord clusterWord = 0;
  return clusterWord;
}

HGCalCluster::ClusterWord HGCalCluster::formatFourthWord( const ClusterAlgoConfig& config ) {
  ClusterWord clusterWord = 0;
  return clusterWord;
}