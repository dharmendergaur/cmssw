#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalCluster_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/Toolkit.h"
// #include <math.h>
// #include <algorithm>

  
using namespace l1thgcfirmware;

#define MEMBERS &HGCalCluster::clock_ , \
                &HGCalCluster::index_ , \
                &HGCalCluster::n_tc_, \
                &HGCalCluster::e_, \
                &HGCalCluster::e_em_, \
                &HGCalCluster::e_em_core_, \
                &HGCalCluster::e_h_early_, \
                &HGCalCluster::w_, \
                &HGCalCluster::n_tc_w_, \
                &HGCalCluster::w2_, \
                &HGCalCluster::wz_, \
                &HGCalCluster::wphi_, \
                &HGCalCluster::wroz_, \
                &HGCalCluster::wz2_, \
                &HGCalCluster::wphi2_, \
                &HGCalCluster::wroz2_, \
                &HGCalCluster::layerbits_, \
                &HGCalCluster::sat_tc_, \
                &HGCalCluster::shapeq_, \
                &HGCalCluster::sortKey_, \
                &HGCalCluster::sortKey2_, \
                &HGCalCluster::L_, \
                &HGCalCluster::R_, \
                &HGCalCluster::X_, \
                &HGCalCluster::lastFrame_ , \
                &HGCalCluster::dataValid_
    // &HGCalCluster::weta_, 
    // &HGCalCluster::weta2_, 

bool HGCalCluster::operator==(const HGCalCluster& other) const { 
  return ::check( *this , other , MEMBERS );
}

std::ostream& operator<< ( std::ostream& aStr , const HGCalCluster& aCell )
{
  return ::print( aStr , aCell ,  MEMBERS );
}

HGCalClusterSAPtrCollection HGCalCluster::ReadDebugFile( const std::string& aFilename )
{
  HGCalClusterSAPtrCollection lRet;
  ::ReadDebugFile( aFilename , lRet , MEMBERS );
  return lRet;
}

#undef MEMBERS







const HGCalCluster& HGCalCluster::operator+=(const HGCalCluster& c) {

  n_tc_ += c.n_tc_;
  e_ += c.e_;
  e_em_ += c.e_em_;
  e_em_core_ += c.e_em_core_;
  e_h_early_ += c.e_h_early_;
  w_ += c.w_;
  // w2_ += c.w2_;
  // wz_ += c.wz_;
  // weta_ += c.weta_;
  // wphi_ += c.wphi_;
  // wroz_ += c.wroz_;
  // wz2_ += c.wz2_;
  // weta2_ += c.weta2_;
  // wphi2_ += c.wphi2_;
  // wroz2_ += c.wroz2_;
  // n_tc_w_ += c.n_tc_w_;
  layerbits_ |= c.layerbits_;
  sat_tc_ |= c.sat_tc_;
  shapeq_ |= c.shapeq_;

  // // if ( w_ <= 52438 and shapeq_ == 1 and c.shapeq_ == 1 ) { // Magic numbers
    // // shapeq_ = 1;
    // w_ = ...
    w2_ += c.w2_;
    wz_ += c.wz_;
    // weta_ += c.weta_;
    wphi_ += c.wphi_;
    wroz_ += c.wroz_;
    wz2_ += c.wz2_;
    // weta2_ += c.weta2_;
    wphi2_ += c.wphi2_;
    wroz2_ += c.wroz2_;
    n_tc_w_ += c.n_tc_w_;      
  // // } else if ( w_ <= c.w_ ) {
    // // shapeq_ = 0;
    // // w_ = c.w_;
    // // w2_ = c.w2_;
    // // wz_ = c.wz_;
    // // // weta_ = c.weta_;
    // // wphi_ = c.wphi_;
    // // wroz_ = c.wroz_;
    // // wz2_ = c.wz2_;
    // // // weta2_ = c.weta2_;
    // // wphi2_ = c.wphi2_;
    // // wroz2_ = c.wroz2_;
    // // n_tc_w_ = c.n_tc_w_;
  // // } else {
    // // shapeq_ = 0;
  // // }

  return *this;
}


void HGCalCluster::saturate()
{
  if( n_tc_.value_ > 0x3FF ) n_tc_.value_ = 0x3FF;
  if( e_.value_ > 0x3FFFFF ) e_.value_ = 0x3FFFFF;
  if( e_em_.value_ > 0x3FFFFF ) e_em_.value_ = 0x3FFFFF;
  if( e_em_core_.value_ > 0x3FFFFF ) e_em_core_.value_ = 0x3FFFFF;
  if( e_h_early_.value_ > 0x3FFFFF ) e_h_early_.value_ = 0x3FFFFF;
  if( w_.value_ > 0xFFFF ) w_.value_ = 0xFFFF;
  if( n_tc_w_.value_ > 0x3FF ) n_tc_w_.value_ = 0x3FF;
  if( w2_.value_ > 0xFFFFFFFF ) w2_.value_ = 0xFFFFFFFF;
  if( wz_.value_ > 0x1FFFFFFF ) wz_.value_ = 0x1FFFFFFF;
  if( wphi_.value_ > 0xFFFFFFF ) wphi_.value_ = 0xFFFFFFF;
  if( wroz_.value_ > 0x1FFFFFFF ) wroz_.value_ = 0x1FFFFFFF;
  if( wz2_.value_ > 0x3FFFFFFFFFF ) wz2_.value_ = 0x3FFFFFFFFFF;
  if( wphi2_.value_ > 0xFFFFFFFFFF ) wphi2_.value_ = 0xFFFFFFFFFF;
  if( wroz2_.value_ > 0x3FFFFFFFFFF ) wroz2_.value_ = 0x3FFFFFFFFFF;      
}

void HGCalCluster::clearClusterSumWords() {
  for ( auto& clusterSumWord : packedData_clustersSums_ ) {
    clusterSumWord = 0;
  }
}

HGCalCluster::ClusterSumWords HGCalCluster::formatClusterSumWords( const ClusterAlgoConfig& config ) {
  clearClusterSumWords();

  const unsigned nb_nTCs = 10;
  const unsigned nb_e = 22;
  const unsigned nb_e_em = 22;
  const unsigned nb_e_em_core = 22;
  const unsigned nb_e_h_early = 22;

  const unsigned nb_w = 16;
  const unsigned nb_n_tc_w = 10;
  const unsigned nb_w2 = 32;
  const unsigned nb_wz = 29;
  // const unsigned nb_weta = 26;
  const unsigned nb_wphi = 28;
  const unsigned nb_wroz = 29;

  const unsigned nb_wz2 = 42;
  // const unsigned nb_weta2 = 36;
  const unsigned nb_wphi2 = 40;
  const unsigned nb_wroz2 = 42;
  const unsigned nb_layerbits = 34;
  const unsigned nb_sat_tc = 1;
  const unsigned nb_shapeq = 1;

  ap_uint<nb_nTCs> hw_nTCs = n_tc().value_;
  ap_uint<nb_e> hw_e = e().value_;
  ap_uint<nb_e_em> hw_e_em = e_em().value_;
  ap_uint<nb_e_em_core> hw_e_em_core = e_em_core().value_;
  ap_uint<nb_e_h_early> hw_e_h_early = e_h_early().value_;

  ap_uint<nb_w> hw_w = w().value_;
  ap_uint<nb_n_tc_w> hw_n_tc_w = n_tc_w().value_;
  ap_uint<nb_w2> hw_w2 = w2().value_;
  ap_uint<nb_wz> hw_wz = wz().value_;
  // ap_uint<nb_weta> hw_weta = 0;
  ap_uint<nb_wphi> hw_wphi = wphi().value_ * 3456. / 1944; // Temp hack (another one...)
  ap_uint<nb_wroz> hw_wroz = wroz().value_;

  ap_uint<nb_wz2> hw_wz2 = wz2().value_;
  // ap_uint<nb_weta2> hw_weta2 = 0;
  ap_uint<nb_wphi2> hw_wphi2 = wphi2().value_;
  ap_uint<nb_wroz2> hw_wroz2 = wroz2().value_;
  ap_uint<nb_layerbits> hw_layerbits = layerbits().value_;
  ap_uint<nb_sat_tc> hw_sat_tc = sat_tc().value_;
  ap_uint<nb_shapeq> hw_shapeq = shapeq().value_;

  const ap_uint<allClusterSumWordsLength> clusterSumRecord = (
    hw_shapeq,
    hw_sat_tc,
    hw_layerbits,
    hw_wroz2,
    hw_wphi2,
    // hw_weta2,
    hw_wz2,
    hw_wroz,
    hw_wphi,
    // hw_weta,
    hw_wz,
    hw_w2,
    hw_n_tc_w,
    hw_w,
    hw_e_h_early,
    hw_e_em_core,
    hw_e_em,
    hw_e,
    hw_nTCs
  );

  for ( unsigned iWord = 0; iWord < nWordsPerClusterSum; ++iWord ) {
    ClusterSumWord word = clusterSumRecord.range((iWord+1)*(clusterSumWordLength)-1,iWord*clusterSumWordLength).to_ulong();
    packedData_clustersSums_[iWord] = word;
  }
  // std::cout << "=== Cluster info ===" << std::endl;
  // std::cout << "NTCs : " << n_tc() << " " << hw_nTCs << " " << hw_nTCs.to_string() << std::endl;
  // std::cout << "E : " << e() << " " << hw_e << " " << hw_e.to_string() << std::endl;
  // std::cout << "w : " << w() << " " << hw_w << " " << hw_w.to_string() << std::endl;
  // std::cout << "n_tc_w : " << n_tc_w() << " " << hw_n_tc_w << " " << hw_n_tc_w.to_string() << std::endl;
  // std::cout << "w2 : " << w2() << " " << hw_w2 << " " << hw_w2.to_string() << std::endl;
  // std::cout << "Cluster sum record : " << clusterSumRecord << " " << clusterSumRecord.to_string()<< std::endl;

  return packedData_clustersSums_;
}
