#ifndef L1Trigger_L1THGCal_HGCalCluster_SA_h
#define L1Trigger_L1THGCal_HGCalCluster_SA_h

#include <vector>
#include <memory>
#include <iostream>
namespace l1thgcfirmware {

  class HGCalCluster {
  public:
    HGCalCluster( unsigned int clock,
                  unsigned int index,
                  bool frameValid,
                  bool dataValid)
                : clock_( clock ),
                  index_(index),
                  layerBits_(0),
                  frameValid_(frameValid),
                  dataValid_(dataValid),
                  n_tc_(0),
                  e_(0),
                  e_em_(0),
                  e_em_core_(0),
                  e_h_early_(0),
                  w_(0),
                  n_tc_w_(0),
                  w2_(0),
                  wz_(0),
                  weta_(0),
                  wphi_(0),
                  wroz_(0),
                  wz2_(0),
                  weta2_(0),
                  wphi2_(0),
                  wroz2_(0),
                  layerbits_(0),
                  sat_tc_(false),
                  shapeq_(1)
                  {}

    ~HGCalCluster(){}

    // Setters
    void setClock( const unsigned int clock ) { clock_ = clock; }
    void setIndex( const unsigned int index ) { index_ = index; }
    void setDataValid( const bool dataValid ) { dataValid_ = dataValid; }

    void set_n_tc( unsigned int n_tc ) { n_tc_ = n_tc; }
    void set_e( unsigned int e ) { e_ = e; }
    void set_e_em( unsigned int e_em ) { e_em_ = e_em; }
    void set_e_em_core( unsigned int e_em_core ) { e_em_core_ = e_em_core; }
    void set_e_h_early( unsigned int e_h_early ) { e_h_early_ = e_h_early; }
    void set_w( unsigned int w ) { w_ = w; }
    void set_n_tc_w( unsigned int n_tc_w ) { n_tc_w_ = n_tc_w; }
    void set_w2( unsigned int w2 ) { w2_ = w2; }
    void set_wz( unsigned int wz ) { wz_ = wz; }
    void set_weta( unsigned int weta ) { weta_ = weta; }
    void set_wphi( unsigned int wphi ) { wphi_ = wphi; }
    void set_wroz( unsigned int wroz ) { wroz_ = wroz; }
    void set_wz2( unsigned int wz2 ) { wz2_ = wz2; }
    void set_weta2( unsigned int weta2 ) { weta2_ = weta2; }
    void set_wphi2( unsigned int wphi2 ) { wphi2_ = wphi2; }
    void set_wroz2( unsigned int wroz2 ) { wroz2_ = wroz2; }
    void set_layerbits( unsigned int layerbits ) { layerbits_ = layerbits; }
    void set_sat_tc( bool sat_tc ) { sat_tc_ = sat_tc; }
    void set_shapeq( unsigned int shapeq ) { shapeq_ = shapeq; }

    // Getters
    unsigned int clock() const { return clock_; }
    unsigned int index() const { return index_; }
    bool frameValid() const { return frameValid_; }
    bool dataValid() const { return dataValid_; }

    unsigned int n_tc() const { return n_tc_; }
    unsigned int e() const { return e_; }
    unsigned int e_em() const { return e_em_; }
    unsigned int e_em_core() const { return e_em_core_; }
    unsigned int e_h_early() const { return e_h_early_; }
    unsigned int w() const { return w_; }
    unsigned int n_tc_w() const { return n_tc_w_; }
    unsigned int w2() const { return w2_; }
    unsigned int wz() const { return wz_; }
    unsigned int weta() const { return weta_; }
    unsigned int wphi() const { return wphi_; }
    unsigned int wroz() const { return wroz_; }
    unsigned int wz2() const { return wz2_; }
    unsigned int weta2() const { return weta2_; }
    unsigned int wphi2() const { return wphi2_; }
    unsigned int wroz2() const { return wroz2_; }
    unsigned int layerbits() const { return layerbits_; }
    bool sat_tc() const { return sat_tc_; }
    unsigned int shapeq() const { return shapeq_; }

    // Operators
    const HGCalCluster& operator+=(const HGCalCluster& hc);

  private:
    unsigned int clock_;
    unsigned int index_;
    unsigned int layerBits_;
    bool frameValid_;
    bool dataValid_;

    unsigned int n_tc_;
    unsigned int e_;
    unsigned int e_em_;
    unsigned int e_em_core_;
    unsigned int e_h_early_;
    unsigned int w_;
    unsigned int n_tc_w_;
    unsigned int w2_;
    unsigned int wz_;
    unsigned int weta_;
    unsigned int wphi_;
    unsigned int wroz_;
    unsigned int wz2_;
    unsigned int weta2_;
    unsigned int wphi2_;
    unsigned int wroz2_;
    unsigned int layerbits_;
    bool sat_tc_;
    unsigned int shapeq_;

  };

  typedef std::vector<HGCalCluster> HGCalClusterSACollection;
  typedef std::shared_ptr<HGCalCluster> HGCalClusterSAPtr;
  typedef std::vector<HGCalClusterSAPtr > HGCalClusterSAPtrCollection;
  typedef std::vector< HGCalClusterSAPtrCollection > HGCalClusterSAPtrCollections;

}  // namespace l1thgcfirmware

#endif
