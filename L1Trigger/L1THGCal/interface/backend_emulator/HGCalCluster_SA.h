#ifndef L1Trigger_L1THGCal_HGCalCluster_SA_h
#define L1Trigger_L1THGCal_HGCalCluster_SA_h

#include "L1Trigger/L1THGCal/interface/backend_emulator/Binary.h"
#include "DataFormats/L1THGCal/interface/HGCalCluster_HW.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringConfig_SA.h"

#include <vector>
#include <memory>


// Forward declare for the sake of our friends
namespace l1thgcfirmware {
  class HGCalCluster;
  typedef std::unique_ptr<HGCalCluster> HGCalClusterSAPtr;
  typedef std::vector<HGCalClusterSAPtr> HGCalClusterSAPtrCollection;
  typedef std::shared_ptr<HGCalCluster> HGCalClusterSAShrPtr;
  typedef std::vector<HGCalClusterSAShrPtr> HGCalClusterSAShrPtrCollection;
}

std::ostream& operator<< ( std::ostream& aStr , const l1thgcfirmware::HGCalCluster& aCell );



namespace l1thgcfirmware {

  class HGCalCluster {
  public:
    HGCalCluster(unsigned int clock, unsigned int index, bool lastFrame, bool dataValid) :
                  clock_(clock),
                  index_(index),
                  n_tc_(0),
                  e_(0),
                  e_em_(0),
                  e_em_core_(0),
                  e_h_early_(0),
                  w_(0),
                  n_tc_w_(0),
                  w2_(0),
                  wz_(0),
                  wphi_(0),
                  wroz_(0),
                  wz2_(0),
                  wphi2_(0),
                  wroz2_(0),
                  layerbits_(0),
                  sat_tc_(false),
                  shapeq_(1),
                  first_layer_(0),
                  last_layer_(0),
                  shower_len_(0),
                  core_shower_len_(0),
                  fractionInCE_E_(0),
                  fractionInCoreCE_E_(0),
                  fractionInEarlyCE_E_(0),
                  sortKey_(0),
                  sortKey2_(0),
                  L_(false), 
                  R_(false), 
                  X_(false), 
                  lastFrame_(lastFrame),
                  dataValid_(dataValid)
                  {
                    hwCluster_.clear();
                  }
    HGCalCluster() : HGCalCluster( 0, 0, 0, 0) {};

    HGCalCluster( const HGCalCluster& aOther /*!< Anonymous argument */ ) = default;
    HGCalCluster& operator = (const HGCalCluster& aOther /*!< Anonymous argument */ ) = default;
    HGCalCluster( HGCalCluster&& aOther /*!< Anonymous argument */ ) = default;
    HGCalCluster& operator = ( HGCalCluster&& aOther /*!< Anonymous argument */ ) = default;

    ~HGCalCluster(){}

    // Firmware representation of cluster sum (input to cluster properties step)
    static constexpr int clusterSumWordLength = 64;
    static constexpr int nWordsPerClusterSum = 8;
    static constexpr int allClusterSumWordsLength = clusterSumWordLength * nWordsPerClusterSum;
    typedef uint64_t ClusterSumWord;
    typedef std::array<ClusterSumWord, nWordsPerClusterSum> ClusterSumWords;

    // Setters
    void setClock( const unsigned int clock ) { clock_ = clock; }
    void setIndex( const unsigned int index ) { index_ = index; }
    void setDataValid( const bool dataValid ) { dataValid_ = dataValid; }

    void set_n_tc( binary n_tc ) { n_tc_ = n_tc; }
    void set_e( binary e ) { e_ = e; }
    void set_e_em( binary e_em ) { e_em_ = e_em; }
    void set_e_em_core( binary e_em_core ) { e_em_core_ = e_em_core; }
    void set_e_h_early( binary e_h_early ) { e_h_early_ = e_h_early; }
    void set_w( binary w ) { w_ = w; }
    void set_n_tc_w( binary n_tc_w ) { n_tc_w_ = n_tc_w; }
    void set_w2( binary w2 ) { w2_ = w2; }
    void set_wz( binary wz ) { wz_ = wz; }
    void set_wphi( binary wphi ) { wphi_ = wphi; }
    void set_wroz( binary wroz ) { wroz_ = wroz; }
    void set_wz2( binary wz2 ) { wz2_ = wz2; }
    void set_wphi2( binary wphi2 ) { wphi2_ = wphi2; }
    void set_wroz2( binary wroz2 ) { wroz2_ = wroz2; }
    void set_layerbits( binary layerbits ) { layerbits_ = layerbits; }
    void set_sat_tc( binary sat_tc ) { sat_tc_ = sat_tc; }
    void set_shapeq( binary shapeq ) { shapeq_ = shapeq; }
    void set_firstLayer(unsigned long int FirstLayer) { first_layer_ = FirstLayer; }
    void set_lastLayer(unsigned long int LastLayer) { last_layer_ = LastLayer; }
    void set_showerLen(unsigned long int ShowerLen) { shower_len_ = ShowerLen; }
    void set_coreShowerLen(unsigned long int CoreShowerLen) { core_shower_len_ = CoreShowerLen; }
    void set_fractionInCE_E(unsigned int fraction ){ fractionInCE_E_ = fraction; }
    void set_fractionInCoreCE_E(unsigned int fraction ){ fractionInCoreCE_E_ = fraction; }
    void set_fractionInEarlyCE_E(unsigned int fraction ){ fractionInEarlyCE_E_ = fraction; }

    // Getters
    unsigned int clock() const { return clock_; }
    unsigned int index() const { return index_; }
    bool lastFrame() const { return lastFrame_; }
    bool dataValid() const { return dataValid_; }

    binary n_tc() const { return n_tc_; }
    binary e() const { return e_; }
    binary e_em() const { return e_em_; }
    binary e_em_core() const { return e_em_core_; }
    binary e_h_early() const { return e_h_early_; }
    binary w() const { return w_; }
    binary n_tc_w() const { return n_tc_w_; }
    binary w2() const { return w2_; }
    binary wz() const { return wz_; }
    binary wphi() const { return wphi_; }
    binary wroz() const { return wroz_; }
    binary wz2() const { return wz2_; }
    binary wphi2() const { return wphi2_; }
    binary wroz2() const { return wroz2_; }
    binary layerbits() const { return layerbits_; }
    binary sat_tc() const { return sat_tc_; }
    binary shapeq() const { return shapeq_; }

    unsigned long int firstLayer() const { return first_layer_; }
    unsigned long int lastLayer() const { return last_layer_; }
    unsigned long int showerLen() const { return shower_len_; }
    unsigned long int coreShowerLen() const { return core_shower_len_; }
    unsigned int fractionInCE_E() const { return fractionInCE_E_; }
    unsigned int fractionInCoreCE_E() const { return fractionInCoreCE_E_; }
    unsigned int fractionInEarlyCE_E() const { return fractionInEarlyCE_E_; }
    // Operators
    const HGCalCluster& operator+=(const HGCalCluster& hc);

    void saturate();


    // Operators
    bool operator==(const HGCalCluster& other) const;

    static HGCalClusterSAPtrCollection ReadDebugFile( const std::string& aFilename );
    
    // Object in firmware representation
    HGCalCluster_HW& hwCluster() { return hwCluster_; }

    // Format data into firmware representation
    void clearClusterSumWords();
    ClusterSumWords formatClusterSumWords( const l1thgcfirmware::ClusterAlgoConfig& config );

  // private:
    unsigned int clock_;
    unsigned int index_;

    binary n_tc_;
    binary e_;
    binary e_em_;
    binary e_em_core_;
    binary e_h_early_;
    binary w_;
    binary n_tc_w_;
    binary w2_;
    binary wz_;
    binary wphi_;
    binary wroz_;
    binary wz2_;
    binary wphi2_;
    binary wroz2_;
    binary layerbits_;
    binary sat_tc_;
    binary shapeq_;

    unsigned long int first_layer_;
    unsigned long int last_layer_;
    unsigned long int shower_len_;
    unsigned long int core_shower_len_;

    unsigned int fractionInCE_E_;
    unsigned int fractionInCoreCE_E_;
    unsigned int fractionInEarlyCE_E_;
    unsigned int sortKey_ , sortKey2_;

    bool L_ , R_ , X_;
    bool lastFrame_;
    bool dataValid_;

    // Firmware representation of cluster as sent on links to L1T
    HGCalCluster_HW hwCluster_;

    // Firmware representation of cluster sum (input to cluster properties)
    ClusterSumWords packedData_clustersSums_;
  };

}  // namespace l1thgcfirmware

#endif
