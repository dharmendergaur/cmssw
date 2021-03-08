#ifndef DataFormats_L1Trigger_HGCalCluster_SA_h
#define DataFormats_L1Trigger_HGCalCluster_SA_h

#include <vector>

namespace l1t {

  class HGCalCluster_SA {
  public:
    HGCalCluster_SA(float x,
                    float y,
                    int zside,
                    unsigned int layer,
                    float eta,
                    float phi,
                    float pt,
                    float mipPt,
                    unsigned int index_cmssw)
        : x_(x),
          y_(y),
          zside_(zside),
          layer_(layer),
          eta_(eta),
          phi_(phi),
          pt_(pt),
          mipPt_(mipPt),
          index_cmssw_(index_cmssw) {}

    ~HGCalCluster_SA(){};

    float x() const { return x_; }
    float y() const { return y_; }
    float zside() const { return zside_; }
    unsigned int layer() const { return layer_; }
    float eta() const { return eta_; }
    float phi() const { return phi_; }
    float pt() const { return pt_; }
    float mipPt() const { return mipPt_; }
    unsigned int index_cmssw() const { return index_cmssw_; }

  private:
    float x_;
    float y_;
    int zside_;
    unsigned int layer_;
    float eta_;
    float phi_;
    float pt_;
    float mipPt_;
    unsigned int index_cmssw_;
  };

  typedef std::vector<HGCalCluster_SA> HGCalClusterSACollection;

}  // namespace l1t

#endif
