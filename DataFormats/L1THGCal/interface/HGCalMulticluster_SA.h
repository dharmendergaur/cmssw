#ifndef DataFormats_L1Trigger_HGCalMulticluster_SA_h
#define DataFormats_L1Trigger_HGCalMulticluster_SA_h

#include "DataFormats/L1THGCal/interface/HGCalCluster_SA.h"

#include <vector>

namespace l1t {

  class HGCalMulticluster_SA {
  public:
    HGCalMulticluster_SA()
        : hOverEValid_(false),
          centre_x_(0),
          centre_y_(0),
          centre_z_(0),
          centreProj_x_(0),
          centreProj_y_(0),
          centreProj_z_(0),
          mipPt_(0),
          sumPt_() {}

    HGCalMulticluster_SA(const l1t::HGCalCluster_SA& tc, float fraction = 1.);

    void addConstituent(const l1t::HGCalCluster_SA& tc, bool updateCentre = true, float fraction = 1.);

    ~HGCalMulticluster_SA(){};

    const std::vector<l1t::HGCalCluster_SA> constituents() const { return constituents_; }

    unsigned size() const { return constituents_.size(); }

    float sumPt() const { return sumPt_; }

  private:
    float hOverE_;
    bool hOverEValid_;

    // Could replace this with own simple implementation of GlobalPoint?
    // Or just a struct?
    float centre_x_;
    float centre_y_;
    float centre_z_;

    float centreProj_x_;
    float centreProj_y_;
    float centreProj_z_;

    float mipPt_;
    float sumPt_;

    std::vector<l1t::HGCalCluster_SA> constituents_;

    void updateP4AndPosition(const l1t::HGCalCluster_SA& tc, bool updateCentre = true, float fraction = 1.);
  };

}  // namespace l1t

#endif
