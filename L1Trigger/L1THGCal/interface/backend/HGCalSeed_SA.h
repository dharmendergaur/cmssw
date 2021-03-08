#ifndef DataFormats_L1Trigger_HGCalSeed_SA_h
#define DataFormats_L1Trigger_HGCalSeed_SA_h

namespace l1t {

  class HGCalSeed_SA {
  public:
    HGCalSeed_SA(float x, float y, float z, float energy) : x_(x), y_(y), z_(z), energy_(energy) {}

    ~HGCalSeed_SA(){};

    float x() const { return x_; }
    float y() const { return y_; }
    float z() const { return z_; }
    float energy() const { return energy_; }

  private:
    float x_;
    float y_;
    float z_;
    float energy_;
  };

  typedef std::vector<HGCalSeed_SA> HGCalSeedSACollection;

}  // namespace l1t

#endif
