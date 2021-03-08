#ifndef __L1Trigger_L1THGCal_HGCalHistoCluteringConfig_SA_h__
#define __L1Trigger_L1THGCal_HGCalHistoCluteringConfig_SA_h__

#include <vector>

namespace l1t {

  class clusterAlgoConfig_SA {
    public:
        clusterAlgoConfig_SA(double kMidRadius,
                            double dr,
                            std::vector<double> dr_byLayer_coefficientA,
                            std::vector<double> dr_byLayer_coefficientB,
                            float ptC3dThreshold)
            : kMidRadius_(kMidRadius),
            dr_(dr),
            dr_byLayer_coefficientA_(dr_byLayer_coefficientA),
            dr_byLayer_coefficientB_(dr_byLayer_coefficientB),
            ptC3dThreshold_(ptC3dThreshold) {}

        void setParameters( const clusterAlgoConfig_SA& newConfig ) {
          kMidRadius_ = newConfig.kMidRadius();
        //   dr_(dr),
        //   dr_byLayer_coefficientA_(dr_byLayer_coefficientA),
        //   dr_byLayer_coefficientB_(dr_byLayer_coefficientB),
        //   ptC3dThreshold_(ptC3dThreshold)
        }
        double kMidRadius() const { return kMidRadius_; }
        double dr() const { return dr_; }
        std::vector<double> dr_byLayer_coefficientA() const { return dr_byLayer_coefficientA_; }
        std::vector<double> dr_byLayer_coefficientB() const { return dr_byLayer_coefficientB_; }
        float ptC3dThreshold() const { return ptC3dThreshold_; }

    private:
        double kMidRadius_;
        double dr_;
        std::vector<double> dr_byLayer_coefficientA_;
        std::vector<double> dr_byLayer_coefficientB_;
        float ptC3dThreshold_;
  };

}  // namespace l1t

#endif
