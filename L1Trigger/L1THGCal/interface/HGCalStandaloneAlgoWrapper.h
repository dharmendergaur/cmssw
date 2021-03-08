#ifndef __L1Trigger_L1THGCal_HGCalStandaloneAlgoWrapper_h__
#define __L1Trigger_L1THGCal_HGCalStandaloneAlgoWrapper_h__

// #include "FWCore/ParameterSet/interface/ParameterSet.h"
// #include "DataFormats/Common/interface/Handle.h"

// #include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"

class HGCalStandaloneAlgoWrapper {
public:
  HGCalStandaloneAlgoWrapper() {}

  virtual ~HGCalStandaloneAlgoWrapper(){}

  virtual void convertCMSSWInputs() const {};
  virtual void convertAlgorithmOutputs() const {};

  // virtual void runAlgorithm(const InputCollection& inputColl, OutputCollection& outColl, const edm::EventSetup& es) = 0;
  virtual void runAlgorithm() const = 0;
};

#endif
