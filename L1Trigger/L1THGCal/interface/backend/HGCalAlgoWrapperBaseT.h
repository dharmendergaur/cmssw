#ifndef __L1Trigger_L1THGCal_HGCalAlgoWrapperBaseT_h__
#define __L1Trigger_L1THGCal_HGCalAlgoWrapperBaseT_h__

template <typename InputCollection, typename OutputCollection, typename Tparam>
class HGCalAlgoWrapperBaseT {
public:
  HGCalAlgoWrapperBaseT() {}

  virtual ~HGCalAlgoWrapperBaseT(){}

  virtual void configure(const Tparam& parameters) = 0;
  virtual void process(const InputCollection& inputCollection, OutputCollection& outputCollection) const = 0;
};

#endif
