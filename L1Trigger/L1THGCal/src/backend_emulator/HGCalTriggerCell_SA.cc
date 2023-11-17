#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalTriggerCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/Toolkit.h"
  
using namespace l1thgcfirmware;



#define MEMBERS &HGCalTriggerCell::clock_ , \
                &HGCalTriggerCell::index_ , \
                &HGCalTriggerCell::rOverZ_ , \
                &HGCalTriggerCell::layer_  , \
                &HGCalTriggerCell::energy_ , \
                &HGCalTriggerCell::phi_ , \
                &HGCalTriggerCell::sortKey_ , \
                &HGCalTriggerCell::deltaR2_ , \
                &HGCalTriggerCell::dX_ , \
                &HGCalTriggerCell::Y_ , \
                &HGCalTriggerCell::lastFrame_ , \
                &HGCalTriggerCell::dataValid_

bool HGCalTriggerCell::operator==(const HGCalTriggerCell& other) const { 
  return ::check( *this , other , MEMBERS );
}

std::ostream& operator<< ( std::ostream& aStr , const HGCalTriggerCell& aCell )
{
  return ::print( aStr , aCell ,  MEMBERS );
}

HGCalTriggerCellSAPtrCollection HGCalTriggerCell::ReadDebugFile( const std::string& aFilename )
{
  HGCalTriggerCellSAPtrCollection lRet;
  ::ReadDebugFile( aFilename , lRet , MEMBERS );
  return lRet;
}

// bool HGCalTriggerCell::operator==(const HGCalTriggerCell& rhs) const {
//   bool isEqual = (this->index() == rhs.index());
//   return isEqual;
// }

// bool HGCalTriggerCell::operator==(const std::shared_ptr<HGCalTriggerCell>& rhs) const {
//   bool isEqual = (this->index() == rhs->index());
//   return isEqual;
// }

#undef MEMBERS