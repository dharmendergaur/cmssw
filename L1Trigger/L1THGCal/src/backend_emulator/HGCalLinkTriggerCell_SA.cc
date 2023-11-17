#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalLinkTriggerCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/Toolkit.h"
  
using namespace l1thgcfirmware;

#define MEMBERS &HGCalLinkTriggerCell::clock_ , \
                &HGCalLinkTriggerCell::index_ , \
                &HGCalLinkTriggerCell::data_ , \
                &HGCalLinkTriggerCell::lastFrame_ , \
                &HGCalLinkTriggerCell::dataValid_

bool HGCalLinkTriggerCell::operator==(const HGCalLinkTriggerCell& other) const { 
  return ::check( *this , other , MEMBERS );
}

std::ostream& operator<< ( std::ostream& aStr , const HGCalLinkTriggerCell& aCell )
{
  return ::print( aStr , aCell ,  MEMBERS );
}

HGCalLinkTriggerCellSAPtrCollection HGCalLinkTriggerCell::ReadDebugFile( const std::string& aFilename )
{
  HGCalLinkTriggerCellSAPtrCollection lRet;
  ::ReadDebugFile( aFilename , lRet , MEMBERS );
  return lRet;
}

#undef MEMBERS