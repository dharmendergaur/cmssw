#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistogramCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/Toolkit.h"

using namespace l1thgcfirmware;

const HGCalHistogramCell& HGCalHistogramCell::operator+=(const HGCalHistogramCell& aOther) {
  S_ += aOther.S_;
  X_ += aOther.X_;
  Y_ += aOther.Y_;
  N_ += aOther.N_;

  return *this;
}

// const HGCalHistogramCell HGCalHistogramCell::operator/(const unsigned int factor) const {
  
  // HGCalHistogramCell aOther(*this);

  // aOther.S_ /= factor;
  // aOther.X_ = 0;
  // aOther.Y_ = 0;
  // aOther.N_ = 0;
  // return aOther;
// }

// const HGCalHistogramCell HGCalHistogramCell::operator+(const HGCalHistogramCell& aOther) const {
  
  // HGCalHistogramCell sum(*this);
  // sum += aOther;
  // return sum;
// }

// const HGCalHistogramCell& HGCalHistogramCell::operator*=(const unsigned int factor) {
  // S_ = int( 1.0*S_ * factor / 262144 ); // Magic numbers

  // return *this;
// }



#define MEMBERS &HGCalHistogramCell::clock_ , \
                &HGCalHistogramCell::index_ , \
                &HGCalHistogramCell::S_ , \
                &HGCalHistogramCell::X_  , \
                &HGCalHistogramCell::Y_ , \
                &HGCalHistogramCell::N_ , \
                &HGCalHistogramCell::sortKey_ , \
                &HGCalHistogramCell::sortKey2_ , \
                &HGCalHistogramCell::maximaOffset_ , \
                &HGCalHistogramCell::left_ , \
                &HGCalHistogramCell::right_ , \
                &HGCalHistogramCell::lastFrame_ , \
                &HGCalHistogramCell::dataValid_

bool HGCalHistogramCell::operator==(const HGCalHistogramCell& other) const { 
  return ::check( *this , other , MEMBERS );
}

std::ostream& operator<< ( std::ostream& aStr , const HGCalHistogramCell& aCell )
{
  return ::print( aStr , aCell ,  MEMBERS );
}

HGCalHistogramCellSAPtrCollection HGCalHistogramCell::ReadDebugFile( const std::string& aFilename )
{
  HGCalHistogramCellSAPtrCollection lRet;
  ::ReadDebugFile( aFilename , lRet , MEMBERS );
  return lRet;
}

#undef MEMBERS