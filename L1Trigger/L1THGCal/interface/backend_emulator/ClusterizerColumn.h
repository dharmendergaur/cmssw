#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalTriggerCell_SA.h"
#include <array>
#include <limits.h>

namespace l1thgcfirmware {

  class HGCalTriggerCell;
  struct ClusterizerColumn
  {
    ClusterizerColumn() : Current( nullptr ) , midr( UINT_MAX ) , counter( -1 )
    {}
    
    std::array< l1thgcfirmware::HGCalTriggerCellSAPtr , 216 > MappedTCs;
    std::list< l1thgcfirmware::HGCalHistogramCellSAPtr > MaximaFifo;
    l1thgcfirmware::HGCalHistogramCellSAPtr Current;
    uint32_t midr;
    int32_t counter;
    
    ClusterizerColumn& pop()
    { 
        if( MaximaFifo.size() ) {
        Current = MaximaFifo.front();
        MaximaFifo.pop_front();
        }else{
        Current = nullptr;
        }
        
        if( Current == nullptr or MaximaFifo.empty() ) midr = UINT_MAX;
        else midr = ( Current->Y_ + MaximaFifo.front()->Y_ ) / 2;
        
        return *this;
    }
  
  };
}