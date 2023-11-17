#ifndef L1Trigger_L1THGCal_HGCalLinkTriggerCell_SA_h
#define L1Trigger_L1THGCal_HGCalLinkTriggerCell_SA_h

#include "L1Trigger/L1THGCal/interface/backend_emulator/Binary.h"

#include <vector>
#include <memory>



// Forward declare for the sake of our friends
namespace l1thgcfirmware {
  class HGCalLinkTriggerCell;
  // typedef std::vector<HGCalLinkTriggerCell> HGCalLinkTriggerCellSACollection;
  typedef std::unique_ptr<HGCalLinkTriggerCell> HGCalLinkTriggerCellSAPtr;
  typedef std::vector<HGCalLinkTriggerCellSAPtr> HGCalLinkTriggerCellSAPtrCollection;
  // typedef std::vector< std::vector<std::shared_ptr<HGCalLinkTriggerCell> > > HGCalLinkTriggerCellSAPtrCollections;
  // typedef std::vector< std::vector< std::vector<std::shared_ptr<HGCalLinkTriggerCell> > > > HGCalLinkTriggerCellSAPtrCollectionss;
}

std::ostream& operator<< ( std::ostream& aStr , const l1thgcfirmware::HGCalLinkTriggerCell& aCell );


namespace l1thgcfirmware {

  class HGCalLinkTriggerCell {
  public:
    HGCalLinkTriggerCell() : clock_(0),
                        index_(0),
                        data_(0),
                        lastFrame_(false),
                        dataValid_(false) {}

    HGCalLinkTriggerCell( const HGCalLinkTriggerCell& aOther /*!< Anonymous argument */ ) = default;
    HGCalLinkTriggerCell& operator = (const HGCalLinkTriggerCell& aOther /*!< Anonymous argument */ ) = default;
    HGCalLinkTriggerCell( HGCalLinkTriggerCell&& aOther /*!< Anonymous argument */ ) = default;
    HGCalLinkTriggerCell& operator = ( HGCalLinkTriggerCell&& aOther /*!< Anonymous argument */ ) = default;
    
    
          // lastFrame_(lastFrame),
          // dataValid_(dataValid) {}

    ~HGCalLinkTriggerCell(){}

    // // Setters
    // bool lastFrame() const { return lastFrame_; }
    // bool dataValid() const { return dataValid_; }

    // Operators
    bool operator==(const HGCalLinkTriggerCell& other) const;

    static HGCalLinkTriggerCellSAPtrCollection ReadDebugFile( const std::string& aFilename );

  // private:
    // friend std::ostream& ::operator<< ( std::ostream& aStr , const HGCalLinkTriggerCell& aCell );
 
    unsigned int clock_;
    unsigned int index_; 
    binary data_;
    bool lastFrame_;
    bool dataValid_;

  };
  
}  // namespace l1thgcfirmware

  
#endif
