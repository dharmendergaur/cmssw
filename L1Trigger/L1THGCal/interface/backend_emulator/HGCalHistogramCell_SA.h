#ifndef L1Trigger_L1THGCal_HGCalHistogramCell_SA_h
#define L1Trigger_L1THGCal_HGCalHistogramCell_SA_h

#include <vector>
#include <memory>
#include <iostream>

namespace l1thgcfirmware {
  class HGCalHistogramCell;
  // typedef std::unique_ptr<HGCalHistogramCell> HGCalHistogramCellSAPtr;  
  typedef std::vector<HGCalHistogramCell> HGCalHistogramCellSACollection;
  // typedef std::vector<HGCalHistogramCellSAPtr> HGCalHistogramCellSAPtrCollection;
  typedef std::shared_ptr<HGCalHistogramCell> HGCalHistogramCellSAPtr;  
  typedef std::vector<HGCalHistogramCellSAPtr> HGCalHistogramCellSAPtrCollection;
}

std::ostream& operator<< ( std::ostream& aStr , const l1thgcfirmware::HGCalHistogramCell& aCell );

namespace l1thgcfirmware {

class HGCalHistogramCell {
  public:
    HGCalHistogramCell( unsigned int clock,
                        unsigned int index,
                        unsigned int S,
                        unsigned int X,
                        unsigned int Y,
                        unsigned int N,
                        unsigned int sortKey,
                        bool lastFrame )
                        : clock_(clock),
                          index_(index),
                          S_(S),
                          X_(X),
                          Y_(Y),
                          N_(N),
                          sortKey_(sortKey),
                          lastFrame_(lastFrame),
                          dataValid_(true) {}

    HGCalHistogramCell() : clock_(0),
                           index_(0),
                           S_(0),
                           X_(0),
                           Y_(0),
                           N_(0),
                           sortKey_(0),
                           lastFrame_(false),
                           dataValid_(false) {}

    HGCalHistogramCell( unsigned int clock,
                        unsigned int index,
                        unsigned int sortKey,
                        bool lastFrame )
                        : HGCalHistogramCell ( clock, index, 0, 0, 0, 0, sortKey, lastFrame) {}

    HGCalHistogramCell( const HGCalHistogramCell& aOther /*!< Anonymous argument */ ) = default;
    HGCalHistogramCell& operator = (const HGCalHistogramCell& aOther /*!< Anonymous argument */ ) = default;
    HGCalHistogramCell( HGCalHistogramCell&& aOther /*!< Anonymous argument */ ) = default;
    HGCalHistogramCell& operator = ( HGCalHistogramCell&& aOther /*!< Anonymous argument */ ) = default;

    ~HGCalHistogramCell(){}

    // Setters
    void setClock( const unsigned int clock ) { clock_ = clock; }
    void addLatency( const unsigned int latency ) { clock_ += latency; }
    void setIndex( const unsigned int index ) { index_ = index; }
    void setSortKey( const unsigned int sortKey ) { sortKey_ = sortKey; }
    void setS( const unsigned int S ) { S_ = S; }
    void setX( const unsigned int X ) { X_ = X; }
    void setY( const unsigned int Y ) { Y_ = Y; }
    void setN( const unsigned int N ) { N_ = N; }

    // Getters
    unsigned int clock() const { return clock_; }
    unsigned int index() const { return index_; }
    unsigned int S() const { return S_; }
    unsigned int X() const { return X_; }
    unsigned int Y() const { return Y_; }
    unsigned int N() const { return N_; }
    unsigned int sortKey() const { return sortKey_; }
    bool lastFrame() const { return lastFrame_; }
    bool dataValid() const { return dataValid_; }

    // Operators
    const HGCalHistogramCell& operator+=(const HGCalHistogramCell& aOther);
    // const HGCalHistogramCell operator/(const unsigned int factor) const;
    // const HGCalHistogramCell operator+(const HGCalHistogramCell& aOther) const;
    // const HGCalHistogramCell& operator*=(const unsigned int factor);

    // Operators
    bool operator==(const HGCalHistogramCell& other) const;

    static HGCalHistogramCellSAPtrCollection ReadDebugFile( const std::string& aFilename );

  // private:
    // friend std::ostream& ::operator<< ( std::ostream& aStr , const HGCalHistogramCell& aCell );
  
    unsigned int clock_;
    unsigned int index_;
    unsigned int S_;
    unsigned int X_;
    unsigned int Y_;
    unsigned int N_;
    unsigned int sortKey_;
    unsigned int sortKey2_;
    unsigned int maximaOffset_;
    bool left_;
    bool right_;
    
    bool lastFrame_;
    bool dataValid_;

  };

}  // namespace l1thgcfirmware

#endif
