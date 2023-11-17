#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoSeeding_SA.h"

#include <cmath>
#include <memory>

using namespace std;
using namespace l1thgcfirmware;

HGCalHistoSeeding::HGCalHistoSeeding(const ClusterAlgoConfig& config) : config_(config) {}

void HGCalHistoSeeding::runSeeding(const HGCalTriggerCellSAPtrCollection& triggerCellsIn,
                                   HGCalHistogramCellSAPtrCollection& histogramOut) const {
  HGCalHistogramCellSAPtrCollection histoCells;

  triggerCellToHistogramCell(triggerCellsIn, histoCells);

  makeHistogram(histoCells, histogramOut);

  // Smearing
  smearHistogram1D(histogramOut);
  interleaving(histogramOut);
  normalizeArea(histogramOut);
  smearHistogram2D(histogramOut);
  deinterleave(histogramOut);

  //Maxima finding
  maximaFinder1D(histogramOut);
  interleaving(histogramOut);
  maximaFinder2D(histogramOut);
  calculateAveragePosition(histogramOut);
  deinterleave(histogramOut);
  maximaFanout(histogramOut);
}

// Convert each trigger cell into a histogram cell
// Phi bin corresponds to the stream the TC is in, r/z bin is calculated from the TC r/z
void HGCalHistoSeeding::triggerCellToHistogramCell(const HGCalTriggerCellSAPtrCollection& triggerCellsIn,
                                                   HGCalHistogramCellSAPtrCollection& histogramOut) const {
  const unsigned int latency = config_.getStepLatency(TcToHc);

  histogramOut.clear();
  for (auto& tc : triggerCellsIn) {
    auto hc = make_unique<HGCalHistogramCell>(tc->clock() + latency,
                                        tc->index(),
                                        tc->energy(),
                                        tc->phi(),
                                        tc->rOverZ(),
                                        1,
                                        int((tc->rOverZ() - config_.rOverZHistOffset()) / config_.rOverZBinSize()),
                                        tc->lastFrame());
    histogramOut.push_back(std::move(hc));
  }
}


// Sum all input histogram cells (each representing one TC) to produce one histogram
void HGCalHistoSeeding::makeHistogram(const HGCalHistogramCellSAPtrCollection& histogramCells,
                                      HGCalHistogramCellSAPtrCollection& histogramOut) const {
  const unsigned int latency = config_.getLatencyUpToAndIncluding(Hist);

  histogramOut.clear();
  const unsigned latencyOffset = 4;
  for (unsigned int iRow = 0; iRow < config_.cRows(); ++iRow) {
    for (unsigned int iColumn = 0; iColumn < config_.cColumns(); ++iColumn) {
      histogramOut.push_back(make_unique<HGCalHistogramCell>(latencyOffset*iRow + latency, iColumn, iRow, (iRow == config_.cRows()-1) ));
    }
  }

  for (const auto& hc : histogramCells) {
    const unsigned int binIndex = config_.cColumns() * hc->sortKey() + hc->index();
    *histogramOut.at(binIndex) += *hc;
  }
}

// Smearing at constant r/z using exponentially falling 2^-n kernel
// Maximum extent of smearing (how many neighbouring bins included) can vary for each r/z bin, and is specified in the configuration
void HGCalHistoSeeding::smearHistogram1D(HGCalHistogramCellSAPtrCollection& histogram) const {
  const unsigned int stepLatency = config_.getStepLatency( Step::Smearing1D );

  HGCalHistogramCellSAPtrCollection lHistogram;
  lHistogram.reserve( histogram.size() );  

  for ( unsigned int iRow = 0; iRow != config_.cRows(); ++iRow ) {   
    for ( unsigned int iColumn = 0; iColumn != config_.cColumns(); ++iColumn ) {
      const unsigned int binIndex = ( config_.cColumns() * iRow ) + iColumn;
      auto hc = std::make_shared<HGCalHistogramCell>( *histogram[binIndex] );
      hc->clock_ += stepLatency;
      
      // We have repeated integer addition and **integer** division, which is non-distributive
      // Each iteration brings in up to two neighbours on each side, and adjusts the scale-factor accordingly
      int scale( 1 ) , width( config_.kernelWidth( iRow ) ), offset( 1 );
      while ( width > 0 ) {

        unsigned int l1(0),l2(0),r1(0),r2(0);

        if ( width >= 2 ) {
          if ( int(iColumn - offset - 1)  >= 0 ) l2 = histogram[binIndex - offset - 1]->S_/4;
          if ( int(iColumn + offset + 1) <= int(config_.cColumns()-1) ) r2 = histogram[binIndex + offset + 1]->S_/4;
        }
        
        if ( int(iColumn - offset)  >= 0 ) l1 = histogram[binIndex - offset]->S_/2;
        if ( int(iColumn + offset)  <= int(config_.cColumns()-1) ) r1 = histogram[binIndex + offset]->S_/2;
        
        hc->S_ += ( ( l2 + l1 ) / scale + ( r2 + r1 ) / scale );
        scale *= 4;
        width -= 2;
        offset += 2;
      }
     
      lHistogram.push_back( hc );
    }
  }
  
  histogram = lHistogram;
}

void HGCalHistoSeeding::interleaving( HGCalHistogramCellSAPtrCollection& histogram ) const
{
  const unsigned int stepLatency = config_.getStepLatency( Step::Interleaving );

  HGCalHistogramCellSAPtrCollection lHistogram;
  lHistogram.reserve( histogram.size() );  
  
  for ( unsigned int iRow = 0; iRow != config_.cRows(); ++iRow ) {  
    for( unsigned int offset = 0; offset != 4; ++offset ) {
      for ( unsigned int iColumn = 0; iColumn != config_.cColumns(); iColumn+=4 ) {
        const unsigned int binIndex = ( config_.cColumns() * iRow ) + iColumn + offset;
        auto hc = std::make_shared<HGCalHistogramCell>( *histogram[binIndex] );

        hc->clock_ += stepLatency + offset;
        hc->sortKey2_ = hc->index_;
        hc->index_ /= 4;
        hc->lastFrame_ = (offset==3) and (iRow == config_.cRows()-1);
        
        lHistogram.push_back( hc );
      }
    }
  }
  
  histogram = lHistogram;  
}

void HGCalHistoSeeding::normalizeArea(HGCalHistogramCellSAPtrCollection& histogram) const {
  const unsigned int stepLatency = config_.getStepLatency(NormArea);
  for (unsigned int iBin = 0; iBin < histogram.size(); ++iBin) {
    HGCalHistogramCell& hc = *histogram.at(iBin);
    hc.addLatency(stepLatency);
    // Check if 131072 is same as (0x1 << (config_.nBitsAreaNormLUT() + 1))
    hc.setS( int( 1.0 * hc.S() * config_.areaNormalization( hc.sortKey() ) / 131072 ) );
    // hc.setS(int(float(hc.S()) * config_.areaNormalization(hc.sortKey()) / (0x1 << (config_.nBitsAreaNormLUT() + 1))));
  }
}

// Smearing at constant phi with exponentially falling 2^-n kernel
// Limited to +/- 1 bin in phi (only consider nearest neighbour), so adding half of the row above and below
void HGCalHistoSeeding::smearHistogram2D(HGCalHistogramCellSAPtrCollection& histogram) const {
    const unsigned int stepLatency = config_.getStepLatency( Step::Smearing2D );
  
    HGCalHistogramCellSAPtrCollection lHistogram;
    lHistogram.reserve( histogram.size() );  

    for ( unsigned int iRow = 0; iRow != config_.cRows(); ++iRow ) {   
        for ( unsigned int iColumn = 0; iColumn != config_.cColumns(); ++iColumn ) {
        const unsigned int binIndex = ( config_.cColumns() * iRow ) + iColumn;
    
        auto hc = std::make_shared<HGCalHistogramCell>( *histogram[binIndex] );
        hc->clock_ += stepLatency;

        if ( iRow != 0 )                 hc->S_ += histogram[ binIndex - config_.cColumns() ]->S_ / 2;
        if ( iRow != config_.cRows()-1 ) hc->S_ += histogram[ binIndex + config_.cColumns() ]->S_ / 2;
        
        lHistogram.push_back( hc );
        }
    }
    
    histogram = lHistogram;
}

void HGCalHistoSeeding::deinterleave( HGCalHistogramCellSAPtrCollection& histogram ) const
{
  const unsigned int stepLatency = config_.getStepLatency( Step::Deinterleaved );
  
  HGCalHistogramCellSAPtrCollection lHistogram;
  lHistogram.reserve( histogram.size() );  
  
  int cQuarterCols = config_.cColumns() / 4;
  auto cClock0 = histogram[0]->clock_;

  for ( unsigned int iRow = 0; iRow != config_.cRows(); ++iRow ) {   
    for ( unsigned int iColumn = 0; iColumn != config_.cColumns(); ++iColumn ) {
  
      auto binIndex  = ( iRow * config_.cColumns() ) + ( (iColumn % 4) * cQuarterCols ) + int( iColumn / 4 );
      auto binIndex2 = ( iRow * config_.cColumns() ) + ( 3 * cQuarterCols );
      
      auto hc = std::make_shared<HGCalHistogramCell>( *histogram[binIndex] );      
      hc->index_ = iColumn; 
      hc->clock_ = cClock0 + ( 4 * iRow ) + stepLatency;
      hc->lastFrame_ = histogram[binIndex2]->lastFrame_;
      lHistogram.push_back( hc );
    }
  }
  histogram = lHistogram;  
}

void HGCalHistoSeeding::maximaFinder1D( HGCalHistogramCellSAPtrCollection& histogram ) const
{  
  const unsigned int stepLatency = config_.getStepLatency( Step::Maxima1D );

  HGCalHistogramCellSAPtrCollection lHistogram;
  lHistogram.reserve( histogram.size() );  
  for( auto& hc : histogram ) lHistogram.push_back( std::make_shared<HGCalHistogramCell>( *hc ) );

  for ( unsigned int iRow = 0; iRow != config_.cRows(); ++iRow ) {  
    auto width = config_.maximaWidths( iRow ); 
    for ( unsigned int iColumn = 0; iColumn != config_.cColumns(); ++iColumn ) {
      const unsigned int binIndex = ( config_.cColumns() * iRow ) + iColumn;

      auto& a = *lHistogram[binIndex];
      a.clock_ += stepLatency;
      
      if ( a.S_ == 0 ) a.sortKey2_ = 0;
      
      for( int i( -width ); i <= (int)width ; ++i ) {
        if( int(iColumn) + i < 0 ) continue;
        if( iColumn + i >= config_.cColumns() ) continue; 
        const unsigned int binIndex2 = binIndex + i;       
        auto& b = *histogram[ binIndex2 ];
     
        if( a.S_ < b.S_ ){
          a.S_ = b.S_;
          a.X_ = b.X_;
          a.Y_ = b.Y_;
          a.N_ = b.N_;
          a.sortKey2_ = b.sortKey2_;
          a.maximaOffset_ = i;          
        }        
      }
      
    }
  }
  
  histogram = lHistogram;  
}

void HGCalHistoSeeding::maximaFinder2D( HGCalHistogramCellSAPtrCollection& histogram ) const
{
  const unsigned int stepLatency = config_.getStepLatency( Step::Maxima2D );
  
  HGCalHistogramCellSAPtrCollection lHistogram;
  lHistogram.reserve( histogram.size() );  
  for( auto& hc : histogram ) lHistogram.push_back( std::make_shared<HGCalHistogramCell>( *hc ) );

  for ( unsigned int iRow = 0; iRow != config_.cRows(); ++iRow ) {   
    for ( unsigned int iColumn = 0; iColumn != config_.cColumns(); ++iColumn ) {
      const unsigned int binIndex = ( config_.cColumns() * iRow ) + iColumn;
  
      auto& i = *lHistogram[binIndex];

      bool criteria0( true ) , criteria1( i.maximaOffset_ == 0 ) , criteria2( false ) , criteria3( i.S_ > config_.thresholdMaxima( i.sortKey_ ) );

      if ( iRow > 0 ) {
        auto& lRef = *histogram[ binIndex - config_.cColumns() ];
        criteria0 = ( i.S_ > lRef.S_ ) or ( ( i.S_ == lRef.S_ ) ) ; // and ( lRef.maximaOffset_ >= 0 ) ); //maximaOffset_ is unsigned
      }

      if ( iRow < (config_.cRows()-1) ) {
        auto& lRef = *histogram[ binIndex + config_.cColumns() ];
        criteria2 = ( i.S_ > lRef.S_ ) or ( ( i.S_ == lRef.S_ ) ) ; // and ( lRef.maximaOffset_ >= 0 ) ); //maximaOffset_ is unsigned
      }
      
      i.clock_ += stepLatency;
      i.maximaOffset_ = ( criteria0 and criteria1 and criteria2 and criteria3 );    
    }
  }

  histogram = lHistogram;   
}

void HGCalHistoSeeding::calculateAveragePosition(HGCalHistogramCellSAPtrCollection& histogram) const {
  const unsigned int stepLatency = config_.getStepLatency(CalcAverage);
  for (auto& hc : histogram) {
    hc->addLatency(stepLatency);
    const unsigned nBits = 17;
    if (hc->N() > 0) {
      unsigned int inv_N = int(round(float(0x1 << nBits) / hc->N()));
      hc->setX((hc->X() * inv_N) >> nBits);
      hc->setY((hc->Y() * inv_N) >> nBits);
    }
  }
}

void HGCalHistoSeeding::maximaFanout( HGCalHistogramCellSAPtrCollection& histogram ) const
{
  const unsigned int stepLatency = config_.getStepLatency( Step::MaximaFanout );
  const std::vector< int > Sign = { -1 , 1 };
  
  
  for ( unsigned int iRow = 0; iRow != config_.cRows(); ++iRow ) {   
    int width = config_.fanoutWidths( iRow ); 
    
    std::vector< unsigned int > lCols;
    
    for ( unsigned int iColumn = 0; iColumn != config_.cColumns(); ++iColumn ) {
      const unsigned int binIndex = ( config_.cColumns() * iRow ) + iColumn;
      auto& i = *histogram[binIndex];
      i.clock_ += stepLatency;
      if( i.maximaOffset_ == 1 ){
        i.maximaOffset_ = width;      
        i.left_ = i.right_ = true;
        lCols.push_back( iColumn );
      }
    }

    for ( int x = 1; x <= width; ++x ) {
      for( auto& col : lCols ) {
        for( auto& sign : Sign ) {
          auto iColumn = int(col) + (sign*x);

          if( iColumn < 0 ) continue;
          if( iColumn >= int(config_.cColumns()) ) continue;
          
          const unsigned int binIndex = ( config_.cColumns() * iRow ) + iColumn;

          auto& i = *histogram[binIndex];

          if( i.maximaOffset_ == 0 )
          {         
            HGCalHistogramCell* l = ( iColumn>0                    ? &*histogram[binIndex-1] : NULL );
            HGCalHistogramCell* r = ( iColumn<int(config_.cColumns()-1) ? &*histogram[binIndex+1] : NULL );        
            HGCalHistogramCell* other( NULL );          
            
            if( l and l->maximaOffset_ and r and r->maximaOffset_ )
              if( l->S_ > r->S_ ) other = l;
              else other = r;
            else if( l and l->maximaOffset_ ) other = l;
            else if( r and r->maximaOffset_ ) other = r;
            
            if( other )
            {
              i.S_ = other->S_;
              i.X_ = other->X_;
              i.Y_ = other->Y_;
              i.N_ = other->N_;
              i.sortKey_ = other->sortKey_;
              i.sortKey2_ = other->sortKey2_;
              i.maximaOffset_ = other->maximaOffset_ - 1;
              if( other == r ){
                i.left_ = other->left_;
                other->left_ = false;
              }
              if( other == l ){
                i.right_ = other->right_;            
                other->right_ = false;            
              }        
            }
            
          }    
        }
      }
    }

    
  }
  
}