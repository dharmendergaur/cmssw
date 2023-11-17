#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClustering_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/ClusterizerColumn.h"

#include <cmath>
#include <algorithm>
#include <memory>

using namespace std;
using namespace l1thgcfirmware;

HGCalHistoClustering::HGCalHistoClustering(const ClusterAlgoConfig& config) : config_(config) {}

void HGCalHistoClustering::runClustering(HGCalTriggerCellSAPtrCollection& triggerCells,
                                         HGCalHistogramCellSAPtrCollection& histogram,
                                         HGCalClusterSAPtrCollection& clustersOut ) const {

  HGCalTriggerCellSAPtrCollection triggerCellRamOut;
  HGCalTriggerCellSAPtrCollection clusteredTriggerCellsOut;
  HGCalHistogramCellSAPtrCollection maximaFifoOut;
  clusterizer(triggerCells, histogram, triggerCellRamOut, maximaFifoOut);
  triggerCellToCluster(triggerCells, histogram, clustersOut);
  clusterAccumulator(clustersOut, histogram);
  clusterTree(clustersOut);
}

void HGCalHistoClustering::clusterizer( HGCalTriggerCellSAPtrCollection& triggerCells, HGCalHistogramCellSAPtrCollection& histogram, HGCalTriggerCellSAPtrCollection& triggerCellsRamOut, HGCalHistogramCellSAPtrCollection& maximaFifoOut ) const
{
  std::array< ClusterizerColumn , 124 > lColumns;

  // Map the TCs into the RAM using the LUT
  auto start = triggerCells.front()->clock();
  for( auto& tc : triggerCells ){
    auto frame = tc->clock() - start;
    auto& lut_out = config_.TriggerCellAddressLUT( ( 216*tc->index() ) + frame );
    lColumns.at( tc->index() ).MappedTCs.at( lut_out ) = std::move(tc);
  }
  triggerCells.clear();
  
  // Map the maxima into the FIFO
  for( auto& i : histogram ){
    if( i->maximaOffset_>0 or i->left_ or i->right_ ) lColumns.at( i->index_ ).MaximaFifo.push_back( i );      
  }
  histogram.clear();

  // Move the first two entries out of the FIFO into "current" and "next"
  // for( auto& lColumn : lColumns ) lColumn.pop().pop();

  HGCalTriggerCellSAPtrCollection triggerCellsOut;
  HGCalHistogramCellSAPtrCollection histogramOut;

  // Read the tcs out sequentially
  for ( unsigned int frame = 0; frame != 216; ++frame ) {  
    for ( unsigned int iColumn = 0; iColumn != config_.cColumns(); ++iColumn ) {      
      auto& col = lColumns.at( iColumn );
      auto& tc = col.MappedTCs.at( frame );
      
                                                        
      // Get the maxima from the FIFO                                                //  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< STILL TRYING TO EMULATE HOW THE FIRMWARE DOES THIS!
      if( col.counter < 0 )
      {
        auto MaximaInFifo = ( ( col.MaximaFifo.size() ) and ( col.MaximaFifo.front()->clock_-306 < frame ) );
        auto NoCurrent     = ( col.Current == nullptr );
        auto TcPastMidpoint = ( tc != nullptr ) and ( tc->rOverZ() > col.midr ); // and ( col.MaximaFifo.front()->clock_-303 < frame ) ) );
        if ( MaximaInFifo and ( NoCurrent or TcPastMidpoint ) ) col.counter = 3;
      }
      else    
      {  
        if( col.counter == 0 ) col.pop();
        if( col.counter >= 0 ) col.counter--;
      }


      if ( col.Current != nullptr )
      {
        // std::cout << *col.Current;        
        maximaFifoOut.push_back( std::make_unique< HGCalHistogramCell >( *col.Current ) );
        auto& hcx = maximaFifoOut.back();
        hcx->clock_ = frame + 289 + 16;     
      }
    
      // Compare the TC against the maxima
      if( tc == nullptr ) continue;
      
      triggerCellsRamOut.push_back( std::make_unique< HGCalTriggerCell >( *tc ) );
      auto& tcx = triggerCellsRamOut.back();
      tcx->setClock( frame + 289 + 20 );
      tcx->setLastFrame( false ); //(frame==215);

      // unsigned int CurrentdR2Cut(5000); // Magic numbers
      // unsigned int CurrentDist = Dist( tc , col.Current , config_ );
      // if( CurrentDist < CurrentdR2Cut )
      // {
        // HGCalHistogramCellSAPtr hc = std::make_shared< HGCalHistogramCell > ( *col.Current );
        // hc->clock_ = frame + 289 + 20;    
        // histogramOut.emplace_back( hc );        

        // tc->clock_ = frame + 289 + 20; 
        // tc->sortKey_ = hc->sortKey_;        
        // triggerCellsOut.push_back( tc );
      // }

      
    }
  }
  triggerCells = move(triggerCellsOut);
  histogram = move(histogramOut);
}

void HGCalHistoClustering::triggerCellToCluster(const HGCalTriggerCellSAPtrCollection& clusteredTriggerCells,
                                                const HGCalHistogramCellSAPtrCollection& histogram,
                                                HGCalClusterSAPtrCollection& clustersOut) const {

    std::map< std::pair< int , int > , HGCalHistogramCellSAPtr > maxima_map;
    for( auto& x : histogram ) maxima_map[ std::make_pair( x->index_ , x->clock_ ) ] = x;
    

    const unsigned int stepLatency = config_.getStepLatency( TriggerCellToCluster );

    for ( const auto& tc : clusteredTriggerCells ) {

        auto cluster = make_unique<HGCalCluster>(
                            tc->clock() + stepLatency,
                            tc->index(),
                            tc->lastFrame(),
                            tc->dataValid()
        );

        auto lIt1( maxima_map.find( std::make_pair( tc->index() , tc->clock() ) ) );
        auto lIt2( maxima_map.find( std::make_pair( tc->index() , tc->clock()+1 ) ) );    
        cluster->L_ = lIt1->second->left_;
        cluster->R_ = lIt1->second->right_;
        cluster->X_ = ( lIt2 == maxima_map.end() ) or ( lIt1->second->sortKey_ != lIt2->second->sortKey_ );
        
        cluster->sortKey_ = lIt1->second->sortKey_;
        cluster->sortKey2_ = lIt1->second->sortKey2_;


        // Cluster from single TC
        // Does this ever happen?
        // if ( tc->deltaR2_ >= 25000 ) { // Magic numbers
        // clusters.push_back( cluster );
        // continue;
        // }

        uint64_t s_TC_W = ( int( tc->energy() / 4 ) == 0 ) ? 1 : tc->energy() / 4;
        uint64_t s_TC_Z = config_.depth( tc->layer() );

        unsigned int triggerLayer = config_.triggerLayer( tc->layer() );
        
        unsigned int s_E_EM = ( (  ( (uint64_t) tc->energy() * config_.layerWeight_E_EM( tc->layer() ) ) + config_.correction() ) >> 18 );
        if ( s_E_EM > config_.saturation() ) s_E_EM = config_.saturation();

        unsigned int s_E_EM_core = ( ( (uint64_t) tc->energy() * config_.layerWeight_E_EM_core( tc->layer() ) + config_.correction() ) >> 18 );
        if ( s_E_EM_core > config_.saturation() ) s_E_EM_core = config_.saturation();

        // Alternative constructor perhaps?
        cluster->set_n_tc(1);
        cluster->set_n_tc_w(1);

        cluster->set_e((config_.layerWeight_E(triggerLayer) == 1) ? tc->energy() : 0);
        cluster->set_e_h_early((config_.layerWeight_E_H_early(triggerLayer) == 1) ? tc->energy() : 0);

        cluster->set_e_em(s_E_EM);
        cluster->set_e_em_core(s_E_EM_core);

        cluster->set_w(s_TC_W);
        cluster->set_w2(s_TC_W * s_TC_W);

        cluster->set_wz(s_TC_W * s_TC_Z);
        cluster->set_wphi(s_TC_W * tc->phi());
        cluster->set_wroz(s_TC_W * tc->rOverZ());

        cluster->set_wz2(s_TC_W * s_TC_Z * s_TC_Z);
        cluster->set_wphi2(s_TC_W * tc->phi() * tc->phi());
        cluster->set_wroz2(s_TC_W * tc->rOverZ() * tc->rOverZ());

        const unsigned nTriggerLayers = 34;  // Should get from config/elsewhere in CMSSW
        cluster->set_layerbits( ( ( (uint64_t) 1) << ( nTriggerLayers - triggerLayer ) ) & 0x3FFFFFFFF);
        cluster->set_sat_tc(cluster->e() == config_.saturation() || cluster->e_em() == config_.saturation());

        cluster->set_shapeq(1);

        clustersOut.push_back( move(cluster) );
    }
}


void HGCalHistoClustering::clusterAccumulator( HGCalClusterSAPtrCollection& clusters, const HGCalHistogramCellSAPtrCollection& histogram ) const
{  
  HGCalClusterSAShrPtrCollection output;
  
  std::map< std::pair< int , int > , HGCalClusterSAShrPtr > cluster_map;

  for( auto& x : clusters ){
    auto lKey = std::make_pair( x->sortKey_ , x->index_ );
    auto lIt = cluster_map.find( lKey );
    if ( lIt == cluster_map.end() ){
      HGCalClusterSAShrPtr lVal = make_shared< HGCalCluster >( *x );
      lVal->X_ = true; // Last entry should always have X_ set
      output.push_back( lVal );
      cluster_map[lKey] = lVal;
    } else {
      *lIt->second += *x;
      lIt->second->L_ = x->L_;
      lIt->second->R_ = x->R_;
      lIt->second->X_ = x->X_;
      lIt->second->sortKey_ = x->sortKey_;
      lIt->second->sortKey2_ = x->sortKey2_;
    }
  }

  for( auto& x : histogram ){  
    auto lIt = cluster_map.find( std::make_pair( x->sortKey_ , x->index_ ) );
    if ( lIt != cluster_map.end() ) lIt->second->clock_ = x->clock_ + 11;     
  }

  for( auto& x : output ) x->saturate();
 
  std::sort( output.begin() , output.end() , []( const HGCalClusterSAShrPtr& a , const HGCalClusterSAShrPtr& b ){ return std::make_pair( a->clock_ , a->index_ ) < std::make_pair( b->clock_ , b->index_ ); } );

  clusters.reserve(output.size());
  for (auto& sharedPtr : output) {
      clusters.push_back(std::make_unique<HGCalCluster>(*sharedPtr));
  }
}

void HGCalHistoClustering::clusterTree( HGCalClusterSAPtrCollection& clusters ) const
{
  HGCalClusterSAShrPtrCollection output;


  // vvvvvvvvvvvvvvvvvv HACK TO VERIFY VALUES
  std::map< std::pair< int , int > , HGCalClusterSAShrPtr > cluster_map;

  for( auto& x : clusters ){
    auto lKey = std::make_pair( x->sortKey_ , x->sortKey2_ );
    auto lIt = cluster_map.find( lKey );
    if ( lIt == cluster_map.end() ){
      auto lVal = make_shared< HGCalCluster >( *x );
      lVal->index_ = 0;
      lVal->X_ = 0;
      output.push_back( lVal );
      cluster_map[lKey] = lVal;
    } else {
      *lIt->second += *x;
      lIt->second->clock_ = max( lIt->second->clock_ , x->clock_ );      
      lIt->second->L_ |= x->L_;
      lIt->second->R_ |= x->R_;
      // lIt->second->X_ = x->X_;
      lIt->second->sortKey_ = x->sortKey_;
      lIt->second->sortKey2_ = x->sortKey2_;
    }
  }
  // ^^^^^^^^^^^^^^^^^ HACK TO VERIFY VALUES
  
  for( auto& x : output ){
    x->saturate();
    x->clock_ += 9;
  }

  std::sort( output.begin() , output.end() , []( const HGCalClusterSAShrPtr& a , const HGCalClusterSAShrPtr& b ){ return std::make_pair( a->clock_ , a->index_ ) < std::make_pair( b->clock_ , b->index_ ); } );

  HGCalClusterSAShrPtr last = nullptr;
  for( auto& x : output ){
    if( last != nullptr and x->clock_ <= last->clock_ ) x->clock_ = last->clock_ + 1;
    last = x;
  }

  std::sort( output.begin() , output.end() , []( const HGCalClusterSAShrPtr& a , const HGCalClusterSAShrPtr& b ){ return std::make_pair( a->clock_ , a->index_ ) < std::make_pair( b->clock_ , b->index_ ); } );

  clusters.reserve(output.size());
  for (auto& sharedPtr : output) {
      clusters.push_back(std::make_unique<HGCalCluster>(*sharedPtr));
  }
}
