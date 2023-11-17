#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringDebug_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/Toolkit.h"

#include <iostream>


using namespace std;
using namespace l1thgcfirmware;

// =======================================================================================================================================================
HGCalHistoClusteringDebugSA::HGCalHistoClusteringDebugSA( ClusterAlgoConfig& config ) : HGCalHistoClusteringImplSA( config ) , cHistSize(0) , cOffset( 7316 ) , DebugIntermediate( true )
{}
// =======================================================================================================================================================

void HGCalHistoClusteringDebugSA::runAlgorithm( const l1thgcfirmware::HGCalLinkTriggerCellSAPtrCollection& linkData )
{
  std::cout << "---------------------------------------" << std::endl;
  config_.printConfiguration();
  std::cout << "---------------------------------------" << std::endl;

  // HACK - Use the HGCalLinkTriggerCell type to hold the full link data
  auto LinksInData = HGCalLinkTriggerCell::ReadDebugFile( "../build/HgcStage2Top/top/DebugFiles/LinksIn.txt" );

  // ===================================
  // Unpack TC from links
  // ===================================

  // linkUnpacking_.runLinkUnpacking( LinksInData, unpackedTCs);
  HGCalTriggerCellSAPtrCollection unpackedTCs;
  HGCalLinkTriggerCellSAPtrCollection triggerCellsWork, towersWork;
  linkUnpacking_.unpackLinks( LinksInData, triggerCellsWork, towersWork);
  if( DebugIntermediate ) check( "UnpackLinks.TriggerCells" , triggerCellsWork , HGCalLinkTriggerCell::ReadDebugFile( "../build/HgcStage2Top/top/DebugFiles/UnpackLinks.TriggerCells.txt" ) );  

  unpackedTCs = linkUnpacking_.triggerCellDistribution( triggerCellsWork );
  if( DebugIntermediate ) check( "TriggerCellDistribution" , unpackedTCs , HGCalTriggerCell::ReadDebugFile( "../build/HgcStage2Top/top/DebugFiles/TriggerCellDistribution.txt" ) );

  linkUnpacking_.unpackTriggerCells( unpackedTCs );
  if( DebugIntermediate ) check( "UnpackTriggerCells" , unpackedTCs , HGCalTriggerCell::ReadDebugFile( "../build/HgcStage2Top/top/DebugFiles/UnpackTriggerCells.txt" ) );  

  // ===================================
  // Histogram and seed finding
  // ===================================

  HGCalHistogramCellSAPtrCollection histogram;
  HGCalHistogramCellSAPtrCollection histogramCells;
  seeding_.triggerCellToHistogramCell(unpackedTCs, histogramCells);
  if( DebugIntermediate ) check( "HistogramCell" , histogramCells , HGCalHistogramCell::ReadDebugFile( "../build/HgcStage2Top/top/DebugFiles/TriggerCellToHistogramCell.txt" ) );

  seeding_.makeHistogram(histogramCells, histogram);
  cHistSize = histogram.size();
  if( DebugIntermediate ) check( "Histogram" , histogram , resize( HGCalHistogramCell::ReadDebugFile( "../build/HgcStage2Top/top/DebugFiles/Histogram.txt" ) , cHistSize , cHistSize ) );

  // Smearing
  seeding_.smearHistogram1D(histogram);
  if( DebugIntermediate ) check( "Smearing1D" , histogram , resize( HGCalHistogramCell::ReadDebugFile( "../build/HgcStage2Top/top/DebugFiles/ExponentialSmearingKernel1D.txt" ) , cOffset , cHistSize ) );  

  seeding_.interleaving(histogram);
  if( DebugIntermediate ) check( "Interleaving1" , histogram , resize( HGCalHistogramCell::ReadDebugFile( "../build/HgcStage2Top/top/DebugFiles/Interleaving.txt" ) , cOffset , cHistSize ) ); 

  seeding_.normalizeArea(histogram);
  if( DebugIntermediate ) check( "AreaNormalization" , histogram , resize( HGCalHistogramCell::ReadDebugFile( "../build/HgcStage2Top/top/DebugFiles/AreaNormalization.txt" ) , cOffset , cHistSize ) );

  seeding_.smearHistogram2D(histogram);
  if( DebugIntermediate ) check( "Smearing2D" , histogram , resize( HGCalHistogramCell::ReadDebugFile( "../build/HgcStage2Top/top/DebugFiles/ExponentialSmearingKernel2D.txt" ) , cOffset , cHistSize ) );

  seeding_.deinterleave(histogram);
  if( DebugIntermediate ) check( "Deinterleaved" , histogram , resize( HGCalHistogramCell::ReadDebugFile( "../build/HgcStage2Top/top/DebugFiles/MaximaFinder1D-Deinterlaced.txt" ) , cOffset , cHistSize ) );

  //Maxima finding
  seeding_.maximaFinder1D(histogram);
  if( DebugIntermediate ) check( "Maxima1D" , histogram , resize( HGCalHistogramCell::ReadDebugFile( "../build/HgcStage2Top/top/DebugFiles/MaximaFinder1D.txt" ) , cOffset , cHistSize ) );

  seeding_.interleaving(histogram);
  if( DebugIntermediate ) check( "Interleaving2" , histogram , resize( HGCalHistogramCell::ReadDebugFile( "../build/HgcStage2Top/top/DebugFiles/Interleaving2.txt" ) , cOffset , cHistSize ) );  

  seeding_.maximaFinder2D(histogram);
  if( DebugIntermediate ) check( "Maxima2D" , histogram , resize( HGCalHistogramCell::ReadDebugFile( "../build/HgcStage2Top/top/DebugFiles/MaximaFinder2D.txt" ) , cOffset , cHistSize ) );

  seeding_.calculateAveragePosition(histogram);
  if( DebugIntermediate ) check( "CalculateAveragePositions" , histogram , resize( HGCalHistogramCell::ReadDebugFile( "../build/HgcStage2Top/top/DebugFiles/CalculateAveragePositions.txt" ) , cOffset , cHistSize ) );

  seeding_.deinterleave(histogram);
  if( DebugIntermediate ) check( "Deinterleaved2" , histogram , resize( HGCalHistogramCell::ReadDebugFile( "../build/HgcStage2Top/top/DebugFiles/MaximaFanout-Deinterlaced.txt" ) , cOffset , cHistSize ) );

  seeding_.maximaFanout(histogram);
  if( DebugIntermediate ) check( "MaximaFanout" , histogram , resize( HGCalHistogramCell::ReadDebugFile( "../build/HgcStage2Top/top/DebugFiles/MaximaFanout.txt" ) , cOffset , cHistSize ) );  

  // ===================================
  // Clustering
  // ===================================

  HGCalTriggerCellSAPtrCollection triggerCellRamOut;
  HGCalTriggerCellSAPtrCollection clusteredTriggerCellsOut;
  HGCalHistogramCellSAPtrCollection maximaFifoOut;
  clustering_.clusterizer(unpackedTCs, histogram, triggerCellRamOut, maximaFifoOut);
  if( DebugIntermediate ) check( "Clusterizer.TriggerCellsRamOut" , triggerCellRamOut , resize( HGCalTriggerCell::ReadDebugFile( "../build/HgcStage2Top/top/DebugFiles/Clusterizer.TriggerCellsRamOut.txt" ) , triggerCellRamOut.size() ) );

  HGCalClusterSAPtrCollection clusters;
  clustering_.triggerCellToCluster(unpackedTCs, histogram, clusters);
  if( DebugIntermediate ) check( "TriggerCellsToClusters" , clusters , HGCalCluster::ReadDebugFile( "../build/HgcStage2Top/top/DebugFiles/TriggerCellToCluster.txt" ) );

  clustering_.clusterAccumulator(clusters, histogram);
  if( DebugIntermediate ) check( "ClusterAccumulator" , clusters , HGCalCluster::ReadDebugFile( "../build/HgcStage2Top/top/DebugFiles/Accumulator.Clusters.txt" ) );

  clustering_.clusterTree(clusters);
  if( DebugIntermediate ) check( "ClusterTreeSum2" , clusters , HGCalCluster::ReadDebugFile( "../build/HgcStage2Top/top/DebugFiles/ClusterTreeSum2.txt" ) );

  // Cluster properties  
  // clusterProperties( clusters );
}

