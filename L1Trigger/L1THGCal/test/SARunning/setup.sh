#!/bin/bash

for f in L1Trigger/L1THGCal/interface/backend/HGCalHistoClusteringImpl_SA.h L1Trigger/L1THGCal/src/backend/HGCalHistoClusteringImpl_SA.cc DataFormats/L1THGCal/interface/HGCalCluster_SA.h DataFormats/L1THGCal/interface/HGCalSeed_SA.h DataFormats/L1THGCal/interface/HGCalMulticluster_SA.h DataFormats/L1THGCal/src/HGCalMulticluster_SA.cc L1Trigger/L1THGCal/test/SARunning/drCoeffA.icc L1Trigger/L1THGCal/test/SARunning/input_seeds.icc L1Trigger/L1THGCal/test/SARunning/input_clusters.icc L1Trigger/L1THGCal/test/SARunning/main.cc;
    do wget -x -nH --cut-dirs 3 https://raw.githubusercontent.com/EmyrClement/cmssw/SADev-Clustering/$f
done

