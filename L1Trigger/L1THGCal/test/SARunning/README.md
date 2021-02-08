Run setup.sh (copy & paste contents outside of CMSSW area, before cmsenv etc.), then compile e.g.
```
g++  -std=c++11 -I./  main.cc L1Trigger/L1THGCal/src/backend/HGCalHistoClusteringImpl_SA.cc DataFormats/L1THGCal/src/HGCalMulticluster_SA.cc
```
Run:
```
./a.out
```
