#include "L1Trigger/L1TNtuples/interface/L1AnalysisPhaseIPFJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "Math/Vector4Dfwd.h"

L1Analysis::L1AnalysisPhaseIPFJet::L1AnalysisPhaseIPFJet()
{
}

L1Analysis::L1AnalysisPhaseIPFJet::~L1AnalysisPhaseIPFJet()
{

}



void L1Analysis::L1AnalysisPhaseIPFJet::SetPhaseIPFJet(const edm::Handle< vector<reco::CaloJet> > phaseIPFJets, unsigned maxL1Extra)
{

  for (unsigned int i=0; i<phaseIPFJets->size() && l1extra_.nPhaseIPFJets<maxL1Extra; i++){
    if (phaseIPFJets->at(i).pt()>0){
      // if(phaseIPFJets->at(i).pt()>30 && abs(phaseIPFJets->at(i).eta())<2.4)
	   // l1extra_.phaseIPFJetHt += phaseIPFJets->at(i).pt();
      l1extra_.phaseIPFJetEt .push_back(phaseIPFJets->at(i).pt());
      l1extra_.phaseIPFJetEta.push_back(phaseIPFJets->at(i).eta());
      l1extra_.phaseIPFJetPhi.push_back(phaseIPFJets->at(i).phi());
      l1extra_.nPhaseIPFJets++;
    }
  }
}

void L1Analysis::L1AnalysisPhaseIPFJet::SetPhaseIPFMET  (const edm::Handle< std::vector<l1t::EtSum> >  phaseIPFSums)
{
  std::vector<l1t::EtSum>::const_iterator sumItr = phaseIPFSums->begin();
  const std::vector<l1t::EtSum>::const_iterator sumItrEnd = phaseIPFSums->end();
  for ( ; sumItr != sumItrEnd; ++sumItr ) {
    l1t::EtSum::EtSumType sumType{ sumItr->getType() };

    // std::cout << "Sum : " << sumItr->getType() << " " << sumItr->hwPt() << " " << sumItr->pt() << std::endl;

    // if ( sumType == l1t::EtSum::kTotalHt ) {
    //   l1extra_.phaseIPFJetHt = sumItr->pt();
    // }
    if ( sumType == l1t::EtSum::kMissingEt ) {
      l1extra_.phaseIPFJetMET = sumItr->pt();
    }
    else if ( sumType == l1t::EtSum::kMissingEtHF ) {
      l1extra_.phaseIPFJetMETHF = sumItr->pt();     
    }
  } 
}

void L1Analysis::L1AnalysisPhaseIPFJet::SetPhaseIPFJetSums  (const edm::Handle< std::vector<l1t::EtSum> >  phaseIPFSums)
{
  std::vector<l1t::EtSum>::const_iterator sumItr = phaseIPFSums->begin();
  const std::vector<l1t::EtSum>::const_iterator sumItrEnd = phaseIPFSums->end();
  for ( ; sumItr != sumItrEnd; ++sumItr ) {
    l1t::EtSum::EtSumType sumType{ sumItr->getType() };

    // std::cout << "Sum : " << sumItr->getType() << " " << sumItr->hwPt() << " " << sumItr->pt() << std::endl;

    if ( sumType == l1t::EtSum::kTotalHt ) {
      l1extra_.phaseIPFJetHt = sumItr->pt();
    }
    // else if ( sumType == l1t::EtSum::kMissingEt ) {
    //   l1extra_.phaseIPFJetMET = sumItr->pt();
    // }
    // else if ( sumType == l1t::EtSum::kMissingEtHF ) {
    //   l1extra_.phaseIPFJetMETHF = sumItr->pt();     
    // }

    else if ( sumType == l1t::EtSum::kMissingHt){
      l1extra_.phaseIPFJetMHT = sumItr->pt();
    }
  } 
}

void L1Analysis::L1AnalysisPhaseIPFJet::SetPFJet(const edm::Handle< vector<l1t::PFJet> > ak4PFJets, unsigned maxL1Extra)
{
  ROOT::Math::PxPyPzEVector myLorentzVector(0, 0, 0, 0);

  for (unsigned int i=0; i<ak4PFJets->size() && l1extra_.nPhaseIPFJets<maxL1Extra; i++){
    if (ak4PFJets->at(i).pt()>0){
      if(ak4PFJets->at(i).pt()>30 && abs(ak4PFJets->at(i).eta())<2.4){
	      l1extra_.ak4PFJetHt += ak4PFJets->at(i).pt();
        
        myLorentzVector += ak4PFJets->at(i).p4();
        
      }
      l1extra_.ak4PFJetEt .push_back(ak4PFJets->at(i).pt());
      l1extra_.ak4PFJetEta.push_back(ak4PFJets->at(i).eta());
      l1extra_.ak4PFJetPhi.push_back(ak4PFJets->at(i).phi());
      l1extra_.nAK4PFJets++;
      
      
    }
  }
  l1extra_.ak4PFJetMHT = myLorentzVector.pt();
}

void L1Analysis::L1AnalysisPhaseIPFJet::SetPFMET(const edm::Handle<reco::PFMETCollection> pfMET)
{
  l1extra_.pfMet = pfMET->at(0).pt();
}


void L1Analysis::L1AnalysisPhaseIPFJet::SetGenJet(const edm::Handle<reco::GenJetCollection> genJets, unsigned maxL1Extra)
{
  ROOT::Math::PxPyPzEVector myLorentzVector(0, 0, 0, 0);

  reco::GenJetCollection::const_iterator genJetItr = genJets->begin();
  reco::GenJetCollection::const_iterator genJetEnd = genJets->end();
  for( ; genJetItr != genJetEnd ; ++genJetItr) {
    
    if(genJetItr->pt()>30 && abs(genJetItr->eta())<2.4){
      l1extra_.genJetHt += genJetItr->pt();  
    
      myLorentzVector += genJetItr->p4();
    }
    l1extra_.genJetPt.push_back( genJetItr->pt() );
    l1extra_.genJetEta.push_back( genJetItr->eta() );
    l1extra_.genJetPhi.push_back( genJetItr->phi() );
    l1extra_.genJetM.push_back( genJetItr->mass() );
    l1extra_.nGenJet++; 
    
  }
  l1extra_.genJetMHT = myLorentzVector.pt(); 
 
}

void L1Analysis::L1AnalysisPhaseIPFJet::SetGenMET(const edm::Handle<reco::GenMETCollection> genMET)
{

  l1extra_.genMet = genMET->at(0).pt();
}


