#include "L1Trigger/Phase2L1ParticleFlow/interface/jetmet/L1SeedConePFJetEmulator.h"

L1SCJetEmu::L1SCJetEmu(bool debug, float coneSize, unsigned nJets)
    : debug_(debug),
      coneSize_(coneSize),
      nJets_(nJets),
      rCone2_(coneSize * coneSize / l1ct::Scales::ETAPHI_LSB / l1ct::Scales::ETAPHI_LSB) {
  init_invert_table<pt_t, inv_pt_t, N_table_inv_pt>(inv_pt_table_);
}

L1SCJetEmu::detaphi_t L1SCJetEmu::deltaPhi(L1SCJetEmu::Particle a, L1SCJetEmu::Particle b) {
  detaphi_t dphi = detaphi_t(a.hwPhi) - detaphi_t(b.hwPhi);
  // phi wrap
  detaphi_t dphi0 =
      dphi > detaphi_t(l1ct::Scales::INTPHI_PI) ? detaphi_t(l1ct::Scales::INTPHI_TWOPI - dphi) : detaphi_t(dphi);
  detaphi_t dphi1 =
      dphi < detaphi_t(-l1ct::Scales::INTPHI_PI) ? detaphi_t(l1ct::Scales::INTPHI_TWOPI + dphi) : detaphi_t(dphi);
  detaphi_t dphiw = dphi > detaphi_t(0) ? dphi0 : dphi1;
  return dphiw;
}

bool L1SCJetEmu::inCone(L1SCJetEmu::Particle seed, L1SCJetEmu::Particle part) const {
  // scale the particle eta, phi to hardware units
  detaphi_t deta = detaphi_t(seed.hwEta) - detaphi_t(part.hwEta);
  detaphi_t dphi = deltaPhi(seed, part);
  bool ret = deta * deta + dphi * dphi < rCone2_;
  //bool ret = r2 < cone2;
  if (debug_) {
    detaphi2_t r2 = detaphi2_t(deta) * detaphi2_t(deta) + detaphi2_t(dphi) * detaphi2_t(dphi);
    dbgCout() << "  part eta, seed eta: " << part.hwEta << ", " << seed.hwEta << std::endl;
    dbgCout() << "  part phi, seed phi: " << part.hwPhi << ", " << seed.hwPhi << std::endl;
    dbgCout() << "  pt, deta, dphi, r2, cone2, lt: " << part.hwPt << ", " << deta << ", " << dphi << ", "
              << deta * deta + dphi * dphi << ", " << rCone2_ << ", " << ret << std::endl;
  }
  return ret;
}

L1SCJetEmu::Jet L1SCJetEmu::makeJet_HW(const std::vector<Particle>& parts, const Particle seed) const {
  // Seed Cone Jet algorithm with ap_fixed types and hardware emulation
  // Particle seed = reduce(parts, op_max);

  // Event with saturation, order of terms doesn't matter since they're all positive
  auto sumpt = [](pt_t(a), const Particle& b) { return a + b.hwPt; };

  // Sum the pt
  pt_t pt = std::accumulate(parts.begin(), parts.end(), pt_t(0), sumpt);
  inv_pt_t inv_pt = invert_with_shift<pt_t, inv_pt_t, N_table_inv_pt>(pt, inv_pt_table_, false);

  // pt weighted d eta
  std::vector<pt_etaphi_t> pt_deta;
  pt_deta.resize(parts.size());
  std::transform(parts.begin(), parts.end(), pt_deta.begin(), [&seed](const Particle& part) {
    // In the firmware we calculate the per-particle pt-weighted deta
    return pt_etaphi_t(part.hwPt * detaphi_t(part.hwEta - seed.hwEta));
  });
  // Accumulate the pt-weighted etas. Init to 0, include seed in accumulation
  pt_etaphi_t sum_pt_eta = std::accumulate(pt_deta.begin(), pt_deta.end(), pt_etaphi_t(0));
  etaphi_t eta = seed.hwEta + etaphi_t(sum_pt_eta * inv_pt);

  // pt weighted d phi
  std::vector<pt_etaphi_t> pt_dphi;
  pt_dphi.resize(parts.size());
  std::transform(parts.begin(), parts.end(), pt_dphi.begin(), [&seed](const Particle& part) {
    // In the firmware we calculate the per-particle pt-weighted dphi
    return pt_etaphi_t(part.hwPt * deltaPhi(part, seed));
  });
  // Accumulate the pt-weighted phis. Init to 0, include seed in accumulation
  pt_etaphi_t sum_pt_phi = std::accumulate(pt_dphi.begin(), pt_dphi.end(), pt_etaphi_t(0));
  etaphi_t phi = seed.hwPhi + etaphi_t(sum_pt_phi * inv_pt);

  Jet jet;
  jet.hwPt = pt;
  jet.hwEta = eta;
  jet.hwPhi = phi;
  jet.constituents = parts;

  if (debug_) {
    std::for_each(pt_dphi.begin(), pt_dphi.end(), [](pt_etaphi_t& x) { dbgCout() << "pt_dphi: " << x << std::endl; });
    std::for_each(pt_deta.begin(), pt_deta.end(), [](pt_etaphi_t& x) { dbgCout() << "pt_deta: " << x << std::endl; });
    dbgCout() << " sum_pt_eta: " << sum_pt_eta << ", 1/pt: " << inv_pt
              << ", sum_pt_eta * 1/pt: " << etaphi_t(sum_pt_eta * inv_pt) << std::endl;
    dbgCout() << " sum_pt_phi: " << sum_pt_phi << ", 1/pt: " << inv_pt
              << ", sum_pt_phi * 1/pt: " << etaphi_t(sum_pt_phi * inv_pt) << std::endl;
    dbgCout() << " uncorr eta: " << seed.hwEta << ", phi: " << seed.hwPhi << std::endl;
    dbgCout() << "   corr eta: " << eta << ", phi: " << phi << std::endl;
    dbgCout() << "         pt: " << pt << std::endl;
  }

  return jet;
}

std::vector<L1SCJetEmu::Jet> L1SCJetEmu::emulateEvent(std::vector<Particle>& parts, std::vector<Particle>& seeds, bool useExternalSeeds, bool allowDoubleCounting) const {
  // The fixed point algorithm emulation
  std::vector<Particle> work;
  work.resize(parts.size());
  std::transform(parts.begin(), parts.end(), work.begin(), [](const Particle& part) { return part; });

  std::vector<Jet> jets;
  jets.reserve(nJets_);
  while ( ( !work.empty() && jets.size() < nJets_ )  ) {
    if ( useExternalSeeds && seeds.size() == 0 ) break; // Can this be combined with previous line?
    // Take the highest pt candidate as a seed
    // Use the firmware reduce function to find the same seed as the firmware
    // in case there are multiple seeds with the same pT
    // ... or use external seed if configured to do so
    Particle seed = (useExternalSeeds) ? seeds.at(0) : reduce(work, op_max);

    // Get the particles within a coneSize_ of the seed
    std::vector<Particle> particlesInCone;
    std::copy_if(work.begin(), work.end(), std::back_inserter(particlesInCone), [&](const Particle& part) {
      return inCone(seed, part);
    });
    if (debug_) {
      dbgCout() << "Seed: " << seed.hwPt << ", " << seed.hwEta << ", " << seed.hwPhi << std::endl;
      std::cout << "N particles : " << particlesInCone.size() << std::endl;
      std::for_each(particlesInCone.begin(), particlesInCone.end(), [&](Particle& part) {
        dbgCout() << "  Part: " << part.hwPt << ", " << part.hwEta << ", " << part.hwPhi << std::endl;
        inCone(seed, part);
      });
    }
    if ( particlesInCone.size() > 0 ) { // Possible hack - some seeds don't have any clustered particles.  Need to understand if real effect (could be) or a bug

     //--------------------------------------- Constituents print  (debug for FW mismatch)  --start ----------------------------
     if(coneSize_ < 0.5){
     std::for_each(particlesInCone.begin(), particlesInCone.end(), [&](Particle& part) {
        std::cout << "  Part (int) --SC4: " << part.hwPt << ", " << part.hwEta << ", " << part.hwPhi << std::endl;
     });
    int count_1=0;
     std::for_each(particlesInCone.begin(), particlesInCone.end(), [&](Particle& part) {
        std::cout << "  Part (float) --SC4: " << part.hwPt << ", " << part.floatEta()  << ", " << part.floatPhi()  << std::endl;
        count_1++;
     });
     std::cout<<"Total particles: "<<count_1<<std::endl;
     }

     if(coneSize_ >0.6){
     std::for_each(particlesInCone.begin(), particlesInCone.end(), [&](Particle& part) {
        std::cout << "  Part (int) --SC8: " << part.hwPt << ", " << part.hwEta << ", " << part.hwPhi << std::endl;
     });
    int count_1=0;
     std::for_each(particlesInCone.begin(), particlesInCone.end(), [&](Particle& part) {
        std::cout << "  Part (float) --SC8: " << part.hwPt << ", " << part.floatEta()  << ", " << part.floatPhi()  << std::endl;
        count_1++;
     });
     std::cout<<"Total particles: "<<count_1<<std::endl;
     }
//--------------------------------------- Constituents print  (debug for FW mismatch)  --end ----------------------------
    // if (debug_){ dbgCout() << "Internal seed used!" << std::endl;}   //debug
     jets.push_back(makeJet_HW(particlesInCone, seed));
     //remove the clustered particles
     if(! allowDoubleCounting){       //if double couting not allowed
    //  if (debug_){ dbgCout() << "Removing candidates!" << std::endl;}   //debug
      work.erase(std::remove_if(work.begin(), work.end(), [&](const Particle& part) { return inCone(seed, part); }),  //erase particles from further jet clustering
               work.end());
     }
   }

    if ( useExternalSeeds ) {
      // if (debug_) {dbgCout() << "External seed used!" << std::endl;}
      seeds.erase(seeds.begin());
      // if (debug_){ dbgCout() << "N seeds remaining and: " << seeds.size() << std::endl;}

    }


  }
  
  //Jet print (without JEC)
  if(coneSize_ < 0.5){
    std::vector<Jet> SC4Jets = jets;
     std::sort(SC4Jets.begin(), SC4Jets.end(), [](Jet seed1, Jet seed2) {
          return seed1.hwPt > seed2.hwPt;
        });
     std::cout << "--------------SC4 (no JEC)----------"<<std::endl;
  std::for_each(SC4Jets.begin(), SC4Jets.end(), [&](const Jet& part) {
      if(part.hwPt >0){
        std::cout << "  Jets : " << part.hwPt << ", " << part.hwEta << ", " << part.hwPhi << std::endl;
      }
    });
  }
  if(coneSize_ > 0.6){
    std::vector<Jet> SC8Jets = jets;
     std::cout << "--------------SC8 (no JEC)----------"<<std::endl;
     std::sort(SC8Jets.begin(), SC8Jets.end(), [](Jet seed1, Jet seed2) {
          return seed1.hwPt > seed2.hwPt;
        });
  std::for_each(SC8Jets.begin(), SC8Jets.end(), [&](const Jet& part) {
      if(part.hwPt >0){
        std::cout << "  Jets : " << part.hwPt << ", " << part.hwEta << ", " << part.hwPhi << std::endl;
      }
    });
  }
  
  return jets;
}
