#include "multififo_regionizer_ref.h"
#ifdef CMSSW_GIT_HASH
#include "../../egamma/pfeginput_ref.h"
#else
#include "../../egamma/l1-input/ref/pfeginput_ref.h"
#endif

#include <iostream>

#include "multififo_regionizer_elements_ref.icc"

#ifdef CMSSW_GIT_HASH
#include "FWCore/ParameterSet/interface/ParameterSet.h"

l1ct::MultififoRegionizerEmulator::MultififoRegionizerEmulator(const edm::ParameterSet& iConfig)
    : MultififoRegionizerEmulator(
          /*nendcaps=*/2,
          iConfig.getParameter<uint32_t>("nClocks"),
          iConfig.getParameter<uint32_t>("nTrack"),
          iConfig.getParameter<uint32_t>("nCalo"),
          iConfig.getParameter<uint32_t>("nEmCalo"),
          iConfig.getParameter<uint32_t>("nMu"),
          /*streaming=*/false,
          /*outii=*/1,
          iConfig.getParameter<bool>("useAlsoVtxCoords")) {
  debug_ = iConfig.getUntrackedParameter<bool>("debug", false);
}
#endif

l1ct::MultififoRegionizerEmulator::MultififoRegionizerEmulator(unsigned int nendcaps,
                                                               unsigned int nclocks,
                                                               unsigned int ntk,
                                                               unsigned int ncalo,
                                                               unsigned int nem,
                                                               unsigned int nmu,
                                                               bool streaming,
                                                               unsigned int outii,
                                                               bool useAlsoVtxCoords)
    : RegionizerEmulator(useAlsoVtxCoords),
      nendcaps_(nendcaps),
      nclocks_(nclocks),
      ntk_(ntk),
      ncalo_(ncalo),
      nem_(nem),
      nmu_(nmu),
      outii_(outii),
      streaming_(streaming),
      emInterceptMode_(noIntercept),
      init_(false),
      tkRegionizer_(ntk, streaming ? (ntk + outii - 1) / outii : ntk, streaming, outii, useAlsoVtxCoords),
      hadCaloRegionizer_(ncalo, streaming ? (ncalo + outii - 1) / outii : ncalo, streaming, outii),
      emCaloRegionizer_(nem, streaming ? (nem + outii - 1) / outii : nem, streaming, outii),
      muRegionizer_(nmu, streaming ? std::max(1u, (nmu + outii - 1) / outii) : nmu, streaming, outii) {
  // now we initialize the routes: track finder
  for (unsigned int ie = 0; ie < nendcaps && ntk > 0; ++ie) {
    for (unsigned int is = 0; is < NTK_SECTORS; ++is) {  // 9 tf sectors
      for (unsigned int il = 0; il < NTK_LINKS; ++il) {  // max tracks per sector per clock
        unsigned int isp = (is + 1) % NTK_SECTORS, ism = (is + NTK_SECTORS - 1) % NTK_SECTORS;
        tkRoutes_.emplace_back(is + NTK_SECTORS * ie, il, is + NTK_SECTORS * ie, il);
        tkRoutes_.emplace_back(is + NTK_SECTORS * ie, il, isp + NTK_SECTORS * ie, il + 2);
        tkRoutes_.emplace_back(is + NTK_SECTORS * ie, il, ism + NTK_SECTORS * ie, il + 4);
      }
    }
  }
  // hgcal
  assert(NCALO_SECTORS == 3 && NTK_SECTORS == 9);  // otherwise math below is broken, but it's hard to make it generic
  for (unsigned int ie = 0; ie < nendcaps; ++ie) {
    for (unsigned int is = 0; is < NCALO_SECTORS; ++is) {  // NCALO_SECTORS sectors
      for (unsigned int il = 0; il < NCALO_LINKS; ++il) {  // max clusters per sector per clock
        for (unsigned int j = 0; j < 3; ++j) {
          caloRoutes_.emplace_back(is + 3 * ie, il, 3 * is + j + 9 * ie, il);
          if (j) {
            caloRoutes_.emplace_back((is + 1) % 3 + 3 * ie, il, 3 * is + j + 9 * ie, il + 4);
          }
        }
      }
    }
  }
  // mu
  for (unsigned int il = 0; il < NMU_LINKS && nmu > 0; ++il) {  // max clusters per sector per clock
    for (unsigned int j = 0; j < NTK_SECTORS * nendcaps; ++j) {
      muRoutes_.emplace_back(0, il, j, il);
    }
  }
}

l1ct::MultififoRegionizerEmulator::~MultififoRegionizerEmulator() {}

void l1ct::MultififoRegionizerEmulator::setEgInterceptMode(bool afterFifo,
                                                           const l1ct::EGInputSelectorEmuConfig& interceptorConfig) {
  emInterceptMode_ = afterFifo ? interceptPostFifo : interceptPreFifo;
  interceptor_.reset(new EGInputSelectorEmulator(interceptorConfig));
}

void l1ct::MultififoRegionizerEmulator::initSectorsAndRegions(const RegionizerDecodedInputs& in,
                                                              const std::vector<PFInputRegion>& out) {
  assert(!init_);
  init_ = true;
  assert(out.size() == NTK_SECTORS * nendcaps_);
  nregions_ = out.size();
  if (ntk_) {
    assert(in.track.size() == NTK_SECTORS * nendcaps_);
    tkRegionizer_.initSectors(in.track);
    tkRegionizer_.initRegions(out);
    tkRegionizer_.initRouting(tkRoutes_);
  }
  if (ncalo_) {
    assert(in.hadcalo.size() == NCALO_SECTORS * nendcaps_);
    hadCaloRegionizer_.initSectors(in.hadcalo);
    hadCaloRegionizer_.initRegions(out);
    hadCaloRegionizer_.initRouting(caloRoutes_);
  }
  if (nem_) {
    assert(in.emcalo.size() == NCALO_SECTORS * nendcaps_);
    emCaloRegionizer_.initSectors(in.emcalo);
    emCaloRegionizer_.initRegions(out);
    emCaloRegionizer_.initRouting(caloRoutes_);
  }
  if (nmu_) {
    muRegionizer_.initSectors(in.muon);
    muRegionizer_.initRegions(out);
    muRegionizer_.initRouting(muRoutes_);
  }
}

// clock-cycle emulation
bool l1ct::MultififoRegionizerEmulator::step(bool newEvent,
                                             const std::vector<l1ct::TkObjEmu>& links,
                                             std::vector<l1ct::TkObjEmu>& out,
                                             bool mux) {
  return ntk_ ? tkRegionizer_.step(newEvent, links, out, mux) : false;
}

bool l1ct::MultififoRegionizerEmulator::step(bool newEvent,
                                             const std::vector<l1ct::EmCaloObjEmu>& links,
                                             std::vector<l1ct::EmCaloObjEmu>& out,
                                             bool mux) {
  assert(emInterceptMode_ == noIntercept);  // otherwise the em & had calo can't be stepped independently
  return nem_ ? emCaloRegionizer_.step(newEvent, links, out, mux) : false;
}

bool l1ct::MultififoRegionizerEmulator::step(bool newEvent,
                                             const std::vector<l1ct::HadCaloObjEmu>& links,
                                             std::vector<l1ct::HadCaloObjEmu>& out,
                                             bool mux) {
  return ncalo_ ? hadCaloRegionizer_.step(newEvent, links, out, mux) : false;
}

bool l1ct::MultififoRegionizerEmulator::step(bool newEvent,
                                             const std::vector<l1ct::MuObjEmu>& links,
                                             std::vector<l1ct::MuObjEmu>& out,
                                             bool mux) {
  return nmu_ ? muRegionizer_.step(newEvent, links, out, mux) : false;
}

bool l1ct::MultififoRegionizerEmulator::step(bool newEvent,
                                             const std::vector<l1ct::TkObjEmu>& links_tk,
                                             const std::vector<l1ct::HadCaloObjEmu>& links_hadCalo,
                                             const std::vector<l1ct::EmCaloObjEmu>& links_emCalo,
                                             const std::vector<l1ct::MuObjEmu>& links_mu,
                                             std::vector<l1ct::TkObjEmu>& out_tk,
                                             std::vector<l1ct::HadCaloObjEmu>& out_hadCalo,
                                             std::vector<l1ct::EmCaloObjEmu>& out_emCalo,
                                             std::vector<l1ct::MuObjEmu>& out_mu,
                                             bool mux) {
  bool ret = false;
  if (ntk_)
    ret = tkRegionizer_.step(newEvent, links_tk, out_tk, mux);
  if (nmu_)
    ret = muRegionizer_.step(newEvent, links_mu, out_mu, mux);
  switch (emInterceptMode_) {
    case noIntercept:
      if (ncalo_)
        ret = hadCaloRegionizer_.step(newEvent, links_hadCalo, out_hadCalo, mux);
      if (nem_)
        ret = emCaloRegionizer_.step(newEvent, links_emCalo, out_emCalo, mux);
      break;
    case interceptPreFifo:
      // we actually intercept at the links, in the software it's equivalent and it's easier
      assert(nem_ > 0 && ncalo_ > 0 && !links_hadCalo.empty() && links_emCalo.empty());
      assert(interceptor_.get());
      {
        std::vector<l1ct::EmCaloObjEmu> intercepted_links;
        interceptor_->select_or_clear(links_hadCalo, intercepted_links);
        ret = hadCaloRegionizer_.step(newEvent, links_hadCalo, out_hadCalo, mux);
        emCaloRegionizer_.step(newEvent, intercepted_links, out_emCalo, mux);
      }
      break;
    case interceptPostFifo:
      assert(nem_ > 0 && ncalo_ > 0 && !links_hadCalo.empty() && links_emCalo.empty());
      assert(interceptor_.get());
      {
        if (mux) {
          std::vector<l1ct::HadCaloObjEmu> hadNoMux;
          hadCaloRegionizer_.step(newEvent, links_hadCalo, hadNoMux, /*mux=*/false);
          std::vector<l1ct::EmCaloObjEmu> emNoMux(hadNoMux.size());
          interceptor_->select_or_clear(hadNoMux, emNoMux);
          ret = hadCaloRegionizer_.muxonly_step(newEvent, /*flush=*/false, hadNoMux, out_hadCalo);
          emCaloRegionizer_.muxonly_step(newEvent, /*flush=*/true, emNoMux, out_emCalo);
        } else {
          ret = hadCaloRegionizer_.step(newEvent, links_hadCalo, out_hadCalo, /*mux=*/false);
          interceptor_->select_or_clear(out_hadCalo, out_emCalo);
        }
      }
      break;
  }
  return ret;
}

void l1ct::MultififoRegionizerEmulator::fillLinks(unsigned int iclock,
                                                  const l1ct::RegionizerDecodedInputs& in,
                                                  std::vector<l1ct::TkObjEmu>& links) {
  if (ntk_ == 0)
    return;
  links.resize(NTK_SECTORS * NTK_LINKS * nendcaps_);
  for (unsigned int is = 0, idx = 0; is < NTK_SECTORS * nendcaps_; ++is) {  // tf sectors
    const l1ct::DetectorSector<l1ct::TkObjEmu>& sec = in.track[is];
    for (unsigned int il = 0; il < NTK_LINKS; ++il, ++idx) {
      unsigned int ioffs = iclock * NTK_LINKS + il;
      if (ioffs < sec.size() && iclock < nclocks_ - 1) {
        links[idx] = sec[ioffs];
      } else {
        links[idx].clear();
      }
    }
  }
}

template <typename T>
void l1ct::MultififoRegionizerEmulator::fillCaloLinks_(unsigned int iclock,
                                                       const std::vector<DetectorSector<T>>& in,
                                                       std::vector<T>& links) {
  links.resize(NCALO_SECTORS * NCALO_LINKS * nendcaps_);
  for (unsigned int is = 0, idx = 0; is < NCALO_SECTORS * nendcaps_; ++is) {
    for (unsigned int il = 0; il < NCALO_LINKS; ++il, ++idx) {
      unsigned int ioffs = iclock * NCALO_LINKS + il;
      if (ioffs < in[is].size() && iclock < nclocks_ - 1) {
        links[idx] = in[is][ioffs];
      } else {
        links[idx].clear();
      }
    }
  }
}

void l1ct::MultififoRegionizerEmulator::fillLinks(unsigned int iclock,
                                                  const l1ct::RegionizerDecodedInputs& in,
                                                  std::vector<l1ct::HadCaloObjEmu>& links) {
  if (ncalo_ == 0)
    return;
  fillCaloLinks_(iclock, in.hadcalo, links);
}

void l1ct::MultififoRegionizerEmulator::fillLinks(unsigned int iclock,
                                                  const l1ct::RegionizerDecodedInputs& in,
                                                  std::vector<l1ct::EmCaloObjEmu>& links) {
  if (nem_ == 0 || emInterceptMode_ != noIntercept)
    return;
  fillCaloLinks_(iclock, in.emcalo, links);
}

void l1ct::MultififoRegionizerEmulator::fillLinks(unsigned int iclock,
                                                  const l1ct::RegionizerDecodedInputs& in,
                                                  std::vector<l1ct::MuObjEmu>& links) {
  if (nmu_ == 0)
    return;
  links.resize(NMU_LINKS);
  // we have 2 muons on odd clock cycles, and 1 muon on even clock cycles.
  assert(NMU_LINKS == 2);
  for (unsigned int il = 0, idx = 0; il < NMU_LINKS; ++il, ++idx) {
    unsigned int ioffs = (iclock * 3) / 2 + il;
    if (ioffs < in.muon.size() && (il == 0 || iclock % 2 == 1) && iclock < nclocks_ - 1) {
      links[idx] = in.muon[ioffs];
    } else {
      links[idx].clear();
    }
  }
}

void l1ct::MultififoRegionizerEmulator::toFirmware(const std::vector<l1ct::TkObjEmu>& emu,
                                                   TkObj fw[NTK_SECTORS][NTK_LINKS]) {
  if (ntk_ == 0)
    return;
  assert(emu.size() == NTK_SECTORS * NTK_LINKS * nendcaps_);
  for (unsigned int is = 0, idx = 0; is < NTK_SECTORS * nendcaps_; ++is) {  // tf sectors
    for (unsigned int il = 0; il < NTK_LINKS; ++il, ++idx) {
      fw[is][il] = emu[idx];
    }
  }
}
void l1ct::MultififoRegionizerEmulator::toFirmware(const std::vector<l1ct::HadCaloObjEmu>& emu,
                                                   HadCaloObj fw[NCALO_SECTORS][NCALO_LINKS]) {
  if (ncalo_ == 0)
    return;
  assert(emu.size() == NCALO_SECTORS * NCALO_LINKS * nendcaps_);
  for (unsigned int is = 0, idx = 0; is < NCALO_SECTORS * nendcaps_; ++is) {  // tf sectors
    for (unsigned int il = 0; il < NCALO_LINKS; ++il, ++idx) {
      fw[is][il] = emu[idx];
    }
  }
}

void l1ct::MultififoRegionizerEmulator::toFirmware(const std::vector<l1ct::EmCaloObjEmu>& emu,
                                                   EmCaloObj fw[NCALO_SECTORS][NCALO_LINKS]) {
  if (nem_ == 0)
    return;
  assert(emu.size() == NCALO_SECTORS * NCALO_LINKS * nendcaps_);
  for (unsigned int is = 0, idx = 0; is < NCALO_SECTORS * nendcaps_; ++is) {  // tf sectors
    for (unsigned int il = 0; il < NCALO_LINKS; ++il, ++idx) {
      fw[is][il] = emu[idx];
    }
  }
}

void l1ct::MultififoRegionizerEmulator::toFirmware(const std::vector<l1ct::MuObjEmu>& emu, MuObj fw[NMU_LINKS]) {
  if (nmu_ == 0)
    return;
  assert(emu.size() == NMU_LINKS);
  for (unsigned int il = 0, idx = 0; il < NMU_LINKS; ++il, ++idx) {
    fw[il] = emu[idx];
  }
}

void l1ct::MultififoRegionizerEmulator::destream(int iclock,
                                                 const std::vector<l1ct::TkObjEmu>& tk_out,
                                                 const std::vector<l1ct::EmCaloObjEmu>& em_out,
                                                 const std::vector<l1ct::HadCaloObjEmu>& calo_out,
                                                 const std::vector<l1ct::MuObjEmu>& mu_out,
                                                 PFInputRegion& out) {
  if (ntk_)
    tkRegionizer_.destream(iclock, tk_out, out.track);
  if (ncalo_)
    hadCaloRegionizer_.destream(iclock, calo_out, out.hadcalo);
  if (nem_)
    emCaloRegionizer_.destream(iclock, em_out, out.emcalo);
  if (nmu_)
    muRegionizer_.destream(iclock, mu_out, out.muon);
}

void l1ct::MultififoRegionizerEmulator::run(const RegionizerDecodedInputs& in, std::vector<PFInputRegion>& out) {
  if (!init_)
    initSectorsAndRegions(in, out);
  tkRegionizer_.reset();
  emCaloRegionizer_.reset();
  hadCaloRegionizer_.reset();
  muRegionizer_.reset();
  std::vector<l1ct::TkObjEmu> tk_links_in, tk_out;
  std::vector<l1ct::EmCaloObjEmu> em_links_in, em_out;
  std::vector<l1ct::HadCaloObjEmu> calo_links_in, calo_out;
  std::vector<l1ct::MuObjEmu> mu_links_in, mu_out;

  // read and sort the inputs
  for (unsigned int iclock = 0; iclock < nclocks_; ++iclock) {
    fillLinks(iclock, in, tk_links_in);
    fillLinks(iclock, in, em_links_in);
    fillLinks(iclock, in, calo_links_in);
    fillLinks(iclock, in, mu_links_in);

    bool newevt = (iclock == 0), mux = true;
    step(newevt, tk_links_in, calo_links_in, em_links_in, mu_links_in, tk_out, calo_out, em_out, mu_out, mux);
  }

  // set up an empty event
  for (auto& l : tk_links_in)
    l.clear();
  for (auto& l : em_links_in)
    l.clear();
  for (auto& l : calo_links_in)
    l.clear();
  for (auto& l : mu_links_in)
    l.clear();

  // read and put the inputs in the regions
  assert(out.size() == nregions_);
  for (unsigned int iclock = 0; iclock < nclocks_; ++iclock) {
    bool newevt = (iclock == 0), mux = true;
    step(newevt, tk_links_in, calo_links_in, em_links_in, mu_links_in, tk_out, calo_out, em_out, mu_out, mux);

    unsigned int ireg = iclock / outii_;
    if (ireg >= nregions_)
      break;

    if (streaming_) {
      destream(iclock, tk_out, em_out, calo_out, mu_out, out[ireg]);
    } else {
      if (iclock % outii_ == 0) {
        out[ireg].track = tk_out;
        out[ireg].emcalo = em_out;
        out[ireg].hadcalo = calo_out;
        out[ireg].muon = mu_out;
      }
    }
  }

  tkRegionizer_.reset();
  emCaloRegionizer_.reset();
  hadCaloRegionizer_.reset();
  muRegionizer_.reset();
}
