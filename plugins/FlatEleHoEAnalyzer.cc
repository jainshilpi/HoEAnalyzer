// -*- C++ -*-
// EleHoEAnalyzer.cc, but flattened out ie no vectors etc
// also, all hcal rechits are saved
// index pointing to ele is also saved if the min_dieta and min_diphi of hit and ele points to the same ele
// dieta and diphi kept signed

// system include files
#include <memory>
#include <cassert>

// user include files
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHadTower.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloEventSetup/plugins/CaloTopologyBuilder.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "TFile.h"
#include "TTree.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/GsfElectronTools.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/GsfElectronAlgo.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include <cmath>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHcalIsolation.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronHcalHelper.h"
#include "Geometry/CaloTopology/interface/CaloTowerConstituentsMap.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EGHcalRecHitSelector.h"

#include "DataFormats/PatCandidates/interface/Electron.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.
class FlatEleHoEAnalyzer : public edm::one::EDAnalyzer<>  {
public:
  explicit FlatEleHoEAnalyzer(const edm::ParameterSet&);
  ~FlatEleHoEAnalyzer();
  
  static edm::ParameterSetDescription makePSetDescription();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  TFile *file;
  TTree *tree;

  int run, lumi_block, event, bunch_crossing, orbit_number, store_number;

  int n_ele;
  std::vector<int>    ele_eb;
  std::vector<int>    ele_ee;
  std::vector<int>    ele_gap_eb_ee;
  std::vector<int>    ele_gap_eb_eta;
  std::vector<int>    ele_gap_eb_phi;
  std::vector<int>    ele_gap_ee_dee;
  std::vector<int>    ele_gap_ee_ring;

  std::vector<int>    ele_golden;
  std::vector<int>    ele_unknown;
  std::vector<int>    ele_bigbrem;
  std::vector<int>    ele_badtrack;
  std::vector<int>    ele_showering;
  std::vector<int>    ele_gap;

  std::vector<int>    ele_genele;
  std::vector<float>  ele_dR_reco_genele;
  std::vector<float>  ele_rpt_reco_genele;

  std::vector<int>    ele_genpho;
  std::vector<float>  ele_dR_reco_genpho;
  std::vector<float>  ele_rpt_reco_genpho;

  std::vector<float>  ele_track_fbrem;
  std::vector<float>  ele_sc_fbrem;
  std::vector<int>    ele_nbrem;
  std::vector<float>  ele_sc_energy;
  std::vector<float>  ele_sc_raw_energy;
  std::vector<float>  ele_ecal_energy;
  std::vector<float>  ele_seed_energy;
  std::vector<float>  ele_seed_corr_energy;
  std::vector<float>  ele_cmssw_hoe;
  std::vector<float>  ele_cmssw_hoe_tower;
  std::vector<float>  ele_cmssw_hoe_5x5;
  std::vector<float>  ele_sc_eta;
  std::vector<float>  ele_sc_phi;
  std::vector<float>  ele_pt;
  std::vector<float>  ele_eta;
  std::vector<float>  ele_phi;
  std::vector<float>  ele_sieie_5x5;
  std::vector<float>  ele_r9_5x5;

  std::vector<float>  ele_pfiso_pho;
  std::vector<float>  ele_pfiso_neu;
  std::vector<float>  ele_pfiso_cha;
  std::vector<float>  ele_pfiso_pu;
  std::vector<float>  ele_pfiso_hcal;
  std::vector<float>  ele_pfiso_ecal;

  std::vector<float>  ele_detiso03_ecalhit;
  std::vector<float>  ele_detiso03_hcaltower1;
  std::vector<float>  ele_detiso03_hcaltower2;
  std::vector<float>  ele_detiso03_trk;
  std::vector<float>  ele_detiso03_trk_heep;

  std::vector<int>    ele_seed_detid;
  std::vector<int>    ele_seed_subdetid;
  std::vector<int>    ele_seed_ieta;
  std::vector<int>    ele_seed_iphi;
  std::vector<float>  ele_seed_eta;
  std::vector<float>  ele_seed_phi;
  std::vector<int>    ele_seed_raw_id;
  std::vector<int>    ele_seed_hcal_ieta;
  std::vector<int>    ele_seed_hcal_iphi;

  std::vector<std::vector<int> > ele_IDbits;
  std::vector<int>  ele_IDVeto;
  std::vector<int>  ele_IDLoose;
  std::vector<int>  ele_IDMedium;
  std::vector<int>  ele_IDTight;
  std::vector<int>  ele_IDMVAiso90;
  std::vector<float>  ele_dEtaIn;
  std::vector<float>  ele_dPhiIn;
  std::vector<float>  ele_1oEm1op;
  std::vector<float>  ele_eSCoP;
  std::vector<float>  ele_eSCoPout;
  std::vector<float>  ele_psEorawE;
  std::vector<float>  ele_gsfTrackChi2;
  //std::vector<float>  ele_convVtxFitProb;
  std::vector<int>  ele_nHit;
  std::vector<int>  ele_missingHit;
  //std::vector<int>  ele_convVeto;
  std::vector<int>  ele_isEcalDriven;
  std::vector<int>  ele_IDHEEP;
  std::vector<int>  ele_IDHEEPBit;

  int n_hcalhit;
  std::vector<int>    hcalhit_ieta;
  std::vector<int>    hcalhit_iphi;
  std::vector<float>  hcalhit_energy;
  std::vector<int>    hcalhit_seed_dieta;
  std::vector<int>    hcalhit_seed_diphi;
  std::vector<int>    hcalhit_raw_id;
  std::vector<int>    hcalhit_depth;
  std::vector<int>    hcalhit_ele_index;
  std::vector<float>  hcalhit_eta;
  std::vector<float>  hcalhit_phi;

  int imin, min_dieta, min_diphi;
  float min_diR2;

  float pu_true;
  int pu_obs;
  float rho;
  float gen_weight;

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  static int calDIEta(int iEta1, int iEta2);
  static int calDIPhi(int iPhi1, int iPhi2);
  void reallocate_setaddress(int n_ele_ = 0, int n_hcalhit_ = 0);
  float getMinEnergyHCAL(HcalDetId id) const;

  int maxDIEta_ = 5;
  int maxDIPhi_ = 5;

  // ----------member data ---------------------------
  //edm::EDGetTokenT<edm::View<reco::GsfElectron> > eleToken_;

  edm::EDGetTokenT<edm::View<pat::Electron> > eleToken_;

  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puCollection_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<GenEventInfoProduct> genEventToken_;
  edm::EDGetTokenT<HBHERecHitCollection> hbhe_rechits_;
  edm::EDGetTokenT<EcalRecHitCollection> ebReducedRecHitCollection_;
  edm::EDGetTokenT<EcalRecHitCollection> eeReducedRecHitCollection_;
  edm::EDGetTokenT<EcalRecHitCollection> esReducedRecHitCollection_;
  edm::ESHandle<CaloGeometry> theCaloGeometry;  
  edm::ESHandle<CaloTowerConstituentsMap> towerMap_;
  edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticlesCollection_;

  /*
  edm::EDGetTokenT<edm::ValueMap<bool> > vidPassToken_;
  edm::EDGetTokenT<edm::ValueMap<unsigned int> > vidBitmapToken_;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> >  vidResultToken_;
  
  
  edm::EDGetTokenT<edm::ValueMap<float> > trkIsoMapToken_; 
  */
  
  std::string output;
  bool Run2_2018 ; // Now two options are supported, Run2_2018 and Run3
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
FlatEleHoEAnalyzer::FlatEleHoEAnalyzer(const edm::ParameterSet& iConfig) :
  //eleToken_(consumes<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("electrons"))),

  eleToken_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"))),

  puCollection_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupCollection"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoSrc"))),
  genEventToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventSrc"))),
  hbhe_rechits_(consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("hbheInput"))),
  ebReducedRecHitCollection_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ebReducedRecHitCollection"))),
  eeReducedRecHitCollection_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("eeReducedRecHitCollection"))),
  esReducedRecHitCollection_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("esReducedRecHitCollection"))),
  genParticlesCollection_(consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleSrc"))),
  output(iConfig.getParameter<std::string>("output_file")),
  Run2_2018(iConfig.getParameter<bool>("Run2_2018_"))
 
{
  //now do what ever initialization is needed

  /*
  vidPassToken_=consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("vid"));
  vidBitmapToken_=consumes<edm::ValueMap<unsigned int> >(iConfig.getParameter<edm::InputTag>("vidBitmap"));
  vidResultToken_=consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("vid"));
  trkIsoMapToken_=consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("trkIsoMap"));
  */


  file = new TFile(output.c_str(), "recreate");
  tree = new TTree("tree", "");
  
  tree->Branch("run", &run, "run/I");
  tree->Branch("lumi_block", &lumi_block, "lumi_block/I");
  tree->Branch("event", &event, "event/I");
  tree->Branch("bunch_crossing", &bunch_crossing, "bunch_crossing/I");
  tree->Branch("orbit_number", &orbit_number, "orbit_number/I");
  tree->Branch("store_number", &store_number, "store_number/I");
  
  tree->Branch("n_ele", &n_ele, "n_ele/I");
  tree->Branch("n_hcalhit", &n_hcalhit, "n_hcalhit/I");
  tree->Branch("pu_true", &pu_true, "pu_true/F");
  tree->Branch("pu_obs", &pu_obs, "pu_obs/I");
  tree->Branch("rho", &rho, "rho/F");
  
  tree->Branch("gen_weight", &gen_weight, "gen_weight/F");
  //}
  
  tree->Branch("ele_eb", &ele_eb);
  tree->Branch("ele_ee", &ele_ee);
  tree->Branch("ele_gap_eb_ee", &ele_gap_eb_ee);
  tree->Branch("ele_gap_eb_eta", &ele_gap_eb_eta);
  tree->Branch("ele_gap_eb_phi", &ele_gap_eb_phi);
  tree->Branch("ele_gap_ee_dee", &ele_gap_ee_dee);
  tree->Branch("ele_gap_ee_ring", &ele_gap_ee_ring);

  tree->Branch("ele_golden", &ele_golden);
  tree->Branch("ele_unknown", &ele_unknown);
  tree->Branch("ele_bigbrem", &ele_bigbrem);
  tree->Branch("ele_gap", &ele_gap);
  tree->Branch("ele_badtrack", &ele_badtrack);
  tree->Branch("ele_showering", &ele_showering);
  
  tree->Branch("ele_track_fbrem", &ele_track_fbrem);
  tree->Branch("ele_sc_fbrem", &ele_sc_fbrem);
  tree->Branch("ele_nbrem", &ele_nbrem);

  tree->Branch("ele_genele", &ele_genele);
  tree->Branch("ele_dR_reco_genele", &ele_dR_reco_genele);
  tree->Branch("ele_rpt_reco_genele", &ele_rpt_reco_genele);
  
  tree->Branch("ele_genpho", &ele_genpho);
  tree->Branch("ele_dR_reco_genpho", &ele_dR_reco_genpho);
  tree->Branch("ele_rpt_reco_genpho", &ele_rpt_reco_genpho);
  
  tree->Branch("ele_sc_energy", &ele_sc_energy);
  tree->Branch("ele_sc_raw_energy", &ele_sc_raw_energy);
  tree->Branch("ele_ecal_energy", &ele_ecal_energy);
  tree->Branch("ele_seed_energy", &ele_seed_energy);
  tree->Branch("ele_seed_corr_energy", &ele_seed_corr_energy);
  tree->Branch("ele_cmssw_hoe", &ele_cmssw_hoe);
  tree->Branch("ele_cmssw_hoe_tower", &ele_cmssw_hoe_tower);
  tree->Branch("ele_cmssw_hoe_5x5", &ele_cmssw_hoe_5x5);
  tree->Branch("ele_sc_eta", &ele_sc_eta);
  tree->Branch("ele_sc_phi", &ele_sc_phi);
  tree->Branch("ele_pt", &ele_pt);
  tree->Branch("ele_eta", &ele_pt);
  tree->Branch("ele_phi", &ele_phi);
  tree->Branch("ele_sieie_5x5", &ele_sieie_5x5);
  tree->Branch("ele_r9_5x5", &ele_r9_5x5);

  tree->Branch("ele_pfiso_pho", &ele_pfiso_pho);
  tree->Branch("ele_pfiso_neu", &ele_pfiso_neu);
  tree->Branch("ele_pfiso_cha", &ele_pfiso_cha);
  tree->Branch("ele_pfiso_pu", &ele_pfiso_pu);
  tree->Branch("ele_pfiso_hcal", &ele_pfiso_hcal);
  tree->Branch("ele_pfiso_ecal", &ele_pfiso_ecal);

  tree->Branch("ele_detiso03_ecalhit", &ele_detiso03_ecalhit);
  tree->Branch("ele_detiso03_hcaltower1", &ele_detiso03_hcaltower1);
  tree->Branch("ele_detiso03_hcaltower2", &ele_detiso03_hcaltower2);
  tree->Branch("ele_detiso03_trk", &ele_detiso03_trk);
  tree->Branch("ele_detiso03_trk_heep", &ele_detiso03_trk_heep);

  tree->Branch("ele_seed_detid", &ele_seed_detid);
  tree->Branch("ele_seed_subdetid", &ele_seed_subdetid);
  tree->Branch("ele_seed_ieta", &ele_seed_ieta);
  tree->Branch("ele_seed_iphi", &ele_seed_iphi);
  tree->Branch("ele_seed_eta", &ele_seed_eta);
  tree->Branch("ele_seed_phi", &ele_seed_phi);
  tree->Branch("ele_seed_raw_id", &ele_seed_raw_id);
  tree->Branch("ele_seed_hcal_ieta", &ele_seed_hcal_ieta);
  tree->Branch("ele_seed_hcal_iphi", &ele_seed_hcal_iphi);
  
  tree->Branch("hcalhit_ieta", &hcalhit_ieta);
  tree->Branch("hcalhit_iphi", &hcalhit_iphi);
  tree->Branch("hcalhit_energy", &hcalhit_energy);
  tree->Branch("hcalhit_seed_dieta", &hcalhit_seed_dieta);
  tree->Branch("hcalhit_seed_diphi", &hcalhit_seed_diphi);
  tree->Branch("hcalhit_raw_id", &hcalhit_raw_id);
  tree->Branch("hcalhit_depth", &hcalhit_depth);
  tree->Branch("hcalhit_ele_index", &hcalhit_ele_index);
  tree->Branch("hcalhit_eta", &hcalhit_eta);
  tree->Branch("hcalhit_phi", &hcalhit_phi);


  tree->Branch("ele_IDbits", &ele_IDbits);
  tree->Branch("ele_IDVeto", &ele_IDVeto);
  tree->Branch("ele_IDLoose", &ele_IDLoose);
  tree->Branch("ele_IDMedium", &ele_IDMedium);
  tree->Branch("ele_IDTight", &ele_IDTight);
  tree->Branch("ele_IDMVAiso90", &ele_IDMVAiso90);
  tree->Branch("ele_IDHEEP", &ele_IDHEEP);
  tree->Branch("ele_dEtaIn", &ele_dEtaIn);
  tree->Branch("ele_dPhiIn", &ele_dPhiIn);
  tree->Branch("ele_psEorawE", &ele_psEorawE);
  tree->Branch("ele_1oEm1op", &ele_1oEm1op);
  tree->Branch("ele_eSCoP", &ele_eSCoP);
  tree->Branch("ele_eSCoPout", &ele_eSCoPout);
  
  //tree->Branch("ele_convVtxFitProb",ele_convVtxFitProb.data(),"ele_convVtxFitProb[n_ele]/I");
  
  tree->Branch("ele_gsfTrackChi2", &ele_gsfTrackChi2);
  tree->Branch("ele_nHit", &ele_nHit);
  tree->Branch("ele_missingHit", &ele_missingHit);
  tree->Branch("ele_isEcalDriven", &ele_isEcalDriven);


}


FlatEleHoEAnalyzer::~FlatEleHoEAnalyzer()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called for each event  ------------
void
FlatEleHoEAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  n_ele = 0;
  ele_eb.clear();
  ele_ee.clear();
  ele_gap_eb_ee.clear();
  ele_gap_eb_eta.clear();
  ele_gap_eb_phi.clear();
  ele_gap_ee_dee.clear();
  ele_gap_ee_ring.clear();

  ele_golden.clear();
  ele_unknown.clear();
  ele_badtrack.clear();
  ele_bigbrem.clear();
  ele_showering.clear();
  ele_gap.clear();

  ele_genele.clear();
  ele_rpt_reco_genele.clear();
  ele_dR_reco_genele.clear();

  ele_genpho.clear();
  ele_rpt_reco_genpho.clear();
  ele_dR_reco_genpho.clear();

  ele_track_fbrem.clear();
  ele_sc_fbrem.clear();
  ele_nbrem.clear();

  ele_sc_energy.clear();
  ele_sc_raw_energy.clear();
  ele_ecal_energy.clear();
  ele_seed_energy.clear();
  ele_seed_corr_energy.clear();
  ele_cmssw_hoe.clear();
  ele_cmssw_hoe_tower.clear();
  ele_cmssw_hoe_5x5.clear();
  ele_sc_eta.clear();
  ele_sc_phi.clear();
  ele_pt.clear();
  ele_eta.clear();
  ele_phi.clear();
  ele_sieie_5x5.clear();
  ele_r9_5x5.clear();

  ele_pfiso_pho.clear();
  ele_pfiso_neu.clear();
  ele_pfiso_cha.clear();
  ele_pfiso_pu.clear();
  ele_pfiso_hcal.clear();
  ele_pfiso_ecal.clear();

  ele_detiso03_ecalhit.clear();
  ele_detiso03_hcaltower1.clear();
  ele_detiso03_hcaltower2.clear();
  ele_detiso03_trk.clear();
  ele_detiso03_trk_heep.clear();

  ele_seed_detid.clear();
  ele_seed_subdetid.clear();
  ele_seed_ieta.clear();
  ele_seed_iphi.clear();
  ele_seed_eta.clear();
  ele_seed_phi.clear();
  ele_seed_raw_id.clear();
  ele_seed_hcal_ieta.clear();
  ele_seed_hcal_iphi.clear();

  ////Electron ID variables
  ele_dEtaIn.clear();
  ele_dPhiIn.clear();
  ele_psEorawE.clear();
  ele_1oEm1op.clear();
  ele_eSCoP.clear();
  ele_eSCoPout.clear();
  //ele_convVtxFitProb.clear();
  ele_gsfTrackChi2.clear();
  ele_nHit.clear();
  ele_missingHit.clear();
  //ele_convVeto.clear();
  ele_IDVeto.clear();
  ele_IDLoose.clear();
  ele_IDMedium.clear();
  ele_IDTight.clear();
  ele_IDMVAiso90.clear();
  ele_IDbits.clear();
  ele_isEcalDriven.clear();
  //ele_TrackIso.clear();
  ele_IDHEEP.clear();
  //ele_IDHEEPBit.clear();
  

  n_hcalhit = 0;
  hcalhit_ieta.clear();
  hcalhit_iphi.clear();
  hcalhit_energy.clear();
  hcalhit_seed_dieta.clear();
  hcalhit_seed_diphi.clear();
  hcalhit_raw_id.clear();
  hcalhit_depth.clear();
  hcalhit_ele_index.clear();
  hcalhit_eta.clear();
  hcalhit_phi.clear();

  pu_true = -999999.f;
  pu_obs = -999999;
  rho = -999999.f;
  gen_weight = 0.f;

  run = iEvent.eventAuxiliary().run();
  lumi_block = iEvent.eventAuxiliary().luminosityBlock();
  event = iEvent.eventAuxiliary().event();
  bunch_crossing = iEvent.eventAuxiliary().bunchCrossing();
  orbit_number = iEvent.eventAuxiliary().orbitNumber();
  store_number = iEvent.eventAuxiliary().storeNumber();

  /*
  ///getHeep rlated handles
  iEvent.getByToken(vidPassToken_,vidPass);
  iEvent.getByToken(vidBitmapToken_,vidBitmap);
  iEvent.getByToken(vidResultToken_,vidResult);
  iEvent.getByToken(trkIsoMapToken_,trkIsoMap);
  */



  edm::Handle<std::vector<PileupSummaryInfo> > genPileupHandle;
  iEvent.getByToken(puCollection_, genPileupHandle);
  
  if (genPileupHandle.isValid()) {
    for (std::vector<PileupSummaryInfo>::const_iterator pu = genPileupHandle->begin(); pu != genPileupHandle->end(); ++pu) {
      if (pu->getBunchCrossing() == 0) {
        pu_true = pu->getTrueNumInteractions();
        pu_obs = pu->getPU_NumInteractions();

        break;
      }
    }
  }

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken_, rhoHandle);
  if (!rhoHandle.failedToGet())
    rho = *(rhoHandle.product());
  else
    rho = -999999.f;

  edm::Handle<GenEventInfoProduct> genEventHandle;
  iEvent.getByToken(genEventToken_, genEventHandle);
  if (!genEventHandle.failedToGet())
    gen_weight = genEventHandle->weight();
  else
    gen_weight = 0.f;

  edm::Handle<HBHERecHitCollection> hbheRechitsHandle;
  iEvent.getByToken(hbhe_rechits_, hbheRechitsHandle);
  iSetup.get<CaloGeometryRecord>().get(theCaloGeometry);
  iSetup.get<CaloGeometryRecord>().get(towerMap_);

  edm::Handle<std::vector<reco::GenParticle> > genParticlesHandle;
  iEvent.getByToken(genParticlesCollection_, genParticlesHandle);

  //if (iEvent.get(eleToken_).size() > ele_golden.capacity())
  //reallocate_setaddress(iEvent.get(eleToken_).size(), 0);

  //std::cout<<"Outside electron loop "<<std::endl;
  
  for (const auto& ele : iEvent.get(eleToken_)) {
    int genele = 0;
    double ele_min_dr2 = 999999.;
  
    //std::cout<<"Inside electron loop "<<std::endl;
    
    double ele_ptR = 999999.;

    int genpho = 0;
    double pho_min_dr2 = 999999.;
    double pho_ptR = 999999.;

    double this_dr2 = 999999.;

    if (genParticlesHandle.isValid()) {
      for (std::vector<reco::GenParticle>::const_iterator ip = genParticlesHandle->begin(); ip != genParticlesHandle->end(); ++ip) {
	const reco::Candidate *p = (const reco::Candidate*)&(*ip);
	//std::cout << " p->pdgId() " << p->pdgId() << std::endl;
        //if ( std::abs(p->pdgId()) == 11 ) { 
	  //std::cout << "-----  p->status() " << p->status() << " p->pdgId() " << p->pdgId() << std::endl;
	//}

	if ( p->status() == 1 ) {
          if (std::abs(p->pdgId()) == 11 or std::abs(p->pdgId()) == 22)
            this_dr2 = reco::deltaR2(ele,*p);
          else 
            continue;

          if (std::abs(p->pdgId()) == 11 and this_dr2 < ele_min_dr2) {
            ele_min_dr2 = this_dr2;
            ele_ptR = ele.pt() / p->pt();
	  }

          if (std::abs(p->pdgId()) == 22 and this_dr2 < pho_min_dr2) {
            pho_min_dr2 = this_dr2;
            pho_ptR = ele.pt() / p->pt();
	  }
	}
      }
    }
  
    // these cuts were decided looking at min_dr and ptR distributions.
    if ( (ele_min_dr2 < 0.0016) and (ele_ptR > 0.7) && (ele_ptR < 1.3) ) 
      genele = 1;
    if ( (pho_min_dr2 < 0.0016) and (pho_ptR > 0.7) && (pho_ptR < 1.3) ) 
      genpho = 1;

    ele_dR_reco_genele.emplace_back( std::sqrt(ele_min_dr2) );
    ele_rpt_reco_genele.emplace_back(ele_ptR);
    ele_genele.emplace_back(genele);

    ele_dR_reco_genpho.emplace_back( std::sqrt(pho_min_dr2) );
    ele_rpt_reco_genpho.emplace_back(pho_ptR);
    ele_genpho.emplace_back(genpho);

    ele_sc_eta.emplace_back(ele.superCluster()->eta());
    ele_sc_phi.emplace_back(ele.superCluster()->phi());
    ele_pt.emplace_back(ele.pt());
    ele_eta.emplace_back(ele.eta());
    ele_phi.emplace_back(ele.phi());
    ele_sieie_5x5.emplace_back(ele.full5x5_sigmaIetaIeta());
    ele_r9_5x5.emplace_back(ele.full5x5_r9());
  
    reco::GsfElectron::PflowIsolationVariables pfIso = ele.pfIsolationVariables();
    ele_pfiso_pho.emplace_back(pfIso.sumPhotonEt);
    ele_pfiso_neu.emplace_back(pfIso.sumNeutralHadronEt);
    ele_pfiso_cha.emplace_back(pfIso.sumChargedHadronPt);
    ele_pfiso_pu.emplace_back(pfIso.sumPUPt);
    ele_pfiso_hcal.emplace_back(ele.hcalPFClusterIso());
    ele_pfiso_ecal.emplace_back(ele.ecalPFClusterIso());

    ele_detiso03_ecalhit.emplace_back(ele.dr03EcalRecHitSumEt());
    ele_detiso03_hcaltower1.emplace_back(ele.dr03HcalDepth1TowerSumEt());
    ele_detiso03_hcaltower2.emplace_back(ele.dr03HcalDepth2TowerSumEt());
    ele_detiso03_trk.emplace_back(ele.dr03TkSumPt());
    ele_detiso03_trk_heep.emplace_back(ele.dr03TkSumPtHEEP());

    ///fill other ID variables
    ele_dEtaIn.push_back(ele.deltaEtaSuperClusterTrackAtVtx());
    ele_dPhiIn.push_back(ele.deltaPhiSuperClusterTrackAtVtx());
    
    bool isPassVeto   = ele.electronID("cutBasedElectronID-Fall17-94X-V2-veto");
    bool isPassLoose   = ele.electronID("cutBasedElectronID-Fall17-94X-V2-loose");
    bool isPassMedium   = ele.electronID("cutBasedElectronID-Fall17-94X-V2-medium");
    bool isPassTight   = ele.electronID("cutBasedElectronID-Fall17-94X-V2-tight");
    bool isPassMVAiso90 = ele.electronID("mvaEleID-Fall17-iso-V2-wp90");

    bool isPassHEEP = ele.electronID("heepElectronID-HEEPV70");

    ele_IDVeto.push_back((int)isPassVeto);
    ele_IDLoose.push_back((int)isPassLoose);
    ele_IDMedium.push_back((int)isPassMedium);
    ele_IDTight.push_back((int)isPassTight);
    ele_IDMVAiso90.push_back((int)isPassMVAiso90);
    ele_IDHEEP.push_back((int)isPassHEEP);
    
    //ele_IDbits.push_back(0);
    
    ele_IDbits.push_back({          
	ele.userInt("cutBasedElectronID-Fall17-94X-V2-veto"),
	  ele.userInt("cutBasedElectronID-Fall17-94X-V2-loose"),
	  ele.userInt("cutBasedElectronID-Fall17-94X-V2-medium"),
	  ele.userInt("cutBasedElectronID-Fall17-94X-V2-tight"),
	  ele.userInt("heepElectronID-HEEPV70")}
	  //ele.userInt("heepElectronID-HEEPV70Bitmap")}
      );
    
    //std::cout<<"n_ele : ID bit for HEEP  : from the array : elePt : "<<n_ele<<" "<<ele.userInt("heepElectronID-HEEPV70")<<" "<<ele_IDbits[n_ele][4]<<" "<<ele_pt[n_ele]<<std::endl;

    ele_missingHit.push_back(ele.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS));
    ele_nHit.push_back(ele.gsfTrack()->hitPattern().trackerLayersWithMeasurement());
    ele_gsfTrackChi2.push_back(ele.gsfTrack()->normalizedChi2());
    ele_eSCoP.push_back(ele.eSuperClusterOverP());
    ele_eSCoPout.push_back(ele.eEleClusterOverPout());
    ele_1oEm1op.push_back( (1.0/ele.ecalEnergy()) - (1.0/ele.trackMomentumAtVtx().r()) );
    ele_psEorawE.push_back( ele.superCluster()->preshowerEnergy()/ele.superCluster()->rawEnergy() );
    ele_isEcalDriven.push_back(ele.ecalDrivenSeed());

  
    //std::cout<<"SC eta : phi : IDveto :  IDbits "<<ele_sc_eta[n_ele]<<" "<<ele_sc_phi[n_ele]<<" "<<ele_IDVeto[n_ele]<<" "<<ele_IDbits[n_ele][4]<<std::endl;
    //ele_convVtxFitProb.push_back();


    EcalClusterLazyTools lazyTool(iEvent, iSetup, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);

    const reco::SuperCluster& superClus = *ele.superCluster();
    const reco::CaloCluster &seedCluster = *superClus.seed();
    DetId seedId = seedCluster.seed() ;
    ele_seed_detid.emplace_back(seedId.det());
    ele_seed_subdetid.emplace_back(seedId.subdetId());

    float var_ele_seed_eta = -999999.f;
    float var_ele_seed_phi = -999999.f;

    DetId seed = (seedCluster.hitsAndFractions())[0].first;
    bool isBarrel = seed.subdetId() == EcalBarrel;
    const EcalRecHitCollection * rechits = isBarrel ? lazyTool.getEcalEBRecHitCollection() : lazyTool.getEcalEERecHitCollection();
    EcalRecHitCollection::const_iterator theSeedHit = rechits->find(seed);
    if (theSeedHit != rechits->end()) {
      if ( (theSeedHit->id().rawId() != 0 ) ) {
	if (theCaloGeometry.product() != nullptr) {
	  const CaloSubdetectorGeometry *ecalgeo = theCaloGeometry->getSubdetectorGeometry(theSeedHit->id());
	  
	  if(ecalgeo->getGeometry(theSeedHit->id()) !=nullptr){
	    const GlobalPoint & ecalrechitPoint = (theCaloGeometry.product())->getPosition(theSeedHit->id());
	    var_ele_seed_eta=ecalrechitPoint.eta();
	    var_ele_seed_phi=ecalrechitPoint.phi();
	  }
	}
      }
    }  
    ele_seed_eta.emplace_back(var_ele_seed_eta);
    ele_seed_phi.emplace_back(var_ele_seed_phi);

    int var_ele_seed_ieta = -999999;
    int var_ele_seed_iphi = -999999;
    int var_ele_seed_raw_id = -999999;
    
    int var_ele_seed_hcal_ieta = -999999;
    int var_ele_seed_hcal_iphi = -999999;

    if ( seedId.det() == DetId::Ecal ) {
      if (seedId.subdetId() == EcalBarrel) {
	EBDetId ebId(seedId);
	var_ele_seed_ieta = ebId.ieta();
	var_ele_seed_iphi = ebId.iphi();
	var_ele_seed_raw_id = ebId.rawId();       
      }
      else if (seedId.subdetId() == EcalEndcap) {
	EEDetId eeId(seedId);
	var_ele_seed_ieta = eeId.ix();
        var_ele_seed_iphi = eeId.iy();
        var_ele_seed_raw_id = eeId.rawId();
      }

      // get hold of the seed hcal behind ele seed
      CaloTowerDetId towerId(towerMap_->towerOf(seedId));       
      var_ele_seed_hcal_ieta = towerId.ieta();
      var_ele_seed_hcal_iphi = towerId.iphi();
    }
  
    ele_track_fbrem.emplace_back(ele.trackFbrem());
    ele_sc_fbrem.emplace_back(ele.superClusterFbrem());
    ele_nbrem.emplace_back(ele.numberOfBrems());

    int var_golden = 0;
    int var_unknown = 0;
    int var_gap = 0;
    int var_badtrack = 0;
    int var_showering = 0;
    int var_bigbrem = 0;

    if (ele.classification() == reco::GsfElectron::GOLDEN)
      var_golden = 1;

    if (ele.classification() == reco::GsfElectron::UNKNOWN)
      var_unknown = 1;

    if (ele.classification() == reco::GsfElectron::BIGBREM)
      var_bigbrem = 1;

    if (ele.classification() == reco::GsfElectron::BADTRACK)
      var_badtrack = 1;

    if (ele.classification() == reco::GsfElectron::SHOWERING)
      var_showering = 1;

    if (ele.classification() == reco::GsfElectron::GAP)
      var_gap = 1;

    ele_golden.emplace_back(var_golden);
    ele_unknown.emplace_back(var_unknown);
    ele_gap.emplace_back(var_gap);
    ele_badtrack.emplace_back(var_badtrack);
    ele_showering.emplace_back(var_showering);
    ele_bigbrem.emplace_back(var_bigbrem);

    ele_sc_energy.emplace_back(superClus.energy());
    ele_sc_raw_energy.emplace_back(superClus.rawEnergy());
    ele_seed_energy.emplace_back(seedCluster.energy());
    ele_seed_corr_energy.emplace_back(seedCluster.correctedEnergy());
    ele_ecal_energy.emplace_back(ele.ecalEnergy());
    ele_cmssw_hoe.emplace_back(ele.hcalOverEcal());
    ele_cmssw_hoe_tower.emplace_back(ele.hcalOverEcalBc());
    ele_cmssw_hoe_5x5.emplace_back(ele.full5x5_hcalOverEcal());
   
    ele_seed_ieta.emplace_back(var_ele_seed_ieta);
    ele_seed_iphi.emplace_back(var_ele_seed_iphi);
    ele_seed_raw_id.emplace_back(var_ele_seed_raw_id);

    ele_seed_hcal_ieta.emplace_back(var_ele_seed_hcal_ieta);
    ele_seed_hcal_iphi.emplace_back(var_ele_seed_hcal_iphi);

    ele_eb.emplace_back(ele.isEB());
    ele_ee.emplace_back(ele.isEE());
    ele_gap_eb_ee.emplace_back(ele.isEBEEGap());
    ele_gap_eb_eta.emplace_back(ele.isEBEtaGap());
    ele_gap_eb_phi.emplace_back(ele.isEBPhiGap());
    ele_gap_ee_dee.emplace_back(ele.isEEDeeGap());
    ele_gap_ee_ring.emplace_back(ele.isEERingGap());

    ++n_ele;
  }

  // given the context, should be ok...
  if (n_ele == 0)
    return;

  // just in case
  assert(((void) "ERROR: ele_eb size doesn't match n_ele!!!", int(ele_eb.size()) == n_ele));
  assert(((void) "ERROR: ele_ee size doesn't match n_ele!!!", int(ele_ee.size()) == n_ele));
  assert(((void) "ERROR: ele_gap_eb_ee size doesn't match n_ele!!!", int(ele_gap_eb_ee.size()) == n_ele));
  assert(((void) "ERROR: ele_gab_eb_eta size doesn't match n_ele!!!", int(ele_gap_eb_eta.size()) == n_ele));
  assert(((void) "ERROR: ele_gap_eb_phi size doesn't match n_ele!!!", int(ele_gap_eb_phi.size()) == n_ele));
  assert(((void) "ERROR: ele_gap_ee_dee size doesn't match n_ele!!!", int(ele_gap_ee_dee.size()) == n_ele));
  assert(((void) "ERROR: ele_gap_ee_ring size doesn't match n_ele!!!", int(ele_gap_ee_ring.size()) == n_ele));

  assert(((void) "ERROR: ele_golden size doesn't match n_ele!!!", int(ele_golden.size()) == n_ele));
  assert(((void) "ERROR: ele_unknown size doesn't match n_ele!!!", int(ele_unknown.size()) == n_ele));
  assert(((void) "ERROR: ele_badtrack size doesn't match n_ele!!!", int(ele_badtrack.size()) == n_ele));
  assert(((void) "ERROR: ele_bigbrem size doesn't match n_ele!!!", int(ele_bigbrem.size()) == n_ele));
  assert(((void) "ERROR: ele_showering size doesn't match n_ele!!!", int(ele_showering.size()) == n_ele));
  assert(((void) "ERROR: ele_gap size doesn't match n_ele!!!", int(ele_gap.size()) == n_ele));

  assert(((void) "ERROR: ele_track_fbrem size doesn't match n_ele!!!", int(ele_track_fbrem.size()) == n_ele));
  assert(((void) "ERROR: ele_sc_fbrem size doesn't match n_ele!!!", int(ele_sc_fbrem.size()) == n_ele));
  assert(((void) "ERROR: ele_nbrem size doesn't match n_ele!!!", int(ele_nbrem.size()) == n_ele));

  assert(((void) "ERROR: ele_genele size doesn't match n_ele!!!", int(ele_genele.size()) == n_ele));
  assert(((void) "ERROR: ele_rpt_reco_genele size doesn't match n_ele!!!", int(ele_rpt_reco_genele.size()) == n_ele));
  assert(((void) "ERROR: ele_dR_reco_genele size doesn't match n_ele!!!", int(ele_dR_reco_genele.size()) == n_ele));

  assert(((void) "ERROR: ele_genpho size doesn't match n_ele!!!", int(ele_genpho.size()) == n_ele));
  assert(((void) "ERROR: ele_rpt_reco_genpho size doesn't match n_ele!!!", int(ele_rpt_reco_genpho.size()) == n_ele));
  assert(((void) "ERROR: ele_dR_reco_genpho size doesn't match n_ele!!!", int(ele_dR_reco_genpho.size()) == n_ele));

  assert(((void) "ERROR: ele_sc_energy size doesn't match n_ele!!!", int(ele_sc_energy.size()) == n_ele));
  assert(((void) "ERROR: ele_sc_raw_energy size doesn't match n_ele!!!", int(ele_sc_raw_energy.size()) == n_ele));
  assert(((void) "ERROR: ele_ecal_energy size doesn't match n_ele!!!", int(ele_ecal_energy.size()) == n_ele));
  assert(((void) "ERROR: ele_seed_energy size doesn't match n_ele!!!", int(ele_seed_energy.size()) == n_ele));
  assert(((void) "ERROR: ele_seed_corr_energy size doesn't match n_ele!!!", int(ele_seed_corr_energy.size()) == n_ele));
  assert(((void) "ERROR: ele_cmssw_hoe size doesn't match n_ele!!!", int(ele_cmssw_hoe.size()) == n_ele));
  assert(((void) "ERROR: ele_cmssw_hoe_tower size doesn't match n_ele!!!", int(ele_cmssw_hoe_tower.size()) == n_ele));
  assert(((void) "ERROR: ele_cmssw_hoe_5x5 size doesn't match n_ele!!!", int(ele_cmssw_hoe_5x5.size()) == n_ele));
  assert(((void) "ERROR: ele_sc_eta size doesn't match n_ele!!!", int(ele_sc_eta.size()) == n_ele));
  assert(((void) "ERROR: ele_sc_phi size doesn't match n_ele!!!", int(ele_sc_phi.size()) == n_ele));
  assert(((void) "ERROR: ele_pt size doesn't match n_ele!!!", int(ele_pt.size()) == n_ele));
  assert(((void) "ERROR: ele_eta size doesn't match n_ele!!!", int(ele_eta.size()) == n_ele));
  assert(((void) "ERROR: ele_phi size doesn't match n_ele!!!", int(ele_phi.size()) == n_ele));
  assert(((void) "ERROR: ele_sieie_5x5 size doesn't match n_ele!!!", int(ele_sieie_5x5.size()) == n_ele));
  assert(((void) "ERROR: ele_r9_5x5 size doesn't match n_ele!!!", int(ele_r9_5x5.size()) == n_ele));

  assert(((void) "ERROR: ele_pfiso_pho size doesn't match n_ele!!!", int(ele_pfiso_pho.size()) == n_ele));
  assert(((void) "ERROR: ele_pfiso_neu size doesn't match n_ele!!!", int(ele_pfiso_neu.size()) == n_ele));
  assert(((void) "ERROR: ele_pfiso_cha size doesn't match n_ele!!!", int(ele_pfiso_cha.size()) == n_ele));
  assert(((void) "ERROR: ele_pfiso_pu size doesn't match n_ele!!!", int(ele_pfiso_pu.size()) == n_ele));
  assert(((void) "ERROR: ele_pfiso_hcal size doesn't match n_ele!!!", int(ele_pfiso_hcal.size()) == n_ele));
  assert(((void) "ERROR: ele_pfiso_ecal size doesn't match n_ele!!!", int(ele_pfiso_ecal.size()) == n_ele));

  assert(((void) "ERROR: ele_detiso03_ecalhit size doesn't match n_ele!!!", int(ele_detiso03_ecalhit.size()) == n_ele));
  assert(((void) "ERROR: ele_detiso03_hcaltower1 size doesn't match n_ele!!!", int(ele_detiso03_hcaltower1.size()) == n_ele));
  assert(((void) "ERROR: ele_detiso03_hcaltower2 size doesn't match n_ele!!!", int(ele_detiso03_hcaltower2.size()) == n_ele));
  assert(((void) "ERROR: ele_detiso03_trk size doesn't match n_ele!!!", int(ele_detiso03_trk.size()) == n_ele));
  assert(((void) "ERROR: ele_detiso03_trk_heep size doesn't match n_ele!!!", int(ele_detiso03_trk_heep.size()) == n_ele));

  assert(((void) "ERROR: ele_seed_detid size doesn't match n_ele!!!", int(ele_seed_detid.size()) == n_ele));
  assert(((void) "ERROR: ele_seed_subdetid size doesn't match n_ele!!!", int(ele_seed_subdetid.size()) == n_ele));
  assert(((void) "ERROR: ele_seed_ieta size doesn't match n_ele!!!", int(ele_seed_ieta.size()) == n_ele));
  assert(((void) "ERROR: ele_seed_iphi size doesn't match n_ele!!!", int(ele_seed_iphi.size()) == n_ele));
  assert(((void) "ERROR: ele_seed_eta size doesn't match n_ele!!!", int(ele_seed_eta.size()) == n_ele));
  assert(((void) "ERROR: ele_seed_phi size doesn't match n_ele!!!", int(ele_seed_phi.size()) == n_ele));
  assert(((void) "ERROR: ele_seed_raw_id size doesn't match n_ele!!!", int(ele_seed_raw_id.size()) == n_ele));
  assert(((void) "ERROR: ele_seed_hcal_ieta size doesn't match n_ele!!!", int(ele_seed_hcal_ieta.size()) == n_ele));
  assert(((void) "ERROR: ele_seed_hcal_iphi size doesn't match n_ele!!!", int(ele_seed_hcal_iphi.size()) == n_ele));

  //if (iEvent.get(hbhe_rechits_).size() > hcalhit_depth.capacity())
  //reallocate_setaddress(0, iEvent.get(hbhe_rechits_).size());

  for (auto& hcalrh : iEvent.get(hbhe_rechits_)) {
    if (hcalrh.energy() < getMinEnergyHCAL(hcalrh.id()))
      continue;
    if (hcalrh.id().depth() < 0 or hcalrh.id().depth() > 100)
      continue;

    hcalhit_ieta.emplace_back(hcalrh.id().ieta());
    hcalhit_iphi.emplace_back(hcalrh.id().iphi());
    hcalhit_energy.emplace_back(hcalrh.energy());

    hcalhit_raw_id.emplace_back(hcalrh.id().rawId());
    hcalhit_depth.emplace_back(hcalrh.id().depth());

    float rechitEta = -999999.f;
    float rechitPhi = -999999.f;
    if (hcalrh.id().rawId() != 0) {
      if (theCaloGeometry.product() != nullptr) {
        const CaloSubdetectorGeometry *geo = theCaloGeometry->getSubdetectorGeometry(hcalrh.id());

        if (geo->getGeometry(hcalrh.id()) != nullptr) {
          const GlobalPoint & rechitPoint = theCaloGeometry.product()->getPosition(hcalrh.id());

          rechitEta=rechitPoint.eta();
          rechitPhi=rechitPoint.phi();	
        }
      }
    }
    hcalhit_eta.emplace_back(rechitEta);
    hcalhit_phi.emplace_back(rechitPhi);

    imin = -1;
    min_dieta = 999999;
    min_diphi = 999999;
    min_diR2 = 999999.f;

    for (int iE = 0; iE < n_ele; ++iE) {
      int dieta = calDIEta(ele_seed_hcal_ieta[iE], hcalhit_ieta.back());
      int diphi = calDIPhi(ele_seed_hcal_iphi[iE], hcalhit_iphi.back());
      float diR2 = float(dieta * dieta) + float(diphi * diphi);

      if (diR2 < min_diR2) {
        min_diR2 = diR2;
        min_dieta = dieta;
        min_diphi = diphi;
        imin = iE;
      }
    }

    hcalhit_seed_dieta.emplace_back(min_dieta);
    hcalhit_seed_diphi.emplace_back(min_diphi);
    hcalhit_ele_index.emplace_back(imin);

    ++n_hcalhit;
  }

  // again, just in case
  assert(((void) "ERROR: hcalhit_ieta size doesn't match n_hcalhit!!!", int(hcalhit_ieta.size()) == n_hcalhit));
  assert(((void) "ERROR: hcalhit_iphi size doesn't match n_hcalhit!!!", int(hcalhit_iphi.size()) == n_hcalhit));
  assert(((void) "ERROR: hcalhit_energy size doesn't match n_hcalhit!!!", int(hcalhit_energy.size()) == n_hcalhit));
  assert(((void) "ERROR: hcalhit_seed_dieta size doesn't match n_hcalhit!!!", int(hcalhit_seed_dieta.size()) == n_hcalhit));
  assert(((void) "ERROR: hcalhit_seed_diphi size doesn't match n_hcalhit!!!", int(hcalhit_seed_diphi.size()) == n_hcalhit));
  assert(((void) "ERROR: hcalhit_raw_id size doesn't match n_hcalhit!!!", int(hcalhit_raw_id.size()) == n_hcalhit));
  assert(((void) "ERROR: hcalhit_depth size doesn't match n_hcalhit!!!", int(hcalhit_depth.size()) == n_hcalhit));
  assert(((void) "ERROR: hcalhit_ele_index size doesn't match n_hcalhit!!!", int(hcalhit_ele_index.size()) == n_hcalhit));
  assert(((void) "ERROR: hcalhit_eta size doesn't match n_hcalhit!!!", int(hcalhit_eta.size()) == n_hcalhit));
  assert(((void) "ERROR: hcalhit_phi size doesn't match n_hcalhit!!!", int(hcalhit_phi.size()) == n_hcalhit));

  tree->Fill();

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}

// doing some blatant copy paste from RecoEgamma/EgammaIsolationAlgos/src/EGHcalRecHitSelector.cc
int FlatEleHoEAnalyzer::calDIPhi(int iPhi1, int iPhi2) {
  int dPhi = iPhi1 - iPhi2;
  if (dPhi > 72 / 2)
    dPhi -= 72;
  else if (dPhi < -72 / 2)
    dPhi += 72;
  return dPhi;
}

int FlatEleHoEAnalyzer::calDIEta(int iEta1, int iEta2) {
  int dEta = iEta1 - iEta2;
  if (iEta1 * iEta2 < 0) {  //-ve to +ve transition and no crystal at zero
    if (dEta < 0)
      ++dEta;
    else
      --dEta;
  }

  return dEta;
}

// HCAL thresholds from here https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/CaloTowersCreator/python/calotowermaker_cfi.py?%21v=CMSSW_10_6_2
// Note: As far as I understood, 
// for 2018, HB threshold is 0.7, and for Run 3 it becomes 0.1 in depth1, 0.2 in depth2, 0.3 in other depths.
// In HE, 2018 and Run3 is same, and it is 0.1 in depth1, and 0.2 in other depths.
// Double check these HCAL thresholds from Sam.
float FlatEleHoEAnalyzer::getMinEnergyHCAL(HcalDetId id) const {
  if ( (id.subdetId() == HcalBarrel)  ) {
    if ( (Run2_2018 == 1) )
      return 0.7f;
    else if ( (Run2_2018 == 0) ) { // means Run3
      if (id.depth() == 1)
	return 0.1f;
      else if (id.depth() == 2)
	return 0.2f;
      else
	return 0.3f;
    }
    else // neither 2018, nor Run3, not supported
      return 999999.f;
  } 
  else if (id.subdetId() == HcalEndcap) {
    if (id.depth() == 1)
      return 0.1f;
    else
      return 0.2f;
  } else
    return 999999.f;
}

void FlatEleHoEAnalyzer::reallocate_setaddress(int n_ele_, int n_hcalhit_)
{

  /*
    static int cap_ele = 8;
  cap_ele = (n_ele_ == 0) ? cap_ele : n_ele_;
 
  ele_eb.reserve(cap_ele);
  ele_ee.reserve(cap_ele);
  ele_gap_eb_ee.reserve(cap_ele);
  ele_gap_eb_eta.reserve(cap_ele);
  ele_gap_eb_phi.reserve(cap_ele);
  ele_gap_ee_dee.reserve(cap_ele);
  ele_gap_ee_ring.reserve(cap_ele);

  ele_golden.reserve(cap_ele);
  ele_unknown.reserve(cap_ele);
  ele_badtrack.reserve(cap_ele);
  ele_bigbrem.reserve(cap_ele);
  ele_showering.reserve(cap_ele);
  ele_gap.reserve(cap_ele);

  ele_track_fbrem.reserve(cap_ele);
  ele_sc_fbrem.reserve(cap_ele);
  ele_nbrem.reserve(cap_ele);

  ele_genele.reserve(cap_ele);
  ele_rpt_reco_genele.reserve(cap_ele);
  ele_dR_reco_genele.reserve(cap_ele);

  ele_genpho.reserve(cap_ele);
  ele_rpt_reco_genpho.reserve(cap_ele);
  ele_dR_reco_genpho.reserve(cap_ele);

  ele_sc_energy.reserve(cap_ele);
  ele_sc_raw_energy.reserve(cap_ele);
  ele_ecal_energy.reserve(cap_ele);
  ele_seed_energy.reserve(cap_ele);
  ele_seed_corr_energy.reserve(cap_ele);
  ele_cmssw_hoe.reserve(cap_ele);
  ele_cmssw_hoe_tower.reserve(cap_ele);
  ele_cmssw_hoe_5x5.reserve(cap_ele);
  ele_sc_eta.reserve(cap_ele);
  ele_sc_phi.reserve(cap_ele);
  ele_pt.reserve(cap_ele);
  ele_eta.reserve(cap_ele);
  ele_phi.reserve(cap_ele);
  ele_sieie_5x5.reserve(cap_ele);
  ele_r9_5x5.reserve(cap_ele);

  ele_pfiso_pho.reserve(cap_ele);
  ele_pfiso_neu.reserve(cap_ele);
  ele_pfiso_cha.reserve(cap_ele);
  ele_pfiso_pu.reserve(cap_ele);
  ele_pfiso_hcal.reserve(cap_ele);
  ele_pfiso_ecal.reserve(cap_ele);

  ele_detiso03_ecalhit.reserve(cap_ele);
  ele_detiso03_hcaltower1.reserve(cap_ele);
  ele_detiso03_hcaltower2.reserve(cap_ele);
  ele_detiso03_trk.reserve(cap_ele);
  ele_detiso03_trk_heep.reserve(cap_ele);

  ele_seed_detid.reserve(cap_ele);
  ele_seed_subdetid.reserve(cap_ele);
  ele_seed_ieta.reserve(cap_ele);
  ele_seed_iphi.reserve(cap_ele);
  ele_seed_eta.reserve(cap_ele);
  ele_seed_phi.reserve(cap_ele);
  ele_seed_raw_id.reserve(cap_ele);
  ele_seed_hcal_ieta.reserve(cap_ele);
  ele_seed_hcal_iphi.reserve(cap_ele);

  ele_dEtaIn.reserve(cap_ele);
  ele_dPhiIn.reserve(cap_ele);
  ele_psEorawE.reserve(cap_ele);
  ele_1oEm1op.reserve(cap_ele);
  ele_eSCoP.reserve(cap_ele);
  ele_eSCoPout.reserve(cap_ele);
  //ele_convVtxFitProb.reserve(cap_ele);
  ele_gsfTrackChi2.reserve(cap_ele);
  ele_nHit.reserve(cap_ele);
  ele_missingHit.reserve(cap_ele);
  ele_IDVeto.reserve(cap_ele);
  ele_IDLoose.reserve(cap_ele);
  ele_IDMedium.reserve(cap_ele);
  ele_IDTight.reserve(cap_ele);
  ele_IDMVAiso90.reserve(cap_ele);
  
  //for(int iele=0; iele<cap_ele; iele++)
  //ele_IDbits[iele].reserve(5);

  ele_isEcalDriven.reserve(cap_ele);
  ele_IDHEEP.reserve(cap_ele);

  

  //static int cap_hcalhit = 128;
  static int cap_hcalhit = 500;
  cap_hcalhit = (n_hcalhit_ == 0) ? cap_hcalhit : n_hcalhit_;
  hcalhit_ieta.reserve(cap_hcalhit);
  hcalhit_iphi.reserve(cap_hcalhit);
  hcalhit_energy.reserve(cap_hcalhit);
  hcalhit_seed_dieta.reserve(cap_hcalhit);
  hcalhit_seed_diphi.reserve(cap_hcalhit);
  hcalhit_raw_id.reserve(cap_hcalhit);
  hcalhit_depth.reserve(cap_hcalhit);
  hcalhit_ele_index.reserve(cap_hcalhit);
  hcalhit_eta.reserve(cap_hcalhit);
  hcalhit_phi.reserve(cap_hcalhit);
  */

  
  //if (n_ele_ == 0 and n_hcalhit_ == 0) {

  /*
  if (n_ele_ != 0) {
  //if(1) {
    //std::cout << "Electron block realloc to " << ele_golden.capacity() << "..." << std::endl;


    //b_ele_IDbits->SetAddress(ele_IDbits.data());
    b_ele_IDVeto->SetAddress(ele_IDVeto.data());
    b_ele_IDLoose->SetAddress(ele_IDLoose.data());
    b_ele_IDMedium->SetAddress(ele_IDMedium.data());
    b_ele_IDTight->SetAddress(ele_IDTight.data());
    b_ele_IDMVAiso90->SetAddress(ele_IDMVAiso90.data());
    b_ele_IDHEEP->SetAddress(ele_IDHEEP.data());
    b_ele_dEtaIn->SetAddress(ele_dEtaIn.data());
    b_ele_dPhiIn->SetAddress(ele_dPhiIn.data());
    b_ele_psEorawE->SetAddress(ele_psEorawE.data());
    b_ele_1oEm1op->SetAddress(ele_1oEm1op.data());
    b_ele_eSCoP->SetAddress(ele_eSCoP.data());
    b_ele_eSCoPout->SetAddress(ele_eSCoPout.data());
    b_ele_gsfTrackChi2->SetAddress(ele_gsfTrackChi2.data());
    b_ele_nHit->SetAddress(ele_nHit.data());
    b_ele_missingHit->SetAddress(ele_missingHit.data());
    b_ele_isEcalDriven->SetAddress(ele_isEcalDriven.data());

    b_ele_eb->SetAddress(ele_eb.data());
    b_ele_ee->SetAddress(ele_ee.data());
    b_ele_gap_eb_ee->SetAddress(ele_gap_eb_ee.data());
    b_ele_gap_eb_eta->SetAddress(ele_gap_eb_eta.data());
    b_ele_gap_eb_phi->SetAddress(ele_gap_eb_phi.data());
    b_ele_gap_ee_dee->SetAddress(ele_gap_ee_dee.data());
    b_ele_gap_ee_ring->SetAddress(ele_gap_ee_ring.data());

    b_ele_golden->SetAddress(ele_golden.data());
    b_ele_unknown->SetAddress(ele_unknown.data());
    b_ele_bigbrem->SetAddress(ele_bigbrem.data());
    b_ele_gap->SetAddress(ele_gap.data());
    b_ele_badtrack->SetAddress(ele_badtrack.data());
    b_ele_showering->SetAddress(ele_showering.data());

    b_ele_track_fbrem->SetAddress(ele_track_fbrem.data());
    b_ele_sc_fbrem->SetAddress(ele_sc_fbrem.data());
    b_ele_nbrem->SetAddress(ele_nbrem.data());

    b_ele_genele->SetAddress(ele_genele.data());
    b_ele_dR_reco_genele->SetAddress(ele_dR_reco_genele.data());
    b_ele_rpt_reco_genele->SetAddress(ele_rpt_reco_genele.data());

    b_ele_genpho->SetAddress(ele_genpho.data());
    b_ele_dR_reco_genpho->SetAddress(ele_dR_reco_genpho.data());
    b_ele_rpt_reco_genpho->SetAddress(ele_rpt_reco_genpho.data());

    b_ele_sc_energy->SetAddress(ele_sc_energy.data());
    b_ele_sc_raw_energy->SetAddress(ele_sc_raw_energy.data());
    b_ele_ecal_energy->SetAddress(ele_ecal_energy.data());
    b_ele_seed_energy->SetAddress(ele_seed_energy.data());
    b_ele_seed_corr_energy->SetAddress(ele_seed_corr_energy.data());
    b_ele_cmssw_hoe->SetAddress(ele_cmssw_hoe.data());
    b_ele_cmssw_hoe_tower->SetAddress(ele_cmssw_hoe_tower.data());
    b_ele_cmssw_hoe_5x5->SetAddress(ele_cmssw_hoe_5x5.data());
    b_ele_sc_eta->SetAddress(ele_sc_eta.data());
    b_ele_sc_phi->SetAddress(ele_sc_phi.data());
    b_ele_pt->SetAddress(ele_pt.data());
    b_ele_eta->SetAddress(ele_eta.data());
    b_ele_phi->SetAddress(ele_phi.data());
    b_ele_sieie_5x5->SetAddress(ele_sieie_5x5.data());
    b_ele_r9_5x5->SetAddress(ele_r9_5x5.data());

    b_ele_pfiso_pho->SetAddress(ele_pfiso_pho.data());
    b_ele_pfiso_neu->SetAddress(ele_pfiso_neu.data());
    b_ele_pfiso_cha->SetAddress(ele_pfiso_cha.data());
    b_ele_pfiso_pu->SetAddress(ele_pfiso_pu.data());
    b_ele_pfiso_hcal->SetAddress(ele_pfiso_hcal.data());
    b_ele_pfiso_ecal->SetAddress(ele_pfiso_ecal.data());

    b_ele_detiso03_ecalhit->SetAddress(ele_detiso03_ecalhit.data());
    b_ele_detiso03_hcaltower1->SetAddress(ele_detiso03_hcaltower1.data());
    b_ele_detiso03_hcaltower2->SetAddress(ele_detiso03_hcaltower2.data());
    b_ele_detiso03_trk->SetAddress(ele_detiso03_trk.data());
    b_ele_detiso03_trk_heep->SetAddress(ele_detiso03_trk_heep.data());

    b_ele_seed_detid->SetAddress(ele_seed_detid.data());
    b_ele_seed_subdetid->SetAddress(ele_seed_subdetid.data());
    b_ele_seed_ieta->SetAddress(ele_seed_ieta.data());
    b_ele_seed_iphi->SetAddress(ele_seed_iphi.data());
    b_ele_seed_eta->SetAddress(ele_seed_eta.data());
    b_ele_seed_phi->SetAddress(ele_seed_phi.data());
    b_ele_seed_raw_id->SetAddress(ele_seed_raw_id.data());
    b_ele_seed_hcal_ieta->SetAddress(ele_seed_hcal_ieta.data());
    b_ele_seed_hcal_iphi->SetAddress(ele_seed_hcal_iphi.data());
  }

  if (n_hcalhit_ != 0) {
    //std::cout << "Hcalhit block realloc to " << hcalhit_depth.capacity() << "..." << std::endl;

    b_hcalhit_ieta->SetAddress(hcalhit_ieta.data());
    b_hcalhit_iphi->SetAddress(hcalhit_iphi.data());
    b_hcalhit_energy->SetAddress(hcalhit_energy.data());
    b_hcalhit_seed_dieta->SetAddress(hcalhit_seed_dieta.data());
    b_hcalhit_seed_diphi->SetAddress(hcalhit_seed_diphi.data());
    b_hcalhit_raw_id->SetAddress(hcalhit_raw_id.data());
    b_hcalhit_depth->SetAddress(hcalhit_depth.data());
    b_hcalhit_ele_index->SetAddress(hcalhit_ele_index.data());
    b_hcalhit_eta->SetAddress(hcalhit_eta.data());
    b_hcalhit_phi->SetAddress(hcalhit_phi.data());
  }
  */
}


// ------------ method called once each job just before starting event loop  ------------
void
FlatEleHoEAnalyzer::beginJob()
{
  //tree->SetAutoSave(0);

  //tree->SetImplicitMT(false);
  //reallocate_setaddress();
}

// ------------ method called once each job just after ending the event loop  ------------
void
FlatEleHoEAnalyzer::endJob()
{
  file->cd();
  tree->Write();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FlatEleHoEAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FlatEleHoEAnalyzer);
