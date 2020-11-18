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
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.
class FlatPhoHoEAnalyzer : public edm::one::EDAnalyzer<>  {
public:
  explicit FlatPhoHoEAnalyzer(const edm::ParameterSet&);
  ~FlatPhoHoEAnalyzer();
  
  static edm::ParameterSetDescription makePSetDescription();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  TFile *file;
  TTree *tree;

  int run, lumi_block, event, bunch_crossing, orbit_number, store_number;

  int n_pho;
  std::vector<int>  pho_haspixseed;

  

  std::vector<int>    pho_genpho;
  std::vector<float>  pho_dR_reco_genpho;
  std::vector<float>  pho_rpt_reco_genpho;

  std::vector<float>  pho_sc_energy;
  std::vector<float>  pho_sc_raw_energy;
  std::vector<float>  pho_seed_energy;
  std::vector<float>  pho_seed_corr_energy;
  std::vector<float>  pho_cmssw_hoe;
  std::vector<float>  pho_sc_eta;
  std::vector<float>  pho_sc_phi;
  std::vector<float>  pho_pt;
  std::vector<float>  pho_eta;
  std::vector<float>  pho_phi;
  std::vector<float>  pho_sieie_5x5;
  std::vector<float>  pho_r9_5x5;

  std::vector<float>  pho_pfiso_pho;
  std::vector<float>  pho_pfiso_neu;
  std::vector<float>  pho_pfiso_cha;
  
  
  std::vector<int>    pho_seed_detid;
  std::vector<int>    pho_seed_subdetid;
  std::vector<int>    pho_seed_ieta;
  std::vector<int>    pho_seed_iphi;
  std::vector<float>  pho_seed_eta;
  std::vector<float>  pho_seed_phi;

  std::vector<int>    pho_nBCs;
  std::vector<std::vector<float> >  pho_seedBC_eta;
  std::vector<std::vector<float> >  pho_seedBC_phi;

  std::vector<int>    pho_seed_raw_id;
  std::vector<int>    pho_seed_hcal_ieta;
  std::vector<int>    pho_seed_hcal_iphi;

  std::vector<int>  pho_IDLoose;
  std::vector<int>  pho_IDMedium;
  std::vector<int>  pho_IDTight;

  
  int n_hcalhit;
  std::vector<int>    hcalhit_ieta;
  std::vector<int>    hcalhit_iphi;
  std::vector<float>  hcalhit_energy;
  std::vector<int>    hcalhit_seed_dieta;
  std::vector<int>    hcalhit_seed_diphi;
  std::vector<int>    hcalhit_raw_id;
  std::vector<int>    hcalhit_depth;
  std::vector<int>    hcalhit_pho_index;
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
  void reallocate_setaddress(int n_pho_ = 0, int n_hcalhit_ = 0);
  float getMinEnergyHCAL(HcalDetId id) const;
  bool GetSeedEtaPhi(bool isBarrel, DetId seedId, EcalClusterLazyTools lazyTool, double &var_pho_seed_eta, double &var_pho_seed_phi);
  

  int maxDIEta_ = 5;
  int maxDIPhi_ = 5;

  // ----------member data ---------------------------
  //edm::EDGetTokenT<edm::View<reco::GsfPhoctron> > phoToken_;

  //edm::EDGetTokenT<edm::View<pat::Electron> > phoToken_;

  edm::EDGetTokenT<edm::View<pat::Photon> > phoToken_;

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

  edm::EDGetTokenT<edm::ValueMap<float> > phoPhoIsoToken_; 
  edm::EDGetTokenT<edm::ValueMap<float> > phoNeuIsoToken_; 
  edm::EDGetTokenT<edm::ValueMap<float> > phoChaIsoToken_; 

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
FlatPhoHoEAnalyzer::FlatPhoHoEAnalyzer(const edm::ParameterSet& iConfig) :
  //phoToken_(consumes<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("electrons"))),

  phoToken_(consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photons"))),

  puCollection_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupCollection"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoSrc"))),
  genEventToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventSrc"))),
  hbhe_rechits_(consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("hbheInput"))),
  ebReducedRecHitCollection_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ebReducedRecHitCollection"))),
  eeReducedRecHitCollection_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("eeReducedRecHitCollection"))),
  esReducedRecHitCollection_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("esReducedRecHitCollection"))),
  genParticlesCollection_(consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleSrc"))),
  phoPhoIsoToken_(consumes <edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("photonIsolation"))),
  phoNeuIsoToken_(consumes <edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("neutralHadronIsolation"))),
  phoChaIsoToken_(consumes <edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("chargedIsolation"))),

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
  
  tree->Branch("n_pho", &n_pho, "n_pho/I");
  tree->Branch("n_hcalhit", &n_hcalhit, "n_hcalhit/I");
  tree->Branch("pu_true", &pu_true, "pu_true/F");
  tree->Branch("pu_obs", &pu_obs, "pu_obs/I");
  tree->Branch("rho", &rho, "rho/F");
  
  tree->Branch("gen_weight", &gen_weight, "gen_weight/F");
  //}
  

  
  tree->Branch("pho_genpho", &pho_genpho);
  tree->Branch("pho_dR_reco_genpho", &pho_dR_reco_genpho);
  tree->Branch("pho_rpt_reco_genpho", &pho_rpt_reco_genpho);
  
  tree->Branch("pho_sc_energy", &pho_sc_energy);
  tree->Branch("pho_sc_raw_energy", &pho_sc_raw_energy);
  tree->Branch("pho_seed_energy", &pho_seed_energy);
  tree->Branch("pho_seed_corr_energy", &pho_seed_corr_energy);
  tree->Branch("pho_cmssw_hoe", &pho_cmssw_hoe);
  tree->Branch("pho_sc_eta", &pho_sc_eta);
  tree->Branch("pho_sc_phi", &pho_sc_phi);
  tree->Branch("pho_pt", &pho_pt);
  tree->Branch("pho_eta", &pho_pt);
  tree->Branch("pho_phi", &pho_phi);
  tree->Branch("pho_sieie_5x5", &pho_sieie_5x5);
  tree->Branch("pho_r9_5x5", &pho_r9_5x5);

  tree->Branch("pho_pfiso_pho", &pho_pfiso_pho);
  tree->Branch("pho_pfiso_neu", &pho_pfiso_neu);
  tree->Branch("pho_pfiso_cha", &pho_pfiso_cha);
  
  tree->Branch("pho_seed_detid", &pho_seed_detid);
  tree->Branch("pho_seed_subdetid", &pho_seed_subdetid);
  tree->Branch("pho_seed_ieta", &pho_seed_ieta);
  tree->Branch("pho_seed_iphi", &pho_seed_iphi);
  tree->Branch("pho_seed_eta", &pho_seed_eta);
  tree->Branch("pho_seed_phi", &pho_seed_phi);
  tree->Branch("pho_seed_raw_id", &pho_seed_raw_id);
  tree->Branch("pho_seed_hcal_ieta", &pho_seed_hcal_ieta);
  tree->Branch("pho_seed_hcal_iphi", &pho_seed_hcal_iphi);
  
  tree->Branch("pho_seedBC_eta", &pho_seedBC_eta);
  tree->Branch("pho_seedBC_phi", &pho_seedBC_phi);

  tree->Branch("pho_nBCs", &pho_nBCs);
  
  tree->Branch("hcalhit_ieta", &hcalhit_ieta);
  tree->Branch("hcalhit_iphi", &hcalhit_iphi);
  tree->Branch("hcalhit_energy", &hcalhit_energy);
  tree->Branch("hcalhit_seed_dieta", &hcalhit_seed_dieta);
  tree->Branch("hcalhit_seed_diphi", &hcalhit_seed_diphi);
  tree->Branch("hcalhit_raw_id", &hcalhit_raw_id);
  tree->Branch("hcalhit_depth", &hcalhit_depth);
  tree->Branch("hcalhit_pho_index", &hcalhit_pho_index);
  tree->Branch("hcalhit_eta", &hcalhit_eta);
  tree->Branch("hcalhit_phi", &hcalhit_phi);


  tree->Branch("pho_IDLoose", &pho_IDLoose);
  tree->Branch("pho_IDMedium", &pho_IDMedium);
  tree->Branch("pho_IDTight", &pho_IDTight);
    

}


FlatPhoHoEAnalyzer::~FlatPhoHoEAnalyzer()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called for each event  ------------
void
FlatPhoHoEAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  n_pho = 0;
  
  pho_genpho.clear();
  pho_rpt_reco_genpho.clear();
  pho_dR_reco_genpho.clear();

  pho_sc_energy.clear();
  pho_sc_raw_energy.clear();
  pho_seed_energy.clear();
  pho_seed_corr_energy.clear();
  pho_cmssw_hoe.clear();
  pho_sc_eta.clear();
  pho_sc_phi.clear();
  pho_pt.clear();
  pho_eta.clear();
  pho_phi.clear();
  pho_sieie_5x5.clear();
  pho_r9_5x5.clear();

  pho_pfiso_pho.clear();
  pho_pfiso_neu.clear();
  pho_pfiso_cha.clear();
  
  pho_seed_detid.clear();
  pho_seed_subdetid.clear();
  pho_seed_ieta.clear();
  pho_seed_iphi.clear();
  pho_seed_eta.clear();
  pho_seed_phi.clear();

  pho_seedBC_eta.clear();
  pho_seedBC_phi.clear();

  pho_nBCs.clear();
  
  pho_seed_raw_id.clear();
  pho_seed_hcal_ieta.clear();
  pho_seed_hcal_iphi.clear();

  pho_IDLoose.clear();
  pho_IDMedium.clear();
  pho_IDTight.clear();
  
  n_hcalhit = 0;
  hcalhit_ieta.clear();
  hcalhit_iphi.clear();
  hcalhit_energy.clear();
  hcalhit_seed_dieta.clear();
  hcalhit_seed_diphi.clear();
  hcalhit_raw_id.clear();
  hcalhit_depth.clear();
  hcalhit_pho_index.clear();
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

  edm::Handle<edm::ValueMap<float> > phoPhoIsoHandle;
  edm::Handle<edm::ValueMap<float> > phoNeuIsoHandle;
  edm::Handle<edm::ValueMap<float> > phoChaIsoHandle;
  iEvent.getByToken(phoPhoIsoToken_, phoPhoIsoHandle);
  iEvent.getByToken(phoNeuIsoToken_, phoNeuIsoHandle);
  iEvent.getByToken(phoChaIsoToken_, phoChaIsoHandle);

  edm::Handle<std::vector<reco::GenParticle> > genParticlesHandle;
  iEvent.getByToken(genParticlesCollection_, genParticlesHandle);

  //if (iEvent.get(phoToken_).size() > pho_golden.capacity())
  //reallocate_setaddress(iEvent.get(phoToken_).size(), 0);

  //std::cout<<"Outside phoctron loop "<<std::endl;
  
  

  for (const auto& pho : iEvent.get(phoToken_)) {

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
            this_dr2 = reco::deltaR2(pho,*p);
          else 
            continue;

          if (std::abs(p->pdgId()) == 11 and this_dr2 < pho_min_dr2) {
            pho_min_dr2 = this_dr2;
            pho_ptR = pho.pt() / p->pt();
	  }
	  
          if (std::abs(p->pdgId()) == 22 and this_dr2 < pho_min_dr2) {
            pho_min_dr2 = this_dr2;
            pho_ptR = pho.pt() / p->pt();
	  }
	}
      }
    }
  
    // these cuts were decided looking at min_dr and ptR distributions.
    if ( (pho_min_dr2 < 0.0016) and (pho_ptR > 0.7) && (pho_ptR < 1.3) ) 
      genpho = 1;

    pho_dR_reco_genpho.emplace_back( std::sqrt(pho_min_dr2) );
    pho_rpt_reco_genpho.emplace_back(pho_ptR);
    pho_genpho.emplace_back(genpho);

    pho_sc_eta.emplace_back(pho.superCluster()->eta());
    pho_sc_phi.emplace_back(pho.superCluster()->phi());
    pho_pt.emplace_back(pho.pt());
    pho_eta.emplace_back(pho.eta());
    pho_phi.emplace_back(pho.phi());
    pho_sieie_5x5.emplace_back(pho.full5x5_sigmaIetaIeta());
    pho_r9_5x5.emplace_back(pho.full5x5_r9());


    pho_pfiso_pho.emplace_back((*phoPhoIsoHandle)[iEvent.get(phoToken_).ptrAt(n_pho)]);
    pho_pfiso_neu.emplace_back((*phoNeuIsoHandle)[iEvent.get(phoToken_).ptrAt(n_pho)]);
    pho_pfiso_cha.emplace_back((*phoChaIsoHandle)[iEvent.get(phoToken_).ptrAt(n_pho)]);
  

    bool isPassLoose = pho.photonID("cutBasedPhotonID-Fall17-94X-V2-loose");
    bool isPassMedium = pho.photonID("cutBasedPhotonID-Fall17-94X-V2-medium");
    bool isPassTight = pho.photonID("cutBasedPhotonID-Fall17-94X-V2-tight");
    

    pho_IDLoose.push_back((int)isPassLoose);
    pho_IDMedium.push_back((int)isPassMedium);
    pho_IDTight.push_back((int)isPassTight);
    

    EcalClusterLazyTools lazyTool(iEvent, iSetup, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
    
    const reco::SuperCluster& superClus = *pho.superCluster();
    const reco::CaloCluster &seedCluster = *superClus.seed();
    DetId seedId = seedCluster.seed() ;
    pho_seed_detid.emplace_back(seedId.det());
    pho_seed_subdetid.emplace_back(seedId.subdetId());

    float var_pho_seed_eta = -999999.f;
    float var_pho_seed_phi = -999999.f;

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
	    var_pho_seed_eta=ecalrechitPoint.eta();
	    var_pho_seed_phi=ecalrechitPoint.phi();
	  }
	}
      }
    }  
    pho_seed_eta.emplace_back(var_pho_seed_eta);
    pho_seed_phi.emplace_back(var_pho_seed_phi);
    

    std::vector<float> tmp_seedbc_eta;
    std::vector<float> tmp_seedbc_phi;
    int tmp_nBCs = 0;
    
    ////add all the basic clusters seed to an array
    for (reco::CaloCluster_iterator cIt = superClus.clustersBegin(); cIt != superClus.clustersEnd(); ++cIt) {
      const reco::CaloClusterPtr cc = *cIt;
      
      DetId seedId = cc->seed() ;

      //DetId seed = (seedID.hitsAndFractions())[0].first;
      bool isBarrel = seedId.subdetId() == EcalBarrel;
      double var_pho_seed_eta=-999, var_pho_seed_phi=-999;
      bool foundHit = GetSeedEtaPhi(isBarrel,seedId, lazyTool,var_pho_seed_eta, var_pho_seed_phi);
      
      if(foundHit){
	tmp_seedbc_eta.push_back(var_pho_seed_eta);
	tmp_seedbc_phi.push_back(var_pho_seed_phi);
      }

      else{
	tmp_seedbc_eta.push_back(-999);
	tmp_seedbc_phi.push_back(-999);
	
      }
      
      tmp_nBCs++;
      //pho_seed_detid.emplace_back(seedId.det());
      //pho_seed_subdetid.emplace_back(seedId.subdetId());
      

      //std::cout<<"iph : BC eta : BC phi : "<<n_pho<<" "<<var_pho_seed_eta<<" "<<var_pho_seed_phi<<std::endl;
    }  // loop on BCs
    
    pho_nBCs.push_back(tmp_nBCs);

    pho_seedBC_eta.push_back(tmp_seedbc_eta);
    pho_seedBC_phi.push_back(tmp_seedbc_phi);

    int var_pho_seed_ieta = -999999;
    int var_pho_seed_iphi = -999999;
    int var_pho_seed_raw_id = -999999;
    
    int var_pho_seed_hcal_ieta = -999999;
    int var_pho_seed_hcal_iphi = -999999;

    if ( seedId.det() == DetId::Ecal ) {
      if (seedId.subdetId() == EcalBarrel) {
	EBDetId ebId(seedId);
	var_pho_seed_ieta = ebId.ieta();
	var_pho_seed_iphi = ebId.iphi();
	var_pho_seed_raw_id = ebId.rawId();       
      }
      else if (seedId.subdetId() == EcalEndcap) {
	EEDetId eeId(seedId);
	var_pho_seed_ieta = eeId.ix();
        var_pho_seed_iphi = eeId.iy();
        var_pho_seed_raw_id = eeId.rawId();
      }

      // get hold of the seed hcal behind pho seed
      CaloTowerDetId towerId(towerMap_->towerOf(seedId));       
      var_pho_seed_hcal_ieta = towerId.ieta();
      var_pho_seed_hcal_iphi = towerId.iphi();
    }
  
    
    pho_sc_energy.emplace_back(superClus.energy());
    pho_sc_raw_energy.emplace_back(superClus.rawEnergy());
    pho_seed_energy.emplace_back(seedCluster.energy());
    pho_seed_corr_energy.emplace_back(seedCluster.correctedEnergy());
    pho_cmssw_hoe.emplace_back(pho.hadTowOverEm());
   
    pho_seed_ieta.emplace_back(var_pho_seed_ieta);
    pho_seed_iphi.emplace_back(var_pho_seed_iphi);
    pho_seed_raw_id.emplace_back(var_pho_seed_raw_id);

    pho_seed_hcal_ieta.emplace_back(var_pho_seed_hcal_ieta);
    pho_seed_hcal_iphi.emplace_back(var_pho_seed_hcal_iphi);

    
    ++n_pho;
  }
  
  // given the context, should be ok...
  if (n_pho == 0)
    return;

  // just in case

  assert(((void) "ERROR: pho_genpho size doesn't match n_pho!!!", int(pho_genpho.size()) == n_pho));
  assert(((void) "ERROR: pho_rpt_reco_genpho size doesn't match n_pho!!!", int(pho_rpt_reco_genpho.size()) == n_pho));
  assert(((void) "ERROR: pho_dR_reco_genpho size doesn't match n_pho!!!", int(pho_dR_reco_genpho.size()) == n_pho));

  assert(((void) "ERROR: pho_sc_energy size doesn't match n_pho!!!", int(pho_sc_energy.size()) == n_pho));
  assert(((void) "ERROR: pho_sc_raw_energy size doesn't match n_pho!!!", int(pho_sc_raw_energy.size()) == n_pho));
  assert(((void) "ERROR: pho_seed_energy size doesn't match n_pho!!!", int(pho_seed_energy.size()) == n_pho));
  assert(((void) "ERROR: pho_seed_corr_energy size doesn't match n_pho!!!", int(pho_seed_corr_energy.size()) == n_pho));
  assert(((void) "ERROR: pho_cmssw_hoe size doesn't match n_pho!!!", int(pho_cmssw_hoe.size()) == n_pho));
  assert(((void) "ERROR: pho_sc_eta size doesn't match n_pho!!!", int(pho_sc_eta.size()) == n_pho));
  assert(((void) "ERROR: pho_sc_phi size doesn't match n_pho!!!", int(pho_sc_phi.size()) == n_pho));
  assert(((void) "ERROR: pho_pt size doesn't match n_pho!!!", int(pho_pt.size()) == n_pho));
  assert(((void) "ERROR: pho_eta size doesn't match n_pho!!!", int(pho_eta.size()) == n_pho));
  assert(((void) "ERROR: pho_phi size doesn't match n_pho!!!", int(pho_phi.size()) == n_pho));
  assert(((void) "ERROR: pho_sieie_5x5 size doesn't match n_pho!!!", int(pho_sieie_5x5.size()) == n_pho));
  assert(((void) "ERROR: pho_r9_5x5 size doesn't match n_pho!!!", int(pho_r9_5x5.size()) == n_pho));

  assert(((void) "ERROR: pho_pfiso_pho size doesn't match n_pho!!!", int(pho_pfiso_pho.size()) == n_pho));
  assert(((void) "ERROR: pho_pfiso_neu size doesn't match n_pho!!!", int(pho_pfiso_neu.size()) == n_pho));
  assert(((void) "ERROR: pho_pfiso_cha size doesn't match n_pho!!!", int(pho_pfiso_cha.size()) == n_pho));
  
  assert(((void) "ERROR: pho_seed_detid size doesn't match n_pho!!!", int(pho_seed_detid.size()) == n_pho));
  assert(((void) "ERROR: pho_seed_subdetid size doesn't match n_pho!!!", int(pho_seed_subdetid.size()) == n_pho));
  assert(((void) "ERROR: pho_seed_ieta size doesn't match n_pho!!!", int(pho_seed_ieta.size()) == n_pho));
  assert(((void) "ERROR: pho_seed_iphi size doesn't match n_pho!!!", int(pho_seed_iphi.size()) == n_pho));
  assert(((void) "ERROR: pho_seed_eta size doesn't match n_pho!!!", int(pho_seed_eta.size()) == n_pho));
  assert(((void) "ERROR: pho_seed_phi size doesn't match n_pho!!!", int(pho_seed_phi.size()) == n_pho));
  assert(((void) "ERROR: pho_seed_raw_id size doesn't match n_pho!!!", int(pho_seed_raw_id.size()) == n_pho));
  assert(((void) "ERROR: pho_seed_hcal_ieta size doesn't match n_pho!!!", int(pho_seed_hcal_ieta.size()) == n_pho));
  assert(((void) "ERROR: pho_seed_hcal_iphi size doesn't match n_pho!!!", int(pho_seed_hcal_iphi.size()) == n_pho));

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

    for (int iE = 0; iE < n_pho; ++iE) {
      int dieta = calDIEta(pho_seed_hcal_ieta[iE], hcalhit_ieta.back());
      int diphi = calDIPhi(pho_seed_hcal_iphi[iE], hcalhit_iphi.back());
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
    hcalhit_pho_index.emplace_back(imin);

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
  assert(((void) "ERROR: hcalhit_pho_index size doesn't match n_hcalhit!!!", int(hcalhit_pho_index.size()) == n_hcalhit));
  assert(((void) "ERROR: hcalhit_eta size doesn't match n_hcalhit!!!", int(hcalhit_eta.size()) == n_hcalhit));
  assert(((void) "ERROR: hcalhit_phi size doesn't match n_hcalhit!!!", int(hcalhit_phi.size()) == n_hcalhit));

  tree->Fill();

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}

// doing some blatant copy paste from RecoEgamma/EgammaIsolationAlgos/src/EGHcalRecHitSphoctor.cc
int FlatPhoHoEAnalyzer::calDIPhi(int iPhi1, int iPhi2) {
  int dPhi = iPhi1 - iPhi2;
  if (dPhi > 72 / 2)
    dPhi -= 72;
  else if (dPhi < -72 / 2)
    dPhi += 72;
  return dPhi;
}

int FlatPhoHoEAnalyzer::calDIEta(int iEta1, int iEta2) {
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
float FlatPhoHoEAnalyzer::getMinEnergyHCAL(HcalDetId id) const {
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

void FlatPhoHoEAnalyzer::reallocate_setaddress(int n_pho_, int n_hcalhit_)
{


}


// ------------ method called once each job just before starting event loop  ------------
void
FlatPhoHoEAnalyzer::beginJob()
{
  //tree->SetAutoSave(0);

  //tree->SetImplicitMT(false);
  //reallocate_setaddress();
}

// ------------ method called once each job just after ending the event loop  ------------
void
FlatPhoHoEAnalyzer::endJob()
{
  file->cd();
  tree->Write();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FlatPhoHoEAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FlatPhoHoEAnalyzer);


bool FlatPhoHoEAnalyzer::GetSeedEtaPhi(bool isBarrel, DetId seed, EcalClusterLazyTools lazyTool, double &var_pho_seed_eta, double &var_pho_seed_phi){
  
  bool found = false;
  const EcalRecHitCollection * rechits = isBarrel ? lazyTool.getEcalEBRecHitCollection() : lazyTool.getEcalEERecHitCollection();
  EcalRecHitCollection::const_iterator theSeedHit = rechits->find(seed);
  if (theSeedHit != rechits->end()) {
    if ( (theSeedHit->id().rawId() != 0 ) ) {
      if (theCaloGeometry.product() != nullptr) {
	const CaloSubdetectorGeometry *ecalgeo = theCaloGeometry->getSubdetectorGeometry(theSeedHit->id());
	
	if(ecalgeo->getGeometry(theSeedHit->id()) !=nullptr){
	  const GlobalPoint & ecalrechitPoint = (theCaloGeometry.product())->getPosition(theSeedHit->id());
	  var_pho_seed_eta=ecalrechitPoint.eta();
	  var_pho_seed_phi=ecalrechitPoint.phi();
	  found = true;
	}
      }
    }
  }  
  
  return found;
}
