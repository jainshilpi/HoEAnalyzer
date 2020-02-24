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

  int n_ele;
  std::vector<int>  ele_golden;
  std::vector<int>  ele_unknown;
  std::vector<int>  ele_bigbrem;
  std::vector<int>  ele_badtrack;
  std::vector<int>  ele_showering;
  std::vector<int>  ele_gap;
  std::vector<float> ele_track_fbrem;
  std::vector<float> ele_sc_fbrem;
  std::vector<int> ele_nbrem;
  std::vector<int>  ele_genmatch;
  std::vector<float>  ele_sc_energy;
  std::vector<float>  ele_dR_reco_gen;
  std::vector<float>  ele_pt_ratio_reco_gen;
  std::vector<float>  ele_sc_raw_energy;
  std::vector<float>  ele_ecal_energy;
  std::vector<float>  ele_seed_energy;
  std::vector<float>  ele_seed_corr_energy;
  std::vector<float>  ele_cmssw_ele_hoe;
  std::vector<float>  ele_cmssw_ele_hoe_5x5;
  std::vector<float>  ele_sc_eta;
  std::vector<float>  ele_pt;
  std::vector<float>  ele_phi;
  std::vector<float>  ele_sieie_5x5;
  std::vector<float>  ele_pfiso_pho;
  std::vector<float>  ele_pfiso_neu;
  std::vector<float>  ele_pfiso_cha;
  std::vector<int>    ele_seed_detid;
  std::vector<int>    ele_seed_subdetid;
  std::vector<int>    ele_seed_ieta;
  std::vector<int>    ele_seed_iphi;
  std::vector<float>  ele_seed_eta;
  std::vector<float>  ele_seed_phi;
  std::vector<int>    ele_seed_raw_id;
  std::vector<int>    ele_seed_hcal_ieta;
  std::vector<int>    ele_seed_hcal_iphi;

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

  int imin_dieta, min_dieta;
  int imin_diphi, min_diphi;

  float pu_true;
  int pu_obs;
  float rho;

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  static int calDIEta(int iEta1, int iEta2);
  static int calDIPhi(int iPhi1, int iPhi2);
  void reallocate_setaddress(int n_ele = 0, int n_hcalhit = 0);
  float getMinEnergyHCAL(HcalDetId id) const;
  
  int maxDIEta_ = 5;
  int maxDIPhi_ = 5;

  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::View<reco::GsfElectron> > eleToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> >     puCollection_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<HBHERecHitCollection> hbhe_rechits_;
  edm::EDGetTokenT<EcalRecHitCollection> ebReducedRecHitCollection_;
  edm::EDGetTokenT<EcalRecHitCollection> eeReducedRecHitCollection_;
  edm::EDGetTokenT<EcalRecHitCollection> esReducedRecHitCollection_;
  edm::ESHandle<CaloGeometry> theCaloGeometry;  
  edm::ESHandle<CaloTowerConstituentsMap> towerMap_;
  edm::EDGetTokenT<std::vector<reco::GenParticle> >     genParticlesCollection_;

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
FlatEleHoEAnalyzer::FlatEleHoEAnalyzer(const edm::ParameterSet& iConfig)
  :
  eleToken_(consumes<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("electrons"))),
  puCollection_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupCollection"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoSrc"))),
  hbhe_rechits_(consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("hbheInput"))),
  ebReducedRecHitCollection_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ebReducedRecHitCollection"))),
  eeReducedRecHitCollection_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("eeReducedRecHitCollection"))),
  esReducedRecHitCollection_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("esReducedRecHitCollection"))),
  genParticlesCollection_(consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleSrc"))),
  output(iConfig.getParameter<std::string>("output_file")),
  Run2_2018(iConfig.getParameter<bool>("Run2_2018_"))
{
  //now do what ever initialization is needed
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
  ele_golden.clear();
  ele_unknown.clear();
  ele_badtrack.clear();
  ele_bigbrem.clear();
  ele_showering.clear();
  ele_gap.clear();

  ele_track_fbrem.clear();
  ele_sc_fbrem.clear();
  ele_nbrem.clear();

  ele_genmatch.clear();
  ele_pt_ratio_reco_gen.clear();
  ele_dR_reco_gen.clear();
  ele_sc_energy.clear();
  ele_sc_raw_energy.clear();
  ele_ecal_energy.clear();
  ele_seed_energy.clear();
  ele_seed_corr_energy.clear();
  ele_cmssw_ele_hoe.clear();
  ele_cmssw_ele_hoe_5x5.clear();
  ele_sc_eta.clear();
  ele_pt.clear();
  ele_phi.clear();
  ele_sieie_5x5.clear();
  ele_pfiso_pho.clear();
  ele_pfiso_neu.clear();
  ele_pfiso_cha.clear();
  ele_seed_detid.clear();
  ele_seed_subdetid.clear();
  ele_seed_ieta.clear();
  ele_seed_iphi.clear();
  ele_seed_eta.clear();
  ele_seed_phi.clear();
  ele_seed_raw_id.clear();
  ele_seed_hcal_ieta.clear();
  ele_seed_hcal_iphi.clear();

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

  pu_true = -999.f;
  pu_obs = -999;
  rho = -999.f;

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
    rho = -999.f;

  edm::Handle<HBHERecHitCollection> hbheRechitsHandle;
  iEvent.getByToken(hbhe_rechits_, hbheRechitsHandle);
  iSetup.get<CaloGeometryRecord>().get(theCaloGeometry);
  iSetup.get<CaloGeometryRecord>().get(towerMap_);

  edm::Handle<std::vector<reco::GenParticle> > genParticlesHandle;
  iEvent.getByToken(genParticlesCollection_, genParticlesHandle);

  for (const auto& ele : iEvent.get(eleToken_)) {
    int genmatch = 0;
    double min_dr = 9999.9;
    double ptR = 9999.9;
     
    if (genParticlesHandle.isValid()) {
      for (std::vector<reco::GenParticle>::const_iterator ip = genParticlesHandle->begin(); ip != genParticlesHandle->end(); ++ip) {
	const reco::Candidate *p = (const reco::Candidate*)&(*ip);
	//std::cout << " p->pdgId() " << p->pdgId() << std::endl;
	if ( (std::abs(p->pdgId())) ==11 ) { 
	  //std::cout << "-----  p->status() " << p->status() << " p->pdgId() " << p->pdgId() << std::endl;
	}
	if ( (p->status()==1) &&  ((std::abs(p->pdgId())) == 11) ) {
	  //std::cout << "checking if the reco ele match with this one" << std::endl;
	  double this_dr=reco::deltaR(ele,*p);
	  //std::cout << "this_dr " << this_dr << std::endl;
	  if (this_dr<min_dr) {
	    min_dr=this_dr;
	    ptR=ele.pt()/p->pt();
	  }
	}  
      }
    }
  
    if ( (min_dr<0.04) && (ptR>0.7) && (ptR<1.3) )  genmatch=1; // these cuts were decided looking at min_dr and ptR distributions.
    ele_dR_reco_gen.emplace_back(min_dr);
    ele_pt_ratio_reco_gen.emplace_back(ptR);
    //  std::cout << "genmatch = " << genmatch <<  " min_dr " << min_dr << " ptR " << ptR   <<  std::endl;    
    ele_genmatch.emplace_back(genmatch);

    ele_sc_eta.emplace_back(ele.superCluster()->eta());
    ele_pt.emplace_back(ele.pt());
    ele_phi.emplace_back(ele.phi());
    ele_sieie_5x5.emplace_back(ele.full5x5_sigmaIetaIeta());
  
    reco::GsfElectron::PflowIsolationVariables pfIso = ele.pfIsolationVariables();
    ele_pfiso_pho.emplace_back(pfIso.sumPhotonEt);
    ele_pfiso_neu.emplace_back(pfIso.sumNeutralHadronEt);
    ele_pfiso_cha.emplace_back(pfIso.sumChargedHadronPt);

    EcalClusterLazyTools lazyTool(iEvent, iSetup, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);

    const reco::SuperCluster& superClus = *ele.superCluster();
    const reco::CaloCluster &seedCluster = *superClus.seed();
    DetId seedId = seedCluster.seed() ;
    ele_seed_detid.emplace_back(seedId.det());
    ele_seed_subdetid.emplace_back(seedId.subdetId());

    float var_ele_seed_eta = -99.f;
    float var_ele_seed_phi = -99.f;

    DetId seed = (seedCluster.hitsAndFractions())[0].first;
    bool isBarrel = seed.subdetId() == EcalBarrel;
    const EcalRecHitCollection * rechits = (isBarrel?lazyTool.getEcalEBRecHitCollection():lazyTool.getEcalEERecHitCollection());
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

    int var_ele_seed_ieta = -999;
    int var_ele_seed_iphi = -999;
    int var_ele_seed_raw_id = -999;
    
    int var_ele_seed_hcal_ieta = -999;
    int var_ele_seed_hcal_iphi = -999;
    
    if ( seedId.det() == DetId::Ecal ) {
      if (seedId.subdetId() == EcalBarrel) {
	EBDetId ebId(seedId);
	var_ele_seed_ieta=ebId.ieta();
	var_ele_seed_iphi=ebId.iphi();
	var_ele_seed_raw_id=ebId.rawId();       
      }
      else if (seedId.subdetId() == EcalEndcap) {
	EEDetId eeId(seedId);
	var_ele_seed_ieta=eeId.ix();
        var_ele_seed_iphi=eeId.iy();
        var_ele_seed_raw_id=eeId.rawId();
      }

      // get hold of the seed hcal behind ele seed
      CaloTowerDetId towerId(towerMap_->towerOf(seedId));       
      var_ele_seed_hcal_ieta=towerId.ieta();
      var_ele_seed_hcal_iphi=towerId.iphi();
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
      var_golden=1;

    if (ele.classification() == reco::GsfElectron::UNKNOWN)
      var_unknown=1;

    if (ele.classification() == reco::GsfElectron::BIGBREM)
      var_bigbrem=1;

    if (ele.classification() == reco::GsfElectron::BADTRACK)
      var_badtrack=1;

    if (ele.classification() == reco::GsfElectron::SHOWERING)
      var_showering=1;

    if (ele.classification() == reco::GsfElectron::GAP)
      var_gap=1;

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
    ele_cmssw_ele_hoe.emplace_back(ele.hcalOverEcal());
    ele_cmssw_ele_hoe_5x5.emplace_back(ele.full5x5_hcalOverEcal());
   
    ele_seed_ieta.emplace_back(var_ele_seed_ieta);
    ele_seed_iphi.emplace_back(var_ele_seed_iphi);
    ele_seed_raw_id.emplace_back(var_ele_seed_raw_id);

    ele_seed_hcal_ieta.emplace_back(var_ele_seed_hcal_ieta);
    ele_seed_hcal_iphi.emplace_back(var_ele_seed_hcal_iphi);

    ++n_ele;
  }

  // just in case
  assert(((void) "ERROR: ele_golden size doesn't match n_ele!!!", int(ele_golden.size()) == n_ele));
  assert(((void) "ERROR: ele_unknown size doesn't match n_ele!!!", int(ele_unknown.size()) == n_ele));
  assert(((void) "ERROR: ele_badtrack size doesn't match n_ele!!!", int(ele_badtrack.size()) == n_ele));
  assert(((void) "ERROR: ele_bigbrem size doesn't match n_ele!!!", int(ele_bigbrem.size()) == n_ele));
  assert(((void) "ERROR: ele_showering size doesn't match n_ele!!!", int(ele_showering.size()) == n_ele));
  assert(((void) "ERROR: ele_gap size doesn't match n_ele!!!", int(ele_gap.size()) == n_ele));

  assert(((void) "ERROR: ele_track_fbrem size doesn't match n_ele!!!", int(ele_track_fbrem.size()) == n_ele));
  assert(((void) "ERROR: ele_sc_fbrem size doesn't match n_ele!!!", int(ele_sc_fbrem.size()) == n_ele));
  assert(((void) "ERROR: ele_nbrem size doesn't match n_ele!!!", int(ele_nbrem.size()) == n_ele));

  assert(((void) "ERROR: ele_genmatch size doesn't match n_ele!!!", int(ele_genmatch.size()) == n_ele));
  assert(((void) "ERROR: ele_pt_ratio_reco_gen size doesn't match n_ele!!!", int(ele_pt_ratio_reco_gen.size()) == n_ele));
  assert(((void) "ERROR: ele_dR_reco_gen size doesn't match n_ele!!!", int(ele_dR_reco_gen.size()) == n_ele));
  assert(((void) "ERROR: ele_sc_energy size doesn't match n_ele!!!", int(ele_sc_energy.size()) == n_ele));
  assert(((void) "ERROR: ele_sc_raw_energy size doesn't match n_ele!!!", int(ele_sc_raw_energy.size()) == n_ele));
  assert(((void) "ERROR: ele_ecal_energy size doesn't match n_ele!!!", int(ele_ecal_energy.size()) == n_ele));
  assert(((void) "ERROR: ele_seed_energy size doesn't match n_ele!!!", int(ele_seed_energy.size()) == n_ele));
  assert(((void) "ERROR: ele_seed_corr_energy size doesn't match n_ele!!!", int(ele_seed_corr_energy.size()) == n_ele));
  assert(((void) "ERROR: ele_cmssw_ele_hoe size doesn't match n_ele!!!", int(ele_cmssw_ele_hoe.size()) == n_ele));
  assert(((void) "ERROR: ele_cmssw_ele_hoe_5x5 size doesn't match n_ele!!!", int(ele_cmssw_ele_hoe_5x5.size()) == n_ele));
  assert(((void) "ERROR: ele_sc_eta size doesn't match n_ele!!!", int(ele_sc_eta.size()) == n_ele));
  assert(((void) "ERROR: ele_pt size doesn't match n_ele!!!", int(ele_pt.size()) == n_ele));
  assert(((void) "ERROR: ele_phi size doesn't match n_ele!!!", int(ele_phi.size()) == n_ele));
  assert(((void) "ERROR: ele_sieie_5x5 size doesn't match n_ele!!!", int(ele_sieie_5x5.size()) == n_ele));
  assert(((void) "ERROR: ele_pfiso_pho size doesn't match n_ele!!!", int(ele_pfiso_pho.size()) == n_ele));
  assert(((void) "ERROR: ele_pfiso_neu size doesn't match n_ele!!!", int(ele_pfiso_neu.size()) == n_ele));
  assert(((void) "ERROR: ele_pfiso_cha size doesn't match n_ele!!!", int(ele_pfiso_cha.size()) == n_ele));
  assert(((void) "ERROR: ele_seed_detid size doesn't match n_ele!!!", int(ele_seed_detid.size()) == n_ele));
  assert(((void) "ERROR: ele_seed_subdetid size doesn't match n_ele!!!", int(ele_seed_subdetid.size()) == n_ele));
  assert(((void) "ERROR: ele_seed_ieta size doesn't match n_ele!!!", int(ele_seed_ieta.size()) == n_ele));
  assert(((void) "ERROR: ele_seed_iphi size doesn't match n_ele!!!", int(ele_seed_iphi.size()) == n_ele));
  assert(((void) "ERROR: ele_seed_eta size doesn't match n_ele!!!", int(ele_seed_eta.size()) == n_ele));
  assert(((void) "ERROR: ele_seed_phi size doesn't match n_ele!!!", int(ele_seed_phi.size()) == n_ele));
  assert(((void) "ERROR: ele_seed_raw_id size doesn't match n_ele!!!", int(ele_seed_raw_id.size()) == n_ele));
  assert(((void) "ERROR: ele_seed_hcal_ieta size doesn't match n_ele!!!", int(ele_seed_hcal_ieta.size()) == n_ele));
  assert(((void) "ERROR: ele_seed_hcal_iphi size doesn't match n_ele!!!", int(ele_seed_hcal_iphi.size()) == n_ele));

  for (auto& hcalrh : iEvent.get(hbhe_rechits_) ) {
    if (hcalrh.energy() < getMinEnergyHCAL(hcalrh.id()))
      continue;

    hcalhit_ieta.emplace_back(hcalrh.id().ieta());
    hcalhit_iphi.emplace_back(hcalrh.id().iphi());
    hcalhit_energy.emplace_back(hcalrh.energy());

    hcalhit_raw_id.emplace_back(hcalrh.id().rawId());
    hcalhit_depth.emplace_back(hcalrh.id().depth());

    float rechitEta = -99.f;
    float rechitPhi = -99.f;
    if (hcalrh.id().rawId() != 0) {
      if (theCaloGeometry.product() != nullptr) {
        const CaloSubdetectorGeometry *geo = theCaloGeometry->getSubdetectorGeometry(hcalrh.id());

        if(geo->getGeometry(hcalrh.id()) != nullptr) {
          const GlobalPoint & rechitPoint = theCaloGeometry.product()->getPosition(hcalrh.id());

          rechitEta=rechitPoint.eta();
          rechitPhi=rechitPoint.phi();	
        }
      }
    }
    hcalhit_eta.emplace_back(rechitEta);
    hcalhit_phi.emplace_back(rechitPhi);

    imin_dieta = -1;
    min_dieta = 999;

    for (int iE = 0; iE < n_ele; ++iE) {
      int dieta = calDIEta(ele_seed_hcal_ieta[iE], hcalhit_ieta.back());

      if (std::abs(dieta) < min_dieta) {
        min_dieta = std::abs(dieta);
        imin_dieta = iE;
      }
    }

    imin_diphi = -1;
    min_diphi = 999;

    for (int iE = 0; iE < n_ele; ++iE) {
      int diphi = calDIEta(ele_seed_hcal_iphi[iE], hcalhit_iphi.back());

      if (std::abs(diphi) < min_diphi) {
        min_diphi = std::abs(diphi);
        imin_diphi = iE;
      }
    }

    if (imin_dieta == imin_diphi) {
      hcalhit_seed_dieta.emplace_back(min_dieta);
      hcalhit_seed_diphi.emplace_back(min_diphi);
      hcalhit_ele_index.emplace_back(imin_dieta);
    }
    else {
      hcalhit_seed_dieta.emplace_back(999);
      hcalhit_seed_diphi.emplace_back(999);
      hcalhit_ele_index.emplace_back(-1);
    }

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

  if (n_ele > int(ele_golden.capacity()))
    reallocate_setaddress(n_ele, 0);
  if (n_hcalhit > int(hcalhit_depth.capacity()))
    reallocate_setaddress(0, n_hcalhit);

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
      return 0.7;
    else if ( (Run2_2018 == 0) ) { // means Run3
      if (id.depth() == 1)
	return 0.1;
      else if (id.depth() == 2)
	return 0.2;
      else
	return 0.3;
    }
    else // neither 2018, nor Run3, not supported
      return 9999.0;
  } 

  else if (id.subdetId() == HcalEndcap) {
    if (id.depth() == 1)
      return 0.1;
    else
      return 0.2;
  } else
    return 9999.0;
}

void FlatEleHoEAnalyzer::reallocate_setaddress(int n_ele, int n_hcalhit)
{
  static int cap_ele = 8;
  cap_ele = (n_ele == 0) ? cap_ele : n_ele;
  ele_golden.reserve(cap_ele);
  ele_unknown.reserve(cap_ele);
  ele_badtrack.reserve(cap_ele);
  ele_bigbrem.reserve(cap_ele);
  ele_showering.reserve(cap_ele);
  ele_gap.reserve(cap_ele);

  ele_track_fbrem.reserve(cap_ele);
  ele_sc_fbrem.reserve(cap_ele);
  ele_nbrem.reserve(cap_ele);

  ele_genmatch.reserve(cap_ele);
  ele_pt_ratio_reco_gen.reserve(cap_ele);
  ele_dR_reco_gen.reserve(cap_ele);
  ele_sc_energy.reserve(cap_ele);
  ele_sc_raw_energy.reserve(cap_ele);
  ele_ecal_energy.reserve(cap_ele);
  ele_seed_energy.reserve(cap_ele);
  ele_seed_corr_energy.reserve(cap_ele);
  ele_cmssw_ele_hoe.reserve(cap_ele);
  ele_cmssw_ele_hoe_5x5.reserve(cap_ele);
  ele_sc_eta.reserve(cap_ele);
  ele_pt.reserve(cap_ele);
  ele_phi.reserve(cap_ele);
  ele_sieie_5x5.reserve(cap_ele);
  ele_pfiso_pho.reserve(cap_ele);
  ele_pfiso_neu.reserve(cap_ele);
  ele_pfiso_cha.reserve(cap_ele);
  ele_seed_detid.reserve(cap_ele);
  ele_seed_subdetid.reserve(cap_ele);
  ele_seed_ieta.reserve(cap_ele);
  ele_seed_iphi.reserve(cap_ele);
  ele_seed_eta.reserve(cap_ele);
  ele_seed_phi.reserve(cap_ele);
  ele_seed_raw_id.reserve(cap_ele);
  ele_seed_hcal_ieta.reserve(cap_ele);
  ele_seed_hcal_iphi.reserve(cap_ele);

  static int cap_hcalhit = 128;
  cap_hcalhit = (n_hcalhit == 0) ? cap_hcalhit : n_hcalhit;
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

  if (n_ele == 0 and n_hcalhit == 0) {
    tree->Branch("n_ele", &n_ele, "n_ele/I");
    tree->Branch("n_hcalhit", &n_hcalhit, "n_hcalhit/I");
    tree->Branch("pu_true", &pu_true, "pu_true/F");
    tree->Branch("pu_obs", &pu_obs, "pu_obs/I");
    tree->Branch("rho", &rho, "rho/F");
  }

  static TBranch *b_ele_golden = tree->Branch("ele_golden", ele_golden.data(), "ele_golden[n_ele]/I");
  static TBranch *b_ele_unknown = tree->Branch("ele_unknown", ele_unknown.data(), "ele_unknown[n_ele]/I");
  static TBranch *b_ele_bigbrem = tree->Branch("ele_bigbrem", ele_bigbrem.data(), "ele_bigbrem[n_ele]/I");
  static TBranch *b_ele_gap = tree->Branch("ele_gap", ele_gap.data(), "ele_gap[n_ele]/I");
  static TBranch *b_ele_badtrack = tree->Branch("ele_badtrack", ele_badtrack.data(), "ele_badtrack[n_ele]/I");
  static TBranch *b_ele_showering = tree->Branch("ele_showering", ele_showering.data(), "ele_showering[n_ele]/I");
  static TBranch *b_ele_track_fbrem = tree->Branch("ele_track_fbrem", ele_track_fbrem.data(), "ele_track_fbrem[n_ele]/F");
  static TBranch *b_ele_sc_fbrem = tree->Branch("ele_sc_fbrem", ele_sc_fbrem.data(), "ele_sc_fbrem[n_ele]/F");
  static TBranch *b_ele_nbrem = tree->Branch("ele_nbrem", ele_nbrem.data(), "ele_nbrem[n_ele]/I");
  static TBranch *b_ele_genmatch = tree->Branch("ele_genmatch", ele_genmatch.data(), "ele_genmatch[n_ele]/I");
  static TBranch *b_ele_dR_reco_gen = tree->Branch("ele_dR_reco_gen", ele_dR_reco_gen.data(), "ele_dR_reco_gen[n_ele]/F");
  static TBranch *b_ele_pt_ratio_reco_gen = tree->Branch("ele_pt_ratio_reco_gen", ele_pt_ratio_reco_gen.data(), "ele_pt_ratio_reco_gen[n_ele]/F");
  static TBranch *b_ele_sc_energy = tree->Branch("ele_sc_energy", ele_sc_energy.data(), "ele_sc_energy[n_ele]/F");
  static TBranch *b_ele_sc_raw_energy = tree->Branch("ele_sc_raw_energy", ele_sc_raw_energy.data(), "ele_sc_raw_energy[n_ele]/F");
  static TBranch *b_ele_ecal_energy = tree->Branch("ele_ecal_energy", ele_ecal_energy.data(), "ele_ecal_energy[n_ele]/F");
  static TBranch *b_ele_seed_energy = tree->Branch("ele_seed_energy", ele_seed_energy.data(), "ele_seed_energy[n_ele]/F");
  static TBranch *b_ele_seed_corr_energy = tree->Branch("ele_seed_corr_energy", ele_seed_corr_energy.data(), "ele_seed_corr_energy[n_ele]/F");
  static TBranch *b_ele_cmssw_ele_hoe = tree->Branch("ele_cmssw_ele_hoe", ele_cmssw_ele_hoe.data(), "cmssw_ele_hoe[n_ele]/F");
  static TBranch *b_ele_cmssw_ele_hoe_5x5 = tree->Branch("ele_cmssw_ele_hoe_5x5", ele_cmssw_ele_hoe_5x5.data(), "ele_cmssw_ele_hoe_5x5[n_ele]/F");
  static TBranch *b_ele_sc_eta = tree->Branch("ele_sc_eta", ele_sc_eta.data(), "ele_sc_eta[n_ele]/F");
  static TBranch *b_ele_pt = tree->Branch("ele_pt", ele_pt.data(), "ele_pt[n_ele]/F");
  static TBranch *b_ele_phi = tree->Branch("ele_phi", ele_phi.data(), "ele_phi[n_ele]/F");
  static TBranch *b_ele_sieie_5x5 = tree->Branch("ele_sieie_5x5", ele_sieie_5x5.data(), "ele_sieie_5x5[n_ele]/F");
  static TBranch *b_ele_pfiso_pho = tree->Branch("ele_pfiso_pho", ele_pfiso_pho.data(), "ele_pfiso_pho[n_ele]/F");
  static TBranch *b_ele_pfiso_neu = tree->Branch("ele_pfiso_neu", ele_pfiso_neu.data(), "ele_pfiso_neu[n_ele]/F");
  static TBranch *b_ele_pfiso_cha = tree->Branch("ele_pfiso_cha", ele_pfiso_cha.data(), "ele_pfiso_cha[n_ele]/F");
  static TBranch *b_ele_seed_detid = tree->Branch("ele_seed_detid", ele_seed_detid.data(), "ele_seed_detid[n_ele]/I");
  static TBranch *b_ele_seed_subdetid = tree->Branch("ele_seed_subdetid", ele_seed_subdetid.data(), "ele_seed_subdetid[n_ele]/I");
  static TBranch *b_ele_seed_ieta = tree->Branch("ele_seed_ieta", ele_seed_ieta.data(), "ele_seed_ieta[n_ele]/I");
  static TBranch *b_ele_seed_iphi = tree->Branch("ele_seed_iphi", ele_seed_iphi.data(), "ele_seed_iphi[n_ele]/I");
  static TBranch *b_ele_seed_eta = tree->Branch("ele_seed_eta", ele_seed_eta.data(), "ele_seed_eta[n_ele]/F");
  static TBranch *b_ele_seed_phi = tree->Branch("ele_seed_phi", ele_seed_phi.data(), "ele_seed_phi[n_ele]/F");
  static TBranch *b_ele_seed_raw_id = tree->Branch("ele_seed_raw_id", ele_seed_raw_id.data(), "ele_seed_raw_id[n_ele]/I");
  static TBranch *b_ele_seed_hcal_ieta = tree->Branch("ele_seed_hcal_ieta", ele_seed_hcal_ieta.data(), "ele_seed_hcal_ieta[n_ele]/I");
  static TBranch *b_ele_seed_hcal_iphi = tree->Branch("ele_seed_hcal_iphi", ele_seed_hcal_iphi.data(), "ele_seed_hcal_iphi[n_ele]/I");

  static TBranch *b_hcalhit_ieta = tree->Branch("hcalhit_ieta", hcalhit_ieta.data(), "hcalhit_ieta[n_hcalhit]/I");
  static TBranch *b_hcalhit_iphi = tree->Branch("hcalhit_iphi", hcalhit_iphi.data(), "hcalhit_iphi[n_hcalhit]/I");
  static TBranch *b_hcalhit_energy = tree->Branch("hcalhit_energy", hcalhit_energy.data(), "hcalhit_energy[n_hcalhit]/F");
  static TBranch *b_hcalhit_seed_dieta = tree->Branch("hcalhit_seed_dieta", hcalhit_seed_dieta.data(), "hcalhit_seed_dieta[n_hcalhit]/I");
  static TBranch *b_hcalhit_seed_diphi = tree->Branch("hcalhit_seed_diphi", hcalhit_seed_diphi.data(), "hcalhit_seed_diphi[n_hcalhit]/I");
  static TBranch *b_hcalhit_raw_id = tree->Branch("hcalhit_raw_id", hcalhit_raw_id.data(), "hcalhit_raw_id[n_hcalhit]/I");
  static TBranch *b_hcalhit_depth = tree->Branch("hcalhit_depth", hcalhit_depth.data(), "hcalhit_depth[n_hcalhit]/I");
  static TBranch *b_hcalhit_ele_index = tree->Branch("hcalhit_ele_index", hcalhit_ele_index.data(), "hcalhit_ele_index[n_hcalhit]/I");
  static TBranch *b_hcalhit_eta = tree->Branch("hcalhit_eta", hcalhit_eta.data(), "hcalhit_eta[n_hcalhit]/F");
  static TBranch *b_hcalhit_phi = tree->Branch("hcalhit_phi", hcalhit_phi.data(), "hcalhit_phi[n_hcalhit]/F");

  if (n_ele != 0) {
    std::cout << "Electron block realloc..." << std::endl;

    b_ele_golden->SetAddress(ele_golden.data());
    b_ele_unknown->SetAddress(ele_unknown.data());
    b_ele_bigbrem->SetAddress(ele_bigbrem.data());
    b_ele_gap->SetAddress(ele_gap.data());
    b_ele_badtrack->SetAddress(ele_badtrack.data());
    b_ele_showering->SetAddress(ele_showering.data());
    b_ele_track_fbrem->SetAddress(ele_track_fbrem.data());
    b_ele_sc_fbrem->SetAddress(ele_sc_fbrem.data());
    b_ele_nbrem->SetAddress(ele_nbrem.data());
    b_ele_genmatch->SetAddress(ele_genmatch.data());
    b_ele_dR_reco_gen->SetAddress(ele_dR_reco_gen.data());
    b_ele_pt_ratio_reco_gen->SetAddress(ele_pt_ratio_reco_gen.data());
    b_ele_sc_energy->SetAddress(ele_sc_energy.data());
    b_ele_sc_raw_energy->SetAddress(ele_sc_raw_energy.data());
    b_ele_ecal_energy->SetAddress(ele_ecal_energy.data());
    b_ele_seed_energy->SetAddress(ele_seed_energy.data());
    b_ele_seed_corr_energy->SetAddress(ele_seed_corr_energy.data());
    b_ele_cmssw_ele_hoe->SetAddress(ele_cmssw_ele_hoe.data());
    b_ele_cmssw_ele_hoe_5x5->SetAddress(ele_cmssw_ele_hoe_5x5.data());
    b_ele_sc_eta->SetAddress(ele_sc_eta.data());
    b_ele_pt->SetAddress(ele_pt.data());
    b_ele_phi->SetAddress(ele_phi.data());
    b_ele_sieie_5x5->SetAddress(ele_sieie_5x5.data());
    b_ele_pfiso_pho->SetAddress(ele_pfiso_pho.data());
    b_ele_pfiso_neu->SetAddress(ele_pfiso_neu.data());
    b_ele_pfiso_cha->SetAddress(ele_pfiso_cha.data());
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

  if (n_hcalhit != 0) {
    std::cout << "Hcalhit block realloc..." << std::endl;

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
}


// ------------ method called once each job just before starting event loop  ------------
void
FlatEleHoEAnalyzer::beginJob()
{
  file = new TFile(output.c_str(), "recreate");
  tree = new TTree("tree", "");
  tree->SetAutoSave(0);
  tree->SetImplicitMT(false);
  reallocate_setaddress();
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
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FlatEleHoEAnalyzer);
