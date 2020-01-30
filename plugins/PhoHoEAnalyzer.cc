// -*- C++ -*-
//
// Package:    HoE/EleHoEAnalyzer
// Class:      EleHoEAnalyzer
//
/**\class EleHoEAnalyzer EleHoEAnalyzer.cc HoE/EleHoEAnalyzer/plugins/EleHoEAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Swagata Mukherjee
//         Created:  Mon, 04 Nov 2019 14:29:29 GMT
//
//


// system include files
#include <memory>
// user include files
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHadTower.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
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
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "TTree.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
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


class PhoHoEAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit PhoHoEAnalyzer(const edm::ParameterSet&);
  ~PhoHoEAnalyzer();
  
  static edm::ParameterSetDescription makePSetDescription();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  edm::Service<TFileService> fs;
  TTree   *tree = fs->make<TTree>("EventTree", "EventData");
  
  std::vector<int>  pho_haspixseed;
  std::vector<int>  pho_ispf;
  std::vector<int>  pho_isstd;

  std::vector<int>  pho_genmatched;
  std::vector<float>  scEn;
  std::vector<float>  dR_recoPho_genPho;
  std::vector<float>  ptRecoPho_by_ptGenPho;
  std::vector<float>  phoSCRawEn;
  std::vector<float>  seedEn;
  std::vector<float>  seedEnCorr;
  std::vector<float>  cmssw_phoHoE;
  std::vector<float>  phoScEta;
  std::vector<float>  phoPt;
  std::vector<float>  phoPhi;
  std::vector<float>  phoSigmaIEtaIEtaFull5x5;
  std::vector<float>  phoR9;
  std::vector<float>  phoR9Full5x5;
  std::vector<float>  phoPFPhoIso;
  std::vector<float>  phoPFNeuIso;
  std::vector<float>  phoPFChIso;
  std::vector<int>    phoSeedDet;
  std::vector<int>    phoSeedSubdet;
  std::vector<int>    phoSeedIeta;
  std::vector<int>    phoSeedIphi;
  std::vector<float>  phoSeedEta;
  std::vector<float>  phoSeedPhi;
  std::vector<int>    phoSeedRawID;
  std::vector<int>    seedHcalIeta;
  std::vector<int>    seedHcalIphi;

  std::vector<std::vector<int>>    hcalRechitIeta;
  std::vector<std::vector<int>>    hcalRechitIphi;
  std::vector<std::vector<float>>  hcalRechitEnergy;
  std::vector<std::vector<int>>    hcalRechitAbsDIetaFromPhoSeed;
  std::vector<std::vector<int>>    hcalRechitAbsDIphiFromPhoSeed;
  std::vector<std::vector<int>>    hcalRechitRawID;
  std::vector<std::vector<int>>    hcalRechitDepth; // mostly for Run 3 //
  std::vector<std::vector<float>>    hcalRechitEta;
  std::vector<std::vector<float>>    hcalRechitPhi;

  std::vector<int>    perPho_hcalRechitIeta;
  std::vector<int>    perPho_hcalRechitIphi;
  std::vector<float>  perPho_hcalRechitEnergy;
  std::vector<int>    perPho_hcalRechitAbsDIetaFromPhoSeed;
  std::vector<int>    perPho_hcalRechitAbsDIphiFromPhoSeed;
  std::vector<int>    perPho_hcalRechitRawID;
  std::vector<int>    perPho_hcalRechitDepth; // mostly for Run 3 //
  std::vector<float>  perPho_hcalRechitEta;
  std::vector<float>  perPho_hcalRechitPhi;

  float  puTrue;
  int puObs;
  float rho;
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  static int calDIEta(int iEta1, int iEta2);
  static int calDIPhi(int iPhi1, int iPhi2);
  float getMinEnergyHCAL_(HcalDetId id) const;
  
  int maxDIEta_=5;
  int maxDIPhi_=5;

  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::View<reco::Photon> > phoToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> >     puCollection_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<HBHERecHitCollection> hbhe_rechits_;
  edm::EDGetTokenT<EcalRecHitCollection> ebReducedRecHitCollection_;
  edm::EDGetTokenT<EcalRecHitCollection> eeReducedRecHitCollection_;
  edm::EDGetTokenT<EcalRecHitCollection> esReducedRecHitCollection_;
  edm::ESHandle<CaloGeometry> theCaloGeometry;  
  edm::ESHandle<CaloTowerConstituentsMap> towerMap_;
  edm::EDGetTokenT<std::vector<reco::GenParticle> >     genParticlesCollection_;
  edm::EDGetTokenT<edm::ValueMap<float> > phoPhoIsoToken_; 
  edm::EDGetTokenT<edm::ValueMap<float> > phoNeuIsoToken_; 
  edm::EDGetTokenT<edm::ValueMap<float> > phoChaIsoToken_; 

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
PhoHoEAnalyzer::PhoHoEAnalyzer(const edm::ParameterSet& iConfig)
  :
  phoToken_(consumes<edm::View<reco::Photon> >(iConfig.getParameter<edm::InputTag>("photons"))),
  puCollection_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupCollection"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoSrc"))),
  hbhe_rechits_(consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("hbheInput"))),
  ebReducedRecHitCollection_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ebReducedRecHitCollection"))),
  eeReducedRecHitCollection_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("eeReducedRecHitCollection"))),
  esReducedRecHitCollection_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("esReducedRecHitCollection"))),
  genParticlesCollection_(consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleSrc"))),
  phoPhoIsoToken_(consumes <edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("photonIsolation"))),
  phoNeuIsoToken_(consumes <edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("neutralHadronIsolation"))),
  phoChaIsoToken_(consumes <edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("chargedIsolation"))),
  Run2_2018(iConfig.getParameter<bool>("Run2_2018_"))
{
  //now do what ever initialization is needed
  tree->Branch("pho_haspixseed_",&pho_haspixseed);
  tree->Branch("pho_ispf_",&pho_ispf);
  tree->Branch("pho_isstd_",&pho_isstd);

  tree->Branch("pho_genmatched_",&pho_genmatched);
  tree->Branch("dR_recoPho_genPho_",&dR_recoPho_genPho);
  tree->Branch("ptRecoPho_by_ptGenPho_",&ptRecoPho_by_ptGenPho);
  tree->Branch("scEn_",&scEn);
  tree->Branch("phoSCRawEn_",&phoSCRawEn);
  tree->Branch("seedEn_",&seedEn);
  tree->Branch("seedEnCorr_",&seedEnCorr);
  tree->Branch("cmssw_phoHoE_",&cmssw_phoHoE);
  tree->Branch("phoScEta_",&phoScEta);
  tree->Branch("phoPt_",&phoPt);
  tree->Branch("phoPhi_",&phoPhi);
  tree->Branch("phoSigmaIEtaIEtaFull5x5_",&phoSigmaIEtaIEtaFull5x5);
  tree->Branch("phoR9_",&phoR9);
  tree->Branch("phoR9Full5x5_",&phoR9Full5x5);
  tree->Branch("phoPFPhoIso_",&phoPFPhoIso);
  tree->Branch("phoPFNeuIso_",&phoPFNeuIso);
  tree->Branch("phoPFChIso_",&phoPFChIso);
  tree->Branch("phoSeedDet_",&phoSeedDet);
  tree->Branch("phoSeedSubdet_",&phoSeedSubdet);
  tree->Branch("phoSeedIeta_",&phoSeedIeta);
  tree->Branch("phoSeedIphi_",&phoSeedIphi);
  tree->Branch("phoSeedEta_",&phoSeedEta);
  tree->Branch("phoSeedPhi_",&phoSeedPhi);
  tree->Branch("phoSeedRawID_",&phoSeedRawID);
  tree->Branch("seedHcalIeta_",&seedHcalIeta);
  tree->Branch("seedHcalIphi_",&seedHcalIphi);
  tree->Branch("hcalRechitIeta_",&hcalRechitIeta);
  tree->Branch("hcalRechitIphi_",&hcalRechitIphi);
  tree->Branch("hcalRechitEnergy_",&hcalRechitEnergy);
  tree->Branch("hcalRechitAbsDIetaFromPhoSeed_",&hcalRechitAbsDIetaFromPhoSeed);
  tree->Branch("hcalRechitAbsDIphiFromPhoSeed_",&hcalRechitAbsDIphiFromPhoSeed);
  tree->Branch("hcalRechitRawID_",&hcalRechitRawID);
  tree->Branch("hcalRechitDepth_",&hcalRechitDepth);
  tree->Branch("hcalRechitEta_",&hcalRechitEta);
  tree->Branch("hcalRechitPhi_",&hcalRechitPhi);
  tree->Branch("puTrue_", &puTrue);
  tree->Branch("puObs_", &puObs);
  tree->Branch("rho_", &rho);
}


PhoHoEAnalyzer::~PhoHoEAnalyzer()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called for each event  ------------
void
PhoHoEAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //   std::cout << " \n ****** new event ... " << std::endl; 
  using namespace edm;

  pho_haspixseed.clear();
  pho_ispf.clear();
  pho_isstd.clear();

  pho_genmatched.clear();
  ptRecoPho_by_ptGenPho.clear();
  dR_recoPho_genPho.clear();
  scEn.clear();
  phoSCRawEn.clear();
  seedEn.clear();
  seedEnCorr.clear();
  cmssw_phoHoE.clear();
  phoScEta.clear();
  phoPt.clear();
  phoPhi.clear();
  phoSigmaIEtaIEtaFull5x5.clear();
  phoR9.clear();
  phoR9Full5x5.clear();
  phoPFPhoIso.clear();
  phoPFNeuIso.clear();
  phoPFChIso.clear();
  phoSeedDet.clear();
  phoSeedSubdet.clear();
  phoSeedIeta.clear();
  phoSeedIphi.clear();
  phoSeedEta.clear();
  phoSeedPhi.clear();
  phoSeedRawID.clear();
  seedHcalIeta.clear();
  seedHcalIphi.clear();

  hcalRechitIeta.clear();
  hcalRechitIphi.clear();
  hcalRechitEnergy.clear();
  hcalRechitAbsDIetaFromPhoSeed.clear();
  hcalRechitAbsDIphiFromPhoSeed.clear();
  hcalRechitRawID.clear();
  hcalRechitDepth.clear();
  hcalRechitEta.clear();
  hcalRechitPhi.clear();

  edm::Handle<std::vector<PileupSummaryInfo> > genPileupHandle;
  iEvent.getByToken(puCollection_, genPileupHandle);
  
  if (genPileupHandle.isValid()) {
    for (std::vector<PileupSummaryInfo>::const_iterator pu = genPileupHandle->begin(); pu != genPileupHandle->end(); ++pu) {
      if (pu->getBunchCrossing() == 0) {
        puTrue = pu->getTrueNumInteractions();
        puObs = pu->getPU_NumInteractions();

        break;
      }
    }
  }

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken_, rhoHandle);
  if (!rhoHandle.failedToGet())
    rho = *(rhoHandle.product());
  else
    rho = 0.f;

  edm::Handle<HBHERecHitCollection> hbheRechitsHandle;
  iEvent.getByToken(hbhe_rechits_, hbheRechitsHandle);
  iSetup.get<CaloGeometryRecord>().get(theCaloGeometry);
  iSetup.get<CaloGeometryRecord>().get(towerMap_);

  edm::Handle<std::vector<reco::GenParticle> > genParticlesHandle;
  iEvent.getByToken(genParticlesCollection_, genParticlesHandle);

  edm::Handle<edm::ValueMap<float> > phoPhoIsoHandle;
  edm::Handle<edm::ValueMap<float> > phoNeuIsoHandle;
  edm::Handle<edm::ValueMap<float> > phoChaIsoHandle;
  iEvent.getByToken(phoPhoIsoToken_, phoPhoIsoHandle);
  iEvent.getByToken(phoNeuIsoToken_, phoNeuIsoHandle);
  iEvent.getByToken(phoChaIsoToken_, phoChaIsoHandle);

  int iPho = 0;
  for(const auto& pho : iEvent.get(phoToken_) ) {
    // std::cout << "\n new phoctron ...\n" ;
    int genmatched=0;
    double min_dr=9999.9;
    double ptR=9999.9;
     
    if (genParticlesHandle.isValid()) {
      //std::cout << "starting gen particle loop \n " ;
      for (std::vector<reco::GenParticle>::const_iterator ip = genParticlesHandle->begin(); ip != genParticlesHandle->end(); ++ip) {
	const reco::Candidate *p = (const reco::Candidate*)&(*ip);
	//std::cout << " p->pdgId() " << p->pdgId() << std::endl;
	if ( p->pdgId() == 22 ) { 
	  //std::cout << "-----  p->status() " << p->status() << " p->pdgId() " << p->pdgId() << std::endl;
	}
	if ( p->status() == 1 &&  p->pdgId() == 22  ) {
	  //std::cout << "checking if the reco pho match with this one" << std::endl;
	  double this_dr=reco::deltaR(pho,*p);
	  //std::cout << "this_dr " << this_dr << std::endl;
	  if (this_dr<min_dr) {
	    min_dr=this_dr;
	    ptR=pho.pt()/p->pt();
	  }
	}  
      }
    }
  
    if ( (min_dr<0.04) && (ptR>0.7) && (ptR<1.3) )  genmatched=1; // these cuts were decided looking at min_dr and ptR distributions.
    dR_recoPho_genPho.push_back(min_dr);
    ptRecoPho_by_ptGenPho.push_back(ptR);
    //  std::cout << "genmatched = " << genmatched <<  " min_dr " << min_dr << " ptR " << ptR   <<  std::endl;    
    pho_genmatched.push_back(genmatched);

    perPho_hcalRechitIeta.clear();
    perPho_hcalRechitIphi.clear();
    perPho_hcalRechitEnergy.clear();
    perPho_hcalRechitAbsDIetaFromPhoSeed.clear();
    perPho_hcalRechitAbsDIphiFromPhoSeed.clear();
    perPho_hcalRechitRawID.clear();
    perPho_hcalRechitDepth.clear();
    perPho_hcalRechitEta.clear();
    perPho_hcalRechitPhi.clear();

    phoScEta.push_back(pho.superCluster()->eta());
    phoPt.push_back(pho.pt());
    phoPhi.push_back(pho.phi());
    phoSigmaIEtaIEtaFull5x5.push_back(pho.full5x5_sigmaIetaIeta());
    phoR9.push_back(pho.r9());
    phoR9Full5x5.push_back(pho.full5x5_r9());

    phoPFPhoIso.push_back((*phoPhoIsoHandle)[iEvent.get(phoToken_).ptrAt(iPho)]);
    phoPFNeuIso.push_back((*phoNeuIsoHandle)[iEvent.get(phoToken_).ptrAt(iPho)]);
    phoPFChIso.push_back((*phoChaIsoHandle)[iEvent.get(phoToken_).ptrAt(iPho)]);

    EcalClusterLazyTools       lazyTool    (iEvent, iSetup, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);

    const reco::SuperCluster& superClus = *pho.superCluster();
    const reco::CaloCluster &seedCluster = *superClus.seed();
    DetId seedId = seedCluster.seed() ;
    phoSeedDet.push_back(seedId.det());
    phoSeedSubdet.push_back(seedId.subdetId());

    float var_phoSeedEta=-99;
    float var_phoSeedPhi=-99;

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
	    var_phoSeedEta=ecalrechitPoint.eta();
	    var_phoSeedPhi=ecalrechitPoint.phi();
	  }
	}
      }
    }
  
    phoSeedEta.push_back(var_phoSeedEta);
    phoSeedPhi.push_back(var_phoSeedPhi);
    
    //     if( seedId.det() == DetId::Forward ) return; // i guess this is not needed for ntuplizing purpose
  
    int var_phoSeedIeta=-999;
    int var_phoSeedIphi=-999;
    int var_phoSeedRawID=-999;
    
    int  var_seedHcalIeta=-999;
    int  var_seedHcalIphi=-999;
    
    if ( seedId.det() == DetId::Ecal ) {
      if (seedId.subdetId() == EcalBarrel) {
	EBDetId ebId(seedId);
	var_phoSeedIeta=ebId.ieta();
	var_phoSeedIphi=ebId.iphi();
	var_phoSeedRawID=ebId.rawId();       
      }
      
      else if (seedId.subdetId() == EcalEndcap) {
	EEDetId eeId(seedId);
	var_phoSeedIeta=eeId.ix();
        var_phoSeedIphi=eeId.iy();
        var_phoSeedRawID=eeId.rawId();

      }

      //Get hold of the seed hcal behind pho seed. Then loop over hcal recHits.
      //Condition to save HCAL recHits is:
      //it's close to pho seed and recHit energy greater than noise.
      //taken from RecoEgamma/EgammaIsolationAlgos/src/EGHcalRecHitSphoctor.cc //
      
      CaloTowerDetId towerId(towerMap_->towerOf(seedId)); 
      int seedHcalIEta = towerId.ieta();
      int seedHcalIPhi = towerId.iphi();
      
      var_seedHcalIeta=seedHcalIEta;
      var_seedHcalIphi=seedHcalIPhi;

      for (auto& hcalrh : iEvent.get(hbhe_rechits_) ) {
	int dIEtaAbs = std::abs(calDIEta(seedHcalIEta, hcalrh.id().ieta()));
	int dIPhiAbs = std::abs(calDIPhi(seedHcalIPhi, hcalrh.id().iphi()));

	if ( (dIEtaAbs <= maxDIEta_) && (dIPhiAbs <= maxDIPhi_) &&  (hcalrh.energy()>getMinEnergyHCAL_(hcalrh.id()) ) ) {
	  // std::cout << "close to pho, save " << std::endl;
	  perPho_hcalRechitIeta.push_back(hcalrh.id().ieta());
	  perPho_hcalRechitIphi.push_back(hcalrh.id().iphi());
	  perPho_hcalRechitEnergy.push_back(hcalrh.energy());
	  perPho_hcalRechitAbsDIetaFromPhoSeed.push_back(dIEtaAbs);
	  perPho_hcalRechitAbsDIphiFromPhoSeed.push_back(dIPhiAbs);
	 
	  perPho_hcalRechitRawID.push_back(hcalrh.id().rawId());
	  perPho_hcalRechitDepth.push_back(hcalrh.id().depth());

	  float rechitEta=-99;
	  float rechitPhi=-99;
	  if ( (hcalrh.id().rawId() != 0 ) ) {
	    if (theCaloGeometry.product() != nullptr) {
	      const CaloSubdetectorGeometry *geo = theCaloGeometry->getSubdetectorGeometry(hcalrh.id());
	      if(geo->getGeometry(hcalrh.id()) !=nullptr){
		const GlobalPoint & rechitPoint = theCaloGeometry.product()->getPosition(hcalrh.id());
		
		rechitEta=rechitPoint.eta();
		rechitPhi=rechitPoint.phi();
		
	      }
	    }
	  }
	  perPho_hcalRechitEta.push_back(rechitEta);
	  perPho_hcalRechitPhi.push_back(rechitPhi);
	}
      }
    }
  
    pho_haspixseed.push_back(pho.hasPixelSeed());
    pho_ispf.push_back(pho.isPFlowPhoton());
    pho_isstd.push_back(pho.isStandardPhoton());

    scEn.push_back(superClus.energy());
    phoSCRawEn.push_back(superClus.rawEnergy());
    seedEn.push_back(seedCluster.energy());
    seedEnCorr.push_back(seedCluster.correctedEnergy());
    cmssw_phoHoE.push_back(pho.hadTowOverEm());
   
    phoSeedIeta.push_back(var_phoSeedIeta);
    phoSeedIphi.push_back(var_phoSeedIphi);
    phoSeedRawID.push_back(var_phoSeedRawID);

    seedHcalIeta.push_back(var_seedHcalIeta);
    seedHcalIphi.push_back(var_seedHcalIphi);
    
    hcalRechitIeta.push_back(perPho_hcalRechitIeta);
    hcalRechitIphi.push_back(perPho_hcalRechitIphi);
    hcalRechitEnergy.push_back(perPho_hcalRechitEnergy);
    hcalRechitAbsDIetaFromPhoSeed.push_back(perPho_hcalRechitAbsDIetaFromPhoSeed);
    hcalRechitAbsDIphiFromPhoSeed.push_back(perPho_hcalRechitAbsDIphiFromPhoSeed);
    hcalRechitRawID.push_back(perPho_hcalRechitRawID);
    hcalRechitDepth.push_back(perPho_hcalRechitDepth);
    hcalRechitEta.push_back(perPho_hcalRechitEta);
    hcalRechitPhi.push_back(perPho_hcalRechitPhi);

    ++iPho;
  }
   
  tree->Fill();
  
  
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}

//doing some blatant copy paste from RecoEgamma/EgammaIsolationAlgos/src/EGHcalRecHitSphoctor.cc //

int PhoHoEAnalyzer::calDIPhi(int iPhi1, int iPhi2) {
  int dPhi = iPhi1 - iPhi2;
  if (dPhi > 72 / 2)
    dPhi -= 72;
  else if (dPhi < -72 / 2)
    dPhi += 72;
  return dPhi;
}

//

int PhoHoEAnalyzer::calDIEta(int iEta1, int iEta2) {
  int dEta = iEta1 - iEta2;
  if (iEta1 * iEta2 < 0) {  //-ve to +ve transistion and no crystal at zero
    if (dEta < 0)
      dEta++;
    else
      dEta--;
  }
  return dEta;
}

//HCAL thresholds from here https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/CaloTowersCreator/python/calotowermaker_cfi.py?%21v=CMSSW_10_6_2
//Note: As far as I understood, 
//for 2018, HB threshold is 0.7, and for Run 3 it becomes 0.1 in depth1, 0.2 in depth2, 0.3 in other depths.
//In HE, 2018 and Run3 is same, and it is 0.1 in depth1, and 0.2 in other depths.
//Double check these HCAL thresholds from Sam.

float PhoHoEAnalyzer::getMinEnergyHCAL_(HcalDetId id) const {
  if ( (id.subdetId() == HcalBarrel)  ) {
    if ( (Run2_2018 == 1) )
      return 0.7;
    else if ( (Run2_2018 == 0) ) { //means Run3
      if (id.depth() == 1)
	return 0.1;
      else if (id.depth() == 2)
	return 0.2;
      else
	return 0.3;
    }
    else //neither 2018 , nor Run3, not supported
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


// ------------ method called once each job just before starting event loop  ------------
void
PhoHoEAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
PhoHoEAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PhoHoEAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhoHoEAnalyzer);
