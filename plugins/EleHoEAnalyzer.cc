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
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
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


class EleHoEAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit EleHoEAnalyzer(const edm::ParameterSet&);
  ~EleHoEAnalyzer();
  
  static edm::ParameterSetDescription makePSetDescription();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  edm::Service<TFileService> fs;
  TTree   *tree = fs->make<TTree>("EventTree", "EventData");
  
  std::vector<int>  ele_golden;
  std::vector<int>  ele_unknown;
  std::vector<int>  ele_bigbrem;
  std::vector<int>  ele_badtrack;
  std::vector<int>  ele_showering;
  std::vector<int>  ele_gap;
  std::vector<float> ele_trackFbrem;
  std::vector<float> ele_superClusterFbrem;
  std::vector<int> ele_numberOfBrems;
  std::vector<int>  ele_genmatched;
  std::vector<float>  scEn;
  std::vector<float>  dR_recoEle_genEle;
  std::vector<float>  ptRecoEle_by_ptGenEle;
  std::vector<float>  eleSCRawEn;
  std::vector<float>  ecalEn;
  std::vector<float>  seedEn;
  std::vector<float>  seedEnCorr;
  std::vector<float>  cmssw_eleHoE;
  std::vector<float>  cmssw_eleHoE_full5x5;
  std::vector<float>  eleScEta;
  std::vector<float>  elePt;
  std::vector<float>  elePhi;
  std::vector<float>  eleSigmaIEtaIEtaFull5x5;
  std::vector<float>  elePFPhoIso;
  std::vector<float>  elePFNeuIso;
  std::vector<float>  elePFChIso;
  std::vector<int>    eleSeedDet;
  std::vector<int>    eleSeedSubdet;
  std::vector<int>    eleSeedIeta;
  std::vector<int>    eleSeedIphi;
  std::vector<float>  eleSeedEta;
  std::vector<float>  eleSeedPhi;
  std::vector<int>    eleSeedRawID;
  std::vector<int>    seedHcalIeta;
  std::vector<int>    seedHcalIphi;

  std::vector<std::vector<int>>    hcalRechitIeta;
  std::vector<std::vector<int>>    hcalRechitIphi;
  std::vector<std::vector<float>>  hcalRechitEnergy;
  std::vector<std::vector<int>>    hcalRechitAbsDIetaFromEleSeed;
  std::vector<std::vector<int>>    hcalRechitAbsDIphiFromEleSeed;
  std::vector<std::vector<int>>    hcalRechitRawID;
  std::vector<std::vector<int>>    hcalRechitDepth; // mostly for Run 3 //
  std::vector<std::vector<float>>    hcalRechitEta;
  std::vector<std::vector<float>>    hcalRechitPhi;

  std::vector<int>    perEle_hcalRechitIeta;
  std::vector<int>    perEle_hcalRechitIphi;
  std::vector<float>  perEle_hcalRechitEnergy;
  std::vector<int>    perEle_hcalRechitAbsDIetaFromEleSeed;
  std::vector<int>    perEle_hcalRechitAbsDIphiFromEleSeed;
  std::vector<int>    perEle_hcalRechitRawID;
  std::vector<int>    perEle_hcalRechitDepth; // mostly for Run 3 //
  std::vector<float>  perEle_hcalRechitEta;
  std::vector<float>  perEle_hcalRechitPhi;

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
EleHoEAnalyzer::EleHoEAnalyzer(const edm::ParameterSet& iConfig)
  :
  eleToken_(consumes<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("electrons"))),
  puCollection_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupCollection"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoSrc"))),
  hbhe_rechits_(consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("hbheInput"))),
  ebReducedRecHitCollection_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ebReducedRecHitCollection"))),
  eeReducedRecHitCollection_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("eeReducedRecHitCollection"))),
  esReducedRecHitCollection_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("esReducedRecHitCollection"))),
  genParticlesCollection_(consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleSrc"))),
  Run2_2018(iConfig.getParameter<bool>("Run2_2018_"))
{
  //now do what ever initialization is needed
  tree->Branch("ele_golden_",&ele_golden);
  tree->Branch("ele_unknown_",&ele_unknown);
  tree->Branch("ele_bigbrem_",&ele_bigbrem);
  tree->Branch("ele_gap_",&ele_gap);
  tree->Branch("ele_badtrack_",&ele_badtrack);
  tree->Branch("ele_showering_",&ele_showering);
  tree->Branch("ele_trackFbrem_",&ele_trackFbrem);
  tree->Branch("ele_superClusterFbrem_",&ele_superClusterFbrem);
  tree->Branch("ele_numberOfBrems_",&ele_numberOfBrems);
  tree->Branch("ele_genmatched_",&ele_genmatched);
  tree->Branch("dR_recoEle_genEle_",&dR_recoEle_genEle);
  tree->Branch("ptRecoEle_by_ptGenEle_",&ptRecoEle_by_ptGenEle);
  tree->Branch("scEn_",&scEn);
  tree->Branch("eleSCRawEn_",&eleSCRawEn);
  tree->Branch("ecalEn_",&ecalEn);
  tree->Branch("seedEn_",&seedEn);
  tree->Branch("seedEnCorr_",&seedEnCorr);
  tree->Branch("cmssw_eleHoE_",&cmssw_eleHoE);
  tree->Branch("cmssw_eleHoE_full5x5_",&cmssw_eleHoE_full5x5);
  tree->Branch("eleScEta_",&eleScEta);
  tree->Branch("elePt_",&elePt);
  tree->Branch("elePhi_",&elePhi);
  tree->Branch("eleSigmaIEtaIEtaFull5x5_",&eleSigmaIEtaIEtaFull5x5);
  tree->Branch("elePFPhoIso_",&elePFPhoIso);
  tree->Branch("elePFNeuIso_",&elePFNeuIso);
  tree->Branch("elePFChIso_",&elePFChIso);
  tree->Branch("eleSeedDet_",&eleSeedDet);
  tree->Branch("eleSeedSubdet_",&eleSeedSubdet);
  tree->Branch("eleSeedIeta_",&eleSeedIeta);
  tree->Branch("eleSeedIphi_",&eleSeedIphi);
  tree->Branch("eleSeedEta_",&eleSeedEta);
  tree->Branch("eleSeedPhi_",&eleSeedPhi);
  tree->Branch("eleSeedRawID_",&eleSeedRawID);
  tree->Branch("seedHcalIeta_",&seedHcalIeta);
  tree->Branch("seedHcalIphi_",&seedHcalIphi);
  tree->Branch("hcalRechitIeta_",&hcalRechitIeta);
  tree->Branch("hcalRechitIphi_",&hcalRechitIphi);
  tree->Branch("hcalRechitEnergy_",&hcalRechitEnergy);
  tree->Branch("hcalRechitAbsDIetaFromEleSeed_",&hcalRechitAbsDIetaFromEleSeed);
  tree->Branch("hcalRechitAbsDIphiFromEleSeed_",&hcalRechitAbsDIphiFromEleSeed);
  tree->Branch("hcalRechitRawID_",&hcalRechitRawID);
  tree->Branch("hcalRechitDepth_",&hcalRechitDepth);
  tree->Branch("hcalRechitEta_",&hcalRechitEta);
  tree->Branch("hcalRechitPhi_",&hcalRechitPhi);
  tree->Branch("puTrue_", &puTrue);
  tree->Branch("puObs_", &puObs);
  tree->Branch("rho_", &rho);
}


EleHoEAnalyzer::~EleHoEAnalyzer()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called for each event  ------------
void
EleHoEAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //   std::cout << " \n ****** new event ... " << std::endl; 
  using namespace edm;
  
  ele_golden.clear();
  ele_unknown.clear();
  ele_badtrack.clear();
  ele_bigbrem.clear();
  ele_showering.clear();
  ele_gap.clear();

  ele_trackFbrem.clear();
  ele_superClusterFbrem.clear();
  ele_numberOfBrems.clear();

  ele_genmatched.clear();
  ptRecoEle_by_ptGenEle.clear();
  dR_recoEle_genEle.clear();
  scEn.clear();
  eleSCRawEn.clear();
  ecalEn.clear();
  seedEn.clear();
  seedEnCorr.clear();
  cmssw_eleHoE.clear();
  cmssw_eleHoE_full5x5.clear();
  eleScEta.clear();
  elePt.clear();
  elePhi.clear();
  eleSigmaIEtaIEtaFull5x5.clear();
  elePFPhoIso.clear();
  elePFNeuIso.clear();
  elePFChIso.clear();
  eleSeedDet.clear();
  eleSeedSubdet.clear();
  eleSeedIeta.clear();
  eleSeedIphi.clear();
  eleSeedEta.clear();
  eleSeedPhi.clear();
  eleSeedRawID.clear();
  seedHcalIeta.clear();
  seedHcalIphi.clear();

  hcalRechitIeta.clear();
  hcalRechitIphi.clear();
  hcalRechitEnergy.clear();
  hcalRechitAbsDIetaFromEleSeed.clear();
  hcalRechitAbsDIphiFromEleSeed.clear();
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

  for(const auto& ele : iEvent.get(eleToken_) ) {
    // std::cout << "\n new electron ...\n" ;
    int genmatched=0;
    double min_dr=9999.9;
    double ptR=9999.9;
     
    if (genParticlesHandle.isValid()) {
      //std::cout << "starting gen particle loop \n " ;
      for (std::vector<reco::GenParticle>::const_iterator ip = genParticlesHandle->begin(); ip != genParticlesHandle->end(); ++ip) {
	const reco::Candidate *p = (const reco::Candidate*)&(*ip);
	//std::cout << " p->pdgId() " << p->pdgId() << std::endl;
	if ( (abs(p->pdgId())) ==11 ) { 
	  //std::cout << "-----  p->status() " << p->status() << " p->pdgId() " << p->pdgId() << std::endl;
	}
	if ( (p->status()==1) &&  ((abs(p->pdgId())) == 11)  ) {
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
  
    if ( (min_dr<0.04) && (ptR>0.7) && (ptR<1.3) )  genmatched=1; // these cuts were decided looking at min_dr and ptR distributions.
    dR_recoEle_genEle.push_back(min_dr);
    ptRecoEle_by_ptGenEle.push_back(ptR);
    //  std::cout << "genmatched = " << genmatched <<  " min_dr " << min_dr << " ptR " << ptR   <<  std::endl;    
    ele_genmatched.push_back(genmatched);

    perEle_hcalRechitIeta.clear();
    perEle_hcalRechitIphi.clear();
    perEle_hcalRechitEnergy.clear();
    perEle_hcalRechitAbsDIetaFromEleSeed.clear();
    perEle_hcalRechitAbsDIphiFromEleSeed.clear();
    perEle_hcalRechitRawID.clear();
    perEle_hcalRechitDepth.clear();
    perEle_hcalRechitEta.clear();
    perEle_hcalRechitPhi.clear();

    eleScEta.push_back(ele.superCluster()->eta());
    elePt.push_back(ele.pt());
    elePhi.push_back(ele.phi());
    eleSigmaIEtaIEtaFull5x5.push_back(ele.full5x5_sigmaIetaIeta());
  
    reco::GsfElectron::PflowIsolationVariables pfIso = ele.pfIsolationVariables();
    elePFPhoIso.push_back(pfIso.sumPhotonEt);
    elePFNeuIso.push_back(pfIso.sumNeutralHadronEt);
    elePFChIso.push_back(pfIso.sumChargedHadronPt);

    EcalClusterLazyTools       lazyTool    (iEvent, iSetup, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);

    const reco::SuperCluster& superClus = *ele.superCluster();
    const reco::CaloCluster &seedCluster = *superClus.seed();
    DetId seedId = seedCluster.seed() ;
    eleSeedDet.push_back(seedId.det());
    eleSeedSubdet.push_back(seedId.subdetId());

    float var_eleSeedEta=-99;
    float var_eleSeedPhi=-99;

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
	    var_eleSeedEta=ecalrechitPoint.eta();
	    var_eleSeedPhi=ecalrechitPoint.phi();
	  }
	}
      }
    }
  
    eleSeedEta.push_back(var_eleSeedEta);
    eleSeedPhi.push_back(var_eleSeedPhi);
    
    //     if( seedId.det() == DetId::Forward ) return; // i guess this is not needed for ntuplizing purpose
  
    int var_eleSeedIeta=-999;
    int var_eleSeedIphi=-999;
    int var_eleSeedRawID=-999;
    
    int  var_seedHcalIeta=-999;
    int  var_seedHcalIphi=-999;
    
    if ( seedId.det() == DetId::Ecal ) {
      if (seedId.subdetId() == EcalBarrel) {
	EBDetId ebId(seedId);
	var_eleSeedIeta=ebId.ieta();
	var_eleSeedIphi=ebId.iphi();
	var_eleSeedRawID=ebId.rawId();       
      }
      
      else if (seedId.subdetId() == EcalEndcap) {
	EEDetId eeId(seedId);
	var_eleSeedIeta=eeId.ix();
        var_eleSeedIphi=eeId.iy();
        var_eleSeedRawID=eeId.rawId();

      }

      //Get hold of the seed hcal behind ele seed. Then loop over hcal recHits.
      //Condition to save HCAL recHits is:
      //it's close to ele seed and recHit energy greater than noise.
      //taken from RecoEgamma/EgammaIsolationAlgos/src/EGHcalRecHitSelector.cc //
      
      CaloTowerDetId towerId(towerMap_->towerOf(seedId)); 
      int seedHcalIEta = towerId.ieta();
      int seedHcalIPhi = towerId.iphi();
      
      var_seedHcalIeta=seedHcalIEta;
      var_seedHcalIphi=seedHcalIPhi;

      for (auto& hcalrh : iEvent.get(hbhe_rechits_) ) {
	int dIEtaAbs = std::abs(calDIEta(seedHcalIEta, hcalrh.id().ieta()));
	int dIPhiAbs = std::abs(calDIPhi(seedHcalIPhi, hcalrh.id().iphi()));

	if ( (dIEtaAbs <= maxDIEta_) && (dIPhiAbs <= maxDIPhi_) &&  (hcalrh.energy()>getMinEnergyHCAL_(hcalrh.id()) ) ) {
	  // std::cout << "close to ele, save " << std::endl;
	  perEle_hcalRechitIeta.push_back(hcalrh.id().ieta());
	  perEle_hcalRechitIphi.push_back(hcalrh.id().iphi());
	  perEle_hcalRechitEnergy.push_back(hcalrh.energy());
	  perEle_hcalRechitAbsDIetaFromEleSeed.push_back(dIEtaAbs);
	  perEle_hcalRechitAbsDIphiFromEleSeed.push_back(dIPhiAbs);
	 
	  perEle_hcalRechitRawID.push_back(hcalrh.id().rawId());
	  perEle_hcalRechitDepth.push_back(hcalrh.id().depth());

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
	  perEle_hcalRechitEta.push_back(rechitEta);
	  perEle_hcalRechitPhi.push_back(rechitPhi);
	}
      }
    }
  
    ele_trackFbrem.push_back(ele.trackFbrem());
    ele_superClusterFbrem.push_back(ele.superClusterFbrem());
    ele_numberOfBrems.push_back(ele.numberOfBrems());
    
    //    std::cout << "track_fbrem=" << ele.trackFbrem() << " sc_fbrem=" << ele.superClusterFbrem() << " nBrem="  <<  ele.numberOfBrems()  << std::endl;
    // std::cout << "classification " << ele.classificationVariables().ClassificationVariables << std::endl;
    int var_golden=0;
    int var_unknown=0;
    int var_gap=0;
    int var_badtrack=0;
    int var_showering=0;
    int var_bigbrem=0;

    if (ele.classification() == reco::GsfElectron::GOLDEN) {
      // std::cout << "golden ele!!!" << std::endl;
      var_golden=1;
    }

    if (ele.classification() == reco::GsfElectron::UNKNOWN) {
      // std::cout << "UNKNOWN ele!!!" << std::endl;
      var_unknown=1;
    }

    if (ele.classification() == reco::GsfElectron::BIGBREM) {
      //std::cout << "BIGBREM ele!!!" << std::endl;
      var_bigbrem=1;
    }

    if (ele.classification() == reco::GsfElectron::BADTRACK) {
      //std::cout << "BADTRACK ele!!!" << std::endl;
      var_badtrack=1;
    }

    if (ele.classification() == reco::GsfElectron::SHOWERING) {
      // std::cout << "SHOWERING ele!!!" << std::endl;
      var_showering=1;
    }

    if (ele.classification() == reco::GsfElectron::GAP) {
      // std::cout << "GAP ele!!!" << std::endl;
      var_gap=1;
    }

    ele_golden.push_back(var_golden);
    ele_unknown.push_back(var_unknown);
    ele_gap.push_back(var_gap);
    ele_badtrack.push_back(var_badtrack);
    ele_showering.push_back(var_showering);
    ele_bigbrem.push_back(var_bigbrem);

    scEn.push_back(superClus.energy());
    eleSCRawEn.push_back(superClus.rawEnergy());
    seedEn.push_back(seedCluster.energy());
    seedEnCorr.push_back(seedCluster.correctedEnergy());
    ecalEn.push_back(ele.ecalEnergy());
    cmssw_eleHoE.push_back(ele.hcalOverEcal());
    cmssw_eleHoE_full5x5.push_back(ele.full5x5_hcalOverEcal());
   
    eleSeedIeta.push_back(var_eleSeedIeta);
    eleSeedIphi.push_back(var_eleSeedIphi);
    eleSeedRawID.push_back(var_eleSeedRawID);

    seedHcalIeta.push_back(var_seedHcalIeta);
    seedHcalIphi.push_back(var_seedHcalIphi);
    
    hcalRechitIeta.push_back(perEle_hcalRechitIeta);
    hcalRechitIphi.push_back(perEle_hcalRechitIphi);
    hcalRechitEnergy.push_back(perEle_hcalRechitEnergy);
    hcalRechitAbsDIetaFromEleSeed.push_back(perEle_hcalRechitAbsDIetaFromEleSeed);
    hcalRechitAbsDIphiFromEleSeed.push_back(perEle_hcalRechitAbsDIphiFromEleSeed);
    hcalRechitRawID.push_back(perEle_hcalRechitRawID);
    hcalRechitDepth.push_back(perEle_hcalRechitDepth);
    hcalRechitEta.push_back(perEle_hcalRechitEta);
    hcalRechitPhi.push_back(perEle_hcalRechitPhi);
  }
   
  tree->Fill();
  
  
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}

//doing some blatant copy paste from RecoEgamma/EgammaIsolationAlgos/src/EGHcalRecHitSelector.cc //

int EleHoEAnalyzer::calDIPhi(int iPhi1, int iPhi2) {
  int dPhi = iPhi1 - iPhi2;
  if (dPhi > 72 / 2)
    dPhi -= 72;
  else if (dPhi < -72 / 2)
    dPhi += 72;
  return dPhi;
}

//

int EleHoEAnalyzer::calDIEta(int iEta1, int iEta2) {
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

float EleHoEAnalyzer::getMinEnergyHCAL_(HcalDetId id) const {
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
EleHoEAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
EleHoEAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
EleHoEAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EleHoEAnalyzer);
