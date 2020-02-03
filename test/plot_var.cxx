// -*- C++ -*-
// root -l -b plot_var.cxx++

//#include "/home/afiqaize/Downloads/eft/EFTFitter/src/PlotterUtil.h"
#include "/afs/cern.ch/work/a/afiqaize/EFTFitter/src/PlotterUtil.h"
#include "TChain.h"

constexpr double etaEB = 1.4442, etaET = 1.566, etaEE = 2.5;

double deltaPhi(const double phi1, const double phi2) {
  return std::acos( std::cos(phi1 - phi2) );
}

double deltaR(const double eta1, const double phi1, const double eta2, const double phi2) {
  return std::sqrt( std::pow(eta1 - eta2, 2.) + std::pow( deltaPhi(phi1, phi2), 2.) );
}

// wrapper to reduce line count
// assumes TH1D... oh well
void assign_plot(const std::string &var, const int nBin, const double min, const double max, 
                 std::vector<std::pair<int, std::array<std::string, 3>>> &vp_plot,
                 std::vector<std::vector<Plot<TH1>>> &v_hist) {
  v_hist.emplace_back( std::vector<Plot<TH1>>(vp_plot.size()) );
  for (int iP = 0; iP < vp_plot.size(); ++iP) {
    if (iP == 0)
      v_hist.back()[iP].plot = std::make_unique<TH1D>(var.c_str(), "", nBin, min, max);
    else
      v_hist.back()[iP].plot = std::make_unique<TH1D>((var + toStr(iP)).c_str(), "", nBin, min, max);

    stylePlot(v_hist.back()[iP].plot.get(), vp_plot[iP].first, 1., 0, 20 + iP, 1.5, 1, 3);

    v_hist.back()[iP].legend_txt = vp_plot[iP].second[0];
    v_hist.back()[iP].legend_opt = vp_plot[iP].second[1];
    v_hist.back()[iP].draw_opt = vp_plot[iP].second[2];
  }
}

// another wrapper to reduce repeating arguments
void std_plot(const std::vector<Plot<TH1>> &v_hist, const std::vector<Plot<TGraph>> &v_graph,
              const std::string &filename, const bool logY,
              const double axy_min, const double axy_max,
              const std::string &axy_txt,
              const double axx_min, const double axx_max,
              const std::string &axx_txt) {
  standard_plot(v_hist, v_graph, filename, // v_hist, v_graph, filename (w/o ext)
                true, true, false, logY, // add_uoflow, normalize_shape, logX, logY
                1, 0, 0, 42, 0.041, "", // legend style, header
                0., 0., 0., 0., // legend position
                axy_min, axy_max, axy_txt, 0.043, 1.13, 0.037, // y axis
                axx_min, axx_max, axx_txt, 0.043, 1.11, 0.037, // x axis
                0.045, 0.11, 0.10, 0.035 // canvas margins
                , ".png" // plot format
                );
}

void plot_var() {
  gROOT->Reset();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2(true);
  TGaxis::SetMaxDigits(3);

  /*/ dump in the root files
  const std::vector<std::string> v_file = {"/home/afiqaize/Downloads/EGM/e_106X/hcal_depth_200122/test_MINIAOD_2023.root", 
                                           "/home/afiqaize/Downloads/EGM/e_106X/hcal_depth_200122/test_MINIAOD_2023.root"};
  const std::vector<double> v_wgt = {1., 1.};
  */

  const std::string idir = "/eos/cms/store/group/phys_egamma/afiqaize/hcal_depth_200129/";
  const std::vector<std::string> v_file = {"DYToEE_M-50_NNPDF31_TuneCP5_14TeV-powheg-pythia8/crab_dy_ee/200129_173513/0000/", 
                                           "GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_14TeV_Pythia8/crab_yjet_pt_40_myy_80/200129_173623/0000/",
                                           "TTToSemiLeptonic_TuneCP5_14TeV-powheg-pythia8/crab_tt_lj/200129_173403/0000/",
                                           "QCD_Pt_15to30_TuneCP5_14TeV_pythia8/crab_qcd_pt_15_30/200129_173422/0000/", 
                                           "QCD_Pt_30to50_TuneCP5_14TeV_pythia8/crab_qcd_pt_30_50/200129_173609/0000/",
                                           "QCD_Pt_50to80_TuneCP5_14TeV_pythia8/crab_qcd_pt_50_80/200129_173350/0000/",
                                           "QCD_Pt_80to120_TuneCP5_14TeV_pythia8/crab_qcd_pt_80_120/200129_173534/0000/", 
                                           "QCD_Pt_120to170_TuneCP5_14TeV_pythia8/crab_qcd_pt_120_170/200129_173542/0000/",
                                           "QCD_Pt_170to300_TuneCP5_14TeV_pythia8/crab_qcd_pt_170_300/200129_173639/0000/",
                                           "QCD_Pt_300to470_TuneCP5_14TeV_pythia8/crab_qcd_pt_300_470/200129_173432/0000/", 
                                           "QCD_Pt_470to600_TuneCP5_14TeV_pythia8/crab_qcd_pt_470_600/200129_173444/0000/",
                                           "QCD_Pt_600to800_TuneCP5_14TeV_pythia8/crab_qcd_pt_600_800/200129_173551/0000/",
                                           "QCD_Pt_800to1000_TuneCP5_14TeV_pythia8/crab_qcd_pt_800_1000/200129_173410/0000/", 
                                           "QCD_Pt_1000to1400_TuneCP5_14TeV_pythia8/crab_qcd_pt_1000_1400/200129_173615/0000/",
                                           "QCD_Pt_1400to1800_TuneCP5_14TeV_pythia8/crab_qcd_pt_1400_1800/200129_173453/0000/",
                                           "QCD_Pt_1800to2400_TuneCP5_14TeV_pythia8/crab_qcd_pt_1800_2400/200129_173504/0000/", 
                                           "QCD_Pt_2400to3200_TuneCP5_14TeV_pythia8/crab_qcd_pt_2400_3200/200129_173523/0000/",
                                           "QCD_Pt_3200toInf_TuneCP5_14TeV_pythia8/crab_qcd_pt_3200/200129_173651/0000/"};
  // ensure same order as above
  // weight = xsec / nevt, xsec just taken off xsecdb, nevt read off files
  // 13 TeV xsec since some aren't available for 14 TeV, but fine
  const std::vector<double> v_xsec = {2137.,
                                      874.2,
                                      831.76 * 4. / 9.,
                                      0., //1.248e9,
                                      0., //1.066e8,
                                      1.568e7,
                                      2.336e6,
                                      4.068e5,
                                      1.033e5,
                                      6826.,
                                      552.6,
                                      156.6,
                                      26.32,
                                      7.5,
                                      0.6479,
                                      0.08715,
                                      0.005242,
                                      1.349e-4};

  // mark in the set of files to go into a single plot
  // key is color, value is legend text, legend type, draw option
  std::vector<std::pair<int, std::array<std::string, 3>>> vp_plot = { {kPink + 9, {"DY", "lp", "hist e2"}},
                                                                      {kPink + 8, {"#gamma + j", "lp", "hist e2"}},
                                                                      {kPink + 1, {"t#bar{t} (s)", "lp", "hist e2"}},
                                                                      {kAzure + 8, {"t#bar{t} (b)", "lp", "hist e2"}},
                                                                      {kAzure + 4, {"QCD", "lp", "hist e2"}},
                                                                      {kPink - 1, {"Signal", "lp", "hist e2"}},
                                                                      {kAzure - 1, {"Background", "lp", "hist e2"}}};

  // list of plots to make
  double xmax = 0.;
  std::vector<std::vector<Plot<TH1>>> v_hist;
  assign_plot("sc_eta_hist", 50, -2.5, 2.5, vp_plot, v_hist); // 0
  for (int iV = 0; iV < 7; ++iV)
    assign_plot("hsum_1depth_" + toStr(iV + 1) + "_eb_hist", 100, 0., 50., vp_plot, v_hist); // 1 - 7
  for (int iV = 0; iV < 7; ++iV)
    assign_plot("hsum_1depth_" + toStr(iV + 1) + "_ee_hist", 100, 0., 200., vp_plot, v_hist); // 8 - 14

  for (int iV = 0; iV < 7; ++iV) {
    if (iV == 1)
      xmax = 2.;
    else if (iV > 3)
      xmax = 0.1;
    else 
      xmax = 1.;

    assign_plot("hoe_1depth_" + toStr(iV + 1) + "_eb_hist", 100, 0., xmax, vp_plot, v_hist); // 15 - 21
  }
  for (int iV = 0; iV < 7; ++iV) {
    if (iV == 3 or iV == 4)
      xmax = 1.;
    else if (iV > 4)
      xmax = 0.4;
    else 
      xmax = 2.;

    assign_plot("hoe_1depth_" + toStr(iV + 1) + "_ee_hist", 100, 0., xmax, vp_plot, v_hist); // 22 - 28
  }

  assign_plot("hsum_7depth_17_eb_hist", 100, 0., 100., vp_plot, v_hist); // 29
  assign_plot("hsum_7depth_17_ee_hist", 100, 0., 400., vp_plot, v_hist); // 30
  assign_plot("hoe_7depth_17_eb_hist", 100, 0., 3., vp_plot, v_hist); // 31
  assign_plot("hoe_7depth_17_ee_hist", 100, 0., 5., vp_plot, v_hist); // 32

  assign_plot("hoe_cmssw_eb_hist", 100, 0., 2., vp_plot, v_hist); // 33
  assign_plot("hoe_cmssw_ee_hist", 100, 0., 2., vp_plot, v_hist); // 34
  assign_plot("hoe_5x5_eb_hist", 100, 0., 2., vp_plot, v_hist); // 35
  assign_plot("hoe_5x5_ee_hist", 100, 0., 2., vp_plot, v_hist); // 36

  for (int iV = 0; iV < 6; ++iV)
    assign_plot("hsum_2depth_" + toStr(iV + 1) + toStr(iV + 2) + "_eb_hist", 100, 0., 100., vp_plot, v_hist); // 37 - 42
  for (int iV = 0; iV < 6; ++iV)
    assign_plot("hsum_2depth_" + toStr(iV + 1) + toStr(iV + 2) + "_ee_hist", 100, 0., 400., vp_plot, v_hist); // 43 - 48

  for (int iV = 0; iV < 6; ++iV) {
    if (iV == 3)
      xmax = 1.;
    else if (iV > 3)
      xmax = 0.1;
    else 
      xmax = 2.;

    assign_plot("hoe_2depth_" + toStr(iV + 1) + toStr(iV + 2) + "_eb_hist", 100, 0., xmax, vp_plot, v_hist); // 49 - 54
  }
  for (int iV = 0; iV < 6; ++iV) {
    if (iV == 3)
      xmax = 2.;
    else if (iV == 4)
      xmax = 1.;
    else if (iV == 5)
      xmax = 0.5;
    else 
      xmax = 4.;

    assign_plot("hoe_2depth_" + toStr(iV + 1) + toStr(iV + 2) + "_ee_hist", 100, 0., xmax, vp_plot, v_hist); // 55 - 60
  }

  for (int iV = 0; iV < 5; ++iV)
    assign_plot("hsum_3depth_" + toStr(iV + 1) + toStr(iV + 3) + "_eb_hist", 100, 0., 100., vp_plot, v_hist); // 61 - 65
  for (int iV = 0; iV < 5; ++iV)
    assign_plot("hsum_3depth_" + toStr(iV + 1) + toStr(iV + 3) + "_ee_hist", 100, 0., 400., vp_plot, v_hist); // 66 - 70

  for (int iV = 0; iV < 5; ++iV) {
    if (iV == 0 or iV == 1)
      xmax = 3.;
    else if (iV == 3)
      xmax = 1.;
    else if (iV == 4)
      xmax = 0.1;
    else 
      xmax = 2.;

    assign_plot("hoe_3depth_" + toStr(iV + 1) + toStr(iV + 3) + "_eb_hist", 100, 0., xmax, vp_plot, v_hist); // 71 - 75
  }
  for (int iV = 0; iV < 5; ++iV) {
    if (iV == 0 or iV == 1)
      xmax = 5.;
    else if (iV == 3)
      xmax = 2.;
    else if (iV == 4)
      xmax = 1.;
    else 
      xmax = 3.;

    assign_plot("hoe_3depth_" + toStr(iV + 1) + toStr(iV + 3) + "_ee_hist", 100, 0., xmax, vp_plot, v_hist); // 76 - 80
  }

  for (int iV = 0; iV < 4; ++iV)
    assign_plot("hsum_4depth_" + toStr(iV + 1) + toStr(iV + 4) + "_eb_hist", 100, 0., 100., vp_plot, v_hist); // 81 - 84
  for (int iV = 0; iV < 4; ++iV)
    assign_plot("hsum_4depth_" + toStr(iV + 1) + toStr(iV + 4) + "_ee_hist", 100, 0., 400., vp_plot, v_hist); // 85 - 88

  for (int iV = 0; iV < 4; ++iV) {
    if (iV == 0 or iV == 1)
      xmax = 3.;
    else if (iV == 3)
      xmax = 1.;
    else 
      xmax = 2.;

    assign_plot("hoe_4depth_" + toStr(iV + 1) + toStr(iV + 4) + "_eb_hist", 100, 0., xmax, vp_plot, v_hist); // 89 - 92
  }
  for (int iV = 0; iV < 4; ++iV) {
    if (iV == 0 or iV == 1)
      xmax = 5.;
    else if (iV == 3)
      xmax = 2.;
    else 
      xmax = 4.;

    assign_plot("hoe_4depth_" + toStr(iV + 1) + toStr(iV + 4) + "_ee_hist", 100, 0., xmax, vp_plot, v_hist); // 93 - 96
  }

  for (int iV = 0; iV < 3; ++iV)
    assign_plot("hsum_5depth_" + toStr(iV + 1) + toStr(iV + 5) + "_eb_hist", 100, 0., 100., vp_plot, v_hist); // 97 - 99
  for (int iV = 0; iV < 3; ++iV)
    assign_plot("hsum_5depth_" + toStr(iV + 1) + toStr(iV + 5) + "_ee_hist", 100, 0., 400., vp_plot, v_hist); // 100 - 102

  for (int iV = 0; iV < 3; ++iV) {
    if (iV == 0 or iV == 1)
      xmax = 3.;
    else 
      xmax = 2.;

    assign_plot("hoe_5depth_" + toStr(iV + 1) + toStr(iV + 5) + "_eb_hist", 100, 0., xmax, vp_plot, v_hist); // 103 - 105
  }
  for (int iV = 0; iV < 3; ++iV) {
    if (iV == 0 or iV == 1)
      xmax = 5.;
    else 
      xmax = 3.;

    assign_plot("hoe_5depth_" + toStr(iV + 1) + toStr(iV + 5) + "_ee_hist", 100, 0., xmax, vp_plot, v_hist); // 106 - 108
  }

  for (int iV = 0; iV < 2; ++iV)
    assign_plot("hsum_6depth_" + toStr(iV + 1) + toStr(iV + 6) + "_eb_hist", 100, 0., 100., vp_plot, v_hist); // 109 - 110
  for (int iV = 0; iV < 2; ++iV)
    assign_plot("hsum_6depth_" + toStr(iV + 1) + toStr(iV + 6) + "_ee_hist", 100, 0., 400., vp_plot, v_hist); // 111 - 112

  for (int iV = 0; iV < 2; ++iV) {
    if (iV == 0 or iV == 1)
      xmax = 3.;
    else 
      xmax = 2.;

    assign_plot("hoe_6depth_" + toStr(iV + 1) + toStr(iV + 6) + "_eb_hist", 100, 0., xmax, vp_plot, v_hist); // 113 - 114
  }
  for (int iV = 0; iV < 2; ++iV) {
    if (iV == 0 or iV == 1)
      xmax = 5.;
    else 
      xmax = 4.;

    assign_plot("hoe_6depth_" + toStr(iV + 1) + toStr(iV + 6) + "_ee_hist", 100, 0., xmax, vp_plot, v_hist); // 115 - 116
  }

  for (int iV = 0; iV < 6; ++iV) {
    for (int iW = iV + 1; iW < 7; ++iW) {
      if (iV == 0 and iW == 1)
        xmax = 20.;
      else if (iV == 0 and (iW == 2 or iW == 3))
        xmax = 15.;
      else if (iV == 1 and iW == 2)
        xmax = 10.;
      else if (((iV == 1 or iV == 2) and iW == 3) or (iV == 3 and iW == 4) or (iW > 4))
        xmax = 8.;
      else 
        xmax = 4.;

      assign_plot("hoh_" + toStr(iW + 1) + "_" + toStr(iV + 1) + "_eb_hist", 100, 0., xmax, vp_plot, v_hist); // 117 - 137
    }
  }
  for (int iV = 0; iV < 6; ++iV) {
    for (int iW = iV + 1; iW < 7; ++iW) {
      if ((iV == 0 and iW == 1) or (iV == 3 and iW == 4))
        xmax = 5.;
      else if ((iV == 0 and iW == 2) or (iV == 2 and iW == 3) or (iV == 3 and iW == 5) or (iV == 4 and iW == 6))
        xmax = 3.;
      else if ((iV == 2 and iW == 5) or (iV == 3 and iW == 6))
        xmax = 2.;
      else if ((iV == 1 and iW == 2) or (iV == 0 and iW == 3) or ((iV == 0 or iV == 2) and iW == 4))
        xmax = 1.5;
      else if ((iV == 1 and iW == 4) or (iV < 2 and iW == 5) or (iV < 3 and iW == 6) or (iV == 1 and iW == 3))
        xmax = 1.;
      else 
        xmax = 4.;

      assign_plot("hoh_" + toStr(iW + 1) + "_" + toStr(iV + 1) + "_ee_hist", 100, 0., xmax, vp_plot, v_hist); // 138 - 158
    }
  }

  for (int iV = 0; iV < 5; ++iV) {
    for (int iW = iV + 2; iW < 6; ++iW) {
      if (iV == 0 and (iW == 2 or iW == 3))
        xmax = 5.;
      else if ((iV == 1 and iW == 3) or (iV == 3 and iW == 5))
        xmax = 4.;
      else if ((iV < 2 and iW == 4) or (iV < 3 and iW == 5))
        xmax = 1.;
      else 
        xmax = 2.;

      const std::string sV = toStr(iV + 1) + toStr(iV + 2), sW = toStr(iW + 1) + toStr(iW + 2);
      assign_plot("hoh_" + sW + "_" + sV + "_eb_hist", 100, 0., xmax, vp_plot, v_hist); // 159 - 168
    }
  }
  for (int iV = 0; iV < 5; ++iV) {
    for (int iW = iV + 2; iW < 6; ++iW) {
      if ((iV == 0 and iW == 2) or (iV == 2 and iW == 4) or (iV == 2 and iW == 5))
        xmax = 1.5;
      else if (((iV == 0 or iV == 1) and iW == 3) or (iV < 2 and iW == 4))
        xmax = 1.;
      else if (iV < 2 and iW == 5)
        xmax = 0.5;
      else 
        xmax = 2.;

      const std::string sV = toStr(iV + 1) + toStr(iV + 2), sW = toStr(iW + 1) + toStr(iW + 2);
      assign_plot("hoh_" + sW + "_" + sV + "_ee_hist", 100, 0., xmax, vp_plot, v_hist); // 169 - 178
    }
  }

  // put the files in chain and remember the cumulative entry count
  std::vector<long long> v_entry;
  auto t1 = std::make_unique<TChain>("demo/EventTree");
  for (auto &f : v_file) {
    //t1->Add(f.c_str());
    t1->Add((idir + f + "test_*.root").c_str());
    v_entry.push_back(t1->GetEntries());
  }

  std::vector<int> *v_match_rg = nullptr;
  t1->SetBranchAddress("ele_genmatched_", &v_match_rg);
  std::vector<float> *v_dR_rg = nullptr;
  t1->SetBranchAddress("dR_recoEle_genEle_", &v_dR_rg);
  std::vector<float> *v_pt_ratio_rg = nullptr;
  t1->SetBranchAddress("ptRecoEle_by_ptGenEle_", &v_pt_ratio_rg);

  std::vector<float> *v_sc_eta = nullptr;
  t1->SetBranchAddress("eleScEta_", &v_sc_eta);
  std::vector<float> *v_pt = nullptr;
  t1->SetBranchAddress("elePt_", &v_pt);
  std::vector<float> *v_phi = nullptr;
  t1->SetBranchAddress("elePhi_", &v_phi);

  std::vector<float> *v_sc_e = nullptr;
  t1->SetBranchAddress("scEn_", &v_sc_e);
  std::vector<float> *v_sc_raw_e = nullptr;
  t1->SetBranchAddress("eleSCRawEn_", &v_sc_raw_e);
  std::vector<float> *v_ecal_e = nullptr;
  t1->SetBranchAddress("ecalEn_", &v_ecal_e);
  std::vector<float> *v_seed_e = nullptr;
  t1->SetBranchAddress("seedEn_", &v_seed_e);
  std::vector<float> *v_seed_corr_e = nullptr;
  t1->SetBranchAddress("seedEnCorr_", &v_seed_corr_e);
  std::vector<float> *v_hoe = nullptr;
  t1->SetBranchAddress("cmssw_eleHoE_", &v_hoe);
  std::vector<float> *v_hoe5x5 = nullptr;
  t1->SetBranchAddress("cmssw_eleHoE_full5x5_", &v_hoe5x5);

  std::vector<float> *v_sieie = nullptr;
  t1->SetBranchAddress("eleSigmaIEtaIEtaFull5x5_", &v_sieie);
  std::vector<float> *v_pfiso_neu = nullptr;
  t1->SetBranchAddress("elePFNeuIso_", &v_pfiso_neu);
  std::vector<float> *v_pfiso_cha = nullptr;
  t1->SetBranchAddress("elePFChIso_", &v_pfiso_cha);

  std::vector<int> *v_seed_det = nullptr;
  t1->SetBranchAddress("eleSeedDet_", &v_seed_det);
  std::vector<int> *v_seed_subdet = nullptr;
  t1->SetBranchAddress("eleSeedSubdet_", &v_seed_subdet);
  std::vector<int> *v_seed_rawid = nullptr;
  t1->SetBranchAddress("eleSeedRawID_", &v_seed_rawid);
  std::vector<int> *v_seed_ieta = nullptr;
  t1->SetBranchAddress("eleSeedIeta_", &v_seed_ieta);
  std::vector<int> *v_seed_iphi = nullptr;
  t1->SetBranchAddress("eleSeedIphi_", &v_seed_iphi);
  std::vector<float> *v_seed_eta = nullptr;
  t1->SetBranchAddress("eleSeedEta_", &v_seed_eta);
  std::vector<float> *v_seed_phi = nullptr;
  t1->SetBranchAddress("eleSeedPhi_", &v_seed_phi);
  std::vector<int> *v_tower_ieta = nullptr;
  t1->SetBranchAddress("seedHcalIeta_", &v_tower_ieta);
  std::vector<int> *v_tower_iphi = nullptr;
  t1->SetBranchAddress("seedHcalIphi_", &v_tower_iphi);

  // the rechits are stored per electron
  // as in start from electron seed, find the hcal tower behind it
  // and then take the rechits around it within some ieta/iphi window
  // thus the vec of vec
  std::vector<std::vector<int>> *v_hcal_rawid = nullptr;
  t1->SetBranchAddress("hcalRechitRawID_", &v_hcal_rawid);
  std::vector<std::vector<int>> *v_hcal_ieta = nullptr;
  t1->SetBranchAddress("hcalRechitIeta_", &v_hcal_ieta);
  std::vector<std::vector<int>> *v_hcal_iphi = nullptr;
  t1->SetBranchAddress("hcalRechitIphi_", &v_hcal_iphi);
  std::vector<std::vector<int>> *v_hcal_depth = nullptr;
  t1->SetBranchAddress("hcalRechitDepth_", &v_hcal_depth);
  std::vector<std::vector<float>> *v_hcal_e = nullptr;
  t1->SetBranchAddress("hcalRechitEnergy_", &v_hcal_e);
  std::vector<std::vector<int>> *v_hcal_adieta = nullptr;
  t1->SetBranchAddress("hcalRechitAbsDIetaFromEleSeed_", &v_hcal_adieta);
  std::vector<std::vector<int>> *v_hcal_adiphi = nullptr;
  t1->SetBranchAddress("hcalRechitAbsDIphiFromEleSeed_", &v_hcal_adiphi);
  std::vector<std::vector<float>> *v_hcal_eta = nullptr;
  t1->SetBranchAddress("hcalRechitEta_", &v_hcal_eta);
  std::vector<std::vector<float>> *v_hcal_phi = nullptr;
  t1->SetBranchAddress("hcalRechitPhi_", &v_hcal_phi);

  float pu_true;
  t1->SetBranchAddress("puTrue_", &pu_true);
  int pu_obs;
  t1->SetBranchAddress("puObs_", &pu_obs);

  // will need these within the loop - let's not assign vecs for every event/ele
  std::vector<int> v_iPlot; // event/electron to plot allocation
  std::vector<int> v_iHR1 = {0, 6, 11, 15, 18, 20}, v_iHR2 = {0, 4, 7, 9, 10}; // how much index to skip for the Hb / Ha plots...
  std::vector<double> v_hcal_sum_e(7, 0.); // hcal energy sums by depth

  const long long reportEvery = std::pow(10, std::floor(std::log10(t1->GetEntries())) - 1);
  for (long long iT = 0; iT < t1->GetEntries(); ++iT) {
    if (iT % reportEvery == 0)
      std::cout << "Processing entry " << iT << "/" << t1->GetEntries() << std::endl;

    // this is to keep track which file we're looking at atm know which histos etc to assign later
    // relies on the files being added to the chain sequentially
    int aF = -1;
    for (int iF = v_entry.size() - 1; iF > -1; --iF) {
      if (iT < v_entry[iF])
        aF = iF;
    }
    const double weight = (aF == 0) ? v_xsec[aF] / v_entry[aF] : v_xsec[aF] / (v_entry[aF] - v_entry[aF - 1]);

    // proceed with the rest of the analysis
    t1->GetEntry(iT);
    const int nEle = v_match_rg->size();

    for (int iE = 0; iE < nEle; ++iE) {
      // cleanup - done first to ensure the vecs are always in the right state prior to filling
      v_iPlot.clear();
      v_hcal_sum_e.assign(7, 0.);

      // determine which plots will the electron be filled into
      // based on gen match and process
      if ((*v_match_rg)[iE]) {
        if (aF < 3) {
          v_iPlot.push_back(5);
          if (aF == 0)
            v_iPlot.push_back(0);
          else if (aF == 1)
            v_iPlot.push_back(1);
          else if (aF == 2)
            v_iPlot.push_back(2);
        }
        else {
          v_iPlot.push_back(6);
          v_iPlot.push_back(4);
        }
      }
      else {
        v_iPlot.push_back(6);
        if (aF == 2)
          v_iPlot.push_back(3);
        else if (aF > 2)
          v_iPlot.push_back(4);
      }

      if ((*v_pt)[iE] < 20.)
        continue;

      const float sc_e = (*v_sc_e)[iE];
      if (sc_e == 0.)
        std::cout << "ERROR: found electron passing pT cut but with 0 energy!!" << std::endl;

      const float &sc_eta = (*v_sc_eta)[iE];
      for (auto &iP : v_iPlot)
        v_hist[0][iP].plot->Fill(sc_eta, weight);

      // simplistic flag denoting the EB - EE offset
      int fEcal = -1;
      if (std::abs(sc_eta) < etaEB)
        fEcal = 0;
      else if (std::abs(sc_eta) > etaET and std::abs(sc_eta) < etaEE)
        fEcal = 7;
      else
        continue;

      // need a simple ID cut, otherwise QCD is too much
      // ... //

      const float &seed_eta = (*v_seed_eta)[iE], &seed_phi = (*v_seed_phi)[iE];

      // sum up the energies by depth
      for (int iH = 0; iH < (*v_hcal_rawid)[iE].size(); ++iH) {
        if ((*v_hcal_depth)[iE][iH] < 1 or (*v_hcal_depth)[iE][iH] > 7) {
          std::cout << "ERROR: found HCAL rechit with depth outside 1 - 7!! Recheck code logic!!!" << std::endl;
          continue;
        }

        const float dR_seed_hcal = deltaR(seed_eta, seed_phi, (*v_hcal_eta)[iE][iH], (*v_hcal_phi)[iE][iH]);
        if (dR_seed_hcal > 0.05)
          continue;

        v_hcal_sum_e[ (*v_hcal_depth)[iE][iH] - 1 ] += (*v_hcal_e)[iE][iH];
      }

      // fill histos
      for (int iV = 0; iV < 7; ++iV) {
        for (auto &iP : v_iPlot) {
          v_hist[1 + iV + fEcal][iP].plot->Fill(v_hcal_sum_e[iV], weight);
          v_hist[15 + iV + fEcal][iP].plot->Fill(v_hcal_sum_e[iV] / sc_e, weight);
        }

        if (iV == 6)
          continue;

        const double hcal_2depth_e = std::accumulate(std::next(std::begin(v_hcal_sum_e), iV), std::next(std::begin(v_hcal_sum_e), iV + 2), 0.);
        for (auto &iP : v_iPlot) {
          v_hist[37 + iV + fEcal - (1 * (fEcal > 0))][iP].plot->Fill(hcal_2depth_e, weight);
          v_hist[49 + iV + fEcal - (1 * (fEcal > 0))][iP].plot->Fill(hcal_2depth_e / sc_e, weight);
        }

        if (iV == 5)
          continue;

        const double hcal_3depth_e = std::accumulate(std::next(std::begin(v_hcal_sum_e), iV), std::next(std::begin(v_hcal_sum_e), iV + 3), 0.);
        for (auto &iP : v_iPlot) {
          v_hist[61 + iV + fEcal - (2 * (fEcal > 0))][iP].plot->Fill(hcal_3depth_e, weight);
          v_hist[71 + iV + fEcal - (2 * (fEcal > 0))][iP].plot->Fill(hcal_3depth_e / sc_e, weight);
        }

        if (iV == 4)
          continue;

        const double hcal_4depth_e = std::accumulate(std::next(std::begin(v_hcal_sum_e), iV), std::next(std::begin(v_hcal_sum_e), iV + 4), 0.);
        for (auto &iP : v_iPlot) {
          v_hist[81 + iV + fEcal - (3 * (fEcal > 0))][iP].plot->Fill(hcal_4depth_e, weight);
          v_hist[89 + iV + fEcal - (3 * (fEcal > 0))][iP].plot->Fill(hcal_4depth_e / sc_e, weight);
        }

        if (iV == 3)
          continue;

        const double hcal_5depth_e = std::accumulate(std::next(std::begin(v_hcal_sum_e), iV), std::next(std::begin(v_hcal_sum_e), iV + 5), 0.);
        for (auto &iP : v_iPlot) {
          v_hist[97 + iV + fEcal - (4 * (fEcal > 0))][iP].plot->Fill(hcal_5depth_e, weight);
          v_hist[103 + iV + fEcal - (4 * (fEcal > 0))][iP].plot->Fill(hcal_5depth_e / sc_e, weight);
        }

        if (iV == 2)
          continue;

        const double hcal_6depth_e = std::accumulate(std::next(std::begin(v_hcal_sum_e), iV), std::next(std::begin(v_hcal_sum_e), iV + 6), 0.);
        for (auto &iP : v_iPlot) {
          v_hist[109 + iV + fEcal - (5 * (fEcal > 0))][iP].plot->Fill(hcal_6depth_e, weight);
          v_hist[113 + iV + fEcal - (5 * (fEcal > 0))][iP].plot->Fill(hcal_6depth_e / sc_e, weight);
        }
      }

      const double hcal_7depth_e = std::accumulate(std::begin(v_hcal_sum_e), std::end(v_hcal_sum_e), 0.);
      if (!fEcal) {
        for (auto &iP : v_iPlot) {
          v_hist[29][iP].plot->Fill(hcal_7depth_e, weight);
          v_hist[31][iP].plot->Fill(hcal_7depth_e / sc_e, weight);
          v_hist[33][iP].plot->Fill((*v_hoe)[iE], weight);
          v_hist[35][iP].plot->Fill((*v_hoe5x5)[iE], weight);
        }
      }
      else {
        for (auto &iP : v_iPlot) {
          v_hist[30][iP].plot->Fill(hcal_7depth_e, weight);
          v_hist[32][iP].plot->Fill(hcal_7depth_e / sc_e, weight);
          v_hist[34][iP].plot->Fill((*v_hoe)[iE], weight);
          v_hist[36][iP].plot->Fill((*v_hoe5x5)[iE], weight);
        }
      }

      // fill up the H ratios
      for (int iV = 0; iV < 6; ++iV) {
        for (int iW = iV + 1; iW < 7; ++iW) {
          for (auto &iP : v_iPlot)
            v_hist[117 + (iW - (iV + 1)) + v_iHR1[iV] + ((3 * (fEcal > 0)) * fEcal)][iP].plot->Fill(v_hcal_sum_e[iW] / v_hcal_sum_e[iV], weight);

          if (iW - iV < 2 or iV > 3 or iW > 5) 
            continue;

          const double hW = v_hcal_sum_e[iW] + v_hcal_sum_e[iW + 1], hV = v_hcal_sum_e[iV] + v_hcal_sum_e[iV + 1];
          for (auto &iP : v_iPlot)
            v_hist[159 + (iW - (iV + 2)) + v_iHR2[iV] + ((3 * (fEcal > 0)) + fEcal)][iP].plot->Fill(hW / hV, weight);
        }
      }
    } // eol ele
  } // eol evt

  // ok we got the plots, draw em
  // args as in the wrapper above
  const std::string odir = "./plot/dR_0p05/";
  std_plot(v_hist[0], {}, odir + v_hist[0][0].plot->GetName(), false, 0.0001, 0.0599, "a. u.", 0., 0., "#eta^{SC}");

  // for dem eff plots
  std::vector<Plot<TGraph>> v_eff(v_hist[0].size());

  // plot the by-depth variables
  for (int iV = 0; iV < 7; ++iV) {
    // first draw the regular plots
    for (int iH : {1, 8})
      std_plot(v_hist[iH + iV], {}, odir + v_hist[iH + iV][0].plot->GetName(), true, 0.0001, 1., "a. u.", 0., 0., "H_{" + toStr(iV + 1) + "} [GeV]");

    for (int iH : {15, 22}) {
      std_plot(v_hist[iH + iV], {}, odir + v_hist[iH + iV][0].plot->GetName(), true, 0.0001, 1., "a. u.", 0., 0., "H_{" + toStr(iV + 1) + "}/E");

      v_eff.clear();
      v_eff.resize(v_hist[iH + iV].size());
      for (int iP = 0; iP < v_eff.size(); ++iP)
        v_eff[iP] = efficiency_profile(v_hist[iH + iV][iP]);

      std_plot({}, v_eff, odir + v_eff[0].plot->GetName(), false, 0.0001, 1.0999, "Efficiency", 
               0., v_eff[0].plot->GetX()[v_eff[0].plot->GetN() - 1], "H_{" + toStr(iV + 1) + "}/E");
    }

    if (iV == 6)
      continue;

    for (int iH : {37, 43}) {
      std_plot(v_hist[iH + iV], {}, odir + v_hist[iH + iV][0].plot->GetName(), true, 0.0001, 1., "a. u.", 
               0., 0., "H_{" + toStr(iV + 1) + toStr(iV + 2) + "} [GeV]");
    }

    for (int iH : {49, 55}) {
      std_plot(v_hist[iH + iV], {}, odir + v_hist[iH + iV][0].plot->GetName(), true, 0.0001, 1., "a. u.", 
               0., 0., "H_{" + toStr(iV + 1) + toStr(iV + 2) + "}/E");

      v_eff.clear();
      v_eff.resize(v_hist[iH + iV].size());
      for (int iP = 0; iP < v_eff.size(); ++iP)
        v_eff[iP] = efficiency_profile(v_hist[iH + iV][iP]);

      std_plot({}, v_eff, odir + v_eff[0].plot->GetName(), false, 0.0001, 1.0999, "Efficiency", 
               0., v_eff[0].plot->GetX()[v_eff[0].plot->GetN() - 1], "H_{" + toStr(iV + 1) + toStr(iV + 2) + "}/E");
    }

    for (int iW = iV + 1; iW < 7; ++iW) {
      for (int iH : {117, 138}) {
        std_plot(v_hist[iH + (iW - (iV + 1)) + v_iHR1[iV]], {}, odir + v_hist[iH + (iW - (iV + 1)) + v_iHR1[iV]][0].plot->GetName(), 
                 true, 0.0001, 1., "a. u.", 0., 0., "H_{" + toStr(iW + 1) + "}/H_{" + toStr(iV + 1) + "}");

        v_eff.clear();
        v_eff.resize(v_hist[iH + (iW - (iV + 1)) + v_iHR1[iV]].size());
        for (int iP = 0; iP < v_eff.size(); ++iP)
          v_eff[iP] = efficiency_profile(v_hist[iH + (iW - (iV + 1)) + v_iHR1[iV]][iP]);

        std_plot({}, v_eff, odir + v_eff[0].plot->GetName(), false, 0.0001, 1.0999, "Efficiency", 
                 0., v_eff[0].plot->GetX()[v_eff[0].plot->GetN() - 1], "H_{" + toStr(iW + 1) + "}/H_{" + toStr(iV + 1) + "}");
      }
    }

    if (iV == 5)
      continue;

    for (int iH : {61, 66}) {
      std_plot(v_hist[iH + iV], {}, odir + v_hist[iH + iV][0].plot->GetName(), true, 0.0001, 1., "a. u.", 
               0., 0., "H_{" + toStr(iV + 1) + toStr(iV + 3) + "} [GeV]");
    }

    for (int iH : {71, 76}) {
      std_plot(v_hist[iH + iV], {}, odir + v_hist[iH + iV][0].plot->GetName(), true, 0.0001, 1., "a. u.", 
               0., 0., "H_{" + toStr(iV + 1) + toStr(iV + 3) + "}/E");

      v_eff.clear();
      v_eff.resize(v_hist[iH + iV].size());
      for (int iP = 0; iP < v_eff.size(); ++iP)
        v_eff[iP] = efficiency_profile(v_hist[iH + iV][iP]);

      std_plot({}, v_eff, odir + v_eff[0].plot->GetName(), false, 0.0001, 1.0999, "Efficiency", 
               0., v_eff[0].plot->GetX()[v_eff[0].plot->GetN() - 1], "H_{" + toStr(iV + 1) + toStr(iV + 3) + "}/E");
    }

    if (iV == 4)
      continue;

    for (int iH : {81, 85}) {
      std_plot(v_hist[iH + iV], {}, odir + v_hist[iH + iV][0].plot->GetName(), true, 0.0001, 1., "a. u.", 
               0., 0., "H_{" + toStr(iV + 1) + toStr(iV + 4) + "} [GeV]");
    }

    for (int iH : {89, 93}) {
      std_plot(v_hist[iH + iV], {}, odir + v_hist[iH + iV][0].plot->GetName(), true, 0.0001, 1., "a. u.", 
               0., 0., "H_{" + toStr(iV + 1) + toStr(iV + 4) + "}/E");

      v_eff.clear();
      v_eff.resize(v_hist[iH + iV].size());
      for (int iP = 0; iP < v_eff.size(); ++iP)
        v_eff[iP] = efficiency_profile(v_hist[iH + iV][iP]);

      std_plot({}, v_eff, odir + v_eff[0].plot->GetName(), false, 0.0001, 1.0999, "Efficiency", 
               0., v_eff[0].plot->GetX()[v_eff[0].plot->GetN() - 1], "H_{" + toStr(iV + 1) + toStr(iV + 4) + "}/E");
    }

    for (int iW = iV + 2; iW < 6; ++iW) {
      const std::string sV = toStr(iV + 1) + toStr(iV + 2), sW = toStr(iW + 1) + toStr(iW + 2);
      for (int iH : {159, 169}) {
        std_plot(v_hist[iH + (iW - (iV + 2)) + v_iHR2[iV]], {}, odir + v_hist[iH + (iW - (iV + 2)) + v_iHR2[iV]][0].plot->GetName(), 
                 true, 0.0001, 1., "a. u.", 0., 0., "H_{" + sW + "}/H_{" + sV + "}");

        v_eff.clear();
        v_eff.resize(v_hist[iH + (iW - (iV + 2)) + v_iHR2[iV]].size());
        for (int iP = 0; iP < v_eff.size(); ++iP)
          v_eff[iP] = efficiency_profile(v_hist[iH + (iW - (iV + 2)) + v_iHR2[iV]][iP]);

        std_plot({}, v_eff, odir + v_eff[0].plot->GetName(), false, 0.0001, 1.0999, "Efficiency", 
                 0., v_eff[0].plot->GetX()[v_eff[0].plot->GetN() - 1], "H_{" + sW + "}/H_{" + sV + "}");
      }
    }

    if (iV == 3)
      continue;

    for (int iH : {97, 100}) {
      std_plot(v_hist[iH + iV], {}, odir + v_hist[iH + iV][0].plot->GetName(), true, 0.0001, 1., "a. u.", 
               0., 0., "H_{" + toStr(iV + 1) + toStr(iV + 5) + "} [GeV]");
    }

    for (int iH : {103, 106}) {
      std_plot(v_hist[iH + iV], {}, odir + v_hist[iH + iV][0].plot->GetName(), true, 0.0001, 1., "a. u.", 
               0., 0., "H_{" + toStr(iV + 1) + toStr(iV + 5) + "}/E");

      v_eff.clear();
      v_eff.resize(v_hist[iH + iV].size());
      for (int iP = 0; iP < v_eff.size(); ++iP)
        v_eff[iP] = efficiency_profile(v_hist[iH + iV][iP]);

      std_plot({}, v_eff, odir + v_eff[0].plot->GetName(), false, 0.0001, 1.0999, "Efficiency", 
               0., v_eff[0].plot->GetX()[v_eff[0].plot->GetN() - 1], "H_{" + toStr(iV + 1) + toStr(iV + 5) + "}/E");
    }

    if (iV == 2)
      continue;

    for (int iH : {109, 111}) {
      std_plot(v_hist[iH + iV], {}, odir + v_hist[iH + iV][0].plot->GetName(), true, 0.0001, 1., "a. u.", 
               0., 0., "H_{" + toStr(iV + 1) + toStr(iV + 6) + "} [GeV]");
    }

    for (int iH : {113, 115}) {
      std_plot(v_hist[iH + iV], {}, odir + v_hist[iH + iV][0].plot->GetName(), true, 0.0001, 1., "a. u.", 
               0., 0., "H_{" + toStr(iV + 1) + toStr(iV + 6) + "}/E");

      v_eff.clear();
      v_eff.resize(v_hist[iH + iV].size());
      for (int iP = 0; iP < v_eff.size(); ++iP)
        v_eff[iP] = efficiency_profile(v_hist[iH + iV][iP]);

      std_plot({}, v_eff, odir + v_eff[0].plot->GetName(), false, 0.0001, 1.0999, "Efficiency", 
               0., v_eff[0].plot->GetX()[v_eff[0].plot->GetN() - 1], "H_{" + toStr(iV + 1) + toStr(iV + 6) + "}/E");
    }
  }

  // plot the sum of all depth variables
  for (int iH : {29, 30})
    std_plot(v_hist[iH], {}, odir + v_hist[iH][0].plot->GetName(), true, 0.0001, 1., "a. u.", 0., 0., "H_{17} [GeV]");

  for (int iH : {31, 32}) {
    std_plot(v_hist[iH], {}, odir + v_hist[iH][0].plot->GetName(), true, 0.0001, 1., "a. u.", 0., 0., "H_{17}/E");

    v_eff.clear();
    v_eff.resize(v_hist[iH].size());
    for (int iP = 0; iP < v_eff.size(); ++iP)
      v_eff[iP] = efficiency_profile(v_hist[iH][iP]);

    std_plot({}, v_eff, odir + v_eff[0].plot->GetName(), false, 0.0001, 1.0999, "Efficiency", 
             0., v_eff[0].plot->GetX()[v_eff[0].plot->GetN() - 1], "H_{17}/E");
  }

  // now dem cmssw vars
  for (int iH : {33, 34}) {
    std_plot(v_hist[iH], {}, odir + v_hist[iH][0].plot->GetName(), true, 0.0001, 1., "a. u.", 0., 0., "CMSSW H/E");

    v_eff.clear();
    v_eff.resize(v_hist[iH].size());
    for (int iP = 0; iP < v_eff.size(); ++iP)
      v_eff[iP] = efficiency_profile(v_hist[iH][iP]);

    std_plot({}, v_eff, odir + v_eff[0].plot->GetName(), false, 0.0001, 1.0999, "Efficiency", 
             0., v_eff[0].plot->GetX()[v_eff[0].plot->GetN() - 1], "CMSSW H/E");
  }

  for (int iH : {35, 36}) {
    std_plot(v_hist[iH], {}, odir + v_hist[iH][0].plot->GetName(), true, 0.0001, 1., "a. u.", 0., 0., "CMSSW 5#times5 H/E");

    v_eff.clear();
    v_eff.resize(v_hist[iH].size());
    for (int iP = 0; iP < v_eff.size(); ++iP)
      v_eff[iP] = efficiency_profile(v_hist[iH][iP]);

    std_plot({}, v_eff, odir + v_eff[0].plot->GetName(), false, 0.0001, 1.0999, "Efficiency", 
             0., v_eff[0].plot->GetX()[v_eff[0].plot->GetN() - 1], "CMSSW 5#times5 H/E");
  }

  // legend, if desired in its own canvas
  legend_canvas(v_hist[0], {}, "./plot/legend", // v_hist, v_graph, filename (w/o ext)
                2, 0, 0, 42, 0.079, "", // legend style, header
                0.05, 1., 0.15, 1. // legend position
                );
  std::cout << std::endl;

  // make the roc, but no need to draw
  const std::vector<int> v_var = {15, 16, 17, 18, 19, 20, 21, 49, 50, 51, 52, 53, 54, 71, 72, 73, 74, 75, 
                                  89, 90, 91, 92, 103, 104, 105, 113, 114, 31,
                                  117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129,
                                  130, 131, 132, 133, 134, 135, 136, 137,
                                  159, 160, 161, 162, 163, 164, 165, 166, 167, 168,
                                  33, 35,

                                  22, 23, 24, 25, 26, 27, 28, 55, 56, 57, 58, 59, 60, 76, 77, 78, 79, 80, 
                                  93, 94, 95, 96, 106, 107, 108, 115, 116, 32,
                                  138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 
                                  151, 152, 153, 154, 155, 156, 157, 158,
                                  169, 170, 171, 172, 173, 174, 175, 176, 177, 178,
                                  34, 36};

  std::vector<Plot<TGraph>> v_roc(v_var.size());
  v_roc.clear();
  v_roc.resize(v_var.size());
  for (int iV = 0; iV < v_var.size(); ++iV)
    v_roc[iV] = discrimination_profile(v_hist[ v_var[iV] ], 5, 6, kBlack, v_hist[ v_var[iV] ][0].plot->GetName(), true);

  gROOT->ProcessLine(".q");
}

/*** old 200122 dirs
  const std::vector<std::string> v_file = {"DYToEE_M-50_NNPDF31_TuneCP5_14TeV-powheg-pythia8/crab_dy_ee/200122_165822/0000/", 
                                           "GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_14TeV_Pythia8/crab_yjet_pt_40_myy_80/200122_165901/0000/",
                                           "TTToSemiLeptonic_TuneCP5_14TeV-powheg-pythia8/crab_tt_lj/200122_165737/0000/",
                                           "QCD_Pt_15to30_TuneCP5_14TeV_pythia8/crab_qcd_pt_15_30/200122_165750/0000/", 
                                           "QCD_Pt_30to50_TuneCP5_14TeV_pythia8/crab_qcd_pt_30_50/200122_165850/0000/",
                                           "QCD_Pt_50to80_TuneCP5_14TeV_pythia8/crab_qcd_pt_50_80/200122_165731/0000/",
                                           "QCD_Pt_80to120_TuneCP5_14TeV_pythia8/crab_qcd_pt_80_120/200122_165833/0000/", 
                                           "QCD_Pt_120to170_TuneCP5_14TeV_pythia8/crab_qcd_pt_120_170/200122_165839/0000/",
                                           "QCD_Pt_170to300_TuneCP5_14TeV_pythia8/crab_qcd_pt_170_300/200122_165906/0000/",
                                           "QCD_Pt_300to470_TuneCP5_14TeV_pythia8/crab_qcd_pt_300_470/200122_165757/0000/", 
                                           "QCD_Pt_470to600_TuneCP5_14TeV_pythia8/crab_qcd_pt_470_600/200122_165804/0000/",
                                           "QCD_Pt_600to800_TuneCP5_14TeV_pythia8/crab_qcd_pt_600_800/200122_165844/0000/",
                                           "QCD_Pt_800to1000_TuneCP5_14TeV_pythia8/crab_qcd_pt_800_1000/200122_165743/0000/", 
                                           "QCD_Pt_1000to1400_TuneCP5_14TeV_pythia8/crab_qcd_pt_1000_1400/200122_165855/0000/",
                                           "QCD_Pt_1400to1800_TuneCP5_14TeV_pythia8/crab_qcd_pt_1400_1800/200122_165811/0000/",
                                           "QCD_Pt_1800to2400_TuneCP5_14TeV_pythia8/crab_qcd_pt_1800_2400/200122_165816/0000/", 
                                           "QCD_Pt_2400to3200_TuneCP5_14TeV_pythia8/crab_qcd_pt_2400_3200/200122_165827/0000/",
                                           "QCD_Pt_3200toInf_TuneCP5_14TeV_pythia8/crab_qcd_pt_3200/200122_165912/0000/"};
***/




