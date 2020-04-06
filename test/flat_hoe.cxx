// -*- C++ -*-
// root -l -b flat_hoe.cxx++

//#include "/home/afiqaize/Downloads/eft/EFTFitter/src/PlotterUtil.h"
#include "/afs/cern.ch/work/a/afiqaize/EFTFitter/src/PlotterUtil.h"
#include "TChain.h"
#include "TH3.h"
#include "TH3F.h"

constexpr float etaEB = 1.4442f, etaET = 1.566f, etaEE = 2.5f;

template <typename K, typename V>
V& ref_to(std::vector<std::pair<K, V>> &vec, const K &key)
{
  for (auto &p : vec) {
    if (key == p.first)
      return p.second;
  }

  return vec[-1].second;
}

float accumulate(const std::vector<float> &vec, int begin, int end) {
  float val = 0.;
  for (int it = begin; it < end; ++it)
    val += vec[it];

  return val;
}

void assign_plot3(const std::string &var, 
                  const int nBin1, const float min1, const float max1,
                  const int nBin2, const float min2, const float max2,
                  const int nBin3, const float min3, const float max3,
                  const std::vector<std::pair<int, std::array<std::string, 4>>> &vp_plot,
                  std::vector<std::pair<std::string, std::vector<Plot<TH3>>>> &v_hist) {
  v_hist.emplace_back( var, std::vector<Plot<TH3>>(vp_plot.size()) );

  for (int iP = 0; iP < vp_plot.size(); ++iP) {
    ref_to(v_hist, var)[iP].plot = std::make_unique<TH3F>((var + vp_plot[iP].second[0]).c_str(), "", 
                                                          nBin1, min1, max1, nBin2, min2, max2, nBin3, min3, max3);

    stylePlot(ref_to(v_hist, var)[iP].plot.get(), vp_plot[iP].first, 1., 0, 20 + iP, 1.5, 1, 3);

    ref_to(v_hist, var)[iP].legend_txt = vp_plot[iP].second[1];
    ref_to(v_hist, var)[iP].legend_opt = vp_plot[iP].second[2];
    ref_to(v_hist, var)[iP].draw_opt = vp_plot[iP].second[3];
  }
}

void flat_hoe() {
  gROOT->Reset();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2(true);
  TGaxis::SetMaxDigits(3);

  const std::string idir = "/eos/cms/store/group/phys_egamma/afiqaize/hcal_depth_200310/";
  const std::vector<std::string> v_file = {"DYToEE_M-50_NNPDF31_TuneCP5_14TeV-powheg-pythia8/crab_dy_ee/200310_173609/0000/",
                                           "GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_14TeV_Pythia8/crab_yjet_pt_40_myy_80/200310_173652/0000/",
                                           "TTToSemiLeptonic_TuneCP5_14TeV-powheg-pythia8/crab_tt_lj/200310_173528/0000/",
                                           "QCD_Pt_15to30_TuneCP5_14TeV_pythia8/crab_qcd_pt_15_30/200310_173540/0000/",
                                           "QCD_Pt_30to50_TuneCP5_14TeV_pythia8/crab_qcd_pt_30_50/200310_173640/0000/",
                                           "QCD_Pt_50to80_TuneCP5_14TeV_pythia8/crab_qcd_pt_50_80/200310_173521/0000/",
                                           "QCD_Pt_80to120_TuneCP5_14TeV_pythia8/crab_qcd_pt_80_120/200310_173622/0000/",
                                           "QCD_Pt_120to170_TuneCP5_14TeV_pythia8/crab_qcd_pt_120_170/200310_173628/0000/",
                                           "QCD_Pt_170to300_TuneCP5_14TeV_pythia8/crab_qcd_pt_170_300/200310_173703/0000/",
                                           "QCD_Pt_300to470_TuneCP5_14TeV_pythia8/crab_qcd_pt_300_470/200310_173546/0000/",
                                           "QCD_Pt_470to600_TuneCP5_14TeV_pythia8/crab_qcd_pt_470_600/200310_173552/0000/",
                                           "QCD_Pt_600to800_TuneCP5_14TeV_pythia8/crab_qcd_pt_600_800/200310_173633/0000/",
                                           "QCD_Pt_800to1000_TuneCP5_14TeV_pythia8/crab_qcd_pt_800_1000/200310_173534/0000/",
                                           "QCD_Pt_1000to1400_TuneCP5_14TeV_pythia8/crab_qcd_pt_1000_1400/200310_173646/0000/",
                                           "QCD_Pt_1400to1800_TuneCP5_14TeV_pythia8/crab_qcd_pt_1400_1800/200310_173558/0000/",
                                           "QCD_Pt_1800to2400_TuneCP5_14TeV_pythia8/crab_qcd_pt_1800_2400/200310_173603/0000/",
                                           "QCD_Pt_2400to3200_TuneCP5_14TeV_pythia8/crab_qcd_pt_2400_3200/200310_173615/0000/",
                                           "QCD_Pt_3200toInf_TuneCP5_14TeV_pythia8/crab_qcd_pt_3200/200310_173710/0000/",
                                           "/afs/cern.ch/work/a/afiqaize/egm/hcal_depth_200122/CMSSW_11_0_1/src/HoE/HoEAnalyzer/test/"};

  // ensure same order as above
  // weight = xsec / nevt, xsec just taken off xsecdb, nevt read off files
  // 13 TeV xsec since some aren't available for 14 TeV, but fine
  const std::vector<float> v_xsec = {2137.f,
                                     874.2f,
                                     831.76f * 4.f / 9.f,
                                     0., //1.248e9f,
                                     0., //1.066e8f,
                                     1.568e7f,
                                     2.336e6f,
                                     4.068e5f,
                                     1.033e5f,
                                     6826.f,
                                     552.6f,
                                     156.6f,
                                     26.32f,
                                     7.5f,
                                     0.6479f,
                                     0.08715f,
                                     0.005242f,
                                     1.349e-4f,
                                     1.};

  // mark in the set of files to go into a single plot
  // key is color, value is sample key, legend text, legend type, draw option
  const std::vector<std::pair<int, std::array<std::string, 4>>> vp_plot = { {kPink + 9, {"_dy", "DY", "lp", "hist e2"}},
                                                                            {kPink + 8, {"_yjet", "#gamma + j", "lp", "hist e2"}},
                                                                            {kPink + 1, {"_ttlj_sig", "t#bar{t} (s)", "lp", "hist e2"}},
                                                                            {kAzure + 8, {"_ttlj_bkg", "t#bar{t} (b)", "lp", "hist e2"}},
                                                                            {kAzure + 4, {"_qcd", "QCD", "lp", "hist e2"}},
                                                                            {kPink - 1, {"_sig", "Signal", "lp", "hist e2"}},
                                                                            {kAzure - 1, {"_bkg", "Background", "lp", "hist e2"}},
                                                                            {kBlack, {"_2ele_pu0", "2e PU0", "lp", "hist e2"}}};
  std::vector<std::pair<std::string, std::vector<Plot<TH3>>>> v_h3;

  const int n_pt = 40, n_rho = 10;
  const float min_pt = 0.f, min_rho = 20.f;
  const float max_pt = 200.f, max_rho = 80.f;

  for (int iV = 0; iV < 7; ++iV) {
    assign_plot3("hse1_" + toStr(iV + 1) + toStr(iV + 1) + "_1p30", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hse1_" + toStr(iV + 1) + toStr(iV + 1) + "_1p48", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hse1_" + toStr(iV + 1) + toStr(iV + 1) + "_2p10", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hse1_" + toStr(iV + 1) + toStr(iV + 1) + "_2p50", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);

    assign_plot3("hoe1_" + toStr(iV + 1) + toStr(iV + 1) + "_1p30", 50, 0.f, 2.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hoe1_" + toStr(iV + 1) + toStr(iV + 1) + "_1p48", 50, 0.f, 2.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hoe1_" + toStr(iV + 1) + toStr(iV + 1) + "_2p10", 50, 0.f, 4.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hoe1_" + toStr(iV + 1) + toStr(iV + 1) + "_2p50", 50, 0.f, 4.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);

    assign_plot3("nxx1_" + toStr(iV + 1) + toStr(iV + 1) + "_1p30", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("nxx1_" + toStr(iV + 1) + toStr(iV + 1) + "_1p48", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("nxx1_" + toStr(iV + 1) + toStr(iV + 1) + "_2p10", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("nxx1_" + toStr(iV + 1) + toStr(iV + 1) + "_2p50", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);

    if (iV == 6) continue;

    assign_plot3("hse2_" + toStr(iV + 1) + toStr(iV + 2) + "_1p30", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hse2_" + toStr(iV + 1) + toStr(iV + 2) + "_1p48", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hse2_" + toStr(iV + 1) + toStr(iV + 2) + "_2p10", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hse2_" + toStr(iV + 1) + toStr(iV + 2) + "_2p50", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);

    assign_plot3("hoe2_" + toStr(iV + 1) + toStr(iV + 2) + "_1p30", 50, 0.f, 2.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hoe2_" + toStr(iV + 1) + toStr(iV + 2) + "_1p48", 50, 0.f, 2.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hoe2_" + toStr(iV + 1) + toStr(iV + 2) + "_2p10", 50, 0.f, 4.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hoe2_" + toStr(iV + 1) + toStr(iV + 2) + "_2p50", 50, 0.f, 4.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);

    assign_plot3("nxx2_" + toStr(iV + 1) + toStr(iV + 2) + "_1p30", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("nxx2_" + toStr(iV + 1) + toStr(iV + 2) + "_1p48", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("nxx2_" + toStr(iV + 1) + toStr(iV + 2) + "_2p10", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("nxx2_" + toStr(iV + 1) + toStr(iV + 2) + "_2p50", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);

    if (iV == 5) continue;

    assign_plot3("hse3_" + toStr(iV + 1) + toStr(iV + 3) + "_1p30", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hse3_" + toStr(iV + 1) + toStr(iV + 3) + "_1p48", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hse3_" + toStr(iV + 1) + toStr(iV + 3) + "_2p10", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hse3_" + toStr(iV + 1) + toStr(iV + 3) + "_2p50", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);

    assign_plot3("hoe3_" + toStr(iV + 1) + toStr(iV + 3) + "_1p30", 50, 0.f, 2.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hoe3_" + toStr(iV + 1) + toStr(iV + 3) + "_1p48", 50, 0.f, 2.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hoe3_" + toStr(iV + 1) + toStr(iV + 3) + "_2p10", 50, 0.f, 4.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hoe3_" + toStr(iV + 1) + toStr(iV + 3) + "_2p50", 50, 0.f, 4.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);

    assign_plot3("nxx3_" + toStr(iV + 1) + toStr(iV + 3) + "_1p30", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("nxx3_" + toStr(iV + 1) + toStr(iV + 3) + "_1p48", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("nxx3_" + toStr(iV + 1) + toStr(iV + 3) + "_2p10", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("nxx3_" + toStr(iV + 1) + toStr(iV + 3) + "_2p50", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);

    if (iV == 4) continue;

    assign_plot3("hse4_" + toStr(iV + 1) + toStr(iV + 4) + "_1p30", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hse4_" + toStr(iV + 1) + toStr(iV + 4) + "_1p48", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hse4_" + toStr(iV + 1) + toStr(iV + 4) + "_2p10", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hse4_" + toStr(iV + 1) + toStr(iV + 4) + "_2p50", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);

    assign_plot3("hoe4_" + toStr(iV + 1) + toStr(iV + 4) + "_1p30", 50, 0.f, 2.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hoe4_" + toStr(iV + 1) + toStr(iV + 4) + "_1p48", 50, 0.f, 2.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hoe4_" + toStr(iV + 1) + toStr(iV + 4) + "_2p10", 50, 0.f, 4.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hoe4_" + toStr(iV + 1) + toStr(iV + 4) + "_2p50", 50, 0.f, 4.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);

    assign_plot3("nxx4_" + toStr(iV + 1) + toStr(iV + 4) + "_1p30", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("nxx4_" + toStr(iV + 1) + toStr(iV + 4) + "_1p48", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("nxx4_" + toStr(iV + 1) + toStr(iV + 4) + "_2p10", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("nxx4_" + toStr(iV + 1) + toStr(iV + 4) + "_2p50", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);

    if (iV == 3) continue;

    assign_plot3("hse5_" + toStr(iV + 1) + toStr(iV + 5) + "_1p30", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hse5_" + toStr(iV + 1) + toStr(iV + 5) + "_1p48", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hse5_" + toStr(iV + 1) + toStr(iV + 5) + "_2p10", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hse5_" + toStr(iV + 1) + toStr(iV + 5) + "_2p50", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);

    assign_plot3("hoe5_" + toStr(iV + 1) + toStr(iV + 5) + "_1p30", 50, 0.f, 2.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hoe5_" + toStr(iV + 1) + toStr(iV + 5) + "_1p48", 50, 0.f, 2.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hoe5_" + toStr(iV + 1) + toStr(iV + 5) + "_2p10", 50, 0.f, 4.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hoe5_" + toStr(iV + 1) + toStr(iV + 5) + "_2p50", 50, 0.f, 4.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);

    assign_plot3("nxx5_" + toStr(iV + 1) + toStr(iV + 5) + "_1p30", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("nxx5_" + toStr(iV + 1) + toStr(iV + 5) + "_1p48", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("nxx5_" + toStr(iV + 1) + toStr(iV + 5) + "_2p10", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("nxx5_" + toStr(iV + 1) + toStr(iV + 5) + "_2p50", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);

    if (iV == 2) continue;

    assign_plot3("hse6_" + toStr(iV + 1) + toStr(iV + 6) + "_1p30", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hse6_" + toStr(iV + 1) + toStr(iV + 6) + "_1p48", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hse6_" + toStr(iV + 1) + toStr(iV + 6) + "_2p10", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hse6_" + toStr(iV + 1) + toStr(iV + 6) + "_2p50", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);

    assign_plot3("hoe6_" + toStr(iV + 1) + toStr(iV + 6) + "_1p30", 50, 0.f, 2.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hoe6_" + toStr(iV + 1) + toStr(iV + 6) + "_1p48", 50, 0.f, 2.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hoe6_" + toStr(iV + 1) + toStr(iV + 6) + "_2p10", 50, 0.f, 4.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hoe6_" + toStr(iV + 1) + toStr(iV + 6) + "_2p50", 50, 0.f, 4.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);

    assign_plot3("nxx6_" + toStr(iV + 1) + toStr(iV + 6) + "_1p30", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("nxx6_" + toStr(iV + 1) + toStr(iV + 6) + "_1p48", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("nxx6_" + toStr(iV + 1) + toStr(iV + 6) + "_2p10", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("nxx6_" + toStr(iV + 1) + toStr(iV + 6) + "_2p50", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);

    if (iV == 1) continue;

    assign_plot3("hse7_" + toStr(iV + 1) + toStr(iV + 7) + "_1p30", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hse7_" + toStr(iV + 1) + toStr(iV + 7) + "_1p48", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hse7_" + toStr(iV + 1) + toStr(iV + 7) + "_2p10", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hse7_" + toStr(iV + 1) + toStr(iV + 7) + "_2p50", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);

    assign_plot3("hoe7_" + toStr(iV + 1) + toStr(iV + 7) + "_1p30", 50, 0.f, 2.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hoe7_" + toStr(iV + 1) + toStr(iV + 7) + "_1p48", 50, 0.f, 2.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hoe7_" + toStr(iV + 1) + toStr(iV + 7) + "_2p10", 50, 0.f, 4.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("hoe7_" + toStr(iV + 1) + toStr(iV + 7) + "_2p50", 50, 0.f, 4.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);

    assign_plot3("nxx7_" + toStr(iV + 1) + toStr(iV + 7) + "_1p30", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("nxx7_" + toStr(iV + 1) + toStr(iV + 7) + "_1p48", 50, 0.f, 100.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("nxx7_" + toStr(iV + 1) + toStr(iV + 7) + "_2p10", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("nxx7_" + toStr(iV + 1) + toStr(iV + 7) + "_2p50", 50, 0.f, 400.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);

    assign_plot3("cmssw_hoe_1p30", 50, 0.f, 4.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("cmssw_hoe_1p48", 50, 0.f, 4.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("cmssw_hoe_2p10", 50, 0.f, 4.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("cmssw_hoe_2p50", 50, 0.f, 4.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);

    assign_plot3("cmssw_hoe_tow_1p30", 50, 0.f, 4.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("cmssw_hoe_tow_1p48", 50, 0.f, 4.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("cmssw_hoe_tow_2p10", 50, 0.f, 4.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("cmssw_hoe_tow_2p50", 50, 0.f, 4.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);

    assign_plot3("cmssw_hoe_5x5_1p30", 50, 0.f, 4.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("cmssw_hoe_5x5_1p48", 50, 0.f, 4.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("cmssw_hoe_5x5_2p10", 50, 0.f, 4.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
    assign_plot3("cmssw_hoe_5x5_2p50", 50, 0.f, 4.f, n_pt, min_pt, max_pt, n_rho, min_rho, max_rho, vp_plot, v_h3);
  }

  // put the files in chain and remember the cumulative entry count
  std::vector<long long> v_entry;
  auto t1 = std::make_unique<TChain>("tree");
  for (int iF = 0; iF < v_file.size(); ++iF) {
    if (iF != v_file.size() - 1)
      t1->Add((idir + v_file[iF] + "flat_electron_*.root").c_str());
    else
      t1->Add((v_file[iF] + "flat_electron.root").c_str());

    v_entry.push_back(t1->GetEntries());
  }
  //t1->SetBranchStatus("*", 0);

  int n_ele;
  t1->SetBranchAddress("n_ele", &n_ele);

  // could do with the realloc trick but it's a pain in small macros
  // just grep the highest size over all files and use that
  const int max_ele = 16;
  int  ele_gap_eb_ee[max_ele]; 
  t1->SetBranchAddress("ele_gap_eb_ee", ele_gap_eb_ee);
  int  ele_gap_eb_eta[max_ele]; 
  t1->SetBranchAddress("ele_gap_eb_eta", ele_gap_eb_eta);
  int  ele_gap_eb_phi[max_ele]; 
  t1->SetBranchAddress("ele_gap_eb_phi", ele_gap_eb_phi);
  int  ele_gap_ee_dee[max_ele]; 
  t1->SetBranchAddress("ele_gap_ee_dee", ele_gap_ee_dee);
  int  ele_gap_ee_ring[max_ele]; 
  t1->SetBranchAddress("ele_gap_ee_ring", ele_gap_ee_ring);

  int  ele_genmatch[max_ele]; 
  t1->SetBranchAddress("ele_genmatch", ele_genmatch);
  float  ele_dR_reco_gen[max_ele]; 
  t1->SetBranchAddress("ele_dR_reco_gen", ele_dR_reco_gen);
  float  ele_pt_ratio_reco_gen[max_ele]; 
  t1->SetBranchAddress("ele_pt_ratio_reco_gen", ele_pt_ratio_reco_gen);

  float  ele_sc_energy[max_ele]; 
  t1->SetBranchAddress("ele_sc_energy", ele_sc_energy);
  float  ele_sc_raw_energy[max_ele]; 
  t1->SetBranchAddress("ele_sc_raw_energy", ele_sc_raw_energy);

  float  ele_cmssw_hoe[max_ele]; 
  t1->SetBranchAddress("ele_cmssw_hoe", ele_cmssw_hoe);
  float  ele_cmssw_hoe_tower[max_ele]; 
  t1->SetBranchAddress("ele_cmssw_hoe_tower", ele_cmssw_hoe_tower);
  float  ele_cmssw_hoe_5x5[max_ele]; 
  t1->SetBranchAddress("ele_cmssw_hoe_5x5", ele_cmssw_hoe_5x5);

  float  ele_sieie_5x5[max_ele]; 
  t1->SetBranchAddress("ele_sieie_5x5", ele_sieie_5x5);
  float  ele_r9_5x5[max_ele]; 
  t1->SetBranchAddress("ele_r9_5x5", ele_r9_5x5);

  float  ele_sc_eta[max_ele]; 
  t1->SetBranchAddress("ele_sc_eta", ele_sc_eta);
  float  ele_pt[max_ele]; 
  t1->SetBranchAddress("ele_pt", ele_pt);
  float  ele_phi[max_ele]; 
  t1->SetBranchAddress("ele_phi", ele_phi);

  int    ele_seed_ieta[max_ele]; 
  t1->SetBranchAddress("ele_seed_ieta", ele_seed_ieta);
  int    ele_seed_iphi[max_ele];
  t1->SetBranchAddress("ele_seed_iphi", ele_seed_iphi);
  float  ele_seed_eta[max_ele]; 
  t1->SetBranchAddress("ele_seed_eta", ele_seed_eta);
  float  ele_seed_phi[max_ele];
  t1->SetBranchAddress("ele_seed_phi", ele_seed_phi);
  int    ele_seed_hcal_ieta[max_ele]; 
  t1->SetBranchAddress("ele_seed_hcal_ieta", ele_seed_hcal_ieta);
  int    ele_seed_hcal_iphi[max_ele];
  t1->SetBranchAddress("ele_seed_hcal_iphi", ele_seed_hcal_iphi);

  int n_hcalhit;
  t1->SetBranchAddress("n_hcalhit", &n_hcalhit);

  const int max_hcalhit = 1024;
  float  hcalhit_energy[max_hcalhit]; 
  t1->SetBranchAddress("hcalhit_energy", hcalhit_energy);
  int    hcalhit_depth[max_hcalhit]; 
  t1->SetBranchAddress("hcalhit_depth", hcalhit_depth);

  int    hcalhit_ieta[max_hcalhit]; 
  t1->SetBranchAddress("hcalhit_ieta", hcalhit_ieta);
  int    hcalhit_iphi[max_hcalhit]; 
  t1->SetBranchAddress("hcalhit_iphi", hcalhit_iphi);

  float  hcalhit_eta[max_hcalhit]; 
  t1->SetBranchAddress("hcalhit_eta", hcalhit_eta);
  float  hcalhit_phi[max_hcalhit]; 
  t1->SetBranchAddress("hcalhit_phi", hcalhit_phi);

  int    hcalhit_seed_dieta[max_hcalhit]; 
  t1->SetBranchAddress("hcalhit_seed_dieta", hcalhit_seed_dieta);
  int    hcalhit_seed_diphi[max_hcalhit]; 
  t1->SetBranchAddress("hcalhit_seed_diphi", hcalhit_seed_diphi);

  int    hcalhit_ele_index[max_hcalhit]; 
  t1->SetBranchAddress("hcalhit_ele_index", hcalhit_ele_index);

  float rho;
  t1->SetBranchAddress("rho", &rho);

  // will need these within the loop - let's not assign vecs for every event/ele
  std::vector<int> v_iPlot; // event/electron to plot allocation
  std::vector<float> v_hse(7, 0.f), v_hoe(7, 0.f), v_nxx(7, 0.f); // hcal energy, h/e and noise sums by depth
  //const float tower_area = 2.5f * 0.0175f;

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
    const float weight = (aF == 0) ? v_xsec[aF] / v_entry[aF] : v_xsec[aF] / (v_entry[aF] - v_entry[aF - 1]);

    // proceed with the rest of the analysis
    t1->GetEntry(iT);

    // set rho of PU0 sample to a fake value
    if (aF == v_entry.size() - 1)
      rho = 30.f;

    for (int iE = 0; iE < n_ele; ++iE) {
      // cleanup - done first to ensure the vecs are always in the right state prior to filling
      v_iPlot.clear();
      v_hse.assign(7, 0.);
      v_hoe.assign(7, 0.);
      v_nxx.assign(7, 0.);

      const bool ele_gap = ele_gap_eb_ee[iE] or ele_gap_eb_eta[iE] or ele_gap_eb_phi[iE] or ele_gap_ee_dee[iE] or ele_gap_ee_ring[iE];
      const bool ele_good = ele_genmatch[iE] and !ele_gap;

      // determine which plots will the electron be filled into
      // based on gen match and process
      if (ele_good) {
        if (aF == v_entry.size() - 1)
          v_iPlot.emplace_back(7);
        else if (aF < 3) {
          v_iPlot.emplace_back(5);
          if (aF == 0)
            v_iPlot.emplace_back(0);
          else if (aF == 1)
            v_iPlot.emplace_back(1);
          else if (aF == 2)
            v_iPlot.emplace_back(2);
        }
        else {
          v_iPlot.emplace_back(6);
          v_iPlot.emplace_back(4);
        }
      }
      else {
        if (aF != v_entry.size() - 1) {
          v_iPlot.emplace_back(6);
          if (aF == 2)
            v_iPlot.emplace_back(3);
          else if (aF > 2)
            v_iPlot.emplace_back(4);
        }
      }

      const float &sc_eta = ele_sc_eta[iE];

      std::string fEta = "_xxx";
      if (std::abs(sc_eta) < 1.3)
        fEta = "_1p30";
      else if (std::abs(sc_eta) < etaEB)
        fEta = "_1p48";
      else if (std::abs(sc_eta) > etaET and std::abs(sc_eta) < 2.1)
        fEta = "_2p10";
      else if (std::abs(sc_eta) > 2.1 and std::abs(sc_eta) < etaEE)
        fEta = "_2p50";
      else
        continue;

      // sum up the energies by depth, if the hit is associated to the electron
      // consider only those within dieta, diphi < 5 that is indexed to the ele
      // dieta and diphi cut might be smaller, read in!
      // hits that are within dieta < 5, 13 < diphi < 23 are tagged as noise   
      for (int iH = 0; iH < n_hcalhit; ++iH) {
        if (hcalhit_depth[iH] < 1 or hcalhit_depth[iH] > 7) {
          std::cout << "ERROR: found HCAL rechit with depth " << hcalhit_depth[iH] << ", outside 1 - 7!! Recheck code logic!!!" << std::endl;
          continue;
        }

        if (hcalhit_ele_index[iH] != iE)
          continue;

        if (std::abs(hcalhit_seed_dieta[iH]) < 5) {
          // regular nxn square arrays
          // ele_ieta_iphi_11x11 --- std::abs(hcalhit_seed_dieta[iH]) < 6 and std::abs(hcalhit_seed_diphi[iH]) < 6
          // ele_ieta_iphi_9x9 --- std::abs(hcalhit_seed_dieta[iH]) < 5 and std::abs(hcalhit_seed_diphi[iH]) < 5
          // ele_ieta_iphi_7x7 --- std::abs(hcalhit_seed_dieta[iH]) < 4 and std::abs(hcalhit_seed_diphi[iH]) < 4
          // ele_ieta_iphi_5x5 --- std::abs(hcalhit_seed_dieta[iH]) < 3 and std::abs(hcalhit_seed_diphi[iH]) < 3
          // ele_ieta_iphi_3x3 --- std::abs(hcalhit_seed_dieta[iH]) < 2 and std::abs(hcalhit_seed_diphi[iH]) < 2

          // these are just the square arrays but with an extra cell when ieta xor iphi == 0 (4 extra cells in total)
          // ele_dieta_plus_diphi_1 --- std::abs(hcalhit_seed_dieta[iH]) + std::abs(hcalhit_seed_diphi[iH]) < 2
          // ele_dieta_plus_diphi_2 --- std::abs(hcalhit_seed_dieta[iH]) + std::abs(hcalhit_seed_diphi[iH]) < 3
          // ele_dieta_plus_diphi_3 --- std::abs(hcalhit_seed_dieta[iH]) + std::abs(hcalhit_seed_diphi[iH]) < 4

          // cone based definitions
          // ele_dR_0p10 --- deltaR(hcalhit_eta[iH], hcalhit_phi[iH], ele_seed_eta[iE], ele_seed_phi[iE]) < 0.1f

          if (deltaR(hcalhit_eta[iH], hcalhit_phi[iH], ele_seed_eta[iE], ele_seed_phi[iE]) < 0.15f) {
            v_hse[ hcalhit_depth[iH] - 1 ] += hcalhit_energy[iH];
            v_hoe[ hcalhit_depth[iH] - 1 ] += hcalhit_energy[iH] / ele_sc_raw_energy[iE];
          }
          else if (std::abs(hcalhit_seed_diphi[iH]) > 13 and std::abs(hcalhit_seed_diphi[iH]) < 23)
            v_nxx[ hcalhit_depth[iH] - 1 ] += hcalhit_energy[iH];
        }
      } // eol hcalhit

      // fill plots
      for (int iV = 0; iV < 7; ++iV) {
        for (auto &iP : v_iPlot) {
          ref_to(v_h3, "hse1_" + toStr(iV + 1) + toStr(iV + 1) + fEta)[iP].plot->Fill( accumulate(v_hse, iV, iV + 1), ele_pt[iE], rho, weight );
          ref_to(v_h3, "hoe1_" + toStr(iV + 1) + toStr(iV + 1) + fEta)[iP].plot->Fill( accumulate(v_hoe, iV, iV + 1), ele_pt[iE], rho, weight );
          ref_to(v_h3, "nxx1_" + toStr(iV + 1) + toStr(iV + 1) + fEta)[iP].plot->Fill( accumulate(v_nxx, iV, iV + 1), ele_pt[iE], rho, weight );
        }

        if (iV == 6) continue;

        for (auto &iP : v_iPlot) {
          ref_to(v_h3, "hse2_" + toStr(iV + 1) + toStr(iV + 2) + fEta)[iP].plot->Fill( accumulate(v_hse, iV, iV + 2), ele_pt[iE], rho, weight );
          ref_to(v_h3, "hoe2_" + toStr(iV + 1) + toStr(iV + 2) + fEta)[iP].plot->Fill( accumulate(v_hoe, iV, iV + 2), ele_pt[iE], rho, weight );
          ref_to(v_h3, "nxx2_" + toStr(iV + 1) + toStr(iV + 2) + fEta)[iP].plot->Fill( accumulate(v_nxx, iV, iV + 2), ele_pt[iE], rho, weight );
        }

        if (iV == 5) continue;

        for (auto &iP : v_iPlot) {
          ref_to(v_h3, "hse3_" + toStr(iV + 1) + toStr(iV + 3) + fEta)[iP].plot->Fill( accumulate(v_hse, iV, iV + 3), ele_pt[iE], rho, weight );
          ref_to(v_h3, "hoe3_" + toStr(iV + 1) + toStr(iV + 3) + fEta)[iP].plot->Fill( accumulate(v_hoe, iV, iV + 3), ele_pt[iE], rho, weight );
          ref_to(v_h3, "nxx3_" + toStr(iV + 1) + toStr(iV + 3) + fEta)[iP].plot->Fill( accumulate(v_nxx, iV, iV + 3), ele_pt[iE], rho, weight );
        }

        if (iV == 4) continue;

        for (auto &iP : v_iPlot) {
          ref_to(v_h3, "hse4_" + toStr(iV + 1) + toStr(iV + 4) + fEta)[iP].plot->Fill( accumulate(v_hse, iV, iV + 4), ele_pt[iE], rho, weight );
          ref_to(v_h3, "hoe4_" + toStr(iV + 1) + toStr(iV + 4) + fEta)[iP].plot->Fill( accumulate(v_hoe, iV, iV + 4), ele_pt[iE], rho, weight );
          ref_to(v_h3, "nxx4_" + toStr(iV + 1) + toStr(iV + 4) + fEta)[iP].plot->Fill( accumulate(v_nxx, iV, iV + 4), ele_pt[iE], rho, weight );
        }

        if (iV == 3) continue;

        for (auto &iP : v_iPlot) {
          ref_to(v_h3, "hse5_" + toStr(iV + 1) + toStr(iV + 5) + fEta)[iP].plot->Fill( accumulate(v_hse, iV, iV + 5), ele_pt[iE], rho, weight );
          ref_to(v_h3, "hoe5_" + toStr(iV + 1) + toStr(iV + 5) + fEta)[iP].plot->Fill( accumulate(v_hoe, iV, iV + 5), ele_pt[iE], rho, weight );
          ref_to(v_h3, "nxx5_" + toStr(iV + 1) + toStr(iV + 5) + fEta)[iP].plot->Fill( accumulate(v_nxx, iV, iV + 5), ele_pt[iE], rho, weight );
        }

        if (iV == 2) continue;

        for (auto &iP : v_iPlot) {
          ref_to(v_h3, "hse6_" + toStr(iV + 1) + toStr(iV + 6) + fEta)[iP].plot->Fill( accumulate(v_hse, iV, iV + 6), ele_pt[iE], rho, weight );
          ref_to(v_h3, "hoe6_" + toStr(iV + 1) + toStr(iV + 6) + fEta)[iP].plot->Fill( accumulate(v_hoe, iV, iV + 6), ele_pt[iE], rho, weight );
          ref_to(v_h3, "nxx6_" + toStr(iV + 1) + toStr(iV + 6) + fEta)[iP].plot->Fill( accumulate(v_nxx, iV, iV + 6), ele_pt[iE], rho, weight );
        }

        if (iV == 1) continue;

        for (auto &iP : v_iPlot) {
          ref_to(v_h3, "hse7_" + toStr(iV + 1) + toStr(iV + 7) + fEta)[iP].plot->Fill( accumulate(v_hse, iV, iV + 7), ele_pt[iE], rho, weight );
          ref_to(v_h3, "hoe7_" + toStr(iV + 1) + toStr(iV + 7) + fEta)[iP].plot->Fill( accumulate(v_hoe, iV, iV + 7), ele_pt[iE], rho, weight );
          ref_to(v_h3, "nxx7_" + toStr(iV + 1) + toStr(iV + 7) + fEta)[iP].plot->Fill( accumulate(v_nxx, iV, iV + 7), ele_pt[iE], rho, weight );

          ref_to(v_h3, "cmssw_hoe" + fEta)[iP].plot->Fill( ele_cmssw_hoe[iE], ele_pt[iE], rho, weight );
          ref_to(v_h3, "cmssw_hoe_tow" + fEta)[iP].plot->Fill( ele_cmssw_hoe_tower[iE], ele_pt[iE], rho, weight );
          ref_to(v_h3, "cmssw_hoe_5x5" + fEta)[iP].plot->Fill( ele_cmssw_hoe_5x5[iE], ele_pt[iE], rho, weight );
        }
      } // eol fill
    } // eol ele
  } // eol evt

  // ok we got the plots, draw em
  // args as in the wrapper above
  const std::string odir = "./plot/";
  auto file = std::make_unique<TFile>((odir + "ele_dR_0p15_raw_hcal_pt_rho.root").c_str(), "recreate");
  file->cd();
  for (auto &vh : v_h3) {
    for (auto &h : vh.second)
      h.plot->Write();
  }

  gROOT->ProcessLine(".q");
}
