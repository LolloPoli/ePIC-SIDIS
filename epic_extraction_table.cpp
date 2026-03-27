//-------------- ePIC extraction table for uTMDs
//--- Authors: Lorenzo Polizzi (lorenzo.polizzi@unife.it), Sara Pucillo (sara.pucillo@cern.ch)

#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include <iostream>
#include <vector>
#include <cmath>
#include "TText.h"
#include <TLatex.h>
#include <TPaveStats.h>
#include <TPaveStatsEditor.h>
#include "TPaletteAxis.h"
#include "TPolyLine.h"
#include "TStyle.h"
#include "TColor.h"
#include "ePIC_style.C"
#include <cstring>  // std::strstr

using namespace std;
namespace fs = std::filesystem;
gROOT->SetBatch(kTRUE);

// rsync -avz -e "ssh -J lpolizzi@login.jlab.org" lpolizzi@ifarm:/lustre24/expphy/volatile/clas12/lpolizzi/sidis/eic/NC_26.03/10x275/minQ2=1/ /Users/lorenzopolizzi/Desktop/PhD/epic/26.03_10x275

// for campaign 26.03 - 10x275 corrupted file: Div1: 0385; Div2: 0275, 1788; Div3: 0274, 0991, 1658; Div4: 0383; Div5: 0216, 1151, 1533

//------------- Generation of logarithmically spaced bin edges between xmin and xmax
vector<double> CreateLogBinning(int nbins, double xmin, double xmax) {
    vector<double> bin_edges(nbins + 1);
    double logxmin = log10(xmin);
    double logxmax = log10(xmax);
    double bin_width = (logxmax - logxmin) / nbins;
    for (int i = 0; i <= nbins; ++i) {
        bin_edges[i] = pow(10, logxmin + i * bin_width);
    }
    return bin_edges;
}

//------------- Adjusts the position of the statistics box on a 2D histogram (top-right)
void SetStatsBox(TH2* hist) {
    TPaveStats* stats0 = (TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
    if (stats0) {
        stats0->SetX1NDC(0.75);  //left
        stats0->SetX2NDC(0.89);  //right
        stats0->SetY1NDC(0.68);  //bottom
        stats0->SetY2NDC(0.88);  //top
    }
}
void SetStatsBox2(TH2* hist) {
    TPaveStats* stats0 = (TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
    if (stats0) {
        stats0->SetX1NDC(0.12);  //left
        stats0->SetX2NDC(0.26);  //right
        stats0->SetY1NDC(0.68);  //bottom
        stats0->SetY2NDC(0.88);  //top
    }
}

//------------- Helper: binning functions in xB–Q2 and z–P_hT
int getBinIndex_xQ2(double xB, double Q2){
    // bin_xB
    //double bin_xB[][2] = {{1e-4, 1e-3}, {1e-3, 3e-3}, {3e-3, 1e-2}, {1e-2, 4e-2}, {4e-2, 1}}; // if 10x100
    double bin_xB[][2] = {{5e-5, 3e-4}, {3e-4, 1e-3}, {1e-3, 3e-3}, {3e-3, 1e-2}, {1e-2, 1}}; // if 18x275
    // bin Q2
    vector<vector<array<double,2>>> binning_Q2_for_xB = {
        {{1,2}, {2,100}},
        {{1,2}, {2,5}, {5,100}},
        {{1,2}, {2,5}, {5,20}, {20,100}},
        {{1,2}, {2,5}, {5,20}, {20,1000}},
        {{1,5}, {5,20}, {20,1000}}
    };
    int binIndex = 1;
    for (int ix = 0; ix < 5; ix++){
        if (xB >= bin_xB[ix][0] && xB <= bin_xB[ix][1]){
            for (size_t iq = 0; iq < binning_Q2_for_xB[ix].size(); iq++){
                double Q2min = binning_Q2_for_xB[ix][iq][0];
                double Q2max = binning_Q2_for_xB[ix][iq][1];
                if (Q2 >= Q2min && Q2 < Q2max){
                    return binIndex;
                }
                binIndex++;
            }
            return -1;
        } else {
            binIndex += binning_Q2_for_xB[ix].size();
        }
    }

    return -1;
}

int getBinIndex_zPt(double z, double Pt){
    // bin_z
    double bin_z[][2] = {{0, 0.1}, {0.1, 0.2}, {0.2, 0.3}, {0.3, 0.5}, {0.5, 0.7}, {0.7, 1}};
    // bin_Pt
    vector<vector<array<double,2>>> binning_Pt_for_z = {
        {{0, 0.25}, {0.25, 0.5}, {0.5, 0.75}, {0.75, 1.25}, {1.25, 2.5}},
        {{0, 0.25}, {0.25, 0.5}, {0.5, 0.75}, {0.75, 1.25}, {1.25, 2.5}},
        {{0, 0.25}, {0.25, 0.5}, {0.5, 0.75}, {0.75, 1.25}, {1.25, 2.5}},
        {{0, 0.25}, {0.25, 0.5}, {0.5, 0.75}, {0.75, 1.25}, {1.25, 2.5}},
        {{0, 0.25}, {0.25, 0.5}, {0.5, 0.75}, {0.75, 1.25}, {1.25, 2.5}},
        {{0, 0.25}, {0.25, 0.5}, {0.5, 0.75}, {0.75, 1.25}, {1.25, 2.5}}
    };
    int binIndex = 1;
    for (int iz = 0; iz < 6; iz++){
        if (z >= bin_z[iz][0] && z <= bin_z[iz][1]){
            for (size_t ip = 0; ip < binning_Pt_for_z[iz].size(); ip++){
                double Ptmin = binning_Pt_for_z[iz][ip][0];
                double Ptmax = binning_Pt_for_z[iz][ip][1];
                if (Pt >= Ptmin && Pt < Ptmax){
                    return binIndex;
                }
                binIndex++;
            }

            return -1;
        } else {
            binIndex += binning_Pt_for_z[iz].size();
        }
    }

    return -1;
}

int getBinIndex_z(double z){
    // bin_z
    double bin_z[][2] = {{0, 0.1}, {0.1, 0.2}, {0.2, 0.3}, {0.3, 0.5}, {0.5, 0.7}, {0.7, 1}};
    int binIndex = 1;
    for (int iz = 0; iz < 6; iz++){
        if (z >= bin_z[iz][0] && z <= bin_z[iz][1]){
            return iz+1;
        }
    }
    return -1;
}
int getBinIndex_Pt(double Pt){
    // bin_Pt
    double bin_Pt[][2] = {{0, 0.25}, {0.25, 0.5}, {0.5, 0.75}, {0.75, 1.25}, {1.25, 2.5}};
    int binIndex = 1;
    for (int ip = 0; ip < 5; ip++){
        if (Pt >= bin_Pt[ip][0] && Pt <= bin_Pt[ip][1]){
            return ip+1;
        }
    }
    return -1;
}



// Gregory binning from Marco

int getBin_xB(double xB){
    vector<double> xB_edges = {0.000630957, 0.001, 0.00158489, 0.00251189, 0.00398107, 0.00630957, 0.01, 0.0158489, 0.0251189, 0.0398107, 0.0630957, 0.1, 0.158489, 0.251189, 0.398107, 0.630957, 1.0};
    for(size_t i = 0; i < xB_edges.size() - 1; i++){
        if(xB >= xB_edges[i] && xB < xB_edges[i+1]){
            return i + 1;  // bin from 1 to ...
        }
    }
    return -1;  // fuori range
}
int getBin_Q2(double Q2){
    vector<double> Q2_edges = {1.0, 1.7782800679308082, 2.3713730200033902, 3.1622776601683795, 4.216965733794858, 5.623415332340303, 7.498939925082745, 10.0, 13.335216533675034, 17.78280067930808, 23.713730200033904}; // chiedere conferma
    for(size_t i = 0; i < Q2_edges.size() - 1; i++){
        if(Q2 >= Q2_edges[i] && Q2 < Q2_edges[i+1]){
            return i + 1;  
        }
    }
    return -1;  // fuori range
}
int getBin_z(double z){
    vector<double> z_edges = {0.05, 0.1, 0.15, 0.2, 0.25, 0.30, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    for(size_t i = 0; i < z_edges.size() - 1; i++){
        if(z >= z_edges[i] && z < z_edges[i+1]){
            return i + 1;  // bin numerati da 1 a 9
        }
    }
    return -1;  // fuori range
}
int getBin_Pt(double Pt){
    vector<double> Pt_edges = {0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8}; // chiedere conferma
    for(size_t i = 0; i < Pt_edges.size() - 1; i++){
        if(Pt >= Pt_edges[i] && Pt < Pt_edges[i+1]){
            return i + 1;  // bin numerati da 1 a 9
        }
    }
    return -1;  // fuori range
}

// for the bin ranges
std::pair<double,double> getBinRange_xB(int bin){
    std::vector<double> xB_edges = {
        0.000630957, 0.001, 0.00158489, 0.00251189, 0.00398107,
        0.00630957, 0.01, 0.0158489, 0.0251189, 0.0398107,
        0.0630957, 0.1, 0.158489, 0.251189,
        0.398107, 0.630957, 1.0
    };
    if (bin < 1 || bin >= xB_edges.size()) {
        return {-1, -1}; // errore
    }
    return {xB_edges[bin-1], xB_edges[bin]};
}
std::pair<double,double> getBinRange_Q2(int bin){
    vector<double> Q2_edges = {1.0, 1.7782800679308082, 2.3713730200033902, 3.1622776601683795, 4.216965733794858, 5.623415332340303, 7.498939925082745, 10.0, 13.335216533675034, 17.78280067930808, 23.713730200033904}; // chiedere conferma
    if (bin < 1 || bin >= Q2_edges.size()) {
        return {-1, -1}; // errore
    }
    return {Q2_edges[bin-1], Q2_edges[bin]};
}
std::pair<double,double> getBinRange_z(int bin){
    vector<double> z_edges = {0.05, 0.1, 0.15, 0.2, 0.25, 0.30, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    if (bin < 1 || bin >= z_edges.size()) {
        return {-1, -1}; // errore
    }
    return {z_edges[bin-1], z_edges[bin]};
}
std::pair<double,double> getBinRange_Pt(int bin){
    vector<double> Pt_edges = {0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0}; // chiedere conferma
    if (bin < 1 || bin >= Pt_edges.size()) {
        return {-1, -1}; // errore
    }
    return {Pt_edges[bin-1], Pt_edges[bin]};
}

//------------- Main function
void epic_extraction_table(int target_pdg = 211, const char* inputDir = "26.03_10x275", int beam_e = 10, int beam_p = 275) {

    //---set_ePIC_style();
    //gROOT->ProcessLine("set_ePIC_style()");
    //gROOT->ForceStyle();

    //gROOT->ForceStyle();

    // --- Hadron selection (PDG code)
    TString tag, label;
    switch (target_pdg) {
        case  211: tag = "pipos"; label = "#pi^{+}"; break;
        case -211: tag = "pineg"; label = "#pi^{-}"; break;
        case  321: tag = "kaonpos";  label = "K^{+}";    break;
        case -321: tag = "kaonneg";  label = "K^{-}";    break;
        default:
            tag   = Form("pdg%d", target_pdg);
            label = Form("PDG %d", target_pdg);
            break;
    }

    // Variables
    // --- electron
    double electron_px, electron_py, electron_pz, electron_mom, electron_Theta, electron_Phi, electron_ThetaDeg, electron_E, electron_W, electron_Q2, electron_ass;
    double electron_eta, electron_y;
    // -- hadron
    double hadron_mom, hadron_Q2, hadron_xB, hadron_xF, hadron_z, hadron_PhT, hadron_Phi_h, hadron_Phi_s, hadron_Phi_lab, hadron_Theta, hadron_eta, hadron_y, hadron_W2, hadron_Mx;
    double helicity, eps, hadron_px, hadron_py, hadron_pz, el_px, el_py, el_pz, el_theta, el_phi, el_eta, el_mom, pr_mom, pr_px, pr_py, pr_pz, pr_phi, pr_theta, pr_eta;
    double rec_pdg, good_PID, el_rec_pdg, rec_pdg_mc;
    double el_ass_rec_pdg, el_ass_px, el_ass_py, el_ass_pz, el_ass_theta, el_ass_phi, el_ass_eta, el_ass_mom;
    double hadron_E;
    //
    double hadron_mom_mc, hadron_Q2_mc, hadron_xB_mc, hadron_xF_mc, hadron_z_mc, hadron_PhT_mc;
    double hadron_Phi_h_mc, hadron_Phi_s_mc, hadron_Phi_lab_mc, hadron_Theta_mc, hadron_eta_mc, hadron_y_mc, hadron_W2_mc, hadron_Mx_mc;
    double hel_mc, eps_mc, hadron_px_mc, hadron_py_mc, hadron_pz_mc;
    double hadron_mc_index, index_mc;
    //
    double hadron_mom_all, hadron_Q2_all, hadron_xB_all, hadron_xF_all, hadron_z_all, hadron_PhT_all;
    double hadron_Phi_h_all, hadron_Phi_s_all, hadron_Phi_lab_all, hadron_Theta_all, hadron_eta_all, hadron_y_all, hadron_W2_all, hadron_Mx_all;
    double hel_all, eps_all, hadron_px_all, hadron_py_all, hadron_pz_all;
    double good_PID_all, pdg_all;
    double hadron_index, index_all;
    // cross section
    double cross_section, cross_section_err_uncorr, cross_section_err_corr;

    //--- input file
    string inputDirStr = inputDir;
    TTree treeHadron(Form("%s_RECO", tag.Data()), Form("RECO %s", label.Data()));
    TTree treeHadron_MC(Form("%s_MC", tag.Data()), Form("MC %s", label.Data()));
    TChain chainElectron("Electron");
    TChain chainHadron_Reco("Hadron Reco");
    TChain chainHadron_MC("Hadron MC");

    int fileCount = 0;

    fs::path inPath(inputDirStr);
    if (fs::exists(inPath) && fs::is_regular_file(inPath) && inPath.extension() == ".root") {
        string filePath = inPath.string();
        chainElectron.Add(Form("%s/ElectronTree_MC", filePath.c_str()));
        chainHadron_MC.Add(Form("%s/HadronTree_MC", filePath.c_str()));
        chainHadron_Reco.Add(Form("%s/HadronTree_RECO", filePath.c_str()));
        fileCount = 1;
    } else {
        for (const auto &entry : fs::directory_iterator(inputDirStr)) {
            if (entry.path().extension() == ".root") {
                string filePath = entry.path().string();
                chainElectron.Add(Form("%s/ElectronTree_MC", filePath.c_str()));
                chainHadron_MC.Add(Form("%s/HadronTree_MC", filePath.c_str()));
                chainHadron_Reco.Add(Form("%s/HadronTree_RECO", filePath.c_str()));
                fileCount++;
            }
        }
    }

    if (fileCount == 0) {
        cerr << "No .root file found in " << inputDirStr << endl;
        return;
    }

    //--- output file
    TString outputFile = Form("epic_2603_%s_%dx%d.root", tag.Data(), beam_e, beam_p);
    TFile outFile(outputFile, "RECREATE");

    // yaml file

    TString yaml_name = Form("%s_%dx%d_26_03.yaml", tag.Data(), beam_e, beam_p);
    std::ofstream yamlFile(yaml_name);

    // csv file

    TString csv_filename = Form("%s_%dx%d_26_03.csv", tag.Data(), beam_e, beam_p);
    std::ofstream csvFile(csv_filename);
    csvFile << "id, s[GeV^{2}], <Q>[GeV], Q_min[GeV], Q_max[GeV], <xB>, xB_min, xB_max, <z>, z_min, z_max, <pT>, pT_min, pT_max, Delta_bin[GeV^3], xSec_diff[fb/GeV^3], xSec_uncorr_error, xSec_corr_error, y_min, y_max, W2_min[GeV^2], W2_max[GeV^2],target_mass[GeV], product_mass[GeV], N_gen, N_reco, abs_stat_error, abs_sys_error \n";

    //--- Here we collect all the variables from the ttree
    if (chainElectron.GetNtrees() > 0) {
        chainElectron.SetBranchAddress("el_px_mc", &electron_px);
        chainElectron.SetBranchAddress("el_py_mc", &electron_py);
        chainElectron.SetBranchAddress("el_pz_mc", &electron_pz);
        chainElectron.SetBranchAddress("el_mom_mc", &electron_mom);
        chainElectron.SetBranchAddress("el_theta_mc", &electron_Theta);
        chainElectron.SetBranchAddress("el_phi_mc", &electron_Phi);
        chainElectron.SetBranchAddress("el_eta_mc", &electron_eta);
        chainElectron.SetBranchAddress("el_y_mc", &electron_y);
    }

    //Hadron reco
    chainHadron_Reco.SetBranchAddress("hadron_index", &hadron_index);
    chainHadron_Reco.SetBranchAddress("hadron_pdg", &rec_pdg);
    chainHadron_Reco.SetBranchAddress("hadron_pdg_mc", &rec_pdg_mc);
    chainHadron_Reco.SetBranchAddress("hadron_good_PID", &good_PID);
    chainHadron_Reco.SetBranchAddress("hadron_px", &hadron_px);
    chainHadron_Reco.SetBranchAddress("hadron_py", &hadron_py);
    chainHadron_Reco.SetBranchAddress("hadron_pz", &hadron_pz);
    chainHadron_Reco.SetBranchAddress("hadron_mom", &hadron_mom);
    chainHadron_Reco.SetBranchAddress("hadron_W2", &hadron_W2);
    chainHadron_Reco.SetBranchAddress("hadron_Q2", &hadron_Q2);
    //chainHadron_Reco.SetBranchAddress("hadron_xF", &hadron_xF);
    chainHadron_Reco.SetBranchAddress("hadron_xB", &hadron_xB);
    chainHadron_Reco.SetBranchAddress("hadron_y", &hadron_y);
    chainHadron_Reco.SetBranchAddress("hadron_z", &hadron_z);
    chainHadron_Reco.SetBranchAddress("hadron_PhT", &hadron_PhT);
    chainHadron_Reco.SetBranchAddress("hadron_Phi_lab", &hadron_Phi_lab);
    chainHadron_Reco.SetBranchAddress("hadron_Theta", &hadron_Theta);
    chainHadron_Reco.SetBranchAddress("hadron_eta", &hadron_eta);
    chainHadron_Reco.SetBranchAddress("hadron_Phi_h", &hadron_Phi_h);
    //chainHadron_Reco.SetBranchAddress("Phi_s", &hadron_Phi_s);
    //chainHadron_Reco.SetBranchAddress("helicity", &helicity);
    //chainHadron_Reco.SetBranchAddress("hadron_Mx", &hadron_Mx);

    //Hadron MC
    chainHadron_MC.SetBranchAddress("hadron_index_mc", &index_mc);
    // Optional: if present, use MC PDG to select the hadron species
    double mc_pdg = 0;
    //bool has_mc_pdg = (chainHadron_MC.SetBranchAddress("hadron_pdg_mc", &mc_pdg) == 0);
    //if (!has_mc_pdg) has_mc_pdg = (chainHadron_MC.SetBranchAddress("hadron_pdg", &mc_pdg) == 0);
    chainHadron_MC.SetBranchAddress("hadron_pdg_mc", &mc_pdg);
    chainHadron_MC.SetBranchAddress("hadron_mom_mc", &hadron_mom_mc);
    chainHadron_MC.SetBranchAddress("hadron_Q2_mc", &hadron_Q2_mc);
    chainHadron_MC.SetBranchAddress("hadron_xB_mc", &hadron_xB_mc);
    //chainHadron_MC.SetBranchAddress("hadron_xF_mc", &hadron_xF_mc);
    chainHadron_MC.SetBranchAddress("hadron_z_mc", &hadron_z_mc);
    chainHadron_MC.SetBranchAddress("hadron_PhT_mc", &hadron_PhT_mc);
    chainHadron_MC.SetBranchAddress("hadron_Phi_lab_mc", &hadron_Phi_lab_mc);
    chainHadron_MC.SetBranchAddress("hadron_Phi_h_mc", &hadron_Phi_h_mc);
    //chainHadron_MC.SetBranchAddress("Phi_s_mc", &hadron_Phi_s_mc);
    chainHadron_MC.SetBranchAddress("hadron_Theta_mc", &hadron_Theta_mc);
    chainHadron_MC.SetBranchAddress("hadron_eta_mc", &hadron_eta_mc);
    chainHadron_MC.SetBranchAddress("hadron_y_mc", &hadron_y_mc);
    //chainHadron_MC.SetBranchAddress("hadron_W2_mc", &hadron_W2_mc);
    //chainHadron_MC.SetBranchAddress("hadron_Mx_mc", &hadron_Mx_mc);
    //chainHadron_MC.SetBranchAddress("helicity_mc", &hel_mc);
    //chainHadron_MC.SetBranchAddress("hadron_epsilon_mc", &eps_mc);
    chainHadron_MC.SetBranchAddress("hadron_px_mc", &hadron_px_mc);
    chainHadron_MC.SetBranchAddress("hadron_py_mc", &hadron_py_mc);
    chainHadron_MC.SetBranchAddress("hadron_pz_mc", &hadron_pz_mc);

    // Interesting hadron --- reco info
    treeHadron.Branch("mc_index", &hadron_index, "mc_index/I");
    treeHadron.Branch("rec_pdg", &rec_pdg, "rec_pdg/D");
    treeHadron.Branch("good_PID", &good_PID, "good_PID/D");
    treeHadron.Branch("px", &hadron_px, "hadron_px/D");
    treeHadron.Branch("py", &hadron_py, "hadron_py/D");
    treeHadron.Branch("pz", &hadron_pz, "hadron_pz/D");
    treeHadron.Branch("E", &hadron_E, "E/D");
    treeHadron.Branch("Mom", &hadron_mom, "Mom/D");
    treeHadron.Branch("Q2", &hadron_Q2, "Q2/D");
    treeHadron.Branch("xB", &hadron_xB, "xB/D");
    treeHadron.Branch("xF", &hadron_xF, "xF/D");
    treeHadron.Branch("z", &hadron_z, "z/D");
    treeHadron.Branch("PhT", &hadron_PhT, "PhT/D");
    treeHadron.Branch("Phi_lab", &hadron_Phi_lab, "Phi_h/D");
    treeHadron.Branch("Phi_h", &hadron_Phi_h, "Phi_h/D");
    treeHadron.Branch("Phi_s", &hadron_Phi_s, "Phi_s/D");
    treeHadron.Branch("theta", &hadron_Theta, "theta/D");
    treeHadron.Branch("eta", &hadron_eta, "eta/D");
    treeHadron.Branch("y", &hadron_y, "y/D");
    treeHadron.Branch("W2", &hadron_W2, "W/D");
    treeHadron.Branch("Mx", &hadron_Mx, "Mx/D");
    treeHadron.Branch("helicity", &helicity, "hel/D");
    treeHadron.Branch("epsilon", &eps, "eps/D");

    // Interesting hadron --- MC info
    treeHadron_MC.Branch("index", &index_mc, "index/I");
    treeHadron_MC.Branch("Mom_mc", &hadron_mom_mc, "Mom_mc/D");
    treeHadron_MC.Branch("Q2_mc", &hadron_Q2_mc, "Q2_mc/D");
    treeHadron_MC.Branch("xB_mc", &hadron_xB_mc, "xB_mc/D");
    treeHadron_MC.Branch("xF_mc", &hadron_xF_mc, "xF_mc/D");
    treeHadron_MC.Branch("z_mc", &hadron_z_mc, "z_mc/D");
    treeHadron_MC.Branch("PhT_mc", &hadron_PhT_mc, "PhT_mc/D");
    treeHadron_MC.Branch("Phi_lab_mc", &hadron_Phi_lab_mc, "Phi_lab_mc/D");
    treeHadron_MC.Branch("Phi_h_mc", &hadron_Phi_h_mc, "Phi_h_mc/D");
    treeHadron_MC.Branch("Phi_s_mc", &hadron_Phi_s_mc, "Phi_s_mc/D");
    treeHadron_MC.Branch("theta_mc", &hadron_Theta_mc, "theta_mc/D");
    treeHadron_MC.Branch("eta_mc", &hadron_eta_mc, "eta_mc/D");
    treeHadron_MC.Branch("y_mc", &hadron_y_mc, "y_mc/D");
    treeHadron_MC.Branch("W_mc", &hadron_W2_mc, "W_mc/D");
    treeHadron_MC.Branch("Mx_mc", &hadron_Mx_mc, "Mx_mc/D");
    treeHadron_MC.Branch("helicity_mc", &hel_mc, "helicity_mc/D");
    treeHadron_MC.Branch("epsilon_mc", &eps_mc, "epsilon_mc/D");
    treeHadron_MC.Branch("pion_px_mc", &hadron_px_mc, "pion_px_mc/D");
    treeHadron_MC.Branch("pion_py_mc", &hadron_py_mc, "pion_py_mc/D");
    treeHadron_MC.Branch("pion_pz_mc", &hadron_pz_mc, "pion_pz_mc/D");

    /*
    treeHadron_all.Branch("index", &index_all, "index/I");
    treeHadron_all.Branch("Mom_all", &hadron_mom_all, "Mom_all/D");
    treeHadron_all.Branch("Q2_all", &hadron_Q2_all, "Q2_all/D");
    treeHadron_all.Branch("xB_all", &hadron_xB_all, "xB_all/D");
    treeHadron_all.Branch("xF_all", &hadron_xF_all, "xF_all/D");
    treeHadron_all.Branch("z_all", &hadron_z_all, "z_all/D");
    treeHadron_all.Branch("PhT_all", &hadron_PhT_all, "PhT_all/D");
    treeHadron_all.Branch("Phi_lab_all", &hadron_Phi_lab_all, "Phi_lab_all/D");
    treeHadron_all.Branch("Phi_h_all", &hadron_Phi_h_all, "Phi_h_all/D");
    treeHadron_all.Branch("Phi_s_all", &hadron_Phi_s_all, "Phi_s_all/D");
    treeHadron_all.Branch("theta_all", &hadron_Theta_all, "theta_all/D");
    treeHadron_all.Branch("eta_all", &hadron_eta_all, "eta_all/D");
    treeHadron_all.Branch("y_all", &hadron_y_all, "y_all/D");
    treeHadron_all.Branch("W_all", &hadron_W2_all, "W_all/D");
    treeHadron_all.Branch("Mx_all", &hadron_Mx_all, "Mx_all/D");
    treeHadron_all.Branch("helicity_all", &hel_all, "helicity_all/D");
    treeHadron_all.Branch("epsilon_all", &eps_all, "epsilon_all/D");
    treeHadron_all.Branch("pion_px_all", &hadron_px_all, "pion_px_all/D");
    treeHadron_all.Branch("pion_py_all", &hadron_py_all, "pion_py_all/D");
    treeHadron_all.Branch("pion_pz_all", &hadron_pz_all, "pion_pz_all/D");
    */

    //
    // --- Distributions
    double bin = 200;
    double chi_min = 0.1;
    double chi_max = 2000;
    const double xmin_xB = 5e-5, xmax_xB = 1;
    const double xmin_Q2 = 1, xmax_Q2 = 2000.;
    auto make_bins = [](int bins, double min, double max){
        return CreateLogBinning(bins, min, max);
    };
    const auto log_chi2 = make_bins(bin, chi_min, chi_max);
    const auto log_bins_Q2 = make_bins(bin, xmin_Q2, xmax_Q2);
    const auto log_bins_xB = make_bins(bin, xmin_xB, xmax_xB);

    TH2D had_Q2VsXb_MC (Form("%s_Q2VsXb_MC", tag.Data()), Form("Correlation Q^{2} vs x_{B}  |  MC %s ; x_{B}; Q^{2} [GeV^{2}]",label.Data()), bin, log_bins_xB.data(), bin, log_bins_Q2.data());
    TH2D had_PhTvsZ_MC (Form("%s_PhTvsZ_MC", tag.Data()), Form("Correlation P_{hT} vs Z  |  MC %s ; z; P_{hT} [GeV]",label.Data()), bin, 0, 1, bin, 0, 3);

    TH1D had_evnt_chi2 (Form("%s_evnt_chi2", tag.Data()), Form("#chi^{2} EventBuilder PID | %s ; only EventBuilder | 1.2 < Mom < 8 GeV ; #chi^{2}; count",label.Data()), bin, -8, 8);
    TH1D had_m (Form("%s_best_mass", tag.Data()), Form("m extracted from #beta | %s ; 1.2 < Mom < 8 GeV ; m [GeV]; count",label.Data()), bin, 0, 1);
    TH1D had_deltaB (Form("%s_delta_beta", tag.Data()), Form("#beta_{meas} - #beta_{th} | %s ; 1.2 < Mom < 8 GeV ; #Delta_{#beta}; count",label.Data()), bin, -0.05, 0.05);
    TH1D had_Mx (Form("%s_missing_mass", tag.Data()), Form("missing mass | %s ; 1.2 < Mom < 8 GeV ; M_{x} [GeV]; count",label.Data()), bin, 0, 60);

    //--- Mom
    TH2D had_MomVsPhT (Form("%s_MomVsPhT", tag.Data()), Form("Correlation Mom vs P_{hT}  |  %s ; P_{hT} [GeV]; Mom [GeV]",label.Data()), bin, 0, 3, bin, 0, 20);
    TH2D had_MomVsXb (Form("%s_MomVsXb", tag.Data()), Form("Correlation Mom vs x_{B}  |  %s ; x_{B}; Mom [GeV]",label.Data()), bin, log_bins_xB.data(), bin, 0, 20);
    TH2D had_MomVsXf (Form("%s_MomVsXf", tag.Data()), Form("Correlation Mom vs x_{F}  |  %s ; x_{F}; Mom [GeV]",label.Data()), bin, -0.5, 0.5, bin, 0, 20);
    TH2D had_MomVsZ (Form("%s_MomVsZ", tag.Data()), Form("Correlation Mom vs Z  |  %s ; z; Mom [GeV]",label.Data()), bin, 0, 1, bin, 0, 20);
    TH2D had_MomVsY (Form("%s_MomVsY", tag.Data()), Form("Correlation Mom vs Y  |  %s ; y; Mom [GeV]",label.Data()), bin, 0.0, 1.0, bin, 0, 20);
    TH2D had_MomVsEta (Form("%s_MomVsEta", tag.Data()), Form("Correlation Mom vs Eta  |  %s ; Eta; Mom [GeV]",label.Data()), bin, -3.5, 3.5, bin, 0, 20);
    TH2D had_MomVsEta_MC (Form("%s_MomVsEta_MC", tag.Data()), Form("Correlation Mom vs Eta MC  |  %s ; Eta; Mom [GeV]",label.Data()), bin, -3.5, 3.5, bin, 0, 20);
    TH2D had_MomVsTheta (Form("%s_MomVsTheta", tag.Data()), Form("Correlation Mom vs Theta  |  %s ; Mom [GeV]; #theta [Rad]",label.Data()), bin, 0, 20, bin, 0, TMath::Pi());
    TH2D had_MomVsPhi_h (Form("%s_MomVsPhi_h", tag.Data()), Form("Correlation Mom vs #Phi_{h}  |  %s ; #Phi_{h} [Rad]; Mom [GeV]",label.Data()), bin, -TMath::Pi(), TMath::Pi(), bin, 0, 20);
    TH2D had_MomVsW2 (Form("%s_MomVsW2", tag.Data()), Form("Correlation Mom vs W^{2} | %s |; Mom [GeV]; W^{2} [GeV^{2}]",label.Data()), bin, 0, 20, bin, 0, 10000);
    //--- Q2
    TH2D had_Q2VsXb (Form("%s_Q2VsXb", tag.Data()), Form("Correlation Q^{2} vs x_{B}  |  %s ; x_{B}; Q^{2} [GeV^{2}]",label.Data()), bin, log_bins_xB.data(), bin, log_bins_Q2.data());
    TH2D had_Q2VsXf (Form("%s_Q2VsXf", tag.Data()), Form("Correlation Q^{2} vs x_{F}  |  %s ; x_{F}; Q^{2} [GeV^{2}]",label.Data()), bin, -0.5, 0.5, bin, log_bins_Q2.data());
    TH2D had_Q2VsMom (Form("%s_Q2VsMom", tag.Data()), Form("Correlation Q^{2} vs Mom  |  %s ; Mom [GeV]; Q^{2} [GeV^{2}]",label.Data()), bin, 0, 20, bin, log_bins_Q2.data());
    TH2D had_Q2VsPhT (Form("%s_Q2VsPhT", tag.Data()), Form("Correlation Q^{2} vs P_{hT}  |  %s ; P_{hT} [GeV]; Q^{2} [GeV^{2}]",label.Data()), bin, 0, 3, bin, log_bins_Q2.data());
    TH2D had_Q2VsZ (Form("%s_Q2VsZ", tag.Data()), Form("Correlation Q^{2} vs Z  |  %s ; z; Q^{2} [GeV^{2}]",label.Data()), bin, 0, 1, bin, log_bins_Q2.data());
    TH2D had_Q2VsY (Form("%s_Q2VsY", tag.Data()), Form("Correlation Q^{2} vs Y  |  %s ; y; Q^{2} [GeV^{2}]",label.Data()), bin, 0.0, 1.0, bin, log_bins_Q2.data());
    TH2D had_Q2VsEta (Form("%s_Q2VsEta", tag.Data()), Form("Correlation Q^{2} vs Eta  |  %s ; Eta; Q^{2} [GeV^{2}]",label.Data()), bin, -3.5, 3.5, bin, log_bins_Q2.data());
    TH2D had_Q2VsPhi_h (Form("%s_Q2VsPhi_h", tag.Data()), Form("Correlation Q^{2} vs #Phi_{h}  |  %s ; #Phi_{h} [Rad]; Q^{2} [GeV^{2}]",label.Data()), bin, -TMath::Pi(), TMath::Pi(), bin, log_bins_Q2.data());
    //--- PhT
    TH2D had_PhTvsZ (Form("%s_PhTvsZ", tag.Data()), Form("Correlation P_{hT} vs Z  |  %s ; z; P_{hT} [GeV]",label.Data()), bin, 0, 1, bin, 0, 3);
    TH2D had_PhTvsXb (Form("%s_PhTvsXb", tag.Data()), Form("Correlation P_{hT} vs x_{B}  |  %s ; x_{B}; P_{hT} [GeV]",label.Data()), bin, log_bins_xB.data(), bin, 0, 3);
    TH2D had_PhTvsEta (Form("%s_PhTvsEta", tag.Data()), Form("Correlation P_{hT} vs Eta  |  %s ; Eta; P_{hT} [GeV]",label.Data()), bin, -3.5, 3.5, bin, 0, 3);
    TH2D had_PhTvsPhi_h (Form("%s_PhTvsPhi_h", tag.Data()), Form("Correlation P_{hT} vs #Phi_{h}  |  %s ; #Phi_{h} [Rad]; P_{hT} [GeV]",label.Data()), bin, -TMath::Pi(), TMath::Pi(), bin, 0, 3);
    //--- z
    TH2D had_zVsXb (Form("%s_zVsXb", tag.Data()), Form("Correlation Z vs x_{B}  |  %s ; x_{B}; z",label.Data()), bin, log_bins_xB.data(), bin, 0, 1);
    TH2D had_zVsXf (Form("%s_zVsXf", tag.Data()), Form("Correlation Z vs x_{F}  |  %s ; x_{F}; z",label.Data()), bin, -0.5, 0.5, bin, 0, 1);
    TH2D had_zVsEta (Form("%s_zVsEta", tag.Data()), Form("Correlation Z vs Eta  |  %s ; Eta; z",label.Data()), bin, -3.5, 3.5, bin, 0, 1);
    TH2D had_zVsPhi_h (Form("%s_zVsPhi_h", tag.Data()), Form("Correlation Z vs #Phi_{h}  |  %s ; #Phi_{h} [Rad]; z",label.Data()), bin, -TMath::Pi(), TMath::Pi(), bin, 0, 1);

    TH2D had_xBvsY (Form("%s_xBvsY", tag.Data()), Form("Correlation y vs x_{B}  |  %s ; x_{B}; y",label.Data()), bin, log_bins_xB.data(), bin, 0.0, 1.0);
    //--- Angles
    TH2D had_ThetaVsPhi_h (Form("%s_ThetaVsPhi_h", tag.Data()), Form("Correlation Theta vs #Phi_{h}  |  %s ; #Phi_{h} [Rad]; #theta [Rad]",label.Data()), bin, -TMath::Pi(), TMath::Pi(), bin, 0, TMath::Pi());
    TH2D had_ThetaVsPhi_Lab (Form("%s_ThetaVsPhi_Lab", tag.Data()), Form("Correlation #theta vs #Phi_{Lab} | %s; #Phi_{Lab} [Rad]; #theta [Rad]",label.Data()), bin, -TMath::Pi(), TMath::Pi(), bin, 0, TMath::Pi());

    int nBin_xQ2 = 16;
    int nBin_zPt = 30;
    vector<vector<double>> had_efficiency_xQ2(nBin_xQ2);
    vector<vector<double>> had_efficiency_zPt(nBin_zPt);
    vector<vector<vector<double>>> had_efficiency_xQ2_zPt(nBin_xQ2, vector<vector<double>> (nBin_zPt));
    vector<vector<vector<double>>> had_efficiency_xQ2_zPt_2(nBin_xQ2, vector<vector<double>> (nBin_zPt));
    vector<vector<vector<double>>> had_efficiency_xQ2_zPt_mc(nBin_xQ2, vector<vector<double>> (nBin_zPt));
    vector<vector<vector<double>>> had_purity_xQ2_zPt_num(nBin_xQ2, vector<vector<double>> (nBin_zPt));
    vector<vector<vector<double>>> had_purity_xQ2_zPt_den(nBin_xQ2, vector<vector<double>> (nBin_zPt));


    // vector for yaml file
    vector<int> v_index;
    vector<double> v_Q2_min;
    // 4D
    vector<double> xB_edges = {0.000630957, 0.001, 0.00158489, 0.00251189, 0.00398107, 0.00630957, 0.01, 0.0158489, 0.0251189, 0.0398107, 0.0630957, 0.1, 0.158489, 0.251189, 0.398107, 0.630957, 1.0};
    vector<double> Q2_edges = {1.0, 1.7782800679308082, 2.3713730200033902, 3.1622776601683795, 4.216965733794858, 5.623415332340303, 7.498939925082745, 10.0, 13.335216533675034, 17.78280067930808, 23.713730200033904};
    vector<double> z_edges = {0.05, 0.1, 0.15, 0.2, 0.25, 0.30, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    vector<double> Pt_edges = {0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0};
    int n_xB = xB_edges.size() - 1;
    int n_Q2 = Q2_edges.size() - 1;
    int n_z = z_edges.size() - 1;
    int n_Pt = Pt_edges.size() - 1;
    //vector<vector<vector<vector<vector<double>>>>> vec_multi_Q2(n_Q2, vector<vector<vector<vector<double>>>> (n_xB, vector<vector<vector<double>>> (n_z, vector<vector<double>> (n_Pt))));
    struct EventData {
        double Q;
        double xB;
        double z;
        double Pt;
    };
    std::map<std::tuple<int,int,int,int>, std::vector<EventData>> data_multi;
    std::map<std::tuple<int,int,int,int>, std::vector<EventData>> data_multi_MC;

    vector<double> total_gen_xB;
    
    // --- Filling of the graphs
    // Hadron of interest
    Long64_t nEntries_had = chainHadron_Reco.GetEntries();
    Long64_t nEntries_hadMC = chainHadron_MC.GetEntries();

    for (Long64_t i = 0; i < nEntries_hadMC; i++) {
        chainHadron_MC.GetEntry(i);
        if (i % 100000 == 0) cout << "MC entry: " << i << "/" << nEntries_hadMC << endl;
        if(hadron_y_mc <= 0.95 && hadron_y_mc >= 0.01 && hadron_Q2_mc >= 1 && hadron_z_mc < 1){
            if (mc_pdg != target_pdg) continue;
            double bin_xQ2 = getBinIndex_xQ2(hadron_xB_mc, hadron_Q2_mc);
            double bin_zPt = getBinIndex_zPt(hadron_z_mc, hadron_PhT_mc);
            if(bin_xQ2 >= 0){
                if(bin_zPt >= 0){
                    had_efficiency_xQ2_zPt_mc[bin_xQ2-1][bin_zPt-1].push_back(hadron_y_mc); // compenso il fatto che i bin partano da 1
                }
            }
            had_Q2VsXb_MC.Fill(hadron_xB_mc, hadron_Q2_mc);
            had_PhTvsZ_MC.Fill(hadron_z_mc, hadron_PhT_mc);
            had_MomVsEta_MC.Fill(hadron_eta_mc, hadron_mom_mc);
            treeHadron_MC.Fill();
            int index_Q2_mc = getBin_Q2(sqrt(hadron_Q2_mc));
            int index_xB_mc = getBin_xB(hadron_xB_mc);
            int index_z_mc = getBin_z(hadron_z_mc);
            int index_Pt_mc = getBin_Pt(hadron_PhT_mc);
            if(index_Q2_mc > 0 && index_xB_mc > 0 && index_z_mc > 0 && index_Pt_mc > 0){
                //vec_multi_Q2[index_Q2-1][index_xB-1][index_z-1][index_Pt-1].push_back(sqrt(hadron_Q2));
                total_gen_xB.push_back(hadron_xB_mc);
                auto key_mc = std::make_tuple(index_Q2_mc, index_xB_mc, index_z_mc, index_Pt_mc);
                EventData ev_mc;
                ev_mc.Q  = sqrt(hadron_Q2_mc);
                ev_mc.xB = hadron_xB_mc;
                ev_mc.z  = hadron_z_mc;
                ev_mc.Pt = hadron_PhT_mc;

                data_multi_MC[key_mc].push_back(ev_mc);
            }
        }
    }

    /*
    for (Long64_t i = 0; i < nEntries_hadAll; i++) {
        chainHadron_all.GetEntry(i);
        if(good_PID_all != 0) continue;
        if(hadron_y_all <= 0.99 && hadron_y_all >= 0.01){ //same request as in epic_studies
            double bin_xQ2 = getBinIndex_xQ2(hadron_xB_all, hadron_Q2_all);
            double bin_zPt = getBinIndex_zPt(hadron_z_all, hadron_PhT_all);
            //
            if(bin_xQ2 >= 0){
                if(bin_zPt >= 0){
                    had_purity_xQ2_zPt_den[bin_xQ2-1][bin_zPt-1].push_back(hadron_y_all);
                    if(pdg_all == 211) had_purity_xQ2_zPt_num[bin_xQ2-1][bin_zPt-1].push_back(hadron_y_all); // compenso il fatto che i bin partano da 1
                }
            }
            //treeHadron_all.Fill();
        }
    }
        */

    for (Long64_t i = 0; i < nEntries_had; i++) {
        chainHadron_Reco.GetEntry(i);
        if (i % 100000 == 0) cout << "RECO entry: " << i << "/" << nEntries_had << endl;
            if(rec_pdg == target_pdg && hadron_z < 1 && hadron_Q2 >= 1){
                if(hadron_y <= 0.99 && hadron_y >= 0.01){ //same request as in epic_studies
                    double bin_xQ2 = getBinIndex_xQ2(hadron_xB, hadron_Q2);
                    double bin_zPt = getBinIndex_zPt(hadron_z, hadron_PhT);
                    if(bin_xQ2 >= 0){
                        if(bin_zPt >= 0){ // here we have reconstructed pion, but we do not ask if they are also MC pion -> contamination
                            had_efficiency_xQ2_zPt_2[bin_xQ2-1][bin_zPt-1].push_back(hadron_y); // compenso il fatto che i bin partano da 1
                            had_purity_xQ2_zPt_den[bin_xQ2-1][bin_zPt-1].push_back(hadron_y_all);
                        }
                    }
                }
            }
            // to understand if we have to apply good_PID == 0 or not!
            if(hadron_y <= 0.95 && hadron_y >= 0.01 && good_PID != -1 && rec_pdg == target_pdg && rec_pdg_mc == target_pdg && hadron_z < 1 && hadron_Q2 >= 1){ //same request as in epic_studies + goodPID + specific pdg
                double bin_xQ2 = getBinIndex_xQ2(hadron_xB, hadron_Q2);
                double bin_zPt = getBinIndex_zPt(hadron_z, hadron_PhT);
                double bin_z = getBinIndex_z(hadron_z);
                double bin_Pt = getBinIndex_Pt(hadron_PhT);
                if(bin_xQ2 >= 0){
                    if(bin_zPt >= 0){
                        had_efficiency_xQ2_zPt[bin_xQ2-1][bin_zPt-1].push_back(hadron_y);
                        had_purity_xQ2_zPt_num[bin_xQ2-1][bin_zPt-1].push_back(hadron_y_all);
                    }
                }
                // Mom
                had_MomVsPhT.Fill(hadron_PhT, hadron_mom);
                had_MomVsEta.Fill(hadron_eta, hadron_mom);
                had_MomVsW2.Fill(hadron_mom, hadron_W2);
                had_MomVsPhi_h.Fill(hadron_Phi_h, hadron_mom);
                had_MomVsTheta.Fill(hadron_mom, hadron_Theta);
                had_MomVsXb.Fill(hadron_xB, hadron_mom);
                had_MomVsXf.Fill(hadron_xF, hadron_mom);
                had_MomVsY.Fill(hadron_y, hadron_mom);
                had_MomVsZ.Fill(hadron_z, hadron_mom);
                // Q2
                had_Q2VsEta.Fill(hadron_eta, hadron_Q2);
                had_Q2VsMom.Fill(hadron_mom, hadron_Q2);
                had_Q2VsPhi_h.Fill(hadron_Phi_h, hadron_Q2);
                had_Q2VsPhT.Fill(hadron_PhT, hadron_Q2);
                had_Q2VsXb.Fill(hadron_xB, hadron_Q2);
                had_Q2VsXf.Fill(hadron_xF, hadron_Q2);
                had_Q2VsY.Fill(hadron_y, hadron_Q2);
                had_Q2VsZ.Fill(hadron_z, hadron_Q2);
                // PhT
                had_PhTvsEta.Fill(hadron_eta, hadron_PhT);
                had_PhTvsXb.Fill(hadron_xB, hadron_PhT);
                had_PhTvsZ.Fill(hadron_z, hadron_PhT);
                had_PhTvsPhi_h.Fill(hadron_Phi_h, hadron_PhT);
                // z
                had_zVsEta.Fill(hadron_eta, hadron_z);
                had_zVsPhi_h.Fill(hadron_Phi_h, hadron_z);
                had_zVsXb.Fill(hadron_xB, hadron_z);
                had_zVsXf.Fill(hadron_xF, hadron_z);
                //
                had_xBvsY.Fill(hadron_xB, hadron_y);
                had_ThetaVsPhi_h.Fill(hadron_Phi_h, hadron_Theta);
                had_ThetaVsPhi_Lab.Fill(hadron_Phi_lab, hadron_Theta);
                treeHadron.Fill();


                // yaml file
                v_index.push_back(i);
                int index_Q2 = getBin_Q2(sqrt(hadron_Q2));
                int index_xB = getBin_xB(hadron_xB);
                int index_z = getBin_z(hadron_z);
                int index_Pt = getBin_Pt(hadron_PhT);
                if(index_Q2 > 0 && index_xB > 0 && index_z > 0 && index_Pt > 0){
                    //vec_multi_Q2[index_Q2-1][index_xB-1][index_z-1][index_Pt-1].push_back(sqrt(hadron_Q2));
                    auto key = std::make_tuple(index_Q2, index_xB, index_z, index_Pt);
                    EventData ev;
                    ev.Q  = sqrt(hadron_Q2);
                    ev.xB = hadron_xB;
                    ev.z  = hadron_z;
                    ev.Pt = hadron_PhT;

                    data_multi[key].push_back(ev);
                }
            }
            
        }

        treeHadron.Write();
        treeHadron_MC.Write();

        vector<TH2D*> hists_had = {
        &had_MomVsPhT, &had_MomVsXb, &had_MomVsXf, &had_MomVsZ, &had_MomVsY, &had_MomVsEta, &had_MomVsEta_MC,
        &had_MomVsTheta, &had_MomVsPhi_h, &had_MomVsW2,
        &had_Q2VsXb, &had_Q2VsXb_MC, &had_Q2VsXf, &had_Q2VsMom, &had_Q2VsPhT, &had_Q2VsZ, &had_Q2VsY, &had_Q2VsEta, &had_Q2VsPhi_h,
        &had_PhTvsZ, &had_PhTvsZ_MC, &had_PhTvsXb, &had_PhTvsEta, &had_PhTvsPhi_h,
        &had_zVsXb, &had_zVsXf, &had_zVsEta, &had_zVsPhi_h,
        &had_xBvsY, &had_ThetaVsPhi_h, &had_ThetaVsPhi_Lab,
    };

    // list to set the statbox2 on the canvas (suffix-based IDs)
    set<string> id_box2 = { "_MomVsXf", "_MomVsZ", "_MomVsY", "_MomVsEta", "_Q2VsXb", "_Q2VsY", "_zVsEta", "_zVsPhi_h", "_MomVsMass_RICH" };
    set<string> id_box3 = { "_MomVsXb", "_Q2VsXb", "_Q2VsXb_MC", "_PhTvsXb", "_zVsXb", "_xBvsY" };
    set<string> id_box4 = { "_Q2VsXb", "_Q2VsXb_MC", "_Q2VsXf", "_Q2VsMom", "_Q2VsPhT", "_Q2VsZ", "_Q2VsY", "_Q2VsEta", "_Q2VsPhi_h" };

    for (size_t i = 0; i < hists_had.size(); ++i) {

        const std::string histName = hists_had[i]->GetName();

        std::string key = histName;
        const std::string prefix = std::string(tag.Data()) + "_";
        if (key.rfind(prefix, 0) == 0) {               // starts_with(prefix)
            key = key.substr(prefix.size() - 1);       // keep leading "_"
        }

        TCanvas *c = new TCanvas(histName.c_str(), histName.c_str(), 800, 600);
        c->SetLogz();

        hists_had[i]->SetStats(0);
        hists_had[i]->SetTitle("");
        hists_had[i]->Draw("colz");

        if (id_box3.count(key)) c->SetLogx();
        if (id_box4.count(key)) c->SetLogy();

        c->Write();
    }


    // cose
    // costruisci insieme di tutte le chiavi
    std::set<std::tuple<int,int,int,int>> all_keys;
    // reco
    for (const auto& [key, _] : data_multi) all_keys.insert(key);
    // mc
    for (const auto& [key, _] : data_multi_MC) all_keys.insert(key);

    yamlFile << "data:\n";
    int n_id = 1;
    double L = 5; // luminositiy, proviamo 5 fb-1
    // for the generated cross-section
    // root -l -b -q /volatile/eic/EPIC/EVGEN/DIS/NC/10x275/minQ2=1/pythia8NCDIS_10x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_2.hepmc3.tree.root -e 'hepmc3_tree->Scan("hepmc3_event.attribute_string", "hepmc3_event.attribute_id==0 && hepmc3_event.attribute_name ==\"GenCrossSection\"", "colsize=30", 1, hepmc3_tree->GetEntries()-1)'
    double total_cc = (6.69559512e+05 + 6.69025860e+05 + 6.69133908e+05 + 6.69334943e+05 + 6.69510172e+05)/5;
    //hiDiv1: 6.69559512e+05, hiDiv2: 6.69025860e+05, hiDiv3: 6.69133908e+05, hiDiv4: 6.69334943e+05, hiDiv5: 6.69510172e+05
    // error cc - errore statistico globale del MC
    //hiDiv1: 4.13015055e+02, hiDiv2: 4.12506870e+02, hiDiv3: 4.12930794e+02, hiDiv4: 4.12961825e+02, hiDiv5: 4.13089841e+02
    double total_gen_count = total_gen_xB.size();
    
    for (const auto& key : all_keys) {

        auto [q2, x, z, pt] = key;

        // -------- RECO --------
        std::vector<EventData> values;
        if (data_multi.count(key))
            values = data_multi[key];

        // -------- MC --------
        int N_gen = data_multi_MC.count(key) ? data_multi_MC[key].size() : 0;

        int N_rec = values.size();
        if(N_rec == 0) continue;

        // -------- medie (RECO) --------
        double sum_Q = 0.0, sum_xB = 0.0, sum_z = 0.0, sum_Pt = 0.0;

        for (const auto& ev : values) {
            sum_Q  += ev.Q;
            sum_xB += ev.xB;
            sum_z  += ev.z;
            sum_Pt += ev.Pt;
        }

        double mean_Q  = (N_rec > 0) ? sum_Q  / N_rec : -1;
        double mean_xB = (N_rec > 0) ? sum_xB / N_rec : -1;
        double mean_z  = (N_rec > 0) ? sum_z  / N_rec : -1;
        double mean_Pt = (N_rec > 0) ? sum_Pt / N_rec : -1;

        // -------- bin edges --------
        auto [Qmin, Qmax]   = getBinRange_Q2(q2);
        auto [xBmin, xBmax] = getBinRange_xB(x);
        auto [zmin, zmax]   = getBinRange_z(z);
        auto [ptmin, ptmax] = getBinRange_Pt(pt);

        double delta = (Qmax-Qmin)*(xBmax-xBmin)*(zmax-zmin)*(ptmax-ptmin);

        // -------- cross section --------
        double sigma = 0.0;
        double sigma_stat = 0.0;
        double sigma_sys = 0.0;

        if (delta > 0 && N_rec > 0 && N_gen > 0) {
            //sigma = N_gen / (L * delta);  // N_rec / (L*eps*delta) with eps = efficiency = N_rec/N_gen
            // new sigma from total cross-section
            sigma = total_cc * (N_gen/total_gen_count) * (1/delta);
            double sigma2 = sigma*sigma;
            // errore statistico
            sigma_stat = sigma / sqrt(N_gen); // approx
            //sigma_stat = sqrt(((2*sigma2)/N_rec) - (sigma2/N_gen)); // not approx? 
        }

        // errore sistematico (capire cosa mettere)
        sigma_sys = sigma * 0.03; // 3%?

        double prod_mass = 0;
        if (std::abs(target_pdg) == 321) prod_mass = 0.497;
        else if (std::abs(target_pdg) == 211) prod_mass = 0.139;

        // YAML
        //yamlFile << "  - bin: [" << q2 << ", " << x << ", " << z << ", " << pt << "]\n";
        yamlFile << "  - bin:\n";
        yamlFile << "      s_GeV2: " << 11000 << "\n";

        yamlFile << "      Q:\n";
        yamlFile << "        Q_mean: " << mean_Q << "\n";
        yamlFile << "        Q_min: " << Qmin << "\n";
        yamlFile << "        Q_max: " << Qmax << "\n";

        yamlFile << "      xB:\n";
        yamlFile << "        xB_mean: "<< mean_xB <<"\n";
        yamlFile << "        xB_min: " << xBmin << "\n";
        yamlFile << "        xB_max: " << xBmax << "\n";

        yamlFile << "      z:\n";
        yamlFile << "        z_mean: " << mean_z << "\n";
        yamlFile << "        z_min: " << zmin << "\n";
        yamlFile << "        z_max: " << zmax << "\n";

        yamlFile << "      Pt:\n";
        yamlFile << "        Pt_mean: " << mean_Pt << "\n";
        yamlFile << "        Pt_min: " << ptmin << "\n";
        yamlFile << "        Pt_max: " << ptmax << "\n";

        yamlFile << "      Delta_bin: "<< delta <<"\n";
        yamlFile << "      xSec_diff: "<< sigma <<"\n";
        yamlFile << "      xSec_uncorr_error: " << sigma_stat << "\n";
        yamlFile << "      xSec_corr_error: " << sigma_sys << "\n";

        yamlFile << "      y_min: " << 0.01 << "\n";
        yamlFile << "      y_max: " << 0.95 << "\n";

        yamlFile << "      W2_min: " << 10.0 << "\n";
        yamlFile << "      W2_max: " << 10000.0 << "\n";

        yamlFile << "      Target_mass_GeV: " << 0.938 << "\n";
        yamlFile << "      Product_mass_GeV: " << prod_mass << "\n";

        yamlFile << "      Num_gen: " << N_gen << "\n";
        yamlFile << "      Num_reco: " << N_rec << "\n";
        yamlFile << "      abs_stat_error: " << 1/sqrt(N_rec) << "\n";
        yamlFile << "      abs_sys_error: " << 0.0 << "\n";


        // CSV
        csvFile << n_id++ << "," << 11000 << "," << mean_Q << "," << Qmin << "," << Qmax << "," << mean_xB << "," << xBmin << "," << xBmax << "," << 
        mean_z << "," << zmin << "," << zmax << "," << mean_Pt << "," << ptmin << "," << ptmax << "," << delta << "," << sigma << "," << sigma_stat << "," << sigma_sys << "," << 
        0.01 << "," << 0.95 << "," << 10.0 << "," << 10000.0 << "," << 0.938 << "," << prod_mass << "," << N_gen << "," << N_rec << "," << 1/sqrt(N_rec) << "," << 0.0 << "\n";

    }

    //outFile.Write();
    outFile.Close();
    //chain.Close();
    TString o_name = Form("%s_%dx%d_26_03", tag.Data(), beam_e, beam_p);
    cout << "-------------------------------------------" << endl;
    cout << "ROOT output file: " << outputFile << endl;
    cout << "yaml and csv file named as: " << o_name << endl;
}