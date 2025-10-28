#include <math.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>
#include <cstring>
#include "TGraph.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TApplication.h"
#include "TMultiGraph.h"
#include "TFeldmanCousins.h"
#include "TGaxis.h"
#include "TLeaf.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TStyle.h>
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TROOT.h"
#include <TText.h>
#include <TLatex.h>
#include <TRandom3.h>
#include "TRandom.h"
#include <TLegend.h>
#include <TParameter.h>
#include <TSpectrum.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <TSystem.h>
#include "TROOT.h"
using namespace std;

std::vector<TH1D*> create_data_spectrum(){
  TFile* file = TFile::Open("crop_file.root");
  TTree* tree = (TTree*)file->Get("Event");
  vector<double>* charge_gamma_finale = nullptr;
  vector<double>* dt_gamma_elec_real_E = nullptr;
  vector<int>* num_om_gamma_f = nullptr;
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("charge_gamma_finale",1);
  tree->SetBranchAddress("charge_gamma_finale", &charge_gamma_finale);
  tree->SetBranchStatus("num_om_gamma_f",1);
  tree->SetBranchAddress("num_om_gamma_f", &num_om_gamma_f);
  
  std::vector<TH1D*> histograms;
  for(int i=0; i<712; i++){   
    histograms.push_back(new TH1D(Form("om_%d",i),Form("om_%d",i),200,0,60000));
  }
  
  for(int entry=0; entry<tree->GetEntries();entry++){
    tree->GetEntry(entry);
    for(size_t k=0; k<charge_gamma_finale->size(); k++){
      histograms[num_om_gamma_f->at(k)]->Fill(charge_gamma_finale->at(k));
    }
  }
  //file->Close();
  return histograms; 
}




std::vector<std::vector<TH1D*>> create_simu_spectrum(int nb_gain_values){
  TFile *fout = new TFile("save_simu_spectra.root", "RECREATE");

  TFile* file = TFile::Open("/sps/nemo/scratch/granjon/magnetic_field/1e_simu/simu_xalbat/calib_gamma/files_root/crop_tot.root");
  TTree* tree = (TTree*)file->Get("Event");
  vector<double>* total_energy_gamma = nullptr;
  vector<double>* gamma_evis = nullptr;
  vector<int>* num_om_gamma_f = nullptr;
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("gamma_evis",1);
  tree->SetBranchAddress("gamma_evis", &gamma_evis);  
  tree->SetBranchStatus("num_om_gamma_f",1);
  tree->SetBranchAddress("num_om_gamma_f", &num_om_gamma_f);
      
  
  std::vector<std::vector<TH1D*>> histo_simu(712);
  for(int i=0; i<712; i++){
    for(int j=0; j<nb_gain_values; j++){
      histo_simu[i].push_back(new TH1D(Form("om_%d_gain_%d",i,j),Form("om_%d_gain_%d",i,j),200,0,2));
    }
  }
  
  for(int entry=0; entry<tree->GetEntries();entry++){
    tree->GetEntry(entry);
    for(size_t k=0; k<gamma_evis->size(); k++){      
      for (int i = 0; i < nb_gain_values; i++) {
	float gain = 0.5 + i * (1.0 / nb_gain_values);
	//if(gamma_evis->at(k)*gain<0.15) continue;
 	histo_simu[num_om_gamma_f->at(k)][i]->Fill(gamma_evis->at(k)*gain);
      }
    }
  }

  fout->cd();
  for(int i=0; i<712; i++){
    for(int j=0; j<nb_gain_values; j++){
      histo_simu[i][j]->Write();
    }
  }
  fout->Close();
  return histo_simu;
}



void fit_chi2_spectra(std::vector<TH1D*> data, int nb_gain_values){
 TFile *fout = new TFile("chi2_results.root", "RECREATE");
  TTree *tree = new TTree("results", "Chi2 comparison results");
  int om=0;
  double chi2=0.0, gain=0.0, hand_chi2=0.0, bin_chi2=0.0;
  TH1D* h_ref1;
  TH1D* h_simu_chi2;
  tree->Branch("om", &om);
  tree->Branch("chi2", &chi2);
  tree->Branch("hand_chi2", &hand_chi2);
  tree->Branch("bin_chi2", &bin_chi2);
  tree->Branch("gain", &gain);

  TFile *f = TFile::Open("save_simu_spectra.root");

  for (int i = 0; i < 712; i++) {// loop on OMs
    om = i;
    TH1D* h_ref = data[i];
    if (!h_ref) continue;
    if(h_ref->GetEntries() == 0) continue;
    for(int k=1; k<16; k++){//16 is because we want to supress 150 keV = 15 bins bc 200 bins = 2MeV
      h_ref->SetBinContent(k, 0);
    }
    TH1D* h_ref1 = new TH1D(Form("h_ref1_%d",i),"Data rebinned",200,0,2);
    h_ref1->SetDirectory(0);
    //project histo charge data in energy range
    for(int k=1; k<=200; k++){
      double bin_content = h_ref->GetBinContent(k);
      h_ref1->SetBinContent(k, bin_content);
    }
    for (int j = 0; j < nb_gain_values; j++) {//loop on gain values
      chi2=0.0;
      hand_chi2=0.0;
      bin_chi2=0.0;
      gain = 0.5+1.0*j/nb_gain_values;
      TH1D* h_simu = (TH1D*)f->Get(Form("om_%d_gain_%d", i, j));
      if (!h_simu) continue;
      int data_entries = h_ref->GetEntries();
      if(data_entries==0){continue;}
      for(int k=1; k<16; k++){
        h_simu->SetBinContent(k, 0);
      }
      TH1D* h_simu_chi2 = new TH1D(Form("h_simu_%d_%d",i,j),"Simu rebinned",200,0,2);
      h_simu_chi2->SetDirectory(0);
      int simu_entries = h_simu->GetEntries();
      if(simu_entries==0){continue;}
      double sum_data = 0;
      double sum_simu = 0;
      for(int k = 1; k <= 200; k++){
        sum_data += h_ref->GetBinContent(k);
        sum_simu += h_simu->GetBinContent(k);
      }
      for(int k=1; k<=200;k++){//normalised simu on data
        double bin_content_simu = h_simu->GetBinContent(k) * (1.0 * sum_data / sum_simu);
        h_simu_chi2->SetBinContent(k,bin_content_simu);
      }
      
      //chi2 by hand
      double s_data=0;
      double s_simu=0;
      for(int k=1; k<=200; k++){
        double d = h_ref1->GetBinContent(k);
        double s = h_simu_chi2->GetBinContent(k);
	s_simu+=s;
	s_data+=d;
        if(s>0 && d>0){
          hand_chi2 += (d - s) * (d - s) / s; // Pearson's chi2	 
          bin_chi2 += 2*(s-d+d*log(d/s));
	}
      }
      chi2 = h_ref1->Chi2Test(h_simu_chi2, "CHI2 UU");
      //Draw part
      // if(i<10){
      // 	TCanvas* canvas = new TCanvas();
      // 	canvas->cd();
      // 	h_ref1->SetLineColor(kRed);
      // 	h_ref1->Draw();
      // 	h_simu_chi2->SetLineColor(kBlue);
      // 	h_simu_chi2->Draw("same");
      // 	canvas->SetLogy();
      // 	canvas->SaveAs(Form("png_fit_save/om_%d_gain_%f.png",i,gain));
      // }
      tree->Fill();
      delete h_simu_chi2;
    }//end j simu gain                                                                               
    delete h_ref1;
  }//end i om
  fout->cd();
  tree->Write();
  fout->Close();
}


void extract_gain_value(){
  float fit_range = 0.1;
  TFile *fin = TFile::Open("chi2_results.root", "READ");
  TTree *tree = (TTree*)fin->Get("results");
  int om;
  double chi2, chi2_square, gain;
  tree->SetBranchAddress("om", &om);
  tree->SetBranchAddress("chi2", &chi2);
  tree->SetBranchAddress("chi2_square", &chi2_square);
  tree->SetBranchAddress("gain", &gain);
  std::map<int, std::vector<std::pair<double,double>>> points_map;
  std::map<int, double> min_chi2_map;
  std::map<int, double> min_gain_map;
  for(Long64_t i=0; i<tree->GetEntries(); ++i){
    tree->GetEntry(i);
    if(chi2==0 || gain ==1 ) continue;
    points_map[om].emplace_back(gain, chi2);
    if(min_chi2_map.find(om) == min_chi2_map.end() || chi2 < min_chi2_map[om]){
      min_chi2_map[om] = chi2;
      min_gain_map[om] = gain;
    }
  }
  TFile *fout = new TFile("fit_results.root","RECREATE");
  TTree *fitTree = new TTree("fitResults","Gain from fit");
  int om_fit;
  double gain_fit, chi2_fit, error_moins, error_plus, gain_fit_extrapolate;
  fitTree->Branch("om",&om_fit);
  fitTree->Branch("gain_fit",&gain_fit);
  fitTree->Branch("chi2_fit",&chi2_fit);
  fitTree->Branch("gain_fit_extrapolate",&gain_fit_extrapolate);
  fitTree->Branch("error_moins",&error_moins);
  fitTree->Branch("error_plus",&error_plus);

  for(int i=0; i<712; i++){
    auto &points = points_map[i];
    int n = points.size();
    TGraph *g = new TGraph(n);
    for(int j=0; j<n; j++){
      g->SetPoint(j, points[j].first, points[j].second);
    }
    double chi2_min = min_chi2_map[i];
    double gain_min = min_gain_map[i];
    double fit_low = gain_min - fit_range;
    double fit_high = gain_min + fit_range;
    TF1 *fit = new TF1(Form("fit_om_%d",i),"pol2",fit_low, fit_high);
    g->Fit(fit,"RQ");
    om_fit = i;
    gain_fit_extrapolate = -fit->GetParameter(1)/(2*fit->GetParameter(2)); // min du pol2 : -b/(2a)  
    // chi2_fit = fit->Eval(gain_fit);                                                               
    chi2_square = fit->GetChisquare();
    gain_fit = fit->GetMinimumX();
    chi2_fit = fit->GetMinimum();
    error_plus = fit->GetX(chi2_fit + 1, gain_fit, gain_fit + 0.05) - gain_fit;
    error_moins = gain_fit - fit->GetX(chi2_fit + 1, gain_fit - 0.05, gain_fit);
    fitTree->Fill();
    TCanvas *c = new TCanvas(Form("om_%d",i), Form("OM %d",i), 800,600);
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kBlack);
    g->Draw("AP");
    fit->Draw("same");
    c->SaveAs(Form("extract_values_png/om_%d.root",i));
  }

  fitTree->Write();
  fout->Close();
  fin->Close();
}





  
//int main(int argc, char const *argv[]){
int main(int argc, char** argv) {
  //TApplication app("app", &argc, argv); // initialise ROOT correctement
  gSystem->Load("libHist");
  gSystem->Load("libGraf");
  gSystem->Load("libGraf3d");
  gSystem->Load("libGpad");
  gSystem->Load("libTree");
  gSystem->Load("libRint");
  gSystem->Load("libPostscript");
  gSystem->Load("libMatrix");
  gSystem->Load("libPhysics");
  gSystem->Load("libMathCore");
  gSystem->Load("libThread");
  gSystem->Load("libCore");
  gSystem->Load("libRIO");
  //
  //first quick scan to get the minimum 
  int nb_gain_scan = 30; //number of gain you want to scan between 0.5 and 1.5
  std::vector<TH1D*> data = create_data_spectrum();
  cout<<"data spectrum created"<<endl;
  //
  //simulation spectra stay the sames, if you want to re-create it uncomment this line
  //std::vector<std::vector<TH1D*>> simu = create_simu_spectrum(nb_gain_scan);
  //cout<<"simulation spectrum created "<<endl;
  //
  fit_chi2_spectra(data, nb_gain_scan);
  cout<<"data and simu chi2 fitted "<<endl;
  extract_gain_value();
  cout<<"gain values extracted "<<endl;

}
