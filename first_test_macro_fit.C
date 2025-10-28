#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TStyle.h"

// Composante Erfc (bord gauche / début de bande)
double erfcStart(double *x, double *par) {
    double A = par[0];
    double x0 = par[1];     // position où l’Erfc commence à descendre
    double sigma = par[2];  // largeur de transition
    return A * TMath::Erfc((x[0]-x0)/sigma);
}

// Composante Gaussienne (fin de bande)
double gaussEnd(double *x, double *par) {
    double A = par[0];
    double mu = par[1];     // centre vers la fin
    double sigma = par[2];  // largeur
    return A * TMath::Exp(-0.5 * TMath::Power((x[0]-mu)/sigma, 2));
}

// Fonction totale : somme des deux composantes sur le même intervalle
double totalFunction(double *x, double *par) {
    double erfcVal = par[0] * TMath::Erfc((x[0]-par[1])/par[2]);
    double gaussVal = par[3] * TMath::Exp(-0.5*TMath::Power((x[0]-par[4])/par[5], 2));
    return erfcVal + gaussVal;
}

void first_test_macro_fit() {
    gStyle->SetOptFit(1111);

    TFile* f = TFile::Open("data_spectrum.root");
    TH1D* h = (TH1D*)f->Get("om_100");
    h->SetLineColor(kBlue);
    h->SetMarkerStyle(20);
    h->Draw();

    // Fonction totale sur le même range
    TF1* ffit = new TF1("ffit", totalFunction, 3000, 13000, 6);
    ffit->SetParName(0,"A_Erfc");
    ffit->SetParName(1,"x0_Erfc");
    ffit->SetParName(2,"sigma_Erfc");
    ffit->SetParName(3,"A_Gauss");
    ffit->SetParName(4,"mu_Gauss");
    ffit->SetParName(5,"sigma_Gauss");

    // Paramètres initiaux réalistes
    ffit->SetParameters(10, 6000, 1000, 20, 12000, 500);

    // Limites pour que l’Erfc ne prenne pas tout et Gauss domine la fin
    ffit->SetParLimits(1, 3000, 9000); // x0_Erfc
    ffit->SetParLimits(4, 9000, 13000); // mu_Gauss

    h->Fit(ffit, "RQ0+");
    ffit->SetLineColor(kRed);
    ffit->Draw("Same");

    // Composante Erfc seule
    TF1* fErfc = new TF1("fErfc", erfcStart, 3000, 13000, 3);
    fErfc->SetParameters(ffit->GetParameter(0), ffit->GetParameter(1),
                         ffit->GetParameter(2));
    fErfc->SetLineColor(kGreen+2);
    fErfc->Draw("Same");

    // Composante Gaussienne seule
    TF1* fGauss = new TF1("fGauss", gaussEnd, 3000, 13000, 3);
    fGauss->SetParameters(ffit->GetParameter(3), ffit->GetParameter(4),
                          ffit->GetParameter(5));
    fGauss->SetLineColor(kMagenta);
    fGauss->Draw("Same");
}
