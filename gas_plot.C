#include <iostream>
#include <vector>
#include <cmath>

#include "Garfield/MediumMagboltz.hh"
#include "TGraph.h"
#include "TCanvas.h"

void gas_plot() {
  using namespace Garfield;

  MediumMagboltz* gas = new MediumMagboltz();
  if (!gas->LoadGasFile("ar_90_ch4_10.gas")) {
    std::cerr << "Error: could not load gas table." << std::endl;
  }

  std::vector<double> redField;    // E/P
  std::vector<double> townsend;   
  std::vector<double> attachment; 
  std::vector<double> fields;
  std::vector<double> vdrift;
  std::vector<double> diffL;
  std::vector<double> diffT;
  
  double bx = 0., by = 0., bz = 0.; // No magnetic field
  double pressureTorr = gas->GetPressure() / 760.; // Torr

 for (double e = 50.; e <= 1000.; e += 50.) {
   double vx, vy, vz;
   if (!gas->ElectronVelocity(e, 0., 0., bx, by, bz, vx, vy, vz)) continue;

   double v = sqrt(vx * vx + vy * vy + vz * vz);

   double dl, dt;
   if (!gas->ElectronDiffusion(e, 0., 0., bx, by, bz, dl, dt)) continue;

   double alpha, eta;
   if (!gas->ElectronTownsend(e, 0., 0., bx, by, bz, alpha)) continue;
   if (!gas->ElectronAttachment(e, 0., 0., bx, by, bz, eta)) continue;

   double ep = e / pressureTorr; 

   fields.push_back(e);
   redField.push_back(ep);
   vdrift.push_back(v);
   diffL.push_back(dl);
   diffT.push_back(dt);
   townsend.push_back(alpha);
   attachment.push_back(eta);

   std::cout << "E = " << e << " V/cm"
	     << "  E/P = " << ep << " V/cm/Torr"
	     << "  vdrift = " << v << " cm/us"
	     << "  DL = " << dl << " cm^0.5"
	     << "  DT = " << dt << " cm^0.5"
	     << "  Alpha = " << alpha << " 1/cm"
	     << "  Eta = " << eta << " 1/cm"
	     << std::endl;
 }

 TCanvas* c1 = new TCanvas("c1", "Transport parameters", 1000, 800);
 c1->Divide(2, 1);

 // Drift velocity
 c1->cd(1);
 TGraph* gV = new TGraph(redField.size(), redField.data(), vdrift.data());
 gV->SetTitle("Drift velocity;E/P [V/cm/Torr];v_{d} [cm/#mus]");
 gV->SetLineColor(kBlue);
 gV->Draw("AL");

 // Diffusion
 c1->cd(2);
 TGraph* gDL = new TGraph(redField.size(), redField.data(), diffL.data());
 gDL->SetTitle("Diffusion coefficients;E/P [V/cm/Torr];Diffusion [cm^{1/2}]");
 gDL->SetLineColor(kRed);
 gDL->Draw("AL");

 TGraph* gDT = new TGraph(redField.size(), redField.data(), diffT.data());
 gDT->SetLineColor(kGreen + 2);
 gDT->Draw("L SAME");

 /*
 // Townsend and attachment
 c1->cd(3);
 TGraph* gAlpha = new TGraph(redField.size(), redField.data(), townsend.data());
 gAlpha->SetTitle("Ionization & Attachment;E/P [V/cm/Torr];Coeff. [1/cm]");
 gAlpha->SetLineColor(kMagenta + 1);
 gAlpha->Draw("AL");

 TGraph* gEta = new TGraph(redField.size(), redField.data(), attachment.data());
 gEta->SetLineColor(kOrange + 2);
 gEta->Draw("L SAME");
 */
 c1->Update();

  
}
