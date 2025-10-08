#include <iostream>
#include <vector>
#include <memory>

#include <TApplication.h>
#include <TCanvas.h>
#include <Garfield/MediumMagboltz.hh>
#include <Garfield/Sensor.hh>
#include <Garfield/ViewField.hh>
#include <Garfield/ComponentAnalyticField.hh>
#include <TH1F.h>
#include "Garfield/ViewSignal.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include <TH2D.h>
#include "Garfield/ViewMedium.hh"
#include "Garfield/ViewCell.hh"
#include <TGraph.h>
#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/DriftLineRKF.hh"

using namespace Garfield;

int main(int argc, char *argv[]) {

  TApplication app("app", &argc, argv);

  
  double mm = 1e-3;
  double kilovolt = 1000.0;
  double volt = 1.0;
  
  // geometry paramters
  // units of cm
  
  const int nStrips = 71;
  const double pitch = 0.3; //1.0 * mm spacing;
  const double stripWidth = 0.2;  //2 mm width;
  const double topY = 1.45; //cm 
  const double midY = 0.0;
  const double botY = -1.45;
  const double leftX = -1.0; 
  const double wireAnodeRadius = 0.002; // 20 um 
  const double wireCathodeRadius = 0.01; //100 um
  // top and bottom plane cathode strip potentials

  const double V_first = -0.5 * kilovolt;
  const double V_fourth = 0.0 * volt;
  const double V_last = -4.0 * kilovolt;

  // Middle plane wires potentials:
  const double V_wire_near = -0.5 * kilovolt; // parallel to first strip
  const double V_wire_anode = 1.5 * kilovolt; // parallel to 4th strip (readout)
  const double V_wire_far = -4.0 * kilovolt;  // parallel to last strip


  MediumMagboltz* gas = new MediumMagboltz();
  gas->LoadGasFile("../ar_90_ch4_10.gas");
  gas->EnableDrift();
  
  
  // Create analytic field component
  ComponentAnalyticField* comp = new ComponentAnalyticField();
  comp->SetMedium(gas);

  // assigning potentials along strips
    std::vector<double> potentials(nStrips, 0.);
    for (int i = 0; i <= 3; ++i) {
        double t = double(i) / 3.0;
        potentials[i] = V_first + t * (V_fourth - V_first);
    }
    for (int i = 3; i < nStrips; ++i) {
        double t = double(i - 3) / double(nStrips - 1 - 3);
        potentials[i] = V_fourth + t * (V_last - V_fourth);
    }

    std::cout << "Strip potentials:" << std::endl;
    for (int i = 0; i < nStrips; ++i) {
      std::cout << "Strip " << i << ": " << potentials[i] << " V" << std::endl;
    }


    for (int i = 0; i < nStrips; ++i) {
    double cx = leftX + i * pitch;
    comp->AddWire(cx, topY, stripWidth/2, potentials[i], "topStrip" + std::to_string(i));
    }
    for (int i = 0; i < nStrips; ++i) {
    double cx = leftX + i * pitch;
    comp->AddWire(cx, botY, stripWidth/2, potentials[i], "botStrip" + std::to_string(i));
    }

    
    // middle plane wires
    
    const double wire1_x = leftX;
    const double wire2_x = leftX + 3 * pitch;
    const double wire3_x = leftX + (nStrips - 1) * pitch;

    comp->AddWire(wire1_x, midY, wireCathodeRadius, V_wire_near, "nearCathode");
    comp->AddWire(wire2_x, midY, wireAnodeRadius, V_wire_anode, "Anode");
    comp->AddWire(wire3_x, midY, wireCathodeRadius, V_wire_far, "farCathode");

    comp->AddReadout("Anode");

  TCanvas* c1 = new TCanvas("c1", "Drift Chamber Geometry", 800, 600);

    comp->PlotCell(c1);

   ViewField* fieldView = new ViewField();
  fieldView->SetComponent(comp);
  fieldView->SetNumberOfContours(50);
  double xminPlot = -2.0; //  leftX - 2 * pitch;
  double xmaxPlot = 20.0 ; //leftX + (nStrips - 1) * pitch + 2 * pitch;
  double yminPlot = botY - 15.;
  double ymaxPlot = topY + 15.;

  fieldView->SetArea(xminPlot, yminPlot, xmaxPlot, ymaxPlot);
  fieldView->PlotContour("V");
  

   
    Sensor* sensor = new Sensor();
    sensor->AddComponent(comp);
    sensor->AddElectrode(comp, "Anode");
    const double tstep = 0.5;
    const double tmin = -0.5 * tstep;
    const unsigned int nbins = 1000;
    sensor->SetTimeWindow(tmin, tstep, nbins);

    // drift lines from a track
    TrackHeed track(sensor);
    track.SetParticle("muon");
    track.SetEnergy(1.e9);
  
    DriftLineRKF drift(sensor);
    const double x0 = 2.;
    const double y0 = 0.5;
    track.NewTrack(x0, y0, 0, 0, 0, 1, 0);

    ViewDrift driftView;
    driftView.SetArea(-5, -5, 5, 5);
    drift.EnablePlotting(&driftView);
    track.EnablePlotting(&driftView);
    

    // Loop over the clusters along the track.
    std::cout << "Number of clusters: " << track.GetClusters().size() << "\n";
    for (const auto& cluster : track.GetClusters()) {
      // Loop over the electrons in the cluster.
      //std::cout<<" electrons associ. w clusters = "<<cluster.electrons.size()<<"\n";
      for (const auto& electron : cluster.electrons) {
	drift.DriftElectron(electron.x, electron.y, electron.z, electron.t);
      }
    }
    
    
   
    TCanvas *c2 = new TCanvas("c2", "", 600, 600);

    ViewCell cellView(comp);
    cellView.SetCanvas(c2);
    cellView.Plot2d();
    driftView.SetCanvas(c2);
    constexpr bool twod = true;
    constexpr bool drawaxis = false;
    driftView.Plot(twod, drawaxis);
    
    TCanvas *cS = new TCanvas("cS", "", 600, 600);
    sensor->PlotSignal("Anode", cS);



  std::cout << "Simulation finished.\n";

  
  app.Run(true);

}
