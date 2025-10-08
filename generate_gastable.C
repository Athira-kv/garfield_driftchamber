#include "Garfield/MediumMagboltz.hh"
#include <iostream>

void generate_gastable() {

    
    Garfield::MediumMagboltz* gas = new Garfield::MediumMagboltz();
    gas->SetComposition("Ar", 90., "CH4", 10.);
    gas->SetTemperature(293.15);  
    gas->SetPressure(760.);       

    gas->EnableDrift();
    // gas->EnableTownsend();
    //gas->EnableAttachment();

    const size_t nE = 50;         // number of field points
    const double emin = 100.;      // V/cm
    const double emax = 100000.;  // V/cm
    constexpr bool useLog = true;
    gas->SetFieldGrid(emin, emax, nE, useLog);

    const int ncoll = 5;          // collisions per field point
    gas->GenerateGasTable(ncoll);

    gas->WriteGasFile("ar_90_ch4_10.gas");

    std::cout << "Gas table generated successfully!\n";
}
