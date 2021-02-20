
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>

int main( int argc, char* argv[] ) {

    TFile *file = new TFile("mdc_O4I_GN_LHV_SIM_PMNS_SHT2_0spin1.M1.root");
    file->cd();

    // save current working area
    TDirectory* currentDir = gDirectory;
    std::string filename = "results";
    TFile *results = new TFile( filename.c_str(), "RECREATE" );

    int numDistances = 5;
    int nBinD = 128;
    float factor, snr;
    double hrss;
    const char* hName;
    std::string histname;

    TTree *mdc = ( TTree* )file->Get( "mdc" );

    mdc->SetBranchAddress( "factor", &factor );
    mdc->SetBranchAddress(   "hrss", &hrss   );
    mdc->SetBranchAddress(    "snr", &snr    );

    std::vector <TH1F*> Vhrss;
    Vhrss.reserve( numDistances );
    for( int i = 0; i < numDistances; ++i ){
        histname = "distance" + std::to_string(i);
        hName = histname.c_str();
        Vhrss.push_back(new TH1F( hName, hName, nBinD, 0, 2.6e-21 ));
    }

    int entr = mdc->GetEntries();
    for( int i = 0; i < entr; ++i ){
        mdc->GetEntry( i );
        for( int j = 1; j < 6; ++j )
            if( factor == j )
                Vhrss[j]->Fill( hrss );
    }

    TCanvas  *ct8 = new TCanvas( "ct8", "ct8", 800, 800);
    ct8->Divide(5,1);
    for( int i = 1; i < 6; ++i ){
        ct8->cd(i);
        Vhrss[i]->Write();
    }

    results->Close();
    delete results;
    // restore working area
    currentDir->cd();
//    ct8->Draw();
    return 0;

}