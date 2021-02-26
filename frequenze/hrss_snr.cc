#include <string>
#include <vector>
#include <iostream>
#include <cmath>

TString path = "/home/aidin/Documenti/Tesi/analisi/frequenze/";
TString head = "mdc_O4I_GN_LHV_SIM_PMNS_";
//std::string EOS  = "SHT2_0spin1";
std::string EOS  = "APR4_q09";
TString tail = ".M1.root";

void hrss_snr(){
    TString namefile = head + EOS.c_str() + tail;
    TString merge_file_name = path + namefile;

    TFile*  filemerge = new TFile( merge_file_name.Data(), "update" );
    TTree* merge_tree =   ( TTree* )    filemerge->Get(       "mdc" );

    const int num_det = 3;
    float factor;
    double hrss[num_det];
    float   snr[num_det];
    int entr = merge_tree->GetEntries();

    const int numDistances = 5;
    const int nBinD = 75;
    const char* hName;
    std::string histname;
    std::string d[numDistances] = { "20000", "10000", "5000", "2500", "1250" };
    int n;

    std::vector <TH1F*> Vhrss;
    Vhrss.reserve( numDistances * num_det );
    std::vector <TH1F*>  Vsnr;
    Vsnr.reserve(  numDistances * num_det );
    // riempio il vettore con un istogramma per ogni distanza e per ogni detector
    for( int i = 0; i < numDistances; ++i )
        for( int j = 0; j < num_det; ++j ){
            histname = "hrssdistance" + std::to_string(i+1) + "_detector" + std::to_string(j+1);
            hName = histname.c_str();
            Vhrss.push_back(new TH1F( hName, hName, nBinD, 0, 1.6e-21 ));
            histname = "snrdistance" + std::to_string(i+1) + "_detector" + std::to_string(j+1);
            hName = histname.c_str();
            Vsnr.push_back(new TH1F(  hName, hName, nBinD, 0,     350 ));
        }

    // associo le variabili ai branch del tree
    merge_tree->SetBranchAddress(    "snr", &snr    );
    merge_tree->SetBranchAddress( "factor", &factor );
    merge_tree->SetBranchAddress(   "hrss", &hrss   );

    // riempio gli istogrammi
    for( int i = 0; i < entr; ++i ){
        merge_tree->GetEntry( i );
            for( int j = 1; j < 6; ++j )
                if( factor == j )
                    for( int k = 0; k < num_det; ++k ){
                        n = 3*( j - 1 ) + k;
                        Vhrss[n]->Fill( hrss[k] );
                        Vsnr[n]->Fill(   snr[k] );
                    }
    }

    std::vector <TCanvas*> Canvas;
    Canvas.reserve( num_det );
    std::string namestr;
    const char* namecanvas;
    // creo tre canvas, uno per ogni detector
    for( int j = 0; j < num_det; ++j ){
    	namestr = "Detector" + std::to_string(j+1);
    	namecanvas = namestr.c_str();
    	Canvas[j] = new TCanvas( namecanvas, namecanvas, 400, 800 );
    }

    std::string nameLeg;
    const char* nameLegend;
    // creo gli oggetti legenda
    std::vector <TLegend*> Legendsnr;
    Legendsnr.reserve( num_det );
    for( int j = 0; j < num_det; ++j ){
    	namestr = "Detector" + std::to_string(j+1);
    	namecanvas = namestr.c_str();
    	Legendsnr[j] = new TLegend( 0.725, 0.800,
                                    0.950, 0.925 );
        Legendsnr[j]->SetHeader("Distances", "C");
        Legendsnr[j]->SetNColumns(2);
    }

    std::vector <TLegend*> Legendhrss;
    Legendhrss.reserve( num_det );
    for( int j = 0; j < num_det; ++j ){
    	namestr = "Detector" + std::to_string(j+1);
    	namecanvas = namestr.c_str();
    	Legendhrss[j] = new TLegend( 0.725, 0.800,
                                     0.950, 0.925 );
        Legendhrss[j]->SetHeader("Distances", "C");
        Legendhrss[j]->SetNColumns(2);
    }

    gStyle->SetOptStat(0); // tolgo le statistiche
    for( int j = 0; j < num_det; ++j ){
        Canvas[j]->Divide(1,2);
        Canvas[j]->cd(1);
        gPad->SetGrid();

        for( int i = 0; i < numDistances * num_det; i+=3 ){
            Vsnr[i]->Draw("SAME");
            Vsnr[i]->SetLineWidth(4);
        }
        Vsnr[j   ]->SetLineColor(kBlue);
        Vsnr[j+3 ]->SetLineColor(kRed);
        Vsnr[j+6 ]->SetLineColor(kGreen);
        Vsnr[j+9 ]->SetLineColor(kBlack);
        Vsnr[j+12]->SetLineColor(kViolet);

        Vsnr[j]->            SetTitle( "simulation snr" );
        Vsnr[j]->GetXaxis()->SetTitle(      "snr"       );
        Vsnr[j]->GetYaxis()->SetTitle(     "counts"     );

        Vsnr[j]->GetXaxis()->CenterTitle();
        Vsnr[j]->GetYaxis()->CenterTitle();

        const char* legNamesnr;
        string lNamesnr;
        for( int i = 0; i < numDistances * num_det; i+=3 ){
            lNamesnr   = d[i/3] + " kPc";
            legNamesnr = lNamesnr.c_str();
            Legendsnr[j]->AddEntry( Vsnr[i], legNamesnr, "l" );
        }

        Legendsnr[j]->Draw("SAME");
        gStyle->SetGridStyle(7);
        gStyle->SetGridWidth(1);
        gPad->SetLogy();


        Canvas[j]->cd(2);
        gPad->SetGrid();
        for( int i = 0; i < numDistances * num_det; i+=3 ){
            Vhrss[i]->Draw("SAME");
            Vhrss[i]->SetLineWidth(4);
        }
        Vhrss[j   ]->SetLineColor(kBlue);
        Vhrss[j+3 ]->SetLineColor(kRed);
        Vhrss[j+6 ]->SetLineColor(kGreen);
        Vhrss[j+9 ]->SetLineColor(kBlack);
        Vhrss[j+12]->SetLineColor(kViolet);

        Vhrss[j]->            SetTitle( "simulation hrss" );
        Vhrss[j]->GetXaxis()->SetTitle(   "hrss"   );
        Vhrss[j]->GetYaxis()->SetTitle(  "counts"  );

        Vhrss[j]->GetXaxis()->CenterTitle();
        Vhrss[j]->GetYaxis()->CenterTitle();

        const char* legNamehrss;
        string lNamehrss;
        for( int i = 0; i < numDistances * num_det; i+=3 ){
            lNamehrss   = d[i/3] + " kPc";
            legNamehrss = lNamehrss.c_str();
            Legendhrss[j]->AddEntry( Vhrss[i], legNamehrss, "l" );
        }

        Legendhrss[j]->Draw("SAME");
        gStyle->SetGridStyle(7);
        gStyle->SetGridWidth(1);
        gPad->SetLogy();
    }

    const char* saveName;
    string saveString;
    for( int j = 0; j < num_det; ++j ){
        saveString = "report/hrss_snr_DistributionFactorsDetector" + std::to_string(j+1) + EOS + ".pdf";
        saveName   = saveString.c_str();
        Canvas[j]->Draw();
        Canvas[j]->Print((saveName));
    }
}