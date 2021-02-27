#include <string>
#include <vector>

TString path = "/home/aidin/Documenti/Tesi/analisi/frequenze/";
TString head = "wave_O4I_GN_LHV_SIM_PMNS_";
//std::string EOS  = "SHT2_0spin1";
std::string EOS  = "APR4_q09";
TString tail = "_cl.M1.root";
const int num_det = 3;

void Overlap(){
    TString namefile = head + EOS.c_str() + tail;

    TString merge_file_name = path + namefile;
    TString       path_save = path + "report/dump";

    TFile*  filemerge = new TFile( merge_file_name.Data(), "update" );
    TTree* merge_tree =   ( TTree* )    filemerge->Get( "waveburst" );

    int numDistances = 5;
    int nBinDx = 8*16;
    int nBinDy = 4*75;
    float  iSNR[num_det],
           oSNR[num_det],
          ioSNR[num_det];
    float  high[num_det],
            low[num_det];
    float factor;
    double   time[6];
    std::string d[5] = { "20", "10", "5", "2.5", "1.25" };

    merge_tree->SetBranchAddress(   "time", &time   );
    merge_tree->SetBranchAddress(   "iSNR", &iSNR   );
    merge_tree->SetBranchAddress(   "oSNR", &oSNR   );
    merge_tree->SetBranchAddress(   "high", &high   );
    merge_tree->SetBranchAddress(    "low", &low    );
    merge_tree->SetBranchAddress(  "ioSNR", &ioSNR  );
    merge_tree->SetBranchAddress( "factor", &factor );

    std::vector <TH2F*>                        VioSNR;
    std::vector <TH2F*>                      Voverlap;
    std::vector <TH2F*>                VLowHighFactor;
    std::vector <TH2F*>                      VLowHigh;

    VioSNR.reserve(          numDistances * num_det );
    Voverlap.reserve(                       num_det );
    VLowHighFactor.reserve(  numDistances * num_det );
    VLowHigh.reserve(                       num_det );

    const char* hName;
    std::string histname;
    for( int j = 0; j < num_det; ++j ){
        histname = "Overlap_detector" + std::to_string(j+1);
        hName = histname.c_str();
        Voverlap.push_back( new TH2F(  hName,  hName,
                                      nBinDx, 0, 250,
                                      nBinDy, 0,   1 ));
        histname = "LowHigh_detector" + std::to_string(j+1);
        hName = histname.c_str();
        VLowHigh.push_back( new TH2F(  hName,     hName,
                                       nBinDx, 0,  4000,
                                       nBinDy, 0,  4500 ));
        for( int i = 0; i < numDistances; ++i ){
            histname = "OverlapDetector" + std::to_string(j+1) + "Distance = " +  d[i] + " Mpc";
            hName = histname.c_str();
            VioSNR.push_back(new TH2F(  hName,  hName,
                                       nBinDx, 0, 250,
                                       nBinDy, 0,   1 ));

            histname = "LHDetector" + std::to_string(j+1) + "Distance = " +  d[i] + " Mpc";
            hName = histname.c_str();
            VLowHighFactor.push_back(new TH2F(  hName,  hName,
                                               nBinDx, 0, 3500,
                                               nBinDy, 0, 4500 ));
        }
    }

    int n;
    int entr = merge_tree->GetEntries();
    for( int i = 0; i < entr; ++i ){
        merge_tree->GetEntry( i );
        for( int k = 0; k < num_det; ++k )
            if( abs( time[k] - time[k+3] ) < 5 ){
                Voverlap[k]->Fill(             sqrt(     iSNR[k]     ),
                                   ioSNR[k] /  sqrt( iSNR[k]*oSNR[k] ));
                VLowHigh[k]->Fill(   low[k],                 high[k]  );
                for( int j = 1; j < 6; ++j ){
                    n = 3*( j - 1 ) + k;
                    if( factor == j ){
                        VioSNR[n]->Fill(               sqrt(     iSNR[k]      ),
                                            ioSNR[k] / sqrt( iSNR[k]*oSNR[k] ));
                        VLowHighFactor[n]->Fill(   low[k],            high[k] );
                    }
                }
            }
    }

    std::vector <TCanvas*>   Canvas_Dist;
    std::vector <TCanvas*>        Canvas;
    std::vector <TCanvas*>      CanvasLH;
    std::vector <TCanvas*> CanvasLH_Dist;
    Canvas_Dist.reserve(   num_det );
    Canvas.reserve(        num_det );
    CanvasLH.reserve(      num_det );
    CanvasLH_Dist.reserve( num_det );
    std::string namestr;
    const char* namecanvas;
    // creo tre canvas, uno per ogni detector
    for( int j = 0; j < num_det; ++j ){
    	namestr             =         "DetectorDistance" + std::to_string(j+1);
    	namecanvas          =                                  namestr.c_str();
    	Canvas_Dist[j]      =  new TCanvas( namecanvas, namecanvas, 800, 600 );

    	namestr             =                 "Detector" + std::to_string(j+1);
    	namecanvas          =                                  namestr.c_str();
    	Canvas[j]           =  new TCanvas( namecanvas, namecanvas, 800, 700 );

    	namestr             =               "LHDetector" + std::to_string(j+1);
    	namecanvas          =                                  namestr.c_str();
    	CanvasLH[j]         = new TCanvas( namecanvas, namecanvas, 1000, 700 );

    	namestr             =       "LHDetectorDistance" + std::to_string(j+1);
    	namecanvas          =                                  namestr.c_str();
    	CanvasLH_Dist[j]    = new TCanvas( namecanvas, namecanvas, 1000, 700 );
    }

    gStyle -> SetOptStat(111111);
    gStyle->SetStatX(0.92);
    gStyle->SetStatY(0.4);

    for( int j = 0; j < num_det; ++j ){
        Canvas[j]->cd();
        Voverlap[j]->SetMarkerSize(10);
        Voverlap[j]->Draw("COLZ");

        Voverlap[j]->SetTitle( "Overlap distribution" );
        Voverlap[j]->GetXaxis()->SetTitle( "Injected SNR" );
        Voverlap[j]->GetYaxis()->SetTitle(    "overlap"   );
        Voverlap[j]->GetXaxis()->CenterTitle();
        Voverlap[j]->GetYaxis()->CenterTitle();

        Canvas_Dist[j]->Divide( 3, 2 );
        for( int i = 0; i < 15; i+=3 ){
            Canvas_Dist[j]->cd(( i + 3 )/3);
            VioSNR[i]->Draw("SAME COLZ");
            VioSNR[i]->SetMarkerSize(10);
        }

        for( int i = 0; i < 15; i+=3 ){
            VioSNR[i]->GetXaxis()->SetTitle( "Injected SNR" );
            VioSNR[i]->GetYaxis()->SetTitle(    "overlap"   );
            VioSNR[i]->GetXaxis()->CenterTitle();
            VioSNR[i]->GetYaxis()->CenterTitle();
        }

        CanvasLH[j]->cd();
        gPad->SetGrid();
        VLowHigh[j]->SetMarkerStyle(20);
        VLowHigh[j]->Draw();

        VLowHigh[j]->SetTitle(   "Max e min frequency distribution" );
        VLowHigh[j]->GetXaxis()->SetTitle( "Minimum frequency [Hz]" );
        VLowHigh[j]->GetYaxis()->SetTitle( "Maximum frequency [Hz]" );
        VLowHigh[j]->GetXaxis()->CenterTitle();
        VLowHigh[j]->GetYaxis()->CenterTitle();

        CanvasLH_Dist[j]->Divide( 3, 2 );
        for( int i = 0; i < 15; i+=3 ){
            CanvasLH_Dist[j]->cd(( i + 3 )/3);
            VLowHighFactor[i]->Draw();
            VLowHighFactor[i]->SetMarkerSize(10);

            if (gPad) gPad->SetLeftMargin(0.15);
            if (gPad) gPad->SetBottomMargin(0.15);
        }

        for( int i = 0; i < 15; i+=3 ){
            histname = "Distance = " +  d[i/3] + " Mpc";
            hName = histname.c_str();
            VLowHighFactor[i]->SetTitle(                                hName );
            VLowHighFactor[i]->GetXaxis()->SetTitle( "Minimum frequency [Hz]" );
            VLowHighFactor[i]->GetYaxis()->SetTitle( "Maximum frequency [Hz]" );
            VLowHighFactor[i]->GetXaxis()->CenterTitle();
            VLowHighFactor[i]->GetYaxis()->CenterTitle();
        }

    }

    const char* saveName;
    string saveString;
    for( int j = 0; j < num_det; ++j ){
        saveString = "report/OverlapDistributionDetector"           + std::to_string(j+1) + EOS + ".pdf";
        saveName = saveString.c_str();
        Canvas[j]->Draw();
        Canvas[j]->Print((saveName));

        saveString = "report/OverlapDistributionFactorDetector"     + std::to_string(j+1) + EOS + ".pdf";
        saveName = saveString.c_str();
        Canvas_Dist[j]->Draw();
        Canvas_Dist[j]->Print((saveName));

        saveString = "report/FrequencyLHDistributionDetector"       + std::to_string(j+1) + EOS + ".pdf";
        saveName = saveString.c_str();
        CanvasLH[j]->Draw();
        CanvasLH[j]->Print((saveName));

        saveString = "report/FrequencyLHDistributionFactorDetector" + std::to_string(j+1) + EOS + ".pdf";
        saveName = saveString.c_str();
        CanvasLH_Dist[j]->Draw();
        CanvasLH_Dist[j]->Print((saveName));
    }
}