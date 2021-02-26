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
    float factor;
    double time[6];
    std::string d[5] = { "20000", "10000", "5000", "2500", "1250" };

    merge_tree->SetBranchAddress(   "time", &time   );
    merge_tree->SetBranchAddress(   "iSNR", &iSNR   );
    merge_tree->SetBranchAddress(   "oSNR", &oSNR   );
    merge_tree->SetBranchAddress(  "ioSNR", &ioSNR  );
    merge_tree->SetBranchAddress( "factor", &factor );

    std::vector <TH2F*> VioSNR;
    VioSNR.reserve( numDistances * num_det );

    std::vector <TH2F*> Voverlap;
    Voverlap.reserve( num_det );

    const char* hName;
    std::string histname;
    for( int j = 0; j < num_det; ++j ){
        histname = "Overlap_detector" + std::to_string(j+1);
        hName = histname.c_str();
        Voverlap.push_back( new TH2F(  hName,  hName,
                                      nBinDx, 0, 250,
                                      nBinDy, 0,   1 ));
        for( int i = 0; i < numDistances; ++i ){
            histname = "OverlapDetector" + std::to_string(j+1) + "Distance = " +  d[i] + " kPc";
            hName = histname.c_str();
            VioSNR.push_back(new TH2F(  hName,  hName,
                                       nBinDx, 0, 250,
                                       nBinDy, 0,   1 ));
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
                for( int j = 1; j < 6; ++j ){
                    n = 3*( j - 1 ) + k;
                    if( factor == j )
                        VioSNR[n]->Fill(            sqrt(     iSNR[k]      ),
                                         ioSNR[k] / sqrt( iSNR[k]*oSNR[k] ));
                }
            }
    }

    std::vector <TCanvas*> Canvas_Dist;
    Canvas_Dist.reserve( num_det );
    std::vector <TCanvas*> Canvas;
    Canvas.reserve( num_det );
    std::string namestr;
    const char* namecanvas;
    // creo tre canvas, uno per ogni detector
    for( int j = 0; j < num_det; ++j ){
    	namestr = "DetectorDistance" + std::to_string(j+1);
    	namecanvas = namestr.c_str();
    	Canvas_Dist[j] = new TCanvas( namecanvas, namecanvas, 800, 600 );

    	namestr = "Detector" + std::to_string(j+1);
    	namecanvas = namestr.c_str();
    	Canvas[j] = new TCanvas( namecanvas, namecanvas, 800, 800 );
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
    }

    const char* saveName;
    string saveString;
    for( int j = 0; j < num_det; ++j ){
        saveString = "report/OverlapDistributionDetector" + std::to_string(j+1) + EOS + ".pdf";
        saveName = saveString.c_str();
        Canvas[j]->Draw();
        Canvas[j]->Print((saveName));
        saveString = "report/OverlapDistributionFactorDetector" + std::to_string(j+1) + EOS + ".pdf";
        saveName = saveString.c_str();
        Canvas_Dist[j]->Draw();
        Canvas_Dist[j]->Print((saveName));
    }

}