#include <string>
#include <vector>

TString path = "/home/aidin/Documenti/Tesi/analisi/frequenze/";
TString head = "wave_O4I_GN_LHV_SIM_PMNS_";
std::string EOS  = "SHT2_0spin1";
//std::string EOS  = "APR4_q09";
TString tail = "_cl.M1.root";
const int num_det = 3;

void Overlap_Distribution(){
    TString namefile = head + EOS.c_str() + tail;

    TString merge_file_name = path + namefile;
    TString       path_save = path + "report/dump";

    TFile*  filemerge = new TFile( merge_file_name.Data(), "update" );
    TTree* merge_tree =   ( TTree* )    filemerge->Get( "waveburst" );

    const int nBinDx = 15;
    const int nBinDy = 35;
    const int larg   = 10;
    float factor;
    float  iSNR[num_det],
           oSNR[num_det],
          ioSNR[num_det];
    double time[6];
    double    **SNR = new double*[num_det];
    double **errSNR = new double*[num_det]; //HA ERRORE L'SNR RICOSTRUITO?
    double    **mea = new double*[num_det];
    double    **rms = new double*[num_det];
    double    **err = new double*[num_det];

    for( int i = 0; i < num_det; ++i ){
        SNR[i]    = new double[nBinDx];
        errSNR[i] = new double[nBinDx];
        rms[i]    = new double[nBinDx];
        mea[i]    = new double[nBinDx];
        err[i]    = new double[nBinDx];
    }

    merge_tree->SetBranchAddress(   "time", &time   );
    merge_tree->SetBranchAddress(   "iSNR", &iSNR   );
    merge_tree->SetBranchAddress(   "oSNR", &oSNR   );
    merge_tree->SetBranchAddress(  "ioSNR", &ioSNR  );
    merge_tree->SetBranchAddress( "factor", &factor );

    std::vector <std::vector <TH1F*>> V(num_det);

    const char* hName;
    std::string histname;
    for( int j = 0; j < num_det; ++j )
        for( int i = 0; i < nBinDx; ++i ){
            SNR[j][i] = ( i * larg + ( i + 1 ) * larg ) / 2;
            histname = "Det" + to_string(j+1) + "SNR" +  to_string( (int)SNR[j][i] );
            hName = histname.c_str();
            V[j].push_back(new TH1F( hName, hName, nBinDy,
                                     0.5+0.5*(1-exp(-0.2*(i-0.2))), 1 ));
        }

    int entr = merge_tree->GetEntries();
    for( int i = 0; i < entr; ++i ){
        merge_tree->GetEntry( i );
        for( int k = 0; k < num_det; ++k )
            if( abs( time[k] - time[k+3] ) < 5 )
                for( int j = 0; j < nBinDx; ++j )
                    if( sqrt(oSNR[k]) >=    j     * larg &&
                        sqrt(oSNR[k]) < ( j + 1 ) * larg )
                        V[k][j]->Fill( ioSNR[k] / sqrt( iSNR[k]*oSNR[k] ));
    }

    for( int j = 0; j < num_det; ++j )
        for( int i = 0; i < nBinDx; ++i ){
            mea[j][i]    = V[j][i]->GetMean();
            rms[j][i]    = V[j][i]->GetStdDev();
            errSNR[j][i] = 1 / sqrt( SNR[j][i] );
        }

    std::vector <TCanvas*> Canvas;
    Canvas.reserve( num_det );
    std::string namestr;
    const char* namecanvas;
    // creo tre canvas, uno per ogni detector
    for( int j = 0; j < num_det; ++j ){
    	namestr = "Detector" + std::to_string(j+1);
    	namecanvas = namestr.c_str();
    	Canvas[j] = new TCanvas( namecanvas, namecanvas, 800, 600 );
    }

    gStyle->SetStatX(0.25);
    gStyle->SetStatY(0.925);
    const char* saveName;
    string saveString;
    for( int j = 0; j < num_det; ++j ){
        Canvas[j]->Divide( 5, 3 );
        for( int i = 0; i < nBinDx; ++i ){
            Canvas[j]->cd(i+1);
            V[j][i]->Draw("");
        }
        saveString = "report/OverlapDistributionsDetector" + std::to_string(j+1) + EOS + ".pdf";
        saveName   = saveString.c_str();
        Canvas[j]->Draw();
        Canvas[j]->Print((saveName));
    }


    std::vector <TCanvas*> Canvas_Tot;
    Canvas_Tot.reserve( num_det );
    // creo tre canvas, uno per ogni detector
    for( int j = 0; j < num_det; ++j ){
    	namestr = "DistributionDetector" + std::to_string(j+1);
    	namecanvas = namestr.c_str();
    	Canvas_Tot[j] = new TCanvas( namecanvas, namecanvas, 800, 600 );
    }

    std::vector <TGraphErrors*> overlaps(num_det);
    for( int j = 0; j < num_det; ++j ){
        overlaps[j] = new TGraphErrors( nBinDx, SNR[j], mea[j],
                                             errSNR[j], rms[j] );
        overlaps[j]->GetXaxis()->SetTitle( "Outputted SNR" );
        overlaps[j]->GetXaxis()->CenterTitle();
        overlaps[j]->GetYaxis()->SetTitle( "averge overlap" );
        overlaps[j]->GetYaxis()->CenterTitle();
        overlaps[j]->SetTitle("Overlaps distribution");
        overlaps[j]->SetMarkerStyle(104);
        gPad->SetGrid();
        Canvas_Tot[j]->cd();
        overlaps[j]->Draw("AP");
    }

    for( int j = 0; j < num_det; ++j ){
        saveString = "report/Overlaps_GrapgDetector" + std::to_string(j+1) + EOS + ".pdf";
        saveName   = saveString.c_str();
        Canvas_Tot[j]->Draw();
        Canvas_Tot[j]->Print((saveName));
    }
}