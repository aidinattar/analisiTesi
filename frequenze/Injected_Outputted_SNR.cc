#include <string>
#include <vector>

TString path = "/home/aidin/Documenti/Tesi/analisi/frequenze/";
TString head = "wave_O4I_GN_LHV_SIM_PMNS_";
std::string EOS  = "SHT2_0spin1";
//std::string EOS  = "APR4_q09";
TString tail = "_cl.M1.root";
const int num_det = 3;

void Injected_Outputted_SNR(){
    TString namefile = head + EOS.c_str() + tail;

    TString merge_file_name = path + namefile;
    TString       path_save = path + "report/dump";

    TFile*  filemerge = new TFile( merge_file_name.Data(), "update" );
    TTree* merge_tree =   ( TTree* )    filemerge->Get( "waveburst" );

    int numDistances = 5;
    int nBinD = 100;
    float factor;
    float oSNR[num_det], iSNR[num_det];
    const char* hName;
    std::string histname;
    double time[6];
    std::string d[5] = { "20000", "10000", "5000", "2500", "1250" };

    // associo le variabili ai branch del tree
    merge_tree->SetBranchAddress(   "time", &time   );
    merge_tree->SetBranchAddress( "factor", &factor );
    merge_tree->SetBranchAddress(   "oSNR", &oSNR   );
    merge_tree->SetBranchAddress(   "iSNR", &iSNR   );

    std::vector <TH1F*> VoSNR;
    VoSNR.reserve( numDistances * num_det );
    std::vector <TH1F*> ViSNR;
    ViSNR.reserve( numDistances * num_det );
    // riempio il vettore con un istogramma per ogni distanza e per ogni detector
    for( int i = 0; i < numDistances; ++i )
        for( int j = 0; j < num_det; ++j ){
            histname = "oSNRdistance" + std::to_string(i+1) + "_detector" + std::to_string(j+1);
            hName = histname.c_str();
            VoSNR.push_back(new TH1F( hName, hName, nBinD, 0, 250 ));
            histname = "iSNRdistance" + std::to_string(i+1) + "_detector" + std::to_string(j+1);
            hName = histname.c_str();
            ViSNR.push_back(new TH1F( hName, hName, nBinD, 0, 250 ));
        }

    int n;
    int entr = merge_tree->GetEntries();
    // riempio gli istogrammi
    for( int i = 0; i < entr; ++i ){
        merge_tree->GetEntry( i );
            for( int j = 1; j < 6; ++j )
                if( factor == j )
                    for( int k = 0; k < num_det; ++k )
                        if( abs( time[k] - time[k+3] ) < 5 ){ // condizione di minimo tempo
                            n = 3*( j - 1 ) + k;
//                            std::cout << n << std::endl;
                            VoSNR[n]->Fill( sqrt(oSNR[k]) );
                            ViSNR[n]->Fill( sqrt(iSNR[k]) );
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
    std::vector <TLegend*> LegendiSNR;
    LegendiSNR.reserve( num_det );
    for( int j = 0; j < num_det; ++j ){
    	namestr = "Detector" + std::to_string(j+1);
    	namecanvas = namestr.c_str();
    	LegendiSNR[j] = new TLegend( 0.725, 0.800,
                                     0.950, 0.925 );
        LegendiSNR[j]->SetHeader("Distances", "C");
        LegendiSNR[j]->SetNColumns(2);
    }

    std::vector <TLegend*> LegendoSNR;
    LegendoSNR.reserve( num_det );
    for( int j = 0; j < num_det; ++j ){
    	namestr = "Detector" + std::to_string(j+1);
    	namecanvas = namestr.c_str();
    	LegendoSNR[j] = new TLegend( 0.725, 0.800,
                                     0.950, 0.925 );
        LegendoSNR[j]->SetHeader("Distances", "C");
        LegendoSNR[j]->SetNColumns(2);
    }

    gStyle->SetOptStat(0); // tolgo le statistiche
    for( int j = 0; j < num_det; ++j ){
        Canvas[j]->Divide(1,2);
        Canvas[j]->cd(1);
        gPad->SetGrid();

        for( int i = 0; i < numDistances * num_det; i+=3 ){
            ViSNR[i]->Draw("SAME");
            ViSNR[i]->SetLineWidth(4);
        }
        ViSNR[j   ]->SetLineColor(kBlue);
        ViSNR[j+3 ]->SetLineColor(kRed);
        ViSNR[j+6 ]->SetLineColor(kGreen);
        ViSNR[j+9 ]->SetLineColor(kBlack);
        ViSNR[j+12]->SetLineColor(kViolet);

        ViSNR[j]->            SetTitle( "Injected" );
        ViSNR[j]->GetXaxis()->SetTitle(   "iSNR"   );
        ViSNR[j]->GetYaxis()->SetTitle(  "counts"  );

        ViSNR[j]->GetXaxis()->CenterTitle();
        ViSNR[j]->GetYaxis()->CenterTitle();

        const char* legNameiSNR;
        string lNameiSNR;
        for( int i = 0; i < numDistances * num_det; i+=3 ){
            lNameiSNR   = d[i/3] + " kPc";
            legNameiSNR = lNameiSNR.c_str();
            LegendiSNR[j]->AddEntry( ViSNR[i], legNameiSNR, "l" );
        }

        LegendiSNR[j]->Draw("SAME");
        gStyle->SetGridStyle(7);
        gStyle->SetGridWidth(1);
        gPad->SetLogy();


        Canvas[j]->cd(2);
        gPad->SetGrid();
        for( int i = 0; i < numDistances * num_det; i+=3 ){
            VoSNR[i]->Draw("SAME");
            VoSNR[i]->SetLineWidth(4);
        }
        VoSNR[j   ]->SetLineColor(kBlue);
        VoSNR[j+3 ]->SetLineColor(kRed);
        VoSNR[j+6 ]->SetLineColor(kGreen);
        VoSNR[j+9 ]->SetLineColor(kBlack);
        VoSNR[j+12]->SetLineColor(kViolet);

        VoSNR[j]->            SetTitle( "Injected" );
        VoSNR[j]->GetXaxis()->SetTitle(   "oSNR"   );
        VoSNR[j]->GetYaxis()->SetTitle(  "counts"  );

        VoSNR[j]->GetXaxis()->CenterTitle();
        VoSNR[j]->GetYaxis()->CenterTitle();

        const char* legNameoSNR;
        string lNameoSNR;
        for( int i = 0; i < numDistances * num_det; i+=3 ){
            lNameoSNR   = d[i/3] + " kPc";
            legNameoSNR = lNameoSNR.c_str();
            LegendoSNR[j]->AddEntry( VoSNR[i], legNameoSNR, "l" );
        }

        LegendoSNR[j]->Draw("SAME");
        gStyle->SetGridStyle(7);
        gStyle->SetGridWidth(1);
        gPad->SetLogy();
    }

    const char* saveName;
    string saveString;
    for( int j = 0; j < num_det; ++j ){
        saveString = "report/SNRDistributionFactorsDetector" + std::to_string(j+1) + EOS + ".pdf";
        saveName   = saveString.c_str();
        Canvas[j]->Draw();
        Canvas[j]->Print((saveName));
    }
}