#include <string>
#include <vector>

TString path = "/home/aidin/Documenti/Tesi/analisi/frequenze/";
TString head = "wave_O4I_GN_LHV_SIM_PMNS_";
//std::string EOS  = "SHT2_0spin1";
std::string EOS  = "APR4_q09";
TString tail = "_cl.M1.root";
const int num_det = 3;

void Frequency_Duration(){
    TString namefile = head + EOS.c_str() + tail;

    TString merge_file_name = path + namefile;
    TString       path_save = path + "report/dump";

    TFile*  filemerge = new TFile( merge_file_name.Data(), "update" );
    TTree* merge_tree =   ( TTree* )    filemerge->Get( "waveburst" );

    int numDistances = 5;
    int nBinD = 50;
    float factor;
    float duration[num_det], frequency[num_det];
    const char* hName;
    std::string histname;
    double time[6];
    std::string d[5] = { "20", "10", "5", "2.5", "1.25" };
    int n;
    int entr = merge_tree->GetEntries();
    double maxtime = 0.1;
    double maxfreq = 3500;

    // associo le variabili ai branch del tree
    merge_tree->SetBranchAddress(   "time", &time   );
    merge_tree->SetBranchAddress( "factor", &factor );
    merge_tree->SetBranchAddress(   "duration", &duration   );
    merge_tree->SetBranchAddress(   "frequency", &frequency   );

    std::vector <TH1F*> Vduration;
    Vduration.reserve( numDistances * num_det );
    std::vector <TH1F*> Vfrequency;
    Vfrequency.reserve( numDistances * num_det );
    // riempio il vettore con un istogramma per ogni distanza e per ogni detector
    for( int i = 0; i < numDistances; ++i )
        for( int j = 0; j < num_det; ++j ){
            histname = "durationdistance" + std::to_string(i+1) + "_detector" + std::to_string(j+1);
            hName = histname.c_str();
            Vduration.push_back(new TH1F( hName, hName, nBinD, 0, maxtime ));
            histname = "frequencydistance" + std::to_string(i+1) + "_detector" + std::to_string(j+1);
            hName = histname.c_str();
            Vfrequency.push_back(new TH1F( hName, hName, nBinD, 0, maxfreq ));
        }

    // riempio gli istogrammi
    for( int i = 0; i < entr; ++i ){
        merge_tree->GetEntry( i );
        for( int j = 1; j < numDistances+1; ++j )
            if( factor == j )
                for( int k = 0; k < num_det; ++k )
                    if( abs( time[k] - time[k+3] ) < 5 ){ // condizione di minimo tempo
                        n = 3*( j - 1 ) + k;
//                        std::cout << n << std::endl;
                        Vduration[n]->Fill( duration[k] );
                        Vfrequency[n]->Fill( frequency[k] );
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
    std::vector <TLegend*> Legendfrequency;
    Legendfrequency.reserve( num_det );
    for( int j = 0; j < num_det; ++j ){
    	namestr = "Detector" + std::to_string(j+1);
    	namecanvas = namestr.c_str();
    	Legendfrequency[j] = new TLegend( 0.725, 0.800,
                                     0.950, 0.925 );
        Legendfrequency[j]->SetHeader("Distances", "C");
        Legendfrequency[j]->SetNColumns(2);
    }

    std::vector <TLegend*> Legendduration;
    Legendduration.reserve( num_det );
    for( int j = 0; j < num_det; ++j ){
    	namestr = "Detector" + std::to_string(j+1);
    	namecanvas = namestr.c_str();
    	Legendduration[j] = new TLegend( 0.725, 0.800,
                                     0.950, 0.925 );
        Legendduration[j]->SetHeader("Distances", "C");
        Legendduration[j]->SetNColumns(2);
    }

    gStyle->SetOptStat(0); // tolgo le statistiche
    for( int j = 0; j < num_det; ++j ){
        Canvas[j]->Divide(1,2);
        Canvas[j]->cd(1);
        gPad->SetGrid();

        for( int i = 0; i < numDistances * num_det; i+=3 ){
            Vfrequency[i]->Draw("SAME");
            Vfrequency[i]->SetLineWidth(4);
        }
        Vfrequency[j   ]->SetLineColor(kBlue);
        Vfrequency[j+3 ]->SetLineColor(kRed);
        Vfrequency[j+6 ]->SetLineColor(kGreen);
        Vfrequency[j+9 ]->SetLineColor(kBlack);
        Vfrequency[j+12]->SetLineColor(kViolet);

        Vfrequency[j]->            SetTitle( "Injected" );
        Vfrequency[j]->GetXaxis()->SetTitle(   "frequency"   );
        Vfrequency[j]->GetYaxis()->SetTitle(  "counts"  );

        Vfrequency[j]->GetXaxis()->CenterTitle();
        Vfrequency[j]->GetYaxis()->CenterTitle();

        const char* legNamefrequency;
        string lNamefrequency;
        for( int i = 0; i < numDistances * num_det; i+=3 ){
            lNamefrequency   = d[i/3] + " Mpc";
            legNamefrequency = lNamefrequency.c_str();
            Legendfrequency[j]->AddEntry( Vfrequency[i], legNamefrequency, "l" );
        }

        Legendfrequency[j]->Draw("SAME");
        gStyle->SetGridStyle(7);
        gStyle->SetGridWidth(1);
        gPad->SetLogy();


        Canvas[j]->cd(2);
        gPad->SetGrid();
        for( int i = 0; i < numDistances * num_det; i+=3 ){
            Vduration[i]->Draw("SAME");
            Vduration[i]->SetLineWidth(4);
        }
        Vduration[j   ]->SetLineColor(kBlue);
        Vduration[j+3 ]->SetLineColor(kRed);
        Vduration[j+6 ]->SetLineColor(kGreen);
        Vduration[j+9 ]->SetLineColor(kBlack);
        Vduration[j+12]->SetLineColor(kViolet);

        Vduration[j]->            SetTitle( "Outputted" );
        Vduration[j]->GetXaxis()->SetTitle(   "duration"    );
        Vduration[j]->GetYaxis()->SetTitle(  "counts"   );

        Vduration[j]->GetXaxis()->CenterTitle();
        Vduration[j]->GetYaxis()->CenterTitle();

        const char* legNameduration;
        std::string lNameduration;
        for( int i = 0; i < numDistances * num_det; i+=3 ){
            lNameduration   = d[i/3] + " Mpc";
            legNameduration = lNameduration.c_str();
            Legendduration[j]->AddEntry( Vduration[i], legNameduration, "l" );
        }

        Legendduration[j]->Draw("SAME");
        gStyle->SetGridStyle(7);
        gStyle->SetGridWidth(1);
        gPad->SetLogy();
    }

    const char* saveName;
    string saveString;
    for( int j = 0; j < num_det; ++j ){
        saveString = "report/FrequencyDistributionFactorsDetector" + std::to_string(j+1) + EOS + ".pdf";
        saveName   = saveString.c_str();
        Canvas[j]->Draw();
        Canvas[j]->Print((saveName));
    }
}