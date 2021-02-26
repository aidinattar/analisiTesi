{
    TString merge_file_name = "mdc_O4I_GN_LHV_SIM_PMNS_SHT2_0spin1.M1.root";
    TFile*  filemerge = new TFile( merge_file_name.Data(), "update" );
    TTree* merge_tree =   ( TTree* )    filemerge->Get(       "mdc" );

    const int num_det = 3;
    float factor;
    double hrss[num_det];
    float  snr[num_det];
    int entr = merge_tree->GetEntries();

    // associo le variabili ai branch del tree
    merge_tree->SetBranchAddress(    "snr", &snr    );
    merge_tree->SetBranchAddress(   "hrss", &hrss   );
    merge_tree->SetBranchAddress( "factor", &factor );

    for( int i = 0; i < entr; ++i ){
      	merge_tree->GetEntry( i );
      	std::cout << factor << "\t" << hrss[] << "\t"
                  << snr[0] << std::endl;
    }
}