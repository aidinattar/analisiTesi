const int num_det = 3;

TString path      = "/home/aidin/Documenti/Tesi/frequenze/";
TString filemerge = "wave_O4I_GN_LHV_SIM_PMNS_SHT2_0spin1_cl.M1.root";




void macro_test_read(){

  TString merge_file_name = path + filemerge;

  TFile*  filemerge = new TFile( merge_file_name.Data(), "update" );
  TTree* merge_tree =   ( TTree* )    filemerge->Get( "waveburst" );

  TH1F* hist_rho = new TH1F( "hist_rho", "hist_rho", 100, 0, 100 );
  TH2F* hist_rho_time[num_det];

  char cuts[4096];
  char name[1028];
  char  sel[1028];

  TCanvas** cx = new TCanvas*[num_det];
  std::string namestr;
  const char* namecanvas;
  for( int j = 0; j < num_det; ++j ){
	namestr = "c" + std::to_string(j);
	namecanvas = namestr.c_str();
	cx[j] = new TCanvas( namecanvas, namecanvas, 400, 400 );
  }

  for( int jj=0; jj<num_det; jj++ ){ // ciclo sui rivelatori
	cx[jj]->cd();
   	sprintf( cuts, "rho<250" );
  	merge_tree->Draw( "rho>>hist_rho", cuts);

  	sprintf( name, "hist_rho_time_%i", jj );
  	hist_rho_time[jj] = new TH2F( name,  name, 1000, 1.1645e+9,
                                      1.1651*1e+9,  100, 0,   250);
  	sprintf( sel, "rho:time[%i]>>hist_rho_time_%i", jj, jj );
	merge_tree->Draw( sel, cuts );
	cx[jj]->Draw("COL");
  }

  int            entr = merge_tree->GetEntries();
  const auto det_cons =		 2 * 	 num_det;
  double time[det_cons];
  float rho;

  merge_tree->SetBranchAddress( "time", &time );
  merge_tree->SetBranchAddress(  "rho", &rho  );

  for( int i=0; i < entr; i++ ){
  	merge_tree->GetEntry( i );
//  	cout << time[0] << "\t" << rho << endl;
  }

  TCanvas* Cplot = new TCanvas( "plot", "plot", 600,
                                   400,    600, 400);
  Cplot->Divide( 2, 1 );
  Cplot->cd(        1 );

  hist_rho->Draw();
  Cplot->cd(2);
  hist_rho_time[0]->Draw();
  char out[1024];

  sprintf( out, "%s/test.png", path_save.Data() );
  Cplot->SaveAs( out );

}
