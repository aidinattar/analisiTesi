#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>


void set_style(){
    // setup graphics
    gStyle -> SetOptTitle   (   kFALSE   );
    gStyle -> SetOptStat    (      0     );
    gStyle -> SetLabelOffset( 0.01 , "x" );
    gStyle -> SetLabelOffset( 0.005, "y" );
    gStyle -> SetTitleOffset( 1.2  , "x" );
    gStyle -> SetTitleOffset( 0.5  , "y" );
}

void PlotFormeOnda( void ){

  TCanvas *c1 = new TCanvas("c1", "My ROOT Plots", 1280, 720);
  c1 -> SetGrid(); //griglia
  set_style();

  gPad   -> SetTopMargin   ( 0.04  );
  gPad   -> SetBottomMargin( 0.09  );
  gPad   -> SetRightMargin ( 0.01  );
  gPad   -> SetLeftMargin  ( 0.05 );

  string title, asx, asy;
  std::cout << "file con la forma d'onda: ";
  std::getline( std::cin, title );
  TGraphErrors *g = new TGraphErrors( title.c_str() );

  std::cout << "asse x: ";
  std::getline( std::cin, asx );
  g -> GetXaxis() -> SetTitle( asx.c_str() );
  g -> GetXaxis() -> CenterTitle();

  std::cout << "asse y: ";
  std::getline( std::cin, asy );
  g -> GetYaxis() -> SetTitle( asy.c_str() );
  g -> GetYaxis() -> CenterTitle();

  g->SetLineColor( kBlue );

  // disegna il grafico a punti (p) in un nuovo frame (a)
  g->Draw("al");

  TLegend *legend = new TLegend();

  legend->AddEntry( g, "experimental data", "l");

  //legend->Draw( "SAME" );

}
