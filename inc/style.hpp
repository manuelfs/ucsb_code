#ifndef H_STYLE
#define H_STYLE

#include "TROOT.h"
#include "TColor.h"
#include "TStyle.h"

void SetStyle(){
  const int myColorNum(11235);
  TColor myColor(myColorNum, 1.0, 1.0, 1.0, "myColor", 0.0);

  gROOT->SetStyle("Modern");
  TStyle *tdrStyle = new TStyle(*gStyle);
  tdrStyle->SetNameTitle("tdr","tdr");
  const int numstops=6;
  double stops[numstops]={0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
  double red[numstops]={1.0, 0.0, 0.0, 0.0, 1.0, 1.0};
  double green[numstops]={0.0, 0.0, 1.0, 1.0, 1.0, 0.0};
  double blue[numstops]={1.0, 1.0, 1.0, 0.0, 0.0, 0.0};

  const int numcolors=999;
  TColor::CreateGradientColorTable(numstops, stops, red, green, blue, numcolors);

  tdrStyle->SetNumberContours(999);
  tdrStyle->SetColorModelPS(1);
  tdrStyle->SetTitleFontSize(0.07);
  tdrStyle->SetTitleSize(0.07, "XYZ");
  tdrStyle->SetPadTopMargin(0.15);
  tdrStyle->SetPadBottomMargin(0.15);
  tdrStyle->SetPadLeftMargin(0.15);
  tdrStyle->SetPadRightMargin(0.05);
  tdrStyle->SetLabelSize(0.05, "XYZ");

  tdrStyle->SetLabelFont(62);
  tdrStyle->SetLegendFont(62);
  tdrStyle->SetStatFont(62);
  tdrStyle->SetTextFont(62);
  tdrStyle->SetTitleFont(62);

  tdrStyle->SetFillColor(myColorNum);
  tdrStyle->SetFillStyle(4000);

  tdrStyle->SetCanvasColor(myColorNum);
  tdrStyle->SetPadColor(myColorNum);
  tdrStyle->SetFrameFillColor(myColorNum);
  tdrStyle->SetStatColor(myColorNum);
  tdrStyle->SetTitleFillColor(myColorNum);
  tdrStyle->SetLegendFillColor(myColorNum);

  tdrStyle->SetFrameFillStyle(4000);

  tdrStyle->SetTitleStyle(0);
  tdrStyle->SetCanvasBorderSize(0);
  tdrStyle->SetFrameBorderSize(0);
  tdrStyle->SetLegendBorderSize(0);
  tdrStyle->SetTitleBorderSize(0);

  gROOT->ForceStyle();
  tdrStyle->cd();
  gROOT->ForceStyle();
  tdrStyle->cd();
}

#endif
