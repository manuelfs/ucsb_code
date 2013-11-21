#include <vector>
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "pu_constants.hpp"

int main(){
  std::vector<float> dataDist(pu::RunsThrough203002, pu::RunsThrough203002+60);
  std::vector<float> MCDist(pu::Summer2012_S10, pu::Summer2012_S10+60);

  TH1D MC("MC", "Number of Primary Vertices;Primary Vertices; Frequency (a.u.)", 60, 0.5, 60.5);
  TH1D data("data", "Number of Primary Vertices;Primary Vertices; Frequency (a.u.)", 60, 0.5, 60.5);

  for(int i(0); i<60; ++i){
    MC.Fill(i+1.0, MCDist[i]);
    data.Fill(i+1.0, dataDist[i]);
  }
  MC.Scale(1.0/MC.Integral("width"));
  data.Scale(1.0/data.Integral("width"));
  MC.SetLineColor(2);
  data.SetLineColor(1);
  MC.SetStats(0);
  data.SetStats(0);

  TCanvas canvas;
  data.Draw("hist");
  MC.Draw("histsame");
  TLegend legend(0.7, 0.8, 0.9, 0.9);
  legend.AddEntry(&MC, "MC", "lpe");
  legend.AddEntry(&data, "CMS data", "lpe");
  legend.Draw("same");
  canvas.Print("NPVCompare.pdf");
}
