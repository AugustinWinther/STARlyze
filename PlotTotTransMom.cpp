// ROOT Includes
#include "TH1D.h"
#include "TCanvas.h"
#include "TLatex.h"

// Local Includes
#include "starlyze.cpp" 

void PlotTotTransMom(const std::string& result_file_path = "slight.out") {
    // Read inn result
    const SimulationResult results = ReadSimulationResults(result_file_path);

    // Create title for plot
    char* title = Form("\\text{STARlight } | \\text{ Pb - Pb } \\sqrt{s_{NN}} = %.2f \\text{ TeV } | \\, %s", 
                       results.sqrt_s_NN/1000, results.decay_latex_str.c_str());

    // Calculate histogram properties
    const double min = std::min_element(results.p_trans_list.begin(), 
                                        results.p_trans_list.end())[0];
    const double max = std::max_element(results.p_trans_list.begin(), 
                                        results.p_trans_list.end())[0];
    const double bin_width = FreedmanDiaconisBinWidth(results.p_trans_list);
    const int n_bins = (max - min)/bin_width;

    // Create histogram object
    TH1D *hist = new TH1D("hist", title, n_bins, min, max);

    // Fill histograms
    for (const double& p_trans : results.p_trans_list) {
        hist->Fill(p_trans);
    }

    // Text information about amount of events
    const char *events_info = Form("\\text{%i events}", results.n_events);
    TLatex *events_info_text = new TLatex(0.70, 0.65, events_info);
    events_info_text->SetNDC();
    events_info_text->SetTextAlign(0);
    events_info_text->SetTextSize(0.04);

    // Create a canvas to draw on
    TCanvas *canvas = new TCanvas("canvas", "", 900, 700);

    // Draw histograms and info texts
    hist->SetStats(kFALSE);
    hist->SetXTitle("\\text{4 Particle Transverse Momentum  [GeV/c]}");
    hist->SetYTitle(Form("\\text{Counts per %.2f [MeV/c]}", bin_width*1000));
    hist->SetAxisRange(0, 0.25);
    hist->SetAxisRange(0, hist->GetBinContent(hist->GetMaximumBin())*1.1, "Y");
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
    hist->GetXaxis()->SetTitleOffset(1.0);
    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->GetXaxis()->SetLabelSize(0.035);
    hist->GetYaxis()->SetLabelSize(0.04);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->SetLineColor(kBlack);
    hist->SetFillColor(kP10Blue);
    hist->Draw();
    events_info_text->Draw();

    // Save plot to TEX file
    std::string file_name = results.repr_str + std::string("_tot_trans_mom.tex");
    canvas->Print(file_name.c_str());
}
