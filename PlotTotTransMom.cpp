// ROOT Includes
#include "TH1D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TFile.h"

// Local Includes
#include "starlyze.cpp" 

void PlotTotTransMom(const std::string& result_file_path = "slight.out") {
    // Read inn result
    const SimulationResult results = ReadSimulationResults(result_file_path);

    // Create ROOT output file before any plotting
    const std::string base_file_name = results.decay_repr_str 
                                     + std::string("_") + std::to_string(results.n_events)
                                     + std::string("_") + std::to_string(results.rnd_seed)
                                     + std::string("_tot_trans_mom");
    const std::string root_file_name = base_file_name + std::string(".root");
    TFile* root_file = new TFile(root_file_name.c_str(), "recreate");

    // Create title for plot
    const char* title = Form("\\text{STARlight } | \\text{ Pb - Pb } \\sqrt{s_{NN}} = %.2f \\text{ TeV } | \\, %s", 
                             results.sqrt_s_NN/1000, results.decay_latex_str.c_str());

    // Create list of all invariant masses
    std::vector<double> p_trans_list;
    for (const Event& event : results.events) {
        p_trans_list.push_back(event.p_trans);
    }

    // Calculate histogram properties
    const double min = std::min_element(p_trans_list.begin(), 
                                        p_trans_list.end())[0];
    const double max = std::max_element(p_trans_list.begin(), 
                                        p_trans_list.end())[0];
    const double bin_width = FreedmanDiaconisBinWidth(p_trans_list);
    const int n_bins = (max - min)/bin_width;

    TH1D* hist = new TH1D("hist", title, n_bins, min, max);

    // Fill histograms
    for (const double& p_trans : p_trans_list) {
        hist->Fill(p_trans);
    }

    // Calculate invariant mass peak and its FWHM
    const int bin_max = hist->GetMaximumBin();
    const double hist_peak = hist->GetXaxis()->GetBinCenter(bin_max);
    const double half_max = hist->GetMaximum() / 2.0;
    int fwhm_left = bin_max;
    int fwhm_right = bin_max;
    while (hist->GetBinContent(fwhm_left) > half_max) fwhm_left--;
    while (hist->GetBinContent(fwhm_right) > half_max) fwhm_right++;
    const double fwhm = hist->GetXaxis()->GetBinCenter(fwhm_right) 
                      - hist->GetXaxis()->GetBinCenter(fwhm_left);

    // Text information about amount of events
    const char* events_info = Form("\\text{%i events}", results.n_events);
    TLatex* events_info_text = new TLatex(0.54, 0.80, events_info);
    events_info_text->SetNDC();

    // Text information about inv. mass. peak
    const char* peak_info = Form("\\text{Peak @ %.1f MeV/c}", hist_peak*1000);
    TLatex* peak_info_text = new TLatex(0.54, 0.75, peak_info);
    peak_info_text->SetNDC();

    // Text information about inv. mass. FWHM
    const char* fwhm_info = Form("\\text{FWHM = %.1f MeV/c}", fwhm*1000);
    TLatex* fwhm_info_text = new TLatex(0.54, 0.70, fwhm_info);
    fwhm_info_text->SetNDC();

    // Create a canvas to draw on
    TCanvas* canvas = new TCanvas("canvas", "", 900, 700);

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
    peak_info_text->Draw();
    fwhm_info_text->Draw();

    // Save plot to TEX file
    const std::string tex_file_name = base_file_name + std::string(".tex");
    canvas->Print(tex_file_name.c_str());

    // Save canvas object to ROOT file
    canvas->Write();
}
