// ROOT Includes
#include "TH1D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TFile.h"

// Local Includes
#include "starlyze.cpp" 

void PlotTotInvMass(const std::string& result_file_path = "slight.out") {
    // Read inn result
    const SimulationResult results = ReadSimulationResults(result_file_path);

    // Create ROOT output file before any plotting
    const std::string base_file_name = results.decay_repr_str 
                                     + std::string("_") + std::to_string(results.n_events)
                                     + std::string("_") + std::to_string(results.rnd_seed)
                                     + std::string("_tot_inv_mass");
    const std::string root_file_name = base_file_name + std::string(".root");
    TFile* root_file = new TFile(root_file_name.c_str(), "recreate");

    // Create title for plot
    const char* title = Form("\\text{STARlight } | \\text{ Pb - Pb } \\sqrt{s_{NN}} = %.2f \\text{ TeV } | \\, %s", 
                             results.sqrt_s_NN/1000, results.decay_latex_str.c_str());

    // Create list of all invariant masses
    std::vector<double> m_inv_list;
    for (const Event& event : results.events) {
        m_inv_list.push_back(event.m_inv);
    }

    // Calculate histogram properties
    const double min = std::min_element(m_inv_list.begin(), 
                                        m_inv_list.end())[0];
    const double max = std::max_element(m_inv_list.begin(), 
                                        m_inv_list.end())[0];
    const double bin_width = FreedmanDiaconisBinWidth(m_inv_list);
    const int n_bins = (max - min)/bin_width;

    // Create histogram object
    TH1D *hist = new TH1D("hist", title, n_bins, min, max);

    // Fill histograms
    for (const double& m_inv : m_inv_list) {
        hist->Fill(m_inv);
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
    const char* peak_info = Form("\\text{Peak @ %.4f GeV/c}^{2}", hist_peak);
    TLatex* peak_info_text = new TLatex(0.54, 0.75, peak_info);
    peak_info_text->SetNDC();

    // Text information about inv. mass. FWHM
    const char* fwhm_info = Form("\\text{FWHM = %.3f keV/c}^{2}", fwhm*1000000);
    TLatex* fwhm_info_text = new TLatex(0.54, 0.70, fwhm_info);
    fwhm_info_text->SetNDC();

    // Create a canvas to draw on
    TCanvas* canvas = new TCanvas("canvas", "", 900, 700);

    // Draw histograms and info texts
    hist->SetStats(kFALSE);
    hist->SetXTitle("\\text{4 Particle Invariant Mass [GeV/c}^{2}\\text{]}");
    hist->SetYTitle(Form("\\text{Counts per %.2f [keV/c}^{2}\\text{]}", bin_width*1000000));
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
    const std::string file_name = base_file_name + std::string(".tex");
    canvas->Print(file_name.c_str());

    // Save canvas object to ROOT file
    canvas->Write();
}


