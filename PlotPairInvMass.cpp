// ROOT Includes
#include "TH1D.h"
#include "TCanvas.h"
#include "TLatex.h"

// Local Includes
#include "starlyze.cpp" 

void PlotPairInvMass(const std::string& result_file_path = "slight.out") {
    // Read inn result
    const SimulationResult results = ReadSimulationResults(result_file_path);

    // Create title for plot
    const char* title = Form("\\text{STARlight } | \\text{ Pb - Pb } \\sqrt{s_{NN}} = %.2f \\text{ TeV } | \\, %s", 
                             results.sqrt_s_NN/1000, results.decay_latex_str.c_str());

    // Collapse m_inv_pair_list to 1D vector
    std::vector<double> inv_masses;
    for (const std::vector<double>& m_inv_pairs : results.m_inv_pairs_list) {
        for (const double& m_inv_pair : m_inv_pairs) {
            inv_masses.push_back(m_inv_pair);
        }
    }

    // Calculate histogram properties
    const double min = std::min_element(inv_masses.begin(), 
                                        inv_masses.end())[0];
    const double max = std::max_element(inv_masses.begin(), 
                                        inv_masses.end())[0];
    const double bin_width = FreedmanDiaconisBinWidth(inv_masses);
    const int n_bins = (max - min)/bin_width;

    // Create histogram object
    TH1D* hist = new TH1D("hist", title, n_bins, min, max);

    // Fill histograms
    for (const double& m_inv : inv_masses) {
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

    // Create a canvas to draw on
    TCanvas* canvas = new TCanvas("canvas", "", 900, 700);

    // Draw histograms and info texts
    hist->SetStats(kFALSE);
    hist->SetXTitle("\\text{2 Particle Invariant Mass [GeV/c}^{2}\\text{]}");
    hist->SetYTitle(Form("\\text{Counts per %.2f [MeV/c}^{2}\\text{]}", bin_width*1000));
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

    // Save plot to TEX file
    const std::string file_name = results.decay_repr_str + std::string("_pair_inv_mass.tex");
    canvas->Print(file_name.c_str());
}
