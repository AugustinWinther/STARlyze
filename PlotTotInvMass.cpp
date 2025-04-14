// ROOT Includes
#include "TH1D.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLatex.h"
#include "TText.h"

// Local Includes
#include "starlyze.cpp" 

void PlotTotInvMass(const std::string& result_file_path = "slight.out") {
    // Read inn result
    const SimulationResult results = ReadSimulationResults(result_file_path);

    // Create title for plot
    char* title = Form("STARlight | Pb-Pb #sqrt{s_{NN}} = %.2f TeV | %s ", 
                       results.sqrt_s_NN/1000, 
                       results.decay_latex_str.c_str());

    // Calculate histogram properties
    const double min = std::min_element(results.m_inv_list.begin(), 
                                              results.m_inv_list.end())[0];
    const double max = std::max_element(results.m_inv_list.begin(), 
                                              results.m_inv_list.end())[0];
    const double bin_width = FreedmanDiaconisBinWidth(results.m_inv_list);
    const int n_bins = (max - min)/bin_width;

    // Create histogram object
    TH1D *hist = new TH1D("hist", title, n_bins, min, max);

    // Fill histograms
    for (const double& m_inv : results.m_inv_list) {
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

    // Text information about the invariant mass
    const char *m_info = Form("#splitline{Peak @ %.4f GeV/c^{2}}{FWHM = %.3f keV/c^{2}}", 
                              hist_peak, fwhm*1000000);
    TLatex *m_info_text = new TLatex(0.13, 0.70, m_info);
    m_info_text->SetNDC();
    m_info_text->SetTextAlign(0);
    m_info_text->SetTextSize(0.04);

    // Text information about the amount of events
    const char *e_info = Form("%i events", results.n_events);
    TText *e_info_text = new TText(0.13, 0.65, e_info);
    e_info_text->SetNDC();
    e_info_text->SetTextSize(0.04);

    // Create a canvas to draw on
    TCanvas *canvas = new TCanvas("canvas", "", 900, 700);

    // Draw histograms and info texts
    hist->SetStats(kFALSE);
    hist->SetXTitle("4 Particle Invariant Mass [GeV/c^{2}]");
    hist->SetYTitle(Form("Counts per %.2f [keV/c^{2}]", bin_width*1000000));
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
    hist->GetXaxis()->SetTitleOffset(0.8);
    hist->GetYaxis()->SetTitleOffset(1.0);
    hist->GetXaxis()->SetLabelSize(0.035);
    hist->GetYaxis()->SetLabelSize(0.04);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->SetLineColor(kBlack);
    hist->SetFillColor(kP10Blue);
    hist->Draw();
    m_info_text->Draw();
    e_info_text->Draw();

    // Save plot to SVG file
    canvas->Print("PlotTotInvMass.svg");
}
