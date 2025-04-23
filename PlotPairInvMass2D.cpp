// ROOT Includes
#include "TH2D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TColor.h"
#include "TFile.h"

// Local Includes
#include "starlyze.cpp" 
#include <iostream>

void PlotPairInvMass2D(const std::string& result_file_path = "slight.out") {
    // Read inn result
    const SimulationResult results = ReadSimulationResults(result_file_path);

    // Create ROOT output file before any plotting
    const std::string base_file_name = results.decay_repr_str 
                                     + std::string("_") + std::to_string(results.n_events)
                                     + std::string("_") + std::to_string(results.rnd_seed)
                                     + std::string("_pair_inv_mass_2d");
    const std::string root_file_name = base_file_name + std::string(".root");
    TFile* root_file = new TFile(root_file_name.c_str(), "recreate");

    // Create title for plot
    const char* title = Form("\\text{STARlight } | \\text{ Pb - Pb } \\sqrt{s_{NN}} = %.2f \\text{ TeV } | \\, %s", 
                             results.sqrt_s_NN/1000, results.decay_latex_str.c_str());

    // Seperate particle pairs invariant masses
    std::vector<double> m_inv_pairs_1;
    std::vector<double> m_inv_pairs_2;
    for (const Event& event : results.events) {
        m_inv_pairs_1.push_back(event.m_inv_pairs[0]);
        m_inv_pairs_2.push_back(event.m_inv_pairs[1]);
    }

    // Create pair inv mass vs pair inv mass histogram
    const double min_1 = std::min_element(m_inv_pairs_1.begin(), m_inv_pairs_1.end())[0];
    const double max_1 = std::max_element(m_inv_pairs_1.begin(), m_inv_pairs_1.end())[0];
    const double width_1 = FreedmanDiaconisBinWidth(m_inv_pairs_1);
    const int nbins_1 = (max_1 - min_1)/width_1;

    const double min_2 = std::min_element(m_inv_pairs_2.begin(), m_inv_pairs_2.end())[0];
    const double max_2 = std::max_element(m_inv_pairs_2.begin(), m_inv_pairs_2.end())[0];
    const double width_2 = FreedmanDiaconisBinWidth(m_inv_pairs_2);
    const int nbins_2 = (max_2 - min_2)/width_2;

    // Create histogram object
    TH2D* hist = new TH2D("hist", title, nbins_1, min_1, max_1, 
                                         nbins_2, min_2, max_2);

    // Fill histogram
    for (int i=0; i < results.n_events/2; i++) {
        hist->Fill(m_inv_pairs_1[i], m_inv_pairs_2[i]);
    }

    // Get peak
    const int bin_max = hist->GetMaximumBin();
    const double hist_peak = hist->GetYaxis()->GetBinCenter(bin_max);
    std::cout << hist_peak;

    // Text information about amount of events
    const char* events_info = Form("\\text{%i events}", results.n_events);
    TLatex* events_info_text = new TLatex(0.54, 0.80, events_info);
    events_info_text->SetNDC();

    // Create a canvas to draw on
    TCanvas* canvas = new TCanvas("canvas", "", 900, 700);

    // Change color pallete
    gStyle->SetPalette(kDeepSea);

    // Draw histograms and info texts
    hist->SetStats(kFALSE);
    hist->SetXTitle("\\text{1. pair invariant mass [GeV/c}^{2}\\text{]}");
    hist->SetYTitle("\\text{2. pair invariant mass [GeV/c}^{2}\\text{]}");
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
    hist->GetXaxis()->SetTitleOffset(1.0);
    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->GetXaxis()->SetLabelSize(0.035);
    hist->GetYaxis()->SetLabelSize(0.04);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->Draw("COLZ");
    events_info_text->Draw();

    // Save plot to TEX file
    const std::string tex_file_name = base_file_name + std::string(".tex");
    canvas->Print(tex_file_name.c_str());

    // Save canvas object to ROOT file
    canvas->Write();
}
