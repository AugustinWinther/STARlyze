// ROOT Includes
#include "TH2D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TColor.h"

// Local Includes
#include "starlyze.cpp" 

void PlotPairInvMass2D(const std::string& result_file_path = "slight.out") {
    // Read inn result
    const SimulationResult results = ReadSimulationResults(result_file_path);

    // Create title for plot
    const char* title = Form("\\text{STARlight } | \\text{ Pb - Pb } \\sqrt{s_{NN}} = %.2f \\text{ TeV } | \\, %s", 
                             results.sqrt_s_NN/1000, results.decay_latex_str.c_str());

    // Seperate particle pairs invariant mases
    std::vector<double> m_inv_pair_1;
    std::vector<double> m_inv_pair_2;
    for (const std::vector<double>& m_inv_pair : results.m_inv_pair_list) {
        m_inv_pair_1.push_back(m_inv_pair[0]);
        m_inv_pair_2.push_back(m_inv_pair[1]);
    }

    // Create pair inv mass vs pair inv mass histogram
    const double min_1 = std::min_element(m_inv_pair_1.begin(), m_inv_pair_1.end())[0];
    const double max_1 = std::max_element(m_inv_pair_1.begin(), m_inv_pair_1.end())[0];
    const double width_1 = FreedmanDiaconisBinWidth(m_inv_pair_1);
    const int nbins_1 = (max_1 - min_1)/width_1;

    const double min_2 = std::min_element(m_inv_pair_2.begin(), m_inv_pair_2.end())[0];
    const double max_2 = std::max_element(m_inv_pair_2.begin(), m_inv_pair_2.end())[0];
    const double width_2 = FreedmanDiaconisBinWidth(m_inv_pair_2);
    const int nbins_2 = (max_2 - min_2)/width_2;

    TH2D *hist = new TH2D("hist", title, nbins_1, min_1, max_1, 
                                         nbins_2, min_2, max_2);

    // Fill histogram
    for (int i=0; i < results.n_events/2; i++) {
        hist->Fill(m_inv_pair_1[i], m_inv_pair_2[i]);
    }

    // Text information about amount of events
    const char *events_info = Form("\\text{%i events}", results.n_events);
    TLatex *events_info_text = new TLatex(0.55, 0.65, events_info);
    events_info_text->SetNDC();
    // events_info_text->SetTextSize(0.04); Leads sto seg fault

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
    const std::string file_name = results.repr_str + std::string("_pair_inv_mass_2d.tex");
    canvas->Print(file_name.c_str());
}