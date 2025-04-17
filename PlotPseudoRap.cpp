// ROOT Includes
#include "TH1D.h"
#include "TCanvas.h"
#include "TLatex.h"

// Local Includes
#include "starlyze.cpp" 

static constexpr double kPSEUDO_RAP_ACCEPT = 0.9;

void PlotPseudoRap(const std::string& result_file_path = "slight.out") {
    // Read inn result
    const SimulationResult results = ReadSimulationResults(result_file_path);

    // Create title for plot
    const char* title = Form("\\text{STARlight } | \\text{ Pb - Pb } \\sqrt{s_{NN}} = %.2f \\text{ TeV } | \\, %s; \\text{Event number}; \\text{Particles Detected}", 
                             results.sqrt_s_NN/1000, results.decay_latex_str.c_str());
    
    // Create bar chart values
    std::string bar_str[5] = {"0","1","2","3","4"};
    int bar_val[5] = {0, 0, 0, 0, 0};;
    int events_detected = 0;
    int particles_detected = 0;
    for (const Event& event : results.events) {
        // Count detected particles in event
        particles_detected = 0;
        for (const double& pseudo_rap : event.pseudo_raps) {
            if (-kPSEUDO_RAP_ACCEPT < pseudo_rap && pseudo_rap < kPSEUDO_RAP_ACCEPT) {
                particles_detected +=1;
            }
        }
        // Add particles detected to bar chart and check if event was detected
        if (particles_detected == 0) {bar_val[0] += 1;}
        if (particles_detected == 1) {bar_val[1] += 1;}
        if (particles_detected == 2) {bar_val[2] += 1;}
        if (particles_detected == 3) {bar_val[3] += 1;}
        if (particles_detected == 4) {bar_val[4] += 1; events_detected += 1;}
    }

    // Create histogram to become barchart
    TH1D* bar = new TH1D("bar",title,5,0,5);

    // Fill bar chart
    for (int i=1; i<6; i++) {
        bar->SetBinContent(i, bar_val[i-1]);
        bar->GetXaxis()->SetBinLabel(i, bar_str[i-1].c_str());
     }

    // Text information about amount of events
    const char* events_info = Form("\\text{%i events}", results.n_events);
    TLatex* events_info_text = new TLatex(0.54, 0.80, events_info);
    events_info_text->SetNDC();

    // Text information about events detected
    const char* detect_info = Form("\\text{where %i fully detected}", events_detected);
    TLatex* detect_info_text = new TLatex(0.54, 0.75, detect_info);
    detect_info_text->SetNDC();

    // Text information acceptence rate
    const char* accept_info = Form("\\text{Acceptence: } |\\eta| < %.1f", kPSEUDO_RAP_ACCEPT);
    TLatex* accept_info_text = new TLatex(0.54, 0.70, accept_info);
    accept_info_text->SetNDC();

    // Create a canvas to draw on
    TCanvas* canvas = new TCanvas("canvas", "", 900, 700);
    canvas->SetGrid();

    // Draw bar chart and info texts
    bar->SetStats(kFALSE);
    bar->SetXTitle("\\text{Number of particles detected}");
    bar->SetYTitle(Form("\\text{Number of events}"));
    bar->SetBarWidth(0.8);
    bar->SetBarOffset(0.1);
    bar->SetMinimum(0);
    bar->GetXaxis()->CenterTitle();
    bar->GetYaxis()->CenterTitle();
    bar->GetXaxis()->SetTitleOffset(1.0);
    bar->GetYaxis()->SetTitleOffset(1.2);
    bar->GetXaxis()->SetLabelSize(0.035);
    bar->GetYaxis()->SetLabelSize(0.04);
    bar->GetXaxis()->SetTitleSize(0.05);
    bar->GetYaxis()->SetTitleSize(0.05);
    bar->SetLineColor(kBlack);
    bar->SetFillColor(kP10Blue);
    bar->Draw("b");
    events_info_text->Draw();
    detect_info_text->Draw();
    accept_info_text->Draw();

    // Save plot to TEX file
    const std::string file_name = results.decay_repr_str + std::string("_pseudo_rap.tex");
    canvas->Print(file_name.c_str());
}
