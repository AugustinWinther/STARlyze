// STD Includes
#include <cmath>     // std::sqrt, std::log, std::pow
#include <random>    // std::random_device, std::mt19937
#include <vector>    // std::vector
#include <fstream>   // std::ifstream
#include <sstream>   // std::istringstream
#include <algorithm> // std::sort, std::shuffle

// Constanta
static constexpr double kPROTON_MASS = 0.93827208816;

// Initialize RNG stuff
auto kRNG = std::default_random_engine {};

double FreedmanDiaconisBinWidth(std::vector<double> data) 
{
    std::sort(data.begin(), data.end());
    const double q1 = data[data.size() / 4];
    const double q3 = data[3 * data.size() / 4];
    const double iqr = q3 - q1;
    const double bin_width = 2 * iqr / std::pow(data.size(), 1.0/3.0);
    return bin_width;
}

std::vector<std::string> SplitStringBy(const std::string& string, 
                                       const char& delimiter) 
{
    std::vector<std::string> string_segments;
    std::istringstream string_stream(string);
    std::string segment;
    
    while (std::getline(string_stream, segment, delimiter)) {
        string_segments.push_back(segment);
    }

    return string_segments;
}

double ParticleIdToMass(const int& mcid) {
    // Returns the mass in GeV of the particle with the given
    // particle ID. Negative numbers are the anti-particle variant
    // and has the same mass, only different charge.
    // IDs follows the PDG mcid format.
    double mass;

    switch (mcid) {
        case   11:  // e-
        case  -11:  // e+
            mass = 0.0005109989499999999; 
            break;
        case   13:  // mu-
        case  -13:  // mu+
            mass = 0.1056583755;    
            break;    
        case  211:  // pi+
        case -211:  // pi-
            mass = 0.139570390983681;   
            break;
        case  321:  // K+
        case -321:  // K-
            mass = 0.49367659945804093; 
            break;
        case  2212:  // p
        case -2212:  // pbar
            mass = 0.93827208816; 
            break;
        default: 
            mass = 0; 
            break;
    }

    return mass;
}

std::string DecayIdToReprStr(const int& decay_id) {
    std::string repr_str;

    switch (decay_id) {
        case 443211:  // J/psi --> pi+pi-pi+pi-
            repr_str = std::string("jpsi_4pi"); 
            break;   
        case 443321211: // J/psi --> K+K-pi+pi-
            repr_str = std::string("jpsi_2K2pi"); 
            break;   
        default: 
            repr_str = std::string("NoReprStrFound"); 
            break;
    }

    return repr_str; 
}

std::string DecayIdToLatexStr(const int& decay_id) {
    std::string latex_str;

    switch (decay_id) {
        case  443011:  // J/psi --> e+e-
            latex_str = std::string("J/#psi #rightarrow e^{+}e^{-}"); 
            break;
        case  443013:  // J/Psi --> mu+mu-
            latex_str = std::string("J/#psi #rightarrow #mu^{+}#mu^{-}"); 
            break;   
        case 443211:  // J/psi --> pi+pi-pi+pi-
            latex_str = std::string("J/\\psi \\rightarrow \\pi^{+}\\pi^{-}\\pi^{+}\\pi^{-}"); 
            break;   
        case 443321211:  // J/psi --> pi+pi-K+K-
            latex_str = std::string("J/\\psi \\rightarrow K^{+}K^{-}\\pi^{+}\\pi^{-}"); 
            break;   
        case 4432212:  // J/psi --> proton anti-proton
            latex_str = std::string("J/#psi #rightarrow p#bar{p}"); 
            break;
        default: 
            latex_str = std::string("NO JETSET ID FOUND"); 
            break;
    }

    return latex_str;
}

class Track {
    public:
    double E, px, py, pz, pseudo_rap;

    Track(const double& px, const double& py, const double& pz, const double& m) {
        const double p_mag = std::sqrt(px*px + py*py + pz*pz);
        this->E = std::sqrt(p_mag*p_mag + m*m);
        this->px = px;
        this->py = py;
        this->pz = pz;
        this->pseudo_rap = 0.5 * std::log((p_mag + pz) / (p_mag - pz));
    }
};

class Event {
    public:
    double m_inv, p_trans;
    std::vector<double> m_inv_pair;

    Event(std::vector<Track> tracks) {
        // In real life, we don't know which particle is which in the detector.
        // Thus we shuffle the list of tracks to remove our knowldege of which
        // track is which particle.
        std::shuffle(tracks.begin(), tracks.end(), kRNG);

        // Calculate invariant mass of the system of the first pair of particles
        const double E1  =  tracks[0].E + tracks[1].E;
        const double px1 = tracks[0].px + tracks[1].px;
        const double py1 = tracks[0].py + tracks[1].py;
        const double pz1 = tracks[0].pz + tracks[1].pz;

        const double m_inv_1 = std::sqrt(E1*E1 - px1*px1 - py1*py1 - pz1*pz1);
        this->m_inv_pair.push_back(m_inv_1);

        // Calculate invariant mass of the system of the second pair of particles
        const double E2  =  tracks[2].E + tracks[3].E;
        const double px2 = tracks[2].px + tracks[3].px;
        const double py2 = tracks[2].py + tracks[3].py;
        const double pz2 = tracks[2].pz + tracks[3].pz;

        const double m_inv_2 = std::sqrt(E2*E2 - px2*px2 - py2*py2 - pz2*pz2);
        this->m_inv_pair.push_back(m_inv_2);

        // Calculate invariant mass of the system of all four particles
        const double  E =  E1 + E2;
        const double px = px1 + px2;
        const double py = py1 + py2;
        const double pz = pz1 + pz2;

        this->m_inv = std::sqrt(E*E - px*px - py*py - pz*pz);

        // Calculate transverse momentum of the system of all four particles
        this->p_trans = std::sqrt(px*px + py*py);
    }
};

class SimulationResult {
    public:
    int n_events;
    double sqrt_s_NN;
    std::string repr_str; 
    std::string decay_latex_str;
    std::vector<double> m_inv_list, p_trans_list;
    std::vector<std::vector<double>> m_inv_pair_list;

    SimulationResult(const std::vector<Event>& events, const int& decay_id,
                     const double& beam_1_gamma, const double& beam_2_gamma) {
        
        // Used to display in plots
        this->n_events = events.size();
        this->repr_str = DecayIdToReprStr(decay_id);
        this->decay_latex_str = DecayIdToLatexStr(decay_id);

        // Energy per nucleon (Only protons are accelerated)
        const double beam_1_E_N = kPROTON_MASS*beam_1_gamma;
        const double beam_2_E_N = kPROTON_MASS*beam_2_gamma;
        this->sqrt_s_NN = beam_1_E_N + beam_2_E_N;
        
        for (Event event : events) {
            this->m_inv_list.push_back(event.m_inv);
            this->m_inv_pair_list.push_back(event.m_inv_pair);
            this->p_trans_list.push_back(event.p_trans);
        }
    }
};

SimulationResult ReadSimulationResults(const std::string& result_file_path) {
    // Variables for track and event values
    double m, px, py, pz, beam_1_gamma, beam_2_gamma;
    int particle_id, decay_id;
    std::vector<Track> tracks;
    std::vector<Event> events;

    // Variables used for parsing
    int tracks_remaining_in_event;
    std::ifstream result_file(result_file_path);
    std::string line;
    std::vector<std::string> line_segments;

    // Go through each line in result_file
    while (std::getline(result_file, line)) {
        line_segments = SplitStringBy(line, ' ');

        if (line_segments[0] == std::string("CONFIG_OPT:")) {
            decay_id = std::stoi(line_segments[2]);
        } 
        else if (line_segments[0] == std::string("BEAM_1:")) {
            beam_1_gamma = std::stod(line_segments[3]);
        } 
        else if (line_segments[0] == std::string("BEAM_2:")) {
            beam_2_gamma = std::stod(line_segments[3]);
        } 
        else if (line_segments[0] == std::string("EVENT:")) {
            tracks_remaining_in_event = std::stoi(line_segments[2]);
        } 
        else if (line_segments[0] == std::string("TRACK:") && tracks_remaining_in_event != 0) {
            px = std::stod(line_segments[3]);
            py = std::stod(line_segments[4]);
            pz = std::stod(line_segments[5]);
            particle_id = std::stoi(line_segments[9]);
            m = ParticleIdToMass(particle_id);
            
            Track track(px, py, pz, m);
            tracks.push_back(track);

            tracks_remaining_in_event -= 1;
        }

        if (tracks_remaining_in_event == 0) {
            Event event(tracks);
            events.push_back(event);
            tracks.clear();
        }
    }

    result_file.close();
    return SimulationResult(events, decay_id, beam_1_gamma, beam_2_gamma);
}
