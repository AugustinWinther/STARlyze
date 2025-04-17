// STD Includes
#include <cmath>     // std::sqrt, std::log, std::pow
#include <random>    // std::random_device, std::mt19937
#include <vector>    // std::vector
#include <fstream>   // std::ifstream
#include <sstream>   // std::istringstream
#include <algorithm> // std::sort, std::shuffle

// Constants 
static constexpr double kELECTRON_MASS = 0.000510998928;
static constexpr double kPROTON_MASS = 0.938272046;
static constexpr double kMUON_MASS = 0.1056583755;
static constexpr double kPION_MASS = 0.13957018;
static constexpr double kKAON_MASS = 0.493677;

// Particle IDs follow the PDG format
enum ParticleId { 
    kELECTRON_ID = 11,
    kPROTON_ID = 2212,
    kMUON_ID = 13,
    kPION_ID = 211,
    kKAON_ID = 321
};

// DecayID (PROD_PID in STARlight) follows STARlight numbering scheme
enum DecayId { 
    kJPSI_2K2PI = 443321211,
    kJPSI_4PI = 443211,
    kJPSI_2MU = 443013,
    kJPSI_2E = 443011,
    kJPSI_2P = 4432212,
};

// Initialize RNG stuff
std::default_random_engine kRNG = std::default_random_engine {};

// Returns optimal bin-width for data. 
// Used for plotting histograms in ROOT Macros
double FreedmanDiaconisBinWidth(std::vector<double> data) {
    std::sort(data.begin(), data.end());
    const double q1 = data[data.size() / 4];
    const double q3 = data[3 * data.size() / 4];
    const double iqr = q3 - q1;
    const double bin_width = 2 * iqr / std::pow(data.size(), 1.0/3.0);
    return bin_width;
}

// Returns vector containg sub-strings seperated by passed delimiter
std::vector<std::string> SplitStringBy(const std::string& string, 
                                       const char& delimiter) {
    std::vector<std::string> string_segments;
    std::istringstream string_stream(string);
    std::string segment;
    
    while (std::getline(string_stream, segment, delimiter)) {
        string_segments.push_back(segment);
    }

    return string_segments;
}

// Returns the mass in GeV of the particle with the given
// particle ID. Negative numbers are the anti-particle variant.
double ParticleIdToMass(const int& particle_id) {
    double mass;

    switch (particle_id) {
        case  kELECTRON_ID:  // e-
        case -kELECTRON_ID:  // e+
            mass = kELECTRON_MASS; 
            break;
        case  kPROTON_ID:  // p
        case -kPROTON_ID:  // pbar
            mass = kPROTON_MASS; 
            break;
        case  kMUON_ID:  // mu-
        case -kMUON_ID:  // mu+
            mass = kMUON_MASS;    
            break;    
        case  kPION_ID:  // pi+
        case -kPION_ID:  // pi-
            mass = kPION_MASS;  
            break;
        case  kKAON_ID:  // K+
        case -kKAON_ID:  // K-
            mass = kKAON_MASS; 
            break;
        default: 
            mass = 0; 
            break;
    }

    return mass;
}

// Returns string reprsentation of the decay
std::string DecayIdToReprStr(const int& decay_id) {
    std::string repr_str;

    switch (decay_id) {
        case kJPSI_2K2PI: // J/psi -> K+ K- pi+ pi-
            repr_str = std::string("jpsi_2K2pi"); 
            break;   
        case kJPSI_4PI:   // J/psi -> pi+ pi- pi+ pi-
            repr_str = std::string("jpsi_4pi"); 
            break;   
        case  kJPSI_2MU:  // J/Psi -> mu+ mu-
            repr_str = std::string("jpsi_2mu"); 
            break;     
        case  kJPSI_2E:   // J/psi -> e+ e-
            repr_str = std::string("jpsi_2e"); 
            break;
        case kJPSI_2P:   // J/psi -> p pbar
            repr_str = std::string("jpsi_2p"); 
            break;
        default: 
            repr_str = std::string("NoReprStrFound"); 
            break;
    }

    return repr_str; 
}

// Returns LaTeX reprsentation of the decay
std::string DecayIdToLatexStr(const int& decay_id) {
    std::string latex_str;

    switch (decay_id) {
        case kJPSI_2K2PI: // J/psi -> K+ K- pi+ pi-
            latex_str = std::string("J/\\psi \\rightarrow K^{+}K^{-}\\pi^{+}\\pi^{-}"); 
            break; 
        case kJPSI_4PI:   // J/psi -> pi+ pi- pi+ pi-
            latex_str = std::string("J/\\psi \\rightarrow \\pi^{+}\\pi^{-}\\pi^{+}\\pi^{-}"); 
            break;   
        case  kJPSI_2MU:  // J/Psi -> mu+ mu-
            latex_str = std::string("J/\\psi \\rightarrow \\mu^{+}\\mu^{-}"); 
            break;     
        case  kJPSI_2E:   // J/psi -> e+ e-
            latex_str = std::string("J/\\psi \\rightarrow e^{+}e^{-}"); 
            break;
        case kJPSI_2P:   // J/psi -> p pbar
            latex_str = std::string("J/\\psi \\rightarrow p\\overline{p}"); 
            break;
        default: 
            latex_str = std::string("NO DECAY ID FOUND"); 
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
    std::vector<double> m_inv_pairs, pseudo_raps;

    Event(std::vector<Track> tracks) {
        // In real life, we don't know which particle is which in the detector.
        // Thus we shuffle the list of tracks to remove our knowldege of which
        // track is which particle.
        std::shuffle(tracks.begin(), tracks.end(), kRNG);

        double E1, E2, px1, px2, py1, py2, pz1, pz2, m_inv_1, m_inv_2;

        // Calculate invariant mass of the system of the first pair of particles
        E1  =  tracks[0].E + tracks[1].E;
        px1 = tracks[0].px + tracks[1].px;
        py1 = tracks[0].py + tracks[1].py;
        pz1 = tracks[0].pz + tracks[1].pz;

        m_inv_1 = std::sqrt(E1*E1 - px1*px1 - py1*py1 - pz1*pz1);
        this->m_inv_pairs.push_back(m_inv_1);

        // Calculate invariant mass of the system of the second pair of particles
        // if there is a second pair
        if (tracks.size() == 4) {
            E2  =  tracks[2].E + tracks[3].E;
            px2 = tracks[2].px + tracks[3].px;
            py2 = tracks[2].py + tracks[3].py;
            pz2 = tracks[2].pz + tracks[3].pz;

            m_inv_2 = std::sqrt(E2*E2 - px2*px2 - py2*py2 - pz2*pz2);
            this->m_inv_pairs.push_back(m_inv_2);
        } else {
            E2  = 0;
            px2 = 0;
            py2 = 0;
            pz2 = 0; 
        }

        // Calculate invariant mass of the system of all particles
        const double  E =  E1 + E2;
        const double px = px1 + px2;
        const double py = py1 + py2;
        const double pz = pz1 + pz2;

        this->m_inv = std::sqrt(E*E - px*px - py*py - pz*pz);

        // Calculate transverse momentum of the system of all four particles
        this->p_trans = std::sqrt(px*px + py*py);

        // Add all pseudo rapidities to pseudo rap. list
        for (const Track& track : tracks){
            this->pseudo_raps.push_back(track.pseudo_rap);
        }
    }
};

class SimulationResult {
    public:
    int n_events;
    double sqrt_s_NN;
    std::string decay_repr_str; 
    std::string decay_latex_str;
    std::vector<Event> events;

    SimulationResult(const std::vector<Event>& events, const int& decay_id,
                     const double& beam_1_gamma, const double& beam_2_gamma) {
        
        // Used to display in plots
        this->n_events = events.size();
        this->decay_repr_str = DecayIdToReprStr(decay_id);
        this->decay_latex_str = DecayIdToLatexStr(decay_id);

        // Energy per nucleon (Only protons are accelerated)
        const double beam_1_E_N = kPROTON_MASS*beam_1_gamma;
        const double beam_2_E_N = kPROTON_MASS*beam_2_gamma;
        this->sqrt_s_NN = beam_1_E_N + beam_2_E_N;
        
        this->events = events;
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
