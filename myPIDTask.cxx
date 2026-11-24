 #include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include <cmath>
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include <vector>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include <CCDB/BasicCCDBManager.h>
#include <fstream>
#include <stdexcept>
#include <string>
#include <TFile.h>
#include <memory>


using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace constants::math;

struct MyPIDTask
{
    Configurable<std::string> cfgReferenceRootFile{"cfgReferenceRootFile", "MeanPtResults.root", "ROOT file containing reference mean pT"};
    ConfigurableAxis cfgAxisCent{"cfgAxisCent", {0.0, 1.0, 5.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0}, "Centrality axis (percentile)"};
    Configurable<int> cfgCentralityChoice{"cfgCentralityChoice", 1, "Centrality estimator: 1=FT0C, 2=FT0M, 3=FDDM, 4=FV0A"};
    Configurable<bool> cfgRunDataFluc{"cfgRunDataFluc", false, "Run DATA fluctuations (C2, subevent)"};
    Configurable<std::string> cfgCCDBurl{"cfgCCDBurl", "https://alice-ccdb.cern.ch", "ccdb url"};
    Configurable<std::string> cfgCCDBUserPath{"cfgCCDBUserPath", "/Users/s/somadutt", "Base CCDB path"};

    static constexpr int KNEta = 9;
    static constexpr int KNCent = 8;
    static constexpr int KNPt = 8;

    inline static constexpr double ptEdges[KNPt + 1] = {0.0, 0.2, 0.3, 0.5, 0.7, 1.0, 2.0, 5.0, 10.0};
    // reference for <[pt]>
    std::unique_ptr<TProfile2D> referenceProfile;
    std::unique_ptr<TProfile> referenceMeanPtA;
    std::unique_ptr<TProfile> referenceMeanPtC;
    // reference for <n(pt)>
    std::unique_ptr<TProfile3D> referenceNPtProfile;
    std::unique_ptr<TProfile2D> referenceNPtA;
    std::unique_ptr<TProfile2D> referenceNPtC;
    
    std::unique_ptr<TProfile2D> referenceMeanPtASpecies;
    std::unique_ptr<TProfile2D> referenceMeanPtCSpecies;
    
    std::unique_ptr<TProfile3D> referenceMeanPtSpecies;
    std::unique_ptr<TProfile3D> referenceNPtASpecies;
    std::unique_ptr<TProfile3D> referenceNPtCSpecies;
    enum PIDIdx {
        kInclusiveIdx = 0,
        kPiMinusIdx,
        kPiPlusIdx,
        kPiAllIdx,
        kKaMinusIdx,
        kKaPlusIdx,
        kKaAllIdx,
        kAntiPrIdx,
        kPrIdx,
        kPrAllIdx,
        KNsp
    };
    enum ECentralityEstimator {
        kCentFT0C = 1,
        kCentFT0M = 2,
        kCentFDDM = 3,
        kCentFV0A = 4
    };
    
    template <typename TCollision>
    float getCentrality(const TCollision& col) const
    {
        if (cfgCentralityChoice.value == kCentFT0C)
            return col.centFT0C();

        if (cfgCentralityChoice.value == kCentFT0M)
            return col.centFT0M();

        if (cfgCentralityChoice.value == kCentFDDM)
            return col.centFDDM();

        if (cfgCentralityChoice.value == kCentFV0A)
            return col.centFV0A();

        return -1.f;
    }
    int getPtBin(double pt) const
    {
        for (int ipt = 0; ipt < KNPt; ++ipt) {

            if (pt >= ptEdges[ipt] &&
                pt < ptEdges[ipt + 1]) {

                return ipt;
            }
        }

        return -1;
    }
    
    // to convert percentages of centrality to an index centrality
    int getCentralityIndex(float cent) const
    {
        const double edges[9] = {
            0., 1., 5., 10., 20.,
            40., 60., 80., 100.
        };

        for (int i = 0; i < KNCent; ++i) {
            if (cent >= edges[i] && cent < edges[i + 1]) {
                return i;
            }
        }

        return -1;
    }

    inline static const std::vector<float> etaLw = {-0.8, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6};
    inline static const std::vector<float> etaUp = {0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8};
    //using MyCollisions = soa::Join<aod::Collisions, aod::CentFT0Cs, aod::PVMults, aod::MultsGlobal>;
    //using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullKa, aod::pidTOFbeta, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::TrackSelection, aod::pidTOFFullPi>;
    using GeneralCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::FT0sCorrected, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFDDMs, aod::CentFV0As, aod::CentNTPVs>;
    Filter collisionFilter = nabs(aod::collision::posZ) < 10.f;
    using AodCollisionsSel = soa::Filtered<GeneralCollisions>;
    using UnfilteredTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullKa, aod::pidTOFbeta, aod::pidTOFFullPi, aod::pidTOFFullPr, aod::pidTOFFullKa>;
    Filter trackFilter = aod::track::pt > 0.2f&&
                         aod::track::pt < 10.0f&&
                         requireGlobalTrackInFilter();
    using AodTracksSel = soa::Filtered<UnfilteredTracks>;
    template <typename TTrack>
    int getNsigmaPID(const TTrack& track)
    {
        // TPC PID selection
        bool isPion = std::abs(track.tpcNSigmaPi()) < 3.0f;
        bool isKaon = std::abs(track.tpcNSigmaKa()) < 3.0f;
        bool isProton = std::abs(track.tpcNSigmaPr()) < 3.0f;

        // TOF threshold
        float ptCutTOF = 0.3f;

        if (track.pt() > ptCutTOF) {

            // Above threshold, require a TOF match
            if (!track.hasTOF()) {
                return 0;
            }
            // Require both TPC and TOF PID
            isPion = isPion && std::abs(track.tofNSigmaPi()) < 3.0f;
            isKaon = isKaon && std::abs(track.tofNSigmaKa()) < 3.0f;
            isProton = isProton && std::abs(track.tofNSigmaPr()) < 3.0f;
        }

        // Reject ambiguous tracks
        if ((isPion && isKaon) || (isPion && isProton) || (isKaon && isProton)) {
            return 0;
        }
        // Return particle ID
        if (isPion) {
            return 1;
        } else if (isKaon) {
            return 2;
        } else if (isProton) {
            return 3;
        }
        return 0;
    }
  HistogramRegistry histos{"histos", {}};
  void init(InitContext const&)
  {
      //const AxisSpec centAxis{cfgAxisCent, "Centrality (%)"};
      const AxisSpec centAxis{std::vector<double>{0., 1., 5., 10., 20., 40., 60., 80., 100.}, "Centrality (%)"};
      const AxisSpec ptAxis{std::vector<double>{0.0, 0.2, 0.3, 0.5, 0.7,1.0, 2.0, 5.0, 10.0},"pT (GeV/c)" };
      // define axes you want to use
  //  const AxisSpec axisCounter{1, 0, +1, ""};
   // const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
   // const AxisSpec axisPt{nBinsPt, 0, 10, "p_{T}"};
    //const AxisSpec axisDeltaPt{100, -1.0, +1.0, "#Delta(p_{T})"};
    histos.add("tpcsignal", "tpcsignal", kTH1F, {{200,0.,10.}});
    histos.add("etaall", "etaall", kTH1F, {{100,-2.,2.}});
    histos.add("eta", "eta", kTH1F, {{100,-2.,2.}});
    histos.add("tracktime", "tracktime", kTH1F, {{200, -50., 50.}});
    // after cuts
    histos.add("tpcSignalVsP", "TPC signal vs momentum;p (GeV/c);TPC dE/dx", kTH2F, {{200, 0., 10.}, {200, 0., 1000.}});
    histos.add("betaVsP","TOF #beta vs momentum;p (GeV/c);#beta", kTH2F, {{200, 0., 5.}, {200, 0., 1.2}});
    // tpc
    histos.add("nSigmaVsPitpc", "nsigmaPitpc;p (GeV/c);nsigmaPi", kTH2F, {{200, 0., 5.}, {200, -10, 10.}});
    histos.add("nSigmaVsPrtpc", "nsigmaPrtpc;p (GeV/c);nsigmaPr", kTH2F, {{200, 0., 5.}, {200, -10, 10.}});
    histos.add("nSigmaVsKatpc", "nsigmaKatpc;p (GeV/c);nsigmaKa", kTH2F, {{200, 0., 5.}, {200, -10, 10.}});
    // tof
    histos.add("nSigmaVsKatof", "nsigmaKatof;p (GeV/c);nsigmaKa", kTH2F, {{200, 0., 5.}, {200, -10, 10.}});
    histos.add("nSigmaVsPitof", "nsigmaPitof;p (GeV/c);nsigmaPi", kTH2F, {{200, 0., 5.}, {200, -10, 10.}});
    histos.add("nSigmaVsPrtof", "nsigmaPrtof;p (GeV/c);nsigmaPr", kTH2F, {{200, 0., 5.}, {200, -10, 10.}});
      
    histos.add("nSigmatofvstpcPi", "nSigmatofvstpcPi; nSigmatpc;nSigmatof", kTH2F, {{500, -10., 10.}, {200, -10, 10.}});
    histos.add("nSigmatofvstpcKa", "nSigmatofvstpcKa; nSigmatpc;nSigmatof", kTH2F, {{500, -10., 10.}, {200, -10, 10.}});
    histos.add("nSigmatofvstpcPr", "nSigmatofvstpcPr; nSigmatpc;nSigmatof", kTH2F, {{500, -10., 10.}, {200, -10, 10.}});
      
    // histograms pt distribution of individual particles
    // pion
    histos.add("ptPi", "Pions;p_{T} (GeV/c); entries", kTH1F, {{100, 0., 10.}});
    histos.add("etaPi", "Pions;#eta;Counts", kTH1F, {{100, -1., 1.}});
    histos.add("phiPi", "Pions;#varphi (rad);entries", kTH1F, {{100, 0., 6.2832}});
    // kaon
    histos.add("ptKa", "Kaons;p_{T} (GeV/c);entries", kTH1F, {{100, 0., 10.}});
    histos.add("etaKa", "Kaons;#eta;Counts", kTH1F, {{100, -1., 1.}});
    histos.add("phiKa", "Kaons;#varphi (rad);entries", kTH1F, {{100, 0., 6.2832}});
    // proton
    histos.add("ptPr", "Protons;p_{T} (GeV/c);entries", kTH1F, {{100, 0., 10.}});
    histos.add("etaPr", "Protons;#eta;entries", kTH1F, {{100, -1., 1.}});
    histos.add("phiPr", "Protons;#varphi (rad);entries", kTH1F, {{100, 0., 6.2832}});
    // centrality for particles
    histos.add("centPi", "Pions;Centrality (%);entries", kTH1F, {{90, 0., 90.}});
    histos.add("centKa", "Kaons;Centrality (%);entries", kTH1F, {{90, 0., 90.}});
    histos.add("centPr", "Protons;Centrality (%);entries", kTH1F, {{90, 0., 90.}});
      
    // histogram of N tracks global on x axis, and Ntracks PV on y axis
    histos.add("tracksize", "tracksize;tracksize;entries", kTH1F, {{4000, 0., 4000.}});
    histos.add("multNtrackPV", "Ntrackpv;Ntrackpv;entries", kTH1F, {{4000, 0., 4000.}});
    histos.add("tracksize_PVtracks","PV tracks vs tracksize; NtrackPV;tracksize", kTH2F, {{4000, 0., 4000.}, {4000, 0., 4000.}});
    // make a plot of tracksize (Ntrackglobal) vs centrality
      histos.add("tracksize_centrality", "N_{tracks} vs centrality;Centrality (%);N_{tracks}", kTH2F, {{101,-0.5,100.5}, {4000, 0., 4000.}});
    // plot mean pt in 1 centrality and 1 eta bin
    histos.add("MeanpT", "Event mean p_{T};[p_{T}] (GeV/c);Entries", kTH1F, {{100, 0., 2.}});
    histos.add("NtracksPV_MeanpT", "Mean p_{T} vs NtracksPV;N_{tracks,PV};[p_{T}] (GeV/c)", kTH2F, {{400, 0., 4000.}, {100, 0., 2.}});
    //histos.add("Prof_Ntrackpv_MeanpT", "Average mean p_{T} vs NtrackPV;N_{tracks,PV};#LT[p_{T}]#GT (GeV/c)", kTProfile, {{400, 0., 4000.}});
    histos.add("Prof_Ntrackpv_MeanpT", "Mean p_{T} vs NtracksPV and particle species;N_{tracks,PV};Particle species;Mean p_{T} (GeV/c)", kTProfile2D, {{1000, 0., 4000.}, {KNsp, -0.5, KNsp - 0.5}});
    // calculate ensemble mean 𝑝t in centrality and 𝜂 bins
      // all particles
      histos.add("Prof_MeanpT_Cent_Eta", "<[pT]> in centrality and in eta bins;Centrality (%);Eta bin;Mean pT (GeV/c)", kTProfile2D, {{101,-0.5,100.5}, {KNEta, -0.5, KNEta - 0.5}});
      // mean pt in eta ranges
      histos.add("Prof_MeanpT_Cent_EtaA", "Mean pT subevent A;Centrality (%);<[pT]_{A}>", kTProfile,
          {{101, -0.5, 100.5}});
      histos.add("Prof_MeanpT_Cent_EtaC", "Mean pT subevent C;Centrality (%);<[pT]_{C}>",
          kTProfile, {{101, -0.5, 100.5}});
      // mean pt per event in centrality and eta bins
      histos.add("hMeanPtEvent_Cent_Eta_inclusive", "Event-wise mean pT inclusive;Centrality (%);Eta bin;Mean pT (GeV/c)", kTH3F, {centAxis, {KNEta, -0.5, KNEta - 0.5}, {100, 0.0, 5.0}});
      histos.add("Prof_MeanpT_Cent_Eta_Species", "Mean pT;Centrality (%);Eta bin;Particle species", kTProfile3D,{{101, -0.5, 100.5}, {KNEta, -0.5, KNEta - 0.5}, {KNsp, -0.5, KNsp - 0.5}});
      
      histos.add("Prof_MeanpT_Cent_EtaA_Species", "Mean pT subevent A;Centrality (%);Particle species",
          kTProfile2D, {{101, -0.5, 100.5}, {KNsp, -0.5, KNsp - 0.5}});
      histos.add("Prof_MeanpT_Cent_EtaC_Species", "Mean pT subevent C;Centrality (%);Particle species",
          kTProfile2D, {{101, -0.5, 100.5}, {KNsp, -0.5, KNsp - 0.5}});
      
      // now we make everything for mean pt versus centrality for different eta bins
      // first a histogram of mean pt versus centrality for the first eta bin
      histos.add("hMeanPtevent_cent_eta0", "mean pT [pt] vs centrality ;centrality (%); [pT] (GeV/c)", kTH2F, {centAxis, {200, 0.0, 2.0}});
      histos.add("hMeanPtevent_cent_eta1", "mean pT [pt] vs centrality ;centrality (%); [pT] (GeV/c)", kTH2F, {centAxis, {200, 0.0, 2.0}});
      histos.add("hMeanPtevent_cent_eta2", "mean pT [pt] vs centrality ;centrality (%); [pT] (GeV/c)", kTH2F, {centAxis, {200, 0.0, 2.0}});
      histos.add("hMeanPtevent_cent_eta3", "mean pT [pt] vs centrality ;centrality (%); [pT] (GeV/c)", kTH2F, {centAxis, {200, 0.0, 2.0}});
      histos.add("hMeanPtevent_cent_eta4", "mean pT [pt] vs centrality ;centrality (%); [pT] (GeV/c)", kTH2F, {centAxis, {200, 0.0, 2.0}});
      // plots for c2
      histos.add("hPbar11", "pbar11 distribution;pbar11 (GeV/c);Entries", kTH1F, {{200, -0.5, 0.5}});
      histos.add("pbar22", "pbar22;pbar22;Entries", kTH1F, {{200, -0.5, 0.5}});
      histos.add("p12", "p12;p12;Entries", kTH1F, {{200, 0., 2.}});
      histos.add("p11", "p11;p11;Entries", kTH1F, {{200, 0., 2.}});
      histos.add("c2", "c2;c2;Entries", kTH1F, {{200, -0.5, 0.5}});
      histos.add("Prof_c2_Cent_Eta", "c2 in centrality and in eta bins;Centrality (%);Eta bin;c2 ", kTProfile2D, {{101,-0.5,100.5}, {KNEta, -0.5, KNEta - 0.5}});
      // c2(centrality, eta, particle)
      histos.add("Prof_c2_Cent_Eta_Species", "c2 for different species;Centrality (%);Eta bin;Particle species", kTProfile3D,{{101, -0.5, 100.5}, {KNEta, -0.5, KNEta - 0.5}, {KNsp, -0.5, KNsp - 0.5}});
      
      histos.add("c2vseta", "c2 (10-20% centrality);Eta bin; c2", kTH2F, {{KNEta, -0.5, KNEta - 0.5}, {200, -1.0, 1.0}});
      // plot delta pt versus centrality
      histos.add("hdelta_ptVScent", "delta pt vs centrality ;centrality (%); delta_pT (GeV/c)", kTH2F, {{101,-0.5,100.5}, {200, -1.0, 1.0}});
      histos.add("Prof_deltaptVScent", "Mean #Delta p_{T} vs centrality;Centrality (%);<#Delta p_{T}> (GeV/c)", kTProfile, {{101,-0.5,100.5}});
      
      // for c2sub
      histos.add("c2subVScent", "c2sub vs centrality ;centrality (%); c2sub", kTH2F, {{101,-0.5,100.5}, {200, -1.0, 1.0}});
      histos.add("Prof_c2sub_Cent_EtaPair", "Subevent c2;Centrality (%);Eta pair index", kTProfile2D, {{101, -0.5, 100.5}, {4, 0.5, 4.5}});
      histos.add("Prof_v0Denominator_Species", "C2sub per species;Centrality (%);Particle species",
          kTProfile2D, {{101, -0.5, 100.5}, {KNsp, -0.5, KNsp - 0.5}});
      
      // for n(pt)
      // nPtVsCent
      // first plot N(pt) 1D histogram
      histos.add("hNpt","pT distribution;pT (GeV/c);N tracks", kTH1F, {ptAxis});
      histos.add("hNpt2", "N(pT);pT bin index;Number of particles", kTH1F, {{KNPt, -0.5, KNPt - 0.5}});
      histos.add("nPtvsPt", "n(Pt) versus pt ; Pt; n(Pt)", kTProfile, {ptAxis});
      histos.add("Prof_NPtVsPt", "N(pT) vs pT;p_{T} (GeV/c);N(p_{T})", kTProfile, {ptAxis});
      histos.add("Prof_NPtVsPtpions", "N(pT) vs pT;p_{T} (GeV/c);N(p_{T})", kTProfile, {ptAxis});
      
      histos.add("Prof_NPt_Cent_Pt_Species", "Mean particle counts;Centrality (%);pT (GeV/c);Particle species", kTProfile3D, {{101, -0.5, 100.5}, ptAxis, {KNsp, -0.5, KNsp - 0.5}});
      
      histos.add("nPtVsCent", "n(Pt) versus cent (0.3-0.5 GeV/c) ; centrality; n(Pt)", kTH2F, {{101,-0.5,100.5}, {200, 0, 10.0}});
      histos.add("Prof_nPtVsCent","Mean n(pT) vs centrality (0.3-0.5 GeV/c);Centrality (%);<n(pT)>",
          kTProfile, {{101, -0.5, 100.5}});
      histos.add("Prof_nPt_Cent_Eta_Pt", "<n(pT)>;Centrality (%);Eta bin;pT (GeV/c)",kTProfile3D, {{101,-0.5,100.5}, {KNEta, -0.5, KNEta - 0.5}, ptAxis});
      histos.add("Prof_nPt_Cent_PtA", "<n_{A}(pT)>;Centrality (%);pT (GeV/c)", kTProfile2D, {{101, -0.5, 100.5}, ptAxis});
      histos.add("Prof_nPt_Cent_PtC", "<n_{C}(pT)>;Centrality (%);pT (GeV/c)", kTProfile2D, {{101, -0.5, 100.5}, ptAxis});
      
      histos.add("Prof_nPt_Cent_PtA_species", "Mean nA(pT) per species;Centrality (%);Particle species;pT (GeV/c)", kTProfile3D, {{101, -0.5, 100.5}, {KNsp, -0.5, KNsp - 0.5}, ptAxis});
      histos.add("Prof_nPt_Cent_PtC_species", "Mean nC(pT) per species;Centrality (%);Particle species;pT (GeV/c)", kTProfile3D, {{101, -0.5, 100.5}, {KNsp, -0.5, KNsp - 0.5}, ptAxis});
      
      histos.add("deltaNPtVsPt"," #delta n(p_{T});p_{T} (GeV/c);#delta n(p_{T})", kTH2F, {ptAxis, {200, -5.0, 5.0}});
      histos.add("deltaNPtVsDeltaPt", "#delta n(p_{T}) vs #delta[p_{T}];#delta[p_{T}] (GeV/c);#delta n(p_{T})", kTH2F, {{200, -0.5, 0.5}, {500, -5.0, 5.0}});
      histos.add("deltaNPtVsDeltaPt2", "#delta n(p_{T}) vs #delta[p_{T}];#delta[p_{T}] (GeV/c);#delta n(p_{T})", kTH2F, {{200, -0.5, 0.5}, {500, -5.0, 5.0}});
      histos.add("deltaNPtVsDeltaPt3", "#delta n(p_{T}) vs #delta[p_{T}];#delta[p_{T}] (GeV/c);#delta n(p_{T})", kTH2F, {{200, -0.5, 0.5}, {500, -5.0, 5.0}});
      histos.add("Prof_DeltaPtDeltaNPt", "Correlation of pT fluctuations;p_{T} (GeV/c);<#delta[p_{T}] #delta n(p_{T})>", kTProfile, {ptAxis});
      
      histos.add("deltaNPtA_vs_deltaPtC", "#delta n_{A}(p_{T}) vs #delta[p_{T}]_{C}; #delta[p_{T}]_{C} (GeV/c);#delta n_{A}(p_{T})", kTH2F,{{100, -0.2, 0.2}, {100, -0.2, 0.2}});
      histos.add("deltaNPtC_vs_deltaPtA", "#delta n_{C}(p_{T}) vs #delta[p_{T}]_{A}; #delta[p_{T}]_{A} (GeV/c);#delta n_{C}(p_{T})", kTH2F,{{100, -0.2, 0.2}, {100, -0.2, 0.2}});
      histos.add("v0vspt", " v0 vs pt; pt (GeV/c);v0", kTH2F,{{200, 0, 5}, {200, -0.5, 0.5}});
      
      // Event-wise particle counts versus pT
      histos.add("ptInclusive", "Inclusive counts;p_{T} (GeV/c);counts", kTH1F, {ptAxis});
      histos.add("ptPion", "Pion counts;p_{T} (GeV/c);counts", kTH1F, {ptAxis});
      histos.add("ptKaon", "Kaon counts;p_{T} (GeV/c);counts", kTH1F, {ptAxis});
      histos.add("ptProton", "Proton counts;p_{T} (GeV/c);counts", kTH1F, {ptAxis});
      
      // for v0 calculation
      
      // Numerator: <delta nA(pT) delta[pT]C>
      histos.add("Prof_v0Numerator", "v0 numerator;Centrality (%);pT (GeV/c);Covariance", kTProfile2D, {{101, -0.5, 100.5}, ptAxis});
      //  <nA(pT)>
      histos.add("Prof_v0MeanNPt", "Mean nA(pT);Centrality (%);pT (GeV/c);<nA(pT)>", kTProfile2D, {{101, -0.5, 100.5}, ptAxis});
      // C2sub/<delta[pT]A delta[pT]C>
      histos.add("Prof_v0Denominator","C2sub;Centrality (%);C2sub (GeV/c)^{2}",kTProfile, {{101, -0.5, 100.5}});
      
      
      histos.add("Prof_v0Numerator_Species", "v0 numerator per species;Centrality (%);Particle species;pT (GeV/c)", kTProfile3D, {{101, -0.5, 100.5}, {KNsp, -0.5, KNsp - 0.5}, ptAxis});
      histos.add("Prof_v0MeanNPt_Species", "Mean nA(pT) per species;Centrality (%);Particle species;pT (GeV/c)", kTProfile3D, {{101, -0.5, 100.5}, {KNsp, -0.5, KNsp - 0.5}, ptAxis});
      
      if (doprocessFluc) {

          // Open the ROOT file from the mean pT calculation
          std::unique_ptr<TFile> referenceFile(
              TFile::Open(cfgReferenceRootFile.value.c_str(), "READ")
          );

          if (!referenceFile || referenceFile->IsZombie()) {
              throw std::runtime_error(
                  "Cannot open MeanPtResults.root"
              );
          }
          auto* meanSpecies = referenceFile->Get<TProfile3D>(
              "my-p-i-d-task/Prof_MeanpT_Cent_Eta_Species"
          );

          if (!meanSpecies) {
              throw std::runtime_error(
                  "Cannot find Prof_MeanpT_Cent_Eta_Species"
              );
          }

          // Check histogram dimensions
          if (meanSpecies->GetNbinsX() != 101 ||
              meanSpecies->GetNbinsY() != KNEta ||
              meanSpecies->GetNbinsZ() != KNsp) {

              throw std::runtime_error(
                  "Incorrect species mean pT reference binning"
              );
          }

          // Clone the profile
          referenceMeanPtSpecies.reset(
              static_cast<TProfile3D*>(
                  meanSpecies->Clone("referenceMeanPtSpecies")
              )
          );

          // Detach it from the ROOT file
          referenceMeanPtSpecies->SetDirectory(nullptr);
          // Retrieve the previously calculated reference profile
          auto* meanProfile = referenceFile->Get<TProfile2D>(
              "my-p-i-d-task/Prof_MeanpT_Cent_Eta"
          );

          if (!meanProfile) {
              throw std::runtime_error(
                  "Cannot find Prof_MeanpT_Cent_Eta"
              );
          }

          // Check that the dimensions are correct
          if (meanProfile->GetNbinsX() != 101 ||
              meanProfile->GetNbinsY() != KNEta) {
              throw std::runtime_error(
                  "Incorrect reference profile binning"
              );
          }

          // Make an independent copy of the histogram
          referenceProfile.reset(
              static_cast<TProfile2D*>(
                  meanProfile->Clone("referenceMeanPt")
              )
          );
          // Detach from the input ROOT file
          referenceProfile->SetDirectory(nullptr);
          
          auto* nPtProfile = referenceFile->Get<TProfile3D>(
                  "my-p-i-d-task/Prof_nPt_Cent_Eta_Pt"
              );

              if (!nPtProfile) {
                  throw std::runtime_error(
                      "Cannot find Prof_nPt_Cent_Eta_Pt"
                  );
              }

              // Check all three dimensions
              if (nPtProfile->GetNbinsX() != 101 ||
                  nPtProfile->GetNbinsY() != KNEta ||
                  nPtProfile->GetNbinsZ() != KNPt) {

                  throw std::runtime_error(
                      "Incorrect n(pT) reference binning"
                  );
              }

              referenceNPtProfile.reset(
                  static_cast<TProfile3D*>(
                      nPtProfile->Clone("referenceMeanNPt")
                  )
              );

              referenceNPtProfile->SetDirectory(nullptr);
          
          auto* meanPtA = referenceFile->Get<TProfile>(
                      "my-p-i-d-task/Prof_MeanpT_Cent_EtaA"
                  );

                  auto* meanPtC = referenceFile->Get<TProfile>(
                      "my-p-i-d-task/Prof_MeanpT_Cent_EtaC"
                  );

                  if (!meanPtA || !meanPtC) {
                      throw std::runtime_error(
                          "Cannot find subevent mean pT references"
                      );
                  }

                  if (meanPtA->GetNbinsX() != 101 ||
                      meanPtC->GetNbinsX() != 101) {
                      throw std::runtime_error(
                          "Incorrect subevent mean pT reference binning"
                      );
                  }

                  referenceMeanPtA.reset(
                      static_cast<TProfile*>(
                          meanPtA->Clone("referenceMeanPtA")
                      )
                  );

                  referenceMeanPtC.reset(
                      static_cast<TProfile*>(
                          meanPtC->Clone("referenceMeanPtC")
                      )
                  );

                  referenceMeanPtA->SetDirectory(nullptr);
                  referenceMeanPtC->SetDirectory(nullptr);

                  auto* nPtA = referenceFile->Get<TProfile2D>(
                      "my-p-i-d-task/Prof_nPt_Cent_PtA"
                  );

                  auto* nPtC = referenceFile->Get<TProfile2D>(
                      "my-p-i-d-task/Prof_nPt_Cent_PtC"
                  );

                  if (!nPtA || !nPtC) {
                      throw std::runtime_error(
                          "Cannot find subevent n(pT) references"
                      );
                  }

                  if (nPtA->GetNbinsX() != 101 ||
                      nPtA->GetNbinsY() != KNPt ||
                      nPtC->GetNbinsX() != 101 ||
                      nPtC->GetNbinsY() != KNPt) {
                      throw std::runtime_error(
                          "Incorrect subevent n(pT) reference binning"
                      );
                  }

                  referenceNPtA.reset(
                      static_cast<TProfile2D*>(
                          nPtA->Clone("referenceNPtA")
                      )
                  );

                  referenceNPtC.reset(
                      static_cast<TProfile2D*>(
                          nPtC->Clone("referenceNPtC")
                      )
                  );

                  referenceNPtA->SetDirectory(nullptr);
                  referenceNPtC->SetDirectory(nullptr);
          
          auto* meanPtASpecies = referenceFile->Get<TProfile2D>(
              "my-p-i-d-task/Prof_MeanpT_Cent_EtaA_Species"
          );

          auto* meanPtCSpecies = referenceFile->Get<TProfile2D>(
              "my-p-i-d-task/Prof_MeanpT_Cent_EtaC_Species"
          );

          if (!meanPtASpecies || !meanPtCSpecies) {
              throw std::runtime_error(
                  "Cannot find species mean pT references"
              );
          }

          referenceMeanPtASpecies.reset(
              static_cast<TProfile2D*>(
                  meanPtASpecies->Clone("referenceMeanPtASpecies")
              )
          );

          referenceMeanPtCSpecies.reset(
              static_cast<TProfile2D*>(
                  meanPtCSpecies->Clone("referenceMeanPtCSpecies")
              )
          );

          referenceMeanPtASpecies->SetDirectory(nullptr);
          referenceMeanPtCSpecies->SetDirectory(nullptr);
          auto* nPtASpecies = referenceFile->Get<TProfile3D>(
              "my-p-i-d-task/Prof_nPt_Cent_PtA_species"
          );

          auto* nPtCSpecies = referenceFile->Get<TProfile3D>(
              "my-p-i-d-task/Prof_nPt_Cent_PtC_species"
          );

          if (!nPtASpecies || !nPtCSpecies) {
              throw std::runtime_error(
                  "Missing species n(pT) reference profiles"
              );
          }

          if (nPtASpecies->GetNbinsX() != 101 ||
              nPtASpecies->GetNbinsY() != KNsp ||
              nPtASpecies->GetNbinsZ() != KNPt ||
              nPtCSpecies->GetNbinsX() != 101 ||
              nPtCSpecies->GetNbinsY() != KNsp ||
              nPtCSpecies->GetNbinsZ() != KNPt) {
              throw std::runtime_error(
                  "Incorrect species n(pT) reference binning"
              );
          }

          referenceNPtASpecies.reset(
              static_cast<TProfile3D*>(
                  nPtASpecies->Clone("referenceNPtASpecies")
              )
          );

          referenceNPtCSpecies.reset(
              static_cast<TProfile3D*>(
                  nPtCSpecies->Clone("referenceNPtCSpecies")
              )
          );

          referenceNPtASpecies->SetDirectory(nullptr);
          referenceNPtCSpecies->SetDirectory(nullptr);

      }
  }
  void processMean(AodCollisionsSel::iterator const& coll, AodTracksSel const& tracks) {
    //double sumWi[KNEta]{}, sumWipti[KNEta]{};
    float cent = getCentrality(coll);
    //double sumWi[KNEta]{}, sumWipti[KNEta]{};
    double sumWi[KNsp][KNEta]{}, sumWipti[KNsp][KNEta]{};
    double sumWi2[KNsp][KNEta]{}, sumWi2pti2[KNsp][KNEta]{};
    double sumWi2pti1[KNsp][KNEta]{};
    // Counts in each pT bin for this collision
    double N[KNsp][KNEta][KNPt]{};
    histos.fill(HIST("tracksize"), tracks.size());
    histos.fill(HIST("multNtrackPV"), coll.multNTracksPV());
    histos.fill(HIST("tracksize_PVtracks"), coll.multNTracksPV(), tracks.size());
      // plot of Ntrackglobal versus centrality
    histos.fill(HIST("tracksize_centrality"), cent, tracks.size());
    for (const auto& track : tracks) {
        float p = track.p();
        double pt = track.pt();
        auto sign = track.sign();
        float eta = track.eta();
        float beta = track.beta();
        double w = 1.0;
        int id = getNsigmaPID(track);
        bool isPi = (id == 1);
        bool isKa = (id == 2);
        bool isPr = (id == 3);
        bool isSpecies[KNsp] = {
            true,
            isPi && sign < 0, isPi && sign > 0, isPi,
            isKa && sign < 0, isKa && sign > 0, isKa,
            isPr && sign < 0, isPr && sign > 0, isPr
        };
        if (cent >= 10.0 && cent < 20.0) {
            histos.fill(HIST("ptInclusive"), pt);
            if (isPi)
                histos.fill(HIST("ptPion"), pt);
            if (isKa)
                histos.fill(HIST("ptKaon"), pt);
            if (isPr)
                histos.fill(HIST("ptProton"), pt);
        }
        
        // calculate the mean pt per event.
        for (int isp = 0; isp < KNsp; ++isp) {
            if (!isSpecies[isp])
                continue;
            for (int ieta = 0; ieta < KNEta; ++ieta) {
                if (eta <= etaLw[ieta] || eta > etaUp[ieta])
                    continue;
                // for p11 (meanpt)
                double pt = track.pt();
                sumWi[isp][ieta] += w;
                sumWipti[isp][ieta] += w * pt;
                // for p22
                sumWi2[isp][ieta] += w * w;
                sumWi2pti2[isp][ieta] += (w * w) * (pt * pt);
                // for p12
                sumWi2pti1[isp][ieta] += (w * w) * pt;
                // count number of particles in pt bins
                int ipt = getPtBin(pt);
                if (ipt >= 0 && ipt < KNPt) {
                    N[isp][ieta][ipt] += 1.0;
                }
            }
        }
      //tpc
        double nSigmaPi = track.tpcNSigmaPi();
        double nSigmaPr = track.tpcNSigmaPr();
        double nSigmaKa = track.tpcNSigmaKa();
      // tof
        double nSigmaKatof = track.tofNSigmaKa();
        double nSigmaPitof = track.tofNSigmaPi();
        double nSigmaPrtof = track.tofNSigmaPr();
        
        float tpcSignal = track.tpcSignal();
        float trackTime = track.trackTime();
    // fill after cuts
        histos.fill(HIST("tracktime"), trackTime);
        histos.fill(HIST("eta"), eta);
        histos.fill(HIST("tpcsignal"), tpcSignal);
        histos.fill(HIST("tpcSignalVsP"), p, tpcSignal);
    // tpc
        histos.fill(HIST("nSigmaVsPitpc"), p, nSigmaPi);
        histos.fill(HIST("nSigmaVsPrtpc"), p, nSigmaPr);
        histos.fill(HIST("nSigmaVsKatpc"), p, nSigmaKa);
    // tof
        if (track.hasTOF()) {

            histos.fill(HIST("betaVsP"), p, track.beta());

            histos.fill(HIST("nSigmaVsKatof"), p, track.tofNSigmaKa());
            histos.fill(HIST("nSigmaVsPitof"), p, track.tofNSigmaPi());
            histos.fill(HIST("nSigmaVsPrtof"), p, track.tofNSigmaPr());

            histos.fill(HIST("nSigmatofvstpcPi"), track.tofNSigmaPi(), track.tpcNSigmaPi());
            histos.fill(HIST("nSigmatofvstpcKa"), track.tofNSigmaKa(), track.tpcNSigmaKa());
            histos.fill(HIST("nSigmatofvstpcPr"), track.tofNSigmaPr(), track.tpcNSigmaPr());
        }
    // Fill histograms for identified particles
        int pid = getNsigmaPID(track);
        if (pid == 1) {
            histos.fill(HIST("ptPi"), track.pt());
            histos.fill(HIST("etaPi"), track.eta());
            histos.fill(HIST("phiPi"), track.phi());
            histos.fill(HIST("centPi"), cent);
        } else if (pid == 2) {
            histos.fill(HIST("ptKa"), track.pt());
            histos.fill(HIST("etaKa"), track.eta());
            histos.fill(HIST("phiKa"), track.phi());
            histos.fill(HIST("centKa"), cent);
        } else if (pid == 3) {
            histos.fill(HIST("ptPr"), track.pt());
            histos.fill(HIST("etaPr"), track.eta());
            histos.fill(HIST("phiPr"), track.phi());
            histos.fill(HIST("centPr"), cent);
        }
    }
      for (int isp = 0; isp < KNsp; ++isp) {
          if (sumWi[isp][0] < 1.0)
              continue;
          // mean pt over ensemble
          histos.fill(HIST("Prof_Ntrackpv_MeanpT"), coll.multNTracksPV(), isp, sumWipti[isp][0] / sumWi[isp][0]);
      }
      // loop over all the eta intervals
      // mean pt per ensemble, that is what the Tprofile does
      for (int ieta = 0; ieta < KNEta; ++ieta) {
          // Total number of particles in each eta interval
          double nTracks = sumWi[kInclusiveIdx][ieta];
          // Skip eta intervals without tracks
          if (nTracks <= 0.0)
              continue;
          // loop over momentum bins
          for (int ipt = 0; ipt < KNPt; ++ipt) {
              double ptCenter = 0.5 * (ptEdges[ipt] + ptEdges[ipt + 1]);
              double ptWidth = ptEdges[ipt + 1] - ptEdges[ipt];
              if (ieta == 0) {
                  for (int isp = 0; isp < KNsp; ++isp) {
                      double totalSpecies = sumWi[isp][ieta];
                      if (totalSpecies <= 0.0)
                          continue;
                      histos.fill(HIST("Prof_NPt_Cent_Pt_Species"), cent, ptCenter, isp, N[isp][ieta][ipt]/(totalSpecies * ptWidth));
                  }
              }
              double nPt = N[kInclusiveIdx][ieta][ipt] / nTracks;
              histos.fill(HIST("Prof_nPt_Cent_Eta_Pt"), cent, ieta, ptCenter, nPt);
              
              if (ieta == 0 && ipt == 2) {
                  histos.fill(HIST("nPtVsCent"), cent, nPt);
                  histos.fill(HIST("Prof_nPtVsCent"), cent, nPt);
              }
              if (cent >= 10.0f && cent < 20.0f && ieta==0) {
                  histos.fill(HIST("nPtvsPt"), ptCenter, nPt);
                  histos.fill(HIST("Prof_NPtVsPt"), ptCenter, N[kInclusiveIdx][ieta][ipt]);
                  histos.fill(HIST("Prof_NPtVsPtpions"), ptCenter, N[kPiAllIdx][ieta][ipt]);
                  histos.fill(HIST("hNpt"),ptCenter, N[kInclusiveIdx][ieta][ipt]);
                  histos.fill(HIST("hNpt2"),ipt, N[kInclusiveIdx][ieta][ipt]);
              }
          }
          // Skip eta intervals without selected tracks
          // this is for all particles
          if (sumWi[kInclusiveIdx][ieta] >= 1.0){
              // Mean pT of this collision in this eta interval meanPt =  p11
              double meanPt = sumWipti[kInclusiveIdx][ieta] / sumWi[kInclusiveIdx][ieta];
              // ensemble mean pt <[pt]> stored in tprofile
              histos.fill(HIST("Prof_MeanpT_Cent_Eta"), cent, ieta, meanPt);
              histos.fill(HIST("hMeanPtEvent_Cent_Eta_inclusive"), cent, ieta, meanPt);
          }
          for (int isp = 0; isp < KNsp; ++isp) {
                  if (sumWi[isp][ieta] <= 0.0)
                      continue;
                  double meanPt = sumWipti[isp][ieta] / sumWi[isp][ieta];
                  histos.fill(HIST("Prof_MeanpT_Cent_Eta_Species"), cent, ieta, isp, meanPt);
              }
      }
      
      // calculate mean n(pt) values for the two different eta ranges for all particle species
      for (int isp = 0; isp < KNsp; ++isp) {
          double nTracksA = 0.0;
          double nTracksC = 0.0;
          // loop over eta A, total number of particles in those eta ranges
          for (int ieta = 1; ieta <= 4; ++ieta) {
              nTracksA += sumWi[isp][ieta];
          }
          // loop over eta C
          for (int ieta = 5; ieta <= 8; ++ieta) {
              nTracksC += sumWi[isp][ieta];
          }
          // calculate in each pt bin the number of particle
          for (int ipt = 0; ipt < KNPt; ++ipt) {
                  double ptCenter = 0.5 * (ptEdges[ipt] + ptEdges[ipt + 1]);
                  double NA = 0.0;
                  double NC = 0.0;
                  // Number of particles in this pT bin, and eta range A
                  for (int ieta = 1; ieta <= 4; ++ieta) {
                      NA += N[isp][ieta][ipt];
                  }
                  // Number of particles in this pT bin
                  for (int ieta = 5; ieta <= 8; ++ieta) {
                      NC += N[isp][ieta][ipt];
                  }
              // save <n(pt)>A and <n(pt)>C
                  if (nTracksA > 0.0) {
                      double nPtA = NA / nTracksA;
                      histos.fill(HIST("Prof_nPt_Cent_PtA_species"), cent, isp, ptCenter, nPtA);
                      if (isp == kInclusiveIdx) {
                          histos.fill(HIST("Prof_nPt_Cent_PtA"), cent, ptCenter, nPtA);
                      }
                  }
                  if (nTracksC > 0.0) {
                      double nPtC = NC / nTracksC;
                      histos.fill(HIST("Prof_nPt_Cent_PtC_species"), cent, isp, ptCenter, nPtC);
                      if (isp == kInclusiveIdx) {
                          histos.fill(HIST("Prof_nPt_Cent_PtC"), cent, ptCenter, nPtC);
                      }
                  }
              }
      }
      
      //calculate mean(pt) values for the two different eta ranges
      for (int isp = 0; isp < KNsp; ++isp) {
          double nTracksA = 0.0;
          double nTracksC = 0.0;

          double sumPtA = 0.0;
          double sumPtC = 0.0;

          // Number of particles and sum of pT in eta range A
          for (int ieta = 1; ieta <= 4; ++ieta) {
              nTracksA += sumWi[isp][ieta];
              sumPtA += sumWipti[isp][ieta];
          }
          // Number of particles and sum of pT in eta range C
          for (int ieta = 5; ieta <= 8; ++ieta) {
              nTracksC += sumWi[isp][ieta];
              sumPtC += sumWipti[isp][ieta];
          }
          if (nTracksA > 0.0) {
              double meanPtA = sumPtA / nTracksA;
              histos.fill(HIST("Prof_MeanpT_Cent_EtaA_Species"), cent, isp, meanPtA);
              if (isp == kInclusiveIdx) {
                  histos.fill(HIST("Prof_MeanpT_Cent_EtaA"), cent, meanPtA);
              }
          }
          if (nTracksC > 0.0) {
              double meanPtC = sumPtC / nTracksC;
              histos.fill(HIST("Prof_MeanpT_Cent_EtaC_Species"), cent, isp, meanPtC);
              if (isp == kInclusiveIdx) {
                  histos.fill(HIST("Prof_MeanpT_Cent_EtaC"), cent, meanPtC);
              }
          }
      }
      
      // make sure there is atleast one track from all particles (inclusiveidx)
      // event wise mean pt for first eta bin
      // for each collision, i fill ( NtracksPV, [pt])
      // in one NtrackPV bin, you have a collection of {[pt1],[pt2],[pt3],}
      if (sumWi[kInclusiveIdx][0] >= 1.0) {
          // mean pt versus centrality for first eta bin only eta 0
          histos.fill(HIST("hMeanPtevent_cent_eta0"), cent, sumWipti[kInclusiveIdx][0] / sumWi[kInclusiveIdx][0]);
      }
      if (sumWi[kInclusiveIdx][1] >= 1.0) {
          // mean pt versus centrality for first eta bin only eta 1
          histos.fill(HIST("hMeanPtevent_cent_eta1"), cent, sumWipti[kInclusiveIdx][1] / sumWi[kInclusiveIdx][1]);
      }
      if (sumWi[kInclusiveIdx][2] >= 1.0) {
          // mean pt versus centrality for first eta bin only eta 2
          histos.fill(HIST("hMeanPtevent_cent_eta2"), cent, sumWipti[kInclusiveIdx][2] / sumWi[kInclusiveIdx][2]);
      }
      if (sumWi[kInclusiveIdx][3] >= 1.0) {
          // mean pt versus centrality for first eta bin only eta 3
          histos.fill(HIST("hMeanPtevent_cent_eta3"), cent, sumWipti[kInclusiveIdx][3] / sumWi[kInclusiveIdx][3]);
      }
      if (sumWi[kInclusiveIdx][4] >= 1.0) {
          // mean pt versus centrality for first eta bin only eta 4
          histos.fill(HIST("hMeanPtevent_cent_eta4"), cent, sumWipti[kInclusiveIdx][4] / sumWi[kInclusiveIdx][4]);
      }
      // plot mean pt in first eta bin and in 1 centrality range.
      if (cent >= 10.0f && cent < 20.0f) {
          histos.fill(HIST("MeanpT"), sumWipti[kInclusiveIdx][0] / sumWi[kInclusiveIdx][0]);
      }
      histos.fill(HIST("NtracksPV_MeanpT"), coll.multNTracksPV(), sumWipti[kInclusiveIdx][0] / sumWi[kInclusiveIdx][0]);
      // mean pt for 1 eta bin, for the first eta bin, this is the distribution of [pt] but then we have 8 eta bins and then the it looks like mean pt eta for 1 centrality for example 10-20 (hMeanPtEvent_Eta_Cent10to20)
      // histos.fill(HIST("MeanpT"), sumWipti[kInclusiveIdx][0] / sumWi[kInclusiveIdx][0]);
  }
    PROCESS_SWITCH(MyPIDTask, processMean, "Calculate reference mean pT", true);
    
    void processFluc(AodCollisionsSel::iterator const& coll, AodTracksSel const& tracks){
        // Find the centrality bin for this collision
        float cent = getCentrality(coll);
        //int icent = getCentralityIndex(cent);
        double sumWi[KNsp][KNEta]{}, sumWipti[KNsp][KNEta]{};
        double sumWi2[KNsp][KNEta]{}, sumWi2pti2[KNsp][KNEta]{};
        double sumWi2pti1[KNsp][KNEta]{};
        // Counts in each pT bin for this collision
        double N[KNsp][KNEta][KNPt]{};
        // store the pbar11 values for each eta interval
        double pbar11Eta[KNEta]{};
        for (const auto& track : tracks) {
            auto sign = track.sign();
            float eta = track.eta();
            double w = 1.0;
            int id = getNsigmaPID(track);
            bool isPi = (id == 1);
            bool isKa = (id == 2);
            bool isPr = (id == 3);
            bool isSpecies[KNsp] = {
                true,
                isPi && sign < 0, isPi && sign > 0, isPi,
                isKa && sign < 0, isKa && sign > 0, isKa,
                isPr && sign < 0, isPr && sign > 0, isPr
            };
            // calculate the mean pt per event. p11 and
            for (int isp = 0; isp < KNsp; ++isp) {
                if (!isSpecies[isp])
                    continue;
                for (int ieta = 0; ieta < KNEta; ++ieta) {
                    if (eta <= etaLw[ieta] || eta > etaUp[ieta])
                        continue;
                    // for p11 (meanpt)
                    double pt = track.pt();
                    sumWi[isp][ieta] += w;
                    sumWipti[isp][ieta] += w * pt;
                    // for p22
                    sumWi2[isp][ieta] += w * w;
                    sumWi2pti2[isp][ieta] += (w * w) * (pt * pt);
                    // for p12
                    sumWi2pti1[isp][ieta] += (w * w) * pt;
                    // count number of particles in pt bins
                    int ipt = getPtBin(pt);
                    if (ipt >= 0 && ipt < KNPt) {
                        N[isp][ieta][ipt] += 1.0;
                    }
                }
            }
        }
        if (cent >= 0.0f && cent < 100.0f) {
            int centBin = referenceProfile->GetXaxis()->FindBin(cent);
            int centBinNPt = referenceNPtProfile->GetXaxis()->FindBin(cent);
            // Loop over eta intervals
            for (int ieta = 0; ieta < KNEta; ++ieta) {
                // Only calculate if this interval has tracks
                if (sumWi[kInclusiveIdx][ieta] < 2.0)
                    continue;
                int bin = referenceProfile->GetBin(centBin, ieta + 1);
                // Check that the reference bin contains entries
                if (referenceProfile->GetBinEntries(bin) <= 0)
                    continue;
                double meanReference = referenceProfile->GetBinContent(bin);
                // mean pT = p11
                double p11 = sumWipti[kInclusiveIdx][ieta] / sumWi[kInclusiveIdx][ieta];
                double p22 = sumWi2pti2[kInclusiveIdx][ieta] / sumWi2[kInclusiveIdx][ieta];
                double p12 = sumWi2pti1[kInclusiveIdx][ieta] / sumWi2[kInclusiveIdx][ieta];
                double tau1 = sumWi2[kInclusiveIdx][ieta] / (sumWi[kInclusiveIdx][ieta] * sumWi[kInclusiveIdx][ieta]);
                // Reference mean pt  [<pt>] for this centrality and eta interval
                // Calculate pbar11 (which is delta_pt)
                double pbar11 = p11 - meanReference;
                pbar11Eta[ieta] = pbar11;
                // for delta_n(pt)
                double totalN = sumWi[kInclusiveIdx][ieta];
                if (totalN <= 0.0)
                    continue;
                // Calculate delta n(pT) for this eta interval
                for (int ipt = 0; ipt < KNPt; ++ipt) {
                    double ptWidth = ptEdges[ipt + 1] - ptEdges[ipt];
                    double ptCenter = 0.5 * (ptEdges[ipt] + ptEdges[ipt + 1]);
                    double nPt = N[kInclusiveIdx][ieta][ipt] /totalN;
                    int binNPt = referenceNPtProfile->GetBin(centBinNPt, ieta + 1, ipt + 1);
                    if (referenceNPtProfile->GetBinEntries(binNPt) <= 0)
                        continue;
                    double meanNPtReference = referenceNPtProfile->GetBinContent(binNPt);
                    double delta_n_Pt = nPt - meanNPtReference;
                    if (cent >= 40.0f && cent < 60.0f && ieta == 0) {
                        histos.fill(HIST("deltaNPtVsPt"), ptCenter, delta_n_Pt);
                        if (ipt == 2) {
                            histos.fill(HIST("deltaNPtVsDeltaPt"), pbar11, delta_n_Pt);
                        }
                        if (ipt == 4){
                            histos.fill(HIST("deltaNPtVsDeltaPt2"), pbar11, delta_n_Pt);
                        }
                        if (ipt == 5){
                            histos.fill(HIST("deltaNPtVsDeltaPt3"), pbar11, delta_n_Pt);
                        }
                        histos.fill(HIST("Prof_DeltaPtDeltaNPt"), ptCenter, pbar11 * delta_n_Pt);
                    }
                    
                }
                // plot delta pt versus centrality for 1 eta bin
                if (ieta==0){
                    histos.fill(HIST("hdelta_ptVScent"), cent, pbar11);
                    histos.fill(HIST("Prof_deltaptVScent"), cent, pbar11);
                }

                // calculate pbar22 and c2
                double pbar22 = p22 - (2.0f * p12 * meanReference) + (meanReference * meanReference);
                double c2nom = (pbar11 * pbar11) - (tau1 * pbar22);
                double c2denom = (1.0f - tau1);
                if (c2denom <= 0.0)
                    continue;
                double c2 = c2nom / c2denom;
                if (!std::isfinite(c2))
                    continue;
                // fill the tprofile thing
                histos.fill(HIST("Prof_c2_Cent_Eta"), cent, ieta, c2);
                // For now, inspect one specific combination: 10-20% centrality, eta bin 0
                if (cent >= 10.0f && cent < 20.0f) {
                    histos.fill(HIST("c2vseta"), ieta, c2);
                    if (ieta ==0){
                        histos.fill(HIST("hPbar11"), pbar11);
                        histos.fill(HIST("pbar22"), pbar22);
                        histos.fill(HIST("c2"), c2);
                    }
                }
            }
            // c2 for particle species
            if (cent >= 0.0f && cent < 100.0f) {
                int centBinSpecies = referenceMeanPtSpecies->GetXaxis()->FindBin(cent);
                for (int ieta = 0; ieta < KNEta; ++ieta) {
                    for (int isp = 0; isp < KNsp; ++isp) {
                        // Need at least two particles of this species
                        if (sumWi[isp][ieta] < 2.0)
                            continue;
                        // Find reference for centrality, eta and species
                        int bin = referenceMeanPtSpecies->GetBin(centBinSpecies, ieta + 1, isp + 1);
                        if (referenceMeanPtSpecies->GetBinEntries(bin) <= 0)
                            continue;
                        double meanReference = referenceMeanPtSpecies->GetBinContent(bin);
                        double p11 = sumWipti[isp][ieta] / sumWi[isp][ieta];
                        double p22 = sumWi2pti2[isp][ieta] / sumWi2[isp][ieta];
                        double p12 = sumWi2pti1[isp][ieta] / sumWi2[isp][ieta];
                        double tau1 = sumWi2[isp][ieta] / (sumWi[isp][ieta] * sumWi[isp][ieta]);
                        double pbar11 = p11 - meanReference;
                        double pbar22 = p22 - 2.0 * p12 * meanReference + meanReference * meanReference;
                        // Two-particle correlator
                        double numerator = pbar11 * pbar11 - tau1 * pbar22;
                        double denominator = 1.0 - tau1;
                        if (denominator <= 0.0)
                            continue;
                        double c2 = numerator / denominator;
                        if (!std::isfinite(c2))
                            continue;
                        histos.fill(HIST("Prof_c2_Cent_Eta_Species"), cent, ieta, isp, c2);
                    }
                }
            }
            // Calculate C2sub for the different particle species
            for (int isp = 0; isp < KNsp; ++isp) {
                double nTracksA = 0.0;
                double nTracksC = 0.0;
                double sumPtA = 0.0;
                double sumPtC = 0.0;

                // Eta range A (-0.8 to 0)
                for (int ieta = 1; ieta <= 4; ++ieta) {
                    nTracksA += sumWi[isp][ieta];
                    sumPtA += sumWipti[isp][ieta];
                }

                // Eta range C (0 to 0.8)
                for (int ieta = 5; ieta <= 8; ++ieta) {
                    nTracksC += sumWi[isp][ieta];
                    sumPtC += sumWipti[isp][ieta];
                }

                if (nTracksA > 0.0 && nTracksC > 0.0) {

                    // Event-wise mean pT
                    double meanPtA = sumPtA / nTracksA;
                    double meanPtC = sumPtC / nTracksC;

                    // Find centrality bins of the reference profiles
                    int centBinA = referenceMeanPtASpecies->GetXaxis()->FindBin(cent);
                    int centBinC = referenceMeanPtCSpecies->GetXaxis()->FindBin(cent);
                    // Find the species bins
                    int binA = referenceMeanPtASpecies->GetBin(centBinA, isp + 1);
                    int binC = referenceMeanPtCSpecies->GetBin(centBinC, isp + 1);

                    // Check that both reference bins contain entries
                    if (referenceMeanPtASpecies->GetBinEntries(binA) <= 0 ||
                        referenceMeanPtCSpecies->GetBinEntries(binC) <= 0) {
                        continue;
                    }

                    // Mean pT references for this species
                    double meanReferenceA = referenceMeanPtASpecies->GetBinContent(binA);
                    double meanReferenceC = referenceMeanPtCSpecies->GetBinContent(binC);
                    // Mean pT fluctuations
                    double deltaPtA = meanPtA - meanReferenceA;
                    double deltaPtC = meanPtC - meanReferenceC;

                    // C2sub for all particle species
                    histos.fill(HIST("Prof_v0Denominator_Species"),cent, isp, deltaPtA * deltaPtC);
                    if (isp == kInclusiveIdx) {
                        histos.fill(HIST("Prof_v0Denominator"),cent, deltaPtA * deltaPtC);
                    }
                    int centBinNPtA = referenceNPtASpecies->GetXaxis()->FindBin(cent);
                    int centBinNPtC = referenceNPtCSpecies->GetXaxis()->FindBin(cent);
                    for (int ipt = 0; ipt < KNPt; ++ipt) {
                        double ptCenter = 0.5 * (ptEdges[ipt] + ptEdges[ipt + 1]);
                        double NA = 0.0;
                        double NC = 0.0;
                        // Tracks in the current pT bin in eta range A
                        for (int ieta = 1; ieta <= 4; ++ieta) {
                            NA += N[isp][ieta][ipt];
                        }
                        // Tracks in the current pT bin in eta range C
                        for (int ieta = 5; ieta <= 8; ++ieta) {
                            NC += N[isp][ieta][ipt];
                        }
                        // Event-wise normalized n(pt)
                        double nPtA = NA / nTracksA;
                        double nPtC = NC / nTracksC;

                        // Reference bins: centrality, species, pT
                        int binNPtA = referenceNPtASpecies->GetBin(centBinNPtA, isp + 1, ipt + 1);
                        int binNPtC = referenceNPtCSpecies->GetBin(centBinNPtC, isp + 1, ipt + 1);
                        // Check that both reference bins have entries
                        if (referenceNPtASpecies->GetBinEntries(binNPtA) <= 0 ||
                            referenceNPtCSpecies->GetBinEntries(binNPtC) <= 0) {
                            continue;
                        }
                        // Ensemble mean n(pT) references
                        double meannPtA = referenceNPtASpecies->GetBinContent(binNPtA);
                        double meannPtC = referenceNPtCSpecies->GetBinContent(binNPtC);
                        // Event-wise n(pT) fluctuations
                        double deltaNPtA = nPtA - meannPtA;
                        double deltaNPtC = nPtC - meannPtC;
                        // Numerator for each species
                        histos.fill(HIST("Prof_v0Numerator_Species"), cent, isp, ptCenter, deltaNPtA * deltaPtC);
                        histos.fill(HIST("Prof_v0MeanNPt_Species"), cent, isp, ptCenter, nPtA);
                        if (isp == kInclusiveIdx) {
                            histos.fill(HIST("Prof_v0Numerator"), cent, ptCenter, deltaNPtA * deltaPtC);
                            histos.fill(HIST("Prof_v0MeanNPt"), cent, ptCenter, nPtA);

                            if (cent >= 10.0f && cent < 20.0f && ipt == 4) {
                                histos.fill(HIST("deltaNPtA_vs_deltaPtC"),deltaPtC, deltaNPtA);
                                histos.fill(HIST("deltaNPtC_vs_deltaPtA"),deltaPtA, deltaNPtC);
                            }
                        }
                    }
                }
            }

            // now we calculate c2sub in two different eta ranges
            for (int ietaA = 1; ietaA <= 4; ++ietaA) {
                int ietaC = KNEta - ietaA;
                double c2sub = pbar11Eta[ietaA] * pbar11Eta[ietaC];
                histos.fill(HIST("c2subVScent"), cent, c2sub);
                histos.fill(HIST("Prof_c2sub_Cent_EtaPair"), cent, ietaA, c2sub);
            }
            if (cent >= 10.0f && cent < 20.0f) {
                // to see if they look the same as meanpt
                histos.fill(HIST("p12"), sumWi2pti1[kInclusiveIdx][0] / sumWi2[kInclusiveIdx][0]);
                histos.fill(HIST("p11"), sumWipti[kInclusiveIdx][0] / sumWi[kInclusiveIdx][0]);
            }
        }
    }
    PROCESS_SWITCH(MyPIDTask, processFluc, "Calculate pT fluctuations", true);
 
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MyPIDTask>(cfgc)};
}
