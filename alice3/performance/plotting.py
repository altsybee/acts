import acts
import acts.examples
import acts.examples.root

#from AliceActsPythonBindings import ResPlotToolConfig

#Residual Plots

#resPlotToolConfig = ResPlotToolConfig()
resPlotToolConfig = acts.examples.root.ResPlotToolConfig()
#binning = resPlotToolConfig.varBinning
#binning["Eta"] = acts.Axis.regular(80, -4,    4,     "#eta")
#binning["Pt"] =  acts.Axis.regular(1000, -0,    50,     "pT [GeV/c]")
#binning["Residual_phi"] =  acts.Axis.regular(200, -0.02,    0.02,     "r_{#phi} [rad]")
#binning["Residual_theta"] =  acts.Axis.regular(200, -0.02,    0.02,     "r_{#theta} [rad]")
#binning["Residual_qop"] =  acts.Axis.regular(200, -0.2,    0.2,     "r_{q/p} [c/GeV]")
#binning["Residual_pt_o_pt"] =  acts.Axis.regular(100, -0.1,    0.1,     "r_{p_{T}}/p_{T}")
#resPlotToolConfig.varBinning = binning


resPlotToolConfig.varBinning = {
    "Eta":             acts.Axis.regular(80,   -4,    4,     "#eta"),
    "Phi":             acts.Axis.regular(70,   -3.5,  3.5,   "#phi"),
    "Pt":              acts.Axis.regular(1000,  0,    50,    "pT [GeV/c]"),
    "Residual_phi":    acts.Axis.regular(200,  -0.02, 0.02,  "r_{#phi} [rad]"),
    "Residual_theta":  acts.Axis.regular(200,  -0.02, 0.02,  "r_{#theta} [rad]"),
    "Residual_qop":    acts.Axis.regular(200,  -0.2,  0.2,   "r_{q/p} [c/GeV]"),
    "Residual_pt_o_pt":acts.Axis.regular(100,  -0.1,  0.1,   "r_{p_{T}}/p_{T}"),
    "Pull":            acts.Axis.regular(100,  -5,    5,     "pull"),
    "Residual_d0":     acts.Axis.regular(100,  -0.5,  0.5,   "r_{d0} [mm]"),
    "Residual_z0":     acts.Axis.regular(100,  -0.5,  0.5,   "r_{z0} [mm]"),
    "Residual_t":      acts.Axis.regular(100,  -1000, 1000,  "r_{t} [s]"),
}



# Duplication Plots
duplicationPlotToolConfig = acts.examples.root.DuplicationPlotToolConfig()
duplicationPlotToolConfig.varBinning = {
    "Eta": acts.Axis.regular(80, -4,    4,     "#eta"),
    "Phi": acts.Axis.regular(100, -3.15, 3.15, "#phi"),
    "Pt":  acts.Axis.regular(1000,  0,    50,   "pT [GeV/c]"),
    "Num": acts.Axis.regular(30, -0.5,  29.5,  "N"),
}

# Efficiency Plots
effPlotToolConfig = acts.examples.root.EffPlotToolConfig()
binning = effPlotToolConfig.varBinning
binning["Eta"] = acts.Axis.regular(80, -4,    4,     "#eta")
binning["Pt"]  = acts.Axis.regular(1000,  0,    50,   "pT [GeV/c]")
effPlotToolConfig.varBinning = binning
effPlotToolConfig.minTruthPt = 0.15

# Fake rate Plots
fakePlotToolConfig = acts.examples.root.FakePlotToolConfig()
binning = fakePlotToolConfig.varBinning
binning["Eta"] = acts.Axis.regular(80, -4,    4,     "#eta")
binning["Pt"] = acts.Axis.regular(1000,  0,    50,   "pT [GeV/c]")
fakePlotToolConfig.varBinning = binning

# TrackQuality Plots
trackQualityPlotToolConfig = acts.examples.root.TrackQualityPlotToolConfig()
binning = trackQualityPlotToolConfig.varBinning
binning["Eta"] = acts.Axis.regular(80, -4,    4,     "#eta")
binning["Pt"] = acts.Axis.regular(1000,  0,    50,   "pT [GeV/c]")
trackQualityPlotToolConfig.varBinning = binning


# Track Summary Plots
trackSummaryPlotToolConfig = acts.examples.root.TrackSummaryPlotToolConfig()
binning = trackSummaryPlotToolConfig.varBinning
binning["Eta"] = acts.Axis.regular(80, -4,    4,     "#eta")
binning["Pt"] = acts.Axis.regular(1000,  0,    50,   "pT [GeV/c]")
trackSummaryPlotToolConfig.varBinning = binning
