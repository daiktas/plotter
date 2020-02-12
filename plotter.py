from style import *
import style
import ROOT
import math
import os 
import json
import yaml

# path to processed nanoAOD ntuples
ntuple_path = "/vols/cms/vc1117/LLP/nanoAOD_friends/HNL/06Feb20/"
lumi = 35.88

def find_xsec(path, xsecs):
    for key, val in xsecs.items():
        if key in path:     
            return val

# Read in event yields and cross-sections for normalisation
with open("eventyields.json") as json_file:
    yields = json.load(json_file)

with open("xsecs.yaml") as yaml_file:
    xsecs = yaml.load(yaml_file, Loader=yaml.FullLoader)

# Define categories by cuts for plotting
categories = {}
categories["mu1mu2jet_"] = "(ntightMuons == 1)*(nlooseMuons == 1)"
categories["mu1mu2jet_dimuoncut_"] = "(ntightMuons == 1)*(nlooseMuons == 1)*(lepJet_deltaR<0.4)*(dimuon_mass < 80.)*(dimuon_mass > 30.)*(dimuon_deltaR>1.)*(dimuon_deltaR<5.)"

# This class is responsible for making the histogram and plotting it for a given variable
class Variable:
    def __init__(self, varexp, name, nbins, xmin, xmax, logy=True):
        self.varexp = varexp
        self.args = (varexp, varexp, nbins, xmin, xmax)
        self.stack = ROOT.THStack(varexp, varexp)
        self.signals = []
        self.name = name
        self.logy = logy
        self.leg = makeLegend(0.4, 0.70, 0.7, 0.88)
        self.leg.SetTextSize(self.leg.GetTextSize()*0.8)
    def Add(self, hist, title, isSignal=False):
        if isSignal:
            self.signals.append(hist)
            hist.SetLineStyle(len(self.signals))
            self.leg.AddEntry(hist, title, "l")
        else:
            self.stack.Add(hist)
            self.leg.AddEntry(hist, title, "f")
    def Draw(self, suffix, opt):
        print ("plotting "+self.varexp)
        self.canvas = makeCanvas(name=self.varexp)
        self.canvas.Draw()
        self.canvas.SetBottomMargin(0.15)
        self.canvas.SetLeftMargin(0.1)
        self.canvas.SetTopMargin(0.08)
        self.stack.Draw(opt)
        self.stack.GetXaxis().SetTitle(self.name)
        self.stack.SetMinimum(1)
        if self.logy:
            self.stack.SetMaximum(self.stack.GetMaximum()*1000)
            self.canvas.SetLogy()
        else:
            self.stack.SetMaximum(self.stack.GetMaximum()*1.6)
        for signal in self.signals:
            signal.SetFillStyle(0)
            signal.Draw("HIST SAME")
        self.leg.Draw("SAME")
        makeCMSText(0.13, 0.97,additionalText="Simulation Preliminary")
        makeText(0.13, 0.78, 0.4, 0.78, "#mu_{1}#mu_{2}")
        makeLumiText(0.7, 0.97)
        self.canvas.Modified()
        self.canvas.SaveAs("plots/"+suffix+self.varexp+".pdf")

# This class prepares a given sample by scaling to int. luminosity
class Sample:
    def __init__(self, name, paths):
        self.paths = paths
        self.name = name
        self.file_list = ROOT.std.vector('string')()
        self.sum_weight = 0
        for path in self.paths:
            for f in os.listdir(os.path.join(ntuple_path, path)):
                self.file_list.push_back(os.path.join(ntuple_path, path, f))
            self.sum_weight += yields[path]
        self.rdf = ROOT.RDataFrame("Friends", self.file_list)
        self.rdf = self.rdf.Define("weightLumi", "genweight*%s*1000.0*%s/%s" %(lumi, find_xsec(path, xsecs), self.sum_weight))
        for category, weight in categories.items():
            self.rdf = self.rdf.Define(category, "weightLumi*%s" %(weight))
        self.hists = []
        print("RDF has entries: "+str(self.rdf.Count().GetValue()))

# A process is a combination of several "Samples" which are all added up internally
class Process:
    def __init__(self, name, title, color):
        self.name = name
        self.title = title
        self.color = color	
        self.hists = []
        self.rdfs = []

    def add(self, *args):
        for arg in args:
            self.rdfs.append(arg.rdf)

    def Histo1D(self, args, varexp, weight):
        for i, rdf in enumerate(self.rdfs):
            if i == 0:
                hist = rdf.Histo1D(args, varexp, weight)
            else:
                tmp_hist = rdf.Histo1D(args, varexp, weight)
                hist.Add(tmp_hist.GetValue())

        hist.GetXaxis().SetTitle(args[1])
        hist.SetLineColor(ROOT.TColor.GetColor(colorscale(self.color, 0.8)))
        hist.SetFillColor(ROOT.TColor.GetColor(colorscale(self.color, 1)))
        hist.SetLineWidth(2)
        self.hists.append(hist.Clone())
        return self.hists[-1]


w0jets = Sample("wjets", ["WToLNu_0J_13TeV-amcatnloFXFX-pythia8-2016"])
w1jets = Sample("w1jets", ["WToLNu_1J_13TeV-amcatnloFXFX-pythia8-2016"])
w2jets = Sample("w2jets", ["WToLNu_2J_13TeV-amcatnloFXFX-pythia8-ext4-2016"])
dy10to50 = Sample("dy10to50",  ["DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-2016"])
dy50 = Sample("dy50",  ["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-ext2-2016"])
ttsemilep = Sample("tt", ["TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8-2016"])
qcd_15to20 = Sample("qcd_15to20", ["QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"]) 
qcd_20to30 = Sample("qcd_20to30", ["QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"])
qcd_30to50 = Sample("qcd_30to50", ["QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"])
qcd_50to80 = Sample("qcd_50to80", ["QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"])
qcd_80to120 = Sample("qcd_80to120", ["QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"])
qcd_120to170 = Sample("qcd_120to170", ["QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"])
qcd_170to300 = Sample("qcd_170to300", ["QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia-2016"])
qcd_170to300 = Sample("qcd_170to300", ["QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia-2016"])
qcd_300to470 = Sample("qcd_300to470", ["QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"])
qcd_470to600 = Sample("qcd_470to600", ["QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"])
qcd_600to800 = Sample("qcd_600to800", ["QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"])
qcd_800to1000 = Sample("qcd_800to1000", ["QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"])
qcd_1000toInf = Sample("qcd_1000toInf", ["QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"])
hnlM4_V0p00183575597507 = Sample("HNL2", ["HeavyNeutrino_lljj_M-4_V-0_00183575597507_mu_Dirac_Moriond17_aug2018_miniAODv3-2016"])
hnlM8_V0p000415932686862 = Sample("HNL", ["HeavyNeutrino_lljj_M-8_V-0_000415932686862_mu_Dirac_Moriond17_aug2018_miniAODv3-2016"])

wjets = Process("W+Jets", "W+jets", "#388e3c")
wjets.add(w0jets, w1jets, w2jets)
dyjets = Process("DY+Jets", "DY+Jets", "#1976d2")
dyjets.add(dy10to50, dy50)
tt = Process("ttbar", "t#bar{t}",  "#ef5350")
tt.add(ttsemilep)
hnl1 = Process("HNL1", "m_{N} = 8 GeV, |V_{#mu}|^{2} = 4.2#times10^{-7}", "#087858")
hnl1.add(hnlM8_V0p000415932686862)
hnl2 = Process("HNL2", "m_{N} = 4 GeV, |V_{#mu}|^{2} = 3.4#times10^{-6}", "#087858")
hnl2.add(hnlM4_V0p00183575597507)
qcd = Process("qcd", "QCD", "#bdbdbd")
qcd.add(qcd_15to20, qcd_20to30, qcd_30to50, qcd_50to80, qcd_80to120, qcd_120to170, qcd_170to300, qcd_300to470, qcd_470to600, qcd_600to800, qcd_800to1000, qcd_1000toInf)

processes = [wjets, dyjets, tt, qcd, hnl1, hnl2]
# add QCD


for suffix, weight in categories.items():
    with open("variables.yaml") as yaml_file:
        variables = yaml.load(yaml_file, Loader=yaml.FullLoader)
        print(variables)
    for variable in variables:
        variable = Variable(*variable)
        print(variable.name)
        if os.path.exists("plots/"+suffix+variable.varexp+".pdf"):
            print("skipping")
            continue
        for process in processes:
            print(process.name)
            if "HNL" in process.name:
                isSignal = True
            else:
                isSignal = False
            variable.Add(process.Histo1D(variable.args, variable.varexp, weight), process.title, isSignal=isSignal)
        variable.Draw(suffix, "hist")
