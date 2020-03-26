from style import *
import style
import ROOT
import math
import os 
import json
import yaml
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-l", "--line", action="store", dest="line")
args = parser.parse_args()
line_number = int(args.line)

with open("variables.yaml") as yaml_file:
    lines = yaml_file.readlines()
    line = lines[line_number]
    print(line)
    variable_infos = yaml.load(line, Loader=yaml.FullLoader)[0]
    print(variable_infos)

# path to processed nanoAOD ntuples
ntuple_path = "/vols/cms/vc1117/LLP/nanoAOD_friends/HNL/25Mar20/"
lumi = {"2016": 35.88, "2017": 41.53, "2018": 58.83}

def find_xsec(path, xsecs):
    for key, val in xsecs.items():
        if key in path:     
            return val


with open("xsecs.yaml") as yaml_file:
    xsecs = yaml.load(yaml_file, Loader=yaml.FullLoader)

# Define categories by cuts for plotting
categories = {}
categories["mu1mu2jet_"] = ["(ntightMuon == 1)*(nlooseMuons == 1)*(nlepJet == 1)", " "]
#categories["mu1mu2jet_met_CR_"] = ["(ntightMuon == 1)*(nlooseMuons == 1)*(lepJet_deltaR<0.4)*(MET_pt > 100.)*(dimuon_mass < 80.)", "p_{T}^{miss} > 100 GeV, m(#mu#mu) < 80 GeV"]
#categories["mu1mu2jet_mt_CR_"] = ["(ntightMuon == 1)*(nlooseMuons == 1)*(lepJet_deltaR<0.4)*(Jet_Muon_MET_Mu_mT > 90.)*(dimuon_mass < 80.)", "m_{T} > 90 GeV, m(#mu#mu) < 80 GeV"]
#categories["mu1mu2jet_ht_CR_"] = ["(ntightMuon == 1)*(nlooseMuons == 1)*(lepJet_deltaR<0.4)*(Jet_Muon_ht > 250.)*(dimuon_mass < 80.)", "H_{T} > 250 GeV, m(#mu#mu) < 80 GeV"]
#categories["mu1mu2jet_DY_CR_"] = ["(ntightMuon == 1)*(nlooseMuons == 1)*(lepJet_deltaR<0.4)*(dimuon_mass > 90.)*(dimuon_mass < 110.)", "90 < m(#mu#mu) < 110 GeV"]

# This class is responsible for making the histogram and plotting it for a given variable
class Variable:
    def __init__(self, varexp, name, nbins, xmin, xmax, logy=True):
        self.varexp = varexp
        self.args = (varexp, varexp, nbins, xmin, xmax)
        self.stack = ROOT.THStack(varexp, varexp)
        self.sumMC = ROOT.TH1F(name, name, nbins, xmin, xmax)
        self.signals = []
        self.data = None
        self.name = name
        self.logy = logy
        #self.leg = makeLegend(0.70, 0.70, 0.90, 0.88)
        self.leg = makeLegend(0.30, 0.70, 0.50, 0.88)
        self.leg.SetTextSize(self.leg.GetTextSize()*0.8)
        self.xmin = xmin
        self.xmax = xmax
    def Add(self, hist, title, isSignal=False, isData=False):
        hist.SetDirectory(0)
        if isSignal:
            self.signals.append(hist)
            hist.SetLineStyle(len(self.signals))
            self.leg.AddEntry(hist, title, "l")
        elif isData:
            self.data = hist
            self.leg.AddEntry(hist, title, "p")
        else:
            self.stack.Add(hist)
            self.sumMC.Add(hist)
            self.leg.AddEntry(hist, title, "f")
    def Draw(self, suffix, opt, draw_text, year="2016"):
        print ("plotting "+self.varexp)
        canvas = makeCanvas(name=self.varexp)
        upperPad = ROOT.TPad("upperPad", "upperPad", 0, 0.33, 1, 1)
        lowerPad = ROOT.TPad("lowerPad", "lowerPad", 0, 0, 1, 0.33)
        upperPad.SetBottomMargin(0.00001)
        upperPad.SetBorderMode(0)
        upperPad.SetTopMargin(0.15)
        lowerPad.SetTopMargin(0.00001)
        lowerPad.SetBottomMargin(0.4)
        lowerPad.SetBorderMode(0)
        canvas.SetBottomMargin(0.2)
        canvas.SetTopMargin(0.1)
        upperPad.Draw()
        lowerPad.Draw()
        upperPad.cd()

        self.stack.Draw(opt)
        self.stack.SetMinimum(1)
        if self.logy:
            self.stack.SetMaximum(self.stack.GetMaximum()*1000)
            upperPad.SetLogy()
        else:
            self.stack.SetMaximum(self.stack.GetMaximum()*1.6)

        for signal in self.signals:
            signal.SetFillStyle(0)
            signal.Draw("HIST SAME")

        if self.data is not None:
            self.data.Draw("P SAME")
            self.hist_ratio = self.data.Clone("ratio histogram")
            self.hist_ratio.Divide(self.sumMC)

        lowerPad.cd()
        axis = self.sumMC.Clone("axis")
        axis.SetMinimum(0.25)
        axis.SetMaximum(1.75)
        axis.GetXaxis().SetTitle(self.name)
        axis.GetXaxis().SetTitleOffset(2.5)
        axis.Draw("AXIS")

        rootObj = []

        line = ROOT.TLine(self.xmin, 1, self.xmax, 1)
        line.Draw("SAME")
        for ibin in range(self.sumMC.GetNbinsX()):
            c = self.sumMC.GetBinCenter(ibin+1)
            w = self.sumMC.GetBinWidth(ibin+1)
            m = self.sumMC.GetBinContent(ibin+1)
            if m > 0.0:
                h = min(self.sumMC.GetBinError(ibin+1)/m, 0.399)
                box = ROOT.TBox(c-0.5*w, 1-h, c+0.5*w, 1+h)
                box.SetFillStyle(3345)
                box.SetLineColor(ROOT.kGray+1)
                box.SetFillColor(ROOT.kGray)
                rootObj.append(box)
                box.Draw("SameF")
                box2 = ROOT.TBox(c-0.5*w, 1-h, c+0.5*w, 1+h)
                box2.SetFillStyle(0)
                box2.SetLineColor(ROOT.kGray+1)
                box2.SetFillColor(ROOT.kGray)
                rootObj.append(box2)
                box2.Draw("SameL")
        canvas.cd()
        self.leg.Draw("SAME")
        makeCMSText(0.13, 0.97,additionalText="Simulation Preliminary")
        makeText(0.12, 0.80, 0.3, 0.80, draw_text)
        makeLumiText(0.8, 0.97, lumi=lumi[year])
        canvas.SaveAs("plots/"+suffix+self.varexp+"_"+year+"_.pdf")

# This class prepares a given sample by scaling to int. luminosity
class Sample:
    def __init__(self, name, paths, isMC=True, year="2016"):
        self.paths = paths
        self.name = name
        self.file_list = ROOT.std.vector('string')()
        self.sum_weight = 0
        self.isMC = isMC
        for path in self.paths:
            for f in os.listdir(os.path.join(ntuple_path, path)):
                self.file_list.push_back(os.path.join(ntuple_path, path, f))
            if self.isMC:
                self.sum_weight += yields[path]["weighted"]
        self.rdf = ROOT.RDataFrame("Friends", self.file_list)
        if self.isMC:
            # print(path, xsecs)
            self.rdf = self.rdf.Define("weightLumi", "IsoMuTrigger_weight_trigger_nominal*MET_filter*tightMuon_weight_iso_nominal*tightMuon_weight_id_nominal*looseMuons_weight_id_nominal*puweight*genweight*%s*1000.0*%s/%s" %(lumi[year], find_xsec(path, xsecs), self.sum_weight))
        else:
            self.rdf = self.rdf.Define("weightLumi", "1")
        for category, properties in categories.items():
            weight = properties[0]
            self.rdf = self.rdf.Define(category, "weightLumi*%s" %(weight))

        self.hists = []
        print("RDF "+name+ " has entries: "+str(self.rdf.Count().GetValue()))

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

# for year in ["2018", "2017", "2016"]:
for year in ["2016"]:
    # Read in event yields and cross-sections for normalisation
    with open("eventyields"+year+".json") as json_file:
        yields = json.load(json_file)

    if year == "2016":
        qcd_15to20 = Sample("qcd_15to20", ["QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-"+year], year=year) 
        qcd_20to30 = Sample("qcd_20to30", ["QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-"+year], year=year)
        qcd_30to50 = Sample("qcd_30to50", ["QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-"+year], year=year)
        qcd_50to80 = Sample("qcd_50to80", ["QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-"+year], year=year)
        qcd_80to120 = Sample("qcd_80to120", ["QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-"+year], year=year)
        qcd_120to170 = Sample("qcd_120to170", ["QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-"+year], year=year)
        qcd_170to300 = Sample("qcd_170to300", ["QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia-"+year], year=year)
        qcd_300to470 = Sample("qcd_300to470", ["QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-"+year], year=year)
        qcd_470to600 = Sample("qcd_470to600", ["QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-"+year], year=year)
        qcd_600to800 = Sample("qcd_600to800", ["QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-"+year])
        qcd_800to1000 = Sample("qcd_800to1000", ["QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-"+year], year=year)
        qcd_1000toInf = Sample("qcd_1000toInf", ["QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-"+year], year=year)
        ttsemilep = Sample("ttsemilep", ["TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8-"+year], year=year)
        ttdilep = Sample("ttdilep", ["TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8-"+year], year=year)

    else:
        qcd_15to20 = Sample("qcd_15to20", ["QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8-"+year], year=year) 
        qcd_20to30 = Sample("qcd_20to30", ["QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8-"+year], year=year)
        qcd_30to50 = Sample("qcd_30to50", ["QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8-"+year], year=year)
        qcd_50to80 = Sample("qcd_50to80", ["QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8-"+year], year=year)
        qcd_80to120 = Sample("qcd_80to120", ["QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8-"+year], year=year)
        qcd_120to170 = Sample("qcd_120to170", ["QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8-"+year], year=year)
        qcd_170to300 = Sample("qcd_170to300", ["QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8-"+year], year=year)
        qcd_300to470 = Sample("qcd_300to470", ["QCD_Pt-300to470_MuEnrichedPt5_TuneCP5_13TeV_pythia8-"+year], year=year)
        qcd_470to600 = Sample("qcd_470to600", ["QCD_Pt-470to600_MuEnrichedPt5_TuneCP5_13TeV_pythia8-"+year], year=year)
        qcd_600to800 = Sample("qcd_600to800", ["QCD_Pt-600to800_MuEnrichedPt5_TuneCP5_13TeV_pythia8-"+year])
        qcd_800to1000 = Sample("qcd_800to1000", ["QCD_Pt-800to1000_MuEnrichedPt5_TuneCP5_13TeV_pythia8-"+year], year=year)
        qcd_1000toInf = Sample("qcd_1000toInf", ["QCD_Pt-1000toInf_MuEnrichedPt5_TuneCP5_13TeV_pythia8-"+year], year=year)

        if year == "2017":
            ttsemilep = Sample("ttsemilep", ["TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_ext1-"+year], year=year)
            ttdilep = Sample("ttdilep", ["TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8-"+year], year=year)

        else:
            ttsemilep = Sample("ttsemilep", ["TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8-"+year], year=year)
            ttdilep = Sample("ttdilep", ["TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8-"+year], year=year)
            
    if year == "2016":
        #runb = Sample("RunBv2", ["SingleMuon_Run2016B_ver2"], isMC=False)
        #rune = Sample("RunE", ["SingleMuon_Run2016E"], isMC=False)
        #runf = Sample("RunF", ["SingleMuon_Run2016F"], isMC=False)
        #rung = Sample("RunG", ["SingleMuon_Run2016G"], isMC=False)
        #runh = Sample("RunH", ["SingleMuon_Run2016H"], isMC=False)
        w0jets = Sample("w0jets", ["WToLNu_0J_13TeV-amcatnloFXFX-pythia8-2016"], year=year)
        w1jets = Sample("w1jets", ["WToLNu_1J_13TeV-amcatnloFXFX-pythia8-2016"], year=year)
        w2jets = Sample("w2jets", ["WToLNu_2J_13TeV-amcatnloFXFX-pythia8-ext4-2016"], year=year)
        dy10to50 = Sample("dy10to50",  ["DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-2016"], year=year)
        dy50 = Sample("dy50",  ["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-ext2-2016"], year=year)

    elif year == "2017":
        #runb = Sample("RunB", ["SingleMuon_Run2017B"], isMC=False)
        #rune = Sample("RunE", ["SingleMuon_Run2017E"], isMC=False)
        #runf = Sample("RunF", ["SingleMuon_Run2017F"], isMC=False)
        w0jets = Sample("w0jets", ["WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2017"], year=year)
        w1jets = Sample("w1jets", ["WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8-ext1-2017"], year=year)
        w2jets = Sample("w2jets", ["WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2017"], year=year)
        dy0jets = Sample("dy0jets", ["DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2017"], year=year)
        dy1jets = Sample("dy1jets", ["DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2017"], year=year)
        dy2jets = Sample("dy2jets", ["DYJetsToLL_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2017"], year=year)

    elif year == "2018":
        #runa = Sample("RunA", ["SingleMuon_Run2018A"], isMC=False)
        #runb = Sample("RunB", ["SingleMuon_Run2017B"], isMC=False)
        w0jets = Sample("w0jets", ["WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2018"], year=year)
        w1jets = Sample("w1jets", ["WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2018"], year=year)
        w2jets = Sample("w2jets", ["WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2018"], year=year)     
        dy0jets = Sample("dy0jets", ["DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2018"], year=year)
        dy1jets = Sample("dy1jets", ["DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2018"], year=year)
        dy2jets = Sample("dy2jets", ["DYJetsToLL_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2018"], year=year)

    '''
    runc = Sample("RunC", ["SingleMuon_Run"+year+"C"], isMC=False)
    rund = Sample("RunD", ["SingleMuon_Run"+year+"D"], isMC=False)
    data = Process("data", "data", "#000000")
    if year == "2016":
        data.add(runb, runc, rund, rune, runf, rung, runh)
    elif year == "2017":
        data.add(runb, runc, rund, rune, runf)
    elif year == "2018":
        data.add(runa, runb, runc, rund)
    #processes = [wjets, dyjets, tt, qcd, data]
    '''

    hnlM4_V0p00183575597507 = Sample("HNL2", ["HeavyNeutrino_lljj_M-4_V-0_00183575597507_mu_Dirac_Moriond17_aug2018_miniAODv3-2016"])
    hnlM8_V0p000415932686862 = Sample("HNL", ["HeavyNeutrino_lljj_M-8_V-0_000415932686862_mu_Dirac_Moriond17_aug2018_miniAODv3-2016"])
    hnl1 = Process("HNL1", "m_{N} = 8 GeV, |V_{#mu}|^{2} = 4.2#times10^{-7}, #times 1000", "#087858")
    hnl1.add(hnlM8_V0p000415932686862)
    hnl2 = Process("HNL2", "m_{N} = 4 GeV, |V_{#mu}|^{2} = 3.4#times10^{-6}, #times 1000", "#087858")
    hnl2.add(hnlM4_V0p00183575597507)

    wjets = Process("W+Jets", "W+jets", "#388e3c")
    wjets.add(w0jets, w1jets, w2jets)
    dyjets = Process("DY+Jets", "DY+Jets", "#1976d2")
    if year == "2016":
        dyjets.add(dy10to50, dy50)
    else:
        dyjets.add(dy0jets, dy1jets, dy2jets)

    tt = Process("ttbar", "t#bar{t}",  "#ef5350")
    tt.add(ttsemilep, ttdilep)
    qcd = Process("qcd", "QCD", "#bdbdbd")
    qcd.add(qcd_15to20, qcd_20to30, qcd_30to50, qcd_50to80, qcd_80to120, qcd_120to170, qcd_170to300, qcd_300to470, qcd_470to600, qcd_600to800, qcd_800to1000, qcd_1000toInf)

    processes = [wjets, dyjets, tt, qcd, hnl1, hnl2]


    for suffix, category in categories.items():
        weight = category[0]
        draw_text = category[1]
        print(weight)
        variable = Variable(*variable_infos)
        print(variable.args, variable.varexp)
        for process in processes:
            print(process.name)
            if "HNL" in process.name:
                isSignal = True
            else:
                isSignal = False
            if process.name == "data":
                isData = True
            else:
                isData = False
            variable.Add(process.Histo1D(variable.args, variable.varexp, suffix), process.title, isSignal=isSignal, isData=isData)
        variable.Draw(suffix, "hist", draw_text, year=year)
