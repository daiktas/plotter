from style import *
import style
import ROOT
import math
import os 
import json
import yaml
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-v", "--variable", action="store", dest="variable")
parser.add_argument("-r", "--region", action="store", dest="region")
parser.add_argument("-y", "--year", action="store", dest="year")
args = parser.parse_args()
variable_number = int(args.variable)
region_number = int(args.region)
year = str(args.year)

with open("variables.yaml") as yaml_file:
    lines = yaml_file.readlines()
    line = lines[variable_number]
    variable_infos = yaml.load(line, Loader=yaml.FullLoader)[0]
    print(variable_infos[0])

with open("regions.yaml") as yaml_file:
    lines = yaml_file.readlines()
    line = lines[region_number]
    category_infos = yaml.load(line, Loader=yaml.FullLoader)[0]
    category = category_infos[0]
    weight = category_infos[1]
    draw_text = category_infos[2]
    print(category, weight, draw_text)

print(year)


# path to processed nanoAOD ntuples
ntuple_path = "/vols/cms/vc1117/LLP/nanoAOD_friends/HNL/26May20_30GeV"
lumi = {"2016": 35.88, "2017": 41.53, "2018": 59.74}

def find_xsec(path, xsecs):
    for key, val in xsecs.items():
        if key in path:     
            return val

with open("/vols/build/cms/LLP/xsec.yaml") as yaml_file:
    xsecs = yaml.load(yaml_file, Loader=yaml.FullLoader)

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
        self.leg = makeLegend(0.70, 0.70, 0.90, 0.88)
        #self.leg = makeLegend(0.30, 0.70, 0.50, 0.88)
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
        self.hist_ratio.Draw("P SAME")
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
        makeText(0.15, 0.75, 0.3, 0.75, "#mu_{1} #mu_{2} (OS)" )
        makeLumiText(0.8, 0.97, lumi[year], year)
        canvas.SaveAs("/vols/cms/vc1117/LLP/plots_30GeV/"+suffix+self.varexp+"_"+year+"_.pdf")
        canvas.SaveAs("/vols/cms/vc1117/LLP/plots_30GeV/"+suffix+self.varexp+"_"+year+"_.png")

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
                self.sum_weight += yields[path]
        self.rdf = ROOT.RDataFrame("Friends", self.file_list) \
                .Filter("nselectedJets_nominal > 0") \
                .Filter("dimuon_mass > 20") \
                .Filter("dimuon_charge == -1") \
                .Filter("nlooseMuons == 1") \
                .Filter("ntightMuon == 1") \
                .Filter("nlepJet_nominal == 1")
        if self.isMC:
            self.rdf = self.rdf.Define("weightLumi", "IsoMuTrigger_weight_trigger_nominal*MET_filter*tightMuon_weight_iso_nominal*tightMuon_weight_id_nominal*looseMuons_weight_id_nominal*puweight*genweight*%s*1000.0*%s/%s" %(lumi[year], find_xsec(path, xsecs), self.sum_weight))
        else:
            self.rdf = self.rdf.Define("weightLumi", "1")
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
        print(args)
        for i, rdf in enumerate(self.rdfs):
            if i == 0:
                rdf = rdf.Define("var", varexp)
                hist = rdf.Histo1D(args, "var", weight)
            else:
                rdf = rdf.Define("var", varexp)
                tmp_hist = rdf.Histo1D(args, "var", weight)
                hist.Add(tmp_hist.GetValue())
        hist.GetXaxis().SetTitle(args[1])
        hist.SetLineColor(ROOT.TColor.GetColor(colorscale(self.color, 0.8)))
        hist.SetFillColor(ROOT.TColor.GetColor(colorscale(self.color, 1)))
        hist.SetLineWidth(2)
        self.hists.append(hist.Clone())
        return self.hists[-1]

with open(os.path.join("/vols/build/cms/LLP/yields_200311", year, "eventyields.json")) as json_file:
    yields = json.load(json_file)

if year == "2016":
    qcd_15to20 = Sample("qcd_15to20", ["QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"], year=year) 
    qcd_20to30 = Sample("qcd_20to30", ["QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"], year=year)
    qcd_30to50 = Sample("qcd_30to50", ["QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"], year=year)
    qcd_50to80 = Sample("qcd_50to80", ["QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"], year=year)
    qcd_80to120 = Sample("qcd_80to120", ["QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"], year=year)
    qcd_120to170 = Sample("qcd_120to170", ["QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"], year=year)
    qcd_170to300 = Sample("qcd_170to300", ["QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"], year=year)
    qcd_300to470 = Sample("qcd_300to470", ["QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"], year=year)
    qcd_470to600 = Sample("qcd_470to600", ["QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"], year=year)
    qcd_600to800 = Sample("qcd_600to800", ["QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"], year=year)
    qcd_800to1000 = Sample("qcd_800to1000", ["QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"], year=year)
    qcd_1000toInf = Sample("qcd_1000toInf", ["QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"], year=year)

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
    qcd_600to800 = Sample("qcd_600to800", ["QCD_Pt-600to800_MuEnrichedPt5_TuneCP5_13TeV_pythia8-"+year], year=year)
    qcd_800to1000 = Sample("qcd_800to1000", ["QCD_Pt-800to1000_MuEnrichedPt5_TuneCP5_13TeV_pythia8-"+year], year=year)
    qcd_1000toInf = Sample("qcd_1000toInf", ["QCD_Pt-1000toInf_MuEnrichedPt5_TuneCP5_13TeV_pythia8-"+year], year=year)

if year == "2016":
    ttsemilep = Sample("ttsemilep", ["TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8-2016"], year=year)
    ttdilep = Sample("ttdilep", ["TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8-2016"], year=year)
    st_t_top = Sample("st_t_top", ["ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1-2016"], year=year)
    st_t_antitop = Sample("st_t_antitop", ["ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8-2016"], year=year)
    st_tW_top = Sample("st_tW_top", ["ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1-2016"], year=year)
    st_tW_antitop = Sample("st_tW_antitop", ["ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1-2016"], year=year)

elif year == "2017":
    ttsemilep = Sample("ttsemilep", ["TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_ext1-2017"], year=year)
    ttdilep = Sample("ttdilep", ["TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8-2017"], year=year)
    st_t_top = Sample("st_t_top", ["ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8-2017"], year=year)
    st_t_antitop = Sample("st_t_antitop", ["ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8-2017"], year=year)
    st_tW_top = Sample("st_tW_top", ["ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8-2017"], year=year)
    st_tW_antitop = Sample("st_tW_antitop", ["ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8-2017"], year=year)

elif year == "2018":
    ttsemilep = Sample("ttsemilep", ["TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8-2018"], year=year)
    ttdilep = Sample("ttdilep", ["TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8-2018"], year=year)
    st_t_top = Sample("st_t_top", ["ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8-2018"], year=year)
    st_t_antitop = Sample("st_t_antitop", ["ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8-2018"], year=year)
    st_tW_top = Sample("st_tW_top" ,["ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8-2018"], year=year)
    st_tW_antitop = Sample("st_tW_antitop", ["ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8-2018"], year=year)

if year == "2016":
    w0jets = Sample("w0jets", ["WToLNu_0J_13TeV-amcatnloFXFX-pythia8-2016"], year=year)
    w1jets = Sample("w1jets", ["WToLNu_1J_13TeV-amcatnloFXFX-pythia8-2016"], year=year)
    w2jets = Sample("w2jets", ["WToLNu_2J_13TeV-amcatnloFXFX-pythia8-ext4-2016"], year=year)
    dy10to50 = Sample("dy10to50",  ["DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-2016"], year=year)
    dy50 = Sample("dy50",  ["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-ext2-2016"], year=year)

elif year == "2017":
    w0jets = Sample("w0jets", ["WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2017"], year=year)
    w1jets = Sample("w1jets", ["WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8-ext1-2017"], year=year)
    w2jets = Sample("w2jets", ["WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2017"], year=year)
    dy10to50 = Sample("dy10to50",  ["DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8-2017"], year=year)
    dy50 = Sample("dy50",  ["DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8-2017"], year=year)

elif year == "2018":
    w0jets = Sample("w0jets", ["WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2018"], year=year)
    w1jets = Sample("w1jets", ["WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2018"], year=year)
    w2jets = Sample("w2jets", ["WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2018"], year=year)     
    dy10to50 = Sample("dy10to50",  ["DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8-2018"], year=year)
    dy50 = Sample("dy50",  ["DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8-2018"], year=year)

wjets = Process("W+Jets", "W+jets", "#388e3c")
wjets.add(w0jets, w1jets, w2jets)

dyjets = Process("DY+Jets", "DY+Jets", "#1976d2")
dyjets.add(dy10to50, dy50)

tt = Process("ttbar", "t#bar{t}/st",  "#ef5350")
tt.add(ttsemilep, ttdilep)
tt.add(st_t_top, st_t_antitop, st_tW_top, st_tW_antitop)

qcd = Process("qcd", "QCD", "#bdbdbd")
qcd.add(qcd_15to20, qcd_20to30, qcd_30to50, qcd_50to80, qcd_80to120, qcd_120to170, qcd_170to300, qcd_300to470, qcd_470to600, qcd_600to800, qcd_800to1000, qcd_1000toInf)

if year == "2016":
    runb = Sample("RunBv2", ["SingleMuon_Run2016B_ver2"], isMC=False)
    runc = Sample("RunC", ["SingleMuon_Run2016C"], isMC=False)
    rund = Sample("RunD", ["SingleMuon_Run2017D"], isMC=False)
    rune = Sample("RunE", ["SingleMuon_Run2016E"], isMC=False)
    runf = Sample("RunF", ["SingleMuon_Run2016F"], isMC=False)
    rung = Sample("RunG", ["SingleMuon_Run2016G"], isMC=False)
    runh = Sample("RunH", ["SingleMuon_Run2016H"], isMC=False)

elif year == "2017":
    runb = Sample("RunB", ["SingleMuon_Run2017B"], isMC=False)
    runc = Sample("RunC", ["SingleMuon_Run2017C"], isMC=False)
    rund = Sample("RunD", ["SingleMuon_Run2017D"], isMC=False)
    rune = Sample("RunE", ["SingleMuon_Run2017E"], isMC=False)
    runf = Sample("RunF", ["SingleMuon_Run2017F"], isMC=False)

elif year == "2018":
    runa = Sample("RunA", ["SingleMuon_Run2018A"], isMC=False)
    runb = Sample("RunB", ["SingleMuon_Run2018B"], isMC=False)
    runc = Sample("RunC", ["SingleMuon_Run2018C"], isMC=False)
    rund = Sample("RunD", ["SingleMuon_Run2018D"], isMC=False)

data = Process("data", "data", "#000000")

if year == "2016":
    data.add(runb, runc, rund, rune, runf, rung, runh)
elif year == "2017":
    data.add(runb, runc, rund, rune, runf)
elif year == "2018":
    data.add(runa, runb, runc, rund)

processes = [qcd, tt, wjets, dyjets, data]

variable = Variable(*variable_infos)
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
    variable.Add(process.Histo1D(variable.args, variable.varexp, category), process.title, isSignal=isSignal, isData=isData)
variable.Draw(category, "hist", draw_text, year=year)
