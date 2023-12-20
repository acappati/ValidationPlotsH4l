#!/bin/env python3
### Example macro for filling standard histograms from H4l nanoAODs.
### Histograms are stored on a file and can then be plotted with

from __future__ import print_function
import math
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ZZAnalysis.NanoAnalysis.tools import getLeptons


pathMC2018 = "/eos/user/n/namapane/H4lnano/220420/" # FIXME: Use 2018 MC for the time being
pathMC2022 = '/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/231214_nano/MC2022/'
pathMC2022EE = '/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/231214_nano/MC2022EE/'
pathDATA = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/231209_nano/Data2022/"

ZmassValue = 91.1876

maxEntriesPerSample = 1e12 # Use only up to this number of events in each MC sample, for quick tests.



ROOT.TH1.SetDefaultSumw2()

####################################
def fillHistos(samplename, filename) :

    ## ZZMass
    h_ZZMass2 = ROOT.TH1F("ZZMass_2GeV_"+samplename,
                          "ZZMass_2GeV_"+samplename,65,70.,200.)
    h_ZZMass2.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h_ZZMass2.GetYaxis().SetTitle("Events / 2 GeV")
    
    h_ZZMass4 = ROOT.TH1F("ZZMass_4GeV_"+samplename,
                          "ZZMass_4GeV_"+samplename,233,70.,1002.)
    h_ZZMass4.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h_ZZMass4.GetYaxis().SetTitle("Events / 4 GeV")

    # h_ZZMass10 = ROOT.TH1F("ZZMass_10GeV_"+samplename,
    #                        "ZZMass_10GeV_"+samplename,93,70.,1000.)
    # h_ZZMass10.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    # h_ZZMass10.GetYaxis().SetTitle("Events / 10 GeV")

    ## Z1 and Z2 masses
    # Z1
    h_Z1Mass = ROOT.TH1F("Z1Mass_2GeV_"+samplename,
                         "Z1Mass_2GeV_"+samplename,40,40.,120.)
    h_Z1Mass.GetXaxis().SetTitle("m_{#it{Z1}} (GeV)")
    h_Z1Mass.GetYaxis().SetTitle("Events / 2 GeV")
    # Z2
    h_Z2Mass = ROOT.TH1F("Z2Mass_2GeV_"+samplename,
                         "Z2Mass_2GeV_"+samplename,60,0.,120.)
    h_Z2Mass.GetXaxis().SetTitle("m_{#it{Z2}} (GeV)")
    h_Z2Mass.GetYaxis().SetTitle("Events / 2 GeV")

    ## KD
    h_KD = ROOT.TH1F("KD_"+samplename,
                     "KD_"+samplename,10,0.,1.)
    h_KD.GetXaxis().SetTitle("#it{D}_{bkg}^{kin}")
    h_KD.GetYaxis().SetTitle("Events / 0.1")



    f = ROOT.TFile.Open(filename)

    event = f.Events
    event.SetBranchStatus("*", 0)
    event.SetBranchStatus("run", 1)
    event.SetBranchStatus("luminosityBlock", 1)
    event.SetBranchStatus("*Muon*", 1)
    event.SetBranchStatus("*Electron*", 1)
    event.SetBranchStatus("*ZZCand*", 1)
    event.SetBranchStatus("bestCandIdx", 1)
    event.SetBranchStatus("HLT_passZZ4l", 1)
    nEntries = event.GetEntries() 

    isMC = False
    if(samplename == "Data"):
        print("Data: sel=", nEntries)
    else:
        isMC = True
        event.SetBranchStatus("overallEventWeight",1)

        # Get sum of weights
        runs = f.Runs
        nRuns = runs.GetEntries()
        iRun = 0
        genEventCount = 0
        genEventSumw = 0.
        while iRun < nRuns and runs.GetEntry(iRun) :
            genEventCount += runs.genEventCount
            genEventSumw += runs.genEventSumw
            iRun +=1
        print (samplename, ": gen=", genEventCount, "sel=",nEntries, "sumw=", genEventSumw)
        #if nEntries>maxEntriesPerSample :
        #    genEventSumw = genEventSumw*maxEntriesPerSample/nEntries
        #    nEntries=maxEntriesPerSample
        #    print("   scaling to:", nEntries, "sumw=", genEventSumw )


    iEntry=0
    printEntries=max(5000,nEntries/10)
    while iEntry<nEntries and event.GetEntry(iEntry):
        iEntry+=1
        if iEntry%printEntries == 0 : print("Processing", iEntry)

        bestCandIdx = event.bestCandIdx

        # Check that the event contains a selected candidate, and that
        # passes the required triggers (which is necessary for samples
        # processed with TRIGPASSTHROUGH=True)
        if(bestCandIdx != -1 and event.HLT_passZZ4l): 
            weight = 1.
            ZZs = Collection(event, 'ZZCand')
            theZZ = ZZs[bestCandIdx]        
            if isMC : 
                weight = (event.overallEventWeight*theZZ.dataMCWeight/genEventSumw)
            ## ZZmass
            m4l=theZZ.mass
            h_ZZMass2.Fill(m4l,weight)
            h_ZZMass4.Fill(m4l,weight)
            # h_ZZMass10.Fill(m4l,weight)
            ## Z1Mass
            mZ1=theZZ.Z1mass
            h_Z1Mass.Fill(mZ1,weight)
            ## Z2Mass
            mZ2=theZZ.Z2mass
            h_Z2Mass.Fill(mZ2,weight)
            ## KD
            KD=theZZ.KD
            h_KD.Fill(KD,weight)
            # Example on how to get the four leptons of the candidates, ordered as
            # [Z1l1, Z2l2, Z2l1, Z2l2]
            #leps = getLeptons(theZZ, event)
            #print(leps[3].pt)
        
    f.Close()
    
    return h_ZZMass2,h_ZZMass4,h_Z1Mass,h_Z2Mass,h_KD


def runMC(outFile): 

    if '2018' in outFile:
        samples = [
            # ggZZ from 2018
            dict(name = "ggTo4mu",filename = pathMC2018+
                        "ggTo4mu_Contin_MCFM701/ZZ4lAnalysis.root"),
            dict(name = "ggTo4e",filename = pathMC2018+
                        "ggTo4e_Contin_MCFM701/ZZ4lAnalysis.root"),
            dict(name = "ggTo4tau",filename = pathMC2018+
                        "ggTo4tau_Contin_MCFM701/ZZ4lAnalysis.root"),
            dict(name = "ggTo2e2mu",filename = pathMC2018+
                        "ggTo2e2mu_Contin_MCFM701/ZZ4lAnalysis.root"),       
            dict(name = "ggTo2e2tau",filename = pathMC2018+
                        "ggTo2e2tau_Contin_MCFM701/ZZ4lAnalysis.root"),
            dict(name = "ggTo2mu2tau",filename = pathMC2018+
                        "ggTo2mu2tau_Contin_MCFM701/ZZ4lAnalysis.root"),

            dict(name = "ZH125",filename = pathMC2018+
                        "ZH125/ZZ4lAnalysis.root"),

            dict(name = "VBFToZZTo4l",filename = pathMC2018 + 
                        "VBFToContinToZZTo4l/ZZ4lAnalysis.root"),
            dict(name = "TTZToLLNuNu",filename = pathMC2018 + 
                        "TTZToLLNuNu_M10ext1/ZZ4lAnalysis.root"),
            dict(name = "TTZJets",filename = pathMC2018 + 
                        "TTZJets_M10_MLMext1/ZZ4lAnalysis.root"),
        ]
    elif '2022EE' in outFile:
        samples = [
            dict(name = "ggH125",filename = pathMC2022EE+
                        "ggH125/ZZ4lAnalysis.root"),
            dict(name = "VBF125",filename = pathMC2022EE+
                        "VBFH125/ZZ4lAnalysis.root"),
            dict(name = "WplusH125",filename = pathMC2022EE+
                        "WplusH125/ZZ4lAnalysis.root"),
            dict(name = "WHminus125",filename = pathMC2022EE+
                        "WHminus125/ZZ4lAnalysis.root"),
            dict(name = "ttH125",filename = pathMC2022EE+
                        "ttH125/ZZ4lAnalysis.root"),

            dict(name = "ZZTo4l",filename = pathMC2022EE+
                        "ZZTo4l/ZZ4lAnalysis.root"),

            dict(name = "WWZ",filename = pathMC2022EE+
                        "WWZ/ZZ4lAnalysis.root"),
            dict(name = "ZZZ",filename = pathMC2022EE+
                        "ZZZ/ZZ4lAnalysis.root"),
            dict(name = "TTWW",filename = pathMC2022EE+
                        "TTWW/ZZ4lAnalysis.root"),
        ]
    elif '2022' in outFile:
        samples = [
            dict(name = "ggH125",filename = pathMC2022+
                        "ggH125/ZZ4lAnalysis.root"),
            dict(name = "VBF125",filename = pathMC2022+
                        "VBFH125/ZZ4lAnalysis.root"),
            dict(name = "WplusH125",filename = pathMC2022+
                        "WplusH125/ZZ4lAnalysis.root"),
            dict(name = "WHminus125",filename = pathMC2022+
                        "WHminus125/ZZ4lAnalysis.root"),
            dict(name = "ttH125",filename = pathMC2022+
                        "ttH125/ZZ4lAnalysis.root"),
            dict(name = "bbH125",filename = pathMC2022+
                        "bbH125/ZZ4lAnalysis.root"),

            dict(name = "ZZTo4l",filename = pathMC2022+
                        "ZZTo4l/ZZ4lAnalysis.root"),

            dict(name = "WWZ",filename = pathMC2022+
                        "WWZ/ZZ4lAnalysis.root"),
            dict(name = "WZZ",filename = pathMC2022+
                        "WZZ/ZZ4lAnalysis.root"),
            dict(name = "ZZZ",filename = pathMC2022+
                        "ZZZ/ZZ4lAnalysis.root"),
            dict(name = "TTWW",filename = pathMC2022+
                        "TTWW/ZZ4lAnalysis.root"),
            dict(name = "TTZZ",filename = pathMC2022+
                        "TTZZ/ZZ4lAnalysis.root"),
        ]


    of = ROOT.TFile.Open(outFile,"recreate") 
    
    for s in samples:
         histos = fillHistos(s["name"], s["filename"])
         for h in histos:
             of.WriteObject(h,h.GetName())

    of.Close()

def runData(outFile):

    of = ROOT.TFile.Open(outFile,"recreate") 
                
    hs_data = fillHistos("Data", pathDATA+ "/ZZ4lAnalysis.root")
    for h in hs_data:
        h.SetBinErrorOption(ROOT.TH1.kPoisson)
        of.WriteObject(h,h.GetName())

    of.Close()

if __name__ == "__main__" :

    # runMC('H4l_MC2018.root')
    # print('Running 2022')
    # runMC('H4l_MC2022.root')
    # print('Running 2022EE')
    # runMC('H4l_MC2022EE.root')
    print('Running data')
    runData('H4l_Data.root')
