from __future__ import print_function
import math
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ZZAnalysis.NanoAnalysis.tools import getLeptons


pathMC = '/eos/user/a/acappati/run3/MC2022/'
pathDATA = '/eos/user/a/acappati/run3/Data2022'

# defaults
ZmassValue = 91.1876
ROOT.TH1.SetDefaultSumw2()



def fillHistos(samplename, filename) :

    # def histo
    h_ZMass = ROOT.TH1F("ZMass_"+samplename,"ZMass_"+samplename,60,60.,120.)
    h_ZMass1 = ROOT.TH1F("ZMass1_"+samplename,"ZMass1_"+samplename,60,60.,120.)

    # open file
    f = ROOT.TFile.Open(filename)

    event = f.Events
    event.SetBranchStatus("*", 0)
    event.SetBranchStatus("run", 1)
    event.SetBranchStatus("luminosityBlock", 1)
    event.SetBranchStatus("*Muon*", 1)
    event.SetBranchStatus("*Electron*", 1)
    event.SetBranchStatus("*ZLLCand*", 1)
    event.SetBranchStatus("*ZCand*", 1)
    event.SetBranchStatus("bestZIdx", 1)
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


    iEntry=0
    printEntries=max(5000,nEntries/10)
    while iEntry<nEntries and event.GetEntry(iEntry):
        iEntry+=1
        if iEntry%printEntries == 0 : print("Processing", iEntry)

        bestCandIdx = event.bestCandIdx
        bestZIdx = event.bestZIdx

        # Check that the event contains a selected candidate, and that
        # passes the required triggers (which is necessary for samples
        # processed with TRIGPASSTHROUGH=True)
        #if(bestCandIdx != -1 and event.HLT_passZZ4l): #FIXME: 4l sel?
        if(bestCandIdx != -1): #FIXME: 4l sel?
        #if(bestZIdx != -1): # sel Z tree? #FIXME: what sel to use?
            weight = 1.
            ZZs = Collection(event, 'ZZCand')
            theZZ = ZZs[bestCandIdx]   
            Zs = Collection(event, 'ZCand')
            theZ = Zs[bestZIdx]        
            if isMC : weight = (event.overallEventWeight*theZZ.dataMCWeight/genEventSumw)
            mZ=theZ.mass
            h_ZMass.Fill(mZ,weight)
            # test
            if (mZ>60. and mZ<120.):
                h_ZMass1.Fill(mZ,weight)
        
    f.Close()
    
    return h_ZMass, h_ZMass1


def runMC():
    outFile = "hist_MC.root" 

    samples = [
        dict(name = "DY",filename = pathMC+'DYJetsToLL/ZZ4lAnalysis.root'),
    ]


    of = ROOT.TFile.Open(outFile,"recreate") 
    
    for s in samples:
         histos = fillHistos(s["name"], s["filename"])
         print('histos ', histos)
         for h in histos:
             of.WriteObject(h,h.GetName()) # works only if there are at least 2 histos

    of.Close()

def runData():
    outFile = "hist_Data.root" 

    of = ROOT.TFile.Open(outFile,"recreate") 
                
    hs_data = fillHistos("Data", pathDATA+ "/ZZ4lAnalysis.root")
    for h in hs_data:
        h.SetBinErrorOption(ROOT.TH1.kPoisson)
        of.WriteObject(h,h.GetName())

    of.Close()

if __name__ == "__main__" :
    runMC()
#    runData()
