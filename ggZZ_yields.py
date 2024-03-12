#!/bin/env python3
### Example macro for filling standard histograms from H4l nanoAODs.
### Histograms are stored on a file and can then be plotted with

# specify as argument the root file where to save yields histograms 
# run with: python3 ggZZ_yields.py ggZZ_2018.root ggZZ_2022EE.root


import argparse
import math
from pathlib import Path
from tabulate import tabulate
from typing import Dict

import ROOT
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ZZAnalysis.NanoAnalysis.tools import getLeptons, get_genEventSumw

ROOT.PyConfig.IgnoreCommandLineOptions = True


pathMC2018 = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/231209_nano/MC2018/" # FIXME: Use 2018 MC for the time being
# pathDATA_CD = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/231209_nano/Data2022_CD/" # data for 2022 period
# pathDATA_EFG = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/231209_nano/Data2022_EFG/" # data for 2022EE period
# pathMC2022 = '/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/231214_nano/MC2022/'
# pathMC2022EE = '/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/231214_nano/MC2022EE/'
pathMC2022EE = '../prod/ggZZ_no0p5_fb_12Mar2024/'


ZmassValue = 91.1876

maxEntriesPerSample = 1e12 # Use only up to this number of events in each MC sample, for quick tests.



ROOT.TH1.SetDefaultSumw2()

####################################
def fillHistos(samplename: str, filename: str, lumi: float) -> Dict[str, ROOT.TH1F] :
    """
    Fill histograms for yields.

    Parameters
    ----------
    samplename : str
        The name of the sample to read.
    filename : str
        The name of the file containing the sample to open.
    lumi : float
        The integrated luminosity

    Returns
    -------
    Dict[str, ROOT.TH1F]
        The dictionary containing the final state as key, and the histogram as values.

    Raises
    ------
    FileNotFoundError
        If the input file is not in the file system.
    ValueError
        If the flavour of the Z bosons is not correct.
    """
    
    ## yields ggZZ
    fs_list = ['4mu', '4e', '2e2mu', '4l']
    h_yield = {fs: ROOT.TH1F(f'h_yield_{fs}_{samplename}', f'h_yield_{fs}_{samplename}', 1, 0.0, 1.0) for fs in fs_list}

    # check if input file exists
    if not Path(filename).is_file():
        raise FileNotFoundError(f'Could not find input file {filename}!')
    
    f = ROOT.TFile.Open(filename)

    # assigning branches
    event = f.Events
    event.SetBranchStatus("*", 0)
    event.SetBranchStatus("run", 1)
    event.SetBranchStatus("luminosityBlock", 1)
    event.SetBranchStatus("event", 1)
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
        genEventSumw = get_genEventSumw(f, maxEntriesPerSample)


    # loop over events
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
                weight = (lumi*1000.* event.overallEventWeight*theZZ.dataMCWeight/genEventSumw)

            h_yield['4l'].Fill(0.5,weight) #yield 4l

            # per final state
            Z1flav = theZZ.Z1flav
            Z2flav = theZZ.Z2flav
            if(Z1flav==-169 and Z2flav==-169):
                currentFinalState = '4mu'
            elif(Z1flav==-121 and Z2flav==-121):
                currentFinalState = '4e'
            elif((Z1flav==-169 and Z2flav==-121) or 
                 (Z1flav==-121 and Z2flav==-169)):
                currentFinalState = '2e2mu'
            else:
                raise ValueError(f'Error in event {event.run}:{event.luminosityBlock}:{event.event}: found Z1flav={Z1flav}, Z2flav={Z2flav}!')

            h_yield[currentFinalState].Fill(0.5,weight) #yield 4l
        
    f.Close()
    

    return h_yield


   

def runMC(outFile): 

    if '2018' in outFile:
        path=pathMC2018
        lumi=59.7 #fb-1
    elif '2022EE' in outFile:
        path=pathMC2022EE
        lumi=27.007 #fb-1
    else:
        print('error: specify the path and lumi')
        
    samples = [
        dict(name = "ggTo4e",     filename = path+"ggTo4e_Contin_MCFM701/ZZ4lAnalysis.root"),
        dict(name = "ggTo4mu",    filename = path+"ggTo4mu_Contin_MCFM701/ZZ4lAnalysis.root"),
        dict(name = "ggTo4tau",   filename = path+"ggTo4tau_Contin_MCFM701/ZZ4lAnalysis.root"),
        dict(name = "ggTo2e2mu",  filename = path+"ggTo2e2mu_Contin_MCFM701/ZZ4lAnalysis.root"),       
        dict(name = "ggTo2e2tau", filename = path+"ggTo2e2tau_Contin_MCFM701/ZZ4lAnalysis.root"),
        dict(name = "ggTo2mu2tau",filename = path+"ggTo2mu2tau_Contin_MCFM701/ZZ4lAnalysis.root"),
    ]

    of = ROOT.TFile.Open(outFile,"recreate") 
    
    for s in samples:
         histos = fillHistos(s["name"], s["filename"], lumi)
         for h in histos.values():
             of.WriteObject(h,h.GetName())
          
    of.Close()



def printYields(inFile: str):
    """
    Print the yields in a pretty way

    Parameters
    ----------
    inFile : str
        File containing the yields histos
    """

    fs_list = ['4mu', '4e', '2e2mu', '4l']
    name_list = ['ggTo4e', 'ggTo4mu', 'ggTo4tau', 'ggTo2e2mu', 'ggTo2e2tau', 'ggTo2mu2tau']

    in_file = ROOT.TFile.Open(inFile, 'READ')

    
    # get histos from input file
    # define dictionary of list of histos
    #   input_histos = {'4mu': ['h_yield_4mu_ggTo4e', 'h_yield_4mu_ggTo4mu', etc.],
    #                   '4e': ['h_yield_4e_ggTo4e', 'h_yield_4e_ggTo4mu', etc.],
    #                    etc.}
    input_histos = {fs: [in_file.Get(f'h_yield_{fs}_{name}') for name in name_list] for fs in fs_list}
    # print yields for check
    print({fs: [h.GetBinContent(1) for h in h_list] for fs, h_list in input_histos.items()})

    # clone the first histogram of the list
    #   output_histos = {'4mu': 'h_yield_4mu_ggTo4e'.Clone('ggZZ_4mu'),
    #                    '4e': 'h_yield_4e_ggTo4e'.Clone('ggZZ_4e'),
    #                     etc.}    
    output_histos = {fs: h_list[0].Clone('ggZZ_{fs}') for fs, h_list in input_histos.items()}
    # print yields for check
    print({fs: h.GetBinContent(1) for fs, h in output_histos.items()})

    # add the remaining histograms
    #   loop over fs and the corresponding list of histos
    #   loop over the histos in the list starting from the second
    #   add to the corresponding output
    for fs, h_list in input_histos.items():
        for h in h_list[1:]:
            output_histos[fs].Add(h, 1.0)

    # print yields for check
    debug_dict = {fs: h.GetBinContent(1) for fs, h in output_histos.items()}
    print(debug_dict)
    debug_left = sum(list(debug_dict.values())[:-1])  # first three fs (4mu, 4e, 2e2mu)
    debug_right = list(debug_dict.values())[-1]  # last fs (4l)
    assert round(debug_left, 3) == round(debug_right, 3), "Yields do not add up!"


    # print yields pretty
    table = []
    for fs, h in output_histos.items():
        table.append([fs, h.GetBinContent(1), '+/-', h.GetBinError(1)])
    table = tabulate(table, headers=['fs', 'yields', '', 'unc.'], tablefmt='pipe', floatfmt='.3f', numalign='right', stralign='left')
    print(table)
         


def main(args: argparse.Namespace) -> int:

    for file in args.input:

        if args.hists:
            print(f'Making histograms for {file}...')
            runMC(file)

        print(f'Printing yields from {file}...')
        printYields(file)

    return 0


if __name__ == "__main__" :

    import sys

    parser = argparse.ArgumentParser(description='Print the yields', epilog='Contact info: Alessandra Cappati <alessandra.cappati@cern.ch>')
    parser.add_argument('input', nargs='+', help='input files')
    parser.add_argument('--hists', action='store_true', help='Remake histograms')
    args = parser.parse_args()

    code = main(args)
    sys.exit(code)

        
