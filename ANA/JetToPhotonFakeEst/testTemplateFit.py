from ROOT import TH1, TH1F, TH2, TH2D, TFile, TCanvas, THStack, TLegend, TMath
from ROOT import gROOT, gStyle
from ROOT import *
import array, glob, ROOT
import io
output=open('FitYields.txt','w')


import optparse
usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)

parser.add_option("--pt",action="store",type="string",dest="pt",default='1520')
parser.add_option("--eta",action="store",type="string",dest="eta",default='EB')
parser.add_option("--SBB",action="store",type="string",dest="SBB",default='')

(options, args) = parser.parse_args()

pt = options.pt
eta = options.eta
SBB = options.SBB


ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTitleSize(0)

def getTemp1DHist(filename, histname):
    file = TFile.Open(filename)    
    temp1DHist = file.Get(histname)
    temp1DHist.SetDirectory(0)
    temp1DHist.SetFillColor(0)

    return temp1DHist

def getWeightedHist(filename, histName, leftCut=None, rightCut=None):

    if leftCut==None and rightCut==None:
        # template is 1D histogram
        sumHist = getTemp1DHist(filename, histName)
        return sumHist
    
def makeFit(varname, varMin, varMax, signalHist, bkgHist, dataHist, MCplot, pteta, plotName):
#def makeFit(varname, varMin, varMax, signalHist, bkgHist, dataHist, pteta, plotName):
    # PhoSCRChHadIso will be input into varname 
    #    print 'variname = ', varname
    # RooFit variables
    argVar = RooRealVar(varname, varname, varMin, varMax)
    argList = RooArgList()
    argList.add(argVar)
    argSet = RooArgSet()
    argSet.add(argVar)

    ############ create PDF files ################
    # signal
    signalDataHist = RooDataHist('signalDataHist', 'signal RooDataHist', argList, signalHist)
    signalPdf = RooHistPdf('signalPdf', varname+' of Signal', argSet, signalDataHist)
    print 'signalPdf = ', signalPdf
    # background
    bkgDataHist = RooDataHist('BackgroundDataHist', 'BackgroundDataHist', argList, bkgHist)
    backgroundPdf = RooHistPdf('backgroundPdf', varname+' of background', argSet, bkgDataHist)
    # extra dy+Zg
    extrabkgDataHist = RooDataHist('ExtraBackgroundDataHist', 'ExtraBackgroundDataHist', argList, MCplot)
    ExtrabackgroundPdf = RooHistPdf('ExtrabackgroundPdf', varname+' of Extra background', argSet, extrabkgDataHist)
    # data
    dataDataHist = RooDataHist('data '+varname, varname+' in Data', argList, dataHist)

    # signal fraction parameter
    fractionSignal = RooRealVar('signal fraction', 'signal fraction', 0.1, 0.0, 1.0)
    #in terms of yields
    rooNTrue = RooRealVar("TrueYield","n True",nTrueStart,0,nMax)
    eTruePdf = RooExtendPdf("eTrue","extended True",signalPdf,rooNTrue)
    rooNFake = RooRealVar("FakeYield","n Fake",nFakeStart,0,nMax)
    eFakePdf = RooExtendPdf("eFake","extended Fake",backgroundPdf,rooNFake)
    rooNExtraFake = RooRealVar("ExtraFakeYield","n Extra Fake",rooNExtraFakeYield)#,0,nMax)
    eExtraFakePdf = RooExtendPdf("eExtraFake","extended extra Fake",ExtrabackgroundPdf,rooNExtraFake)
   # sumPDF =RooAddPdf("sumPdf","sig and bkg pdf",RooArgList(eTruePdf,eFakePdf))
    sumPDF =RooAddPdf("sumPdf","sig and bkg pdf",RooArgList(eTruePdf,eFakePdf,eExtraFakePdf))
    
    # total pdf of signal and background
    # model(x) = fsig*sig(x) + (1-fsig)*bkg(x)
 #   sumPDF = RooAddPdf('sumPdf', 'sig and bkg pdf', signalPdf, backgroundPdf, fractionSignal)

    ########## Now do the fit and plot the result #########
    # fit
    sumPDF.fitTo(dataDataHist, RooFit.SumW2Error(kFALSE))#, RooFit.PrintLevel(-1))

    # plot PDF and toy data overlaid
    if plotName!='':
        leg = TLegend(0.30, 0.80, 0.50, 0.95)
        leg.SetFillColor(0)
        leg.SetBorderSize(0)
        c1 = TCanvas('c1', 'c1', 700, 600)
        plotter = RooPlot(argVar, varMin, varMax, 20) # nBins = 20 is dummy
        dataDataHist.plotOn(plotter, RooFit.Name('Data'))
        sumPDF.plotOn(plotter, RooFit.Name('sumPdf'), RooFit.LineColor(kOrange+10))
        sumPDF.plotOn(plotter, RooFit.Components('signalPdf'), RooFit.Name('signal'), RooFit.LineColor(kGreen-4))
        sumPDF.plotOn(plotter, RooFit.Components('backgroundPdf'), RooFit.Name('background'), RooFit.LineColor(kAzure))
        sumPDF.plotOn(plotter, RooFit.Components('ExtrabackgroundPdf'), RooFit.Name('extra background'), RooFit.LineColor(kMagenta))
        #plotter.Draw()
        sumPDF.paramOn(plotter)
        plotter.Draw()
        MCplot.SetLineColor(kCyan)
        MCplot.SetFillColor(kCyan)
        MCplot.SetFillStyle(3001)
        #MCplot.Draw('samehist,e1')
        plotter.Draw("same")
        chi2 = RooChi2Var("chi2","chi2",sumPDF,dataDataHist)
       # print 'plotter->chiSquare: ', chi2.getVal(), pseudo_data_templ.GetNbinsX()-2, chi2.getVal()/(pseudo_data_templ.GetNbinsX()-2)
        print 'FakeYield', rooNFake.getVal(), '+-', rooNFake.getError()
        print 'TrueYield', rooNTrue.getVal(), '+-', rooNTrue.getError()
        output.write('FakeYield')

        print 'Photon_pT_'+pt+'_eta_'+eta
        print 'plotter->chiSquare: ', round(plotter.chiSquare('sumPdf','Data'),2)
        chi2 = 'chi2/ndf: '+str(round(plotter.chiSquare('sumPdf','Data'),2))
        DeRR = sqrt((rooNFake.getError()*rooNFake.getError())+(rooNTrue.getError()*rooNTrue.getError()))
        NeRR = rooNTrue.getError()
        eRR = (rooNTrue.getVal()/(rooNTrue.getVal()+rooNFake.getVal()))*sqrt((DeRR/(rooNTrue.getVal()+rooNFake.getVal()))**2+(NeRR/rooNTrue.getVal())**2)
        print '========================='
        print 'signal fraction', round(rooNTrue.getVal()/(rooNTrue.getVal()+rooNFake.getVal()),2), ' +- ', round(eRR,2)
        sigFrac = 'signal fraction: '+str(round(rooNTrue.getVal()/(rooNTrue.getVal()+rooNFake.getVal()),2))+'  +-  '+str(round(eRR,2))
        print '========================='
        leg.AddEntry(sumPdf, 'Fit', 'l')
        leg.AddEntry(signal, 'Signal', 'l')
        leg.AddEntry(background, 'Background', 'l')
        leg.Draw('same')
        box = TPaveText(0.60,0.70,.95,.80,"TRNDC")
        box.Draw(); 
        box.SetFillColor(0)
        box.SetBorderSize(0)
        #box.AddText(str(sigFrac));
        box1 = TPaveText(0.60,0.65,.95,.70,"TRNDC")
        box1.Draw(); 
        box1.SetFillColor(0)
        box1.SetBorderSize(0)
        #box1.AddText(str(chi2));
        box2 = TPaveText(0.60,0.55,.95,.60,"TRNDC")
        box2.Draw(); 
        box2.SetFillColor(0)
        box2.SetBorderSize(0)
        box2.AddText(str(pteta));
        c1.SaveAs(plotName)

    print 'fit returned value ',(rooNTrue.getVal()/(rooNTrue.getVal()+rooNFake.getVal())), ' +- ',eRR
    return rooNTrue.getVal(),rooNFake.getVal()

def optimizeBinBoundaries(templ, minAccumulate, firstBinValue = -9999, endBinValue = 9999):
    nbins = templ.GetNbinsX()
    print "-------------------------------"
    print "nbins", nbins
    accumulate = 0
    binlist = []

    # ignore empty bins in the beginning
    for firstIndex in xrange(1, nbins+1):
        #print "1", templ.GetBinLowEdge(firstIndex)
        if templ.GetBinLowEdge(firstIndex)>=firstBinValue:
            binlist.append(templ.GetBinLowEdge(firstIndex))
            break

    for Index in xrange(firstIndex, nbins+1):
        #print "2", templ.GetBinLowEdge(firstIndex)
        if templ.GetBinLowEdge(Index)>=endBinValue:
            binlist.append(endBinValue)
            break
        if accumulate > minAccumulate:
            accumulate=0
            #print "3", templ.GetBinLowEdge(firstIndex)
            binlist.append(templ.GetBinLowEdge(Index))
            #print "4", templ.GetBinContent(Index)
        accumulate+=templ.GetBinContent(Index)
        print "accumulate:", accumulate
    print "BinLowEdge:", templ.GetBinLowEdge(nbins+1)
    binlist.append(templ.GetBinLowEdge(nbins+1)) # bin end

    newarray = array.array('d') 
    newarray.fromlist(binlist)
    print "-------------------------------"
    return newarray

def makeUniformHist(hist):
    nbins = hist.GetNbinsX()
    print 'nbins = ', nbins

    outHist = TH1F(hist.GetName()+'_uni',hist.GetName()+'_uni',nbins,1,nbins+1)
    for index in xrange(1, nbins+1):
        outHist.SetBinContent(index, hist.GetBinContent(index))
        outHist.SetBinError(index, hist.GetBinError(index))
    return outHist

#################### MAIN CODE: Do The Fitting #######################
#InputFile0 = 'Temp_pseudodata_W.root'
InputData = 'Hist_miniTreewSFsNoSieieIsoMegPSVwMt_mrgd_data_vetoEIDonly.root'
MCDY = 'Hist_miniTreewSFsNoSieieIsoMegPSVwMt_job_summer12_DYJetsToLL_skimmed_wg_vetoEIDonly.root'
MCZg = 'Hist_miniTreewSFsNoSieieIsoMegPSVwMt_job_summer12_Zg_skimmed_wg_vetoEIDonly.root'

#inFake =     'ZjetsISOTemp_mrgd_BkgISOTemp_job_2electron_data.root'
#inReal =     'ZgFSRISOTemp_mrgd_SigISOTemp_job_2electron_data.root'
inFake = 'CorrZjetsISOTemp_mrgd_BkgISOTemp_job_2electron_data.root'
inReal = 'CorrZgFSRISOTemp_mrgd_SigISOTemp_job_2electron_data.root'

#for photon pt plot
#Sig_temp = "h_Pho_ISO_"+eta
#Bkg_temp = "h_ISO_"+eta+"_"+pt
#All_temp = "h_ChHadIso_"+eta+"_"+pt
#fine_data_temp = "hfine_ChHadIso_"+eta+"_"+pt
#fine_real_temp = "hfine_Pho_ISO_"+eta

#for others
#Sig_temp = "h_Pho_ISO_"+eta
#Bkg_temp = "h_ISO_"+eta
#All_temp = "h_ChHadIso"+eta+"_Zveto_MegGT50"
#fine_data_temp = "hfine_ChHadIso"+eta+"_Zveto_MegGT50"
#fine_real_temp = "hfine_Pho_ISO_"+eta

#for W_MT plot, but fake factors from each pt(g) bin
Sig_temp = "hfine_Pho_ISO_"+eta
Bkg_temp = "hfine_ISO_"+eta+"_"+pt
if (eta == 'EE'):
    Bkg_temp = "hfine_ISO_"+eta
All_temp =       "hfine_ChHadIso_"+eta+"_Zveto_MegGT50_"+pt#+"_WMTGT40"
fine_data_temp = "hfine_ChHadIso_"+eta+"_Zveto_MegGT50_"+pt#+"_WMTGT40"
fine_real_temp = "hfine_Pho_ISO_"+eta
plotMC= "hfine_ChHadIso_"+eta+"_Zveto_MegGT50_"+pt#+"_WMTGT40"

#Sig_temp = "hfine_Pho_ISO_"+eta
#Bkg_temp = "hfine_ISO_"+eta

#if (eta == 'EE' and (pt=='4') ):
 #   Bkg_temp = "h_ISO_"+eta+"_3"
#if (pt=='4'):
#    Bkg_temp = "h_ISO_"+eta+"_3"
#if (pt=='6' or pt=='7' or pt=='8' or pt=='9' or pt=='10' or pt=='11' or pt=='12'):
 #   Bkg_temp = "h_ISO_"+eta+"_6"

print Sig_temp, Bkg_temp

#try with nominal stuff
#All_temp = "h_ChHadIso_"+eta+"_"+pt
#Sig_temp = "h_RandConeChHadIso_"+eta+"_"+pt
#Bkg_temp = "h_ChHadIso_"+eta+"_SB"+SBB+"_"+pt

FitVariable = 'ChargedHadIsolationSB'+SBB
print Sig_temp, Bkg_temp, All_temp

# Fit range
LowFitRange = -1
HighFitRange = 20

# signal template
sig_templ = getWeightedHist(inReal, Sig_temp)
# background template
bkg_templ = getWeightedHist(inFake, Bkg_temp)
# pseudo data
pseudo_data_templ = getWeightedHist(InputData, All_temp)
fine_real_templ = getWeightedHist(inReal, fine_real_temp)
#signal MC_truth (gen-matched)
#MC_truth_plot= getWeightedHist(InputFile, MCtruth_plot)
#bkg MC_truth (gen-matched)
#MC_truth_plot= getWeightedHist(InputFile1, MCtruth_plot)
hMCplotDY= getWeightedHist(MCDY, plotMC)
hMCplotZg= getWeightedHist(MCZg, plotMC)
# ##########
# pseudo_sig = getWeightedHist(InputFile, trueSig_hist2d_temp, 0.0, SihihSelCut)
# Lbin = pseudo_data_templ.FindBin(LowFitRange)
# Rbin = pseudo_data_templ.FindBin(HighFitRange)
# MCtrueSelSignalFraction = pseudo_sig.Integral(Lbin, Rbin)/pseudo_data_templ.Integral(Lbin, Rbin)
# ##########
print "data",pseudo_data_templ.Integral()
#pseudo_data_templ.Add(hMCplot,-1)
#print "data",pseudo_data_templ.Integral()

errS = ROOT.Double()
errB = ROOT.Double()
errD = ROOT.Double()
print 'Signal Integral = ',sig_templ.IntegralAndError(-1,43,errS),errS
print 'Background Integral = ',bkg_templ.IntegralAndError(-1,43,errB),errB
print 'Data Integral = ',pseudo_data_templ.IntegralAndError(-1,43,errD),errD
#print 'Pseudo Data Integral = ',MC_truth_plot.Integral()

nMax = pseudo_data_templ.Integral()
if(eta == 'EB' ):
    nTrueStart = 0.5*nMax
    nFakeStart = 0.5*nMax
if(eta == 'EE' ):
    nTrueStart = 0.5*nMax
    nFakeStart = 0.5*nMax

DYscalePtbinned_EB = [1.561,1.556,1.475,1.463]
DYscalePtbinned_EE = [1.109,1.094,1.024,0.994]
#DYscalePtbinned_EB = [1.561,1.556,1.475,1.467,1.386,1.603,1.870,1.530,1.194,1.274,1.894,1.679]
#DYscalePtbinned_EE = [1.109,1.094,1.024,0.903,0.948,1.214,1.313,1.037,0.790,1.104,0.901,0.898]
rooNExtraFakeYield = 1#MCintegralDY*1.215
MCintegralZg = hMCplotZg.Integral()
MCintegralDY = hMCplotDY.Integral()
print MCintegralZg, MCintegralDY, MCintegralZg+MCintegralDY

if(eta == 'EB' ):
    hMCplotDY.Scale(DYscalePtbinned_EB[3])
if(eta == 'EE' ):
    hMCplotDY.Scale(DYscalePtbinned_EE[0])

hMCplotDY.Add(hMCplotZg)
MCintegral = hMCplotDY.Integral()
print MCintegralZg, MCintegralDY, MCintegralZg+MCintegralDY, MCintegral
rooNExtraFakeYield = MCintegral
#pseudo_data_templ.Add(hMCplotDY,-1)

print "DY+Zg",MCintegralDY, DYscalePtbinned_EB[0], rooNExtraFakeYield
#hMCplot = hMCplotDY.Add(hMCplotZg);
# Rebin the histograms to make equal bins
#if(eta == 'EB' ):
#    sig_templ.Rebin(2)
#    bkg_templ.Rebin(2)
#    pseudo_data_templ.Rebin(2)
#if(eta == 'EE' ):
#    sig_templ.Rebin(4)
#    bkg_templ.Rebin(4)
#    pseudo_data_templ.Rebin(4)
#MC_truth_plot.Rebin(4)

boundaries = optimizeBinBoundaries(bkg_templ, 0., LowFitRange, HighFitRange)

#binlist = [-1.0,0.0,2.0,5.0, 8.,13.]
#binlist = [-1.0, 1.5,3., 5., 7.,10., 14.0,21.0][-1.0, 0.0, 1., 2.0,5.0, 10.0, 15.0, 20]
if(eta == 'EB' ):
    binlist = [-1,0,1.5,2.5,3.5,5.0,7.0,10.0,14.0,20.0]
   # binlist = [-1,0.,1.,2,3,4,5.0,6,7.0,8,9,10.0,11,12,13,14.0,15,16,17,18,19,20.0]
if(eta == 'EE' ):
    binlist = [-1.0,0.,1.2,2.8,5,7,10,14.,20.0]
print 'using custom bins for Data:',binlist
#boundaries = array.array('d')
#boundaries.fromlist(binlist)

# need to set the bin boundaries manually for fitting data. 
#(read for Rebin: https://root.cern.ch/doc/master/classTH1.html#aff6520fdae026334bf34fa1800946790)
# Rebin (number of variable size bins, 
print 'new number of boundaries = ',len(boundaries)

pseudo_data_templRbin = pseudo_data_templ.Clone()
bkg_templRbin = bkg_templ.Clone()
sig_templRbin = sig_templ.Clone()
hMCplotDYRbin = hMCplotDY.Clone()
#pseudo_data_templRbin = pseudo_data_templ.Rebin(len(boundaries)-1, pseudo_data_templ.GetName()+'rebin', boundaries)
#bkg_templRbin = bkg_templ.Rebin(len(boundaries)-1, bkg_templ.GetName()+'rebin', boundaries)
#sig_templRbin = sig_templ.Rebin(len(boundaries)-1, sig_templ.GetName()+'rebin', boundaries)
#hMCplotDYRbin = hMCplotDY.Rebin(len(boundaries)-1, hMCplotDY.GetName()+'rebin', boundaries)

# to avoid signal going below fit range in the plot
print 'sig_templRbin underflow = ', sig_templRbin.GetBinContent(0)

sig_templRbin.SetBinContent(0, 0)

pseudo_data_templUni = makeUniformHist(pseudo_data_templRbin)
bkg_templUni = makeUniformHist(bkg_templRbin)
sig_templUni = makeUniformHist(sig_templRbin)

# using LowFitRange and HighFitRange to fo the fit doesn't work well
# to get the good fit result, do fit for fixed size histograms
LowUniFitRange = pseudo_data_templRbin.FindBin(LowFitRange)
if LowUniFitRange == 0:
#    print 'lower bin is underflow, setting to 1'
    LowUniFitRange = 1

HighUniFitRange = pseudo_data_templRbin.FindBin(HighFitRange) + 1
if pseudo_data_templRbin.GetBinLowEdge(HighUniFitRange-1) == HighFitRange:
#    print 'upper bin in on the border of fit range, reducing'
    HighUniFitRange -= 1

print 'fitting in the bin range = ',LowUniFitRange,': ', HighUniFitRange


# Fit results: overlap signal, bkg, sum of sig and bkg, and pseudo data
(nsig,nbkg) = makeFit(FitVariable, LowFitRange, HighFitRange, sig_templRbin, bkg_templRbin, pseudo_data_templRbin, hMCplotDYRbin, eta+', '+pt+'-bin', 'fit_'+FitVariable+'_'+eta+'_'+pt+'.png')


#############################
if(eta == 'EB' ):
    leftbin = -1#LowFitRange
    rightbin = 43#HighFitRange
if(eta == 'EE' ):
    leftbin = -1#LowFitRange
    rightbin = 106#HighFitRange

c2 = TCanvas('c2', 'c2', 700, 600)

leg = TLegend(0.55, 0.65, 0.85, 0.85)
leg.SetFillColor(0)
leg.SetBorderSize(0)

print 'bins for normalized = ', leftbin, ' : ', rightbin
# signal template
sig_DDScale = getWeightedHist(inReal, Sig_temp)
# background template
bkg_DDScale = getWeightedHist(inFake, Bkg_temp)

sig_DDScale.Scale(nsig/sig_DDScale.Integral(leftbin,rightbin))
#fitSigFrac*pseudo_data_templRbin.Integral(leftbin,rightbin)
bkg_DDScale.Scale(nbkg/bkg_DDScale.Integral(leftbin,rightbin))

fine_data_templ = getWeightedHist(InputData, fine_data_temp)

err = ROOT.Double()
errBkg = ROOT.Double()
if(eta == 'EB' ):
    print "------------"
#    print "data in sig region: ",fine_data_templ.Integral(-1,5)
#    print "effi of real temp: ",fine_real_templ.Integral(-1,5)/fine_real_templ.Integral(-1,43)
#    print "real yield in sig region: ",nsig*fine_real_templ.Integral(-1,5)/fine_real_templ.Integral(-1,43)
#    print "fake in sig region: ", fine_data_templ.Integral(-1,5) - (nsig*fine_real_templ.Integral(-1,5)/fine_real_templ.Integral(-1,43))
    print 'Signal & Err  = ',sig_DDScale.IntegralAndError(-1,5,err),err
    print 'Background  = ',bkg_DDScale.IntegralAndError(-1,5,errBkg),errBkg
    print 'Data  = ',pseudo_data_templRbin.Integral(-1,5)
    print "------------"
if(eta == 'EE' ):
    print "------------"
#    print "data in sig region: ",fine_data_templ.Integral(-1,11)
#    print "effi of real temp: ",fine_real_templ.Integral(-1,11)/fine_real_templ.Integral(-1,107)
#    print "real yield in sig region: ",nsig*fine_real_templ.Integral(-1,11)/fine_real_templ.Integral(-1,107)
#    print "fake in sig region: ", fine_data_templ.Integral(-1,11) - (nsig*fine_real_templ.Integral(-1,11)/fine_real_templ.Integral(-1,107))
    print 'Signal & Err  = ',sig_DDScale.IntegralAndError(-1,11,err),err
    print 'Background  = ',bkg_DDScale.IntegralAndError(-1,11,errBkg),errBkg
    print 'Data  = ',pseudo_data_templRbin.Integral(-1,11)
    print "------------"

#print 'Pseudo Data Integral = ',MC_truth_plot.Integral()


sig_DDScale.SetLineColor(kGreen-4)
#sig_DDScale.SetFillColor(kGreen-4)
sig_DDScale.GetXaxis().SetRangeUser(LowFitRange, HighFitRange)
sig_DDScale.Sumw2()
#sig_DDScale.Draw()
#c2.SaveAs("sig_DDScale.png")

bkg_DDScale.SetLineColor(kAzure)
#bkg_DDScale.SetFillColor(kAzure)
bkg_DDScale.GetXaxis().SetRangeUser(LowFitRange, HighFitRange)
bkg_DDScale.Sumw2()
bkg_DDScale.Draw("text")
c2.SaveAs("bkg_DDScale.png")

# stack = THStack(FitVariable+'_stack', FitVariable+'_stack')
# stack.Add(bkg_DDScale)
# stack.Add(sig_DDScale)
# stack.Draw('hist')
# #stack.SetMaximum(120)
# stack.GetXaxis().SetTitle('photon '+FitVariable + ' [GeV]')

pseudo_data_templRbin.SetMarkerColor(1)
pseudo_data_templRbin.SetMarkerStyle(8)
pseudo_data_templRbin.Sumw2()
#pseudo_data_templRbin.Draw()
#c2.SaveAs("pseudo_data_templRbin.png")


#stack = THStack(FitVariable+'_stack', FitVariable+'_stack')
#stack.Add(bkg_templRbin)
#stack.Add(sig_templRbin)
#stack.Draw('hist')
#stack.SetMaximum(120)
#stack.GetXaxis().SetTitle('photon '+FitVariable + ' [GeV]')

leg.AddEntry(pseudo_data_templRbin, 'Pseudo Data', 'lp')
leg.AddEntry(sig_DDScale, 'Signal', 'f')
leg.AddEntry(bkg_DDScale, 'Background', 'f')
#leg.Draw('same')
#pseudo_data_templRbin.Draw("sames")



output.close()


