from ROOT import TH1, TH1F, TH2, TH2D, TFile, TCanvas, THStack, TLegend, TMath, TList
from ROOT import gROOT, gStyle
from ROOT import *
import array, glob, ROOT
import sys, os, re, uuid, math, copy, imp, random, collections, pickle, time
from array import array
from uncertainties import ufloat
from uncertainties import unumpy

ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)

def format_hist( hist, color ) :
    hist.SetLineColor( color )
    hist.SetMarkerColor( color )
    hist.SetLineWidth(2)
    hist.SetMarkerSize(0)
    hist.SetStats(0)

def get_default_cuts(var) :
    
    if var == 'sigmaIetaIeta' :
        return { 'EB' : { 'tight' : ( 0, 0.01099 ), 'loose' : ( 0.01101, 0.0299 ) },
                 'EE' : { 'tight' : ( 0, 0.03299 ), 'loose' : ( 0.033, 0.099 ) } 
               }

def get1PhoHist(filename, histname) :
    file = TFile.Open(filename)
    phoHist = file.Get(histname)
    #   phoHist.SetDirectory(0)
    #    phoHist.SetFillColor(0)

    return phoHist

def getTemp1DHist(filename, histname) :
    file = TFile.Open(filename)
    temp1DHist = file.Get(histname)
    #    temp1DHist.SetDirectory(0)
    #    temp1DHist.SetFillColor(0)

    return temp1DHist

def getWeightedHist(filename, histName):

    # template is 1D histogram
    sumHist = getTemp1DHist(filename, histName)
    return sumHist

# fit on a single photon
def run_photon_fit(templates, gg_hist, lead_eta, var, outputDir=None) :

    accept_reg = ['EB', 'EE']
    if lead_eta not in accept_reg :
        print 'Lead region does not make sense'
        return

    # cut point for sigmaIetaIeta
    sigmaIetaIetaCut = get_default_cuts(var)

    bins_lead_tight = (gg_hist.GetXaxis().FindBin(sigmaIetaIetaCut[lead_eta]['tight'][0]), gg_hist.GetXaxis().FindBin(sigmaIetaIetaCut[lead_eta]['tight'][1]))
    bins_lead_loose = (gg_hist.GetXaxis().FindBin(sigmaIetaIetaCut[lead_eta]['loose'][0]), gg_hist.GetXaxis().FindBin(sigmaIetaIetaCut[lead_eta]['loose'][1]))

    print 'cut, bins lead, tight = %f-%f (%d - %d) ' %(sigmaIetaIetaCut[lead_eta]['tight'][0], sigmaIetaIetaCut[lead_eta]['tight'][1], bins_lead_tight[0], bins_lead_tight[1])
    print 'cut, bins lead, loose = %f-%f (%d - %d) ' %(sigmaIetaIetaCut[lead_eta]['loose'][0], sigmaIetaIetaCut[lead_eta]['loose'][1], bins_lead_loose[0], bins_lead_loose[1])

    Ndata_T = gg_hist.Integral(bins_lead_tight[0], bins_lead_tight[1])
    Ndata_L = gg_hist.Integral(bins_lead_loose[0], bins_lead_loose[1])

    Ndata = {}
    Ndata['T'] = ufloat( Ndata_T, math.sqrt(Ndata_T ), 'Ndata_T' )
    Ndata['L'] = ufloat( Ndata_L, math.sqrt(Ndata_L ), 'Ndata_L' )

    print 'N data T = ', Ndata['T']
    print 'N data L = ', Ndata['L']

    eff_cuts = {}
    eff_cuts['lead'] = {}
    eff_cuts['lead']['tight'] = sigmaIetaIetaCut[lead_eta]['tight']
    eff_cuts['lead']['loose'] = sigmaIetaIetaCut[lead_eta]['loose']

    (eff_1d_stat, eff_1d_syst) = generate_1d_efficiencies(templates, eff_cuts, lead_eta)
    print 'EFF_1D'
    print eff_1d_stat

    results_stat = run_fit({'T': Ndata['T'], 'L' : Ndata['L']}, eff_1d_stat)
    results_syst = run_fit( {'T': Ndata['T'], 'L' : Ndata['L']}, eff_1d_syst )

    print '**********************************'
    print 'Real Norm = ', results_stat.item(0)
    print 'Fake Norm = ', results_stat.item(1)
    print '**********************************'

    hist_temp_lead_real = templates['lead']['real'].Clone( 'hist_temp_lead_real' )
    hist_temp_lead_fake = templates['lead']['real'].Clone( 'hist_temp_lead_fake' )

    p_R_T = results_stat.item(0)*eff_1d_stat['eff_R_T_lead']
    p_R_L = results_stat.item(0)*eff_1d_stat['eff_R_L_lead']
    p_F_T = results_stat.item(1)*eff_1d_stat['eff_F_T_lead']
    p_F_L = results_stat.item(1)*eff_1d_stat['eff_F_L_lead']

    print '====================================='
    print 'nPred Real Tight = ', p_R_T
    print 'nPred Real Loose = ', p_R_L
    print 'nPred Fake Tight = ', p_F_T
    print 'nPred Fake Loose = ', p_F_L
    print '====================================='


    text_results_stat = collect_results( results_stat, Ndata, eff_1d_stat, templates, bins_lead_loose, bins_lead_tight)
    text_results_syst = collect_results( results_syst, Ndata, eff_1d_stat, templates, bins_lead_loose, bins_lead_tight)
    return text_results_stat, text_results_syst


def collect_results(results, data, efficiencies, templates, bins_lead_loose, bins_lead_tight) :

    text_results = collections.OrderedDict()

    for key, val in efficiencies.iteritems() :
        text_results[key] = val

    text_results['Ndata_T'] = data['T']
    text_results['Ndata_L'] = data['L']

    text_results['alpha_R'] = results.item(0)
    text_results['alpha_F'] = results.item(1)

    text_results['Npred_R_T'] = text_results['alpha_R']*text_results['eff_R_T_lead']
    text_results['Npred_R_L'] = text_results['alpha_R']*text_results['eff_R_L_lead']

    text_results['Npred_F_T'] = text_results['alpha_F']*text_results['eff_F_T_lead']
    text_results['Npred_F_L'] = text_results['alpha_F']*text_results['eff_F_L_lead']

    int_lead_real_loose = get_integral_and_error(templates['lead']['real'], bins_lead_loose )
    int_lead_real_tight = get_integral_and_error(templates['lead']['real'], bins_lead_tight )
    int_lead_fake_loose = get_integral_and_error(templates['lead']['fake'], bins_lead_loose )
    int_lead_fake_tight = get_integral_and_error(templates['lead']['fake'], bins_lead_tight )

    text_results['template_int_lead_real_loose'] = int_lead_real_loose
    text_results['template_int_lead_real_tight'] = int_lead_real_tight
    text_results['template_int_lead_fake_loose'] = int_lead_fake_loose
    text_results['template_int_lead_fake_tight'] = int_lead_fake_tight

    return text_results

def generate_1d_efficiencies( templates, eff_cuts, lead_eta) :

    (int_stat, int_syst) = get_template_integrals( templates, eff_cuts, lead_eta)

    print 'Template integral'
    print int_stat

    (eff_1d_stat, eff_1d_syst) = get_1d_loose_efficiencies( int_stat, int_syst, lead_eta)

    return eff_1d_stat, eff_1d_syst

def get_1d_loose_efficiencies( int_stat, int_syst, lead_eta) :
    eff_stat = {}
    eff_syst = {}

    if int_stat['lead']['real']['loose'].n == 0 :
        eff_stat['eff_R_T_lead'] = ufloat( 1.0, int_stat['lead']['real']['tight'].s/int_stat['lead']['real']['tight'].n )
    else :
        eff_stat['eff_R_T_lead'] = int_stat['lead']['real']['tight'] / (int_stat['lead']['real']['tight']+int_stat['lead']['real']['loose'])
    eff_stat['eff_F_T_lead'] = int_stat['lead']['fake']['tight'] / (int_stat['lead']['fake']['tight']+int_stat['lead']['fake']['loose'])
    eff_stat['eff_R_L_lead'] = int_stat['lead']['real']['loose'] / (int_stat['lead']['real']['tight']+int_stat['lead']['real']['loose'])
    eff_stat['eff_F_L_lead'] = int_stat['lead']['fake']['loose'] / (int_stat['lead']['fake']['tight']+int_stat['lead']['fake']['loose'])

    eff_syst['eff_R_T_lead'] = int_syst['lead']['real']['tight'] / (int_syst['lead']['real']['tight']+int_syst['lead']['real']['loose'])
    eff_syst['eff_F_T_lead'] = int_syst['lead']['fake']['tight'] / (int_syst['lead']['fake']['tight']+int_syst['lead']['fake']['loose'])
    eff_syst['eff_R_L_lead'] = int_syst['lead']['real']['loose'] / (int_syst['lead']['real']['tight']+int_syst['lead']['real']['loose'])
    eff_syst['eff_F_L_lead'] = int_syst['lead']['fake']['loose'] / (int_syst['lead']['fake']['tight']+int_syst['lead']['fake']['loose'])

    # Do systematics
    # the integrals may already have some systematics
    # that are propagated to the eff_*
    # therefore, don't overwrite the existing 
    # systematics, but make a ufloat with a
    # zero value, and non-zero syst
    eff_syst['eff_R_L_lead'] = eff_syst['eff_R_L_lead'] + ufloat( 0.0, math.fabs(eff_syst['eff_R_L_lead'].n)*0.0, 'Template_lead_real_loose')
    eff_syst['eff_F_L_lead'] = eff_syst['eff_F_L_lead'] + ufloat( 0.0, math.fabs(eff_syst['eff_F_L_lead'].n)*0.0, 'Template_lead_fake_loose' )
    eff_syst['eff_R_T_lead'] = eff_syst['eff_R_T_lead'] + ufloat( 0.0, math.fabs(eff_syst['eff_R_T_lead'].n)*0.0, 'Template_lead_real_loose')
    eff_syst['eff_F_T_lead'] = eff_syst['eff_F_T_lead'] + ufloat( 0.0, math.fabs(eff_syst['eff_F_T_lead'].n)*0.0, 'Template_lead_fake_loose' )

    return eff_stat, eff_syst

def get_integral_and_error( hist, bins=None, name='' ) :

    err = ROOT.Double()
    if bins is None :
        val = hist.IntegralAndError( 1, hist.GetNbinsX(), err )
    else :
        if bins[1] is None :
            val = hist.IntegralAndError( bins[0], hist.GetNbinsX(), err )
        else :
            val = hist.IntegralAndError( bins[0], bins[1], err )

    return ufloat( val, err, name )

def get_template_integrals( templates, eff_cuts, lead_eta) :

    int_stat = {}
    int_stat['lead']={}
    int_stat['lead']['real']={}
    int_stat['lead']['fake']={}

    int_syst = {}
    int_syst['lead']={}
    int_syst['lead']['real']={}
    int_syst['lead']['fake']={}

    bins_lead_real_tight = ( templates['lead']['real'].GetXaxis().FindBin( eff_cuts['lead']['tight'][0] ), templates['lead']['real'].GetXaxis().FindBin( eff_cuts['lead']['tight'][1] ) )
    bins_lead_real_loose = ( templates['lead']['real'].GetXaxis().FindBin( eff_cuts['lead']['loose'][0] ), templates['lead']['real'].GetXaxis().FindBin( eff_cuts['lead']['loose'][1] ) )
    bins_lead_fake_tight = ( templates['lead']['fake'].GetXaxis().FindBin( eff_cuts['lead']['tight'][0] ), templates['lead']['fake'].GetXaxis().FindBin( eff_cuts['lead']['tight'][1] ) )
    bins_lead_fake_loose = ( templates['lead']['fake'].GetXaxis().FindBin( eff_cuts['lead']['loose'][0] ), templates['lead']['fake'].GetXaxis().FindBin( eff_cuts['lead']['loose'][1] ) )

    int_stat['lead']['real']['tight'] = get_integral_and_error( templates['lead']['real'], bins_lead_real_tight, 'Data_lead_real_tight' )
    int_stat['lead']['real']['loose'] = get_integral_and_error( templates['lead']['real'], bins_lead_real_loose, 'Data_lead_real_loose' )
    int_stat['lead']['fake']['tight'] = get_integral_and_error( templates['lead']['fake'], bins_lead_fake_tight, 'Data_lead_fake_tight' )
    int_stat['lead']['fake']['loose'] = get_integral_and_error( templates['lead']['fake'], bins_lead_fake_loose, 'Data_lead_fake_loose' )

    # If running with systematics, set the data systs to zero
    # May need to implement non-zero systematics for data in the future
    # The overall template systematics should not be set here
    int_syst['lead']['real']['tight'] = ufloat(int_stat['lead']['real']['tight'].n, 0.0 , 'Data_lead_real_tight' )
    int_syst['lead']['real']['loose'] = ufloat(int_stat['lead']['real']['loose'].n, 0.0 , 'Data_lead_real_loose' )
    int_syst['lead']['fake']['tight'] = ufloat(int_stat['lead']['fake']['tight'].n, 0.0 , 'Data_lead_fake_tight' )
    int_syst['lead']['fake']['loose'] = ufloat(int_stat['lead']['fake']['loose'].n, 0.0 , 'Data_lead_fake_loose' )

    return int_stat, int_syst

def run_fit( data, efficiencies ) :

    # make the matrix
    matrix = generate_eff_matrix( efficiencies )
    print matrix

    #do the fit!  Invert the matrix and multiply the by counts vectors
    results = solve_matrix_eq( matrix, [data['T'], data['L']] )

    return results 

def generate_eff_matrix( eff_dic) :

    eff_matrix = [ [ eff_dic['eff_R_T_lead'], eff_dic['eff_F_T_lead'] ],
                   [ eff_dic['eff_R_L_lead'], eff_dic['eff_F_L_lead'] ] ]
    
    return eff_matrix

def solve_matrix_eq( matrix_ntries, vector_entries ) :

    ms = []
    mn = []
    for row in matrix_ntries :
        ms_row = []
        mn_row = []
        for col in row :
            ms_row.append( col.s )
            mn_row.append( col.n )
        ms.append( ms_row )
        mn.append( mn_row )

    matrix = unumpy.umatrix( mn, ms )

    print matrix

    vs = []
    vn = []
    for row in vector_entries :
        vn.append( [ row.n ] )
        vs.append( [ row.s ] )

    vector = unumpy.umatrix( vn, vs )
    
    inv_matrix = None
    try :
        inv_matrix = matrix.getI()
    except :
        print 'Failed to invert matrix, aborting'
        raise
        return None

    print inv_matrix
    print vector

    return inv_matrix*vector

def run_fit_manual(data, eff) :

    alpha_rf = ( (1.0/( (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl'] 
                      + (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                      + (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl']
                      - (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                      - (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl']
                      - (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl'] ) )
              * ( (  eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl']
                   - eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl'] )*data['TL']
                + (  (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*eff['eff_R_L_subl']
                   - (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_F_L_lead']*eff['eff_F_L_subl'] )*data['LT']
                + (   (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])
                    - (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl']) )*data['LL']
              ) )

    alpha_fr = ( (1.0/( (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl']     
                    + (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                    + (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl']
                    - (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                    - (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl']
                    - (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl'] ) )
            * ( (   eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                  - eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl'] )*data['TL'] 
              + (   (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*eff['eff_F_L_subl']
                  - (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_R_L_lead']*eff['eff_F_L_subl'] )*data['LT']
              + (   (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])
                  - (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl']) ) *data['LL']
              ) )
                
    alpha_ff = ( (1.0/( (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl']    
                    + (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                    + (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl']
                    - (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                    - (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl']
                    - (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl'] ) )
            * ( (   eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl']
                  - eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl'] )*data['TL']
              + (  (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                 - (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*eff['eff_R_L_subl'] ) * data['LT']
              + (   (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])
                  - (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl']) )*data['LL']
              ) )


    nPred_RF_TT = ( (1-eff['eff_R_L_lead'])*(1-eff['eff_F_L_subl'])* 
                      ( (1.0/( (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl'] 
                      + (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                      + (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl']
                      - (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                      - (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl']
                      - (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl'] ) )
              * ( (  eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl']
                   - eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl'] )*data['TL']
                + (  (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*eff['eff_R_L_subl']
                   - (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_F_L_lead']*eff['eff_F_L_subl'] )*data['LT']
                + (   (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])
                    - (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl']) )*data['LL']
              ) ) )

    nPred_FR_TT = ( (1-eff['eff_F_L_lead'])*(1-eff['eff_R_L_subl'])* 
                   ( (1.0/( (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl']     
                    + (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                    + (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl']
                    - (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                    - (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl']
                    - (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl'] ) )
            * ( (   eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                  - eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl'] )*data['TL'] 
              + (   (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*eff['eff_F_L_subl']
                  - (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_R_L_lead']*eff['eff_F_L_subl'] )*data['LT']
              + (   (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])
                  - (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl']) ) *data['LL']
              ) ) )

    nPred_FF_TT = ( (1-eff['eff_F_L_lead'])*(1-eff['eff_F_L_subl'])* 
                  ( (1.0/( (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl']    
                    + (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                    + (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl']
                    - (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                    - (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl']
                    - (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl'] ) )
            * ( (   eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl']
                  - eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl'] )*data['TL']
              + (  (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                 - (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*eff['eff_R_L_subl'] ) * data['LT']
              + (   (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])
                  - (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl']) )*data['LL']
              ) ) )
       
    return {'alpha_RF' : alpha_rf, 'alpha_FR' : alpha_fr, 'alpha_FF' : alpha_ff, 'pred_RF_TT' : nPred_RF_TT, 'pred_FR_TT' : nPred_FR_TT, 'pred_FF_TT' : nPred_FF_TT }


def save_results( results, outputDir, namePostfix='' ) :
    if outputDir is None :
        return

    fname = outputDir + '/results%s.pickle' %namePostfix

    if not os.path.isdir( os.path.dirname( fname ) ) :
        os.makedirs( os.path.dirname( fname ) )
        
    file = open( fname, 'w' )
    pickle.dump( results, file )
    file.close()


####################################
## do_nominal_fit##
InputFile0 = 'TEMP_ROOT/SieieTemp_miniTree_mrgd_data.root' # data
InputFile1 = 'TEMP_ROOT/sigTempSieie_SieieSigTemp_job_2electron_mrgd.root' # signal
InputFile2 = 'TEMP_ROOT/bkgTempSieie_BkgSieieSigTemp_job_2electron_mrgd.root' # bkg

regions = ['EB', 'EE']
for reg in regions :

    # templates
    true_lead_hist = 'h_Sieie_'+reg+'_1'
    fake_lead_hist = 'h_Sieie_'+reg+'_SB_1'
    
    gg_hist_name = 'h_Sieie_'+reg+'_1'
    gg_hist = get1PhoHist(InputFile0, gg_hist_name)
    
    templates = {}
    templates['lead'] = {}
    templates['lead']['real'] = get1PhoHist(InputFile1, true_lead_hist)
    templates['lead']['fake'] = get1PhoHist(InputFile2, fake_lead_hist)
    
    namePostfix = '_%s' %( reg)
    outputDir = 'output/'
    
    (results_stat, results_syst) = run_photon_fit(templates, gg_hist, reg, 'sigmaIetaIeta', outputDir)
    
    save_results(results_stat, outputDir, namePostfix)
    save_results(results_syst, outputDir, namePostfix)


