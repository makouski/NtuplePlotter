from distribution_mod import distribution
import ROOT

# initialize variables, assign values later
WJetsSF = 1.0
TopSF = 1.0
QCDSF = 0.0

#import array
#binarray = array.array('d')
#binarray.fromlist([0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,300])

# load cross-sections and N gen as global variables
execfile('SF.py')

# function definitions ####################################################

def saveTemplatesToFile(templateList, varlist, outFileName):
	outfile = ROOT.TFile(outFileName,'RECREATE')
	for template in templateList:
		for var in varlist:
			template.histList[var].SetDirectory(outfile.GetDirectory(''))
			template.histList[var].Write()
			template.histList[var].SetDirectory(0)
	outfile.Close()

def plotTemplates(dataTemplate, MCTemplateList, SignalTemplateZoomList, varlist, outDirName):
	canvas = ROOT.TCanvas('c1','c1',640,800)
	
	for var in varlist:
		legend = ROOT.TLegend(0.7, 1 - 0.05*(1 + len(MCTemplateList) + len(SignalTemplateZoomList)), 0.99, 1.00)
		legend.SetBorderSize(0)
		legend.SetFillColor(10)
		
		if dataTemplate is not None:
			legend.AddEntry(dataTemplate.histList[var], dataTemplate.name, 'pl')
		
		# MC templates listed in the order they appear in legend
		for mc in MCTemplateList:
			mcHist = mc.histList[var]
			legend.AddEntry(mcHist, mc.name, 'f')
		
		stack = ROOT.THStack('stack_'+var,var)
		# reverse order for stack to be consistent with legend
		MCTemplateList.reverse()
		for mc in MCTemplateList:
			mcHist = mc.histList[var]
			#if var == 'nub':
			#	mcHist.Rebin(5)
			stack.Add(mcHist)
		MCTemplateList.reverse()

		if dataTemplate is not None:
			#if var == 'nub':
			#	dataTemplate.histList[var].Rebin(5)
			if dataTemplate.histList[var].GetMaximum() > stack.GetMaximum():
				stack.SetMaximum(dataTemplate.histList[var].GetMaximum())

		if 'cut_flow' in var: # or 'MET' in var:
			canvas.SetLogy(1)
			stack.SetMinimum(100)
		else:
			canvas.SetLogy(0)
		
		stack.Draw('HIST')
		
		if 'barrel' in outDirName and 'photon1SigmaIEtaIEta' in var:
			stack.GetXaxis().SetRangeUser(0.0,0.025)
		
		if dataTemplate is not None:
			stack.GetXaxis().SetTitle(dataTemplate.histList[var].GetXaxis().GetTitle())
			stack.GetYaxis().SetTitle(dataTemplate.histList[var].GetYaxis().GetTitle())
		stack.SetTitle('')

		if dataTemplate is not None:
			dataTemplate.histList[var].Draw('ESAME')
					
		for signal,zoom in SignalTemplateZoomList:
			sigHist = signal.histList[var].Clone()
			#sigHist.SetFillStyle(3244)
			sigHist.Scale(zoom)
			sigHist.Draw('HISTSAME')
			if zoom != 1:
				legend.AddEntry(sigHist, signal.name + ' x ' + str(zoom), 'f')
			else:
				legend.AddEntry(sigHist, signal.name, 'f')
		if 'cut_flow' not in var:
			legend.Draw()
		canvas.SaveAs(outDirName+'/'+var+'.png')
		
def loadDataTemplate(varlist, inputDir, prefix):
	templPrefix = inputDir+prefix
	DataTempl = distribution('Data', [
		(templPrefix+'Data_a.root', 1),
		(templPrefix+'Data_b.root', 1),
		(templPrefix+'Data_c.root', 1),
		(templPrefix+'Data_d.root', 1),
		], varlist)
	return DataTempl


def loadQCDTemplate(varlist, inputDir, prefix):
	templPrefix = inputDir+prefix
	QCD_sf = QCDSF
	QCDTempl = distribution('QCD', [
		(templPrefix+'Data_a.root', QCD_sf),
		(templPrefix+'Data_b.root', QCD_sf),
		(templPrefix+'Data_c.root', QCD_sf),
		(templPrefix+'Data_d.root', QCD_sf),
		(templPrefix+'TTJets1l.root', -1 * QCD_sf * gSF * TTJets1l_xs/TTJets1l_num),
		(templPrefix+'TTJets2l.root', -1 * QCD_sf * gSF * TTJets2l_xs/TTJets2l_num),
		(templPrefix+'WJets.root', -1 * QCD_sf * gSF * WJets_xs/WJets_num),
	], varlist, 46)
	return QCDTempl


def loadMCTemplates(varList, inputDir, prefix, titleSuffix, fillStyle):
	templPrefix = inputDir+prefix
	
	MCtemplates = {}
	
	#MCtemplates['mst_510_200'] = distribution('mst_510_200'+titleSuffix,[(templPrefix+'mst_510_M3_5050_M1_200.root', gSF*0.0751004/15000)],varList, 620,3244)
	
	MCtemplates['WHIZARD'] = distribution('WHIZARD'+titleSuffix, [(templPrefix+'WHIZARD.root', TopSF*gSF*TTgamma_xs/WHIZARD_num)], varList, 98, fillStyle)
	
	MCtemplates['TTJets'] = distribution('TTJets'+titleSuffix, [
		(templPrefix+'TTJets1l.root', TopSF*gSF*TTJets1l_xs/TTJets1l_num),
		(templPrefix+'TTJets2l.root', TopSF*gSF*TTJets2l_xs/TTJets2l_num),
		], varList ,11, fillStyle)
	
	###################################
	#return MCtemplates
	###################################
	
	MCtemplates['Vgamma'] = distribution('Vgamma'+titleSuffix, [
        (templPrefix+'Zgamma.root', gSF*Zgamma_xs/Zgamma_num),
        (templPrefix+'Wgamma.root', gSF*Wgamma_xs/Wgamma_num),
    #    (templPrefix+'WWgamma.root', gSF*WWgamma_xs/WWgamma_num),
        ], varList, 90, fillStyle)

	MCtemplates['SingleTop'] = distribution('SingleTop'+titleSuffix, [
		(templPrefix+'SingleT_t.root',      gSF*SingTopT_xs/SingTopT_num),
        (templPrefix+'SingleT_s.root',      gSF*SingTopS_xs/SingTopS_num),
        (templPrefix+'SingleT_tw.root',     gSF*SingToptW_xs/SingToptW_num),
        (templPrefix+'SingleTbar_t.root',   gSF*SingTopbarT_xs/SingTopbarT_num),
        (templPrefix+'SingleTbar_s.root',   gSF*SingTopbarS_xs/SingTopbarS_num),
        (templPrefix+'SingleTbar_tw.root',  gSF*SingTopbartW_xs/SingTopbartW_num),
		], varList, 8, fillStyle)
	
	MCtemplates['WJets'] = distribution('WJets'+titleSuffix, [
        (templPrefix+'WJets.root', WJetsSF*gSF*WJets_xs/WJets_num)], varList, 7, fillStyle)
		
	MCtemplates['ZJets'] = distribution('ZJets'+titleSuffix, [
		(templPrefix+'ZJets.root', gSF*ZJets_xs/ZJets_num)], varList, 9, fillStyle)

	
	MCtemplates['Other'] = distribution('Other'+titleSuffix, [

        (templPrefix+'WZ_3lnu.root', gSF*WZ_3lnu_xs/WZ_3lnu_num),
        (templPrefix+'WZ_2l2q.root', gSF*WZ_2l2q_xs/WZ_2l2q_num),
        
        (templPrefix+'ZZ_2e2mu.root', gSF*ZZ_2e2mu_xs/ZZ_2e2mu_num),
        (templPrefix+'ZZ_2e2tau.root', gSF*ZZ_2e2tau_xs/ZZ_2e2tau_num),
        (templPrefix+'ZZ_2mu2tau.root', gSF*ZZ_2mu2tau_xs/ZZ_2mu2tau_num),
        (templPrefix+'ZZ_4e.root', gSF*ZZ_4e_xs/ZZ_4e_num),
        (templPrefix+'ZZ_4mu.root', gSF*ZZ_4mu_xs/ZZ_4mu_num),
        (templPrefix+'ZZ_4tau.root', gSF*ZZ_4tau_xs/ZZ_4tau_num),
        
        (templPrefix+'WW_2l2nu.root', gSF*WW_2l2nu_xs/WW_2l2nu_num),

          #(templPrefix+'TTW.root', gSF*TTW_xs/TTW_num),
          #(templPrefix+'TTZ.root', gSF*TTZ_xs/TTZ_num),
		], varList, 7, fillStyle)

	return MCtemplates


def saveNoMETTemplates(inputDir, outFileName):
	varList = ['MET','M3']
	DataTempl = loadDataTemplate(varList, inputDir, 'hist_1phoNoMET_top_')
	MCTemplDict = loadMCTemplates(varList, inputDir, 'hist_1phoNoMET_top_','',1001)
	MCTempl = []
	MCTempl.append(MCTemplDict['WHIZARD'])
	MCTempl.append(MCTemplDict['TTJets'])
	MCTempl.append(MCTemplDict['Vgamma'])
	MCTempl.append(MCTemplDict['SingleTop'])
	MCTempl.append(MCTemplDict['WJets'])
	MCTempl.append(MCTemplDict['ZJets'])
	MCTempl.append(MCTemplDict['Other'])
	saveTemplatesToFile([DataTempl] + MCTempl, varList, outFileName)


def saveBarrelFitTemplates(inputDir, outFileName):
	varList = ['MET','M3','photon1ChHadSCRIso', 'photon1ChHadRandIso', 'photon1_Sigma_ChSCRIso']
	DataTempl_b = loadDataTemplate(varList, inputDir, 'hist_1pho_barrel_top_')
	
	MCTempl_b = loadMCTemplates(varList, inputDir, 'hist_1pho_barrel_top_','',1001)	
	MCTempl_rs_b = loadMCTemplates(varList, inputDir, 'hist_1pho_rs_barrel_top_', '_signal', 1001)
	MCTempl_fe_b = loadMCTemplates(varList, inputDir, 'hist_1pho_fe_barrel_top_', '_electron', 3005)
	MCTempl_fjrb_b = loadMCTemplates(varList, inputDir, 'hist_1pho_fjrb_barrel_top_', '_fake', 3005)
	
	saveTemplatesToFile([DataTempl_b] + [
		MCTempl_b['WHIZARD'],MCTempl_b['TTJets'],
		MCTempl_rs_b['WHIZARD'],MCTempl_rs_b['TTJets'],
		MCTempl_fe_b['WHIZARD'],MCTempl_fe_b['TTJets'],
		MCTempl_fjrb_b['WHIZARD'],MCTempl_fjrb_b['TTJets'],
	], varList, outFileName)

def savePreselTemplates(inputDir, outFileName):
	if WJetsSF != 1.0 or TopSF != 1.0:
		print 'We want to save templates for M3 fit, but the SFs are not 1.0'
		print 'exiting'
		return
	
	varList = ['MET','M3']
	DataTempl = loadDataTemplate(varList, inputDir, 'hist_1pho_top_')
	if QCDSF > 0.0001:
		QCDTempl = loadQCDTemplate(varList, inputDir[:-1]+'_qcd/', 'hist_1pho_top_')
	else:
		print 'The purpose of this function is to save templates for M3 fit, without QCD it is useless'
	MCTemplDict = loadMCTemplates(varList, inputDir, 'hist_1pho_top_','',1001)
	MCTempl = []
	MCTempl.append(MCTemplDict['WHIZARD'])
	MCTempl.append(MCTemplDict['TTJets'])
	MCTempl.append(MCTemplDict['Vgamma'])
	MCTempl.append(MCTemplDict['SingleTop'])
	MCTempl.append(MCTemplDict['WJets'])
	MCTempl.append(MCTemplDict['ZJets'])
	MCTempl.append(MCTemplDict['Other'])
	if QCDSF > 0.0001:
		MCTempl.append(QCDTempl)
	saveTemplatesToFile([DataTempl] + MCTempl, varList, outFileName)


def makeAllPlots(varList, inputDir, outDirName):
	shortVarList = varList[:]
	shortVarList.remove('cut_flow')

	# load templates PreSel
	
	DataTempl = loadDataTemplate(varList, inputDir, 'hist_1pho_top_')
	if QCDSF > 0.0001:
		QCDTempl = loadQCDTemplate(varList, inputDir[:-1]+'_qcd/', 'hist_1pho_top_')
	MCTemplDict = loadMCTemplates(varList, inputDir, 'hist_1pho_top_','',1001)
	MCTempl = []
	#MCTempl.append(MCTemplDict['mst_510_200'])
	MCTempl.append(MCTemplDict['WHIZARD'])
	MCTempl.append(MCTemplDict['TTJets'])
	MCTempl.append(MCTemplDict['Vgamma'])
	MCTempl.append(MCTemplDict['SingleTop'])
	MCTempl.append(MCTemplDict['WJets'])
	MCTempl.append(MCTemplDict['ZJets'])
	MCTempl.append(MCTemplDict['Other'])
	if QCDSF > 0.0001:
		MCTempl.append(QCDTempl)
	
	if WJetsSF == 1.0 and TopSF == 1.0:
		pass
	else:
		# save final templates, exactly as they are on the plots
		saveTemplatesToFile([DataTempl] + MCTempl, ['MET','M3','WtransMass'], 'templates_presel_scaled.root')
		
	plotTemplates( DataTempl, MCTempl, [], varList, outDirName+'/presel')
	
	region = 'barrel'
	# load templates
	DataTempl_b = loadDataTemplate(shortVarList, inputDir, 'hist_1pho_'+region+'_top_')
	if QCDSF > 0.0001:
		QCDTempl_b = loadQCDTemplate(shortVarList, inputDir[:-1]+'_qcd/', 'hist_1pho_'+region+'_top_')
	MCTemplDict_b = loadMCTemplates(shortVarList, inputDir, 'hist_1pho_'+region+'_top_','',1001)
	MCTempl_b = []
	#MCTempl_b.append(MCTemplDict_b['mst_510_200'])
	MCTempl_b.append(MCTemplDict_b['WHIZARD'])
	MCTempl_b.append(MCTemplDict_b['TTJets'])
	MCTempl_b.append(MCTemplDict_b['Vgamma'])
	MCTempl_b.append(MCTemplDict_b['SingleTop'])
	MCTempl_b.append(MCTemplDict_b['WJets'])
	MCTempl_b.append(MCTemplDict_b['ZJets'])
	MCTempl_b.append(MCTemplDict_b['Other'])
	if QCDSF > 0.0001:
		MCTempl_b.append(QCDTempl_b)
	
	#susysignal = [
	#		(distribution('mst_510_M3_5050_M1_200', [('/Users/makouski/dis/plotting_trees/Hist/hist_1pho_barrel_top_mst_510_M3_5050_M1_200.root', gSF*0.0751004/15000)],shortVarList,620,3244),1),
	#		#(,1),
	#		#(,1)
	#		]
	if WJetsSF != 1.0 or TopSF != 1.0:
		# save final templates, exactly as they are on the plots
		saveTemplatesToFile([DataTempl_b] + MCTempl_b, ['MET','M3','WtransMass'], 'templates_barrel_scaled.root')
	plotTemplates( DataTempl_b, MCTempl_b, [], shortVarList, outDirName+'/'+region+'_samples')
	
	############################
	return
	############################
	
	MCTempl_rs_c = loadMCTemplates(shortVarList, inputDir, 'hist_1pho_rs_'+region+'_top_', '_signal', 1001)
	MCTempl_fe_c = loadMCTemplates(shortVarList, inputDir, 'hist_1pho_fe_'+region+'_top_', '_electron', 3005)
	MCTempl_fjrb_c = loadMCTemplates(shortVarList, inputDir, 'hist_1pho_fjrb_'+region+'_top_', '_fake', 3005)
	
	# combine MC categories
	TTphoton = MCTempl_rs_c['WHIZARD']
	TTphoton.mergeWith(MCTempl_rs_c['TTJets'])
	TTphoton.name = 'TopPair+photon'
	
	TTfake = MCTempl_fe_c['WHIZARD']
	TTfake.mergeWith(MCTempl_fjrb_c['WHIZARD'])
	TTfake.mergeWith(MCTempl_fe_c['TTJets'])
	TTfake.mergeWith(MCTempl_fjrb_c['TTJets'])
	TTfake.name = 'TopPair+fake'
	
	BGphoton = MCTempl_rs_c['Vgamma']
	BGphoton.mergeWith(MCTempl_rs_c['SingleTop'])
	#BGphoton.mergeWith(MCTempl_rs_c['Other'])
	BGphoton.mergeWith(MCTempl_rs_c['WJets'])
	BGphoton.mergeWith(MCTempl_rs_c['ZJets'])
	BGphoton.name = 'BG+photon'
	
	BGfake = MCTempl_fe_c['Vgamma']
	BGfake.mergeWith(MCTempl_fjrb_c['Vgamma'])
	BGfake.mergeWith(MCTempl_fe_c['SingleTop'])
	BGfake.mergeWith(MCTempl_fjrb_c['SingleTop'])
	#BGfake.mergeWith(MCTempl_fe_c['Other'])
	#BGfake.mergeWith(MCTempl_fjrb_c['Other'])
	BGfake.mergeWith(MCTempl_fe_c['WJets'])
	BGfake.mergeWith(MCTempl_fjrb_c['WJets'])
	BGfake.mergeWith(MCTempl_fe_c['ZJets'])
	BGfake.mergeWith(MCTempl_fjrb_c['ZJets'])
	BGfake.name = 'BG+fake'

	
	MCTempl_c = []
	MCTempl_c.append(TTphoton)
	MCTempl_c.append(TTfake)
	MCTempl_c.append(BGphoton)
	MCTempl_c.append(BGfake)
	
	plotTemplates( DataTempl_b, MCTempl_c, [], shortVarList, outDirName+'/'+region+'_categories')

# main part ##############################################################################################
# templates for data drive fit or closure test. No rescaling necessary
#saveBarrelFitTemplates('/Users/makouski/dis/plotting_trees/hist/', 'templates_barrel.root')

# for MET fit. No rescaling
#if WJetsSF == 1.0 and TopSF == 1.0:
#	saveNoMETTemplates('/Users/makouski/dis/plotting_trees/hist/', 'templates_presel_nomet.root')
#	saveNoMETTemplates('/Users/makouski/dis/plotting_trees/hist_qcd/', 'templates_presel_nomet_qcd.root')

# after doing MET fit, update the QCDSF

# ==3j 
#QCDSF = 0.14

# >=4j
QCDSF = 0.1048

# save templates for M3 fit
savePreselTemplates('/Users/makouski/dis/plotting_trees/hist/', 'templates_presel.root')
# do M3 fit, update SF for Top and WJets

# ==3j 
#TopSF = 1.15
#WJetsSF = 1.10

# >=4j
TopSF = 0.936
WJetsSF = 1.59

varList_all = ['nVtx',
			'MET','Ht','WtransMass','M3','M3first','minM3','enub','eb','nub','enubjj',
			'ele1Pt','ele1Eta','ele1RelIso',
			'ele1D0','ele1MVA','ele1Dz',
			'ele2Pt','ele2RelIso',
			'ele1sigmaIetaIeta','ele1EoverP',
			'ele1DrJet','ele1pho1Mass',
			'looseEleDrGenPho',
			'cut_flow',
			'nJets',
			'jet1Pt','jet2Pt','jet3Pt','jet4Pt','jet1Eta','jet2Eta','jet3Eta','jet4Eta',
			'photon1Et','photon1Eta','photon1HoverE','photon1SigmaIEtaIEta',
			'photon1DrElectron','photon1DrJet',
			'photon1ChHadIso','photon1NeuHadIso','photon1PhoIso',
			'photon1ChHadSCRIso','photon1PhoSCRIso',
			'photon1ChHadRandIso','photon1PhoRandIso',
			'photon1MotherID','photon1GMotherID','photon1DrMCbquark','GenPhotonEt',
			#'photon1_Sigma_ChSCRIso'
			]

makeAllPlots(varList_all, '/Users/makouski/dis/plotting_trees/hist/', 'plots')

