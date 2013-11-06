from distribution_mod import distribution
import ROOT

# load cross-sections and N gen as global variables
execfile('SF.py')

# function definitions ####################################################

def saveTemplatesToFile(templateList, varList, outFileName):
	outfile = ROOT.TFile(outFileName,'RECREATE')
	for template in templateList:
		for var in varList:
			template.histList[var].SetDirectory(outfile.GetDirectory(''))
			template.histList[var].Write()
			template.histList[var].SetDirectory(0)
	outfile.Close()

def plotTemplates(dataTemplate, MCTemplateList, SignalTemplateZoomList, varList, outDirName):
	canvas = ROOT.TCanvas('c1','c1',640,800)
	
	for var in varList:
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
			stack.Add(mcHist)
		MCTemplateList.reverse()

		if dataTemplate is not None:
			if dataTemplate.histList[var].GetMaximum() > stack.GetMaximum():
				stack.SetMaximum(dataTemplate.histList[var].GetMaximum())

		if 'cut_flow' in var:
			canvas.SetLogy(1)
			stack.SetMinimum(100)
		else:
			canvas.SetLogy(0)
		
		stack.Draw('HIST')
		
		if 'barrel' in outDirName and 'photon1SigmaIEtaIEta' in var:
			stack.GetXaxis().SetRangeUser(0.0,0.025)
		
		stack.GetXaxis().SetTitle(dataTemplate.histList[var].GetXaxis().GetTitle())
		stack.GetYaxis().SetTitle(dataTemplate.histList[var].GetYaxis().GetTitle())
		stack.SetTitle('')

		if dataTemplate is not None:
			#dataTemplate.histList[var].SetMarkerStyle(8)
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


def loadDataTemplate(varList, inputDir, prefix):
	templPrefix = inputDir+prefix
	DataTempl = distribution('Data', [
		(templPrefix+'Data_a_Aug6.root', 1),
		(templPrefix+'Data_a_Jul13.root', 1),
		(templPrefix+'Data_b_Jul13.root', 1),
		(templPrefix+'Data_c_Aug24.root', 1),
		(templPrefix+'Data_c_Dec11.root', 1),
		(templPrefix+'Data_c_PRv2.root', 1),
		(templPrefix+'Data_d_PRv1_part1.root', 1),
		(templPrefix+'Data_d_PRv1_part2.root', 1),
		(templPrefix+'Data_d_PRv1_part3.root', 1),
		(templPrefix+'Data_d_PRv1_part4.root', 1),
		(templPrefix+'Data_d_PRv1_part5.root', 1),
		], varList)
	return DataTempl


def loadMCTemplates(varList, inputDir, prefix, titleSuffix, fillStyle):
	templPrefix = inputDir+prefix
	
	WHIZARDTempl = distribution('WHIZARD'+titleSuffix, [(templPrefix+'WHIZARD.root', gSF*TTgamma_xs/WHIZARD_num)], varList, 98, fillStyle)
	
	TTJetsTempl = distribution('TTJets'+titleSuffix, [
		(templPrefix+'TTJets1l.root', gSF*TTJets1l_xs/TTJets1l_num),
		(templPrefix+'TTJets2l.root', gSF*TTJets2l_xs/TTJets2l_num),
		], varList ,11, fillStyle)
	
	#VgammaTempl = distribution('Vgamma'+titleSuffix, [
    #    (templPrefix+'Zgamma.root', gSF*Zgamma_xs/Zgamma_num),
    #    (templPrefix+'Wgamma.root', gSF*Wgamma_xs/Wgamma_num),
    #    (templPrefix+'WWgamma.root', gSF*WWgamma_xs/WWgamma_num),
    #    ], varList, 90, fillStyle)

	OtherTempl = distribution('Other'+titleSuffix, [
        (templPrefix+'WJets1.root', gSF*WJets_xs/(WJets1_num+WJets2_num)),
        (templPrefix+'WJets2.root', gSF*WJets_xs/(WJets1_num+WJets2_num)),

		(templPrefix+'ZJets.root', gSF*ZJets_xs/ZJets_num),

        (templPrefix+'Zgamma.root', gSF*Zgamma_xs/Zgamma_num),
        (templPrefix+'Wgamma.root', gSF*Wgamma_xs/Wgamma_num),
        (templPrefix+'WWgamma.root', gSF*WWgamma_xs/WWgamma_num),

		(templPrefix+'SingleT_t.root',      gSF*SingTopT_xs/SingTopT_num),
        (templPrefix+'SingleT_s.root',      gSF*SingTopS_xs/SingTopS_num),
        (templPrefix+'SingleT_tw.root',     gSF*SingToptW_xs/SingToptW_num),
        (templPrefix+'SingleTbar_t.root',   gSF*SingTopbarT_xs/SingTopbarT_num),
        (templPrefix+'SingleTbar_s.root',   gSF*SingTopbarS_xs/SingTopbarS_num),
        (templPrefix+'SingleTbar_tw.root',  gSF*SingTopbartW_xs/SingTopbartW_num),
		
        (templPrefix+'WZ_3lnu.root', gSF*WZ_3lnu_xs/WZ_3lnu_num),
        (templPrefix+'WZ_2l2q.root', gSF*WZ_2l2q_xs/WZ_2l2q_num),
        
        (templPrefix+'ZZ_2e2mu.root', gSF*ZZ_2e2mu_xs/ZZ_2e2mu_num),
        (templPrefix+'ZZ_2e2tau.root', gSF*ZZ_2e2tau_xs/ZZ_2e2tau_num),
        (templPrefix+'ZZ_2mu2tau.root', gSF*ZZ_2mu2tau_xs/ZZ_2mu2tau_num),
        (templPrefix+'ZZ_4e.root', gSF*ZZ_4e_xs/ZZ_4e_num),
        (templPrefix+'ZZ_4mu.root', gSF*ZZ_4mu_xs/ZZ_4mu_num),
        (templPrefix+'ZZ_4tau.root', gSF*ZZ_4tau_xs/ZZ_4tau_num),
        
        (templPrefix+'WW_2l2nu.root', gSF*WW_2l2nu_xs/WW_2l2nu_num),

        (templPrefix+'TTW.root', gSF*TTW_xs/TTW_num),
        (templPrefix+'TTZ.root', gSF*TTZ_xs/TTZ_num),
		], varList, 7, fillStyle)
	
	MCtemplates = []
	MCtemplates.append(WHIZARDTempl)
	MCtemplates.append(TTJetsTempl)
	#MCtemplates.append(VgammaTempl)
	MCtemplates.append(OtherTempl)
	return MCtemplates

def makeAllPlots(varList, inputDir, outDirName):
	shortVarList = varList[:]
	shortVarList.remove('cut_flow')
	shortVarList.remove('Nele_before_cut')

	
	# load templates (barrel+endcap)
	
	DataTempl = loadDataTemplate(varList, inputDir, 'hist_1pho_top_')
	
	MCTempl_rs = loadMCTemplates(varList, inputDir, 'hist_1pho_rs_top_', '(signal)', 3001)
	#MCTempl_rb = loadMCTemplates(varList, inputDir, 'hist_1pho_rb_top_', '(rb)', 3003)
	MCTempl_fe = loadMCTemplates(varList, inputDir, 'hist_1pho_fe_top_', '(electron)', 3006)
	MCTempl_fjrb = loadMCTemplates(varList, inputDir, 'hist_1pho_fjrb_top_', '(fake)', 3007)
	# combine MC categories in correct order
	MCTempl = []
	for templInd in xrange(len(MCTempl_rs)):
		MCTempl.append(MCTempl_rs[templInd])
		#MCTempl.append(MCTempl_rb[templInd])
		MCTempl.append(MCTempl_fe[templInd])
		MCTempl.append(MCTempl_fjrb[templInd])
	
	plotTemplates( DataTempl, MCTempl, [], varList, outDirName)
	
	
	# load templates (barrel)
	
	DataTempl_b = loadDataTemplate(shortVarList, inputDir, 'hist_1pho_barrel_top_')
	
	MCTempl_rs_b = loadMCTemplates(shortVarList, inputDir, 'hist_1pho_rs_barrel_top_', '(signal)', 3001)
	#MCTempl_rb_b = loadMCTemplates(shortVarList, inputDir, 'hist_1pho_rb_barrel_top_', '(rb)', 3003)
	MCTempl_fe_b = loadMCTemplates(shortVarList, inputDir, 'hist_1pho_fe_barrel_top_', '(electron)', 3006)
	MCTempl_fjrb_b = loadMCTemplates(shortVarList, inputDir, 'hist_1pho_fjrb_barrel_top_', '(fake)', 3007)
	# combine MC categories in correct order
	MCTempl_b = []
	for templInd in xrange(len(MCTempl_rs_b)):
		MCTempl_b.append(MCTempl_rs_b[templInd])
		#MCTempl_b.append(MCTempl_rb_b[templInd])
		MCTempl_b.append(MCTempl_fe_b[templInd])
		MCTempl_b.append(MCTempl_fjrb_b[templInd])
	
	plotTemplates( DataTempl_b, MCTempl_b, [], shortVarList, outDirName+'/barrel')
	
	# save barrel sihih templates for fit
	saveTemplatesToFile([DataTempl_b]+MCTempl_b, ['photon1SigmaIEtaIEta'], 'templates_sihih_barrel.root')

	# load templates (endcap)
	
	DataTempl_e = loadDataTemplate(shortVarList, inputDir, 'hist_1pho_endcap_top_')
	
	MCTempl_rs_e = loadMCTemplates(shortVarList, inputDir, 'hist_1pho_rs_endcap_top_', '(signal)', 3001)
	#MCTempl_rb_e = loadMCTemplates(shortVarList, inputDir, 'hist_1pho_rb_endcap_top_', '(rb)', 3003)
	MCTempl_fe_e = loadMCTemplates(shortVarList, inputDir, 'hist_1pho_fe_endcap_top_', '(electron)', 3006)
	MCTempl_fjrb_e = loadMCTemplates(shortVarList, inputDir, 'hist_1pho_fjrb_endcap_top_', '(fake)', 3007)
	# combine MC categories in correct order
	MCTempl_e = []
	for templInd in xrange(len(MCTempl_rs_e)):
		MCTempl_e.append(MCTempl_rs_e[templInd])
		#MCTempl_e.append(MCTempl_rb_e[templInd])
		MCTempl_e.append(MCTempl_fe_e[templInd])
		MCTempl_e.append(MCTempl_fjrb_e[templInd])
	
	plotTemplates( DataTempl_e, MCTempl_e, [], shortVarList, outDirName+'/endcap')

    
# main part ##############################################################################################

varList = ['nVtx',
			'MET','Ht','WtransMass','WtransMass_photon',

			'ele1Pt','ele1Eta','ele1RelIso',
			'ele1D0','ele1MVA','ele1Dz',
			'ele2Pt','ele2RelIso',
			'ele1sigmaIetaIeta','ele1EoverP',
			'ele1DrJet','ele1pho1Mass',
			'looseEleDrGenPho',

			'cut_flow',
			'Nele_before_cut',
			'jjg_inv_mass','jje_inv_mass','bjjg_inv_mass','bjje_inv_mass','ele_met_g_trans_mass','Njets_vetoed',
           
			'nJets',
			'jet1Pt','jet1CHF','jet1CEF','jet1NHF','jet1NEF','jet1NCharged','jet1NConstituents',
			'jet2Pt','jet3Pt','jet4Pt','jet1Eta','jet2Eta','jet3Eta','jet4Eta',

			'photon1Et','photon1Eta','photon1HoverE','photon1SigmaIEtaIEta',
			'photon1DrElectron','photon1DrJet',
			'photon1RelIso','photon1ChHadIso','photon1NeuHadIso','photon1PhoIso',
			'photon1MotherID','photon1GMotherID','photon1DrMCbquark','GenPhotonEt',
			]

makeAllPlots(varList, '/Users/makouski/dis/plotting_trees/Hist/', 'plots')
#makeAllPlots(varList, '/uscms_data/d2/makouski/ttgamma/plotting/Hist/', 'plots')

