import ROOT
import templateFits

templateFits.InputFilename = 'templates_barrel.root'
templateFits.fitData = False ## to do closure test
templateFits.NpseudoExp = 0

#templateFits.datasetsToMix = ['TTGamma','TTJets']
#templateFits.datasetsToMix = ['MGTTGamma','TTJets']

scales = [0.1, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 5.0]
#scales = [1.0]
graph = ROOT.TGraphErrors(len(scales))
gpoint = 0

for s in scales:
	templateFits.signalSF = s
	phoPurity,phoPurityError,MCtruth = templateFits.doTheFit()
	
	graph.SetPoint(gpoint, MCtruth, phoPurity)
	graph.SetPointError(gpoint, 0.0, phoPurityError)
	gpoint+=1
	
canvas = ROOT.TCanvas('c1','c1',640,800)
graph.Draw('Ap')
graph.SetMarkerStyle(8)

line = ROOT.TLine(0.0, 0.0, 1.0, 1.0)
line.Draw('SAME')

canvas.SaveAs('linearityNewMG.png')