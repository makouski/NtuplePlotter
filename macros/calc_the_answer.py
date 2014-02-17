import ROOT

def integral(file,histname):
	err = ROOT.Double(0.0)
	hist = file.Get(histname)
	integr = hist.IntegralAndError(0,hist.GetNbinsX() + 1, err)
	print integr,err
	return integr,err

## >=4j
photnPurity = 0.5706
photnPurityErr = 0.073

QCDSF = 0.1048
QCDSFErr = 0.0086

M3TopSF = 0.933
M3TopSFErr = 0.0066

M3WJetsSF = 1.61
M3WJetsSFErr = 0.064

## >=3j



print 'input data'
print 'photon Purity: ',photnPurity,' +-',photnPurityErr
print
print 'Opening file with pre-selection histograms'
preselfile = ROOT.TFile('templates_presel_scaled.root','READ')
print 'File was last modified: ',preselfile.GetModificationDate().AsString()
print
topPreselInt,topPreselErr = integral(preselfile,'TTJets_MET')
print 'TTJets presel events expected ',topPreselInt,' +-',topPreselErr,'(MC stat)'
topPreselErr = (topPreselErr**2 + (topPreselInt*M3TopSFErr/M3TopSF)**2 )**0.5
print 'taking into account fit error: ',topPreselInt,' +-',topPreselErr,'(MC stat+fit)'
print
whizPreselInt,whizPreselErr = integral(preselfile,'WHIZARD_MET')
print 'WHIZARD presel events expected ',whizPreselInt,' +-',whizPreselErr,'(MC stat)'
print
print 'Opening file with photon in barrel'
barrelfile = ROOT.TFile('templates_barrel_scaled.root','READ')
print 'File was last modified: ',barrelfile.GetModificationDate().AsString()
print
whizBarrInt,whizBarrErr = integral(barrelfile,'WHIZARD_MET')
print 'WHIZARD barrel events expected ',whizBarrInt, ' +-',whizBarrErr,'(MC stat)'
print

phoAcc = whizBarrInt/whizPreselInt
# These errors are not completely independent. Treat them if they are for now.
phoAccErr = phoAcc * ((whizBarrErr/whizBarrInt)**2 + (whizPreselErr/whizPreselInt)**2)**0.5
print 'signal photon acceptance',phoAcc,' +-',phoAccErr

print

DataBarrInt,DataBarrErr = integral(barrelfile,'Data_MET')
print 'Data events',DataBarrInt,' +-',DataBarrErr

print 'Now getting expected background events from MC'
VgammaBarrInt,VgammaBarrErr = integral(barrelfile,'Vgamma_MET')
SingleTopBarrInt,SingleTopBarrErr = integral(barrelfile,'SingleTop_MET')
WJetsBarrInt,WJetsBarrErr = integral(barrelfile,'WJets_MET')
ZJetsBarrInt,ZJetsBarrErr = integral(barrelfile,'ZJets_MET')
OtherBarrInt,OtherBarrErr = integral(barrelfile,'Other_MET')
QCDBarrInt,QCDBarrErr = integral(barrelfile,'QCD_MET')
bckgExp = VgammaBarrInt+SingleTopBarrInt+WJetsBarrInt+ZJetsBarrInt+OtherBarrInt+QCDBarrInt
SqBckgErr = 0.0
for bc,bcerr in [(VgammaBarrInt,VgammaBarrErr),
					(SingleTopBarrInt,SingleTopBarrErr),
					(WJetsBarrInt, (WJetsBarrErr**2 + (WJetsBarrInt*M3WJetsSFErr/M3WJetsSF)**2)**0.5 ),
					(ZJetsBarrInt,ZJetsBarrErr),
					(OtherBarrInt,OtherBarrErr),
					(QCDBarrInt, (QCDBarrErr**2 + (QCDBarrInt*QCDSFErr/QCDSF)**2)**0.5 )
				]:
	SqBckgErr += (bcerr)**2
bckgExpErr = (SqBckgErr)**0.5
print 'total background expected',bckgExp,' +-',bckgExpErr
print
xsRatio = (DataBarrInt - bckgExp) * photnPurity / phoAcc / topPreselInt
xsRatioRelErr = ( (DataBarrErr**2 + bckgExpErr**2)/(DataBarrInt - bckgExp)**2 + (photnPurityErr/photnPurity)**2 + (phoAccErr/phoAcc)**2 + (topPreselErr/topPreselInt)**2 )**0.5
print xsRatio,' +-',xsRatio*xsRatioRelErr
