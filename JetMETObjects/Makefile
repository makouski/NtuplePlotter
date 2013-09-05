all: FactorizedJetCorrector.o JetCorrectorParameters.o SimpleJetCorrector.o

clean:
	rm FactorizedJetCorrector.o JetCorrectorParameters.o SimpleJetCorrector.o

FactorizedJetCorrector.o: FactorizedJetCorrector.cc FactorizedJetCorrector.h JetCorrectorParameters.h SimpleJetCorrector.h Utilities.cc
	g++ -c `root-config --libs` -I`root-config --incdir` FactorizedJetCorrector.cc

JetCorrectorParameters.o: JetCorrectorParameters.cc JetCorrectorParameters.h 
	g++ -c `root-config --libs` -I`root-config --incdir` JetCorrectorParameters.cc

SimpleJetCorrector.o: JetCorrectorParameters.cc JetCorrectorParameters.h SimpleJetCorrector.h
	g++ -c `root-config --libs` -I`root-config --incdir` SimpleJetCorrector.cc 
 
