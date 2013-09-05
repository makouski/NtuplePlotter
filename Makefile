all: makeTemplates makeSkim

makeTemplates: Histogrammer.o Selector.o EventTree.o makeTemplates.cpp
	g++ -o makeTemplates `root-config --libs` -I`root-config --incdir` EventTree.o Selector.o Histogrammer.o makeTemplates.cpp
##  JetMETObjects/FactorizedJetCorrector.o JetMETObjects/JetCorrectorParameters.o JetMETObjects/SimpleJetCorrector.o

makeSkim: Selector.o EventTree.o makeSkim.cpp
	g++ -o makeSkim `root-config --libs` -I`root-config --incdir` EventTree.o Selector.o makeSkim.cpp

EventTree.o: EventTree.cpp EventTree.h
	g++ -c -I`root-config --incdir` EventTree.cpp

Selector.o: Selector.cpp Selector.h
	g++ -c -I`root-config --incdir` Selector.cpp

Histogrammer.o: Histogrammer.cpp Histogrammer.h
	g++ -c -I`root-config --incdir` Histogrammer.cpp

clean:
	rm EventTree.o Selector.o Histogrammer.o
