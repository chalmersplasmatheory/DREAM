
TEX=latexmk
TEXF=-pdf -pdflatex="pdflatex" -outdir="pdf"
BIB=biber

all: adaptive_timestep ampere_faraday discretisation dreicer frozen-current heat-transport magnetic_compression psi0psi1evaluation SPIcoordinates SPIDeltaSource SvenssonTransport theory trappedbc

adaptive_timestep: adaptive_timestep.pdf

ampere_faraday: ampere_faraday.pdf

discretisation: discretisation.pdf

dreicer: dreicer.pdf

frozen-current: frozen-current.pdf

heat-transport: heat-transport.pdf

magnetic_compression: magnetic_compression.pdf

psi0psi1evaluation: psi0psi1evaluation.pdf

SPIcoordinates: SPIcoordinates.pdf

SPIDeltaSource: SPIDeltaSource.pdf

SvenssonTransport: SvenssonTransport.pdf

theory: theory.pdf

trappedbc: trappedbc.pdf

%.pdf: %.tex
	$(TEX) $(TEXF) -use-make $<
	
clean:
	rm -f *.fdb_latexmk *.fls
	rm -rf pdf/*
