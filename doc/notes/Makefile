
TEX=latexmk
TEXF=-pdf -pdflatex="pdflatex" -outdir="pdf"
BIB=biber

all: adaptive_timestep ampere_faraday charge_states discretisation dreicer frozen-current heat-transport kinionz_approx magnetic_compression psi0psi1evaluation RunawayTransportRechesterRosenbluth SPIcoordinates SPIDeltaSource SvenssonTransport theory trappedbc

adaptive_timestep: adaptive_timestep.pdf

ampere_faraday: ampere_faraday.pdf

charge_states: charge_states.pdf

discretisation: discretisation.pdf

dreicer: dreicer.pdf

frozen-current: frozen-current.pdf

heat-transport: heat-transport.pdf

kinioniz_approx: kinioniz_approx.pdf

magnetic_compression: magnetic_compression.pdf

psi0psi1evaluation: psi0psi1evaluation.pdf

RunawayTransportRechesterRosenbluth: RunawayTransportRechesterRosenbluth.pdf

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
