MAKEFLAGS += --no-print-directory
SHELL := /bin/bash

particles = 10000000
energy = 100
start = 1
stop = 25
interval = 1
repeats = 10
vis = 0



default :
	@echo "Call something specific, like 'make files particles=50000 energy=200'"


files :
	@cp -f topastemplate.txt output/neutrons.topas
	@julia writer.jl $(particles) $(energy) $(vis) >> log.txt

run :
	make files
	topas output/neutrons.topas #>> log.txt
	touch output/results.csv
	julia -t auto matcher.jl #>> log.txt
	julia -t auto efficienter.jl $(particles) 1
	
many :
	make clean
	for i in $$(seq 1 1 $(repeats)) ; do \
		make files ; \
		echo "run $$i / $(repeats)" ; \
		topas output/neutrons.topas >> log.txt ; \
		julia -t auto matcher.jl >> log.txt ; \
	done
	julia -t auto efficienter.jl $(particles) $(repeats)

curve :
	make clean
	rm -f geometryresults.txt
	touch geometryresults.txt
	echo "Scintillators, matchrate (match / (neutron/cm²))" >> geometryresults.txt

	for E in $$(seq $(start) $(interval) $(stop)) ; do \
		for i in $$(seq 1 1 $(repeats)) ; do \
			echo "Scintillators: $$E / $(stop), particles= $(particles), measurement $$i / $(repeats)." ; \
			julia geometrywriter.jl $$E ; \
			make files energy=$(energy) ; \
			topas output/neutrons.topas >> log.txt ; \
			julia -t auto matcher.jl >> log.txt ; \
		done ; \
		echo -n "$$E, " >> geometryresults.txt ; \
		julia -t auto efficienter.jl $(particles) $(repeats) >> log.txt ; \
	done
	
clean :
	rm -f log.txt
	rm -f output/*
	rm -f plots/*
	rm -f results.csv


test :
	for i in $$(seq 2 10); do \
		echo "$$i"; \
	done
