MAKEFLAGS += --no-print-directory
SHELL := /bin/bash

particles = 10000000
energy = 230
start = 25
stop = 250
interval = 25
repeats = 10
vis = 0
multi = 5
dist = 0.1


default :
	@echo "Call something specific, like 'make files particles=50000 energy=200'"


files :
	@cp -f topastemplate.txt output/neutrons.topas
	@julia writer.jl $(particles) $(energy) $(vis) $(dist) >> log.txt

run :
	make files
	topas output/neutrons.topas #>> log.txt
	touch output/results.csv
	julia -t auto matcher.jl #>> log.txt
	julia -t auto plotter.jl
	
many :
	make clean
	for i in $$(seq 1 1 $(repeats)) ; do \
		make files ; \
		echo "run $$i / $(repeats)" ; \
		topas output/neutrons.topas >> log.txt ; \
		julia -t auto matcher.jl >> log.txt ; \
	done
	julia -t auto plotter.jl

clean :
	rm -f log.txt
	rm -rf output/*
	mkdir output/PhaseSpaces
	rm -f plots/*
	rm -f results.csv


test :
	for i in $$(seq 2 10); do \
		echo "$$i"; \
	done
