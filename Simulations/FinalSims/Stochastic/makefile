MAKEFLAGS += --no-print-directory
SHELL := /bin/bash

particles = 10000000
energy = 100
start = 1
stop = 25
interval = 1
repeats = 10
vis = 0
z = 25
dets = 24
time = 1000000

default :
	@echo "Call something specific, like 'make files particles=50000 energy=200'"


files :
	@cp -f topastemplate.txt output/neutrons.topas
	@julia writer.jl $(particles) $(energy) $(vis) >> log.txt

run :
	make files
	julia -t auto geometrywriter.jl $(dets) $(z)
	topas output/neutrons.topas #>> log.txt
	touch output/results.csv
	julia -t auto matcher.jl $(time) #>> log.txt
	julia -t auto efficienter.jl $(particles) 1 0
	
many :
	make clean
	for i in $$(seq 1 1 $(repeats)) ; do \
		make files ; \
		echo "run $$i / $(repeats)" ; \
		topas output/neutrons.topas >> log.txt ; \
		julia -t auto matcher.jl $(time)  >> log.txt ; \
	done
	julia -t auto efficienter.jl $(particles) $(repeats) 0

stochastic :
	make clean
	for i in $$(seq 1 1 $(repeats)) ; do \
		make files ; \
		echo "run $$i / $(repeats)" ; \
		topas output/neutrons.topas >> log.txt ; \
		julia -t auto stochasticmatcher.jl >> log.txt ; \
	done
	julia -t auto stochasticefficienter.jl $(particles) $(repeats) 0
	
clean :
	rm -f log.txt
	rm -f output/*
	rm -f plots/*
	rm -f results.csv


test :
	for i in $$(seq 2 10); do \
		echo "$$i"; \
	done
