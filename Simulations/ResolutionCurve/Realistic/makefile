MAKEFLAGS += --no-print-directory

particles = 10000000
energy = 100
start = 10
stop = 250
interval = 10
repeats = 10
vis = 0



default :
	@echo "Call something specific, like 'make files particles=50000 energy=200'"


files :
	cp -f topastemplate.txt output/neutrons.topas
	julia writer.jl $(particles) $(energy) $(vis) >> log.txt

run :
	make files
	topas output/neutrons.topas #>> log.txt
	touch output/results.csv
	julia matcher.jl #>> log.txt
	julia plotter.jl
	
many :
	make clean
	for i in $(shell seq 1 1 $(repeats)) ; do \
		make files ; \
		echo "run $$i / $(repeats)" ; \
		topas output/neutrons.topas >> log.txt ; \
		julia matcher.jl >> log.txt ; \
	done
	julia plotter.jl
	
curve :
	make clean
	for E in $(shell seq $(start) $(interval) $(stop)) ; do \
		for i in $(shell seq 1 1 $(repeats)) ; do \
			echo "Energy: $$E / $(stop) MeV, particles= $(particles), measurement $$i / $(repeats)." ; \
			make files energy=$$E ; \
			topas output/neutrons.topas >> log.txt ; \
			julia matcher.jl >> log.txt ; \
		done ; \
		julia plotter.jl ; \
	done
	
clean :
	rm -f log.txt
	rm -f output/*
	rm -f plots/*
	rm -f results.csv
	

