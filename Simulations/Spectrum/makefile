particles = 10000000
energy = 100
start = 50
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
	topas output/neutrons.topas >> log.txt
	touch output/results.csv
	julia matcher.jl #>> log.txt
	
curve :
	make clean
	rm -f output/results.csv
	touch output/results.csv
	echo -n "Energy, Mean, STD \n" >> output/results.csv
	for E in $(shell seq $(start) $(interval) $(stop)) ; do \
		for i in $(shell seq 1 $(repeats)) ; do \
			echo -n "$$E, " >> output/results.csv ; \
			echo "Energy: $$E / $(stop) MeV, particles= $(particles), measurement $$i / $(repeats)." ; \
			make -s run energy=$$E ; \
			echo -n "\n" >> output/results.csv ; \
		done ; \
	done
	julia plotter.jl >> log.txt
	
	
clean :
	rm -f log.txt
	rm -f output/*
	

