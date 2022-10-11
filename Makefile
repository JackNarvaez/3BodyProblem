Body2 = 1 		#Row's number in data of Body 2 
Body3 = 2		#Row's number in data of Body 3
jump = 100		# Stepsize of writing
k = 10000		# Number of orbital periods of outer body.
n = 10000		# Iterations in one orbital period of outer body.

all: Main

Main: Main.cpp Integrator.cpp Integrator.h conversion.cpp conversion.h
	g++ -std=c++17 $< Integrator.cpp conversion.cpp -o Evolution.x;\
	./Evolution.x ${Body2} ${Body3} ${k} ${jump} ${n};\

COEs_Evolve: orb_elements.py
	python3 $< ${Body3}

Phase_Space: Phase_space.py
	python3 $< ${Body3}

Energy: Energy.cpp Integrator.cpp Integrator.h
	g++ -std=c++17 $< Integrator.cpp Integrator.h -o Energy.x;\
	./Energy.x ${Body2} ${Body3} ${k} ${jump} ${n};\

Plot_Energy: EnergyPlot.py
	python3 $< ${Body2} ${Body3}

Profiling_gprof: Main.cpp Integrator.cpp Integrator.h conversion.h conversion.h
	g++ -std=c++17 -g -pg -Wall  $< Integrator.cpp conversion.cpp -o Evolution.x;\
	./Evolution.x ${Body2} ${Body3} 1000 10 1000;\
	gprof ./Evolution.x gmon.out > report.txt;\

#Profiling_cache: Main.cpp Integrator.cpp Integrator.h conversion.h conversion.h
#	g++ -std=c++17 -g -pg -Wall  $< Integrator.cpp conversion.cpp -o Evolution.x;\
#	valgrind --tool=cachegrind ./Evolution.x  ${Body2} ${Body3} 100 10 100;\
#	cg_annotate --auto=yes cachegrind.out.<pid> > report2.txt;\

clean:
	rm -f *.x *.txt *.png *.out *.gif
	rm -f ./Files/* *.txt
