Body2 = 1 		#Row's number in data of Body 2 
Body3 = 2		#Row's number in data of Body 3
jump = 100		# Stepsize of writing
k = 10000		# Number of orbital periods of outer body.
n = 10000		# Iterations in one orbital period of outer body.

all: Main

Main: Main.cpp Integrator.cpp Integrator.h conversion.h conversion.h
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

clean:
	rm -f *.x *.txt *.png *.out *.gif
	rm -f ./Files/* *.txt
