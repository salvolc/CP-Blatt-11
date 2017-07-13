data_a1 = 	a
plots_a1 = 	test
data_a2 = 	
plots_a2 = 	
movies = animation.mp4
movies2 = 

all: $(plots_a1) $(plots_a2) 

animation: $(plots_a1) $(plots_a2) animation.py $(movies)

aufgabe1: aufgabe1.cpp functions.hpp
	g++ -std=c++11 aufgabe1.cpp -o aufgabe1 -O2 -fopenmp 

aufgabe2: aufgabe2.cpp functions.hpp
	g++ -std=c++11 aufgabe2.cpp -o aufgabe2 -O2 -fopenmp 


$(movies): animation.py $(data_a1)
	python3.5 animation.py

$(movies2): animation_2.py $(data_a2)
	python3.5 animation_2.py

$(data_a1): aufgabe1
	./aufgabe1

$(data_a2): aufgabe2
	./aufgabe2


$(plots_a1): $(data_a1) aufgabe1.py
	python3.5 aufgabe1.py

$(plots_a2): $(data_a2) aufgabe2.py
	python3.5 aufgabe2.py


clean:
	@rm -rf $(data_a1) $(data_a2) $(plots_a1) $(plots_a2) $(movies) aufgabe1 aufgabe2 functions 
