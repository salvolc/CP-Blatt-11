#ifndef FUNCTIONS
#define FUNCTIONS

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "../../eigen/Eigen/Dense"
#include <functional>

using std::complex;
using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;
using std::ofstream;
using std::cout;
using std::endl;

void ladebalken(int i, int max);
void init_random_one_d(Eigen::VectorXd &r, double min, double max);
void speichere(std::string name, MatrixXd data);
void speichere(std::string name, VectorXd data);

void init_linspace_f(Eigen::VectorXd &r, double min, double max){
	int size = r.size();
	double inter = max - min;
	for (int i = 0; i < size; ++i)
	{
		r[i] = min + (inter/(size-1)) * i;
	}
	return;
}

void ladebalken(int i, int max){
	double progress = (1.*i)/(1.*max)*100;
	std::cout << "\rSchritt " << i << " von " << max << " geschafft! " << "Also " << progress << " %";
	if(i == max-1 || i==max){std::cout << "\rSchritt " << max << " von " << max << " geschafft. Fertig!" << std::endl;}
	
	return;
}

void init_random_one_d(Eigen::VectorXd &r, double min, double max){
	int size = r.size();
	double inter = max - min;
	for (int i = 0; i < size; ++i)
	{
		r[i] = min + (rand()%11)/10*inter;
	}
	return;
}
void speichere(std::string name, MatrixXd data){
	cout << ("\nSpeichere Datei " +name+ ".txt ab:").c_str() << endl;
	ofstream dat((name+".txt").c_str());
	dat.is_open();
	dat << data << "\n";
	dat.close();
	cout << "Datei " +name+ ".txt abgespeichert!\n" << endl;
	return;
}

void speichere(std::string name, VectorXd data){
	cout << ("\nSpeichere Datei " +name+ ".txt ab:").c_str() << endl;
	ofstream dat((name+".txt").c_str());
	dat.is_open();
	dat << data << "\n";
	dat.close();
	cout << "Datei " +name+ ".txt abgespeichert!\n" << endl;
	return;
}




#endif 