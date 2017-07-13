//#include "../eigen/Eigen/Dense"
//#include <iostream>
//#include <fstream>
//#include <complex>
#include "functions.hpp"
// #include "oldfunctions.hpp"
#include <random>
#include <cmath>
//using Eigen::MatrixXd;
//using Eigen::MatrixXcd;
//using Eigen::VectorXd;
//using Eigen::VectorXcd;
//using std::cout;
//using std::endl;


int main(){
	
	std::random_device rd;
    std::mt19937 gen(rd());
 	double alpha = 2;
    std::exponential_distribution<> d(1);

    d.min();

    int n = 1e5;

	VectorXd z(n),zexp(n);

	for (int i = 0; i < int(n/2); ++i)
	{
		zexp[i] = alpha*d(gen);
		zexp[int(n/2)+i] = -alpha*d(gen);
	}
	cout << d(gen) << endl;//gen.max() << endl;*.
	speichere("test",zexp);

	return 0;

}
