#include "../../eigen/Eigen/Dense"
#include <iostream>
#include <fstream>
#include <math.h> 
#include <cmath>
#include <functional>
#include <random>
#include "functions.hpp"

using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;
using std::cout;
using std::endl;
using std::ofstream;
using std::complex;
using std::function;

VectorXd parallelZufall(int n);
double pi_MC(int n, int dim, int threads);

int main()
{
	std::random_device rd;  
   	std::mt19937 gen(rd()); 
   	std::uniform_real_distribution<> dis(0, 1);

   	//int n = 1e9;
   	//VectorXd z;
   	//z = parallelZufall(n);

   	double pi=pi_MC(1e9,3,8);

   	cout << pi << endl;

	return 0;
}



VectorXd parallelZufall(int n){
	std::random_device rd;  
   	std::mt19937 gen(rd()); 
   	std::uniform_real_distribution<> dis(0, 1);

	VectorXd x(n);
	int n4 = int((n/4)+0.5);

	MatrixXd x4(4,n4);


	#pragma omp parallel for
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < n4; ++j)
		{
			x[j+i*n4] = dis(gen);
		}
	}
	return x;
}

int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

double pi_MC(int n, int dim=3, int threads=4){
	int NKreis=0;
	int punkteEcht=0;

	int n4 = int((n/threads)+0.5);


	MatrixXd punkt(threads,dim);

	std::random_device rd;  
   	std::mt19937 gen(rd()); 
   	std::uniform_real_distribution<> dis(-1, 1);

	#pragma omp parallel for
	for (int t = 0; t < threads; ++t)
	{
		for (int i = 0; i < n4; ++i)
		{
			for (int j = 0; j < dim; ++j)
			{
				punkt(t,j) = dis(gen);
			}
			if (VectorXd(punkt.row(t)).squaredNorm() < 1){
				#pragma omp atomic
				NKreis++;
			}
			#pragma omp atomic
			punkteEcht++;
		}
	}

	if(punkteEcht != n){
		cout << "Es gab wohl probleme!" << endl;
		for (int i = 0; i < n - punkteEcht ; ++i)
		{
			for (int j = 0; j < dim; ++j)
			{
				punkt(0,j) = dis(gen);
			}
			if (VectorXd(punkt.row(0)).squaredNorm() < 1){
				NKreis++;
			}	
		}
	}

	if(dim == 3)
	return ((1.*NKreis)/(1.*n))*(3./4.)*pow(2,3);
	if(dim == 2)
	return ((1.*NKreis)/(1.*n))*pow(2,2);
	if(dim%2==1)
	return pow(tgamma(((1.*dim)/2.)+1)*((1.*NKreis)/(1.*n))*pow(2,dim),2./(1.*dim));
	if(dim%2==0)
	return pow(factorial(((1.*dim)/2.))*((1.*NKreis)/(1.*n))*pow(2,dim),2./(1.*dim));
}





























// //Der Lineare Kongruenz Zufallszahlengenerator
// VectorXd linkon(long long int r0, long long int a, long long int c, long long int m, long long int n);
// //Die Funktion sin^4 aus c)
// double p(double x);
// //Neumann verfahren mit übergabe von Zufallszahlen ACHTUNG! erzeugt nicht immer gleich viele soll es auch nicht s.u.
// VectorXd neumann(VectorXd zx, VectorXd zy, function < double(double)  > f);


// int main(){
// 	//Mit long long erzeugt man 64bit Ints

// 	//Paramatersatz für den Linearen Kongruenz Generator i)
// 	long long int r01,a1,c1,n1,m1;
// 	r01 = 1234;
// 	a1 = 20;
// 	c1 = 120;
// 	m1 = 6075;
// 	n1 = 1e5;

// 	//Paramatersatz für den Linearen Kongruenz Generator ii)
// 	long long int r02,a2,c2,n2,m2;
// 	r02 = 123456789;
// 	a2 = 65539;
// 	c2 = 0;
// 	m2 = 2147483648;
// 	n2 = 1e5;


// 	//Ziehe die Zufallzahlen aud dem Linearen Kongruenz Generator
// 	VectorXd r1;
// 	r1 = linkon(r01,a1,c1,m1,n1);
// 	speichere("Aufgabe1i",r1);

// 	VectorXd r2;
// 	r2 = linkon(r02,a2,c2,m2,n2);
// 	speichere("Aufgabe1ii",r2);

// 	//Implentiere den Mersenne Twister
// 	std::random_device rd;  
//     std::mt19937 gen(rd()); 
//     std::uniform_real_distribution<> dis(0, 1);

//     VectorXd rM(n2);

//     for (int i = 0; i < n2; ++i)
//     {
//     	rM[i] = dis(gen);			//dis(gen) erzeugt die Zufallszahlen durch den Mersenne Twister, speichere die im Vektor
//     }

//     speichere("Aufgabe1Mars",rM);

//     //#########################################################
//     //####################Aufgabe 1c###########################
//     //#########################################################

//     //Teile die Zufallszahlen in 2 Teile (einen für den x und einen für das y für Neumann)
//     VectorXd zx1,zy1,zx2,zy2,zx3,zy3;
//     zx1 = r1.segment(0,int(n1/2)); zy1 = r1.segment(int(n1/2),int(n1/2));
//     zx2 = r2.segment(0,int(n2/2)); zy2 = r2.segment(int(n2/2),int(n2/2));
//     zx3 = rM.segment(0,int(n2/2)); zy3 = rM.segment(int(n2/2),int(n2/2));

//     //Ziehe mit Neumann und speichern
//     VectorXd z1,z2,z3;
//     z1 = neumann(zx1,zy1,p);
//     z2 = neumann(zx2,zy2,p);
//     z3 = neumann(zx3,zy3,p);
//     speichere("Aufgabe1c1",z1);
//     speichere("Aufgabe1c2",z2);
//     speichere("Aufgabe1c3",z3);

//     //Hier kommt die Analytische Berechnung von p(x)
//     int n = 1e5;
//     VectorXd x(n),y(n);

//     for (int i = 0; i < n; ++i)
//     {
//     	x[i] = i*(1/(1.*n));	//Das hier lasst x von 0 bis 1 laufen
//     	y[i] = p(x[i]);			//Und berechnet y
//     }

//     speichere("Aufgabe1xAN",x);
//     speichere("Aufgabe1yAN",y);

//     //Die Spektraltests sind einfache Plots und daher in Python da alles mit 
//     //plt.plot(zufallszahl[:-1],zufallszahl[1:]) abgefrühstückt ist

// 	return 0;
// } 


// VectorXd linkon(long long int r0, long long int a, long long int c, long long int m, long long int n){
// 	VectorXd r(n);
// 	r[0] = r0;							//Setze Startwert
// 	for (int i = 0; i < n-1; ++i)
// 	{
// 		r[i+1] = fmod((a*r[i]+c),m);	//Berechne nach Vorschrift die nächsten Zahlenwerte
// 	}
// 	return r*(1./(1.*m));				//Teile durch m um floats zu erhalten
// }


// double p(double x){
// 	return (3./8.)*pow(sin(M_PI*x),4);	//Ganz simpel die Funktion aus c)
// }

// VectorXd neumann(VectorXd zx, VectorXd zy, function < double(double)  > f){
// 	long long int n = zx.size();

// 	VectorXd z(n);


// 	for (int i = 0; i < n; ++i)
// 	{
// 		if(zy[i] <= f(zx[i]))	//Liegt der Wert von y in f(x) nimm die Zufallszahl sonst NAN
// 		{						//Dass Neumann nicht n=1e5 Zufallszahlen erzeugt ist absicht damit
// 			z[i] = zx[i];		//Schlechte Zufallszahlgeneratoren auffällig schlechter sind
// 		}else{					//Bei 2 ist der Algorithmus der n Zahlen erzeugt drin :)
// 			z[i] = NAN;
// 		}
// 	}

// 	return z;
// }