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



int main()
{
	
	return 0;
}
























// //Box-Muller-Algorithmus
// VectorXd boxm(int n);
// //Zentraler Grenzwertsatz Funktion flexible sigma Eingabe
// VectorXd zgw(int n, double sigma);
// //Normaler Gauß halt
// double gaus(double x, double mu, double sigma);
// //Neumann Algorithmus um eine Gauß verteilung zu erhalten, dabei wird die Fläche angepasst siehe unten
// VectorXd neumannGAUS(int n, double mu, double sigma, double alpha, function < double(double,double,double)  > f);
// //Transformationsverfahren aus Kierfeld Sktipt
// VectorXd inversionp(int n);
// //Transformationsverfahren aus SMD
// VectorXd inversionpSMD(int n);
// //Die gesuchte Verteilung für d)
// double p(double x);

// int main(){
// 	int n = 1e5;	//Wir wollen immer 10^5 Zahlen

// 	VectorXd zg = boxm(n);		//Box-Muller für a)
// 	speichere("Aufgabe2a",zg);

// 	VectorXd Zzgw = zgw(n,1);	//Zentraler Grenzwertsatz für b)
// 	speichere("Aufgabe2b",Zzgw);

// 	VectorXd xAN(n),yAN(n);		//Anlytischer Gauß für Plot
// 	for (int i = 0; i < n; ++i)
// 	{
// 		xAN[i] = ((1.*i)/(1.*n))*10-5;	//Benutze den Gauß zwischen -5 und 5 (sollte hier reichen)
// 		yAN[i] = gaus(xAN[i],0,1);
// 	}
// 	speichere("Aufgabe2bxAN",xAN);
// 	speichere("Aufgabe2byAN",yAN);

// 	//Das hier ist nur zur Kontrolle gewesen mit einem alpha = 1.5 (wie im skript)
// 	VectorXd zNEU(n);
// 	double alpha = 1.5;
// 	zNEU = neumannGAUS(n,0,1,alpha,gaus);	
// 	speichere("Aufgabe2c",zNEU);

// 	//Berechne die Gaußverteilung mit Neumann für verschiedene alphas 
// 	MatrixXd zaNEU(6,n);
// 	VectorXd alphaV(6);

// 	alphaV << 0.01, 0.1, 1.0, 1.5, 10, 100;

// 	for (int i = 0; i < 6; ++i)
// 	{
// 		zaNEU.row(i) = neumannGAUS(n,0,1,alphaV[i],gaus);
// 	}
// 	speichere("Aufgabe2c2",zaNEU);



// 	VectorXd zinv(n);		//Inversionsmethode aus Skript
// 	zinv = inversionp(n);
// 	speichere("Aufgabe2d",zinv);

// 	//Analytische Verteilung aus d)
// 	VectorXd xANd(n),yANd(n);
// 	for (int i = 0; i < n; ++i)
// 	{
// 		xANd[i] = ((1.*i)/(1.*n))*20-10;	//Zwischen -5 und 5 Werte für p(x)
// 		yANd[i] = p(xANd[i]);
// 	}
// 	speichere("Aufgabe2dxAN",xANd);
// 	speichere("Aufgabe2dyAN",yANd);

// 	VectorXd zinvSMD(n);			//Andere Inversionsmethode
// 	zinvSMD = inversionpSMD(n);
// 	speichere("Aufgabe2dSMD",zinvSMD);

// 	return 0;
// }

// //Box-Muller Algorithmus
// VectorXd boxm(int n){
// 	VectorXd x1(n),x2(n),z1(n),z2(n);

// 	std::random_device rd;  
//     std::mt19937 gen(rd()); 
//     std::uniform_real_distribution<> dis(0, 1);

//     //Ziehe und speichere Zunächst 2 Zufallszahlen (jeweils n Stück)
//     for (int i = 0; i < n; ++i)
//     {
//     	x1[i] = dis(gen);
//     	x2[i] = dis(gen);
//     }
//     //Nutze Box-Muller zum Erzeugen der Gauß Verteilungen
//     for (int i = 0; i < n; ++i)
//     {
//     	z1[i] = sqrt(-2*log(x1[i]))*cos(2*M_PI*x2[i]);
//     	z2[i] = sqrt(-2*log(x1[i]))*sin(2*M_PI*x2[i]);
//     }
//     return z1;	//Übergebe einen der beiden (evtl wäre Mischen effizienter aber die Zahlen könnten korreliert sein...)
// }

// //Gauß mit Zentralen Grenzwertsatz
// VectorXd zgw(int n, double sigma){

// 	//Berechne zunächst wie viele Zahlen aufsummiert werden müssen
// 	//Das ist von Sigma abhängig und wird aufgerundet abgeschnitten für Int
// 	int N = int((pow(sigma,2)*12)+0.5); 

// 	//z sind später die Zufallszahlen und y hilfsvariable zum summieren
// 	VectorXd z(n);
// 	double y;

// 	std::random_device rd;  
//     std::mt19937 gen(rd()); 
//     std::uniform_real_distribution<> dis(0, 1);

//     //Gehe durch n (Anzahl an gewollten Zufallszahlen)
//     for (int i = 0; i < n; ++i)
//     {
//     	y = 0;							//Setze zunächst y = 0
//     	for (int j = 0; j < N; ++j)
//     	{
//     		y += dis(gen);				//Summiere N Zufallszahlen [0,1] auf
//     	}
//     	y -= N/2.;						//Ziehe N/2 ab (s.h. Skript)
//     	z[i] = y;						//Abspeichern
//     }
//     return z;							//Zurückgeben
// }

// //Üblicher Gauß mit x,mu,sigma als Input
// double gaus(double x, double mu, double sigma){
// 	return 1./(sqrt(2*M_PI)*sigma) * exp((-1.*pow((x-mu),2))/(2*pow(sigma,2)));
// }


// //Gauß mit Neumann Ziehen, hier kann theoretisch auch eine andere Funktion der form double(double,double,double)
// //erzeugt werden, jedoch ist zu beachten dass die Fläche in der Zahlen gezogen Werden sehr angepasst an gauß ist 
// VectorXd neumannGAUS(int n, double mu, double sigma, double alpha, function < double(double,double,double)  > f){
	
// 	std::random_device rd;
//     std::mt19937 gen(rd());
 	
//  	//Ziehe hier exponential verteilte Zahlen hier kommen zahlen der form 1 * exp(-1*x) heraus
//     std::exponential_distribution<> d(1);
// 	std::uniform_real_distribution<> dis(0, 1);
//     //z zum Abspeichern und zexp und zy als Hilfsvariablen
// 	VectorXd z(n);
// 	double zexp,zy;

// 	//i_n Zählt wie viele Zahlen wir durch Neumann haben damit wir auch 10^5 Stück bekommen
// 	int i_n = 0;
// 	do
// 	{
// 		//Ziehe zunächst einen Exponentialverteilung für x Wert
// 		//dabei wird eine weitere Zahl gezogen um auch den negativen Bereich abzudecken
// 		zexp =  (dis(gen) > 0.5)?(-1)*d(gen):d(gen); 
// 		//Ziehe in y Richtung gleichverteilte Zahlen gewichtet mit der Funktion die gefordert war um den Definitionsberech abzudeken
// 		zy = dis(gen)*alpha*exp(-fabs(zexp));
// 		//wenn die y Zahl in dem Gauß liegt akzeptieren sonst gehts weiter...
// 		if(zy <= f(zexp,mu,sigma))
// 		{
// 			z[i_n] = zexp;
// 			if(i_n < n){i_n++;}
// 		}
// 	}while(i_n < n);

// 	return z;
// }

// //Die gesuchte Verteilung
// double p(double x){
// 	return (1/M_PI)*(1/(1+pow(x,2)));
// }

// //Die inversen der gesuchten Verteilung
// double pi(double y){
// 	return sqrt((1/(M_PI*y))-1);
// }

// //Die Ableitung der inversen Verteilung
// double dpi(double y){
// 	return -1/(2*M_PI*pow(y,2)*sqrt((1/(M_PI*y))-1));
// }

// //Inversionsmethode aus dem Kierfeld Skript
// VectorXd inversionp(int n){
// 	std::random_device rd;  
//     std::mt19937 gen(rd()); 
//     std::uniform_real_distribution<> dis(0, 1);

//     VectorXd z(n),zN(n);
//     double y;

//     for (int i = 0; i < n; ++i)
//     {
//     	//Ziehe zufällig gleichverteilt im Intervall 0 bis 1/pi
//     	y = dis(gen)*((1/M_PI));
//     	zN[i] = (dis(gen) > 0.5)?((-1)*fabs(dpi(y))):(fabs(dpi(y)));
//     }
//     return zN;
// }

// //Inversionsmethode aus SMD
// VectorXd inversionpSMD(int n){
// 	std::random_device rd;  
//     std::mt19937 gen(rd()); 
//     std::uniform_real_distribution<> dis(0, 1);

//     VectorXd z(n),zN(n);
//     double y;

//     for (int i = 0; i < n; ++i)
//     {
//     	//Ziehe Zufallszahlen nach (z-1/2) also zwischen -1/2 und 1/2 und einsetzten -> s.h. getechte Formel
//     	y = dis(gen)-0.5;
//     	zN[i] = tan(M_PI*y);
//     }
//     return zN;
// }