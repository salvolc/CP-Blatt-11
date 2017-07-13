#ifndef FUNCTIONS
#define FUNCTIONS

#include <cmath>
#include <iostream>
#include <iomanip>
#include "../eigen/Eigen/Dense"
#include <functional>

using std::complex;
using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;
using std::cout;
using std::endl;


double lennard_jones(double& r);	
	



//Der Verlet Algorithmus bekommt Anfangsbedingungen, masse,
//mögliche Parameter für die Kräfte, eine Kraftfunktion (am besten eine die periodische RB hat da sonst alles keinen
//Sinn hat) und die periodischen RB selbst, um die Vektoren zu beschränken, dazu optional Thermostat Einstellungen
Eigen::MatrixXd* verlet_period(
		int steps,
		double h,
		Eigen::MatrixXd y0,
		Eigen::VectorXd m,
		Eigen::VectorXd params,
		std::function < Eigen::VectorXd(Eigen::MatrixXd,Eigen::VectorXd,int,Eigen::VectorXd)  > force,
		//Force Bekommt alle Teilchen diesmal auf jede Fall, eine ID für das Teilchen und massen und Parameter
		Eigen::VectorXd L,
		//Der L Vektor für die RB
		bool thermostat,
		double T
		//Thermostat Optionen (optional)
	);

//Diese Funktion soll einen Vektor r in die Schranken von L und 0 weisen also nix größer L und nix kleiner 0 mehr
//L ist die Randbedingung und r der zu prüfende Vektor
void test_period(Eigen::VectorXd& r,Eigen::VectorXd L); 

//Einfach die Kraft zwischen x und y berechnen wirkt auf x (hoffentlich?)
Eigen::VectorXd lennard_jones_kraft_naiv(Eigen::VectorXd x,Eigen::VectorXd y);
//Kraft auf ein Teilchen mit der particle_id wird berechnet, es werden periodische Randbedingungen angenommen
//welche mit Params übergeben werden können, die Masse wird aus konsistenz zu anderen Kraftfunktionen mitübergeben
Eigen::VectorXd lennard_jones_kraft_period(Eigen::MatrixXd y,Eigen::VectorXd m, int particle_id, Eigen::VectorXd params);


Eigen::MatrixXd schwerpunkts_geschwindigkeit(Eigen::MatrixXd* y,int steps);
Eigen::VectorXd E_Pot_period_par(Eigen::MatrixXd* y,int steps, Eigen::VectorXd L);

double E_kin_single(Eigen::VectorXd v,double m);
double E_kin_N_part_sum(Eigen::MatrixXd v_N, Eigen::VectorXd m);
Eigen::VectorXd E_kin_N_step(Eigen::MatrixXd* y,int steps, Eigen::VectorXd m);

Eigen::VectorXd temp_step(Eigen::MatrixXd* y,int steps,Eigen::VectorXd m);

Eigen::MatrixXd paarcorr(Eigen::MatrixXd *y, int steps, int steps_anf, int bins, Eigen::VectorXd L);

//Thermostat, gibt die angepassten Geschwindigkeiten für alle Teilchen zurück
Eigen::MatrixXd setze_temp(Eigen::MatrixXd y, double T);


//Initialisierungfunktionen, verändern den eingegebenen Vektor bzw. Matrix
void init_linspace_f(Eigen::VectorXd &r, double min, double max);
//Initialisiere Geschwindigkeit auf gewünschte Temperatur
void init_v0_T(Eigen::MatrixXd& v0,double T);
//Initialisiere Orte in einem Quadratischen Gitter 
void init_r0_period_2D_f(Eigen::MatrixXd& r0,Eigen::VectorXd L);



//Hilfsfunktionen und ein Ladebalken damit man weiß ob das Programm noch läuft
std::string toStrMaxDecimals(double value, int decimals);
void ladebalken(int i, int max);


//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//#####################################Lennard Jones Potential##################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################

//Simples Lennard Jones Abstand rein, Potential raus...
double lennard_jones(double& r){
	double val;
	val = 4*( (1/pow(r,12)) - (1/pow(r,6)) );
	return val;
}


//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//#####################################Lennard Jones Kraft######################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################


Eigen::VectorXd lennard_jones_kraft_naiv(Eigen::VectorXd x,Eigen::VectorXd y){
	double r = (x-y).norm();							//Berechne den Abstand der Teilchen
	Eigen::VectorXd val,r_vek;							//Erzeuge Ausgabe und Vektor der von x auf y zeigt
	//Falls der Abstand null ist üben die Teilchen keine Kraft aus.. das sollte eigentlich nicht vorkommen deshalb der 
	//cerr, ist vermutlich nicht das beste Handeling aber man merkt ob was nicht stimmt
	if(r == 0){std::cerr << "Die Teilchen sind an einem Ort!";val.resize(x.size());for (int i = 0; i < x.size(); ++i){val[i]=0;}return val;} 
	r_vek = x-y;										//Berechne den Vektor von x nach y
	val = 24*((2/pow(r,14)) - (1/pow(r,8)))*r_vek;		//Lennard Jones Kraft hoffentlich richtig
	return val;											//Gib den Kraftvektor zurück
}

//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//#####################################Verlet Algorithmus Sachen ###############################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################


//Der Verlet Algorithmus mit periodischen Randbedingungen
Eigen::MatrixXd* verlet_period(
		//Steps und h halt
		int steps,
		double h,
		//Anfangsbedingungen
		Eigen::MatrixXd y0,
		Eigen::VectorXd m,
		Eigen::VectorXd params,
		std::function < Eigen::VectorXd(Eigen::MatrixXd,Eigen::VectorXd,int,Eigen::VectorXd)  > force,
		Eigen::VectorXd L,
		bool thermostat = false,
		double T = 1.
		//Force Bekommt alle Teilchen diesmal auf jede Fall, eine ID für das Teilchen und massen und Parameter
	)
{
	if(y0.size()==0){std::cerr << "Startwerte unbekannt, Problem nicht bestimmbar!" << std::endl;exit(1);}

	//Wichtige oder nützliche Größen
	int space_dim = int(y0.row(0).size()/2);		//Anhand der Startwerte die Dimension und Teilchenzahl bestimmen
	int n_particles = y0.col(0).size();				//Die Startwerte sollten die Form y_start[n_particles][r,v] haben

	//Diese Vektoren helfen als Zwischenspeicher für die Verlet Berechnungen damit nicht immer y mit Rattenschwanz
	//geschrieben werden muss...
	Eigen::VectorXd a_pre(space_dim),a_post(space_dim),r_pre(space_dim),v_pre(space_dim),r_post(space_dim);

	//Diese Matrizen mit allen Beschleunigungen ist praktisch damit diese nicht mehrmals Berechnet werden muss
	//da die Kraftberechnung das aufwendigste ist
	Eigen::MatrixXd a_pre_all(n_particles,space_dim);
	Eigen::MatrixXd a_post_all(n_particles,space_dim);

	//Ausgabepointer bzw. Array mit allen Daten die wichtig sind
	Eigen::MatrixXd* y;
	y = new Eigen::MatrixXd[steps];

	//Gib den steps großen Matrix Array die Dimenisonen für die Matrizen
	//y hat letztlich die Gestalt von steps vielen Matrizen die Teilchenzahl viele Zeilen hat
	//in denen jeweils der Ort und die Geschwindigkeit steht
	for (int init = 0; init < steps; ++init)
	{
		y[init].resize(n_particles,space_dim*2);
	}

	//Anfangswerte setzten
	y[0] = y0;

	//Berechne die Beschleunigungen für den ersten schritt, so kann für i=0 bereits die Übergabe a_pre = a_post genutzt werden
	for (int particle_id = 0; particle_id < n_particles; ++particle_id)
	{
		a_post_all.row(particle_id) = force(y[0],m,particle_id,L)/m[particle_id];	
	}

	for (int i = 0; i < steps-1; ++i)
	{
		//Für den Verlet Algorithmus werden drei Rechenoperationen durchgefürht
		//Zunächst wird r[n+1] bestimmt, aus diesem wird die Kraft bestimmt
		//a[n+1] welches die force(r[n+1],t[n+1]) ist wobei t hier irrelavant ist
		//zuletzt wird v[n+1] berechnet Dies geschieht für alle Teilchen in einem Schritt
		//Doch alles nochmal der Reihe nach
		a_pre_all = a_post_all;	//Die jetztigen Beschleunigungen sind die i+1 Beschleunigung der letzten Iteration
		for (int particle_id = 0; particle_id < n_particles; ++particle_id)
		{
			a_pre = a_pre_all.row(particle_id);							//a_n für das Teilchen	
			r_pre = y[i].row(particle_id).segment(0,space_dim);			//r_n für das Teilchen
			v_pre = y[i].row(particle_id).segment(space_dim,space_dim);	//v_n für das Teilchen

			r_post = r_pre + v_pre*h + 0.5 * pow(h,2) * a_pre;			//Bestimme r_n+1

			//Hier kommen Randbedingungen ins Spiel mit dieser Funktion wird r wieder in Kasten gebracht
			test_period(r_post,L);										

			//Ist sichergestellt dass r noch im Kasten ist wird der Teilchenort in y gespeichert
			y[i+1].row(particle_id).segment(0,space_dim) = r_post;
		}
		//Der Teilchenort ist zunächst für alle Teilchen berechnet worden da hier bei a_post...
		for (int particle_id = 0; particle_id < n_particles; ++particle_id)
		{
			//.. y[i+1] genutzt wird um die neuen Beschleunigungen zu berechnen (hier sollte noch durch m geteilt werden)
			a_post = force(y[i+1],m,particle_id,L)/m[particle_id];
			a_post_all.row(particle_id) = a_post;

			//a_pre ist nicht mehr das für das Teilchen da die Schleife oben das a_pre vom letzten Teilchen drin hatte
			//somit muss unser in weiser vorraussicht angelegter Speicher für die Beschleunigungen genutzt werden
			//damit sparen wir eine mühsame Kraftrechnung an dieser Stelle
			a_pre = a_pre_all.row(particle_id);

			//Fast geschaft wegen der oberen Schleife ist v_n auch noch das v_n des allerletzten Teilchen
			//das ist schlecht also holen wir und das v_n mit unserem y[i]
			v_pre = y[i].row(particle_id).segment(space_dim,space_dim);


			//Jetzt muss nur noch die Geschwindigkeit mit dem Verlet Algorithmus berechnet werden und das nochmal
			//für alle steps und das wars auch schon
			v_pre = v_pre+0.5*(a_pre+a_post)*h;
			y[i+1].row(particle_id).segment(space_dim,space_dim) = v_pre;
		}
		//Thermostat fall notwendig, passt die Geschwindigkeiten an die Temperatur an..
		if(thermostat)y[i+1].block(0,space_dim,n_particles,space_dim) = setze_temp(y[i+1],T);
		ladebalken(i,steps-1);
	}
	return y;	//Puh das war anstrengend aber jetzt kann der zeiger auf das Array wieder zur Main zurück
}


//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//#############################Dinge für periodische Randbedingungen############################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################



//Hier kommt der Test ob der Vektor L in dem Kasten liegt
void test_period(Eigen::VectorXd& r,Eigen::VectorXd L){ //L ist die Randbedingung und r der zu prüfende Vektor
	int dim = r.size();
	for (int i = 0; i < dim; ++i) //Teste jede Dimension einzeln
	{
		if(r[i] >= L[i] || r[i] <= 0){
			r[i] = r[i] - L[i] * floor(r[i]/L[i]);	//Aus dem Skript das Zurücksetzten der Teilchen
		}
	}
	return;
}

//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//###############################Die Kraftfunktion periodisch###################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################

Eigen::VectorXd lennard_jones_kraft_period(Eigen::MatrixXd y,Eigen::VectorXd m, int particle_id, Eigen::VectorXd params){
	Eigen::VectorXd L; 												//Das sind die Begrenzungen für die Teilchen	
	L = params;														//L = {L,L} ist wie in der Aufgabe so L/2
	int n_particles = y.col(0).size();								//Teilchen und Dimension wie auch sonst immer über y
	int space_dim = int(y.row(0).size()/2.);						
	double r_max = (L[0]+L[1])*0.25;								//r_max für gleiche Begrenzungen in den Richtungen 
	Eigen::MatrixXd r_mat = y.block(0,0,n_particles,space_dim);		//Die Orte aller Teilchen in r_mat damit sind die v's weg

	Eigen::VectorXd kraft(space_dim),r_particle(space_dim);			//Die letztliche Kraft auf ein Teilchen r_particle
	r_particle = r_mat.row(particle_id);							//Speichere das Teilchen um das es geht in r_particle
	for (int init = 0; init < space_dim; ++init)					//ist sehr nützlich
	{
		kraft[init]=0;												//Initialisiere Kraft damit wenn nix gerechnet wird  
	}																//auch null rauskommt

	double dist;
	//Hilfsvektoren r_ref ist das Teilchen zu dem die Kraft Bestimmt wird, r_ref_bild findet das nächste Teilchen bzw.
	//Bildteilchen (s.h. unten)
	Eigen::VectorXd r_diff(space_dim),r_ref(space_dim),r_ref_bild(space_dim);

	for (int i_particle = 0; i_particle < n_particles; ++i_particle)
	{
		if(particle_id == i_particle)continue;				//Teilchen wirkt nicht aus sich selbst
		r_ref = r_mat.row(i_particle);						//Speichere den Ort des Teilches zu dem die Kraft bestimmt wird
		r_diff = r_ref-r_particle;							//Zünachst Diffrenzvektor zum nicht Bildteilchen
		r_ref_bild = r_ref;									
		for (int i_dim = 0; i_dim < space_dim; ++i_dim)
		{
			//Falls es ein Bildteilchen gibt welches näher an unserem Teilchen ist dann berechne die Kraft zu
			//dem Bildteilchen, hier wird also der Abstand in den einzelnen Komponenten minimiert indem 
			//Bildteilchen berücksichtigt werden, somit erhält man das nächste (Bild)Teilchen zu unserem Teilchen
			//Es ist nur wichtig das nächste zu betrachten weil es nicht 2 geben kann die in L/2 radius liegen
			if(sqrt(pow(r_diff[i_dim],2)) > sqrt(pow(r_ref[i_dim]-r_particle[i_dim]+L[i_dim],2)) ) 
				r_ref_bild[i_dim]=r_ref[i_dim]+L[i_dim];
			if(sqrt(pow(r_diff[i_dim],2)) > sqrt(pow(r_ref[i_dim]-r_particle[i_dim]-L[i_dim],2)) ) 
				r_ref_bild[i_dim]=r_ref[i_dim]-L[i_dim];
		}
		r_diff = r_ref_bild-r_particle;		//Abstandsvektor zum nächsten (Bild)Teilchen und unserem Teilchen
		dist = r_diff.norm();
		//Ist der Abstand im CutOff drin dann addiere den Beitrag der Kraft auf
		if(dist <= r_max)kraft += lennard_jones_kraft_naiv(r_particle,r_ref_bild);
	}
	return kraft;													
}

//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//###############################Wichtige Größen des Systems####################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################

//Berechne den Schwerpunktsgeschwindigkeitsvektor für alle steps
Eigen::MatrixXd schwerpunkts_geschwindigkeit(Eigen::MatrixXd* y,int steps){
	int n_particles = y[0].col(0).size();			//Teilchenzahl und Dimension wie oben auch
	int space_dim 	= int(y[0].row(0).size()/2);

	Eigen::VectorXd v_sum(space_dim),v_i(space_dim);//in V_sum kommt pro step die 
	Eigen::MatrixXd v_s(steps,space_dim);			//Hier kommt die Schwerpunktsgeschwindigkeit rein
	std::cout << "\nBerechne die Schwerpunkts Geschwindigkeit:" << std::endl;
	for (int i = 0; i < steps; ++i)
	{	
		for (int init = 0; init < space_dim; ++init)
		{
			v_sum[init]=0;
		}
		for (int particle_id = 0; particle_id < n_particles; ++particle_id)
		{
			v_i = y[i].row(particle_id).segment(space_dim,space_dim);		//in y sind die Geschwindigkeiten in den Spalten
			v_sum += v_i;													//Summiere die Beiträge aller Teilchen
		}
		v_s.row(i) = v_sum;													//Speichere die Summe 
		ladebalken(i,steps);
	}
	std::cout << "Die Schwerpunkts Geschwindigkeit ist fertig berechnet!" << std::endl;
	return v_s;
}


Eigen::VectorXd E_Pot_period_par(Eigen::MatrixXd* y,int steps, Eigen::VectorXd L){
	int n_particles = y[0].col(0).size();							//L = {L,L} ist wie in der Aufgabe so L/2
	int space_dim = int(y[0].row(0).size()/2.);						//Teilchen und Dimension wie auch sonst immer über y
	double r_max = (L[0]+L[1])*0.25;								//r_max für gleiche Begrenzungen in den Richtungen 
	Eigen::MatrixXd r_particle(steps,space_dim);					

	Eigen::VectorXd dist(steps),pot(steps),pot_v(steps);	//In pot wird für jedes Teilchen aufsummiert
															//das ergebnis kommt in pot_v und wird zurückgegeben
	for (int init = 0; init < steps; ++init)
	{
		pot_v[init] = 0;
	}

	//Hilfs"vektoren" für die einzelnen Steps
	Eigen::MatrixXd r_diff(steps,space_dim),r_ref(steps,space_dim),r_ref_bild(steps,space_dim);

	std::cout << "\nBerechne die potentielle Energie:" << std::endl;
	//Durch das Multithreading ist die dimension der Hilfvektoren höher damit gleichzeite Threads sich diese nicht 
	//gegenseitig überschreiben..
	#pragma omp parallel for
	for (int i = 0; i < steps; ++i)
	{
		pot[i] = 0;
		//Identisch zur Kraft Berechnung..
		for (int j_particle = 0; j_particle < n_particles; ++j_particle)
		{
			for (int i_particle = j_particle; i_particle < n_particles; ++i_particle)
			{
				if(j_particle == i_particle)continue;
				r_ref.row(i) = Eigen::VectorXd(y[i].row(i_particle)).segment(0,space_dim);
				r_particle.row(i) = Eigen::VectorXd(y[i].row(j_particle)).segment(0,space_dim);

				r_diff.row(i) = r_ref.row(i)-r_particle.row(i);
				r_ref_bild.row(i) = r_ref.row(i);
				for (int i_dim = 0; i_dim < space_dim; ++i_dim)
				{
					if(sqrt(pow(r_diff(i,i_dim),2)) > sqrt(pow(r_ref(i,i_dim)-r_particle(i,i_dim)+L[i_dim],2)) ) 
						r_ref_bild(i,i_dim)=r_ref(i,i_dim)+L[i_dim];
					if(sqrt(pow(r_diff(i,i_dim),2)) > sqrt(pow(r_ref(i,i_dim)-r_particle(i,i_dim)-L[i_dim],2)) ) 
						r_ref_bild(i,i_dim)=r_ref(i,i_dim)-L[i_dim];
				}
				r_diff.row(i) = r_ref_bild.row(i)-r_particle.row(i);
				dist[i] = r_diff.row(i).norm();
				//Einziger Unterschies es wird die Potentielle Energie Gespeichert
				if(dist[i] <= r_max)pot[i] += lennard_jones(dist[i]);
			}
		}
		pot_v[i] = (pot[i]- lennard_jones(r_max));
		#pragma omp critical	//damit die Ausgabe nicht kaputt geht
		ladebalken(i,steps);

	}
	std::cout << "Die potentielle Energie ist fertig berechnet!" << std::endl;
	return pot_v;			 										
}


//Simples Ekin ausrechnen mit v und m sollte klar sein
double E_kin_single(Eigen::VectorXd v,double m){
	double ekin,vbetrag;
	vbetrag = v.norm();
	ekin = 0.5*m*pow(vbetrag,2);
	return ekin;
}

//Berechne Ekin für viele Teilchen (aufsummiert)
double E_kin_N_part_sum(Eigen::MatrixXd v_N, Eigen::VectorXd m){
	Eigen::VectorXd v;
	int n_particles = v_N.col(0).size();	//Wie immer
	double E_Kin = 0;						//Hier drin wird die energie pro Teilchen aufsummiert
	for (int i_particle = 0; i_particle < n_particles; ++i_particle)
	{
		v = v_N.row(i_particle);
		E_Kin += E_kin_single(v,m[i_particle]); //Addiere die Einzelbeiträge
	}
	return E_Kin;
}

//Berechne Ekin für alle Steps für alle Teilchen gib für jeden Step eine Energie in einem Vektor zurück
Eigen::VectorXd E_kin_N_step(Eigen::MatrixXd* y,int steps, Eigen::VectorXd m){
	int n_particles = y[0].col(0).size();
	int space_dim = int(y[0].row(0).size()/2);
	Eigen::VectorXd E_Kin(steps);
	Eigen::MatrixXd v;
	for (int i = 0; i < steps; ++i)
	{
		v = y[i].block(0,space_dim,n_particles,space_dim);
		E_Kin[i] = E_kin_N_part_sum(v,m);	
	}
	return E_Kin;
}

//Berechne die Temperatur pro Step
Eigen::VectorXd temp_step(Eigen::MatrixXd* y,int steps,Eigen::VectorXd m){
	int n_particles = y[0].col(0).size();
	int dim = int(y[0].row(0).size()/2);
	double E_Kin;
	double kb = 1.;
	double Nf = dim*n_particles-dim;	//Anzahl der Freiheitsgrade wie im Skript
	Eigen::VectorXd T(steps);
	for (int i = 0; i < steps; ++i)
	{
		E_Kin = E_kin_N_part_sum(y[i].block(0,dim,n_particles,dim),m);
		T[i] = (2*E_Kin)/(kb*Nf);		//Formal für T aus Kierfeld Skript 
	}
	return T;
}

//Die Paarcorrleation..., hier wird eine Matrix ausgeben die nur zwei zeilen hat
//Eine mit den Binmitten für die Abstände und eine in der das gemittelte Auftrittsvorkommen der Teilchen steht
//Das steps_anf steht für den Schritt ab dem Gerechnet werden also ab wann ist das alles äquilibriert
Eigen::MatrixXd paarcorr(Eigen::MatrixXd *y, int steps, int steps_anf, int bins, Eigen::VectorXd L)
{
	double r_max = L[0]*0.5;						//Es ist nur sinvoll bis L/2 zu rechenen s.h. Kierfeld Skript

	int n_particles = y[0].col(0).size();			//Wie auch schon oben
	int space_dim = int(y[0].row(0).size()/2.);

	Eigen::MatrixXd paare(2,bins);					//Für die Ausgabe
	Eigen::VectorXd bingrenzen(bins+1),binmitten(bins),paar(bins);
	double binbreite = r_max/(1.*bins);				//Die Binbreite um die Binnummer für einen Abstanf zu finden

	init_linspace_f(bingrenzen,0,r_max);			//Initialisiere die Bingrenzen

	for (int i_bin = 0; i_bin < bins; ++i_bin)		//und berechne dann die Binmitten die um einen kleiner sind
	{
		binmitten[i_bin] = (bingrenzen[i_bin]+bingrenzen[i_bin+1])/2.;
	}
	paare.row(0) = binmitten;						//Wie oben in die erste Zeile kommen die Binmitten zum plotten

	Eigen::VectorXd r_diff,r_ref,r_ref_bild,r_particle;	//Hilfsvektoren
	double dist;

	int binnummer,n_steps,y_step;		//Dieses y_step ist damit der richtige step ausgewählt wird also nach dem 
	n_steps = steps - steps_anf;		//äquilibrieren bis zum schluss

	//Der erste teil hier ist wieder gleich der Kraft und der Potentiellen Energie
	for (int i = 0; i < n_steps; ++i)
	{
		for (int i_particle = 0; i_particle < n_particles; ++i_particle)
		{
			for (int j_particle = i_particle; j_particle < n_particles; ++j_particle)
			{
				y_step = steps_anf + i;
				r_ref = y[y_step].row(i_particle).segment(0,space_dim);
				r_particle = y[y_step].row(j_particle).segment(0,space_dim);

				r_diff = r_ref-r_particle;
				r_ref_bild = r_ref;

				for (int i_dim = 0; i_dim < space_dim; ++i_dim)
				{
					if(sqrt(pow(r_diff[i_dim],2)) > sqrt(pow(r_ref[i_dim]-r_particle[i_dim]+L[i_dim],2)) ) 
						r_ref_bild[i_dim]=r_ref[i_dim]+L[i_dim];
					if(sqrt(pow(r_diff[i_dim],2)) > sqrt(pow(r_ref[i_dim]-r_particle[i_dim]-L[i_dim],2)) ) 
						r_ref_bild[i_dim]=r_ref[i_dim]-L[i_dim];
				}
				r_diff = r_ref_bild-r_particle;
				dist = r_diff.norm();
				//Hier wird es nochmal interessant die Binnummer ergibt sich durch aufrunden nachdem der Abstand
				//Durch die Binbreite geteilt worden ist, dies lässt sich durch überlegen und testen herausfinden
				binnummer = int((dist/binbreite)+0.5);
				//Diese letzten Kritieren Prüfen ob der Abstand kleiner war als r_max 
				if(binnummer < bins && binnummer > 0){paar[binnummer] += 1;}
			}
		}
		ladebalken(i,n_steps);
	}

	Eigen::VectorXd flaeche(bins);
	for (int i_bin = 0; i_bin < bins; ++i_bin)
	{
		//Fläche der Kreisringe mit denen die Bins gewichtet werden müssen
		flaeche[i_bin] = M_PI*(pow(bingrenzen[i_bin+1],2)-pow(bingrenzen[i_bin],2)); 
		//Hier geschieht der letzte Schritt mit der Korrektur der Fläche, der Dichte und der Teilchenzahl
		//Der Faktor 2 bei n_particles/2 sorgt dafür dass die anzahl an Einträgen in den Bins verdoppelt wird
		//Dies hat den Grund, dass in der Berechnung im Skipt das Paar i,j doppelt gezählt wird
		//Da dies nicht effizient wäre wird hier nur die nötige Hälfte gerechnet und dann mit 2 multipliziert
		paar[i_bin] = paar[i_bin]/(n_steps*flaeche[i_bin]*(n_particles/2)*(n_particles/(L[0]*L[1]))); //(flaeche[r_i]*n_steps*n_particles*(n_particles/(L[0]*L[1])));
	}
	paare.row(1) = paar;


	return paare;
}

//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//###############################Initialisierungs Funktionen####################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################

//Initialisiere die Anfangsgeschwindigkeit
void init_v0_T(Eigen::MatrixXd& v0,double T){
	int n_particles = v0.col(0).size();		//Dinge die weiter oben schon geklärt sind
	int dim = v0.row(0).size();
	double kb = 1;
	double Nf = dim*n_particles-dim;
	double vbetrag=0;

	Eigen::VectorXd schwer_v0(dim);	

	for (int init_dim = 0; init_dim < dim; ++init_dim)
	{
		schwer_v0[init_dim]=0;
	}
	for (int init_p = 0; init_p < n_particles; ++init_p)
	{	
		for (int init_dim = 0; init_dim < dim; ++init_dim)
		{
			v0(init_p,init_dim)=((rand()%21)/10.-1);
		} 
		//Berechne die Schwepunktsgeschwindigkeit
		schwer_v0 += (1./(1.*n_particles))*Eigen::VectorXd(v0.row(init_p));
	}
	for (int init = 0; init < n_particles; ++init)
	{
		//Ziehe die Schwepunktsgeschwindigkeit ab und berechne die Quadartische norm die nötig ist um die Temp zu skalieren
		v0.row(init) = (Eigen::VectorXd(v0.row(init))-schwer_v0);
		vbetrag += v0.row(init).squaredNorm();
	}
	//Skalierung für die Temperatur die Formel stammt umgestellt aus dem Kierfeld Skript und muss an alle Geschwindigkeiten dran
	double t_weight = sqrt((T*Nf*kb)/(vbetrag));
	for (int init = 0; init < n_particles; ++init)
	{
		v0.row(init) = Eigen::VectorXd(v0.row(init))*t_weight;
	}
	return;
}

//Diese Funktion erzeugt einen Vektor mit gleichen Abständen von min bis max
void init_linspace_f(Eigen::VectorXd &r, double min, double max){
	int size = r.size();
	double inter = max - min;
	for (int i = 0; i < size; ++i)
	{
		r[i] = min + (inter/(size-1)) * i;
	}
	return;
}

//Anfangsbedingungen für die Orte
void init_r0_period_2D_f(Eigen::MatrixXd& r0,Eigen::VectorXd L){
	int n_particles = r0.col(0).size();
	int dim 		= r0.row(0).size();
	double Lx = L[0];
	double Ly = L[1]; 
	int rows = 1;
	while(pow(rows,2) < n_particles){rows++;}	//Ein Quadartisches Gitter bis alle Punkte draufpassen
	Eigen::VectorXd x(rows),y(rows);
	init_linspace_f(x,1,Lx); init_linspace_f(y,1,Ly);
	double i_x=0;double i_y=0;
	for (int i = 0; i < n_particles; ++i)
	{

		r0(i,0) = x[i_x];
		r0(i,1) = y[i_y];
		i_y++;
		if(i_y > y.size()-1){i_x++;i_y=0;}
	}
	return;
}

//Setzte die Temperatur, dies ist ein Isokinetisches Thermostat
Eigen::MatrixXd setze_temp(Eigen::MatrixXd y, double T){
	Eigen::MatrixXd v;
	int n_particles = y.col(0).size();
	int dim = int(y.row(0).size()/2);
	v = y.block(0,dim,n_particles,dim);

	double kb = 1;
	double Nf = dim*n_particles-dim;
	double vbetrag = 0;

	for (int i_particle = 0; i_particle < n_particles; ++i_particle)
	{	
		vbetrag += v.row(i_particle).squaredNorm();
	}
	//Der Gewichtungsfaktor ist äquivalent zu dem bei der Initialisierung
	double t_weight = sqrt(T*Nf*kb/vbetrag);
	v *= t_weight;

	return v;
}




//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//###########################################Hilfsfunktionen####################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################

//Zahl und nachkommastellen rein und string kommt raus
std::string toStrMaxDecimals(double value, int decimals)
{
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(decimals) << value;
    std::string s = ss.str();
    if(decimals > 0 && s[s.find_last_not_of('0')] == '.') {
        s.erase(s.size() - decimals + 1);
    }
    return s;
}

//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//#############################Ein kleiner Ladebalken###########################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################
//##############################################################################################################

//Ganz Praktisch
void ladebalken(int i, int max){
	double progress = (1.*i)/(1.*max)*100;
	std::cout << "\rSchritt " << i << " von " << max << " geschafft! " << "Also " << progress << " %";
	if(i == max-1 || i==max){std::cout << "\rSchritt " << max << " von " << max << " geschafft. Fertig!" << std::endl;}
	
	return;
}



int delta(int n,int m){if(n == m) return 1; else return 0;}     //Implementierung der benötigten deltafunktiom delta_nm
int theta(double x){if(x < 0) return 0; else return 1;}			//Implementierung der benötigten Thetafunktion theta(b/2 - |xi_n|) 


//################################################################################################################################
//##################################### Funktion für Komponenten des diskr. Hamilton  ############################################
//################################################################################################################################

complex<double> h_nm(int n,int m,double V0, double b, double delXi, double Ximin)  
{
	complex<double> H; //H_nm als komplexe Einträge für H
	double Hr;         //Realteil von H_nm  
 
	Hr = -(1/(pow(delXi,2)))*(delta(n,m+1)+delta(n,m-1)-2.*delta(n,m)) + V0*theta((b/2.)-fabs(Ximin+n*delXi))*delta(n,m); //Berechnung von Realteil von H_nm
																														  //mit auf dem Zettel angegebener Formel
	H.real(Hr);        //Lösung wird in den Realteil von H_nm geschrieben
	return H;		   //komplexes Matrixelement H_nm wird zurückgegeben
}

//################################################################################################################################
//####################################### Funktion für Matrix des diskr. Hamilton  ###############################################
//################################################################################################################################

MatrixXcd Ham(double maxXi, double minXi, double delXi, double V0, double b){
	int n = int((maxXi-minXi)/delXi);                    // Dimension der Matrix für den Hamilton bestimmen
	int m = n;						                     // Hamilton ist quadratisch

	MatrixXcd H(n,m);
	for (int i_n = 0; i_n < n; ++i_n)
	{
		for (int i_m = 0; i_m < m; ++i_m)
		{
			H(i_n,i_m) = h_nm(i_n,i_m,V0,b,delXi,minXi); //Matrix für alle n und m aus nxm mit den aus h_nm() bestimmten
		}                                                //Elementen füllen.
	}
	return H;                                            // komplexe Matrix des Hamiltonoperators wird zurückgegbeben
}


//################################################################################################################################
//############################## Funktion für Matrix des diskr. Zeitentwicklungsoperators  #######################################
//################################################################################################################################


MatrixXcd ZeitEnt(double delT, MatrixXcd H){
	MatrixXcd ZeitEnt;  				  // Initialisierung des Zeitentwicklungsoperator als komplexe Matrix
	int n = H.col(0).size();   			  // Dimension des Hamilton bestimmen
	int m = H.row(0).size();        

	MatrixXd E = MatrixXd::Identity(n,m); // Erzeugung einer Einheitsmatrix der Dimension nxm (Dimension des Hamiton)

	ZeitEnt = (E+complex<double>(0,0.5)*H*delT).inverse() * (E-complex<double>(0,0.5)*H*delT); //Bestimmung des Zeitentwicklungsoperators mit
 																							   //auf Zettel angegebener Formel
	return ZeitEnt;						  // komplexe Matrix des Zeitentwicklungsoperators wird zurückgegbeben	
}


//################################################################################################################################
//############################ Funktion für komplexen diskr. Zustandsvektor (normiertes Gauß-Paket)  #############################
//################################################################################################################################

VectorXcd gaus_P_dis(double maxXi, double minXi, double delXi, double Xi0, double sigma, double k){
	int n = int((maxXi-minXi)/delXi);     // Dimension für den Zustandsvektor bestimmen
	VectorXcd gaus(n);					  // Initialisierung des Zustandsvektors der Dimension nx1	

	for (int i_n = 0; i_n < n; ++i_n)
	{
		gaus(i_n) = pow((1/(2*M_PI*sigma)),0.25)*exp(-pow(((minXi+delXi*i_n)-Xi0),2) /(4*sigma))*exp(complex<double>(0,1)*k*(minXi+i_n*delXi)); 
		                                  // Berechnung des Zustandsvektors mit auf Zettel angegebener Formel, jedoch diskretisiert
	}
	return gaus;                          // komplexer Vektor für Zustandsvektor wird zurückgegeben
}


//################################################################################################################################
//###################################### Zeitentwicklung eines Zustandsvektors  ##################################################
//################################################################################################################################


MatrixXcd schroed(MatrixXcd Sh, VectorXcd psi0, int steps){

	int dim = psi0.size();    
	MatrixXcd psi(steps,dim);                     //Matrix für alle Zustandsvektoren in Zeitintervallen deltatau als Zeilen
	psi.row(0) = psi0;					          //Anfangszustand ist erste Zeile der Matrix
	for (int i = 1; i < steps; ++i)

	{
		psi.row(i) = Sh*VectorXcd(psi.row(i-1));  //Wiederholte Anwendung des Zeitentwicklungsoperators auf Anfangszustand, 
	}											  //jeder Zwischenschritt wird in MatrixZeile geschrieben.
												  //Steps Iterationen, bis Zeit tau = steps * deltatau erreicht		

	return psi;
}


//################################################################################################################################
//###################################### Bestimmung der Transmissionskoeffizienten  ##############################################
//################################################################################################################################


VectorXd trans(MatrixXcd psi, double limXi, double delXi, double minXi){ // Übergebe Matrix mit Zustandsvektoren für Zeitschritte
	int tn = psi.col(0).size();				//Anzahl der Zeitschritte als Größe der Zeilen
	int n = psi.row(0).size();				//Dimension der 1xn Zustandsvektoren
	VectorXd trans = VectorXd::Zero(tn);	//Initialisierung des Vektors für Transmissionskoeffizienten für t_n Zeitschritte

	for (int i_t = 0; i_t < tn; ++i_t)		//Iteration durch Zeitschritte
	{
		for (int i_n = 0; i_n < n; ++i_n)	//Iteration durch Ortschritte
		{
			if(minXi+delXi*i_n > limXi)		//Bedingung Zettel xi > 0
				trans[i_t] += delXi*norm(psi(i_t,i_n));  //Berechnung Transmissionskoeffizienten mit Formel von Zettel
		}
	}
	return trans;
}


//################################################################################################################################
//################################################# Gauß-Seidel Iteration  #######################################################
//################################################################################################################################


//####################################RHO ist Platzhalter#########################################
MatrixXd gaus_seidel(MatrixXd phi, MatrixXd rho, double del){ // Übergabe Anfangsbedingungen, Ladungsdichte, Schrittweite
	int j_max = phi.col(0).size();							  // Maximales j für x Iteration
	int l_max = phi.row(0).size();							  // Maximales l für y Iteration	 
	MatrixXd phi_post=phi;                                    // Initialisierung diskr. Phi als Matrix der Dimension j_max x l_max
															  // Anfangsbedingungen zunächst in Phi-Matrix schreiben
	for (int j = 1; j < j_max-1; ++j)
	{
		for (int l = 1; l < l_max-1; ++l)
		{
			phi_post(j,l) = 0.25*(phi(j+1,l)+phi_post(j-1,l)+phi(j,l+1)+phi_post(j,l-1))+0.25*pow(del,2)*rho(j,l);
			//Iteration durch alle l, j und berechnung von Phi mit Gau-Seidel Iteration, Rand wird ausgespart
		}
	}
	return phi_post;
}


//################################################################################################################################
//#################################################### E-Feld Berechnung #########################################################
//################################################################################################################################


MatrixXd E_Feld(MatrixXd phi,double del){
	int j_max = phi.col(0).size(); 
	int l_max = phi.row(0).size();
	MatrixXd E_Feld(j_max,l_max);           //E-Feld Dimension ist gleich der Dimension von Phi
	for (int j = 1; j < j_max-1; ++j)
	{
		for (int l = 1; l < l_max-1; ++l)
		{
			E_Feld(j,l) = sqrt(pow(((phi(j+1,l)-phi(j-1,l))/(2*del)),2)+pow(((phi(j,l+1)-phi(j,l-1))/(2*del)),2));
			//Ableitung des Potentials mit symmetrischem Differenzenquotient für x und y und Betrag des E-Feldvektors aus sqrt(x^2+y^2)
		}
	}
	return E_Feld;

}

//################################################################################################################################
//################################################ Analytische Lösung für Psi ####################################################
//################################################################################################################################


double ana_d(double x,double y,int maxit){
	double erg_pre=0,erg_post=0; //Erst Werte auf Null setzen damit aufaddiert werden kann
	double n = 1;
	do
	{
		erg_pre = erg_post;
		erg_post = erg_pre + ( (2*(1-cos(n*M_PI))) / (n*M_PI*sinh(n*M_PI)) ) * sin(n*M_PI*x) * sinh(n*M_PI*y);
		// Summe über n für die Formel der analytischen Lösung von dem Zettel für ein übergebenes x und y
		n++;
	}
	while(n < maxit); // Auf jeden Fall mehr als 10 Iterationen
	return erg_post;
}



//################################################################################################################################
//####################################### diskr. Analytische Lösung für Psi als Matrix ###########################################
//################################################################################################################################


MatrixXd ana_dis(VectorXd L, double del) //Systemgröße und Schrittweite übergeben
{
	int j_max,l_max;		
	j_max = int(L[0]/del)+1; // Maximales j für x-Iteration
	l_max = int(L[1]/del)+1; // Maximales l für y-Iteration

	MatrixXd ana_dis(j_max,l_max); //Initialisierung Matrix für diskr. analytische Lösung, gleiche Dimension wie Phi_numerisch

	for (int j = 0; j < j_max; ++j)
	{
		for (int l = 0; l < l_max; ++l)
		{
			ana_dis(j,l) = ana_d(j*del,l*del,200);
			//Iteration durch alle j,l, bzw. x,y und berechnen der Matrixelemente mit ana_d für bestimmten x bzw y
		}
	}
	return ana_dis;
}

double cube(double r, double x){
	return r*x-pow(x,3);
}
VectorXd cube_it(double r, int steps, double x0){
	VectorXd val = VectorXd::Zero(steps);
	val[0] = x0;
	for (int i = 0; i < steps-1; ++i)
	{
		val[i+1] = cube(r,val[i]);
	}
	return val;
}
MatrixXd cube_it_M(double r, int steps, int stepsX){
	MatrixXd val = MatrixXd::Zero(stepsX,steps);
	VectorXd x0 = VectorXd::Zero(stepsX); 
	//init_linspace_f(x0,-sqrt(1+r),sqrt(1+r));
	//x0[int(x0.size()/2)] = 0;
	init_random_one_d(x0,-sqrt(1+r),sqrt(1+r));

	for (int i_x = 0; i_x < stepsX; ++i_x)
	{
		val.row(i_x) = cube_it(r,steps,x0[i_x]);
	}
	return val;
}


double logist(double r, double x){
	return r*x*(1-x);
}

VectorXd logist_it(double r, int steps, double x0){
	VectorXd val = VectorXd::Zero(steps);
	val[0] = x0;
	for (int i = 0; i < steps-1; ++i)
	{
		val[i+1] = logist(r,val[i]);
	}
	return val;
}

MatrixXd logist_itM(double r, int steps, int stepsX){
	MatrixXd val = MatrixXd::Zero(stepsX,steps);
	VectorXd x0 = VectorXd::Zero(stepsX); 
	init_linspace_f(x0,0.,1.);

	for (int i_x = 0; i_x < stepsX; ++i_x)
	{
		val.row(i_x) = logist_it(r,steps,x0[i_x]);
	}
	return val;
}

#endif 











/*#####################OBSOLETE##########################################################
Eigen::VectorXd temp_step_with_Ekin(Eigen::MatrixXd* y,int steps,Eigen::VectorXd E_Kin){
	int n_particles = y[0].col(0).size();
	int dim = int(y[0].row(0).size()/2);
	double kb = 1.;
	double Nf = dim*n_particles-dim;
	std::cout << Nf << std::endl;
	Eigen::VectorXd T(steps);
	for (int i = 0; i < steps; ++i)
	{
		T[i] = (2*E_Kin[i])/(kb*Nf);	
	}
	return T;
}
*/
