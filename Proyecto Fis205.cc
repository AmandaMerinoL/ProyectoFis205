#include "Pythia8/Pythia.h"
#include <iostream>
#include <fstream>
#include <list>
#include <vector>

using namespace Pythia8;


//Funciones que me pueden servir
//Para calcular factoriales
int fact(int num) {
	int factorial=1;
	for (int a=1;a<=num;a++) {
		factorial=factorial*a;
	}
	return factorial;
	}
//Para calcular combinatoria de a sobre 2 (se forman pares de particulas)
int comb(int a) {
	int combi=a*(a-1)/2;
	return combi;
	}


int main() {
	
	//Defino vectores que voy a usar mas adelante. vec donde guardo las particulas. FinJ donde guardo las particulas en el jet final
	//PsJ donde guardo las particulas combinadas en un pseudo jet. rep es un valor para "eliminar" una particula de vec
	std::vector<std::vector<double> > vec;
	std::vector<std::vector<double> > FinJ;
	std::vector<double> aux;
	std::vector<std::vector<double> > PsJ;
	std::vector<std::vector<double> > Part;



  // Numero de eventos
  int nEvent    = 1;
  int nListJets = 5;

  // Generator. LHC process and output selection. Initialization.
  Pythia pythia;
  pythia.readString("Beams:eCM = 14000.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 200.");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.init();

  // Common parameters for the two jet finders.
  double etaMax   = 4.;
  double radius   = 0.7;
  double pTjetMin = 10.;
  // Exclude neutrinos (and other invisible) from study.
  int    nSel     = 2;


  // Set up SlowJet jet finder, with anti-kT clustering
  // and pion mass assumed for non-photons..
  SlowJet slowJet( -1, radius, pTjetMin, etaMax, nSel, 1);

  
  std::list<double>iev;
  //Abro archivo para anotar los datos necesarios
  std::ofstream nev;
  nev.open ("nev.txt");
  nev<<"El numero de eventos nEvents es: "<<nEvent<<".\n"<<"\n";
  
  std::cout<<"El numero de eventos nEvents es: "<<nEvent<<".\n"<<"\n";

 //Se genera el evento, se salta si hay error.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;
    
    
    iev.push_back(iEvent);
    nev <<"Evento numero: "<< iEvent+1<<".\n";
    
    std::vector<std::vector<double> > vec;


    //Se analizan las propiedades del jet.
    slowJet. analyze( pythia.event );
    if (iEvent < nListJets) slowJet.list();
    
    nev <<"size jet: "<< slowJet.sizeJet()<<". Obtendremos "<<comb(slowJet.sizeJet())<<" distancias dij"<<".\n";
    
    std::list<double> dist={};
    //Vector donde se guarda la distancia LdiB con el eje.
	std::vector<double> LdiB;
	
	
	//De las propiedades se extrae pT, y, phi, y se guarda en el vector de particulas Part
	for(int i=0;i<slowJet.sizeJet();i++){
		
		double I=i;
		Part.push_back({I,slowJet.pT(i),slowJet.phi(i),slowJet.y(i)});
	}
	//Imprimo el vector Part en la terminal y en un archivo de texto
	std::cout<<Part.size()<<std::endl;
		std::cout<<"Part"<<std::endl;
	for (int i = 0; i < Part.size(); i++) {
	for (int j = 0; j < Part[i].size(); j++)
		std::cout << Part[i][j] << " ";
		std::cout << std::endl;
	}
	nev<<Part.size()<<std::endl;
		nev<<"Part"<<std::endl;
	for (int i = 0; i < Part.size(); i++) {
	for (int j = 0; j < Part[i].size(); j++)
		nev << Part[i][j] << " ";
		nev << std::endl;}
		
	/*Mi intencion era que los pasos se repitan hasta repasar todas las particulas, pero el while es el que me da ese error*/	
    
	//while(Part.size()>0){	
		
		//Se rellena y se imprime LdiB
    	for (int i = 0; i < slowJet.sizeJet(); ++i){
			double diB=Part[i][1]*Part[i][1];
			LdiB.push_back(diB);
			
		}
	
		std::cout<<"LdiB"<<std::endl;
	for (int i = 0; i < LdiB.size(); i++) {
		std::cout << LdiB[i] << " ";
		std::cout << std::endl;
	}
    
    	//Comienzo a calcular parametros dij y diB.
    	for (int i = 0; i < (Part.size() - 1); ++i)
    	for (int j = i +1; j < Part.size(); ++j) {
    		double dEta = Part[i][3] - Part[j][3];
    		
    		double dPhi = abs( Part[i][2] - Part[j][2] );
    		double dR = sqrt( pow2(dEta) + pow2(dPhi) );
    		double diB=Part[i][1]*Part[i][1];
    	
    		std::vector<double> v1;
    
  	
  	//Se calcula el resto de parametros para guardar un vector de particula ( i , j , dij , diB ). El vector se agrega al vector de vectores vec.
			if(Part[i][1]*Part[i][1] < Part[j][1]*Part[j][1]){
				double dij=Part[i][1]*Part[i][1]*dR*dR/radius;
				dist.push_back(dij);
				nev << "i: "<<i<<" j: "<<j<<" dij: "<< dij <<" diB: "<<diB <<".\n";
				v1.push_back(i);
				v1.push_back(j);
				v1.push_back(dij);
				v1.push_back(diB);
		
				vec.push_back(v1);
			
			} else{
		
		 		double dij=Part[j][1]*Part[j][1]*dR*dR/radius;
		 		dist.push_back(dij);
				nev << "i: "<<i<<" j: "<<j<<" dij: "<< dij <<" diB: "<<diB << ".\n";
				v1.push_back(i);
				v1.push_back(j);
				v1.push_back(dij);
				v1.push_back(diB);
		
		
				vec.push_back(v1);
			}
		
		}
		
		std::vector<double> aux3;	
		std::vector<double> Ldij;
	

		
	//Creo vectores de las distancias dij y diB, para ordenarlas y buscar el minimo valor
		for (int k = 0; k < vec.size(); k++) {
			Ldij.push_back(vec[k][2]);
		}
	
	
  		//Se ordenan los vectores par comparar minimos
		sort(Ldij.begin(),Ldij.end());
		sort(LdiB.begin(),LdiB.end());
	
	
	
	//Segun el menor valor, se mueven las particulas de lista
		if(Ldij[0]<LdiB[0]){
		
			for (int l = 0; l < vec.size(); l++) {
				double Nnum=0;
				double Npt=0;
				double Ny=0;
				double Nphi=0;
				//Se combinan i y j en un pseudo jet
				if(vec[l][2]==Ldij[0]){
					aux=vec[l];
					PsJ.push_back(aux);
					for (int m = 0; m < vec.size(); m++) {
						if(vec[m][0]==aux[0]){
							Nnum=Part[m][0];
							Npt=Npt+Part[m][1];
							Nphi=Nphi+Part[m][2];
							Ny=Ny+Part[m][3];
							Part.erase(vec.begin()+m);
						}
						if(vec[m][1]==aux[1]){
							Nnum=Part[m][0];
							Npt=Npt+Part[m][1];
							Nphi=Nphi+Part[m][2];
							Ny=Ny+Part[m][3];
							Part.erase(vec.begin()+m);
						}
					}
				}
			}
	
		}else{
			//La particula se elimina de Part y se lleva a FinJ (final jet)
			int aa=0;
			for (int i = 0; i < vec.size(); i++) {
				if(vec[i][3]==LdiB[0]){
					aa+=1;
					aux=vec[i];
					
					Part.erase(Part.begin()+i);
				}
			}
			if(aa==0){
				for(int n=0;n<vec.size();n++){
					int len=Part.size()-1.00;
					if(vec[n][1]==len){
						aux=vec[n];
						Part.erase(Part.begin()+n);
						
					}
				}
			
			}
			FinJ.push_back(aux);
			aa=0;
			std::cout<<Part.size()<<std::endl;
		}
		std::cout<<"part size "<<Part.size()<<std::endl;
		//}
		
		


	
		nev<<"Part"<<std::endl;
	for (int i = 0; i < Part.size(); i++) {
	for (int j = 0; j < Part[i].size(); j++)
		nev<< Part[i][j] << " ";
		nev << std::endl;}
		
			std::cout<<"Part"<<std::endl;
	for (int i = 0; i < Part.size(); i++) {
	for (int j = 0; j < Part[i].size(); j++)
		std::cout<< Part[i][j] << " ";
		std::cout << std::endl;}
		
		nev<<"FinJ"<<std::endl;
	for (int i = 0; i < FinJ.size(); i++) {
	for (int j = 0; j < FinJ[i].size(); j++)
		nev<< FinJ[i][j] << " ";
		nev<< std::endl;}
		
			std::cout<<"FinJ"<<std::endl;
	for (int i = 0; i < FinJ.size(); i++) {
	for (int j = 0; j < FinJ[i].size(); j++)
		std::cout<< FinJ[i][j] << " ";
		std::cout<< std::endl;
	}

	nev<<std::endl;

  }
  
  
	
	nev.close();
       
  
  pythia.stat();
  


  // Done.
  return 0;
}
