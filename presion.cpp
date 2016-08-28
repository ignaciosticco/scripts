#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <string>
using namespace std;



int main(int argc, char const *argv[])

{
	// El archivo input debe tener presion x y (en columnas)
	string archivoDatos;
	string archivoGuardar;
	cout<<"ingrese el nombre del archivo con los datos: "<<endl;
	cin>>archivoDatos;
	cout<<"ingrese el nombre del archivo donde guardar: "<<endl;
	cin>>archivoGuardar;
	string line;
	ofstream fileOut(archivoGuardar.c_str());
	ifstream fileIn(archivoDatos.c_str());

	int xstep=1;
	int xlong=20; 	// Longitud del recinto
	float ymin= 9.5; // min y-position to take into account
	float ymax= 10.5; // Max y-position to take into account
	vector<int> positionx;
	vector<float> acum_press;
	vector<int> N;
	vector<float> social_pressure;
	int sizepositionx=(xlong/xstep);

	for (int i = 0; i < sizepositionx+1; i++){

		positionx.push_back(xlong-i*xstep);
		acum_press.push_back(0);
		N.push_back(0);
		social_pressure.push_back(0);

	}
/* */


		while(fileIn.good()){
			int cantAtoms;
			int num_iteracion;
			float tiempo;
			getline(fileIn,line,'P');
			getline(fileIn,line,'I');
			num_iteracion=atoi(line.c_str());
			tiempo=num_iteracion*0.0001; 	// Iteraciones* paso temporal
			getline(fileIn,line,'S');
			getline(fileIn,line,'I');
			cantAtoms=atoi(line.c_str());

			int iterador=0;
			float pressure_data;
			float x_data;
			float y_data;
			getline(fileIn,line,'I');
			getline(fileIn,line,'y');
			getline(fileIn,line,' ');

			while(iterador<cantAtoms){
				getline(fileIn,line,' ');
				pressure_data=atof(line.c_str());
				getline(fileIn,line,' ');
				x_data=atof(line.c_str());
				getline(fileIn,line,' ');
				y_data=atof(line.c_str());

				if(0<=x_data<= xlong && ymin<=y_data<ymax){
					/*  POR QUE ESTA MAL ESTO??:
					int i=0;
					while(abs(x_data - positionx[i])> xstep/2){
						i++;
					}

					acum_press[i]=acum_press[i]+pressure_data;
					N[i]=N[i]+1;
					*/
					for (int i=0; i<sizepositionx+1;i++){
						if (abs(x_data-positionx[i])<0.5){
							acum_press[i]=acum_press[i]+pressure_data;
							N[i]=N[i]+1;

						}
					}
				}
			iterador++;
			}

		}

		for (int i=1;  i<sizepositionx+1; i++){
			social_pressure[i]=acum_press[i]/N[i];
		}


	fileIn.close();
	for (int i = 0; i<sizepositionx+1; i++)
	{
		fileOut<< positionx[i]<<" "<< social_pressure[i]<<endl;
	}
	cout<<"End"<<endl;
	return 0;

}
