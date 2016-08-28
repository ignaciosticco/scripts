#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <string>
using namespace std;



int main(int argc, char const *argv[])

{
	// El archivo input debe tener presion x y (en columnas)
	// El output es una tabla de radio vs presion

	string archivoDatos;
	string archivoGuardar;
	cout<<"ingrese el nombre del archivo con los datos: "<<endl;
	cin>>archivoDatos;
	cout<<"ingrese el nombre del archivo donde guardar: "<<endl;
	cin>>archivoGuardar;
	string line;
	ofstream fileOut(archivoGuardar.c_str());
	ifstream fileIn(archivoDatos.c_str());

	float rstep=0.3; // paso del vector R
	float r_delta= 0.1; // radio "diferencial"
	float r_first= 0.0; // primer radio
	vector<float> R;
	vector<double> acum_press;
	vector<double> N;
	vector<double> social_pressure;
	int size_R= 24;
	float door_x= 20.0;	// posicion del centro de la puerta en x
	float door_y= 10.0;	// posicion del centro de la puerta en y

	for (int i = 0; i < size_R+1; i++){

		R.push_back(r_first+rstep*i);
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
			tiempo=num_iteracion*0.0001; 	// Iteraciones*paso temporal
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



					for (int i=0; i<size_R+1;i++){

						float rsq_inf = (R[i]-r_delta)*(R[i]-r_delta);
						float rsq_sup = (R[i]+r_delta)*(R[i]+r_delta);
						float rsq = (x_data-door_x)*(x_data-door_x)+(y_data-door_y)*(y_data-door_y);

						if (rsq_inf < rsq && rsq < rsq_sup && x_data<=door_x){
							acum_press[i]=acum_press[i]+pressure_data;			// presion asociada al iesimo radio
							N[i]=N[i]+1;

						}
					}

			iterador++;
			}

		}

		for (int i=1;  i<size_R+1; i++){

			social_pressure[i]=acum_press[i]/N[i];
		}


	fileIn.close();
	for (int i = 0; i<size_R+1; i++)
	{
		fileOut<< R[i]<<" "<< social_pressure[i]<<endl;
	}
	cout<<"End"<<endl;
	return 0;

}
