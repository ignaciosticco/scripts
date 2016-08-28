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

/* */


		while(fileIn.good()){
			int cantAtoms;
			int num_iteracion;

			getline(fileIn,line,'P');
			getline(fileIn,line,'I');
			num_iteracion=atoi(line.c_str());

			getline(fileIn,line,'S');
			getline(fileIn,line,'I');
			cantAtoms=atoi(line.c_str());

			int iterador=0;
			float pressure_data;
			float x_data;
			float y_data;
			float speed_data;
			getline(fileIn,line,'I');
			getline(fileIn,line,'d');
			getline(fileIn,line,' ');

			while(iterador<cantAtoms){
				getline(fileIn,line,' ');
				x_data = atof(line.c_str());
				getline(fileIn,line,' ');
				y_data = atof(line.c_str());
				getline(fileIn,line,' ');
				pressure_data = atof(line.c_str());
				getline(fileIn,line,' ');
				speed_data = atof(line.c_str());
				if (num_iteracion>50000){				// Ignotamos los primeros 5s
		    fileOut<< x_data<<" "<< y_data<<" "<< pressure_data<<" "<< speed_data<<endl;
        //cout<< x_data<<" "<< y_data<<" "<< pressure_data<<endl;
				}
        iterador++;

		}


  }
  	fileIn.close();

}
