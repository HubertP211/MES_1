#include <iostream>
#include <fstream>
using namespace std;


int liczba_elementow;
int liczba_wezlow;
double global_l;
double global_k;
double global_s;
double global_temp_otoczenia;
double global_alfa;
double global_q;


ifstream plik;


struct wezel
{
	double temp;
	int warunek_brzegowy;
	double x;
};

struct element
{
	int nod1;
	int nod2;
	double s_element;
	double k_element;
	double l_element;
	double h_lokalne[2][2];
	double p_lokalne[2];

};


bool loadData(string nazwaPliku)
{
	
	plik.open(nazwaPliku.c_str());
	if(!plik.good())
	{
		return false;
	}

	while(true)
	{
		plik >> liczba_elementow >>  global_l >> global_k >> global_s >> global_temp_otoczenia >> global_alfa >> global_q;
		break;
	}
	
	liczba_wezlow = liczba_elementow + 1;

}

bool loadData_element(string nazwaPliku, element *tablica_elementow, wezel *tablica_wezlow)
{
	cout<<"//////"<<endl;
	cout<<"laduje elementy"<<endl;
	cout<<"//////"<<endl;
	if(!plik.good())
	{
		return false;
	}

	for(int i = 0; i<liczba_elementow; i++)
	{
		while(true)
		{
			plik >> tablica_elementow[i].nod1;
			plik >> tablica_elementow[i].nod2;
			plik >> tablica_elementow[i].s_element;
			plik >> tablica_elementow[i].k_element;
			tablica_elementow[i].l_element = tablica_wezlow[i+1].x - tablica_wezlow[i].x;
			//tablica_elementow[i].s_element = global_s;
			break;
		}
	}
}

bool loadData_wezel(string nazwaPliku, wezel *tablica_wezlow)
{
	cout<<"//////"<<endl;
	cout<<"laduje wezly"<<endl;
	cout<<"//////"<<endl;
	if(!plik.good())
	{
		return false;
	}

	for(int i = 0; i<liczba_wezlow; i++)
	{
		while(true)
		{
			plik >> tablica_wezlow[i].temp;
			plik >> tablica_wezlow[i].warunek_brzegowy;
			plik >> tablica_wezlow[i].x;
			break;
		}
	}
}


void generateLocalMatrixForElement(element *tablica_elementow, wezel *tablica_wezlow)
{
	for(int i = 0; i<liczba_elementow; i++)
	{
		tablica_elementow[i].h_lokalne[0][0] = (tablica_elementow[i].s_element * tablica_elementow[i].k_element) / tablica_elementow[i].l_element;
		tablica_elementow[i].h_lokalne[0][1] = -((tablica_elementow[i].s_element * tablica_elementow[i].k_element) / tablica_elementow[i].l_element);
		tablica_elementow[i].h_lokalne[1][0] = -((tablica_elementow[i].s_element * tablica_elementow[i].k_element) / tablica_elementow[i].l_element);
		tablica_elementow[i].h_lokalne[1][1] = (tablica_elementow[i].s_element * tablica_elementow[i].k_element) / tablica_elementow[i].l_element;

	}

	for(int i = 0; i<liczba_elementow; i++)
	{
		if(i==0 && tablica_wezlow[tablica_elementow[i].nod1].warunek_brzegowy==1)
		{
			tablica_elementow[i].p_lokalne[0] = global_q * global_s;
			tablica_elementow[i].p_lokalne[1] = 0;
		}
		else if(i==0 && tablica_wezlow[tablica_elementow[i].nod1].warunek_brzegowy==2 )
		{
			tablica_elementow[i].p_lokalne[0] = -(global_alfa * global_temp_otoczenia * global_s);
			tablica_elementow[i].p_lokalne[1] = 0;
		}
		else if(i==liczba_elementow-1 && tablica_wezlow[tablica_elementow[i].nod2].warunek_brzegowy==1 )
		{
			tablica_elementow[i].p_lokalne[0] = 0;
			tablica_elementow[i].p_lokalne[1] = global_q * global_s;
		}
		else if(i==liczba_elementow-1 && tablica_wezlow[tablica_elementow[i].nod2].warunek_brzegowy==2 )
		{
			tablica_elementow[i].p_lokalne[0] = 0;
			tablica_elementow[i].p_lokalne[1] = -(global_alfa * global_temp_otoczenia * global_s);
		}
		else
		{
			tablica_elementow[i].p_lokalne[0] = 0;
			tablica_elementow[i].p_lokalne[1] = 0;
		}
	}
//////////////////////////////////////////////////////////////////////
	cout<<"//////"<<endl;
	cout<<"sprawdzam macierze lokalne elementow"<<endl;
	for(int k = 0; k<liczba_elementow; k++)
	{
	for(int i = 0; i<2; i++)
	{
		for(int j = 0; j<2; j++)
		{
			cout<<tablica_elementow[k].h_lokalne[i][j]<<" ";
		}
		cout<<endl;
	}
	cout<<"///////////////"<<endl;
	}

	cout<<"//////"<<endl;
	cout<<"sprawdzam macierze lokalne warunkow brzegowych elementow"<<endl;
	
	for(int i = 0; i<liczba_elementow; i++)
	{
		for(int j = 0; j<2; j++)
		{
			cout<<tablica_elementow[i].p_lokalne[j]<<endl;
		}
		cout<<endl;
	}
	cout<<"///////////////"<<endl;
	

}


void checkElement(element *tablica_elementow)
{
	cout<<"//////"<<endl;
	cout<<"sprawdzam elementy "<<endl;
	cout<<"//////"<<endl;

	for(int i = 0; i<liczba_elementow; i++)
	{
		cout<<"nod1: "<<tablica_elementow[i].nod1<<endl;
		cout<<"nod2: "<<tablica_elementow[i].nod2<<endl;
		cout<<"s element: "<<tablica_elementow[i].s_element<<endl;
		cout<<"k element: "<<tablica_elementow[i].k_element<<endl;
		cout<<"l element: "<<tablica_elementow[i].l_element<<endl;
		cout<<endl;
	}
}

void checkWezel(wezel *tablica_wezlow)
{
	cout<<"//////"<<endl;
	cout<<"sprawdzam wezly"<<endl;
	cout<<"//////"<<endl;

	for(int i = 0; i<liczba_wezlow; i++)
	{
		cout<<"temp: "<<tablica_wezlow[i].temp<<endl;
		cout<<"warunek brzegowy: "<<tablica_wezlow[i].warunek_brzegowy<<endl;
		cout<<"x: "<<tablica_wezlow[i].x<<endl;
		cout<<endl;

	}
}

void checkData()
{
	cout<<"//////"<<endl;
	cout<<"liczba elementow: "<<liczba_elementow<<endl;
	cout<<"liczba wezlow: "<<liczba_wezlow<<endl;
	cout<<"global_l: "<<global_l<<endl;
	cout<<"global_k: "<<global_k<<endl;
	cout<<"global s: "<<global_s<<endl;
	cout<<"temperatura otoczenia: "<<global_temp_otoczenia<<endl;
	cout<<"global alfa: "<<global_alfa<<endl;
	cout<<"global q: "<<global_q<<endl;
	cout<<"//////"<<endl;
}

void generateFEM_grid(element *tablica_elementow, wezel *tablica_wezlow)
{

	
	if(!loadData_wezel("dane.txt", tablica_wezlow)) cout<<"Blad przy otwieraniu pliku"<<endl;
	checkWezel(tablica_wezlow);

	if(!loadData_element("dane.txt", tablica_elementow, tablica_wezlow)) cout<<"Blad przy otwieraniu pliku"<<endl;
	checkElement(tablica_elementow);

	generateLocalMatrixForElement(tablica_elementow, tablica_wezlow);

}

void calculateMatrixH(double *matrix_h[], element *tablica_elementow, double *warunki)
{
	
	for(int i = 0; i<liczba_elementow+1; i++)
	{
		for(int j = 0; j<liczba_elementow+1; j++)
		{
			matrix_h[i][j] = 0;
		}
	}

	for(int k = 0; k<liczba_elementow; k++)
	{
		for(int i = 0; i<2; i++)
		{
			for(int j = 0; j<2; j++)
			{
				matrix_h[i+k][j+k] += tablica_elementow[k].h_lokalne[i][j];
				cout<<"k: "<<k<<" i: "<<i<<" j: "<<j<<endl;				
			}			
		}
		cout<<k<<" macierz dodana do macierzy globalnej"<<endl;
		
	}
	cout<<"/////////////////"<<endl;
	matrix_h[liczba_elementow][liczba_elementow] += global_alfa * global_s;

	for(int i = 0; i<liczba_elementow+1; i++)
	{
		for(int j = 0; j<liczba_elementow+1; j++)
		{
			cout<<matrix_h[i][j]<<" ";
		}
		cout<<endl;
	}

	for(int i = 0; i<liczba_wezlow; i++)
	{
		warunki[i] = 0;
	}

	warunki[0] = tablica_elementow[0].p_lokalne[0];
	warunki[liczba_wezlow-1] = tablica_elementow[liczba_elementow-1].p_lokalne[1];
	cout<<"////////////////"<<endl;
	for(int i = 0; i<liczba_wezlow; i++)
	{
		cout<<warunki[i]<<endl;
	}
	

}

void calculate(double *matrix_h[], element *tablica_elementow, wezel *tablica_wezlow, double *warunki)
{
	double **gauss = new double *[liczba_elementow+1];
	for(int i = 0; i<liczba_elementow+1; i++)
	{
		gauss[i] = new double[liczba_elementow+2];
	}
	cout<<"utworzono"<<endl;
	////////////////////////////
	cout<<"zaczynam kopiowanie"<<endl;

	for(int i = 0; i<liczba_elementow+1; i++)
	{
		for(int j = 0; j<liczba_elementow+1; j++)
		{
			gauss[i][j] = matrix_h[i][j];
		}
	}
	cout<<"przekopiowano matrix_h"<<endl;
/////////////////////////////////////////////////
	for(int i = 0; i<liczba_elementow+1; i++)
	{
		if(warunki[i]==0) gauss[i][liczba_elementow+1] = warunki[i];
		else gauss[i][liczba_elementow+1] = (-1)*warunki[i];
	}
	cout<<"przekopiowano warunki"<<endl;
	cout<<"sprawdzam poprawnosc ukladu"<<endl;
	for(int i = 0; i<liczba_elementow+1; i++)
	{
		for(int j = 0; j<liczba_elementow+2; j++)
		{
			cout<<gauss[i][j]<<" ";
		}
		cout<<endl;
	}

//////////////////////////////

	cout<<"liczenie"<<endl;

	for(int i = 1; i<liczba_wezlow; i++)
	{
		for(int k=i; k<liczba_wezlow; k++)
		{
			for(int j=i; j<liczba_wezlow+1; j++)
			{
				gauss[k][j] = gauss[k][j] - (gauss[i-1][j]*(gauss[k][i-1]/gauss[i-1][i-1]));
			}
			gauss[k][i-1]=0;
		}
	}
	

	///////////////////

	for(int i = 0; i<liczba_elementow+1; i++)
	{
		for(int j = 0; j<liczba_elementow+2; j++)
		{
			cout<<gauss[i][j]<<" ";
		}
		cout<<endl;
	}


	////////////////
	cout<<"liczenie temperatur"<<endl;

	tablica_wezlow[liczba_wezlow-1].temp = gauss[liczba_wezlow-1][liczba_wezlow] / gauss[liczba_wezlow-1][liczba_wezlow-1];
	
	for(int i = liczba_wezlow-2; i>=0; i--)
	{
		for(int j = 1; j<liczba_wezlow-1; j++)
		{
			if(gauss[i][i] != 0) tablica_wezlow[i].temp += ((-1)*gauss[i][j]*tablica_wezlow[j].temp) / gauss[i][i];
		}
		if(gauss[i][i] != 0) tablica_wezlow[i].temp += gauss[i][liczba_wezlow] / gauss[i][i];


	}

	for(int i = 0; i<liczba_wezlow-1; i++)
	{
		tablica_wezlow[i].temp += tablica_wezlow[liczba_wezlow-1].temp;

	}

	cout<<"wyswietlam temperatury"<<endl;
	for(int i = 0; i<liczba_wezlow; i++)
	{
		cout<<"t"<<i+1<<": "<<tablica_wezlow[i].temp<<" K -> "<<tablica_wezlow[i].temp - 273.15<<" C"<<endl;
	}
}



int main()
{
	
	if(!loadData("dane.txt")) cout<<"Blad przy otwieraniu pliku"<<endl;
	checkData();

	element *tablica_elementow = new element[liczba_elementow];
	wezel *tablica_wezlow = new wezel[liczba_elementow+1];
	double *warunki = new double[liczba_wezlow];

	double **matrix_h = new double * [liczba_elementow+1];
	for(int i = 0; i<liczba_elementow+1; i++)
	{
		matrix_h[i] =new double[liczba_elementow+1];
	}


	generateFEM_grid(tablica_elementow, tablica_wezlow);	
	calculateMatrixH(matrix_h, tablica_elementow, warunki);

	calculate(matrix_h, tablica_elementow, tablica_wezlow, warunki);

	system("pause");


}