#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;
void init(double mass[], double wing_square[], double traction[], double Cd[], double Cl[], double& v, double& h);
void mini(double t[], double &h);

int main()
{
	double mass[3], wing_square[3], traction[3], Cd[3], Cl[3];
	double p, g, v, h;
	cout << fixed << setprecision(2);
	p = 0.17; g = 9.81;
	init(mass, wing_square, traction, Cd, Cl, v, h);
	double L[3];
	for (int i = 0; i < 3; i++) {
		L[i] = 0.5 * p * v * v * wing_square[i] * Cl[i];
	}

	double D[3];
	for (int i = 0; i < 3; i++) {
		D[i] = 0.5 * p * v * v * wing_square[i] * Cd[i];
	}

	double ay[3];
	for (int i = 0; i < 3; i++) {
		ay[i] = (L[i] - mass[i] * g) / mass[i];
	}

	double t[3];
	for (int i = 0; i < 3; i++) {
		t[i] = sqrt(2 * h / ay[i]);
	}
	mini(t, h);

	return 0;
}



void init(double mass[], double wing_square[], double traction[], double Cd[], double Cl[], double& v, double& h)
{

	for (int i = 0; i <= 2; i++)//ìàññà
	{
		cout << "Enter the mass of the aircraft " << i + 1 << " (kg): ";
		cin >> mass[i];
		while (mass[i] <= 0 or cin.fail())
		{
			cin.clear();
			cin.ignore(1000, '\n');
			cout << "Error! Enter a valid positive value: ";
			cin >> mass[i];
		}
		cin.ignore(1000, '\n');
	}

	for (int i = 0; i <= 2; i++)//ïëîùàäü êðûëüåâ
	{
		cout << "Enter the wing square of the aircraft " << i + 1 << " (m^2): ";
		cin >> wing_square[i];
		while (wing_square[i] <= 0 or cin.fail())
		{
			cin.clear();
			cin.ignore(1000, '\n');
			cout << "Error! Enter a valid positive value: ";
			cin >> wing_square[i];
		}
		cin.ignore(1000, '\n');
	}

	for (int i = 0; i <= 2; i++)//òÿãà
	{
		cout << "Enter the traction of the aircraft " << i + 1 << " (N): ";
		cin >> traction[i];
		while (traction[i] <= 0 or cin.fail())
		{
			cin.clear();
			cin.ignore(1000, '\n');
			cout << "Error! Enter a valid positive value: ";
			cin >> traction[i];
		}
		cin.ignore(1000, '\n');
	}

	for (int i = 0; i <= 2; i++)//Cd
	{
		cout << "Enter Cd of the aircraft " << i + 1 << ": ";
		cin >> Cd[i];
		while (Cd[i] <= 0 or cin.fail())
		{
			cin.clear();
			cin.ignore(1000, '\n');
			cout << "Error! Enter a valid positive value: ";
			cin >> Cd[i];
		}
		cin.ignore(1000, '\n');
	}

	for (int i = 0; i <= 2; i++)//Cl
	{
		cout << "Enter Cl of the aircraft " << i + 1 << ": ";
		cin >> Cl[i];
		while (Cl[i] <= 0 or cin.fail())
		{
			cin.clear();
			cin.ignore(1000, '\n');
			cout << "Error! Enter a valid positive value: ";
			cin >> Cl[i];
		}
		cin.ignore(1000, '\n');
	}
	cout << "Enter velocity of the aircrafts (m/s): ";
	cin >> v;
	while (v <= 0 or cin.fail())
	{
		cin.clear();
		cin.ignore(1000, '\n');
		cout << "Error! Enter a valid positive value: ";
		cin >> v;
	}
	cin.ignore(1000, '\n');
	cout << "Enter height (m): ";
	cin >> h;
	while (h <= 0 or cin.fail())
	{
		cin.clear();
		cin.ignore(1000, '\n');
		cout << "Error! Enter a valid positive value: ";
		cin >> h;
	}
	cin.ignore(1000, '\n');
}

void mini(double t[],double &h)
{
	double mik = t[0];
	int n = 0;

		for (int i = 0; i < 3; i++)
		{
			if (mik > t[i])
			{
				mik = t[i];
			}
		}
		int air1 = 0;
		int air2 = 0;
		int air3 = 0;
		for (int i = 0; i < 3; i++)
		{
			if (t[i] == mik)
			{
				if (i == 0)
				{
					air1 = 1;
				}
				if (i == 1) 
				{
					air2 = 2;
				}
				if (i == 2) 
				{
					air3 = 3;
				}
			}
		}
		if (air1 + air2 + air3 == 6)
		{
			cout << "Airplanes gain heighr in the same amount of time"<< endl;
		}
		else if (air1 + air2 + air3 == 3)
		{
			if (air3 != 3)
			{
				cout << "Airplanes 1 and 2 gained height " << h << " meters faster in " << mik << " seconds" << endl;
			}
			if (air3 == 3)
			{
				cout << "Airplan 3 gained height " << h << " meters faster in " << mik << " seconds" << endl;
			}
		}
		else if (air1 + air2 + air3 == 4)
		{
			cout << "Airplanes 1 and 3 gained height " << h << " meters faster in " << mik << " seconds" << endl;
		}
		else if (air1 + air2 + air3 == 5)
		{
			cout << "Airplanes 2 and 3 gained height " << h << " meters faster in " << mik << " seconds" << endl;
		}

		else if (air1 + air2 + air3 == 1)
		{
			cout << "Airplan 1 gained height " << h << " meters faster in " << mik << " seconds" << endl;
		}

		else if (air1 + air2 + air3 == 2)
		{
			cout << "Airplan 2 gained height " << h << " meters faster in " << mik << " seconds" << endl;
		}

		


}
