#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

void init(double& m, double& L, double& D, double& T);

int main()
{
	double m, L, D, T, a, ay, g;
	cout << fixed << setprecision(2);
	g = 9.81;
	init(m, L, D, T);
	a = (T - D) / m;
	ay = (L - m * g) / m;
	cout << "\nYour parametrs:" << endl;
	cout << "Mass: " << m << " kg; Lift: " << L << " N; Drag: " << D << " N; Thrust: " << T << " N" << endl;
	cout << "\nResult:" << endl;
	cout << "Horizontal acceleration: " << a << " m/s^2" << endl;
	cout << "Vertical acceleration: " << ay << " m/s^2" << endl;
	return 0;
}


void init(double& m, double& L, double& D, double& T)
{
	cout << "Enter mass (kg): ";
	cin >> m;
	while (m <= 0 or cin.fail())
	{
		cin.clear();
		cin.ignore(1000, '\n');
		cout << "Error! Enter a valid positive number: ";
		cin >> m;
	}

	cout << "Enter lifting force (N): ";
	cin >> L;
	while (L < 0||cin.fail())
	{
		cin.clear();
		cin.ignore(1000, '\n');
		cout << "Error! Enter a valid non-negative number: ";
		cin >> L;
	}

	cout << "Enter drag (N): ";
	cin >> D;
	while (D < 0 || cin.fail())
	{
		cin.clear();
		cin.ignore(1000, '\n');
		cout << "Error! Enter a valid non-negative number: ";
		cin >> D;
	}

	cout << "Enter traction (N): ";
	cin >> T;
	while (T < 0 || cin.fail())
	{
		cin.clear();
		cin.ignore(1000, '\n');
		cout << "Error! Enter a valid non-negative number: ";
		cin >> T;
	}
}
