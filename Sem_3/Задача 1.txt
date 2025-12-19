#include <iostream>
#include <iomanip>
using std::cin;
using std::cout;
using std::endl;
using std::setw;

int main()
{
	double  s, v, p, Cl;
	cout << "Enter squier (m): ";
	cin >> s;
	while (s < 0 or cin.fail())
	{
		cin.clear();
		cin.ignore(1000, '\n');
		cout << "Error! Enter a valid non-negative number: ";
		cin >> s;
	}
	cout << "Enter velocity (m/s): ";
	cin >> v;
	while (v < 0 or cin.fail())
	{
		cin.clear();
		cin.ignore(1000, '\n');
		cout << "Error! Enter a valid non-negative number: ";
		cin >> v;
	}

	cout << "Enter density (kg/m^3):";
	cin >> p;
	while (p < 0 or cin.fail())
	{
		cin.clear();
		cin.ignore(1000, '\n');
		cout << "Error! Enter a valid non-negative number: ";
		cin >> p;
	}

	cout << "Enter lifting force coefficient: ";
	cin >> Cl;

	while (Cl < 0 or cin.fail())
	{
		cin.clear();
		cin.ignore(1000, '\n');
		cout << "Error! Enter a valid non-negative number: ";
		cin >> Cl;
	}
	double L = 0.5 * v * v * p * s * Cl;
	cout << "Result L = " << L << " N" << endl;
	return 0;
}

