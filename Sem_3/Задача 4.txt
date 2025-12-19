#include<iostream>
#include<iomanip>
#include<cmath>
using namespace std;
void init(double& h, double& ay);

int main()
{
	double h, ay, t;
	cout << fixed << setprecision(2);
	init(h, ay);
	t = sqrt(2 * h / ay);
	cout << "Your parameters: \n" << "Height: " << h << "m\n" << "Horizontal acceleration: " << ay << "m/s^2" << endl;
	cout << "Result: \n" << "Aircraft time to height: " << t << "s" << endl;
	return 0;
}

void init(double& h, double& ay)
{
	cout << "Enter horizontal acceleration (m/s^2): ";
	cin >> ay;
	while (ay <= 0 or cin.fail())
	{
		cin.clear();
		cin.ignore(1000, '\n');
		cout << "Error! Enter a valid positive number: ";
		cin >> ay;
	}
	cin.ignore(1000, '\n');
	cout << "Enter height (m): ";
	cin >> h;
	while (h <= 0 or cin.fail())
	{
		cin.clear();
		cin.ignore(1000, '\n');
		cout << "Error! Enter a valid positive number: ";
		cin >> h;
	}
	cin.ignore(1000, '\n');

}
