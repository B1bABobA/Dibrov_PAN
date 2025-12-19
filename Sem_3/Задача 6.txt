#include <iostream>
#include <iomanip>
#include <cmath>
#include <windows.h>
using namespace std;
void init(double* v, double* p, int i, double& Sw, double& Cl);
void lift_force(double* L, double* v, double* p, int i, double& Sw, double& Cl);
void print_results(double* v, double* p, double* L, int i, double Sw, double Cl);

int main() {
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	int i;
	cout << "Введите количество измерений: ";
	cin >> i;
	while (i <= 0 or cin.fail())
	{
		cin.clear();
		cin.ignore(1000, '\n');
		cout << "Ошибка! Введите положительное целое число: ";
		cin >> i;
	}
	cin.ignore(1000, '\n');

	double Sw, Cl;
	double* v = new double[i];
	double* p = new double[i];
	init(v, p, i, Sw, Cl);

	double* L = new double[i];
	lift_force(L, v, p, i, Sw, Cl);
	print_results(v, p, L, i, Sw, Cl);
	delete[]  v;
	delete[] p;
	delete[] L;

	return 0;
}
void init(double* v, double* p, int i, double& Sw, double& Cl)
{
	cout << "Введите площадь крыла (м^2): ";
	cin >> Sw;
	while (Sw <= 0 or cin.fail())
	{
		cin.clear();
		cin.ignore(1000, '\n');
		cout << "Ошибка! Введите положительное число: ";
		cin >> Sw;
	}
	cin.ignore(1000, '\n');

	cout << "Введите коэффициент подъёмной силы: ";
	cin >> Cl;
	while (Cl <= 0 or cin.fail())
	{
		cin.clear();
		cin.ignore(1000, '\n');
		cout << "Ошибка! Введите положительное число: ";
		cin >> Cl;
	}
	cin.ignore(1000, '\n');

	for (int n = 0; n < i; n++)
	{
		cout << "Введите скорость на шаге " <<n+1 << " (м/с^2): ";
		cin >> v[n];
		while (v[n] <= 0 or cin.fail())
		{
			cin.clear();
			cin.ignore(1000, '\n');
			cout << "Ошибка! Введите положительное число: ";
			cin >> v[n];
		}
		cin.ignore(1000, '\n');
	}

	for (int n = 0; n < i; n++)
	{
		cout << "Введите плотность воздуха на шаге " << n + 1 << " (кг/м^3): ";
		cin >> p[n];
		while (p[n] <= 0 or cin.fail())
		{
			cin.clear();
			cin.ignore(1000, '\n');
			cout << "Ошибка! Введите положительное число: ";
			cin >> p[n];
		}
		cin.ignore(1000, '\n');
	}
}

void lift_force(double* L, double* v, double* p, int i, double& Sw, double& Cl)
{
	for (int n = 0; n < i; n++)
	{
		L[n] = 0.5 * p[n] * v[n] * v[n] * Sw * Cl;
	}
}

void print_results(double* v, double* p, double* L, int i, double Sw, double Cl) {
	cout << "Результаты расчёта подъёмной силы" << endl;
	cout << "Введённые параметры:" << endl;
	cout << "Площадь крыла: " << Sw << " м^2" << endl;
	cout << "Коэффициент подъемной силы: " << Cl << endl;
	cout << setw(5) << "Шаг" << setw(20) << "Скорость(м/с)" << setw(22) << "Плотность(кг/м^3)" << setw(9) << "Сила(Н)" << endl;

	for (int n = 0; n < i; n++) {
		cout << setw(5) << n + 1
			<< setw(15) << fixed << setprecision(2) << v[n]
			<< setw(20) << fixed << setprecision(4) << p[n]
			<< setw(15) << fixed << setprecision(2) << L[n] << endl;
	}
}
