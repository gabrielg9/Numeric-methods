#include<iostream>
#include<cmath>

using namespace std;

double funkcja(double x, double y)
{
	return y - x * x + 1;
}

void Forward_Euler(double xn, double yn, double dx)
{
	double yn1 = 0;
	for (int i = 0; i < 10; i++)
	{
		yn1 = yn + dx * funkcja(xn, yn);
		xn += dx;
		yn = yn1;
		cout <<"x"<<i+1<<"="<<xn<<"  y"<<i+1<<"="<< yn1 << endl;
	}
}

void Backward_Euler(double xn, double yn, double dx)
{
	double yn1 = 0;
	for (int i = 0; i < 10; i++)
	{
		xn += dx;
		yn = yn1;
		yn1 = yn + dx * funkcja(xn, yn);
		cout << "x" << i + 1 << "=" << xn << "  y" << i + 1 << "=" << yn1 << endl;
	}
}

void Crank_nicolson(double xn, double yn, double dx)
{
	double yn1 = 0;
	for (int i = 0; i < 10; i++)
	{
		double x_next = xn + dx;
		double y_next = yn1;
		yn1 = yn + dx * ((funkcja(xn, yn) + funkcja(x_next, y_next)) / 2);
		cout << "x" << i + 1 << "=" << x_next << "  y" << i + 1 << "=" << yn1 << endl;
		xn += dx;
		yn = yn1;
	}
}

void Extrapolation_Richardson(double xn, double yn, double dx)
{
	double yn1 = 0;
	double x_I = xn;
	double y_I = yn;

	double tab1[10];
	double tab2[20];
	for (int i = 0; i < 10; i++)
	{
		yn1 = y_I + dx * funkcja(x_I, y_I);
		tab1[i] = yn1;
		x_I += dx;
		y_I = yn1;
		//cout << tab1[i] << endl;
	}
	cout << endl;
	dx /= 2.0;
	yn1 = 0;
	x_I = xn;
	y_I = yn;
	for (int i = 0; i < 20; i++)
	{
		yn1 = y_I + dx * funkcja(x_I, y_I);
		tab2[i] = yn1;
		x_I += dx;
		y_I = yn1;
		//cout << tab2[i] << endl;
	}
	dx *= 2.0;
	for (int i = 0; i < 10; i++)
	{
		xn += dx;
		cout << "x" << i + 1 << "=" << xn << "  y" << i + 1 << "=" << ((4 * tab2[2 * i + 1] - tab1[i]) / 3.0) << endl;
	}
}

void Direct_Iteration(double xn, double yn, double dx)
{
	double yn1 = 0;
	double yp = 0;
	double yp1 = 0;
	for (int i = 0; i < 10; i++)
	{
		xn += dx;

		yn1 = yn + dx * funkcja(xn, yn);
		cout << endl << "Predyktor " << yn1 << endl;
		//cout << "x" << i + 1 << "=" << xn << "  y" << i + 1 << "=" << yn1 << endl;
		yp = yn + dx * funkcja(xn, yn1);
		cout << "Korektor " << endl;
		for (int i = 0; ; i++)
		{
			yp1 = yn + dx * funkcja(xn, yp);
			cout << yp << endl;
			if (abs((yp1 - yp) / yn1) < 0.01)
			{
				cout << yp1 << endl;
				break;
			}
			yp = yp1;
		}
		yn = yp1;
	}
}

void Improved_Euler(double xn, double yn, double dx)
{
	double yn1 = 0;
	double yk1 = 0;
	double ysr = 0;
	for (int i = 0; i < 10; i++)
	{
		yn1 = yn + dx * funkcja(xn, yn);
		xn += dx;
		yk1 = yn + dx * funkcja(xn, yn1);
		ysr = (yn1 + yk1) / 2.0;
		cout << "x" << i + 1 << "=" << xn << "  y" << i + 1 << "="<<ysr << endl;
		yn = ysr;
	}
}

void Knutty_4_rzedu(double xn, double yn, double dx)
{
	double yn1 = 0;
	double k1 = 0, k2 = 0, k3 = 0, k4 = 0;
	for (int i = 0; i < 10; i++)
	{
		k1 = dx * funkcja(xn, yn);
		k2 = dx * funkcja(xn + dx * 1.0 / 2.0, yn + k1 * 1.0 / 2.0);
		k3 = dx * funkcja(xn + dx * 1.0 / 2.0, yn + k2 * 1.0 / 2.0);
		k4 = dx * funkcja(xn + dx, yn + k3);
		yn1 = yn + ((k1 + 2 * k2 + 2 * k3 + k4) / 6.0);
		xn += dx;
		cout << "x" << i + 1 << "=" << xn << "  y" << i + 1 << "=" << yn1 << endl;
		yn = yn1;
	}
}

int main()
{
	cout << "Podaj x0, y0, oraz krok dx" << endl;
	double x, y, dx;
	cin >> x >> y >> dx;
	cout << "Forward Euler" << endl;
	Forward_Euler(x, y, dx);
	cout << endl << "Backward Euler" << endl;
	Backward_Euler(x, y, dx);
	cout << endl << "Crank - Nicolson "<< endl;
	Crank_nicolson(x, y, dx);
	cout << endl << "Extralopation Richardson" << endl;
	Extrapolation_Richardson(x, y, dx);
	cout << endl << "Direct Iteration" << endl;
	Direct_Iteration(x, y, dx);
	cout << endl << "Improvd Euler" << endl;
	Improved_Euler(x, y, dx);
	cout << endl << "Knutty 4 rzedu" << endl;
	Knutty_4_rzedu(x, y, dx);

	system("pause");
	return 0;
}