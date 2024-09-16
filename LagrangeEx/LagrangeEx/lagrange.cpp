#include "lagrange.h"
#include <iostream>

using namespace std;
/**
 * Как полагается, в хедере только объявляем
 * а в .cpp уже расписываем.
 **/

lagrange::lagrange ()
//конструктор
{	
	int n1;
	cout << "\nCreating an object...\n\n";
	/**
	 * Вводим количество точек, которые введём
	 **/
	cout << "Number of dots: ";
	for (;;)
	{
		cin >> n1;
		if (n1 < 1)
		/**
		 * Для одной точки Лагранж с математической т.зр. ещё работает,
		 * а вот для 0 или меньше - уже нет, поэтому следует проверять.
		 * Тем более, что для 0 или меньше ячеек память не выделяеца =)
		 **/
			cout << "N must be more than 1, enter it again: ";
		else break;
	}
	cout << endl;
	n=n1;
	
	x=new double[n];
	y=new double[n];

	for (int i=0; i<n; i++)
	/**
	 * Заполняем сразу парами, каждая пара - и есть координата точки,
	 * например (0; 0), (1; 1), (1.5; 2.25) и т.д.
	 **/
	{
		cout << "X[" << i+1 << "]=";
		cin >> x[i];
		cout << "Y[" << i+1 << "]=";
		cin >> y[i];
	}
}
	
lagrange::~lagrange()
/**
 * Использование в этой программе деконструктора совсем не обязательно
 * поэтому его можно было бы вообще не объявлять и не писать, ни в 
 * lagrange.h, ни здесь.
 **/
{
	delete [] x; delete [] y;
}
	
void lagrange::show()
/**
 * Метод для вывода на экран тех точек, которые мы ввели, заполняя
 * данный объект.
 * **/
{
	cout << "\nThe Table:\n\n";
	for (int i=0; i<n; i++)
	{
		cout << x[i] << "\t";
		cout << y[i] << "\n";
	}
}

double lagrange::polynomial (double argx)
/**
 * Сам полином, реализация формулы, которую можно найти где угодно.
 **/
{
	double c, s=0;
	for (int i=0; i<n; i++)
	{
		c=1;
		for (int j=0; j<n; j++)
			if (i!=j)
				c*=(argx-x[j])/(x[i]-x[j]);
		s+=c*y[i];
	}
	return s;
}
