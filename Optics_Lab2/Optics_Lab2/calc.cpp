//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Численные методы
// Главный расчет
// calc.cpp
// 
// Кафедра Прикладной и компьютерной оптики, http://aco.ifmo.ru
// Университет ИТМО
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "integration.h"
#include "sample.h"
#include <iostream>
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main()
{
	// Первое задание
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Файлы вывода
	std::ofstream NewtonCotes, Gausse;
	NewtonCotes.open("D:\\Documents\\Учёба\\VS19projects\\Optics_Lab2\\NewtonCotes.txt");
	Gausse.open("D:\\Documents\\Учёба\\VS19projects\\Optics_Lab2\\Gausse.txt");
	// Конструкторы
	NodeRate nodes;
	NodeRate* nd = &nodes;
	Function func;
	Function* f = &func;
	Integration integ(f, nd);
	// Границы интегр и другие переменные
	double a = 0.0, b = 3 * PI / 4, I_an, I_nc1, I_nc2, I_nc3, I_nc4, I_nc5, I_g1, I_g2, I_g3, I_g4, I_g5, I_g6;
	unsigned int count;
	// Аналитический интеграл
	I_an = func.ProbaIntegral(a, b);
	// Цикл расчета и вывода результатов
	for (count = 1; count <= 100; count++)
	{
		I_nc1 = integ.NewtonCotes(a, b, count, 1);
		I_nc2 = integ.NewtonCotes(a, b, count, 2);
		I_nc3 = integ.NewtonCotes(a, b, count, 3);
		I_nc4 = integ.NewtonCotes(a, b, count, 4);
		I_nc5 = integ.NewtonCotes(a, b, count, 5);

		NewtonCotes << count << " " << I_an - I_nc1 << " " << I_an - I_nc2 << " " << I_an - I_nc3 << " " << I_an - I_nc4 << " " << I_an - I_nc5 << std::endl;

		nodes.NodeCalc(a, b, count);
	
		I_g1 = integ.Gausse(a, b, count, 1);
		I_g2 = integ.Gausse(a, b, count, 2);
		I_g3 = integ.Gausse(a, b, count, 3);
		I_g4 = integ.Gausse(a, b, count, 4);
		I_g5 = integ.Gausse(a, b, count, 5);
		I_g6 = integ.Gausse(a, b, count, 6);

		Gausse << count << " " << I_an - I_g1 << " " << I_an - I_g2 << " " << I_an - I_g3 << " " << I_an - I_g4 << " " << I_an - I_g5 << " " << I_an - I_g6 << std::endl;

		nodes.NodeReset();
	}

	NewtonCotes.close();
	Gausse.close();
	// Второе задание
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Конструктор
	Sample<double> data(256, 256);
	// Ввод ФРТ из файла
	std::ifstream frtin("D:\\Documents\\Учёба\\VS19projects\\Optics_Lab2\\lab2frt.txt");
	frtin >> data;

	std::ofstream fkeout("D:\\Documents\\Учёба\\VS19projects\\Optics_Lab2\\lab2frtout.txt");
	
	double suminf = 0.0, sum;
	int size = data.GetSizeX();

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			suminf += data(i, j);
		}
	}

	for (int i = 0; i < size/2; i ++)
	{
		sum = 0.0;
		for (int j = size / 2 - i; j <= size / 2 + i; j++)
		{
			for (int k = size / 2 - i; k <= size / 2 + i; k++)
			{
				sum += data(j, k);
			}
		}
		fkeout << sum / suminf << std::endl;
	}

	frtin.close();
	fkeout.close();

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	return 0;
}