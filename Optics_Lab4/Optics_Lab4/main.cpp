/////////////////////////////////////////////////////////////////////////////
// Численные методы
// Вычисление показателя преломления для заданной длины волны 
// по формуле Герцбергера методом наименьших квадратов 
// (см.статью Родионов С.А., Резник В. Г. "Структура и параметры машинного каталога стекла")
// main.cpp
// 
// Кафедра Прикладной и компьютерной оптики, http://aco.ifmo.ru
// Университет ИТМО
/////////////////////////////////////////////////////////////////////////////

#include "medium.h"
#include <fstream>
#include <iostream>

/////////////////////////////////////////////////////////////////////////////
int main()
{
	setlocale(LC_ALL, "Russian");
	std::vector<double> lambda_list(6), n_list(6);
	std::ifstream input_n("D:\\Documents\\VS19projects\\Optics_Lab4\\input_n.txt");
	std::ifstream input_l("D:\\Documents\\VS19projects\\Optics_Lab4\\input_l.txt");

	for (int i = 0; i < 6; i++)
	{
		input_n >> n_list[i];
		input_l >> lambda_list[i];
	}

	input_n.close();
	input_l.close();

	Medium med(lambda_list, n_list);

	std::cout << "ng = " << med.GetN(0.43583) << " nD = " << med.GetN(0.5893) << " n_спец = " << med.GetN(0.8302);

	return 0;
}