//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Численные методы
// Пересчет узловых точек и их весовых коэффициентов
// NodeRate.h
// 
// Кафедра Прикладной и компьютерной оптики, http://aco.ifmo.ru
// Университет ИТМО
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#if !defined NODERATE_H
#define NODERATE_H

#include <fstream>
#include <cmath>
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class NodeRate
{
public:
	// Массивы узлов и весовых коэфф
	double Xi[6][6], Ci[6][6], width;
	// Конструктор
	NodeRate();
	// Деструктор
	~NodeRate();

	//пересчет узловых точек
	void NodeCalc(double a, double b, int count);

	//пересчет весовых коэффициентов
	//void RateCalc();

	// Сброс до начальных значений
	void NodeReset();
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Инициализация массивов в конструкторе
NodeRate::NodeRate() :width(0.0), Ci{}, Xi{}
{
	std::ifstream input("D:\\Documents\\Учёба\\VS19projects\\Optics_Lab2\\lab2Xi.txt");
	int n, i;
	for (int j = 0; j < 21; j++)
	{
		input >> n >> i;
		input >> Xi[n][i] >> Ci[n][i];
	}
	input.close();
}
	// Деструктор
NodeRate::~NodeRate() {}
	// Расчет узлов
void NodeRate::NodeCalc(double a, double b, int count)
{
	width = (b - a) / count;
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			Xi[i][j] = a + ((Xi[i][j] + 1) * (width) / 2);
		}
	}
}

/*void NodeRate::RateCalc()
{
	for (int i = 0; i < 6; i++)
	{
		switch (i+1)
		{
		case 1:
			{
			for (int j = 0; j <= i; j++)
				{
				Ci[i][j] = 2 / ((1 - pow(Xi[i][j], 2)) * pow(1, 2));
				}
			break;
			}
		case 2:
			{
			for (int j = 0; j <= i; j++)
			{
				Ci[i][j] = 2 / ((1 - pow(Xi[i][j], 2)) * pow((-(2 * Xi[i][j] - Xi[i][j] * (3 * pow(Xi[i][j], 2) - 1)) / (pow(Xi[i][j], 2) - 1)), 2));
			}
			break;
			}
		case 3:
			{
			for (int j = 0; j <= i; j++)
			{
				Ci[i][j] = 2 / ((1 - pow(Xi[i][j], 2)) * pow(((3 * Xi[i][j] * (pow(Xi[i][j], 3) +
					(3 * Xi[i][j] * (pow(Xi[i][j], 2) - 1)) / 2) - (9 * pow(Xi[i][j], 2) - 3) / 2) / (pow(Xi[i][j], 2) - 1)), 2));
			}
			break;
			}
		case 4:
			{
			for (int j = 0; j <= i; j++)
			{
				Ci[i][j] = 2 / ((1 - pow(Xi[i][j], 2)) * pow(-(4 * pow(Xi[i][j], 3) + 6 * Xi[i][j] * (pow(Xi[i][j], 2) - 1)
					- 4 * Xi[i][j] * (3 * pow(Xi[i][j], 2) * (pow(Xi[i][j], 2) - 1) + (3 * pow(pow(Xi[i][j], 2) - 1, 2) / 8)
						+ pow(Xi[i][j], 4))) / (pow(Xi[i][j], 2) - 1), 2));
			}
			break;
			}
		case 5:
			{
			for (int j = 0; j <= i; j++)
			{
				Ci[i][j] = 2 / ((1 - pow(Xi[i][j], 2)) * pow(-(15 * pow(Xi[i][j], 2) * (pow(Xi[i][j], 2) - 1)
					- 5 * Xi[i][j] * ((15 * Xi[i][j] * pow(pow(Xi[i][j], 2) - 1, 2)) / 8 + 5 * pow(Xi[i][j], 3) * (pow(Xi[i][j], 2) - 1)
						+ pow(Xi[i][j], 5)) + (15 * pow(pow(Xi[i][j], 2) - 1, 2)) / 8 + 5 * pow(Xi[i][j], 4)) / (pow(Xi[i][j], 2) - 1), 2));
			}
			break;
			}
		case 6:
			{
			for (int j = 0; j <= i; j++)
			{
				Ci[i][j] = 2 / ((1 - pow(Xi[i][j], 2)) * pow(-((45 * Xi[i][j] * pow(pow(Xi[i][j], 2) - 1, 2)) / 4
					- 6 * Xi[i][j] * ((15 * pow(Xi[i][j], 4) * (pow(Xi[i][j], 2) - 1)) / 2 + (5 * pow(pow(Xi[i][j], 2) - 1, 3)) / 16 + pow(Xi[i][j], 6)
						+ (45 * pow(Xi[i][j], 2) * pow(pow(Xi[i][j], 2) - 1, 2)) / 8) + 30 * pow(Xi[i][j], 3) * (pow(Xi[i][j], 2) - 1)
							+ 6 * pow(Xi[i][j], 5)) / (pow(Xi[i][j], 2) - 1), 2));
			}
			break;
			}
		}
		
	}
}*/
	// Сброс
void NodeRate::NodeReset()
{
	std::ifstream input("D:\\Documents\\Учёба\\VS19projects\\Optics_Lab2\\lab2Xi.txt");
	int n, i;
	for (int j = 0; j < 21; j++)
	{
		input >> n >> i;
		input >> Xi[n][i] >> Ci[n][i];
	}
	input.close();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif //NODERATE_H