//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Численные методы
// Результаты дифференцирования
// main.cpp
// 
// Кафедра Прикладной и компьютерной оптики, http://aco.ifmo.ru
// Университет ИТМО
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "diff.h"
#include <fstream>

int main()
{
	std::ofstream fout_dw, fout_eps, fout_step;
	fout_dw.open("D:\\Documents\\VS19projects\\Optics_Lab1\\fout_dw.txt");
	fout_eps.open("D:\\Documents\\VS19projects\\Optics_Lab1\\fout_eps.txt");
	fout_step.open("D:\\Documents\\VS19projects\\Optics_Lab1\\fout_step.txt");

	Function f;
	Function* fm = &f;
	Diff diff(fm);

	double ro, step, dw_dx, dw_dx_1, dw_dx_2, dw_dy, dw_dy_1, dw_dy_2;
	step = 0.05;

	for (ro = -1.0; ro <= 1.0; ro += step)
	{
		dw_dx = f.ProbaDiffRoX(ro, 0.);  // аналитически вычисленная производная 
		dw_dx_1 = diff.DiffOnesideX(ro, 0., step);    // односторонняя разность 
		dw_dx_2 = diff.DiffDoublesideX(ro, 0., step); // двусторонняя разность 

		dw_dy = f.ProbaDiffRoY(0., ro);  // аналитически вычисленная производная 
		dw_dy_1 = diff.DiffOnesideY(0., ro, step);    // односторонняя разность 
		dw_dy_2 = diff.DiffDoublesideY(0., ro, step); // двусторонняя разность 

		// вывод данных для графиков поперечных аберраций 
		fout_dw << ro << " " << dw_dx_1 << " " << dw_dx_2 << " " << dw_dx << " " << dw_dy_1 << " " << dw_dy_2 << " " << dw_dy << std::endl;
		// вывод данных для графиков погрешности от ro 
		fout_eps << ro << " " << dw_dx - dw_dx_1 << " " << dw_dx - dw_dx_2 << " " << dw_dy - dw_dy_1 << " " << dw_dy - dw_dy_2 << std::endl;
	}
	
	ro = 0.9;
	for (step = 0.01; step <= 0.5; step += 0.01)
	{
		dw_dx = f.ProbaDiffRoX(ro, 0.);  // аналитически вычисленная производная
		dw_dx_1 = diff.DiffOnesideX(ro, 0., step);    // односторонняя разность 
		dw_dx_2 = diff.DiffDoublesideX(ro, 0., step); // двусторонняя разность 

		dw_dy = f.ProbaDiffRoY(0., ro);  // аналитически вычисленная производная
		dw_dy_1 = diff.DiffOnesideY(0., ro, step);    // односторонняя разность 
		dw_dy_2 = diff.DiffDoublesideY(0., ro, step); // двусторонняя разность 

		// вывод данных для графиков погрешности от шага 
		fout_step << step << " " << dw_dx - dw_dx_1 << " " << dw_dx - dw_dx_2 << " " << dw_dy - dw_dy_1 << " " << dw_dy - dw_dy_2 << std::endl;
	}
	
	fout_dw.close();
	fout_eps.close();
	fout_step.close();

	return 0;
}