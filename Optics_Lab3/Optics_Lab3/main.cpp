/////////////////////////////////////////////////////////////////////////////
// Численные методы
// Интерполяция
// main.cpp
// 
// Кафедра Прикладной и компьютерной оптики, http://aco.ifmo.ru
// Университет ИТМО
/////////////////////////////////////////////////////////////////////////////

#include "interpolation.h"

/////////////////////////////////////////////////////////////////////////////
int main()
{
	std::ofstream line("D:\\Documents\\Учёба\\VS19projects\\Optics_Lab3\\ResLineInterp.txt");
	std::ofstream square("D:\\Documents\\Учёба\\VS19projects\\Optics_Lab3\\ResSquareInterp.txt");
	std::ofstream newton("D:\\Documents\\Учёба\\VS19projects\\Optics_Lab3\\ResNewtonInterp.txt");
	std::ofstream lagrange("D:\\Documents\\Учёба\\VS19projects\\Optics_Lab3\\ResLagrangeInterp.txt");
	std::ofstream Delta("D:\\Documents\\Учёба\\VS19projects\\Optics_Lab3\\ResDelta.txt");

	Interpolation inter_line(-1.5 * PI, 10, 0.333 * PI), inter_square(-1.5 * PI, 10, 0.333 * PI), inter_newton(-1.5 * PI, 10, 0.333 * PI), inter_lagrange(-1.5 * PI, 10, 0.333 * PI);

	inter_line.CalcLine(100);
	inter_square.CalcSquare(100);
	inter_newton.CalcNewton(100);
	inter_lagrange.CalcLagrange(100);

	line << inter_line;
	square << inter_square;
	newton << inter_newton;
	lagrange << inter_lagrange;
	Delta << inter_line.CalcSKW() << std::endl << inter_square.CalcSKW() << std::endl << inter_newton.CalcSKW() << std::endl << inter_lagrange.CalcSKW();

	line.close();
	square.close();
	newton.close();
	lagrange.close();
	Delta.close();

	std::ifstream func_in1("D:\\Documents\\Учёба\\VS19projects\\Optics_Lab3\\func_in1.txt");
	std::ifstream func_in2("D:\\Documents\\Учёба\\VS19projects\\Optics_Lab3\\func_in2.txt");
	std::ofstream func_outLI("D:\\Documents\\Учёба\\VS19projects\\Optics_Lab3\\func_outLI.txt");
	std::ofstream func_outSQ("D:\\Documents\\Учёба\\VS19projects\\Optics_Lab3\\func_outSQ.txt");

	Interpolation inter_512line, inter_512square;

	func_in1 >> inter_512line;
	func_in2 >> inter_512square;

	func_in1.close();
	func_in2.close();

	inter_512line.CalcLine(512);
	inter_512square.CalcSquare(512);

	func_outLI << inter_512line;
	func_outSQ << inter_512square;

	func_outLI.close();
	func_outSQ.close();

	return 0;
}