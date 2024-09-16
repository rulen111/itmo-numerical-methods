//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Численные методы
// Дифференцирование методами односторонней и двусторонней разности
// diff.h
// 
// Кафедра Прикладной и компьютерной оптики, http://aco.ifmo.ru
// Университет ИТМО
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "function.h"

#if !defined DIFF_H
#define DIFF_H
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Diff
{
public:
	Function* m_f;

	// конструктор
	Diff(Function* f);
	// деструктор
	~Diff();

	// дифференцирование методом односторонней разности по x: f'(x,y)/dx
	double DiffOnesideX(double x, double y, double dx);
	// дифференцирование методом односторонней разности по y: f'(x,y)/dy
	double DiffOnesideY(double x, double y, double dy);
	// дифференцирование методом двухсторонней разности по x: f'(x,y)/dx
	double DiffDoublesideX(double x, double y, double dx);
	// дифференцирование методом двухсторонней разности по y: f'(x,y)/dy
	double DiffDoublesideY(double x, double y, double dy);
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Реализация методов

//Конструктор
Diff::Diff(Function* f) : m_f(f) {}

// деструктор
Diff::~Diff() {}

// дифференцирование методом односторонней разности по x: f'(x,y)/dx
double Diff::DiffOnesideX(double x, double y, double dx)
{
	return ((m_f->ProbaFunction(x + dx, y) - m_f->ProbaFunction(x, y)) / dx);
}

// дифференцирование методом односторонней разности по y: f'(x,y)/dy
double Diff::DiffOnesideY(double x, double y, double dy)
{
	return ((m_f->ProbaFunction(x, y + dy) - m_f->ProbaFunction(x, y)) / dy);
}

// дифференцирование методом двухсторонней разности по x: f'(x,y)/dx
double Diff::DiffDoublesideX(double x, double y, double dx)
{
	return ((m_f->ProbaFunction(x + dx, y) - m_f->ProbaFunction(x - dx, y)) / (2 * dx));
}

// дифференцирование методом двухсторонней разности по y: f'(x,y)/dy
double Diff::DiffDoublesideY(double x, double y, double dy)
{
	return ((m_f->ProbaFunction(x, y + dy) - m_f->ProbaFunction(x, y - dy)) / (2 * dy));
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif //defined DIFF_H