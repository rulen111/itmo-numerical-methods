//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Численные методы
// Функция: проба функции и проба интеграла
// function.h
// 
// Кафедра Прикладной и компьютерной оптики, http://aco.ifmo.ru
// Университет ИТМО
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#if !defined FUNCTION_H
#define FUNCTION_H

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define PI 3.14159265358979
#define E  2.71828182845905
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <cmath>

class Function
{
public:
    // номер варианта
    int m_variant;
    // конструктор
    Function();
    // деструктор
    ~Function();

    // проба функции f(x)
    double Proba(double x);
    // проба интеграла (вычисленного аналитически) функции f(x) на интервале [a,b] 
    double ProbaIntegral(double a, double b);
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//реализация методов
Function::Function() :m_variant(13) {}

Function::~Function() {}

double Function::Proba(double x)
{
    return (pow(sin(x), 2));
}

double Function::ProbaIntegral(double a, double b)
{
    return ((b / 2 - sin(2 * b) / 4) - (a / 2 - sin(2 * a) / 4));
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif //FUNCTION_H