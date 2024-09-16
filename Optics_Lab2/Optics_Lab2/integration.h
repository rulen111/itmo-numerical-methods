//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Численные методы
// Класс численного интегрирования
// integration.h
// 
// Кафедра Прикладной и компьютерной оптики, http://aco.ifmo.ru
// Университет ИТМО
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#if !defined INTEGRATION_H
#define INTEGRATION_H

#include "function.h"
#include "NodeRate.h"
#include <cmath>
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Integration
{
    Function* m_f; // интегрируемая функция 
    NodeRate* m_n; // узлы
    int Cn[6][6];  // вес
    int C[6]{ 1,2,6,8,90,288 };

public:
    // конструктор
    Integration(Function* f, NodeRate* nd);
    // деструктор
    ~Integration();

    // интегрирование методом Ньютона-Котеса на интервале [a,b]
    // count - число подинтервалов
    // order - порядок метода
    double NewtonCotes(double a, double b, unsigned int count, unsigned int order);

    // интегрирование методом Гаусса на интервале [a,b]
    // count - число шагов
    // order - порядок метода
    double Gausse(double a, double b, unsigned int count, unsigned int order);
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Инициализация
Integration::Integration(Function* f, NodeRate* nd) :m_f(f), m_n(nd)
{
    Cn[0][0] = 1; Cn[1][0] = 1; Cn[1][1] = 1; Cn[2][0] = 1; Cn[2][1] = 4; Cn[2][2] = 1; Cn[3][0] = 1; Cn[3][1] = 3; Cn[3][2] = 3; Cn[3][3] = 1;
        Cn[4][0] = 7; Cn[4][1] = 32; Cn[4][2] = 12; Cn[4][3] = 32; Cn[4][4] = 7; Cn[5][0] = 19; Cn[5][1] = 75; Cn[5][2] = 50; Cn[5][3] = 50;
            Cn[5][4] = 75; Cn[5][5] = 19;
}
    // Деструктор
Integration::~Integration() {}
    // Метод Ньютона-Котеса
double Integration::NewtonCotes(double a, double b, unsigned int count, unsigned int order)
{
    double h = ((b - a) / count) / order, sum1 = 0.0, x0 = a, x = 0.0;
    for (unsigned int j = 1; j <= count; j++)
    {
        for (unsigned int i = 0; i <= order; i++)
        {
            x = x0 + i * h;
            sum1 += Cn[order][i] * m_f->Proba(x);
            
        }
        x0 = x;
    }
    return ((order * h / C[order]) * sum1);
}
    // Метод Гаусса
double Integration::Gausse(double a, double b, unsigned int count, unsigned int order)
{
    double sum3 = 0.0;
    for (unsigned int j = 1; j <= count; j++)
    {
        for (unsigned int i = 0; i <= order-1; i++)
        {
            sum3 += m_n->Ci[order - 1][i] * m_f->Proba(m_n->Xi[order - 1][i] + m_n->width * (j - 1.0));
        }
    }
    return (((b - a) / (2.0 * count)) * sum3);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif // INTEGRATION_H