/////////////////////////////////////////////////////////////////////////////
// Численные методы
// Класс интерполяции
// interpolation.h
// 
// Кафедра Прикладной и компьютерной оптики, http://aco.ifmo.ru
// Университет ИТМО
/////////////////////////////////////////////////////////////////////////////
#if !defined  INTERPOLATION_H
#define INTERPOLATION_H

#include "boost/numeric/ublas/matrix.hpp"
#include "sample1d.h"
#include <cmath>
/////////////////////////////////////////////////////////////////////////////
class Interpolation
{
public:
    Sample1D m_func; // начальная выборка
    Sample1D m_inter;// полученная выборка

    Interpolation();
    //конструктор
    Interpolation(double x1, int n, double dx);
    //деструктор
    ~Interpolation();

    // Кусочно-линейная интерполяция
    // n - новое количество точек в выборке 
    void CalcLine(int n);
    // Кусочно-квадратичная интерполяция
    // n - новое количество точек в выборке 
    void CalcSquare(int n);

    // интерполяция многочленами Лагранжа
    // y[i]=SUM{y[j]*L[j]}, i=0,n, j=0,k
    // n - новое количество точек в выборке 
    void CalcLagrange(int n);
    // интерполяция многочленами Ньютона
    // y[i]=SUM{c[j]*N[j]}, j=0,n
    // n - новое количество точек в выборке 
    void CalcNewton(int n);

    // вычисление среднеквадратического отклонения восстановленной функции от аналитически заданной 
    double CalcSKW();

    friend std::ostream& operator<<(std::ostream& out, Interpolation& one);

    friend std::istream& operator>>(std::istream& in, Interpolation& one);
};
/////////////////////////////////////////////////////////////////////////////
Interpolation::Interpolation() {}

Interpolation::Interpolation(double x1, int n, double dx)
{
    m_func.SetSample(x1, n, dx);
}

Interpolation::~Interpolation() {}

void Interpolation::CalcLine(int n)
{
    //m_inter.ReSize(n);
    int n0 = m_func.GetSize(), q = 0;
    double n1 = floor(n / (n0 - 1));
    double dx = (m_func.GetX(1) - m_func.GetX(0)) * (n0 * 1.0 - 1) / n, a0, a1, x0, x1;

    for (int i = 1; i < n0; i++)
    {
        x0 = m_func.GetX(i - 1);
        x1 = m_func.GetX(i);
        a1 = (m_func.GetY(i) - m_func.GetY(i - 1)) / (x1 - x0);
        a0 = m_func.GetY(i - 1) - a1 * x0;
        for (int j = 0; j < n1; j++)
        {
            m_inter.AddX(x0);
            m_inter.AddY(a0 + a1 * x0);
            x0 += dx;
            q++;
        }
        if (i == n0 - 1)
        {
            dx = (m_func.GetX(n0 - 1) - x0) / (n * 1.0 - q * 1.0);
            for (int k = 1; k <= n - q; k++)
            {
                m_inter.AddX(x0 + k * dx);
                m_inter.AddY(a0 + a1 * (x0 + k * dx));
            }
        }
    }
}

void Interpolation::CalcSquare(int n)
{
    //m_inter.ReSize(n);
    int n0 = m_func.GetSize(), q = 0;
    double n1 = floor(n / (n0 - 1));
    double dx = (m_func.GetX(1) - m_func.GetX(0)) * (n0 * 1.0 - 1) / n, a0, a1, a2, x0, x1, x2;

    for (int i = 1; i < n0 - 1; i += 2)
    {
        if (i == n0 - 3)
        {
            x0 = m_func.GetX(i - 1);
            x1 = m_func.GetX(i);
            x2 = m_func.GetX(i + 1);
            a2 = ((m_func.GetY(i + 1) - m_func.GetY(i - 1)) / ((x2 - x0) * (x2 - x1))) - ((m_func.GetY(i) - m_func.GetY(i - 1)) / ((x1 - x0) * (x2 - x1)));
            a1 = ((m_func.GetY(i) - m_func.GetY(i - 1)) / (x1 - x0)) - a2 * (x1 + x0);
            a0 = m_func.GetY(i - 1) - a1 * x0 - a2 * pow(x0, 2);
            for (int j = 0; j < n1; j++)
            {
                m_inter.AddX(x0);
                m_inter.AddY(a0 + a1 * x0 + a2 * pow(x0, 2));
                x0 += dx;
                q++;
            }
            i -= 1;
        }
        else
        {
            x0 = m_func.GetX(i - 1);
            x1 = m_func.GetX(i);
            x2 = m_func.GetX(i + 1);
            a2 = ((m_func.GetY(i + 1) - m_func.GetY(i - 1)) / ((x2 - x0) * (x2 - x1))) - ((m_func.GetY(i) - m_func.GetY(i - 1)) / ((x1 - x0) * (x2 - x1)));
            a1 = ((m_func.GetY(i) - m_func.GetY(i - 1)) / (x1 - x0)) - a2 * (x1 + x0);
            a0 = m_func.GetY(i - 1) - a1 * x0 - a2 * pow(x0, 2);
            for (int j = 0; j < 2 * n1; j++)
            {
                m_inter.AddX(x0);
                m_inter.AddY(a0 + a1 * x0 + a2 * pow(x0, 2));
                x0 += dx;
                q++;
            }
            if (i == n0 - 2)
            {
                dx = (m_func.GetX(n0 - 1) - x0) / (n * 1.0 - q * 1.0);
                for (int k = 1; k <= n - q; k++)
                {
                    m_inter.AddX(x0 + k * dx);
                    m_inter.AddY(a0 + a1 * (x0 + k * dx) + a2 * pow((x0 + k * dx), 2));
                }
            }
        }
    }
}

void Interpolation::CalcLagrange(int n)
{
    //m_inter.ReSize(n);
    int n0 = m_func.GetSize();
    double dx = (m_func.GetX(n0 - 1) - m_func.GetX(0)) / (n * 1.0 - 1), x, result, P;
    x = m_func.GetX(0);
    for (int i = 0; i < n; i++)
    {
        result = 0.0;
        for (int j = 0; j < n0; j++)
        {
            P = 1.0;
            for (int k = 0; k < n0; k++)
            {
                if (k != j)
                {
                    P *= (x - m_func.GetX(k)) / (m_func.GetX(j) - m_func.GetX(k));
                }
            }
            result += P * m_func.GetY(j);
        }
        m_inter.AddX(x);
        m_inter.AddY(result);
        x += dx;
    }
}

void Interpolation::CalcNewton(int n)
{
    //m_inter.ReSize(n);
    int n0 = m_func.GetSize();
    double dx = (m_func.GetX(n0 - 1) - m_func.GetX(0)) / (n * 1.0 - 1), P, x = m_func.GetX(0), Xx;
    boost::numeric::ublas::matrix<double> m_B;
    m_B.resize(n0 - 1, n0 - 1);
    for (int i = 0; i < n0 - 1; i++)
    {
        m_B(i, 0) = (m_func.GetY(i + 1) - m_func.GetY(i)) / (m_func.GetX(i + 1) - m_func.GetX(i));
    }
    for (int i = 1; i < n0 - 1; i++)
    {
        for (int j = 0; j < n0 - i - 1; j++)
        {
            m_B(j, i) = (m_B(j + 1, i - 1) - m_B(j, i - 1)) / (m_func.GetX(j + i + 1) - m_func.GetX(j));
        }
    }
    for (int i = 0; i < n; i++)
    {
        Xx = 1.0;
        P = m_func.GetY(0);
        for (int j = 0; j < n0 - 1; j++)
        {
            Xx *= x - m_func.GetX(j);
            P += m_B(0, j) * Xx;
        }
        m_inter.AddX(x);
        m_inter.AddY(P);
        x += dx;
    }
}

double Interpolation::CalcSKW()
{
    int N = m_inter.GetSize();
    double sum1 = 0.0, sum2 = 0.0, xin;
    for (int i = 0; i < N; i++)
    {
        xin = m_inter.GetX(i);
        sum1 += m_inter.GetY(i) - m_func.Y(xin);
        sum2 += pow(m_inter.GetY(i) - m_func.Y(xin), 2);
    }
    sum1 = sum1 / N;
    sum2 = sum2 / N;
    return (sqrt(sum2 - pow(sum1, 2)));
}

std::ostream& operator<<(std::ostream& out, Interpolation& one)
{
    out << one.m_inter;
    return out;
}

std::istream& operator>>(std::istream& in, Interpolation& one)
{
    in >> one.m_func;
    return in;
}
/////////////////////////////////////////////////////////////////////////////
#endif // INTERPOLATION_H