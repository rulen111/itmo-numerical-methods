/////////////////////////////////////////////////////////////////////////////
// Численные методы
// Комплексная двумерная или одномерная выборка - наследник от шаблона Sample<PAR>
// sample_complex.h
// 
// Кафедра Прикладной и компьютерной оптики, http://aco.ifmo.ru
// Университет ИТМО
/////////////////////////////////////////////////////////////////////////////
#if !defined  SAMPLE_COMPLEX_H
#define SAMPLE_COMPLEX_H

#include "sample.h"
#include <complex>

#define PI 3.14159265358979
/////////////////////////////////////////////////////////////////////////////
class SampleComplex : public Sample<std::complex<double> >
{
public:
    // полный конструктор
    SampleComplex(int size_x, int size_y);
    // деструктор
    virtual ~SampleComplex();

    // вычисление интенсивности выборки (квадрат модуля) 
    // записывается в вещественную часть, мнимая часть обнуляется
    void CalcIntens();

    // вывод в файл только вещественной части выборки
    std::ostream& PrintReal(std::ostream& out, SampleComplex& other);
    // вывод в файл только мнимой части выборки
    std::ostream& PrintImag(std::ostream& out, SampleComplex& other);
    // вывод в файл модуля выборки
    std::ostream& PrintModule(std::ostream& out, SampleComplex& other);

    // создание одномерной выборки-функции (dx - шаг по выборке)            x = 1
    void Create1DFunc(double dx);

    // создание выборки-функции circ(ro) 
    void CreateCirc(double dx);
    // создание выборки зрачковой функции 
    void CreatePupil(double dx);

};
/////////////////////////////////////////////////////////////////////////////
SampleComplex::SampleComplex(int size_x, int size_y) : Sample(size_x, size_y) {}

SampleComplex::~SampleComplex() {}

void SampleComplex::CalcIntens()
{
    double z;
    for (int i = 0; i < m_data.size(); i++)
    {
        z = pow((abs(m_data[i])), 2);
        m_data[i].real(z / pow(PI, 2));
        m_data[i].imag(0.0);
    }
}

std::ostream& SampleComplex::PrintReal(std::ostream& out, SampleComplex& other)
{
    for (int i = 0; i < other.m_size_x; i++)
    {
        for (int j = 0; j < other.m_size_y; j++)
        {
            out << std::real(other(i, j)) << " ";
        }
        out << std::endl;
    }
    return out;
}

std::ostream& SampleComplex::PrintImag(std::ostream& out, SampleComplex& other)
{
    for (int i = 0; i < other.m_size_x; i++)
    {
        for (int j = 0; j < other.m_size_y; j++)
        {
            out << std::imag(other(i, j)) << " ";
        }
        out << std::endl;
    }
    return out;
}

std::ostream& SampleComplex::PrintModule(std::ostream& out, SampleComplex& other)
{
    for (int i = 0; i < other.m_size_x; i++)
    {
        for (int j = 0; j < other.m_size_y; j++)
        {
            out << abs(other(i, j)) << " ";
        }
        out << std::endl;
    }
    return out;
}

void SampleComplex::Create1DFunc(double dx)
{
    double n = 1 / pow(dx, 2);
    m_size_x = n;
    m_size_y = 1;
    m_data.resize(m_size_x * m_size_y);
    for (int i = 0; i < n; i++)
    {
        m_data[i] = 1.;
    }
}

void SampleComplex::CreateCirc(double dx)
{
    double n = 1 / pow(dx, 2);
    m_size_x = n;
    m_size_y = n;
    m_data.resize(m_size_x * m_size_y);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if ((pow((i - n / 2) * dx, 2) + pow((j - n / 2) * dx, 2)) <= 1)
            {
                m_data[j * n + i] = 1.;
            }
            else
            {
                m_data[j * n + i] = 0.;
            }
        }
    }
}

void SampleComplex::CreatePupil(double dx)
{
    double n = 1 / pow(dx, 2);
    m_size_x = n;
    m_size_y = n;
    m_data.resize(m_size_x * m_size_y);
    std::complex<double> li;
    li.real(0.0);
    li.imag(1.0);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if ((pow((i - n / 2) * dx, 2) + pow((j - n / 2) * dx, 2)) <= 1)
            {
                m_data[j * n + i] = exp(2 * PI * li * (0.15 * (6 * pow(sqrt(pow(((i - n / 2) * dx), 2) + pow(((j - n / 2) * dx), 2)), 4)
                    - 6 * pow(sqrt(pow(((i - n / 2) * dx), 2) + pow(((j - n / 2) * dx), 2)), 2) + 1)));
            }
            else
            {
                m_data[j * n + i] = 0;
            }
        }
    }
}
/////////////////////////////////////////////////////////////////////////////
#endif //SAMPLE_COMPLEX_H