/////////////////////////////////////////////////////////////////////////////
// Численные методы
// Класс одномерной выборки
// sample1d.h
// 
// Кафедра Прикладной и компьютерной оптики, http://aco.ifmo.ru
// Университет ИТМО
/////////////////////////////////////////////////////////////////////////////
#if !defined  SAMPLE1D_H
#define SAMPLE1D_H

#include <vector>
#include <fstream>
/////////////////////////////////////////////////////////////////////////////
#define PI 3.14159265358979
#define E  2.71828182845905
/////////////////////////////////////////////////////////////////////////////
class Sample1D
{
public:
    // значения координат x
    std::vector<double> m_x;
    // значения функции 
    std::vector<double> m_y;
    //номер варианта
    int m_variant;              //13       512 elements

    //конструктор
    Sample1D();
    //деструктор
    ~Sample1D();

    // генерация выборки по функции           //f(x)=sin(x)*x         [-1.5PI;1.5PI]  
    // x1 - значение x первой точки
    // n - количество точек
    // dx - шаг по x 
    void SetSample(double x1, int n, double dx);

    void AddX(double x);

    void AddY(double y);

    int GetSize();

    //void ReSize(int n);

    double GetX(int i);

    double GetY(int i);

    //  вычисление значения функции в точке x
    double Y(double x);

    // присваивание выборки
    Sample1D& operator= (const Sample1D& temp);

    // функция ввода и вывода выборки
    friend std::ostream& operator<<(std::ostream& out, Sample1D& one);
    friend std::istream& operator>>(std::istream& out, Sample1D& one);
};
/////////////////////////////////////////////////////////////////////////////
Sample1D::Sample1D() :m_variant(13) {}

Sample1D::~Sample1D() {}

void Sample1D::SetSample(double x1, int n, double dx)
{
    for (int i = 0; i < n; i++)
    {
        m_x.push_back(x1 + i * dx);
        m_y.push_back(sin(x1 + i * dx) * (x1 + i * dx));
    }
}

void Sample1D::AddX(double x)
{
    m_x.push_back(x);
}

void Sample1D::AddY(double y)
{
    m_y.push_back(y);
}

int Sample1D::GetSize()
{
    return m_x.size();
}

/*void Sample1D::ReSize(int n)
{
    m_x.resize(n);
    m_y.resize(n);
}*/

double Sample1D::GetX(int i)
{
    return m_x[i];
}

double Sample1D::GetY(int i)
{
    return m_y[i];
}

double Sample1D::Y(double x)
{
    return (sin(x) * x);
}

Sample1D& Sample1D::operator= (const Sample1D& temp)
{
    m_x = temp.m_x;
    m_y = temp.m_y;
    return *this;
}
/////////////////////////////////////////////////////////////////////////////
// Функция вывода выборки
std::ostream& operator<<(std::ostream& out, Sample1D& one)
{
    for (size_t i = 0; i < one.m_y.size(); i++)
    {
        out << one.m_x[i] << " " << one.m_y[i] << std::endl;
    }
    return out;
}
/////////////////////////////////////////////////////////////////////////////
// Функция чтения выборки
std::istream& operator>>(std::istream& in, Sample1D& one)
{
    double x, y;
    // цикл будет выполняться до тех пор, пока не встретиться конец файла
    while (!in.eof())
    {
        // читаем данные из файла
        in >> x >> y;
        // записываем прочитанное число в конец контейнеров
        one.m_x.push_back(x);
        one.m_y.push_back(y);
    }
    return in;
}
/////////////////////////////////////////////////////////////////////////////
#endif //defined SAMPLE1D_H