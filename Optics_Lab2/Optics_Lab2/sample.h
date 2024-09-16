//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Численные методы
// Шаблон двумерной выборки. 
// sample.h
// 
// Кафедра Прикладной и компьютерной оптики, http://aco.ifmo.ru
// Университет ИТМО
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#if !defined  SAMPLE_H
#define SAMPLE_H

#include <iostream>
#include <vector>
#include <cmath>
#include <exception>
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class PAR>
class Sample
{
protected:
    // число элементов выборки по x и y
    int m_size_x, m_size_y;
    // указатель на динамический массив данных
    std::vector<PAR> m_data;

public:
    // пустой конструктор запрещен
    explicit Sample();
    // полный конструктор
    Sample(int size_x, int size_y);
    // деструктор
    ~Sample();

    // размер выборки по X
    int GetSizeX() const;
    // размер выборки по Y
    int GetSizeY() const;
    // оператор получения значения элемента с номером i,j
    PAR& operator()(int i, int j);
    // оператор получения const значения элемента с номером i,j
    const PAR& operator()(int i, int j) const;

    // получение указателя на всю выборку
    PAR* GetPointer();

    // инициализация всех элементов выборки одинаковым значением
    void SetValue(PAR value);

    // операции с выборкой и числом
    // сложение всех элементов выборки и числа
    Sample<PAR>& operator+=(PAR value);
    // умножение всех элементов выборки на число
    Sample<PAR>& operator*=(PAR value);
    // деление всех элементов выборки на число
    Sample<PAR>& operator/=(PAR value);

    // операции с двумя выборками
    // присваивание
    Sample<PAR>& operator= (const Sample<PAR>& temp);
    // поэлементное сложение двух выборок
    Sample<PAR>& operator+=(const Sample<PAR>& temp);
    // поэлементное перемножение двух выборок
    Sample<PAR>& operator*=(const Sample<PAR>& temp);

    // перегрузка оператора вывода выборки в файл
    template <class PAR>
    friend std::ostream& operator<<(std::ostream& out, Sample<PAR>& other);
    // перегрузка оператора чтения выборки из файла (размерности по X и Y считаем одинаковыми)
    template <class PAR>
    friend std::istream& operator>>(std::istream& out, Sample<PAR>& other);

    // вычисление реальной координаты точки по x
    // i - номер точки
    // dx - шаг по выборке
    double CalcCoordX(int i, double dx = 1.) const;
    // вычисление реальной координаты точки по y
    // i - номер точки
    // dy - шаг по выборке
    double CalcCoordY(int i, double dy = 1.) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// полный конструктор
template <class PAR>
Sample<PAR>::Sample(int size_x, int size_y)
{
    m_size_x = size_x;
    m_size_y = size_y;
    m_data.resize(m_size_x * m_size_y);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// деструктор
template <class PAR>
Sample<PAR>::~Sample()
{
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// размер выборки по X
template <class PAR>
int Sample<PAR>::GetSizeX() const
{
    return(m_size_x);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// размер выборки по Y
template <class PAR>
int Sample<PAR>::GetSizeY() const
{
    return(m_size_y);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// оператор получения значения элемента с номером i,j
template <class PAR>
PAR& Sample<PAR>::operator()(int i, int j)
{
    return (m_data[j * m_size_x + i]);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// оператор получения const значения элемента с номером i,j
template <class PAR>
const PAR& Sample<PAR>::operator()(int i, int j) const
{
    return (m_data[j * m_size_x + i]);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// получение указателя на всю выборку
template <class PAR>
PAR* Sample<PAR>::GetPointer()
{
    return &m_data[0];
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// инициализация всех элементов выборки одинаковым значением
template <class PAR>
void Sample<PAR>::SetValue(PAR value)
{
    for (size_t i = 0; i < m_data.size(); i++)
    {
        m_data[i] = value;
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// сложение всех элементов выборки и числа
template <class PAR>
Sample<PAR>& Sample<PAR>::operator+=(PAR value)
{
    for (size_t i = 0; i < m_data.size(); i++)
    {
        m_data[i] += value;
    }
    return *this;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// умножение всех элементов выборки на число
template <class PAR>
Sample<PAR>& Sample<PAR>::operator*=(PAR value)
{
    for (size_t i = 0; i < m_data.size(); i++)
    {
        m_data[i] *= value;
    }
    return *this;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// деление всех элементов выборки на число
template <class PAR>
Sample<PAR>& Sample<PAR>::operator/=(PAR value)
{
    for (size_t i = 0; i < m_data.size(); i++)
    {
        m_data[i] /= value;
    }
    return *this;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// присваивание
template <class PAR>
Sample<PAR>& Sample<PAR>::operator=(const Sample<PAR>& temp)
{
    m_size_x = temp.m_size_x;
    m_size_y = temp.m_size_y;
    m_data = temp.m_data;

    return *this;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// поэлементное сложение двух выборок
template <class PAR>
Sample<PAR>& Sample<PAR>::operator+=(const Sample<PAR>& temp)
{
    if (m_size_x != temp.m_size_x || m_size_y != temp.m_size_y)
        throw("Sample::operator+= - samples sizes not valid");

    for (size_t i = 0; i < m_data.size(); i++)
    {
        m_data[i] += temp.m_data[i];
    }
    return *this;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// поэлементное перемножение двух выборок
template <class PAR>
Sample<PAR>& Sample<PAR>::operator*=(const Sample<PAR>& temp)
{
    if (m_size_x != temp.m_size_x || m_size_y != temp.m_size_y)
        throw("Sample::operator+= - samples sizes not valid");

    for (size_t i = 0; i < m_data.size(); i++)
    {
        m_data[i] *= temp.m_data[i];
    }
    return *this;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// вычисление реальной координаты точки по X
// i - номер точки
// dx - шаг по выборке
template <class PAR>
double Sample<PAR>::CalcCoordX(int i, double dx) const
{
    return (i - m_size_x / 2.) * dx;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// вычисление реальной координаты точки по Y
// i - номер точки
// dy - шаг по выборке
template <class PAR>
double Sample<PAR>::CalcCoordY(int i, double dy) const
{
    return (i - m_size_y / 2.) * dy;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// перегрузка оператора вывода выборки в файл
template <class PAR>
std::ostream& operator<<(std::ostream& out, Sample<PAR>& other)
{
    for (size_t i = 0; i < other.m_size_x; i++)
    {
        for (size_t j = 0; j < other.m_size_y; j++)
        {
            out << other(i, j) << " ";
        }
        out << std::endl;
    }
    return out;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// перегрузка оператора чтения выборки из файла (размерности по X и Y считаем одинаковыми)
template <class PAR>
std::istream& operator>>(std::istream& out, Sample<PAR>& other)
{
    // чтение временного конейнера из файла
    std::vector<PAR> data;
    PAR value;
    // цикл будет выполняться до тех пор, пока не встретиться конец файла
    while (!out.eof())
    {
        out >> value;           // читаем число из файла
        data.push_back(value);// записываем прочитанное число в конец контейнера
    }

    // проверка размерности прочитанного контейнера
    size_t size = size_t(std::sqrt(double(data.size())));
    if (size * size != data.size())
        throw std::runtime_error("Wrong file size: impossible create square array.");

    // создание выборки вычисленного размера
    other.m_size_x = other.m_size_y = size;
    other.m_data = data;
    return out;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif //defined SAMPLE_H