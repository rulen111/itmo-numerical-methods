/////////////////////////////////////////////////////////////////////////////
// Численные методы
// Класс "Оптическая среда" вычисляет показатель преломления для заданной длины волны 
// по формуле Герцбергера методом наименьших квадратов 
// (см.статью Родионов С.А., Резник В. Г. "Структура и параметры машинного каталога стекла")
// medium.h
// 
// Кафедра Прикладной и компьютерной оптики, http://aco.ifmo.ru
// Университет ИТМО
/////////////////////////////////////////////////////////////////////////////
#if !defined MEDIUM_H
#define MEDIUM_H
/////////////////////////////////////////////////////////////////////////////

#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/vector.hpp"
#include "inverse.h"
#include <cmath>

#include <iostream>

using namespace boost::numeric::ublas;
/////////////////////////////////////////////////////////////////////////////
class Medium
{
private:
    // Матрица L коэффициентов при мю. (Размер 6х6)
    boost::numeric::ublas::matrix<double> m_mL;
    // Матрица весов Q (Размер 6х6)
    boost::numeric::ublas::matrix<double> m_mQ;
    // Вектор-столбец показателей преломления (размер 6)
    boost::numeric::ublas::vector<double> m_vN;
    // вычисленный вектор коэффициентов формулы Герцбергера  (размер 6)
    // (заполняется первый раз при вычислении показателя преломления)
    boost::numeric::ublas::vector<double> m_vMu;

public:

    // Конструктор. Инициализация матриц и векторов.
    // lambda_list - список длин волн, мкм (размер 6)
    // n_list - список известных показателей преломления, соответствующих длинам волн (размер 6)
    Medium(std::vector<double>& lambda_list, std::vector<double>& n_list);
    // Деструктор
    virtual ~Medium();

    // Вычисление показателя преломления оптической среды для произвольной длины волны
    // medium - среда, показатель преломления которой нужно вычислить
    // lambda - длина волны (в мкм)
    // Возвращаемое значение: показатель преломления 
    double GetN(double lambda);
};
/////////////////////////////////////////////////////////////////////////////
Medium::Medium(std::vector<double>& lambda_list, std::vector<double>& n_list)
{
    m_mL.resize(6, 6);
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            switch (j + 1)
            {
            case 1:
                {
                m_mL(i, j) = 1.0;
                break;
                }
            case 2:
                {
                m_mL(i, j) = pow(lambda_list[i], 2);
                break;
                }
            case 3:
                {
                m_mL(i, j) = pow(lambda_list[i], 4);
                break;
                }
            case 4:
                {
                m_mL(i, j) = pow((pow(lambda_list[i], 2) - 0.028), -1);
                break;
                }
            case 5:
                {
                m_mL(i, j) = pow((pow(lambda_list[i], 2) - 0.028), -2);
                break;
                }
            case 6:
                {
                m_mL(i, j) = pow((pow(lambda_list[i], 2) - 0.028), -3);
                break;
                }
            }
        }
    }
    /////////////////////////////////////////////////////////////////////////////
    m_mQ.resize(6, 6);
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            if (lambda_list[i] >= 0.365 && lambda_list[i] <= 0.405)
            {
                if (j == i)
                {
                    m_mQ(i, j) = 1.0;
                }
                else
                {
                    m_mQ(i, j) = 0.0;
                }
            }
            else if (lambda_list[i] >= 0.434 && lambda_list[i] <= 0.657)
            {
                if (j == i)
                {
                    m_mQ(i, j) = 10.0;
                }
                else
                {
                    m_mQ(i, j) = 0.0;
                }
            }
            else if (lambda_list[i] >= 0.7 && lambda_list[i] <= 1.4)
            {
                if (j == i)
                {
                    m_mQ(i, j) = 1.0;
                }
                else
                {
                    m_mQ(i, j) = 0.0;
                }
            }
            else
            {
                if (j == i)
                {
                    m_mQ(i, j) = 0.1;
                }
                else
                {
                    m_mQ(i, j) = 0.0;
                }
            }
        }
    }
    /////////////////////////////////////////////////////////////////////////////
    m_vN.resize(6);
    for (int i = 0; i < 6; i++)
    {
        m_vN[i] = n_list[i];
    }
    /////////////////////////////////////////////////////////////////////////////
    m_vMu.resize(6);
    boost::numeric::ublas::matrix<double> QQ = prod(m_mQ, m_mQ);
    boost::numeric::ublas::matrix<double> LT = trans(m_mL);
    boost::numeric::ublas::matrix<double> LTQQ = prod(LT, QQ);
    boost::numeric::ublas::matrix<double> LTQQL = prod(LTQQ, m_mL);
    boost::numeric::ublas::matrix<double> LQQ(6, 6);
    bool flag = false;
    LQQ = gjinverse<double>(LTQQL, flag);
    boost::numeric::ublas::matrix<double> LQQLT = prod(LQQ, LT);
    boost::numeric::ublas::matrix<double> LQQLTQQ = prod(LQQLT, QQ);
    boost::numeric::ublas::vector<double> N(6);
    for (int i = 0; i < 6; i++)
    {
        N[i] = n_list[i];
    }
    m_vMu = prod(LQQLTQQ, N);
}

Medium::~Medium() {}

double Medium::GetN(double lambda)
{
    double N = 0.0;
    boost::numeric::ublas::vector<double> L(6);
    for (int i = 0; i < 6; i++)
    {
        switch (i + 1)
        {
        case 1:
            {
            L[i] = 1.0;
            break;
            }
        case 2:
            {
            L[i] = pow(lambda, 2);
            break;
            }
        case 3:
            {
            L[i] = pow(lambda, 4);
            break;
            }
        case 4:
            {
            L[i] = pow((pow(lambda, 2) - 0.028), -1);
            break;
            }
        case 5:
            {
            L[i] = pow((pow(lambda, 2) - 0.028), -2);
            break;
            }
        case 6:
            {
            L[i] = pow((pow(lambda, 2) - 0.028), -3);
            break;
            }
        }

        N += m_vMu[i] * L[i];
    }
    return N;
}
/////////////////////////////////////////////////////////////////////////////
#endif  // MEDIUM_H