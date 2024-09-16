/////////////////////////////////////////////////////////////////////////////
// Численные методы
// Комплексная двумерная или одномерная выборка - наследник от шаблона Sample<PAR>
// sample_complex.h
// 
// Кафедра Прикладной и компьютерной оптики, http://aco.ifmo.ru
// Университет ИТМО
/////////////////////////////////////////////////////////////////////////////
#if !defined FFT_H
#define FFT_H

#include <fftw3.h>
#include "sample_complex.h"
/////////////////////////////////////////////////////////////////////////////
class FFT
{
public:
    //-----------------------
    // направление преобразования Фурье:
    enum class EDirection
    {
        eDirect,   // прямое 
        eReverse   // обратное 
    };
    // размерность преобразования Фурье:
    enum class EDimentions
    {
        e1D,      // одномерное 
        e2D       // двумерное 
    };
    //-----------------------

private:
    // план вычислений преобразования Фурье
    fftw_plan m_plan = 0;
    // направление (прямое | обратное)
    EDirection m_direction;
    // размерность (одномерное | двумерное)
    EDimentions m_dimension;
    // размер выборки по X, Y
    int m_size_x, m_size_y;


public:
    //конструктор 
    // size_x, size_y - размеры выборки, по умолчанию size_y=1 (одномерная выборка)
    // direction - направление (прямое | обратное), по умолчанию прямое
    // dimension - размерность (одномерное | двумерное), по умолчанию одномерное

    FFT(int size_x, int size_y = 1, EDirection direction = EDirection::eDirect, EDimentions dimension = EDimentions::e1D);
    //деструктор
    ~FFT();

    // функция полного сдвигового преобразования Фурье 
    //		(сдвиг выборки и спектра и учет коэф-та 1/sqrt(n*m) )
    // sample - исходная выборка и результат преобразования 
    void Run(SampleComplex& sample);

private:
    // сдвиг выборки или спектра на size/2 
    // sample - выборка 
    // is_sample: true - sample shift, false - spectre shift
    void ShiftSample(SampleComplex& sample, bool is_sample);
};
/////////////////////////////////////////////////////////////////////////////
FFT::FFT(int size_x, int size_y, EDirection direction, EDimentions dimension)
{
    m_size_x = size_x;
    m_size_y = size_y;
    m_direction = direction;
    m_dimension = dimension;
}

FFT::~FFT() {}

void FFT::Run(SampleComplex& sample)
{
    switch (m_dimension)
    {
    case EDimentions::e1D:
        {
        if (m_direction == EDirection::eDirect)
            {
            m_plan = fftw_plan_dft_1d(m_size_x, (fftw_complex*)sample.GetPointer(), (fftw_complex*)sample.GetPointer(), FFTW_FORWARD, FFTW_ESTIMATE);
            }
        else
            {
            m_plan = fftw_plan_dft_1d(m_size_x, (fftw_complex*)sample.GetPointer(), (fftw_complex*)sample.GetPointer(), FFTW_BACKWARD, FFTW_ESTIMATE);
            }
        break;
        }
    case EDimentions::e2D:
        {
        if (m_direction == EDirection::eDirect)
            {
            m_plan = fftw_plan_dft_2d(m_size_x, m_size_y, (fftw_complex*)sample.GetPointer(), (fftw_complex*)sample.GetPointer(), FFTW_FORWARD, FFTW_ESTIMATE);
            }
        else
            {
            m_plan = fftw_plan_dft_2d(m_size_x, m_size_y, (fftw_complex*)sample.GetPointer(), (fftw_complex*)sample.GetPointer(), FFTW_BACKWARD, FFTW_ESTIMATE);
            }
        break;
        }
    }
    std::complex<double> sz;
    sz.real(1 / sqrt(sample.GetSizeX() * sample.GetSizeY()));
    sz.imag(0.);
    switch (m_direction)
    {
    case EDirection::eDirect:
        {
        ShiftSample(sample, true);
        fftw_execute(m_plan);
        for (int i = 0; i < sample.GetSizeX(); i++)
        {
            for (int j = 0; j < sample.GetSizeY(); j++)
            {
                sample(i, j) *= sz;
            }
        }
        ShiftSample(sample, false);
        break;
        }
    case EDirection::eReverse:
        {
        ShiftSample(sample, false);
        fftw_execute(m_plan);
        for (int i = 0; i < sample.GetSizeX(); i++)
        {
            for (int j = 0; j < sample.GetSizeY(); j++)
            {
                sample(i, j) *= sz;
            }
        }
        ShiftSample(sample, true);
        break;
        }
    }
    fftw_destroy_plan(m_plan);
}

void FFT::ShiftSample(SampleComplex& sample, bool is_sample)
{
    int Nx = sample.GetSizeX(), Ny = sample.GetSizeY();
    double ks = Nx / 2, ms = Nx / 2, ls = Ny / 2, ns = Ny / 2, K, K1, K2;
    std::complex<double> li;
    li.real(0.0);
    li.imag(1.0);
    switch (m_dimension)
    {
    case EDimentions::e1D:
        {
        if (is_sample == true)
        {
            for (int i = 0; i < Nx; i++)
            {
                K = i * ms / Nx;
                sample(i, 0) *= exp(-2 * PI * li * K);
            }
        }
        else
        {
            for (int i = 0; i < Nx; i++)
            {
                K = i * ks / Nx;
                sample(i, 0) *= exp(-2 * PI * li * K);
            }
        }
        break;
        }
    case EDimentions::e2D:
        {
        if (is_sample == true)
        {
            for (int i = 0; i < Nx; i++)
            {
                K1 = i * ms / Nx;
                for (int j = 0; j < Ny; j++)
                {
                    K2 = j * ns / Ny;
                    sample(i, j) *= exp(-2 * PI * li * K1) * exp(-2 * PI * li * K2);
                }
            }
        }
        else
        {
            for (int i = 0; i < Nx; i++)
            {
                K1 = i * ks / Nx;
                for (int j = 0; j < Ny; j++)
                {
                    K2 = j * ls / Ny;
                    sample(i, j) *= exp(-2 * PI * li * K1) * exp(-2 * PI * li * K2);
                }
            }
        }
        break;
        }
    }
}
/////////////////////////////////////////////////////////////////////////////
#endif //defined FFT_H