/////////////////////////////////////////////////////////////////////////////
// Численные методы
// Выполнение заданий
// main.cpp
// 
// Кафедра Прикладной и компьютерной оптики, http://aco.ifmo.ru
// Университет ИТМО
/////////////////////////////////////////////////////////////////////////////

#include "fft.h"

/////////////////////////////////////////////////////////////////////////////
int main()
{
	std::ofstream samp1("D:\\Documents\\VS19projects\\Optics_Lab55\\sample1.txt");
	std::ofstream samp2("D:\\Documents\\VS19projects\\Optics_Lab55\\sample2.txt");
	std::ofstream samp2d("D:\\Documents\\VS19projects\\Optics_Lab55\\sample2d.txt");
	std::ofstream samp2d_ref("D:\\Documents\\VS19projects\\Optics_Lab55\\sample2d_ref.txt");
	std::ofstream samp1_im("D:\\Documents\\VS19projects\\Optics_Lab55\\spectr_im1.txt");
	std::ofstream samp2_im("D:\\Documents\\VS19projects\\Optics_Lab55\\spectr_im2.txt");
	std::ofstream samp1_re("D:\\Documents\\VS19projects\\Optics_Lab55\\spectr_re1.txt");
	std::ofstream samp2_re("D:\\Documents\\VS19projects\\Optics_Lab55\\spectr_re2.txt");
	std::ofstream spectr2d("D:\\Documents\\VS19projects\\Optics_Lab55\\spectr2d.txt");
	std::ofstream spectr2d_ref("D:\\Documents\\VS19projects\\Optics_Lab55\\spectr2d_ref.txt");

	SampleComplex sample1(1, 256);
	sample1.Create1DFunc(0.0625);

	SampleComplex sample2(256, 256);
	sample2.CreateCirc(0.0625);

	SampleComplex sample2d(256, 256);
	sample2d.CreatePupil(0.0625);

	SampleComplex sample2d_ref(256, 256);
	sample2d_ref.CreateCirc(0.0625);

	sample1.PrintReal(samp1, sample1);
	sample2.PrintReal(samp2, sample2);
	sample2d.PrintReal(samp2d, sample2d);
	//samp2d << sample2d;
	sample2d_ref.PrintReal(samp2d_ref, sample2d_ref);

	samp1.close();
	samp2.close();
	samp2d.close();
	samp2d_ref.close();

	FFT fft1(256, 1, FFT::EDirection::eDirect, FFT::EDimentions::e1D);
	FFT fft2(256, 256, FFT::EDirection::eDirect, FFT::EDimentions::e2D);
	FFT fft3(256, 256, FFT::EDirection::eReverse, FFT::EDimentions::e2D);
	FFT fft4(256, 256, FFT::EDirection::eReverse, FFT::EDimentions::e2D);

	fft1.Run(sample1);
	fft2.Run(sample2);
	fft3.Run(sample2d);
	fft4.Run(sample2d_ref);

	sample2d.CalcIntens();
	sample2d_ref.CalcIntens();

	sample1.PrintReal(samp1_re, sample1);
	sample1.PrintImag(samp1_im, sample1);
	sample2.PrintReal(samp2_re, sample2);
	sample2.PrintImag(samp2_im, sample2);
	
	sample2d.PrintReal(spectr2d, sample2d);
	sample2d_ref.PrintReal(spectr2d_ref, sample2d_ref);

	samp1_im.close();
	samp2_im.close();
	samp1_re.close();
	samp2_re.close();
	spectr2d.close();
	spectr2d_ref.close();

	return 0;
}