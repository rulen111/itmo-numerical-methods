#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>

#include <vector>

/**
 * Обращение матрицы с помощью алгоритма Гаусса-Жордана (PARTIAL PIVOT)
 *
 *	Параметры:
 * m -  Матрица для обращения. Д.б. квадратной.
 * singular - Если матрица передана сингулярная, то установите эту
 *            переменную в true, иначе - false.
 *
 * Возвращаемое значение. Если singular - false, то возвращается обратная матрица.
 *						   Если нет, то возвращается матрица со случайными элементами.
 *
 *	Автор метода
 *		Yi Wang, wangy01@mails.tsinghua.edu.cn
 * Взято с сайта
 *		http://www.crystalclearsoftware.com/cgi-bin/boost_wiki/wiki.pl?BOOST_WIKI
 */
template<class T>
boost::numeric::ublas::matrix<T>
gjinverse(const boost::numeric::ublas::matrix<T>& m,
    bool& singular)
{
    using namespace boost::numeric::ublas;

    const int size = m.size1();

    // Cannot invert if non-square matrix or 0x0 matrix.
    // Report it as singular in these cases, and return 
    // a 0x0 matrix.
    if (size != m.size2() || size == 0)
    {
        singular = true;
        matrix<T> A(0, 0);
        return A;
    }

    // Handle 1x1 matrix edge case as general purpose 
    // inverter below requires 2x2 to function properly.
    if (size == 1)
    {
        matrix<T> A(1, 1);
        if (m(0, 0) == 0.0)
        {
            singular = true;
            return A;
        }
        singular = false;
        A(0, 0) = 1 / m(0, 0);
        return A;
    }

    // Create an augmented matrix A to invert. Assign the
    // matrix to be inverted to the left hand side and an
    // identity matrix to the right hand side.
    matrix<T> A(size, 2 * size);
    matrix_range<matrix<T> > Aleft(A,
        range(0, size),
        range(0, size));
    Aleft = m;
    matrix_range<matrix<T> > Aright(A,
        range(0, size),
        range(size, 2 * size));
    Aright = identity_matrix<T>(size);

    // Swap rows to eliminate zero diagonal elements.
    for (int k = 0; k < size; k++)
    {
        if (A(k, k) == 0) // XXX: test for "small" instead
        {
            // Find a row(l) to swap with row(k)
            int l = -1;
            for (int i = k + 1; i < size; i++)
            {
                if (A(i, k) != 0)
                {
                    l = i;
                    break;
                }
            }

            // Swap the rows if found
            if (l < 0)
            {
                singular = true;
                return Aleft;
            }
            else
            {
                matrix_row<matrix<T> > rowk(A, k);
                matrix_row<matrix<T> > rowl(A, l);
                rowk.swap(rowl);
            }
        }
    }

    // Doing partial pivot
    for (int k = 0; k < size; k++)
    {
        // normalize the current row
        for (int j = k + 1; j < 2 * size; j++)
            A(k, j) /= A(k, k);
        A(k, k) = 1;

        // normalize other rows
        for (int i = 0; i < size; i++)
        {
            if (i != k)  // other rows  
            {
                if (A(i, k) != 0)
                {
                    for (int j = k + 1; j < 2 * size; j++)
                        A(i, j) -= A(k, j) * A(i, k);
                    A(i, k) = 0;
                }
            }
        }
    }

    singular = false;
    return Aright;
}