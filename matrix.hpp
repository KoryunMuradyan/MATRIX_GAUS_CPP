#ifndef __MATRIX_HPP__
#define __MATRIX_HPP__

#include <boost/range/algorithm_ext/for_each.hpp>
#include <iterator>
#include <vector>

using namespace std;
template <typename T>
using vec_type = std::vector<T>; 
template <typename T>
using matrix_type = std::vector<std::vector<T>>;

namespace math 
{

template <typename T>
class Matrix
{
	private:
		matrix_type<T> raw_matrix_;
	public:
		Matrix()
		{}
		Matrix(matrix_type<T>&);
		matrix_type<T> get_matrix();
}; // class Matrix

template <typename T>
Matrix<T>::Matrix(matrix_type<T>& arg_vec) 
	: raw_matrix_(arg_vec)
{}

template <typename T>
matrix_type<T> Matrix<T>::get_matrix()
{
	return raw_matrix_;
}

} // namespace math

#endif // __MATRIX_HPP__
