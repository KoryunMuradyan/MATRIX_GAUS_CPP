#ifndef GAUS_JACOBI_HPP
#define GAUS_JACOBI_HPP

#include <numeric>
#include <map>
#include <iterator>
#include <algorithm>
using namespace math;

template <typename T> bool Not_Zero(T& arg)
{
	return arg != T(NULL);
}

template <typename T> T find_GCD(T arg_num_1, T arg_num_2)
{
	if (arg_num_1 == arg_num_2) {
		return arg_num_1;
	if (arg_num_1 > arg_num_2)
		return find_GCD(arg_num_1-arg_num_2,arg_num_2);
	return find_GCD(arg_num_1, arg_num_2-arg_num_1);
	}
	if (arg_num_1 > arg_num_2) {
		return find_GCD(arg_num_1-arg_num_2,arg_num_2);
	}
	return find_GCD(arg_num_1, arg_num_2-arg_num_1);
}

template <typename T> T find_LCM(const T& arg_num_1, const T& arg_num_2)
{
	T gcd = find_GCD(arg_num_1, arg_num_2);
	return (arg_num_1*arg_num_2)/gcd; 
}

template <typename T> void sort_rows_with_zero(matrix_type<T>& raw_matrix_)
{
	multimap<int, vec_type<T>> tmp_mltmap;
	int row_size = int(raw_matrix_[0].size() - 1);
	for_each(raw_matrix_.begin(), raw_matrix_.end(), 
			[&row_size, &tmp_mltmap](vec_type<T> i) {
			auto front_zero_num = find_if(i.begin(), i.end(), 
					Not_Zero<T>);
			int pos = int(front_zero_num - i.begin());
			if (pos == row_size) {
			pos = -1;
			}
			tmp_mltmap.insert(pair<int, vec_type<T>>(pos, i));
			}
		);
	auto it = raw_matrix_.begin();
	for_each(tmp_mltmap.begin(), tmp_mltmap.end(), 
			[&it](auto i) {
			*it++ = i.second;
			}
		);
}

template <typename T> void define_each_X(std::vector<T>& row, int&& pos, 
		std::map<std::string, T>& _variables)
{
	string str = "x" + std::to_string(pos + 1);
	auto tmp_num = accumulate(row.begin() + pos + 1, row.end() - 1, T(0));
	row[pos] = (row.back()-tmp_num)/(row[pos]);
	_variables.insert(pair<std::string, T>(str, row[pos]));
}

template <typename T>
void put_known_vars(typename matrix_type<T>::iterator& begin, 
		typename matrix_type<T>::iterator& end, int& pos, T& var)
{
	auto a = begin;
	for_each(begin, end, [&](auto& row) {row[pos] *= var;});
}

template <typename T> 
std::map<std::string, T> defineVariables(matrix_type<T>& raw_matrix_)
{
	matrix_type<T> arg_vec(raw_matrix_);
	*(arg_vec[0].end() - 2) = arg_vec[0].back()/(*(arg_vec[0].end() - 2));
	std::map<std::string, T> _variables;
	_variables.insert(pair<std::string, T>("x0", *(arg_vec[0].end() - 2)));
	auto ptr_vars = &_variables;
	auto it_begin = arg_vec.begin() + 1;
	auto it_end = arg_vec.end();
	for_each(arg_vec.begin(), arg_vec.end() - 1, [&](auto& row) {
			auto front_not_zero_num = find_if(row.begin(), 
					row.end(), Not_Zero<T>);
			int pos = int(front_not_zero_num - row.begin());
			put_known_vars(it_begin, it_end, pos, row[pos]);
			define_each_X(*it_begin, std::move(pos - 1), 
					_variables);
			it_begin++;
			}
		);
	return _variables;
}

template <typename T> 
void Transform(std::vector<T>& i, std::vector<T>& tmp_v, 
		int& pos1, std::vector<T> it)
{
	T mult = find_LCM(abs(i[pos1]), abs(tmp_v[pos1]));
	T mult_i = mult / i[pos1];
	T mult_it = mult/tmp_v[pos1];
	transform(i.begin(), i.end(), i.begin(),
			bind1st(std::multiplies<T>(), mult_i));
	transform(it.begin(), it.end(), tmp_v.begin(),
			bind1st(std::multiplies<T>(), mult_it));
	if (mult_i*i[pos1] == -(mult_it*tmp_v[pos1])) {
		transform(i.begin(), i.end(), tmp_v.begin(), i.begin(), 
				plus<T>());
	} else {
		transform(i.begin(), i.end(),
				tmp_v.begin(), i.begin(), minus<T>());
	}

}

template <typename T> void gausHelper(matrix_type<T>& raw_matrix_)
{	
	auto it = raw_matrix_.begin();
	for_each(raw_matrix_.begin(), raw_matrix_.end() - 1, [&](auto& i) {
		it++;
		auto No_zero_num1 = find_if(i.begin(), i.end(), Not_Zero<T>);
		auto No_zero_num2 = find_if(it->begin(), it->end(), Not_Zero<T>);
		int pos1 = int(No_zero_num1 - i.begin());
		int pos2 = int(No_zero_num2 - it->begin());
		vec_type<T> tmp_v = *it;
		if (pos1 != int(i.size())) { 
			if (pos1 == pos2) {
				Transform(i, tmp_v, pos1, *it);
			}
		}
		} );
	if(*((raw_matrix_.begin())->end() - 3) != 0) gausHelper(raw_matrix_);
	else return;
}

template <typename T> 
std::map<std::string, T> gaus_solve(Matrix<T>& arg_matrix)
{	
	matrix_type<T> raw_matrix_ = arg_matrix.get_matrix();
	sort_rows_with_zero(raw_matrix_);
	reverse(raw_matrix_.begin(), raw_matrix_.end());
	gausHelper(raw_matrix_);	
	std::map<std::string, T> _variables = defineVariables(raw_matrix_);
	return _variables;
}

// Gaus helper functions end

#endif
