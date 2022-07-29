#include <iostream>
#include "matrix.hpp"
#include "Read_from_write_into_file.hpp"
#include "Gaus.hpp"

int main(int argc, char** argv)
{
	matrix_type<double> init_vec = MatrixRead<double>(argv[1]);
	math::Matrix<double> obj(init_vec);
	std::map<std::string, double> variables = gaus_solve(obj);
	Generate_Output_File(variables);
	return 0;
}

