#include <algorithm>
#include <set>
#include <sstream>
#include <fstream>
#include <iterator>
#include <iostream>
#include <cassert>


std::multiset<std::string> read_from_file(const std::string filename)
{
	std::ifstream file_to_read(filename);
	const std::string str_file((std::istreambuf_iterator<char>(file_to_read)),
					std::istreambuf_iterator<char>());
	file_to_read.close();
	std::istringstream iss(str_file);
	std::multiset<std::string> str_set; 
	while(iss){
		std::string s;
		iss >> s;
		if(std::any_of(s.begin(), s.end(), [](char ch){return std::isdigit(ch);}))
			str_set.insert(s);
	}
	return str_set;
}

void add_to_str(std::string str1, std::string str2, std::string& result_str)
{
	if(str1 == str2){
		str1.resize(15, ' ');
		str2.resize(11, ' ');
		result_str += (str1 + str2 + "   passed\n");
	} else {	
		str1.resize(15, ' ');
		str2.resize(11, ' ');
		result_str += (str1 + str2 + "   failed\n");
	}
}

void generate_result_file(std::multiset<std::string>& set1, std::multiset<std::string>& set2)
{
	assert(set1.size() == set2.size());
	std::string result_str = "expected      output\n\n";
	std::multiset<std::string>::iterator it1 = set1.begin(), it2 = set2.begin();
	while(set1.end() != it1) {
		add_to_str(*it1++, *it2++, result_str);
	}
	std::ofstream result_file("result.txt");
	result_file << result_str;
	result_file.close();
}

void test()
{
	std::multiset<std::string> gold_result, output_result;
	gold_result = read_from_file(std::string("golden.txt"));
	output_result = read_from_file(std::string("output.txt"));
	generate_result_file(gold_result, output_result);

}

int main()
{
	test();
}
