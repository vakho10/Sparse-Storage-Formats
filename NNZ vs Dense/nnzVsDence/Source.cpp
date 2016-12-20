#include<iostream>
#include<random>
#include<chrono>
#include"string"
#include <filesystem>
#include"Make_data.h"
#include "cg_matrix.h"
#include"cg_jnzSubmatrix.h"
#include <fstream>
#include <sstream>
#include <Windows.h>
//
using namespace std::experimental::filesystem::v1;
using namespace std;

vector<string> split(const string &s, char delim) {
	vector<string> elems;
	stringstream ss(s);
	string item;
	while (getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

int main()
{
	ofstream ofs("results.txt");

	string patternPath = "matricesX3new\\";
	path current_dir(patternPath);

	int count = 0;
	for (auto i = directory_iterator(current_dir); i != directory_iterator(); i++)
	{
		if (!is_directory(i->path())) //we eliminate directories
		{
		
			string fileName = i->path().string();
			cout << fileName << "\t";
			vector<string> x = split(fileName, '\\');

			if (count % 5 == 0)
				ofs << "DIMENSION = " << split(x[1], '_')[0] << endl;

			ofs << split(split(x[1], '_')[1], '.')[0] << "%" << '\t' << '\t';

			long long time = getMinimal1(fileName);

			ofs << "CG_JNZ_SUBMATRIX TIME : " << time << endl;
			cout << "CG_JNZ_SUBMATRIX TIME : " << time << endl;

			count++;
		}
		else continue;
	}
} 
