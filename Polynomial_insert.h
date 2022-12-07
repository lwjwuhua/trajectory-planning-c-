#pragma once

#ifndef POLYNOMIAL_INSERT_
#define POLYNOMIAL_INSERT_


#include <cmath>
#include <ctime>
#include <iostream>  
#include<string>
#include <vector>
#include <map>
#include <math.h>
#include "Eigen/Dense"
using Eigen::MatrixXd;
using Eigen::VectorXd;

// #include<graphics.h>
#include<vector>
using namespace std;


class PolyNomial
{
public:
	PolyNomial();
	~PolyNomial();

	vector <VectorXd> pos;
	
	vector<VectorXd> catMullRomSpline(vector<VectorXd> inputPoints);



private:

};






#endif 
