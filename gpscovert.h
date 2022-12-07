#pragma once
#ifndef GPSCOVERT_H_
#define GPSCOVERT_H_
#include <iostream>
#include "Eigen/Dense"
using namespace std;

class GpsTo {
public:
	GpsTo();
	virtual ~GpsTo();
	
    double epsilon ;
    double pi ;
    double d2r;
    double r2d ;

    double a ;		//���򳤰���
    double f_inverse ;			//���ʵ���
    double b ;
	//const double b = 6356752.314245;			//����̰���

    double e ;

    double maximumLatitude;

    double lontitude ;
    double latitude  ;

    double mercatorAngleToGeodeticLatitude(double mercatorAngle);
    double geodeticLatitudeToMercatorAngle(double latitude);
	void Blh2Wmc(double& x, double& y, double& z);
	void Wmc2Blh(double& x, double& y, double& z);

    void Gcj02_To_Gps84(double lat, double lon);


};

#endif 
