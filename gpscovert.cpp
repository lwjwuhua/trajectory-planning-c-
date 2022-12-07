#include"gpscovert.h"
#include<iostream>
#include <cmath>
using namespace std;

GpsTo::GpsTo() {
	 epsilon = 0.000000000000001;
	 pi = 3.14159265358979323846;
	 d2r = pi / 180;
	 r2d = 180 / pi;

	a = 6378137.0;		//ÍÖÇò³¤°ëÖá
	 f_inverse = 298.257223563;			//±âÂÊµ¹Êý
	 b = a - a / f_inverse;
	//const double b = 6356752.314245;			//ÍÖÇò¶Ì°ëÖá

	e = sqrt(a * a - b * b) / a;
	 maximumLatitude = mercatorAngleToGeodeticLatitude(pi);
	





};
GpsTo::~GpsTo() {};

 double GpsTo:: mercatorAngleToGeodeticLatitude(double mercatorAngle)
{
	return pi / 2.0 - (2.0 * atan(exp(-mercatorAngle)));
	//return 2.0 * atan(exp(mercatorAngle)) - pi / 2.0;
}
 


 double GpsTo::geodeticLatitudeToMercatorAngle(double latitude) {

	 // Clamp the latitude coordinate to the valid Mercator bounds.
	 if (latitude > maximumLatitude)
	 {
		 latitude = maximumLatitude;
	 }
	 else if (latitude < -maximumLatitude)
	 {
		 latitude = -maximumLatitude;
	 }
	 double sinLatitude = sin(latitude);
	 return 0.5 * log((1.0 + sinLatitude) / (1.0 - sinLatitude));

  }

 void GpsTo:: Blh2Wmc(double& x, double& y, double& z)
 {
	 x = x * d2r * a;
	 y = geodeticLatitudeToMercatorAngle(y * d2r) * a;
 }

 void GpsTo:: Wmc2Blh(double& x, double& y, double& z)
 {
	 //var oneOverEarthSemimajorAxis = this._oneOverSemimajorAxis;
	 x = x / a * r2d;
	 y = mercatorAngleToGeodeticLatitude(y / a) * r2d;
 }


 void GpsTo::Gcj02_To_Gps84(double lon, double lat) {

	  double pi = 3.1415926535897932384626;
	  double x_pi = 3.14159265358979324 * 3000.0 / 180.0;
	  double a = 6378245.0;
	//  double a = 6378137.0;
	  double ee = 0.00669342162296594323;

		
		 double ret = -100.0 + 2.0 * (lon - 105.0) + 3.0 * (lat - 35.0) + 0.2 * (lat - 35.0) * (lat - 35.0) + 0.1 * (lon - 105.0) * (lat - 35.0)
			 + 0.2 * sqrt( abs(lon - 105.0));
		 ret += (20.0 * sin(6.0 * (lon - 105.0) * pi) + 20.0 * sin(2.0 * (lon - 105.0) * pi)) * 2.0 / 3.0;
		 ret += (20.0 * sin((lat - 35.0) * pi) + 40.0 * sin((lat - 35.0) / 3.0 * pi)) * 2.0 / 3.0;
		 ret += (160.0 * sin((lat - 35.0) / 12.0 * pi) + 320 * sin((lat - 35.0) * pi / 30.0)) * 2.0 / 3.0;
		 double dLat = lat;

		
		 double ret_ = 300.0 + (lon - 105.0) + 2.0 * (lat - 35.0) + 0.1 * (lon - 105.0) * (lon - 105.0) + 0.1 * (lon - 105.0) * (lat - 35.0) + 0.1
			 * sqrt( abs(lon - 105.0));
		 ret_ += (20.0 * sin(6.0 * (lon - 105.0) * pi) + 20.0 * sin(2.0 * (lon - 105.0) * pi)) * 2.0 / 3.0;
		 ret_ += (20.0 * sin((lon - 105.0) * pi) + 40.0 * sin((lon - 105.0) / 3.0 * pi)) * 2.0 / 3.0;
		 ret_ += (150.0 * sin((lon - 105.0) / 12.0 * pi) + 300.0 * sin((lon - 105.0) / 30.0
			 * pi)) * 2.0 / 3.0;
		 double dLon = ret_;
		


		 double radLat = lat / 180.0 * pi;
		 double magic =  sin(radLat);
		 magic = 1 - ee * magic * magic;
		 double sqrtMagic = sqrt(magic);
		 dLat = (dLat * 180.0) / ((a * (1 - ee)) / (magic * sqrtMagic) * pi);
		 dLon = (dLon * 180.0) / (a / sqrtMagic * cos(radLat) * pi);
		 double mgLat = lat + dLat;
		 double mgLon = lon + dLon;
		



	    lontitude = lon * 2 - mgLon;
	    latitude = lat * 2 - mgLat;
	 

 }