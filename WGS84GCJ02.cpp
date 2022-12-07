#define _WGSTOGCJ_API
#include "WGS84GCJ02.h"
#include <iostream>
#include <cmath>
#include "WGS84GCJ02.h"
#include <ceres/ceres.h>

constexpr double kKRASOVSKY_A = 6378245.0;				// equatorial radius [unit: meter]
constexpr double kKRASOVSKY_B = 6356863.0187730473;	// polar radius
constexpr double kKRASOVSKY_ECCSQ = 6.6934216229659332e-3; // first eccentricity squared
constexpr double kKRASOVSKY_ECC2SQ = 6.7385254146834811e-3; // second eccentricity squared
constexpr double PI = 3.14159265358979323846;   //¦Ð

constexpr double kDEG2RAD = PI / 180.0;
constexpr double kRAD2DEG = 180.0 / PI;

constexpr inline double Deg2Rad(const double deg) {
	return deg * kDEG2RAD;
}

constexpr inline double Rad2Deg(const double rad) {
	return rad * kRAD2DEG;

}

std::pair<double, double> GetGeodeticOffset(const double& wgs84lon, const double& wgs84lat)
{
	//get geodetic offset relative to 'center china'
	double lon0 = wgs84lon - 105.0;
	double lat0 = wgs84lat - 35.0;

	//generate an pair offset roughly in meters
	double lon1 = 300.0 + lon0 + 2.0 * lat0 + 0.1 * lon0 * lon0 + 0.1 * lon0 * lat0 + 0.1 * sqrt(fabs(lon0));
	lon1 = lon1 + (20.0 * sin(6.0 * lon0 * PI) + 20.0 * sin(2.0 * lon0 * PI)) * 2.0 / 3.0;
	lon1 = lon1 + (20.0 * sin(lon0 * PI) + 40.0 * sin(lon0 / 3.0 * PI)) * 2.0 / 3.0;
	lon1 = lon1 + (150.0 * sin(lon0 / 12.0 * PI) + 300.0 * sin(lon0 * PI / 30.0)) * 2.0 / 3.0;
	double lat1 = -100.0 + 2.0 * lon0 + 3.0 * lat0 + 0.2 * lat0 * lat0 + 0.1 * lon0 * lat0 + 0.2 * sqrt(fabs(lon0));
	lat1 = lat1 + (20.0 * sin(6.0 * lon0 * PI) + 20.0 * sin(2.0 * lon0 * PI)) * 2.0 / 3.0;
	lat1 = lat1 + (20.0 * sin(lat0 * PI) + 40.0 * sin(lat0 / 3.0 * PI)) * 2.0 / 3.0;
	lat1 = lat1 + (160.0 * sin(lat0 / 12.0 * PI) + 320.0 * sin(lat0 * PI / 30.0)) * 2.0 / 3.0;

	//latitude in radian
	double B = Deg2Rad(wgs84lat);
	double sinB = sin(B), cosB = cos(B);
	double W = sqrt(1 - kKRASOVSKY_ECCSQ * sinB * sinB);
	double N = kKRASOVSKY_A / W;

	//geodetic offset used by GCJ-02
	double lon2 = Rad2Deg(lon1 / (N * cosB));
	double lat2 = Rad2Deg(lat1 * W * W / (N * (1 - kKRASOVSKY_ECCSQ)));
	return { lon2, lat2 };


}


WGSTOGCJ_API_EXPORTS std::pair<double, double> Wgs2Gcj(const double& wgs84lon, const double& wgs84lat)
{
	std::pair<double, double> dlonlat = GetGeodeticOffset(wgs84lon, wgs84lat);
	double gcj02lon = wgs84lon + dlonlat.first;
	double gcj02lat = wgs84lat + dlonlat.second;
	return { gcj02lon, gcj02lat };
		

}

std::pair<double, double> Gcj2Wgs_SimpleIteration(const double& gcj02lon, const double& gcj02lat)
{
	std::pair<double, double> lonlat = Wgs2Gcj(gcj02lon, gcj02lat);
	double lon0 = lonlat.first;
	double lat0 = lonlat.second;
	int iterCounts = 0;
	while (++iterCounts < 1000)
	{
		std::pair<double, double> lonlat1 = Wgs2Gcj(lon0, lat0);
		double lon1 = lonlat1.first;
		double lat1 = lonlat1.second;
		double dlon = lon1 - gcj02lon;
		double dlat = lat1 - gcj02lat;
		lon1 = lon0 - dlon;
		lat1 = lat0 - dlat;
		//1.0e-9 degree corresponding to 0.1mm
		if (fabs(dlon) < 1.0e-9 && fabs(dlat) < 1.0e-9)
			break;
		lon0 = lon1;
		lat0 = lat1;
	}
	return { lon0 , lat0 };
	

}

class AutoDiffCostFunc
{
public:
	AutoDiffCostFunc(const double lon, const double lat) :
		mlonGcj(lon), mlatGcj(lat) {}
	template <typename T>
	bool operator() (const T* const lonWgs, const T* const latWgs, T* residuals) const
	{
		//get geodetic offset relative to 'center china'
		T lon0 = lonWgs[0] - T(105.0);
		T lat0 = latWgs[0] - T(35.0);

		//generate an pair offset roughly in meters
		T lon1 = T(300.0) + lon0 + T(2.0) * lat0 + T(0.1) * lon0 * lon0 + T(0.1) * lon0 * lat0
			+ T(0.1) * ceres::sqrt(ceres::abs(lon0));
		lon1 = lon1 + (T(20.0) * ceres::sin(T(6.0) * lon0 * T(PI)) + T(20.0) * ceres::sin(T(2.0) * lon0 * T(PI))) * T(2.0) / T(3.0);
		lon1 = lon1 + (T(20.0) * ceres::sin(lon0 * T(PI)) + T(40.0) * ceres::sin(lon0 / T(3.0) * T(PI))) * T(2.0) / T(3.0);
		lon1 = lon1 + (T(150.0) * ceres::sin(lon0 / T(12.0) * T(PI)) + T(300.0) * ceres::sin(lon0 * T(PI) / T(30.0))) * T(2.0) / T(3.0);
		T lat1 = T(-100.0) + T(2.0) * lon0 + T(3.0) * lat0 + T(0.2) * lat0 * lat0 + T(0.1) * lon0 * lat0
			+ T(0.2) * ceres::sqrt(ceres::abs(lon0));
		lat1 = lat1 + (T(20.0) * ceres::sin(T(6.0) * lon0 * T(PI)) + T(20.0) * ceres::sin(T(2.0) * lon0 * T(PI))) * T(2.0) / T(3.0);
		lat1 = lat1 + (T(20.0) * ceres::sin(lat0 * T(PI)) + T(40.0) * ceres::sin(lat0 / T(3.0) * T(PI))) * T(2.0) / T(3.0);
		lat1 = lat1 + (T(160.0) * ceres::sin(lat0 / T(12.0) * T(PI)) + T(320.0) * ceres::sin(lat0 * T(PI) / T(30.0))) * T(2.0) / T(3.0);

		//latitude in radian
		T B = latWgs[0] * T(kDEG2RAD);
		T sinB = ceres::sin(B), cosB = ceres::cos(B);
		T W = ceres::sqrt(T(1) - T(kKRASOVSKY_ECCSQ) * sinB * sinB);
		T N = T(kKRASOVSKY_A) / W;

		//geodetic offset used by GCJ-02
		T lon2 = T(kRAD2DEG) * lon1 / (N * cosB);
		T lat2 = T(kRAD2DEG) * (lat1 * W * W / (N * (1 - kKRASOVSKY_ECCSQ)));

		//residuals
		residuals[0] = lonWgs[0] + lon2 - mlonGcj;
		residuals[1] = latWgs[0] + lat2 - mlatGcj;
		return true;
	}

private:
	double mlonGcj;
	double mlatGcj;
};



WGSTOGCJ_API_EXPORTS std::pair<double, double> Gcj2Wgs_AutoDiff(const double& gcj02lon, const double& gcj02lat)
{
	ceres::Problem* poProblem = new ceres::Problem;
	AutoDiffCostFunc* pCostFunc = new AutoDiffCostFunc(gcj02lon, gcj02lat);

	double wgslon = gcj02lon, wgslat = gcj02lat;
	poProblem->AddResidualBlock(new ceres::AutoDiffCostFunction<AutoDiffCostFunc, 2, 1, 1>(pCostFunc),
		nullptr,
		&wgslon,
		&wgslat);

	ceres::Solver::Options options;
	options.max_num_iterations = 30;
	options.linear_solver_type = ceres::DENSE_QR;
	options.minimizer_progress_to_stdout = false;
	options.gradient_tolerance = 1e-16;
	options.function_tolerance = 1e-12;
	options.parameter_tolerance = 1e-14;
	ceres::Solver::Summary summary;
	ceres::Solve(options, poProblem, &summary);
	delete poProblem;		//auto free memory of cost function "pCostFunc"
	return { wgslon, wgslat };
	


}

std::pair<double, double> GetPartialDerivative_Lon(const double& wgs84lon, const double& wgs84lat, const double& dlon)
{
	double lonBk = wgs84lon + dlon;
	double lonFw = wgs84lon - dlon;
	std::pair<double, double > gcjlatlonBk = Wgs2Gcj(lonBk, wgs84lat);
	double gcjlonBk = gcjlatlonBk.first;
	double gcjlatBk = gcjlatlonBk.second;

	std::pair<double, double > gcjlatlonFw = Wgs2Gcj(lonFw, wgs84lat);
	double gcjlonFw = gcjlatlonFw.first;
	double gcjlatFw = gcjlatlonFw.second;
	double dlongcj_dlonwgs = (gcjlonBk - gcjlonFw) / (dlon * 2.0);
	double dlatgcj_dlonwgs = (gcjlatBk - gcjlatFw) / (dlon * 2.0);
	return { dlongcj_dlonwgs , dlatgcj_dlonwgs };
}

std::pair<double, double> GetPartialDerivative_Lat(const double& wgs84lon, const double& wgs84lat, const double& dlat)
{
	double latBk = wgs84lat + dlat;
	double latFw = wgs84lat - dlat;
	std::pair<double, double> gcjlonlatBk = Wgs2Gcj(wgs84lon, latBk);
	double gcjlonBk = gcjlonlatBk.first;
	double gcjlatBk = gcjlonlatBk.second;
	std::pair<double, double> gcjlonlatFw = Wgs2Gcj(wgs84lon, latFw);
	double gcjlonFw = gcjlonlatFw.first;
	double gcjlatFw = gcjlonlatFw.second;
	double dlongcj_dlatwgs = (gcjlonBk - gcjlonFw) / (dlat * 2.0);
	double dlatgcj_dlatwgs = (gcjlatBk - gcjlatFw) / (dlat * 2.0);
	return { dlongcj_dlatwgs , dlatgcj_dlatwgs };
}




WGSTOGCJ_API_EXPORTS std::pair<double, double> Gcj2Wgs_NumbericDiff(const double& gcj02lon, const double& gcj02lat)
{
	
	double wgs84lon = gcj02lon, wgs84lat = gcj02lat;
	int nIterCount = 0;
	double tol = 1e-9;
	while (++nIterCount < 1000)
	{
		std::pair<double, double > data1 = GetPartialDerivative_Lon(wgs84lon, wgs84lat, tol);
		double dlongcj_dlonwgs = data1.first;
		double dlatgcj_dlonwgs = data1.second;
		std::pair<double, double > data2 = GetPartialDerivative_Lat(wgs84lon, wgs84lat, tol);
		double dlongcj_dlatwgs = data2.first;
		double dlatgcj_dlatwgs = data2.second;
		std::pair<double, double > data3 = Wgs2Gcj(wgs84lon, wgs84lat);
		double gcj02lonEst = data3.first;
		double gcj02latEst = data3.second;
		double l_lon = gcj02lon - gcj02lonEst;
		double l_lat = gcj02lat - gcj02latEst;
		double d_latwgs = (l_lon * dlatgcj_dlonwgs - l_lat * dlongcj_dlonwgs) /
			(dlongcj_dlatwgs * dlatgcj_dlonwgs - dlatgcj_dlatwgs * dlongcj_dlonwgs);
		double d_lonwgs = (l_lon - dlongcj_dlatwgs * d_latwgs) / dlongcj_dlonwgs;

		if (fabs(d_latwgs) < tol && fabs(d_lonwgs) < tol)
			break;
		wgs84lon = wgs84lon + d_lonwgs;
		wgs84lat = wgs84lat + d_latwgs;
	}
	return { wgs84lon, wgs84lat };
		

}

WGSTOGCJ_API_EXPORTS std::pair<double, double> Gcj2Wgs_AnalyticDiff(const double& gcj02lon, const double& gcj02lat)
{
	double wgs84lon = gcj02lon, wgs84lat = gcj02lat;
	int nIterCount = 0;
	while (++nIterCount < 1000)
	{
		//get geodetic offset relative to 'center china'
		double lon0 = wgs84lon - 105.0;
		double lat0 = wgs84lat - 35.0;

		//generate an pair offset roughly in meters
		double lon1 = 300.0 + lon0 + 2.0 * lat0 + 0.1 * lon0 * lon0 + 0.1 * lon0 * lat0 + 0.1 * sqrt(fabs(lon0));
		lon1 = lon1 + (20.0 * sin(6.0 * lon0 * PI) + 20.0 * sin(2.0 * lon0 * PI)) * 2.0 / 3.0;
		lon1 = lon1 + (20.0 * sin(lon0 * PI) + 40.0 * sin(lon0 / 3.0 * PI)) * 2.0 / 3.0;
		lon1 = lon1 + (150.0 * sin(lon0 / 12.0 * PI) + 300.0 * sin(lon0 * PI / 30.0)) * 2.0 / 3.0;
		double lat1 = -100.0 + 2.0 * lon0 + 3.0 * lat0 + 0.2 * lat0 * lat0 + 0.1 * lon0 * lat0 + 0.2 * sqrt(fabs(lon0));
		lat1 = lat1 + (20.0 * sin(6.0 * lon0 * PI) + 20.0 * sin(2.0 * lon0 * PI)) * 2.0 / 3.0;
		lat1 = lat1 + (20.0 * sin(lat0 * PI) + 40.0 * sin(lat0 / 3.0 * PI)) * 2.0 / 3.0;
		lat1 = lat1 + (160.0 * sin(lat0 / 12.0 * PI) + 320.0 * sin(lat0 * PI / 30.0)) * 2.0 / 3.0;

		double g_lon0 = 0;
		if (lon0 > 0)
			g_lon0 = 0.05 / sqrt(lon0);
		else
			if (lon0 < 0)
				g_lon0 = -0.05 / sqrt(-lon0);
			else
				g_lon0 = 0;

		double PIlon0 = PI * lon0, PIlat0 = PI * lat0;
		double dlon1_dlonwgs = 1 + 0.2 * lon0 + 0.1 * lat0 + g_lon0
			+ ((120 * PI * cos(6 * PIlon0) + 40 * PI * cos(2 * PIlon0))
				+ (20 * PI * cos(PIlon0) + 40 * PI / 3.0 * cos(PIlon0 / 3.0))
				+ (12.5 * PI * cos(PIlon0 / 12.0) + 10 * PI * cos(PIlon0 / 30.0))) * 2.0 / 3.0;
		double dlon1_dlatwgs = 2 + 0.1 * lon0;

		double dlat1_dlonwgs = 2 + 0.1 * lat0 + 2 * g_lon0
			+ (120 * PI * cos(6 * PIlon0) + 40 * PI * cos(2 * PIlon0)) * 2.0 / 3.0;
		double dlat1_dlatwgs = 3 + 0.4 * lat0 + 0.1 * lon0
			+ ((20 * PI * cos(PIlat0) + 40.0 * PI / 3.0 * cos(PIlat0 / 3.0))
				+ (40 * PI / 3.0 * cos(PIlat0 / 12.0) + 32.0 * PI / 3.0 * cos(PIlat0 / 30.0))) * 2.0 / 3.0;

		//latitude in radian
		double B = Deg2Rad(wgs84lat);
		double sinB = sin(B), cosB = cos(B);
		double WSQ = 1 - kKRASOVSKY_ECCSQ * sinB * sinB;
		double W = sqrt(WSQ);
		double N = kKRASOVSKY_A / W;

		double dW_dlatwgs = -PI * kKRASOVSKY_ECCSQ * sinB * cosB / (180.0 * W);
		double dN_dlatwgs = -kKRASOVSKY_A * dW_dlatwgs / WSQ;

		double PIxNxCosB = PI * N * cosB;
		double dlongcj_dlonwgs = 1.0 + 180.0 * dlon1_dlonwgs / PIxNxCosB;
		double dlongcj_dlatwgs = 180 * dlon1_dlatwgs / PIxNxCosB -
			180 * lon1 * PI * (dN_dlatwgs * cosB - PI * N * sinB / 180.0) / (PIxNxCosB * PIxNxCosB);

		double PIxNxSubECCSQ = PI * N * (1 - kKRASOVSKY_ECCSQ);
		double dlatgcj_dlonwgs = 180 * WSQ * dlat1_dlonwgs / PIxNxSubECCSQ;
		double dlatgcj_dlatwgs = 1.0 + 180 * (N * (dlat1_dlatwgs * WSQ + 2.0 * lat1 * W * dW_dlatwgs) - lat1 * WSQ * dN_dlatwgs) /
			(N * PIxNxSubECCSQ);

		std::pair<double, double> gcj02lonlatEst = Wgs2Gcj(wgs84lon, wgs84lat);
		double gcj02lonEst = gcj02lonlatEst.first;
		double gcj02latEst = gcj02lonlatEst.second;
		double l_lon = gcj02lon - gcj02lonEst;
		double l_lat = gcj02lat - gcj02latEst;

		double d_latwgs = (l_lon * dlatgcj_dlonwgs - l_lat * dlongcj_dlonwgs) /
			(dlongcj_dlatwgs * dlatgcj_dlonwgs - dlatgcj_dlatwgs * dlongcj_dlonwgs);
		double d_lonwgs = (l_lon - dlongcj_dlatwgs * d_latwgs) / dlongcj_dlonwgs;

		if (fabs(d_latwgs) < 1.0e-9 && fabs(d_lonwgs) < 1.0e-9)
			break;
		wgs84lon = wgs84lon + d_lonwgs;
		wgs84lat = wgs84lat + d_latwgs;
	}
	return { wgs84lon, wgs84lat };
	

}

void Test(double lon, double lat)
{
	double lonWgs = lon;
	double latWgs = lat;
	std::cout.precision(16);
	std::cout << "WGS84 Point: (" << latWgs << ", " << lonWgs << ")\n";


	std::pair<double, double> lonlatGcj = Wgs2Gcj(lonWgs, latWgs);
	double lonGcj = lonlatGcj.first;
	double latGcj = lonlatGcj.second;
	std::cout << "GCJ-02 Point [Wgs2Gcj]: (" << latGcj << ", " << lonGcj << ")\n";


	//simple linear iteration 
	std::pair<double, double> lonlatWgsNS = Gcj2Wgs_SimpleIteration(lonGcj, latGcj);
	double lonWgsNS = lonlatWgsNS.first;
	double latWgsNS = lonlatWgsNS.second;
	std::cout << "WGS84 Point [simple linear iteration]: (" << latWgsNS << ", " << lonWgsNS << ")\n";

	//numerically differentiated cost function

	std::pair<double, double> lonlatWgsND = Gcj2Wgs_NumbericDiff(lonGcj, latGcj);
	double lonWgsND = lonlatWgsND.first;
	double latWgsND = lonlatWgsND.second;
	std::cout << "WGS84 Point [numerically derivation]: (" << latWgsND << ", " << lonWgsND << ")\n";

	//analytical differentiated cost function

	std::pair<double, double> lonlatWgsNAn = Gcj2Wgs_AnalyticDiff(lonGcj, latGcj);
	double lonWgsNAn = lonlatWgsNAn.first;
	double latWgsNAn = lonlatWgsNAn.second;
	std::cout << "WGS84 Point [analytical derivation]: (" << latWgsNAn << ", " << lonWgsNAn << ")\n";

	//auto differentiable, use ceres
	std::pair<double, double> lonlatWgsNA = Gcj2Wgs_AutoDiff(lonGcj, latGcj);
	double lonWgsNA = lonlatWgsNA.first;
	double latWgsNA = lonlatWgsNA.second;
	std::cout << "WGS84 Point [auto differentiable]: (" << latWgsNA << ", " << lonWgsNA << ")\n";

	
}

std::pair<double, double> Wsg2Gcj(double Wlon, double Wlat)
{
	double lonWgs = Wlon;
	double latWgs = Wlat;

	std::cout.precision(16);
	std::cout << "Ô­Ê¼WGS84 Point: (" << lonWgs << ", " << latWgs << ")\n";


	std::pair<double, double> lonlatGcj = Wgs2Gcj(lonWgs, latWgs);
	double lonGcj = lonlatGcj.first;
	double latGcj = lonlatGcj.second;
	std::cout << "GCJ-02 Point [Wgs2Gcj]: (" << lonGcj << ", " << latGcj << ")\n";

	return { lonGcj, latGcj };

}



std::pair<double,double> Gcj2Wsg(double Glon, double Glat)
{
	double lonGcj = Glon;
	double latGcj = Glat;



	////simple linear iteration 
	//std::pair<double, double> lonlatWgsNS = Gcj2Wgs_SimpleIteration(lonGcj, latGcj);
	//double lonWgsNS = lonlatWgsNS.first;
	//double latWgsNS = lonlatWgsNS.second;
	//std::cout << "WGS84 Point [simple linear iteration]: (" << latWgsNS << ", " << lonWgsNS << ")\n";

	////numerically differentiated cost function

	//std::pair<double, double> lonlatWgsND = Gcj2Wgs_NumbericDiff(lonGcj, latGcj);
	//double lonWgsND = lonlatWgsND.first;
	//double latWgsND = lonlatWgsND.second;
	//std::cout << "WGS84 Point [numerically derivation]: (" << latWgsND << ", " << lonWgsND << ")\n";

	////analytical differentiated cost function

	//std::pair<double, double> lonlatWgsNAn = Gcj2Wgs_AnalyticDiff(lonGcj, latGcj);
	//double lonWgsNAn = lonlatWgsNAn.first;
	//double latWgsNAn = lonlatWgsNAn.second;
	//std::cout << "WGS84 Point [analytical derivation]: (" << latWgsNAn << ", " << lonWgsNAn << ")\n";



	//auto differentiable, use ceres
	std::pair<double, double> lonlatWgsNA = Gcj2Wgs_AutoDiff(lonGcj, latGcj);
	double lonWgsNA = lonlatWgsNA.first;
	double latWgsNA = lonlatWgsNA.second;
	std::cout << "WGS84 Point [auto differentiable]:===== (" << lonWgsNA << ", " << latWgsNA << ")\n";

	
	return { lonWgsNA, latWgsNA };

}



