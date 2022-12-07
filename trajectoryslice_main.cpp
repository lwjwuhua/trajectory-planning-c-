#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include "Eigen/Dense"
#include"gpscovert.h"
#include"Polynomial_insert.h"
#include <cmath>
#include <iomanip>
#include <iterator>
// #include <ceres/ceres.h>
#include"WGS84GCJ02.h"
#include"trajectory_guidelist.h"
#include"RSU.h"  
#include"Spline.h"

using Eigen::MatrixXd;
using Eigen::VectorXf;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::string;
using std::vector;
using namespace SplineSpace;

int main() {

	GpsTo gps2;
 	PolyNomial polynomi;

	string out_file_name_ = "4-output.txt";
	


	/* ========================================== 前处理  =====================================================*/



	/* 1,读取txt文件里的高德坐标点 GCJ-02坐标系 */

	string in_file_name_ = "2.txt";
	ifstream in_file(in_file_name_.c_str(), ifstream::in);
	ofstream out_file(out_file_name_.c_str(), ofstream::out);

	if (!in_file.is_open()) {
		cout << "Cannot open input file: " << in_file_name_ << endl;
	}
    string line;
    int i = 0;

	VectorXd gd_pos (2);
	vector <VectorXd> gd_pos_list;

	while (getline(in_file, line) && (i <= 600)) {


		istringstream iss(line);
		double pos_x;
		iss >> pos_x; // reads first element from the current line

		char c;
		iss >> c;

		double pos_y;
		iss >> pos_y;

		gd_pos << pos_x, pos_y;
	 	gd_pos_list.push_back(gd_pos);

	}

     
	for (i = 0; i < gd_pos_list.size(); i++) {
		cout << fixed << setprecision(6)<< gd_pos_list[i][0]<<"\t"<< gd_pos_list[i][1] << endl;
	}

	/* 实例化RSU，并添加相应属性 */
	RSU rsu;

	vector<RSU> rsu_;
	rsu.id = 3;
	rsu.lon = 106.4840548;   rsu.lat = 29.6635501;
	rsu.r = 500;
	rsu_.push_back(rsu);
	rsu.id = 2;
	rsu.lon = 106.4841978;   rsu.lat = 29.6619316;
	rsu.r = 500;
	rsu_.push_back(rsu);
	rsu.id = 1;
	rsu.lon = 106.4850102;   rsu.lat = 29.6618484;
	rsu.r = 500;
	rsu_.push_back(rsu);


//	Trajectory_list trajec_list_; //轨迹引导列表
	vector<Trajectory_list> trajec_list;
	Path_list path_list;  
	vector <Path_list> path_step;  //路段轨迹列表



	/* ==========================分段插值添加， 逐个读取每一段的经纬轨迹，并分别插值==============================  */
	int N = gd_pos_list.size();
	
	int count = 0;
	int n=0;
	vector<VectorXd> Output;

	for (size_t K = 0; K < N; K++)         
	{
		Trajectory_list trajec_list_;
		int i = K;
		int e = 0;
		vector<VectorXd> Output_temp;
		Output_temp.clear();

		VectorXd gd_pos_(2);
		vector <VectorXd> gd_pos_list_;
		gd_pos_list_.clear();

		for ( i ; (e < 200) && (i < N); i++)
		{

			
			if (i != K)
			{
				e = gd_pos_list[i][0] - gd_pos_list[i - 1][0];
				gd_pos_ << gd_pos_list[i - 1][0], gd_pos_list[i - 1][1];
				gd_pos_list_.push_back(gd_pos_);   // 用于后续轨迹插值

			}

		

		}

		cout << "========================第" << n << "段路的高德经纬轨迹=========================" << endl;
		for (size_t i = 0; i < gd_pos_list_.size(); i++)
		{
			cout << "第" << i << "个点的高德经纬轨迹====" << fixed << setprecision(6)<< gd_pos_list_[i] << endl;
		
		}

		/*  ======================================   对该段进行轨迹插值切片  ======================================   */


 

		/* 2 ，坐标投影 */


		VectorXd proj_pos(2);
		vector <VectorXd> proj_pos_list;


		for (size_t i = 0; i < gd_pos_list_.size(); i++)
		{
			double x = gd_pos_list_[i][0];
			double y = gd_pos_list_[i][1];
			double z = 100;

			printf("原大地经纬度坐标：%.10lf\t%.10lf\t%.10lf\n", x, y, z);


			gps2.Blh2Wmc(x, y, z);
			printf("Web墨卡托坐标：%.10lf\t%.10lf\t%.10lf\n", x, y, z);

			//	cout << fixed << setprecision(6) << x << endl;
			proj_pos << x, y;
			proj_pos_list.push_back(proj_pos);

			/*	gps2.Wmc2Blh(x, y, z);
				printf("转回大地经纬度坐标：%.10lf\t%.10lf\t%.10lf\n", x, y, z);*/

		}

		for (i = 0; i < proj_pos_list.size(); i++) {
			cout << fixed << setprecision(6) << proj_pos_list[i][0] << "\t";

			cout << fixed << setprecision(6) << proj_pos_list[i][1] << "\n";
		}


		/* 3, 轨迹插值切片 */

		// vector<VectorXd> Output;
	   
		double x[10];	//插值点
		double y[10];
		VectorXd proj_temp(2);
		if (  proj_pos_list.size()<4  )   
		{
			VectorXd proj_pos_list_temp_x(proj_pos_list.size());
			VectorXd proj_pos_list_temp_y(proj_pos_list.size());
			
			double x0[100];		//已知的数据点
			double y0[100] ;

			for (size_t l = 0; l < proj_pos_list.size(); l++)     //已知的数据点
			{

				x0[l] = proj_pos_list[l][0];
				y0[l]=  proj_pos_list[l][1];
			}

			//for (size_t l_ = 0; l_ < size(proj_pos_list); l_++)     //已知的数据点
			//{
			//	
			//	cout<<"==========" << x0[l_] << "\t" << y0[l_] << endl;

			//}
			cout << "===================================" << endl;

			//double x[10];	//插值点
			//double y[10];
			// VectorXd proj_temp(2);

			try
			{

				SplineInterface* sp = new Spline(x0, y0, proj_pos_list.size());	//使用接口，且使用默认边界条件
				sp->AutoInterp(10, x, y);			//求x的插值结果y

				for (int i = 0; i < 10; i++) {
					cout << x[i] << "\t" << y[i] << endl;
				/*	proj_temp << x[i], y[i];
					Output.push_back(proj_temp);*/

				}


			}
			catch (SplineFailure sf)
			{
			//	cout << sf.GetMessage() << endl;
			}

		}
		else {
		
			Output_temp = polynomi.catMullRomSpline(proj_pos_list);
		
		}

	


		/*数组倒序,高德序列与插值序列相反的时候需要倒序*/
		if (proj_pos_list[0][0] != x[0])     // 数组倒序,高德序列与插值序列相反的时候需要倒序
		{


			int i = 0;  //循环变量1, i的值为数组第一个元素的下标

			int j = ( sizeof(x)/sizeof(x[0]) ) - 1;  //循环变量2, j的值为数组最后一个元素的下标

			double buf;  //互换时的中间存储变量

			for (; i < j; ++i, --j)  /*因为i和j已经初始化过了, 所以表达式1可以省略, 但表达式1后面的分号不能省。*/

			{

				buf = x[i];

				x[i] = x[j];

				x[j] = buf;

			}
			int i_ = 0;  

			int j_ = ( sizeof(y)/sizeof(y[0]) ) - 1;
			double buf_;
			for (; i_ < j_; ++i_, --j_)

			{
				

				buf_ = y[i_];

				y[i_] = y[j_];

				y[j_] = buf_;

			}



		}

	

		for (int i = 0; i < ( sizeof(x)/sizeof(x[0]) ); i++) {
			cout << x[i] << "//t" << y[i] << endl;
			/*	proj_temp << x[i], y[i];
				Output.push_back(proj_temp);*/

		}




		/* 4， 转换回经纬坐标  */


		for (size_t i = 0; i < ( sizeof(x)/sizeof(x[0]) ); i++)
		{
			double z_ = 100;

			gps2.Wmc2Blh(x[i], y[i], z_);

			gd_pos_ << x[i], y[i];
			Output_temp.push_back(gd_pos_);   //当前路段轨迹序列

			proj_temp << x[i], y[i];
			Output.push_back(proj_temp);
			
		}




		/* 5，  GCJ-02坐标系 转换gps坐标  */



		/* ========================================== 路径分发  =====================================================*/

		/* 6，路径分发      */

		 /* 6.1，创建轨迹引导列表， 通过高德坐标点选取RSU并排序   */

		//Trajectory_list trajec_list_;
		//vector<Trajectory_list> trajec_list;

		///* 实例化RSU，并添加相应属性 */
		//RSU rsu;

		///*;
		// rsu0.id = 0;
		// rsu0.lon = 106.4840548;   rsu0.lat = 29.6635501;  //34
		// rsu1.id = 1;
		// rsu1.lon = 106.4841978;   rsu1.lat = 29.6619316; //1109
		// rsu2.id = 2;
		// rsu2.lon = 106.4850102;   rsu2.lat = 29.6618484; //1639
		// */

		//vector<RSU> rsu_;
		//rsu.id = 0;
		//rsu.lon = 106.4840548;   rsu.lat = 29.6635501;
		//rsu.r = 500;
		//rsu_.push_back(rsu);
		//rsu.id = 1;
		//rsu.lon = 106.4841978;   rsu.lat = 29.6619316;
		//rsu.r = 500;
		//rsu_.push_back(rsu);
		//rsu.id = 2;
		//rsu.lon = 106.4850102;   rsu.lat = 29.6618484;
		//rsu.r = 500;
		//rsu_.push_back(rsu);
 


		/* 6.2,   轨迹按路段切片分配给对应RSU */

	   //  int traject_secqcount = trajec_list.size();

	//	Path_list path_list;

	//	vector <Path_list> path_step;  //路段轨迹列表

			path_list.path_step_id = n;
			path_list.path_step_trajectory = Output_temp;
			
			path_step.push_back(path_list);   //带有当前路段编号的轨迹序列


		
		//cout << "recent_path_pos_count=====" << size(path_step) << endl;


		//for (size_t p = 0; p < size(path_step); p++)
		//{


		//	for (size_t p_ = 0; p_ < size(path_step[p].path_step_trajectory); p_++)
		//	{

		//		cout << "recent_path_id===" << path_step[p].path_step_id << '\n' << "recent_path_pos_list======" << path_step[p].path_step_trajectory[p_] << endl;
		//	}

		//}

			cout << "recent_path_id===" << path_step[n].path_step_id << endl;
			for (size_t p_ = 0; p_ < path_step[n].path_step_trajectory.size(); p_++)
			{
				cout  << "recent_path_pos_list======" << path_step[n].path_step_trajectory[p_] << endl;

			}


		/* 判断当前段轨迹所属RSU  (目前暂时用每个路段的起始坐标来确立RSU)       */
		int path_count = path_step.size();

		double z_=100;
		gps2.Blh2Wmc(rsu_[0].lon, rsu_[0].lat, z_);
		gps2.Blh2Wmc(rsu_[1].lon, rsu_[1].lat, z_);
		gps2.Blh2Wmc(rsu_[2].lon, rsu_[2].lat, z_);

		cout << "b2w====" << rsu_[0].lon << '\t' << rsu_[0].lat << endl;
		cout << "b2w====" << rsu_[1].lon << '\t' << rsu_[1].lat << endl;
		cout << "b2w====" << rsu_[2].lon << '\t' << rsu_[2].lat << endl;

		
		int traject_secqcount = trajec_list.size();

		double tra_lon = path_step[n].path_step_trajectory[0][0];
		double tra_lat = path_step[n].path_step_trajectory[0][1];
		gps2.Blh2Wmc(tra_lon, tra_lat, z_);

		double r_e;

		// trajec_list_.rsu_secq = k;

		bool recentpath_step_in_rsu = true;
		for (size_t i = 0; i < rsu_.size(); i++)
		{

			if (recentpath_step_in_rsu)
			{
				r_e = sqrt((rsu_[i].lon - tra_lon) * (rsu_[i].lon - tra_lon) + (rsu_[i].lat - tra_lat) * (rsu_[i].lat - tra_lat));
				cout << "r_e====" << r_e << endl;

				if (r_e < 500)   // 判断在哪个RSU的覆盖半径内 ,如果在，就把rsu id 添加进引导列表
				{

					trajec_list_.rsu_id = rsu_[i].id;

					recentpath_step_in_rsu = false;

				}
				else
				{
					cout << "不在RSU范围内" << endl;
				}

			}



		}

		gps2.Wmc2Blh(rsu_[0].lon, rsu_[0].lat, z_);
		gps2.Wmc2Blh(rsu_[1].lon, rsu_[1].lat, z_);
		gps2.Wmc2Blh(rsu_[2].lon, rsu_[2].lat, z_);


		trajec_list_.rsu_secq = traject_secqcount;
		cout << "trajeclist_count======" << traject_secqcount << endl;

		if ((trajec_list.size() == 0) || (trajec_list_.rsu_id != trajec_list[traject_secqcount - 1].rsu_id))
		{


			//for (size_t ii = 0; ii < size(path_step[n].path_step_trajectory); ii++)
			//{
			//	Trajectory_withlane a;
			//	a.traject_poslonlat = path_step[n].path_step_trajectory[ii];
			//	trajec_list_.tarject_withlane_list.push_back(a);

			//}

			Trajectory_withlane b;
			b.traject_poslonlat_withid = path_step[n];
			b.path_step_id = path_step[n].path_step_id;
			trajec_list_.tarject_withlane_list.push_back(b);
			trajec_list.push_back(trajec_list_);


		}
		else
		{
			/* 仅添加经纬序列  */

			//for (size_t ii = 0; ii < size(path_step[n].path_step_trajectory); ii++)
			//{
			//	Trajectory_withlane c;
			//	c.path_step_id = n;
			//	c.traject_poslonlat = path_step[n].path_step_trajectory[ii];

			//

			//	trajec_list[traject_secqcount - 1].tarject_withlane_list.push_back(c);

			//}

			Trajectory_withlane c;
			c.path_step_id = n;
			
			c.traject_poslonlat_withid = path_step[n];

			trajec_list[traject_secqcount - 1].tarject_withlane_list.push_back(c);

		}





		
		n = n + 1;
		K = K + gd_pos_list_.size();

	}


	
	cout << "==============================列表结果=====================================" << endl;

	for (size_t i = 0; i < trajec_list.size(); i++)
	{

		cout << "+++++++++++++++++++++++++++++++++++++++" << endl;
		cout << "traclist_id=====" << trajec_list[i].rsu_secq << "\t" << "rsu_id=======" << trajec_list[i].rsu_id << endl;

		
	//	cout << "tarject_withlane_list.size===" << trajec_list[i].tarject_withlane_list.size() << endl;

		for (size_t i_ = 0; i_ < trajec_list[i].tarject_withlane_list.size(); i_++)
		{
			cout << "path_step_id===" << trajec_list[i].tarject_withlane_list[i_].path_step_id << endl;

			for (size_t p = 0; p < trajec_list[i].tarject_withlane_list[i_].traject_poslonlat_withid.path_step_trajectory.size(); p++)
			{
				cout << "tarject_poslonlat====" << trajec_list[i].tarject_withlane_list[i_].traject_poslonlat_withid.path_step_trajectory[p][0]
					<<'\t'<< trajec_list[i].tarject_withlane_list[i_].traject_poslonlat_withid.path_step_trajectory[p][1] << endl;
				
				/* 用于最后结果回放 */
				gd_pos << trajec_list[i].tarject_withlane_list[i_].traject_poslonlat_withid.path_step_trajectory[p][0], trajec_list[i].tarject_withlane_list[i_].traject_poslonlat_withid.path_step_trajectory[p][1];
				Output.push_back(gd_pos);

			}


		}

		
	}

	


	/* ============================结果输出 ============================= */

	

	Output.erase(unique(Output.begin(), Output.end()), Output.end());
	for (size_t i = 0; i < Output.size(); i++)
	{
		cout << Output[i][0] << "=====" << Output[i][1] << endl;

	}















	 /*   结果导出  */

	for (size_t k = 0; k < Output.size(); ++k) {

		// out_file << "[" ;
		// out_file << fixed << setprecision(6)  << Output[k][0] << "\t";
		// out_file << fixed << setprecision(6) << Output[k][1] << "\n";
		out_file << fixed << setprecision(6) <<"["<< Output[k][0]<<","<<Output[k][1]<<"]" <<"," << "\n";
		//out_file << fixed << setprecision(6)  << Output[k][0] << "\t" << Output[k][1]   << "\n";

	}


	if (in_file.is_open()) {
		in_file.close();
	}


	return 0;
}
