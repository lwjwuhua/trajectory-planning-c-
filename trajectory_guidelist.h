#pragma once
#ifndef TRAJECTORY_GUIDELIST_H_
#define TRAJECTORY_GUIDELIST_H_

#include "eigen3/Eigen/Dense"




class Path_list
{
public:

	int path_step_id;
	Eigen::VectorXd path_list;
	vector <VectorXd> path_step_trajectory;

private:

};




class Trajectory_withlane
{
public:
	
	Eigen::VectorXd traject_poslonlat;
	 int laneid;
	 int path_step_id;

	 Path_list traject_poslonlat_withid;

private:

};




class Trajectory_list
{
public:
	
	int rsu_secq;
	int rsu_id;
	Eigen::VectorXd Rsu_ori_des_list;
	
	vector <Trajectory_withlane> tarject_withlane_list;

	
private:

};










#endif // !TRAJECTORY_GUIDELIST_H_
