#include "Polynomial_insert.h"
#include <string>
using namespace std;
using std::vector;

PolyNomial::PolyNomial()
{
}

PolyNomial::~PolyNomial()
{
}


double tj(double ti, VectorXd Pi, VectorXd Pj) {
    double t = sqrt(sqrt(pow((Pj[0] - Pi[0]), 2) + pow((Pj[1] - Pi[1]), 2))) + ti;
    return t;
}



vector<VectorXd> PolyNomial::catMullRomSpline(vector<VectorXd> inputPoints)
{
    // ÿ�������м������ٸ���
    int numSpace = 232;

    int numPoints = inputPoints.size();
    cout <<"num =="<< numPoints << endl;

    VectorXd c(2);
    vector<VectorXd> C;

    if (numPoints ==4 )
    {
        /*  P0-P1֮��Ĳ�ֵ  */
        VectorXd P0(2), P1(2), P2(2), P3(2), P_temp(2), P_temp_(2);

        P0[0] = inputPoints[0][0];
        P0[1] = inputPoints[0][1];
        P1[0] = inputPoints[1][0];
        P1[1] = inputPoints[1][1];
        P2[0] = inputPoints[2][0];
        P2[1] = inputPoints[2][1];
        P3[0] = inputPoints[3][0];
        P3[1] = inputPoints[3][1];

        P_temp_ = P3;
        P_temp = 2 * P0 - P1;

        P3 = P2;
        P2 = P1;
        P1 = P0;
        P0 = P_temp;


        double t0 = 0;
        double t1 = tj(t0, P0, P1);
        double t2 = tj(t1, P1, P2);
        double t3 = tj(t2, P2, P3);

        //  cout << "t2-t1==="<<t2-t1 << endl;

         // �������Ϊ�����֮��ļ��
        double linespace = (t2 - t1) / numSpace;
        // cout << linespace << endl;
        double t = t1;


        VectorXd c1(2);
        vector<VectorXd> C2;
        while (t <= t2) {
            double A1_x = (t1 - t) / (t1 - t0) * P0[0] + (t - t0) / (t1 - t0) * P1[0];
            double A1_y = (t1 - t) / (t1 - t0) * P0[1] + (t - t0) / (t1 - t0) * P1[1];
            double A2_x = (t2 - t) / (t2 - t1) * P1[0] + (t - t1) / (t2 - t1) * P2[0];
            double A2_y = (t2 - t) / (t2 - t1) * P1[1] + (t - t1) / (t2 - t1) * P2[1];
            double A3_x = (t3 - t) / (t3 - t2) * P2[0] + (t - t2) / (t3 - t2) * P3[0];
            double A3_y = (t3 - t) / (t3 - t2) * P2[1] + (t - t2) / (t3 - t2) * P3[1];
            double B1_x = (t2 - t) / (t2 - t0) * A1_x + (t - t0) / (t2 - t0) * A2_x;
            double B1_y = (t2 - t) / (t2 - t0) * A1_y + (t - t0) / (t2 - t0) * A2_y;
            double B2_x = (t3 - t) / (t3 - t1) * A2_x + (t - t1) / (t3 - t1) * A3_x;
            double B2_y = (t3 - t) / (t3 - t1) * A2_y + (t - t1) / (t3 - t1) * A3_y;
            double C_x = (t2 - t) / (t2 - t1) * B1_x + (t - t1) / (t2 - t1) * B2_x;
            double C_y = (t2 - t) / (t2 - t1) * B1_y + (t - t1) / (t2 - t1) * B2_y;
            /* C_x = floor(C_x);
             C_y = floor(C_y);*/


             //c2 << C_x, C_y;
             //C2.push_back(c2);

            c1 << C_x, C_y;
            C.push_back(c1);

            t = t + linespace;
        }
        // C ������4�����У���2����͵�3����֮��Ĳ�ֵ��




        /*  P1-P2֮��Ĳ�ֵ  */

        P0 = P1;
        P1 = P2;
        P2 = P3;
        P3 = P_temp_;



        t0 = 0;
        t1 = tj(t0, P0, P1);
        t2 = tj(t1, P1, P2);
        t3 = tj(t2, P2, P3);
        linespace = (t2 - t1) / numSpace;
        t = t1;


        VectorXd c2(2);
        vector<VectorXd> C1;

        while (t <= t2) {
            double A1_x = (t1 - t) / (t1 - t0) * P0[0] + (t - t0) / (t1 - t0) * P1[0];
            double A1_y = (t1 - t) / (t1 - t0) * P0[1] + (t - t0) / (t1 - t0) * P1[1];
            double A2_x = (t2 - t) / (t2 - t1) * P1[0] + (t - t1) / (t2 - t1) * P2[0];
            double A2_y = (t2 - t) / (t2 - t1) * P1[1] + (t - t1) / (t2 - t1) * P2[1];
            double A3_x = (t3 - t) / (t3 - t2) * P2[0] + (t - t2) / (t3 - t2) * P3[0];
            double A3_y = (t3 - t) / (t3 - t2) * P2[1] + (t - t2) / (t3 - t2) * P3[1];
            double B1_x = (t2 - t) / (t2 - t0) * A1_x + (t - t0) / (t2 - t0) * A2_x;
            double B1_y = (t2 - t) / (t2 - t0) * A1_y + (t - t0) / (t2 - t0) * A2_y;
            double B2_x = (t3 - t) / (t3 - t1) * A2_x + (t - t1) / (t3 - t1) * A3_x;
            double B2_y = (t3 - t) / (t3 - t1) * A2_y + (t - t1) / (t3 - t1) * A3_y;
            double C_x = (t2 - t) / (t2 - t1) * B1_x + (t - t1) / (t2 - t1) * B2_x;
            double C_y = (t2 - t) / (t2 - t1) * B1_y + (t - t1) / (t2 - t1) * B2_y;

            //c1 << C_x, C_y;
            //C1.push_back(c1);


            c2 << C_x, C_y;
            C.push_back(c2);

            t = t + linespace;

        }


        /*  P2-P3֮��Ĳ�ֵ  */

        P_temp_ = P3;
        P_temp = 2 * P3 - P2;


        P0 = P1;
        P1 = P2;
        P2 = P3;
        P3 = P_temp;


        t0 = 0;
        t1 = tj(t0, P0, P1);
        t2 = tj(t1, P1, P2);
        t3 = tj(t2, P2, P3);
        linespace = (t2 - t1) / numSpace;
        t = t1;


        VectorXd c3(2);


        while (t <= t2) {
            double A1_x = (t1 - t) / (t1 - t0) * P0[0] + (t - t0) / (t1 - t0) * P1[0];
            double A1_y = (t1 - t) / (t1 - t0) * P0[1] + (t - t0) / (t1 - t0) * P1[1];
            double A2_x = (t2 - t) / (t2 - t1) * P1[0] + (t - t1) / (t2 - t1) * P2[0];
            double A2_y = (t2 - t) / (t2 - t1) * P1[1] + (t - t1) / (t2 - t1) * P2[1];
            double A3_x = (t3 - t) / (t3 - t2) * P2[0] + (t - t2) / (t3 - t2) * P3[0];
            double A3_y = (t3 - t) / (t3 - t2) * P2[1] + (t - t2) / (t3 - t2) * P3[1];
            double B1_x = (t2 - t) / (t2 - t0) * A1_x + (t - t0) / (t2 - t0) * A2_x;
            double B1_y = (t2 - t) / (t2 - t0) * A1_y + (t - t0) / (t2 - t0) * A2_y;
            double B2_x = (t3 - t) / (t3 - t1) * A2_x + (t - t1) / (t3 - t1) * A3_x;
            double B2_y = (t3 - t) / (t3 - t1) * A2_y + (t - t1) / (t3 - t1) * A3_y;
            double C_x = (t2 - t) / (t2 - t1) * B1_x + (t - t1) / (t2 - t1) * B2_x;
            double C_y = (t2 - t) / (t2 - t1) * B1_y + (t - t1) / (t2 - t1) * B2_y;

            //c1 << C_x, C_y;
            //C1.push_back(c1);


            c3 << C_x, C_y;
            C.push_back(c3);

            t = t + linespace;

        }

        return C;

    }
   
     else if (numPoints > 4) {


        int numSplines = numPoints - 3;

        std::vector<VectorXd> result;
        std::vector<VectorXd> spline;

        int last = inputPoints.size()-1;
	

        VectorXd P0(2), P1(2), P2(2), P3(2), P_temp(2), P_temp_(2);

        /* ��һ����͵ڶ�����Ĳ�ֵ  */

        P0[0] = inputPoints[0][0];
        P0[1] = inputPoints[0][1];
        P1[0] = inputPoints[1][0];
        P1[1] = inputPoints[1][1];
        P2[0] = inputPoints[2][0];
        P2[1] = inputPoints[2][1];
        P3[0] = inputPoints[3][0];
        P3[1] = inputPoints[3][1];

        P_temp_ = P3;
        P_temp = 2 * P0 - P1;

        P3 = P2;
        P2 = P1;
        P1 = P0;
        P0 = P_temp;


        double t0 = 0;
        double t1 = tj(t0, P0, P1);
        double t2 = tj(t1, P1, P2);
        double t3 = tj(t2, P2, P3);

        //  cout << "t2-t1==="<<t2-t1 << endl;

         // �������Ϊ�����֮��ļ��
        double linespace = (t2 - t1) / numSpace;
        // cout << linespace << endl;
         double t = t1;


        VectorXd c1(2);
        vector<VectorXd> C2;
        while (t <= t2) {
            double A1_x = (t1 - t) / (t1 - t0) * P0[0] + (t - t0) / (t1 - t0) * P1[0];
            double A1_y = (t1 - t) / (t1 - t0) * P0[1] + (t - t0) / (t1 - t0) * P1[1];
            double A2_x = (t2 - t) / (t2 - t1) * P1[0] + (t - t1) / (t2 - t1) * P2[0];
            double A2_y = (t2 - t) / (t2 - t1) * P1[1] + (t - t1) / (t2 - t1) * P2[1];
            double A3_x = (t3 - t) / (t3 - t2) * P2[0] + (t - t2) / (t3 - t2) * P3[0];
            double A3_y = (t3 - t) / (t3 - t2) * P2[1] + (t - t2) / (t3 - t2) * P3[1];
            double B1_x = (t2 - t) / (t2 - t0) * A1_x + (t - t0) / (t2 - t0) * A2_x;
            double B1_y = (t2 - t) / (t2 - t0) * A1_y + (t - t0) / (t2 - t0) * A2_y;
            double B2_x = (t3 - t) / (t3 - t1) * A2_x + (t - t1) / (t3 - t1) * A3_x;
            double B2_y = (t3 - t) / (t3 - t1) * A2_y + (t - t1) / (t3 - t1) * A3_y;
            double C_x = (t2 - t) / (t2 - t1) * B1_x + (t - t1) / (t2 - t1) * B2_x;
            double C_y = (t2 - t) / (t2 - t1) * B1_y + (t - t1) / (t2 - t1) * B2_y;
            /* C_x = floor(C_x);
             C_y = floor(C_y);*/


             //c2 << C_x, C_y;
             //C2.push_back(c2);

            c1 << C_x, C_y;
            C.push_back(c1);

            t = t + linespace;
        }


        /*  ��ȥ��β�㣬�м��Ĳ�ֵ  */

        for (int i = 0; i < numSplines; i++) {
            // �ӵ�0���㿪ʼ��ÿ��ȡ4����   
            std::vector<VectorXd>::const_iterator first = inputPoints.begin() + i;
            std::vector<VectorXd>::const_iterator last = inputPoints.begin() + i + 4;
            std::vector<VectorXd> inputPoints_part(first, last);

            // ��ʱ�������Ϊ4�����µ��ú���������else��֧
         //   spline = catMullRomSpline(inputPoints_part);   // ע�͵�

           

            P0[0] = inputPoints_part[0][0];
            P0[1] = inputPoints_part[0][1];
            P1[0] = inputPoints_part[1][0];
            P1[1] = inputPoints_part[1][1];
            P2[0] = inputPoints_part[2][0];
            P2[1] = inputPoints_part[2][1];
            P3[0] = inputPoints_part[3][0];
            P3[1] = inputPoints_part[3][1];


             t0 = 0;
             t1 = tj(t0, P0, P1);
             t2 = tj(t1, P1, P2);
             t3 = tj(t2, P2, P3);

            //  cout << "t2-t1==="<<t2-t1 << endl;

             // �������Ϊ�����֮��ļ��
              linespace = (t2 - t1) / numSpace;
            // cout << linespace << endl;
              t = t1;


            VectorXd c2(2);
            vector<VectorXd> C2;
            while (t <= t2) {
                double A1_x = (t1 - t) / (t1 - t0) * P0[0] + (t - t0) / (t1 - t0) * P1[0];
                double A1_y = (t1 - t) / (t1 - t0) * P0[1] + (t - t0) / (t1 - t0) * P1[1];
                double A2_x = (t2 - t) / (t2 - t1) * P1[0] + (t - t1) / (t2 - t1) * P2[0];
                double A2_y = (t2 - t) / (t2 - t1) * P1[1] + (t - t1) / (t2 - t1) * P2[1];
                double A3_x = (t3 - t) / (t3 - t2) * P2[0] + (t - t2) / (t3 - t2) * P3[0];
                double A3_y = (t3 - t) / (t3 - t2) * P2[1] + (t - t2) / (t3 - t2) * P3[1];
                double B1_x = (t2 - t) / (t2 - t0) * A1_x + (t - t0) / (t2 - t0) * A2_x;
                double B1_y = (t2 - t) / (t2 - t0) * A1_y + (t - t0) / (t2 - t0) * A2_y;
                double B2_x = (t3 - t) / (t3 - t1) * A2_x + (t - t1) / (t3 - t1) * A3_x;
                double B2_y = (t3 - t) / (t3 - t1) * A2_y + (t - t1) / (t3 - t1) * A3_y;
                double C_x = (t2 - t) / (t2 - t1) * B1_x + (t - t1) / (t2 - t1) * B2_x;
                double C_y = (t2 - t) / (t2 - t1) * B1_y + (t - t1) / (t2 - t1) * B2_y;
                /* C_x = floor(C_x);
                 C_y = floor(C_y);*/


                c2 << C_x, C_y;
                C.push_back(c2);

                t = t + linespace;
            }



        }


        /*   ����������֮��Ĳ�ֵ   */


        P0[0] = inputPoints[last-3][0];
        P0[1] = inputPoints[last-3][1];
        P1[0] = inputPoints[last-2][0];
        P1[1] = inputPoints[last-2][1];
        P2[0] = inputPoints[last-1][0];
        P2[1] = inputPoints[last-1][1];
        P3[0] = inputPoints[last][0];
        P3[1] = inputPoints[last][1];

        P_temp_ = P3;
        P_temp = 2 * P3 - P2;


        P0 = P1;
        P1 = P2;
        P2 = P3;
        P3 = P_temp;

         t0 = 0;
         t1 = tj(t0, P0, P1);
         t2 = tj(t1, P1, P2);
         t3 = tj(t2, P2, P3);

        //  cout << "t2-t1==="<<t2-t1 << endl;

         // �������Ϊ�����֮��ļ��
         linespace = (t2 - t1) / numSpace;
        // cout << linespace << endl;
         t = t1;


        VectorXd c3(2);
       //  vector<VectorXd> C2;
        while (t <= t2) {
            double A1_x = (t1 - t) / (t1 - t0) * P0[0] + (t - t0) / (t1 - t0) * P1[0];
            double A1_y = (t1 - t) / (t1 - t0) * P0[1] + (t - t0) / (t1 - t0) * P1[1];
            double A2_x = (t2 - t) / (t2 - t1) * P1[0] + (t - t1) / (t2 - t1) * P2[0];
            double A2_y = (t2 - t) / (t2 - t1) * P1[1] + (t - t1) / (t2 - t1) * P2[1];
            double A3_x = (t3 - t) / (t3 - t2) * P2[0] + (t - t2) / (t3 - t2) * P3[0];
            double A3_y = (t3 - t) / (t3 - t2) * P2[1] + (t - t2) / (t3 - t2) * P3[1];
            double B1_x = (t2 - t) / (t2 - t0) * A1_x + (t - t0) / (t2 - t0) * A2_x;
            double B1_y = (t2 - t) / (t2 - t0) * A1_y + (t - t0) / (t2 - t0) * A2_y;
            double B2_x = (t3 - t) / (t3 - t1) * A2_x + (t - t1) / (t3 - t1) * A3_x;
            double B2_y = (t3 - t) / (t3 - t1) * A2_y + (t - t1) / (t3 - t1) * A3_y;
            double C_x = (t2 - t) / (t2 - t1) * B1_x + (t - t1) / (t2 - t1) * B2_x;
            double C_y = (t2 - t) / (t2 - t1) * B1_y + (t - t1) / (t2 - t1) * B2_y;
            /* C_x = floor(C_x);
             C_y = floor(C_y);*/


             //c2 << C_x, C_y;
             //C2.push_back(c2);

            c3 << C_x, C_y;
            C.push_back(c3);

            t = t + linespace;
        }




        // �������յĽ��

        return C;
    }
    // �����������4��ʱ
   
    else {
        cout << "The number of input points must > 4" << endl;
        vector<VectorXd> nopoint;
        nopoint.push_back(VectorXd(0, 0));
        return nopoint;
    }
    
    
//    return  inputPoints;




}


