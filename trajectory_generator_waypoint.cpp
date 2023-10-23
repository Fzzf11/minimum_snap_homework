#include "trajectory_generator_waypoint.h"
#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;    
using namespace Eigen;

TrajectoryGeneratorWaypoint::TrajectoryGeneratorWaypoint(){}
TrajectoryGeneratorWaypoint::~TrajectoryGeneratorWaypoint(){}

//define factorial function, input i, output i!
int TrajectoryGeneratorWaypoint::Factorial(int x)
{
    int fac = 1;
    for(int i = x; i > 0; i--)
        fac = fac * i;
    return fac;
}
/*

    STEP 2: Learn the "Closed-form solution to minimum snap" in L5, then finish this PolyQPGeneration function

    variable declaration: input       const int d_order,                    // the order of derivative
                                      const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
                                      const Eigen::MatrixXd &Vel,           // boundary velocity
                                      const Eigen::MatrixXd &Acc,           // boundary acceleration
                                      const Eigen::VectorXd &Time)          // time allocation in each segment
                          output      MatrixXd PolyCoeff(m, 3 * p_num1d);   // position(x,y,z), so we need (3 * p_num1d) coefficients

*/

Eigen::MatrixXd TrajectoryGeneratorWaypoint::getM(const int n_seg,const int d_order,const int p_numld,const Eigen::VectorXd& time)
{
    // get M
    /*
        n_seg: the number of segments
        d_order: the deferential order----->4
        p_numld: the number of variables in each segment
        time: size: m*1
        M: size: p_numld*m x p_numld*m
    */
    int M_row = p_numld*n_seg;
    int M_col = p_numld*n_seg;
    double t;

    MatrixXd M = MatrixXd::Zero(M_row,M_col);

    MatrixXd cof(d_order,p_numld);
    if(d_order == 4){
        cof <<  1,1,1,1,1,1,1,1,
                0,1,2,3,4,5,6,7,
                0,0,2,6,12,20,30,42,
                0,0,0,6,24,60,120,210;
    }else if(d_order == 3){
        cof <<  1,1,1,1,1,1,
                0,1,2,3,4,5,
                0,0,2,6,12,20;
    }
    
    for(int k=0;k<n_seg;++k){
        MatrixXd M_k = MatrixXd::Zero(p_numld,p_numld);
        t = time(k);
        for(int i=0;i<p_numld;i++){
            if(i<=3){
                M_k(i,i) = cof(i,i);
            }
            else{
                M_k(i,i-4) = cof(i-4,i-4);
                for(int j=i-3;j<p_numld;j++){
                    M_k(i,j) = cof(i-4,j)*pow(t,j-(i-4));
                }
            }
        }
        M.block(k*p_numld,k*p_numld,p_numld,p_numld) = M_k;
    }

    cout << M << endl;
    return M;
}

Eigen::MatrixXd TrajectoryGeneratorWaypoint::getQ(const int n_seg,const int d_order,const int p_numld,const Eigen::VectorXd& time)
{
    // get Q
    /*
        n_seg: the number of segments
        d_order: the deferential order----->4
        p_numld: the number of variables in each segment
        time: size: m*1
    return:
        Q:n_seg*p_numld x n_seg*p_numld
    */
    int Q_row = n_seg*p_numld;
    int Q_col = n_seg*p_numld;
    MatrixXd Q = MatrixXd::Zero(Q_row,Q_col);
    double t;

    for(int k=0;k<n_seg;k++){
        t = time(k);
        MatrixXd Q_k = MatrixXd::Zero(p_numld,p_numld);
        for(int i=0;i<p_numld;i++){
            for(int j=0;j<p_numld;j++){
                if(i<d_order || j<d_order){
                    continue;
                }
                Q_k(i,j)=(Factorial(i)*Factorial(j)*pow(t,i+j-7))/(Factorial(i-4)*Factorial(j-4)*(i+j-7));
            }
        }
        Q.block(k*p_numld,k*p_numld,p_numld,p_numld) = Q_k;
    }

    cout << Q << endl;

    return Q;
}

Eigen::MatrixXd TrajectoryGeneratorWaypoint::getCt(const int n_seg,const int d_order)
{
    // get Ct
    /*
        n_seg:the number of segments
        d_order:the order of derivative
    return:
        Ct:n_seg*2*d_order x (n_seg*2*d_order-(n_seg-1)*d_order)
    */
    int Ct_row = n_seg*2*d_order;
    int Ct_col = n_seg*2*d_order-(n_seg-1)*d_order;
    MatrixXd Ct = MatrixXd::Zero(Ct_row,Ct_col);
    //dF
    for(int k=0;k<d_order;k++){
        Ct(k,k) = 1;
    }

    for(int k=1;k<n_seg;k++){
        Ct(d_order+(k-1)*2*d_order,d_order+(k-1)*1) = 1;
        Ct(2*d_order+(k-1)*2*d_order,d_order+(k-1)*1) = 1;
    }

    for(int k=0;k<d_order;k++){
        Ct(d_order+2*d_order*(n_seg-1)+k,d_order+(n_seg-1)*1+k) = 1;
    }

    for(int k=1;k<n_seg;k++){
        for(int i=1;i<d_order;i++){
            Ct((k-1)*2*d_order+d_order+i,2*d_order+n_seg-1+(k-1)*(d_order-1)+i-1) = 1;
            Ct((k-1)*2*d_order+d_order+i+d_order,2*d_order+n_seg-1+(k-1)*(d_order-1)+i-1) = 1;
        }
    }

    cout << Ct << endl;

    return Ct;
}

Eigen::MatrixXd TrajectoryGeneratorWaypoint::PolyQPGeneration(
            const int d_order,                    // the order of derivative ----->4
            const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
            const Eigen::MatrixXd &Vel,           // boundary velocity
            const Eigen::MatrixXd &Acc,           // boundary acceleration
            const Eigen::VectorXd &Time)          // time allocation in each segment
{
    // enforce initial and final velocity and accleration, for higher order derivatives, just assume them be 0;
    int p_order   = 2 * d_order - 1;              // the order of polynomial
    int p_num1d   = p_order + 1;                  // the number of variables in each segment

    int m = Time.size();                          // the number of segments
    MatrixXd PolyCoeff = MatrixXd::Zero(m, 3 * p_num1d);           // position(x,y,z), so we need (3 * p_num1d) coefficients
    VectorXd Px(p_num1d * m), Py(p_num1d * m), Pz(p_num1d * m);

    /*   Produce Mapping Matrix A to the entire trajectory, A is a mapping matrix that maps polynomial coefficients to derivatives.   */
    //get Q
    MatrixXd Q = getQ(m,d_order,p_num1d,Time);

    //get M
    MatrixXd M = getM(m,d_order,p_num1d,Time);

    //get Ct
    MatrixXd Ct = getCt(m,d_order);

    MatrixXd M_inv = M.inverse();
    MatrixXd M_inv_T = M_inv.transpose();

    MatrixXd Ct_T = Ct.transpose(); 

    MatrixXd R = Ct_T*M_inv_T*Q*M_inv*Ct;

    int num_dF = 2*d_order+m-1;
    int num_dP = (d_order-1)*(m-1);

    MatrixXd Rpp = R.bottomRightCorner(num_dP,num_dP);
    MatrixXd Rfp = R.topRightCorner(num_dF,num_dP);



    /*   Produce the dereivatives in X, Y and Z axis directly.  */
    for(int dim=0;dim<3;dim++){
        VectorXd waypoint = Path.col(dim);
        VectorXd dF = VectorXd::Zero(num_dF);

        dF(0) = waypoint(0);

        for(int i=0;i<m-1;++i){
            dF(d_order+i) = waypoint(i+1);
        }

        dF(d_order+(m-1)) = waypoint(m);

        VectorXd dP = (-1)*Rpp.inverse()*Rfp.transpose()*dF;
        VectorXd dtotal(num_dF+num_dP);
        dtotal << dF,dP;

        VectorXd poly_coef_ld = M_inv*Ct*dtotal;

        MatrixXd poly_coef_ld_T = poly_coef_ld.transpose();
        for(int k=0;k<m;k++){
                PolyCoeff.block(k,dim*p_num1d,1,p_num1d) = poly_coef_ld_T.block(0,k*p_num1d,1,p_num1d);
        }
    }    
    cout << PolyCoeff << endl;
    /*   Produce the Minimum Snap cost function, the Hessian Matrix   */
    return PolyCoeff;
}
