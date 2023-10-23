#ifndef _TRAJECTORY_GENERATOR_WAYPOINT_H_
#define _TRAJECTORY_GENERATOR_WAYPOINT_H_

#include <Eigen/Eigen>
#include <vector>

class TrajectoryGeneratorWaypoint {
    private:
		double _qp_cost;
		Eigen::MatrixXd _Q;
		Eigen::VectorXd _Px, _Py, _Pz;
    public:
        TrajectoryGeneratorWaypoint();

        ~TrajectoryGeneratorWaypoint();

        Eigen::MatrixXd PolyQPGeneration(
            const int order,
            const Eigen::MatrixXd &Path,
            const Eigen::MatrixXd &Vel,
            const Eigen::MatrixXd &Acc,
            const Eigen::VectorXd &Time);
        
        int Factorial(int x);

        Eigen::MatrixXd getM(const int n_seg,const int d_order,const int p_numld,const Eigen::VectorXd& time);
        Eigen::MatrixXd getQ(const int n_seg,const int d_order,const int p_numld,const Eigen::VectorXd& time);
        Eigen::MatrixXd getCt(const int n_seg,const int d_order);
};
        

#endif
