// Wrapper to C++ micromagnetics simulation to be called in julia

#include "LandauLifshitz.cpp" // Include Femeko C++ FEM functions

extern "C"{
	
	// Get the demagnetizing field on the mesh nodes
	void LL(double* m_in,
			double* p_in,
		    int* t_in,
		    int* surfaceT_in,
		    double* normal_in,
		    double* AE, double* VE,
		    int nv, int nt, int ne,
		    double* M_avg_out,
		    int maxAtt)
	{

		// Node coordinates
		Eigen::Map<Eigen::MatrixXd> p(p_in,3,nv);
		
		// Element connectivity
		Eigen::Map<Eigen::MatrixXi> t(t_in,4,nt);

		// Surface element node connectivity
		Eigen::Map<Eigen::MatrixXi> surfaceT(surfaceT_in,4,ne);
		
		// Surface normals
		Eigen::Map<Eigen::MatrixXd> normal(normal_in,3,ne);

		// Magnetization
		Eigen::Map<Eigen::MatrixXd> m(m_in,3,nv);

		// Run
		Eigen::MatrixXd M_avg = LandauLifshitz(
										  	m, p, t,
										  	surfaceT,
										  	normal,
										  	AE, VE,
										  	maxAtt);

		// Copy results to Julia output
        Eigen::Map<Eigen::MatrixXd>(M_avg_out, 3, maxAtt) = M_avg;

	} // Wrapper to C++ demag field

} // extern C