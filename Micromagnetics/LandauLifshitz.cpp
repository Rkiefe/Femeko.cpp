#include "../src/BEMc.cpp"

void LandauLifshitz(
	Eigen::Ref<Eigen::MatrixXd> m,
	Eigen::Ref<Eigen::MatrixXd> p,
	Eigen::Ref<Eigen::MatrixXi> t,
	Eigen::Ref<Eigen::MatrixXi> surfaceT,
	Eigen::Ref<Eigen::MatrixXd> normal,
	double* AE,
	double* VE)
{

	// Constants
	double mu0 = pi*4e-7; 		  // Vaccum magnetic permeability
	double giro = 2.210173e5 /mu0;// Gyromagnetic ratio (rad T-1 s-1)

	double scl = 1e-9; 		// Scale of the model (nm)

	double Ms = 860e3;
	double Aexc = 13e-12; 	// Exchange (J/m)
	double Aan = 0.0; 		// Anisotropy (J/m3)
	
	double totalTime = 0.4*1e-9*giro; 	// giro * seconds | total time of simulation
	double dt = 0.0117; 				// giro seconds

	Eigen::Vector3d uan(0.0,0.0,0.0); 		// Easy axis direction
	Eigen::Vector3d Hext(0.0,50e3,0.0); 	// A/m | Applied field

	double maxTorque = 1e-14;   // Maximum <|dm/dt|> allowed
	int maxAtt = 15000; 		// Maximum number of time steps

	// Process the mesh input data
	int nv = p.cols();
	int nt = t.cols();
	int ne = surfaceT.cols();

	std::vector<double> Vn(nv,0.0);
	std::vector<double> Lagrange(nv,0.0);
	for(int k = 0; k<nt; k++)
	{
		for(int i = 0; i<4; i++)
		{
			Vn[t(i,k)] 		 += VE[k];
			Lagrange[t(i,k)] += VE[k]/4;
		}
	}

	// Make the FEM-BEM matrices
	Eigen::MatrixXd A = denseStiffnessMatrix(p,t,VE); 		// nv by nv
	Eigen::MatrixXd B = Bmatrix(p,surfaceT,AE); 			// nv by ne
	Eigen::MatrixXd C = Cmatrix(p,surfaceT,normal,AE); 	// ne by nv
	Eigen::MatrixXd D = Dmatrix(p,surfaceT,AE); 			// ne by ne

	// Extend the matrix
	// LHS = [-A B; C D]
	Eigen::MatrixXd LHS(nv+ne,nv+ne);
	LHS.block(0,0,nv,nv) 		= -A;
	LHS.block(0, nv, nv, ne) 	= B;
	LHS.block(nv, 0, ne, nv) 	= C;
	LHS.block(nv, nv, ne, ne) 	= D;


	// Calculate initial magnetic field
	Eigen::MatrixXd Heff = Eigen::MatrixXd::Zero(3,nv);
	for(int i = 0; i<nv; i++)
	{
		Heff.col(i) += Hext;
	}

	// Demagnetizing field
	Eigen::MatrixXd Hd = Eigen::MatrixXd::Zero(3,nv);
	BEMdmag(Hd, p, t, surfaceT, normal, AE, VE, Vn.data(), LHS, m);
	Hd *= mu0*Ms;

	// Exchange field
	Eigen::MatrixXd Hexc(3,nv);
	for(int i = 0; i<3; i++)
	{
		Hexc.row(i) = -2*Aexc * (A*(m.row(i).transpose())).transpose();
	}

	// Anisotropy field
	Eigen::MatrixXd Han(3,nv);
	for(int i = 0; i<nv; i++)
	{
		Han.col(i) = 2*Aan/Ms *uan.dot(m.col(i)) * uan;
	}

	// Effective field
	Heff += Hd + Hexc + Han;

	// double time = 0.0;
	// int att = 0;
	// while (time < totalTime && att < maxAtt)
	// {
	// 	att += 1;
	// 	time += dt;


	// }

}

