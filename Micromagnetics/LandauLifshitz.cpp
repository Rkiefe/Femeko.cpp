#include "../src/BEMc.cpp"


// Find new magnetization
void timeStep(
	Eigen::Ref<Eigen::Vector3d> m,
	Eigen::Ref<Eigen::Vector3d> H,
	Eigen::Ref<Eigen::Vector3d> Hold,
	Eigen::Ref<Eigen::Vector3d> Heff,
	double dt,
	double giro,
	double damp,
	double precession)
{
	double d = dt*giro/2.0;
	
	// m (n+1)
	Eigen::Vector3d mNew, m12, H12;

	// 1) Initial guess of the new magnetic field
	H12 = 1.5 *H - 0.5 *Hold;

	// Find new m(n+1) until it converges
	Eigen::Vector3d oldGuess = m;
	Eigen::PartialPivLU<Eigen::Matrix3d> luSolver;

	double err = 1.0;
	int att = 0;
	Eigen::Matrix3d mat;
	while(err > 1e-6 && att < 1000)
	{
		att++;

		// 2) m (n+1) from m(n) and H(n+1/2)
		mat << 1.0,      d*H12(2),  -d*H12(1),
	          	-d*H12(2),  1.0,     d*H12(0),
	           	d*H12(1), -d*H12(0),  1.0;

	    luSolver.compute(mat);

	    mNew = luSolver.solve(m - d * m.cross(H12));

	    // mNew = mat.lu().solve(m - d * m.cross(H12));  // Best for 3x3 systems

	    // 3) m (n+1/2)
	    m12 = 0.5*(m + mNew);

	    // 4) Get H(n+1/2) from m (n+1/2)
	    H12 = precession*Heff + damp*m12.cross(Heff);

	    // Max difference between mNew and oldGuess
	    for(int i = 0; i<3; i++)
	    {
	    	err = std::max(err,std::abs(mNew(i) - oldGuess(i)));
	    }

	    // Update oldGuess
	    oldGuess = mNew;
	}

	// Update input m to the new m
	m = mNew;
}

Eigen::MatrixXd LandauLifshitz(
	Eigen::Ref<Eigen::MatrixXd> m,
	Eigen::Ref<Eigen::MatrixXd> p,
	Eigen::Ref<Eigen::MatrixXi> t,
	Eigen::Ref<Eigen::MatrixXi> surfaceT,
	Eigen::Ref<Eigen::MatrixXd> normal,
	double* AE,
	double* VE,
	int maxAtt)
{

	// Constants
	double mu0 = pi*4e-7; 		  // Vaccum magnetic permeability
	double giro = 2.210173e5 /mu0;// Gyromagnetic ratio (rad T-1 s-1)

	double scl = 1e-9; 		// Scale of the model (nm)

	double Ms = 860e3;
	double Aexc = 13e-12; 	// Exchange (J/m)
	double Aan = 0.0; 		// Anisotropy (J/m3)

	double damp = 0.1; 		// Damping [0,1]
	double precession = 1.0;// 0.0 or 1.0 | Consider or not the precession term
	
	double totalTime = 0.4*1e-9*giro; 	// giro * seconds | total time of simulation
	double dt = 0.0117; 				// giro seconds

	Eigen::Vector3d uan(0.0,0.0,0.0); 		// Easy axis direction
	Eigen::Vector3d Hext(0.0,50e3,0.0); 	// A/m | Applied field

	double maxTorque = 1e-14;   // Maximum <|dm/dt|> allowed
	// maxAtt = 10; 			// Maximum number of time steps

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

	// Factorize LHS once
	Eigen::PartialPivLU<Eigen::MatrixXd> lu(LHS);


	// Calculate initial magnetic field
	Eigen::MatrixXd Heff = Eigen::MatrixXd::Zero(3,nv);
	for(int i = 0; i<nv; i++)
	{
		Heff.col(i) += mu0 *Hext;
	}

	// Demagnetizing field
	Eigen::MatrixXd Hd = Eigen::MatrixXd::Zero(3,nv);
	BEMdmag(Hd, p, t, surfaceT, normal, AE, VE, Vn.data(), lu, m);
	Hd *= mu0*Ms;

	// Exchange field
	Eigen::MatrixXd Hexc(3,nv);
	for(int i = 0; i<3; i++)
	{
		// Update the entire xyz row
		Hexc.row(i) = -2*Aexc * (A*(m.row(i).transpose())).transpose();
		
		// Update each node to be at the correct units
		for(int nd = 0; nd<nv; nd++)
		{
			Hexc(i,nd) /= Ms*pow(scl,2)*Lagrange[nd];
		} 
	}

	// Anisotropy field
	Eigen::MatrixXd Han(3,nv);
	for(int i = 0; i<nv; i++)
	{
		Han.col(i) = 2*Aan/Ms *uan.dot(m.col(i)) * uan;
	}

	// Effective field
	Heff += Hd + Hexc + Han;

	// H = Heff + damp M cross Heff
	Eigen::MatrixXd H(3,nv);
	for(int i = 0; i<nv; i++)
	{
		// Eigen::Vector3d aux = m.col(i).head<3>().cross(Heff.col(i).head<3>());
		H.col(i) = precession * Heff.col(i) + damp * m.col(i).head<3>().cross(Heff.col(i).head<3>());
	}

	// Time iteration
	Eigen::MatrixXd Hold = H;
	
	Eigen::MatrixXd M_avg = Eigen::MatrixXd::Zero(3,maxAtt);

	double time = 0.0;
	int att = 0;
	while (time < totalTime && att < maxAtt)
	{
		// New magnetization
		#pragma omp parallel for
		for(int i = 0; i<nv; i++)
		{
			timeStep(m.col(i).head<3>(),
					 H.col(i),
					 Hold.col(i),
					 Heff.col(i),
					 dt,
					 1.0,
					 damp,
					 precession); // Considering time was normalized by giro
		

		} // New magnetization

		// -- New magnetic field --
		Hold = H; // Store the old one

		// Reset effective field
		for(int i = 0; i<nv; i++)
		{
			Heff.col(i) = mu0 *Hext;
		}

		// Demagnetizing field
		Hd.setZero();
		BEMdmag(Hd, p, t, surfaceT, normal, AE, VE, Vn.data(), lu, m);
		Hd *= mu0*Ms;

		// Exchange field
		for(int i = 0; i<3; i++)
		{
			// Update the entire xyz row
			Hexc.row(i) = -2*Aexc * (A*(m.row(i).transpose())).transpose();
			
			// Update each node to be at the correct units
			for(int nd = 0; nd<nv; nd++)
			{
				Hexc(i,nd) /= Ms*pow(scl,2)*Lagrange[nd];
			}
		}

		// Anisotropy field
		for(int i = 0; i<nv; i++)
		{
			Han.col(i) = 2*Aan/Ms *uan.dot(m.col(i)) * uan;
		}

		// Add all field contributions
		Heff += Hd + Hexc + Han;

		// New H = Heff + damp M cross Heff
		for(int i = 0; i<nv; i++)
		{
			// Eigen::Vector3d aux = m.col(i).head<3>().cross(Heff.col(i).head<3>());
			H.col(i) = precession * Heff.col(i) + damp * m.col(i).head<3>().cross(Heff.col(i).head<3>());
		}

		// Average magnetization
		M_avg.col(att) = m.rowwise().mean();

		att++;
		time += dt;

		// Log results
		if(!(att%100))
		{
			std::cout << att << std::endl;
		}

	} // Time iteration

	return M_avg;

} // Landau Lifshitz solver

