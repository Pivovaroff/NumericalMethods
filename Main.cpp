#include "Funct.h"
//#include <iomanip>


datatype a, b;
funct f, g, phi, psi;
funct2 f22, g2, phi2, psi2,  AnalSol;
funct3 T2;
//functsys fs;
dsys df;
bigexample fantasy;
tVector X0 = nullptr;

int main()
{
	setlocale(LC_ALL, "rus");
	size_t key = 10;
	//		tMatrix A = nullptr;
	////		tMatrix LU = nullptr;
	//		int dim;
	//		dim = Dimension();
	//		CreatingArray(A, dim);
	//		EnteringArray(A, dim);
	////		A[0][0] = 1;
	////		A[0][1] = 2;
	////		A[0][2] = 1;
	////		A[1][0] = 2;
	////		A[1][1] = 1;
	////		A[1][2] = 1;
	////		A[2][0] = 1;
	////		A[2][1] = -1;
	////		A[2][2] = 2;
	////		LU = LUDecomposition(A, dim);
	//		LUDecomposition_pr(A, dim);
	//		PrintingArray(A, dim);
	//
	////		tVector f = nullptr;
	////		CreatingVector(f, dim);
	////		EnteringVector(f, dim);
	////	//	f[0] = 1;
	////	//	f[1] = 0;
	////	//	f[2] = 0;
	//////		LUSolve(LU, f, dim);
	////
	////		DeletingArray(A, dim);
	////		if (LU != nullptr) { DeletingArray(LU, dim); }


		//std::cout<<std::setprecision(20)<<1.234e-15<<endl;
	while (key != 0)
	{
		clear();
		key = PrintMENU();
		switch (key)
		{
		case 1:
		{
			AnalSol = HeatAnal;
			//HeatTransfer_Example1(0.5);
			//HeatTransfer_Example2(0.5);
			//HeatTransfer_Example1_Anal(0, AnalSol);
			//HeatTransfer_Example1_Quasilinear(1);
			//HeatTransfer_Example1_Quasilinear_Iterations(1);
			break;
		}
		case 3:
		{
			f = f_w2;
			g = g_w;
			phi = phi_w;
			psi = psi_w;
			Wave(f, g, phi, psi);
			break;
		}
		case 2:
		{
			CreatingVector(X0, 2);
			X0[0] = 1;
			X0[1] = 0.1;
			df = dsys2;
			//Pendulum_Euler(df, X0, 2);
			//Pendulum_SymPlan(df, X0, 2);
			//ImplicitEuler(df, X0, 2);
			//Runge_Kutt(df, X0, 2);
			fantasy = BigExample;
			SolveBigExampleRK(fantasy, X0);
			DeletingVector(X0);
			break;
			//break;
		}
		case 4:
		{
			f = f_w;
			g = g_w;
			phi = phi_w;
			psi = psi_w;
			Waveanal(f, g, phi, psi);
			break;
		}
		case 5:
		{
			f22 = f_fe;
			g2 = g_fe;
			phi2 = phi_fe;
			psi2 = psi_fe;
			T2 = T_fe;
			datatype d, ksi;
			cout << "D=" ;
			cin >> d;
			cout << endl;
			cout << "ksi=";
			cin >> ksi;
			cout << endl;
			
			//Finite_Element_BYTHEDICK( ksi, d);
			break;
		}
		case 6:
		{
			CreatingVector(X0, 2);
			X0[0] = 1;
			X0[1] = 1;
			//df = dsys2;
			df = Runge1;
			Rozenbrok_Runge(df, X0, 2);
			//Rozenbrok(df, X0, 2);
			ClearRAM();
			break;
		}
		case 0:
		{
			cin.get();
			return 0;
			break;
		}
		}
	}
	//////////////////////////////////////////////////
//		switch (key)
//		{
//		case 1:
//		{
//			clear();
//			CreatingVector(X0, 2);
//			X0[0] = 0;
//			X0[1] = 0;
//			df = dsys1;
//			Runge_Kutt(df, X0, 2);
//			Predictor_Corrector(df, X0, 2);
//			Enclosed_Runge_Kutt(df, X0, 2);
//			Adams(df, X0, 2);
//			ImplicitEuler(df, X0, 2);
//			DeletingVector(X0);
//			break;
//		}
//		case 2:
//		{
//			clear();
//			CreatingVector(X0, 2);
//			X0[0] = 0;
//			X0[1] = 0;
//			df = dsys2;
//			Runge_Kutt(df, X0, 2);
//			Predictor_Corrector(df, X0, 2);
//			Enclosed_Runge_Kutt(df, X0, 2);
//			Adams(df, X0, 2);
//			DeletingVector(X0);
//			break;
//		}
//		case 3:
//		{
//			clear();
//			CreatingVector(X0, 3);
//			X0[0] = 2;
//			X0[1] = 0;
//			X0[2] = 0;
//			df = dsys3;
//			Runge_Kutt(df, X0, 3);
//			Predictor_Corrector(df, X0, 3);
//			Enclosed_Runge_Kutt(df, X0, 3);
//			Adams(df, X0, 3);
//			DeletingVector(X0);
//			break;
//		}
//		case 4:
//		{
//			clear();
//			CreatingVector(X0, 2);
//			X0[0] = 1;
//			X0[1] = 0;
//			df = Oscill;
//			Pendulum_Euler(df, X0, 2);
//			Pendulum_SymPlan(df, X0, 2);
//			ImplicitEuler(df, X0, 2);
//			Runge_Kutt2(df, X0, 2);
//			DeletingVector(X0);
//			break;
//		}
//		case 5:
//		{
//			clear();
//			CreatingVector(X0, 2);
//			X0[0] = 1;
//			X0[1] = 0.1;
//			fantasy = BigExample;
//			SolveBigExampleRK(fantasy, X0);
//			DeletingVector(X0);
//			break;
//		}
//		case 6:
//		{
//			clear();
//			//cout << "Runge-Kutt vs Analytical Solve" << endl;
//			CreatingVector(X0, 2);
//			X0[0] = 6;
//			X0[1] = 5;
//			df = dsysAnal;
//			RKvsAS(df, X0, 2);
//
//			break;
//		}
//		case 0:
//		{
//			//DeletingVector(Standings);
//			ClearRAM();
//			cin.get();
//			return 0;
//			break;
//		}
//		}
//	}


	cin.get();
	return 0;
}
