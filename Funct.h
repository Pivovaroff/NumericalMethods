#pragma once
#define _USE_MATH_DEFINES 
#include <iostream>
#include <locale>
#include <cmath>
#include <fstream>
#include <ctime>
#include <cstdio>
#include <stdio.h>
#include <math.h>
using namespace std;

typedef double datatype;
typedef datatype* tVector;
typedef datatype** tMatrix;

typedef datatype(*funct)(datatype);
typedef datatype(*funct2)(datatype, datatype);
typedef datatype(*funct3)(datatype, datatype, datatype);
typedef tVector(*functsys) (datatype, datatype);
typedef tVector(*dsys)(tVector&, int);
typedef tVector(*bigexample)(tVector&, datatype);
//typedef datatype(*HeatStream)(datatype);

//datatype eps = 1.e-10;
//int N = 0;

void CreatingArray(tMatrix &mass, size_t n);//выделение памяти матрицы

void EnteringArray(tMatrix &mass, size_t n); //ввод значений матрицы

void PrintingArray(tMatrix &mass, size_t n);//вывод матрицы

void DeletingArray(tMatrix &Matrix, size_t N);//высвобождение памяти матрицы

void DeletingArray(tMatrix &Matrix);//высвобождение памяти матрицы


tMatrix Summ(tMatrix &mass1, tMatrix &mass2, size_t n);//матричная сумма

tMatrix Difference(tMatrix &mass1, tMatrix &mass2, size_t n);//матричная разность

tMatrix Multiplication(tMatrix &mass1, tMatrix &mass2, size_t n);//матричное перемножение

tMatrix MatrCopy(tMatrix &mass, size_t n);

tMatrix Transpose(tMatrix &A, size_t n);

void DeletingVector(tVector &vect);//высвобождение памяти вектора

void CreatingVector(tVector &vect, size_t n);//выделение памяти для вектора

void EnteringVector(tVector &vect, size_t n);//ввод вектора

void PrintingVector(tVector &vect);//вывод вектора

void PrintingVector(tVector &vect, int n);

tVector Summ(tVector &vect1, tVector &vect2, int n);

tVector Difference(tVector &vect1, tVector &vect2, size_t n);

tVector NumberVectorMultiplication(datatype num, tVector &mass, size_t n);

datatype Scalar(tVector &vect1, tVector &vect2);

tVector VectCopy(tVector &vect, size_t n);

tVector Gauss(tMatrix &A, tVector &b);

tVector Gauss(tMatrix &A, tVector &b, int N);

datatype Determinant(tMatrix &A1, tVector &B, size_t n); //приведение матрицы к верхнегреулольному виду

datatype NormOfVector(tVector &vect1, size_t n);

datatype NormOfVectorC(tVector &vect1, size_t n);

tVector MVMultiplication(tMatrix &mass1, tVector &vect, int n);//перемножение вектора на матрицу (слева)

void ResultsSave(tMatrix &A, tVector &B, tVector &X, size_t n, datatype iter); //сохранение результата в файл

void ResultsSavetest(tMatrix &A, size_t n);

tMatrix Inverse(tMatrix &A, size_t n);

tVector QRmethod(tMatrix &mass1, tVector &vect);//QR-разложение

datatype ConditionNumberOne(tMatrix &A, size_t n);  //max via columns

datatype ConditionNumberInf(tMatrix &A, size_t n); //max via rows

datatype NormOfDifferenceB(tMatrix &A, tVector &B, tVector &X, size_t n);

void StabilityOfSolve(tMatrix &A, tVector &B, tVector &X1, size_t n);

void EstimateOfCond(tMatrix &A, tVector &B, tVector &X, size_t n);


void EnteringFromFile(tMatrix &A, tVector &B);//вывод из файла матрицы А и вектора B.

datatype NormOfMatrix(tMatrix &A);

tMatrix NumberMatrixMultiplication(datatype num, tMatrix &A);

tMatrix NumberMatrixMultiplication(datatype num, tMatrix &A, int n);

void NMMultiplication(datatype num, tMatrix &A, int n);

tVector NumberVectorMultiplication(datatype num, tVector &mass);

tVector SimpleIterationMethod(tMatrix &a, tVector &b);

tVector JacobiMethod(tMatrix &a, tVector &b);

tVector JacobiMethod2(tMatrix &a, tVector &b);

tVector ZeydelMethod(tMatrix &a, tVector &b);

tVector Relax(tMatrix &a, tVector &b);

void CreatingDataFile();

tVector Relax3Diag(tMatrix &a, tVector &b, datatype omega);

tVector Relax3Diag(tMatrix &a, tVector &b, datatype omega, size_t N);

tVector Shuttle(tVector &l, tVector &d, tVector &u, tVector &b, int m);

void EnteringFromFile(tMatrix &A);

tMatrix Hessenberg(tMatrix &A);

void PrintingArray(tMatrix &A);

tMatrix QR_T(tMatrix &A);

tMatrix Transpose(tMatrix &A);

tMatrix Multiplication(tMatrix &mass1, tMatrix &mass2);

bool TrueQ(tMatrix &A);

tMatrix One(tMatrix &A);

tMatrix EugineVectors(tMatrix &mass1, tVector &lambda);

tMatrix Reley(tMatrix &mass1);

void CreatingVector(tVector &vect);

tVector MainDiag(tMatrix &A, tVector &diag1);

tMatrix MatrCopy(tMatrix &mass);

datatype f1(datatype x);

datatype f2(datatype x);

datatype f3(datatype x);

datatype f4(datatype x);

datatype f5(datatype x);

datatype f_w(datatype x);

datatype f_w2(datatype x);

datatype g_w(datatype x);

datatype phi_w(datatype x);

datatype psi_w(datatype x);

datatype Kx(datatype x);

datatype Kxpr(datatype x);

datatype P_Zero(datatype x);

datatype P_Heat(datatype x);

datatype K10(datatype x);

datatype P1(datatype x);

datatype P2(datatype x);

datatype P3(datatype x);

datatype P4(datatype x);

size_t PrintMENU();

void clear();

tVector FindRoot_WilO4KA(funct &f, datatype a, datatype b);

tVector FindRoot_Neuton(funct &f, datatype a, datatype b);

tVector FindRoot_Neuton_Anal(funct &f, datatype a, datatype b);

tVector FindRoot_Neuton_Num(funct &f, datatype a, datatype b);

tVector dsys1(tVector &X0, int n);

tVector dsys2(tVector &X0, int n);

tVector dsys3(tVector &X0, int n);

tVector dsys4(tVector &X0, int n);

tVector dsysAnal(tVector &X0, int n);

tVector dsysAnalSolve(datatype t);

tVector dsysRozenbrok1(tVector &X0, int n);

tVector Oscill(tVector &X0, int n);

tVector fsys2(datatype x, datatype y);

tVector fsys3(datatype x, datatype y);

tVector FindRoot_WilO4KA_CQ5(funct &f);

tMatrix FindRoot_Sys_Num(functsys &fsys, datatype a, datatype b);

tMatrix FindRoot_Sys2(functsys &fsys, datatype a, datatype b);

tMatrix FindRoot_Sys1(functsys &fsys, datatype a, datatype b);

void Runge_Kutt(dsys &fsys, tVector &X0, int n);

void Predictor_Corrector(dsys &fsys, tVector &X, int n);

tVector Summ(tVector &X1, tVector &X2, tVector &X3, tVector &X4, int n);

tVector Summ(tVector &X1, tVector &X2, tVector &X3, tVector &X4, tVector&X5, int n);

void Pendulum_Euler(dsys &fsys, tVector &X, int n);

tVector FindRoot_SymPlan(dsys &fsys, tVector &X, int n);

void Pendulum_SymPlan(dsys &fsys, tVector &X, int n);

void Enclosed_Runge_Kutt(dsys &fsys, tVector &X, int n);

void Adams(dsys &fsys, tVector &X, int n);

tVector FindRoot_ImplicitEuler(dsys &fsys, tVector &X, int n);

void ImplicitEuler(dsys &fsys, tVector &X, int n);

tVector BigExample(tVector& X0, datatype t);

void SolveBigExampleRK(bigexample &fsys, tVector &X);

void RKvsAS(dsys &fsys, tVector &X, int n);

//tMatrix LUDecomposition(tMatrix &A, int n);

tMatrix LUDecomposition(tMatrix &A, int n);

void LUDecomposition_pr(tMatrix &A, int n);

tVector LUSolve(tMatrix &LU, tVector &f, int n);

int Dimension();

void ClearRAM();

void ADvsAS(dsys &fsys, tVector &X, int n);

void PKvsAS(dsys &fsys, tVector &X, int n);

void Runge_Kutt2(dsys &fsys, tVector &X, int n);

datatype Integrate(funct &f, datatype a, datatype b);

void HeatTransfer(datatype sigma);

void HeatTransfer_Example1(datatype sigma);

void HeatTransfer_Example1_Quasilinear(datatype sigma);

void HeatTransfer_Example1_Quasilinear_Iterations(datatype sigma);

void HeatTransfer_Example2(datatype sigma);

datatype HeatAnal(datatype x, datatype t);

datatype P100(datatype x);

void HeatTransfer_Example2_Anal(datatype sigma, funct2 &AnalSolve);

void HeatTransfer_Example1_Anal(datatype sigma, funct2 &AnalSolve);

void HeatTransfer_Example2_Quasilinear(datatype sigma);

datatype K_Example3(datatype u);

datatype AnalSolve_Example3(datatype x, datatype t);

void HeatTransfer_Example3(datatype sigma, funct2& AnalSolve_Example3);

void HeatTransfer_Example3_Anal(datatype sigma, funct2& AnalSolve_Example3);

//void HeatTransfer_Example3_Iter(datatype sigma, funct2& AnalSolve);

void HeatTransfer_Example3_Iter_Anal(datatype sigma, funct2& AnalSolve);

void HeatTransfer_Example3_Iter(datatype sigma, funct2& AnalSolve);

void Wave(funct &f, funct &g, funct &phi, funct &psi);

void Waveanal(funct &f, funct &g, funct &phi, funct &psi);

void NVMultiplication(datatype num, tVector &v, int n);

void Rozenbrok(dsys &fsys, tVector &X, int n);

void Rozenbrok_Runge(dsys &fsys, tVector &X, int n);

tVector Runge1(tVector &X0, int n);

datatype f_fe(datatype x, datatype t);

datatype g_fe(datatype x, datatype t);

datatype phi_fe(datatype x, datatype t);

datatype psi_fe(datatype x, datatype t);

datatype T_fe(datatype x, datatype t);

datatype T_fe(datatype x, datatype t, datatype betta);

void Finite_Element(funct2 &f, funct2 &g, funct2 &phi, funct2 &psi, funct3 &T,datatype ksi, datatype d, datatype betta);

void Lab3_Ex1();
