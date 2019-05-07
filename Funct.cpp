#include <iostream>
#include <locale>
#include <cmath>
#include <fstream>
#include <ctime>
#include <cstdio>
#include <stdio.h>
#include "Funct.h"
#include <iomanip>
using namespace std;

typedef double datatype;
typedef datatype* tVector;
typedef datatype** tMatrix;

datatype eps = 1.e-3;
datatype eps7 = 1.e-7;
datatype eps5 = 1.e-5;
tVector Standings = nullptr;
const datatype tao = 0.0100;
const datatype q = 1.e-3;
int N = 0;

const datatype h = 0.01;


/*  data type      byte          max value
bool               =  1         255.00
char               =  1         255.00
short int          =  2         32767.00
unsigned short int =  2         65535.00
int                =  4         2147483647.00
unsigned int       =  4         4294967295.00
long int           =  4         2147483647.00
unsigned long int  =  4         4294967295.00
float              =  4         2147483647.00
long float         =  8         9223372036854775800.00
double             =  8         9223372036854775800.00  */




void CreatingArray(tMatrix &mass, size_t n)//выделение памяти матрицы
{
	mass = new tVector[n];	//массив указателей на ячейки памяти
	for (size_t i = 0; i < n; i++)
		mass[i] = new datatype[n];
}
void EnteringArray(tMatrix &mass, size_t n) //ввод значений матрицы
{
	//ofstream IN("Data.txt");
	//if (mass = nullptr) { CreatingArray(mass, n); }
	for (int row = 0; row < n; row++)  //запись значений в массив
	{
		for (int col = 0; col < n; col++)
		{
			cout << "Введите элемент " << row + 1 << " строки " << col + 1 << " столбца: ";
			cin >> mass[row][col];
			//		IN << mass[row][col] << " ";
		}
		cout << endl;
		//	IN << endl;
	}
	cout << endl;
	//IN.close();
}
void PrintingArray(tMatrix &mass, size_t n)//вывод матрицы
{
	for (size_t row = 0; row < n; row++)  //вывод значений массива
	{
		for (size_t col = 0; col < n; col++)
		{
			cout << mass[row][col] << "  ";
		}
		cout << endl;
	}
	cout << endl;
}
void DeletingArray(tMatrix &Matrix, size_t N)//высвобождение памяти матрицы
{
	for (int count = 0; count < N; count++)  //высвобождение памяти
	{
		delete[] Matrix[count];
	}
	delete[] Matrix;
}
void DeletingArray(tMatrix &Matrix)//высвобождение памяти матрицы
{
	for (int count = 0; count < N; count++)  //высвобождение памяти
	{
		delete[] Matrix[count];
	}
	delete[] Matrix;
}
tMatrix Summ(tMatrix &mass1, tMatrix &mass2, size_t n)//матричная сумма
{
	tMatrix a = nullptr;
	CreatingArray(a, n);
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			a[i][j] = mass1[i][j] + mass2[i][j];
		}
	}
	return a;
}
tMatrix Difference(tMatrix &mass1, tMatrix &mass2, size_t n)//матричная разность
{

	tMatrix a = nullptr;
	CreatingArray(a, n);
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n; j++)

		{
			a[i][j] = mass1[i][j] - mass2[i][j];
		}
	}
	return a;
}
tMatrix Multiplication(tMatrix &mass1, tMatrix &mass2, size_t n)//матричное перемножение
{
	tMatrix a = nullptr;
	CreatingArray(a, n);
	datatype result;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			result = 0;
			for (int k = 0; k < n; k++)
			{
				result = result + mass1[i][k] * mass2[k][j];
			}
			a[i][j] = result;
		}
	}

	return a;
}
tMatrix Multiplication(tMatrix &mass1, tMatrix &mass2)//матричное перемножение
{
	tMatrix a = nullptr;
	CreatingArray(a, N);
	datatype result;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			result = 0;
			for (int k = 0; k < N; k++)
			{
				result = result + mass1[i][k] * mass2[k][j];
			}
			a[i][j] = result;
		}
	}

	return a;
}
tMatrix MatrCopy(tMatrix &mass, size_t n)
{
	tMatrix A = nullptr;
	CreatingArray(A, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			A[i][j] = mass[i][j];
		}
	return A;
}
tMatrix MatrCopy(tMatrix &mass)
{
	tMatrix A = nullptr;
	CreatingArray(A, N);
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		{
			A[i][j] = mass[i][j];
		}
	return A;
}
tMatrix Transpose(tMatrix &A, size_t n)
{
	tMatrix Copy = new tVector[n];
	for (int j = 0; j < n; j++)
	{
		Copy[j] = new datatype[n];
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Copy[i][j] = A[j][i];
		}
	}
	return Copy;
}
tMatrix Transpose(tMatrix &A)
{
	tMatrix Copy = new tVector[N];
	for (int j = 0; j < N; j++)
	{
		Copy[j] = new datatype[N];
	}
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			Copy[i][j] = A[j][i];
		}
	}
	return Copy;
}
void DeletingVector(tVector &vect)//высвобождение памяти вектора
{
	//cout<<"+"<<endl;
	delete[] vect;
	//cout<<"+"<<endl;
}
void CreatingVector(tVector &vect, size_t n)//выделение памяти для вектора
{
	vect = new datatype[n];

}
void CreatingVector(tVector &vect)//выделение памяти для вектора
{
	vect = new datatype[N];
}
void EnteringVector(tVector &vect, size_t n)//ввод вектора
{

	for (int row = 0; row < n; row++)
	{
		cout << "Введите " << row + 1 << " элемент столбца B: ";
		cin >> vect[row];
	}

	cout << endl;
}
void PrintingVector(tVector &vect)//вывод вектора
{

	for (size_t row = 0; row < N; row++)
	{
		cout << vect[row] << endl;
	}

	cout << endl;
}

void PrintingVector(tVector &vect, int n)
{
	for (size_t row = 0; row < n; row++)
		{
			cout << vect[row] << endl;
		}

		cout << endl;
}
tVector Summ(tVector &vect1, tVector &vect2, int n)
{
	tVector result = new datatype[n];
//	CreatingVector(result, n);
	for (int i = 0; i < n; i++)
	{
		result[i] = vect1[i] + vect2[i];
	}
	return result;
}
tVector Difference(tVector &vect1, tVector &vect2, size_t n)
{
	tVector result = nullptr;
	CreatingVector(result, n);
	for (size_t i = 0; i < n; i++)
	{
		result[i] = vect1[i] - vect2[i];
	}
	return result;
}
datatype Scalar(tVector &vect1, tVector &vect2)
{
	datatype result = 0;
	//CreatingVector(result, n);
	for (size_t i = 0; i < N; i++)
	{
		result = result + (vect1[i] * vect2[i]);
	}
	return result;
}
tVector VectCopy(tVector &vect, size_t n)
{
	tVector A = nullptr;
	CreatingVector(A, n);
	for (int i = 0; i < n; i++)
	{
		A[i] = vect[i];
	}

	return A;
}
tVector Gauss(tMatrix &A, tVector &b)
{
	tMatrix AA = MatrCopy(A, N);
	tVector B = VectCopy(b, N);
	tVector Sol = nullptr;
	//int row_max;
	for (int j = 0; j < N; j++)
	{
		int row_max = j;
		datatype max = abs(AA[j][j]);
		for (int i = j; i < N; i++)
		{
			if (abs(AA[i][j]) > max)
			{
				max = abs(AA[i][j]);
				row_max = i;
			}
			//else row_max = j;
		}//поиск наибольшего по модулю значения в столбце


		if (row_max != j)
		{
			datatype boofer;
			tVector show = nullptr;
			show = AA[j];
			AA[j] = AA[row_max];
			AA[row_max] = show;
			boofer = B[j];
			B[j] = B[row_max];
			B[row_max] = boofer;
		}

		if (abs(AA[j][j]) < eps7)
		{
			AA[j][j] = 0;
			cout << "Матрица системы вырожденная, система несовместна" << endl;
			tVector nul = nullptr;
			return nul;
		}
		else
		{
			for (int i = j + 1; i < N; i++)
			{
				if (abs(AA[i][j]) > eps)
				{
					datatype c;
					c = AA[i][j] / AA[j][j];
					for (int k = j; k < N; k++) //работа со строками ниже топовой
					{
						AA[i][k] = AA[i][k] - c * AA[j][k];
					}
					B[i] = B[i] - c * B[j];
				}
				else AA[i][j] = 0;
			}//создание нужных строк, в которых на месте эл-та A[i][j] стоит "0"

			/*datatype koeff = AA[j][j];
			B[j] = B[j] / koeff;
			for (int i = j; i < n; i++)
			{
			AA[j][i] = AA[j][i] / koeff;
			}//получение строки с элеметом "1" на главной диагонали */
		}
	}
	//PrintingArray(AA, n);
	Sol = new datatype[N];
	Sol[N - 1] = B[N - 1] / AA[N - 1][N - 1];
	for (int i = N - 2; i >= 0; i--)
	{
		datatype summa = 0;
		for (size_t j = N - 1; j > i; j--)
		{
			summa = summa + AA[i][j] * Sol[j];
		}
		Sol[i] = (B[i] - summa) / AA[i][i];
	}
	return Sol;
	delete[] Sol;
}

tVector Gauss(tMatrix &A, tVector &b, int N)
{
	tMatrix AA = MatrCopy(A, N);
	tVector B = VectCopy(b, N);
	tVector Sol = nullptr;
	//int row_max;
	for (int j = 0; j < N; j++)
	{
		int row_max = j;
		datatype max = abs(AA[j][j]);
		for (int i = j; i < N; i++)
		{
			if (abs(AA[i][j]) > max)
			{
				max = abs(AA[i][j]);
				row_max = i;
			}
			//else row_max = j;
		}//поиск наибольшего по модулю значения в столбце


		if (row_max != j)
		{
			datatype boofer;
			tVector show = nullptr;
			show = AA[j];
			AA[j] = AA[row_max];
			AA[row_max] = show;
			boofer = B[j];
			B[j] = B[row_max];
			B[row_max] = boofer;
		}

		if (abs(AA[j][j]) < eps7)
		{
			AA[j][j] = 0;
			cout << "Матрица системы вырожденная, система несовместна" << endl;
			tVector nul = nullptr;
			return nul;
		}
		else
		{
			for (int i = j + 1; i < N; i++)
			{
				if (abs(AA[i][j]) > eps)
				{
					datatype c;
					c = AA[i][j] / AA[j][j];
					for (int k = j; k < N; k++) //работа со строками ниже топовой
					{
						AA[i][k] = AA[i][k] - c * AA[j][k];
					}
					B[i] = B[i] - c * B[j];
				}
				else AA[i][j] = 0;
			}//создание нужных строк, в которых на месте эл-та A[i][j] стоит "0"

			/*datatype koeff = AA[j][j];
			B[j] = B[j] / koeff;
			for (int i = j; i < n; i++)
			{
			AA[j][i] = AA[j][i] / koeff;
			}//получение строки с элеметом "1" на главной диагонали */
		}
	}
	//PrintingArray(AA, n);
	Sol = new datatype[N];
	Sol[N - 1] = B[N - 1] / AA[N - 1][N - 1];
	for (int i = N - 2; i >= 0; i--)
	{
		datatype summa = 0;
		for (size_t j = N - 1; j > i; j--)
		{
			summa = summa + AA[i][j] * Sol[j];
		}
		Sol[i] = (B[i] - summa) / AA[i][i];
	}
	return Sol;
	delete[] Sol;
}

/*for (int i = 0; i < n; i++)
{
Sol[i] = 0;
}*/
datatype Determinant(tMatrix &A1, tVector &B, size_t n) //приведение матрицы к верхнегреулольному виду
{														//+ изменение вектора В
	datatype det = 1;
	size_t beta = 0;//счетчик числа перестановок строк


	for (int j = 0; j < n; j++)
	{
		int row_max = j;
		datatype max = abs(A1[j][j]);
		for (int i = j; i < n; i++)
		{
			if (abs(A1[i][j]) > max)
			{
				max = abs(A1[i][j]);
				row_max = i;
			}
			//else row_max = j;
		}//поиск наибольшего по модулю значения в столбце


		if (row_max != j)
		{
			beta = beta + 1;
			datatype boofer;
			tVector show = nullptr;
			show = A1[j];
			A1[j] = A1[row_max];
			A1[row_max] = show;
			boofer = B[j];
			B[j] = B[row_max];
			B[row_max] = boofer;
		}//перетаскивание строки с наибольшим элементом на топ место

		if (abs(A1[j][j]) < eps)
		{
			A1[j][j] = 0;
			cout << "Матрица системы вырожденная, система несовместна" << endl;
			/*tVector nul = nullptr;
			return nul;*/
		}
		else
		{
			for (int i = j + 1; i < n; i++)
			{
				datatype c;
				c = A1[i][j] / A1[j][j];
				for (int k = j; k < n; k++)  //работа со строками ниже топовой
				{
					A1[i][k] = A1[i][k] - c * A1[j][k];
				}
				B[i] = B[i] - c * B[j];
			}//создание нужных строк, в которых на месте эл-та A[i][j] стоит "0"
		}//исходно было то, что закомменчено снизу(для справки, если что-то не будет правильно работать)

	}

	for (int j = 0; j < n; j++)
	{
		det = det * A1[j][j];
	}
	if (beta % 2 == 1)
	{
		det = -1 * det;
	}


	return det;
}
datatype NormOfVector(tVector &vect1, size_t n) //||X||_1
{
	datatype Norm = 0;
	for (size_t i = 0; i < n; i++)
	{
		Norm = Norm + abs(vect1[i])*abs(vect1[i]);
	}
	Norm = sqrt(Norm);
	//cout << "norm=" << Norm << endl;
	return Norm;
}
datatype NormOfVectorC(tVector &vect1, size_t n)
{
	datatype norm = abs(vect1[0]);
	for (int i = 1; i<n; i++)
	{
		datatype a= abs(vect1[i]);
		if (a>norm) norm =a;
	}
	return norm;
}
tVector MVMultiplication(tMatrix &mass1, tVector &vect, int n)//перемножение вектора на матрицу (слева)
{
	tVector a = nullptr;
	CreatingVector(a, n);
	datatype result=0;
	for (int i = 0; i < n; i++)
	{
		result = 0;
		for (int k = 0; k < n; k++)
		{
			result = result + mass1[i][k] * vect[k];
		}
		a[i] = result;
	}

	return a;
}
void ResultsSave(tMatrix &A, tVector &B, tVector &X, size_t n, datatype iter) //сохранение результата в файл
{
	char *FileName = new char[50];
	cout << "Введите имя файла, в который запишется результат решения СЛАУ: ";
	cin.getline(FileName, 49);
	ofstream fout(FileName);
	fout << "A:";
	for (size_t row = 0; row < n; row++)  //вывод значений массива
	{
		for (size_t col = 0; col < n; col++)
		{
			if (row > 0 && col == 0)
				fout << "  " << A[row][col] << "  ";
			else fout << A[row][col] << "  ";
		}
		fout << endl;
	}
	fout << endl;
	fout << "x:";
	for (size_t row = 0; row < n; row++)
	{
		if (row > 0)fout << "  " << X[row] << endl;
		else fout << X[row] << endl;
	}

	fout << endl;
	fout << "b:";
	for (size_t row = 0; row < n; row++)
	{
		if (row > 0) fout << "  " << B[row] << endl;
		else fout << B[row] << endl;
	}

	fout << endl;

	tVector B1 = MVMultiplication(A, X, n);
	fout << "Ax=: ";
	for (size_t row = 0; row < n; row++)
	{
		if (row > 0) fout << "   " << B1[row] << endl;
		else fout << B1[row] << endl;
	}
	fout << "Невязка: ";
	for (size_t row = 0; row < n; row++)
	{
		if (row > 0) fout << "         " << B[row] - B1[row] << endl;
		else fout << B[row] - B1[row] << endl;
	}
	fout << "Точность: " << eps << endl;
	fout << "Шагов для сходимости: " << iter << endl;

	fout << endl;
	fout.close();
	delete[] FileName;
	DeletingVector(B1);
}
void ResultsSavetest(tMatrix &A, size_t n)
{
	char *FileName = new char[50];
	cout << "Введите имя файла, в который запишется результат : ";
	cin.getline(FileName, 49);
	ofstream fout(FileName);
	fout << "A:";
	for (size_t row = 0; row < n; row++)  //вывод значений массива
	{
		for (size_t col = 0; col < n; col++)
		{
			if (row > 0 && col == 0)
				fout << "  " << A[row][col] << "  ";
			else fout << A[row][col] << "  ";
		}
		fout << endl;
	}

	fout.close();
	delete[] FileName;



}
tMatrix Inverse(tMatrix &A, size_t n)
{
	tMatrix InvA = nullptr;
	tMatrix AA = MatrCopy(A, n);

	CreatingArray(InvA, n);


	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			if (i == j) InvA[i][j] = 1; else InvA[i][j] = 0;
		}

	}//инициализация единичной матрицы (its working)


	for (int j = 0; j < n; j++)
	{
		int row_max = j;
		datatype max = abs(AA[j][j]);
		for (int i = j; i < n; i++)
		{
			if (abs(AA[i][j]) > max)
			{
				max = abs(AA[i][j]);
				row_max = i;
			}
		}//поиск наибольшего по модулю значения в столбце


		if (row_max != j)
		{
			tVector show = nullptr;
			show = AA[j];
			AA[j] = AA[row_max];
			AA[row_max] = show;

			show = nullptr;
			show = InvA[j];
			InvA[j] = InvA[row_max];
			InvA[row_max] = show;
			show = nullptr;
			DeletingVector(show);
		}//перетаскивание строки с наибольшим элементом на топ место

		for (int i = j + 1; i < n; i++)
		{
			datatype c;
			c = AA[i][j] / AA[j][j];
			for (int k = 0; k < n; k++)  //работа со строками ниже топовой
			{

				AA[i][k] = AA[i][k] - c * AA[j][k];
				InvA[i][k] = InvA[i][k] - c * InvA[j][k];
			}//cout << "==================================" << endl;
		}//создание нужных строк, в которых на месте эл-та A[i][j] стоит "0"



		if (abs(AA[j][j]) < eps7)
		{
			AA[j][j] = 0;
			//cout << "Матрица системы вырожденная, обратной не существует" << endl;
			tMatrix nul = nullptr;
			return nul;

		}
		else
		{
			datatype koeff = AA[j][j];
			for (int i = j; i < n; i++)
			{
				AA[j][i] = AA[j][i] / koeff;
			}//получение строки с элеметом "1" на главной диагонали
			for (int i = 0; i < n; i++)
			{
				InvA[j][i] = InvA[j][i] / koeff;
			}
		}
		//cout << "Iteratoin " << j + 1 << " : " << endl;
		//PrintingArray(AA, n);
		//cout << endl;
		//	PrintingArray(InvA, n);
	}
	//cout << "левая матрица после спуска вниз" << endl;
	//PrintingArray(AA, n);
	//cout << "правая матрицы после спуска вниз" << endl;
	//PrintingArray(InvA, n);


	for (int j = n - 1; j >= 1; j--)
	{
		for (int i = j - 1; i >= 0; i--)
		{
			datatype c = AA[i][j];
			AA[i][j] = 0;
			for (int k = 0; k < n; k++)
			{
				InvA[i][k] = InvA[i][k] - c * InvA[j][k];
			}
		}

	}
	DeletingArray(AA, n);
	//ResultsSavetest(InvA, n);
	return InvA;
}

tVector QRmethod(tMatrix &mass1, tVector &vect) //QR-разложение
{
	//cout << "Поехали" << endl;
	tVector Sol = nullptr;
	tVector B = nullptr;
	tMatrix A = nullptr;
	tMatrix T = nullptr;
	datatype InvariantI = 0;
	datatype InvariantF = 0;

	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			InvariantI = InvariantI + mass1[i][j] * mass1[i][j];
		}
	}

	CreatingArray(T, N);//Матрица Т
	A = MatrCopy(mass1, N);//Определяем А как копию исходной матрицы
	B = VectCopy(vect, N);//Определяем B как копию исходного вектора правой части
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (i == j) T[i][j] = 1;
			else T[i][j] = 0;
		}//Изначально матрица T имеет вид единичной диагональной
	}

	for (int j = 0; j < N - 1; j++)
	{

		for (int i = j + 1; i < N; i++)
		{
			if (abs(A[i][j]) < eps)
			{
				A[i][j] = 0;
			}
			else
			{
				/*cout << "Поехали6" << endl;*/
				datatype c = A[j][j] / (sqrt(A[j][j] * A[j][j] + A[i][j] * A[i][j])); //Коэффициент с
				datatype s = A[i][j] / (sqrt(A[j][j] * A[j][j] + A[i][j] * A[i][j]));//Коэффициент s
				tMatrix T1 = nullptr;
				CreatingArray(T1, N);//Cоздаем промежуточную матрицу матрицу T1 в виде даиагональной единичной
				for (int k = 0; k < N; k++)
				{
					for (int l = 0; l < N; l++)
					{
						if (k == l) T1[k][l] = 1.;
						else T1[k][l] = 0.;
						/*cout << "Поехали7" << endl;*/
					}
				}//Изначально матрица T имеет вид единичной диагональной
				/*PrintingArray(T1, n);*/

				T1[j][j] = c;
				T1[j][i] = s;
				T1[i][j] = -s;
				T1[i][i] = c;
				/*cout << "Поехали8" << endl;
				PrintingArray(T1, n);*/

				//A = Multiplication(T1, A, n);//домножаем А на T1 слева
				//for (int k = 0; k < n; k++)//домножаем А на T1 слева
				//for (int k = j; k < n; k++)
				for (int k = 0; k < N; k++)
				{
					datatype IK = A[i][k];
					datatype JK = A[j][k];
					A[j][k] = T1[j][j] * JK + T1[j][i] * IK;
					A[i][k] = T1[i][j] * JK + T1[i][i] * IK;
				}
				B = MVMultiplication(T1, B, N);//домножаем В на T1 слева
				//T = Multiplication(T1, T, n);//домножаем Т на T1 слева
				//for (int k = j; k < n; k++)//домножаем T на T1 слева
				for (int k = 0; k < N; k++)
				{
					datatype IK = T[i][k];
					datatype JK = T[j][k];
					T[j][k] = T1[j][j] * JK + T1[j][i] * IK;
					T[i][k] = T1[i][j] * JK + T1[i][i] * IK;
				}
				//PrintingArray(A, n);
				DeletingArray(T1, N);//удаляем промежуточную Т1
			}
		}
	}

	//ОБРАТНЫЙ ХОД
	//PrintingArray(A);
	for (int y = 0; y < N; y++)
	{
		if (abs(A[y][y]) < eps7)
		{
			A[y][y] = 0;
			cout << "Матрица системы вырожденная, система несовместна" << endl;
			tVector nul = nullptr;
			return nul;
		}
	}//проверка на вырожденность

	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			if (i != j)
			{
				if (abs(A[i][j]) < eps) A[i][j] = 0;
			}

		}
	}
	//PrintingArray(A);
	CreatingVector(Sol, N); //матрица А = R
	datatype opa;				//матрица В = b*
	Sol[N - 1] = B[N - 1] / A[N - 1][N - 1];
	for (int l = (N - 1); l >= 1; l--)
	{
		opa = B[l - 1];
		for (int k = (l + 1); k <= N; k++) opa = opa - Sol[k - 1] * A[l - 1][k - 1];
		Sol[l - 1] = opa / A[l - 1][l - 1];
	}

	/*	+
		+
		+
		+
													 НА ЭТОМ РЕШЕНИЕ ОКОНЧЕНО!! ДАЛЬШЕ ДОП ВОПРОСЫ!!
		+
		+
		+
		+
	*/

	/*cout << "Матрица R (TA матрица)" << endl;
	PrintingArray(A, N);*/

	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			InvariantF = InvariantF + A[i][j] * A[i][j];
		}
	}

	//cout << "Матрица Q (транспонированная матрица T):" << endl;
	tMatrix TT = nullptr;;
	TT = Transpose(T, N);
	tMatrix QR = nullptr;
	tMatrix P = nullptr;


	QR = Multiplication(TT, A, N);
	/*cout << "QR:" << endl;
	PrintingArray(QR, N);
	cout << "Difference of Matrix QR:" << endl;*/
	P = Difference(mass1, QR, N);
	//PrintingArray(P, N);
	DeletingArray(QR, N);
	DeletingArray(P, N);
	//cout << "Матрица Q (транспонированная матрица T):" << endl;
	//PrintingArray(TT, N);
	DeletingArray(TT, N);
	//PrintingArray(T, N);
	/*if (abs(InvariantF - InvariantI)<0.01)
	{
		cout << "Начало: " << InvariantI << endl;
		cout << "Конец: " << InvariantF << endl;
		cout << "Сумма квадратов элементов матрицы - инвариант" << endl;
	}*/

	//PrintingVector(Sol);
	return Sol;
}
datatype ConditionNumberOne(tMatrix &A, size_t n)  //max via columns
{
	tMatrix AI = Inverse(A, n);
	datatype cond1 = -1;
	datatype cond2 = -1;
	datatype sum;

	for (int j = 0; j < n; j++)  // ||A||_1
	{
		sum = 0;
		for (size_t i = 0; i < n; i++)
		{
			sum = sum + abs(A[i][j]);
		}
		if (sum > cond1) cond1 = sum;
	}

	for (int j = 0; j < n; j++)  // ||AInversed||_1
	{
		sum = 0;
		for (size_t i = 0; i < n; i++)
		{
			sum = sum + abs(AI[i][j]);
		}
		if (sum > cond2) cond2 = sum;
	}

	/*cout << "condA(1)=" << cond1*cond2 << endl;*/
	DeletingArray(AI, n);
	return cond1 * cond2;
}
datatype ConditionNumberInf(tMatrix &A, size_t n)  //max via rows   //ГОВНО ЕБАНОЕ, НАДО ДЕЛАТЬ ПО-ДРУГОМУ
{
	tMatrix AI = Inverse(A, n);
	datatype cond1 = -1;
	datatype cond2 = -1;
	datatype sum;

	for (int i = 0; i < n; i++)  // ||A||_1
	{
		sum = 0;
		for (size_t j = 0; j < n; j++)
		{
			sum = sum + abs(A[i][j]);
		}
		if (sum > cond1) cond1 = sum;
	}

	for (int i = 0; i < n; i++)  // ||AInversed||_1
	{
		sum = 0;
		for (size_t j = 0; j < n; j++)
		{
			sum = sum + abs(AI[i][j]);
		}
		if (sum > cond2) cond2 = sum;
	}

	/*cout << "condA(Inf)=" << cond1 * cond2 << endl;*/
	DeletingArray(AI, n);
	return cond1 * cond2;
}
datatype NormOfDifferenceB(tMatrix &A, tVector &B, tVector &X, size_t n)
{
	tVector B1 = MVMultiplication(A, X, n);
	datatype NOD = 0;
	for (size_t i = 0; i < n; i++)
	{
		NOD = NOD + (B1[i] - B[i])*(B1[i] - B[i]);
	}
	NOD = sqrt(NOD);
	DeletingVector(B1);
	return NOD;
}
void StabilityOfSolve(tMatrix &A, tVector &B, tVector &X1, size_t n)
{
	tVector B1 = VectCopy(B, n);
	for (size_t i = 0; i < n; i++)
	{
		B1[i] = B1[i] - 0.01;
	}
	tVector X2 = Gauss(A, B1);
	DeletingVector(B1);
	tVector DX = nullptr;
	DX = Difference(X1, X2, n);
	if (NormOfVector(DX, n) < 0.01)
	{
		cout << "Решение устойчиво к погрешностям исходных данных" << endl;
	}
	else
	{
		cout << "Решение НЕустойчиво к погрешностям исходных данных" << endl;
	}
} //
void EstimateOfCond(tMatrix &A, tVector &B, tVector &X, size_t n)
{//ОЦЕНКА СНИЗУ
	datatype NormB = NormOfVector(B, n);
	datatype NormX = NormOfVector(X, n);
	datatype NormDX = 0;
	datatype const delta = 0.001;
	datatype EOC = 1e+10;

	tVector B1 = nullptr;
	CreatingVector(B1, n);

	tVector DB = nullptr;
	CreatingVector(DB, n);

	for (size_t j = 0; j < 50; j++)  //50 итераций
	{
		srand(time(0));

		for (size_t i = 0; i < n; i++)
		{
			DB[i] = delta * (5 + rand() % 6); //генерация погрешности от 0,005 до 0,01 вектора В
		}
		datatype NormDB = NormOfVector(DB, n); //Вычисление нормы вектора погрешности
		for (size_t i = 0; i < n; i++)
		{
			B1[i] = B[i] + DB[i];		//создание нового вектора В
		}

		tVector Solve = Gauss(A, B1); //получаем решение, отличающееся от Х
		tVector DX = Difference(Solve, X, n);  //вектор погрешности вычисления Х

		DeletingVector(Solve);
		NormDX = NormOfVector(DX, n);
		DeletingVector(DX);

		datatype bred = (NormDX / NormX) / (NormDB / NormB);
		if (bred < EOC)
		{
			EOC = bred;;
		}
	}
	DeletingVector(B1);
	DeletingVector(DB);
	cout << "Condition number of this matrix >= " << EOC << endl;
}

void EnteringFromFile(tMatrix &A, tVector &B)//вывод из файла матрицы А и вектора B.
{
	char *FileName = new char[50];
	cout << "Введите имя файла, из которого будет производиться запись массива А " << endl << "и вектора-столбца В: ";
	//gets_s(FileName, 49);

	int i = 0;
	char *str = new char[501];
	ifstream OUT(FileName);
	while (!OUT.eof())
	{
		OUT.getline(str, 500, '\n');
		i++;
	}
	OUT.close();
	OUT.open(FileName);
	N = i;
	cout << "Размерность массива: " << N << endl << endl;
	CreatingArray(A, N);
	CreatingVector(B, N);
	for (int k = 0; k < N; k++)
	{
		for (int j = 0; j < N; j++)
		{
			OUT >> A[k][j];
		}
		OUT >> B[k];
	}
	OUT.close();
	delete[] str;
	delete[] FileName;
}


void CreatingDataFile()  //работает с косяком: последняя и предпоследняя строчки заполняются неправильно
{
	int n = 213;
	ofstream IN("DATA213.txt");

	IN << 4 << " " << 1 << " ";
	for (size_t j = 2; j < n; j++)
	{
		IN << 0 << " ";
	}
	IN << 6 << endl; //инициальзация первой строки

	for (int i = 1; i < n; i++)//строка
	{
		for (int j = 0; j < n; j++)//столбец
		{
			if (abs(i - j) == 1 & i != n - 1 & j != n - 1)
			{
				IN << 1 << " ";
			}
			else
			{
				if (i == j)
				{
					IN << 4 << " ";
				}
				else
				{
					IN << 0 << " ";
				}
			}

		}
		if (i != n - 1)
		{
			IN << 10 - 2 * ((i + 1) % 2) << endl;
		}
		else
		{
			IN << 6;
		}

	}
	IN.close();
}



datatype NormOfMatrix(tMatrix &A) //||A||_1
{
	datatype NormMatrix = 0;


	for (size_t j = 0; j < N; j++)
	{
		datatype sum = 0;
		for (size_t i = 0; i < N; i++)
		{
			sum = sum + abs(A[i][j]);
		}
		if (NormMatrix < sum)
		{
			NormMatrix = sum;
		}
	}
	return NormMatrix;
}


tMatrix NumberMatrixMultiplication(datatype num, tMatrix &A)
{
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			A[i][j] = A[i][j] * num;
		}
	}
	return A;
}

tMatrix NumberMatrixMultiplication(datatype num, tMatrix &A, int n)
{
	tMatrix AA=nullptr;
	CreatingArray(AA, n);
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			AA[i][j] = A[i][j] * num;
		}
	}
	return AA;
}

void NMMultiplication(datatype num, tMatrix &A, int n)
{
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			A[i][j] = A[i][j] * num;
		}
	}
}


tVector NumberVectorMultiplication(datatype num, tVector &mass)
{
	tVector A = VectCopy(mass, N);
	for (size_t j = 0; j < N; j++)
	{
		A[j] = A[j] * num;
	}

	return A;
}

tVector NumberVectorMultiplication(datatype num, tVector &mass, int n)
{
	//tVector A = VectCopy(mass, n);
	tVector A = new datatype[n];
	for (int j = 0; j < n; j++)
	{
		A[j] = mass[j] * num;
	}

	return A;
}


tVector SimpleIterationMethod(tMatrix &a, tVector &b)
{
	tMatrix A = MatrCopy(a, N);
	tVector B = VectCopy(b, N);

	//PrintingArray(A, N);

	datatype LambdaMax = 0;
	for (size_t i = 0; i < N; i++)
	{
		datatype RowSum = 0;
		for (size_t j = 0; j < N; j++)
		{
			RowSum = RowSum + abs(A[i][j]);
		}
		if (A[i][i] < 0)
		{
			RowSum = RowSum + A[i][i];
		}
		if (LambdaMax <= RowSum)
		{
			LambdaMax = RowSum;
		}
	}//Оценка наибольшего собственного числа матрицы через круги Гершгорина
	//cout << LambdaMax;
	datatype t = 1 / LambdaMax; //при этом параметре итерационный метод сойдётся


	tMatrix C, E = nullptr;
	tVector X1 = nullptr;
	tVector X0 = nullptr;

	CreatingArray(C, N);
	CreatingArray(E, N);


	CreatingVector(X0, N);
	for (size_t i = 0; i < N; i++)
	{
		X0[i] = 1;
	}//инициальзиация вектора Х0 единицами


	B = NumberVectorMultiplication(t, B);//задание tB для итераций

	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			if (i == j)
			{
				E[i][j] = 1;
			}
			else
			{
				E[i][j] = 0;
			}
		}
	}//initialization E

	//C = -(t*A - E);
	A = NumberMatrixMultiplication(t, A);//tA
	C = Difference(E, A, N);// -(tA-E)

	bool flag = 1;
	tVector CX0 = nullptr;
	tVector Ax_k = nullptr;
	tVector r_0 = nullptr;
	int iterations = 0;


	//X1 = C * X0 + B;
	while (flag)
	{
		iterations++;
		CX0 = MVMultiplication(C, X0, N);//C*X0
		X1 = Summ(CX0, B, N); //C*X0+t*B

		Ax_k = MVMultiplication(a, X1, N);
		r_0 = Difference(Ax_k, b, N);        //критерий остановы |Ax-b|<eps

		if (NormOfVector(r_0, N) > eps)
		{
			X0 = X1;
		}
		else
		{
			flag = 0;
			DeletingVector(CX0);
		}
	}

	ResultsSave(a, b, X1, N, iterations);

	cout << "Метод простой итерации сошелся за " << iterations << " итераций" << endl;
	DeletingArray(E, N);
	DeletingArray(C, N);
	DeletingArray(A, N);
	DeletingVector(B);
	DeletingVector(Ax_k);
	DeletingVector(r_0);
	DeletingVector(X0);
	return X1;


}


tVector JacobiMethod(tMatrix &a, tVector &b)   //ГАВНО, 2 вариант рабочий
{
	tMatrix A = MatrCopy(a, N);
	tVector B = VectCopy(b, N);
	tMatrix C, D = nullptr;
	tVector X1 = nullptr;
	tVector X0 = nullptr;

	CreatingArray(C, N);
	CreatingArray(D, N);

	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			if (i == j)
			{
				D[i][i] = -1 / a[i][i];
			}
			else
			{
				D[i][j] = 0;
			}
		}
	}//задаём матрицу -D^{-1}




	CreatingVector(X0, N);
	for (size_t i = 0; i < N; i++)
	{
		X0[i] = 1;
	}//инициальзиация вектора Х0 единицами


	for (size_t i = 0; i < N; i++)
	{
		B[i] = -1 * B[i] * D[i][i];
	}//задание D^{-1}B для итераций

	for (size_t i = 0; i < N; i++)
	{
		A[i][i] = 0;
	}//L+U


	//C = -D^{-1}(L+U);
	C = Multiplication(D, A, N);// -D^{-1}(L+U)



	bool flag = 1;
	tVector CX0 = nullptr;
	tVector Ax_k = nullptr;
	tVector r_0 = nullptr;
	int iterations = 0;


	//X1 = C * X0 + B;
	while (flag)
	{
		iterations++;
		CX0 = MVMultiplication(C, X0, N);//C*X0
		X1 = Summ(CX0, B, N); //C*X0+t*B
		//cout << "CHECK" << endl;

		Ax_k = MVMultiplication(a, X1, N);
		r_0 = Difference(Ax_k, b, N);        //критерий остановы |Ax-b|<eps

		if (NormOfVector(r_0, N) > eps)
		{
			X0 = X1;
		}
		else
		{
			flag = 0;
			DeletingVector(CX0);
		}
	}

	ResultsSave(a, b, X1, N, iterations);

	cout << "Метод Якоби сошелся за " << iterations << " итераций" << endl;
	DeletingArray(D, N);
	DeletingArray(C, N);
	DeletingArray(A, N);
	DeletingVector(B);
	DeletingVector(Ax_k);
	DeletingVector(r_0);
	DeletingVector(X0);
	return X1;
}

tVector JacobiMethod2(tMatrix &a, tVector &b)
{
	tMatrix A = MatrCopy(a, N);
	tVector B = VectCopy(b, N);
	tMatrix C = nullptr;
	tVector D = nullptr;
	tVector X1 = nullptr;
	tVector X0 = nullptr;

	CreatingArray(C, N);
	CreatingVector(D, N);


	CreatingVector(X0, N);
	for (size_t i = 0; i < N; i++)
	{
		X0[i] = 1;
	}//инициальзиация вектора Х0 единицами


	for (size_t i = 0; i < N; i++)
	{
		B[i] = B[i] / A[i][i];
	}//задание D^{-1}B для итераций



	for (size_t i = 0; i < N; i++)//C = -D ^ { -1 }(L + U);
		for (size_t j = 0; j < N; j++)
		{
			if (i == j)
			{
				C[i][j] = 0;
			}
			else
			{
				C[i][j] = -A[i][j] / A[i][i];
			}
		}


	bool flag = 1;
	tVector CX0 = nullptr;
	tVector Ax_k = nullptr;
	tVector r_0 = nullptr;
	//tVector DX = nullptr;
	int iterations = 0;
	//datatype Epsilon = ((1 - NormOfMatrix(C))*eps) / NormOfMatrix(C);


	//X1 = C * X0 + B;
	while (flag)
	{
		iterations++;
		CX0 = MVMultiplication(C, X0, N);//C*X0
		X1 = Summ(CX0, B, N); //C*X0+t*B
		//cout << "CHECK" << endl;

		Ax_k = MVMultiplication(a, X1, N);
		r_0 = Difference(Ax_k, b, N);        //критерий остановы |Ax-b|<eps
		//DX = Difference(X1, X0, N);
		//PrintingVector(DX, N);
		//cout << endl;
		//PrintingVector(X1, N);

		if (NormOfVector(r_0, N) > eps)
		{
			X0 = X1;
		}
		else
		{
			flag = 0;
			DeletingVector(CX0);
		}
		//cin.get();
	}

	ResultsSave(a, b, X1, N, iterations);

	cout << "Метод Якоби сошелся за " << iterations << " итераций" << endl;
	DeletingVector(D);
	DeletingArray(C, N);
	DeletingArray(A, N);
	DeletingVector(B);
	DeletingVector(Ax_k);
	DeletingVector(r_0);
	//DeletingVector(DX);
	DeletingVector(X0);
	return X1;
}

tVector Relax(tMatrix &a, tVector &b)
{
	tMatrix A = MatrCopy(a, N);
	tVector B = VectCopy(b, N);
	//tMatrix C, D = nullptr;
	tVector X1 = nullptr;
	tVector X0 = nullptr;
	datatype omega = 0.5;//допустим
	CreatingVector(X0, N);
	for (size_t i = 0; i < N; i++)
	{
		X0[i] = 1;
	}//инициальзиация вектора Х0 единицами

	/*for (int i = 0; i < N; i++)
	{
		B[i] =  (b[i] / a[i][i]);
	}*///Это значит у нас ебать вектор у.
	CreatingVector(X1, N);
	int iterations = 0;
	bool flag = 1;
	while (flag)
	{
		iterations++;
		//cout << "iterations="<<iterations <<endl;
		tVector Ax_k = nullptr;
		tVector r_0 = nullptr;

		for (int i = 0; i < N; i++)
		{
			//cout << "i=" << i << endl;
			datatype penis1 = 0.0;
			datatype penis2 = 0.0;
			//int k = 0;
			for (int j = 0; j < i; j++)
			{
				//k++;
				penis1 += A[i][j] * X1[j];
			}
			//cout << "p1=" << penis1 << endl;
			for (int j = i + 1; j < N; j++)
			{
				penis2 += A[i][j] * X0[j];
			}
			//cout << "p2=" << penis2 << endl;
			X1[i] = omega * (B[i] - penis1 - penis2) / A[i][i] + (1 - omega)*X0[i];
			//cout << "x1[" <<i<<"]="<< X1[i]<< endl;
			//cout << endl;
		}
		//cout << "x1=";
		//PrintingVector(X1, N);
		//cout << endl;
		Ax_k = MVMultiplication(a, X1, N);
		r_0 = Difference(Ax_k, b, N);        //критерий остановы |Ax-b|<eps
		//cout << "r_0="<< NormOfVector(r_0, N) <<endl;
		if (NormOfVector(r_0, N) > eps)
		{
			X0 = X1;
		}
		else
		{
			flag = 0;
		}

		DeletingVector(Ax_k);
		DeletingVector(r_0);
		//cin.get();
	}



	//вычисление матрицы С=-(L+D)^{-1}*U
	tMatrix LPlusD = nullptr;
	LPlusD = MatrCopy(a, N);

	for (size_t i = 0; i < N - 1; i++)
	{
		for (size_t j = i + 1; j < N; j++)
		{
			LPlusD[i][j] = 0;
		}
	} //задание L+D матрицы
	/*cout << "L+D matrix: " << endl;
	PrintingArray(LPlusD, N);*/

	for (size_t i = 1; i < N; i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			LPlusD[i][j] = LPlusD[i][j] * omega;
		}
	} //задание wL+D матрицы
	/*cout << "wL+D matrix: " << endl;
	PrintingArray(LPlusD, N);*/

	LPlusD = Inverse(LPlusD, N);
	/*cout << "Inv(wL+D) matrix: " << endl;
	PrintingArray(LPlusD, N);*/

	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < i + 1; j++)
		{
			LPlusD[i][j] = -LPlusD[i][j] * omega;
		}
	} //задание -w(wL+D)^{-1} матрицы
	/*cout << "Inv(wL+D)*(-w) matrix: " << endl;
	PrintingArray(LPlusD, N);*/

	LPlusD = Multiplication(LPlusD, a, N);
	/*cout << "Inv(wL+D)*(-w)A matrix: " << endl;
	PrintingArray(LPlusD, N);*/


	for (size_t i = 0; i < N; i++)
	{
		LPlusD[i][i] = LPlusD[i][i] + 1;
	} //матрица С
	/*cout << "Inv(wL+D)*(-w)A+E matrix: " << endl;
	PrintingArray(LPlusD, N);*/

	/*cout << "L+D matrix: " << endl;
	PrintingArray(LPlusD, N);*/
	cout << "Матрица С: " << endl;
	PrintingArray(LPlusD, N);
	cout << "Норма ||C||_1 = " << NormOfMatrix(LPlusD) << endl;



	DeletingArray(LPlusD, N);



	ResultsSave(a, b, X1, N, iterations);
	DeletingArray(A, N);
	DeletingVector(B);
	DeletingVector(X0);

	return X1;

}





tVector ZeydelMethod(tMatrix &a, tVector &b)
{
	tMatrix A = MatrCopy(a, N);
	tVector B = VectCopy(b, N);
	//tMatrix C, D = nullptr;
	tVector X1 = nullptr;
	tVector X0 = nullptr;
	datatype omega = 1.0;//допустим


	CreatingVector(X0, N);
	for (size_t i = 0; i < N; i++)
	{
		X0[i] = 1;
	}//инициальзиация вектора Х0 единицами

	/*for (int i = 0; i < N; i++)
	{
		B[i] =  (b[i] / a[i][i]);
	}*///Это значит у нас ебать вектор у.
	CreatingVector(X1, N);
	int iterations = 0;
	bool flag = 1;
	while (flag)
	{
		iterations++;
		//cout << "iterations=" << iterations << endl;
		tVector Ax_k = nullptr;
		tVector r_0 = nullptr;

		for (int i = 0; i < N; i++)
		{
			//cout << "i=" << i << endl;
			datatype penis1 = 0.0;
			datatype penis2 = 0.0;
			//int k = 0;
			for (int j = 0; j < i; j++)
			{
				//k++;
				penis1 += A[i][j] * X1[j];
			}
			//cout << "p1=" << penis1 << endl;
			for (int j = i + 1; j < N; j++)
			{
				penis2 += A[i][j] * X0[j];
			}
			//cout << "p2=" << penis2 << endl;
			X1[i] = omega * (B[i] - penis1 - penis2) / A[i][i] + (1 - omega)*X0[i];
			//cout << "x1[" << i << "]=" << X1[i] << endl;
			//cout << endl;
		}
		//cout << "x1=";
		//PrintingVector(X1, N);
		//cout << endl;
		Ax_k = MVMultiplication(a, X1, N);
		r_0 = Difference(Ax_k, b, N);        //критерий остановы |Ax-b|<eps
		//cout << "r_0=" << NormOfVector(r_0, N) << endl;
		if (NormOfVector(r_0, N) > eps)
		{
			X0 = X1;
		}
		else
		{
			flag = 0;
		}

		DeletingVector(Ax_k);
		DeletingVector(r_0);
		//cin.get();
	}
	ResultsSave(a, b, X1, N, iterations);


	//вычисление матрицы С=-(L+D)^{-1}*U
	tMatrix LPlusD, U = nullptr;
	LPlusD = MatrCopy(a, N);

	U = MatrCopy(a, N);
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < i + 1; j++)
		{
			U[i][j] = 0;
		}
	}
	/*cout << "U matrix: " << endl;
	PrintingArray(U, N);*/

	for (size_t i = 0; i < N - 1; i++)
	{
		for (size_t j = i + 1; j < N; j++)
		{
			LPlusD[i][j] = 0;
		}
	} //задание L+D матрицы
	/*cout << "L+D matrix: " << endl;
	PrintingArray(LPlusD, N);*/

	LPlusD = Inverse(LPlusD, N);
	for (size_t i = 0; i < N - 1; i++)
	{
		for (size_t j = 0; j < N - i - 1; i++)
		{
			LPlusD[i][j] = -LPlusD[i][j];
		}
	} //задание -(L+D)^{-1} матрицы
	tMatrix S = MatrCopy(LPlusD, N);
	LPlusD = Multiplication(S, U, N);//LPlusD=C
	cout << "Матрица С: " << endl;
	PrintingArray(LPlusD, N);
	cout << "Норма ||C||_1 = " << NormOfMatrix(LPlusD) << endl;


	DeletingArray(S, N);
	DeletingArray(U, N);
	DeletingArray(LPlusD, N);
	DeletingArray(A, N);
	DeletingVector(B);
	DeletingVector(X0);

	return X1;
}



tVector Relax3Diag(tMatrix &a, tVector &b, datatype omega)
{
	tMatrix A = MatrCopy(a, N);
	tVector B = VectCopy(b, N);

	tVector X1 = nullptr;
	tVector X0 = nullptr;
	tVector Diag, LDiag, UDiag = nullptr;
	CreatingVector(Diag, N);
	CreatingVector(LDiag, N - 1);
	CreatingVector(UDiag, N - 1);

	for (size_t i = 0; i < N; i++)
	{
		Diag[i] = A[i][i];
	}
	for (size_t i = 0; i < N - 1; i++)
	{
		LDiag[i] = A[i + 1][i];
	}
	for (size_t i = 0; i < N - 1; i++)
	{
		UDiag[i] = A[i][i + 1];
	}

	CreatingVector(X0, N);
	for (size_t i = 0; i < N; i++)
	{
		X0[i] = 1;
	}//инициальзиация вектора Х0 единицами

	CreatingVector(X1, N);


	tVector Ax_k = nullptr;
	tVector r_0 = nullptr;
	CreatingVector(Ax_k, N);
	size_t iterations = 0;
	bool flag = 1;
	while (flag)
	{
		iterations++;


		for (int i = 0; i < N; i++)
		{

			datatype p1 = 0.0;
			datatype p2 = 0.0;
			if (i == 0)
			{
				p1 = 0;
			}
			else
			{
				p1 = LDiag[i - 1] * X1[i - 1];
			}

			if (i == N - 1)
			{
				p2 = 0;
			}
			else
			{
				p2 = UDiag[i] * X0[i + 1];
			}

			X1[i] = omega * (B[i] - p1 - p2) / A[i][i] + (1 - omega)*X0[i];

		}

		Ax_k[0] = Diag[0] * X1[0] + UDiag[0] * X1[1];
		for (size_t i = 1; i < N - 1; i++)
		{
			Ax_k[i] = X1[i - 1] * LDiag[i - 1] + Diag[i] * X1[i] + X1[i + 1] * UDiag[i];
		}

		Ax_k[N - 1] = LDiag[N - 2] * X1[N - 2] + Diag[N - 1] * X1[N - 1];


		r_0 = Difference(Ax_k, b, N);//критерий остановы |Ax-b|<eps


		if (NormOfVector(r_0, N) > eps)
		{
			X0 = X1;
		}
		else
		{
			flag = 0;
		}


		//cin.get();
	}
	DeletingVector(Ax_k);
	DeletingVector(r_0);
	ResultsSave(a, b, X1, N, iterations);
	cout << "Метод сошёлся за " << iterations << " итераций" << endl;


	//вычисление матрицы С=-(L+D)^{-1}*U
	tMatrix LPlusD = nullptr;
	LPlusD = MatrCopy(a, N);

	for (size_t i = 0; i < N - 1; i++)
	{
		for (size_t j = i + 1; j < N; j++)
		{
			LPlusD[i][j] = 0;
		}
	} //задание L+D матрицы
	/*cout << "L+D matrix: " << endl;
	PrintingArray(LPlusD, N);*/

	for (size_t i = 1; i < N; i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			LPlusD[i][j] = LPlusD[i][j] * omega;
		}
	} //задание wL+D матрицы
	/*cout << "wL+D matrix: " << endl;
	PrintingArray(LPlusD, N);*/

	LPlusD = Inverse(LPlusD, N);
	/*cout << "Inv(wL+D) matrix: " << endl;
	PrintingArray(LPlusD, N);*/

	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < i + 1; j++)
		{
			LPlusD[i][j] = -LPlusD[i][j] * omega;
		}
	} //задание -w(wL+D)^{-1} матрицы
	/*cout << "Inv(wL+D)*(-w) matrix: " << endl;
	PrintingArray(LPlusD, N);*/

	LPlusD = Multiplication(LPlusD, a, N);
	/*cout << "Inv(wL+D)*(-w)A matrix: " << endl;
	PrintingArray(LPlusD, N);*/


	for (size_t i = 0; i < N; i++)
	{
		LPlusD[i][i] = LPlusD[i][i] + 1;
	} //матрица С

	cout << "Норма ||C||_1 = " << NormOfMatrix(LPlusD) << endl;






	DeletingArray(LPlusD, N);
	DeletingVector(Diag);
	DeletingVector(LDiag);
	DeletingVector(UDiag);
	DeletingArray(A, N);
	DeletingVector(B);
	DeletingVector(X0);

	return X1;

}

tVector Relax3Diag(tMatrix &a, tVector &b, datatype omega, size_t N)
{
	tMatrix A = MatrCopy(a, N);
	tVector B = VectCopy(b, N);

	tVector X1 = nullptr;
	tVector X0 = nullptr;
	tVector Diag, LDiag, UDiag = nullptr;
	CreatingVector(Diag, N);
	CreatingVector(LDiag, N - 1);
	CreatingVector(UDiag, N - 1);

	for (size_t i = 0; i < N; i++)
	{
		Diag[i] = A[i][i];
	}
	for (size_t i = 0; i < N - 1; i++)
	{
		LDiag[i] = A[i + 1][i];
	}
	for (size_t i = 0; i < N - 1; i++)
	{
		UDiag[i] = A[i][i + 1];
	}

	CreatingVector(X0, N);
	for (size_t i = 0; i < N; i++)
	{
		X0[i] = 1;
	}//инициальзиация вектора Х0 единицами

	CreatingVector(X1, N);


	tVector Ax_k = nullptr;
	tVector r_0 = nullptr;
	CreatingVector(Ax_k, N);
	size_t iterations = 0;
	bool flag = 1;
	while (flag)
	{
		iterations++;


		for (int i = 0; i < N; i++)
		{

			datatype p1 = 0.0;
			datatype p2 = 0.0;
			if (i == 0)
			{
				p1 = 0;
			}
			else
			{
				p1 = LDiag[i - 1] * X1[i - 1];
			}

			if (i == N - 1)
			{
				p2 = 0;
			}
			else
			{
				p2 = UDiag[i] * X0[i + 1];
			}

			X1[i] = omega * (B[i] - p1 - p2) / A[i][i] + (1 - omega)*X0[i];

		}

		Ax_k[0] = Diag[0] * X1[0] + UDiag[0] * X1[1];
		for (size_t i = 1; i < N - 1; i++)
		{
			Ax_k[i] = X1[i - 1] * LDiag[i - 1] + Diag[i] * X1[i] + X1[i + 1] * UDiag[i];
		}

		Ax_k[N - 1] = LDiag[N - 2] * X1[N - 2] + Diag[N - 1] * X1[N - 1];


		r_0 = Difference(Ax_k, b, N);//критерий остановы |Ax-b|<eps


		if (NormOfVector(r_0, N) > eps)
		{
			X0 = X1;
		}
		else
		{
			flag = 0;
		}


		//cin.get();
	}
	DeletingVector(Ax_k);
	DeletingVector(r_0);
	ResultsSave(a, b, X1, N, iterations);
	cout << "Метод сошёлся за " << iterations << " итераций" << endl;


	//вычисление матрицы С=-(L+D)^{-1}*U
	tMatrix LPlusD = nullptr;
	LPlusD = MatrCopy(a, N);

	for (size_t i = 0; i < N - 1; i++)
	{
		for (size_t j = i + 1; j < N; j++)
		{
			LPlusD[i][j] = 0;
		}
	} //задание L+D матрицы
	/*cout << "L+D matrix: " << endl;
	PrintingArray(LPlusD, N);*/

	for (size_t i = 1; i < N; i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			LPlusD[i][j] = LPlusD[i][j] * omega;
		}
	} //задание wL+D матрицы
	/*cout << "wL+D matrix: " << endl;
	PrintingArray(LPlusD, N);*/

	LPlusD = Inverse(LPlusD, N);
	/*cout << "Inv(wL+D) matrix: " << endl;
	PrintingArray(LPlusD, N);*/

	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < i + 1; j++)
		{
			LPlusD[i][j] = -LPlusD[i][j] * omega;
		}
	} //задание -w(wL+D)^{-1} матрицы
	/*cout << "Inv(wL+D)*(-w) matrix: " << endl;
	PrintingArray(LPlusD, N);*/

	LPlusD = Multiplication(LPlusD, a, N);
	/*cout << "Inv(wL+D)*(-w)A matrix: " << endl;
	PrintingArray(LPlusD, N);*/


	for (size_t i = 0; i < N; i++)
	{
		LPlusD[i][i] = LPlusD[i][i] + 1;
	} //матрица С

	cout << "Норма ||C||_1 = " << NormOfMatrix(LPlusD) << endl;






	DeletingArray(LPlusD, N);
	DeletingVector(Diag);
	DeletingVector(LDiag);
	DeletingVector(UDiag);
	DeletingArray(A, N);
	DeletingVector(B);
	DeletingVector(X0);

	return X1;

}


tVector Shuttle(tVector &lowdiagg, tVector &diagg, tVector &updiagg, tVector &b, int m)
{
	tVector X = nullptr;
	tVector alpha = nullptr;
	tVector beta = nullptr; // коэффициенты прогонки
	tVector B = VectCopy(b, m);
	tVector L = VectCopy(lowdiagg, m - 1);//нижняя диагональ
	tVector Diag = VectCopy(diagg, m);//главная диагональ
	tVector U = VectCopy(updiagg, m - 1);//верхняя диагональ
	CreatingVector(X, m);
	CreatingVector(alpha, m - 1);
	CreatingVector(beta, m);

	alpha[0] = -U[0] / Diag[0];

	for (int i = 1; i < m - 1; i++)
	{
		alpha[i] = -U[i] / (L[i - 1] * alpha[i - 1] + Diag[i]); //
	}

	beta[0] = B[0] / Diag[0];//Задаем первый элемент бета

	for (int i = 1; i < m; i++)
	{
		beta[i] = (B[i] - L[i - 1] * beta[i - 1]) / (L[i - 1] * alpha[i - 1] + Diag[i]); //все остальные
	}

	X[m - 1] = beta[m - 1];

	for (int i = m - 2; i >= 0; i--)
	{
		X[i] = alpha[i] * X[i + 1] + beta[i]; //теперь все
	}


	DeletingVector(alpha);
	DeletingVector(beta);
	DeletingVector(B);
	DeletingVector(L);
	DeletingVector(Diag);
	DeletingVector(U);


	return X;
}




tMatrix Hessenberg(tMatrix &mass)
{
	tMatrix A = nullptr;
	A = MatrCopy(mass, N);
	PrintingArray(A);

	for (int k = 1; k < N; k++)
	{

		for (int l = k + 1; l < N; l++)
		{
			if (abs(A[l][k - 1]) < eps)
			{
				A[l][k - 1] = 0;
			}
			else
			{

				datatype alpha = A[k][k - 1] / (sqrt(A[k][k - 1] * A[k][k - 1] + A[l][k - 1] * A[l][k - 1])); //Коэффициент альфа
				datatype betta = A[l][k - 1] / (sqrt(A[k][k - 1] * A[k][k - 1] + A[l][k - 1] * A[l][k - 1]));//Коэффициент бета

				for (int p = 0; p < N; p++)
				{
					datatype IK = A[l][p];
					datatype JK = A[k][p];
					A[k][p] = alpha * JK + betta * IK;
					A[l][p] = -betta * JK + alpha * IK;
				}
				//PrintingArray(A);
				for (int p = 0; p < N; p++)
				{
					datatype IK = A[p][l];
					datatype JK = A[p][k];
					A[p][k] = alpha * JK + betta * IK;
					A[p][l] = -betta * JK + alpha * IK;
				}
				/*PrintingArray(A);
				cin.get();*/
			}
		}
		//cin.get();
	}
	for (int k = 0; k < N; k++)
	{
		for (int l = 0; l < N; l++)
		{
			if (abs(A[l][k]) < eps)
			{
				A[l][k] = 0;
			}
		}
	}




	PrintingArray(A);
	return A;
}

void EnteringFromFile(tMatrix &A)//вывод из файла матрицы А и вектора B.
{
	char FileName[] = { "1.TXT" };
	/*cout << "Введите имя файла, из которого будет производиться запись массива А: ";
	gets_s(FileName, 49);*/

	int i = 0;
	char *str = new char[501];
	ifstream OUT(FileName);
	while (!OUT.eof())
	{
		OUT.getline(str, 500, '\n');
		i++;
	}
	OUT.close();
	OUT.open(FileName);
	N = i;
	cout << "Размерность массива: " << N << endl << endl;
	CreatingArray(A, N);
	for (int k = 0; k < N; k++)
	{
		for (int j = 0; j < N; j++)
		{
			OUT >> A[k][j];
		}
	}
	OUT.close();
	delete[] str;
	//delete[] FileName;
}


void PrintingArray(tMatrix &mass)//вывод матрицы
{
	for (size_t row = 0; row < N; row++)  //вывод значений массива
	{
		for (size_t col = 0; col < N; col++)
		{
			cout << mass[row][col] << "  ";
		}
		cout << endl;
	}
	cout << endl;
}


//tMatrix QR_T(tMatrix &mass1)
//{
//	tMatrix T = nullptr;
//	tMatrix A = MatrCopy(mass1,N);
//	CreatingArray(T, N);//Матрица Т
//
//	for (int i = 0; i < N; i++)
//	{
//		for (int j = 0; j < N; j++)
//		{
//			if (i == j) T[i][j] = 1;
//			else T[i][j] = 0;
//		}//Изначально матрица T имеет вид единичной диагональной
//	}
//
//	for (int j = 0; j < N - 1; j++)
//	{
//
//		for (int i = j + 1; i < N; i++)
//		{
//			if (abs(A[i][j]) < eps)
//			{
//				A[i][j] = 0;
//			}
//			else
//			{
//				/*cout << "Поехали6" << endl;*/
//				datatype c = A[j][j] / (sqrt(A[j][j] * A[j][j] + A[i][j] * A[i][j])); //Коэффициент с
//				datatype s = A[i][j] / (sqrt(A[j][j] * A[j][j] + A[i][j] * A[i][j]));//Коэффициент s
//				tMatrix T1 = nullptr;
//				CreatingArray(T1, N);//Cоздаем промежуточную матрицу матрицу T1 в виде даиагональной единичной
//				for (int k = 0; k < N; k++)
//				{
//					for (int l = 0; l < N; l++)
//					{
//						if (k == l) T1[k][l] = 1.;
//						else T1[k][l] = 0.;
//						/*cout << "Поехали7" << endl;*/
//					}
//				}//Изначально матрица T имеет вид единичной диагональной
//				/*PrintingArray(T1, n);*/
//
//				T1[j][j] = c;
//				T1[j][i] = s;
//				T1[i][j] = -s;
//				T1[i][i] = c;
//				/*cout << "Поехали8" << endl;
//				PrintingArray(T1, n);*/
//
//				//A = Multiplication(T1, A, n);//домножаем А на T1 слева
//				//for (int k = 0; k < n; k++)//домножаем А на T1 слева
//				//for (int k = j; k < n; k++)
//				for (int k = 0; k < N; k++)
//				{
//					datatype IK = A[i][k];
//					datatype JK = A[j][k];
//					A[j][k] = T1[j][j] * JK + T1[j][i] * IK;
//					A[i][k] = T1[i][j] * JK + T1[i][i] * IK;
//				}
//				//B = MVMultiplication(T1, B, n);//домножаем В на T1 слева
//				//T = Multiplication(T1, T, n);//домножаем Т на T1 слева
//				//for (int k = j; k < n; k++)//домножаем T на T1 слева
//				for (int k = 0; k < N; k++)
//				{
//					datatype IK = T[i][k];
//					datatype JK = T[j][k];
//					T[j][k] = T1[j][j] * JK + T1[j][i] * IK;
//					T[i][k] = T1[i][j] * JK + T1[i][i] * IK;
//				}
//				//PrintingArray(A, n);
//				DeletingArray(T1, N);//удаляем промежуточную Т1
//			}
//		}
//	}
//	return T;
//}





bool TrueQ(tMatrix &A)
{
	for (size_t j = 0; j < N - 1; j++)
	{
		for (size_t i = j + 1; i < N; i++)
		{
			if (abs(A[i][j]) > eps)
			{
				return 1;
			}
		}
	}
	return 0;
}

tMatrix One(tMatrix &A)
{

	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			if (i == j)
			{
				A[i][i] = 1;
			}
			else
			{
				A[i][j] = 0;
			}
		}
	}
	return A;
}

tMatrix One(int n)
{
	tMatrix A=nullptr;
	CreatingArray(A,n);
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			if (i == j)
			{
				A[i][i] = 1;
			}
			else
			{
				A[i][j] = 0;
			}
		}
	}
	return A;
}

//tMatrix QR_T(tMatrix &mass1)
//{
//	//+++++++++++++++++++++ATTENTION!!!++++++++++++++++++++++++++++++++++
//	// СМОТРИ СТР. 40!
//	//tMatrix T = nullptr;
//	tMatrix A = MatrCopy(mass1, N);
//	//CreatingArray(T, N);//Матрица Т
//	bool flagIter = 1;
//	int smesh = -1;//порядок смещения(стр.40)
//	while (flagIter)
//	{
//		cout << "hui1" << endl;
//		PrintingArray(A);
//		smesh=smesh+1;
//		if (smesh >= 4) return A;
//		cout << smesh << endl;
//		datatype sigma = A[N - 1 - smesh][N - 1 - smesh];//задаем смещение
//		cout << sigma << endl;
//		for (int opa = 0; opa < N; opa++)
//		{
//			A[opa][opa] = A[opa][opa] - sigma;
//		}
//		cout << "hui2" << endl;
//		PrintingArray(A);
//		bool flag = true;
//		while (flag)
//		{
//		for (int k = 0; k < N; k++) //A:=Q^(-1)AQ
//		{
//			for (int l = k + 1; l < N; l++)
//			{
//				if (abs(A[l][k]) < eps)
//				{
//					A[l][k] = 0;
//				}
//				else
//				{
//
//					datatype alpha = A[k][k] / (sqrt(A[k][k] * A[k][k] + A[l][k] * A[l][k])); //Коэффициент c
//					datatype betta = A[l][k] / (sqrt(A[k][k] * A[k][k] + A[l][k] * A[l][k]));//Коэффициент s
//
//					for (int p = 0; p < N; p++)
//					{
//						datatype IK = A[l][p];
//						datatype JK = A[k][p];
//						A[k][p] = alpha * JK + betta * IK;
//						A[l][p] = alpha * IK - (betta)* JK;
//					}
//					for (int p = 0; p < N; p++)
//					{
//						datatype IK = A[p][l];
//						datatype JK = A[p][k];
//						A[p][k] = alpha * JK - betta * IK;
//						A[p][l] = betta * JK + alpha * IK;
//					}
//
//				}
//			}
//
//		}
//		}
//
//		cout << "hui3" << endl;
//		PrintingArray(A);
//		for (int opa = 1; opa < N; opa++)//обратное смещение
//		{
//			A[opa][opa] = A[opa][opa] + sigma;
//		}
//		flagIter = TrueQ(A);
//		cout << "hui4" << endl;
//		PrintingArray(A);
//	}
//
//	return A;
//}


tMatrix QR_T(tMatrix &mass1)
{
	//+++++++++++++++++++++ATTENTION!!!++++++++++++++++++++++++++++++++++
	// СМОТРИ СТР. 40! ПРЕЖДЕ ЧЕМ ТУПО БЕЗДУМНО ИСПРАВЛЯТЬ ПОСМОТРИ И ПОЙМИ ЧТО ВСЕ СУКА ПРАВИЛЬНО
	//tMatrix T = nullptr;
	tMatrix A = MatrCopy(mass1, N);
	int s4et4ik = 0;
	int s4et4ik2 = 0;
	for (int smesh = 0; smesh < N - 1; smesh++)
	{
		bool flag = 1;
		/*cout << "hui1" << endl;
		PrintingArray(A);*/

		while (flag)//A:=Q^(-1)AQ
		{
			datatype sigma = A[N - 1 - smesh][N - 1 - smesh];//задаем смещение
			cout << sigma << endl;
			for (int opa = 0; opa < N; opa++)
			{
				A[opa][opa] = A[opa][opa] - sigma;
			}

			for (int k = 0; k < N - smesh; k++)
			{
				for (int l = k + 1; l < N - smesh; l++)
				{
					if (abs(A[l][k]) < eps)
					{
						A[l][k] = 0;
					}
					else
					{
						/*cout << "hui3" << endl;
						PrintingArray(A);*/
						s4et4ik2++;
						datatype alpha = A[k][k] / (sqrt(A[k][k] * A[k][k] + A[l][k] * A[l][k])); //Коэффициент c
						datatype betta = A[l][k] / (sqrt(A[k][k] * A[k][k] + A[l][k] * A[l][k]));//Коэффициент s

						for (int p = 0; p < N; p++)
						{
							datatype IK = A[l][p];
							datatype JK = A[k][p];
							A[k][p] = alpha * JK + betta * IK;
							A[l][p] = -betta * JK + alpha * IK;
						}
						for (int p = 0; p < N; p++)
						{
							datatype IK = A[p][l];
							datatype JK = A[p][k];
							A[p][k] = alpha * JK + betta * IK;
							A[p][l] = -betta * JK + alpha * IK;
						}
					}
				}

			} // В итоге Q^(-1)*A*Q
			/*cout << "hui4" << endl;
			PrintingArray(A);*/

			for (int opa = 0; opa < N; opa++)//смещаем назад!
			{
				A[opa][opa] = A[opa][opa] + sigma;
			}

			s4et4ik++;

			datatype suma = 0;
			for (int i = 0; i < N - smesh - 1; i++)// проверка на то что внедиагональные элементы в последней строке нашего минора близки к 0
			{
				suma = suma + abs(A[N - 1 - smesh][i]);
			}
			if (suma < eps7)
			{
				flag = false;
			}
			else flag = true;


		}

	}
	cout << "количество иераций=" << s4et4ik << endl;
	cout << "количество опреаций QR=" << s4et4ik2 << endl;
	return A;
}

tMatrix EugineVectors(tMatrix &mass1, tVector &lambda)
{
	tMatrix A, EVectors = nullptr;
	tVector X0 = nullptr;
	CreatingVector(X0, N);
	CreatingArray(EVectors, N);
	A = MatrCopy(mass1, N);

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			EVectors[i][j] = 0;
		}
	}


	for (size_t i = 0; i < N; i++)//для каждого собственного значения
	{

		tVector X = nullptr;
		tMatrix E = nullptr;
		CreatingArray(E, N);
		E = One(E);
		//E = NumberMatrixMultiplication(lambda[i],E);
		for (int j = 0; j < N; j++)
		{
			E[j][j] = E[j][j] * lambda[i];
		}
		//tMatrix C = Difference(A,E,N);
		tMatrix C = MatrCopy(A);
		cout << "Matrix C:" << endl;
		PrintingArray(C);
		for (int j = 0; j < N; j++)
		{
			C[j][j] = C[j][j] - E[j][j];
		}
		/////////////////////////////
		cout << "Matrix C:" << endl;
		PrintingArray(C);
		cout << endl;



		for (size_t j = 0; j < N; j++)
		{
			X0[j] = 0;
		}
		X0[0] = 1; //инициализировали вектор Х0

		/////////////////////////////
		PrintingVector(X0);
		/////////////////////////////
		cin.get();
		//double NormVectDX = 1;
		bool flaq = true;
		while (flaq)//может выдавать ошибку из-за медленной сходимоcти
		{

			cout << "1(X0): " << endl;
			PrintingVector(X0);
			cout << "lambda = " << lambda[i] << endl;
			double norm;
			//X = MVMultiplication(C, X0, N);//получаем вектор y из методы
			X = QRmethod(C, X0);//получаем вектор y из методы
			//X = Gauss(C, X0);
			/////////////////////////////
			PrintingVector(X);
			cout << endl;
			/////////////////////////////

			norm = NormOfVector(X, N);
			cout << "2(X): " << endl;
			X = NumberVectorMultiplication(1 / norm, X);//нормируем его

			/////////////////////////////
			PrintingVector(X);
			NormOfVector(X, N);
			cout << endl;
			/////////////////////////////
			cin.get();
			cout << "3(DX): " << endl;
			/*tVector DX = Difference(X, X0, N);
			PrintingVector(DX);
			cout << endl;
			NormVectDX = NormOfVector(DX, N);*/
			/*NormVectDX = abs(abs(Scalar(X0, X)) - 1);
			cout << NormVectDX << endl;*/
			tVector AX = MVMultiplication(A, X, N);
			tVector LX = NumberVectorMultiplication(lambda[i], X);
			tVector HUY = Difference(AX, LX, N);
			if (NormOfVector(HUY, N) < eps)
			{
				flaq = false;
			}
			X0 = VectCopy(X, N);
			cout << "So1: " << endl;
			PrintingVector(X0);
		}

		cout << "So: " << endl;
		PrintingVector(X0);
		for (int j = 0; j < N; j++)
		{
			EVectors[i][j] = X0[j];
		}

		PrintingArray(EVectors);
		/*cout << lambda[i] << endl;
		PrintingVector(X0,N);*/
		DeletingVector(X);
		cin.get();
	}

	PrintingArray(EVectors);

	for (size_t i = 0; i < N; i++)
	{
		cout << "lambda = " << lambda[i] << "  Vect: ";
		for (size_t j = 0; j < N; j++)
		{
			cout << EVectors[i][j] << "  ";
		}
		cout << endl;
	}

	DeletingArray(A, N);
	DeletingVector(X0);
	return EVectors;
}
























//ЧТОБЫ НЕ ПУТАТЬ





























tMatrix Reley(tMatrix &mass1)
{
	tMatrix A, EVectors = nullptr;
	tVector X0, lambda = nullptr;
	CreatingVector(X0, N);
	CreatingVector(lambda, N);
	CreatingArray(EVectors, N);
	A = MatrCopy(mass1, N);

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			EVectors[i][j] = 0;
		}
	}


	for (size_t i = 0; i < N; i++)//для каждого собственного значения
	{
		tVector X = nullptr;
		tMatrix E = nullptr;

		//for (size_t j = 0; j < N; j++)
		//{
		//	X0[j] = i-j;
		//}
		//X0[i] = -1 -i; //инициализировали вектор Х0

		if (i == 0)
		{
			X0[0] = -0.8;
			X0[1] = 0;
			X0[2] = -0.2;
			X0[3] = -0.4;
		}

		if (i == 1)
		{
			X0[0] = 0;
			X0[1] = 0.7;
			X0[2] = -0.6;
			X0[3] = 0.35;
		}

		if (i == 2)
		{
			X0[0] = 0.5;
			X0[1] = 0;
			X0[2] = -0.43;
			X0[3] = -0.75;
		}

		if (i == 3)
		{
			X0[0] = 0;
			X0[1] = 0.7;
			X0[2] = 0.6;
			X0[3] = -0.35;
		}

		//		tVector X = nullptr;
		bool flaq = true;

		while (flaq)//может выдавать ошибку из-за медленной сходимоcти
		{
			tVector AX0 = MVMultiplication(A, X0, N);

			lambda[i] = Scalar(AX0, X0);
			cout << "lambda[" << i << "]=" << lambda[i] << endl;
			CreatingArray(E, N);

			E = One(E);

			for (int j = 0; j < N; j++)
			{
				E[j][j] = E[j][j] * lambda[i];
			}
			cout << "Matrix E:" << endl;
			PrintingArray(E);

			tMatrix C = MatrCopy(A);
			cout << "Matrix C:" << endl;
			PrintingArray(C);
			for (int j = 0; j < N; j++)
			{
				C[j][j] = C[j][j] - E[j][j];
			}

			cout << "1(X0): " << endl;
			PrintingVector(X0);
			cout << "lambda = " << lambda[i] << endl;
			double norm;

			X = QRmethod(C, X0);

			PrintingVector(X);
			cout << endl;


			norm = NormOfVector(X, N);
			cout << "2(X): " << endl;
			X = NumberVectorMultiplication(1 / norm, X);//нормируем его


			PrintingVector(X);
			NormOfVector(X, N);
			cout << endl;

			cin.get();
			cout << "3(DX): " << endl;
			tVector AX = MVMultiplication(A, X, N);
			tVector LX = NumberVectorMultiplication(lambda[i], X);
			tVector HUY = Difference(AX, LX, N);
			if (NormOfVector(HUY, N) < eps)
			{
				flaq = false;
			}
			X0 = VectCopy(X, N);
			cout << "So1: " << endl;
			PrintingVector(X0);
		}

		cout << "So: " << endl;
		PrintingVector(X0);
		for (int j = 0; j < N; j++)
		{
			EVectors[i][j] = X0[j];
		}

		PrintingArray(EVectors);
		/*cout << lambda[i] << endl;
		PrintingVector(X0,N);*/
		DeletingVector(X);
		cin.get();
	}

	PrintingArray(EVectors);

	for (size_t i = 0; i < N; i++)
	{
		cout << "lambda = " << lambda[i] << "  Vect: ";
		for (size_t j = 0; j < N; j++)
		{
			cout << EVectors[i][j] << "  ";
		}
		cout << endl;
	}

	DeletingArray(A, N);
	DeletingVector(X0);
	return EVectors;
}


//tMatrix Reley(tMatrix &mass1)
//{
//	tMatrix A, EVectors = nullptr;
//	tVector X0, lambda = nullptr;
//	CreatingVector(X0, N);
//	CreatingVector(lambda, N);
//	CreatingArray(EVectors, N);
//	A = MatrCopy(mass1, N);
//
//	for (int i = 0; i < N; i++)
//	{
//		for (int j = 0; j < N; j++)
//		{
//			EVectors[i][j] = 0;
//		}
//	}
//
//
//	for (size_t i = 0; i < N; i++)//для каждого собственного значения
//	{
//		tVector X = nullptr;
//		tMatrix E = nullptr;
//
//		for (size_t j = 0; j < N; j++)
//		{
//			X0[j] = 0;
//		}
//		X0[1] = 1; //инициализировали вектор Х0
//
//
//
//		//		tVector X = nullptr;
//		bool flaq = true;
//
//		while (flaq)//может выдавать ошибку из-за медленной сходимоcти
//		{
//			tVector AX0 = MVMultiplication(A, X0, N);
//
//			lambda[i] = Scalar(AX0, X0);
//			cout << "lambda[" << i << "]=" << lambda[i] << endl;
//			CreatingArray(E, N);
//
//			E = One(E);
//
//			for (int j = 0; j < N; j++)
//			{
//				E[j][j] = E[j][j] * lambda[i];
//			}
//			cout << "Matrix E:" << endl;
//			PrintingArray(E);
//
//			tMatrix C = MatrCopy(A);
//			cout << "Matrix C:" << endl;
//			PrintingArray(C);
//			for (int j = 0; j < N; j++)
//			{
//				C[j][j] = C[j][j] - E[j][j];
//			}
//
//			cout << "1(X0): " << endl;
//			PrintingVector(X0);
//			cout << "lambda = " << lambda[i] << endl;
//			double norm;
//
//			X = QRmethod(C, X0);
//
//			PrintingVector(X);
//			cout << endl;
//
//
//			norm = NormOfVector(X, N);
//			cout << "2(X): " << endl;
//			X = NumberVectorMultiplication(1 / norm, X);//нормируем его
//
//
//			PrintingVector(X);
//			NormOfVector(X, N);
//			cout << endl;
//
//			cin.get();
//			cout << "3(DX): " << endl;
//			tVector AX = MVMultiplication(A, X, N);
//			tVector LX = NumberVectorMultiplication(lambda[i], X);
//			tVector HUY = Difference(AX, LX, N);
//			if (NormOfVector(HUY, N) < eps)
//			{
//				flaq = false;
//			}
//			X0 = VectCopy(X, N);
//			cout << "So1: " << endl;
//			PrintingVector(X0);
//		}
//
//		cout << "So: " << endl;
//		PrintingVector(X0);
//		for (int j = 0; j < N; j++)
//		{
//			EVectors[i][j] = X0[j];
//		}
//
//		PrintingArray(EVectors);
//		/*cout << lambda[i] << endl;
//		PrintingVector(X0,N);*/
//		DeletingVector(X);
//		cin.get();
//	}
//
//	PrintingArray(EVectors);
//
//	for (size_t i = 0; i < N; i++)
//	{
//		cout << "lambda = " << lambda[i] << "  Vect: ";
//		for (size_t j = 0; j < N; j++)
//		{
//			cout << EVectors[i][j] << "  ";
//		}
//		cout << endl;
//	}
//
//	DeletingArray(A, N);
//	DeletingVector(X0);
//	return EVectors;
//}



tVector MainDiag(tMatrix &A, tVector &diag1)
{
	tVector diag = nullptr;
	CreatingVector(diag, N);
	for (size_t i = 0; i < N; i++)
	{
		diag[i] = A[i][i];
	}
	return diag;
}


datatype f1(datatype x)
{
	return (x - 0.1)*(x - 0.22)*(x - 0.55)*(x - 0.7)*(x - 0.75);
}

datatype f2(datatype x)
{
	return sqrt(x + 1) - 1;
}

datatype f3(datatype x)
{
	return 35 * x*x*x - 67 * x*x - 3 * x + 3;
}

datatype f4(datatype x)
{
	return x * x;
}

datatype f5(datatype x)
{
	return -q * (x - 1) / 2;
}

datatype f_w(datatype x)
{
	return sin(M_PI*x);
}

datatype f_w2(datatype x)
{
	return x*(1-x);
}

datatype g_w(datatype x)
{
	return 0;
}

datatype phi_w(datatype x)
{
	return 0;
}

datatype psi_w(datatype x)
{
	return 0;
}

datatype Kxpr(datatype x)
{
	if ( 0 <= x && x <= 0.25)
	{
		return 10;
	}
	if ( 0.25 < x && x < 0.5)
	{
		return 1./ (-0.4 * x + 0.2 + 6 * x - 1.5);
	}
	return 0.66666666666666666;
}

datatype Kx(datatype x)
{
	if ( 0 <= x && x <= 0.25)
	{
		return 0.1;
	}
	if ( 0.25 < x && x < 0.5)
	{
		return (-0.4 * x + 0.2 + 6 * x - 1.5);
	}
	return 1.5;
}

datatype P_Zero(datatype x)
{
	return 0;
}

datatype P_Heat(datatype x)
{
	return (5 + 0.1*x*x*x*x);
}

datatype K10(datatype x)
{
	return 10.;
}

datatype K1(datatype x)
{
	return 1.;
}

datatype P1(datatype x)
{
	if ( 0 <= x && x < 0.5)
	{
		return 10;
	}
//	if ( 0.5 <= x )
//	{
//		return 0;
//	}
	return 0;
}

datatype P2(datatype x)
{
	if ( 0 <= x && x < 0.5)
	{
		return 20 * x;
	}
//	if ( 0.5 <= x )
//	{
//		return 0;
//	}
	return 0;
}

datatype P3(datatype x)
{
	if ( 0 <= x && x < 0.5)
	{
		return 20 * (0.5 - x);
	}
//	if ( 0.5 <= x )
//	{
//		return 0;
//	}
	return 0;
}

datatype P4(datatype x)
{
	if ( 0 <= x && x <= 0.25 )
	{
		return 20 * x;
	}
	if ( 0.25 < x && x < 0.5 )
	{
		return 20 * (0.5 - x);
	}
//	if ( 0.5 <= x)
//	{
//		return 0;
//	}
	return 0;
}


size_t PrintMENU()
{
	size_t s = 0;
	cout << endl;
	cout << "Choose your destiny: " << endl;
	cout << "1 - 1 example " << endl;
	cout << "2 - 2 example " << endl;
	cout << "3 - 3 example  " << endl;
	cout << "4 - Pendeulum( sym )" << endl;
	/*cout << "5 - 5-й пример ( система 4 решения )" << endl;*/
	//cout << "6 - y=x*x" << endl;
	cout << "5 - Big Example" << endl;
	cout << "6 - Analytical Solve" << endl;
	cout << "0 - Exit" << endl;
	cout << endl;
	cout << "Your choise: ";
	cin >> s;
	cout << endl;
	return s;
}

void clear()
{
	for (size_t i = 0; i < 6; i++)
	{
		cout << endl;
	}
}

tVector FindRoot_WilO4KA_CQ5(funct &f)
{
	cout << "МЕТОД БИСЕКЦИИ КОНТРОЛЬНЫЙ ВОПРОС 5" << endl;
	size_t n = 0;
	datatype fi, ff; //f initial, f final
	tMatrix AB = nullptr;
	datatype a = 1 / (1 + q);
	datatype b = 4;
	while (a < b)
	{
		fi = a;
		ff = a + h;
		if (f(fi)*f(ff) < 0)
		{
			n++;
			if (n == 1)
			{
				AB = new tVector[1];
				AB[0] = new datatype[2];
				AB[0][0] = fi;
				AB[0][1] = ff;
			}
			else
			{
				tMatrix boofAB = AB;
				AB = nullptr;
				AB = new tVector[n];
				for (size_t i = 0; i < n; i++)
				{
					AB[i] = new datatype[2];
				}
				for (size_t i = 0; i < n - 1; i++)
				{
					for (size_t j = 0; j < 2; j++)
					{
						AB[i][j] = boofAB[i][j];
					}
					delete[] boofAB[i];  //ОПАСНО ИМХО
				}
				delete boofAB;
				AB[n - 1][0] = fi;
				AB[n - 1][1] = ff;
			}
		}
		a += h;
	}
	if (AB == nullptr)
	{
		cout << "Корней нет!" << endl;
		return nullptr;
	}
	else
	{
		cout << "Всего корней: " << n << endl;
	}

	tVector Root = new datatype[n];
	size_t count_iter = 0;
	for (size_t i = 0; i < n; i++)
	{
		a = AB[i][0];
		b = AB[i][1];

		//cout <<"Begin: "<< a << "  " << b << endl;
		//cout << f(a)*f((b-a)/2) << endl;
		//cin.get();
		while ((b - a) > 2 * eps7)
		{
			datatype c = (b + a) / 2;
			if (abs(f(c)) < eps7) cout << "Срабатывание метода f(c) < Eps" << endl;
			if ((f(a)*f(c)) < 0)
			{
				b = c;
			}
			else
			{
				a = c;
			}
			//cout << a << "  " << b << endl;
			//cin.get();
			count_iter++;
		}
		Root[i] = (a + b) / 2;
		if (abs(Root[i]) < 1.e-7)
		{
			Root[i] = 0;
		}
	}
	for (size_t i = 0; i < n; i++)
	{
		cout << i + 1 << ")  " << Root[i] << endl;
	}

	cout << "Метод бисекции сошёлся за " << count_iter << " итераций" << endl;
	cout << endl;
	return Root;
}

tVector FindRoot_WilO4KA(funct &f, datatype a, datatype b)
{
	cout << "МЕТОД БИСЕКЦИИ" << endl;
	size_t n = 0;
	datatype fi, ff; //f initial, f final
	tMatrix AB = nullptr;
	while (a < b)
	{
		fi = a;
		ff = a + h;
		if (f(fi)*f(ff) < 0)
		{
			n++;
			if (n == 1)
			{
				AB = new tVector[1];
				AB[0] = new datatype[2];
				AB[0][0] = fi;
				AB[0][1] = ff;
			}
			else
			{
				tMatrix boofAB = AB;
				AB = nullptr;
				AB = new tVector[n];
				for (size_t i = 0; i < n; i++)
				{
					AB[i] = new datatype[2];
				}
				for (size_t i = 0; i < n - 1; i++)
				{
					for (size_t j = 0; j < 2; j++)
					{
						AB[i][j] = boofAB[i][j];
					}
					delete[] boofAB[i];  //ОПАСНО ИМХО
				}
				delete boofAB;
				AB[n - 1][0] = fi;
				AB[n - 1][1] = ff;
			}
		}
		a += h;
	}
	if (AB == nullptr)
	{
		cout << "Корней нет!" << endl;
		return nullptr;
	}
	else
	{
		cout << "Всего корней: " << n << endl;
	}
	////////////////////////////////////////
	cout << endl;
	cout << "Отрезки локализации: " << endl;
	for (size_t i = 0; i < n; i++)
	{
		cout << "[ " << AB[i][0] << " , " << AB[i][1] << " ]" << endl;
	}
	cout << endl;
	////////////////////////////////////////
	tVector Root = new datatype[n];
	size_t count_iter = 0;
	for (size_t i = 0; i < n; i++)
	{
		a = AB[i][0];
		b = AB[i][1];

		//cout <<"Begin: "<< a << "  " << b << endl;
		//cout << f(a)*f((b-a)/2) << endl;
		//cin.get();
		while ((b - a) > 2 * eps5)
		{
			datatype c = (b + a) / 2;
			cout << "Очередное приближение х = " << c << endl;
			if ((f(a)*f(c)) < 0)
			{
				b = c;
			}
			else
			{
				a = c;
			}
			//cout << a << "  " << b << endl;
			//cin.get();
			count_iter++;
		}
		Root[i] = (a + b) / 2;
		if (abs(Root[i]) < eps7)
		{
			Root[i] = 0;
		}
	}
	for (size_t i = 0; i < n; i++)
	{
		cout << i + 1 << ")  " << Root[i] << endl;
	}

	cout << "Метод бисекции сошёлся за " << count_iter << " итераций" << endl;
	cout << endl;
	return Root;
}


tVector FindRoot_Neuton_Num(funct &f, datatype a, datatype b)
{
	cout << "МЕТОД НЬЮТОНА ЧИСЛЕННЫЙ" << endl;
	size_t count_iter = 0;
	size_t n = 0;
	datatype fi, ff; //f initial, f final
	tMatrix AB = nullptr;
	while (a < b)
	{
		fi = a;
		ff = a + h;
		if (f(fi)*f(ff) < 0)
		{
			n++;
			if (n == 1)
			{
				AB = new tVector[1];
				AB[0] = new datatype[2];
				AB[0][0] = fi;
				AB[0][1] = ff;
			}
			else
			{
				tMatrix boofAB = AB;
				AB = nullptr;
				AB = new tVector[n];
				for (size_t i = 0; i < n; i++)
				{
					AB[i] = new datatype[2];
				}
				for (size_t i = 0; i < n - 1; i++)
				{
					for (size_t j = 0; j < 2; j++)
					{
						AB[i][j] = boofAB[i][j];
					}
					delete[] boofAB[i];  //ОПАСНО ИМХО
				}
				delete boofAB;
				AB[n - 1][0] = fi;
				AB[n - 1][1] = ff;
			}
		}
		a += h;
	}
	if (AB == nullptr)
	{
		cout << "Корней нет!" << endl;
		return nullptr;
	}
	else
	{
		cout << "Всего корней: " << n << endl;
	}
	////////////////////////////////////////
	cout << endl;
	cout << "Отрезки локализации: " << endl;
	for (size_t i = 0; i < n; i++)
	{
		cout << "[ " << AB[i][0] << " , " << AB[i][1] << " ]" << endl;
	}
	cout << endl;
	////////////////////////////////////////

	tVector Root = new datatype[n];
	for (size_t i = 0; i < n; i++)
	{
		datatype x0, xk, fprime;
		datatype difference = 1;
		x0 = AB[i][0];
		//xk = x0 + eps;
		while (difference > eps5)
		{
			fprime = (f(x0 + eps) - f(x0)) / eps;
			xk = x0 - (f(x0) / fprime);
			difference = abs(x0 - xk);
			x0 = xk;
			if (xk<AB[i][0] || xk>AB[i][1])
			{
				x0 = AB[i][1];
			}
			cout << "Очередное приближение х = " << x0 << endl;
			count_iter++;
		}
		Root[i] = xk;
		if (abs(Root[i]) < eps7)
		{
			Root[i] = 0;
		}
		cout << "Корень найден" << endl;
	}
	for (size_t i = 0; i < n; i++)
	{
		cout << i + 1 << ")  " << Root[i] << endl;
	}
	cout << "Численный метод Ньютона сошёлся за " << count_iter << " итераций" << endl;
	cout << endl;
	delete[] AB;
	return Root;
}

tVector FindRoot_Neuton_Anal(funct &f, datatype a, datatype b)
{
	cout << "МЕТОД НЬЮТОНА АНАЛИТИЧЕСКИЙ" << endl;
	size_t count_iter = 0;
	size_t n = 0;
	datatype fi, ff; //f initial, f final
	tMatrix AB = nullptr;
	while (a < b)
	{
		fi = a;
		ff = a + h;
		if (f(fi)*f(ff) < 0)
		{
			n++;
			if (n == 1)
			{
				AB = new tVector[1];
				AB[0] = new datatype[2];
				AB[0][0] = fi;
				AB[0][1] = ff;
			}
			else
			{
				tMatrix boofAB = AB;
				AB = nullptr;
				AB = new tVector[n];
				for (size_t i = 0; i < n; i++)
				{
					AB[i] = new datatype[2];
				}
				for (size_t i = 0; i < n - 1; i++)
				{
					for (size_t j = 0; j < 2; j++)
					{
						AB[i][j] = boofAB[i][j];
					}
					delete[] boofAB[i];  //ОПАСНО ИМХО
				}
				delete boofAB;
				AB[n - 1][0] = fi;
				AB[n - 1][1] = ff;
			}
		}
		a += h;
	}
	if (AB == nullptr)
	{
		cout << "Корней нет!" << endl;
		return nullptr;
	}
	else
	{
		cout << "Всего корней: " << n << endl;
	}
	//////////////////////////////////////////
	//cout << endl;
	//cout << "Отрезки локализации: " << endl;
	//for (size_t i = 0; i < n; i++)
	//{
	//	cout << "[ " << AB[i][0] << " , " << AB[i][1] << " ]" << endl;
	//}
	//cout << endl;
	//////////////////////////////////////////

	tVector Root = new datatype[n];
	for (size_t i = 0; i < n; i++)
	{
		datatype x0, xk, fprime;
		datatype difference = 1;
		x0 = AB[i][0];
		xk = x0 + eps;
		while (difference > eps5)
		{
			fprime = 5 * x0*x0*x0*x0 - 232 * x0*x0*x0 / 25 + 11907 * x0*x0 / 2000 - 15119 * x0 / 10000 + 24299 / 200000;
			xk = x0 - (f(x0) / fprime);
			difference = abs(x0 - xk);
			x0 = xk;
			if (xk<AB[i][0] || xk>AB[i][1])
			{
				x0 = AB[i][1];
			}
			cout << "Очередное приближение х = " << x0 << endl;
			count_iter++;
		}
		Root[i] = xk;
		if (abs(Root[i]) < eps5)
		{
			Root[i] = 0;
		}
		cout << "Корень найден" << endl;
	}
	for (size_t i = 0; i < n; i++)
	{
		cout << i + 1 << ")  " << Root[i] << endl;
	}
	cout << "Аналитический метод Ньютона сошёлся за " << count_iter << " итераций" << endl;
	cout << endl;
	delete[] AB;
	return Root;
}


tVector FindRoot_Neuton(funct &f, datatype a, datatype b)
{
	cout << "МЕТОД НЬЮТОНА" << endl;
	size_t count_iter = 0;
	size_t n = 0;
	datatype fi, ff; //f initial, f final
	tMatrix AB = nullptr;
	while (a < b)
	{
		fi = a;
		ff = a + h;
		if (f(fi)*f(ff) < 0)
		{
			n++;
			if (n == 1)
			{
				AB = new tVector[1];
				AB[0] = new datatype[2];
				AB[0][0] = fi;
				AB[0][1] = ff;
			}
			else
			{
				tMatrix boofAB = AB;
				AB = nullptr;
				AB = new tVector[n];
				for (size_t i = 0; i < n; i++)
				{
					AB[i] = new datatype[2];
				}
				for (size_t i = 0; i < n - 1; i++)
				{
					for (size_t j = 0; j < 2; j++)
					{
						AB[i][j] = boofAB[i][j];
					}
					delete[] boofAB[i];  //ОПАСНО ИМХО
				}
				delete boofAB;
				AB[n - 1][0] = fi;
				AB[n - 1][1] = ff;
			}
		}
		a += h;
	}
	if (AB == nullptr)
	{
		cout << "Корней нет!" << endl;
		return nullptr;
	}
	else
	{
		cout << "Всего корней: " << n << endl;
	}
	//////////////////////////////////////////
	//cout << endl;
	//cout << "Отрезки локализации: " << endl;
	//for (size_t i = 0; i < n; i++)
	//{
	//	cout << "[ " << AB[i][0] << " , " << AB[i][1] << " ]" << endl;
	//}
	//cout << endl;
	//////////////////////////////////////////

	tVector Root = new datatype[n];
	for (size_t i = 0; i < n; i++)
	{
		datatype x0, xk, fprime;
		datatype difference = 1;
		x0 = AB[i][0];
		xk = x0 + eps;
		while (difference > eps5)
		{
			fprime = (f(x0 + eps) - f(x0)) / eps;
			xk = x0 - (f(x0) / fprime);
			difference = abs(x0 - xk);
			x0 = xk;
			//cout << "Очередное приближение х = " << x0 << endl;
			count_iter++;
		}
		Root[i] = xk;
		if (abs(Root[i]) < eps5)
		{
			Root[i] = 0;
		}
	}
	for (size_t i = 0; i < n; i++)
	{
		cout << i + 1 << ")  " << Root[i] << endl;
	}
	cout << "Метод Ньютона сошёлся за " << count_iter << " итераций" << endl;
	cout << endl;
	delete[] AB;
	return Root;
}

tVector dsys1(tVector &X0, int n)
{
	tVector sol = new datatype[n];
	sol[0] = (2 * X0[0] + X0[1] * X0[1] - 1);
	sol[1] = (6 * X0[0] - X0[1] * X0[1] + 1);
	return sol;
}

tVector dsys2(tVector &X0, int n)
{
	tVector sol = new datatype[n];
	sol[0] = (1 - X0[0] * X0[0] - X0[1] * X0[1]);
	sol[1] = (2 * X0[0]);
	return sol;
}

tVector dsys3(tVector &X0, int n)
{
	tVector sol = new datatype[n];
	sol[0] = 10 * (X0[1] - X0[0]);
	sol[1] = X0[0] * (28 - X0[2]) - X0[1];
	sol[2] = X0[0] * X0[1] - (8 / 3)*X0[2];
	return sol;
}

tVector dsys4(tVector &X0, int n)
{
	tVector sol = new datatype[n];
	sol[0] = -X0[0];
	sol[1] = -2 * X0[1];
	return sol;
}

tVector dsysAnal(tVector &X0, int n)
{
	tVector sol = new datatype[n];
	sol[0] = 2 * X0[0] - 5 * X0[1] + 3;
	sol[1] = 5 * X0[0] - 6 * X0[1] + 1;
	return sol;
}

tVector dsysAnalSolve(datatype t)
{
	tVector sol = new datatype[2];
	sol[0] = 5 * (exp(-2 * t))*cos(3 * t) + 1;
	sol[1] = (exp(-2 * t))*(4 * (cos(3 * t)) + 3 * (sin(3 * t))) + 1;
	return sol;
}

tVector dsysRozenbrok1(tVector &X0, int n)
{
	tVector sol = new datatype[n];
	sol[0] = -50*(X0[0] - cos(X0[1]));
	sol[1] = 1;
	return sol;
}

tVector Oscill(tVector &X0, int n)
{
	tVector sol = new datatype[n];
	sol[0] = X0[1];
	sol[1] = -(200 / 3) * X0[0];
	return sol;
}

tVector fsys2(datatype x, datatype y)
{
	tVector sol = new datatype[2];
	sol[0] = x * x + y * y + x + y - 8;
	sol[1] = x * x + y * y + x * y - 7;
	return sol;
}

tVector fsys3(datatype x, datatype y)
{
	tVector sol = new datatype[2];
	sol[0] = 3 * x + 5 * y - 1;
	sol[1] = x + 3 * y + 2;
	return sol;
}

tMatrix FindRoot_Sys_Num(functsys &fsys, datatype a, datatype b)
{
	char FileName[] = "Data.txt";
	ofstream OUT(FileName);
	OUT << "{";
	size_t sumiter = 0;
	cout << "НЬЮТОН ЧИСЛЕННЫЙ" << endl;
	datatype hh;
	size_t n = 0;   //число решений
	tMatrix Root = nullptr; //массив решений. Координаты записаны по строкам, как х и у соответственно
	datatype x, y = -a;
	cout << "Введите шаг разбиения: ";
	cin >> hh;
	cout << endl << endl;
	datatype alliter = (2 * a) / (hh)-1; //всего узлов на квадрате со стороной a и шагом h
	for (size_t i = 0; i < alliter; i++) //ходим по оси у
	{

		y += hh;
		x = -b;
		for (size_t j = 0; j < alliter; j++) //ходим по оси х по узлам
		{
			bool flag = 0; //изначально 0 - нет совпадений  с решениями в массиве
			bool nullJ = 0; //флаг перехода к следующему узлу из-за вырожденности матрицы Якоби
			x += hh;
			if (abs(x) < eps7) x = 0;
			size_t counter = 0; //счетчик числа итераций метода Ньютона
			datatype normdif = 10; //для нормы вектора между итерациями
			tVector x0 = new datatype[2];
			tVector xk = new datatype[2];
			tVector jac = new datatype[2];
			tVector jacboof = new datatype[2];
			tVector Fxk = new datatype[2];
			tVector dx = new datatype[2];
			tMatrix Jacobi = nullptr;
			CreatingArray(Jacobi, 2);
			x0[0] = x;					//x0 - вектор {x,y}
			x0[1] = y;					//x0 - вектор {x,y}

			while (normdif > eps7 && counter < 31)  //если не сходится за 30 итераций, то выход из цикла
			{									//и переход к следующему узлу
				sumiter++;
				counter++;
				//строим матрицу Якоби в точке x0
				//jac = fsys[(x + eps, y) - fsys(x, y)] / eps;					//сначала производные по х
				jac = fsys(x0[0] + eps, x0[1]);
				jacboof = fsys(x0[0], x0[1]);
				jac = Difference(jac, jacboof, 2);
				jac = NumberVectorMultiplication(1 / eps, jac, 2);
				Jacobi[0][0] = jac[0];
				Jacobi[1][0] = jac[1];
				//jac = fsys[(x, y + eps) - fsys(x, y)] / eps;					//и по у
				jac = fsys(x0[0], x0[1] + eps);
				jacboof = fsys(x0[0], x0[1]);
				jac = Difference(jac, jacboof, 2);
				jac = NumberVectorMultiplication(1 / eps, jac, 2);
				Jacobi[0][1] = jac[0];
				Jacobi[1][1] = jac[1];
				/*cout << "Якоби: " << endl;
				PrintingArray(Jacobi, 2);*/
				//готовая матрица частных производных
				Jacobi = Inverse(Jacobi, 2);
				if (Jacobi == nullptr)
				{
					nullJ = 1; goto newiter;
				}
				Fxk = fsys(x0[0], x0[1]);
				Fxk = MVMultiplication(Jacobi, Fxk, 2);
				xk = Difference(x0, Fxk, 2);
				/*cout << "Вектор xk: ";
				cout << xk[0] << "  " << xk[1] << endl;
				cout << endl;*/
				dx = Difference(x0, xk, 2);
				normdif = NormOfVector(dx, 2);
				//cout << "Норма вектора dx: " << normdif << endl;
				x0 = xk;
			} //выход из цикла, либо если counter = 30, либо итерационный процесс сошёлся
			if (counter < 30)//если процесс сошёлся
			{
				if (n == 0) //если решение первое
				{
					Root = new tVector[1];
					Root[0] = new datatype[2];
					Root[0][0] = xk[0];
					Root[0][1] = xk[1];
					n++;
				}
				else //иначе проверяем сходство с данными в массиве
				{
					for (size_t ii = 0; ii < n; ii++)
					{
						if (abs(Root[ii][0] - xk[0]) < eps  && abs(Root[ii][1] - xk[1]) < eps) //значения х и у отличаются
						{																    	// более, чем на eps
							flag = 1; //есть совпадение
						}
					}
					if (flag == 0)   //если совпадений не было, то добавляем размерность массива, иначе ничего не делаем
					{
						n++;
						tMatrix boofRoot = Root;
						Root = nullptr;
						Root = new tVector[n];
						for (size_t k = 0; k < n; k++)
						{
							Root[k] = new datatype[2];
						}
						for (size_t k = 0; k < n - 1; k++)
						{
							for (size_t jj = 0; jj < 2; jj++)
							{
								Root[k][jj] = boofRoot[k][jj];
							}
							delete[] boofRoot[k];  //ОПАСНО ИМХО
						}
						delete boofRoot;
						Root[n - 1][0] = xk[0];
						Root[n - 1][1] = xk[1];
						/*cout << "Записанный вектор xk в массив: " << endl;
						cout<< Root[n - 1][0] <<"  "<< Root[n - 1][1] <<endl;*/
					}
				}
			}//а если не сошёлся, то переходим к следующему узлу
			//cin.get();
		newiter:
			if (nullJ)//если матрица Якоби вырождена
			{
				OUT << "{" << x << "," << y << "," << 30 << "}";
				if (!(i + 1 >= alliter && j + 1 >= alliter))//если узел НЕ последний
				{
					OUT << ",";
				}
			}
			else
			{
				OUT << "{" << x << "," << y << "," << counter << "}";
				if (!(i + 1 >= alliter && j + 1 >= alliter))//если узел НЕ последний
				{
					OUT << ",";
				}
			}


			if (Jacobi != nullptr)
			{
				DeletingArray(Jacobi, 2);
			}

			delete[] x0;
			delete[] jac;
			delete[] jacboof;
			delete[] Fxk;
			delete[] dx;
			/*cout << xk[0] << "  " << xk[1] << endl;
			delete[] xk;*/
		}

	}
	OUT << "}";
	OUT.close();
	if (Root == nullptr)
	{
		cout << "Решений нет!" << endl;
		return nullptr;
	}
	else
	{
		cout << "Всего решений: " << n << endl;
	}
	for (size_t i = 0; i < n; i++)
	{
		cout << i + 1 << ") x = " << Root[i][0] << " , y = " << Root[i][1] << endl;
	}
	cout << "Всего произведено " << sumiter << " итераций" << endl;
	cout << endl;
	return Root;
}


tMatrix FindRoot_Sys1(functsys &fsys, datatype a, datatype b)
{
	size_t sumiter = 0;
	cout << "НЬЮТОН АНАЛИТИЧЕСКИЙ" << endl;
	datatype hh;
	size_t n = 0;   //число решений
	tMatrix Root = nullptr; //массив решений. Координаты записаны по строкам, как х и у соответственно
	datatype x, y = -a;
	cout << "Введите шаг разбиения: ";
	cin >> hh;
	cout << endl << endl;
	datatype alliter = (2 * a) / (hh)-1; //всего узлов на квадрате со стороной a и шагом 0.1
	for (size_t i = 0; i < alliter; i++) //ходим по оси у
	{

		y += hh;
		x = -10;
		for (size_t j = 0; j < alliter; j++) //ходим по оси х по узлам
		{
			bool flag = 0; //изначально 0 - нет совпадений
			x += hh;
			size_t counter = 0; //счетчик числа итераций метода Ньютона
			datatype normdif = 10; //для нормы вектора между итерациями
			tVector x0 = new datatype[2];
			tVector xk = new datatype[2];
			//tVector jac = new datatype[2];
			//tVector jacboof = new datatype[2];
			tVector Fxk = new datatype[2];
			tVector dx = new datatype[2];
			tMatrix Jacobi = nullptr;
			CreatingArray(Jacobi, 2);
			x0[0] = x;					//x0 - вектор {x,y}
			x0[1] = y;					//x0 - вектор {x,y}

			while (normdif > eps7 && counter < 31)  //если не сходится за 30 итераций, то выход из цикла
			{									//и переход к следующему узлу
				sumiter++;
				counter++;
				//строим матрицу Якоби в точке x0
				//jac = fsys[(x + eps, y) - fsys(x, y)] / eps;					//сначала производные по х
				/*jac = fsys(x0[0] + eps, x0[1]);
				jacboof = fsys(x0[0], x0[1]);
				jac = Difference(jac, jacboof, 2);
				jac = NumberVectorMultiplication(1 / eps, jac, 2);
				Jacobi[0][0] = jac[0];
				Jacobi[1][0] = jac[1];*/
				Jacobi[0][0] = 2 * x0[0];
				Jacobi[1][0] = -2 * x0[1];
				//jac = fsys[(x, y + eps) - fsys(x, y)] / eps;					//и по у
				/*jac = fsys(x0[0], x0[1] + eps);
				jacboof = fsys(x0[0], x0[1]);
				jac = Difference(jac, jacboof, 2);
				jac = NumberVectorMultiplication(1 / eps, jac, 2);
				Jacobi[0][1] = jac[0];
				Jacobi[1][1] = jac[1];*/
				Jacobi[0][1] = x0[1];
				Jacobi[1][1] = x0[0];
				/*cout << "Якоби: " << endl;
				PrintingArray(Jacobi, 2);*/
				//готовая матрица частных производных
				Jacobi = Inverse(Jacobi, 2);
				if (Jacobi == nullptr) goto newiter;
				Fxk = fsys(x0[0], x0[1]);
				Fxk = MVMultiplication(Jacobi, Fxk, 2);
				xk = Difference(x0, Fxk, 2);
				/*cout << "Вектор xk: ";
				cout << xk[0] << "  " << xk[1] << endl;
				cout << endl;*/
				dx = Difference(x0, xk, 2);
				normdif = NormOfVector(dx, 2);
				//cout << "Норма вектора dx: " << normdif << endl;
				x0 = xk;
			} //выход из цикла, либо если counter = 30, либо итерационный процесс сошёлся

			if (counter < 30)//если процесс сошёлся
			{
				if (n == 0) //если решение первое
				{
					Root = new tVector[1];
					Root[0] = new datatype[2];
					Root[0][0] = xk[0];
					Root[0][1] = xk[1];
					n++;
				}
				else //иначе проверяем сходство с данными в массиве
				{
					for (size_t ii = 0; ii < n; ii++)
					{
						if (abs(Root[ii][0] - xk[0]) < eps  && abs(Root[ii][1] - xk[1]) < eps) //значения х и у отличаются
						{																    	// более, чем на eps
							flag = 1; //есть совпадение
						}
					}
					if (flag == 0)   //если совпадений не было, то добавляем размерность массива, иначе ничего не делаем
					{
						n++;
						tMatrix boofRoot = Root;
						Root = nullptr;
						Root = new tVector[n];
						for (size_t k = 0; k < n; k++)
						{
							Root[k] = new datatype[2];
						}
						for (size_t k = 0; k < n - 1; k++)
						{
							for (size_t jj = 0; jj < 2; jj++)
							{
								Root[k][jj] = boofRoot[k][jj];
							}
							delete[] boofRoot[k];  //ОПАСНО ИМХО
						}
						delete boofRoot;
						Root[n - 1][0] = xk[0];
						Root[n - 1][1] = xk[1];
						/*cout << "Записанный вектор xk в массив: " << endl;
						cout<< Root[n - 1][0] <<"  "<< Root[n - 1][1] <<endl;*/
					}
				}
			}//а если не сошёлся, то переходим к следующему узлу
			//cin.get();
		newiter:
			if (Jacobi != nullptr)
			{
				DeletingArray(Jacobi, 2);
			}

			delete[] x0;
			//delete[] jac;
			//delete[] jacboof;
			delete[] Fxk;
			delete[] dx;
			/*cout << xk[0] << "  " << xk[1] << endl;
			delete[] xk;*/
		}

	}
	if (Root == nullptr)
	{
		cout << "Решений нет!" << endl;
		return nullptr;
	}
	else
	{
		cout << "Всего решений: " << n << endl;
	}
	for (size_t i = 0; i < n; i++)
	{
		cout << i + 1 << ") x = " << Root[i][0] << " , y = " << Root[i][1] << endl;
	}

	cout << "Всего произведено " << sumiter << " итераций" << endl;
	cout << endl;
	return Root;
}


tMatrix FindRoot_Sys2(functsys &fsys, datatype a, datatype b)
{
	size_t sumiter = 0;
	cout << "НЬЮТОН АНАЛИТИЧЕСКИЙ" << endl;
	datatype hh;
	size_t n = 0;   //число решений
	tMatrix Root = nullptr; //массив решений. Координаты записаны по строкам, как х и у соответственно
	datatype x, y = -a;
	cout << "Введите шаг разбиения: ";
	cin >> hh;
	cout << endl << endl;
	datatype alliter = (2 * a) / (hh)-1; //всего узлов на квадрате со стороной a и шагом 0.1
	for (size_t i = 0; i < alliter; i++) //ходим по оси у
	{

		y += hh;
		x = -10;
		for (size_t j = 0; j < alliter; j++) //ходим по оси х по узлам
		{
			bool flag = 0; //изначально 0 - нет совпадений
			x += hh;
			size_t counter = 0; //счетчик числа итераций метода Ньютона
			datatype normdif = 10; //для нормы вектора между итерациями
			tVector x0 = new datatype[2];
			tVector xk = new datatype[2];
			//tVector jac = new datatype[2];
			//tVector jacboof = new datatype[2];
			tVector Fxk = new datatype[2];
			tVector dx = new datatype[2];
			tMatrix Jacobi = nullptr;
			CreatingArray(Jacobi, 2);
			x0[0] = x;					//x0 - вектор {x,y}
			x0[1] = y;					//x0 - вектор {x,y}

			while (normdif > eps7 && counter < 31)  //если не сходится за 30 итераций, то выход из цикла
			{									//и переход к следующему узлу
				sumiter++;
				counter++;
				//строим матрицу Якоби в точке x0
				//jac = fsys[(x + eps, y) - fsys(x, y)] / eps;					//сначала производные по х
				/*jac = fsys(x0[0] + eps, x0[1]);
				jacboof = fsys(x0[0], x0[1]);
				jac = Difference(jac, jacboof, 2);
				jac = NumberVectorMultiplication(1 / eps, jac, 2);
				Jacobi[0][0] = jac[0];
				Jacobi[1][0] = jac[1];*/
				Jacobi[0][0] = 2 * x0[0] + 1;
				Jacobi[1][0] = 2 * x0[1] + 1;
				//jac = fsys[(x, y + eps) - fsys(x, y)] / eps;					//и по у
				/*jac = fsys(x0[0], x0[1] + eps);
				jacboof = fsys(x0[0], x0[1]);
				jac = Difference(jac, jacboof, 2);
				jac = NumberVectorMultiplication(1 / eps, jac, 2);
				Jacobi[0][1] = jac[0];
				Jacobi[1][1] = jac[1];*/
				Jacobi[0][1] = 2 * x0[0] + x0[1];
				Jacobi[1][1] = 2 * x0[1] + x0[0];
				/*cout << "Якоби: " << endl;
				PrintingArray(Jacobi, 2);*/
				//готовая матрица частных производных
				Jacobi = Inverse(Jacobi, 2);
				if (Jacobi == nullptr) goto newiter;
				Fxk = fsys(x0[0], x0[1]);
				Fxk = MVMultiplication(Jacobi, Fxk, 2);
				xk = Difference(x0, Fxk, 2);
				/*cout << "Вектор xk: ";
				cout << xk[0] << "  " << xk[1] << endl;
				cout << endl;*/
				dx = Difference(x0, xk, 2);
				normdif = NormOfVector(dx, 2);
				//cout << "Норма вектора dx: " << normdif << endl;
				x0 = xk;
			} //выход из цикла, либо если counter = 30, либо итерационный процесс сошёлся

			if (counter < 30)//если процесс сошёлся
			{
				if (n == 0) //если решение первое
				{
					Root = new tVector[1];
					Root[0] = new datatype[2];
					Root[0][0] = xk[0];
					Root[0][1] = xk[1];
					n++;
				}
				else //иначе проверяем сходство с данными в массиве
				{
					for (size_t ii = 0; ii < n; ii++)
					{
						if (abs(Root[ii][0] - xk[0]) < eps  && abs(Root[ii][1] - xk[1]) < eps) //значения х и у отличаются
						{																    	// более, чем на eps
							flag = 1; //есть совпадение
						}
					}
					if (flag == 0)   //если совпадений не было, то добавляем размерность массива, иначе ничего не делаем
					{
						n++;
						tMatrix boofRoot = Root;
						Root = nullptr;
						Root = new tVector[n];
						for (size_t k = 0; k < n; k++)
						{
							Root[k] = new datatype[2];
						}
						for (size_t k = 0; k < n - 1; k++)
						{
							for (size_t jj = 0; jj < 2; jj++)
							{
								Root[k][jj] = boofRoot[k][jj];
							}
							delete[] boofRoot[k];  //ОПАСНО ИМХО
						}
						delete boofRoot;
						Root[n - 1][0] = xk[0];
						Root[n - 1][1] = xk[1];
						/*cout << "Записанный вектор xk в массив: " << endl;
						cout<< Root[n - 1][0] <<"  "<< Root[n - 1][1] <<endl;*/
					}
				}
			}//а если не сошёлся, то переходим к следующему узлу
			//cin.get();
		newiter:
			if (Jacobi != nullptr)
			{
				DeletingArray(Jacobi, 2);
			}

			delete[] x0;
			//delete[] jac;
			//delete[] jacboof;
			delete[] Fxk;
			delete[] dx;
			/*cout << xk[0] << "  " << xk[1] << endl;
			delete[] xk;*/
		}

	}
	if (Root == nullptr)
	{
		cout << "Решений нет!" << endl;
		return nullptr;
	}
	else
	{
		cout << "Всего решений: " << n << endl;
	}
	for (size_t i = 0; i < n; i++)
	{
		cout << i + 1 << ") x = " << Root[i][0] << " , y = " << Root[i][1] << endl;
	}
	cout << "Всего произведено " << sumiter << " итераций" << endl;
	cout << endl;
	return Root;
}


void Runge_Kutt(dsys &fsys, tVector &X, int n)
{
	cout << "Runge_Kutt" << endl;
	tVector X0 = VectCopy(X, n);
	tVector X0pr = nullptr;
	unsigned int key = 0;
	datatype tao1 = tao;
	datatype tj = 0;

	char FileName[] = "Data_Runge-Kutt.txt";
	char FileName2[] = "RK_tao n.txt";
	char FileName3[] = "RK_eps n.txt";
	ofstream OUT(FileName);
	ofstream OUT2(FileName2);
	ofstream OUT3(FileName3);
	OUT.setf(ios::fixed);
	OUT2.setf(ios::fixed);
	OUT3.setf(ios::fixed);
	OUT << "{";
	OUT2 << "{";
	OUT3 << "{";

	//OUT << "{" << X0[0] << "," << X0[1] << "}"<<",";

	OUT << "{";
	for (size_t i = 0; i < n; i++)
	{
		OUT << X[i];
		if (i != n - 1)
		{
			OUT << ",";
		}
	}
	OUT << "},";

	cout << "Control step?" << endl << "YOUR CHOISE : ";
	cin >> key;
	cout << endl;
	switch (key)
	{
	case 0:
	{
		tVector k1 = nullptr;
		tVector k2 = nullptr;
		tVector k3 = nullptr;
		tVector k4 = nullptr;
		tVector Xn = nullptr;
		/*CreatingVector(k1, n);
		CreatingVector(k2, n);
		CreatingVector(k3, n);
		CreatingVector(k4, n);*/



		for (size_t i = 1; i < 500; i++)
		{
			k1 = fsys(X0, n);
			tVector k1pr = NumberVectorMultiplication(tao / 2, k1, n);
			tVector k1prpr = Summ(X0, k1pr, n);
			k2 = fsys(k1prpr, n);
			tVector k2pr = NumberVectorMultiplication(tao / 2, k2, n);
			tVector k2prpr = Summ(X0, k2pr, n);
			k3 = fsys(k2prpr, n);
			tVector k3pr = NumberVectorMultiplication(tao, k3, n);
			tVector k3prpr = Summ(X0, k3pr, n);
			k4 = fsys(k3prpr, n);

			DeletingVector(k1pr);
			DeletingVector(k1prpr);
			DeletingVector(k2pr);
			DeletingVector(k2prpr);
			DeletingVector(k3pr);
			DeletingVector(k3prpr);

			k1pr = NumberVectorMultiplication(tao / 6, k1, n);
			k2pr = NumberVectorMultiplication(tao / 3, k2, n);
			k3pr = NumberVectorMultiplication(tao / 3, k3, n);
			tVector k4pr = NumberVectorMultiplication(tao / 6, k4, n);
			tVector K = Summ(k1pr, k2pr, k3pr, k4pr, n);
			DeletingVector(k1);
			DeletingVector(k2);
			DeletingVector(k3);
			DeletingVector(k4);
			//K = Summ(K, k3pr, n);
			//K = Summ(K, k4pr, n);
			//K = NumberVectorMultiplication(tao, K, n);
			Xn = Summ(X0, K, n);
			//OUT << "{" << Xn[0] << "," << Xn[1] << "}";
			OUT << "{";
			for (size_t i = 0; i < n; i++)
			{
				OUT << Xn[i];
				if (i != n - 1)
				{
					OUT << ",";
				}
			}
			OUT << "}";
			if (i % 499 != 0) OUT << ",";
			DeletingVector(X0);
			X0 = nullptr;
			X0 = VectCopy(Xn, n);
			DeletingVector(K);
			DeletingVector(Xn);
			DeletingVector(k1pr);
			DeletingVector(k2pr);
			DeletingVector(k3pr);
			DeletingVector(k4pr);

		}
		OUT << "}";
		cout << "End Runge_Kutt" << endl;
		OUT.close();


		DeletingVector(X0);
		break;
	}
	case 1:
	{
		tVector k1 = nullptr;
		tVector k2 = nullptr;
		tVector k3 = nullptr;
		tVector k4 = nullptr;
		tVector Xn = nullptr;
		tVector Xn1 = nullptr;
		tVector DX = nullptr;
		datatype norm = 0;

		/*CreatingVector(k1, n);
		CreatingVector(k2, n);
		CreatingVector(k3, n);
		CreatingVector(k4, n);*/


	beginning:
		while (tj < 4)
		{
			k1 = fsys(X0, n);
			tVector k1pr = NumberVectorMultiplication(tao1 / 2, k1, n);
			tVector k1prpr = Summ(X0, k1pr, n);
			k2 = fsys(k1prpr, n);
			tVector k2pr = NumberVectorMultiplication(tao1 / 2, k2, n);
			tVector k2prpr = Summ(X0, k2pr, n);
			k3 = fsys(k2prpr, n);
			tVector k3pr = NumberVectorMultiplication(tao1, k3, n);
			tVector k3prpr = Summ(X0, k3pr, n);
			k4 = fsys(k3prpr, n);

			DeletingVector(k1pr);
			DeletingVector(k1prpr);
			DeletingVector(k2pr);
			DeletingVector(k2prpr);
			DeletingVector(k3pr);
			DeletingVector(k3prpr);

			k1pr = NumberVectorMultiplication(tao1 / 6, k1, n);
			k2pr = NumberVectorMultiplication(tao1 / 3, k2, n);
			k3pr = NumberVectorMultiplication(tao1 / 3, k3, n);
			tVector k4pr = NumberVectorMultiplication(tao1 / 6, k4, n);
			tVector K = Summ(k1pr, k2pr, k3pr, k4pr, n);
			DeletingVector(k1);
			DeletingVector(k2);
			DeletingVector(k3);
			DeletingVector(k4);

			Xn = Summ(X0, K, n); //вычислено приближение при шаге tao
			//DeletingVector(X0);
			//X0 = nullptr;
			//X0 = VectCopy(Xn,n);
			DeletingVector(K);
			//DeletingVector(Xn);
			DeletingVector(k1pr);
			DeletingVector(k2pr);
			DeletingVector(k3pr);
			DeletingVector(k4pr);



			X0pr = VectCopy(X0, n); //вектор начальных координат
			tao1 = tao1 / 2; //подробили шаг
			for (size_t i = 0; i < 2; i++)
			{
				k1 = fsys(X0pr, n);
				tVector k1pr = NumberVectorMultiplication(tao1 / 2, k1, n);
				tVector k1prpr = Summ(X0pr, k1pr, n);
				k2 = fsys(k1prpr, n);
				tVector k2pr = NumberVectorMultiplication(tao1 / 2, k2, n);
				tVector k2prpr = Summ(X0pr, k2pr, n);
				k3 = fsys(k2prpr, n);
				tVector k3pr = NumberVectorMultiplication(tao1, k3, n);
				tVector k3prpr = Summ(X0pr, k3pr, n);
				k4 = fsys(k3prpr, n);

				DeletingVector(k1pr);
				DeletingVector(k1prpr);
				DeletingVector(k2pr);
				DeletingVector(k2prpr);
				DeletingVector(k3pr);
				DeletingVector(k3prpr);

				k1pr = NumberVectorMultiplication(tao1 / 6, k1, n);
				k2pr = NumberVectorMultiplication(tao1 / 3, k2, n);
				k3pr = NumberVectorMultiplication(tao1 / 3, k3, n);
				tVector k4pr = NumberVectorMultiplication(tao1 / 6, k4, n);
				tVector K = Summ(k1pr, k2pr, k3pr, k4pr, n);
				DeletingVector(k1);
				DeletingVector(k2);
				DeletingVector(k3);
				DeletingVector(k4);

				Xn1 = Summ(X0pr, K, n);
				DeletingVector(X0pr);
				//X0pr = nullptr;
				X0pr = VectCopy(Xn1, n);    //вычислено полуприближение при шаге tao/2 при i=0
				//и нормальное приближение при i=1
				DeletingVector(K);
				DeletingVector(Xn1);
				DeletingVector(k1pr);
				DeletingVector(k2pr);
				DeletingVector(k3pr);
				DeletingVector(k4pr);

			} //имеем 2 приближения: 1-но шажное и 2-х шажное (НЕ ЗАБЫТЬ УДАЛИТЬ!!)
			tao1 = tao1 * 2; //вернул исходный шаг

			cout<<"step is "<<tao1<<endl;

			DX = Difference(X0pr, Xn, n);
			norm = NormOfVector(DX, n);
			DeletingVector(Xn);
			norm = norm / 15;   //страница 16 методички
			cout<<norm<<endl;
			DeletingVector(DX);
			//OUT2 << tao1 << ",";
			//OUT3 << norm << ",";
			if (norm < 1.e-3 && norm > 1.e-6)  // Делаю как сам вижу, в методе непонятно написано
			{
				//tao1 = 2 * tao1; //возвращаем исходный шаг
				OUT << "{";
				for (size_t i = 0; i < n; i++)
				{
					OUT << X0pr[i];
					if (i != n - 1)
					{
						OUT << ",";
					}
				}
				OUT << "},";
				OUT2 << tao1 << ",";
				OUT3 << norm << ",";
				tj = tj + tao1;
				DeletingVector(X0);
				X0 = VectCopy(X0pr, n);
				DeletingVector(X0pr);
				goto beginning;
			}
			else
				if (norm < 1.e-6)  //если решение очень точное
				{
					tao1 = 3 * tao1; //удлиняем шаг в 3 раза
					DeletingVector(X0pr);
					goto beginning;
				}
				else //если решение очень неточное
				{
					tao1 = tao1 / 2;
					DeletingVector(X0pr);
					goto beginning;
				}

		}

		OUT << "{";
		OUT2 << tao1 << "}";
		OUT3 << norm << "}";
		for (size_t i = 0; i < n; i++) //дублирую запись последней точки в массив значений
		{							// потому что почему бы и нет, чтоб не мучиться с последней запятой
			OUT << X0[i];
			if (i != n - 1)
			{
				OUT << ",";
			}
		}
		OUT << "}";
		OUT << "}";
		cout << "End Runge_Kutt" << endl;
		OUT.close();
		OUT2.close();
		OUT3.close();


		DeletingVector(X0);
		break;
	}
	}


}

tVector Summ(tVector &X1, tVector &X2, tVector &X3, int n)
{
	tVector X = nullptr;
	CreatingVector(X, n);
	for (size_t i = 0; i < n; i++)
	{
		X[i] = X1[i] + X2[i] + X3[i];
	}
	return X;
}

tVector Summ(tVector &X1, tVector &X2, tVector &X3, tVector &X4, int n)
{
	tVector X = nullptr;
	CreatingVector(X, n);
	for (size_t i = 0; i < n; i++)
	{
		X[i] = X1[i] + X2[i] + X3[i] + X4[i];
	}
	return X;
}

tVector Summ(tVector &X1, tVector &X2, tVector &X3, tVector &X4, tVector&X5, int n)
{
	tVector X = nullptr;
	CreatingVector(X, n);
	for (size_t i = 0; i < n; i++)
	{
		X[i] = X1[i] + X2[i] + X3[i] + X4[i] + X5[i];
	}
	return X;
}

void Predictor_Corrector(dsys &fsys, tVector &X, int n)
{
	cout << "Predictor-Corrector method" << endl;
	tVector X0 = VectCopy(X, n);

	char FileName[] = "Data_Predictor-Corrector.txt";
	ofstream OUT(FileName);
	OUT << "{{";
	for (size_t i = 0; i < n; i++)
	{
		OUT << X[i];
		if (i != n - 1)
		{
			OUT << ",";
		}
	}
	OUT << "},";

	tVector k1 = nullptr;
	tVector k2 = nullptr;
	tVector k3 = nullptr;
	tVector k4 = nullptr;
	tVector Xn = nullptr;

	tMatrix Xinitials = new tVector[4];//выделяю память под 4 вектора, нужных для начала итерационного процесса
	//for (size_t i = 0; i < 4; i++)
	//{
	//	//Xinitials[i] = new datatype[2];
	//	Xinitials[i] = VectCopy(X, n);
	//}
	//PrintingArray(Xinitials, n);
	Xinitials[0] = fsys(X0, n);

	for (size_t i = 1; i < 4; i++) //по м-ду Р-К вычисляю 4 раза вектор y, записываю в массив Xinitials
	{
		k1 = fsys(X0, n);
		tVector k1pr = NumberVectorMultiplication(tao / 2, k1, n);
		tVector k1prpr = Summ(X0, k1pr, n);
		k2 = fsys(k1prpr, n);
		tVector k2pr = NumberVectorMultiplication(tao / 2, k2, n);
		tVector k2prpr = Summ(X0, k2pr, n);
		k3 = fsys(k2prpr, n);
		tVector k3pr = NumberVectorMultiplication(tao, k3, n);
		tVector k3prpr = Summ(X0, k3pr, n);
		k4 = fsys(k3prpr, n);

		DeletingVector(k1pr);
		DeletingVector(k1prpr);
		DeletingVector(k2pr);
		DeletingVector(k2prpr);
		DeletingVector(k3pr);
		DeletingVector(k3prpr);

		k1pr = NumberVectorMultiplication(tao / 6, k1, n);
		k2pr = NumberVectorMultiplication(tao / 3, k2, n);
		k3pr = NumberVectorMultiplication(tao / 3, k3, n);
		tVector k4pr = NumberVectorMultiplication(tao / 6, k4, n);
		DeletingVector(k1);
		DeletingVector(k2);
		DeletingVector(k3);
		DeletingVector(k4);
		tVector K = Summ(k1pr, k2pr, k3pr, k4pr, n);
		//K = Summ(K, k3pr, n);
		//K = Summ(K, k4pr, n);
		//K = NumberVectorMultiplication(tao, K, n);
		Xn = Summ(X0, K, n);
		DeletingVector(K);
		DeletingVector(X0);
		Xinitials[i] = VectCopy(Xn, n); //копирую в массив
		//X0 = nullptr;
		X0 = VectCopy(Xn, n);
		DeletingVector(Xn);
		DeletingVector(k1pr);
		DeletingVector(k2pr);
		DeletingVector(k3pr);
		DeletingVector(k4pr);
	}

	/*for (size_t i = 0; i < 4; i++)
	{
		cout << Xinitials[i][0] << " " << Xinitials[i][1] << endl;
	}*/

	//реализую формулы 1.18, имея 4 начальных вектора(f заменяю на k, нумерую слева направо)
	for (size_t i = 4; i < 500; i++)
	{
		Xn = nullptr;
		k1 = fsys(Xinitials[3], n); //fn
		//PrintingVector(k1);
		k2 = fsys(Xinitials[2], n); //fn-1
		//PrintingVector(k2);
		k3 = fsys(Xinitials[1], n); //fn-2
		//PrintingVector(k3);
		k4 = fsys(Xinitials[0], n); //fn-3
		//PrintingVector(k4);
		tVector k1prime = NumberVectorMultiplication(55 * tao / 24, k1, n);
		tVector k2prime = NumberVectorMultiplication(-59 * tao / 24, k2, n);
		tVector k3prime = NumberVectorMultiplication(37 * tao / 24, k3, n);
		tVector k4prime = NumberVectorMultiplication(-9 * tao / 24, k4, n);
		/*DeletingVector(k1);
		DeletingVector(k2);
		DeletingVector(k3);*/
		DeletingVector(k4);
		tVector K = Summ(k1prime, k2prime, k3prime, k4prime, n);
		DeletingVector(k1prime);
		DeletingVector(k2prime);
		DeletingVector(k3prime);
		DeletingVector(k4prime);
		tVector Yzero = Summ(Xinitials[3], K, n);
		DeletingVector(K);
		tVector Fzero = fsys(Yzero, n);
		DeletingVector(Yzero);
		k1prime = NumberVectorMultiplication(9 * tao / 24, Fzero, n);
		DeletingVector(Fzero);
		k2prime = NumberVectorMultiplication(19 * tao / 24, k1, n);
		k3prime = NumberVectorMultiplication(-5 * tao / 24, k2, n);
		k4prime = NumberVectorMultiplication(1 * tao / 24, k3, n);
		DeletingVector(k1);
		DeletingVector(k2);
		DeletingVector(k3);
		K = Summ(k1prime, k2prime, k3prime, k4prime, n);
		DeletingVector(k1prime);
		DeletingVector(k2prime);
		DeletingVector(k3prime);
		DeletingVector(k4prime);
		Xn = Summ(Xinitials[3], K, n);
		DeletingVector(K);
		OUT << "{";
		for (size_t j = 0; j < n; j++)
		{
			OUT << Xn[j];
			if (j != n - 1)
			{
				OUT << ",";
			}
		}
		OUT << "}";
		if (i % 499 != 0) OUT << ",";

		DeletingVector(Xinitials[0]);
		Xinitials[0] = Xinitials[1];
		Xinitials[1] = Xinitials[2];
		Xinitials[2] = Xinitials[3];
		Xinitials[3] = VectCopy(Xn, n);
		DeletingVector(Xn);

	}
	OUT << "}";
	cout << "End of Predictor-Corrector method" << endl;
	OUT.close();

	for (size_t i = 0; i < 4; i++)
	{
		delete[] Xinitials[i];
	}
	delete[] Xinitials;


	/*DeletingVector(k1prime);
	DeletingVector(k2prime);
	DeletingVector(k3prime);
	DeletingVector(k4prime);*/


}

void Pendulum_Euler(dsys &fsys, tVector &X, int n)
{
	cout << "Oscillator ( Euler method )" << endl;
	tVector X0 = VectCopy(X, n);

	char FileName[] = "Oscillator Euler.txt";
	ofstream OUT(FileName);
	OUT << "{{";
	for (size_t i = 0; i < n; i++)
	{
		OUT << X[i];
		if (i != n - 1)
		{
			OUT << ",";
		}
	}
	OUT << "},";
	tMatrix C = nullptr;//задаю матрицу итераций для явного м-да Эйлера
	CreatingArray(C, n);
	C[0][0] = 1; C[0][1] = tao; C[1][0] = -tao * 200 / 3; C[1][1] = 1;
	tVector Xn = nullptr;
	for (size_t i = 1; i < 501; i++)
	{
		OUT << "{";
		Xn = MVMultiplication(C, X0, n);
		for (size_t j = 0; j < n; j++)
		{
			OUT << Xn[j];
			if (j != n - 1)
			{
				OUT << ",";
			}
		}
		OUT << "}";
		if (i % 500 != 0) OUT << ",";

		DeletingVector(X0);
		X0 = VectCopy(Xn, n);
		DeletingVector(Xn);
	}
	OUT << "}";
	cout << "End ofOscillator ( Euler method )" << endl;
	OUT.close();
	DeletingVector(X0);
	DeletingArray(C, n);
}



tVector FindRoot_ImplicitEuler(dsys &fsys, tVector &X, int n)
{
	datatype norm = 10;
	tVector X0 = VectCopy(X, n); //в-р начальных значений
	//tVector Solve = nullptr;
	tVector Xn = nullptr;
	tMatrix Jacobi = nullptr;
	tMatrix JacobiInv = nullptr;
	//CreatingArray(JacobiPrime, n);
	tVector dfx = nullptr;
	tVector C = VectCopy(X, n);
	tVector jac1 = nullptr;
	tVector jac2 = nullptr;
	tVector jacboof = nullptr;

	tVector diff = nullptr;
	while (norm > 1.e-3)
	{
		CreatingArray(Jacobi, n);
		jac1 = fsys(X0, n); //хочу столбец частных производных по x
		X0[0] = X0[0] + eps;
		jac2 = fsys(X0, n);
		X0[0] = X0[0] - eps;
		jacboof = Difference(jac2, jac1, n);
		dfx = NumberVectorMultiplication(-tao / eps, jacboof, n); //  1/eps
		Jacobi[0][0] = dfx[0] + 1;
		Jacobi[1][0] = dfx[1];
		DeletingVector(jac1);
		DeletingVector(jac2);
		DeletingVector(jacboof);
		DeletingVector(dfx);
		jac1 = fsys(X0, n); //хочу столбец частных производных по y
		X0[1] = X0[1] + eps;
		jac2 = fsys(X0, n);
		X0[1] = X0[1] - eps;
		jacboof = Difference(jac2, jac1, n);
		//dfx = NumberVectorMultiplication(1 / eps, jacboof, n);
		dfx = NumberVectorMultiplication(-tao / eps, jacboof, n);
		Jacobi[0][1] = dfx[0];
		Jacobi[1][1] = dfx[1] + 1;
		//DeletingVector(jac1);
		DeletingVector(jac2);
		DeletingVector(jacboof);
		DeletingVector(dfx);
		DeletingVector(jac1);
		JacobiInv = Inverse(Jacobi, n);
		if (JacobiInv == nullptr)
		{
			cout << "Matrix of Jacobi is degenerate, calculating is impossible" << endl;
			return X0;
		}
		DeletingArray(Jacobi, n);  //на этом этапе есть только обратная матрица в S-ой точке внутренних итераций
		jac1 = VectCopy(X0, n);
		jac2 = fsys(X0, n);
		jacboof = NumberVectorMultiplication(-tao, jac2, n);
		DeletingVector(jac2);
		jac2 = Summ(jac1, jacboof, n);
		DeletingVector(jac1);
		DeletingVector(jacboof);
		jac1 = Difference(jac2, C, n);  //это выражение F(X_k), которое нужно для итерационного процесса
		DeletingVector(jac2);

		jac2 = MVMultiplication(JacobiInv, jac1, n);
		DeletingArray(JacobiInv, n);
		DeletingVector(jac1);
		Xn = Difference(X0, jac2, n); // (S+1) приближение к точному решению
		DeletingVector(jac2);



		//Xn = X0 - JacobiInv * fsys(X0, n); надо вычислить немного другие величины в fsys(X0), а именно -2тао*fsys(X0)
		//C = NumberVectorMultiplication(2 * tao, jac1, n);

		//jac2 = MVMultiplication(JacobiInv, C, n);//вычисление не привязано к матрице якоби,
		//DeletingVector(C);
		//DeletingArray(JacobiInv, n);//я так сделал, чтобы не вводить лишние переменные. Их и так много
		//Xn = Summ(X0, jac2, n);
		//DeletingVector(jac2);

		diff = Difference(Xn, X0, n);
		norm = NormOfVector(diff, n);
		cout << "Norm of the difference of vectors: " << norm << endl;
		DeletingVector(diff);
		DeletingVector(X0);
		X0 = VectCopy(Xn, n);
		DeletingVector(Xn);

	}
	//Solve = VectCopy(X0, n);
	//DeletingVector(X0);
	DeletingVector(C);
	//DeletingVector(Xn);
	//DeletingVector(dfx);

	//DeletingArray(Jacobi, n);
	return X0;

}

void ImplicitEuler(dsys &fsys, tVector &X, int n)
{
	cout << "Implicit Euler" << endl;
	tVector X0 = VectCopy(X, n);

	char FileName[] = "Implicit Euler.txt";
	ofstream OUT(FileName);
	OUT << "{{";
	for (size_t i = 0; i < n; i++)
	{
		OUT << X[i];
		if (i != n - 1)
		{
			OUT << ",";
		}
	}
	OUT << "},";

	tVector Xn = nullptr;
	for (size_t i = 1; i < 10001; i++)
	{
		cout << "Calculating of nonlinear equation" << endl;
		OUT << "{";
		Xn = FindRoot_ImplicitEuler(fsys, X0, n);
		for (size_t j = 0; j < n; j++)
		{
			OUT << Xn[j];
			if (j != n - 1)
			{
				OUT << ",";
			}
		}
		OUT << "}";
		if (i % 10000 != 0) OUT << ",";

		DeletingVector(X0);
		X0 = VectCopy(Xn, n);
		DeletingVector(Xn);
	}
	OUT << "}";
	cout << "End of Implicit Euler " << endl;
	OUT.close();
	DeletingVector(X0);
}


//tVector FindRoot_SymPlan_Prime(dsys &fsys, tVector &X, int n)//решение нелинейного уравнения для симметричной схемы (работает онли для размерности 2)
//{
//	datatype norm = 10;
//	tVector X0 = VectCopy(X, n); //в-р начальных значений
//	tVector Solve = nullptr;
//	tVector Xn = nullptr;
//	tMatrix Jacobi = nullptr;
//	tMatrix JacobiInv = nullptr;
//	//CreatingArray(JacobiPrime, n);
//	tVector dfx = nullptr;
//	tVector C = nullptr;
//	tVector jac1 = nullptr;
//	tVector jac2 = nullptr;
//	tVector jacboof = nullptr;
//	tVector diff = nullptr;
//	while (norm > 1.e-1)
//	{
//		CreatingArray(Jacobi, n);
//
//		jac1 = fsys(X0, 2); //хочу столбец частных производных по y
//	    //(1/(4+((200*tao*tao)/3)))
//		//(1 / (-(50 / 3) + (1 / (tao*tao))))
//		CreatingArray(JacobiInv, 2);
//		JacobiInv[0][0] = (1 / (4 + ((200 * tao*tao) / 3)))*(-2);
//		JacobiInv[1][0] = (1 / (4 + ((200 * tao*tao) / 3)))*(-200 * tao / 3);
//		JacobiInv[0][1] = (1 / (4 + ((200 * tao*tao) / 3)))*(tao);
//		JacobiInv[1][1] = (1 / (4 + ((200 * tao*tao) / 3)))*(-2);
//		if (JacobiInv == nullptr)
//		{
//			cout << "Matrix of Jacobi is degenerate, calculating is impossible" << endl;
//			return X0;
//		}
//		DeletingArray(Jacobi, n);
//		//Xn = X0 - JacobiInv * fsys(X0, n); надо вычислить немного другие величины в fsys(X0), а именно -2тао*fsys(X0)
//		C = NumberVectorMultiplication(-2*tao, jac1, n);
//		DeletingVector(jac1);
//		jac2 = MVMultiplication(JacobiInv, C, n);//вычисление не привязано к матрице якоби,
//		DeletingVector(C);
//		DeletingArray(JacobiInv, n);//я так сделал, чтобы не вводить лишние переменные. Их и так много
//		Xn = Difference(X0, jac2, n);
//		DeletingVector(jac2);
//
//		diff = Difference(Xn, X0, n);
//		norm = NormOfVector(diff, n);
//		//cout << "Norm of the difference of vectors: " << norm << endl;
//		DeletingVector(diff);
//		DeletingVector(X0);
//		X0 = VectCopy(Xn, n);
//		DeletingVector(Xn);
//
//	}
//	Solve = VectCopy(X0, n);
//	DeletingVector(X0);
//	//DeletingVector(Xn);
//	//DeletingVector(dfx);
//
//	//DeletingArray(Jacobi, n);
//	return Solve;
//}


//tVector FindRoot_SymPlan_Prime(dsys &fsys, tVector &X, int n)//решение нелинейного уравнения для симметричной схемы (работает онли для размерности 2)
//{
//	datatype norm = 10;
//	tVector X0 = VectCopy(X, n); //в-р начальных значений
//	tVector X0prime = VectCopy(X, n); //в-р начальных значений неизменяемый
//	tVector Solve = nullptr;
//	tVector Xn = nullptr;
//	tMatrix Jacobi = nullptr;
//	tMatrix JacobiInv = nullptr;
//	//CreatingArray(JacobiPrime, n);
//	tVector dfx = nullptr;
//	tVector C = nullptr;
//	tVector jac1 = nullptr;
//	tVector jac2 = nullptr;
//	tVector jac3 = nullptr;
//	tVector jac4 = nullptr;
//	tVector jacboof = nullptr;
//	tVector diff = nullptr;
//	while (norm > 0.7)
//	{
//		CreatingArray(Jacobi, n);
//
//		jac1 = fsys(X0, 2); //хочу столбец промежуточных решений
//		jac2 = fsys(X0prime, 2);
//		jac1= NumberVectorMultiplication(-tao, jac1, n);
//		jac2 = NumberVectorMultiplication(-tao, jac2, n);
//		//(1/(4+((200*tao*tao)/3)))
//		//(1 / (-(50 / 3) + (1 / (tao*tao))))
//		CreatingArray(JacobiInv, 2);
//		JacobiInv[0][0] = (1 / (4 + ((200 * tao*tao) / 3)))*(-2);
//		JacobiInv[1][0] = (1 / (4 + ((200 * tao*tao) / 3)))*(-200 * tao / 3);
//		JacobiInv[0][1] = (1 / (4 + ((200 * tao*tao) / 3)))*(tao);
//		JacobiInv[1][1] = (1 / (4 + ((200 * tao*tao) / 3)))*(-2);
//		if (JacobiInv == nullptr)
//		{
//			cout << "Matrix of Jacobi is degenerate, calculating is impossible" << endl;
//			return X0;
//		}
//		DeletingArray(Jacobi, n);
//		//Xn = X0 - JacobiInv * fsys(X0, n); надо вычислить немного другие величины в fsys(X0), а именно -2тао*fsys(X0)
//		C = Summ(jac2, jac1, n);
//	    jac3= NumberVectorMultiplication(2, X0, n);
//		jac4= NumberVectorMultiplication(-2, X0prime, n);
//		C=Summ(C, jac3, n);
//		C = Summ(C, jac4, n);
//		DeletingVector(jac1);
//		jac2 = MVMultiplication(JacobiInv, C, n);//вычисление не привязано к матрице якоби,
//
//		/*for (int i = 0; i < 2; i++)
//		{
//			cout << jac2[i] << endl;
//		}*/
//		DeletingVector(C);
//		DeletingArray(JacobiInv, n);//я так сделал, чтобы не вводить лишние переменные. Их и так много
//		Xn = Difference(X0, jac2, n);
//		DeletingVector(jac2);
//		DeletingVector(jac3);
//		DeletingVector(jac4);
//		diff = Difference(Xn, X0, n);
//		norm = NormOfVector(diff, n);
//		//cout << "Norm of the difference of vectors: " << norm << endl;
//		DeletingVector(diff);
//		DeletingVector(X0);
//		X0 = VectCopy(Xn, n);
//		DeletingVector(Xn);
//
//	}
//	Solve = VectCopy(X0, n);
//	DeletingVector(X0);
//	//DeletingVector(Xn);
//	//DeletingVector(dfx);
//
//	//DeletingArray(Jacobi, n);
//	return Solve;
//}


//void Pendulum_SymPlan_Prime(dsys &fsys, tVector &X, int n)
//{
//	cout << "Pendulum ( sym )" << endl;
//	tVector X0 = VectCopy(X, n);
//
//	char FileName[] = "Oscillator sym.txt";
//	ofstream OUT(FileName);
//	OUT << "{{";
//	for (size_t i = 0; i < n; i++)
//	{
//		OUT << X[i];
//		if (i != n - 1)
//		{
//			OUT << ",";
//		}
//	}
//	OUT << "},";
//
//	tVector Xn = nullptr;
//	for (size_t i = 0; i < 100; i++)
//	{
//		//cout << "Calculating of nonlinear equation" << endl;
//		OUT << "{";
//		Xn = FindRoot_SymPlan_Prime(fsys, X0, n);
//		for (size_t j = 0; j < n; j++)
//		{
//			OUT << Xn[j];
//			if (j != n - 1)
//			{
//				OUT << ",";
//			}
//		}
//		OUT << "}";
//		if (i % 99 != 0) OUT << ",";
//
//		DeletingVector(X0);
//		X0 = VectCopy(Xn, n);
//		DeletingVector(Xn);
//	}
//	OUT << "}";
//	cout << "End of Pendulum ( sym )" << endl;
//	OUT.close();
//	DeletingVector(X0);
//}

//void Pendulum_SymPlan_Open(dsys &fsys, tVector &X, int n)
//{
//	cout << "Pendulum ( sym )" << endl;
//	tVector X0 = VectCopy(X, n);
//	tVector Xn = VectCopy(X, n);
//	datatype k;
//	k = 200 / 3;
//	char FileName[] = "Oscillator_sym_open.txt";
//	ofstream OUT(FileName);
//	OUT << "{{";
//	for (size_t i = 0; i < n; i++)
//	{
//		OUT << X[i];
//		if (i != n - 1)
//		{
//			OUT << ",";
//		}
//	}
//	OUT << "},";
//
//
//	//tVector Xn = nullptr;
//	//Xn = VectCopy(X0, n);
//	for (size_t i = 0; i < 100; i++)
//	{
//		datatype x = X0[0];
//		datatype y = X0[1];
//		datatype xyi1= (2 * y - 4 * k*tao*x - k * tao*tao*y) / (2 + k * tao*tao);
//		Xn[1] = xyi1;
//		datatype y_prime = Xn[1];
//		datatype xyi2 = x + tao * ((y + y_prime) / (2));
//		Xn[0] = xyi2;
//		OUT << "{";
//		for (size_t j = 0; j < n; j++)
//		{
//			OUT << Xn[j];
//			if (j != n - 1)
//			{
//				OUT << ",";
//			}
//		}
//		OUT << "}";
//		if (i % 99 != 0) OUT << ",";
//
//		DeletingVector(X0);
//		X0 = VectCopy(Xn, n);
//		DeletingVector(Xn);
//	}
//	OUT << "}";
//	cout << "End of Pendulum" << endl;
//	OUT.close();
//	DeletingVector(X0);
//}


tVector FindRoot_SymPlan(dsys &fsys, tVector &X, int n)//решение нелинейного уравнения для симметричной схемы (работает онли для размерности 2)
{
	datatype norm = 10;
	tVector X0 = VectCopy(X, n); //в-р начальных значений
	//tVector Solve = nullptr;
	tVector Xn = nullptr;
	tMatrix Jacobi = nullptr;
	tMatrix JacobiInv = nullptr;
	//CreatingArray(JacobiPrime, n);
	tVector dfx = nullptr;
	tVector C = nullptr;
	tVector jac1 = nullptr;
	tVector jac2 = nullptr;
	tVector jacboof = nullptr;
	jac1 = NumberVectorMultiplication(-2, X, n); //-2*y_n
	jac2 = fsys(X, n);
	jacboof = NumberVectorMultiplication(-tao, jac2, n);
	DeletingVector(jac2);
	C = Summ(jac1, jacboof, n); //определил постоянную часть выражения F(X)=0;
	DeletingVector(jac1);
	DeletingVector(jacboof);
	tVector diff = nullptr;
	while (norm > 1.e-3)
	{
		CreatingArray(Jacobi, n);
		jac1 = fsys(X0, n); //хочу столбец частных производных по x
		X0[0] = X0[0] + eps;
		jac2 = fsys(X0, n);
		X0[0] = X0[0] - eps;
		jacboof = Difference(jac2, jac1, n);
		dfx = NumberVectorMultiplication(-tao / eps, jacboof, n); //  1/eps
		Jacobi[0][0] = dfx[0] + 2;
		Jacobi[1][0] = dfx[1];
		DeletingVector(jac1);
		DeletingVector(jac2);
		DeletingVector(jacboof);
		DeletingVector(dfx);
		jac1 = fsys(X0, n); //хочу столбец частных производных по y
		X0[1] = X0[1] + eps;
		jac2 = fsys(X0, n);
		X0[1] = X0[1] - eps;
		jacboof = Difference(jac2, jac1, n);
		dfx = NumberVectorMultiplication(-tao / eps, jacboof, n);
		Jacobi[0][1] = dfx[0];
		Jacobi[1][1] = dfx[1] + 2;
		DeletingVector(jac2);
		DeletingVector(jacboof);
		DeletingVector(dfx);
		DeletingVector(jac1);
		JacobiInv = Inverse(Jacobi, n);
		if (JacobiInv == nullptr)
		{
			cout << "Matrix of Jacobi is degenerate, calculating is impossible" << endl;
			return X0;
		}
		DeletingArray(Jacobi, n);  //на этом этапе есть только обратная матрица в S-ой точке внутренних итераций
		jac1 = NumberVectorMultiplication(2, X0, n);
		jac2 = fsys(X0, n);
		jacboof = NumberVectorMultiplication(-tao, jac2, n);
		DeletingVector(jac2);
		jac2 = Summ(jac1, jacboof, n);
		DeletingVector(jac1);
		DeletingVector(jacboof);
		jac1 = Summ(jac2, C, n);  //это выражение F(X_k), которое нужно для итерационного процесса
		DeletingVector(jac2);

		jac2 = MVMultiplication(JacobiInv, jac1, n);
		DeletingArray(JacobiInv, n);
		DeletingVector(jac1);
		Xn = Difference(X0, jac2, n); // (S+1) приближение к точному решению
		DeletingVector(jac2);

		diff = Difference(Xn, X0, n);
		norm = NormOfVector(diff, n);
		//cout << "Norm of the difference of vectors: " << norm << endl;
		DeletingVector(diff);
		DeletingVector(X0);
		X0 = VectCopy(Xn, n);
		DeletingVector(Xn);

	}

	DeletingVector(C);

	return X0;
}

void Pendulum_SymPlan(dsys &fsys, tVector &X, int n)
{
	cout << "Pendulum ( sym )" << endl;
	tVector X0 = VectCopy(X, n);

	char FileName[] = "Oscillator SymPlan.txt";
	ofstream OUT(FileName);
	OUT << "{{";
	for (size_t i = 0; i < n; i++)
	{
		OUT << X[i];
		if (i != n - 1)
		{
			OUT << ",";
		}
	}
	OUT << "},";

	tVector Xn = nullptr;
	for (size_t i = 1; i < 1000001; i++)
	{
		//cout << "Calculating of nonlinear equation" << endl;
		OUT << "{";
		Xn = FindRoot_SymPlan(fsys, X0, n);
		for (size_t j = 0; j < n; j++)
		{
			OUT << Xn[j];
			if (j != n - 1)
			{
				OUT << ",";
			}
		}
		OUT << "}";
		if (i % 1000000 != 0) OUT << ",";

		DeletingVector(X0);
		X0 = VectCopy(Xn, n);
		DeletingVector(Xn);
	}
	OUT << "}";
	cout << "End of Pendulum" << endl;
	OUT.close();
	DeletingVector(X0);
}

void Enclosed_Runge_Kutt(dsys &fsys, tVector &X, int n)
{
	cout << "Enclosed_Runge_Kutt method" << endl;
	tVector X0 = VectCopy(X, n);
	datatype tj = 0;
	datatype normDX = 0;
	datatype tao1 = tao;

	char FileName[] = "Data_Enclosed_Runge-Kutt.txt";
	char FileName2[] = "ERK_tao n.txt";
	char FileName3[] = "ERK_eps n.txt";
	ofstream OUT(FileName);
	ofstream OUT2(FileName2);
	ofstream OUT3(FileName3);
	OUT << "{{";
	OUT2 << "{";
	OUT3 << "{";
	for (size_t i = 0; i < n; i++)
	{
		OUT << X[i];
		if (i != n - 1)
		{
			OUT << ",";
		}
	}
	OUT << "},";

	tVector k1 = nullptr;
	tVector k2 = nullptr;
	tVector k3 = nullptr;
	tVector k4 = nullptr;
	tVector Xn = nullptr;
	tVector Xn1 = nullptr;
	tVector DX = nullptr;
	while (tj < 2)
	{
		beginning:
		k1 = fsys(X0, n);
		tVector k1pr = NumberVectorMultiplication(tao1 / 2, k1, n);
		tVector k1prpr = Summ(X0, k1pr, n);
		k2 = fsys(k1prpr, n);
		tVector k2pr = NumberVectorMultiplication(tao1 / 2, k2, n);
		tVector k2prpr = Summ(X0, k2pr, n);
		k3 = fsys(k2prpr, n);
		tVector k3pr = NumberVectorMultiplication(tao1, k3, n);
		tVector k3prpr = Summ(X0, k3pr, n);
		k4 = fsys(k3prpr, n);

		DeletingVector(k1pr);
		DeletingVector(k1prpr);
		DeletingVector(k2pr);
		DeletingVector(k2prpr);
		DeletingVector(k3pr);
		DeletingVector(k3prpr);

		k1pr = NumberVectorMultiplication(tao1 / 6, k1, n);
		k2pr = NumberVectorMultiplication(tao1 / 3, k2, n);
		k3pr = NumberVectorMultiplication(tao1 / 3, k3, n);
		tVector k4pr = NumberVectorMultiplication(tao1 / 6, k4, n);
		tVector K = Summ(k1pr, k2pr, k3pr, k4pr, n);
		DeletingVector(k4pr);
		Xn = Summ(X0, K, n); //вычислено значение y_n+1
		DeletingVector(K);

		k4pr = NumberVectorMultiplication(tao1 / 15, k4, n);
		tVector k5pr = NumberVectorMultiplication(tao1 / 10, Xn, n);
		K = Summ(k1pr, k2pr, k3pr, k4pr, k5pr, n);
		DeletingVector(k1);
		DeletingVector(k2);
		DeletingVector(k3);
		DeletingVector(k4);
		Xn1 = Summ(X0, K, n);
		DeletingVector(K);
		DeletingVector(k1pr);
		DeletingVector(k2pr);
		DeletingVector(k3pr);
		DeletingVector(k4pr);
		DeletingVector(k5pr);

		DX = Difference(Xn, Xn1, n);
		normDX = NormOfVector(DX, n);
		DeletingVector(DX);
		DeletingVector(Xn1);

		if (normDX<1.e-3 && normDX>1.e-5)
		{
			OUT << "{";
			OUT2 << tao1 << ",";
			OUT3 << normDX << ",";
			for (size_t i = 0; i < n; i++)
			{
				OUT << Xn[i];
				if (i != n - 1)
				{
					OUT << ",";
				}
			}
			OUT << "},";
			tj = tj + tao1; //!!!!!!
			DeletingVector(X0);
			X0 = nullptr;
			X0 = VectCopy(Xn, n);
			DeletingVector(Xn);
		}
		else if (normDX < 1.e-5)
		{
			tao1 = 3 * tao1; //удлиняем шаг в 3 раза
			DeletingVector(Xn);
			goto beginning;
		}
		else
		{
			tao1 = tao1 / 2; //уменьшаем шаг в 2 раза
			DeletingVector(Xn);
			goto beginning;
		}

	}

	tao1 = tao1 - tj + 1; //шаг последней итерации
	//cout<<"Last step is "<<tao1<<endl;
	k1 = fsys(X0, n);
	tVector k1pr = NumberVectorMultiplication(tao1 / 2, k1, n);
	tVector k1prpr = Summ(X0, k1pr, n);
	k2 = fsys(k1prpr, n);
	tVector k2pr = NumberVectorMultiplication(tao1 / 2, k2, n);
	tVector k2prpr = Summ(X0, k2pr, n);
	k3 = fsys(k2prpr, n);
	tVector k3pr = NumberVectorMultiplication(tao1, k3, n);
	tVector k3prpr = Summ(X0, k3pr, n);
	k4 = fsys(k3prpr, n);

	DeletingVector(k1pr);
	DeletingVector(k1prpr);
	DeletingVector(k2pr);
	DeletingVector(k2prpr);
	DeletingVector(k3pr);
	DeletingVector(k3prpr);

	k1pr = NumberVectorMultiplication(tao1 / 6, k1, n);
	k2pr = NumberVectorMultiplication(tao1 / 3, k2, n);
	k3pr = NumberVectorMultiplication(tao1 / 3, k3, n);
	tVector k4pr = NumberVectorMultiplication(tao1 / 6, k4, n);
	tVector K = Summ(k1pr, k2pr, k3pr, k4pr, n);
	DeletingVector(k1);
	DeletingVector(k2);
	DeletingVector(k3);
	DeletingVector(k4);
	Xn = Summ(X0, K, n);
	OUT << "{";
	for (size_t i = 0; i < n; i++)
	{
		OUT << Xn[i];
		if (i != n - 1)
		{
			OUT << ",";
		}
	}
	OUT << "}}";
	OUT2 << tao1 << "}";
	OUT3 << normDX << "}";
	DeletingVector(X0);
	DeletingVector(K);
	DeletingVector(Xn);
	DeletingVector(k1pr);
	DeletingVector(k2pr);
	DeletingVector(k3pr);
	DeletingVector(k4pr);

	cout << "End Enclosed_Runge_Kutt method" << endl;
	OUT.close();
	OUT2.close();
	OUT3.close();
}


void Adams(dsys &fsys, tVector &X, int n)
{
	cout << "Adams-Bushfort method(extrapolation)" << endl;
	tVector X0 = VectCopy(X, n);

	char FileName[] = "Data_Adams-Bushford.txt";
	ofstream OUT(FileName);
	OUT << "{{";
	for (size_t i = 0; i < n; i++)
	{
		OUT << X0[i];
		if (i != n - 1)
		{
			OUT << ",";
		}
	}
	OUT << "},";

	tVector k1 = nullptr;
	tVector k2 = nullptr;
	tVector k3 = nullptr;
	tVector k4 = nullptr;
	tVector Xn = nullptr;

	tMatrix Xinitials = new tVector[4];//выделяю память под 4 вектора, нужных для начала итерационного процесса

	Xinitials[0] = VectCopy(X0, n);
	for (size_t i = 1; i < 4; i++) //по м-ду Р-К вычисляю 4 раза вектор y, записываю в массив Xinitials
	{
		k1 = fsys(X0, n);
		tVector k1pr = NumberVectorMultiplication(tao / 2, k1, n);
		tVector k1prpr = Summ(X0, k1pr, n);
		k2 = fsys(k1prpr, n);
		tVector k2pr = NumberVectorMultiplication(tao / 2, k2, n);
		tVector k2prpr = Summ(X0, k2pr, n);
		k3 = fsys(k2prpr, n);
		tVector k3pr = NumberVectorMultiplication(tao, k3, n);
		tVector k3prpr = Summ(X0, k3pr, n);
		k4 = fsys(k3prpr, n);

		DeletingVector(k1pr);
		DeletingVector(k1prpr);
		DeletingVector(k2pr);
		DeletingVector(k2prpr);
		DeletingVector(k3pr);
		DeletingVector(k3prpr); // k1-k4

		k1pr = NumberVectorMultiplication(tao / 6, k1, n);
		k2pr = NumberVectorMultiplication(tao / 3, k2, n);
		k3pr = NumberVectorMultiplication(tao / 3, k3, n);
		tVector k4pr = NumberVectorMultiplication(tao / 6, k4, n);
		tVector K = Summ(k1pr, k2pr, k3pr, k4pr, n);
		DeletingVector(k1pr);
		DeletingVector(k2pr);
		DeletingVector(k3pr);
		DeletingVector(k4pr);
		//K = Summ(K, k3pr, n);
		//K = Summ(K, k4pr, n);
		//K = NumberVectorMultiplication(tao, K, n);
		Xn = Summ(X0, K, n);
		OUT << "{";
		for (size_t j = 0; j < n; j++)
		{
			OUT << Xn[j];
			if (j != n - 1)
			{
				OUT << ",";
			}
		}
		OUT << "}";
		OUT << ",";
		Xinitials[i] = VectCopy(Xn, n); //копирую в массив
		DeletingVector(X0);
		//X0 = nullptr;
		X0 = VectCopy(Xn, n);
		DeletingVector(K);
		DeletingVector(Xn);


		DeletingVector(k1);
		DeletingVector(k2);
		DeletingVector(k3);
		DeletingVector(k4);

	}
	DeletingVector(X0);  //allright


	//реализую формулы 1.18, имея 4 начальных вектора(f заменяю на k, нумерую слева направо)  ЗДЕСЬ ЧТО_ТО НЕПРАВИЛЬНО
	for (size_t i = 4; i < 500; i++)
	{
		//Xn = nullptr;
		k1 = fsys(Xinitials[3], n); //fn

		k2 = fsys(Xinitials[2], n); //fn-1

		k3 = fsys(Xinitials[1], n); //fn-2

		k4 = fsys(Xinitials[0], n); //fn-3

		tVector k1prime = NumberVectorMultiplication((55 * tao) / 24, k1, n);
		tVector k2prime = NumberVectorMultiplication(-(59 * tao) / 24, k2, n);
		tVector k3prime = NumberVectorMultiplication((37 * tao) / 24, k3, n);
		tVector k4prime = NumberVectorMultiplication(-(9 * tao) / 24, k4, n);


		tVector K = Summ(k1prime, k2prime, k3prime, k4prime, n);
		DeletingVector(k1prime);
		DeletingVector(k2prime);
		DeletingVector(k3prime);
		DeletingVector(k4prime);


		//tVector Yzero = Summ(Xinitials[3], K, n);
		Xn = Summ(Xinitials[3], K, n);
		DeletingVector(K);

		OUT << "{";
		for (size_t j = 0; j < n; j++)
		{
			OUT << Xn[j];
			if (j != n - 1)
			{
				OUT << ",";
			}
		}
		OUT << "}";
		if (i % 499 != 0) OUT << ",";
		//

		DeletingVector(k1);
		DeletingVector(k2);
		DeletingVector(k3);
		DeletingVector(k4);




		DeletingVector(Xinitials[0]);
		Xinitials[0] = Xinitials[1];
		Xinitials[1] = Xinitials[2];
		Xinitials[2] = Xinitials[3];
		Xinitials[3] = VectCopy(Xn, n);
		DeletingVector(Xn);
		//DeletingVector(Yzero);
	}


	OUT << "}";
	cout << "End of Adams-Bushfort method(extrapolation)" << endl;
	OUT.close();


	for (size_t i = 0; i < 4; i++)
	{
		delete[] Xinitials[i];
	}
	delete[] Xinitials;



}


tVector BigExample(tVector& X0, datatype t)
{
	tVector sol = new datatype[2];
	sol[0] = (4 * (X0[0] * X0[0] + X0[1] * X0[1])*X0[1] + 0.16*X0[1]) - 5 * ((X0[0] * X0[0] + X0[1] * X0[1])*(X0[0] * X0[0] + X0[1] * X0[1]) - 0.08*(X0[0] * X0[0] - X0[1] * X0[1]))* (4 * (X0[0] * X0[0] + X0[1] * X0[1])*X0[0] - 0.16*X0[0]) + 0.015*X0[1] * sin(5 * t);
	sol[1] = -(4 * (X0[0] * X0[0] + X0[1] * X0[1])*X0[0] - 0.16*X0[0]) - 5 * ((X0[0] * X0[0] + X0[1] * X0[1])*(X0[0] * X0[0] + X0[1] * X0[1]) - 0.08*(X0[0] * X0[0] - X0[1] * X0[1]))*(4 * (X0[0] * X0[0] + X0[1] * X0[1])*X0[1] + 0.16*X0[1]) + 0.015*X0[1] * sin(5 * t);
	return sol;
}

void SolveBigExampleRK(bigexample &fsys, tVector &X)
{
	cout << "Big Example" << endl;
	tVector X0 = VectCopy(X, 2);
	//tVector X0pr = nullptr;
	datatype dt = tao;
	datatype tj = 0;
	datatype norm = 0;
	int key = 0;

	char FileName[] = "Big Example.txt";
	char FileName2[]="BE tao.txt";
	char FileName3[]="BE eps.txt";
	ofstream OUT(FileName);
	ofstream OUTt(FileName2);
	ofstream OUTe(FileName3);

	OUT.setf(ios::fixed);
	OUTt.setf(ios::fixed);
	OUTe.setf(ios::fixed);


	OUT << "{{";
	for (size_t i = 0; i < 2; i++)
	{
		OUT << X[i];
		if (i != 1)
		{
			OUT << ",";
		}
	}
	OUT << "},";
	OUTt<<"{";
	OUTe<<"{";
	cout<<"Control step?"<<endl;
	cout<<"Your choise: ";
	cin>>key;
	cout<<endl;
	switch (key)
	{
	case 0: //временной интервал 10 секунд
		{
			tVector k1 = nullptr;
			tVector k2 = nullptr;
			tVector k3 = nullptr;
			tVector k4 = nullptr;
			tVector Xn = nullptr;

			for (size_t i = 1; i < 50001; i++)
			{
				k1 = fsys(X0, tj);
				tVector k1pr = NumberVectorMultiplication(tao / 2, k1, 2);
				tVector k1prpr = Summ(X0, k1pr, 2);
				k2 = fsys(k1prpr, tj + tao / 2);
				tVector k2pr = NumberVectorMultiplication(tao / 2, k2, 2);
				tVector k2prpr = Summ(X0, k2pr, 2);
				k3 = fsys(k2prpr, tj + tao / 2);
				tVector k3pr = NumberVectorMultiplication(tao, k3, 2);
				tVector k3prpr = Summ(X0, k3pr, 2);
				k4 = fsys(k3prpr, tj + tao);
				tj = tj + tao;

				DeletingVector(k1pr);
				DeletingVector(k1prpr);
				DeletingVector(k2pr);
				DeletingVector(k2prpr);
				DeletingVector(k3pr);
				DeletingVector(k3prpr);

				k1pr = NumberVectorMultiplication(tao / 6, k1, 2);
				k2pr = NumberVectorMultiplication(tao / 3, k2, 2);
				k3pr = NumberVectorMultiplication(tao / 3, k3, 2);
				tVector k4pr = NumberVectorMultiplication(tao / 6, k4, 2);
				tVector K = Summ(k1pr, k2pr, k3pr, k4pr, 2);
				DeletingVector(k1);
				DeletingVector(k2);
				DeletingVector(k3);
				DeletingVector(k4);
				Xn = Summ(X0, K, 2);
				OUT << "{";
				for (size_t i = 0; i < 2; i++)
				{
					OUT << Xn[i];
					if (i != 1)
					{
						OUT << ",";
					}
				}
				OUT << "}";
				if (i % 50000 != 0) OUT << ",";
				DeletingVector(X0);
				X0 = nullptr;
				X0 = VectCopy(Xn, 2);
				DeletingVector(K);
				DeletingVector(Xn);
				DeletingVector(k1pr);
				DeletingVector(k2pr);
				DeletingVector(k3pr);
				DeletingVector(k4pr);

			}
			OUT << "}";
			cout << "End of Big Example." << endl;
			OUT.close();
			DeletingVector(X0);
			break;
		}
	case 1:
		{
			while (tj<50)
			{
				tVector k1 = nullptr;
				tVector k2 = nullptr;
				tVector k3 = nullptr;
				tVector k4 = nullptr;

				k1=fsys(X0,tj);
				tVector kpr = NumberVectorMultiplication(dt/2, k1, 2);
				tVector kprpr = Summ(X0, kpr, 2);
				k2=fsys(kprpr, tj + dt/2);
				DeletingVector(kpr);
				DeletingVector(kprpr);
				kpr = NumberVectorMultiplication(dt/2, k2, 2);
				kprpr = Summ(X0, kpr, 2);
				k3=fsys(kprpr, tj + dt/2);
				DeletingVector(kpr);
				DeletingVector(kprpr);
				kpr = NumberVectorMultiplication(dt, k3, 2);
				kprpr = Summ(X0, kpr, 2);
				k4=fsys(kprpr, tj + dt);
				DeletingVector(kpr);
				DeletingVector(kprpr);
				tVector k1pr = NumberVectorMultiplication(dt/6, k1, 2);
				tVector k2pr = NumberVectorMultiplication(dt/3, k2, 2);
				tVector k3pr = NumberVectorMultiplication(dt/3, k3, 2);
				tVector k4pr = NumberVectorMultiplication(dt/6, k4, 2);
				DeletingVector(k1);
				DeletingVector(k2);
				DeletingVector(k3);
				DeletingVector(k4);
				tVector K = Summ(k1pr, k2pr, k3pr, k4pr, 2);
				DeletingVector(k1pr);
				DeletingVector(k2pr);
				DeletingVector(k3pr);
				DeletingVector(k4pr);
				tVector Xn = Summ(X0, K, 2); //решение на следующем слое за один шаг

				tVector X0pr = VectCopy(X0, 2);
				dt=dt/2;
				for (int i=0; i<2; i++)
				{
					k1=fsys(X0pr,tj);
					kpr = NumberVectorMultiplication(dt/2, k1, 2);
					kprpr = Summ(X0pr, kpr, 2);
					k2=fsys(kprpr, tj + dt/2);
					DeletingVector(kpr);
					DeletingVector(kprpr);
					kpr = NumberVectorMultiplication(dt/2, k2, 2);
					kprpr = Summ(X0pr, kpr, 2);
					k3=fsys(kprpr, tj + dt/2);
					DeletingVector(kpr);
					DeletingVector(kprpr);
					kpr = NumberVectorMultiplication(dt, k3, 2);
					kprpr = Summ(X0pr, kpr, 2);
					k4=fsys(kprpr, tj + dt);
					DeletingVector(kpr);
					DeletingVector(kprpr);
					k1pr = NumberVectorMultiplication(dt/6, k1, 2);
					k2pr = NumberVectorMultiplication(dt/3, k2, 2);
					k3pr = NumberVectorMultiplication(dt/3, k3, 2);
					k4pr = NumberVectorMultiplication(dt/6, k4, 2);
					DeletingVector(k1);
					DeletingVector(k2);
					DeletingVector(k3);
					DeletingVector(k4);
					K = Summ(k1pr, k2pr, k3pr, k4pr, 2);
					DeletingVector(k1pr);
					DeletingVector(k2pr);
					DeletingVector(k3pr);
					DeletingVector(k4pr);
					tVector Xnpr = Summ(X0pr, K, 2);//решение на слое 1/2
					DeletingVector(X0pr);
					X0pr = VectCopy(Xnpr, 2); //"переход" к следующему слою
					DeletingVector(Xnpr);
					tj = tj + dt;
				}//в X0pr сидит значение на новом слое
				dt = 2*dt;
				tj = tj - dt;
				tVector DX = Difference(X0pr, Xn, 2);
				norm = NormOfVectorC(DX, 2);

				if (norm < 0.1e-8 && norm > 0.99e-8)
				{
					DeletingVector(X0pr);
					OUT<<"{"<<Xn[0]<<","<<Xn[1]<<"},";
					OUTt<<dt<<",";
					OUTe<<fixed<<setprecision(15)<<norm<<",";
					DeletingVector(X0);
					X0 = VectCopy(Xn, 2);
					DeletingVector(Xn);
					tj = tj + dt;
				}else
				if (norm <= 0.99e-8)
				{
					DeletingVector(X0pr);
					OUT<<"{"<<Xn[0]<<","<<Xn[1]<<"},";
					OUTt<<dt<<",";
					OUTe<<fixed<<setprecision(15)<<norm<<",";
					DeletingVector(X0);
					X0 = VectCopy(Xn, 2);
					DeletingVector(Xn);
					//dt=dt*2;
					tj = tj + dt;
					dt = dt*2;
				}else
				if (norm >= 0.1e-8)
				{
					DeletingVector(X0pr);
					DeletingVector(Xn);
					dt=dt/2.5;
//					OUTt<<dt<<",";
//					OUTe<<norm<<",";
				}
			}
			OUTt<<dt<<"}";
			OUTt.close();
			OUTe<<norm<<"}";
			OUTe.close();
			OUT<<X0[0]<<","<<X0[1]<<"}";
			OUT.close();
			cout << "End of Big Example with step control" << endl;
			break;
		}
	}

}


void RKvsAS(dsys &fsys, tVector &X, int n)
{


	//tVector X0pr = nullptr;
	tVector AS = nullptr; //от слова Analytical Solve
	datatype tao1 = 15*tao;
	for (size_t u = 0; u < 20; u++)
	{
		datatype max = 0;
		tVector X0 = VectCopy(X, n);
		cout << "Runge-Kutt vs Analytical Solve ,  h = " << tao1 << endl;
		datatype tj = 0;
		/*char FileName[] = "Data_Runge-Kutt.txt";
		ofstream OUT(FileName);
		OUT << "{{";
		for (size_t i = 0; i < n; i++)
		{
			OUT << X[i];
			if (i != n - 1)
			{
				OUT << ",";
			}
		}
		OUT << "},";*/

		tVector k1 = nullptr;
		tVector k2 = nullptr;
		tVector k3 = nullptr;
		tVector k4 = nullptr;
		tVector Xn = nullptr;

		while (tj < 5)
		{
			k1 = fsys(X0, n);
			tVector k1pr = NumberVectorMultiplication(tao1 / 2, k1, n);
			tVector k1prpr = Summ(X0, k1pr, n);
			k2 = fsys(k1prpr, n);
			tVector k2pr = NumberVectorMultiplication(tao1 / 2, k2, n);
			tVector k2prpr = Summ(X0, k2pr, n);
			k3 = fsys(k2prpr, n);
			tVector k3pr = NumberVectorMultiplication(tao1, k3, n);
			tVector k3prpr = Summ(X0, k3pr, n);
			k4 = fsys(k3prpr, n);

			DeletingVector(k1pr);
			DeletingVector(k1prpr);
			DeletingVector(k2pr);
			DeletingVector(k2prpr);
			DeletingVector(k3pr);
			DeletingVector(k3prpr);

			k1pr = NumberVectorMultiplication(tao1 / 6, k1, n);
			k2pr = NumberVectorMultiplication(tao1 / 3, k2, n);
			k3pr = NumberVectorMultiplication(tao1 / 3, k3, n);
			tVector k4pr = NumberVectorMultiplication(tao1 / 6, k4, n);
			tVector K = Summ(k1pr, k2pr, k3pr, k4pr, n);
			DeletingVector(k1);
			DeletingVector(k2);
			DeletingVector(k3);
			DeletingVector(k4);

			Xn = Summ(X0, K, n);
			AS = dsysAnalSolve(tj + tao1);
			datatype a = abs(Xn[0] - AS[0]);
			datatype b = abs(Xn[1] - AS[1]);
			if (a > b && max < a)
			{
				max = a;
			}
			else if (a <= b && max < b)
			{
				max = b;
			}
			DeletingVector(AS);

			/*OUT << "{";
			for (size_t i = 0; i < n; i++)
			{
				OUT << Xn[i];
				if (i != n - 1)
				{
					OUT << ",";
				}
			}
			OUT << "}";
			if (i % 499 != 0) OUT << ",";*/
			DeletingVector(X0);
			X0 = nullptr;
			X0 = VectCopy(Xn, n);
			DeletingVector(K);
			DeletingVector(Xn);
			DeletingVector(k1pr);
			DeletingVector(k2pr);
			DeletingVector(k3pr);
			DeletingVector(k4pr);
			tj = tj + tao1;
		}
		/*OUT << "}";
		cout << "End of RK method." << endl;
		OUT.close();*/
		DeletingVector(X0);
		cout << "max| y - U | = " << max << endl;
		cout << endl;
		tao1 = tao1 * 0.7;
	}
}


void LUDecomposition_pr(tMatrix &A, int n)
{
	datatype sum1;
	int row_max = 0;
	datatype max = 0;
	tVector boof = nullptr;
	//tVector Standings = nullptr; //в St-gs[i] записан номер строки, которую нужно поставить на место i-й
	//tVector b = nullptr;
	//tVector y = nullptr;
	//tMatrix LU = nullptr;
	//tMatrix Solve = nullptr;
	//CreatingArray(Solve, n);
	//CreatingArray(LU, n); //выделил память под итоговую матрицу
	CreatingVector(Standings, n);
	for (size_t k = 0; k < n; k++)
	{
		Standings[k] = k;
	}

	for (size_t k = 0; k < n; k++) //счетчик стр и стб, с которыми ведётся работа
	{
		row_max = k;
		max = abs(A[k][k]);
		for	(size_t i = k + 1; i < n; i++)
		{
			if (abs(A[i][k])>max)
			{
				row_max = i;
				max = abs(A[i][k]);
			}
		}//ищу максимальный элемент в столбце
		if (row_max != k)
		{
			boof = A[k];
			A[k] = A[row_max];
			A[row_max] = boof;

			boof = A[k];
			A[k] = A[row_max];
			A[row_max] = boof;

			Standings[k] = row_max;
			Standings[row_max] = k;
		}
		for (size_t j = k; j < n; j++) //работа с частями U
		{
			sum1 = 0;
			for (size_t i = 0; i < k; i++)
			{
				sum1 = sum1 + A[k][i] * A[i][j];
			}
			A[k][j] = A[k][j] - sum1;
			if (A[k][k] == 0)
			{
				cout << "LU decomposition is impossible ( det = 0 )" << endl;
				//cout << "Returning nullptr..." << endl;
				//return nullptr;
			}
		}
		for (size_t j = k+1; j < n; j++) //работа с частями L
		{
			sum1 = 0;
			for (size_t i = 0; i < k; i++)
			{
				sum1 = sum1 + A[j][i] * A[i][k];
			}
			A[j][k] = (A[j][k] - sum1) / A[k][k];
		}
	}
}

tMatrix LUDecomposition(tMatrix &A, int n)
{
	datatype sum1;
	int row_max = 0;
	datatype max = 0;
	tVector boof = nullptr;
	//tVector Standings = nullptr; //в St-gs[i] записан номер строки, которую нужно поставить на место i-й
	//tVector b = nullptr;
	//tVector y = nullptr;
	tMatrix LU = nullptr;
	//tMatrix Solve = nullptr;
	//CreatingArray(Solve, n);
	CreatingArray(LU, n); //выделил память под итоговую матрицу
	CreatingVector(Standings, n);
	for (size_t k = 0; k < n; k++)
	{
		Standings[k] = k;
	}

	for (size_t k = 0; k < n; k++) //счетчик стр и стб, с которыми ведётся работа
	{
		row_max = k;
		max = abs(A[k][k]);
		for	(size_t i = k + 1; i < n; i++)
		{
			if (abs(A[i][k])>max)
			{
				row_max = i;
				max = abs(A[i][k]);
			}
		}//ищу максимальный элемент в столбце
		if (row_max != k)
		{
			boof = A[k];
			A[k] = A[row_max];
			A[row_max] = boof;

			boof = LU[k];
			LU[k] = LU[row_max];
			LU[row_max] = boof;

			Standings[k] = row_max;
			Standings[row_max] = k;
		}
		for (size_t j = k; j < n; j++) //работа с частями U
		{
			sum1 = 0;
			for (size_t i = 0; i < k; i++)
			{
				sum1 = sum1 + LU[k][i] * LU[i][j];
			}
			LU[k][j] = A[k][j] - sum1;
			if (LU[k][k] == 0)
			{
				cout << "LU decomposition is impossible ( det = 0 )" << endl;
				cout << "Returning nullptr..." << endl;
				return nullptr;
			}
		}
		for (size_t j = k+1; j < n; j++) //работа с частями L
		{
			sum1 = 0;
			for (size_t i = 0; i < k; i++)
			{
				sum1 = sum1 + LU[j][i] * LU[i][k];
			}
			LU[j][k] = (A[j][k] - sum1) / LU[k][k];
		}
	}
	//PrintingArray(LU, n); //это я проверял работоспособность на тестовых примерах
	//PrintingVector(Standings, n);

//	for (size_t k = 0; k < n; k++) //векторы в массиве В записаны по строкам
//	{
//		for (size_t i = 0; i < n; i++)
//		{
//			b[i] = B[k][Standings[i]];
//		}
//
//	}
//
//	DeletingVector(b);
//	DeletingVector(y);
//	DeletingVector(Standings);
//	DeletingArray(LU, n);
	return LU;
}

tVector LUSolve(tMatrix &LU, tVector &f, int n)
{
	size_t num = 0;
	datatype sum = 0;
	tVector b = nullptr;
	tVector y = nullptr;
	tVector Solve = nullptr;
	CreatingVector(b, n);
	CreatingVector(y, n);
	CreatingVector(Solve, n);
	for (size_t i = 0; i < n; i++)
	{
		num = Standings[i];
		b[i] = f[num];
	} //запись вектора правой части в нужном порядке

	for (size_t i = 0; i < n; i++)
	{
		sum = 0;
		for (size_t k = 0; k < i; k++)
		{
			sum = sum + LU[i][k] * y[k];
		}
		y[i] = b[i] - sum;
	}

	for (int i = n - 1; i > -1; i--)
	{
		sum = 0;
		for (size_t k = i + 1; k < n; k++)
		{
			sum = sum + LU[i][k] * Solve[k];
		}
		Solve[i] = (y[i] - sum) / LU[i][i];
	}
	//PrintingVector(Solve, n);
	DeletingVector(y);
	DeletingVector(b);
	return Solve;
}


int Dimension()
{
	int dim = 0;
	cout << "Enter the Dimension: ";
	cin >> dim;
	cout << endl;
	return dim;
}

void ClearRAM()
{
	DeletingVector(Standings);
}

void ADvsAS(dsys &fsys, tVector &X, int n)
{

		cout << "ADvsAS" << endl;
		//tVector X0 = VectCopy(X, n);

		char FileName[] = "ADvsAS.txt";
		ofstream OUT(FileName);
		OUT << "{";


		//tVector X0pr = nullptr;
		tVector AS = nullptr; //от слова Analytical Solve
		datatype tao1 = 15 * tao;
		for (size_t u = 0; u < 20; u++)
		{
			datatype max = 0;
			tVector X0 = VectCopy(X, n);
			cout << "Adams vs Analytical Solve ,  h = " << tao1 << endl;
			datatype tj = 0;


			tVector k1 = nullptr;
			tVector k2 = nullptr;
			tVector k3 = nullptr;
			tVector k4 = nullptr;
			tVector Xn = nullptr;

			tMatrix Xinitials = new tVector[4];//выделяю память под 4 вектора, нужных для начала итерационного процесса

			Xinitials[0] = VectCopy(X0, n);

				for (size_t i = 1; i < 4; i++) //по м-ду Р-К вычисляю 4 раза вектор y, записываю в массив Xinitials
				{
					k1 = fsys(X0, n);
					tVector k1pr = NumberVectorMultiplication(tao1 / 2, k1, n);
					tVector k1prpr = Summ(X0, k1pr, n);
					k2 = fsys(k1prpr, n);
					tVector k2pr = NumberVectorMultiplication(tao1 / 2, k2, n);
					tVector k2prpr = Summ(X0, k2pr, n);
					k3 = fsys(k2prpr, n);
					tVector k3pr = NumberVectorMultiplication(tao1, k3, n);
					tVector k3prpr = Summ(X0, k3pr, n);
					k4 = fsys(k3prpr, n);

					DeletingVector(k1pr);
					DeletingVector(k1prpr);
					DeletingVector(k2pr);
					DeletingVector(k2prpr);
					DeletingVector(k3pr);
					DeletingVector(k3prpr); // k1-k4

					k1pr = NumberVectorMultiplication(tao1 / 6, k1, n);
					k2pr = NumberVectorMultiplication(tao1 / 3, k2, n);
					k3pr = NumberVectorMultiplication(tao1 / 3, k3, n);
					tVector k4pr = NumberVectorMultiplication(tao1 / 6, k4, n);
					tVector K = Summ(k1pr, k2pr, k3pr, k4pr, n);
					DeletingVector(k1pr);
					DeletingVector(k2pr);
					DeletingVector(k3pr);
					DeletingVector(k4pr);
					//K = Summ(K, k3pr, n);
					//K = Summ(K, k4pr, n);
					//K = NumberVectorMultiplication(tao, K, n);
					Xn = Summ(X0, K, n);
					/*OUT << "{";
					for (size_t j = 0; j < n; j++)
					{
						OUT << Xn[j];
						if (j != n - 1)
						{
							OUT << ",";
						}
					}
					OUT << "}";
					OUT << ",";*/
					Xinitials[i] = VectCopy(Xn, n); //копирую в массив
					DeletingVector(X0);
					//X0 = nullptr;
					X0 = VectCopy(Xn, n);
					DeletingVector(K);
					DeletingVector(Xn);


					DeletingVector(k1);
					DeletingVector(k2);
					DeletingVector(k3);
					DeletingVector(k4);

					tj = tj + tao1;
				//	cout << "tj " << tj << endl;
				}
				DeletingVector(X0);  //allright
				while (tj < 5)
				{
					//cout << "tj " << tj << endl;
				//реализую формулы 1.18, имея 4 начальных вектора(f заменяю на k, нумерую слева направо)  ЗДЕСЬ ЧТО_ТО НЕПРАВИЛЬНО
				/*for (size_t i = 4; i < 500; i++)
				{*/
					//Xn = nullptr;
					k1 = fsys(Xinitials[3], n); //fn

					k2 = fsys(Xinitials[2], n); //fn-1

					k3 = fsys(Xinitials[1], n); //fn-2

					k4 = fsys(Xinitials[0], n); //fn-3

					tVector k1prime = NumberVectorMultiplication((55 * tao1) / 24, k1, n);
					tVector k2prime = NumberVectorMultiplication(-(59 * tao1) / 24, k2, n);
					tVector k3prime = NumberVectorMultiplication((37 * tao1) / 24, k3, n);
					tVector k4prime = NumberVectorMultiplication(-(9 * tao1) / 24, k4, n);


					tVector K = Summ(k1prime, k2prime, k3prime, k4prime, n);
					DeletingVector(k1prime);
					DeletingVector(k2prime);
					DeletingVector(k3prime);
					DeletingVector(k4prime);



					//tVector Yzero = Summ(Xinitials[3], K, n);
					Xn = Summ(Xinitials[3], K, n);
					DeletingVector(K);
                    AS = dsysAnalSolve(tj+tao1);


					datatype a = abs(Xn[0] - AS[0]);
					datatype b = abs(Xn[1] - AS[1]);
					/*cout << "tj " << tj << endl;
					cout << "Xn[0] " << Xn[0] << endl;
					cout << "AS[0] " << AS[0] << endl;
					cout << "Xn[1] " << Xn[1] << endl;
					cout << "AS[1] " << AS[1] << endl;
					cout << "a " << a << endl;
					cout << "b " << b << endl;*/
					if (a > b && max < a)
					{
						if (a > max) { max = a; }
					}
					else if (a <= b && max < b)
					{
						if (b > max) { max = b; }
					}
					DeletingVector(AS);


					//if (i % 499 != 0) OUT << ",";
					//

					DeletingVector(k1);
					DeletingVector(k2);
					DeletingVector(k3);
					DeletingVector(k4);




					DeletingVector(Xinitials[0]);
					Xinitials[0] = Xinitials[1];
					Xinitials[1] = Xinitials[2];
					Xinitials[2] = Xinitials[3];
					Xinitials[3] = VectCopy(Xn, n);
					DeletingVector(Xn);
					//DeletingVector(Yzero);
				    //}
					tj = tj + tao1;
			    }


			//OUT << "}";
			//cout << "End of Adams-Bushfort method(extrapolation)" << endl;
			//OUT.close();

			//DeletingVector(X0);
			cout << "max| y - U | = " << max << endl;
			cout << endl;
			OUT << "{";

			OUT << tao1;
			OUT << ",";
			OUT << max;

			OUT << "}";
			if (u != 19) OUT << ",";

			tao1 = tao1 * 0.7;

			for (size_t i = 0; i < 4; i++)
			{
				delete[] Xinitials[i];
			}
			delete[] Xinitials;


		}
		OUT << "}";
		cout << "End of htaoanal method." << endl;
		OUT.close();
}





void PKvsAS(dsys &fsys, tVector &X, int n)
{

	cout << "PKvsAS" << endl;
	//tVector X0 = VectCopy(X, n);

	char FileName[] = "PKvsAS.txt";
	ofstream OUT(FileName);
	OUT << "{";


	//tVector X0pr = nullptr;
	tVector AS = nullptr; //от слова Analytical Solve
	datatype tao1 = 15 * tao;
	for (size_t u = 0; u < 20; u++)
	{
		datatype max = 0;
		tVector X0 = VectCopy(X, n);
		cout << "Adams vs Analytical Solve ,  h = " << tao1 << endl;
		datatype tj = 0;


		tVector k1 = nullptr;
		tVector k2 = nullptr;
		tVector k3 = nullptr;
		tVector k4 = nullptr;
		tVector Xn = nullptr;

		tMatrix Xinitials = new tVector[4];//выделяю память под 4 вектора, нужных для начала итерационного процесса

		Xinitials[0] = VectCopy(X0, n);

		for (size_t i = 1; i < 4; i++) //по м-ду Р-К вычисляю 4 раза вектор y, записываю в массив Xinitials
	{
		k1 = fsys(X0, n);
		tVector k1pr = NumberVectorMultiplication(tao1 / 2, k1, n);
		tVector k1prpr = Summ(X0, k1pr, n);
		k2 = fsys(k1prpr, n);
		tVector k2pr = NumberVectorMultiplication(tao1 / 2, k2, n);
		tVector k2prpr = Summ(X0, k2pr, n);
		k3 = fsys(k2prpr, n);
		tVector k3pr = NumberVectorMultiplication(tao1, k3, n);
		tVector k3prpr = Summ(X0, k3pr, n);
		k4 = fsys(k3prpr, n);

		DeletingVector(k1pr);
		DeletingVector(k1prpr);
		DeletingVector(k2pr);
		DeletingVector(k2prpr);
		DeletingVector(k3pr);
		DeletingVector(k3prpr);

		k1pr = NumberVectorMultiplication(tao1 / 6, k1, n);
		k2pr = NumberVectorMultiplication(tao1 / 3, k2, n);
		k3pr = NumberVectorMultiplication(tao1 / 3, k3, n);
		tVector k4pr = NumberVectorMultiplication(tao1 / 6, k4, n);
		DeletingVector(k1);
		DeletingVector(k2);
		DeletingVector(k3);
		DeletingVector(k4);
		tVector K = Summ(k1pr, k2pr, k3pr, k4pr, n);
		//K = Summ(K, k3pr, n);
		//K = Summ(K, k4pr, n);
		//K = NumberVectorMultiplication(tao, K, n);
		Xn = Summ(X0, K, n);
		DeletingVector(K);
		DeletingVector(X0);
		Xinitials[i] = VectCopy(Xn, n); //копирую в массив
		//X0 = nullptr;
		X0 = VectCopy(Xn, n);
		DeletingVector(Xn);
		DeletingVector(k1pr);
		DeletingVector(k2pr);
		DeletingVector(k3pr);
		DeletingVector(k4pr);
		tj = tj + tao1;
	}
		DeletingVector(X0);  //allright
		while (tj < 5)
		{
						Xn = nullptr;
						k1 = fsys(Xinitials[3], n); //fn
						//PrintingVector(k1);
						k2 = fsys(Xinitials[2], n); //fn-1
						//PrintingVector(k2);
						k3 = fsys(Xinitials[1], n); //fn-2
						//PrintingVector(k3);
						k4 = fsys(Xinitials[0], n); //fn-3
						//PrintingVector(k4);
						tVector k1prime = NumberVectorMultiplication(55 * tao1 / 24, k1, n);
						tVector k2prime = NumberVectorMultiplication(-59 * tao1 / 24, k2, n);
						tVector k3prime = NumberVectorMultiplication(37 * tao1 / 24, k3, n);
						tVector k4prime = NumberVectorMultiplication(-9 * tao1 / 24, k4, n);
						/*DeletingVector(k1);
						DeletingVector(k2);
						DeletingVector(k3);*/
						DeletingVector(k4);
						tVector K = Summ(k1prime, k2prime, k3prime, k4prime, n);
						DeletingVector(k1prime);
						DeletingVector(k2prime);
						DeletingVector(k3prime);
						DeletingVector(k4prime);
						tVector Yzero = Summ(Xinitials[3], K, n);
						DeletingVector(K);
						tVector Fzero = fsys(Yzero, n);
						DeletingVector(Yzero);
						k1prime = NumberVectorMultiplication(9 * tao1 / 24, Fzero, n);
						DeletingVector(Fzero);
						k2prime = NumberVectorMultiplication(19 * tao1 / 24, k1, n);
						k3prime = NumberVectorMultiplication(-5 * tao1 / 24, k2, n);
						k4prime = NumberVectorMultiplication(1 * tao1 / 24, k3, n);
						DeletingVector(k1);
						DeletingVector(k2);
						DeletingVector(k3);
						K = Summ(k1prime, k2prime, k3prime, k4prime, n);
						DeletingVector(k1prime);
						DeletingVector(k2prime);
						DeletingVector(k3prime);
						DeletingVector(k4prime);
					    Xn = Summ(Xinitials[3], K, n);


						DeletingVector(K);

						AS = dsysAnalSolve(tj + tao1);


						datatype a = abs(Xn[0] - AS[0]);
						datatype b = abs(Xn[1] - AS[1]);
						/*cout << "tj " << tj << endl;
						cout << "Xn[0] " << Xn[0] << endl;
						cout << "AS[0] " << AS[0] << endl;
						cout << "Xn[1] " << Xn[1] << endl;
						cout << "AS[1] " << AS[1] << endl;
						cout << "a " << a << endl;
						cout << "b " << b << endl;*/
						if (a > b && max < a)
						{
							if (a > max) { max = a; }
						}
						else if (a <= b && max < b)
						{
							if (b > max) { max = b; }
						}
						DeletingVector(AS);

						tj = tj + tao1;
						DeletingVector(Xinitials[0]);
						Xinitials[0] = Xinitials[1];
						Xinitials[1] = Xinitials[2];
						Xinitials[2] = Xinitials[3];
						Xinitials[3] = VectCopy(Xn, n);
						DeletingVector(Xn);

		}


		//OUT << "}";
		//cout << "End of Adams-Bushfort method(extrapolation)" << endl;
		//OUT.close();

		//DeletingVector(X0);
		cout << "max| y - U | = " << max << endl;
		cout << endl;
		OUT << "{";

		OUT << tao1;
		OUT << ",";
		OUT << max;

		OUT << "}";
		if (u != 19) OUT << ",";

		tao1 = tao1 * 0.7;

		for (size_t i = 0; i < 4; i++)
		{
			delete[] Xinitials[i];
		}
		delete[] Xinitials;


	}
	OUT << "}";
	cout << "End of htaoanal method." << endl;
	OUT.close();

}


void Runge_Kutt2(dsys &fsys, tVector &X, int n)
{
	cout << "Runge-Kutt Method" << endl;
	tVector X0 = VectCopy(X, n);
	tVector X0pr = nullptr;
	unsigned int key = 0;
	datatype tao1 = tao;
	datatype tj = 0;

	char FileName[] = "Data_Runge-Kutt2.txt";
	ofstream OUT(FileName);
	OUT << "{";

	//OUT << "{" << X0[0] << "," << X0[1] << "}"<<",";

	OUT << "{";
	for (size_t i = 0; i < n; i++)
	{
		OUT << X[i];
		if (i != n - 1)
		{
			OUT << ",";
		}
	}
	OUT << "},";

	cout << "Control the step?" << endl << "YOUR CHOISE(0-no) : ";
	cin >> key;
	cout << endl;
	switch (key)
	{
	case 0:
	{
		tVector k1 = nullptr;
		tVector k2 = nullptr;

		tVector Xn = nullptr;
		/*CreatingVector(k1, n);
		CreatingVector(k2, n);
		CreatingVector(k3, n);
		CreatingVector(k4, n);*/



		for (size_t i = 1; i < 100001; i++)
		{
			k1 = fsys(X0, n);
			tVector k1pr = NumberVectorMultiplication(tao / 2, k1, n);
			tVector k1prpr = Summ(X0, k1pr, n);
			k2 = fsys(k1prpr, n);


			DeletingVector(k1pr);
			DeletingVector(k1prpr);



			tVector k2pr = NumberVectorMultiplication(tao , k2, n);


			//K = Summ(K, k3pr, n);
			//K = Summ(K, k4pr, n);
			//K = NumberVectorMultiplication(tao, K, n);
			Xn = Summ(X0, k2pr, n);
			//OUT << "{" << Xn[0] << "," << Xn[1] << "}";
			OUT << "{";
			for (size_t i = 0; i < n; i++)
			{
				OUT << Xn[i];
				if (i != n - 1)
				{
					OUT << ",";
				}
			}
			OUT << "}";
			if (i % 100000 != 0) OUT << ",";
			DeletingVector(X0);
			X0 = nullptr;
			X0 = VectCopy(Xn, n);
			//DeletingVector(K);
			DeletingVector(Xn);
			//DeletingVector(k1pr);
			DeletingVector(k2pr);


		}
		OUT << "}";
		cout << "End of RK2 method." << endl;
		OUT.close();


		DeletingVector(X0);
		break;
	}
	case 1:
	{
		tVector k1 = nullptr;
		tVector k2 = nullptr;
		tVector k3 = nullptr;
		tVector k4 = nullptr;
		tVector Xn = nullptr;
		tVector Xn1 = nullptr;
		tVector DX = nullptr;
		datatype norm = 0;

		/*CreatingVector(k1, n);
		CreatingVector(k2, n);
		CreatingVector(k3, n);
		CreatingVector(k4, n);*/


	beginning:
		while (tj < 1)
		{
			k1 = fsys(X0, n);
			tVector k1pr = NumberVectorMultiplication(tao1 / 2, k1, n);
			tVector k1prpr = Summ(X0, k1pr, n);
			k2 = fsys(k1prpr, n);
			tVector k2pr = NumberVectorMultiplication(tao1 / 2, k2, n);
			tVector k2prpr = Summ(X0, k2pr, n);
			k3 = fsys(k2prpr, n);
			tVector k3pr = NumberVectorMultiplication(tao1, k3, n);
			tVector k3prpr = Summ(X0, k3pr, n);
			k4 = fsys(k3prpr, n);

			DeletingVector(k1pr);
			DeletingVector(k1prpr);
			DeletingVector(k2pr);
			DeletingVector(k2prpr);
			DeletingVector(k3pr);
			DeletingVector(k3prpr);

			k1pr = NumberVectorMultiplication(tao1 / 6, k1, n);
			k2pr = NumberVectorMultiplication(tao1 / 3, k2, n);
			k3pr = NumberVectorMultiplication(tao1 / 3, k3, n);
			tVector k4pr = NumberVectorMultiplication(tao1 / 6, k4, n);
			tVector K = Summ(k1pr, k2pr, k3pr, k4pr, n);
			DeletingVector(k1);
			DeletingVector(k2);
			DeletingVector(k3);
			DeletingVector(k4);

			Xn = Summ(X0, K, n); //вычислено приближение при шаге tao
			//DeletingVector(X0);
			//X0 = nullptr;
			//X0 = VectCopy(Xn,n);
			DeletingVector(K);
			//DeletingVector(Xn);
			DeletingVector(k1pr);
			DeletingVector(k2pr);
			DeletingVector(k3pr);
			DeletingVector(k4pr);



			X0pr = VectCopy(X0, n);
			tao1 = tao1 / 2; //подробили шаг
			for (size_t i = 0; i < 2; i++)
			{
				k1 = fsys(X0pr, n);
				tVector k1pr = NumberVectorMultiplication(tao1 / 2, k1, n);
				tVector k1prpr = Summ(X0pr, k1pr, n);
				k2 = fsys(k1prpr, n);
				tVector k2pr = NumberVectorMultiplication(tao1 / 2, k2, n);
				tVector k2prpr = Summ(X0pr, k2pr, n);
				k3 = fsys(k2prpr, n);
				tVector k3pr = NumberVectorMultiplication(tao1, k3, n);
				tVector k3prpr = Summ(X0pr, k3pr, n);
				k4 = fsys(k3prpr, n);

				DeletingVector(k1pr);
				DeletingVector(k1prpr);
				DeletingVector(k2pr);
				DeletingVector(k2prpr);
				DeletingVector(k3pr);
				DeletingVector(k3prpr);

				k1pr = NumberVectorMultiplication(tao1 / 6, k1, n);
				k2pr = NumberVectorMultiplication(tao1 / 3, k2, n);
				k3pr = NumberVectorMultiplication(tao1 / 3, k3, n);
				tVector k4pr = NumberVectorMultiplication(tao1 / 6, k4, n);
				tVector K = Summ(k1pr, k2pr, k3pr, k4pr, n);
				DeletingVector(k1);
				DeletingVector(k2);
				DeletingVector(k3);
				DeletingVector(k4);

				Xn1 = Summ(X0pr, K, n);
				DeletingVector(X0pr);
				X0pr = nullptr;
				X0pr = VectCopy(Xn1, n);    //вычислено полуприближение при шаге tao/2 при i=0
				//и нормальное приближение при i=1
				DeletingVector(K);
				DeletingVector(Xn1);
				DeletingVector(k1pr);
				DeletingVector(k2pr);
				DeletingVector(k3pr);
				DeletingVector(k4pr);

			} //имеем 2 приближения: 1-но шажное и 2-х шажное (НЕ ЗАБЫТЬ УДАЛИТЬ!!)


			DX = Difference(X0pr, Xn, n);
			norm = NormOfVector(DX, n);
			DeletingVector(Xn);
			norm = norm / 15;   //страница 16 методички
			DeletingVector(DX);
			if (norm < 1.e-6 && norm > 1.e-10)  // Делаю как сам вижу, в методе непонятно написано
			{
				tao1 = 2 * tao1; //возвращаем исходный шаг
				OUT << "{";
				for (size_t i = 0; i < n; i++)
				{
					OUT << X0pr[i];
					if (i != n - 1)
					{
						OUT << ",";
					}
				}
				OUT << "},";
				tj = tj + tao1;
				DeletingVector(X0);
				X0 = VectCopy(X0pr, n);
				DeletingVector(X0pr);
				goto beginning;
			}
			else
				if (norm < 1.e-10)
				{
					tao1 = 4 * tao1; //удлиняем шаг в 2 раза
					DeletingVector(X0pr);
					goto beginning;
				}
				else
				{
					DeletingVector(X0pr);
					goto beginning;
				}

		}
		OUT << "{";
		for (size_t i = 0; i < n; i++) //дублирую запись последней точки в массив значений
		{							// потому что почему бы и нет, чтоб не мучиться
			OUT << X0[i];
			if (i != n - 1)
			{
				OUT << ",";
			}
		}
		OUT << "}";
		OUT << "}";
		cout << "End of RK method" << endl;
		OUT.close();


		DeletingVector(X0);
		break;
	}
	}

}

//datatype Integrate(funct &f, datatype a, datatype b)
//{
//
//	return (( (b - a) / 6) * (f(a) + 4 * f((a + b) / 2) + f(b)));
//}

datatype Integrate(funct &f, datatype a, datatype b)
{
	cout << "(a + b) / 2   " << (a + b) / 2 << endl;
	cout << "f((a + b) / 2)    "<< f((a + b) / 2) << endl;
	cout << "(b - a)   " << (b - a) << endl;
	cout<<"itogo   " <<f((a + b) / 2)*(b - a) <<endl;
	return (f((a + b) / 2)*(b - a));
}

//void HeatTransfer(datatype sigma)
//{
//	char FileName[] = "x.txt";
//	char FileName2[] = "t.txt";
//	char FileName3[] = "u.txt";
//	ofstream OUTx(FileName);
//	ofstream OUTt(FileName2);
//	ofstream OUTu(FileName3);
//
//	funct f = Kx; //для вычисления коэффициентов a
//	funct P = P3; // граничное условие потока варианта 15
//	//предполагаю решать задачу на временном промежутке от 0 до 0.7с  и с длиной стержня 1
//	datatype dt = 0.0007;
//	datatype dh = 0.05;
//	int n = 21; //всего узлов на отрезке от 0 до 1 с шагом h  0---20
//	int m = 1001; // временных слоёв на указанной мною сетке   0---1000
//
//	datatype current_t = 0; //переменная текущего времени
//	
//	OUTu << "{" << endl;
//	OUTu << "{";
//	tVector solve = new datatype[n]; //массив точек температур в узлах сетки на текущем временном слое
//	for (int i = 0; i < n-1; i++)
//	{
//		solve[i] = 0.04;  // начальные условия варианта 15
//		OUTu << 0.04 << ",";
//	}
//	solve[n - 1] = 0.04;
//	OUTu << 0.04 << "}," << endl; //запись первой строчки температур в файл
//
//	OUTx << "{";
//	tVector x = new datatype[n]; //массив пространственной сетки 0---20
//	for (int i = 0; i < n-1; i++)
//	{
//		x[i] = i * dh;
//		OUTx << x[i] << ",";
//	}
//	x[n - 1] = (n - 1)*dh;
//	OUTx << x[n - 1] << "}"; //запись простанственной сетки в файл
//	OUTx.close();
//
//	OUTt << "{";
//	tVector t = new datatype[m];  // массив временной сетки 0---1000
//	for (int i = 0; i < m-1; i++)
//	{
//		t[i] = i * dt;
//		OUTt << t[i] << ",";
//	}
//	t[m - 1] = (m - 1)*dt;
//	OUTt << t[m - 1] << "}"; //запись временнОй сетки в файл
//	OUTt.close();
//
//	tVector a = new datatype[n - 1]; //массив маленьких а для вычисления А, В, С (n-1 есть число промежутков интегрирования)
//	for (int i = 0; i < n - 1; i++)
//	{
//		a[i] =dh / Integrate(f, x[i], x[i + 1]); //обрабатываются все отрезочки
//	}
//	cout << "a:" << endl;
//	PrintingVector(a, n-1);
//	cin.get();
//
//	//в терминах методички N = 20
//
//	tVector A = new datatype[n - 2];  //0---18
//	for (int i = 0; i < n - 2; i++)
//	{
//		A[i] = sigma * a[i] / dh;
//	}
//	cout << "A:" << endl;
//	PrintingVector(A, n-2);
//	cin.get();
//
//	tVector B = new datatype[n - 2];  //0---18
//	for (int i = 0; i < n - 2; i++)
//	{
//		B[i] = sigma * a[i + 1] / dh;
//	}
//	PrintingVector(B, n - 2);
//	cin.get();
//
//	tVector C = new datatype[n - 2];  //0---18
//	for (int i = 0; i < n - 2; i++)
//	{
//		C[i] = A[i] + B[i] + dh/dt;
//	}
//	PrintingVector(C, n - 2);
//	cin.get();
//
//	datatype kappa = (sigma * a[n - 2] / dh) / (dh / (2 * dt) + sigma * a[n - 2] / dh);
//	cout << kappa;
//	cin.get();
//
//	while (current_t <= 0.7)
//	{
//	// задал все статичные параметры, осталось задать F и mu. Приступим:
//
//		tVector F = new datatype[n - 2];
//		for (int i = 0; i < n - 2; i++)  //старался сделать без ошибки, но мб она тут есть
//		{
//			F[i] = (dh / dt) * solve[i + 1] + (1 - sigma) * (a[i + 1] * (solve[i + 2] - solve[i + 1]) - a[i] * (solve[i + 1] - solve[i])) / dh;
//		}
//		cout << "F:" << endl;
//		PrintingVector(F, n - 2);
//		cin.get();
//		datatype mu = (solve[n - 1] * dh / (2 * dt) + sigma * P(current_t + dt) + (1 - sigma) * (P(current_t) - a[n - 2] * (solve[n - 1] - solve[n - 2]) / h))
//		/ (dh / (2 * dt) + sigma * a[n - 2] / dh);
//		cout << "mu "<<mu << endl;
//		cin.get();
//		//теперь коэффициенты трехдиагональной матрицы:
//
//		tVector ldiag = new datatype[n - 3];  //0---17
//		for (int i = 0; i < n - 3; i++)
//		{
//			ldiag[i] = A[i + 1];
//		}
//
//		tVector diag = new datatype[n - 2];  //записываю как последний ублюдок
//		for (int i = 0; i < n - 3; i++)      //в терминах методички с минусами
//		{
//			diag[i] = -C[i];
//		}
//		diag[n - 3] = -C[n - 3] + kappa * B[n - 3];
//
//		tVector udiag = new datatype[n - 3];
//		for (int i = 0; i < n - 3; i++)
//		{
//			udiag[i] = B[i];
//		}
//
//
//
//		//СМОТРИ ЕБАТЬ НАХУЙ ПОХОДУ ОДОГО КОЭФФИЦИЕНТА ДЛЯ А Б  С И Ф НЕ ХВАТАЕТ ОЛООООО
//		//СМОТРИ ЕБАТЬ НАХУЙ ПОХОДУ ОДОГО КОЭФФИЦИЕНТА ДЛЯ А Б  С И Ф НЕ ХВАТАЕТ ОЛООООО
//		//СМОТРИ ЕБАТЬ НАХУЙ ПОХОДУ ОДОГО КОЭФФИЦИЕНТА ДЛЯ А Б  С И Ф НЕ ХВАТАЕТ ОЛООООО
//		//СМОТРИ ЕБАТЬ НАХУЙ ПОХОДУ ОДОГО КОЭФФИЦИЕНТА ДЛЯ А Б  С И Ф НЕ ХВАТАЕТ ОЛООООО
//		//СМОТРИ ЕБАТЬ НАХУЙ ПОХОДУ ОДОГО КОЭФФИЦИЕНТА ДЛЯ А Б  С И Ф НЕ ХВАТАЕТ ОЛООООО
//		//СМОТРИ ЕБАТЬ НАХУЙ ПОХОДУ ОДОГО КОЭФФИЦИЕНТА ДЛЯ А Б  С И Ф НЕ ХВАТАЕТ ОЛООООО
//		//СМОТРИ ЕБАТЬ НАХУЙ ПОХОДУ ОДОГО КОЭФФИЦИЕНТА ДЛЯ А Б  С И Ф НЕ ХВАТАЕТ ОЛООООО
//		//СМОТРИ ЕБАТЬ НАХУЙ ПОХОДУ ОДОГО КОЭФФИЦИЕНТА ДЛЯ А Б  С И Ф НЕ ХВАТАЕТ ОЛООООО
//		//СМОТРИ ЕБАТЬ НАХУЙ ПОХОДУ ОДОГО КОЭФФИЦИЕНТА ДЛЯ А Б  С И Ф НЕ ХВАТАЕТ ОЛООООО
//		//СМОТРИ ЕБАТЬ НАХУЙ ПОХОДУ ОДОГО КОЭФФИЦИЕНТА ДЛЯ А Б  С И Ф НЕ ХВАТАЕТ ОЛООООО
//
//
//		//Осталось записать только вектор правой части:
//
//		tVector b = new datatype[n - 2]; // 0---18
//		b[0] = -F[0] - 0.04 * A[0];
//		for (int i = 1; i < n - 4; i++)
//		{
//			b[i] = -F[i];
//		}
//		b[n - 3] = -F[n - 3] - mu * B[n - 3];
//
//		////////////////////////// Все готово для нахождения решения на новом слое.
//		////////////////////////// Не уверен только в правильности индексации, но все должно быть ок.
//		////////////////////////// Сопоставлял индексы в уме. Ты если будешь править ошибки, то лучше распиши на бумаге
//		////////////////////////// В один момент щёлнет, как надо правильно индексировать.
//		////////////////////////// я сосал меня ебали
//		tVector solveboof = Shuttle(ldiag, diag, udiag, b, n - 2);
//		tVector solvenew = new datatype[n]; // вектор решения на новом временном слое
//
//		solvenew[0] = 0.04;
//		OUTu << "{" << solvenew[0] << ",";
//		for (int i = 1; i < n - 1; i++)
//		{
//			solvenew[i] = solveboof[i - 1];
//			OUTu << solvenew[i] << ",";
//		}
//		solvenew[n - 1] = kappa * solvenew[n - 2] + mu;
//		OUTu << solvenew[n - 1] << "}," << endl;
////		for (int i = 0; i < n ; i++)
////		{
////			cout << solvenew[i] << endl;
////		}
//		current_t = current_t + dt;
//		cout<<current_t<<endl;
//
//		DeletingVector(F);
//		DeletingVector(ldiag);
//		DeletingVector(diag);
//		DeletingVector(udiag);
//		DeletingVector(b);
//		DeletingVector(solveboof);
//		DeletingVector(solve);
//		solve = VectCopy(solvenew, n);
//		DeletingVector(solvenew);
//		// добавь запись этого вектора в файл. Файл не такой, как аньше, а другой.
//		// Я рисовал на лабе, как должен выглядеть этот файл. Если что, напиши, я напомню.
//
//	}
//
//	DeletingVector(x);
//	DeletingVector(t);
//	DeletingVector(a);
//	DeletingVector(A);
//	DeletingVector(B);
//	DeletingVector(C);
//	//PrintingVector(solve, n);
//	DeletingVector(solve);
//
////	// задал все статичные параметры, осталось задать F и mu. Приступим:
////
////	tVector F = new datatype[n - 2];
////	for (int i = 0; i < n - 2; i++)  //старался сделать без ошибки, но мб она тут есть
////	{
////		F[i] = (dh / dt) * solve[i + 1] + (1 - sigma) * (a[i + 1] * (solve[i + 2] - solve[i + 1]) - a[i] * (solve[i + 1] - solve[i])) / dh;
////	}
////
////	datatype mu = (solve[n - 1] * dh / (2 * dt) + sigma * P(current_t + dt) + (1 - sigma) * (P(current_t) - a[n - 2] * (solve[n - 1] - solve[n - 2]) / h))
////	/ (dh / (2 * dt) + sigma * a[n - 2] / dh);
////
////	//теперь коэффициенты трехдиагональной матрицы:
////
////	tVector ldiag = new datatype[n - 3];  //0---17
////	for (int i = 0; i < n - 3; i++)
////	{
////		ldiag[i] = A[i + 1];
////	}
////
////	tVector diag = new datatype[n - 2];  //записываю как последний ублюдок
////	for (int i = 0; i < n - 3; i++)      //в терминах методички с минусами
////	{
////		diag[i] = -C[i];
////	}
////	diag[n - 3] = -C[n - 3] + kappa * B[n - 3];
////
////	tVector udiag = new datatype[n - 3];
////	for (int i = 0; i < n - 3; i++)
////	{
////		udiag[i] = B[i];
////	}
////
////	//Осталось записать только вектор правой части:
////
////	tVector b = new datatype[n - 2]; // 0---18
////	b[0] = -F[0] - 0.04 * A[0];
////	for (int i = 1; i < n - 4; i++)
////	{
////		b[i] = -F[i];
////	}
////	b[n - 3] = -F[n - 3] - mu * B[n - 3];
////
////	////////////////////////// Все готово для нахождения решения на новом слое.
////	////////////////////////// Не уверен только в правильности индексации, но все должно быть ок.
////	////////////////////////// Сопоставлял индексы в уме. Ты если будешь править ошибки, то лучше распиши на бумаге
////	////////////////////////// В один момент щёлнет, как надо правильно индексировать.
////	tVector solveboof = Shuttle(ldiag, diag, udiag, b, n - 2);
//////	for (int i = 0; i < n - 2; i++)
//////	{
//////		cout << solvenew[i] << endl;
//////	}
////	solvenew[0] = 0.04;
////	for (int i = 1; i < n - 1; i++)
////	{
////		solvenew[i] = solveboof[i - 1];
////	}
////	solvenew[n - 1] = kappa * solvenew[n - 2] + mu;
////	for (int i = 0; i < n ; i++)
////	{
////		cout << solvenew[i] << endl;
////	}
//}


void HeatTransfer(datatype sigma)
{
	char FileName[] = "x.txt";
	char FileName2[] = "t.txt";
	char FileName3[] = "u.txt";
	ofstream OUTx(FileName);
	ofstream OUTt(FileName2);
	ofstream OUTu(FileName3);

	funct f = Kx; //для вычисления коэффициентов a
	funct P = P3; // граничное условие потока варианта 15
	//предполагаю решать задачу на временном промежутке от 0 до 0.7с  и с длиной стержня 1
	datatype dt = 0.0001;
	datatype dh = 0.05;
	int n = 21; //всего узлов на отрезке от 0 до 1 с шагом h  0---20
	int m = 7001; // временных слоёв на указанной мною сетке   0---1000

	datatype current_t = 0; //переменная текущего времени
	OUTu << "{" << endl;
	OUTu << "{";
	tVector solve = new datatype[n]; //массив точек температур в узлах сетки на текущем временном слое 0---20
	for (int i = 0; i < n - 1; i++)
	{
		solve[i] = 0.04;  // начальные условия варианта 15
		OUTu << 0.04 << ",";
	}
	solve[n - 1] = 0.04;
	OUTu << 0.04 << "}," << endl; //запись первой строчки температур в файл

	OUTx << "{";
	tVector x = new datatype[n]; //массив пространственной сетки 0---20
	for (int i = 0; i < n - 1; i++)
	{
		x[i] = i * dh;
		OUTx << x[i] << ",";
	}
	x[n - 1] = (n - 1)*dh;
	OUTx << x[n - 1] << "}"; //запись простанственной сетки в файл
	OUTx.close();

	OUTt << "{";
	tVector t = new datatype[m];  // массив временной сетки 0---1000
	for (int i = 0; i < m - 1; i++)
	{
		t[i] = i * dt;
		OUTt << t[i] << ",";
	}
	t[m - 1] = (m - 1)*dt;
	OUTt << t[m - 1] << "}"; //запись временнОй сетки в файл
	OUTt.close();

	tVector a = new datatype[n - 1]; //массив маленьких а для вычисления А, В, С (n-1 есть число промежутков интегрирования)1---20(0---19)
	for (int i = 0; i < n - 1; i++)
	{
		a[i] = dh / Integrate(f, x[i], x[i + 1]); //обрабатываются все отрезочки
	}
//	cout << "a:" << endl;
//	PrintingVector(a, n - 1);
//	cin.get();

	//в терминах методички N = 20

	tVector A = new datatype[n - 2];  //1---19(0---18)
	for (int i = 0; i < n - 2; i++)
	{
		A[i] = sigma * a[i] / dh;
	}
//	cout << "A:" << endl;
//	PrintingVector(A, n - 2);
//	cin.get();

	tVector B = new datatype[n - 2];  //1---19(0---18)
	for (int i = 0; i < n - 2; i++)
	{
		B[i] = sigma * a[i + 1] / dh;
	}
//	PrintingVector(B, n - 2);
//	cin.get();

	tVector C = new datatype[n - 2];  //1---19(0---18)
	for (int i = 0; i < n - 2; i++)
	{
		C[i] = A[i] + B[i] + dh / dt;
	}
//	PrintingVector(C, n - 2);
//	cin.get();

	datatype kappa = (sigma * a[n - 2] / dh) / (dh / (2 * dt) + sigma * a[n - 2] / dh);
//	cout << kappa;
//	cin.get();

	while (current_t < 0.6999)
	{
		// задал все статичные параметры, осталось задать F и mu. Приступим:

		tVector F = new datatype[n - 2]; //1---19(0---18)
		for (int i = 0; i < n - 2; i++)  //старался сделать без ошибки, но мб она тут есть
		{
			F[i] = (dh / dt) * solve[i + 1] + (1 - sigma) * (a[i + 1] * (solve[i + 2] - solve[i + 1]) - a[i] * (solve[i + 1] - solve[i])) / dh;
		}

		datatype mu = (solve[n - 1] * dh / (2 * dt) + sigma * P(current_t + dt) + (1 - sigma) * (P(current_t) - a[n - 2] * (solve[n - 1] - solve[n - 2]) / h))
			/ (dh / (2 * dt) + sigma * a[n - 2] / dh);


		//теперь коэффициенты трехдиагональной матрицы:

		tVector ldiag = new datatype[n - 3];  //0---17 (1--18)
		for (int i = 0; i < n - 3; i++)
		{
			ldiag[i] = A[i + 1];
		}

		tVector diag = new datatype[n - 2];  //0---18 (1--19)//записываю как последний ублюдок
		for (int i = 0; i < n - 3; i++)      //в терминах методички с минусами(похуй)
		{
			diag[i] = -C[i];
		}
		diag[n - 3] = -C[n - 3] + kappa * B[n - 3];

		tVector udiag = new datatype[n - 3];//0---17 (1--18)
		for (int i = 0; i < n - 3; i++)
		{
			udiag[i] = B[i];
		}


		//Осталось записать только вектор правой части:

		tVector b = new datatype[n - 2]; // 0---18
		b[0] = -F[0] - 0.04 * A[0];
		for (int i = 1; i < n - 3; i++)
		{
			b[i] = -F[i];
		}
		b[n - 3] = -F[n - 3] - mu * B[n - 3];

		//PrintingVector(b, n - 2);
		////////////////////////// Все готово для нахождения решения на новом слое.
		////////////////////////// Не уверен только в правильности индексации, но все должно быть ок.
		////////////////////////// Сопоставлял индексы в уме. Ты если будешь править ошибки, то лучше распиши на бумаге
		////////////////////////// В один момент щёлнет, как надо правильно индексировать.
		////////////////////////// я сосал меня ебали
		tVector solveboof = Shuttle(ldiag, diag, udiag, b, n - 2);
		tVector solvenew = new datatype[n]; // вектор решения на новом временном слое


		solvenew[0] = 0.04;
		OUTu << "{" << solvenew[0] << ",";
		for (int i = 1; i < n - 1; i++)
		{
			solvenew[i] = solveboof[i - 1];
			OUTu << solvenew[i] << ",";
		}
		solvenew[n - 1] = kappa * solvenew[n - 2] + mu;
		OUTu << solvenew[n - 1] << "}," << endl;
		current_t = current_t + dt;
		cout << current_t << endl;

		DeletingVector(F);
		DeletingVector(ldiag);
		DeletingVector(diag);
		DeletingVector(udiag);
		DeletingVector(b);
		DeletingVector(solveboof);
		DeletingVector(solve);
		solve = VectCopy(solvenew, n);
		DeletingVector(solvenew);
		// добавь запись этого вектора в файл. Файл не такой, как pаньше, а другой.
		// Я рисовал на лабе, как должен выглядеть этот файл. Если что, напиши, я напомню.

	}
	tVector F = new datatype[n - 2]; //1---19(0---18)
	for (int i = 0; i < n - 2; i++)  //старался сделать без ошибки, но мб она тут есть
	{
		F[i] = (dh / dt) * solve[i + 1] + (1 - sigma) * (a[i + 1] * (solve[i + 2] - solve[i + 1]) - a[i] * (solve[i + 1] - solve[i])) / dh;
	}

	datatype mu = (solve[n - 1] * dh / (2 * dt) + sigma * P(current_t + dt) + (1 - sigma) * (P(current_t) - a[n - 2] * (solve[n - 1] - solve[n - 2]) / h))
		/ (dh / (2 * dt) + sigma * a[n - 2] / dh);


	//теперь коэффициенты трехдиагональной матрицы:

	tVector ldiag = new datatype[n - 3];  //0---17 (1--18)
	for (int i = 0; i < n - 3; i++)
	{
		ldiag[i] = A[i + 1];
	}

	tVector diag = new datatype[n - 2];  //0---18 (1--19)//записываю как последний ублюдок
	for (int i = 0; i < n - 3; i++)      //в терминах методички с минусами(похуй)
	{
		diag[i] = -C[i];
	}
	diag[n - 3] = -C[n - 3] + kappa * B[n - 3];

	tVector udiag = new datatype[n - 3];//0---17 (1--18)
	for (int i = 0; i < n - 3; i++)
	{
		udiag[i] = B[i];
	}


	//Осталось записать только вектор правой части:

	tVector b = new datatype[n - 2]; // 0---18
	b[0] = -F[0] - 0.04 * A[0];
	for (int i = 1; i < n - 3; i++)
	{
		b[i] = -F[i];
	}
	b[n - 3] = -F[n - 3] - mu * B[n - 3];

	//PrintingVector(b, n - 2);
	////////////////////////// Все готово для нахождения решения на новом слое.
	////////////////////////// Не уверен только в правильности индексации, но все должно быть ок.
	////////////////////////// Сопоставлял индексы в уме. Ты если будешь править ошибки, то лучше распиши на бумаге
	////////////////////////// В один момент щёлнет, как надо правильно индексировать.
	////////////////////////// я сосал меня ебали
	tVector solveboof = Shuttle(ldiag, diag, udiag, b, n - 2);
	tVector solvenew = new datatype[n]; // вектор решения на новом временном слое


	solvenew[0] = 0.04;
	OUTu << "{" << solvenew[0] << ",";
	for (int i = 1; i < n - 1; i++)
	{
		solvenew[i] = solveboof[i - 1];
		OUTu << solvenew[i] << ",";
	}
	solvenew[n - 1] = kappa * solvenew[n - 2] + mu;
	OUTu << solvenew[n - 1] << "}}" << endl;
	current_t = current_t + dt;
	cout << current_t << endl;

	DeletingVector(F);
	DeletingVector(ldiag);
	DeletingVector(diag);
	DeletingVector(udiag);
	DeletingVector(b);
	DeletingVector(solveboof);
	DeletingVector(solve);
	solve = VectCopy(solvenew, n);
	DeletingVector(solvenew);
	DeletingVector(x);
	DeletingVector(t);
	DeletingVector(a);
	DeletingVector(A);
	DeletingVector(B);
	DeletingVector(C);
	//PrintingVector(solve, n);
	DeletingVector(solve);
}


void HeatTransfer_Example1(datatype sigma)
{
	char FileName[] = "x.txt";
	char FileName2[] = "t.txt";
	char FileName3[] = "1.txt";
	ofstream OUTx(FileName);
	ofstream OUTt(FileName2);
	ofstream OUTu(FileName3);

	datatype T = 0;
	datatype dt = 0;
	cout<<"Enter dt step:";
	cin>>dt;
	cout<<endl;
	cout<<"Enter time: ";
	cin>>T;
	cout<<endl;
	int m = T/dt;
	m++;
	funct f = Kx; //для вычисления коэффициентов a
	//funct P = P_Zero; // граничное условие потока
	//предполагаю решать задачу на временном промежутке от 0 до 0.7с  и с длиной стержня 1
	//datatype dt = 0.0001;
	datatype dh = 0.05;
	int n = 21; //всего узлов на отрезке от 0 до 1 с шагом h  0---20
	//int m = 7001; // временных слоёв на указанной мною сетке   0---7000

	datatype current_t = 0; //переменная текущего времени
	OUTu << "{" << endl;
	OUTu << "{";
	tVector solve = new datatype[n]; //массив точек температур в узлах сетки на текущем временном слое 0---20
	for (int i = 0; i < n - 1; i++)
	{
		solve[i] = 0.04 + i * dh * (1 - i*dh);  // начальные условия варианта 15
		OUTu << solve[i] << ",";
	}
	solve[n - 1] = 0.04;
	OUTu << 0.04 << "}," << endl; //запись первой строчки температур в файл

	OUTx << "{";
	tVector x = new datatype[n]; //массив пространственной сетки 0---20
	for (int i = 0; i < n - 1; i++)
	{
		x[i] = i * dh;
		OUTx << x[i] << ",";
	}
	x[n - 1] = (n - 1)*dh;
	OUTx << x[n - 1] << "}"; //запись простанственной сетки в файл
	OUTx.close();

	OUTt << "{";
	tVector t = new datatype[m];  // массив временной сетки 0---1000
	for (int i = 0; i < m - 1; i++)
	{
		t[i] = i * dt;
		OUTt << t[i] << ",";
	}
	t[m - 1] = (m - 1)*dt;
	OUTt << t[m - 1] << "}"; //запись временнОй сетки в файл
	OUTt.close();
	m--;

	tVector a = new datatype[n - 1]; //массив маленьких а для вычисления А, В, С (n-1 есть число промежутков интегрирования)1---20(0---19)
	for (int i = 0; i < n - 1; i++)
	{
		a[i] = dh / Integrate(f, x[i], x[i + 1]); //обрабатываются все отрезочки
	}
//	cout << "a:" << endl;
//	PrintingVector(a, n - 1);
//	cin.get();

	//в терминах методички N = 20

	tVector A = new datatype[n - 2];  //1---19(0---18)
	for (int i = 0; i < n - 2; i++)
	{
		A[i] = sigma * a[i] / dh;
	}
//	cout << "A:" << endl;
//	PrintingVector(A, n - 2);
//	cin.get();

	tVector B = new datatype[n - 2];  //1---19(0---18)
	for (int i = 0; i < n - 2; i++)
	{
		B[i] = sigma * a[i + 1] / dh;
	}
//	PrintingVector(B, n - 2);
//	cin.get();

	tVector C = new datatype[n - 2];  //1---19(0---18)
	for (int i = 0; i < n - 2; i++)
	{
		C[i] = A[i] + B[i] + dh / dt;
	}
//	PrintingVector(C, n - 2);
//	cin.get();

//	datatype kappa = (sigma * a[n - 2] / dh) / (dh / (2 * dt) + sigma * a[n - 2] / dh);
//	cout << kappa;
//	cin.get();

	for (int vr = 0; vr < m-1; vr++)
	{
		// задал все статичные параметры, осталось задать F и mu. Приступим:

		tVector F = new datatype[n - 2]; //1---19(0---18)
		for (int i = 0; i < n - 2; i++)  //старался сделать без ошибки, но мб она тут есть
		{
			F[i] = (dh / dt) * solve[i + 1] + (1 - sigma) * (a[i + 1] * (solve[i + 2] - solve[i + 1]) - a[i] * (solve[i + 1] - solve[i])) / dh;
		}

		//datatype mu = (solve[n - 1] * dh / (2 * dt) + sigma * P(current_t + dt) + (1 - sigma) * (P(current_t) - a[n - 2] * (solve[n - 1] - solve[n - 2]) / h))
		//	/ (dh / (2 * dt) + sigma * a[n - 2] / dh);


		//теперь коэффициенты трехдиагональной матрицы:

		tVector ldiag = new datatype[n - 3];  //0---17 (1--18)
		for (int i = 0; i < n - 3; i++)
		{
			ldiag[i] = A[i + 1];
		}

		tVector diag = new datatype[n - 2];  //0---18 (1--19)//записываю как последний ублюдок
		for (int i = 0; i < n - 3; i++)      //в терминах методички с минусами
		{
			diag[i] = -C[i];
		}
		//diag[n - 3] = -C[n - 3] + kappa * B[n - 3];
		diag[n - 3] = -C[n - 3];

		tVector udiag = new datatype[n - 3];//0---17 (1--18)
		for (int i = 0; i < n - 3; i++)
		{
			udiag[i] = B[i];
		}


		//Осталось записать только вектор правой части:

		tVector b = new datatype[n - 2]; // 0---18
		b[0] = -F[0] - 0.04 * A[0];
		for (int i = 1; i < n - 3; i++)
		{
			b[i] = -F[i];
		}
		b[n - 3] = -F[n - 3] - 0.04 * B[n - 3];

		//PrintingVector(b, n - 2);
		////////////////////////// Все готово для нахождения решения на новом слое.
		////////////////////////// Не уверен только в правильности индексации, но все должно быть ок.
		////////////////////////// Сопоставлял индексы в уме. Ты если будешь править ошибки, то лучше распиши на бумаге
		////////////////////////// В один момент щёлнет, как надо правильно индексировать.
		tVector solveboof = Shuttle(ldiag, diag, udiag, b, n - 2);
		tVector solvenew = new datatype[n]; // вектор решения на новом временном слое


		solvenew[0] = 0.04;
		OUTu << "{" << solvenew[0] << ",";
		for (int i = 1; i < n - 1; i++)
		{
			solvenew[i] = solveboof[i - 1];
			OUTu << solvenew[i] << ",";
		}
		solvenew[n - 1] = 0.04;
		OUTu << solvenew[n - 1] << "}," << endl;
		current_t = current_t + dt;
		cout << current_t << endl;

		DeletingVector(F);
		DeletingVector(ldiag);
		DeletingVector(diag);
		DeletingVector(udiag);
		DeletingVector(b);
		DeletingVector(solveboof);
		DeletingVector(solve);
		solve = VectCopy(solvenew, n);
		DeletingVector(solvenew);
		// добавь запись этого вектора в файл. Файл не такой, как pаньше, а другой.
		// Я рисовал на лабе, как должен выглядеть этот файл. Если что, напиши, я напомню.

	}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// задал все статичные параметры, осталось задать F и mu. Приступим:

			tVector F = new datatype[n - 2]; //1---19(0---18)
			for (int i = 0; i < n - 2; i++)  //старался сделать без ошибки, но мб она тут есть
			{
				F[i] = (dh / dt) * solve[i + 1] + (1 - sigma) * (a[i + 1] * (solve[i + 2] - solve[i + 1]) - a[i] * (solve[i + 1] - solve[i])) / dh;
			}

			//datatype mu = (solve[n - 1] * dh / (2 * dt) + sigma * P(current_t + dt) + (1 - sigma) * (P(current_t) - a[n - 2] * (solve[n - 1] - solve[n - 2]) / h))
			//	/ (dh / (2 * dt) + sigma * a[n - 2] / dh);


			//теперь коэффициенты трехдиагональной матрицы:

			tVector ldiag = new datatype[n - 3];  //0---17 (1--18)
			for (int i = 0; i < n - 3; i++)
			{
				ldiag[i] = A[i + 1];
			}

			tVector diag = new datatype[n - 2];  //0---18 (1--19)//записываю как последний ублюдок
			for (int i = 0; i < n - 3; i++)      //в терминах методички с минусами(похуй)
			{
				diag[i] = -C[i];
			}
			//diag[n - 3] = -C[n - 3] + kappa * B[n - 3];
			diag[n - 3] = -C[n - 3];

			tVector udiag = new datatype[n - 3];//0---17 (1--18)
			for (int i = 0; i < n - 3; i++)
			{
				udiag[i] = B[i];
			}


			//Осталось записать только вектор правой части:

			tVector b = new datatype[n - 2]; // 0---18
			b[0] = -F[0] - 0.04 * A[0];
			for (int i = 1; i < n - 3; i++)
			{
				b[i] = -F[i];
			}
			b[n - 3] = -F[n - 3] - 0.04 * B[n - 3];

			//PrintingVector(b, n - 2);
			////////////////////////// Все готово для нахождения решения на новом слое.
			////////////////////////// Не уверен только в правильности индексации, но все должно быть ок.
			////////////////////////// Сопоставлял индексы в уме. Ты если будешь править ошибки, то лучше распиши на бумаге
			////////////////////////// В один момент щёлнет, как надо правильно индексировать.
			////////////////////////// я сосал меня ебали
			tVector solveboof = Shuttle(ldiag, diag, udiag, b, n - 2);
			tVector solvenew = new datatype[n]; // вектор решения на новом временном слое


			solvenew[0] = 0.04;
			OUTu << "{" << solvenew[0] << ",";
			for (int i = 1; i < n - 1; i++)
			{
				solvenew[i] = solveboof[i - 1];
				OUTu << solvenew[i] << ",";
			}
			solvenew[n - 1] = 0.04;
			OUTu << solvenew[n - 1] << "}}" << endl;
			current_t = current_t + dt;
			cout << current_t << endl;

			DeletingVector(F);
			DeletingVector(ldiag);
			DeletingVector(diag);
			DeletingVector(udiag);
			DeletingVector(b);
			DeletingVector(solveboof);
			DeletingVector(solve);
			solve = VectCopy(solvenew, n);
			DeletingVector(solvenew);
			// добавь запись этого вектора в файл. Файл не такой, как pаньше, а другой.
			// Я рисовал на лабе, как должен выглядеть этот файл. Если что, напиши, я напомню.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//OUTu<<"}";
	DeletingVector(x);
	DeletingVector(t);
	DeletingVector(a);
	DeletingVector(A);
	DeletingVector(B);
	DeletingVector(C);
	//PrintingVector(solve, n);
	DeletingVector(solve);
}

void HeatTransfer_Example1_Anal(datatype sigma, funct2 &AnalSolve)
{
//	char FileName[] = "x.txt";
//	char FileName2[] = "t.txt";
//	char FileName3[] = "1.txt";
//	ofstream OUTx(FileName);
//	ofstream OUTt(FileName2);
//	ofstream OUTu(FileName3);
//	OUTu.setf(ios::fixed);
	datatype T = 0;
	datatype dt = 0;
	cout<<"Enter dt step:";
	cin>>dt;
	cout<<endl;
	cout<<"Enter time: ";
	cin>>T;
	cout<<endl;
	tVector Errors = new datatype[5];
	datatype difference;
	datatype maxError;

	//funct f = Kx; //для вычисления коэффициентов a
	funct f = K1;
	//funct P = P_Zero; // граничное условие потока
	//предполагаю решать задачу на временном промежутке от 0 до 0.7с  и с длиной стержня 1
	//datatype dt = 0.0001;
	datatype dh = 0.05;

for (int count = 0; count <5; count++)
{
	maxError = 0;
	difference = 0;
	int n = 1./dh+1; //всего узлов на отрезке от 0 до 1 с шагом h  0---20
	int m = T/dt +1; // временных слоёв на указанной мною сетке   0---7000

	datatype current_t = 0; //переменная текущего времени
	//OUTu << "{" << endl;
	//OUTu << "{";
	tVector solve = new datatype[n]; //массив точек температур в узлах сетки на текущем временном слое 0---20
	for (int i = 0; i < n - 1; i++)
	{
		solve[i] = sin(M_PI*i*dh);  // начальные условия варианта 15
		//OUTu << setprecision(17)<<solve[i] << ",";
	}
	solve[n - 1]=sin(M_PI*(n-1)*dh);
	//OUTu << setprecision(17)<<solve[n-1] << "}," << endl; //запись первой строчки температур в файл

	//OUTx << "{";
	tVector x = new datatype[n]; //массив пространственной сетки 0---20
	for (int i = 0; i < n - 1; i++)
	{
		x[i] = i * dh;
		//OUTx << x[i] << ",";
	}
	x[n - 1] = (n - 1)*dh;
	//OUTx << x[n - 1] << "}"; //запись простанственной сетки в файл
	//OUTx.close();

	//OUTt << "{";
	tVector t = new datatype[m];  // массив временной сетки 0---1000
	for (int i = 0; i < m - 1; i++)
	{
		t[i] = i * dt;
		//OUTt << t[i] << ",";
	}
	t[m - 1] = (m - 1)*dt;
	//OUTt << t[m - 1] << "}"; //запись временнОй сетки в файл
	//OUTt.close();
	m--;

	tVector a = new datatype[n - 1]; //массив маленьких а для вычисления А, В, С (n-1 есть число промежутков интегрирования)1---20(0---19)
	for (int i = 0; i < n - 1; i++)
	{
		//a[i] = dh / Integrate(f, x[i], x[i + 1]); //обрабатываются все отрезочки
		a[i]=0.5*(f(x[i])+f(x[i+1]));
	}
//	cout << "a:" << endl;
//	PrintingVector(a, n - 1);
//	cin.get();

	//в терминах методички N = 20

	tVector A = new datatype[n - 2];  //1---19(0---18)
	for (int i = 0; i < n - 2; i++)
	{
		A[i] = sigma * a[i] / dh;
	}
//	cout << "A:" << endl;
//	PrintingVector(A, n - 2);
//	cin.get();

	tVector B = new datatype[n - 2];  //1---19(0---18)
	for (int i = 0; i < n - 2; i++)
	{
		B[i] = sigma * a[i + 1] / dh;
	}
//	PrintingVector(B, n - 2);
//	cin.get();

	tVector C = new datatype[n - 2];  //1---19(0---18)
	for (int i = 0; i < n - 2; i++)
	{
		C[i] = A[i] + B[i] + dh / dt;
	}
//	PrintingVector(C, n - 2);
//	cin.get();

//	datatype kappa = (sigma * a[n - 2] / dh) / (dh / (2 * dt) + sigma * a[n - 2] / dh);
//	cout << kappa;
//	cin.get();

	for (int vr = 0; vr < m-1; vr++)
	{
		// задал все статичные параметры, осталось задать F и mu. Приступим:

		tVector F = new datatype[n - 2]; //1---19(0---18)
		for (int i = 0; i < n - 2; i++)  //старался сделать без ошибки, но мб она тут есть
		{
			F[i] = (dh / dt) * solve[i + 1] + (1 - sigma) * (a[i + 1] * (solve[i + 2] - solve[i + 1]) - a[i] * (solve[i + 1] - solve[i])) / dh;
		}

		//datatype mu = (solve[n - 1] * dh / (2 * dt) + sigma * P(current_t + dt) + (1 - sigma) * (P(current_t) - a[n - 2] * (solve[n - 1] - solve[n - 2]) / h))
		//	/ (dh / (2 * dt) + sigma * a[n - 2] / dh);


		//теперь коэффициенты трехдиагональной матрицы:

		tVector ldiag = new datatype[n - 3];  //0---17 (1--18)
		for (int i = 0; i < n - 3; i++)
		{
			ldiag[i] = A[i + 1];
		}

		tVector diag = new datatype[n - 2];  //0---18 (1--19)//записываю как последний ублюдок
		for (int i = 0; i < n - 3; i++)      //в терминах методички с минусами
		{
			diag[i] = -C[i];
		}
		//diag[n - 3] = -C[n - 3] + kappa * B[n - 3];
		diag[n - 3] = -C[n - 3];

		tVector udiag = new datatype[n - 3];//0---17 (1--18)
		for (int i = 0; i < n - 3; i++)
		{
			udiag[i] = B[i];
		}


		//Осталось записать только вектор правой части:

		tVector b = new datatype[n - 2]; // 0---18
		b[0] = -F[0] - solve[0] * A[0];
		for (int i = 1; i < n - 3; i++)
		{
			b[i] = -F[i];
		}
		b[n - 3] = -F[n - 3] - solve[n-1] * B[n - 3];

		//PrintingVector(b, n - 2);
		////////////////////////// Все готово для нахождения решения на новом слое.
		////////////////////////// Не уверен только в правильности индексации, но все должно быть ок.
		////////////////////////// Сопоставлял индексы в уме. Ты если будешь править ошибки, то лучше распиши на бумаге
		////////////////////////// В один момент щёлнет, как надо правильно индексировать.
		tVector solveboof = Shuttle(ldiag, diag, udiag, b, n - 2);
		tVector solvenew = new datatype[n]; // вектор решения на новом временном слое


		solvenew[0] = 0;
		//OUTu << "{" << solvenew[0] << ",";
		for (int i = 1; i < n - 1; i++)
		{
			solvenew[i] = solveboof[i - 1];
			//OUTu << setprecision(17) << solvenew[i] << ",";
		}
		solvenew[n - 1] = 0;
		//OUTu << setprecision(17)<< solvenew[n - 1] << "}," << endl;

		for(int z=0; z<n; z++)
		{
			difference = abs(solvenew[z]-AnalSolve(z*dh,current_t + dt));
			if (difference>maxError)
			{
				maxError = difference;
			}
		}

		current_t = current_t + dt;
		//cout << current_t << endl;

		DeletingVector(F);
		DeletingVector(ldiag);
		DeletingVector(diag);
		DeletingVector(udiag);
		DeletingVector(b);
		DeletingVector(solveboof);
		DeletingVector(solve);
		solve = VectCopy(solvenew, n);
		DeletingVector(solvenew);
		// добавь запись этого вектора в файл. Файл не такой, как pаньше, а другой.
		// Я рисовал на лабе, как должен выглядеть этот файл. Если что, напиши, я напомню.

	}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// задал все статичные параметры, осталось задать F и mu. Приступим:

			tVector F = new datatype[n - 2]; //1---19(0---18)
			for (int i = 0; i < n - 2; i++)  //старался сделать без ошибки, но мб она тут есть
			{
				F[i] = (dh / dt) * solve[i + 1] + (1 - sigma) * (a[i + 1] * (solve[i + 2] - solve[i + 1]) - a[i] * (solve[i + 1] - solve[i])) / dh;
			}

			//datatype mu = (solve[n - 1] * dh / (2 * dt) + sigma * P(current_t + dt) + (1 - sigma) * (P(current_t) - a[n - 2] * (solve[n - 1] - solve[n - 2]) / h))
			//	/ (dh / (2 * dt) + sigma * a[n - 2] / dh);


			//теперь коэффициенты трехдиагональной матрицы:

			tVector ldiag = new datatype[n - 3];  //0---17 (1--18)
			for (int i = 0; i < n - 3; i++)
			{
				ldiag[i] = A[i + 1];
			}

			tVector diag = new datatype[n - 2];  //0---18 (1--19)//записываю как последний ублюдок
			for (int i = 0; i < n - 3; i++)      //в терминах методички с минусами(похуй)
			{
				diag[i] = -C[i];
			}
			//diag[n - 3] = -C[n - 3] + kappa * B[n - 3];
			diag[n - 3] = -C[n - 3];

			tVector udiag = new datatype[n - 3];//0---17 (1--18)
			for (int i = 0; i < n - 3; i++)
			{
				udiag[i] = B[i];
			}


			//Осталось записать только вектор правой части:

			tVector b = new datatype[n - 2]; // 0---18
			b[0] = -F[0] - solve[0] * A[0];
			for (int i = 1; i < n - 3; i++)
			{
				b[i] = -F[i];
			}
			b[n - 3] = -F[n - 3] - solve[n-1] * B[n - 3];

			//PrintingVector(b, n - 2);
			////////////////////////// Все готово для нахождения решения на новом слое.
			////////////////////////// Не уверен только в правильности индексации, но все должно быть ок.
			////////////////////////// Сопоставлял индексы в уме. Ты если будешь править ошибки, то лучше распиши на бумаге
			////////////////////////// В один момент щёлнет, как надо правильно индексировать.
			////////////////////////// я сосал меня ебали
			tVector solveboof = Shuttle(ldiag, diag, udiag, b, n - 2);
			tVector solvenew = new datatype[n]; // вектор решения на новом временном слое


			solvenew[0] = 0;
			//OUTu << "{" << solvenew[0] << ",";
			for (int i = 1; i < n - 1; i++)
			{
				solvenew[i] = solveboof[i - 1];
				//OUTu << solvenew[i] << ",";
			}
			solvenew[n - 1] = 0;
			//OUTu << solvenew[n - 1] << "}}" << endl;
			//OUTu.close();
			for(int z=0; z<n; z++)
			{
				difference = abs(solvenew[z]-AnalSolve(z*dh,current_t + dt));
				if (difference>maxError)
				{
					maxError = difference;
				}
			}

			current_t = current_t + dt;
			//cout << current_t << endl;

			DeletingVector(F);
			DeletingVector(ldiag);
			DeletingVector(diag);
			DeletingVector(udiag);
			DeletingVector(b);
			DeletingVector(solveboof);
			DeletingVector(solve);
			solve = VectCopy(solvenew, n);
			DeletingVector(solvenew);
			// добавь запись этого вектора в файл. Файл не такой, как pаньше, а другой.
			// Я рисовал на лабе, как должен выглядеть этот файл. Если что, напиши, я напомню.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//OUTu<<"}";
	DeletingVector(x);
	DeletingVector(t);
	DeletingVector(a);
	DeletingVector(A);
	DeletingVector(B);
	DeletingVector(C);
	//PrintingVector(solve, n);
	DeletingVector(solve);
	Errors[count]=maxError;
	//cout<<setprecision(17)<<maxError<<endl;
	dt=dt/4;
	dh=dh/2;
}//end for определения максимальной ошибки
cout<<"       sigma=0, dt=dt/4, dh=dh/2"<<endl;
cout<<"Погрешность на шаге:         Отношение ошибок: "<<endl;
for (int i=0; i<4; i++)
{
	cout<<Errors[i]<<"                   "<<Errors[i]/Errors[i+1]<<endl;
}
cout<<Errors[4]<<"                     ---"<<endl;
DeletingVector(Errors);
}

void HeatTransfer_Example1_Quasilinear(datatype sigma)
{
	char FileName[] = "x.txt";
	char FileName2[] = "t.txt";
	char FileName3[] = "1q.txt";
	ofstream OUTx(FileName);
	ofstream OUTt(FileName2);
	ofstream OUTu(FileName3);

	funct f = P_Heat; //для вычисления коэффициентов a
	//funct P = P_Zero; // граничное условие потока
	//предполагаю решать задачу на временном промежутке от 0 до 0.7с  и с длиной стержня 1
	datatype dt = 0.0001;
	datatype dh = 0.05;
	int n = 21; //всего узлов на отрезке от 0 до 1 с шагом h  0---20
	int m = 1001; // временных слоёв на указанной мною сетке   0---7000

	datatype current_t = 0; //переменная текущего времени
	OUTu << "{" << endl;
	OUTu << "{";
	tVector solve = new datatype[n]; //массив точек температур в узлах сетки на текущем временном слое 0---20
	for (int i = 0; i < n - 1; i++)
	{
		solve[i] = 0.04 + i * dh * (1 - i*dh);  // начальные условия варианта 15
		OUTu << solve[i] << ",";
	}
	solve[n - 1] = 0.04;
	OUTu << 0.04 << "}," << endl; //запись первой строчки температур в файл

	OUTx << "{";
	tVector x = new datatype[n]; //массив пространственной сетки 0---20
	for (int i = 0; i < n - 1; i++)
	{
		x[i] = i * dh;
		OUTx << x[i] << ",";
	}
	x[n - 1] = (n - 1)*dh;
	OUTx << x[n - 1] << "}"; //запись простанственной сетки в файл
	OUTx.close();

	OUTt << "{";
	tVector t = new datatype[m];  // массив временной сетки 0---1000
	for (int i = 0; i < m - 1; i++)
	{
		t[i] = i * dt;
		OUTt << t[i] << ",";
	}
	t[m - 1] = (m - 1)*dt;
	OUTt << t[m - 1] << "}"; //запись временнОй сетки в файл
	OUTt.close();

	tVector a = new datatype[n - 1]; //массив маленьких а для вычисления А, В, С (n-1 есть число промежутков интегрирования)1---20(0---19)
	for (int i = 0; i < n - 1; i++)
	{
		//a[i] = dh / Integrate(f, x[i], x[i + 1]); //обрабатываются все отрезочки
		a[i] = 0.5*(f(solve[i+1]) + f(solve[i]));
	}
//	cout << "a:" << endl;
//	PrintingVector(a, n - 1);
//	cin.get();

	//в терминах методички N = 20

	tVector A = new datatype[n - 2];  //1---19(0---18)
	for (int i = 0; i < n - 2; i++)
	{
		A[i] = sigma * a[i] / dh;
	}
//	cout << "A:" << endl;
//	PrintingVector(A, n - 2);
//	cin.get();

	tVector B = new datatype[n - 2];  //1---19(0---18)
	for (int i = 0; i < n - 2; i++)
	{
		B[i] = sigma * a[i + 1] / dh;
	}
//	PrintingVector(B, n - 2);
//	cin.get();

	tVector C = new datatype[n - 2];  //1---19(0---18)
	for (int i = 0; i < n - 2; i++)
	{
		C[i] = A[i] + B[i] + dh / dt;
	}
//	PrintingVector(C, n - 2);
//	cin.get();

//	datatype kappa = (sigma * a[n - 2] / dh) / (dh / (2 * dt) + sigma * a[n - 2] / dh);
//	cout << kappa;
//	cin.get();

	while (current_t < 0.0999)
	{
		// задал все статичные параметры, осталось задать F и mu. Приступим:

		tVector F = new datatype[n - 2]; //1---19(0---18)
		for (int i = 0; i < n - 2; i++)  //старался сделать без ошибки, но мб она тут есть
		{
			F[i] = (dh / dt) * solve[i + 1] + (1 - sigma) * (a[i + 1] * (solve[i + 2] - solve[i + 1]) - a[i] * (solve[i + 1] - solve[i])) / dh;
		}

		//datatype mu = (solve[n - 1] * dh / (2 * dt) + sigma * P(current_t + dt) + (1 - sigma) * (P(current_t) - a[n - 2] * (solve[n - 1] - solve[n - 2]) / h))
		//	/ (dh / (2 * dt) + sigma * a[n - 2] / dh);


		//теперь коэффициенты трехдиагональной матрицы:

		tVector ldiag = new datatype[n - 3];  //0---17 (1--18)
		for (int i = 0; i < n - 3; i++)
		{
			ldiag[i] = A[i + 1];
		}

		tVector diag = new datatype[n - 2];  //0---18 (1--19)//записываю как последний ублюдок
		for (int i = 0; i < n - 3; i++)      //в терминах методички с минусами(похуй)
		{
			diag[i] = -C[i];
		}
		//diag[n - 3] = -C[n - 3] + kappa * B[n - 3];
		diag[n - 3] = -C[n - 3];

		tVector udiag = new datatype[n - 3];//0---17 (1--18)
		for (int i = 0; i < n - 3; i++)
		{
			udiag[i] = B[i];
		}


		//Осталось записать только вектор правой части:

		tVector b = new datatype[n - 2]; // 0---18
		b[0] = -F[0] - 0.04 * A[0];
		for (int i = 1; i < n - 3; i++)
		{
			b[i] = -F[i];
		}
		b[n - 3] = -F[n - 3] - 0.04 * B[n - 3];

		//PrintingVector(b, n - 2);
		////////////////////////// Все готово для нахождения решения на новом слое.
		////////////////////////// Не уверен только в правильности индексации, но все должно быть ок.
		////////////////////////// Сопоставлял индексы в уме. Ты если будешь править ошибки, то лучше распиши на бумаге
		////////////////////////// В один момент щёлнет, как надо правильно индексировать.
		////////////////////////// я сосал меня ебали
		tVector solveboof = Shuttle(ldiag, diag, udiag, b, n - 2);
		tVector solvenew = new datatype[n]; // вектор решения на новом временном слое


		solvenew[0] = 0.04;
		OUTu << "{" << solvenew[0] << ",";
		for (int i = 1; i < n - 1; i++)
		{
			solvenew[i] = solveboof[i - 1];
			OUTu << solvenew[i] << ",";
		}
		solvenew[n - 1] = 0.04;
		OUTu << solvenew[n - 1] << "}," << endl;
		current_t = current_t + dt;
		cout << current_t << endl;

		DeletingVector(F);
		DeletingVector(ldiag);
		DeletingVector(diag);
		DeletingVector(udiag);
		DeletingVector(b);
		DeletingVector(solveboof);
		DeletingVector(solve);
		solve = VectCopy(solvenew, n);
		for (int i = 0; i < n - 1; i++)
			{
				//a[i] = dh / Integrate(f, x[i], x[i + 1]); //обрабатываются все отрезочки
				a[i] = 0.5*(f(solve[i+1]) + f(solve[i]));  //записываются новые значения
			}
		DeletingVector(solvenew);
		// добавь запись этого вектора в файл. Файл не такой, как pаньше, а другой.
		// Я рисовал на лабе, как должен выглядеть этот файл. Если что, напиши, я напомню.

	}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// задал все статичные параметры, осталось задать F и mu. Приступим:

			tVector F = new datatype[n - 2]; //1---19(0---18)
			for (int i = 0; i < n - 2; i++)  //старался сделать без ошибки, но мб она тут есть
			{
				F[i] = (dh / dt) * solve[i + 1] + (1 - sigma) * (a[i + 1] * (solve[i + 2] - solve[i + 1]) - a[i] * (solve[i + 1] - solve[i])) / dh;
			}

			//datatype mu = (solve[n - 1] * dh / (2 * dt) + sigma * P(current_t + dt) + (1 - sigma) * (P(current_t) - a[n - 2] * (solve[n - 1] - solve[n - 2]) / h))
			//	/ (dh / (2 * dt) + sigma * a[n - 2] / dh);


			//теперь коэффициенты трехдиагональной матрицы:

			tVector ldiag = new datatype[n - 3];  //0---17 (1--18)
			for (int i = 0; i < n - 3; i++)
			{
				ldiag[i] = A[i + 1];
			}

			tVector diag = new datatype[n - 2];  //0---18 (1--19)//записываю как последний ублюдок
			for (int i = 0; i < n - 3; i++)      //в терминах методички с минусами(похуй)
			{
				diag[i] = -C[i];
			}
			//diag[n - 3] = -C[n - 3] + kappa * B[n - 3];
			diag[n - 3] = -C[n - 3];

			tVector udiag = new datatype[n - 3];//0---17 (1--18)
			for (int i = 0; i < n - 3; i++)
			{
				udiag[i] = B[i];
			}


			//Осталось записать только вектор правой части:

			tVector b = new datatype[n - 2]; // 0---18
			b[0] = -F[0] - 0.04 * A[0];
			for (int i = 1; i < n - 3; i++)
			{
				b[i] = -F[i];
			}
			b[n - 3] = -F[n - 3] - 0.04 * B[n - 3];

			//PrintingVector(b, n - 2);
			////////////////////////// Все готово для нахождения решения на новом слое.
			////////////////////////// Не уверен только в правильности индексации, но все должно быть ок.
			////////////////////////// Сопоставлял индексы в уме. Ты если будешь править ошибки, то лучше распиши на бумаге
			////////////////////////// В один момент щёлнет, как надо правильно индексировать.
			////////////////////////// я сосал меня ебали
			tVector solveboof = Shuttle(ldiag, diag, udiag, b, n - 2);
			tVector solvenew = new datatype[n]; // вектор решения на новом временном слое


			solvenew[0] = 0.04;
			OUTu << "{" << solvenew[0] << ",";
			for (int i = 1; i < n - 1; i++)
			{
				solvenew[i] = solveboof[i - 1];
				OUTu << solvenew[i] << ",";
			}
			solvenew[n - 1] = 0.04;
			OUTu << solvenew[n - 1] << "}}" << endl;
			current_t = current_t + dt;
			cout << current_t << endl;

			DeletingVector(F);
			DeletingVector(ldiag);
			DeletingVector(diag);
			DeletingVector(udiag);
			DeletingVector(b);
			DeletingVector(solveboof);
			DeletingVector(solve);
			solve = VectCopy(solvenew, n);
			DeletingVector(solvenew);
			// добавь запись этого вектора в файл. Файл не такой, как pаньше, а другой.
			// Я рисовал на лабе, как должен выглядеть этот файл. Если что, напиши, я напомню.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//OUTu<<"}";
	DeletingVector(x);
	DeletingVector(t);
	DeletingVector(a);
	DeletingVector(A);
	DeletingVector(B);
	DeletingVector(C);
	//PrintingVector(solve, n);
	DeletingVector(solve);
}


void HeatTransfer_Example1_Quasilinear_Iterations(datatype sigma)
{
	char FileName[] = "x.txt";
	char FileName2[] = "t.txt";
	char FileName3[] = "1q.txt";
	ofstream OUTx(FileName);
	ofstream OUTt(FileName2);
	ofstream OUTu(FileName3);
	size_t counter = 0;
	datatype norm = 100;
	funct f = P_Heat; //для вычисления коэффициентов a
	//funct P = P_Zero; // граничное условие потока
	//предполагаю решать задачу на временном промежутке от 0 до 0.7с  и с длиной стержня 1
	datatype dt = 0.0001;
	datatype dh = 0.05;
	int n = 21; //всего узлов на отрезке от 0 до 1 с шагом h  0---20
	int m = 1001; // временных слоёв на указанной мною сетке   0---7000

	datatype current_t = 0; //переменная текущего времени
	OUTu << "{" << endl;
	OUTu << "{";
	tVector solve = new datatype[n]; //массив точек температур в узлах сетки на текущем временном слое 0---20
	for (int i = 0; i < n - 1; i++)
	{
		solve[i] = 0.04 + i * dh * (1 - i*dh);  // начальные условия варианта 15
		OUTu << solve[i] << ",";
	}
	solve[n - 1] = 0.04;
	OUTu << 0.04 << "}," << endl; //запись первой строчки температур в файл

	OUTx << "{";
	tVector x = new datatype[n]; //массив пространственной сетки 0---20
	for (int i = 0; i < n - 1; i++)
	{
		x[i] = i * dh;
		OUTx << x[i] << ",";
	}
	x[n - 1] = (n - 1)*dh;
	OUTx << x[n - 1] << "}"; //запись простанственной сетки в файл
	OUTx.close();

	OUTt << "{";
	tVector t = new datatype[m];  // массив временной сетки 0---7000
	for (int i = 0; i < m - 1; i++)
	{
		t[i] = i * dt;
		OUTt << t[i] << ",";
	}
	t[m - 1] = (m - 1)*dt;
	OUTt << t[m - 1] << "}"; //запись временнОй сетки в файл
	OUTt.close();

	tVector a = new datatype[n - 1]; //массив маленьких а для вычисления А, В, С (n-1 есть число промежутков интегрирования)1---20(0---19)


	//в терминах методички N = 20

	tVector A = new datatype[n - 2];  //1---19(0---18)


	tVector B = new datatype[n - 2];  //1---19(0---18)


	tVector C = new datatype[n - 2];  //1---19(0---18)


	while (current_t < 0.0999)
	{
		counter = 0;
		norm=100;
		tVector constantU = VectCopy(solve, n);
		while (norm>eps5)
		{
			for (int i = 0; i < n - 1; i++)
				{
					a[i] = 0.5*(f(solve[i+1]) + f(solve[i]));
				}
			for (int i = 0; i < n - 2; i++)
				{
					A[i] = sigma * a[i] / dh;
				}
			for (int i = 0; i < n - 2; i++)
				{
					B[i] = sigma * a[i + 1] / dh;
				}
			for (int i = 0; i < n - 2; i++)
				{
					C[i] = A[i] + B[i] + dh / dt;
				}
			tVector F = new datatype[n - 2]; //1---19(0---18)
			for (int i = 0; i < n - 2; i++)  //старался сделать без ошибки, но мб она тут есть
			{
				F[i] = (dh / dt) * constantU[i + 1] + (1 - sigma) * (a[i + 1] * (solve[i + 2] - solve[i + 1]) - a[i] * (solve[i + 1] - solve[i])) / dh;
			}

			//теперь коэффициенты трехдиагональной матрицы:

			tVector ldiag = new datatype[n - 3];  //0---17 (1--18)
			for (int i = 0; i < n - 3; i++)
			{
				ldiag[i] = A[i + 1];
			}

			tVector diag = new datatype[n - 2];  //0---18 (1--19)//записываю как последний ублюдок
			for (int i = 0; i < n - 3; i++)      //в терминах методички с минусами(похуй)
			{
				diag[i] = -C[i];
			}
			//diag[n - 3] = -C[n - 3] + kappa * B[n - 3];
			diag[n - 3] = -C[n - 3];

			tVector udiag = new datatype[n - 3];//0---17 (1--18)
			for (int i = 0; i < n - 3; i++)
			{
				udiag[i] = B[i];
			}


			//Осталось записать только вектор правой части:

			tVector b = new datatype[n - 2]; // 0---18
			b[0] = -F[0] - 0.04 * A[0];
			for (int i = 1; i < n - 3; i++)
			{
				b[i] = -F[i];
			}
			b[n - 3] = -F[n - 3] - 0.04 * B[n - 3];

			tVector solveboof = Shuttle(ldiag, diag, udiag, b, n - 2);
			tVector solvenew = new datatype[n]; // вектор решения на новом временном слое


			solvenew[0] = 0.04;
			for (int i = 1; i < n - 1; i++)
			{
				solvenew[i] = solveboof[i - 1];
			}
			solvenew[n - 1] = 0.04;

			tVector DX = Difference(solvenew, solve, n);
			norm = NormOfVectorC(DX, n);
			DeletingVector(DX);
			counter++;

			DeletingVector(F);
			DeletingVector(ldiag);
			DeletingVector(diag);
			DeletingVector(udiag);
			DeletingVector(b);
			DeletingVector(solveboof);
			DeletingVector(solve);
			solve = VectCopy(solvenew, n); //запись s-ой итерации решения
			DeletingVector(solvenew);
			cout<<norm<<endl;
		}
		cout<<"Iterations before convergence: "<<counter<<endl;

		OUTu << "{" << solve[0] << ",";
		for (int i = 1; i < n - 1; i++)
			{
				OUTu << solve[i] << ",";
			}
		OUTu << solve[n - 1] << "}," << endl;
		current_t = current_t + dt;
		//cout << current_t <<endl;
		DeletingVector(constantU);
	}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// задал все статичные параметры, осталось задать F и mu. Приступим:
	counter=0;
	norm=100;
	tVector constantU = VectCopy(solve, n);
	while (norm>eps5)
	{
		for (int i = 0; i < n - 1; i++)
			{
				a[i] = 0.5*(f(solve[i+1]) + f(solve[i]));
			}
		for (int i = 0; i < n - 2; i++)
			{
				A[i] = sigma * a[i] / dh;
			}
		for (int i = 0; i < n - 2; i++)
			{
				B[i] = sigma * a[i + 1] / dh;
			}
		for (int i = 0; i < n - 2; i++)
			{
				C[i] = A[i] + B[i] + dh / dt;
			}
		tVector F = new datatype[n - 2]; //1---19(0---18)
		for (int i = 0; i < n - 2; i++)  //старался сделать без ошибки, но мб она тут есть
		{
			F[i] = (dh / dt) * solve[i + 1] + (1 - sigma) * (a[i + 1] * (solve[i + 2] - solve[i + 1]) - a[i] * (solve[i + 1] - solve[i])) / dh;
		}

		//datatype mu = (solve[n - 1] * dh / (2 * dt) + sigma * P(current_t + dt) + (1 - sigma) * (P(current_t) - a[n - 2] * (solve[n - 1] - solve[n - 2]) / h))
		//	/ (dh / (2 * dt) + sigma * a[n - 2] / dh);


		//теперь коэффициенты трехдиагональной матрицы:

		tVector ldiag = new datatype[n - 3];  //0---17 (1--18)
		for (int i = 0; i < n - 3; i++)
		{
			ldiag[i] = A[i + 1];
		}

		tVector diag = new datatype[n - 2];  //0---18 (1--19)//записываю как последний ублюдок
		for (int i = 0; i < n - 3; i++)      //в терминах методички с минусами(похуй)
		{
			diag[i] = -C[i];
		}
		//diag[n - 3] = -C[n - 3] + kappa * B[n - 3];
		diag[n - 3] = -C[n - 3];

		tVector udiag = new datatype[n - 3];//0---17 (1--18)
		for (int i = 0; i < n - 3; i++)
		{
			udiag[i] = B[i];
		}


		//Осталось записать только вектор правой части:

		tVector b = new datatype[n - 2]; // 0---18
		b[0] = -F[0] - 0.04 * A[0];
		for (int i = 1; i < n - 3; i++)
		{
			b[i] = -F[i];
		}
		b[n - 3] = -F[n - 3] - 0.04 * B[n - 3];


		tVector solveboof = Shuttle(ldiag, diag, udiag, b, n - 2);
		tVector solvenew = new datatype[n]; // вектор решения на новом временном слое


		solvenew[0] = 0.04;
		for (int i = 1; i < n - 1; i++)
		{
			solvenew[i] = solveboof[i - 1];
		}
		solvenew[n - 1] = 0.04;

		tVector DX=Difference(solve, solvenew, n);
		norm = NormOfVectorC(DX, n);
		DeletingVector(DX);
		counter++;

		DeletingVector(F);
		DeletingVector(ldiag);
		DeletingVector(diag);
		DeletingVector(udiag);
		DeletingVector(b);
		DeletingVector(solveboof);
		DeletingVector(solve);
		solve = VectCopy(solvenew, n);
		DeletingVector(solvenew);
	}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	cout<<"Iterations before convergence: "<<counter<<endl;
	OUTu<<"{"<<solve[0] << ",";
	for (int i = 1; i < n - 1; i++)
		{
			OUTu << solve[i] << ",";
		}
	OUTu << solve[n - 1] << "}}" << endl;
	//OUTu<<"}";
	DeletingVector(x);
	DeletingVector(t);
	DeletingVector(a);
	DeletingVector(A);
	DeletingVector(B);
	DeletingVector(C);
	//PrintingVector(solve, n);
	DeletingVector(solve);
	DeletingVector(constantU);
}

datatype P100(datatype x)
{
	return 100.;
}

void HeatTransfer_Example2(datatype sigma)
{
	char FileName[] = "x.txt";
	char FileName2[] = "t.txt";
	char FileName3[] = "2.txt";
	ofstream OUTx(FileName);
	ofstream OUTt(FileName2);
	ofstream OUTu(FileName3);

	datatype T = 0;
	datatype dt = 0;
	cout<<"Enter dt step:";
	cin>>dt;
	cout<<endl;
	cout<<"Enter time: ";
	cin>>T;
	cout<<endl;
	int m = T/dt;
	m++;

	//funct f = Kxpr; //для вычисления коэффициентов a
	//funct K=Kx;
	funct K=K1;
	funct Pright = P_Zero; // граничное условие потока
	funct Pleft = P_Zero;
	//предполагаю решать задачу на временном промежутке от 0 до 0.7с  и с длиной стержня 1
	//datatype dt = 0.0001;
	datatype dh = 0.05;
	int n = 21; //всего узлов на отрезке от 0 до 1 с шагом h  0---20
	//int m = 7001; // временных слоёв на указанной мною сетке   0---1000

	datatype current_t = 0; //переменная текущего времени
	OUTu << "{" << endl;
	OUTu << "{";
	tVector solve = new datatype[n]; //массив точек температур в узлах сетки на текущем временном слое 0---20
	//datatype Energy_i = 0;
	for (int i = 0; i < n - 1; i++)
	{
		solve[i] = 0.04 + i * dh * (1 - i*dh);  // начальные условия варианта 15
		OUTu << solve[i] << ",";
	}
	solve[n - 1] = 0.04;
	OUTu << 0.04 << "}," << endl; //запись первой строчки температур в файл

	OUTx << "{";
	tVector x = new datatype[n]; //массив пространственной сетки 0---20
	for (int i = 0; i < n - 1; i++)
	{
		x[i] = i * dh;
		OUTx << x[i] << ",";
	}
	x[n - 1] = (n - 1)*dh;
	OUTx << x[n - 1] << "}"; //запись простанственной сетки в файл
	OUTx.close();

	OUTt << "{";
	tVector t = new datatype[m];  // массив временной сетки 0---1000
	for (int i = 0; i < m - 1; i++)
	{
		t[i] = i * dt;
		OUTt << t[i] << ",";
	}
	t[m - 1] = (m - 1)*dt;
	OUTt << t[m - 1] << "}"; //запись временнОй сетки в файл
	OUTt.close();
	m--;

	tVector a = new datatype[n - 1]; //массив маленьких а для вычисления А, В, С (n-1 есть число промежутков интегрирования)1---20(0---19)
	for (int i = 0; i < n - 1; i++)
	{
		//a[i] = dh / Integrate(f, x[i], x[i + 1]); //обрабатываются все отрезочки
		a[i]=0.5*(K(x[i])+K(x[i+1]));
	}
//	cout << "a:" << endl;
//	PrintingVector(a, n - 1);
//	cin.get();

	//в терминах методички N = 20

	tVector A = new datatype[n - 2];  //1---19(0---18)
	for (int i = 0; i < n - 2; i++)
	{
		A[i] = sigma * a[i] / dh;
	}
//	cout << "A:" << endl;
//	PrintingVector(A, n - 2);
//	cin.get();

	tVector B = new datatype[n - 2];  //1---19(0---18)
	for (int i = 0; i < n - 2; i++)
	{
		B[i] = sigma * a[i + 1] / dh;
	}
//	PrintingVector(B, n - 2);
//	cin.get();

	tVector C = new datatype[n - 2];  //1---19(0---18)
	for (int i = 0; i < n - 2; i++)
	{
		C[i] = A[i] + B[i] + dh / dt;
	}
//	PrintingVector(C, n - 2);
//	cin.get();

	datatype kappa_r = (sigma * a[n - 2] / dh) / (dh / (2 * dt) + sigma * a[n - 2] / dh);
	//datatype kappa_l = (sigma*dt*a[0])/(dt*sigma*a[0] + 0.5*dh*dh);
	datatype kappa_l = (sigma*dt*a[0])/(dh*(dh/2+dt*sigma*a[0]/dh));

//	cout << kappa;
//	cin.get();

	for (int vr=0; vr<m-1; vr++)
	{
		// задал все статичные параметры, осталось задать F и mu. Приступим:

		tVector F = new datatype[n - 2]; //1---19(0---18)
		for (int i = 0; i < n - 2; i++)  //старался сделать без ошибки, но мб она тут есть
		{
			F[i] = (dh / dt) * solve[i + 1] + (1 - sigma) * (a[i + 1] * (solve[i + 2] - solve[i + 1]) - a[i] * (solve[i + 1] - solve[i])) / dh;
		}

		datatype mu_r = (solve[n - 1] * dh / (2 * dt) + sigma * Pright(current_t + dt) + (1 - sigma) * (Pright(current_t) - a[n - 2] * (solve[n - 1] - solve[n - 2]) / dh))
			/ (dh / (2 * dt) + sigma * a[n - 2] / dh);
//		datatype mu_l = (dt*(1-sigma)*a[0]*(solve[1]-solve[0]) + 0.5*dh*dh*solve[0])
//				/(0.5*dh*dh + dt*sigma*a[0]);
		datatype mu_l = (dh*solve[0]/2 + dt*sigma*Pleft(current_t + dt) + a[0]*dt*(1-sigma)/dh*(solve[1]-solve[0]) +dt*(1-sigma)*Pleft(current_t))/(dh/2+dt*sigma*a[0]/dh);


		//теперь коэффициенты трехдиагональной матрицы:

		tVector ldiag = new datatype[n - 3];  //0---17 (1--18)
		for (int i = 0; i < n - 3; i++)
		{
			ldiag[i] = A[i + 1];
		}

		tVector diag = new datatype[n - 2];  //0---18 (1--19)//записываю как последний ублюдок
		diag[0] = -C[0] + A[0]*kappa_l;
		for (int i = 1; i < n - 3; i++)      //в терминах методички с минусами(похуй)
		{
			diag[i] = -C[i];
		}
		diag[n - 3] = -C[n - 3] + kappa_r * B[n - 3];

		tVector udiag = new datatype[n - 3];//0---17 (1--18)
		for (int i = 0; i < n - 3; i++)
		{
			udiag[i] = B[i];
		}


		//Осталось записать только вектор правой части:

		tVector b = new datatype[n - 2]; // 0---18
		b[0] = -F[0] - mu_l * A[0];
		for (int i = 1; i < n - 3; i++)
		{
			b[i] = -F[i];
		}
		b[n - 3] = -F[n - 3] - mu_r * B[n - 3];

		//PrintingVector(b, n - 2);
		////////////////////////// Все готово для нахождения решения на новом слое.
		////////////////////////// Не уверен только в правильности индексации, но все должно быть ок.
		////////////////////////// Сопоставлял индексы в уме. Ты если будешь править ошибки, то лучше распиши на бумаге
		////////////////////////// В один момент щёлнет, как надо правильно индексировать.
		////////////////////////// я сосал меня ебали
		tVector solveboof = Shuttle(ldiag, diag, udiag, b, n - 2);
		tVector solvenew = new datatype[n]; // вектор решения на новом временном слое


		solvenew[0] = solveboof[0]*kappa_l + mu_l;
		OUTu << "{" << solvenew[0] << ",";
		for (int i = 1; i < n - 1; i++)
		{
			solvenew[i] = solveboof[i - 1];
			OUTu << solvenew[i] << ",";
		}
		solvenew[n - 1] = solveboof[n - 3]*kappa_r + mu_r;
		OUTu << solvenew[n - 1] << "}," << endl;
		current_t = current_t + dt;
		cout << current_t << endl;

		DeletingVector(F);
		DeletingVector(ldiag);
		DeletingVector(diag);
		DeletingVector(udiag);
		DeletingVector(b);
		DeletingVector(solveboof);
		DeletingVector(solve);
		solve = VectCopy(solvenew, n);
		DeletingVector(solvenew);
		// добавь запись этого вектора в файл. Файл не такой, как pаньше, а другой.
		// Я рисовал на лабе, как должен выглядеть этот файл. Если что, напиши, я напомню.

	}
	tVector F = new datatype[n - 2]; //1---19(0---18)
	for (int i = 0; i < n - 2; i++)  //старался сделать без ошибки, но мб она тут есть
	{
		F[i] = (dh / dt) * solve[i + 1] + (1 - sigma) * (a[i + 1] * (solve[i + 2] - solve[i + 1]) - a[i] * (solve[i + 1] - solve[i])) / dh;
	}

	datatype mu_r = (solve[n - 1] * dh / (2 * dt) + sigma * Pright(current_t + dt) + (1 - sigma) * (Pright(current_t) - a[n - 2] * (solve[n - 1] - solve[n - 2]) / dh))
		/ (dh / (2 * dt) + sigma * a[n - 2] / dh);
	datatype mu_l = (dt*(1-sigma)*a[0]*(solve[1]-solve[0]) + 0.5*dh*dh*solve[0])
			/(0.5*dh*dh + dt*sigma*a[0]);


	//теперь коэффициенты трехдиагональной матрицы:

	tVector ldiag = new datatype[n - 3];  //0---17 (1--18)
	for (int i = 0; i < n - 3; i++)
	{
		ldiag[i] = A[i + 1];
	}

	tVector diag = new datatype[n - 2];  //0---18 (1--19)//записываю как последний ублюдок
	diag[0] = -C[0] + A[0]*kappa_l;
	for (int i = 1; i < n - 3; i++)      //в терминах методички с минусами(похуй)
	{
		diag[i] = -C[i];
	}
	diag[n - 3] = -C[n - 3] + kappa_r * B[n - 3];

	tVector udiag = new datatype[n - 3];//0---17 (1--18)
	for (int i = 0; i < n - 3; i++)
	{
		udiag[i] = B[i];
	}


	//Осталось записать только вектор правой части:

	tVector b = new datatype[n - 2]; // 0---18
	b[0] = -F[0] - mu_l * A[0];
	for (int i = 1; i < n - 3; i++)
	{
		b[i] = -F[i];
	}
	b[n - 3] = -F[n - 3] - mu_r * B[n - 3];

	//PrintingVector(b, n - 2);
	////////////////////////// Все готово для нахождения решения на новом слое.
	////////////////////////// Не уверен только в правильности индексации, но все должно быть ок.
	////////////////////////// Сопоставлял индексы в уме. Ты если будешь править ошибки, то лучше распиши на бумаге
	////////////////////////// В один момент щёлнет, как надо правильно индексировать.
	////////////////////////// я сосал меня ебали
	tVector solveboof = Shuttle(ldiag, diag, udiag, b, n - 2);
	tVector solvenew = new datatype[n]; // вектор решения на новом временном слое


	solvenew[0] = solveboof[0]*kappa_l + mu_l;
	OUTu << "{" << solvenew[0] << ",";
	for (int i = 1; i < n - 1; i++)
	{
		solvenew[i] = solveboof[i - 1];
		OUTu << solvenew[i] << ",";
	}
	solvenew[n - 1] = solveboof[n - 3]*kappa_r + mu_r;
	OUTu << solvenew[n - 1] << "}}" << endl;
	current_t = current_t + dt;
	cout << current_t << endl;

	DeletingVector(F);
	DeletingVector(ldiag);
	DeletingVector(diag);
	DeletingVector(udiag);
	DeletingVector(b);
	DeletingVector(solveboof);
	DeletingVector(solve);
	solve = VectCopy(solvenew, n);
	DeletingVector(solvenew);
	//OUTu<<"}";
	DeletingVector(x);
	DeletingVector(t);
	DeletingVector(a);
	DeletingVector(A);
	DeletingVector(B);
	DeletingVector(C);
	//PrintingVector(solve, n);
	DeletingVector(solve);
}

datatype HeatAnal(datatype x, datatype t)
{
	return sin(M_PI*x)*exp(-t*M_PI*M_PI);
}

void HeatTransfer_Example2_Anal(datatype sigma, funct2 &AnalSolve)
{

	char FileName[] = "x.txt";
	char FileName2[] = "t.txt";
	char FileName3[] = "2Anal.txt";
	ofstream OUTx(FileName);
	ofstream OUTt(FileName2);
	ofstream OUTu(FileName3);

	datatype T = 0;
	datatype dt = 0;
	cout<<"Enter dt step:";
	cin>>dt;
	cout<<endl;
	cout<<"Enter time: ";
	cin>>T;
	cout<<endl;
	tVector Errors = new datatype[6];
	datatype difference=0;

	//funct f = Kxpr; //для вычисления коэффициентов a
	//funct K=Kx;
	funct K=K1;
	funct P = P_Zero; // граничное условие потока
	//предполагаю решать задачу на временном промежутке от 0 до 0.7с  и с длиной стержня 1
	//datatype dt = 0.0001;
	datatype dh = 0.1;
for (int count = 0; count <6; count++)
{
	datatype maxError = 0;
	int n = 1./dh + 1; //всего узлов на отрезке от 0 до 1 с шагом h  0---20
	int m = T/dt +1;
	//int m = 7001; // временных слоёв на указанной мною сетке   0---1000

	datatype current_t = 0; //переменная текущего времени
	OUTu << "{" << endl;
	OUTu << "{";
	tVector solve = new datatype[n]; //массив точек температур в узлах сетки на текущем временном слое 0---20
	//datatype Energy_i = 0;
	for (int i = 0; i < n; i++)
	{
		solve[i] = sin(i*dh);  // начальные условия варианта 15
		OUTu << solve[i] << ",";
	}
	OUTu << solve[n - 1] << "}," << endl; //запись первой строчки температур в файл

	OUTx << "{";
	tVector x = new datatype[n]; //массив пространственной сетки 0---20
	for (int i = 0; i < n - 1; i++)
	{
		x[i] = i * dh;
		OUTx << x[i] << ",";
	}
	x[n - 1] = (n - 1)*dh;
	OUTx << x[n - 1] << "}"; //запись простанственной сетки в файл
	OUTx.close();

	OUTt << "{";
	tVector t = new datatype[m];  // массив временной сетки 0---1000
	for (int i = 0; i < m - 1; i++)
	{
		t[i] = i * dt;
		OUTt << t[i] << ",";
	}
	t[m - 1] = (m - 1)*dt;
	OUTt << t[m - 1] << "}"; //запись временнОй сетки в файл
	OUTt.close();
	m--;

	tVector a = new datatype[n - 1]; //массив маленьких а для вычисления А, В, С (n-1 есть число промежутков интегрирования)1---20(0---19)
	for (int i = 0; i < n - 1; i++)
	{
		//a[i] = dh / Integrate(f, x[i], x[i + 1]); //обрабатываются все отрезочки
		a[i]=0.5*(K(x[i])+K(x[i+1]));
	}
//	cout << "a:" << endl;
//	PrintingVector(a, n - 1);
//	cin.get();

	//в терминах методички N = 20

	tVector A = new datatype[n - 2];  //1---19(0---18)
	for (int i = 0; i < n - 2; i++)
	{
		A[i] = sigma * a[i] / dh;
	}
//	cout << "A:" << endl;
//	PrintingVector(A, n - 2);
//	cin.get();

	tVector B = new datatype[n - 2];  //1---19(0---18)
	for (int i = 0; i < n - 2; i++)
	{
		B[i] = sigma * a[i + 1] / dh;
	}
//	PrintingVector(B, n - 2);
//	cin.get();

	tVector C = new datatype[n - 2];  //1---19(0---18)
	for (int i = 0; i < n - 2; i++)
	{
		C[i] = A[i] + B[i] + dh / dt;
	}
//	PrintingVector(C, n - 2);
//	cin.get();

	datatype kappa_r = (sigma * a[n - 2] / dh) / (dh / (2 * dt) + sigma * a[n - 2] / dh);
	datatype kappa_l = (sigma*dt*a[0])/(dt*sigma*a[0] + 0.5*dh*dh);
//	cout << kappa;
//	cin.get();

	for (int vr=0; vr<m-1; vr++)
	{
		// задал все статичные параметры, осталось задать F и mu. Приступим:

		tVector F = new datatype[n - 2]; //1---19(0---18)
		for (int i = 0; i < n - 2; i++)  //старался сделать без ошибки, но мб она тут есть
		{
			F[i] = (dh / dt) * solve[i + 1] + (1 - sigma) * (a[i + 1] * (solve[i + 2] - solve[i + 1]) - a[i] * (solve[i + 1] - solve[i])) / dh;
		}

		datatype mu_r = (solve[n - 1] * dh / (2 * dt) + sigma * P(current_t + dt) + (1 - sigma) * (P(current_t) - a[n - 2] * (solve[n - 1] - solve[n - 2]) / dh))
			/ (dh / (2 * dt) + sigma * a[n - 2] / dh);
		datatype mu_l = (dt*(1-sigma)*a[0]*(solve[1]-solve[0]) + 0.5*dh*dh*solve[0])
				/(0.5*dh*dh + dt*sigma*a[0]);


		//теперь коэффициенты трехдиагональной матрицы:

		tVector ldiag = new datatype[n - 3];  //0---17 (1--18)
		for (int i = 0; i < n - 3; i++)
		{
			ldiag[i] = A[i + 1];
		}

		tVector diag = new datatype[n - 2];  //0---18 (1--19)//записываю как последний ублюдок
		diag[0] = -C[0] + A[0]*kappa_l;
		for (int i = 1; i < n - 3; i++)      //в терминах методички с минусами(похуй)
		{
			diag[i] = -C[i];
		}
		diag[n - 3] = -C[n - 3] + kappa_r * B[n - 3];

		tVector udiag = new datatype[n - 3];//0---17 (1--18)
		for (int i = 0; i < n - 3; i++)
		{
			udiag[i] = B[i];
		}


		//Осталось записать только вектор правой части:

		tVector b = new datatype[n - 2]; // 0---18
		b[0] = -F[0] - mu_l * A[0];
		for (int i = 1; i < n - 3; i++)
		{
			b[i] = -F[i];
		}
		b[n - 3] = -F[n - 3] - mu_r * B[n - 3];

		//PrintingVector(b, n - 2);
		////////////////////////// Все готово для нахождения решения на новом слое.
		////////////////////////// Не уверен только в правильности индексации, но все должно быть ок.
		////////////////////////// Сопоставлял индексы в уме. Ты если будешь править ошибки, то лучше распиши на бумаге
		////////////////////////// В один момент щёлнет, как надо правильно индексировать.
		////////////////////////// я сосал меня ебали
		tVector solveboof = Shuttle(ldiag, diag, udiag, b, n - 2);
		tVector solvenew = new datatype[n]; // вектор решения на новом временном слое


		solvenew[0] = solveboof[0]*kappa_l + mu_l;
		OUTu << "{" << solvenew[0] << ",";
		for (int i = 1; i < n - 1; i++)
		{
			solvenew[i] = solveboof[i - 1];
			OUTu << solvenew[i] << ",";
		}
		solvenew[n - 1] = solveboof[n - 3]*kappa_r + mu_r;
		OUTu << solvenew[n - 1] << "}," << endl;

		for(int z=0; z<n; z++)
		{
			difference = abs(solvenew[z]-AnalSolve(z*dh,current_t));
			if (difference>maxError)
			{
				maxError = difference;
			}
		}

		current_t = current_t + dt;
		cout << current_t << endl;

		DeletingVector(F);
		DeletingVector(ldiag);
		DeletingVector(diag);
		DeletingVector(udiag);
		DeletingVector(b);
		DeletingVector(solveboof);
		DeletingVector(solve);
		solve = VectCopy(solvenew, n);
		DeletingVector(solvenew);
		// добавь запись этого вектора в файл. Файл не такой, как pаньше, а другой.
		// Я рисовал на лабе, как должен выглядеть этот файл. Если что, напиши, я напомню.

	}
	tVector F = new datatype[n - 2]; //1---19(0---18)
	for (int i = 0; i < n - 2; i++)  //старался сделать без ошибки, но мб она тут есть
	{
		F[i] = (dh / dt) * solve[i + 1] + (1 - sigma) * (a[i + 1] * (solve[i + 2] - solve[i + 1]) - a[i] * (solve[i + 1] - solve[i])) / dh;
	}

	datatype mu_r = (solve[n - 1] * dh / (2 * dt) + sigma * P(current_t + dt) + (1 - sigma) * (P(current_t) - a[n - 2] * (solve[n - 1] - solve[n - 2]) / dh))
		/ (dh / (2 * dt) + sigma * a[n - 2] / dh);
	datatype mu_l = (dt*(1-sigma)*a[0]*(solve[1]-solve[0]) + 0.5*dh*dh*solve[0])
			/(0.5*dh*dh + dt*sigma*a[0]);


	//теперь коэффициенты трехдиагональной матрицы:

	tVector ldiag = new datatype[n - 3];  //0---17 (1--18)
	for (int i = 0; i < n - 3; i++)
	{
		ldiag[i] = A[i + 1];
	}

	tVector diag = new datatype[n - 2];  //0---18 (1--19)//записываю как последний ублюдок
	diag[0] = -C[0] + A[0]*kappa_l;
	for (int i = 1; i < n - 3; i++)      //в терминах методички с минусами(похуй)
	{
		diag[i] = -C[i];
	}
	diag[n - 3] = -C[n - 3] + kappa_r * B[n - 3];

	tVector udiag = new datatype[n - 3];//0---17 (1--18)
	for (int i = 0; i < n - 3; i++)
	{
		udiag[i] = B[i];
	}


	//Осталось записать только вектор правой части:

	tVector b = new datatype[n - 2]; // 0---18
	b[0] = -F[0] - mu_l * A[0];
	for (int i = 1; i < n - 3; i++)
	{
		b[i] = -F[i];
	}
	b[n - 3] = -F[n - 3] - mu_r * B[n - 3];

	//PrintingVector(b, n - 2);
	////////////////////////// Все готово для нахождения решения на новом слое.
	////////////////////////// Не уверен только в правильности индексации, но все должно быть ок.
	////////////////////////// Сопоставлял индексы в уме. Ты если будешь править ошибки, то лучше распиши на бумаге
	////////////////////////// В один момент щёлнет, как надо правильно индексировать.
	////////////////////////// я сосал меня ебали
	tVector solveboof = Shuttle(ldiag, diag, udiag, b, n - 2);
	tVector solvenew = new datatype[n]; // вектор решения на новом временном слое


	solvenew[0] = solveboof[0]*kappa_l + mu_l;
	OUTu << "{" << solvenew[0] << ",";
	for (int i = 1; i < n - 1; i++)
	{
		solvenew[i] = solveboof[i - 1];
		OUTu << solvenew[i] << ",";
	}
	solvenew[n - 1] = solveboof[n - 3]*kappa_r + mu_r;
	OUTu << solvenew[n - 1] << "}}" << endl;

	for(int z=0; z<n; z++)
	{
		difference = abs(solvenew[z]-AnalSolve(z*dh,current_t));
		if (difference>maxError)
		{
			maxError = difference;
		}
	}

	current_t = current_t + dt;
	cout << current_t << endl;

	DeletingVector(F);
	DeletingVector(ldiag);
	DeletingVector(diag);
	DeletingVector(udiag);
	DeletingVector(b);
	DeletingVector(solveboof);
	DeletingVector(solve);
	solve = VectCopy(solvenew, n);
	DeletingVector(solvenew);
	//OUTu<<"}";
	DeletingVector(x);
	DeletingVector(t);
	DeletingVector(a);
	DeletingVector(A);
	DeletingVector(B);
	DeletingVector(C);
	//PrintingVector(solve, n);
	DeletingVector(solve);

Errors[count]=maxError;
dt=dt/2;
dh=dh/2;
}//end for определения максимальной ошибки
cout<<"Отношение ошибок: "<<endl;
for (int i=0; i<5; i++)
{
	cout<<Errors[i]/Errors[i+1]<<endl;
}
DeletingVector(Errors);
}

void HeatTransfer_Example2_Quasilinear(datatype sigma)
{
	char FileName[] = "x.txt";
	char FileName2[] = "t.txt";
	char FileName3[] = "2Q.txt";
	ofstream OUTx(FileName);
	ofstream OUTt(FileName2);
	ofstream OUTu(FileName3);

	funct f = P_Heat; //для вычисления коэффициентов a
	funct P = P_Zero; // граничное условие потока
	//предполагаю решать задачу на временном промежутке от 0 до 0.7с  и с длиной стержня 1
	datatype dt = 0.0001;
	datatype dh = 0.05;
	int n = 21; //всего узлов на отрезке от 0 до 1 с шагом h  0---20
	int m = 7001; // временных слоёв на указанной мною сетке   0---1000

	datatype current_t = 0; //переменная текущего времени
	OUTu << "{" << endl;
	OUTu << "{";
	tVector solve = new datatype[n]; //массив точек температур в узлах сетки на текущем временном слое 0---20
	for (int i = 0; i < n - 1; i++)
	{
		solve[i] = 0.04 + i * dh * (1 - i*dh);  // начальные условия варианта 15
		OUTu << solve[i] << ",";
	}
	solve[n - 1] = 0.04;
	OUTu << 0.04 << "}," << endl; //запись первой строчки температур в файл

	OUTx << "{";
	tVector x = new datatype[n]; //массив пространственной сетки 0---20
	for (int i = 0; i < n - 1; i++)
	{
		x[i] = i * dh;
		OUTx << x[i] << ",";
	}
	x[n - 1] = (n - 1)*dh;
	OUTx << x[n - 1] << "}"; //запись простанственной сетки в файл
	OUTx.close();

	OUTt << "{";
	tVector t = new datatype[m];  // массив временной сетки 0---1000
	for (int i = 0; i < m - 1; i++)
	{
		t[i] = i * dt;
		OUTt << t[i] << ",";
	}
	t[m - 1] = (m - 1)*dt;
	OUTt << t[m - 1] << "}"; //запись временнОй сетки в файл
	OUTt.close();

	tVector a = new datatype[n - 1]; //массив маленьких а для вычисления А, В, С (n-1 есть число промежутков интегрирования)1---20(0---19)
	for (int i = 0; i < n - 1; i++)
	{
		//a[i] = dh / Integrate(f, x[i], x[i + 1]); //обрабатываются все отрезочки
		a[i] = 0.5*(f(solve[i+1]) + f(solve[i]));
	}
//	cout << "a:" << endl;
//	PrintingVector(a, n - 1);
//	cin.get();

	//в терминах методички N = 20

	tVector A = new datatype[n - 2];  //1---19(0---18)
	for (int i = 0; i < n - 2; i++)
	{
		A[i] = sigma * a[i] / dh;
	}
//	cout << "A:" << endl;
//	PrintingVector(A, n - 2);
//	cin.get();

	tVector B = new datatype[n - 2];  //1---19(0---18)
	for (int i = 0; i < n - 2; i++)
	{
		B[i] = sigma * a[i + 1] / dh;
	}
//	PrintingVector(B, n - 2);
//	cin.get();

	tVector C = new datatype[n - 2];  //1---19(0---18)
	for (int i = 0; i < n - 2; i++)
	{
		C[i] = A[i] + B[i] + dh / dt;
	}
//	PrintingVector(C, n - 2);
//	cin.get();

	datatype kappa_r = (sigma * a[n - 2] / dh) / (dh / (2 * dt) + sigma * a[n - 2] / dh);
	datatype kappa_l = (sigma*dt*a[0])/(dt*sigma*a[0] + 0.5*dh*dh);
//	cout << kappa;
//	cin.get();

	while (current_t < 0.6999)
	{
		// задал все статичные параметры, осталось задать F и mu. Приступим:

		tVector F = new datatype[n - 2]; //1---19(0---18)
		for (int i = 0; i < n - 2; i++)  //старался сделать без ошибки, но мб она тут есть
		{
			F[i] = (dh / dt) * solve[i + 1] + (1 - sigma) * (a[i + 1] * (solve[i + 2] - solve[i + 1]) - a[i] * (solve[i + 1] - solve[i])) / dh;
		}
		
		datatype mu_r = (solve[n - 1] * dh / (2 * dt) + sigma * P(current_t + dt) + (1 - sigma) * (P(current_t) - a[n - 2] * (solve[n - 1] - solve[n - 2]) / dh))
			/ (dh / (2 * dt) + sigma * a[n - 2] / dh);
		datatype mu_l = (dt*(1-sigma)*a[0]*(solve[1]-solve[0]) + 0.5*dh*dh*solve[0])
				/(0.5*dh*dh + dt*sigma*a[0]);


		//теперь коэффициенты трехдиагональной матрицы:

		tVector ldiag = new datatype[n - 3];  //0---17 (1--18)
		for (int i = 0; i < n - 3; i++)
		{
			ldiag[i] = A[i + 1];
		}

		tVector diag = new datatype[n - 2];  //0---18 (1--19)//записываю как последний ублюдок
		diag[0] = -C[0] + A[0]*kappa_l;
		for (int i = 1; i < n - 3; i++)      //в терминах методички с минусами(похуй)
		{
			diag[i] = -C[i];
		}
		diag[n - 3] = -C[n - 3] + kappa_r * B[n - 3];

		tVector udiag = new datatype[n - 3];//0---17 (1--18)
		for (int i = 0; i < n - 3; i++)
		{
			udiag[i] = B[i];
		}


		//Осталось записать только вектор правой части:

		tVector b = new datatype[n - 2]; // 0---18
		b[0] = -F[0] - mu_l * A[0];
		for (int i = 1; i < n - 3; i++)
		{
			b[i] = -F[i];
		}
		b[n - 3] = -F[n - 3] - mu_r * B[n - 3];

		//PrintingVector(b, n - 2);
		////////////////////////// Все готово для нахождения решения на новом слое.
		////////////////////////// Не уверен только в правильности индексации, но все должно быть ок.
		////////////////////////// Сопоставлял индексы в уме. Ты если будешь править ошибки, то лучше распиши на бумаге
		////////////////////////// В один момент щёлнет, как надо правильно индексировать.
		////////////////////////// я сосал меня ебали
		tVector solveboof = Shuttle(ldiag, diag, udiag, b, n - 2);
		tVector solvenew = new datatype[n]; // вектор решения на новом временном слое

		
		solvenew[0] = solveboof[0]*kappa_l + mu_l;
		OUTu << "{" << solvenew[0] << ",";
		for (int i = 1; i < n - 1; i++)
		{
			solvenew[i] = solveboof[i - 1];
			OUTu << solvenew[i] << ",";
		}
		solvenew[n - 1] = solveboof[n - 3]*kappa_r + mu_r;
		OUTu << solvenew[n - 1] << "}," << endl;
		current_t = current_t + dt;
		cout << current_t << endl;

		DeletingVector(F);
		DeletingVector(ldiag);
		DeletingVector(diag);
		DeletingVector(udiag);
		DeletingVector(b);
		DeletingVector(solveboof);
		DeletingVector(solve);
		solve = VectCopy(solvenew, n);
		for (int i = 0; i < n - 1; i++)
		{
			//a[i] = dh / Integrate(f, x[i], x[i + 1]); //обрабатываются все отрезочки
			a[i] = 0.5*(f(solve[i+1]) + f(solve[i]));
		}
		DeletingVector(solvenew);
		// добавь запись этого вектора в файл. Файл не такой, как pаньше, а другой.
		// Я рисовал на лабе, как должен выглядеть этот файл. Если что, напиши, я напомню.

	}
	tVector F = new datatype[n - 2]; //1---19(0---18)
	for (int i = 0; i < n - 2; i++)  //старался сделать без ошибки, но мб она тут есть
	{
		F[i] = (dh / dt) * solve[i + 1] + (1 - sigma) * (a[i + 1] * (solve[i + 2] - solve[i + 1]) - a[i] * (solve[i + 1] - solve[i])) / dh;
	}

	datatype mu_r = (solve[n - 1] * dh / (2 * dt) + sigma * P(current_t + dt) + (1 - sigma) * (P(current_t) - a[n - 2] * (solve[n - 1] - solve[n - 2]) / dh))
		/ (dh / (2 * dt) + sigma * a[n - 2] / dh);
	datatype mu_l = (dt*(1-sigma)*a[0]*(solve[1]-solve[0]) + 0.5*dh*dh*solve[0])
			/(0.5*dh*dh + dt*sigma*a[0]);


	//теперь коэффициенты трехдиагональной матрицы:

	tVector ldiag = new datatype[n - 3];  //0---17 (1--18)
	for (int i = 0; i < n - 3; i++)
	{
		ldiag[i] = A[i + 1];
	}

	tVector diag = new datatype[n - 2];  //0---18 (1--19)//записываю как последний ублюдок
	diag[0] = -C[0] + A[0]*kappa_l;
	for (int i = 1; i < n - 3; i++)      //в терминах методички с минусами(похуй)
	{
		diag[i] = -C[i];
	}
	diag[n - 3] = -C[n - 3] + kappa_r * B[n - 3];

	tVector udiag = new datatype[n - 3];//0---17 (1--18)
	for (int i = 0; i < n - 3; i++)
	{
		udiag[i] = B[i];
	}


	//Осталось записать только вектор правой части:

	tVector b = new datatype[n - 2]; // 0---18
	b[0] = -F[0] - mu_l * A[0];
	for (int i = 1; i < n - 3; i++)
	{
		b[i] = -F[i];
	}
	b[n - 3] = -F[n - 3] - mu_r * B[n - 3];

	//PrintingVector(b, n - 2);
	////////////////////////// Все готово для нахождения решения на новом слое.
	////////////////////////// Не уверен только в правильности индексации, но все должно быть ок.
	////////////////////////// Сопоставлял индексы в уме. Ты если будешь править ошибки, то лучше распиши на бумаге
	////////////////////////// В один момент щёлнет, как надо правильно индексировать.
	////////////////////////// я сосал меня ебали
	tVector solveboof = Shuttle(ldiag, diag, udiag, b, n - 2);
	tVector solvenew = new datatype[n]; // вектор решения на новом временном слое


	solvenew[0] = solveboof[0]*kappa_l + mu_l;
	OUTu << "{" << solvenew[0] << ",";
	for (int i = 1; i < n - 1; i++)
	{
		solvenew[i] = solveboof[i - 1];
		OUTu << solvenew[i] << ",";
	}
	solvenew[n - 1] = solveboof[n - 3]*kappa_r + mu_r;
	OUTu << solvenew[n - 1] << "}}" << endl;
	current_t = current_t + dt;
	cout << current_t << endl;

	DeletingVector(F);
	DeletingVector(ldiag);
	DeletingVector(diag);
	DeletingVector(udiag);
	DeletingVector(b);
	DeletingVector(solveboof);
	DeletingVector(solve);
	solve = VectCopy(solvenew, n);
	DeletingVector(solvenew);
	//OUTu<<"}";
	DeletingVector(x);
	DeletingVector(t);
	DeletingVector(a);
	DeletingVector(A);
	DeletingVector(B);
	DeletingVector(C);
	//PrintingVector(solve, n);
	DeletingVector(solve);
}

datatype AnalSolve_Example3(datatype x, datatype t)
{
	if (x<=5*t)
	{
		return sqrt(100*t - 20*x);
	}
	else
	{
		return 0;
	}
}

datatype K_Example3(datatype u)
{
	return 0.5*u*u;
}

void HeatTransfer_Example3(datatype sigma, funct2& AnalSolve)
{
		char FileName[] = "x.txt";
		char FileName2[] = "t.txt";
		char FileName3[] = "3.txt";
		ofstream OUTx(FileName);
		ofstream OUTt(FileName2);
		ofstream OUTu(FileName3);
		OUTu.setf(ios::fixed);

		datatype T = 0;
		datatype dt = 0;
		datatype dh=0;
		cout<<"Enter dt step:";
		cin>>dt;
		cout<<endl;
		cout<<"Enter time: ";
		cin>>T;
		cout<<"Enter dh step:";
		cin>>dh;
		cout<<endl;
		int m = T/dt;
		m++;

		//funct f = Kxpr; //для вычисления коэффициентов a
		//funct K=Kx;
		funct K=K_Example3;
		funct Pright = P_Zero; // граничное условие потока
		//funct Pleft = P_Zero;
		//предполагаю решать задачу на временном промежутке от 0 до 0.7с  и с длиной стержня 1
		//datatype dt = 0.0001;
		//datatype dh = 0.05;
		int n = 1./dh +1; //всего узлов на отрезке от 0 до 1 с шагом h  0---20
		//int m = 7001; // временных слоёв на указанной мною сетке   0---1000

		datatype current_t = 0; //переменная текущего времени
		OUTu << "{" << endl;
		OUTu << "{";
		tVector solve = new datatype[n]; //массив точек температур в узлах сетки на текущем временном слое 0---20
		//datatype Energy_i = 0;
		for (int i = 0; i < n - 1; i++)
		{
			solve[i] = 0;  // начальные условия варианта 15
			OUTu << 0 << ",";
		}
		solve[n - 1] = 0;
		OUTu << 0 << "}," << endl; //запись первой строчки температур в файл

		OUTx << "{";
		tVector x = new datatype[n]; //массив пространственной сетки 0---20
		for (int i = 0; i < n - 1; i++)
		{
			x[i] = i * dh;
			OUTx << x[i] << ",";
		}
		x[n - 1] = (n - 1)*dh;
		OUTx << x[n - 1] << "}"; //запись простанственной сетки в файл
		OUTx.close();

		OUTt << "{";
		tVector t = new datatype[m];  // массив временной сетки 0---1000
		for (int i = 0; i < m - 1; i++)
		{
			t[i] = i * dt;
			OUTt << t[i] << ",";
		}
		t[m - 1] = (m - 1)*dt;
		OUTt << t[m - 1] << "}"; //запись временнОй сетки в файл
		OUTt.close();
		m--;

		tVector a = new datatype[n - 1]; //массив маленьких а для вычисления А, В, С (n-1 есть число промежутков интегрирования)1---20(0---19)
		for (int i = 0; i < n - 1; i++)
		{
			//a[i] = dh / Integrate(f, x[i], x[i + 1]); //обрабатываются все отрезочки
			a[i]=0.5*K(solve[i])+K(solve[i+1]);
		}
	//	cout << "a:" << endl;
	//	PrintingVector(a, n - 1);
	//	cin.get();

		//в терминах методички N = 20

		tVector A = new datatype[n - 2];  //1---19(0---18)
		for (int i = 0; i < n - 2; i++)
		{
			A[i] = sigma * a[i] / dh;
		}
	//	cout << "A:" << endl;
	//	PrintingVector(A, n - 2);
	//	cin.get();

		tVector B = new datatype[n - 2];  //1---19(0---18)
		for (int i = 0; i < n - 2; i++)
		{
			B[i] = sigma * a[i + 1] / dh;
		}
	//	PrintingVector(B, n - 2);
	//	cin.get();

		tVector C = new datatype[n - 2];  //1---19(0---18)
		for (int i = 0; i < n - 2; i++)
		{
			C[i] = A[i] + B[i] + dh / dt;
		}

		datatype kappa_r = (sigma * a[n - 2] / dh) / (dh / (2 * dt) + sigma * a[n - 2] / dh);
		datatype mu_r = 0;


		for (int vr=0; vr<m-1; vr++)
		{
			// задал все статичные параметры, осталось задать F и mu. Приступим:

			tVector F = new datatype[n - 2]; //1---19(0---18)
			for (int i = 0; i < n - 2; i++)  //старался сделать без ошибки, но мб она тут есть
			{
				F[i] = (dh / dt) * solve[i + 1] + (1 - sigma) * (a[i + 1] * (solve[i + 2] - solve[i + 1]) - a[i] * (solve[i + 1] - solve[i])) / dh;
			}

			mu_r = (solve[n - 1] * dh / (2 * dt) + sigma * Pright(current_t + dt) + (1 - sigma) * (Pright(current_t) - a[n - 2] * (solve[n - 1] - solve[n - 2]) / dh))
				/ (dh / (2 * dt) + sigma * a[n - 2] / dh);
	//		datatype mu_l = (dt*(1-sigma)*a[0]*(solve[1]-solve[0]) + 0.5*dh*dh*solve[0])
	//				/(0.5*dh*dh + dt*sigma*a[0]);
			//datatype mu_l = (dh*solve[0]/2 + dt*sigma*Pleft(current_t + dt) + a[0]*dt*(1-sigma)/dh*(solve[1]-solve[0]) +dt*(1-sigma)*Pleft(current_t))/(dh/2+dt*sigma*a[0]/dh);


			//теперь коэффициенты трехдиагональной матрицы:

			tVector ldiag = new datatype[n - 3];  //0---17 (1--18)
			for (int i = 0; i < n - 3; i++)
			{
				ldiag[i] = A[i + 1];
			}

			tVector diag = new datatype[n - 2];  //0---18 (1--19)//записываю как последний ублюдок
			//diag[0] = -C[0] + A[0]*kappa_l;
			for (int i = 0; i < n - 3; i++)      //в терминах методички с минусами(похуй)
			{
				diag[i] = -C[i];
			}
			diag[n - 3] = -C[n - 3] + kappa_r * B[n - 3];

			tVector udiag = new datatype[n - 3];//0---17 (1--18)
			for (int i = 0; i < n - 3; i++)
			{
				udiag[i] = B[i];
			}


			//Осталось записать только вектор правой части:

			tVector b = new datatype[n - 2]; // 0---18
			b[0] = -F[0] - AnalSolve(0, current_t) * A[0];
			for (int i = 1; i < n - 3; i++)
			{
				b[i] = -F[i];
			}
			b[n - 3] = -F[n - 3] - mu_r * B[n - 3];

			tVector solveboof = Shuttle(ldiag, diag, udiag, b, n - 2);
			tVector solvenew = new datatype[n]; // вектор решения на новом временном слое


			solvenew[0] = AnalSolve(0, current_t + dt);
			OUTu << "{" << setprecision(20)<< solvenew[0] << ",";
			for (int i = 1; i < n - 1; i++)
			{
				solvenew[i] = solveboof[i - 1];
				OUTu <<setprecision(20)<< solvenew[i] << ",";
			}
			solvenew[n - 1] = solveboof[n - 3]*kappa_r + mu_r;
			OUTu <<setprecision(20)<< solvenew[n - 1] << "}," << endl;
			current_t = current_t + dt;
			cout << current_t << endl;

			DeletingVector(F);
			DeletingVector(ldiag);
			DeletingVector(diag);
			DeletingVector(udiag);
			DeletingVector(b);
			DeletingVector(solveboof);
			DeletingVector(solve);
			solve = VectCopy(solvenew, n);
			for (int i = 0; i < n - 1; i++)
			{
				//a[i] = dh / Integrate(f, x[i], x[i + 1]); //обрабатываются все отрезочки
				a[i] = 0.5*(K(solve[i+1]) + K(solve[i]));
			}
//			cout<<"a:"<<endl;
//			PrintingVector(a, n-1);
			for (int i = 0; i < n - 2; i++)
			{
				A[i] = sigma * a[i] / dh;
			}
			for (int i = 0; i < n - 2; i++)
			{
				B[i] = sigma * a[i + 1] / dh;
			}
			for (int i = 0; i < n - 2; i++)
			{
				C[i] = A[i] + B[i] + dh / dt;
			}
			kappa_r = (sigma * a[n - 2] / dh) / (dh / (2 * dt) + sigma * a[n - 2] / dh);
			mu_r = (solve[n - 1] * dh / (2 * dt) + sigma * Pright(current_t + dt) + (1 - sigma) * (Pright(current_t) - a[n - 2] * (solve[n - 1] - solve[n - 2]) / dh))
							/ (dh / (2 * dt) + sigma * a[n - 2] / dh);
			DeletingVector(solvenew);
			// добавь запись этого вектора в файл. Файл не такой, как pаньше, а другой.
			// Я рисовал на лабе, как должен выглядеть этот файл. Если что, напиши, я напомню.

		}
		tVector F = new datatype[n - 2]; //1---19(0---18)
		for (int i = 0; i < n - 2; i++)  //старался сделать без ошибки, но мб она тут есть
		{
			F[i] = (dh / dt) * solve[i + 1] + (1 - sigma) * (a[i + 1] * (solve[i + 2] - solve[i + 1]) - a[i] * (solve[i + 1] - solve[i])) / dh;
		}

		mu_r = (solve[n - 1] * dh / (2 * dt) + sigma * Pright(current_t + dt) + (1 - sigma) * (Pright(current_t) - a[n - 2] * (solve[n - 1] - solve[n - 2]) / dh))
			/ (dh / (2 * dt) + sigma * a[n - 2] / dh);
		//datatype mu_l = (dt*(1-sigma)*a[0]*(solve[1]-solve[0]) + 0.5*dh*dh*solve[0])
		//		/(0.5*dh*dh + dt*sigma*a[0]);


		//теперь коэффициенты трехдиагональной матрицы:

		tVector ldiag = new datatype[n - 3];  //0---17 (1--18)
		for (int i = 0; i < n - 3; i++)
		{
			ldiag[i] = A[i + 1];
		}

		tVector diag = new datatype[n - 2];  //0---18 (1--19)//записываю как последний ублюдок
		//diag[0] = -C[0] + A[0]*kappa_l;
		for (int i = 0; i < n - 3; i++)      //в терминах методички с минусами(похуй)
		{
			diag[i] = -C[i];
		}
		diag[n - 3] = -C[n - 3] + kappa_r * B[n - 3];

		tVector udiag = new datatype[n - 3];//0---17 (1--18)
		for (int i = 0; i < n - 3; i++)
		{
			udiag[i] = B[i];
		}


		//Осталось записать только вектор правой части:

		tVector b = new datatype[n - 2]; // 0---18
		b[0] = -F[0] - AnalSolve(0, current_t) * A[0];
		for (int i = 1; i < n - 3; i++)
		{
			b[i] = -F[i];
		}
		b[n - 3] = -F[n - 3] - mu_r * B[n - 3];

		//PrintingVector(b, n - 2);
		////////////////////////// Все готово для нахождения решения на новом слое.
		////////////////////////// Не уверен только в правильности индексации, но все должно быть ок.
		////////////////////////// Сопоставлял индексы в уме. Ты если будешь править ошибки, то лучше распиши на бумаге
		////////////////////////// В один момент щёлнет, как надо правильно индексировать.
		tVector solveboof = Shuttle(ldiag, diag, udiag, b, n - 2);
		tVector solvenew = new datatype[n]; // вектор решения на новом временном слое


		solvenew[0] = AnalSolve(0, current_t + dt);
		OUTu << "{" <<setprecision(20)<< solvenew[0] << ",";
		for (int i = 1; i < n - 1; i++)
		{
			solvenew[i] = solveboof[i - 1];
			OUTu <<setprecision(20)<< solvenew[i] << ",";
		}
		solvenew[n - 1] = solveboof[n - 3]*kappa_r + mu_r;
		OUTu <<setprecision(20)<< solvenew[n - 1] << "}}" << endl;
		current_t = current_t + dt;
		cout << current_t << endl;

		DeletingVector(F);
		DeletingVector(ldiag);
		DeletingVector(diag);
		DeletingVector(udiag);
		DeletingVector(b);
		DeletingVector(solveboof);
		DeletingVector(solve);
		solve = VectCopy(solvenew, n);
		DeletingVector(solvenew);
		//OUTu<<"}";
		DeletingVector(x);
		DeletingVector(t);
		DeletingVector(a);
		DeletingVector(A);
		DeletingVector(B);
		DeletingVector(C);
		//PrintingVector(solve, n);
		DeletingVector(solve);
}

void HeatTransfer_Example3_Anal(datatype sigma, funct2& AnalSolve)
{
//		char FileName[] = "x.txt";
//		char FileName2[] = "t.txt";
//		char FileName3[] = "3.txt";
//		ofstream OUTx(FileName);
//		ofstream OUTt(FileName2);
//		ofstream OUTu(FileName3);
//		OUTu.setf(ios::fixed);

		datatype T = 0;
		datatype dt = 0;
		cout<<"Enter dt step:";
		cin>>dt;
		cout<<endl;
		cout<<"Enter time: ";
		cin>>T;
		cout<<endl;
		int m;
		tVector Errors = new datatype[5];
		datatype difference;
		datatype maxError;

		//funct f = Kxpr; //для вычисления коэффициентов a
		//funct K=Kx;
		funct K=K_Example3;
		funct Pright = P_Zero; // граничное условие потока
		//funct Pleft = P_Zero;
		//предполагаю решать задачу на временном промежутке от 0 до 0.7с  и с длиной стержня 1
		//datatype dt = 0.0001;
		datatype dh=0.05;
		int n; //всего узлов на отрезке от 0 до 1 с шагом h  0---20
		//int m = 7001; // временных слоёв на указанной мною сетке   0---1000
for (int count=0; count<5; count++)
{
	difference =0;
	maxError = 0;
	m=T/dt + 1;
	n=1./dh + 1;
		datatype current_t = 0; //переменная текущего времени
//		OUTu << "{" << endl;
//		OUTu << "{";
		tVector solve = new datatype[n]; //массив точек температур в узлах сетки на текущем временном слое 0---20
		//datatype Energy_i = 0;
		for (int i = 0; i < n - 1; i++)
		{
			solve[i] = 0;  // начальные условия варианта 15
			//OUTu << 0 << ",";
		}
		solve[n - 1] = 0;
		//OUTu << 0 << "}," << endl; //запись первой строчки температур в файл

		//OUTx << "{";
		tVector x = new datatype[n]; //массив пространственной сетки 0---20
		for (int i = 0; i < n - 1; i++)
		{
			x[i] = i * dh;
			//OUTx << x[i] << ",";
		}
		x[n - 1] = (n - 1)*dh;
	//	OUTx << x[n - 1] << "}"; //запись простанственной сетки в файл
		//OUTx.close();

		//OUTt << "{";
		tVector t = new datatype[m];  // массив временной сетки 0---1000
		for (int i = 0; i < m - 1; i++)
		{
			t[i] = i * dt;
		//	OUTt << t[i] << ",";
		}
		t[m - 1] = (m - 1)*dt;
		//OUTt << t[m - 1] << "}"; //запись временнОй сетки в файл
		//OUTt.close();
		m--;

		tVector a = new datatype[n - 1]; //массив маленьких а для вычисления А, В, С (n-1 есть число промежутков интегрирования)1---20(0---19)
		for (int i = 0; i < n - 1; i++)
		{
			//a[i] = dh / Integrate(f, x[i], x[i + 1]); //обрабатываются все отрезочки
			a[i]=0.5*K(solve[i])+K(solve[i+1]);
		}
	//	cout << "a:" << endl;
	//	PrintingVector(a, n - 1);
	//	cin.get();

		//в терминах методички N = 20

		tVector A = new datatype[n - 2];  //1---19(0---18)
		for (int i = 0; i < n - 2; i++)
		{
			A[i] = sigma * a[i] / dh;
		}
	//	cout << "A:" << endl;
	//	PrintingVector(A, n - 2);
	//	cin.get();

		tVector B = new datatype[n - 2];  //1---19(0---18)
		for (int i = 0; i < n - 2; i++)
		{
			B[i] = sigma * a[i + 1] / dh;
		}
	//	PrintingVector(B, n - 2);
	//	cin.get();

		tVector C = new datatype[n - 2];  //1---19(0---18)
		for (int i = 0; i < n - 2; i++)
		{
			C[i] = A[i] + B[i] + dh / dt;
		}

		datatype kappa_r = (sigma * a[n - 2] / dh) / (dh / (2 * dt) + sigma * a[n - 2] / dh);
		datatype mu_r = 0;


		for (int vr=0; vr<m-1; vr++)
		{
			// задал все статичные параметры, осталось задать F и mu. Приступим:

			tVector F = new datatype[n - 2]; //1---19(0---18)
			for (int i = 0; i < n - 2; i++)  //старался сделать без ошибки, но мб она тут есть
			{
				F[i] = (dh / dt) * solve[i + 1] + (1 - sigma) * (a[i + 1] * (solve[i + 2] - solve[i + 1]) - a[i] * (solve[i + 1] - solve[i])) / dh;
			}

			mu_r = (solve[n - 1] * dh / (2 * dt) + sigma * Pright(current_t + dt) + (1 - sigma) * (Pright(current_t) - a[n - 2] * (solve[n - 1] - solve[n - 2]) / dh))
				/ (dh / (2 * dt) + sigma * a[n - 2] / dh);
	//		datatype mu_l = (dt*(1-sigma)*a[0]*(solve[1]-solve[0]) + 0.5*dh*dh*solve[0])
	//				/(0.5*dh*dh + dt*sigma*a[0]);
			//datatype mu_l = (dh*solve[0]/2 + dt*sigma*Pleft(current_t + dt) + a[0]*dt*(1-sigma)/dh*(solve[1]-solve[0]) +dt*(1-sigma)*Pleft(current_t))/(dh/2+dt*sigma*a[0]/dh);


			//теперь коэффициенты трехдиагональной матрицы:

			tVector ldiag = new datatype[n - 3];  //0---17 (1--18)
			for (int i = 0; i < n - 3; i++)
			{
				ldiag[i] = A[i + 1];
			}

			tVector diag = new datatype[n - 2];  //0---18 (1--19)//записываю как последний ублюдок
			//diag[0] = -C[0] + A[0]*kappa_l;
			for (int i = 0; i < n - 3; i++)
			{
				diag[i] = -C[i];
			}
			diag[n - 3] = -C[n - 3] + kappa_r * B[n - 3];

			tVector udiag = new datatype[n - 3];//0---17 (1--18)
			for (int i = 0; i < n - 3; i++)
			{
				udiag[i] = B[i];
			}


			//Осталось записать только вектор правой части:

			tVector b = new datatype[n - 2]; // 0---18
			b[0] = -F[0] - AnalSolve(0, current_t) * A[0];
			for (int i = 1; i < n - 3; i++)
			{
				b[i] = -F[i];
			}
			b[n - 3] = -F[n - 3] - mu_r * B[n - 3];

			tVector solveboof = Shuttle(ldiag, diag, udiag, b, n - 2);
			tVector solvenew = new datatype[n]; // вектор решения на новом временном слое


			solvenew[0] = AnalSolve(0, current_t + dt);
			//OUTu << "{" << setprecision(20)<< solvenew[0] << ",";
			for (int i = 1; i < n - 1; i++)
			{
				solvenew[i] = solveboof[i - 1];
			//	OUTu <<setprecision(20)<< solvenew[i] << ",";
			}
			solvenew[n - 1] = solveboof[n - 3]*kappa_r + mu_r;
			//OUTu <<setprecision(20)<< solvenew[n - 1] << "}," << endl;
			for (int i = 1; i < n; i++)
			{
				difference = abs(solvenew[i]-AnalSolve(i*dh, current_t + dt));
				if (difference > maxError) maxError = difference;
			}


			current_t = current_t + dt;
			cout << current_t << endl;

			DeletingVector(F);
			DeletingVector(ldiag);
			DeletingVector(diag);
			DeletingVector(udiag);
			DeletingVector(b);
			DeletingVector(solveboof);
			DeletingVector(solve);
			solve = VectCopy(solvenew, n);
			for (int i = 0; i < n - 1; i++)
			{
				//a[i] = dh / Integrate(f, x[i], x[i + 1]); //обрабатываются все отрезочки
				a[i] = 0.5*(K(solve[i+1]) + K(solve[i]));
			}
//			cout<<"a:"<<endl;
//			PrintingVector(a, n-1);
			for (int i = 0; i < n - 2; i++)
			{
				A[i] = sigma * a[i] / dh;
			}
			for (int i = 0; i < n - 2; i++)
			{
				B[i] = sigma * a[i + 1] / dh;
			}
			for (int i = 0; i < n - 2; i++)
			{
				C[i] = A[i] + B[i] + dh / dt;
			}
			kappa_r = (sigma * a[n - 2] / dh) / (dh / (2 * dt) + sigma * a[n - 2] / dh);
			mu_r = (solve[n - 1] * dh / (2 * dt) + sigma * Pright(current_t + dt) + (1 - sigma) * (Pright(current_t) - a[n - 2] * (solve[n - 1] - solve[n - 2]) / dh))
							/ (dh / (2 * dt) + sigma * a[n - 2] / dh);
			DeletingVector(solvenew);
			// добавь запись этого вектора в файл. Файл не такой, как pаньше, а другой.
			// Я рисовал на лабе, как должен выглядеть этот файл. Если что, напиши, я напомню.

		}


		tVector F = new datatype[n - 2]; //1---19(0---18)
		for (int i = 0; i < n - 2; i++)  //старался сделать без ошибки, но мб она тут есть
		{
			F[i] = (dh / dt) * solve[i + 1] + (1 - sigma) * (a[i + 1] * (solve[i + 2] - solve[i + 1]) - a[i] * (solve[i + 1] - solve[i])) / dh;
		}

		mu_r = (solve[n - 1] * dh / (2 * dt) + sigma * Pright(current_t + dt) + (1 - sigma) * (Pright(current_t) - a[n - 2] * (solve[n - 1] - solve[n - 2]) / dh))
			/ (dh / (2 * dt) + sigma * a[n - 2] / dh);
		//datatype mu_l = (dt*(1-sigma)*a[0]*(solve[1]-solve[0]) + 0.5*dh*dh*solve[0])
		//		/(0.5*dh*dh + dt*sigma*a[0]);


		//теперь коэффициенты трехдиагональной матрицы:

		tVector ldiag = new datatype[n - 3];  //0---17 (1--18)
		for (int i = 0; i < n - 3; i++)
		{
			ldiag[i] = A[i + 1];
		}

		tVector diag = new datatype[n - 2];  //0---18 (1--19)//записываю как последний ублюдок
		//diag[0] = -C[0] + A[0]*kappa_l;
		for (int i = 0; i < n - 3; i++)      //в терминах методички с минусами(похуй)
		{
			diag[i] = -C[i];
		}
		diag[n - 3] = -C[n - 3] + kappa_r * B[n - 3];

		tVector udiag = new datatype[n - 3];//0---17 (1--18)
		for (int i = 0; i < n - 3; i++)
		{
			udiag[i] = B[i];
		}


		//Осталось записать только вектор правой части:

		tVector b = new datatype[n - 2]; // 0---18
		b[0] = -F[0] - AnalSolve(0, current_t) * A[0];
		for (int i = 1; i < n - 3; i++)
		{
			b[i] = -F[i];
		}
		b[n - 3] = -F[n - 3] - mu_r * B[n - 3];

		//PrintingVector(b, n - 2);
		////////////////////////// Все готово для нахождения решения на новом слое.
		////////////////////////// Не уверен только в правильности индексации, но все должно быть ок.
		////////////////////////// Сопоставлял индексы в уме. Ты если будешь править ошибки, то лучше распиши на бумаге
		////////////////////////// В один момент щёлнет, как надо правильно индексировать.
		////////////////////////// я сосал меня ебали
		tVector solveboof = Shuttle(ldiag, diag, udiag, b, n - 2);
		tVector solvenew = new datatype[n]; // вектор решения на новом временном слое


		solvenew[0] = AnalSolve(0, current_t + dt);
		//OUTu << "{" <<setprecision(20)<< solvenew[0] << ",";
		for (int i = 1; i < n - 1; i++)
		{
			solvenew[i] = solveboof[i - 1];
		//	OUTu <<setprecision(20)<< solvenew[i] << ",";
		}
		solvenew[n - 1] = solveboof[n - 3]*kappa_r + mu_r;
		//OUTu <<setprecision(20)<< solvenew[n - 1] << "}}" << endl;
		for (int i = 1; i < n; i++)
		{
			difference = abs(solvenew[i]-AnalSolve(i*dh, current_t + dt));
			if (difference > maxError) maxError = difference;
		}
		current_t = current_t + dt;
		cout << current_t << endl;

		DeletingVector(F);
		DeletingVector(ldiag);
		DeletingVector(diag);
		DeletingVector(udiag);
		DeletingVector(b);
		DeletingVector(solveboof);
		DeletingVector(solve);
		solve = VectCopy(solvenew, n);
		DeletingVector(solvenew);
		//OUTu<<"}";
		DeletingVector(x);
		DeletingVector(t);
		DeletingVector(a);
		DeletingVector(A);
		DeletingVector(B);
		DeletingVector(C);
		//PrintingVector(solve, n);
		DeletingVector(solve);
Errors[count]=maxError;
dh=dh/2;
dt=dt/4;
}
cout<<"Погрешность на шаге:         Отношение ошибок: "<<endl;
for (int i=0; i<4; i++)
{
	cout<<Errors[i]<<"                   "<<Errors[i]/Errors[i+1]<<endl;
}
cout<<Errors[4]<<"                     ---"<<endl;
DeletingVector(Errors);
}


void HeatTransfer_Example3_Iter(datatype sigma, funct2 &AnalSolve)
{
	char FileName[] = "x.txt";
	char FileName2[] = "t.txt";
	char FileName3[] = "3_Iter.txt";
	ofstream OUTx(FileName);
	ofstream OUTt(FileName2);
	ofstream OUTu(FileName3);
	OUTu.setf(ios::fixed);
	datatype T, dt, dh;
	cout<< "Enter Time: ";
	cin>>T;
	cout<<endl;
	cout<< "Enter dt step: ";
	cin>>dt;
	cout<<endl;
	cout<< "Enter dh step: ";
	cin>>dh;
	cout<<endl;
	int n = 1./dh +1; //число узлов
	int m = T/dt; //число шагов по времени
	funct K=K_Example3;
	funct Pright= P_Zero;
	datatype current_t = 0;

	OUTx<<"{";
	for (int i=0; i<n-1; i++)
	{
		OUTx<<i*dh<<",";
	}
	OUTx<<(n-1)*dh<<"}";
	OUTx.close();

	OUTt<<"{";
	for (int i=0; i<m; i++)
	{
		OUTt<<i*dt<<",";
	}
	OUTt<<m*dt<<"}";
	OUTt.close();

	tVector solve=new datatype[n];
	OUTu<<"{{";
	for(int i=0; i<n-1; i++)
	{
		solve[i]=0;
		OUTu<< 0 <<",";
	}
	solve[n-1]=0;
	OUTu<< 0 <<"},";

	tVector A = new datatype[n - 2];

	tVector B = new datatype[n - 2];

	tVector C = new datatype[n - 2];

	tVector F = new datatype[n - 2];

	tVector a = new datatype[n - 1];

	tVector ud = new datatype[n - 3];

	tVector d = new datatype[n - 2];

	tVector ld = new datatype[n - 3];

	tVector b = new datatype[n - 2];

	for(int time=0; time<m-1; time++)
	{
		datatype norm =100;
		tVector constantU=VectCopy(solve, n);
		while(norm > eps5)
		{
			for(int i=0; i<n-1; i++)
			{
				a[i]=0.5*(K(solve[i])+K(solve[i+1]));
			}
			for(int i=0; i<n-2; i++)
			{
				A[i]=sigma*a[i]/dh;
			}
			for(int i=0; i<n-2; i++)
			{
				B[i]=sigma*a[i+1]/dh;
			}
			for(int i=0; i<n-2; i++)
			{
				C[i]=A[i]+B[i]+dh/dt;
			}
			for(int i=0; i<n-2; i++)
			{
				F[i]=dh*constantU[i+1]/dt + (1-sigma)*(a[i+1]*(constantU[i+2]-constantU[i+1])/dh - a[i]*(constantU[i+1]-constantU[i])/dh);
			}
			datatype mu_r = (constantU[n-1]*dh/(2*dt) + (sigma-1)*a[n-2]*(constantU[n-1]-constantU[n-2])/dh)/(dh/(2*dt)+sigma*a[n-2]/dh);
			datatype kappa_r = (sigma*a[n-2]/dh)/(dh/(2*dt) + sigma*a[n-2]/dh);

			for (int i=0; i<n-3; i++)
			{
				d[i]=-C[i];
			}
			d[n-3]=-C[n-3]+kappa_r*B[n-3];

			for(int i=0; i<n-3; i++)
			{
				ld[i]=A[i+1];
			}

			for(int i=0; i<n-3; i++)
			{
				ud[i]=B[i];
			}

			b[0]=-F[0]-A[0]*AnalSolve(0, current_t+dt); //правильно ли тут записываю? это же s-ый шаг
			for(int i=1; i<n-3; i++)
			{
				b[i]=-F[i];
			}
			b[n-3]=-F[n-3]-mu_r*B[n-3];

			tVector solveboof = Shuttle(ld, d, ud, b, n-2);
			tVector solvenew = new datatype[n];

			solvenew[0]= AnalSolve(0, current_t + dt);
			for(int i=1; i<n-1;i++)
			{
				solvenew[i]=solveboof[i-1];
			}
			solvenew[n-1]= kappa_r*solveboof[n-2]+mu_r;

			tVector DX = Difference(solve, solvenew, n);
			norm = NormOfVectorC(DX, n);
			DeletingVector(solveboof);
			DeletingVector(solve);
			solve=VectCopy(solvenew, n);
			DeletingVector(solvenew);
		}
		OUTu<<"{";
		for(int i=0; i<n-1; i++)
		{
			OUTu<<setprecision(17)<<solve[i]<<",";
		}
		OUTu<<solve[n-1]<<"},";
		current_t=current_t+dt;
	}
	datatype norm =100;
	tVector constantU=VectCopy(solve, n);
	while(norm > eps5)
	{
		for(int i=0; i<n-1; i++)
		{
			a[i]=0.5*(K(solve[i])+K(solve[i+1]));
		}
		for(int i=0; i<n-2; i++)
		{
			A[i]=sigma*a[i]/dh;
		}
		for(int i=0; i<n-2; i++)
		{
			B[i]=sigma*a[i+1]/dh;
		}
		for(int i=0; i<n-2; i++)
		{
			C[i]=A[i]+B[i]+dh/dt;
		}
		for(int i=0; i<n-2; i++)
		{
			F[i]=dh*constantU[i+1]/dt + (1-sigma)*(a[i+1]*(constantU[i+2]-constantU[i+1])/dh - a[i]*(constantU[i+1]-constantU[i])/dh);
		}
		datatype mu_r = (constantU[n-1]*dh/(2*dt) + (sigma-1)*a[n-2]*(constantU[n-1]-constantU[n-2])/dh)/(dh/(2*dt)+sigma*a[n-2]/dh);
		datatype kappa_r = (sigma*a[n-2]/dh)/(dh/(2*dt) + sigma*a[n-2]/dh);

		for (int i=0; i<n-3; i++)
		{
			d[i]=-C[i];
		}
		d[n-3]=-C[n-3]+kappa_r*B[n-3];

		for(int i=0; i<n-3; i++)
		{
			ld[i]=A[i+1];
		}

		for(int i=0; i<n-3; i++)
		{
			ud[i]=B[i];
		}

		b[0]=-F[0]-A[0]*AnalSolve(0, current_t+dt); //правильно ли тут записываю? это же s-ый шаг
		for(int i=1; i<n-3; i++)
		{
			b[i]=-F[i];
		}
		b[n-3]=-F[n-3]-mu_r*B[n-3];

		tVector solveboof = Shuttle(ld, d, ud, b, n-2);
		tVector solvenew = new datatype[n];

		solvenew[0]= AnalSolve(0, current_t + dt);
		for(int i=1; i<n-1;i++)
		{
			solvenew[i]=solveboof[i-1];
		}
		solvenew[n-1]= kappa_r*solveboof[n-2]+mu_r;

		tVector DX = Difference(solve, solvenew, n);
		norm = NormOfVectorC(DX, n);
		DeletingVector(solveboof);
		DeletingVector(solve);
		solve=VectCopy(solvenew, n);
		DeletingVector(solvenew);
	}
	OUTu<<"{";
	for(int i=0; i<n-1; i++)
	{
		OUTu<<setprecision(17)<<solve[i]<<",";
	}
	OUTu<<solve[n-1]<<"}}";
	OUTu.close();
	DeletingVector(solve);
}

void HeatTransfer_Example3_Iter_Anal(datatype sigma, funct2 &AnalSolve)
{
//	char FileName[] = "x.txt";
//	char FileName2[] = "t.txt";
//	char FileName3[] = "3_Iter.txt";
//	ofstream OUTx(FileName);
//	ofstream OUTt(FileName2);
//	ofstream OUTu(FileName3);
//	OUTu.setf(ios::fixed);
	datatype T, dt, dh;
	cout<< "Enter Time: ";
	cin>>T;
	cout<<endl;
	cout<< "Enter dt step: ";
	cin>>dt;
	cout<<endl;
	cout<< "Enter dh step: ";
	cin>>dh;
	cout<<endl;
	tVector Error=new datatype[5];
for(int count=0; count<5; count++)
{
	int n = 1./dh +1; //число узлов
	int m = T/dt; //число шагов по времени
	funct K=K_Example3;
	funct Pright= P_Zero;
	datatype current_t = 0;

//	OUTx<<"{";
//	for (int i=0; i<n-1; i++)
//	{
//		OUTx<<i*dh<<",";
//	}
//	OUTx<<(n-1)*dh<<"}";
//	OUTx.close();

//	OUTt<<"{";
//	for (int i=0; i<m; i++)
//	{
//		OUTt<<i*dt<<",";
//	}
//	OUTt<<m*dt<<"}";
//	OUTt.close();

	tVector solve=new datatype[n];
	//OUTu<<"{{";
	for(int i=0; i<n-1; i++)
	{
		solve[i]=0;
		//OUTu<< 0 <<",";
	}
	solve[n-1]=0;
	//OUTu<< 0 <<"},";

	tVector A = new datatype[n - 2];

	tVector B = new datatype[n - 2];

	tVector C = new datatype[n - 2];

	tVector F = new datatype[n - 2];

	tVector a = new datatype[n - 1];

	tVector ud = new datatype[n - 3];

	tVector d = new datatype[n - 2];

	tVector ld = new datatype[n - 3];

	tVector b = new datatype[n - 2];

	for(int time=0; time<m-1; time++)
	{
		datatype norm =100;
		tVector constantU=VectCopy(solve, n);
		while(norm > eps5)
		{
			for(int i=0; i<n-1; i++)
			{
				a[i]=0.5*(K(solve[i])+K(solve[i+1]));
			}
			for(int i=0; i<n-2; i++)
			{
				A[i]=sigma*a[i]/dh;
			}
			for(int i=0; i<n-2; i++)
			{
				B[i]=sigma*a[i+1]/dh;
			}
			for(int i=0; i<n-2; i++)
			{
				C[i]=A[i]+B[i]+dh/dt;
			}
			for(int i=0; i<n-2; i++)
			{
				F[i]=dh*constantU[i+1]/dt + (1-sigma)*(a[i+1]*(constantU[i+2]-constantU[i+1])/dh - a[i]*(constantU[i+1]-constantU[i])/dh);
			}
			datatype mu_r = (constantU[n-1]*dh/(2*dt) + (sigma-1)*a[n-2]*(constantU[n-1]-constantU[n-2])/dh)/(dh/(2*dt)+sigma*a[n-2]/dh);
			datatype kappa_r = (sigma*a[n-2]/dh)/(dh/(2*dt) + sigma*a[n-2]/dh);

			for (int i=0; i<n-3; i++)
			{
				d[i]=-C[i];
			}
			d[n-3]=-C[n-3]+kappa_r*B[n-3];

			for(int i=0; i<n-3; i++)
			{
				ld[i]=A[i+1];
			}

			for(int i=0; i<n-3; i++)
			{
				ud[i]=B[i];
			}

			b[0]=-F[0]-A[0]*AnalSolve(0, current_t+dt); //правильно ли тут записываю? это же s-ый шаг
			for(int i=1; i<n-3; i++)
			{
				b[i]=-F[i];
			}
			b[n-3]=-F[n-3]-mu_r*B[n-3];

			tVector solveboof = Shuttle(ld, d, ud, b, n-2);
			tVector solvenew = new datatype[n];

			solvenew[0]= AnalSolve(0, current_t + dt);
			for(int i=1; i<n-1;i++)
			{
				solvenew[i]=solveboof[i-1];
			}
			solvenew[n-1]= kappa_r*solveboof[n-2]+mu_r;

			tVector DX = Difference(solve, solvenew, n);
			norm = NormOfVectorC(DX, n);
			DeletingVector(solveboof);
			DeletingVector(solve);
			solve=VectCopy(solvenew, n);
			DeletingVector(solvenew);
		}
		current_t=current_t+dt;
	}
	datatype norm =100;
	tVector constantU=VectCopy(solve, n);
	while(norm > eps5)
	{
		for(int i=0; i<n-1; i++)
		{
			a[i]=0.5*(K(solve[i])+K(solve[i+1]));
		}
		for(int i=0; i<n-2; i++)
		{
			A[i]=sigma*a[i]/dh;
		}
		for(int i=0; i<n-2; i++)
		{
			B[i]=sigma*a[i+1]/dh;
		}
		for(int i=0; i<n-2; i++)
		{
			C[i]=A[i]+B[i]+dh/dt;
		}
		for(int i=0; i<n-2; i++)
		{
			F[i]=dh*constantU[i+1]/dt + (1-sigma)*(a[i+1]*(constantU[i+2]-constantU[i+1])/dh - a[i]*(constantU[i+1]-constantU[i])/dh);
		}
		datatype mu_r = (constantU[n-1]*dh/(2*dt) + (sigma-1)*a[n-2]*(constantU[n-1]-constantU[n-2])/dh)/(dh/(2*dt)+sigma*a[n-2]/dh);
		datatype kappa_r = (sigma*a[n-2]/dh)/(dh/(2*dt) + sigma*a[n-2]/dh);

		for (int i=0; i<n-3; i++)
		{
			d[i]=-C[i];
		}
		d[n-3]=-C[n-3]+kappa_r*B[n-3];

		for(int i=0; i<n-3; i++)
		{
			ld[i]=A[i+1];
		}

		for(int i=0; i<n-3; i++)
		{
			ud[i]=B[i];
		}

		b[0]=-F[0]-A[0]*AnalSolve(0, current_t+dt); //правильно ли тут записываю? это же s-ый шаг
		for(int i=1; i<n-3; i++)
		{
			b[i]=-F[i];
		}
		b[n-3]=-F[n-3]-mu_r*B[n-3];

		tVector solveboof = Shuttle(ld, d, ud, b, n-2);
		tVector solvenew = new datatype[n];

		solvenew[0]= AnalSolve(0, current_t + dt);
		for(int i=1; i<n-1;i++)
		{
			solvenew[i]=solveboof[i-1];
		}
		solvenew[n-1]= kappa_r*solveboof[n-2]+mu_r;

		tVector DX = Difference(solve, solvenew, n);
		norm = NormOfVectorC(DX, n);
		DeletingVector(solveboof);
		DeletingVector(solve);
		solve=VectCopy(solvenew, n);
		DeletingVector(solvenew);
	}
	Error[count]=abs(solve[(n-1)/2]-AnalSolve(dh*(n-1)/2, current_t+dt));
	DeletingVector(solve);
	dt=dt/4;
	dh=dh/2;
}
cout<<"Погрешность на шаге:         Отношение ошибок: "<<endl;
for (int i=0; i<4; i++)
{
	cout<<Error[i]<<"                   "<<Error[i]/Error[i+1]<<endl;
}
cout<<Error[4]<<"                     ---"<<endl;
DeletingVector(Error);
}
//void HeatTransfer_Example3_Iter(datatype sigma, funct2& AnalSolve)
//{
//		char FileName[] = "x.txt";
//		char FileName2[] = "t.txt";
//		char FileName3[] = "3_Iter.txt";
//		ofstream OUTx(FileName);
//		ofstream OUTt(FileName2);
//		ofstream OUTu(FileName3);
//		OUTu.setf(ios::fixed);
//
//		datatype T = 0;
//		datatype dt = 0;
//		datatype dh=0;
//		cout<<"Enter dt step:";
//		cin>>dt;
//		cout<<endl;
//		cout<<"Enter time: ";
//		cin>>T;
//		cout<<endl;
//		cout<<"Enter dh step:";
//		cin>>dh;
//		cout<<endl;
//		int m;
////		tVector Errors = new datatype[5];
////		datatype difference;
////		datatype maxError;
//
//		//funct f = Kxpr; //для вычисления коэффициентов a
//		//funct K=Kx;
//		funct K=K_Example3;
//		funct Pright = P_Zero; // граничное условие потока
//		//funct Pleft = P_Zero;
//		//предполагаю решать задачу на временном промежутке от 0 до 0.7с  и с длиной стержня 1
//		//datatype dt = 0.0001;
//		//datatype dh=0.05;
//		int n; //всего узлов на отрезке от 0 до 1 с шагом h  0---20
//		//int m = 7001; // временных слоёв на указанной мною сетке   0---1000
////for (int count=0; count<5; count++)
////{
////	difference =0;
////	maxError = 0;
//	m=T/dt + 1;
//	n=1./dh + 1;
//		datatype current_t = 0; //переменная текущего времени
//		OUTu << "{" << endl;
//		OUTu << "{";
//		tVector solve = new datatype[n]; //массив точек температур в узлах сетки на текущем временном слое 0---20
//		//datatype Energy_i = 0;
//		for (int i = 0; i < n - 1; i++)
//		{
//			solve[i] = 0;  // начальные условия варианта 15
//			OUTu << 0 << ",";
//		}
//		solve[n - 1] = 0;
//		OUTu << 0 << "}," << endl; //запись первой строчки температур в файл
//
//		OUTx << "{";
//		tVector x = new datatype[n]; //массив пространственной сетки 0---20
//		for (int i = 0; i < n - 1; i++)
//		{
//			x[i] = i * dh;
//			OUTx << x[i] << ",";
//		}
//		x[n - 1] = (n - 1)*dh;
//		OUTx << x[n - 1] << "}"; //запись простанственной сетки в файл
//		OUTx.close();
//
//		OUTt << "{";
//		tVector t = new datatype[m];  // массив временной сетки 0---1000
//		for (int i = 0; i < m - 1; i++)
//		{
//			t[i] = i * dt;
//			OUTt << t[i] << ",";
//		}
//		t[m - 1] = (m - 1)*dt;
//		OUTt << t[m - 1] << "}"; //запись временнОй сетки в файл
//		OUTt.close();
//		m--;
//
//		tVector A = new datatype[n - 2];  //1---19(0---18)
//
//		tVector B = new datatype[n - 2];  //1---19(0---18)
//
//		tVector C = new datatype[n - 2];  //1---19(0---18)
//
//		tVector a = new datatype[n - 1];
//
//
//	for (int vr=0; vr<m-1; vr++)
//	{
//			datatype norm =100;
//			tVector constantU=VectCopy(solve, n);
//			int counter=0;
//			while(norm>eps5)
//		{
//
//			for (int i = 0; i < n - 1; i++)
//			{
//
//				a[i]=0.5*K(solve[i])+K(solve[i+1]);
//			}
//
//			for (int i = 0; i < n - 2; i++)
//			{
//				A[i] = sigma * a[i] / dh;
//			}
//
//
//			for (int i = 0; i < n - 2; i++)
//			{
//				B[i] = sigma * a[i + 1] / dh;
//			}
//
//			for (int i = 0; i < n - 2; i++)
//			{
//				C[i] = A[i] + B[i] + dh / dt;
//			}
//
//			datatype kappa_r = (sigma * a[n - 2] / dh) / (dh / (2 * dt) + sigma * a[n - 2] / dh);
//			datatype mu_r = 0;
//			// задал все статичные параметры, осталось задать F и mu. Приступим:
//
//			tVector F = new datatype[n - 2]; //1---19(0---18)
//			for (int i = 0; i < n - 2; i++)  //старался сделать без ошибки, но мб она тут есть
//			{
//				F[i] = (dh / dt) * constantU[i + 1] + (1 - sigma) * (a[i + 1] * (constantU[i + 2] - constantU[i + 1]) - a[i] * (constantU[i + 1] - constantU[i])) / dh;
//			}
//
//			mu_r = (constantU[n - 1] * dh / (2 * dt) + sigma * Pright(current_t + dt) + (1 - sigma) * (Pright(current_t) - a[n - 2] * (constantU[n - 1] - constantU[n - 2]) / dh))
//				/ (dh / (2 * dt) + sigma * a[n - 2] / dh);
//
//			//теперь коэффициенты трехдиагональной матрицы:
//
//			tVector ldiag = new datatype[n - 3];  //0---17 (1--18)
//			for (int i = 0; i < n - 3; i++)
//			{
//				ldiag[i] = A[i + 1];
//			}
//
//			tVector diag = new datatype[n - 2];  //0---18 (1--19)//записываю как последний ублюдок
//			//diag[0] = -C[0] + A[0]*kappa_l;
//			for (int i = 0; i < n - 3; i++)      //в терминах методички с минусами(похуй)
//			{
//				diag[i] = -C[i];
//			}
//			diag[n - 3] = -C[n - 3] + kappa_r * B[n - 3];
//
//			tVector udiag = new datatype[n - 3];//0---17 (1--18)
//			for (int i = 0; i < n - 3; i++)
//			{
//				udiag[i] = B[i];
//			}
//
//
//			//Осталось записать только вектор правой части:
//
//			tVector b = new datatype[n - 2]; // 0---18
//			b[0] = -F[0] - AnalSolve(0, current_t) * A[0];
//			for (int i = 1; i < n - 3; i++)
//			{
//				b[i] = -F[i];
//			}
//			b[n - 3] = -F[n - 3] - mu_r * B[n - 3];
//
//			tVector solveboof = Shuttle(ldiag, diag, udiag, b, n - 2);
//			tVector solvenew = new datatype[n]; // вектор решения на новом временном слое
//
//
//			solvenew[0] = AnalSolve(0, current_t + dt);
//			for (int i = 1; i < n - 1; i++)
//			{
//				solvenew[i] = solveboof[i - 1];
//			}
//			solvenew[n - 1] = solveboof[n - 3]*kappa_r + mu_r;
//
//			tVector DX=Difference(solvenew, solve, n);
//			norm=NormOfVectorC(DX,n);
//			DeletingVector(DX);
//
//			DeletingVector(F);
//			DeletingVector(ldiag);
//			DeletingVector(diag);
//			DeletingVector(udiag);
//			DeletingVector(b);
//			DeletingVector(solveboof);
//			DeletingVector(solve);
//			solve = VectCopy(solvenew, n);
//			DeletingVector(solvenew);
//			counter++;
//		}
//			cout<<"Итераций до сходимости: "<<counter<<endl;
//			OUTu << "{" << solve[0] << ",";
//			for (int i = 1; i < n - 1; i++)
//				{
//					OUTu << solve[i] << ",";
//				}
//			OUTu << solve[n - 1] << "}," << endl;
//			current_t = current_t + dt;
//			DeletingVector(constantU);
//	}
/////////////////////////////////////////////////
//	datatype norm =100;
//		tVector constantU=VectCopy(solve, n);
//		int counter=0;
//		while(norm>eps5)
//	{
//
//		for (int i = 0; i < n - 1; i++)
//		{
//
//			a[i]=0.5*K(solve[i])+K(solve[i+1]);
//		}
//
//		for (int i = 0; i < n - 2; i++)
//		{
//			A[i] = sigma * a[i] / dh;
//		}
//
//
//		for (int i = 0; i < n - 2; i++)
//		{
//			B[i] = sigma * a[i + 1] / dh;
//		}
//
//		for (int i = 0; i < n - 2; i++)
//		{
//			C[i] = A[i] + B[i] + dh / dt;
//		}
//
//		datatype kappa_r = (sigma * a[n - 2] / dh) / (dh / (2 * dt) + sigma * a[n - 2] / dh);
//		datatype mu_r = 0;
//		// задал все статичные параметры, осталось задать F и mu. Приступим:
//
//		tVector F = new datatype[n - 2]; //1---19(0---18)
//		for (int i = 0; i < n - 2; i++)  //старался сделать без ошибки, но мб она тут есть
//		{
//			F[i] = (dh / dt) * constantU[i + 1] + (1 - sigma) * (a[i + 1] * (constantU[i + 2] - constantU[i + 1]) - a[i] * (constantU[i + 1] - constantU[i])) / dh;
//		}
//
//		mu_r = (constantU[n - 1] * dh / (2 * dt) + sigma * Pright(current_t + dt) + (1 - sigma) * (Pright(current_t) - a[n - 2] * (constantU[n - 1] - constantU[n - 2]) / dh))
//			/ (dh / (2 * dt) + sigma * a[n - 2] / dh);
////		datatype mu_l = (dt*(1-sigma)*a[0]*(solve[1]-solve[0]) + 0.5*dh*dh*solve[0])
////				/(0.5*dh*dh + dt*sigma*a[0]);
//		//datatype mu_l = (dh*solve[0]/2 + dt*sigma*Pleft(current_t + dt) + a[0]*dt*(1-sigma)/dh*(solve[1]-solve[0]) +dt*(1-sigma)*Pleft(current_t))/(dh/2+dt*sigma*a[0]/dh);
//
//
//		//теперь коэффициенты трехдиагональной матрицы:
//
//		tVector ldiag = new datatype[n - 3];  //0---17 (1--18)
//		for (int i = 0; i < n - 3; i++)
//		{
//			ldiag[i] = A[i + 1];
//		}
//
//		tVector diag = new datatype[n - 2];  //0---18 (1--19)//записываю как последний ублюдок
//		//diag[0] = -C[0] + A[0]*kappa_l;
//		for (int i = 0; i < n - 3; i++)      //в терминах методички с минусами(похуй)
//		{
//			diag[i] = -C[i];
//		}
//		diag[n - 3] = -C[n - 3] + kappa_r * B[n - 3];
//
//		tVector udiag = new datatype[n - 3];//0---17 (1--18)
//		for (int i = 0; i < n - 3; i++)
//		{
//			udiag[i] = B[i];
//		}
//
//
//		//Осталось записать только вектор правой части:
//
//		tVector b = new datatype[n - 2]; // 0---18
//		b[0] = -F[0] - AnalSolve(0, current_t) * A[0];
//		for (int i = 1; i < n - 3; i++)
//		{
//			b[i] = -F[i];
//		}
//		b[n - 3] = -F[n - 3] - mu_r * B[n - 3];
//
//		tVector solveboof = Shuttle(ldiag, diag, udiag, b, n - 2);
//		tVector solvenew = new datatype[n]; // вектор решения на новом временном слое
//
//
//		solvenew[0] = AnalSolve(0, current_t + dt);
//		for (int i = 1; i < n - 1; i++)
//		{
//			solvenew[i] = solveboof[i - 1];
//		}
//		solvenew[n - 1] = solveboof[n - 3]*kappa_r + mu_r;
//
//		tVector DX=Difference(solvenew, solve, n);
//		norm=NormOfVectorC(DX,n);
//		DeletingVector(DX);
//
//		DeletingVector(F);
//		DeletingVector(ldiag);
//		DeletingVector(diag);
//		DeletingVector(udiag);
//		DeletingVector(b);
//		DeletingVector(solveboof);
//		DeletingVector(solve);
//		solve = VectCopy(solvenew, n);
//		DeletingVector(solvenew);
//		counter++;
//	}
//		cout<<"Итераций до сходимости: "<<counter<<endl;
//		OUTu << "{" << solve[0] << ",";
//		for (int i = 1; i < n - 1; i++)
//			{
//				OUTu << solve[i] << ",";
//			}
//		OUTu << solve[n - 1] << "}}" << endl;
//		OUTu.close();
//		DeletingVector(constantU);
//		DeletingVector(solve);
//}


void Wave(funct &f, funct &g, funct &phi, funct &psi)//соответственно начальные условия задаем как входные параметры
{
	// задам х от 0 до 1(почему бы и нет), возьму наверное 50 кусочков, то есть 51 точку 
	// вторую производную функции f для аппроксимации второго слоя буду считать численно
	// также будем считать а=1
	char FileName[] = "x_Wave.txt";
	char FileName2[] = "t_Wave.txt";
	char FileName3[] = "u_Wave.txt";
	ofstream OUTx(FileName);
	ofstream OUTt(FileName2);
	ofstream OUTu(FileName3);

	//хуярим в файлы в нормальном виде без ебаных е кстати
	OUTx.setf(ios::fixed);
	OUTt.setf(ios::fixed);
	OUTu.setf(ios::fixed);
	
	datatype dt = 0.01;
	datatype dh = 0.01;//условие Куранта |a*tao/h|<=1 выполняется
	int n = 101; //всего узлов на отрезке от 0 до 1 с шагом dh  0---5000
	int m = 3001; // временных слоёв на указанной мною сетке   0---7000

	//хуярим на 7 секунд ебать
	
	OUTu << "{" << endl;
	OUTu << "{";
	tVector solve0 = new datatype[n]; //массив точек точек струны на предыдущем слое.
	tVector solve1 = new datatype[n]; //массив точек точек струны на текущем временном слое.
	tVector solve2 = new datatype[n]; //массив точек точек струны наследующем временном слое.
	for (int i = 0; i < n - 1; i++)
	{
		solve0[i] = f( i * dh);  // начальные условия, форма струны при t=0
		OUTu << solve0[i] << ",";
	}
	solve0[n - 1] = f((n-1) * dh);
	OUTu << solve0[n - 1] << "}," << endl; //запись первой строчки формы в файл

	OUTx << "{";
	tVector x = new datatype[n]; //массив пространственной сетки 0---5000
	for (int i = 0; i < n - 1; i++)
	{
		x[i] = i * dh;
		OUTx << x[i] << ",";
	}
	x[n - 1] = (n - 1)*dh;
	OUTx << x[n - 1] << "}"; //запись простанственной сетки в файл
	OUTx.close();

	OUTt << "{";
	tVector t = new datatype[m];  // массив временной сетки 0---7000
	for (int i = 0; i < m - 1; i++)
	{
		t[i] = i * dt;
		OUTt << t[i] << ",";
	}
	t[m - 1] = (m - 1)*dt;
	OUTt << t[m - 1] << "}"; //запись временнОй сетки в файл
	OUTt.close();

	//Находим форму струны на втором временном слое 
	OUTu << "{";
	for (int i = 0; i < n - 1; i++) //не забываем что для f берем вторую производную
	{
		solve1[i] = solve0[i] + dt * g(x[i]) + (dt*dt / 2)*((1 / (h*h))*(f(x[i] - dh) - 2 * f(x[i]) + f(x[i] + dh)));
		OUTu << solve1[i] << ",";
	}
	//solve1[n - 1] = solve0[n - 1] + dt * g(x[n - 1]) + (dt*dt / 2)*((1 / (h*h))*(f(x[n - 1] - dh) - 2 * f(x[n - 1]) + f(x[n - 1] + dh)));
	solve1[n - 1] = f((n - 1) * dh);
	OUTu << solve1[n - 1] << "}," << endl;

	//имеем решение на первых двух слоях
	//все готово для поиска решения на остальных слоях
	//хуярим итерационный процесс, не забываем про условия на концах.

	datatype current_t = dt; //так как уже на двух слоях мы уже нашли

	while (current_t<=dt*(m-1))
	{
		current_t += dt;
		cout << "current_t: " << current_t << endl;
        //первое и последнее решение хуярим отдельно, остальные в цикле по пространству.

        OUTu << "{";
		solve2[0] = phi(current_t);
		OUTu << solve2[0] << ","; //записали решение на даноом временном слое в точке х=0
		for (int i = 1; i < n - 1; i++)
		{
			solve2[i] = ((dt*dt) / (dh*dh))*(solve1[i + 1] - 2 * solve1[i] + solve1[i - 1]) + 2 * solve1[i] - solve0[i];
			OUTu << solve2[i] << ",";
		}
		solve2[n - 1] = psi(current_t);
		//solve2[n - 1] = f((n - 1) * dh);
		if (current_t< dt*(m - 1))
		{
			OUTu << solve0[n - 1] << "}," << endl;//если наше время меньше 7 то будет еще запись
		}
		else if (current_t >= dt * (m - 1))
		{
			OUTu << solve0[n - 1] << "}}";//если время уже 7 или больше то записи больше не будет и все пиздец
		}
		

		//оп решение есть.
		//теперь надо сделать шаг вперед по времени на шаблоне
		//то есть сделать solve1=>solve0 и solve2=>solve1
		//не пробуем через буфер(НО СУКА СТОИТ ПОПРОБОВАТЬ)
		//сделаем покомпонентно
		//ATTENTION
		//тут может быть насрано в виде проебов памяти и прочего говна.
		for (int i = 0; i <= n - 1; i++)
		{
			solve0[i] = solve1[i];
		}
		for (int i = 0; i <= n - 1; i++)
		{
			solve1[i] = solve2[i];
		}

	}
	//Бля там в конце запятая но меня не ебет
	OUTu.close();
	DeletingVector(solve0);
	DeletingVector(solve1);
	DeletingVector(solve2);
	DeletingVector(x);
	DeletingVector(t);
}

void Waveanal(funct &f, funct &g, funct &phi, funct &psi)//соответственно начальные условия задаем как входные параметры
{
	// задам х от 0 до 1(почему бы и нет), возьму наверное 50 кусочков, то есть 51 точку 
	// вторую производную функции f для аппроксимации второго слоя буду считать численно
	// также будем считать а=1


	datatype dt = 0.01;
	datatype dh = 0.01;//условие Куранта |a*tao/h|<=1 выполняется
	int n = 101; //всего узлов на отрезке от 0 до 1 с шагом dh  0---5000
	int m = 10; // временных слоёв на указанной мною сетке   0---7000
	datatype raznost;

	for (int step=0;step<4;step++)
	{ 
	tVector solve0 = new datatype[n]; //массив точек точек струны на предыдущем слое.
	tVector solve1 = new datatype[n]; //массив точек точек струны на текущем временном слое.
	tVector solve2 = new datatype[n]; //массив точек точек струны наследующем временном слое.
	for (int i = 0; i < n - 1; i++)
	{
		solve0[i] = f(i * dh);  // начальные условия, форма струны при t=0
	}
	solve0[n - 1] = f((n - 1) * dh);
	//datatype current_t = dt;
	datatype current_t = 2 * dt;
	tVector x = new datatype[n]; //массив пространственной сетки 0---5000
	for (int i = 0; i < n - 1; i++)
	{
		x[i] = i * dh;
	//	OUTx << x[i] << ",";
	}
	x[n - 1] = (n - 1)*dh;

	tVector t = new datatype[m];  // массив временной сетки 0---7000
	for (int i = 0; i < m - 1; i++)
	{
		t[i] = i * dt;
	}
	t[m - 1] = (m - 1)*dt;

	solve1[0] = psi(current_t);
	for (int i = 0; i < n - 1; i++) //не забываем что для f берем вторую производную
	{
		solve1[i] = solve0[i] + dt * g(x[i]) + (dt*dt / 2)*((1 / (h*h))*(f(x[i] - dh) - 2 * f(x[i]) + f(x[i] + dh)));
}
	//solve1[n - 1] = solve0[n - 1] + dt * g(x[n - 1]) + (dt*dt / 2)*((1 / (h*h))*(f(x[n - 1] - dh) - 2 * f(x[n - 1]) + f(x[n - 1] + dh)));
	solve1[n - 1] = f((n - 1) * dh);

	

	while (current_t <= dt * m)
	{
		
		//cout << "current_t: " << current_t << endl;
		//первое и последнее решение хуярим отдельно, остальные в цикле по пространству.

	
		solve2[0] = phi(current_t);
		
		for (int i = 1; i < n - 1; i++)
		{
			solve2[i] = ((dt*dt) / (dh*dh))*(solve1[i + 1] - 2 * solve1[i] + solve1[i - 1]) + 2 * solve1[i] - solve0[i];
			
		}
		solve2[n - 1] = psi(current_t);
		
		for (int i = 0; i <= n - 1; i++)
		{
			solve0[i] = solve1[i];
		}
		for (int i = 0; i <= n - 1; i++)
		{
			solve1[i] = solve2[i];
		}
		
		current_t += dt;//raznost = abs(solve1[(n - 1) / 2] - sin(M_PI*((n - 1) / 2)*dh)*cos(M_PI*current_t));
	}
	raznost = abs(solve1[(n - 1) / 2] - sin(M_PI*((n - 1) / 2)*dh)*cos(M_PI*current_t));
	cout << "raznost"<< raznost << endl;
	dt = dt / 10;
	dh = dh / 10;
	m = m * 10;
	n = ((n-1) * 10)+1;

	DeletingVector(solve0);
	DeletingVector(solve1);
	DeletingVector(solve2);
	DeletingVector(x);
	DeletingVector(t);
	}
}

void NVMultiplication(datatype num, tVector &v, int n)
{
	for(int i=0; i<n; i++)
	{
		v[i]=num*v[i];
	}
}

//void Rozenbrok_Runge(dsys &fsys, tVector &X, int n) //реализация L-устойчивого метода с модификациями
//{
//	char FileName[]="Rozenbrok.txt";
//	ofstream File(FileName);
//	File.setf(ios::fixed);
//	File<<"{";
//	File<<"{";
//	for (int count = 0; count < n - 1; count++)
//	{
//		File<<setprecision(17)<<X[count]<<",";
//	}
//	File<<setprecision(17)<<X[n-1];
//	File<<"}";
//
//	cout<<"Rozenbrok method starts here"<<endl;
//	tVector X0=VectCopy(X,n);
//	datatype dh=0;
//	datatype T=0;
//	cout<<"Enter step: ";
//	cin>>dh;
//	cout<<endl;
//	cout<<"Enter the range: ";
//	cin>>T;
//	cout<<endl;
//	datatype currentT=0;
//
//	datatype gamma=0.57281606248213; //условие метода
//	datatype alpha21=2*gamma;  //условие метода
//	datatype alpha3=(0.2 - 0.25*alpha21)/(0.25 - alpha21/3);// условие метода
//	//datatype b3=0; //условие метода
////	datatype b2=0.5; //предположение  its wrong
////	datatype b4=(1/4 - 1/3 - b2*(8*gamma*gamma*gamma - 4*gamma*gamma))/(alpha3*alpha3*alpha3-alpha3*alpha3); its wrong
//	datatype b4=(8*gamma -3)/(24*gamma*alpha3*alpha3 - 12*alpha3*alpha3*alpha3);
//	datatype b2=(0.25 - b4*alpha3*alpha3*alpha3)/(8*gamma*gamma*gamma);
//	datatype b1=1-b2-b4;
//
//	datatype a=((1./24) - gamma/2 + 1.5*gamma*gamma - gamma*gamma*gamma)/b4;
//	datatype c=((1./6) - gamma + gamma*gamma)/b4;
//
//	//cout<<gamma<<" "<<alpha21<<" "<<alpha3<<" "<<b4<<" "<<b1<<" "<<a<<" "<<c<<endl;
//
//////решим систему для определения beta prime
//	tVector B=new datatype[3];
//	B[0]=0;
//	B[1]=0.5 - gamma;
//	B[2]=-alpha21*alpha21*((1./6) - gamma + gamma*gamma);
//	tMatrix AMatr=nullptr;
//	CreatingArray(AMatr,3);
//	AMatr[0][0]=(a-c)*alpha3*alpha3;
//	AMatr[0][1]=c*alpha21*alpha21;
//	AMatr[0][2]=-a*alpha21*alpha21;
//	AMatr[1][0]=b2;
//	AMatr[1][1]=0;
//	AMatr[1][2]=b4;
//	AMatr[2][0]=b4*alpha3*alpha3 - (1./12) + gamma/3;
//	AMatr[2][1]=-b4*alpha21*alpha21;
//	AMatr[2][2]=0;
//	tVector betaprime=Gauss(AMatr, B, 3);
//	DeletingVector(B);
//	DeletingArray(AMatr, 3);
//////have defined betaprime
//	datatype beta32=a/betaprime[0];
//	datatype beta42=(c-betaprime[1])/(betaprime[0]);
//	datatype alpha32=((1./8) - gamma/3)/(b4*alpha3*betaprime[0]); //equal alpha42
//	alpha3=alpha3 - alpha32; //now its alpha31, equal alpha41
//	datatype gamma21=betaprime[0] - alpha21;
//	datatype gamma32=beta32-alpha32;
//	datatype gamma31=betaprime[1]-beta32-alpha3;
//	datatype gamma42=beta42-alpha32;
//	datatype gamma41=betaprime[3]-beta42-1-alpha3;
//	//gamma43 = 1
//	cout<<"Проверка условий порядка"<<endl;
//	cout<<"(7.15a) "<<b1+b2+b4 - 1<<endl;
//	cout<<"(7.15b) "<<b2*betaprime[0] + b4*betaprime[2] - 0.5 + gamma<<endl;
//	cout<<"(7.15c) "<<b2*alpha21*alpha21 + b4*(alpha3+alpha32)*(alpha3+alpha32) - (1./3)<<endl;
//	//cout<<(1./3)<<endl;
//	cout<<"(7.15d) "<<b4*(beta42*betaprime[0] + betaprime[1]) - (1./6) + gamma - gamma*gamma<<endl;
//	cout<<"(7.15e) "<<b2*alpha21*alpha21*alpha21 + b4*(alpha3+alpha32)*(alpha3+alpha32)*(alpha3+alpha32) -(1./4)<<endl;
//	cout<<"(7.15f) "<<b4*(alpha3+alpha32)*alpha32*betaprime[0]-(1./8) +gamma/3<<endl;
//	cout<<"(7.15g) "<<b4*(beta42*alpha21*alpha21 + (alpha3+alpha32)*(alpha3+alpha32)) - (1./12) + gamma/3<<endl;
//	cout<<"(7.15h) "<<b4*beta32*betaprime[0]-(1./24) + gamma/2 - 1.5*gamma*gamma + gamma*gamma*gamma<<endl;
//	DeletingVector(betaprime);
//	cin.get();
//	//cout<<beta32<<beta42<<alpha3<<alpha32<<gamma21<<gamma31<<gamma41<<gamma42<<endl;
//	//now all coeffs are defined for starting calculating
//
//	tMatrix Jacobi=nullptr;
//	while(currentT<T)
//	{
//			//cout<<"основной"<<endl;
//			CreatingArray(Jacobi,n);
//			for(int i=0; i<n; i++)
//			{ //eps=1.e-3
//				X0[i]=X0[i]-eps;
//				tVector jacminus=fsys(X0,n);
//				X0[i]=X0[i]+2*eps;
//				tVector jacplus=fsys(X0,n);
//				X0[i]=X0[i]-eps;
//				tVector djac=Difference(jacplus, jacminus, n);
//				DeletingVector(jacminus);
//				DeletingVector(jacplus);
//				tVector jac=NumberVectorMultiplication(1/(2*eps),djac, n);
//				DeletingVector(djac);
//				for(int j=0; j<n; j++)
//				{
//					Jacobi[j][i]=jac[j];
//				}
//				DeletingVector(jac);
//			}//calculated Jacobi matrix
//			//PrintingArray(Jacobi, n);
//			tMatrix LU = NumberMatrixMultiplication(-dh*gamma, Jacobi, n);
//			for(int i=0;i<n;i++)
//			{
//				LU[i][i]=LU[i][i]+1;
//			}
//			LUDecomposition_pr(LU, n);
//
//		//calculating k1
//			tVector F1 = fsys(X0,n);
//			NVMultiplication(dh, F1, n);
//			tVector k1=LUSolve(LU, F1, n);
//			DeletingVector(F1);
//			//cout<<"k1"<<endl;
//		//calculating k2
//			tVector X1_2=NumberVectorMultiplication(alpha21, k1, n);
//			tVector X_2=Summ(X0, X1_2, n);
//			DeletingVector(X1_2);
//			tVector F2=fsys(X_2, n);
//			DeletingVector(X_2);
//
//			tVector gk_2 = NumberVectorMultiplication(gamma21, k1, n);;
//			tVector F2pr=MVMultiplication(Jacobi, gk_2, n);
//			DeletingVector(gk_2);
//			tVector RP_2 = Summ(F2, F2pr, n);
//			DeletingVector(F2);
//			DeletingVector(F2pr);
//			NumberVectorMultiplication(dh, RP_2, n);
//			tVector k2=LUSolve(LU, RP_2, n);
//			DeletingVector(RP_2);
//			//cout<<"k2"<<endl;
//		//calculating k3
//			tVector X1_3 = NumberVectorMultiplication(alpha3, k1, n); //alpha3=alpha31
//			tVector X2_3 = NumberVectorMultiplication(alpha32, k2, n);
//			tVector X_3=Summ(X0, X1_3, X2_3, n);
//			DeletingVector(X1_3);
//			DeletingVector(X2_3);
//			tVector F3 = fsys(X_3, n); //equal F4
//			DeletingVector(X_3);
//
//			tVector gk1_3 = NumberVectorMultiplication(gamma31, k1, n);
//			tVector gk2_3 = NumberVectorMultiplication(gamma32, k2, n);
//			tVector gk_3 = Summ(gk1_3, gk2_3, n);
//			DeletingVector(gk1_3);
//			DeletingVector(gk2_3);
//			tVector F3pr = MVMultiplication(Jacobi, gk_3, n);
//			DeletingVector(gk_3);
//
//			tVector RP_3 = Summ(F3, F3pr, n);
//			//DeletingVector(F3); пригодится для k4
//			DeletingVector(F3pr);
//			NumberVectorMultiplication(dh, RP_3, n);
//			tVector k3=LUSolve(LU, RP_3, n);
//			DeletingVector(RP_3);
//
//		//calculating k4
//			tVector gk1_4 = NumberVectorMultiplication(gamma41, k1, n);
//			tVector gk2_4 = NumberVectorMultiplication(gamma42, k2, n);
//			//tVector gammak3 = NumberVectorMultiplication(gamma43, k3, n); gamma43=1
//			tVector gk_4 = Summ(gk1_4, gk2_4, k3, n);
//
//			//cout<<"summ gamma*k"<<endl;
//			//PrintingVector(gk, n);
//
//			DeletingVector(gk1_4);
//			DeletingVector(gk2_4);
//			tVector F4pr = MVMultiplication(Jacobi, gk_4, n);
//			DeletingVector(gk_4);
//
//			tVector RP_4 = Summ(F3, F4pr, n);
//			NVMultiplication(dh, RP_4, n);
//			DeletingVector(F3);
//			DeletingVector(F4pr);
//
//			tVector k4 = LUSolve(LU, RP_4, n);
//			DeletingVector(RP_4);
//			DeletingArray(LU, n);
//			DeletingArray(Jacobi, n);
//			//cout<<"k4"<<endl;
//
//			NVMultiplication(b1, k1, n);
//			NVMultiplication(b2, k2, n);
//			//NVMultiplication(b3, k3, n); equals zero
//			NVMultiplication(b4, k4, n);
//
//			tVector Xn = Summ(X0, k1, k2, k4, n);
//			DeletingVector(k1);
//			DeletingVector(k2);
//			DeletingVector(k3);
//			DeletingVector(k4);
//			//DeletingVector(X0);
//
//			tVector X0pr = VectCopy(X0, n);
//			tVector Xn1 = nullptr;
//			dh=dh/2;
//			for (int i=0; i<2; i++)
//				{
//					//cout<<"вложенный"<<endl;
//					CreatingArray(Jacobi,n);
//					for(int i=0; i<n; i++)
//						{ //eps=1.e-3
//							X0pr[i]=X0pr[i]-eps;
//							tVector jacminus=fsys(X0pr,n);
//							X0pr[i]=X0pr[i]+2*eps;
//							tVector jacplus=fsys(X0pr,n);
//							X0pr[i]=X0pr[i]-eps;
//							tVector djac=Difference(jacplus, jacminus, n);
//							DeletingVector(jacminus);
//							DeletingVector(jacplus);
//							tVector jac=NumberVectorMultiplication(1/(2*eps),djac, n);
//							DeletingVector(djac);
//							for(int j=0; j<n; j++)
//								{
//									Jacobi[j][i]=jac[j];
//								}
//							DeletingVector(jac);
//						}//calculated Jacobi matrix
//					//PrintingArray(Jacobi, n);
//					tMatrix LU = NumberMatrixMultiplication(-dh*gamma, Jacobi, n);
//					for(int i=0;i<n;i++)
//						{
//							LU[i][i]=LU[i][i]+1;
//						}
//					LUDecomposition_pr(LU, n);
//
//				//calculating k1
//					tVector F1 = fsys(X0pr,n);
//					NVMultiplication(dh, F1, n);
//					tVector k1=LUSolve(LU, F1, n);
//					DeletingVector(F1);
//					//cout<<"k1"<<endl;
//				//calculating k2
//					tVector X1_2=NumberVectorMultiplication(alpha21, k1, n);
//					tVector X_2=Summ(X0pr, X1_2, n);
//					DeletingVector(X1_2);
//					tVector F2=fsys(X_2, n);
//					DeletingVector(X_2);
//
//					tVector gk_2 = NumberVectorMultiplication(gamma21, k1, n);;
//					tVector F2pr=MVMultiplication(Jacobi, gk_2, n);
//					DeletingVector(gk_2);
//					tVector RP_2 = Summ(F2, F2pr, n);
//					DeletingVector(F2);
//					DeletingVector(F2pr);
//					//PrintingVector(RPpr, n);
//					NumberVectorMultiplication(dh, RP_2, n);
//					tVector k2=LUSolve(LU, RP_2, n);
//					DeletingVector(RP_2);
//					//cout<<"k2"<<endl;
//				//calculating k3
//					tVector X1_3 = NumberVectorMultiplication(alpha3, k1, n); //alpha3=alpha31
//					tVector X2_3 = NumberVectorMultiplication(alpha32, k2, n);
//					tVector X_3=Summ(X0pr, X1_3, X2_3, n);
//					DeletingVector(X1_3);
//					DeletingVector(X2_3);
//					tVector F3 = fsys(X_3, n); //equal F4
//					DeletingVector(X_3);
//
//					tVector gk1_3 = NumberVectorMultiplication(gamma31, k1, n);
//					tVector gk2_3 = NumberVectorMultiplication(gamma32, k2, n);
//					tVector gk_3 = Summ(gk1_3, gk2_3, n);
//					DeletingVector(gk1_3);
//					DeletingVector(gk2_3);
//					tVector F3pr = MVMultiplication(Jacobi, gk_3, n);
//					DeletingVector(gk_3);
//
//					tVector RP_3 = Summ(F3, F3pr, n);
//					//DeletingVector(F3); пригодится для k4
//					DeletingVector(F3pr);
//					NumberVectorMultiplication(dh, RP_3, n);
//					tVector k3=LUSolve(LU, RP_3, n);
//					DeletingVector(RP_3);
//
//				//calculating k4
//					tVector gk1_4 = NumberVectorMultiplication(gamma41, k1, n);
//					tVector gk2_4 = NumberVectorMultiplication(gamma42, k2, n);
//					//tVector gammak3 = NumberVectorMultiplication(gamma43, k3, n); gamma43=1
//					tVector gk_4 = Summ(gk1_4, gk2_4, k3, n);
//
//					//cout<<"summ gamma*k"<<endl;
//					//PrintingVector(gk, n);
//
//					DeletingVector(gk1_4);
//					DeletingVector(gk2_4);
//					tVector F4pr = MVMultiplication(Jacobi, gk_4, n);
//					DeletingVector(gk_4);
//
//					tVector RP_4 = Summ(F3, F4pr, n);
//					NVMultiplication(dh, RP_4, n);
//					DeletingVector(F3);
//					DeletingVector(F4pr);
//
//					tVector k4 = LUSolve(LU, RP_4, n);
//					DeletingVector(RP_4);
//					DeletingArray(LU, n);
//					DeletingArray(Jacobi, n);
//					//cout<<"k4"<<endl;
//
//					NVMultiplication(b1, k1, n);
//					NVMultiplication(b2, k2, n);
//					//NVMultiplication(b3, k3, n); equals zero
//					NVMultiplication(b4, k4, n);
//
//					Xn1 = Summ(X0pr, k1, k2, k4, n);
//					DeletingVector(k1);
//					DeletingVector(k2);
//					DeletingVector(k3);
//					DeletingVector(k4);
//					DeletingVector(X0pr);
//					X0pr = VectCopy(Xn1, n);
//					DeletingVector(Xn1);
//			}//X0pr - solve. Конец двух шагов по 0.5
//
////			cout<<"X0"<<endl;
////			PrintingVector(X0, n);
////			cout<<"X0pr"<<endl;
////			PrintingVector(X0pr, n);
////			cout<<"Xn"<<endl;
////			PrintingVector(Xn, n);
//
//
//			dh=2*dh;
//			tVector DX=	Difference(Xn, X0pr, n);
//			datatype norm=NormOfVectorC(DX,n);
//
////			cout<<"Vector DX:"<<endl;
////			PrintingVector(DX, n);
//
//			DeletingVector(DX);
//			if (norm<=1.e-5 && norm >=1.e-7)
//			{
//				cout<<"Запись в файл"<<endl;
//				File<<"{";
//				for (int count = 0; count < n - 1; count++)
//				{
//					File<<setprecision(17)<<Xn[count]<<",";
//				}
//				File<<setprecision(17)<<Xn[n-1];
//				File<<"},";
//				DeletingVector(X0);
//				X0=VectCopy(Xn, n);
//				DeletingVector(Xn);
//				DeletingVector(X0pr);
//				currentT = currentT + dh;
//			}
//			if ( norm<1.e-7 )
//			{
//				cout<<"Запись в файл"<<endl;
//				File<<"{";
//				for (int count = 0; count < n - 1; count++)
//				{
//					File<<setprecision(17)<<Xn[count]<<",";
//				}
//				File<<setprecision(17)<<Xn[n-1];
//				File<<"},";
//				DeletingVector(X0);
//				X0=VectCopy(Xn, n);
//				DeletingVector(Xn);
//				DeletingVector(X0pr);
//				currentT = currentT + dh;
//				dh=2*dh;
//			}
//			if (norm>1.e-5)
//			{
//				cout<<"Норма большая"<<endl;
//				dh=dh/2;
//				cout<<dh<<endl;
//				DeletingVector(X0pr);
//				DeletingVector(Xn);
//				//cin.get();
//
//			}
//	}//end while
//
///////////
////	while(true)
////	{
////			cout<<"WHILE"<<endl;
////			CreatingArray(Jacobi,n);
////			for(int i=0; i<n; i++)
////			{ //eps=1.e-3
////				X0[i]=X0[i]-eps;
////				tVector jacminus=fsys(X0,n);
////				X0[i]=X0[i]+2*eps;
////				tVector jacplus=fsys(X0,n);
////				X0[i]=X0[i]-eps;
////				tVector djac=Difference(jacplus, jacminus, n);
////				DeletingVector(jacminus);
////				DeletingVector(jacplus);
////				tVector jac=NumberVectorMultiplication(1/(2*eps),djac, n);
////				DeletingVector(djac);
////				for(int j=0; j<n; j++)
////				{
////					Jacobi[j][i]=jac[j];
////				}
////				DeletingVector(jac);
////			}//calculated Jacobi matrix
////			tMatrix LU = NumberMatrixMultiplication(-h*gamma, Jacobi, n);
////			//PrintingArray(A, n);
////			for(int i=0;i<n;i++)
////			{
////				LU[i][i]=LU[i][i]+1;
////			}
////			LUDecomposition_pr(LU, n);
////
////		//calculating k1
////			tVector F1 = fsys(X0,n);
////			NVMultiplication(h, F1, n);
////			tVector k1=LUSolve(LU, F1, n);
////			DeletingVector(F1);
////			//cout<<"k1"<<endl;
////		//calculating k2
////			tVector X1_2=NumberVectorMultiplication(alpha21, k1, n);
////			tVector X_2=Summ(X0, X1_2, n);
////			DeletingVector(X1_2);
////			tVector F2=fsys(X_2, n);
////			DeletingVector(X_2);
////
////			tVector gk_2 = NumberVectorMultiplication(gamma21, k1, n);
////			tVector F2pr=MVMultiplication(Jacobi, gk_2, n);
////			DeletingVector(gk_2);
////			tVector RP_2 = Summ(F2, F2pr, n);
////			DeletingVector(F2);
////			DeletingVector(F2pr);
////			NumberVectorMultiplication(h, RP_2, n);
////			tVector k2=LUSolve(LU, RP_2, n);
////			DeletingVector(RP_2);
////			//cout<<"k2"<<endl;
////		//calculating k3
////			tVector X1_3 = NumberVectorMultiplication(alpha3, k1, n); //alpha3=alpha31
////
////			tVector X2_3 = NumberVectorMultiplication(alpha32, k2, n);
////			tVector X_3=Summ(X0, X1_3, X2_3, n);
////			DeletingVector(X1_3);
////			DeletingVector(X2_3);
////			tVector F3 = fsys(X_3, n); //equal F4
////			DeletingVector(X_3);
////
////			tVector gk1_3 = NumberVectorMultiplication(gamma31, k1, n);
////			tVector gk2_3 = NumberVectorMultiplication(gamma32, k2, n);
////			tVector gk_3 = Summ(gk1_3, gk2_3, n);
////			DeletingVector(gk1_3);
////			DeletingVector(gk2_3);
////			tVector F3pr = MVMultiplication(Jacobi, gk_3, n);
////			DeletingVector(gk_3);
////
////			tVector RP_3 = Summ(F3, F3pr, n);
////			//DeletingVector(F3); пригодится для k4
////			DeletingVector(F3pr);
////			NumberVectorMultiplication(h, RP_3, n);
////			tVector k3=LUSolve(LU, RP_3, n);
////			DeletingVector(RP_3);
////
////		//calculating k4
////			tVector gk1_4 = NumberVectorMultiplication(gamma41, k1, n);
////			tVector gk2_4 = NumberVectorMultiplication(gamma42, k2, n);
////			//tVector gammak3 = NumberVectorMultiplication(gamma43, k3, n); gamma43=1
////			tVector gk_4 = Summ(gk1_4, gk2_4, k3, n);
////
////			DeletingVector(gk1_4);
////			DeletingVector(gk2_4);
////			tVector F4pr = MVMultiplication(Jacobi, gk_4, n);
////			DeletingVector(gk_4);
////
////			tVector RP_4 = Summ(F3, F4pr, n);
////			NVMultiplication(h, RP_4, n);
////			DeletingVector(F3);
////			DeletingVector(F4pr);
////
////			tVector k4 = LUSolve(LU, RP_4, n);
////			DeletingVector(RP_4);
////			DeletingArray(LU, n);
////			DeletingArray(Jacobi, n);
////			//cout<<"k4"<<endl;
////
////			NVMultiplication(b1, k1, n);
////			NVMultiplication(b2, k2, n);
////			//NVMultiplication(b3, k3, n); equals zero
////			NVMultiplication(b4, k4, n);
////
////			tVector Xn = Summ(X0, k1, k2, k4, n);
////			DeletingVector(k1);
////			DeletingVector(k2);
////			DeletingVector(k3);
////			DeletingVector(k4);
////			//DeletingVector(X0);
////
////			tVector X0pr = VectCopy(X0, n);
////			tVector Xn1 = nullptr;
////			h=h/2;
////			for (int i=0; i<2; i++)
////			{
////				CreatingArray(Jacobi,n);
////				for(int i=0; i<n; i++)
////				{ //eps=1.e-3
////					X0pr[i]=X0pr[i]-eps;
////					tVector jacminus=fsys(X0pr,n);
////					X0pr[i]=X0pr[i]+2*eps;
////					tVector jacplus=fsys(X0pr,n);
////					X0pr[i]=X0pr[i]-eps;
////					tVector djac=Difference(jacplus, jacminus, n);
////					DeletingVector(jacminus);
////					DeletingVector(jacplus);
////					tVector jac=NumberVectorMultiplication(1/(2*eps),djac, n);
////					DeletingVector(djac);
////					for(int j=0; j<n; j++)
////					{
////						Jacobi[j][i]=jac[j];
////					}
////					DeletingVector(jac);
////				}//calculated Jacobi matrix
////				tMatrix LU = NumberMatrixMultiplication(-h*gamma, Jacobi, n);
////				for(int i=0;i<n;i++)
////				{
////					LU[i][i]=LU[i][i]+1;
////				}
////				LUDecomposition_pr(LU, n);
////
////			//calculating k1
////				tVector F1 = fsys(X0pr,n);
////				NVMultiplication(h, F1, n);
////				tVector k1=LUSolve(LU, F1, n);
////				DeletingVector(F1);
////				//cout<<"k1"<<endl;
////			//calculating k2
////				tVector X1_2=NumberVectorMultiplication(alpha21, k1, n);
////				tVector X_2=Summ(X0pr, X1_2, n);
////				DeletingVector(X1_2);
////				tVector F2=fsys(X_2, n);
////				DeletingVector(X_2);
////
////				tVector gk_2 = NumberVectorMultiplication(gamma21, k1, n);;
////				tVector F2pr=MVMultiplication(Jacobi, gk_2, n);
////				DeletingVector(gk_2);
////				tVector RP_2 = Summ(F2, F2pr, n);
////				DeletingVector(F2);
////				DeletingVector(F2pr);
////				//PrintingVector(RPpr, n);
////				NumberVectorMultiplication(h, RP_2, n);
////				tVector k2=LUSolve(LU, RP_2, n);
////				DeletingVector(RP_2);
////				//cout<<"k2"<<endl;
////			//calculating k3
////				tVector X1_3 = NumberVectorMultiplication(alpha3, k1, n); //alpha3=alpha31
////				tVector X2_3 = NumberVectorMultiplication(alpha32, k2, n);
////				tVector X_3=Summ(X0pr, X1_3, X2_3, n);
////				DeletingVector(X1_3);
////				DeletingVector(X2_3);
////				tVector F3 = fsys(X_3, n); //equal F4
////				DeletingVector(X_3);
////
////				tVector gk1_3 = NumberVectorMultiplication(gamma31, k1, n);
////				tVector gk2_3 = NumberVectorMultiplication(gamma32, k2, n);
////				tVector gk_3 = Summ(gk1_3, gk2_3, n);
////				DeletingVector(gk1_3);
////				DeletingVector(gk2_3);
////				tVector F3pr = MVMultiplication(Jacobi, gk_3, n);
////				DeletingVector(gk_3);
////
////				tVector RP_3 = Summ(F3, F3pr, n);
////				//DeletingVector(F3); пригодится для k4
////				DeletingVector(F3pr);
////				NumberVectorMultiplication(h, RP_3, n);
////				tVector k3=LUSolve(LU, RP_3, n);
////				DeletingVector(RP_3);
////
////			//calculating k4
////				tVector gk1_4 = NumberVectorMultiplication(gamma41, k1, n);
////				tVector gk2_4 = NumberVectorMultiplication(gamma42, k2, n);
////				//tVector gammak3 = NumberVectorMultiplication(gamma43, k3, n); gamma43=1
////				tVector gk_4 = Summ(gk1_4, gk2_4, k3, n);
////
////				//cout<<"summ gamma*k"<<endl;
////				//PrintingVector(gk, n);
////
////				DeletingVector(gk1_4);
////				DeletingVector(gk2_4);
////				tVector F4pr = MVMultiplication(Jacobi, gk_4, n);
////				DeletingVector(gk_4);
////
////				tVector RP_4 = Summ(F3, F4pr, n);
////				NVMultiplication(h, RP_4, n);
////				DeletingVector(F3);
////				DeletingVector(F4pr);
////
////				tVector k4 = LUSolve(LU, RP_4, n);
////				DeletingVector(RP_4);
////				DeletingArray(LU, n);
////				DeletingArray(Jacobi, n);
////				//cout<<"k4"<<endl;
////
////				NVMultiplication(b1, k1, n);
////				NVMultiplication(b2, k2, n);
////				//NVMultiplication(b3, k3, n); equals zero
////				NVMultiplication(b4, k4, n);
////
////				Xn1 = Summ(X0pr, k1, k2, k4, n);
////				DeletingVector(k1);
////				DeletingVector(k2);
////				DeletingVector(k3);
////				DeletingVector(k4);
////				DeletingVector(X0pr);
////				X0pr = VectCopy(Xn1, n);
////				DeletingVector(Xn1);
////			}//X0pr - solve
////		h=2*h;
////		tVector DX=	Difference(Xn, X0pr, n);
////		datatype norm=NormOfVectorC(DX,n);
////		DeletingVector(DX);
////		if (norm<=1.e-6 && norm >=1.e-8)
////		{
////			File<<"{";
////			for (int count = 0; count < n - 1; count++)
////			{
////				File<<setprecision(17)<<Xn[count]<<",";
////			}
////			File<<setprecision(17)<<Xn[n-1];
////			File<<"},";
////			DeletingVector(X0);
////			//X0=VectCopy(Xn, n);
////			DeletingVector(Xn);
////			DeletingVector(X0pr);
////			//currentT = currentT + h;
////			break;
////		}
////		if ( norm<1.e-8 )
////		{
////			File<<"{";
////			for (int count = 0; count < n - 1; count++)
////			{
////				File<<setprecision(17)<<Xn[count]<<",";
////			}
////			File<<setprecision(17)<<Xn[n-1];
////			File<<"},";
////			DeletingVector(X0);
////			//X0=VectCopy(Xn, n);
////			DeletingVector(Xn);
////			DeletingVector(X0pr);
////			break;
//////			currentT = currentT + h;
//////			h=3*h;
////		}
////		if (norm>1.e-6)
////		{
////			h=h/2;
////			DeletingVector(X0pr);
////			DeletingVector(Xn);
////
////		}
////	}
//	File<<"}";
//	File.close();
//	//}
//}

void Rozenbrok(dsys &fsys, tVector &X, int n) //реализация L-устойчивого метода с модификациями
{
	char FileName[]="Rozenbrok.txt";
	ofstream File(FileName);
	File.setf(ios::fixed);
	File<<"{";
	File<<"{";
	for (int count = 0; count < n - 1; count++)
	{
		File<<setprecision(17)<<X[count]<<",";
	}
	File<<setprecision(17)<<X[n-1];
	File<<"}";

	cout<<"Rozenbrok method starts here"<<endl;
	tVector X0=VectCopy(X,n);
	datatype dh=0;
	datatype T=0;
	cout<<"Enter step: ";
	cin>>dh;
	cout<<endl;
	cout<<"Enter the range: ";
	cin>>T;
	cout<<endl;
	datatype currentT=0;

	datatype gamma=0.57281606248213; //условие метода
	datatype alpha21=2*gamma;  //условие метода
	datatype alpha3=(0.2 - 0.25*alpha21)/(0.25 - alpha21/3);// условие метода
	//datatype b3=0; //условие метода
//	datatype b2=0.5; //предположение  its wrong
//	datatype b4=(1/4 - 1/3 - b2*(8*gamma*gamma*gamma - 4*gamma*gamma))/(alpha3*alpha3*alpha3-alpha3*alpha3); its wrong
	datatype b4=(8*gamma -3)/(24*gamma*alpha3*alpha3 - 12*alpha3*alpha3*alpha3);
	datatype b2=(0.25 - b4*alpha3*alpha3*alpha3)/(8*gamma*gamma*gamma);
	datatype b1=1-b2-b4;

	datatype a=(1./24 - gamma/2 + 1.5*gamma*gamma - gamma*gamma*gamma)/b4;
	datatype c=(1./6 - gamma + gamma*gamma)/b4;

	//cout<<gamma<<" "<<alpha21<<" "<<alpha3<<" "<<b4<<" "<<b1<<" "<<a<<" "<<c<<endl;

////решим систему для определения beta prime
	tVector B=new datatype[3];
	B[0]=0;
	B[1]=0.5 - gamma;
	B[2]=-alpha21*alpha21*(1./6 - gamma + gamma*gamma);
	tMatrix AMatr=nullptr;
	CreatingArray(AMatr,3);
	AMatr[0][0]=(a-c)*alpha3*alpha3;
	AMatr[0][1]=c*alpha21*alpha21;
	AMatr[0][2]=-a*alpha21*alpha21;
	AMatr[1][0]=b2;
	AMatr[1][1]=0;
	AMatr[1][2]=b4;
	AMatr[2][0]=b4*alpha3*alpha3 - 1./12 + gamma/3;
	AMatr[2][1]=-b4*alpha21*alpha21;
	AMatr[2][2]=0;
	tVector betaprime=Gauss(AMatr, B, 3);
	DeletingVector(B);
	DeletingArray(AMatr, 3);
////have defined betaprime
	datatype beta32=a/betaprime[0];
	datatype beta42=(c-betaprime[1])/(betaprime[0]);
	datatype alpha32=(1./8 - gamma/3)/(b4*alpha3*betaprime[0]); //equal alpha42
	alpha3=alpha3 - alpha32; //now its alpha31, equal alpha41
	datatype gamma21=betaprime[0] - alpha21;
	datatype gamma32=beta32-alpha32;
	datatype gamma31=betaprime[1]-beta32-alpha3;
	datatype gamma42=beta42-alpha32;
	datatype gamma41=betaprime[3]-beta42-1-alpha3;
	//gamma43 = 1
	cout<<"Проверка условий порядка"<<endl;
	cout<<"(7.15a) "<<b1+b2+b4 - 1<<endl;
	cout<<"(7.15b) "<<b2*betaprime[0] + b4*betaprime[2] - 0.5 + gamma<<endl;
	cout<<"(7.15c) "<<b2*alpha21*alpha21 + b4*(alpha3+alpha32)*(alpha3+alpha32) - (1./3)<<endl;
	//cout<<(1./3)<<endl;
	cout<<"(7.15d) "<<b4*(beta42*betaprime[0] + betaprime[1]) - (1./6) + gamma - gamma*gamma<<endl;
	cout<<"(7.15e) "<<b2*alpha21*alpha21*alpha21 + b4*(alpha3+alpha32)*(alpha3+alpha32)*(alpha3+alpha32) -(1./4)<<endl;
	cout<<"(7.15f) "<<b4*(alpha3+alpha32)*alpha32*betaprime[0]-(1./8) +gamma/3<<endl;
	cout<<"(7.15g) "<<b4*(beta42*alpha21*alpha21 + (alpha3+alpha32)*(alpha3+alpha32)) - (1./12) + gamma/3<<endl;
	cout<<"(7.15h) "<<b4*beta32*betaprime[0]-(1./24) + gamma/2 - 1.5*gamma*gamma + gamma*gamma*gamma<<endl;
	DeletingVector(betaprime);
	//cout<<beta32<<beta42<<alpha3<<alpha32<<gamma21<<gamma31<<gamma41<<gamma42<<endl;
	//now all coeffs are defined for starting calculating

	tMatrix Jacobi=nullptr;
	while(currentT<T)
	{
			//cout<<"основной"<<endl;
			CreatingArray(Jacobi,n);
			for(int i=0; i<n; i++)
			{ //eps=1.e-3
				X0[i]=X0[i]-eps;
				tVector jacminus=fsys(X0,n);
				X0[i]=X0[i]+2*eps;
				tVector jacplus=fsys(X0,n);
				X0[i]=X0[i]-eps;
				tVector djac=Difference(jacplus, jacminus, n);
				DeletingVector(jacminus);
				DeletingVector(jacplus);
				tVector jac=NumberVectorMultiplication(1/(2*eps),djac, n);
				DeletingVector(djac);
				for(int j=0; j<n; j++)
				{
					Jacobi[j][i]=jac[j];
				}
				DeletingVector(jac);
			}//calculated Jacobi matrix
			tMatrix LU = NumberMatrixMultiplication(-dh*gamma, Jacobi, n);
			for(int i=0;i<n;i++)
			{
				LU[i][i]=LU[i][i]+1;
			}
			LUDecomposition_pr(LU, n);

		//calculating k1
			tVector F1 = fsys(X0,n);
			NVMultiplication(dh, F1, n);
			tVector k1=LUSolve(LU, F1, n);
			DeletingVector(F1);
			//cout<<"k1"<<endl;
		//calculating k2
			tVector X1_2=NumberVectorMultiplication(alpha21, k1, n);
			tVector X_2=Summ(X0, X1_2, n);
			DeletingVector(X1_2);
			tVector F2=fsys(X_2, n);
			DeletingVector(X_2);

			tVector gk_2 = NumberVectorMultiplication(gamma21, k1, n);;
			tVector F2pr=MVMultiplication(Jacobi, gk_2, n);
			DeletingVector(gk_2);
			tVector RP_2 = Summ(F2, F2pr, n);
			DeletingVector(F2);
			DeletingVector(F2pr);
			//PrintingVector(RPpr, n);
			NumberVectorMultiplication(dh, RP_2, n);
			tVector k2=LUSolve(LU, RP_2, n);
			DeletingVector(RP_2);
			//cout<<"k2"<<endl;
		//calculating k3
			tVector X1_3 = NumberVectorMultiplication(alpha3, k1, n); //alpha3=alpha31
			tVector X2_3 = NumberVectorMultiplication(alpha32, k2, n);
			tVector X_3=Summ(X0, X1_3, X2_3, n);
			DeletingVector(X1_3);
			DeletingVector(X2_3);
			tVector F3 = fsys(X_3, n); //equal F4
			DeletingVector(X_3);

			tVector gk1_3 = NumberVectorMultiplication(gamma31, k1, n);
			tVector gk2_3 = NumberVectorMultiplication(gamma32, k2, n);
			tVector gk_3 = Summ(gk1_3, gk2_3, n);
			DeletingVector(gk1_3);
			DeletingVector(gk2_3);
			tVector F3pr = MVMultiplication(Jacobi, gk_3, n);
			DeletingVector(gk_3);

			tVector RP_3 = Summ(F3, F3pr, n);
			//DeletingVector(F3); пригодится для k4
			DeletingVector(F3pr);
			NumberVectorMultiplication(dh, RP_3, n);
			tVector k3=LUSolve(LU, RP_3, n);
			DeletingVector(RP_3);

		//calculating k4
			tVector gk1_4 = NumberVectorMultiplication(gamma41, k1, n);
			tVector gk2_4 = NumberVectorMultiplication(gamma42, k2, n);
			//tVector gammak3 = NumberVectorMultiplication(gamma43, k3, n); gamma43=1
			tVector gk_4 = Summ(gk1_4, gk2_4, k3, n);

			DeletingVector(gk1_4);
			DeletingVector(gk2_4);
			tVector F4pr = MVMultiplication(Jacobi, gk_4, n);
			DeletingVector(gk_4);

			tVector RP_4 = Summ(F3, F4pr, n);
			NVMultiplication(dh, RP_4, n);
			DeletingVector(F3);
			DeletingVector(F4pr);

			tVector k4 = LUSolve(LU, RP_4, n);
			DeletingVector(RP_4);
			DeletingArray(LU, n);
			DeletingArray(Jacobi, n);
			//cout<<"k4"<<endl;

			NVMultiplication(b1, k1, n);
			NVMultiplication(b2, k2, n);
			//NVMultiplication(b3, k3, n); equals zero
			NVMultiplication(b4, k4, n);

			tVector Xn = Summ(X0, k1, k2, k4, n);
			DeletingVector(k1);
			DeletingVector(k2);
			DeletingVector(k3);
			DeletingVector(k4);
			File<<"{";
			for (int i=0; i<n-1; i++)
			{
				File<<Xn[i]<<",";
			}
			File<< setprecision(17) <<Xn[n-1]<<"},";
			DeletingVector(X0);
			X0 = VectCopy(Xn, n);
			DeletingVector(Xn);
			currentT=currentT+dh;
	}//end while(currentT<T)


	CreatingArray(Jacobi,n);
	for(int i=0; i<n; i++)
	{ //eps=1.e-3
		X0[i]=X0[i]-eps;
		tVector jacminus=fsys(X0,n);
		X0[i]=X0[i]+2*eps;
		tVector jacplus=fsys(X0,n);
		X0[i]=X0[i]-eps;
		tVector djac=Difference(jacplus, jacminus, n);
		DeletingVector(jacminus);
		DeletingVector(jacplus);
		tVector jac=NumberVectorMultiplication(1/(2*eps),djac, n);
		DeletingVector(djac);
		for(int j=0; j<n; j++)
		{
			Jacobi[j][i]=jac[j];
		}
		DeletingVector(jac);
	}//calculated Jacobi matrix
	tMatrix LU = NumberMatrixMultiplication(-dh*gamma, Jacobi, n);
	for(int i=0;i<n;i++)
	{
		LU[i][i]=LU[i][i]+1;
	}
	LUDecomposition_pr(LU, n);

//calculating k1
	tVector F1 = fsys(X0,n);
	NVMultiplication(dh, F1, n);
	tVector k1=LUSolve(LU, F1, n);
	DeletingVector(F1);
	//cout<<"k1"<<endl;
//calculating k2
	tVector X1_2=NumberVectorMultiplication(alpha21, k1, n);
	tVector X_2=Summ(X0, X1_2, n);
	DeletingVector(X1_2);
	tVector F2=fsys(X_2, n);
	DeletingVector(X_2);

	tVector gk_2 = NumberVectorMultiplication(gamma21, k1, n);;
	tVector F2pr=MVMultiplication(Jacobi, gk_2, n);
	DeletingVector(gk_2);
	tVector RP_2 = Summ(F2, F2pr, n);
	DeletingVector(F2);
	DeletingVector(F2pr);
	//PrintingVector(RPpr, n);
	NumberVectorMultiplication(dh, RP_2, n);
	tVector k2=LUSolve(LU, RP_2, n);
	DeletingVector(RP_2);
	//cout<<"k2"<<endl;
//calculating k3
	tVector X1_3 = NumberVectorMultiplication(alpha3, k1, n); //alpha3=alpha31
	tVector X2_3 = NumberVectorMultiplication(alpha32, k2, n);
	tVector X_3=Summ(X0, X1_3, X2_3, n);
	DeletingVector(X1_3);
	DeletingVector(X2_3);
	tVector F3 = fsys(X_3, n); //equal F4
	DeletingVector(X_3);

	tVector gk1_3 = NumberVectorMultiplication(gamma31, k1, n);
	tVector gk2_3 = NumberVectorMultiplication(gamma32, k2, n);
	tVector gk_3 = Summ(gk1_3, gk2_3, n);
	DeletingVector(gk1_3);
	DeletingVector(gk2_3);
	tVector F3pr = MVMultiplication(Jacobi, gk_3, n);
	DeletingVector(gk_3);

	tVector RP_3 = Summ(F3, F3pr, n);
	//DeletingVector(F3); пригодится для k4
	DeletingVector(F3pr);
	NumberVectorMultiplication(dh, RP_3, n);
	tVector k3=LUSolve(LU, RP_3, n);
	DeletingVector(RP_3);

//calculating k4
	tVector gk1_4 = NumberVectorMultiplication(gamma41, k1, n);
	tVector gk2_4 = NumberVectorMultiplication(gamma42, k2, n);
	//tVector gammak3 = NumberVectorMultiplication(gamma43, k3, n); gamma43=1
	tVector gk_4 = Summ(gk1_4, gk2_4, k3, n);

	//cout<<"summ gamma*k"<<endl;
	//PrintingVector(gk, n);

	DeletingVector(gk1_4);
	DeletingVector(gk2_4);
	tVector F4pr = MVMultiplication(Jacobi, gk_4, n);
	DeletingVector(gk_4);

	tVector RP_4 = Summ(F3, F4pr, n);
	NVMultiplication(dh, RP_4, n);
	DeletingVector(F3);
	DeletingVector(F4pr);

	tVector k4 = LUSolve(LU, RP_4, n);
	DeletingVector(RP_4);
	DeletingArray(LU, n);
	DeletingArray(Jacobi, n);
	//cout<<"k4"<<endl;

	NVMultiplication(b1, k1, n);
	NVMultiplication(b2, k2, n);
	//NVMultiplication(b3, k3, n); equals zero
	NVMultiplication(b4, k4, n);

	tVector Xn = Summ(X0, k1, k2, k4, n);
	DeletingVector(k1);
	DeletingVector(k2);
	DeletingVector(k3);
	DeletingVector(k4);
	File<<"{";
	for (int i=0; i<n-1; i++)
	{
		File<<Xn[i]<<",";
	}
	File<< setprecision(17) <<Xn[n-1]<<"}";

	File<<"}";
	File.close();
	DeletingVector(X0);
	DeletingVector(Xn);
}

void Rozenbrok_Runge(dsys& fsys, tVector &X, int n)
{
	char FileName[]="Rozenbrok.txt";
	ofstream File(FileName);
	char FileName2[]="RozenbrokStep.txt";
	File.setf(ios::fixed);
	ofstream FileStep(FileName2);
	FileStep.setf(ios::fixed);
	FileStep<<"{";
	File<<"{";
	File<<"{";
	for (int count = 0; count < n - 1; count++)
	{
		File<<setprecision(17)<<X[count]<<",";
	}
	File<<setprecision(17)<<X[n-1];
	File<<"},";

	cout<<"Rozenbrok method starts here"<<endl;
	tVector X0=VectCopy(X,n);
	datatype dh=0;
	datatype T=0;
	cout<<"Enter step: ";
	cin>>dh;
	cout<<endl;
	cout<<"Enter the range: ";
	cin>>T;
	cout<<endl;
	datatype currentT=0;

	datatype gamma=0.57281606248213; //условие метода
	datatype alpha21=2*gamma;  //условие метода
	datatype alpha3=(0.2 - 0.25*alpha21)/(0.25 - alpha21/3);// условие метода
	//datatype b3=0; //условие метода
//	datatype b2=0.5; //предположение  its wrong
//	datatype b4=(1/4 - 1/3 - b2*(8*gamma*gamma*gamma - 4*gamma*gamma))/(alpha3*alpha3*alpha3-alpha3*alpha3); its wrong
	datatype b4=(8*gamma -3)/(24*gamma*alpha3*alpha3 - 12*alpha3*alpha3*alpha3);
	datatype b2=(0.25 - b4*alpha3*alpha3*alpha3)/(8*gamma*gamma*gamma);
	datatype b1=1-b2-b4;

	datatype a=(1./24 - gamma/2 + 1.5*gamma*gamma - gamma*gamma*gamma)/b4;
	datatype c=(1./6 - gamma + gamma*gamma)/b4;

	//cout<<gamma<<" "<<alpha21<<" "<<alpha3<<" "<<b4<<" "<<b1<<" "<<a<<" "<<c<<endl;

////решим систему для определения beta prime
	tVector B=new datatype[3];
	B[0]=0;
	B[1]=0.5 - gamma;
	B[2]=-alpha21*alpha21*(1./6 - gamma + gamma*gamma);
	tMatrix AMatr=nullptr;
	CreatingArray(AMatr,3);
	AMatr[0][0]=(a-c)*alpha3*alpha3;
	AMatr[0][1]=c*alpha21*alpha21;
	AMatr[0][2]=-a*alpha21*alpha21;
	AMatr[1][0]=b2;
	AMatr[1][1]=0;
	AMatr[1][2]=b4;
	AMatr[2][0]=b4*alpha3*alpha3 - 1./12 + gamma/3;
	AMatr[2][1]=-b4*alpha21*alpha21;
	AMatr[2][2]=0;
	tVector betaprime=Gauss(AMatr, B, 3);
	DeletingVector(B);
	DeletingArray(AMatr, 3);
////have defined betaprime
	datatype beta32=a/betaprime[0];
	datatype beta42=(c-betaprime[1])/(betaprime[0]);
	datatype alpha32=(1./8 - gamma/3)/(b4*alpha3*betaprime[0]); //equal alpha42
	alpha3=alpha3 - alpha32; //now its alpha31, equal alpha41
	datatype gamma21=betaprime[0] - alpha21;
	datatype gamma32=beta32-alpha32;
	datatype gamma31=betaprime[1]-beta32-alpha3;
	datatype gamma42=beta42-alpha32;
	datatype gamma41=betaprime[3]-beta42-1-alpha3;
	//gamma43 = 1
	cout<<"Проверка условий порядка"<<endl;
	cout<<"(7.15a) "<<b1+b2+b4 - 1<<endl;
	cout<<"(7.15b) "<<b2*betaprime[0] + b4*betaprime[2] - 0.5 + gamma<<endl;
	cout<<"(7.15c) "<<b2*alpha21*alpha21 + b4*(alpha3+alpha32)*(alpha3+alpha32) - (1./3)<<endl;
	//cout<<(1./3)<<endl;
	cout<<"(7.15d) "<<b4*(beta42*betaprime[0] + betaprime[1]) - (1./6) + gamma - gamma*gamma<<endl;
	cout<<"(7.15e) "<<b2*alpha21*alpha21*alpha21 + b4*(alpha3+alpha32)*(alpha3+alpha32)*(alpha3+alpha32) -(1./4)<<endl;
	cout<<"(7.15f) "<<b4*(alpha3+alpha32)*alpha32*betaprime[0]-(1./8) +gamma/3<<endl;
	cout<<"(7.15g) "<<b4*(beta42*alpha21*alpha21 + (alpha3+alpha32)*(alpha3+alpha32)) - (1./12) + gamma/3<<endl;
	cout<<"(7.15h) "<<b4*beta32*betaprime[0]-(1./24) + gamma/2 - 1.5*gamma*gamma + gamma*gamma*gamma<<endl;
	DeletingVector(betaprime);
	//cout<<beta32<<beta42<<alpha3<<alpha32<<gamma21<<gamma31<<gamma41<<gamma42<<endl;
	//now all coeffs are defined for starting calculating

	tMatrix Jacobi_initial = nullptr;
	CreatingArray(Jacobi_initial, n);
	for (int i=0; i<n; i++)
	{
		X0[i]=X0[i]+eps5;
		tVector jacplus=fsys(X0, n);
		X0[i]=X0[i]-2*eps5;
		tVector jacminus=fsys(X0, n);
		X0[i]=X0[i]+eps5;
		tVector jac = Difference(jacplus, jacminus, n);
		NVMultiplication(1/(2*eps5), jac, n);
		for (int j=0; j<n ; j++)
		{
			Jacobi_initial[j][i]=jac[j];
		}
		DeletingVector(jacplus);
		DeletingVector(jacminus);
		DeletingVector(jac);
	}//end Jacobi calculating
	tMatrix LU_initial = NumberMatrixMultiplication(-dh*gamma, Jacobi_initial, n);
	for (int j=0; j<n ; j++)
	{
		LU_initial[j][j]++;
	}
	LUDecomposition_pr(LU_initial, n);

	tMatrix LU_middle = nullptr;
	tMatrix Jacobi_middle = nullptr;
	CreatingArray(Jacobi_middle, n);
	while (currentT<T)
	{
//		for (int i=0; i<n; i++)
//		{
//			X0[i]=X0[i]+eps5;
//			tVector jacplus=fsys(X0, n);
//			X0[i]=X0[i]-2*eps5;
//			tVector jacminus=fsys(X0, n);
//			X0[i]=X0[i]+eps5;
//			tVector jac = Difference(jacplus, jacminus, n);
//			NVMultiplication(1/(2*eps5), jac, n);
//			for (int j=0; j<n ; j++)
//			{
//				Jacobi_initial[j][i]=jac[j];
//			}
//			DeletingVector(jacplus);
//			DeletingVector(jacminus);
//			DeletingVector(jac);
//		}//end Jacobi calculating
//		tMatrix LU_initial = NumberMatrixMultiplication(-dh*gamma, Jacobi_initial, n);
//		for (int j=0; j<n ; j++)
//		{
//			LU_initial[j][j]++;
//		}
//		LUDecomposition_pr(LU_initial, n);

	//k1
		tVector F1 = fsys(X0, n);
		NVMultiplication(dh, F1, n);
		tVector k1=LUSolve(LU_initial, F1, n);
		DeletingVector(F1);
	//k2
		tVector X2_1 = NumberVectorMultiplication(alpha21, k1, n);
		tVector X2 = Summ(X0, X2_1, n);
		tVector F2 = fsys(X2, n);
		tVector gk2 = NumberVectorMultiplication(gamma21, k1, n);
		tVector F2pr = MVMultiplication(Jacobi_initial, gk2, n);
		tVector RP_2=Summ(F2, F2pr, n);
		NVMultiplication(dh, RP_2, n);
		tVector k2 = LUSolve(LU_initial, RP_2, n);
		DeletingVector(X2_1);
		DeletingVector(X2);
		DeletingVector(F2);
		DeletingVector(gk2);
		DeletingVector(F2pr);
		DeletingVector(RP_2);
	//k3
		tVector X3_1 = NumberVectorMultiplication(alpha3, k1, n);
		tVector X3_2 = NumberVectorMultiplication(alpha32, k2, n);
		tVector X3 = Summ(X0, X3_1, X3_2, n);
		tVector F3 = fsys(X3, n);
		tVector gk3_1 = NumberVectorMultiplication(gamma31, k1, n);
		tVector gk3_2 = NumberVectorMultiplication(gamma32, k2, n);
		tVector gk3 = Summ(gk3_1, gk3_2, n);
		tVector F3pr = MVMultiplication(Jacobi_initial, gk3, n);
		tVector RP_3=Summ(F3, F3pr, n);
		NVMultiplication(dh, RP_3, n);
		tVector k3 = LUSolve(LU_initial, RP_3, n);
		DeletingVector(X3_1);
		DeletingVector(X3_2);
		DeletingVector(X3);
		DeletingVector(gk3_1);
		DeletingVector(gk3_2);
		DeletingVector(gk3);
		DeletingVector(F3pr);
		DeletingVector(RP_3);
	//k4
		tVector gk4_1 = NumberVectorMultiplication(gamma41, k1, n);
		tVector gk4_2 = NumberVectorMultiplication(gamma42, k2, n);
		tVector gk4 = Summ(gk4_1, gk4_2, k3, n); //gamma43 = 1
		tVector F4pr = MVMultiplication(Jacobi_initial, gk4, n);
		tVector RP_4 = Summ(F3, F4pr, n);
		NVMultiplication(dh, RP_4, n);
		tVector k4 = LUSolve(LU_initial, RP_4, n);
		DeletingVector(gk4_1);
		DeletingVector(gk4_2);
		DeletingVector(gk4);
		DeletingVector(F4pr);
		DeletingVector(RP_4);
		DeletingVector(F3);
	//X1
		NVMultiplication(dh, k1, n);
		NVMultiplication(dh, k2, n);
		NVMultiplication(dh, k4, n);
		tVector X1 = Summ(X0, k1, k2, k4, n);
		DeletingVector(k1);
		DeletingVector(k2);
		DeletingVector(k3);
		DeletingVector(k4);
/////////////////////////////////////   X_middle   /////////////////////////////////////////////////
		dh=dh/2;
		tVector X0pr = VectCopy(X0, n);
	//k1
		F1 = fsys(X0pr, n);
		NVMultiplication(dh, F1, n);
		k1=LUSolve(LU_initial, F1, n);
		DeletingVector(F1);
	//k2
		X2_1 = NumberVectorMultiplication(alpha21, k1, n);
		X2 = Summ(X0pr, X2_1, n);
		F2 = fsys(X2, n);
		gk2 = NumberVectorMultiplication(gamma21, k1, n);
		F2pr = MVMultiplication(Jacobi_initial, gk2, n);
		RP_2=Summ(F2, F2pr, n);
		NVMultiplication(dh, RP_2, n);
		k2 = LUSolve(LU_initial, RP_2, n);
		DeletingVector(X2_1);
		DeletingVector(X2);
		DeletingVector(F2);
		DeletingVector(gk2);
		DeletingVector(F2pr);
		DeletingVector(RP_2);
	//k3
		X3_1 = NumberVectorMultiplication(alpha3, k1, n);
		X3_2 = NumberVectorMultiplication(alpha32, k2, n);
		X3 = Summ(X0pr, X3_1, X3_2, n);
		F3 = fsys(X3, n);
		gk3_1 = NumberVectorMultiplication(gamma31, k1, n);
		gk3_2 = NumberVectorMultiplication(gamma32, k2, n);
		gk3 = Summ(gk3_1, gk3_2, n);
		F3pr = MVMultiplication(Jacobi_initial, gk3, n);
		RP_3=Summ(F3, F3pr, n);
		NVMultiplication(dh, RP_3, n);
		k3 = LUSolve(LU_initial, RP_3, n);
		DeletingVector(X3_1);
		DeletingVector(X3_2);
		DeletingVector(X3);
		DeletingVector(gk3_1);
		DeletingVector(gk3_2);
		DeletingVector(gk3);
		DeletingVector(F3pr);
		DeletingVector(RP_3);
	//k4
		gk4_1 = NumberVectorMultiplication(gamma41, k1, n);
		gk4_2 = NumberVectorMultiplication(gamma42, k2, n);
		gk4 = Summ(gk4_1, gk4_2, k3, n); //gamma43 = 1
		F4pr = MVMultiplication(Jacobi_initial, gk4, n);
		RP_4 = Summ(F3, F4pr, n);
		NVMultiplication(dh, RP_4, n);
		k4 = LUSolve(LU_initial, RP_4, n);
		DeletingVector(gk4_1);
		DeletingVector(gk4_2);
		DeletingVector(gk4);
		DeletingVector(F4pr);
		DeletingVector(RP_4);
		DeletingVector(F3);
	//X_middle
		NVMultiplication(dh, k1, n);
		NVMultiplication(dh, k2, n);
		NVMultiplication(dh, k4, n);
		tVector X1pr = Summ(X0pr, k1, k2, k4, n);
		DeletingVector(k1);
		DeletingVector(k2);
		DeletingVector(k3);
		DeletingVector(k4);
		DeletingVector(X0pr);
		X0pr = VectCopy(X1pr, n);
		DeletingVector(X1pr);
/////////////////////////////////////   X1pr   /////////////////////////////////////////////////
		for (int i=0; i<n; i++)
		{
			X0pr[i]=X0pr[i]+eps5;
			tVector jacplus=fsys(X0pr, n);
			X0pr[i]=X0pr[i]-2*eps5;
			tVector jacminus=fsys(X0pr, n);
			X0pr[i]=X0pr[i]+eps5;
			tVector jac = Difference(jacplus, jacminus, n);
			NVMultiplication(1/(2*eps5), jac, n);
			for (int j=0; j<n ; j++)
			{
				Jacobi_middle[j][i]=jac[j];
			}
			DeletingVector(jacplus);
			DeletingVector(jacminus);
			DeletingVector(jac);
		}//end Jacobi_middle calculating
		LU_middle = NumberMatrixMultiplication(-dh*gamma, Jacobi_middle, n);
		for (int j=0; j<n ; j++)
		{
			LU_middle[j][j]++;
		}
		LUDecomposition_pr(LU_middle, n);

		//k1
			F1 = fsys(X0pr, n);
			NVMultiplication(dh, F1, n);
			k1=LUSolve(LU_middle, F1, n);
			DeletingVector(F1);
		//k2
			X2_1 = NumberVectorMultiplication(alpha21, k1, n);
			X2 = Summ(X0pr, X2_1, n);
			F2 = fsys(X2, n);
			gk2 = NumberVectorMultiplication(gamma21, k1, n);
			F2pr = MVMultiplication(Jacobi_middle, gk2, n);
			RP_2=Summ(F2, F2pr, n);
			NVMultiplication(dh, RP_2, n);
			k2 = LUSolve(LU_middle, RP_2, n);
			DeletingVector(X2_1);
			DeletingVector(X2);
			DeletingVector(F2);
			DeletingVector(gk2);
			DeletingVector(F2pr);
			DeletingVector(RP_2);
		//k3
			X3_1 = NumberVectorMultiplication(alpha3, k1, n);
			X3_2 = NumberVectorMultiplication(alpha32, k2, n);
			X3 = Summ(X0pr, X3_1, X3_2, n);
			F3 = fsys(X3, n);
			gk3_1 = NumberVectorMultiplication(gamma31, k1, n);
			gk3_2 = NumberVectorMultiplication(gamma32, k2, n);
			gk3 = Summ(gk3_1, gk3_2, n);
			F3pr = MVMultiplication(Jacobi_middle, gk3, n);
			RP_3=Summ(F3, F3pr, n);
			NVMultiplication(dh, RP_3, n);
			k3 = LUSolve(LU_middle, RP_3, n);
			DeletingVector(X3_1);
			DeletingVector(X3_2);
			DeletingVector(X3);
			DeletingVector(gk3_1);
			DeletingVector(gk3_2);
			DeletingVector(gk3);
			DeletingVector(F3pr);
			DeletingVector(RP_3);
		//k4
			gk4_1 = NumberVectorMultiplication(gamma41, k1, n);
			gk4_2 = NumberVectorMultiplication(gamma42, k2, n);
			gk4 = Summ(gk4_1, gk4_2, k3, n); //gamma43 = 1
			F4pr = MVMultiplication(Jacobi_middle, gk4, n);
			RP_4 = Summ(F3, F4pr, n);
			NVMultiplication(dh, RP_4, n);
			k4 = LUSolve(LU_middle, RP_4, n);
			DeletingVector(gk4_1);
			DeletingVector(gk4_2);
			DeletingVector(gk4);
			DeletingVector(F4pr);
			DeletingVector(RP_4);
			DeletingVector(F3);
		//X1pr
			NVMultiplication(dh, k1, n);
			NVMultiplication(dh, k2, n);
			NVMultiplication(dh, k4, n);
			X1pr = Summ(X0pr, k1, k2, k4, n);
			DeletingVector(k1);
			DeletingVector(k2);
			DeletingVector(k3);
			DeletingVector(k4);
			DeletingVector(X0pr);
			DeletingArray(LU_middle, n);

			dh=dh*2;
			tVector DX = Difference(X1, X1pr, n);
			datatype norm = NormOfVectorC(DX, n);
			norm = norm/15;
			if (norm <= 1.e-6 && norm >=1.e-8)
			{
				File<<"{";
				FileStep<<setprecision(20)<<dh<<",";
				for (int i=0; i<n-1; i++)
				{
					File<<setprecision(17)<<X1[i]<<",";
				}
				File<< setprecision(17) <<X1[n-1]<<"},";
				DeletingVector(X0);
				X0 = VectCopy(X1, n);
				DeletingVector(X1);
				DeletingVector(X1pr);

				for (int i=0; i<n; i++)
				{
					X0[i]=X0[i]+eps5;
					tVector jacplus=fsys(X0, n);
					X0[i]=X0[i]-2*eps5;
					tVector jacminus=fsys(X0, n);
					X0[i]=X0[i]+eps5;
					tVector jac = Difference(jacplus, jacminus, n);
					NVMultiplication(1/(2*eps5), jac, n);
					for (int j=0; j<n ; j++)
					{
						Jacobi_initial[j][i]=jac[j];
					}
					DeletingVector(jacplus);
					DeletingVector(jacminus);
					DeletingVector(jac);
				}//end Jacobi calculating
				DeletingArray(LU_initial, n);
				LU_initial = NumberMatrixMultiplication(-dh*gamma, Jacobi_initial, n);
				for (int j=0; j<n ; j++)
				{
					LU_initial[j][j]++;
				}
				LUDecomposition_pr(LU_initial, n);

				DeletingArray(Jacobi_middle, n);
				CreatingArray(Jacobi_middle, n);

				currentT=currentT+dh;
			}
			if (norm < 1.e-8)
			{
				File<<"{";
				FileStep<<setprecision(20)<<dh<<",";
				for (int i=0; i<n-1; i++)
				{
					File<<setprecision(17)<<X1[i]<<",";
				}
				File<< setprecision(17) <<X1[n-1]<<"},";
				DeletingVector(X0);
				X0 = VectCopy(X1, n);
				DeletingVector(X1);
				DeletingVector(X1pr);

				for (int i=0; i<n; i++)
				{
					X0[i]=X0[i]+eps5;
					tVector jacplus=fsys(X0, n);
					X0[i]=X0[i]-2*eps5;
					tVector jacminus=fsys(X0, n);
					X0[i]=X0[i]+eps5;
					tVector jac = Difference(jacplus, jacminus, n);
					NVMultiplication(1/(2*eps5), jac, n);
					for (int j=0; j<n ; j++)
					{
						Jacobi_initial[j][i]=jac[j];
					}
					DeletingVector(jacplus);
					DeletingVector(jacminus);
					DeletingVector(jac);
				}//end Jacobi calculating
				DeletingArray(LU_initial, n);
				LU_initial = NumberMatrixMultiplication(-dh*gamma, Jacobi_initial, n);
				for (int j=0; j<n ; j++)
				{
					LU_initial[j][j]++;
				}
				LUDecomposition_pr(LU_initial, n);

				DeletingArray(Jacobi_middle, n);
				CreatingArray(Jacobi_middle, n);

				currentT=currentT+dh;
				dh=2*dh;
			}
			if (norm > 1.e-6)
			{
				dh=dh/2;
				DeletingVector(X1);
				DeletingVector(X1pr);
				DeletingArray(Jacobi_middle, n);
				CreatingArray(Jacobi_middle, n);
			}

	}//end while (currentT<T)
	File<<"}";
	FileStep<<"}";
	FileStep.close();
	File.close();
}

tVector Runge1(tVector &X0, int n)
{
	tVector sol = new datatype[n];
	sol[0] = -X0[0];
	sol[1] = -1000000*X0[1];
	return sol;
}

datatype f_fe(datatype x, datatype t)
{
	return 0;
}

datatype g_fe(datatype x, datatype t)
{
	return 0;
}

datatype phi_fe(datatype x, datatype t)
{
	return 0;
}

datatype psi_fe(datatype x, datatype t)
{
	return 0;
}

datatype T_fe(datatype x, datatype t, datatype betta)
{
	return betta * (exp(-t-x*betta)*(exp(t)-1-t)) / (-1 + exp(betta * 10));
	//	return 0;
}

//void Finite_Element(funct2 &f, funct2 &g, funct2 &phi, funct2 &psi, funct2 T)
//{
//	char FileName[] = "x_FE.txt";
//	char FileName2[] = "t_FE.txt";
//	char FileName3[] = "u_FE.txt";
//	ofstream OUTx(FileName);
//	ofstream OUTt(FileName2);
//	ofstream OUTu(FileName3);
//
//	//хуярим в файлы в нормальном виде без ебаных е кстати
//	OUTx.setf(ios::fixed);
//	OUTt.setf(ios::fixed);
//	OUTu.setf(ios::fixed);
//
//	//datatype L = 10.;
//
//	datatype dt = 0.01;
//	datatype dh = 0.2;//условие Куранта |a*tao/h|<=1 выполняется
//	int n = 51; //всего узлов на отрезке от 0 до 1 с шагом dh  0---5000
//	int m = 701; // временных слоёв на указанной мною сетке   0---7000
//
//	//хуярим на 7 секунд ебать
//
//	OUTx << "{";
//	tVector x = new datatype[n]; //массив пространственной сетки 0---5000
//	for (int i = 0; i < n - 1; i++)
//	{
//		x[i] = i * dh;
//		OUTx << x[i] << ",";
//	}
//	x[n - 1] = (n - 1)*dh;
//	OUTx << x[n - 1] << "}"; //запись простанственной сетки в файл
//	OUTx.close();
//
//	OUTt << "{";
//	tVector t = new datatype[m];  // массив временной сетки 0---7000
//	for (int i = 0; i < m - 1; i++)
//	{
//		t[i] = i * dt;
//		OUTt << t[i] << ",";
//	}
//	t[m - 1] = (m - 1)*dt;
//	OUTt << t[m - 1] << "}"; //запись временнОй сетки в файл
//	OUTt.close();
//
//
//	OUTu << "{" << endl;
//	OUTu << "{";
//	tVector solve0 = new datatype[n]; //массив точек точек струны на предыдущем слое.
//	tVector solve1 = new datatype[n]; //массив точек точек струны на текущем временном слое.
//	tVector solve2 = new datatype[n]; //массив точек точек струны наследующем временном слое.
//	for (int i = 0; i < n - 1; i++)
//	{
//		solve0[i] = f(i * dh, 0);  // начальные условия, форма струны при t=0
//		OUTu << solve0[i] << ",";
//	}
//	solve0[n - 1] = f((n - 1) * dh, 0);
//	OUTu << solve0[n - 1] << "}," << endl; //запись первой строчки формы в файл
//
//
//
//	datatype lambda = 1;
//	datatype mu = 1;
//	datatype rho = 1000;
////	datatype a1 = 1;
////	datatype a2 = 1;
//	datatype alpha = 1;
//
//	//Находим форму струны на втором временном слое
//	//OUTu << "{";
//	//for (int i = 0; i < n - 1; i++) //не забываем что для f берем вторую производную
//	//{
//	//	solve1[i] = -(dt*dt/2)*(3 * lambda + 2 * mu)*((T_fe(x[i], 0 - dt) - 2 * T_fe(x[i], 0) + T_fe(x[i], 0 + dt) / (dt*dt)));
//	//	OUTu << solve1[i] << ",";
//	//}
//	//solve1[n - 1] = -(3 * lambda + 2 * mu)*((T_fe(x[n-1], 0 - dt) - 2 * T_fe(x[n - 1], 0) + T_fe(x[n - 1], 0 + dt) / (dt*dt)));
//	//OUTu << solve1[n - 1] << "}," << endl;
//
//	//Находим форму струны на втором временном слое
//	OUTu << "{";
//	for (int i = 0; i < n - 1; i++) //не забываем что для f берем вторую производную
//	{
//		solve1[i] = -(dt*dt / 2)*(3 * lambda + 2 * mu)*(-sin(x[i] * M_PI / 10)*cos(0));
//		OUTu << solve1[i] << ",";
//	}
//	solve1[n - 1] = -(dt*dt / 2)*(3 * lambda + 2 * mu)*(-sin(x[n - 1] * M_PI / 10)*cos(0));
//	OUTu << solve1[n - 1] << "}," << endl;
//
//
//	//сборка матриц
//	//матрица масс диагональная
//	tVector M = new datatype[n];
//	for (int i = 0; i < n; i++)
//	{
//		M[i] = rho * dh / 2;
//	}
//	//вроде ок
//
//
//	//жесткость
//	tVector K_l = new datatype[n]; //поддиагональ
//	for (int i = 0; i < n - 1; i++)
//	{
//		K_l[i] = -(lambda + 2 * mu)*(-1 / dh);
//	}
//	tVector K = new datatype[n]; //диагональ
//	for (int i = 1; i < n - 1; i++)
//	{
//		K[i] = -(lambda + 2 * mu) * 2 / dh;
//	}
//	K[0] = 1 / dh;
//	K[n - 1] = 1 / dh;
//
//	tVector K_u = new datatype[n]; //наддиагональ
//	for (int i = 0; i < n - 1; i++)
//	{
//		K_u[i] = -(lambda + 2 * mu)*(-1 / dh);
//	}
//	//проверяй
//
//
//
//
//
//
//	//имеем решение на первых двух слоях
//	//все готово для поиска решения на остальных слоях
//	//хуярим итерационный процесс, не забываем про условия на концах.
//
//	datatype current_t = dt; //так как уже на двух слоях мы уже нашли
//
//	while (current_t <= 7.01)
//	{
//		//вектор нагрузки(считаем на каждом временном шаге так как иеем косинус от т)
//		tVector F = new datatype[n];
//		F[0] = alpha * rho*(3 * lambda + 2 * mu)*(10 / (dh*M_PI* M_PI))*(cos(current_t)*(-dh * M_PI + 10 * sin((dh*M_PI) / 10)));
//		F[n - 1] = alpha * rho*(3 * lambda + 2 * mu)*(10 / (dh*M_PI* M_PI))*(cos(current_t)*(-dh * M_PI + 10 * sin((dh*M_PI) / 10)));
//		for (int i = 1; i < n - 1; i++)
//		{
//			F[i] = alpha * rho*(3 * lambda + 2 * mu)*(2 * 10 * 10 / (dh*M_PI* M_PI))*((cos((dh*M_PI) / 10) - 1)*cos(current_t)*sin((x[i] * M_PI) / 10));
//		}
//
//		current_t += dt;
//		cout << "current_t: " << current_t << endl;
//		//первое и последнее решение хуярим отдельно, остальные в цикле по пространству.
//
//		OUTu << "{";
//		solve2[0] = phi(current_t, 0);
//		OUTu << solve2[0] << ","; //записали решение на даноом временном слое в точке х=0
//		for (int i = 1; i < n - 1; i++)
//		{
//			solve2[i] = (dt*dt / M[i])*(K_l[i - 1] * solve1[i - 1] + K[i] * solve1[i] + K_u[i] * solve1[i + 1]) + 2 * solve1[i] - solve0[i] - F[i];
//			OUTu << solve2[i] << ",";
//		}
//		solve2[n - 1] = psi(current_t, 0);
//		if (current_t < 7.01)
//		{
//			OUTu << solve0[n - 1] << "}," << endl;//если наше время меньше 7 то будет еще запись
//		}
//		else if (current_t >= 7.01)
//		{
//			OUTu << solve0[n - 1] << "}}";//если время уже 7 или больше то записи больше не будет и все пиздец
//		}
//		DeletingVector(F);
//
//		//оп решение есть.
//		//теперь надо сделать шаг вперед по времени на шаблоне
//		//то есть сделать solve1=>solve0 и solve2=>solve1
//		//не пробуем через буфер(НО СУКА СТОИТ ПОПРОБОВАТЬ)
//		//сделаем покомпонентно
//		//ATTENTION
//		//тут может быть насрано в виде проебов памяти и прочего говна.
//		for (int i = 0; i <= n - 1; i++)
//		{
//			solve0[i] = solve1[i];
//		}
//		for (int i = 0; i <= n - 1; i++)
//		{
//			solve1[i] = solve2[i];
//		}
//
//	}
//	//Бля там в конце запятая но меня не ебет
//	OUTu.close();
//	DeletingVector(solve0);
//	DeletingVector(solve1);
//	DeletingVector(solve2);
//	DeletingVector(K);
//	DeletingVector(K_l);
//	DeletingVector(K_u);
//	DeletingVector(M);
//	DeletingVector(x);
//	DeletingVector(t);
//}

void Finite_Element(funct2 &f, funct2 &g, funct2 &phi, funct2 &psi, funct3 &T,datatype ksi, datatype d, datatype betta)
{
	char FileName[] = "x_FE.txt";
	char FileName2[] = "t_FE.txt";
	char FileName3[] = "u_FE.txt";
	ofstream OUTx(FileName);
	ofstream OUTt(FileName2);
	ofstream OUTu(FileName3);

	//хуярим в файлы в нормальном виде без ебаных е кстати
	OUTx.setf(ios::fixed);
	OUTt.setf(ios::fixed);
	OUTu.setf(ios::fixed);

	//datatype L = 10.;

	datatype dt = 0.1;
	datatype dh = 0.1;//условие Куранта |a*tao/h|<=1 выполняется
	int n = 101; //всего узлов на отрезке от 0 до 1 с шагом dh  0---5000
	int m = 101; // временных слоёв на указанной мною сетке   0---7000

	//хуярим на 7 секунд ебать

	OUTx << "{";
	tVector x = new datatype[n]; //массив пространственной сетки 0---5000
	for (int i = 0; i < n - 1; i++)
	{
		x[i] = i * dh;
		OUTx << x[i] << ",";
	}
	x[n - 1] = (n - 1)*dh;
	OUTx << x[n - 1] << "}"; //запись простанственной сетки в файл
	OUTx.close();

	OUTt << "{";
	tVector t = new datatype[m];  // массив временной сетки 0---7000
	for (int i = 0; i < m - 1; i++)
	{
		t[i] = i * dt;
		OUTt << t[i] << ",";
	}
	t[m - 1] = (m - 1)*dt;
	OUTt << t[m - 1] << "}"; //запись временнОй сетки в файл
	OUTt.close();


	OUTu << "{" << endl;
	OUTu << "{";
	tVector solve0 = new datatype[n]; //массив точек точек струны на предыдущем слое.
	tVector solve1 = new datatype[n]; //массив точек точек струны на текущем временном слое.
	tVector solve2 = new datatype[n]; //массив точек точек струны наследующем временном слое.
	for (int i = 0; i < n - 1; i++)
	{
		solve0[i] = f(i * dh, 0);  // начальные условия, форма струны при t=0
		OUTu << solve0[i] << ",";
	}
	solve0[n - 1] = f((n - 1) * dh, 0);
	OUTu << solve0[n - 1] << "}," << endl; //запись первой строчки формы в файл




	OUTu << "{";
	solve1[0] = 0;
	OUTu << solve1[0] << ",";
	for (int i = 1; i < n - 1; i++) //не забываем что для f берем вторую производную
	{
		solve1[i] = 0;
		OUTu << solve1[i] << ",";
	}
	solve1[n - 1] = 0;
	OUTu << solve1[n - 1] << "}," << endl;



	//сборка матриц
	//матрица масс диагональная
	tVector M = new datatype[n];
	for (int i = 0; i < n; i++)
	{
		M[i] =  dh / 2;
	}
	//вроде ок


	//жесткость
	tVector K_l = new datatype[n]; //поддиагональ
	for (int i = 0; i < n - 1; i++)
	{
		K_l[i] = (-1 / dh);
	}
	tVector K = new datatype[n]; //диагональ
	for (int i = 1; i < n - 1; i++)
	{
		K[i] = 2 / dh;
	}
	K[0] = 1 / dh;
	K[n - 1] = 1 / dh;

	tVector K_u = new datatype[n]; //наддиагональ
	for (int i = 0; i < n - 1; i++)
	{
		K_u[i] = (-1 / dh);
	}
	//проверяй






	//имеем решение на первых двух слоях
	//все готово для поиска решения на остальных слоях
	//хуярим итерационный процесс, не забываем про условия на концах.

	datatype current_t = dt; //так как уже на двух слоях мы уже нашли

	while (current_t <= 7.01)
	{
		//вектор нагрузки(считаем на каждом временном шаге так как иеем косинус от т)
		tVector F = new datatype[n];
		F[0] = (exp(-current_t-betta*(dh-10))*(current_t-1)*(exp(dh*betta)*(1-dh*betta)-1))/(dh*betta*(exp(10*betta)-1));
		F[n - 1] = (exp(-current_t )*(current_t - 1)*(1-exp(dh*betta)+dh*betta)) / (dh*betta*(exp(10 * betta) - 1));
		for (int i = 1; i < n - 1; i++)
		{
			F[i] = (exp(-current_t-betta*(x[i]+dh-10))*(exp(dh * betta) - 1)*(exp(dh * betta) - 1)*(current_t -1))/(dh*betta*(exp(10 * betta) - 1));
		}

		current_t += dt;
		cout << "current_t: " << current_t << endl;
		//первое и последнее решение хуярим отдельно, остальные в цикле по пространству.

		OUTu << "{";
		solve2[0] = phi(current_t, 0);
		OUTu << solve2[0] << ","; //записали решение на даноом временном слое в точке х=0
		for (int i = 1; i < n - 1; i++)
		{
			solve2[i] = (dt*dt / M[i])*(K_l[i - 1] * solve1[i - 1] + K[i] * solve1[i] + K_u[i] * solve1[i + 1]) + 2 * solve1[i] - solve0[i] - F[i];
			OUTu << solve2[i] << ",";
		}
		solve2[n - 1] = psi(current_t, 0);
		if (current_t < 7.01)
		{
			OUTu << solve0[n - 1] << "}," << endl;//если наше время меньше 7 то будет еще запись
		}
		else if (current_t >= 7.01)
		{
			OUTu << solve0[n - 1] << "}}";//если время уже 7 или больше то записи больше не будет
		}
		DeletingVector(F);

		for (int i = 0; i <= n - 1; i++)
		{
			solve0[i] = solve1[i];
		}
		for (int i = 0; i <= n - 1; i++)
		{
			solve1[i] = solve2[i];
		}

	}
	OUTu.close();
	DeletingVector(solve0);
	DeletingVector(solve1);
	DeletingVector(solve2);
	DeletingVector(K);
	DeletingVector(K_l);
	DeletingVector(K_u);
	DeletingVector(M);
	DeletingVector(x);
	DeletingVector(t);
}

void Lab3_Ex1()
{
	datatype err=100;
	datatype dt=0.01;
	datatype h1, h2;
	int N1, N2;
	cout<<"Enter h1 step: ";
	cin>>h1;
	cout<<endl<<"Enter h2 step: "; //                   x2 | . . . . . .
	cin>>h2;						//					   | . . . . . .
	cout<<endl;						//					   | . . . . . .
									//  				   |___________
	N1=1./h1 +1;					//								  x1
	N2=1./h2 +1;
	tMatrix solve = new tVector[N1];
	tMatrix solvemid = new tVector[N1];
//	tMatrix solvenew = new tVector[N1];
	for(int i=0; i<N1; i++)
	{
		solve[i]=new double[N2];
		solvemid[i]= new double[N2];
//		solvenew[i]=new double[N2];
	}
	for(int i=0; i<N1; i++)
	{
		for(int j=0 ; j<N2; j++)
		{
			solve[i][j]=0;
		}
	}
	for(int i=0; i<N2; i++)
	{
		solve[0][i]=1;
		solve[N1-1][i]=1;
//		solvenew[0][i]=1;
//		solvenew[N1-1][i]=1;
		solvemid[0][i]=1;
		solvemid[N1-1][i]=1;
	}
	for(int i=0; i<N1; i++)
	{
		solve[i][0]=1;
		solve[i][N2-1]=1;
//		solvenew[i][0]=1;
//		solvenew[i][N2-1]=1;
		solvemid[i][0]=1;
		solvemid[i][N2-1]=1;
	}
	for(int i=0; i<N1; i++)
	{
		for(int j=0; j<N2; j++)
		{
			cout<<solve[i][j]<<" | ";
		}
		cout<<endl;
	}
	tVector solveboof = nullptr;
	tVector F = new double[N1 - 2];
	tVector Fpr = new double[N2 - 2];
	tVector d = new double[N1 - 2];
	tVector ld = new double[N1 - 3];
	tVector ud = new double[N1 - 3];
//	while(err>eps5)
//	{
	for(int count=0; count<500; count++)
	{
		err=0;
		for(int i=0; i<N1 - 3; i++)
		{
			ld[i]=1./(h1*h1);
			ud[i]=1./(h1*h1);
			d[i]=-2*(1./(h1*h1) + 1./dt);
		}
		d[N1-3]=-2*(1./(h1*h1) + 1./dt);
		for(int j=1; j < N2 - 1; j++)
		{
			F[0]=-(2*solve[1][j]/dt + (solve[1][j+1]-2*solve[1][j]+solve[1][j-1])/(h2*h2))-1./(h1*h1);
			for(int i=1; i<N1-3; i++)
			{
				F[i]=-(2*solve[i+1][j]/dt + (solve[i+1][j+1]-2*solve[i+1][j]+solve[i+1][j-1])/(h2*h2));
			}
			F[N1-3]=-(2*solve[N1-2][j]/dt + (solve[N1-2][j+1]-2*solve[N1-2][j]+solve[N1-2][j-1])/(h2*h2))-1./(h1*h1);
			solveboof = Shuttle(ld, d, ud, F, N1-2);
			for(int i=0; i<N1-2; i++)
			{
				solvemid[i+1][j]=solveboof[i];
			}
			DeletingVector(solveboof);
		}//найден промежуточный слой
		for(int i=0; i<N1 - 3; i++)
		{
			ld[i]=1./(h2*h2);
			ud[i]=1./(h2*h2);
			d[i]=-2*(1./(h2*h2) + 1./dt);
		}
		d[N1-3]=-2*(1./(h2*h2) + 1./dt);
		for(int i=1; i<N1-1; i++)
		{
			Fpr[0]=-(2*solvemid[i][1]/dt + (solvemid[i+1][1]-2*solvemid[i][1]+solvemid[i-1][1])/(h1*h1))-1./(h2*h2);
			for(int j=1; j<N2 -3; j++)
			{
				Fpr[j]=-(2*solvemid[i][j+1]/dt + (solvemid[i+1][j+1]-2*solvemid[i][j+1]+solvemid[i-1][j+1])/(h1*h1));
			}
			Fpr[N2-3]=-(2*solvemid[i][N2-2]/dt + (solvemid[i+1][N2-2]-2*solvemid[i][N2-2]+solvemid[i-1][N2-2])/(h1*h1))-1./(h2*h2);
			solveboof = Shuttle(ld,d,ud,Fpr,N2 - 2);
			for(int j=0; j<N2-2; j++)
			{
				solve[i][j+1]=solveboof[j];
			}
			DeletingVector(solveboof);
		}

	}
	for(int i=0; i<N1; i++)
	{
		for(int j=0; j<N2; j++)
		{
			cout<<solve[i][j]<<" | ";
		}
		cout<<endl;
	}

	for(int i=0; i<N1; i++)
	{
		delete[] solve[i];
		delete[] solvemid[i];
	}
	delete[] solve;
	delete[] solvemid;
	DeletingVector(F);
	DeletingVector(Fpr);
	DeletingVector(ld);
	DeletingVector(d);
	DeletingVector(ud);
//	}
}















































