// ConsoleApplication4.cpp : Defines the entry point for the console application.

// Name:  Thomas Chuna

// Date: 7/29/16

// Description: 4-D Lattice Implementing Metropolis Algorithm

// The code loops through each link and updates it using the metropolis algorithm. The metropolis algorithm uses the boundary conditions and it's action

// after updating is complete you have a configuration. This configuration is then used to calculate the 

// To Do List:	Work out units in this code, fix Save NPV algorithm for 4-D cube so you can use first 64 links.

//		Create anisotropic lattice and Neighboring position algorithm use vector ideas.



//Libraries

//#include "stdafx.h"
#include <string>
#include <vector>
#include <stdio.h>    /*printf*/
#include <math.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h> 

#include <iostream>		/* For Reading and writing to txt files*/
#include <fstream>



//Quaternion Data Structure

struct Quaternion {
	double Coord[4];
} Plaquette1, Plaquette2, q, q0, q1, q2, q3, q4, q5, line, StapleTotal;

//Cruetz Algorithm
double CalculateNewA0(double B, double k);
Quaternion Complete_A0_with_Random_Avector(double a0);

//Finding Neighbooring Links
void FindNeighboorsPositions(int i, int j, int size, std::vector<std::vector<int> > & arr, int d);
void BoundaryCondition(int i, int j, int size, std::vector<std::vector<int> > & arr, int d);

//SU(2) Operations
Quaternion RandomQuaternion();
void ZeroQuaternion(Quaternion q1);

//Multiplication
Quaternion HamiltonProduct(Quaternion q1, Quaternion q2);
Quaternion HamiltonProduct(Quaternion q1, Quaternion q2, Quaternion q3);

//Addition
Quaternion QuaternionAddition(Quaternion q1, Quaternion q2);
Quaternion QuaternionAddition(Quaternion q1, Quaternion q2, Quaternion q3);
Quaternion QuaternionAddition(Quaternion q1, Quaternion q2, Quaternion q3, Quaternion q4);
Quaternion QuaternionAddition(Quaternion q1, Quaternion q2, Quaternion q3, Quaternion q4, Quaternion q5);
Quaternion QuaternionAddition(Quaternion q1, Quaternion q2, Quaternion q3, Quaternion q4, Quaternion q5, Quaternion q6);

double QuaternionDeterminant(Quaternion q1);

int main( int argc, char *argv[])
{
	//Initializations
	//Constants & Initailize RNG
	const double Pi = 3.14159265359;		// pi
	const double boltz = 1.38064e-23;		// Boltzmann Constant
//	const double B = 1 / (2.27e-6);			// "temperature" of the system
	const double B = atof(argv[1])/100.0; 		// "temperature" of the system
	const int a = 40;				// iterations per dipole
	const int size = 4;				// number of sites per latice dimension
	const int d = 4;				// number of dimensions
	srand(1);

	//Updating Links
	double Action1, Action2, oldA0, newA0, k;
	Action1 = Action2 = oldA0 = newA0 = 0;

	//Neighboor Position Vector
	std::vector < std::vector < int > >NPV;
	NPV.assign(6 * (d - 1), std::vector<int>(2));

	//2 Dimensional Vector containing 2-D Lattice Links
	std::vector < std::vector < Quaternion > >s;
	s.assign(d, std::vector< Quaternion >(static_cast<int>(pow(size, d))));

	//Initialize Links With SU(2) Element (i.e. quaternion)
	for (int i = 0; i < d; i++) {
		for (int j = 0; j < static_cast<int>(pow(size, d)); j++) {
			s[i][j] = RandomQuaternion();
		}
	}

	//open output file
	std::ofstream out1("SelectedA0.txt");
	std::ofstream app("PolyakovLoop.txt", std::ofstream::app);

	//Metropolis Algorithm: # iterations is "a" iterations per lattice link
	// i loops over possible directions
	// j loops over the sites
	for (int n = 0; n < a; n++) {
		for (int i = 0; i < d; i++) {
			for (int j = 0; j < static_cast<int>(pow(size, d) - 2); j++) {
				//Reset Parameters for good measure
				int m = 0; k = 0;
				StapleTotal.Coord[0] = StapleTotal.Coord[1] = StapleTotal.Coord[2] = StapleTotal.Coord[3] = 0;

				//Locate Neighboors
				FindNeighboorsPositions(i, j, size, NPV, d);

				//Calculate Alpha,Beta,Gamma Staple
				//need to include conjugation for directed link this has not been implemented yet.

				for (int m = 0; m < 2 * d - 3; m++) {
					ZeroQuaternion(q);
					q = HamiltonProduct(s[NPV[3 * m][0]][NPV[3 * m][1]], s[NPV[3 * m + 1][0]][NPV[3 * m + 1][1]], s[NPV[3 * m + 2][0]][NPV[3 * m + 2][1]]);
					switch (m)
					{
					case 0:	q0 = q;	case 1:	q1 = q;
					case 2: q2 = q;	case 3: q3 = q;
					case 4:	q4 = q;	case 5:	q5 = q;
					}
				}

				//Use Cruetz Algorithm to Update Link
				StapleTotal = QuaternionAddition(q0, q1, q2, q3, q4, q5);
				k = sqrt(QuaternionDeterminant(StapleTotal));
				newA0 = CalculateNewA0(B, k);
				out1 << newA0 << std::endl;
				s[i][j] = Complete_A0_with_Random_Avector(newA0);


				//Let world know you're still running
//				if (j == size*size - 1 && i == 1) {
//					printf("Completed %dth iteration \n", n);
//				}
			}
		}
		if (n % a == 0) {
			double PolyakovLoop = 0;
			for (int i = 0; i < d; i++) {
				line = s[i][0];
//				for (int m = 0; m < 4; m++) { printf("line.Coord[%d] is : %f \n ", m, line.Coord[m]); }

				for (int j = 0; j < static_cast<int>(pow(size, d) - 2); j++) {
					line = HamiltonProduct(line, s[i][j]);
					for (int m = 0; m < 6 * (d - 1); m++) {
						//why is nothing here?
					}
				}
//				printf("Polyakov Loop value in %d direction %f \n", i, 2 * line.Coord[0]);
				PolyakovLoop = PolyakovLoop + (2.0 * line.Coord[0] / (static_cast<double>(d)));
			}
//			printf("Average Polyakov Loop: %f \n", PolyakovLoop);
			app << 1/B << " " <<  PolyakovLoop << std::endl;

		}
	}
	out1.close();
	app.close();
}

// Note: For methods always declare the size because the algorithm which calculates vector size is time intensive

//Cruetz Algorithm
double CalculateNewA0(double B, double k) {
	bool test = false;
	double x = 0; double newA0 = 0; double r = 0;

	while (test == false) {
		x = exp(-2 * B*k) + (1 - exp(-2 * B*k))*(double)rand() / (double)RAND_MAX;
		newA0 = 1 + log(x) / (B*k);

		r = (double)rand() / (double)RAND_MAX;
		if (r < sqrt(1 - newA0*newA0))
		{
			test = true;
		}
	}
	return newA0;
}
Quaternion Complete_A0_with_Random_Avector(double a0)
{
	//This is the algorithm written on the chalkboard
	//during ALexei and I's Metropolis Algorithm discussion 7/25/16
	Quaternion Q;
	const double Pi = 3.14159265359;
	double s1, theta, phi;

	theta = 2 * Pi * (double)rand() / (double)RAND_MAX;
	phi = acos(2 * ((double)rand() / (double)RAND_MAX) - 1);
	s1 = sqrt(1 - a0*a0);

	Q.Coord[0] = a0;
	Q.Coord[1] = s1*sin(theta)*sin(phi);
	Q.Coord[2] = s1*cos(theta)*cos(phi);
	Q.Coord[3] = s1*cos(theta);
	return Q;
}
//FindingNeighboorPositions
void FindNeighboorsPositions(int i, int j, int size, std::vector<std::vector<int> >& arr, int d) {
	//Instance Variable Declaration
	int x = 1; int y = size; int z = size*size; int t = size*size*size;
	int LinkDirection0 = 0; int LinkDirection1 = 0; int dummy0 = 0; int dummy1 = 0; int counter = 0;

	for (int m = 0; m < d; m++) {
		if ((j % static_cast<int>(pow(size, m))) / pow(size, m - 1) == 0 || (j % static_cast<int>(pow(size, m))) / pow(size, m - 1) == (size - 1)) {
			BoundaryCondition(i, j, size, arr, d);
			return;
		}
	}

	//Main Program
	if (i == 0) { LinkDirection0 = 0; LinkDirection1 = x; }
	else if (i == 1) { LinkDirection0 = 1; LinkDirection1 = y; }
	else if (i == 2) { LinkDirection0 = 2;  LinkDirection1 = z; }
	else { LinkDirection0 = 3; LinkDirection1 = t; }
	for (int m = 0; m < d; m++) {
		if (m != i) {
			if (m == 0) { dummy0 = 0; dummy1 = x; }
			if (m == 1) { dummy0 = 1; dummy1 = y; }
			if (m == 2) { dummy0 = 2; dummy1 = z; }
			if (m == 3) { dummy0 = 3; dummy1 = t; }

			arr[0 + 6 * counter][0] = dummy0;			arr[0 + 6 * counter][1] = j + LinkDirection1;
			arr[1 + 6 * counter][0] = LinkDirection0;	arr[1 + 6 * counter][1] = j + dummy1;
			arr[2 + 6 * counter][0] = dummy0;			arr[2 + 6 * counter][1] = j;
			arr[3 + 6 * counter][0] = dummy0;			arr[3 + 6 * counter][1] = j - dummy1 + LinkDirection1;
			arr[4 + 6 * counter][0] = LinkDirection0;	arr[4 + 6 * counter][1] = j - dummy1;
			arr[5 + 6 * counter][0] = dummy0;			arr[5 + 6 * counter][1] = j - dummy1;
			counter++;
		}
	}
}
void BoundaryCondition(int i, int j, int size, std::vector<std::vector<int> >& arr, int d) {
	//declarations
	int x = 1; int y = size; int z = size*size; int t = size*size*size;
	int LinkDirection0 = 0; int LinkDirection1 = 0; int dummy0 = 0; int dummy1 = 0; int counter = 0;
	int z1 = 0; int z2 = 0;

	//
	if (i == 0) { LinkDirection0 = 0; LinkDirection1 = x; }
	else if (i == 1) { LinkDirection0 = 1; LinkDirection1 = y; }
	else if (i == 2) { LinkDirection0 = 2;  LinkDirection1 = z; }
	else { LinkDirection0 = 3; LinkDirection1 = t; }
	for (int m = 0; m < d; m++) {
		//"m != i" is there becasue you can't have a staple pointing in the link direction (only 3 loops needed)
		if (m != i) {
			if (m == 0) { dummy0 = 0; dummy1 = x; }
			if (m == 1) { dummy0 = 1; dummy1 = y; }
			if (m == 2) { dummy0 = 2; dummy1 = z; }
			if (m == 3) { dummy0 = 3; dummy1 = t; }


			arr[0 + 6 * counter][1] = ((j + LinkDirection1) % (LinkDirection1 * size)) + ((j / (LinkDirection1 * size))*LinkDirection1 * size);
			arr[1 + 6 * counter][1] = ((j + dummy1) % (dummy1 * size)) + ((j / (dummy1*size))*dummy1*size);
			arr[2 + 6 * counter][1] = j;
			z1 = j - dummy1;
			if (j - dummy1 < 0) { z1 = dummy1 * size + (j - dummy1); }
			z2 = ((z1) % (dummy1 * size)) + ((j / (dummy1 * size))*dummy1 * size);
			arr[3 + 6 * counter][1] = ((z2 + LinkDirection1) % (LinkDirection1 * size)) + ((j / (LinkDirection1 * size))*LinkDirection1 * size);
			arr[4 + 6 * counter][1] = z2;
			arr[5 + 6 * counter][1] = z2;
			counter++;
		}
	}
}

//Quaternion Operations
Quaternion RandomQuaternion()
{
	//This is the randomization algorithm proposed by
	// K. Shoemaker. "Uniform Random Rotations".
	//D. Kirk, editor, Graphics Gems 3 pg 124-132
	Quaternion Q;
	const double Pi = 3.14159265359;
	double t, tt, a1, a2, s1, s2;
	t = (double)rand() / (double)RAND_MAX;
	tt = 1.0 - t;
	a1 = 2 * Pi * (double)rand() / (double)RAND_MAX;
	a2 = 2 * Pi * (double)rand() / (double)RAND_MAX;
	s1 = sqrt(tt);
	s2 = sqrt(t);
	Q.Coord[0] = s1*cos(a1);
	Q.Coord[1] = s1*sin(a1);
	Q.Coord[2] = s2*cos(a2);
	Q.Coord[3] = s2*sin(a2);
	return Q;
}
void ZeroQuaternion(Quaternion q1) {
	q1.Coord[0] = 0;
	q1.Coord[1] = 0;
	q1.Coord[2] = 0;
	q1.Coord[3] = 0;
}
//Quaternion "Multiplication"
Quaternion HamiltonProduct(Quaternion q1, Quaternion q2, Quaternion q3)
{
	return HamiltonProduct(HamiltonProduct(q1, q2), q3);
}
Quaternion HamiltonProduct(Quaternion q1, Quaternion q2)
{
	Quaternion q3;
	q3.Coord[0] = q1.Coord[0] * q2.Coord[0] - q1.Coord[1] * q2.Coord[1] - q1.Coord[2] * q2.Coord[2] - q1.Coord[3] * q2.Coord[3];
	q3.Coord[1] = q1.Coord[0] * q2.Coord[1] + q1.Coord[1] * q2.Coord[0] + q1.Coord[2] * q2.Coord[3] - q1.Coord[3] * q2.Coord[2];
	q3.Coord[2] = q1.Coord[0] * q2.Coord[2] - q1.Coord[1] * q2.Coord[3] + q1.Coord[2] * q2.Coord[0] + q1.Coord[3] * q2.Coord[1];
	q3.Coord[3] = q1.Coord[0] * q2.Coord[3] + q1.Coord[1] * q2.Coord[2] - q1.Coord[2] * q2.Coord[1] + q1.Coord[3] * q2.Coord[0];
	return q3;
}
//Quaternion Addition
Quaternion QuaternionAddition(Quaternion q1, Quaternion q2)
{
	Quaternion Q;
	Q.Coord[0] = q1.Coord[0] + q2.Coord[0];
	Q.Coord[1] = q1.Coord[1] + q2.Coord[1];
	Q.Coord[2] = q1.Coord[2] + q2.Coord[2];
	Q.Coord[3] = q1.Coord[3] + q2.Coord[3];
	return Q;
}
Quaternion QuaternionAddition(Quaternion q1, Quaternion q2, Quaternion q3) {
	return QuaternionAddition(QuaternionAddition(q1, q2), q3);
}
Quaternion QuaternionAddition(Quaternion q1, Quaternion q2, Quaternion q3, Quaternion q4) {
	return QuaternionAddition(QuaternionAddition(QuaternionAddition(q1, q2), q3), q4);
}
Quaternion QuaternionAddition(Quaternion q1, Quaternion q2, Quaternion q3, Quaternion q4, Quaternion q5) {
	return QuaternionAddition(QuaternionAddition(QuaternionAddition(QuaternionAddition(q1, q2), q3), q4), q5);
}
Quaternion QuaternionAddition(Quaternion q1, Quaternion q2, Quaternion q3, Quaternion q4, Quaternion q5, Quaternion q6) {
	return QuaternionAddition(QuaternionAddition(QuaternionAddition(QuaternionAddition(QuaternionAddition(q1, q2), q3), q4), q5), q6);
}

double QuaternionDeterminant(Quaternion q1)
{
	return pow(q1.Coord[0], 2) + pow(q1.Coord[1], 2) + pow(q1.Coord[2], 2) + pow(q1.Coord[3], 2);
}
