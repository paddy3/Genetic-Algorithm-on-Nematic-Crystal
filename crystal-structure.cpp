#include<iostream>
#include <iomanip> 
#include<cmath>
#include<fstream>
#include<cassert>
#include<list>

#include<Eigen/Dense>


using namespace std;
using namespace Eigen;

double const PI = 4*atan(1);

struct Axis {
	double value;
	Vector3f axis;
} ;
bool operator< (Axis in1, Axis in2) {return in1.value < in2.value;}

double GridDiff(MatrixXf,MatrixXf,int dim);
Axis Optimal_Axis(const MatrixXf, int dim, Matrix3f (*Operation)(double,Vector3f), double);

Matrix3f Rotation(double ang, Vector3f ax) {
	Matrix3f result;
	result=AngleAxisf(ang,ax);
return result;
}
int PM(int val) {if (val<=0) return -1; return 1;}

int main(int argc, char* argv[]) {
	
	double x,y,phi,psi,theta;
	int rotN[]={6,4,3,2};
	double r_max=0;

	
	if (argc != 7)
		return 1;

	int gridN=atoi(argv[6]);
cout << "gridNumber: " << gridN << endl;
	
	x=atof(argv[1]);
	y=atof(argv[2]);
	phi=atof(argv[3]);
	psi=atof(argv[4]);
	theta=atof(argv[5]);

	r_max=pow(.75/PI*gridN*pow(x,2)*y*sin(phi)*sin(theta),1./3);
cout << "rmax: " << r_max << endl;
	
	
	// part where grid points are selected
	double x2[] = {x*abs(cos(phi)),x*sin(phi)};
	double x3[] = {x*y*sin(psi)*abs(cos(theta)),
				    x*y*cos(psi)*abs(cos(theta)),x*y*sin(theta)};

	double r=0,rx=0,ry=0,rz=0;
	double rt=0;
	double xt=0,yt=0;
	double xm=0,dx=0,bj=0,rj=0;
	
	int imax=0, jmax=0, kmax=0;

	kmax=r_max/x3[2];

	gridN=0;
	for(int  k=-kmax;k<=kmax;k++) {
		rt=pow(r_max,2)-pow(k*x3[2],2);
		if (rt<0)
			continue;
		rt=sqrt(rt);

		xt= k*x3[0];
		yt= k*x3[1];

		xm=yt/tan(phi)-xt;
		dx=rt/sin(phi);

		imax=(xm+dx)+1;

		for (int i=xm-dx;i<imax;i++) {
			bj= -(x2[0]*(i+xt)+x2[1]*yt)/pow(x,2);
			rj= -(pow(i,2)+pow(k*x*y,2)+2*i*xt-pow(rt,2))/pow(x,2);

			if (rj<0)
				continue;
			
			jmax=bj+sqrt(rj);

			for (int j=bj-sqrt(rj);j<=jmax;j++) {
				rx=i+j*x2[0]+xt;
				ry=j*x2[1]+yt;
				rz=k*x3[2];

				r=sqrt(pow(rx,2)+pow(ry,2)+pow(rz,2));

				if (r>r_max || r==0)
					continue;

				gridN++;				
			}// for j (y)
		}// for i (x)		
	}// for k (z)


cout << "GridNumber: " << gridN << endl;
	
	int matC=0;
	MatrixXf ausConf(3,gridN), resConf(3,gridN);
	
	for(int  k=-kmax;k<=kmax;k++) {
		rt=pow(r_max,2)-pow(k*x3[2],2);
		if (rt<0)
			continue;
		rt=sqrt(rt);

		xt= k*x3[0];
		yt= k*x3[1];

		xm=yt/tan(phi)-xt;
		dx=rt/sin(phi);

		imax=(xm+dx)+1;

		for (int i=xm-dx;i<imax;i++) {
			bj= -(x2[0]*(i+xt)+x2[1]*yt)/pow(x,2);
			rj= -(pow(i,2)+pow(k*x*y,2)+2*i*xt-pow(rt,2))/pow(x,2);

			if (rj<0)
				continue;
			
			jmax=bj+sqrt(rj);

			for (int j=bj-sqrt(rj);j<=jmax;j++) {
				rx=i+j*x2[0]+xt;
				ry=j*x2[1]+yt;
				rz=k*x3[2];

				r=sqrt(pow(rx,2)+pow(ry,2)+pow(rz,2));

				if (r>r_max || r==0)
					continue;

				ausConf(0,matC)=rx;
				ausConf(1,matC)=ry;
				ausConf(2,matC)=rz;
				
				matC++;			
			}// for j (y)
		}// for i (x)		
	}// for k (z)
	// end grid point selection
	
	Axis axisOpti;
	
	axisOpti=Optimal_Axis(ausConf, gridN, Rotation, 2*PI/2);
	
	cout << "optimal axis: " << axisOpti.value << endl;
	cout << axisOpti.axis.transpose() << endl;
	
	
	
return 0;
}

Axis Optimal_Axis(const MatrixXf ausM, int dim, Matrix3f (*Operation) (double,Vector3f), double angle) {
	list<Axis> axisL;
	list<Axis>::iterator axisI;
	
	Axis axisNew;
	Matrix3f operM;
	MatrixXf resM(3,dim);
	Vector3f lastAxis;
	
	double diffSum=0;
	int loopC=0;
	
	axisNew.axis=Vector3f::UnitX();
		operM=(*Operation)(angle,axisNew.axis);
		axisNew.value=GridDiff(ausM,operM*ausM,dim);
		axisL.push_back(axisNew);		

	axisNew.axis=Vector3f::UnitY();
		operM=(*Operation)(angle,axisNew.axis);
		axisNew.value=GridDiff(ausM,operM*ausM,dim);
		axisL.push_back(axisNew);	
	
	axisNew.axis=Vector3f::UnitZ();
		operM=(*Operation)(angle,axisNew.axis);
		axisNew.value=GridDiff(ausM,operM*ausM,dim);
		axisL.push_back(axisNew);	

	
	axisL.sort();
	loopC=0;
	do {
		axisL.resize(3);
		diffSum=0;
		for (axisI=axisL.begin(); axisI != axisL.end(); axisI++)
			diffSum+=axisI->value;
			
for (axisI=axisL.begin(); axisI != axisL.end(); axisI++)
	cout << axisI->axis.transpose() << "  : " << axisI->value << endl;
	
		for (int pmC=0;pmC<4;pmC++) {
			axisI=axisL.begin();
				axisNew.axis=(1-axisI->value/diffSum)*axisI->axis;
			axisI++;
				axisNew.axis+=PM(pmC%2)*(1-axisI->value/diffSum)*axisI->axis;
			axisI++;
				axisNew.axis+=PM(pmC%3)*(1-axisI->value/diffSum)*axisI->axis;
			
			axisNew.axis/=axisNew.axis.norm();
			
			operM=(*Operation)(angle,axisNew.axis);
			axisNew.value=GridDiff(ausM,operM*ausM,dim);
			axisL.push_back(axisNew);
			
cout << "   val pmC " << pmC << ": " << axisNew.value << endl;
cout << axisNew.axis << endl;			
		} // for pmC
		
		axisL.sort();
		axisI=axisL.begin();
cout << "value loopC " << loopC << ": " << axisI->value << endl;
cout << axisI->axis << endl;

	loopC++;
	} while (loopC < 4 );
	
return *axisL.begin();		
}

double GridDiff(const MatrixXf ausM, const MatrixXf resM, int dim) {
	
	double diffSum=0,diffMin=0,diff=0;
	double diffSet=pow(10,7);
	double ausLen=0;
	
	for (int ausC=0; ausC<dim; ausC++) {
		diffMin=diffSet;
		for (int resC=0; resC<dim; resC++) {
			diff=0;
			for (int i=0; i<3; i++)
				diff+= pow(ausM(i,ausC)-resM(i,resC),2);
			diff=sqrt(diff);
			
			if (diff<diffMin)
				diffMin=diff;
		}//for resC
	ausLen=0;
	for (int i=0; i<3; i++)
		ausLen+=pow(ausM(i,ausC),2);
		
	diffSum+=diffMin/ausLen;
	}// for ausC

return diffSum;	
}
