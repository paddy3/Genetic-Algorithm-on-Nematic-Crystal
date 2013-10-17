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

double GridDiff(MatrixXf,MatrixXf);
Axis Optimal_Axis(const MatrixXf, Matrix3f (*Operation)(double,Vector3f), double, double *val =0);

Axis Optimal_Axis_Climb(const MatrixXf, Matrix3f (*Operation)(double,Vector3f), double);

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
	
	double *val=new double[5];
	for (int i=0;i<5;i++)
		val[i]=atof(argv[i+1]);
		
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
	
	Axis axisOpti,axisOptiClimb;

	for (int rotC=0; rotC<4; rotC++) {
	cout << endl << "  Rotation N: " << rotN[rotC] << endl;
		
		axisOpti=Optimal_Axis(ausConf, Rotation, 2*PI/rotN[rotC],val);
		
		cout << endl;
		cout << "optimal axis: " << axisOpti.value << endl;
		cout << axisOpti.axis.transpose() << endl;
		
		axisOptiClimb=Optimal_Axis_Climb(ausConf, Rotation, 2*PI/rotN[rotC]);
		
		cout << endl;
		cout << "optimal axis climb: " << axisOptiClimb.value << endl;
		cout << axisOptiClimb.axis.transpose() << endl;
	}	
cout << endl;
return 0;
}

Axis Optimal_Axis(const MatrixXf ausM, Matrix3f (*Operation) (double,Vector3f), double angle, double * val) {
	list<Axis> axisL;
	list<Axis>::iterator axisI;
	
	Axis axisNew;
	Matrix3f operM;
	Vector3f lastAxis;
	
	double diffSum=0;
	int loopC=0;
	
	for (int i=0;i<3;i++) {
		axisNew.axis=Vector3f::Unit(i);
			operM=(*Operation)(angle,axisNew.axis);
			axisNew.value=GridDiff(ausM,operM*ausM);
			axisL.push_back(axisNew);	
	}
	
	if (val) {
		axisNew.axis(0)=cos(val[2]);
		axisNew.axis(1)=sin(val[2]);
		axisNew.axis(2)=0;
			operM=(*Operation)(angle,axisNew.axis);
			axisNew.value=GridDiff(ausM,operM*ausM);
			axisL.push_back(axisNew);	
		
		axisNew.axis(0)=sin(val[3])*cos(val[4]);
		axisNew.axis(1)=cos(val[3])*cos(val[4]);
		axisNew.axis(2)=sin(val[4]);
			operM=(*Operation)(angle,axisNew.axis);
			axisNew.value=GridDiff(ausM,operM*ausM);
			axisL.push_back(axisNew);
	}
	
	axisL.sort();			
//for (axisI=axisL.begin(); axisI != axisL.end(); axisI++)
//	cout << axisI->axis.transpose() << "  : " << axisI->value << endl;

	loopC=0;
	do {
		axisL.resize(3);
		diffSum=0;
		for (axisI=axisL.begin(); axisI != axisL.end(); axisI++)
			diffSum+=axisI->value;
			
//for (axisI=axisL.begin(); axisI != axisL.end(); axisI++)
//	cout << axisI->axis.transpose() << "  : " << axisI->value << endl;
	
		for (int pmC=0;pmC<4;pmC++) {
			axisI=axisL.begin();
				axisNew.axis=(1-axisI->value/diffSum)*axisI->axis;
			axisI++;
				axisNew.axis+=PM(pmC%2)*(1-axisI->value/diffSum)*axisI->axis;
			axisI++;
				axisNew.axis+=PM(pmC%3)*(1-axisI->value/diffSum)*axisI->axis;
			
			axisNew.axis/=axisNew.axis.norm();
			
			operM=(*Operation)(angle,axisNew.axis);
			axisNew.value=GridDiff(ausM,operM*ausM);
			axisL.push_back(axisNew);
			
//cout << "   val pmC " << pmC << ": " << axisNew.value << endl;
//cout << axisNew.axis << endl;			
		} // for pmC
		
		axisL.sort();
		axisI=axisL.begin();
//cout << "value loopC " << loopC << ": " << axisI->value << endl;
//cout << axisI->axis << endl;

	loopC++;
	} while (loopC < 4 );
	
return *axisL.begin();		
}

Axis Optimal_Axis_Climb(const MatrixXf ausM, Matrix3f (*Operation) (double,Vector3f), double angle) {
	Matrix3f presM,minM=Matrix3f::Identity(),bufM;
	Matrix3f operM;
	
	double minDiff=pow(10,7),bufDiff=0;
	double rotAng=PI/30,minRotAng=PI/pow(10,10);
	bool diff;
	
	Matrix3f rotY,rotZ;
	rotY=AngleAxisf(rotAng, Vector3f::UnitY());
	rotZ=AngleAxisf(rotAng, Vector3f::UnitZ());
	
	Matrix3f matZ=Matrix3f::Identity();
	

	for (double angZ=0; angZ<PI/2; angZ+=rotAng) {
		bufM=matZ;
		for (double angY=0; angY<2*PI; angY+=rotAng) {
			operM=(*Operation)(angle,bufM.col(0));
			bufDiff=GridDiff(ausM,operM*ausM);
			
			if (bufDiff<minDiff) {
				minM=bufM;
				minDiff=bufDiff;
			}
		bufM=rotY*bufM;
		}//for angY
	matZ=rotZ*matZ;
	}//for angZ	
//cout << endl << "  axis pre climb: " << minDiff << endl;
//cout << minM.col(0).transpose() << endl;

rotAng/=2;

	int loopC=0;
	while(rotAng>minRotAng) {
		diff=false;
		presM=minM;
		for(int col=1;col<3;col++) {
			for(int sig=0;sig<=1;sig++) {
				bufM=AngleAxisf(pow(-1,sig)*rotAng,presM.col(col))*presM;
				
				operM=(*Operation)(angle,bufM.col(0));
				bufDiff=GridDiff(ausM,operM*ausM);
				
				if(bufDiff < minDiff) {
					minM=bufM;;
					minDiff=bufDiff;
					diff=true;
//cout << "Climb " << rotAng <<  " new min: " << minDiff << endl;
//cout << minM << endl;
				}
			}//for sig
		}//for col
	if(! diff)
		rotAng/=2;	
	loopC++;
	}// while climb rotAng
	
	Axis res;
	res.value=minDiff;
	res.axis=minM.col(0);
return res;
}

double GridDiff(const MatrixXf ausM, const MatrixXf resM) {
	
	double diffSum=0,diffMin=0,diff=0;
	double diffSet=pow(10,7);
	double ausLen=0;
	int dim=ausM.cols();
	
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
