#include<iostream>
#include <iomanip> 
#include<cmath>
#include<random>
#include<list>
#include<fstream>
#include<cassert>

#include"Indi_1cr_1al_origin.h"
#include"Evaluation.h"

using namespace std;


void Disp_Mat(double ** mat) {
	for (int i=0;i<3;i++) {
			for (int j=0;j<3;j++) 
				cout << mat[i][j] <<"  " ;
		cout << endl;
	}			
}

int main(int argc, char* argv[])
{
	Indi_1cr_1al_origin indi;
	
	Evaluation<Indi_1cr_1al_origin> eval(0.05,3);
	
	double *bas=0,*var=new double[5];
	double a=0, x=0, y=0, phi=0, psi=0, theta=0;
	int bas_len=0;
	
	switch (argc) {
		case 6:
			for (int i=0;i<5;i++)	
				var[i]=atof(argv[i+1]);
		break;
		case 10: {
			double argf[9];
			for (int i=0;i<9;i++)
				argf[i]=atof(argv[i+1]);
				
			double vec1[] = {argf[0],argf[1],argf[2]};
			double vec2[] = {argf[3],argf[4],argf[5]};
			double vec3[] = {argf[6],argf[7],argf[8]};
			
			double *MatV[]={vec1,vec2,vec3};
			
			Disp_Mat(MatV);
			cout << endl;
			
			double surf = indi.Surface_min(MatV);
			cout << "surMini: " << surf << endl;
			Disp_Mat(MatV);
			cout << endl;
				
			
			
			int sort = indi.Sort(MatV,3);
			cout << "sort: " << sort << endl;
			Disp_Mat(MatV);
			cout << endl;
			
	
			indi.Givens_rotate(MatV,0,0);
			cout << "Givens rotated" << endl;
			Disp_Mat(MatV);
			cout << endl;
			
			if (MatV[1][0] < 0 ) {
				for (int mi=0; mi<2; mi++)
					MatV[1][mi] = -MatV[1][mi];
			}
			if (MatV[1][1] < 0 ) { // rotate 180° um x
				MatV[1][1] *= -1;
				MatV[2][1] *= -1;
				MatV[2][2] *= -1;
			}

			if (MatV[2][0] < 0 )
				indi.Vec_sub(MatV[2],0,MatV[2]);
			MatV[2][2] = abs(MatV[2][2]);
	
			Disp_Mat(MatV);
			cout << endl;
			
			cout << "MatFinal" << endl;
			for (int i=0;i<3;i++) {
				for (int j=0;j<3;j++) 
					cout << MatV[i][j]/MatV[0][0] <<"  " ;
			cout << endl;
			}
			cout << endl;
			
			double a=MatV[0][0];
			
			x=sqrt(pow(MatV[1][0],2)+pow(MatV[1][1],2))/a;
			phi=atan(MatV[1][1]/MatV[1][0]);
			
			y=sqrt(pow(MatV[2][0],2)+pow(MatV[2][1],2)+pow(MatV[2][2],2))/(a*x);
			psi=abs(atan2(MatV[2][0],MatV[2][1]));
			theta=asin(MatV[2][2]/(a*x*y));
			
			var[0]=x; var[1]=y; var[2]=phi; var[3]=psi; var[4]=theta;	
		}
		break;
		default:
			cout << "Invalid input " << argc << endl;
			return 1;
	}
			
	
	eval.Disp_Array(var,5);
			
	indi.Unique(bas,bas_len,var);
	eval.Disp_Array(var,5);

	indi.Unique(bas,bas_len,var);
	eval.Disp_Array(var,5);
	
	cout << " Pot: " << eval.Calc_Pot(var,bas,bas_len) << endl;		

return 0;
}
