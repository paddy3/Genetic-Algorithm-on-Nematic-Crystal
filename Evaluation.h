#include<iostream>
#include<cmath>
#include<cstdio> //only for testing, need printf
#include<cassert>

using namespace std;

template <class T>
class Evaluation {
	double rho,f,s,l;
	double fcc_pot;
	double r_min, maxPot;
	
	
/*	
double Star_Pot(double r, double ct){
if (r==0) return 0;
	double rf=sqrt(f);
	if (r<=s)
		return 5./18*pow(rf,3)*(1./(1+rf/2)-log(r/s));

return 5./18*pow(rf,3)*(s/r/(1+rf/2)*exp(-rf*(r-s)/2*s));	
}	
*/

double Star_Pot(double r, double ct){
	if(r<=s) 
		return -log(r/s)+1/(1+sqrt(l)/2);

return s/(1+sqrt(f)/2)*exp(-sqrt(f)/2*(r-s)/s)/r;
}

double Pot_Sum(double * arr, double len, int n,int i) {
	if (i >= n || i <0 || len <= n*(n+1)/2)
		return -1;
		
	int count=0,pos=0;
	int k=0;
	double pot=0;
	
	if (i==0) {
		for (k=0;k<n;k++)
			pot += arr[k];
		return pot;	
	}
	
	for (k=0;k<i;k++) {
		pot += arr[k*(2*n-k-3)/2+i];
		count++;
	}
	
	if (count >= n-1)
		return pot;
		
	pot += arr[pos=k*(2*n-k-3)/2+1+i];
	count++;
	
	for (;count<n-1;count++) {
		pot += arr[++pos];
	}

return pot;	
}

public:

Evaluation(double r=1 , double la=1, double maxPot=1./0.) : rho(r), f(64), s(1), l(la), maxPot(maxPot) {
	double *var= new double[5];
	double epsilon=pow(10,-20);
	
	r_min = 7;

//cout << "r_min: " << r_min << endl;
//cout << "pot(r_min): " << Star_Pot(r_min) <<endl;

	var[0]=var[1]=1;
	var[2]=1.57;
	var[3]=0.785;
	var[4]=0.785;
	fcc_pot=Calc_Pot(var,0,0,0);
//fcc_pot=104;
}

double Rho() {return rho;}
double Rho(double r) {return rho=r;}

double Fcc() {return fcc_pot;}

template<typename A>
void Disp_Array(A * arr,int len) {
	
	if (len <=0) return;
	assert(arr != 0);
	for (int i=0;i<len;i++) 
		cout << arr[i] << "  ";
	cout << endl;
}

template<typename A>
void Disp_Array(A *arr,int len, std::ofstream & ofs) {
	if (len <=0) return;
	assert(arr != 0);
	for (int i=0;i<len;i++) 
		ofs << arr[i] << "  ";
}

double Pot_between(const double *var, const double *basis,int blen,const double *angle=0,int b=0,int B=0) {
	if (! var) return -1;
	if (blen < 0) return -1;			
	int n = blen/3+1;
	
	double angx,angy,angz,cTheta=0;
	if (angle) {
		angx=sin(angle[0])*cos(angle[1]);
		angy=sin(angle[0])*sin(angle[1]);
		angz=cos(angle[0]);
	}
	else {
		angz=1;
		angx=angy=0;
	}
		
	double x=var[0], y=var[1], phi=var[2], psi=var[3], theta=var[4];

if ( ! (0 < phi && phi<PI))
cout << "Phi value error " << phi << endl;
	double a = pow(n*1./(rho*x*x*y*sin(phi)*sin(theta)),1./3);
	double x2[] = {a*x*abs(cos(phi)),a*x*sin(phi)};
	//double x3[] = {a*x*y*sin(psi)*abs(cos(theta)),
	//			    a*x*y*cos(psi)*abs(cos(theta)),a*x*y*sin(theta)};
	double x3[] = {a*x*y*cos(psi)*abs(cos(theta)),
				    a*x*y*sin(psi)*abs(cos(theta)),a*x*y*sin(theta)};
			
	double b1=0, b2=0, b3=0, by=0, bx=0;	  
	double B1=0, B2=0, B3=0, By=0, Bx=0;	  
		  
	double pot=0,r=0;
	double dx=0,xt=0, yt=0, xm=0, rt=0;
	int imax=0, jmax=0, kmax=0;
	double bj=0, rj=0, sqrj=0;
	
	double rx=0, ry=0, rz=0;
	
////////////
	
	if (B==0)
			B1=B2=B3=0;
	else {
			B1=basis[3*(B-1)];
			B2=basis[3*(B-1)+1];
			B3=basis[3*(B-1)+2];
	}

	Bx=B1*a+B2*x2[0]+B3*x3[0];
	By=B2*x2[1]+B3*x3[1];
	
	
	if (b==0)
		b1=b2=b3=0;
	else {
		b1=basis[3*(b-1)];
		b2=basis[3*(b-1)+1];
		b3=basis[3*(b-1)+2];
	}
	
	bx=b1*a+b2*x2[0]+b3*x3[0];
	by=b2*x2[1]+b3*x3[1];
////////////

	kmax=r_min/x3[2];
	
	for(int k=-kmax-1;k<=kmax;k++) {
		rt=pow(r_min,2)-pow((k+b3-B3)*x3[2],2);
		if(rt<0)
			continue;
		rt=sqrt(rt);
		
		xt= (k+b3)*x3[0]+b2*x2[0]+b1*a-Bx;
		yt= (k+b3)*x3[1]+b2*x2[1]-By;
		 
		xm= yt/tan(phi)-xt;
		dx= rt/sin(phi);
		
		imax= (xm+dx)/a;
			
		for(int i=(xm-dx)/a;i<=imax;i++) {
			bj= -(x2[0]*(i*a+k*x3[0]+bx-Bx)+x2[1]*(k*x3[1]+by-By))/pow(a*x,2);
			
			rj= pow(bj,2)-(pow(i*a,2)+pow(k*a*x*y,2)+pow((b3-B3)*x3[2],2)
					+2*i*k*a*x3[0]+2*(bx-Bx)*(i*a+k*x3[0])+2*(by-By)*k*x3[1]
					+pow((bx-Bx),2)+pow((by-By),2)-pow(rt,2))/pow(a*x,2);
					
			if (rj<0)
				continue;
			sqrj=sqrt(rj);
			
			jmax=bj+sqrj;
			

			for(int j=bj-sqrj;j<=jmax;j++) {
				rx=(i+b1)*a+(j+b2)*x2[0]+(k+b3)*x3[0]-Bx;
				ry=(j+b2)*x2[1]+(k+b3)*x3[1]-By;
				rz=(k+b3-B3)*x3[2];
					
				r=sqrt(pow(rx,2)+pow(ry,2)+pow(rz,2));
				
				cTheta=(rx*angx+ry*angy+rz*angz)/r;
				pot += Star_Pot(r,cTheta);	
if (pot > maxPot) return pot;
//cout << "pot " << pot << endl;										
			} // for j
		} // for i		
	} // for k	
//cout << "potBetween " << pot << endl;	
return pot;
}

double Calc_Pot(const double *var, const double *basis,int blen,const double *angle=0) {

	if (! var) return -1;
	if (blen < 0) return -1;			
	int n = blen/3+1;
	
	double pot=0;
		
	for(int B=0;B<n;B++) {
		for(int b=B+1;b<n;b++) {
				pot += Pot_between(var,basis,blen,angle,b,B);
			} //for basis
	} //for Basis
	
	pot /= n;
	pot += Pot_between(var,basis,blen,angle,0,0)/2;
	
//cout << "pot: " << pot << endl;
return pot;
}

double Eval(T &indi, int gen=1) {
	double * basis=0,*var=0,*ang=0;
	int blen=0;
	double pot=0,fitness=0;
	
	blen=indi.Basis(basis);
	indi.Variables(var);
	indi.Angles(ang);
	
	pot = Calc_Pot(var,basis,blen,ang);
	fitness = exp(1-pow(pot/fcc_pot,1+gen/10));
	indi.Energy(pot);
	indi.Fitness(fitness);
return fitness;
}

double Climb(double *var, double *basis, int blen,double *ang) {
	if (! var) return -1;
	if (blen < 0) return -1;	
	
	T indi;
	int n=blen/3+1;
	double min_fak=10;
	double eps_len=pow(10,-20), eps_arc=pow(10,-20);
	double *step=0, *pot_arr=0;
	int 	steplen=0, potlen=0;
	double *new_var=new double[5],*pot_buf=new double[n-1],*new_ang=new double[2];
	double pot=0,newpot=0, newvalue=0,pot_k[]={0,0};
	double var_border[]={1,1,PI/2,PI,PI/2};
	
	bool eps=false,nochange=false;
	
	steplen= 7+blen;
	step = new double[steplen];
	step[0]=step[1]= min_fak*pow(2,-indi.LEN_GENEL);
	step[2]=step[4]= PI/2*min_fak*pow(2,-indi.ARC_GENEL);
	step[3]=PI*min_fak*pow(2,-indi.ARC_GENEL);
	step[5]=step[6]=step[2];
	
	for (int i=steplen-blen;i<steplen;i++)
		step[i]= min_fak*pow(2,-indi.BAS_GENEL);
	
	potlen=(n*(n-1))/2+1;
	pot_arr=new double[potlen];
	
	pot = Calc_Pot(var,basis,blen,ang);
int count=0;	
	do {

		if (nochange)
			for (int i=0;i<steplen;i++)
				step[i] /= 2;
		nochange=true;		
	
		for (int i=0;i<5;i++) {
			for (int j=0;j<5;j++)
				new_var[j]=var[j];
			
			for (int k=0;k<=1;k++) {
				new_var[i]=Wrap(var[i]+pow(-1,k)*step[i],var_border[i]);
				if (new_var[i] == -1) {
					pot_k[k]=pot;
					continue;
				}	
				pot_k[k]=Calc_Pot(new_var,basis,blen,ang);
			}
			
			if (pot_k[0] < pot_k[1] && pot_k[0]<pot) {
				pot = pot_k[0];
				var[i] = Wrap(var[i]+step[i],var_border[i]);
				nochange = false;
			}
			else if (pot_k[1] < pot) {
				pot = pot_k[1];
				var[i] = Wrap(var[i]-step[i],var_border[i]);
				nochange = false;
			}
		} // first opt.
			
		for (int i=0;i<2;i++) {
			for (int j=0;j<2;j++)
				new_ang[j]=ang[j];
			
			for (int k=0;k<=1;k++) {
				new_ang[i]=Wrap(ang[i]+pow(-1,k)*step[i],PI/2);
				if (new_ang[i] == -1) {
					pot_k[k] = pot;
					continue;
				}	
				pot_k[k]=Calc_Pot(var,basis,blen,new_ang);
			}
				
			if (pot_k[0] < pot_k[1] && pot_k[0]<pot) {
				pot = pot_k[0];
				ang[i] = Wrap(ang[i]+step[i],PI/2);
				nochange = false;
			}
			else if (pot_k[1] < pot) {
				pot = pot_k[1];
				ang[i] = Wrap(ang[i]-step[i],PI/2);
				nochange = false;
			}
		} // second opt.
		
				
	if (step[1] > eps_len)
		continue;
	if (max(step[2],step[3]) > eps_arc)
		continue;
		
	eps=true;	
	} while (! eps && ++count < 5000);
cout << count << endl;	
}

double Wrap(double val, double border) {
	if (val>border)
		return border;
	if (val<=0)
		return -1;
return val;
}

};
