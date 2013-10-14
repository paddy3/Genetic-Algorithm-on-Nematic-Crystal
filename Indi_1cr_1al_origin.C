#include"Indi_1cr_1al_origin.h"

#include<iostream>
#include<cmath>
#include<random>

#include<cassert>

using namespace std;



// std Constructor
Indi_1cr_1al_origin::Indi_1cr_1al_origin() : 
	bp_chromosome(0), chromosome_length(0), dp_basis(0), basis_n(0), 
	pot_energy(0), fitness(0) {}
	
void   Indi_1cr_1al_origin::Initiate(int basis) {Initialise_Chromosome(basis);}

// Constructor
Indi_1cr_1al_origin::Indi_1cr_1al_origin(int basis) : 
	bp_chromosome(0), dp_basis(0)
{
	Initialise_Chromosome(basis);
}

Indi_1cr_1al_origin::~Indi_1cr_1al_origin()
{
	delete bp_chromosome;
		bp_chromosome = 0;
	delete dp_basis;
		dp_basis = 0;
}

Indi_1cr_1al_origin::Indi_1cr_1al_origin(const Indi_1cr_1al_origin &cSource) : 
	bp_chromosome(0), dp_basis(0)  {
	Copy(cSource);	
}

Indi_1cr_1al_origin & Indi_1cr_1al_origin::operator=(const Indi_1cr_1al_origin &cSource) {
	if (this == &cSource)
		return *this;
	Copy(cSource);
return *this;
}

void Indi_1cr_1al_origin::Copy(const Indi_1cr_1al_origin &cSource) {
	chromosome_length = cSource.chromosome_length;
	
	delete bp_chromosome;
	bp_chromosome = new bool [chromosome_length];
	
	for (int i=0;i<chromosome_length;i++)
		bp_chromosome[i] = cSource.bp_chromosome[i];	
		
	for (int i=0;i<5;i++)
		variables[i] = cSource.variables[i];
	
	for (int i=0;i<2;i++)
		angles[i] = cSource.angles[i];
		
	basis_n=cSource.basis_n;
	delete dp_basis;
	dp_basis = 0;
	if (basis_n > 1) {
		dp_basis = new double [3*(basis_n-1)];
		for (int i=0;i<3*(basis_n-1);i++)
			dp_basis[i] = cSource.dp_basis[i];
	}
	
	pot_energy = cSource.pot_energy;
	fitness = cSource.fitness;

}

bool Indi_1cr_1al_origin::operator< (const Indi_1cr_1al_origin & rh) const {
	return fitness < rh.fitness;
}

std::ostream& Indi_1cr_1al_origin::Disp(std::ostream &os,int flag=0) {
	
	os << angles[0] << "  " << angles[1]<< endl;
	
	for (int i=0;i<5;i++)
		os<<variables[i]<<"  ";
	os << endl;
	
	if (flag != 1) {
		for (int i=0;i<basis_n-1;i++)
			os<<dp_basis[3*i+1] <<"  "<<dp_basis[3*i+1]<< "  " << dp_basis[3*i+1] << endl;
	}
	
	os << "Energy: " << pot_energy << " ; Fitness: " << fitness << endl;
return os;
}

void Indi_1cr_1al_origin::Initialise_Chromosome(int basis=1) {
	random_device dev;
	bernoulli_distribution bdist(0.5);
	
	if (basis <=0)
		basis = 1;
	basis_n = basis;
	
	chromosome_length = 2*LEN_GENEL+5*ARC_GENEL+(basis_n-1)*3*BAS_GENEL;
	
	delete bp_chromosome;
	bp_chromosome = new bool [chromosome_length];

	for (int i=0;i<chromosome_length;i++)
		bp_chromosome[i]= bdist(dev);

	Unique();
}

void Indi_1cr_1al_origin::Breed(const Indi_1cr_1al_origin & par, Indi_1cr_1al_origin & off1, Indi_1cr_1al_origin & off2){
	int len_par1 = chromosome_length;
	int len_par2 = par.chromosome_length;
	
	if (len_par1 != len_par2) {
			cout << "len1: " << len_par1 << endl;
			cout << "len2: " << len_par2 << endl;
	}
	
assert(len_par1 == len_par2);

	random_device rd;
	uniform_int_distribution<int> distribution(0,min(len_par1,len_par2));
	
	int cut = distribution(rd);
//cout << "cut: " << cut << endl;
		
	delete off1.bp_chromosome; off1.bp_chromosome = 0;
	delete off2.bp_chromosome; off2.bp_chromosome = 0;
	
	//using 1x from now on:
	off1.bp_chromosome = new bool[len_par2];
	off1.chromosome_length = len_par2;
	
	off2.bp_chromosome = new bool[len_par1];
	off2.chromosome_length = len_par2;
	
	off1.basis_n = par.basis_n;
	off2.basis_n = basis_n;
	
	off1.fitness = off1.pot_energy = 0;
	off2.fitness = off2.pot_energy = 0;

	for(int i=0;i<cut;i++) {
			off1.bp_chromosome[i] = bp_chromosome[i];
			off2.bp_chromosome[i] = par.bp_chromosome[i];
	}	
	
	for(int i=cut;i<len_par2;i++)
			off1.bp_chromosome[i] = par.bp_chromosome[i];			
	for(int i=cut;i<len_par2;i++)
			off2.bp_chromosome[i] = bp_chromosome[i];
			
	
	// now the general end:		
	if (! off1.Mutate(1./off1.chromosome_length))
		off1.Unique();
				
	if (! off2.Mutate(1./off2.chromosome_length))
		off2.Unique();

	if (off1.chromosome_length != len_par1 || off2.chromosome_length != len_par2) {
			cout << "len1 " << len_par1 << " lenO1 " << off1.chromosome_length;
			cout << " lenO2 " << off2.chromosome_length << endl;
	}
}

// Based on  Anderson (2000)
double Indi_1cr_1al_origin::Givens_coeff(const double *vec,int k, int l, double &c, double &s) {
	double a=vec[k], b=vec[l];
	double r=0,t=0,u=0;
	
	if (b == 0) {c = copysign(1,a); s = 0; r = abs(a);}
	else if (a == 0) {c = 0; s = -copysign(1,b); r = abs(b);}
	else if (abs(b) > abs(a)) {
	  t = a/b;
	  u = copysign(sqrt(1+t*t),b);
	  s = -1/u;
	  c = -s*t;
	  r = b*u;
	}  
	else {
	  t = b/a;
	  u = copysign(sqrt(1+t*t),a);
	  c = 1/u;
	  s = -c*t;
	  r = a*u;
	}
return r;	
}

// Givens rotation on a 3x3 Matrix
void Indi_1cr_1al_origin::Givens_rotate(double **Mat, double **MatB=0, int lenB=0) {
double	prod=0;
double c=0,s=0;

	for (int vec=0;vec<2;vec++) {
		for (int l=2;l>vec;l--) {
			if (Mat[vec][l] == 0)
				continue;
			Mat[vec][vec] = Givens_coeff(Mat[vec],vec,l,c,s);	
			Mat[vec][l] = 0;
			
			for (int rem_vec=vec+1;rem_vec<3;rem_vec++) {
				prod = c*Mat[rem_vec][vec]-s*Mat[rem_vec][l];
				Mat[rem_vec][l] = s*Mat[rem_vec][vec]+c*Mat[rem_vec][l];
				Mat[rem_vec][vec] = prod;
			}// for rem_vec		
			
			if (! (MatB && lenB) ) continue;
			for (int vecB=0; vecB<lenB;vecB++) {
				prod = c*MatB[vecB][vec]-s*MatB[vecB][l];
				MatB[vecB][l] = s*MatB[vecB][vec]+c*MatB[vecB][l];
				MatB[vecB][vec] = prod;
			}
		} // for l
	} // for vec
}

double Indi_1cr_1al_origin::Kreuz_betrag(const double *vec1,const double *vec2) {
	double sum=0;
	for (int i=0;i<2;i++)
		for (int j=i+1;j<=2;j++)
			 sum += pow(vec1[i]*vec2[j]-vec1[j]*vec2[i],2);
	
return sqrt(sum);
}

// Adds two vectors and saves result in result-vector
void Indi_1cr_1al_origin::Vec_add(double *result, double *vec1, double *vec2=0) {
	if ( ! result)
		return;
	if 	(! vec2) {
		for (int i=0;i<3;i++)
			result[i]=vec1[i];
		return;
	}
	if (! vec1) {
		for (int i=0;i<3;i++)
			result[i]=vec2[i];
		return;
	}
	for (int i=0;i<3;i++)
		result[i] = vec1[i] + vec2[i];
}

// Subtracts two vectors (vec1 - vec2) and saves result in result-vector
void Indi_1cr_1al_origin::Vec_sub(double *result, double *vec1, double *vec2) {
	if ( !(result && vec2))
		return;
	if (! vec1) {
		for (int i=0;i<3;i++)
			result[i]=-vec2[i];
		return;
	}		
	for (int i=0;i<3;i++)
		result[i] = vec1[i] - vec2[i];
}

/*
int Indi_1cr_1al_origin::Surface_min(double **Mat) {
	int nr_chg=0;
	double vec_buf[3],surf_buf=0,surface[3];
	bool change=false;
	
	// initiate surface
	for (int i=0;i<3;i++) {
		surface[i]=0;
		for (int j=0;j<3;j++)
			if (i != j)
				surface[i] += Kreuz_betrag(Mat[i],Mat[j]);
	}
	
	do {
	change = false;
		for (int vec=0;vec<3;vec++) {
			for (int add_vec=0;add_vec<3;add_vec++) {
				if (add_vec==vec) continue;
				for (int b=0;b<2;b++) {
					if (b)
						Vec_add(vec_buf,Mat[vec],Mat[add_vec]);
					else
						Vec_sub(vec_buf,Mat[vec],Mat[add_vec]);
					
					surf_buf=0;
					for (int i=0;i<3;i++)
						if (i != vec)
							surf_buf += Kreuz_betrag(vec_buf,Mat[i]);
							
					if (surf_buf < surface[vec]) {					
						change = true;
						nr_chg++;
						Vec_add(Mat[vec],vec_buf);
						surface[vec] = surf_buf;
					}		
				} // for plus/minus
			} // for over vector to add
		} // for over vectors		
	} while (change);

return nr_chg;	
}
*/

double Indi_1cr_1al_origin::Surface_min(double **Mat) {
	int nr_chg=0,ivecC;
	double vec_buf[3],dvecC[3],surf_buf=0,surface=0;
	double *stVec, *ndVec;
	bool change=false;
	
	// initiate surface
	for (int i=0;i<3;i++) {
		for (int j=i+1;j<3;j++)
				surface += Kreuz_betrag(Mat[i],Mat[j]);
	}
	
	do {
	change = false;
		for (int vec=0;vec<3;vec++) {
			for (int add_vec=0;add_vec<3;add_vec++) {
				if (add_vec==vec) continue;
				for (int b=0;b<2;b++) {
					if (b)
						Vec_add(vec_buf,Mat[vec],Mat[add_vec]);
					else
						Vec_sub(vec_buf,Mat[vec],Mat[add_vec]);
					
					surf_buf=0;
					for (int i=0;i<3;i++) {
						if (i==vec)
							stVec=vec_buf;
						else
							stVec=Mat[i];
							
						for (int j=i+1;j<3;j++) {
							if (j==vec)
								ndVec=vec_buf;
							else
								ndVec=Mat[j];
								
							surf_buf += Kreuz_betrag(stVec,ndVec);
						}	
					}
							
					if (surf_buf < surface) {					
						change = true;
						surface = surf_buf;
						ivecC=vec;
						Vec_add(dvecC,vec_buf);
					}		
				} // for plus/minus
			} // for over vector to add
		} // for over vectors
		if (change) {
			Vec_add(Mat[ivecC],dvecC);
			nr_chg++;
		}		
	} while (change);

return surface;	
}


int Indi_1cr_1al_origin::Sort(double **Mat,int len=3) {
	if (! Mat) return 0;
	
	double *vec_buf=0;
	double length[len], length_buf=0;
	int chg=0, sort=0;
	
	for (int i=0;i<len;i++) {
		length[i] =0;
		for (int j=0;j<3;j++) 
			length[i] += Mat[i][j]*Mat[i][j];			
	}
	
	for (int i=1;i<len;i++) {
		for (chg=i-1;chg>=0;chg--) 
			if (length[i] <= length[chg]) break;
		if (i == ++chg) continue;
		sort++;
		length_buf = length[i];
		vec_buf = Mat[i];
		for (int j=i;j>chg;j--) {
			length[j] = length[j-1];
			Mat[j] = Mat[j-1];
		}
		length[chg] = length_buf;
		Mat[chg] = vec_buf;
	}
return sort;
}

void Indi_1cr_1al_origin::VtC(double **MatV,double **MatB, int lenB) {
	if (! (MatB && MatV)) return;
	if (! lenB) return;
	double sum=0;
	for (int i=0;i<lenB;i++) {
		for (int j=0;j<3;j++) {
			sum=0;
			for (int s=0;s<3;s++)
				sum += MatV[s][j] * MatB[i][s];
			MatB[i][j] = sum;
		}
	}
}

void Indi_1cr_1al_origin::CtV(double **MatV,double **MatB, int lenB) {
	if (! (MatB && MatV)) return;
	if (! lenB) return;
	double sum =0,newB=0;
	for (int vecB=0;vecB<lenB;vecB++) {
		for (int coorB=2;coorB>=0;coorB--) {
			sum = 0;
			for (int s=2;s>coorB;s--) 
				sum += MatB[vecB][s] * MatV[s][coorB];
			newB = (MatB[vecB][coorB]-sum) /MatV[coorB][coorB];
			if (newB!=1) {
				if (abs(newB) > 1) newB -= (int) newB;
				if (newB < 0) newB +=1;
			}	
			MatB[vecB][coorB] = newB;
		}
	}
}
	
double ** Indi_1cr_1al_origin::BastMat(double *Basis,int blen,int &lenB) {
	if (! (Basis && blen)) return 0;
	lenB = blen/3;
	double *buf=0;
	double **Mat=new double *[lenB];
	
	for (int i=0;i<lenB;i++) {
		buf = new double[3];
		for (int j=0;j<3;j++)
			buf[j] = Basis[3*i+j];
		Mat[i] = buf;
		buf = 0;
	}
return Mat;
}
	
double * Indi_1cr_1al_origin::MattBas(double **Mat, int lenB, int &blen) {
	if (! (Mat && lenB)) return 0;
	blen = lenB*3;
	
	double *basis = new double[blen];
	
	for (int i=0;i<lenB;i++)
		for (int j=0;j<3;j++)
			basis[3*i+j] = Mat[i][j];
return basis;
}

double Indi_1cr_1al_origin::BTD(const bool * array, int start, int end) {

	if (!array) {
		cout << "BTD:: 0array übergeben" << endl;
		return 0;
	}
	if (start > end) {
		cout << "BTD:: start > end" << endl;
		return 0;
	}
	
//cout << "start: " << start << " end: " << end << endl;
	double dbl = 0;
	long int pow=1;		
	for (int i=end; i>=start; i--) {
		if (array[i]) {
//cout << "add " << i <<endl;
			dbl+= pow;
		}
		pow *= 2;
	}	

//return dbl;		
return (dbl + 1)/pow;
}

double Indi_1cr_1al_origin::BTD(int start, int end) {
	return BTD(bp_chromosome,start,end);
}

void Indi_1cr_1al_origin::ITB(int value,int start,int len,bool * chrom) {
	for (int i=start+len-1;i>=start;i--) {
		if (value) {
			chrom[i] = value % 2;
			value /= 2;
		}
		else
			chrom[i] = 0;
	}
}

void Indi_1cr_1al_origin::ITB(int value,int start,int len) {
	return ITB(value,start,len,bp_chromosome);
}

void Indi_1cr_1al_origin::I1TB(int start,int len,bool *chrom) {
	for (int i=start;i<start+len;i++)
		chrom[i] = 1;
}

void Indi_1cr_1al_origin::I1TB(int start,int len) {
	return I1TB(start,len,bp_chromosome);
}

int Indi_1cr_1al_origin::Write_Chromosome(double *var,double *bas=0,int blen=0,double *ang=0){
	double border[]={1,1,PI/2,PI,PI/2}, pw=pow(2,BAS_GENEL);
	double length[]={LEN_GENEL,LEN_GENEL,ARC_GENEL,ARC_GENEL,ARC_GENEL};
	int perm[] ={0,2,1,3,4},per=0,point=0;
	
	point=0;
	for (int i=0;i<5;i++) {
		per = perm[i];
		if (var[per] >= border[per])
			I1TB(point,length[per]);
		else
			ITB(var[per]/border[per]*(pow(2,length[per])+1)-.5,point,length[per]);
		point += length[per];
	}
	
	if (ang) {
		for (int i=0;i<2;i++) {
			if (ang[i]>=PI/2)
				I1TB(point,BAS_GENEL);
			else
				ITB(ang[i]/(PI/2)*(pow(2,ARC_GENEL)+1)-.5,point,ARC_GENEL);
			point += ARC_GENEL;
		}			
	}
	else
		point += 2*ARC_GENEL;

	if (! (blen && bas))
		return point;
	for (int i=0;i<blen;i++) {
		if (bas[i]>=1)
			I1TB(point,BAS_GENEL);
		else
			ITB(bas[i]*(pw+1)-.5,point,BAS_GENEL);
		point += BAS_GENEL;
	}
return point;	
}

int Indi_1cr_1al_origin::Read_Chromosome(double *var,double *&bas,double *ang=0) {
	double border[]={1,1,PI/2,PI,PI/2};
	double length[]={LEN_GENEL,LEN_GENEL,ARC_GENEL,ARC_GENEL,ARC_GENEL};
	int perm[] ={0,2,1,3,4},point=0,blen=0;
	
	for (int i=0;i<5;i++) {
		var[perm[i]] = border[perm[i]] * BTD(point,point+length[perm[i]]-1);
		point += length[perm[i]];
	}
	
	if (ang) {
		for (int i=0;i<2;i++) {
			ang[i] = PI/2 * BTD(point,point+ARC_GENEL-1);
			point += ARC_GENEL;
		}		
	}
	else
		point += 2*ARC_GENEL;
		
	blen = (chromosome_length-point)/BAS_GENEL;
	if (! (blen && bas)) {
		delete bas; bas=0;
		return 0;
	}
	
	delete bas;
	bas = new double[blen];
	
	for (int i=0;i<blen;i++) {
		bas[i] = BTD(point,point+BAS_GENEL-1);
		point += BAS_GENEL;
	}
	
return blen;	
}

void Indi_1cr_1al_origin::Unique() {
	int blen = Read_Chromosome(variables,dp_basis,angles);

	Unique(dp_basis,blen,variables);
	
	int chrl=Write_Chromosome(variables,dp_basis,blen);
	if (chromosome_length != chrl) {
		cout << "Unique:: Chromosomelänge stimmt nicht überein!  ";
		cout << chromosome_length << "  " << chrl << endl;
	}	
	basis_n = blen/3+1;
}

void Indi_1cr_1al_origin::Unique(double *bas,int &blen,double *var) {

//cout << "unique" << endl;	
	
	double x=var[0], y=var[1], phi=var[2], psi=var[3], theta=var[4];

	double vec1[] = {1,0,0};
	double vec2[] = {x*abs(cos(phi)),x*sin(phi),0};
	//double vec3[] = {x*y*cos(psi)*abs(cos(theta)),
	//			    x*y*sin(psi)*abs(cos(theta)),x*y*sin(theta)};
	double vec3[] = {x*y*sin(psi)*abs(cos(theta)),
				    x*y*cos(psi)*abs(cos(theta)),x*y*sin(theta)};
	double *MatV[]={vec1,vec2,vec3};

/*
cout << "mat_roh" << endl;
for (int i=0;i<3;i++) {
	for (int j=0;j<3;j++) 
		cout << MatV[i][j] <<"  " ;
cout << endl;
}
cout << endl;
*/			
	double **MatB=0;
	int lenB=0;
	MatB = BastMat(bas,blen,lenB);
	VtC(MatV,MatB,lenB);
	
	double surf = Surface_min(MatV);
	if (! Sort(MatV) && ! surf)
		return;
		
	Givens_rotate(MatV,MatB,lenB);
	
	if (MatV[1][0] < 0 )
		Vec_sub(MatV[1],0,MatV[1]);
	if (MatV[1][1] < 0 ) { // rotate 180° um x
		MatV[1][1] *= -1;
		MatV[2][1] *= -1;
		MatV[2][2] *= -1;
	}
	
	if (MatV[2][0] < 0 )
		Vec_sub(MatV[2],0,MatV[2]);
	MatV[2][2] = abs(MatV[2][2]);

	
	CtV(MatV,MatB,lenB);
	Sort(MatB,lenB);
	bas=MattBas(MatB,lenB,blen);
	
	double a=MatV[0][0];
	
	x=sqrt(pow(MatV[1][0],2)+pow(MatV[1][1],2));
	//phi=abs(asin(MatV[1][1]/(a*x)));
	phi=atan(MatV[1][1]/MatV[1][0]);

//cout << "phi: " << phi << endl;
//assert(phi<=PI/2);

if ( ! (0 <= phi && phi<=PI/2)) {
	cout << "Unique: Phi value error " << phi << endl;
	cout << "sort: " << surf << endl;
	for (int i=0;i<3;i++) {
		for (int j=0;j<3;j++) 
			cout << MatV[i][j] <<"  " ;
	cout << endl;
	}
	cout << endl;
}
assert(0 <= phi && phi<=PI/2);

	y=sqrt(pow(MatV[2][0],2)+pow(MatV[2][1],2)+pow(MatV[2][2],2));
	psi=atan2(MatV[2][0],MatV[2][1]);
	
if ( ! (0<=psi && psi <= PI) ) {
	cout << "Unique: Psi value error " << psi << endl;
	for (int i=0;i<3;i++) {
		for (int j=0;j<3;j++) 
			cout << MatV[i][j] <<"  " ;
	cout << endl;
	}
	cout << endl;
}
assert(	0<=psi && psi <= PI);


	theta=asin(MatV[2][2]/(y));
	
if ( ! (0<=theta && theta <= PI/2) ) {
	cout << "Unique: Theta value error " << theta << endl;
	for (int i=0;i<3;i++) {
		for (int j=0;j<3;j++) 
			cout << MatV[i][j] <<"  " ;
	cout << endl;
	}
	cout << endl;
}
assert(	0<=theta && theta <= PI/2);
	
	var[0]=x/a; var[1]=y/x; var[2]=phi; var[3]=psi; var[4]=theta;	
}

int Indi_1cr_1al_origin::Chrom(bool * &chrom) {	
	int len=chromosome_length;
	
	delete chrom; chrom=0;
	if (len <= 0) return 0;
	
	chrom = new bool [len];
	for (int i=0;i<len;i++)
		chrom[i] = bp_chromosome[i];
return len;		
}

int Indi_1cr_1al_origin::Basis(double * &basis) {	
	delete basis; basis=0;
	
	if (basis_n == 1) return 0;
	if (basis_n < 1) return -1;
	
	int len = 3*(basis_n-1);
	basis = new double [len];
	
	for (int i=0;i<len;i++)
		basis[i] = dp_basis[i];
return len;		
}

int Indi_1cr_1al_origin::Basis_n() {return basis_n;}

void Indi_1cr_1al_origin::Variables(double *& variable) {
	delete variable; variable=0;
//	if (chromosome_length <= 0) return;
	variable = new double[5];
	for (int i=0;i<5;i++)
		variable[i] = variables[i];
}

void Indi_1cr_1al_origin::Angles(double *& ang) {
	delete ang; ang=0;
	ang = new double[2];
	for (int i=0;i<2;i++)
		ang[i] = angles[i];
}

double Indi_1cr_1al_origin::Fitness(double fit) {return fitness=fit;}
double Indi_1cr_1al_origin::Fitness() {return fitness;}
double Indi_1cr_1al_origin::Energy(double energy) {return pot_energy=energy;}
double Indi_1cr_1al_origin::Energy() {return pot_energy;}

int Indi_1cr_1al_origin::Mutate(double pm) {
	int mutation_c=0;
	if (pm>1 || pm<0) pm=0.02;
//	random_device dev;
	default_random_engine dev;
	bernoulli_distribution bdist(pm);
	
	for(int i=0;i<chromosome_length;i++) {
		if (bdist(dev)) {
			bp_chromosome[i] = !bp_chromosome[i];
			mutation_c++;
		}	
	}
	
	if (mutation_c) {
		Unique();
	}	
return mutation_c;
}

