#ifndef Indi_1cr_1al_origin_H
#define Indi_1cr_1al_origin_H

#include<cmath>
#include<iostream>

static double const PI = 4*atan(1);

class Indi_1cr_1al_origin
{
protected:

bool *bp_chromosome;
int chromosome_length;

double variables[5], angles[2], *dp_basis;
int basis_n;

double pot_energy, fitness;

int Write_Chromosome(double*,double*,int,double *);
int  Read_Chromosome(double *,double*&,double *);


void Initialise_Chromosome(int);
double BTD(const bool *,int start, int end);
double BTD(int start, int end);
void ITB(int value,int start,int len,bool *chrom);
void ITB(int value,int start,int len);
void I1TB(int start,int len,bool *chrom);
void I1TB(int start,int len);


void Unique();
double Givens_coeff(const double *vec,int k, int l, double &c, double &s);
//Givens_rotate
double Kreuz_betrag(const double *vec1,const double *vec2);
void Vec_add(double *result, double *vec1, double *vec2);
//Vec_sub
//Sort
void VtC(double **MatV,double **MatB, int lenB);
//CtV

double ** BastMat(double *Basis,int blen,int &lenB);
double * MattBas(double **Mat, int lenB, int &blen);
//SurfMin

public:
void CtV(double **MatV,double **MatB, int lenB);	//back to protected
void Givens_rotate(double **Mat, double **MatB, int lenB); //back to protected
void Vec_sub(double *result, double *vec1, double *vec2); //back to protected
int Sort(double **Mat,int len); //back to protected
//int Surface_min(double **Mat); //back to protected
double Surface_min(double **Mat); //original is with int

static int const LEN_GENEL=31;
static int const ARC_GENEL=31;
static int const BAS_GENEL=31;

Indi_1cr_1al_origin();
Indi_1cr_1al_origin(int basis_n);
Indi_1cr_1al_origin(const Indi_1cr_1al_origin &);
~Indi_1cr_1al_origin();
Indi_1cr_1al_origin & operator=(const Indi_1cr_1al_origin &);

std::ostream& Disp(std::ostream &os,int);

int Mutate(double);

void Copy(const Indi_1cr_1al_origin &);

void Initiate(int);

void Breed(const Indi_1cr_1al_origin &, Indi_1cr_1al_origin &, Indi_1cr_1al_origin &);

double Fitness();
double Fitness(double);
double Energy();
double Energy(double);

int Basis_n();
int Basis(double * &);
int Chrom(bool * &);

void Variables(double *&);

void Angles(double *&);

bool operator<(const Indi_1cr_1al_origin &) const;

void Unique(double *bas,int &blen,double *var) ;

};
#endif
