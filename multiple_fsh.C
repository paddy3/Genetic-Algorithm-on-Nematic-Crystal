#include<iostream>
#include <iomanip> 
#include<cmath>
#include<random>
#include<list>
#include<fstream>
#include<cassert>

#include<sstream>
#include<string>

#include"Indi_1cr_1al_origin.h"
#include"Evaluation.h"

using namespace std;

int main(int argc, char* argv[]) //input: r la nrun maxPot
{
	
	double r= atof(argv[1]), la= atof(argv[2]), maxPot=atof(argv[4]);
	int nrun= atoi(argv[3]);
	
	list<Indi_1cr_1al_origin> * mylist=new list<Indi_1cr_1al_origin>;
	list<Indi_1cr_1al_origin>::iterator it, ifit;
	
	Evaluation<Indi_1cr_1al_origin> eval(r,la,maxPot);
	
	int pop_size=100;
	if (pop_size%2 == 0) 
		pop_size++;
cout << "popSize: " << pop_size << endl;	
	
	double stdErr=0.005;
	
	double *ps = new double[pop_size];
	int ps_count=0;
	double fitSum=0, fit =0, bestfit=0;
	
	double *var=0,*bas=0, *ang=0;
	int bas_len=0;

	stringstream file_sstream;
	file_sstream << "Population-" << pop_size;
	file_sstream << "_R-" << r << "_la-" << la;
	file_sstream << ".dat";
	
	FILE *gp = popen("gnuplot -persist", "w");
	ofstream datFile;
	
	const string tmp = file_sstream.str();
	const char* fileName = tmp.c_str();
	
	datFile.open(fileName,ofstream::app);
	
int min_baslen=0;
double new_pot=0, min_pot=pow(10,32), min_var[5], min_ang[2], *min_bas=0;

for (int run=0;run<nrun;run++) {	
	mylist->resize(pop_size);

	for (it=mylist->begin();it!=mylist->end();it++) {
		it->Initiate(1);
		eval.Eval(*it,1);
	}
	mylist->sort();
		
	for (it=mylist->begin();it!=mylist->end();it++)
		it->Disp(cout,0);
		
///////////////

int generation=0;
for(generation=1;generation<=240;generation++) {	
cout << "\tGeneration " << generation << endl;
	
	fitSum=0; ps_count=0; bestfit=0;
	for (it=mylist->begin();it!=mylist->end();it++) {
//		if (! generation)
//			it->Initiate(1);
//cout << "calculate generation " << generation << " indi "<< ps_count+1;
//cout << " popSize " << pop_size << endl;

		if (generation > 1)
			ps[ps_count++]= fitSum+= fit =eval.Eval(*it,generation);
		else
			ps[ps_count++]= fitSum+= fit =it->Fitness();	
		if (fit > bestfit) {
			bestfit = fit;
			ifit = it;
		}
	}
assert(ps_count == pop_size);
	
	for (int i=0; i<ps_count; i++)
		ps[i] /= fitSum;
// output //
	{
		double val, sum, qsum, mini, maxi, mittelW, stdAbw;
		int n;
		n=sum=qsum=maxi= 0;
		mini=1./0.;
		for (it=mylist->begin();it!=mylist->end();) {
			val = it->Energy();
			if (val != val || val<=0) {
					it=mylist->erase(it);
cout << "deleted list element  " << val << endl;
					continue;
			}
			sum += val;
			qsum += val*val;
			maxi=max(maxi,val);
			mini=min(mini,val);
			n++;	
		it++;	
		}
		mittelW = sum/n;

cout << "qsum-sum*mittelW: " << qsum-sum*mittelW << " n: " << n<<" sum " << sum << " qsum " << qsum << endl;
assert((qsum-sum*mittelW) >= 0);
		
		stdAbw=sqrt((qsum-sum*mittelW)/(n+1));
cout <<"  --StdAbw: "<< stdAbw<< endl;		
		cout <<generation<<"\t"<<mittelW<<"\t"<<mini;
		cout <<"\t"<<mittelW-stdAbw<<"\t"<<mittelW+stdAbw<<"\t"<<maxi<<endl;

//		datFile <<generation<<"\t"<<mittelW<<"\t"<<mini;
//		datFile <<"\t"<<mittelW-stdAbw<<"\t"<<mittelW+stdAbw<<"\t"<<maxi<<endl;
	
	int genlen=1000;	
	if (generation%genlen==0) { 	
fclose(gp);
gp = popen("gnuplot -persist", "w");	
		fprintf(gp, "set bars 4.0\n");
		fprintf(gp, "set xrange [%f:%f]\n",max(0.,generation - 1.5* genlen), generation+.5*genlen);
//		fprintf(gp, "plot '%s' using 1:4:3:6:5 with candlestick whiskerbars ",fileName);
//		fprintf(gp, "lt -1 title 'evolution of population', \\\n");
//		fprintf(gp, "'' using 1:2:2:2:2 with candlestick lt -1 notitle, \\\n");
		fprintf(gp, "plot '%s'using 1:2 w l notitle, \\\n",fileName);
//		fprintf(gp, "'' using 1:2 w l notitle, \\\n");
		fprintf(gp, "'' using 1:3 w l notitle \n");
	}
		
if (generation >= 240) break;		
	} //output
	
cout << "  newGeneration" << endl;	
	{ // breed new generation
		list<Indi_1cr_1al_origin>	*buf = new list<Indi_1cr_1al_origin>;
		buf->resize(pop_size);
		list<Indi_1cr_1al_origin>::iterator ilist,ioff1, ioff2, ipar1,ipar2;
		
		random_device rd;
		uniform_real_distribution<double> distri(0.0,1.0);
		
		double ppar1=0,ppar2=0;
cout << "  **FITTEST**" << endl;

		// print fittest
		ifit->Variables(var);
		bas_len=ifit->Basis(bas);
		eval.Disp_Array(var,5);
		eval.Disp_Array(bas,bas_len);
		cout << "  " << ifit->Energy() << "  " << ifit->Fitness() << endl;		

		// overtake the fittest
		ioff1 = ioff2 = buf->begin();			
		ioff1->Copy(*ifit);

cout << "generate offspr." << endl;		
		for (int i=0;i<ps_count/2; i++) {
			if (i)
				ioff1++; 
			ioff1++;
			ioff2++; ioff2++;
			do {
				ppar1=distri(rd);
				ppar2=distri(rd);
				ipar1 = ipar2 = mylist->begin();	
				for (int j=0;j<pop_size;j++) {
					if (ppar1 > ps[j])
						ipar1++;
					if (ppar2 > ps[j])
						ipar2++;
				} //for find parents
			} while (ipar1 == ipar2);
			
			ipar1->Breed(*ipar2,*ioff1,*ioff2);
		} //for breed
		delete mylist;
		mylist = buf;		
		buf = 0;
	} //code to breed
	

}
	
	mylist->sort();
					
	for (it=mylist->begin();it!=mylist->end();it++)
		it->Disp(cout,0);	
	
cout << "optimiere..." << endl;
	Indi_1cr_1al_origin best;
	
	best=mylist->back();
	best.Variables(var);
	best.Angles(ang);
	bas_len=best.Basis(bas);
	
	eval.Climb(var,bas,bas_len,ang);
	
cout << "optimiert: " << endl;
	
	best.Unique(bas,bas_len,var);

	eval.Disp_Array(ang,2);
	eval.Disp_Array(var,5);
	eval.Disp_Array(bas,bas_len);
cout << endl;

	new_pot=eval.Calc_Pot(var,bas,bas_len,ang);	
	cout << "energy " << new_pot << endl;
	
//########	


	if (new_pot<min_pot) {
		min_pot=new_pot;
		for (int i=0;i<2;i++)
			min_ang[i] = ang[i];
		for (int i=0;i<5;i++)
			min_var[i] = var[i];
		
		min_baslen = bas_len;
		
		delete min_bas;
		min_bas = new double[min_baslen];
		for (int i=0;i<bas_len;i++)	
			min_bas[i] = bas[i];
	}
} // for runs

	cout << endl << "Gewonnene Energie:" << endl;
	eval.Disp_Array(min_ang,2);
	eval.Disp_Array(min_var,5);
	eval.Disp_Array(min_bas,min_baslen);
	cout<< setprecision (20) << "energy " << min_pot << endl;
	
	datFile << min_pot << " ";
	eval.Disp_Array(min_ang,2,datFile);
	eval.Disp_Array(min_var,5,datFile);
	eval.Disp_Array(min_bas,min_baslen,datFile);
	datFile << endl;
	
	double *basis=0;
	int blen=0;
	
	var[0]=1;
	var[1]=1;
	var[2]=1.0472;
	var[3]=1.0472;
	var[4]=0.955317;

	eval.Disp_Array(var,5);
	cout << "  " << eval.Calc_Pot(var,basis,blen) << endl;		
/*
	var[0]=1./4;
	var[1]=2./3;
	var[2]=0.02454;
	var[3]=1.37445;
	var[4]=0.0490874;
				
	eval.Disp_Array(var,5);
	cout << "  " << eval.Calc_Pot(var,basis,blen) << endl;	
*/
fclose(gp);
datFile.close();
return 0;
}
