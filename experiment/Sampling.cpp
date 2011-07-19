/*
 * Sampling.cpp
 *
 *  Created on: Jul 7, 2011
 *      Author: fenrisulfr
 */
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <fstream>
#include <iostream>
#include <cstring>

using namespace std;

double samplehkl (int h, int k, int l, double size, int pkpts, double wv);
double Samplemax (double size, int pkPts, double wv);
void interp (double *xs,int nxs,double *inx, double *iny, int ninvalues, double *outvalues);
void splinterp (double *xs,int nxs,double *inx, double *iny, int nin, double *outvalues);
void linterp (double *xs,int nxs,double *inx, double *iny, int nin, double *outvalues);
bool addelement(double**myarray, double newelement, unsigned int &elements);
int sparse (double *xs, double *ys, int elements, int sparse, double * xsparse, double *ysparse);
long double squares (double *ys, int nxs, double *yinterp);

double samplehkl (int h, int k, int l, double size, int pkpts, double wv)
/*no need to make thiings more robust, but it might be fun*/
{
	return 3.1;
}
double Samplemax (double size, double wv,  int pkPts)
/*
 * For a cubic particle, gives smallest peak size, well smallest
 */
{

	double ib = .83*wv/size;
	double sigma = sqrt(ib/(2*M_PI));
	return 6*sigma/pkPts;
}
/*
 * TODO need a function to find center so we can define the start end end of peaks, but to what use that
 * that would be put I don't know
 */
double Peakcenter (double wv, double a, int h, int k, int l)
/*
 * a non-cubic version just needs more math, not hard
 */
{
	double d = sqrt ((h*h+k*k+l*l)/(a*a));
	return wv/d;
}
void interp (double *xs,int nxs,double *inx, double *iny, int ninvalues, double **outvalues)
{
	unsigned int nixes = 0;
	for (int i = 0; i<nxs; i++)
	{
		int adv = 0;
		for ( ; inx[adv]<xs[i]; adv++) ;
		addelement(outvalues,iny[adv-1]+(iny[adv]-iny[adv-1])*(xs[i]-inx[adv-1])/(inx[adv]-inx[adv-1]),nixes);
	}
	if(nixes != (unsigned int) nxs)
		cerr << "interpolation fail" << endl;
	return;
}

void splinterp (double *xs,int nxs,double *inx, double *iny, int nin, double **outvalues)
{
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *splointer = gsl_spline_alloc (gsl_interp_cspline, nin);
	gsl_spline_init(splointer, inx, iny, nin);
	for (int i = 0; i<nxs; i++)
		{
		(*outvalues)[i]=gsl_spline_eval(splointer,xs[i],acc);
		}
	gsl_spline_free(splointer);
	gsl_interp_accel_free(acc);
	return;
}
void linterp (double *xs,int nxs,double *inx, double *iny, int nin, double **outvalues)
{
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *splointer = gsl_spline_alloc (gsl_interp_linear, nin);
	gsl_spline_init(splointer, inx, iny, nin);
	for (int i = 0; i<nxs; i++)
		{
		(*outvalues)[i]=gsl_spline_eval(splointer,xs[i],acc);
		}
	gsl_spline_free(splointer);
	gsl_interp_accel_free(acc);
	return;
}
void akimterp (double *xs,int nxs,double *inx, double *iny, int nin, double **outvalues)
{
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *splointer = gsl_spline_alloc (gsl_interp_akima, nin);
	gsl_spline_init(splointer, inx, iny, nin);
	for (int i = 0; i<nxs; i++)
		{
		(*outvalues)[i]=gsl_spline_eval(splointer,xs[i],acc);
		}
	gsl_spline_free(splointer);
	gsl_interp_accel_free(acc);
	return;
}
void superakimterp (double *xs,int nxs,double *inx, double *iny, int nin)
{
	ofstream akimasuper;
	akimasuper.open("./akimas.xy");
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *splointer = gsl_spline_alloc (gsl_interp_akima, nin);
	gsl_spline_init(splointer, inx, iny, nin);
	for (int i = 0; i<10000; i++)
		{
		double x = i * (xs[nxs-1]-xs[0])/10000 + xs[0];
		akimasuper << x << "	" << gsl_spline_eval(splointer,x,acc) << endl;
		}
	gsl_spline_free(splointer);
	gsl_interp_accel_free(acc);
	akimasuper.close();
	return;
}
bool addelement(double **myarray, double newelement, unsigned int &elements)
{
	double *newarray = new double[elements+1];
	if (!newarray) return false;
	//if (elements>0)
	for(unsigned int i=0; i < elements; i++)
		newarray[i]=(*myarray)[i];

	delete [] *myarray;

	*myarray=newarray;
	(*myarray)[elements]=newelement;
//	newarray=NULL;
	elements++;
	return true;
}
int sparse (double *xs, double *ys, int elements, int spars, double **xsparse, double **ysparse)
{
	unsigned int sparsee=0;
	for (int i=0; spars*i<elements; i++)
	{
		addelement(xsparse,xs[spars*i], sparsee);
		sparsee--;
		addelement(ysparse,ys[spars*i], sparsee);
	}
	if ( elements % spars != 0)
	{
		addelement(xsparse,xs[elements-1], sparsee);
		sparsee--;
		addelement(ysparse,ys[elements-1], sparsee);
	}
	return sparsee;
}
long double squares (double *ys, int nxs, double *yinterp)
{
	long double wss=0;
	long double rts=0;
	for (int i=0; i<nxs; i++)
	{
		rts = (yinterp[i]-ys[i])*(yinterp[i]-ys[i]);
		wss += rts/ys[i];
	}
	return wss;
}
int main ()
{
	double x, y, z, a, b, c, d;
	double *xs=NULL;
	double *ys=NULL;
	double *xsparse=NULL, *ysparse=NULL, *yinterp=NULL;
	unsigned int elements = 0;
	unsigned int sparsee = 0;
	ifstream datain;

	datain.open("./test.dint");
	if(!datain)
	{
	// file couldn't be opened
	      cerr << "Error: file could not be opened" << endl;
	      return 1;
	}
	datain.ignore(250, '\n');
	datain.ignore(250, '\n');
	while (!datain.eof())
	{
		datain >> x >> z >> a >> b >> c >> d >> y;
		addelement(&xs, x, elements);

		elements--;
		addelement(&ys, y, elements);
	}
	ofstream interpout;
	cout << "readin complete" << endl;
	datain.close();
	sparsee=sparse(xs, ys, elements, 2, &xsparse, &ysparse);
	cout << "sparse complete" << endl;
	interp(xs, elements, xsparse, ysparse, sparsee, &yinterp);
	interpout.open("./linpterp.xy");
	for (int i=0; i<elements; i++)
	{
		interpout << xs[i] << "	" << yinterp[i] << endl;
	}
	interpout.close();
	cout << "lin-man " << squares(ys, elements, yinterp) << endl;
	splinterp(xs, elements, xsparse, ysparse, sparsee, &yinterp);
	interpout.open("./splinpterp.xy");
	for (int i=0; i<elements; i++)
		{
			interpout << xs[i] << "	" << yinterp[i] << endl;
		}
	interpout.close();
	cout << "spline " << squares(ys, elements, yinterp) << endl;
	linterp(xs, elements, xsparse, ysparse, sparsee, &yinterp);
	cout << "gsl-lin " << squares(ys, elements, yinterp) << endl;
	akimterp(xs, elements, xsparse, ysparse, sparsee, &yinterp);
	cout << "gsl-akima " << squares(ys, elements, yinterp) << endl;
	interpout.open("./akimterp.xy");
	for (int i=0; i<elements; i++)
		{
			interpout << xs[i] << "	" << yinterp[i] << endl;
		}
	interpout.close();
	interpout.open("./orig.xy");
	for (int i=0; i<elements; i++)
		{
			interpout << xs[i] << "	" << ys[i] << endl;
		}
	interpout.close();
	superakimterp(xs, elements, xsparse, ysparse, sparsee);
	delete [] xs;
	delete [] ys;
	delete [] xsparse;
	delete [] ysparse;
	delete [] yinterp;

	return 0;
}
