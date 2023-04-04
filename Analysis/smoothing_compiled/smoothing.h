#ifndef MYSMOOTHING
#define MYSMOOTHING

//#include <TMath.h>
//#include <Math/ProbFuncMathCore.h>
#include <iostream>
#include <cmath>
#include <limits>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
//#include "helper.h"
#include <eigen3/Eigen/Dense>

using namespace std;
namespace py = pybind11;


	template <typename VEC>
double smoothautoquad(double x, const VEC& xs, const VEC& xws, const VEC& ys, const VEC& ws, double pend)
{

	int it = 0;
	int bin = -1;
	double range = xs[xs.size()-1]-xs[0];
	double dmin = range;
	for(decltype(xs.size()) i = 0 ; i < xs.size() ; ++i)
	{   
		if(std::abs(x-xs[i]) < dmin ) {dmin = std::abs(x-xs[i]); bin = i;}
	}

	double sq = xws[bin]*2.;

	double pold = 1.;
	double step = sq*0.2;
	double u = 1.;
	Eigen::Matrix3d M;
	Eigen::Vector3d V;
	while(true)
	{
		it++;
		M.array() = 0.;
		V.array() = 0.;
		for(decltype(xs.size()) i = 0 ; i < xs.size() ; ++i)
		{
			if(ws[i] == 0.) {continue;}

			double d = std::exp(-0.5*pow((xs[i]-x)/sq, 2));
	
			double dwq = d/(ws[i]*ws[i]);	
			double xsq = xs[i]*xs[i];
			M(0,0) += dwq;
			M(0,1) += dwq*xs[i];
			M(0,2) += dwq*xsq;
			M(1,0) += dwq*xs[i];
			M(1,1) += dwq*xsq;
			M(1,2) += dwq*xsq*xs[i];
			M(2,0) += dwq*xsq;
			M(2,1) += dwq*xsq*xs[i];
			M(2,2) += dwq*xsq*xsq;

			V(0) += dwq*ys[i];
			V(1) += dwq*ys[i]*xs[i];
			V(2) += dwq*ys[i]*xsq;

		}

		Eigen::Vector3d R = M.colPivHouseholderQr().solve(V);

		double chi = 0.;
		double ndf = 0.;
		for(decltype(xs.size()) i = 0 ; i < xs.size() ; ++i)
		{
			if(ws[i] == 0.) {continue;}
			double d = std::exp(-0.5*pow((xs[i]-x)/sq, 2));
			chi += pow((R(0)+ R(1)*xs[i] + R(2)*xs[i]*xs[i] - ys[i])/ws[i], 2)*d;
			ndf+=d;
		}
        //double p = 0.0;
		double p = chi/(std::max(ys.size()/2., ndf-2.));
		//double p = ROOT::Math::chisquared_cdf_c(chi, std::max(ys.size()/2., ndf-2.));
		//pout("SMquad", it++, sq, abs(p - pend));
		if(std::abs(p - pend) < 0.0001 || sq > 5*range || sq < xws[bin] || it > 1000) 
		{
			return (R(0) +R(1)*x+ R(2)*x*x);
		}
		if((p-pend)*(pold-pend) < 0.) {step*=0.5;}
		if((p-pend)*(pold-pend) > 0.) {step*=1.2;}

		if(p-pend > 0.) {sq+=step;}
		else {sq-=step;}

		if(sq < 0.) {u++; sq = xws[bin]*0.5/u; step = sq*0.2;}

		pold = p;
	}

	return -1.;
}


	template <typename XVEC, typename VEC>
VEC smoothhist(const XVEC& bins,const VEC& y, const VEC& yw, double pend)
{
//cout << y.transpose() << endl;
//cout << yw.transpose() << endl;
    VEC x;
	x.resize(y.size());
    VEC xw;
	xw.resize(y.size());
	for(decltype(y.size()) b = 0 ; b < y.size() ; ++b)
	{
		x[b] = 0.5*(bins[b+1]+bins[b]);
		xw[b] = (bins[b+1]-bins[b]);
	}

    VEC res;
	res.resize(y.size());
	for(decltype(y.size()) b = 0 ; b < y.size() ; ++b)
	{
		if(yw[b] == 0.) {res[b] = y[b];}
		else
		{
			//res[b] = smoothauto(x[b], x, xw, y, yw, pend);
			res[b] = smoothautoquad(x[b], x, xw, y, yw, pend);
		}
	}

//cout << res.transpose() << endl;
//for(int b = 0  ; b < res.size() ; ++b) {if(!isnormal(res(b))){cout << res[b] << ": !!!!!!!!!!!!!!!!!!!!!NAN!!!!!!!!!!!!!!!!!!!!!"<< endl;}}

    return res;
}

//using py_smoothhist = smoothist < py::array, py::array >;
void py_smoothhist(const py::array& py_bins, const py::array& py_y, const py::array& py_yw, double pend, py::array& out){
    std::vector<double> bins(py_bins.size());
    std::vector<double> y(py_y.size());
    std::vector<double> yw(py_yw.size());
	for(size_t b = 0 ; b < bins.size() ; ++b) {bins[b] = static_cast<const double*>(py_bins.data())[b];}
	for(size_t b = 0 ; b < y.size() ; ++b) {
        y[b] = static_cast<const double*>(py_y.data())[b];
        yw[b] = static_cast<const double*>(py_yw.data())[b];
    }
    std::vector<double> result = smoothhist<std::vector<double>, std::vector<double> >(bins, y, yw, pend);

    //std::assert (out.size() == result.size(), "Input and output array sizes are not the same!");

    double* out_data = (double*)(void*) out.data();
	for(size_t b = 0 ; b < result.size() ; ++b) {out_data[b] = result[b];}
}


#endif
