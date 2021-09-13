#ifndef LJSOLVER
#define LJSOLVER

#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

using namespace std;
namespace py = pybind11;

struct Soln2d {
  Eigen::Vector3d vec;
  double k;
};


class LorentzVector {
private:
  double x_, y_, z_, e_;

public:

  LorentzVector():
    x_(0), 
    y_(0), 
    z_(0), 
    e_(0) {}
  
  LorentzVector(double x, double y, double z, double e):
    x_(x), 
    y_(y), 
    z_(z), 
    e_(e) {}

  LorentzVector(const py::array& arr) {
    double* npdata = (double*)(void*) arr.data();
    x_ = npdata[0]; 
    y_ = npdata[1]; 
    z_ = npdata[2]; 
    e_ = npdata[3];
  }

  void print() const {
    std::cout << "LV(" << x_ << ", " << y_  << ", " << z_  << ", " << e_ << ")" << std::endl;
  }
  double x() const {return x_;}
  double y() const {return y_;}
  double z() const {return z_;}
  double energy() const {return e_;}
  double mass() const {
    return std::sqrt(
	    std::max(
	     std::pow(e_, 2.) - (
	      std::pow(x_, 2.) +
	      std::pow(y_, 2.) +
	      std::pow(z_, 2.)
	     ),
	     0.d
	    )
      );
  }
};

class NuSolveLJ {
private:
  const double Mt_ = 172.5;
  const double Mw_ = 80.4;
  double Ml_ = 0.;
  double Mb_ = 0.;
  const double Mn_ = 0.;

  bool ERROR_ = false;
  bool WEIRD_ = false;
  Eigen::Matrix3d H_;
  Eigen::Vector2d MET_;
  Eigen::Matrix2d VM_;
  Eigen::Matrix3d V0_;
  Eigen::Matrix3d SIGm2_;

  LorentzVector nu_;

  bool error_ =false;
  bool solved_ = false;

  const Eigen::Matrix3d Unit_;
  static Eigen::Matrix3d make_unit() {
    Eigen::Matrix3d ret = Eigen::Matrix3d::Identity();
    ret(2,2) = -1;
    return ret;
  }
  const Eigen::Matrix3d Deriv_;
  static Eigen::Matrix3d make_deriv() {
    Eigen::Matrix3d ret;
    ret << 0.d, -1.d, 0.d,
      1.d, 0.d , 0.d,
      0.d, 0.d , 0.d;
    return ret;
  }

  inline double Cofactor(const Eigen::Matrix3d& M, int row, int col);
  Eigen::Matrix3d MatCross(const Eigen::Vector3d& L, const Eigen::Matrix3d& M);
  inline int swapI(int i);
  Eigen::Matrix3d SwapMat(const Eigen::Matrix3d& G);
  Eigen::Vector3d SwapLine(const Eigen::Vector3d& L);
  std::vector<Eigen::Vector3d> Lfactor(Eigen::Matrix3d M);
  std::vector<Soln2d> IntersectLineEllipse(const Eigen::Vector3d& line, const Eigen::Matrix3d& ellipse);
  std::vector<Soln2d> IntersectEllipses(const Eigen::Matrix3d& ellipse1, const Eigen::Matrix3d& ellipse2);

  double nu_chisq_ = numeric_limits<double>::max();

  void Solve();
  void Reset() {
    ERROR_ = false;
    WEIRD_ = false;
    solved_ = false;
    Ml_ = 0.;
    Mb_ = 0.;
    //Mn_ = 0.;

    H_ = Eigen::Matrix3d::Zero();
    MET_ = Eigen::Vector2d::Zero();
    VM_ = Eigen::Matrix2d::Zero();
    V0_ = Eigen::Matrix3d::Zero();
    SIGm2_ = Eigen::Matrix3d::Zero();

    nu_chisq_ = numeric_limits<double>::max();
  }

  bool Error() {
    return ERROR_;
  }

  double Chi2(const Eigen::Matrix3d& X, const Eigen::Vector3d& soln) {
    return (soln.transpose() * X * soln)(0,0);
  }

  void SetMET(double metx, double mety, double metxerr, double metyerr, double metxyrho);
  void BuildH(const LorentzVector& lep, const LorentzVector& bjet);
  int eps_[3][3][3];

public:

  double NuChi2() {
    if( !solved_ ) std::cout << "solved_ is false! Check if Solve() is being called!" << std::endl;
    return nu_chisq_;
  }

  LorentzVector* Nu(){
    if(solved_){
      return &nu_;
    }
    else{
      std::cout << "solved_ is false! Check if Solve() is being called!" << std::endl;
      return nullptr;
    }
  }

  void MakeSolution(const LorentzVector& bj, const LorentzVector& lep, double met_x, double met_y, double metunc_x, double metunc_y, double metxy_rho);

  NuSolveLJ(); //configuration
  ~NuSolveLJ();

};

void run_nu_solver(const py::array& lep, const py::array& jet, const py::array& met, py::array& out) {
  NuSolveLJ solver;
  LorentzVector lepton(lep);
  LorentzVector b_jet(jet);
#ifdef DEBUG_NuSolveLJ
  lepton.print();
  b_jet.print();
#endif
  double* met_data = (double*)(void*) met.data();
#ifdef DEBUG_NuSolveLJ
  std::cout << "MET(" << met_data[0] << ", " << met_data[1] << ")" << std::endl;
#endif
  solver.MakeSolution(b_jet, lepton, met_data[0], met_data[1], 1, 1, 0);
  double* out_data = (double*)(void*) out.data();
  out_data[0] = (*solver.Nu()).x();
  out_data[1] = (*solver.Nu()).y();
  out_data[2] = (*solver.Nu()).z();
  out_data[3] = solver.NuChi2();
#ifdef DEBUG_NuSolveLJ
  std::cout << "OUT(" << (*solver.Nu()).x() << ", " << (*solver.Nu()).y() << ", "
	    << (*solver.Nu()).z() << solver.NuChi2() << ")"	<< std::endl;
#endif
}

#endif
