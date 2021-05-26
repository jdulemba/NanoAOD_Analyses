#include "NuSolveLJ.h"
//#define DEBUG_NuSolveLJ // Uncomment to unleash debugging

//init.
NuSolveLJ::NuSolveLJ():
    Unit_(make_unit()),
    Deriv_(make_deriv()) {

        for (size_t i = 0; i < 3; i++) {
            for (size_t j = 0; j < 3; j++) {
                for (size_t k = 0; k < 3; k++) {
                    if (i == j || j == k || k == i) {
                        eps_[i][j][k]= 0;
                    }
                    else if (i == 0) {
                        if (j==1) {eps_[i][j][k]= 1;}
                        else {eps_[i][j][k]= -1;}
                    }
                    else if (i == 1) {
                        if (j == 2) {eps_[i][j][k]= 1;}
                        else {eps_[i][j][k]= -1;}
                    }
                    else if (i == 2) {
                        if (j == 0) {eps_[i][j][k]= 1;}
                        else {eps_[i][j][k]= -1;}
                    }
                }
            }
        }
    }


NuSolveLJ::~NuSolveLJ()
{}

void NuSolveLJ::MakeSolution(const LorentzVector& bj, const LorentzVector& lep, double met_x, double met_y, double metunc_x, double metunc_y, double metxy_rho)
{
    Reset();
    BuildH(bj, lep);
#ifdef DEBUG_NuSolveLJ
    std::cout << endl << "H Built: " << endl << H_ << endl;
#endif
    SetMET(met_x,met_y,metunc_x,metunc_y,metxy_rho);
    Solve();
#ifdef DEBUG_NuSolveLJ
    std::cout << endl << "MET: " << endl << MET_ << endl
        << "VM_: " << endl << VM_  << endl
        << "V0_: " << endl << V0_ << endl
        << "SIGm2_: " << endl << SIGm2_ << endl
        << "nu_chisq_: " << nu_chisq_ << endl;
#endif
}


//builds the matrix H and H_perp from Burt's paper

void NuSolveLJ::BuildH(const LorentzVector& bjet, const LorentzVector& lep)
{
    Ml_ = lep.mass();
    Mb_ = bjet.mass();
    //Mn_ = 0.;
    Eigen::Vector3d l(lep.x(), lep.y(), lep.z());
    Eigen::Vector3d b(bjet.x(), bjet.y(), bjet.z());

#ifdef DEBUG_NuSolveLJ
    std::cout << endl << "L: " << endl << l << endl << "Mass: " << Ml_ << endl;
    std::cout << endl << "B: " << endl << b << endl << "Mass: " << Mb_ << endl;
#endif

    double pl = l.norm();
    double pb = b.norm();
    double El = std::sqrt(Ml_*Ml_ + pl*pl);
    double Eb = std::sqrt(Mb_*Mb_ + pb*pb);

    double cosbl = (l.dot(b)) / (pl * pb);
    double sinbl = std::sqrt(1. - cosbl * cosbl);

    double betal = pl/El;
    double betab = pb/Eb;
    double gammali = Ml_/El;

    double x0 = -0.5/El*(Mw_*Mw_ - Ml_*Ml_ - Mn_*Mn_);
    double x0p = -0.5/Eb*(Mt_*Mt_ - Mw_*Mw_ - Mb_*Mb_);
    double epsilon = (Mw_*Mw_ - Mn_*Mn_)*gammali*gammali;
    double Sx = (x0*betal - pl*gammali*gammali)/(betal*betal);
    double Sy = 1./sinbl*(x0p/betab - cosbl*Sx);

    double omega = 1./sinbl*(betal/betab - cosbl);

    double OmegaS = omega*omega - gammali*gammali;
    double Omega = std::sqrt(OmegaS);
    double x1 = Sx - (Sx + omega*Sy)/OmegaS;
    double y1 = Sy - (Sx + omega*Sy)*omega/OmegaS;
    double ZS = x1*x1*OmegaS - (Sy - omega*Sx)*(Sy - omega*Sx) - Mw_*Mw_ + x0*x0 + epsilon*epsilon;


    if(ZS < 0)
    {
        ERROR_ = true;
        solved_ =true;
        nu_ = LorentzVector();
        return;
    }

    double Z = std::sqrt(ZS);

    Eigen::Matrix3d Ht;
    Ht << Z/Omega      , 0.d, x1-pl,
       Z*omega/Omega, 0.d, y1,
       0.d          , Z  , 0.d;

    double w1 = std::atan2(l[1], l[0]); // TODO: CHECK!!!
    Eigen::AngleAxisd rot(-1.*w1, Eigen::Vector3d::UnitZ());
    Eigen::Vector3d ln = rot * l;
    Eigen::Vector3d bn = rot * b;

    double w2 = std::atan2(ln[2], ln[0]);
    rot = Eigen::AngleAxisd(w2, Eigen::Vector3d::UnitY());
    bn = rot * bn;

    double w3 = std::atan2(bn[2], bn[1]);

    Eigen::Matrix3d mat_transform;
    mat_transform = Eigen::AngleAxisd(w1, Eigen::Vector3d::UnitZ()) *
        Eigen::AngleAxisd(-1*w2, Eigen::Vector3d::UnitY()) *
        Eigen::AngleAxisd(w3, Eigen::Vector3d::UnitX());
    H_ = mat_transform * Ht;
}

//define a matrix cofactor

inline double NuSolveLJ::Cofactor(const Eigen::Matrix3d& M, int row, int col)
{
    return(M((row + 1) % 3 , (col + 1) % 3) * M((row + 2) % 3,(col + 2) % 3) - M((row + 1) % 3,(col + 2) % 3) * M((row + 2) % 3,(col + 1) % 3));
}

//Define the "cross product" of a line with a matrix, used to find intersections of line with ellipse

Eigen::Matrix3d NuSolveLJ::MatCross(const Eigen::Vector3d& L, const Eigen::Matrix3d& M)
{
    Eigen::Matrix3d D = Eigen::Matrix3d::Zero();
    for (int i = 0; i < 3; ++i)
    {
        for (int m = 0; m < 3; ++m)
        {
            for (int j = 0; j < 3; ++j)
            {
                for (int k = 0; k < 3; ++k)
                {
                    D(i,m)+=eps_[i][j][k]*L(j)*M(k,m);
                }
            }
        }
    }
    return D;
}

//tool to swap first and second indices, which sometimes you need to do
inline int NuSolveLJ::swapI(int i){
    switch(i){
        case 0 : return 1;
        case 1 : return 0;
        case 2 : return 2;
        default : return 3;
    }
}

//returns matrix with swapped first and second indices
Eigen::Matrix3d NuSolveLJ::SwapMat(const Eigen::Matrix3d& G)
{
    Eigen::Matrix3d Gs; // TODO, try to make it better
    for (int x = 0; x<3; ++x)
    {
        for (int y = 0; y < 3; ++y)
        {
            Gs(x, y) = G(swapI(x), swapI(y));
        }
    }
    return Gs;
}

//returns line with swapped first and second indices

Eigen::Vector3d NuSolveLJ::SwapLine(const Eigen::Vector3d& L)
{
    Eigen::Vector3d Ls;
    for (int x = 0; x<3; ++x)
    {
        Ls(x)=L(swapI(x));
    }
    return Ls;
}

//factors degenrate conic into outer product of two lines. These lines can then be intersected with an ellipse to get the intersection of the conic and ellipse

std::vector<Eigen::Vector3d> NuSolveLJ::Lfactor(Eigen::Matrix3d M)
{
    Eigen::Vector3d Lp;
    Eigen::Vector3d Lm;
    std::vector<Eigen::Vector3d> Ls;
    double small = 10e-10;
    bool swapQ = abs(M(0, 0)) > abs(M(1, 1));
    if (swapQ) {M = SwapMat(M);}
    M *= (1/M(1,1));
    if( abs(M(0, 0)) < small && abs(M(1, 1)) < small )
    {
        Lp = Eigen::Vector3d( M(0,1), 0, M(1,2));
        Lm = Eigen::Vector3d( 0, M(0,1), M(0,2)-M(1,2));
    }
    else if( abs( Cofactor(M, 2, 2) ) < small && abs(M(1, 1)) > small )
    {
        if( Cofactor(M, 0, 0) > 0) { return Ls; }
        Lp = Eigen::Vector3d(M(0,1), M(1,1), M(1,2) + std::sqrt( -Cofactor(M, 0, 0) ));
        Lm = Eigen::Vector3d(M(0,1), M(1,1), M(1,2) - std::sqrt( -Cofactor(M, 0, 0) ));
    }
    else if(abs(M(1, 1)) > small)
    {
        if(Cofactor(M, 2, 2) > 0) {return Ls;}
        Lp = Eigen::Vector3d( M(0,1) + std::sqrt( -Cofactor(M, 2, 2 )), M(1, 1),
                -1*(Cofactor(M, 1, 2)/Cofactor(M, 2, 2) + (Cofactor(M, 0, 2)/Cofactor(M, 2, 2)) * (M(0, 1)+std::sqrt(-Cofactor(M, 2, 2)))));
        Lm = Eigen::Vector3d(M(0, 1) - std::sqrt( -Cofactor(M, 2, 2)), M(1, 1),
                -1*(Cofactor(M, 1, 2)/Cofactor(M, 2, 2) + (Cofactor(M, 0, 2)/Cofactor(M, 2, 2)) * (M(0, 1)-std::sqrt(-Cofactor(M, 2, 2)))));
    }
    if (swapQ) { 
        Ls.push_back(SwapLine(Lp)); 
        Ls.push_back(SwapLine(Lm));
    }
    else {
        Ls.push_back(Lp); 
        Ls.push_back(Lm);
    }

    return Ls;
}

//the intersections of a line with ellipse is found by looking at eigenvectors of the 'matrix cross product'
//the parameter k has been somewhat problematic here. Burt says it should be zero for valid solutions, but sometimes it is of order 1e-3 for valid solutions

std::vector<Soln2d> NuSolveLJ::IntersectLineEllipse(const Eigen::Vector3d& L, const Eigen::Matrix3d& ellipse)
{
#ifdef DEBUG_NuSolveLJ
    std::cout<<endl<<"NuSolveLJ::IntersectLineEllipse(L, M) " << endl
        << L << endl 
        << ellipse << endl;
    std::cout << "-----" << endl;
#endif

    std::vector<Soln2d> list;
    Eigen::Matrix3d M = 0.5 * (ellipse + ellipse.transpose());

    Eigen::Matrix3d Mx = MatCross(L, M);
    Eigen::ComplexEigenSolver<Eigen::Matrix3d> eigensolver(Mx);
    auto eivals = eigensolver.eigenvalues();
    auto evecs = eigensolver.eigenvectors();

#ifdef DEBUG_NuSolveLJ
    std::cout << "Mx:" << endl << Mx << endl;
    std::cout << "eigenvalues: " << endl;
    std::cout << eivals << endl;
#endif

    for (int i = 0; i < 3; ++i)
    {
        if(abs(eivals(i).imag()) > 1e-6 || abs(eivals(i).real()) < 1e-9){ continue; }
        Eigen::Vector3d v = evecs.col(i).real();
        double k = pow(L.dot(v), 2) + pow(v.dot(M*v), 2);
#ifdef DEBUG_NuSolveLJ
        std::cout<<endl<<"k= "<<k;
#endif
        if(abs(k)>1e-2){continue;}
        else {
            struct Soln2d soln;
            soln.vec = (1/v(2))*v;
            //			soln.k = k;
            soln.k = abs(eivals(i).real());
            list.push_back(soln);
        }
    }

#ifdef DEBUG_NuSolveLJ
    std::cout<<endl<<"N_EL: "<<list.size()<<" ";
#endif
    return list;
}

//Form the Pencil, pick a value of lambda for which it is degenerate. then factor and intersect line with ellipse.

std::vector<Soln2d> NuSolveLJ::IntersectEllipses(const Eigen::Matrix3d& ellipse1, const Eigen::Matrix3d& ellipse2)
{
#ifdef DEBUG_NuSolveLJ
    std::cout << "NuSolveLJ::IntersectEllipses(e1, e2). E1: " << endl
        << ellipse1 << endl << "E2: " << ellipse2 << endl;
#endif
    Eigen::Matrix3d A = 0.5 * (ellipse1.transpose() + ellipse1);
    Eigen::Matrix3d B = 0.5 * (ellipse2 + ellipse2.transpose());
    if (abs(A.determinant()) > abs(B.determinant()))
    {
        Eigen::Matrix3d M = A;
        A = B;
        B = M;
    }

    Eigen::Matrix3d Mat = A.inverse()*B;
    auto eivals = Mat.eigenvalues();
    std::vector<Soln2d> pointlist;

    for(int evindex=0 ; evindex<3; evindex++) {
        //constexpr double minsize	= 0.001;

        if(abs(eivals(evindex).imag()) > 10e-10){ continue; }
        Eigen::Matrix3d DegMat = B - eivals(evindex).real() * A;
        std::vector<Eigen::Vector3d> lines = Lfactor(DegMat);
        for(const auto& line : lines) {			

          std::vector<Soln2d> points = IntersectLineEllipse(line, A);
            for(const auto& point : points) {
                // bool newmin =true;
                // for(size_t n=0; n<pointlist.size(); n++) {
                // 	if ((point.vec - pointlist[n].vec).Mag() < minsize) {
                // 		newmin=false;
                // 		if (points[j].k<pointlist[n].k) {
                // 			pointlist[n]=points[j];
                // 		}
                // 		continue;
                // 	}
                // }
                //if (newmin) {pointlist.push_back(points[j]);}
                pointlist.push_back(point);
            }
        }
    }
#ifdef DEBUG_NuSolveLJ
    std::cout<<endl<<"N_EE= "<<pointlist.size()<<" ";
#endif
    return pointlist;
}

//actually find the solutions

void NuSolveLJ::Solve()
{
    if (ERROR_) {
        return;
    }

    Eigen::Matrix3d Lambda = V0_ - H_;
    Eigen::Matrix3d X = Lambda.transpose() * SIGm2_ * Lambda;
    Eigen::Matrix3d M = (X * Deriv_).transpose() + X * Deriv_;

#ifdef DEBUG_NuSolveLJ
    std::cout << "NuSolveLJ::Solve()" << endl << "Lambda: " << endl << Lambda << endl;
    std::cout << "X:" << endl << X << endl;
    std::cout << "M:" << endl << M << endl;
    std::cout << "Deriv_" << endl << Deriv_ << endl;
#endif

    std::vector<Soln2d>  solns2d  = IntersectEllipses(M, Unit_);
    Eigen::Vector3d soln;

    if (solns2d.size()==0) {
        ERROR_ = true;
        WEIRD_ = true;
        nu_ = LorentzVector();
        solved_ = true;
        return;
    }
    auto best = min_element(
            solns2d.begin(), solns2d.end(),
            [X, this ](const Soln2d & a, const Soln2d & b) -> bool {
            return Chi2(X,a.vec) < Chi2(X,b.vec);
            }
            );

    nu_chisq_ = Chi2(X, best->vec);
    soln = H_ * best->vec;
    nu_ = LorentzVector(
            soln(0), soln(1), 
            soln(2), soln.norm()
            );
    solved_ = true;
}

//set the MET

void NuSolveLJ::SetMET(double metx, double mety, double metxerr, double metyerr, double metxyrho)
{
    MET_(0,0) = metx;
    MET_(1,0) = mety;

    V0_(0,2) = metx;
    V0_(1,2) = mety;

    VM_(0,0) = metxerr*metxerr;
    VM_(1,1) = metyerr*metyerr;
    VM_(1,0) = metxerr*metyerr*metxyrho;
    VM_(0,1) = metxerr*metyerr*metxyrho;

    VM_.inverse();

    //Sigm2 is the 3x3 matrix SIGMA^-2 from Burt's paper

    SIGm2_(0,0) = VM_(0,0);
    SIGm2_(0,1) = VM_(0,1);
    SIGm2_(1,0) = VM_(1,0);
    SIGm2_(1,1) = VM_(1,1);

}
