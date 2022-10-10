#include "libXBase/XMat3D.h"
#include "libXBase/XStringTools.h"

//-----------------------------------------------------------------------------
// Constructeurs
//-----------------------------------------------------------------------------
XMat3D::XMat3D(XPt3D axiateur)
{
    A = XPt3D(0, -axiateur.Z, axiateur.Y);
    B = XPt3D(axiateur.Z, 0, -axiateur.X);
    C = XPt3D(-axiateur.Y, axiateur.X, 0);
}

XMat3D::XMat3D(double omega,double phi,double kapa, TypeRotation type)
{
    double sin_kapa = sin(kapa);
    double sin_omega= sin(omega);
    double sin_phi	= sin(phi);
    double cos_kapa = cos(kapa);
    double cos_omega= cos(omega);
    double cos_phi	= cos(phi);

    switch(type) {
    case KapaPhiOmega:
        A=XPt3D(	cos_phi * cos_kapa,
                    cos_omega * sin_kapa + sin_phi * sin_omega * cos_kapa,
                    sin_omega * sin_kapa - sin_phi * cos_omega * cos_kapa);
        B=XPt3D(	cos_phi * sin_kapa * (-1.0),
                    cos_omega * cos_kapa - sin_phi * sin_omega * sin_kapa,
                    sin_omega * cos_kapa + sin_phi * cos_omega * sin_kapa);
        C=XPt3D(	sin_phi,
                    cos_phi * sin_omega * (-1.0),
                    cos_phi * cos_omega);
        break;
    case KapaOmegaPhi:
        A=XPt3D(	cos_phi * cos_kapa + sin_phi * sin_omega * sin_kapa,
                    cos_omega * sin_kapa,
                    sin_omega * cos_phi * sin_kapa - sin_phi * cos_kapa);
        B=XPt3D(	sin_phi * sin_omega * cos_kapa - cos_phi * sin_kapa,
                    cos_omega * cos_kapa,
                    sin_phi * sin_kapa + sin_omega * cos_phi * cos_kapa);
        C=XPt3D(	sin_phi * cos_omega,
                    sin_omega * (-1.0),
                    cos_phi * cos_omega);
        break;
    case OPK:
        A=XPt3D(	cos_phi * cos_kapa,
                    cos_phi * sin_kapa * (-1.0),
                    sin_phi);
        B=XPt3D(	cos_omega * sin_kapa + sin_phi * sin_omega * cos_kapa,
                    cos_omega * cos_kapa - sin_phi * sin_omega * sin_kapa,
                    cos_phi * sin_omega * (-1.0));
        C=XPt3D(	sin_omega * sin_kapa - sin_phi * cos_omega * cos_kapa,
                    sin_omega * cos_kapa + sin_phi * cos_omega * sin_kapa,
                    cos_phi * cos_omega);
        break;
    case KPO:
        A=XPt3D( cos_kapa * cos_phi,
                 cos_kapa * sin_phi * sin_omega - sin_kapa * cos_omega,
                 sin_kapa * sin_omega + cos_kapa * sin_phi * cos_omega);
        B=XPt3D( sin_kapa * cos_phi,
                 cos_kapa * cos_omega + sin_kapa * sin_phi * sin_omega,
                 sin_kapa * sin_phi * cos_omega - cos_kapa * sin_omega);
        C=XPt3D( sin_phi * (-1.0),
                 cos_phi * sin_omega,
                 cos_phi * cos_omega);
        break;
    }
}

//-----------------------------------------------------------------------------
// Recuperation des lignes et des colonnes
//-----------------------------------------------------------------------------
XPt3D XMat3D::lig(int n) const	// n : numero de ligne de 1 a 3
{
    switch(n) {
    case 1:
        return A;
    case 2:
        return B;
    case 3:
        return C;
    }
    return XPt3D(0,0,0);
}

XPt3D XMat3D::col(int n) const	// n : numero de colonne de 1 a 3
{
    switch(n) {
    case 1:
        return XPt3D(A.X, B.X, C.X);
    case 2:
        return XPt3D(A.Y, B.Y, C.Y);
    case 3:
        return XPt3D(A.Z, B.Z, C.Z);
    }
    return XPt3D(0,0,0);
}

//-----------------------------------------------------------------------------
// Operateurs de calcul
//-----------------------------------------------------------------------------
XMat3D& XMat3D::operator+=(XMat3D M)
{
    A += M.A;	B += M.B;	C += M.C;
    return *this;
}

XMat3D& XMat3D::operator-=(XMat3D M)
{
    A -= M.A;	B -= M.B;	C -= M.C;
    return *this;
}

XMat3D& XMat3D::operator*=(XMat3D M)
{
    XMat3D N = M.Trn();
    A = XPt3D(prodScal(A, N.A),prodScal(A, N.B),prodScal(A, N.C));
    B = XPt3D(prodScal(B, N.A),prodScal(B, N.B),prodScal(B, N.C));
    C = XPt3D(prodScal(C, N.A),prodScal(C, N.B),prodScal(C, N.C));
    return *this;
}

XMat3D& XMat3D::operator*=(double k)
{
    A *= k;		B *= k;		C *= k;
    return *this;
}

XMat3D& XMat3D::operator/=(double k)
{
    A /= k;		B /= k;		C /= k;
    return *this;
}

double XMat3D::Det() const		// Determinant
{
    return prodMixt(A,B,C);
}

XMat3D XMat3D::Trn()	const		// Transposee
{
    return XMat3D(A.X, B.X, C.X,
                  A.Y, B.Y, C.Y,
                  A.Z, B.Z, C.Z);
}

//-----------------------------------------------------------------------------
// Calcul des angles d'une matrice rotation
//-----------------------------------------------------------------------------
double XMat3D::Omega(TypeRotation type) const
{
    switch(type) {
    case KapaPhiOmega:
        return asin(C.Y * (-1.0) / cos(Phi(type)));
    case KapaOmegaPhi:
        return asin((-1.0) * C.Y);
    case OPK:
        return asin(B.Z * (-1.0) / cos(Phi(type)));
    case KPO:
        return asin(C.Y / cos(Phi(type)));
    }
    return 0.0;
}

double XMat3D::Phi(TypeRotation type) const
{
    switch(type) {
    case KapaPhiOmega:
        return asin(C.X);
    case KapaOmegaPhi:
        return asin(C.X / cos(Omega(KapaOmegaPhi)));
    case OPK:
        return asin(A.Z);
    case KPO:
        return asin(C.X*(-1.0));
    }
    return 0.0;
}

double XMat3D::Kapa(TypeRotation type) const
{
    double phi, omega;
    switch(type) {
    case KapaPhiOmega:
        phi = Phi(type);
        if (asin((-1.0)*B.X / cos(phi)) > 0)
            return acos(A.X / cos(phi));
        else
            return (-1.0) * acos(A.X / cos(phi));
    case KapaOmegaPhi:
        omega = Omega(KapaOmegaPhi);
        if (asin(A.Y / cos(omega)) > 0)
        {
            double val = B.Y / cos(omega);
            return acos(val);
        }
        else
            return (-1.0) * acos(B.Y / cos(omega));
    case OPK:
        phi = Phi(type);
        if (asin((-1.0)*A.Y / cos(phi)) > 0)
            return acos(A.X / cos(phi));
        else
            return (-1.0) * acos(A.X / cos(phi));
    case KPO:
        phi = Phi(type);
        if (asin(A.Y / cos(phi)) > 0)
            return acos(A.X / cos(phi));
        else
            return (-1.0) * acos(A.X / cos(phi));
    }
    return 0.0;
}
//-----------------------------------------------------------------------------
// Ecriture dans un fichier XML (image)
//------------------------------------------------------------------------
bool XMat3D::XmlWrite(std::ostream* out)
{
    *out << "<mat3d> " << std::endl;
    *out << "<l1> " << std::endl;
    A.XmlWrite(out);
    *out << "</l1> " << std::endl;
    *out << "<l2> " << std::endl;
    B.XmlWrite(out);
    *out << "</l2> " << std::endl;
    *out << "<l3> " << std::endl;
    C.XmlWrite(out);
    *out << "</l3> " << std::endl;
    *out << "</mat3d>" << std::endl;
    return out->good();
}
//-----------------------------------------------------------------------------
// Ecriture dans un fichier texte
//------------------------------------------------------------------------
bool XMat3D::WriteTxt(std::ostream* out)
{
    //ancienne version exportait le omega phi kappa (utilisations ????)
    //*out << Omega()*180/M_PI << "  " << Phi()*180/M_PI << "  " << Kapa()*180/M_PI <<std::endl;
    //*out << std::endl;
    *out << A << '\n';
    *out << B << '\n';
    *out << C << '\n';

    return out->good();
}
//------------------------------------------------------------------------
bool XMat3D::ReadTxt(std::ifstream* in)
{
    uint16 count =0;
    XStringTools ST;
    A = XPt3D();
    B = XPt3D();
    C = XPt3D();
    while( (in->good()) && (!in->eof()))
    {
        char ligne[1024];
        in->getline(ligne,1023);
        std::string strLigne = ligne;
        if(ST.IdentifyLineAsComment(strLigne))
            continue;
        std::vector<std::string> tokens;
        ST.Tokenize(strLigne,tokens,std::string(" \t,"));
        if(tokens.size() != 3)
            continue;
        count++;
        if(count==1)
        {
            A.X = atof(tokens[0].c_str());
            A.Y = atof(tokens[1].c_str());
            A.Z = atof(tokens[2].c_str());
            continue;
        }
        if(count==2)
        {
            B.X = atof(tokens[0].c_str());
            B.Y = atof(tokens[1].c_str());
            B.Z = atof(tokens[2].c_str());
            continue;
        }
        if(count==3)
        {
            C.X = atof(tokens[0].c_str());
            C.Y = atof(tokens[1].c_str());
            C.Z = atof(tokens[2].c_str());
            return true;
        }
    }
    return false;
}
//-----------------------------------------------------------------------------
// Operations entre matrices
//-----------------------------------------------------------------------------
XMat3D operator+(XMat3D M, XMat3D P)
{
    XMat3D Q = M;
    return Q += P;
}

XMat3D operator-(XMat3D M, XMat3D P)
{
    XMat3D Q = M;
    return Q -= P;
}

XMat3D operator*(XMat3D M, XMat3D P)
{
    XMat3D Q = M;
    return Q*= P;
}

XMat3D operator*(XMat3D M, double k)
{
    XMat3D Q = M;
    return Q *= k;
}

XMat3D operator*(double k, XMat3D M)
{
    XMat3D Q = M;
    return Q *= k;
}

XPt3D operator*(XMat3D M, XPt3D P)
{
    return XPt3D(prodScal(M.lig(1),P),prodScal(M.lig(2),P),prodScal(M.lig(3),P));
}

XPt3D operator*(XPt3D P, XMat3D M)
{
    return XPt3D(prodScal(M.col(1),P),prodScal(M.col(2),P),prodScal(M.col(3),P));
}

bool operator==(XMat3D M, XMat3D P)
{
    return M.lig(1)==P.lig(1)&& M.lig(2)==P.lig(2)&& M.lig(3)==P.lig(3);
}

bool operator!=(XMat3D M, XMat3D P)
{
    return !(M==P);
}

//-----------------------------------------------------------------------------
// Operateurs d'entree/sortie
//-----------------------------------------------------------------------------
std::istream& operator>>(std::istream& s, XMat3D& M)
{
    XPt3D A, B, C;
    s >> A >> B >> C;
    M = XMat3D(A, B, C);
    return s;
}
//------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& s, XMat3D M)
{
    s << M.lig(1) << std::endl << M.lig(2) << std::endl << M.lig(3);
    return s;
}


//Calcul Axe et angle a partir de R
//------------------------------------------------------------------------
void XMat3D::R2Axe_Angle(XPt3D* p,double* angle)const{
    XMat3D tmp=*this - this->Trn();
    double d=sqrt(tmp.B.Z*tmp.B.Z+tmp.A.Z*tmp.A.Z+tmp.A.Y*tmp.A.Y);
    if(d>0.00000001){
        p->X=-tmp.B.Z/d;
        p->Y= tmp.A.Z/d;
        p->Z=-tmp.A.Y/d;
        *angle=asin(d/2);
    }
    else{
        p->X= 0;
        p->Y= 0;
        p->Z=1;
        *angle=0;

    }

}


// Fonctions vectorielles
//------------------------------------------------------------------------
XMat3D prod_Pa_tPb(XPt3D Pa, XPt3D Pb){
    XMat3D out;
    out.A.X=Pa.X*Pb.X;out.A.Y=Pa.X*Pb.Y;out.A.Z=Pa.X*Pb.Z;
    out.B.X=Pa.Y*Pb.X;out.B.Y=Pa.Y*Pb.Y;out.B.Z=Pa.Y*Pb.Z;
    out.C.X=Pa.Z*Pb.X;out.C.Y=Pa.Z*Pb.Y;out.C.Z=Pa.Z*Pb.Z;

    return out;
}

//------------------------------------------------------------------------
XMat3D  XMat3D::R_plus_dR(double da,double db,double dc){
    /*XMat3D out=*this;
  XPt3D var_axe(da,db,dc);
  XMat3D delta(var_axe);
  XMat3D dR=out*delta;
  XMat3D mat_apres=out+dR;
  XMat3D mat_apres_debug=mat_apres;
  mat_apres.Normalise();
  */
    XPt3D axe=XPt3D(da,db,dc);
    axe.Normalise();
    XMat3D dr=Axe_Angle2R(axe,XPt3D(da,db,dc).Norme());
    XMat3D out=(*this)*dr;
    return out;
}

//Calcul R a partir de Axe et angle  
//------------------------------------------------------------------------
XMat3D XMat3D::Axe_Angle2R(XPt3D p,double angle){
    XMat3D axiateur=  XMat3D(p);
    XMat3D tmp= Identite() + axiateur * sin(angle)+ axiateur * axiateur * (1-cos(angle));
    return tmp;
}
//------------------------------------------------------------------------
void XMat3D::BinaryRead(std::ifstream& in){
    A.BinaryRead(in);
    B.BinaryRead(in);
    C.BinaryRead(in);
}
//------------------------------------------------------------------------
void XMat3D::BinaryWrite(std::ofstream& out) const {
    A.BinaryWrite(out);
    B.BinaryWrite(out);
    C.BinaryWrite(out);
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
#if USE_GSL
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
void  XMat3D::Normalise(){
    gsl_matrix * this_gsl=gsl_matrix_alloc (3, 3);
    gsl_matrix * V=gsl_matrix_alloc (3, 3);
    gsl_vector *work = gsl_vector_alloc (3);
    gsl_vector *S= gsl_vector_alloc (3);
    gsl_matrix_set (this_gsl,0,0,A.X);
    gsl_matrix_set (this_gsl,0,1,A.Y);
    gsl_matrix_set (this_gsl,0,2,A.Z);
    gsl_matrix_set (this_gsl,1,0,B.X);
    gsl_matrix_set (this_gsl,1,1,B.Y);
    gsl_matrix_set (this_gsl,1,2,B.Z);
    gsl_matrix_set (this_gsl,2,0,C.X);
    gsl_matrix_set (this_gsl,2,1,C.Y);
    gsl_matrix_set (this_gsl,2,2,C.Z);
    gsl_linalg_SV_decomp (this_gsl,V,S,work);
    gsl_matrix_transpose (V);
    A.X=gsl_matrix_get(this_gsl,0,0)*gsl_matrix_get(V,0,0)+gsl_matrix_get(this_gsl,0,1)*gsl_matrix_get(V,1,0)+gsl_matrix_get(this_gsl,0,2)*gsl_matrix_get(V,2,0);
    A.Y=gsl_matrix_get(this_gsl,0,0)*gsl_matrix_get(V,0,1)+gsl_matrix_get(this_gsl,0,1)*gsl_matrix_get(V,1,1)+gsl_matrix_get(this_gsl,0,2)*gsl_matrix_get(V,2,1);
    A.Z=gsl_matrix_get(this_gsl,0,0)*gsl_matrix_get(V,0,2)+gsl_matrix_get(this_gsl,0,1)*gsl_matrix_get(V,1,2)+gsl_matrix_get(this_gsl,0,2)*gsl_matrix_get(V,2,2);
    B.X=gsl_matrix_get(this_gsl,1,0)*gsl_matrix_get(V,0,0)+gsl_matrix_get(this_gsl,1,1)*gsl_matrix_get(V,1,0)+gsl_matrix_get(this_gsl,1,2)*gsl_matrix_get(V,2,0);
    B.Y=gsl_matrix_get(this_gsl,1,0)*gsl_matrix_get(V,0,1)+gsl_matrix_get(this_gsl,1,1)*gsl_matrix_get(V,1,1)+gsl_matrix_get(this_gsl,1,2)*gsl_matrix_get(V,2,1);
    B.Z=gsl_matrix_get(this_gsl,1,0)*gsl_matrix_get(V,0,2)+gsl_matrix_get(this_gsl,1,1)*gsl_matrix_get(V,1,2)+gsl_matrix_get(this_gsl,1,2)*gsl_matrix_get(V,2,2);
    C.X=gsl_matrix_get(this_gsl,2,0)*gsl_matrix_get(V,0,0)+gsl_matrix_get(this_gsl,2,1)*gsl_matrix_get(V,1,0)+gsl_matrix_get(this_gsl,2,2)*gsl_matrix_get(V,2,0);;
    C.Y=gsl_matrix_get(this_gsl,2,0)*gsl_matrix_get(V,0,1)+gsl_matrix_get(this_gsl,2,1)*gsl_matrix_get(V,1,1)+gsl_matrix_get(this_gsl,2,2)*gsl_matrix_get(V,2,1);;
    C.Z=gsl_matrix_get(this_gsl,2,0)*gsl_matrix_get(V,0,2)+gsl_matrix_get(this_gsl,2,1)*gsl_matrix_get(V,1,2)+gsl_matrix_get(this_gsl,2,2)*gsl_matrix_get(V,2,2);;
    gsl_vector_free(work);
    gsl_vector_free(S);
    gsl_matrix_free (this_gsl);
    gsl_matrix_free (V);
}	
#else
void  XMat3D::Normalise(){std::cerr<<"La matrice n'est pas normalisee"<<std::endl;}
#endif

