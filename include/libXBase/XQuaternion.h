//XQuaternion
//Convertion adaptation de la lib Matis

#ifndef __XQUATERNION_H__
#define __XQUATERNION_H__

#include "libXBase/XPt3D.h"
#include "libXBase/XMat3D.h"

//// Classe de quaternions, tres utiles pour representer les rotations
class XQuaternion
{
	public :
		XQuaternion();

		//// Constructeur a partir de 4 valeurs : x,y,z,w ou q1,q2,q3,q0
		XQuaternion(double x, double y, double z, double w);

		//// Constructeur a partir de 4 valeurs : x,y,z,w ou q1,q2,q3,q0
		XQuaternion(XPt3D const &P, double w = 0);

		//// Constructeur par recopie
		XQuaternion(XQuaternion const &q);

		//// Constructeur a partir d'une matrice rotation
		XQuaternion(XMat3D const &Rot);

		//// Constructeur a des angles d'Euler format Helico MAP
		XQuaternion(double roll, double pitch, double yaw);

		bool XmlWrite(std::ostream* out);

		double Norme() const;
		void Normalise();

		//// Operateurs mathematiques
		XQuaternion Conjugue();
		XQuaternion operator* (const XQuaternion & q) const;
		XQuaternion operator+ (const XQuaternion & q) const;
		XQuaternion operator- (const XQuaternion & q) const;
		void operator += (const XQuaternion &q);
		void operator -= (const XQuaternion &q);

		XQuaternion & operator=(const XQuaternion &q);

		//// Constructeur a partir d'une matrice rotation
		bool CalculateFromRotationMatrix(XMat3D const &R);

		//// Calcul de la matrice rotation a partir du quaternion
		void GetRotationMatrix(XMat3D &Rot) const;
		//// Calcul de la matrice rotation a partir du quaternion
		XMat3D GetRotationMatrix() const;

		//// Calcul de la derivee de la matrice rotation pa rapport a w (q0)
		XMat3D GetdR_dw() const;
		//// Calcul de la derivee de la matrice rotation pa rapport a x (q1)
		XMat3D GetdR_dx() const;
		//// Calcul de la derivee de la matrice rotation pa rapport a y (q2)
		XMat3D GetdR_dy() const;
		//// Calcul de la derivee de la matrice rotation pa rapport a z (q3)
		XMat3D GetdR_dz() const;

		inline double q0() const {return m_w;}
		inline double q1() const {return m_x;}
		inline double q2() const {return m_y;}
		inline double q3() const {return m_z;}

		inline double w() const {return m_w;}
		inline double x() const {return m_x;}
		inline double y() const {return m_y;}
		inline double z() const {return m_z;}

		//// Calcul de l'angle de la rotation a partir du quaternion
		double GetAngle() const;

		//// Calcul de l'axe de la rotation a partir du quaternion
		XPt3D GetAxis() const;

		virtual bool WriteTxt(std::ostream* out);

		XPt3D OmegaPhiKappaTopo();
		XPt3D OmegaPhiKappa();
		XPt3D OmegaPhiKappaBis();
		XPt3D OmegaPhiKappaTer();
		XPt3D OmegaPhiKappaStar();

		void w(double val) {m_w = val;}
		void x(double val) {m_x = val;}
		void y(double val) {m_y = val;}
		void z(double val) {m_z = val;}

	protected :
		double m_w,m_x,m_y,m_z;
};

#endif
