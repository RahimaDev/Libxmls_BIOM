// Classe extraite de la lib Matis IGN

#include "libXBase/XQuaternion.h"

//-----------------------------------------------------------------------------
XQuaternion::XQuaternion()
{;}

//-----------------------------------------------------------------------------
XQuaternion::XQuaternion(double x, double y, double z, double w)
{
  m_x = x;
  m_y = y;
  m_z = z;
  m_w = w;
}

//-----------------------------------------------------------------------------
XQuaternion::XQuaternion(XPt3D const &P, double w)
{
  m_x = P.X;
  m_y = P.Y;
  m_z = P.Z;
  m_w = w;
}

//-----------------------------------------------------------------------------
XQuaternion::XQuaternion(XMat3D const &R)
{
  CalculateFromRotationMatrix(R);
}

//-----------------------------------------------------------------------------
XQuaternion::XQuaternion(XQuaternion const &q)
{
  (*this) = q;
}

//-----------------------------------------------------------------------------
//// Constructeur a des angles d'Euler format Helico MAP

// CODE INUTILISe

//-----------------------------------------------------------------------------
XQuaternion::XQuaternion(double roll, double pitch, double yaw)
{
	//Construction d'un quaternion depuis roll/pitch/yaw
	m_x = cos(roll/2.)*cos(pitch/2.)*cos(yaw/2.) + sin(roll/2.)*sin(pitch/2.)*sin(yaw/2.);
	m_y = sin(roll/2.)*cos(pitch/2.)*cos(yaw/2.) - cos(roll/2.)*sin(pitch/2.)*sin(yaw/2.);
	m_z = cos(roll/2.)*sin(pitch/2.)*cos(yaw/2.) + sin(roll/2.)*cos(pitch/2.)*sin(yaw/2.);
	m_w = cos(roll/2.)*cos(pitch/2.)*sin(yaw/2.) - sin(roll/2.)*sin(pitch/2.)*cos(yaw/2.);
	Normalise();
    //double rollQ  = atan2(2*(m_z*m_w + m_x*m_y),m_x*m_x - m_y*m_y - m_z*m_z + m_w*m_w)*180/M_PI;
    //double pitchQ = asin(-2*(m_y*m_w - m_x*m_z)) *180/M_PI;
    //double yawQ   = atan2(2*(m_y*m_z + m_x*m_w),m_x*m_x + m_y*m_y - m_z*m_z - m_w*m_w)*180/M_PI;


	//Construction d'une Matrice depuis roll/pich/yaw
	XPt3D a,b,c;
	a.X = cos(pitch)*cos(yaw);
	a.Y = cos(pitch)*sin(yaw);
	a.Z = -sin(pitch);

	b.X = sin(roll)*sin(pitch)*cos(yaw) - cos(roll)*sin(yaw);
	b.Y = sin(roll)*sin(pitch)*sin(yaw) + cos(roll)*cos(yaw);
	b.Z = sin(roll)*cos(pitch);

	c.X = cos(roll)*sin(pitch)*cos(yaw) + sin(roll)*sin(yaw);
	c.Y = cos(roll)*sin(pitch)*sin(yaw) - sin(roll)*cos(yaw);
	c.Z = cos(roll)*cos(pitch);
	XMat3D DCM(a,b,c);

	//>>> roll//pitch/yaw depuis matrice
    //double rollDCM  = atan2(b.Z, c.Z)*180/M_PI;
    //double pitchDCM = asin(-a.Z) *180/M_PI;
    //double yawDCM   = atan2(a.Y, a.X)*180/M_PI;


	//matrice depuis quaternion
	XPt3D a1,b1,c1;
	a1.X = m_x*m_x + m_y*m_y - m_z*m_z - m_w*m_w;
	a1.Y = 2*(m_y*m_z + m_x*m_w);
	a1.Z = 2*(m_y*m_w - m_x*m_z);

	b1.X = 2*(m_y*m_z - m_x*m_w);
	b1.Y = m_x*m_x - m_y*m_y + m_z*m_z - m_w*m_w;
	b1.Z = 2*(m_z*m_w + m_x*m_y);

	c1.X = 2*(m_y*m_w + m_x*m_z);
	c1.Y = 2*(m_z*m_w - m_x*m_y);
	c1.Z = m_x*m_x - m_y*m_y - m_z*m_z + m_w*m_w;
	XMat3D DCM1(a1,b1,c1);
	//>>roll//pitch/yaw depuis quaternion

/*  Rot.A.X = m_w*m_w + m_x*m_x - m_y*m_y - m_z*m_z;
  Rot.A.Y = 2*(-m_w*m_z + m_x*m_y);
  Rot.A.Z = 2*( m_w*m_y + m_x*m_z);

  Rot.B.X = 2*( m_w*m_z + m_y*m_x);
  Rot.B.Y = m_w*m_w + m_y*m_y - m_x*m_x - m_z*m_z;
  Rot.B.Z = 2*(-m_w*m_x + m_y*m_z);

  Rot.C.X = 2*(-m_w*m_y + m_z*m_x);
  Rot.C.Y = 2*( m_w*m_x + m_z*m_y);
  Rot.C.Z = m_w*m_w + m_z*m_z - m_x*m_x - m_y*m_y;

  Rot /= m_w*m_w + m_x*m_x + m_y*m_y + m_z*m_z;
*/
    //double rollDCM1  = atan2(b1.Z, c1.Z)*180/M_PI;
    //double pitchDCM1 = asin(-a1.Z)*180/M_PI;
    //double yawDCM1   = atan2(a1.Y, a1.X)*180/M_PI;

}
//-----------------------------------------------------------------------------
double XQuaternion::Norme() const
{
  return sqrt((m_w*m_w) + (m_x*m_x) + (m_y*m_y) + (m_z*m_z));
}

//-----------------------------------------------------------------------------
void XQuaternion::Normalise()
{
  double n=Norme();
  if (n != 0)
    {
      m_w /= n;
      m_x /= n;
      m_y /= n;
      m_z /= n;
    }
}

//-----------------------------------------------------------------------------
XQuaternion XQuaternion::Conjugue()
{
  return XQuaternion(-m_x,-m_y,-m_z,m_w);
}

//-----------------------------------------------------------------------------
XQuaternion XQuaternion::operator* (const XQuaternion & r) const
{
  XQuaternion res;
  res.m_w = m_w*r.m_w - m_x*r.m_x - m_y*r.m_y - m_z*r.m_z;
  res.m_x = m_w*r.m_x + m_x*r.m_w + m_y*r.m_z - m_z*r.m_y;
  res.m_y = m_w*r.m_y + m_y*r.m_w + m_z*r.m_x - m_x*r.m_z;
  res.m_z = m_w*r.m_z + m_z*r.m_w + m_x*r.m_y - m_y*r.m_x;
  return res;
}

//-----------------------------------------------------------------------------
XQuaternion XQuaternion::operator+ (const XQuaternion & q) const
{
  XQuaternion res;
  res.m_w = m_w + q.m_w;
  res.m_x = m_x + q.m_x;
  res.m_y = m_y + q.m_y;
  res.m_z = m_z + q.m_z;
  return res;
}

//-----------------------------------------------------------------------------
XQuaternion XQuaternion::operator- (const XQuaternion & q) const
{
  XQuaternion res;
  res.m_w = m_w - q.m_w;
  res.m_x = m_x - q.m_x;
  res.m_y = m_y - q.m_y;
  res.m_z = m_z - q.m_z;
  return res;
}

//-----------------------------------------------------------------------------
void XQuaternion::operator+= (const XQuaternion &q)
{
  m_w += q.m_w;
  m_x += q.m_x;
  m_y += q.m_y;
  m_z += q.m_z;
}

//-----------------------------------------------------------------------------
void XQuaternion::operator-= (const XQuaternion &q)
{
  m_w -= q.m_w;
  m_x -= q.m_x;
  m_y -= q.m_y;
  m_z -= q.m_z;
}

//-----------------------------------------------------------------------------
XQuaternion & XQuaternion::operator=(const XQuaternion &q)
{
  m_w = q.m_w;
  m_x = q.m_x;
  m_y = q.m_y;
  m_z = q.m_z;
  return (*this);
}

//-----------------------------------------------------------------------------
XMat3D XQuaternion::GetRotationMatrix() const
{
  XMat3D rot;
  GetRotationMatrix(rot);
  return rot;
}

//-----------------------------------------------------------------------------
void XQuaternion::GetRotationMatrix(XMat3D &Rot) const
{
  Rot.A.X = m_w*m_w + m_x*m_x - m_y*m_y - m_z*m_z;
  Rot.A.Y = 2*(-m_w*m_z + m_x*m_y);
  Rot.A.Z = 2*( m_w*m_y + m_x*m_z);

  Rot.B.X = 2*( m_w*m_z + m_y*m_x);
  Rot.B.Y = m_w*m_w + m_y*m_y - m_x*m_x - m_z*m_z;
  Rot.B.Z = 2*(-m_w*m_x + m_y*m_z);

  Rot.C.X = 2*(-m_w*m_y + m_z*m_x);
  Rot.C.Y = 2*( m_w*m_x + m_z*m_y);
  Rot.C.Z = m_w*m_w + m_z*m_z - m_x*m_x - m_y*m_y;

  Rot /= m_w*m_w + m_x*m_x + m_y*m_y + m_z*m_z;
}

//-----------------------------------------------------------------------------
XMat3D XQuaternion::GetdR_dw() const
{
  XMat3D dR;
  dR.A.X = 2*m_w;
  dR.A.Y =-2*m_z;
  dR.A.Z = 2*m_y;

  dR.B.X = 2*m_z;
  dR.B.Y = 2*m_w;
  dR.B.Z =-2*m_x;

  dR.C.X =-2*m_y;
  dR.C.Y = 2*m_x;
  dR.C.Z = 2*m_w;

  XMat3D Rot;
  GetRotationMatrix(Rot);
  dR -= Rot * (2 * m_w) ;

  dR /= m_w*m_w + m_x*m_x + m_y*m_y + m_z*m_z;

  return dR;
}

//-----------------------------------------------------------------------------
XMat3D XQuaternion::GetdR_dx() const
{
  XMat3D dR;
  dR.A.X = 2*m_x;
  dR.A.Y = 2*m_y;
  dR.A.Z = 2*m_z;

  dR.B.X = 2*m_y;
  dR.B.Y =-2*m_x;
  dR.B.Z =-2*m_w;

  dR.C.X = 2*m_z;
  dR.C.Y = 2*m_w;
  dR.C.Z =-2*m_x;

  XMat3D Rot;
  GetRotationMatrix(Rot);
  dR -= Rot * (2 * m_x) ;

  dR /= m_w*m_w + m_x*m_x + m_y*m_y + m_z*m_z;

  return dR;
}

//-----------------------------------------------------------------------------
XMat3D XQuaternion::GetdR_dy() const
{
  XMat3D dR;
  dR.A.X =-2*m_y;
  dR.A.Y = 2*m_x;
  dR.A.Z = 2*m_w;

  dR.B.X = 2*m_x;
  dR.B.Y = 2*m_y;
  dR.B.Z = 2*m_z;

  dR.C.X =-2*m_w;
  dR.C.Y = 2*m_z;
  dR.C.Z =-2*m_y;

  XMat3D Rot;
  GetRotationMatrix(Rot);
  dR -= Rot * (2 * m_y) ;

  dR /= m_w*m_w + m_x*m_x + m_y*m_y + m_z*m_z;

  return dR;
}

//-----------------------------------------------------------------------------
XMat3D XQuaternion::GetdR_dz() const
{
  XMat3D dR;
  dR.A.X =-2*m_z;
  dR.A.Y =-2*m_w;
  dR.A.Z = 2*m_x;

  dR.B.X = 2*m_w;
  dR.B.Y =-2*m_z;
  dR.B.Z = 2*m_y;

  dR.C.X = 2*m_x;
  dR.C.Y = 2*m_y;
  dR.C.Z = 2*m_z;

  XMat3D Rot;
  GetRotationMatrix(Rot);
  dR -= Rot * (2 * m_z) ;

  dR /= m_w*m_w + m_x*m_x + m_y*m_y + m_z*m_z;

  return dR;
}

//-----------------------------------------------------------------------------
bool XQuaternion::CalculateFromRotationMatrix(XMat3D const &R)
{
  int max = -1;
  double val, valmax = 0;
  val = 0.25 * (1 + R.A.X + R.B.Y + R.C.Z);
  if (val > valmax)
    {
      max = 0;
      valmax = val;
    }
  val = 0.25 * (1 + R.A.X - R.B.Y - R.C.Z);
  if (val > valmax)
    {
      max = 1;
      valmax = val;
    }
  val = 0.25 * (1 - R.A.X + R.B.Y - R.C.Z);
  if (val > valmax)
    {
      max = 2;
      valmax = val;
    }
  val = 0.25 * (1 - R.A.X - R.B.Y + R.C.Z);
  if (val > valmax)
    {
      max = 3;
      valmax = val;
    }
  switch (max)
    {
    case 0 :
      {
	m_w = sqrt(valmax);
	m_x = 0.25 * (R.C.Y - R.B.Z) / m_w;
	m_y = 0.25 * (R.A.Z - R.C.X) / m_w;
	m_z = 0.25 * (R.B.X - R.A.Y) / m_w;
	break;
      }
    case 1 :
      {
	m_x = sqrt(valmax);
	m_w = 0.25 * (R.C.Y - R.B.Z) / m_x;
	if (m_w < 0)
	  {
	    m_w = -m_w;
	    m_x = -m_x;
	  }
	m_y = 0.25 * (R.A.Y + R.B.X) / m_x;
	m_z = 0.25 * (R.A.Z + R.C.X) / m_x;
	break;
      }
    case 2 :
      {
	m_y = sqrt(valmax);
	m_w = 0.25 * (R.A.Z - R.C.X) / m_y;
	if (m_w < 0)
	  {
	    m_w = -m_w;
	    m_y = -m_y;
	  }
	m_x = 0.25 * (R.A.Y + R.B.X) / m_y;
	m_z = 0.25 * (R.B.Z + R.C.Y) / m_y;
	break;
      }
    case 3 :
      {
	m_z = sqrt(valmax);
	m_w = 0.25 * (R.B.X - R.A.Y) / m_z;
	if (m_w < 0)
	  {
	    m_w = -m_w;
	    m_z = -m_z;
	  }
	m_x = 0.25 * (R.A.Z + R.C.X) / m_z;
	m_y = 0.25 * (R.B.Z + R.C.Y) / m_z;
	break;
      }
    default :
    return false;
//      THROWMESS(Erreur, "Quaternion:CalculateFromRotationMatrix() : Cas impossible...\n");
    };
    return true;
}

//-----------------------------------------------------------------------------
double XQuaternion::GetAngle() const
{
  return 2.*atan2(sqrt((m_x*m_x)+(m_y*m_y)+(m_z*m_z)),m_w);
}

//-----------------------------------------------------------------------------
XPt3D XQuaternion::GetAxis() const
{
  double teta=GetAngle();
  if (teta==0.0)
    return XPt3D(1,0,0);
  XPt3D axe(m_x/sin(teta/2.),m_y/sin(teta/2.),m_z/sin(teta/2.));
  double norme = sqrt((axe.X*axe.X) + (axe.Y*axe.Y) +(axe.Z*axe.Z));
  return axe/norme;
 // return axe;
}

/*Nostream& operator << (Nostream& sortie, Quaternion const & p)
{
  sortie <<"[(" << p.m_x << "," << p.m_y << "," << p.m_z << "),"<< p.m_w<<"]";
  return sortie;
}*/


//-----------------------------------------------------------------------------
// Ecriture dans un fichier Texte
//-----------------------------------------------------------------------------
bool XQuaternion::WriteTxt(std::ostream* out)
{
	*out <<  "q0: "<< q0() << "  q1: "<< q1()  << "  q2: " << q2() << "  q3: " << q3() <<std::endl;
	XPt3D axe = GetAxis();
	*out <<  "axe:  x = "<< axe.X << "  y = "<< axe.Y  << "  z = " << axe.Z  <<std::endl;
	*out << "angle = "<< GetAngle() << std::endl;
	*out << std::endl;
	return out->good();
}

//-----------------------------------------------------------------------------
XPt3D XQuaternion::OmegaPhiKappaBis()
{
	Normalise();

	//resultats tres bizarres avec cet algo
	//permutation pour test a partir de matrice apx
	//permutation des elements
	double w = m_x;
	double x = m_y;
	double y = m_z;
	double z = m_w;
	//fin permutation

    double sqw = w*w;
	double sqx = x*x;
	double sqy = y*y;
	double sqz = z*z;
	double unit = sqx + sqy + sqz + sqw; // if normalised is one, otherwise is correction factor
	double test = x*y + z*w;

	double attitude,heading,bank;

	if (test > 0.499*unit) { // singularity at north pole
		heading = 2 * atan2(x,w);
		attitude = M_PI/2;
		bank = 0;
		return XPt3D(attitude,heading,bank);
	}
	if (test < -0.499*unit) { // singularity at south pole
		heading = -2 * atan2(x,w);
		attitude = -M_PI/2;
		bank = 0;
		return XPt3D(attitude,heading,bank);
	}
	heading = atan2(2*y*w-2*x*z , sqx - sqy - sqz + sqw);
	attitude = asin(2*test/unit);
	bank = atan2(2*x*w-2*y*z , -sqx + sqy - sqz + sqw);
	return XPt3D(attitude,heading,bank);

}
//-----------------------------------------------------------------------------
XPt3D XQuaternion::OmegaPhiKappaTopo()
{
	double rollQ  = atan2(2*(m_z*m_w + m_x*m_y),m_x*m_x - m_y*m_y - m_z*m_z + m_w*m_w);
	double pitchQ = asin(-2*(m_y*m_w - m_x*m_z));
	double yawQ   = atan2(2*(m_y*m_z + m_x*m_w),m_x*m_x + m_y*m_y - m_z*m_z - m_w*m_w);

	return XPt3D(yawQ, pitchQ, 2*M_PI-rollQ);
}
//-----------------------------------------------------------------------------
XPt3D XQuaternion::OmegaPhiKappa()
{

	double rollQ  = atan2(2*(m_z*m_w + m_x*m_y),m_x*m_x - m_y*m_y - m_z*m_z + m_w*m_w);
	double pitchQ = asin(-2*(m_y*m_w - m_x*m_z));
	double yawQ   = atan2(2*(m_y*m_z + m_x*m_w),m_x*m_x + m_y*m_y - m_z*m_z - m_w*m_w);
//	return XPt3D(rollQ,pitchQ,yawQ);
	return XPt3D(M_PI-yawQ,-pitchQ,rollQ);
}

//-----------------------------------------------------------------------------
XPt3D XQuaternion::OmegaPhiKappaTer()
{
//************ CODE INUTILISe  en test

	double rollQ  = atan2(2*(m_x*m_y + m_z*m_w), 1-2*(m_y*m_y + m_z*m_z));
	double pitchQ = asin(2*(m_x*m_z - m_w*m_y));
	double yawQ   = atan2(2*(m_x*m_w + m_y*m_z), 1-2*(m_z*m_z + m_w*m_w));
	return XPt3D(rollQ,pitchQ,yawQ);
}
//-----------------------------------------------------------------------------
XPt3D XQuaternion::OmegaPhiKappaStar()
{
	/* Recuperation de l'angle de rotation */
	double angle = acos(m_w) * 2;

	/* Recuperation des composantes de l'axe de rotation */
	double vx = m_x;
	double vy = m_y;
	double vz = m_z;

	/* Normalisation de l'axe de rotation */
	double norm = sqrt(vx * vx + vy * vy + vz * vz);
	if (norm > 0.0005)
	{
		vx /= norm;
		vy /= norm;
		vz /= norm;
	}

	/* Calcul de la latitude */
	double latitude = -asin(vy);

	double longitude;
	/* Calcul de la longitude */
	if (vx * vx + vz * vz < 0.0005)
		longitude = 0;
	else
		longitude = atan2(vx, vz);

	/* Si la longitude est negative, on la ramene du côte positif */
	if (longitude < 0)
		longitude += 2 * M_PI;

	return XPt3D(latitude,longitude,angle);
}
//-----------------------------------------------------------------------------
// Ecriture dans un fichier XML
//-----------------------------------------------------------------------------
bool XQuaternion::XmlWrite(std::ostream* out)
{
	*out << "<quaternion> " << std::endl;
	*out << "<x> " << m_x << " </x>" << std::endl;
	*out << "<y> " << m_y << " </y>" << std::endl;
	*out << "<z> " << m_z << " </z>" << std::endl;
	*out << "<w> " << m_w << " </w>" << std::endl;
	*out << "</quaternion>" << std::endl;
	return out->good();
}
