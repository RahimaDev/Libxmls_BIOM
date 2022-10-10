#ifndef _XCOLORTOOLS_H
#define _XCOLORTOOLS_H

#include "XPt3D.h"

class XColorTools 
{
protected:
		
	float Magic(float rm1, float rm2, float rh)
	{
	  float retval = rm1;
	  if (rh > 360.0)
		rh -= 360.0;
	  if (rh < 0.0)
		rh += 360.0;
	  if (rh < 60.0)
		retval = (float)(rm1 + (rm2 - rm1) * rh / 60.0);
	  else if (rh < 180.0)
		retval = rm2;
	  else if (rh < 240.0)
		retval = (float)(rm1 + (rm2 - rm1) * (240.0 - rh) / 60.0);
	  return (float)ceil (retval * 255);
	}

public:
	// Constructeur
	XColorTools() {;}

	std::string HSL2HexaKML(float h,float s,float l, byte transparency = 255)
	{
		XPt3D RGB = HSL2RGB(h,s,l);
		char str[1024];
		sprintf(str,"%02X%02X%02X%02X",transparency,(byte)RGB.Z,(byte)RGB.Y,(byte)RGB.X);
		return std::string(str);
	}
	static XPt3D RGBtoHSL(XPt3D colorRGB)
	{
	    float r, g, b, h, s, l; //this function works with floats between 0 and 1
	    r = (float)(colorRGB.X / 256.0);
		g = (float)(colorRGB.Y / 256.0);
		b = (float)(colorRGB.Z / 256.0);

	 //Then, minColor and maxColor are defined. Mincolor is the value of the color component with the smallest value,
	 //while maxColor is the value of the color component with the largest value.
	 //These two variables are needed because the Lightness is defined as (minColor + maxColor) / 2.
	    float maxColor = XMax(r, XMax(g, b));
	    float minColor = XMin(r, XMin(g, b));

	//If minColor equals maxColor, we know that R=G=B and thus the color is a shade of gray.
	//This is a trivial case, hue can be set to anything, saturation has to be set to 0 because only then it's a shade of gray, and lightness is set to R=G=B, the shade of the gray.
	//R == G == B, so it's a shade of gray
	    if(minColor == maxColor)
	    {
	        h = 0.0; //it doesn't matter what value it has
	        s = 0.0;
	        l = r; //doesn't matter if you pick r, g, or b
	    }

	//If minColor is not equal to maxColor, we have a real color instead of a shade of gray, so more calculations are needed:
	//Lightness (l) is now set to it's definition of (minColor + maxColor)/2.
	//Saturation (s) is then calculated with a different formula depending if light is in the first half of the second half. This is because the HSL model can be represented as a double cone, the first cone has a black tip and corresponds to the first half of lightness values, the second cone has a white tip and contains the second half of lightness values.
	//Hue (h) is calculated with a different formula depending on which of the 3 color components is the dominating one, and then normalized to a number between 0 and 1.

	    else
	    {
	        l = (minColor + maxColor) / 2;

	        if(l < 0.5) 
				s = (maxColor - minColor) / (maxColor + minColor);
	        else 
				s = (float)((maxColor - minColor) / (2.0 - maxColor - minColor));

	        if(r == maxColor) 
				h = (g - b) / (maxColor - minColor);
			else if(g == maxColor) 
				h = (float)(2.0 + (b - r) / (maxColor - minColor));
			else h = (float)(4.0 + (r - g) / (maxColor - minColor));

	        h /= 6; //to bring it to a number between 0 and 1
	        if(h < 0) 
				h ++;
	    }

	//Finally, H, S and L are calculated out of h,s and l as integers between 0 and 255 and "returned" as the result.
	//Returned, because H, S and L were passed by reference to the function.

	    XPt3D colorHSL;
	    colorHSL.X = int(h * 255.0);
	    colorHSL.Y = int(s * 255.0);
	    colorHSL.Z = int(l * 255.0);

    	return colorHSL;
	}



	static XPt3D HSLToRGB(XPt3D colorHSL)
	{
    float r, g, b, h, s, l; //this function works with floats between 0 and 1
    float temp1, temp2, tempr, tempg, tempb;
	h = (float)(colorHSL.X / 256.0);
	s = (float)(colorHSL.Y / 256.0);
	l = (float)(colorHSL.Z / 256.0);

    //Then follows a trivial case: if the saturation is 0,
    //the color will be a grayscale color, and the calculation is then very simple: r, g and b are all set to the lightness.

    //If saturation is 0, the color is a shade of gray
    if(s == 0) r = g = b = l;

	//If the saturation is higher than 0, more calculations are needed again. red, green and blue are calculated
	//with the formulas defined in the code.

    //If saturation > 0, more complex calculations are needed
    else
    {
        //Set the temporary values
        if(l < 0.5) temp2 = l * (1 + s);
        else temp2 = (l + s) - (l * s);
        temp1 = 2 * l - temp2;
		tempr = (float)(h + 1.0 / 3.0);
        if(tempr > 1) tempr--;
        tempg = h;
		tempb =(float)( h - 1.0 / 3.0);
        if(tempb < 0) tempb++;

        //Red
		if(tempr < 1.0 / 6.0) r = (float)(temp1 + (temp2 - temp1) * 6.0 * tempr);
        else if(tempr < 0.5) r = temp2;
		else if(tempr < 2.0 / 3.0) r = (float)(temp1 + (temp2 - temp1) * ((2.0 / 3.0) - tempr) * 6.0);
        else r = temp1;

        //Green
		if(tempg < 1.0 / 6.0) g = (float)(temp1 + (temp2 - temp1) * 6.0 * tempg);
        else if(tempg < 0.5) g = temp2;
		else if(tempg < 2.0 / 3.0) g = (float)(temp1 + (temp2 - temp1) * ((2.0 / 3.0) - tempg) * 6.0);
        else g = temp1;

        //Blue
		if(tempb < 1.0 / 6.0) b = (float)(temp1 + (temp2 - temp1) * 6.0 * tempb);
        else if(tempb < 0.5) b = temp2;
		else if(tempb < 2.0 / 3.0) b = (float)(temp1 + (temp2 - temp1) * ((2.0 / 3.0) - tempb) * 6.0);
        else b = temp1;
    }

	//And finally, the results are returned as integers between 0 and 255.

    XPt3D  colorRGB;
    colorRGB.X = int(r * 255.0);
    colorRGB.Y = int(g * 255.0);
    colorRGB.Z = int(b * 255.0);
    return colorRGB;
    }



	XPt3D HSL2RGB (float h,float s,float l)
	{
	  float hue = (float)(h * 360.0 / 255.0);
	  float saturation = (float)(s / 255.0);
	  float luminance = (float)(l / 255.0);
	  float red;
	  float green;
	  float blue;
	  if (saturation == 0.0)
	  {
		  red = green = blue = (float)(ceil (luminance * 255.0));
	  }
	  else
	  {
		float rm1;
		float rm2;
		if (luminance <= 0.5)
		  rm2 = luminance + luminance * saturation;
		else
		  rm2 = luminance + saturation - luminance * saturation;
		rm1 = (float)(2.0 * luminance - rm2);
		red = Magic (rm1, rm2, (float)(hue + 120.0));
		green = Magic (rm1, rm2, hue);
		blue = Magic (rm1, rm2, (float)(hue - 120.0));
	  }
	  return XPt3D(red, green, blue);
	}
};
template <class T>
inline T Borne(float n)
{
  if(n>T(0xffffffff))
    return T(0xffffffff);
  if(n<0)
    return 0;
  else
    return T(n);
}
template<class T>
inline void RVB_ITS2(const T r,const T v,const T b,T &i,T &t,T &s){

   double di, dt, ds;

   double dr = double(r);
   double dv = double(v);
   double db = double(b);

   // Normalisation [0;255] -> [0;1]
   dr /= T(0xffffffff);
   dv /= T(0xffffffff);
   db /= T(0xffffffff);

   // Calcul de l'intensite
   di = (dr+dv+db) / 3.;


   // Si r=v=b la couleur est une nuance de gris :
   // la saturation est nulle
   // la teinte n'est pas definie (on la prend egale a 0)
   // cette convention conforme aux formules de retranscription
   // ITS->RVB n'empeche pas la resversibilite de la transfo
   if(r==v && r==b){
      ds = 0;
      dt = 0;
   }
   else{
      // Calcul de la saturation
      double min = dr;
      if(dv < min)
         min = dv;
      if(db < min)
         min = db;

      ds = 1 - min/di;  // di=0 <=> R=V=B : cas ecarte !

      // Calcul de la teinte
      dt = acos(              ((dr-dv)+(dr-db)) 
                 / (2*sqrt((dr-dv)*(dr-dv)+(dr-db)*(dv-db)))  );

      if (db > dv)
         dt = 2*M_PI - dt;

      dt /= 2*M_PI;       // Normalisation -> [0;1]
   }

   // Retranscription dans le type template T
   // et passage [0;1]->[0;255]
   i = T(di * T(0xffffffff));
   t = T(dt * T(0xffffffff));
   s = T(ds * T(0xffffffff));

}

template<class T>
inline void ITS_RVB2(const T i,const T t,const T s,T &r,T &v,T &b){

   double dr=0.0,dv=0.0,db=0.0;

   double di = double(i);

   // Normalisation [0;255] -> [0;1]
   di /= T(0xffffffff);

   double dt = double(t);
   double ds = double(s);

   // Normalisation [0;255] -> [0;1]
   dt /= T(0xffffffff);
   ds /= T(0xffffffff);

   // retraduction de la teinte en terme d'angle
   dt *= 2*M_PI;

   if( dt>=0. && dt<=2*M_PI/3. ){
      db = (1-ds) / 3.;
      dr = ( 1 + ( ds*cos(dt) / cos(M_PI/3.-dt) ) ) / 3.;
      dv = 1. - (db+dr);
   }

   if( dt>2*M_PI/3. && dt<=4*M_PI/3. ){
      dt -= 2*M_PI/3.;
      dr = (1-ds) / 3.;
      dv = ( 1 + ( ds*cos(dt) / cos(M_PI/3.-dt) ) ) / 3.;
      db = 1. - (dr+dv);
   }

   if( dt>4*M_PI/3. && dt<=2*M_PI ){
      dt -= 4*M_PI/3.;
      dv = (1-ds) / 3.;
      db = ( 1 + ( ds*cos(dt) / cos(M_PI/3.-dt)) ) / 3.;
      dr = 1. - (dv+db);
   }

   // jusque la, dr dv et db tels que dr+dv+db=1
   // denormalisation sachant que I=(R+G+B)/3 
   // et passage [0;1]->[0;255]

   dr *= 3*di*T(0xffffffff);
   dv *= 3*di*T(0xffffffff);
   db *= 3*di*T(0xffffffff);

   // Retranscription dans le type template T

 /*  r = T(dr);
   v = T(dv);
   b = T(db);*/
   r = Borne<T>(dr);
   v = Borne<T>(dv);
   b = Borne<T>(db);

}

#endif //_XCOLORTOOLS_H
/*

ColorHSL RGBtoHSL(ColorRGB colorRGB)
{
    float r, g, b, h, s, l; //this function works with floats between 0 and 1
    r = colorRGB.r / 256.0;
    g = colorRGB.g / 256.0;
    b = colorRGB.b / 256.0;



Then, minColor and maxColor are defined. Mincolor is the value of the color component with the smallest value, while maxColor is the value of the color component with the largest value. These two variables are needed because the Lightness is defined as (minColor + maxColor) / 2.


    float maxColor = max(r, max(g, b));
    float minColor = min(r, min(g, b));



If minColor equals maxColor, we know that R=G=B and thus the color is a shade of gray. This is a trivial case, hue can be set to anything, saturation has to be set to 0 because only then it's a shade of gray, and lightness is set to R=G=B, the shade of the gray.


    //R == G == B, so it's a shade of gray
    {
        h = 0.0; //it doesn't matter what value it has
        s = 0.0;
        l = r; //doesn't matter if you pick r, g, or b
    }



If minColor is not equal to maxColor, we have a real color instead of a shade of gray, so more calculations are needed:

Lightness (l) is now set to it's definition of (minColor + maxColor)/2.
Saturation (s) is then calculated with a different formula depending if light is in the first half of the second half. This is because the HSL model can be represented as a double cone, the first cone has a black tip and corresponds to the first half of lightness values, the second cone has a white tip and contains the second half of lightness values.
Hue (h) is calculated with a different formula depending on which of the 3 color components is the dominating one, and then normalized to a number between 0 and 1.



    else
    {
        l = (minColor + maxColor) / 2;

        if(l < 0.5) s = (maxColor - minColor) / (maxColor + minColor);
        else s = (maxColor - minColor) / (2.0 - maxColor - minColor);

        if(r == maxColor) h = (g - b) / (maxColor - minColor);
        else if(g == maxColor) h = 2.0 + (b - r) / (maxColor - minColor);
        else h = 4.0 + (r - g) / (maxColor - minColor);

        h /= 6; //to bring it to a number between 0 and 1
        if(h < 0) h ++;
    }



Finally, H, S and L are calculated out of h,s and l as integers between 0 and 255 and "returned" as the result. Returned, because H, S and L were passed by reference to the function.


    ColorHSL colorHSL;
    colorHSL.h = int(h * 255.0);
    colorHSL.s = int(s * 255.0);
    colorHSL.l = int(l * 255.0);
    return colorHSL;
}




HSL to RGB

This is the opposite conversion, so this function will calculate the inverse of the RGBtoHSL function.

First, internally the variables with small letters are defined as floating point numbers between 0 and 1 again. Some temporary values for the calculations are also declared.


ColorRGB HSLtoRGB(ColorHSL colorHSL)
{
    float r, g, b, h, s, l; //this function works with floats between 0 and 1
    float temp1, temp2, tempr, tempg, tempb;
    h = colorHSL.h / 256.0;
    s = colorHSL.s / 256.0;
    l = colorHSL.l / 256.0;



Then follows a trivial case: if the saturation is 0, the color will be a grayscale color, and the calculation is then very simple: r, g and b are all set to the lightness.


    //If saturation is 0, the color is a shade of gray
    if(s == 0) r = g = b = l;



If the saturation is higher than 0, more calculations are needed again. red, green and blue are calculated with the formulas defined in the code.


    //If saturation > 0, more complex calculations are needed
    else
    {
        //Set the temporary values
        if(l < 0.5) temp2 = l * (1 + s);
        else temp2 = (l + s) - (l * s);
        temp1 = 2 * l - temp2;
        tempr = h + 1.0 / 3.0;
        if(tempr > 1) tempr--;
        tempg = h;
        tempb = h - 1.0 / 3.0;
        if(tempb < 0) tempb++;

        //Red
        if(tempr < 1.0 / 6.0) r = temp1 + (temp2 - temp1) * 6.0 * tempr;
        else if(tempr < 0.5) r = temp2;
        else if(tempr < 2.0 / 3.0) r = temp1 + (temp2 - temp1) * ((2.0 / 3.0) - tempr) * 6.0;
        else r = temp1;

        //Green
        if(tempg < 1.0 / 6.0) g = temp1 + (temp2 - temp1) * 6.0 * tempg;
        else if(tempg < 0.5) g = temp2;
        else if(tempg < 2.0 / 3.0) g = temp1 + (temp2 - temp1) * ((2.0 / 3.0) - tempg) * 6.0;
        else g = temp1;

        //Blue
        if(tempb < 1.0 / 6.0) b = temp1 + (temp2 - temp1) * 6.0 * tempb;
        else if(tempb < 0.5) b = temp2;
        else if(tempb < 2.0 / 3.0) b = temp1 + (temp2 - temp1) * ((2.0 / 3.0) - tempb) * 6.0;
        else b = temp1;
    }



And finally, the results are returned as integers between 0 and 255.


    ColorRGB colorRGB;
    colorRGB.r = int(r * 255.0);
    colorRGB.g = int(g * 255.0);
    colorRGB.b = int(b * 255.0);
    return colorRGB;
}




RGB to HSV
The function RGBtoHSV works very similar as the RGBtoHSL function, the only difference is that now the variable V (Value) instead of L (Lightness) is used, and Value is defined as maxColor. This can immediately be calculated at the beginning of the function:


ColorHSV RGBtoHSV(ColorRGB colorRGB)
{
    float r, g, b, h, s, v; //this function works with floats between 0 and 1
    r = colorRGB.r / 256.0;
    g = colorRGB.g / 256.0;
    b = colorRGB.b / 256.0;
    float maxColor = max(r, max(g, b));
    float minColor = min(r, min(g, b));
    v = maxColor;



Then, the saturation is calculated. If the color is black, the value of saturation doesn't matter so it can be set to 0. This has to be done to avoid a division by zero.


    if(maxColor == 0)//avoid division by zero when the color is black
    {
        s = 0;
    }
    else
    {
        s = (maxColor - minColor) / maxColor;
    }



Finally, the hue is calculated. If saturation is 0, the color is gray so hue doesn't matter. Again this case is handled separately to avoid divisions by zero.


    if(s == 0)
    {
        h = 0; //it doesn't matter what value it has
    }
    else
    {
        if(r == maxColor) h = (g - b) / (maxColor-minColor);
        else if(g == maxColor) h = 2.0 + (b - r) / (maxColor - minColor);
        else h = 4.0 + (r - g) / (maxColor - minColor);
        h /= 6.0; //to bring it to a number between 0 and 1
        if (h < 0) h++;
    }



And finally, the results are returned as integers between 0 and 255.


    ColorHSV colorHSV;
    colorHSV.h = int(h * 255.0);
    colorHSV.s = int(s * 255.0);
    colorHSV.v = int(v * 255.0);
    return colorHSV;
}




HSV to RGB

First the floating point numbers between 0 and 1 are declared again:


ColorRGB HSVtoRGB(ColorHSV colorHSV)
{
    float r, g, b, h, s, v; //this function works with floats between 0 and 1
    h = colorHSV.h / 256.0;
    s = colorHSV.s / 256.0;
    v = colorHSV.v / 256.0;



The trivial case for saturation = zero is handled:


    //If saturation is 0, the color is a shade of gray
    if(s == 0) r = g = b = v;



The HSV model can be presented on a cone with hexagonal shape. For each of the sides of the hexagon, a separate case is calculated:


    //If saturation > 0, more complex calculations are needed
    else
    {
        float f, p, q, t;
        int i;
        h *= 6; //to bring hue to a number between 0 and 6, better for the calculations
        i = int(floor(h));  //e.g. 2.7 becomes 2 and 3.01 becomes 3 or 4.9999 becomes 4
        f = h - i;  //the fractional part of h
        p = v * (1 - s);
        q = v * (1 - (s * f));
        t = v * (1 - (s * (1 - f)));
        switch(i)
        {
            case 0: r = v; g = t; b = p; break;
            case 1: r = q; g = v; b = p; break;
            case 2: r = p; g = v; b = t; break;
            case 3: r = p; g = q; b = v; break;
            case 4: r = t; g = p; b = v; break;
            case 5: r = v; g = p; b = q; break;
        }
    }



And again, the results are "returned" as integers between 0 and 255.


    ColorRGB colorRGB;
    colorRGB.r = int(r * 255.0);
    colorRGB.g = int(g * 255.0);
    colorRGB.b = int(b * 255.0);
    return colorRGB;
}



*/

/*
/ JScript source code
//Red : 0..255
//Green : 0..255
//Blue : 0..255
//Hue : 0,0..360,0<=>0..255
//Lum : 0,0..1,0<=>0..255
//Sat : 0,0..1,0<=>0..255

//Retourne un tableau de 3 valeurs : H,S,L
function RGB2HSL (r, g, b)
{
  red = Math.round (r);
  green = Math.round (g);
  blue = Math.round (b);
  var minval = Math.min (red, Math.min (green, blue));
  var maxval = Math.max (red, Math.max (green, blue));
  var mdiff = maxval - minval + 0.0;
  var msum = maxval + minval + 0.0;
  var luminance = msum / 510.0;
  var saturation;
  var hue;
  if (maxval == minval)
  {
    saturation = 0.0;
    hue = 0.0;
  }
  else
  {
    var rnorm = (maxval - red) / mdiff;
    var gnorm = (maxval - green) / mdiff;
    var bnorm = (maxval - blue) / mdiff;
    saturation = (luminance <= 0.5) ? (mdiff / msum) : (mdiff / (510.0 - msum));
    if (red == maxval)
      hue = 60.0 * (6.0 + bnorm - gnorm);
    if (green == maxval)
      hue = 60.0 * (2.0 + rnorm - bnorm);
    if (blue == maxval)
      hue = 60.0 * (4.0 + gnorm - rnorm);
    if (hue > 360.0)
      hue -= 360.0;
  }
  return new Array (Math.round (hue * 255.0 / 360.0), Math.round (saturation * 255.0), Math.round (luminance * 255.0));
}

function Magic (rm1, rm2, rh)
{
  var retval = rm1;
  if (rh > 360.0)
    rh -= 360.0;
  if (rh < 0.0)
    rh += 360.0;
  if (rh < 60.0)
    retval = rm1 + (rm2 - rm1) * rh / 60.0;
  else if (rh < 180.0)
    retval = rm2;
  else if (rh < 240.0)
    retval = rm1 + (rm2 - rm1) * (240.0 - rh) / 60.0;
  return Math.round (retval * 255);
}

//Retourne un tableau de 3 valeurs : R,G,B
function HSL2RGB (h, s, l)
{
  var hue = h * 360.0 / 255.0;
  var saturation = s / 255.0;
  var luminance = l / 255.0;
  var red;
  var green;
  var blue;
  if (saturation == 0.0)
  {
    red = green = blue = Math.round (luminance * 255.0);
  }
  else
  {
    var rm1;
    var rm2;
    if (luminance <= 0.5)
      rm2 = luminance + luminance * saturation;
    else
      rm2 = luminance + saturation - luminance * saturation;
    rm1 = 2.0 * luminance - rm2;
    red = Magic (rm1, rm2, hue + 120.0);
    green = Magic (rm1, rm2, hue);
    blue = Magic (rm1, rm2, hue - 120.0);
  }
  return new Array (red, green, blue);
}
*/

