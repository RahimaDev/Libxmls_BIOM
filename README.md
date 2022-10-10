# libXMLS README
Mobile Laser Scanning library by Bruno Vallet for IGN-France

Required dependencies: TinyXML, Boost (date_time system filesystem)

Optional dependencies: Ori CGAL

Test:

Unzip data.zip at the root

Assuming your build dir is next to the root, go to it and:

./bin/XMlsInfo ../libxmls/data/traj/ 20140616 ../libxmls/data/Calib.xml ../libxmls/data/test.ept

the output should end with:

Range in [2.124,446.8]

Amplitude in [0.56,70.98]

Reflectance in [-20.29,38.2]

Deviation in [0,255]

numEcho in [0,3]

X in [-30.14,103.6]

Y in [-446.5,413.4]

Z in [-2.096,1.175e-38]

Xw in [6.508e+05,6.516e+05]

Yw in [6.861e+06,6.862e+06]

Zw in [-56.51,75.76]

#./bin/XMlsInfo ../data/traj/20140616 ../data/Calib.xml ../data/test.ept


# Detecting points inside buildings in the outside scan (by Rahima Djahel)

Download the necessary data available at:




./bin/Polygon ../data/trajecto/sbet_EMS-160624-LC_1.out     ../data/CalibRiegl_clone.xml     ../data/EMS-20160624_0755-01-00003.ept     25859 0

with:

sbet_EMS-160624-LC_1.out= sbet file
CalibRiegl_clone.xml= calibration file for the laser
EMS-20160624_0755-01-00003.ept= folder containing the echo pulse tables (generated with EptExport)
25859=id of the block to process (look at info.txt in ept_folder for valid range)
0=parameter to enable/disable a contrario method (0 disable/1 enable)

# Openings detection (by Rahima Djahel)

./bin/Openings ../data/trajecto/sbet_EMS-160624-LC_1.out     ../data/CalibRiegl_clone.xml     ../data/EMS-20160624_0755-01-00003.ept     25727 1

# Related papers:
Openings detection:

Djahel, R., Vallet, B., & Monasse, P. (2022). DETECTING OPENINGS FOR INDOOR/OUTDOOR REGISTRATION. The International Archives of Photogrammetry, Remote Sensing and Spatial Information Sciences, 43, 177-184.

Detecting points inside buildings in the outside scan:

Djahel, R., Vallet, B., & Monasse, P. (2021). Towards Efficient Indoor/Outdoor Registration Using Planar Polygons. ISPRS annals of the photogrammetry, remote sensing and spatial information sciences, 2, 51-58.
