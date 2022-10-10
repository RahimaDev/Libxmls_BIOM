
#define WHAT "EptReader: example usage of echo/pulse tables (ept)"

#include <ctime>
#include <iostream>
#include "libXMls/XEchoPulseTables.h"

using namespace std;

int main(int argc, char **argv)
{
    cout << WHAT << endl;
    if(argc != 2)
    {
        cout << "Usage: " << argv[0] << " ept_basename" << endl;
        return 1;
    }
    cout<<"Reading "<<argv[1]<<endl;

    clock_t start = clock();
    XEchoPulseTable ep_table(argv[1]);
    cout << ep_table.NPulseAttrib() << "/" << ep_table.NEchoAttrib() << " pulse/echo attributes found" << endl;
    cout << ep_table.NTotalPulse() << "/" << ep_table.NTotalEcho() << " pulse/echo attributes found" << endl;

    // Select all blocks = create blocks structure but do not load them, only required attributes are accessible
    ep_table.SelectAll();
    // you can also select only a certain time interval with ep_table.Select(XSecond i_start, XSecond i_end)
    cout << ep_table.NBlock() << " blocks selected" << endl;

    // this is the simplest way to make an existing attribute accessible
    XFloatAttrib * p_refl_attrib = ep_table.GetEchoAttrib<XFloatAttrib>("reflectance");
    XFloatAttrib * p_ampl_attrib = ep_table.GetEchoAttrib<XFloatAttrib>("amplitude");
    XUCharAttrib * p_dev_attrib = ep_table.GetEchoAttrib<XUCharAttrib>("deviation");
    cout << "p_refl_attrib->m_ept_path=" << p_refl_attrib->EptPath() << endl;

    // this is the simplest way to create a new (accessible) attribute and to access it
    XUCharAttrib * p_lum_attrib = ep_table.AddEchoAttrib<XUCharAttrib>("lum");
    cout << "p_lum_attrib->m_ept_path=" << p_lum_attrib->EptPath() << endl;

    // display info on (accessible) attributes and (selected) blocks
    for(auto& it:ep_table.AttribList())
        cout << it->m_object << " attrib " << it->m_attribname << " type "
             << it->Typename() << " folder " << it->AttributeFolder() <<
                " accessible " << it->m_is_accessible << " NBlocks " << it->NBlock() << endl;

    double time_min=1.e30, time_max=0.;
    float range_min=1.e30, range_max=0.;
    float refl_min=1.e30, refl_max=-1.e30;
    float ampl_min=1.e30, ampl_max=-1.e30;
    unsigned char dev_min=255, dev_max=0;
    XPt3D G;
    int n_min=8, n_max=0;

    // example iteration on blocks
    for(XBlockIndex block_idx=0; block_idx<ep_table.NBlock(); block_idx++)
    {
        cout << "Processing block " << block_idx << " sec " << ep_table.mts_echo.Second(block_idx) << endl;
        ep_table.Load(block_idx);

        // example iteration on pulses
        for(XPulseIndex i_pulse=0; i_pulse<ep_table.NPulse(block_idx); i_pulse++)
        {
            double t=ep_table.Time(block_idx, i_pulse);
            if(t < time_min) time_min=t;
            if(t > time_max) time_max=t;
        }

        // example iteration on echos
        for(XEchoIndex i_echo=0; i_echo<ep_table.NEcho(block_idx); i_echo++)
        {
            double r = ep_table.Range(block_idx, i_echo);
            if(r < range_min) range_min=r;
            if(r > range_max) range_max=r;
            float refl = p_refl_attrib->at(block_idx)[i_echo];
            float ampl = p_ampl_attrib->at(block_idx)[i_echo];
            unsigned char dev = p_dev_attrib->at(block_idx)[i_echo];
            p_lum_attrib->at(block_idx)[i_echo] = (refl+20.)*12.75;
            if(refl < refl_min) refl_min=refl;
            if(refl > refl_max) refl_max=refl;
            if(ampl < ampl_min) ampl_min=ampl;
            if(ampl > ampl_max) ampl_max=ampl;
            if(dev < dev_min) dev_min=dev;
            if(dev > dev_max) dev_max=dev;
            int n = ep_table.NumEcho(block_idx, i_echo);
            if(n < n_min) n_min=n;
            if(n > n_max) n_max=n;
            G += ep_table.P(block_idx, i_echo);
            if(i_echo % 1000000 == 0) cout << "P=" << ep_table.P(block_idx, i_echo) << endl;
        }
        p_lum_attrib->Save(block_idx); // save the new attrib before freeing a block (else its lost)
        ep_table.Free(block_idx);
    }
    cout << "Time " << time_min << "-" << time_max << endl;  
    cout << "Range " << range_min << "-" << range_max << endl;
    cout << "Reflectance " << refl_min << "-" << refl_max << endl;
    cout << "Amplitude " << ampl_min << "-" << ampl_max << endl;
    cout << "Deviation " << dev_min << "-" << dev_max << endl;
    cout << "N echo " << n_min << "-" << n_max << endl;
    cout << "Barycenter " << G / ep_table.NTotalEcho() << endl;
    cout << "Done in " << (double)(clock() - start) / (double) CLOCKS_PER_SEC << "s" << endl;
    return 0;
}
