#pragma once

#include <string>
#include <vector>
#include <fstream>
#include "proj_api.h"
#include "libXBase/XArchiGeoref.h"

/**
 * @brief POSPac post-processed (position,speed,acceleration) file (sbet_xxx.out) format
 * @note Time is expressed as seconds of week in GPS time system. GPS week starts
 * at midnight between saturday and sunday.
 */
struct SbetEvent
{
    double m_time;//seconds of week in GPS time system
    double m_latitude;//radians
    double m_longitude;//radians
    double m_altitude;//meters
    double m_xVelocity;//meters/second
    double m_yVelocity;//meters/second
    double m_zVelocity;//meters/second
    double m_roll;//radians
    double m_pitch;//radians
    double m_plateformHeading;//radians
    double m_wanderAngle;//radians
    double m_xAcceleration;//meters/second2
    double m_yAcceleration;//meters/second2
    double m_zAcceleration;//meters/second2
    double m_xBodyAngularRate;//radians/second
    double m_yBodyAngularRate;//radians/second
    double m_zBodyAngularRate;//radians/second
};
/// linear interpolation in time (only positional info)
SbetEvent InterpolPosOnly(const SbetEvent & s0, SbetEvent & s1, double time);
/// linear interpolation in time
SbetEvent Interpol(const SbetEvent & s0, SbetEvent & s1, double time);
std::ostream& operator<<(std::ostream& out, const SbetEvent& s);


/**
 * @brief POSPac post-processed accuracy file (smrmsg_xxx.out where xxx
 * is the processing Kernel) format
 * @note Time is expressed as seconds of week in GPS time system. GPS week starts
 * at midnight between saturday and sunday.
 */
struct AccuracyEvent
{
    double m_time;//seconds of week in GPS time system
    double m_northPositionRMSError;//meters
    double m_eastPositionRMSError;//meters
    double m_downPositionRMSError;//meters
    double m_northVelocityRMSError;//meters/second
    double m_eastVelocityRMSError;//meters/second
    double m_downVelocityRMSError;//meters/second
    double m_RollRMSError;//arc-minutes/second
    double m_PitchRMSError;//arc-minutes/second
    double m_headingRMSError;//arc-minutes/second
};
// linear interpolation between a0 (alpha=0) and a1 (alpha=1)
AccuracyEvent Interpol(const AccuracyEvent & a0, AccuracyEvent & a1, double time);
std::ostream& operator<<(std::ostream& out, const AccuracyEvent& a);

/// mutualize code for reading/saving sbet (position,speed acceleration) and smrmsg (accuracy)
template <class T> class TEventSeries
{
protected:
    std::string m_filename;
    std::vector<T> mv_event;
    double m_start_time, m_end_time;
    unsigned int m_n_event; //, m_nb_jour;
    bool m_meta_info_present;
    int m_last_index; // for optimization in sequential access

public:
    bool HasFile(){return m_filename.size()>0;}
    TEventSeries():m_filename(""), m_start_time(0.), m_end_time(0.),
        m_n_event(0), m_meta_info_present(false),m_last_index(0){}
    TEventSeries(std::string filename):m_filename(filename),m_meta_info_present(false),m_last_index(0){if(HasFile()) GetMetaInfo();}

    void SetFilename(std::string filename){m_filename=filename;GetMetaInfo();}
    std::string Filename(){return m_filename;}

    int Nevent(){return m_n_event;}
    //int NbJour(){return m_nb_jour;}
    std::vector<T> & Events(){return mv_event;}
    T & Event(int i){return mv_event.at(i);}
    /// At construction: of the whole file, after a load: of the loaded interval
    double StartTime(){return m_start_time;}
    /// At construction: of the whole file, after a load: of the loaded interval
    double EndTime(){return m_end_time;}
    void Display(int n_lines=7)
    {
        std::cout.precision(16);
        for(int i=0; i<mv_event.size(); i+=mv_event.size()/n_lines)
            std::cout << mv_event[i] << std::endl;
    }

    /// get index of first block for which t < time
    int PrevIndex(double time)
    {
        if(!m_meta_info_present) return 0;
        if(time < m_start_time) {
            std::cout.precision(16); std::cout << "Time: " << time << "<" << m_start_time << std::endl; return -1;
        }
        if(time > m_end_time) {
            std::cout.precision(16); std::cout << "Time: " << time << ">" << m_end_time << std::endl; return m_n_event-1;
        }
        // optimization for sequential access: check if the last index commputed is good
        if(time > mv_event[m_last_index].m_time && time < mv_event[m_last_index+1].m_time )
            return m_last_index;
        return PrevIndex(time, 0, m_n_event-1);
    }

    /// recursive part, assumes time is in the i_min, i_max interval
    int PrevIndex(double time, int i_min, int i_max)
    {
        if(!m_meta_info_present) return 0;
        if(i_max <= i_min+1)
        {
            m_last_index = i_min;
            return i_min; // end recursion
        }
        int i_med = (i_min+i_max)/2;
        if(time < mv_event[i_med].m_time)
            return PrevIndex(time, i_min, i_med);
        else return PrevIndex(time, i_med, i_max);
    }

    /// returns false if the time was outside bounds (in this case block is the first or last)
    bool Interpol_event(T & event, double time)
    {
        if(!m_meta_info_present) {return false;}
        int index = PrevIndex(time);
        if(index < 0) {event = mv_event.front(); return false;}
        if((unsigned int) index+1 >= m_n_event) {event = mv_event.back(); return false;}
        event = Interpol(mv_event[index], mv_event[index+1], time);
        //std::cout << "time " << time << " index " << index << std::endl;
        return true;
    }

    /// get meta info (number of events, start and end time) without loading everything
    bool GetMetaInfo()
    {
        if(m_meta_info_present || !HasFile()) return true;
        std::ifstream ifs(m_filename.c_str(), std::ios::binary | std::ios_base::in | std::ios_base::ate);
        if(!ifs.good() || !ifs.is_open() || ifs.fail())
        {
            std::cout << "XTrajecto::GetMetaInfo() Failed to open " << m_filename << std::endl;
            return false;
        }

        // number of event is infered from file size
        m_n_event = ifs.tellg()/sizeof(T);

        // read start and end time in the file
        ifs.seekg(-(int)sizeof(T), std::ios::end);
        ifs.read((char*)&m_end_time, sizeof(double));
        ifs.seekg(0, std::ios::beg);
        ifs.read((char*)&m_start_time, sizeof(double));
        //m_nb_jour=(int)(m_start_time/24.0/3600.0);
        m_meta_info_present = true;
        return true;
    }

    bool Load(double t_min=0., double t_max=0.)
    {
        if(!GetMetaInfo()) return false;
        std::ifstream ifs(m_filename.c_str(), std::ios::binary | std::ios_base::in);
        if(!ifs.good() || !ifs.is_open() || ifs.fail())
        {
            std::cout << "XTrajecto::Load() Failed to open " << m_filename << std::endl;
            return false;
        }
        if(t_max < t_min) // swap
        {
            double temp = t_min; t_min=t_max; t_max=temp;
        }
        int start_id = 0;
        int end_id = m_n_event-1;
        if((t_min == 0. && t_max == 0) || (t_min < m_start_time && t_max > m_end_time)) // load all trajecto
        {
            std::cout << "Reading " << m_n_event << " events" << std::endl;
            mv_event.resize(m_n_event);
            ifs.read((char*)&mv_event[0], sizeof(T)*mv_event.size());
        }
        else // partial load mode
        {
            // add a margin because linear estimation is not perfect
            t_min = t_min-1.; t_max=t_max+1;
            if(t_max < m_start_time || t_min > m_end_time)
            {
                std::cout << "Required interval [" << t_min << " " << t_max <<
                             "] is not in trajectory interval [" << m_start_time << " " << m_end_time <<
                             "], check your time interval" << std::endl;
                return false;
            }
            if(t_min < m_start_time) t_min = m_start_time;
            if(t_max > m_end_time) t_max = m_end_time;
            start_id = (t_min - m_start_time)*m_n_event/(m_end_time-m_start_time);
            end_id = (t_max - m_start_time)*m_n_event/(m_end_time-m_start_time)+1; // upper bound
            m_n_event = end_id - start_id + 1; // start == end => 1 event
        }
        mv_event.resize(m_n_event);
        ifs.seekg(sizeof(T)*start_id, std::ios::beg);
        ifs.read((char*)&mv_event[0], sizeof(T)*mv_event.size());
        m_start_time = mv_event.front().m_time;
        m_end_time = mv_event.back().m_time;
        //std::cout << "Loaded " << m_n_event << " events from index " << start_id << "(time " << m_start_time << ") to "
        //          << end_id << "(time " << m_end_time << ")" << std::endl;
        ifs.close();
        return true;
    }

    //-----------------------------------------------------------------------------
    bool Unload()
    {
        mv_event.clear();
        return true;
    }
    bool Save(std::string filename)
    {
        if(mv_event.empty()) return false;
        std::ofstream ofs(filename.c_str(), std::ios::binary | std::ios_base::out);
        if(!ofs.good() || !ofs.is_open() || ofs.fail())
        {
            std::cout << "XTrajecto::Save() Failed to open " << filename << std::endl;
            return false;
        }
        ofs.write((char*)&mv_event[0], sizeof(T)*mv_event.size());
        ofs.close();
        return true;
    }
};

/// adding georeferencing functionalities to sbet events
class XGeorefEventSeries:public TEventSeries<SbetEvent>
{
    projPJ m_proj_in, m_proj_out;
    std::vector<XArchiGeoref> mv_georef;

public:
    XGeorefEventSeries(){}
    XGeorefEventSeries(std::string filename,
                       std::string strProj4_out="+init=IGNF:LAMB93",
                       std::string strProj4_in="+init=IGNF:RGF93G");

    /// Set the input projection (the one used in the raw trajecto file)
    void SetInputProjection(std::string strProj4_in="+init=IGNF:LAMB93");

    /// Set the output projection (the one used to answer all geometric queries)
    void SetOutputProjection(std::string strProj4_in="+init=IGNF:LAMB93");

    /// Get the meridian convergence angle
    double GetConvMeridien(SbetEvent event);

    /// Precompute the georef for (required for GetGeoref_precomputed)
    bool PreComputeGeoref();

    /// Create a georeference from an event
    void GetGeoref(const SbetEvent & event, XArchiGeoref & georef);

    /// direct access once precomputed
    inline XArchiGeoref & GetGeoref(unsigned int i){return mv_georef.at(i);}

    /// Get the translation only
    XPt3D GetTranslation(const SbetEvent & event);
    inline XPt3D  GetTranslation(int i){return GetTranslation(Event(i));}
    bool GetTranslation(double time, XPt3D & P, bool precomputed=true);

    /// Get the rotation only
    XMat3D GetRotation(const SbetEvent & event);
    inline XMat3D  GetRotation(int i){return GetRotation(Event(i));}

    /// Create a georeference by linear interpolation in time between events then georeferencing
    /** @param time seconds of week in GPS time */
    bool GetGeoref(double time, XArchiGeoref & georef);

    /// Create a georeference by linear interpolation in time between precomputed georefs
    /** @param time seconds of week in GPS time */
    bool GetGeoref_precomputed(double time, XArchiGeoref & georef);

    /// Apply the precomputed georef by linear interpolation between the point with previous and next sampled georef
    /** @param time seconds of week in GPS time */
    XPt3D ApplyGeoref_precomputed(double time, const XPt3D & P);

    /// export INS info in a ply file
    void ExportPly(std::string ply_filename);
    /// export INS info in a ply file, keep only some events (keep size should be n_event)
    void ExportPly(std::string ply_filename, std::vector<bool> & keep);

    /// apply a translation to the georefevent
    void ApplyTranslation(unsigned int i, const XPt3D & P);
};

/// Top level classes are prefixed with X
typedef TEventSeries<AccuracyEvent> XAccuracy;

/**
 * @brief POSPac post-processed solution (position,speed,acceleration) + accuracies
 */
class XTrajecto:public XGeorefEventSeries
{  
protected:
    XAccuracy m_accuracy;

public:
    XTrajecto(std::string sbetFile, std::string accFile="",
              std::string strProj4_out="+init=IGNF:LAMB93",
              std::string strProj4_in="+init=IGNF:RGF93G");

    void SetParams(std::string sbetFile, std::string accFile);

    std::string AccuracyFile(){return m_accuracy.Filename();}
    bool HasAcc(){return m_accuracy.HasFile();}

    XAccuracy & Accuracy(){return m_accuracy;}

    /// get meta info (number of events, start and end time) without loading everything
    bool GetMetaInfo();
    /// loads all trajecto (default values) or a time interval (given from day or week start)
    bool Load(double t_min=0., double t_max=0.);
    bool Unload();
    bool Save(std::string pos_filename, std::string acc_filename="");

    /// returns false if the time was outside bounds (in this case block is the first or last)
    bool Interpol(SbetEvent & sbet_event, AccuracyEvent & acc_event, double time);
};
