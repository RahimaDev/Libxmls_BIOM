
#define WHAT "XTrajectoSample: sample usage of XTrajecto"

#include "libXMls/XTrajecto.h"

#include <iostream>
#include <fstream>
#include <ctime>

using namespace std;

int main(int argc, char* argv[])
{
    cout << WHAT << endl;
    if(argc < 4)
    {
        cout << "Usage: " << argv[0] << " sbet_name smrmsg_name n_queries [ply_name]" << endl;
        return 1;
    }
    cout.precision(12);
    XTrajecto trajecto(argv[1], argv[2]);
    int n_queries = atoi(argv[3]);
    string ply_name;
    if(argc > 4) ply_name = string(argv[4]);
    clock_t start = clock();
    trajecto.Load();
    cout << endl << trajecto.Nevent() << " events from " <<
            trajecto.StartTime() << " to " << trajecto.EndTime()  << " read in "
         << (double)(clock() - start) / (double) CLOCKS_PER_SEC << "s" << endl;
    start = clock();
    trajecto.PreComputeGeoref();
    cout << trajecto.Nevent() << " geotransforms in "
         << (double)(clock() - start) / (double) CLOCKS_PER_SEC << "s" << endl;

    // test queries
    if(n_queries>0)
    {
        vector<double> time;

        int progressbar_size=10;
        int progress_step = n_queries/progressbar_size+1;
        double start_time=trajecto.StartTime()+2.;
        double end_time=trajecto.EndTime()-2.;
        double dt=(end_time-start_time)/n_queries;
        for(int i=0; i<n_queries; i++)
            time.push_back(start_time+i*dt);

        cout << '|';
        for(int i=1;i<progressbar_size-1;i++) cout << '-';
        cout << '|' << endl;
        start = clock();
        for(unsigned int i=0;i<time.size();i++)
        {
            //SolutionEvent sol_event;
            //bool OK = trajecto.Interpol_event_TimeFromWeek(sol_event, time.at(i));
            XArchiGeoref georef;
            bool OK = trajecto.GetGeoref_precomputed(time.at(i), georef); //
            if(OK && i%progress_step==0)
            {
                //cout << '.' << flush;
                cout << "Solution Event: "<< time.at(i) << (OK?" OK":" KO") << endl;
                cout << georef.Translation() << endl << georef.Rotation() << endl;
            }
        }
        //trajecto.Display();
        cout << endl << n_queries << " queries in " << (double)(std::clock() - start) / (double) CLOCKS_PER_SEC << "s" << endl;
    }
    if(ply_name.size()>0) trajecto.ExportPly(ply_name);
    return 0;
}

