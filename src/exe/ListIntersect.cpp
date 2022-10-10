
#define WHAT "Intersects two lists (id_pt, id_mes, time)"

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <set>

using namespace std;

void RobustGetLine(ifstream & ifs, string & line)
{
    getline(ifs, line);
    if(line[line.size()-1] == '\r')
        line.erase(line.size()-1);
}

struct Point
{
    int id;
    double x, y, z;
    Point():id(0),x(0),y(0),z(0){}
};

bool operator<(const Point & p1, const Point & p2)
{
    return p1.id < p2.id;
}

struct Entry
{
    Point pt;
    int id_match;
    float time;
    Entry(): id_match(0), time(0){}
};

std::vector<Entry> Read(string filename)
{
    vector<Entry> ret;
    cout.precision(9);

    // open file
    ifstream ifs(filename.c_str());
    if(!ifs.good()) cout << "ListIntersect::Read() Failed to open " << filename << endl;
    string line;
    int i_line=0, id_match = 0, cur_id_pt=0;
    bool end_reached = false;
    while(!ifs.eof() && !end_reached && i_line++<100000)
    {
        RobustGetLine(ifs, line);
        if(!line.empty())
        {
            istringstream iss(line);
            Entry e;
            iss >> e.pt.id >> e.time >> e.pt.x >> e.pt.y >> e.pt.z;
            cout << e.pt.id << " " << e.time << endl;
            if(e.pt.id != cur_id_pt) // new point
            {
                cur_id_pt = e.pt.id;
                id_match = 0;
            } else id_match++;
            e.id_match = id_match;
            ret.push_back(e);
        }
    }
    return ret;
}

//-----------------------------------------------------------------------------
int main(int argc, char **argv)
{
    cout << WHAT << endl;
    cout << "argc: " << argc << endl;
    if(argc < 2)
    {
        cout << "Usage: " << argv[0] << "  list1.txt list2.txt [time_thr=5]" << endl;
        cout << "time_thr: threshold to make matches base on time (default=5s)" << endl;
        return 0;
    }

    int i_arg=1;
    // required
    string list1(argv[i_arg++]);
    cout << "list1: " << list1 << endl;
    string list2(argv[i_arg++]);
    cout << "list2: " << list2 << endl;

    // optional
    float time_thr=0.5;
    if(i_arg < argc) time_thr = atof(argv[i_arg++]);

    clock_t start = clock();
    vector<Entry> v1=Read(list1), v2=Read(list2);
    int n_match = 0;
    set<Point> pt_set, pt1_set, pt2_set;
    cout << "id_pt\tid1\tid2\ttime1\ttime2" << endl;
    for(auto & it1:v1) for(auto & it2:v2)
    {
        if(it1.pt.id == it2.pt.id && fabs(it1.time-it2.time) < time_thr)
        {
            cout << it1.pt.id << '\t' << it1.id_match << '\t' << it2.id_match <<
                    '\t' << it1.time << '\t' << it2.time << endl;
            n_match++;
            pt_set.insert(it1.pt);
        }
    }
    for(auto & it1:v1) pt1_set.insert(it1.pt);
    for(auto & it2:v2) pt2_set.insert(it2.pt);

    cout << "Visible" << endl;
    for(auto & it:pt_set) cout << it.id << " " << it.x << " " << it.y << " " << it.z << endl;

    cout << v1.size() << "," << v2.size() << "->" << n_match << " common matches " <<
            pt1_set.size() << "," << pt2_set.size() << "->" << pt_set.size() << " common points in " << (double)(clock() - start) / (double) CLOCKS_PER_SEC << "s" << endl;
    return 0;
}


