#include <iostream>
#include <iomanip>
#include <ctime>
#include <sstream>
#include <string>
#include <list>
#include <vector>
#include <map>
#include <algorithm>
#include <cstdlib>
#include <iterator>
#include <set>
#include <fstream>
#include <chrono>
using namespace std;
//graph class declaration
class Graph
    {
    public:
        set<int> V; //only vertices with positive degree are stored
        map<int, float> strength;
        map<int, map<int, float> > weight;
        void read_edgelist(string&, bool, bool);
        inline int order(){ return V.size(); }
        int ecount();
        float get_weight(int, int);
        void print_graph();
        float proximity(int, int);
        float get_sn_nbrs_proximity(int, map<int, float>&);
        void memberships_to_coms(map<int, set<int> >&, map<int, set<int> >&);
        friend float proximity(Graph&, int, int);
        friend float overlap(Graph&, set<int>&, set<int>&);
    };
//other function declarations
void set_parameters(int, char*[], bool &, float &, float &);
void merge_communities(Graph&, map<int, set<int> >&, map<int, set<int> >&, float, int);
void write_seeds(set<int>&, string&);
float overlap(Graph& g, set<int>&, set<int>&);
void usage(void);
void print_set(set<int>&);

//setting the seed for random number generator. This generates new sequence at each run of the program.
void rand_seed()
{
    int seed = static_cast<int>(time(0));
    srand(seed);
}

int main(int argc, char* argv[])
{
    rand_seed();
    std::chrono::time_point<std::chrono::system_clock> start_time, end_time;
    start_time = std::chrono::system_clock::now();  //time starts
    bool weighted = false;
    float max_overlap = 0.4, rho_0 = 0.4;
    string network_file = argv[1];
    set_parameters(argc, argv, weighted, max_overlap, rho_0);
    cout<<"-------------------------------------------------------------------"<<endl;
    cout<<"Neighbourhood Proximity based Community Detection "<<endl;
    cout<<"Authors:  Pawan Kumar and Ravins Dohare"<<endl;
    cout<<"Email: pkumariitd@gmail.com, ravinsdohare@gmail.com"<<endl;
    cout<<"-------------------------------------------------------------------"<<endl;
    Graph g;
    map<int, set<int> > coms;
    g.read_edgelist(network_file, false, weighted);
    //g.print_graph();
    map<int, set<int> >::iterator ci;
    set<int> uncovered = g.V;
    map<int, set<int> > membership;
    set<int>::iterator si, sj;
    int com_count = 1;
    while(!uncovered.empty())
    {
        si = uncovered.begin();
        int r = rand()%uncovered.size();
        for(int i = 1; i <= r; ++i)
            ++si;
        map<int, float> prox;
        float max_prox;
        max_prox = g.get_sn_nbrs_proximity(*si, prox);
        map<int, float>::iterator mi;
        for(mi = prox.begin(); mi != prox.end(); ++mi)
        {
            if(mi->second > rho_0 || mi->second >= max_prox)
                membership[*si].insert(membership[mi->first].begin(), membership[mi->first].end());
        }
        if(membership[*si].empty())
        {
            membership[*si].insert(com_count);
            com_count++;
        }
        for(mi = prox.begin(); mi != prox.end(); ++mi)
             if(mi->second > rho_0 && membership[mi->first].empty())
             {
                 membership[mi->first].insert(membership[*si].begin(), membership[*si].end());
                 uncovered.erase(mi->first);
             }
        uncovered.erase(*si);
    }
    //membership to communities
    g.memberships_to_coms(membership, coms);   //convert memberships into communities
    int avg_csize = 0;
    for(ci = coms.begin(); ci != coms.end(); ++ci)
        {   //cout<<ci->second.size()<<", ";
            avg_csize += ci->second.size();
        }
    avg_csize /= coms.size();
    merge_communities(g, coms, membership, max_overlap, avg_csize);
    ostringstream comfile;
    comfile<<"./coms-npcd.txt";
    ofstream fout(comfile.str());
    if(!fout.is_open())
    {
        cout<<"Destination file for communities could not be opened.";
        exit(1);
    }
    for(ci = coms.begin(); ci != coms.end(); ++ci)
    {
        copy(ci->second.begin(), ci->second.end(), ostream_iterator<int>(fout, " "));
        fout<<endl;
    }
    fout.close();
    end_time = std::chrono::system_clock::now();  //time ends
    std::chrono::duration<double> elapsed_seconds = end_time-start_time;
    cout.setf(ios::left, ios::adjustfield);
    cout<<"-------------------------------------------------------------------"<<endl;
    cout<<setw(33)<<"Network file"<<"= "<<setw(20)<<network_file<<endl;
    cout<<setw(33)<<"Network order"<<"= "<<setw(20)<<g.order()<<endl;
    cout<<setw(33)<<"No. of edges"<<"= "<<setw(20)<<g.ecount()<<endl;
    cout<<setw(33)<<"Total non-singleton communities"<<"= "<<coms.size()<<endl;
    cout<<setw(33)<<"Community file"<<"= "<<setw(20)<<comfile.str()<<endl;
    cout<<setw(33)<<"Time elapsed"<<"= "<<elapsed_seconds.count()<<"s\n";
    cout<<"-------------------------------------------------------------------"<<endl;
    return 0;
}

//here are the definitions of all the functions used in the code above
void set_parameters(int argc, char* argv[], bool &weighted, float &max_overlap, float &rho_0)
{
    if(argc < 2 || argc > 7)
        usage();
    if(argv[1][0] == '-')
        usage();
    int i=2;
    while(i < argc)
    {
        string arg = string(argv[i]);
        if(arg == "-w")
        {
            weighted = true;
            i++;
        }
        else
            if(arg == "-ov")
            {
                istringstream is(argv[i+1]);
                is>>max_overlap;
                if( max_overlap < 0 || max_overlap > 0.5)
                    usage();
                i += 2;
            }
            else
                if(arg == "-rh")
                {
                    istringstream is(argv[i+1]);
                    is>>rho_0;
                    if( rho_0 < 0 || rho_0 >= 1)
                        usage();
                    i += 2;
                }
                else
                    usage();
    }
}

void Graph::memberships_to_coms(map<int, set<int> >&membership, map<int, set<int> >&coms)
{
   for(auto si = V.begin(); si != V.end(); ++si)
        for(auto sj = membership[*si].begin(); sj != membership[*si].end(); ++sj)
            coms[*sj].insert(*si);
}

void Graph::read_edgelist(string& edgefile, bool directed, bool weighted)
{
    ifstream fin;
    fin.open(edgefile);
    if(!fin.is_open())
    {
        cout<<"The file containing the edgelist could not be opened."<<endl;
        exit(1);
    }
    string line;
    while ( getline(fin, line ) )
    {
        istringstream is( line );
        int u, v;
        is>>u;
        is>>v;
        if(u == v || weight[u].find(v) != weight[u].end())
            continue;
        float w;
        if(weighted == true)
        {
            if(is.str() == " " || is.eof())
            {   cout<<endl<<"The edge list has missing weights."<<endl;
                exit(1);
            }
            is>>w;
        }
        else
            w = 1;
        weight[u][v] = w;
        if(strength.find(u) == strength.end())
            strength[u] = w;
        else
            strength[u] = strength[u]+w;
        if(directed == false)
        {   weight[v][u] = w;
            if(strength.find(v) == strength.end())
                strength[v] = w;
            else
                strength[v] = strength[v]+w;
        }
        V.insert(u);
        V.insert(v);
    }
}
int Graph::ecount()
{
    map<int, map<int, float> >::iterator mi;
    int degree_sum = 0;
    for(mi = weight.begin(); mi != weight.end(); ++mi)
        degree_sum += mi->second.size();
    return degree_sum/2;
}
float Graph::get_weight(int u, int v)
{
    if(weight[u].find(v) != weight[u].end())
        return weight[u][v];
    else
        return 0;
}
float Graph::get_sn_nbrs_proximity(int u, map<int, float>& prox)
{
    map<int, float>::iterator mi, mj;
    float max_prox = 0;
    for(mi = weight[u].begin(); mi != weight[u].end(); ++mi)
    {
        prox[mi->first] = proximity(u, mi->first);
        if(prox[mi->first] > max_prox)
            max_prox = prox[mi->first];
        for(mj = weight[mi->first].begin(); mj != weight[mi->first].end(); ++mj)
            if(mj->first != u && prox.find(mj->first) == prox.end())
            {
                prox[mj->first] = proximity(u, mj->first);
                if(prox[mj->first] > max_prox)
                    max_prox = prox[mj->first];
            }
    }
    return max_prox;
}

float Graph::proximity(int u, int v)
{
    float min_weight_sum = 0, common_nbrs_strength = 0, weight_sum = 0;
    map<int, float>::iterator i, j;
    if(weight[u].size() > weight[v].size()) //to ensure that next loop runs minimum times.
    {
        int temp = u;
        u = v;
        v = temp;
    }
    set<int> common_nbrs;
    set<int>::iterator si;
    for(i = weight[u].begin(); i != weight[u].end(); ++i)
    {
        j = weight[v].find(i->first);
        if(j != weight[v].end())
        {
            min_weight_sum += min(i->second, j->second);
            weight_sum += i->second + j->second;
            common_nbrs_strength += strength[i->first];
            for(si = common_nbrs.begin(); si != common_nbrs.end(); ++si)
                weight_sum += 2*get_weight(i->first, *si);
            common_nbrs.insert(i->first);
        }
    }
    float prox;
    if(common_nbrs_strength == 0)
        prox = get_weight(u,v)/min(strength[u], strength[v]);
    else
        prox = (min_weight_sum + get_weight(u,v))*weight_sum/(min(strength[u], strength[v])*common_nbrs_strength);
    return prox;
}

void merge_communities(Graph& g, map<int, set<int> >& coms, map<int, set<int> >& members, float given_max_ov, int minc)
{
    map<int, set<int> >::iterator ci, cj;
    set<int>::iterator si, sj, sk;
    map<int, float>::iterator mi;
    set<int> labels;
    for(ci = coms.begin(); ci != coms.end(); ++ci)
        labels.insert(ci->first);
    do
    {
        si = labels.begin();
        int r = rand()%labels.size();
        for(int i = 1; i<=r; ++i)
            ++si;
        set<int> neighboring_com_labels;
        for(sj = coms[*si].begin(); sj != coms[*si].end(); ++sj)
            for(mi = g.weight[*sj].begin(); mi != g.weight[*sj].end(); ++mi)
                if(coms[*si].find(mi->first) == coms[*si].end())
                {
                    for(sk = members[mi->first].begin(); sk != members[mi->first].end(); ++sk)
                        if(coms[*sk].size() >= coms[*si].size())  //neighboring coms that are bigger than coms[*si]
                            neighboring_com_labels.insert(*sk);
                }
        float max_ov = 0;
        set<int> high_overlapping_coms;
        for(sj = neighboring_com_labels.begin(); sj != neighboring_com_labels.end(); ++sj)
        {
            float curr_ov = overlap(g, coms[*si], coms[*sj]);
            if(coms[*si].size() < minc)
            {
                if(curr_ov > max_ov)
                {
                    max_ov = curr_ov;
                    high_overlapping_coms.clear();
                    high_overlapping_coms.insert(*sj);
                }
                else
                    if(curr_ov == max_ov)
                        high_overlapping_coms.insert(*sj);
            }
            else
            {
                if( curr_ov >= given_max_ov)
                   high_overlapping_coms.insert(*sj);
            }
        }
        for(sj = high_overlapping_coms.begin(); sj != high_overlapping_coms.end(); ++sj)
        {
            for(sk = coms[*si].begin(); sk != coms[*si].end(); ++sk)
            {
                coms[*sj].insert(*sk);
                members[*sk].erase(*si);
                members[*sk].insert(*sj);
            }
        }
        if(!high_overlapping_coms.empty())
            coms.erase(*si);
        labels.erase(*si);
    }while(!labels.empty());
}

float overlap(Graph& g, set<int>& C1, set<int>& C2)
{
    set<int>::iterator si, sj;
    float intercom_weight = 0, deg1 = 0, deg2 = 0;
    for(si = C1.begin(); si != C1.end(); ++si)
    {
        for(sj = C2.begin(); sj != C2.end(); ++sj)
            intercom_weight += g.get_weight(*si, *sj);
        deg1 += g.strength[*si];
    }
    for(sj = C2.begin(); sj != C2.end(); ++sj)
        deg2 += g.strength[*sj];
    return intercom_weight/min(deg1, deg2);
}
void usage(void)
{
    cout<<endl<<"Call program as shown below:"<<endl;
    int i;
    for(i = 0; i < 22; ++i)
        cout<<"---";
    cout<<endl<<"./npcd network_file -w -ov [option] -rh [option]"<<endl;
    for(i = 0; i < 22; ++i)
        cout<<"---";
    cout<<endl<<"where:"<<endl;
    cout<<setw(15)<<"network_file"<<" --> The file containing the network in edge list format"<<endl;
    cout<<setw(15)<<"w"<<" --> A flag which if supplied indicates that the network is weighted"<<endl;
    cout<<setw(15)<<"ov"<<" --> Maximum allowable overlap between communities (between 0 and 1)"<<endl;
    cout<<setw(15)<<"rh"<<" --> Minimum neighborhood proximity for generation of initial communities(between 0 and 1)"<<endl;
    cout<<endl<<"The minimal syntax is:"<<endl;
    for(i = 0; i < 22; ++i)
        cout<<"---";
    cout<<endl<<"./npcd network_file"<<endl;
    for(i = 0; i < 22; ++i)
        cout<<"---";
    cout<<endl<<"In this case, the network is treated as unweighted and ov = 0.20 and rh = 0.4."<<endl;
    exit(1);
}
