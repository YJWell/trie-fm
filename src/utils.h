#ifndef UTIL_H
#define UTIL_H


#include <tr1/unordered_map>
#include <stdio.h>
#include <string>
#include <iostream>
#include <vector>
#include <time.h>
#include <algorithm>
#include <fstream>
#include <stack>
#include <math.h>

using namespace std;


#define CHILDNUM 4 //# of children for each node
#define MAXREADLEN 200

#define KMER 17
#define BIT_SET(a,b) ((a) |= (1<<(b)))
#define BIT_CLEAR(a,b) ((a) &= ~(1<<(b)))
#define BIT_FLIP(a,b) ((a) ^= (1<<(b)))
#define BIT_CHECK(a,b) ((a) & (1<<(b)))
typedef int32_t saint64_t;
typedef u_int32_t fmint_t;


//int gettimeofday(struct  timeval*tv,struct  timezone *tz );

void complementSequence(string &seq)
{
    for(auto it=seq.begin();it<seq.end();it++)
    {
        if(*it=='A' || *it=='a')
            *it='T';
        else if(*it=='T' || *it=='t')
            *it='A';
        else if(*it=='C' || *it=='c')
            *it='G';
        else if(*it=='G' || *it=='g')
            *it='C';
    }
}

u_int8_t encode(char &c)
{
    u_int8_t i;
    switch(c)
    {
        case 'a':
        case 'A': i=0;break;
        case 'c':
        case 'C': i=1;break;
        case 'g':
        case 'G': i=2;break;
        case 't':
        case 'T': i=3;break;
        default: i=0;break;
    }
    return i;
}

char decode(size_t i)
{
    char c;
    switch(i)
    {
        case 0: c='A';break;
        case 1: c='C';break;
        case 2: c='G';break;
        case 3: c='T';break;
        default: c='N';break;
    }
    return c;
}

class FastaHead
{
public:
    string name;
    fmint_t length;
    fmint_t start; // if fasta is reversed then (start,end],
    fmint_t end;   // else, [start,end)

    FastaHead(string s,fmint_t l)
    {
        name=s;
        length=l;
        start=-1;
        end=-1;
    }
    /*
private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & name;
        ar & length;
        ar & start;
        ar & end;
    }

    friend class boost::serialization::access;//serialization
    */
};

class RealPosition
{
public:
    string name;
    fmint_t pos;
};

class ChromIndex
{
public:
    ChromIndex()
    {
        READY=false;
    }

    void addhead(FastaHead head)
    {
        heads.push_back(head);
    }

    void adjust()
    {
        fmint_t tmp_start=0;
        for(size_t i=heads.size()-1;i>0;i--)
        {
            heads[i].start=tmp_start;
            heads[i].end=tmp_start+heads[i].length;
            tmp_start=heads[i].end+1;
        }
        heads[0].start=tmp_start;
        heads[0].end=tmp_start+heads[0].length;
        READY=true;
    }

    RealPosition calRealPos(fmint_t &p,size_t l)
    {
        if(!READY)
            adjust();
        for(size_t i=0;i<heads.size();i++)
        {
            if(p>heads[i].start)
            {
                if(p<=heads[i].end)
                {
                    rpos.pos=heads[i].length-p+heads[i].start-l+2;
                    rpos.name=heads[i].name;
                }
                else
                {
                    //cout<<"read wrong, because it matches to 'end marker'"<<endl;
                    rpos.pos=p;
                    rpos.name="WrongRead";
                }
                break;
            }
        }
        return rpos;
    }

private:
    /*
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & heads;
        ar & rpos;
        ar & READY;
    }*/

    vector<FastaHead> heads;
    RealPosition rpos;
    bool READY;
    //friend class boost::serialization::access;//serialization
};




/*
void storeIndiceToDisk(FMIndex &fmidx,ChromIndex &chrindex,string fname)
{
    string out_name1=fname+".1.fm";
    ofstream ofs1(out_name1);
    boost::archive::binary_oarchive oa1(ofs1);
    oa1<<fmidx;
    string out_name2=fname+".2.fm";
    ofstream ofs2(out_name2);
    boost::archive::binary_oarchive oa2(ofs2);
    oa2<<chrindex;
}

FMIndex readFMIndexFromDisk(string fname)
{
    string in_name=fname+".1.fm";
    FMIndex fmindex;
    ifstream ifs(in_name);
    boost::archive::binary_iarchive ia(ifs);
    ia>>fmindex;
    return fmindex;
}


ChromIndex readChrindexxFromDisk(string fname)
{
    string in_name=fname+".2.fm";
    ChromIndex chrindex;
    ifstream ifs(in_name);
    boost::archive::binary_iarchive ia(ifs);
    ia>>chrindex;
    return chrindex;
}
*/



#endif // UTIL_H
