///Author:Yujia Wu
///Created on Dec 20 2014
///Can build bwt for a sequence within 2^32-1 bp
///



#ifndef BWT_H
#define BWT_H

#include <vector>
#include <tr1/unordered_map>
#include <algorithm>
#include <string>
#include <array>
#include <utility>
#include <stack>
#include <iostream>

#include "divsufsort.h"

using namespace std;


typedef int32_t saint64_t;
typedef u_int32_t fmint_t;


class FMIndex
{
public:
    ///__init__
    FMIndex(fmint_t STEP=128,fmint_t GAP=16)
    {
        step=STEP>1?STEP:1;
        gap=GAP>1?GAP:1;
        EOS='$';
        EOC='%';
        best_bound_found=0;
        num_alignments_output=10;
        matches.reserve(1000);
        inex_bounds.reserve(500);
        best_bound.fill(0);
        inex_start_stat.resize(MAXREADLEN,0);
        bound_stack.fill(0);
        entries_of_bound_stack=(CHILDNUM+2)*2;
    }


    /*----------------------------------------------------------
     * For building bwt and fmindex
     * ---------------------------------------------------------*/

    ///transform string to bwt string
    void transform(string *orig_str)
    {
        orig_str->shrink_to_fit();
        size_t found=orig_str->find(EOS);
        if(found!=string::npos)
        {
            cout<<"original string contains $"<<endl;
            exit (EXIT_FAILURE);
        }
        for(string::iterator it=orig_str->begin();it<orig_str->end();it++)
        {
            if(*it=='N' || *it=='n')
            {
                *it='A';
            }
        }
        //add end mark to string
        orig_str->push_back(EOS);
        len=orig_str->size();
        full_sa=new saint64_t[len];
        divsufsort((unsigned char*)orig_str->c_str(),full_sa,len);
        //bwt_str.reserve(len);
        encoded_bwt_str.reserve(len);
        for(fmint_t i=0;i<len;i++)
        {
            switch(full_sa[i])
            {
                case 0: encoded_bwt_str.push_back(encode(EOS)); break;
                default: encoded_bwt_str.push_back(encode(orig_str->at(full_sa[i]-1))); break;
            }
        }
        encoded_bwt_str.shrink_to_fit();
    }


    ///build the index
    void buildIndex()
    {
        _calcCheckpointingSA(full_sa);
        _calcOccAndCheckpoints(encoded_bwt_str);
    }

    /*----------------------------------------------------------
     * For searching with trie
     * ---------------------------------------------------------*/

    void findBound(vector<u_int8_t> &q,fmint_t &top,fmint_t &bot)
    {
        vector<u_int8_t>::iterator qc;
        for(qc=q.begin();qc!=q.end();++qc)
        {
            top=_lf(top,*qc);
            bot=_lf(bot,*qc);
            if(top==bot)
            {
                //top=-1;
                //bot=-1;
                break;
            }
        }

    }


    void findBound(vector<u_int8_t> &q,const fmint_t start_pos,fmint_t &top,fmint_t &bot)
    {
        vector<u_int8_t>::iterator qc;
        for(qc=q.begin()+start_pos;qc!=q.end();++qc)
        {
            top=_lf(top,*qc);
            bot=_lf(bot,*qc);
            if(top==bot)
            {
                break;
            }
        }
    }


    vector<u_int8_t>::iterator findBound1(vector<u_int8_t> &q,fmint_t start_pos,fmint_t &top,fmint_t &bot)
    {
        fmint_t new_top,new_bot;
        vector<u_int8_t>::iterator qc;
        for(qc=q.begin()+start_pos;qc!=q.end();++qc)
        {
            new_top=_lf(top,*qc);
            new_bot=_lf(bot,*qc);
            if(new_top<new_bot)
            {
                top=new_top;
                bot=new_bot;
            }
            else
            {
                break;
            }
        }
        return qc;
    }

    ///trie inexact search the positions of query q,start from truncated position in leaf node
    vector<fmint_t>& trieInexSearch(vector<u_int8_t> &q,const fmint_t start_pos,fmint_t &top,fmint_t &bot,const int z)
    {
        //find the suffixes for the query
        match_status=0;
        matches.clear();
        vector<u_int8_t>::iterator qc=findBound1(q,start_pos,top,bot);//exact matched part
        inex_start_stat[qc-q.begin()]++;
        //cout<<"first index of failed matching: "<<qc-q.begin()<<endl;
        if(qc==q.end())//exact matched
        {
            if((bot-top)>num_alignments_output)
                bot=top+num_alignments_output;
            for(fmint_t i=top;i<bot;i++)
            {
                fmint_t pos=_walk(i);
                matches.push_back(pos);
            }
            match_status=1;
        }
        else
        {
            inex_bounds.clear();
            //_findInexBounds(q,qc,z,top,bot);
            _findInexBoundsOnceForAllLetters(q,qc,z,top,bot,0,false,false);
            for(auto inex_bound=inex_bounds.begin();inex_bound!=inex_bounds.end();inex_bound++)
            {
                if(inex_bound->second-inex_bound->first>num_alignments_output)
                    inex_bound->second=inex_bound->first+num_alignments_output;
                //cout<<inex_bound->first<<" "<<inex_bound->second<<endl;
                for(fmint_t i=inex_bound->first;i<inex_bound->second;i++)
                {
                    fmint_t pos=_walk(i);
                    matches.push_back(pos);
                }
            }
        }
        return matches;
    }

    ///trie inexact search the positions of query q, start traverse node->key, find all possible inexact matches
    vector<fmint_t>& trieInexSearchFromLeave(vector<u_int8_t> &q,const fmint_t start_pos,fmint_t &top,fmint_t &bot,const int z)
    {
        //find the suffixes for the query
        match_status=0;
        matches.clear();
        fmint_t leaf_top=top,leaf_bot=bot;
        vector<u_int8_t>::iterator qc=findBound1(q,start_pos,top,bot);//exact matched part
        inex_start_stat[qc-q.begin()]++;
        //cout<<"first index of failed matching: "<<qc-q.begin()<<endl;
        if(qc==q.end())//exact matched
        {
            if((bot-top)>num_alignments_output)
                bot=top+num_alignments_output;
            for(fmint_t i=top;i<bot;i++)
            {
                fmint_t pos=_walk(i);
                matches.push_back(pos);
            }
            match_status=1;
        }
        else
        {
            inex_bounds.clear();
            //_findInexBounds(q,qc,z,top,bot);
            _findInexBoundsOnceForAllLetters(q,q.begin()+start_pos,z,leaf_top,leaf_bot,0,false,false);
            for(auto inex_bound=inex_bounds.begin();inex_bound!=inex_bounds.end();inex_bound++)
            {
                if(inex_bound->second-inex_bound->first>num_alignments_output)
                    inex_bound->second=inex_bound->first+num_alignments_output;
                //cout<<inex_bound->first<<" "<<inex_bound->second<<endl;
                for(fmint_t i=inex_bound->first;i<inex_bound->second;i++)
                {
                    fmint_t pos=_walk(i);
                    matches.push_back(pos);
                }
            }
        }
        return matches;
    }

    ///trie inexact search the positions of query q,start from leave, find the first exact or inexact bound
    vector<fmint_t>& trieInexSearchFromLeaveFindBest(vector<u_int8_t> &q,const fmint_t start_pos,const fmint_t &top,const fmint_t &bot,const int z)
    {
        best_bound[0]=0;
        best_bound[1]=0;
        best_bound_found=0;
        _findBestBound(q,q.begin()+start_pos,z,top,bot);
        findAlignments(best_bound[0],best_bound[1]);
        return matches;
    }

   ///Do lf(occ+c) for all letters
   template<size_t size>
   void trieLFOnceForAllChildren(const fmint_t idx,const u_int8_t &children_flag,array<fmint_t,size> &trie_bound_stack,const int stack_level,const u_int8_t tb)
   {
       //find the nearest checkpoint for idx
       fmint_t check=(idx+(step/2))/step;
       if(check>=C_len)
           check=C_len-1;
       fmint_t pos=check*step;

       //range between pos and idx
       fmint_t range[2];
       if(pos>idx)
       {
           range[0]=idx;
           range[1]=pos;
       }
       else
       {
           range[0]=pos;
           range[1]=idx;
       }

       //count of all children letters between pos and idx
       u_int8_t letter;
       for(fmint_t i=range[0];i<range[1];i++)
       {
           letter=encoded_bwt_str[i];
           if(BIT_CHECK(children_flag,letter))
           {
               trie_bound_stack[stack_level*9+tb+letter]+=1;
           }
       }

       //calculate the letter count upto idx (not included)
       if(pos>idx)
       {
           for(size_t j=0;j<CHILDNUM;j++)
           {
               if(BIT_CHECK(children_flag,j))
               {
                    trie_bound_stack[stack_level*9+tb+j]=C[check][j]-trie_bound_stack[stack_level*9+tb+j]+occ[j];
               }
           }
       }
       else
       {
           for(size_t j=0;j<CHILDNUM;j++)
           {
               if(BIT_CHECK(children_flag,j))
               {
                    trie_bound_stack[stack_level*9+tb+j]+=C[check][j]+occ[j];
               }
           }
       }

   }

    vector<fmint_t>& findAlignments(const fmint_t top,const fmint_t bot)
    {
        matches.clear();
        if((bot-top)<num_alignments_output)
        {
            for(fmint_t i=top;i<bot;i++)
            {
                fmint_t pos=_walk(i);
                matches.push_back(pos);
            }
        }
        else
        {
            for(fmint_t i=top;i<top+num_alignments_output;i++)
            {
                fmint_t pos=_walk(i);
                matches.push_back(pos);
            }
        }
        //sort(matches.begin(),matches.end());
        return matches;
    }
    

    /*----------------------------------------------------------
     * For direct search
     * ---------------------------------------------------------*/
    ///search the positions of query q

    vector<fmint_t>& search(vector<u_int8_t> &q)
    {
        //find the suffixes for the query
        fmint_t top=0,bot=len;
        vector<u_int8_t>::iterator qc=_findBound(q,top,bot);//common process should use _findBound()
        // find the location of the suffixes
        // by walking the reverse text from that position
        // with lf mapping
        //cout<<"first index of failed matching: "<<qc-q.begin()<<endl;
        matches.clear();
        //inex_start_stat[qc-q.begin()]++;
        if(qc==q.end())
        {
            if((bot-top)>num_alignments_output)
                bot=top+num_alignments_output;
            for(fmint_t i=top;i<bot;i++)
            {
                fmint_t pos=_walk(i);
                matches.push_back(pos);
            }
            //sort(matches.begin(),matches.end());
        }
        return matches;
    }

    ///inexact search the positions of query q
    vector<fmint_t>& inexSearch(vector<u_int8_t> &q,const int z)
    {
        //find the suffixes for the query
        match_status=0;
        matches.clear();
        fmint_t top=0,bot=len;
        vector<u_int8_t>::iterator qc=_findBound(q,top,bot);//exact matched part
        inex_start_stat[qc-q.begin()]++;
        //cout<<"first index of failed matching: "<<qc-q.begin()<<endl;
        if(qc==q.end())//exact matched
        {
            if((bot-top)>num_alignments_output)
                bot=top+num_alignments_output;
            for(fmint_t i=top;i<bot;i++)
            {
                fmint_t pos=_walk(i);
                matches.push_back(pos);
            }
            match_status=1;
        }
        else
        {           
            inex_bounds.clear();
            _findInexBounds(q,qc,z,top,bot);
            //_findInexBoundsOnceForAllLetters(q,qc,z,top,bot,0,false,false);
            for(auto inex_bound=inex_bounds.begin();inex_bound!=inex_bounds.end();inex_bound++)
            {
                if(inex_bound->second-inex_bound->first>num_alignments_output)
                    inex_bound->second=inex_bound->first+num_alignments_output;
                //cout<<inex_bound->first<<" "<<inex_bound->second<<endl;
                for(fmint_t i=inex_bound->first;i<inex_bound->second;i++)
                {
                    fmint_t pos=_walk(i);
                    matches.push_back(pos);
                }
            }
        }
        return matches;
    }

    ///inexact search the positions of query q
    vector<fmint_t>& inexSearchFromStart(vector<u_int8_t> &q,const int z)
    {
        //find the suffixes for the query
        match_status=0;
        matches.clear();
        fmint_t top=0,bot=len;
        vector<u_int8_t>::iterator qc=_findBound(q,top,bot);//exact matched part
        inex_start_stat[qc-q.begin()]++;
        //cout<<"first index of failed matching: "<<qc-q.begin()<<endl;
        if(qc==q.end())//exact matched
        {
            if((bot-top)>num_alignments_output)
                bot=top+num_alignments_output;
            for(fmint_t i=top;i<bot;i++)
            {
                fmint_t pos=_walk(i);
                matches.push_back(pos);
            }
            match_status=1;
        }
        else
        {
            inex_bounds.clear();
            //_findInexBounds(q,q.begin(),z,0,len);
            _findInexBoundsOnceForAllLetters(q,q.begin(),z,0,len,0,false,false);
            for(auto inex_bound=inex_bounds.begin();inex_bound!=inex_bounds.end();inex_bound++)
            {
                if(inex_bound->second-inex_bound->first>num_alignments_output)
                    inex_bound->second=inex_bound->first+num_alignments_output;
                for(fmint_t i=inex_bound->first;i<inex_bound->second;i++)
                {
                    fmint_t pos=_walk(i);
                    matches.push_back(pos);
                }
            }
        }
        return matches;
    }

    ///inexact search the positions of query q, start from trunc, find the first exact or inexact bound
    vector<fmint_t>& inexSearchFromTruncFindBest(vector<u_int8_t> &q,const int z)
    {
        //find the suffixes for the query
        match_status=0;
        matches.clear();
        fmint_t top=0,bot=len;
        vector<u_int8_t>::iterator qc=_findBound(q,top,bot);//exact matched part
        inex_start_stat[qc-q.begin()]++;
        //cout<<"first index of failed matching: "<<qc-q.begin()<<endl;
        if(qc==q.end())//exact matched
        {
            if((bot-top)>num_alignments_output)
                bot=top+num_alignments_output;
            for(fmint_t i=top;i<bot;i++)
            {
                fmint_t pos=_walk(i);
                matches.push_back(pos);
            }
            match_status=1;
        }
        else
        {
            best_bound[0]=0;
            best_bound[1]=0;
            best_bound_found=0;
            _findBestBound(q,qc,z,top,bot);
            //_findInexBounds(q,qc,z,top,bot);
            //_findInexBoundsOnceForAllLetters(q,qc,z,top,bot,0,false,false);
            findAlignments(best_bound[0],best_bound[1]);
        }
        return matches;
    }

    ///inexact search the positions of query q, start from specified base position, find the first exact or inexact bound
    vector<fmint_t>& inexSearchFromSpecifiedPosFindBest(vector<u_int8_t> &q,const int z)
    {
        //find the suffixes for the query
        matches.clear();
        fmint_t top=0,bot=len;
        const uint pos=10;
        vector<u_int8_t>::iterator qc=_findBound(q,top,bot,pos);//exact matched part
        inex_start_stat[qc-q.begin()]++;
        //cout<<"first index of failed matching: "<<qc-q.begin()<<endl;
        best_bound[0]=0;
        best_bound[1]=0;
        best_bound_found=0;
        _findBestBound(q,qc,z,top,bot);
        //_findInexBounds(q,qc,z,top,bot);
        //_findInexBoundsOnceForAllLetters(q,qc,z,top,bot,0,false,false);
        findAlignments(best_bound[0],best_bound[1]);
        return matches;
    }
    ///inexact search the positions of query q, find the first exact or inexact bound
    vector<fmint_t>& inexSearchFindBest(vector<u_int8_t> &q,const int z)
    {
        best_bound[0]=0;
        best_bound[1]=0;
        best_bound_found=0;
        _findBestBound(q,q.begin(),z,0,len);
        findAlignments(best_bound[0],best_bound[1]);
        return matches;
    }


    int returnMatchStatus()
    {
        return match_status;
    }

    ///return transformed string
    vector<u_int8_t>& returnBWTString()
    {
        return encoded_bwt_str;
    }

    array<fmint_t,CHILDNUM+2>& returnOcc()
    {
        return occ;
    }

    vector< array<fmint_t,CHILDNUM+2> >& returnCheckpoints()
    {
        return C;
    }

    vector<fmint_t>& returnSA()
    {
        return sa;
    }

    size_t returnBWTStringLength()
    {
        return len;
    }

    size_t returnCheckpointsLength()
    {
        return C_len;
    }

    char returnEOC()
    {
        return this->EOC;
    }

    char returnEOS()
    {
        return this->EOS;
    }

    void resetInexStartStat()
    {
        inex_start_stat.resize(inex_start_stat.size(),0);
    }

    void printInexStartStat()
    {
        for(auto it=inex_start_stat.begin();it!=inex_start_stat.end();it++)
        {
            cout<<it-inex_start_stat.begin()<<": "<<*it<<endl;
        }
    }

private:


    /*----------------------------------------------------------
     * For building bwt and fmindex
     * ---------------------------------------------------------*/


    ///calculate the first occurance of a letter in sorted string s or say C[c]
    ///count the number of letters for each step and return list of the counts
    /// s is the bwt transformed string
    void _calcOccAndCheckpoints(vector<u_int8_t> &encoded_s)
    {
        /*
        tr1::unordered_map<u_int8_t,fmint_t> A;//letter count
        fmint_t i=0;
        for(vector<u_int8_t>::iterator encoded_c=encoded_s.begin();encoded_c!=encoded_s.end();++encoded_c,i++)
        {
            if(i%step==0)
                C.push_back(A);
            if(A.find(*encoded_c)==A.end())//dont find c in A
                A[*encoded_c]=1;
            else// find c in A
                A[*encoded_c]=A[*encoded_c]+1;
        }
        C.shrink_to_fit();
        vector<u_int8_t> letters;
        for(tr1::unordered_map<u_int8_t,fmint_t>::iterator it_A=A.begin();it_A!=A.end();++it_A)
            if(it_A->first<5)
                letters.push_back(it_A->first);
        sort(letters.begin(),letters.end());//sort letters, if not sorted by map/unorded_map EOS,EOC not included
        fmint_t idx=0;
        occ[encode(EOS)]=idx;
        idx=idx+A[encode(EOS)];
        occ[encode(EOC)]=idx;
        idx=idx+A[encode(EOC)];
        for(vector<u_int8_t>::iterator encoded_c=letters.begin();encoded_c!=letters.end();++encoded_c)
        {
            occ[*encoded_c]=idx;//first index of letters
            idx=idx+A[*encoded_c];
        }
        C_len=C.size();
        */
        array<fmint_t,CHILDNUM+2> A;
        A.fill(0);
        fmint_t i;
        C.reserve(len/step+1);
        for(i=0;i<len;i++)
        {
            if(i%step==0)
                C.push_back(A);
            A[encoded_s[i]]++;
        }

        //occ.resize(CHILDNUM+2);
        fmint_t idx=0;
        occ[CHILDNUM]=idx;
        idx+=A[CHILDNUM];
        occ[CHILDNUM+1]=idx;
        idx+=A[CHILDNUM+1];
        for(i=0;i<CHILDNUM;i++)
        {
            occ[i]=idx;//first index of letters
            idx+=A[i];
        }
        C_len=C.size();
    }


    ///remain only one entry for full_sa[] for each gap return checkpoints of full_sa
    void _calcCheckpointingSA(saint64_t *full_sa)
    {
        sa.reserve(len/gap+1); //get ceiling of len/gap
        for(size_t i=0;i<len;i+=gap)
            sa.push_back((fmint_t)full_sa[i]);
        delete [] full_sa;
    }



    /*----------------------------------------------------------
     * For query search
     * ---------------------------------------------------------*/

    /*
    ///find the first and last suffix positions for query q
    void _findBound(vector<u_int8_t> &q,fmint_t &top,fmint_t &bot)
    {
        vector<u_int8_t>::iterator qc;
        for(qc=q.begin();qc!=q.end();++qc)
        {
            top=_lf(top,*qc);
            bot=_lf(bot,*qc);
            if(top==bot)
            {
                //top=-1;
                //bot=-1;
                break;
            }
        }
    }*/

    ///find the first and last suffix positions for query q
    vector<u_int8_t>::iterator _findBound(vector<u_int8_t> &q,fmint_t &top,fmint_t &bot,uint pos)
    {
        fmint_t new_top,new_bot;
        vector<u_int8_t>::iterator qc;
        for(qc=q.begin();qc!=q.begin()+pos;++qc)
        {
            new_top=_lf(top,*qc);
            new_bot=_lf(bot,*qc);
            if(new_top<new_bot)
            {
                top=new_top;
                bot=new_bot;
            }
            else
            {
                //qc=qc-1;
                break;
            }
        }
        return qc;
    }
    
    ///find the first and last suffix positions for query q
    vector<u_int8_t>::iterator _findBound(vector<u_int8_t> &q,fmint_t &top,fmint_t &bot)
    {
        fmint_t new_top,new_bot;
        vector<u_int8_t>::iterator qc;
        for(qc=q.begin();qc!=q.end();++qc)
        {
            new_top=_lf(top,*qc);
            new_bot=_lf(bot,*qc);
            if(new_top<new_bot)
            {
                top=new_top;
                bot=new_bot;
            }
            else
            {
                //qc=qc-1;
                break;
            }
        }
        return qc;
    }

    ///find the first and last suffix positions for query q
    template <typename Iter>
    void _findInexBoundsOriginal(const vector<u_int8_t> &q,const Iter qc,const int z,const fmint_t top,const fmint_t bot)
    {
        if(z<0)//exceed number allowed of mismatches
        {
            return;
        }
        if(qc==q.end())//complete matching
        {
            //cout<<top<<" "<<bot<<endl;
            inex_bounds.push_back(make_pair(top,bot));
            return;
        }
        //_findInexBoundsOriginal(q,qc+1,z-1,top,bot);//skip one letter of q,insertion
        fmint_t new_top,new_bot;
        for(u_int8_t c=0;c<CHILDNUM;c++)
        {
            new_top=_lf(top,c);
            new_bot=_lf(bot,c);
            if(new_top<new_bot)
            {
                //_findInexBoundsOriginal(q,qc,z-1,new_top,new_bot);//add one letter for q,deletion
                if(c==*qc)
                {
                    _findInexBoundsOriginal(q,qc+1,z,new_top,new_bot);
                }
                else
                {
                    _findInexBoundsOriginal(q,qc+1,z-1,new_top,new_bot);
                }
            }
        }
    }

    ///find the first and last suffix positions for query q
    template <typename Iter>
    void _findInexBounds(const vector<u_int8_t> &q,Iter qc,const int z,const fmint_t top,const fmint_t bot)
    {
        if(z<0)//exceed number allowed of mismatches
        {
            return;
        }
        if(qc==q.end())//complete matching
        {
            //cout<<top<<" "<<bot<<endl;
            inex_bounds.push_back(make_pair(top,bot));
            return;
        }
        if(z)
        {
            //_findInexBounds(q,qc+1,z-2,top,bot);//skip one letter of q,insertion
            fmint_t new_top,new_bot;
            for(u_int8_t c=0;c<CHILDNUM;c++)
            {
                new_top=_lf(top,c);
                new_bot=_lf(bot,c);
                if(new_top<new_bot)
                {
                    //_findInexBounds(q,qc,z-2,new_top,new_bot);//add one letter for q,deletion
                    if(c==*qc)
                    {
                        _findInexBounds(q,qc+1,z,new_top,new_bot);
                    }
                    else
                    {
                        _findInexBounds(q,qc+1,z-1,new_top,new_bot);
                    }
                }
            }
        }
        else//if z==0
        {
            fmint_t new_top=top,new_bot=bot;
            for(;qc!=q.end();++qc)
            {
                new_top=_lf(new_top,*qc);
                new_bot=_lf(new_bot,*qc);
                if(new_top==new_bot)
                {   break;}
            }
            if(qc==q.end())
            {
                inex_bounds.push_back(make_pair(new_top,new_bot));
                //cout<<top<<" "<<bot<<endl;
            }
        }
    }
    
    ///find the first and last suffix positions for query q
    template <typename Iter>
    void _findInexBoundsOnceForAllLetters(const vector<u_int8_t> &q,Iter qc,const int z,fmint_t top,fmint_t bot,const int stack_level,bool IN=false,bool DEL=false)
    {
        if(z<0)//exceed number allowed of mismatches
        {
            return;
        }
        if(qc==q.end())//complete matching
        {
            inex_bounds.push_back(make_pair(top,bot));
            return;
        }
        if(z)
        {
            if(DEL)
            {
                _findInexBoundsOnceForAllLetters(q,qc+1,-1,top,bot,stack_level,true,false);//skip one letter of q,insertion
            }
            else
            {
                _findInexBoundsOnceForAllLetters(q,qc+1,z-1,top,bot,stack_level,true,false);//skip one letter of q,insertion
            }
            _lfAllLetters(top,stack_level,0);
            _lfAllLetters(bot,stack_level,1);
            for(u_int8_t c=0;c<CHILDNUM;c++)
            {
                top=bound_stack[entries_of_bound_stack*stack_level+c*2];
                bot=bound_stack[entries_of_bound_stack*stack_level+c*2+1];
                if(top<bot)
                {
                    if(IN)
                    {
                        _findInexBoundsOnceForAllLetters(q,qc,-1,top,bot,stack_level+1,false,true);//add one letter for q,deletion
                    }
                    else
                    {
                        _findInexBoundsOnceForAllLetters(q,qc,z-1,top,bot,stack_level+1,false,true);//add one letter for q,deletion
                    }
                    if(c==*qc)
                    {
                        _findInexBoundsOnceForAllLetters(q,qc+1,z,top,bot,stack_level+1,false,false);
                    }
                    else
                    {
                        _findInexBoundsOnceForAllLetters(q,qc+1,z-1,top,bot,stack_level+1,false,false);
                    }
                }
            }
            for(uint i=entries_of_bound_stack*stack_level;i<entries_of_bound_stack*stack_level+CHILDNUM*2;i++)
                bound_stack[i]=0;
        }
        else//if z==0
        {
            fmint_t new_top=top,new_bot=bot;
            for(;qc!=q.end();++qc)
            {
                new_top=_lf(new_top,*qc);
                new_bot=_lf(new_bot,*qc);
                if(new_top==new_bot)
                {   break;}
            }
            if(qc==q.end())
            {
                inex_bounds.push_back(make_pair(new_top,new_bot));
            }
        }
    }

    ///find the first bound met, whatever it is exact bound or inexact bound
    template <typename Iter>
    void _findBestBound(const vector<u_int8_t> &q,Iter qc,const int z,const fmint_t top,const fmint_t bot)
    {
        if(z<0)//exceed number allowed of mismatches
        {
            return;
        }
        if(qc==q.end())//complete matching
        {
            //cout<<top<<" "<<bot<<endl;
            best_bound[0]=top;
            best_bound[1]=bot;
            best_bound_found=1;
            return;
        }

        fmint_t new_top,new_bot;
        if(z)
        {
            new_top=_lf(top,*qc);
            new_bot=_lf(bot,*qc);
            if(new_top<new_bot)
            {
                _findBestBound(q,qc+1,z,new_top,new_bot);
            }
            for(u_int8_t i=1;i<4;i++)
            {
                if(best_bound_found)
                    break;
                new_top=_lf(top,((*qc)+i)%4);
                new_bot=_lf(bot,((*qc)+i)%4);
                if(new_top<new_bot)
                {
                    _findBestBound(q,qc+1,z-1,new_top,new_bot);
                }
            }
        }
        else//z==0
        {
            new_top=top;
            new_bot=bot;
            for(;qc!=q.end();++qc)
            {
                new_top=_lf(new_top,*qc);
                new_bot=_lf(new_bot,*qc);
                if(new_top==new_bot)
                {   break;  }
            }
            if(qc==q.end())
            {
                best_bound[0]=top;
                best_bound[1]=bot;
                best_bound_found=1;
                //cout<<top<<" "<<bot<<endl;
            }
        }
    }

    ///find the first bound met, whatever it is exact bound or inexact bound
    template <typename Iter>
    void _findBestBound(const vector<u_int8_t> &q,Iter qc,const int z,const fmint_t top,const fmint_t bot,bool IN,bool DEL)
    {
        if(z<0)//exceed number allowed of mismatches
        {
            return;
        }
        if(qc==q.end())//complete matching
        {
            //cout<<top<<" "<<bot<<endl;
            best_bound[0]=top;
            best_bound[1]=bot;
            best_bound_found=1;
            return;
        }

        fmint_t new_top,new_bot;
        if(z)
        {
            new_top=_lf(top,*qc);
            new_bot=_lf(bot,*qc);
            if(new_top<new_bot)
            {
                _findBestBound(q,qc+1,z,new_top,new_bot);
                if(IN)
                {
                    _findBestBound(q,qc,-1,new_top,new_bot,false,true);//add one letter for q,deletion
                }
                else
                {
                    _findBestBound(q,qc,z-1,new_top,new_bot,false,true);//add one letter for q,deletion
                }
            }
            for(u_int8_t i=1;i<4;i++)
            {
                if(best_bound_found)
                    break;
                new_top=_lf(top,((*qc)+i)%4);
                new_bot=_lf(bot,((*qc)+i)%4);
                if(new_top<new_bot)
                {
                    _findBestBound(q,qc+1,z-1,new_top,new_bot);
                    if(IN)
                    {
                        _findBestBound(q,qc,-1,new_top,new_bot,false,true);//add one letter for q,deletion
                    }
                    else
                    {
                        _findBestBound(q,qc,z-1,new_top,new_bot,false,true);//add one letter for q,deletion
                    }
                }
            }
            if(DEL)
            {
                _findBestBound(q,qc+1,-1,top,bot,true,false);//skip one letter of q,insertion
            }
            else
            {
                _findBestBound(q,qc+1,z-1,top,bot,true,false);//skip one letter of q,insertion
            }
        }
        else//z==0
        {
            new_top=top;
            new_bot=bot;
            for(;qc!=q.end();++qc)
            {
                new_top=_lf(new_top,*qc);
                new_bot=_lf(new_bot,*qc);
                if(new_top==new_bot)
                {   break;  }
            }
            if(qc==q.end())
            {
                best_bound[0]=top;
                best_bound[1]=bot;
                best_bound_found=1;
                //cout<<top<<" "<<bot<<endl;
            }
        }
    }


    void _lfAllLetters(const fmint_t idx,const int stack_level,const int tb)
    {
        //find the nearest checkpoint for idx
        fmint_t check=(idx+(step/2))/step;
        if(check>=C_len)
            check=C_len-1;
        fmint_t pos=check*step;

        //range between pos and idx
        fmint_t range[2];
        if(pos>idx)
        {
            range[0]=idx;
            range[1]=pos;
        }
        else
        {
            range[0]=pos;
            range[1]=idx;
        }

        //count of all children letters between pos and idx
        for(fmint_t i=range[0];i<range[1];i++)
        {
            //bound_stack[stack_level][A/T/C/G][top/bot]++
            bound_stack[entries_of_bound_stack*stack_level+encoded_bwt_str[i]*2+tb]+=1;
        }

        //calculate the letter count upto idx (not included)
        /*
        tr1::unordered_map<u_int8_t,fmint_t>::iterator Ccheck_it,occ_it;
        if(pos>idx)
        {
            for(u_int8_t j=0;j<CHILDNUM;j++)
            {
                Ccheck_it=C[check].find(j);
                if(Ccheck_it!=C[check].end())
                    bound_stack[entries_of_bound_stack*stack_level+j*2+tb]=Ccheck_it->second-bound_stack[entries_of_bound_stack*stack_level+j*2+tb];
                else
                    bound_stack[entries_of_bound_stack*stack_level+j*2+tb]=Ccheck_it->second;
                occ_it=occ.find(j);
                if(occ_it!=occ.end())
                    bound_stack[entries_of_bound_stack*stack_level+j*2+tb]+=occ_it->second;
            }
        }
        else
        {
            for(u_int8_t j=0;j<CHILDNUM;j++)
            {
                Ccheck_it=C[check].find(j);
                if(Ccheck_it!=C[check].end())
                    bound_stack[entries_of_bound_stack*stack_level+j*2+tb]+=Ccheck_it->second;
                occ_it=occ.find(j);
                if(occ_it!=occ.end())
                    bound_stack[entries_of_bound_stack*stack_level+j*2+tb]+=occ_it->second;
            }
        }*/


        if(pos>idx)
        {
            for(u_int8_t j=0;j<CHILDNUM;j++)
            {
                bound_stack[entries_of_bound_stack*stack_level+j*2+tb]=C[check][j]-bound_stack[entries_of_bound_stack*stack_level+j*2+tb]+occ[j];
            }
        }
        else
        {
            for(u_int8_t j=0;j<CHILDNUM;j++)
            {
                bound_stack[entries_of_bound_stack*stack_level+j*2+tb]+=C[check][j]+occ[j];
            }
        }

    }


    ///Count the number of a letter upto idx in s using checkpoints
    ///
    /// Arguments:
    /// idx     -- count upto this position
    /// letter  -- count for this letter
    /// C       -- is the list of checkpoints
    /// step    -- is the step of the checkpoints
    /// bwt_str -- the transformed string
    fmint_t _countLetterWithCheckpoints(const fmint_t idx,const u_int8_t letter)
    {
        //find the nearest checkpoint for idx
        fmint_t check=(idx+(step/2))/step;
        if(check>=C_len)
        {
            check=C_len-1;
        }
        fmint_t pos=check*step;

        //count of the letter s[idx] upto pos (not included)
        fmint_t count;/*
        tr1::unordered_map<u_int8_t,fmint_t>::iterator Ccheck_it= C[check].find(letter);
        if(Ccheck_it==C[check].end())//cannot find
        {
            count=0;
        }
        else
        {
            count=Ccheck_it->second;
        }*/
        count=C[check][letter];
        //range between pos and idx
        fmint_t range[2];
        if(pos<idx)
        {
            range[0]=pos;
            range[1]=idx;
        }
        else
        {
            range[0]=idx;
            range[1]=pos;
        }

        //count of letters between pos and idx
        fmint_t i,k=0;
        for(i=range[0];i<range[1];i++)
        {
            if(letter==encoded_bwt_str[i])
            {
                k=k+1;
            }
        }

        //calculate the letter count upto idx (not included)
        if(pos<idx)
        {
            count=count+k;
        }
        else
        {
            count=count-k;
        }

        return count;
    }

    ///find the offset in position idx of transformed string from the beginning
    fmint_t _walk(const fmint_t idx)
    {
        //walk to the a checkpoint using lf mapping
        fmint_t r=0,i=idx;
        //bool quit_early=false;
        while(i%gap)
        {
            r=r+1;
            i=_lf(i,encoded_bwt_str[i]);
        }
        //usage of sa for faster searches
        //if(!quit_early)
        //r=(sa[i/gap]+r)%len;
        return (sa[i/gap]+r)%len;
    }
    /*
    ///get the first occurance of letter qc in left-column
    fmint_t _occ(const u_int8_t qc)
    {
        fmint_t c;
        tr1::unordered_map<u_int8_t,fmint_t>::iterator occ_it=occ.find(qc);
        if(occ_it==occ.end())//cannot find qc
        {
            c=0;
        }
        else
        {
            c=occ_it->second;
        }
        return c;
    }*/

    ///count the occurances of letter qc (rank of qc) upto position idx
    fmint_t _count(const fmint_t idx,const u_int8_t qc)
    {
        //fmint_t count=_countLetterWithCheckpoints(idx,qc);
        //find the nearest checkpoint for idx
        fmint_t check=(idx+(step/2))/step;
        if(check>=C_len)
        {
            check=C_len-1;
        }
        fmint_t pos=check*step;

        //count of the letter s[idx] upto pos (not included)
        fmint_t count;/*
        tr1::unordered_map<u_int8_t,fmint_t>::iterator Ccheck_it= C[check].find(letter);
        if(Ccheck_it==C[check].end())//cannot find
        {
            count=0;
        }
        else
        {
            count=Ccheck_it->second;
        }*/
        count=C[check][qc];

        //range between pos and idx
        fmint_t range[2];
        if(pos<idx)
        {
            range[0]=pos;
            range[1]=idx;
        }
        else
        {
            range[0]=idx;
            range[1]=pos;
        }

        //count of letters between pos and idx
        fmint_t i,k=0;
        for(i=range[0];i<range[1];i++)
        {
            if(qc==encoded_bwt_str[i])
            {
                k=k+1;
            }
        }

        //calculate the letter count upto idx (not included)
        if(pos<idx)
        {
            count=count+k;
        }
        else
        {
            count=count-k;
        }

        return count;
    }


    ///get the nearest lf mapping for letter qc at position idx
    fmint_t _lf(const fmint_t idx,const u_int8_t qc)
    {
        //fmint_t o,c;
        //o=_occ(qc);
        //c=_count(idx,qc);
        //return o+c;
        //return _occ(qc)+_count(idx,qc);
        return occ[qc]+_countLetterWithCheckpoints(idx,qc);
    }

    ///encode string
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
            default: i=CHILDNUM+1;break;
        }

        if(i==CHILDNUM+1)
        {
            if(c==EOS)
                i=CHILDNUM;
        }
        return (u_int8_t)i;
    }
    /*----------------------------------------------------------
     * For serialization
     * ---------------------------------------------------------*/

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & step;
        ar & gap;
        ar & encoded_bwt_str;
        ar & sa;
        ar & occ;
        ar & C;
        ar & len;
        ar & C_len;
    }


    fmint_t step;//step of the checkpoints C
    fmint_t gap;//gap of the checkpoints sa
    char EOS;//end mark for whold string
    char EOC;//end mark for each chromosome substring
    vector<u_int8_t> encoded_bwt_str;//transformed bwt string
    saint64_t *full_sa;//suffix array
    vector<fmint_t> sa;//checkpoints of suffix array
    //tr1::unordered_map<u_int8_t,fmint_t> occ;//first occurance of letters
    array<fmint_t,CHILDNUM+2> occ;//first occurance of letters
    //vector< tr1::unordered_map<u_int8_t,fmint_t> > C;//list of checkpoints
    vector< array<fmint_t,CHILDNUM+2> > C;
    size_t len;//length of transformed bwt string
    size_t C_len;//length of checkpoints
    vector< pair<fmint_t,fmint_t> > inex_bounds;//container holding inexact matched top,bot pairs
    array<fmint_t,2> best_bound;
    array<fmint_t,(CHILDNUM+2)*2*MAXREADLEN> bound_stack; //collapse a 3d stack array into 1d bound_stack[level][letter][top/bot]
    uint num_alignments_output;
    int entries_of_bound_stack;
    int best_bound_found;
    friend class boost::serialization::access;//serialization

public:
    int match_status;//indicate a query is exactly matched or inexactly matched, 1 = exact, 0 = inexact
    vector<fmint_t> matches;//query results
    vector<size_t> inex_start_stat;

};

#endif // BWT_H

