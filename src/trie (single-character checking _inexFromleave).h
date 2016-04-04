///Author:Yujia Wu
///Created on Feb 2 2015
///


#ifndef TRIE_H
#define TRIE_H



#include <vector>
#include <string>
#include <time.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <stack>
#include <math.h>
#include <utility>

#include "utils.h"
#include "bwt.h"





class BoundInfo
{
public:
    BoundInfo(fmint_t start_pos,fmint_t top,fmint_t bot)
    {
        this->start_pos=start_pos;
        this->top=top;
        this->bot=bot;
    }

    fmint_t start_pos,top,bot;

};

class CNode
{

public:

    CNode(const vector<u_int8_t> &seq)
    {
        key=seq;
        child=NULL;
        sibling=NULL;
        nameList=NULL;
        junc_flag=false;
    }

    template <typename Iter>
    CNode(Iter it,Iter end)
    {
        key.assign(it,end);
        child=NULL;
        sibling=NULL;
        nameList=NULL;
        junc_flag=false;
    }

    CNode(const vector<u_int8_t> &seq,const string &name)
    {
        key=seq;
        child=NULL;
        sibling=NULL;
        junc_flag=false;
        nameList=new vector<string>;
        nameList->push_back(std::move(name));
    }

    template <typename Iter>
    CNode(Iter it,Iter end,const string &name)
    {
        key.assign(it,end);
        child=NULL;
        sibling=NULL;
        junc_flag=false;
        nameList=new vector<string>;
        nameList->push_back(std::move(name));
    }

    ~CNode()
    {
        if(child)
            delete child;
        if(sibling)
            delete sibling;
        if(nameList)
            delete nameList;
    }

    void insert(const vector<u_int8_t> &seq,const string &name,CNode* parent)
    {
        u_int16_t k=this->commPrefix(seq);
        if(k==0)
        {
            if(!this->sibling)
            {
                this->sibling=new CNode(seq,name);
                //if a new sibling is created, then its parent node must be a junction
                if(parent)
                    parent->junc_flag=true;
                return;
            }
            this->sibling->insert(seq,name,parent);
        }
        else if(k<seq.size())
        {
            if(k<key.size())
                this->split(k);
            vector<u_int8_t> sub_seq(seq.begin()+k,seq.end());
            if(!this->child)
            {
                this->child=new CNode(sub_seq,name);
                return;
            }
            this->child->insert(sub_seq,name,this);
        }
        //deal with repeated sequence
        else if(k==seq.size())
        {
            if(k<key.size()) //e.g. insert 'meet' when 'meeting' exists
                this->split(k);
            else if(k==key.size()) //the inserting read has already existed e.g. insert 'meet' when 'meet' exists
            {
                if(!this->nameList)
                    this->nameList=new vector<string>;
                this->nameList->push_back(std::move(name));
            }

        }
    }

    template <typename Iter>
    void insertByIter(Iter it,Iter begin,Iter end,const string &name, CNode* parent)
    {
        this->commPrefixByIter(it,end);
        if(it==begin)
        {
            if(!this->sibling)
            {
                this->sibling=new CNode(begin,end,name);
                //if a new sibling is created, then its parent node must be a junction
                if(parent)
                {
                    parent->junc_flag=true;
                }
                return;
            }
            this->sibling->insertByIter(begin,begin,end,name,parent);
        }
        else if(it!=end)
        {
            if((it-begin)<key.size())
            {
                Iter key_it=key.begin()+(it-begin);
                this->splitByIter(key_it);
            }
            //vector<u_int8_t> sub_seq(begin+k,seq.end());
            if(!this->child)
            {
                this->child=new CNode(it,end,name);
                return;
            }
            this->child->insertByIter(it,it,end,name,this);
        }
        //deal with repeated sequence
        else if(it==end)
        {
            if((it-begin)<key.size()) //e.g. insert 'meet' when 'meeting' exists
            {
                Iter key_it=key.begin()+(it-begin);
                this->splitByIter(key_it);
            }
            else if((it-begin)==key.size()) //the inserting read has already existed e.g. insert 'meet' when 'meet' exists
            {
                if(!this->nameList)
                    this->nameList=new vector<string>;
                this->nameList->push_back(std::move(name));
            }
        }
    }

private:

    u_int16_t commPrefix(const vector<u_int8_t> &seq)
    {
        for(size_t i=0;i<seq.size();i++)
            if(i==key.size() || seq[i]!=key[i])
                return i;
        return seq.size();
    }

    template <typename Iter>
    void commPrefixByIter(Iter &seq_it,Iter seq_end)
    {
        Iter key_it=key.begin();
        for(;seq_it!=seq_end;seq_it++,key_it++)
            if(key_it==key.end() || *seq_it!=*key_it)
                return;
    }

    void split(u_int16_t k)
    {
        vector<u_int8_t> s(this->key.begin()+k,this->key.end());
        CNode* newNode=new CNode(s);
        newNode->child=this->child;
        this->child=newNode;
        newNode->nameList=this->nameList;
        this->nameList=NULL;
        newNode->junc_flag=this->junc_flag;
        this->junc_flag=false;
        vector<u_int8_t> sub_key(this->key.begin(),this->key.begin()+k);
        this->key=sub_key;
    }

    template <typename Iter>
    void splitByIter(Iter key_it)
    {
        CNode* newNode=new CNode(key_it,this->key.end());
        newNode->child=this->child;
        this->child=newNode;
        newNode->nameList=this->nameList;
        this->nameList=NULL;
        newNode->junc_flag=this->junc_flag;
        this->junc_flag=false;
        //this->key.assign(this->key.begin(),key_it);
        this->key.resize(key_it-this->key.begin());
    }

public:
    vector<string>* nameList;
    vector<u_int8_t> key;
    CNode* child;
    CNode* sibling;
    bool junc_flag;

};



class CTrie
{
public:
    CTrie()
    {
        root=NULL;
        empty_flag=false;
        bound_stack.fill(0);
        stack_level=-1;
        seq_stack.reserve(MAXREADLEN);
        seq_out=new char[MAXREADLEN];
        resetStatistics();
    }
    ~CTrie(){delete root;}

    void resetStatistics()
    {
        count_matched=0;
        count_exact_matched=0;
        count_unmatched=0;
        count_node=0;
        count_junc_node=0;
        count_repeated_node=0;
        count_path=0;
        count_node_on_this_path=0;
        count_path_node=0;
        count_comparison=0;
        two_nodes_junc_count=0;
        three_nodes_junc_count=0;
        four_nodes_junc_count=0;
        five_nodes_junc_count=0;
        more_nodes_junc_count=0;
        unqiue_suffix_len=0;
        leaf_node_len=0;
        unique_shortest_suffix_len=MAXREADLEN;
        unique_longest_suffix_len=0;
    }
    void pushBound(const fmint_t start_pos,const fmint_t top,const fmint_t bot)
    {
        stack_level+=1;
        bound_stack[3*stack_level]=start_pos;
        bound_stack[3*stack_level+1]=top;
        bound_stack[3*stack_level+2]=bot;
    }
    inline void popBound()
    {
        //bound_stack[3*stack_level]=0;
        //bound_stack[3*stack_level+1]=0;
        //bound_stack[3*stack_level+2]=0;
        stack_level-=1;
    }

    //point to a fmindex
    void specifyIndex(FMIndex *fmidx,ChromIndex *chridx)
    {
        this->fmindex=fmidx;
        this->chrindex=chridx;
        //BoundInfo bound(0,0,(fmint_t)fmindex->returnBWTString().size());
        //bound_stack.push(bound);
        pushBound(0,0,(fmint_t)fmindex->returnBWTString().size());
    }

    ///output statistics info
    void summary()
    {
        float avgPathLen=(float)count_path_node/count_path;
        size_t num_bwt_array_scanned_in_junction_by_ITFM=two_nodes_junc_count+three_nodes_junc_count+four_nodes_junc_count+five_nodes_junc_count;
        size_t num_bwt_array_scanned_in_junction_by_TFM=two_nodes_junc_count*2+three_nodes_junc_count*3+four_nodes_junc_count*4+five_nodes_junc_count*5;
        cout<<"# of nodes: "<<count_node<<endl;
        cout<<"# of junction nodes: "<<count_junc_node<<endl;
        cout<<"# of repeated reads: "<<count_repeated_node<<endl;
        cout<<"# of path: "<<count_path<<endl;
        cout<<"average path length: "<<avgPathLen<<" nodes/path"<<endl;
        cout<<"# of junction nodes that have two children: "<<two_nodes_junc_count<<endl;
        cout<<"# of junction nodes that have three children: "<<three_nodes_junc_count<<endl;
        cout<<"# of junction nodes that have four children: "<<four_nodes_junc_count<<endl;
        cout<<"# of junction nodes that have five children: "<<five_nodes_junc_count<<endl;
        cout<<"# of junction nodes that have more children: "<<more_nodes_junc_count<<endl;
        //cout<<"# of invoking bwt searching by TFM: "<<count_comparison<<endl;
        //cout<<"# of invoking bwt searching by ITFM: "<<count_comparison-(num_bwt_array_scanned_in_junction_by_TFM-num_bwt_array_scanned_in_junction_by_ITFM)<<endl;
        //cout<<"# of invoking bwt searching in junction by TFM: "<<num_bwt_array_scanned_in_junction_by_TFM<<endl;
        //cout<<"# of invoking bwt searching in junction by ITFM: "<<num_bwt_array_scanned_in_junction_by_ITFM<<endl;
        //cout<<"# of saved bwt array scan by ITFM: "<<num_bwt_array_scanned_in_junction_by_TFM-num_bwt_array_scanned_in_junction_by_ITFM<<endl;
        //cout<<"performance increase TFM vs ITFM: "<<(float)((float)count_comparison/(float)(count_comparison-(num_bwt_array_scanned_in_junction_by_TFM-num_bwt_array_scanned_in_junction_by_ITFM)))<<endl;
        cout<<"total unqiue suffix length:"<<unqiue_suffix_len<<endl;
        //cout<<"total leaf nodes length: "<<leaf_node_len<<endl;
        //cout<<"average leaf node length: "<<(float)leaf_node_len/(float)(count_node-count_junc_node)<<endl;
        cout<<"shortest unqiue suffix length: "<<unique_shortest_suffix_len<<endl;
        cout<<"longest unqiue suffix length: "<<unique_longest_suffix_len<<endl;

    }
/*
    void appendRead(const string &seq,const string name)
    {
        if(!root)
        {
            root=new CNode(seq,name);
            return;
        }
        root->insert(seq,name,NULL);
    }count_comparison
*/
    void appendRead(vector<u_int8_t> &seq,const string name)
    {
        vector<u_int8_t>::iterator begin=seq.begin(),end=seq.end();
        if(!root)
        {
            root=new CNode(begin,end,name);
            return;
        }
        root->insertByIter(begin,begin,end,name,NULL);
    }

    template <typename Iter>
    void appendRead(Iter begin,Iter end,const string name)
    {
        if(!root)
        {
            root=new CNode(begin,end,name);
            return;
        }
        root->insertByIter(begin,begin,end,name,NULL);
    }

    inline void traverse()
    {
        _traverse(root);
    }

    void run(string out_path)
    {
        string unmatched_out_path;
        unmatched_out_path.append("unmatched.bin");
        ctrie_out.open(out_path,ios::out|ios::app);
        ctrie_unmatched_out.open(unmatched_out_path,ios::out|ios::app|ios::binary);
        _run(root);
        ctrie_out.close();
        ctrie_unmatched_out.close();
        cout<<"count_exact_matched: "<<count_exact_matched<<endl;
        cout<<"count_unmatched: "<<count_unmatched<<endl;
    }

    void inexRun(string out_path,int mismatch_allowed)
    {
        this->mismatch_allowed=mismatch_allowed;
        string unmatched_out_path;
        unmatched_out_path.append("unmatched.txt");
        ctrie_out.open(out_path,ios::out|ios::app);
        ctrie_unmatched_out.open(unmatched_out_path,ios::out|ios::app|ios::binary);
        _inexRun(root);
        ctrie_out.close();
        ctrie_unmatched_out.close();
        cout<<"count_matched: "<<count_matched<<endl;
        cout<<"count_exact_matched: "<<count_exact_matched<<endl;
        cout<<"count_INexact_matched: "<<(int)count_matched-(int)count_exact_matched<<endl;
        cout<<"count_unmatched: "<<count_unmatched<<endl;
        cout<<"mapping rate: "<<(float)count_matched/(float)(count_unmatched+count_matched)<<endl;
    }

private:

    void _traverse(CNode* pnode)
    {
        if(!pnode)
            return;

        count_node++;
        count_node_on_this_path++;
        size_t key_len=pnode->key.size();
        for(size_t i=0;i<key_len;i++)
            seq_stack.push_back(pnode->key[i]);
        count_comparison+=key_len;
        _traverse(pnode->child);
        //run the following when pnode->child is null
        if(pnode->junc_flag)
        {
            count_junc_node++;
            CNode* childNode=pnode->child;
            int i=1;
            while((childNode=childNode->sibling))
            {
               i++;
            }
            switch(i)
            {
               case 2: two_nodes_junc_count++; break;
               case 3: three_nodes_junc_count++; break;
               case 4: four_nodes_junc_count++; break;
               case 5: five_nodes_junc_count++; break;
               default: more_nodes_junc_count++; break;
            }
        }
        vector<string>* nameList=pnode->nameList;
        if(nameList)
        {
            size_t size=nameList->size();
            count_repeated_node=count_repeated_node+size-1;
            if(!pnode->child)
            {
               count_path++;
               count_path_node+=count_node_on_this_path;
               leaf_node_len+=pnode->key.size();
               if(size==1)
               {
                   if(pnode->key.size()<unique_shortest_suffix_len)
                        unique_shortest_suffix_len=pnode->key.size();
                   if(pnode->key.size()>unique_longest_suffix_len)
                       unique_longest_suffix_len=pnode->key.size();
                   unqiue_suffix_len+=pnode->key.size();
               }
            }
        }
        for(size_t i=0;i<key_len;i++)
            seq_stack.pop_back();
        _traverse(pnode->sibling);

        count_node_on_this_path--;
    }

    void _run(CNode* pnode)
    {
        if(!pnode)
            return;
        size_t i,key_len=pnode->key.size();
        for(i=0;i<key_len;i++)
            seq_stack.push_back(pnode->key[i]);
        if(pnode->junc_flag)//junction node
        {
            fmint_t seq_start=bound_stack[3*stack_level];
            fmint_t top=bound_stack[3*stack_level+1];
            fmint_t bot=bound_stack[3*stack_level+2];
            fmindex->findBound(seq_stack,seq_start,top,bot);
            if(top==bot)
            {
                empty_flag=true;
            }
            else
            {
                pushBound(seq_stack.size(),top,bot);
            }
        }
        vector<string>* nameList=pnode->nameList;
        if(nameList)
        {
            size_t size=nameList->size();
            fmint_t seq_start=bound_stack[3*stack_level];
            fmint_t top=bound_stack[3*stack_level+1],bot=bound_stack[3*stack_level+2];
            size_t seq_len=seq_stack.size();
            size_t j;
            //if(seq_start<(fmint_t)seq_stack.size())//if this node is only an end not a junction
            if(!pnode->junc_flag)
            {
                fmindex->findBound(seq_stack,seq_start,top,bot);
                fmindex->findAlignments(top,bot);
                if(fmindex->matches.size())
                {
                    /*for(j=0;j<seq_len;j++)
                        seq_out[j]=decode(seq_stack[j]);*/
                    count_exact_matched+=size;
                    for(i=0;i<size;i++)
                    {
                        ctrie_out<<(*nameList)[i]<<" ";
                        for(j=0;j<seq_len;j++)
                        {
                            ctrie_out<<seq_stack[j];
                        }
                        ctrie_out<<" "<<fmindex->matches.size()<<"  ";
                        for(j=0;j<fmindex->matches.size();j++)
                        {
                            RealPosition rpos=chrindex->calRealPos(fmindex->matches[j],seq_len);
                            ctrie_out<<rpos.name<<","<<rpos.pos<<" ";
                        }
                        ctrie_out<<endl;
                        //memset(seq_out,0,seq_len);
                    }
                }
                else
                {
                    count_unmatched+=size;
                    ctrie_unmatched_out.write((const char*)&seq_stack[0], seq_len);
                    ctrie_unmatched_out<<" "<<size<<endl;
                }
            }
            else
            {
                if(empty_flag)
                {
                    count_unmatched+=size;
                    ctrie_unmatched_out.write((const char*)&seq_stack[0], seq_len);
                    ctrie_unmatched_out<<" "<<size<<endl;
                }
                else
                {
                    fmindex->findAlignments(top,bot);
                    /*for(j=0;j<seq_len;j++)
                        seq_out[j]=decode(seq_stack[j]);*/
                    count_exact_matched+=size;
                    for(i=0;i<size;i++)
                    {
                        ctrie_out<<(*nameList)[i]<<" ";
                        for(j=0;j<seq_len;j++)
                        {
                            ctrie_out<<seq_stack[j];
                        }
                        ctrie_out<<" "<<fmindex->matches.size()<<"  ";
                        for(j=0;j<fmindex->matches.size();j++)
                        {
                            RealPosition rpos=chrindex->calRealPos(fmindex->matches[j],seq_len);
                            ctrie_out<<rpos.name<<","<<rpos.pos<<" ";
                        }
                        ctrie_out<<endl;
                        //memset(seq_out,0,seq_len);
                    }
                }
            }
        }
        if(empty_flag)
        {
            _run1(pnode->child);
            empty_flag=false;
        }
        else
        {
            _run(pnode->child);
            //if(pnode->junc_flag)
            //    bound_stack.pop();
            if(pnode->junc_flag)
                  popBound();
        }

        //run the following when pnode->child is null such as popping stack
        for(i=0;i<key_len;i++)
            seq_stack.pop_back();

        _run(pnode->sibling);
    }

    ///exact matching terminated earlier for mismatch in non leaf node
    void _run1(CNode* pnode)
    {
        if(!pnode)
            return;
        size_t i,key_len=pnode->key.size();
        for(i=0;i<key_len;i++)
            seq_stack.push_back(pnode->key[i]);
        vector<string>* nameList=pnode->nameList;
        if(nameList)
        {
            size_t size=nameList->size();
            count_unmatched+=size;
            ctrie_unmatched_out.write((const char*)&seq_stack[0], seq_stack.size());
            ctrie_unmatched_out<<" "<<size<<endl;
            /*
            for(i=0;i<size;i++)
            {
                ctrie_out<<(*nameList)[i]<<" ";
                for(size_t j=0;j<seq_stack.size();j++)
                {
                    ctrie_out<<seq_stack[j];
                }
                ctrie_out<<" "<<0<<endl;
            }*/
        }
        _run1(pnode->child);
        //run the following when pnode->child is null such as popping stack
        for(i=0;i<key_len;i++)
            seq_stack.pop_back();

        _run1(pnode->sibling);
    }


    void _inexRun(CNode* pnode)
    {
        if(!pnode)
            return;
        size_t i,key_len=pnode->key.size();
        for(i=0;i<key_len;i++)
            seq_stack.push_back(pnode->key[i]);
        fmint_t seq_start,top,bot;
        if(pnode->junc_flag)//junction node
        {
            seq_start=bound_stack[3*stack_level];
            top=bound_stack[3*stack_level+1];
            bot=bound_stack[3*stack_level+2];
            fmindex->findBound(seq_stack,seq_start,top,bot);
            if(top==bot)
            {
                empty_flag=true;
            }
            else
            {
                pushBound(seq_stack.size(),top,bot);
            }
        }
        vector<string>* nameList=pnode->nameList;
        if(nameList)
        {
            size_t size=nameList->size();
            seq_start=bound_stack[3*stack_level];
            top=bound_stack[3*stack_level+1];
            bot=bound_stack[3*stack_level+2];
            size_t seq_len=seq_stack.size();
            size_t j;
            //if(seq_start<(fmint_t)seq_stack.size())//if this node is only an end not a junction
            if(!pnode->junc_flag)
            {
                //string sub_seq(seq_stack.begin()+seq_start,seq_stack.end());
                fmindex->trieInexSearchFromLeaveFindBest(seq_stack,seq_start,top,bot,mismatch_allowed);
                if(fmindex->matches.size())
                {
                    /*for(j=0;j<seq_len;j++)
                        seq_out[j]=decode(seq_stack[j]);*/
                    if(fmindex->returnMatchStatus())
                        count_exact_matched+=size;
                    count_matched+=size;
                    for(i=0;i<size;i++)
                    {
                        ctrie_out<<(*nameList)[i]<<" ";
                        for(j=0;j<seq_len;j++)
                        {
                            ctrie_out<<seq_stack[j];
                        }
                        ctrie_out<<" "<<fmindex->matches.size()<<"  ";
                        for(j=0;j<fmindex->matches.size();j++)
                        {
                            RealPosition rpos=chrindex->calRealPos(fmindex->matches[j],seq_len);
                            ctrie_out<<rpos.name<<","<<rpos.pos<<" ";
                        }
                        ctrie_out<<endl;
                        //memset(seq_out,0,seq_len);
                    }
                }
                else
                {
                    count_unmatched+=size;
                    ctrie_unmatched_out.write((const char*)&seq_stack[0], seq_len);
                    ctrie_unmatched_out<<" "<<size<<endl;
                }
            }
            else
            {
                if(empty_flag)
                 {
                     count_unmatched+=size;
                     ctrie_unmatched_out.write((const char*)&seq_stack[0], seq_len);
                     ctrie_unmatched_out<<" "<<size<<endl;
                 }
                else
                {
                    fmindex->findAlignments(top,bot);
                    /*for(j=0;j<seq_len;j++)
                        seq_out[j]=decode(seq_stack[j]);*/
                    count_exact_matched+=size;
                    count_matched+=size;
                    for(i=0;i<size;i++)
                    {
                        ctrie_out<<(*nameList)[i]<<" ";
                        for(j=0;j<seq_len;j++)
                        {
                            ctrie_out<<seq_stack[j];
                        }
                        ctrie_out<<" "<<fmindex->matches.size()<<"  ";
                        for(j=0;j<fmindex->matches.size();j++)
                        {
                            RealPosition rpos=chrindex->calRealPos(fmindex->matches[j],seq_len);
                            ctrie_out<<rpos.name<<","<<rpos.pos<<" ";
                        }
                        ctrie_out<<endl;
                        //memset(seq_out,0,seq_len);
                    }
                }

            }
        }
        if(empty_flag)
        {
            _inexRun1(pnode->child);
            empty_flag=false;
        }
        else
        {
            _inexRun(pnode->child);
            //if(pnode->junc_flag)
            //    bound_stack.pop();
            if(pnode->junc_flag)
                  popBound();
        }
        //run the following when pnode->child is null such as popping stack
        for(i=0;i<key_len;i++)
            seq_stack.pop_back();

        _inexRun(pnode->sibling);
    }

    ///exact matching terminated earlier for mismatch in non leaf node
    void _inexRun1(CNode* pnode)
    {
        if(!pnode)
            return;
        size_t i,key_len=pnode->key.size();
        for(i=0;i<key_len;i++)
            seq_stack.push_back(pnode->key[i]);
        vector<string>* nameList=pnode->nameList;
        if(nameList)
        {
            size_t size=nameList->size();
            fmint_t seq_start=bound_stack[3*stack_level];
            fmint_t top=bound_stack[3*stack_level+1];
            fmint_t bot=bound_stack[3*stack_level+2];
            size_t j,seq_len=seq_stack.size();
            fmindex->trieInexSearchFromLeaveFindBest(seq_stack,seq_start,top,bot,mismatch_allowed);
            if(fmindex->matches.size())
            {
               /*for(j=0;j<seq_len;j++)
                   seq_out[j]=decode(seq_stack[j]);*/
               if(fmindex->returnMatchStatus())
                   count_exact_matched+=size;
               count_matched+=size;
               for(i=0;i<size;i++)
               {
                   ctrie_out<<(*nameList)[i]<<" ";
                   for(j=0;j<seq_len;j++)
                   {
                       ctrie_out<<seq_stack[j];
                   }
                   ctrie_out<<" "<<fmindex->matches.size()<<"  ";
                   for(j=0;j<fmindex->matches.size();j++)
                   {
                       RealPosition rpos=chrindex->calRealPos(fmindex->matches[j],seq_len);
                       ctrie_out<<rpos.name<<","<<rpos.pos<<" ";
                   }
                   ctrie_out<<endl;
                   //memset(seq_out,0,seq_len);
               }
           }
           else
           {
               count_unmatched+=size;
               ctrie_unmatched_out.write((const char*)&seq_stack[0], seq_len);
               ctrie_unmatched_out<<" "<<size<<endl;
           }
        }
        _inexRun1(pnode->child);
        //run the following when pnode->child is null such as popping stack
        for(i=0;i<key_len;i++)
            seq_stack.pop_back();

        _inexRun1(pnode->sibling);
    }



    void _run2(CNode* pnode)
    {
        if(!pnode)
            return;
        size_t i,key_len=pnode->key.size();
        for(i=0;i<key_len;i++)
            seq_stack.push_back(pnode->key[i]);
        if(pnode->junc_flag)//junction node
        {/*
            BoundInfo top_bound=bound_stack.top();
            fmint_t seq_start=top_bound.start_pos;
            fmint_t top=top_bound.top;
            fmint_t bot=top_bound.bot;
            //string sub_seq(seq_stack.begin()+seq_start,seq_stack.end());
            //fmindex->findBound(sub_seq,top,bot);
            fmindex->findBound(seq_stack,seq_start,top,bot);
            BoundInfo new_bound(seq_stack.size(),top,bot);
            bound_stack.push(new_bound);
            */
            fmint_t seq_start=bound_stack[3*stack_level];
            fmint_t top=bound_stack[3*stack_level+1];
            fmint_t bot=bound_stack[3*stack_level+2];
            fmindex->findBound(seq_stack,seq_start,top,bot);
            pushBound(seq_stack.size(),top,bot);
        }
        vector<string>* nameList=pnode->nameList;
        if(nameList)
        {
            size_t size=nameList->size();/*
            fmint_t seq_start=bound_stack.top().start_pos;
            fmint_t top=bound_stack.top().top,bot=bound_stack.top().bot;
            */
            fmint_t seq_start=bound_stack[3*stack_level];
            fmint_t top=bound_stack[3*stack_level+1],bot=bound_stack[3*stack_level+2];

            //if(seq_start<(fmint_t)seq_stack.size())//if this node is only an end not a junction
            if(!pnode->junc_flag)
            {
                //string sub_seq(seq_stack.begin()+seq_start,seq_stack.end());
                fmindex->findBound(seq_stack,seq_start,top,bot);
            }
            fmindex->findAlignments(top,bot);
            size_t seq_len=seq_stack.size();
            if(fmindex->matches.size())
            {
                size_t j;
                /*for(j=0;j<seq_len;j++)
                    seq_out[j]=decode(seq_stack[j]);*/
                count_exact_matched+=size;
                for(i=0;i<size;i++)
                {
                    ctrie_out<<(*nameList)[i]<<" ";
                    for(j=0;j<seq_len;j++)
                    {
                        ctrie_out<<seq_stack[j];
                    }
                    ctrie_out<<" "<<fmindex->matches.size()<<"  ";
                    for(j=0;j<fmindex->matches.size();j++)
                    {
                        RealPosition rpos=chrindex->calRealPos(fmindex->matches[j],seq_len);
                        ctrie_out<<rpos.name<<","<<rpos.pos<<" ";
                    }
                    ctrie_out<<endl;
                    //memset(seq_out,0,seq_len);
                }
            }
            else
            {
                count_unmatched+=size;
                ctrie_unmatched_out.write((const char*)&seq_stack[0], seq_len);
                ctrie_unmatched_out<<" "<<size<<endl;
            }

        }
        _run(pnode->child);
        //run the following when pnode->child is null such as popping stack
        for(i=0;i<key_len;i++)
            seq_stack.pop_back();
        //if(pnode->junc_flag)
        //    bound_stack.pop();
        if(pnode->junc_flag)
              popBound();
        _run(pnode->sibling);
    }

    void _inexRun2(CNode* pnode)
    {
        if(!pnode)
            return;
        size_t i,key_len=pnode->key.size();
        for(i=0;i<key_len;i++)
            seq_stack.push_back(pnode->key[i]);
        fmint_t seq_start,top,bot;
        if(pnode->junc_flag)//junction node
        {

            seq_start=bound_stack[3*stack_level];
            top=bound_stack[3*stack_level+1];
            bot=bound_stack[3*stack_level+2];
            fmindex->findBound(seq_stack,seq_start,top,bot);
            pushBound(seq_stack.size(),top,bot);
        }
        vector<string>* nameList=pnode->nameList;
        if(nameList)
        {
            size_t size=nameList->size();/*
            fmint_t seq_start=bound_stack.top().start_pos;
            fmint_t top=bound_stack.top().top,bot=bound_stack.top().bot;*/

            seq_start=bound_stack[3*stack_level];
            top=bound_stack[3*stack_level+1];
            bot=bound_stack[3*stack_level+2];
            //if(seq_start<(fmint_t)seq_stack.size())//if this node is only an end not a junction
            if(!pnode->junc_flag)
            {
                //string sub_seq(seq_stack.begin()+seq_start,seq_stack.end());
                fmindex->trieInexSearchFromLeaveFindBest(seq_stack,seq_start,top,bot,mismatch_allowed);
            }
            else
            {
                fmindex->findAlignments(top,bot);
            }
            //vector<fmint_t> matches=fmindex->matches;
            size_t seq_len=seq_stack.size();
            if(fmindex->matches.size())
            {
                size_t j;
                /*for(j=0;j<seq_len;j++)
                    seq_out[j]=decode(seq_stack[j]);*/
                if(fmindex->returnMatchStatus())
                    count_exact_matched+=size;
                count_matched+=size;
                for(i=0;i<size;i++)
                {
                    ctrie_out<<(*nameList)[i]<<" ";
                    for(j=0;j<seq_len;j++)
                    {
                        ctrie_out<<seq_stack[j];
                    }
                    ctrie_out<<" "<<fmindex->matches.size()<<"  ";
                    for(j=0;j<fmindex->matches.size();j++)
                    {
                        RealPosition rpos=chrindex->calRealPos(fmindex->matches[j],seq_len);
                        ctrie_out<<rpos.name<<","<<rpos.pos<<" ";
                    }
                    ctrie_out<<endl;
                    //memset(seq_out,0,seq_len);
                }

            }
            else
            {
                count_unmatched+=size;
                ctrie_unmatched_out.write((const char*)&seq_stack[0], seq_len);
                ctrie_unmatched_out<<" "<<size<<endl;
            }

        }
        _inexRun(pnode->child);
        //run the following when pnode->child is null such as popping stack
        for(i=0;i<key_len;i++)
            seq_stack.pop_back();
        //if(pnode->junc_flag)
        //    bound_stack.pop();
        if(pnode->junc_flag)
              popBound();

        _inexRun(pnode->sibling);
    }

    CNode* root;
    vector<u_int8_t> seq_stack;
    char *seq_out;
    bool empty_flag;
    //stack<BoundInfo> bound_stack;
    array<fmint_t,MAXREADLEN*3> bound_stack;//bound_stack[3*i+0]=start_pos,bound_stack[3*i+1]=top,bound_stack[3*i+2]=bot
    int stack_level;
    FMIndex* fmindex;
    ChromIndex* chrindex;
    int mismatch_allowed;
    ofstream ctrie_out,ctrie_unmatched_out;
    size_t count_exact_matched,count_matched,count_unmatched;
    size_t count_node,count_junc_node,count_repeated_node,
            count_path,count_node_on_this_path,count_path_node,
            two_nodes_junc_count,three_nodes_junc_count,
            four_nodes_junc_count,five_nodes_junc_count,
            more_nodes_junc_count,count_comparison,unqiue_suffix_len,
            unique_shortest_suffix_len,unique_longest_suffix_len,
            leaf_node_len;
};




#endif // TRIE_H
