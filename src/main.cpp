///Author:Yujia Wu
///Created on Dec 20 2014


#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <vector>
#include <array>
#include <string>
#include <tr1/unordered_map>
//#include <boost/archive/binary_oarchive.hpp>
//#include <boost/archive/binary_iarchive.hpp>

#include <zlib.h>
#include "kseq.h"
#include "trie.h"



KSEQ_INIT(gzFile, gzread)
using namespace std;



//int main(int argc,char *argv[])
int main()
{
    struct timeval start,finish;
    int l,count=0;
    double duration;
    fmint_t bwt_str_size=0,sa_size=0,occ_size=0,C_size=0;


    //output files
    string trie_out="ctrie_out.txt";
    string direct_out="out.txt";



    //genome
    string fasta_file="../test_data/genome/Cyanidioschyzon_merolae.ASM9120v1.28.dna.toplevel.fa";

    //reads
    string fastq_file="../test_data/reads/5M_simulated_35_long_reads_on_Cmerolae.fq";





    /*----------------------------------------------------------
     * simple test
     * ---------------------------------------------------------*/



/*
    //string str="GTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTACATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGTATGGCATCGCTGAAATGCAGCGGGCGAGTAGGCATCGCTGAAATGCAGCGGGCGAGTA";
    string str="AATG%ATGTTGGG";
    l=str.length();
    cout<<"length: "<<l<<endl;
    gettimeofday(&start,NULL);
    FMIndex fmindex;
    fmindex.transform(&str);//pass pointer of s
    fmindex.buildIndex();
    string bw=fmindex.returnBWTString();
    gettimeofday(&finish,NULL);
    duration=finish.tv_sec-start.tv_sec+(finish.tv_usec-start.tv_usec)/1000000.0;
    cout << bw << endl;
    cout<<"construction time: "<<duration<<" s"<<endl;


    string q="ATG";
    gettimeofday(&start,NULL);
    vector<fmint_t> results=fmindex.search(q);
    gettimeofday(&finish,NULL);
    for(vector<fmint_t>::iterator r_it=results.begin();r_it!=results.end();++r_it)
    {
        cout<<*r_it<<" ";
    }
    cout<<endl;
    duration=finish.tv_sec-start.tv_sec+(finish.tv_usec-start.tv_usec)/1000000.0;
    cout<<"searching time: "<<duration<<" s"<<endl;

*/



    /*----------------------------------------------------------
     * build fm-index and chrindex
     * ---------------------------------------------------------*/

    ChromIndex chrindex;
    FMIndex fmindex;
    {
        gzFile fp_fa;
        kseq_t *seq;
        fp_fa = gzopen(fasta_file.c_str(), "r");
        seq = kseq_init(fp_fa);
        string fasta_seq;
        while ((l = kseq_read(seq)) >= 0)
        {
            FastaHead head(seq->name.s,seq->seq.l);
            chrindex.addhead(head);
            fasta_seq.append(seq->seq.s);
            fasta_seq.push_back(fmindex.returnEOC());
        }
        kseq_destroy(seq);
        gzclose(fp_fa);

        fp_fa = gzopen(fasta_file.c_str(), "r");
        seq = kseq_init(fp_fa);
        string fasta_complement_seq;
        while ((l = kseq_read(seq)) >= 0)
        {
            string chrname(seq->name.s);
            FastaHead head(chrname+"complement",seq->seq.l);
            chrindex.addhead(head);
            fasta_complement_seq.append(seq->seq.s);
            fasta_complement_seq.push_back(fmindex.returnEOC());
        }
        kseq_destroy(seq);
        gzclose(fp_fa);

        chrindex.adjust();
        gettimeofday(&start,NULL);
        std::reverse(fasta_seq.begin(),fasta_seq.end());
        complementSequence(fasta_complement_seq);

        cout<<"fasta length: "<<fasta_seq.length()<<endl;

        fasta_seq.append(fasta_complement_seq);

        fmindex.transform(&fasta_seq);//pass pointer of s
    }
    fmindex.buildIndex();
    //string bw=fmindex.returnBWTString();

    gettimeofday(&finish,NULL);
    duration=finish.tv_sec-start.tv_sec+(finish.tv_usec-start.tv_usec)/1000000.0;
    //cout<<"bwt length: "<<bw.length()<<endl;
    cout<<"construction time: "<<duration<<" s"<<endl;





    bwt_str_size=0,bwt_str_capacity=0,sa_size=0,sa_capacity=0,occ_size=0,C_size=0,C_capacity=0;
    bwt_str_size=sizeof(char)*fmindex.returnBWTString().size();
    bwt_str_capacity=sizeof(char)*fmindex.returnBWTString().capacity();
    sa_size=sizeof(fmint_t)*fmindex.returnSA().size();
    sa_capacity=sizeof(fmint_t)*fmindex.returnSA().capacity();
    occ_size=sizeof(tr1::unordered_map<char,fmint_t>::value_type)*fmindex.returnOcc().size();
    for(fmint_t i=0;i<(fmint_t)fmindex.returnCheckpoints().size();i++)
    {
        C_size+=sizeof(tr1::unordered_map<char,fmint_t>::value_type)*fmindex.returnCheckpoints()[i].size();
    }



    cout<<"bwt_str_size: "<<bwt_str_size<<" bytes"<<endl;
    cout<<"sa_size: "<<sa_size<<" bytes"<<endl;
    cout<<"occ_size: "<<occ_size<<" bytes"<<endl;
    cout<<"C_size: "<<C_size<<" bytes"<<endl;

    /*----------------------------------------------------------
     * test single read
     * ---------------------------------------------------------*/

    /*
        //AGACGGTGCAAACACAAGCCTTGGTTTATCCTGAAATCGTTTTCCTGGGA
        //ATCATACCTAGCTTTGTGTGTGTGGGAATGTGTAGAATACACACACACACACACACACACACACACACACACACA
        string q_str="TCTCCGGGTGCCGCTCGGTGGATCCGTTGAGAATCAGCGCCATGAACGGC";
        vector<u_int8_t> q;
        for(size_t j=0;j<q_str.size();j++)
            q.push_back(encode(q_str[j]));
        gettimeofday(&start,NULL);
        vector<fmint_t> results=fmindex.search(q);
        gettimeofday(&finish,NULL);
        cout<<"positions: ";
        for(auto r_it=results.begin();r_it!=results.end();++r_it)
        {
            cout<<*r_it<<" ";
        }
        cout<<endl;
        cout<<"exact matched number of positions: "<<results.size()<<endl;
        duration=finish.tv_sec-start.tv_sec+(finish.tv_usec-start.tv_usec)/1000000.0;
        cout<<"exact searching time: "<<duration<<" s"<<endl<<endl;


        int mismatch_allowed=3;
        gettimeofday(&start,NULL);
        for(int i=0;i<100;i++)
        {
            results.clear();
            results=fmindex.inexSearchFromSpecifiedPosFindBest(q,mismatch_allowed);
            //results=fmindex.inexSearchFromStart(q,mismatch_allowed);
        }
        gettimeofday(&finish,NULL);
        cout<<"positions: ";
        for(auto r_it=results.begin();r_it!=results.end();++r_it)
        {
            cout<<*r_it<<" ";
        }
        cout<<endl;
        duration=finish.tv_sec-start.tv_sec+(finish.tv_usec-start.tv_usec)/1000000.0;

        if(fmindex.returnMatchStatus())
        {
            cout<<"exact matched number of positions: "<<results.size()<<endl;
            cout<<"exact searching time: "<<duration<<" s"<<endl;
        }
        else
        {
            cout<<"INexact matched number of positions: "<<results.size()<<endl;
            cout<<"INexact searching time: "<<duration<<" s"<<endl;
        }

    */



    /*----------------------------------------------------------
     * build compact trie for reads
     * ---------------------------------------------------------*/



    CTrie ctrie;

    gzFile fp_fq;
    kseq_t *read;

    fp_fq = gzopen(fastq_file.c_str(), "r");
    read = kseq_init(fp_fq);
    count=0;
    cout<<"constructing compact trie......"<<endl;
    gettimeofday(&start,NULL);
    vector<u_int8_t> q;
    q.reserve(130);
    size_t j;
    vector<u_int8_t>::iterator begin,end;
    while ((l = kseq_read(read)) >= 0)
    {
        count++;
        //ctrie.appendRead(read->seq.s,read->name.s);
        for(j=0;j<read->seq.l;j++)
            q.push_back(encode(read->seq.s[j]));
        begin=q.begin();
        end=q.end();
        ctrie.appendRead(begin,end,read->name.s);
        q.clear();
    }
    gettimeofday(&finish,NULL);
    cout<<"read amount: "<<count<<endl;
    kseq_destroy(read);
    gzclose(fp_fq);
    duration=finish.tv_sec-start.tv_sec+(finish.tv_usec-start.tv_usec)/1000000.0;
    cout<<"compact trie's construction time: "<<duration<<" s"<<endl;


    gettimeofday(&start,NULL);
    ctrie.traverse();
    gettimeofday(&finish,NULL);
    duration=finish.tv_sec-start.tv_sec+(finish.tv_usec-start.tv_usec)/1000000.0;
    cout<<"traversing time: "<<duration<<" s"<<endl;
    ctrie.summary();

    ///search by compact trie

    FMIndex *fmidx;
    fmidx=&fmindex;
    ChromIndex *chridx;
    chridx=&chrindex;
    ctrie.specifyIndex(fmidx,chridx);
    int mismatch_allowed=3;
    gettimeofday(&start,NULL);
    ctrie.run(trie_out);
    //ctrie.inexRun(trie_out,mismatch_allowed);
    gettimeofday(&finish,NULL);
    duration=finish.tv_sec-start.tv_sec+(finish.tv_usec-start.tv_usec)/1000000.0;
    cout<<"trie running time: "<<duration<<" s"<<endl;





    /*----------------------------------------------------------
     * search without trie
     * ---------------------------------------------------------*/




    //gzFile fp_fq;
    //kseq_t *read;

    size_t count_matched=0,count_unmatched=0;
    ofstream out_file,out_file1;
    out_file.open(direct_out,ios::out);
    out_file1.open(direct_out.append("_unmatched"),ios::out);
    fp_fq = gzopen(fastq_file.c_str(), "r");
    read = kseq_init(fp_fq);
    count=0;
    //vector<u_int8_t> q;
    //q.reserve(200);
    cout<<"align reads directly......"<<endl;
    gettimeofday(&start,NULL);

    while ((l = kseq_read(read)) >= 0)
    {
        count++;
        //string q(read->seq.s);
        for(j=0;j<read->seq.l;j++)
            q.push_back(encode(read->seq.s[j]));
        vector<fmint_t> results=fmindex.search(q);
        q.clear();
        if(results.size())
        {
            count_matched++;
            out_file<<read->name.s<<" ";
            for(j=0;j<read->seq.l;j++)
            {
               out_file<<q[0];
            }
            //out_file<<read->name.s<<" "<<read->seq.s<<" "<<results.size()<<"  ";
            for(j=0;j<results.size();j++)
            {
                RealPosition rpos=chrindex.calRealPos(results[j],q.size());
                out_file<<rpos.name<<","<<rpos.pos<<" ";
            }
            out_file<<endl;
        }
        else
        {
            count_unmatched++;
            out_file1<<"1"<<" "<<read->seq.s<<endl;
        }
    }
    gettimeofday(&finish,NULL);
    cout<<"read amount: "<<count<<endl;
    kseq_destroy(read);
    gzclose(fp_fq);
    out_file.close();
    duration=finish.tv_sec-start.tv_sec+(finish.tv_usec-start.tv_usec)/1000000.0;
    cout<<"searching time: "<<duration<<" s"<<endl;
    cout<<"count_matched: "<<count_matched<<endl;
    cout<<"count_unmatched: "<<count_unmatched<<endl;



    /*----------------------------------------------------------
     * inexact search without trie
     * ---------------------------------------------------------*/
/*

    gzFile fp_fq;
    kseq_t *read;
    fmindex.resetInexStartStat();
    size_t count_exact_matched=0,count_matched=0,count_unmatched=0;
    ofstream out_file,out_file1;
    out_file.open(direct_out,ios::out);
    out_file1.open(direct_out.append("_unmatched"),ios::out);
    fp_fq = gzopen(fastq_file.c_str(), "r");
    read = kseq_init(fp_fq);
    count=0;
    vector<u_int8_t> q;
    q.reserve(200);
    int mismatch_allowed=3;
    cout<<"align reads directly......"<<endl;
    gettimeofday(&start,NULL);
    size_t j;
    while ((l = kseq_read(read)) >= 0)
    {
        count++;
        //string q(read->seq.s);
        for(j=0;j<read->seq.l;j++)
            q.push_back(encode(read->seq.s[j]));
        vector<fmint_t> results=fmindex.inexSearchFromSpecifiedPosFindBest(q,mismatch_allowed);
        q.clear();
        if(fmindex.returnMatchStatus())//exact match
        {
            count_exact_matched++;
        }
        if(results.size())
        {
            count_matched++;
            out_file<<read->name.s<<" ";
            for(j=0;j<read->seq.l;j++)
            {
               out_file<<q[0];
            }
            out_file<<" "<<results.size()<<endl;
            //out_file<<read->name.s<<" "<<read->seq.s<<" "<<results.size()<<"  ";
            for(j=0;j<results.size();j++)
            {
                RealPosition rpos=chrindex.calRealPos(results[j],q.size());
                out_file<<rpos.name<<","<<rpos.pos<<" ";
            }
            out_file<<endl;
        }
        else
        {
            count_unmatched++;
            out_file1<<"1"<<" "<<read->seq.s<<endl;
        }
    }
    gettimeofday(&finish,NULL);
    cout<<"read amount: "<<count<<endl;
    kseq_destroy(read);
    gzclose(fp_fq);
    out_file.close();
    duration=finish.tv_sec-start.tv_sec+(finish.tv_usec-start.tv_usec)/1000000.0;
    cout<<"searching time: "<<duration<<" s"<<endl;
    cout<<"count_matched: "<<count_matched<<endl;
    cout<<"count_exact_matched: "<<count_exact_matched<<endl;
    cout<<"count_INexact_matched: "<<(int)(count_matched-count_exact_matched)<<endl;
    cout<<"count_unmatched: "<<count_unmatched<<endl;
    cout<<"mapping rate: "<<(float)count_matched/(float)count<<endl;
    //fmindex.printInexStartStat();

*/




    return 0;
}



