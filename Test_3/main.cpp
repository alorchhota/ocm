#include<bits/stdc++.h>
#include"include/ocm.h"
using namespace std;
void my_binary(int64_t n);
int64_t cal(string str_k_mer);

string sample_kmer="AAACCGGTCTGTTGGGACCACT";  // len = 22
int counter=0;

template <typename ccmbase_obj> void update_count_from_file(string file, int len_k_mer, ccmbase_obj &ccmbase_obj_name);
template <typename ccmbase_obj> void update_collision_from_file(string file, int len_k_mer, ccmbase_obj &ccmbase_obj_name, int round);
template <typename ccmbase_obj> void update_count_collision_from_file(string file, int len_k_mer, ccmbase_obj &ccmbase_obj_name, int round, int total_round);


int main(){

    sketch::ocm::ccmbase <int32_t, sketch::hash::WangHash > cms_obj(9,10,137,false);   // np = 5, nh = 10, seed = 137,
    sketch::ocm::ccmbase <int32_t, sketch::hash::WangHash > ccms_obj(9,10,137,true);

    update_count_from_file("rymv.sim.fa",22, cms_obj);
    //cout<<"True Count: "<<counter<<endl;
    update_count_from_file("rymv.sim.fa",22, ccms_obj);


    //estimation
    //cout<<"CMS Estimate count: "<<cms_obj.est_count(cal(sample_kmer))<<endl; //estimation of the sample k-mer
    //cout<<"CCMS Estimate count: "<<ccms_obj.est_count(cal(sample_kmer))<<endl;


    //construct OCM
    sketch::ocm::ocmbase <uint64_t, sketch::hash::WangHash > sketch1(9,10,137);
    for(int r = 1; r<= 5; r++){
        if (r > 1){
             // for all kmers update collision
            update_collision_from_file("rymv.sim.fa",22,sketch1,r);
        }
        sketch1.clear_core();
        // for all kmer update count.
        update_count_from_file("rymv.sim.fa",22,sketch1);
        cout<<"Sketch 1: Round "<<r<<" is complete\n";
    }
    //end construct OCM

    //cout<<"Offline CMS Estimate count: "<<sketch1.est_count(cal(sample_kmer))<<endl;


    //construct OCCM
    sketch::ocm::ocmbase <uint64_t, sketch::hash::WangHash > sketch2(9,10,137);
    for(int r = 1; r<= 5; r++){
        if (r > 1){
            // for all kmers update collision
            update_collision_from_file("rymv.sim.fa",22,sketch2,r);
        }
        sketch2.clear_core();
        // for all kmer update count.
        update_count_collision_from_file("rymv.sim.fa",22,sketch2,r,5);
        cout<<"Sketch 2: Round "<<r<<" is complete\n";
    }
    //end construct OCCM
    //cout<<"Offline CCMS Estimate count: "<<sketch2.est_count(cal(sample_kmer))<<endl;

    //query using original data:
    ifstream qifile("rymv.sim.22mer.counts.txt");
    ofstream ofile("QueryOutput.txt");
    string kmer_str; int file_count;
    while(!qifile.eof()){
        qifile>>kmer_str>>file_count;
        uint64_t kmer_cal = cal(kmer_str);
        ofile<<kmer_str<<" File Count: "<<file_count<<" CMS-count: "<<cms_obj.est_count(kmer_cal)<<" CCMS-count: "<<ccms_obj.est_count(kmer_cal)<<" OCMS-count: "<<sketch1.est_count(kmer_cal)<<" OCCMS-count: "<<sketch2.est_count(kmer_cal)<<" \n";
        //cout<<kmer_str<<" File Count: "<<file_count<<"  CMS-count: "<<cms_obj.est_count(kmer_cal)<<" CCMS-count: "<<ccms_obj.est_count(kmer_cal)<<" OCMS-count: "<<sketch1.est_count(kmer_cal)<<" OCCMS-count: "<<sketch2.est_count(kmer_cal)<<" \n";
        }

}

void my_binary(int64_t n) //converting k-mer to binary representation
{
    int a[100];
    string str="";
    int i;
    for(i=0; i<64; i++)
    {
        if(n>0)
        {
            a[i]=n%2;
            n= n/2;
        }
        else a[i]=0;

    }
    cout<<endl<<"Binary: ";
    for(int j=i-1 ;j>=0 ;j--)
    {
        cout<<a[j];
        //str+=to_string(a[i]);
    }
    cout<<endl;
}




int64_t cal(string str_k_mer)
{
    const char* char_array = str_k_mer.c_str();
    int64_t k_mer = 0;
    for(int j=0; char_array[j]!='\0'; j++)
    {
        switch(char_array[j])
        {
            case 'A':
                    k_mer = k_mer<<2; // A=00
                    break;
            case 'T':
                    k_mer = k_mer<<2;
                    k_mer = k_mer | 1;   //T=01
                    break;
            case 'G':
                    k_mer = k_mer<<2;  //G=10
                    k_mer = k_mer | 2;
                    break;
            case 'C':
                    k_mer = k_mer<<2; //C=11
                    k_mer = k_mer | 3;
                    break;
        }
    }
    //cout<<"cal value: "<<k_mer;
    //my_binary(k_mer);
    return k_mer;

}


template <typename ccmbase_obj> void update_count_from_file(string file, int len_k_mer, ccmbase_obj &ccmbase_obj_name)
{
    std::ifstream input(file);
    if(!input.good()){cout<<"Can't open"<<endl;}

    std::string line, name, content;
    while( std::getline( input, line ).good()){ //reading 1 line from input

        if( line.empty() || line[0] == '>' ){
            if( !name.empty() ){
                int l=content.length();
                for(int i=0; i<=l-len_k_mer; i++) //extracting k-mers from the line
                {
                    string part = content.substr(i, len_k_mer);
                    int64_t k_mer = cal(part);
                    /////////Updating part.
                    ccmbase_obj_name.update_count(k_mer);
                    if(k_mer == cal(sample_kmer))counter++;
                }
                name.clear();
            }
            if( !line.empty() ){name = line.substr(1);}
            content.clear();
        }
        else if( !name.empty() ){
            if( line.find(' ') != std::string::npos ){
                name.clear();
                content.clear();
            }
            else {
                content += line;
            }
        }
    }
    if( !name.empty() ){
        //std::cout << name << " : " << content << std::endl;
        int l=content.length();
        for(int i=0; i<=l-len_k_mer; i++)
        {
            string part = content.substr(i, len_k_mer);
            int64_t k_mer = cal(part);

            //// updating part.
            ccmbase_obj_name.update_count(k_mer);
            if(k_mer == cal(sample_kmer))counter++; //manual count of the sample k-mer

        }
    }
}


template <typename ccmbase_obj> void update_collision_from_file(string file, int len_k_mer, ccmbase_obj &ccmbase_obj_name, int round)
{
    std::ifstream input(file);
    if(!input.good()){cout<<"Can't open"<<endl;}

    std::string line, name, content;
    while( std::getline( input, line ).good()){ //reading 1 line from input

        if( line.empty() || line[0] == '>' ){
            if( !name.empty() ){
                int l=content.length();
                for(int i=0; i<=l-len_k_mer; i++) //extracting k-mers from the line
                {
                    string part = content.substr(i, len_k_mer);
                    int64_t k_mer = cal(part);
                    /////////Updating part.
                    ccmbase_obj_name.update_collision(k_mer, round);
                }
                name.clear();
            }
            if( !line.empty() ){name = line.substr(1);}
            content.clear();
        }
        else if( !name.empty() ){
            if( line.find(' ') != std::string::npos ){
                name.clear();
                content.clear();
            }
            else {
                content += line;
            }
        }
    }
    if( !name.empty() ){
        //std::cout << name << " : " << content << std::endl;
        int l=content.length();
        for(int i=0; i<=l-len_k_mer; i++)
        {
            string part = content.substr(i, len_k_mer);
            int64_t k_mer = cal(part);

            //// updating part.
            ccmbase_obj_name.update_collision(k_mer, round);

        }
    }
}

template <typename ccmbase_obj> void update_count_collision_from_file(string file, int len_k_mer, ccmbase_obj &ccmbase_obj_name, int round, int total_round)
{
    std::ifstream input(file);
    if(!input.good()){cout<<"Can't open"<<endl;}

    std::string line, name, content;
    while( std::getline( input, line ).good()){ //reading 1 line from input

        if( line.empty() || line[0] == '>' ){
            if( !name.empty() ){
                int l=content.length();
                for(int i=0; i<=l-len_k_mer; i++) //extracting k-mers from the line
                {
                    string part = content.substr(i, len_k_mer);
                    int64_t k_mer = cal(part);
                    /////////Updating part.
                    ccmbase_obj_name.update_count_collision(k_mer, round, total_round);
                }
                name.clear();
            }
            if( !line.empty() ){name = line.substr(1);}
            content.clear();
        }
        else if( !name.empty() ){
            if( line.find(' ') != std::string::npos ){
                name.clear();
                content.clear();
            }
            else {
                content += line;
            }
        }
    }
    if( !name.empty() ){
        //std::cout << name << " : " << content << std::endl;
        int l=content.length();
        for(int i=0; i<=l-len_k_mer; i++)
        {
            string part = content.substr(i, len_k_mer);
            int64_t k_mer = cal(part);

            //// updating part.
            ccmbase_obj_name.update_count_collision(k_mer, round, total_round);

        }
    }
}
