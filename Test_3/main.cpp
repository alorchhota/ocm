#include<bits/stdc++.h>
#include"include/ocm.h"
#include <string>
#include <cmath>

using namespace std;
void my_binary(int64_t n);
int64_t cal(string str_k_mer);

string sample_kmer="CCCAGGAATTTCACCCGGGTCG";  // len = 22

//int kmer_len = 22;
int counter=0;

template <typename ccmbase_obj> void update_count_from_file(string file, int len_k_mer, ccmbase_obj &ccmbase_obj_name);
template <typename ccmbase_obj> void update_collision_from_file(string file, int len_k_mer, ccmbase_obj &ccmbase_obj_name, int round);
template <typename ccmbase_obj> void update_count_collision_from_file(string file, int len_k_mer, ccmbase_obj &ccmbase_obj_name, int round, int total_round);


int main(int argc, char *argv[]){
    int kmer_length, counter_h, counter_w, num_round, num_threads, start;
    int np, nh;
    string out_ocm_file, in_fasta_file;
    bool canonicalize;
    string mode(argv[1]);
    string test(argv[2]);
    if(test=="-C")
    {
        canonicalize = false;
        start = 3;

    }
    else{
        canonicalize = true;
        start = 2;
    }
    for(int i=start; i<argc-1; i++)
    {
        string arg(argv[i]);
        string param(argv[++i]);
        if(arg=="-k") kmer_length=stoi(param);
        else if(arg=="-h") counter_h=stoi(param);
        else if(arg=="-w") counter_w=stoi(param);
        else if(arg=="-r") num_round=stoi(param);
        else if(arg=="-t") num_threads=stoi(param);
        else if(arg=="-o") out_ocm_file=param;
        else if(arg=="-fa") in_fasta_file=param;
    }
    nh=counter_h;
    np=log2(counter_w);
//        for(int i=start; i<argc; i++)
//    {
//        cout<<i<<" "<<argv[i]<<endl;
//        cout<<argc<<endl;
//    }
    cout<<kmer_length<<" "<<counter_h<<" "<<counter_w<<" "<<num_round<<" "<<num_threads<<" "<<in_fasta_file<<endl;




     sketch::ocm::ocmbase <uint64_t, sketch::hash::WangHash > sketch1(out_ocm_file, np,nh,137);
    //sketch::ocm::ocmbase <uint64_t, sketch::hash::WangHash > sketch2(12,10,137);

    if(mode=="count")
    {

    //cout<<cal("AAACCGGTCTGTTGGGACCACT")<<endl;  //268931867597
    sketch::ocm::ccmbase <int32_t, sketch::hash::WangHash > cms_obj(np,nh,137,false);   // np = 5, nh = 10, seed = 137,
    sketch::ocm::ccmbase <int32_t, sketch::hash::WangHash > ccms_obj(np,nh,137,true);

    update_count_from_file(in_fasta_file,kmer_length, cms_obj);
    //cout<<"True Count: "<<counter<<endl;
    update_count_from_file(in_fasta_file,kmer_length, ccms_obj);


    //estimation
    //cout<<"CMS Estimate count: "<<cms_obj.est_count(cal(sample_kmer))<<endl; //estimation of the sample k-mer
    //cout<<"CCMS Estimate count: "<<ccms_obj.est_count(cal(sample_kmer))<<endl;


    //construct OCM

    for(int r = 1; r<= num_round; r++){
        clock_t start_time = clock();
        if (r > 1){
             // for all kmers update collision
            update_collision_from_file(in_fasta_file,kmer_length,sketch1,r);
        }
        sketch1.clear_core();
        // for all kmer update count.
        update_count_from_file(in_fasta_file,kmer_length,sketch1);
        clock_t end_time = clock();
        cout<<"Sketch 1: Round "<<r<<" is completed in "<<(double)(end_time - start_time)/CLOCKS_PER_SEC<<" seconds\n";
        //sketch1.showCounters(cal(sample_kmer));
    }
    //end construct OCM
    sketch1.save_core_and_collision();
    cout<<"Offline CMS Estimate count: "<<sketch1.est_count(cal(sample_kmer))<<endl;


    //construct OCCM

//    for(int r = 1; r<= 5; r++){
//        clock_t start_time = clock();
//        if (r > 1){
//            // for all kmers update collision
//            update_collision_from_file(srcfilename,kmer_len,sketch2,r);
//        }
//        sketch2.clear_core();
//        // for all kmer update count.
//        update_count_collision_from_file(srcfilename,kmer_len,sketch2,r,5);
//        clock_t end_time = clock();
//        cout<<"Sketch 2: Round "<<r<<" is completed in "<<(double)(end_time - start_time)/CLOCKS_PER_SEC<<" seconds\n";
//        //sketch2.showCounters(cal(sample_kmer));
//    }
//    //end construct OCCM
//    sketch2.save_core_and_collision();
//    cout<<"Offline CCMS Estimate count: "<<sketch2.est_count(cal(sample_kmer))<<endl;



    //query using original data:
    ifstream qifile("rymv.sim.22mer.counts.txt");
    ofstream ofile("QueryOutput.txt");
    ofstream ocsvfile("Report.csv");
    ocsvfile<<"kmer,true_count,cm,ccm,ocm,occm\n";
    string kmer_str; int file_count;
    while(!qifile.eof()){
        qifile>>kmer_str>>file_count;
        uint64_t kmer_cal = cal(kmer_str);
        //ofile<<kmer_str<<" File Count: "<<file_count<<" CMS-count: "<<cms_obj.est_count(kmer_cal)<<" CCMS-count: "<<ccms_obj.est_count(kmer_cal)<<" OCMS-count: "<<sketch1.est_count(kmer_cal)<<" OCCMS-count: "<<sketch2.est_count(kmer_cal)<<" \n";
        ofile<<kmer_str<<" File Count: "<<file_count<<" CMS-count: "<<cms_obj.est_count(kmer_cal)<<" CCMS-count: "<<ccms_obj.est_count(kmer_cal)<<" OCMS-count: "<<sketch1.est_count(kmer_cal)<<" OCCMS-count: "<<" \n";
        //ocsvfile<<kmer_str<<","<<file_count<<","<<cms_obj.est_count(kmer_cal)<<","<<ccms_obj.est_count(kmer_cal)<<","<<sketch1.est_count(kmer_cal)<<","<<sketch2.est_count(kmer_cal)<<"\n";
        ocsvfile<<kmer_str<<","<<file_count<<","<<cms_obj.est_count(kmer_cal)<<","<<ccms_obj.est_count(kmer_cal)<<","<<sketch1.est_count(kmer_cal)<<","<<"\n";
        //cout<<kmer_str<<" File Count: "<<file_count<<"  CMS-count: "<<cms_obj.est_count(kmer_cal)<<" CCMS-count: "<<ccms_obj.est_count(kmer_cal)<<" OCMS-count: "<<sketch1.est_count(kmer_cal)<<" OCCMS-count: "<<sketch2.est_count(kmer_cal)<<" \n";
        }


        //cout<<sketch1.get_core_size()<<"          "<<sketch1.get_collision_size()<<endl;
        //sketch1.show_data();
        //sketch1.show_core(20);
    }




    else if(mode=="query")
    {
        string input_ocm_file;
        string query_file;
        string query_result_file;

        if(test=="-C")
        {
            canonicalize = false;
            start = 3;

        }
        else{
            canonicalize = true;
            start = 2;
        }
        for(int i=start; i<argc-1; i++)
        {
            string arg(argv[i]);
            string param(argv[++i]);
            if(arg=="-f") input_ocm_file=param;
            else if(arg=="-q") query_file=param;
            else if(arg=="-o") query_result_file=param;
        }
    cout<<"query"<<endl;
    //get_min_collision(5);


        sketch::ocm::ocmbase <uint64_t, sketch::hash::WangHash > sketch3(input_ocm_file, 2,2,137);
        sketch3.load_from_sketch(input_ocm_file);

        cout<<"loaded"<<endl;
        //cout<<"Estimated count using saved ocm file: "<<sketch3.est_count(cal(sample_kmer))<<endl;
        //sketch3.show_core(20);



//        cout<<"in main"<<endl;
//        get_data(input_file);
//        cout<<"after data"<<endl;
//        for(int i=0; i<10; i++)
//            get_nth_core(input_file,i);
//        cout<<"collision"<<endl;
//        for(int i=0; i<10; i++)
//            get_nth_collision(input_file,i);
       // estimate e;
        //e.est_count(6);
       // e.est_count(10);

        ifstream infile(query_file);
        string s;
        ofstream myfile;
        while (infile >> s)
        {

                   cout<<s<<": "<<sketch3.est_count(cal(s))<<endl;


                  myfile.open (query_result_file, std::ios_base::app);
                 myfile <<s<<": "<<sketch3.est_count(cal(s))<<endl;
                  myfile.close();

                   //cout<<"Offline CCMS Estimate count: "<<sketch3.est_count(cal(sample_kmer))<<endl<<endl;
        }






//        cout<<"Offline CMS Estimate count: "<<sketch::ocm::est_count(cal(sample_kmer), nps, nhs, seedseeds, conservatives)<<endl;
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
