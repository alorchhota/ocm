#include<iostream>
#include<ctime>
#include"include/ocm.h"
using namespace std;
uint64_t cal(string str_k_mer);
clock_t start_time,end_time,total_time;

string sample_kmer="AAACCGGTCTGTTGGGACCACT";  // len = 22
string INPUT_FASTA_FILE = "",OUTPUT_SKETCH_FILE="";   //  "rymv.sim.fa"; //
unsigned int kmer_len = 0;
int counter=0;
unsigned int NP = 20, NH = 7, TOTAL_ROUND = 4,SEED=137, counter_w,num_threads;
bool CANONICALIZE = true, CONSERVATIVE = false;

int main(int argc, char *argv[]){
    string mode(argv[1]);
    if( mode == "count"){
        for(int i= 2; i<argc-1; i++){
            string arg(argv[i]);
            //string param(argv[++i]);
            if(arg=="-k") kmer_len=stoi(argv[++i]);
            else if(arg=="-h") NH =stoi(argv[++i]);
            else if(arg=="-w") NP = log2(stoi(argv[++i]));
            else if(arg=="-n") TOTAL_ROUND=stoi(argv[++i]);
            else if(arg=="-t") num_threads=stoi(argv[++i]);
            else if(arg=="-o") OUTPUT_SKETCH_FILE = argv[++i];
            else if(arg=="-fa") INPUT_FASTA_FILE=argv[++i];
            else if(arg=="-r") CANONICALIZE = false;
            else if(arg=="-c") CONSERVATIVE = true;
        }
        if(kmer_len <= 0 || INPUT_FASTA_FILE=="" || OUTPUT_SKETCH_FILE==""){
            cout<<"INVALID INPUT"; return 0;
        }
        //end inouting
        //cout<<NP<<" "<<NH<<" "<<TOTAL_ROUND<<" "<<OUTPUT_SKETCH_FILE<<" "<<INPUT_FASTA_FILE<<" "<<CANONICALIZE<<" ";
        if(!CONSERVATIVE){
            //construct OCM
            total_time = clock();
            sketch::ocm::ocmbase <uint64_t, sketch::hash::WangHash,2> sketch1(NP,NH,137);
            for(int r = 0; r< TOTAL_ROUND; r++){
                sketch1.count_function = &sketch::ocm::ocmbase<uint64_t, sketch::hash::WangHash,2>::update_collision;
                if (r > 0){
                    // for all kmers update collision
                    start_time = clock();
                    sketch1.update_from_file(INPUT_FASTA_FILE, kmer_len, r, CANONICALIZE, 0);
                    end_time = clock();
                    cout<<"Updating collision for ocms round "<<r<<" is completed in "<<(double)(end_time - start_time)/CLOCKS_PER_SEC<<" seconds\n";
                }

                sketch1.clear_core();
                // for all kmer update count.
                sketch1.count_function = &sketch::ocm::ocmbase<uint64_t, sketch::hash::WangHash,2>::update_count;
                start_time = clock();
                sketch1.update_from_file(INPUT_FASTA_FILE, kmer_len, r, CANONICALIZE, 0);
                end_time = clock();
                cout<<"Updating count for ocms round "<<r<<" is completed in "<<(double)(end_time - start_time)/CLOCKS_PER_SEC<<" seconds\n";
                //sketch1.showCounters(cal(sample_kmer));
            }
            cout<<"Constructing ocms is completed in :"<<(double)(clock()-total_time)/CLOCKS_PER_SEC<<" seconds\n";
            //end construct OCM
            sketch1.save_sketch(OUTPUT_SKETCH_FILE);
        }
        else{
        //construct occm
        total_time = clock();
        sketch::ocm::ocmbase <uint64_t, sketch::hash::WangHash,2> sketch2(NP,NH,137);
        sketch2.count_function = &sketch::ocm::ocmbase<uint64_t, sketch::hash::WangHash,2>::update_collision;
        for(int r = 0; r< TOTAL_ROUND; r++){
            if (r > 0){
                // for all kmers update collision
                clock_t start_time = clock();
                sketch2.update_from_file(INPUT_FASTA_FILE, kmer_len, r, CANONICALIZE, 0);
                clock_t end_time = clock();
                cout<<"Updating collision for occms round "<<r<<" is completed in "<<(double)(end_time - start_time)/CLOCKS_PER_SEC<<" seconds\n";
            }
            sketch2.clear_core();
            sketch2.count_function = &sketch::ocm::ocmbase<uint64_t, sketch::hash::WangHash,2>::update_count_collision;
            // for all kmer update count collision
            start_time = clock();
            //for(auto kmer: kmers_vec) sketch2.update_count_collision(kmer,r,TOTAL_ROUND);
            sketch2.update_from_file(INPUT_FASTA_FILE, kmer_len, r, CANONICALIZE, TOTAL_ROUND);
            end_time = clock();
            cout<<"Updating count-collision for occms round "<<r<<" is completed in "<<(double)(end_time - start_time)/CLOCKS_PER_SEC<<" seconds\n";
            //sketch2.showCounters(cal(sample_kmer));

        }
        //end occm
        cout<<"Constructing occm is completed in :"<<(double)(clock()-total_time)/CLOCKS_PER_SEC<<" seconds\n";
        sketch2.save_sketch(OUTPUT_SKETCH_FILE);
        }
    }
    if (mode == "query"){
        string input_sketch_name;
        string query_file_name;
        string query_result_file_name;

        for(int i=2; i<argc-1; i++)
        {
            string arg(argv[i]);
            string param(argv[++i]);
            if(arg=="-f") input_sketch_name=param;
            else if(arg=="-q") query_file_name=param;
            else if(arg=="-o") query_result_file_name=param;
        }
        cout<<"query"<<endl;
        ifstream input_sketch_file(input_sketch_name, std::ios::in | std::ios::binary);
        if(!input_sketch_file.is_open()){
            cout<<"Can't open sketch file\n";
            return 0;
        }
        input_sketch_file.read(reinterpret_cast<char *>(&NP), sizeof(NP));
        input_sketch_file.read(reinterpret_cast<char *>(&NH), sizeof(NH));
        input_sketch_file.read(reinterpret_cast<char *>(&SEED), sizeof(SEED));
        input_sketch_file.close();

        //cout<<"Read input sketch file "<<NP<<" "<<NH<<" "<<SEED<<"\n";

        sketch::ocm::ocmbase <uint64_t, sketch::hash::WangHash,2> query_sketch(NP,NH,SEED);
        query_sketch.load_from_sketch(input_sketch_name);
        ifstream infile(query_file_name);
        if(!infile.good()){cout<<"Couldn't open input query file\n";}
        string kmer;
        int true_count;
        ofstream query_result_file;
        query_result_file.open(query_result_file_name, ios::out);
        if(!query_result_file.good()){cout<<"Couldn't open output query file file\n";}
        query_result_file<<"kmer,true_count,estimated_count\n";
        while (infile >> kmer >> true_count){
            //cout<<kmer<<","<<true_count<<","<<query_sketch.est_count(cal(kmer))<<endl;
            query_result_file <<kmer<<","<<true_count<<","<<query_sketch.est_count(cal(kmer))<<endl;
        }
    }
    return 0;
}


uint64_t cal(string str_k_mer)
{
    const char* char_array = str_k_mer.c_str();
    uint64_t k_mer = 0;
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
