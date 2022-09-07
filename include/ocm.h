#include <iostream>
#include <vector>
#include <fstream>
#include "hash.h"
#include "compact_vector/compact_vector.hpp"

namespace sketch {
inline namespace ocm{
using std::allocator;

//OLD work starts

template<typename CounterType=int32_t , typename HashStruct = WangHash >
class ccmbase{
    std::vector<CounterType, allocator<CounterType>> core_;     //resisters of the hash table
    uint32_t np_;                                               // no of column (W) is 2^np_
    uint32_t nh_;                                               // number of hash functions
    uint64_t mask_;                                             // and-ing a number with mask will give X mod W
    uint64_t seedseed_;
    const bool conservative_;
    const HashStruct hf_;
    std::vector<uint64_t, allocator<uint64_t>> seeds_;
    CounterType       *data()       {return core_.data();}     //data is a pointer to a function
    const CounterType *data() const {return core_.data();}
    size_t size() const {return core_.size();}

public:
    bool valid_kmer=true;
    ccmbase(unsigned np, unsigned nh=10, unsigned seedseed=137, bool conservative = false):
        np_(np),
        nh_(nh),
        mask_((1ull << np_) - 1),
        seedseed_(seedseed),
        //hf_(nh_, seedseed),
        conservative_(conservative)
        {
            //assert(hf_.size() == nh_);
            nh_ += (nh % 2 == 0);
            core_.resize(nh_ << np_);
            //POST_REQ(core_.size() == (nh_ << np_), "core must be properly sized");            //for throwing exception
            std::mt19937_64 mt(seedseed_ + 4);
            while(seeds_.size() < static_cast<unsigned>(nh_)) seeds_.emplace_back(mt());
            
        }

    void update_count(uint64_t val) {
        std::vector<uint64_t> pos(nh_);                     //cptr points to the beginning of counts_ vector
        for(unsigned added = 0; added < nh_; added++){
            CounterType hv = hf_(val ^ seeds_[added]);
            //cptr[added] = hv;               //counts vector now contains hash values
            pos[added] = (hv & mask_) + (added << np_);   // exact positions where we will increase the counter by one.
        }

        if(conservative_ == false){
            for(unsigned added = 0; added < nh_; added++) core_[pos[added]]++;
        }
        else if(conservative_ == true){

            CounterType min_count = std::numeric_limits<CounterType>::max();
            for(unsigned added = 0; added < nh_; added++){
                min_count = (std::min<CounterType>)(core_[pos[added]], min_count);
            }

            for(unsigned added = 0; added < nh_; added++){
                if(core_[pos[added]] == min_count) core_[pos[added]]++;
            }
        }
        
    }

    CounterType est_count(uint64_t val) const {
        std::vector<uint64_t> pos(nh_);
        for(unsigned added = 0; added < nh_; added++){
            CounterType hv = hf_(val ^ seeds_[added]);
            pos[added] = (hv & mask_) + (added << np_);   // exact positions where we will increase the counter by one.
        }

        CounterType min_count = std::numeric_limits<CounterType>::max();
        for(int i=0; i< nh_; i++){
            min_count = (std::min<CounterType>)(core_[pos[i]], min_count);
        }
        return min_count;
    }
    //support functions
    uint64_t addChar(uint64_t k_mer, char ch){
        switch(ch)
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
        return k_mer;
    }

    int64_t reverse_compliment(uint64_t cal_kmer, int kmer_length)
    {
        uint64_t k_mer = 0;
        uint64_t mask = 3;

        for(int i=0; i<kmer_length; i++)
        {
            switch(cal_kmer & mask)
            {

            case 0:
                k_mer = k_mer<<2;
                k_mer = k_mer | 1;   //A=00->T=01
                break;
            case 1:
                k_mer = k_mer<<2; //T=01->A=00
                break;
            case 2:
                k_mer = k_mer<<2; //G=10->C=11
                k_mer = k_mer | 3;
                break;
            case 3:
                k_mer = k_mer<<2;  //C=11->G=10
                k_mer = k_mer | 2;
                break;

            }
            cal_kmer=cal_kmer>>2;
        }
        return k_mer;
    }
    
    //update count from file
    void update_count_from_file(std::string filename, int len_kmer, bool canonicalize){
        std::ifstream fasta_file(filename);
        int64_t currentKmer = 0; int current_len = 0;
        int chunk_size = 1000;
        char arr_chunk[chunk_size];
        bool isInHeader = false;
        uint64_t MASK = (1ull << (2*len_kmer)) -1;
        // std::cout<<"value of mask: "<<MASK<<std::endl;
        // int chunk_count = 0; int line_no =1;
        while(!fasta_file.eof())
        {
            // chunk_count++;
            fasta_file.read(arr_chunk, chunk_size);
            int i = 0;

            while(i<chunk_size)
            {
                char ch = arr_chunk[i];
                if(ch==EOF)break;            
                // if (ch == '\n')line_no++;
                // std::cout<<"char found "<<ch<<std::endl;
                if(ch == '>'){
                    isInHeader = true;
                    currentKmer = 0; current_len = 0;
                    i++; continue;
                }
                else if (isInHeader==true && ch=='\n'){
                    isInHeader = false;
                    i++; continue;
                }
                if(isInHeader){i+=1;continue;}
                if(ch=='\n' || ch=='\r' || ch==' '){i+=1;continue;}
                if(ch=='N'){
                    currentKmer = 0;current_len = 0;
                    i+=1;continue;
                }
                else
                {
                    // kmer_count++;
                    if(current_len < len_kmer){
                        currentKmer = addChar(currentKmer, ch);
                        current_len++;
                    }
                    else{
                        currentKmer = addChar(currentKmer, ch) & MASK;
                    }
                    // std::cout<<currentKmer<<" "<<current_len<<std::endl;

                    if(current_len==len_kmer)
                    {
                        //GOT KMER -- do the necessary things
                        update_count(currentKmer);
                        if(canonicalize) update_count(reverse_compliment(currentKmer, len_kmer));
                    }
                }
                i+=1;
            }
        }
    }

    void load_from_sketch(std::string input_file_name){
        std::cout<<"loading data from sketch file------\n";
        std::ifstream input_file;
        input_file.open(input_file_name, std::ios::in | std::ios::binary);
        if(input_file.is_open()){
            input_file.seekg(sizeof(uint32_t)*2 + sizeof(uint64_t), std::ios::beg);
            for(uint32_t i=0; i< (nh_<<np_) ; i++){
                input_file.read(reinterpret_cast<char *>(&core_[i]), sizeof(core_[i]));
                //std::cout<<core_[i]<<std::endl;
            }
            std::cout<<"core is read from sketch file------"<<std::endl;
        }

    }

    // output starts
    void save_sketch(std::string output_file_name){
        std::cout<<"saving sketch------\n";
        std::ofstream outputfile;
        outputfile.open(output_file_name, std::ios::out | std::ios::binary);
        // Write Binary File //
        if(outputfile.is_open()){
            outputfile.write(reinterpret_cast<char*>(&np_), sizeof(np_));
            outputfile.write(reinterpret_cast<char*>(&nh_), sizeof(nh_));
            outputfile.write(reinterpret_cast<char*>(&seedseed_), sizeof(seedseed_));


            for(unsigned i = 0; i < core_.size(); i++){
                outputfile.write(reinterpret_cast<char*>(&core_[i]), sizeof(core_[i]));
            }
            outputfile.close();
        }
        std::cout<<"Sketch is saved------\n";
    }
    // output ends
    //support functions ends
};


//old work ends

// Base for Offline Count min Sketch
template<typename CounterType=int32_t , typename HashStruct = WangHash , unsigned int BitSize = 4>
class ocmbase{
    std::vector<CounterType, allocator<CounterType>> core_;     //resisters of the hash table
    //compact::vector<unsigned int, BitSize> collision_;                                // will keep track of collision after each round
    std::vector<unsigned int> collision_;
    uint32_t np_;                                               // no of column (W) is 2^np_
    uint32_t nh_;                                               // number of hash functions
    uint64_t mask_;                                             // and-ing a number with mask will give X mod W
    uint64_t seedseed_;
    const HashStruct hf_;
    std::vector<uint64_t, allocator<uint64_t>> seeds_;
    CounterType       *data()       {return core_.data();}     //data is a pointer to a function
    const CounterType *data() const {return core_.data();}
    size_t size() const {return core_.size();}

public:
    void (ocmbase::*count_function)(uint64_t val, int round, int total_round);
    bool valid_kmer=true;
    ocmbase(unsigned np, unsigned nh=10, unsigned seedseed=137):
        np_(np),
        nh_(nh),
        mask_((1ull << np_) - 1),
        seedseed_(seedseed)
        {
            //assert(hf_.size() == nh_);
            nh_ += (nh % 2 == 0);
            core_.resize(nh_ << np_);
            //POST_REQ(core_.size() == (nh_ << np_), "core must be properly sized");            //for throwing exception
            collision_.resize(nh_ << np_);
            //std::cout<<"Collision DEBUG: \nBitSize= "<<BitSize<<" "<<sizeof(collision_[0])<<"\n";
            for(int i=0;i<collision_.size(); i++) collision_[i] = 0;
            std::mt19937_64 mt(seedseed_ + 4);
            while(seeds_.size() < static_cast<unsigned>(nh_)) seeds_.emplace_back(mt());
        }

    void clear_core(){
        core_.clear();
        core_.shrink_to_fit();
        core_.resize(nh_ << np_);
    }

    void update_count(uint64_t val, int round, int total_round) {                       // update count function
        int min_collision = std::numeric_limits<int>::max();
        //std::vector<CounterType> counts(nh_);
        std::vector<uint64_t> pos(nh_);
        //auto cptr = counts.data();                       //cptr points to the beginning of counts_ vector
        // map kmers to counters + find min collision
        for(unsigned added = 0; added < nh_; added++){
            CounterType hv = hf_(val ^ seeds_[added]);
            //cptr[added] = hv;               //counts vector now contains hash values
            pos[added] = (hv & mask_) + (added << np_);   // exact positions where we will increase the counter by one.
            min_collision = std::min(min_collision, (int)collision_[pos[added]]);
        }

        for(unsigned added = 0; added < nh_; added++){
            if( collision_[pos[added]] == min_collision) core_[pos[added]]++;
        }
    }

    void update_count_collision(uint64_t val, int round, int total_round){
        //std::vector<CounterType> counts(nh_);
        int min_collision = std::numeric_limits<int>::max();
        std::vector<uint64_t> pos(nh_);
        //auto cptr = counts.data();                       //cptr points to the beginning of counts_ vector
        for(unsigned added = 0; added < nh_; added++){
            CounterType hv = hf_(val ^ seeds_[added]);
            //cptr[added] = hv;               //counts vector now contains hash values
            pos[added] = (hv & mask_) + (added << np_);   // exact positions where we will increase the counter by one.
            min_collision = std::min(min_collision, (int)collision_[pos[added]]);
        }

        if(min_collision < round - 1){
            // #>=1 cell without collision in prev round
            CounterType min_count = std::numeric_limits<int>::max();
            for(unsigned added=0; added< nh_; added++){
                if(collision_[pos[added]] == min_collision) min_count = std::min(min_count, core_[pos[added]]);
            }

            for(unsigned added=0; added< nh_; added++){
                if(collision_[pos[added]] == min_collision){
                    // c[i,j] = min(C[i,j]+1, min_count)        Changed code from paper
                    if(core_[pos[added]] == min_count){
                        core_[pos[added]] = min_count+ 1;
                    }
                }
            }
        }

        else{
            // every cell has a collision in the prev round
            CounterType min_count = std::numeric_limits<int>::max();
            for(unsigned added=0; added< nh_; added++){
                min_count = std::min(min_count, core_[pos[added]]);
            }

            for(unsigned added=0; added< nh_; added++){
                if(round < total_round && core_[pos[added]] > min_count){
                    collision_[pos[added]] = round;
                }
                // Changed Code from paper
                //core_[pos[added]] = std::min( core_[pos[added]] + 1, min_count);
                if(core_[pos[added]] == min_count){
                    core_[pos[added]] = min_count + 1;
                }
            }
        }
    }

    void update_collision(uint64_t val, int round, int total_round){
        //std::vector<CounterType> counts(nh_);
        // map kmer to counter and find min collision
        int min_collision = std::numeric_limits<int>::max();
        std::vector<uint64_t> pos(nh_);
        //auto cptr = counts.data();                       //cptr points to the beginning of counts_ vector
        for(unsigned added = 0; added < nh_; added++){
            CounterType hv = hf_(val ^ seeds_[added]);
            //cptr[added] = hv;                             //counts vector now contains hash values
            pos[added] = (hv & mask_) + (added << np_);   // exact positions where we will increase the counter by one.
            min_collision = std::min(min_collision, (int)collision_[pos[added]]);
        }


        if (min_collision >= round-2){
            // find min-count
            CounterType min_count = std::numeric_limits<CounterType>::max();
            for(unsigned added = 0; added < nh_; added++){
                min_count = (std::min<CounterType>)(core_[pos[added]], min_count);
            }

            for(unsigned added = 0; added < nh_; added++){
                if (core_[pos[added]] > min_count){
                    collision_[pos[added]] = round - 1;
                }
            }
        }
    }

    CounterType est_count(uint64_t val) const {
        //std::vector<CounterType> counts(nh_);
        int min_collision = std::numeric_limits<int>::max();
        std::vector<uint64_t> pos(nh_);
        //auto cptr = counts.data();                       //cptr points to the beginning of counts_ vector
        // map kmers and find min collision
        for(unsigned added = 0; added < nh_; added++){
            CounterType hv = hf_(val ^ seeds_[added]);
            //cptr[added] = hv;               //counts vector now contains hash values
            pos[added] = (hv & mask_) + (added << np_);   // exact positions where we will increase the counter by one.
            min_collision = std::min(min_collision, (int)collision_[pos[added]]);
        }

        CounterType min_count = std::numeric_limits<CounterType>::max();

        // min count with smallest collision number
        for(unsigned added = 0; added < nh_; added++){
            if( collision_[pos[added]] == min_collision) min_count = std::min(min_count,core_[pos[added]]);
        }
        return min_count;
    }

    void showSeeds(){
        for(auto seed : seeds_) std::cout<<seed<<" ";
        std::cout<<std::endl;
    }

    void showCounters(uint64_t val){
        std::vector<CounterType> counts(nh_);
        std::vector<u_int64_t> pos(nh_);
        auto cptr = counts.data();                       //cptr points to the beginning of counts_ vector
        std::cout<<"Showing counters for value "<<val<<" : \n";
        for(unsigned added = 0; added < nh_; added++){
            CounterType hv = hf_(val ^ seeds_[added]);
            cptr[added] = hv;               //counts vector now contains hash values
            pos[added] = (hv & mask_) + (added << np_);   // exact positions where we will increase the counter by one.
            std::cout<<"hf: "<<added<<" core: "<<core_[pos[added]]<<" Colision: "<<collision_[pos[added]]<<std::endl;
        }
        std::cout<<std::endl;
    }

    void load_from_sketch(std::string input_file_name){
        std::cout<<"loading data from sketch file------\n";
        std::ifstream input_file;
        input_file.open(input_file_name, std::ios::in | std::ios::binary);
        if(input_file.is_open()){
            input_file.seekg(sizeof(uint32_t)*2 + sizeof(uint64_t), std::ios::beg);            //future work
            for(uint32_t i=0; i< (nh_<<np_) ; i++){
                input_file.read(reinterpret_cast<char *>(&core_[i]), sizeof(core_[i]));
                //std::cout<<core_[i]<<std::endl;
            }
            std::cout<<"core is read from sketch file------"<<std::endl;
            for(uint32_t i=0; i< (nh_<<np_); i++){
                int temp;
                input_file.read(reinterpret_cast<char *>(&temp), sizeof(temp));
                collision_[i]=temp;
                //input_file.read(reinterpret_cast<char *>(&collision_[i]), sizeof(collision_[i]));
            }
        }
        std::cout<<"loading from sketch file is complete------\n";

    }

    // output starts
    void save_sketch(std::string output_file_name){
        std::cout<<"saving sketch------\n";
        std::ofstream outputfile;
        outputfile.open(output_file_name, std::ios::out | std::ios::binary);
        // Write Binary File //
        if(outputfile.is_open()){
            outputfile.write(reinterpret_cast<char*>(&np_), sizeof(np_));
            outputfile.write(reinterpret_cast<char*>(&nh_), sizeof(nh_));
            outputfile.write(reinterpret_cast<char*>(&seedseed_), sizeof(seedseed_));


            for(unsigned i = 0; i < core_.size(); i++){
                outputfile.write(reinterpret_cast<char*>(&core_[i]), sizeof(core_[i]));
            }
            std::cout<<"core is saved------\n";
            for(unsigned i = 0; i < collision_.size(); i++){
                int temp = collision_[i];
                outputfile.write(reinterpret_cast<char*>(&temp), sizeof(temp));
                //outputfile.write(reinterpret_cast<char*>(&collision_[i]), sizeof(collision_[i]));
            }
            outputfile.close();
        }
        std::cout<<"Sketch is saved------\n";
    }
    // output ends


    //support functions
    int64_t addChar(int64_t k_mer, char ch){
        switch(ch)
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
        return k_mer;
    }

    int64_t reverse_compliment(uint64_t cal_kmer, int kmer_length)
    {
        uint64_t k_mer = 0;
        uint64_t mask = 3;

        for(int i=0; i<kmer_length; i++)
        {
            switch(cal_kmer & mask)
            {

            case 0:
                k_mer = k_mer<<2;
                k_mer = k_mer | 1;   //A=00->T=01
                break;
            case 1:
                k_mer = k_mer<<2; //T=01->A=00
                break;
            case 2:
                k_mer = k_mer<<2; //G=10->C=11
                k_mer = k_mer | 3;
                break;
            case 3:
                k_mer = k_mer<<2;  //C=11->G=10
                k_mer = k_mer | 2;
                break;

            }
            cal_kmer=cal_kmer>>2;
        }
        return k_mer;
    }

    // update count from file 
    void update_from_file(std::string filename, int len_kmer, int round, bool canonicalize, int total_round){
        std::ifstream fasta_file(filename);
        int64_t currentKmer = 0; int current_len = 0;
        int chunk_size = 100000;
        char arr_chunk[chunk_size];
        bool isInHeader = false;
        uint64_t MASK = (1ull << (2*len_kmer)) -1;
        // std::cout<<"value of mask: "<<MASK<<std::endl;
        // int chunk_count = 0; int line_no =1;
        while(!fasta_file.eof())
        {
            // chunk_count++;
            fasta_file.read(arr_chunk, chunk_size);
            int i = 0;

            while(i<chunk_size)
            {
                char ch = arr_chunk[i];
                if(ch==EOF)break;            
                // if (ch == '\n')line_no++;
                // std::cout<<"char found "<<ch<<std::endl;
                if(ch == '>'){
                    isInHeader = true;
                    currentKmer = 0; current_len = 0;
                    i++; continue;
                }
                else if (isInHeader==true && ch=='\n'){
                    isInHeader = false;
                    i++; continue;
                }
                if(isInHeader){i+=1;continue;}
                if(ch=='\n' || ch=='\r' || ch==' '){i+=1;continue;}
                if(ch=='N'){
                    currentKmer = 0;current_len = 0;
                    i+=1;continue;
                }
                else
                {
                    // kmer_count++;
                    if(current_len < len_kmer){
                        currentKmer = addChar(currentKmer, ch);
                        current_len++;
                    }
                    else{
                        currentKmer = addChar(currentKmer, ch) & MASK;
                    }
                    // std::cout<<currentKmer<<" "<<current_len<<std::endl;

                    if(current_len==len_kmer)
                    {
                        //GOT KMER -- do the necessary things
                        (this->*(this->count_function))(currentKmer, round, total_round);
                        if(canonicalize)(this->*(this->count_function))(reverse_compliment(currentKmer,len_kmer),round,total_round);
                    }
                }
                i+=1;
            }
        }
    }

};

} }
