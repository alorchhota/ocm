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
        //std::vector<CounterType> counts(nh_);
        std::vector<uint64_t> pos(nh_);
        //auto cptr = counts.data();                       //cptr points to the beginning of counts_ vector
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

    CounterType est_count(uint64_t val, int round) const {
        //std::vector<CounterType> counts(nh_);
        std::vector<uint64_t> pos(nh_);
        //auto cptr = counts.data();                       //cptr points to the beginning of counts_ vector
        for(unsigned added = 0; added < nh_; added++){
            CounterType hv = hf_(val ^ seeds_[added]);
            //cptr[added] = hv;               //counts vector now contains hash values
            pos[added] = (hv & mask_) + (added << np_);   // exact positions where we will increase the counter by one.
        }

        CounterType min_count = std::numeric_limits<CounterType>::max();
        for(int i=0; i< nh_; i++){
            min_count = (std::min<CounterType>)(core_[pos[i]], min_count);
        }
        return min_count;
        }
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


    int64_t process_first_kmer(int& start, int len_kmer, char arr_kmer[])
    {
        int64_t k_mer = 0;
        size_t len = strlen(arr_kmer);
        if(start+len_kmer>len)
        {
            valid_kmer=false;
            return k_mer;

        }
        int t=0;

        for(int l=0; l<len_kmer; l++)     //iterating over 1 kmer
        {
            switch(arr_kmer[start+l])
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
                default:
                    k_mer=0;
                    start+=(l+1);
                    l=-1;

            }

        }
        return k_mer;



    }

void process_one_line(char arr_kmer[], int arr_len, int len_kmer,int round, bool canonicalize, int total_round)
{
    int64_t mask_1, mask_2;
    if(2*len_kmer>30)
    {

        mask_1 = (1 << 30);
        mask_1 = ~(mask_1 << (2*(len_kmer-1)-30));
        mask_2 = (1 << 30);
        mask_2 = ~(mask_2 << ((2*(len_kmer-1)-30)+1));

    }
    else{
        mask_1 = ~(1<<(2*(len_kmer-1)));
        mask_2 = ~(1<<((2*(len_kmer-1))+1));
    }
    int start=0;
    int64_t k_mer = process_first_kmer(start, len_kmer, arr_kmer);
    if(!valid_kmer) return;
    (this->*(this->count_function))(k_mer, round, total_round);
    if(canonicalize)(this->*(this->count_function))(reverse_compliment(k_mer,len_kmer),round,total_round);



    for(int l=start+1; arr_kmer[l+len_kmer-1]!='\0'; l++)     //iterating over 1 kmer
    {
        k_mer = k_mer & mask_1;
        k_mer = k_mer & mask_2;
            switch(arr_kmer[l+len_kmer-1])
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
                default:
                    start=l+len_kmer;
                    if(arr_kmer[start]=='\0') return;
                    k_mer = process_first_kmer(start, len_kmer, arr_kmer);
                    l=start+1;
                    if(!valid_kmer) return;
                    (this->*(this->count_function))(k_mer, round, total_round);
                    if(canonicalize)(this->*(this->count_function))(reverse_compliment(k_mer,len_kmer),round,total_round);



            }

        (this->*(this->count_function))(k_mer, round, total_round);
        if(canonicalize)(this->*(this->count_function))(reverse_compliment(k_mer,len_kmer),round,total_round);

        }
        k_mer = 0;
}

    // update count from file muti-threaded
    void update_from_file(std::string filename, int len_kmer, int round, bool canonicalize, int total_round){
        int char_size = 500;

        int64_t k_mer = 0;

        std::ifstream file(filename);
        char arr_kmer[char_size];
        int global_file_pos=0;
        int global_currect_pos=0;
        int flag=0;
        file.seekg(0, std::ios_base::end);
        int length = file.tellg();
        file.seekg(0);
        bool is_last_chunk=false;

        while(!file.eof())          //reading n lines at a time
        {
            if(global_currect_pos>=length) return;
            flag=0;
            arr_kmer[0]='\0';
            int read_length;
            if(length-global_currect_pos>char_size)
            {
                read_length=char_size-1;

            }
            else{
                is_last_chunk=true;
                read_length=length-global_currect_pos;
            }
            file.read(arr_kmer, read_length);

            arr_kmer[read_length]='\0';
            int cur_idx=0;
            if(arr_kmer[cur_idx]=='\n') cur_idx++;
            int one_line_len=0;
            char one_line[char_size];
            while(true)       //extracting 1 line at a time
            {
                if(global_currect_pos>=length) return;
                int id_length=0;
                if(arr_kmer[cur_idx]=='>')
                {
                    while(arr_kmer[cur_idx++]!='\n')
                    {
                        id_length++;

                        if(cur_idx==char_size-1){
                            file.seekg(global_currect_pos);
                            flag=1;
                            cur_idx=0;
                            break;

                        }
                    }
                    if(flag==0) global_currect_pos+=(id_length+1);
                }

                int tmp=0;
                if(flag==1)
                {
                    one_line[0]='\0';
                    one_line_len=0;
                    break;
                }
                int num_of_lines=-1;

                while(arr_kmer[cur_idx]!='>')
                {


                    if(cur_idx==char_size-1){
                        if(is_last_chunk)
                        {
                            global_currect_pos+=(one_line_len);
                            global_currect_pos+=num_of_lines;
                            file.seekg(global_currect_pos);
                            break;


                        }

                        file.seekg(global_currect_pos);
                        flag=1;
                        cur_idx=0;
                        break;
                    }
                    if(arr_kmer[cur_idx]!='\n')
                        one_line[one_line_len++]=arr_kmer[cur_idx++];
                    else{
                        num_of_lines++;
                        cur_idx++;
                    }
                }
                if(flag==0)
                {
                    one_line[one_line_len]='\0';
                    global_currect_pos+=one_line_len+1;
                    global_currect_pos+=(num_of_lines);


                    //print_array(one_line);
                    valid_kmer=true;
                    process_one_line(one_line, one_line_len, len_kmer, round, canonicalize, total_round);


                    one_line[0]='\0';
                    one_line_len=0;

                }
                else{
                    break;

                }

            }

        }







    }
    //support functions


    uint64_t reverse_compliment(uint64_t cal_kmer, int kmer_length)
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

    //support functions ends

};

} }
