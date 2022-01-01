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

    CounterType est_count(uint64_t val) const {
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

    void update_count(uint64_t val) {                       // update count function
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

    void update_count_collision(uint64_t val, int current_round, int total_round){
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

        if(min_collision < current_round - 1){
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
                if(current_round < total_round && core_[pos[added]] > min_count){
                    collision_[pos[added]] = current_round;
                }
                // Changed Code from paper
                //core_[pos[added]] = std::min( core_[pos[added]] + 1, min_count);
                if(core_[pos[added]] == min_count){
                    core_[pos[added]] = min_count + 1;
                }
            }
        }
    }

    void update_collision(uint64_t val, int round){
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
        //this->load_data(input_file);
        //this->load_core(input_file);
        std::ifstream input_file;
        input_file.open(input_file_name, std::ios::in | std::ios::binary);
        if(input_file.is_open()){
            input_file.seekg(sizeof(uint32_t)*2 + sizeof(uint64_t), std::ios::beg);            //future work
            for(uint32_t i=0; i< (nh_<<np_) ; i++){
                input_file.read(reinterpret_cast<char *>(&core_[i]), sizeof(core_[i]));
                //std::cout<<core_[i]<<std::endl;
            }
            for(uint32_t i=0; i< (nh_<<np_); i++){
                input_file.read(reinterpret_cast<char *>(&collision_[i]), sizeof(collision_[i]));
               // std::cout<<collision_[i]<<std::endl;
            }
        }

    }

    // output starts
    void save_sketch(std::string output_file_name){
        std::ofstream outputfile;
        outputfile.open(output_file_name, std::ios::out | std::ios::binary);

        // Write Binary File //
        if(outputfile.is_open()){
            // std::cout<<"np:   "<<nh_<<std::endl;
            //char true_= '1', false_ = '0';
            outputfile.write(reinterpret_cast<char*>(&np_), sizeof(np_));
            outputfile.write(reinterpret_cast<char*>(&nh_), sizeof(nh_));
            outputfile.write(reinterpret_cast<char*>(&seedseed_), sizeof(seedseed_));
//            if(conservative_==true) outputFile.write(reinterpret_cast<char*>(&true_), sizeof(true_));
//            else outputFile.write(reinterpret_cast<char*>(&false_), sizeof(false_));
            for(unsigned i = 0; i < core_.size(); i++){
                outputfile.write(reinterpret_cast<char*>(&core_[i]), sizeof(core_[i]));
            }
            for(unsigned i = 0; i < collision_.size(); i++){
                outputfile.write(reinterpret_cast<char*>(&collision_[i]), sizeof(collision_[i]));
            }
            outputfile.close();
        }
    }
    // output ends

};

} }
