#include <bits/stdc++.h>
#include "hash.h"



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
            //clear();
        }

    void update_count(uint64_t val) {
        std::vector<CounterType> counts(nh_);
        std::vector<u_int64_t> pos(nh_);
        auto cptr = counts.data();                       //cptr points to the beginning of counts_ vector
        for(unsigned added = 0; added < nh_; added++){
            CounterType hv = hf_(val ^ seeds_[added]);
            cptr[added] = hv;               //counts vector now contains hash values
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
        std::vector<CounterType> counts(nh_);
        std::vector<u_int64_t> pos(nh_);
        auto cptr = counts.data();                       //cptr points to the beginning of counts_ vector
        for(unsigned added = 0; added < nh_; added++){
            CounterType hv = hf_(val ^ seeds_[added]);
            cptr[added] = hv;               //counts vector now contains hash values
            pos[added] = (hv & mask_) + (added << np_);   // exact positions where we will increase the counter by one.
        }

        CounterType min_count = std::numeric_limits<CounterType>::max();
        //std::cout<<"Counters for val "<<val<<": ";
        for(int i=0; i< nh_; i++){
            //std::cout<<core_[pos[i]]<<" ";
            min_count = (std::min<CounterType>)(core_[pos[i]], min_count);
        }
        //std::cout<<std::endl;
        return min_count;
        }
};


//old work ends

// Base for Offline Count min Sketch
template<typename CounterType=int32_t , typename HashStruct = WangHash >
class ocmbase{
    std::vector<CounterType, allocator<CounterType>> core_;     //resisters of the hash table
    std::vector<int> collision_;                                // will keep track of collision after each round
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
    ocmbase(unsigned np, unsigned nh=10, unsigned seedseed=137, bool conservative = false):
        np_(np),
        nh_(nh),
        mask_((1ull << np_) - 1),
        seedseed_(seedseed),
        hf_(nh_, seedseed),
        conservative_(conservative)
        {
            //assert(hf_.size() == nh_);
            nh_ += (nh % 2 == 0);
            core_.resize(nh_ << np_);
            //POST_REQ(core_.size() == (nh_ << np_), "core must be properly sized");            //for throwing exception
            collision_.resize(nh_ << np_);
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
        std::vector<CounterType> counts(nh_);
        std::vector<u_int64_t> pos(nh_);
        auto cptr = counts.data();                       //cptr points to the beginning of counts_ vector
        for(unsigned added = 0; added < nh_; added++){
            CounterType hv = hf_(val ^ seeds_[added]);
            cptr[added] = hv;               //counts vector now contains hash values
            pos[added] = (hv & mask_) + (added << np_);   // exact positions where we will increase the counter by one.
        }

        int min_collision = std::numeric_limits<int>::max();
        for(unsigned added = 0; added < nh_; added++){
            min_collision = std::min(min_collision, collision_[pos[added]]);
        }
        //std::cout<<"DEBUG "<<"min collision is : "<<min_collision<<std::endl;
        for(unsigned added = 0; added < nh_; added++){
            if( collision_[pos[added]] == min_collision) core_[pos[added]]++;
            //else std::cout<<"Collision Deceted while updating"<<std::endl;
        }
    }

    void update_count_collision(uint64_t val, int current_round, int total_round){
        std::vector<CounterType> counts(nh_);
        std::vector<u_int64_t> pos(nh_);
        auto cptr = counts.data();                       //cptr points to the beginning of counts_ vector
        for(unsigned added = 0; added < nh_; added++){
            CounterType hv = hf_(val ^ seeds_[added]);
            cptr[added] = hv;               //counts vector now contains hash values
            pos[added] = (hv & mask_) + (added << np_);   // exact positions where we will increase the counter by one.
        }

        // get smallest collision round among counters
        int min_collision = std::numeric_limits<int>::max();
        for(unsigned added = 0; added < nh_; added++){
            min_collision = std::min(min_collision, collision_[pos[added]]);
        }
        //std::cout<<"DEBUG "<<"min collision is : "<<min_collision<<std::endl;
        if(min_collision < current_round - 1){
            // # >= 1 cell without collision in prev round
            CounterType min_count = std::numeric_limits<int>::max();
            for(unsigned added=0; added< nh_; added++){
                if(collision_[pos[added]] == min_collision) min_count = std::min(min_count, core_[pos[added]]);
            }

            for(unsigned added=0; added< nh_; added++){
                if(collision_[pos[added]] == min_collision){
                    // c[i,j] = min(C[i,j]+1, min_count)        Changed code from paper
                    core_[pos[added]] = min_count+ 1;
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
                core_[pos[added]] = min_count + 1;
            }

        }
    }

    void update_collision(int64_t val, int round){
        std::vector<CounterType> counts(nh_);
        std::vector<u_int64_t> pos(nh_);
        auto cptr = counts.data();                       //cptr points to the beginning of counts_ vector
        for(unsigned added = 0; added < nh_; added++){
            CounterType hv = hf_(val ^ seeds_[added]);
            cptr[added] = hv;                             //counts vector now contains hash values
            pos[added] = (hv & mask_) + (added << np_);   // exact positions where we will increase the counter by one.
        }

        int min_collision = std::numeric_limits<int>::max();
        for(unsigned added = 0; added < nh_; added++){
            min_collision = std::min(min_collision, collision_[pos[added]]);
        }
        //std::cout<<"DEBUG "<<"min collision is : "<<min_collision<<std::endl;

        if (min_collision >= round-2){
            //std::cout<<"Came to update collision"<<std::endl;
            // find min-count
            CounterType min_count = std::numeric_limits<CounterType>::max();
            for(unsigned added = 0; added < nh_; added++){
                min_count = (std::min<CounterType>)(core_[pos[added]], min_count);
            }

            for(unsigned added = 0; added < nh_; added++){
                if (core_[pos[added]] > min_count){
                    collision_[pos[added]] = round - 1;
                    //std::cout<<"Collision Detected"<<std::endl;
                }

            }
        }
    }

    CounterType est_count(uint64_t val) const {
        std::vector<CounterType> counts(nh_);
        std::vector<u_int64_t> pos(nh_);
        auto cptr = counts.data();                       //cptr points to the beginning of counts_ vector
        for(unsigned added = 0; added < nh_; added++){
            CounterType hv = hf_(val ^ seeds_[added]);
            cptr[added] = hv;               //counts vector now contains hash values
            pos[added] = (hv & mask_) + (added << np_);   // exact positions where we will increase the counter by one.
        }

        CounterType min_count = std::numeric_limits<CounterType>::max();

        int min_collision = std::numeric_limits<int>::max();
        for(unsigned added = 0; added < nh_; added++){
            min_collision = std::min(min_collision, collision_[pos[added]]);
        }
        //std::cout<<"DEBUG "<<"min collision is : "<<min_collision<<std::endl;
        for(unsigned added = 0; added < nh_; added++){
            if( collision_[pos[added]] == min_collision) min_count = core_[pos[added]];
        }
        return min_count;
    }

    void showSeeds(){
        for(auto seed : seeds_) std::cout<<seed<<" ";
        std::cout<<std::endl;
    }

};

} }
