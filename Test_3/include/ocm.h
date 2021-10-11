#include <bits/stdc++.h>
#include "hash.h"



namespace sketch {
inline namespace ocm{
using std::allocator;

//OLD work starts

template<typename CounterType=int32_t , typename HasherSetType= HasherSet<WangHash> >
class ccmbase{
    std::vector<CounterType, allocator<CounterType>> core_;     //resisters of the hash table
    uint32_t np_;                                               // no of column (W) is 2^np_
    uint32_t nh_;                                               // number of hash functions
    uint64_t mask_;                                             // and-ing a number with mask will give X mod W
    uint64_t seedseed_;
    const bool conservative_;
    const HasherSetType hf_;
    CounterType       *data()       {return core_.data();}     //data is a pointer to a function
    const CounterType *data() const {return core_.data();}
    size_t size() const {return core_.size();}

public:
    ccmbase(unsigned np, unsigned nh=10, unsigned seedseed=137, bool conservative = false):
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
    }

    void addh_val(uint64_t val) {
        std::vector<CounterType> counts(nh_);
        std::vector<u_int64_t> pos(nh_);
        auto cptr = counts.data();                       //cptr points to the beginning of counts_ vector
        for(unsigned added = 0; added < nh_; added++){
            CounterType hv = hf_(val, added);
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
            CounterType hv = hf_(val, added);
            cptr[added] = hv;               //counts vector now contains hash values
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

}
}
