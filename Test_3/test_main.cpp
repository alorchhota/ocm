#include<bits/stdc++.h>
#include"include/hash.h"
bool conservative_update = false;
using namespace std;
int64_t kmar_to_int(string str_k_mer);



// Thomas Wang hash
// Original site down, available at https://naml.us/blog/tag/thomas-wang
// This is our core 64-bit hash.
// It a bijection within [0,1<<64)
// and can be inverted with irving_inv_hash.
struct WangHash {
    template<typename...Args> WangHash(Args &&...) {}
    INLINE auto operator()(uint64_t key) const {
          key = (~key) + (key << 21); // key = (key << 21) - key - 1;
          key = key ^ (key >> 24);
          key = (key + (key << 3)) + (key << 8); // key * 265
          key = key ^ (key >> 14);
          key = (key + (key << 2)) + (key << 4); // key * 21
          key = key ^ (key >> 28);
          key = key + (key << 31);
          return key;
    }
    static constexpr auto hash(uint64_t key) {
          key = (~key) + (key << 21); // key = (key << 21) - key - 1;
          key = key ^ (key >> 24);
          key = (key + (key << 3)) + (key << 8); // key * 265
          key = key ^ (key >> 14);
          key = (key + (key << 2)) + (key << 4); // key * 21
          key = key ^ (key >> 28);
          key = key + (key << 31);
          return key;
    }
    INLINE auto operator()(int64_t key) const {return operator()(uint64_t(key));}
    INLINE uint32_t operator()(uint32_t key) const {
        key += ~(key << 15);
        key ^=  (key >> 10);
        key +=  (key << 3);
        key ^=  (key >> 6);
        key += ~(key << 11);
        key ^=  (key >> 16);
        return key;
    }
    INLINE auto operator()(int32_t key) const {return operator()(uint32_t(key));}
#ifdef _VEC_H__
    INLINE Type operator()(Type element) const {
        VType key = Space::add(Space::slli(element, 21), ~element); // key = (~key) + (key << 21);
        key = Space::srli(key.simd_, 24) ^ key.simd_; //key ^ (key >> 24)
        key = Space::add(Space::add(Space::slli(key.simd_, 3), Space::slli(key.simd_, 8)), key.simd_); // (key + (key << 3)) + (key << 8);
        key = key.simd_ ^ Space::srli(key.simd_, 14);  // key ^ (key >> 14);
        key = Space::add(Space::add(Space::slli(key.simd_, 2), Space::slli(key.simd_, 4)), key.simd_); // (key + (key << 2)) + (key << 4); // key * 21
        key = key.simd_ ^ Space::srli(key.simd_, 28); // key ^ (key >> 28);
        key = Space::add(Space::slli(key.simd_, 31), key.simd_);    // key + (key << 31);
        return key.simd_;
    }
#endif
#if VECTOR_WIDTH > 16
    INLINE auto operator()(__m128i key) const {
        key = _mm_add_epi64(~key, _mm_slli_epi64(key, 21)); // key = (key << 21) - key - 1;
        key ^= _mm_srli_epi64(key, 24);
        key = _mm_add_epi64(key, _mm_add_epi64(_mm_slli_epi64(key, 3), _mm_slli_epi64(key, 8)));
        key ^= _mm_srli_epi64(key, 14);
        key = _mm_add_epi64(_mm_add_epi64(key, _mm_slli_epi64(key, 2)), _mm_slli_epi64(key, 4));
        key ^= _mm_srli_epi64(key, 28);
        key = _mm_add_epi64(_mm_slli_epi64(key, 31), key);
        return key;
    }
#endif
    INLINE uint64_t inverse(uint64_t key) const {
        // https://naml.us/blog/tag/thomas-wang
        uint64_t tmp;
        // Invert key = key + (key << 31)
        tmp = key-(key<<31);
        key = key-(tmp<<31);
        // Invert key = key ^ (key >> 28)
        tmp = key^key>>28;
        key = key^tmp>>28;
        // Invert key *= 21
        key *= UINT64_C(14933078535860113213);
        // Invert key = key ^ (key >> 14)
        tmp = key^key>>14;
        tmp = key^tmp>>14;
        tmp = key^tmp>>14;
        key = key^tmp>>14;
        // Invert key *= 265
        key *= UINT64_C(15244667743933553977);
        // Invert key = key ^ (key >> 24)
        tmp = key^key>>24;
        key = key^tmp>>24;
        // Invert key = (~key) + (key << 21)
        tmp = ~key;
        tmp = ~(key-(tmp<<21));
        tmp = ~(key-(tmp<<21));
        key = ~(key-(tmp<<21));
        return key;
    }
    template<typename T>
    INLINE uint64_t operator()(const T&x) const {return this->operator()(static_cast<uint64_t>(x));}
};

template<typename BaseHash>
struct SeededHash: public BaseHash {
    const uint64_t seed_;
    template<typename...A>
    SeededHash(uint64_t seed, A &&...a): seed_(seed), BaseHash(std::forward<A>(a)...) {}
    auto operator()(uint64_t item) const {return BaseHash::operator()(item ^ seed_);}
};


template<typename Hasher=SeededHash<WangHash>>
struct HasherSet {
    std::vector<Hasher> hashers_;
    HasherSet(size_t nh, uint64_t seedseed=137) {
        std::mt19937_64 mt(seedseed);
        while(hashers_.size() < nh)
            hashers_.emplace_back(mt());
    }
    size_t size() const {return hashers_.size();}
    uint64_t operator()(uint64_t v, unsigned ind) const {
        return hashers_[ind](v);
    }
    //uint64_t operator()(uint64_t v) const {throw std::runtime_error("Should not be called.");}
};



template<typename CounterType=int32_t , typename HasherSetType= HasherSet<WangHash> >
class ccmbase{
    std::vector<CounterType, allocator<CounterType>> core_;     //resisters of the hash table
    uint32_t np_;                                               // no of column (W) is 2^np_
    uint32_t nh_;                                               // number of hash functions
    uint64_t mask_;                                             // and-ing a number with mask will give X mod W
    uint64_t seedseed_;
    const HasherSetType hf_;
    CounterType       *data()       {return core_.data();}     //data is a pointer to a function
    const CounterType *data() const {return core_.data();}
    size_t size() const {return core_.size();}

public:
    ccmbase(unsigned np, unsigned nh=10, unsigned seedseed=137):
        np_(np),
        nh_(nh),
        mask_((1ull << np_) - 1),
        seedseed_(seedseed),
        hf_(nh_, seedseed)
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

        if(conservative_update == false){
            for(unsigned added = 0; added < nh_; added++) core_[pos[added]]++;
        }
    }

    CounterType est_count(uint64_t val) const {
        std::vector<CounterType> counts(nh_);
        std::vector<u_int64_t> pos(nh_);
        auto cptr = counts.data();                         //cptr points to the beginning of counts_ vector
        for(unsigned added = 0; added < nh_; added++){
            CounterType hv = hf_(val, added);
            cptr[added] = hv;                             //counts vector now contains hash values
            pos[added] = (hv & mask_) + (added << np_);   // exact positions where we will increase the counter by one.
        }

        CounterType min_count = std::numeric_limits<CounterType>::max();

        for(int i=0; i< nh_; i++){
            min_count =min(core_[pos[i]], min_count);
        }
        return min_count;
        }
};


int main(){

    ccmbase <int32_t, HasherSet<WangHash> >c1(5,10);

    //adding for test
    for(int i=0; i<90; i++){
        c1.addh_val(108);
    }

    cout<<c1.est_count(108)<<endl;


//    //read_file("rymv.sim.fa",len_k_mer);
}


int64_t kmar_to_int(string str_k_mer){
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
    return k_mer;

}


