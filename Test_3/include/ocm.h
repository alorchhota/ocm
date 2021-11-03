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
    std::vector<int> collision_;
    int core_size;                              // will keep track of collision after each round
    uint32_t np_;                                               // no of column (W) is 2^np_
    uint32_t nh_;                                               // number of hash functions
    uint64_t mask_;                                             // and-ing a number with mask will give X mod W
    uint64_t seedseed_;
    bool conservative_;
    const HashStruct hf_;
    std::string output_file;
    std::vector<uint64_t, allocator<uint64_t>> seeds_;
    CounterType       *data()       {return core_.data();}     //data is a pointer to a function
    const CounterType *data() const {return core_.data();}
    size_t size() const {return core_.size();}

public:
    ocmbase(std::string file, unsigned np, unsigned nh=10, unsigned seedseed=137, int core_s=0, bool conservative = false):
        np_(np),
        nh_(nh),
        core_size(core_s),
        mask_((1ull << np_) - 1),
        seedseed_(seedseed),
        hf_(nh_, seedseed),
        conservative_(conservative)
        {
           //std::cout<<"nh  "<<nh_<<std::endl;
            //assert(hf_.size() == nh_);
            output_file=file;
            nh_ += (nh % 2 == 0);
            //std::cout<<"nh  "<<nh_<<std::endl;
            core_.resize(nh_ << np_);
            core_size = core_.size();
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


    void save_core_and_collision()
    {
    //std::cout<<"sameeeee"<<core_size<<" "<<collision_size<<std::endl;;
        //std::cout<<std::endl<<core_.size()<<std::endl<<std::endl;
            std::ofstream outputFile;
            outputFile.open(output_file, std::ios::out | std::ios::binary);

       /// Write Binary File ////
       if(outputFile.is_open())
       {
          // std::cout<<"np:   "<<nh_<<std::endl;
       char true_= '1', false_ = '0';
            outputFile.write(reinterpret_cast<char*>(&np_), sizeof(np_));
           outputFile.write(reinterpret_cast<char*>(&nh_), sizeof(nh_));
         //  outputFile.write(reinterpret_cast<char*>(&mask_), sizeof(mask_));
           outputFile.write(reinterpret_cast<char*>(&seedseed_), sizeof(seedseed_));
           if(conservative_==true) outputFile.write(reinterpret_cast<char*>(&true_), sizeof(true_));
           else outputFile.write(reinterpret_cast<char*>(&false_), sizeof(false_));
            int size_core_ = core_.size();
            outputFile.write(reinterpret_cast<char*>(&size_core_), sizeof(size_core_));
            //int size_collision_ = collision_.size();
            //outputFile.write(reinterpret_cast<char*>(&size_collision_), sizeof(size_collision_));
            for(unsigned i = 0; i < core_.size(); i++){
                    //outputFile.write(reinterpret_cast<char*>(&core_[i]), sizeof(CounterType));
                    //outputFile<<core_[i];
                    //std::cout<<sizeof(int32_t)<<std::endl<<std::endl;
                    outputFile.write(reinterpret_cast<char*>(&core_[i]), sizeof(core_[i]));
                  // if(i<=15) std::cout<<core_[i]<<std::endl;
            }
            for(unsigned i = 0; i < collision_.size(); i++){
                    //std::cout<<sizeof(int)<<std::endl<<std::endl;
                   // outputFile.write(reinterpret_cast<char*>(&collision_[i]), sizeof(CounterType));
                   // std::cout<<collision_[i]<<std::endl;
                    outputFile.write(reinterpret_cast<char*>(&collision_[i]), sizeof(collision_[i]));

            }
            for(unsigned i = 0; i < seeds_.size(); i++){
                    //std::cout<<sizeof(int)<<std::endl<<std::endl;
                   // outputFile.write(reinterpret_cast<char*>(&collision_[i]), sizeof(CounterType));
                  // std::cout<<sizeof(seeds_[i])<<std::endl;
                    outputFile.write(reinterpret_cast<char*>(&seeds_[i]), sizeof(seeds_[i]));

            }

            outputFile.close();
       }



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
//        // DeBUG CODE
//        if(val == 268931867597u){
//            std::cout<<"\nBefore updating: \n";
//            std::cout<<core_[pos[0]]<<" "<<core_[pos[1]]<<" "<<core_[pos[2]]<<" "<<core_[pos[3]]<<" "<<core_[pos[4]]<<" "<<core_[pos[5]]<<" "<<core_[pos[6]]<<" "<<core_[pos[7]]<<" "<<core_[pos[8]]<<" "<<core_[pos[9]]<<" "<<core_[pos[10]]<<" \n";
//            std::cout<<collision_[pos[0]]<<" "<<collision_[pos[1]]<<" "<<collision_[pos[2]]<<" "<<collision_[pos[3]]<<" "<<collision_[pos[4]]<<" "<<collision_[pos[5]]<<" "<<collision_[pos[6]]<<" "<<collision_[pos[7]]<<" "<<collision_[pos[8]]<<" "<<collision_[pos[9]]<<" "<<collision_[pos[10]]<<" \n";
//            std::cout<<"DEBUG round: "<<round<<" min collision is : "<<min_collision<<std::endl;
//        }
//        // END DEBUG CODE

        if(min_collision < current_round - 1){
            //if(val == 268931867597u) std::cout<<" # >= 1 cell without collision in prev round\n";
            CounterType min_count = std::numeric_limits<int>::max();
            for(unsigned added=0; added< nh_; added++){
                if(collision_[pos[added]] == min_collision) min_count = std::min(min_count, core_[pos[added]]);
            }
            //if(val == 268931867597u) std::cout<<"min-count is "<<min_count<<std::endl;
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
            //if(val == 268931867597u) std::cout<<" every cell has a collision in the prev round\n";
            CounterType min_count = std::numeric_limits<int>::max();
            for(unsigned added=0; added< nh_; added++){
                min_count = std::min(min_count, core_[pos[added]]);
            }
//            if(val == 268931867597u){
//                std::cout<<"min-count is "<<min_count<<std::endl;
//            }

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
//        // DeBUG CODE
//        if(val == 268931867597u){
//            std::cout<<"\nAfter updating: \n";
//            std::cout<<core_[pos[0]]<<" "<<core_[pos[1]]<<" "<<core_[pos[2]]<<" "<<core_[pos[3]]<<" "<<core_[pos[4]]<<" "<<core_[pos[5]]<<" "<<core_[pos[6]]<<" "<<core_[pos[7]]<<" "<<core_[pos[8]]<<" "<<core_[pos[9]]<<" "<<core_[pos[10]]<<" \n";
//            std::cout<<collision_[pos[0]]<<" "<<collision_[pos[1]]<<" "<<collision_[pos[2]]<<" "<<collision_[pos[3]]<<" "<<collision_[pos[4]]<<" "<<collision_[pos[5]]<<" "<<collision_[pos[6]]<<" "<<collision_[pos[7]]<<" "<<collision_[pos[8]]<<" "<<collision_[pos[9]]<<" "<<collision_[pos[10]]<<" \n";
//        }
//        // END DEBUG CODE

    }

    void update_collision(uint64_t val, int round){
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
        //std::cout<<"Estimate function is called\n";
        //this->showCounters(val);
        for(unsigned added = 0; added < nh_; added++){
            if( collision_[pos[added]] == min_collision) min_count = core_[pos[added]];
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
    void show_data()
    {
        std::cout<<np_<<std::endl<<nh_<<std::endl<<mask_<<std::endl<<seedseed_<<std::endl<<conservative_<<std::endl<<core_size<<std::endl;;
        //for(int i=0; i<10; i++)
          //  std::cout<<this->get_nth_collision(output_file,i)<<std::endl;

    }



    void load_core(std::string output_file)
    {
    //std::cout<<"sameeeee"<<core_size<<" "<<collision_size<<std::endl;;
        std::ifstream outputFile;
        outputFile.open(output_file, std::ios::in | std::ios::binary);
        if(outputFile.is_open())
       {
       outputFile.seekg(sizeof(uint32_t)*2 + sizeof(uint64_t) + sizeof(char) + sizeof(int), std::ios::beg);
           for(uint32_t i=0; i<core_size; i++)
           {
                outputFile.read(reinterpret_cast<char *>(&core_[i]), sizeof(core_[i]));
                //std::cout<<core_[i]<<std::endl;

           }

        outputFile.close();
       }

    }
    void load_collision(std::string output_file)
    {
        std::ifstream outputFile;
        outputFile.open(output_file, std::ios::in | std::ios::binary);
        if(outputFile.is_open())
       {
       outputFile.seekg(sizeof(uint32_t)*2 + sizeof(uint64_t) + sizeof(char) + sizeof(int) + 8*core_size, std::ios::beg);
           for(uint32_t i=0; i<core_size; i++)
           {
                outputFile.read(reinterpret_cast<char *>(&collision_[i]), sizeof(collision_[i]));
               // std::cout<<collision_[i]<<std::endl;

           }

        outputFile.close();
       }

    }
    void load_seed(std::string output_file)
    {
        std::mt19937_64 mt(seedseed_ + 4);
        while(seeds_.size() < static_cast<unsigned>(nh_)) seeds_.emplace_back(mt());

    }
    void show_core(uint32_t n) const
    {

           for(uint32_t i=0; i<n; i++)
           {
                std::cout<<core_[i]<<std::endl;

           }


    }
    void load_from_sketch(std::string input_file)
    {
        //this->load_data(input_file);
        this->load_core(input_file);
        this->load_collision(input_file);
        this->load_seed(input_file);

    }
    uint32_t get_core_size(){return core_.size();}
    uint32_t get_collision_size(){return collision_.size();}



};

} }
