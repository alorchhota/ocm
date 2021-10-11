#include<bits/stdc++.h>
#include"include/ocm.h"
using namespace std;


int main(){
    cout<<"--TEST 3 --"<<endl;
    vector<uint64_t> kmers;

    unsigned int k;
    ifstream myfile("randints.txt");
    for(int i=0;i<1000; i++){
        myfile >> k;
        kmers.push_back(k);
    }

    //for(int i=0;i<kmers.size(); i++) cout<<kmers[i]<<endl;

    sketch::ocm::ccmbase <int32_t, sketch::hash::HasherSet<sketch::hash::WangHash> > cms_obj(5,10,137,false);   // np = 5, nh = 10, seed = 137,
    sketch::ocm::ccmbase <int32_t, sketch::hash::HasherSet<sketch::hash::WangHash> > ccms_obj(5,10,137,true);

    for(int i=0;i<kmers.size(); i++){
        cms_obj.addh_val(kmers[i]);
        ccms_obj.addh_val(kmers[i]);
    }




    //construct OCM
    sketch::ocm::ocmbase <uint64_t, sketch::hash::HasherSet<sketch::hash::WangHash> > sketch1(5,10);
    for(int r = 1; r<= 3; r++){
        if (r > 1){
            // for all kmers update collision
            for(auto kmer: kmers) sketch1.update_collision(kmer, r);
        }
        sketch1.clear_core();
        // for all kmer update count.
        for(auto kmer: kmers) sketch1.addh_val(kmer);
    }
    //end construct OCM


    for(int i=0;i<kmers.size(); i++) cout<<"kmer:    "<<kmers[i]<<"\tCMS_est:  "<<cms_obj.est_count(kmers[i])<<"\t Conservative Estimate: "<<ccms_obj.est_count(kmers[i])<<"\t Offline Cont: "<<sketch1.est_count(kmers[i]) <<endl;


}

