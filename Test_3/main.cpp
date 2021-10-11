#include<bits/stdc++.h>
#include"include/ocm.h"
using namespace std;


int main(){

    sketch::ocm::ccmbase <int32_t, sketch::hash::HasherSet<sketch::hash::WangHash> >c1(5,10,137,true);   // np = 5, nh = 10, seed = 137,

    //adding for test
    for(int i=0; i<90; i++){
        c1.addh_val(108);
    }

    cout<<c1.est_count(108)<<endl;


//    //read_file("rymv.sim.fa",len_k_mer);
}
