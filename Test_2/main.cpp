#include"include/ccm.h"
#include<iostream>
#include<fstream>
#include <stdio.h>
#include<time.h>
#include <ctime>
#define ccm_flag true
#define cm_flag false



using namespace std;
void my_binary(int64_t n) //converting k-mer to binary representation
{
    int a[100];
    string str="";
    int i;
    for(i=0; i<64; i++)
    {
        if(n>0)
        {
            a[i]=n%2;
            n= n/2;
        }
        else a[i]=0;

    }
    cout<<endl<<"Binary: ";
    for(int j=i-1 ;j>=0 ;j--)
    {
        cout<<a[j];
        //str+=to_string(a[i]);
    }
    cout<<endl;
}
int64_t cal(string str_k_mer)
{
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
    //cout<<"cal value: "<<k_mer;
    //my_binary(k_mer);
    return k_mer;

}
void read_file(string file, int len_k_mer)
{
    string sample="TGAAATTCCTGGGTGCCTCCAA";
    int counter=0;
    size_t x =10;
    sketch::cm::cs4wbase_t<> c2(16, 10); //or csbase_t
    sketch::cm::cs4wbase_t<>&& c = std::move(c2);
    sketch::cm::SlidingWindow<sketch::cm::cs4w_t> sw(x, std::move(c));

    std::ifstream input(file);
    if(!input.good()){
        cout<<"Can't open"<<endl;
    }

    std::string line, name, content;
    while( std::getline( input, line ).good()){ //reading 1 line from input

        if( line.empty() || line[0] == '>' ){
            if( !name.empty() ){
                int l=content.length();
                for(int i=0; i<=l-len_k_mer; i++) //extracting k-mers from the line
                {
                    //int64_t k_mer = 0;

                    string part = content.substr(i, len_k_mer);
                    int64_t k_mer = cal(part);

                    std::cout<<k_mer  << " : " << content.substr(i, len_k_mer) <<" : "<< std::endl;
                    my_binary(k_mer);
                    //std::cout << name << " : " << content << std::endl;


                    sw.addh(k_mer);

                    if(k_mer == cal(sample))
                    {
                        counter++;

                    }
                }
                //std::cout << name << " : " << content << std::endl;
                name.clear();
            }
            if( !line.empty() ){
                name = line.substr(1);
            }
            content.clear();
        } else if( !name.empty() ){
            if( line.find(' ') != std::string::npos ){
                name.clear();
                content.clear();
            } else {
                content += line;
            }
        }
    }
    if( !name.empty() ){
        //std::cout << name << " : " << content << std::endl;
        int l=content.length();
        for(int i=0; i<=l-len_k_mer; i++)
        {
            //int64_t k_mer = 0;

            string part = content.substr(i, len_k_mer);
            int64_t k_mer = cal(part);

            std::cout<<k_mer  << " : " << content.substr(i, len_k_mer) <<" : "<< std::endl;
            my_binary(k_mer); //binary representation of the k-mers
            //std::cout << name << " : " << content << std::endl;


            sw.addh(k_mer); //inserting the k-mers

            if(k_mer == cal(sample))
            {
                counter++; //manual count of the sample k-mer

            }
        }
}
    cout<<"Estimate count: "<<sw.sketch().est_count(cal(sample))<<endl; //estimation of the sample k-mer
    //cout<<"Estimate count: "<<c2.est_count(cal(sample))<<endl;
    cout<<counter<<endl;



}


//rand generator function
//not relevant now
//vector<int> my_vec;
//vector<int> counter;
//int my_random()
//{
//    string s="";
//    string my_arr[]={"65", "66", "67"};
//    string flag="not found";
//
//    for(int i=0; i<3; i++)
//    {
//        s=s+my_arr[rand()%3];
//    }
//
//    int x=stoi(s);
//
//
//    for(int i=0; i<my_vec.size(); i++)
//    {
//        if(my_vec[i]==x)
//        {
//            counter[i]=counter[i]+1;
//            flag="found";
//            break;
//        }
//    }
//    if(flag=="not found")
//    {
//        my_vec.push_back(x);
//        counter.push_back(1);
//    }
//    return x;
//}

int main()
{
//    //srand(time(0));
    int len_k_mer =22;
//    string sample = "GCGTCGCTGTGGAGCGAGCCTG"; //74
//    //bool ccm_flag = true, cm_flag=false;
//
//
//    //sketch::cm::ccmbase_t<sketch::update::Increment,sketch::cm::DefaultCompactVectorType,sketch::hash::WangHash,false> *c = new sketch::cm::ccmbase_t<sketch::update::Increment,sketch::cm::DefaultCompactVectorType,
//    //sketch::hash::WangHash,false>(4, 10);
////    int A=65, B=66, C=67;
////    for(int i=0; i<300; i++)
////    {
////        int test = my_random();
////        cout<<i+1<<"ith kmer: "<<test<<endl;
////        c->addh(test);
////    }
////
////    int a= rand()%my_vec.size();
////    cout<<endl<<"kmer: "<<my_vec[a]<<endl;
////    cout<<endl<<"Real count: "<<counter[a]<<endl;
////
//
//    //read_file("rymv.sim.fa",len_k_mer);
//    //read_file("test.fa", 2, c);
//    //cout<<"Estimate count: "<<c->est_count(cal("GCGT"))<<endl;
////
//
//
//    //cout<<"Estimate count: "<<c->est_count(185)<<endl;



//using ccmbase_t checking different values of nbits, nhashes...only nbits=5 gives correct estimate

//    int number=50;
//
//    for(unsigned nhashes=10; nhashes<=25; nhashes+=5)
//    {
//
//        for(unsigned nbits=64; nbits<100; nbits++)
//        {
//            sketch::cm::ccmbase_t<sketch::update::Increment,sketch::cm::DefaultCompactVectorType,sketch::hash::WangHash,false> *c = new sketch::cm::ccmbase_t<sketch::update::Increment,sketch::cm::DefaultCompactVectorType,
//            sketch::hash::WangHash,false>(nbits, nhashes);
//
//            for(int i=0; i<100; i++)
//            {
//                c->addh_val(number);
//
//            }
//            cout<<"nbits: "<<nbits<<"nhashes: "<<nhashes<<"\tEstimate count: "<<c->est_count(number)<<endl;
//        }
//
//    }



//using cs4wbase_t and sliding window
//inserting a demo number to test
    size_t x =50000;
    sketch::cm::cs4wbase_t<> c2(5, 10);
    sketch::cm::cs4wbase_t<>&& c = std::move(c2);

    sketch::cm::SlidingWindow<sketch::cm::cs4w_t> sw(x, std::move(c));
    for(int i=0; i<67; i++)
    {
        sw.addh(108);
    }

    cout<<sw.sketch().est_count(108)<<endl;



//using cs4wbase_t and sliding window
//reading original k-mers from the file and estimating the count
    read_file("rymv.sim.fa",len_k_mer);







	return 0;
}
