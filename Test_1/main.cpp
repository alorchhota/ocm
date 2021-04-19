#include"include/ccm.h"
#include<iostream>
#include<fstream>
#include <stdio.h>
#include<time.h>
#include <ctime>



using namespace std;

void read_file(string file)
{
    std::ifstream input(file);
    if(!input.good()){
        cout<<"Can't open"<<endl;
    }

    std::string line, name, content;
    while( std::getline( input, line ).good() ){
        if( line.empty() || line[0] == '>' ){
            if( !name.empty() ){
                std::cout << name << " : " << content << std::endl;
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
        std::cout << name << " : " << content << std::endl;
    }



}
vector<int> my_vec;
vector<int> counter;
int my_random()
{
    string s="";
    string my_arr[]={"65", "66", "67"};
    string flag="not found";

    for(int i=0; i<3; i++)
    {
        s=s+my_arr[rand()%3];
    }

    int x=stoi(s);


    for(int i=0; i<my_vec.size(); i++)
    {
        if(my_vec[i]==x)
        {
            counter[i]=counter[i]+1;
            flag="found";
            break;
        }
    }
    if(flag=="not found")
    {
        my_vec.push_back(x);
        counter.push_back(1);
    }
    return x;
}

int main()
{
    srand(time(0));

    sketch::cm::ccmbase_t<sketch::update::Increment,sketch::cm::DefaultCompactVectorType,sketch::hash::WangHash,false> *c = new sketch::cm::ccmbase_t<sketch::update::Increment,sketch::cm::DefaultCompactVectorType,
    sketch::hash::WangHash,false>(5, 10);
    int A=65, B=66, C=67;
    for(int i=0; i<300; i++)
    {
        int test = my_random();
        cout<<i+1<<"ith kmer: "<<test<<endl;
        c->addh(test);
    }

    int a= rand()%my_vec.size();
    cout<<endl<<"kmer: "<<my_vec[a]<<endl;
    cout<<endl<<"Real count: "<<counter[a]<<endl;

    cout<<"Estimate count: "<<c->est_count(my_vec[a])<<endl;
    //read_file("rymv.sim.fa");




	return 0;
}
