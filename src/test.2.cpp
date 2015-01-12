#include <iostream>
#include <sstream>
using namespace std;
    

int main () {
            for (int residueNumber1 = 0; residueNumber1 < 7000; residueNumber1++) {
            stringstream ss3a(stringstream::in | stringstream::out);
            //cout<<"check 0 :"<<endl;
            ss3a.clear();
            //cout<<ss3a.str()<<endl;
            //cout <<residueNumber1<<" , "<<myLeontisWesthofBondRow.residue1Atom[0]<<endl;
            //cout<<residueNumber1<<"/"<<myLeontisWesthofBondRow.residue1Atom[0]<<endl;;
         
            //cout<<"check 1 :"<<residueNumber1<<endl;
            //ss3a<<"anything";
            //cout<<ss3a.str()<<endl;
            ss3a.str("");
            //ss3a<<residueNumber1;
            cout<<"after check 1"<<endl;
            //ss3a<<"/";
            //cout<<"/"<<endl;
            //ss3a<<myLeontisWesthofBondRow.residue1Atom[0];
            cout<<"check 1.5:"<<residueNumber1<<"/"<<"CBlah"<<endl;;
            ss3a<<residueNumber1<<"/"<<"Cblah";
            //cout<<"check 2 :"<<endl;
            cout<<"check 1.7 :"<<ss3a.str()<<endl;
            //cout<<"check 3 :"<<endl;

            stringstream ss4(stringstream::in | stringstream::out);//(""); 
            /// ss4<<" ";
         
            //cout<<ss4.str()<<endl;
            //cout<<"check 3.5 :"<<endl;
            ss4.str("");
            //cout<<"check 3.6 :"<<endl;
            ss4.clear();
            ss4<<residueNumber1<<"/"<<residueNumber1;
            //cout<<ss4.str()<<endl;
            //cout<<"check 4 :"<<endl;
          }
};
