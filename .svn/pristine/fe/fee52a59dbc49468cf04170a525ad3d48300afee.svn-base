#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

int  main()
{
	
	ifstream number("../include/resources/torqueThetaSI.csv");

	vector<vector<double> > torqueTheta;	// store data
	vector<double> gridRow;	// row of numberGrid

	const int maxChars =100;
	int row;
	int col;
	double data;
    char * v;
	v=new char[maxChars];
	stringstream u;
	string stringData;
	
	if (number.is_open())
	{
		for (row = 0; row < 360; row++)
		{
			gridRow.clear();
			for (col = 0; col < 360; col++)
			{
                                cout<<"reading column: "<<col<<endl;
				number.getline(v,maxChars,','); 
                                cout<<"v =   "<<v  <<endl;
                                
				u.clear();
				u.str("");
				u<<string(v);
                                cout<<"read    data: "<<u   <<endl;
				u>>stringData;
				data = atof(stringData.c_str());
                                cout<<"pushing data: "<<data<<endl;
                 
				gridRow.push_back(data);
			}
			torqueTheta.push_back(gridRow);
		}
	}
	
	char run = 'y';
	int getRow, getCol;
	
	cout << "loadTorqueTheta.cpp" << endl;

	while(run == 'y')
	{
		cout << "Enter indices to retrieve:  ";
		cin >> getRow >> getCol;
		cout << "row = " << getRow << "and col = " << getCol;
		cout << torqueTheta[getRow][getCol] << endl;
		
		cout << "Would you like to do another calculation? (y/n)  ";
		cin >> run;
		cout << "run:  " << run << endl;
	}
    return 0;
}
