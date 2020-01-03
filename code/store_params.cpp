// FUNCTION FOR STORING FERTILITY RATES

vector<vector<double>> store_rates(const string &name){
	string line;
    ifstream myfile;
    myfile.open (name.c_str());
	if(myfile.fail()){ cout<<"rates input didn't work"<<endl;}
	int cols,rows=0;
	std::getline(myfile, line);
	++rows;
	istringstream iss(line);
	vector<string> tokens;
	copy(istream_iterator<string>(iss),
     istream_iterator<string>(),
     back_inserter(tokens));
	cols=tokens.size();
	while (std::getline(myfile, line))
        ++rows;
	vector<vector<double>>  rates(rows,vector<double>(cols));
	myfile.clear() ;
	myfile.seekg(0, ios::beg) ;
	for(int i=0; i<rows; i++){
        for(int j=0;j<cols;j++) {
        myfile>> rates[i][j];
        }
	}
	myfile.close();
return rates;
}

