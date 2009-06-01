LDFLAGS  =   -L /usr/local/lib/ -lcondor  -lm
CXXFLAGS = -I /usr/local/include -I ..


t-condor : t-condor.cc
	$(CXX) $(CXXFLAGS) -o t-condor t-condor.cc $(LDFLAGS)  
	
