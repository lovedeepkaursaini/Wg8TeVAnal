CPPFLAGS=`root-config --cflags` -I/usr/local/root/include
LDFLAGS=`root-config --libs` -lMinuit
skimmer_wg: skimmer_wg.cc
	g++ $(CPPFLAGS) $(LDFLAGS) -o $@ $^

