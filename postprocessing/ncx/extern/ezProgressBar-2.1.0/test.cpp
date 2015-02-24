/* 
CHANGELOG
20110522 rsz Created.
*/

#include "ezProgressBar.hpp"
#include "ezETAProgressBar.hpp"
#include <assert.h>

int main() {
	ez::ezProgressBar p(3);
	p.start();
	while(p.cur != p.n) ++p;
	std::cout << std::endl;

	ez::ezETAProgressBar eta(3);
	// seconds = d*24*60*60 + h*60*60 + m*60 + s
	// d=8,h=0,m=0,s=0
	//std::cout << eta.secondsToString(691200) << std::endl;
	assert(eta.secondsToString(691200).compare("8d ") == 0);
	// d=10,h=0,m=0,s=0
	//std::cout << eta.secondsToString(864000) << std::endl;
	assert(eta.secondsToString(864000).compare("10d ") == 0);
	// d=123,h=0,m=0,s=0
	//std::cout << eta.secondsToString(10627200) << std::endl;
	assert(eta.secondsToString(10627200).compare("123d ") == 0);
	// d=1,h=0,m=0,s=0
	//std::cout << eta.secondsToString(86400) << std::endl;
	assert(eta.secondsToString(86400).compare("1d ") == 0);
	// d=0,h=23,m=0,s=0
	//std::cout << eta.secondsToString(82800) << std::endl;
	assert(eta.secondsToString(82800).compare("23h ") == 0);
	// d=0,h=12,m=0,s=0
	//std::cout << eta.secondsToString(43200) << std::endl;
	assert(eta.secondsToString(43200).compare("12h ") == 0);
	// d=0,h=1,m=0,s=0
	//std::cout << eta.secondsToString(3600) << std::endl;
	assert(eta.secondsToString(3600).compare("1h ") == 0);
	// d=0,h=0,m=41,s=0
	//std::cout << eta.secondsToString(2460) << std::endl;
	assert(eta.secondsToString(2460).compare("41m ") == 0);
	// d=0,h=0,m=1,s=0
	//std::cout << eta.secondsToString(60) << std::endl;
	assert(eta.secondsToString(60).compare("1m ") == 0);
	// d=0,h=0,m=0,s=29
	//std::cout << eta.secondsToString(29) << std::endl;
	assert(eta.secondsToString(29).compare("29s") == 0);
	// d=13,h=10,m=0,s=0
	//std::cout << eta.secondsToString(1159200) << std::endl;
	assert(eta.secondsToString(1159200).compare("13d 10h ") == 0);
	// d=13,h=10,m=56,s=0
	//std::cout << eta.secondsToString(1162560) << std::endl;
	assert(eta.secondsToString(1162560).compare("13d 10h 56m ") == 0);
	// d=12,h=34,m=56,s=43
	//std::cout << eta.secondsToString(1162603) << std::endl;
	assert(eta.secondsToString(1162603).compare("13d 10h 56m 43s") == 0);
	// d=0,h=13,m=56,s=43
	//std::cout << eta.secondsToString(50203) << std::endl;
	assert(eta.secondsToString(50203).compare("13h 56m 43s") == 0);
	// d=0,h=0,m=56,s=43
	//std::cout << eta.secondsToString(3403) << std::endl;
	assert(eta.secondsToString(3403).compare("56m 43s") == 0);
	
	eta.start();
	while(eta.cur != eta.n) ++eta;

	return 0;
}
