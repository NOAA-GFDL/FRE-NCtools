#include "ezRateProgressBar.hpp"

#ifdef WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

int main() {
	int i;
	int n = 1000;
	
	ez::ezRateProgressBar<int> rate(n);
	rate.units = "MB";
	rate.start();
	for(i=0; i <= n; i+=100) {
		rate.update(i);
	#ifdef WIN32
		Sleep(1000);
	#else
		sleep(1);
	#endif
	}

	std::cout << std::endl;
	
	rate.reset();
	n = 10;
	int delta = 10000;
	rate.n = delta*n;
	rate.start();
	for(i=0; i <= n; ++i) {
	#ifdef WIN32
		Sleep(1000);
	#else
		sleep(1);
	#endif
		rate+=delta;
	}

	return 0;
}