#include "ezETAProgressBar.hpp"

#ifdef WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

int main() {
	int i;
	int n = 30;
	ez::ezETAProgressBar eta(n);

	eta.start();
	for(i=0; i <= n; ++i, ++eta) {
	#ifdef WIN32
		Sleep(1000);
	#else
		sleep(1);
	#endif
	}

	n = 99999;
	eta.n = n;
	eta.reset();
	eta.start();
	for(i=0; i <= 5; ++i, ++eta) {
	#ifdef WIN32
		Sleep(1000);
	#else
		sleep(1);
	#endif
	}

	return 0;
}