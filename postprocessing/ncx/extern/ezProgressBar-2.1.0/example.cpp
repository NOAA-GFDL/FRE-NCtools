/* 
g++ example.cpp -o example

CHANGELOG
20110512 rsz Created.
*/

#include "ezProgressBar.hpp"

#ifdef WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

int main() {
	int n = 100;
	ez::ezProgressBar p1(n);
	int i;

	p1.start();
	for(i=0; i <= n; ) {
	#ifdef WIN32
		Sleep(1000);
	#else
		sleep(1);
	#endif
		++p1;
		++p1;
		i += 2;
	}

	std::cout << std::endl;
	
	p1.reset();
	n = 20;
	p1.n = n;
	p1.start();
	for(i=0; i <= n; ++i) {
	#ifdef WIN32
		Sleep(1000);
	#else
		sleep(1);
	#endif
		++p1;
	}

	std::cout << std::endl;

	return 0;
}