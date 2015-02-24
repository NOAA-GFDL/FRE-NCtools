#include "ezOdometer.hpp"
#include <assert.h>

int main() {
  ez::ezOdometer<int> o1(3);
  o1.setrange(0,0,1);
  o1.setrange(1,0,1);
  o1.setrange(2,0,1);
  std::cout << "o1="; o1.print(); std::cout << "\n";
  ++o1;
  assert( o1.counters[0]==0 );
  assert( o1.counters[1]==0 );
  assert( o1.counters[2]==1 );
  
  while(!o1.isDone()) {
    std::cout << "o1="; o1.print(); std::cout << "\n";
    ++o1;
  };
  assert( o1.counters[0]==1 );
  assert( o1.counters[1]==1 );
  assert( o1.counters[2]==1 );
  
  std::cout << "\n";
  o1.reset();
  o1.lock(1, 0);
  while(!o1.isDone()) {
    std::cout << "o1="; o1.print(); std::cout << "\n";
    ++o1;
  };
  assert( o1.counters[0]==1 );
  assert( o1.counters[1]==0 );
  assert( o1.counters[2]==1 );

  std::cout << "\n";
  o1.reset();
  o1.lock(1, 0);
  ++o1;
  assert( o1.atLimit(0) == false );
  assert( o1.atLimit(2) );
  
  while(!o1.isDone()) {
    std::cout << "o1="; o1.print(); std::cout << "\n";
    ++o1;
  };
  assert( o1.counters[0]==1 );
  assert( o1.counters[1]==0 );
  assert( o1.counters[2]==1 );
  assert( o1.isDone() );
  assert( o1.atLimit(0) );
  assert( o1.atLimit(2) );
  std::cout << "o1="; o1.print(); std::cout << "\n";

  std::cout << "\n";
  o1.reset();
  o1.setrange(0,-1,1);
  o1.setrange(2,2,10,3);
  assert( o1.counters[0]==-1 );
  assert( o1.counters[1]==0 );
  assert( o1.counters[2]==2 );
  while(!o1.isDone()) {
    std::cout << "o1="; o1.print(','); std::cout << "\n";
    ++o1;
  };
  assert( o1.counters[0]==1 );
  assert( o1.counters[1]==0 );
  assert( o1.counters[2]==8 );
  std::cout << "o1="; o1.print(','); std::cout << "\n";

  /// Test custom lists.
  ez::ezCustomOdometer<int> o2(4);
  o2.setrange(0,-1,5,3);
  std::vector<int> values;
  values.push_back(-5);
  values.push_back(10);
  values.push_back(123);
  o2.setCustomValues(1,values);
  values.clear();
  values.push_back(9);
  values.push_back(99);  
  o2.setCustomValues(2,values);
  o2.setCustomValue(3,2);

  assert( o2.digits[0]==-1 );
  assert( o2.digits[1]==-5 );
  assert( o2.digits[2]==9 );
  assert( o2.digits[3]==2 );

  std::cout << "\n";
  std::cout << "o2="; o2.print(','); std::cout << "\n";
  ++o2;
  assert( o2.digits[0]==-1 );
  assert( o2.digits[1]==-5 );
  assert( o2.digits[2]==99 );
  assert( o2.digits[3]==2 );
  
  while(!o2.isDone()) {
    std::cout << "o2="; o2.print(','); std::cout << "\n";
    ++o2;
  };
  assert( o2.digits[0]==5 );
  assert( o2.digits[1]==123 );
  assert( o2.digits[2]==99 );
  assert( o2.digits[3]==2 );
  std::cout << "o2="; o2.print(','); std::cout << "\n";
  
  /// Test custom lists with a lock.
  o2.clear();
  o2.init(4);
  o2.setrange(0,-1,5,3);
  values.clear();
  values.push_back(-5);
  values.push_back(10);
  values.push_back(123);
  o2.setCustomValues(1,values);
  o2.setCustomValue(2,1);
  o2.setCustomValue(3,2);

  std::cout << "\n";
  std::cout << "o2="; o2.print(','); std::cout << "\n";

  ++o2;
  assert( o2.atLimit(0) == false );
  assert( o2.atLimit(3) );
  
  while(!o2.isDone()) {
    std::cout << "o2="; o2.print(','); std::cout << "\n";
    ++o2;
  };
  assert( o2.digits[0]==5 );
  assert( o2.digits[1]==123 );
  assert( o2.digits[2]==1 );
  assert( o2.digits[3]==2 );
  assert( o2.isDone() );
  assert( o2.atLimit(0) );
  assert( o2.atLimit(2) );
  std::cout << "o2="; o2.print(','); std::cout << "\n";

  /// Test simple odometer normalized.
  o1.clear();
  o1.init(3);
  o1.setrange(0,-1,5,3);
  o1.setrange(1,5,10,2);
  o1.setrange(2,-2,0);
  assert( o1.counters[0]==-1 );
  assert( o1.counters[1]==5 );
  assert( o1.counters[2]==-2 );
  
  o1.normalize();
  assert( o1.counters[0]==0 );
  assert( o1.counters[1]==0 );
  assert( o1.counters[2]==0 );
  assert( o1.maxs[0]==2 );
  assert( o1.maxs[1]==2 );
  assert( o1.maxs[2]==2 );
  
  while(!o1.isDone()) {
    std::cout << "o1="; o1.print(','); std::cout << "\n";
    ++o1;
  };
  assert( o1.counters[0]==2 );
  assert( o1.counters[1]==2 );
  assert( o1.counters[2]==2 );
  std::cout << "o1="; o1.print(','); std::cout << "\n";

  /// Test custom lists that are then normalized.
  o2.clear();
  o2.init(4);
  o2.setrange(0,-1,5,3);
  values.clear();
  values.push_back(-5);
  values.push_back(10);
  values.push_back(123);
  o2.setCustomValues(1,values);
  o2.setCustomValue(2,1);
  o2.setCustomValue(3,2);
  
  ez::ezCustomOdometer<int> o3(4);
  o3 = o2;
  assert( o3.digits[0]==-1 );
  assert( o3.digits[1]==-5 );
  assert( o3.digits[2]==1 );
  assert( o3.digits[3]==2 );

  o3.normalize();
  assert( o3.digits[0]==0 );
  assert( o3.digits[1]==0 );
  assert( o3.digits[2]==0 );
  assert( o3.digits[3]==0 );
  assert( o3.maxs[0]==2 );
  assert( o3.maxs[1]==2 );
  assert( o3.maxs[2]==0 );
  assert( o3.maxs[3]==0 );
  
  while(!o3.isDone()) {
    std::cout << "o3="; o3.print(','); std::cout << "\n";
    ++o3;
  };
  assert( o3.digits[0]==2 );
  assert( o3.digits[1]==2 );
  assert( o3.digits[2]==0 );
  assert( o3.digits[3]==0 );
  assert( o3.isDone() );
  std::cout << "o3="; o3.print(','); std::cout << "\n";

  // Number of combinations reported.
  o1.clear();
  assert( o1.getNumCombinations() == 0 );
  
  o1.clear();
  o1.init(1);
  o1.setrange(0,0,0);
  assert( o1.getNumCombinations() == 1 );

  o1.clear();
  o1.init(3);
  o1.setrange(0,0,1);
  o1.setrange(1,0,1);
  o1.setrange(2,0,1);
  assert( o1.getNumCombinations() == 2*2*2 );

  o1.clear();
  o1.init(3);
  o1.setrange(0,0,1);
  o1.setrange(1,0,0);
  o1.setrange(2,0,1);
  assert( o1.getNumCombinations() == 2*1*2 );

  o1.clear();
  o1.init(3);
  o1.setrange(0,1,10);
  o1.setrange(1,1,10,2);
  o1.setrange(2,-5,-1,3);
  assert( o1.getNumCombinations() == 10*5*2 );

  o2.clear();
  o2.init(4);
  o2.setrange(0,-1,5,3);
  values.clear();
  values.push_back(-5);
  values.push_back(10);
  values.push_back(123);
  o2.setCustomValues(1,values);
  o2.setCustomValue(2,1);
  o2.setCustomValue(3,2);
  assert( o2.getNumCombinations() == 3*3*1*1 );

  // Reversed sequences.
  o1.clear();
  o1.init(3);
  o1.setrange(0,10,8,-1);
  o1.setrange(1,10,5,-2);
  o1.setrange(2,-1,-5,-3);
  assert( o1.getNumCombinations() == 3*3*2 );
  assert( o1.counters[0]==10 );
  assert( o1.counters[1]==10 );
  assert( o1.counters[2]==-1 );

  while(!o1.isDone()) {
    std::cout << "o1="; o1.print(','); std::cout << "\n";
    ++o1;
  };
  assert( o1.counters[0]==8 );
  assert( o1.counters[1]==6 );
  assert( o1.counters[2]==-4 );

  o3.clear();
  o3 = o2;
  o3.normalize();
  assert( o3.getNumCombinations() == 3*3*1*1 );

  return 0;
}