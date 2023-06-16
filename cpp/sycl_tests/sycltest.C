#include <iostream>
#include <cmath>

#include <CL/sycl.hpp>


static int custom_device_selector(const sycl::device &d) {
       using namespace cl::sycl;
      int rating = 0;
      if (d.is_gpu() && (d.get_info<info::device::name>().find("Nvidia") != std::string::npos)) {
        rating = 3;
      } else if (d.is_gpu()) {
        rating = 2;
      } else if (d.is_cpu()) {
        rating = 1;
      }
      return rating;
   }

  static void list_devices(){
    for (auto platform : sycl::platform::get_platforms())
    {
        std::cout << "Platform: "
                  << platform.get_info<sycl::info::platform::name>()
                  << std::endl;

        for (auto device : platform.get_devices())
        {
            std::cout << "\tDevice: "
                      << device.get_info<sycl::info::device::name>()
                      << std::endl;
        }
    }
  }


using namespace sycl;
int main(int argc, char* argv[])
{
  list_devices();

  sycl::device preferred_device { custom_device_selector };

   queue q(preferred_device);
 // sycl::queue q;

  std::cout << "Using device: "
            << q.get_device().get_info<sycl::info::device::name>()
            << std::endl;;

  // Compute the first n_items values in a well known sequence
  constexpr int n_items = 16;
  int *items = sycl::malloc_shared<int>(n_items, q);
  q.parallel_for(sycl::range<1>(n_items), [items] (sycl::id<1> i) {
      double x1 = pow((1.0 + sqrt(5.0))/2, i);
      double x2 = pow((1.0 - sqrt(5.0))/2, i);
      items[i] = round((x1 - x2)/sqrt(5));
  }).wait();

  for(int i = 0 ; i < n_items ; ++i) {
    std::cout << items[i] << std::endl;
  }
  free(items, q);

  return 0;
}
