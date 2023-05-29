#ifndef HeterogeneousTest_ROCmWrapper_interface_DeviceAdditionWrapper_h
#define HeterogeneousTest_ROCmWrapper_interface_DeviceAdditionWrapper_h

#include <cstddef>

namespace cms::rocmtest {

  void wrapper_add_vectors_f(const float* __restrict__ in1,
                             const float* __restrict__ in2,
                             float* __restrict__ out,
                             size_t size);

  void wrapper_add_vectors_d(const double* __restrict__ in1,
                             const double* __restrict__ in2,
                             double* __restrict__ out,
                             size_t size);

}  // namespace cms::rocmtest

#endif  // HeterogeneousTest_ROCmWrapper_interface_DeviceAdditionWrapper_h
