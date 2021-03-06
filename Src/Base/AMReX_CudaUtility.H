#ifndef AMREX_CUDA_UTILITY_H_
#define AMREX_CUDA_UTILITY_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_CudaDevice.H>
#include <AMReX_Extension.H>
#include <AMReX_REAL.H>
#include <iostream>
#include <cmath>

namespace amrex {
namespace Cuda {

    inline bool isManaged (void const* p) {
#ifdef AMREX_USE_CUDA
        cudaPointerAttributes ptr_attr;
        cudaPointerGetAttributes(&ptr_attr, p);
        cudaError_t err = cudaGetLastError();
        return (err != cudaErrorInvalidValue) && ptr_attr.isManaged;
#else
        return false;
#endif
    }

    inline bool isDevicePtr (void const* p) {
#ifdef AMREX_USE_CUDA
        cudaPointerAttributes ptr_attr;
        cudaPointerGetAttributes(&ptr_attr, p);
        cudaError_t err = cudaGetLastError();
        return (err != cudaErrorInvalidValue) && (ptr_attr.memoryType == cudaMemoryTypeDevice);
#else
        return false;
#endif
    }

    inline bool isHostPtr (void const* p) {
#ifdef AMREX_USE_CUDA
        cudaPointerAttributes ptr_attr;
        cudaPointerGetAttributes(&ptr_attr, p);
        cudaError_t err = cudaGetLastError();
        return (err != cudaErrorInvalidValue) && (ptr_attr.memoryType == cudaMemoryTypeHost);
#else
        return false;
#endif
    }

    inline bool isCudaPtr (void const* p) {
#ifdef AMREX_USE_CUDA
        cudaPointerAttributes ptr_attr;
        cudaPointerGetAttributes(&ptr_attr, p);
        cudaError_t err = cudaGetLastError();
        return (err != cudaErrorInvalidValue);
#else
        return false;
#endif
    }

    namespace Atomic {

#ifdef AMREX_USE_CUDA
        namespace detail {

            AMREX_GPU_DEVICE AMREX_INLINE
            float atomicMax(float* address, float val)
            {
                int* address_as_i = (int*) address;
                int old = *address_as_i, assumed;
                do {
                    assumed = old;
                    old = atomicCAS(address_as_i, assumed,
                                    __float_as_int(fmaxf(val, __int_as_float(assumed))));
                } while (assumed != old);
                return __int_as_float(old);
            }
            
            AMREX_GPU_DEVICE AMREX_INLINE
            double atomicMax(double* address, double val)
            {
                unsigned long long int* address_as_ull = 
                    (unsigned long long int*) address;
                unsigned long long int old = *address_as_ull, assumed;
                do {
                    assumed = old;
                    old = atomicCAS(address_as_ull, assumed,
                                    __double_as_longlong(fmax(val, __longlong_as_double(assumed))));
                } while (assumed != old);
                return __longlong_as_double(old);
            }

            AMREX_GPU_DEVICE AMREX_INLINE
            float atomicMin(float* address, float val)
            {
                int* address_as_i = (int*) address;
                int old = *address_as_i, assumed;
                do {
                    assumed = old;
                    old = atomicCAS(address_as_i, assumed,
                                    __float_as_int(fminf(val, __int_as_float(assumed))));
                } while (assumed != old);
                return __int_as_float(old);
            }

            AMREX_GPU_DEVICE AMREX_INLINE
            double atomicMin(double* address, double val)
            {
                unsigned long long int* address_as_ull = 
                    (unsigned long long int*) address;
                unsigned long long int old = *address_as_ull, assumed;
                do {
                    assumed = old;
                    old = atomicCAS(address_as_ull, assumed,
                                    __double_as_longlong(fmin(val, __longlong_as_double(assumed))));
                } while (assumed != old);
                return __longlong_as_double(old);
            }
        } // namespace detail
#endif  // AMREX_USE_CUDA

        template<class T>
        AMREX_GPU_HOST_DEVICE AMREX_INLINE
        void Add (T* sum, T value)
        {
#if defined(__CUDA_ARCH__)
            atomicAdd(sum, value);
#else
            *sum += value;
#endif
        }
        
        template<class T>
        AMREX_GPU_HOST_DEVICE AMREX_INLINE
        void Min (T* m, T value)
        {
#if defined(__CUDA_ARCH__)
            detail::atomicMin(m, value);
#else
            *m = (*m) < value ? (*m) : value;
#endif
        }
        
        template<class T>
        AMREX_GPU_HOST_DEVICE AMREX_INLINE
        void Max (T* m, T value)
        {
#if defined(__CUDA_ARCH__)
            detail::atomicMax(m, value);
#else
            *m = (*m) > value ? (*m) : value;
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_INLINE
        void Or (int* m, int value)
        {
#if defined(__CUDA_ARCH__)
            atomicOr(m, value);
#else
            *m = (*m) || value; 
#endif
        }

    } // namespace Atomic

    template <class T>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    bool isnan(T m)
    {
#if defined(__CUDA_ARCH__)
        return isnan(m);
#else
        return std::isnan(m);
#endif
    }

    template <class T>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    bool isinf(T m)
    {
#ifdef __CUDA_ARCH__
        return isinf(m);
#else
        return std::isinf(m);
#endif
    }

    class StreamIter
    {
    public:
        StreamIter (const int n, bool is_thread_safe=true);
        ~StreamIter ();

        StreamIter (StreamIter const&) = delete;
        StreamIter (StreamIter &&) = delete;
        void operator= (StreamIter const&) = delete;
        void operator= (StreamIter &&) = delete;

        int operator() () const { return m_i; }

        bool isValid () const { return m_i < m_n; }

#if !defined(AMREX_USE_CUDA)
        void operator++ () { ++m_i; }
#else
        void operator++ ();
#endif

    private:
        int m_n;
        int m_i;
        bool m_threadsafe;
    };

} // namespace Cuda

#if AMREX_USE_CUDA
std::ostream& operator<< (std::ostream& os, const dim3& d);
#endif

using Cuda::isnan;
using Cuda::isinf;

} // namespace amrex

#endif
