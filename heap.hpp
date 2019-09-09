#include <algorithm>

#if __cplusplus >= 202000
#include <version>
#else

#include <ciso646>

#endif

namespace std
{

#if defined(_LIBCPP_VERSION)

    template<class _RandomAccessIterator, class _Compare>
    inline _LIBCPP_INLINE_VISIBILITY
    void
    update_heap(_RandomAccessIterator __first, _RandomAccessIterator __last, _RandomAccessIterator __elem,
                _Compare __comp)
    {
        typedef typename __comp_ref_type<_Compare>::type _Comp_ref;

        if((__first + 1) >= __last)
        {
            return;
        }

        if(__elem == __first)
        {
            return __sift_down<_Comp_ref>(__first, __last, __comp, __last - __first, __elem);
        }

        typename iterator_traits<_RandomAccessIterator>::difference_type __parent = (__elem - __first - 1) / 2;

        if(__comp(*(__first + __parent), *__elem))
        {
            return __sift_up<_Comp_ref>(__first, __elem + 1, __comp, __elem + 1 - __first);
        }
        else
        {
            return __sift_down<_Comp_ref>(__first, __last, __comp, __last - __first, __elem);
        }
    }

    template<class _RandomAccessIterator>
    inline _LIBCPP_INLINE_VISIBILITY
    void
    update_heap(_RandomAccessIterator __first, _RandomAccessIterator __last, _RandomAccessIterator __elem)
    {
        _VSTD::update_heap(__first, __last, __elem,
                           __less<typename iterator_traits<_RandomAccessIterator>::value_type>());
    }

#elif defined(__GLIBCXX__) || (__GLIBCPP__)

    /**
     *  @brief  Update an element on a heap.
     *  @param  __first  Start of heap.
     *  @param  __last   End of heap.
     *  @param  __elem   Updated element.
     *  @ingroup heap_algorithms
     *
     *  This operation potentially updates the position of
     *  the element to restore [__first,__last) to be a valid
     *  heap. The element must be in [__first,__last).
     */
    template<typename _RandomAccessIterator>
    inline void
    update_heap(_RandomAccessIterator __first, _RandomAccessIterator __last, _RandomAccessIterator __elem)
    {
        typedef typename iterator_traits<_RandomAccessIterator>::value_type _ValueType;
        typedef typename iterator_traits<_RandomAccessIterator>::difference_type _DistanceType;

        // concept requirements
        __glibcxx_function_requires(_Mutable_RandomAccessIteratorConcept<_RandomAccessIterator>)
        __glibcxx_function_requires(_LessThanComparableConcept<_ValueType>)
        __glibcxx_requires_non_empty_range(__first, __last);
        __glibcxx_requires_valid_range(__first, __last);
        __glibcxx_requires_valid_range(__first, __elem);
        __glibcxx_requires_non_empty_range(__elem, __last);
        __glibcxx_requires_valid_range(__elem, __last);
        __glibcxx_requires_irreflexive(__first, __last);
        __glibcxx_requires_heap(__first, __last);

        if ((__first + 1) >= __last)
        {
            return;
        }

        __gnu_cxx::__ops::_Iter_less_iter __compIter;
        __gnu_cxx::__ops::_Iter_less_val __compVal;

        _ValueType __value = _GLIBCXX_MOVE(*__elem);

        if (__elem == __first)
        {
            return std::__adjust_heap(__first, _DistanceType(0), _DistanceType(__last - __first), _GLIBCXX_MOVE(__value), __compIter);
        }

        _DistanceType __parent = (__elem - __first - 1) / 2;

        if (__compVal(__first + __parent, __value))
        {
            return std::__push_heap(__first, _DistanceType(__elem - __first), _DistanceType(0), _GLIBCXX_MOVE(__value), __compVal);
        } else {
            return std::__adjust_heap(__first, _DistanceType(__elem - __first), _DistanceType(__last - __first), _GLIBCXX_MOVE(__value), __compIter);
        }
    }

    /**
     *  @brief  Update an element on a heap.
     *  @param  __first  Start of heap.
     *  @param  __last   End of heap.
     *  @param  __elem   Updated element.
     *  @param  __comp   Comparison functor.
     *  @ingroup heap_algorithms
     *
     *  This operation potentially updates the position of
     *  the element to restore [__first,__last) to be a valid
     *  heap. The element must be in [__first,__last).
     *  Compare operations are performed using comp.
     */
    template<typename _RandomAccessIterator, typename _Compare>
    inline void
    update_heap(_RandomAccessIterator __first, _RandomAccessIterator __last, _RandomAccessIterator __elem, _Compare __comp)
    {
        typedef typename iterator_traits<_RandomAccessIterator>::value_type _ValueType;
        typedef typename iterator_traits<_RandomAccessIterator>::difference_type _DistanceType;

        // concept requirements
        __glibcxx_function_requires(_Mutable_RandomAccessIteratorConcept<_RandomAccessIterator>)
        __glibcxx_requires_non_empty_range(__first, __last);
        __glibcxx_requires_valid_range(__first, __last);
        __glibcxx_requires_valid_range(__first, __elem);
        __glibcxx_requires_non_empty_range(__elem, __last);
        __glibcxx_requires_valid_range(__elem, __last);
        __glibcxx_requires_irreflexive_pred(__first, __last, __comp);
        __glibcxx_requires_heap_pred(__first, __last, __comp);

        if ((__first + 1) >= __last)
        {
            return;
        }

        __decltype(__gnu_cxx::__ops::__iter_comp_iter(_GLIBCXX_MOVE(__comp))) __compIter(_GLIBCXX_MOVE(__comp));
        __decltype(__gnu_cxx::__ops::__iter_comp_val(_GLIBCXX_MOVE(__comp))) __compVal(_GLIBCXX_MOVE(__comp));

        _ValueType __value = _GLIBCXX_MOVE(*__elem);

        if (__elem == __first)
        {
            return std::__adjust_heap(__first, _DistanceType(0), _DistanceType(__last - __first), _GLIBCXX_MOVE(__value), __compIter);
        }

        _DistanceType __parent = (__elem - __first - 1) / 2;

        if (__compVal((__first + __parent), __value))
        {
            return std::__push_heap(__first, _DistanceType(__elem - __first), _DistanceType(0), _GLIBCXX_MOVE(__value), __compVal);
        } else {
            return std::__adjust_heap(__first, _DistanceType(__elem - __first), _DistanceType(__last - __first), _GLIBCXX_MOVE(__value), __compIter);
        }
    }

#elif defined(_CPPLIB_VER) || defined(__INTEL_CXXLIB_ICC)

    template<class _RanIt, class _Pr>
    inline void
    update_heap(_RanIt _First, _RanIt _Last, _RanIt _Elem, _Pr _Pred)
    {
        typedef _Iter_diff_t<_RanIt> _Diff;

        // Commented out since debug macros are not stable across versions in dinumkware
        //_Adl_verify_range(_First, _Last);
        //_Adl_verify_range(_Elem, _Last);

        if ((_First + 1) >= _Last)
        {
            return;
        }

        _Iter_value_t<_RanIt> _Val = _STD move(*_Elem);

        if (_Elem == _First)
        {
            return _Pop_heap_hole_by_index(_First, _Diff(0), _Diff(_Last - _First), _STD move(_Val), _Pass_fn(_Pred));
        }

        _Diff _Parent = (_Elem - _First - 1) / 2;

        if (_DEBUG_LT_PRED(_Pred, *(_First + _Parent), _Val))
        {
            return _Push_heap_by_index(_First, _Diff(_Elem - _First), _Diff(0), _STD move(_Val), _Pass_fn(_Pred));
        } else {
            return _Pop_heap_hole_by_index(_First, _Diff(_Elem - _First), _Diff(_Last - _First), _STD move(_Val), _Pass_fn(_Pred));
        }
    }

    template<class _RanIt>
    inline void
    update_heap(_RanIt _First, _RanIt _Last, _RanIt _Elem)
    {
        _STD update_heap(_First, _Last, _Elem, less<>());
    }

#else

#error "You are using a currently unsupported stdlib (libc++ (clang, icc), libstdc++ (gcc, icc) and dinkumware (MSVC, icc) are supported). Please contact the project maintainer."

#endif

}
