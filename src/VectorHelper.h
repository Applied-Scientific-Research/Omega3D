#pragma once

#include <vector>

#ifdef USE_VC
#include <Vc/Vc>

// Must use Vc's allocator to ensure memory alignment
template <class S> using Vector = std::vector<S, Vc::Allocator<S>>;
// now we can use Vector<float> in code

//using VectorF = std::vector<float, Vc::Allocator<float>>;
//using VectorD = std::vector<double, Vc::Allocator<double>>;

// here's a way to make other vectors with the same number of entries as float_v
//using Vc::float_v;
//typedef Vc::SimdArray<double, float_v::size()> double_v;
//typedef Vc::SimdArray<std::uint16_t, float_v::size()> uint16_v;
//typedef Vc::SimdArray<std::uint32_t, float_v::size()> uint32_v;

// use this to safely convert a std::vector<S> to an array of Vc::Vector<S>
template <class S>
inline const Vc::Memory<Vc::Vector<S>> stdvec_to_vcvec (const Vector<S>& in, const S defaultval) {
    Vc::Memory<Vc::Vector<S>> out(in.size());
    for (size_t i=0; i<in.size(); ++i) out[i] = in[i];
    for (size_t i=in.size(); i<out.entriesCount(); ++i) out[i] = defaultval;
    return out;
}

#else	// no Vc present, use stdlib instead

template <class S> using Vector = std::vector<S>;

//using VectorF = std::vector<float>;
//using VectorD = std::vector<double>;

#endif
