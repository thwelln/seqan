// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Implementation of the Packed String class.
// ==========================================================================

#ifndef SEQAN_SEQUENCE_STRING_PACKED_H_
#define SEQAN_SEQUENCE_STRING_PACKED_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename T>
struct HostIterator;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// --------------------------------------------------------------------------
// Specialization Packed String
// --------------------------------------------------------------------------

/*!
 * @class PackedString Packed String
 * @extends String
 * @headerfile <seqan/sequence.h>
 * @brief A string that stores as many values in one machine word as possible.
 * 
 * @signature template <typename TValue, typename THostSpec>
 *            class String<TValue, Packed<THostSpec> >;
 * 
 * @tparam TValue The value type, that is the type of the items/characters
 *                stored in the string.Use @link Value @endlink to get the value
 *                type for a given class.
 * @tparam THostSpec The specializing type.This is the specialization of the
 *                   host string that is used for storing the packed values.
 *                   Default: @link AllocString @endlink
 */

template <typename THostspec = Alloc<> >
struct Packed;

// --------------------------------------------------------------------------
// Metafunction PackedTraits_
// --------------------------------------------------------------------------

template <typename TPackedString>
struct PackedTraits_
{
    typedef typename Size<TPackedString>::Type                          TSize;
    typedef typename Value<TPackedString>::Type                         TValue;
    typedef typename Value<typename Host<TPackedString>::Type>::Type    THostValue;

    enum
    {
        BITS_PER_VALUE = BitsPerValue<TValue>::VALUE,
        VALUES_PER_HOST_VALUE = LENGTH<THostValue>::VALUE,
        MAX_BIT_POS = (VALUES_PER_HOST_VALUE - 1) * BITS_PER_VALUE
    };

    static inline
    typename Size<typename Host<TPackedString>::Type>::Type
    toHostLength(typename Size<TPackedString>::Type len)
    {
        return (len + VALUES_PER_HOST_VALUE - 1) / VALUES_PER_HOST_VALUE;
    }
};

/**
.Spec.Packed String:
..cat:Strings
..general:Class.String
..summary:A string that stores as many values in one machine word as possible.
..signature:String<TValue, Packed<THostspec> >
..param.TValue:The value type, that is the type of the items/characters stored in the string.
...remarks:Use @Metafunction.Value@ to get the value type for a given class.
..param.THostspec:The specializing type.
...remarks:This is the specialization of the host string that is used for storing the packed values.
...default:@Spec.Alloc String.Alloc<>@
..include:seqan/sequence.h
*/

/*???TODO Optimierungsm�glichkeiten:
- _clearSpace kopiert Zeichenweise im Packed-String, und nicht im Host-String
- _clearSpace verwendet resize, um den Host zu vergr��ern, d.h. der Inhalt wird eventuell doppelt kopiert.
*/

template <typename TValue, typename THostspec>
class String<TValue, Packed<THostspec> >
{
public:
    typedef typename Host<String>::Type THost;
    typedef typename Size<String>::Type TSize;
    typedef PackedTraits_<String>       TTraits;

    THost data_host;
    TSize data_length;

    String():
        data_length(0)
    {
    }

    template <typename TSource>
    String(TSource & source):
        data_length(0)
    {
        assign(*this, source);
    }
    template <typename TSource>
    String(TSource const & source):
        data_length(0)
    {
        assign(*this, source);
    }
    String(String const & source):
        data_length(0)
    {
        assign(*this, source);
    }

    template <typename TSource>
    String & operator =(TSource const & source)
    {
        assign(*this, source);
        return *this;
    }
    String & operator =(String const & source)
    {
        assign(*this, source);
        return *this;
    }


    // ----------------------------------------------------------------------
    // Subscription operators; have to be defined in class def.
    // ----------------------------------------------------------------------

    template <typename TPos>
    inline typename Reference<String>::Type
    operator[](TPos pos)
    {
        return value(*this, pos);
    }

    template <typename TPos>
    inline typename Reference<String const>::Type 
    operator[](TPos pos) const
    {
        return data_host[pos / TTraits::VALUES_PER_HOST_VALUE][pos % TTraits::VALUES_PER_HOST_VALUE];
    }
};

// --------------------------------------------------------------------------
// Specialization Packed String Iter
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
class Iter<TPackedString, Packed<THostspec> >
{
public:
    typedef typename HostIterator<Iter>::Type THostIterator;
    typedef typename Position<TPackedString>::Type TPosition;

    THostIterator data_iterator;
    unsigned char data_localpos;

    Iter(THostIterator host_begin):
          data_iterator(host_begin),
          data_localpos(0)
    {
    }

    Iter(THostIterator host_begin, TPosition pos_):
          data_iterator(host_begin),
          data_localpos(0)
    {
        *this += pos_;
    }

    Iter(TPackedString &container):
          data_iterator(begin(host(container), Standard())),
          data_localpos(0)
    {
    }

    Iter(TPackedString &container, TPosition pos_):
          data_iterator(begin(host(container), Standard())),
          data_localpos(0)
    {
        *this += pos_;
    }

//    inline
//    Iter const & 
//    operator=(Iter const & other_)
//    {
//        data_iterator = other_.data_iterator;
//        data_localpos = other_.data_localpos;
//        return *this;
//    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction DefaultOverflowImplicit
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
struct DefaultOverflowImplicit<String<TValue, Packed<THostspec> > >
        : DefaultOverflowImplicit<typename Host<String<TValue, Packed<THostspec> > >::Type>
{};

template <typename TValue, typename THostspec>
struct DefaultOverflowImplicit<String<TValue, Packed<THostspec> > const>
        : DefaultOverflowImplicit<typename Host<String<TValue, Packed<THostspec> > const>::Type>
{};

// --------------------------------------------------------------------------
// Metafunction DefaultOverflowExplicit
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
struct DefaultOverflowExplicit<String<TValue, Packed<THostspec> > >
        : DefaultOverflowExplicit<typename Host<String<TValue, Packed<THostspec> > >::Type>
{};

template <typename TValue, typename THostspec>
struct DefaultOverflowExplicit<String<TValue, Packed<THostspec> > const>
        : DefaultOverflowExplicit<typename Host<String<TValue, Packed<THostspec> > const>::Type>
{};

// --------------------------------------------------------------------------
// Metafunction IsContiguous
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
struct IsContiguous<String<TValue, Packed<THostspec> > >
{
    typedef False Type;
    enum { VALUE = false };
};

// --------------------------------------------------------------------------
// Metafunction PackedHostValue_
// --------------------------------------------------------------------------

template <typename TString>
struct PackedHostValue_
{
    typedef typename Value<TString>::Type TValue;
    typedef Tuple<TValue, 64 / BitsPerValue<TValue>::VALUE, BitPacked<> > Type; // use 64bit words
};

// --------------------------------------------------------------------------
// Metafunction Host
// --------------------------------------------------------------------------

///.Metafunction.Host.param.T.type:Spec.Packed String
template <typename TValue, typename THostspec>
struct Host<String<TValue, Packed<THostspec> > >
{
    typedef typename PackedHostValue_<String<TValue, Packed<THostspec> > >::Type TInternalValue;
    typedef String<TInternalValue, THostspec> Type;
};

template <typename TValue, typename THostspec>
struct Host<String<TValue, Packed<THostspec> > const>
{
    typedef typename PackedHostValue_<String<TValue, Packed<THostspec> > >::Type TInternalValue;
    typedef String<TInternalValue, THostspec> const Type;
};

// --------------------------------------------------------------------------
// Metafunction GetValue
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
struct GetValue<String<TValue, Packed<THostspec> > > :
    public Value<String<TValue, Packed<THostspec> > > {};

template <typename TValue, typename THostspec>
struct GetValue<String<TValue, Packed<THostspec> > const> :
    public Value<String<TValue, Packed<THostspec> > const> {};

// --------------------------------------------------------------------------
// Metafunction Reference
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
struct Reference<String<TValue, Packed<THostspec> > >
{
    typedef typename Iterator<String<TValue, Packed<THostspec> >, Standard>::Type TIterator;
    typedef Proxy<IteratorProxy<TIterator> > Type;
};

template <typename TValue, typename THostspec>
struct Reference<String<TValue, Packed<THostspec> > const> :
    public GetValue<String<TValue, Packed<THostspec> > const> {};

// --------------------------------------------------------------------------
// Metafunction Size
// --------------------------------------------------------------------------

/*
template <typename TValue, typename THostspec>
struct Size<String<TValue, Packed<THostspec> > >
{
    typedef __int64 Type;
};
template <typename TValue, typename THostspec>
struct Size<String<TValue, Packed<THostspec> > const>
{
    typedef __int64 Type;
};
*/

// --------------------------------------------------------------------------
// Metafunction Iterator
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec, typename TSpec>
struct Iterator<String<TValue, Packed<THostspec> >, TSpec>
{
    typedef Iter<String<TValue, Packed<THostspec> >, Packed<THostspec> > Type;
};

template <typename TValue, typename THostspec, typename TSpec>
struct Iterator<String<TValue, Packed<THostspec> > const, TSpec>
{
    typedef Iter<String<TValue, Packed<THostspec> > const, Packed<THostspec> > Type;
};

// --------------------------------------------------------------------------
// Metafunction HostIterator
// --------------------------------------------------------------------------

// TODO(holtgrew): Actually is internal, mark so here and rename to HostIterator_
template <typename T>
struct HostIterator;

template <typename TPackedString, typename THostspec>
struct HostIterator<Iter<TPackedString, Packed<THostspec> > > :
    public Iterator<typename Host<TPackedString>::Type, Standard> {};

template <typename TPackedString, typename THostspec>
struct HostIterator<Iter<TPackedString, Packed<THostspec> > const>
{
    typedef typename Host<TPackedString>::Type THost_;
    typedef typename Iterator<THost_, Standard>::Type const Type;
};

// --------------------------------------------------------------------------
// Internal Metafunction TempCopy_
// --------------------------------------------------------------------------

// Note: this works only, if the copy assignment is done without using TempCopy_.
template <typename TValue, typename THostspec>
struct TempCopy_<String<TValue, Packed<THostspec> > >
{
    typedef String<TValue, Packed<THostspec> > Type;
};

// ============================================================================
// Functions
// ============================================================================

// ****************************************************************************
// Functions for Packed String
// ****************************************************************************

// --------------------------------------------------------------------------
// Function host
// --------------------------------------------------------------------------

///.Function.host.param.object.type:Spec.Packed String
///.Function.host.class:Spec.Packed String

template <typename TValue, typename THostspec>
inline typename Host<String<TValue, Packed<THostspec> > >::Type &
host(String<TValue, Packed<THostspec> > & me)
{
    return me.data_host;
}

template <typename TValue, typename THostspec>
inline typename Host<String<TValue, Packed<THostspec> > const>::Type const &
host(String<TValue, Packed<THostspec> > const & me)
{
    return me.data_host;
}

// --------------------------------------------------------------------------
// Function length
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
inline typename Size<String<TValue, Packed<THostspec> > >::Type
length(String<TValue, Packed<THostspec> > & me) 
{
    return me.data_length;
}

template <typename TValue, typename THostspec>
inline typename Size<String<TValue, Packed<THostspec> > const>::Type
length(String<TValue, Packed<THostspec> > const & me) 
{
    return me.data_length;
}

// --------------------------------------------------------------------------
// Function _setLength
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec, typename TSize>
inline void 
_setLength(
    String<TValue, Packed<THostspec> > & me, 
    TSize new_length)
{
    typedef String<TValue, Packed<THostspec> > TString;
    me.data_length = new_length;
    _setLength(host(me), PackedTraits_<TString>::toHostLength(new_length));
}

// --------------------------------------------------------------------------
// Function assign()
//
// Helpers: _assignCopyPackedString()
// --------------------------------------------------------------------------

// optimized variant for copy assignment. The host sequence is copied instead of
// copying the packed string value by value.
template <typename TTarget, typename TSource, typename TTag>
inline void 
_assignCopyPackedString(TTarget & target,
                           TSource & source,
                           Tag<TTag> const & tag)
{
    typedef typename Size<TTarget>::Type TSize2;

    assign(host(target), host(source), tag);
    TSize2 new_length_limit = length(host(target)) * PackedTraits_<TTarget>::VALUES_PER_HOST_VALUE;
    TSize2 new_length = length(source);
    if (new_length > new_length_limit)
    {
        new_length = new_length_limit;
    }
    _setLength(target, new_length);
}

template <typename TTarget, typename TSource, typename TSize, typename TTag>
inline void 
_assignCopyPackedString(TTarget & target,
                        TSource & source,
                        TSize limit,
                        Tag<TTag> const & tag)
{
    typedef typename Size<TTarget>::Type TSize2;

    TSize2 host_limit = PackedTraits_<TTarget>::toHostLength(limit);
    assign(host(target), host(source), host_limit, tag);
    TSize2 new_length_limit = length(host(target)) * PackedTraits_<TTarget>::VALUES_PER_HOST_VALUE;
    TSize2 new_length = length(source);
    if (new_length > new_length_limit)
    {
        new_length = new_length_limit;
    }
    if (new_length > limit)
    {
        new_length = limit;
    }
    _setLength(target, new_length);
}

template <typename TValue, typename THostspec, typename TTag>
inline void 
assign(String<TValue, Packed<THostspec> > & target,
       String<TValue, Packed<THostspec> > & source,
       Tag<TTag> const & tag)
{
    _assignCopyPackedString(target, source, tag);
}

template <typename TValue, typename THostspec, typename TTag>
inline void 
assign(String<TValue, Packed<THostspec> > & target,
       String<TValue, Packed<THostspec> > const & source,
       Tag<TTag> const & tag)
{
    _assignCopyPackedString(target, source, tag);
}

template <typename TValue, typename THostspec, typename TSize, typename TTag>
void assign(String<TValue, Packed<THostspec> > & target,
            String<TValue, Packed<THostspec> > & source,
            TSize limit,
            Tag<TTag> const & tag)
{
    _assignCopyPackedString(target, source, limit, tag);
}
template <typename TValue, typename THostspec, typename TSize, typename TTag>
void assign(String<TValue, Packed<THostspec> > & target,
            String<TValue, Packed<THostspec> > const & source,
            TSize limit,
            Tag<TTag> const & tag)
{
    _assignCopyPackedString(target, source, limit, tag);
}

// --------------------------------------------------------------------------
// Function getObjectId()
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
inline void const * 
getObjectId(String<TValue, Packed<THostspec> > const & me)
{
    return getObjectId(host(me));
}

// --------------------------------------------------------------------------
// Function iter()
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec, typename TPos, typename TTag>
inline typename Iterator<String<TValue, Packed<THostspec> >, Tag<TTag> const>::Type 
iter(String<TValue, Packed<THostspec> > & me,
     TPos pos_,
     Tag<TTag> const &)
{
    typedef typename Iterator<String<TValue, Packed<THostspec> >, Tag<TTag> const>::Type TIterator;
    return TIterator(me, pos_);
}

template <typename TValue, typename THostspec, typename TPos, typename TTag>
inline typename Iterator<String<TValue, Packed<THostspec> > const, Tag<TTag> const>::Type 
iter(String<TValue, Packed<THostspec> > const & me,
     TPos pos_,
     Tag<TTag> const &)
{
    typedef typename Iterator<String<TValue, Packed<THostspec> > const, Tag<TTag> const>::Type TIterator;
    return TIterator(me, pos_);
}

// --------------------------------------------------------------------------
// Function begin()
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec, typename TTag>
inline typename Iterator<String<TValue, Packed<THostspec> >, Tag<TTag> const>::Type 
begin(String<TValue, Packed<THostspec> > & me,
      Tag<TTag> const &)
{
    typedef typename Iterator<String<TValue, Packed<THostspec> >, Tag<TTag> const>::Type TIterator;
    return TIterator(me);
}

template <typename TValue, typename THostspec, typename TTag>
inline typename Iterator<String<TValue, Packed<THostspec> > const, Tag<TTag> const>::Type 
begin(String<TValue, Packed<THostspec> > const & me,
      Tag<TTag> const &)
{
    typedef typename Iterator<String<TValue, Packed<THostspec> > const, Tag<TTag> const>::Type TIterator;
    return TIterator(me);
}

// --------------------------------------------------------------------------
// Function end()
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec, typename TTag>
inline typename Iterator<String<TValue, Packed<THostspec> >, Tag<TTag> const>::Type
end(String<TValue, Packed<THostspec> > & me,
    Tag<TTag> const & tag_)
{
    return iter(me, length(me), tag_);
}

template <typename TValue, typename THostspec, typename TTag>
inline typename Iterator<String<TValue, Packed<THostspec> > const, Tag<TTag> const>::Type 
end(String<TValue, Packed<THostspec> > const & me,
    Tag<TTag> const & tag_)
{
    return iter(me, length(me), tag_);
}

// --------------------------------------------------------------------------
// Function value()
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec, typename TPos>
inline typename Reference<String<TValue, Packed<THostspec> > >::Type
value(String<TValue, Packed<THostspec> > & me, 
      TPos pos)
{
    return *iter(me, pos, Standard());
} 

template <typename TValue, typename THostspec, typename TPos>
inline typename Reference<String<TValue, Packed<THostspec> > const>::Type
value(String<TValue, Packed<THostspec> > const & me, 
      TPos pos)
{
    typedef String<TValue, Packed<THostspec> > TPackedString;
    typedef PackedTraits_<TPackedString> TTraits;
    return me.data_host[pos / TTraits::VALUES_PER_HOST_VALUE][pos % TTraits::VALUES_PER_HOST_VALUE];
} 

// --------------------------------------------------------------------------
// Function capacity()
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
inline typename Size<String<TValue, Packed<THostspec> > const>::Type
capacity(String<TValue, Packed<THostspec> > const & me)
{
    typedef String<TValue, Packed<THostspec> > TPackedString;
    typedef PackedTraits_<TPackedString> TTraits;
    typedef typename Size<TPackedString>::Type TSize;
    
    return capacity(host(me)) * (TSize)TTraits::VALUES_PER_HOST_VALUE;
}

// --------------------------------------------------------------------------
// Function clear()
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
inline void 
clear(String<TValue, Packed<THostspec> > & me)
{
    clear(host(me));
    _setLength(me, 0);
}

// --------------------------------------------------------------------------
// Function _clearSpace()
//
// Helper struct ClearSpaceStringPacked_.
// --------------------------------------------------------------------------

//implementation for all expand tags other than "limit"
template <typename TExpand>
struct ClearSpaceStringPacked_
{
    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq, 
        typename Size<T>::Type size)
    {
        typedef typename Size<T>::Type TSize;
        TSize wanted_host_length = PackedTraits_<T>::toHostLength(size);
        TSize new_host_length = resize(host(seq), wanted_host_length, TExpand());
        if (new_host_length < wanted_host_length)
        {
            size = new_host_length * PackedTraits_<T>::VALUES_PER_HOST_VALUE;
        }
        _setLength(seq, size);
        return size;
    }

    template <typename T>
    static inline typename Size<T>::Type 
    _clearSpace_(
        T & seq, 
        typename Size<T>::Type size,
        typename Size<T>::Type limit)
    {
        if (limit < size)
        {
            size = limit;
        }
        return _clearSpace_(seq, limit);
    }

    template <typename T>
    static inline typename Size<T>::Type 
    _clearSpace_(
        T & seq, 
        typename Size<T>::Type size, 
        typename Size<T>::Type start, 
        typename Size<T>::Type end)
    {
        return _clearSpace_(seq, size, start, end, maxValue<typename Size<T>::Type >());
    }

    template <typename T>
    static typename Size<T>::Type 
    _clearSpace_(
        T & seq, 
        typename Size<T>::Type size, 
        typename Size<T>::Type start, 
        typename Size<T>::Type end, 
        typename Size<T>::Type limit)
    {
//??? TODO: This function can be accelerated this way: 
//              - move values in host
//              - avoid double moving of the rest-part if "resize" allocates a new block

        typedef typename Size<T>::Type TSize;

        TSize old_length = length(seq);
        TSize old_size = end - start;
        TSize wanted_new_length = old_length + size - old_size;

        if (wanted_new_length > limit)
        {
            wanted_new_length = limit;
        }

        TSize wanted_host_length = PackedTraits_<T>::toHostLength(wanted_new_length);
        TSize new_host_length = resize(host(seq), wanted_host_length, TExpand());

        TSize new_length;
        if (new_host_length < wanted_host_length)
        {
            new_length = new_host_length * PackedTraits_<T>::VALUES_PER_HOST_VALUE;
            if (new_length <= start + size)
            {
                goto FINISH;
            }
            old_length = new_length - size + old_size;
        }
        else
        {
            new_length = wanted_new_length;
        }
/*
        //move [end:right_end] to [start + size:..]
        if (old_size > size)
        {//move rest to left
            ::std::copy(iter(seq, end, Standard()), iter(seq, old_length, Standard()), iter(seq, end + size - old_size, Standard()));
        }
        else
        {//move rest to right
            ::std::copy_backward(iter(seq, end, Standard()), iter(seq, old_length, Standard()), iter(seq,  new_length, Standard()));
        }
*/
        if (old_size > size)
        {
            arrayMoveForward(iter(seq, end, Standard()), iter(seq, old_length, Standard()), iter(seq, end + size - old_size, Standard()));
        }
        else
        {
            arrayMoveBackward(iter(seq, end, Standard()), iter(seq, old_length, Standard()), iter(seq, end + size - old_size, Standard()));
        }
FINISH:
        _setLength(seq, new_length);
        return size;
    }
/*
    template <typename T>
    static inline typename Size<T>::Type 
    _clearSpace_(
        T & seq, 
        typename Size<T>::Type size, 
        typename Iterator<T>::Type start, 
        typename Iterator<T>::Type end)
    {
        typename Iterator<T>::Type seq_begin = begin(seq);
        return _clearSpace(seq, size, start - seq_begin, end - seq_begin, Insist());
    }

    template <typename T>
    static inline typename Size<T>::Type 
    _clearSpace_(
        T & seq, 
        typename Size<T>::Type size,  
        typename Iterator<T>::Type start,
        typename Iterator<T>::Type end,
        typename Size<T>::Type limit) 
    {
        typename Iterator<T>::Type seq_begin = begin(seq);
        return _clearSpace(seq, size, start - seq_begin, end - seq_begin, limit, Insist());
    }
*/
};

template<typename TValue, typename THostspec, typename TExpand>
inline typename Size< String<TValue, Packed<THostspec> > >::Type 
_clearSpace(String<TValue, Packed<THostspec> > & me, 
        typename Size< String<TValue, Packed<THostspec> > >::Type size, 
        Tag<TExpand>)
{
    return ClearSpaceStringPacked_<Tag<TExpand> >::_clearSpace_(me, size);
}

template<typename TValue, typename THostspec, typename TExpand>
inline typename Size< String<TValue, Packed<THostspec> > >::Type 
_clearSpace(String<TValue, Packed<THostspec> > & me, 
        typename Size< String<TValue, Packed<THostspec> > >::Type size, 
        typename Size< String<TValue, Packed<THostspec> > >::Type limit, 
        Tag<TExpand>)
{
    return ClearSpaceStringPacked_<Tag<TExpand> >::_clearSpace_(me, size, limit);
}

template<typename TValue, typename THostspec, typename TPosition, typename TExpand>
inline typename Size< String<TValue, Packed<THostspec> > >::Type 
_clearSpace(String<TValue, Packed<THostspec> > & me, 
            typename Size< String<TValue, Packed<THostspec> > >::Type size, 
            TPosition pos_begin, 
            TPosition pos_end, 
            Tag<TExpand>)
{
    return ClearSpaceStringPacked_<Tag<TExpand> >::_clearSpace_(me, size, pos_begin, pos_end);
}

template<typename TValue, typename THostspec, typename TPosition, typename TExpand>
inline typename Size< String<TValue, Packed<THostspec> > >::Type 
_clearSpace(String<TValue, Packed<THostspec> > & me, 
            typename Size< String<TValue, Packed<THostspec> > >::Type size, 
            TPosition pos_begin, 
            TPosition pos_end, 
            typename Size< String<TValue, Packed<THostspec> > >::Type limit, 
            Tag<TExpand>)
{
    return ClearSpaceStringPacked_<Tag<TExpand> >::_clearSpace_(me, size, pos_begin, pos_end, limit);
}

// --------------------------------------------------------------------------
// Function reserve()
// --------------------------------------------------------------------------

///.Function.reserve.param.object.type:Spec.Packed String

template <typename TValue, typename TSpec, typename TSize_, typename TExpand>
inline typename Size< String<TValue, Packed<TSpec> > >::Type
reserve(
    String<TValue, Packed<TSpec> > & seq, 
    TSize_ new_capacity,
    Tag<TExpand> tag)
{

    typedef String<TValue, Packed<TSpec> > TString;
    typedef typename Size<TString>::Type TSize;
    TSize ret_value = reserve(host(seq), PackedTraits_<TString>::toHostLength(new_capacity), tag);
    return ret_value * PackedTraits_<TString>::VALUES_PER_HOST_VALUE;
}

// ****************************************************************************
// Functions for Packed String Iter
// ****************************************************************************

// --------------------------------------------------------------------------
// Function hostIterator()
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
inline typename HostIterator<Iter<TPackedString, Packed<THostspec> > >::Type &
hostIterator(Iter<TPackedString, Packed<THostspec> > & me)
{
    return me.data_iterator;
}

template <typename TPackedString, typename THostspec>
inline typename HostIterator<Iter<TPackedString, Packed<THostspec> > const>::Type  &
hostIterator(Iter<TPackedString, Packed<THostspec> > const & me)
{
    return me.data_iterator;
}


// --------------------------------------------------------------------------
// Function value()
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
inline typename Reference<Iter<TPackedString, Packed<THostspec> > >::Type 
value(Iter<TPackedString, Packed<THostspec> > & me)
{
    return typename Reference<Iter<TPackedString, Packed<THostspec> > >::Type(me);
}

template <typename TPackedString, typename THostspec>
inline typename Reference<Iter<TPackedString, Packed<THostspec> > const>::Type 
value(Iter<TPackedString, Packed<THostspec> > const & me)
{
    return typename Reference<Iter<TPackedString, Packed<THostspec> > const>::Type(me);
}

template <typename TPackedString, typename THostspec>
inline typename Reference<Iter<TPackedString const, Packed<THostspec> > >::Type
value(Iter<TPackedString const, Packed<THostspec> > & me)
{
    return getValue(hostIterator(me))[me.data_localpos];
}

template <typename TPackedString, typename THostspec>
inline typename Reference<Iter<TPackedString const, Packed<THostspec> > const>::Type
value(Iter<TPackedString const, Packed<THostspec> > const & me)
{
    return getValue(hostIterator(me))[me.data_localpos];
}

// --------------------------------------------------------------------------
// Function getValue()
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
inline typename GetValue<Iter<TPackedString, Packed<THostspec> > >::Type 
getValue(Iter<TPackedString, Packed<THostspec> > & me)
{
    return getValue(hostIterator(me))[me.data_localpos];
}

template <typename TPackedString, typename THostspec>
inline typename GetValue<Iter<TPackedString, Packed<THostspec> > const>::Type 
getValue(Iter<TPackedString, Packed<THostspec> > const & me)
{
    return getValue(hostIterator(me))[me.data_localpos];
}

// --------------------------------------------------------------------------
// Function assignValue()
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec, typename TValue>
inline void
assignValue(Iter<TPackedString, Packed<THostspec> > const & me,
            TValue const & _value)
{
    typedef Iter<TPackedString, Packed<THostspec> > const TIterator;
    assignValue(value(hostIterator(me)), me.data_localpos, (typename Value<TIterator>::Type)_value);
}

template <typename TPackedString, typename THostspec, typename TValue>
inline void
assignValue(Iter<TPackedString, Packed<THostspec> > & me,
            TValue const & _value)
{
    typedef Iter<TPackedString, Packed<THostspec> > TIterator;
    assignValue(value(hostIterator(me)), me.data_localpos, (typename Value<TIterator>::Type)_value);
}

// --------------------------------------------------------------------------
// Function moveValue()
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec, typename TValue>
inline void
moveValue(Iter<TPackedString, Packed<THostspec> > & me,
          TValue const & _value)
{
    assignValue(me, _value);
}

template <typename TPackedString, typename THostspec, typename TValue>
inline void
moveValue(Iter<TPackedString, Packed<THostspec> > const & me,
          TValue const & _value)
{
    assignValue(me, _value);
}

// --------------------------------------------------------------------------
// Function valueConstruct()
// --------------------------------------------------------------------------

//emulate construction and destruction 

template <typename TPackedString, typename THostspec>
inline void
valueConstruct(Iter<TPackedString, Packed<THostspec> > const & /*it*/)
{
    // TODO(holtgrew): Why not assign default-constructed?
}

template <typename TPackedString, typename THostspec, typename TParam>
inline void
valueConstruct(Iter<TPackedString, Packed<THostspec> > const & it,
               TParam const & param_)
{
    assignValue(it, param_);
}

template <typename TPackedString, typename THostspec, typename TParam>
inline void
valueConstruct(Iter<TPackedString, Packed<THostspec> > const & it,
               TParam const & param_,
               Move const & /*tag*/)
{
    moveValue(it, param_);
}

// --------------------------------------------------------------------------
// Function valueDestruct()
// --------------------------------------------------------------------------

// Packed strings cannot contain non-POD data types.

template <typename TPackedString, typename THostspec>
inline void
valueDestruct(Iter<TPackedString, Packed<THostspec> > const & /*it*/)
{
}

// --------------------------------------------------------------------------
// Function operator==()
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
inline bool 
operator==(Iter<TPackedString, Packed<THostspec> > const & left,
           Iter<TPackedString, Packed<THostspec> > const & right)
{
    return (hostIterator(left) == hostIterator(right)) && (left.data_localpos == right.data_localpos);
}

// --------------------------------------------------------------------------
// Function operator!=()
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
inline bool 
operator!=(Iter<TPackedString, Packed<THostspec> > const & left,
           Iter<TPackedString, Packed<THostspec> > const & right)
{
    return (hostIterator(left) != hostIterator(right)) || (left.data_localpos != right.data_localpos);
}

// --------------------------------------------------------------------------
// Function operator>()
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
inline bool 
operator>(Iter<TPackedString, Packed<THostspec> > const & left,
          Iter<TPackedString, Packed<THostspec> > const & right)
{
    return (hostIterator(left) > hostIterator(right)) || ((hostIterator(left) == hostIterator(right)) && (left.data_localpos > right.data_localpos));
}

// --------------------------------------------------------------------------
// Function operator>=()
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
inline bool 
operator>=(Iter<TPackedString, Packed<THostspec> > const & left,
           Iter<TPackedString, Packed<THostspec> > const & right)
{
    return (hostIterator(left) > hostIterator(right)) || ((hostIterator(left) == hostIterator(right)) && (left.data_localpos >= right.data_localpos));
}

// --------------------------------------------------------------------------
// Function operator<()
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
inline bool 
operator<(Iter<TPackedString, Packed<THostspec> > const & left,
          Iter<TPackedString, Packed<THostspec> > const & right)
{
    return (hostIterator(left) < hostIterator(right)) || ((hostIterator(left) == hostIterator(right)) && (left.data_localpos < right.data_localpos));
}

// --------------------------------------------------------------------------
// Function operator<=()
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
inline bool 
operator <= (Iter<TPackedString, Packed<THostspec> > const & left,
             Iter<TPackedString, Packed<THostspec> > const & right)
{
    return (hostIterator(left) < hostIterator(right)) || ((hostIterator(left) == hostIterator(right)) && (left.data_localpos <= right.data_localpos));
}

// --------------------------------------------------------------------------
// Function operator++()
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
inline Iter<TPackedString, Packed<THostspec> > &
operator++(Iter<TPackedString, Packed<THostspec> > & me)
{
    if (++me.data_localpos == PackedTraits_<TPackedString>::VALUES_PER_HOST_VALUE)
    {
        me.data_localpos = 0;
        ++hostIterator(me);
    }
    return me;
}

// --------------------------------------------------------------------------
// Function operator--()
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
inline Iter<TPackedString, Packed<THostspec> > &
operator--(Iter<TPackedString, Packed<THostspec> > & me)
{
    if (me.data_localpos-- == 0)
    {
        me.data_localpos = PackedTraits_<TPackedString>::VALUES_PER_HOST_VALUE - 1;
        --hostIterator(me);
    }
    return me;
}

// ----------------------------------------------------------------------------
// Helper Function _helperIsNegative()
// ----------------------------------------------------------------------------

// to remove '... < 0 is always false' warning
template <typename T>
inline bool
_isNegative(T, False)
{
    return false;
}

template <typename T>
inline bool
_isNegative(T t, True)
{
    return t < 0;
}

template <typename T>
inline bool
_isNegative(T t)
{
    return _isNegative(t, typename IsSameType<T, typename MakeSigned_<T>::Type>::Type());
}

// --------------------------------------------------------------------------
// Function operator+() for (iter, integral)
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec, typename TIntegral>
inline Iter<TPackedString, Packed<THostspec> >  
operator+(Iter<TPackedString, Packed<THostspec> > const & iter,
          TIntegral const & delta)
{
    typedef PackedTraits_<TPackedString> TTraits;

    if (_isNegative(delta))
        return iter - (typename MakeUnsigned<TIntegral>::Type)(-delta);

    TIntegral ofs = (TIntegral)iter.data_localpos + delta;
    return Iter<TPackedString, Packed<THostspec> >(
        hostIterator(iter) + ofs / TTraits::VALUES_PER_HOST_VALUE,
        ofs % TTraits::VALUES_PER_HOST_VALUE);
}

template <typename TPackedString, typename THostspec, typename TIntegral>
inline Iter<TPackedString, Packed<THostspec> >  
operator+(TIntegral const & delta,
          Iter<TPackedString, Packed<THostspec> > const & iter)
{
    typedef PackedTraits_<TPackedString> TTraits;

    if (_isNegative(delta))
        return iter - (typename MakeUnsigned<TIntegral>::Type)(-delta);

    TIntegral ofs = (TIntegral)iter.data_localpos + delta;
    return Iter<TPackedString, Packed<THostspec> >(
        hostIterator(iter) + ofs / TTraits::VALUES_PER_HOST_VALUE,
        ofs % TTraits::VALUES_PER_HOST_VALUE);
}

// --------------------------------------------------------------------------
// Function operator+=() for (iter, integral)
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec, typename TIntegral>
inline Iter<TPackedString, Packed<THostspec> > &
operator+=(Iter<TPackedString, Packed<THostspec> > & iter,
           TIntegral const & delta)
{
    typedef PackedTraits_<TPackedString> TTraits;

    if (_isNegative(delta))
        return iter -= (typename MakeUnsigned<TIntegral>::Type)(-delta);

    TIntegral ofs = (TIntegral)iter.data_localpos + delta;
    hostIterator(iter) += ofs / TTraits::VALUES_PER_HOST_VALUE;
    iter.data_localpos = ofs % TTraits::VALUES_PER_HOST_VALUE;
    return iter;
}

// --------------------------------------------------------------------------
// Function operator-() for (iter, integral)
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec, typename TIntegral>
inline Iter<TPackedString, Packed<THostspec> >  
operator-(Iter<TPackedString, Packed<THostspec> > const & iter,
          TIntegral const & delta)
{
    typedef PackedTraits_<TPackedString> TTraits;

    if (_isNegative(delta))
        return iter + (typename MakeUnsigned<TIntegral>::Type)(-delta);

    TIntegral ofs = delta + (TIntegral)(TTraits::VALUES_PER_HOST_VALUE - 1) - (TIntegral)iter.data_localpos;
    return Iter<TPackedString, Packed<THostspec> >(
        hostIterator(iter) - ofs / TTraits::VALUES_PER_HOST_VALUE,
        (TTraits::VALUES_PER_HOST_VALUE - 1) - (ofs % TTraits::VALUES_PER_HOST_VALUE));
}

// --------------------------------------------------------------------------
// Function operator-=() for (iter, integral)
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec, typename TIntegral>
inline Iter<TPackedString, Packed<THostspec> > &
operator-=(Iter<TPackedString, Packed<THostspec> > & iter,
           TIntegral const & delta)
{
    typedef PackedTraits_<TPackedString> TTraits;

    if (_isNegative(delta))
        return iter += (typename MakeUnsigned<TIntegral>::Type)(-delta);

    TIntegral ofs = delta + (TIntegral)(TTraits::VALUES_PER_HOST_VALUE - 1) - (TIntegral)iter.data_localpos;
    hostIterator(iter) -= ofs / TTraits::VALUES_PER_HOST_VALUE;
    iter.data_localpos = (TTraits::VALUES_PER_HOST_VALUE - 1) - (ofs % TTraits::VALUES_PER_HOST_VALUE);
    return iter;
}

// --------------------------------------------------------------------------
// Function operator-() for (iter, iter)
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
inline typename Difference<Iter<TPackedString, Packed<THostspec> > >::Type  
operator-(Iter<TPackedString, Packed<THostspec> > const & iterLeft,
          Iter<TPackedString, Packed<THostspec> > const & iterRight)
{
    typedef PackedTraits_<TPackedString> TTraits;
    typedef typename Difference<Iter<TPackedString, Packed<THostspec> > >::Type TDiff;
    return (TDiff)(hostIterator(iterLeft) - hostIterator(iterRight)) * (TDiff)TTraits::VALUES_PER_HOST_VALUE +
           (TDiff)iterLeft.data_localpos - (TDiff)iterRight.data_localpos;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQUENCE_STRING_PACKED_H_
