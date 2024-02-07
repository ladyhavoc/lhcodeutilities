/*
Copyright (c) 2024 Ashley Rose Hale (LadyHavoc)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

/// This is a single-header Predictive Partial Matching + Range Coder library
/// intended for use on network packets, it can also tell you if a packet is
/// corrupted in transit, taking the place of a regular CRC.

#pragma once

#ifndef LHNETCODEC_H
#define LHNETCODEC_H

#include <string.h>
#include <stdint.h>

typedef enum LHNETCODEC_enum
{
	LHNETCODEC_STATUS_UNINITIALIZED = 0,
	LHNETCODEC_STATUS_OK = 0x1000,
	LHNETCODEC_STATUS_INVALID_PARAMETER = 0x2000,
	LHNETCODEC_STATUS_ENCODE_FULL = 0x3000,
	LHNETCODEC_STATUS_DECODE_EOF = 0x4000,
	LHNETCODEC_STATUS_DECODE_CORRUPT = 0x4001,

	/// number of independent contexts for bits that are tracked by the system.
	/// this limit only affects memory usage, and has no ongoing cpu cost.
	/// 
	/// example table of context ranges:
	/// * encode int64 as firstcontext=0, numcontexts=64
	/// * encode float (unioned to int32) as firstcontext=64, numcontexts=32
	/// * encode UTF-8 text bytes as firstcontext=96, numcontexts=8
	/// * encode ULEB128/SLEB128 bytes as firstcontext=104, numcontexts=8
	///   (ULEB128 = a variable length integer encoding, best known for use in
	///    protobuf)
	/// * encode protocol enums as firstcontext=112, numcontexts=8
	///   (or however many you need)
	/// 
	/// the above example is 120 contexts, which is below this limit
	/// 
	/// notes:
	/// * ULEB128/SLEB128 should work very well with the AdaptiveCoder, it is
	///   quite reasonable to use them instead of encoding int64, while the
	///   AdaptiveCoder can make bits cost far less than a bit, it is still best
	///   to encode the fewest bits you can get away with.
#ifdef LHNETCODEC_SETTINGS_MAXCONTENTRANGE
	LHNETCODEC_MAXCONTEXTRANGE = LHNETCODEC_SETTINGS_MAXCONTENTRANGE,
#else
	LHNETCODEC_MAXCONTEXTRANGE = 256,
#endif
	/// maximum number of bits of context supported - a bit is predicted based
	/// on the same bit position in this many previous numbers
	/// 
	/// NOTE: this increases size of the histogram by powers of 2, tune
	/// carefully.  Values above 10 don't seem to provide much benefit, and have
	/// very significant memory footprint.
#ifdef LHNETCODEC_SETTINGS_MAXCONTEXTLENGTH
	LHNETCODEC_MAXCONTEXTLENGTH = LHNETCODEC_SETTINGS_MAXCONTEXTLENGTH,
#else
	LHNETCODEC_MAXCONTEXTLENGTH = 4,
#endif
	/// number of possible contexts for a bit position
	LHNETCODEC_MAXCONTEXTSIZE = (1 << LHNETCODEC_MAXCONTEXTLENGTH),
	/// the largest histogram that we may have to use, depending on settings
	/// this can be a very big number, if you get a stack overflow it is
	/// probably this limit
	LHNETCODEC_MAXHISTOGRAMCONTEXTS = LHNETCODEC_MAXCONTEXTSIZE * LHNETCODEC_MAXCONTEXTRANGE,
	/// the largest histogram that may have to use, depending on settings
	/// note that this is a very big number, if you get a stack overflow it is
	/// probably this limit
	LHNETCODEC_MAXHISTOGRAMSIZE = LHNETCODEC_MAXHISTOGRAMCONTEXTS * 2,

#if LHNETCODEC_ENABLE_ADAPTIVECODER
	/// the context length is configurable, higher values allow more specialized
	/// adaptation but may result in lower compression ratios on short streams
	/// in the adaptive coder.
	LHNETCODEC_ADAPTIVECODER_DEFAULT_CONTEXTLENGTH = 4,
	/// the initial value of the histogram pretends that all contexts had this
	/// many occurrences for each symbol.
	LHNETCODEC_ADAPTIVECODER_DEFAULT_PROBABILITY_BASE = 1,
	/// the probability of a bit value is incremented by this much, in relation
	/// to base this can be used to make the probabilities adapt faster (higher
	/// compression if correct) or slower (lower cost of misprediction).
	LHNETCODEC_ADAPTIVECODER_DEFAULT_PROBABILITY_STEP = 1,
	/// the maximum bias that the adaptive coder will allow - when a probability
	/// exceeds this value, both probabilities are shifted down
	LHNETCODEC_ADAPTIVECODER_DEFAULT_PROBABILITY_LIMIT = 64,
	/// shift down the probabilities by this much when limit is reached
	LHNETCODEC_ADAPTIVECODER_DEFAULT_PROBABILITY_SHIFT = 1,
#endif

#if LHNETCODEC_ENABLE_BLOCKCODER
	// the first version of the BlockCoder codec
	LHNETCODEC_BLOCKCODER_VERSION_0 = 0x4000,
	// the first version of the BlockCoderFixedPoint codec
	LHNETCODEC_BLOCKCODERFIXEDPOINT_VERSION_0 = 0x4001,

	/// the context size is configurable, higher values allow higher compression
	/// ratios in the block coder, but significantly increase the header size in
	/// the block coder, which may be counter productive on small blocks.
	LHNETCODEC_BLOCKCODER_DEFAULT_CONTEXTLENGTH = 1,
	/// number of bits written for each symbol probability, after the extend bit
	LHNETCODEC_BLOCKCODER_TABLEBITS = 4,
	/// the block coder needs to replay the values being encoded each time it
	/// flushes, so this is the maximum queue it allows for one block
	LHNETCODEC_BLOCKCODER_REPLAY_NUMBERS = 65536,

	/// number of bits used to represent fixed point splits in the histogram
	LHNETCODEC_BLOCKCODERFIXEDPOINT_SPLITBITS = 12,
#endif
}
LHNETCODEC_enum;

// macros that are used in several places in inner loops

// reduce probabilities if they are getting so high as to make the
// divide potentially unsafe (split must not be equal to low or high)
#define MAKE_PROBABILITY_SAFE do{\
	while ((probability_0 | probability_1) > 0x10000)\
	{\
		if (probability_0 > 0)\
		{\
			if (probability_0 >= 2)\
				probability_0 >>= 1;\
			else\
				probability_0 = 1;\
		}\
		if (probability_1 > 0)\
		{\
			if (probability_1 >= 2)\
				probability_1 >>= 1;\
			else\
				probability_1 = 1;\
		}\
	}\
} while(0)
// when the low and high have the same significant byte, emit a byte (this is the difference between range coders which emit bytes and arithmetic coders which emit individual bits)
// this also checks if the range has collapsed (low == high), which is probably
// impossible, but we check anyway.
//
// we adjust cursor_trimmed only when we write a non-zero byte
#define ENCODE_CHECKRANGE do {\
	while (((low ^ high) & 0xFF0000) == 0)\
	{\
		if (cursor >= maxbytes)\
		{\
			state->status = LHNETCODEC_STATUS_ENCODE_FULL;\
			break;\
		}\
		unsigned int b = (low >> 16) & 0xFF;\
		if (coded)\
			coded[cursor] = b;\
		cursor++;\
		if (b)\
			cursor_trimmed = cursor;\
		low = (low << 8) & 0xFFFFFF;\
		high = (high << 8) & 0xFFFFFF;\
		if (low == high)\
		{\
			low = 0;\
			high = 0xFFFFFF;\
		}\
	}\
} while (0)
// encode one number with a specified range of possible values
#define ENCODE_NUMBER(code_number, code_range) do {\
	unsigned int code_low = low + (code_number) * (high + 1 - low) / (code_range);\
	unsigned int code_high = low + ((code_number) + 1) * (high + 1 - low) / (code_range) - 1;\
	low = code_low;\
	high = code_high;\
	ENCODE_CHECKRANGE;\
} while(0)
// check if we need to read another byte from the stream
#define DECODE_CHECKRANGE do {\
	while (((low ^ high) & 0xFF0000) == 0)\
	{\
		if (cursor >= maxbytes)\
		{\
			state->status = LHNETCODEC_STATUS_DECODE_CORRUPT;\
			break;\
		}\
		stream = ((stream << 8) | (coded[cursor++] & 0xFF)) & 0xFFFFFF;\
		low = (low << 8) & 0xFFFFFF;\
		high = (high << 8) & 0xFFFFFF;\
		if (low == high)\
		{\
			low = 0;\
			high = 0xFFFFFF;\
		}\
	}\
} while (0)
// decode one number with a specified range of possible values
#define DECODE_NUMBER(assign_number, code_range) do{\
	unsigned int code_number = (stream - low) * (code_range) / (high + 1 - low);\
	unsigned int code_low = low + code_number * (high + 1 - low) / (code_range);\
	unsigned int code_high = low + (code_number + 1) * (high + 1 - low) / (code_range) - 1;\
	low = code_low;\
	high = code_high;\
	DECODE_CHECKRANGE;\
	(assign_number) = code_number;\
} while(0)

typedef struct LHNETCODEC_AdaptiveCoder_State
{
	/// counts of how many times a 0 or 1 has occurred so far after each
	/// possible context, for each possible bit position, this increments when
	/// encoding each bit, so it adapts to the numbers actually encoded in the
	/// stream, but has an initial phase where it is a poor fit (i.e. completely
	/// uncompressed).
	unsigned char histogram[LHNETCODEC_MAXHISTOGRAMSIZE];
	/// store the previous N bits for each of the bits in a 64bit number
	/// being encoded/decoded, reading/writing smaller sizes will only advance
	/// the contexts for the number of bits specified by the caller.
	unsigned char contexts[LHNETCODEC_MAXCONTEXTRANGE];
	/// context length in bits
	unsigned int contextlength;
	/// context_mask = (1u << contextlength) - 1;
	unsigned int contextlength_mask;
	/// minimum value of the histogram probabilities, probabilities can not go
	/// below this value (in actuality the histogram[] bytes are relative to
	/// this value)
	/// setting this above 65536 might cause math errors.
	/// range: 1 <= base <= 65536,
	unsigned int probability_base;
	/// each time a bit is encoded, the probability of that bit in this context
	/// is increased by this much, the ratio of base to step allows faster
	/// or slower adaptation to the data stream.
	/// setting this above 65536 might cause math errors.
	/// range: 0 <= step <= 255,
	unsigned int probability_step;
	/// maximum value of the histogram probabilities - if this is reached,
	/// probability_maximum_shift is applied to both.
	/// setting this above 65536 might cause math errors.
	/// range: 1 <= limit <= 256.
	unsigned int probability_limit;
	/// how much to shift down the probabilities when maximum is reached, this
	/// is performed on the probability values relative to base, they can never
	/// go below 0, relative to base, so base is the lowest effective value used
	/// in computing the range split.  if this is set to more than the size of
	/// histogram values, the shift down is effectively a reset for the context.
	/// range: 1 <= shift <= 32.
	unsigned char probability_shift;
	/// current coded data buffer which stores the resulting stream, can be NULL
	/// if merely being used for simulations.
	unsigned char* coded_start;
	/// size of the coded working buffer in bytes, when coded_cursor
	/// reaches this value, no more bits can be appended and an error status is
	/// returned.
	size_t coded_maxbytes;
	/// current position to write the next byte in the coded data buffer, this
	/// is for any byte value we may want to put there
	size_t coded_cursor_any_byte;
	/// the trimmed length only goes to the last non-zero byte
	size_t coded_cursor_trimmed;
	/// current operating range in the range codec
	unsigned int range_low;
	/// current operating range in the range codec
	unsigned int range_high;
	/// (decoder only) this is the most recent bits read from the coded buffer,
	/// which is basically a really long binary fraction that guides us to the
	/// range choices that the encoder made when writing the stream.
	unsigned int range_stream;
	/// current operational status of the coder
	LHNETCODEC_enum status;
}
LHNETCODEC_AdaptiveCoder_State;

LHNETCODEC_enum LHNETCODEC_AdaptiveCoder_Encode_Begin(LHNETCODEC_AdaptiveCoder_State* state, unsigned char *coded_start, size_t coded_maxbytes, unsigned int contextlength, unsigned int probability_base, unsigned int probability_step, unsigned int probability_limit, unsigned int probability_shift)
{
	if (contextlength > LHNETCODEC_MAXCONTEXTLENGTH)
		contextlength = LHNETCODEC_MAXCONTEXTLENGTH;
	if (probability_base < 1)
		probability_base = 1;
	if (probability_limit < 1)
		probability_limit = 1;
	if (probability_shift < 1)
		probability_shift = 1;
	if (probability_shift > 32)
		probability_shift = 32;
	memset(state, 0, sizeof(*state));
	state->contextlength = contextlength;
	state->contextlength_mask = (1u << state->contextlength) - 1;
	state->probability_base = probability_base;
	state->probability_step = probability_step;
	state->probability_limit = probability_limit;
	state->probability_shift = probability_shift;
	state->coded_start = coded_start;
	state->coded_maxbytes = coded_maxbytes;
	state->coded_cursor_any_byte = 0;
	state->coded_cursor_trimmed = 0;
	state->range_low = 0;
	state->range_high = 0xFFFFFF;
	state->status = LHNETCODEC_STATUS_OK;
	return state->status;
}

LHNETCODEC_enum LHNETCODEC_AdaptiveCoder_Decode_Begin(LHNETCODEC_AdaptiveCoder_State* state, unsigned char* coded_start, size_t coded_maxbytes, unsigned int contextlength, unsigned int probability_base, unsigned int probability_step, unsigned int probability_limit, unsigned int probability_shift)
{
	if (contextlength > LHNETCODEC_MAXCONTEXTLENGTH)
		contextlength = LHNETCODEC_MAXCONTEXTLENGTH;
	if (probability_base < 1)
		probability_base = 1;
	if (probability_limit < 1)
		probability_limit = 1;
	if (probability_shift < 1)
		probability_shift = 1;
	if (probability_shift > 32)
		probability_shift = 32;
	memset(state, 0, sizeof(*state));
	state->contextlength = contextlength;
	state->contextlength_mask = (1u << state->contextlength) - 1;
	state->probability_base = probability_base;
	state->probability_step = probability_step;
	state->probability_limit = probability_limit;
	state->probability_shift = probability_shift;
	unsigned int contextsize2 = 2u << state->contextlength;
	unsigned int stream = 0;
	const unsigned char* coded = coded_start;
	size_t maxbytes = coded_maxbytes;
	size_t cursor = 0;
	state->status = LHNETCODEC_STATUS_DECODE_CORRUPT;
	// when starting to decode the stream we need to read enough bytes to be
	// able to make decisions
	if (coded)
	{
		stream = (stream << 8) & 0xFFFF00; if (cursor < maxbytes) stream |= coded[cursor] & 0xFF; cursor++;
		stream = (stream << 8) & 0xFFFF00; if (cursor < maxbytes) stream |= coded[cursor] & 0xFF; cursor++;
		stream = (stream << 8) & 0xFFFF00; if (cursor < maxbytes) stream |= coded[cursor] & 0xFF; cursor++;
	}
	state->range_low = 0;
	state->range_high = 0xFFFFFF;
	state->range_stream = stream;
	state->coded_start = coded_start;
	state->coded_maxbytes = maxbytes;
	state->coded_cursor_any_byte = cursor;
	state->status = LHNETCODEC_STATUS_OK;
	return state->status;
}

LHNETCODEC_enum LHNETCODEC_AdaptiveCoder_Encode_Number(LHNETCODEC_AdaptiveCoder_State* state, unsigned int firstcontext, unsigned int numcontexts, uint64_t number)
{
	if ((unsigned int)firstcontext + numcontexts > LHNETCODEC_MAXCONTEXTRANGE)
	{
		state->status = LHNETCODEC_STATUS_INVALID_PARAMETER;
		return state->status;
	}
	int position;
	unsigned int low = state->range_low;
	unsigned int high = state->range_high;
	unsigned char* contexts = state->contexts;
	unsigned int probability_base = state->probability_base;
	unsigned int probability_step = state->probability_step;
	unsigned int probability_limit = state->probability_limit;
	unsigned int probability_shift = state->probability_shift;
	unsigned int probability_base2 = probability_base * 2;
	unsigned int contextlength = state->contextlength;
	unsigned int contextsize2 = 2u << contextlength;
	unsigned int contextmask = (1u << contextlength) - 1;
	size_t cursor = state->coded_cursor_any_byte;
	size_t cursor_trimmed = state->coded_cursor_trimmed;
	unsigned char* coded = state->coded_start;
	unsigned char* histogram = state->histogram;
	size_t maxbytes = state->coded_maxbytes;
	if (state->status != LHNETCODEC_STATUS_OK)
		return state->status;
	unsigned char* h;
	position = numcontexts - 1;
	h = histogram + (firstcontext + position) * contextsize2;
	for (; position >= 0; position--, h -= contextsize2)
	{
		unsigned int c = firstcontext + position;
		unsigned int context = contexts[c] * 2;
		unsigned int probability_0 = h[context];
		unsigned int probability_1 = h[context + 1];
		unsigned int split = low + (unsigned int)
			((uint64_t)(high - low) * (probability_0 + probability_base) /
			(probability_0 + probability_1 + probability_base2));
		unsigned int bit = (number >> position) & 1;
		if (bit)
		{
			low = split;
			probability_1 += probability_step;
			// if the probability limit is reached, apply the shift to keep the
			// coder from becoming too specialized, and to keep the histogram
			// from overflowing
			if (probability_1 >= probability_limit)
			{
				probability_0 >>= probability_shift;
				probability_1 >>= probability_shift;
				h[context] = probability_0;
			}
			h[context + 1] = probability_1;
			context++;
		}
		else
		{
			high = split - 1;
			probability_0 += probability_step;
			// if the probability limit is reached, apply the shift to keep the
			// coder from becoming too specialized, and to keep the histogram
			// from overflowing
			if (probability_0 >= probability_limit)
			{
				probability_0 >>= probability_shift;
				probability_1 >>= probability_shift;
				h[context + 1] = probability_1;
			}
			h[context] = probability_0;
		}
		context &= contextmask;
		contexts[c] = context;
		ENCODE_CHECKRANGE;
	}
	state->range_low = low;
	state->range_high = high;
	state->coded_cursor_any_byte = cursor;
	state->coded_cursor_trimmed = cursor_trimmed;
	return state->status;
}

LHNETCODEC_enum LHNETCODEC_AdaptiveCoder_Encode_Flush(LHNETCODEC_AdaptiveCoder_State* state)
{
	unsigned int low = state->range_low;
	size_t cursor = state->coded_cursor_any_byte;
	size_t cursor_trimmed = state->coded_cursor_trimmed;
	unsigned char* coded = state->coded_start;
	size_t maxbytes = state->coded_maxbytes;
	if (state->status != LHNETCODEC_STATUS_OK)
		return state->status;
	// flush the last non-zero bytes of low to the stream
	while (low)
	{
		if (cursor >= maxbytes)
		{
			state->status = LHNETCODEC_STATUS_ENCODE_FULL;
			break;
		}
		unsigned int b = (low >> 16) & 0xFF;
		if (coded)
			coded[cursor] = b;
		cursor++;
		if (b)
			cursor_trimmed = cursor;
		low = (low << 8) & 0xFFFFFF;
	}
	state->coded_cursor_any_byte = cursor;
	state->coded_cursor_trimmed = cursor_trimmed;
	return state->status;
}

uint64_t LHNETCODEC_AdaptiveCoder_Decode_Number(LHNETCODEC_AdaptiveCoder_State* state, unsigned int firstcontext, unsigned int numcontexts)
{
	if ((unsigned int)firstcontext + numcontexts > LHNETCODEC_MAXCONTEXTRANGE)
	{
		state->status = LHNETCODEC_STATUS_INVALID_PARAMETER;
		return 0;
	}
	if (state->status != LHNETCODEC_STATUS_OK || !state->coded_start)
		return 0;
	uint64_t number = 0;
	int position;
	unsigned int low = state->range_low;
	unsigned int high = state->range_high;
	unsigned int stream = state->range_stream;
	unsigned char* contexts = state->contexts;
	unsigned int probability_base = state->probability_base;
	unsigned int probability_step = state->probability_step;
	unsigned int probability_limit = state->probability_limit;
	unsigned int probability_shift = state->probability_shift;
	unsigned int probability_base2 = probability_base * 2;
	unsigned int contextlength = state->contextlength + 1;
	unsigned int contextsize2 = 2u << contextlength;
	unsigned int contextmask = state->contextlength_mask;
	size_t cursor = state->coded_cursor_any_byte;
	unsigned char* coded = state->coded_start;
	unsigned char* histogram = state->histogram;
	size_t maxbytes = state->coded_maxbytes;
	if (state->status != LHNETCODEC_STATUS_OK)
		return state->status;
	unsigned char* h;
	for (position = numcontexts - 1, h = histogram + position * contextsize2; position >= 0; position--, h -= contextsize2)
	{
		unsigned int c = firstcontext + position;
		unsigned int context = contexts[c] * 2;
		unsigned int probability_0 = h[context];
		unsigned int probability_1 = h[context + 1];
		unsigned int split = low +
			(unsigned int)(
				(uint64_t)(high - low) *
				(probability_0 + probability_base
					) / (
						probability_0 + probability_1 + probability_base2
						));
		unsigned int bit = stream >= split ? 1 : 0;
		if (bit)
		{
			low = split;
			probability_1 += probability_step;
			// if the probability limit is reached, apply the shift to keep the
			// coder from becoming too specialized, and to keep the histogram
			// from overflowing
			if (probability_1 >= probability_limit)
			{
				probability_0 >>= probability_shift;
				probability_1 >>= probability_shift;
				h[context] = probability_0;
			}
			h[context + 1] = probability_1;
			context++;
			number |= 1ull << position;
		}
		else
		{
			high = split - 1;
			probability_0 += probability_step;
			// if the probability limit is reached, apply the shift to keep the
			// coder from becoming too specialized, and to keep the histogram
			// from overflowing
			if (probability_0 >= probability_limit)
			{
				probability_0 >>= probability_shift;
				probability_1 >>= probability_shift;
				h[context + 1] = probability_1;
			}
			h[context] = probability_0;
		}
		context &= contextmask;
		contexts[c] = context;
		DECODE_CHECKRANGE;
	}
	state->range_low = low;
	state->range_high = high;
	state->range_stream = stream;
	state->coded_cursor_any_byte = cursor;
	return state->status;
}



/// when encoding using this block coder (where the histogram is stored
/// explicitly in the stream before the main text) we need to build a histogram
/// before we can actually encode the numbers, but we also must abide by a limit
/// on the size of the encoded text, so we have to flush the pending numbers to
/// the encoder to get that exact size, and estimate new pending numbers as very
/// poorly compressed, as the buffer fills up we will flush more frequently and
/// get an exact fit on the last several numbers, at a higher cost in CPU toward
/// the end of course, the caller can adjust the flush threshold.
/// 
/// we also have to support rollbacks of the encoder state because it is often
/// desired to only included a data structure in its entirety or defer it to the
/// next packet, so if writes failed partway through the structure, the caller
/// can use the rollback functionality and try appending other data or simply
/// finish building the packet and send it, next time the data structure may be
/// able to fit.
///
/// the advantages of using the block coder are the following:
/// * free CRC functionality - the histogram is written at the start of the
///   stream, and is decremented as the decoding proceeds, so if reading numbers
///   from the stream goes past the number counted in the histogram, the stream
///   is clearly corrupt - unfortunately you won't know this until you decode
///   that far, so if you can not handle corrupt numbers being read, it won't be
///   very useful to know this by that point.
/// * known length - the histogram contains an exact count of every bit value so
///   when it reaches zero we know to stop parsing at that moment.
/// * efficient encoding - the PPM is prone to mispredictions because it has no
///   prior knowledge of what is to come, block coding knows in advance what
///   values will appear, so compression ratios are very good from the start,
///   and become near-perfect at the end of the block when the contexts start to
///   have zero-probability of some values and no data needs to be read to
///   finish decoding those cases.  this comes at a clear cost for storing the
///   histogram in the packet, but if you already needed to include a CRC this
///   is similar in size to a regular hash.
/// 
/// disadvantages:
/// * potentially very high CPU cost when encoding - CheckFullness is cheap on a
///   block that is mostly empty, but has exponentially higher cost when called
///   repeatedly on a mostly-full block, as it has to call Flush frequently to
///   determine if there is still a little bit of space left; using a higher
///   value for fullness_threshold greatly improves its CPU cost.
typedef struct LHNETCODEC_BlockCoder_Encode_State
{
	/// counts of how many times a 0 or 1 occurs after each possible context,
	/// this decrements when decoding each bit, so the buffer is fully decoded
	/// when remaining_bits == 0 and remaining_bits is a sum of histogram[].
	unsigned int histogram[LHNETCODEC_MAXHISTOGRAMSIZE];
	/// the first time a new histogram context is encountered, the values for
	/// it are encoded, this avoids the need for a very large header in practice
	/// and avoids all zero values (since if the context is never reached, it is
	/// never stored)
	unsigned char histogramcontextsused[LHNETCODEC_MAXHISTOGRAMCONTEXTS];
	/// store the previous N bits for each of the 64 bits in a 64bit number
	/// being encoded/decoded, reading/writing smaller sizes will only advance
	/// the contexts for the used bits.
	unsigned char contexts[LHNETCODEC_MAXCONTEXTRANGE];
	/// context length in bits
	unsigned int contextlength;
	/// context_mask = (1u << contextlength) - 1;
	unsigned int contextlength_mask;
	/// current coded data buffer which stores the resulting stream, can be NULL
	/// if merely being used for simulations.
	unsigned char* coded_start;
	/// size of the coded working buffer in bytes, when coded_cursor_byte
	/// reaches this value, no more bits can be appended and an error status is
	/// returned.
	size_t coded_maxbytes;
	/// current position to write the next byte in the coded data buffer, this
	/// is for any byte value we may want to put there
	size_t coded_cursor_any_byte;
	/// the trimmed length only goes to the last non-zero byte
	size_t coded_cursor_trimmed;
	// this indicates a rollback has occurred but the encode has not been run,
	// again yet, so coded_cursor is a valid number for estimating if appending
	// numbers, but the actual data stream needs to be regenerated on flush.
	unsigned int finish_needs_flush;
	/// which version of the format we are encoding
	unsigned int coded_version;
	/// indicates the caller wants the versioned magic number to be included
	unsigned int coded_store_magic_number;
	/// indicates the value of cursor at the end of the header
	/// (after versioned magic number, contextlength, histogram)
	size_t coded_header_length;
	/// current operating range in the range coder
	unsigned int range_low;
	/// current operating range in the range coder
	unsigned int range_high;
	// firstcontext of each number to encode
	uint8_t numberfirstcontext[LHNETCODEC_BLOCKCODER_REPLAY_NUMBERS];
	// numcontexts of each number to encode
	uint8_t numbernumcontexts[LHNETCODEC_BLOCKCODER_REPLAY_NUMBERS];
	// value of each number
	uint64_t numbers[LHNETCODEC_BLOCKCODER_REPLAY_NUMBERS];
	// number of bits used in the buffer so far
	unsigned int number_totalwidth;
	// how many numbers are stored
	unsigned int number_count;
	// value of number_totalwidth at the time of last encode attempt, the delta
	// between number_totalwidth and this is the outstanding (at risk) data that
	// has not been encoded yet, so we use a conservative estimate on the size
	// of that data.
	unsigned int number_encodedwidth;

	/// current operational status of the coder
	LHNETCODEC_enum status;
}
LHNETCODEC_BlockCoder_Encode_State;

/// Initializes the state for a new encoding session.
void LHNETCODEC_BlockCoder_Encode_Init(LHNETCODEC_BlockCoder_Encode_State* state, unsigned char *coded_start, size_t coded_maxbytes, LHNETCODEC_enum version, unsigned int include_versioned_magic_number, unsigned int contextlength)
{
	memset(state, 0, sizeof(*state));
	state->coded_start = coded_start;
	state->coded_maxbytes = coded_maxbytes;
	state->coded_cursor_any_byte = 0;
	state->coded_cursor_trimmed = 0;
	state->coded_version = version;
	// we use this just to guarantee that Flush is called by Finish if Number
	// does not trigger a Flush
	state->finish_needs_flush = 1;
	state->status = LHNETCODEC_STATUS_UNINITIALIZED;
	switch (version) {
	case LHNETCODEC_BLOCKCODER_VERSION_0:
		if (contextlength > LHNETCODEC_MAXCONTEXTLENGTH)
			contextlength = LHNETCODEC_MAXCONTEXTLENGTH;
		state->contextlength = contextlength;
		state->contextlength_mask = (1u << contextlength) - 1;
		// just bump the cursor so we don't underestimate the bytes we'll
		// write in Flush.
		state->coded_store_magic_number = include_versioned_magic_number;
		if (state->coded_store_magic_number)
			state->coded_cursor_any_byte += 4;
		// add a byte for contextlength
		state->coded_cursor_any_byte++;
		// add some bytes for at least one histogram context
		state->coded_cursor_any_byte += 2;
		// this is only an estimate - Flush will determine the real length
		state->coded_header_length = state->coded_cursor_any_byte;
		state->coded_version = LHNETCODEC_BLOCKCODER_VERSION_0;
		state->status = LHNETCODEC_STATUS_OK;
		break;
	// no default case because we want the compiler to warn if a new version is
	// not handled here
	}
}

void LHNETCODEC_BlockCoder_GenerateHistogram(LHNETCODEC_BlockCoder_Encode_State* state)
{
	unsigned int* histogram = state->histogram;
	unsigned char* contexts = state->contexts;
	unsigned int contextlength = state->contextlength;
	unsigned int contextsize2 = 2u << contextlength;
	unsigned int contextmask = (1u << contextlength) - 1;
	unsigned char* numberfirstcontext = state->numberfirstcontext;
	unsigned char* numbernumcontexts = state->numbernumcontexts;
	uint64_t* numbers = state->numbers;
	unsigned int n;
	memset(state->histogram, 0, sizeof(state->histogram));
	memset(state->histogramcontextsused, 0, sizeof(state->histogramcontextsused));
	memset(state->contexts, 0, sizeof(state->contexts));
	for (n = 0; n < state->number_count; n++)
	{
		unsigned int position;
		uint64_t number = numbers[n];
		unsigned int firstcontext = numberfirstcontext[n];
		unsigned int numcontexts = numbernumcontexts[n];
		unsigned int* h;
		// we walk the position in opposite order here because it is slightly
		// easier and position doesn't matter for this purpose
		for (position = 0, h = histogram; position < numcontexts; position++, h += contextsize2)
		{
			unsigned int c = firstcontext + position;
			unsigned int context = contexts[c] * 2 + ((number >> position) & 1);
			h[context]++;
			contexts[c] = context & contextmask;
		}
	}
}

/// attempts to encode the buffered numbers, this is not expected to fail when
/// called by AppendNumber which uses conservative estimates
LHNETCODEC_enum LHNETCODEC_BlockCoder_Encode_Flush(LHNETCODEC_BlockCoder_Encode_State* state)
{
	// reset all encoder state except the numbers buffer
	state->number_encodedwidth = 0;
	state->range_low = 0;
	state->range_high = 0xFFFFFF;
	state->coded_cursor_any_byte = 0;
	state->coded_cursor_trimmed = 0;
	state->finish_needs_flush = 0;
	// keep local copies of the variables so the compiler knows the values are
	// completely under our control here
	unsigned char* numberfirstcontext = state->numberfirstcontext;
	unsigned char* numbernumcontexts = state->numbernumcontexts;
	uint64_t* numbers = state->numbers;
	unsigned int low = state->range_low;
	unsigned int high = state->range_high;
	unsigned char* contexts = state->contexts;
	unsigned int contextlength = state->contextlength;
	unsigned int contextlength1 = contextlength + 1;
	unsigned int contextsize2 = 2 << contextlength;
	unsigned int histogramsize = LHNETCODEC_MAXCONTEXTRANGE * contextsize2;
	unsigned int contextmask = (1 << contextlength) - 1;
	size_t cursor = state->coded_cursor_any_byte;
	size_t cursor_trimmed = state->coded_cursor_trimmed;
	unsigned char* coded = state->coded_start;
	unsigned int* histogram = state->histogram;
	size_t maxbytes = state->coded_maxbytes;
	// calculate the histogram before encoding with it
	LHNETCODEC_BlockCoder_GenerateHistogram(state);
	// reset contexts again because we need to start fresh
	memset(state->contexts, 0, sizeof(state->contexts));
	// write the magic number if the caller wants to
	if (state->coded_store_magic_number)
	{
		switch (state->coded_version)
		{
		case LHNETCODEC_BLOCKCODER_VERSION_0:
			if (cursor + 4 > maxbytes)
			{
				state->status = LHNETCODEC_STATUS_ENCODE_FULL;
				return state->status;
			}
			if (coded)
				memcpy(coded + cursor, "LHN0", 4);
			cursor += 4;
			cursor_trimmed = cursor;
			break;
		}
	}
	unsigned int flags = state->contextlength;
	if (coded)
		coded[cursor] = flags;
	cursor++;
	// if the flags is not zero, include it in the trimmed length
	if (flags)
		cursor_trimmed = cursor;
	state->coded_header_length = cursor_trimmed;
	// now encode the numbers with the diminishing histogram
	unsigned int n;
	for (n = 0; n < state->number_count; n++)
	{
		int position;
		uint64_t number = numbers[n];
		unsigned int firstcontext = numberfirstcontext[n];
		unsigned int numcontexts = numbernumcontexts[n];
		unsigned int* h;
		state->number_encodedwidth += numcontexts;
		position = numcontexts - 1;
		h = histogram + (firstcontext + position) * contextsize2;
		for (; position >= 0; position--, h -= contextsize2)
		{
			unsigned int c = firstcontext + position;
			unsigned int context = contexts[c];
			if (!state->histogramcontextsused[(c << contextlength) + context])
			{
				// first time we've seen this context, so we need to code the
				// probabilities
				int side;
				state->histogramcontextsused[(c << contextlength) + context] = 1;
				for (side = 0; side < 2; side++)
				{
					unsigned int count = h[context * 2 + side];
					unsigned int c = count;
					for (;;c >>= LHNETCODEC_BLOCKCODER_TABLEBITS)
					{
						unsigned int value = c;
						unsigned int extend = value > 0 ? 1 : 0;
						ENCODE_NUMBER(extend, 2);
						if (!extend)
							break;
						value &= (1 << LHNETCODEC_BLOCKCODER_TABLEBITS) - 1;
						ENCODE_NUMBER(value, 1 << LHNETCODEC_BLOCKCODER_TABLEBITS);
					}
				}
			}
			context *= 2;
			unsigned int probability_0 = h[context];
			unsigned int probability_1 = h[context + 1];
			MAKE_PROBABILITY_SAFE;
			unsigned int split = low + (unsigned int)((uint64_t)(high - low) * probability_0 / (probability_0 + probability_1));
			unsigned int bit = (number >> position) & 1;
			if (bit)
			{
				low = split;
				context++;
			}
			else
				high = split - 1;
			h[context]--;
			context &= contextmask;
			contexts[c] = context;
			ENCODE_CHECKRANGE;
		}
	}
	// now flush the remaining bytes until low becomes zero, this may write up
	// to 3 bytes
	while (low)
	{
		unsigned int b = (low >> 16) & 0xFF;
		if (coded)
			coded[cursor] = b;
		cursor++;
		if (b)
			cursor_trimmed = cursor;
		low = (low << 8) & 0xFFFFFF;
	}
	// at this point, the entire block has been encoded, and histogram[] is
	// entirely zero (as we decremented while encoding), encodedwidth is equal
	// to total_width
	state->coded_cursor_any_byte = cursor;
	state->coded_cursor_trimmed = cursor_trimmed;
	return state->status;
}

/// append a single number with the specified number of bits to the block, this
/// can fail if too many numbers are written to the block, or if flush fails to
/// fit the compressed data.
/// 
/// never_flush is advisable when writing groups of several numbers if you
/// inte if appending several numbers at once and relying on
/// rollback to revert all of them as a unit if the encoded size is too big.
LHNETCODEC_enum LHNETCODEC_BlockCoder_Encode_CheckFullness(LHNETCODEC_BlockCoder_Encode_State* state, unsigned int fullness_threshold)
{
	if (state->status != LHNETCODEC_STATUS_OK)
		return state->status;
	// finishing_bytes has to account for up to 3 flush bytes at the end, and
	// an increase in number of prefix bytes (at least 1), and one more
	// 64bit number with the worst imaginable compression ratio, so this is a
	// conservative estimate as typically numbers compress at least a little,
	// but we want to be absolutely sure it fits.
	const unsigned int finishing_bytes = 12;
	// if the caller wants us to say the buffer is full a bit early, do so
	fullness_threshold += finishing_bytes;
	if (state->coded_cursor_any_byte + fullness_threshold > state->coded_maxbytes)
	{
		state->status = LHNETCODEC_STATUS_ENCODE_FULL;
		return state->status;
	}
	// this causes Flush to be called more often as the block fills up,
	// which refines our estimated remaining space, but at high cost.
	unsigned int outstanding_bytes = finishing_bytes
		+ ((state->number_totalwidth - state->number_encodedwidth + 7) >> 3);
	if (state->coded_cursor_any_byte + outstanding_bytes > state->coded_maxbytes)
		LHNETCODEC_BlockCoder_Encode_Flush(state);
	return state->status;
}

/// append a single number with the specified number of bits to the block,
/// without checking if they will fit in the encoded block size limit (see
/// CheckFullness and Rollback).
/// 
/// this can fail if the maximum number of replayable numbers has been reached
/// in the block, by which point you have probably greatly exceeded the encoded
/// size limit and should have called CheckFullness already.
LHNETCODEC_enum LHNETCODEC_BlockCoder_Encode_Number(LHNETCODEC_BlockCoder_Encode_State* state, unsigned int firstcontext, unsigned int numcontexts, uint64_t number)
{
	if (state->status != LHNETCODEC_STATUS_OK)
		return state->status;
	if (state->number_count >= LHNETCODEC_BLOCKCODER_REPLAY_NUMBERS)
	{
		state->status = LHNETCODEC_STATUS_ENCODE_FULL;
		return state->status;
	}
	state->numberfirstcontext[state->number_count] = firstcontext;
	state->numbernumcontexts[state->number_count] = numcontexts;
	state->numbers[state->number_count] = number;
	state->number_totalwidth += numcontexts;
	state->number_count++;
	return state->status;
}

typedef struct LHNETCODEC_BlockCoder_Encode_RollbackPosition
{
	unsigned int number_count;
	size_t cursor;
}
LHNETCODEC_BlockCoder_Encode_RollbackPosition;

/// gets a rollback position for undoing Number operations when a group of
/// several related numbers will only partially fit - for example a data
/// structure that you want to be included in its entirety or not at all.
LHNETCODEC_BlockCoder_Encode_RollbackPosition LHNETCODEC_BlockCoder_Encode_GetRollbackPosition(const LHNETCODEC_BlockCoder_Encode_State* state)
{
	LHNETCODEC_BlockCoder_Encode_RollbackPosition p;
	p.number_count = state->number_count;
	p.cursor = state->coded_cursor_any_byte;
	return p;
}

LHNETCODEC_enum LHNETCODEC_BlockCoder_Encode_Rollback(LHNETCODEC_BlockCoder_Encode_State* state, LHNETCODEC_BlockCoder_Encode_RollbackPosition p)
{
	state->status = LHNETCODEC_STATUS_OK;
	while (state->number_count > p.number_count)
		state->number_totalwidth -= state->numbernumcontexts[--state->number_count];
	state->coded_cursor_any_byte = p.cursor;
	// while coded_cursor is sufficient for estimating further Number calls, the
	// coded data is corrupted, so Finish must call Flush
	state->finish_needs_flush = 1;
	return state->status;
}

LHNETCODEC_enum LHNETCODEC_BlockCoder_Encode_Finish(LHNETCODEC_BlockCoder_Encode_State* state, unsigned char **coded_data, uint64_t *coded_length)
{
	if (state->finish_needs_flush)
		LHNETCODEC_BlockCoder_Encode_Flush(state);
	if (coded_data)
		*coded_data = state->coded_start;
	if (coded_length)
		*coded_length = state->coded_cursor_trimmed;
	return state->status;
}

typedef struct LHNETCODEC_BlockCoder_Decode_State
{
	/// counts of how many times a 0 or 1 occurs after each possible context,
	/// this decrements when decoding each bit, so the buffer is fully decoded
	/// when remaining_bits == 0 and remaining_bits is a sum of histogram[].
	unsigned int histogram[LHNETCODEC_MAXHISTOGRAMSIZE];
	/// the first time a new histogram context is encountered, the values for
	/// it are encoded, this avoids the need for a very large header in practice
	/// and avoids all zero values (since if the context is never reached, it is
	/// never stored)
	unsigned char histogramcontextsused[LHNETCODEC_MAXHISTOGRAMCONTEXTS];
	/// store the previous N bits for each of the 64 bits in a 64bit number
	/// being encoded/decoded, reading/writing smaller sizes will only advance
	/// the contexts for the used bits.
	unsigned char contexts[LHNETCODEC_MAXCONTEXTRANGE];
	/// context size in bits
	unsigned int contextlength;
	/// context_mask = (1u << contextlength) - 1;
	unsigned int context_mask;
	/// current coded data buffer which stores the resulting stream, can be NULL
	/// if merely being used for simulations.
	const unsigned char* coded_start;
	/// size of the coded working buffer in bytes, when coded_cursor
	/// reaches this value, no more bits can be appended and an error status is
	/// returned.
	size_t coded_maxbytes;
	/// current position in the coded data buffer, in bytes.
	size_t coded_cursor;
	/// which version of the format we are decoding
	unsigned int coded_version;
	/// indicates the block has a magic number included
	unsigned int coded_store_magic_number;
	/// indicates the value of cursor at the end of the header
	/// (after versioned magic number, contextlength, histogram)
	size_t coded_header_length;
	/// current operating range in the range coder
	unsigned int range_low;
	/// current operating range in the range coder
	unsigned int range_high;
	/// (decode only) fractional value being streamed in
	unsigned int range_stream;
	/// sum of histogram[], used for detecting corrupt blocks
	unsigned int total_width;

	/// current operational status of the coder
	LHNETCODEC_enum status;
}
LHNETCODEC_BlockCoder_Decode_State;

/// Initializes the state for a new encoding session.
LHNETCODEC_enum LHNETCODEC_BlockCoder_Decode_Init(LHNETCODEC_BlockCoder_Decode_State* state, const unsigned char* coded_start, size_t coded_length, LHNETCODEC_enum version_hint, unsigned int read_versioned_magic_number)
{
	const unsigned char* coded = coded_start;
	size_t maxbytes = coded_length;
	size_t cursor = 0;
	memset(state, 0, sizeof(*state));
	// default to this status until everything looks ok
	state->status = LHNETCODEC_STATUS_DECODE_CORRUPT;
	state->coded_version = version_hint;
	if (read_versioned_magic_number)
	{
		if (cursor + 4 <= maxbytes
			&& !memcmp(coded + cursor, "LHN0", 4))
		{
			state->coded_version = LHNETCODEC_BLOCKCODER_VERSION_0;
			cursor += 4;
			state->coded_store_magic_number = 1;
		}
	}
	switch (state->coded_version)
	{
	case LHNETCODEC_BLOCKCODER_VERSION_0:
		{
			unsigned int total_width = 0;
			unsigned int stream = 0;
			unsigned int flags = cursor < maxbytes ? coded[cursor] : 0;
			cursor++;
			state->contextlength = flags & 0x1F;
			if (state->contextlength > LHNETCODEC_MAXCONTEXTLENGTH)
			{
				state->status = LHNETCODEC_STATUS_DECODE_CORRUPT;
				return state->status;
			}
			state->context_mask = (1u << state->contextlength) - 1;
			state->total_width = total_width;
			state->coded_header_length = cursor;
			// when starting to decode the stream we need to read enough bytes to be
			// able to make decisions
			stream = 0;
			stream = (stream << 8) & 0xFFFF00; if (cursor < maxbytes) stream |= coded[cursor] & 0xFF; cursor++;
			stream = (stream << 8) & 0xFFFF00; if (cursor < maxbytes) stream |= coded[cursor] & 0xFF; cursor++;
			stream = (stream << 8) & 0xFFFF00; if (cursor < maxbytes) stream |= coded[cursor] & 0xFF; cursor++;
			state->coded_start = coded;
			state->coded_maxbytes = maxbytes;
			state->coded_cursor = cursor;
			state->range_low = 0;
			state->range_high = 0xFFFFFF;
			state->range_stream = stream;
			state->status = LHNETCODEC_STATUS_OK;
		}
		break;
	}
	return state->status;
}

/// attempts to decode one number as a set of bits
uint64_t LHNETCODEC_BlockCoder_Decode_Number(LHNETCODEC_BlockCoder_Decode_State* state, unsigned int firstcontext, unsigned int numcontexts)
{
	uint64_t number = 0;
	int position;
	unsigned int low = state->range_low;
	unsigned int high = state->range_high;
	unsigned int stream = state->range_stream;
	unsigned char* contexts = state->contexts;
	unsigned int contextlength = state->contextlength;
	unsigned int contextlength1 = contextlength;
	unsigned int contextsize2 = 1 << contextlength1;
	unsigned int contextmask = state->context_mask;
	size_t cursor = state->coded_cursor;
	const unsigned char* coded = state->coded_start;
	unsigned int* histogram = state->histogram;
	size_t maxbytes = state->coded_maxbytes;
	unsigned int total_width = state->total_width;
	if (state->status != LHNETCODEC_STATUS_OK)
		return 0;
	if (state->total_width == 0)
	{
		state->status = LHNETCODEC_STATUS_DECODE_EOF;
		return 0;
	}
	// if the histogram is exhausted partway through a number, the stream is
	// corrupt
	if (total_width < numcontexts)
	{
		state->status = LHNETCODEC_STATUS_DECODE_CORRUPT;
		return 0;
	}
	unsigned int* h;
	position = numcontexts - 1;
	h = histogram + (firstcontext + position) * contextsize2;
	for (; position >= 0; position--, h -= contextsize2)
	{
		unsigned int c = firstcontext + position;
		unsigned int context = contexts[c];
		if (!state->histogramcontextsused[(c << contextlength) + context])
		{
			// first time we've seen this context, so we need to code the
			// probabilities
			int side;
			state->histogramcontextsused[(c << contextlength) + context] = 1;
			for (side = 0; side < 2; side++)
			{
				unsigned int count = 0;
				unsigned int shift = 0;
				for (;;)
				{
					unsigned int extend;
					DECODE_NUMBER(extend, 2);
					if (!extend)
						break;
					unsigned int value;
					DECODE_NUMBER(value, 1 << LHNETCODEC_BLOCKCODER_TABLEBITS);
					count += value << shift;
					shift += LHNETCODEC_BLOCKCODER_TABLEBITS;
				}
				h[context * 2 + side] = count;
			}
		}
		context *= 2;
		unsigned int probability_0 = h[context];
		unsigned int probability_1 = h[context + 1];
		MAKE_PROBABILITY_SAFE;
		unsigned int split = low + (unsigned int)((uint64_t)(high - low) * probability_0 / (probability_0 + probability_1));
		unsigned int bit = stream >= split ? 1 : 0;
		if (bit)
		{
			low = split;
			context++;
			number |= 1ull << position;
		}
		else
			high = split - 1;
		h[context]--;
		context &= contextmask;
		contexts[c] = context;
		DECODE_CHECKRANGE;
	}
	state->range_low = low;
	state->range_high = high;
	state->range_stream = stream;
	state->coded_cursor = cursor;
	state->total_width = total_width;
	return state->status;
}




/// this is a variant of the above BlockCoder that uses a smaller encoding for
/// the histogram but can not take advantage of diminishing probabilities, it
/// is faster because it does not use integer division after histogram build.
typedef struct LHNETCODEC_BlockCoderFixedPoint_Encode_State
{
	/// counts of how many times a 0 or 1 occurs after each possible context,
	/// this decrements when decoding each bit, so the buffer is fully decoded
	/// when remaining_bits == 0 and remaining_bits is a sum of histogram[].
	unsigned int histogram[LHNETCODEC_MAXHISTOGRAMSIZE];
	/// the first time a new histogram context is encountered, the values for
	/// it are encoded, this avoids the need for a very large header in practice
	/// and avoids all zero values (since if the context is never reached, it is
	/// never stored)
	unsigned char histogramcontextsused[LHNETCODEC_MAXHISTOGRAMCONTEXTS];
	/// each value here is a fixed point representation of the ratio between the
	/// two context probabilities, which is what we actually write to the stream
	/// instead of the raw numbers
	unsigned short probability_splits[LHNETCODEC_MAXHISTOGRAMCONTEXTS];
	/// store the previous N bits for each of the 64 bits in a 64bit number
	/// being encoded/decoded, reading/writing smaller sizes will only advance
	/// the contexts for the used bits.
	unsigned char contexts[LHNETCODEC_MAXCONTEXTRANGE];
	/// context size in bits
	unsigned int contextlength;
	/// context_mask = (1u << contextlength) - 1;
	unsigned int context_mask;
	/// current coded data buffer which stores the resulting stream, can be NULL
	/// if merely being used for simulations.
	unsigned char* coded_start;
	/// size of the coded working buffer in bytes, when coded_cursor_byte
	/// reaches this value, no more bits can be appended and an error status is
	/// returned.
	size_t coded_maxbytes;
	/// current position to write the next byte in the coded data buffer, this
	/// is for any byte value we may want to put there
	size_t coded_cursor_any_byte;
	/// the trimmed length only goes to the last non-zero byte
	size_t coded_cursor_trimmed;
	// this indicates a rollback has occurred but the encode has not been run,
	// again yet, so coded_cursor is a valid number for estimating if appending
	// numbers, but the actual data stream needs to be regenerated on flush.
	unsigned int finish_needs_flush;
	/// which version of the format we are encoding
	unsigned int coded_version;
	/// indicates the caller wants the versioned magic number to be included
	unsigned int coded_store_magic_number;
	/// indicates the value of cursor at the end of the header
	/// (after versioned magic number, contextlength, histogram)
	size_t coded_header_length;
	/// current operating range in the range coder
	unsigned int range_low;
	/// current operating range in the range coder
	unsigned int range_high;
	// firstcontext of each number to encode
	uint8_t numberfirstcontext[LHNETCODEC_BLOCKCODER_REPLAY_NUMBERS];
	// numcontexts of each number to encode
	uint8_t numbernumcontexts[LHNETCODEC_BLOCKCODER_REPLAY_NUMBERS];
	// value of each number
	uint64_t numbers[LHNETCODEC_BLOCKCODER_REPLAY_NUMBERS];
	// number of bits used in the buffer so far
	// 0 <= count_width <= LHNETCODEC_BLOCKCODER_BLOCK_BITS
	unsigned int number_totalwidth;
	// how many numbers are stored
	// 0 <= count_values <= LHNETCODEC_BLOCKCODER_BLOCK_BITS
	unsigned int number_count;
	// value of number_totalwidth at the time of last encode attempt, the delta
	// between number_totalwidth and this is the outstanding (at risk) data that
	// has not been encoded yet, so we use a conservative estimate on the size
	// of that data.
	unsigned int number_encodedwidth;

	/// current operational status of the coder
	LHNETCODEC_enum status;
}
LHNETCODEC_BlockCoderFixedPoint_Encode_State;

/// Initializes the state for a new encoding session.
void LHNETCODEC_BlockCoderFixedPoint_Encode_Init(LHNETCODEC_BlockCoderFixedPoint_Encode_State* state, unsigned char* coded_start, size_t coded_maxbytes, LHNETCODEC_enum version, unsigned int include_versioned_magic_number, unsigned int contextlength)
{
	memset(state, 0, sizeof(*state));
	state->coded_start = coded_start;
	state->coded_maxbytes = coded_maxbytes;
	state->coded_cursor_any_byte = 0;
	state->coded_cursor_trimmed = 0;
	state->coded_version = version;
	// we use this just to guarantee that Flush is called by Finish if Number
	// does not trigger a Flush
	state->finish_needs_flush = 1;
	state->status = LHNETCODEC_STATUS_UNINITIALIZED;
	switch (version) {
	case LHNETCODEC_BLOCKCODERFIXEDPOINT_VERSION_0:
		if (contextlength > LHNETCODEC_MAXCONTEXTLENGTH)
			contextlength = LHNETCODEC_MAXCONTEXTLENGTH;
		state->contextlength = contextlength;
		state->context_mask = (1u << contextlength) - 1;
		// just bump the cursor so we don't underestimate the bytes we'll
		// write in Flush.
		state->coded_store_magic_number = include_versioned_magic_number;
		if (state->coded_store_magic_number)
			state->coded_cursor_any_byte += 4;
		// add a byte for contextlength
		state->coded_cursor_any_byte++;
		// add some bytes for at least one histogram context
		state->coded_cursor_any_byte += 2;
		// this is only an estimate - Flush will determine the real length
		state->coded_header_length = state->coded_cursor_any_byte;
		state->coded_version = LHNETCODEC_BLOCKCODERFIXEDPOINT_VERSION_0;
		state->status = LHNETCODEC_STATUS_OK;
		break;
		// no default case because we want the compiler to warn if a new version is
		// not handled here
	}
}

void LHNETCODEC_BlockCoderFixedPoint_GenerateHistogram(LHNETCODEC_BlockCoderFixedPoint_Encode_State* state)
{
	unsigned int* histogram = state->histogram;
	unsigned char* contexts = state->contexts;
	unsigned int contextlength = state->contextlength;
	unsigned int contextlength1 = contextlength + 1;
	unsigned int contextsize2 = 1 << contextlength1;
	unsigned int histogramcontexts = LHNETCODEC_MAXCONTEXTRANGE << contextlength;
	unsigned int contextmask = state->context_mask;
	unsigned char* numberfirstcontext = state->numberfirstcontext;
	unsigned char* numbernumcontexts = state->numbernumcontexts;
	uint64_t* numbers = state->numbers;
	unsigned int n;
	memset(state->histogram, 0, sizeof(state->histogram));
	memset(state->probability_splits, 0, sizeof(state->probability_splits));
	memset(state->contexts, 0, sizeof(state->contexts));
	for (n = 0; n < state->number_count; n++)
	{
		unsigned int position;
		uint64_t number = numbers[n];
		unsigned int firstcontext = numberfirstcontext[n];
		unsigned int numcontexts = numbernumcontexts[n];
		unsigned int* h;
		// we walk the position in opposite order here because it is slightly
		// easier and position doesn't matter for this purpose
		for (position = 0, h = histogram; position < numcontexts; position++, h += contextsize2)
		{
			unsigned int c = firstcontext + position;
			unsigned int context = contexts[c] * 2 + ((number >> position) & 1);
			h[context]++;
			contexts[c] = context & contextmask;
		}
	}
	unsigned int* h;
	for (n = 0, h = state->histogram; n < histogramcontexts / 2; n++, h += 2)
	{
		unsigned int probability_0 = h[0];
		unsigned int probability_1 = h[1];
		MAKE_PROBABILITY_SAFE;
		unsigned int split;
		if (probability_0 + probability_1 == 0)
			split = 1 << (LHNETCODEC_BLOCKCODERFIXEDPOINT_SPLITBITS - 1);
		else
		{
			split = (probability_0 << LHNETCODEC_BLOCKCODERFIXEDPOINT_SPLITBITS) / (probability_0 + probability_1);
			if (split == (1 << LHNETCODEC_BLOCKCODERFIXEDPOINT_SPLITBITS))
				split--;
		}
		state->probability_splits[n] = split;
	}
}

/// attempts to encode the buffered numbers, this is not expected to fail when
/// called by AppendNumber which uses conservative estimates
LHNETCODEC_enum LHNETCODEC_BlockCoderFixedPoint_Encode_Flush(LHNETCODEC_BlockCoderFixedPoint_Encode_State* state)
{
	// reset all encoder state except the numbers buffer
	state->number_encodedwidth = 0;
	state->range_low = 0;
	state->range_high = 0xFFFFFF;
	state->coded_cursor_any_byte = 0;
	state->coded_cursor_trimmed = 0;
	state->finish_needs_flush = 0;
	// keep local copies of the variables so the compiler knows the values are
	// completely under our control here
	unsigned char* numberfirstcontext = state->numberfirstcontext;
	unsigned char* numbernumcontexts = state->numbernumcontexts;
	uint64_t* numbers = state->numbers;
	unsigned int low = state->range_low;
	unsigned int high = state->range_high;
	unsigned char* contexts = state->contexts;
	unsigned int contextlength = state->contextlength;
	unsigned int histogramsize = (LHNETCODEC_MAXCONTEXTRANGE * 2) << contextlength;
	unsigned int histogramcontexts = LHNETCODEC_MAXCONTEXTRANGE << contextlength;
	unsigned int contextmask = state->context_mask;
	size_t cursor = state->coded_cursor_any_byte;
	size_t cursor_trimmed = state->coded_cursor_trimmed;
	unsigned char* coded = state->coded_start;
	unsigned int* histogram = state->histogram;
	size_t maxbytes = state->coded_maxbytes;
	// calculate the histogram before encoding with it
	LHNETCODEC_BlockCoderFixedPoint_GenerateHistogram(state);
	// reset contexts again because we need to start fresh
	memset(state->contexts, 0, sizeof(state->contexts));
	// write the magic number if the caller wants to
	if (state->coded_store_magic_number)
	{
		switch (state->coded_version)
		{
		case LHNETCODEC_BLOCKCODERFIXEDPOINT_VERSION_0:
			if (cursor + 4 > maxbytes)
			{
				state->status = LHNETCODEC_STATUS_ENCODE_FULL;
				return state->status;
			}
			if (coded)
				memcpy(coded + cursor, "LHF0", 4);
			cursor += 4;
			cursor_trimmed = cursor;
			break;
		}
	}
	unsigned int flags = state->contextlength;
	if (coded)
		coded[cursor] = flags;
	cursor++;
	// if the flags is not zero, include it in the trimmed length
	if (flags)
		cursor_trimmed = cursor;
	state->coded_header_length = cursor_trimmed;
	// now encode the numbers with the diminishing histogram
	unsigned int n;
	for (n = 0; n < state->number_count; n++)
	{
		int position;
		uint64_t number = numbers[n];
		unsigned int firstcontext = numberfirstcontext[n];
		unsigned int numcontexts = numbernumcontexts[n];
		state->number_encodedwidth += numcontexts;
		for (position = numcontexts - 1; position >= 0; position--)
		{
			unsigned int c = firstcontext + position;
			unsigned int context = contexts[c];
			unsigned int splitfixedpoint = state->probability_splits[(c << contextlength) + context];
			if (!state->histogramcontextsused[(c << contextlength) + context])
			{
				// first time we've seen this context, so we need to code the
				// probabilities
				state->histogramcontextsused[(c << contextlength) + context] = 1;
				ENCODE_NUMBER(splitfixedpoint, 1 << LHNETCODEC_BLOCKCODERFIXEDPOINT_SPLITBITS);
			}
			context *= 2;
			unsigned int split = low + (unsigned int)(((uint64_t)(high - low) * splitfixedpoint) >> LHNETCODEC_BLOCKCODERFIXEDPOINT_SPLITBITS);
			unsigned int bit = (number >> position) & 1;
			if (bit)
			{
				low = split;
				context++;
			}
			else
				high = split - 1;
			context &= contextmask;
			contexts[c] = context;
			ENCODE_CHECKRANGE;
		}
	}
	// now flush the remaining bytes until low becomes zero, this may write up
	// to 3 bytes
	while (low)
	{
		unsigned int b = (low >> 16) & 0xFF;
		if (coded)
			coded[cursor] = b;
		cursor++;
		if (b)
			cursor_trimmed = cursor;
		low = (low << 8) & 0xFFFFFF;
	}
	// at this point, the entire block has been encoded, and histogram[] is
	// entirely zero (as we decrmeneted while encoding), encodedwidth is equal
	// to total_width
	state->coded_cursor_any_byte = cursor;
	state->coded_cursor_trimmed = cursor_trimmed;
	return state->status;
}

/// append a single number with the specified number of bits to the block, this
/// can fail if too many numbers are written to the block, or if flush fails to
/// fit the compressed data.
/// 
/// never_flush is advisable when writing groups of several numbers if you
/// inte if appending several numbers at once and relying on
/// rollback to revert all of them as a unit if the encoded size is too big.
LHNETCODEC_enum LHNETCODEC_BlockCoderFixedPoint_Encode_CheckFullness(LHNETCODEC_BlockCoderFixedPoint_Encode_State* state, unsigned int fullness_threshold)
{
	if (state->status != LHNETCODEC_STATUS_OK)
		return state->status;
	// finishing_bytes has to account for up to 3 flush bytes at the end, and
	// an increase in number of prefix bytes (at least 1), and one more
	// 64bit number with the worst imaginable compression ratio, so this is a
	// conservative estimate as typically numbers compress at least a little,
	// but we want to be absolutely sure it fits.
	const unsigned int finishing_bytes = 12;
	// if the caller wants us to say the buffer is full a bit early, do so
	fullness_threshold += finishing_bytes;
	if (state->coded_cursor_any_byte + fullness_threshold > state->coded_maxbytes)
	{
		state->status = LHNETCODEC_STATUS_ENCODE_FULL;
		return state->status;
	}
	// this causes Flush to be called more often as the block fills up,
	// which refines our estimated remaining space, but at high cost.
	unsigned int outstanding_bytes = finishing_bytes
		+ ((state->number_totalwidth - state->number_encodedwidth + 7) >> 3);
	if (state->coded_cursor_any_byte + outstanding_bytes > state->coded_maxbytes)
		LHNETCODEC_BlockCoderFixedPoint_Encode_Flush(state);
	return state->status;
}

/// append a single number with the specified number of bits to the block,
/// without checking if they will fit in the encoded block size limit (see
/// CheckFullness and Rollback).
/// 
/// this can fail if the maximum number of replayable numbers has been reached
/// in the block, by which point you have probably greatly exceeded the encoded
/// size limit and should have called CheckFullness already.
LHNETCODEC_enum LHNETCODEC_BlockCoderFixedPoint_Encode_Number(LHNETCODEC_BlockCoderFixedPoint_Encode_State* state, unsigned int firstcontext, unsigned int numcontexts, uint64_t number)
{
	if (state->status != LHNETCODEC_STATUS_OK)
		return state->status;
	if (state->number_count >= LHNETCODEC_BLOCKCODER_REPLAY_NUMBERS)
	{
		state->status = LHNETCODEC_STATUS_ENCODE_FULL;
		return state->status;
	}
	state->numberfirstcontext[state->number_count] = firstcontext;
	state->numbernumcontexts[state->number_count] = numcontexts;
	state->numbers[state->number_count] = number;
	state->number_totalwidth += numcontexts;
	state->number_count++;
	return state->status;
}

typedef struct LHNETCODEC_BlockCoderFixedPoint_Encode_RollbackPosition
{
	size_t number_count;
	size_t cursor;
}
LHNETCODEC_BlockCoderFixedPoint_Encode_RollbackPosition;

/// gets a rollback position for undoing Number operations when a group of
/// several related numbers will only partially fit - for example a data
/// structure that you want to be included in its entirety or not at all.
LHNETCODEC_BlockCoderFixedPoint_Encode_RollbackPosition LHNETCODEC_BlockCoderFixedPoint_Encode_GetRollbackPosition(const LHNETCODEC_BlockCoderFixedPoint_Encode_State* state)
{
	LHNETCODEC_BlockCoderFixedPoint_Encode_RollbackPosition p;
	p.number_count = state->number_count;
	p.cursor = state->coded_cursor_any_byte;
	return p;
}

LHNETCODEC_enum LHNETCODEC_BlockCoderFixedPoint_Encode_Rollback(LHNETCODEC_BlockCoderFixedPoint_Encode_State* state, LHNETCODEC_BlockCoderFixedPoint_Encode_RollbackPosition p)
{
	state->status = LHNETCODEC_STATUS_OK;
	while (state->number_count > p.number_count)
		state->number_totalwidth -= state->numbernumcontexts[--state->number_count];
	state->coded_cursor_any_byte = p.cursor;
	// while coded_cursor is sufficient for estimating further Number calls, the
	// coded data is corrupted, so Finish must call Flush
	state->finish_needs_flush = 1;
	return state->status;
}

LHNETCODEC_enum LHNETCODEC_BlockCoderFixedPoint_Encode_Finish(LHNETCODEC_BlockCoderFixedPoint_Encode_State* state, unsigned char** coded_data, uint64_t* coded_length)
{
	if (state->finish_needs_flush)
		LHNETCODEC_BlockCoderFixedPoint_Encode_Flush(state);
	if (coded_data)
		*coded_data = state->coded_start;
	if (coded_length)
		*coded_length = state->coded_cursor_trimmed;
	return state->status;
}

typedef struct LHNETCODEC_BlockCoderFixedPoint_Decode_State
{
	/// counts of how many times a 0 or 1 occurs after each possible context,
	/// this decrements when decoding each bit, so the buffer is fully decoded
	/// when remaining_bits == 0 and remaining_bits is a sum of histogram[].
	unsigned int histogram[LHNETCODEC_MAXHISTOGRAMSIZE];
	/// the first time a new histogram context is encountered, the values for
	/// it are encoded, this avoids the need for a very large header in practice
	/// and avoids all zero values (since if the context is never reached, it is
	/// never stored)
	unsigned char histogramcontextsused[LHNETCODEC_MAXHISTOGRAMCONTEXTS];
	/// each value here is a fixed point representation of the ratio between the
	/// two context probabilities, which is what we actually write to the stream
	/// instead of the raw numbers
	unsigned short probability_splits[LHNETCODEC_MAXHISTOGRAMCONTEXTS];
	/// store the previous N bits for each of the 64 bits in a 64bit number
	/// being encoded/decoded, reading/writing smaller sizes will only advance
	/// the contexts for the used bits.
	unsigned char contexts[LHNETCODEC_MAXCONTEXTRANGE];
	/// context size in bits
	unsigned int contextlength;
	/// context_mask = (1u << contextlength) - 1;
	unsigned int context_mask;
	/// current coded data buffer which stores the resulting stream, can be NULL
	/// if merely being used for simulations.
	const unsigned char* coded_start;
	/// size of the coded working buffer in bytes, when coded_cursor
	/// reaches this value, no more bits can be appended and an error status is
	/// returned.
	size_t coded_maxbytes;
	/// current position in the coded data buffer, in bytes.
	size_t coded_cursor;
	/// which version of the format we are decoding
	unsigned int coded_version;
	/// indicates the block has a magic number included
	unsigned int coded_store_magic_number;
	/// indicates the value of cursor at the end of the header
	/// (after versioned magic number, contextlength, histogram)
	size_t coded_header_length;
	/// current operating range in the range coder
	unsigned int range_low;
	/// current operating range in the range coder
	unsigned int range_high;
	/// (decode only) fractional value being streamed in
	unsigned int range_stream;
	/// sum of histogram[], used for detecting corrupt blocks
	unsigned int total_width;

	/// current operational status of the coder
	LHNETCODEC_enum status;
}
LHNETCODEC_BlockCoderFixedPoint_Decode_State;

/// Initializes the state for a new encoding session.
LHNETCODEC_enum LHNETCODEC_BlockCoderFixedPoint_Decode_Init(LHNETCODEC_BlockCoderFixedPoint_Decode_State* state, const unsigned char* coded_start, size_t coded_length, LHNETCODEC_enum version_hint, unsigned int read_versioned_magic_number)
{
	const unsigned char* coded = coded_start;
	size_t maxbytes = coded_length;
	size_t cursor = 0;
	memset(state, 0, sizeof(*state));
	// default to this status until everything looks ok
	state->status = LHNETCODEC_STATUS_DECODE_CORRUPT;
	state->coded_version = version_hint;
	if (read_versioned_magic_number)
	{
		if (cursor + 4 <= maxbytes
			&& !memcmp(coded + cursor, "LHF0", 4))
		{
			state->coded_version = LHNETCODEC_BLOCKCODERFIXEDPOINT_VERSION_0;
			cursor += 4;
			state->coded_store_magic_number = 1;
		}
	}
	switch (state->coded_version)
	{
	case LHNETCODEC_BLOCKCODERFIXEDPOINT_VERSION_0:
	{
		unsigned int total_width = 0;
		unsigned int stream = 0;
		unsigned int flags = cursor < maxbytes ? coded[cursor] : 0;
		cursor++;
		state->contextlength = flags & 0x1F;
		if (state->contextlength > LHNETCODEC_MAXCONTEXTLENGTH)
		{
			state->status = LHNETCODEC_STATUS_DECODE_CORRUPT;
			return state->status;
		}
		state->context_mask = (1u << state->contextlength) - 1;
		state->total_width = total_width;
		state->coded_header_length = cursor;
		// when starting to decode the stream we need to read enough bytes to be
		// able to make decisions
		stream = 0;
		stream = (stream << 8) & 0xFFFF00; if (cursor < maxbytes) stream |= coded[cursor] & 0xFF; cursor++;
		stream = (stream << 8) & 0xFFFF00; if (cursor < maxbytes) stream |= coded[cursor] & 0xFF; cursor++;
		stream = (stream << 8) & 0xFFFF00; if (cursor < maxbytes) stream |= coded[cursor] & 0xFF; cursor++;
		state->coded_start = coded;
		state->coded_maxbytes = maxbytes;
		state->coded_cursor = cursor;
		state->range_low = 0;
		state->range_high = 0xFFFFFF;
		state->range_stream = stream;
		state->status = LHNETCODEC_STATUS_OK;
	}
	break;
	}
	return state->status;
}

/// attempts to decode one number as a set of bits
uint64_t LHNETCODEC_BlockCoderFixedPoint_Decode_Number(LHNETCODEC_BlockCoderFixedPoint_Decode_State* state, unsigned int firstcontext, unsigned int numcontexts)
{
	uint64_t number = 0;
	int position;
	unsigned int low = state->range_low;
	unsigned int high = state->range_high;
	unsigned int stream = state->range_stream;
	unsigned char* contexts = state->contexts;
	unsigned int contextlength = state->contextlength;
	unsigned int contextlength1 = contextlength;
	unsigned int contextsize2 = 1 << contextlength1;
	unsigned int contextmask = state->context_mask;
	size_t cursor = state->coded_cursor;
	const unsigned char* coded = state->coded_start;
	unsigned int* histogram = state->histogram;
	size_t maxbytes = state->coded_maxbytes;
	unsigned int total_width = state->total_width;
	if (state->status != LHNETCODEC_STATUS_OK)
		return 0;
	if (state->total_width == 0)
	{
		state->status = LHNETCODEC_STATUS_DECODE_EOF;
		return 0;
	}
	// if the histogram is exhausted partway through a number, the stream is
	// corrupt
	if (total_width < numcontexts)
	{
		state->status = LHNETCODEC_STATUS_DECODE_CORRUPT;
		return 0;
	}
	for (position = numcontexts - 1; position >= 0; position--)
	{
		unsigned int c = firstcontext + position;
		unsigned int context = contexts[c];
		unsigned int splitfixedpoint = state->probability_splits[(c << contextlength) + context];
		if (!state->histogramcontextsused[(c << contextlength) + context])
		{
			// first time we've seen this context, so we need to code the
			// probabilities
			state->histogramcontextsused[(position << contextlength) + context] = 1;
			DECODE_NUMBER(splitfixedpoint, 1 << LHNETCODEC_BLOCKCODERFIXEDPOINT_SPLITBITS);
			state->probability_splits[(position << contextlength) + context] = splitfixedpoint;
		}
		context *= 2;
		unsigned int split = low + (unsigned int)(((uint64_t)(high - low) * splitfixedpoint) >> LHNETCODEC_BLOCKCODERFIXEDPOINT_SPLITBITS);
		unsigned int bit = stream >= split ? 1 : 0;
		if (bit)
		{
			low = split;
			context++;
			number |= 1ull << position;
		}
		else
			high = split - 1;
		context &= contextmask;
		contexts[c] = context;
		DECODE_CHECKRANGE;
	}
	state->range_low = low;
	state->range_high = high;
	state->range_stream = stream;
	state->coded_cursor = cursor;
	state->total_width = total_width;
	return state->status;
}

#endif
