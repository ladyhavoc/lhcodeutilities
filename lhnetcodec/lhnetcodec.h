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
	LHNETCODEC_STATUS_OK = 0x1001,
	LHNETCODEC_STATUS_ENCODE_FULL = 0x2002,
	LHNETCODEC_STATUS_DECODE_EOF = 0x3002,
	LHNETCODEC_STATUS_DECODE_CORRUPT = 0x3003,

	// the first version of the BlockCoder codec
	LHNETCODEC_BLOCKCODER_VERSION_0 = 0x4000,

	/// maximum number of bits of context supported - a bit is predicted based
	/// on the same bit position in this many previous numbers
	LHNETCODEC_MAXCONTEXTBITS = 8,
	LHNETCODEC_MAXCONTEXTSIZE = (1u << LHNETCODEC_MAXCONTEXTBITS),
	LHNETCODEC_MAXHISTOGRAMSIZE = LHNETCODEC_MAXCONTEXTSIZE * 2,

	/// the context size is configurable, higher values allow more specialized
	/// adaptation but may result in lower compression ratios on short streams
	/// in the adaptive coder (in contrast to the block coder).
	LHNETCODEC_ADAPTIVECODER_DEFAULTCONTEXTBITS = 8,
	/// the initial value of the histogram pretends that all contexts had this
	/// many occurrences for each symbol.
	LHNETCODEC_ADAPTIVECODER_HISTOGRAM_INITIAL_VALUE = 1,

	/// the context size is configurable, higher values allow higher compression
	/// ratios in the block coder, but significantly increase the header size in
	/// the block coder, which may be counter productive on small blocks.
	LHNETCODEC_BLOCKCODER_DEFAULTCONTEXTBITS = 6,
	/// the block coder needs to replay the values being encoded each time it
	/// flushes, so this is the maximum queue it allows for one block
	LHNETCODEC_BLOCKCODER_REPLAY_NUMBERS = 65536,
}
LHNETCODEC_enum;

typedef struct LHNETCODEC_AdaptiveCoder_State
{
	/// counts of how many times a 0 or 1 has occurred so far after each
	/// possible context, this increments when encoding each bit, so it adapts
	/// to the numbers actually encoded in the stream, but has an initial phase
	/// where it is a poor fit (i.e. completely uncompressed).
	unsigned int histogram[LHNETCODEC_MAXHISTOGRAMSIZE];
	/// store the previous N bits for each of the 64 bits in a 64bit number
	/// being encoded/decoded, reading/writing smaller sizes will only advance
	/// the contexts for the number of bits specified by the caller.
	unsigned char contexts[64];
	/// context size in bits
	unsigned int context_bits;
	/// context_mask = (1u << context_bits) - 1;
	unsigned int context_mask;
	/// current coded data buffer which stores the resulting stream, can be NULL
	/// if merely being used for simulations.
	unsigned char* coded_start;
	/// size of the coded working buffer in bytes, when coded_cursor
	/// reaches this value, no more bits can be appended and an error status is
	/// returned.
	unsigned int coded_maxbytes;
	/// current position in the coded data buffer, in bytes.
	unsigned int coded_cursor;
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

LHNETCODEC_enum LHNETCODEC_AdaptiveCoder_Encode_Begin(LHNETCODEC_AdaptiveCoder_State* state, unsigned char *coded_start, unsigned int coded_maxbytes, unsigned int coded_cursor, unsigned int contextbits)
{
	unsigned int i;
	memset(state, 0, sizeof(*state));
	for (i = 0; i < LHNETCODEC_MAXHISTOGRAMSIZE; i++)
		state->histogram[i] = LHNETCODEC_ADAPTIVECODER_HISTOGRAM_INITIAL_VALUE;
	if (contextbits > LHNETCODEC_MAXCONTEXTBITS)
		contextbits = LHNETCODEC_MAXCONTEXTBITS;
	state->context_bits = contextbits;
	state->context_mask = (1u << contextbits) - 1;
	state->coded_start = coded_start;
	state->coded_maxbytes = coded_maxbytes;
	state->coded_cursor = coded_cursor;
	state->range_low = 0;
	state->range_high = 0xFFFFFF;
	state->status = LHNETCODEC_STATUS_OK;
	return state->status;
}

LHNETCODEC_enum LHNETCODEC_AdaptiveCoder_Decode_Begin(LHNETCODEC_AdaptiveCoder_State* state, unsigned char* coded_start, unsigned int coded_maxbytes, unsigned int coded_cursor, unsigned int contextbits)
{
	unsigned int i;
	memset(state, 0, sizeof(*state));
	for (i = 0; i < LHNETCODEC_MAXHISTOGRAMSIZE; i++)
		state->histogram[i] = LHNETCODEC_ADAPTIVECODER_HISTOGRAM_INITIAL_VALUE;
	if (contextbits > LHNETCODEC_MAXCONTEXTBITS)
		contextbits = LHNETCODEC_MAXCONTEXTBITS;
	state->context_bits = contextbits;
	state->context_mask = (1u << contextbits) - 1;
	unsigned int stream = 0;
	const unsigned char* coded = coded_start;
	unsigned int maxbytes = coded_maxbytes;
	unsigned int cursor = coded_cursor;
	state->status = LHNETCODEC_STATUS_DECODE_CORRUPT;
	// when starting to decode the stream we need to read enough bytes to be
	// able to make decisions
	if (coded)
	{
		stream = (stream << 8) & 0xFFFF00; if (cursor < maxbytes) stream |= coded[cursor] & 0xFF; cursor++;
		stream = (stream << 8) & 0xFFFF00; if (cursor < maxbytes) stream |= coded[cursor] & 0xFF; cursor++;
		stream = (stream << 8) & 0xFFFF00; if (cursor < maxbytes) stream |= coded[cursor] & 0xFF; cursor++;
		state->status = LHNETCODEC_STATUS_OK;
	}
	state->range_low = 0;
	state->range_high = 0xFFFFFF;
	state->range_stream = stream;
	state->coded_start = coded_start;
	state->coded_maxbytes = maxbytes;
	state->coded_cursor = cursor;
	return state->status;
}

LHNETCODEC_enum LHNETCODEC_AdaptiveCoder_Encode_Number(LHNETCODEC_AdaptiveCoder_State* state, uint64_t number, unsigned int width)
{
	int position;
	unsigned int low = state->range_low;
	unsigned int high = state->range_high;
	unsigned char* contexts = state->contexts;
	unsigned int contextmask = state->context_mask;
	unsigned int cursor = state->coded_cursor;
	unsigned char* coded = state->coded_start;
	unsigned int* histogram = state->histogram;
	unsigned int maxbytes = state->coded_maxbytes;
	if (state->status != LHNETCODEC_STATUS_OK)
		return state->status;
	for (position = width-1; position >= 0; position--)
	{
		unsigned int context = contexts[position] * 2;
		unsigned int probability_0 = histogram[context + 0];
		unsigned int probability_1 = histogram[context + 1];
		// reduce probabilities if they are getting so high as to make the
		// divide potentially unsafe (split must not be equal to low or high)
		while ((probability_0 | probability_1) > 0x10000)
		{
			if (probability_0 > 0)
			{
				if (probability_0 >= 2)
					probability_0 >>= 1;
				else
					probability_0 = 1;
			}
			if (probability_1 > 0)
			{
				if (probability_1 >= 2)
					probability_1 >>= 1;
				else
					probability_1 = 1;
			}
		}
		unsigned int split = low + (unsigned int)((uint64_t)(high - low) * probability_0 / (probability_0 + probability_1));
		unsigned int bit = (number >> position) & 1;
		if (bit)
		{
			low = split;
			context++;
		}
		else
			high = split - 1;
		histogram[context]++;
		context &= contextmask;
		contexts[position] = context;
		// when the low and high have the same significant byte, emit a byte
		// (this is the difference between range coders which emit bytes and
		//  arithmetic coders which emit individual bits)
		while (((low ^ high) & 0xFF0000) == 0)
		{
			if (cursor >= maxbytes)
			{
				state->status = LHNETCODEC_STATUS_ENCODE_FULL;
				break;
			}
			if (coded)
				coded[cursor] = (low >> 16) & 0xFF;
			cursor++;
			low = (low << 8) & 0xFFFFFF;
			high = (high << 8) & 0xFFFFFF;
			// This is probably impossible, but just in case...
			if (low == high)
			{
				low = 0;
				high = 0xFFFFFF;
			}
		}
	}
	state->range_low = low;
	state->range_high = high;
	state->coded_cursor = cursor;
	return state->status;
}

LHNETCODEC_enum LHNETCODEC_AdaptiveCoder_Encode_Flush(LHNETCODEC_AdaptiveCoder_State* state)
{
	unsigned int low = state->range_low;
	unsigned int cursor = state->coded_cursor;
	unsigned char* coded = state->coded_start;
	unsigned int maxbytes = state->coded_maxbytes;
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
		if (coded)
			coded[cursor] = (low >> 16) & 0xFF;
		cursor++;
		low = (low << 8) & 0xFFFFFF;
	}
	// we can trim trailing zero bytes here because the decoder will assume any
	// missing bytes are zero
	if (coded)
		while (cursor > 0 && coded[cursor - 1] == 0)
			cursor--;
	state->coded_cursor = cursor;
	return state->status;
}

uint64_t LHNETCODEC_AdaptiveCoder_Decode_Number(LHNETCODEC_AdaptiveCoder_State* state, unsigned int width)
{
	uint64_t number = 0;
	int position;
	unsigned int low = state->range_low;
	unsigned int high = state->range_high;
	unsigned int stream = state->range_stream;
	unsigned char* contexts = state->contexts;
	unsigned int contextmask = state->context_mask;
	unsigned int cursor = state->coded_cursor;
	unsigned char* coded = state->coded_start;
	unsigned int* histogram = state->histogram;
	unsigned int maxbytes = state->coded_maxbytes;
	if (state->status != LHNETCODEC_STATUS_OK || !state->coded_start)
		return 0;
	for (position = width - 1; position >= 0; position--)
	{
		unsigned int context = contexts[position] * 2;
		unsigned int probability_0 = histogram[context + 0];
		unsigned int probability_1 = histogram[context + 1];
		// reduce probabilities if they are getting so high as to make the
		// divide potentially unsafe (split must not be equal to low or high)
		while ((probability_0 | probability_1) > 0x10000)
		{
			if (probability_0 > 0)
			{
				if (probability_0 >= 2)
					probability_0 >>= 1;
				else
					probability_0 = 1;
			}
			if (probability_1 > 0)
			{
				if (probability_1 >= 2)
					probability_1 >>= 1;
				else
					probability_1 = 1;
			}
		}
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
		histogram[context]++;
		context &= contextmask;
		contexts[position] = context;
		// when the low and high have the same significant byte, read a byte
		// (this is the difference between range coders which read bytes and
		//  arithmetic coders which read individual bits)
		while (((low ^ high) & 0xFF0000) == 0)
		{
			if (cursor >= maxbytes)
			{
				state->status = LHNETCODEC_STATUS_DECODE_CORRUPT;
				break;
			}
			stream = ((stream << 8) | (coded[cursor++] & 0xFF)) & 0xFFFFFF;
			low = (low << 8) & 0xFFFFFF;
			high = (high << 8) & 0xFFFFFF;
			// This is probably impossible, but just in case...
			if (low == high)
			{
				low = 0;
				high = 0xFFFFFF;
			}
		}
	}
	state->range_low = low;
	state->range_high = high;
	state->range_stream = stream;
	state->coded_cursor = cursor;
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
	/// store the previous N bits for each of the 64 bits in a 64bit number
	/// being encoded/decoded, reading/writing smaller sizes will only advance
	/// the contexts for the used bits.
	unsigned char contexts[64];
	/// context size in bits
	unsigned int context_bits;
	/// context_mask = (1u << context_bits) - 1;
	unsigned int context_mask;
	/// current coded data buffer which stores the resulting stream, can be NULL
	/// if merely being used for simulations.
	unsigned char* coded_start;
	/// size of the coded working buffer in bytes, when coded_cursor_byte
	/// reaches this value, no more bits can be appended and an error status is
	/// returned.
	unsigned int coded_maxbytes;
	/// current position in the coded data buffer, in bytes.
	unsigned int coded_cursor;
	// this indicates a rollback has occurred but the encode has not been run,
	// again yet, so coded_cursor is a valid number for estimating if appending
	// numbers, but the actual data stream needs to be regenerated on flush.
	unsigned int finish_needs_flush;
	/// which version of the format we are encoding
	unsigned int coded_version;
	/// indicates the caller wants the versioned magic number to be included
	unsigned int coded_store_magic_number;
	/// indicates the value of cursor at the end of the header
	/// (after versioned magic number, contextbits, histogram)
	unsigned int coded_header_length;
	/// current operating range in the range coder
	unsigned int range_low;
	/// current operating range in the range coder
	unsigned int range_high;
	// width of each number in bits
	uint8_t numberwidths[LHNETCODEC_BLOCKCODER_REPLAY_NUMBERS];
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
LHNETCODEC_BlockCoder_Encode_State;

/// Initializes the state for a new encoding session.
void LHNETCODEC_BlockCoder_Encode_Init(LHNETCODEC_BlockCoder_Encode_State* state, unsigned char *coded_start, unsigned int coded_maxbytes, LHNETCODEC_enum version, unsigned int include_versioned_magic_number, unsigned int contextbits)
{
	memset(state, 0, sizeof(*state));
	state->coded_start = coded_start;
	state->coded_maxbytes = coded_maxbytes;
	state->coded_cursor = 0;
	state->coded_version = version;
	// we use this just to guarantee that Flush is called by Finish if Number
	// does not trigger a Flush
	state->finish_needs_flush = 1;
	state->status = LHNETCODEC_STATUS_UNINITIALIZED;
	switch (version) {
	case LHNETCODEC_BLOCKCODER_VERSION_0:
		if (contextbits > LHNETCODEC_MAXCONTEXTBITS)
			contextbits = LHNETCODEC_MAXCONTEXTBITS;
		state->context_bits = contextbits;
		state->context_mask = (1u << contextbits) - 1;
		// just bump the cursor so we don't underestimate the bytes we'll
		// write in Flush.
		state->coded_store_magic_number = include_versioned_magic_number;
		if (state->coded_store_magic_number)
			state->coded_cursor += 4;
		// add a byte for contextbits
		state->coded_cursor++;
		// add bytes for histogram
		state->coded_cursor += 2 << contextbits;
		// this is only an estimate - Flush will determine the real length
		state->coded_header_length = state->coded_cursor;
		state->coded_version = LHNETCODEC_BLOCKCODER_VERSION_0;
		state->status = LHNETCODEC_STATUS_OK;
		break;
	// no default case because we want the compiler to warn if a new version is
	// not handled here
	}
}

/// attempts to encode the buffered numbers, this is not expected to fail when
/// called by AppendNumber which uses conservative estimates
LHNETCODEC_enum LHNETCODEC_BlockCoder_Encode_Flush(LHNETCODEC_BlockCoder_Encode_State* state)
{
	// reset all encoder state except the numbers buffer
	memset(state->histogram, 0, sizeof(state->histogram));
	memset(state->contexts, 0, sizeof(state->contexts));
	state->number_encodedwidth = 0;
	state->range_low = 0;
	state->range_high = 0xFFFFFF;
	state->coded_cursor = 0;
	state->finish_needs_flush = 0;
	// keep local copies of the variables so the compiler knows the values are
	// completely under our control here
	unsigned char* numberwidths = state->numberwidths;
	uint64_t* numbers = state->numbers;
	unsigned int low = state->range_low;
	unsigned int high = state->range_high;
	unsigned char* contexts = state->contexts;
	unsigned int contextmask = state->context_mask;
	unsigned int cursor = state->coded_cursor;
	unsigned char* coded = state->coded_start;
	unsigned int* histogram = state->histogram;
	unsigned int maxbytes = state->coded_maxbytes;
	// calculate the histogram before encoding with it
	unsigned int n;
	for (n = 0; n < state->number_count; n++)
	{
		int position;
		uint64_t number = numbers[n];
		unsigned int width = numberwidths[n];
		for (position = width - 1; position >= 0;position--)
		{
			unsigned int context = contexts[position] * 2 + ((number >> position) & 1);
			histogram[context]++;
			contexts[position] = context & contextmask;
		}
	}
	memset(state->contexts, 0, sizeof(state->contexts));
	// write the magic number if the caller wants to
	if (state->coded_store_magic_number)
	{
		switch (state->coded_version)
		{
		case LHNETCODEC_BLOCKCODER_VERSION_0:
			if (cursor + 5 > maxbytes)
			{
				state->status = LHNETCODEC_STATUS_ENCODE_FULL;
				return state->status;
			}
			if (coded)
			{
				memcpy(coded + cursor, "LHN0", 4);
				coded[cursor + 4] = state->context_bits & 0xFF;
			}
			cursor += 5;
			break;
		}
	}
	// now encode the histogram before the compressed stream
	unsigned int histogramsize = 2u << state->context_bits;
	for (n = 0;n < histogramsize;n++)
	{
		// encode the number as ULEB128
		unsigned int number = histogram[n];
		unsigned int b;
		do
		{
			b = (number & 0x7F);
			if (number > b)
				b |= 0x80;
			if (cursor >= maxbytes)
			{
				state->status = LHNETCODEC_STATUS_ENCODE_FULL;
				return state->status;
			}
			if (coded)
				coded[cursor] = b;
			cursor++;
			number >>= 7;
		}
		while (number);
	}
	state->coded_header_length = cursor;
	// now encode the numbers with the diminishing histogram
	for (n = 0; n < state->number_count; n++)
	{
		int position;
		uint64_t number = numbers[n];
		unsigned int width = numberwidths[n];
		state->number_encodedwidth += width;
		for (position = width - 1; position >= 0; position--)
		{
			unsigned int context = contexts[position] * 2;
			unsigned int probability_0 = histogram[context + 0];
			unsigned int probability_1 = histogram[context + 1];
			// reduce probabilities if they are getting so high as to make the
			// divide potentially unsafe (split must not be equal to low or high)
			while ((probability_0 | probability_1) > 0x10000)
			{
				if (probability_0 > 0)
				{
					if (probability_0 >= 2)
						probability_0 >>= 1;
					else
						probability_0 = 1;
				}
				if (probability_1 > 0)
				{
					if (probability_1 >= 2)
						probability_1 >>= 1;
					else
						probability_1 = 1;
				}
			}
			unsigned int split = low + (unsigned int)((uint64_t)(high - low) * probability_0 / (probability_0 + probability_1));
			unsigned int bit = (number >> position) & 1;
			if (bit)
			{
				low = split;
				context++;
			}
			else
				high = split - 1;
			histogram[context]--;
			context &= contextmask;
			contexts[position] = context;
			// when the low and high have the same significant byte, emit a byte
			// (this is the difference between range coders which emit bytes and
			//  arithmetic coders which emit individual bits)
			while (((low ^ high) & 0xFF0000) == 0)
			{
				if (cursor >= maxbytes)
				{
					state->status = LHNETCODEC_STATUS_ENCODE_FULL;
					break;
				}
				if (coded)
					coded[cursor] = (low >> 16) & 0xFF;
				cursor++;
				low = (low << 8) & 0xFFFFFF;
				high = (high << 8) & 0xFFFFFF;
				// This is probably impossible, but just in case...
				if (low == high)
				{
					low = 0;
					high = 0xFFFFFF;
				}
			}
		}
	}
	// now flush the remaining bytes until low becomes zero, this may write up
	// to 3 bytes
	while (low)
	{
		if (coded)
			coded[cursor] = (low >> 16) & 0xFF;
		cursor++;
		low = (low << 8) & 0xFFFFFF;
	}
	// we can trim trailing zero bytes here because the decoder will assume any
	// missing bytes are zero
	if (coded)
		while (cursor > state->coded_header_length && coded[cursor - 1] == 0)
			cursor--;
	// at this point, the entire block has been encoded, and histogram[] is
	// entirely zero (as we decrmeneted while encoding), encodedwidth is equal
	// to total_width
	state->coded_cursor = cursor;
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
	if (state->coded_cursor + fullness_threshold > state->coded_maxbytes)
	{
		state->status = LHNETCODEC_STATUS_ENCODE_FULL;
		return state->status;
	}
	// this causes Flush to be called more often as the block fills up,
	// which refines our estimated remaining space, but at high cost.
	unsigned int outstanding_bytes = finishing_bytes
		+ ((state->number_totalwidth - state->number_encodedwidth + 7) >> 3);
	if (state->coded_cursor + outstanding_bytes > state->coded_maxbytes)
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
LHNETCODEC_enum LHNETCODEC_BlockCoder_Encode_Number(LHNETCODEC_BlockCoder_Encode_State* state, uint64_t number, unsigned int width)
{
	if (state->status != LHNETCODEC_STATUS_OK)
		return state->status;
	if (state->number_count >= LHNETCODEC_BLOCKCODER_REPLAY_NUMBERS)
	{
		state->status = LHNETCODEC_STATUS_ENCODE_FULL;
		return state->status;
	}
	state->numberwidths[state->number_count] = width;
	state->numbers[state->number_count] = number;
	state->number_totalwidth += width;
	state->number_count++;
	return state->status;
}

typedef struct LHNETCODEC_BlockCoder_Encode_RollbackPosition
{
	unsigned int number_count;
	unsigned int cursor;
}
LHNETCODEC_BlockCoder_Encode_RollbackPosition;

/// gets a rollback position for undoing Number operations when a group of
/// several related numbers will only partially fit - for example a data
/// structure that you want to be included in its entirety or not at all.
LHNETCODEC_BlockCoder_Encode_RollbackPosition LHNETCODEC_BlockCoder_Encode_GetRollbackPosition(const LHNETCODEC_BlockCoder_Encode_State* state)
{
	LHNETCODEC_BlockCoder_Encode_RollbackPosition p;
	p.number_count = state->number_count;
	p.cursor = state->coded_cursor;
	return p;
}

LHNETCODEC_enum LHNETCODEC_BlockCoder_Encode_Rollback(LHNETCODEC_BlockCoder_Encode_State* state, LHNETCODEC_BlockCoder_Encode_RollbackPosition p)
{
	state->status = LHNETCODEC_STATUS_OK;
	while (state->number_count > p.number_count)
		state->number_totalwidth -= state->numberwidths[--state->number_count];
	state->coded_cursor = p.cursor;
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
		*coded_length = state->coded_cursor;
	return state->status;
}

typedef struct LHNETCODEC_BlockCoder_Decode_State
{
	/// counts of how many times a 0 or 1 occurs after each possible context,
	/// this decrements when decoding each bit, so the buffer is fully decoded
	/// when remaining_bits == 0 and remaining_bits is a sum of histogram[].
	unsigned int histogram[LHNETCODEC_MAXHISTOGRAMSIZE];
	/// sum of histogram[], used for detecting corrupt blocks
	unsigned int total_width;
	/// store the previous N bits for each of the 64 bits in a 64bit number
	/// being encoded/decoded, reading/writing smaller sizes will only advance
	/// the contexts for the used bits.
	unsigned char contexts[64];
	/// context size in bits
	unsigned int context_bits;
	/// context_mask = (1u << context_bits) - 1;
	unsigned int context_mask;
	/// current coded data buffer which stores the resulting stream, can be NULL
	/// if merely being used for simulations.
	const unsigned char* coded_start;
	/// size of the coded working buffer in bytes, when coded_cursor
	/// reaches this value, no more bits can be appended and an error status is
	/// returned.
	unsigned int coded_maxbytes;
	/// current position in the coded data buffer, in bytes.
	unsigned int coded_cursor;
	/// which version of the format we are decoding
	unsigned int coded_version;
	/// indicates the block has a magic number included
	unsigned int coded_store_magic_number;
	/// indicates the value of cursor at the end of the header
	/// (after versioned magic number, contextbits, histogram)
	unsigned int coded_header_length;
	/// current operating range in the range coder
	unsigned int range_low;
	/// current operating range in the range coder
	unsigned int range_high;
	/// (decode only) fractional value being streamed in
	unsigned int range_stream;

	/// current operational status of the coder
	LHNETCODEC_enum status;
}
LHNETCODEC_BlockCoder_Decode_State;

/// Initializes the state for a new encoding session.
LHNETCODEC_enum LHNETCODEC_BlockCoder_Decode_Init(LHNETCODEC_BlockCoder_Decode_State* state, const unsigned char* coded_start, unsigned int coded_length, LHNETCODEC_enum version_hint, unsigned int read_versioned_magic_number)
{
	const unsigned char* coded = coded_start;
	unsigned int maxbytes = coded_length;
	unsigned int cursor = 0;
	memset(state, 0, sizeof(*state));
	// default to this status until everything looks ok
	state->status = LHNETCODEC_STATUS_DECODE_CORRUPT;
	state->coded_version = version_hint;
	if (read_versioned_magic_number)
	{
		if (cursor + 5 <= maxbytes
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
			// read the histogram
			unsigned int n;
			unsigned int total_width = 0;
			unsigned int stream = 0;
			if (cursor >= maxbytes)
			{
				state->status = LHNETCODEC_STATUS_DECODE_CORRUPT;
				return state->status;
			}
			state->context_bits = coded[cursor++];
			if (state->context_bits > LHNETCODEC_MAXCONTEXTBITS)
			{
				state->status = LHNETCODEC_STATUS_DECODE_CORRUPT;
				return state->status;
			}
			state->context_mask = (1u << state->context_bits) - 1;
			unsigned int histogramsize = 2u << state->context_bits;
			for (n = 0; n < histogramsize; n++)
			{
				unsigned int number = 0;
				unsigned int b;
				unsigned int shift = 0;
				do
				{
					if (cursor >= maxbytes)
					{
						state->status = LHNETCODEC_STATUS_DECODE_CORRUPT;
						return state->status;
					}
					b = coded[cursor];
					cursor++;
					number |= ((b & 0x7F) << shift);
					shift += 7;
				} while (b & 0x80);
				// we don't really care if this is accurate if the data turns
				// out to be corrupt because we'll discover that later (the
				// histogram acts as a CRC after all)
				state->histogram[n] = number;
				total_width += number;
			}
			state->total_width = total_width;
			state->coded_header_length = cursor;
			// when starting to decode the stream we need to read enough bytes to be
			// able to make decisions
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
uint64_t LHNETCODEC_BlockCoder_Decode_Number(LHNETCODEC_BlockCoder_Decode_State* state, unsigned int width)
{
	uint64_t number = 0;
	int position;
	unsigned int low = state->range_low;
	unsigned int high = state->range_high;
	unsigned int stream = state->range_stream;
	unsigned char* contexts = state->contexts;
	unsigned int contextmask = state->context_mask;
	unsigned int cursor = state->coded_cursor;
	const unsigned char* coded = state->coded_start;
	unsigned int* histogram = state->histogram;
	unsigned int maxbytes = state->coded_maxbytes;
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
	if (total_width < width)
	{
		state->status = LHNETCODEC_STATUS_DECODE_CORRUPT;
		return 0;
	}
	for (position = width - 1; position >= 0; position--)
	{
		unsigned int context = contexts[position] * 2;
		unsigned int probability_0 = histogram[context];
		unsigned int probability_1 = histogram[context + 1];
		// reduce probabilities if they are getting so high as to make the
		// divide potentially unsafe (split must not be equal to low or high)
		while ((probability_0 | probability_1) > 0x10000)
		{
			if (probability_0 > 0)
			{
				if (probability_0 >= 2)
					probability_0 >>= 1;
				else
					probability_0 = 1;
			}
			if (probability_1 > 0)
			{
				if (probability_1 >= 2)
					probability_1 >>= 1;
				else
					probability_1 = 1;
			}
		}
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
		histogram[context]--;
		total_width--;
		context &= contextmask;
		contexts[position] = context;
		// when the low and high have the same significant byte, read a byte
		// (this is the difference between range coders which read bytes and
		//  arithmetic coders which read individual bits)
		while (((low ^ high) & 0xFF0000) == 0)
		{
			// if we run out of bytes, assume there are just trailing zeros, the
			// encoder trims trailing zero bytes
			unsigned char b = (cursor < maxbytes) ? (coded[cursor] & 0xFF) : 0;
			cursor++;
			stream = ((stream << 8) | b) & 0xFFFFFF;
			low = (low << 8) & 0xFFFFFF;
			high = (high << 8) & 0xFFFFFF;
			// This is probably impossible, but just in case...
			if (low == high)
			{
				low = 0;
				high = 0xFFFFFF;
			}
		}
	}
	state->range_low = low;
	state->range_high = high;
	state->range_stream = stream;
	state->coded_cursor = cursor;
	state->total_width = total_width;
	return state->status;
}

#endif
