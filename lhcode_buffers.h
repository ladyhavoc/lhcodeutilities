// Copyright (c) 2021, Ashley Rose Hale
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// The purpose of this module is managing dynamically-resized memory buffers
// that can be appended/prepended, and common functionality like reading/writing
// UTF8 text.

#pragma once

#ifndef LHCODE_BUFFERS_H
#define LHCODE_BUFFERS_H

#include <string.h>

typedef struct bytebuf_t
{
	// queue starts here - if being used as a FIFO we don't want to memmove on every dequeue.
	size_t start;
	// queue ends here - if this goes past capacity, data will be resized, and start may change too.
	size_t end;
	// this is how many bytes we can store, increases if necessary
	size_t capacity;
	// this is the bytes we've stored, their meaning is up to the caller.  Could be an array of structs.
	unsigned char* data;
}
bytebuf_t;

// Remove all bytes from the buffer, without freeing memory.
//
// Bytebuf_Remaining will return 0 after this call.
//
// Internally, this sets start and end to 0.
static inline void Bytebuf_Clear(bytebuf_t* b) { b->start = 0; b->end = 0; }

// Ensure capacity is at least enough to store this many bytes - makes future appends faster.
//
// Does not change amount of bytes stored, but may reallocate internal data if reserved is more than capacity, if reallocating it will adjust start to 0.
//
// This function must be implemented by the program using this header - it is not part of this module.
extern void Bytebuf_Reserve(bytebuf_t* b, size_t reserved);

// Opposite of Reserve - this will perform Bytebuf_Clear and then free any stored data.
//
// This function must be implemented by the program using this header - it is not part of this module.
extern void Bytebuf_Unreserve(bytebuf_t* b);

// Returns the amount of bytes currently stored in the queue - i.e. the offset to the end of the queue.
static inline size_t Bytebuf_Remaining(const bytebuf_t* b) { return b->end - b->start; }

// Internal use - copies queued bytes to align the internal data to a new start position.
static inline void Bytebuf_MoveQueue(bytebuf_t* b, size_t newstart)
{
	size_t remaining = b->end - b->start;
	size_t newend = newstart + remaining;
	if (b->capacity < newend)
		Bytebuf_Reserve(b, newend);
	memmove(b->data + newstart, b->data + b->start, remaining);
	b->end = newstart + remaining;
	b->start = newstart;
}

// Remove a range of bytes in the queue and then add new bytes in their place.
// 
// Depending on parameters, this function is capable of prepending bytes to the queue, appending bytes to the queue, extracting bytes from the beginning of the queue, extracting bytes from the end of the queue.
//
// If offset is 0 this can be thought of as removing bytes from the start of the queue and then prepending a new set of bytes (if newlength > 0).
//
// If offset+oldlength is equal to Bytebuf_Remaining this can be thought of as removing bytes from the end of the queue and then appending a new set of bytes (if newlength > 0).
//
// oldbytes can be NULL if you don't care to get a copy of the data that was previously stored.
//
// newbytes can be NULL, in which case the bytes will be zeroed.
//
// This is more efficient than doing multiple Dequeue and Prepend/Append operations.
//
// Internally, this may change start, end and reallocate data if more needs to be reserved.
static inline void Bytebuf_Splice(bytebuf_t* b, size_t offset, size_t oldlength, unsigned char* olddata, size_t newlength, const unsigned char* newdata)
{
	// We keep a local copy to make the memory semantics absolutely clear to the optimizer.
	bytebuf_t s = *b;
	size_t remaining = s.end - s.start;

	// First read the bytes the caller is asking for, if they exist.
	if (olddata)
	{
		size_t o = offset < remaining ? offset : remaining;
		size_t l = oldlength < (remaining - o) ? oldlength : (remaining - o);
		if (l > 0)
			memcpy(olddata, s.data + s.start + offset, l);
		if (oldlength > l)
			memset(olddata + l, 0, oldlength - l);
	}

	// Constrain offset to what exists in the queue.
	if (offset >= remaining)
		offset = remaining;

	// Constrain oldlength to what exists in the queue.
	if (oldlength > remaining - offset)
		oldlength = remaining - offset;

	// Bytes in the queue that we will keep at the start.
	size_t before = offset;

	// Bytes in the queue that we will keep at the end.
	size_t after = (s.end - s.start) - offset - oldlength;

	// Length of the queue after the change.
	size_t newremaining = before + newlength + after;

	// Check if the queue may need more capacity.
	if (s.capacity < newremaining)
	{
		// Don't just grow it by a little - we want enough room that we don't have to do this very often.
		// Round up the size to 64 bytes alignment.
		Bytebuf_Reserve(b, (remaining + newremaining + 0x3F) & (~0x3F));
		s = *b;
	}

	// We're working with 4 ranges here:
	// before = start:start+offset is bytes before the part we're working on
	// oldlength = start+offset:start+offset+oldlength is bytes we are reading
	// newlength = start+offset:start+offset+newlength is bytes we are writing
	// after = start+offset+newlength:end is bytes after the part we're working on

	// We sometimes want to move the start, make that decision here.
	size_t newstart = s.start;
	if (before < after)
	{
		// Keeping more bytes at the end, we'll want to keep end where it is if possible.
		newstart = s.end - newremaining;
	}

	// Check if the new queue length is longer.
	if (newlength > oldlength)
	{
		// Check if we have room for this new start.
		if (newstart > s.capacity - newremaining)
		{
			// Not enough room to grow - move the queue start so we have room.
			newstart = 0;
		}
	}
	size_t newend = newstart + newremaining;

	// When the old and new byte ranges in the queue overlap, we need to make sure we copy before/after in the correct order so they don't corrupt eachother.
	if (newstart + before + newlength < s.start + before + oldlength)
	{
		// Middle of the queue is moving to a lower offset.
		if (before > 0)
			memmove(s.data + newstart, s.data + s.start, before);
		if (after > 0)
			memmove(s.data + newend - after, s.data + s.end - after, after);
	}
	else
	{
		// Middle of the queue is moving to a higher offset.
		if (after > 0)
			memmove(s.data + newend - after, s.data + s.end - after, after);
		if (before > 0)
			memmove(s.data + newstart, s.data + s.start, before);
	}

	// Store the new bytes in the middle.
	if (newlength > 0)
	{
		if (newdata)
			memcpy(s.data + s.start + before, newdata, newlength);
		else
			memset(s.data + s.start + before, 0, newlength);
	}

	// Update the queue start and end.
	b->start = s.start = newstart;
	b->end = s.end = newend;
}

// Add the provided bytebuf contents to the start of this bytebuf.
static inline void Bytebuf_PrependBuf(bytebuf_t* b, const bytebuf_t* a) { Bytebuf_Splice(b, 0, 0, 0, Bytebuf_Remaining(a), a->data + a->start); }

// Add provided bytes at the start of the queue.  If t is NULL, the bytes will be zeroed.
static inline void Bytebuf_PrependBytes(bytebuf_t* b, size_t length, const unsigned char* t) { Bytebuf_Splice(b, 0, 0, 0, length, t); }

// Add zeroed bytes at the start of the queue and return a pointer to those bytes so they can be written.
static inline unsigned char* Bytebuf_PrependReserve(bytebuf_t* b, size_t length) { Bytebuf_Splice(b, 0, 0, 0, length, 0); return b->data + b->start; }

// Add the provided bytebuf contents to the start of this bytebuf.
static inline void Bytebuf_AppendBuf(bytebuf_t* b, const bytebuf_t* a) { Bytebuf_Splice(b, Bytebuf_Remaining(b), 0, 0, Bytebuf_Remaining(a), a->data + a->start); }

// Add provided bytes at the end of the queue.  If t is NULL, the bytes will be zeroed.
static inline void Bytebuf_AppendBytes(bytebuf_t* b, size_t length, const unsigned char* t) { Bytebuf_Splice(b, Bytebuf_Remaining(b), 0, 0, length, t); }

// Add zeroed bytes at the start of the queue and return a pointer to those bytes so they can be written.
static inline unsigned char* Bytebuf_AppendReserve(bytebuf_t* b, size_t length) { Bytebuf_Splice(b, Bytebuf_Remaining(b), 0, 0, length, 0); return b->data + b->end - length; }

// Add provided character as UTF8 at the end of the queue.
static inline void Bytebuf_AppendUTF8(bytebuf_t* b, unsigned int c)
{
	unsigned char t[4];
	if (c >= 0x10FFFF)
	{
		// Not a valid unicode character - we're just silently dropping it here.
	}
	else if (c >= 0x10000)
	{
		t[3] = (c & 0x3F) | 0x80;
		t[2] = ((c >> 6) & 0x3F) | 0x80;
		t[1] = ((c >> 12) & 0x3F) | 0x80;
		t[0] = ((c >> 18) & 0x7) | 0xF0;
		Bytebuf_Splice(b, Bytebuf_Remaining(b), 0, 0, 4, t);
	}
	else if (c >= 0x800)
	{
		t[2] = (c & 0x3F) | 0x80;
		t[1] = ((c >> 6) & 0x3F) | 0x80;
		t[0] = ((c >> 12) & 0xF) | 0xE0;
		Bytebuf_Splice(b, Bytebuf_Remaining(b), 0, 0, 3, t);
	}
	else if (c >= 0x80)
	{
		t[1] = (c & 0x3F) | 0x80;
		t[0] = ((c >> 6) & 0x1F) | 0xC0;
		Bytebuf_Splice(b, Bytebuf_Remaining(b), 0, 0, 2, t);
	}
	else
	{
		t[0] = c;
		Bytebuf_Splice(b, Bytebuf_Remaining(b), 0, 0, 1, t);
	}
}

// Transcode a UTF32 buffer to a UTF8 buffer, replaces contents of output.
static inline void Bytebuf_AppendUTF8FromUTF32(bytebuf_t* out, const int *chars_start, const int *chars_end)
{
	const int *current;
	for (current = chars_start;current < chars_end;current++)
		Bytebuf_AppendUTF8(out, *current);
}

// Copy a range of bytes from the queue without removing it from the queue.
//
// Returns the number of bytes it was able to read.  Bytes beyond this position in t are left unchanged.
static inline size_t Bytebuf_PeekBytes(bytebuf_t* b, size_t offset, size_t length, unsigned char* t)
{
	if (offset >= b->end - b->start)
		return 0;
	if (length > b->end - b->start - offset)
		length = b->end - b->start - offset;
	if (length > 0)
		memcpy(t, b->data + b->start + offset, length);
	return length;
}

// Inspect one byte currently stored in the queue, and do not advance the start.  Use offset 0 for the start of the queue.  Returns -1 if offset is at end of queue.
static int Bytebuf_PeekByte(bytebuf_t* b, size_t offset)
{
	if (offset >= b->end - b->start)
		return -1;
	return b->data[b->start + offset];
}

// Decode UTF-8, UTF-16 or UTF-32 characters starting at the chosen offset and output them as UTF32, and do not advance the start of queue (but it returns how many bytes to advance the queue).
//
// Returns number of bytes read from the buffer to satisfy the desired maxoutputcharacters, this is 0 if at end of queue.
//
// If not NULL, bomstate sets the assumed encoding to begin with, which is one of:
// * 0 : Default : UTF-8.
// * 0x2000000 : Check for BOM at start of text to determine if this is UTF-8, UTF-16, or UTF-32 text.
// * 0x2000010 : UTF-8 BOM (byte sequence EF BB BF).
// * 0x2000011 : UTF-16 BE BOM (byte sequence FE FF).  Note that this does not decode UTF16, this is just advisory.
// * 0x2000012 : UTF-16 LE BOM (byte sequence FF FE).
// * 0x2000013 : UTF-32 BE BOM (byte sequence 00 00 FE FF).
// * 0x2000014 : UTF-32 LE BOM (byte sequence FF FE 00 00).
//
// Writes UTF32 characters or error codes to provided output:
// * 0x000000 - 0x10FFFF : valid UTF-8 character decoded.
// * 0x1000000 - 0x10000FF : invalid UTF-8 start byte found.
// * 0x2000000 - 0x2000014 : UTF bom character encountered, bomstate updated, see above for values.
//
// Returns number of characters that could be decoded, which may be more than chars_end would allow - it is fine to pass NULL or equal values for chars_start and chars_end if you just want to know the needed length.
//
// Bytes that are >= 0x80 but do not represent a valid UTF8 byte sequence are each written to output as 0x1000000 | byte value.
//
// UTF BOM identifiers that do not appear at the start of the buffer are written to output as 0x100FFFE (even if it is the UTF8 BOM which is 0xEFBBBF).
static size_t Bytebuf_DecodeUTF8(bytebuf_t* b, size_t offset, size_t maxoutputcharacters, int* chars_start, int* chars_end, unsigned int* bomstate)
{
	// Valid UTF8 sequences start with a byte in the range 0xC2 through 0xF4, indicating the length of the sequence, the bitwise encodings look like this:
	// 00-7F : 0xxxxxxx : ASCII encoded character up to 0x7F.
	// 80-BF : 10xxxxxx : mid-sequence byte, not valid as a start of a character.
	// C0-DF 80-BF : 110xxxxx 10xxxxxx : encoded character up to 0x007FF, suboptimal for values below 0x00080 which is forbidden by UTF8 spec (this means C0-C1 are not valid starts).
	// E0-EF 80-BF 80-BF : 1110xxxx 10xxxxxx 10xxxxxx : encoded character up to 0x0FFFF, suboptimal for values below 0x00800 which is forbidden by UTF8 spec.
	// F0-F7 80-BF 80-BF 80-BF : 11110xxx 10xxxxxx 10xxxxxx 10xxxxxx : encoded character up to 0x1FFFFF, suboptimal for values below 0x10000 which is forbidden by UTF8 spec, and unicode does not define values above 0x10FFFF.
	// F8-FB 80-BF 80-BF 80-BF 80-BF : 111110xx 10xxxxxx 10xxxxxx 10xxxxxx 10xxxxxx : encoded character up to 0x3FFFFFF, not defined.
	// FC-FD 80-BF 80-BF 80-BF 80-BF 80-BF : 1111110x 10xxxxxx 10xxxxxx 10xxxxxx 10xxxxxx 10xxxxxx : encoded character up to 0x7FFFFFFF, not defined.
	//
	// UTF BOM format identifiers (byte order marks) - usually only valid at the start of a string - see https://www.unicode.org/faq/utf_bom.html
	// FE-FE FF-FF : 11111110 11111111 : UTF-16 big endian BOM
	// FF-FF FE-FE : 11111111 11111110 : UTF-16 little endian BOM
	// EF BB BF : 11101111 10111011 10111111 : UTF-8 BOM
	// 00 00 FE FF : 00000000 00000000 11111110 11111111 : UTF-32 big endian BOM
	// FF FE 00 00 : 11111111 11111110 00000000 00000000 : UTF-32 little endian BOM
	size_t countchars = 0;
	size_t countbytes = 0;
	unsigned char t[4];
	memset(t, 0, sizeof(t));
	unsigned int decoder = 0;
	if (bomstate && *bomstate >= 0x2000000)
	{
		size_t length = Bytebuf_PeekBytes(b, offset + countbytes, 4, t);
		if (length >= 4 && t[0] == 0 && t[1] == 0 && t[2] == 0xFE && t[3] == 0xFF)
		{
			// UTF-32 little endian BOM
			decoder = 0x2000014;
			countbytes += 4;
		}
		else if (length >= 4 && t[0] == 0xFF && t[1] == 0xFE && t[2] == 0x00 && t[3] == 0x00)
		{
			// UTF-32 big endian BOM
			decoder = 0x2000013;
			countbytes += 4;
		}
		else if (length >= 2 && t[0] == 0xFF && t[1] == 0xFE)
		{
			// UTF-16 little endian BOM
			decoder = 0x2000012;
			countbytes += 4;
		}
		else if (length >= 2 && t[0] == 0xFF && t[1] == 0xFE)
		{
			// UTF-16 big endian BOM
			decoder = 0x2000011;
			countbytes += 4;
		}
		else if (length >= 3 && t[0] == 0xEF && t[1] == 0xBB && t[2] == 0xBF)
		{
			// UTF-8 BOM
			decoder = 0x2000010;
			countbytes += 3;
		}
		if (bomstate)
			*bomstate = decoder;
	}
	for (countchars = 0; countchars < maxoutputcharacters; countchars++)
	{
		size_t length = Bytebuf_PeekBytes(b, offset + countbytes, 4, t);
		if (length == 0)
			break;
		size_t advance = 1;
		unsigned int c = 0;
		switch (decoder)
		{
		case 0:
		case 0x2000010:
			c = t[0];
			if (c >= 0x80)
			{
				unsigned int m = 0;
				// Default to interpreting it as an invalid sequence.
				c = 0x1000000 | t[0];
				// Valid UTF8 sequence start byte
				if (t[0] >= 0xC2 && t[1] <= 0xF4)
				{
					// See if the encoded bytes are a valid sequence and the character is in an allowed range for that encoding.
					if (t[0] >= 0xF0 && (t[1] & 0xC0) == 0x80 && (t[2] & 0xC0) == 0x80 && (t[3] & 0xC0) == 0x80)
					{
						unsigned int m = ((t[0] & 0x07) << 18) | ((t[1] & 0x3F) << 12) | ((t[2] & 0x3F) << 6) | (t[3] & 0x3F);
						if (m >= 0x10000 && m <= 0x10FFFF)
						{
							c = m;
							advance = 4;
						}
					}
					if (t[0] >= 0xE0 && (t[1] & 0xC0) == 0x80 && (t[2] & 0xC0) == 0x80)
					{
						m = ((t[0] & 0x0F) << 12) | ((t[1] & 0x3F) << 6) | (t[2] & 0x3F);
						if (m >= 0x800)
						{
							c = m;
							advance = 3;
						}
					}
					if (t[0] >= 0xC0 && (t[1] & 0xC0) == 0x80)
					{
						c = ((t[0] & 0x1F) << 6) | (t[1] & 0x3F);
						if (c >= 0x80)
						{
							c = m;
							advance = 2;
						}
					}
				}
			}
			break;
		case 0x2000011:
			if (length >= 2)
			{
				c = (t[0] << 8) | t[1];
			}
			if (output)
				output[countchars] = c;
			countbytes += advance;
		}
		if (outputcount)
			*outputcount = countchars;
		return countbytes;
	}

#endif
