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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define LHNETCODEC_ENABLE_ADAPTIVECODER 1
#define LHNETCODEC_ENABLE_BLOCKCODER 1

#include "lhnetcodec.h"

int findflag(int argc, char** argv, const char* flagname, int minargs, int maxargs, int *value, int defaultvalue, char const **strvalue, const char *defaultstringvalue)
{
	int i;
	if (value)
		*value = defaultvalue;
	if (strvalue)
		*strvalue = defaultstringvalue;
	for (i = 1; i < argc;i++)
	{
		if (!strcmp(argv[i], flagname))
		{
			int j;
			for (j = i + 1; j < argc; j++)
				if (argv[j][0] == '-' && argv[j][1] != '-' && argv[j][1] != 0)
					break;
			if (j - i < minargs || j - i > maxargs)
			{
				if (minargs == maxargs)
					fprintf(stderr, "flag %s takes %i args\n", flagname, minargs);
				else
					fprintf(stderr, "flag %s takes between %i and %i args\n", flagname, minargs, maxargs);
				return -1;
			}
			const char* first = (i < argc && argv[i + 1]) ? argv[i + 1] : "";
			int v = atoi(first);
			if (value && v != 0 || first[0] == '0')
				*value = v;
			if (strvalue)
				*strvalue = first;
			return i + 1;
		}
	}
	return 0;
}

int refillbuffer(FILE* file, unsigned char* buffer, size_t* bufferused, size_t buffercapacity)
{
	size_t used = *bufferused;
	size_t want = buffercapacity - used;
	if (want > 0x1000000)
		want = 0x1000000;
	size_t r = fread(buffer + used, 1, want, file);
	if (r <= want)
		*bufferused = used + r;
	return (int)r;
}

void printusage(int argc, char** argv)
{
	printf(
		"usage: %s command [-flags]\n"
		"about: a tool by LadyHavoc to test lhnetcodec.h\n"
		"flags:\n"
		"-v 1 - verbose messages\n"
		"-v 2 - verbose messages, show context histogram\n"
		"-contextlength # - set the number of context bits to encode with/decode with\n"
		"-pbase # - set the probability_base value for the AdaptiveCoder\n"
		"-pstep # - set the probability_step value for the AdaptiveCoder\n"
		"-plimit # - set the probability_limit value for the AdaptiveCoder\n"
		"-pshift # - set the probability_shift value for the AdaptiveCoder\n"
		"-infile abc - set the input file, use - for stdin, otherwise a test string is used\n"
		"-outfile abc - set the output file, use - for stdout, otherwise only stats are printed\n"
		"- encode_ac - encode file using AdaptiveCoder\n"
		"- decode_ac - decode file using AdaptiveCoder\n"
		"- encode_bc - encode file using BlockCoder\n"
		"- decode_bc - decode file using BlockCoder\n"
		"- encode_bcfp - encode file using BlockCoderFixedPoint\n"
		"- decode_bcfp - decode file using BlockCoderFixedPoint\n"
		, argv[0] ? argv[0] : "lhnetcodec_test_util"
	);
}

int main(int argc, char** argv)
{
	if (argc == 1)
	{
		printusage(argc, argv);
		return 0;
	}
	int verbose; if (findflag(argc, argv, "-v", 1, 1, &verbose, 2, NULL, NULL) < 0) return -1;
	int contextlength; if (findflag(argc, argv, "-contextlength", 1, 1, &contextlength, -1, NULL, NULL) < 0) return -1;
	int probability_base; if (findflag(argc, argv, "-pbase", 1, 1, &probability_base, LHNETCODEC_ADAPTIVECODER_DEFAULT_PROBABILITY_BASE, NULL, NULL) < 0) return -1;
	int probability_step; if (findflag(argc, argv, "-pstep", 1, 1, &probability_step, LHNETCODEC_ADAPTIVECODER_DEFAULT_PROBABILITY_STEP, NULL, NULL) < 0) return -1;
	int probability_limit; if (findflag(argc, argv, "-plimit", 1, 1, &probability_limit, LHNETCODEC_ADAPTIVECODER_DEFAULT_PROBABILITY_LIMIT, NULL, NULL) < 0) return -1;
	int probability_shift; if (findflag(argc, argv, "-pshift", 1, 1, &probability_shift, LHNETCODEC_ADAPTIVECODER_DEFAULT_PROBABILITY_SHIFT, NULL, NULL) < 0) return -1;
	int write_versioned_magic_number; if (findflag(argc, argv, "-use_magic", 1, 1, &write_versioned_magic_number, 1, NULL, NULL) < 0) return -1;
	const char* infilename; if (findflag(argc, argv, "-infile", 1, 1, NULL, 0, &infilename, "") < 0) return -1;
	const char* outfilename; if (findflag(argc, argv, "-outfile", 1, 1, NULL, 0, &outfilename, "") < 0) return -1;
	FILE* infile = infilename[0] != 0 ? (!strcmp(infilename, "-") ? stdin : fopen(infilename, "rb")) : NULL;
	FILE* outfile = outfilename[0] != 0 ? (!strcmp(outfilename, "-") ? stdout : fopen(outfilename, "wb")) : NULL;
	static unsigned char textbuffer[65536];
	static unsigned char codedbuffer[65536];
	const char* testtext = "TEST TEST The quick brown fox jumps over the lazy dog\nHello World!\nTEST TEST The quick brown fox jumps over the lazy dog\nHello World!\nTEST TEST The quick brown fox jumps over the lazy dog\nHello World!\nTEST TEST The quick brown fox jumps over the lazy dog\nHello World!\n";
	memcpy(textbuffer, testtext, strlen(testtext));
	size_t textbufferlength = strlen(testtext);
	size_t codedbufferlength = 0;
	if (contextlength == -1 && (!strcmp(argv[1], "encode_ac") || !strcmp(argv[1], "decode_ac")))
		contextlength = LHNETCODEC_ADAPTIVECODER_DEFAULT_CONTEXTLENGTH;
	if (contextlength == -1 && (!strcmp(argv[1], "encode_bc") || !strcmp(argv[1], "decode_bc")))
		contextlength = LHNETCODEC_BLOCKCODER_DEFAULT_CONTEXTLENGTH;
	if (contextlength == -1 && (!strcmp(argv[1], "encode_bcfp") || !strcmp(argv[1], "decode_bcfp")))
		contextlength = LHNETCODEC_BLOCKCODER_DEFAULT_CONTEXTLENGTH;
	if (infile && (!strcmp(argv[1], "encode_ac") || !strcmp(argv[1], "encode_bc") || !strcmp(argv[1], "encode_bcfp")))
	{
		textbufferlength = 0;
		for (;;)
		{
			int r = refillbuffer(infile, textbuffer, &textbufferlength, sizeof(textbuffer));
			if (r < 0)
			{
				fprintf(stderr, "fread returned %i", r);
				return -1;
			}
			if (r > 0)
				textbufferlength += r;
			if (r == 0 || errno == EOF || textbufferlength == sizeof(textbuffer))
				break;
		}
	}
	if (infile && (!strcmp(argv[1], "decode_ac") || !strcmp(argv[1], "decode_bc") || !strcmp(argv[1], "decode_bcfp")))
	{
		codedbufferlength = 0;
		for (;;)
		{
			int r = refillbuffer(infile, codedbuffer, &codedbufferlength, sizeof(codedbuffer));
			if (r < 0)
			{
				fprintf(stderr, "fread returned %i", r);
				return -1;
			}
			if (r > 0)
				codedbufferlength += r;
			if (r == 0 || errno == EOF || codedbufferlength == sizeof(codedbuffer))
				break;
		}
	}
	if (!strcmp(argv[1], "encode_ac"))
	{
		unsigned int p;
		static LHNETCODEC_AdaptiveCoder_State state;
		LHNETCODEC_AdaptiveCoder_Encode_Begin(
			&state, codedbuffer, sizeof(codedbuffer), contextlength,
			probability_base, probability_step, probability_limit,
			probability_shift);
		for (p = 0; p < textbufferlength; p++)
		{
			uint64_t number = textbuffer[p];
			LHNETCODEC_AdaptiveCoder_Encode_Number(&state, 0, 8, number);
		}
		LHNETCODEC_enum status = LHNETCODEC_AdaptiveCoder_Encode_Flush(&state);
		if (status != LHNETCODEC_STATUS_OK)
		{
			fprintf(stderr, "LHNETCODEC_AdaptiveCoder_Encode_Flush returned %i", (int)status);
			return -1;
		}
		if (verbose > 0)
		{
			fprintf(stderr,
				"AdaptiveCoder stats:\n"
				"bytes in: %u\n"
				"cursor: %u\n"
				"contextlength: %u\n"
				, (unsigned int)textbufferlength
				, (unsigned int)state.coded_cursor_trimmed
				, (unsigned int)state.contextlength
			);
			if (verbose > 1)
			{
				fprintf(stderr, "histogram:\n");
				char s[256];
				unsigned int context;
				unsigned int contextsize = 1u << state.contextlength;
				unsigned int position;
				const unsigned char* h = state.histogram;
				for (position = 0; position < LHNETCODEC_MAXCONTEXTRANGE; position++)
				{
					for (context = 0; context < contextsize; context++, h += 2)
					{
						if (h[0] + h[1] == 0)
							continue;
						unsigned int bit;
						unsigned int p = 0;
						s[p++] = '0';
						s[p++] = 'b';
						for (bit = contextsize >> 1; bit; bit >>= 1)
							s[p++] = '0' + ((context & bit) ? 1 : 0);
						s[p] = 0;
						fprintf(stderr, "%i:%s %i %i\n", position, s, h[0], h[1]);
					}
				}
			}
		}
	}
	else if (!strcmp(argv[1], "decode_ac"))
	{
		return 1;
	}
	else if (!strcmp(argv[1], "encode_bc"))
	{
		unsigned int p;
		LHNETCODEC_enum version = LHNETCODEC_BLOCKCODER_VERSION_0;
		static LHNETCODEC_BlockCoder_Encode_State state;
		LHNETCODEC_BlockCoder_Encode_Init(&state, codedbuffer, sizeof(codedbuffer), version, 1, contextlength);
		for (p = 0; p < textbufferlength; p++)
		{
			uint64_t number = textbuffer[p];
			LHNETCODEC_BlockCoder_Encode_Number(&state, 0, 8, number);
		}
		uint64_t coded_length = 0;
		LHNETCODEC_enum status = LHNETCODEC_BlockCoder_Encode_Finish(&state, NULL, &coded_length);
		if (status != LHNETCODEC_STATUS_OK)
		{
			fprintf(stderr, "LHNETCODEC_BlockCoder_Encode_Flush returned %i", (int)status);
			return -1;
		}
		if (verbose > 0)
		{
			static LHNETCODEC_BlockCoder_Decode_State dstate;
			LHNETCODEC_BlockCoder_Decode_Init(&dstate, codedbuffer, (unsigned int)coded_length, version, write_versioned_magic_number);
			for (p = 0; p < textbufferlength; p++)
				LHNETCODEC_BlockCoder_Decode_Number(&dstate, 0, 8);
			fprintf(stderr,
				"BlockCoder stats:\n"
				"bytes in: %u\n"
				"cursor: %u\n"
				"contextlength: %u\n"
				"store_magic_number: %u\n"
				"version: %u\n"
				"header length: %u\n"
				, (unsigned int)textbufferlength
				, (unsigned int)coded_length
				, (unsigned int)dstate.contextlength
				, (unsigned int)dstate.coded_store_magic_number
				, (unsigned int)dstate.coded_version
				, (unsigned int)dstate.coded_header_length
			);
			if (verbose > 1)
			{
				// we need to build the histogram again if we want to print it
				static LHNETCODEC_BlockCoder_Encode_State hstate;
				LHNETCODEC_BlockCoder_Encode_Init(&hstate, codedbuffer, sizeof(codedbuffer), version, 1, contextlength);
				for (p = 0; p < textbufferlength; p++)
				{
					uint64_t number = textbuffer[p];
					LHNETCODEC_BlockCoder_Encode_Number(&hstate, 0, 8, number);
				}
				LHNETCODEC_BlockCoder_GenerateHistogram(&hstate);
				fprintf(stderr, "histogram:\n");
				char s[256];
				unsigned int context;
				unsigned int contextsize = 1u << hstate.contextlength;
				unsigned int position;
				const unsigned int* h = hstate.histogram;
				for (position = 0; position < LHNETCODEC_MAXCONTEXTRANGE; position++)
				{
					for (context = 0; context < contextsize; context++, h += 2)
					{
						if (h[0] + h[1] == 0)
							continue;
						unsigned int bit;
						unsigned int p = 0;
						s[p++] = '0';
						s[p++] = 'b';
						for (bit = contextsize >> 1; bit; bit >>= 1)
							s[p++] = '0' + ((context & bit) ? 1 : 0);
						s[p] = 0;
						fprintf(stderr, "%i:%s %i %i\n", position, s, h[0], h[1]);
					}
				}
			}
		}
	}
	else if (!strcmp(argv[1], "decode_bc"))
	{
		return 1;
	}
	else if (!strcmp(argv[1], "encode_bcfp"))
	{
		unsigned int p;
		LHNETCODEC_enum version = LHNETCODEC_BLOCKCODERFIXEDPOINT_VERSION_0;
		static LHNETCODEC_BlockCoderFixedPoint_Encode_State state;
		LHNETCODEC_BlockCoderFixedPoint_Encode_Init(&state, codedbuffer, sizeof(codedbuffer), version, 1, contextlength);
		for (p = 0; p < textbufferlength; p++)
		{
			uint64_t number = textbuffer[p];
			LHNETCODEC_BlockCoderFixedPoint_Encode_Number(&state, 0, 8, number);
		}
		uint64_t coded_length = 0;
		LHNETCODEC_enum status = LHNETCODEC_BlockCoderFixedPoint_Encode_Finish(&state, NULL, &coded_length);
		if (status != LHNETCODEC_STATUS_OK)
		{
			fprintf(stderr, "LHNETCODEC_BlockCoderFixedPoint_Encode_Flush returned %i", (int)status);
			return -1;
		}
		if (verbose > 0)
		{
			static LHNETCODEC_BlockCoderFixedPoint_Decode_State dstate;
			LHNETCODEC_BlockCoderFixedPoint_Decode_Init(&dstate, codedbuffer, (unsigned int)coded_length, version, write_versioned_magic_number);
			for (p = 0; p < textbufferlength; p++)
				LHNETCODEC_BlockCoderFixedPoint_Decode_Number(&dstate, 0, 8);
			fprintf(stderr,
				"BlockCoderFixedPoint stats:\n"
				"bytes in: %u\n"
				"cursor: %u\n"
				"contextlength: %u\n"
				"store_magic_number: %u\n"
				"version: %u\n"
				"header length: %u\n"
				, (unsigned int)textbufferlength
				, (unsigned int)coded_length
				, (unsigned int)dstate.contextlength
				, (unsigned int)dstate.coded_store_magic_number
				, (unsigned int)dstate.coded_version
				, (unsigned int)dstate.coded_header_length
			);
			if (verbose > 1)
			{
				// we need to build the histogram again if we want to print it
				static LHNETCODEC_BlockCoderFixedPoint_Encode_State hstate;
				LHNETCODEC_BlockCoderFixedPoint_Encode_Init(&hstate, codedbuffer, sizeof(codedbuffer), version, 1, contextlength);
				for (p = 0; p < textbufferlength; p++)
				{
					uint64_t number = textbuffer[p];
					LHNETCODEC_BlockCoderFixedPoint_Encode_Number(&hstate, 0, 8, number);
				}
				LHNETCODEC_BlockCoderFixedPoint_GenerateHistogram(&hstate);
				fprintf(stderr, "histogram:\n");
				char s[256];
				unsigned int context;
				unsigned int contextsize = 1u << hstate.contextlength;
				unsigned int position;
				const unsigned int* h = hstate.histogram;
				for (position = 0; position < LHNETCODEC_MAXCONTEXTRANGE; position++)
				{
					for (context = 0; context < contextsize; context++, h += 2)
					{
						if (h[0] + h[1] == 0)
							continue;
						unsigned int bit;
						unsigned int p = 0;
						s[p++] = '0';
						s[p++] = 'b';
						for (bit = contextsize >> 1; bit; bit >>= 1)
							s[p++] = '0' + ((context & bit) ? 1 : 0);
						s[p] = 0;
						fprintf(stderr, "%i:%s %i %i\n", position, s, h[0], h[1]);
					}
				}
			}
		}
	}
	else if (!strcmp(argv[1], "decode_bcfp"))
	{
		return 1;
	}
	else
	{
		printusage(argc, argv);
		return -1;
	}
	return 0;
}


