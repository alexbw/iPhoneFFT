//
//  OouraFFT.h
//  oScope
//
//  Created by Alex Wiltschko on 10/14/09.
//  This software is free because I love you.
//	But remember: Prof Ooura did the hard work.
//	A metaphor: Ooura raised the cow, slaughtered it, ground the beef,
//	made the beef sausages, and delivered it to my fridge.
//	I took it out of the fridge, microwaved it, 
//	slathered it with mustard and then called it dinner.
//
// Algorithm from fft4g.c
//
// Quirks and to-do's...
// * Rewrite to use single-precision. 
//	 I don't know what's going on in the bit reversal array, which might depend
//	 on a particular word length.
// * Rewrite to use UInt32 as length specifiers. It's silly we can't even analyze a whole second of data. Ah well.

// data[2*k] = R[k], 0<=k<n/2, real frequency data
// data[2*k+1] = I[k], 0<k<n/2, imaginary frequency data
// data[1] = R[n/2]



#import <Foundation/Foundation.h>
#import "fft4g.h"

// Define integer values for various windows
// NOTE: not using any windows that require parameter input.
#define kWindow_Rectangular	1
#define kWindow_Hamming		2
#define kWindow_Hann		3
#define kWindow_Bartlett	4
#define kWindow_Triangular	5



@interface OouraFFT : NSObject {
	double *inputData;
	
	int dataLength;
	int numFrequencies;
	BOOL dataIsFrequency;


	double *window;
	int windowType;
	int numWindows;
	int oldestDataIndex;
	
	double **allWindowedFFTs;
	double *spectrumData;
	// TODO: hold onto real and imaginary parts? no? maybe?
	
}

@property double *inputData;

@property double *window;
@property int windowType;
@property int numWindows;

@property int numFrequencies;
@property int dataLength;
@property double **allWindowedFFTs;
@property double *spectrumData;

- (void)doFFT;
- (void)doIFFT;

- (void)calculateWelchPeriodogramWithNewSignalSegment;
- (void)windowSignalSegment;

- (id)initForSignalsOfLength:(int)numPoints andNumWindows:(int)numberOfWindows;
- (BOOL)checkDataLength:(int)inDataLength;

- (void)setWindowType:(int)inWindowType;

@end