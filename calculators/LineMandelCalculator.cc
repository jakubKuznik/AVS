/**
 * @file LineMandelCalculator.cc
 * @author FULL NAME <xlogin00@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over lines
 * @date DATE
 */
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <stdlib.h>
#include <cstring>

#include "LineMandelCalculator.h"

#include <immintrin.h>


LineMandelCalculator::LineMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "LineMandelCalculator")
{
	data = (int *)(_mm_malloc(height * width * sizeof(int), 64));
	real = (double *)(_mm_malloc(width * sizeof(double), 64));
	imag = (double *)(_mm_malloc(width * sizeof(double), 64));
	curr = (int *)(_mm_malloc(width * sizeof(int), 64));
}

LineMandelCalculator::~LineMandelCalculator() {
	_mm_free(data);
	data = NULL;
	_mm_free(real);
	real = NULL;
	_mm_free(imag);
	imag = NULL;
	_mm_free(curr);
	curr = NULL;
}

int * LineMandelCalculator::calculateMandelbrot () {

	int half = height/2;
	double *re = real;
	double *im = imag;
	int *cu = curr; 
	double x, y_s, y, x2, y2, x_s;

	int sum = 0;

	for (int i = 0; i < half; i++){
		y_s = y_start + i * dy; // current imaginary value
		// INIT VALUES

		std::uninitialized_fill(curr, curr+width, limit);
		std::uninitialized_fill(imag, imag+width, y_s);

		#pragma omp simd 
		for (int p = 0; p < width; p++){
			*re = x_start + p * dx;
			re++;
		}
		re = real;
		im = imag;

		for (int k = 0; (k < limit); k++){
			sum = 0; 	
			#pragma omp simd 
			for (int j = 0; j < width; j++){
				x_s = x_start + j*dx;

				x2 = (*re)*(*re);
				y2 = (*im)*(*im);
				if (x2 + y2 > 4.0){
					if (*(cu) == limit){
						*(cu) = k;
						sum++;
					}
				}
				*im = 2.0f * (*re) * (*im) + y_s;
				*re = x2 - y2 + x_s;
				im++;
				re++;
				cu++;
			}
			cu = curr;
			re = real;
			im = imag;
			if (sum == width)
				break;
		}
		memcpy(data+(i*width), curr, width*sizeof(int));
		memcpy(data+((width*(height-i-1))), curr, width*sizeof(int));
	}
	return data;
}


