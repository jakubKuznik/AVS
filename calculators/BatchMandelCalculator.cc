/**
 * @file BatchMandelCalculator.cc
 * @author FULL NAME <xlogin00@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over small batches
 * @date DATE
 */

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <stdlib.h>
#include <stdexcept>

#include "BatchMandelCalculator.h"

#include <immintrin.h>

#define MS 64

BatchMandelCalculator::BatchMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "BatchMandelCalculator")
{
	data = (int *)(_mm_malloc(height * width * sizeof(int), 64));
}

BatchMandelCalculator::~BatchMandelCalculator() {
	_mm_free(data);
	data = NULL;
}


int * BatchMandelCalculator::calculateMandelbrot () {
	
	int half = height/2;
	int matrixSize = MS*MS;
	int numberOfMatrix = (height*width)/matrixSize;
	
	for (int m = 0; m < numberOfMatrix/2; m++){
		int row = MS*(int)((MS*(m))/width);
		int column = (MS*m)%width;
		
		for (int i = 0; i < MS; i++)
		{
			float x, y_s, y, x2, y2;
			y_s = y_start + (row+i) * dy; // current imaginary value
			#pragma omp simd simdlen(64) 
			for (int j = 0; j < MS; j++)
			{
				float x_s = x_start + (j+column) *dx;
				x = x_s; 
				y = y_s; // current imaginary value

				int value = limit;
				for (int k = 0;(k < limit); k++)
				{
					x2 = x*x;
					y2 = y*y;
					y = 2.0f * x * y + y_s;
					x = x2 - y2 + x_s; 

					if (x2 + y2 > 4.0f){
						value = k;
						break;
					}
				}	
				data[((row+i)*width)+(j+column)] = value;
				data[(((height-1)*width)-((row+i)*width)) + (j+column)] = value;
			}
		}
	}
	return data;
}


