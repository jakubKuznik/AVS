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


#include "LineMandelCalculator.h"


LineMandelCalculator::LineMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "LineMandelCalculator")
{
	data = (int *)(malloc(height * width * sizeof(int)));
}

LineMandelCalculator::~LineMandelCalculator() {
	free(data);
	data = NULL;
}

int * LineMandelCalculator::calculateMandelbrot () {

	int half = height/2;
	int *fData = data; 
	int *bData = data;
	for (int j = 0; j < width; j++)
	{
		float x_s = x_start + j *dx;
		float x, y_s, y, x2, y2;
		bData = &data[((height-1)*width)+j];
		fData = &data[j];
		#pragma omp simd 
		for (int i = 0; i < half; i++)
		{
			int p = i*width;
			x = x_s; 
			y_s = y_start + i * dy; // current imaginary value
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

			/// inverse 
			*(fData + p) = *(bData - p) = value;
			//data[((height-i-1)*width)+j] = data[i*width+j] = value;  
		}
	}
	return data;
	
}


