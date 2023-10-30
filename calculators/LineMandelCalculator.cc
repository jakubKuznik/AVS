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

	// H == 2, W == 3 

	int half = height/2;
	int *pdata = data;
	for (int j = 0; j < width; j++)
	{
		float x_s = x_start + j *dx;
		#pragma omp simd 
		for (int i = 0; i < half; i++)
		{
			float x = x_s; 
			float y_s = y_start + i * dy; // current imaginary value
			float y = y_s; // current imaginary value

			int value = limit;
			for (int k=0;(k < limit);k++)
			{
				float x2 = x*x;
				float y2 = y*y;
				y = 2.0f * x * y + y_s;
				x = x2 - y2 + x_s; 
				if (x2 + y2 > 4.0f){
					value = k;
					break;
				}
			}	

			/// inverse 
			pdata[((height-i-1)*width)+j] = pdata[i*width+j] = value;  
		}
	}
	return data;

	
}


