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

BatchMandelCalculator::BatchMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "BatchMandelCalculator")
{
	data = (int *)(malloc(height * width * sizeof(int)));
}

BatchMandelCalculator::~BatchMandelCalculator() {
	free(data);
	data = NULL;
}


int * BatchMandelCalculator::calculateMandelbrot () {


	int half = height/2;
	int *fData = data; 
	int *bData = &data[(height)*width];
	

	int matrixSize = 16*16;
	int numberOfMatrix = (height*width)/matrixSize;
	
	for (int m = 0; m < numberOfMatrix; m++){
		int row = 16*(int)((16*(m))/width);
		int column = (16*m)%width;
		int mx = column + row*width;
		int *fData = &data[mx]; 
		for (int i = 0; i < 16; i++)
		{
			float x, y_s, y, x2, y2;
			y_s = y_start + (row+i) * dy; // current imaginary value
			#pragma omp simd 
			for (int j = 0; j < 16; j++)
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
				*(fData++) = value;
				// *(fData++) = *((bData--)-(width-(j*2)))  = value;	
			}
			fData += width - 16;
		}
	}
		return data;

}


