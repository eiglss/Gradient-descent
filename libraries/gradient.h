/******************************************************************************/
/**
*      @file : header.h
*   @version : 0
*      @date : month dd, yyyy
*    @author : Enzo IGLESIS
* @copyright : Copyright (c) 2019 Enzo IGLESIS
*              MIT License (MIT)
*******************************************************************************/

#ifndef HEADER_H_
#define HEADER_H_

/*******************************   LIBRARIES    *******************************/
#include <stdio.h>
#include <stdlib.h>

/*******************************     MACROS     *******************************/
/**** SYSTEM ****/
#define SYSTEM_SIZE 3

/**** GRADIENT DESCENT PARAMETERS *****/
#define GAMMA       0.0001
#define MAX_ITER    100000
#define FUNC_TOL    0.1

#define GD_VERBOSE   /* comment the line to not display information
                        of gradient descent */

//#define PTHREAD

/*******************************     TYPES      *******************************/
typedef struct
{
    double x[3];
}state_t;

/*******************************   CONSTANTES   *******************************/

/*******************************    VARIABLES   *******************************/

/*******************************   FUNCTIONS    *******************************/
void gradient_descent_loop(state_t * x);

#endif /* HEADER_H_ */
