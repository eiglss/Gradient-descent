/******************************************************************************/
/**
*      @file c.c
*   @version 0
*      @date month dd, yyyy
*    @author Enzo IGLESIS
* @copyright Copyright (c) 2019 Enzo IGLESIS
*              MIT License (MIT)
*******************************************************************************/

/*******************************    LIBRARIES    *******************************/
#include "gradient.h"
#include <math.h>

/*******************************     MACROS     *******************************/

/*******************************      TYPES     *******************************/

/*******************************    CONSTANTS   *******************************/

/*******************************    VARIABLES   *******************************/

/*******************************    FUNCTIONS   *******************************/
/**** SYSTEM DEPENDENT ****/

/**** function ****/
static inline double function_0(state_t * x)
{
    return 3.*x->x[0]-cos(x->x[1]*x->x[2])-(3./2.);
}

static inline double function_1(state_t * x)
{
    return 4.*(x->x[0]*x->x[0])-625.*(x->x[1]*x->x[1])+2.*x->x[1]-1.;
}

static inline double function_2(state_t * x)
{
    return exp(-x->x[0]*x->x[1])+20.*x->x[2]+(10.*M_PI-3.)/3.;
}

static double (*function[SYSTEM_SIZE])(state_t * x) = {function_0,\
                                                       function_1,\
                                                       function_2};

/**** jacobian ****/
/* row 1 column 1 */
static inline double jacobian_0_0(state_t * x)
{
    return 3.;
}

/* row 1 column 2 */
static inline double jacobian_0_1(state_t * x)
{
    return sin(x->x[1]*x->x[2])*x->x[2];
}

/* row 1 column 3 */
static inline double jacobian_0_2(state_t * x)
{
    return sin(x->x[1]*x->x[2])*x->x[1];
}

/* row 2 column 1 */
static inline double jacobian_1_0(state_t * x)
{
    return 8.*x->x[0];
}

/* row 2 column 2 */
static inline double jacobian_1_1(state_t * x)
{
    return -1250.*x->x[1]+2.;
}

/* row 2 column 3 */
static inline double jacobian_1_2(state_t * x)
{
    return 0.;
}

/* row 3 column 1 */
static inline double jacobian_2_0(state_t * x)
{
    return -x->x[1]*exp(-x->x[0]*x->x[1]);
}

/* row 3 column 2 */
static inline double jacobian_2_1(state_t * x)
{
    return -x->x[0]*exp(-x->x[0]*x->x[1]);
}

/* row 3 column 3 */
static inline double jacobian_2_2(state_t * x)
{
    return 20.;
}

static double (*jacobian[SYSTEM_SIZE][SYSTEM_SIZE])(state_t * x) = {{jacobian_0_0,\
                                                                     jacobian_0_1,\
                                                                     jacobian_0_2},\
                                                                    {jacobian_1_0,\
                                                                     jacobian_1_1,\
                                                                     jacobian_1_2},\
                                                                    {jacobian_2_0,\
                                                                     jacobian_2_1,\
                                                                     jacobian_2_2}};

/**** gradient of the function ****/
static inline double gradient_function_0(state_t * x)
{
    return jacobian[0][0](x)*function[0](x)+\
           jacobian[1][0](x)*function[1](x)+\
           jacobian[2][0](x)*function[2](x);
}

static inline double gradient_function_1(state_t * x)
{
    return jacobian[0][1](x)*function[0](x)+\
           jacobian[1][1](x)*function[1](x)+\
           jacobian[2][1](x)*function[2](x);
}

static inline double gradient_function_2(state_t * x)
{
    return jacobian[0][2](x)*function[0](x)+\
           jacobian[1][2](x)*function[1](x)+\
           jacobian[2][2](x)*function[2](x);
}

static double (*gradient_function[SYSTEM_SIZE])(state_t * x) = {gradient_function_0,\
                                                                gradient_function_1,\
                                                                gradient_function_2};

/**** objective function ****/
static inline double objective_function(state_t * x)
{
    /* declaration */
    size_t i;
    double sq_y[SYSTEM_SIZE];
    double sum = 0;
    /* program */
    for(i = 0; i < SYSTEM_SIZE; i++) /* multi-threading possible */
    {
        sq_y[i] = function[i](x)*function[i](x);
    }
    for(i = 0; i < SYSTEM_SIZE; i++) /* multi-threading impossible */
    {
        sum += sq_y[i];
    }
    return 0.5*sum;
}

/**** gradient descent ****/
static inline double gradient_descent_0(state_t * x)
{
    return x->x[0] = x->x[0]-GAMMA*gradient_function[0](x);
}

static inline double gradient_descent_1(state_t * x)
{
    return x->x[1] = x->x[1]-GAMMA*gradient_function[1](x);
}

static inline double gradient_descent_2(state_t * x)
{
    return x->x[2] = x->x[2]-GAMMA*gradient_function[2](x);
}

static double (*gradient_descent[SYSTEM_SIZE])(state_t * x) = {gradient_descent_0,\
                                                               gradient_descent_1,\
                                                               gradient_descent_2};

void gradient_descent_loop(state_t * x)
{
    /* declaration */
    size_t iter = 0;
    size_t i;
    /* program */
    do
    {
        #ifdef GD_VERBOSE
            printf("[%04lu]: x[%f, %f, %f]; F(x) = %f\r\n", iter, x->x[0], x->x[1], x->x[2], objective_function(x));
        #endif
        for(i = 0; i < SYSTEM_SIZE; i++) /* multi-threading possible */
        {
            gradient_descent[i](x);
        }
    }
    while(++iter < MAX_ITER && objective_function(x) > FUNC_TOL); /* multi-threading impossible */
    #ifdef GD_VERBOSE
        printf("[%04lu]: x[%f, %f, %f]; F(x) = %f\r\n", iter, x->x[0], x->x[1], x->x[2], objective_function(x));
    #endif
}
