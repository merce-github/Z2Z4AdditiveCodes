/************************************************************/
/*                                                          */
/* Package or project name: Z2Z4AdditiveCodes package       */
/* Regression test file name: Z2Z4AdditiveCodes_R_test.m    */
/*                                                          */
/* Comments:                                                */
/*                                                          */
/* Authors: J. Pujol and M. Villanueva                      */
/*                                                          */
/* Revision version and last date: v1.0, 2018/02/23         */
/*                                 v1.1, 2022/07/26         */
/*                                                          */
/************************************************************/

SetAssertions(true);
Alarm(30*60);
SetQuitOnError(false);

/************************************************************/
/*                                                          */
/* BLACK-BOX TESTS                                          */
/*                                                          */
/************************************************************/
load "Z2Z4AdditiveCode_BB_test.m";
load "Z2Z4StandardForm_BB_test.m";
load "Z2Z4NumInvariants_BB_test.m";
load "Z2Z4GrayMap_BB_test.m";
load "Z2Z4SpanKernelCodes_BB_test.m";
load "Z2Z4CyclicCodes_BB_test.m";
load "Z2Z4MinimumWeight_BB_test.m";
load "Z2Z4Decoding_BB_test.m";
load "Z2Z4Families_BB_test.m";
load "Z2Z4CoveringRadius_BB_test.m";



