#ifndef _GLOBALS_H //This check insures that C-preprocessor only includes
#define _GLOBALS_H //the following text once (not multiple times)


// The following variable defines the number of dimensions the 
// simulation runs in (usually 3).

//static const int g_dim = 3;

extern int g_dim;

static const int ERR_SYNTAX = 1;
static const int ERR_INPUT = 2;


#endif //#ifndef _GLOBALS_H
