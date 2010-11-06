#ifndef AOS_SYS_PRIVATE_H
#define AOS_SYS_PRIVAYE_H
//Contains definitions that should not be accessed by the end user.
#if defined(__CYGWIN__)
void GetTempPath(long, char*);
#endif
void init_path(void);
#endif
