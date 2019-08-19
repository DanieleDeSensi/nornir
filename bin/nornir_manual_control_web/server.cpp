#include <assert.h>
#include <stdio.h>
#include "../../src/external/riff/src/external/cppnanomsg/nn.hpp"
#include "../../src/external/riff/src/external/nanomsg/src/pair.h"


int main (const int argc, const char **argv){
  int sock = nn_socket (AF_SP, NN_PAIR);
  assert(sock >= 0);
  int rbind = nn_connect (sock, "ws://0.0.0.0:3001");
  printf("nn_bind returned %d\n", rbind);
  assert(rbind >= 0);
  while(1){
      char *buf = NULL;
      int bytes = nn_recv (sock, &buf, NN_MSG, 0);
      assert (bytes >= 0);
      printf ("NODE0: RECEIVED \"%s\"\n", buf);
      buf[0] = 0;
      nn_freemsg (buf);
    }
}