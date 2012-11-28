#include "../lib/aos.h"
#include "../cuda/gpu.h"

int main(){
#define ngpu 2
    int gpus[ngpu]={2, 3};
    mvm_iwfs("mvm1.bin", "mvm2.bin", "grad1.bin", "grad2.bin", gpus, ngpu, 200000);
}
