#ifndef __NVME_CTRL_H__
#define __NVME_CTRL_H__
#ifdef __cplusplus
extern "C" {
#endif

#include "types.h"
#include "memory.h"
#include <stddef.h>
#include <stdint.h>


/*
 * Reset the NVM controller and initialize the controller structure.
 *
 * ctrl             pointer to controller handle
 *
 * ioctl_fd         open file descriptor to kernel module
 *                  used to call kernel functions for page-locking and aquiring
 *                  physical address of queues
 *
 * register_ptr     mmap'd BAR0 address space of the controller
 *
 * Returns 0 on success, or an error code on failure.
 *
 */
int nvm_init(nvm_ctrl_t* ctrl, int ioctl_fd, volatile void* register_ptr);


/*
 * Free allocated resources.
 */
void nvm_free(nvm_ctrl_t* ctrl, int ioctl_fd);


/*
 * Prepare queue handle and create IO completion queue (CQ)
 *
 * At the moment, only page-sized queues are supported.
 */
int nvm_create_cq(nvm_queue_t* cq, uint16_t no, nvm_ctrl_t* ctrl, void* vaddr, uint64_t paddr);


/* 
 * Prepare queue handle and create IO submission queue (SQ)
 *
 * Corresponding CQ must have been created first.
 * At the moment, only page-sized queues are supported.
 */
int nvm_create_sq(nvm_queue_t* sq, uint16_t cq_no, nvm_ctrl_t* ctrl, void* vaddr, uint64_t paddr);


#ifdef __cplusplus
}
#endif
#endif