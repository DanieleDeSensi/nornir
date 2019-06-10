.PHONY = all

all:
# Build queues
ifeq ("$(wildcard ./queues/hzdptr.o)","")
	$(MAKE) -C queues msqueue
endif
ifeq ("$(wildcard ./queues/xxhash.o)","")
	$(MAKE) -C queues msqueue
endif
