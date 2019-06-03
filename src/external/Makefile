.PHONY = all

all:
# TODO: Checkout correct versions of external libraries. After clone:
	if [ ! -d "./mammut" ]; then  \
		git clone https://github.com/DanieleDeSensi/mammut.git; \
	fi

	if [ ! -d "./riff" ]; then  \
		git clone https://github.com/DanieleDeSensi/riff.git; \
	fi

# Build mammut
#	$(MAKE) -C mammut
# Build queues
ifeq ("$(wildcard ./queues/hzdptr.o)","")
	$(MAKE) -C queues msqueue
endif
ifeq ("$(wildcard ./queues/xxhash.o)","")
	$(MAKE) -C queues msqueue
endif
# Build riff
#	$(MAKE) -C riff
# Build graphite client
	$(MAKE) -C graphite-c-client all
clean:
	$(MAKE) -C leo clean
cleanall:
	$(MAKE) -C leo cleanall

