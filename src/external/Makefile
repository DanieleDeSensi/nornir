.PHONY = all

all:
# TODO: Checkout correct versions of external libraries. After clone:
#&& cd mammut && git checkout 47ef67d33a998279345ea6f8741b7fa9ec9847b7
#&& cd riff  && git checkout aae8e80398a6a516705f0b5472db1e2881dae396
	if [ ! -d "./mammut" ]; then  \
		git clone https://github.com/DanieleDeSensi/mammut.git; \
	fi

	if [ ! -d "./riff" ]; then  \
		git clone https://github.com/DanieleDeSensi/riff.git; \
	fi

# Build mammut
	$(MAKE) -C mammut
# Build leo
ifeq ("$(wildcard ./leo/leo.o)","")
	$(MAKE) -C leo
endif
# Build queues
ifeq ("$(wildcard ./queues/hzdptr.o)","")
	$(MAKE) -C queues msqueue
endif
ifeq ("$(wildcard ./queues/xxhash.o)","")
	$(MAKE) -C queues msqueue
endif
# Build riff
	$(MAKE) -C riff
# Build graphite client
	$(MAKE) -C graphite-c-client all
clean:
	$(MAKE) -C leo clean
cleanall:
	$(MAKE) -C leo cleanall

