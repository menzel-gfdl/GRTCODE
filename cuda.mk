CUDA_COMPILE = $(NVCC) -x cu -dc -ccbin $(CC) --compiler-options="$(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS)"
#CUDA_LIBTOOL_COMPILE = $(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile $(CUDA_COMPILE) -prefer-non-pic
#CUDA_LINK = $(NVCC) --compiler-options="$(AM_CFLAGS) $(CFLAGS)" --linker-options="$(AM_LDFLAGS) $(LDFLAGS)"
#CUDA_LIBTOOL_LINK = $(AM_V_CCLD) $(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=link $(CUDA_LINK)

.c.o:
	$(CUDA_COMPILE) -o $@ -c $<

.cu.o:
	$(CUDA_COMPILE) -o $@ -c $<

$(gpu_object): $(objects)
	$(NVCC) -dlink -o $@ $^


#.c.lo:
#	$(CUDA_LIBTOOL_COMPILE) -o $@ -c $<

#.cu.lo:
#	$(CUDA_LIBTOOL_COMPILE) -o $@ -c $<
