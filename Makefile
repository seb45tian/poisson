BUILDDIR := build

SRCDIR := src

HDRS := $(wildcard $(SRCDIR)/*.h)

SRCS := $(wildcard $(SRCDIR)/*.c)

CC=gcc-5
CFLAGS=-Wall -O3

LD := $(CC)

all: poisson

poisson: $(BUILDDIR)/alloc.o $(BUILDDIR)/simulation.o $(BUILDDIR)/poisson.o 
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BUILDDIR)/%.o : $(SRCDIR)/%.c
	@echo compiling $<
	$(maketargetdir)
	$(CC) $(CFLAGS) -c -o $@ $<

define maketargetdir
	-@mkdir -p $(dir $@) > /dev/null 2>&1
endef

clean:
	rm -f poisson $(OBJS)
	rm -rf $(BUILDDIR)
