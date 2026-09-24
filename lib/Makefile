USR := $(shell whoami)
CC=gcc

# Portable GSL detection via gsl-config (works on Linux, macOS/Homebrew, etc.)
GSL_CFLAGS := $(shell gsl-config --cflags)
GSL_LIBS   := $(shell gsl-config --libs)

CFLAGS = -Wall $(GSL_CFLAGS) -std=gnu99 -m64 -O3 -fpic
LDFLAGS = $(GSL_LIBS) -DHCUBATURE -lm

SOURCES=main.c hcubature.c source/halo.c source/effFormFact.c source/difCrosSec.c source/difRateGen.c source/WriteData.c source/coeffs_eft.c source/phys_consts.c source/binning_general.c source/efficiency_curve.c

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=main

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $(LDFLAGS) -o $@

.c:
	$(CC) $(CFLAGS) $< -c $@

clean:
	rm -f *.o $(EXECUTABLE)