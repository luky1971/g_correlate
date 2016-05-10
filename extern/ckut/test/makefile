CC = gcc
CFLAGS += -std=c99

INCLUDE = ../include
SRC = .
BUILD = build

$(BUILD)/test_string: $(SRC)/test_string.c $(INCLUDE)/ckut_string.h
	$(CC) $(CFLAGS) -o $(BUILD)/test_string $(SRC)/test_string.c -I$(INCLUDE)

clean:
	rm -f $(BUILD)/*.o $(BUILD)/test_string
