SUBDIRS = 	buchberger

all: $(SUBDIRS)

$(SUBDIRS):
		$(MAKE) -C $@

run: 
	./bin/benchmark_buchberger

.PHONY: all $(SUBDIRS) run
